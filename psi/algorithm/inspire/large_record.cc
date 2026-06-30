#include "psi/algorithm/inspire/large_record.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <vector>

#include "yacl/base/exception.h"

#include "psi/algorithm/inspire/internal_primitives.h"
#include "psi/algorithm/inspire/large_record_he.h"
#include "psi/algorithm/ypir/legacy/ypir_util.h"

using namespace psi::inspire::lrhe;

namespace psi::inspire::internal {

namespace {

using psi::ypir::ypir_internal::Cdks21Lwe2RlweInplace;
using psi::ypir::ypir_internal::EltwiseAddMod;
using psi::ypir::ypir_internal::EltwiseFMAMod;
using psi::ypir::ypir_internal::EltwiseMultMod;
using psi::ypir::ypir_internal::EltwiseSubMod;
using psi::ypir::ypir_internal::IsPowerOfTwo;
using psi::ypir::ypir_internal::MatrixMultiplicationFlatU16;
using psi::ypir::ypir_internal::MatrixVectorMultiplicationU16;
using psi::ypir::ypir_internal::ModInverse;
using psi::ypir::ypir_internal::PirParams;
using psi::ypir::ypir_internal::PseudorandomMatrixGenerate;
using psi::ypir::ypir_internal::Secret;
using psi::ypir::ypir_internal::YpirHexlNtt;
using psi::ypir::ypir_internal::kThirdDimensionSeed;

using Record = std::array<uint8_t, LargeRecordConfig::kRecordBytes>;

uint64_t ModAdd(uint64_t lhs, uint64_t rhs, uint64_t mod) {
  uint64_t sum = lhs + rhs;
  if (sum >= mod) {
    sum -= mod;
  }
  return sum;
}

uint64_t ModMul(uint64_t lhs, uint64_t rhs, uint64_t mod) {
  return static_cast<uint64_t>((static_cast<__uint128_t>(lhs) * rhs) % mod);
}

void NegacyclicMonomialMul(const std::vector<uint64_t>& input, uint64_t shift,
                           uint64_t mod, std::vector<uint64_t>& output) {
  const uint64_t degree = input.size();
  std::fill(output.begin(), output.end(), 0);
  for (uint64_t i = 0; i < degree; ++i) {
    uint64_t exponent = i + shift;
    bool neg = false;
    while (exponent >= degree) {
      exponent -= degree;
      neg = !neg;
    }
    uint64_t value = input[i] % mod;
    if (!neg) {
      output[exponent] = ModAdd(output[exponent], value, mod);
    } else if (value != 0) {
      output[exponent] = (output[exponent] + mod - value) % mod;
    }
  }
}

std::vector<std::vector<uint64_t>> InterpolateColumn(
    const std::vector<std::vector<uint64_t>>& values,
    const LargeRecordConfig& cfg) {
  const uint64_t row = cfg.folding;
  const uint64_t degree = cfg.degree;
  const uint64_t mod = cfg.pmod;
  const uint64_t step = (degree << 1) / row;
  const uint64_t row_inv = static_cast<uint64_t>(
      ModInverse(static_cast<int64_t>(row), static_cast<int64_t>(mod)));

  std::vector<std::vector<uint64_t>> coeffs(
      row, std::vector<uint64_t>(degree, 0));
  std::vector<uint64_t> rotated(degree, 0);

  for (uint64_t j = 0; j < row; ++j) {
    for (uint64_t r = 0; r < row; ++r) {
      uint64_t power = (j * r) % row;
      uint64_t shift = ((row - power) % row) * step;
      NegacyclicMonomialMul(values[r], shift, mod, rotated);
      for (uint64_t k = 0; k < degree; ++k) {
        coeffs[j][k] = ModAdd(coeffs[j][k], rotated[k], mod);
      }
    }
    for (uint64_t k = 0; k < degree; ++k) {
      coeffs[j][k] = ModMul(coeffs[j][k], row_inv, mod);
    }
  }
  return coeffs;
}

uint64_t GetBundleFactor(const LargeRecordConfig& cfg) {
  return cfg.packed_coeffs / LargeRecordConfig::kRecordBytes;
}

Record GenerateRecord(uint64_t record_idx) {
  Record record{};
  for (uint64_t b = 0; b < LargeRecordConfig::kRecordBytes; ++b) {
    record[b] = static_cast<uint8_t>((17 * record_idx + b) & 0xFFULL);
  }
  return record;
}

std::vector<uint64_t> BuildBundledEntry(uint64_t bundle_idx,
                                        const LargeRecordConfig& cfg) {
  const uint64_t bundle_factor = GetBundleFactor(cfg);
  std::vector<uint64_t> entry(cfg.degree, 0);
  for (uint64_t slot = 0; slot < bundle_factor; ++slot) {
    Record record = GenerateRecord(bundle_idx * bundle_factor + slot);
    uint64_t base = slot * LargeRecordConfig::kRecordBytes;
    for (uint64_t b = 0; b < LargeRecordConfig::kRecordBytes; ++b) {
      entry[base + b] = static_cast<uint64_t>(record[b]);
    }
  }
  return entry;
}

void EncodeDatabaseLargeRecord(const LargeRecordConfig& cfg,
                               std::vector<std::vector<uint16_t>>& encoded_db) {
  const uint64_t bundle_factor = GetBundleFactor(cfg);
  const uint64_t bundled_entries = cfg.total_records / bundle_factor;
  const uint64_t t = cfg.folding;
  const uint64_t cols = bundled_entries / t;
  const uint64_t rows = t * cfg.packed_coeffs;

  encoded_db.assign(rows, std::vector<uint16_t>(cols, 0));

#pragma omp parallel for
  for (uint64_t col = 0; col < cols; ++col) {
    std::vector<std::vector<uint64_t>> values(
        t, std::vector<uint64_t>(cfg.degree, 0));
    for (uint64_t row = 0; row < t; ++row) {
      values[row] = BuildBundledEntry(col * t + row, cfg);
    }
    std::vector<std::vector<uint64_t>> coeffs = InterpolateColumn(values, cfg);
    for (uint64_t k = 0; k < t; ++k) {
      for (uint64_t i = 0; i < cfg.packed_coeffs; ++i) {
        encoded_db[k * cfg.packed_coeffs + i][col] =
            static_cast<uint16_t>(coeffs[k][i] % cfg.pmod);
      }
    }
  }
}

void BuildPaperPackKeyMaterial(
    Secret& sk, AESCTR_PRNG& prng,
    const psi::ypir::ypir_internal::FheParams& fparm,
    PaperPackKeyMaterial& key_material) {
  const uint64_t degree = fparm.get_poly_degree();
  const uint64_t cmod = fparm.get_rlwe_cmod();
  const uint64_t t_auto = fparm.get_t_auto();
  const uint64_t half = degree / 2;
  const uint64_t g = 5;
  const uint64_t h = (degree << 1) - 1;

  std::vector<std::vector<uint64_t>> wg(
      t_auto, std::vector<uint64_t>(degree, 0));
  std::vector<std::vector<uint64_t>> yg(
      t_auto, std::vector<uint64_t>(degree, 0));
  std::vector<std::vector<uint64_t>> wh(
      t_auto, std::vector<uint64_t>(degree, 0));
  std::vector<std::vector<uint64_t>> yh(
      t_auto, std::vector<uint64_t>(degree, 0));

  prng.refresh(kThirdDimensionSeed);
  PseudorandomMatrixGenerate(wg, cmod, prng);
  PseudorandomMatrixGenerate(wh, cmod, prng);

  std::vector<uint64_t> sk_copy = sk.data;
  std::vector<uint64_t> g_key(degree, 0);
  std::vector<uint64_t> h_key(degree, 0);
  ApplyAutoCoefForm(g_key, sk_copy, static_cast<int32_t>(g), cmod);
  ApplyAutoCoefForm(h_key, sk_copy, static_cast<int32_t>(h), cmod);
  GadgetEncrypt(sk, wg, g_key, yg, fparm);
  GadgetEncrypt(sk, wh, h_key, yh, fparm);

  key_material.g_automaps.assign(half, std::vector<uint32_t>(degree, 0));
  key_material.gh_automaps.assign(half, std::vector<uint32_t>(degree, 0));
  uint64_t g_power = 1;
  for (uint64_t i = 0; i < half; ++i) {
    psi::ypir::ypir_internal::PrecomputeAutomap(
        static_cast<uint32_t>(degree), static_cast<uint32_t>(g_power),
        key_material.g_automaps[i]);
    psi::ypir::ypir_internal::PrecomputeAutomap(
        static_cast<uint32_t>(degree),
        static_cast<uint32_t>(
            (static_cast<__uint128_t>(h) * g_power) % (degree << 1)),
        key_material.gh_automaps[i]);
    g_power = static_cast<uint64_t>(
        (static_cast<__uint128_t>(g_power) * g) % (degree << 1));
  }

  key_material.kg_a_powers.assign(
      half - 1, std::vector<std::vector<uint64_t>>(
                    t_auto, std::vector<uint64_t>(degree, 0)));
  key_material.kg_b_powers.assign(
      half - 1, std::vector<std::vector<uint64_t>>(
                    t_auto, std::vector<uint64_t>(degree, 0)));
  key_material.gh_kg_a_powers.assign(
      half - 1, std::vector<std::vector<uint64_t>>(
                    t_auto, std::vector<uint64_t>(degree, 0)));
  key_material.gh_kg_b_powers.assign(
      half - 1, std::vector<std::vector<uint64_t>>(
                    t_auto, std::vector<uint64_t>(degree, 0)));
  for (uint64_t i = 0; i + 1 < half; ++i) {
    RotateKskRows(wg, key_material.g_automaps[i],
                  key_material.kg_a_powers[i]);
    RotateKskRows(yg, key_material.g_automaps[i],
                  key_material.kg_b_powers[i]);
    RotateKskRows(wg, key_material.gh_automaps[i],
                  key_material.gh_kg_a_powers[i]);
    RotateKskRows(yg, key_material.gh_automaps[i],
                  key_material.gh_kg_b_powers[i]);
  }

  key_material.kh_a = std::move(wh);
  key_material.kh_b = std::move(yh);
  BuildNttMonomials(degree, degree, fparm.get_ntt(),
                    key_material.monomials_ntt);
}

void AggregatePaperTransformHalves(
    const std::vector<std::vector<uint64_t>>& block_rows,
    const PaperPackKeyMaterial& key_material,
    const psi::ypir::ypir_internal::FheParams& fparm,
    std::vector<std::vector<uint64_t>>& first_half,
    std::vector<std::vector<uint64_t>>& second_half) {
  const uint64_t degree = fparm.get_poly_degree();
  const uint64_t modulus = fparm.get_rlwe_cmod();
  const uint64_t half = degree / 2;
  const uint64_t degree_inv = static_cast<uint64_t>(
      ModInverse(static_cast<int64_t>(degree), static_cast<int64_t>(modulus)));
  YpirHexlNtt& ntts = fparm.get_ntt();

  first_half.assign(half, std::vector<uint64_t>(degree, 0));
  second_half.assign(half, std::vector<uint64_t>(degree, 0));

  std::vector<uint64_t> a_tilde(degree, 0);
  std::vector<uint64_t> rotated(degree, 0);
  std::vector<uint64_t> shifted(degree, 0);

  for (uint64_t k = 0; k < degree; ++k) {
    a_tilde = block_rows[k];
    Cdks21Lwe2RlweInplace(a_tilde.data(), degree, modulus, ntts);
    EltwiseFMAMod(a_tilde.data(), a_tilde.data(), degree_inv, nullptr, degree,
                  modulus);

    for (uint64_t j = 0; j < half; ++j) {
      ApplyAutoNttFormGeneral(a_tilde, rotated, key_material.g_automaps[j]);
      EltwiseMultMod(shifted.data(), rotated.data(),
                     key_material.monomials_ntt[k].data(), degree, modulus);
      EltwiseAddMod(first_half[j].data(), first_half[j].data(),
                    shifted.data(), degree, modulus);

      ApplyAutoNttFormGeneral(a_tilde, rotated, key_material.gh_automaps[j]);
      EltwiseMultMod(shifted.data(), rotated.data(),
                     key_material.monomials_ntt[k].data(), degree, modulus);
      EltwiseAddMod(second_half[j].data(), second_half[j].data(),
                    shifted.data(), degree, modulus);
    }
  }
}

void CollapseOnePreprocessPaper(
    std::vector<std::vector<uint64_t>>& current,
    const std::vector<std::vector<uint64_t>>& ksk_a,
    std::vector<std::vector<std::vector<uint64_t>>>& decomp_steps,
    const psi::ypir::ypir_internal::FheParams& fparm) {
  const uint64_t degree = fparm.get_poly_degree();
  const uint64_t modulus = fparm.get_rlwe_cmod();

  std::vector<uint64_t> switched(degree, 0);
  std::vector<std::vector<std::vector<uint64_t>>> local_decomp;
  std::vector<uint64_t> last = current.back();
  KeyswitchPreprocess(ksk_a, last, switched, local_decomp, fparm);
  current.pop_back();
  EltwiseAddMod(current.back().data(), current.back().data(),
                switched.data(), degree, modulus);
  decomp_steps.push_back(std::move(local_decomp[0]));
}

void CollapseHalfPreprocessPaper(
    std::vector<std::vector<uint64_t>> current,
    const std::vector<std::vector<std::vector<uint64_t>>>& rotated_ksk_a,
    std::vector<std::vector<std::vector<uint64_t>>>& decomp_steps,
    std::vector<uint64_t>& final_a_ntt,
    const psi::ypir::ypir_internal::FheParams& fparm) {
  decomp_steps.clear();
  while (current.size() > 1) {
    CollapseOnePreprocessPaper(current, rotated_ksk_a[current.size() - 2],
                               decomp_steps, fparm);
  }
  final_a_ntt = current[0];
}

void CollapseOneOnlinePaper(
    std::vector<uint64_t>& current_b_ntt,
    const std::vector<std::vector<uint64_t>>& ksk_b,
    std::vector<std::vector<uint64_t>>& decomp_step,
    const psi::ypir::ypir_internal::FheParams& fparm) {
  std::vector<uint64_t> next_b(current_b_ntt.size(), 0);
  KeyswitchOnline(ksk_b, current_b_ntt, next_b, decomp_step, fparm);
  current_b_ntt.swap(next_b);
}

void CollapseHalfOnlinePaper(
    std::vector<uint64_t>& current_b_ntt,
    const std::vector<std::vector<std::vector<uint64_t>>>& rotated_ksk_b,
    const std::vector<std::vector<std::vector<uint64_t>>>& decomp_steps,
    const psi::ypir::ypir_internal::FheParams& fparm) {
  const uint64_t num_steps = decomp_steps.size();
  for (uint64_t i = 0; i < num_steps; ++i) {
    CollapseOneOnlinePaper(
        current_b_ntt, rotated_ksk_b[num_steps - 1 - i],
        const_cast<std::vector<std::vector<uint64_t>>&>(decomp_steps[i]),
        fparm);
  }
}

void PreprocessPackedHintsPaper(
    const std::vector<std::vector<uint64_t>>& hint,
    const PaperPackKeyMaterial& key_material,
    std::vector<PaperPackedHintBlock>& packed_hint_blocks,
    const psi::ypir::ypir_internal::FheParams& fparm,
    const LargeRecordConfig& cfg) {
  packed_hint_blocks.resize(cfg.folding);

  for (uint64_t block = 0; block < cfg.folding; ++block) {
    std::vector<std::vector<uint64_t>> block_rows(
        cfg.packed_coeffs, std::vector<uint64_t>(cfg.degree, 0));
    for (uint64_t i = 0; i < cfg.packed_coeffs; ++i) {
      block_rows[i] = hint[block * cfg.packed_coeffs + i];
    }

    std::vector<std::vector<uint64_t>> first_half;
    std::vector<std::vector<uint64_t>> second_half;
    AggregatePaperTransformHalves(block_rows, key_material, fparm, first_half,
                                  second_half);

    std::vector<uint64_t> a1_ntt;
    std::vector<uint64_t> a2_ntt;
    CollapseHalfPreprocessPaper(
        std::move(first_half), key_material.kg_a_powers,
        packed_hint_blocks[block].decomp_half1, a1_ntt, fparm);
    CollapseHalfPreprocessPaper(
        std::move(second_half), key_material.gh_kg_a_powers,
        packed_hint_blocks[block].decomp_half2, a2_ntt, fparm);

    std::vector<std::vector<uint64_t>> final_pair;
    final_pair.push_back(std::move(a1_ntt));
    final_pair.push_back(std::move(a2_ntt));
    std::vector<std::vector<std::vector<uint64_t>>> final_decomp_steps;
    CollapseOnePreprocessPaper(final_pair, key_material.kh_a,
                               final_decomp_steps, fparm);
    packed_hint_blocks[block].packed_a_ntt = std::move(final_pair[0]);
    packed_hint_blocks[block].decomp_final = std::move(final_decomp_steps[0]);
  }
}

RlweCiphertext PackSelectedBlockOnlinePaper(
    uint64_t block_idx, const std::vector<uint64_t>& selected_column_response,
    const PaperPackKeyMaterial& key_material,
    const std::vector<PaperPackedHintBlock>& packed_hint_blocks,
    const psi::ypir::ypir_internal::FheParams& fparm,
    const LargeRecordConfig& cfg) {
  std::vector<uint64_t> b_agg(cfg.degree, 0);
  for (uint64_t i = 0; i < cfg.packed_coeffs; ++i) {
    b_agg[i] =
        selected_column_response[block_idx * cfg.packed_coeffs + i] % cfg.qmod;
  }
  fparm.get_ntt().Forward(b_agg.data(), cfg.degree);

  CollapseHalfOnlinePaper(b_agg, key_material.kg_b_powers,
                          packed_hint_blocks[block_idx].decomp_half1, fparm);
  CollapseHalfOnlinePaper(b_agg, key_material.gh_kg_b_powers,
                          packed_hint_blocks[block_idx].decomp_half2, fparm);
  CollapseOneOnlinePaper(
      b_agg, key_material.kh_b,
      const_cast<std::vector<std::vector<uint64_t>>&>(
          packed_hint_blocks[block_idx].decomp_final),
      fparm);

  RlweCiphertext ct(cfg.degree, cfg.qmod, true);
  ct.a = Poly(packed_hint_blocks[block_idx].packed_a_ntt, cfg.qmod, true);
  ct.b = Poly(b_agg, cfg.qmod, true);
  return ct;
}

RlweCiphertext AnswerLargeRecordOnlinePaper(
    const std::vector<uint64_t>& selected_column_response,
    const PaperPackKeyMaterial& key_material,
    const std::vector<PaperPackedHintBlock>& packed_hint_blocks,
    ApproximateRgswCiphertext& point_ct, intel::hexl::NTT& ntts,
    const psi::ypir::ypir_internal::FheParams& fparm,
    const LargeRecordConfig& cfg) {
  std::vector<RlweCiphertext> coeff_cts;
  coeff_cts.reserve(cfg.folding);
  for (uint64_t k = 0; k < cfg.folding; ++k) {
    coeff_cts.push_back(PackSelectedBlockOnlinePaper(
        k, selected_column_response, key_material, packed_hint_blocks, fparm,
        cfg));
  }

  RlweCiphertext acc = coeff_cts.back();
  for (int64_t j = static_cast<int64_t>(coeff_cts.size()) - 2; j >= 0; --j) {
    RlweCiphertext mul_in = acc;
    RlweCiphertext mul_out(acc.get_degree(), acc.get_modulus(), true);
    ExternalProduct(point_ct, mul_in, mul_out, ntts);
    RlweBfvAdd(
        mul_out, const_cast<RlweCiphertext&>(coeff_cts[static_cast<size_t>(j)]),
        acc);
  }
  return acc;
}

ApproximateRgswCiphertext EncryptPointQuery(uint64_t row_idx,
                                            lrhe::Secret& rlwe_sk,
                                            const LargeRecordConfig& cfg) {
  std::vector<uint64_t> point_poly(cfg.degree, 0);
  uint64_t shift = (row_idx % cfg.folding) * ((cfg.degree << 1) / cfg.folding);
  uint64_t exponent = shift;
  bool neg = false;
  while (exponent >= cfg.degree) {
    exponent -= cfg.degree;
    neg = !neg;
  }
  point_poly[exponent] = neg ? (cfg.qmod - 1) : 1;

  Plaintext point_pt(point_poly, cfg.qmod, false);
  ApproximateRgswCiphertext point_ct(cfg.degree, cfg.qmod, cfg.rgsw_b,
                                     cfg.rgsw_z, cfg.rgsw_t);
  RgswEncode(point_pt, point_ct, rlwe_sk, cfg.sigma);
  return point_ct;
}

Record ExtractRecordFromBundle(const std::vector<uint64_t>& bundle,
                               uint64_t local_offset) {
  uint64_t base = local_offset * LargeRecordConfig::kRecordBytes;
  Record record{};
  for (uint64_t i = 0; i < LargeRecordConfig::kRecordBytes; ++i) {
    record[i] = static_cast<uint8_t>(bundle[base + i] & 0xFFULL);
  }
  return record;
}

void SerializeRgswCt(
    const ApproximateRgswCiphertext& ct,
    std::vector<std::vector<std::vector<uint64_t>>>& ct_m_out,
    std::vector<std::vector<std::vector<uint64_t>>>& ct_sm_out) {
  ct_m_out.resize(ct.get_t());
  ct_sm_out.resize(ct.get_t());
  for (uint64_t i = 0; i < ct.get_t(); ++i) {
    ct_m_out[i] = {ct.ct_m[i].a.get_data(), ct.ct_m[i].b.get_data()};
    ct_sm_out[i] = {ct.ct_sm[i].a.get_data(), ct.ct_sm[i].b.get_data()};
  }
}

ApproximateRgswCiphertext DeserializeRgswCt(
    const std::vector<std::vector<std::vector<uint64_t>>>& ct_m_ser,
    const std::vector<std::vector<std::vector<uint64_t>>>& ct_sm_ser,
    const LargeRecordConfig& cfg) {
  ApproximateRgswCiphertext ct(cfg.degree, cfg.qmod, cfg.rgsw_b, cfg.rgsw_z,
                               cfg.rgsw_t);
  for (uint64_t i = 0; i < ct.get_t(); ++i) {
    ct.ct_m[i].a = Poly(ct_m_ser[i][0], cfg.qmod, true);
    ct.ct_m[i].b = Poly(ct_m_ser[i][1], cfg.qmod, true);
    ct.ct_sm[i].a = Poly(ct_sm_ser[i][0], cfg.qmod, true);
    ct.ct_sm[i].b = Poly(ct_sm_ser[i][1], cfg.qmod, true);
  }
  return ct;
}

}  // namespace

LargeRecordPrecomputedState PrepareLargeRecordState(
    const LargeRecordConfig& cfg) {
  const uint64_t bundle_factor = GetBundleFactor(cfg);
  const uint64_t bundled_entries = cfg.total_records / bundle_factor;
  const uint64_t cols = bundled_entries / cfg.folding;
  const uint64_t encoded_rows = cfg.folding * cfg.packed_coeffs;

  psi::ypir::ypir_internal::FheParams fparm(
      cfg.degree, cfg.qmod, cfg.pmod, cfg.lwe_dimension, cfg.lwe_cmod,
      cfg.lwe_pmod, 0.0, 3.19, {20, 18, 2}, {16, 16, 1});
  [[maybe_unused]] PirParams pparm(encoded_rows, cols);

  const uint8_t key[16] = "I am xct's son";
  AESCTR_PRNG prng(key);

  srand(123);
  std::vector<int64_t> sk_vec_signed(cfg.degree, 0);
  std::vector<uint64_t> sk_vec(cfg.degree, 0);
  for (uint64_t i = 0; i < cfg.degree; ++i) {
    sk_vec_signed[i] = rand() & 1;
    sk_vec[i] = static_cast<uint64_t>(sk_vec_signed[i]);
  }

  Secret packing_sk(cfg.degree, cfg.qmod);
  packing_sk.data = sk_vec;

  LargeRecordPrecomputedState state;

  EncodeDatabaseLargeRecord(cfg, state.encoded_db);
  fparm.set_persudo_matrix_simplepir(cols);

  state.hint.assign(encoded_rows,
                    std::vector<uint64_t>(cfg.lwe_dimension, 0));
  MatrixMultiplicationFlatU16(
      state.hint, state.encoded_db, fparm.get_persudo_matrix_simplepir_flat(),
      cfg.lwe_dimension, cfg.lwe_cmod);

  BuildPaperPackKeyMaterial(packing_sk, prng, fparm, state.key_material);
  PreprocessPackedHintsPaper(state.hint, state.key_material,
                             state.packed_hint_blocks, fparm, cfg);
  return state;
}

LargeRecordQuery GenerateLargeRecordQuery(
    uint64_t target_index, const LargeRecordConfig& cfg,
    LargeRecordClientSecrets& secrets) {
  const uint64_t bundle_factor = GetBundleFactor(cfg);
  const uint64_t bundled_entries = cfg.total_records / bundle_factor;
  const uint64_t cols = bundled_entries / cfg.folding;
  const uint64_t encoded_rows = cfg.folding * cfg.packed_coeffs;

  psi::ypir::ypir_internal::FheParams fparm(
      cfg.degree, cfg.qmod, cfg.pmod, cfg.lwe_dimension, cfg.lwe_cmod,
      cfg.lwe_pmod, 0.0, 3.19, {20, 18, 2}, {16, 16, 1});
  [[maybe_unused]] PirParams pparm(encoded_rows, cols);

  const uint64_t bundle_idx = target_index / bundle_factor;
  const uint64_t col_idx = bundle_idx / cfg.folding;
  const uint64_t row_idx = bundle_idx % cfg.folding;

  srand(123);
  secrets.sk_vec_signed.resize(cfg.degree, 0);
  secrets.sk_vec.resize(cfg.degree, 0);
  for (uint64_t i = 0; i < cfg.degree; ++i) {
    secrets.sk_vec_signed[i] = rand() & 1;
    secrets.sk_vec[i] = static_cast<uint64_t>(secrets.sk_vec_signed[i]);
  }

  Secret simple_sk(cfg.lwe_dimension, cfg.lwe_cmod);
  simple_sk.data = secrets.sk_vec;

  fparm.set_persudo_matrix_simplepir(cols);

  LargeRecordQuery query;
  SimplepirQuery(col_idx, simple_sk, query.column_query, fparm, pparm);

  lrhe::Secret rlwe_sk(secrets.sk_vec_signed, cfg.qmod, false);
  auto point_ct = EncryptPointQuery(row_idx, rlwe_sk, cfg);
  SerializeRgswCt(point_ct, query.rgsw_ct_m, query.rgsw_ct_sm);

  return query;
}

LargeRecordResponse ProcessLargeRecordQuery(
    const LargeRecordQuery& query, const LargeRecordPrecomputedState& state,
    const LargeRecordConfig& cfg) {
  psi::ypir::ypir_internal::FheParams fparm(
      cfg.degree, cfg.qmod, cfg.pmod, cfg.lwe_dimension, cfg.lwe_cmod,
      cfg.lwe_pmod, 0.0, 3.19, {20, 18, 2}, {16, 16, 1});

  std::vector<uint64_t> selected_column_response;
  MatrixVectorMultiplicationU16(selected_column_response, state.encoded_db,
                                query.column_query, cfg.lwe_cmod);

  auto point_ct = DeserializeRgswCt(query.rgsw_ct_m, query.rgsw_ct_sm, cfg);

  auto response_ct =
      AnswerLargeRecordOnlinePaper(selected_column_response, state.key_material,
                                   state.packed_hint_blocks, point_ct,
                                   fparm.get_ntt().Raw(), fparm, cfg);

  LargeRecordResponse response;
  response.ct_a = response_ct.a.get_data();
  response.ct_b = response_ct.b.get_data();
  return response;
}

std::vector<uint8_t> RecoverLargeRecordResponse(
    const LargeRecordResponse& response, uint64_t target_index,
    const LargeRecordConfig& cfg, const LargeRecordClientSecrets& secrets) {
  const uint64_t bundle_factor = GetBundleFactor(cfg);
  const uint64_t local_offset = target_index % bundle_factor;

  lrhe::Secret rlwe_sk(
      const_cast<std::vector<int64_t>&>(secrets.sk_vec_signed), cfg.qmod,
      false);

  RlweCiphertext ct(cfg.degree, cfg.qmod, true);
  ct.a = Poly(response.ct_a, cfg.qmod, true);
  ct.b = Poly(response.ct_b, cfg.qmod, true);

  Plaintext dec_pt(cfg.degree, cfg.pmod);
  RlweBfvDecrypt(ct, dec_pt, rlwe_sk);

  Record recovered = ExtractRecordFromBundle(dec_pt.pt.get_data(), local_offset);
  return std::vector<uint8_t>(recovered.begin(), recovered.end());
}

}  // namespace psi::inspire::internal
