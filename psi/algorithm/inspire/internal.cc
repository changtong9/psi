#include "psi/algorithm/inspire/internal.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include "yacl/base/exception.h"

#include "psi/algorithm/ypir/legacy/client.h"
#include "psi/algorithm/ypir/legacy/hexl.h"
#include "psi/algorithm/ypir/legacy/server.h"
#include "psi/algorithm/ypir/legacy/ypir_params.h"
#include "psi/algorithm/ypir/legacy/ypir_util.h"
#include "psi/algorithm/ypir/util.h"

namespace psi::inspire::internal {
namespace {

using psi::ypir::YpirMode;
using psi::ypir::ypir_internal::AESCTR_PRNG;
using psi::ypir::ypir_internal::Cdks21Lwe2RlweInplace;
using psi::ypir::ypir_internal::EltwiseAddMod;
using psi::ypir::ypir_internal::EltwiseFMAMod;
using psi::ypir::ypir_internal::EltwiseMultMod;
using psi::ypir::ypir_internal::EltwiseSubMod;
using psi::ypir::ypir_internal::FheParams;
using psi::ypir::ypir_internal::GetLog2;
using psi::ypir::ypir_internal::IsPowerOfTwo;
using psi::ypir::ypir_internal::kSecondDimensionSeed;
using psi::ypir::ypir_internal::kThirdDimensionSeed;
using psi::ypir::ypir_internal::MatrixMultiplicationFlat;
using psi::ypir::ypir_internal::MatrixMultiplicationFlatU16;
using psi::ypir::ypir_internal::MatrixRowDecompose;
using psi::ypir::ypir_internal::MatrixTranspose;
using psi::ypir::ypir_internal::MatrixVectorFirstDimension;
using psi::ypir::ypir_internal::MatrixVectorMultiplicationU16;
using psi::ypir::ypir_internal::MatVecU8U32Mod2p32;
using psi::ypir::ypir_internal::ModInverse;
using psi::ypir::ypir_internal::PirParams;
using psi::ypir::ypir_internal::PrecomputeAutomap;
using psi::ypir::ypir_internal::PseudorandomMatrixGenerate;
using psi::ypir::ypir_internal::SampleGauss;
using psi::ypir::ypir_internal::Secret;
using psi::ypir::ypir_internal::VectorColDecompose;
using psi::ypir::ypir_internal::YpirHexlNtt;

uint64_t MakeSeedMaterial() {
  const uint64_t time_seed = static_cast<uint64_t>(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());
  std::random_device rd;
  return time_seed ^ (static_cast<uint64_t>(rd()) << 32) ^ rd();
}

std::mt19937_64 MakePrng() {
  return std::mt19937_64(MakeSeedMaterial());
}

uint64_t GetBase(uint64_t b, uint64_t z, uint64_t ti) {
  return 1ULL << (b + ti * z);
}

uint64_t GetInspiringGenerator(uint64_t gamma, uint64_t degree) {
  if (gamma == degree) {
    return 5;
  }
  return ((degree << 1) / gamma) + 1;
}

void PrecomputeAutomaps(uint32_t length, const std::vector<uint64_t>& idx,
                        std::vector<std::vector<uint32_t>>& automaps) {
  automaps.assign(idx.size(), std::vector<uint32_t>(length, 0));
  for (uint64_t i = 0; i < idx.size(); ++i) {
    PrecomputeAutomap(length, static_cast<uint32_t>(idx[i]), automaps[i]);
  }
}

void ApplyAutoNttFormGeneral(const std::vector<uint64_t>& vec,
                             std::vector<uint64_t>& result,
                             const std::vector<uint32_t>& automap) {
  const uint64_t length = vec.size();
  if (result.size() != length) {
    result.resize(length);
  }
  if (vec.data() == result.data()) {
    std::vector<uint64_t> scratch(length, 0);
    for (uint64_t i = 0; i < length; ++i) {
      scratch[i] = vec[automap[i]];
    }
    result.swap(scratch);
    return;
  }
  for (uint64_t i = 0; i < length; ++i) {
    result[i] = vec[automap[i]];
  }
}

void RotateKskRows(const std::vector<std::vector<uint64_t>>& in,
                   const std::vector<uint32_t>& automap,
                   std::vector<std::vector<uint64_t>>& out) {
  const uint64_t rows = in.size();
  const uint64_t degree = rows == 0 ? 0 : in[0].size();
  out.assign(rows, std::vector<uint64_t>(degree, 0));
  for (uint64_t i = 0; i < rows; ++i) {
    ApplyAutoNttFormGeneral(in[i], out[i], automap);
  }
}

void BuildNttMonomials(uint64_t degree, uint64_t count, YpirHexlNtt& ntt,
                       std::vector<std::vector<uint64_t>>& monomials_ntt) {
  monomials_ntt.assign(count, std::vector<uint64_t>(degree, 0));
  for (uint64_t i = 0; i < count; ++i) {
    monomials_ntt[i][i] = 1;
    ntt.Forward(monomials_ntt[i].data(), degree);
  }
}

void ApplyAutoCoefForm(std::vector<uint64_t>& result,
                       const std::vector<uint64_t>& input, int32_t index,
                       uint64_t modulus) {
  const uint64_t length = input.size();
  result.assign(length, 0);
  for (uint64_t i = 0; i < length; ++i) {
    uint64_t destination = (i * index) % (2 * length);
    if (destination >= length) {
      result[destination - length] = (modulus - input[i]) % modulus;
    } else {
      result[destination] = input[i];
    }
  }
}

void ApproximateGadgetDecomp(const std::vector<uint64_t>& poly,
                             std::vector<std::vector<uint64_t>>& mat,
                             uint64_t b, uint64_t z, uint64_t t) {
  const uint64_t degree = poly.size();
  const uint64_t mask = (1ULL << z) - 1;
  mat.assign(t, std::vector<uint64_t>(degree, 0));
  for (uint64_t i = 0; i < degree; ++i) {
    uint64_t val = poly[i] >> b;
    for (uint64_t j = 0; j < t; ++j) {
      mat[j][i] = val & mask;
      val >>= z;
    }
  }
}

void KeyswitchPreprocess(
    const std::vector<std::vector<uint64_t>>& ksk_a,
    std::vector<uint64_t>& a_in, std::vector<uint64_t>& a_out,
    std::vector<std::vector<std::vector<uint64_t>>>& decomp_buf,
    const FheParams& fparm) {
  const uint64_t length = fparm.get_poly_degree();
  const uint64_t modulus = fparm.get_rlwe_cmod();
  const uint64_t t = fparm.get_t_auto();
  YpirHexlNtt& ntt = fparm.get_ntt();

  std::vector<std::vector<uint64_t>> decomp_a;
  ntt.Inverse(a_in.data(), length);
  ApproximateGadgetDecomp(a_in, decomp_a, fparm.get_b_auto(),
                          fparm.get_z_auto(), t);
  for (uint64_t i = 0; i < t; ++i) {
    ntt.Forward(decomp_a[i].data(), length);
  }

  decomp_buf.push_back(std::move(decomp_a));
  a_out.assign(length, 0);
  std::vector<uint64_t> tmp(length, 0);
  for (uint64_t i = 0; i < t; ++i) {
    EltwiseMultMod(tmp.data(), ksk_a[i].data(), decomp_buf.back()[i].data(),
                   length, modulus);
    EltwiseSubMod(a_out.data(), a_out.data(), tmp.data(), length, modulus);
  }
}

void KeyswitchOnline(const std::vector<std::vector<uint64_t>>& ksk_b,
                     std::vector<uint64_t>& b_in, std::vector<uint64_t>& b_out,
                     const std::vector<std::vector<uint64_t>>& decomp_buf,
                     const FheParams& fparm) {
  const uint64_t length = fparm.get_poly_degree();
  const uint64_t modulus = fparm.get_rlwe_cmod();
  const uint64_t t = fparm.get_t_auto();
  b_out = b_in;
  std::vector<uint64_t> tmp(length, 0);
  for (uint64_t i = 0; i < t; ++i) {
    EltwiseMultMod(tmp.data(), ksk_b[i].data(), decomp_buf[i].data(), length,
                   modulus);
    EltwiseSubMod(b_out.data(), b_out.data(), tmp.data(), length, modulus);
  }
}

void LweEncrypt(Secret& sk, const std::vector<uint64_t>& a, uint64_t message,
                uint64_t& b, uint64_t pmod, double sig) {
  const uint64_t lwe_dimension = sk.get_len();
  const uint64_t cmod = sk.get_mod();
  const long double delta =
      static_cast<long double>(cmod) / static_cast<long double>(pmod);

  auto rng = MakePrng();
  const uint64_t e = SampleGauss(sig, cmod, rng);
  std::vector<uint64_t> tmp(lwe_dimension, 0);
  EltwiseMultMod(tmp.data(), a.data(), sk.data.data(), lwe_dimension, cmod);
  b = 0;
  for (uint64_t i = 0; i < lwe_dimension; ++i) {
    b = (b + tmp[i]) % cmod;
  }
  b = (b + e) % cmod;
  b = (b + static_cast<uint64_t>(message * delta)) % cmod;
}

void RlweEncode(Secret& sk, const std::vector<uint64_t>& a,
                std::vector<uint64_t>& message, std::vector<uint64_t>& b,
                const FheParams& fparm) {
  const uint64_t poly_degree = fparm.get_poly_degree();
  const uint64_t rlwe_cmod = fparm.get_rlwe_cmod();
  YpirHexlNtt& ntt = fparm.get_ntt();
  if (!sk.get_ntt_form()) {
    ntt.Forward(sk.data.data(), poly_degree);
    sk.switch_ntt_format();
  }

  std::vector<uint64_t> err(poly_degree, 0);
  auto rng = MakePrng();
  SampleGauss(err, fparm.get_sig_ring(), rlwe_cmod, rng);
  EltwiseFMAMod(err.data(), message.data(), 1, err.data(), poly_degree,
                rlwe_cmod);
  ntt.Forward(err.data(), poly_degree);
  EltwiseMultMod(b.data(), sk.data.data(), a.data(), poly_degree, rlwe_cmod);
  EltwiseAddMod(b.data(), err.data(), b.data(), poly_degree, rlwe_cmod);
}

void GadgetEncrypt(Secret& sk, const std::vector<std::vector<uint64_t>>& a,
                   std::vector<uint64_t>& message,
                   std::vector<std::vector<uint64_t>>& b,
                   const FheParams& fparm) {
  const uint64_t t_auto = fparm.get_t_auto();
  const uint64_t poly_degree = fparm.get_poly_degree();
  const uint64_t rlwe_cmod = fparm.get_rlwe_cmod();
  std::vector<uint64_t> message_tmp(poly_degree, 0);
  for (uint64_t i = 0; i < t_auto; ++i) {
    const uint64_t base =
        GetBase(fparm.get_b_auto(), fparm.get_z_auto(), i) % rlwe_cmod;
    EltwiseFMAMod(message_tmp.data(), message.data(), base, nullptr,
                  poly_degree, rlwe_cmod);
    RlweEncode(sk, a[i], message_tmp, b[i], fparm);
  }
}

void GenerateInspiringPartialKey(
    Secret& sk, uint64_t gamma, std::vector<std::vector<uint64_t>>& ksk_b,
    AESCTR_PRNG& prng, const FheParams& fparm) {
  const uint64_t degree = fparm.get_poly_degree();
  const uint64_t cmod = fparm.get_rlwe_cmod();
  const uint64_t t_auto = fparm.get_t_auto();
  const uint64_t generator = GetInspiringGenerator(gamma, degree);
  YpirHexlNtt& ntt = fparm.get_ntt();

  if (sk.get_ntt_form()) {
    ntt.Inverse(sk.data.data(), degree);
    sk.switch_ntt_format();
  }

  std::vector<std::vector<uint64_t>> ksk_a(
      t_auto, std::vector<uint64_t>(degree, 0));
  ksk_b.assign(t_auto, std::vector<uint64_t>(degree, 0));
  prng.refresh(kThirdDimensionSeed);
  PseudorandomMatrixGenerate(ksk_a, cmod, prng);

  std::vector<uint64_t> newkey(degree, 0);
  const std::vector<uint64_t> copykey = sk.data;
  ApplyAutoCoefForm(newkey, copykey, static_cast<int32_t>(generator), cmod);
  GadgetEncrypt(sk, ksk_a, newkey, ksk_b, fparm);
}

void DoublepirQuery(uint64_t c_idx, uint64_t r_idx, Secret& simple_sk,
                    Secret& double_sk, std::vector<uint64_t>& qu0,
                    std::vector<uint64_t>& qu1, AESCTR_PRNG& prng,
                    const FheParams& fparm, const PirParams& pparm) {
  const uint64_t cols = pparm.get_col();
  qu0.resize(cols);
  const auto& matrix0 = fparm.get_persudo_matrix_simplepir();
  for (uint64_t i = 0; i < cols; ++i) {
    LweEncrypt(simple_sk, matrix0[i], i == c_idx ? 1 : 0, qu0[i],
               fparm.get_lwe_pmod(), fparm.get_sig());
  }

  const uint64_t rows = pparm.get_row();
  const uint64_t degree = fparm.get_poly_degree();
  qu1.resize(rows);
  std::vector<std::vector<uint64_t>> matrix1(rows,
                                             std::vector<uint64_t>(degree, 0));
  prng.refresh(kSecondDimensionSeed);
  PseudorandomMatrixGenerate(matrix1, fparm.get_rlwe_cmod(), prng);
  for (uint64_t i = 0; i < rows; ++i) {
    LweEncrypt(double_sk, matrix1[i], i == r_idx ? 1 : 0, qu1[i],
               fparm.get_rlwe_pmod(), fparm.get_sig_ring());
  }
}

std::vector<uint64_t> InspiringPartialPreprocess(
    const std::vector<std::vector<uint64_t>>& a_ct_tilde, uint64_t gamma,
    const std::vector<std::vector<uint64_t>>& ksk_a,
    std::vector<std::vector<std::vector<uint64_t>>>& decomp_buf,
    const FheParams& fparm) {
  const uint64_t degree = fparm.get_poly_degree();
  const uint64_t modulus = fparm.get_rlwe_cmod();
  YpirHexlNtt& ntt = fparm.get_ntt();
  if (gamma == 1) {
    return a_ct_tilde[0];
  }
  YACL_ENFORCE_LE(gamma, degree >> 1);
  YACL_ENFORCE(IsPowerOfTwo(gamma));

  const uint64_t generator = GetInspiringGenerator(gamma, degree);
  std::vector<uint64_t> generator_powers(gamma, 1);
  for (uint64_t i = 1; i < gamma; ++i) {
    generator_powers[i] =
        static_cast<uint64_t>((static_cast<unsigned __int128>(
                                   generator_powers[i - 1]) *
                               generator) %
                              (degree << 1));
  }

  std::vector<std::vector<uint32_t>> automaps;
  PrecomputeAutomaps(static_cast<uint32_t>(degree), generator_powers, automaps);

  std::vector<std::vector<uint64_t>> monomials_ntt;
  BuildNttMonomials(degree, gamma, ntt, monomials_ntt);

  const uint64_t mod_inv = ModInverse(static_cast<int64_t>(gamma),
                                      static_cast<int64_t>(modulus));
  std::vector<std::vector<uint64_t>> scaled_a = a_ct_tilde;
  for (uint64_t i = 0; i < gamma; ++i) {
    EltwiseFMAMod(scaled_a[i].data(), scaled_a[i].data(), mod_inv, nullptr,
                  degree, modulus);
  }

  std::vector<std::vector<uint64_t>> aggregated(
      gamma, std::vector<uint64_t>(degree, 0));
  std::vector<uint64_t> rotated(degree, 0);
  std::vector<uint64_t> shifted(degree, 0);
  for (uint64_t k = 0; k < gamma; ++k) {
    for (uint64_t j = 0; j < gamma; ++j) {
      ApplyAutoNttFormGeneral(scaled_a[k], rotated, automaps[j]);
      EltwiseMultMod(shifted.data(), rotated.data(), monomials_ntt[k].data(),
                     degree, modulus);
      EltwiseAddMod(aggregated[j].data(), aggregated[j].data(), shifted.data(),
                    degree, modulus);
    }
  }

  std::vector<std::vector<uint64_t>> current = std::move(aggregated);
  std::vector<std::vector<uint64_t>> rotated_ksk_a;
  for (uint64_t level = gamma - 1; level > 0; --level) {
    RotateKskRows(ksk_a, automaps[level - 1], rotated_ksk_a);
    std::vector<uint64_t> source = current[level];
    std::vector<uint64_t> switched(degree, 0);
    KeyswitchPreprocess(rotated_ksk_a, source, switched, decomp_buf, fparm);
    EltwiseAddMod(current[level - 1].data(), current[level - 1].data(),
                  switched.data(), degree, modulus);
  }
  return current[0];
}

std::vector<uint64_t> InspiringPartialOnline(
    const std::vector<uint64_t>& b_values, uint64_t gamma,
    const std::vector<std::vector<uint64_t>>& ksk_b,
    const std::vector<std::vector<std::vector<uint64_t>>>& decomp_buf,
    const FheParams& fparm) {
  const uint64_t degree = fparm.get_poly_degree();
  const uint64_t modulus = fparm.get_rlwe_cmod();
  YpirHexlNtt& ntt = fparm.get_ntt();

  std::vector<uint64_t> packed_b(degree, 0);
  for (uint64_t i = 0; i < gamma; ++i) {
    packed_b[i] = b_values[i] % modulus;
  }
  ntt.Forward(packed_b.data(), degree);
  if (gamma == 1) {
    return packed_b;
  }
  YACL_ENFORCE_LE(gamma, degree >> 1);
  YACL_ENFORCE(IsPowerOfTwo(gamma));
  YACL_ENFORCE_EQ(decomp_buf.size(), gamma - 1);

  const uint64_t generator = GetInspiringGenerator(gamma, degree);
  std::vector<uint64_t> generator_powers(gamma, 1);
  for (uint64_t i = 1; i < gamma; ++i) {
    generator_powers[i] =
        static_cast<uint64_t>((static_cast<unsigned __int128>(
                                   generator_powers[i - 1]) *
                               generator) %
                              (degree << 1));
  }
  std::vector<std::vector<uint32_t>> automaps;
  PrecomputeAutomaps(static_cast<uint32_t>(degree), generator_powers, automaps);

  std::vector<std::vector<uint64_t>> rotated_ksk_b;
  for (uint64_t level = gamma - 1, step = 0; level > 0; --level, ++step) {
    RotateKskRows(ksk_b, automaps[level - 1], rotated_ksk_b);
    std::vector<uint64_t> switched(degree, 0);
    KeyswitchOnline(rotated_ksk_b, packed_b, switched, decomp_buf[step],
                    fparm);
    packed_b.swap(switched);
  }
  return packed_b;
}

void DoublepirHintGenerate(
    std::vector<std::vector<uint64_t>>& db,
    std::vector<std::vector<uint64_t>>& double_server_hint,
    std::vector<std::vector<uint64_t>>& double_client_hint,
    const FheParams& fparm) {
  const uint64_t lwe_dimension = fparm.get_lwe_dimension();
  const uint64_t poly_degree = fparm.get_poly_degree();
  const uint64_t rlwe_cmod = fparm.get_rlwe_cmod();

  const uint64_t* matrix_flat = fparm.get_persudo_matrix_simplepir_flat();
  std::vector<std::vector<uint64_t>> simple_hint(
      db.size(), std::vector<uint64_t>(lwe_dimension, 0));
  MatrixMultiplicationFlat(simple_hint, db, matrix_flat, lwe_dimension,
                           fparm.get_lwe_cmod());

  std::vector<std::vector<uint64_t>> matrix_decomp;
  MatrixRowDecompose(simple_hint, matrix_decomp, fparm.get_b_decomp(),
                     fparm.get_z_decomp(), fparm.get_t_decomp());
  MatrixTranspose(matrix_decomp, double_server_hint);

  const uint64_t* matrix_d2_flat = fparm.get_persudo_matrix_doublepir_flat();
  MatrixMultiplicationFlat(double_client_hint, double_server_hint,
                           matrix_d2_flat, poly_degree, rlwe_cmod);
}

void DoublepirAnswer(const uint8_t* db, uint32_t* qu0,
                     std::vector<uint64_t>& qu1,
                     const std::vector<std::vector<uint64_t>>& server_hint,
                     std::vector<uint64_t>& h2,
                     std::vector<std::vector<uint64_t>>& h3,
                     std::vector<uint64_t>& h4, const FheParams& fparm,
                     const PirParams& pparm) {
  std::vector<uint32_t> simple_res(pparm.get_row(), 0);
  MatVecU8U32Mod2p32(db, qu0, simple_res.data(), pparm.get_row(),
                     pparm.get_col());

  std::vector<std::vector<uint16_t>> trans_simple_res_decomp;
  VectorColDecompose(simple_res, trans_simple_res_decomp,
                     fparm.get_b_decomp(), fparm.get_z_decomp(),
                     fparm.get_t_decomp());

  MatrixVectorFirstDimension(h2, server_hint, qu1, fparm.get_rlwe_cmod());
  MatrixMultiplicationFlatU16(
      h3, trans_simple_res_decomp,
      fparm.get_persudo_matrix_doublepir_flat(), fparm.get_poly_degree(),
      fparm.get_rlwe_cmod());
  MatrixVectorMultiplicationU16(h4, trans_simple_res_decomp, qu1,
                                fparm.get_rlwe_cmod());
}

std::vector<std::vector<uint64_t>> ExpandDatabase(const std::vector<uint8_t>& db,
                                                  uint64_t rows,
                                                  uint64_t cols) {
  YACL_ENFORCE_EQ(db.size(), rows * cols);
  std::vector<std::vector<uint64_t>> out(rows, std::vector<uint64_t>(cols, 0));
  for (uint64_t row = 0; row < rows; ++row) {
    for (uint64_t col = 0; col < cols; ++col) {
      out[row][col] = db[row * cols + col];
    }
  }
  return out;
}

}  // namespace

InspireQuery GenerateQuery(uint64_t raw_idx,
                           const psi::ypir::YpirParameters& params,
                           ClientSecrets& secrets,
                           const psi::ypir::internal::ypir::Context& context) {
  YACL_ENFORCE(params.mode == YpirMode::kDoublepir);
  YACL_ENFORCE_LT(raw_idx, params.NumItems());

  const uint64_t row_idx = raw_idx / params.db_cols;
  const uint64_t col_idx = raw_idx % params.db_cols;
  secrets.simple_secret =
      Secret(context.fhe_params->get_lwe_dimension(),
             context.fhe_params->get_lwe_cmod());
  secrets.double_secret =
      Secret(context.fhe_params->get_poly_degree(),
             context.fhe_params->get_rlwe_cmod());
  secrets.initialized = true;

  InspireQuery query;
  DoublepirQuery(col_idx, row_idx, secrets.simple_secret, secrets.double_secret,
                 query.qu0, query.qu1, *context.prng, *context.fhe_params,
                 *context.pir_params);
  GenerateInspiringPartialKey(
      secrets.double_secret,
      context.fhe_params->get_lwe_dimension() *
          context.fhe_params->get_t_decomp(),
      query.ksk_b, *context.prng, *context.fhe_params);
  return query;
}

std::vector<uint8_t> RecoverResponse(
    const InspireResponse& response, const psi::ypir::YpirParameters& params,
    const ClientSecrets& secrets,
    const psi::ypir::internal::ypir::Context& context) {
  YACL_ENFORCE(params.mode == YpirMode::kDoublepir);
  YACL_ENFORCE(secrets.initialized,
               "GenerateQuery must be called before decode");
  auto simple_secret = secrets.simple_secret;
  auto double_secret = secrets.double_secret;
  auto result = response.doublepir_response;
  uint64_t message = 0;
  psi::ypir::ypir_internal::YpirRecover(
      simple_secret, double_secret, result, message, *context.fhe_params,
      *context.pir_params);
  return psi::ypir::EncodeIntegerValue(message, params.value_bytes);
}

InspirePrecomputedState PrepareOfflineState(
    const std::vector<uint8_t>& db, const psi::ypir::YpirParameters& params,
    const psi::ypir::internal::ypir::Context& context) {
  YACL_ENFORCE(params.mode == YpirMode::kDoublepir);
  auto db_matrix = ExpandDatabase(db, params.db_rows, params.db_cols);

  InspirePrecomputedState state;
  std::vector<std::vector<uint64_t>> double_server_hint;
  const uint64_t pack_num =
      context.fhe_params->get_lwe_dimension() *
      context.fhe_params->get_t_decomp();
  std::vector<std::vector<uint64_t>> double_client_hint(
      pack_num,
      std::vector<uint64_t>(context.fhe_params->get_poly_degree(), 0));
  DoublepirHintGenerate(db_matrix, double_server_hint, double_client_hint,
                        *context.fhe_params);

  std::vector<std::vector<uint64_t>> ksk_a(
      context.fhe_params->get_t_auto(),
      std::vector<uint64_t>(context.fhe_params->get_poly_degree(), 0));
  context.prng->refresh(kThirdDimensionSeed);
  PseudorandomMatrixGenerate(ksk_a, context.fhe_params->get_rlwe_cmod(),
                             *context.prng);

  YpirHexlNtt& ntt = context.fhe_params->get_ntt();
  for (uint64_t i = 0; i < pack_num; ++i) {
    Cdks21Lwe2RlweInplace(double_client_hint[i].data(),
                          context.fhe_params->get_poly_degree(),
                          context.fhe_params->get_rlwe_cmod(), ntt);
  }

  state.decomp_buf.clear();
  state.hint_0 = InspiringPartialPreprocess(
      double_client_hint, pack_num, ksk_a, state.decomp_buf,
      *context.fhe_params);
  MatrixTranspose(double_server_hint, state.server_hint);
  return state;
}

InspireResponse ProcessQuery(
    const std::vector<uint8_t>& db, const InspireQuery& query,
    const InspirePrecomputedState& state,
    const psi::ypir::YpirParameters& params,
    const psi::ypir::internal::ypir::Context& context) {
  YACL_ENFORCE(params.mode == YpirMode::kDoublepir);
  YACL_ENFORCE_EQ(context.fhe_params->get_t_decomp(), 1ULL,
                  "Inspire currently requires t_decomp == 1");

  std::vector<uint32_t> qu0(query.qu0.size(), 0);
  for (size_t i = 0; i < query.qu0.size(); ++i) {
    qu0[i] = static_cast<uint32_t>(query.qu0[i]);
  }

  auto qu1 = query.qu1;
  std::vector<std::vector<uint64_t>> result;
  result.push_back(state.hint_0);

  const uint64_t lwe_dimension = context.fhe_params->get_lwe_dimension();
  const uint64_t poly_degree = context.fhe_params->get_poly_degree();
  const uint64_t pack_num =
      lwe_dimension * context.fhe_params->get_t_decomp();
  std::vector<std::vector<uint64_t>> h3(
      context.fhe_params->get_t_decomp(),
      std::vector<uint64_t>(poly_degree, 0));
  std::vector<uint64_t> h2(pack_num, 0);
  std::vector<uint64_t> h4(context.fhe_params->get_t_decomp(), 0);
  DoublepirAnswer(db.data(), qu0.data(), qu1, state.server_hint, h2, h3, h4,
                  *context.fhe_params, *context.pir_params);

  result.push_back(InspiringPartialOnline(h2, pack_num, query.ksk_b,
                                          state.decomp_buf,
                                          *context.fhe_params));

  std::vector<uint64_t> pack_b(poly_degree, 0);
  pack_b[0] = h4[0];
  context.fhe_params->get_ntt().Forward(pack_b.data(), poly_degree);
  Cdks21Lwe2RlweInplace(h3[0].data(), poly_degree,
                        context.fhe_params->get_rlwe_cmod(),
                        context.fhe_params->get_ntt());
  result.push_back(h3[0]);
  result.push_back(pack_b);

  InspireResponse response;
  response.doublepir_response = std::move(result);
  return response;
}

}  // namespace psi::inspire::internal
