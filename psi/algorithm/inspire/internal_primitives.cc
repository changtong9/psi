#include "psi/algorithm/inspire/internal_primitives.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include "yacl/base/exception.h"

#include "psi/algorithm/ypir/legacy/ypir_util.h"

namespace psi::inspire::internal {

namespace {

using psi::ypir::ypir_internal::kThirdDimensionSeed;
using psi::ypir::ypir_internal::PseudorandomMatrixGenerate;

uint64_t MakeSeedMaterial() {
  const uint64_t time_seed = static_cast<uint64_t>(
      std::chrono::high_resolution_clock::now().time_since_epoch().count());
  std::random_device rd;
  return time_seed ^ (static_cast<uint64_t>(rd()) << 32) ^ rd();
}

std::mt19937_64 MakePrng() { return std::mt19937_64(MakeSeedMaterial()); }

}  // namespace

uint64_t GetBase(uint64_t b, uint64_t z, uint64_t ti) {
  return 1ULL << (b + ti * z);
}

uint64_t GetInspiringGenerator(uint64_t gamma, uint64_t degree) {
  if (gamma == degree) {
    return 5;
  }
  return ((degree << 1) / gamma) + 1;
}

void PrecomputeAutomaps(uint32_t length,
                        const std::vector<uint64_t>& idx,
                        std::vector<std::vector<uint32_t>>& automaps) {
  automaps.assign(idx.size(), std::vector<uint32_t>(length, 0));
  for (uint64_t i = 0; i < idx.size(); ++i) {
    psi::ypir::ypir_internal::PrecomputeAutomap(
        length, static_cast<uint32_t>(idx[i]), automaps[i]);
  }
}

void ApplyAutoCoefForm(std::vector<uint64_t>& result,
                       const std::vector<uint64_t>& input, int32_t index,
                       uint64_t modulus) {
  const uint64_t length = input.size();
  result.assign(length, 0);
  for (uint64_t i = 0; i < length; ++i) {
    uint64_t destination = (i * static_cast<uint64_t>(index)) % (2 * length);
    if (destination >= length) {
      result[destination - length] = (modulus - input[i]) % modulus;
    } else {
      result[destination] = input[i];
    }
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
                     std::vector<uint64_t>& b_in,
                     std::vector<uint64_t>& b_out,
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

void GadgetEncrypt(Secret& sk,
                   const std::vector<std::vector<uint64_t>>& a,
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

void SimplepirQuery(uint64_t idx, Secret& sk, std::vector<uint64_t>& qu,
                    const FheParams& fparm, const PirParams& pparm) {
  const uint64_t length = pparm.get_col();
  const uint64_t lwe_pmod = fparm.get_lwe_pmod();
  const float sig = fparm.get_sig();
  qu.resize(length, 0);

  const auto& matrix = fparm.get_persudo_matrix_simplepir();
  for (uint64_t i = 0; i < length; ++i) {
    LweEncrypt(sk, matrix[i], i == idx ? 1 : 0, qu[i], lwe_pmod, sig);
  }
}

}  // namespace psi::inspire::internal
