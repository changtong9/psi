#pragma once

#include <cstdint>
#include <random>
#include <vector>

#include "psi/algorithm/ypir/legacy/hexl.h"
#include "psi/algorithm/ypir/legacy/ypir_params.h"

namespace psi::inspire::internal {

using psi::ypir::ypir_internal::AESCTR_PRNG;
using psi::ypir::ypir_internal::EltwiseAddMod;
using psi::ypir::ypir_internal::EltwiseFMAMod;
using psi::ypir::ypir_internal::EltwiseMultMod;
using psi::ypir::ypir_internal::EltwiseSubMod;
using psi::ypir::ypir_internal::FheParams;
using psi::ypir::ypir_internal::PirParams;
using psi::ypir::ypir_internal::SampleGauss;
using psi::ypir::ypir_internal::Secret;
using psi::ypir::ypir_internal::YpirHexlNtt;

uint64_t GetBase(uint64_t b, uint64_t z, uint64_t ti);

uint64_t GetInspiringGenerator(uint64_t gamma, uint64_t degree);

void PrecomputeAutomaps(uint32_t length,
                        const std::vector<uint64_t>& idx,
                        std::vector<std::vector<uint32_t>>& automaps);

void ApplyAutoCoefForm(std::vector<uint64_t>& result,
                       const std::vector<uint64_t>& input, int32_t index,
                       uint64_t modulus);

void ApplyAutoNttFormGeneral(const std::vector<uint64_t>& vec,
                             std::vector<uint64_t>& result,
                             const std::vector<uint32_t>& automap);

void RotateKskRows(const std::vector<std::vector<uint64_t>>& in,
                   const std::vector<uint32_t>& automap,
                   std::vector<std::vector<uint64_t>>& out);

void BuildNttMonomials(uint64_t degree, uint64_t count, YpirHexlNtt& ntt,
                       std::vector<std::vector<uint64_t>>& monomials_ntt);

void ApproximateGadgetDecomp(const std::vector<uint64_t>& poly,
                             std::vector<std::vector<uint64_t>>& mat,
                             uint64_t b, uint64_t z, uint64_t t);

void KeyswitchPreprocess(
    const std::vector<std::vector<uint64_t>>& ksk_a,
    std::vector<uint64_t>& a_in, std::vector<uint64_t>& a_out,
    std::vector<std::vector<std::vector<uint64_t>>>& decomp_buf,
    const FheParams& fparm);

void KeyswitchOnline(const std::vector<std::vector<uint64_t>>& ksk_b,
                     std::vector<uint64_t>& b_in,
                     std::vector<uint64_t>& b_out,
                     const std::vector<std::vector<uint64_t>>& decomp_buf,
                     const FheParams& fparm);

void LweEncrypt(Secret& sk, const std::vector<uint64_t>& a, uint64_t message,
                uint64_t& b, uint64_t pmod, double sig);

void RlweEncode(Secret& sk, const std::vector<uint64_t>& a,
                std::vector<uint64_t>& message, std::vector<uint64_t>& b,
                const FheParams& fparm);

void GadgetEncrypt(Secret& sk,
                   const std::vector<std::vector<uint64_t>>& a,
                   std::vector<uint64_t>& message,
                   std::vector<std::vector<uint64_t>>& b,
                   const FheParams& fparm);

void SimplepirQuery(uint64_t idx, Secret& sk, std::vector<uint64_t>& qu,
                    const FheParams& fparm, const PirParams& pparm);

}  // namespace psi::inspire::internal
