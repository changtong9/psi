#pragma once

#include <vector>

#include "psi/algorithm/inspire/types.h"
#include "psi/algorithm/ypir/params.h"
#include "psi/algorithm/ypir/ypir_internal_params.h"

namespace psi::inspire::internal {

struct ClientSecrets {
  psi::ypir::ypir_internal::Secret simple_secret;
  psi::ypir::ypir_internal::Secret double_secret;
  bool initialized = false;
};

InspireQuery GenerateQuery(uint64_t raw_idx,
                           const psi::ypir::YpirParameters& params,
                           ClientSecrets& secrets,
                           const psi::ypir::internal::ypir::Context& context);

std::vector<uint8_t> RecoverResponse(
    const InspireResponse& response, const psi::ypir::YpirParameters& params,
    const ClientSecrets& secrets,
    const psi::ypir::internal::ypir::Context& context);

InspirePrecomputedState PrepareOfflineState(
    const std::vector<uint8_t>& db, const psi::ypir::YpirParameters& params,
    const psi::ypir::internal::ypir::Context& context);

InspireResponse ProcessQuery(
    const std::vector<uint8_t>& db, const InspireQuery& query,
    const InspirePrecomputedState& state,
    const psi::ypir::YpirParameters& params,
    const psi::ypir::internal::ypir::Context& context);

}  // namespace psi::inspire::internal
