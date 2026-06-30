#pragma once

#include <cstdint>
#include <vector>

#include "psi/algorithm/inspire/types.h"

namespace psi::inspire::internal {

struct LargeRecordClientSecrets {
  std::vector<int64_t> sk_vec_signed;
  std::vector<uint64_t> sk_vec;
};

LargeRecordPrecomputedState PrepareLargeRecordState(
    const LargeRecordConfig& cfg);

LargeRecordQuery GenerateLargeRecordQuery(
    uint64_t target_index, const LargeRecordConfig& cfg,
    LargeRecordClientSecrets& secrets);

LargeRecordResponse ProcessLargeRecordQuery(
    const LargeRecordQuery& query,
    const LargeRecordPrecomputedState& state,
    const LargeRecordConfig& cfg);

std::vector<uint8_t> RecoverLargeRecordResponse(
    const LargeRecordResponse& response, uint64_t target_index,
    const LargeRecordConfig& cfg,
    const LargeRecordClientSecrets& secrets);

}  // namespace psi::inspire::internal
