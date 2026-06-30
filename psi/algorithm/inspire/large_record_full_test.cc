#include <cstdio>
#include <cstdint>
#include <vector>

#include "psi/algorithm/inspire/large_record.h"
#include "psi/algorithm/inspire/types.h"

int main() {
  psi::inspire::LargeRecordConfig cfg;
  cfg.total_records = 1ULL << 20;
  cfg.folding = 8;
  cfg.degree = 2048;
  cfg.packed_coeffs = 2048;

  const uint64_t target_index = 123456;
  const uint64_t kRecordBytes = psi::inspire::LargeRecordConfig::kRecordBytes;

  std::fprintf(stderr, "=== Large Record Paper - Full 2^20 Test ===\n");
  std::fprintf(stderr, "total_records: %lu\n", cfg.total_records);
  std::fprintf(stderr, "target_index: %lu\n", target_index);

  std::vector<uint8_t> expected(kRecordBytes, 0);
  for (uint64_t b = 0; b < kRecordBytes; ++b) {
    expected[b] = static_cast<uint8_t>((17 * target_index + b) & 0xFFULL);
  }

  auto state = psi::inspire::internal::PrepareLargeRecordState(cfg);
  std::fprintf(stderr, "setup done\n");

  psi::inspire::internal::LargeRecordClientSecrets secrets;
  auto query = psi::inspire::internal::GenerateLargeRecordQuery(
      target_index, cfg, secrets);
  std::fprintf(stderr, "query done\n");

  auto response = psi::inspire::internal::ProcessLargeRecordQuery(
      query, state, cfg);
  std::fprintf(stderr, "answer done\n");

  auto recovered = psi::inspire::internal::RecoverLargeRecordResponse(
      response, target_index, cfg, secrets);
  std::fprintf(stderr, "recover done\n");

  bool ok = recovered.size() == kRecordBytes;
  if (ok) {
    for (uint64_t i = 0; i < kRecordBytes; ++i) {
      if (recovered[i] != expected[i]) {
        ok = false;
        std::fprintf(stderr, "mismatch at byte %lu: got %u expected %u\n", i,
                     recovered[i], expected[i]);
        break;
      }
    }
  }
  std::fprintf(stderr, "full PIR recover check: %s\n", ok ? "ok" : "FAIL");
  return ok ? 0 : 1;
}
