#include <cstdio>
#include <cstdint>
#include <vector>

#include "psi/algorithm/inspire/large_record.h"
#include "psi/algorithm/inspire/types.h"

int main() {
  psi::inspire::LargeRecordConfig cfg;
  cfg.total_records = 1ULL << 16;
  cfg.folding = 8;
  cfg.degree = 2048;
  cfg.packed_coeffs = 2048;
  cfg.qmod = psi::inspire::LargeRecordConfig::kCrtMod;
  cfg.pmod = 257;
  cfg.sigma = 0.41f;
  cfg.rgsw_b = 0;
  cfg.rgsw_z = 8;
  cfg.rgsw_t = 7;
  cfg.lwe_dimension = 2048;
  cfg.lwe_cmod = psi::inspire::LargeRecordConfig::kCrtMod;
  cfg.lwe_pmod = 257;

  const uint64_t target_index = 12345;
  const uint64_t kRecordBytes = psi::inspire::LargeRecordConfig::kRecordBytes;

  std::fprintf(stderr, "large_record_paper smoke test\n");
  std::fprintf(stderr, "total_records: %lu\n", cfg.total_records);
  std::fprintf(stderr, "record_bytes: %lu\n", kRecordBytes);
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
  } else {
    std::fprintf(stderr, "size mismatch: got %zu expected %lu\n",
                 recovered.size(), kRecordBytes);
  }

  std::fprintf(stderr, "full PIR recover check: %s\n", ok ? "ok" : "FAIL");
  return ok ? 0 : 1;
}
