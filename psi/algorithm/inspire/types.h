#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace psi::inspire {

struct InspirePrecomputedState {
  std::vector<uint64_t> hint_0;
  std::vector<std::vector<uint64_t>> server_hint;
  std::vector<std::vector<std::vector<uint64_t>>> decomp_buf;
};

struct InspireQuery {
  std::vector<uint64_t> qu0;
  std::vector<uint64_t> qu1;
  std::vector<std::vector<uint64_t>> ksk_b;
};

struct InspireResponse {
  std::vector<std::vector<uint64_t>> doublepir_response;
};

struct LargeRecordConfig {
  static constexpr uint64_t kCrtQ1 = 268369921ULL;
  static constexpr uint64_t kCrtQ2 = 249561089ULL;
  static constexpr uint64_t kCrtMod = kCrtQ1 * kCrtQ2;
  static constexpr size_t kRecordBytes = 256;

  uint64_t total_records = 1ULL << 20;
  uint64_t folding = 8;
  uint64_t degree = 2048;
  uint64_t packed_coeffs = 2048;
  uint64_t qmod = kCrtMod;
  uint64_t pmod = 257;
  float sigma = 0.41f;
  uint64_t rgsw_b = 0;
  uint64_t rgsw_z = 8;
  uint64_t rgsw_t = 7;
  uint64_t lwe_dimension = 2048;
  uint64_t lwe_cmod = kCrtMod;
  uint64_t lwe_pmod = 257;
};

struct PaperPackKeyMaterial {
  std::vector<std::vector<std::vector<uint64_t>>> kg_a_powers;
  std::vector<std::vector<std::vector<uint64_t>>> kg_b_powers;
  std::vector<std::vector<std::vector<uint64_t>>> gh_kg_a_powers;
  std::vector<std::vector<std::vector<uint64_t>>> gh_kg_b_powers;
  std::vector<std::vector<uint64_t>> kh_a;
  std::vector<std::vector<uint64_t>> kh_b;
  std::vector<std::vector<uint32_t>> g_automaps;
  std::vector<std::vector<uint32_t>> gh_automaps;
  std::vector<std::vector<uint64_t>> monomials_ntt;
};

struct PaperPackedHintBlock {
  std::vector<uint64_t> packed_a_ntt;
  std::vector<std::vector<std::vector<uint64_t>>> decomp_half1;
  std::vector<std::vector<std::vector<uint64_t>>> decomp_half2;
  std::vector<std::vector<uint64_t>> decomp_final;
};

struct LargeRecordPrecomputedState {
  std::vector<std::vector<uint16_t>> encoded_db;
  std::vector<std::vector<uint64_t>> hint;
  std::vector<PaperPackedHintBlock> packed_hint_blocks;
  PaperPackKeyMaterial key_material;
};

struct LargeRecordQuery {
  std::vector<uint64_t> column_query;
  std::vector<std::vector<std::vector<uint64_t>>> rgsw_ct_m;
  std::vector<std::vector<std::vector<uint64_t>>> rgsw_ct_sm;
};

struct LargeRecordResponse {
  std::vector<uint64_t> ct_a;
  std::vector<uint64_t> ct_b;
};

}  // namespace psi::inspire
