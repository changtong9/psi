#include <cstdio>
#include <vector>

#include "psi/algorithm/inspire/internal.h"
#include "psi/algorithm/inspire/params.h"
#include "psi/algorithm/ypir/ypir_internal_params.h"

int main() {
  auto params = psi::inspire::CreateSmallTestParams();
  auto client_context = psi::ypir::internal::ypir::CreateContext(params);
  auto server_context = psi::ypir::internal::ypir::CreateContext(params);
  psi::inspire::internal::ClientSecrets secrets;

  std::vector<uint8_t> db(params.db_rows * params.db_cols, 0);
  for (uint64_t row = 0; row < params.db_rows; ++row) {
    for (uint64_t col = 0; col < params.db_cols; ++col) {
      db[row * params.db_cols + col] = static_cast<uint8_t>((row + col) % 251);
    }
  }

  const auto state =
      psi::inspire::internal::PrepareOfflineState(db, params, server_context);
  const uint64_t row = 111;
  const uint64_t col = 222;
  const uint64_t raw_idx = row * params.db_cols + col;

  const auto query = psi::inspire::internal::GenerateQuery(
      raw_idx, params, secrets, client_context);
  const auto response = psi::inspire::internal::ProcessQuery(
      db, query, state, params, server_context);
  const auto decoded = psi::inspire::internal::RecoverResponse(
      response, params, secrets, client_context);

  const uint8_t expected = static_cast<uint8_t>((row + col) % 251);
  if (decoded.size() != 1 || decoded[0] != expected) {
    std::fprintf(stderr, "decoded mismatch: got %u expected %u\n",
                 decoded.empty() ? 0 : decoded[0], expected);
    return 1;
  }
  return 0;
}
