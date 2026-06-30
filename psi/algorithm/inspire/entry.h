#pragma once

#include <cstdint>
#include <memory>
#include <string>

#include "yacl/link/context.h"

namespace psi::inspire {

struct InspireSenderOptions {
  uint64_t db_rows = 0;
  uint64_t db_cols = 0;
  uint64_t item_size_bits = 0;
  std::string db_file;
};

struct InspireReceiverOptions {
  uint64_t db_rows = 0;
  uint64_t db_cols = 0;
  uint64_t item_size_bits = 0;
  std::string query_file;
  std::string output_file;
};

int SenderOnline(const InspireSenderOptions& options,
                 std::shared_ptr<yacl::link::Context> lctx);

int ReceiverOnline(const InspireReceiverOptions& options,
                   std::shared_ptr<yacl::link::Context> lctx);

}  // namespace psi::inspire
