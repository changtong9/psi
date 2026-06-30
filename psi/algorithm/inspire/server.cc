#include "psi/algorithm/inspire/server.h"

#include <cstring>
#include <ostream>
#include <string>
#include <utility>

#include "yacl/base/exception.h"

#include "psi/algorithm/inspire/serialize.h"
#include "psi/algorithm/ypir/ypir_internal_params.h"

namespace psi::inspire {

InspireServer::InspireServer(InspireParameters params)
    : psi::pir::IndexPirDataBase(psi::pir::PirType::YPIR_PIR),
      params_(std::move(params)),
      context_(std::make_unique<psi::ypir::internal::ypir::Context>(
          psi::ypir::internal::ypir::CreateContext(params_))) {
  YACL_ENFORCE(params_.mode == psi::ypir::YpirMode::kDoublepir);
  YACL_ENFORCE_EQ(params_.value_bytes, 1U,
                  "Inspire currently expects uint8_t database values");
}

void InspireServer::GenerateFromRawData(
    const psi::pir::RawDatabase& raw_database) {
  const bool item_layout = raw_database.Rows() <= params_.NumItems() &&
                           raw_database.RowByteLen() == params_.value_bytes;
  const bool matrix_layout = raw_database.Rows() == params_.db_rows &&
                             raw_database.RowByteLen() == params_.db_cols;
  YACL_ENFORCE(item_layout || matrix_layout,
               "raw database shape does not match Inspire parameters");

  db_row_major_.assign(params_.db_rows * params_.db_cols, 0);
  if (item_layout) {
    for (uint64_t raw_idx = 0; raw_idx < raw_database.Rows(); ++raw_idx) {
      db_row_major_[raw_idx] = raw_database.At(raw_idx)[0];
    }
  } else {
    for (uint64_t row = 0; row < params_.db_rows; ++row) {
      const auto& row_bytes = raw_database.At(row);
      std::memcpy(db_row_major_.data() + row * params_.db_cols,
                  row_bytes.data(), row_bytes.size());
    }
  }
  db_set_ = true;
}

void InspireServer::GenerateFromSimpleHashTable(
    const psi::pir::RawDatabase& raw_database) {
  GenerateFromRawData(raw_database);
}

void InspireServer::Dump(std::ostream& out_stream) const {
  out_stream << "InspireServer{db_rows=" << params_.db_rows
             << ", db_cols=" << params_.db_cols
             << ", value_bytes=" << params_.value_bytes
             << ", db_set=" << db_set_ << "}";
}

InspirePrecomputedState InspireServer::PerformOfflinePrecomputation() const {
  YACL_ENFORCE(db_set_, "database must be loaded before precomputation");
  return internal::PrepareOfflineState(db_row_major_, params_, *context_);
}

InspireResponse InspireServer::ProcessQuery(const InspireQuery& query) const {
  return ProcessQuery(query, PerformOfflinePrecomputation());
}

InspireResponse InspireServer::ProcessQuery(
    const InspireQuery& query, const InspirePrecomputedState& state) const {
  YACL_ENFORCE(db_set_, "database must be loaded before query processing");
  return internal::ProcessQuery(db_row_major_, query, state, params_, *context_);
}

yacl::Buffer InspireServer::Response(
    const yacl::ByteContainerView& query_buffer) const {
  return SerializeResponse(ProcessQuery(DeserializeQuery(query_buffer)));
}

yacl::Buffer InspireServer::Response(
    const yacl::ByteContainerView& query_buffer,
    const yacl::Buffer& /*pks_buffer*/) const {
  return Response(query_buffer);
}

std::string InspireServer::Response(const yacl::ByteContainerView& query_buffer,
                                    const std::string& /*pks_buffer*/) const {
  auto buffer = Response(query_buffer);
  return std::string(static_cast<std::string_view>(buffer));
}

}  // namespace psi::inspire
