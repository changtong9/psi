#pragma once

#include <cstdint>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "yacl/base/buffer.h"
#include "yacl/base/byte_container_view.h"

#include "psi/algorithm/inspire/internal.h"
#include "psi/algorithm/inspire/params.h"
#include "psi/algorithm/inspire/types.h"
#include "psi/algorithm/pir_interface/pir_db.h"

namespace psi::inspire {

class InspireServer : public psi::pir::IndexPirDataBase {
 public:
  explicit InspireServer(InspireParameters params);

  void GenerateFromRawData(const psi::pir::RawDatabase& raw_database) override;
  void GenerateFromSimpleHashTable(
      const psi::pir::RawDatabase& raw_database) override;
  void Dump(std::ostream& out_stream) const override;
  std::size_t MaxElementsOfOnePt() const override { return 1; }
  bool DbSeted() const override { return db_set_; }

  [[nodiscard]] bool DbSet() const { return db_set_; }
  [[nodiscard]] const InspireParameters& GetParameters() const {
    return params_;
  }

  InspirePrecomputedState PerformOfflinePrecomputation() const;
  InspireResponse ProcessQuery(const InspireQuery& query) const;
  InspireResponse ProcessQuery(const InspireQuery& query,
                               const InspirePrecomputedState& state) const;
  yacl::Buffer Response(const yacl::ByteContainerView& query_buffer) const;
  yacl::Buffer Response(const yacl::ByteContainerView& query_buffer,
                        const yacl::Buffer& pks_buffer) const override;
  std::string Response(const yacl::ByteContainerView& query_buffer,
                       const std::string& pks_buffer) const override;

 private:
  InspireParameters params_;
  bool db_set_ = false;
  std::vector<uint8_t> db_row_major_;
  std::unique_ptr<psi::ypir::internal::ypir::Context> context_;
};

}  // namespace psi::inspire
