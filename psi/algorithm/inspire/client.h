#pragma once

#include <memory>
#include <string>
#include <vector>

#include "yacl/base/buffer.h"
#include "yacl/base/byte_container_view.h"

#include "psi/algorithm/inspire/internal.h"
#include "psi/algorithm/inspire/params.h"
#include "psi/algorithm/inspire/types.h"
#include "psi/algorithm/pir_interface/index_pir.h"

namespace psi::inspire {

class InspireClient : public psi::pir::IndexPirClient {
 public:
  explicit InspireClient(InspireParameters params);

  const InspireParameters& GetParameters() const { return params_; }

  pir::PirType GetPirType() const override { return pir::PirType::YPIR_PIR; }

  yacl::Buffer GeneratePksBuffer() const override;
  std::string GeneratePksString() const override;

  InspireQuery GenerateQuery(uint64_t raw_idx) const;
  yacl::Buffer GenerateQueryBuffer(uint64_t raw_idx) const;
  yacl::Buffer GenerateIndexQuery(uint64_t raw_idx) const override;
  std::string GenerateIndexQueryStr(uint64_t raw_idx) const override;

  std::vector<uint8_t> DecodeResponse(const InspireResponse& response,
                                      uint64_t raw_idx) const;
  std::vector<uint8_t> DecodeResponseBuffer(
      const yacl::ByteContainerView& response_buffer, uint64_t raw_idx) const;
  std::vector<uint8_t> DecodeIndexResponse(
      const yacl::ByteContainerView& response_buffer,
      uint64_t raw_idx) const override;

 private:
  InspireParameters params_;
  std::unique_ptr<psi::ypir::internal::ypir::Context> context_;
  mutable internal::ClientSecrets client_secrets_;
};

}  // namespace psi::inspire
