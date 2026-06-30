#include "psi/algorithm/inspire/client.h"

#include <string>
#include <utility>

#include "yacl/base/exception.h"

#include "psi/algorithm/inspire/serialize.h"
#include "psi/algorithm/ypir/ypir_internal_params.h"

namespace psi::inspire {

InspireClient::InspireClient(InspireParameters params)
    : params_(std::move(params)),
      context_(std::make_unique<psi::ypir::internal::ypir::Context>(
          psi::ypir::internal::ypir::CreateContext(params_))) {
  YACL_ENFORCE(params_.mode == psi::ypir::YpirMode::kDoublepir);
}

yacl::Buffer InspireClient::GeneratePksBuffer() const {
  return yacl::Buffer();
}

std::string InspireClient::GeneratePksString() const { return {}; }

InspireQuery InspireClient::GenerateQuery(uint64_t raw_idx) const {
  YACL_ENFORCE_LT(raw_idx, params_.NumItems());
  return internal::GenerateQuery(raw_idx, params_, client_secrets_, *context_);
}

yacl::Buffer InspireClient::GenerateQueryBuffer(uint64_t raw_idx) const {
  return SerializeQuery(GenerateQuery(raw_idx));
}

yacl::Buffer InspireClient::GenerateIndexQuery(uint64_t raw_idx) const {
  return GenerateQueryBuffer(raw_idx);
}

std::string InspireClient::GenerateIndexQueryStr(uint64_t raw_idx) const {
  auto buffer = GenerateQueryBuffer(raw_idx);
  return std::string(static_cast<std::string_view>(buffer));
}

std::vector<uint8_t> InspireClient::DecodeResponse(
    const InspireResponse& response, uint64_t raw_idx) const {
  YACL_ENFORCE_LT(raw_idx, params_.NumItems());
  return internal::RecoverResponse(response, params_, client_secrets_,
                                   *context_);
}

std::vector<uint8_t> InspireClient::DecodeResponseBuffer(
    const yacl::ByteContainerView& response_buffer, uint64_t raw_idx) const {
  return DecodeResponse(DeserializeResponse(response_buffer), raw_idx);
}

std::vector<uint8_t> InspireClient::DecodeIndexResponse(
    const yacl::ByteContainerView& response_buffer, uint64_t raw_idx) const {
  return DecodeResponseBuffer(response_buffer, raw_idx);
}

}  // namespace psi::inspire
