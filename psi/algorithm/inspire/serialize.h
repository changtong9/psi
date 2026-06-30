#pragma once

#include "yacl/base/buffer.h"
#include "yacl/base/byte_container_view.h"

#include "psi/algorithm/inspire/types.h"

namespace psi::inspire {

yacl::Buffer SerializeQuery(const InspireQuery& query);
InspireQuery DeserializeQuery(const yacl::ByteContainerView& buffer);

yacl::Buffer SerializeResponse(const InspireResponse& response);
InspireResponse DeserializeResponse(const yacl::ByteContainerView& buffer);

}  // namespace psi::inspire
