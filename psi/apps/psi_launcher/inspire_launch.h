#pragma once

#include <memory>

#include "yacl/link/context.h"

#include "psi/proto/pir.pb.h"

namespace psi {

PirResultReport RunPir(const InspireReceiverConfig& inspire_receiver_config,
                       const std::shared_ptr<yacl::link::Context>& lctx);

PirResultReport RunPir(const InspireSenderConfig& inspire_sender_config,
                       const std::shared_ptr<yacl::link::Context>& lctx);

}  // namespace psi
