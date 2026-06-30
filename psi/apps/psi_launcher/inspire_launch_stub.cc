#include "yacl/base/exception.h"

#include "psi/apps/psi_launcher/inspire_launch.h"

namespace psi {

PirResultReport RunPir(const InspireReceiverConfig&,
                       const std::shared_ptr<yacl::link::Context>&) {
  YACL_THROW("Inspire is only supported on x86_64");
}

PirResultReport RunPir(const InspireSenderConfig&,
                       const std::shared_ptr<yacl::link::Context>&) {
  YACL_THROW("Inspire is only supported on x86_64");
}

}  // namespace psi
