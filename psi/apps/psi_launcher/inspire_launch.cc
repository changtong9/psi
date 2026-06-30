#include "psi/apps/psi_launcher/inspire_launch.h"

#include "yacl/base/exception.h"

#include "psi/algorithm/inspire/entry.h"

namespace psi {

PirResultReport RunPir(const InspireReceiverConfig& inspire_receiver_config,
                       const std::shared_ptr<yacl::link::Context>& lctx) {
  psi::inspire::InspireReceiverOptions options;
  options.db_rows = inspire_receiver_config.db_rows();
  options.db_cols = inspire_receiver_config.db_cols();
  options.item_size_bits = inspire_receiver_config.item_size_bits();
  options.query_file = inspire_receiver_config.query_file();
  options.output_file = inspire_receiver_config.output_file();

  YACL_ENFORCE_EQ(psi::inspire::ReceiverOnline(options, lctx), 0);
  return PirResultReport();
}

PirResultReport RunPir(const InspireSenderConfig& inspire_sender_config,
                       const std::shared_ptr<yacl::link::Context>& lctx) {
  psi::inspire::InspireSenderOptions options;
  options.db_rows = inspire_sender_config.db_rows();
  options.db_cols = inspire_sender_config.db_cols();
  options.item_size_bits = inspire_sender_config.item_size_bits();
  options.db_file = inspire_sender_config.db_file();

  YACL_ENFORCE_EQ(psi::inspire::SenderOnline(options, lctx), 0);
  return PirResultReport();
}

}  // namespace psi
