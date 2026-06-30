#pragma once

#include "psi/algorithm/ypir/params.h"

namespace psi::inspire {

using InspireParameters = psi::ypir::YpirParameters;

inline InspireParameters CreateParamsForScenario(uint64_t num_items,
                                                 uint64_t item_size_bits) {
  return psi::ypir::CreateParamsForScenarioDoublePIR(num_items,
                                                     item_size_bits);
}

inline InspireParameters CreateParamsForShape(uint64_t db_rows,
                                              uint64_t db_cols,
                                              uint64_t item_size_bits) {
  return psi::ypir::CreateParamsForShapeDoublePIR(db_rows, db_cols,
                                                  item_size_bits);
}

inline InspireParameters CreateSmallTestParams() {
  return psi::ypir::CreateSmallTestParamsDoublePIR();
}

}  // namespace psi::inspire
