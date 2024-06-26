/*
 * (C) Copyright 2021- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_PSEUDOMODELSTATESET_H_
#define OOPS_GENERIC_PSEUDOMODELSTATESET_H_

#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/base/StateSet.h"
#include "oops/generic/ModelBase.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

/// Generic implementation of the pseudo model initialized with 4D State
/// (steps through time by stepping through states in 4D state)
template <typename MODEL>
class PseudoModelStateSet : public ModelBase<MODEL> {
  typedef Geometry<MODEL>          Geometry_;
  typedef ModelAuxControl<MODEL>   ModelAux_;
  typedef State<MODEL>             State_;
  typedef StateSet<MODEL>           StateSet_;

 public:
  static const std::string classname() {return "oops::PseudoModelStateSet";}

  /// Initialize pseudo model with \p stateset - 4D state to loop through
  /// in the model run and \p tstep - time resolution of the model
  PseudoModelStateSet(const StateSet_ & stateset,
                     const util::Duration & tstep = util::Duration(0));

  /// initialize forecast
  void initialize(State_ &) const override;
  /// one forecast step
  void step(State_ &, const ModelAux_ &) const override;
  /// finalize forecast
  void finalize(State_ &) const override;

  /// model time step
  const util::Duration & timeResolution() const override {return tstep_;}

 private:
  void print(std::ostream &) const override;

  /// Reference to 4D state that is used in the model
  const StateSet_ & stateset_;
  /// Model's time resolution
  util::Duration   tstep_;
  /// Index of the current state
  mutable size_t currentstate_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
PseudoModelStateSet<MODEL>::PseudoModelStateSet(const StateSet_ & stateset,
                                              const util::Duration & tstep)
  : stateset_(stateset), tstep_(tstep) {
  const std::vector<util::DateTime> validTimes = stateset_.validTimes();
  if (validTimes.size() > 1) tstep_ = validTimes[1] - validTimes[0];
  Log::trace() << "PseudoModelStateSet<MODEL>::PseudoModelStateSet done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModelStateSet<MODEL>::initialize(State_ & xx) const {
  currentstate_ = 0;
  xx = stateset_[currentstate_];
  Log::trace() << "PseudoModelStateSet<MODEL>::initialize done currentstate is " << currentstate_ << std::endl;
  Log::trace() << "PseudoModelStateSet<MODEL>::initialize done size is " << stateset_[currentstate_] << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModelStateSet<MODEL>::step(State_ & xx, const ModelAux_ & merr) const {
  Log::trace() << "PseudoModelStateSet<MODEL>:step Starting " << currentstate_ << std::endl;
  xx = stateset_[currentstate_++];
  Log::trace() << "PseudoModelStateSet<MODEL>::step done" << currentstate_ << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModelStateSet<MODEL>::finalize(State_ & xx) const {
  Log::trace() << "PseudoModelStateSet<MODEL>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void PseudoModelStateSet<MODEL>::print(std::ostream & os) const {
  os << "Pseudo model stepping through 4D state with " << tstep_ << " time resolution";
}

}  // namespace oops

#endif  // OOPS_GENERIC_PSEUDOMODELSTATESET_H_
