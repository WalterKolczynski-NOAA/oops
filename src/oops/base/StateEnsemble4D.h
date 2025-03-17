/*
 * (C) Copyright 2019-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_STATEENSEMBLE4D_H_
#define OOPS_BASE_STATEENSEMBLE4D_H_

#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Accumulator.h"
#include "oops/base/Geometry.h"
#include "oops/base/State.h"
#include "oops/base/StateSet.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/Logger.h"

namespace oops {

class Variables;

// -----------------------------------------------------------------------------

/// \brief Ensemble of 4D states
template<typename MODEL> class StateEnsemble4D {
  typedef Geometry<MODEL>      Geometry_;
  typedef StateSet<MODEL>      StateSet_;
  typedef State<MODEL>         State_;

 public:
  /// Create ensemble of 4D states
  StateEnsemble4D(const Geometry_ &, const eckit::Configuration &);

  StateEnsemble4D(const Geometry_ & resol,
                                     const eckit::Configuration & config,
                                     const Variables & vars,
                                     const std::vector<util::DateTime> & times,
                                     const eckit::mpi::Comm & commTime,
                                     const std::vector<int> & ens,
                                     const eckit::mpi::Comm & commEns,
                                     const int mymember);

  /// Create ensemble of 4D states
  StateEnsemble4D(const Geometry_ &, const eckit::Configuration &,
                  StateSet_ & stateSet);

  explicit StateEnsemble4D(std::vector<StateSet_> & stateSetVec, const int ensNum);

  /// calculate ensemble mean
  StateSet_ mean() const;

  /// Accessors
  unsigned int size() const { return states_.size(); }
  StateSet_ & operator[](const size_t ii) { return states_[ii]; }
  const StateSet_ & operator[](const size_t ii) const { return states_[ii]; }

  /// Information
  const Variables & variables() const {return states_[0].variables();}

 private:
  void getMembers(const eckit::Configuration &);
  std::vector<eckit::LocalConfiguration> membersConfig;
  // TODO(mpotts) remove std::vector<StateSet_> and replace it with a singular StateSet_
  // variable states_
  std::vector<StateSet_> states_;
};

// ====================================================================================

template<typename MODEL>
StateEnsemble4D<MODEL>::StateEnsemble4D(std::vector<StateSet_> & stateSetVec, const int ensNum) :
  states_(stateSetVec) {
  Log::trace() << "StateEnsemble4D:contructor done" << std::endl;
}

// ====================================================================================

template<typename MODEL>
StateEnsemble4D<MODEL>::StateEnsemble4D(const Geometry_ & resol,
      const eckit::Configuration & config, StateSet_ & stateSet) : states_() {
  Log::trace() << "StateEnsemble4D:contructor starting" << std::endl;
  states_.emplace_back(stateSet);

  Log::trace() << "StateEnsemble4D:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StateEnsemble4D<MODEL>::StateEnsemble4D(const Geometry_ & resol,
                                        const eckit::Configuration & config,
                                        const Variables & vars,
                                        const std::vector<util::DateTime> & times,
                                        const eckit::mpi::Comm & commTime,
                                        const std::vector<int> & ens,
                                        const eckit::mpi::Comm & commEns,
                                        const int mymember)
  : states_() {
  // Abort if both "members" and "members from template" are specified
  getMembers(config);

  // Reserve memory to hold ensemble
  states_.reserve(times.size());

  // read in ensemble members on appropriate communicator
  states_.emplace_back(StateSet_(resol, membersConfig[mymember-1]));

  Log::trace() << "StateEnsemble4D:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL>
StateEnsemble4D<MODEL>::StateEnsemble4D(const Geometry_ & resol,
                                        const eckit::Configuration & config)
  : states_() {
  // Abort if both "members" and "members from template" are specified
  getMembers(config);

  // Reserve memory to hold ensemble
  states_.reserve(membersConfig.size());

  // Loop over all ensemble members
  for (size_t jj = 0; jj < membersConfig.size(); ++jj) {
    states_.emplace_back(StateSet_(resol, membersConfig[jj]));
  }
  Log::trace() << "StateEnsemble4D:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StateSet<MODEL> StateEnsemble4D<MODEL>::mean() const {
  // Compute ensemble mean
  Accumulator<MODEL, StateSet_, StateSet_> ensmean(states_[0]);

  const double rr = 1.0/static_cast<double>(states_.size());
  Log::info() << "calculating mean with states_.size() of " << states_.size() << std::endl;
  for (size_t iens = 0; iens < states_.size(); ++iens) {
    ensmean.accumul(rr, states_[iens]);
  }

  Log::trace() << "StateEnsemble4D::mean done" << std::endl;
  return std::move(ensmean);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void StateEnsemble4D<MODEL>::getMembers(const eckit::Configuration & config) {
  if (config.has("members") && config.has("members from template"))
    ABORT("StateEnsemble4D:constructor: both members and members from template are specified");

  if (config.has("members")) {
    // Explicit members
    config.get("members", membersConfig);
  } else if (config.has("members from template")) {
    // Templated members
    eckit::LocalConfiguration templateConfig;
    config.get("members from template", templateConfig);
    eckit::LocalConfiguration membersTemplate;
    templateConfig.get("template", membersTemplate);
    std::string pattern;
    templateConfig.get("pattern", pattern);
    int ne;
    templateConfig.get("nmembers", ne);
    int start = 1;
    if (templateConfig.has("start")) {
      templateConfig.get("start", start);
    }
    std::vector<int> except;
    if (templateConfig.has("except")) {
      templateConfig.get("except", except);
    }
    int zpad = 0;
    if (templateConfig.has("zero padding")) {
      templateConfig.get("zero padding", zpad);
    }
    int count = start;
    for (int ie=0; ie < ne; ++ie) {
      while (std::count(except.begin(), except.end(), count)) {
        count += 1;
      }
      eckit::LocalConfiguration memberConfig(membersTemplate);
      util::seekAndReplace(memberConfig, pattern, count, zpad);
      membersConfig.push_back(memberConfig);
      count += 1;
    }
  } else {
    ABORT("StateEnsemble4D: ensemble not specified");
  }
}

}  // namespace oops

#endif  // OOPS_BASE_STATEENSEMBLE4D_H_
