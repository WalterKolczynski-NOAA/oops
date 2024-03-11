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
                  StateSet_ & stateSet );

  /// calculate ensemble mean
  StateSet_ mean() const;

  /// Accessors
  unsigned int size() const { return states_.size(); }
  State_ & operator()(const int ii) { return (stateSet_)[ii]; }
  const State_ & operator()(const int ii) const { return (stateSet_)[ii]; }
  StateSet_ & operator[](const int ii) { return states_[ii]; }
  const StateSet_ & operator[](const int ii) const { return states_[ii]; }

  /// Information
  const Variables & variables() const {return states_[0].variables();}

 private:
  std::vector<StateSet_> states_;
  StateSet_ stateSet_;
};

// ====================================================================================

template<typename MODEL>
StateEnsemble4D<MODEL>::StateEnsemble4D(const Geometry_ &, const eckit::Configuration &,
                  StateSet_ & stateSet ): states_(), stateSet_(stateSet) {

  // copy stateSet into class
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
  : states_(), stateSet_(resol,vars,times,commTime,ens,commEns) {
  // Abort if both "members" and "members from template" are specified
  if (config.has("members") && config.has("members from template"))
    ABORT("StateEnsemble4D:constructor: both members and members from template are specified");

  std::vector<eckit::LocalConfiguration> membersConfig;
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

  // Reserve memory to hold ensemble
  states_.reserve(membersConfig.size());

  // Loop over all ensemble members
  for (size_t jj = 0; jj < membersConfig.size(); ++jj) {
    states_.emplace_back(StateSet_(resol, membersConfig[jj]));
  }
  for( size_t jt = 0; jt < times.size(); ++jt) { // FIX THIS. 
//    stateSet_[jt] = State(resol,membersConfig[mymember]);
    stateSet_[jt] = (states_[mymember])[jt];
  }
  Log::trace() << "StateEnsemble4D:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL>
StateEnsemble4D<MODEL>::StateEnsemble4D(const Geometry_ & resol,
                                        const eckit::Configuration & config)
  : states_(), stateSet_(resol,config) {
  // Abort if both "members" and "members from template" are specified
  if (config.has("members") && config.has("members from template"))
    ABORT("StateEnsemble4D:constructor: both members and members from template are specified");

  std::vector<eckit::LocalConfiguration> membersConfig;
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

  return states_[0].ens_mean();
/*
  // Compute ensemble mean
  Accumulator<MODEL, StateSet_, StateSet_> ensmean(states_[0]);

  const double rr = 1.0/static_cast<double>(states_.size());
  for (size_t iens = 0; iens < states_.size(); ++iens) {
    ensmean.accumul(rr, states_[iens]);
  }

  Log::trace() << "StateEnsemble4D::mean done" << std::endl;
  return std::move(ensmean);
*/
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_STATEENSEMBLE4D_H_
