/*
 * (C) Copyright 2019-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_STATEENSEMBLESET_H_
#define OOPS_BASE_STATEENSEMBLESET_H_

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

/// \brief Ensemble of stateSets
template<typename MODEL> class StateEnsembleSet {
  typedef Geometry<MODEL>      Geometry_;
  typedef StateSet<MODEL>      StateSet_;
  typedef State<MODEL>         State_;

 public:
  /// Create ensemble of stateSets
  StateEnsembleSet(const Geometry_ &, const eckit::Configuration &);

  StateEnsembleSet(const Geometry_ & resol,
                                     const eckit::Configuration & config,
                                     const Variables & vars,
                                     const std::vector<util::DateTime> & times,
                                     const eckit::mpi::Comm & commTime,
                                     const std::vector<int> & ens,
                                     const eckit::mpi::Comm & commEns,
                                     const int mymember);

  /// Create ensemble of stateSets
  StateEnsembleSet(const Geometry_ &, const eckit::Configuration &,
                  StateSet_ & stateSet );

  /// calculate ensemble mean
  StateSet_ mean() const;

  /// Accessors
  unsigned int size() const { return stateSet_.size(); }
  State_ & operator()(const int ii) { return (stateSet_)[ii]; }
  const State_ & operator()(const int ii) const { return (stateSet_)[ii]; }
//  StateSet_ & operator[](const int ii) { return states_[ii]; }
//  const StateSet_ & operator[](const int ii) const { return states_[ii]; }
  /// Information
  const Variables & variables() const {return stateSet_.variables();}
  const StateSet_ & stateSet() const {return stateSet_;}

 private:
  StateSet_ stateSet_;
};

// ====================================================================================

template<typename MODEL>
StateEnsembleSet<MODEL>::StateEnsembleSet(const Geometry_ &, const eckit::Configuration &,
                  StateSet_ & stateSet ): stateSet_(stateSet) {

  Log::trace() << "StateEnsembleSet:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StateEnsembleSet<MODEL>::StateEnsembleSet(const Geometry_ & resol,
                                        const eckit::Configuration & config,
                                        const Variables & vars,
                                        const std::vector<util::DateTime> & times,
                                        const eckit::mpi::Comm & commTime,
                                        const std::vector<int> & ens,
                                        const eckit::mpi::Comm & commEns,
                                        const int mymember)
  : stateSet_(resol,vars,times,commTime,ens,commEns) {
  // Abort if both "members" and "members from template" are specified
  if (config.has("members") && config.has("members from template"))
    ABORT("StateEnsembleSet:constructor: both members and members from template are specified");

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
    ABORT("StateEnsembleSet: ensemble not specified");
  }

  std::vector<StateSet_> states_;
  // Reserve memory to hold ensemble
  states_.reserve(times.size());

  // read in ensemble members on appropriate communicator
  states_.emplace_back(StateSet_(resol, membersConfig[mymember-1]));

  for (size_t jj = 0; jj < stateSet_.local_ens_size(); ++jj) {
    for( size_t jt = 0; jt < times.size(); ++jt) { // FIX THIS. 
      stateSet_[jt] = (states_[jj])[jt];
    }
  }
  Log::trace() << "StateEnsembleSet:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL>
StateEnsembleSet<MODEL>::StateEnsembleSet(const Geometry_ & resol,
                                        const eckit::Configuration & config)
  : stateSet_(resol,config) {
  // Abort if both "members" and "members from template" are specified
  if (config.has("members") && config.has("members from template"))
    ABORT("StateEnsembleSet:constructor: both members and members from template are specified");

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
    ABORT("StateEnsembleSet: ensemble not specified");
  }

  // Reserve memory to hold ensemble
  states_.reserve(membersConfig.size());

  // Loop over all ensemble members
  for (size_t jj = 0; jj < membersConfig.size(); ++jj) {
    states_.emplace_back(StateSet_(resol, membersConfig[jj]));
  }
  Log::trace() << "StateEnsembleSet:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StateSet<MODEL> StateEnsembleSet<MODEL>::mean() const {

  return stateSet_.ens_mean();
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_STATEENSEMBLESET_H_
