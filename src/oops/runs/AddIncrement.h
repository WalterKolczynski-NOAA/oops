/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ADDINCREMENT_H_
#define OOPS_RUNS_ADDINCREMENT_H_

#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

/// Application that adds an increment to a state and writes the sum to a file.
///
/// The increment may optionally be multiplied by a scaling factor and have a different resolution
/// than the state.
template <typename MODEL> class AddIncrement : public Application {
  typedef Geometry<MODEL>  Geometry_;
  typedef State<MODEL>     State_;
  typedef Increment<MODEL> Increment_;

 public:
// -----------------------------------------------------------------------------
  explicit AddIncrement(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~AddIncrement() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Setup resolution
    const Geometry_ stateResol(eckit::LocalConfiguration(fullConfig, "state geometry"),
                               this->getComm());

    const Geometry_ incResol(eckit::LocalConfiguration(fullConfig, "increment geometry"),
                             this->getComm());

//  Read state
    State_ xx(stateResol, eckit::LocalConfiguration(fullConfig, "state"));
    Log::test() << "State: " << xx << std::endl;

//  Read increment
    const eckit::LocalConfiguration incParams(fullConfig, "increment");
    oops::Variables addedVars(incParams, "added variables");
    Increment_ dx(incResol, addedVars, xx.validTime());
    dx.read(incParams);
    Log::test() << "Increment: " << dx << std::endl;

//  Scale increment
    if (incParams.has("scaling factor")) {
      dx *= incParams.getDouble("scaling factor");
      Log::test() << "Scaled the increment: " << dx << std::endl;
    }

//  Assertions on state versus increment
    ASSERT(xx.validTime() == dx.validTime());

//  Add increment to state
    xx += dx;

//  Write state
    xx.write(eckit::LocalConfiguration(fullConfig, "output"));

    Log::test() << "State plus increment: " << xx << std::endl;

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::AddIncrement<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_ADDINCREMENT_H_
