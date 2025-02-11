/*
 * (C) Copyright 2018-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_CONVERTINCREMENT_H_
#define OOPS_RUNS_CONVERTINCREMENT_H_

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/interface/LinearVariableChange.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

/// Options describing each increment read by the ConvertIncrement application.
template <typename MODEL> class IncrementParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(IncrementParameters, Parameters);

 public:
  RequiredParameter<util::DateTime> date{"date", this};
  RequiredParameter<oops::Variables> inputVariables{"input variables", this};
  RequiredParameter<eckit::LocalConfiguration> input{"input", this};
  RequiredParameter<eckit::LocalConfiguration> output{"output", this};
  RequiredParameter<eckit::LocalConfiguration> trajectory{"trajectory", this};
};


/// Top-level options taken by the ConvertIncrement application.
template <typename MODEL> class ConvertIncrementParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(ConvertIncrementParameters, ApplicationParameters);

 public:
  /// Input geometry parameters.
  RequiredParameter<eckit::LocalConfiguration> inputGeometry{"input geometry", this};

  /// Output geometry parameters.
  RequiredParameter<eckit::LocalConfiguration> outputGeometry{"output geometry", this};

  /// Linear variable change.
  OptionalParameter<eckit::LocalConfiguration> linearVarChange{"linear variable change", this};

  /// List of increments.
  RequiredParameter<std::vector<eckit::LocalConfiguration>> increments{"increments", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class ConvertIncrement : public Application {
  typedef Geometry<MODEL>                    Geometry_;
  typedef Increment<MODEL>                   Increment_;
  typedef State<MODEL>                       State_;
  typedef LinearVariableChange<MODEL>        LinearVariableChange_;

  typedef ConvertIncrementParameters<MODEL>  ConvertIncrementParameters_;

 public:
// -------------------------------------------------------------------------------------------------
  explicit ConvertIncrement(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {}
// -------------------------------------------------------------------------------------------------
  virtual ~ConvertIncrement() {}
// -------------------------------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Deserialize parameters
    ConvertIncrementParameters_ params;
    params.deserialize(fullConfig);

//  Setup resolution for intput and output
    const Geometry_ resol1(params.inputGeometry, this->getComm());
    const Geometry_ resol2(params.outputGeometry, this->getComm());

// Check if there is a change of variable defined in Parameters
    bool lvcDefined = false;
    auto linVarChangeParams = params.linearVarChange.value();
    if (linVarChangeParams != boost::none) {
        if (linVarChangeParams.value().has("output variables")) {
            lvcDefined = true;
        }
    }

//  List of input and output increments
    const std::vector<eckit::LocalConfiguration>& incrementParams = params.increments;
    const int nincrements = incrementParams.size();

//  Loop over increments
    for (int jm = 0; jm < nincrements; ++jm) {
//    Print output
      Log::info() << "Converting increment " << jm+1 << " of " << nincrements << std::endl;

//    Datetime for increment
      const util::DateTime incdatetime(incrementParams[jm].getString("date"));

//    Variables for input increment
      const Variables incvars(incrementParams[jm], "input variables");

//    Read input
      const eckit::LocalConfiguration inputParams(incrementParams[jm], "input");
      Increment_ dxi(resol1, incvars, incdatetime);
      dxi.read(inputParams);
      Log::test() << "Input increment: " << dxi << std::endl;

//    Copy and change resolution
      Increment_ dx(resol2, dxi);

//    Variable transform
      if (lvcDefined) {
        const eckit::LocalConfiguration trajConf(incrementParams[jm], "trajectory");
        State_ xTrajBg(resol1, trajConf);
        ASSERT(xTrajBg.validTime() == dx.validTime());  // Check time is consistent
        Log::test() << "Trajectory state: " << xTrajBg << std::endl;

        // Create variable change
        LinearVariableChange_ lvc(resol2, linVarChangeParams.value());
        Variables varout(linVarChangeParams.value(), "output variables");
        lvc.changeVarTraj(xTrajBg, varout);
        if (linVarChangeParams.value().getBool("do inverse", false)) {
          lvc.changeVarInverseTL(dx, varout);
        } else {
          lvc.changeVarTL(dx, varout);
        }
      }

//    Write state
      const eckit::LocalConfiguration outputParams(incrementParams[jm], "output");
      dx.write(outputParams);

      Log::test() << "Output increment: " << dx << std::endl;
    }
    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::ConvertIncrement<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_CONVERTINCREMENT_H_
