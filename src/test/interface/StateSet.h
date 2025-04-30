/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef TEST_INTERFACE_STATESET_H_
#define TEST_INTERFACE_STATESET_H_

#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Geometry.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/State.h"
#include "oops/base/StateSet.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/IgnoreOtherParameters.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"
#include "test/interface/State.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------
/// \brief tests transpose and reverseTranspose
///
template <typename MODEL>
 void testStateSetTranspose() {
    typedef oops::Geometry<MODEL>  Geometry_;
    typedef oops::State<MODEL>  State_;
    typedef oops::StateSet<MODEL>  StateSet_;
    typedef oops::IncrementSet<MODEL>  IncrementSet_;
    const eckit::Configuration & config = test::TestEnvironment::config();
    int nEns;  // Number of ensemble members

    config.get("ensemble members", nEns);
    const eckit::LocalConfiguration daGeomConfig(TestEnvironment::config(), "da geometry");
    eckit::LocalConfiguration fcGeomConfig(TestEnvironment::config(), "fc geometry");
    const eckit::LocalConfiguration fcstConfig(TestEnvironment::config(), "fcst");
    // DA geom uses all mpi tasks on MPI_COMM_WORLD
    // Get DA geometry configuration from config_
    // Get full configuration

    // FC geometry is created on each ensemble communicator
    // Get FC geometry configuration from config_

    // Create DA communicator
    const eckit::mpi::Comm & worldComm = oops::mpi::world();
    const int ntasks = worldComm.size();
    const int mytask = worldComm.rank();
    const int tasks_per_member = ntasks / nEns;
    const int mymember = mytask / tasks_per_member + 1;

    // Create member communicator
    std::string commNameStr = "comm_member_" + std::to_string(mymember);
    char const *commName = commNameStr.c_str();
    eckit::mpi::Comm & commMember = worldComm.split(mymember, commName);
    const int subrank = commMember.rank();

    // Create patch communicator
    std::string patchNameStr = "patch_member_" + std::to_string(subrank);
    char const *patchName = patchNameStr.c_str();
    eckit::mpi::Comm & patchMember = worldComm.split(subrank, patchName);

    oops::Log::info() << "size of patchMember/ENS comm is " << patchMember.size() << std::endl;
    // Create geometries using appropriate communicators
    Geometry_ daGeom(daGeomConfig, worldComm);
    fcGeomConfig.set("member_number", mymember);
    Geometry_ fcGeom(fcGeomConfig, commMember);

    //  Setup times
    oops::Log::info() << "setting up times" << std::endl;
    eckit::LocalConfiguration model = fcstConfig.getSubConfiguration("model");
    const util::Duration tstep(model.getString("tstep"));
    eckit::LocalConfiguration ic = fcstConfig.getSubConfiguration("initial condition");
    const util::DateTime bgndate(ic.getString("datetime"));
    const util::Duration fclength(fcstConfig.getString("forecast length"));
    const util::DateTime enddate(bgndate + fclength);
    std::vector<util::DateTime> times;
    const oops::Variables vars(ic, "state variables");
    // Create ensemble members vector
    std::vector<int> ensMembers;
    for (int i = 1; i <= nEns; i++) {
      ensMembers.push_back(i);
    }
    for (util::DateTime ii=(bgndate+tstep); ii <= enddate; ii=ii+tstep) {
       oops::Log::info() << "pushing back time " << ii << std::endl;
       times.push_back(ii);
    }
    oops::mpi::world().barrier();

    // Create StateSet for FC geometry
    StateSet_ fcStateSet(fcGeom, vars, times, oops::mpi::myself(),
                              ensMembers, patchMember);

    std::random_device rd;  // Used to obtain a seed for the random number engine
    //  std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::mt19937 gen(123+oops::mpi::world().rank());

    // Create distribution for the range you want (e.g., between 0.0 and 1.0)
    std::uniform_real_distribution<double> dis(0.0d, 1.0d);
    std::vector<std::vector<double> > zz(fcStateSet.local_ens_size());

    size_t dataSize = fcStateSet(0, 0).serialSize() -3;


  // Put random numbers in each State
    for (size_t jm = 0; jm < fcStateSet.local_ens_size(); ++jm) {
      fcStateSet(0, jm).serialize(zz[jm]);
      for (size_t i = 0; i < dataSize; ++i) {
      // Generate a random double
          zz[jm][i] = dis(gen);
      }
    }

  // deserialize back to stateSet
    for (size_t jt = 0; jt < fcStateSet.local_time_size(); ++jt) {
      // Put States from each local ensemble member in a vector
      for (size_t jm = 0; jm < fcStateSet.local_ens_size(); ++jm) {
        // serialize local ensembles
        size_t indx = 0;
        fcStateSet(jt, jm).deserialize(zz[jm], indx);
      }
    }
  // Put random numbers in each State
    for (size_t jm = 0; jm < fcStateSet.local_ens_size(); ++jm) {
      fcStateSet(0, jm).serialize(zz[jm]);
    }
    oops::Log::info() << "before transpose fc stateset is " << fcStateSet << std::endl;
    // Test transpose functionality between geometries
    // daStateSet is on the daGeom with both ensemble members on each MPI proc
    std::vector<StateSet_> daStateSet = fcStateSet.transpose(worldComm, daGeom, mymember);
    oops::mpi::world().barrier();
    std::vector<State_> states_;
    // after reverseTranspose, states are back to distributed across ensemble ranks
    for (size_t ens=0; ens < daStateSet.size(); ++ens) {
        states_.emplace_back(((daStateSet[ens]).reverseTranspose(worldComm, fcGeom,
          mymember, ens)));
    }
    oops::Log::info() << "state[0] after reverseTranspose is " << states_[0] << std::endl;
    StateSet_ *newFCStateSet = new StateSet_(states_, mymember - 1, times, oops::mpi::myself(),
                   ensMembers, patchMember);
    IncrementSet_ newState(fcGeom, vars, times, oops::mpi::myself(), ensMembers, patchMember);
    newState.diff(fcStateSet, *newFCStateSet);
    // This should be zero for all variables 
    oops::Log::info() << "diff between stateSets is " << newState << std::endl;
    // Verify that the norm is zero
    EXPECT(newState[0].norm() == 0.0);

    delete newFCStateSet;
}

// -----------------------------------------------------------------------------
template <typename MODEL>
class StateSet : public oops::Test {
 public:
  using oops::Test::Test;

 private:
  std::string testid() const override {return "test::StateSet<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/StateSet/testStateSetTranspose")
      { testStateSetTranspose<MODEL>(); });
  }

  void clear() const override {}
};

}  // namespace test

#endif  // TEST_INTERFACE_STATESET_H_
