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

#if 0
// -----------------------------------------------------------------------------
/// Configuration of the state set test.
template <typename MODEL>
class StateSetTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(StateSetTestParameters, Parameters)

 public:
  /// Configuration of the first state file to load
  oops::RequiredParameter<eckit::LocalConfiguration> statefile1{"statefile1", this};
  /// Configuration of the second state file to load
  oops::RequiredParameter<eckit::LocalConfiguration> statefile2{"statefile2", this};
  /// Validity time for states
  oops::RequiredParameter<util::DateTime> date{"date", this};
  /// Number of ensemble members (determines resolution ratio)
  oops::RequiredParameter<int> nens{"number of members", this};
};

// -----------------------------------------------------------------------------
/// Top-level test parameters.
template <typename MODEL>
class TopTestParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(TopTestParameters, Parameters)

  typedef oops::Geometry<MODEL>           Geometry_;
  typedef StateSetTestParameters<MODEL>   StateSetTestParameters_;

 public:
//  oops::RequiredParameter<StateSetTestParameters_> stateSetTest{"state set test", this};
  oops::RequiredParameter<eckit::LocalConfiguration> daGeomConfig{"da geometry", this};
  oops::RequiredParameter<eckit::LocalConfiguration> fcGeomConfig{"fc geometry", this};
  oops::RequiredParameter<eckit::LocalConfiguration> fcstConfig{"fcst", this};
  oops::RequiredParameter<eckit::LocalConfiguration> ensemble{"ensemble members", this};
  oops::IgnoreOtherParameters ignore{this};
};
#endif

// -----------------------------------------------------------------------------
template <typename MODEL> class StateSetFixture : private boost::noncopyable {
 public:
  typedef oops::Geometry<MODEL>      Geometry_;
//  typedef StateSetTestParameters<MODEL> StateSetTestParameters_;

//  static const StateSetTestParameters_ & test()  {return *getInstance().test_;}
//  static const Geometry_ & daGeom() {return *getInstance().daGeom_;}
//  static const Geometry_ & fcGeom() {return *getInstance().fcGeom_;}
  static void reset() {
/*
    getInstance().daGeom_.reset();
    getInstance().fcGeom_.reset();
    getInstance().test_.reset();
*/
  }
 private:
  static StateSetFixture<MODEL>& getInstance() {
    static StateSetFixture<MODEL> theStateSetFixture;
    return theStateSetFixture;
  }

  StateSetFixture() {
    std::cout << "in statesetfixture" << std::endl;
    /*
    const eckit::LocalConfiguration ensembleConfig(TestEnvironment::config(), "ensemble members"));
    const eckit::LocalConfiguration daGeomConfig(TestEnvironment::config(), "da geometry");
    const eckit::LocalConfiguration fcGeomConfig(TestEnvironment::config(), "fc geometry");
    const eckit::LocalConfiguration fcstConfig(TestEnvironment::config(), "fcst");
    */

/*
    TopTestParameters<MODEL> parameters;
    parameters.validateAndDeserialize(TestEnvironment::config());
    std::cout << "starting dageom" << std::endl;
    daGeom_ = std::make_unique<Geometry_>(parameters.daGeometry, oops::mpi::world());
    std::cout << "starting fcgeom" << std::endl;
    fcGeom_ = std::make_unique<Geometry_>(parameters.fcGeometry, oops::mpi::world());
    */
  }

  ~StateSetFixture<MODEL>() {}

//  std::unique_ptr<StateSetTestParameters_> test_;
//  std::unique_ptr<Geometry_> daGeom_;
//  std::unique_ptr<Geometry_> fcGeom_;
};

// -----------------------------------------------------------------------------
/// \brief tests constructors
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
    std::cout << "da geom" << daGeomConfig << std::endl;
    eckit::LocalConfiguration fcGeomConfig(TestEnvironment::config(), "fc geometry");
    std::cout << "fc geom" << fcGeomConfig << std::endl;
    const eckit::LocalConfiguration fcstConfig(TestEnvironment::config(), "fcst");
    std::cout << "we have nEns = " << nEns << std::endl;
    // DA geom uses all mpi tasks on MPI_COMM_WORLD
    // Get DA geometry configuration from config_
    // Get full configuration

    //FC geometry is created on each ensemble communicator
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
    std::cout << "Hey, myrank is " << mytask << " and mymember is " << mymember << std::endl;
    fcGeomConfig.set("member_number",mymember);
    std::cout << "fcgeomconf is " << fcGeomConfig << std::endl;
    std::cout << "commMember size is " << commMember.size() << std::endl;
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

    fcStateSet.random();
    oops::Log::trace() << "before transpose fc stateset is " << fcStateSet << std::endl;
    // Test transpose functionality between geometries
    // daStateSet is on the daGeom with both ensemble members on each MPI proc
    std::vector<StateSet_> daStateSet = fcStateSet.transpose(worldComm, daGeom, mymember);
    oops::mpi::world().barrier();
    oops::Log::trace() << "size of daStateSet is " << daStateSet.size() << std::endl;
    oops::Log::trace() << "localvec[0] is " << daStateSet[0] << std::endl;
    oops::Log::trace() << "localvec[1] is " << daStateSet[1] << std::endl;
    std::vector<State_> states_;
    // after Rtranspose, states are back to distributed across ensemble ranks
    for(size_t ens=0; ens < daStateSet.size(); ++ens){
//      if((mymember - 1) == ens) { //we only want our ensemble member
        states_.emplace_back(((daStateSet[ens]).Rtranspose(worldComm, fcGeom,
          mymember,ens)));
//      }
    }
    oops::Log::trace() << "size of states_ is " << states_.size() << std::endl;
    oops::Log::trace() << "state after Rtranspose is " << states_[0] << std::endl;
    oops::Log::trace() << "state[1] after Rtranspose is " << states_[1] << std::endl;
    oops::Log::trace() << "times[0] is " << times[0] << std::endl;
    StateSet_ *newFCStateSet = new StateSet_(states_, mymember - 1, times, oops::mpi::myself(),
		                       ensMembers, patchMember);
    oops::Log::trace() << "newFCState after Rtranspose is " << *newFCStateSet << std::endl;
    IncrementSet_ newState(fcGeom, vars, times, oops::mpi::myself(), ensMembers, patchMember);
    newState.diff(fcStateSet,*newFCStateSet);
    oops::Log::trace() << "diff between stateSets is " << newState << std::endl;
    delete newFCStateSet;
    /*

    // Verify dimensions
    EXPECT(daStateSet.size() == nEns);
    EXPECT(daStateSet.size() == nEns * times.size());
    EXPECT(fcStateSet.size() == times.size());
    */
}

template <typename MODEL> void testStateSetConstructors() {
  typedef StateSetFixture<MODEL>     Test_;
  typedef oops::State<MODEL>      State_;
  typedef oops::StateSet<MODEL>   StateSet_;
  typedef oops::Geometry<MODEL>   Geometry_;
  typedef oops::GeometryIterator<MODEL> GeometryIterator_;

#if 0 
  const util::DateTime vt(Test_::test().date);
  const int nens = Test_::test().nens;
  
  // Get the MPI partition
  const int ntasks = oops::mpi::world().size();
  const int mytask = oops::mpi::world().rank();
  const int tasks_per_set = ntasks / 2;  // Split into 2 groups
  const int myset = mytask / tasks_per_set + 1;

  std::cout << "splitting communicators" << std::endl;
  // Create split communicators for each set
  std::string commNameStr = "comm_set_" + std::to_string(myset);
  char const *commName = commNameStr.c_str();
  eckit::mpi::Comm & commSet = oops::mpi::world().split(myset, commName);

  // Create states on DA geometry (higher resolution)
  std::cout << "creating DA states" << std::endl;
  std::unique_ptr<State_> state1(new State_(Test_::daGeom(), Test_::test().statefile1));
  std::unique_ptr<State_> state2(new State_(Test_::daGeom(), Test_::test().statefile2));
  
  std::vector<State_> states;
  states.push_back(*state1);
  states.push_back(*state2);

  // Test main constructor, using split communicators like in LocalEnsembleDA.h
  std::vector<util::DateTime> times = {vt, vt};
  // Create ss1 on first set of processors using DA geometry
  std::unique_ptr<StateSet_> ss1;
  if (myset == 1) {
    ss1.reset(new StateSet_(Test_::daGeom(), state1->variables(), times, commSet));
    EXPECT(ss1.get());
    EXPECT(ss1->size() == 2);
    oops::oops::Log::test() << "Printing DA StateSet on set 1: " << *ss1 << std::endl;
  }

  // Create ss2 on second set of processors using FC geometry
  std::unique_ptr<StateSet_> ss2;
  if (myset == 2) {
    ss2.reset(new StateSet_(Test_::fcGeom(), state1->variables(), times, commSet));
    EXPECT(ss2.get());
    if (ss1) {  // Only check size if ss1 exists on this processor
      EXPECT(ss2->size() == ss1->size());
    }
    oops::oops::Log::test() << "Printing FC StateSet on set 2: " << *ss2 << std::endl;
  }

  // Test geometry resolution ratio
  if (myset == 1) {
    // Get resolution info from both geometries
    GeometryIterator_ daIt = Test_::daGeom().begin();
    GeometryIterator_ fcIt = Test_::fcGeom().begin();
    
    int daPoints = 0;
    int fcPoints = 0;
    
    // Count grid points in each geometry
    while (daIt != Test_::daGeom().end()) {
      ++daPoints;
      ++daIt;
    }
    while (fcIt != Test_::fcGeom().end()) {
      ++fcPoints;
      ++fcIt;
    }
    
    // Check that DA has N times more points than FC
    const double ratio = static_cast<double>(daPoints) / static_cast<double>(fcPoints);
    EXPECT(std::abs(ratio - nens) < 0.1);  // Allow for small rounding differences
    
    oops::oops::Log::test() << "DA points: " << daPoints << ", FC points: " << fcPoints 
                      << ", Ratio: " << ratio << " (expected " << nens << ")" << std::endl;
  }

  // Cleanup
  ss1.reset();
  ss2.reset();
  EXPECT(!ss1.get());
  EXPECT(!ss2.get());
#endif
}

// -----------------------------------------------------------------------------
template <typename MODEL>
class StateSet : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~StateSet() {StateSetFixture<MODEL>::reset();}

 private:
  std::string testid() const override {return "test::StateSet<" + MODEL::name() + ">";}

  void register_tests() const override {
    std::cout << "at beginning of test" << std::endl;
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();
/*
    ts.emplace_back(CASE("interface/StateSet/testStateSetConstructors")
      { testStateSetConstructors<MODEL>(); });
*/
    ts.emplace_back(CASE("interface/StateSet/testStateSetTranspose")
      { testStateSetTranspose<MODEL>(); });
  }

  void clear() const override {}
};

}  // namespace test

#endif  // TEST_INTERFACE_STATESET_H_
