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
  oops::RequiredParameter<StateSetTestParameters_> stateSetTest{"state set test", this};
  oops::RequiredParameter<eckit::LocalConfiguration> daGeometry{"da geometry", this};
  oops::RequiredParameter<eckit::LocalConfiguration> fcGeometry{"fc geometry", this};
  oops::IgnoreOtherParameters ignore{this};
};

// -----------------------------------------------------------------------------
template <typename MODEL> class StateSetFixture : private boost::noncopyable {
 public:
  typedef oops::Geometry<MODEL>      Geometry_;
  typedef StateSetTestParameters<MODEL> StateSetTestParameters_;

  static const StateSetTestParameters_ & test()  {return *getInstance().test_;}
  static const Geometry_ & daGeom() {return *getInstance().daGeom_;}
  static const Geometry_ & fcGeom() {return *getInstance().fcGeom_;}
  static void reset() {
    getInstance().daGeom_.reset();
    getInstance().fcGeom_.reset();
    getInstance().test_.reset();
  }

 private:
  static StateSetFixture<MODEL>& getInstance() {
    static StateSetFixture<MODEL> theStateSetFixture;
    return theStateSetFixture;
  }

  StateSetFixture() {
    TopTestParameters<MODEL> parameters;
    parameters.validateAndDeserialize(TestEnvironment::config());
    std::cout << "in statesetfixture" << std::endl;

    test_ = std::make_unique<StateSetTestParameters_>(parameters.stateSetTest);
    std::cout << "starting dageom" << std::endl;
    daGeom_ = std::make_unique<Geometry_>(parameters.daGeometry, oops::mpi::world());
    std::cout << "starting fcgeom" << std::endl;
    fcGeom_ = std::make_unique<Geometry_>(parameters.fcGeometry, oops::mpi::world());
  }

  ~StateSetFixture<MODEL>() {}

  std::unique_ptr<StateSetTestParameters_> test_;
  std::unique_ptr<Geometry_> daGeom_;
  std::unique_ptr<Geometry_> fcGeom_;
};

// -----------------------------------------------------------------------------
/// \brief tests constructors
template <typename MODEL> void testStateSetConstructors() {
  typedef StateSetFixture<MODEL>     Test_;
  typedef oops::State<MODEL>      State_;
  typedef oops::StateSet<MODEL>   StateSet_;
  typedef oops::Geometry<MODEL>   Geometry_;
  typedef oops::GeometryIterator<MODEL> GeometryIterator_;

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
#if 0 
  // Create ss1 on first set of processors using DA geometry
  std::unique_ptr<StateSet_> ss1;
  if (myset == 1) {
    ss1.reset(new StateSet_(Test_::daGeom(), state1->variables(), times, commSet));
    EXPECT(ss1.get());
    EXPECT(ss1->size() == 2);
    oops::Log::test() << "Printing DA StateSet on set 1: " << *ss1 << std::endl;
  }

  // Create ss2 on second set of processors using FC geometry
  std::unique_ptr<StateSet_> ss2;
  if (myset == 2) {
    ss2.reset(new StateSet_(Test_::fcGeom(), state1->variables(), times, commSet));
    EXPECT(ss2.get());
    if (ss1) {  // Only check size if ss1 exists on this processor
      EXPECT(ss2->size() == ss1->size());
    }
    oops::Log::test() << "Printing FC StateSet on set 2: " << *ss2 << std::endl;
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
    
    oops::Log::test() << "DA points: " << daPoints << ", FC points: " << fcPoints 
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

    ts.emplace_back(CASE("interface/StateSet/testStateSetConstructors")
      { testStateSetConstructors<MODEL>(); });
  }

  void clear() const override {}
};

}  // namespace test

#endif  // TEST_INTERFACE_STATESET_H_
