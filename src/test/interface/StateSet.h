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
  oops::RequiredParameter<eckit::LocalConfiguration> geometry{"geometry", this};
  oops::IgnoreOtherParameters ignore{this};
};


// -----------------------------------------------------------------------------

template <typename MODEL> class StateSetFixture : private boost::noncopyable {
 public:
  typedef oops::Geometry<MODEL>      Geometry_;
  typedef StateSetTestParameters<MODEL> StateSetTestParameters_;

  static const StateSetTestParameters_ & test()  {return *getInstance().test_;}
  static const Geometry_            & resol() {return *getInstance().resol_;}
  static void reset() {
    getInstance().resol_.reset();
    getInstance().test_.reset();
  }

 private:
  static StateSetFixture<MODEL>& getInstance() {
    static StateSetFixture<MODEL> theStateSetFixture;
    return theStateSetFixture;
  }

  StateSetFixture<MODEL>() {
    TopTestParameters<MODEL> parameters;
    parameters.validateAndDeserialize(TestEnvironment::config());

    test_ = std::make_unique<StateSetTestParameters_>(parameters.stateSetTest);
    resol_ = std::make_unique<Geometry_>(parameters.geometry,
                                         oops::mpi::world(), oops::mpi::myself());
  }

  ~StateSetFixture<MODEL>() {}

  std::unique_ptr<StateSetTestParameters_> test_;
  std::unique_ptr<Geometry_>            resol_;
};

// -----------------------------------------------------------------------------
/// \brief tests constructors
template <typename MODEL> void testStateSetConstructors() {
  typedef StateSetFixture<MODEL>     Test_;
  typedef oops::State<MODEL>      State_;
  typedef oops::StateSet<MODEL>   StateSet_;

  const util::DateTime vt(Test_::test().date);

  // Test constructor from individual states
  std::unique_ptr<State_> state1(new State_(Test_::resol(), Test_::test().statefile1));
  std::unique_ptr<State_> state2(new State_(Test_::resol(), Test_::test().statefile2));
  
  std::vector<State_> states;
  states.push_back(*state1);
  states.push_back(*state2);

  // Test main constructor
  std::vector<util::DateTime> times = {vt, vt};
  std::unique_ptr<StateSet_> ss1(new StateSet_(Test_::resol(), state1->variables(), times, oops::mpi::world()));
  EXPECT(ss1.get());
  EXPECT(ss1->size() == 2);
  oops::Log::test() << "Printing StateSet: " << *ss1 << std::endl;

  // Test copy constructor
  std::unique_ptr<StateSet_> ss2(new StateSet_(*ss1));
  EXPECT(ss2.get());
  EXPECT(ss2->size() == ss1->size());

  // Destruct copy
  ss2.reset();
  EXPECT(!ss2.get());

  // Test empty constructor
//  StateSet_ ss3;
//  EXPECT(ss3.size() == 0);
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
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("interface/StateSet/testStateSetConstructors")
      { testStateSetConstructors<MODEL>(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_STATESET_H_
