#ifndef TEST_BASE_STATESET_H_TEST
#define TEST_BASE_STATESET_H_TEST

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/StateSet.h"
#include "oops/base/Geometry.h"
#include "oops/util/DateTime.h"
#include "test/TestEnvironment.h"

namespace oops {
//namespace test {

template <typename MODEL>
class StateSetMPITest : public oops::Application {
  typedef Geometry<MODEL>  Geometry_;
//  typedef Model<MODEL>  Model_;
  typedef State<MODEL>  State_;
  typedef StateSet<MODEL>  StateSet_;
  typedef IncrementSet<MODEL>  IncrementSet_;
//  const eckit::Configuration & config2 = test::TestEnvironment::config();
 public:
  explicit StateSetMPITest(const eckit::Configuration & config) 
    : Application(oops::mpi::world()), config_(config) {}
  virtual ~StateSetMPITest() = default;

  // Required override for Application
  std::string appname() const override { return "StateSetMPITest"; }

  // Required override for Application
  int execute(const eckit::Configuration & config) const override {
    testDifferentGeometries();
    return 0;
  }

 private:
  void testDifferentGeometries() const {
//    const eckit::Configuration & config = test::TestEnvironment::config();  
    eckit::PathName confPath("testinput/test_stateset.yaml");
    const eckit::YAMLConfiguration config2(confPath);
    int nEns;  // Number of ensemble members
  
    config2.get("ensemble members", nEns); 
    // DA geom uses all mpi tasks on MPI_COMM_WORLD
    // Get DA geometry configuration from config_
    eckit::LocalConfiguration daGeomConfig = config2.getSubConfiguration("da geometry");
    std::cout << "da geom" << daGeomConfig << std::endl;
    // Get full configuration
//    config_.get("da geometry", daGeomConfig);

    //FC geometry is created on each ensemble communicator
    // Get FC geometry configuration from config_
    eckit::LocalConfiguration fcGeomConfig;
    config2.get("fc geometry", fcGeomConfig);
    std::cout << "fc geom" << fcGeomConfig << std::endl;

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

    // Create geometries using appropriate communicators
    Geometry_ daGeom(daGeomConfig, worldComm);
    std::cout << "Hey, myrank is " << mytask << " and mymember is " << mymember << std::endl;
    fcGeomConfig.set("member_number",mymember);
    std::cout << "fcgeomconf is " << fcGeomConfig << std::endl;
    std::cout << "commMember size is " << commMember.size() << std::endl;
    Geometry<MODEL> fcGeom(fcGeomConfig, commMember);

    //  Setup times
    Log::info() << "setting up times" << std::endl;
    eckit::LocalConfiguration fcstparams = config2.getSubConfiguration("fcst");
    eckit::LocalConfiguration model = fcstparams.getSubConfiguration("model");
    const util::Duration tstep(model.getString("tstep"));
    eckit::LocalConfiguration ic = fcstparams.getSubConfiguration("initial condition");
    const util::DateTime bgndate(ic.getString("datetime"));
    const util::Duration fclength(fcstparams.getString("forecast length"));
    const util::DateTime enddate(bgndate + fclength);
    std::vector<util::DateTime> times;
    const Variables vars(ic, "state variables");
    // Create ensemble members vector
    std::vector<int> ensMembers;
    for (int i = 1; i <= nEns; i++) {
      ensMembers.push_back(i);
    }
    for (util::DateTime ii=(bgndate+tstep); ii <= enddate; ii=ii+tstep) {
       Log::info() << "pushing back time " << ii << std::endl;
       times.push_back(ii);
    }
    oops::mpi::world().barrier();

    // Create StateSet for FC geometry
    StateSet<MODEL> fcStateSet(fcGeom, vars, times, oops::mpi::myself(),
                              ensMembers, patchMember);

    fcStateSet.random();
    Log::trace() << "before transpose fc stateset is " << fcStateSet << std::endl;
    // Test transpose functionality between geometries
    // daStateSet is on the daGeom with both ensemble members on each MPI proc
    std::vector<StateSet<MODEL>> daStateSet = fcStateSet.transpose(worldComm, daGeom, mymember);
    oops::mpi::world().barrier();
    Log::trace() << "size of daStateSet is " << daStateSet.size() << std::endl;
    Log::trace() << "localvec[0] is " << daStateSet[0] << std::endl;
    Log::trace() << "localvec[1] is " << daStateSet[1] << std::endl;
    std::vector<State_> states_;
    // after Rtranspose, states are back to distributed across ensemble ranks
    for(size_t ens=0; ens < daStateSet.size(); ++ens){
//      if((mymember - 1) == ens) { //we only want our ensemble member
        states_.emplace_back(((daStateSet[ens]).Rtranspose(this->getComm(), fcGeom,
          mymember,ens)));
//      }
    }
    Log::trace() << "size of states_ is " << states_.size() << std::endl;
    Log::trace() << "state after Rtranspose is " << states_[0] << std::endl;
    Log::trace() << "state[1] after Rtranspose is " << states_[1] << std::endl;
    Log::trace() << "times[0] is " << times[0] << std::endl;
    StateSet_ *newFCStateSet = new StateSet_(states_, mymember - 1, times, oops::mpi::myself(),
		                       ensMembers, patchMember);
    Log::trace() << "newFCState after Rtranspose is " << *newFCStateSet << std::endl;
    IncrementSet_ newState(fcGeom, vars, times, oops::mpi::myself(), ensMembers, patchMember);
    newState.diff(fcStateSet,*newFCStateSet);
    Log::trace() << "diff between stateSets is " << newState << std::endl;
    delete newFCStateSet;
    /*

    // Verify dimensions
    EXPECT(daStateSet.size() == nEns);
    EXPECT(daStateSet.size() == nEns * times.size());
    EXPECT(fcStateSet.size() == times.size());
    */
  }

  const eckit::Configuration & config_;
};

//}  // namespace test
}  // namespace oops

#endif  // TEST_BASE_STATESET_H_TEST
