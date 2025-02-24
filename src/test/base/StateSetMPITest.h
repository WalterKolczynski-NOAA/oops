#ifndef TEST_BASE_STATESET_H_TEST
#define TEST_BASE_STATESET_H_TEST

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/base/StateSet.h"
#include "oops/base/Geometry.h"
#include "oops/util/DateTime.h"
#include "test/TestEnvironment.h"

namespace oops {
namespace test {

template <typename MODEL>
class StateSetMPITest : public oops::Application {
  typedef Geometry<MODEL>  Geometry_;
//  typedef Model<MODEL>  Model_;
  typedef State<MODEL>  State_;
  typedef StateSet<MODEL>  StateSet_;
 public:
  explicit StateSetMPITest(const eckit::Configuration & config) 
    : Application(oops::mpi::world()), config_(config) {}
  virtual ~StateSetMPITest() = default;

  // Required override for Application
  std::string appname() const override { return "StateSetMPITest"; }

  // Required override for Application
  int execute(const eckit::Configuration & config) const override {
    std::cout << "starting testDiffGeom " << std::endl;
    testDifferentGeometries();
    return 0;
  }

 private:
  void testDifferentGeometries() const {
    const int nEns = 2;  // Number of ensemble members
    
    std::cout << "testDiffGeom 1 " << std::endl;
    // DA geom uses all mpi tasks on MPI_COMM_WORLD
    // Get DA geometry configuration from config_
    eckit::LocalConfiguration daGeomConfig = config_.getSubConfiguration("da geometry");
    // Get full configuration
//    config_.get("da geometry", daGeomConfig);

    //FC geometry is created on each ensemble communicator
    std::cout << "testDiffGeom 2 " << std::endl;
    // Get FC geometry configuration from config_
    eckit::LocalConfiguration fcGeomConfig;
    config_.get("fc geometry", fcGeomConfig);

    std::cout << "testDiffGeom 3 " << std::endl;
    // Create DA communicator
    const eckit::mpi::Comm & worldComm = oops::mpi::world();
    const int ntasks = worldComm.size();
    const int mytask = worldComm.rank();
    const int tasks_per_member = ntasks / nEns;
    const int mymember = mytask / tasks_per_member + 1;

    std::cout << "testDiffGeom 4 -- " << ntasks << " " << nEns << " " << tasks_per_member << " " << mymember << " " << mytask << " " << std::endl;
    // Create member communicator
    std::string commNameStr = "comm_member_" + std::to_string(mymember);
    char const *commName = commNameStr.c_str();
    eckit::mpi::Comm & commMember = worldComm.split(mymember, commName);
    const int subrank = commMember.rank();

    std::cout << "my rank is " << mytask << " and subrank is " << subrank << std::endl;
    // Create patch communicator
    std::string patchNameStr = "patch_member_" + std::to_string(subrank);
    char const *patchName = patchNameStr.c_str();
    eckit::mpi::Comm & patchMember = worldComm.split(subrank, patchName);

    std::cout << "testDiffGeom 6 " << std::endl;
    // Create geometries using appropriate communicators
    std::cout << "daGeom is " << daGeomConfig << std::endl;
    Geometry_ daGeom(daGeomConfig, worldComm);
    std::cout << "daGeom constructed" << std::endl;
    Geometry<MODEL> fcGeom(fcGeomConfig, commMember);
    std::cout << "fcGeom constructed" << std::endl;


    // Create variables and times for StateSet
    Variables vars;
    vars.push_back(Variable("ua"));
    vars.push_back(Variable("va"));
    vars.push_back(Variable("t"));
    vars.push_back(Variable("delp"));
    vars.push_back(Variable("sphum"));
    
    std::vector<util::DateTime> times;
    times.push_back(util::DateTime("2021-03-23T06:00:00Z"));

    // Create ensemble members vector
    std::vector<int> ensMembers;
    for (int i = 1; i <= nEns; i++) {
      ensMembers.push_back(i);
    }

    // Create StateSet for DA geometry
    /*
    StateSet<MODEL> daStateSet(daGeom, vars, times, oops::mpi::myself(), 
                              ensMembers, patchMember);
    std::cout << "about to randomize daStatSet" << std::endl;
    daStateSet.random();
    */
    // Create StateSet for FC geometry
    StateSet<MODEL> fcStateSet(fcGeom, vars, times, oops::mpi::myself(),
                              ensMembers, patchMember);

    fcStateSet.random();
    std::cout << "fc stateset size is " << fcStateSet.size() << std::endl;
    std::cout << "fc stateset is " << fcStateSet << std::endl;
    // Test transpose functionality between geometries
    std::vector<StateSet<MODEL>> localVec = fcStateSet.transpose(worldComm, daGeom, mymember);
    oops::mpi::world().barrier();
    /*
    for(size_t ii=0; ii < localVec.size(); ++ii) {
	std::cout << "localVec[" << ii << "] is " << localVec[ii] << std::endl;
    }
    */
    /*

    // Verify dimensions
    EXPECT(localVec.size() == nEns);
    EXPECT(daStateSet.size() == nEns * times.size());
    EXPECT(fcStateSet.size() == times.size());
    */
  }

  const eckit::Configuration & config_;
};

}  // namespace test
}  // namespace oops

#endif  // TEST_BASE_STATESET_H_TEST
