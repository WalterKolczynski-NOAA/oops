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
    const int nEns = 4;  // Number of ensemble members
    
    // Create DA geometry with nx=4, ny=4 layout (16 total points)
    eckit::LocalConfiguration daGeomConfig;
    daGeomConfig.set("nx", 4);
    daGeomConfig.set("ny", 4);
    
    // Create FC geometry with ax=2, ay=2 layout (4 points, 16/nEns)
    eckit::LocalConfiguration fcGeomConfig;
    fcGeomConfig.set("nx", 2); 
    fcGeomConfig.set("ny", 2);

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
    Geometry<MODEL> daGeom(daGeomConfig, worldComm);
    Geometry<MODEL> fcGeom(fcGeomConfig, commMember);

    // Create variables and times for StateSet
    Variables vars;
    vars.push_back(Variable("var1"));
    std::vector<util::DateTime> times;
    times.push_back(util::DateTime("2021-01-01T00:00:00Z"));

    // Create ensemble members vector
    std::vector<int> ensMembers;
    for (int i = 1; i <= nEns; i++) {
      ensMembers.push_back(i);
    }

    // Create StateSet for DA geometry
    StateSet<MODEL> daStateSet(daGeom, vars, times, oops::mpi::myself(), 
                              ensMembers, patchMember);

    // Create StateSet for FC geometry
    StateSet<MODEL> fcStateSet(fcGeom, vars, times, oops::mpi::myself(),
                              ensMembers, patchMember);

    // Test transpose functionality between geometries
    std::vector<StateSet<MODEL>> localVec = fcStateSet.transpose(worldComm, daGeom, mymember);

    // Verify dimensions
    EXPECT(localVec.size() == nEns);
    EXPECT(daStateSet.size() == nEns * times.size());
    EXPECT(fcStateSet.size() == times.size());
  }

  const eckit::Configuration & config_;
};

}  // namespace test
}  // namespace oops

#endif  // TEST_BASE_STATESET_H_TEST