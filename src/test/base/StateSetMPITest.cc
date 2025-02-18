#include "oops/runs/Run.h"
#include "oops/test/base/StateSetMPITest.h"
#include "test/interface/ModelTraits.h"
#include "test/TestEnvironment.h"

int main(int argc, char ** argv) {
  oops::Run run(argc, argv);
  eckit::LocalConfiguration config;
  test::TestEnvironment::getInstance().setup(config);
  oops::test::StateSetMPITest<test::TestModelTraits> test(config);
  return run.execute(test);
}