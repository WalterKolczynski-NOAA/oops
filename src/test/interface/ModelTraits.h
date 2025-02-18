/*
 * (C) Copyright 2024 JEDI
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_MODELTRAITS_H_
#define TEST_INTERFACE_MODELTRAITS_H_

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
#include "test/base/Variables.h"

namespace test {

/// Test Model Traits for unit testing
struct TestModelTraits {
  // Basic model trait type definitions
  typedef int                            GeometryValueType;
  typedef double                         StateValueType;
  typedef oops::Geometry<TestModelTraits>     Geometry;
  typedef oops::State<TestModelTraits>        State;
  typedef oops::Variables    Variables;
  typedef oops::GeoVaLs<TestModelTraits>      GeoVaLs;

  static const std::string classname() {return "TestModel";}
};

// Model-specific configuration for test geometry
struct TestGeometryData {
  TestGeometryData() = default;
  TestGeometryData(const eckit::Configuration & config) {
    config.get("nx", nx);
    config.get("ny", ny);
    npoints = nx * ny;
  }
  int nx;
  int ny;
  int npoints;
};

// Model-specific configuration for test state
struct TestStateData {
  TestStateData() = default;
  explicit TestStateData(const TestGeometryData & geom) : values(geom.npoints, 0.0) {}
  std::vector<double> values;
};

}  // namespace test

#endif  // TEST_INTERFACE_MODELTRAITS_H_
