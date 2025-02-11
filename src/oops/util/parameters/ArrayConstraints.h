/*
 * (C) Copyright 2024 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_PARAMETERS_ARRAYCONSTRAINTS_H_
#define OOPS_UTIL_PARAMETERS_ARRAYCONSTRAINTS_H_

#include <exception>
#include <memory>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "oops/util/parameters/ParameterConstraint.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

/// \brief Constrains a parameter whose JSON representation is an array to contain no less than
/// the specified number of items.
///
/// \tparam T Type of the constrained parameter.
template <typename T>
class MinItemsConstraint : public ParameterConstraint<T> {
 public:
  explicit MinItemsConstraint(size_t minItems) : minItems_(minItems) {}

  void checkValue(const std::string &path, const T &value) const override {
    if (std::size(value) < minItems_)
      throw eckit::BadValue(path + ": Array length is less than " + std::to_string(minItems_),
                            Here());
  }

  PropertyJsonSchema jsonSchema() const override {
    return {{"minItems", std::to_string(minItems_)}};
  }

 private:
  size_t minItems_;
};

template <typename T>
std::shared_ptr<MinItemsConstraint<T>> minItemsConstraint(size_t minItems) {
  return std::make_shared<MinItemsConstraint<T>>(minItems);
}

template <typename T>
std::shared_ptr<MinItemsConstraint<T>> nonEmptyConstraint() {
  return std::make_shared<MinItemsConstraint<T>>(1);
}

/// \brief Constrains a parameter whose JSON representation is an array to contain no more than
/// the specified number of items.
///
/// \tparam T Type of the constrained parameter.
template <typename T>
class MaxItemsConstraint : public ParameterConstraint<T> {
 public:
  explicit MaxItemsConstraint(size_t maxItems) : maxItems_(maxItems) {}

  void checkValue(const std::string &path, const T &value) const override {
    if (std::size(value) > maxItems_)
      throw eckit::BadValue(path + ": Array length is greater than " + std::to_string(maxItems_),
                            Here());
  }

  PropertyJsonSchema jsonSchema() const override {
    return {{"maxItems", std::to_string(maxItems_)}};
  }

 private:
  size_t maxItems_;
};

template <typename T>
std::shared_ptr<MaxItemsConstraint<T>> maxItemsConstraint(size_t maxItems) {
  return std::make_shared<MaxItemsConstraint<T>>(maxItems);
}

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_ARRAYCONSTRAINTS_H_
