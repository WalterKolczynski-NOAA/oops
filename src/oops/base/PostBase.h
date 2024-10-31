/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_POSTBASE_H_
#define OOPS_BASE_POSTBASE_H_

#include <boost/noncopyable.hpp>

#include "oops/base/PostTimer.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace eckit {
  class Configuration;
}

namespace oops {

class PostTimerParameters;

// -----------------------------------------------------------------------------

/// Handles post-processing of model fields.
/*!
 *  PostBase is the base class for all state post processors, it
 *  is mostly used so that PostProcessor can hold a vector of such
 *  processors.
 *  By default processing is performed on every call.
 */

template <typename FLDS> class PostBase : private boost::noncopyable {
 public:
/// Constructors and basic operators
  PostBase() : timer_() {}
  explicit PostBase(const eckit::Configuration & config) : timer_(config) {}
  explicit PostBase(const PostTimerParameters & timerParams) : timer_(timerParams) {}
  PostBase(const util::DateTime & start, const util::DateTime & finish,
           const util::Duration & freq = util::Duration(0))
    : timer_(start, finish, freq) {}

  virtual ~PostBase() {}

/// Setup
  void initialize(const FLDS & xx, const util::DateTime & end,
                  const util::Duration & tstep) {
    std::cout << "in Postbase, initializing timer_ " << xx.validTime() << std::endl;
    std::cout << "in Postbase, initializing timer_ end is " << end << std::endl;
    timer_.initialize(xx.validTime(), end);
//    std::cout << "in Postbase, timer(validtime) is " << timer_.itIsTime(xx.validTime()) << std::endl;
    this->doInitialize(xx, end, tstep);
  }

/// Process state or increment
  void process(const FLDS & xx) {
    std::cout << "in Postbase process, validtime is " << xx.validTime() << std::endl;
//    std::cout << "in Postbase process, timer(validtime) is " << timer_.itIsTime(xx.validTime()) << std::endl;
    if (timer_.itIsTime(xx.validTime())) {
        std::cout << "in Postbase, xx is " << xx << std::endl;
        this->doProcessing(xx);  
    } else {
      std::cout << "in Postbase, DID NOT PROCESS " << std::endl;
    }
    std::cout << "in Postbase, leaving process " << xx.validTime() << std::endl;
  }

/// Final
  void finalize(const FLDS & xx) {
    std::cout << "in Postbase, starting finalize" << xx.validTime() << std::endl;
    this->doFinalize(xx);
    std::cout << "in Postbase, done finalize" << xx.validTime() << std::endl;
  }

 private:
  PostTimer timer_;

/// Actual processing
  virtual void doProcessing(const FLDS &) = 0;
  virtual void doInitialize(const FLDS &, const util::DateTime &,
                            const util::Duration &) {}
  virtual void doFinalize(const FLDS &) {}
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_POSTBASE_H_
