/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_INCREMENTENSEMBLESET_H_
#define OOPS_BASE_INCREMENTENSEMBLESET_H_

#include <Eigen/Dense>
#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/LocalIncrement.h"
#include "oops/base/StateSet.h"
#include "oops/base/StateEnsembleSet.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Ensemble of 4D increments
template<typename MODEL> class IncrementEnsembleSet {
  typedef Geometry<MODEL>            Geometry_;
  typedef GeometryIterator<MODEL>    GeometryIterator_;
  typedef StateSet<MODEL>             StateSet_;
  typedef StateEnsembleSet<MODEL>     StateEnsembleSet_;
  typedef IncrementSet<MODEL>         IncrementSet_;

 public:
  /// Constructor
  IncrementEnsembleSet(const Geometry_ & resol,
                      const Variables & vars,
                      const std::vector<util::DateTime> &,
                      const int rank);
  /// \brief construct ensemble of perturbations as \p ens - \p mean; holding
  //         \p vars variables
  IncrementEnsembleSet(const StateEnsembleSet_ & ens, const StateSet_ & mean,
                      const Variables & vars);

  IncrementEnsembleSet(const StateSet_ & ensemble,
                      const StateSet_ & mean, 
                      const Geometry_ & resol,
                      const Variables & vars);
  /// Accessors
  size_t size() const {return ensemblePerturbs_.size();}
  IncrementSet_ incrementSet() const {return ensemblePerturbsSet_;}

 private:
  IncrementSet_ ensemblePerturbsSet_;
};

// ====================================================================================

template<typename MODEL>
IncrementEnsembleSet<MODEL>::IncrementEnsembleSet(const Geometry_ & resol, const Variables & vars,
                                                const std::vector<util::DateTime> & timeslots,
                                                const int rank)
  : ensemblePerturbsSet_(resol, vars, timeslots, oops::mpi::myself())
{
  Log::trace() << "IncrementEnsembleSet:contructor done" << std::endl;
}

// ====================================================================================
template<typename MODEL>
IncrementEnsembleSet<MODEL>::IncrementEnsembleSet(const StateSet_ & ensemble,
                                                const StateSet_ & mean, 
                                                const Geometry_ & resol,
                                                const Variables & vars)
  : ensemblePerturbsSet_(resol, vars, ensemble)     
{
  ensemblePerturbsSet_.diff(ensemble,mean);
  Log::trace() << "IncrementEnsembleSet:contructor(StateEnsembleSet) done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
IncrementEnsembleSet<MODEL>::IncrementEnsembleSet(const StateEnsembleSet_ & ensemble,
                                                const StateSet_ & mean, const Variables & vars)
  : ensemblePerturbsSet_(mean[0].geometry(), vars, mean.times(), oops::mpi::myself())
{
  Log::trace() << "IncrementEnsembleSet:contructor(StateEnsembleSet) done" << std::endl;
}

}  // namespace oops

#endif  // OOPS_BASE_INCREMENTENSEMBLESET_H_
