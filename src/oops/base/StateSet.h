/*
 * (C) Copyright 2023 UCAR
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/DataSetBase.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/interface/GeometryIterator.h"   
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template<typename MODEL>
class StateSet : public DataSetBase< State<MODEL>, Geometry<MODEL> > {
  typedef Geometry<MODEL>                 Geometry_;
  typedef GeometryIterator<MODEL>         GeometryIterator_;
  typedef Increment<MODEL>                Increment_;
  typedef State<MODEL>                    State_;

 public:
  StateSet(const Geometry_ &, const Variables &,
           const std::vector<util::DateTime> &, const eckit::mpi::Comm &,
           const std::vector<int> & ens = {0},
           const eckit::mpi::Comm & commEns = oops::mpi::myself());
  StateSet(const Geometry_ &, const eckit::Configuration &,
           const eckit::mpi::Comm & commTime = oops::mpi::myself(),
           const eckit::mpi::Comm & commEns = oops::mpi::myself());
  // create a StateSet variable from a std::vector of State variables distributed
  // across communicators
  StateSet(const Geometry_ &, const StateSet &, const int);
  StateSet(const Geometry_ &, const StateSet &);
  StateSet(const StateSet &) = default;
  // Calculate the ensemble mean and return a new StateSet variable
  StateSet ens_mean() const;
  // Collect distributed states and return a local subset
  std::unique_ptr<StateSet> get_local(const eckit::mpi::Comm &, const Geometry_ &, const int &, const int &) const;
  /// Zero
  void zero();
  /// Accumulator
  void accumul(const double &, const StateSet &);
  virtual ~StateSet() = default;

 private:
  std::string classname() const {return "StateSet";}
};

// -----------------------------------------------------------------------------

template<typename MODEL>
StateSet<MODEL>::StateSet(const Geometry_ & resol,
                          const Variables & vars,
                          const std::vector<util::DateTime> & times,
                          const eckit::mpi::Comm & commTime,
                          const std::vector<int> & ens,
                          const eckit::mpi::Comm & commEns)
  : DataSetBase<State_, Geometry_>(times, commTime, ens, commEns)
{
  size_t mytime = this->local_time_size() * commTime.rank();
  for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
    for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
      this->dataset().emplace_back(new State_(resol, vars, times[mytime + jt]));
    }
  }
  this->check_consistency();
  Log::trace() << "StateSet::StateSet" << std::endl;
  Log::info() << "StateSet::StateSet done" << std::endl;
}


// -----------------------------------------------------------------------------
template<typename MODEL>
StateSet<MODEL>::StateSet(const Geometry_ & resol, const eckit::Configuration & config,
                          const eckit::mpi::Comm & commTime, const eckit::mpi::Comm & commEns)
  : DataSetBase<State_, Geometry_>(commTime, commEns)
{
  Log::trace() << "StateSet::StateSet read start " << config << std::endl;
  Log::info() << "StateSet::StateSet read start " << config << std::endl;

// get vector of local configurations
  std::vector<eckit::LocalConfiguration> locals = this->configure(config);

  size_t indx = 0;
  for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
    for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
      this->dataset().emplace_back(new State_(resol, locals.at(indx)));
      ++indx;
    }
  }

  this->sync_times();
  this->check_consistency();

  Log::trace() << "StateSet::StateSet read done" << std::endl;
  Log::info() << "StateSet::StateSet read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StateSet<MODEL>::StateSet(const Geometry_ & resol, const StateSet & other)
  : DataSetBase<State_, Geometry_>(other.times(), other.commTime(),
                                   other.members(), other.commEns())
{
  Log::trace() << "StateSet::StateSet chres start" << std::endl;
  for (size_t jj = 0; jj < other.size(); ++jj) {
    this->dataset().emplace_back(std::make_unique<State_>(resol, other[jj]));
  }
  Log::trace() << "StateSet::StateSet chres done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StateSet<MODEL>::StateSet(const Geometry_ & resol, const StateSet & other, 
     const int local_ens_size )
  : DataSetBase<State_, Geometry_>(other.commTime(), oops::mpi::myself())
{
  Log::trace() << "StateSet::StateSet redist start" << std::endl;
/*
  for (size_t jj = 0; jj < other.size(); ++jj) {
    this->dataset().emplace_back(std::make_unique<State_>(resol, other[jj]));
  }
*/
  Log::trace() << "StateSet::StateSet redist done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
std::unique_ptr<StateSet<MODEL> > StateSet<MODEL>::get_local(const eckit::mpi::Comm & global, 
           const Geometry_ & subgeom, const int & mytask, const int & ensNum) const 
{
  std::unique_ptr<StateSet<MODEL> > local;
  std::vector<int> local_ens;
  size_t dataSize = (*this)(0, 0).serialSize()-3;  // would be good to make this a method
  std::vector<double> zz;
  for(int i = 1; i <= this->ens_size(); ++i) { local_ens.push_back(i); }

  local = std::unique_ptr<StateSet<MODEL> >(new StateSet(subgeom, this->variables(), this->times(), 
                  this->commTime(), local_ens, oops::mpi::myself()));
  std::vector<int> buf(2);
  std::vector<int> global_indices;
  this->geometry().get_indices(global_indices);
  int nxg = global_indices[1] - global_indices[0] + 1;
  int nyg = global_indices[3] - global_indices[2] + 1;
  int nvars = this->variables().size();
  std::cout << "HEY, rank size is " << global.size() << " and mytask is " << mytask << std::endl;
  std::cout << "HEY, nxg is " << nxg << " and nyg " << nyg << std::endl;
  for(int i = 0; i < global.size(); ++i){
    if(i == mytask ) {
      buf[0] = this->geometry().tileNum();  // The tile number about to be sent
      buf[1] = ensNum;
    }
    std::cout << "HEY, starting broadcasts " << i << " on task " << mytask << std::endl;
    global.broadcast(buf, i);
    //broadcast the state
    (*this)(0,0).serialize(zz); //serialize so that we don't overwrite local state
    global.broadcast(zz, i);
    std::cout << "HEY, broadcasts are done on task " << mytask << std::endl;
    //If the incoming tile number matches what this task needs, copy it into local ens
    if(subgeom.tileNum() == buf[0]) { // we need part of the zz buffer
      std::vector<double> yy;
      int ist, iend, jst, jend, npz;
      std::vector<int> indices;
      subgeom.get_indices(indices);
      ist = indices[0];
      iend = indices[1];
      jst = indices[2];
      jend = indices[3];
      std::cout << "calling ssect from " << mytask << " " << ist << "," << iend << "," << jst << "," << jend << std::endl;
      (*local)(0,ensNum-1).serializeSect(yy,ist,iend,jst,jend);
      std::cout << "returned fromssect from " << mytask << " " << yy.size() << std::endl;
      std::cout << "serializing subgeom tileNum " << subgeom.tileNum() << std::endl;
      size_t indx = 0;
      (*local)(0,ensNum-1).deserialize(yy,indx);
      std::cout << "DONE serializing subgeom tileNum " << subgeom.tileNum() << std::endl;
    } else {  
      std::cout << "dont have this tile on " << mytask << " " << std::endl;
    }
#if 0
/*
      (*local)(0,i).serialize(yy);
      int ist, iend, jst, jend, npz;
      std::vector<int> indices;
      subgeom.get_indices(indices);
      ist = indices[0];
      iend = indices[1];
      jst = indices[2];
      jend = indices[3];
      npz = indices[6];
      std::cout << "HEY, ist-- is " << ist <<","<<iend<<","<<jst<<","<<jend<<","<<npz << std::endl;
      std::cout << "HEY, size of yy is " << yy.size() << " zz is " << zz.size() << std::endl;
      // yy is the vector for the subgeom, so we will copy values from zz into it       
      // zz is serialized by variable, zlevels, j, i
      // we will take strides of nvars*zlevels*ny for i values from ist to iend and ny = jend - jst
      int nx = iend - ist + 1;
      int ny = jend - jst + 1;
      int stride = nvars*npz*ny; 
      int index;
      int local_ind = 0;
      for( int v = 0; v < nvars; ++v) {
        for( int kk = 0; kk < npz; ++kk) {
          for( int jj = (jst -1); jj < jend; ++jj) {
            for( int ii = (ist - 1); ii < iend; ++ii) {
              index = v*npz*nyg*nxg + kk * nyg * nxg + jj * nxg + ii; 
              std::cout << mytask << " pushing back index " << local_ind << " " << index <<" "<< ii << " " << jj << " " << kk << " " << v << std::endl;
              yy[local_ind]=zz[index];
              local_ind++;
            }
          }
        }
      }
*/
    }
#endif
//    oops::mpi::world()::barrier();
  }
  return(std::move(local));
}
// -----------------------------------------------------------------------------
template<typename MODEL>
StateSet<MODEL> StateSet<MODEL>::ens_mean() const {
  Log::trace() << "StateSet::ens_mean start" << std::endl;
  StateSet<MODEL> mean = StateSet<MODEL>(this->geometry(), (*this));

  const double fact = 1.0 / static_cast<double>(this->ens_size());

  for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
    size_t dataSize = (*this)(jt, 0).serialSize()-3;  // would be good to make this a method
    // Put States from each local ensemble member in a vector
    std::vector<std::vector<double> > zz(this->local_ens_size());
    for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
      // serialize local ensembles
      (*this)(jt, jm).serialize(zz[jm]);
    }

// add up all the state values on the local communicator and put them in zz[0][:]
    for (size_t jm = 1; jm < this->local_ens_size(); ++jm) {
      for ( int i = 0; i < dataSize; ++i) zz[0][i] += zz[jm][i];
    }
    if (this->commEns().size() > 1) {
      // if commEns > 1, then sum up across commEns communicators
      this->commEns().allReduceInPlace(&(zz[0].front()), dataSize, eckit::mpi::Operation::SUM);
    }
    // Divide by total number of members to get average
    for ( int i = 0; i < dataSize; ++i) zz[0][i] *= fact;

// deserialize back to stateSet
    size_t indx = 0;
    mean[jt].deserialize(zz[0], indx);
  }
  Log::trace() << "StateSet::ens_mean done" << std::endl;
  return mean;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void StateSet<MODEL>::zero() {
  Log::trace() << "StateSet<MODEL>::zero starting" << std::endl;
  for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
    for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
      (*this)(jt, jm).zero();
    }
  }
  Log::trace() << "StateSet<MODEL>::zero done" << std::endl;
}
// -----------------------------------------------------------------------------

template<typename MODEL>
void StateSet<MODEL>::accumul(const double & zz, const StateSet & xx) {
  Log::trace() << "StateSet<MODEL>::accumul starting" << std::endl;
  for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
    for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
      (*this)(jt, jm).accumul(zz, xx(jt, jm));
    }
  }
  Log::trace() << "StateSet<MODEL>::accumul done" << std::endl;
}

}  // namespace oops
