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
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/DataSetBase.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------

template<typename MODEL>
class StateSet : public DataSetBase< State<MODEL>, Geometry<MODEL> > {
  typedef Geometry<MODEL>                 Geometry_;
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
  StateSet(const Geometry_ &, const StateSet &);
  StateSet(const StateSet &) = default;
  StateSet(const StateSet &, const int);
  // Calculate the ensemble mean and return a new StateSet variable
  StateSet ens_mean() const;
  // Collect distributed states and return a local subset
  std::unique_ptr<StateSet> localize(const eckit::mpi::Comm &, const Geometry_ &, const int &,
       const int &) const;
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
//  size_t mytime = this->local_time_size() * commTime.rank();
  util::DateTime localtime = times[0];
  std::cout << "MYDBG setting time to be " << localtime << " in ctr " << std::endl;
  for (size_t jm = 0; jm < this->local_ens_size(); ++jm) {
    for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
      this->dataset().emplace_back(new State_(resol, vars, times[jt]));
    }
  }
  this->sync_times();
  this->check_consistency();
  std::cout << "MYDBG time is " << this->times()[0] << " at end of ctr " << std::endl;
  Log::trace() << "StateSet::StateSet" << std::endl;
}


// -----------------------------------------------------------------------------
template<typename MODEL>
StateSet<MODEL>::StateSet(const Geometry_ & resol, const eckit::Configuration & config,
                          const eckit::mpi::Comm & commTime, const eckit::mpi::Comm & commEns)
  : DataSetBase<State_, Geometry_>(commTime, commEns)
{
  Log::trace() << "StateSet::StateSet read start " << config << std::endl;

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
StateSet<MODEL>::StateSet(const StateSet & other, const int ensNum)
  : DataSetBase<State_, Geometry_>(other.times(), other.commTime(),
                                   other.members(), other.commEns())
{
  std::vector<double> zz;
  size_t indx = 0;
  for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
    this->dataset().emplace_back(new State_(other.geometry(), other(jt, ensNum)));
  }
}
// -----------------------------------------------------------------------------

template<typename MODEL>
std::unique_ptr<StateSet<MODEL> > StateSet<MODEL>::localize(const eckit::mpi::Comm & global,
           const Geometry_ & DAgeometry, const int & mytask, const int & ensNum) const
{
/* This method collects parts of the distributed StateSet and places all ensemble
   member states in a smaller patch (1/N the size of Forecast geometry) of a StateSet 
   held in the local_ensemble. It is essentially a transpose of a distributed StateSet
   to a locally held StateSet. The DAgeometry should be have a decomposition that is
   spread across N (number of ensemble members) times the number of MPI tasks that the 
   forecast geometry decomposition. In other words, if the forecast geometry has a 
   layout of [4,4] and there are 9 ensemble members, the DA geometry should have a 
   layout that multiplies to 4*4*9 or something like 12,12. Note that the resolution
   of both geometries is the same (e.g. C48, C96, etc.). Just the decomposition
   is different between the geometries.
*/
  int ist_fc, iend_fc, jst_fc, jend_fc, kst_fc, kend_fc, npz_fc;
  int ist_da, iend_da, jst_da, jend_da, kst_da, kend_da, npz_da;
  int ist_rcv, iend_rcv, jst_rcv, jend_rcv, kst_rcv, kend_rcv, npz_rcv;
  std::unique_ptr<StateSet<MODEL> > local;
  std::vector<int> local_ens;
  size_t dataSize = (*this)(0, 0).serialSize()-3;  // would be good to make this a method
  std::vector<double> zz;
  for (int i = 1; i <= this->ens_size(); ++i) { local_ens.push_back(i); }

  local = std::unique_ptr<StateSet<MODEL> >(new StateSet(DAgeometry, this->variables(), this->times(),
                  this->commTime(), local_ens, oops::mpi::myself()));
  std::cout << "MYDBG in localize, the time is " << this->times()[0] << std::endl;
  std::cout << "MYDBG in localize, the time local has times " << local->times()[0] << std::endl;
  std::vector<int> buf(11);
  std::vector<int> recipients;  // This will contain list of mpi tasks where local tile will be sent
  std::vector<int> senders;  // This will contain list of mpi tasks which will be sending data to me
  std::vector<int> tileEnsNum;  // This will contain list of ensemble numbers that I am receiving
  int mytile = this->geometry().tileNum();
  std::vector<int> global_indices = this->geometry().get_indices();  // pull from this geom and put
                                                                     // into DAgeometry
  ist_fc = global_indices[0];
  iend_fc = global_indices[1];
  jst_fc = global_indices[2];
  jend_fc = global_indices[3];
  kst_fc = global_indices[4];
  kend_fc = global_indices[5];
  npz_fc = global_indices[6];
  int nxg = iend_fc - ist_fc + 1;
  int nyg = jend_fc - jst_fc + 1;
  int nvars = this->variables().size();  // number of variable state

  std::vector<int> indices = DAgeometry.get_indices();
  ist_da = indices[0];
  iend_da = indices[1];
  jst_da = indices[2];
  jend_da = indices[3];
  kst_da = indices[4];
  kend_da = indices[5];
  npz_da = indices[6];

//std::cout << "geom vars broadcasting are " << ist_fc << ", " << iend_fc << ", " << jst_fc << ", " << jend_fc << ", " << std::endl;
//std::cout << "geom vars needed are " << ist_da << ", " << iend_da << ", " << jst_da << ", " << jend_da << ", " << std::endl;

  for (int i = 0; i < global.size(); ++i) {
    if (i == mytask) {  // mytask is global rank
      buf[0] = mytile;   // The tile number that this rank holds
      buf[1] = DAgeometry.tileNum();  // The tile number that I need
      buf[2] = ensNum;   // the ensemble number this tile belongs to
      buf[3] = ist_fc;   // the start of my broadcast domain decomp in i
      buf[4] = iend_fc;  // the start of my broadcast domain decomp in i
      buf[5] = jst_fc;   // the start of my broadcast domain decomp in j
      buf[6] = jend_fc;  // the start of my broadcast domain decomp in j
      buf[7] = ist_da;   // the start of my i domain decomp I NEED
      buf[8] = iend_da;  // the start of my i domain decomp I NEED
      buf[9] = jst_da;   // the start of my j domain decomp I NEED
      buf[10] = jend_da;  // the start of my j domain decomp I NEED
    }
    global.broadcast(buf, i);                 // This is to figure out who is sending domain I NEED
    if ((buf[0] == DAgeometry.tileNum()) &&   // *_fc indices will have a larger span than *_da indices
      ((buf[3] <= ist_da) && (iend_da <= buf[4] )) &&  // *_da indices must be within *_fc indices
      ((buf[5] <= jst_da) && (jend_da <= buf[6] ))) {  //  if the tile, ist, and jst that the sender
                                              // has matches what I need, this is one of my senders
      senders.push_back(i);
      ist_rcv = buf[3];
      iend_rcv = buf[4];
      jst_rcv = buf[5];
      jend_rcv = buf[6];
//    std::cout << "proc " << i << " is sending " << buf[3] << ", " << buf[4] << ", " << buf[5] << ", " << buf[6] << ", " << std::endl;
      tileEnsNum.push_back(buf[2]);
    }
    if ((buf[1] == mytile) &&   // buf here contains indices of domain that is NEEDED by the other processor
       ((ist_fc <= buf[7]) && (buf[8] <= iend_fc)) &&  // NEEDED domain must be within my indices
       ((jst_fc <= buf[9]) && (buf[10] <= jend_fc)) ) {  //  if the DAgeometryetry tile needed
                                     // matches the tile I have, this is who I will send it to
//    std::cout << "proc " << i << " will recive " << ist_fc << ", " << iend_fc<< ", " << jst_fc << ", " <<jend_fc << ", " << std::endl;
      recipients.push_back(i);
    }
  }

// ---- now  send and collect messages
  std::vector<eckit::mpi::Request> send_req_;
  std::vector<eckit::mpi::Request> recv_req_;
  std::vector<size_t> recv_tasks_;
  size_t indx = 0;
  std::vector<double> yy;

  (*this)(0, 0).serialize(zz);  // serialize the forecast state in time 0 and local_ens_number 0
  std::vector<std::vector<double> > zz_recv(senders.size());  // vector of vectors to receive
                                                              // each buffer
  for ( int j = 0; j < recipients.size(); ++j ) {  // reserve space for receiving
     for ( int k = 0; k < zz.size(); ++k ) {  // fill up recv buffers with zeros
       zz_recv[j].push_back(0.0);
     }
  }

  for ( int j = 0; j < recipients.size(); ++j ) {  // loop through list of rcpts/sndrs and send/recv
    if (recipients[j] != mytask) {  // dont send anything to myself
      send_req_.push_back(global.iSend(&zz.front(), zz.size(), recipients[j], ensNum));
    }
    if (senders[j] != mytask) {  // dont need to receive from myself
        recv_req_.push_back(global.iReceive(&zz_recv[j][0], zz.size(), senders[j], tileEnsNum[j]));
        recv_tasks_.push_back(tileEnsNum[j]);
    } else {  // I already have this forecast state
      // copy from my local version
      size_t itask = ensNum-1;
      zz_recv[itask] = zz;
      indx = 0;
  //    int size_fld = (*local)(0, 0).serialSize() - 3;  // get the serialsize of the local tile
      int size_fld = zz_recv[itask].size();  // get the serialsize of the local tile
      (*local)(0, itask).deserializeSection(zz_recv[itask], size_fld, ist_rcv, iend_rcv,
              jst_rcv, jend_rcv, ist_da, iend_da, jst_da, jend_da, indx);  // deserialize state section
    }
  }

// Start looking for messages
  for (size_t r = 0; r < recv_req_.size(); ++r) {
    int ireq = -1;
    eckit::mpi::Status rst = global.waitAny(recv_req_, ireq);
    ASSERT(rst.error() == 0);
    size_t itask = recv_tasks_[ireq] - 1;
    indx = 0;
    //int size_fld = (*local)(0, 0).serialSize() - 3;  // get the serialsize of the local tile
    int size_fld = zz_recv[itask].size();  // get the serialsize of the local tile
    (*local)(0, itask).deserializeSection(zz_recv[itask], size_fld, ist_rcv, iend_rcv,
             jst_rcv, jend_rcv, ist_da, iend_da, jst_da, jend_da, indx);  // deserialize state section
  }
//  (*local).times()[0] = this->times()[0];
//  const std::vector<util::DateTime> times = (*local).validTimes();
//  std::cout << "MYDBG times[0] at end of localize is " << times[0] << std::endl;
  oops::mpi::world().barrier();
  local->sync_times();
  std::cout << "MYDBG times[0] after sync_times is " << local->times()[0] << std::endl;
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
      for (int i = 0; i < dataSize; ++i) {
          zz[0][i] += zz[jm][i]; }
    }
    if (this->local_ens_size() != this->ens_size()) {
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
