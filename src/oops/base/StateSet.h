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
  StateSet(const StateSet &, const int );
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

template<typename MODEL>
StateSet<MODEL>::StateSet(const StateSet & other, const int ensNum)
  : DataSetBase<State_, Geometry_>(other.times(), other.commTime(),
                                   other.members(), other.commEns())
{

  std::vector<double> zz;  
  size_t indx = 0;
  for (size_t jt = 0; jt < this->local_time_size(); ++jt) {
    std::cout << "pushing back ens number" << ensNum << std::endl;
//    std::cout << "other(jt,ensNum)" << other(jt,ensNum) << std::endl;
    this->dataset().emplace_back(new State_(other.geometry(), other(jt,ensNum)));
    std::cout << "done pushing back ens number" << ensNum << std::endl;
/*
    other(jt,ensNum).serialize(zz); //serialize 
    std::cout << "zz is " << zz[0] << "," << zz[1] << "," << zz[2] << std::endl;
    (*this)[jt].deserialize(zz,indx);
    std::cout << "done deserializing" << ensNum << std::endl;
*/
/*
*/
  } 
    
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
  std::vector<int> buf(3);
  std::vector<int> recipients; // This will contain list of mpi tasks where local tile should be sent
  std::vector<int> senders; // This will contain list of mpi tasks which will be sending data to me
  std::vector<int> tileEnsNum; // This will contain list of ensemble numbers that I am receiving
  std::vector<int> global_indices;
  int mytile = this->geometry().tileNum();
  this->geometry().get_indices(global_indices); // returns ist, iend, jst, jend, kst, kend, npz
  int nxg = global_indices[1] - global_indices[0] + 1;
  int nyg = global_indices[3] - global_indices[2] + 1;
  int nvars = this->variables().size(); //number of variable state
  std::cout << "HEY, rank size is " << global.size() << " and mytask is " << mytask << std::endl;
  std::cout << "HEY, nxg is " << nxg << " and nyg " << nyg << std::endl;
  std::cout << "before bcast state looks like this " << (*this)(0,0) << std::endl;
  StateSet<MODEL> tmpState = StateSet<MODEL>(*this);  //ths tmpState is full size
  for(int i = 0; i < global.size(); ++i){
    if(i == mytask ) { //mytask is global rank
      buf[0] = mytile;  // The tile number that this rank holds
      buf[1] = subgeom.tileNum();  // The tile number that I need
      buf[2] = ensNum;  // the ensemble number this tile belongs to
    }
    std::cout << "HEY, starting broadcasts " << i << " on task " << mytask << std::endl;
    global.broadcast(buf, i);
    std::cout << "buf on task " << i << " is " << buf[0] << " " << buf[1] << std::endl;
    if(buf[0] == subgeom.tileNum()) { //  if the tile that the sender has matches what I need, this is one of my senders
      senders.push_back(i);
      tileEnsNum.push_back(buf[2]);
    }
    if(buf[1] == mytile) { //  if the subgeometry tile needed matches the tile I have, this is who I will send it to
      recipients.push_back(i);
    }
#if 0
    //broadcast the state
    const std::vector<util::DateTime> times = (*this).validTimes(); 
    std::cout << "before serialize, times are " << times[0] << " and " << times.size() << std::endl;
    (*this)(0,0).serialize(zz);  //serialize the state in time 0 and local_ens_number 0
    global.broadcast(zz, i);  //broadcast from root of rank i
    //If the incoming tile number matches what this task needs, copy it into local ens
    if(subgeom.tileNum() == buf[0]) { // we need part of the zz buffer
      size_t indx = 0;
      tmpState[0].deserialize(zz,indx);  
      Log::info() << "Start of zz is " << zz[0] << " " << zz[10] << std::endl;
      std::cout << "HEY, broadcasts are done on task " << mytask << std::endl;
      int size_fld = (*local)(0,0).serialSize() - 3;  // get the serialsize of the smaller tile (local)
      std::vector<double> yy;
      int ist, iend, jst, jend, npz;
      std::vector<int> indices;
      subgeom.get_indices(indices);
      ist = indices[0];
      iend = indices[1];
      jst = indices[2];
      jend = indices[3];
      std::cout << "calling ssect from " << mytask << " " << ist << "," << iend << "," << jst << "," << jend << std::endl;
      tmpState[0].serializeSect(yy,ist,iend,jst,jend,size_fld);
      Log::info() << "Start of yy is " << yy[0] << " " << yy[10] << std::endl;
      std::cout << "returned fromssect from " << mytask << " " << yy.size() << std::endl;
      std::cout << "serializing subgeom tileNum " << subgeom.tileNum() << "," << ensNum << std::endl;
      std::cout << "Start of yy is " << yy[0] << "," << yy[10] << std::endl;
      indx = 0;
      (*local)(0,ensNum-1).deserialize(yy,indx);
      std::cout << "DONE deserializing subgeom tileNum " << subgeom.tileNum() << std::endl;
    } else {  
      std::cout << "dont have this tile on " << mytask << " " << std::endl;
    }
  std::cout << "writing out tmpState[0]" << tmpState[0] << std::endl;
#endif
  }
  std::cout << "From task " << mytask << " my recipients are " << recipients[0] << " " << recipients[1] << "|||" << recipients.size() << std::endl;
  std::cout << "From task " << mytask << " my senders are " << senders[0] << " " << senders[1] << "|||" << senders.size() << std::endl;
  std::vector<eckit::mpi::Request> send_req_;
  std::vector<eckit::mpi::Request> recv_req_;
  std::vector<size_t> recv_tasks_;
  size_t indx = 0;

  (*this)(0,0).serialize(zz);  //serialize the state in time 0 and local_ens_number 0
  std::vector<std::vector<double> > zz_recv(senders.size());  // vector of vectors to receive each buffer
  for( int j = 0; j < recipients.size(); ++j ) {  // reserve space for receiving
     for( int k = 0; k < zz.size(); ++k ) {
       zz_recv[j].push_back(0.0);
     }
  }
  std::cout << "From task " << mytask << " zz size is " << zz.size() << " and  " << zz_recv[0].size() << "|||" << zz_recv[1].size() << std::endl;
  for( int j = 0; j < recipients.size(); ++j ) {  // loop through list of recipients/senders and send/recv
    if(recipients[j] != mytask) {  // dont send anything to myself
      std::cout << "on mytask " << mytask << " sending to " << recipients[j] << std::endl;
      send_req_.push_back(global.iSend(&zz.front(), zz.size(),recipients[j],ensNum));

    }
    if(senders[j] != mytask) { // dont need to receive from myself
        std::cout << "on mytask " << mytask << " looking for mesg from " << senders[j] << std::endl;
        recv_req_.push_back(global.iReceive(&zz_recv[j][0], zz.size(),senders[j],tileEnsNum[j]));
        recv_tasks_.push_back(tileEnsNum[j]);
        std::cout << "on task " << mytask << " pushing back recv_tasks " << tileEnsNum[j] << std::endl;
    } else { 
      // copy from my local version
      size_t itask = ensNum-1;
      zz_recv[itask] = zz;
      indx = 0;
      std::cout << "size of zz_recv[" << itask << "] is " << zz_recv[itask].size() << std::endl;
      tmpState[0].deserialize(zz_recv[itask],indx); // put the full zz recv vector in a state
      int size_fld = (*local)(0,0).serialSize() - 3;  // get the serialsize of the smaller tile (local)
      std::vector<double> yy;
      int ist, iend, jst, jend, npz;
      std::vector<int> indices;
      subgeom.get_indices(indices);
      ist = indices[0];
      iend = indices[1];
      jst = indices[2];
      jend = indices[3];
      std::cout << "calling ssect from " << mytask << " " << ist << "," << iend << "," << jst << "," << jend << std::endl;
      tmpState[0].serializeSect(yy,ist,iend,jst,jend,size_fld);  // copy out on the part of the state we need
      Log::info() << "Start of yy is " << yy[0] << " " << yy[10] << std::endl;
      std::cout << "returned fromssect from " << mytask << " " << yy.size() << std::endl;
      std::cout << "serializing subgeom tileNum " << subgeom.tileNum() << "," << ensNum << std::endl;
      std::cout << "Start of yy is " << yy[0] << "," << yy[10] << std::endl;
      indx = 0;
      (*local)(0,itask).deserialize(yy,indx); // assuming itask is the tag number, put this state in ens number itask
      std::cout << "DONE deserializing subgeom tileNum " << subgeom.tileNum() << std::endl;
    } 
  }
  std::cout << "starting to look for messages. I need " << recv_req_.size() << std::endl;
//  tmpState[0].deserialize(zz,indx);  
  for( size_t r = 0; r < recv_req_.size(); ++r) {
    int ireq = -1;
    eckit::mpi::Status rst = global.waitAny(recv_req_, ireq);
    std::cout << "On rank " << mytask << " got a message with ireq = " << ireq << " " << recv_tasks_[ireq] << std::endl;
    ASSERT(rst.error() == 0);
    size_t itask = recv_tasks_[ireq] - 1;
      indx = 0;
      std::cout << "size of zz_recv[" << itask << "] is " << zz_recv[itask].size() << std::endl;
      tmpState[0].deserialize(zz_recv[itask],indx); // put the full zz recv vector in a state
      int size_fld = (*local)(0,0).serialSize() - 3;  // get the serialsize of the smaller tile (local)
      std::vector<double> yy;
      int ist, iend, jst, jend, npz;
      std::vector<int> indices;
      subgeom.get_indices(indices);
      ist = indices[0];
      iend = indices[1];
      jst = indices[2];
      jend = indices[3];
      std::cout << "calling ssect from " << mytask << " " << ist << "," << iend << "," << jst << "," << jend << std::endl;
      tmpState[0].serializeSect(yy,ist,iend,jst,jend,size_fld);  // copy out on the part of the state we need
      Log::info() << "Start of yy is " << yy[0] << " " << yy[10] << std::endl;
      std::cout << "returned fromssect from " << mytask << " " << yy.size() << std::endl;
      std::cout << "serializing subgeom tileNum " << subgeom.tileNum() << "," << ensNum << std::endl;
      std::cout << "Start of yy is " << yy[0] << "," << yy[10] << std::endl;
      indx = 0;
      (*local)(0,itask).deserialize(yy,indx); // assuming itask is the tag number, put this state in ens number itask
      std::cout << "DONE deserializing subgeom tileNum " << subgeom.tileNum() << std::endl;

  }
  std::cout << "DONE deserializing subgeom state 0 is " << (*local)(0,0) << std::endl;
  std::cout << "DONE deserializing subgeom state 1 is " << (*local)(0,1) << std::endl;
  const std::vector<util::DateTime> times = (*local).validTimes(); 
  return(std::move(local));
}
// -----------------------------------------------------------------------------
template<typename MODEL>
StateSet<MODEL> StateSet<MODEL>::ens_mean() const {
  Log::trace() << "StateSet::ens_mean start" << std::endl;
  StateSet<MODEL> mean = StateSet<MODEL>(this->geometry(), (*this));

  const double fact = 1.0 / static_cast<double>(this->ens_size());

  Log::info() << "in ens_mean, ens_size, local_ens_size, local_time_size are " << this->ens_size() << " " << this->local_ens_size() << " " << this->local_time_size() << std::endl;
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
//      Log::info() << "in ens_mean, summing local ensmebles 0 and " << jm << " fact is " << fact << std::endl;
      for ( int i = 0; i < dataSize; ++i) { 
//          Log::info() << zz[0][i] << " " << zz[jm][i] << std::endl;
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
