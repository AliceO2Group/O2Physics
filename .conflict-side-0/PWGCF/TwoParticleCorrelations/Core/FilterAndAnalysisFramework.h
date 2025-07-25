// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#ifndef PWGCF_TWOPARTICLECORRELATIONS_CORE_FILTERANDANALYSISFRAMEWORK_H_
#define PWGCF_TWOPARTICLECORRELATIONS_CORE_FILTERANDANALYSISFRAMEWORK_H_

#include "PWGCF/TwoParticleCorrelations/Core/EventSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/PIDSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/TrackSelectionFilterAndAnalysis.h"

#include <CCDB/BasicCCDBManager.h>

#include <TList.h>
#include <TNamed.h>
#include <TObject.h>
#include <TString.h>

#include <Rtypes.h>

#include <string>
#include <vector>

namespace o2
{
namespace analysis
{
namespace PWGCF
{

/// \brief Base class for filter and selection once filetered
class FilterAndAnalysisFramework : public TNamed
{
  friend void registerConfiguration(FilterAndAnalysisFramework*);

 public:
  FilterAndAnalysisFramework();
  FilterAndAnalysisFramework(std::string& url, std::string& path, std::string& date, bool force = false)
    : ccdburl(url), ccdbpath(path), filterdate(date), forceupdate(force) {}

  void SetConfiguration(EventSelectionConfigurable&, TrackSelectionConfigurable&, PIDSelectionConfigurable&, SelectionFilterAndAnalysis::selmodes mode);

 protected:
  void RegisterConfiguration();

 public:
  void Init();

  /// \brief get the valid (armed) mask associated to the configuration of the different objects
  /// \return the corresponding armed mask
  uint64_t getCollisionMask() { return fEventFilter->getMask(); }
  uint64_t getTrackMask() { return fTrackFilter->getMask(); }
  uint64_t getPIDMask() { return fPIDFilter->getMask(); }
  /// \brief get the valid (armed) optional part mask associated to the configuration of the different objects
  /// \return the corresponding armed optional mask
  /// A clear example of the optional part mask is the mask of the multiplicity classes
  /// where only one of the available in the whole mask will be flagged
  std::vector<uint64_t>& getCollisionOptMask() { return fEventFilter->getOptMask(); }
  std::vector<uint64_t>& getTrackOptMask() { return fTrackFilter->getOptMask(); }
  std::vector<uint64_t>& getPIDOptMask() { return fPIDFilter->getOptMask(); }
  /// \brief provides the optional part mask of the different objects in a printable shape
  TString printCollisionOptionalMasks() const { return fEventFilter->printOptionalMasks(); }
  TString printTrackOptionalMasks() const { return fTrackFilter->printOptionalMasks(); }
  TString printPIDOptionalMasks() const { return fPIDFilter->printOptionalMasks(); }
  /// \brief get the valid (armed) mandatory part mask associated to the configuration of the different objects
  /// \return the corresponding armed forced mask
  /// A clear example of the mandatory part mask is the mask of the zvertex and their
  /// alternatives where only a concrete one is required to be flagged
  uint64_t getCollisionForcedMask() { return fEventFilter->getForcedMask(); }
  uint64_t getTrackForcedMask() { return fTrackFilter->getForcedMask(); }
  uint64_t getPIDForcedMask() { return fPIDFilter->getForcedMask(); }
  /// \brief filter the different objects
  template <typename CollisionToFilter, typename AssociatedTracks>
  uint64_t FilterCollision(CollisionToFilter const& col, AssociatedTracks const& trks, int bfield)
  {
    return fEventFilter->Filter(col, trks, bfield);
  }
  template <typename TrackToFilter>
  uint64_t FilterTrack(TrackToFilter const& track)
  {
    return fTrackFilter->Filter(track);
  }
  template <typename TrackToFilter>
  uint64_t FilterTrackPID(TrackToFilter const& track)
  {
    return fPIDFilter->Filter(track);
  }
  /// \brief get the event multiplicities
  std::vector<float> GetCollisionMultiplicities() { return fEventFilter->GetMultiplicities(); }
  /// \brief Gets the index of the active multiplicity value within the multiplicities array
  int getCollisionMultiplicityIndex() { return fEventFilter->getMultiplicityIndex(); }

  /// \brief get the cut string signatures
  const TString& getEventFilterCutStringSignature() { return fEventFilter->getCutStringSignature(); }
  const TString& getTrackFilterCutStringSignature() { return fTrackFilter->getCutStringSignature(); }
  const TString& getPIDFilterCutStringSignature() { return fPIDFilter->getCutStringSignature(); }

 private:
  std::string ccdburl = "";
  std::string ccdbpath = "";
  std::string filterdate = "";
  bool forceupdate = false;
  PWGCF::TrackSelectionFilterAndAnalysis* fTrackFilter = nullptr;  /// the track filter
  PWGCF::EventSelectionFilterAndAnalysis* fEventFilter = nullptr;  /// the event filter
  PWGCF::PIDSelectionFilterAndAnalysis* fPIDFilter = nullptr;      /// the PID filter
  PWGCF::TrackSelectionFilterAndAnalysis* _fTrackFilter = nullptr; /// the track filter temporal storage until initialized
  PWGCF::EventSelectionFilterAndAnalysis* _fEventFilter = nullptr; /// the event filter temporal storage until initialized
  PWGCF::PIDSelectionFilterAndAnalysis* _fPIDFilter = nullptr;     /// the PID filter temporal storage until initialized

  ClassDefNV(FilterAndAnalysisFramework, 1)
};

extern void registerConfiguration(o2::analysis::PWGCF::FilterAndAnalysisFramework*);

} // namespace PWGCF
} // namespace analysis
} // namespace o2

#endif // PWGCF_TWOPARTICLECORRELATIONS_CORE_FILTERANDANALYSISFRAMEWORK_H_
