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

/// \file utilsTpcSkimsTableCreator.h
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>
/// \author Christian Sonnabend <christian.sonnabend@cern.ch>
/// \author Jeremy Wilkinson <jeremy.wilkinson@cern.ch>
/// \author Ana Marin <ana.marin@cern.ch>
/// \author Oleksii Lubynets <oleksii.lubynets@cern.ch>
/// \brief  Helper functions used both in V0 and TOF structs of tpcSkimsTableCreator.cxx

#ifndef DPG_TASKS_TPC_UTILSTPCSKIMSTABLECREATOR_H_
#define DPG_TASKS_TPC_UTILSTPCSKIMSTABLECREATOR_H_

#include "tpcSkimsTableCreator.h"

#include "Common/Core/CollisionTypeHelper.h"
#include "Common/DataModel/OccupancyTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/Logger.h>
#include <ReconstructionDataFormats/PID.h>

#include <TFile.h>
#include <TH2.h>
#include <TRandom3.h>

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <map>
#include <string>

namespace o2::dpg_tpcskimstablecreator
{
enum {
  TrackSelectionNoCut = 0,
  TrackSelectionGlobalTrack,
  TrackSelectionTrackWoPtEta,
  TrackSelectionGlobalTrackWoDCA,
  TrackSelectionQualityTracks,
  TrackSelectionInAcceptanceTracks
};

enum {
  EventSelectionNo = 0,
  EventSelectionRun2,
  EventSelectionRun3
};

/// Event selection
template <typename CollisionType>
bool isEventSelected(const CollisionType& collision, const int applyEvSel)
{
  if ((applyEvSel == EventSelectionRun2 && !collision.sel7()) || (applyEvSel == EventSelectionRun3 && !collision.sel8())) {
    return false;
  } else {
    return true;
  }
};

// Track selection
template <typename TrackType>
bool isTrackSelected(const TrackType& track, const int trackSelection)
{
  bool isSelected{false};
  isSelected |= trackSelection == TrackSelectionNoCut;
  isSelected |= (trackSelection == TrackSelectionGlobalTrack) && track.isGlobalTrack();
  isSelected |= (trackSelection == TrackSelectionTrackWoPtEta) && track.isGlobalTrackWoPtEta();
  isSelected |= (trackSelection == TrackSelectionGlobalTrackWoDCA) && track.isGlobalTrackWoDCA();
  isSelected |= (trackSelection == TrackSelectionQualityTracks) && track.isQualityTrack();
  isSelected |= (trackSelection == TrackSelectionInAcceptanceTracks) && track.isInAcceptanceTrack();

  return isSelected;
}

/// Evaluate tpcSignal with or without dEdx correction
template <bool IsCorrectedDeDx, typename TrkType>
double tpcSignalGeneric(const TrkType& track)
{
  if constexpr (IsCorrectedDeDx) {
    return track.tpcSignalCorrected();
  } else {
    return track.tpcSignal();
  }
}

/// Determine interaction rate source from CCDB
std::string evaluateIrSource(const o2::framework::Service<o2::ccdb::BasicCCDBManager>& ccdb, const std::string& ccdbPathGrpLhcIf, const uint64_t timestamp)
{
  std::string irSource{""};
  o2::parameters::GRPLHCIFData* genRunParams = ccdb->template getForTimeStamp<o2::parameters::GRPLHCIFData>(ccdbPathGrpLhcIf, timestamp);
  if (genRunParams != nullptr) {
    const auto collSys = CollisionSystemType::getCollisionTypeFromGrp(genRunParams);
    if (collSys == CollisionSystemType::kCollSyspp) {
      irSource = "T0VTX";
    } else {
      irSource = "ZNC hadronic";
    }
    LOG(info) << "irSource determined from General Run Parameters: " << irSource;
  } else {
    LOG(warning) << "No General Run Parameters object found. irSource will remain undefined";
  }
  return irSource;
}

struct OccupancyValues {
  float tmoPrimUnfm80{UndefValueFloat};
  float tmoFV0AUnfm80{UndefValueFloat};
  float tmoFT0AUnfm80{UndefValueFloat};
  float tmoFT0CUnfm80{UndefValueFloat};
  float tmoRT0V0PrimUnfm80{UndefValueFloat};
  float twmoPrimUnfm80{UndefValueFloat};
  float twmoFV0AUnfm80{UndefValueFloat};
  float twmoFT0AUnfm80{UndefValueFloat};
  float twmoFT0CUnfm80{UndefValueFloat};
  float twmoRT0V0PrimUnfm80{UndefValueFloat};
};

using TrackMeanOccs = soa::Join<aod::TmoTrackIds, aod::TmoToTrackQA, aod::TmoPrim, aod::TmoT0V0, aod::TmoRT0V0Prim, aod::TwmoPrim, aod::TwmoT0V0, aod::TwmoRT0V0Prim>;

/// Evaluate occupancy-related variables
template <typename TrkType>
void evaluateOccupancyVariables(const TrkType& track, OccupancyValues& occValues)
{
  if (track.tmoId() == -1) {
    return;
  }
  const auto tmoFromTrack = track.template tmo_as<TrackMeanOccs>();
  occValues.tmoPrimUnfm80 = tmoFromTrack.tmoPrimUnfm80();
  occValues.tmoFV0AUnfm80 = tmoFromTrack.tmoFV0AUnfm80();
  occValues.tmoFT0AUnfm80 = tmoFromTrack.tmoFT0AUnfm80();
  occValues.tmoFT0CUnfm80 = tmoFromTrack.tmoFT0CUnfm80();
  occValues.tmoRT0V0PrimUnfm80 = tmoFromTrack.tmoRobustT0V0PrimUnfm80();
  occValues.twmoPrimUnfm80 = tmoFromTrack.twmoPrimUnfm80();
  occValues.twmoFV0AUnfm80 = tmoFromTrack.twmoFV0AUnfm80();
  occValues.twmoFT0AUnfm80 = tmoFromTrack.twmoFT0AUnfm80();
  occValues.twmoFT0CUnfm80 = tmoFromTrack.twmoFT0CUnfm80();
  occValues.twmoRT0V0PrimUnfm80 = tmoFromTrack.twmoRobustT0V0PrimUnfm80();
}

inline void initAdaptDwnsplHistograms(const std::string& ccdbPathForAdaptDwnsplFile, const std::string& adaptDwnsplFileName, const std::string& adaptDwnsplHistoName, TH2*& histo)
{
  o2::ccdb::CcdbApi ccdbApi;
  ccdbApi.init("http://alice-ccdb.cern.ch");
  const std::map<std::string, std::string> metadata;
  if (!std::filesystem::exists(adaptDwnsplFileName.c_str())) {
    // Retrieve the file from CCDB
    const bool isFileAvailable = ccdbApi.retrieveBlob(ccdbPathForAdaptDwnsplFile, ".", metadata, -1, false, adaptDwnsplFileName);
    if (!isFileAvailable) {
      LOGF(fatal, "Failed to retrieve expected-momentum histograms file from CCDB: %s/%s", ccdbPathForAdaptDwnsplFile.c_str(), adaptDwnsplFileName.c_str());
      return;
    }
  } else {
    LOGF(fatal, "The file %s already exists. Make sure that 'adaptDwnsplFileName' has different values for V0 and TPC-TOF", adaptDwnsplFileName.c_str());
    return;
  }
  TFile* adaptDwsplFile = TFile::Open(adaptDwnsplFileName.c_str(), "read");
  if (adaptDwsplFile == nullptr || adaptDwsplFile->IsZombie()) {
    LOGF(fatal, "Failed to open expected-momentum histograms file from CCDB: %s/%s", ccdbPathForAdaptDwnsplFile.c_str(), adaptDwnsplFileName.c_str());
    return;
  }
  histo = dynamic_cast<TH2*>(adaptDwsplFile->Get<TH2>(adaptDwnsplHistoName.c_str())->Clone());
  if (histo == nullptr) {
    LOGF(fatal, "Histogram %s not found in expected-momentum histograms file!", adaptDwnsplHistoName.c_str());
    return;
  }
  if (histo->GetXaxis()->IsVariableBinSize()) {
    LOGF(fatal, "Histogram %s has variable bin size, which is not allowed!", adaptDwnsplHistoName.c_str());
    return;
  }
  histo->SetDirectory(nullptr);
  adaptDwsplFile->Close();
  delete adaptDwsplFile;
  adaptDwsplFile = nullptr;
}

inline bool performAdaptiveDownsampling(TRandom3* fRndm, const double p, const o2::track::PID::ID pid, const double targetYieldPerBin, const TH2* histoExpected, const bool isPositive)
{
  const int pidShift = isPositive ? 0 : o2::track::PID::NIDs;
  const double yieldExpected = histoExpected->GetBinContent(histoExpected->GetXaxis()->FindBin(std::log10(p)), histoExpected->GetYaxis()->FindBin(pidShift + pid));
  if (yieldExpected < targetYieldPerBin) {
    return true;
  }
  if (fRndm->Uniform(1) > (targetYieldPerBin / yieldExpected)) {
    return false;
  }
  return true;
}

} // namespace o2::dpg_tpcskimstablecreator
#endif // DPG_TASKS_TPC_UTILSTPCSKIMSTABLECREATOR_H_
