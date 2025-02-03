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

/// \file mchAlignRecord.cxx
/// \brief A useful task to get MillePede record for alignment and apply different geometries
///
/// \author Chi ZHANG, CEA-Saclay, chi.zhang@cern.ch

#include <gsl/span>
#include <cmath>
#include <string>
#include <unordered_map>
#include <vector>
#include <iostream>

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TParameter.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TChain.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TSystem.h>

#include "CommonConstants/LHCConstants.h"
#include "CommonUtils/NameConf.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGDQ/Core/VarManager.h"

#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/GRPGeomHelper.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/Logger.h"
#include "Framework/CallbackService.h"
#include "CCDB/BasicCCDBManager.h"

#include "MCHGeometryTransformer/Transformations.h"
#include "DataFormatsMCH/Cluster.h"
#include "DataFormatsMCH/TrackMCH.h"
#include "MCHTracking/Track.h"
#include "MCHTracking/TrackExtrap.h"
#include "MCHTracking/TrackParam.h"
#include "MCHTracking/TrackFitter.h"
#include "MCHBase/TrackerParam.h"
#include "ReconstructionDataFormats/TrackMCHMID.h"
#include "MCHAlign/Aligner.h"
#include "DetectorsCommonDataFormats/AlignParam.h"
#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsCommonDataFormats/DetectorNameConf.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

using namespace std;
using std::cout;
using std::endl;

const int fgNCh = 10;
const int fgNDetElemCh[fgNCh] = {4, 4, 4, 4, 18, 18, 26, 26, 26, 26};
const int fgSNDetElemCh[fgNCh + 1] = {0, 4, 8, 12, 16, 34, 52, 78, 104, 130, 156};

// using MyEvents = soa::Join<aod::Collisions, aod::EvSels>;
using MyEvents = aod::Collisions;

struct mchAlignRecordTask {

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  parameters::GRPMagField* grpmag;
  TGeoManager* geo;

  mch::TrackFitter trackFitter;
  mch::Aligner mAlign{};
  double mImproveCutChi2{};
  Double_t weightRecord{1.0};

  int fCurrentRun;

  map<int, math_utils::Transform3D> transformOld;
  map<int, math_utils::Transform3D> transformNew;
  mch::geo::TransformationCreator transformation;

  Configurable<string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<string> fFixChamber{"fix-chamber", "", "Fixing chamber"};
  Configurable<bool> fDoNewGeo{"do-realign", false, "Transform to a given new geometry"};
  Configurable<bool> fDoEvaluation{"do-evaluation", false, "Enable storage of residuals"};
  Configurable<string> fConfigNewGeoFile{"new-geo", "o2sim_geometry-aligned.root", "New geometry for transformation"};
  Configurable<double> fAllowedVarX{"variation-x", 2.0, "Allowed variation for x axis in cm"};
  Configurable<double> fAllowedVarY{"variation-y", 0.3, "Allowed variation for y axis in cm"};
  Configurable<double> fAllowedVarPhi{"variation-phi", 0.002, "Allowed variation for phi axis in rad"};
  Configurable<double> fAllowedVarZ{"variation-z", 2.0, "Allowed variation for z axis in cm"};
  Configurable<double> cfgSigmaX{"cfgSigmaX", 1000., "Sigma cut along X"};
  Configurable<double> cfgSigmaY{"cfgSigmaY", 1000., "Sigma cut along Y"};
  Configurable<double> cfgChamberResolutionX{"cfgChamberResolutionX", 0.04, "Chamber resolution along X configuration for refit"}; // 0.4cm pp, 0.2cm PbPb
  Configurable<double> cfgChamberResolutionY{"cfgChamberResolutionY", 0.04, "Chamber resolution along Y configuration for refit"}; // 0.4cm pp, 0.2cm PbPb
  Configurable<double> cfgSigmaCutImprove{"cfgSigmaCutImprove", 6., "Sigma cut for track improvement"};
  struct : ConfigurableGroup {
    Configurable<std::vector<int>> cfgDetElem{"cfgDetElem",
                                              {},
                                              "List of DetElem to be fixed"};
    Configurable<std::vector<int>> cfgParMask{"cfgParMask",
                                              {},
                                              "List of param mask for d.o.f to be fixed"};
  } fFixDetElem;

  void init(InitContext& ic)
  {

    // Load field and geometry informations here
    fCCDB->setURL(fConfigCcdbUrl);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();

    // Configuration for alignment object
    mAlign.SetDoEvaluation(fDoEvaluation.value);
    mAlign.SetAllowedVariation(0, fAllowedVarX.value);
    mAlign.SetAllowedVariation(1, fAllowedVarY.value);
    mAlign.SetAllowedVariation(2, fAllowedVarPhi.value);
    mAlign.SetAllowedVariation(3, fAllowedVarZ.value);
    mAlign.SetSigmaXY(cfgSigmaX.value, cfgSigmaY.value);

    // Configuration for track fitter
    const auto& trackerParam = o2::mch::TrackerParam::Instance();
    trackFitter.setBendingVertexDispersion(trackerParam.bendingVertexDispersion);
    trackFitter.setChamberResolution(cfgChamberResolutionX.value, cfgChamberResolutionY.value);
    trackFitter.smoothTracks(true);
    trackFitter.useChamberResolution();
    mImproveCutChi2 = 2. * cfgSigmaCutImprove.value * cfgSigmaCutImprove.value;

    // Configuration for chamber fixing
    TString chambersString = fFixChamber.value;
    std::unique_ptr<TObjArray> objArray(chambersString.Tokenize(","));
    if (objArray->GetEntries() > 0) {
      for (int iVar = 0; iVar < objArray->GetEntries(); ++iVar) {
        LOG(info) << Form("%s%d", "Fixing chamber: ", std::stoi(objArray->At(iVar)->GetName()));
        mAlign.FixChamber(std::stoi(objArray->At(iVar)->GetName()));
      }
    }

    // Configuration for DE fixing with given axes
    auto DEs = fFixDetElem.cfgDetElem.value;
    auto Masks = fFixDetElem.cfgParMask.value;
    if (DEs.size() > 0) {
      if (DEs.size() != Masks.size()) {
        LOG(fatal) << "Inconsistent size of mask list.";
      }
      for (int i = 0; i < static_cast<int>(DEs.size()); i++) {
        LOG(info) << Form("%s%d%s%d", "Fixing DE: ", DEs.at(i), " with mask: ", Masks.at(i));
        mAlign.FixDetElem(DEs.at(i), Masks.at(i));
      }
    }

    // Init for output saving
    mAlign.init();

    ic.services().get<CallbackService>().set<CallbackService::Id::Stop>([this]() {
      LOG(info) << "Saving records into ROOT file";
      mAlign.terminate();
    });
  }

  //_________________________________________________________________________________________________
  Int_t GetDetElemId(Int_t iDetElemNumber)
  {
    // make sure detector number is valid
    if (!(iDetElemNumber >= fgSNDetElemCh[0] &&
          iDetElemNumber < fgSNDetElemCh[fgNCh])) {
      LOG(fatal) << "Invalid detector element number: " << iDetElemNumber;
    }
    /// get det element number from ID
    // get chamber and element number in chamber
    int iCh = 0;
    int iDet = 0;
    for (int i = 1; i <= fgNCh; i++) {
      if (iDetElemNumber < fgSNDetElemCh[i]) {
        iCh = i;
        iDet = iDetElemNumber - fgSNDetElemCh[i - 1];
        break;
      }
    }

    // make sure detector index is valid
    if (!(iCh > 0 && iCh <= fgNCh && iDet < fgNDetElemCh[iCh - 1])) {
      LOG(fatal) << "Invalid detector element id: " << 100 * iCh + iDet;
    }

    // add number of detectors up to this chamber
    return 100 * iCh + iDet;
  }

  //_________________________________________________________________________________________________
  bool RemoveTrack(mch::Track& track)
  {
    bool removeTrack = false;

    try {
      trackFitter.fit(track, false);
    } catch (exception const& e) {
      removeTrack = true;
      return removeTrack;
    }

    auto itStartingParam = std::prev(track.rend());

    while (true) {
      try {
        trackFitter.fit(track, true, false, (itStartingParam == track.rbegin()) ? nullptr : &itStartingParam);
      } catch (exception const&) {
        removeTrack = true;
        break;
      }

      double worstLocalChi2 = -1.0;

      track.tagRemovableClusters(0x1F, false);
      auto itWorstParam = track.end();

      for (auto itParam = track.begin(); itParam != track.end(); ++itParam) {
        if (itParam->getLocalChi2() > worstLocalChi2) {
          worstLocalChi2 = itParam->getLocalChi2();
          itWorstParam = itParam;
        }
      }

      if (worstLocalChi2 < mImproveCutChi2)
        break;

      if (!itWorstParam->isRemovable()) {
        removeTrack = true;
        track.removable();
        break;
      }

      auto itNextParam = track.removeParamAtCluster(itWorstParam);
      auto itNextToNextParam = (itNextParam == track.end()) ? itNextParam : std::next(itNextParam);
      itStartingParam = track.rbegin();

      if (track.getNClusters() < 10) {
        removeTrack = true;
        break;
      } else {
        while (itNextToNextParam != track.end()) {
          if (itNextToNextParam->getClusterPtr()->getChamberId() != itNextParam->getClusterPtr()->getChamberId()) {
            itStartingParam = std::make_reverse_iterator(++itNextParam);
            break;
          }
          ++itNextToNextParam;
        }
      }
    }

    if (!removeTrack) {
      for (auto& param : track) {
        param.setParameters(param.getSmoothParameters());
        param.setCovariances(param.getSmoothCovariances());
      }
    }

    return removeTrack;
  }

  template <typename TEvent, typename TTracks, typename TClusters>
  void runProcessTracks(TEvent const& collision, aod::BCsWithTimestamps const&, TTracks const& tracks, TClusters const& clusters)
  {

    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();

    if (fCurrentRun != bc.runNumber()) {

      grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
      if (grpmag != nullptr) {
        o2::base::Propagator::initFieldFromGRP(grpmag);
      }
      // Configuration for magnetic field and track extrapolation
      VarManager::SetupMuonMagField();
      mch::TrackExtrap::useExtrapV2();
      mAlign.SetBFieldOn(mch::TrackExtrap::isFieldON());
      trackFitter.smoothTracks(true);

      // Load reference geometry
      geo = fCCDB->getForTimeStamp<TGeoManager>(geoPath, bc.timestamp());
      transformation = mch::geo::transformationFromTGeoManager(*geo);
      for (int i = 0; i < 156; i++) {
        int iDEN = GetDetElemId(i);
        transformOld[iDEN] = transformation(iDEN);
      }

      if (fDoNewGeo.value) {
        // Load new geometry with which we want to check
        base::GeometryManager::loadGeometry(fConfigNewGeoFile.value);
        transformation = mch::geo::transformationFromTGeoManager(*gGeoManager);
        for (int i = 0; i < 156; i++) {
          int iDEN = GetDetElemId(i);
          transformNew[iDEN] = transformation(iDEN);
        }
      }

      fCurrentRun = bc.runNumber();
    }

    // Loop over forward tracks
    for (auto const& track : tracks) {

      int clIndex = -1;
      mch::Track convertedTrack;

      // Use only MCH-MID matched tracks
      if (static_cast<int>(track.trackType()) != 3) {
        continue;
      }

      // Loop over attached clusters
      for (auto const& cluster : clusters) {

        if (cluster.template fwdtrack_as<TTracks>() != track) {
          continue;
        }

        clIndex += 1;

        mch::Cluster* mch_cluster = new mch::Cluster();
        mch_cluster->x = cluster.x();
        mch_cluster->y = cluster.y();
        mch_cluster->z = cluster.z();

        if (fDoNewGeo.value) {
          math_utils::Point3D<double> local;
          math_utils::Point3D<double> master;

          master.SetXYZ(cluster.x(), cluster.y(), cluster.z());

          transformOld[cluster.deId()].MasterToLocal(master, local);
          transformNew[cluster.deId()].LocalToMaster(local, master);

          mch_cluster->x = master.x();
          mch_cluster->y = master.y();
          mch_cluster->z = master.z();
        }

        uint32_t ClUId = mch::Cluster::buildUniqueId(static_cast<int>(cluster.deId() / 100) - 1, cluster.deId(), clIndex);
        mch_cluster->uid = ClUId;

        mch_cluster->ex = cluster.isGoodX() ? 0.2 : 10.0;
        mch_cluster->ey = cluster.isGoodY() ? 0.2 : 10.0;

        convertedTrack.createParamAtCluster(*mch_cluster);
      }

      if (convertedTrack.getNClusters() > 9) {
        // Erase removable track
        if (!RemoveTrack(convertedTrack)) {
          mAlign.ProcessTrack(convertedTrack, transformation, true, weightRecord);
        }
      }
    }
  }

  void processTracks(MyEvents::iterator const& collision, aod::BCsWithTimestamps const& bcs, aod::FwdTracks const& tracks, aod::FwdTrkCls const& clusters)
  {
    runProcessTracks(collision, bcs, tracks, clusters);
  }

  PROCESS_SWITCH(mchAlignRecordTask, processTracks, "Process tracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<mchAlignRecordTask>(cfgc)};
}
