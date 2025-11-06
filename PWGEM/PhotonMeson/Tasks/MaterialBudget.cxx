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

/// \file MaterialBudget.cxx
/// \brief Task to analyse and calculate the material budget weights
/// \author S. Mrozinski, smrozins@cern.ch

#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Core/CutsLibrary.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/Core/HistogramsLibrary.h"
#include "PWGEM/PhotonMeson/Core/PairCut.h"
#include "PWGEM/PhotonMeson/Core/V0PhotonCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/MCUtilities.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TPCVDriftManager.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/McCollisionExtra.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/runDataProcessing.h>

#include <TMCProcess.h>
#include <TPDGCode.h>
#include <TProfile2D.h>

#include <algorithm>
#include <deque>
#include <iterator>
#include <map>
#include <ranges>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photonmeson::photonpair;
using namespace o2::aod::pwgem::photon;
using namespace o2::aod::pwgem::dilepton::utils::emtrackutil;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;
using namespace o2::aod::pwgem::dilepton::utils;
using o2::constants::math::TwoPI;

using MyCollisions = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent>;
using MyCollision = MyCollisions::iterator;

using MyCollisionsMC = soa::Join<aod::EMEvents, aod::EMMCEventLabels>;
using MyCollisionMC = MyCollisionsMC::iterator;

using MyMCCollisions = soa::Join<aod::EMMCEvents, aod::BinnedGenPts, aod::MostProbableEMEventIdsInMC>;
using MyMCCollision = MyMCCollisions::iterator;

using MyV0Photons = soa::Join<aod::V0PhotonsKF, aod::V0KFEMEventIds>;
using MyV0Photon = MyV0Photons::iterator;

using MyPrimaryElectrons = soa::Filtered<soa::Join<aod::EMPrimaryElectronsFromDalitz, aod::EMPrimaryElectronEMEventIds, aod::EMPrimaryElectronsPrefilterBitDerived>>;
using MyPrimaryElectron = MyPrimaryElectrons::iterator;

using MyMCV0Legs = soa::Join<aod::V0Legs, aod::V0LegMCLabels>;
using MyMCV0Leg = MyMCV0Legs::iterator;

using MyTracks = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::McTrackLabels>;
struct MaterialBudget {

  Configurable<int> cfgCentEstimator{"cfgCentEstimator", 2, "FT0M:0, FT0A:1, FT0C:2"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0, "min. centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 999.f, "max. centrality"};

  struct SimplePhoton {
    float pt;
    float eta;
    float phi;
    float vz;
    float r;
  };

  static constexpr int MaxMixEvents = 5;
  std::map<std::tuple<int, int>, std::deque<std::vector<SimplePhoton>>> mixingPools;

  std::vector<float> centBinEdges = {0, 10, 30, 50, 90, 200};
  std::vector<float> zvtxBinEdges = {-10, -5, 0, 5, 10};

  std::tuple<int, int> getPoolBin(float cent, float zvtx)
  {
    int centbin = std::lower_bound(centBinEdges.begin(), centBinEdges.end(), cent) - centBinEdges.begin() - 1;
    centbin = std::max(0, std::min(centbin, static_cast<int>(centBinEdges.size()) - 2));

    int zbin = std::lower_bound(zvtxBinEdges.begin(), zvtxBinEdges.end(), zvtx) - zvtxBinEdges.begin() - 1;
    zbin = std::max(0, std::min(zbin, static_cast<int>(zvtxBinEdges.size()) - 2));
    return {centbin, zbin};
  }

  HistogramRegistry registry{"output", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  static constexpr std::string_view ItsClsNames[] = {
    "ITSCls0", "ITSCls1", "ITSCls2", "ITSCls3",
    "ITSCls4", "ITSCls5", "ITSCls6", "ITSCls7"};

  static constexpr std::string_view EventTypes[2] = {"before/", "after/"};

  Configurable<bool> cfgPlotBremsstrahlung{"cfgPlotBremsstrahlung", false, "produce plots to study Bremsstrahlung"};
  Configurable<bool> cfgPlotResolution{"cfgPlotResolution", false, "produce plots to study resolution"};
  Configurable<bool> cfgPlotPurity{"cfgPlotPurity", false, "produce plots to study purity"};
  Configurable<bool> cfgPlotMBDetailed{"cfgPlotMBDetailed", false, "produce plots to study material distribution distribution"};
  Configurable<bool> cfgPlotMBGeneral{"cfgPlotMBGeneral", true, "produce plots to study material distribution"};
  Configurable<bool> cfgPlotMBCollisions{"cfgPlotMBCollisions", false, "produce plots to study material distribution if collision association is wrong"};

  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0}, ""};
  const AxisSpec axisPt{binsPt, "#it{p}_{T} [GeV/c]"};

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", true, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -2, "min. occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -2, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStandard{"cfgRequireNoCollInTimeRangeStandard", false, "require no collision in time range standard"};
    Configurable<bool> cfgRequireNoCollInTimeRangeStrict{"cfgRequireNoCollInTimeRangeStrict", false, "require no collision in time range strict"};
    Configurable<bool> cfgRequireNoCollInITSROFStandard{"cfgRequireNoCollInITSROFStandard", false, "require no collision in time range standard"};
    Configurable<bool> cfgRequireNoCollInITSROFStrict{"cfgRequireNoCollInITSROFStrict", false, "require no collision in time range strict"};
    Configurable<bool> cfgRequireNoHighMultCollInPrevRof{"cfgRequireNoHighMultCollInPrevRof", false, "require no HM collision in previous ITS ROF"};
  } eventcuts;

  DalitzEECut fDileptonCut;
  struct : ConfigurableGroup {
    std::string prefix = "dileptoncut_group";
    Configurable<float> cfgMinMass{"cfgMinMass", 0.0, "min mass"};
    Configurable<float> cfgMaxMass{"cfgMaxMass", 0.1, "max mass"};
    Configurable<bool> cfgApplyPhiv{"cfgApplyPhiv", true, "flag to apply phiv cut"};
    Configurable<bool> cfgRequireItsibAny{"cfgRequireItsibAny", false, "flag to require ITS ib any hits"};
    Configurable<bool> cfgRequireItsib1st{"cfgRequireItsib1st", true, "flag to require ITS ib 1st hit"};
    Configurable<float> cfgPhivSlope{"cfgPhivSlope", 0.0185, "slope for m vs. phiv"};
    Configurable<float> cfgPhivIntercept{"cfgPhivIntercept", -0.0280, "intercept for m vs. phiv"};

    Configurable<float> cfgMinPtTrack{"cfgMinPtTrack", 0.1, "min pT for single track"};
    Configurable<float> cfgMaxEtaTrack{"cfgMaxEtaTrack", 0.8, "max eta for single track"};
    Configurable<int> cfgMinnclusterTPC{"cfgMinnclusterTPC", 0, "min ncluster tpc"};
    Configurable<int> cfgMinnclusterITS{"cfgMinnclusterITS", 5, "min ncluster its"};
    Configurable<int> cfgMinncrossedrows{"cfgMinncrossedrows", 70, "min ncrossed rows"};
    Configurable<float> cfgMaxchi2TPC{"cfgMaxchi2TPC", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfgMaxchi2ITS{"cfgMaxchi2ITS", 36.0, "max chi2/NclsITS"};
    Configurable<float> cfgMaxDCAxy{"cfgMaxDCAxy", 0.05, "max dca XY for single track in cm"};
    Configurable<float> cfgMaxDCAz{"cfgMaxDCAz", 0.05, "max dca Z for single track in cm"};
    Configurable<float> cfgMaxDCA3dsigmaTrack{"cfgMaxDCA3dsigmaTrack", 1.5, "max DCA 3D in sigma"};
    Configurable<float> cfgMaxFracSharedClustersTPC{"cfgMaxFracSharedClustersTPC", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<bool> cfgApplyCutsFromPrefilterDerived{"cfgApplyCutsFromPrefilterDerived", false, "flag to apply prefilter to electron"};
    Configurable<bool> includeITSsa{"includeITSsa", false, "Flag to enable ITSsa tracks"};
    Configurable<float> cfgMaxpttrackITSsa{"cfgMaxpttrackITSsa", 0.15, "max pt for ITSsa tracks"};

    Configurable<int> cfgPIDscheme{"cfgPIDscheme", static_cast<int>(DalitzEECut::PIDSchemes::kTOFif), "pid scheme [kTOFif : 0, kTPConly : 1]"};
    Configurable<float> cfgMinTPCNsigmaEl{"cfgMinTPCNsigmaEl", -2.0, "min. TPC n sigma for electron inclusion"};
    Configurable<float> cfgMaxTPCNsigmaEl{"cfgMaxTPCNsigmaEl", +3.0, "max. TPC n sigma for electron inclusion"};
    Configurable<float> cfgMinTPCNsigmaPi{"cfgMinTPCNsigmaPi", -0.0, "min. TPC n sigma for pion exclusion"};
    Configurable<float> cfgMaxTPCNsigmaPi{"cfgMaxTPCNsigmaPi", +0.0, "max. TPC n sigma for pion exclusion"};
    Configurable<float> cfgMinTOFNsigmaEl{"cfgMinTOFNsigmaEl", -3.0, "min. TOF n sigma for electron inclusion"};
    Configurable<float> cfgMaxTOFNsigmaEl{"cfgMaxTOFNsigmaEl", +3.0, "max. TOF n sigma for electron inclusion"};
  } dileptoncuts;

  V0PhotonCut fV0PhotonCut;
  struct : ConfigurableGroup {
    std::string prefix = "pcmcut_group";
    Configurable<bool> cfgRequireV0WithITSTPC{"cfgRequireV0WithITSTPC", false, "flag to select V0s with ITS-TPC matched tracks"};
    Configurable<bool> cfgRequireV0WithITSOnly{"cfgRequireV0WithITSOnly", false, "flag to select V0s with ITSonly tracks"};
    Configurable<bool> cfgRequireV0WithTPCOnly{"cfgRequireV0WithTPCOnly", false, "flag to select V0s with TPConly tracks"};
    Configurable<float> cfgMinPtV0{"cfgMinPtV0", 0.1, "min pT for v0 photons at PV"};
    Configurable<float> cfgMaxPtV0{"cfgMaxPtV0", 1e+10, "max pT for v0 photons at PV"};
    Configurable<float> cfgMinEtaV0{"cfgMinEtaV0", -0.8, "min eta for v0 photons at PV"};
    Configurable<float> cfgMaxEtaV0{"cfgMaxEtaV0", +0.8, "max eta for v0 photons at PV"};
    Configurable<float> cfgMinV0Radius{"cfgMinV0Radius", 4.0, "min v0 radius"};
    Configurable<float> cfgMaxV0Radius{"cfgMaxV0Radius", 90.0, "max v0 radius"};
    Configurable<float> cfgMaxAlphaAp{"cfgMaxAlphaAp", 0.95, "max alpha for AP cut"};
    Configurable<float> cfgMaxQtAp{"cfgMaxQtAp", 0.01, "max qT for AP cut"};
    Configurable<float> cfgMinCospa{"cfgMinCospa", 0.999, "min V0 CosPA"};
    Configurable<float> cfgMaxPca{"cfgMaxPca", 1.5, "max distance btween 2 legs"};
    Configurable<float> cfgMaxChi2kf{"cfgMaxChi2kf", 1e+10, "max chi2/ndf with KF"};
    Configurable<bool> cfgRejectV0OnITSib{"cfgRejectV0OnITSib", true, "flag to reject V0s on ITSib"};
    Configurable<int> cfgMinNclusterTPC{"cfgMinNclusterTPC", 0, "min ncluster tpc"};
    Configurable<int> cfgMinNcrossedrows{"cfgMinNcrossedrows", 40, "min ncrossed rows"};
    Configurable<float> cfgMaxFracSharedClustersTPC{"cfgMaxFracSharedClustersTPC", 999.f, "max fraction of shared clusters in TPC"};
    Configurable<float> cfgMaxChi2TPC{"cfgMaxChi2TPC", 4.0, "max chi2/NclsTPC"};
    Configurable<float> cfgMaxChi2ITS{"cfgMaxChi2ITS", 36.0, "max chi2/NclsITS"};
    Configurable<float> cfgMinTPCNsigmaEl{"cfgMinTPCNsigmaEl", -3.0, "min. TPC n sigma for electron"};
    Configurable<float> cfgMaxTPCNsigmaEl{"cfgMaxTPCNsigmaEl", +3.0, "max. TPC n sigma for electron"};
    Configurable<bool> cfgDisableITSonlytrack{"cfgDisableITSonlytrack", false, "flag to disable ITSonly tracks"};
    Configurable<bool> cfgDisableTPCOnlytrack{"cfgDisableTPCOnlytrack", false, "flag to disable TPConly tracks"};
  } pcmcuts;

  void init(InitContext&)
  {

    defineEMEventCut();
    defineDileptonCut();
    definePCMCut();
    addhistograms();
  }

  void defineEMEventCut()
  {
    fEMEventCut = EMPhotonEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(-eventcuts.cfgZvtxMax, +eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
    fEMEventCut.SetRequireNoCollInTimeRangeStandard(eventcuts.cfgRequireNoCollInTimeRangeStandard);
    fEMEventCut.SetRequireNoCollInTimeRangeStrict(eventcuts.cfgRequireNoCollInTimeRangeStrict);
    fEMEventCut.SetRequireNoCollInITSROFStandard(eventcuts.cfgRequireNoCollInITSROFStandard);
    fEMEventCut.SetRequireNoCollInITSROFStrict(eventcuts.cfgRequireNoCollInITSROFStrict);
    fEMEventCut.SetRequireNoHighMultCollInPrevRof(eventcuts.cfgRequireNoHighMultCollInPrevRof);
  }

  void defineDileptonCut()
  {
    fDileptonCut = DalitzEECut("fDileptonCut", "fDileptonCut");

    // for pair
    fDileptonCut.SetMeeRange(dileptoncuts.cfgMinMass, dileptoncuts.cfgMaxMass);
    fDileptonCut.SetMaxPhivPairMeeDep([&](float mll) { return (mll - dileptoncuts.cfgPhivIntercept) / dileptoncuts.cfgPhivSlope; });
    fDileptonCut.ApplyPhiV(dileptoncuts.cfgApplyPhiv);
    fDileptonCut.RequireITSibAny(dileptoncuts.cfgRequireItsibAny);
    fDileptonCut.RequireITSib1st(dileptoncuts.cfgRequireItsib1st);

    // for tracks
    fDileptonCut.SetTrackPtRange(dileptoncuts.cfgMinPtTrack, 1e+10f);
    fDileptonCut.SetTrackEtaRange(-dileptoncuts.cfgMaxEtaTrack, +dileptoncuts.cfgMaxEtaTrack);
    fDileptonCut.SetMinNClustersTPC(dileptoncuts.cfgMinnclusterTPC);
    fDileptonCut.SetMinNCrossedRowsTPC(dileptoncuts.cfgMinncrossedrows);
    fDileptonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fDileptonCut.SetMaxFracSharedClustersTPC(dileptoncuts.cfgMaxFracSharedClustersTPC);
    fDileptonCut.SetChi2PerClusterTPC(0.0, dileptoncuts.cfgMaxchi2TPC);
    fDileptonCut.SetChi2PerClusterITS(0.0, dileptoncuts.cfgMaxchi2ITS);
    fDileptonCut.SetNClustersITS(dileptoncuts.cfgMinnclusterITS, 7);
    fDileptonCut.SetMaxDcaXY(dileptoncuts.cfgMaxDCAxy);
    fDileptonCut.SetMaxDcaZ(dileptoncuts.cfgMaxDCAz);
    fDileptonCut.SetTrackDca3DRange(0.f, dileptoncuts.cfgMaxDCA3dsigmaTrack); // in sigma
    fDileptonCut.IncludeITSsa(dileptoncuts.includeITSsa, dileptoncuts.cfgMaxpttrackITSsa);

    // for eID
    fDileptonCut.SetPIDScheme(dileptoncuts.cfgPIDscheme);
    fDileptonCut.SetTPCNsigmaElRange(dileptoncuts.cfgMinTPCNsigmaEl, dileptoncuts.cfgMaxTPCNsigmaEl);
    fDileptonCut.SetTPCNsigmaPiRange(dileptoncuts.cfgMinTPCNsigmaPi, dileptoncuts.cfgMaxTPCNsigmaPi);
    fDileptonCut.SetTOFNsigmaElRange(dileptoncuts.cfgMinTOFNsigmaEl, dileptoncuts.cfgMaxTOFNsigmaEl);
  }

  void definePCMCut()
  {
    fV0PhotonCut = V0PhotonCut("fV0PhotonCut", "fV0PhotonCut");

    // for v0
    fV0PhotonCut.SetV0PtRange(pcmcuts.cfgMinPtV0, pcmcuts.cfgMaxPtV0);
    fV0PhotonCut.SetV0EtaRange(pcmcuts.cfgMinEtaV0, pcmcuts.cfgMaxEtaV0);
    fV0PhotonCut.SetMinCosPA(pcmcuts.cfgMinCospa);
    fV0PhotonCut.SetMaxPCA(pcmcuts.cfgMaxPca);
    fV0PhotonCut.SetMaxChi2KF(pcmcuts.cfgMaxChi2kf);
    fV0PhotonCut.SetRxyRange(pcmcuts.cfgMinV0Radius, pcmcuts.cfgMaxV0Radius);
    fV0PhotonCut.SetAPRange(pcmcuts.cfgMaxAlphaAp, pcmcuts.cfgMaxQtAp);
    fV0PhotonCut.RejectITSib(pcmcuts.cfgRejectV0OnITSib);

    // for track
    fV0PhotonCut.SetMinNClustersTPC(pcmcuts.cfgMinNclusterTPC);
    fV0PhotonCut.SetMinNCrossedRowsTPC(pcmcuts.cfgMinNcrossedrows);
    fV0PhotonCut.SetMinNCrossedRowsOverFindableClustersTPC(0.8);
    fV0PhotonCut.SetMaxFracSharedClustersTPC(pcmcuts.cfgMaxFracSharedClustersTPC);
    fV0PhotonCut.SetChi2PerClusterTPC(0.0, pcmcuts.cfgMaxChi2TPC);
    fV0PhotonCut.SetTPCNsigmaElRange(pcmcuts.cfgMinTPCNsigmaEl, pcmcuts.cfgMaxTPCNsigmaEl);
    fV0PhotonCut.SetChi2PerClusterITS(-1e+10, pcmcuts.cfgMaxChi2ITS);
    fV0PhotonCut.SetNClustersITS(0, 7);
    fV0PhotonCut.SetMeanClusterSizeITSob(0.0, 16.0);
    fV0PhotonCut.SetDisableITSonly(pcmcuts.cfgDisableITSonlytrack);
    fV0PhotonCut.SetDisableTPConly(pcmcuts.cfgDisableTPCOnlytrack);
    fV0PhotonCut.SetRequireITSTPC(pcmcuts.cfgRequireV0WithITSTPC);
    fV0PhotonCut.SetRequireITSonly(pcmcuts.cfgRequireV0WithITSOnly);
    fV0PhotonCut.SetRequireTPConly(pcmcuts.cfgRequireV0WithTPCOnly);
  }

  void addhistograms()
  {

    auto hCollisionCounter = registry.add<TH1>("Event/before/hCollisionCounter", "collision counter;;Number of events", kTH1F, {{10, 0.5, 10.5}}, false);
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "No TF border");
    hCollisionCounter->GetXaxis()->SetBinLabel(3, "No ITS ROF border");
    hCollisionCounter->GetXaxis()->SetBinLabel(4, "No Same Bunch Pileup");
    hCollisionCounter->GetXaxis()->SetBinLabel(5, "Is Vertex ITSTPC");
    hCollisionCounter->GetXaxis()->SetBinLabel(6, "Is Good Zvtx FT0vsPV");
    hCollisionCounter->GetXaxis()->SetBinLabel(7, "FT0AND");
    hCollisionCounter->GetXaxis()->SetBinLabel(8, "sel8");
    hCollisionCounter->GetXaxis()->SetBinLabel(9, "|Z_{vtx}| < 10 cm");
    hCollisionCounter->GetXaxis()->SetBinLabel(10, "accepted");

    registry.add("Event/before/hZvtx", "vertex z; Z_{vtx} (cm)", kTH1F, {{100, -50, +50}}, false);
    registry.add("Event/before/hMultNTracksPV", "hMultNTracksPV; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
    registry.add("Event/before/hMultNTracksPVeta1", "hMultNTracksPVeta1; N_{track} to PV", kTH1F, {{6001, -0.5, 6000.5}}, false);
    registry.add("Event/before/hMultFT0", "hMultFT0;mult. FT0A;mult. FT0C", kTH2F, {{300, 0, 6000}, {300, 0, 6000}}, false);
    registry.add("Event/before/hCentFT0A", "hCentFT0A;centrality FT0A (%)", kTH1F, {{110, 0, 110}}, false);
    registry.add("Event/before/hCentFT0C", "hCentFT0C;centrality FT0C (%)", kTH1F, {{110, 0, 110}}, false);
    registry.add("Event/before/hCentFT0M", "hCentFT0M;centrality FT0M (%)", kTH1F, {{110, 0, 110}}, false);
    registry.add("Event/before/hCentFT0MvsMultNTracksPV", "hCentFT0MvsMultNTracksPV;centrality FT0M (%);N_{track} to PV", kTH2F, {{110, 0, 110}, {600, 0, 6000}}, false);
    registry.add("Event/before/hMultFT0MvsMultNTracksPV", "hMultFT0MvsMultNTracksPV;mult. FT0M;N_{track} to PV", kTH2F, {{600, 0, 6000}, {600, 0, 6000}}, false);
    registry.addClone("Event/before/", "Event/after/");

    if (cfgPlotMBGeneral) {

      registry.add("MBGeneral", ";z_{conv} (cm); R_{conv} (cm); #eta; #varphi (rad); p_{T}", kTHnSparseF,
                   {{200, -100.f, 100.f},
                    {100, 0.f, 100.f},
                    {80, -2, +2},
                    {90, 0, o2::constants::math::TwoPI},
                    {1000, 0, 10}}, // pT 5
                   true);

      registry.add("Pi0/Same",
                   "Pi0;  m_{#gamma#gamma}; #it{p}_{T}; R_{conv} (cm); R_{conv} (cm)",
                   kTHnSparseF,
                   {{400, 0.f, 0.3},   // x 0
                    {500, 0, 10},      // pT 1
                    {100, 0.f, 100.f}, // R of photon 1
                    {100, 0.f, 100.f}},
                   true);

      registry.add("Dalitz/Same",
                   "Dalitz; m_{eeγ}; p_{T}; R_{γ} (cm); R_{e+}; R_{e-}",
                   kTHnSparseF,
                   {{400, 0.f, 0.5f},
                    {500, 0.f, 10.f},
                    {100, 0.f, 100.f}},
                   true);

      registry.add("Pi0/Mix", "Pi0 mixed-event;m_{#gamma#gamma};p_{T}", kTHnSparseF,
                   {{400, 0, 0.3}, {500, 0, 10}, {100, 0, 100}, {100, 0, 100}}, true);
    }

    if (cfgPlotMBDetailed) {

      registry.add("MBStudiesWireLeft",
                   "photon characteristics for MB studies; x (cm); y (cm); z (cm); #varphi (rad.);  R_{conv} (cm)",
                   kTHnSparseF,
                   {
                     {200, -20.f, 20.f}, // x 0
                     {200, -20.f, 20.f}, // y 1
                     {80, -20.f, 20.f},  // z 2
                     {40, 3.15, 3.4},    // phi  3
                     {80, 0.f, 20.f},    // conversion radius  4
                     {1000, 0, 10}       // pT 5
                   },
                   true);

      registry.add("MBStudiesWireRight",
                   "photon characteristics for MB studies; x (cm); y (cm); z (cm); #varphi (rad.);  R_{conv} (cm)",
                   kTHnSparseF,
                   {
                     {200, -20.f, 20.f}, // x 0
                     {200, -20.f, 20.f}, // y 1
                     {80, -20.f, 20.f},  // z 2
                     {40, 6.00, 6.15},   // phi  3
                     {80, 0.f, 20.f},    // conversion radius  4
                     {1000, 0, 10}       // pT 5
                   },
                   true);

      registry.add("MBStudiesWireITS",
                   "photon characteristics for MB studies; x (cm); y (cm); z (cm); #varphi (rad.);  R_{conv} (cm)",
                   kTHnSparseF,
                   {
                     {160, -40.f, 40.f}, // x 0
                     {160, -40.f, 40.f}, // y 1
                     {80, -40.f, 40.f},  // z 2
                     {90, 0, TwoPI},     // phi  3
                     {150, 0.f, 60.f},   // conversion radius  4
                     {1000, 0, 10}       // pT 5
                   },
                   true);

      registry.add("MBStudiesMFT",
                   "photon characteristics for MB studies; x (cm); y (cm); z (cm); #varphi (rad.);  R_{conv} (cm)",
                   kTHnSparseF,
                   {
                     {160, -80.f, 80.f}, // x 0
                     {160, -80.f, 80.f}, // y 1
                     {40, -40.f, 40.f},  // z 2
                     {90, 0, TwoPI},     // phi  3
                     {40, 40.f, 60.f},   // conversion radius  4
                     {500, 0, 10}        // pT 5
                   },
                   true);

      registry.add("MBStudiesTPC",
                   "photon characteristics for MB studies; x (cm); y (cm); z (cm); #varphi (rad.);  R_{conv} (cm)",
                   kTHnSparseF,
                   {
                     {160, -90.f, 90.f}, // x 0
                     {160, -90.f, 90.f}, // y 1
                     {80, -40.f, 40.f},  // z 2
                     {90, 0, TwoPI},     // phi  3
                     {80, 60.f, 80.f},   // conversion radius  4
                     {1000, 0, 10}       // pT 5
                   },
                   true);

      for (int nCls = 0; nCls <= 7; nCls++) { // o2-linter: disable=magic-number (just numbers for ITS cluster)
        registry.add(Form("ITSHits/Neg/ITSCls%d", nCls),
                     Form("Conversion radius NEG legs with %d ITS clusters;R_{conv} (cm);Entries", nCls),
                     HistType::kTH1F,
                     {AxisSpec{100, 0.f, 100.f}});

        registry.add(Form("ITSHits/Pos/ITSCls%d", nCls),
                     Form("Conversion radius POS legs with %d ITS clusters;R_{conv} (cm);Entries", nCls),
                     HistType::kTH1F,
                     {AxisSpec{100, 0.f, 100.f}});
      }
    }

    if (cfgPlotResolution) {

      registry.add("ResolutionGen/Z_res",
                   "Z resolution (gen ref);z_gen (cm);R_gen (cm);phi_gen;eta_gen;pT_gen;Δz (cm)",
                   kTHnSparseF,
                   {{200, -100, 100}, {200, 0, 100}, {90, 0, TwoPI}, {200, -1.0f, 1.0f}, {500, 0, 10}, {120, -30, 30}}, false);

      registry.add("ResolutionGen/R_res",
                   "R resolution (gen ref);z_gen (cm);R_gen (cm);phi_gen;eta_gen;pT_gen;ΔR (cm)",
                   kTHnSparseF,
                   {{200, -100, 100}, {200, 0, 100}, {90, 0, TwoPI}, {200, -1.0f, 1.0f}, {500, 0, 10}, {120, -30, 30}}, false);

      registry.add("ResolutionGen/Phi_res",
                   "Phi resolution (gen ref);z_gen (cm);R_gen (cm);phi_gen;eta_gen;pT_gen;Δφ",
                   kTHnSparseF,
                   {{200, -100, 100}, {200, 0, 100}, {90, 0, TwoPI}, {200, -1.0f, 1.0f}, {500, 0, 10}, {100, -0.2f, 0.2f}}, false);

      registry.add("ResolutionGen/Pt_res",
                   "pT resolution (gen ref);z_gen (cm);R_gen (cm);phi_gen;eta_gen;pT_gen;ΔpT/pT_gen",
                   kTHnSparseF,
                   {{200, -100, 100}, {200, 0, 100}, {90, 0, TwoPI}, {200, -1.0f, 1.0f}, {500, 0, 10}, {200, -1.0f, 1.0f}}, false);

      registry.add("ResolutionGen/Eta_res",
                   "Eta resolution (gen ref);z_gen (cm);R_gen (cm);phi_gen;eta_gen;pT_gen;Δη",
                   kTHnSparseF,
                   {{200, -100, 100}, {200, 0, 100}, {90, 0, TwoPI}, {200, -1.0f, 1.0f}, {500, 0, 10}, {100, -0.5f, 0.5f}}, false);

      registry.add("ResolutionRec/Z_res",
                   "Z resolution (rec ref);z_rec (cm);R_rec (cm);phi_rec;eta_rec;pT_rec;Δz (cm)",
                   kTHnSparseF,
                   {{200, -100, 100}, {200, 0, 100}, {90, 0, TwoPI}, {200, -1.0f, 1.0f}, {500, 0, 10}, {120, -30, 30}}, false);

      registry.add("ResolutionRec/R_res",
                   "R resolution (rec ref);z_rec (cm);R_rec (cm);phi_rec;eta_rec;pT_rec;ΔR (cm)",
                   kTHnSparseF,
                   {{200, -100, 100}, {200, 0, 100}, {90, 0, TwoPI}, {200, -1.0f, 1.0f}, {500, 0, 10}, {120, -30, 30}}, false);

      registry.add("ResolutionRec/Phi_res",
                   "Phi resolution (rec ref);z_rec (cm);R_rec (cm);phi_rec;eta_rec;pT_rec;Δφ",
                   kTHnSparseF,
                   {{200, -100, 100}, {200, 0, 100}, {90, 0, TwoPI}, {200, -1.0f, 1.0f}, {500, 0, 10}, {100, -0.2f, 0.2f}}, false);

      registry.add("ResolutionRec/Pt_res",
                   "pT resolution (rec ref);z_rec (cm);R_rec (cm);phi_rec;eta_rec;pT_rec;ΔpT/pT_gen",
                   kTHnSparseF,
                   {{200, -100, 100}, {200, 0, 100}, {90, 0, TwoPI}, {200, -1.0f, 1.0f}, {500, 0, 10}, {200, -1.0f, 1.0f}}, false);

      registry.add("ResolutionRec/Eta_res",
                   "Eta resolution (rec ref);z_rec (cm);R_rec (cm);phi_rec;eta_rec;pT_rec;Δη",
                   kTHnSparseF,
                   {{200, -100, 100}, {200, 0, 100}, {90, 0, TwoPI}, {200, -1.0f, 1.0f}, {500, 0, 10}, {100, -0.5f, 0.5f}}, false);
    }

    if (cfgPlotPurity) {

      registry.add("PhotonPurity",
                   "Photon purity classification;p_{T} (GeV/c);R_{conv} (cm);#eta;#varphi;category",
                   kTHnSparseF,
                   {
                     {500, 0, 10},                        // pT
                     {200, 0, 100},                       // R
                     {100, -1.0, 1.0},                    // eta
                     {90, 0, o2::constants::math::TwoPI}, // phi
                     {12, 0.5, 12.5}                      // purity categories
                   },
                   true);

      auto hPurity = registry.get<THnSparse>(HIST("PhotonPurity"));

      auto axis = hPurity->GetAxis(4);

      axis->SetBinLabel(1, "e^{+}e^{-}");
      axis->SetBinLabel(2, "e^{-} + #pi^{+}");
      axis->SetBinLabel(3, "#pi^{-} + e^{+}");
      axis->SetBinLabel(4, "K^{-} + e^{+}");
      axis->SetBinLabel(5, "e^{-} + K^{+}");
      axis->SetBinLabel(6, "#bar{p} + p");
      axis->SetBinLabel(7, "K^{-} + K^{+}");
      axis->SetBinLabel(8, "#pi^{-} + p");
      axis->SetBinLabel(9, "#bar{p} + #pi^{+}");
      axis->SetBinLabel(10, "#mu^{-} + #mu^{+}");
    }

    if (cfgPlotBremsstrahlung) {

      registry.add<TH2>("Bremsstrahlung/EnergyLossVsR", ";R (cm);E_{loss} (GeV)",
                        kTH2F, {{50, 0, 100}, {100, 0, 5}});
      registry.add<TH2>("Bremsstrahlung/EnergyLossXY", ";X (cm);Y (cm)",
                        kTH2F, {{100, -100, 100}, {100, -100, 100}});
      registry.add<TH2>("Bremsstrahlung/EnergyLossXYWeigh", ";X (cm);Y (cm)",
                        kTH2F, {{100, -100, 100}, {100, -100, 100}});
      registry.add<TH2>("Bremsstrahlung/EnergyLossVsZ", ";Z (cm);E_{loss} (GeV)",
                        kTH2F, {{100, -100, 100}, {200, 0, 5}});

      registry.add<TH1>("Bremsstrahlung/NBrem", ";N_{Brem} per e^{-};Counts",
                        kTH1I, {{10, 0, 10}});

      registry.add<TH2>("Bremsstrahlung/DeltaPoverPvsR", ";R (cm);(p_{Reco}-p_{MC})/p_{MC}",
                        kTH2F, {{50, 0, 100}, {200, -1, 1}});
      registry.add<TH2>("Bremsstrahlung/DeltaPtvsR", ";R (cm);(p_{T,Reco}-p_{T,MC})/p_{T,MC}",
                        kTH2F, {{50, 0, 100}, {300, -2, 2}});
      registry.add<TH2>("Bremsstrahlung/DeltaPhivsR", ";R (cm);#Delta#varphi (rad)",
                        kTH2F, {{50, 0, 100}, {200, -0.1, 0.1}});
      registry.add<TH2>("Bremsstrahlung/DeltaEtavsR", ";R (cm);#Delta#eta",
                        kTH2F, {{50, 0, 100}, {200, -0.1, 0.1}});

      registry.add<TH2>("Bremsstrahlung/Sigma1PtVsR", ";R (cm);#sigma(1/p_{T}) [c/GeV]",
                        kTH2F, {{50, 0, 100}, {100, 0, 0.1}});

      registry.add<TH2>("Bremsstrahlung/SigmaYVsR", ";R (cm);#sigma(y) [cm]",
                        kTH2F, {{50, 0, 100}, {100, 0, 0.5}});
      registry.add<TH2>("Bremsstrahlung/SigmaZVsR", ";R (cm);#sigma(z) [cm]",
                        kTH2F, {{50, 0, 100}, {100, 0, 1.0}});

      registry.add<TProfile>("Bremsstrahlung/relativeResoPtWBrems",
                             "relative #it{p}_{T} resolution; #it{p}_{T}; #sigma(#it{p}_{T})/#it{p}_{T}",
                             kTProfile,
                             {axisPt});

      registry.add<TProfile>("Bremsstrahlung/relativeResoPtWOBrems",
                             "relative #it{p}_{T} resolution; #it{p}_{T}; #sigma(#it{p}_{T})/#it{p}_{T}",
                             kTProfile,
                             {axisPt});
    }

    if (cfgPlotMBCollisions) {

      registry.add("RightCollisions/SparseDeltaCol",
                   "RightCollisions;R_{rec} (cm);#it{p}_{T} (GeV/c);#Delta #it{p}_{T}/#it{p}_{T,gen};#Delta z (cm);#Delta R (cm);#Delta col ID",
                   kTHnSparseF,
                   {
                     {100, 0.f, 100.f},  // R_rec
                     {200, 0.f, 10.f},   // pT
                     {200, -1.f, 1.f},   // ΔpT/pT
                     {200, -20.f, 20.f}, // Δz
                     {200, -8.f, 8.f},   // ΔR
                     {21, -10.5, 10.5}   // ΔCollisionID
                   },
                   true);

      registry.add("WrongCollisions/SparseDeltaCol",
                   "WrongCollisions;R_{rec} (cm);#it{p}_{T} (GeV/c);#Delta #it{p}_{T}/#it{p}_{T,gen};#Delta z (cm);#Delta R (cm);#Delta col ID",
                   kTHnSparseF,
                   {{100, 0.f, 100.f},
                    {200, 0.f, 10.f},
                    {200, -1.f, 1.f},
                    {200, -20.f, 20.f},
                    {200, -8.f, 8.f},
                    {21, -10.5, 10.5}},
                   true);
    }
  }

  template <int N>
  void fillITSClsNeg(HistogramRegistry& reg, float value)
  {
    if constexpr (N == 0) // o2-linter: disable=magic-number (just numbers for ITS cluster)
      reg.fill(HIST("ITSHits/Neg/ITSCls0"), value);
    if constexpr (N == 1) // o2-linter: disable=magic-number (just numbers for ITS cluster)
      reg.fill(HIST("ITSHits/Neg/ITSCls1"), value);
    if constexpr (N == 2) // o2-linter: disable=magic-number (just numbers for ITS cluster)
      reg.fill(HIST("ITSHits/Neg/ITSCls2"), value);
    if constexpr (N == 3) // o2-linter: disable=magic-number (just numbers for ITS cluster)
      reg.fill(HIST("ITSHits/Neg/ITSCls3"), value);
    if constexpr (N == 4) // o2-linter: disable=magic-number (just numbers for ITS cluster)
      reg.fill(HIST("ITSHits/Neg/ITSCls4"), value);
    if constexpr (N == 5) // o2-linter: disable=magic-number (just numbers for ITS cluster)
      reg.fill(HIST("ITSHits/Neg/ITSCls5"), value);
    if constexpr (N == 6) // o2-linter: disable=magic-number (just numbers for ITS cluster)
      reg.fill(HIST("ITSHits/Neg/ITSCls6"), value);
    if constexpr (N == 7) // o2-linter: disable=magic-number (just numbers for ITS cluster)
      reg.fill(HIST("ITSHits/Neg/ITSCls7"), value);
  }

  template <int N>
  void fillITSClsPos(HistogramRegistry& reg, float value)
  {
    if constexpr (N == 0) // o2-linter: disable=magic-number (just numbers for ITS cluster)
      reg.fill(HIST("ITSHits/Pos/ITSCls0"), value);
    if constexpr (N == 1) // o2-linter: disable=magic-number (just numbers for ITS cluster)
      reg.fill(HIST("ITSHits/Pos/ITSCls1"), value);
    if constexpr (N == 2) // o2-linter: disable=magic-number (just numbers for ITS cluster)
      reg.fill(HIST("ITSHits/Pos/ITSCls2"), value);
    if constexpr (N == 3) // o2-linter: disable=magic-number (just numbers for ITS cluster)
      reg.fill(HIST("ITSHits/Pos/ITSCls3"), value);
    if constexpr (N == 4) // o2-linter: disable=magic-number (just numbers for ITS cluster)
      reg.fill(HIST("ITSHits/Pos/ITSCls4"), value);
    if constexpr (N == 5) // o2-linter: disable=magic-number (just numbers for ITS cluster)
      reg.fill(HIST("ITSHits/Pos/ITSCls5"), value);
    if constexpr (N == 6) // o2-linter: disable=magic-number (just numbers for ITS cluster)
      reg.fill(HIST("ITSHits/Pos/ITSCls6"), value);
    if constexpr (N == 7) // o2-linter: disable=magic-number (just numbers for ITS cluster)
      reg.fill(HIST("ITSHits/Pos/ITSCls7"), value);
  }

  template <typename TV0s>
  void reconstructMesons(const TV0s& v0s)
  {
    for (auto i = 0; i < v0s.size(); i++) {
      auto g1 = v0s.iteratorAt(i);

      for (auto j = i + 1; j < v0s.size(); j++) {
        auto g2 = v0s.iteratorAt(j);

        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

        registry.fill(HIST("Pi0/Same"), v12.M(), v12.Pt(),
                      g1.v0radius(), g2.v0radius());
      }
    }
  }

  void fillClusterNeg(HistogramRegistry& reg, int nCls, float value)
  {
    switch (nCls) {
      case 0: // o2-linter: disable=magic-number (just numbers for ITS cluster)
        fillITSClsNeg<0>(reg, value);
        break;
      case 1: // o2-linter: disable=magic-number (just numbers for ITS cluster)
        fillITSClsNeg<1>(reg, value);
        break;
      case 2: // o2-linter: disable=magic-number (just numbers for ITS cluster)
        fillITSClsNeg<2>(reg, value);
        break;
      case 3: // o2-linter: disable=magic-number (just numbers for ITS cluster)
        fillITSClsNeg<3>(reg, value);
        break;
      case 4: // o2-linter: disable=magic-number (just numbers for ITS cluster)
        fillITSClsNeg<4>(reg, value);
        break;
      case 5: // o2-linter: disable=magic-number (just numbers for ITS cluster)
        fillITSClsNeg<5>(reg, value);
        break;
      case 6: // o2-linter: disable=magic-number (just numbers for ITS cluster)
        fillITSClsNeg<6>(reg, value);
        break;
      case 7: // o2-linter: disable=magic-number (just numbers for ITS cluster)
        fillITSClsNeg<7>(reg, value);
        break;
      default:
        break;
    }
  }

  void fillClusterPos(HistogramRegistry& reg, int nCls, float value)
  {
    switch (nCls) {
      case 0:
        fillITSClsPos<0>(reg, value);
        break;
      case 1:
        fillITSClsPos<1>(reg, value);
        break;
      case 2:
        fillITSClsPos<2>(reg, value);
        break;
      case 3:
        fillITSClsPos<3>(reg, value);
        break;
      case 4:
        fillITSClsPos<4>(reg, value);
        break;
      case 5:
        fillITSClsPos<5>(reg, value);
        break;
      case 6:
        fillITSClsPos<6>(reg, value);
        break;
      case 7:
        fillITSClsPos<7>(reg, value);
        break;
      default:
        break;
    }
  }
  template <const int evID, typename TCollision>
  void fillEventInfo(TCollision const& collision, const float /*weight*/ = 1.f)
  {
    registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hCollisionCounter"), 1.0);
    if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hCollisionCounter"), 2.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hCollisionCounter"), 3.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hCollisionCounter"), 4.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hCollisionCounter"), 5.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hCollisionCounter"), 6.0);
    }
    if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hCollisionCounter"), 7.0);
    }
    if (collision.sel8()) {
      registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hCollisionCounter"), 8.0);
    }
    if (std::fabs(collision.posZ()) < 10.0) { // o2-linter: disable=magic-number (vertex position)
      registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hCollisionCounter"), 9.0);
    }

    registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hZvtx"), collision.posZ());

    registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hMultNTracksPV"), collision.multNTracksPV());
    registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hMultNTracksPVeta1"), collision.multNTracksPVeta1());
    registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hMultFT0"), collision.multFT0A(), collision.multFT0C());
    registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hCentFT0A"), collision.centFT0A());
    registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hCentFT0C"), collision.centFT0C());
    registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hCentFT0M"), collision.centFT0M());
    registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hCentFT0MvsMultNTracksPV"), collision.centFT0M(), collision.multNTracksPV());
    registry.fill(HIST("Event/") + HIST(EventTypes[evID]) + HIST("hMultFT0MvsMultNTracksPV"), collision.multFT0A() + collision.multFT0C(), collision.multNTracksPV());
  }
  template <typename TV0s>
  void reconstructMesonsMixed(const TV0s& current,
                              const std::deque<std::vector<SimplePhoton>>& pool)
  {
    for (const auto& prev : pool) {
      for (const auto& g1 : current) {
        if (!fV0PhotonCut.IsSelected<aod::V0Legs>(g1)) {
          continue;
        }
        for (const auto& g2 : prev) {
          ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
          ROOT::Math::PtEtaPhiMVector v2(g2.pt, g2.eta, g2.phi, 0.);
          ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

          registry.fill(HIST("Pi0/Mix"), v12.M(), v12.Pt(),
                        g1.v0radius(), g2.r);
        }
      }
    }
  }

  template <typename TV0s, typename TPositrons, typename TElectrons>
  void reconstructDalitz(const TV0s& v0s,
                         const TPositrons& positronsPerCollision,
                         const TElectrons& electronsPerCollision,
                         MyCollision const& collision)
  {

    for (const auto& g1 : v0s) {
      if (!fV0PhotonCut.IsSelected<aod::V0Legs>(g1)) {
        continue;
      }
      ROOT::Math::PtEtaPhiMVector vGamma(g1.pt(), g1.eta(), g1.phi(), 0.);

      for (const auto& [pos, ele] :
           combinations(CombinationsFullIndexPolicy(positronsPerCollision,
                                                    electronsPerCollision))) {

        if (pos.trackId() == ele.trackId()) {
          continue;
        }

        if (!fDileptonCut.template IsSelectedTrack<false>(pos, collision) ||
            !fDileptonCut.template IsSelectedTrack<false>(ele, collision)) {
          continue; // Track-Cuts
        }

        ROOT::Math::PtEtaPhiMVector vPos(pos.pt(), pos.eta(), pos.phi(),
                                         o2::constants::physics::MassElectron);
        ROOT::Math::PtEtaPhiMVector vEle(ele.pt(), ele.eta(), ele.phi(),
                                         o2::constants::physics::MassElectron);

        auto vee = vPos + vEle;
        auto veGamma = vGamma + vee;

        registry.fill(HIST("Dalitz/Same"),
                      veGamma.M(), vee.Pt(),
                      g1.v0radius());
      }
    }
  }

  SliceCache cache;
  Preslice<MyV0Photons> perCollision = aod::v0photonkf::emeventId;
  Preslice<MyPrimaryElectrons> perCollisionElectron = aod::emprimaryelectron::emeventId;

  Partition<MyPrimaryElectrons> positrons = o2::aod::emprimaryelectron::sign > int8_t(0) && dileptoncuts.cfgMinPtTrack < o2::aod::track::pt&& nabs(o2::aod::track::eta) < dileptoncuts.cfgMaxEtaTrack;
  Partition<MyPrimaryElectrons> electrons = o2::aod::emprimaryelectron::sign < int8_t(0) && dileptoncuts.cfgMinPtTrack < o2::aod::track::pt && nabs(o2::aod::track::eta) < dileptoncuts.cfgMaxEtaTrack;

  Filter prefilterPrimaryelectron = ifnode(dileptoncuts.cfgApplyCutsFromPrefilterDerived.node(), o2::aod::emprimaryelectron::pfbderived == static_cast<uint16_t>(0), true);

  void processRec(MyCollisions const& collisions,
                  MyV0Photons const& v0photons,
                  aod::V0Legs const&,
                  MyPrimaryElectrons const& electrons)
  {

    for (const auto& collision : collisions) {
      const float centralities[3] = {collision.centFT0M(), collision.centFT0A(), collision.centFT0C()};
      if (centralities[cfgCentEstimator] < cfgCentMin || cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      fillEventInfo<0>(collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }
      fillEventInfo<1>(collision);
      registry.fill(HIST("Event/before/hCollisionCounter"), 10.0); // accepted
      registry.fill(HIST("Event/after/hCollisionCounter"), 10.0);  // accepted

      auto v0sThisCollision = v0photons.sliceBy(perCollision, collision.globalIndex());

      for (const auto& v0 : v0sThisCollision) {

        if (!fV0PhotonCut.IsSelected<aod::V0Legs>(v0)) {
          continue;
        }

        auto negLeg = v0.template negTrack_as<aod::V0Legs>();
        auto posLeg = v0.template posTrack_as<aod::V0Legs>();

        int nClsNeg = 0;
        int nClsPos = 0;

        if constexpr (requires { negLeg.itsNCls(); }) {
          nClsNeg = negLeg.itsNCls();
          nClsPos = posLeg.itsNCls();
        } else if constexpr (requires { negLeg.itsClusterMap(); }) {
          auto countBits = [](uint8_t map) {
            return std::bitset<8>(map).count();
          };
          nClsNeg = countBits(negLeg.itsClusterMap());
          nClsPos = countBits(posLeg.itsClusterMap());
        }
        if (nClsNeg >= 0 && nClsNeg <= 7) { // o2-linter: disable=magic-number (just numbers for ITS cluster)
          fillClusterNeg(registry, nClsNeg, v0.v0radius());
        }
        if (nClsPos >= 0 && nClsPos <= 7) { // o2-linter: disable=magic-number (just numbers for ITS cluster)
          fillClusterPos(registry, nClsPos, v0.v0radius());
        }

        if (cfgPlotMBGeneral) {

          registry.fill(
            HIST("MBGeneral"),
            v0.vz(),       // 0
            v0.v0radius(), // 1
            v0.eta(),      // 2
            v0.phi(),      // 3
            v0.pt());      // 4
        }

        if (cfgPlotMBDetailed) {

          registry.fill(
            HIST("MBStudiesWireLeft"),
            v0.vx(),       // 0
            v0.vy(),       // 1
            v0.vz(),       // 2
            v0.phi(),      // 3
            v0.v0radius(), // 4
            v0.pt());      // 5

          registry.fill(
            HIST("MBStudiesWireRight"),
            v0.vx(),       // 0
            v0.vy(),       // 1
            v0.vz(),       // 2
            v0.phi(),      // 3
            v0.v0radius(), // 4
            v0.pt());      // 5

          registry.fill(
            HIST("MBStudiesWireITS"),
            v0.vx(),       // 0
            v0.vy(),       // 1
            v0.vz(),       // 2
            v0.phi(),      // 3
            v0.v0radius(), // 4
            v0.pt());      // 5

          registry.fill(
            HIST("MBStudiesMFT"),
            v0.vx(),       // 0
            v0.vy(),       // 1
            v0.vz(),       // 2
            v0.phi(),      // 3
            v0.v0radius(), // 4
            v0.pt());      // 5

          registry.fill(
            HIST("MBStudiesTPC"),
            v0.vx(),       // 0
            v0.vy(),       // 1
            v0.vz(),       // 2
            v0.phi(),      // 3
            v0.v0radius(), // 4
            v0.pt());      // 5
        }
      }

      reconstructMesons(v0sThisCollision);

      auto electronsPerCollision = electrons.sliceBy(perCollisionElectron, collision.globalIndex());
      auto positronsPerCollision = positrons.sliceBy(perCollisionElectron, collision.globalIndex());

      reconstructDalitz(v0sThisCollision, positronsPerCollision, electronsPerCollision, collision);

      auto key = getPoolBin(centralities[cfgCentEstimator], collision.posZ());
      auto& pool = mixingPools[key];
      reconstructMesonsMixed(v0sThisCollision, pool);

      std::vector<SimplePhoton> eventCopy;
      eventCopy.reserve(v0sThisCollision.size());
      for (const auto& v0 : v0sThisCollision) {
        if (!fV0PhotonCut.IsSelected<aod::V0Legs>(v0)) {
          continue;
        }
        eventCopy.push_back({v0.pt(), v0.eta(), v0.phi(), v0.vz(), v0.v0radius()});
      }
      pool.emplace_back(std::move(eventCopy));

      if (pool.size() > MaxMixEvents) {
        pool.pop_front();
      }
    }
  }

  int classifyPurity(int pdgNeg, int pdgPos)
  {
    if (pdgNeg == kElectron && pdgPos == kPositron) {
      return 1;
    }
    if (pdgNeg == kElectron && pdgPos == kPiPlus)
      return 2;
    if (pdgNeg == kPiMinus && pdgPos == kPositron)
      return 3;
    if (pdgNeg == kKMinus && pdgPos == kPositron)
      return 4;
    if (pdgNeg == kElectron && pdgPos == kKPlus)
      return 5;
    if (pdgNeg == kProtonBar && pdgPos == kProton)
      return 6;
    if (pdgNeg == kKMinus && pdgPos == kKPlus)
      return 7;
    if (pdgNeg == kPiMinus && pdgPos == kProton)
      return 8;
    if (pdgNeg == kProtonBar && pdgPos == kPiPlus)
      return 9;
    if (pdgNeg == kMuonMinus && pdgPos == kMuonPlus)
      return 10;
    return 11;
  }

  Partition<MyCollisions> groupedCollisions = cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax;                                          // this goes to same event.
  Filter collisionFilterCommon = nabs(o2::aod::collision::posZ) < 10.f && o2::aod::collision::numContrib > static_cast<uint16_t>(0) && o2::aod::evsel::sel8 == true; // o2-linter: disable=magic-number (vertex position)
  Filter collisionFilterCentrality = (cfgCentMin < o2::aod::cent::centFT0M && o2::aod::cent::centFT0M < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0A && o2::aod::cent::centFT0A < cfgCentMax) || (cfgCentMin < o2::aod::cent::centFT0C && o2::aod::cent::centFT0C < cfgCentMax);
  using MyFilteredCollisions = soa::Filtered<MyCollisions>; // this goes to same event pairing.

  void processTrue(MyCollisions const&, MyFilteredCollisions const& filteredCollisions,
                   MyV0Photons const& v0photons,
                   MyMCV0Legs const&,
                   aod::EMMCParticles const& mcparticles,
                   aod::EMMCEvents const&)
  {
    for (const auto& collision : filteredCollisions) {

      const float centralities[3] = {
        collision.centFT0M(),
        collision.centFT0A(),
        collision.centFT0C()};

      if (centralities[cfgCentEstimator] < cfgCentMin ||
          cfgCentMax < centralities[cfgCentEstimator]) {
        continue;
      }

      fillEventInfo<0>(collision);
      if (!fEMEventCut.IsSelected(collision)) {
        continue;
      }

      fillEventInfo<1>(collision);

      registry.fill(HIST("Event/before/hCollisionCounter"), 10.0);
      registry.fill(HIST("Event/after/hCollisionCounter"), 10.0);

      auto v0PhotonsColl = v0photons.sliceBy(perCollision, collision.globalIndex());

      for (const auto& v0 : v0PhotonsColl) {

        auto pos = v0.posTrack_as<MyMCV0Legs>();
        auto ele = v0.negTrack_as<MyMCV0Legs>();

        auto posmc = pos.template emmcparticle_as<aod::EMMCParticles>();
        auto elemc = ele.template emmcparticle_as<aod::EMMCParticles>();

        if (!fV0PhotonCut.IsSelected<MyMCV0Legs>(v0)) {
          continue;
        }

        int purityCat = 12; // default: unmatched

        purityCat = classifyPurity(elemc.pdgCode(), posmc.pdgCode());

        float rConv = v0.v0radius();
        registry.fill(HIST("PhotonPurity"),
                      v0.pt(), rConv, v0.eta(), v0.phi(), purityCat);

        int photonid = FindCommonMotherFrom2Prongs(posmc, elemc, -11, 11, 22, mcparticles);
        if (photonid < 0) {
          continue;
        }

        auto mcphoton = mcparticles.iteratorAt(photonid);

        float rRec = v0.v0radius();
        float rGen = std::sqrt(mcphoton.vx() * mcphoton.vx() +
                               mcphoton.vy() * mcphoton.vy());
        float deltaR = rRec - rGen;
        float deltaZ = v0.vz() - mcphoton.vz();
        float deltaPhi = v0.phi() - mcphoton.phi();
        float deltaEta = v0.eta() - mcphoton.eta();
        float deltaPt = v0.pt() - mcphoton.pt();

        registry.fill(HIST("ResolutionGen/Z_res"),
                      mcphoton.vz(), rGen, mcphoton.phi(),
                      mcphoton.eta(), mcphoton.pt(), deltaZ);

        registry.fill(HIST("ResolutionGen/R_res"),
                      mcphoton.vz(), rGen, mcphoton.phi(),
                      mcphoton.eta(), mcphoton.pt(), deltaR);

        registry.fill(HIST("ResolutionGen/Phi_res"),
                      mcphoton.vz(), rGen, mcphoton.phi(),
                      mcphoton.eta(), mcphoton.pt(), deltaPhi);

        registry.fill(HIST("ResolutionGen/Eta_res"),
                      mcphoton.vz(), rGen, mcphoton.phi(),
                      mcphoton.eta(), mcphoton.pt(), deltaEta);

        registry.fill(HIST("ResolutionGen/Pt_res"),
                      mcphoton.vz(), rGen, mcphoton.phi(),
                      mcphoton.eta(), mcphoton.pt(),
                      deltaPt / mcphoton.pt());

        registry.fill(HIST("ResolutionRec/Z_res"),
                      v0.vz(), rRec, v0.phi(), v0.eta(),
                      v0.pt(), deltaZ);

        registry.fill(HIST("ResolutionRec/R_res"),
                      v0.vz(), rRec, v0.phi(), v0.eta(),
                      v0.pt(), deltaR);

        registry.fill(HIST("ResolutionRec/Phi_res"),
                      v0.vz(), rRec, v0.phi(), v0.eta(),
                      v0.pt(), deltaPhi);

        registry.fill(HIST("ResolutionRec/Eta_res"),
                      v0.vz(), rRec, v0.phi(), v0.eta(),
                      v0.pt(), deltaEta);

        registry.fill(HIST("ResolutionRec/Pt_res"),
                      v0.vz(), rRec, v0.phi(), v0.eta(),
                      v0.pt(), deltaPt / mcphoton.pt());
      }
    }
  }

  void processCollAssoc(MyCollisionsMC const& collisions,
                        MyMCCollisions const&,
                        aod::EMMCParticles const& mcparticles,
                        MyV0Photons const& v0photons,
                        MyMCV0Legs const&)
  {
    for (auto const& col : collisions) {

      auto mccollision = col.emmcevent_as<MyMCCollisions>();

      auto v0s = v0photons.sliceBy(perCollision, col.globalIndex());

      int mcColIdDerived = mccollision.globalIndex(); // MC event (derived) global index

      for (auto const& v0 : v0s) {

        if (!fV0PhotonCut.IsSelected<MyMCV0Legs>(v0)) {
          continue;
        }

        auto negLeg = v0.negTrack_as<MyMCV0Legs>();
        auto posLeg = v0.posTrack_as<MyMCV0Legs>();

        auto negMC = negLeg.template emmcparticle_as<aod::EMMCParticles>();
        auto posMC = posLeg.template emmcparticle_as<aod::EMMCParticles>();

        // find the photon mother in MC
        int photonId = FindCommonMotherFrom2Prongs(posMC, negMC, -11, 11, 22, mcparticles);
        if (photonId < 0) {
          continue;
        }

        auto mcPhoton = mcparticles.iteratorAt(photonId);
        auto trueMcCol = mcPhoton.emmcevent_as<MyMCCollisions>(); // MC event where the photon was generated
        int trueMcColIndex = trueMcCol.globalIndex();

        int deltaCol = mcColIdDerived - trueMcColIndex;

        float rRec = v0.v0radius();
        float rGen = std::sqrt(mcPhoton.vx() * mcPhoton.vx() + mcPhoton.vy() * mcPhoton.vy());
        float deltaR = rRec - rGen;
        float deltaZ = v0.vz() - mcPhoton.vz();
        float deltaPt = v0.pt() - mcPhoton.pt();
        float relPtRes = (mcPhoton.pt() > 0.f) ? deltaPt / mcPhoton.pt() : 0.f;

        // Fill base debug
        registry.fill(HIST("DeltaRecDerived"), deltaCol);
        registry.fill(HIST("AllR"), rRec);

        if (deltaCol == 0) {
          registry.fill(HIST("RightCollisions/SparseDeltaCol"),
                        rRec,
                        v0.pt(),
                        relPtRes,
                        deltaZ,
                        deltaR,
                        deltaCol);
        } else {
          registry.fill(HIST("WrongCollisions/SparseDeltaCol"),
                        rRec,
                        v0.pt(),
                        relPtRes,
                        deltaZ,
                        deltaR,
                        deltaCol);
        }

        if (deltaCol == 0) {
          registry.fill(HIST("RightCollisions/DeltaColvsZ"), v0.vz());
          registry.fill(HIST("RightCollisions/DeltaColvsR"), rRec);
          registry.fill(HIST("RightCollisions/DeltaColvspT"), v0.pt());
          registry.fill(HIST("RightCollisions/DeltaColvsZRes"), deltaZ);
          registry.fill(HIST("RightCollisions/DeltaColvsRRes"), deltaR);
          registry.fill(HIST("RightCollisions/DeltaColvspTRes"), relPtRes);
        } else {
          registry.fill(HIST("WrongCollisions/DeltaColvsZ"), v0.vz());
          registry.fill(HIST("WrongCollisions/DeltaColvsR"), rRec);
          registry.fill(HIST("WrongCollisions/DeltaColvspT"), v0.pt());
          registry.fill(HIST("WrongCollisions/DeltaColvsZRes"), deltaZ);
          registry.fill(HIST("WrongCollisions/DeltaColvsRRes"), deltaR);
          registry.fill(HIST("WrongCollisions/DeltaColvspTRes"), relPtRes);
        }
      }
    }
  }

  void processBremsstrahlung(MyTracks const& tracks, aod::McParticles const& mcparticles)
  {

    for (const auto& trk : tracks) {
      if (!trk.has_mcParticle()) {
        continue;
      }
      auto mc = trk.mcParticle();
      if (std::abs(mc.pdgCode()) != kElectron) {
        continue;
      }

      if (trk.pt() < dileptoncuts.cfgMinPtTrack) {
        continue;
      }

      double pMC = mc.p();
      double ptMC = mc.pt();
      double phiMC = mc.phi();
      double etaMC = mc.eta();

      double pReco = trk.p();
      double ptReco = trk.pt();
      double phiReco = trk.phi();
      double etaReco = trk.eta();

      int nBrem = 0; // o2-linter: disable=magic-number

      registry.fill(HIST("Bremsstrahlung/relativeResoPtWOBrems"), trk.pt(), trk.pt() * std::sqrt(trk.c1Pt21Pt2()));

      bool isFirst = true;
      for (const auto& dId : mc.daughtersIds()) {
        if (dId < 0 || dId >= mcparticles.size()) {
          continue;
        }
        auto daughter = mcparticles.iteratorAt(dId);
        if (daughter.getProcess() != kPBrem) {
          continue;
        }

        if (std::abs(daughter.eta()) > pcmcuts.cfgMaxEtaV0) {
          continue;
        }

        double r = std::hypot(daughter.vx(), daughter.vy());
        double z = daughter.vz();

        nBrem++;

        registry.fill(HIST("Bremsstrahlung/EnergyLossVsR"), r, daughter.e());
        registry.fill(HIST("Bremsstrahlung/EnergyLossVsZ"), z, daughter.e());
        registry.fill(HIST("Bremsstrahlung/EnergyLossXY"), daughter.vx(), daughter.vy());
        registry.fill(HIST("Bremsstrahlung/EnergyLossXYWeigh"), daughter.vx(), daughter.vy(), daughter.e());

        if (isFirst) {
          registry.fill(HIST("Bremsstrahlung/relativeResoPtWBrems"), trk.pt(), trk.pt() * std::sqrt(trk.c1Pt21Pt2()));
        }

        registry.fill(HIST("Bremsstrahlung/Sigma1PtVsR"), r, trk.sigma1Pt());

        registry.fill(HIST("Bremsstrahlung/SigmaYVsR"), r, trk.sigmaY());
        registry.fill(HIST("Bremsstrahlung/SigmaZVsR"), r, trk.sigmaZ());

        if (pMC > 0) {
          registry.fill(HIST("Bremsstrahlung/DeltaPoverPvsR"), r, (pReco - pMC) / pMC);
        }
        if (ptMC > 0) {
          registry.fill(HIST("Bremsstrahlung/DeltaPtvsR"), r, (ptReco - ptMC) / ptMC);
        }
        registry.fill(HIST("Bremsstrahlung/DeltaPhivsR"), r, phiReco - phiMC);
        registry.fill(HIST("Bremsstrahlung/DeltaEtavsR"), r, etaReco - etaMC);
      }

      registry.fill(HIST("Bremsstrahlung/NBrem"), nBrem);
      isFirst = false;
    }
  }

  PROCESS_SWITCH(MaterialBudget, processRec, "process material budget", true);
  PROCESS_SWITCH(MaterialBudget, processTrue, "process material budget", false);
  PROCESS_SWITCH(MaterialBudget, processCollAssoc, "process material budget", false);
  PROCESS_SWITCH(MaterialBudget, processBremsstrahlung, "process material budget", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MaterialBudget>(cfgc)};
}
