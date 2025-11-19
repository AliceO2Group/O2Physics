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

/// \file phosElId.cxx
/// \struct PHOS electron id analysis
/// \brief Task for calculating electron identification parameters
///
/// \author Yeghishe Hambardzumyan, MIPT
/// \since Apr, 2024

#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CaloClusters.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "PHOSBase/Geometry.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

#include "TF1.h"
#include "TPDGCode.h"

#include <climits>
#include <cmath>
#include <map>
#include <memory>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::aod::evsel;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace phos_match
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_INDEX_COLUMN(Track, track);
DECLARE_SOA_INDEX_COLUMN(CaloCluster, caloCluster);
DECLARE_SOA_COLUMN(SignalCharge, signalCharge, float);
DECLARE_SOA_COLUMN(SignalZ, signalZ, float);
DECLARE_SOA_COLUMN(SignalPx, signalPx, double);
DECLARE_SOA_COLUMN(SignalPy, signalPy, double);
DECLARE_SOA_COLUMN(SignalPz, signalPz, double);
DECLARE_SOA_COLUMN(SignalE, signalE, double);
} // namespace phos_match

DECLARE_SOA_TABLE(PHOSMatchindexTable, "AOD", "PHSMTCH",                                                      //!
                  o2::soa::Index<>, phos_match::CollisionId, phos_match::CaloClusterId, phos_match::TrackId); //!

} // namespace o2::aod

// globalized estimator names for centrality
enum CentEstimators { FV0A,
                      FT0M,
                      FT0A,
                      FT0C,
                      FDDM,
                      NTPV };

bool testLambda(float pt, float l1, float l2, float cutThreshold, bool useNegativeCrossTerm)
{
  float l2Mean = 1.53126f + 9.50835e+06f / (1.f + 1.08728e+07f * pt + 1.73420e+06f * pt * pt);
  float l1Mean = 1.12365f + 0.123770f * std::exp(-pt * 0.246551f) + 5.30000e-03f * pt;
  float l2Sigma = 6.48260e-02f + 7.60261e+10f / (1.f + 1.53012e+11f * pt + 5.01265e+05f * pt * pt) + 9.00000e-03f * pt;
  float l1Sigma = 4.44719e-04f + 6.99839e-01f / (1.f + 1.22497e+00f * pt + 6.78604e-07f * pt * pt) + 9.00000e-03f * pt;
  float c = -0.35f - 0.550f * std::exp(-0.390730f * pt);
  if (l1Sigma == 0.f || l2Sigma == 0.f)
    return false;

  float term1 = 0.5f * (l1 - l1Mean) * (l1 - l1Mean) / (l1Sigma * l1Sigma);
  float term2 = 0.5f * (l2 - l2Mean) * (l2 - l2Mean) / (l2Sigma * l2Sigma);
  float crossTerm = 0.5f * c * (l1 - l1Mean) * (l2 - l2Mean) / (l1Sigma * l2Sigma);

  float rSquared;
  if (useNegativeCrossTerm) {
    rSquared = term1 + term2 - crossTerm;
  } else {
    rSquared = term1 + term2 + crossTerm;
  }

  return rSquared < cutThreshold;
}

struct PhosElId {

  Produces<o2::aod::PHOSMatchindexTable> phosMatch;

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                             aod::pidTOFFullEl, aod::pidTPCFullEl, aod::pidTPCFullPi,
                             aod::pidTPCFullKa, aod::pidTPCFullPr>;
  Configurable<bool> isMC{"isMC", true, "Enable MC analysis"},
    isSel8{"isSel8", 1, "check if event is Single Event Latch-up 8"},
    mSwapM20M02ForTestLambda{"mSwapM20M02ForTestLambda", false, "Swap m20 and m02 arguments for testLambda (false for note's correct order, true for swapped/original incorrect order)"},
    mUseNegativeCrossTerm{"mUseNegativeCrossTerm", true, "Use negative sign for the cross-term in testLambda (true for analysis note version, false for old version)"},
    mFillSingleLoopHistos{"mFillSingleLoopHistos", false, "Fill single loop histograms"};

  Configurable<float> mColMaxZ{"mColMaxZ", 10.f, "maximum z accepted in analysis"},
    mMinCluE{"mMinCluE", 0.3, "Minimum cluster energy for analysis"},
    mMinCluTime{"minCluTime", -25.e-9, "Min. cluster time"},
    mMaxCluTime{"mMaxCluTime", 25.e-9, "Max. cluster time"},
    mCluTimeAxisMin{"mCluTimeAxisMin", -100, "lower axis limit for cluster time in nanoseconds"},
    mCluTimeAxisMax{"mCluTimeAxisMax", 100, "upper axis limit for cluster time in nanoseconds"},
    mDeltaXmin{"mDeltaXmin", -100., "Min for track and cluster coordinate delta"},
    mDeltaXmax{"mDeltaXmax", 100., "Max for track and cluster coordinate delta"},
    mDeltaZmin{"mDeltaZmin", -100., "Min for track and cluster coordinate delta"},
    mDeltaZmax{"mDeltaZmax", 100., "Max for track and cluster coordinate delta"},
    mEpmin{"mEpmin", -1., "Min for E/p histograms"},
    mEpmax{"mEpmax", 3., "Max for E/p histograms"},
    EtaMax{"EtaMax", {0.8f}, "eta ranges"},
    PtMin{"PtMin", {0.2f}, "pt min"},
    PtMax{"PtMax", {20.f}, "pt max"},
    DCAxyMax{"DCAxyMax", {3.f}, "dcaxy max"},
    DCAzMax{"DCAzMax", {3.f}, "dcaz max"},
    ITSchi2Max{"ITSchi2Max", {5.f}, "its chi2 max"},
    ITSnclsMin{"ITSnclsMin", {4.5f}, "min number of ITS clusters"},
    ITSnclsMax{"ITSnclsMax", {7.5f}, "max number of ITS clusters"},
    TPCchi2Max{"TPCchi2Max", {4.f}, "tpc chi2 max"},
    TPCnclsMin{"TPCnclsMin", {90.f}, "min number of TPC clusters"},
    TPCnclsMax{"TPCnclsMax", {170.f}, "max number of TPC clusters"},
    TPCnclsCRMin{"TPCnclsCRMin", {80.f}, "min number of TPC crossed rows"},
    TPCnclsCRMax{"TPCnclsCRMax", {161.f}, "max number of TPC crossed rows"},
    TPCNSigmaElMin{"TPCNSigmaElMin", {-3.f}, "min TPC nsigma e for inclusion"},
    TPCNSigmaElMax{"TPCNSigmaElMax", {2.f}, "max TPC nsigma e for inclusion"},
    TPCNSigmaPiMin{"TPCNSigmaPiMin", {-3.f}, "min TPC nsigma pion for exclusion"},
    TPCNSigmaPiMax{"TPCNSigmaPiMax", {3.5f}, "max TPC nsigma pion for exclusion"},
    TPCNSigmaPrMin{"TPCNSigmaPrMin", {-3.f}, "min TPC nsigma proton for exclusion"},
    TPCNSigmaPrMax{"TPCNSigmaPrMax", {4.f}, "max TPC nsigma proton for exclusion"},
    TPCNSigmaKaMin{"TPCNSigmaKaMin", {-3.f}, "min TPC nsigma kaon for exclusion"},
    TPCNSigmaKaMax{"TPCNSigmaKaMax", {4.f}, "max TPC nsigma kaon for exclusion"},
    TOFNSigmaElMin{"TOFNSigmaElMin", {-3.f}, "min TOF nsigma e for inclusion"},
    TOFNSigmaElMax{"TOFNSigmaElMax", {3.f}, "max TOF nsigma e for inclusion"},
    NsigmaTrackMatch{"NsigmaTrackMatch", {2.f}, "PHOS Track Matching Nsigma for inclusion"},
    mShowerShapeCutValue{"mShowerShapeCutValue", 4.f, "Cut threshold for testLambda shower shape"};

  Configurable<int> mEvSelTrig{"mEvSelTrig", kTVXinPHOS, "Select events with this trigger"},
    mAmountOfModules{"mAmountOfModules", 4, "amount of modules for PHOS"},
    mMinCluNcell{"minCluNcell", 3, "min cells in cluster"},
    nBinsDeltaX{"nBinsDeltaX", 500, "N bins for track and cluster coordinate delta"},
    nBinsDeltaZ{"nBinsDeltaZ", 500, "N bins for track and cluster coordinate delta"},
    nBinsEp{"nBinsEp", 400, "N bins for E/p histograms"};

  Configurable<std::vector<float>> pSigmadZ{"pSigmadZ", {0.642, 0., 1.77, 2.725, 0.}, "parameters for sigma dz function"},
    pSigmadX{"pSigmadX", {2.17769, 1.60275, 2.24136}, "parameters for sigma dx function"},
    pPhosShiftZ{"pPhosShiftZ", {4.78838, 2.75138, 1.40825, 2.28735}, "Phos coordinate centering Z per module"},
    pPhosShiftX{"pPhosShiftX", {2.158702, -1.526772, -0.814658, -1.852678}, "Phos coordinate centering X per module"},
    pMeandXPosMod{"pMeandXPosMod", {-10.57, -0.42, 1.06, -8.1, -0.42, 1.14, -8.34, -0.42, 1.04, -7.38, -0.42, 1.17}, "parameters for mean dx function for positive tracks"},
    pMeandXNegMod{"pMeandXNegMod", {9.92, -0.42, 1.29, 7.82, -0.4, 1.34, 8.45, -0.33, 1.5, 7.5, -0.42, 1.25}, "parameters for mean dx function for negative tracks"};

  Filter ptFilter = (aod::track::pt > PtMin) && (aod::track::pt < PtMax),
         etafilter = nabs(aod::track::eta) < EtaMax,
         dcaxyfilter = nabs(aod::track::dcaXY) < DCAxyMax,
         dcazfilter = nabs(aod::track::dcaZ) < DCAzMax,
         itschi2filter = aod::track::itsChi2NCl < ITSchi2Max,
         tpcchi2filter = aod::track::tpcChi2NCl < TPCchi2Max,
         mapfilter = (aod::track::itsClusterMap & uint8_t(1)) > 0;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  std::unique_ptr<o2::phos::Geometry> geomPHOS;
  double bz{0.}; // magnetic field
  int runNumber{0};

  HistogramRegistry mHistManager{"PhosElIdHistograms"};
  TF1 *fSigma_dz, *fSigma_dx;
  std::array<TF1*, 4> fMeandXPosMod, fMeandXNegMod;

  void init(InitContext const&)
  {
    LOG(info) << "Initializing PHOS electron identification analysis task ...";

    std::vector<double> momentumBinning = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                                           1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                                           4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 10.};

    const AxisSpec axisCounter{3, 0, +3, ""},
      axisP{momentumBinning, "p (GeV/c)"},
      axisPt{momentumBinning, "p_{T} (GeV/c)"},
      axisEta{200, -0.2, 0.2, "#eta"},
      axisPhi{80, 240, 320, "#varphi"},
      axisE{200, 0, 10, "E (GeV)", "E (GeV)"},
      axisMassSpectrum{4000, 0, 4, "M (GeV/c^{2})", "Mass e^{+}e^{-} (GeV/c^{2})"},
      axisEp{nBinsEp, mEpmin, mEpmax, "E/p", "E_{cluster}/p_{track}"},
      axisdX{nBinsDeltaX, mDeltaXmin, mDeltaXmax, "x_{tr}-x_{clu} (cm)", "x_{tr}-x_{clu} (cm)"},
      axisdZ{nBinsDeltaZ, mDeltaZmin, mDeltaZmax, "z_{tr}-z_{clu} (cm)", "z_{tr}-z_{clu} (cm)"},
      axisCells{20, 0., 20., "number of cells", "number of cells"},
      axisTime{200, mCluTimeAxisMin, mCluTimeAxisMax, "time (ns)", "time (nanoseconds)"},
      axisModes{4, 1., 5., "module", "module"},
      axisX{150, -75., 75., "x (cm)", "x (cm)"},
      axisZ{150, -75., 75., "z (cm)", "z (cm)"},
      axisDCATrackXY{400, -.1, .1, "DCA XY (cm)", "DCA XY (cm)"},
      axisDCATrackZ{400, -.1, .1, "DCA Z (cm)", "DCA Z (cm)"},
      axisVColX{400, -.1, .1, "colision vertex x (cm)", "colision vertex x (cm)"},
      axisVColY{400, -.1, .1, "colision vertex y (cm)", "colision vertex y (cm)"},
      axisVColZ{400, -20., 20., "colision vertex z (cm)", "colision vertex z (cm)"}, // should look like gauss
      axisVTrackX{400, -10., 10., "track vertex x (cm)", "track vertex x (cm)"},
      axisVTrackY{400, -10., 10., "track vertex y (cm)", "track vertex y (cm)"},
      axisVTrackZ{400, -10., 10., "track vertex z (cm)", "track vertex z (cm)"};

    mHistManager.add("eventCounter", "eventCounter", kTH1F, {axisCounter});

    mHistManager.add("tracks/hTrackPtEtaPhi", "Track pt vs eta vs phi", HistType::kTH3F, {axisPt, axisEta, axisPhi});
    mHistManager.add("tracks/hTrackPtEtaPhi_Phos", "Track pt vs eta vs phi on Phos surface", HistType::kTH3F, {axisPt, axisEta, axisPhi});
    mHistManager.add("tracks/hTrackDCA", "Track DCA info", HistType::kTH2F, {axisDCATrackXY, axisDCATrackZ});
    mHistManager.add("tracks/hTrackVX", "Track vertex coordinate X", HistType::kTH1F, {axisVTrackX});
    mHistManager.add("tracks/hTrackVY", "Track vertex coordinate Y", HistType::kTH1F, {axisVTrackY});
    mHistManager.add("tracks/hTrackVZ", "Track vertex coordinate Z", HistType::kTH1F, {axisVTrackZ});
    mHistManager.add("collision/hColVX", "Collision vertex coordinate X", HistType::kTH1F, {axisVColX});
    mHistManager.add("collision/hColVY", "Collision vertex coordinate Y", HistType::kTH1F, {axisVColY});
    mHistManager.add("collision/hColVZ", "Collision vertex coordinate Z", HistType::kTH1F, {axisVColZ});
    mHistManager.add("tracks/hTrackPhosProjMod", "Track projection coordinates on PHOS modules", HistType::kTH3F, {axisX, axisZ, axisModes});

    mHistManager.add("clusterSpectra/hCluE_v_mod_v_time", "Cluster energy spectrum (E > 0.3 GeV) vs time per module", HistType::kTH3F, {axisE, axisTime, axisModes});
    mHistManager.add("clusterSpectra/hCluE_mod_energy_cut", "Cluster energy spectrum (E > 0.3 GeV) per module", HistType::kTH2F, {axisE, axisModes});
    mHistManager.add("clusterSpectra/hCluE_mod_time_cut", "Cluster energy spectrum (E > 0.3 GeV)(time +-25 ns) per module", HistType::kTH2F, {axisE, axisModes});
    mHistManager.add("clusterSpectra/hCluE_mod_cell_cut", "Cluster energy spectrum (E > 0.3 GeV)(time +-25 ns)(ncells > 3) per module", HistType::kTH2F, {axisE, axisModes});
    mHistManager.add("clusterSpectra/hCluE_mod_disp", "Cluster energy spectrum OK dispersion and (E > 0.3 GeV)(time +-25 ns)(ncells > 3) per module", HistType::kTH2F, {axisE, axisModes});
    mHistManager.add("clusterSpectra/hCluE_ncells_mod", "Cluster energy spectrum vs cell ammount per module", HistType::kTH3F, {axisE, axisCells, axisModes});

    mHistManager.add("coordinateMatching/hCluXZ_mod", "Local cluster X Z per module", HistType::kTH3F, {axisX, axisZ, axisModes});
    mHistManager.add("coordinateMatching/hdZpmod", "dz,p_{tr},module", HistType::kTH3F, {axisdZ, axisPt, axisModes});
    mHistManager.add("coordinateMatching/hdZpmod_pos", "dz,p_{tr},module positive tracks", HistType::kTH3F, {axisdZ, axisPt, axisModes});
    mHistManager.add("coordinateMatching/hdZpmod_neg", "dz,p_{tr},module negative tracks", HistType::kTH3F, {axisdZ, axisPt, axisModes});
    mHistManager.add("coordinateMatching/hdXpmod", "dx,p_{tr},module", HistType::kTH3F, {axisdX, axisPt, axisModes});
    mHistManager.add("coordinateMatching/hdXpmod_pos", "dx,p_{tr},module positive tracks", HistType::kTH3F, {axisdX, axisPt, axisModes});
    mHistManager.add("coordinateMatching/hdXpmod_neg", "dx,p_{tr},module negative tracks", HistType::kTH3F, {axisdX, axisPt, axisModes});

    mHistManager.add("clusterSpectra/hCluE_v_pt_disp", "Cluster energy vs p | OK dispersion", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("clusterSpectra/hCluE_v_pt_Nsigma", "Cluster energy vs p within trackmatch Nsigma", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("clusterSpectra/hCluE_v_pt_Nsigma_disp", "Cluster energy vs p within trackmatch Nsigma | OK dispersion", HistType::kTH3F, {axisE, axisPt, axisModes});

    mHistManager.add("energyMomentumRatio/hEp_v_pt_disp", "E/p ratio vs p | OK dispersion", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("energyMomentumRatio/hEp_v_pt_Nsigma", "E/p ratio vs p within trackmatch Nsigma", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("energyMomentumRatio/hEp_v_pt_Nsigma_disp", "E/p ratio vs p within trackmatch Nsigma | OK dispersion", HistType::kTH3F, {axisEp, axisPt, axisModes});

    mHistManager.add("energyMomentumRatio/hEp_v_E_disp", "E/p ratio vs cluster E | OK dispersion", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("energyMomentumRatio/hEp_v_E_Nsigma", "E/p ratio vs cluster E within trackmatch Nsigma", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("energyMomentumRatio/hEp_v_E_Nsigma_disp", "E/p ratio vs cluster E within trackmatch Nsigma | OK dispersion", HistType::kTH3F, {axisEp, axisE, axisModes});

    mHistManager.add("doubleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma", "Cluster energy vs p within trackmatch Nsigma", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("doubleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma_disp", "Cluster energy vs p within trackmatch Nsigma | OK dispersion", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("doubleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma_TPCel", "Cluster energy vs p within trackmatch Nsigma", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("doubleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma_disp_TPCel", "Cluster energy vs p within trackmatch Nsigma | OK dispersion", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("doubleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma", "E/p ratio vs p within trackmatch Nsigma", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("doubleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma_disp", "E/p ratio vs p within trackmatch Nsigma | OK dispersion", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("doubleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma_TPCel", "E/p ratio vs p within trackmatch Nsigma", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("doubleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma_disp_TPCel", "E/p ratio vs p within trackmatch Nsigma | OK dispersion", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("doubleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma", "E/p ratio vs cluster E within trackmatch Nsigma", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("doubleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma_disp", "E/p ratio vs cluster E within trackmatch Nsigma | OK dispersion", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("doubleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma_TPCel", "E/p ratio vs cluster E within trackmatch Nsigma", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("doubleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma_disp_TPCel", "E/p ratio vs cluster E within trackmatch Nsigma | OK dispersion", HistType::kTH3F, {axisEp, axisE, axisModes});

    const char* commonHistos[][2] = {
      {"hCluE_v_pt_disp", "Cluster energy vs p | OK dispersion"},
      {"hCluE_v_pt_Nsigma", "Cluster energy vs p within trackmatch Nsigma"},
      {"hCluE_v_pt_Nsigma_disp", "Cluster energy vs p within trackmatch Nsigma | OK dispersion"},
      {"hEp_v_pt_disp", "E/p ratio vs p | OK dispersion"},
      {"hEp_v_pt_Nsigma", "E/p ratio vs p within trackmatch Nsigma"},
      {"hEp_v_pt_Nsigma_disp", "E/p ratio vs p within trackmatch Nsigma | OK dispersion"},
      {"hEp_v_E_disp", "E/p ratio vs cluster E | OK dispersion"},
      {"hEp_v_E_Nsigma", "E/p ratio vs cluster E within trackmatch Nsigma"},
      {"hEp_v_E_Nsigma_disp", "E/p ratio vs cluster E within trackmatch Nsigma | OK dispersion"}};

    for (auto const& histo : commonHistos) {
      AxisSpec axis1 = (TString(histo[0]).Contains("hCluE")) ? axisE : axisEp;
      AxisSpec axis2 = (TString(histo[0]).Contains("_v_E_")) ? axisE : axisPt;
      const char* subdir = (TString(histo[0]).Contains("hCluE")) ? "clusterSpectra" : "energyMomentumRatio";
      // Histograms for TPC/TOF identified particles
      mHistManager.add(Form("PID/%s/%s_TPCel", subdir, histo[0]), Form("%s | TPCel", histo[1]), HistType::kTH3F, {axis1, axis2, axisModes});
    }

    if (isMC) {
      mHistManager.add("TrueEl/hTrueElInPhos", "True electrons in PHOS acceptance", HistType::kTH2F, {axisPt, axisModes});
      mHistManager.add("TrueEl/hTrueElWithCluster", "True electrons with a cluster", HistType::kTH2F, {axisPt, axisModes});

      for (auto const& histo : commonHistos) {
        AxisSpec axis1 = (TString(histo[0]).Contains("hCluE")) ? axisE : axisEp;
        AxisSpec axis2 = (TString(histo[0]).Contains("_v_E_")) ? axisE : axisPt;
        const char* subdir = (TString(histo[0]).Contains("hCluE")) ? "clusterSpectra" : "energyMomentumRatio";

        // Histograms for all true electrons
        mHistManager.add(Form("TrueEl/%s/%s", subdir, histo[0]), Form("%s | TrueEl", histo[1]), HistType::kTH3F, {axis1, axis2, axisModes});

        // Histograms for true electrons that also pass TPC/TOF PID
        mHistManager.add(Form("TrueEl_after_PID/%s/%s_TPCel", subdir, histo[0]), Form("%s | TrueEl+TPCel", histo[1]), HistType::kTH3F, {axis1, axis2, axisModes});
      }
    }

    geomPHOS = std::make_unique<o2::phos::Geometry>("PHOS");

    fSigma_dz = new TF1("fSigma_dz", "[0]/(x+[1])^[2]+pol1(3)", 0.3, 10);
    fSigma_dz->SetParameters(((std::vector<float>)pSigmadZ).at(0), ((std::vector<float>)pSigmadZ).at(1), ((std::vector<float>)pSigmadZ).at(2), ((std::vector<float>)pSigmadZ).at(3), ((std::vector<float>)pSigmadZ).at(4));

    fSigma_dx = new TF1("fSigma_dx", "[0]/x^[1]+[2]", 0.1, 10);
    fSigma_dx->SetParameters(((std::vector<float>)pSigmadX).at(0), ((std::vector<float>)pSigmadX).at(1), ((std::vector<float>)pSigmadX).at(2));

    for (int i = 0; i < mAmountOfModules; i++) {
      fMeandXPosMod[i] = new TF1(Form("funcMeandx_pos_mod%i", i + 1), "[0]/(x+[1])^[2]", 0.1, 10);
      fMeandXPosMod[i]->SetParameters(pMeandXPosMod->at(3 * i), pMeandXPosMod->at(3 * i + 1), pMeandXPosMod->at(3 * i + 2));

      fMeandXNegMod[i] = new TF1(Form("funcMeandx_neg_mod%i", i + 1), "[0]/(x+[1])^[2]", 0.1, 10);
      fMeandXNegMod[i]->SetParameters(pMeandXNegMod->at(3 * i), pMeandXNegMod->at(3 * i + 1), pMeandXNegMod->at(3 * i + 2));
    }
  }

  void processData(SelCollisions::iterator const& collision,
                   aod::CaloClusters const& clusters,
                   soa::Filtered<MyTracks> const& tracks,
                   aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
      o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo == nullptr) {
        LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      bz = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      runNumber = bc.runNumber();
    }
    if (std::fabs(collision.posZ()) > mColMaxZ)
      return;
    mHistManager.fill(HIST("eventCounter"), 0.5);
    if (!isMC && !collision.alias_bit(mEvSelTrig))
      return;
    mHistManager.fill(HIST("eventCounter"), 1.5);
    if (isSel8) {
      if (!collision.sel8())
        return;
      mHistManager.fill(HIST("eventCounter"), 2.5);
    }
    if (clusters.size() == 0)
      return; // Nothing to process
    mHistManager.fill(HIST("collision/hColVX"), collision.posX());
    mHistManager.fill(HIST("collision/hColVY"), collision.posY());
    mHistManager.fill(HIST("collision/hColVZ"), collision.posZ());

    for (auto const& track : tracks) {

      if (!track.has_collision() || !track.hasTPC())
        continue;
      if (track.itsNCls() < ITSnclsMin || track.itsNCls() > ITSnclsMax || !((track.itsClusterMap() & uint8_t(1)) > 0))
        continue;
      if (track.tpcNClsFound() < TPCnclsMin || track.tpcNClsFound() > TPCnclsMax)
        continue;
      if (track.tpcNClsCrossedRows() < TPCnclsCRMin || track.tpcNClsCrossedRows() > TPCnclsCRMax)
        continue;

      mHistManager.fill(HIST("tracks/hTrackVX"), track.x());
      mHistManager.fill(HIST("tracks/hTrackVY"), track.y());
      mHistManager.fill(HIST("tracks/hTrackVZ"), track.z());

      int16_t module;
      float trackX = 999., trackZ = 999.;

      auto trackPar = getTrackPar(track);
      if (!impactOnPHOS(trackPar, track.trackEtaEmcal(), track.trackPhiEmcal(), track.collision_as<SelCollisions>().posZ(), module, trackX, trackZ))
        continue;

      float trackMom = track.p();
      float trackPT = track.pt();

      bool isElectron = false;
      if (track.hasTPC()) {
        float nsigmaTPCEl = track.tpcNSigmaEl();
        float nsigmaTOFEl = track.tofNSigmaEl();
        bool isTPCElectron = nsigmaTPCEl > TPCNSigmaElMin && nsigmaTPCEl < TPCNSigmaElMax;
        bool isTOFElectron = nsigmaTOFEl > TOFNSigmaElMin && nsigmaTOFEl < TOFNSigmaElMax;
        isElectron = isTPCElectron || isTOFElectron;

        float nsigmaTPCPi = track.tpcNSigmaPi();
        float nsigmaTPCKa = track.tpcNSigmaKa();
        float nsigmaTPCPr = track.tpcNSigmaPr();
        bool isPion = nsigmaTPCPi > TPCNSigmaPiMin && nsigmaTPCPi < TPCNSigmaPiMax;
        bool isKaon = nsigmaTPCKa > TPCNSigmaKaMin && nsigmaTPCKa < TPCNSigmaKaMax;
        bool isProton = nsigmaTPCPr > TPCNSigmaPrMin && nsigmaTPCPr < TPCNSigmaPrMax;
        if (isElectron && !(isPion || isKaon || isProton))
          isElectron = true;
      }

      bool posTrack = track.sign() * bz > 0;
      for (auto const& clu : clusters) {
        if (module != clu.mod())
          continue;
        double cluE = clu.e();

        if (cluE < mMinCluE ||
            clu.ncell() < mMinCluNcell ||
            clu.time() > mMaxCluTime || clu.time() < mMinCluTime)
          continue;

        bool isDispOK = false;
        if (mSwapM20M02ForTestLambda)
          isDispOK = testLambda(cluE, clu.m02(), clu.m20(), mShowerShapeCutValue, mUseNegativeCrossTerm);
        else
          isDispOK = testLambda(cluE, clu.m20(), clu.m02(), mShowerShapeCutValue, mUseNegativeCrossTerm);
        float posX = clu.x(), posZ = clu.z(), dX = trackX - posX, dZ = trackZ - posZ, Ep = cluE / trackMom;

        mHistManager.fill(HIST("coordinateMatching/hdZpmod"), dZ, trackPT, module);
        mHistManager.fill(HIST("coordinateMatching/hdXpmod"), dX, trackPT, module);
        if (posTrack) {
          mHistManager.fill(HIST("coordinateMatching/hdZpmod_pos"), dZ, trackPT, module);
          mHistManager.fill(HIST("coordinateMatching/hdXpmod_pos"), dX, trackPT, module);
        } else {
          mHistManager.fill(HIST("coordinateMatching/hdZpmod_neg"), dZ, trackPT, module);
          mHistManager.fill(HIST("coordinateMatching/hdXpmod_neg"), dX, trackPT, module);
        }

        if (isDispOK) {
          mHistManager.fill(HIST("clusterSpectra/hCluE_v_pt_disp"), cluE, trackPT, module);
          mHistManager.fill(HIST("energyMomentumRatio/hEp_v_pt_disp"), Ep, trackPT, module);
          mHistManager.fill(HIST("energyMomentumRatio/hEp_v_E_disp"), Ep, cluE, module);
          if (isElectron) {
            mHistManager.fill(HIST("PID/clusterSpectra/hCluE_v_pt_disp_TPCel"), cluE, trackPT, module);
            mHistManager.fill(HIST("PID/energyMomentumRatio/hEp_v_pt_disp_TPCel"), Ep, trackPT, module);
            mHistManager.fill(HIST("PID/energyMomentumRatio/hEp_v_E_disp_TPCel"), Ep, cluE, module);
          }
        }
        if (clu.trackdist() < NsigmaTrackMatch) {
          mHistManager.fill(HIST("doubleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma"), cluE, trackPT, module);
          mHistManager.fill(HIST("doubleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma"), Ep, trackPT, module);
          mHistManager.fill(HIST("doubleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma"), Ep, cluE, module);
          if (isElectron) {
            mHistManager.fill(HIST("doubleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma_TPCel"), cluE, trackPT, module);
            mHistManager.fill(HIST("doubleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma_TPCel"), Ep, trackPT, module);
            mHistManager.fill(HIST("doubleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma_TPCel"), Ep, cluE, module);
          }
          if (isDispOK) {
            mHistManager.fill(HIST("doubleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma_disp"), cluE, trackPT, module);
            mHistManager.fill(HIST("doubleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma_disp"), Ep, trackPT, module);
            mHistManager.fill(HIST("doubleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma_disp"), Ep, cluE, module);
            if (isElectron) {
              mHistManager.fill(HIST("doubleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma_disp_TPCel"), cluE, trackPT, module);
              mHistManager.fill(HIST("doubleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma_disp_TPCel"), Ep, trackPT, module);
              mHistManager.fill(HIST("doubleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma_disp_TPCel"), Ep, cluE, module);
            }
          }
        }
        if (!isWithinNSigma(module, trackPT, dZ, dX, posTrack))
          continue;
        mHistManager.fill(HIST("clusterSpectra/hCluE_v_pt_Nsigma"), cluE, trackPT, module);
        mHistManager.fill(HIST("energyMomentumRatio/hEp_v_pt_Nsigma"), Ep, trackPT, module);
        mHistManager.fill(HIST("energyMomentumRatio/hEp_v_E_Nsigma"), Ep, cluE, module);
        if (isElectron) {
          mHistManager.fill(HIST("PID/clusterSpectra/hCluE_v_pt_Nsigma_TPCel"), cluE, trackPT, module);
          mHistManager.fill(HIST("PID/energyMomentumRatio/hEp_v_pt_Nsigma_TPCel"), Ep, trackPT, module);
          mHistManager.fill(HIST("PID/energyMomentumRatio/hEp_v_E_Nsigma_TPCel"), Ep, cluE, module);
        }
        if (isDispOK) {
          mHistManager.fill(HIST("clusterSpectra/hCluE_v_pt_Nsigma_disp"), cluE, trackPT, module);
          mHistManager.fill(HIST("energyMomentumRatio/hEp_v_pt_Nsigma_disp"), Ep, trackPT, module);
          mHistManager.fill(HIST("energyMomentumRatio/hEp_v_E_Nsigma_disp"), Ep, cluE, module);
          if (isElectron) {
            mHistManager.fill(HIST("PID/clusterSpectra/hCluE_v_pt_Nsigma_disp_TPCel"), cluE, trackPT, module);
            mHistManager.fill(HIST("PID/energyMomentumRatio/hEp_v_pt_Nsigma_disp_TPCel"), Ep, trackPT, module);
            mHistManager.fill(HIST("PID/energyMomentumRatio/hEp_v_E_Nsigma_disp_TPCel"), Ep, cluE, module);
          }
          phosMatch(collision.index(), clu.index(), track.index());
        }
      }

      mHistManager.fill(HIST("tracks/hTrackPtEtaPhi"), track.pt(), track.eta(), track.phi() * TMath::RadToDeg());
      mHistManager.fill(HIST("tracks/hTrackPtEtaPhi_Phos"), track.pt(), track.trackEtaEmcal(), track.trackPhiEmcal() * TMath::RadToDeg());
      mHistManager.fill(HIST("tracks/hTrackDCA"), track.dcaXY(), track.dcaZ());
      mHistManager.fill(HIST("tracks/hTrackPhosProjMod"), trackX, trackZ, module);
    } // end of double loop

    for (auto const& clu : clusters) {
      double cluE = clu.e(), cluTime = clu.time();
      int mod = clu.mod();
      bool isDispOK = false;
      if (mSwapM20M02ForTestLambda)
        isDispOK = testLambda(cluE, clu.m02(), clu.m20(), mShowerShapeCutValue, mUseNegativeCrossTerm);
      else
        isDispOK = testLambda(cluE, clu.m20(), clu.m02(), mShowerShapeCutValue, mUseNegativeCrossTerm);
      if (cluE > mMinCluE) {
        mHistManager.fill(HIST("clusterSpectra/hCluE_mod_energy_cut"), cluE, mod);
        mHistManager.fill(HIST("clusterSpectra/hCluE_v_mod_v_time"), cluE, cluTime * 1e9, mod);
        if (cluTime < mMaxCluTime && cluTime > mMinCluTime) {
          mHistManager.fill(HIST("clusterSpectra/hCluE_mod_time_cut"), cluE, mod);
          if (clu.ncell() >= mMinCluNcell) {
            mHistManager.fill(HIST("clusterSpectra/hCluE_mod_cell_cut"), cluE, mod);
            mHistManager.fill(HIST("coordinateMatching/hCluXZ_mod"), clu.x(), clu.z(), mod);
            mHistManager.fill(HIST("clusterSpectra/hCluE_ncells_mod"), cluE, clu.ncell(), mod);
            if (isDispOK)
              mHistManager.fill(HIST("clusterSpectra/hCluE_mod_disp"), cluE, mod);
          }
        }
      }

      if (cluE < mMinCluE ||
          clu.ncell() < mMinCluNcell ||
          clu.time() > mMaxCluTime || clu.time() < mMinCluTime)
        continue;

      // The following block is disabled by default as it appears to be broken as of now.
      if (mFillSingleLoopHistos) {
        if (clu.trackdist() > NsigmaTrackMatch)
          continue;
        auto matchedTrack = tracks.iteratorAt(clu.trackIndex());
        if (!matchedTrack.has_collision() || !matchedTrack.hasTPC())
          continue;

        if (matchedTrack.itsNCls() < ITSnclsMin || matchedTrack.itsNCls() > ITSnclsMax || !((matchedTrack.itsClusterMap() & uint8_t(1)) > 0))
          continue;
        if (matchedTrack.tpcNClsFound() < TPCnclsMin || matchedTrack.tpcNClsFound() > TPCnclsMax)
          continue;
        if (matchedTrack.tpcNClsCrossedRows() < TPCnclsCRMin || matchedTrack.tpcNClsCrossedRows() > TPCnclsCRMax)
          continue;

        mHistManager.fill(HIST("singleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma"), cluE, matchedTrack.pt(), mod);
        mHistManager.fill(HIST("singleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma"), cluE / matchedTrack.p(), matchedTrack.pt(), mod);
        mHistManager.fill(HIST("singleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma"), cluE / matchedTrack.p(), cluE, mod);
        if (isDispOK) {
          mHistManager.fill(HIST("singleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma_disp"), cluE, matchedTrack.pt(), mod);
          mHistManager.fill(HIST("singleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma_disp"), cluE / matchedTrack.p(), matchedTrack.pt(), mod);
          mHistManager.fill(HIST("singleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma_disp"), cluE / matchedTrack.p(), cluE, mod);
        }
        bool isElectron = false;
        if (matchedTrack.hasTPC()) {
          float nsigmaTPCEl = matchedTrack.tpcNSigmaEl();
          float nsigmaTOFEl = matchedTrack.tofNSigmaEl();
          bool isTPCElectron = nsigmaTPCEl > TPCNSigmaElMin && nsigmaTPCEl < TPCNSigmaElMax;
          bool isTOFElectron = nsigmaTOFEl > TOFNSigmaElMin && nsigmaTOFEl < TOFNSigmaElMax;
          isElectron = isTPCElectron || isTOFElectron;

          float nsigmaTPCPi = matchedTrack.tpcNSigmaPi();
          float nsigmaTPCKa = matchedTrack.tpcNSigmaKa();
          float nsigmaTPCPr = matchedTrack.tpcNSigmaPr();
          bool isPion = nsigmaTPCPi > TPCNSigmaPiMin && nsigmaTPCPi < TPCNSigmaPiMax;
          bool isKaon = nsigmaTPCKa > TPCNSigmaKaMin && nsigmaTPCKa < TPCNSigmaKaMax;
          bool isProton = nsigmaTPCPr > TPCNSigmaPrMin && nsigmaTPCPr < TPCNSigmaPrMax;
          if (isElectron && !(isPion || isKaon || isProton))
            isElectron = true;
        }
        if (isElectron) {
          mHistManager.fill(HIST("singleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma_TPCel"), cluE, matchedTrack.pt(), mod);
          mHistManager.fill(HIST("singleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma_TPCel"), cluE / matchedTrack.p(), matchedTrack.pt(), mod);
          mHistManager.fill(HIST("singleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma_TPCel"), cluE / matchedTrack.p(), cluE, mod);
          if (isDispOK) {
            mHistManager.fill(HIST("singleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma_disp_TPCel"), cluE, matchedTrack.pt(), mod);
            mHistManager.fill(HIST("singleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma_disp_TPCel"), cluE / matchedTrack.p(), matchedTrack.pt(), mod);
            mHistManager.fill(HIST("singleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma_disp_TPCel"), cluE / matchedTrack.p(), cluE, mod);
          }
        }
      }
    } // end of cluster loop
  }
  PROCESS_SWITCH(PhosElId, processData, "process data", false);

  void processMC(SelCollisions::iterator const& collision,
                 aod::CaloClusters const& clusters,
                 soa::Join<soa::Filtered<MyTracks>, aod::McTrackLabels> const& tracks,
                 aod::McParticles const& mcParticles,
                 aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();

    if (runNumber != bc.runNumber()) {
      LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
      o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo == nullptr) {
        LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      bz = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      runNumber = bc.runNumber();
    }
    if (std::fabs(collision.posZ()) > mColMaxZ)
      return;
    mHistManager.fill(HIST("eventCounter"), 0.5);
    if (!isMC && !collision.alias_bit(mEvSelTrig))
      return;
    mHistManager.fill(HIST("eventCounter"), 1.5);
    if (isSel8) {
      if (!collision.sel8())
        return;
      mHistManager.fill(HIST("eventCounter"), 2.5);
    }
    if (clusters.size() == 0)
      return; // Nothing to process
    mHistManager.fill(HIST("collision/hColVX"), collision.posX());
    mHistManager.fill(HIST("collision/hColVY"), collision.posY());
    mHistManager.fill(HIST("collision/hColVZ"), collision.posZ());

    for (auto const& track : tracks) {

      if (!track.has_collision() || !track.hasTPC())
        continue;
      if (track.itsNCls() < ITSnclsMin || track.itsNCls() > ITSnclsMax || !((track.itsClusterMap() & uint8_t(1)) > 0))
        continue;
      if (track.tpcNClsFound() < TPCnclsMin || track.tpcNClsFound() > TPCnclsMax)
        continue;
      if (track.tpcNClsCrossedRows() < TPCnclsCRMin || track.tpcNClsCrossedRows() > TPCnclsCRMax)
        continue;

      mHistManager.fill(HIST("tracks/hTrackVX"), track.x());
      mHistManager.fill(HIST("tracks/hTrackVY"), track.y());
      mHistManager.fill(HIST("tracks/hTrackVZ"), track.z());

      int16_t module;
      float trackX = 999., trackZ = 999.;

      auto trackPar = getTrackPar(track);
      if (!impactOnPHOS(trackPar, track.trackEtaEmcal(), track.trackPhiEmcal(), track.collision_as<SelCollisions>().posZ(), module, trackX, trackZ))
        continue;

      float trackMom = track.p();
      float trackPT = track.pt();

      bool isElectron = false;
      if (track.hasTPC()) {
        float nsigmaTPCEl = track.tpcNSigmaEl();
        float nsigmaTOFEl = track.tofNSigmaEl();
        bool isTPCElectron = nsigmaTPCEl > TPCNSigmaElMin && nsigmaTPCEl < TPCNSigmaElMax;
        bool isTOFElectron = nsigmaTOFEl > TOFNSigmaElMin && nsigmaTOFEl < TOFNSigmaElMax;
        isElectron = isTPCElectron || isTOFElectron;

        float nsigmaTPCPi = track.tpcNSigmaPi();
        float nsigmaTPCKa = track.tpcNSigmaKa();
        float nsigmaTPCPr = track.tpcNSigmaPr();
        bool isPion = nsigmaTPCPi > TPCNSigmaPiMin && nsigmaTPCPi < TPCNSigmaPiMax;
        bool isKaon = nsigmaTPCKa > TPCNSigmaKaMin && nsigmaTPCKa < TPCNSigmaKaMax;
        bool isProton = nsigmaTPCPr > TPCNSigmaPrMin && nsigmaTPCPr < TPCNSigmaPrMax;
        if (isElectron && !(isPion || isKaon || isProton))
          isElectron = true;
      }

      bool isTrueElectron = false;
      auto mcLabel = track.mcParticleId();
      if (mcLabel > -1 && mcLabel < mcParticles.size()) {
        auto mcpart = mcParticles.iteratorAt(mcLabel);
        if (std::abs(mcpart.pdgCode()) == PDG_t::kElectron) {
          isTrueElectron = true;
        }
      }

      if (isTrueElectron) {
        mHistManager.fill(HIST("TrueEl/hTrueElInPhos"), trackPT, module);
      }

      bool posTrack = track.sign() * bz > 0;
      for (auto const& clu : clusters) {
        if (module != clu.mod())
          continue;
        double cluE = clu.e();

        if (cluE < mMinCluE ||
            clu.ncell() < mMinCluNcell ||
            clu.time() > mMaxCluTime || clu.time() < mMinCluTime)
          continue;

        if (isTrueElectron) {
          mHistManager.fill(HIST("TrueEl/hTrueElWithCluster"), trackPT, module);
        }

        bool isDispOK = false;
        if (mSwapM20M02ForTestLambda)
          isDispOK = testLambda(cluE, clu.m02(), clu.m20(), mShowerShapeCutValue, mUseNegativeCrossTerm);
        else
          isDispOK = testLambda(cluE, clu.m20(), clu.m02(), mShowerShapeCutValue, mUseNegativeCrossTerm);
        float posX = clu.x(), posZ = clu.z(), dX = trackX - posX, dZ = trackZ - posZ, Ep = cluE / trackMom;

        mHistManager.fill(HIST("coordinateMatching/hdZpmod"), dZ, trackPT, module);
        mHistManager.fill(HIST("coordinateMatching/hdXpmod"), dX, trackPT, module);
        if (posTrack) {
          mHistManager.fill(HIST("coordinateMatching/hdZpmod_pos"), dZ, trackPT, module);
          mHistManager.fill(HIST("coordinateMatching/hdXpmod_pos"), dX, trackPT, module);
        } else {
          mHistManager.fill(HIST("coordinateMatching/hdZpmod_neg"), dZ, trackPT, module);
          mHistManager.fill(HIST("coordinateMatching/hdXpmod_neg"), dX, trackPT, module);
        }

        if (isDispOK) {
          mHistManager.fill(HIST("clusterSpectra/hCluE_v_pt_disp"), cluE, trackPT, module);
          mHistManager.fill(HIST("energyMomentumRatio/hEp_v_pt_disp"), Ep, trackPT, module);
          mHistManager.fill(HIST("energyMomentumRatio/hEp_v_E_disp"), Ep, cluE, module);
          if (isElectron) {
            mHistManager.fill(HIST("PID/clusterSpectra/hCluE_v_pt_disp_TPCel"), cluE, trackPT, module);
            mHistManager.fill(HIST("PID/energyMomentumRatio/hEp_v_pt_disp_TPCel"), Ep, trackPT, module);
            mHistManager.fill(HIST("PID/energyMomentumRatio/hEp_v_E_disp_TPCel"), Ep, cluE, module);
            if (isTrueElectron) {
              mHistManager.fill(HIST("TrueEl_after_PID/clusterSpectra/hCluE_v_pt_disp_TPCel"), cluE, trackPT, module);
              mHistManager.fill(HIST("TrueEl_after_PID/energyMomentumRatio/hEp_v_pt_disp_TPCel"), Ep, trackPT, module);
              mHistManager.fill(HIST("TrueEl_after_PID/energyMomentumRatio/hEp_v_E_disp_TPCel"), Ep, cluE, module);
            }
          }
          if (isTrueElectron) {
            mHistManager.fill(HIST("TrueEl/clusterSpectra/hCluE_v_pt_disp"), cluE, trackPT, module);
            mHistManager.fill(HIST("TrueEl/energyMomentumRatio/hEp_v_pt_disp"), Ep, trackPT, module);
            mHistManager.fill(HIST("TrueEl/energyMomentumRatio/hEp_v_E_disp"), Ep, cluE, module);
          }
        }
        if (clu.trackdist() < NsigmaTrackMatch) {
          mHistManager.fill(HIST("doubleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma"), cluE, trackPT, module);
          mHistManager.fill(HIST("doubleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma"), Ep, trackPT, module);
          mHistManager.fill(HIST("doubleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma"), Ep, cluE, module);
          if (isElectron) {
            mHistManager.fill(HIST("doubleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma_TPCel"), cluE, trackPT, module);
            mHistManager.fill(HIST("doubleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma_TPCel"), Ep, trackPT, module);
            mHistManager.fill(HIST("doubleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma_TPCel"), Ep, cluE, module);
          }
          if (isDispOK) {
            mHistManager.fill(HIST("doubleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma_disp"), cluE, trackPT, module);
            mHistManager.fill(HIST("doubleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma_disp"), Ep, trackPT, module);
            mHistManager.fill(HIST("doubleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma_disp"), Ep, cluE, module);
            if (isElectron) {
              mHistManager.fill(HIST("doubleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma_disp_TPCel"), cluE, trackPT, module);
              mHistManager.fill(HIST("doubleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma_disp_TPCel"), Ep, trackPT, module);
              mHistManager.fill(HIST("doubleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma_disp_TPCel"), Ep, cluE, module);
            }
          }
        }
        if (!isWithinNSigma(module, trackPT, dZ, dX, posTrack))
          continue;
        mHistManager.fill(HIST("clusterSpectra/hCluE_v_pt_Nsigma"), cluE, trackPT, module);
        mHistManager.fill(HIST("energyMomentumRatio/hEp_v_pt_Nsigma"), Ep, trackPT, module);
        mHistManager.fill(HIST("energyMomentumRatio/hEp_v_E_Nsigma"), Ep, cluE, module);
        if (isElectron) {
          mHistManager.fill(HIST("PID/clusterSpectra/hCluE_v_pt_Nsigma_TPCel"), cluE, trackPT, module);
          mHistManager.fill(HIST("PID/energyMomentumRatio/hEp_v_pt_Nsigma_TPCel"), Ep, trackPT, module);
          mHistManager.fill(HIST("PID/energyMomentumRatio/hEp_v_E_Nsigma_TPCel"), Ep, cluE, module);
          if (isTrueElectron) {
            mHistManager.fill(HIST("TrueEl_after_PID/clusterSpectra/hCluE_v_pt_Nsigma_TPCel"), cluE, trackPT, module);
            mHistManager.fill(HIST("TrueEl_after_PID/energyMomentumRatio/hEp_v_pt_Nsigma_TPCel"), Ep, trackPT, module);
            mHistManager.fill(HIST("TrueEl_after_PID/energyMomentumRatio/hEp_v_E_Nsigma_TPCel"), Ep, cluE, module);
          }
        }
        if (isTrueElectron) {
          mHistManager.fill(HIST("TrueEl/clusterSpectra/hCluE_v_pt_Nsigma"), cluE, trackPT, module);
          mHistManager.fill(HIST("TrueEl/energyMomentumRatio/hEp_v_pt_Nsigma"), Ep, trackPT, module);
          mHistManager.fill(HIST("TrueEl/energyMomentumRatio/hEp_v_E_Nsigma"), Ep, cluE, module);
        }
        if (isDispOK) {
          mHistManager.fill(HIST("clusterSpectra/hCluE_v_pt_Nsigma_disp"), cluE, trackPT, module);
          mHistManager.fill(HIST("energyMomentumRatio/hEp_v_pt_Nsigma_disp"), Ep, trackPT, module);
          mHistManager.fill(HIST("energyMomentumRatio/hEp_v_E_Nsigma_disp"), Ep, cluE, module);
          if (isElectron) {
            mHistManager.fill(HIST("PID/clusterSpectra/hCluE_v_pt_Nsigma_disp_TPCel"), cluE, trackPT, module);
            mHistManager.fill(HIST("PID/energyMomentumRatio/hEp_v_pt_Nsigma_disp_TPCel"), Ep, trackPT, module);
            mHistManager.fill(HIST("PID/energyMomentumRatio/hEp_v_E_Nsigma_disp_TPCel"), Ep, cluE, module);
            if (isTrueElectron) {
              mHistManager.fill(HIST("TrueEl_after_PID/clusterSpectra/hCluE_v_pt_Nsigma_disp_TPCel"), cluE, trackPT, module);
              mHistManager.fill(HIST("TrueEl_after_PID/energyMomentumRatio/hEp_v_pt_Nsigma_disp_TPCel"), Ep, trackPT, module);
              mHistManager.fill(HIST("TrueEl_after_PID/energyMomentumRatio/hEp_v_E_Nsigma_disp_TPCel"), Ep, cluE, module);
            }
          }
          if (isTrueElectron) {
            mHistManager.fill(HIST("TrueEl/clusterSpectra/hCluE_v_pt_Nsigma_disp"), cluE, trackPT, module);
            mHistManager.fill(HIST("TrueEl/energyMomentumRatio/hEp_v_pt_Nsigma_disp"), Ep, trackPT, module);
            mHistManager.fill(HIST("TrueEl/energyMomentumRatio/hEp_v_E_Nsigma_disp"), Ep, cluE, module);
          }
          phosMatch(collision.index(), clu.index(), track.index());
        }
      }

      mHistManager.fill(HIST("tracks/hTrackPtEtaPhi"), track.pt(), track.eta(), track.phi() * TMath::RadToDeg());
      mHistManager.fill(HIST("tracks/hTrackPtEtaPhi_Phos"), track.pt(), track.trackEtaEmcal(), track.trackPhiEmcal() * TMath::RadToDeg());
      mHistManager.fill(HIST("tracks/hTrackDCA"), track.dcaXY(), track.dcaZ());
      mHistManager.fill(HIST("tracks/hTrackPhosProjMod"), trackX, trackZ, module);
    } // end of double loop

    for (auto const& clu : clusters) {
      double cluE = clu.e(), cluTime = clu.time();
      int mod = clu.mod();
      bool isDispOK = false;
      if (mSwapM20M02ForTestLambda)
        isDispOK = testLambda(cluE, clu.m02(), clu.m20(), mShowerShapeCutValue, mUseNegativeCrossTerm);
      else
        isDispOK = testLambda(cluE, clu.m20(), clu.m02(), mShowerShapeCutValue, mUseNegativeCrossTerm);
      if (cluE > mMinCluE) {
        mHistManager.fill(HIST("clusterSpectra/hCluE_mod_energy_cut"), cluE, mod);
        mHistManager.fill(HIST("clusterSpectra/hCluE_v_mod_v_time"), cluE, cluTime * 1e9, mod);
        if (cluTime < mMaxCluTime && cluTime > mMinCluTime) {
          mHistManager.fill(HIST("clusterSpectra/hCluE_mod_time_cut"), cluE, mod);
          if (clu.ncell() >= mMinCluNcell) {
            mHistManager.fill(HIST("clusterSpectra/hCluE_mod_cell_cut"), cluE, mod);
            mHistManager.fill(HIST("coordinateMatching/hCluXZ_mod"), clu.x(), clu.z(), mod);
            mHistManager.fill(HIST("clusterSpectra/hCluE_ncells_mod"), cluE, clu.ncell(), mod);
            if (isDispOK)
              mHistManager.fill(HIST("clusterSpectra/hCluE_mod_disp"), cluE, mod);
          }
        }
      }

      if (cluE < mMinCluE ||
          clu.ncell() < mMinCluNcell ||
          clu.time() > mMaxCluTime || clu.time() < mMinCluTime)
        continue;

      // The following block is disabled by default because the track and cluster labels are not consistent, causing crashes when trying to access tracks using cluster information.
      if (mFillSingleLoopHistos) {
        if (clu.trackdist() > NsigmaTrackMatch)
          continue;
        if (clu.trackIndex() == UCHAR_MAX) {
          continue;
        }
        auto matchedTrack = tracks.iteratorAt(clu.trackIndex());

        if (!matchedTrack.has_collision() || !matchedTrack.hasTPC())
          continue;

        if (matchedTrack.itsNCls() < ITSnclsMin || matchedTrack.itsNCls() > ITSnclsMax || !((matchedTrack.itsClusterMap() & uint8_t(1)) > 0))
          continue;
        if (matchedTrack.tpcNClsFound() < TPCnclsMin || matchedTrack.tpcNClsFound() > TPCnclsMax)
          continue;
        if (matchedTrack.tpcNClsCrossedRows() < TPCnclsCRMin || matchedTrack.tpcNClsCrossedRows() > TPCnclsCRMax)
          continue;

        mHistManager.fill(HIST("singleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma"), cluE, matchedTrack.pt(), mod);
        mHistManager.fill(HIST("singleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma"), cluE / matchedTrack.p(), matchedTrack.pt(), mod);
        mHistManager.fill(HIST("singleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma"), cluE / matchedTrack.p(), cluE, mod);
        if (isDispOK) {
          mHistManager.fill(HIST("singleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma_disp"), cluE, matchedTrack.pt(), mod);
          mHistManager.fill(HIST("singleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma_disp"), cluE / matchedTrack.p(), matchedTrack.pt(), mod);
          mHistManager.fill(HIST("singleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma_disp"), cluE / matchedTrack.p(), cluE, mod);
        }
        bool isElectron = false;
        if (matchedTrack.hasTPC()) {
          float nsigmaTPCEl = matchedTrack.tpcNSigmaEl();
          float nsigmaTOFEl = matchedTrack.tofNSigmaEl();
          bool isTPCElectron = nsigmaTPCEl > TPCNSigmaElMin && nsigmaTPCEl < TPCNSigmaElMax;
          bool isTOFElectron = nsigmaTOFEl > TOFNSigmaElMin && nsigmaTOFEl < TOFNSigmaElMax;
          isElectron = isTPCElectron || isTOFElectron;

          float nsigmaTPCPi = matchedTrack.tpcNSigmaPi();
          float nsigmaTPCKa = matchedTrack.tpcNSigmaKa();
          float nsigmaTPCPr = matchedTrack.tpcNSigmaPr();
          bool isPion = nsigmaTPCPi > TPCNSigmaPiMin && nsigmaTPCPi < TPCNSigmaPiMax;
          bool isKaon = nsigmaTPCKa > TPCNSigmaKaMin && nsigmaTPCKa < TPCNSigmaKaMax;
          bool isProton = nsigmaTPCPr > TPCNSigmaPrMin && nsigmaTPCPr < TPCNSigmaPrMax;
          if (isElectron && !(isPion || isKaon || isProton))
            isElectron = true;
        }
        if (isElectron) {
          mHistManager.fill(HIST("singleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma_TPCel"), cluE, matchedTrack.pt(), mod);
          mHistManager.fill(HIST("singleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma_TPCel"), cluE / matchedTrack.p(), matchedTrack.pt(), mod);
          mHistManager.fill(HIST("singleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma_TPCel"), cluE / matchedTrack.p(), cluE, mod);
          if (isDispOK) {
            mHistManager.fill(HIST("singleLoop/trackdist/clusterSpectra/hCluE_v_pt_Nsigma_disp_TPCel"), cluE, matchedTrack.pt(), mod);
            mHistManager.fill(HIST("singleLoop/trackdist/energyMomentumRatio/hEp_v_pt_Nsigma_disp_TPCel"), cluE / matchedTrack.p(), matchedTrack.pt(), mod);
            mHistManager.fill(HIST("singleLoop/trackdist/energyMomentumRatio/hEp_v_E_Nsigma_disp_TPCel"), cluE / matchedTrack.p(), cluE, mod);
          }
        }
      }
    } // end of cluster loop
  }
  PROCESS_SWITCH(PhosElId, processMC, "process mc", false);

  void processDummy(SelCollisions::iterator const&) {}
  PROCESS_SWITCH(PhosElId, processDummy, "Dummy process", true);

  bool isWithinNSigma(int16_t& mod, float p, float deltaZ, float deltaX, bool positiveCharge)
  {
    int modMinus1 = mod - 1;
    if (std::fabs(deltaZ - ((std::vector<float>)pPhosShiftZ).at(modMinus1)) > NsigmaTrackMatch * fSigma_dz->Eval(p))
      return false;
    if (positiveCharge) {
      if (std::fabs(deltaX - fMeandXPosMod[modMinus1]->Eval(p) - ((std::vector<float>)pPhosShiftX).at(modMinus1)) > NsigmaTrackMatch * fSigma_dx->Eval(p))
        return false;
    } else {
      if (std::fabs(deltaX - fMeandXNegMod[modMinus1]->Eval(p) - ((std::vector<float>)pPhosShiftX).at(modMinus1)) > NsigmaTrackMatch * fSigma_dx->Eval(p))
        return false;
    }
    return true;
  }
  /////////////////////////////////////////////////////////////////////// taken from PHOSAlign
  bool impactOnPHOS(o2::track::TrackParametrization<float>& trackPar, float trackEta, float trackPhi, float /*zvtx*/, int16_t& module, float& trackX, float& trackZ)
  {
    // eta,phi was calculated at EMCAL radius.
    // Extrapolate to PHOS assuming zeroB and current vertex
    // Check if direction in PHOS acceptance+20cm and return phos module number and coordinates in PHOS module plane
    const float phiMin = 240. * 0.017453293; // degToRad
    const float phiMax = 323. * 0.017453293; // PHOS+20 cm * degToRad
    const float etaMax = 0.178266;
    if (trackPhi < phiMin || trackPhi > phiMax || std::abs(trackEta) > etaMax) {
      return false;
    }
    const float dphi = 20. * 0.017453293;
    trackPhi = RecoDecay::constrainAngle<float, float>(trackPhi);

    module = 1 + static_cast<int16_t>((trackPhi - phiMin) / dphi);
    if (module < 1) {
      module = 1;
    }
    if (module > mAmountOfModules) { // > 4
      module = mAmountOfModules;     // = 4
    }

    // get PHOS radius
    constexpr float ShiftY = -1.26;    // Depth-optimized
    double posL[3] = {0., 0., ShiftY}; // local position at the center of module
    double posG[3] = {0};
    geomPHOS->getAlignmentMatrix(module)->LocalToMaster(posL, posG);
    double rPHOS = std::sqrt(posG[0] * posG[0] + posG[1] * posG[1]);
    double alpha = (230. + 20. * module) * 0.017453293;

    // During main reconstruction track was propagated to radius 460 cm with accounting material
    // now material is not available. Therefore, start from main rec. position and extrapoate to actual radius without material
    // Get track parameters at point where main reconstruction stop
    float xPHOS = 460.f, xtrg = 0.f;
    if (!trackPar.getXatLabR(xPHOS, xtrg, bz, o2::track::DirType::DirOutward)) {
      return false;
    }
    auto prop = o2::base::Propagator::Instance();
    if (!trackPar.rotate(alpha) ||
        !prop->PropagateToXBxByBz(trackPar, xtrg, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrNONE)) {
      return false;
    }
    // calculate xyz from old (Phi, eta) and new r
    float r = std::sqrt(trackPar.getX() * trackPar.getX() + trackPar.getY() * trackPar.getY());
    trackPar.setX(r * std::cos(trackPhi - alpha));
    trackPar.setY(r * std::sin(trackPhi - alpha));
    trackPar.setZ(r / std::tan(2. * std::atan(std::exp(-trackEta))));

    if (!prop->PropagateToXBxByBz(trackPar, rPHOS, 0.95, 10, o2::base::Propagator::MatCorrType::USEMatCorrNONE)) {
      return false;
    }
    alpha = trackPar.getAlpha();
    double ca = std::cos(alpha), sa = std::sin(alpha);
    posG[0] = trackPar.getX() * ca - trackPar.getY() * sa;
    posG[1] = trackPar.getY() * ca + trackPar.getX() * sa;
    posG[2] = trackPar.getZ();

    geomPHOS->getAlignmentMatrix(module)->MasterToLocal(posG, posL);
    trackX = posL[0];
    trackZ = posL[1];
    return true;
  }
};

struct MassSpectra {

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults,
                                  aod::FT0sCorrected, aod::CentFT0Ms,
                                  aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As,
                                  aod::CentFDDMs, aod::CentNTPVs>;
  using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                             aod::pidTOFFullEl, aod::pidTPCFullEl, aod::pidTPCFullPi,
                             aod::pidTPCFullKa, aod::pidTPCFullPr>;
  Configurable<bool> isMC{"isMC", false, "Enable MC analysis"},
    isSel8{"isSel8", 1, "check if event is Single Event Latch-up 8"};
  Configurable<int> mEvSelTrig{"mEvSelTrig", kTVXinPHOS, "Select events with this trigger"},
    MassBinning{"MassBinning", 1000, "Binning for mass"},
    EnergyBinning{"EnergyBinning", 100, "Binning for energy"},
    EpRatioBinning{"EpRatioBinning", 200, "Binning for energy to momentum ratio"},
    CentBinning{"CentBinning", 10, "Binning for centrality"},
    CentEst{"CentEst", 1, "Centrality estimator, 0: FV0A, 1: FT0M, 2: FT0A, 3: FT0C, 4: FDDM, 5: NTPV"};

  Configurable<float> mColMaxZ{"mColMaxZ", 10.f, "maximum z accepted in analysis"},
    fEtaMax{"fEtaMax", {0.8f}, "eta ranges"},
    fEtaMaxPhos{"fEtaMaxPhos", {0.15f}, "eta ranges of phos"},
    fPtMin{"fPtMin", {0.2f}, "pt min"},
    fPtMax{"fPtMax", {20.f}, "pt max"},
    fMassSpectraMin{"fMassSpectraMin", {2.5f}, "mass spectra min for e+e-"},
    fMassSpectraMax{"fMassSpectraMax", {3.5f}, "mass spcetra max for e+e-"},
    fDCAxyMax{"fDCAxyMax", {3.f}, "dcaxy max"},
    fDCAzMax{"fDCAzMax", {3.f}, "dcaz max"},
    fITSchi2Max{"fITSchi2Max", {5.f}, "its chi2 max"},
    fITSnclsMin{"fITSnclsMin", {4.5f}, "min number of ITS clusters"},
    fITSnclsMax{"fITSnclsMax", {7.5f}, "max number of ITS clusters"},
    fTPCchi2Max{"fTPCchi2Max", {4.f}, "tpc chi2 max"},
    fTPCnclsMin{"fTPCnclsMin", {90.f}, "min number of TPC clusters"},
    fTPCnclsMax{"fTPCnclsMax", {170.f}, "max number of TPC clusters"},
    fTPCnclsCRMin{"fTPCnclsCRMin", {80.f}, "min number of TPC crossed rows"},
    fTPCnclsCRMax{"fTPCnclsCRMax", {161.f}, "max number of TPC crossed rows"},
    fTPCNSigmaElMin{"fTPCNSigmaElMin", {-3.f}, "min TPC nsigma e for inclusion"},
    fTPCNSigmaElMax{"fTPCNSigmaElMax", {2.f}, "max TPC nsigma e for inclusion"},
    fTPCNSigmaPiMin{"fTPCNSigmaPiMin", {-3.f}, "min TPC nsigma pion for exclusion"},
    fTPCNSigmaPiMax{"fTPCNSigmaPiMax", {3.5f}, "max TPC nsigma pion for exclusion"},
    fTPCNSigmaPrMin{"fTPCNSigmaPrMin", {-3.f}, "min TPC nsigma proton for exclusion"},
    fTPCNSigmaPrMax{"fTPCNSigmaPrMax", {4.f}, "max TPC nsigma proton for exclusion"},
    fTPCNSigmaKaMin{"fTPCNSigmaKaMin", {-3.f}, "min TPC nsigma kaon for exclusion"},
    fTPCNSigmaKaMax{"fTPCNSigmaKaMax", {4.f}, "max TPC nsigma kaon for exclusion"},
    fTOFNSigmaElMin{"fTOFNSigmaElMin", {-3.f}, "min TOF nsigma e for inclusion"},
    fTOFNSigmaElMax{"fTOFNSigmaElMax", {3.f}, "max TOF nsigma e for inclusion"},
    fShiftEp{"fShiftEp", {0.055f}, "PHOS E/p shift for electrons"},
    fNsigmaEp{"fNsigmaEp", {2.f}, "PHOS E/p nsigma for inclusion"};

  Configurable<std::vector<float>> fEpSigmaPars{"fEpSigmaPars", {1.3e-02, 1.9e-02, 1.1e-02, 3.e-02}, "E/p sigma function parameters (from alice 3 mc tests + const)"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  std::unique_ptr<o2::phos::Geometry> geomPHOS;
  double bz{0.}; // magnetic field
  int runNumber{0};

  HistogramRegistry mHistManager{"MassSpectraHistograms"};
  TF1* fEpSigmaPhos;

  void init(InitContext const&)
  {
    LOG(info) << "Initializing PHOS electron identification analysis task ...";

    std::vector<double> momentumBinning = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                                           1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                                           4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 10.};
    const AxisSpec axisCounter{3, 0, +3, ""},
      axisCent{CentBinning, 0, 100, "centrality percentage"},
      axisPt{momentumBinning, "p_{T} (GeV/c)"},
      axisEp{EpRatioBinning, 0., 2., "E/p", "E_{cluster}/p_{track}"},
      axisE{EnergyBinning, 0, 10, "E (GeV)", "E (GeV)"},
      axisMassSpectrum{MassBinning, fMassSpectraMin, fMassSpectraMax, "M (GeV/c^{2})", "Mass e^{+}e^{-} (GeV/c^{2})"};

    mHistManager.add("eventCounter", "eventCounter", kTH1F, {axisCounter});

    mHistManager.add("h_eh_pp_mass_spectra_v_pt_v_cent", "Mass e^{+}h^{+} vs momentum e^{+}h^{+}", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("h_ee_pp_mass_spectra_v_pt_v_cent", "Mass e^{+}e^{+} vs momentum e^{+}e^{+}", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("h_eh_mm_mass_spectra_v_pt_v_cent", "Mass e^{-}h^{-} vs momentum e^{-}h^{-}", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("h_ee_mm_mass_spectra_v_pt_v_cent", "Mass e^{-}e^{-} vs momentum e^{-}e^{-}", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});

    mHistManager.add("h_eh_pp_mass_spectra_v_E_v_cent", "Mass e^{+}h^{+} vs cluster E e^{+}h^{+}", HistType::kTH3F, {axisMassSpectrum, axisE, axisCent});
    mHistManager.add("h_ee_pp_mass_spectra_v_E_v_cent", "Mass e^{+}e^{+} vs cluster E e^{+}e^{+}", HistType::kTH3F, {axisMassSpectrum, axisE, axisCent});
    mHistManager.add("h_eh_mm_mass_spectra_v_E_v_cent", "Mass e^{-}h^{-} vs cluster E e^{-}h^{-}", HistType::kTH3F, {axisMassSpectrum, axisE, axisCent});
    mHistManager.add("h_ee_mm_mass_spectra_v_E_v_cent", "Mass e^{-}e^{-} vs cluster E e^{-}e^{-}", HistType::kTH3F, {axisMassSpectrum, axisE, axisCent});

    mHistManager.add("h_eh_mp_mass_spectra_v_pt_v_cent", "Mass e^{#pm}h^{#mp} vs momentum e^{#pm}h^{#mp}", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("h_ee_mp_mass_spectra_v_pt_v_cent", "Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("h_eh_mp_mass_spectra_v_E_v_cent", "Mass e^{#pm}h^{#mp} vs cluster E e^{#pm}h^{#mp}", HistType::kTH3F, {axisMassSpectrum, axisE, axisCent});
    mHistManager.add("h_ee_mp_mass_spectra_v_E_v_cent", "Mass e^{#pm}e^{#mp} vs cluster E e^{#pm}e^{#mp}", HistType::kTH3F, {axisMassSpectrum, axisE, axisCent});

    mHistManager.add("hEp_v_E_v_cent", "E/p ratio vs cluster E", HistType::kTH3F, {axisEp, axisE, axisCent});
    mHistManager.add("hEp_v_E_v_cent_cutEp", "E/p ratio vs cluster E within nSigma corridor", HistType::kTH3F, {axisEp, axisE, axisCent});

    geomPHOS = std::make_unique<o2::phos::Geometry>("PHOS");

    std::vector<float> epSigmaPars = fEpSigmaPars;
    fEpSigmaPhos = new TF1("fEpSigmaPhos", "sqrt([0]*[0]/x/x+[1]*[1]/x+[2]*[2])+[3]", 0.01, 10);
    fEpSigmaPhos->SetParameters(epSigmaPars.at(0), epSigmaPars.at(1), epSigmaPars.at(2), epSigmaPars.at(3));
  }

  void process(SelCollisions::iterator const& collision,
               aod::CaloClusters const& clusters,
               MyTracks const& tracks,
               o2::aod::PHOSMatchindexTable const& matches,
               aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
      o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo == nullptr) {
        LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      bz = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      runNumber = bc.runNumber();
    }
    if (std::fabs(collision.posZ()) > mColMaxZ)
      return;
    mHistManager.fill(HIST("eventCounter"), 0.5);
    if (!isMC && !collision.alias_bit(mEvSelTrig))
      return;
    mHistManager.fill(HIST("eventCounter"), 1.5);
    if (isSel8) {
      if (!collision.sel8())
        return;
      mHistManager.fill(HIST("eventCounter"), 2.5);
    }

    if (clusters.size() == 0)
      return; // Nothing to process

    float cent = -1.;
    switch (CentEst) {
      case FV0A:
        cent = collision.centFV0A();
        break;
      case FT0M:
        cent = collision.centFT0M();
        break;
      case FT0A:
        cent = collision.centFT0A();
        break;
      case FT0C:
        cent = collision.centFT0C();
        break;
      case FDDM:
        cent = collision.centFDDM();
        break;
      case NTPV:
        cent = collision.centNTPV();
        break;
    }

    for (auto const& TPCel : tracks) {
      if (!TPCel.has_collision() || std::fabs(TPCel.dcaXY()) > fDCAxyMax || std::fabs(TPCel.dcaZ()) > fDCAzMax || !TPCel.hasTPC() || std::fabs(TPCel.eta()) > fEtaMaxPhos)
        continue;
      if (TPCel.pt() < fPtMin || TPCel.pt() > fPtMax)
        continue;
      if (TPCel.itsChi2NCl() > fITSchi2Max)
        continue;
      if (TPCel.itsNCls() < fITSnclsMin || TPCel.itsNCls() > fITSnclsMax || !((TPCel.itsClusterMap() & uint8_t(1)) > 0))
        continue;
      if (TPCel.tpcChi2NCl() > fTPCchi2Max)
        continue;
      if (TPCel.tpcNClsFound() < fTPCnclsMin || TPCel.tpcNClsFound() > fTPCnclsMax)
        continue;
      if (TPCel.tpcNClsCrossedRows() < fTPCnclsCRMin || TPCel.tpcNClsCrossedRows() > fTPCnclsCRMax)
        continue;

      bool isElectron = false;
      float nsigmaTPCEl = TPCel.tpcNSigmaEl();
      float nsigmaTOFEl = TPCel.tofNSigmaEl();
      bool isTPCElectron = nsigmaTPCEl > fTPCNSigmaElMin && nsigmaTPCEl < fTPCNSigmaElMax;
      bool isTOFElectron = nsigmaTOFEl > fTOFNSigmaElMin && nsigmaTOFEl < fTOFNSigmaElMax;
      isElectron = isTPCElectron || isTOFElectron;

      float nsigmaTPCPi = TPCel.tpcNSigmaPi();
      float nsigmaTPCKa = TPCel.tpcNSigmaKa();
      float nsigmaTPCPr = TPCel.tpcNSigmaPr();
      bool isPion = nsigmaTPCPi > fTPCNSigmaPiMin && nsigmaTPCPi < fTPCNSigmaPiMax;
      bool isKaon = nsigmaTPCKa > fTPCNSigmaKaMin && nsigmaTPCKa < fTPCNSigmaKaMax;
      bool isProton = nsigmaTPCPr > fTPCNSigmaPrMin && nsigmaTPCPr < fTPCNSigmaPrMax;
      if (isElectron && !(isPion || isKaon || isProton))
        isElectron = true;
      if (!isElectron)
        continue;
      bool posTrack = TPCel.sign() * bz > 0;

      for (auto const& match : matches) {
        auto clust2 = clusters.iteratorAt(match.caloClusterId());
        auto track2 = tracks.iteratorAt(match.trackId());

        if (TPCel.collisionId() != track2.collisionId())
          continue;
        if (TPCel.index() >= track2.index())
          break;

        float mass2Tracks = 0, momProbeTrack = track2.pt(), cluE = clust2.e();
        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> fourVectorP1, fourVectorP2;
        fourVectorP1.SetPxPyPzE(TPCel.px(), TPCel.py(), TPCel.pz(), TPCel.energy(0));
        fourVectorP2.SetPxPyPzE(track2.px(), track2.py(), track2.pz(), track2.energy(0));
        mass2Tracks = (fourVectorP1 + fourVectorP2).M();
        bool elCandidate = (std::fabs(cluE / track2.p() - fShiftEp - 1) < fNsigmaEp * fEpSigmaPhos->Eval(cluE));

        if (TPCel.sign() == track2.sign()) {
          if (posTrack) {
            mHistManager.fill(HIST("h_eh_pp_mass_spectra_v_pt_v_cent"), mass2Tracks, momProbeTrack, cent);
            mHistManager.fill(HIST("h_eh_pp_mass_spectra_v_E_v_cent"), mass2Tracks, cluE, cent);
            if (elCandidate) {
              mHistManager.fill(HIST("h_ee_pp_mass_spectra_v_pt_v_cent"), mass2Tracks, momProbeTrack, cent);
              mHistManager.fill(HIST("h_ee_pp_mass_spectra_v_E_v_cent"), mass2Tracks, cluE, cent);
            }
          } else {
            mHistManager.fill(HIST("h_eh_mm_mass_spectra_v_pt_v_cent"), mass2Tracks, momProbeTrack, cent);
            mHistManager.fill(HIST("h_eh_mm_mass_spectra_v_E_v_cent"), mass2Tracks, cluE, cent);
            if (elCandidate) {
              mHistManager.fill(HIST("h_ee_mm_mass_spectra_v_pt_v_cent"), mass2Tracks, momProbeTrack, cent);
              mHistManager.fill(HIST("h_ee_mm_mass_spectra_v_E_v_cent"), mass2Tracks, cluE, cent);
            }
          }
        } else {
          mHistManager.fill(HIST("h_eh_mp_mass_spectra_v_pt_v_cent"), mass2Tracks, momProbeTrack, cent);
          mHistManager.fill(HIST("h_eh_mp_mass_spectra_v_E_v_cent"), mass2Tracks, cluE, cent);
          if (elCandidate) {
            mHistManager.fill(HIST("h_ee_mp_mass_spectra_v_pt_v_cent"), mass2Tracks, momProbeTrack, cent);
            mHistManager.fill(HIST("h_ee_mp_mass_spectra_v_E_v_cent"), mass2Tracks, cluE, cent);
          }
        }
      }
    } // end of double loop

    for (auto const& match : matches) {
      auto clust = clusters.iteratorAt(match.caloClusterId());
      auto track = tracks.iteratorAt(match.trackId());
      float cluE = clust.e();
      float epRatio = cluE / track.p();
      mHistManager.fill(HIST("hEp_v_E_v_cent"), epRatio, cluE, cent);
      bool elCandidate = (std::fabs(epRatio - fShiftEp - 1) < fNsigmaEp * fEpSigmaPhos->Eval(cluE));
      if (elCandidate)
        mHistManager.fill(HIST("hEp_v_E_v_cent_cutEp"), epRatio, cluE, cent);
    }
  }
  PROCESS_SWITCH(MassSpectra, process, "process", false);

  void processDummy(SelCollisions::iterator const&) {}
  PROCESS_SWITCH(MassSpectra, processDummy, "Dummy process", true);
};

struct TpcElIdMassSpectrum {

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults,
                                  aod::FT0sCorrected, aod::CentFT0Ms,
                                  aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As,
                                  aod::CentFDDMs, aod::CentNTPVs>;
  using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                             aod::pidTOFFullEl, aod::pidTPCFullEl, aod::pidTPCFullPi,
                             aod::pidTPCFullKa, aod::pidTPCFullPr>;
  Configurable<bool> isSel8{"isSel8", 1, "check if event is Single Event Latch-up 8"},
    mSwapM20M02ForTestLambda{"mSwapM20M02ForTestLambda", false, "Swap m20 and m02 arguments for testLambda (false for note's correct order, true for swapped/original incorrect order)"},
    mUseNegativeCrossTerm{"mUseNegativeCrossTerm", true, "Use negative sign for the cross-term in testLambda (true for analysis note version, false for old version)"};

  Configurable<bool> isMC{"isMC", true, "Enable MC analysis"};
  Configurable<float> mColMaxZ{"mColMaxZ", 10.f, "maximum z accepted in analysis"},
    mMinCluE{"mMinCluE", 0.1, "Minimum cluster energy for photons in the analysis"},
    mCutMIPCluE{"mCutMIPCluE", 0.3, "Min cluster energy to reject MIPs in the analysis"},
    mMaxCluE{"mMaxCluE", 1., "Maximum cluster energy for photons in the analysis"},
    mMinCluTime{"minCluTime", -25.e-9, "Min. cluster time"},
    mMaxCluTime{"mMaxCluTime", 25.e-9, "Max. cluster time"},
    EtaMax{"EtaMax", {0.8f}, "eta ranges"},
    PtMin{"PtMin", {0.2f}, "pt min"},
    PtMax{"PtMax", {20.f}, "pt max"},
    MassSpectraJpsiMin{"MassSpectraJpsiMin", {0.5f}, "mass spectra min for Jpsi region"},
    MassSpectraJpsiMax{"MassSpectraJpsiMax", {3.5f}, "mass spcetra max for Jpsi region"},
    MassSpectraChicMin{"MassSpectraChicMin", {3.f}, "mass spectra min Chic region"},
    MassSpectraChicMax{"MassSpectraChicMax", {4.f}, "mass spcetra max Chic region"},
    DCAxyMax{"DCAxyMax", {3.f}, "dcaxy max"},
    DCAzMax{"DCAzMax", {3.f}, "dcaz max"},
    ITSchi2Max{"ITSchi2Max", {5.f}, "its chi2 max"},
    ITSnclsMin{"ITSnclsMin", {4.5f}, "min number of ITS clusters"},
    ITSnclsMax{"ITSnclsMax", {7.5f}, "max number of ITS clusters"},
    TPCchi2Max{"TPCchi2Max", {4.f}, "tpc chi2 max"},
    TPCnclsMin{"TPCnclsMin", {90.f}, "min number of TPC clusters"},
    TPCnclsMax{"TPCnclsMax", {170.f}, "max number of TPC clusters"},
    TPCnclsCRMin{"TPCnclsCRMin", {80.f}, "min number of TPC crossed rows"},
    TPCnclsCRMax{"TPCnclsCRMax", {161.f}, "max number of TPC crossed rows"},
    TPCNSigmaElMin{"TPCNSigmaElMin", {-3.f}, "min TPC nsigma e for inclusion"},
    TPCNSigmaElMax{"TPCNSigmaElMax", {2.f}, "max TPC nsigma e for inclusion"},
    TPCNSigmaPiMin{"TPCNSigmaPiMin", {-3.f}, "min TPC nsigma pion for exclusion"},
    TPCNSigmaPiMax{"TPCNSigmaPiMax", {3.5f}, "max TPC nsigma pion for exclusion"},
    TPCNSigmaPrMin{"TPCNSigmaPrMin", {-3.f}, "min TPC nsigma proton for exclusion"},
    TPCNSigmaPrMax{"TPCNSigmaPrMax", {4.f}, "max TPC nsigma proton for exclusion"},
    TPCNSigmaKaMin{"TPCNSigmaKaMin", {-3.f}, "min TPC nsigma kaon for exclusion"},
    TPCNSigmaKaMax{"TPCNSigmaKaMax", {4.f}, "max TPC nsigma kaon for exclusion"},
    TOFNSigmaElMin{"TOFNSigmaElMin", {-3.f}, "min TOF nsigma e for inclusion"},
    TOFNSigmaElMax{"TOFNSigmaElMax", {3.f}, "max TOF nsigma e for inclusion"},
    PhosRangeEta{"PhosRangeEta", {0.12f}, "Phos range definition plus minus eta"},
    PhosRangePhiMin{"PhosRangePhiMin", {230.f}, "Phos range angle phi min"},
    PhosRangePhiMax{"PhosRangePhiMax", {330.f}, "Phos range angle phi max"},
    eeMassMin{"eeMassMin", {2.9f}, "J/psi(e+e-) Mass corridor lower limit (for Chic selection)"},
    eeMassMax{"eeMassMax", {3.3f}, "J/psi(e+e-) Mass corridor upper limit (for Chic selection)"},
    JpsiMass{"JpsiMass", {3.097f}, "J/psi Mass constant"},
    mMassSpectrumLowerCutoff{"mMassSpectrumLowerCutoff", {0.01f}, "Used to exclude 0+0 masses"},
    mShowerShapeCutValue{"mShowerShapeCutValue", 4.f, "Cut threshold for testLambda shower shape"};

  Configurable<int> mEvSelTrig{"mEvSelTrig", kTVXinPHOS, "Select events with this trigger"},
    CentBinning{"CentBinning", 10, "Binning for centrality"},
    CentEst{"CentEst", 1, "Centrality estimator, 0: FV0A, 1: FT0M, 2: FT0A, 3: FT0C, 4: FDDM, 5: NTPV"},
    MassBinning{"MassBinning", 1000, "Binning for mass"},
    EnergyBinning{"EnergyBinning", 100, "Binning for energy"},
    mMinCluNcell{"minCluNcell", 3, "min cells in cluster"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  double bz{0.};
  int runNumber{0};

  HistogramRegistry mHistManager{"tpcElIdHistograms"};

  void init(InitContext const&)
  {
    LOG(info) << "Initializing ee mass spectrum via TPC electron identification analysis task ...";

    std::vector<double> momentumBinning = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                                           1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                                           4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 10.};
    const AxisSpec axisCounter{3, 0, +3, ""},
      axisCent{CentBinning, 0, 100, "centrality percentage"},
      axisVTrackX{400, -5., 5., "track vertex x (cm)", "track vertex x (cm)"},
      axisVTrackY{400, -5., 5., "track vertex y (cm)", "track vertex y (cm)"},
      axisVTrackZ{400, -20., 20., "track vertex z (cm)", "track vertex z (cm)"},
      axisE{EnergyBinning, 0, 10, "E (GeV)", "E (GeV)"},
      axisMassSpectrum{MassBinning, MassSpectraJpsiMin, MassSpectraJpsiMax, "M (GeV/c^{2})", "Mass e^{+}e^{-} (GeV/c^{2})"},
      axisMassSpectrumChiC{MassBinning, MassSpectraChicMin, MassSpectraChicMax, "M (GeV/c^{2})", "Mass e^{+}e^{-}#gamma (GeV/c^{2})"},
      axisMassSpectrumChiCNoJpsiErrors{MassBinning, MassSpectraChicMin, MassSpectraChicMax, "M (GeV/c^{2})", "Mass e^{+}e^{-}#gamma - Mass e^{+}e^{-} + Mass J/#psi (GeV/c^{2})"},
      axisMassSpectrumgammagamma{MassBinning, 0, 0.3, "M (GeV/c^{2})", "Mass #gamma#gamma (GeV/c^{2})"},
      axisTPC{1000, 0, 200, "TPC signal (dE/dx)"},
      axisPt{momentumBinning, "p_{T} (GeV/c)"},
      axisPtProbe{momentumBinning, "Probe p_{T} (GeV/c)"},
      axisPtBig{2000, 0, 20, "p_{T} (GeV/c)"},
      axisEta{600, -3., 3., "#eta"};

    mHistManager.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    mHistManager.add("centCounter", "centCounter", kTH1F, {axisCent});
    mHistManager.add("hTPCspectra", "TPC dE/dx spectra", HistType::kTH2F, {axisPt, axisTPC});
    mHistManager.add("hTPCspectra_isElectronRej", "isElectron with rejection | TPC dE/dx spectra", HistType::kTH2F, {axisPt, axisTPC});

    mHistManager.add("TPCee/h_MS_mp_v_pt_v_cent", "Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp} (from TPC candidates)", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("TPCee/h_MS_mm_v_pt_v_cent", "Mass e^{-}e^{-} vs momentum e^{-}e^{-} (from TPC candidates)", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("TPCee/h_MS_pp_v_pt_v_cent", "Mass e^{+}e^{+} vs momentum e^{+}e^{+} (from TPC candidates)", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});

    mHistManager.add("TPCee/h_MS_mp_kTVXinPHOS_v_pt_v_cent", "TVXinPHOS | Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp} (from TPC candidates)", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("TPCee/h_MS_mm_kTVXinPHOS_v_pt_v_cent", "TVXinPHOS | Mass e^{-}e^{-} vs momentum e^{-}e^{-} (from TPC candidates)", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("TPCee/h_MS_pp_kTVXinPHOS_v_pt_v_cent", "TVXinPHOS | Mass e^{+}e^{+} vs momentum e^{+}e^{+} (from TPC candidates)", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});

    mHistManager.add("TPCee/h_MS_mp_phosRange_v_pt_v_cent", "Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp} (from TPC candidates) with one e in phos acceptance range", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("TPCee/h_MS_mm_phosRange_v_pt_v_cent", "Mass e^{-}e^{-} vs momentum e^{-}e^{-} (from TPC candidates) with one e in phos acceptance range", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("TPCee/h_MS_pp_phosRange_v_pt_v_cent", "Mass e^{+}e^{+} vs momentum e^{+}e^{+} (from TPC candidates) with one e in phos acceptance range", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});

    mHistManager.add("TPCee/h_MS_mp_phosRange_kTVXinPHOS_v_pt_v_cent", "TVXinPHOS | Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp} (from TPC candidates) with one e in phos acceptance range", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("TPCee/h_MS_mm_phosRange_kTVXinPHOS_v_pt_v_cent", "TVXinPHOS | Mass e^{-}e^{-} vs momentum e^{-}e^{-} (from TPC candidates) with one e in phos acceptance range", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("TPCee/h_MS_pp_phosRange_kTVXinPHOS_v_pt_v_cent", "TVXinPHOS | Mass e^{+}e^{+} vs momentum e^{+}e^{+} (from TPC candidates) with one e in phos acceptance range", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});

    mHistManager.add("TPCeePhosGamma/h_MS_noMatches_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("TPCeePhosGamma/h_MS_noMatches_aroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("TPCeePhosGamma/h_MS_noMatches_aroundJpsi_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | DispOK", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("TPCeePhosGamma/h_MS_noMatches_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma | DispOK", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("TPCeePhosGamma/h_MS_aroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("TPCeePhosGamma/h_MS_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma | DispOK", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});

    mHistManager.add("TPCeePhosGamma/h_MS_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma (TPC candidates + Phos cluster)", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("TPCeePhosGamma/h_MS_v_cluE_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs cluster Energy left by the photon", HistType::kTH3F, {axisMassSpectrumChiC, axisE, axisCent});

    mHistManager.add("TPCeePhosGamma/h_minusee_MS_noMatches_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}#gamma", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("TPCeePhosGamma/h_minusee_MS_noMatches_aroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("TPCeePhosGamma/h_minusee_MS_noMatches_aroundJpsi_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | DispOK", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("TPCeePhosGamma/h_minusee_MS_noMatches_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}#gamma | DispOK", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("TPCeePhosGamma/h_minusee_MS_aroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("TPCeePhosGamma/h_minusee_MS_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}#gamma | DispOK", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});

    mHistManager.add("TPCeePhosGamma/h_minusee_MS_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp}  - Mass e^{#pm}e^{#mp} + Mass J/#psi vs momentum e^{#pm}e^{#mp}#gamma (TPC candidates + Phos cluster)", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("TPCeePhosGamma/h_minusee_MS_v_cluE_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp}  - Mass e^{#pm}e^{#mp} + Mass J/#psi vs cluster Energy left by the photon", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisE, axisCent});

    mHistManager.add("twoPhoton/MS_noCuts", "Mass vs Transverse Momentum for #gamma#gamma", HistType::kTH3F, {axisMassSpectrumgammagamma, axisPt, axisCent});
    mHistManager.add("twoPhoton/MS_noMatches", "Mass vs Transverse Momentum for #gamma#gamma excluding trackmatched clusters", HistType::kTH3F, {axisMassSpectrumgammagamma, axisPt, axisCent});

    mHistManager.add("TPCeff/h_eh_pp_mass_spectra_v_pt_v_cent", "Mass e^{+}h^{+} vs momentum e^{+}h^{+}", HistType::kTH3F, {axisMassSpectrum, axisPtProbe, axisCent});
    mHistManager.add("TPCeff/h_ee_pp_mass_spectra_v_pt_v_cent", "Mass e^{+}e^{+} vs momentum e^{+}e^{+}", HistType::kTH3F, {axisMassSpectrum, axisPtProbe, axisCent});
    mHistManager.add("TPCeff/h_eh_mm_mass_spectra_v_pt_v_cent", "Mass e^{-}h^{-} vs momentum e^{-}h^{-}", HistType::kTH3F, {axisMassSpectrum, axisPtProbe, axisCent});
    mHistManager.add("TPCeff/h_ee_mm_mass_spectra_v_pt_v_cent", "Mass e^{-}e^{-} vs momentum e^{-}e^{-}", HistType::kTH3F, {axisMassSpectrum, axisPtProbe, axisCent});
    mHistManager.add("TPCeff/h_eh_mp_mass_spectra_v_pt_v_cent", "Mass e^{#pm}h^{#mp} vs momentum e^{#pm}h^{#mp}", HistType::kTH3F, {axisMassSpectrum, axisPtProbe, axisCent});
    mHistManager.add("TPCeff/h_ee_mp_mass_spectra_v_pt_v_cent", "Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}", HistType::kTH3F, {axisMassSpectrum, axisPtProbe, axisCent});

    mHistManager.add("hTrackVX", "Track vertex coordinate X", HistType::kTH1F, {axisVTrackX});
    mHistManager.add("hTrackVY", "Track vertex coordinate Y", HistType::kTH1F, {axisVTrackY});
    mHistManager.add("hTrackVZ", "Track vertex coordinate Z", HistType::kTH1F, {axisVTrackZ});
    mHistManager.add("hTrackVX_Cut", "Track vertex coordinate X after cut", HistType::kTH1F, {axisVTrackX});
    mHistManager.add("hTrackVY_Cut", "Track vertex coordinate Y after cut", HistType::kTH1F, {axisVTrackY});
    mHistManager.add("hTrackVZ_Cut", "Track vertex coordinate Z after cut", HistType::kTH1F, {axisVTrackZ});

    mHistManager.add("hTrackPt", "Track pt", HistType::kTH1F, {axisPtBig});
    mHistManager.add("hTrackPt_Cut", "Track pt after cut", HistType::kTH1F, {axisPtBig});
    mHistManager.add("hTrackEta", "Track eta", HistType::kTH1F, {axisEta});
    mHistManager.add("hTrackEta_Cut", "Track eta after cut", HistType::kTH1F, {axisEta});

    if (isMC) {
      mHistManager.add("True/hTrackPt", "True Electron Track pt", HistType::kTH1F, {axisPtBig});
      mHistManager.add("True/hTPCspectra", "True Electron TPC dE/dx spectra", HistType::kTH2F, {axisPt, axisTPC});
      mHistManager.add("True/hTPCspectra_isElectronRej", "True Electron isElectron with rejection | TPC dE/dx spectra", HistType::kTH2F, {axisPt, axisTPC});
      mHistManager.add("True/TPCee/h_MS_mp_v_pt_v_cent", "Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp} (from TPC candidates)", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
      mHistManager.add("True/TPCee/h_MS_mm_v_pt_v_cent", "Mass e^{-}e^{-} vs momentum e^{-}e^{-} (from TPC candidates)", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
      mHistManager.add("True/TPCee/h_MS_pp_v_pt_v_cent", "Mass e^{+}e^{+} vs momentum e^{+}e^{+} (from TPC candidates)", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
      mHistManager.add("True/TPCeff/h_eh_pp_mass_spectra_v_pt_v_cent", "Mass e^{+}h^{+} vs momentum e^{+}h^{+}", HistType::kTH3F, {axisMassSpectrum, axisPtProbe, axisCent});
      mHistManager.add("True/TPCeff/h_ee_pp_mass_spectra_v_pt_v_cent", "Mass e^{+}e^{+} vs momentum e^{+}e^{+}", HistType::kTH3F, {axisMassSpectrum, axisPtProbe, axisCent});
      mHistManager.add("True/TPCeff/h_eh_mm_mass_spectra_v_pt_v_cent", "Mass e^{-}h^{-} vs momentum e^{-}h^{-}", HistType::kTH3F, {axisMassSpectrum, axisPtProbe, axisCent});
      mHistManager.add("True/TPCeff/h_ee_mm_mass_spectra_v_pt_v_cent", "Mass e^{-}e^{-} vs momentum e^{-}e^{-}", HistType::kTH3F, {axisMassSpectrum, axisPtProbe, axisCent});
      mHistManager.add("True/TPCeff/h_eh_mp_mass_spectra_v_pt_v_cent", "Mass e^{#pm}h^{#mp} vs momentum e^{#pm}h^{#mp}", HistType::kTH3F, {axisMassSpectrum, axisPtProbe, axisCent});
      mHistManager.add("True/TPCeff/h_ee_mp_mass_spectra_v_pt_v_cent", "Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}", HistType::kTH3F, {axisMassSpectrum, axisPtProbe, axisCent});
    }
  }

  void processData(SelCollisions::iterator const& collision,
                   aod::CaloClusters const& clusters,
                   MyTracks const& tracks,
                   o2::aod::PHOSMatchindexTable const& matches,
                   aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
      o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo == nullptr) {
        LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      bz = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      runNumber = bc.runNumber();
    }
    if (std::fabs(collision.posZ()) > mColMaxZ)
      return;

    float cent = -1.;
    switch (CentEst) {
      case FV0A:
        cent = collision.centFV0A();
        break;
      case FT0M:
        cent = collision.centFT0M();
        break;
      case FT0A:
        cent = collision.centFT0A();
        break;
      case FT0C:
        cent = collision.centFT0C();
        break;
      case FDDM:
        cent = collision.centFDDM();
        break;
      case NTPV:
        cent = collision.centNTPV();
        break;
    }
    mHistManager.fill(HIST("eventCounter"), 0.5);
    mHistManager.fill(HIST("centCounter"), cent);
    if ((isMC || collision.alias_bit(mEvSelTrig))) {
      mHistManager.fill(HIST("eventCounter"), 1.5);
    }
    if (isSel8) {
      if (!collision.sel8())
        return;
      mHistManager.fill(HIST("eventCounter"), 2.5);
    }

    auto isGoodElectronForSignal = [&](const MyTracks::iterator& track) -> bool {
      if (!track.has_collision() || !track.hasTPC())
        return false;
      if (track.pt() <= PtMin || track.pt() >= PtMax)
        return false;
      if (std::fabs(track.eta()) >= EtaMax)
        return false;
      if (std::fabs(track.dcaXY()) >= DCAxyMax)
        return false;
      if (std::fabs(track.dcaZ()) >= DCAzMax)
        return false;
      if (track.itsChi2NCl() >= ITSchi2Max)
        return false;
      if (track.tpcChi2NCl() >= TPCchi2Max)
        return false;
      if (!((track.itsClusterMap() & uint8_t(1)) > 0))
        return false;
      if (track.itsNCls() < ITSnclsMin || track.itsNCls() > ITSnclsMax)
        return false;
      if (track.tpcNClsFound() < TPCnclsMin || track.tpcNClsFound() > TPCnclsMax)
        return false;
      if (track.tpcNClsCrossedRows() < TPCnclsCRMin || track.tpcNClsCrossedRows() > TPCnclsCRMax)
        return false;

      bool isTPCElectron = (track.tpcNSigmaEl() > TPCNSigmaElMin) && (track.tpcNSigmaEl() < TPCNSigmaElMax);
      bool isTOFElectron = (track.tofNSigmaEl() > TOFNSigmaElMin) && (track.tofNSigmaEl() < TOFNSigmaElMax);
      if (!isTPCElectron && !isTOFElectron)
        return false;

      bool isPion = (track.tpcNSigmaPi() >= TPCNSigmaPiMin && track.tpcNSigmaPi() <= TPCNSigmaPiMax);
      bool isKaon = (track.tpcNSigmaKa() >= TPCNSigmaKaMin && track.tpcNSigmaKa() <= TPCNSigmaKaMax);
      bool isProton = (track.tpcNSigmaPr() >= TPCNSigmaPrMin && track.tpcNSigmaPr() <= TPCNSigmaPrMax);
      if (isPion || isKaon || isProton)
        return false;
      return true;
    };

    auto isGoodTagElectron = [&](const MyTracks::iterator& track) -> bool {
      if (!track.has_collision() || !track.hasTPC())
        return false;
      if (!((track.itsClusterMap() & uint8_t(1)) > 0))
        return false;
      if (track.itsChi2NCl() > ITSchi2Max || track.tpcChi2NCl() > TPCchi2Max)
        return false;
      if (track.itsNCls() < ITSnclsMin || track.itsNCls() > ITSnclsMax)
        return false;
      if (track.tpcNClsFound() < TPCnclsMin || track.tpcNClsFound() > TPCnclsMax)
        return false;
      if (track.tpcNClsCrossedRows() < TPCnclsCRMin || track.tpcNClsCrossedRows() > TPCnclsCRMax)
        return false;
      if (std::fabs(track.eta()) >= EtaMax)
        return false;
      if (std::fabs(track.dcaXY()) >= DCAxyMax)
        return false;
      if (std::fabs(track.dcaZ()) >= DCAzMax)
        return false;

      bool isTPCElectron = (track.tpcNSigmaEl() > TPCNSigmaElMin) && (track.tpcNSigmaEl() < TPCNSigmaElMax);
      bool isTOFElectron = (track.tofNSigmaEl() > TOFNSigmaElMin) && (track.tofNSigmaEl() < TOFNSigmaElMax);
      if (!isTPCElectron && !isTOFElectron)
        return false;

      bool isPionSignal = (track.tpcNSigmaPi() >= TPCNSigmaPiMin && track.tpcNSigmaPi() <= TPCNSigmaPiMax);
      bool isKaonSignal = (track.tpcNSigmaKa() >= TPCNSigmaKaMin && track.tpcNSigmaKa() <= TPCNSigmaKaMax);
      bool isProtonSignal = (track.tpcNSigmaPr() >= TPCNSigmaPrMin && track.tpcNSigmaPr() <= TPCNSigmaPrMax);
      if (isPionSignal || isKaonSignal || isProtonSignal)
        return false;
      return true;
    };

    auto isGoodProbeBaseTrack = [&](const MyTracks::iterator& track) -> bool {
      if (!track.has_collision() || !track.hasTPC())
        return false;
      if (!((track.itsClusterMap() & uint8_t(1)) > 0))
        return false;
      if (track.itsChi2NCl() > ITSchi2Max || track.tpcChi2NCl() > TPCchi2Max)
        return false;
      if (track.itsNCls() < ITSnclsMin || track.itsNCls() > ITSnclsMax)
        return false;
      if (track.tpcNClsFound() < TPCnclsMin || track.tpcNClsFound() > TPCnclsMax)
        return false;
      if (track.tpcNClsCrossedRows() < TPCnclsCRMin || track.tpcNClsCrossedRows() > TPCnclsCRMax)
        return false;
      if (std::fabs(track.dcaXY()) > DCAxyMax || std::fabs(track.dcaZ()) > DCAzMax)
        return false;
      if (std::fabs(track.eta()) >= EtaMax)
        return false;
      return true;
    };

    auto isProbeIdentifiedAsElectron = [&](const MyTracks::iterator& track) -> bool {
      if (!track.hasTPC())
        return false;
      bool isTPCElectron = (track.tpcNSigmaEl() > TPCNSigmaElMin) && (track.tpcNSigmaEl() < TPCNSigmaElMax);
      bool isTOFElectron = (track.tofNSigmaEl() > TOFNSigmaElMin) && (track.tofNSigmaEl() < TOFNSigmaElMax);
      if (!isTPCElectron && !isTOFElectron)
        return false;

      bool isPionSignal = (track.tpcNSigmaPi() >= TPCNSigmaPiMin && track.tpcNSigmaPi() <= TPCNSigmaPiMax);
      bool isKaonSignal = (track.tpcNSigmaKa() >= TPCNSigmaKaMin && track.tpcNSigmaKa() <= TPCNSigmaKaMax);
      bool isProtonSignal = (track.tpcNSigmaPr() >= TPCNSigmaPrMin && track.tpcNSigmaPr() <= TPCNSigmaPrMax);
      if (isPionSignal || isKaonSignal || isProtonSignal)
        return false;
      return true;
    };

    for (auto const& [track1_iterator, track2_iterator] : combinations(CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {
      if (track1_iterator.collisionId() != track2_iterator.collisionId()) {
        continue;
      }

      bool track1IsSignalE = isGoodElectronForSignal(track1_iterator);
      bool track2IsSignalE = isGoodElectronForSignal(track2_iterator);

      if (track1IsSignalE && track2IsSignalE) {
        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> fourVectorP1, fourVectorP2;
        fourVectorP1.SetPxPyPzE(track1_iterator.px(), track1_iterator.py(), track1_iterator.pz(), track1_iterator.energy(0));
        fourVectorP2.SetPxPyPzE(track2_iterator.px(), track2_iterator.py(), track2_iterator.pz(), track2_iterator.energy(0));

        bool inPhosEtaRange1 = std::fabs(track1_iterator.eta()) < PhosRangeEta;
        bool inPhosEtaRange2 = std::fabs(track2_iterator.eta()) < PhosRangeEta;
        bool inPhosPhiRange1 = (track1_iterator.phi() * TMath::RadToDeg() > PhosRangePhiMin && track1_iterator.phi() * TMath::RadToDeg() < PhosRangePhiMax);
        bool inPhosPhiRange2 = (track2_iterator.phi() * TMath::RadToDeg() > PhosRangePhiMin && track2_iterator.phi() * TMath::RadToDeg() < PhosRangePhiMax);
        bool inPhosRange = (inPhosEtaRange1 && inPhosPhiRange1) || (inPhosEtaRange2 && inPhosPhiRange2);

        double pairMass = (fourVectorP1 + fourVectorP2).M(), pairPt = (fourVectorP1 + fourVectorP2).Pt();

        if (track1_iterator.sign() == track2_iterator.sign()) {
          bool track1IsPositive = track1_iterator.sign() * bz > 0;
          if (track1IsPositive) {
            mHistManager.fill(HIST("TPCee/h_MS_pp_v_pt_v_cent"), pairMass, pairPt, cent);
            if ((isMC || collision.alias_bit(mEvSelTrig)))
              mHistManager.fill(HIST("TPCee/h_MS_pp_kTVXinPHOS_v_pt_v_cent"), pairMass, pairPt, cent);
            if (inPhosRange) {
              mHistManager.fill(HIST("TPCee/h_MS_pp_phosRange_v_pt_v_cent"), pairMass, pairPt, cent);
              if ((isMC || collision.alias_bit(mEvSelTrig)))
                mHistManager.fill(HIST("TPCee/h_MS_pp_phosRange_kTVXinPHOS_v_pt_v_cent"), pairMass, pairPt, cent);
            }
          } else {
            mHistManager.fill(HIST("TPCee/h_MS_mm_v_pt_v_cent"), pairMass, pairPt, cent);
            if ((isMC || collision.alias_bit(mEvSelTrig)))
              mHistManager.fill(HIST("TPCee/h_MS_mm_kTVXinPHOS_v_pt_v_cent"), pairMass, pairPt, cent);
            if (inPhosRange) {
              mHistManager.fill(HIST("TPCee/h_MS_mm_phosRange_v_pt_v_cent"), pairMass, pairPt, cent);
              if ((isMC || collision.alias_bit(mEvSelTrig)))
                mHistManager.fill(HIST("TPCee/h_MS_mm_phosRange_kTVXinPHOS_v_pt_v_cent"), pairMass, pairPt, cent);
            }
          }
        } else {
          mHistManager.fill(HIST("TPCee/h_MS_mp_v_pt_v_cent"), pairMass, pairPt, cent);
          if ((isMC || collision.alias_bit(mEvSelTrig)))
            mHistManager.fill(HIST("TPCee/h_MS_mp_kTVXinPHOS_v_pt_v_cent"), pairMass, pairPt, cent);
          if (inPhosRange) {
            mHistManager.fill(HIST("TPCee/h_MS_mp_phosRange_v_pt_v_cent"), pairMass, pairPt, cent);
            if ((isMC || collision.alias_bit(mEvSelTrig)))
              mHistManager.fill(HIST("TPCee/h_MS_mp_phosRange_kTVXinPHOS_v_pt_v_cent"), pairMass, pairPt, cent);
          }

          if ((isMC || collision.alias_bit(mEvSelTrig)) && clusters.size() != 0) {
            for (auto const& gamma : clusters) {
              float cluE = gamma.e();
              if (cluE < mMinCluE || cluE > mMaxCluE || gamma.ncell() < mMinCluNcell || gamma.time() > mMaxCluTime || gamma.time() < mMinCluTime)
                continue;
              bool matchFlag = false;
              bool isJpsi = (pairMass > eeMassMin && pairMass < eeMassMax);
              bool isDispOK = false;
              if (mSwapM20M02ForTestLambda)
                isDispOK = testLambda(cluE, gamma.m02(), gamma.m20(), mShowerShapeCutValue, mUseNegativeCrossTerm);
              else
                isDispOK = testLambda(cluE, gamma.m20(), gamma.m02(), mShowerShapeCutValue, mUseNegativeCrossTerm);
              for (auto const& match : matches) {
                if (gamma.index() == match.caloClusterId()) {
                  matchFlag = true;
                  break;
                }
              }
              ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> fourVectorP3;
              fourVectorP3.SetPxPyPzE(gamma.px(), gamma.py(), gamma.pz(), cluE);
              double tripletMass = (fourVectorP1 + fourVectorP2 + fourVectorP3).M();
              double tripletPt = (fourVectorP1 + fourVectorP2 + fourVectorP3).Pt();
              double tripletMinusPairPlusJpsiMass = tripletMass - pairMass + JpsiMass;

              mHistManager.fill(HIST("TPCeePhosGamma/h_MS_v_3pt_v_cent"), tripletMass, tripletPt, cent);
              mHistManager.fill(HIST("TPCeePhosGamma/h_minusee_MS_v_3pt_v_cent"), tripletMinusPairPlusJpsiMass, tripletPt, cent);
              mHistManager.fill(HIST("TPCeePhosGamma/h_MS_v_cluE_v_cent"), tripletMass, cluE, cent);
              mHistManager.fill(HIST("TPCeePhosGamma/h_minusee_MS_v_cluE_v_cent"), tripletMinusPairPlusJpsiMass, cluE, cent);

              if (!matchFlag) {
                mHistManager.fill(HIST("TPCeePhosGamma/h_MS_noMatches_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                mHistManager.fill(HIST("TPCeePhosGamma/h_minusee_MS_noMatches_v_3pt_v_cent"), tripletMinusPairPlusJpsiMass, tripletPt, cent);
                if (isJpsi) {
                  mHistManager.fill(HIST("TPCeePhosGamma/h_MS_noMatches_aroundJpsi_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                  mHistManager.fill(HIST("TPCeePhosGamma/h_minusee_MS_noMatches_aroundJpsi_v_3pt_v_cent"), tripletMinusPairPlusJpsiMass, tripletPt, cent);
                  if (isDispOK) {
                    mHistManager.fill(HIST("TPCeePhosGamma/h_MS_noMatches_aroundJpsi_DispOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                    mHistManager.fill(HIST("TPCeePhosGamma/h_minusee_MS_noMatches_aroundJpsi_DispOK_v_3pt_v_cent"), tripletMinusPairPlusJpsiMass, tripletPt, cent);
                  }
                }
                if (isDispOK) {
                  mHistManager.fill(HIST("TPCeePhosGamma/h_MS_noMatches_DispOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                  mHistManager.fill(HIST("TPCeePhosGamma/h_minusee_MS_noMatches_DispOK_v_3pt_v_cent"), tripletMinusPairPlusJpsiMass, tripletPt, cent);
                }
              }
              if (isJpsi) {
                mHistManager.fill(HIST("TPCeePhosGamma/h_MS_aroundJpsi_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                mHistManager.fill(HIST("TPCeePhosGamma/h_minusee_MS_aroundJpsi_v_3pt_v_cent"), tripletMinusPairPlusJpsiMass, tripletPt, cent);
              }
              if (isDispOK) {
                mHistManager.fill(HIST("TPCeePhosGamma/h_MS_DispOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                mHistManager.fill(HIST("TPCeePhosGamma/h_minusee_MS_DispOK_v_3pt_v_cent"), tripletMinusPairPlusJpsiMass, tripletPt, cent);
              }
            }
          }
        }
      }

      if (isGoodTagElectron(track1_iterator) && isGoodProbeBaseTrack(track2_iterator)) {
        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pTag1, pProbe2;
        pTag1.SetPxPyPzE(track1_iterator.px(), track1_iterator.py(), track1_iterator.pz(), track1_iterator.energy(0));
        pProbe2.SetPxPyPzE(track2_iterator.px(), track2_iterator.py(), track2_iterator.pz(), track2_iterator.energy(0));
        float massTag1Probe2 = (pTag1 + pProbe2).M();
        float ptProbe2 = track2_iterator.pt();
        bool tag1IsPositive = track1_iterator.sign() * bz > 0;

        if (track1_iterator.sign() == track2_iterator.sign()) {
          if (tag1IsPositive) {
            mHistManager.fill(HIST("TPCeff/h_eh_pp_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
          } else {
            mHistManager.fill(HIST("TPCeff/h_eh_mm_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
          }
        } else {
          mHistManager.fill(HIST("TPCeff/h_eh_mp_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
        }
        if (isProbeIdentifiedAsElectron(track2_iterator)) {
          if (track1_iterator.sign() == track2_iterator.sign()) {
            if (tag1IsPositive) {
              mHistManager.fill(HIST("TPCeff/h_ee_pp_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
            } else {
              mHistManager.fill(HIST("TPCeff/h_ee_mm_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
            }
          } else {
            mHistManager.fill(HIST("TPCeff/h_ee_mp_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
          }
        }
      }

      if (isGoodTagElectron(track2_iterator) && isGoodProbeBaseTrack(track1_iterator)) {
        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pTag2, pProbe1;
        pTag2.SetPxPyPzE(track2_iterator.px(), track2_iterator.py(), track2_iterator.pz(), track2_iterator.energy(0));
        pProbe1.SetPxPyPzE(track1_iterator.px(), track1_iterator.py(), track1_iterator.pz(), track1_iterator.energy(0));
        float massTag2Probe1 = (pTag2 + pProbe1).M();
        float ptProbe1 = track1_iterator.pt();
        bool tag2IsPositive = track2_iterator.sign() * bz > 0;

        if (track2_iterator.sign() == track1_iterator.sign()) {
          if (tag2IsPositive) {
            mHistManager.fill(HIST("TPCeff/h_eh_pp_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
          } else {
            mHistManager.fill(HIST("TPCeff/h_eh_mm_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
          }
        } else {
          mHistManager.fill(HIST("TPCeff/h_eh_mp_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
        }
        if (isProbeIdentifiedAsElectron(track1_iterator)) {
          if (track2_iterator.sign() == track1_iterator.sign()) {
            if (tag2IsPositive) {
              mHistManager.fill(HIST("TPCeff/h_ee_pp_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
            } else {
              mHistManager.fill(HIST("TPCeff/h_ee_mm_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
            }
          } else {
            mHistManager.fill(HIST("TPCeff/h_ee_mp_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
          }
        }
      }
    }

    for (auto const& gamma1 : clusters) {
      float cluE1 = gamma1.e();
      if (cluE1 < mMinCluE || gamma1.ncell() < mMinCluNcell || gamma1.time() > mMaxCluTime || gamma1.time() < mMinCluTime)
        continue;
      bool matchFlag1 = false;

      bool isDispOKClu1 = false;
      if (mSwapM20M02ForTestLambda)
        isDispOKClu1 = testLambda(cluE1, gamma1.m02(), gamma1.m20(), mShowerShapeCutValue, mUseNegativeCrossTerm);
      else
        isDispOKClu1 = testLambda(cluE1, gamma1.m20(), gamma1.m02(), mShowerShapeCutValue, mUseNegativeCrossTerm);
      if (!isDispOKClu1)
        continue;
      for (auto const& match : matches) {
        if (gamma1.index() == match.caloClusterId()) {
          matchFlag1 = true;
          break;
        }
      }
      for (auto const& gamma2 : clusters) {
        if (gamma1.index() >= gamma2.index())
          continue;
        float cluE2 = gamma2.e();
        if (cluE2 < mMinCluE || gamma2.ncell() < mMinCluNcell || gamma2.time() > mMaxCluTime || gamma2.time() < mMinCluTime)
          continue;
        bool isDispOKClu2 = false;
        if (mSwapM20M02ForTestLambda)
          isDispOKClu2 = testLambda(cluE2, gamma2.m02(), gamma2.m20(), mShowerShapeCutValue, mUseNegativeCrossTerm);
        else
          isDispOKClu2 = testLambda(cluE2, gamma2.m20(), gamma2.m02(), mShowerShapeCutValue, mUseNegativeCrossTerm);
        if (!isDispOKClu2)
          continue;
        bool matchFlag2 = false;
        for (auto const& match : matches) {
          if (gamma2.index() == match.caloClusterId()) {
            matchFlag2 = true;
            break;
          }
        }
        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> fourVectorG1, fourVectorG2;
        fourVectorG1.SetPxPyPzE(gamma1.px(), gamma1.py(), gamma1.pz(), cluE1);
        fourVectorG2.SetPxPyPzE(gamma2.px(), gamma2.py(), gamma2.pz(), cluE2);
        double pairMassGG = (fourVectorG1 + fourVectorG2).M();
        double pairPtGG = (fourVectorG1 + fourVectorG2).Pt();

        if (pairMassGG < mMassSpectrumLowerCutoff)
          continue;

        mHistManager.fill(HIST("twoPhoton/MS_noCuts"), pairMassGG, pairPtGG, cent);
        if (matchFlag1 || matchFlag2)
          continue;
        mHistManager.fill(HIST("twoPhoton/MS_noMatches"), pairMassGG, pairPtGG, cent);
      }
    }

    for (auto const& track : tracks) {
      mHistManager.fill(HIST("hTrackPt"), track.pt());
      mHistManager.fill(HIST("hTrackEta"), track.eta());
      mHistManager.fill(HIST("hTrackVX"), track.x());
      mHistManager.fill(HIST("hTrackVY"), track.y());
      mHistManager.fill(HIST("hTrackVZ"), track.z());
      mHistManager.fill(HIST("hTPCspectra"), track.pt(), track.tpcSignal());

      if (isGoodElectronForSignal(track)) {
        mHistManager.fill(HIST("hTPCspectra_isElectronRej"), track.pt(), track.tpcSignal());
        mHistManager.fill(HIST("hTrackPt_Cut"), track.pt());
        mHistManager.fill(HIST("hTrackEta_Cut"), track.eta());
        mHistManager.fill(HIST("hTrackVX_Cut"), track.x());
        mHistManager.fill(HIST("hTrackVY_Cut"), track.y());
        mHistManager.fill(HIST("hTrackVZ_Cut"), track.z());
      }
    }
  }
  PROCESS_SWITCH(TpcElIdMassSpectrum, processData, "process data", false);

  void processMC(SelCollisions::iterator const& collision,
                 aod::CaloClusters const& clusters,
                 soa::Join<MyTracks, aod::McTrackLabels> const& tracks,
                 o2::aod::PHOSMatchindexTable const& matches,
                 aod::McParticles const& mcParticles,
                 aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      LOG(info) << ">>>>>>>>>>>> Current run number: " << runNumber;
      o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
      if (grpo == nullptr) {
        LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      o2::base::Propagator::initFieldFromGRP(grpo);
      bz = o2::base::Propagator::Instance()->getNominalBz();
      LOG(info) << ">>>>>>>>>>>> Magnetic field: " << bz;
      runNumber = bc.runNumber();
    }
    if (std::fabs(collision.posZ()) > mColMaxZ)
      return;

    float cent = -1.;
    switch (CentEst) {
      case FV0A:
        cent = collision.centFV0A();
        break;
      case FT0M:
        cent = collision.centFT0M();
        break;
      case FT0A:
        cent = collision.centFT0A();
        break;
      case FT0C:
        cent = collision.centFT0C();
        break;
      case FDDM:
        cent = collision.centFDDM();
        break;
      case NTPV:
        cent = collision.centNTPV();
        break;
    }
    mHistManager.fill(HIST("eventCounter"), 0.5);
    mHistManager.fill(HIST("centCounter"), cent);
    if ((isMC || collision.alias_bit(mEvSelTrig))) {
      mHistManager.fill(HIST("eventCounter"), 1.5);
    }
    if (isSel8) {
      if (!collision.sel8())
        return;
      mHistManager.fill(HIST("eventCounter"), 2.5);
    }

    auto isGoodElectronForSignal = [&](auto const& track) -> bool {
      if (!track.has_collision() || !track.hasTPC())
        return false;
      if (track.pt() <= PtMin || track.pt() >= PtMax)
        return false;
      if (std::fabs(track.eta()) >= EtaMax)
        return false;
      if (std::fabs(track.dcaXY()) >= DCAxyMax)
        return false;
      if (std::fabs(track.dcaZ()) >= DCAzMax)
        return false;
      if (track.itsChi2NCl() >= ITSchi2Max)
        return false;
      if (track.tpcChi2NCl() >= TPCchi2Max)
        return false;
      if (!((track.itsClusterMap() & uint8_t(1)) > 0))
        return false;
      if (track.itsNCls() < ITSnclsMin || track.itsNCls() > ITSnclsMax)
        return false;
      if (track.tpcNClsFound() < TPCnclsMin || track.tpcNClsFound() > TPCnclsMax)
        return false;
      if (track.tpcNClsCrossedRows() < TPCnclsCRMin || track.tpcNClsCrossedRows() > TPCnclsCRMax)
        return false;

      bool isTPCElectron = (track.tpcNSigmaEl() > TPCNSigmaElMin) && (track.tpcNSigmaEl() < TPCNSigmaElMax);
      bool isTOFElectron = (track.tofNSigmaEl() > TOFNSigmaElMin) && (track.tofNSigmaEl() < TOFNSigmaElMax);
      if (!isTPCElectron && !isTOFElectron)
        return false;

      bool isPion = (track.tpcNSigmaPi() >= TPCNSigmaPiMin && track.tpcNSigmaPi() <= TPCNSigmaPiMax);
      bool isKaon = (track.tpcNSigmaKa() >= TPCNSigmaKaMin && track.tpcNSigmaKa() <= TPCNSigmaKaMax);
      bool isProton = (track.tpcNSigmaPr() >= TPCNSigmaPrMin && track.tpcNSigmaPr() <= TPCNSigmaPrMax);
      if (isPion || isKaon || isProton)
        return false;
      return true;
    };

    auto isGoodTagElectron = [&](auto const& track) -> bool {
      if (!track.has_collision() || !track.hasTPC())
        return false;
      if (!((track.itsClusterMap() & uint8_t(1)) > 0))
        return false;
      if (track.itsChi2NCl() > ITSchi2Max || track.tpcChi2NCl() > TPCchi2Max)
        return false;
      if (track.itsNCls() < ITSnclsMin || track.itsNCls() > ITSnclsMax)
        return false;
      if (track.tpcNClsFound() < TPCnclsMin || track.tpcNClsFound() > TPCnclsMax)
        return false;
      if (track.tpcNClsCrossedRows() < TPCnclsCRMin || track.tpcNClsCrossedRows() > TPCnclsCRMax)
        return false;
      if (std::fabs(track.eta()) >= EtaMax)
        return false;
      if (std::fabs(track.dcaXY()) >= DCAxyMax)
        return false;
      if (std::fabs(track.dcaZ()) >= DCAzMax)
        return false;

      bool isTPCElectron = (track.tpcNSigmaEl() > TPCNSigmaElMin) && (track.tpcNSigmaEl() < TPCNSigmaElMax);
      bool isTOFElectron = (track.tofNSigmaEl() > TOFNSigmaElMin) && (track.tofNSigmaEl() < TOFNSigmaElMax);
      if (!isTPCElectron && !isTOFElectron)
        return false;

      bool isPionSignal = (track.tpcNSigmaPi() >= TPCNSigmaPiMin && track.tpcNSigmaPi() <= TPCNSigmaPiMax);
      bool isKaonSignal = (track.tpcNSigmaKa() >= TPCNSigmaKaMin && track.tpcNSigmaKa() <= TPCNSigmaKaMax);
      bool isProtonSignal = (track.tpcNSigmaPr() >= TPCNSigmaPrMin && track.tpcNSigmaPr() <= TPCNSigmaPrMax);
      if (isPionSignal || isKaonSignal || isProtonSignal)
        return false;
      return true;
    };

    auto isGoodProbeBaseTrack = [&](auto const& track) -> bool {
      if (!track.has_collision() || !track.hasTPC())
        return false;
      if (!((track.itsClusterMap() & uint8_t(1)) > 0))
        return false;
      if (track.itsChi2NCl() > ITSchi2Max || track.tpcChi2NCl() > TPCchi2Max)
        return false;
      if (track.itsNCls() < ITSnclsMin || track.itsNCls() > ITSnclsMax)
        return false;
      if (track.tpcNClsFound() < TPCnclsMin || track.tpcNClsFound() > TPCnclsMax)
        return false;
      if (track.tpcNClsCrossedRows() < TPCnclsCRMin || track.tpcNClsCrossedRows() > TPCnclsCRMax)
        return false;
      if (std::fabs(track.dcaXY()) > DCAxyMax || std::fabs(track.dcaZ()) > DCAzMax)
        return false;
      if (std::fabs(track.eta()) >= EtaMax)
        return false;
      return true;
    };

    auto isProbeIdentifiedAsElectron = [&](auto const& track) -> bool {
      if (!track.hasTPC())
        return false;
      bool isTPCElectron = (track.tpcNSigmaEl() > TPCNSigmaElMin) && (track.tpcNSigmaEl() < TPCNSigmaElMax);
      bool isTOFElectron = (track.tofNSigmaEl() > TOFNSigmaElMin) && (track.tofNSigmaEl() < TOFNSigmaElMax);
      if (!isTPCElectron && !isTOFElectron)
        return false;

      bool isPionSignal = (track.tpcNSigmaPi() >= TPCNSigmaPiMin && track.tpcNSigmaPi() <= TPCNSigmaPiMax);
      bool isKaonSignal = (track.tpcNSigmaKa() >= TPCNSigmaKaMin && track.tpcNSigmaKa() <= TPCNSigmaKaMax);
      bool isProtonSignal = (track.tpcNSigmaPr() >= TPCNSigmaPrMin && track.tpcNSigmaPr() <= TPCNSigmaPrMax);
      if (isPionSignal || isKaonSignal || isProtonSignal)
        return false;
      return true;
    };

    for (auto const& [track1_iterator, track2_iterator] : combinations(CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {
      if (track1_iterator.collisionId() != track2_iterator.collisionId()) {
        continue;
      }

      bool track1IsTrueE = false;
      auto mcLabel1 = track1_iterator.mcParticleId();
      if (mcLabel1 > -1 && mcLabel1 < mcParticles.size()) {
        auto mcpart = mcParticles.iteratorAt(mcLabel1);
        if (std::abs(mcpart.pdgCode()) == PDG_t::kElectron) {
          track1IsTrueE = true;
        }
      }
      bool track2IsTrueE = false;
      auto mcLabel2 = track2_iterator.mcParticleId();
      if (mcLabel2 > -1 && mcLabel2 < mcParticles.size()) {
        auto mcpart = mcParticles.iteratorAt(mcLabel2);
        if (std::abs(mcpart.pdgCode()) == PDG_t::kElectron) {
          track2IsTrueE = true;
        }
      }

      bool track1IsSignalE = isGoodElectronForSignal(track1_iterator);
      bool track2IsSignalE = isGoodElectronForSignal(track2_iterator);

      if (track1IsSignalE && track2IsSignalE) {
        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> fourVectorP1, fourVectorP2;
        fourVectorP1.SetPxPyPzE(track1_iterator.px(), track1_iterator.py(), track1_iterator.pz(), track1_iterator.energy(0));
        fourVectorP2.SetPxPyPzE(track2_iterator.px(), track2_iterator.py(), track2_iterator.pz(), track2_iterator.energy(0));

        bool inPhosEtaRange1 = std::fabs(track1_iterator.eta()) < PhosRangeEta;
        bool inPhosEtaRange2 = std::fabs(track2_iterator.eta()) < PhosRangeEta;
        bool inPhosPhiRange1 = (track1_iterator.phi() * TMath::RadToDeg() > PhosRangePhiMin && track1_iterator.phi() * TMath::RadToDeg() < PhosRangePhiMax);
        bool inPhosPhiRange2 = (track2_iterator.phi() * TMath::RadToDeg() > PhosRangePhiMin && track2_iterator.phi() * TMath::RadToDeg() < PhosRangePhiMax);
        bool inPhosRange = (inPhosEtaRange1 && inPhosPhiRange1) || (inPhosEtaRange2 && inPhosPhiRange2);

        double pairMass = (fourVectorP1 + fourVectorP2).M(), pairPt = (fourVectorP1 + fourVectorP2).Pt();

        if (track1_iterator.sign() == track2_iterator.sign()) {
          bool track1IsPositive = track1_iterator.sign() * bz > 0;
          if (track1IsPositive) {
            mHistManager.fill(HIST("TPCee/h_MS_pp_v_pt_v_cent"), pairMass, pairPt, cent);
            if ((isMC || collision.alias_bit(mEvSelTrig)))
              mHistManager.fill(HIST("TPCee/h_MS_pp_kTVXinPHOS_v_pt_v_cent"), pairMass, pairPt, cent);
            if (inPhosRange) {
              mHistManager.fill(HIST("TPCee/h_MS_pp_phosRange_v_pt_v_cent"), pairMass, pairPt, cent);
              if ((isMC || collision.alias_bit(mEvSelTrig)))
                mHistManager.fill(HIST("TPCee/h_MS_pp_phosRange_kTVXinPHOS_v_pt_v_cent"), pairMass, pairPt, cent);
            }
          } else {
            mHistManager.fill(HIST("TPCee/h_MS_mm_v_pt_v_cent"), pairMass, pairPt, cent);
            if ((isMC || collision.alias_bit(mEvSelTrig)))
              mHistManager.fill(HIST("TPCee/h_MS_mm_kTVXinPHOS_v_pt_v_cent"), pairMass, pairPt, cent);
            if (inPhosRange) {
              mHistManager.fill(HIST("TPCee/h_MS_mm_phosRange_v_pt_v_cent"), pairMass, pairPt, cent);
              if ((isMC || collision.alias_bit(mEvSelTrig)))
                mHistManager.fill(HIST("TPCee/h_MS_mm_phosRange_kTVXinPHOS_v_pt_v_cent"), pairMass, pairPt, cent);
            }
          }
        } else {
          mHistManager.fill(HIST("TPCee/h_MS_mp_v_pt_v_cent"), pairMass, pairPt, cent);
          if ((isMC || collision.alias_bit(mEvSelTrig)))
            mHistManager.fill(HIST("TPCee/h_MS_mp_kTVXinPHOS_v_pt_v_cent"), pairMass, pairPt, cent);
          if (inPhosRange) {
            mHistManager.fill(HIST("TPCee/h_MS_mp_phosRange_v_pt_v_cent"), pairMass, pairPt, cent);
            if ((isMC || collision.alias_bit(mEvSelTrig)))
              mHistManager.fill(HIST("TPCee/h_MS_mp_phosRange_kTVXinPHOS_v_pt_v_cent"), pairMass, pairPt, cent);
          }

          if ((isMC || collision.alias_bit(mEvSelTrig)) && clusters.size() != 0) {
            for (auto const& gamma : clusters) {
              float cluE = gamma.e();
              if (cluE < mMinCluE || cluE > mMaxCluE || gamma.ncell() < mMinCluNcell || gamma.time() > mMaxCluTime || gamma.time() < mMinCluTime)
                continue;
              bool matchFlag = false;
              bool isJpsi = (pairMass > eeMassMin && pairMass < eeMassMax);
              bool isDispOK = false;
              if (mSwapM20M02ForTestLambda)
                isDispOK = testLambda(cluE, gamma.m02(), gamma.m20(), mShowerShapeCutValue, mUseNegativeCrossTerm);
              else
                isDispOK = testLambda(cluE, gamma.m20(), gamma.m02(), mShowerShapeCutValue, mUseNegativeCrossTerm);
              for (auto const& match : matches) {
                if (gamma.index() == match.caloClusterId()) {
                  matchFlag = true;
                  break;
                }
              }
              ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> fourVectorP3;
              fourVectorP3.SetPxPyPzE(gamma.px(), gamma.py(), gamma.pz(), cluE);
              double tripletMass = (fourVectorP1 + fourVectorP2 + fourVectorP3).M();
              double tripletPt = (fourVectorP1 + fourVectorP2 + fourVectorP3).Pt();
              double tripletMinusPairPlusJpsiMass = tripletMass - pairMass + JpsiMass;

              mHistManager.fill(HIST("TPCeePhosGamma/h_MS_v_3pt_v_cent"), tripletMass, tripletPt, cent);
              mHistManager.fill(HIST("TPCeePhosGamma/h_minusee_MS_v_3pt_v_cent"), tripletMinusPairPlusJpsiMass, tripletPt, cent);
              mHistManager.fill(HIST("TPCeePhosGamma/h_MS_v_cluE_v_cent"), tripletMass, cluE, cent);
              mHistManager.fill(HIST("TPCeePhosGamma/h_minusee_MS_v_cluE_v_cent"), tripletMinusPairPlusJpsiMass, cluE, cent);

              if (!matchFlag) {
                mHistManager.fill(HIST("TPCeePhosGamma/h_MS_noMatches_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                mHistManager.fill(HIST("TPCeePhosGamma/h_minusee_MS_noMatches_v_3pt_v_cent"), tripletMinusPairPlusJpsiMass, tripletPt, cent);
                if (isJpsi) {
                  mHistManager.fill(HIST("TPCeePhosGamma/h_MS_noMatches_aroundJpsi_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                  mHistManager.fill(HIST("TPCeePhosGamma/h_minusee_MS_noMatches_aroundJpsi_v_3pt_v_cent"), tripletMinusPairPlusJpsiMass, tripletPt, cent);
                  if (isDispOK) {
                    mHistManager.fill(HIST("TPCeePhosGamma/h_MS_noMatches_aroundJpsi_DispOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                    mHistManager.fill(HIST("TPCeePhosGamma/h_minusee_MS_noMatches_aroundJpsi_DispOK_v_3pt_v_cent"), tripletMinusPairPlusJpsiMass, tripletPt, cent);
                  }
                }
                if (isDispOK) {
                  mHistManager.fill(HIST("TPCeePhosGamma/h_MS_noMatches_DispOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                  mHistManager.fill(HIST("TPCeePhosGamma/h_minusee_MS_noMatches_DispOK_v_3pt_v_cent"), tripletMinusPairPlusJpsiMass, tripletPt, cent);
                }
              }
              if (isJpsi) {
                mHistManager.fill(HIST("TPCeePhosGamma/h_MS_aroundJpsi_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                mHistManager.fill(HIST("TPCeePhosGamma/h_minusee_MS_aroundJpsi_v_3pt_v_cent"), tripletMinusPairPlusJpsiMass, tripletPt, cent);
              }
              if (isDispOK) {
                mHistManager.fill(HIST("TPCeePhosGamma/h_MS_DispOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                mHistManager.fill(HIST("TPCeePhosGamma/h_minusee_MS_DispOK_v_3pt_v_cent"), tripletMinusPairPlusJpsiMass, tripletPt, cent);
              }
            }
          }
        }

        if (track1IsTrueE && track2IsTrueE) {
          if (track1_iterator.sign() == track2_iterator.sign()) {
            bool track1IsPositive = track1_iterator.sign() * bz > 0;
            if (track1IsPositive) {
              mHistManager.fill(HIST("True/TPCee/h_MS_pp_v_pt_v_cent"), pairMass, pairPt, cent);
            } else {
              mHistManager.fill(HIST("True/TPCee/h_MS_mm_v_pt_v_cent"), pairMass, pairPt, cent);
            }
          } else {
            mHistManager.fill(HIST("True/TPCee/h_MS_mp_v_pt_v_cent"), pairMass, pairPt, cent);
          }
        }
      }

      if (isGoodTagElectron(track1_iterator) && isGoodProbeBaseTrack(track2_iterator)) {
        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pTag1, pProbe2;
        pTag1.SetPxPyPzE(track1_iterator.px(), track1_iterator.py(), track1_iterator.pz(), track1_iterator.energy(0));
        pProbe2.SetPxPyPzE(track2_iterator.px(), track2_iterator.py(), track2_iterator.pz(), track2_iterator.energy(0));
        float massTag1Probe2 = (pTag1 + pProbe2).M();
        float ptProbe2 = track2_iterator.pt();
        bool tag1IsPositive = track1_iterator.sign() * bz > 0;

        if (track1_iterator.sign() == track2_iterator.sign()) {
          if (tag1IsPositive) {
            mHistManager.fill(HIST("TPCeff/h_eh_pp_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
          } else {
            mHistManager.fill(HIST("TPCeff/h_eh_mm_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
          }
        } else {
          mHistManager.fill(HIST("TPCeff/h_eh_mp_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
        }
        if (isProbeIdentifiedAsElectron(track2_iterator)) {
          if (track1_iterator.sign() == track2_iterator.sign()) {
            if (tag1IsPositive) {
              mHistManager.fill(HIST("TPCeff/h_ee_pp_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
            } else {
              mHistManager.fill(HIST("TPCeff/h_ee_mm_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
            }
          } else {
            mHistManager.fill(HIST("TPCeff/h_ee_mp_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
          }
        }
        // Fill MC Truth Efficiency Histograms
        if (track1IsTrueE && track2IsTrueE) {
          if (track1_iterator.sign() == track2_iterator.sign()) {
            if (tag1IsPositive) {
              mHistManager.fill(HIST("True/TPCeff/h_eh_pp_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
            } else {
              mHistManager.fill(HIST("True/TPCeff/h_eh_mm_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
            }
          } else {
            mHistManager.fill(HIST("True/TPCeff/h_eh_mp_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
          }
          if (isProbeIdentifiedAsElectron(track2_iterator)) {
            if (track1_iterator.sign() == track2_iterator.sign()) {
              if (tag1IsPositive) {
                mHistManager.fill(HIST("True/TPCeff/h_ee_pp_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
              } else {
                mHistManager.fill(HIST("True/TPCeff/h_ee_mm_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
              }
            } else {
              mHistManager.fill(HIST("True/TPCeff/h_ee_mp_mass_spectra_v_pt_v_cent"), massTag1Probe2, ptProbe2, cent);
            }
          }
        }
      }

      if (isGoodTagElectron(track2_iterator) && isGoodProbeBaseTrack(track1_iterator)) {
        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> pTag2, pProbe1;
        pTag2.SetPxPyPzE(track2_iterator.px(), track2_iterator.py(), track2_iterator.pz(), track2_iterator.energy(0));
        pProbe1.SetPxPyPzE(track1_iterator.px(), track1_iterator.py(), track1_iterator.pz(), track1_iterator.energy(0));
        float massTag2Probe1 = (pTag2 + pProbe1).M();
        float ptProbe1 = track1_iterator.pt();
        bool tag2IsPositive = track2_iterator.sign() * bz > 0;

        if (track2_iterator.sign() == track1_iterator.sign()) {
          if (tag2IsPositive) {
            mHistManager.fill(HIST("TPCeff/h_eh_pp_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
          } else {
            mHistManager.fill(HIST("TPCeff/h_eh_mm_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
          }
        } else {
          mHistManager.fill(HIST("TPCeff/h_eh_mp_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
        }
        if (isProbeIdentifiedAsElectron(track1_iterator)) {
          if (track2_iterator.sign() == track1_iterator.sign()) {
            if (tag2IsPositive) {
              mHistManager.fill(HIST("TPCeff/h_ee_pp_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
            } else {
              mHistManager.fill(HIST("TPCeff/h_ee_mm_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
            }
          } else {
            mHistManager.fill(HIST("TPCeff/h_ee_mp_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
          }
        }
        // Fill MC Truth Efficiency Histograms (Tag 2, Probe 1)
        if (track1IsTrueE && track2IsTrueE) {
          if (track2_iterator.sign() == track1_iterator.sign()) {
            if (tag2IsPositive) {
              mHistManager.fill(HIST("True/TPCeff/h_eh_pp_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
            } else {
              mHistManager.fill(HIST("True/TPCeff/h_eh_mm_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
            }
          } else {
            mHistManager.fill(HIST("True/TPCeff/h_eh_mp_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
          }
          if (isProbeIdentifiedAsElectron(track1_iterator)) {
            if (track2_iterator.sign() == track1_iterator.sign()) {
              if (tag2IsPositive) {
                mHistManager.fill(HIST("True/TPCeff/h_ee_pp_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
              } else {
                mHistManager.fill(HIST("True/TPCeff/h_ee_mm_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
              }
            } else {
              mHistManager.fill(HIST("True/TPCeff/h_ee_mp_mass_spectra_v_pt_v_cent"), massTag2Probe1, ptProbe1, cent);
            }
          }
        }
      }
    }

    for (auto const& gamma1 : clusters) {
      float cluE1 = gamma1.e();
      if (cluE1 < mMinCluE || gamma1.ncell() < mMinCluNcell || gamma1.time() > mMaxCluTime || gamma1.time() < mMinCluTime)
        continue;
      bool matchFlag1 = false;

      bool isDispOKClu1 = false;
      if (mSwapM20M02ForTestLambda)
        isDispOKClu1 = testLambda(cluE1, gamma1.m02(), gamma1.m20(), mShowerShapeCutValue, mUseNegativeCrossTerm);
      else
        isDispOKClu1 = testLambda(cluE1, gamma1.m20(), gamma1.m02(), mShowerShapeCutValue, mUseNegativeCrossTerm);
      if (!isDispOKClu1)
        continue;
      for (auto const& match : matches) {
        if (gamma1.index() == match.caloClusterId()) {
          matchFlag1 = true;
          break;
        }
      }
      for (auto const& gamma2 : clusters) {
        if (gamma1.index() >= gamma2.index())
          continue;
        float cluE2 = gamma2.e();
        if (cluE2 < mMinCluE || gamma2.ncell() < mMinCluNcell || gamma2.time() > mMaxCluTime || gamma2.time() < mMinCluTime)
          continue;
        bool isDispOKClu2 = false;
        if (mSwapM20M02ForTestLambda)
          isDispOKClu2 = testLambda(cluE2, gamma2.m02(), gamma2.m20(), mShowerShapeCutValue, mUseNegativeCrossTerm);
        else
          isDispOKClu2 = testLambda(cluE2, gamma2.m20(), gamma2.m02(), mShowerShapeCutValue, mUseNegativeCrossTerm);
        if (!isDispOKClu2)
          continue;
        bool matchFlag2 = false;
        for (auto const& match : matches) {
          if (gamma2.index() == match.caloClusterId()) {
            matchFlag2 = true;
            break;
          }
        }
        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> fourVectorG1, fourVectorG2;
        fourVectorG1.SetPxPyPzE(gamma1.px(), gamma1.py(), gamma1.pz(), cluE1);
        fourVectorG2.SetPxPyPzE(gamma2.px(), gamma2.py(), gamma2.pz(), cluE2);
        double pairMassGG = (fourVectorG1 + fourVectorG2).M();
        double pairPtGG = (fourVectorG1 + fourVectorG2).Pt();

        if (pairMassGG < mMassSpectrumLowerCutoff)
          continue;

        mHistManager.fill(HIST("twoPhoton/MS_noCuts"), pairMassGG, pairPtGG, cent);
        if (matchFlag1 || matchFlag2)
          continue;
        mHistManager.fill(HIST("twoPhoton/MS_noMatches"), pairMassGG, pairPtGG, cent);
      }
    }

    for (auto const& track : tracks) {
      bool isTrueElectron = false;
      auto mcLabel = track.mcParticleId();
      if (mcLabel > -1 && mcLabel < mcParticles.size()) {
        auto mcpart = mcParticles.iteratorAt(mcLabel);
        if (std::abs(mcpart.pdgCode()) == PDG_t::kElectron) {
          isTrueElectron = true;
        }
      }

      mHistManager.fill(HIST("hTrackPt"), track.pt());
      mHistManager.fill(HIST("hTrackEta"), track.eta());
      mHistManager.fill(HIST("hTrackVX"), track.x());
      mHistManager.fill(HIST("hTrackVY"), track.y());
      mHistManager.fill(HIST("hTrackVZ"), track.z());
      mHistManager.fill(HIST("hTPCspectra"), track.pt(), track.tpcSignal());
      if (isTrueElectron) {
        mHistManager.fill(HIST("True/hTrackPt"), track.pt());
        mHistManager.fill(HIST("True/hTPCspectra"), track.pt(), track.tpcSignal());
      }

      if (isGoodElectronForSignal(track)) {
        mHistManager.fill(HIST("hTPCspectra_isElectronRej"), track.pt(), track.tpcSignal());
        if (isTrueElectron) {
          mHistManager.fill(HIST("True/hTPCspectra_isElectronRej"), track.pt(), track.tpcSignal());
        }
        mHistManager.fill(HIST("hTrackPt_Cut"), track.pt());
        mHistManager.fill(HIST("hTrackEta_Cut"), track.eta());
        mHistManager.fill(HIST("hTrackVX_Cut"), track.x());
        mHistManager.fill(HIST("hTrackVY_Cut"), track.y());
        mHistManager.fill(HIST("hTrackVZ_Cut"), track.z());
      }
    }
  }

  PROCESS_SWITCH(TpcElIdMassSpectrum, processMC, "process mc", false);
  void processDummy(SelCollisions::iterator const&)
  {
  }
  PROCESS_SWITCH(TpcElIdMassSpectrum, processDummy, "Dummy process", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{
    adaptAnalysisTask<PhosElId>(cfgc),
    adaptAnalysisTask<MassSpectra>(cfgc),
    adaptAnalysisTask<TpcElIdMassSpectrum>(cfgc)};
  return workflow;
}
