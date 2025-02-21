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

#include <climits>
#include <cstdlib>
#include <map>
#include <memory>
#include <vector>
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/CaloClusters.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "ReconstructionDataFormats/TrackParametrization.h"
#include "Framework/ConfigParamSpec.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"
#include "PHOSBase/Geometry.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CommonDataFormat/InteractionRecord.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPLHCIFData.h"
#include "DetectorsBase/Propagator.h"
#include "TF1.h"
#include "TLorentzVector.h"

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

struct PhosElId {

  Produces<o2::aod::PHOSMatchindexTable> phosMatch;

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                             aod::TracksDCACov, aod::pidTOFFullEl, aod::pidTPCFullEl,
                             aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
  Configurable<float> mMinCluE{"mMinCluE", 0.3, "Minimum cluster energy for analysis"},
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
    cfgEtaMax{"cfgEtaMax", {0.8f}, "eta ranges"},
    cfgPtMin{"cfgPtMin", {0.2f}, "pt min"},
    cfgPtMax{"cfgPtMax", {20.f}, "pt max"},
    cfgDCAxyMax{"cfgDCAxyMax", {3.f}, "dcaxy max"},
    cfgDCAzMax{"cfgDCAzMax", {3.f}, "dcaz max"},
    cfgITSchi2Max{"cfgITSchi2Max", {5.f}, "its chi2 max"},
    cfgITSnclsMin{"cfgITSnclsMin", {4.5f}, "min number of ITS clusters"},
    cfgITSnclsMax{"cfgITSnclsMax", {7.5f}, "max number of ITS clusters"},
    cfgTPCchi2Max{"cfgTPCchi2Max", {4.f}, "tpc chi2 max"},
    cfgTPCnclsMin{"cfgTPCnclsMin", {90.f}, "min number of TPC clusters"},
    cfgTPCnclsMax{"cfgTPCnclsMax", {170.f}, "max number of TPC clusters"},
    cfgTPCnclsCRMin{"cfgTPCnclsCRMin", {80.f}, "min number of TPC crossed rows"},
    cfgTPCnclsCRMax{"cfgTPCnclsCRMax", {161.f}, "max number of TPC crossed rows"},
    cfgTPCNSigmaElMin{"cfgTPCNSigmaElMin", {-3.f}, "min TPC nsigma e for inclusion"},
    cfgTPCNSigmaElMax{"cfgTPCNSigmaElMax", {2.f}, "max TPC nsigma e for inclusion"},
    cfgTPCNSigmaPiMin{"cfgTPCNSigmaPiMin", {-3.f}, "min TPC nsigma pion for exclusion"},
    cfgTPCNSigmaPiMax{"cfgTPCNSigmaPiMax", {3.5f}, "max TPC nsigma pion for exclusion"},
    cfgTPCNSigmaPrMin{"cfgTPCNSigmaPrMin", {-3.f}, "min TPC nsigma proton for exclusion"},
    cfgTPCNSigmaPrMax{"cfgTPCNSigmaPrMax", {4.f}, "max TPC nsigma proton for exclusion"},
    cfgTPCNSigmaKaMin{"cfgTPCNSigmaKaMin", {-3.f}, "min TPC nsigma kaon for exclusion"},
    cfgTPCNSigmaKaMax{"cfgTPCNSigmaKaMax", {4.f}, "max TPC nsigma kaon for exclusion"},
    cfgTOFNSigmaElMin{"cfgTOFNSigmaElMin", {-3.f}, "min TOF nsigma e for inclusion"},
    cfgTOFNSigmaElMax{"cfgTOFNSigmaElMax", {3.f}, "max TOF nsigma e for inclusion"},
    cfgNsigmaTrackMatch{"cfgNsigmaTrackMatch", {2.f}, "PHOS Track Matching Nsigma for inclusion"};

  Configurable<int> mEvSelTrig{"mEvSelTrig", kTVXinPHOS, "Select events with this trigger"},
    mMinCluNcell{"minCluNcell", 3, "min cells in cluster"},
    nBinsDeltaX{"nBinsDeltaX", 500, "N bins for track and cluster coordinate delta"},
    nBinsDeltaZ{"nBinsDeltaZ", 500, "N bins for track and cluster coordinate delta"},
    nBinsEp{"nBinsEp", 400, "N bins for E/p histograms"};

  Configurable<std::vector<float>> pSigmadZ{"pSigmadZ", {0.642, 0., 1.77, 2.725, 0.}, "parameters for sigma dz function"},
    pSigmadX{"pSigmadX", {2.17769, 1.60275, 2.24136}, "parameters for sigma dx function"},
    pPhosShiftZ{"pPhosShiftZ", {4.78838, 2.75138, 1.40825, 2.28735}, "Phos coordinate centering Z per module"},
    pPhosShiftX{"pPhosShiftX", {2.158702, -1.526772, -0.814658, -1.852678}, "Phos coordinate centering X per module"},
    pMeandXPosMod1{"pMeandXPosMod1", {-10.57, -0.42, 1.06}, "parameters for mean dx function on module 1 for positive tracks"},
    pMeandXPosMod2{"pMeandXPosMod2", {-8.1, -0.42, 1.14}, "parameters for mean dx function on module 2 for positive tracks"},
    pMeandXPosMod3{"pMeandXPosMod3", {-8.34, -0.42, 1.04}, "parameters for mean dx function on module 3 for positive tracks"},
    pMeandXPosMod4{"pMeandXPosMod4", {-7.38, -0.42, 1.17}, "parameters for mean dx function on module 4 for positive tracks"},
    pMeandXNegMod1{"pMeandXNegMod1", {9.92, -0.42, 1.29}, "parameters for mean dx function on module 1 for negative tracks"},
    pMeandXNegMod2{"pMeandXNegMod2", {7.82, -0.4, 1.34}, "parameters for mean dx function on module 2 for negative tracks"},
    pMeandXNegMod3{"pMeandXNegMod3", {8.45, -0.33, 1.5}, "parameters for mean dx function on module 3 for negative tracks"},
    pMeandXNegMod4{"pMeandXNegMod4", {7.5, -0.42, 1.25}, "parameters for mean dx function on module 4 for negative tracks"};

  Filter ptFilter = (aod::track::pt > cfgPtMin) && (aod::track::pt < cfgPtMax);
  Filter etafilter = nabs(aod::track::eta) < cfgEtaMax;
  Filter dcaxyfilter = nabs(aod::track::dcaXY) < cfgDCAxyMax;
  Filter dcazfilter = nabs(aod::track::dcaZ) < cfgDCAzMax;
  Filter itschi2filter = aod::track::itsChi2NCl < cfgITSchi2Max;
  Filter tpcchi2filter = aod::track::tpcChi2NCl < cfgTPCchi2Max;
  Filter mapfilter = (aod::track::itsClusterMap & uint8_t(1)) > 0;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  std::unique_ptr<o2::phos::Geometry> geomPHOS;
  double bz{0.}; // magnetic field
  int runNumber{0};

  HistogramRegistry mHistManager{"PhosElIdHistograms"};
  TF1 *fSigma_dz, *fSigma_dx;
  float *PhosShiftX, *PhosShiftZ;

  TF1* fMeandXPosMod1;
  TF1* fMeandXNegMod1;
  TF1* fMeandXPosMod2;
  TF1* fMeandXNegMod2;
  TF1* fMeandXPosMod3;
  TF1* fMeandXNegMod3;
  TF1* fMeandXPosMod4;
  TF1* fMeandXNegMod4;

  void init(InitContext const&)
  {
    LOG(info) << "Initializing PHOS electron identification analysis task ...";

    std::vector<double> momentumBinning = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                                           1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                                           4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 10.};
    std::vector<float> parametersSigmadZ = pSigmadZ;
    std::vector<float> parametersSigmadX = pSigmadX;
    std::vector<float> vPhosShiftX = pPhosShiftX;
    std::vector<float> vPhosShiftZ = pPhosShiftZ;
    std::vector<float> meandXPosMod1 = pMeandXPosMod1;
    std::vector<float> meandXPosMod2 = pMeandXPosMod2;
    std::vector<float> meandXPosMod3 = pMeandXPosMod3;
    std::vector<float> meandXPosMod4 = pMeandXPosMod4;
    std::vector<float> meandXNegMod1 = pMeandXNegMod1;
    std::vector<float> meandXNegMod2 = pMeandXNegMod2;
    std::vector<float> meandXNegMod3 = pMeandXNegMod3;
    std::vector<float> meandXNegMod4 = pMeandXNegMod4;

    PhosShiftX = new float[4];
    PhosShiftX[0] = vPhosShiftX.at(0);
    PhosShiftX[1] = vPhosShiftX.at(1);
    PhosShiftX[2] = vPhosShiftX.at(2);
    PhosShiftX[3] = vPhosShiftX.at(3);

    PhosShiftZ = new float[4];
    PhosShiftZ[0] = vPhosShiftZ.at(0);
    PhosShiftZ[1] = vPhosShiftZ.at(1);
    PhosShiftZ[2] = vPhosShiftZ.at(2);
    PhosShiftZ[3] = vPhosShiftZ.at(3);

    const AxisSpec axisCounter{1, 0, +1, ""},
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
    mHistManager.add("TVXinPHOSCounter", "TVXinPHOSCounter", kTH1F, {axisCounter});

    mHistManager.add("hTrackPtEtaPhi", "Track pt vs eta vs phi", HistType::kTH3F, {axisPt, axisEta, axisPhi});
    mHistManager.add("hTrackPtEtaPhi_Phos", "Track pt vs eta vs phi on Phos surface", HistType::kTH3F, {axisPt, axisEta, axisPhi});
    mHistManager.add("hTrackDCA", "Track DCA info", HistType::kTH2F, {axisDCATrackXY, axisDCATrackZ});
    mHistManager.add("hTrackVX", "Track vertex coordinate X", HistType::kTH1F, {axisVTrackX});
    mHistManager.add("hTrackVY", "Track vertex coordinate Y", HistType::kTH1F, {axisVTrackY});
    mHistManager.add("hTrackVZ", "Track vertex coordinate Z", HistType::kTH1F, {axisVTrackZ});
    mHistManager.add("hColVX", "Collision vertex coordinate X", HistType::kTH1F, {axisVColX});
    mHistManager.add("hColVY", "Collision vertex coordinate Y", HistType::kTH1F, {axisVColY});
    mHistManager.add("hColVZ", "Collision vertex coordinate Z", HistType::kTH1F, {axisVColZ});
    mHistManager.add("hTrackPhosProjMod", "Track projection coordinates on PHOS modules", HistType::kTH3F, {axisX, axisZ, axisModes});

    mHistManager.add("hCluE_v_mod_v_time", "Cluster energy spectrum (E > 0.3 GeV) vs time per module", HistType::kTH3F, {axisE, axisTime, axisModes});

    mHistManager.add("hCluE_mod_energy_cut", "Cluster energy spectrum (E > 0.3 GeV) per module", HistType::kTH2F, {axisE, axisModes});
    mHistManager.add("hCluE_mod_time_cut", "Cluster energy spectrum (E > 0.3 GeV)(time +-25 ns) per module", HistType::kTH2F, {axisE, axisModes});
    mHistManager.add("hCluE_mod_cell_cut", "Cluster energy spectrum (E > 0.3 GeV)(time +-25 ns)(ncells > 3) per module", HistType::kTH2F, {axisE, axisModes});
    mHistManager.add("hCluE_mod_disp", "Cluster energy spectrum OK dispersion and (E > 0.3 GeV)(time +-25 ns)(ncells > 3) per module", HistType::kTH2F, {axisE, axisModes});

    mHistManager.add("hCluE_ncells_mod", "Cluster energy spectrum per module", HistType::kTH3F, {axisE, axisCells, axisModes});
    mHistManager.add("hCluXZ_mod", "Local cluster X Z per module", HistType::kTH3F, {axisX, axisZ, axisModes});

    mHistManager.add("hdZpmod", "dz,p_{tr},module", HistType::kTH3F, {axisdZ, axisPt, axisModes});
    mHistManager.add("hdZpmod_pos", "dz,p_{tr},module positive tracks", HistType::kTH3F, {axisdZ, axisPt, axisModes});
    mHistManager.add("hdZpmod_neg", "dz,p_{tr},module negative tracks", HistType::kTH3F, {axisdZ, axisPt, axisModes});
    mHistManager.add("hdXpmod", "dx,p_{tr},module", HistType::kTH3F, {axisdX, axisPt, axisModes});
    mHistManager.add("hdXpmod_pos", "dx,p_{tr},module positive tracks", HistType::kTH3F, {axisdX, axisPt, axisModes});
    mHistManager.add("hdXpmod_neg", "dx,p_{tr},module negative tracks", HistType::kTH3F, {axisdX, axisPt, axisModes});

    mHistManager.add("hCluE_v_pt_disp", "Cluster energy vs p | OK dispersion", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("hCluE_v_pt_Nsigma", "Cluster energy vs p within trackmatch Nsigma", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("hCluE_v_pt_Nsigma_disp", "Cluster energy vs p within trackmatch Nsigma | OK dispersion", HistType::kTH3F, {axisE, axisPt, axisModes});

    mHistManager.add("hCluE_v_pt_disp_TPC", "Cluster energy vs p | OK dispersion", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("hCluE_v_pt_Nsigma_TPC", "Cluster energy vs p within trackmatch Nsigma", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("hCluE_v_pt_Nsigma_disp_TPC", "Cluster energy vs p within trackmatch Nsigma | OK dispersion", HistType::kTH3F, {axisE, axisPt, axisModes});

    mHistManager.add("hEp_v_pt_disp", "E/p ratio vs p | OK dispersion", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("hEp_v_pt_Nsigma", "E/p ratio vs p within trackmatch Nsigma", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("hEp_v_pt_Nsigma_disp", "E/p ratio vs p within trackmatch Nsigma | OK dispersion", HistType::kTH3F, {axisEp, axisPt, axisModes});

    mHistManager.add("hEp_v_E_disp", "E/p ratio vs cluster E | OK dispersion", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("hEp_v_E_Nsigma", "E/p ratio vs cluster E within trackmatch Nsigma", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("hEp_v_E_Nsigma_disp", "E/p ratio vs cluster E within trackmatch Nsigma | OK dispersion", HistType::kTH3F, {axisEp, axisE, axisModes});

    mHistManager.add("hEp_v_pt_disp_TPC", "E/p ratio vs p | OK dispersion", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("hEp_v_pt_Nsigma_TPC", "E/p ratio vs p within trackmatch Nsigma", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("hEp_v_pt_Nsigma_disp_TPC", "E/p ratio vs p within trackmatch Nsigma | OK dispersion", HistType::kTH3F, {axisEp, axisPt, axisModes});

    mHistManager.add("hEp_v_E_disp_TPC", "E/p ratio vs cluster E | OK dispersion", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("hEp_v_E_Nsigma_TPC", "E/p ratio vs cluster E within trackmatch Nsigma", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("hEp_v_E_Nsigma_disp_TPC", "E/p ratio vs cluster E within trackmatch Nsigma | OK dispersion", HistType::kTH3F, {axisEp, axisE, axisModes});

    geomPHOS = std::make_unique<o2::phos::Geometry>("PHOS");
    fSigma_dz = new TF1("fSigma_dz", "[0]/(x+[1])^[2]+pol1(3)", 0.3, 10);
    fSigma_dz->SetParameters(parametersSigmadZ.at(0), parametersSigmadZ.at(1), parametersSigmadZ.at(2), parametersSigmadZ.at(3), parametersSigmadZ.at(4));

    fSigma_dx = new TF1("fSigma_dx", "[0]/x^[1]+[2]", 0.1, 10);
    fSigma_dx->SetParameters(parametersSigmadX.at(0), parametersSigmadX.at(1), parametersSigmadX.at(2));

    fMeandXPosMod1 = new TF1("funcMeandx_pos_mod1", "[0]/(x+[1])^[2]", 0.1, 10);
    fMeandXNegMod1 = new TF1("funcMeandx_neg_mod1", "[0]/(x+[1])^[2]", 0.1, 10);
    fMeandXPosMod2 = new TF1("funcMeandx_pos_mod2", "[0]/(x+[1])^[2]", 0.1, 10);
    fMeandXNegMod2 = new TF1("funcMeandx_neg_mod2", "[0]/(x+[1])^[2]", 0.1, 10);
    fMeandXPosMod3 = new TF1("funcMeandx_pos_mod3", "[0]/(x+[1])^[2]", 0.1, 10);
    fMeandXNegMod3 = new TF1("funcMeandx_neg_mod3", "[0]/(x+[1])^[2]", 0.1, 10);
    fMeandXPosMod4 = new TF1("funcMeandx_pos_mod4", "[0]/(x+[1])^[2]", 0.1, 10);
    fMeandXNegMod4 = new TF1("funcMeandx_neg_mod4", "[0]/(x+[1])^[2]", 0.1, 10);

    fMeandXPosMod1->SetParameters(meandXPosMod1.at(0), meandXPosMod1.at(1), meandXPosMod1.at(2));
    fMeandXPosMod2->SetParameters(meandXPosMod2.at(0), meandXPosMod2.at(1), meandXPosMod2.at(2));
    fMeandXPosMod3->SetParameters(meandXPosMod3.at(0), meandXPosMod3.at(1), meandXPosMod3.at(2));
    fMeandXPosMod4->SetParameters(meandXPosMod4.at(0), meandXPosMod4.at(1), meandXPosMod4.at(2));

    fMeandXNegMod1->SetParameters(meandXNegMod1.at(0), meandXNegMod1.at(1), meandXNegMod1.at(2));
    fMeandXNegMod2->SetParameters(meandXNegMod2.at(0), meandXNegMod2.at(1), meandXNegMod2.at(2));
    fMeandXNegMod3->SetParameters(meandXNegMod3.at(0), meandXNegMod3.at(1), meandXNegMod3.at(2));
    fMeandXNegMod4->SetParameters(meandXNegMod4.at(0), meandXNegMod4.at(1), meandXNegMod4.at(2));
  }
  void process(SelCollisions::iterator const& collision,
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
    if (std::fabs(collision.posZ()) > 10.f)
      return;
    mHistManager.fill(HIST("eventCounter"), 0.5);
    if (!collision.alias_bit(mEvSelTrig))
      return;
    mHistManager.fill(HIST("TVXinPHOSCounter"), 0.5);

    if (clusters.size() == 0)
      return; // Nothing to process
    mHistManager.fill(HIST("hColVX"), collision.posX());
    mHistManager.fill(HIST("hColVY"), collision.posY());
    mHistManager.fill(HIST("hColVZ"), collision.posZ());

    for (auto const& track : tracks) {

      if (!track.has_collision() || std::fabs(track.dcaXY()) > cfgDCAxyMax || std::fabs(track.dcaZ()) > cfgDCAzMax || !track.hasTPC() || std::fabs(track.eta()) > 0.15)
        continue;
      if (track.pt() < cfgPtMin || track.pt() > cfgPtMax)
        continue;
      if (track.itsChi2NCl() > cfgITSchi2Max)
        continue;
      if (track.itsNCls() < cfgITSnclsMin || track.itsNCls() > cfgITSnclsMax || !((track.itsClusterMap() & uint8_t(1)) > 0))
        continue;
      if (track.tpcChi2NCl() > cfgTPCchi2Max)
        continue;
      if (track.tpcNClsFound() < cfgTPCnclsMin || track.tpcNClsFound() > cfgTPCnclsMax)
        continue;
      if (track.tpcNClsCrossedRows() < cfgTPCnclsCRMin || track.tpcNClsCrossedRows() > cfgTPCnclsCRMax)
        continue;

      mHistManager.fill(HIST("hTrackVX"), track.x());
      mHistManager.fill(HIST("hTrackVY"), track.y());
      mHistManager.fill(HIST("hTrackVZ"), track.z());

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
        bool isTPCElectron = nsigmaTPCEl > cfgTPCNSigmaElMin && nsigmaTPCEl < cfgTPCNSigmaElMax;
        bool isTOFElectron = nsigmaTOFEl > cfgTOFNSigmaElMin && nsigmaTOFEl < cfgTOFNSigmaElMax;
        isElectron = isTPCElectron || isTOFElectron;

        float nsigmaTPCPi = track.tpcNSigmaPi();
        float nsigmaTPCKa = track.tpcNSigmaKa();
        float nsigmaTPCPr = track.tpcNSigmaPr();
        bool isPion = nsigmaTPCPi > cfgTPCNSigmaPiMin && nsigmaTPCPi < cfgTPCNSigmaPiMax;
        bool isKaon = nsigmaTPCKa > cfgTPCNSigmaKaMin && nsigmaTPCKa < cfgTPCNSigmaKaMax;
        bool isProton = nsigmaTPCPr > cfgTPCNSigmaPrMin && nsigmaTPCPr < cfgTPCNSigmaPrMax;
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

        bool isDispOK = testLambda(cluE, clu.m02(), clu.m20());
        float posX = clu.x(), posZ = clu.z(), dX = trackX - posX, dZ = trackZ - posZ, Ep = cluE / trackMom;

        mHistManager.fill(HIST("hdZpmod"), dZ, trackPT, module);
        mHistManager.fill(HIST("hdXpmod"), dX, trackPT, module);
        if (posTrack) {
          mHistManager.fill(HIST("hdZpmod_pos"), dZ, trackPT, module);
          mHistManager.fill(HIST("hdXpmod_pos"), dX, trackPT, module);
        } else {
          mHistManager.fill(HIST("hdZpmod_neg"), dZ, trackPT, module);
          mHistManager.fill(HIST("hdXpmod_neg"), dX, trackPT, module);
        }

        if (isDispOK) {
          mHistManager.fill(HIST("hCluE_v_pt_disp"), cluE, trackPT, module);
          mHistManager.fill(HIST("hEp_v_pt_disp"), Ep, trackPT, module);
          mHistManager.fill(HIST("hEp_v_E_disp"), Ep, cluE, module);
          if (isElectron) {
            mHistManager.fill(HIST("hCluE_v_pt_disp_TPC"), cluE, trackPT, module);
            mHistManager.fill(HIST("hEp_v_pt_disp_TPC"), Ep, trackPT, module);
            mHistManager.fill(HIST("hEp_v_E_disp_TPC"), Ep, cluE, module);
          }
        }
        if (!isWithinNSigma(module, trackPT, dZ, dX))
          continue;
        mHistManager.fill(HIST("hCluE_v_pt_Nsigma"), cluE, trackPT, module);
        mHistManager.fill(HIST("hEp_v_pt_Nsigma"), Ep, trackPT, module);
        mHistManager.fill(HIST("hEp_v_E_Nsigma"), Ep, cluE, module);
        if (isElectron) {
          mHistManager.fill(HIST("hCluE_v_pt_Nsigma_TPC"), cluE, trackPT, module);
          mHistManager.fill(HIST("hEp_v_pt_Nsigma_TPC"), Ep, trackPT, module);
          mHistManager.fill(HIST("hEp_v_E_Nsigma_TPC"), Ep, cluE, module);
        }
        if (isDispOK) {
          mHistManager.fill(HIST("hCluE_v_pt_Nsigma_disp"), cluE, trackPT, module);
          mHistManager.fill(HIST("hEp_v_pt_Nsigma_disp"), Ep, trackPT, module);
          mHistManager.fill(HIST("hEp_v_E_Nsigma_disp"), Ep, cluE, module);
          if (isElectron) {
            mHistManager.fill(HIST("hCluE_v_pt_Nsigma_disp_TPC"), cluE, trackPT, module);
            mHistManager.fill(HIST("hEp_v_pt_Nsigma_disp_TPC"), Ep, trackPT, module);
            mHistManager.fill(HIST("hEp_v_E_Nsigma_disp_TPC"), Ep, cluE, module);
          }
          phosMatch(collision.index(), clu.index(), track.index());
        }
      }

      mHistManager.fill(HIST("hTrackPtEtaPhi"), track.pt(), track.eta(), track.phi() * TMath::RadToDeg());
      mHistManager.fill(HIST("hTrackPtEtaPhi_Phos"), track.pt(), track.trackEtaEmcal(), track.trackPhiEmcal() * TMath::RadToDeg());
      mHistManager.fill(HIST("hTrackDCA"), track.dcaXY(), track.dcaZ());
      mHistManager.fill(HIST("hTrackPhosProjMod"), trackX, trackZ, module);
    } // end of double loop

    for (auto const& clu : clusters) {
      double cluE = clu.e(), cluTime = clu.time();
      int mod = clu.mod();
      if (cluE > mMinCluE) {
        mHistManager.fill(HIST("hCluE_mod_energy_cut"), cluE, mod);
        mHistManager.fill(HIST("hCluE_v_mod_v_time"), cluE, cluTime * 1e9, mod);
        if (cluTime < mMaxCluTime && cluTime > mMinCluTime) {
          mHistManager.fill(HIST("hCluE_mod_time_cut"), cluE, mod);
          if (clu.ncell() >= mMinCluNcell) {
            mHistManager.fill(HIST("hCluE_mod_cell_cut"), cluE, mod);
            mHistManager.fill(HIST("hCluXZ_mod"), clu.x(), clu.z(), mod);
            mHistManager.fill(HIST("hCluE_ncells_mod"), cluE, clu.ncell(), mod);
            if (testLambda(cluE, clu.m02(), clu.m20()))
              mHistManager.fill(HIST("hCluE_mod_disp"), cluE, mod);
          }
        }
      }
    } // end of cluster loop
  }

  bool isWithinNSigma(int16_t& mod, float p, float deltaZ, float deltaX)
  {
    if (mod == 1) {
      if (std::fabs(deltaZ - PhosShiftZ[0]) > cfgNsigmaTrackMatch * fSigma_dz->Eval(p))
        return false;
      if (std::fabs(deltaX - fMeandXPosMod1->Eval(p) + PhosShiftX[0]) > cfgNsigmaTrackMatch * fSigma_dx->Eval(p))
        return false;
    } else if (mod == 2) {
      if (std::fabs(deltaZ - PhosShiftZ[1]) > cfgNsigmaTrackMatch * fSigma_dz->Eval(p))
        return false;
      if (std::fabs(deltaX - fMeandXPosMod2->Eval(p) + PhosShiftX[1]) > cfgNsigmaTrackMatch * fSigma_dx->Eval(p))
        return false;
    } else if (mod == 3) {
      if (std::fabs(deltaZ - PhosShiftZ[2]) > cfgNsigmaTrackMatch * fSigma_dz->Eval(p))
        return false;
      if (std::fabs(deltaX - fMeandXPosMod3->Eval(p) + PhosShiftX[2]) > cfgNsigmaTrackMatch * fSigma_dx->Eval(p))
        return false;
    } else if (mod == 4) {
      if (std::fabs(deltaZ - PhosShiftZ[3]) > cfgNsigmaTrackMatch * fSigma_dz->Eval(p))
        return false;
      if (std::fabs(deltaX - fMeandXPosMod4->Eval(p) + PhosShiftX[3]) > cfgNsigmaTrackMatch * fSigma_dx->Eval(p))
        return false;
    }
    return true;
  }

  ///////////////////////////////////////////////////////////////////////
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
    module = 1 + static_cast<int16_t>((trackPhi - phiMin) / dphi);
    if (module < 1) {
      module = 1;
    }
    if (module > 4) {
      module = 4;
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
  //_____________________________________________________________________________
  bool testLambda(float pt, float l1, float l2)
  {
    // Parameterization for full dispersion
    float l2Mean = 1.53126 + 9.50835e+06 / (1. + 1.08728e+07 * pt + 1.73420e+06 * pt * pt);
    float l1Mean = 1.12365 + 0.123770 * std::exp(-pt * 0.246551) + 5.30000e-03 * pt;
    float l2Sigma = 6.48260e-02 + 7.60261e+10 / (1. + 1.53012e+11 * pt + 5.01265e+05 * pt * pt) + 9.00000e-03 * pt;
    float l1Sigma = 4.44719e-04 + 6.99839e-01 / (1. + 1.22497e+00 * pt + 6.78604e-07 * pt * pt) + 9.00000e-03 * pt;
    float c = -0.35 - 0.550 * std::exp(-0.390730 * pt);

    return 0.5 * (l1 - l1Mean) * (l1 - l1Mean) / l1Sigma / l1Sigma +
             0.5 * (l2 - l2Mean) * (l2 - l2Mean) / l2Sigma / l2Sigma +
             0.5 * c * (l1 - l1Mean) * (l2 - l2Mean) / l1Sigma / l2Sigma <
           4.;
  }
};

struct MassSpectra {

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults,
                                  aod::FT0sCorrected, aod::CentFT0Ms,
                                  aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As,
                                  aod::CentFDDMs, aod::CentNTPVs>;
  using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                             aod::TracksDCACov, aod::pidTOFFullEl, aod::pidTPCFullEl,
                             aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
  Configurable<int> mEvSelTrig{"mEvSelTrig", kTVXinPHOS, "Select events with this trigger"},
    cfgMassBinning{"cfgMassBinning", 1000, "Binning for mass"},
    cfgEnergyBinning{"cfgEnergyBinning", 100, "Binning for energy"},
    cfgEpRatioBinning{"cfgEpRatioBinning", 200, "Binning for energy to momentum ratio"},
    cfgCentBinning{"cfgCentBinning", 10, "Binning for centrality"},
    cfgCentEst{"cfgCentEst", 1, "Centrality estimator, 0: FV0A, 1: FT0M, 2: FT0A, 3: FT0C, 4: FDDM, 5: NTPV"};

  Configurable<float> cfgEtaMax{"cfgEtaMax", {0.8f}, "eta ranges"},
    cfgPtMin{"cfgPtMin", {0.2f}, "pt min"},
    cfgPtMax{"cfgPtMax", {20.f}, "pt max"},
    cfgMassSpectraMin{"cfgMassSpectraMin", {2.5f}, "mass spectra min for e+e-"},
    cfgMassSpectraMax{"cfgMassSpectraMax", {3.5f}, "mass spcetra max for e+e-"},
    cfgDCAxyMax{"cfgDCAxyMax", {3.f}, "dcaxy max"},
    cfgDCAzMax{"cfgDCAzMax", {3.f}, "dcaz max"},
    cfgITSchi2Max{"cfgITSchi2Max", {5.f}, "its chi2 max"},
    cfgITSnclsMin{"cfgITSnclsMin", {4.5f}, "min number of ITS clusters"},
    cfgITSnclsMax{"cfgITSnclsMax", {7.5f}, "max number of ITS clusters"},
    cfgTPCchi2Max{"cfgTPCchi2Max", {4.f}, "tpc chi2 max"},
    cfgTPCnclsMin{"cfgTPCnclsMin", {90.f}, "min number of TPC clusters"},
    cfgTPCnclsMax{"cfgTPCnclsMax", {170.f}, "max number of TPC clusters"},
    cfgTPCnclsCRMin{"cfgTPCnclsCRMin", {80.f}, "min number of TPC crossed rows"},
    cfgTPCnclsCRMax{"cfgTPCnclsCRMax", {161.f}, "max number of TPC crossed rows"},
    cfgTPCNSigmaElMin{"cfgTPCNSigmaElMin", {-3.f}, "min TPC nsigma e for inclusion"},
    cfgTPCNSigmaElMax{"cfgTPCNSigmaElMax", {2.f}, "max TPC nsigma e for inclusion"},
    cfgTPCNSigmaPiMin{"cfgTPCNSigmaPiMin", {-3.f}, "min TPC nsigma pion for exclusion"},
    cfgTPCNSigmaPiMax{"cfgTPCNSigmaPiMax", {3.5f}, "max TPC nsigma pion for exclusion"},
    cfgTPCNSigmaPrMin{"cfgTPCNSigmaPrMin", {-3.f}, "min TPC nsigma proton for exclusion"},
    cfgTPCNSigmaPrMax{"cfgTPCNSigmaPrMax", {4.f}, "max TPC nsigma proton for exclusion"},
    cfgTPCNSigmaKaMin{"cfgTPCNSigmaKaMin", {-3.f}, "min TPC nsigma kaon for exclusion"},
    cfgTPCNSigmaKaMax{"cfgTPCNSigmaKaMax", {4.f}, "max TPC nsigma kaon for exclusion"},
    cfgTOFNSigmaElMin{"cfgTOFNSigmaElMin", {-3.f}, "min TOF nsigma e for inclusion"},
    cfgTOFNSigmaElMax{"cfgTOFNSigmaElMax", {3.f}, "max TOF nsigma e for inclusion"},
    cfgShiftEp{"cfgShiftEp", {0.055f}, "PHOS E/p shift for electrons"},
    cfgNsigmaEp{"cfgNsigmaEp", {2.f}, "PHOS E/p nsigma for inclusion"};

  Configurable<std::vector<float>> cfgEpSigmaPars{"cfgEpSigmaPars", {1.3e-02, 1.9e-02, 1.1e-02, 3.e-02}, "E/p sigma function parameters (from alice 3 mc tests + const)"};

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
    const AxisSpec axisCounter{1, 0, +1, ""},
      axisCent{cfgCentBinning, 0, 100, "centrality percentage"},
      axisPt{momentumBinning, "p_{T} (GeV/c)"},
      axisEp{cfgEpRatioBinning, 0., 2., "E/p", "E_{cluster}/p_{track}"},
      axisE{cfgEnergyBinning, 0, 10, "E (GeV)", "E (GeV)"},
      axisMassSpectrum{cfgMassBinning, cfgMassSpectraMin, cfgMassSpectraMax, "M (GeV/c^{2})", "Mass e^{+}e^{-} (GeV/c^{2})"};

    mHistManager.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    mHistManager.add("TVXinPHOSCounter", "TVXinPHOSCounter", kTH1F, {axisCounter});

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

    std::vector<float> epSigmaPars = cfgEpSigmaPars;
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
    if (std::fabs(collision.posZ()) > 10.f)
      return;
    mHistManager.fill(HIST("eventCounter"), 0.5);
    if (!collision.alias_bit(mEvSelTrig))
      return;
    mHistManager.fill(HIST("TVXinPHOSCounter"), 0.5);

    if (clusters.size() == 0)
      return; // Nothing to process

    float cent = -1.;
    switch (cfgCentEst) {
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

      if (!TPCel.has_collision() || std::fabs(TPCel.dcaXY()) > cfgDCAxyMax || std::fabs(TPCel.dcaZ()) > cfgDCAzMax || !TPCel.hasTPC() || std::fabs(TPCel.eta()) > 0.15)
        continue;
      if (TPCel.pt() < cfgPtMin || TPCel.pt() > cfgPtMax)
        continue;
      if (TPCel.itsChi2NCl() > cfgITSchi2Max)
        continue;
      if (TPCel.itsNCls() < cfgITSnclsMin || TPCel.itsNCls() > cfgITSnclsMax || !((TPCel.itsClusterMap() & uint8_t(1)) > 0))
        continue;
      if (TPCel.tpcChi2NCl() > cfgTPCchi2Max)
        continue;
      if (TPCel.tpcNClsFound() < cfgTPCnclsMin || TPCel.tpcNClsFound() > cfgTPCnclsMax)
        continue;
      if (TPCel.tpcNClsCrossedRows() < cfgTPCnclsCRMin || TPCel.tpcNClsCrossedRows() > cfgTPCnclsCRMax)
        continue;

      bool isElectron = false;
      float nsigmaTPCEl = TPCel.tpcNSigmaEl();
      float nsigmaTOFEl = TPCel.tofNSigmaEl();
      bool isTPCElectron = nsigmaTPCEl > cfgTPCNSigmaElMin && nsigmaTPCEl < cfgTPCNSigmaElMax;
      bool isTOFElectron = nsigmaTOFEl > cfgTOFNSigmaElMin && nsigmaTOFEl < cfgTOFNSigmaElMax;
      isElectron = isTPCElectron || isTOFElectron;

      float nsigmaTPCPi = TPCel.tpcNSigmaPi();
      float nsigmaTPCKa = TPCel.tpcNSigmaKa();
      float nsigmaTPCPr = TPCel.tpcNSigmaPr();
      bool isPion = nsigmaTPCPi > cfgTPCNSigmaPiMin && nsigmaTPCPi < cfgTPCNSigmaPiMax;
      bool isKaon = nsigmaTPCKa > cfgTPCNSigmaKaMin && nsigmaTPCKa < cfgTPCNSigmaKaMax;
      bool isProton = nsigmaTPCPr > cfgTPCNSigmaPrMin && nsigmaTPCPr < cfgTPCNSigmaPrMax;
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

        float mass2Tracks = 0, mom2Tracks = 0, cluE = clust2.e();
        TLorentzVector fourVectorP1, fourVectorP2;
        fourVectorP1.SetPxPyPzE(TPCel.px(), TPCel.py(), TPCel.pz(), TPCel.energy(0));
        fourVectorP2.SetPxPyPzE(track2.px(), track2.py(), track2.pz(), track2.energy(0));
        mom2Tracks = (fourVectorP1 + fourVectorP2).Pt();
        mass2Tracks = (fourVectorP1 + fourVectorP2).M();
        bool elCandidate = (std::fabs(cluE / track2.p() - cfgShiftEp - 1) < cfgNsigmaEp * fEpSigmaPhos->Eval(cluE));

        if (TPCel.sign() == track2.sign()) {
          if (posTrack) {
            mHistManager.fill(HIST("h_eh_pp_mass_spectra_v_pt_v_cent"), mass2Tracks, mom2Tracks, cent);
            mHistManager.fill(HIST("h_eh_pp_mass_spectra_v_E_v_cent"), mass2Tracks, cluE, cent);
            if (elCandidate) {
              mHistManager.fill(HIST("h_ee_pp_mass_spectra_v_pt_v_cent"), mass2Tracks, mom2Tracks, cent);
              mHistManager.fill(HIST("h_ee_pp_mass_spectra_v_E_v_cent"), mass2Tracks, cluE, cent);
            }
          } else {
            mHistManager.fill(HIST("h_eh_mm_mass_spectra_v_pt_v_cent"), mass2Tracks, mom2Tracks, cent);
            mHistManager.fill(HIST("h_eh_mm_mass_spectra_v_E_v_cent"), mass2Tracks, cluE, cent);
            if (elCandidate) {
              mHistManager.fill(HIST("h_ee_mm_mass_spectra_v_pt_v_cent"), mass2Tracks, mom2Tracks, cent);
              mHistManager.fill(HIST("h_ee_mm_mass_spectra_v_E_v_cent"), mass2Tracks, cluE, cent);
            }
          }
        } else {
          mHistManager.fill(HIST("h_eh_mp_mass_spectra_v_pt_v_cent"), mass2Tracks, mom2Tracks, cent);
          mHistManager.fill(HIST("h_eh_mp_mass_spectra_v_E_v_cent"), mass2Tracks, cluE, cent);
          if (elCandidate) {
            mHistManager.fill(HIST("h_ee_mp_mass_spectra_v_pt_v_cent"), mass2Tracks, mom2Tracks, cent);
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
      bool elCandidate = (std::fabs(epRatio - cfgShiftEp - 1) < cfgNsigmaEp * fEpSigmaPhos->Eval(cluE));
      if (elCandidate)
        mHistManager.fill(HIST("hEp_v_E_v_cent_cutEp"), epRatio, cluE, cent);
    }
  }
};

struct TpcElIdMassSpectrum {

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::Mults,
                                  aod::FT0sCorrected, aod::CentFT0Ms,
                                  aod::CentFT0As, aod::CentFT0Cs, aod::CentFV0As,
                                  aod::CentFDDMs, aod::CentNTPVs>;
  using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                             aod::TracksDCACov, aod::pidTOFFullEl, aod::pidTPCFullEl,
                             aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;
  Configurable<float> mMinCluE{"mMinCluE", 0.1, "Minimum cluster energy for photons in the analysis"},
    mCutMIPCluE{"mCutMIPCluE", 0.3, "Min cluster energy to reject MIPs in the analysis"},
    mMaxCluE{"mMaxCluE", 1., "Maximum cluster energy for photons in the analysis"},
    mMinCluTime{"minCluTime", -25.e-9, "Min. cluster time"},
    mMaxCluTime{"mMaxCluTime", 25.e-9, "Max. cluster time"},
    cfgEtaMax{"cfgEtaMax", {0.8f}, "eta ranges"},
    cfgPtMin{"cfgPtMin", {0.2f}, "pt min"},
    cfgPtMax{"cfgPtMax", {20.f}, "pt max"},
    cfgMassSpectraJpsiMin{"cfgMassSpectraJpsiMin", {2.5f}, "mass spectra min for Jpsi region"},
    cfgMassSpectraJpsiMax{"cfgMassSpectraJpsiMax", {3.5f}, "mass spcetra max for Jpsi region"},
    cfgMassSpectraChicMin{"cfgMassSpectraChicMin", {3.f}, "mass spectra min Chic region"},
    cfgMassSpectraChicMax{"cfgMassSpectraChicMax", {4.f}, "mass spcetra max Chic region"},
    cfgDCAxyMax{"cfgDCAxyMax", {3.f}, "dcaxy max"},
    cfgDCAzMax{"cfgDCAzMax", {3.f}, "dcaz max"},
    cfgITSchi2Max{"cfgITSchi2Max", {5.f}, "its chi2 max"},
    cfgITSnclsMin{"cfgITSnclsMin", {4.5f}, "min number of ITS clusters"},
    cfgITSnclsMax{"cfgITSnclsMax", {7.5f}, "max number of ITS clusters"},
    cfgTPCchi2Max{"cfgTPCchi2Max", {4.f}, "tpc chi2 max"},
    cfgTPCnclsMin{"cfgTPCnclsMin", {90.f}, "min number of TPC clusters"},
    cfgTPCnclsMax{"cfgTPCnclsMax", {170.f}, "max number of TPC clusters"},
    cfgTPCnclsCRMin{"cfgTPCnclsCRMin", {80.f}, "min number of TPC crossed rows"},
    cfgTPCnclsCRMax{"cfgTPCnclsCRMax", {161.f}, "max number of TPC crossed rows"},
    cfgTPCNSigmaElMin{"cfgTPCNSigmaElMin", {-3.f}, "min TPC nsigma e for inclusion"},
    cfgTPCNSigmaElMax{"cfgTPCNSigmaElMax", {2.f}, "max TPC nsigma e for inclusion"},
    cfgTPCNSigmaPiMin{"cfgTPCNSigmaPiMin", {-3.f}, "min TPC nsigma pion for exclusion"},
    cfgTPCNSigmaPiMax{"cfgTPCNSigmaPiMax", {3.5f}, "max TPC nsigma pion for exclusion"},
    cfgTPCNSigmaPrMin{"cfgTPCNSigmaPrMin", {-3.f}, "min TPC nsigma proton for exclusion"},
    cfgTPCNSigmaPrMax{"cfgTPCNSigmaPrMax", {4.f}, "max TPC nsigma proton for exclusion"},
    cfgTPCNSigmaKaMin{"cfgTPCNSigmaKaMin", {-3.f}, "min TPC nsigma kaon for exclusion"},
    cfgTPCNSigmaKaMax{"cfgTPCNSigmaKaMax", {4.f}, "max TPC nsigma kaon for exclusion"},
    cfgTOFNSigmaElMin{"cfgTOFNSigmaElMin", {-3.f}, "min TOF nsigma e for inclusion"},
    cfgTOFNSigmaElMax{"cfgTOFNSigmaElMax", {3.f}, "max TOF nsigma e for inclusion"},
    cfgPhosRangeEta{"cfgPhosRangeEta", {0.12f}, "Phos range definition plus minus eta"},
    cfgPhosRangePhiMin{"cfgPhosRangePhiMin", {230.f}, "Phos range angle phi min"},
    cfgPhosRangePhiMax{"cfgPhosRangePhiMax", {330.f}, "Phos range angle phi max"},
    cfgeeMassMin{"cfgeeMassMin", {2.9f}, "J/psi(e+e-) Mass corridor lower limit"},
    cfgeeMassMax{"cfgeeMassMax", {3.3f}, "J/psi(e+e-) Mass corridor upper limit"},
    cfgJpsiMass{"cfgJpsiMass", {3.097f}, "J/psi Mass constant"};

  Configurable<int> mEvSelTrig{"mEvSelTrig", kTVXinPHOS, "Select events with this trigger"},
    cfgCentBinning{"cfgCentBinning", 10, "Binning for centrality"},
    cfgCentEst{"cfgCentEst", 1, "Centrality estimator, 0: FV0A, 1: FT0M, 2: FT0A, 3: FT0C, 4: FDDM, 5: NTPV"},
    cfgMassBinning{"cfgMassBinning", 1000, "Binning for mass"},
    cfgEnergyBinning{"cfgEnergyBinning", 100, "Binning for energy"},
    mMinCluNcell{"minCluNcell", 3, "min cells in cluster"};

  Filter ptFilter = (aod::track::pt > cfgPtMin) && (aod::track::pt < cfgPtMax),
         etaFilter = nabs(aod::track::eta) < cfgEtaMax,
         dcaxyFilter = nabs(aod::track::dcaXY) < cfgDCAxyMax,
         dcazFilter = nabs(aod::track::dcaZ) < cfgDCAzMax;

  Filter tpctofEl = ((aod::pidtpc::tpcNSigmaEl > cfgTPCNSigmaElMin) && (aod::pidtpc::tpcNSigmaEl < cfgTPCNSigmaElMax)) ||
                    ((aod::pidtof::tofNSigmaEl > cfgTOFNSigmaElMin) && (aod::pidtof::tofNSigmaEl < cfgTOFNSigmaElMax)),
         tpcPiRej = (aod::pidtpc::tpcNSigmaPi < cfgTPCNSigmaPiMin) || (aod::pidtpc::tpcNSigmaPi > cfgTPCNSigmaPiMax),
         tpcKaRej = (aod::pidtpc::tpcNSigmaKa < cfgTPCNSigmaKaMin) || (aod::pidtpc::tpcNSigmaKa > cfgTPCNSigmaPrMax),
         tpcPrRej = (aod::pidtpc::tpcNSigmaPr < cfgTPCNSigmaPrMin) || (aod::pidtpc::tpcNSigmaPr > cfgTPCNSigmaPrMax);

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  double bz{0.}; // magnetic field
  int runNumber{0};

  HistogramRegistry mHistManager{"tpcElIdHistograms"};

  void init(InitContext const&)
  {
    LOG(info) << "Initializing ee mass spectrum via TPC electron identification analysis task ...";

    std::vector<double> momentumBinning = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                                           1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                                           4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 10.};
    const AxisSpec axisCounter{1, 0, +1, ""},
      axisCent{cfgCentBinning, 0, 100, "centrality percentage"},
      axisVTrackX{400, -5., 5., "track vertex x (cm)", "track vertex x (cm)"},
      axisVTrackY{400, -5., 5., "track vertex y (cm)", "track vertex y (cm)"},
      axisVTrackZ{400, -20., 20., "track vertex z (cm)", "track vertex z (cm)"},
      axisE{cfgEnergyBinning, 0, 10, "E (GeV)", "E (GeV)"},
      axisMassSpectrum{cfgMassBinning, cfgMassSpectraJpsiMin, cfgMassSpectraJpsiMax, "M (GeV/c^{2})", "Mass e^{+}e^{-} (GeV/c^{2})"},
      axisMassSpectrumChiC{cfgMassBinning, cfgMassSpectraChicMin, cfgMassSpectraChicMax, "M (GeV/c^{2})", "Mass e^{+}e^{-}#gamma (GeV/c^{2})"},
      axisMassSpectrumChiCNoJpsiErrors{cfgMassBinning, cfgMassSpectraChicMin, cfgMassSpectraChicMax, "M (GeV/c^{2})", "Mass e^{+}e^{-}#gamma - Mass e^{+}e^{-} + Mass J/#psi (GeV/c^{2})"},
      axisTPC{1000, 0, 200, "TPC signal (dE/dx)"},
      axisPt{momentumBinning, "p_{T} (GeV/c)"},
      axisPtBig{2000, 0, 20, "p_{T} (GeV/c)"},
      axisEta{600, -3., 3., "#eta"};

    mHistManager.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    mHistManager.add("hTPCspectra", "TPC dE/dx spectra", HistType::kTH2F, {axisPt, axisTPC});
    mHistManager.add("hTPCspectra_isElectronRej", "isElectron with rejection | TPC dE/dx spectra", HistType::kTH2F, {axisPt, axisTPC});

    mHistManager.add("h_TPCee_MS_mp_v_pt_v_cent", "Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp} (from TPC candidates)", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("h_TPCee_MS_mm_v_pt_v_cent", "Mass e^{-}e^{-} vs momentum e^{-}e^{-} (from TPC candidates)", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("h_TPCee_MS_pp_v_pt_v_cent", "Mass e^{+}e^{+} vs momentum e^{+}e^{+} (from TPC candidates)", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});

    mHistManager.add("h_TPCee_MS_mp_phosRange_v_pt_v_cent", "Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp} (from TPC candidates) with one e in phos acceptance range", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("h_TPCee_MS_mm_phosRange_v_pt_v_cent", "Mass e^{-}e^{-} vs momentum e^{-}e^{-} (from TPC candidates) with one e in phos acceptance range", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("h_TPCee_MS_pp_phosRange_v_pt_v_cent", "Mass e^{+}e^{+} vs momentum e^{+}e^{+} (from TPC candidates) with one e in phos acceptance range", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});

    mHistManager.add("h_TPCee_MS_mp_phosRange_kTVXinPHOS_v_pt_v_cent", "Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp} (from TPC candidates) with one e in phos acceptance range", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("h_TPCee_MS_mm_phosRange_kTVXinPHOS_v_pt_v_cent", "Mass e^{-}e^{-} vs momentum e^{-}e^{-} (from TPC candidates) with one e in phos acceptance range", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});
    mHistManager.add("h_TPCee_MS_pp_phosRange_kTVXinPHOS_v_pt_v_cent", "Mass e^{+}e^{+} vs momentum e^{+}e^{+} (from TPC candidates) with one e in phos acceptance range", HistType::kTH3F, {axisMassSpectrum, axisPt, axisCent});

    mHistManager.add("h_TPCeePhosGamma_MS_withMatches_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}Track vs momentum e^{#pm}e^{#mp}Track", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_MS_noMatches_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_MS_noMatches_aroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_MS_noMatches_aroundJpsi_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | DispOK", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_MS_noMatches_aroundJpsi_DispNotOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | DispNotOK", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_MS_noMatches_notAroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma (not around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_MS_noMatches_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma | DispOK", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_MS_noMatches_DispNotOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma | DispNotOK", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_MS_aroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_MS_notAroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma (not around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_MS_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma | DispOK", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_MS_DispNotOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma | DispNotOK", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_MS_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_MS_withMatches_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}Track vs momentum e^{#pm}e^{#mp}Track | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_MS_noMatches_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_MS_noMatches_aroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_MS_noMatches_aroundJpsi_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | DispOK | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_MS_noMatches_aroundJpsi_DispNotOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | DispNotOK | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_MS_noMatches_notAroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma (not around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_MS_noMatches_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma | DispOK | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_MS_noMatches_DispNotOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma | DispNotOK | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_MS_aroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_MS_notAroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma (not around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_MS_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma | DispOK | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_MS_DispNotOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma | DispNotOK | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_isMIP_MS_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma | cluE < E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});

    mHistManager.add("h_TPCeePhosGamma_MS_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma vs momentum e^{#pm}e^{#mp}#gamma (TPC candidates + Phos cluster)", HistType::kTH3F, {axisMassSpectrumChiC, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_MS_v_cluE", "Mass e^{#pm}e^{#mp}#gamma vs cluster Energy left by the photon", HistType::kTH3F, {axisMassSpectrumChiC, axisE, axisCent});

    mHistManager.add("h_TPCeePhosGamma_minusee_MS_withMatches_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}Track - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}Track", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_minusee_MS_noMatches_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}#gamma", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_minusee_MS_noMatches_aroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_minusee_MS_noMatches_aroundJpsi_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | DispOK", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_minusee_MS_noMatches_aroundJpsi_DispNotOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | DispNotOK", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_minusee_MS_noMatches_notAroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp}  (not around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_minusee_MS_noMatches_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}#gamma | DispOK", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_minusee_MS_noMatches_DispNotOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}#gamma | DispNotOK", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_minusee_MS_aroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_minusee_MS_notAroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp}  (not around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_minusee_MS_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}#gamma | DispOK", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_minusee_MS_DispNotOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}#gamma | DispNotOK", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_minusee_MS_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}#gamma | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_minusee_MS_withMatches_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}Track - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}Track | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_minusee_MS_noMatches_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}#gamma | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_minusee_MS_noMatches_aroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_minusee_MS_noMatches_aroundJpsi_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | DispOK | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_minusee_MS_noMatches_aroundJpsi_DispNotOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | DispNotOK | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_minusee_MS_noMatches_notAroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp}  (not around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_minusee_MS_noMatches_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}#gamma | DispOK | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_minusee_MS_noMatches_DispNotOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}#gamma | DispNotOK | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_minusee_MS_aroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} (around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_minusee_MS_notAroundJpsi_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp}  (not around J/#psi) vs momentum e^{#pm}e^{#mp}#gamma | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_minusee_MS_DispOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}#gamma | DispOK | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_notMIP_minusee_MS_DispNotOK_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}#gamma | DispNotOK | cluE > E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_isMIP_minusee_MS_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}#gamma | cluE < E_{MIP}", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});

    mHistManager.add("h_TPCeePhosGamma_minusee_MS_v_3pt_v_cent", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp}  - Mass e^{#pm}e^{#mp} + Mass J/#psi vs momentum e^{#pm}e^{#mp}#gamma (TPC candidates + Phos cluster)", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisPt, axisCent});
    mHistManager.add("h_TPCeePhosGamma_minusee_MS_v_cluE", "Mass e^{#pm}e^{#mp}#gamma - Mass e^{#pm}e^{#mp}  - Mass e^{#pm}e^{#mp} + Mass J/#psi vs cluster Energy left by the photon", HistType::kTH3F, {axisMassSpectrumChiCNoJpsiErrors, axisE, axisCent});

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
  }
  void process(SelCollisions::iterator const& collision,
               aod::CaloClusters const& clusters,
               MyTracks const& tracks,
               soa::Filtered<MyTracks> const& filteredTracks,
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
    mHistManager.fill(HIST("eventCounter"), 0.5);
    if (std::fabs(collision.posZ()) > 10.f)
      return;

    float cent = -1.;
    switch (cfgCentEst) {
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

    for (auto const& [track1, track2] : combinations(CombinationsStrictlyUpperIndexPolicy(filteredTracks, filteredTracks))) {
      if (!track1.has_collision() || !track1.hasTPC())
        continue;
      if (!track2.has_collision() || !track2.hasTPC())
        continue;
      if (track1.collisionId() != track2.collisionId())
        continue;
      if (!((track1.itsClusterMap() & uint8_t(1)) > 0) || !((track2.itsClusterMap() & uint8_t(1)) > 0))
        continue;
      if (track1.itsChi2NCl() > cfgITSchi2Max || track2.itsChi2NCl() > cfgITSchi2Max)
        continue;
      if (track1.tpcChi2NCl() > cfgTPCchi2Max || track2.tpcChi2NCl() > cfgTPCchi2Max)
        continue;
      if (track1.itsNCls() < cfgITSnclsMin || track2.itsNCls() < cfgITSnclsMin)
        continue;
      if (track1.itsNCls() > cfgITSnclsMax || track2.itsNCls() > cfgITSnclsMax)
        continue;
      if (track1.tpcNClsFound() < cfgTPCnclsMin || track2.tpcNClsFound() < cfgTPCnclsMin)
        continue;
      if (track1.tpcNClsFound() > cfgTPCnclsMax || track2.tpcNClsFound() > cfgTPCnclsMax)
        continue;
      if (track1.tpcNClsCrossedRows() < cfgTPCnclsCRMin || track2.tpcNClsCrossedRows() < cfgTPCnclsCRMin)
        continue;
      if (track1.tpcNClsCrossedRows() > cfgTPCnclsCRMax || track2.tpcNClsCrossedRows() > cfgTPCnclsCRMax)
        continue;

      TLorentzVector fourVectorP1, fourVectorP2;
      fourVectorP1.SetPxPyPzE(track1.px(), track1.py(), track1.pz(), track1.energy(0));
      fourVectorP2.SetPxPyPzE(track2.px(), track2.py(), track2.pz(), track2.energy(0));

      bool inPhosEtaRange1 = std::fabs(track1.eta()) < cfgPhosRangeEta;
      bool inPhosEtaRange2 = std::fabs(track2.eta()) < cfgPhosRangeEta;
      bool inPhosPhiRange1 = (track1.phi() * TMath::RadToDeg() > cfgPhosRangePhiMin && track1.phi() * TMath::RadToDeg() < cfgPhosRangePhiMax);
      bool inPhosPhiRange2 = (track2.phi() * TMath::RadToDeg() > cfgPhosRangePhiMin && track2.phi() * TMath::RadToDeg() < cfgPhosRangePhiMax);
      bool inPhosRange = (inPhosEtaRange1 && inPhosPhiRange1) || (inPhosEtaRange2 && inPhosPhiRange2);
      bool posTrack = track1.sign() * bz > 0;

      double pairMass = (fourVectorP1 + fourVectorP2).M(), pairPt = (fourVectorP1 + fourVectorP2).Pt();

      if (track1.sign() == track2.sign()) {
        if (posTrack) {
          mHistManager.fill(HIST("h_TPCee_MS_pp_v_pt_v_cent"), pairMass, pairPt, cent);
          if (inPhosRange) {
            mHistManager.fill(HIST("h_TPCee_MS_pp_phosRange_v_pt_v_cent"), pairMass, pairPt, cent);
            if (collision.alias_bit(mEvSelTrig))
              mHistManager.fill(HIST("h_TPCee_MS_pp_phosRange_kTVXinPHOS_v_pt_v_cent"), pairMass, pairPt, cent);
          }
        } else {
          mHistManager.fill(HIST("h_TPCee_MS_mm_v_pt_v_cent"), pairMass, pairPt, cent);
          if (inPhosRange) {
            mHistManager.fill(HIST("h_TPCee_MS_mm_phosRange_v_pt_v_cent"), pairMass, pairPt, cent);
            if (collision.alias_bit(mEvSelTrig))
              mHistManager.fill(HIST("h_TPCee_MS_mm_phosRange_kTVXinPHOS_v_pt_v_cent"), pairMass, pairPt, cent);
          }
        }
      } else {
        mHistManager.fill(HIST("h_TPCee_MS_mp_v_pt_v_cent"), pairMass, pairPt, cent);
        if (inPhosRange) {
          mHistManager.fill(HIST("h_TPCee_MS_mp_phosRange_v_pt_v_cent"), pairMass, pairPt, cent);
          if (collision.alias_bit(mEvSelTrig))
            mHistManager.fill(HIST("h_TPCee_MS_mp_phosRange_kTVXinPHOS_v_pt_v_cent"), pairMass, pairPt, cent);
        }

        if (collision.alias_bit(mEvSelTrig) && clusters.size() != 0) {
          for (auto const& gamma : clusters) {
            float cluE = gamma.e();

            if (cluE < mMinCluE || cluE > mMaxCluE ||
                gamma.ncell() < mMinCluNcell ||
                gamma.time() > mMaxCluTime || gamma.time() < mMinCluTime)
              continue;

            bool matchFlag = 0,
                 isJpsi = 0,
                 isNotMIP = cluE > mCutMIPCluE,
                 isDispOK = testLambda(cluE, gamma.m02(), gamma.m20());

            if (pairMass > cfgeeMassMin && pairMass < cfgeeMassMax)
              isJpsi = 1;

            for (auto const& match : matches) {
              if (gamma.index() == match.caloClusterId()) {
                matchFlag = 1;
                break;
              }
            }

            TLorentzVector fourVectorP3;
            fourVectorP3.SetPxPyPzE(gamma.px(), gamma.py(), gamma.pz(), cluE);
            double tripletMass = (fourVectorP1 + fourVectorP2 + fourVectorP3).M(), tripletPt = (fourVectorP1 + fourVectorP2 + fourVectorP3).Pt();

            mHistManager.fill(HIST("h_TPCeePhosGamma_MS_v_3pt_v_cent"), tripletMass, tripletPt, cent);
            mHistManager.fill(HIST("h_TPCeePhosGamma_minusee_MS_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);

            mHistManager.fill(HIST("h_TPCeePhosGamma_MS_v_cluE"), tripletMass, cluE, cent);
            mHistManager.fill(HIST("h_TPCeePhosGamma_minusee_MS_v_cluE"), tripletMass - pairMass + cfgJpsiMass, cluE, cent);

            if (matchFlag) {
              mHistManager.fill(HIST("h_TPCeePhosGamma_MS_withMatches_v_3pt_v_cent"), tripletMass, tripletPt, cent);
              mHistManager.fill(HIST("h_TPCeePhosGamma_minusee_MS_withMatches_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
            } else {
              mHistManager.fill(HIST("h_TPCeePhosGamma_MS_noMatches_v_3pt_v_cent"), tripletMass, tripletPt, cent);
              mHistManager.fill(HIST("h_TPCeePhosGamma_minusee_MS_noMatches_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
              if (isJpsi) {
                mHistManager.fill(HIST("h_TPCeePhosGamma_MS_noMatches_aroundJpsi_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                mHistManager.fill(HIST("h_TPCeePhosGamma_minusee_MS_noMatches_aroundJpsi_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
                if (isDispOK) {
                  mHistManager.fill(HIST("h_TPCeePhosGamma_MS_noMatches_aroundJpsi_DispOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                  mHistManager.fill(HIST("h_TPCeePhosGamma_minusee_MS_noMatches_aroundJpsi_DispOK_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
                } else {
                  mHistManager.fill(HIST("h_TPCeePhosGamma_MS_noMatches_aroundJpsi_DispNotOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                  mHistManager.fill(HIST("h_TPCeePhosGamma_minusee_MS_noMatches_aroundJpsi_DispNotOK_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
                }
              } else {
                mHistManager.fill(HIST("h_TPCeePhosGamma_MS_noMatches_notAroundJpsi_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                mHistManager.fill(HIST("h_TPCeePhosGamma_minusee_MS_noMatches_notAroundJpsi_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
              }
              if (isDispOK) {
                mHistManager.fill(HIST("h_TPCeePhosGamma_MS_noMatches_DispOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                mHistManager.fill(HIST("h_TPCeePhosGamma_minusee_MS_noMatches_DispOK_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
              } else {
                mHistManager.fill(HIST("h_TPCeePhosGamma_MS_noMatches_DispNotOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                mHistManager.fill(HIST("h_TPCeePhosGamma_minusee_MS_noMatches_DispNotOK_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
              }
            }

            if (isJpsi) {
              mHistManager.fill(HIST("h_TPCeePhosGamma_MS_aroundJpsi_v_3pt_v_cent"), tripletMass, tripletPt, cent);
              mHistManager.fill(HIST("h_TPCeePhosGamma_minusee_MS_aroundJpsi_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
            } else {
              mHistManager.fill(HIST("h_TPCeePhosGamma_MS_notAroundJpsi_v_3pt_v_cent"), tripletMass, tripletPt, cent);
              mHistManager.fill(HIST("h_TPCeePhosGamma_minusee_MS_notAroundJpsi_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
            }

            if (isDispOK) {
              mHistManager.fill(HIST("h_TPCeePhosGamma_MS_DispOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
              mHistManager.fill(HIST("h_TPCeePhosGamma_minusee_MS_DispOK_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
            } else {
              mHistManager.fill(HIST("h_TPCeePhosGamma_MS_DispNotOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
              mHistManager.fill(HIST("h_TPCeePhosGamma_minusee_MS_DispNotOK_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
            }

            if (isNotMIP) {
              mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_MS_v_3pt_v_cent"), tripletMass, tripletPt, cent);
              mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_minusee_MS_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
              if (matchFlag) {
                mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_MS_withMatches_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_minusee_MS_withMatches_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
              } else {
                mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_MS_noMatches_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_minusee_MS_noMatches_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
                if (isJpsi) {
                  mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_MS_noMatches_aroundJpsi_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                  mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_minusee_MS_noMatches_aroundJpsi_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
                  if (isDispOK) {
                    mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_MS_noMatches_aroundJpsi_DispOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                    mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_minusee_MS_noMatches_aroundJpsi_DispOK_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
                  } else {
                    mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_MS_noMatches_aroundJpsi_DispNotOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                    mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_minusee_MS_noMatches_aroundJpsi_DispNotOK_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
                  }
                } else {
                  mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_MS_noMatches_notAroundJpsi_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                  mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_minusee_MS_noMatches_notAroundJpsi_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
                }
                if (isDispOK) {
                  mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_MS_noMatches_DispOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                  mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_minusee_MS_noMatches_DispOK_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
                } else {
                  mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_MS_noMatches_DispNotOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                  mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_minusee_MS_noMatches_DispNotOK_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
                }
              }

              if (isJpsi) {
                mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_MS_aroundJpsi_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_minusee_MS_aroundJpsi_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
              } else {
                mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_MS_notAroundJpsi_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_minusee_MS_notAroundJpsi_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
              }

              if (isDispOK) {
                mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_MS_DispOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_minusee_MS_DispOK_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
              } else {
                mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_MS_DispNotOK_v_3pt_v_cent"), tripletMass, tripletPt, cent);
                mHistManager.fill(HIST("h_TPCeePhosGamma_notMIP_minusee_MS_DispNotOK_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
              }
            } else {
              mHistManager.fill(HIST("h_TPCeePhosGamma_isMIP_MS_v_3pt_v_cent"), tripletMass, tripletPt, cent);
              mHistManager.fill(HIST("h_TPCeePhosGamma_isMIP_minusee_MS_v_3pt_v_cent"), tripletMass - pairMass + cfgJpsiMass, tripletPt, cent);
            }
          }
        }
      }
    }

    for (auto const& track : tracks) {
      mHistManager.fill(HIST("hTrackPt"), track.pt());
      mHistManager.fill(HIST("hTrackEta"), track.eta());
      mHistManager.fill(HIST("hTrackVX"), track.x());
      mHistManager.fill(HIST("hTrackVY"), track.y());
      mHistManager.fill(HIST("hTrackVZ"), track.z());
      mHistManager.fill(HIST("hTPCspectra"), track.pt(), track.tpcSignal());
    }
    for (auto const& track : filteredTracks) {
      if (!track.has_collision() || !track.hasTPC())
        continue;
      if (track.itsChi2NCl() > cfgITSchi2Max || track.tpcChi2NCl() > cfgTPCchi2Max)
        continue;
      if (track.itsNCls() < cfgITSnclsMin || track.itsNCls() > cfgITSnclsMax || !((track.itsClusterMap() & uint8_t(1)) > 0))
        continue;
      if (track.tpcNClsFound() < cfgTPCnclsMin || track.tpcNClsFound() > cfgTPCnclsMax)
        continue;
      if (track.tpcNClsCrossedRows() < cfgTPCnclsCRMin || track.tpcNClsCrossedRows() > cfgTPCnclsCRMax)
        continue;
      mHistManager.fill(HIST("hTPCspectra_isElectronRej"), track.pt(), track.tpcSignal());
      mHistManager.fill(HIST("hTrackPt_Cut"), track.pt());
      mHistManager.fill(HIST("hTrackEta_Cut"), track.eta());
      mHistManager.fill(HIST("hTrackVX_Cut"), track.x());
      mHistManager.fill(HIST("hTrackVY_Cut"), track.y());
      mHistManager.fill(HIST("hTrackVZ_Cut"), track.z());
    }
  }
  //_____________________________________________________________________________
  bool testLambda(float pt, float l1, float l2)
  {
    // Parameterization for full dispersion
    float l2Mean = 1.53126 + 9.50835e+06 / (1. + 1.08728e+07 * pt + 1.73420e+06 * pt * pt);
    float l1Mean = 1.12365 + 0.123770 * std::exp(-pt * 0.246551) + 5.30000e-03 * pt;
    float l2Sigma = 6.48260e-02 + 7.60261e+10 / (1. + 1.53012e+11 * pt + 5.01265e+05 * pt * pt) + 9.00000e-03 * pt;
    float l1Sigma = 4.44719e-04 + 6.99839e-01 / (1. + 1.22497e+00 * pt + 6.78604e-07 * pt * pt) + 9.00000e-03 * pt;
    float c = -0.35 - 0.550 * std::exp(-0.390730 * pt);

    return 0.5 * (l1 - l1Mean) * (l1 - l1Mean) / l1Sigma / l1Sigma +
             0.5 * (l2 - l2Mean) * (l2 - l2Mean) / l2Sigma / l2Sigma +
             0.5 * c * (l1 - l1Mean) * (l2 - l2Mean) / l1Sigma / l2Sigma <
           4.;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{
    adaptAnalysisTask<PhosElId>(cfgc),
    adaptAnalysisTask<MassSpectra>(cfgc),
    adaptAnalysisTask<TpcElIdMassSpectrum>(cfgc)};
  return workflow;
}
