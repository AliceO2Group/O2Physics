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

/// \struct PHOS electron id analysis
/// \brief Task for calculating electron identification parameters
/// \author Yeghishe Hambardzumyan, MIPT
/// \since Apr, 2024
/// @note Inherits functions and variables from phosAlign & phosPi0.
/// @note Results will be used for candidate table producing task.
///

using namespace o2;
using namespace o2::soa;
using namespace o2::aod::evsel;
using namespace o2::framework;
using namespace o2::framework::expressions;

namespace o2::aod
{
namespace PHOSMatch
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
} // namespace PHOSMatch

DECLARE_SOA_TABLE(PHOSMatchindexTable, "AOD", "PHSMTCH",                                                   //!
                  o2::soa::Index<>, PHOSMatch::CollisionId, PHOSMatch::CaloClusterId, PHOSMatch::TrackId); //!

} // namespace o2::aod

struct phosElId {

  Produces<o2::aod::PHOSMatchindexTable> phosMatch;

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using myTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::pidTOFFullEl, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;

  Configurable<float> mMinCluE{"mMinCluE", 0.3, "Minimum cluster energy for analysis"};
  Configurable<float> mMinCluTime{"minCluTime", -25.e-9, "Min. cluster time"};
  Configurable<float> mMaxCluTime{"mMaxCluTime", 25.e-9, "Max. cluster time"};
  Configurable<int> mEvSelTrig{"mEvSelTrig", kTVXinPHOS, "Select events with this trigger"};
  Configurable<int> mMinCluNcell{"minCluNcell", 3, "min cells in cluster"};

  Configurable<int> nBinsDeltaX{"nBinsDeltaX", 500, "N bins for track and cluster coordinate delta"};
  Configurable<float> mDeltaXmin{"mDeltaXmin", -100., "Min for track and cluster coordinate delta"};
  Configurable<float> mDeltaXmax{"mDeltaXmax", 100., "Max for track and cluster coordinate delta"};

  Configurable<int> nBinsDeltaZ{"nBinsDeltaZ", 500, "N bins for track and cluster coordinate delta"};
  Configurable<float> mDeltaZmin{"mDeltaZmin", -100., "Min for track and cluster coordinate delta"};
  Configurable<float> mDeltaZmax{"mDeltaZmax", 100., "Max for track and cluster coordinate delta"};

  Configurable<int> nBinsEp{"nBinsEp", 400, "N bins for E/p histograms"};
  Configurable<float> mEpmin{"mEpmin", -1., "Min for E/p histograms"};
  Configurable<float> mEpmax{"mEpmax", 3., "Max for E/p histograms"};

  Configurable<std::vector<float>> pSigma_dz{"pSigma_dz", {0.642, 0., 1.77, 2.725, 0.}, "parameters for sigma dz function"};
  Configurable<std::vector<float>> pSigma_dx{"pSigma_dx", {2.17769, 1.60275, 2.24136}, "parameters for sigma dx function"};

  Configurable<std::vector<float>> pPhosShiftZ{"pPhosShiftZ", {4.78838, 2.75138, 1.40825, 2.28735}, "Phos coordinate centering Z per module"};
  Configurable<std::vector<float>> pPhosShiftX{"pPhosShiftX", {2.158702, -1.526772, -0.814658, -1.852678}, "Phos coordinate centering X per module"};

  Configurable<std::vector<float>> pMean_dx_pos_mod1{"pMean_dx_pos_mod1", {-10.57, -0.42, 1.06}, "parameters for mean dx function on module 1 for positive tracks"};
  Configurable<std::vector<float>> pMean_dx_pos_mod2{"pMean_dx_pos_mod2", {-8.1, -0.42, 1.14}, "parameters for mean dx function on module 2 for positive tracks"};
  Configurable<std::vector<float>> pMean_dx_pos_mod3{"pMean_dx_pos_mod3", {-8.34, -0.42, 1.04}, "parameters for mean dx function on module 3 for positive tracks"};
  Configurable<std::vector<float>> pMean_dx_pos_mod4{"pMean_dx_pos_mod4", {-7.38, -0.42, 1.17}, "parameters for mean dx function on module 4 for positive tracks"};

  Configurable<std::vector<float>> pMean_dx_neg_mod1{"pMean_dx_neg_mod1", {9.92, -0.42, 1.29}, "parameters for mean dx function on module 1 for negative tracks"};
  Configurable<std::vector<float>> pMean_dx_neg_mod2{"pMean_dx_neg_mod2", {7.82, -0.4, 1.34}, "parameters for mean dx function on module 2 for negative tracks"};
  Configurable<std::vector<float>> pMean_dx_neg_mod3{"pMean_dx_neg_mod3", {8.45, -0.33, 1.5}, "parameters for mean dx function on module 3 for negative tracks"};
  Configurable<std::vector<float>> pMean_dx_neg_mod4{"pMean_dx_neg_mod4", {7.5, -0.42, 1.25}, "parameters for mean dx function on module 4 for negative tracks"};

  Configurable<float> cfgEtaMax{"cfgEtaMax", {0.8f}, "eta ranges"};
  Configurable<float> cfgPtMin{"cfgPtMin", {0.2f}, "pt min"};
  Configurable<float> cfgPtMax{"cfgPtMax", {20.f}, "pt max"};
  Configurable<float> cfgDCAxyMax{"cfgDCAxyMax", {3.f}, "dcaxy max"};
  Configurable<float> cfgDCAzMax{"cfgDCAzMax", {3.f}, "dcaz max"};
  Configurable<float> cfgITSchi2Max{"cfgITSchi2Max", {5.f}, "its chi2 max"};
  Configurable<float> cfgITSnclsMin{"cfgITSnclsMin", {4.5f}, "min number of ITS clusters"};
  Configurable<float> cfgITSnclsMax{"cfgITSnclsMax", {7.5f}, "max number of ITS clusters"};
  Configurable<float> cfgTPCchi2Max{"cfgTPCchi2Max", {4.f}, "tpc chi2 max"};
  Configurable<float> cfgTPCnclsMin{"cfgTPCnclsMin", {90.f}, "min number of TPC clusters"};
  Configurable<float> cfgTPCnclsMax{"cfgTPCnclsMax", {170.f}, "max number of TPC clusters"};
  Configurable<float> cfgTPCnclsCRMin{"cfgTPCnclsCRMin", {80.f}, "min number of TPC crossed rows"};
  Configurable<float> cfgTPCnclsCRMax{"cfgTPCnclsCRMax", {161.f}, "max number of TPC crossed rows"};
  Configurable<float> cfgTPCNSigmaElMin{"cfgTPCNSigmaElMin", {-3.f}, "min TPC nsigma e for inclusion"};
  Configurable<float> cfgTPCNSigmaElMax{"cfgTPCNSigmaElMax", {2.f}, "max TPC nsigma e for inclusion"};
  Configurable<float> cfgTPCNSigmaPiMin{"cfgTPCNSigmaPiMin", {-3.f}, "min TPC nsigma pion for exclusion"};
  Configurable<float> cfgTPCNSigmaPiMax{"cfgTPCNSigmaPiMax", {3.5f}, "max TPC nsigma pion for exclusion"};
  Configurable<float> cfgTPCNSigmaPrMin{"cfgTPCNSigmaPrMin", {-3.f}, "min TPC nsigma proton for exclusion"};
  Configurable<float> cfgTPCNSigmaPrMax{"cfgTPCNSigmaPrMax", {4.f}, "max TPC nsigma proton for exclusion"};
  Configurable<float> cfgTPCNSigmaKaMin{"cfgTPCNSigmaKaMin", {-3.f}, "min TPC nsigma kaon for exclusion"};
  Configurable<float> cfgTPCNSigmaKaMax{"cfgTPCNSigmaKaMax", {4.f}, "max TPC nsigma kaon for exclusion"};
  Configurable<float> cfgTOFNSigmaElMin{"cfgTOFNSigmaElMin", {-3.f}, "min TOF nsigma e for inclusion"};
  Configurable<float> cfgTOFNSigmaElMax{"cfgTOFNSigmaElMax", {3.f}, "max TOF nsigma e for inclusion"};

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

  HistogramRegistry mHistManager{"phosElIdHistograms"};
  TF1 *fSigma_dz, *fSigma_dx;
  float *PhosShiftX, *PhosShiftZ;

  TF1* fMean_dx_pos_mod1;
  TF1* fMean_dx_neg_mod1;
  TF1* fMean_dx_pos_mod2;
  TF1* fMean_dx_neg_mod2;
  TF1* fMean_dx_pos_mod3;
  TF1* fMean_dx_neg_mod3;
  TF1* fMean_dx_pos_mod4;
  TF1* fMean_dx_neg_mod4;

  void init(InitContext const&)
  {
    LOG(info) << "Initializing PHOS electron identification analysis task ...";

    std::vector<double> momentum_binning = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                                            1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                                            4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 10.};
    std::vector<float> parameters_sigma_dz = pSigma_dz;
    std::vector<float> parameters_sigma_dx = pSigma_dx;
    std::vector<float> vPhosShiftX = pPhosShiftX;
    std::vector<float> vPhosShiftZ = pPhosShiftZ;
    std::vector<float> Mean_dx_pos_mod1 = pMean_dx_pos_mod1;
    std::vector<float> Mean_dx_pos_mod2 = pMean_dx_pos_mod2;
    std::vector<float> Mean_dx_pos_mod3 = pMean_dx_pos_mod3;
    std::vector<float> Mean_dx_pos_mod4 = pMean_dx_pos_mod4;
    std::vector<float> Mean_dx_neg_mod1 = pMean_dx_neg_mod1;
    std::vector<float> Mean_dx_neg_mod2 = pMean_dx_neg_mod2;
    std::vector<float> Mean_dx_neg_mod3 = pMean_dx_neg_mod3;
    std::vector<float> Mean_dx_neg_mod4 = pMean_dx_neg_mod4;

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

    const AxisSpec
      axisCounter{1, 0, +1, ""},
      axisP{momentum_binning, "p (GeV/c)"},
      axisPt{momentum_binning, "p_{T} (GeV/c)"},
      axisEta{200, -0.2, 0.2, "#eta"},
      axisPhi{80, 240, 320, "#varphi"},
      axisE{200, 0, 10, "E (GeV)", "E (GeV)"},
      axisMassSpectrum{4000, 0, 4, "M (GeV/c^{2})", "Mass e^{+}e^{-} (GeV/c^{2})"},
      axisEp{nBinsEp, mEpmin, mEpmax, "E/p", "E_{cluster}/p_{track}"},
      axisdX{nBinsDeltaX, mDeltaXmin, mDeltaXmax, "x_{tr}-x_{clu} (cm)", "x_{tr}-x_{clu} (cm)"},
      axisdZ{nBinsDeltaZ, mDeltaZmin, mDeltaZmax, "z_{tr}-z_{clu} (cm)", "z_{tr}-z_{clu} (cm)"},
      axisCells{20, 0., 20., "number of cells", "number of cells"},
      axisTime{100, 2e9 * mMinCluTime, 2e9 * mMaxCluTime, "time (ns)", "time (nanoseconds)"},
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
    mHistManager.add("hCluE_v_pt_1sigma", "Cluster energy vs p within trackmatch 1sigma", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("hCluE_v_pt_1sigma_disp", "Cluster energy vs p within trackmatch 1sigma | OK dispersion", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("hCluE_v_pt_2sigma", "Cluster energy vs p within trackmatch 2sigma", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("hCluE_v_pt_2sigma_disp", "Cluster energy vs p within trackmatch 2sigma | OK dispersion", HistType::kTH3F, {axisE, axisPt, axisModes});

    mHistManager.add("hCluE_v_pt_disp_TPC", "Cluster energy vs p | OK dispersion", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("hCluE_v_pt_1sigma_TPC", "Cluster energy vs p within trackmatch 1sigma", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("hCluE_v_pt_1sigma_disp_TPC", "Cluster energy vs p within trackmatch 1sigma | OK dispersion", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("hCluE_v_pt_2sigma_TPC", "Cluster energy vs p within trackmatch 2sigma", HistType::kTH3F, {axisE, axisPt, axisModes});
    mHistManager.add("hCluE_v_pt_2sigma_disp_TPC", "Cluster energy vs p within trackmatch 2sigma | OK dispersion", HistType::kTH3F, {axisE, axisPt, axisModes});

    mHistManager.add("hEp_v_pt_disp", "E/p ratio vs p | OK dispersion", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("hEp_v_pt_1sigma", "E/p ratio vs p within trackmatch 1sigma", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("hEp_v_pt_1sigma_disp", "E/p ratio vs p within trackmatch 1sigma | OK dispersion", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("hEp_v_pt_2sigma", "E/p ratio vs p within trackmatch 2sigma", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("hEp_v_pt_2sigma_disp", "E/p ratio vs p within trackmatch 2sigma | OK dispersion", HistType::kTH3F, {axisEp, axisPt, axisModes});

    mHistManager.add("hEp_v_E_disp", "E/p ratio vs cluster E | OK dispersion", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("hEp_v_E_1sigma", "E/p ratio vs cluster E within trackmatch 1sigma", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("hEp_v_E_1sigma_disp", "E/p ratio vs cluster E within trackmatch 1sigma | OK dispersion", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("hEp_v_E_2sigma", "E/p ratio vs cluster E within trackmatch 2sigma", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("hEp_v_E_2sigma_disp", "E/p ratio vs cluster E within trackmatch 2sigma | OK dispersion", HistType::kTH3F, {axisEp, axisE, axisModes});

    mHistManager.add("hEp_v_pt_disp_TPC", "E/p ratio vs p | OK dispersion", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("hEp_v_pt_1sigma_TPC", "E/p ratio vs p within trackmatch 1sigma", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("hEp_v_pt_1sigma_disp_TPC", "E/p ratio vs p within trackmatch 1sigma | OK dispersion", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("hEp_v_pt_2sigma_TPC", "E/p ratio vs p within trackmatch 2sigma", HistType::kTH3F, {axisEp, axisPt, axisModes});
    mHistManager.add("hEp_v_pt_2sigma_disp_TPC", "E/p ratio vs p within trackmatch 2sigma | OK dispersion", HistType::kTH3F, {axisEp, axisPt, axisModes});

    mHistManager.add("hEp_v_E_disp_TPC", "E/p ratio vs cluster E | OK dispersion", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("hEp_v_E_1sigma_TPC", "E/p ratio vs cluster E within trackmatch 1sigma", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("hEp_v_E_1sigma_disp_TPC", "E/p ratio vs cluster E within trackmatch 1sigma | OK dispersion", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("hEp_v_E_2sigma_TPC", "E/p ratio vs cluster E within trackmatch 2sigma", HistType::kTH3F, {axisEp, axisE, axisModes});
    mHistManager.add("hEp_v_E_2sigma_disp_TPC", "E/p ratio vs cluster E within trackmatch 2sigma | OK dispersion", HistType::kTH3F, {axisEp, axisE, axisModes});

    geomPHOS = std::make_unique<o2::phos::Geometry>("PHOS");
    fSigma_dz = new TF1("fSigma_dz", "[0]/(x+[1])^[2]+pol1(3)", 0.3, 10);
    fSigma_dz->SetParameters(parameters_sigma_dz.at(0), parameters_sigma_dz.at(1), parameters_sigma_dz.at(2), parameters_sigma_dz.at(3), parameters_sigma_dz.at(4));

    fSigma_dx = new TF1("fSigma_dx", "[0]/x^[1]+[2]", 0.1, 10);
    fSigma_dx->SetParameters(parameters_sigma_dx.at(0), parameters_sigma_dx.at(1), parameters_sigma_dx.at(2));

    fMean_dx_pos_mod1 = new TF1("funcMeandx_pos_mod1", "[0]/(x+[1])^[2]", 0.1, 10);
    fMean_dx_neg_mod1 = new TF1("funcMeandx_neg_mod1", "[0]/(x+[1])^[2]", 0.1, 10);
    fMean_dx_pos_mod2 = new TF1("funcMeandx_pos_mod2", "[0]/(x+[1])^[2]", 0.1, 10);
    fMean_dx_neg_mod2 = new TF1("funcMeandx_neg_mod2", "[0]/(x+[1])^[2]", 0.1, 10);
    fMean_dx_pos_mod3 = new TF1("funcMeandx_pos_mod3", "[0]/(x+[1])^[2]", 0.1, 10);
    fMean_dx_neg_mod3 = new TF1("funcMeandx_neg_mod3", "[0]/(x+[1])^[2]", 0.1, 10);
    fMean_dx_pos_mod4 = new TF1("funcMeandx_pos_mod4", "[0]/(x+[1])^[2]", 0.1, 10);
    fMean_dx_neg_mod4 = new TF1("funcMeandx_neg_mod4", "[0]/(x+[1])^[2]", 0.1, 10);

    fMean_dx_pos_mod1->SetParameters(Mean_dx_pos_mod1.at(0), Mean_dx_pos_mod1.at(1), Mean_dx_pos_mod1.at(2));
    fMean_dx_pos_mod2->SetParameters(Mean_dx_pos_mod2.at(0), Mean_dx_pos_mod2.at(1), Mean_dx_pos_mod2.at(2));
    fMean_dx_pos_mod3->SetParameters(Mean_dx_pos_mod3.at(0), Mean_dx_pos_mod3.at(1), Mean_dx_pos_mod3.at(2));
    fMean_dx_pos_mod4->SetParameters(Mean_dx_pos_mod4.at(0), Mean_dx_pos_mod4.at(1), Mean_dx_pos_mod4.at(2));

    fMean_dx_neg_mod1->SetParameters(Mean_dx_neg_mod1.at(0), Mean_dx_neg_mod1.at(1), Mean_dx_neg_mod1.at(2));
    fMean_dx_neg_mod2->SetParameters(Mean_dx_neg_mod2.at(0), Mean_dx_neg_mod2.at(1), Mean_dx_neg_mod2.at(2));
    fMean_dx_neg_mod3->SetParameters(Mean_dx_neg_mod3.at(0), Mean_dx_neg_mod3.at(1), Mean_dx_neg_mod3.at(2));
    fMean_dx_neg_mod4->SetParameters(Mean_dx_neg_mod4.at(0), Mean_dx_neg_mod4.at(1), Mean_dx_neg_mod4.at(2));
  }
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
               aod::CaloClusters& clusters,
               soa::Filtered<myTracks> const& tracks,
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
    if (fabs(collision.posZ()) > 10.f)
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

      if (!track.has_collision() || fabs(track.dcaXY()) > cfgDCAxyMax || fabs(track.dcaZ()) > cfgDCAzMax || !track.hasTPC() || fabs(track.eta()) > 0.15)
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
        bool isTPC_electron = nsigmaTPCEl > cfgTPCNSigmaElMin && nsigmaTPCEl < cfgTPCNSigmaElMax;
        bool isTOF_electron = nsigmaTOFEl > cfgTOFNSigmaElMin && nsigmaTOFEl < cfgTOFNSigmaElMax;
        isElectron = isTPC_electron || isTOF_electron;

        float nsigmaTPCPi = track.tpcNSigmaPi();
        float nsigmaTPCKa = track.tpcNSigmaKa();
        float nsigmaTPCPr = track.tpcNSigmaPr();
        bool isPion = nsigmaTPCPi > cfgTPCNSigmaPiMin && nsigmaTPCPi < cfgTPCNSigmaPiMax;
        bool isKaon = nsigmaTPCKa > cfgTPCNSigmaKaMin && nsigmaTPCKa < cfgTPCNSigmaKaMax;
        bool isProton = nsigmaTPCPr > cfgTPCNSigmaPrMin && nsigmaTPCPr < cfgTPCNSigmaPrMax;
        if (isElectron && !(isPion || isKaon || isProton))
          isElectron = true;
      }

      bool posTrack = (track.sign() > 0 && bz > 0) || (track.sign() < 0 && bz < 0);
      for (auto const& clu : clusters) {
        if (module != clu.mod())
          continue;
        double cluE = clu.e();
        mHistManager.fill(HIST("hCluE_ncells_mod"), cluE, clu.ncell(), module);

        if (cluE < mMinCluE ||
            clu.ncell() < mMinCluNcell ||
            clu.time() > mMaxCluTime || clu.time() < mMinCluTime)
          continue;

        bool isDispOK = testLambda(cluE, clu.m02(), clu.m20());

        float posX = clu.x(), posZ = clu.z(), dX = trackX - posX, dZ = trackZ - posZ, Ep = cluE / trackMom;

        mHistManager.fill(HIST("hCluXZ_mod"), posX, posZ, module);

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
        if (!isWithin2Sigma(module, trackPT, dZ, dX))
          continue;
        mHistManager.fill(HIST("hCluE_v_pt_2sigma"), cluE, trackPT, module);
        mHistManager.fill(HIST("hEp_v_pt_2sigma"), Ep, trackPT, module);
        mHistManager.fill(HIST("hEp_v_E_2sigma"), Ep, cluE, module);
        if (isElectron) {
          mHistManager.fill(HIST("hCluE_v_pt_2sigma_TPC"), cluE, trackPT, module);
          mHistManager.fill(HIST("hEp_v_pt_2sigma_TPC"), Ep, trackPT, module);
          mHistManager.fill(HIST("hEp_v_E_2sigma_TPC"), Ep, cluE, module);
        }
        if (isDispOK) {
          mHistManager.fill(HIST("hCluE_v_pt_2sigma_disp"), cluE, trackPT, module);
          mHistManager.fill(HIST("hEp_v_pt_2sigma_disp"), Ep, trackPT, module);
          mHistManager.fill(HIST("hEp_v_E_2sigma_disp"), Ep, cluE, module);
          if (isElectron) {
            mHistManager.fill(HIST("hCluE_v_pt_2sigma_disp_TPC"), cluE, trackPT, module);
            mHistManager.fill(HIST("hEp_v_pt_2sigma_disp_TPC"), Ep, trackPT, module);
            mHistManager.fill(HIST("hEp_v_E_2sigma_disp_TPC"), Ep, cluE, module);
          }
          phosMatch(collision.index(), clu.index(), track.index());
        }
        if (isWithin1Sigma(module, trackPT, dZ, dX)) {
          mHistManager.fill(HIST("hCluE_v_pt_1sigma"), cluE, trackPT, module);
          mHistManager.fill(HIST("hEp_v_pt_1sigma"), Ep, trackPT, module);
          mHistManager.fill(HIST("hEp_v_E_1sigma"), Ep, cluE, module);
          if (isElectron) {
            mHistManager.fill(HIST("hCluE_v_pt_1sigma_TPC"), cluE, trackPT, module);
            mHistManager.fill(HIST("hEp_v_pt_1sigma_TPC"), Ep, trackPT, module);
            mHistManager.fill(HIST("hEp_v_E_1sigma_TPC"), Ep, cluE, module);
          }
          if (isDispOK) {
            mHistManager.fill(HIST("hCluE_v_pt_1sigma_disp"), cluE, trackPT, module);
            mHistManager.fill(HIST("hEp_v_pt_1sigma_disp"), Ep, trackPT, module);
            mHistManager.fill(HIST("hEp_v_E_1sigma_disp"), Ep, cluE, module);
            if (isElectron) {
              mHistManager.fill(HIST("hCluE_v_pt_1sigma_disp_TPC"), cluE, trackPT, module);
              mHistManager.fill(HIST("hEp_v_pt_1sigma_disp_TPC"), Ep, trackPT, module);
              mHistManager.fill(HIST("hEp_v_E_1sigma_disp_TPC"), Ep, cluE, module);
            }
          }
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
            if (testLambda(cluE, clu.m02(), clu.m20()))
              mHistManager.fill(HIST("hCluE_mod_disp"), cluE, mod);
          }
        }
      }
    } // end of cluster loop
  }

  bool isWithin2Sigma(int16_t& mod, float p, float deltaZ, float deltaX)
  {
    if (mod == 1) {
      if (fabs(deltaZ - PhosShiftZ[0]) > 2 * fSigma_dz->Eval(p))
        return false;
      if (fabs(deltaX - fMean_dx_pos_mod1->Eval(p) + PhosShiftX[0]) > 2 * fSigma_dx->Eval(p))
        return false;
    } else if (mod == 2) {
      if (fabs(deltaZ - PhosShiftZ[1]) > 2 * fSigma_dz->Eval(p))
        return false;
      if (fabs(deltaX - fMean_dx_pos_mod2->Eval(p) + PhosShiftX[1]) > 2 * fSigma_dx->Eval(p))
        return false;
    } else if (mod == 3) {
      if (fabs(deltaZ - PhosShiftZ[2]) > 2 * fSigma_dz->Eval(p))
        return false;
      if (fabs(deltaX - fMean_dx_pos_mod3->Eval(p) + PhosShiftX[2]) > 2 * fSigma_dx->Eval(p))
        return false;
    } else if (mod == 4) {
      if (fabs(deltaZ - PhosShiftZ[3]) > 2 * fSigma_dz->Eval(p))
        return false;
      if (fabs(deltaX - fMean_dx_pos_mod4->Eval(p) + PhosShiftX[3]) > 2 * fSigma_dx->Eval(p))
        return false;
    }
    return true;
  }

  bool isWithin1Sigma(int16_t& mod, float p, float deltaZ, float deltaX)
  {
    if (mod == 1) {
      if (fabs(deltaZ - PhosShiftZ[0]) > fSigma_dz->Eval(p))
        return false;
      if (fabs(deltaX - fMean_dx_pos_mod1->Eval(p) + PhosShiftX[0]) > fSigma_dx->Eval(p))
        return false;
    } else if (mod == 2) {
      if (fabs(deltaZ - PhosShiftZ[1]) > fSigma_dz->Eval(p))
        return false;
      if (fabs(deltaX - fMean_dx_pos_mod2->Eval(p) + PhosShiftX[1]) > fSigma_dx->Eval(p))
        return false;
    } else if (mod == 3) {
      if (fabs(deltaZ - PhosShiftZ[2]) > fSigma_dz->Eval(p))
        return false;
      if (fabs(deltaX - fMean_dx_pos_mod3->Eval(p) + PhosShiftX[2]) > fSigma_dx->Eval(p))
        return false;
    } else if (mod == 4) {
      if (fabs(deltaZ - PhosShiftZ[3]) > fSigma_dz->Eval(p))
        return false;
      if (fabs(deltaX - fMean_dx_pos_mod4->Eval(p) + PhosShiftX[3]) > fSigma_dx->Eval(p))
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
    if (trackPhi < phiMin || trackPhi > phiMax || abs(trackEta) > etaMax) {
      return false;
    }

    const float dphi = 20. * 0.017453293;
    if (trackPhi < 0.) {
      trackPhi += TMath::TwoPi();
    }
    if (trackPhi > TMath::TwoPi()) {
      trackPhi -= TMath::TwoPi();
    }
    module = 1 + static_cast<int16_t>((trackPhi - phiMin) / dphi);
    if (module < 1) {
      module = 1;
    }
    if (module > 4) {
      module = 4;
    }

    // get PHOS radius
    constexpr float shiftY = -1.26;    // Depth-optimized
    double posL[3] = {0., 0., shiftY}; // local position at the center of module
    double posG[3] = {0};
    geomPHOS->getAlignmentMatrix(module)->LocalToMaster(posL, posG);
    double rPHOS = sqrt(posG[0] * posG[0] + posG[1] * posG[1]);
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
    double ca = cos(alpha), sa = sin(alpha);
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
    float l1Mean = 1.12365 + 0.123770 * TMath::Exp(-pt * 0.246551) + 5.30000e-03 * pt;
    float l2Sigma = 6.48260e-02 + 7.60261e+10 / (1. + 1.53012e+11 * pt + 5.01265e+05 * pt * pt) + 9.00000e-03 * pt;
    float l1Sigma = 4.44719e-04 + 6.99839e-01 / (1. + 1.22497e+00 * pt + 6.78604e-07 * pt * pt) + 9.00000e-03 * pt;
    float c = -0.35 - 0.550 * TMath::Exp(-0.390730 * pt);

    return 0.5 * (l1 - l1Mean) * (l1 - l1Mean) / l1Sigma / l1Sigma +
             0.5 * (l2 - l2Mean) * (l2 - l2Mean) / l2Sigma / l2Sigma +
             0.5 * c * (l1 - l1Mean) * (l2 - l2Mean) / l1Sigma / l2Sigma <
           4.;
  }
};

struct massSpectra {

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using myTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::pidTOFFullEl, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;

  Configurable<int> mEvSelTrig{"mEvSelTrig", kTVXinPHOS, "Select events with this trigger"};

  Configurable<float> cfgEtaMax{"cfgEtaMax", {0.8f}, "eta ranges"};
  Configurable<float> cfgPtMin{"cfgPtMin", {0.2f}, "pt min"};
  Configurable<float> cfgPtMax{"cfgPtMax", {20.f}, "pt max"};
  Configurable<float> cfgDCAxyMax{"cfgDCAxyMax", {3.f}, "dcaxy max"};
  Configurable<float> cfgDCAzMax{"cfgDCAzMax", {3.f}, "dcaz max"};
  Configurable<float> cfgITSchi2Max{"cfgITSchi2Max", {5.f}, "its chi2 max"};
  Configurable<float> cfgITSnclsMin{"cfgITSnclsMin", {4.5f}, "min number of ITS clusters"};
  Configurable<float> cfgITSnclsMax{"cfgITSnclsMax", {7.5f}, "max number of ITS clusters"};
  Configurable<float> cfgTPCchi2Max{"cfgTPCchi2Max", {4.f}, "tpc chi2 max"};
  Configurable<float> cfgTPCnclsMin{"cfgTPCnclsMin", {90.f}, "min number of TPC clusters"};
  Configurable<float> cfgTPCnclsMax{"cfgTPCnclsMax", {170.f}, "max number of TPC clusters"};
  Configurable<float> cfgTPCnclsCRMin{"cfgTPCnclsCRMin", {80.f}, "min number of TPC crossed rows"};
  Configurable<float> cfgTPCnclsCRMax{"cfgTPCnclsCRMax", {161.f}, "max number of TPC crossed rows"};

  Configurable<float> cfgTPCNSigmaElMin{"cfgTPCNSigmaElMin", {-3.f}, "min TPC nsigma e for inclusion"};
  Configurable<float> cfgTPCNSigmaElMax{"cfgTPCNSigmaElMax", {2.f}, "max TPC nsigma e for inclusion"};
  Configurable<float> cfgTPCNSigmaPiMin{"cfgTPCNSigmaPiMin", {-3.f}, "min TPC nsigma pion for exclusion"};
  Configurable<float> cfgTPCNSigmaPiMax{"cfgTPCNSigmaPiMax", {3.5f}, "max TPC nsigma pion for exclusion"};
  Configurable<float> cfgTPCNSigmaPrMin{"cfgTPCNSigmaPrMin", {-3.f}, "min TPC nsigma proton for exclusion"};
  Configurable<float> cfgTPCNSigmaPrMax{"cfgTPCNSigmaPrMax", {4.f}, "max TPC nsigma proton for exclusion"};
  Configurable<float> cfgTPCNSigmaKaMin{"cfgTPCNSigmaKaMin", {-3.f}, "min TPC nsigma kaon for exclusion"};
  Configurable<float> cfgTPCNSigmaKaMax{"cfgTPCNSigmaKaMax", {4.f}, "max TPC nsigma kaon for exclusion"};
  Configurable<float> cfgTOFNSigmaElMin{"cfgTOFNSigmaElMin", {-3.f}, "min TOF nsigma e for inclusion"};
  Configurable<float> cfgTOFNSigmaElMax{"cfgTOFNSigmaElMax", {3.f}, "max TOF nsigma e for inclusion"};

  Configurable<std::vector<float>> cfgEpSigmaPars{"cfgEpSigmaPars", {1.3e-02, 1.9e-02, 1.1e-02, 3.e-02}, "E/p sigma function parameters (from alice 3 mc tests + const)"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  std::unique_ptr<o2::phos::Geometry> geomPHOS;
  double bz{0.}; // magnetic field
  int runNumber{0};

  HistogramRegistry mHistManager{"MassSpectraHistograms"};
  TF1* fEp_sigma_phos;

  void init(InitContext const&)
  {
    LOG(info) << "Initializing PHOS electron identification analysis task ...";

    std::vector<double> momentum_binning = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                                            1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                                            4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 10.};

    const AxisSpec
      axisCounter{1, 0, +1, ""},
      axisPt{momentum_binning, "p_{T} (GeV/c)"},
      axisE{200, 0, 10, "E (GeV)", "E (GeV)"},
      axisMassSpectrum{4000, 0, 4, "M (GeV/c^{2})", "Mass e^{+}e^{-} (GeV/c^{2})"};

    mHistManager.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    mHistManager.add("TVXinPHOSCounter", "TVXinPHOSCounter", kTH1F, {axisCounter});

    mHistManager.add("h_eh_pp_mass_spectra_v_Pt", "Mass e^{+}h^{+} vs momentum e^{+}h^{+}", HistType::kTH2F, {axisMassSpectrum, axisPt});
    mHistManager.add("h_ee_pp_mass_spectra_v_Pt", "Mass e^{+}e^{+} vs momentum e^{+}e^{+}", HistType::kTH2F, {axisMassSpectrum, axisPt});
    mHistManager.add("h_eh_mm_mass_spectra_v_Pt", "Mass e^{-}h^{-} vs momentum e^{-}h^{-}", HistType::kTH2F, {axisMassSpectrum, axisPt});
    mHistManager.add("h_ee_mm_mass_spectra_v_Pt", "Mass e^{-}e^{-} vs momentum e^{-}e^{-}", HistType::kTH2F, {axisMassSpectrum, axisPt});

    mHistManager.add("h_eh_pp_mass_spectra_v_E", "Mass e^{+}h^{+} vs cluster E e^{+}h^{+}", HistType::kTH2F, {axisMassSpectrum, axisE});
    mHistManager.add("h_ee_pp_mass_spectra_v_E", "Mass e^{+}e^{+} vs cluster E e^{+}e^{+}", HistType::kTH2F, {axisMassSpectrum, axisE});
    mHistManager.add("h_eh_mm_mass_spectra_v_E", "Mass e^{-}h^{-} vs cluster E e^{-}h^{-}", HistType::kTH2F, {axisMassSpectrum, axisE});
    mHistManager.add("h_ee_mm_mass_spectra_v_E", "Mass e^{-}e^{-} vs cluster E e^{-}e^{-}", HistType::kTH2F, {axisMassSpectrum, axisE});

    mHistManager.add("h_eh_mp_mass_spectra_v_Pt", "Mass e^{#pm}h^{#mp} vs momentum e^{#pm}h^{#mp}", HistType::kTH2F, {axisMassSpectrum, axisPt});
    mHistManager.add("h_ee_mp_mass_spectra_v_Pt", "Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp}", HistType::kTH2F, {axisMassSpectrum, axisPt});
    mHistManager.add("h_eh_mp_mass_spectra_v_E", "Mass e^{#pm}h^{#mp} vs cluster E e^{#pm}h^{#mp}", HistType::kTH2F, {axisMassSpectrum, axisE});
    mHistManager.add("h_ee_mp_mass_spectra_v_E", "Mass e^{#pm}e^{#mp} vs cluster E e^{#pm}e^{#mp}", HistType::kTH2F, {axisMassSpectrum, axisE});

    geomPHOS = std::make_unique<o2::phos::Geometry>("PHOS");

    std::vector<float> EpSigmaPars = cfgEpSigmaPars;
    fEp_sigma_phos = new TF1("fEp_sigma_phos", "sqrt([0]*[0]/x/x+[1]*[1]/x+[2]*[2])+[3]", 0.01, 10);
    fEp_sigma_phos->SetParameters(EpSigmaPars.at(0), EpSigmaPars.at(1), EpSigmaPars.at(2), EpSigmaPars.at(3));
  }

  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
               aod::CaloClusters const& clusters,
               myTracks const& tracks,
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
    if (fabs(collision.posZ()) > 10.f)
      return;
    mHistManager.fill(HIST("eventCounter"), 0.5);
    if (!collision.alias_bit(mEvSelTrig))
      return;
    mHistManager.fill(HIST("TVXinPHOSCounter"), 0.5);

    if (clusters.size() == 0)
      return; // Nothing to process

    for (auto const& TPCel : tracks) {

      if (!TPCel.has_collision() || fabs(TPCel.dcaXY()) > cfgDCAxyMax || fabs(TPCel.dcaZ()) > cfgDCAzMax || !TPCel.hasTPC() || fabs(TPCel.eta()) > 0.15)
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
      bool isTPC_electron = nsigmaTPCEl > cfgTPCNSigmaElMin && nsigmaTPCEl < cfgTPCNSigmaElMax;
      bool isTOF_electron = nsigmaTOFEl > cfgTOFNSigmaElMin && nsigmaTOFEl < cfgTOFNSigmaElMax;
      isElectron = isTPC_electron || isTOF_electron;

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

      for (auto const& match : matches) {
        auto clust2 = clusters.iteratorAt(match.caloClusterId());
        auto track2 = tracks.iteratorAt(match.trackId());

        if (TPCel.collisionId() != track2.collisionId())
          continue;
        if (TPCel.index() >= track2.index())
          break;

        float mass_2tracks = 0, mom_2tracks = 0, cluE = clust2.e();
        TLorentzVector P1, P2;
        P1.SetPxPyPzE(TPCel.px(), TPCel.py(), TPCel.pz(), TPCel.energy(0));
        P2.SetPxPyPzE(track2.px(), track2.py(), track2.pz(), track2.energy(0));
        mom_2tracks = (P1 + P2).Pt();
        mass_2tracks = (P1 + P2).M();
        bool posTrack = (TPCel.sign() > 0 && bz > 0) || (TPCel.sign() < 0 && bz < 0);
        bool elCandidate = (fabs(clust2.e() / track2.p() - 1) < 2 * fEp_sigma_phos->Eval(mom_2tracks));

        if (TPCel.sign() == track2.sign()) {
          if (posTrack) {
            mHistManager.fill(HIST("h_eh_pp_mass_spectra_v_Pt"), mass_2tracks, mom_2tracks);
            mHistManager.fill(HIST("h_eh_pp_mass_spectra_v_E"), mass_2tracks, cluE);
            if (elCandidate) {
              mHistManager.fill(HIST("h_ee_pp_mass_spectra_v_Pt"), mass_2tracks, mom_2tracks);
              mHistManager.fill(HIST("h_ee_pp_mass_spectra_v_E"), mass_2tracks, cluE);
            }
          } else {
            mHistManager.fill(HIST("h_eh_mm_mass_spectra_v_Pt"), mass_2tracks, mom_2tracks);
            mHistManager.fill(HIST("h_eh_mm_mass_spectra_v_E"), mass_2tracks, cluE);
            if (elCandidate) {
              mHistManager.fill(HIST("h_ee_mm_mass_spectra_v_Pt"), mass_2tracks, mom_2tracks);
              mHistManager.fill(HIST("h_ee_mm_mass_spectra_v_E"), mass_2tracks, cluE);
            }
          }
        } else {
          mHistManager.fill(HIST("h_eh_mp_mass_spectra_v_Pt"), mass_2tracks, mom_2tracks);
          mHistManager.fill(HIST("h_eh_mp_mass_spectra_v_E"), mass_2tracks, cluE);
          if (elCandidate) {
            mHistManager.fill(HIST("h_ee_mp_mass_spectra_v_Pt"), mass_2tracks, mom_2tracks);
            mHistManager.fill(HIST("h_ee_mp_mass_spectra_v_E"), mass_2tracks, cluE);
          }
        }
      }
    } // end of double loop
  }
};

struct tpcElIdMassSpectrum {

  using SelCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using myTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksDCACov, aod::pidTOFFullEl, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr>;

  Configurable<float> cfgEtaMax{"cfgEtaMax", {0.8f}, "eta ranges"};
  Configurable<float> cfgPtMin{"cfgPtMin", {0.2f}, "pt min"};
  Configurable<float> cfgPtMax{"cfgPtMax", {20.f}, "pt max"};
  Configurable<float> cfgDCAxyMax{"cfgDCAxyMax", {3.f}, "dcaxy max"};
  Configurable<float> cfgDCAzMax{"cfgDCAzMax", {3.f}, "dcaz max"};
  Configurable<float> cfgITSchi2Max{"cfgITSchi2Max", {5.f}, "its chi2 max"};
  Configurable<float> cfgITSnclsMin{"cfgITSnclsMin", {4.5f}, "min number of ITS clusters"};
  Configurable<float> cfgITSnclsMax{"cfgITSnclsMax", {7.5f}, "max number of ITS clusters"};
  Configurable<float> cfgTPCchi2Max{"cfgTPCchi2Max", {4.f}, "tpc chi2 max"};
  Configurable<float> cfgTPCnclsMin{"cfgTPCnclsMin", {90.f}, "min number of TPC clusters"};
  Configurable<float> cfgTPCnclsMax{"cfgTPCnclsMax", {170.f}, "max number of TPC clusters"};
  Configurable<float> cfgTPCnclsCRMin{"cfgTPCnclsCRMin", {80.f}, "min number of TPC crossed rows"};
  Configurable<float> cfgTPCnclsCRMax{"cfgTPCnclsCRMax", {161.f}, "max number of TPC crossed rows"};
  Configurable<float> cfgTPCNSigmaElMin{"cfgTPCNSigmaElMin", {-3.f}, "min TPC nsigma e for inclusion"};
  Configurable<float> cfgTPCNSigmaElMax{"cfgTPCNSigmaElMax", {2.f}, "max TPC nsigma e for inclusion"};
  Configurable<float> cfgTPCNSigmaPiMin{"cfgTPCNSigmaPiMin", {-3.f}, "min TPC nsigma pion for exclusion"};
  Configurable<float> cfgTPCNSigmaPiMax{"cfgTPCNSigmaPiMax", {3.5f}, "max TPC nsigma pion for exclusion"};
  Configurable<float> cfgTPCNSigmaPrMin{"cfgTPCNSigmaPrMin", {-3.f}, "min TPC nsigma proton for exclusion"};
  Configurable<float> cfgTPCNSigmaPrMax{"cfgTPCNSigmaPrMax", {4.f}, "max TPC nsigma proton for exclusion"};
  Configurable<float> cfgTPCNSigmaKaMin{"cfgTPCNSigmaKaMin", {-3.f}, "min TPC nsigma kaon for exclusion"};
  Configurable<float> cfgTPCNSigmaKaMax{"cfgTPCNSigmaKaMax", {4.f}, "max TPC nsigma kaon for exclusion"};
  Configurable<float> cfgTOFNSigmaElMin{"cfgTOFNSigmaElMin", {-3.f}, "min TOF nsigma e for inclusion"};
  Configurable<float> cfgTOFNSigmaElMax{"cfgTOFNSigmaElMax", {3.f}, "max TOF nsigma e for inclusion"};

  Filter ptFilter = (aod::track::pt > cfgPtMin) && (aod::track::pt < cfgPtMax);
  Filter etafilter = nabs(aod::track::eta) < cfgEtaMax;
  Filter dcaxyfilter = nabs(aod::track::dcaXY) < cfgDCAxyMax;
  Filter dcazfilter = nabs(aod::track::dcaZ) < cfgDCAzMax;

  Filter tpcEl = ((aod::pidtpc::tpcNSigmaEl > cfgTPCNSigmaElMin) && (aod::pidtpc::tpcNSigmaEl < cfgTPCNSigmaElMax)) || ((aod::pidtof::tofNSigmaEl > cfgTOFNSigmaElMin) && (aod::pidtof::tofNSigmaEl < cfgTOFNSigmaElMax));
  Filter tpcPi_rej = (aod::pidtpc::tpcNSigmaPi < cfgTPCNSigmaPiMin) || (aod::pidtpc::tpcNSigmaPi > cfgTPCNSigmaPiMax);
  Filter tpcKa_rej = (aod::pidtpc::tpcNSigmaKa < cfgTPCNSigmaKaMin) || (aod::pidtpc::tpcNSigmaKa > cfgTPCNSigmaPrMax);
  Filter tpcPr_rej = (aod::pidtpc::tpcNSigmaPr < cfgTPCNSigmaPrMin) || (aod::pidtpc::tpcNSigmaPr > cfgTPCNSigmaPrMax);

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  double bz{0.}; // magnetic field
  int runNumber{0};

  HistogramRegistry mHistManager{"tpcElIdHistograms"};

  void init(InitContext const&)
  {
    LOG(info) << "Initializing ee mass spectrum via TPC electron identification analysis task ...";

    std::vector<double> momentum_binning = {0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0,
                                            1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,
                                            4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 10.};
    const AxisSpec
      axisCounter{1, 0, +1, ""},
      axisVTrackX{400, -5., 5., "track vertex x (cm)", "track vertex x (cm)"},
      axisVTrackY{400, -5., 5., "track vertex y (cm)", "track vertex y (cm)"},
      axisVTrackZ{400, -20., 20., "track vertex z (cm)", "track vertex z (cm)"},
      axisMassSpectrum{4000, 0, 4, "M (GeV/c^{2})", "Mass e^{+}e^{-} (GeV/c^{2})"},
      axisTPC{1000, 0, 200, "TPC signal (dE/dx)"},
      axisPt{momentum_binning, "p_{T} (GeV/c)"},
      axisPtBig{2000, 0, 20, "p_{T} (GeV/c)"},
      axisEta{600, -3., 3., "#eta"};

    mHistManager.add("eventCounter", "eventCounter", kTH1F, {axisCounter});
    mHistManager.add("hTPCspectra_isElectronRej", "isElectron with rejection | TPC dE/dx spectra", HistType::kTH2F, {axisPt, axisTPC});

    mHistManager.add("h_TPCee_MS_mp", "Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp} (from TPC candidates) vs pair pt", HistType::kTH2F, {axisMassSpectrum, axisPt});
    mHistManager.add("h_TPCee_MS_mm", "Mass e^{-}e^{-} vs momentum e^{-}e^{-} (from TPC candidates) vs pair pt", HistType::kTH2F, {axisMassSpectrum, axisPt});
    mHistManager.add("h_TPCee_MS_pp", "Mass e^{+}e^{+} vs momentum e^{+}e^{+} (from TPC candidates) vs pair pt", HistType::kTH2F, {axisMassSpectrum, axisPt});

    mHistManager.add("h_TPCee_MS_mp_phosRange", "Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp} (from TPC candidates) vs pair pt with one e in phos acceptance range", HistType::kTH2F, {axisMassSpectrum, axisPt});
    mHistManager.add("h_TPCee_MS_mm_phosRange", "Mass e^{-}e^{-} vs momentum e^{-}e^{-} (from TPC candidates) vs pair pt with one e in phos acceptance range", HistType::kTH2F, {axisMassSpectrum, axisPt});
    mHistManager.add("h_TPCee_MS_pp_phosRange", "Mass e^{+}e^{+} vs momentum e^{+}e^{+} (from TPC candidates) vs pair pt with one e in phos acceptance range", HistType::kTH2F, {axisMassSpectrum, axisPt});

    mHistManager.add("h_TPCee_MS_mp_phosRange_kTVXinPHOS", "Mass e^{#pm}e^{#mp} vs momentum e^{#pm}e^{#mp} (from TPC candidates) vs pair pt with one e in phos acceptance range", HistType::kTH2F, {axisMassSpectrum, axisPt});
    mHistManager.add("h_TPCee_MS_mm_phosRange_kTVXinPHOS", "Mass e^{-}e^{-} vs momentum e^{-}e^{-} (from TPC candidates) vs pair pt with one e in phos acceptance range", HistType::kTH2F, {axisMassSpectrum, axisPt});
    mHistManager.add("h_TPCee_MS_pp_phosRange_kTVXinPHOS", "Mass e^{+}e^{+} vs momentum e^{+}e^{+} (from TPC candidates) vs pair pt with one e in phos acceptance range", HistType::kTH2F, {axisMassSpectrum, axisPt});

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
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision,
               soa::Filtered<myTracks> const& tracks,
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
    if (fabs(collision.posZ()) > 10.f)
      return;

    for (auto& [track1, track2] : combinations(CombinationsStrictlyUpperIndexPolicy(tracks, tracks))) {
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

      TLorentzVector P1, P2;
      P1.SetPxPyPzE(track1.px(), track1.py(), track1.pz(), track1.energy(0));
      P2.SetPxPyPzE(track2.px(), track2.py(), track2.pz(), track2.energy(0));

      bool inPhosEtaRange1 = fabs(track1.eta()) < 0.12;
      bool inPhosEtaRange2 = fabs(track2.eta()) < 0.12;
      bool inPhosPhiRange1 = (track1.phi() * TMath::RadToDeg() > 250 && track1.phi() * TMath::RadToDeg() < 320);
      bool inPhosPhiRange2 = (track2.phi() * TMath::RadToDeg() > 250 && track2.phi() * TMath::RadToDeg() < 320);
      bool inPhosRange = (inPhosEtaRange1 && inPhosPhiRange1) || (inPhosEtaRange2 && inPhosPhiRange2);
      bool posTrack = (track1.sign() > 0 && bz > 0) || (track1.sign() < 0 && bz < 0);

      if (track1.sign() == track2.sign()) {
        if (posTrack) {
          mHistManager.fill(HIST("h_TPCee_MS_pp"), (P1 + P2).M(), (P1 + P2).Pt());
          if (inPhosRange) {
            mHistManager.fill(HIST("h_TPCee_MS_pp_phosRange"), (P1 + P2).M(), (P1 + P2).Pt());
            if (collision.alias_bit(kTVXinPHOS))
              mHistManager.fill(HIST("h_TPCee_MS_pp_phosRange_kTVXinPHOS"), (P1 + P2).M(), (P1 + P2).Pt());
          }
        } else {
          mHistManager.fill(HIST("h_TPCee_MS_mm"), (P1 + P2).M(), (P1 + P2).Pt());
          if (inPhosRange) {
            mHistManager.fill(HIST("h_TPCee_MS_mm_phosRange"), (P1 + P2).M(), (P1 + P2).Pt());
            if (collision.alias_bit(kTVXinPHOS))
              mHistManager.fill(HIST("h_TPCee_MS_mm_phosRange_kTVXinPHOS"), (P1 + P2).M(), (P1 + P2).Pt());
          }
        }
      } else {
        mHistManager.fill(HIST("h_TPCee_MS_mp"), (P1 + P2).M(), (P1 + P2).Pt());
        if (inPhosRange) {
          mHistManager.fill(HIST("h_TPCee_MS_mp_phosRange"), (P1 + P2).M(), (P1 + P2).Pt());
          if (collision.alias_bit(kTVXinPHOS))
            mHistManager.fill(HIST("h_TPCee_MS_mp_phosRange_kTVXinPHOS"), (P1 + P2).M(), (P1 + P2).Pt());
        }
      }
    }

    for (auto const& track1 : tracks) {
      mHistManager.fill(HIST("hTrackPt"), track1.pt());
      mHistManager.fill(HIST("hTrackEta"), track1.eta());
      mHistManager.fill(HIST("hTrackVX"), track1.x());
      mHistManager.fill(HIST("hTrackVY"), track1.y());
      mHistManager.fill(HIST("hTrackVZ"), track1.z());

      if (!track1.has_collision() || !track1.hasTPC())
        continue;
      if (track1.itsChi2NCl() > cfgITSchi2Max || track1.tpcChi2NCl() > cfgTPCchi2Max)
        continue;
      if (track1.itsNCls() < cfgITSnclsMin || track1.itsNCls() > cfgITSnclsMax || !((track1.itsClusterMap() & uint8_t(1)) > 0))
        continue;
      if (track1.tpcNClsFound() < cfgTPCnclsMin || track1.tpcNClsFound() > cfgTPCnclsMax)
        continue;
      if (track1.tpcNClsCrossedRows() < cfgTPCnclsCRMin || track1.tpcNClsCrossedRows() > cfgTPCnclsCRMax)
        continue;

      mHistManager.fill(HIST("hTPCspectra_isElectronRej"), track1.pt(), track1.tpcSignal());

      mHistManager.fill(HIST("hTrackPt_Cut"), track1.pt());
      mHistManager.fill(HIST("hTrackEta_Cut"), track1.eta());
      mHistManager.fill(HIST("hTrackVX_Cut"), track1.x());
      mHistManager.fill(HIST("hTrackVY_Cut"), track1.y());
      mHistManager.fill(HIST("hTrackVZ_Cut"), track1.z());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{
    adaptAnalysisTask<phosElId>(cfgc),
    adaptAnalysisTask<massSpectra>(cfgc),
    adaptAnalysisTask<tpcElIdMassSpectrum>(cfgc)};
  return workflow;
}
