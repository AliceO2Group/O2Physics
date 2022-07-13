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
/// \author Mattia Faggin <mattia.faggin@cern.ch>, Padova University and INFN
///
/// Event selection: o2-analysis-timestamp --aod-file AO2D.root -b | o2-analysis-event-selection -b | ---> not working with Run2 converted data/MC
/// Track selection: o2-analysis-trackextension | o2-analysis-trackselection | ---> add --isRun3 1 with Run 3 data/MC (then global track selection works)
/// PID: o2-analysis-pid-tpc-full | o2-analysis-pid-tof-full
/// Working configuration (2021 Oct 20th): o2-analysis-trackextension -b --aod-file ./AO2D.root | o2-analysis-trackselection -b --isRun3 1 | o2-analysis-pid-tpc-full -b | o2-analysis-pid-tof-full -b | o2-analysis-pp-qa-impact-parameter -b

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/DCA.h"
#include "Common/Core/trackUtilities.h" // for propagation to primary vertex

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/PIDResponse.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonUtils/NameConf.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/TrackSelection.h"
#include "DetectorsVertexing/PVertexer.h"
#include "ReconstructionDataFormats/Vertex.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPObject.h"

#include "iostream"
#include "vector"
#include "set"

using namespace o2::framework;
using namespace o2::framework::expressions;

// void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
//{
//   ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, false, {"Fill MC histograms."}};
//   workflowOptions.push_back(optionDoMC);
// }

#include "Framework/runDataProcessing.h"

/// QA task for impact parameter distribution monitoring
struct QaImpactPar {

  /// Input parameters
  Configurable<bool> fDebug{"fDebug", false, "Debug flag enabling outputs"};
  ConfigurableAxis binningImpPar{"binningImpPar", {200, -500.f, 500.f}, "Impact parameter binning"};
  Configurable<bool> keepOnlyPhysPrimary{"keepOnlyPhysPrimary", false, "Consider only phys. primary particles (MC)"};
  // Configurable<int> numberContributorsMin{"numberContributorsMin", 0, "Minimum number of contributors for the primary vertex"};
  Configurable<bool> useTriggerkINT7{"useTriggerkINT7", false, "Use kINT7 trigger"};
  Configurable<bool> usesel8{"usesel8", true, "Use or not the sel8() (T0A & T0C) event selection"};
  Configurable<bool> useIsGlobalTrackFilter{"useIsGlobalTrackFilter", true, "Use or not 'isGlobalTrack' filter"};
  Configurable<float> zVtxMax{"zVtxMax", 10.f, "Maximum value for |z_vtx|"};
  // Configurable<int> keepOnlyGlobalTracks{"keepOnlyGlobalTracks", 1, "Keep only global tracks or not"};
  Configurable<float> ptMin{"ptMin", 0.1f, "Minimum track pt [GeV/c]"};
  Configurable<float> nSigmaTPCPionMin{"nSigmaTPCPionMin", -99999.f, "Minimum nSigma value in TPC, pion hypothesis"};
  Configurable<float> nSigmaTPCPionMax{"nSigmaTPCPionMax", 99999.f, "Maximum nSigma value in TPC, pion hypothesis"};
  Configurable<float> nSigmaTPCKaonMin{"nSigmaTPCKaonMin", -99999.f, "Minimum nSigma value in TPC, kaon hypothesis"};
  Configurable<float> nSigmaTPCKaonMax{"nSigmaTPCKaonMax", 99999.f, "Maximum nSigma value in TPC, kaon hypothesis"};
  Configurable<float> nSigmaTPCProtonMin{"nSigmaTPCProtonMin", -99999.f, "Minimum nSigma value in TPC, proton hypothesis"};
  Configurable<float> nSigmaTPCProtonMax{"nSigmaTPCProtonMax", 99999.f, "Maximum nSigma value in TPC, proton hypothesis"};
  Configurable<float> nSigmaTOFPionMin{"nSigmaTOFPionMin", -99999.f, "Minimum nSigma value in TOF, pion hypothesis"};
  Configurable<float> nSigmaTOFPionMax{"nSigmaTOFPionMax", 99999.f, "Maximum nSigma value in TOF, pion hypothesis"};
  Configurable<float> nSigmaTOFKaonMin{"nSigmaTOFKaonMin", -99999.f, "Minimum nSigma value in TOF, kaon hypothesis"};
  Configurable<float> nSigmaTOFKaonMax{"nSigmaTOFKaonMax", 99999.f, "Maximum nSigma value in TOF, kaon hypothesis"};
  Configurable<float> nSigmaTOFProtonMin{"nSigmaTOFProtonMin", -99999.f, "Minimum nSigma value in TOF, proton hypothesis"};
  Configurable<float> nSigmaTOFProtonMax{"nSigmaTOFProtonMax", 99999.f, "Maximum nSigma value in TOF, proton hypothesis"};
  Configurable<bool> doPVrefit{"doPVrefit", true, "Do PV refit"};
  Configurable<int> nBins_DeltaX_PVrefit{"nBins_DeltaX_PVrefit", 1000, "Number of bins of DeltaX for PV refit"};
  Configurable<int> nBins_DeltaY_PVrefit{"nBins_DeltaY_PVrefit", 1000, "Number of bins of DeltaY for PV refit"};
  Configurable<int> nBins_DeltaZ_PVrefit{"nBins_DeltaZ_PVrefit", 1000, "Number of bins of DeltaZ for PV refit"};
  Configurable<float> minDeltaX_PVrefit{"minDeltaX_PVrefit", -0.5, "Min. DeltaX value for PV refit (cm)"};
  Configurable<float> maxDeltaX_PVrefit{"maxDeltaX_PVrefit", 0.5, "Max. DeltaX value for PV refit (cm)"};
  Configurable<float> minDeltaY_PVrefit{"minDeltaY_PVrefit", -0.5, "Min. DeltaY value for PV refit (cm)"};
  Configurable<float> maxDeltaY_PVrefit{"maxDeltaY_PVrefit", 0.5, "Max. DeltaY value for PV refit (cm)"};
  Configurable<float> minDeltaZ_PVrefit{"minDeltaZ_PVrefit", -0.5, "Min. DeltaZ value for PV refit (cm)"};
  Configurable<float> maxDeltaZ_PVrefit{"maxDeltaZ_PVrefit", 0.5, "Max. DeltaZ value for PV refit (cm)"};
  Configurable<uint16_t> minPVcontrib{"minPVcontrib", 0, "Minimum number of PV contributors"};
  Configurable<uint16_t> maxPVcontrib{"maxPVcontrib", 10000, "Maximum number of PV contributors"};
  Configurable<bool> removeDiamondConstraint{"removeDiamondConstraint", true, "Remove the diamond constraint for the PV refit"};
  Configurable<bool> keepAllTracksPVrefit{"keepAllTracksPVrefit", false, "Keep all tracks for PV refit (for debug)"};
  Configurable<bool> use_customITSHitMap{"use_customITSHitMap", false, "Use custom ITS hitmap selection"};
  Configurable<int> customITShitmap{"customITShitmap", 0, "Custom ITS hitmap (consider the binary representation)"};
  Configurable<int> n_customMinITShits{"n_customMinITShits", 0, "Minimum number of layers crossed by a track among those in \"customITShitmap\""};
  Configurable<bool> custom_forceITSTPCmatching{"custom_forceITSTPCmatching", false, "Consider or not only ITS-TPC macthed tracks when using custom ITS hitmap"};

  /// Custom cut selection objects
  TrackSelection selector_ITShitmap;

  /// Selections with Filter (from o2::framework::expressions)
  // Primary vertex |z_vtx|<XXX cm
  Filter collisionZVtxFilter = (nabs(o2::aod::collision::posZ) < zVtxMax);
  Filter collisionNumContribPV = (minPVcontrib <= o2::aod::collision::numContrib) && (o2::aod::collision::numContrib < maxPVcontrib);
  // Global tracks
  // with Run 3 Reco/MC enable '--isRun3 1' option
  Filter globalTrackFilter = (useIsGlobalTrackFilter.node() == false) || (requireGlobalTrackInFilter()); /// filterbit 4 track selections + tight DCA cuts
  // Pt selection
  Filter ptMinFilter = o2::aod::track::pt > ptMin;

  /// Histogram registry (from o2::framework)
  HistogramRegistry histograms{"HistogramsImpParQA"};
  bool isPIDPionApplied = ((nSigmaTPCPionMin > -10.001 && nSigmaTPCPionMax < 10.001) || (nSigmaTOFPionMin > -10.001 && nSigmaTOFPionMax < 10.001));
  bool isPIDKaonApplied = ((nSigmaTPCKaonMin > -10.001 && nSigmaTPCKaonMax < 10.001) || (nSigmaTOFKaonMin > -10.001 && nSigmaTOFKaonMax < 10.001));
  bool isPIDProtonApplied = ((nSigmaTPCProtonMin > -10.001 && nSigmaTPCProtonMax < 10.001) || (nSigmaTOFProtonMin > -10.001 && nSigmaTOFProtonMax < 10.001));

  // Needed for PV refitting
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  const char* ccdbpath_lut = "GLO/Param/MatLUT";
  const char* ccdbpath_geo = "GLO/Config/GeometryAligned";
  const char* ccdbpath_grp = "GLO/GRP/GRP";
  const char* ccdburl = "http://alice-ccdb.cern.ch";
  // o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int mRunNumber;

  /////////////////////////////////////////////////////////////
  ///                   Process functions                   ///
  /////////////////////////////////////////////////////////////

  /// Data
  using collisionRecoTable = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
  using trackTable = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::TracksExtra>;
  using trackFullTable = o2::soa::Join<o2::aod::Tracks, o2::aod::TrackSelection, o2::aod::TracksCov, o2::aod::TracksExtra, o2::aod::TracksDCA,
                                       o2::aod::pidTPCFullPi, o2::aod::pidTPCFullKa, o2::aod::pidTPCFullPr,
                                       o2::aod::pidTOFFullPi, o2::aod::pidTOFFullKa, o2::aod::pidTOFFullPr>;
  void processData(o2::soa::Filtered<collisionRecoTable>::iterator& collision,
                   const trackTable& tracksUnfiltered,
                   const o2::soa::Filtered<trackFullTable>& tracks,
                   o2::aod::BCsWithTimestamps const& bcs)
  {
    /// here call the template processReco function
    auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
    processReco<false>(collision, tracksUnfiltered, tracks, 0, bc);
  }
  PROCESS_SWITCH(QaImpactPar, processData, "process data", true);

  /// MC
  using collisionMCRecoTable = o2::soa::Join<collisionRecoTable, o2::aod::McCollisionLabels>;
  using trackMCFullTable = o2::soa::Join<trackFullTable, o2::aod::McTrackLabels>;
  void processMC(o2::soa::Filtered<collisionMCRecoTable>::iterator& collision,
                 trackTable const& tracksUnfiltered,
                 o2::soa::Filtered<trackMCFullTable> const& tracks,
                 const o2::aod::McParticles& mcParticles,
                 const o2::aod::McCollisions&,
                 o2::aod::BCsWithTimestamps const& bcs)
  {
    /// here call the template processReco function
    auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
    processReco<true>(collision, tracksUnfiltered, tracks, mcParticles, bc);
  }
  PROCESS_SWITCH(QaImpactPar, processMC, "process MC", false);

  /// core template process function
  /// template<bool IS_MC, typename C, typename T, typename T_MC>
  /// void processReco(const C& collision, const trackTable& unfilteredTracks, const T& tracks,
  ///                 const T_MC& mcParticles,
  ///                 o2::aod::BCsWithTimestamps const& bcs);

  /////////////////////////////////////////////////////////////

  /// init function - declare and define histograms
  void init(InitContext&)
  {
    // Primary vertex
    const AxisSpec collisionXAxis{100, -20.f, 20.f, "X (cm)"};
    const AxisSpec collisionYAxis{100, -20.f, 20.f, "Y (cm)"};
    const AxisSpec collisionZAxis{100, -20.f, 20.f, "Z (cm)"};
    const AxisSpec collisionXOrigAxis{1000, -20.f, 20.f, "X original PV (cm)"};
    const AxisSpec collisionYOrigAxis{1000, -20.f, 20.f, "Y original PV (cm)"};
    const AxisSpec collisionZOrigAxis{1000, -20.f, 20.f, "Z original PV (cm)"};
    const AxisSpec collisionNumberContributorAxis{1000, 0, 1000, "Number of contributors"};
    const AxisSpec collisionDeltaX_PVrefit{nBins_DeltaX_PVrefit, minDeltaX_PVrefit, maxDeltaX_PVrefit, "#Delta x_{PV} (cm)"};
    const AxisSpec collisionDeltaY_PVrefit{nBins_DeltaY_PVrefit, minDeltaY_PVrefit, maxDeltaY_PVrefit, "#Delta y_{PV} (cm)"};
    const AxisSpec collisionDeltaZ_PVrefit{nBins_DeltaZ_PVrefit, minDeltaZ_PVrefit, maxDeltaZ_PVrefit, "#Delta z_{PV} (cm)"};

    histograms.add("Reco/vertices", "", kTH1D, {{2, 0.5f, 2.5f, ""}});
    histograms.get<TH1>(HIST("Reco/vertices"))->GetXaxis()->SetBinLabel(1, "All PV");
    histograms.get<TH1>(HIST("Reco/vertices"))->GetXaxis()->SetBinLabel(2, "PV refit doable");
    histograms.add("Reco/vertices_perTrack", "", kTH1D, {{3, 0.5f, 3.5f, ""}});
    histograms.get<TH1>(HIST("Reco/vertices_perTrack"))->GetXaxis()->SetBinLabel(1, "All PV");
    histograms.get<TH1>(HIST("Reco/vertices_perTrack"))->GetXaxis()->SetBinLabel(2, "PV refit doable");
    histograms.get<TH1>(HIST("Reco/vertices_perTrack"))->GetXaxis()->SetBinLabel(3, "PV refit #chi^{2}!=-1");
    histograms.add("Reco/vertexZ", "", kTH1D, {collisionZAxis});
    histograms.add("Reco/numberContributors", "", kTH1D, {collisionNumberContributorAxis});
    if (doPVrefit) {
      histograms.add("Reco/nContrib_vs_DeltaX_PVrefit", "", kTH2D, {collisionNumberContributorAxis, collisionDeltaX_PVrefit});
      histograms.add("Reco/nContrib_vs_DeltaY_PVrefit", "", kTH2D, {collisionNumberContributorAxis, collisionDeltaY_PVrefit});
      histograms.add("Reco/nContrib_vs_DeltaZ_PVrefit", "", kTH2D, {collisionNumberContributorAxis, collisionDeltaZ_PVrefit});
      histograms.add("Reco/nContrib_vs_Chi2PVrefit", "", kTH2D, {collisionNumberContributorAxis, {102, -1.5, 100.5, "#chi^{2} PV refit"}});
      histograms.add("Reco/X_PVrefitChi2minus1", "PV refit with #chi^{2}==-1", kTH2D, {collisionXAxis, collisionXOrigAxis});
      histograms.add("Reco/Y_PVrefitChi2minus1", "PV refit with #chi^{2}==-1", kTH2D, {collisionYAxis, collisionYOrigAxis});
      histograms.add("Reco/Z_PVrefitChi2minus1", "PV refit with #chi^{2}==-1", kTH2D, {collisionZAxis, collisionZOrigAxis});
      histograms.add("Reco/nContrib_PVrefitNotDoable", "N. contributors for PV refit not doable", kTH1D, {collisionNumberContributorAxis});
      histograms.add("Reco/nContrib_PVrefitChi2minus1", "N. contributors orginal PV for PV refit #chi^{2}==-1", kTH1D, {collisionNumberContributorAxis});
    }

    // Needed for PV refitting
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbpath_lut));
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(ccdbpath_geo);
    }
    mRunNumber = 0;

    /// Custom cut selection objects
    std::set<uint8_t> set_customITShitmap; // = {};
    if (use_customITSHitMap) {
      for (int index_ITSlayer = 0; index_ITSlayer < 7; index_ITSlayer++) {
        if ((customITShitmap & (1 << index_ITSlayer)) > 0) {
          set_customITShitmap.insert(static_cast<uint8_t>(index_ITSlayer));
        }
      }
      LOG(info) << "### customITShitmap: " << customITShitmap;
      LOG(info) << "### n_customMinITShits: " << n_customMinITShits;
      LOG(info) << "### set_customITShitmap.size(): " << set_customITShitmap.size();
      LOG(info) << "### Custom ITS hitmap checked: ";
      for (std::set<uint8_t>::iterator it = set_customITShitmap.begin(); it != set_customITShitmap.end(); it++) {
        LOG(info) << "Layer " << (int)(*it) << " ";
      }
      LOG(info) << "############";

      selector_ITShitmap.SetRequireHitsInITSLayers(n_customMinITShits, set_customITShitmap);
    }

    // tracks
    const AxisSpec trackPtAxis{100, 0.f, 10.f, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec trackEtaAxis{40, -2.f, 2.f, "#it{#eta}"};
    const AxisSpec trackPhiAxis{24, 0.f, TMath::TwoPi(), "#varphi"};
    const AxisSpec trackImpParRPhiAxis{binningImpPar, "#it{d}_{r#it{#varphi}} (#mum)"};
    const AxisSpec trackImpParZAxis{binningImpPar, "#it{d}_{z} (#mum)"};
    const AxisSpec trackNSigmaTPCPionAxis{20, -10.f, 10.f, "Number of #sigma TPC #pi^{#pm}"};
    const AxisSpec trackNSigmaTPCKaonAxis{20, -10.f, 10.f, "Number of #sigma TPC K^{#pm}"};
    const AxisSpec trackNSigmaTPCProtonAxis{20, -10.f, 10.f, "Number of #sigma TPC proton"};
    const AxisSpec trackNSigmaTOFPionAxis{20, -10.f, 10.f, "Number of #sigma TOF #pi^{#pm}"};
    const AxisSpec trackNSigmaTOFKaonAxis{20, -10.f, 10.f, "Number of #sigma TOF K^{#pm}"};
    const AxisSpec trackNSigmaTOFProtonAxis{20, -10.f, 10.f, "Number of #sigma TOF proton"};
    const AxisSpec trackPDGAxis{5, -1.5f, 3.5f, "species (-1: not matched, 0: unknown, 1: pi, 2: K, 3: p)"};

    histograms.add("Reco/pt", "", kTH1D, {trackPtAxis});
    histograms.add("Reco/itsHits", "Number of hits vs ITS layer;layer ITS", kTH2D, {{8, -1.5, 6.5, "ITS layer"}, {8, -0.5, 7.5, "Number of hits"}});
    histograms.add("Reco/refitRun3", "Refit conditions for tracks", kTH1D, {{5, 0.5, 5.5, ""}});
    histograms.get<TH1>(HIST("Reco/refitRun3"))->GetXaxis()->SetBinLabel(1, "hasTPC");
    histograms.get<TH1>(HIST("Reco/refitRun3"))->GetXaxis()->SetBinLabel(2, "hasITS");
    histograms.get<TH1>(HIST("Reco/refitRun3"))->GetXaxis()->SetBinLabel(3, "hasTPC && !hasITS");
    histograms.get<TH1>(HIST("Reco/refitRun3"))->GetXaxis()->SetBinLabel(4, "!hasTPC && hasITS");
    histograms.get<TH1>(HIST("Reco/refitRun3"))->GetXaxis()->SetBinLabel(5, "hasTPC && hasITS");
    histograms.add("Reco/h4ImpPar", "", kTHnSparseD, {trackPtAxis, trackImpParRPhiAxis, trackEtaAxis, trackPhiAxis, trackPDGAxis});
    histograms.add("Reco/h4ImpParZ", "", kTHnSparseD, {trackPtAxis, trackImpParZAxis, trackEtaAxis, trackPhiAxis, trackPDGAxis});
    if (isPIDPionApplied) {
      histograms.add("Reco/h4ImpPar_Pion", "", kTHnSparseD, {trackPtAxis, trackImpParRPhiAxis, trackEtaAxis, trackPhiAxis, trackPDGAxis});
      histograms.add("Reco/h4ImpParZ_Pion", "", kTHnSparseD, {trackPtAxis, trackImpParZAxis, trackEtaAxis, trackPhiAxis, trackPDGAxis});
    }
    if (isPIDKaonApplied) {
      histograms.add("Reco/h4ImpPar_Kaon", "", kTHnSparseD, {trackPtAxis, trackImpParRPhiAxis, trackEtaAxis, trackPhiAxis, trackPDGAxis});
      histograms.add("Reco/h4ImpParZ_Kaon", "", kTHnSparseD, {trackPtAxis, trackImpParZAxis, trackEtaAxis, trackPhiAxis, trackPDGAxis});
    }
    if (isPIDProtonApplied) {
      histograms.add("Reco/h4ImpPar_Proton", "", kTHnSparseD, {trackPtAxis, trackImpParRPhiAxis, trackEtaAxis, trackPhiAxis, trackPDGAxis});
      histograms.add("Reco/h4ImpParZ_Proton", "", kTHnSparseD, {trackPtAxis, trackImpParZAxis, trackEtaAxis, trackPhiAxis, trackPDGAxis});
    }
    histograms.add("Reco/hNSigmaTPCPion", "", kTH2D, {trackPtAxis, trackNSigmaTPCPionAxis});
    histograms.add("Reco/hNSigmaTPCKaon", "", kTH2D, {trackPtAxis, trackNSigmaTPCKaonAxis});
    histograms.add("Reco/hNSigmaTPCProton", "", kTH2D, {trackPtAxis, trackNSigmaTPCProtonAxis});
    histograms.add("Reco/hNSigmaTOFPion", "", kTH2D, {trackPtAxis, trackNSigmaTOFPionAxis});
    histograms.add("Reco/hNSigmaTOFKaon", "", kTH2D, {trackPtAxis, trackNSigmaTOFKaonAxis});
    histograms.add("Reco/hNSigmaTOFProton", "", kTH2D, {trackPtAxis, trackNSigmaTOFProtonAxis});
    histograms.add("Reco/hNSigmaTPCPion_afterPID", "", kTH2D, {trackPtAxis, trackNSigmaTPCPionAxis});
    histograms.add("Reco/hNSigmaTPCKaon_afterPID", "", kTH2D, {trackPtAxis, trackNSigmaTPCKaonAxis});
    histograms.add("Reco/hNSigmaTPCProton_afterPID", "", kTH2D, {trackPtAxis, trackNSigmaTPCProtonAxis});
    histograms.add("Reco/hNSigmaTOFPion_afterPID", "", kTH2D, {trackPtAxis, trackNSigmaTOFPionAxis});
    histograms.add("Reco/hNSigmaTOFKaon_afterPID", "", kTH2D, {trackPtAxis, trackNSigmaTOFKaonAxis});
    histograms.add("Reco/hNSigmaTOFProton_afterPID", "", kTH2D, {trackPtAxis, trackNSigmaTOFProtonAxis});

    histograms.add("MC/vertexZ_MCColl", "", kTH1D, {collisionZAxis});
    histograms.add("MC/ptMC", "", kTH1D, {trackPtAxis});
  }

  /// core template process function
  template <bool IS_MC, typename C, typename T, typename T_MC>
  void processReco(const C& collision, const trackTable& unfilteredTracks, const T& tracks,
                   const T_MC& mcParticles,
                   o2::aod::BCsWithTimestamps::iterator const& bc)
  {
    constexpr float toMicrometers = 10000.f; // Conversion from [cm] to [mum]

    auto PDGtoIndex = [](const int pdg) {
      switch (pdg) {
        case 211: // pion
          return 1;
        case 321: // kaon
          return 2;
        case 2212: // proton
          return 3;
        default: // not identified
          return 0;
      }
    };

    /// trigger selection
    if (useTriggerkINT7) {
      // from Tutorial/src/multiplicityEventTrackSelection.cxx
      if (!collision.alias()[kINT7]) {
        return;
      }
    }
    /// offline event selections
    if (usesel8 && !collision.sel8()) {
      return;
    }

    histograms.fill(HIST("Reco/vertices"), 1);
    histograms.fill(HIST("Reco/vertexZ"), collision.posZ());
    histograms.fill(HIST("Reco/numberContributors"), collision.numContrib());
    if constexpr (IS_MC) {
      if (collision.has_mcCollision()) {
        histograms.fill(HIST("MC/vertexZ_MCColl"), collision.mcCollision().posZ());
      }
    }

    ///////////////////////////////////
    ///       For PV refit          ///
    ///////////////////////////////////
    /// retrieve the tracks contributing to the primary vertex fitting
    std::vector<int64_t> vec_globID_contr = {};
    std::vector<o2::track::TrackParCov> vec_TrkContributos = {};
    if (fDebug) {
      LOG(info) << "\n === New collision";
    }
    const int nTrk = unfilteredTracks.size();
    int nContrib = 0;
    int nNonContrib = 0;
    for (const auto& unfilteredTrack : unfilteredTracks) {
      if (!unfilteredTrack.isPVContributor()) {
        /// the track di not contribute to fit the primary vertex
        nNonContrib++;
        continue;
      }
      vec_globID_contr.push_back(unfilteredTrack.globalIndex());
      vec_TrkContributos.push_back(getTrackParCov(unfilteredTrack));
      nContrib++;
      if (fDebug) {
        LOG(info) << "---> a contributor! stuff saved";
        LOG(info) << "vec_contrib size: " << vec_TrkContributos.size() << ", nContrib: " << nContrib;
      }
    }
    if (fDebug) {
      LOG(info) << "===> nTrk: " << nTrk << ",   nContrib: " << nContrib << ",   nNonContrib: " << nNonContrib;
    }

    if (vec_TrkContributos.size() != collision.numContrib()) {
      LOG(info) << "!!! something wrong in the number of contributor tracks for PV fit !!! " << vec_TrkContributos.size() << " vs. " << collision.numContrib();
      return;
    }

    std::vector<bool> vec_useTrk_PVrefit(vec_globID_contr.size(), true);

    /// Prepare the vertex refitting
    // Get the magnetic field for the Propagator
    o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    // auto bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
    if (mRunNumber != bc.runNumber()) {
      auto grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbpath_grp, bc.timestamp());
      if (grpo != nullptr) {
        o2::base::Propagator::initFieldFromGRP(grpo);
        o2::base::Propagator::Instance()->setMatLUT(lut);
        LOGF(info, "Setting magnetic field to %d kG for run %d from its GRP CCDB object", grpo->getNominalL3Field(), bc.runNumber());
      } else {
        LOGF(fatal, "GRP object is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
      }
      mRunNumber = bc.runNumber();
    }
    // build the VertexBase to initialize the vertexer
    o2::dataformats::VertexBase Pvtx;
    Pvtx.setX(collision.posX());
    Pvtx.setY(collision.posY());
    Pvtx.setZ(collision.posZ());
    Pvtx.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    // configure PVertexer
    o2::vertexing::PVertexer vertexer;
    if (removeDiamondConstraint) {
      o2::conf::ConfigurableParam::updateFromString("pvertexer.useMeanVertexConstraint=false"); // we want to refit w/o MeanVertex constraint
    }
    vertexer.init();
    bool PVrefit_doable = vertexer.prepareVertexRefit(vec_TrkContributos, Pvtx);
    if (!PVrefit_doable) {
      LOG(info) << "Not enough tracks accepted for the refit";
      if (doPVrefit) {
        histograms.fill(HIST("Reco/nContrib_PVrefitNotDoable"), collision.numContrib());
      }
    } else {
      histograms.fill(HIST("Reco/vertices"), 2);
    }

    if (fDebug) {
      LOG(info) << "prepareVertexRefit = " << PVrefit_doable << " Ncontrib= " << vec_TrkContributos.size() << " Ntracks= " << collision.numContrib() << " Vtx= " << Pvtx.asString();
    }
    ///////////////////////////////////
    ///////////////////////////////////

    /// loop over tracks
    float pt = -999.f;
    float impParRPhi = -999.f;
    float impParZ = -999.f;
    float tpcNSigmaPion = -999.f;
    float tpcNSigmaKaon = -999.f;
    float tpcNSigmaProton = -999.f;
    float tofNSigmaPion = -999.f;
    float tofNSigmaKaon = -999.f;
    float tofNSigmaProton = -999.f;
    int ntr = tracks.size();
    int cnt = 0;
    for (const auto& track : tracks) {

      /// Specific MC selections
      int pdgIndex = -1;
      if constexpr (IS_MC) {
        if (keepOnlyPhysPrimary) {
          /// we want only physical primary particles
          if (!track.has_mcParticle()) {
            continue;
          }
          auto particle = track.mcParticle();
          if (keepOnlyPhysPrimary && particle.isPhysicalPrimary()) {
            continue;
          }
          pdgIndex = PDGtoIndex(std::abs(particle.pdgCode()));
          histograms.fill(HIST("MC/ptMC"), particle.pt());
        } else {
          if (track.has_mcParticle()) {
            auto particle = track.mcParticle();
            pdgIndex = PDGtoIndex(std::abs(particle.pdgCode()));
            histograms.fill(HIST("MC/ptMC"), particle.pt());
          }
        }
      }

      /// Using the Filter instead
      /// if ((keepOnlyGlobalTracks) && (!track.isGlobalTrack())) {
      ///  /// not a global track (FB 4 with tight DCA cuts)
      ///  continue;
      ///}

      /// apply custom ITS hitmap selections, if asked
      if (use_customITSHitMap && !selector_ITShitmap.IsSelected(track, TrackSelection::TrackCuts::kITSHits)) {
        /// skip this track and go on, because it does not satisfy the ITS hit requirements
        continue;
      }
      if (use_customITSHitMap && custom_forceITSTPCmatching && (!track.hasITS() || !track.hasTPC())) {
        // if (use_customITSHitMap && custom_forceITSTPCmatching && track.hasITS()) {  ///ATTEMPT: REMOVE TRACKS WITH ITS
        /// skip this track because it is not global (no matching ITS-TPC)
        continue;
      }
      int itsNhits = 0;
      for (unsigned int i = 0; i < 7; i++) {
        if (track.itsClusterMap() & (1 << i)) {
          itsNhits += 1;
        }
      }
      bool trkHasITS = false;
      for (unsigned int i = 0; i < 7; i++) {
        if (track.itsClusterMap() & (1 << i)) {
          trkHasITS = true;
          histograms.fill(HIST("Reco/itsHits"), i, itsNhits);
        }
      }
      if (!trkHasITS) {
        histograms.fill(HIST("Reco/itsHits"), -1, itsNhits);
      }
      /// Test refit conditions
      if (track.hasTPC()) { // hasTPC
        histograms.fill(HIST("Reco/refitRun3"), 1);
      }
      if (track.hasITS()) { // hasITS
        histograms.fill(HIST("Reco/refitRun3"), 2);
      }
      if (track.hasTPC() && !track.hasITS()) { // hasTPC && !hasITS
        histograms.fill(HIST("Reco/refitRun3"), 3);
      }
      if (!track.hasTPC() && track.hasITS()) { // !hasTPC && hasITS
        histograms.fill(HIST("Reco/refitRun3"), 4);
      }
      if (track.hasTPC() && track.hasITS()) { // hasTPC && hasITS
        histograms.fill(HIST("Reco/refitRun3"), 5);
      }

      pt = track.pt();
      tpcNSigmaPion = track.tpcNSigmaPi();
      tpcNSigmaKaon = track.tpcNSigmaKa();
      tpcNSigmaProton = track.tpcNSigmaPr();
      tofNSigmaPion = track.tofNSigmaPi();
      tofNSigmaKaon = track.tofNSigmaKa();
      tofNSigmaProton = track.tofNSigmaPr();

      histograms.fill(HIST("Reco/pt"), pt);
      histograms.fill(HIST("Reco/hNSigmaTPCPion"), pt, tpcNSigmaPion);
      histograms.fill(HIST("Reco/hNSigmaTPCKaon"), pt, tpcNSigmaKaon);
      histograms.fill(HIST("Reco/hNSigmaTPCProton"), pt, tpcNSigmaProton);
      histograms.fill(HIST("Reco/hNSigmaTOFPion"), pt, tofNSigmaPion);
      histograms.fill(HIST("Reco/hNSigmaTOFKaon"), pt, tofNSigmaKaon);
      histograms.fill(HIST("Reco/hNSigmaTOFProton"), pt, tofNSigmaProton);

      histograms.fill(HIST("Reco/vertices_perTrack"), 1);
      if (PVrefit_doable) {
        histograms.fill(HIST("Reco/vertices_perTrack"), 2);
      }
      /// PV refitting, if the tracks contributed to this at the beginning
      o2::dataformats::VertexBase PVbase_recalculated;
      bool recalc_imppar = false;
      if (doPVrefit && PVrefit_doable) {
        recalc_imppar = true;
        auto it_trk = std::find(vec_globID_contr.begin(), vec_globID_contr.end(), track.globalIndex()); /// track global index
        // if( it_trk==vec_globID_contr.end() ) {
        //   /// not found: this track did not contribute to the initial PV fitting
        //   continue;
        // }
        if (it_trk != vec_globID_contr.end()) {
          /// this track contributed to the PV fit: let's do the refit without it
          const int entry = std::distance(vec_globID_contr.begin(), it_trk);
          if (!keepAllTracksPVrefit) {
            vec_useTrk_PVrefit[entry] = false; /// remove the track from the PV refitting
          }
          auto Pvtx_refitted = vertexer.refitVertex(vec_useTrk_PVrefit, Pvtx); // vertex refit
          if (fDebug) {
            LOG(info) << "refit " << cnt << "/" << ntr << " result = " << Pvtx_refitted.asString();
          }

          if (Pvtx_refitted.getChi2() < 0) {
            LOG(info) << "---> Refitted vertex has bad chi2 = " << Pvtx_refitted.getChi2();
            histograms.fill(HIST("Reco/X_PVrefitChi2minus1"), Pvtx_refitted.getX(), collision.posX());
            histograms.fill(HIST("Reco/Y_PVrefitChi2minus1"), Pvtx_refitted.getY(), collision.posY());
            histograms.fill(HIST("Reco/Z_PVrefitChi2minus1"), Pvtx_refitted.getZ(), collision.posZ());
            histograms.fill(HIST("Reco/nContrib_PVrefitChi2minus1"), collision.numContrib());
            recalc_imppar = false;
          } else {
            histograms.fill(HIST("Reco/vertices_perTrack"), 3);
          }
          // histograms.fill(HIST("Reco/nContrib_vs_Chi2PVrefit"), /*Pvtx_refitted.getNContributors()*/collision.numContrib()-1, Pvtx_refitted.getChi2());
          histograms.fill(HIST("Reco/nContrib_vs_Chi2PVrefit"), vec_useTrk_PVrefit.size() - 1, Pvtx_refitted.getChi2());

          vec_useTrk_PVrefit[entry] = true; /// restore the track for the next PV refitting

          if (recalc_imppar) {
            // fill the histograms for refitted PV with good Chi2
            const double DeltaX = Pvtx.getX() - Pvtx_refitted.getX();
            const double DeltaY = Pvtx.getY() - Pvtx_refitted.getY();
            const double DeltaZ = Pvtx.getZ() - Pvtx_refitted.getZ();
            histograms.fill(HIST("Reco/nContrib_vs_DeltaX_PVrefit"), collision.numContrib(), DeltaX);
            histograms.fill(HIST("Reco/nContrib_vs_DeltaY_PVrefit"), collision.numContrib(), DeltaY);
            histograms.fill(HIST("Reco/nContrib_vs_DeltaZ_PVrefit"), collision.numContrib(), DeltaZ);

            // fill the newly calculated PV
            PVbase_recalculated.setX(Pvtx_refitted.getX());
            PVbase_recalculated.setY(Pvtx_refitted.getY());
            PVbase_recalculated.setZ(Pvtx_refitted.getZ());
            PVbase_recalculated.setCov(Pvtx_refitted.getSigmaX2(), Pvtx_refitted.getSigmaXY(), Pvtx_refitted.getSigmaY2(), Pvtx_refitted.getSigmaXZ(), Pvtx_refitted.getSigmaYZ(), Pvtx_refitted.getSigmaZ2());
          }

          cnt++;
        }
      } /// end 'if (doPVrefit && PVrefit_doable)'

      /// impact parameter to the PV
      // value calculated wrt global PV (not recalculated) ---> coming from trackextension workflow
      impParRPhi = toMicrometers * track.dcaXY(); // dca.getY();
      impParZ = toMicrometers * track.dcaZ();     // dca.getY();
      // updated value after PV recalculation
      if (recalc_imppar) {
        auto trackPar = getTrackPar(track);
        o2::gpu::gpustd::array<float, 2> dcaInfo{-999., -999.};
        if (o2::base::Propagator::Instance()->propagateToDCABxByBz({PVbase_recalculated.getX(), PVbase_recalculated.getY(), PVbase_recalculated.getZ()}, trackPar, 2.f, matCorr, &dcaInfo)) {
          // if (o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, matCorr, &dcaInfo)) {
          impParRPhi = dcaInfo[0] * toMicrometers;
          impParZ = dcaInfo[1] * toMicrometers;
        }
      }

      /// all tracks
      histograms.fill(HIST("Reco/h4ImpPar"), pt, impParRPhi, track.eta(), track.phi(), pdgIndex);
      histograms.fill(HIST("Reco/h4ImpParZ"), pt, impParZ, track.eta(), track.phi(), pdgIndex);

      if (isPIDPionApplied && nSigmaTPCPionMin < tpcNSigmaPion && tpcNSigmaPion < nSigmaTPCPionMax && nSigmaTOFPionMin < tofNSigmaPion && tofNSigmaPion < nSigmaTOFPionMax) {
        /// PID selected pions
        histograms.fill(HIST("Reco/h4ImpPar_Pion"), pt, impParRPhi, track.eta(), track.phi(), pdgIndex);
        histograms.fill(HIST("Reco/h4ImpParZ_Pion"), pt, impParZ, track.eta(), track.phi(), pdgIndex);
        histograms.fill(HIST("Reco/hNSigmaTPCPion_afterPID"), pt, tpcNSigmaPion);
        histograms.fill(HIST("Reco/hNSigmaTOFPion_afterPID"), pt, tofNSigmaPion);
      }
      if (isPIDKaonApplied && nSigmaTPCKaonMin < tpcNSigmaKaon && tpcNSigmaKaon < nSigmaTPCKaonMax && nSigmaTOFKaonMin < tofNSigmaKaon && tofNSigmaKaon < nSigmaTOFKaonMax) {
        /// PID selected kaons
        histograms.fill(HIST("Reco/h4ImpPar_Kaon"), pt, impParRPhi, track.eta(), track.phi(), pdgIndex);
        histograms.fill(HIST("Reco/h4ImpParZ_Kaon"), pt, impParZ, track.eta(), track.phi(), pdgIndex);
        histograms.fill(HIST("Reco/hNSigmaTPCKaon_afterPID"), pt, tpcNSigmaKaon);
        histograms.fill(HIST("Reco/hNSigmaTOFKaon_afterPID"), pt, tofNSigmaKaon);
      }
      if (isPIDProtonApplied && nSigmaTPCProtonMin < tpcNSigmaProton && tpcNSigmaProton < nSigmaTPCProtonMax && nSigmaTOFProtonMin < tofNSigmaProton && tofNSigmaProton < nSigmaTOFProtonMax) {
        /// PID selected Protons
        histograms.fill(HIST("Reco/h4ImpPar_Proton"), pt, impParRPhi, track.eta(), track.phi(), pdgIndex);
        histograms.fill(HIST("Reco/h4ImpParZ_Proton"), pt, impParZ, track.eta(), track.phi(), pdgIndex);
        histograms.fill(HIST("Reco/hNSigmaTPCProton_afterPID"), pt, tpcNSigmaProton);
        histograms.fill(HIST("Reco/hNSigmaTOFProton_afterPID"), pt, tofNSigmaProton);
      }
    }
  } /// end processReco
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec w{
    adaptAnalysisTask<QaImpactPar>(cfgc, TaskName{"qa-impact-par"})};
  return w;
}
