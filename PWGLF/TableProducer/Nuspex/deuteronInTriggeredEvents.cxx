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
///
/// \file deuteronInTriggeredEvents.cxx
///
/// \brief (Anti-)nuclei spectra analysis task in jet-triggered events
/// \author Cristian Moscatelli (cristian.moscatelli@cern.ch)
///
/// Based on PWGLF/TableProducer/Nuspex/nucleiSpectra.cxx
/// \since 05/2026
//
// ================
// Executable + dependencies:
//
// data (run3):
// o2-analysis-lf-deuteron-in-triggered-events, o2-analysis-event-selection-service
// o2-analysis-propagationservice, o2-analysis-trackselection, o2-analysis-track-extra-v002-converter
// o2-analysis-pid-tof-merge, o2-analysis-pid-tpc-service, o2-analysis-ft0-corrected-table
//
//  mc:
// same as Data and o2-analysis-mccollision-converter

#include "PWGLF/DataModel/LFSlimNucleiTables.h"
#include "PWGLF/Utils/inelGt.h"
//
#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/JetUtilities.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/TableProducer/PID/pidTOFBase.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/BetheBlochAleph.h>
#include <ReconstructionDataFormats/Track.h>

#include <Math/Vector4D.h>
#include <TMCProcess.h>
#include <TPDGCode.h> // for PDG codes
#include <TRandom3.h>

#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/PseudoJet.hh>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct NucleusCandidate {
  int globalIndex = 0;
  int collTrackIndex = 0;
  float pt = 0.f;
  float eta = 0.f;
  float phi = 0.f;
  float tpcInnerParam = 0.f;
  float beta = 0.f;
  float zVertex = 0.f;
  int nContrib = 0;
  float DCAxy = 0.f;
  float DCAz = 0.f;
  float TPCsignal = 0.f;
  float ITSchi2 = 0.f;
  float TPCchi2 = 0.f;
  float TOFchi2 = 0.f;
  std::array<float, 5> nSigmaTPC{};
  std::array<float, 5> tofMasses{};
  bool fillTree = false;
  bool fillDCAHist = false;
  bool correctPV = false;
  bool isSecondary = false;
  bool fromWeakDecay = false;
  uint16_t flags = 0;
  uint8_t TPCfindableCls = 0;
  uint8_t TPCcrossedRows = 0;
  uint8_t ITSclsMap = 0;
  uint8_t TPCnCls = 0;
  uint8_t TPCnClsShared = 0;
  uint8_t ITSnCls = 0;
  uint32_t clusterSizesITS = 0;
};

namespace nuclei
{
constexpr double bbMomScalingDefault[5][2]{
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.}};
constexpr double betheBlochDefault[5][6]{
  {-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09},
  {-136.71, 0.441, 0.2269, 1.347, 0.8035, 0.09},
  {-239.99, 1.155, 1.099, 1.137, 1.006, 0.09},
  {-321.34, 0.6539, 1.591, 0.8225, 2.363, 0.09},
  {-586.66, 1.859, 4.435, 0.282, 3.201, 0.09}};
constexpr double nSigmaTPCdefault[5][2]{
  {-5., 5.},
  {-5., 5.},
  {-5., 5.},
  {-5., 5.},
  {-5., 5.}};
constexpr double DCAcutDefault[5][2]{
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.},
  {1., 1.}};
constexpr int TreeConfigDefault[5][2]{
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0}};
constexpr int DCAHistDefault[5][2]{
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0},
  {0, 0}};
constexpr double DownscalingDefault[5][1]{
  {1.},
  {1.},
  {1.},
  {1.},
  {1.}};
constexpr int species{5};                                                                    // Number of analysed species
constexpr int codes[5]{2212, 1000010020, 1000010030, 1000020030, 1000020040};                // PGD code of these particles
constexpr float charges[5]{1.f, 1.f, 1.f, 2.f, 2.f};                                         // Their charge
constexpr float masses[5]{MassProton, MassDeuteron, MassTriton, MassHelium3, MassAlpha};     // Their masses
static const std::vector<std::string> matter{"M", "A"};                                      // Type of particles (matter or antimatter)
static const std::vector<std::string> pidName{"TPC", "TOF"};                                 // Type of PID
static const std::vector<int> hfMothCodes{511, 521, 531, 541, 5122};                         // b-mesons + Lambda_b
static const std::vector<std::string> names{"proton", "deuteron", "triton", "He3", "alpha"}; // Particles name
static const std::vector<std::string> treeConfigNames{"Filter trees", "Use TOF selection"};
static const std::vector<std::string> DCAConfigNames{"Save DCA hist", "Matter/Antimatter"};
static const std::vector<std::string> nSigmaConfigName{"nsigma_min", "nsigma_max"};                   // label for nsigmaTPC selection
static const std::vector<std::string> nDCAConfigName{"max DCAxy", "max DCAz"};                        // label for DCA selection
static const std::vector<std::string> DownscalingConfigName{"Fraction of kept candidates"};           // Fraction of candidates to be kept
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"}; // BB paramteters
static const std::vector<std::string> binnedVariableNames{"DCAxy", "DCAz", "TPCnsigma", "TOFnsigma", "TOFmass"};
static const std::vector<std::string> chargeLabelNames{"Positive", "Negative"}; // Particles charge (pos or neg)

float pidCuts[2][5][2];
// This is a set of 2x5x2 3D histograms
// described by smart pointer.
// Keep in mind: first index: TPC or TOF. Second index: type of particle. Third index: charge.
std::shared_ptr<TH2> hNsigma[2][5][2];
std::shared_ptr<TH2> hTOFmass[5][2];
std::shared_ptr<TH1> hGenNuclei[5][2];
std::shared_ptr<TH2> hMomRes[5][2];
std::shared_ptr<TH3> hNsigmaEta[2][5][2];
std::shared_ptr<TH3> hTOFmassEta[5][2];
std::shared_ptr<TH2> hDCAxy[2][5][2];
std::shared_ptr<TH2> hDCAz[2][5][2];
std::shared_ptr<TH2> hGloTOFtracks[2];
std::shared_ptr<TH2> hDeltaP[2][5];
std::shared_ptr<THnSparse> hDCAHists[2][5];
o2::base::MatLayerCylSet* lut = nullptr; // describe detector material

std::vector<NucleusCandidate> candidates;

enum evSel {
  kTVX = 0,
  kTFBorder, // Here we can substitute with sel8
  kITSROFborder,
  kZvtx,
  kNoSameBunchPileup,
  kIsGoodZvtxFT0vsPV,
  kIsGoodITSLayersAll,
  kIsEPtriggered,
  kINELgt0,
  kIsJetTriggered, // Check: evSel for event with at least one jet over pT-threshold
  kNevSels
};

static const std::vector<std::string> eventSelectionTitle{"Event selections"};
static const std::vector<std::string> eventSelectionLabels{"TVX", "TF border", "ITS ROF border", "Z vtx", "No same-bunch pile-up", "kIsGoodZvtxFT0vsPV", "isGoodITSLayersAll", "isEPtriggered", "IsINELgt0", "IsJetTriggered"};

constexpr int EvSelDefault[10][1]{
  {1},
  {1},
  {1}, // Default: sel8
  {0},
  {0},
  {0},
  {0},
  {0},
  {0},  // INEL > 0
  {1}}; // Triggered on jets

enum evGenSel : uint8_t {
  kGenIsINELgt0 = 1 << 0,
  kGenIsJetTriggered = 1 << 1,
  kGenHasRecoEv = 1 << 2
};

enum triggerListName {
  fChJetLowPt = 0,
  fChJetHighPt = 1
};
} // namespace nuclei

struct DeuteronInTriggeredEvents {
  enum {
    kProton = BIT(0),
    kDeuteron = BIT(1),
    kTriton = BIT(2),
    kHe3 = BIT(3),
    kHe4 = BIT(4),
    kHasTOF = BIT(5),
    kHasTRD = BIT(6),
    kIsAmbiguous = BIT(7), /// just a placeholder now
    kITSrof = BIT(8),
    kIsPhysicalPrimary = BIT(9), /// MC flags starting from the second half of the short
    kIsSecondaryFromMaterial = BIT(10),
    kIsSecondaryFromWeakDecay = BIT(11) /// the last 4 bits are reserved for the PID in tracking
  };

  Produces<o2::aod::NucleiTable> nucleiTable;                       // For data
  Produces<o2::aod::NucleiTableMCExtension> nucleiTableMCExtension; // For MC analysis
  Produces<o2::aod::GenEventMCSel> GenEventMCSel;                   // For MC reco events
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdgDB; // For INELgt0 gen MC selection
  Zorro zorro;                                 // Definition of Zorro: helpful for skimmed data
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  Configurable<bool> cfgCompensatePIDinTracking{"cfgCompensatePIDinTracking", false, "If true, divide tpcInnerParam by the electric charge"};

  struct : o2::framework::ConfigurableGroup {
    std::string prefix{"cfgTrackCut"};
    Configurable<LabeledArray<double>> dcaMax{"dcaMax", {nuclei::DCAcutDefault[0], 5, 2, nuclei::names, nuclei::nDCAConfigName}, "Max DCAxy and DCAz for light nuclei"};
    Configurable<float> etaMax{"etaMax", 0.8f, "Max Eta for tracks"};
    Configurable<int> itsNClusMin{"itsNClusMin", 5, "Minimum number of ITS clusters"};
    Configurable<float> itsChi2ClusMax{"itsChi2ClusMax", 36.f, "Max ITS Chi2 per cluster"};
    Configurable<float> rapidityMax{"rapidityMax", 1.f, "Maximum rapidity for tracks"};
    Configurable<float> rapidityMin{"rapidityMin", -1.f, "Minimum rapidity for tracks"};
    Configurable<bool> rapidityToggle{"rapidityToggle", false, "If true, use rapidity cuts"};
    Configurable<float> tpcChi2ClusMax{"tpcChi2ClusMax", 4.f, "Max TPC Chi2 per cluster"};
    Configurable<int> tpcNCrossedRowsMin{"tpcNCrossedRowsMin", 70, "Minimum number of TPC crossed rows"};
    Configurable<float> tpcNCrossedRowsOverFindableMin{"tpcNCrossedRowsOverFindableMin", 0.8f, "Minimum ratio of crossed rows over findable clusters"};
    Configurable<int> tpcNClsMin{"tpcNClsMin", 80, "Minimum number of TPC clusters"};
    Configurable<float> tpcRigidityMin{"tpcRigidityMin", 0.5f, "Minimum TPC rigidity for tracks"};
    Configurable<LabeledArray<double>> tpcNSigmaMax{"tpcNSigmaMax", {nuclei::nSigmaTPCdefault[0], 5, 2, nuclei::names, nuclei::nSigmaConfigName}, "TPC nsigma selection for light nuclei"};
  } cfgTrackCut;

  Configurable<float> cfgCutPtMinTree{"cfgCutPtMinTree", 0.2f, "Minimum track transverse momentum for tree saving"};
  Configurable<float> cfgCutPtMaxTree{"cfgCutPtMaxTree", 15.0f, "Maximum track transverse momentum for tree saving"};

  // Event selections
  Configurable<LabeledArray<int>> cfgEventSelections{"cfgEventSelections", {nuclei::EvSelDefault[0], 10, 1, nuclei::eventSelectionLabels, nuclei::eventSelectionTitle}, "Event selections"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};

  Configurable<LabeledArray<double>> cfgMomentumScalingBetheBloch{"cfgMomentumScalingBetheBloch", {nuclei::bbMomScalingDefault[0], 5, 2, nuclei::names, nuclei::chargeLabelNames}, "TPC Bethe-Bloch momentum scaling for light nuclei"};
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {nuclei::betheBlochDefault[0], 5, 6, nuclei::names, nuclei::betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgDownscaling{"cfgDownscaling", {nuclei::DownscalingDefault[0], 5, 1, nuclei::names, nuclei::DownscalingConfigName}, "Fraction of kept candidates for light nuclei"};
  Configurable<LabeledArray<int>> cfgTreeConfig{"cfgTreeConfig", {nuclei::TreeConfigDefault[0], 5, 2, nuclei::names, nuclei::treeConfigNames}, "Filtered trees configuration"};
  Configurable<int> cfgFillGenSecondaries{"cfgFillGenSecondaries", 0, "Fill generated secondaries (0: no, 1: only weak decays, 2: all of them)"};
  Configurable<LabeledArray<int>> cfgDCAHists{"cfgDCAHists", {nuclei::DCAHistDefault[0], 5, 2, nuclei::names, nuclei::DCAConfigNames}, "DCA hist configuration"};

  Configurable<double> cfgNsigmaTPCcutDCAhists{"cfgNsigmaTPCcutDCAhists", 3., "TPC nsigma cut for DCA hists"};
  Configurable<double> cfgDeltaTOFmassCutDCAhists{"cfgDeltaTOFmassCutDCAhists", 0.2, "Delta TOF mass cut for DCA hists"};
  ConfigurableAxis cfgDCAxyBinsProtons{"cfgDCAxyBinsProtons", {1500, -1.5f, 1.5f}, "DCAxy binning for Protons"};
  ConfigurableAxis cfgDCAxyBinsDeuterons{"cfgDCAxyBinsDeuterons", {1500, -1.5f, 1.5f}, "DCAxy binning for Deuterons"};
  ConfigurableAxis cfgDCAxyBinsTritons{"cfgDCAxyBinsTritons", {1500, -1.5f, 1.5f}, "DCAxy binning for Tritons"};
  ConfigurableAxis cfgDCAxyBinsHe3{"cfgDCAxyBinsHe3", {1500, -1.5f, 1.5f}, "DCAxy binning for He3"};
  ConfigurableAxis cfgDCAxyBinsAlpha{"cfgDCAxyBinsAlpha", {1500, -1.5f, 1.5f}, "DCAxy binning for Alpha"};

  ConfigurableAxis cfgPtBinsProtons{"cfgPtBinsProtons", {100, 0., 10.}, "Pt binning for Protons"};
  ConfigurableAxis cfgPtBinsDeuterons{"cfgPtBinsDeuterons", {100, 0., 10.}, "Pt binning for Deuterons"};
  ConfigurableAxis cfgPtBinsTritons{"cfgPtBinsTritons", {100, 0., 10.}, "Pt binning for Tritons"};
  ConfigurableAxis cfgPtBinsHe3{"cfgPtBinsHe3", {100, 0., 10.}, "Pt binning for He3"};
  ConfigurableAxis cfgPtBinsAlpha{"cfgPtBinsAlpha", {100, 0., 10.}, "Pt binning for Alpha"};

  Configurable<float> cfgCMrapidity{"cfgCMrapidity", 0.f, "Rapidity of the center of mass (only for p-Pb)"};
  ConfigurableAxis cfgMomResBins{"cfgMomResBins", {200, -1., 1.}, "Momentum resolution binning"};
  ConfigurableAxis cfgNsigmaTPCbins{"cfgNsigmaTPCbins", {100, -5., 5.}, "nsigma_TPC binning"};
  ConfigurableAxis cfgNsigmaTOFbins{"cfgNsigmaTOFbins", {100, -5., 5.}, "nsigma_TOF binning"};
  ConfigurableAxis cfgTOFmassBins{"cfgTOFmassBins", {200, -5., 5.}, "TOF mass binning"};
  ConfigurableAxis cfgNITSClusBins{"cfgNITSClusBins", {3, 4.5, 7.5}, "N ITS clusters binning"};
  ConfigurableAxis cfgNTPCClusBins{"cfgNTPCClusBins", {3, 89.5, 159.5}, "N TPC clusters binning"};

  // Configurable for working with skimmed data
  Configurable<bool> cfgApplyMCEvSel{"cfgApplyMCEvSel", false, "If true, apply jet-trigger selection on gen events"};
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "Skimmed dataset processing"};
  Configurable<int> cfgTriggerList{"cfgTriggerList", 0, "0 : Low jet-pT thr, 1 : High jet-pT thr"};

  // Configurable for jet identification
  Configurable<double> cfgRJet{"cfgRJet", 0.6, "R_jet"};
  Configurable<double> cfgEtaJet{"cfgEtaJet", 0.9, "Jet acceptance"};

  // running variables for track tuner
  o2::dataformats::DCA mDcaInfoCov;
  o2::track::TrackParametrizationWithError<float> mTrackParCov;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

  // CCDB options
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrLUT), "Type of material correction"};
  Configurable<std::string> cfgCCDBurl{"cfgCCDBurl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> cfgZorroCCDBpath{"cfgZorroCCDBpath", "EventFiltering/Zorro/", "path to the zorro ccdb objects"};
  int mRunNumber = 0;
  float mBz = 0.f;

  // Utility object for jet background subtraction methods
  JetBkgSubUtils backgroundSub;

  using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksCovIU, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime, aod::Tracks, aod::TrackSelection, aod::TrackSelectionExtension, aod::TracksDCA>;

  HistogramRegistry spectra{"spectra", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  double computeAbsoDecL(aod::McParticles::iterator particle)
  {
    if (!particle.has_daughters())
      return -1.f;

    float mothVtx[3]{particle.vx(), particle.vy(), particle.vz()};
    float dauVtx[3]{0.f, 0.f, 0.f};
    auto daughters = particle.daughters_as<aod::McParticles>();
    for (const auto& dau : daughters) {
      if (std::abs(dau.pdgCode()) != PDG_t::kGamma && std::abs(dau.pdgCode()) != PDG_t::kElectron) {
        dauVtx[0] = dau.vx();
        dauVtx[1] = dau.vy();
        dauVtx[2] = dau.vz();
        break;
      }
    }
    return std::hypot(mothVtx[0] - dauVtx[0], mothVtx[1] - dauVtx[1], mothVtx[2] - dauVtx[2]);
  }

  template <class Tcoll>
  bool eventSelectionWithHisto(Tcoll& collision)
  {
    spectra.fill(HIST("hEventSelections"), 0);

    if (cfgEventSelections->get(nuclei::evSel::kTVX) && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    spectra.fill(HIST("hEventSelections"), nuclei::evSel::kTVX + 1);

    if (cfgEventSelections->get(nuclei::evSel::kTFBorder) && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    spectra.fill(HIST("hEventSelections"), nuclei::evSel::kTFBorder + 1);

    if (cfgEventSelections->get(nuclei::evSel::kITSROFborder) && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    spectra.fill(HIST("hEventSelections"), nuclei::evSel::kITSROFborder + 1);

    if (cfgEventSelections->get(nuclei::evSel::kZvtx) && std::abs(collision.posZ()) > cfgCutVertex) {
      return false;
    }
    spectra.fill(HIST("hEventSelections"), nuclei::evSel::kZvtx + 1);

    if (cfgEventSelections->get(nuclei::evSel::kNoSameBunchPileup) && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    spectra.fill(HIST("hEventSelections"), nuclei::evSel::kNoSameBunchPileup + 1);

    if (cfgEventSelections->get(nuclei::evSel::kIsGoodZvtxFT0vsPV) && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    spectra.fill(HIST("hEventSelections"), nuclei::evSel::kIsGoodZvtxFT0vsPV + 1);

    if (cfgEventSelections->get(nuclei::evSel::kIsGoodITSLayersAll) && !collision.selection_bit(aod::evsel::kIsGoodITSLayersAll)) {
      return false;
    }
    spectra.fill(HIST("hEventSelections"), nuclei::evSel::kIsGoodITSLayersAll + 1);

    if constexpr (
      requires {
        collision.triggereventep();
      }) {
      if (cfgEventSelections->get(nuclei::evSel::kIsEPtriggered) && !collision.triggereventep()) {
        return false;
      }
      spectra.fill(HIST("hEventSelections"), nuclei::evSel::kIsEPtriggered + 1);
    }

    if (cfgEventSelections->get(nuclei::evSel::kINELgt0) && !collision.selection_bit(aod::kINELgtZERO)) {
      return false;
    }
    spectra.fill(HIST("hEventSelections"), nuclei::evSel::kINELgt0 + 1);

    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);

    gRandom->SetSeed(bc.timestamp());

    if (cfgSkimmedProcessing) {
      bool isTriggered = zorro.isSelected(bc.globalBC()); /// Just let Zorro do the accounting
      if (!isTriggered)
        return false;
    }

    spectra.fill(HIST("hEventSelections"), nuclei::evSel::kIsJetTriggered + 1);

    return true;
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), (cfgTriggerList.value == nuclei::fChJetLowPt) ? "fJetChLowPt" : (cfgTriggerList.value == nuclei::fChJetHighPt) ? "fJetChHighPt"
                                                                                                                                                                                  : throw std::runtime_error("Invalid TriggerList value"));
      zorro.populateHistRegistry(spectra, bc.runNumber());
    }
    auto timestamp = bc.timestamp();
    mRunNumber = bc.runNumber();

    o2::parameters::GRPMagField* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    o2::base::Propagator::Instance()->setMatLUT(nuclei::lut);
    mBz = static_cast<float>(grpmag->getNominalL3Field());
    LOGF(info, "Retrieved GRP for timestamp %ull (%i) with magnetic field of %1.2f kZG", timestamp, mRunNumber, mBz);
  }

  void init(o2::framework::InitContext&)
  {
    zorroSummary.setObject(zorro.getZorroSummary());
    zorro.setBaseCCDBPath(cfgZorroCCDBpath.value);
    ccdb->setURL(cfgCCDBurl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    const AxisSpec nSigmaAxes[2]{{cfgNsigmaTPCbins, "n#sigma_{TPC}"}, {cfgNsigmaTOFbins, "n#sigma_{TOF}"}};
    const AxisSpec tofMassAxis{cfgTOFmassBins, "TOF mass - PDG mass"};
    const AxisSpec ptResAxis{cfgMomResBins, "#Delta#it{p}_{T}/#it{p}_{T}"};
    const AxisSpec nITSClusAxis{cfgNITSClusBins, "N ITS clusters"};
    const AxisSpec nTPCClusAxis{cfgNTPCClusBins, "N TPC clusters"};
    const AxisSpec hasTRDAxis{2, -0.5, 1.5, "Has TRD"};
    const AxisSpec correctPVAxis{2, -0.5, 1.5, "Correct PV"};
    const AxisSpec isSecondaryAxis{2, -0.5, 1.5, "Is secondary"};
    const AxisSpec fromWeakDecayAxis{2, -0.5, 1.5, "From weak decay"};

    const AxisSpec ptAxes[5]{
      {cfgPtBinsProtons, "#it{p}_{T} (GeV/#it{c})"},
      {cfgPtBinsDeuterons, "#it{p}_{T} (GeV/#it{c})"},
      {cfgPtBinsTritons, "#it{p}_{T} (GeV/#it{c})"},
      {cfgPtBinsHe3, "#it{p}_{T} (GeV/#it{c})"},
      {cfgPtBinsAlpha, "#it{p}_{T} (GeV/#it{c})"}};
    const AxisSpec dcaxyAxes[5]{
      {cfgDCAxyBinsProtons, "DCA_{xy} (cm)"},
      {cfgDCAxyBinsDeuterons, "DCA_{xy} (cm)"},
      {cfgDCAxyBinsTritons, "DCA_{xy} (cm)"},
      {cfgDCAxyBinsHe3, "DCA_{xy} (cm)"},
      {cfgDCAxyBinsAlpha, "DCA_{xy} (cm)"}};
    const AxisSpec dcazAxes[5]{
      {cfgDCAxyBinsProtons, "DCA_{z} (cm)"},
      {cfgDCAxyBinsDeuterons, "DCA_{z} (cm)"},
      {cfgDCAxyBinsTritons, "DCA_{z} (cm)"},
      {cfgDCAxyBinsHe3, "DCA_{z} (cm)"},
      {cfgDCAxyBinsAlpha, "DCA_{z} (cm)"}};
    const AxisSpec etaAxis{40, -1., 1., "#eta"};

    spectra.add("hEventSelections", "hEventSelections", {HistType::kTH1D, {{nuclei::evSel::kNevSels + 1, -0.5f, static_cast<float>(nuclei::evSel::kNevSels) + 0.5f}}});
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(1, "all");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kTVX + 2, "TVX");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kTFBorder + 2, "TFborder");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kITSROFborder + 2, "ITSROFborder");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kZvtx + 2, "Zvtx");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kNoSameBunchPileup + 2, "kNoSameBunchPileup");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kIsGoodZvtxFT0vsPV + 2, "isGoodZvtxFT0vsPV");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kIsGoodITSLayersAll + 2, "IsGoodITSLayersAll");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kIsEPtriggered + 2, "IsEPtriggered");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kINELgt0 + 2, "IsINELgt0");
    spectra.get<TH1>(HIST("hEventSelections"))->GetXaxis()->SetBinLabel(nuclei::evSel::kIsJetTriggered + 2, "IsJetTriggered");

    // Distribution of z-vertex of selected events
    spectra.add("hRecVtxZData", "collision z position", HistType::kTH1F, {{200, -20., +20., "z position (cm)"}});
    // Distribution of z-vertex of generated events
    if (doprocessMC) {
      spectra.add("hGenVtxZ", " generated collision z position", HistType::kTH1F, {{200, -20., +20., "z position (cm)"}});
    }
    spectra.add("hTpcSignalData", "Specific energy loss", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    spectra.add("hTpcSignalDataSelected", "Specific energy loss for selected particles", HistType::kTH2F, {{600, -6., 6., "#it{p} (GeV/#it{c})"}, {1400, 0, 1400, "d#it{E} / d#it{X} (a. u.)"}});
    spectra.add("hTofSignalData", "TOF beta", HistType::kTH2F, {{500, 0., 5., "#it{p} (GeV/#it{c})"}, {750, 0, 1.5, "TOF #beta"}});

    for (unsigned int iC{0}; iC < nuclei::matter.size(); ++iC) {
      nuclei::hGloTOFtracks[iC] = spectra.add<TH2>(fmt::format("hTPCTOFtracks{}", nuclei::matter[iC]).data(), fmt::format("Global vs TOF matched {} tracks in a collision", nuclei::chargeLabelNames[iC]).data(), HistType::kTH2D, {{300, -0.5, 300.5, "Number of global tracks"}, {300, -0.5, 300.5, "Number of TOF matched tracks"}});

      for (int iS{0}; iS < nuclei::species; ++iS) {
        nuclei::hNsigma[0][iS][iC] = spectra.add<TH2>(fmt::format("h{}nsigma{}_{}", nuclei::pidName[0], nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("n#sigma_{{}} {} {}", nuclei::pidName[0], nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTH2D, {ptAxes[iS], nSigmaAxes[0]});
        nuclei::hNsigmaEta[0][iS][iC] = spectra.add<TH3>(fmt::format("h{}nsigmaEta{}_{}", nuclei::pidName[0], nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("n#sigma_{{}} {} {} vs #eta", nuclei::pidName[0], nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTH3D, {etaAxis, ptAxes[iS], nSigmaAxes[0]});

        for (unsigned int iPID{0}; iPID < nuclei::matter.size(); ++iPID) {
          nuclei::hDCAxy[iPID][iS][iC] = spectra.add<TH2>(fmt::format("hDCAxy{}_{}_{}", nuclei::pidName[iPID], nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("DCAxy {} {} {}", nuclei::pidName[iPID], nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTH2D, {ptAxes[iS], dcaxyAxes[iS]});
          nuclei::hDCAz[iPID][iS][iC] = spectra.add<TH2>(fmt::format("hDCAz{}_{}_{}", nuclei::pidName[iPID], nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("DCAz {} {} {}", nuclei::pidName[iPID], nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTH2D, {ptAxes[iS], dcazAxes[iS]});
        }
        nuclei::hTOFmass[iS][iC] = spectra.add<TH2>(fmt::format("h{}TOFmass{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("TOF mass - {}  PDG mass", nuclei::names[iS]).data(), HistType::kTH2D, {ptAxes[iS], tofMassAxis});
        nuclei::hTOFmassEta[iS][iC] = spectra.add<TH3>(fmt::format("h{}TOFmassEta{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("TOF mass - {}  PDG mass", nuclei::names[iS]).data(), HistType::kTH3D, {etaAxis, ptAxes[iS], tofMassAxis});
        nuclei::hDeltaP[iC][iS] = spectra.add<TH2>(fmt::format("hDeltaP{}_{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("#Delta#it{{p}}/#it{{p}} {} {}", nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTH2D, {{232, 0.2, 6., "#it{p} (GeV/#it{{c}})"}, {200, -1.0, 1.0, "(#it{p}_{IU} - #it{p}_{TPC}) / #it{p}"}});
        if (doprocessMC) {
          nuclei::hMomRes[iS][iC] = spectra.add<TH2>(fmt::format("h{}MomRes{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("Momentum resolution {}", nuclei::names[iS]).data(), HistType::kTH2D, {ptAxes[iS], ptResAxis});
          nuclei::hGenNuclei[iS][iC] = spectra.add<TH1>(fmt::format("h{}Gen{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("Generated {}", nuclei::names[iS]).data(), HistType::kTH1D, {ptAxes[iS]});
          nuclei::hDCAHists[iC][iS] = spectra.add<THnSparse>(fmt::format("hDCAHists{}_{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("DCA histograms {} {}", nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTHnSparseF, {ptAxes[iS], dcaxyAxes[iS], dcazAxes[iS], nSigmaAxes[0], tofMassAxis, nITSClusAxis, nTPCClusAxis, correctPVAxis, isSecondaryAxis, fromWeakDecayAxis});
        } else {
          nuclei::hDCAHists[iC][iS] = spectra.add<THnSparse>(fmt::format("hDCAHists{}_{}", nuclei::matter[iC], nuclei::names[iS]).data(), fmt::format("DCA histograms {} {}", nuclei::matter[iC], nuclei::names[iS]).data(), HistType::kTHnSparseF, {ptAxes[iS], dcaxyAxes[iS], dcazAxes[iS], nSigmaAxes[0], tofMassAxis, nITSClusAxis, nTPCClusAxis});
        }
      }
    }

    // Setting the nsigma interval in which the the candidate is kept
    // The first index is set to 0 because it is for TPC
    for (int iS{0}; iS < nuclei::species; ++iS) {
      for (unsigned int iMax{0}; iMax < nuclei::pidName.size(); ++iMax) {
        nuclei::pidCuts[0][iS][iMax] = cfgTrackCut.tpcNSigmaMax->get(iS, iMax);
      }
    }

    nuclei::lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
  }

  // Find hit on ITS layer
  template <typename TrackIts>
  bool hasHitITS(const TrackIts& track, int layer)
  {
    int ibit = layer - 1;
    return (track.itsClusterMap() & (1 << ibit));
  }

  template <typename Ttrks>
  bool isJetTriggered(Ttrks const& tracks, nuclei::triggerListName triggerCondition)
  {
    // Defining trigger condition
    double jetPtThreshold(0.0);

    switch (triggerCondition) {
      case nuclei::fChJetLowPt:
        jetPtThreshold = 30.0;
        break;
      case nuclei::fChJetHighPt:
        jetPtThreshold = 55.0;
        break;
      default:
        return false; // Non-valid trigger
    }

    // Loop over tracks
    int id(-1);
    std::vector<fastjet::PseudoJet> fjParticles;
    fjParticles.clear();
    for (auto const& track : tracks) {
      id++;
      if (!passedTrackSelectionForJetReconstruction(track))
        continue;

      // 4-momentum representation of a particle
      fastjet::PseudoJet fourMomentum(track.px(), track.py(), track.pz(), track.energy(MassPionCharged));
      fourMomentum.set_user_index(id);
      fjParticles.emplace_back(fourMomentum);
    }

    if (fjParticles.empty())
      return false; // reject 0-particles events

    // Cluster particles using the the anti-kT algorithm
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, cfgRJet);                    // Defining the algorithm to cluster, and the jet radius
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0)); // Activate area evaluation, and set area evaluation method
    fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);                        // Declare the will of applying clustering algorithm with area evaluation to the selected candidates
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());

    if (jets.empty())
      return false; // skip events with no reconstructed jets
    auto [rhoPerp, rhoMPerp] = jetutilities::estimateRhoPerpCone(fjParticles, jets[0], cfgRJet);

    // Analyse only leading jet
    const auto& jet = jets[0];
    if (std::fabs(jet.eta()) > cfgEtaJet)
      return false; // Jet must be fully contained in the acceptance

    // Jet pt must be larger than threshold
    auto jetForSub = jet;
    fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
    if (jetMinusBkg.pt() < jetPtThreshold)
      return false;

    return true;
  }

  template <typename McParts>
  bool isMCJetTriggered(McParts const& McParticles, aod::McParticles const& particlesMC, nuclei::triggerListName triggerCondition)
  {
    // Defining trigger condition
    double jetPtThreshold(0.0);

    switch (triggerCondition) {
      case nuclei::fChJetLowPt:
        jetPtThreshold = 30.0;
        break;
      case nuclei::fChJetHighPt:
        jetPtThreshold = 55.0;
        break;
      default:
        return false; // Non-valid trigger
    }

    std::vector<fastjet::PseudoJet> fjParticles;
    fjParticles.clear();

    for (auto const& particle : McParticles) {

      // MC particles per collision
      if (!isPhysicalPrimaryOrFromHF(particle, particlesMC))
        continue;

      double minPtParticle = 0.1;
      if (std::fabs(particle.eta()) > cfgEtaJet || particle.pt() < minPtParticle)
        continue;

      // Build 4-momentum assuming charged pion mass
      static constexpr float kMassPionChargedSquared = o2::constants::physics::MassPionCharged * o2::constants::physics::MassPionCharged;
      const double energy = std::sqrt(particle.p() * particle.p() + kMassPionChargedSquared);
      fastjet::PseudoJet fourMomentum(particle.px(), particle.py(), particle.pz(), energy);
      fourMomentum.set_user_index(particle.pdgCode());
      fjParticles.emplace_back(fourMomentum);
    }

    if (fjParticles.empty())
      return false; // reject 0-particles events

    // Cluster particles using the the anti-kT algorithm
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, cfgRJet);                    // Defining the algorithm to cluster, and the jet radius
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0)); // Activate area evaluation, and set area evaluation method
    fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);                        // Declare the will of applying clustering algorithm with area evaluation to the selected candidates
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());

    if (jets.empty())
      return false; // skip events with no reconstructed jets
    auto [rhoPerp, rhoMPerp] = jetutilities::estimateRhoPerpCone(fjParticles, jets[0], cfgRJet);

    // Analyse only leading jet
    const auto& jet = jets[0];
    if (std::fabs(jet.eta()) > cfgEtaJet)
      return false; // Jet must be fully contained in the acceptance

    // Jet pt must be larger than threshold
    auto jetForSub = jet;
    fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
    if (jetMinusBkg.pt() < jetPtThreshold)
      return false;

    return true;
  }

  // Check if particle is a physical primary or a decay product of a heavy-flavor hadron
  bool isPhysicalPrimaryOrFromHF(aod::McParticle const& particle, aod::McParticles const& mcParticles)
  {
    // Keep only pi, K, p, e, mu
    int pdg = std::abs(particle.pdgCode());
    if (!(pdg == PDG_t::kPiPlus || pdg == PDG_t::kKPlus || pdg == PDG_t::kProton || pdg == PDG_t::kElectron || pdg == PDG_t::kMuonMinus))
      return false;

    // Constants for identifying heavy-flavor (charm and bottom) content from PDG codes
    static constexpr int kCharmQuark = 4;
    static constexpr int kBottomQuark = 5;
    static constexpr int hundreds = 100;
    static constexpr int thousands = 1000;

    // Check if particle is from heavy-flavor decay
    bool fromHF = false;
    if (particle.has_mothers()) {
      auto mother = mcParticles.iteratorAt(particle.mothersIds()[0]);
      int motherPdg = std::abs(mother.pdgCode());
      fromHF = (motherPdg / hundreds == kCharmQuark || motherPdg / hundreds == kBottomQuark || motherPdg / thousands == kCharmQuark || motherPdg / thousands == kBottomQuark);
    }

    // Select only physical primary particles or from heavy-flavor
    return (particle.isPhysicalPrimary() || fromHF);
  }

  // Track selection criteria for creating jets
  template <typename JetTrack>
  bool passedTrackSelectionForJetReconstruction(const JetTrack& track)
  {
    static constexpr int MinTpcCr = 70;
    static constexpr double MaxChi2Tpc = 4.0;
    static constexpr double MaxChi2Its = 36.0;
    static constexpr double MinPtTrack = 0.1;
    static constexpr double DcaxyMaxTrackPar0 = 0.0105;
    static constexpr double DcaxyMaxTrackPar1 = 0.035;
    static constexpr double DcaxyMaxTrackPar2 = 1.1;
    static constexpr double DcazMaxTrack = 2.0;

    // General part
    if (std::fabs(track.eta()) > cfgEtaJet)
      return false;
    if (track.pt() < MinPtTrack)
      return false;
    if (std::fabs(track.dcaXY()) > (DcaxyMaxTrackPar0 + DcaxyMaxTrackPar1 / std::pow(track.pt(), DcaxyMaxTrackPar2)))
      return false; // DCAxy cut
    if (std::fabs(track.dcaZ()) > DcazMaxTrack)
      return false; // DCAz cut
    // Part relative to ITS
    if (!track.hasITS())
      return false;
    if ((!hasHitITS(track, 1)) && (!hasHitITS(track, 2)) && (!hasHitITS(track, 3)))
      return false; // Has Inner Barrel hit
    if (track.itsChi2NCl() >= MaxChi2Its)
      return false;
    // Part relative to TPC
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < MinTpcCr)
      return false;
    if (track.tpcChi2NCl() >= MaxChi2Tpc)
      return false;

    return true;
  }

  template <typename Tcoll, typename Ttrks>
  void fillDataInfo(Tcoll const& collision, Ttrks const& tracks)
  {
    spectra.fill(HIST("hRecVtxZData"), collision.posZ());

    const o2::math_utils::Point3D<float> collVtx{collision.posX(), collision.posY(), collision.posZ()};

    const double bgScalings[5][2]{
      {nuclei::charges[0] * cfgMomentumScalingBetheBloch->get(0u, 0u) / nuclei::masses[0], nuclei::charges[0] * cfgMomentumScalingBetheBloch->get(0u, 1u) / nuclei::masses[0]},
      {nuclei::charges[1] * cfgMomentumScalingBetheBloch->get(1u, 0u) / nuclei::masses[1], nuclei::charges[1] * cfgMomentumScalingBetheBloch->get(1u, 1u) / nuclei::masses[1]},
      {nuclei::charges[2] * cfgMomentumScalingBetheBloch->get(2u, 0u) / nuclei::masses[2], nuclei::charges[2] * cfgMomentumScalingBetheBloch->get(2u, 1u) / nuclei::masses[2]},
      {nuclei::charges[3] * cfgMomentumScalingBetheBloch->get(3u, 0u) / nuclei::masses[3], nuclei::charges[3] * cfgMomentumScalingBetheBloch->get(3u, 1u) / nuclei::masses[3]},
      {nuclei::charges[4] * cfgMomentumScalingBetheBloch->get(4u, 0u) / nuclei::masses[4], nuclei::charges[4] * cfgMomentumScalingBetheBloch->get(4u, 1u) / nuclei::masses[4]}};

    int nGloTracks[2]{0, 0}, nTOFTracks[2]{0, 0};

    constexpr float MinMomentumForDeltaP = 0.2f;

    for (const auto& track : tracks) { // start loop over tracks
      if (std::abs(track.eta()) > cfgTrackCut.etaMax ||
          track.tpcInnerParam() < cfgTrackCut.tpcRigidityMin ||
          track.itsNCls() < cfgTrackCut.itsNClusMin ||
          track.tpcNClsFound() < cfgTrackCut.tpcNClsMin ||
          track.tpcNClsCrossedRows() < cfgTrackCut.tpcNCrossedRowsMin ||
          track.tpcNClsCrossedRows() < cfgTrackCut.tpcNCrossedRowsOverFindableMin * track.tpcNClsFindable() ||
          track.tpcChi2NCl() > cfgTrackCut.tpcChi2ClusMax ||
          track.itsChi2NCl() > cfgTrackCut.itsChi2ClusMax) {
        continue;
      }
      // temporary fix: tpcInnerParam() returns the momentum in all the software tags before --- CM: what does it mean
      bool heliumPID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
      float correctedTpcInnerParam = (heliumPID && cfgCompensatePIDinTracking) ? track.tpcInnerParam() / 2 : track.tpcInnerParam();

      spectra.fill(HIST("hTpcSignalData"), correctedTpcInnerParam * track.sign(), track.tpcSignal());
      float nSigma[2][5]{
        {-10., -10., -10., -10., -10.},
        {0.f, 0.f, 0.f, 0.f, 0.f}}; /// then we will calibrate the TOF mass for the He3 and Alpha
      const int iC{track.sign() < 0};

      /// Checking if we have outliers in the TPC-TOF correlation
      nGloTracks[iC]++;
      if (track.hasTOF()) {
        nTOFTracks[iC]++;
      }

      bool selectedTPC[5]{false}, goodToAnalyse{false};
      std::array<float, 5> nSigmaTPC;
      // This part apply a selection on tracks:
      // if the nsigma of the track is less than a certain value, it is accepted
      for (int iS{0}; iS < nuclei::species; ++iS) { // The first index here is set to 0 because it is TPC analysis
        double expBethe{common::BetheBlochAleph(static_cast<double>(correctedTpcInnerParam * bgScalings[iS][iC]), cfgBetheBlochParams->get(iS, 0u), cfgBetheBlochParams->get(iS, 1u), cfgBetheBlochParams->get(iS, 2u), cfgBetheBlochParams->get(iS, 3u), cfgBetheBlochParams->get(iS, 4u))};
        double expSigma{expBethe * cfgBetheBlochParams->get(iS, 5u)};
        nSigma[0][iS] = static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
        nSigmaTPC[iS] = nSigma[0][iS];
        selectedTPC[iS] = (nSigma[0][iS] > nuclei::pidCuts[0][iS][0] && nSigma[0][iS] < nuclei::pidCuts[0][iS][1]);
        goodToAnalyse = goodToAnalyse || selectedTPC[iS];
        if (selectedTPC[iS] && track.p() > MinMomentumForDeltaP) {
          nuclei::hDeltaP[iC][iS]->Fill(track.p(), 1 - correctedTpcInnerParam / track.p()); // Important: see of there is a shift between track momentum and tpc innerparam
        }
      }
      if (!goodToAnalyse) {
        continue;
      }

      mDcaInfoCov.set(999, 999, 999, 999, 999);
      setTrackParCov(track, mTrackParCov);
      mTrackParCov.setPID(track.pidForTracking());
      std::array<float, 2> dcaInfo;
      o2::base::Propagator::Instance()->propagateToDCA(collVtx, mTrackParCov, mBz, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfo);

      float beta{o2::pid::tof::Beta::GetBeta(track)};
      spectra.fill(HIST("hTpcSignalDataSelected"), correctedTpcInnerParam * track.sign(), track.tpcSignal());
      spectra.fill(HIST("hTofSignalData"), correctedTpcInnerParam, beta);
      beta = std::min(1.f - 1.e-6f, std::max(1.e-4f, beta));                       /// sometimes beta > 1 or < 0
      uint16_t flag = static_cast<uint16_t>((track.pidForTracking() & 0xF) << 12); // Store the PID of the track in the 12-15 bit
      std::array<float, 5> tofMasses{-3.f, -3.f, -3.f, -3.f, -3.f};
      bool fillTree{false};
      bool fillDCAHist{false};
      bool correctPV{false};
      bool isSecondary{false};
      bool fromWeakDecay{false};

      if (track.hasTOF()) {
        flag |= kHasTOF;
      }
      if (track.hasTRD()) {
        flag |= kHasTRD;
      }
      if (!collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        flag |= kITSrof;
      }

      for (int iS{0}; iS < nuclei::species; ++iS) {
        bool selectedTOF{false};
        if (std::abs(dcaInfo[1]) > cfgTrackCut.dcaMax->get(iS, 1)) {
          continue;
        }
        ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>> fvector{mTrackParCov.getPt() * nuclei::charges[iS], mTrackParCov.getEta(), mTrackParCov.getPhi(), nuclei::masses[iS]};
        float y{fvector.Rapidity() + cfgCMrapidity};
        for (unsigned int iPID{0}; iPID < nuclei::pidName.size(); ++iPID) { /// 0 TPC, 1 TOF
          if (selectedTPC[iS]) {
            if (iPID && !track.hasTOF()) {
              continue;
            } else if (iPID) {
              selectedTOF = true; /// temporarly skipped
              float charge{nuclei::charges[iS]};
              tofMasses[iS] = correctedTpcInnerParam * charge * std::sqrt(1.f / (beta * beta) - 1.f) - nuclei::masses[iS];
            }
            if (!cfgTrackCut.rapidityToggle || (y > cfgTrackCut.rapidityMin && y < cfgTrackCut.rapidityMax)) {
              if (std::abs(nSigmaTPC[iS]) < cfgNsigmaTPCcutDCAhists && (!iPID || std::abs(tofMasses[iS]) < cfgDeltaTOFmassCutDCAhists)) {
                nuclei::hDCAxy[iPID][iS][iC]->Fill(fvector.pt(), dcaInfo[0]);
                nuclei::hDCAz[iPID][iS][iC]->Fill(fvector.pt(), dcaInfo[1]);
              }
              if (std::abs(dcaInfo[0]) < cfgTrackCut.dcaMax->get(iS, 0u)) {
                if (!iPID) { /// temporary exclusion of the TOF nsigma PID for the He3 and Alpha
                  nuclei::hNsigma[iPID][iS][iC]->Fill(fvector.pt(), nSigma[iPID][iS]);
                  nuclei::hNsigmaEta[iPID][iS][iC]->Fill(fvector.eta(), fvector.pt(), nSigma[iPID][iS]);
                }
                if (iPID) {
                  nuclei::hTOFmass[iS][iC]->Fill(fvector.pt(), tofMasses[iS]);
                  nuclei::hTOFmassEta[iS][iC]->Fill(fvector.eta(), fvector.pt(), tofMasses[iS]);
                }
              }
            }
          }
        }
        if (selectedTPC[iS]) {
          if (cfgTreeConfig->get(iS, 1u) && !selectedTOF) {
            continue;
          }
          !fillTree && cfgTreeConfig->get(iS, 0u) ? fillTree = true : fillTree;
          !fillDCAHist && cfgDCAHists->get(iS, iC) ? fillDCAHist = true : fillDCAHist;
          bool setPartFlag = cfgTreeConfig->get(iS, 0u) || cfgDCAHists->get(iS, iC);
          if (setPartFlag) {
            if (cfgDownscaling->get(iS) < 1. && gRandom->Rndm() > cfgDownscaling->get(iS)) {
              continue;
            }
            flag |= BIT(iS);
          }
        }
      }
      if (flag & (kProton | kDeuteron | kTriton | kHe3 | kHe4) || doprocessMC) { /// ignore PID pre-selections for the MC
        if (fillTree) {
          if (flag & kTriton) {
            if (track.pt() < cfgCutPtMinTree || track.pt() > cfgCutPtMaxTree || track.sign() > 0)
              continue;
          }
        }
        nuclei::candidates.emplace_back(NucleusCandidate{
          static_cast<int>(track.globalIndex()), static_cast<int>(track.collisionId()), (1 - 2 * iC) * mTrackParCov.getPt(), mTrackParCov.getEta(), mTrackParCov.getPhi(),
          correctedTpcInnerParam, beta, collision.posZ(), collision.numContrib(), dcaInfo[0], dcaInfo[1], track.tpcSignal(), track.itsChi2NCl(), track.tpcChi2NCl(), track.tofChi2(),
          nSigmaTPC, tofMasses, fillTree, fillDCAHist, correctPV, isSecondary, fromWeakDecay, flag, track.tpcNClsFindable(), static_cast<uint8_t>(track.tpcNClsCrossedRows()), track.itsClusterMap(),
          static_cast<uint8_t>(track.tpcNClsFound()), static_cast<uint8_t>(track.tpcNClsShared()), static_cast<uint8_t>(track.itsNCls()), static_cast<uint32_t>(track.itsClusterSizes())});
      }
    } // end loop over tracks

    nuclei::hGloTOFtracks[0]->Fill(nGloTracks[0], nTOFTracks[0]);
    nuclei::hGloTOFtracks[1]->Fill(nGloTracks[1], nTOFTracks[1]);
  }

  void processData(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, TrackCandidates const& tracks, aod::BCsWithTimestamps const&)
  {
    nuclei::candidates.clear();
    if (!eventSelectionWithHisto(collision)) {
      return;
    }

    fillDataInfo(collision, tracks); // Apply selection on track candidates and saving in structure nuclei

    for (size_t i1{0}; i1 < nuclei::candidates.size(); ++i1) {
      auto& c1 = nuclei::candidates[i1];
      if (c1.fillTree) {
        nucleiTable(c1.pt, c1.eta, c1.phi, c1.tpcInnerParam, c1.beta, c1.zVertex, c1.nContrib, c1.DCAxy, c1.DCAz, c1.TPCsignal, c1.ITSchi2, c1.TPCchi2, c1.TOFchi2, c1.flags, c1.TPCfindableCls, c1.TPCcrossedRows, c1.ITSclsMap, c1.TPCnCls, c1.TPCnClsShared, c1.clusterSizesITS);
      }
      if (c1.fillDCAHist) {
        for (int iS{0}; iS < nuclei::species; ++iS) {
          if (c1.flags & BIT(iS)) {
            nuclei::hDCAHists[c1.pt < 0][iS]->Fill(std::abs(c1.pt), c1.DCAxy, c1.DCAz, c1.nSigmaTPC[iS], c1.tofMasses[iS], c1.ITSnCls, c1.TPCnCls);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(DeuteronInTriggeredEvents, processData, "Data analysis", true);

  Preslice<TrackCandidates> tracksPerCollisions = aod::track::collisionId;
  Preslice<o2::aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  void processMC(soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels> const& collisions, aod::McCollisions const& mcCollisions, soa::Join<TrackCandidates, aod::McTrackLabels> const& tracks, aod::McParticles const& particlesMC, aod::BCsWithTimestamps const&)
  {
    nuclei::candidates.clear();

    std::vector<bool> goodCollisions(mcCollisions.size(), false);
    std::vector<uint8_t> eventMask(mcCollisions.size(), 0);

    // Jet trigger condition
    auto trigger = static_cast<nuclei::triggerListName>(cfgTriggerList.value);

    for (const auto& c : mcCollisions) {

      spectra.fill(HIST("hGenVtxZ"), c.posZ());

      auto mcParticlesPerColl = particlesMC.sliceBy(perMcCollision, c.globalIndex());
      auto& mask = eventMask[c.globalIndex()];

      if (o2::pwglf::isINELgt0mc(mcParticlesPerColl, pdgDB))
        mask |= nuclei::kGenIsINELgt0;

      if (isMCJetTriggered(mcParticlesPerColl, particlesMC, trigger))
        mask |= nuclei::kGenIsJetTriggered;
    }

    for (const auto& collision : collisions) {
      if (!eventSelectionWithHisto(collision)) {
        continue;
      }

      // Avoid unwanted memory leaks
      if (!collision.has_mcCollision())
        continue;

      int mcId = collision.mcCollisionId();
      if (mcId < 0 || mcId >= static_cast<int>(eventMask.size()))
        continue;

      const auto& slicedTracks = tracks.sliceBy(tracksPerCollisions, collision.globalIndex());

      if (cfgApplyMCEvSel) {
        if (!isJetTriggered(slicedTracks, trigger))
          continue;
      }

      goodCollisions[mcId] = true;
      auto& mask = eventMask[mcId];
      mask |= nuclei::kGenHasRecoEv;

      fillDataInfo(collision, slicedTracks);
    }

    for (const auto& c : mcCollisions) {
      GenEventMCSel(eventMask[c.globalIndex()]);
    }

    std::vector<bool> isReconstructed(particlesMC.size(), false);
    for (auto& c : nuclei::candidates) {
      auto label = tracks.iteratorAt(c.globalIndex);
      if (label.mcParticleId() < -1 || label.mcParticleId() >= particlesMC.size()) {
        continue;
      }
      auto particle = particlesMC.iteratorAt(label.mcParticleId());
      bool storeIt{false};
      for (int iS{0}; iS < nuclei::species; ++iS) {
        if (std::abs(particle.pdgCode()) == nuclei::codes[iS]) {
          if (c.fillTree && !storeIt) {
            nuclei::hMomRes[iS][particle.pdgCode() < 0]->Fill(1., std::abs(c.pt * nuclei::charges[iS]), 1. - std::abs(c.pt * nuclei::charges[iS]) / particle.pt());
            storeIt = cfgTreeConfig->get(iS, 0u) || cfgTreeConfig->get(iS, 1u); /// store only the particles of interest
          }
          auto coll = collisions.iteratorAt(c.collTrackIndex);
          int collMCGlobId = coll.mcCollisionId();
          if (particle.mcCollisionId() == collMCGlobId) {
            c.correctPV = true;
          }
          if (!c.correctPV) {
            c.flags |= kIsAmbiguous;
          }
          if (c.fillDCAHist && cfgDCAHists->get(iS, c.pt < 0)) {
            nuclei::hDCAHists[c.pt < 0][iS]->Fill(std::abs(c.pt), c.DCAxy, c.DCAz, c.nSigmaTPC[iS], c.tofMasses[iS], c.ITSnCls, c.TPCnCls, c.correctPV, c.isSecondary, c.fromWeakDecay);
          }
        }
      }
      if (!storeIt) {
        continue;
      }
      if (particle.y() < cfgTrackCut.rapidityMin || particle.y() > cfgTrackCut.rapidityMax) {
        continue;
      }

      int motherPdgCode = 0;
      float motherDecRadius = -1;
      isReconstructed[particle.globalIndex()] = true;
      if (particle.isPhysicalPrimary()) {
        c.flags |= kIsPhysicalPrimary;
        if (particle.has_mothers()) {
          for (const auto& motherparticle : particle.mothers_as<aod::McParticles>()) {
            if (std::find(nuclei::hfMothCodes.begin(), nuclei::hfMothCodes.end(), std::abs(motherparticle.pdgCode())) != nuclei::hfMothCodes.end()) {
              c.flags |= kIsSecondaryFromWeakDecay;
              motherPdgCode = motherparticle.pdgCode();
              motherDecRadius = std::hypot(particle.vx() - motherparticle.vx(), particle.vy() - motherparticle.vy());
              break;
            }
          }
        }
      } else if (particle.has_mothers()) {
        c.flags |= kIsSecondaryFromWeakDecay;
        for (const auto& motherparticle : particle.mothers_as<aod::McParticles>()) {
          motherPdgCode = motherparticle.pdgCode();
          motherDecRadius = std::hypot(particle.vx() - motherparticle.vx(), particle.vy() - motherparticle.vy());
        }
      } else {
        c.flags |= kIsSecondaryFromMaterial;
      }

      isReconstructed[particle.globalIndex()] = true;
      float absoDecL = computeAbsoDecL(particle);
      nucleiTableMCExtension(c.pt, c.eta, c.phi, c.tpcInnerParam, c.beta, c.zVertex, c.nContrib, c.DCAxy, c.DCAz, c.TPCsignal, c.ITSchi2, c.TPCchi2, c.TOFchi2, c.flags, c.TPCfindableCls, c.TPCcrossedRows, c.ITSclsMap, c.TPCnCls, c.TPCnClsShared, c.clusterSizesITS, goodCollisions[particle.mcCollisionId()], particle.pt(), particle.eta(), particle.phi(), particle.pdgCode(), motherPdgCode, motherDecRadius, absoDecL, eventMask[particle.mcCollisionId()]);
    }

    int index{0};
    for (const auto& particle : particlesMC) {
      int pdg{std::abs(particle.pdgCode())};
      for (int iS{0}; iS < nuclei::species; ++iS) {
        if (pdg != nuclei::codes[iS]) {
          continue;
        }
        if (particle.y() < cfgTrackCut.rapidityMin || particle.y() > cfgTrackCut.rapidityMax) {
          continue;
        }

        uint16_t flags = 0;
        int motherPdgCode = 0;
        float motherDecRadius = -1;
        if (particle.isPhysicalPrimary()) {
          flags |= kIsPhysicalPrimary;
          nuclei::hGenNuclei[iS][particle.pdgCode() < 0]->Fill(1., particle.pt());
          // antinuclei from B hadrons are classified as physical primaries
          if (particle.has_mothers()) {
            for (const auto& motherparticle : particle.mothers_as<aod::McParticles>()) {
              if (std::find(nuclei::hfMothCodes.begin(), nuclei::hfMothCodes.end(), std::abs(motherparticle.pdgCode())) != nuclei::hfMothCodes.end()) {
                flags |= kIsSecondaryFromWeakDecay;
                motherPdgCode = motherparticle.pdgCode();
                motherDecRadius = std::hypot(particle.vx() - motherparticle.vx(), particle.vy() - motherparticle.vy());
                break;
              }
            }
          }
        } else if (particle.getProcess() == TMCProcess::kPDecay) {
          if (!particle.has_mothers()) {
            continue; // skip secondaries from weak decay without mothers
          }
          flags |= kIsSecondaryFromWeakDecay;
          for (const auto& motherparticle : particle.mothers_as<aod::McParticles>()) {
            motherPdgCode = motherparticle.pdgCode();
            motherDecRadius = std::hypot(particle.vx() - motherparticle.vx(), particle.vy() - motherparticle.vy());
          }
        } else {
          flags |= kIsSecondaryFromMaterial;
        }

        if (!isReconstructed[index] && (cfgTreeConfig->get(iS, 0u) || cfgTreeConfig->get(iS, 1u))) {
          if ((flags & kIsPhysicalPrimary) == 0 && cfgFillGenSecondaries == 0) {
            continue; // skip secondaries if not requested
          }
          if ((flags & (kIsPhysicalPrimary | kIsSecondaryFromWeakDecay)) == 0 && cfgFillGenSecondaries == 1) {
            continue; // skip secondaries from material if not requested
          }
          float absDecL = computeAbsoDecL(particle);
          nucleiTableMCExtension(999., 999., 999., 0., 0., 999., -1, 999., 999., -1, -1, -1, -1, flags, 0, 0, 0, 0, 0, 0, goodCollisions[particle.mcCollisionId()], particle.pt(), particle.eta(), particle.phi(), particle.pdgCode(), motherPdgCode, motherDecRadius, absDecL, eventMask[particle.mcCollisionId()]);
        }
        break;
      }
      index++;
    }
  }
  PROCESS_SWITCH(DeuteronInTriggeredEvents, processMC, "MC analysis", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<DeuteronInTriggeredEvents>(cfgc)};
}
