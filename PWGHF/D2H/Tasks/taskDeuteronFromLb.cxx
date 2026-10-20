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
//
/// \file taskDeuteronFromLb.cxx
/// \brief A filter task for non prompt deuterons from beauty-hadron decays
/// \author Marta Razza <marta.razza@cern.ch>
/// \author Francesca Ercolessi <francesca.ercolessi@cern.ch>
/// \since May 25, 2026

#include "Common/Core/Zorro.h"
#include "Common/Core/ZorroSummary.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>

#include <array>
#include <chrono>
#include <cmath>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct HfTaskDeuteronFromLb {

  Zorro zorro;
  o2::base::Propagator::MatCorrType noMatCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> applySkimming{"applySkimming", false, "Skimmed dataset processing"};
  Configurable<std::string> cfgSkimming{"cfgSkimming", "fH2fromLb", "Configurable for skimming"};
  Configurable<bool> sel8{"sel8", true, "sel8 event selection"};
  Configurable<float> cfgEta{"cfgEta", 0.8f, "Track eta selection"};
  Configurable<float> cfgTPCNclsFound{"cfgTPCNclsFound", 100, "Minimum TPC clusters found"};
  Configurable<float> cfgTPCChi2Ncl{"cfgTPCChi2Ncl", 4.0f, "Maximum TPC chi2 per N clusters"};
  Configurable<float> cfgITSChi2Ncl{"cfgITSChi2Ncl", 36.0f, "Maximum ITS chi2 per N clusters"};
  Configurable<float> cfgITScls{"cfgITScls", 2, "Minimum ITS clusters"};
  Configurable<float> cfgMaxPt{"cfgMaxPt", 5.0f, "Maximum pT cut"};
  Configurable<float> cfgMinPt{"cfgMinPt", 0.5f, "Minimum pT cut"};
  Configurable<float> cfgTPCNsigma{"cfgTPCNsigma", 3.0f, "TPC n sigma for deuteron PID"};
  Configurable<float> cfgTofNsigma{"cfgTofNsigma", 3.0f, "TOF n sigma for deuteron PID"};
  Configurable<float> cfgITSNsigma{"cfgITSNsigma", -2.0f, "ITS n sigma for deuteron PID"};
  Configurable<float> cfgDCAmin{"cfgDCAmin", 0.05f, "Minimum DCA for deuteron PID"};
  Configurable<float> cfgDCAmax{"cfgDCAmax", 1000.0f, "Maximum DCA for deuteron PID"};
  Configurable<float> rapidityCut{"rapidityCut", 0.5f, "Rapidity cut"};
  Configurable<float> ptThresholdforPID{"ptThresholdforPID", 1.0f, "pT threshold for PID, above this value TPC+TOF PID is used, below TPC+ITS PID is used"};
  // PDG codes
  Configurable<int> pdgCodeBeautyMeson{"pdgCodeBeautyMeson", -521, "PDG code of the beauty meson mother particle (default: B-)"};
  Configurable<int> pdgCodeBeautyBaryon{"pdgCodeBeautyBaryon", -5122, "PDG code of the beauty baryon mother particle (default: anti-Lambda_b)"};
  Configurable<int> pdgCodeDaughter{"pdgCodeDaughter", -1000010020, "PDG code of the daughter particle (default: anti-deuteron)"};
  Configurable<int> pdgCodeHypertriton{"pdgCodeHypertriton", -1010010030, "PDG code of the hypertriton mother particle (default: anti-hypertriton)"};

  int mRunNumber = 0;
  int mCurrentRun = -1;

  framework::Service<ccdb::BasicCCDBManager> ccdb{};

  using CollisionCandidates = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
  using TrackCandidates = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::TracksExtra, o2::aod::TracksDCA, o2::aod::TrackSelection, o2::aod::pidTPCFullDe, o2::aod::pidTOFFullDe, o2::aod::pidTPCFullEl, o2::aod::pidTPCFullPi, o2::aod::pidTPCFullKa, o2::aod::pidTPCFullPr, o2::aod::pidTOFFullEl, o2::aod::pidTOFFullPi, o2::aod::pidTOFFullKa, o2::aod::pidTOFFullPr>;
  using MCTrackCandidates = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksCov, o2::aod::TracksExtra, o2::aod::TracksDCA, o2::aod::McTrackLabels, o2::aod::pidTPCFullDe, o2::aod::pidTOFFullDe, o2::aod::pidTPCFullEl, o2::aod::pidTPCFullPi, o2::aod::pidTPCFullKa, o2::aod::pidTPCFullPr, o2::aod::pidTOFFullEl, o2::aod::pidTOFFullPi, o2::aod::pidTOFFullKa, o2::aod::pidTOFFullPr>;
  using MCCollisionCandidates = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels, o2::aod::McCollisionLabels>;

  Preslice<o2::aod::TrackAssoc> trackIndicesPerCollision = o2::aod::track_association::collisionId;
  Preslice<o2::aod::McParticles> mcParticlesPerCollision = o2::aod::mcparticle::mcCollisionId;

  ConfigurableAxis ptAxis{"ptAxis", {100, 0.f, 10.f}, "p_{T} GeV/c"};
  ConfigurableAxis nSigmaAxis{"nSigmaAxis", {200, -10.f, 10.f}, "nSigma"};
  ConfigurableAxis dcaXyAxis{"dcaXyAxis", {4000, -2.0f, 2.0f}, "DCA xy (cm)"};
  ConfigurableAxis dcaZAxis{"dcaZAxis", {4000, -2.0f, 2.0f}, "DCA z (cm)"};
  ConfigurableAxis ctauHypertritonAxis{"ctauHypertritonAxis", {100, 0.f, 30.f}, "c#tau (cm)"};

  HistogramRegistry qaHistos{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};
  OutputObj<TH1D> hProcessedEvents{TH1D("hProcessedEvents", "Event filtered;; Number of events", 4, 0., 4.)};

  bool isSelectedMother(int pdgCode)
  {
    return pdgCode == pdgCodeBeautyMeson || pdgCode == pdgCodeBeautyBaryon || pdgCode == pdgCodeHypertriton;
  }

  double getMotherMass(int pdgCode)
  {
    if (pdgCode == pdgCodeBeautyMeson) {
      return o2::constants::physics::MassBPlus;
    }

    if (pdgCode == pdgCodeBeautyBaryon) {
      return o2::constants::physics::MassLambdaB0;
    }
    if (pdgCode == pdgCodeHypertriton) {
      return o2::constants::physics::MassHyperTriton;
    }
    return -1.;
  }

  void init(framework::InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());

    if (applySkimming) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }

    qaHistos.add("Data/hDCAxyVsPt", "DCAxy #bar{d} vs p_{T}", {HistType::kTH2D, {ptAxis, dcaXyAxis}});
    qaHistos.add("Data/hDCAzVsPt", "DCAz #bar{d} vs p_{T}", {HistType::kTH2D, {ptAxis, dcaZAxis}});
    qaHistos.add("Data/hnSigmaTPCVsPtFinalSel", "n#sigma TPC vs p_{T} for #bar{d} hypothesis; p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2D, {ptAxis, nSigmaAxis}});
    qaHistos.add("Data/hnSigmaTOFVsPtFinalSel", "n#sigma TOF vs p_{T} for #bar{d} hypothesis; p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2D, {ptAxis, nSigmaAxis}});
    qaHistos.add("Data/ptAntiDeuteron", "ptAntiDeuteron", {HistType::kTH1F, {ptAxis}});
    qaHistos.add("Data/etaAntideuteron", "etaAntideuteron", {HistType::kTH1F, {{100, -1.0f, 1.0f, "eta #bar{d}"}}});
    qaHistos.add("Data/hVtxZ", "Z-Vertex distribution after selection;Z (cm)", HistType::kTH1F, {{100, -50, 50}});
    qaHistos.add("Data/hnSigmaTOFVsPtPurityNoVeto", "n#sigma TOF vs p_{T} for #bar{d} hypothesis for Data-driven purity check; p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2D, {ptAxis, nSigmaAxis}});
    qaHistos.add("Data/hnSigmaTOFVsPtPurityVeto", "n#sigma TOF vs p_{T} for #bar{d} hypothesis for Data-driven purity check with TPC+TOF veto for e,pi,K,p; p_{T} (GeV/c); n#sigma TOF", {HistType::kTH2D, {ptAxis, nSigmaAxis}});
    qaHistos.add("Data/hnSigmaTPCVsPtPurity", "n#sigma TPC vs p_{T} for #bar{d} hypothesis for Data-driven purity check; p_{T} (GeV/c); n#sigma TPC", {HistType::kTH2D, {ptAxis, nSigmaAxis}});
    qaHistos.add("Data/hnSigmaITSVsPt", "n#sigma ITS vs p_{T} for #bar{d} hypothesis, used to clean up the TPC/TOF purity sample; p_{T} (GeV/c); n#sigma ITS", {HistType::kTH2D, {ptAxis, nSigmaAxis}});
    // MC generated-level histograms
    qaHistos.add("MCGen/ptGeneratedBminus", "p_{T} generated B^{-};p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCGen/ptGeneratedAntiLambdaB", "p_{T} generated #bar{#Lambda}_{b};p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCGen/ptAntiDeuteronPrimary", "p_{T} #bar{d} primary gen;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCGen/hMotherPdgCode", "PDG code of mother, gen level;PDG code;Counts", HistType::kTH1I, {{12000, -6000, 6000}});
    qaHistos.add("MCGen/ctauBminus", "ctau of B^{-}, gen level;ctau (#mu m);Counts", HistType::kTH1F, {{25, 0., 2000.f}});
    qaHistos.add("MCGen/ctauAntiLambdaB", "ctau of #bar{#Lambda}_{b}, gen level;ctau (#mu m);Counts", HistType::kTH1F, {{25, 0., 2000.f}});
    qaHistos.add("MCGen/ptAntiDeuteronFromBminus", "p_{T} #bar{d} from B^{-} gen;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCGen/ptAntiDeuteronFromAntiLambdaB", "p_{T} #bar{d} from #bar{#Lambda}_{b} gen;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCGen/ptAntiDeuteronAll", "p_{T} #bar{d} all gen;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCGen/hEvents", "Events MC gen;;Counts", HistType::kTH1F, {{2, 0., 2.}});
    qaHistos.add("MCGen/ptGeneratedAntiHypertriton", "p_{T} generated ^{3}_{#bar{#Lambda}}#bar{H};p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCGen/ctauAntiHypertriton", "ctau of ^{3}_{#bar{#Lambda}}#bar{H}, gen level;c#tau (cm);Counts", HistType::kTH1F, {ctauHypertritonAxis});
    qaHistos.add("MCGen/ptAntiDeuteronFromAntiHypertriton", "p_{T} #bar{d} from ^{3}_{#bar{#Lambda}}#bar{H} gen;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    // MC reco/MC-anchored histograms
    qaHistos.add("MCReco/hDCAxyVsPtDeFromBminus", "DCAxy #bar{d} from B^{-} vs p_{T};p_{T} (GeV/c);DCA_{xy} (cm)", {HistType::kTH2D, {ptAxis, dcaXyAxis}});
    qaHistos.add("MCReco/hDCAxyVsPtDeFromAntiLambdaB", "DCAxy #bar{d} from #bar{#Lambda}_{b} vs p_{T};p_{T} (GeV/c);DCA_{xy} (cm)", {HistType::kTH2D, {ptAxis, dcaXyAxis}});
    qaHistos.add("MCReco/hDCAxyVsPtDePrimary", "DCAxy #bar{d} primary vs p_{T};p_{T} (GeV/c);DCA_{xy} (cm)", {HistType::kTH2D, {ptAxis, dcaXyAxis}});
    qaHistos.add("MCReco/hDCAzVsPtDeFromBminus", "DCAz #bar{d} from B^{-} vs p_{T};p_{T} (GeV/c);DCA_{z} (cm)", {HistType::kTH2D, {ptAxis, dcaZAxis}});
    qaHistos.add("MCReco/hDCAzVsPtDeFromAntiLambdaB", "DCAz #bar{d} from #bar{#Lambda}_{b} vs p_{T};p_{T} (GeV/c);DCA_{z} (cm)", {HistType::kTH2D, {ptAxis, dcaZAxis}});
    qaHistos.add("MCReco/hDCAzVsPtDePrimary", "DCAz #bar{d} primary vs p_{T};p_{T} (GeV/c);DCA_{z} (cm)", {HistType::kTH2D, {ptAxis, dcaZAxis}});
    qaHistos.add("MCReco/hDCAxyVsPtDeFromBminusSameMC", "DCAxy from B^{-}, same MC collision", {HistType::kTH2D, {ptAxis, dcaXyAxis}});
    qaHistos.add("MCReco/hDCAxyVsPtDeFromAntiLambdaBSameMC", "DCAxy from anti-Lambda_b, same MC collision", {HistType::kTH2D, {ptAxis, dcaXyAxis}});
    qaHistos.add("MCReco/hDCAxyVsPtDeFromAntiHypertritonSameMC", "DCAxy from anti-hypertriton, same MC collision", {HistType::kTH2D, {ptAxis, dcaXyAxis}});
    qaHistos.add("MCReco/hDCAxyVsPtDePrimarySameMC", "DCAxy primary, same MC collision", {HistType::kTH2D, {ptAxis, dcaXyAxis}});
    qaHistos.add("MCReco/ptAntiDeuteronFromBminus", "p_{T} #bar{d} from B^{-} reco/MC anchored;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCReco/ptAntiDeuteronFromAntiLambdaB", "p_{T} #bar{d} from #bar{#Lambda}_{b} reco/MC anchored;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCReco/ptAntiDeuteronPrimary", "p_{T} #bar{d} primary reco/MC anchored;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCReco/hEvents", "Events MC reco;;Counts", HistType::kTH1F, {{3, 0., 3.}});
    qaHistos.add("MCReco/hDCAxy-Primary", "DCAxy primary reco/MC anchored;DCA xy (cm);Counts", {HistType::kTH1D, {dcaXyAxis}});
    qaHistos.add("MCReco/hDCAxy-FromBMeson", "DCAxy from B^{-} reco/MC anchored;DCA xy (cm);Counts", {HistType::kTH1D, {dcaXyAxis}});
    qaHistos.add("MCReco/hDCAxy-FromAntiLambdaB", "DCAxy from #bar{#Lambda}_{b} reco/MC anchored;DCA xy (cm);Counts", {HistType::kTH1D, {dcaXyAxis}});
    qaHistos.add("MCReco/hMotherPdgCode", "PDG code of mother, reco/MC anchored;PDG code;Counts", HistType::kTH1I, {{12000, -6000, 6000}});
    qaHistos.add("MCReco/ctauBminus", "ctau of B^{-}, reco/MC anchored;ctau (#mu m);Counts", HistType::kTH1F, {{25, 0., 2000.f}});
    qaHistos.add("MCReco/ctauAntiLambdaB", "ctau of #bar{#Lambda}_{b}, reco/MC anchored;ctau (#mu m);Counts", HistType::kTH1F, {{25, 0., 2000.f}});
    qaHistos.add("MCReco/ptAntiDeuteronPIDSelected", "p_{T} #bar{d} reco/MC anchored PID selection;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCReco/ptAntiDeuteronPIDSelectedTrue", "p_{T} #bar{d} reco/MC anchored PID and PDG;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCReco/ptAntiDeuteronTPCOnlyPIDSelected", "p_{T} #bar{d} reco/MC anchored TPC only PID selection;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCReco/ptAntiDeuteronTPCOnlyPIDSelectedTrue", "p_{T} #bar{d} reco/MC anchored TPC only PID and PDG;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCReco/ptAntiDeuteronTPCorTPCTOFPIDSelected", "p_{T} #bar{d} reco/MC anchored TPC or TPC+TOF PID selection;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCReco/ptAntiDeuteronTPCorTPCTOFPIDSelectedTrue", "p_{T} #bar{d} reco/MC anchored TPC or TPC+TOF PID and PDG;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCReco/ptAntiDeuteronGenFromBminus", "p_{T} #bar{d} from B^{-} gen in reco events;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCReco/ptAntiDeuteronGenFromAntiLambdaB", "p_{T} #bar{d} from #bar{#Lambda}_{b} gen in reco events;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCReco/hDCAxyVsPtDeFromAntiHypertriton", "DCAxy #bar{d} from ^{3}_{#bar{#Lambda}}#bar{H} vs p_{T};p_{T} (GeV/c);DCA_{xy} (cm)", {HistType::kTH2D, {ptAxis, dcaXyAxis}});
    qaHistos.add("MCReco/hDCAzVsPtDeFromAntiHypertriton", "DCAz #bar{d} from ^{3}_{#bar{#Lambda}}#bar{H} vs p_{T};p_{T} (GeV/c);DCA_{z} (cm)", {HistType::kTH2D, {ptAxis, dcaZAxis}});
    qaHistos.add("MCReco/ptAntiDeuteronFromAntiHypertriton", "p_{T} #bar{d} from ^{3}_{#bar{#Lambda}}#bar{H} reco/MC anchored;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCReco/hDCAxy-FromAntiHypertriton", "DCAxy from ^{3}_{#bar{#Lambda}}#bar{H} reco/MC anchored;DCA xy (cm);Counts", {HistType::kTH1D, {dcaXyAxis}});
    qaHistos.add("MCReco/ctauAntiHypertriton", "ctau of ^{3}_{#bar{#Lambda}}#bar{H}, reco/MC anchored;c#tau (cm);Counts", HistType::kTH1F, {ctauHypertritonAxis});
    qaHistos.add("MCReco/ptAntiDeuteronGenFromAntiHypertriton", "p_{T} #bar{d} from ^{3}_{#bar{#Lambda}}#bar{H} gen in reco events;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});
    qaHistos.add("MCReco/ptAntiDeuteronGenPrimary", "p_{T} #bar{d} primary gen in reco events;p_{T} (GeV/c);Counts", HistType::kTH1F, {ptAxis});

    hProcessedEvents->GetXaxis()->SetBinLabel(1, "Events processed");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "ZORRO");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, "sel8");
    hProcessedEvents->GetXaxis()->SetBinLabel(4, "z-vertex");
    auto hEvents = qaHistos.get<TH1>(HIST("MCReco/hEvents"));
    hEvents->GetXaxis()->SetBinLabel(1, "All");
    hEvents->GetXaxis()->SetBinLabel(2, "sel8");
    hEvents->GetXaxis()->SetBinLabel(3, "z-vertex");
    auto hEventsGen = qaHistos.get<TH1>(HIST("MCGen/hEvents"));
    hEventsGen->GetXaxis()->SetBinLabel(1, "All");
    hEventsGen->GetXaxis()->SetBinLabel(2, "z-vertex");

    qaHistos.add("MCReco/hMotherCategory", "Mother of #bar{d};;Counts", HistType::kTH1F, {{5, 0., 5.}});
    auto hMotherCategory = qaHistos.get<TH1>(HIST("MCReco/hMotherCategory"));
    hMotherCategory->GetXaxis()->SetBinLabel(1, "B^{-}");
    hMotherCategory->GetXaxis()->SetBinLabel(2, "#bar{#Lambda}_{b}");
    hMotherCategory->GetXaxis()->SetBinLabel(3, "^{3}_{#bar{#Lambda}}#bar{H}");
    hMotherCategory->GetXaxis()->SetBinLabel(4, "primary");
    hMotherCategory->GetXaxis()->SetBinLabel(5, "other");
  }

  void initCCDB(o2::aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    if (applySkimming) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), cfgSkimming.value);
      zorro.populateHistRegistry(registry, bc.runNumber());
    }
    mRunNumber = bc.runNumber();
  }

  template <typename T1>
  bool passedSingleTrackSelection(const T1& track)
  {
    if (std::abs(track.eta()) > cfgEta) {
      return false;
    }
    if (std::abs(track.rapidity(o2::constants::physics::MassDeuteron)) > rapidityCut) {
      return false;
    }
    if (!track.hasITS()) {
      return false;
    }
    if (!track.hasTPC()) {
      return false;
    }
    if (track.tpcNClsFound() < cfgTPCNclsFound) {
      return false;
    }
    if (track.tpcChi2NCl() > cfgTPCChi2Ncl) {
      return false;
    }
    if (track.itsChi2NCl() > cfgITSChi2Ncl) {
      return false;
    }
    if (track.itsNCls() < cfgITScls) {
      return false;
    }
    if (track.pt() > cfgMaxPt) {
      return false;
    }
    if (track.pt() < cfgMinPt) {
      return false;
    }
    if (track.sign() > 0) {
      return false;
    }

    return true;
  }

  void processData(CollisionCandidates const& collisions,
                   o2::aod::TrackAssoc const& trackIndices,
                   TrackCandidates const& tracks,
                   o2::aod::BCsWithTimestamps const&)
  {
    auto tracksWithItsPid = o2::soa::Attach<TrackCandidates, o2::aod::pidits::ITSNSigmaDe>(tracks);
    for (const auto& collision : collisions) {
      if (mCurrentRun != collision.bc_as<o2::aod::BCsWithTimestamps>().runNumber()) {
        auto* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", collision.bc_as<o2::aod::BCsWithTimestamps>().timestamp());
        o2::base::Propagator::initFieldFromGRP(grpo);
        mCurrentRun = collision.bc_as<o2::aod::BCsWithTimestamps>().runNumber();
      }

      const auto& bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
      initCCDB(bc);
      hProcessedEvents->Fill(0.5);
      if (applySkimming) {
        if (!zorro.isSelected(bc.globalBC())) {
          continue;
        }
      }
      hProcessedEvents->Fill(1.5);
      if (sel8 && !collision.sel8()) {
        continue;
      }
      hProcessedEvents->Fill(2.5);
      if (std::abs(collision.posZ()) > cutzvertex) {
        continue;
      }
      hProcessedEvents->Fill(3.5);
      qaHistos.fill(HIST("Data/hVtxZ"), collision.posZ());

      const auto& trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());

      for (const auto& trackId : trackIdsThisCollision) {
        const auto& track = tracksWithItsPid.rawIteratorAt(trackId.trackId());
        std::array<float, 2> dca{track.dcaXY(), track.dcaZ()};

        if (track.collisionId() != collision.globalIndex()) {
          auto trackPar = getTrackParCov(track);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, noMatCorr, &dca);
        }

        if (std::abs(dca[0]) < cfgDCAmin || std::abs(dca[0]) > cfgDCAmax) {
          continue;
        }

        if (!passedSingleTrackSelection(track)) {
          continue;
        }

        const bool isTPCDe = std::abs(track.tpcNSigmaDe()) < cfgTPCNsigma;
        const bool isTOFDe = std::abs(track.tofNSigmaDe()) < cfgTofNsigma;
        const bool isITSDe = track.itsNSigmaDe() > cfgITSNsigma;

        qaHistos.fill(HIST("Data/hnSigmaITSVsPt"), track.pt(), track.itsNSigmaDe());

        const bool isElectron = std::abs(track.tpcNSigmaEl()) < cfgTPCNsigma && std::abs(track.tofNSigmaEl()) < cfgTofNsigma;
        const bool isNotElectron = !isElectron;
        const bool isPion = std::abs(track.tpcNSigmaPi()) < cfgTPCNsigma && std::abs(track.tofNSigmaPi()) < cfgTofNsigma;
        const bool isNotPion = !isPion;
        const bool isKaon = std::abs(track.tpcNSigmaKa()) < cfgTPCNsigma && std::abs(track.tofNSigmaKa()) < cfgTofNsigma;
        const bool isNotKaon = !isKaon;
        const bool isProton = std::abs(track.tpcNSigmaPr()) < cfgTPCNsigma && std::abs(track.tofNSigmaPr()) < cfgTofNsigma;
        const bool isNotProton = !isProton;

        const bool isTPCTOFDeClean = isTPCDe && isNotElectron && isNotPion && isNotKaon && isNotProton;

        if (track.pt() < ptThresholdforPID) {
          if (isITSDe) {
            qaHistos.fill(HIST("Data/hnSigmaTPCVsPtPurity"), track.pt(), track.tpcNSigmaDe());
          }
        } else {
          if (isTPCDe && track.hasTOF()) {
            qaHistos.fill(HIST("Data/hnSigmaTOFVsPtPurityNoVeto"), track.pt(), track.tofNSigmaDe());
          }
          if (isTPCTOFDeClean && track.hasTOF()) {
            qaHistos.fill(HIST("Data/hnSigmaTOFVsPtPurityVeto"), track.pt(), track.tofNSigmaDe());
          }
        }

        const bool passFinalSel = (track.pt() < ptThresholdforPID && isTPCDe && isITSDe) ||
                                  (track.pt() >= ptThresholdforPID && isTPCTOFDeClean && isTOFDe && track.hasTOF());

        if (passFinalSel) {
          qaHistos.fill(HIST("Data/ptAntiDeuteron"), track.pt());
          qaHistos.fill(HIST("Data/etaAntideuteron"), track.eta());
          qaHistos.fill(HIST("Data/hDCAxyVsPt"), track.pt(), dca[0]);
          qaHistos.fill(HIST("Data/hDCAzVsPt"), track.pt(), dca[1]);
          qaHistos.fill(HIST("Data/hnSigmaTPCVsPtFinalSel"), track.pt(), track.tpcNSigmaDe());
          if (track.pt() >= ptThresholdforPID) {
            qaHistos.fill(HIST("Data/hnSigmaTOFVsPtFinalSel"), track.pt(), track.tofNSigmaDe());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskDeuteronFromLb, processData, "processData", true);

  void processMC(MCCollisionCandidates const& collisions,
                 o2::aod::TrackAssoc const& trackIndices,
                 MCTrackCandidates const& tracks,
                 o2::aod::McParticles const& mcParticles,
                 o2::aod::BCsWithTimestamps const&)
  {
    auto tracksWithItsPid = o2::soa::Attach<MCTrackCandidates, o2::aod::pidits::ITSNSigmaDe>(tracks);
    tracksWithItsPid.bindExternalIndices(&mcParticles);

    for (const auto& collision : collisions) {
      const auto& bc = collision.bc_as<o2::aod::BCsWithTimestamps>();
      if (mCurrentRun != bc.runNumber()) {
        auto* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", bc.timestamp());
        o2::base::Propagator::initFieldFromGRP(grpo);
        mCurrentRun = bc.runNumber();
      }

      qaHistos.fill(HIST("MCReco/hEvents"), 0.5);
      if (sel8 && !collision.sel8()) {
        continue;
      }
      qaHistos.fill(HIST("MCReco/hEvents"), 1.5);
      if (std::abs(collision.posZ()) > cutzvertex) {
        continue;
      }
      qaHistos.fill(HIST("MCReco/hEvents"), 2.5);

      if (collision.has_mcCollision()) {
        const auto mcParticlesThisEvent = mcParticles.sliceBy(mcParticlesPerCollision, collision.mcCollisionId());

        for (const auto& mcParticle : mcParticlesThisEvent) {

          if (mcParticle.pdgCode() != pdgCodeDaughter) {
            continue;
          }
          if (std::abs(mcParticle.eta()) > cfgEta || std::abs(mcParticle.y()) > rapidityCut) {
            continue;
          }
          if (mcParticle.pt() < cfgMinPt || mcParticle.pt() > cfgMaxPt) {
            continue;
          }
          bool fromMother = false;
          for (const auto& mom : mcParticle.mothers_as<o2::aod::McParticles>()) {
            if (mom.pdgCode() == pdgCodeBeautyMeson) {
              qaHistos.fill(HIST("MCReco/ptAntiDeuteronGenFromBminus"), mcParticle.pt());
              fromMother = true;
              break;
            }
            if (mom.pdgCode() == pdgCodeBeautyBaryon) {
              qaHistos.fill(HIST("MCReco/ptAntiDeuteronGenFromAntiLambdaB"), mcParticle.pt());
              fromMother = true;
              break;
            }
            if (mom.pdgCode() == pdgCodeHypertriton) {
              qaHistos.fill(HIST("MCReco/ptAntiDeuteronGenFromAntiHypertriton"), mcParticle.pt());
              fromMother = true;
              break;
            }
          }
          if (!fromMother && mcParticle.isPhysicalPrimary()) {
            qaHistos.fill(HIST("MCReco/ptAntiDeuteronGenPrimary"), mcParticle.pt());
          }
        }
      }

      const auto& trackIdsThisCollision = trackIndices.sliceBy(trackIndicesPerCollision, collision.globalIndex());

      for (const auto& trackId : trackIdsThisCollision) {
        const auto& track = tracksWithItsPid.rawIteratorAt(trackId.trackId());
        std::array<float, 2> dca{track.dcaXY(), track.dcaZ()};

        if (track.collisionId() != collision.globalIndex()) {
          auto trackPar = getTrackParCov(track);
          o2::base::Propagator::Instance()->propagateToDCABxByBz({collision.posX(), collision.posY(), collision.posZ()}, trackPar, 2.f, noMatCorr, &dca);
        }

        if (std::abs(dca[0]) < cfgDCAmin || std::abs(dca[0]) > cfgDCAmax) {
          continue;
        }
        if (!passedSingleTrackSelection(track)) {
          continue;
        }

        if (!track.has_mcParticle()) {
          continue;
        }

        auto mcParticle = track.mcParticle();
        const bool isSameMcCollision = collision.has_mcCollision() && (mcParticle.mcCollisionId() == collision.mcCollisionId());
        const bool isTrueDeuteron = (mcParticle.pdgCode() == pdgCodeDaughter);
        const bool isTPCDe = std::abs(track.tpcNSigmaDe()) < cfgTPCNsigma;
        const bool isTOFDe = std::abs(track.tofNSigmaDe()) < cfgTofNsigma;
        const bool isITSDe = track.itsNSigmaDe() > cfgITSNsigma;

        const bool selTpcOrTpcTof = isTPCDe && (!track.hasTOF() || isTOFDe);

        if (isTPCDe) {
          qaHistos.fill(HIST("MCReco/ptAntiDeuteronTPCOnlyPIDSelected"), track.pt());
          if (isTrueDeuteron) {
            qaHistos.fill(HIST("MCReco/ptAntiDeuteronTPCOnlyPIDSelectedTrue"), track.pt());
          }
        }

        if (selTpcOrTpcTof) {
          qaHistos.fill(HIST("MCReco/ptAntiDeuteronTPCorTPCTOFPIDSelected"), track.pt());
          if (isTrueDeuteron) {
            qaHistos.fill(HIST("MCReco/ptAntiDeuteronTPCorTPCTOFPIDSelectedTrue"), track.pt());
          }
        }
        // Purity study: TPC && TOF
        if (isTPCDe && isTOFDe && track.hasTOF()) {
          qaHistos.fill(HIST("MCReco/ptAntiDeuteronPIDSelected"), track.pt());
          if (isTrueDeuteron) {
            qaHistos.fill(HIST("MCReco/ptAntiDeuteronPIDSelectedTrue"), track.pt());
          }
        }

        const bool isElectron = std::abs(track.tpcNSigmaEl()) < cfgTPCNsigma && std::abs(track.tofNSigmaEl()) < cfgTofNsigma;
        const bool isNotElectron = !isElectron;
        const bool isPion = std::abs(track.tpcNSigmaPi()) < cfgTPCNsigma && std::abs(track.tofNSigmaPi()) < cfgTofNsigma;
        const bool isNotPion = !isPion;
        const bool isKaon = std::abs(track.tpcNSigmaKa()) < cfgTPCNsigma && std::abs(track.tofNSigmaKa()) < cfgTofNsigma;
        const bool isNotKaon = !isKaon;
        const bool isProton = std::abs(track.tpcNSigmaPr()) < cfgTPCNsigma && std::abs(track.tofNSigmaPr()) < cfgTofNsigma;
        const bool isNotProton = !isProton;

        const bool isTPCTOFDeClean = isTPCDe && isNotElectron && isNotPion && isNotKaon && isNotProton;

        if (isTrueDeuteron) {
          int motherPdg = 0;
          for (const auto& mom : mcParticle.mothers_as<o2::aod::McParticles>()) {
            qaHistos.fill(HIST("MCReco/hMotherPdgCode"), mom.pdgCode());
            if (mom.pdgCode() == pdgCodeBeautyMeson || mom.pdgCode() == pdgCodeBeautyBaryon || mom.pdgCode() == pdgCodeHypertriton) {
              motherPdg = mom.pdgCode();
              double dx = mcParticle.vx() - mom.vx();
              double dy = mcParticle.vy() - mom.vy();
              double dz = mcParticle.vz() - mom.vz();

              double flightDistance = std::sqrt(dx * dx + dy * dy + dz * dz);
              double massBeauty = getMotherMass(mom.pdgCode());

              if (massBeauty > 0.) {
                double ctauBeauty = flightDistance / mom.p() * massBeauty;

                if (motherPdg == pdgCodeBeautyMeson) {
                  qaHistos.fill(HIST("MCReco/ctauBminus"), ctauBeauty * 1.e4);
                } else if (motherPdg == pdgCodeBeautyBaryon) {
                  qaHistos.fill(HIST("MCReco/ctauAntiLambdaB"), ctauBeauty * 1.e4);
                } else if (motherPdg == pdgCodeHypertriton) {
                  qaHistos.fill(HIST("MCReco/ctauAntiHypertriton"), ctauBeauty); // cm
                }
              }
              break;
            }
          }

          if (motherPdg == pdgCodeBeautyMeson) {
            qaHistos.fill(HIST("MCReco/hMotherCategory"), 0.5);
          } else if (motherPdg == pdgCodeBeautyBaryon) {
            qaHistos.fill(HIST("MCReco/hMotherCategory"), 1.5);
          } else if (motherPdg == pdgCodeHypertriton) {
            qaHistos.fill(HIST("MCReco/hMotherCategory"), 2.5);
          } else if (mcParticle.isPhysicalPrimary()) {
            qaHistos.fill(HIST("MCReco/hMotherCategory"), 3.5);
          } else {
            qaHistos.fill(HIST("MCReco/hMotherCategory"), 4.5);
          }

          // same final selection as in data
          const bool passFinalSel = (track.pt() < ptThresholdforPID && isTPCDe && isITSDe) ||
                                    (track.pt() >= ptThresholdforPID && isTPCTOFDeClean && isTOFDe && track.hasTOF());
          if (passFinalSel) {
            if (motherPdg == pdgCodeBeautyMeson) {
              qaHistos.fill(HIST("MCReco/hDCAxyVsPtDeFromBminus"), track.pt(), dca[0]);
              if (isSameMcCollision) {
                qaHistos.fill(HIST("MCReco/hDCAxyVsPtDeFromBminusSameMC"), track.pt(), dca[0]);
              }
              qaHistos.fill(HIST("MCReco/hDCAzVsPtDeFromBminus"), track.pt(), dca[1]);
              qaHistos.fill(HIST("MCReco/hDCAxy-FromBMeson"), dca[0]);
              qaHistos.fill(HIST("MCReco/ptAntiDeuteronFromBminus"), track.pt());
            } else if (motherPdg == pdgCodeBeautyBaryon) {
              qaHistos.fill(HIST("MCReco/hDCAxyVsPtDeFromAntiLambdaB"), track.pt(), dca[0]);
              if (isSameMcCollision) {
                qaHistos.fill(HIST("MCReco/hDCAxyVsPtDeFromAntiLambdaBSameMC"), track.pt(), dca[0]);
              }
              qaHistos.fill(HIST("MCReco/hDCAzVsPtDeFromAntiLambdaB"), track.pt(), dca[1]);
              qaHistos.fill(HIST("MCReco/hDCAxy-FromAntiLambdaB"), dca[0]);
              qaHistos.fill(HIST("MCReco/ptAntiDeuteronFromAntiLambdaB"), track.pt());
            } else if (motherPdg == pdgCodeHypertriton) {
              qaHistos.fill(HIST("MCReco/hDCAxyVsPtDeFromAntiHypertriton"), track.pt(), dca[0]);
              if (isSameMcCollision) {
                qaHistos.fill(HIST("MCReco/hDCAxyVsPtDeFromAntiHypertritonSameMC"), track.pt(), dca[0]);
              }
              qaHistos.fill(HIST("MCReco/hDCAzVsPtDeFromAntiHypertriton"), track.pt(), dca[1]);
              qaHistos.fill(HIST("MCReco/hDCAxy-FromAntiHypertriton"), dca[0]);
              qaHistos.fill(HIST("MCReco/ptAntiDeuteronFromAntiHypertriton"), track.pt());
            } else if (mcParticle.isPhysicalPrimary()) {
              qaHistos.fill(HIST("MCReco/hDCAxyVsPtDePrimary"), track.pt(), dca[0]);
              qaHistos.fill(HIST("MCReco/hDCAzVsPtDePrimary"), track.pt(), dca[1]);
              qaHistos.fill(HIST("MCReco/hDCAxy-Primary"), dca[0]);
              if (isSameMcCollision) {
                qaHistos.fill(HIST("MCReco/hDCAxyVsPtDePrimarySameMC"), track.pt(), dca[0]);
              }
              qaHistos.fill(HIST("MCReco/ptAntiDeuteronPrimary"), track.pt());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(HfTaskDeuteronFromLb, processMC, "processMC", false);

  void processGen(o2::aod::McCollision const& collision, o2::aod::McParticles const& mcParticles)
  {
    qaHistos.fill(HIST("MCGen/hEvents"), 0.5);
    if (std::abs(collision.posZ()) > cutzvertex) {
      return;
    }
    qaHistos.fill(HIST("MCGen/hEvents"), 1.5);

    for (const auto& mcParticle : mcParticles) {
      if (isSelectedMother(mcParticle.pdgCode())) {
        if (mcParticle.pdgCode() == pdgCodeBeautyMeson) {
          qaHistos.fill(HIST("MCGen/ptGeneratedBminus"), mcParticle.pt());
        } else if (mcParticle.pdgCode() == pdgCodeBeautyBaryon) {
          qaHistos.fill(HIST("MCGen/ptGeneratedAntiLambdaB"), mcParticle.pt());
        } else if (mcParticle.pdgCode() == pdgCodeHypertriton) {
          qaHistos.fill(HIST("MCGen/ptGeneratedAntiHypertriton"), mcParticle.pt());
        }
        if (mcParticle.has_daughters()) {

          for (const auto& daughter : mcParticle.daughters_as<o2::aod::McParticles>()) {

            double dx = daughter.vx() - mcParticle.vx();
            double dy = daughter.vy() - mcParticle.vy();
            double dz = daughter.vz() - mcParticle.vz();
            double flightDistance = std::sqrt(dx * dx + dy * dy + dz * dz);
            double massBeauty = getMotherMass(mcParticle.pdgCode());
            if (massBeauty > 0.) {
              double ctauBeauty = flightDistance / mcParticle.p() * massBeauty;

              if (mcParticle.pdgCode() == pdgCodeBeautyMeson) {
                qaHistos.fill(HIST("MCGen/ctauBminus"), ctauBeauty * 1.e4);
              } else if (mcParticle.pdgCode() == pdgCodeBeautyBaryon) {
                qaHistos.fill(HIST("MCGen/ctauAntiLambdaB"), ctauBeauty * 1.e4);
              } else if (mcParticle.pdgCode() == pdgCodeHypertriton) {
                qaHistos.fill(HIST("MCGen/ctauAntiHypertriton"), ctauBeauty); // cm
              }
            }
            break;
          }
        }
      }

      if (mcParticle.pdgCode() == pdgCodeDaughter) {

        if (std::abs(mcParticle.eta()) > cfgEta) {
          continue;
        }

        if (std::abs(mcParticle.y()) > rapidityCut) {
          continue;
        }

        qaHistos.fill(HIST("MCGen/ptAntiDeuteronAll"), mcParticle.pt());

        int motherPdg = 0;
        if (mcParticle.has_mothers()) {
          for (const auto& mom : mcParticle.mothers_as<o2::aod::McParticles>()) {

            qaHistos.fill(HIST("MCGen/hMotherPdgCode"), mom.pdgCode());

            if (mom.pdgCode() == pdgCodeBeautyMeson) {
              motherPdg = mom.pdgCode();
              break;
            }

            if (mom.pdgCode() == pdgCodeBeautyBaryon) {
              motherPdg = mom.pdgCode();
              break;
            }

            if (mom.pdgCode() == pdgCodeHypertriton) {
              motherPdg = mom.pdgCode();
              break;
            }
          }
        }

        if (motherPdg == pdgCodeBeautyMeson) {
          qaHistos.fill(HIST("MCGen/ptAntiDeuteronFromBminus"), mcParticle.pt());
        } else if (motherPdg == pdgCodeBeautyBaryon) {
          qaHistos.fill(HIST("MCGen/ptAntiDeuteronFromAntiLambdaB"), mcParticle.pt());
        } else if (motherPdg == pdgCodeHypertriton) {
          qaHistos.fill(HIST("MCGen/ptAntiDeuteronFromAntiHypertriton"), mcParticle.pt());
        } else if (mcParticle.isPhysicalPrimary()) {
          qaHistos.fill(HIST("MCGen/ptAntiDeuteronPrimary"), mcParticle.pt());
        }
      }
    }
  }

  PROCESS_SWITCH(HfTaskDeuteronFromLb, processGen, "processGen", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTaskDeuteronFromLb>(cfgc)};
}
