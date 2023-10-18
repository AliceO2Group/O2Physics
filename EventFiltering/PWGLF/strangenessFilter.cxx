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
/// \brief A filter task for strangeness triggers
/// \author Chiara De Martin (chiara.de.martin@cern.ch)
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \since June 1, 2021

#include <cmath>
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "PWGHF/Core/PDG.h"
#include "../filterTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

struct strangenessFilter {

  // Recall the output table
  Produces<aod::StrangenessFilters> strgtable;
  TrackSelection mTrackSelector;

  // Define a histograms and registries
  HistogramRegistry QAHistos{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry QAHistosTopologicalVariables{"QAHistosTopologicalVariables", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry QAHistosTriggerParticles{"QAHistosTriggerParticles", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry QAHistosStrangenessTracking{"QAHistosStrangenessTracking", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry EventsvsMultiplicity{"EventsvsMultiplicity", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  OutputObj<TH1F> hProcessedEvents{TH1F("hProcessedEvents", "Strangeness - event filtered; Event counter; Number of events", 14, 0., 14.)};
  OutputObj<TH1F> hCandidate{TH1F("hCandidate", "; Candidate pass selection; Number of events", 30, 0., 30.)};
  OutputObj<TH1F> hEvtvshMinPt{TH1F("hEvtvshMinPt", " Number of h-Xi events with pT_h higher than thrd; hadrons with p_{T}>bincenter (GeV/c); Number of events", 11, 0., 11.)};
  OutputObj<TH1F> hhXiPairsvsPt{TH1F("hhXiPairsvsPt", "pt distributions of Xi in events with a trigger particle; #it{p}_{T} (GeV/c); Number of Xi", 100, 0., 10.)};

  // Selection criteria for cascades
  Configurable<bool> doextraQA{"doextraQA", 1, "do extra QA"};
  Configurable<float> cutzvertex{"cutzvertex", 100.0f, "Accepted z-vertex range"};
  Configurable<float> v0cospa{"v0cospa", 0.95, "V0 CosPA"};
  Configurable<float> casccospaxi{"casccospaxi", 0.95, "Casc CosPA"};
  Configurable<float> casccospaomega{"casccospaomega", 0.95, "Casc CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 2.0, "DCA V0 Daughters"};
  Configurable<float> dcacascdau{"dcacascdau", 2.0, "DCA Casc Daughters"};
  Configurable<float> dcamesontopv{"dcamesontopv", 0.05, "DCA Meson To PV"};
  Configurable<float> dcabaryontopv{"dcabaryontopv", 0.05, "DCA Baryon To PV"};
  Configurable<float> dcabachtopv{"dcabachtopv", 0.05, "DCA Bach To PV"};
  Configurable<float> dcav0topv{"dcav0topv", 0.0, "DCA V0 To PV"};
  Configurable<float> v0radius{"v0radius", 1.0, "V0 Radius"};
  Configurable<float> cascradius{"cascradius", 0.6, "cascradius"};
  Configurable<float> rapidity{"rapidity", 2, "rapidity"};
  Configurable<float> eta{"eta", 2, "Eta"};
  Configurable<float> minpt{"minpt", 0.5, "minpt"};
  Configurable<float> etadau{"etadau", 0.9, "EtaDaughters"};
  Configurable<float> masslambdalimit{"masslambdalimit", 0.02, "masslambdalimit"};
  Configurable<float> omegarej{"omegarej", 0.005, "omegarej"};
  Configurable<float> xirej{"xirej", 0.008, "xirej"}; // merge the two rejection variables into one?
  Configurable<float> ximasswindow{"ximasswindow", 0.075, "Xi Mass Window"};
  Configurable<float> omegamasswindow{"omegamasswindow", 0.075, "Omega Mass Window"}; // merge the two windows variables into one?
  Configurable<int> properlifetimefactor{"properlifetimefactor", 5, "Proper Lifetime cut"};
  Configurable<float> lowerradiusXiYN{"lowerradiusXiYN", 24.39, "Cascade lower radius for single Xi trigger"};
  Configurable<float> lowerradiusOmega{"lowerradiusOmega", 19.0, "Omega lower radius for high radius Omega trigger"};
  Configurable<float> upperradiusOmega{"upperradiusOmega", 19.0, "Omega upper radius for low radius Omega trigger"};
  Configurable<float> nsigmatpcpi{"nsigmatpcpi", 6, "N Sigmas TPC pi"};
  Configurable<float> nsigmatpcka{"nsigmatpcka", 6, "N Sigmas TPC ka"};
  Configurable<float> nsigmatpcpr{"nsigmatpcpr", 6, "N Sigmas TPC pr"};
  Configurable<bool> hastof{"hastof", 1, "Has TOF (OOB condition)"};
  Configurable<float> ptthrtof{"ptthrtof", 1.0, "Pt threshold to apply TOF condition"};
  Configurable<bool> kint7{"kint7", 0, "Apply kINT7 event selection"};
  Configurable<bool> sel7{"sel7", 0, "Apply sel7 event selection"};
  Configurable<bool> sel8{"sel8", 0, "Apply sel8 event selection"};

  // Selections criteria for tracks
  Configurable<float> hEta{"hEta", 0.9f, "Eta range for trigger particles"};
  Configurable<float> hMinPt{"hMinPt", 1.0f, "Min pt for trigger particles"};
  Configurable<bool> isTrackFilter{"isTrackFilter", true, "Apply track myTrackSelections"};

  // Settings for strangeness tracking filter
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpMagPath{"grpMagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<double> bz{"bz", -5., "magnetic field"};
  Configurable<float> minPtTrackedCascade{"minPtTrackedCascade", 0., "Min. pt for tracked cascades"};
  Configurable<float> massWindowTrackedOmega{"massWindowTrackedOmega", 0.05, "Inv. mass window for tracked Omega-"};
  Configurable<float> massWindowXiExclTrackedOmega{"massWindowXiExclTrackedOmega", 0.005, "Inv. mass window for exclusion of Xi for tracked Omega-"};
  Configurable<float> massWindowTrackedXi{"massWindowTrackedXi", 0.05, "Inv. mass window for tracked Xi-"};
  Configurable<float> massWindowLambda{"massWindowLambda", 0.05, "Inv. mass window for Lambda (ST)"};
  Configurable<float> maxMatchingChi2TrackedCascade{"maxMatchingChi2TrackedCascade", 2000., "Max matching chi2 for tracked cascades"};

  Configurable<float> maxNSigmaBachelorTrackedXi{"maxNSigmaBachelorTrackedXi", 3., "Max Nsigma for bachelor of tracked Xi (pi)"};
  Configurable<float> maxNSigmaBachelorTrackedOmega{"maxNSigmaBachelorTrackedOmega", 3., "Max Nsigma for bachelor of tracked Xi (Ka)"};
  Configurable<float> maxNSigmaV0PrTrackedCascade{"maxNSigmaV0PrTrackedCascade", 3., "Max Nsigma for proton from V0 fromtracked Xi"};
  Configurable<float> maxNSigmaV0PiTrackedCascade{"maxNSigmaV0PiTrackedCascade", 3., "Max Nsigma for pion from V0 fromtracked Xi"};
  Configurable<float> minPtTrackedV0{"minPtTrackedV0", 0., "Min. pt for tracked V0"};
  Configurable<float> minPtTracked3Body{"minPtTracked3Body", 0., "Min. pt for tracked 3Body"};

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    mTrackSelector.SetTrackType(o2::aod::track::TrackTypeEnum::Track);
    mTrackSelector.SetPtRange(hMinPt, 1e10f);
    mTrackSelector.SetEtaRange(-hEta, hEta);
    mTrackSelector.SetRequireITSRefit(true);
    mTrackSelector.SetRequireTPCRefit(true);
    mTrackSelector.SetRequireGoldenChi2(false);
    mTrackSelector.SetMinNCrossedRowsTPC(70);
    mTrackSelector.SetMinNCrossedRowsOverFindableClustersTPC(0.8f);
    mTrackSelector.SetMaxChi2PerClusterTPC(4.f);
    mTrackSelector.SetRequireHitsInITSLayers(1, {0, 1, 2}); // one hit in any of the first three layers of IB
    mTrackSelector.SetMaxChi2PerClusterITS(36.f);
    // mTrackSelector.SetMaxDcaXYPtDep([](float pt) { return 0.0105f + 0.0350f / pow(pt, 1.1f); });
    mTrackSelector.SetMaxDcaXY(1.f);
    mTrackSelector.SetMaxDcaZ(2.f);

    hProcessedEvents->GetXaxis()->SetBinLabel(1, "Events processed");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "Events w/ high-#it{p}_{T} hadron");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, "#Omega");
    hProcessedEvents->GetXaxis()->SetBinLabel(4, "high-#it{p}_{T} hadron - #Xi");
    hProcessedEvents->GetXaxis()->SetBinLabel(5, "2#Xi");
    hProcessedEvents->GetXaxis()->SetBinLabel(6, "3#Xi");
    hProcessedEvents->GetXaxis()->SetBinLabel(7, "4#Xi");
    hProcessedEvents->GetXaxis()->SetBinLabel(8, "#Xi-YN");
    hProcessedEvents->GetXaxis()->SetBinLabel(9, "#Omega high radius");
    hProcessedEvents->GetXaxis()->SetBinLabel(10, "#Xi");
    hProcessedEvents->GetXaxis()->SetBinLabel(11, "trk. #Xi");
    hProcessedEvents->GetXaxis()->SetBinLabel(12, "trk. #Omega");
    hProcessedEvents->GetXaxis()->SetBinLabel(13, "trk. V^{0}");
    hProcessedEvents->GetXaxis()->SetBinLabel(14, "trk. 3body");

    hCandidate->GetXaxis()->SetBinLabel(1, "All");
    hCandidate->GetXaxis()->SetBinLabel(2, "Has_V0");
    hCandidate->GetXaxis()->SetBinLabel(3, "DCA_meson");
    hCandidate->GetXaxis()->SetBinLabel(4, "DCA_baryon");
    hCandidate->GetXaxis()->SetBinLabel(5, "TPCNsigma_pion");
    hCandidate->GetXaxis()->SetBinLabel(6, "TPCNsigma_proton");
    hCandidate->GetXaxis()->SetBinLabel(7, "Eta_dau");
    hCandidate->GetXaxis()->SetBinLabel(8, "DCABachToPV");
    hCandidate->GetXaxis()->SetBinLabel(9, "V0Radius");
    hCandidate->GetXaxis()->SetBinLabel(10, "CascRadius");
    hCandidate->GetXaxis()->SetBinLabel(11, "V0CosPA");
    hCandidate->GetXaxis()->SetBinLabel(12, "DCAV0Dau");
    hCandidate->GetXaxis()->SetBinLabel(13, "DCACascDau");
    hCandidate->GetXaxis()->SetBinLabel(14, "MassLambdaLimit");
    hCandidate->GetXaxis()->SetBinLabel(15, "Eta");
    hCandidate->GetXaxis()->SetBinLabel(16, "HasTOFOneLeg");
    hCandidate->GetXaxis()->SetBinLabel(17, "CascCosPA");
    hCandidate->GetXaxis()->SetBinLabel(18, "DCAV0ToPV");
    hCandidate->GetXaxis()->SetBinLabel(19, "ProperLifeTime");

    std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};
    AxisSpec multAxisNTPV = {100, 0.0f, 100.0f, "N. tracks PV estimator"};
    AxisSpec multAxisT0M = {600, 0.0f, 6000.0f, "T0M multiplicity estimator"};
    AxisSpec multAxisV0A = {500, 0.0f, 25000.0f, "V0A multiplicity estimator"};
    AxisSpec ximassAxis = {200, 1.28f, 1.36f};
    AxisSpec omegamassAxis = {200, 1.59f, 1.75f};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pTPCAxis = {100, 0.0f, 10.0f, "#it{p} TPC (GeV/#it{c})"};
    AxisSpec etaAxis = {200, -2.0f, 2.0f, "#eta"};
    AxisSpec phiAxis = {100, -TMath::Pi() / 2, 3. * TMath::Pi() / 2, "#varphi"};
    AxisSpec ptTriggAxis = {150, 0.0f, 15.0f, "#it{p}_{T} (GeV/#it{c})"};

    // general QA histograms
    QAHistos.add("hVtxZ", "Z-Vertex distribution after selection;Z (cm)", HistType::kTH1F, {{100, -50, 50}});
    QAHistos.add("hMassXiBefSelvsPt", "hMassXiBefSelvsPt", HistType::kTH2F, {ximassAxis, ptAxis});
    QAHistos.add("hMassOmegaBefSelvsPt", "hMassOmegaBefSelvsPt", HistType::kTH2F, {omegamassAxis, ptAxis});
    QAHistos.add("hMassXiAfterSelvsPt", "hMassXiAfterSelvsPt", HistType::kTH2F, {ximassAxis, ptAxis});
    QAHistos.add("hMassOmegaAfterSelvsPt", "hMassOmegaAfterSelvsPt", HistType::kTH2F, {omegamassAxis, ptAxis});
    QAHistos.add("hPtXi", "pt distribution of selected Xi candidates", HistType::kTH1F, {ptAxis});
    QAHistos.add("hPtOmega", "pt distribution of selected Omega candidates", HistType::kTH1F, {ptAxis});
    QAHistos.add("hEtaXi", "eta distribution of selected Xi candidates", HistType::kTH1F, {etaAxis});
    QAHistos.add("hEtaOmega", "eta distribution of selected Omega candidates", HistType::kTH1F, {etaAxis});

    // topological variables distributions
    QAHistosTopologicalVariables.add("hCascCosPAXi", "hCascCosPAXi", HistType::kTH1F, {{350, 0.65f, 1.0f}});
    QAHistosTopologicalVariables.add("hV0CosPAXi", "hV0CosPAXi", HistType::kTH1F, {{250, 0.75f, 1.0f}});
    QAHistosTopologicalVariables.add("hCascRadiusXi", "hCascRadiusXi", HistType::kTH1F, {{500, 0.0f, 50.0f}});
    QAHistosTopologicalVariables.add("hV0RadiusXi", "hV0RadiusXi", HistType::kTH1F, {{500, 0.0f, 50.0f}});
    QAHistosTopologicalVariables.add("hDCAV0DaughtersXi", "hDCAV0DaughtersXi", HistType::kTH1F, {{110, 0.0f, 2.2f}});
    QAHistosTopologicalVariables.add("hDCACascDaughtersXi", "hDCACascDaughtersXi", HistType::kTH1F, {{110, 0.0f, 2.2f}});
    QAHistosTopologicalVariables.add("hDCAV0ToPVXi", "hDCAV0ToPVXi", HistType::kTH1F, {{220, 0.0f, 2.2f}});
    QAHistosTopologicalVariables.add("hDCABachToPVXi", "|hDCABachToPVXi|", HistType::kTH1F, {{400, 0.0f, 2.0f}});
    QAHistosTopologicalVariables.add("hDCAPosToPVXi", "|hDCAPosToPVXi|", HistType::kTH1F, {{400, 0.0f, 2.0f}});
    QAHistosTopologicalVariables.add("hDCANegToPVXi", "|hDCANegToPVXi|", HistType::kTH1F, {{400, 0.0f, 2.0f}});
    QAHistosTopologicalVariables.add("hInvMassLambdaXi", "InvMassLambdaXi", HistType::kTH1F, {{200, 1.07f, 1.17f}});
    QAHistosTopologicalVariables.add("hProperLifetimeXi", "Proper Lifetime Xi", HistType::kTH1F, {{50, 0, 50}});
    //
    QAHistosTopologicalVariables.add("hCascCosPAOmega", "hCascCosPAOmega", HistType::kTH1F, {{350, 0.65f, 1.0f}});
    QAHistosTopologicalVariables.add("hV0CosPAOmega", "hV0CosPAOmega", HistType::kTH1F, {{250, 0.75f, 1.0f}});
    QAHistosTopologicalVariables.add("hCascRadiusOmega", "hCascRadiusOmega", HistType::kTH1F, {{500, 0.0f, 50.0f}});
    QAHistosTopologicalVariables.add("hV0RadiusOmega", "hV0RadiusOmega", HistType::kTH1F, {{500, 0.0f, 50.0f}});
    QAHistosTopologicalVariables.add("hDCAV0DaughtersOmega", "hDCAV0DaughtersOmega", HistType::kTH1F, {{110, 0.0f, 2.2f}});
    QAHistosTopologicalVariables.add("hDCACascDaughtersOmega", "hDCACascDaughtersOmega", HistType::kTH1F, {{110, 0.0f, 2.2f}});
    QAHistosTopologicalVariables.add("hDCAV0ToPVOmega", "hDCAV0ToPVOmega", HistType::kTH1F, {{220, 0.0f, 2.2f}});
    QAHistosTopologicalVariables.add("hDCABachToPVOmega", "hDCABachToPVOmega", HistType::kTH1F, {{400, 0.0f, 2.0f}});
    QAHistosTopologicalVariables.add("hDCAPosToPVOmega", "hDCAPosToPVOmega", HistType::kTH1F, {{400, 0.0f, 2.0f}});
    QAHistosTopologicalVariables.add("hDCANegToPVOmega", "hDCANegToPVOmega", HistType::kTH1F, {{400, 0.0f, 2.0f}});
    QAHistosTopologicalVariables.add("hInvMassLambdaOmega", "InvMassLambdaOmega", HistType::kTH1F, {{200, 1.07f, 1.17f}});
    QAHistosTopologicalVariables.add("hProperLifetimeOmega", "Proper Lifetime Omega", HistType::kTH1F, {{50, 0, 50}});
    QAHistosTopologicalVariables.add("hCascRadiusOmegaLargeR", "hCascRadiusOmegaLargeR", HistType::kTH1F, {{500, 0.0f, 50.0f}});
    QAHistosTopologicalVariables.add("hCascRadiusXiYN", "hCascRadiusXiYN", HistType::kTH1F, {{500, 0.0f, 50.0f}});

    // trigger particles QA
    QAHistosTriggerParticles.add("hTriggeredParticlesAllEv", "Distribution of #tracks w/ pt > pt,trigg,min", HistType::kTH1F, {{20, 0.5, 20.5, "Trigger counter"}});
    QAHistosTriggerParticles.add("hTriggeredParticlesSelEv", "Distribution of #tracks w/ pt > pt,trigg,min (after sel)", HistType::kTH1F, {{20, 0.5, 20.5, "Trigger counter"}});
    QAHistosTriggerParticles.add("hPtTriggerAllEv", "hPtTriggerAllEv", HistType::kTH1F, {{300, 0, 30, "Pt of trigger particles"}});
    QAHistosTriggerParticles.add("hPtTriggerSelEv", "hPtTriggerSelEv", HistType::kTH1F, {{300, 0, 30, "Pt of trigger particles after selections"}});
    QAHistosTriggerParticles.add("hEtaTriggerAllEv", "hEtaTriggerAllEv", HistType::kTH2F, {{180, -1.4, 1.4, "Eta of trigger particles"}, {ptTriggAxis}});
    QAHistosTriggerParticles.add("hPhiTriggerAllEv", "hPhiTriggerAllEv", HistType::kTH2F, {{100, 0, 2 * TMath::Pi(), "Phi of trigger particles"}, {ptTriggAxis}});
    QAHistosTriggerParticles.add("hDCAxyTriggerAllEv", "hDCAxyTriggerAllEv", HistType::kTH2F, {{400, -0.2, 0.2, "DCAxy of trigger particles"}, {ptTriggAxis}});
    QAHistosTriggerParticles.add("hDCAzTriggerAllEv", "hDCAzTriggerAllEv", HistType::kTH2F, {{400, -0.2, 0.2, "DCAz of trigger particles"}, {ptTriggAxis}});

    if (doextraQA) {
      EventsvsMultiplicity.add("AllEventsvsMultiplicityZeqV0A", "ZeqV0A distribution of all events", HistType::kTH1F, {multAxisV0A});
      EventsvsMultiplicity.add("hadEventsvsMultiplicityZeqV0A", "ZeqV0A distribution of events with hight pT hadron", HistType::kTH1F, {multAxisV0A});
      EventsvsMultiplicity.add("hadEventsvsMultiplicityZeqV0AvsPt", "ZeqV0A distribution of events with hight pT hadron", HistType::kTH2F, {{multAxisV0A}, {11, 0, 11}});

      EventsvsMultiplicity.add("AllEventsvsMultiplicityZeqT0M", "ZeqT0M distribution of all events", HistType::kTH1F, {multAxisT0M});
      EventsvsMultiplicity.add("hadEventsvsMultiplicityZeqT0M", "ZeqT0M distribution of events with hight pT hadron", HistType::kTH1F, {multAxisT0M});
      EventsvsMultiplicity.add("hadEventsvsMultiplicityZeqT0MvsPt", "ZeqT0M distribution of events with hight pT hadron", HistType::kTH2F, {{multAxisT0M}, {11, 0, 11}});

      EventsvsMultiplicity.add("AllEventsvsMultiplicityZeqNTracksPV", "ZeqNTracksPV distribution of all events", HistType::kTH1F, {multAxisNTPV});
      EventsvsMultiplicity.add("hadEventsvsMultiplicityZeqNTracksPV", "ZeqNTracksPV distribution of events with hight pT hadron", HistType::kTH1F, {multAxisNTPV});
      EventsvsMultiplicity.add("hadEventsvsMultiplicityZeqNTracksPVvsPt", "ZeqNTracksPV distribution of events with hight pT hadron", HistType::kTH2F, {{multAxisNTPV}, {11, 0, 11}});

      // additional QA histos
      QAHistos.add("hTPCNsigmaXiBachPiPlus", "nsigma TPC distribution bachelor pion+", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaXiV0PiPlus", "nsigma TPC distribution pi+", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaXiV0Proton", "nsigma TPC distribution proton", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaXiBachPiMinus", "nsigma TPC distribution bachelor pion-", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaXiV0PiMinus", "nsigma TPC distribution pi-", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaXiV0AntiProton", "nsigma TPC distribution antiproton", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaOmegaBachKaPlus", "nsigma TPC distribution bachelor kaon+", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaOmegaV0PiPlus", "nsigma TPC distribution pi+", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaOmegaV0Proton", "nsigma TPC distribution proton", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaOmegaBachKaMinus", "nsigma TPC distribution bachelor kaon-", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaOmegaV0PiMinus", "nsigma TPC distribution pi-", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaOmegaV0AntiProton", "nsigma TPC distribution antiproton", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hHasTOFBachKa", "bachelor kaon has TOF", HistType::kTH2F, {{2, 0, 2}, {ptAxis}});
      QAHistos.add("hHasTOFBachPi", "bachelor pi has TOF", HistType::kTH2F, {{2, 0, 2}, {ptAxis}});
      QAHistos.add("hHasTOFPr", "pr dau has TOF", HistType::kTH2F, {{2, 0, 2}, {ptAxis}});
      QAHistos.add("hHasTOFPi", "pi dau has TOF", HistType::kTH2F, {{2, 0, 2}, {ptAxis}});
      QAHistos.add("hRapXi", "Rap Xi", HistType::kTH1F, {{100, -1, 1}});
      QAHistos.add("hRapOmega", "Rap Omega", HistType::kTH1F, {{100, -1, 1}});

      QAHistosStrangenessTracking.add("hStRVsPtTrkCasc", "Tracked cascades;p_{T} (GeV/#it{c});R (cm)", HistType::kTH2D, {{200, 0., 10.}, {200, 0., 50}});
      QAHistosStrangenessTracking.add("hMassOmegaTrkCasc", "Tracked cascades;m_{#Omega} (GeV/#it{c}^{2})", HistType::kTH1D, {{1000, 1., 3.}});
      QAHistosStrangenessTracking.add("hMassXiTrkCasc", "Tracked cascades;m_{#Xi} (GeV/#it{c}^{2})", HistType::kTH1D, {{1000, 1., 3.}});
      QAHistosStrangenessTracking.add("hMassV0TrkCasc", "Tracked cascades;m_{V^{0}} (GeV/#it{c}^{2})", HistType::kTH1D, {{1000, 1., 3.}});
      QAHistosStrangenessTracking.add("hMatchChi2TrkCasc", "Tracked cascades;#chi^{2}", HistType::kTH1D, {{1000, 0., 2000.}});
      QAHistosStrangenessTracking.add("hMassOmegaVsMatchChi2TrkCasc", "Tracked cascades;m_{#Omega} (GeV/#it{c}^{2});#chi^{2}", HistType::kTH2D, {{1000, 1., 3.}, {1000, 0., 2000.}});
      QAHistosStrangenessTracking.add("hMassXiVsMatchChi2TrkCasc", "Tracked cascades;m_{#Xi} (GeV/#it{c}^{2});#chi^{2}", HistType::kTH2D, {{1000, 1., 3.}, {1000, 0., 2000.}});
      QAHistosStrangenessTracking.add("hMassOmegaVsTopChi2TrkCasc", "Tracked cascades;m_{#Omega} (GeV/#it{c}^{2});#chi^{2}", HistType::kTH2D, {{1000, 1., 3.}, {1000, 0., 2000.}});
      QAHistosStrangenessTracking.add("hMassXiVsTopChi2TrkCasc", "Tracked cascades;m_{#Xi} (GeV/#it{c}^{2});#chi^{2}", HistType::kTH2D, {{1000, 1., 3.}, {1000, 0., 2000.}});
      QAHistosStrangenessTracking.add("hNSigmaTpcPiTrkCascBachelor", "Tracked cascades;N_{#sigma, #pi}", HistType::kTH1D, {{100, -5., 5.}});
      QAHistosStrangenessTracking.add("hNSigmaTpcKaTrkCascBachelor", "Tracked cascades;N_{#sigma, K}", HistType::kTH1D, {{100, -5., 5.}});
      QAHistosStrangenessTracking.add("hNSigmaTpcPrTrkCascV0", "Tracked cascades;N_{#sigma, p}", HistType::kTH1D, {{100, -5., 5.}});
      QAHistosStrangenessTracking.add("hNSigmaTpcPiTrkCascV0", "Tracked cascades;N_{#sigma, #pi}", HistType::kTH1D, {{100, -5., 5.}});
      QAHistosStrangenessTracking.add("hMassH3LTrkV0", "Tracked V0;m_{H3L} (GeV/#it{c}^{2})", HistType::kTH1D, {{1000, 2.8, 3.8}});
      QAHistosStrangenessTracking.add("hMassH4LTrkV0", "Tracked V0;m_{H4L} (GeV/#it{c}^{2})", HistType::kTH1D, {{1000, 3.8, 4.8}});
      QAHistosStrangenessTracking.add("hMassH3LTrk3body", "Tracked 3body;m_{H3L} (GeV/#it{c}^{2})", HistType::kTH1D, {{200, 0., 10.}});
      QAHistosStrangenessTracking.add("hMassHe4LTrk3body", "Tracked 3body;m_{He4L} (GeV/#it{c}^{2})", HistType::kTH1D, {{200, 0., 10.}});
      QAHistosStrangenessTracking.add("hDcaXY", "DCA;DCA_{xy} (cm)", HistType::kTH1D, {{200, -.5, .5}});
      QAHistosStrangenessTracking.add("hDcaXYVsPt", "DCA;p_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2D, {{200, 0., 10.}, {200, -.5, .5}});
      QAHistosStrangenessTracking.add("hDcaZ", "DCA;DCA_{z} (cm)", HistType::kTH1D, {{200, -.5, .5}});
      QAHistosStrangenessTracking.add("hDcaZVsPt", "DCA;p_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2D, {{200, 0., 10.}, {200, -.5, .5}});
      QAHistosStrangenessTracking.add("hDcaVsPt", "DCA;DCA (cm);p_{T} (GeV/#it{c})", HistType::kTH2D, {{200, 0., .5}, {200, 0., 10.}});
      QAHistosStrangenessTracking.add("hDcaVsR", "DCA;DCA (cm);R (cm)", HistType::kTH2D, {{200, 0., .5}, {200, 0., 10.}});
      QAHistosStrangenessTracking.add("hPtCascCand", "cascades;p_{T} (GeV/#it{c})", HistType::kTH1D, {{200, 0., 10.}});
      QAHistosStrangenessTracking.add("hPtCascTracked", "tracked cascades;p_{T} (GeV/#it{c})", HistType::kTH1D, {{200, 0., 10.}});
    }
  }

  // Filters
  Filter trackFilter = (nabs(aod::track::eta) < hEta) && (aod::track::pt > hMinPt);

  // Tables
  using CollisionCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>::iterator;
  // using CollisionCandidatesRun3 = soa::Join<aod::Collisions, aod::EvSels>::iterator;
  using CollisionCandidatesRun3 = soa::Join<aod::Collisions, aod::EvSels, aod::MultZeqs>::iterator;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA>>;
  using DaughterTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCLfFullPi, aod::pidTPCLfFullPr, aod::pidTPCLfFullKa>;
  using Cascades = aod::CascDataExt;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int runNumber;

  ////////////////////////////////////////////////////////
  ////////// Strangeness Filter - Run 2 conv /////////////
  ////////////////////////////////////////////////////////

  void fillTriggerTable(bool keepEvent[])
  {
    strgtable(keepEvent[0], keepEvent[1], keepEvent[2], keepEvent[3], keepEvent[4], keepEvent[5], keepEvent[6], keepEvent[7], keepEvent[8], keepEvent[9], keepEvent[10]);
  }

  void processRun2(CollisionCandidates const& collision, TrackCandidates const& tracks, Cascades const& fullCasc, aod::V0sLinked const&, aod::V0Datas const& v0data, DaughterTracks& dtracks)
  {
    // Is event good? [0] = Omega, [1] = high-pT hadron + Xi, [2] = 2Xi, [3] = 3Xi, [4] = 4Xi, [5] single-Xi, [6] Omega with high radius
    // [7] tracked Xi, [8] tracked Omega, [9] tracked V0, [10] tracked 3Body
    bool keepEvent[11]{}; // explicitly zero-initialised

    if (kint7 && !collision.alias_bit(kINT7)) {
      fillTriggerTable(keepEvent);
      return;
    }
    if (sel7 && !collision.sel7()) {
      fillTriggerTable(keepEvent);
      return;
    }
    if (sel8 && !collision.sel8()) {
      fillTriggerTable(keepEvent);
      return;
    }

    if (TMath::Abs(collision.posZ()) > cutzvertex) {
      fillTriggerTable(keepEvent);
      return;
    }

    if (doextraQA) {
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicity"), collision.centRun2V0M());
      QAHistos.fill(HIST("hCentrality"), collision.centRun2V0M());
    }
    hProcessedEvents->Fill(0.5);

    // constants
    const float ctauxi = 4.91;     // from PDG
    const float ctauomega = 2.461; // from PDG

    // variables
    float xipos = -1.;
    float xiproperlifetime = -1.;
    float omegaproperlifetime = -1.;
    float xiptotmom = -1.;
    int xicounter = 0;
    int xicounterYN = 0;
    int omegacounter = 0;
    int omegalargeRcounter = 0;
    int triggcounterForEstimates = 0;
    int triggcounter = 0;

    for (auto& casc : fullCasc) { // loop over cascades
      triggcounterForEstimates = 0;
      auto v0index = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0index.has_v0Data())) {
        continue; // skip those cascades for which V0 doesn't exist
      }
      auto v0 = v0index.v0Data(); // de-reference index to correct v0data in case it exists
      auto bachelor = casc.bachelor_as<DaughterTracks>();
      auto posdau = v0.posTrack_as<DaughterTracks>();
      auto negdau = v0.negTrack_as<DaughterTracks>();

      bool isXi = false;
      bool isXiYN = false;
      bool isOmega = false;
      bool isOmegalargeR = false;

      // Position
      xipos = std::hypot(casc.x() - collision.posX(), casc.y() - collision.posY(), casc.z() - collision.posZ());
      // Total momentum
      xiptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
      // Proper lifetime
      xiproperlifetime = o2::analysis::pdg::MassXiMinus * xipos / (xiptotmom + 1e-13);
      omegaproperlifetime = o2::analysis::pdg::MassOmegaMinus * xipos / (xiptotmom + 1e-13);

      if (casc.sign() == 1) {
        if (TMath::Abs(casc.dcapostopv()) < dcamesontopv) {
          continue;
        }
        if (TMath::Abs(casc.dcanegtopv()) < dcabaryontopv) {
          continue;
        }
        if (TMath::Abs(posdau.tpcNSigmaPi()) > nsigmatpcpi) {
          continue;
        }
        if (TMath::Abs(negdau.tpcNSigmaPr()) > nsigmatpcpr) {
          continue;
        }
      } else {
        if (TMath::Abs(casc.dcanegtopv()) < dcamesontopv) {
          continue;
        }
        if (TMath::Abs(casc.dcapostopv()) < dcabaryontopv) {
          continue;
        }
        if (TMath::Abs(posdau.tpcNSigmaPr()) > nsigmatpcpr) {
          continue;
        }
        if (TMath::Abs(negdau.tpcNSigmaPi()) > nsigmatpcpi) {
          continue;
        }
      }
      // these selection differ for Xi and Omegas:
      if (TMath::Abs(posdau.eta()) > etadau) {
        continue;
      }
      if (TMath::Abs(negdau.eta()) > etadau) {
        continue;
      }
      if (TMath::Abs(bachelor.eta()) > etadau) {
        continue;
      }
      if (TMath::Abs(casc.dcabachtopv()) < dcabachtopv) {
        continue;
      }
      if (casc.v0radius() < v0radius) {
        continue;
      }
      if (casc.cascradius() < cascradius) {
        continue;
      }
      if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa) {
        continue;
      }
      if (casc.dcaV0daughters() > dcav0dau) {
        continue;
      }
      if (casc.dcacascdaughters() > dcacascdau) {
        continue;
      }
      if (TMath::Abs(casc.mLambda() - constants::physics::MassLambda) > masslambdalimit) {
        continue;
      }
      if (TMath::Abs(casc.eta()) > eta) {
        continue;
      }

      isXi = (TMath::Abs(bachelor.tpcNSigmaPi()) < nsigmatpcpi) &&
             (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > casccospaxi) &&
             (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv) &&
             (TMath::Abs(casc.mXi() - o2::analysis::pdg::MassXiMinus) < ximasswindow) &&
             (TMath::Abs(casc.mOmega() - o2::analysis::pdg::MassOmegaMinus) > omegarej) &&
             (xiproperlifetime < properlifetimefactor * ctauxi) &&
             (TMath::Abs(casc.yXi()) < rapidity); // add PID on bachelor
      isXiYN = (TMath::Abs(bachelor.tpcNSigmaPi()) < nsigmatpcpi) &&
               (casc.cascradius() > lowerradiusXiYN) &&
               (TMath::Abs(casc.mXi() - o2::analysis::pdg::MassXiMinus) < ximasswindow) &&
               (TMath::Abs(casc.mOmega() - o2::analysis::pdg::MassOmegaMinus) > omegarej) &&
               (xiproperlifetime < properlifetimefactor * ctauxi) &&
               (TMath::Abs(casc.yXi()) < rapidity); // add PID on bachelor
      isOmega = (TMath::Abs(bachelor.tpcNSigmaKa()) < nsigmatpcka) &&
                (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > casccospaomega) &&
                (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv) &&
                (casc.cascradius() < upperradiusOmega) &&
                (TMath::Abs(casc.mOmega() - o2::analysis::pdg::MassOmegaMinus) < omegamasswindow) &&
                (TMath::Abs(casc.mXi() - o2::analysis::pdg::MassXiMinus) > xirej) &&
                (omegaproperlifetime < properlifetimefactor * ctauomega) &&
                (TMath::Abs(casc.yOmega()) < rapidity); // add PID on bachelor
      isOmegalargeR = (TMath::Abs(bachelor.tpcNSigmaKa()) < nsigmatpcka) &&
                      (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > casccospaomega) &&
                      (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv) &&
                      (casc.cascradius() > lowerradiusOmega) &&
                      (TMath::Abs(casc.mOmega() - o2::analysis::pdg::MassOmegaMinus) < omegamasswindow) &&
                      (TMath::Abs(casc.mXi() - o2::analysis::pdg::MassXiMinus) > xirej) &&
                      (omegaproperlifetime < properlifetimefactor * ctauomega) &&
                      (TMath::Abs(casc.yOmega()) < rapidity); // add PID on bachelor

      if (isXi) {
        // Count number of Xi candidates
        xicounter++;

        // Plot for estimates
        if (tracks.size() > 0)
          triggcounterForEstimates = 1;
        if (triggcounterForEstimates && (TMath::Abs(casc.mXi() - o2::analysis::pdg::MassXiMinus) < 0.01))
          hhXiPairsvsPt->Fill(casc.pt()); // Fill the histogram with all the Xis produced in events with a trigger particle
        // End plot for estimates
      }
      if (isXiYN) {
        // Xis for YN interactions
        xicounterYN++;
      }
      if (isOmega) {
        // Count number of Omega candidates
        omegacounter++;
      }
      if (isOmegalargeR) {
        // Count number of Omega candidates with high radius
        omegalargeRcounter++;
      }
    } // end loop over cascades

    // Omega trigger definition
    if (omegacounter > 0) {
      keepEvent[0] = true;
    }

    // High-pT hadron + Xi trigger definition
    if (xicounter > 0) {
      for (auto track : tracks) { // start loop over tracks
        if (isTrackFilter && !mTrackSelector.IsSelected(track)) {
          continue;
        }
        triggcounter++;
        keepEvent[1] = true;
      } // end loop over tracks
    }

    // 2Xi trigger definition
    if (xicounter > 1) {
      keepEvent[2] = true;
    }

    // 3Xi trigger definition
    if (xicounter > 2) {
      keepEvent[3] = true;
    }

    // 4Xi trigger definition
    if (xicounter > 3) {
      keepEvent[4] = true;
    }

    // Single-Xi (YN) trigger definition
    if (xicounterYN > 0) {
      keepEvent[5] = true;
    }

    // Omega with high radius trigger definition
    if (omegalargeRcounter > 0) {
      keepEvent[6] = true;
    }

    // Fill centrality dependent histos
    if (keepEvent[0]) {
      hProcessedEvents->Fill(2.5);
    }
    if (keepEvent[1]) {
      hProcessedEvents->Fill(3.5);
    }
    if (keepEvent[2]) {
      hProcessedEvents->Fill(4.5);
    }
    if (keepEvent[3]) {
      hProcessedEvents->Fill(5.5);
    }
    if (keepEvent[4]) {
      hProcessedEvents->Fill(6.5);
    }
    if (keepEvent[5]) {
      hProcessedEvents->Fill(7.5);
    }
    if (keepEvent[6]) {
      hProcessedEvents->Fill(8.5);
    }

    // Filling the table
    fillTriggerTable(keepEvent);
  }
  //
  PROCESS_SWITCH(strangenessFilter, processRun2, "Process data Run2", true);

  //////////////////////////////////////////////////////
  ////////// Strangeness Filter - Run 3 MC /////////////
  //////////////////////////////////////////////////////

  void processRun3(CollisionCandidatesRun3 const& collision, TrackCandidates const& tracks, Cascades const& fullCasc, aod::V0sLinked const&, aod::V0Datas const& v0data, DaughterTracks& dtracks,
                   aod::AssignedTrackedCascades const& trackedCascades, aod::Cascades const& cascades, aod::AssignedTrackedV0s const& trackedV0s, aod::AssignedTracked3Bodys const& tracked3Bodys, aod::BCsWithTimestamps const&)
  {
    // Is event good? [0] = Omega, [1] = high-pT hadron + Xi, [2] = 2Xi, [3] = 3Xi, [4] = 4Xi, [5] single-Xi, [6] Omega with high radius
    // [7] tracked Xi, [8] tracked Omega, [9] tracked V0, [10] tracked 3Body
    bool keepEvent[11]{}; // explicitly zero-initialised

    if (sel8 && !collision.sel8()) {
      fillTriggerTable(keepEvent);
      return;
    }
    // all processed events after event selection
    hProcessedEvents->Fill(0.5);

    if (TMath::Abs(collision.posZ()) > cutzvertex) {
      fillTriggerTable(keepEvent);
      return;
    }
    QAHistos.fill(HIST("hVtxZ"), collision.posZ());
    if (doextraQA) {
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityZeqV0A"), collision.multZeqFV0A());
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityZeqT0M"), collision.multZeqFT0A() + collision.multZeqFT0C());
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityZeqNTracksPV"), collision.multZeqNTracksPV());
    }

    // constants
    const float ctauxi = 4.91;     // from PDG
    const float ctauomega = 2.461; // from PDG

    // variables
    float xipos = -1.;
    float xiproperlifetime = -1.;
    float omegaproperlifetime = -1.;
    float xiptotmom = -1.;
    int xicounter = 0;
    int xicounterYN = 0;
    int omegacounter = 0;
    int omegalargeRcounter = 0;
    int triggcounter = 0;
    int triggcounterAllEv = 0;
    int triggcounterForEstimates = 0;

    for (auto& casc : fullCasc) { // loop over cascades
      triggcounterForEstimates = 0;

      hCandidate->Fill(0.5); // All candidates

      auto v0index = casc.v0_as<o2::aod::V0sLinked>();
      if (!(v0index.has_v0Data())) {
        continue; // skip those cascades for which V0 doesn't exist
      }
      hCandidate->Fill(1.5);      // V0 exists
      auto v0 = v0index.v0Data(); // de-reference index to correct v0data in case it exists
      auto bachelor = casc.bachelor_as<DaughterTracks>();
      auto posdau = v0.posTrack_as<DaughterTracks>();
      auto negdau = v0.negTrack_as<DaughterTracks>();

      bool isXi = false;
      bool isXiYN = false;
      bool isOmega = false;
      bool isOmegalargeR = false;

      // QA
      QAHistos.fill(HIST("hMassXiBefSelvsPt"), casc.mXi(), casc.pt());
      QAHistos.fill(HIST("hMassOmegaBefSelvsPt"), casc.mOmega(), casc.pt());

      // Position
      xipos = std::hypot(casc.x() - collision.posX(), casc.y() - collision.posY(), casc.z() - collision.posZ());
      // Total momentum
      xiptotmom = std::hypot(casc.px(), casc.py(), casc.pz());
      // Proper lifetime
      xiproperlifetime = o2::analysis::pdg::MassXiMinus * xipos / (xiptotmom + 1e-13);
      omegaproperlifetime = o2::analysis::pdg::MassOmegaMinus * xipos / (xiptotmom + 1e-13);

      if (casc.sign() > 0) {
        if (TMath::Abs(casc.dcapostopv()) < dcamesontopv) {
          continue;
        }
        hCandidate->Fill(2.5);
        if (TMath::Abs(casc.dcanegtopv()) < dcabaryontopv) {
          continue;
        }
        hCandidate->Fill(3.5);
        if (TMath::Abs(posdau.tpcNSigmaPi()) > nsigmatpcpi) {
          continue;
        }
        hCandidate->Fill(4.5);
        if (TMath::Abs(negdau.tpcNSigmaPr()) > nsigmatpcpr) {
          continue;
        }
        hCandidate->Fill(5.5);
      } else if (casc.sign() < 0) {
        if (TMath::Abs(casc.dcanegtopv()) < dcamesontopv) {
          continue;
        }
        hCandidate->Fill(2.5);
        if (TMath::Abs(casc.dcapostopv()) < dcabaryontopv) {
          continue;
        }
        hCandidate->Fill(3.5);
        if (TMath::Abs(negdau.tpcNSigmaPi()) > nsigmatpcpi) {
          continue;
        }
        hCandidate->Fill(4.5);
        if (TMath::Abs(posdau.tpcNSigmaPr()) > nsigmatpcpr) {
          continue;
        }
        hCandidate->Fill(5.5);
      }
      if (TMath::Abs(posdau.eta()) > etadau) {
        continue;
      }
      if (TMath::Abs(negdau.eta()) > etadau) {
        continue;
      }
      if (TMath::Abs(bachelor.eta()) > etadau) {
        continue;
      }
      hCandidate->Fill(6.5);
      if (TMath::Abs(casc.dcabachtopv()) < dcabachtopv) {
        continue;
      }
      hCandidate->Fill(7.5);
      if (casc.v0radius() < v0radius) {
        continue;
      }
      hCandidate->Fill(8.5);
      if (casc.cascradius() < cascradius) {
        continue;
      }
      hCandidate->Fill(9.5);
      if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < v0cospa) {
        continue;
      }
      hCandidate->Fill(10.5);
      if (casc.dcaV0daughters() > dcav0dau) {
        continue;
      }
      hCandidate->Fill(11.5);
      if (casc.dcacascdaughters() > dcacascdau) {
        continue;
      }
      hCandidate->Fill(12.5);
      if (TMath::Abs(casc.mLambda() - constants::physics::MassLambda) > masslambdalimit) {
        continue;
      }
      hCandidate->Fill(13.5);
      if (TMath::Abs(casc.eta()) > eta) {
        continue;
      }
      hCandidate->Fill(14.5);
      if (hastof &&
          (!posdau.hasTOF() && posdau.pt() > ptthrtof) &&
          (!negdau.hasTOF() && negdau.pt() > ptthrtof) &&
          (!bachelor.hasTOF() && bachelor.pt() > ptthrtof)) {
        continue;
      }
      hCandidate->Fill(15.5);

      // Fill selections QA for XiMinus
      if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > casccospaxi) {
        hCandidate->Fill(16.5);
        if (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv) {
          hCandidate->Fill(17.5);
          if (xiproperlifetime < properlifetimefactor * ctauxi) {
            hCandidate->Fill(18.5);
            if (TMath::Abs(casc.yXi()) < rapidity) {
              hCandidate->Fill(19.5);
            }
          }
        }
      }

      isXi = (TMath::Abs(bachelor.tpcNSigmaPi()) < nsigmatpcpi) &&
             (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > casccospaxi) &&
             (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv) &&
             (TMath::Abs(casc.mXi() - o2::analysis::pdg::MassXiMinus) < ximasswindow) &&
             (TMath::Abs(casc.mOmega() - o2::analysis::pdg::MassOmegaMinus) > omegarej) &&
             (xiproperlifetime < properlifetimefactor * ctauxi) &&
             (TMath::Abs(casc.yXi()) < rapidity);
      isXiYN = (TMath::Abs(bachelor.tpcNSigmaPi()) < nsigmatpcpi) &&
               (casc.cascradius() > lowerradiusXiYN) &&
               (TMath::Abs(casc.mXi() - o2::analysis::pdg::MassXiMinus) < ximasswindow) &&
               (TMath::Abs(casc.mOmega() - o2::analysis::pdg::MassOmegaMinus) > omegarej) &&
               (xiproperlifetime < properlifetimefactor * ctauxi) &&
               (TMath::Abs(casc.yXi()) < rapidity);
      isOmega = (TMath::Abs(bachelor.tpcNSigmaKa()) < nsigmatpcka) &&
                (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > casccospaomega) &&
                (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv) &&
                (TMath::Abs(casc.mOmega() - o2::analysis::pdg::MassOmegaMinus) < omegamasswindow) &&
                (TMath::Abs(casc.mXi() - o2::analysis::pdg::MassXiMinus) > xirej) &&
                (casc.cascradius() < upperradiusOmega) &&
                (omegaproperlifetime < properlifetimefactor * ctauomega) &&
                (TMath::Abs(casc.yOmega()) < rapidity);
      isOmegalargeR = (TMath::Abs(bachelor.tpcNSigmaKa()) < nsigmatpcka) &&
                      (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) > casccospaomega) &&
                      (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) > dcav0topv) &&
                      (casc.cascradius() > lowerradiusOmega) &&
                      (TMath::Abs(casc.mOmega() - o2::analysis::pdg::MassOmegaMinus) < omegamasswindow) &&
                      (TMath::Abs(casc.mXi() - o2::analysis::pdg::MassXiMinus) > xirej) &&
                      (omegaproperlifetime < properlifetimefactor * ctauomega) &&
                      (TMath::Abs(casc.yOmega()) < rapidity);

      if (isXi) {
        QAHistos.fill(HIST("hMassXiAfterSelvsPt"), casc.mXi(), casc.pt());
        QAHistos.fill(HIST("hPtXi"), casc.pt());
        QAHistos.fill(HIST("hEtaXi"), casc.eta());
        QAHistosTopologicalVariables.fill(HIST("hProperLifetimeXi"), xiproperlifetime);
        QAHistosTopologicalVariables.fill(HIST("hCascCosPAXi"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
        QAHistosTopologicalVariables.fill(HIST("hV0CosPAXi"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
        QAHistosTopologicalVariables.fill(HIST("hCascRadiusXi"), casc.cascradius());
        QAHistosTopologicalVariables.fill(HIST("hV0RadiusXi"), casc.v0radius());
        QAHistosTopologicalVariables.fill(HIST("hDCAV0ToPVXi"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
        QAHistosTopologicalVariables.fill(HIST("hDCAV0DaughtersXi"), casc.dcaV0daughters());
        QAHistosTopologicalVariables.fill(HIST("hDCACascDaughtersXi"), casc.dcacascdaughters());
        QAHistosTopologicalVariables.fill(HIST("hDCABachToPVXi"), TMath::Abs(casc.dcabachtopv()));
        QAHistosTopologicalVariables.fill(HIST("hDCAPosToPVXi"), TMath::Abs(casc.dcapostopv()));
        QAHistosTopologicalVariables.fill(HIST("hDCANegToPVXi"), TMath::Abs(casc.dcanegtopv()));
        QAHistosTopologicalVariables.fill(HIST("hInvMassLambdaXi"), casc.mLambda());

        if (doextraQA) {

          QAHistos.fill(HIST("hHasTOFBachPi"), bachelor.hasTOF(), bachelor.pt());
          // QA PID
          if (casc.sign() > 0) {
            QAHistos.fill(HIST("hTPCNsigmaXiBachPiPlus"), bachelor.tpcNSigmaPi(), bachelor.tpcInnerParam());
            QAHistos.fill(HIST("hTPCNsigmaXiV0PiPlus"), posdau.tpcNSigmaPi(), posdau.tpcInnerParam());
            QAHistos.fill(HIST("hTPCNsigmaXiV0AntiProton"), negdau.tpcNSigmaPr(), negdau.tpcInnerParam());
            QAHistos.fill(HIST("hHasTOFPi"), posdau.hasTOF(), posdau.pt());
            QAHistos.fill(HIST("hHasTOFPr"), negdau.hasTOF(), negdau.pt());
          } else {
            QAHistos.fill(HIST("hTPCNsigmaXiBachPiMinus"), bachelor.tpcNSigmaPi(), bachelor.tpcInnerParam());
            QAHistos.fill(HIST("hTPCNsigmaXiV0Proton"), posdau.tpcNSigmaPr(), posdau.tpcInnerParam());
            QAHistos.fill(HIST("hTPCNsigmaXiV0PiMinus"), negdau.tpcNSigmaPi(), negdau.tpcInnerParam());
            QAHistos.fill(HIST("hHasTOFPr"), posdau.hasTOF(), posdau.pt());
            QAHistos.fill(HIST("hHasTOFPi"), negdau.hasTOF(), negdau.pt());
          }
          QAHistos.fill(HIST("hRapXi"), casc.yXi());
        }

        // Count number of Xi candidates
        xicounter++;

        // Plot for estimates
        for (auto track : tracks) { // start loop over tracks
          if (isTrackFilter && !mTrackSelector.IsSelected(track)) {
            continue;
          }
          triggcounterForEstimates++;
          if (triggcounterForEstimates > 0)
            break;
        }
        if (triggcounterForEstimates && (TMath::Abs(casc.mXi() - o2::analysis::pdg::MassXiMinus) < 0.01))
          hhXiPairsvsPt->Fill(casc.pt()); // Fill the histogram with all the Xis produced in events with a trigger particle
        // End plot for estimates
      }
      if (isXiYN) {
        // Xis for YN interactions
        xicounterYN++;
        QAHistosTopologicalVariables.fill(HIST("hCascRadiusXiYN"), casc.cascradius());
      }
      if (isOmega) {
        QAHistos.fill(HIST("hMassOmegaAfterSelvsPt"), casc.mOmega(), casc.pt());
        QAHistos.fill(HIST("hPtOmega"), casc.pt());
        QAHistos.fill(HIST("hEtaOmega"), casc.eta());
        QAHistosTopologicalVariables.fill(HIST("hProperLifetimeOmega"), omegaproperlifetime);
        QAHistosTopologicalVariables.fill(HIST("hCascCosPAOmega"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()));
        QAHistosTopologicalVariables.fill(HIST("hV0CosPAOmega"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()));
        QAHistosTopologicalVariables.fill(HIST("hCascRadiusOmega"), casc.cascradius());
        QAHistosTopologicalVariables.fill(HIST("hV0RadiusOmega"), casc.v0radius());
        QAHistosTopologicalVariables.fill(HIST("hDCAV0ToPVOmega"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
        QAHistosTopologicalVariables.fill(HIST("hDCAV0DaughtersOmega"), casc.dcaV0daughters());
        QAHistosTopologicalVariables.fill(HIST("hDCACascDaughtersOmega"), casc.dcacascdaughters());
        QAHistosTopologicalVariables.fill(HIST("hDCABachToPVOmega"), TMath::Abs(casc.dcabachtopv()));
        QAHistosTopologicalVariables.fill(HIST("hDCAPosToPVOmega"), TMath::Abs(casc.dcapostopv()));
        QAHistosTopologicalVariables.fill(HIST("hDCANegToPVOmega"), TMath::Abs(casc.dcanegtopv()));
        QAHistosTopologicalVariables.fill(HIST("hInvMassLambdaOmega"), casc.mLambda());

        if (doextraQA) {

          // QA PID
          if (casc.sign() > 0) {
            QAHistos.fill(HIST("hTPCNsigmaOmegaBachKaPlus"), bachelor.tpcNSigmaKa(), bachelor.tpcInnerParam());
            QAHistos.fill(HIST("hTPCNsigmaOmegaV0PiPlus"), posdau.tpcNSigmaPi(), posdau.tpcInnerParam());
            QAHistos.fill(HIST("hTPCNsigmaOmegaV0AntiProton"), negdau.tpcNSigmaPr(), negdau.tpcInnerParam());
            QAHistos.fill(HIST("hHasTOFPi"), posdau.hasTOF(), posdau.pt());
            QAHistos.fill(HIST("hHasTOFPr"), negdau.hasTOF(), negdau.pt());
          } else {
            QAHistos.fill(HIST("hTPCNsigmaOmegaBachKaMinus"), bachelor.tpcNSigmaKa(), bachelor.tpcInnerParam());
            QAHistos.fill(HIST("hTPCNsigmaOmegaV0Proton"), posdau.tpcNSigmaPr(), posdau.tpcInnerParam());
            QAHistos.fill(HIST("hTPCNsigmaOmegaV0PiMinus"), negdau.tpcNSigmaPi(), negdau.tpcInnerParam());
            QAHistos.fill(HIST("hHasTOFPr"), posdau.hasTOF(), posdau.pt());
            QAHistos.fill(HIST("hHasTOFPi"), negdau.hasTOF(), negdau.pt());
          }
          QAHistos.fill(HIST("hHasTOFBachKa"), bachelor.hasTOF(), bachelor.pt());
          QAHistos.fill(HIST("hRapOmega"), casc.yOmega());
        }

        // Count number of Omega candidates
        omegacounter++;
      }
      if (isOmegalargeR) {
        omegalargeRcounter++;
        QAHistosTopologicalVariables.fill(HIST("hCascRadiusOmegaLargeR"), casc.cascradius());
      }
    } // end loop over cascades

    // Omega trigger definition
    if (omegacounter > 0) {
      keepEvent[0] = true;
    }

    bool EvtwhMinPt[11];
    bool EvtwhMinPtXi[11];
    float ThrdPt[11];
    for (int i = 0; i < 11; i++) {
      EvtwhMinPt[i] = 0.;
      EvtwhMinPtXi[i] = 0.;
      ThrdPt[i] = static_cast<float>(i);
    }

    // QA tracks
    for (auto track : tracks) { // start loop over tracks
      if (isTrackFilter && !mTrackSelector.IsSelected(track)) {
        continue;
      }
      triggcounterAllEv++;
      QAHistosTriggerParticles.fill(HIST("hPtTriggerAllEv"), track.pt());
      QAHistosTriggerParticles.fill(HIST("hPhiTriggerAllEv"), track.phi(), track.pt());
      QAHistosTriggerParticles.fill(HIST("hEtaTriggerAllEv"), track.eta(), track.pt());
      QAHistosTriggerParticles.fill(HIST("hDCAxyTriggerAllEv"), track.dcaXY(), track.pt());
      QAHistosTriggerParticles.fill(HIST("hDCAzTriggerAllEv"), track.dcaZ(), track.pt());
      for (int i = 0; i < 11; i++) {
        if (track.pt() > ThrdPt[i])
          EvtwhMinPt[i] = 1;
      }
    } // end loop over tracks
    for (int i = 0; i < 11; i++) {
      if (EvtwhMinPt[i]) {
        EventsvsMultiplicity.fill(HIST("hadEventsvsMultiplicityZeqV0AvsPt"), collision.multZeqFV0A(), i + 0.5);
        EventsvsMultiplicity.fill(HIST("hadEventsvsMultiplicityZeqT0MvsPt"), collision.multZeqFT0A() + collision.multZeqFT0C(), i + 0.5);
        EventsvsMultiplicity.fill(HIST("hadEventsvsMultiplicityZeqNTracksPVvsPt"), collision.multZeqNTracksPV(), i + 0.5);
      }
    }
    if (triggcounterAllEv > 0) {
      hProcessedEvents->Fill(1.5);
      if (doextraQA) {
        EventsvsMultiplicity.fill(HIST("hadEventsvsMultiplicityZeqV0A"), collision.multZeqFV0A());
        EventsvsMultiplicity.fill(HIST("hadEventsvsMultiplicityZeqT0M"), collision.multZeqFT0A() + collision.multZeqFT0C());
        EventsvsMultiplicity.fill(HIST("hadEventsvsMultiplicityZeqNTracksPV"), collision.multZeqNTracksPV());
      }
    }
    QAHistosTriggerParticles.fill(HIST("hTriggeredParticlesAllEv"), triggcounterAllEv);

    // High-pT hadron + Xi trigger definition
    if (xicounter > 0) {
      for (auto track : tracks) { // start loop over tracks
        if (isTrackFilter && !mTrackSelector.IsSelected(track)) {
          continue;
        }
        triggcounter++;
        QAHistosTriggerParticles.fill(HIST("hPtTriggerSelEv"), track.pt());
        for (int i = 0; i < 11; i++) {
          if (track.pt() > ThrdPt[i])
            EvtwhMinPtXi[i] = 1;
        }
        keepEvent[1] = true;
      } // end loop over tracks
      QAHistosTriggerParticles.fill(HIST("hTriggeredParticlesSelEv"), triggcounter);
    }

    for (int i = 0; i < 11; i++) {
      if (EvtwhMinPtXi[i])
        hEvtvshMinPt->Fill(i + 0.5);
    }

    // 2Xi trigger definition
    if (xicounter > 1) {
      keepEvent[2] = true;
    }

    // 3Xi trigger definition
    if (xicounter > 2) {
      keepEvent[3] = true;
    }

    // 4Xi trigger definition
    if (xicounter > 3) {
      keepEvent[4] = true;
    }

    // Single-Xi (YN) trigger definition
    if (xicounterYN > 0) {
      keepEvent[5] = true;
    }

    // Omega with high radius trigger definition
    if (omegalargeRcounter > 0) {
      keepEvent[6] = true;
    }

    // strangeness tracking selection
    const auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      runNumber = bc.runNumber();
      auto timestamp = bc.timestamp();

      if (o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, timestamp)) {
        o2::base::Propagator::initFieldFromGRP(grpo);
      } else if (o2::parameters::GRPMagField* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpMagPath, timestamp)) {
        o2::base::Propagator::initFieldFromGRP(grpmag);
      } else {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpMagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << timestamp;
      }
    }

    const auto primaryVertex = getPrimaryVertex(collision);
    o2::dataformats::DCA impactParameterTrk;

    for (const auto& casc : fullCasc) {
      QAHistosStrangenessTracking.fill(HIST("hPtCascCand"), casc.pt());
    }

    for (const auto& trackedCascade : trackedCascades) {
      const auto trackCasc = trackedCascade.track_as<DaughterTracks>();
      QAHistosStrangenessTracking.fill(HIST("hPtCascTracked"), trackCasc.pt());
      QAHistosStrangenessTracking.fill(HIST("hStRVsPtTrkCasc"), trackCasc.pt(), RecoDecay::sqrtSumOfSquares(trackCasc.x(), trackCasc.y()));
      // QAHistosStrangenessTracking.fill(HIST("hMassOmegaTrkCasc"), trackedCascade.omegaMass());
      // QAHistosStrangenessTracking.fill(HIST("hMassXiTrkCasc"), trackedCascade.xiMass());
      QAHistosStrangenessTracking.fill(HIST("hMatchChi2TrkCasc"), trackedCascade.matchingChi2());
      QAHistosStrangenessTracking.fill(HIST("hMassOmegaVsMatchChi2TrkCasc"), trackedCascade.omegaMass(), trackedCascade.matchingChi2());
      QAHistosStrangenessTracking.fill(HIST("hMassXiVsMatchChi2TrkCasc"), trackedCascade.xiMass(), trackedCascade.matchingChi2());
      QAHistosStrangenessTracking.fill(HIST("hMassOmegaVsTopChi2TrkCasc"), trackedCascade.omegaMass(), trackedCascade.topologyChi2());
      QAHistosStrangenessTracking.fill(HIST("hMassXiVsTopChi2TrkCasc"), trackedCascade.xiMass(), trackedCascade.topologyChi2());

      auto trackParCovTrk = getTrackParCov(trackCasc);
      o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackParCovTrk, bz, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrNONE, &impactParameterTrk);

      QAHistosStrangenessTracking.fill(HIST("hDcaXY"), impactParameterTrk.getY());
      QAHistosStrangenessTracking.fill(HIST("hDcaXYVsPt"), trackParCovTrk.getPt(), impactParameterTrk.getY());
      QAHistosStrangenessTracking.fill(HIST("hDcaZ"), impactParameterTrk.getZ());
      QAHistosStrangenessTracking.fill(HIST("hDcaZVsPt"), trackParCovTrk.getPt(), impactParameterTrk.getZ());
      QAHistosStrangenessTracking.fill(HIST("hDcaVsPt"), impactParameterTrk.getY(), trackCasc.pt());
      QAHistosStrangenessTracking.fill(HIST("hDcaVsR"), impactParameterTrk.getY(), RecoDecay::sqrtSumOfSquares(trackCasc.x(), trackCasc.y()));

      // const auto itsTrack = trackedCascade.itsTrack();
      const auto cascade = trackedCascade.cascade();
      const auto bachelor = cascade.bachelor_as<DaughterTracks>();
      const auto v0 = cascade.v0_as<o2::aod::V0sLinked>();
      const auto negTrack = v0.negTrack_as<DaughterTracks>();
      const auto posTrack = v0.posTrack_as<DaughterTracks>();

      std::array<double, 2> masses{o2::analysis::pdg::MassProton, o2::analysis::pdg::MassPiMinus};
      std::array<std::array<float, 3>, 2> momenta;
      std::array<double, 2> nsigma;
      if (trackCasc.sign() < 0) {
        // Omega-, Xi-
        momenta[0] = {posTrack.px(), posTrack.py(), posTrack.pz()};
        momenta[1] = {negTrack.px(), negTrack.py(), negTrack.pz()};
        nsigma[0] = posTrack.tpcNSigmaPr();
        nsigma[1] = negTrack.tpcNSigmaPi();
      } else {
        // Omega+, Xi+
        momenta[0] = {negTrack.px(), negTrack.py(), negTrack.pz()};
        momenta[1] = {posTrack.px(), posTrack.py(), posTrack.pz()};
        nsigma[0] = negTrack.tpcNSigmaPr();
        nsigma[1] = posTrack.tpcNSigmaPi();
      }

      const auto v0mass = RecoDecay::m(momenta, masses);
      QAHistosStrangenessTracking.fill(HIST("hMassV0TrkCasc"), v0mass);
      QAHistosStrangenessTracking.fill(HIST("hNSigmaTpcPrTrkCascV0"), nsigma[0]);
      QAHistosStrangenessTracking.fill(HIST("hNSigmaTpcPiTrkCascV0"), nsigma[1]);
      QAHistosStrangenessTracking.fill(HIST("hNSigmaTpcPiTrkCascBachelor"), bachelor.tpcNSigmaPi());
      QAHistosStrangenessTracking.fill(HIST("hNSigmaTpcKaTrkCascBachelor"), bachelor.tpcNSigmaKa());

      momenta[0] = {posTrack.px() + negTrack.px(), posTrack.py() + negTrack.py(), posTrack.pz() + negTrack.pz()};
      momenta[1] = {bachelor.px(), bachelor.py(), bachelor.pz()};
      masses = {o2::analysis::pdg::MassLambda0, o2::analysis::pdg::MassK0};
      const auto massOmega = RecoDecay::m(momenta, masses);
      if (posTrack.hasTPC() && negTrack.hasTPC()) {
        QAHistosStrangenessTracking.fill(HIST("hMassOmegaTrkCasc"), massOmega);
      }
      masses = {o2::analysis::pdg::MassLambda0, o2::analysis::pdg::MassPi0};
      const auto massXi = RecoDecay::m(momenta, masses);
      if (posTrack.hasTPC() && negTrack.hasTPC()) {
        QAHistosStrangenessTracking.fill(HIST("hMassXiTrkCasc"), massXi);
      }

      if ((trackCasc.pt() > minPtTrackedCascade) &&
          (trackedCascade.matchingChi2() < maxMatchingChi2TrackedCascade) &&
          (std::abs(v0mass - o2::analysis::pdg::MassLambda0) < massWindowLambda) &&
          (std::abs(nsigma[0]) < maxNSigmaV0PrTrackedCascade) &&
          (std::abs(nsigma[1]) < maxNSigmaV0PiTrackedCascade)) {
        // Xi
        if ((std::abs(massXi - o2::analysis::pdg::MassXiMinus) < massWindowTrackedXi) &&
            (std::abs(bachelor.tpcNSigmaPi()) < maxNSigmaBachelorTrackedXi)) {
          keepEvent[7] = true;
        }
        // Omega
        if ((std::abs(massOmega - o2::analysis::pdg::MassOmegaMinus) < massWindowTrackedOmega) &&
            (std::abs(massXi - o2::analysis::pdg::MassXiMinus) >= massWindowXiExclTrackedOmega) &&
            (std::abs(bachelor.tpcNSigmaKa()) < maxNSigmaBachelorTrackedOmega)) {
          keepEvent[8] = true;
        }
      }
    }

    for (const auto& trackedV0 : trackedV0s) {
      const auto trackV0 = trackedV0.track_as<DaughterTracks>();
      QAHistosStrangenessTracking.fill(HIST("hMassH3LTrkV0"), trackedV0.h3Lmass());
      QAHistosStrangenessTracking.fill(HIST("hMassH4LTrkV0"), trackedV0.h4Lmass());
      if (trackV0.pt() > minPtTrackedV0) {
        keepEvent[9] = true;
      }
    }

    for (const auto& tracked3Body : tracked3Bodys) {
      const auto track3Body = tracked3Body.track_as<DaughterTracks>();
      QAHistosStrangenessTracking.fill(HIST("hMassH3LTrk3body"), tracked3Body.h3Lmass());
      QAHistosStrangenessTracking.fill(HIST("hMassHe4LTrk3body"), tracked3Body.he4Lmass());
      if (track3Body.pt() > minPtTracked3Body) {
        keepEvent[10] = true;
      }
    }

    // Fill centrality dependent histos
    if (keepEvent[0]) {
      hProcessedEvents->Fill(2.5);
    }
    if (keepEvent[1]) {
      hProcessedEvents->Fill(3.5);
    }
    if (keepEvent[2]) {
      hProcessedEvents->Fill(4.5);
    }
    if (keepEvent[3]) {
      hProcessedEvents->Fill(5.5);
    }
    if (keepEvent[4]) {
      hProcessedEvents->Fill(6.5);
    }
    if (keepEvent[5]) {
      hProcessedEvents->Fill(7.5);
    }
    if (keepEvent[6]) {
      hProcessedEvents->Fill(8.5);
    }
    if (xicounter > 0) {
      hProcessedEvents->Fill(9.5);
    }
    if (keepEvent[7]) {
      hProcessedEvents->Fill(10.5);
    }
    if (keepEvent[8]) {
      hProcessedEvents->Fill(11.5);
    }
    if (keepEvent[9]) {
      hProcessedEvents->Fill(12.5);
    }
    if (keepEvent[10]) {
      hProcessedEvents->Fill(13.5);
    }

    // Filling the table
    fillTriggerTable(keepEvent);
  }
  //
  PROCESS_SWITCH(strangenessFilter, processRun3, "Process Run3", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangenessFilter>(cfgc, TaskName{"lf-strangeness-filter"})};
}
