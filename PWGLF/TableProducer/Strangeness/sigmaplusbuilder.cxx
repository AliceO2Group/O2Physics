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

/// \file   sigmaplusbuilder.cxx
/// \brief  Task for Sigma+ -> p + pi0 reconstruction, pi0 reconstructed via one converted photon (PCM)
/// \author Henrik Fribert (TUM)

#include "PWGLF/DataModel/LFKinkDecayTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/PID.h>
#include <ReconstructionDataFormats/Track.h>

#include <TF1.h>
#include <TMCProcess.h>
#include <TPDGCode.h>
#include <TVector3.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA,
                             aod::pidTPCFullEl, aod::pidTPCFullPr, aod::pidTOFFullPr>;
using TracksFullMC = soa::Join<TracksFull, aod::McTrackLabels>;
using CollisionsFull = aod::Collisions;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels>;

struct Sigmaplusbuilder {

  // photon (PCM) selection
  Configurable<float> photonMaxMass{"photonMaxMass", 0.20, "Max photon mass (GeV/c^2)"};
  Configurable<float> photonMinRapidity{"photonMinRapidity", -0.8, "Min photon rapidity"};
  Configurable<float> photonMaxRapidity{"photonMaxRapidity", 0.8, "Max photon rapidity"};
  Configurable<float> photonDauEtaMin{"photonDauEtaMin", -0.8, "Min eta of photon daughter tracks"};
  Configurable<float> photonDauEtaMax{"photonDauEtaMax", 0.8, "Max eta of photon daughter tracks"};
  Configurable<float> photonMinRadius{"photonMinRadius", 3.0, "Min photon conversion radius (cm)"};
  Configurable<float> photonMaxRadius{"photonMaxRadius", 115., "Max photon conversion radius (cm)"};
  Configurable<float> photonMinV0cospa{"photonMinV0cospa", 0.80, "Min V0 CosPA"};
  Configurable<float> photonMaxDCAV0Dau{"photonMaxDCAV0Dau", 3.5, "Max DCA between photon daughters (cm)"};
  Configurable<float> photonMaxQt{"photonMaxQt", 0.15, "Max Armenteros qT for photons (GeV/c)"};
  Configurable<float> photonMaxAlpha{"photonMaxAlpha", 1.0, "Max |Armenteros alpha| for photons"};
  Configurable<float> photonMaxTPCNSigmaEl{"photonMaxTPCNSigmaEl", 15, "Max |TPC nSigma_el| for photon daughters"};

  // proton selection
  Configurable<float> protonMinPt{"protonMinPt", 0.3, "Minimum proton pT (GeV/c)"};
  Configurable<float> protonMaxEta{"protonMaxEta", 0.9, "Maximum |eta| for proton track"};
  Configurable<float> protonMaxTPCNSigma{"protonMaxTPCNSigma", 4, "Max |TPC nSigma_pr| for proton"};
  Configurable<float> protonMaxTOFNSigma{"protonMaxTOFNSigma", 4, "Max |TOF nSigma_pr| for proton, if TOF available"};
  Configurable<float> protonPtMinRequireTOF{"protonPtMinRequireTOF", 0.75, "Above this pT, require TOF PID for proton"};

  // proton-photon candidate selection
  Configurable<float> candMaxDcaProtonGamma{"candMaxDcaProtonGamma", 5.0, "Max DCA between proton and photon at the fitted vertex (cm)"};
  Configurable<float> candMinRadius{"candMinRadius", 1.0, "Min candidate decay radius (cm)"};
  Configurable<float> candMaxRadius{"candMaxRadius", 100., "Max candidate decay radius (cm)"};

  // missing-photon discriminant retry
  // resolution-shaped function and tuned parameters from Run 2 used
  Configurable<int> discrRetryMaxIter{"discrRetryMaxIter", 10, "Max attempts to recover a negative discriminant by perturbing the flight direction"};
  Configurable<float> discrRetryThetaRange{"discrRetryThetaRange", 0.02, "Max |delta-theta| perturbation of the flight direction per retry (rad)"};
  Configurable<float> discrRetryThetaPar0{"discrRetryThetaPar0", 0.000133299, "par0 of the delta-theta resolution function"};
  Configurable<float> discrRetryThetaPar1{"discrRetryThetaPar1", 0.0016761, "par1 of the delta-theta resolution function"};
  Configurable<float> discrRetryPhiRange{"discrRetryPhiRange", 0.02, "Max |delta-phi| perturbation of the flight direction per retry (rad)"};
  Configurable<float> discrRetryPhiPar0{"discrRetryPhiPar0", 0.000129845, "par0 of the delta-phi resolution function"};
  Configurable<float> discrRetryPhiPar1{"discrRetryPhiPar1", 0.00199688, "par1 of the delta-phi resolution function"};

  Configurable<std::string> ccdbPath{"ccdbPath", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  Produces<aod::SigmaPlusCands> sigmaPlusCands;
  Produces<aod::SigmaPlusCandsMC> sigmaPlusCandsMC;

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::vertexing::DCAFitterN<2> fitter;
  int mRunNumber = 0;
  float mBz = 0;
  TF1 mThetaResoFunc;
  TF1 mPhiResoFunc;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Preslice<aod::V0Datas> v0PerCollision = aod::v0data::collisionId;
  Preslice<TracksFull> tracksPerCollision = aod::track::collisionId;
  Preslice<TracksFullMC> tracksPerCollisionMC = aod::track::collisionId;

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbPath);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);

    mThetaResoFunc = TF1("thetaResoFunc", "[0]/(abs(x)+[1])", -discrRetryThetaRange, discrRetryThetaRange);
    mThetaResoFunc.SetParameters(discrRetryThetaPar0.value, discrRetryThetaPar1.value);
    mPhiResoFunc = TF1("phiResoFunc", "[0]/(abs(x)+[1])", -discrRetryPhiRange, discrRetryPhiRange);
    mPhiResoFunc.SetParameters(discrRetryPhiPar0.value, discrRetryPhiPar1.value);

    const AxisSpec axisVertexZ{100, -15., 15., "vrtx_{Z} (cm)"};

    const AxisSpec axisPhotonSel{11, -0.5, 10.5, "selection step"};
    const AxisSpec axisPhotonMass{200, 0., 0.3, "m_{#gamma} (GeV/c^{2})"};
    const AxisSpec axisPhotonPt{100, 0., 5., "#it{p}_{T,#gamma} (GeV/c)"};
    const AxisSpec axisPhotonRadius{200, 0., 200., "R_{conv} (cm)"};
    const AxisSpec axisAlpha{100, -1., 1., "#alpha_{AP}"};
    const AxisSpec axisQt{100, 0., 0.3, "q_{T,AP} (GeV/c)"};
    const AxisSpec axisConvXY{200, -100., 100., "conv. point (cm)"};

    const AxisSpec axisProtonSel{5, -0.5, 4.5, "selection step"};
    const AxisSpec axisProtonPt{100, 0., 5., "#it{p}_{T,p} (GeV/c)"};
    const AxisSpec axisNSigma{100, -5., 5., "n#sigma"};

    histos.add("hVertexZ", "hVertexZ", kTH1F, {axisVertexZ});

    histos.add("Photon/hSelectionCounter", "Photon/hSelectionCounter", kTH1F, {axisPhotonSel});
    auto hPhotonSel = histos.get<TH1>(HIST("Photon/hSelectionCounter"));
    hPhotonSel->GetXaxis()->SetBinLabel(1, "All");
    hPhotonSel->GetXaxis()->SetBinLabel(2, "Mass");
    hPhotonSel->GetXaxis()->SetBinLabel(3, "Rapidity");
    hPhotonSel->GetXaxis()->SetBinLabel(4, "Neg eta");
    hPhotonSel->GetXaxis()->SetBinLabel(5, "Pos eta");
    hPhotonSel->GetXaxis()->SetBinLabel(6, "DCA daughters");
    hPhotonSel->GetXaxis()->SetBinLabel(7, "Radius");
    hPhotonSel->GetXaxis()->SetBinLabel(8, "CosPA");
    hPhotonSel->GetXaxis()->SetBinLabel(9, "Qt");
    hPhotonSel->GetXaxis()->SetBinLabel(10, "Alpha");
    hPhotonSel->GetXaxis()->SetBinLabel(11, "TPC nSigma el");

    histos.add("Photon/hMass", "Photon/hMass", kTH1F, {axisPhotonMass});
    histos.add("Photon/hPt", "Photon/hPt", kTH1F, {axisPhotonPt});
    histos.add("Photon/hRadius", "Photon/hRadius", kTH1F, {axisPhotonRadius});
    histos.add("Photon/h2ArmenterosPodolanski", "Photon/h2ArmenterosPodolanski", kTH2F, {axisAlpha, axisQt});
    histos.add("Photon/h2ConvPointXY", "Photon/h2ConvPointXY", kTH2F, {axisConvXY, axisConvXY});

    histos.add("Proton/hSelectionCounter", "Proton/hSelectionCounter", kTH1F, {axisProtonSel});
    auto hProtonSel = histos.get<TH1>(HIST("Proton/hSelectionCounter"));
    hProtonSel->GetXaxis()->SetBinLabel(1, "All");
    hProtonSel->GetXaxis()->SetBinLabel(2, "Pt");
    hProtonSel->GetXaxis()->SetBinLabel(3, "Eta");
    hProtonSel->GetXaxis()->SetBinLabel(4, "TPC nSigma");
    hProtonSel->GetXaxis()->SetBinLabel(5, "TOF nSigma");

    histos.add("Proton/hPt", "Proton/hPt", kTH1F, {axisProtonPt});
    histos.add("Proton/h2TPCNSigmaVsPt", "Proton/h2TPCNSigmaVsPt", kTH2F, {axisProtonPt, axisNSigma});
    histos.add("Proton/h2TOFNSigmaVsPt", "Proton/h2TOFNSigmaVsPt", kTH2F, {axisProtonPt, axisNSigma});

    const AxisSpec axisCandSel{7, -0.5, 6.5, "selection step"};
    const AxisSpec axisDca{100, 0., 10., "DCA(p,#gamma) (cm)"};
    const AxisSpec axisCandRadius{200, 0., 200., "R_{dec} (cm)"};
    const AxisSpec axisDiscriminant{200, -1., 1., "discriminant (GeV^{4}/#it{c}^{4})"};
    const AxisSpec axisMassSigma{200, 1.0, 1.4, "m_{p#gamma#gamma} (GeV/#it{c}^{2})"};
    const AxisSpec axisSigmaPt{100, 0., 6., "#it{p}_{T,#Sigma^{+}} (GeV/#it{c})"};
    const AxisSpec axisMomentum{100, 0., 10., "#it{p} (GeV/#it{c})"};

    histos.add("Candidate/hSelectionCounter", "Candidate/hSelectionCounter", kTH1F, {axisCandSel});
    auto hCandSel = histos.get<TH1>(HIST("Candidate/hSelectionCounter"));
    hCandSel->GetXaxis()->SetBinLabel(1, "All pairs");
    hCandSel->GetXaxis()->SetBinLabel(2, "Vertex fit");
    hCandSel->GetXaxis()->SetBinLabel(3, "DCA(p,#gamma)");
    hCandSel->GetXaxis()->SetBinLabel(4, "Radius");
    hCandSel->GetXaxis()->SetBinLabel(5, "Real root");
    hCandSel->GetXaxis()->SetBinLabel(6, "Valid root");
    hCandSel->GetXaxis()->SetBinLabel(7, "Filled");

    histos.add("Candidate/hDcaProtonGamma", "Candidate/hDcaProtonGamma", kTH1F, {axisDca});
    histos.add("Candidate/hRadius", "Candidate/hRadius", kTH1F, {axisCandRadius});
    histos.add("Candidate/hDiscriminant", "Candidate/hDiscriminant", kTH1F, {axisDiscriminant});
    histos.add("Candidate/hMassSigmaPlus", "Candidate/hMassSigmaPlus", kTH1F, {axisMassSigma});
    histos.add("Candidate/h2MassVsPt", "Candidate/h2MassVsPt", kTH2F, {axisSigmaPt, axisMassSigma});

    const AxisSpec axisDiscrIter{discrRetryMaxIter + 2, -0.5, discrRetryMaxIter + 1.5, "discriminant retry iteration"};
    histos.add("Candidate/hDiscriminantRetryIter", "Candidate/hDiscriminantRetryIter", kTH1F, {axisDiscrIter});

    if (doprocessMc) {
      histos.add("Photon/True/hSelectionCounter", "Photon/True/hSelectionCounter", kTH1F, {axisPhotonSel});
      auto hPhotonSelSignal = histos.get<TH1>(HIST("Photon/True/hSelectionCounter"));
      hPhotonSelSignal->GetXaxis()->SetBinLabel(1, "All");
      hPhotonSelSignal->GetXaxis()->SetBinLabel(2, "Mass");
      hPhotonSelSignal->GetXaxis()->SetBinLabel(3, "Rapidity");
      hPhotonSelSignal->GetXaxis()->SetBinLabel(4, "Neg eta");
      hPhotonSelSignal->GetXaxis()->SetBinLabel(5, "Pos eta");
      hPhotonSelSignal->GetXaxis()->SetBinLabel(6, "DCA daughters");
      hPhotonSelSignal->GetXaxis()->SetBinLabel(7, "Radius");
      hPhotonSelSignal->GetXaxis()->SetBinLabel(8, "CosPA");
      hPhotonSelSignal->GetXaxis()->SetBinLabel(9, "Qt");
      hPhotonSelSignal->GetXaxis()->SetBinLabel(10, "Alpha");
      hPhotonSelSignal->GetXaxis()->SetBinLabel(11, "TPC nSigma el");

      histos.add("Photon/True/hMass", "Photon/True/hMass", kTH1F, {axisPhotonMass});
      histos.add("Photon/True/hPt", "Photon/True/hPt", kTH1F, {axisPhotonPt});
      histos.add("Photon/True/hRadius", "Photon/True/hRadius", kTH1F, {axisPhotonRadius});
      histos.add("Photon/True/h2ArmenterosPodolanski", "Photon/True/h2ArmenterosPodolanski", kTH2F, {axisAlpha, axisQt});
      histos.add("Photon/True/h2ConvPointXY", "Photon/True/h2ConvPointXY", kTH2F, {axisConvXY, axisConvXY});

      histos.add("Proton/True/hPt", "Proton/True/hPt", kTH1F, {axisProtonPt});
      histos.add("Proton/True/h2TPCNSigmaVsPt", "Proton/True/h2TPCNSigmaVsPt", kTH2F, {axisProtonPt, axisNSigma});
      histos.add("Proton/True/h2TOFNSigmaVsPt", "Proton/True/h2TOFNSigmaVsPt", kTH2F, {axisProtonPt, axisNSigma});

      histos.add("Candidate/True/hSelectionCounter", "Candidate/True/hSelectionCounter", kTH1F, {axisCandSel});
      auto hCandSelSignal = histos.get<TH1>(HIST("Candidate/True/hSelectionCounter"));
      hCandSelSignal->GetXaxis()->SetBinLabel(1, "All pairs");
      hCandSelSignal->GetXaxis()->SetBinLabel(2, "Vertex fit");
      hCandSelSignal->GetXaxis()->SetBinLabel(3, "DCA(p,#gamma)");
      hCandSelSignal->GetXaxis()->SetBinLabel(4, "Radius");
      hCandSelSignal->GetXaxis()->SetBinLabel(5, "Real root");
      hCandSelSignal->GetXaxis()->SetBinLabel(6, "Valid root");
      hCandSelSignal->GetXaxis()->SetBinLabel(7, "Filled");

      histos.add("Candidate/True/hDcaProtonGamma", "Candidate/True/hDcaProtonGamma", kTH1F, {axisDca});
      histos.add("Candidate/True/hRadius", "Candidate/True/hRadius", kTH1F, {axisCandRadius});
      const AxisSpec axisVtxRes{200, 0., 20., "|vtx_{fit} - vtx_{MC}| (cm)"};
      histos.add("Candidate/True/hVertexResFromMcTruth", "Candidate/True/hVertexResFromMcTruth", kTH1F, {axisVtxRes});
      const AxisSpec axisMomRes{200, -1., 1., "(p_{fit} - p_{MC}) / p_{MC}"};
      histos.add("Candidate/True/hProtonMomResFromMcTruth", "Candidate/True/hProtonMomResFromMcTruth", kTH1F, {axisMomRes});
      histos.add("Candidate/True/hPhotonMomResFromMcTruth", "Candidate/True/hPhotonMomResFromMcTruth", kTH1F, {axisMomRes});
      const AxisSpec axisProtonFlightAngle{180, 0., 3.15, "angle(p_{proton}, n) (rad)"};
      histos.add("Candidate/True/hProtonFlightAngle", "Candidate/True/hProtonFlightAngle", kTH1F, {axisProtonFlightAngle});
      histos.add("Candidate/True/hDiscriminant", "Candidate/True/hDiscriminant", kTH1F, {axisDiscriminant});
      histos.add("Candidate/True/hDiscriminantRetryIter", "Candidate/True/hDiscriminantRetryIter", kTH1F, {axisDiscrIter});
      histos.add("Candidate/True/hMassSigmaPlus", "Candidate/True/hMassSigmaPlus", kTH1F, {axisMassSigma});
      histos.add("Candidate/True/h2MassVsPt", "Candidate/True/h2MassVsPt", kTH2F, {axisSigmaPt, axisMassSigma});

      histos.add("MC/hGenSigmaPlusPt", "MC/hGenSigmaPlusPt", kTH1F, {axisSigmaPt});

      const AxisSpec axisProtonTruthQA{3, -0.5, 2.5, "step"};
      histos.add("MC/hProtonTruthQA", "MC/hProtonTruthQA", kTH1F, {axisProtonTruthQA});
      auto hProtonTruthQA = histos.get<TH1>(HIST("MC/hProtonTruthQA"));
      hProtonTruthQA->GetXaxis()->SetBinLabel(1, "Accepted");
      hProtonTruthQA->GetXaxis()->SetBinLabel(2, "Has mcParticle");
      hProtonTruthQA->GetXaxis()->SetBinLabel(3, "Sigma+ mother");

      const AxisSpec axisPhotonTruthQA{5, -0.5, 4.5, "step"};
      histos.add("MC/hPhotonTruthQA", "MC/hPhotonTruthQA", kTH1F, {axisPhotonTruthQA});
      auto hPhotonTruthQA = histos.get<TH1>(HIST("MC/hPhotonTruthQA"));
      hPhotonTruthQA->GetXaxis()->SetBinLabel(1, "Accepted");
      hPhotonTruthQA->GetXaxis()->SetBinLabel(2, "Both legs have mcParticle");
      hPhotonTruthQA->GetXaxis()->SetBinLabel(3, "Shared real gamma mother");
      hPhotonTruthQA->GetXaxis()->SetBinLabel(4, "Gamma's mother is pi0");
      hPhotonTruthQA->GetXaxis()->SetBinLabel(5, "Pi0's mother is Sigma+");
    }

    if (doprocessFindable) {
      const AxisSpec axisDetectorPresence{3, -0.5, 2.5, "detector"};
      const AxisSpec axisDuplicateTrack{2, -0.5, 1.5, "track"};
      const AxisSpec axisV0Presence{2, -0.5, 1.5, "V0 match"};
      const AxisSpec axisTPCClusters{160, -0.5, 159.5, "TPC clusters"};
      const AxisSpec axisPhotonMomResolution{200, -1., 1., "(#it{p}_{reco,#gamma} - #it{p}_{MC,#gamma})/#it{p}_{MC,#gamma}"};
      const AxisSpec axisPhotonCosPA{200, -1., 1., "cosPA_{#gamma}"};
      const AxisSpec axisPairMass{200, 0., 0.1, "m_{e^{+}e^{-}} (GeV/#it{c}^{2})"};
      histos.add("Findable/hSigmaPlusPt", "Findable/hSigmaPlusPt", kTH1F, {axisSigmaPt});
      histos.add("Findable/h2ProtonPtVsSigmaPlusPt", "Findable/h2ProtonPtVsSigmaPlusPt", kTH2F, {axisSigmaPt, axisMomentum});
      histos.add("Findable/hElectronPt", "Findable/hElectronPt", kTH1F, {axisMomentum});
      histos.add("Findable/hPositronPt", "Findable/hPositronPt", kTH1F, {axisMomentum});
      histos.add("Findable/hElectronPositronMass", "Findable/hElectronPositronMass", kTH1F, {axisPairMass});
      histos.add("Findable/hPhotonConversionRadius", "Findable/hPhotonConversionRadius", kTH1F, {axisPhotonRadius});
      histos.add("Findable/hPhotonMomentumResolution", "Findable/hPhotonMomentumResolution", kTH1F, {axisPhotonMomResolution});
      histos.add("Findable/hPhotonCosPA", "Findable/hPhotonCosPA", kTH1F, {axisPhotonCosPA});
      histos.add("Findable/hConversionPairV0Presence", "Findable/hConversionPairV0Presence", kTH1F, {axisV0Presence});
      histos.add("Findable/hDuplicateConversionTrackCounter", "Findable/hDuplicateConversionTrackCounter", kTH1F, {axisDuplicateTrack});
      histos.add("Findable/hDuplicateElectronTPCNClsFound", "Findable/hDuplicateElectronTPCNClsFound", kTH1F, {axisTPCClusters});
      histos.add("Findable/hDuplicatePositronTPCNClsFound", "Findable/hDuplicatePositronTPCNClsFound", kTH1F, {axisTPCClusters});
      histos.add("Findable/hElectronDetectorPresence", "Findable/hElectronDetectorPresence", kTH1F, {axisDetectorPresence});
      histos.add("Findable/hPositronDetectorPresence", "Findable/hPositronDetectorPresence", kTH1F, {axisDetectorPresence});
      auto hConversionPairV0Presence = histos.get<TH1>(HIST("Findable/hConversionPairV0Presence"));
      hConversionPairV0Presence->GetXaxis()->SetBinLabel(1, "valid pair");
      hConversionPairV0Presence->GetXaxis()->SetBinLabel(2, "in V0");
      auto hDuplicateConversionTrackCounter = histos.get<TH1>(HIST("Findable/hDuplicateConversionTrackCounter"));
      hDuplicateConversionTrackCounter->GetXaxis()->SetBinLabel(1, "e^{-}");
      hDuplicateConversionTrackCounter->GetXaxis()->SetBinLabel(2, "e^{+}");
      auto hElectronDetectorPresence = histos.get<TH1>(HIST("Findable/hElectronDetectorPresence"));
      auto hPositronDetectorPresence = histos.get<TH1>(HIST("Findable/hPositronDetectorPresence"));
      hElectronDetectorPresence->GetXaxis()->SetBinLabel(1, "ITS");
      hElectronDetectorPresence->GetXaxis()->SetBinLabel(2, "TPC");
      hElectronDetectorPresence->GetXaxis()->SetBinLabel(3, "TOF");
      hPositronDetectorPresence->GetXaxis()->SetBinLabel(1, "ITS");
      hPositronDetectorPresence->GetXaxis()->SetBinLabel(2, "TPC");
      hPositronDetectorPresence->GetXaxis()->SetBinLabel(3, "TOF");
    }
  }

  // photon (PCM) candidate selection
  template <bool IsMC, typename TTracks, typename TV0>
  bool selectPhoton(const TV0& v0)
  {
    auto posTrack = v0.template posTrack_as<TTracks>();
    auto negTrack = v0.template negTrack_as<TTracks>();

    bool isSignal = false;
    if constexpr (IsMC) {
      if (posTrack.has_mcParticle() && negTrack.has_mcParticle()) {
        auto mcPos = posTrack.template mcParticle_as<aod::McParticles>();
        auto mcNeg = negTrack.template mcParticle_as<aod::McParticles>();
        isSignal = findSigmaPlusMotherOfPhoton(mcPos, mcNeg) >= 0;
      }
    }

    auto fillPhotonStep = [&](int step) {
      histos.fill(HIST("Photon/hSelectionCounter"), step);
      if constexpr (IsMC) {
        if (isSignal) {
          histos.fill(HIST("Photon/True/hSelectionCounter"), step);
        }
      }
    };
    fillPhotonStep(0);

    if (v0.mGamma() < 0 || v0.mGamma() > photonMaxMass) {
      return false;
    }
    fillPhotonStep(1);

    float photonY = RecoDecay::y(std::array{v0.px(), v0.py(), v0.pz()}, o2::constants::physics::MassGamma);
    if (photonY < photonMinRapidity || photonY > photonMaxRapidity) {
      return false;
    }
    fillPhotonStep(2);

    if (v0.negativeeta() < photonDauEtaMin || v0.negativeeta() > photonDauEtaMax) {
      return false;
    }
    fillPhotonStep(3);

    if (v0.positiveeta() < photonDauEtaMin || v0.positiveeta() > photonDauEtaMax) {
      return false;
    }
    fillPhotonStep(4);

    if (std::abs(v0.dcaV0daughters()) > photonMaxDCAV0Dau) {
      return false;
    }
    fillPhotonStep(5);

    if (v0.v0radius() < photonMinRadius || v0.v0radius() > photonMaxRadius) {
      return false;
    }
    fillPhotonStep(6);

    if (v0.v0cosPA() < photonMinV0cospa) {
      return false;
    }
    fillPhotonStep(7);

    if (v0.qtarm() > photonMaxQt) {
      return false;
    }
    fillPhotonStep(8);

    if (std::abs(v0.alpha()) > photonMaxAlpha) {
      return false;
    }
    fillPhotonStep(9);

    if (std::abs(posTrack.tpcNSigmaEl()) > photonMaxTPCNSigmaEl || std::abs(negTrack.tpcNSigmaEl()) > photonMaxTPCNSigmaEl) {
      return false;
    }
    fillPhotonStep(10);

    histos.fill(HIST("Photon/hMass"), v0.mGamma());
    histos.fill(HIST("Photon/hPt"), v0.pt());
    histos.fill(HIST("Photon/hRadius"), v0.v0radius());
    histos.fill(HIST("Photon/h2ArmenterosPodolanski"), v0.alpha(), v0.qtarm());
    histos.fill(HIST("Photon/h2ConvPointXY"), v0.x(), v0.y());

    return true;
  }

  // proton candidate selection
  template <typename TTrack>
  bool selectProton(const TTrack& track)
  {
    histos.fill(HIST("Proton/hSelectionCounter"), 0);

    if (track.pt() < protonMinPt) {
      return false;
    }
    histos.fill(HIST("Proton/hSelectionCounter"), 1);

    if (std::abs(track.eta()) > protonMaxEta) {
      return false;
    }
    histos.fill(HIST("Proton/hSelectionCounter"), 2);

    if (std::abs(track.tpcNSigmaPr()) > protonMaxTPCNSigma) {
      return false;
    }
    histos.fill(HIST("Proton/hSelectionCounter"), 3);

    if (track.pt() > protonPtMinRequireTOF) {
      if (!track.hasTOF() || std::abs(track.tofNSigmaPr()) > protonMaxTOFNSigma) {
        return false;
      }
    }
    histos.fill(HIST("Proton/hSelectionCounter"), 4);

    histos.fill(HIST("Proton/hPt"), track.pt());
    histos.fill(HIST("Proton/h2TPCNSigmaVsPt"), track.pt(), track.tpcNSigmaPr());
    if (track.hasTOF()) {
      histos.fill(HIST("Proton/h2TOFNSigmaVsPt"), track.pt(), track.tofNSigmaPr());
    }

    return true;
  }

  // small vector helpers used in the Sigma+ reconstruction
  static float dot3(const std::array<float, 3>& a, const std::array<float, 3>& b)
  {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  }
  static std::array<float, 3> cross3(const std::array<float, 3>& a, const std::array<float, 3>& b)
  {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
  }
  static std::array<float, 3> normalize3(const std::array<float, 3>& a)
  {
    float norm = std::sqrt(dot3(a, a));
    return {a[0] / norm, a[1] / norm, a[2] / norm};
  }

  // MC truth
  template <typename TMcPart>
  int findSigmaPlusMotherOfProton(const TMcPart& mcProton)
  {
    if (std::abs(mcProton.pdgCode()) != PDG_t::kProton) {
      return -1;
    }
    auto const& mothers = mcProton.template mothers_as<aod::McParticles>();
    if (mothers.empty() || std::abs(mothers.front().pdgCode()) != PDG_t::kSigmaPlus) {
      return -1;
    }
    return mothers.front().globalIndex();
  }

  template <typename TMcPart>
  int findSigmaPlusMotherOfPhoton(const TMcPart& mcPos, const TMcPart& mcNeg)
  {
    auto const& posMothers = mcPos.template mothers_as<aod::McParticles>();
    auto const& negMothers = mcNeg.template mothers_as<aod::McParticles>();
    if (posMothers.empty() || negMothers.empty()) {
      return -1;
    }
    auto mcGamma = posMothers.front();
    if (mcGamma.globalIndex() != negMothers.front().globalIndex() || mcGamma.pdgCode() != PDG_t::kGamma) {
      return -1;
    }

    auto const& pi0Mothers = mcGamma.template mothers_as<aod::McParticles>();
    if (pi0Mothers.empty() || std::abs(pi0Mothers.front().pdgCode()) != PDG_t::kPi0) {
      return -1;
    }

    auto const& sigmaMothers = pi0Mothers.front().template mothers_as<aod::McParticles>();
    if (sigmaMothers.empty() || std::abs(sigmaMothers.front().pdgCode()) != PDG_t::kSigmaPlus) {
      return -1;
    }
    return sigmaMothers.front().globalIndex();
  }

  template <typename TMcPart>
  bool isSigmaPlusToProtonPi0(const TMcPart& mcPart)
  {
    if (std::abs(mcPart.pdgCode()) != PDG_t::kSigmaPlus) {
      return false;
    }
    int pdgProton = mcPart.pdgCode() > 0 ? PDG_t::kProton : PDG_t::kProtonBar;
    bool hasProton = false;
    bool hasPi0 = false;
    for (const auto& daughter : mcPart.template daughters_as<aod::McParticles>()) {
      hasProton |= (daughter.pdgCode() == pdgProton);
      hasPi0 |= (std::abs(daughter.pdgCode()) == PDG_t::kPi0);
    }
    return hasProton && hasPi0;
  }

  template <bool IsElectron, typename TTrack>
  void fillFindableTrackDetectors(const TTrack& track)
  {
    auto fillDetector = [&](int detector) {
      if constexpr (IsElectron) {
        histos.fill(HIST("Findable/hElectronDetectorPresence"), detector);
      } else {
        histos.fill(HIST("Findable/hPositronDetectorPresence"), detector);
      }
    };

    if (track.hasITS()) {
      fillDetector(0);
    }
    if (track.hasTPC()) {
      fillDetector(1);
    }
    if (track.hasTOF()) {
      fillDetector(2);
    }
  }

  template <typename TMcPart>
  int findSigmaPlusMotherOfConversionElectron(const TMcPart& mcElectron, int& gammaIndex, std::array<float, 3>& conversionVertex)
  {
    if (std::abs(mcElectron.pdgCode()) != PDG_t::kElectron || mcElectron.getProcess() != TMCProcess::kPPair) {
      return -1;
    }

    auto const& gammaMothers = mcElectron.template mothers_as<aod::McParticles>();
    if (gammaMothers.empty() || gammaMothers.front().pdgCode() != PDG_t::kGamma) {
      return -1;
    }

    auto mcGamma = gammaMothers.front();
    if (mcGamma.getProcess() != TMCProcess::kPDecay) {
      return -1;
    }

    auto const& pi0Mothers = mcGamma.template mothers_as<aod::McParticles>();
    if (pi0Mothers.empty() || std::abs(pi0Mothers.front().pdgCode()) != PDG_t::kPi0) {
      return -1;
    }

    auto mcPi0 = pi0Mothers.front();
    auto const& sigmaMothers = mcPi0.template mothers_as<aod::McParticles>();
    if (sigmaMothers.empty() || std::abs(sigmaMothers.front().pdgCode()) != PDG_t::kSigmaPlus) {
      return -1;
    }

    gammaIndex = mcGamma.globalIndex();
    conversionVertex = {mcElectron.vx(), mcElectron.vy(), mcElectron.vz()};
    return sigmaMothers.front().globalIndex();
  }

  // Build a Sigma+ -> p pi0 candidate from a proton track and a PCM photon
  template <bool IsMC, typename TTracks, typename TTrack, typename TV0>
  void buildSigmaPlusCandidate(const TTrack& protonTrack, const TV0& photon, const std::array<float, 3>& pv)
  {
    auto posTrack = photon.template posTrack_as<TTracks>();
    auto negTrack = photon.template negTrack_as<TTracks>();

    bool isSignal = false;
    std::array<float, 3> mcTrueVtx{};       // Sigma+ decay vertex
    std::array<float, 3> mcTrueMomProton{}; // true MC proton momentum
    std::array<float, 3> mcTrueMomGamma{};  // true MC momentum of the measured photon
    int protonPdgCode = 0;
    int protonMotherPdgCode = 0;
    int gammaPdgCode = 0;
    int gammaMotherPdgCode = 0;
    int gammaGMotherPdgCode = 0;
    if constexpr (IsMC) {
      if (protonTrack.has_mcParticle()) {
        auto mcProton = protonTrack.template mcParticle_as<aod::McParticles>();
        protonPdgCode = mcProton.pdgCode();
        auto const& protonMothers = mcProton.template mothers_as<aod::McParticles>();
        if (!protonMothers.empty()) {
          protonMotherPdgCode = protonMothers.front().pdgCode();
        }

        int protonSigmaIdx = findSigmaPlusMotherOfProton(mcProton);

        if (posTrack.has_mcParticle() && negTrack.has_mcParticle()) {
          auto mcPos = posTrack.template mcParticle_as<aod::McParticles>();
          auto mcNeg = negTrack.template mcParticle_as<aod::McParticles>();

          auto const& posMothers = mcPos.template mothers_as<aod::McParticles>();
          if (!posMothers.empty()) {
            auto mcGamma = posMothers.front();
            gammaPdgCode = mcGamma.pdgCode();
            mcTrueMomGamma = {mcGamma.px(), mcGamma.py(), mcGamma.pz()};

            auto const& gammaMothers = mcGamma.template mothers_as<aod::McParticles>();
            if (!gammaMothers.empty()) {
              auto mcPi0 = gammaMothers.front();
              gammaMotherPdgCode = mcPi0.pdgCode();

              auto const& pi0Mothers = mcPi0.template mothers_as<aod::McParticles>();
              if (!pi0Mothers.empty()) {
                gammaGMotherPdgCode = pi0Mothers.front().pdgCode();
              }
            }
          }

          int photonSigmaIdx = findSigmaPlusMotherOfPhoton(mcPos, mcNeg);
          isSignal = (protonSigmaIdx >= 0 && photonSigmaIdx >= 0 && photonSigmaIdx == protonSigmaIdx);
        }

        if (isSignal) {
          mcTrueVtx = {mcProton.vx(), mcProton.vy(), mcProton.vz()};
          mcTrueMomProton = {mcProton.px(), mcProton.py(), mcProton.pz()};
        }
      }
    }

    auto fillCandStep = [&](int step) {
      histos.fill(HIST("Candidate/hSelectionCounter"), step);
      if constexpr (IsMC) {
        if (isSignal) {
          histos.fill(HIST("Candidate/True/hSelectionCounter"), step);
        }
      }
    };
    fillCandStep(0); // all pairs

    auto protonTrackParCov = getTrackParCov(protonTrack);

    std::array<float, 21> zeroCov{};
    auto photonTrackParCov = o2::track::TrackParCov({photon.x(), photon.y(), photon.z()}, {photon.px(), photon.py(), photon.pz()}, zeroCov, 0, true);
    photonTrackParCov.setAbsCharge(0);
    photonTrackParCov.setPID(o2::track::PID::Photon);

    int nCand = 0;
    try {
      nCand = fitter.process(protonTrackParCov, photonTrackParCov);
    } catch (...) {
      return;
    }
    if (nCand == 0 || !fitter.propagateTracksToVertex()) {
      return;
    }
    fillCandStep(1); // Vertex fit

    float fitChi2 = fitter.getChi2AtPCACandidate();
    float dcaProtonGamma = std::sqrt(fitChi2);
    histos.fill(HIST("Candidate/hDcaProtonGamma"), dcaProtonGamma);
    if constexpr (IsMC) {
      if (isSignal) {
        histos.fill(HIST("Candidate/True/hDcaProtonGamma"), dcaProtonGamma);
      }
    }
    if (dcaProtonGamma > candMaxDcaProtonGamma) {
      return;
    }
    fillCandStep(2); // DCA(p,gamma)

    std::array<float, 3> secVtx = fitter.getPCACandidatePos();
    float radius = std::hypot(secVtx[0], secVtx[1]);
    histos.fill(HIST("Candidate/hRadius"), radius);
    if constexpr (IsMC) {
      if (isSignal) {
        histos.fill(HIST("Candidate/True/hRadius"), radius);
        float vtxResFromMc = std::hypot(secVtx[0] - mcTrueVtx[0], secVtx[1] - mcTrueVtx[1], secVtx[2] - mcTrueVtx[2]);
        histos.fill(HIST("Candidate/True/hVertexResFromMcTruth"), vtxResFromMc);
      }
    }
    if (radius < candMinRadius || radius > candMaxRadius) {
      return;
    }
    fillCandStep(3); // radius

    // flight direction n and the decay-plane basis n, eIn, eOut
    std::array<float, 3> flightVec{secVtx[0] - pv[0], secVtx[1] - pv[1], secVtx[2] - pv[2]};
    std::array<float, 3> nHat = normalize3(flightVec);
    float flightDistance = std::sqrt(dot3(flightVec, flightVec));

    std::array<float, 3> pProton;
    std::array<float, 3> pGamma1;
    fitter.getTrack(0).getPxPyPzGlo(pProton);
    fitter.getTrack(1).getPxPyPzGlo(pGamma1);

    if constexpr (IsMC) {
      if (isSignal) {
        float trueProtonP = std::sqrt(dot3(mcTrueMomProton, mcTrueMomProton));
        float fitProtonP = std::sqrt(dot3(pProton, pProton));
        histos.fill(HIST("Candidate/True/hProtonMomResFromMcTruth"), (fitProtonP - trueProtonP) / trueProtonP);

        float trueGammaP = std::sqrt(dot3(mcTrueMomGamma, mcTrueMomGamma));
        float fitGammaP = std::sqrt(dot3(pGamma1, pGamma1));
        histos.fill(HIST("Candidate/True/hPhotonMomResFromMcTruth"), (fitGammaP - trueGammaP) / trueGammaP);
      }
    }

    float e1 = std::sqrt(dot3(pGamma1, pGamma1)); // |p_gamma1|
    float massPi0 = o2::constants::physics::MassPionNeutral;

    std::array<float, 3> eOut{};
    std::array<float, 3> eIn{};
    float coefA = 0.f, coefB = 0.f, coefK = 0.f, coefL = 0.f;
    float pGamma2In = 0.f, pGamma2Out = 0.f, tPerp2 = 0.f;
    float discriminant = -1.f;
    std::array<float, 3> nHatOrig = nHat;             // each retry perturbs around this
    constexpr float CoefADegenerateThreshold = 1e-6f; // guards against dividing by a near-zero quadratic coefficient

    if constexpr (IsMC) {
      if (isSignal) {
        float cosProtonFlight = dot3(pProton, nHatOrig) / std::sqrt(dot3(pProton, pProton));
        float protonFlightAngle = std::acos(std::clamp(cosProtonFlight, -1.f, 1.f));
        histos.fill(HIST("Candidate/True/hProtonFlightAngle"), protonFlightAngle);
      }
    }

    int discrIter = 0;
    for (; discrIter <= discrRetryMaxIter; ++discrIter) {
      if (discrIter > 0) {
        TVector3 nHatVec(nHatOrig[0], nHatOrig[1], nHatOrig[2]);
        float dTheta = mThetaResoFunc.GetRandom();
        float dPhi = mPhiResoFunc.GetRandom();
        nHatVec.SetMagThetaPhi(1., nHatVec.Theta() + dTheta, nHatVec.Phi() + dPhi);
        nHat = {static_cast<float>(nHatVec.X()), static_cast<float>(nHatVec.Y()), static_cast<float>(nHatVec.Z())};
      }

      eOut = normalize3(cross3(pProton, nHat)); // normal to the decay plane
      eIn = cross3(eOut, nHat);                 // in-plane, transverse to n

      float pProtonIn = dot3(pProton, eIn); // proton has no out-of-plane component, by construction
      float pGamma1Long = dot3(pGamma1, nHat);
      float pGamma1In = dot3(pGamma1, eIn);
      float pGamma1Out = dot3(pGamma1, eOut);

      // in-plane transverse momentum cancels between p, gamma1, gamma2
      pGamma2In = -(pProtonIn + pGamma1In);
      // out-of-plane momentum cancels between gamma1 and gamma2 alone
      pGamma2Out = -pGamma1Out;
      tPerp2 = pGamma2In * pGamma2In + pGamma2Out * pGamma2Out;

      // m_pi0^2 = 2*(|p_gamma1||p_gamma2| - p_gamma1.p_gamma2), solved for
      // the remaining unknown (p_gamma2 along n) -> quadratic coefA*x^2+coefB*x+coefC=0
      float dot1 = pGamma1In * pGamma2In + pGamma1Out * pGamma2Out;
      coefK = massPi0 * massPi0 + 2.f * dot1;
      coefL = 2.f * pGamma1Long;

      coefA = 4.f * e1 * e1 - coefL * coefL;
      if (std::abs(coefA) < CoefADegenerateThreshold) {
        continue;
      }
      coefB = -2.f * coefK * coefL;
      float coefC = 4.f * e1 * e1 * tPerp2 - coefK * coefK;

      discriminant = coefB * coefB - 4.f * coefA * coefC;
      if (discriminant >= 0.f) {
        break;
      }
    }

    histos.fill(HIST("Candidate/hDiscriminantRetryIter"), discrIter);
    histos.fill(HIST("Candidate/hDiscriminant"), discriminant);
    if constexpr (IsMC) {
      if (isSignal) {
        histos.fill(HIST("Candidate/True/hDiscriminantRetryIter"), discrIter);
        histos.fill(HIST("Candidate/True/hDiscriminant"), discriminant);
      }
    }
    if (discriminant < 0.f) {
      return;
    }
    fillCandStep(4); // real root

    // two roots from the quadratic, among both we keep the mass closest to the nominal Sigma+ mass
    float sqrtDisc = std::sqrt(discriminant);
    std::array<float, 2> roots{(-coefB + sqrtDisc) / (2.f * coefA), (-coefB - sqrtDisc) / (2.f * coefA)};

    bool haveCandidate = false;
    float bestMass = -999.f;
    std::array<float, 3> bestMomGamma2{};
    for (const float& pGamma2Long : roots) {
      std::array<float, 3> pGamma2{
        eIn[0] * pGamma2In + eOut[0] * pGamma2Out + nHat[0] * pGamma2Long,
        eIn[1] * pGamma2In + eOut[1] * pGamma2Out + nHat[1] * pGamma2Long,
        eIn[2] * pGamma2In + eOut[2] * pGamma2Out + nHat[2] * pGamma2Long};
      float eGamma2 = std::sqrt(tPerp2 + pGamma2Long * pGamma2Long);
      float eProton = std::sqrt(dot3(pProton, pProton) + o2::constants::physics::MassProton * o2::constants::physics::MassProton);
      std::array<float, 3> pTotal{pProton[0] + pGamma1[0] + pGamma2[0], pProton[1] + pGamma1[1] + pGamma2[1], pProton[2] + pGamma1[2] + pGamma2[2]};
      float eTotal = eProton + e1 + eGamma2;
      float mass2 = eTotal * eTotal - dot3(pTotal, pTotal);
      if (mass2 < 0.f) {
        continue;
      }
      float mass = std::sqrt(mass2);
      if (!haveCandidate || std::abs(mass - o2::constants::physics::MassSigmaPlus) < std::abs(bestMass - o2::constants::physics::MassSigmaPlus)) {
        haveCandidate = true;
        bestMass = mass;
        bestMomGamma2 = pGamma2;
      }
    }
    if (!haveCandidate) {
      return;
    }
    fillCandStep(5); // valid root

    std::array<float, 3> pSigma{pProton[0] + pGamma1[0] + bestMomGamma2[0], pProton[1] + pGamma1[1] + bestMomGamma2[1], pProton[2] + pGamma1[2] + bestMomGamma2[2]};
    float ptSigma = std::hypot(pSigma[0], pSigma[1]);

    histos.fill(HIST("Candidate/hMassSigmaPlus"), bestMass);
    histos.fill(HIST("Candidate/h2MassVsPt"), ptSigma, bestMass);
    if constexpr (IsMC) {
      if (isSignal) {
        histos.fill(HIST("Candidate/True/hMassSigmaPlus"), bestMass);
        histos.fill(HIST("Candidate/True/h2MassVsPt"), ptSigma, bestMass);
      }
    }
    fillCandStep(6); // filled

    // photon (V0) opening angle: angle between the e+/e- daughter momenta at their own reference point
    std::array<float, 3> pPosDau{posTrack.px(), posTrack.py(), posTrack.pz()};
    std::array<float, 3> pNegDau{negTrack.px(), negTrack.py(), negTrack.pz()};
    float photonOpeningAngle = std::acos(std::clamp(dot3(pPosDau, pNegDau) / std::sqrt(dot3(pPosDau, pPosDau) * dot3(pNegDau, pNegDau)), -1.f, 1.f));

    // photon pointing angle: angle between the fitted photon momentum and the line from its conversion point to the p-gamma decay vertex
    std::array<float, 3> convToDecVtx{secVtx[0] - photon.x(), secVtx[1] - photon.y(), secVtx[2] - photon.z()};
    float photonPointingAngle = std::acos(std::clamp(dot3(pGamma1, convToDecVtx) / std::sqrt(dot3(pGamma1, pGamma1) * dot3(convToDecVtx, convToDecVtx)), -1.f, 1.f));

    // photon DCA to PV: distance from the PV to the line through the conversion point along the photon momentum direction
    std::array<float, 3> convPoint{photon.x(), photon.y(), photon.z()};
    std::array<float, 3> photonDir = normalize3({photon.px(), photon.py(), photon.pz()});
    std::array<float, 3> pvToConv{pv[0] - convPoint[0], pv[1] - convPoint[1], pv[2] - convPoint[2]};
    std::array<float, 3> pvToConvCrossDir = cross3(pvToConv, photonDir);
    float photonDcaToPV = std::sqrt(dot3(pvToConvCrossDir, pvToConvCrossDir));

    if constexpr (IsMC) {
      sigmaPlusCandsMC(secVtx[0], secVtx[1], secVtx[2],
                       radius, flightDistance, dcaProtonGamma, fitChi2,
                       pProton[0], pProton[1], pProton[2],
                       pGamma1[0], pGamma1[1], pGamma1[2],
                       bestMomGamma2[0], bestMomGamma2[1], bestMomGamma2[2],
                       protonTrack.tpcNSigmaPr(), protonTrack.tofNSigmaPr(),
                       posTrack.tpcNSigmaEl(), negTrack.tpcNSigmaEl(),
                       photon.mGamma(), photon.alpha(), photon.qtarm(), photon.v0radius(),
                       photonOpeningAngle, photonPointingAngle, photonDcaToPV,
                       protonTrack.itsNCls(), protonTrack.tpcNClsFound(), protonTrack.dcaXY(), protonTrack.dcaZ(),
                       posTrack.itsNCls(), posTrack.tpcNClsFound(), negTrack.itsNCls(), negTrack.tpcNClsFound(),
                       isSignal,
                       protonPdgCode, protonMotherPdgCode,
                       gammaPdgCode, gammaMotherPdgCode, gammaGMotherPdgCode,
                       mcTrueVtx[0], mcTrueVtx[1], mcTrueVtx[2],
                       mcTrueMomProton[0], mcTrueMomProton[1], mcTrueMomProton[2],
                       mcTrueMomGamma[0], mcTrueMomGamma[1], mcTrueMomGamma[2]);
    } else {
      sigmaPlusCands(secVtx[0], secVtx[1], secVtx[2],
                     radius, flightDistance, dcaProtonGamma, fitChi2,
                     pProton[0], pProton[1], pProton[2],
                     pGamma1[0], pGamma1[1], pGamma1[2],
                     bestMomGamma2[0], bestMomGamma2[1], bestMomGamma2[2],
                     protonTrack.tpcNSigmaPr(), protonTrack.tofNSigmaPr(),
                     posTrack.tpcNSigmaEl(), negTrack.tpcNSigmaEl(),
                     photon.mGamma(), photon.alpha(), photon.qtarm(), photon.v0radius(),
                     photonOpeningAngle, photonPointingAngle, photonDcaToPV,
                     protonTrack.itsNCls(), protonTrack.tpcNClsFound(), protonTrack.dcaXY(), protonTrack.dcaZ(),
                     posTrack.itsNCls(), posTrack.tpcNClsFound(), negTrack.itsNCls(), negTrack.tpcNClsFound());
    }
  }

  void initCCDB(aod::BCs::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    o2::parameters::GRPMagField* grpmag = ccdb->getForRun<o2::parameters::GRPMagField>(grpmagPath, mRunNumber);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    mBz = grpmag->getNominalL3Field();
    fitter.setBz(mBz);
    LOG(info) << "Task initialized for run " << mRunNumber << " with magnetic field " << mBz << " kZG";
  }

  void processData(CollisionsFull const& collisions, aod::V0Datas const& v0s, TracksFull const& tracks, aod::BCs const&)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision.bc_as<aod::BCs>());
      histos.fill(HIST("hVertexZ"), collision.posZ());
      std::array<float, 3> pv{collision.posX(), collision.posY(), collision.posZ()};

      auto v0sThisCollision = v0s.sliceBy(v0PerCollision, collision.globalIndex());
      std::vector<aod::V0Datas::iterator> acceptedPhotons;
      for (const auto& v0 : v0sThisCollision) {
        if (selectPhoton<false, TracksFull>(v0)) {
          acceptedPhotons.push_back(v0);
        }
      }

      auto tracksThisCollision = tracks.sliceBy(tracksPerCollision, collision.globalIndex());
      std::vector<TracksFull::iterator> acceptedProtons;
      for (const auto& track : tracksThisCollision) {
        if (selectProton(track)) {
          acceptedProtons.push_back(track);
        }
      }

      for (const auto& photon : acceptedPhotons) {
        for (const auto& proton : acceptedProtons) {
          buildSigmaPlusCandidate<false, TracksFull>(proton, photon, pv);
        }
      }
    }
  }
  PROCESS_SWITCH(Sigmaplusbuilder, processData, "Process data", true);

  void processMc(CollisionsFullMC const& collisions, aod::V0Datas const& v0s, TracksFullMC const& tracks, aod::BCs const&, aod::McParticles const& mcParticles)
  {
    for (const auto& collision : collisions) {
      initCCDB(collision.bc_as<aod::BCs>());
      histos.fill(HIST("hVertexZ"), collision.posZ());
      std::array<float, 3> pv{collision.posX(), collision.posY(), collision.posZ()};

      auto v0sThisCollision = v0s.sliceBy(v0PerCollision, collision.globalIndex());
      std::vector<aod::V0Datas::iterator> acceptedPhotons;
      for (const auto& v0 : v0sThisCollision) {
        if (selectPhoton<true, TracksFullMC>(v0)) {
          acceptedPhotons.push_back(v0);

          histos.fill(HIST("MC/hPhotonTruthQA"), 0);
          auto posTrack = v0.template posTrack_as<TracksFullMC>();
          auto negTrack = v0.template negTrack_as<TracksFullMC>();
          if (posTrack.has_mcParticle() && negTrack.has_mcParticle()) {
            histos.fill(HIST("MC/hPhotonTruthQA"), 1);
            auto mcPos = posTrack.template mcParticle_as<aod::McParticles>();
            auto mcNeg = negTrack.template mcParticle_as<aod::McParticles>();

            auto const& posMothers = mcPos.template mothers_as<aod::McParticles>();
            auto const& negMothers = mcNeg.template mothers_as<aod::McParticles>();
            if (!posMothers.empty() && !negMothers.empty()) {
              auto mcGamma = posMothers.front();
              if (mcGamma.globalIndex() == negMothers.front().globalIndex() && mcGamma.pdgCode() == PDG_t::kGamma) {
                histos.fill(HIST("MC/hPhotonTruthQA"), 2);
                auto const& pi0Mothers = mcGamma.template mothers_as<aod::McParticles>();
                if (!pi0Mothers.empty() && std::abs(pi0Mothers.front().pdgCode()) == PDG_t::kPi0) {
                  histos.fill(HIST("MC/hPhotonTruthQA"), 3);
                  auto const& sigmaMothers = pi0Mothers.front().template mothers_as<aod::McParticles>();
                  if (!sigmaMothers.empty() && std::abs(sigmaMothers.front().pdgCode()) == PDG_t::kSigmaPlus) {
                    histos.fill(HIST("MC/hPhotonTruthQA"), 4);

                    histos.fill(HIST("Photon/True/hMass"), v0.mGamma());
                    histos.fill(HIST("Photon/True/hPt"), v0.pt());
                    histos.fill(HIST("Photon/True/hRadius"), v0.v0radius());
                    histos.fill(HIST("Photon/True/h2ArmenterosPodolanski"), v0.alpha(), v0.qtarm());
                    histos.fill(HIST("Photon/True/h2ConvPointXY"), v0.x(), v0.y());
                  }
                }
              }
            }
          }
        }
      }

      auto tracksThisCollision = tracks.sliceBy(tracksPerCollisionMC, collision.globalIndex());
      std::vector<TracksFullMC::iterator> acceptedProtons;
      for (const auto& track : tracksThisCollision) {
        if (selectProton(track)) {
          acceptedProtons.push_back(track);

          histos.fill(HIST("MC/hProtonTruthQA"), 0);
          if (track.has_mcParticle()) {
            histos.fill(HIST("MC/hProtonTruthQA"), 1);
            auto mcProton = track.template mcParticle_as<aod::McParticles>();
            if (findSigmaPlusMotherOfProton(mcProton) >= 0) {
              histos.fill(HIST("MC/hProtonTruthQA"), 2);

              histos.fill(HIST("Proton/True/hPt"), track.pt());
              histos.fill(HIST("Proton/True/h2TPCNSigmaVsPt"), track.pt(), track.tpcNSigmaPr());
              if (track.hasTOF()) {
                histos.fill(HIST("Proton/True/h2TOFNSigmaVsPt"), track.pt(), track.tofNSigmaPr());
              }
            }
          }
        }
      }

      for (const auto& photon : acceptedPhotons) {
        for (const auto& proton : acceptedProtons) {
          buildSigmaPlusCandidate<true, TracksFullMC>(proton, photon, pv);
        }
      }
    }

    // all generated Sigma+ -> p pi0 decays, regardless of reconstruction
    for (const auto& mcPart : mcParticles) {
      if (isSigmaPlusToProtonPi0(mcPart)) {
        histos.fill(HIST("MC/hGenSigmaPlusPt"), mcPart.pt());
      }
    }
  }
  PROCESS_SWITCH(Sigmaplusbuilder, processMc, "Process MC", false);

  void processFindable(CollisionsFullMC const& collisions, aod::V0Datas const& v0s, TracksFullMC const& tracks, aod::McParticles const&)
  {
    for (const auto& collision : collisions) {
      auto tracksThisCollision = tracks.sliceBy(tracksPerCollisionMC, collision.globalIndex());
      auto v0sThisCollision = v0s.sliceBy(v0PerCollision, collision.globalIndex());

      for (const auto& protonTrack : tracksThisCollision) {
        if (!protonTrack.has_mcParticle()) {
          continue;
        }

        auto mcProton = protonTrack.template mcParticle_as<aod::McParticles>();
        int sigmaIndex = findSigmaPlusMotherOfProton(mcProton);
        if (sigmaIndex < 0) {
          continue;
        }

        auto const& protonMothers = mcProton.template mothers_as<aod::McParticles>();
        if (protonMothers.empty()) {
          continue;
        }
        auto mcSigmaPlus = protonMothers.front();
        if (std::abs(mcSigmaPlus.y()) > 1) {
          continue;
        }
        histos.fill(HIST("Findable/h2ProtonPtVsSigmaPlusPt"), mcSigmaPlus.pt(), protonTrack.pt());
        std::vector<int> seenElectronMcIds;
        std::vector<int> seenElectronTrackIds;
        std::vector<int> seenPositronMcIds;
        std::vector<int> seenPositronTrackIds;

        for (const auto& electronTrack : tracksThisCollision) {
          if (!electronTrack.has_mcParticle()) {
            continue;
          }

          if (electronTrack.tpcNClsFound() < 90 || electronTrack.sign() > 0) {
            continue;
          }

          auto mcElectron = electronTrack.template mcParticle_as<aod::McParticles>();
          if (mcElectron.pdgCode() != PDG_t::kElectron) {
            continue;
          }

          int electronGammaIndex = -1;
          std::array<float, 3> electronVertex{};
          if (findSigmaPlusMotherOfConversionElectron(mcElectron, electronGammaIndex, electronVertex) != sigmaIndex) {
            continue;
          }

          for (const auto& positronTrack : tracksThisCollision) {
            if (!positronTrack.has_mcParticle()) {
              continue;
            }

            if (positronTrack.tpcNClsFound() < 90 || positronTrack.sign() < 0) {
              continue;
            }

            auto mcPositron = positronTrack.template mcParticle_as<aod::McParticles>();
            if (mcPositron.pdgCode() != PDG_t::kPositron) {
              continue;
            }

            int positronGammaIndex = -1;
            std::array<float, 3> positronVertex{};
            if (findSigmaPlusMotherOfConversionElectron(mcPositron, positronGammaIndex, positronVertex) != sigmaIndex || positronGammaIndex != electronGammaIndex) {
              continue;
            }

            bool sameConversionPoint = electronVertex[0] == positronVertex[0] &&
                                       electronVertex[1] == positronVertex[1] &&
                                       electronVertex[2] == positronVertex[2];
            if (!sameConversionPoint) {
              continue;
            }

            int electronMcId = mcElectron.globalIndex();
            int electronTrackId = electronTrack.globalIndex();
            int positronMcId = mcPositron.globalIndex();
            int positronTrackId = positronTrack.globalIndex();

            bool seenElectronMc = false;
            bool duplicateElectronTrack = false;
            for (int iSeen = 0; iSeen < static_cast<int>(seenElectronMcIds.size()); ++iSeen) {
              if (seenElectronMcIds[iSeen] == electronMcId) {
                seenElectronMc = true;
                duplicateElectronTrack |= (seenElectronTrackIds[iSeen] != electronTrackId);
              }
            }
            if (duplicateElectronTrack) {
              histos.fill(HIST("Findable/hDuplicateConversionTrackCounter"), 0);
              histos.fill(HIST("Findable/hDuplicateElectronTPCNClsFound"), electronTrack.tpcNClsFound());
            }
            if (!seenElectronMc) {
              seenElectronMcIds.push_back(electronMcId);
              seenElectronTrackIds.push_back(electronTrackId);
            }

            bool seenPositronMc = false;
            bool duplicatePositronTrack = false;
            for (int iSeen = 0; iSeen < static_cast<int>(seenPositronMcIds.size()); ++iSeen) {
              if (seenPositronMcIds[iSeen] == positronMcId) {
                seenPositronMc = true;
                duplicatePositronTrack |= (seenPositronTrackIds[iSeen] != positronTrackId);
              }
            }
            if (duplicatePositronTrack) {
              histos.fill(HIST("Findable/hDuplicateConversionTrackCounter"), 1);
              histos.fill(HIST("Findable/hDuplicatePositronTPCNClsFound"), positronTrack.tpcNClsFound());
            }
            if (!seenPositronMc) {
              seenPositronMcIds.push_back(positronMcId);
              seenPositronTrackIds.push_back(positronTrackId);
            }
            if (duplicateElectronTrack || duplicatePositronTrack) {
              continue;
            }

            std::array<float, 3> recoGammaMom{electronTrack.px() + positronTrack.px(), electronTrack.py() + positronTrack.py(), electronTrack.pz() + positronTrack.pz()};
            float recoGammaP = std::sqrt(dot3(recoGammaMom, recoGammaMom));
            constexpr float electronMass = o2::constants::physics::MassElectron;
            float electronP2 = electronTrack.px() * electronTrack.px() + electronTrack.py() * electronTrack.py() + electronTrack.pz() * electronTrack.pz();
            float positronP2 = positronTrack.px() * positronTrack.px() + positronTrack.py() * positronTrack.py() + positronTrack.pz() * positronTrack.pz();
            float pairEnergy = std::sqrt(electronP2 + electronMass * electronMass) + std::sqrt(positronP2 + electronMass * electronMass);
            float pairMass2 = pairEnergy * pairEnergy - recoGammaP * recoGammaP;
            histos.fill(HIST("Findable/hElectronPositronMass"), std::sqrt(std::max(pairMass2, 0.f)));

            auto const& gammaMothers = mcElectron.template mothers_as<aod::McParticles>();
            if (!gammaMothers.empty()) {
              auto mcGamma = gammaMothers.front();
              std::array<float, 3> trueGammaMom{mcGamma.px(), mcGamma.py(), mcGamma.pz()};
              float trueGammaP = std::sqrt(dot3(trueGammaMom, trueGammaMom));
              if (trueGammaP > 0.f) {
                histos.fill(HIST("Findable/hPhotonMomentumResolution"), (recoGammaP - trueGammaP) / trueGammaP);
              }
            }

            std::array<float, 3> photonFlightVec{electronVertex[0] - collision.posX(), electronVertex[1] - collision.posY(), electronVertex[2] - collision.posZ()};
            float flightNorm = std::sqrt(dot3(photonFlightVec, photonFlightVec));
            if (flightNorm > 0.f && recoGammaP > 0.f) {
              histos.fill(HIST("Findable/hPhotonCosPA"), dot3(photonFlightVec, recoGammaMom) / (flightNorm * recoGammaP));
            }

            histos.fill(HIST("Findable/hSigmaPlusPt"), mcSigmaPlus.pt());
            histos.fill(HIST("Findable/hConversionPairV0Presence"), 0);
            for (const auto& v0 : v0sThisCollision) {
              auto posTrack = v0.template posTrack_as<TracksFullMC>();
              auto negTrack = v0.template negTrack_as<TracksFullMC>();
              bool sameChargeMatched = posTrack.globalIndex() == positronTrack.globalIndex() && negTrack.globalIndex() == electronTrack.globalIndex();
              if (sameChargeMatched) {
                histos.fill(HIST("Findable/hConversionPairV0Presence"), 1);
                break;
              }
            }
            histos.fill(HIST("Findable/hElectronPt"), electronTrack.pt());
            histos.fill(HIST("Findable/hPositronPt"), positronTrack.pt());
            histos.fill(HIST("Findable/hPhotonConversionRadius"), std::hypot(electronVertex[0], electronVertex[1]));
            fillFindableTrackDetectors<true>(electronTrack);
            fillFindableTrackDetectors<false>(positronTrack);
          }
        }
      }
    }
  }
  PROCESS_SWITCH(Sigmaplusbuilder, processFindable, "Process findable MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<Sigmaplusbuilder>(cfgc)};
}
