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

#include "PWGEM/Dilepton/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"
#include "PWGLF/DataModel/LFKinkDecayTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/svPoolCreator.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TPCVDriftManager.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/MatLayerCylSet.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/Primitive2D.h>
#include <ReconstructionDataFormats/PID.h>
#include <ReconstructionDataFormats/Track.h>

#include <TAxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TMCProcess.h>
#include <TPDGCode.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;

using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA,
                             aod::pidTPCFullEl, aod::pidTPCFullPr, aod::pidTOFFullPr>;
using TracksFullMC = soa::Join<TracksFull, aod::McTrackLabels>;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels>;

struct Sigmaplusbuilder {

  // event selection
  Configurable<float> cutZVertex{"cutZVertex", 10.0f, "Accepted z-vertex range (cm)"};

  // photon (PCM) selection
  Configurable<float> photonMaxMass{"photonMaxMass", 0.05, "Max photon mass (GeV/c^2)"};
  Configurable<float> photonMinRapidity{"photonMinRapidity", -1.0, "Min photon rapidity"};
  Configurable<float> photonMaxRapidity{"photonMaxRapidity", 1.0, "Max photon rapidity"};
  Configurable<float> photonDauEtaMin{"photonDauEtaMin", -1.0, "Min eta of photon daughter tracks"};
  Configurable<float> photonDauEtaMax{"photonDauEtaMax", 1.0, "Max eta of photon daughter tracks"};
  Configurable<float> photonMinRadius{"photonMinRadius", 3.0, "Min photon conversion radius (cm)"};
  Configurable<float> photonMaxRadius{"photonMaxRadius", 115., "Max photon conversion radius (cm)"};
  Configurable<float> photonMinV0cospa{"photonMinV0cospa", 0.9, "Min photon cosine of pointing angle to the PV"};
  Configurable<float> photonMaxDCAV0Dau{"photonMaxDCAV0Dau", 1.5, "Max DCA between photon daughters (cm)"};
  Configurable<float> photonMaxOpeningAngle{"photonMaxOpeningAngle", 0.4, "Max opening angle between the photon's e+/e- daughter momenta (rad)"};
  Configurable<float> photonMaxDeltaTheta{"photonMaxDeltaTheta", 0.15, "Max |theta_pos - theta_neg| of the photon's daughter tracks (rad)"};
  Configurable<float> photonMaxQt{"photonMaxQt", 0.05, "Max Armenteros qT for photons (GeV/c)"};
  Configurable<float> photonMaxAlpha{"photonMaxAlpha", 0.95, "Max |Armenteros alpha| for photons"};
  Configurable<float> photonDauMinTPCNSigmaEl{"photonDauMinTPCNSigmaEl", -3., "Min TPC nSigma_el of the photon daughters"};
  Configurable<float> photonDauMaxTPCNSigmaEl{"photonDauMaxTPCNSigmaEl", 3., "Max TPC nSigma_el of the photon daughters"};
  Configurable<float> photonDauMaxPt{"photonDauMaxPt", 0.5, "Max pT of the photon daughters (GeV/c)"};
  Configurable<float> photonDauMinTpcNCls{"photonDauMinTpcNCls", 30, "Min number of found TPC clusters for the photon daughter tracks"};

  // photon daughters paired by this task instead of taken from V0Datas (V0Datas assumes photons from primary vertex)
  Configurable<bool> useCustomVertexer{"useCustomVertexer", false, "Pair the photon's e+/e- daughter tracks in this task instead of reading V0Datas"};
  Configurable<bool> photonSkipAmbiTracks{"photonSkipAmbiTracks", false, "Skip ambiguous tracks when pairing photon daughters"};
  Configurable<float> photonPoolTimeMarginNS{"photonPoolTimeMarginNS", 800., "Time margin (ns) added to a daughter track's time range when matching it to collisions"};
  Configurable<float> photonMaxDXYIni{"photonMaxDXYIni", 4., "Max xy distance (cm) between the two daughter tracks at the start of the photon vertex fit"};
  Configurable<float> photonMaxCircleTouchDist{"photonMaxCircleTouchDist", 4., "Max |circle centre distance - (R1 + R2)| (cm) of the two daughter tracks, applied before the fit"};

  // proton selection
  Configurable<float> protonMinPt{"protonMinPt", 0.4, "Minimum proton pT (GeV/c)"};
  Configurable<float> protonMaxEta{"protonMaxEta", 0.9, "Maximum |eta| for proton track"};
  Configurable<float> protonMinTpcNCls{"protonMinTpcNCls", 80, "Min number of found TPC clusters for the proton track"};
  Configurable<float> protonMaxTPCNSigma{"protonMaxTPCNSigma", 4, "Max |TPC nSigma_pr| for proton"};
  Configurable<float> protonMaxTOFNSigma{"protonMaxTOFNSigma", 4, "Max |TOF nSigma_pr| for proton, if TOF available"};
  Configurable<float> protonPtMinRequireTOF{"protonPtMinRequireTOF", 0.75, "Above this pT, require TOF PID for proton"};
  Configurable<bool> protonRequireTofHit{"protonRequireTofHit", false, "Above protonPtMinRequireTOF, reject a proton with no TOF hit at all"};
  Configurable<float> protonMinDcaToPV{"protonMinDcaToPV", 0.005, "Min DCAxy of the proton track to the PV (cm)"};
  Configurable<float> protonMaxDcaToPV{"protonMaxDcaToPV", 5.0, "Max DCAxy of the proton track to the PV (cm)"};

  // proton-photon candidate selection
  Configurable<float> candVertexProtonWeight{"candVertexProtonWeight", 1.0, "Sigma+ vertex between the proton's (1) and the photon's (0) point of closest approach"};
  Configurable<float> candMaxDcaProtonGamma{"candMaxDcaProtonGamma", 1.5, "Max DCA between proton and photon at the fitted vertex (cm)"};
  Configurable<float> candMaxDcaToPV{"candMaxDcaToPV", 20., "Max DCA of the candidate's total (reconstructed) momentum line to the PV (cm)"};
  Configurable<bool> candRejectNegRootCenter{"candRejectNegRootCenter", true, "Reject candidates with rootCenter<0"};
  Configurable<float> candMaxRootCenter{"candMaxRootCenter", 15, "Max rootCenter=-coefB/(2*coefA) (GeV/c)"};
  Configurable<float> candMinAntiSigmaPointingAngle{"candMinAntiSigmaPointingAngle", 0.0, "Min AntiSigmaPointingAngle (rad)"};
  Configurable<float> candMaxAntiSigmaPointingAngle{"candMaxAntiSigmaPointingAngle", 10., "Max AntiSigmaPointingAngle (rad)"};
  Configurable<float> candMaxSigmaMass{"candMaxSigmaMass", 1.35, "Max reconstructed Sigma+ candidate mass (GeV/c^2)"};
  Configurable<float> candMaxRapidity{"candMaxRapidity", 0.9, "Max |rapidity| of the reconstructed Sigma+ candidate"};
  Configurable<float> candMinRadius{"candMinRadius", 1.0, "Min candidate decay radius (cm)"};
  Configurable<float> candMaxRadius{"candMaxRadius", 100., "Max candidate decay radius (cm)"};
  Configurable<float> candMinFlightDistance{"candMinFlightDistance", 0.0, "Min 3D distance from PV to candidate decay vertex (cm)"};
  Configurable<float> candMaxFlightDistance{"candMaxFlightDistance", 250.0, "Max 3D distance from PV to candidate decay vertex (cm)"};
  Configurable<float> candMaxPhotonOpeningAngle{"candMaxPhotonOpeningAngle", 3.15, "Max photon opening angle (rad), recomputed from the daughters' track momenta"};
  Configurable<float> candMaxPhotonPointingAngle{"candMaxPhotonPointingAngle", 0.5, "Max angle between the photon momentum and the decay-vertex-to-conversion-point line (rad)"};
  Configurable<float> candMaxPhotonDcaToPV{"candMaxPhotonDcaToPV", 250.0, "Max DCA of the photon's flight line to the PV (cm)"};
  Configurable<bool> candDeduplicatePhotons{"candDeduplicatePhotons", true, "Per timeframe, write only the candidate with the smallest DCA to PV among candidates sharing the same photon"};

  // tilt search: if the measured flight direction gives no real solution for the missing photon,
  // tilted flight directions are tried starting with the smallest tilt
  Configurable<float> candMaxTilt{"candMaxTilt", 0.3, "Max tilt of the flight direction searched for a real root (rad)"};
  Configurable<float> candTiltStep{"candTiltStep", 0.005, "Tilt step between the rings of tried flight directions (rad)"};
  Configurable<int> candTiltNAzimuth{"candTiltNAzimuth", 36, "Flight directions tried per ring"};

  // MC
  Configurable<float> cutRapMotherMC{"cutRapMotherMC", 1.0f, "Rapidity cut for generated mother Sigma+ in MC"};
  Configurable<float> cutPtGenMC{"cutPtGenMC", 0.5f, "Minimum pT for generated Sigma+ in MC"};

  Configurable<bool> fillSlimTables{"fillSlimTables", false, "write the slim candidate tables instead of the full ones"};

  Configurable<std::string> ccdbPath{"ccdbPath", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "CCDB path of the material budget LUT"};

  Produces<aod::SigmaPlusCands> sigmaPlusCands;
  Produces<aod::SigmaPlusCandsMC> sigmaPlusCandsMC;
  Produces<aod::SlimSigmaPlusCands> slimSigmaPlusCands;
  Produces<aod::SlimSigmaPlusCandsMC> slimSigmaPlusCandsMC;

  Service<o2::ccdb::BasicCCDBManager> ccdb{};
  o2::vertexing::DCAFitterN<2> fitter;
  int mRunNumber = 0;
  float mBz = 0;
  o2::base::MatLayerCylSet* mLut = nullptr;
  o2::aod::common::TPCVDriftManager mVDriftMgr;

  // daughter pairing (electrons and positrons with the range of collisions each is compatible with)
  svPoolCreator svPhotonPoolCreator{PDG_t::kElectron, PDG_t::kPositron};
  std::vector<TrackCand> mElectronPool;
  std::vector<TrackCand> mPositronPool;
  std::vector<bool> mGoodCollision; // collisions passing the event selection

  // values of a candidate, kept until the end of the timeframe and then written to the tables
  struct SigmaPlusCandidate {
    uint64_t photonId = 0;
    bool isSignal = false;
    int photonLegsWithoutIts = 0;
    int matchedSigmaId = -1;
    std::array<float, 3> decVtx{};
    float radius = 0.f;
    float flightDistance = 0.f;
    float dcaProtonGamma = 0.f;
    float dcaToPV = 0.f;
    float tiltAngle = 0.f;
    float rootCenter = 0.f;
    float antiSigmaPointingAngle = 0.f;
    std::array<float, 3> pProton{};
    std::array<float, 3> pGamma1{};
    std::array<float, 3> pGamma2{};
    float protonTpcNSigma = 0.f;
    float protonTofNSigma = 0.f;
    int protonSign = 0;
    uint8_t protonItsNCls = 0;
    int16_t protonTpcNCls = 0;
    float protonDcaXY = 0.f;
    float protonDcaZ = 0.f;
    float photonMass = 0.f;
    float photonAlpha = 0.f;
    float photonQt = 0.f;
    float photonRadius = 0.f;
    float photonOpeningAngle = 0.f;
    float photonPointingAngle = 0.f;
    float photonDcaDau = 0.f;
    float photonCosPAToPV = 0.f;
    float photonDcaXYToPV = 0.f;
    float photonDcaZToPV = 0.f;
    float photonPsiPair = 0.f;
    float posTpcNSigmaEl = 0.f;
    float negTpcNSigmaEl = 0.f;
    uint8_t posItsNCls = 0;
    uint8_t negItsNCls = 0;
    int16_t posTpcNCls = 0;
    int16_t negTpcNCls = 0;
    int16_t posTpcNClsFindable = 0;
    int16_t negTpcNClsFindable = 0;
    float posTpcChi2NCl = 0.f;
    float negTpcChi2NCl = 0.f;
    float posDcaXY = 0.f;
    float posDcaZ = 0.f;
    float negDcaXY = 0.f;
    float negDcaZ = 0.f;
    // MC
    bool collisionIdCheck = false;
    bool protonIsSignal = false;
    bool photonIsSignal = false;
    std::array<float, 3> decVtxMC{};
    std::array<float, 3> pProtonMC{};
    std::array<float, 3> pGammaMC{};
    std::array<float, 3> pSigmaPlusMC{};
    float decayRadiusMC = -999.f;
    float massMC = -999.f;
  };
  std::vector<SigmaPlusCandidate> mCandidatesOfTimeframe;

  // MC: generated Sigma+ -> p pi0 in acceptance of current timeframe and whether a true candidate of it was written
  std::unordered_map<int, bool> mGenSigmaWritten;

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
    fitter.setMaxDXYIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);

    mVDriftMgr.init(&ccdb->instance());
    svPhotonPoolCreator.setTimeMargin(photonPoolTimeMarginNS);
    if (photonSkipAmbiTracks) {
      svPhotonPoolCreator.setSkipAmbiTracks();
    }

    const AxisSpec axisVertexZ{100, -15., 15., "vrtx_{Z} (cm)"};

    const AxisSpec axisPhotonSel{14, -0.5, 13.5, "selection step"};
    const AxisSpec axisPhotonMass{200, 0., 0.3, "m_{#gamma} (GeV/c^{2})"};
    const AxisSpec axisPhotonPt{100, 0., 5., "#it{p}_{T,#gamma} (GeV/c)"};
    const AxisSpec axisPhotonRadius{200, 0., 200., "R_{conv} (cm)"};
    const AxisSpec axisAlpha{100, -1., 1., "#alpha_{AP}"};
    const AxisSpec axisQt{100, 0., 0.3, "q_{T,AP} (GeV/c)"};
    const AxisSpec axisConvXY{200, -100., 100., "conv. point (cm)"};
    const AxisSpec axisNSigmaEl{100, -10., 10., "n#sigma_{el}"};
    const AxisSpec axisTpcNCls{160, -0.5, 159.5, "TPC clusters"};
    const AxisSpec axisPhotonOpeningAngle{180, 0., 3.15, "opening angle (rad)"};
    const AxisSpec axisPhotonDeltaTheta{200, -1., 1., "#Delta#theta (rad)"};
    const AxisSpec axisPhotonDcaV0Dau{200, 0., 10., "DCA between photon daughters (cm)"};
    const AxisSpec axisCollisionsPerPhoton{100, 0.5, 100.5, "collisions a photon is kept in"};
    const AxisSpec axisTpcTimeRangeNColl{200, 0.5, 200.5, "collisions in the TPC time range of a TPC-only daughter"};

    const AxisSpec axisProtonSel{7, -0.5, 6.5, "selection step"};
    const AxisSpec axisProtonPt{100, 0., 5., "#it{p}_{T,p} (GeV/c)"};
    const AxisSpec axisNSigma{100, -5., 5., "n#sigma"};
    const AxisSpec axisProtonDcaToPV{1000, 0., 5., "DCA_{xy,p} to PV (cm)"};

    const AxisSpec axisCandSel{15, -0.5, 14.5, "selection step"};
    const AxisSpec axisPhotonLegType{3, -0.5, 2.5, "photon daughters without ITS"};
    const AxisSpec axisCandPhotonDcaToPV{250, 0., 250., "DCA_{#gamma-line} to PV (cm)"};
    const AxisSpec axisDca{100, 0., 10., "DCA(p,#gamma) (cm)"};
    const AxisSpec axisDcaToPV{500, 0., 0.5, "DCA_{cand} to PV (cm)"};
    const AxisSpec axisCandRadius{200, 0., 200., "R_{dec} (cm)"};
    const AxisSpec axisFlightDistance{250, 0., 250., "|SV-PV| (cm)"};
    const AxisSpec axisTilt{200, 0., 0.5, "flight direction tilt for a real root (rad)"};
    const AxisSpec axisRootCenter{400, -10., 10., "-b/(2a) (GeV/#it{c})"};
    const AxisSpec axisAntiSigmaPA{180, 0., 3.15, "AntiPA (rad)"};
    const AxisSpec axisMassSigma{200, 1.0, 1.4, "m_{p#gamma#gamma} (GeV/#it{c}^{2})"};
    const AxisSpec axisRapidity{200, -2., 2., "y_{#Sigma^{+}}"};
    const AxisSpec axisSigmaPt{100, 0., 6., "#it{p}_{T,#Sigma^{+}} (GeV/#it{c})"};
    const AxisSpec axisCandidatesPerPhoton{50, 0.5, 50.5, "candidates sharing one photon"};

    const AxisSpec axisSigmaSel{6, -0.5, 5.5, "step"};

    const std::vector<std::string> photonSteps{"All", "Neg eta", "Pos eta", "TPC clusters", "TPC nSigma_{el}, p_{T}", "DCA daughters", "Radius",
                                               "Opening angle", "Delta theta", "CosPA", "Rapidity", "Qt", "Alpha", "Mass"};
    const std::vector<std::string> protonSteps{"All", "Pt", "Eta", "TPC clusters", "TPC nSigma", "TOF nSigma", "DCA to PV"};
    const std::vector<std::string> candSteps{"All pairs", "Autocorrelation", "Vertex fit", "DCA(p,#gamma)", "Radius", "Flight distance", "Real root",
                                             "Valid root", "DCA to PV", "Mass", "Rapidity", "Photon opening angle", "Photon pointing angle", "Photon DCA to PV", "Written"};
    const std::vector<std::string> sigmaSteps{"Generated", "Proton reconstructed", "Proton selected", "Photon daughters reconstructed", "Photon selected", "Candidate written"};
    auto setBinLabels = [](TAxis* axis, const std::vector<std::string>& labels) {
      for (size_t i = 0; i < labels.size(); ++i) {
        axis->SetBinLabel(i + 1, labels[i].c_str());
      }
    };
    histos.add("hVertexZ", "hVertexZ", kTH1F, {axisVertexZ});

    setBinLabels(histos.add<TH1>("Photon/Inclusive/hSelectionCounter", "hSelectionCounter", kTH1D, {axisPhotonSel})->GetXaxis(), photonSteps);
    histos.add("Photon/Inclusive/h2TPCNClsPosVsPt", "h2TPCNClsPosVsPt", kTH2F, {axisPhotonPt, axisTpcNCls});
    histos.add("Photon/Inclusive/h2TPCNClsNegVsPt", "h2TPCNClsNegVsPt", kTH2F, {axisPhotonPt, axisTpcNCls});
    histos.add("Photon/Inclusive/h2TPCNSigmaElPosVsPt", "h2TPCNSigmaElPosVsPt", kTH2F, {axisPhotonPt, axisNSigmaEl});
    histos.add("Photon/Inclusive/h2TPCNSigmaElNegVsPt", "h2TPCNSigmaElNegVsPt", kTH2F, {axisPhotonPt, axisNSigmaEl});
    histos.add("Photon/Inclusive/hDcaV0Daughters", "hDcaV0Daughters", kTH1F, {axisPhotonDcaV0Dau});
    histos.add("Photon/Inclusive/hOpeningAngle", "hOpeningAngle", kTH1F, {axisPhotonOpeningAngle});
    histos.add("Photon/Inclusive/hDeltaTheta", "hDeltaTheta", kTH1F, {axisPhotonDeltaTheta});
    histos.add("Photon/Inclusive/hMass", "hMass", kTH1F, {axisPhotonMass});
    histos.add("Photon/Inclusive/hPt", "hPt", kTH1F, {axisPhotonPt});
    histos.add("Photon/Inclusive/hRadius", "hRadius", kTH1F, {axisPhotonRadius});
    histos.add("Photon/Inclusive/h2ArmenterosPodolanski", "h2ArmenterosPodolanski", kTH2F, {axisAlpha, axisQt});
    histos.add("Photon/Inclusive/h2ConvPointXY", "h2ConvPointXY", kTH2F, {axisConvXY, axisConvXY});
    histos.add("Photon/Inclusive/hCollisionsPerPhoton", "hCollisionsPerPhoton", kTH1F, {axisCollisionsPerPhoton});

    setBinLabels(histos.add<TH1>("Proton/Inclusive/hSelectionCounter", "hSelectionCounter", kTH1D, {axisProtonSel})->GetXaxis(), protonSteps);
    histos.add("Proton/Inclusive/hPt", "hPt", kTH1F, {axisProtonPt});
    histos.add("Proton/Inclusive/h2TPCNSigmaVsPt", "h2TPCNSigmaVsPt", kTH2F, {axisProtonPt, axisNSigma});
    histos.add("Proton/Inclusive/h2TOFNSigmaVsPt", "h2TOFNSigmaVsPt", kTH2F, {axisProtonPt, axisNSigma});
    histos.add("Proton/Inclusive/h2TPCNClsVsPt", "h2TPCNClsVsPt", kTH2F, {axisProtonPt, axisTpcNCls});
    histos.add("Proton/Inclusive/h2DcaToPVVsPt", "h2DcaToPVVsPt", kTH2F, {axisProtonPt, axisProtonDcaToPV});

    setBinLabels(histos.add<TH1>("Candidate/Inclusive/hSelectionCounter", "hSelectionCounter", kTH1D, {axisCandSel})->GetXaxis(), candSteps);
    setBinLabels(histos.add<TH2>("Candidate/Inclusive/h2SelectionCounterVsLegType", "h2SelectionCounterVsLegType", kTH2D, {axisCandSel, axisPhotonLegType})->GetXaxis(), candSteps);
    histos.add("Candidate/Inclusive/hDcaProtonGamma", "hDcaProtonGamma", kTH1F, {axisDca});
    histos.add("Candidate/Inclusive/hRadius", "hRadius", kTH1F, {axisCandRadius});
    histos.add("Candidate/Inclusive/hFlightDistance", "hFlightDistance", kTH1F, {axisFlightDistance});
    histos.add("Candidate/Inclusive/hTiltAngle", "hTiltAngle", kTH1F, {axisTilt});
    histos.add("Candidate/Inclusive/hRootCenter", "hRootCenter", kTH1F, {axisRootCenter});
    histos.add("Candidate/Inclusive/hDcaToPV", "hDcaToPV", kTH1F, {axisDcaToPV});
    histos.add("Candidate/Inclusive/hAntiSigmaPointingAngle", "hAntiSigmaPointingAngle", kTH1F, {axisAntiSigmaPA});
    histos.add("Candidate/Inclusive/hMassSigmaPlus", "hMassSigmaPlus", kTH1F, {axisMassSigma});
    histos.add("Candidate/Inclusive/h2MassVsPt", "h2MassVsPt", kTH2F, {axisSigmaPt, axisMassSigma});
    histos.add("Candidate/Inclusive/h2MassVsTilt", "h2MassVsTilt", kTH2F, {axisTilt, axisMassSigma});
    histos.add("Candidate/Inclusive/hRapidity", "hRapidity", kTH1F, {axisRapidity});
    histos.add("Candidate/Inclusive/hPhotonOpeningAngle", "hPhotonOpeningAngle", kTH1F, {axisPhotonOpeningAngle});
    histos.add("Candidate/Inclusive/hPhotonPointingAngle", "hPhotonPointingAngle", kTH1F, {axisPhotonOpeningAngle});
    histos.add("Candidate/Inclusive/hPhotonDcaToPV", "hPhotonDcaToPV", kTH1F, {axisCandPhotonDcaToPV});

    if (doprocessMc || doprocessFindable) {
      histos.addClone("Photon/Inclusive/", "Photon/True/");
    }
    if (doprocessMc) {
      histos.addClone("Proton/Inclusive/", "Proton/True/");
      histos.addClone("Candidate/Inclusive/", "Candidate/True/");
    }

    histos.add("Photon/Inclusive/hTpcTimeRangeNColl", "hTpcTimeRangeNColl", kTH1F, {axisTpcTimeRangeNColl});
    histos.add("Candidate/Inclusive/hCandidatesPerPhoton", "hCandidatesPerPhoton", kTH1F, {axisCandidatesPerPhoton});

    if (doprocessMc) {
      histos.add("MC/hGenSigmaPlusPt", "MC/hGenSigmaPlusPt", kTH1F, {axisSigmaPt});
      setBinLabels(histos.add<TH1>("MC/hSigmaPlusCounter", "MC/hSigmaPlusCounter", kTH1D, {axisSigmaSel})->GetXaxis(), sigmaSteps);
    }

    if (doprocessFindable) {
      const AxisSpec axisMomentum{100, 0., 10., "#it{p} (GeV/#it{c})"};
      const AxisSpec axisDetectorPresence{3, -0.5, 2.5, "detector"};
      const AxisSpec axisDuplicateTrack{2, -0.5, 1.5, "track"};
      const AxisSpec axisTPCClusters{160, -0.5, 159.5, "TPC clusters"};
      const AxisSpec axisPhotonMomResolution{200, -1., 1., "(#it{p}_{reco,#gamma} - #it{p}_{MC,#gamma})/#it{p}_{MC,#gamma}"};
      const AxisSpec axisPhotonCosPA{200, -1., 1., "cosPA_{#gamma}"};
      const AxisSpec axisPairMass{200, 0., 0.1, "m_{e^{+}e^{-}} (GeV/#it{c}^{2})"};
      const AxisSpec axisV0Presence{2, -0.5, 1.5, "V0 match"};
      histos.add("Findable/hSigmaPlusPt", "Findable/hSigmaPlusPt", kTH1F, {axisSigmaPt});
      histos.add("Findable/h2ProtonPtVsSigmaPlusPt", "Findable/h2ProtonPtVsSigmaPlusPt", kTH2F, {axisSigmaPt, axisMomentum});
      histos.add("Findable/hElectronPt", "Findable/hElectronPt", kTH1F, {axisMomentum});
      histos.add("Findable/hPositronPt", "Findable/hPositronPt", kTH1F, {axisMomentum});
      histos.add("Findable/hElectronPositronMass", "Findable/hElectronPositronMass", kTH1F, {axisPairMass});
      histos.add("Findable/hPhotonConversionRadius", "Findable/hPhotonConversionRadius", kTH1F, {axisPhotonRadius});
      histos.add("Findable/hPhotonMomentumResolution", "Findable/hPhotonMomentumResolution", kTH1F, {axisPhotonMomResolution});
      histos.add("Findable/hPhotonCosPA", "Findable/hPhotonCosPA", kTH1F, {axisPhotonCosPA});
      histos.add("Findable/hConversionPairV0Presence", "Findable/hConversionPairV0Presence", kTH1F, {axisV0Presence});
      histos.add("Findable/hPhotonSearchPresence", "Findable/hPhotonSearchPresence", kTH1F, {axisV0Presence});
      histos.add("Findable/hDuplicateConversionTrackCounter", "Findable/hDuplicateConversionTrackCounter", kTH1F, {axisDuplicateTrack});
      histos.add("Findable/hDuplicateElectronTPCNClsFound", "Findable/hDuplicateElectronTPCNClsFound", kTH1F, {axisTPCClusters});
      histos.add("Findable/hDuplicatePositronTPCNClsFound", "Findable/hDuplicatePositronTPCNClsFound", kTH1F, {axisTPCClusters});
      histos.add("Findable/hElectronDetectorPresence", "Findable/hElectronDetectorPresence", kTH1F, {axisDetectorPresence});
      histos.add("Findable/hPositronDetectorPresence", "Findable/hPositronDetectorPresence", kTH1F, {axisDetectorPresence});

      setBinLabels(histos.get<TH1>(HIST("Findable/hConversionPairV0Presence"))->GetXaxis(), {"valid pair", "in V0"});
      setBinLabels(histos.get<TH1>(HIST("Findable/hPhotonSearchPresence"))->GetXaxis(), {"valid pair", "passed photon selection"});
      setBinLabels(histos.get<TH1>(HIST("Findable/hDuplicateConversionTrackCounter"))->GetXaxis(), {"e^{-}", "e^{+}"});
      setBinLabels(histos.get<TH1>(HIST("Findable/hElectronDetectorPresence"))->GetXaxis(), {"ITS", "TPC", "TOF"});
      setBinLabels(histos.get<TH1>(HIST("Findable/hPositronDetectorPresence"))->GetXaxis(), {"ITS", "TPC", "TOF"});
    }
  }

  // photon (electron/positron conversion pair) candidate
  template <typename TTrack>
  struct PhotonCand {
    float x = 0.f, y = 0.f, z = 0.f;
    float px = 0.f, py = 0.f, pz = 0.f;
    float mGamma = 0.f;
    float alpha = 0.f;
    float qtarm = 0.f;
    float radius = 0.f;
    float dcaDau = 0.f;
    float cosPAToPV = 0.f;
    TTrack negTrack;
    TTrack posTrack;
  };

  // MC counter of generated Sigma+ -> p pi0 in acceptance
  void fillSigmaCounter(int sigmaId, int step)
  {
    if (mGenSigmaWritten.contains(sigmaId)) {
      histos.fill(HIST("MC/hSigmaPlusCounter"), step);
    }
  }

  template <typename TName, typename... Ts>
  void fillPhotonHist(const TName& name, bool isSignal, Ts... values)
  {
    histos.fill(HIST("Photon/Inclusive/") + name, values...);
    if (isSignal) {
      histos.fill(HIST("Photon/True/") + name, values...);
    }
  }

  template <typename TName, typename... Ts>
  void fillProtonHist(const TName& name, bool isSignal, Ts... values)
  {
    histos.fill(HIST("Proton/Inclusive/") + name, values...);
    if (isSignal) {
      histos.fill(HIST("Proton/True/") + name, values...);
    }
  }

  template <typename TName, typename... Ts>
  void fillCandHist(const TName& name, bool isSignal, Ts... values)
  {
    histos.fill(HIST("Candidate/Inclusive/") + name, values...);
    if (isSignal) {
      histos.fill(HIST("Candidate/True/") + name, values...);
    }
  }

  // photon candidates from V0Datas
  template <bool IsMC, typename TV0s, typename TTracks>
  std::vector<PhotonCand<typename TTracks::iterator>> findPhotonsFromV0s(const TV0s& v0s, const TTracks&, const std::array<float, 3>& pv)
  {
    std::vector<PhotonCand<typename TTracks::iterator>> photons;

    for (const auto& v0 : v0s) {
      auto posTrack = v0.template posTrack_as<TTracks>();
      auto negTrack = v0.template negTrack_as<TTracks>();

      int sigmaId = -1;
      if constexpr (IsMC) {
        if (posTrack.has_mcParticle() && negTrack.has_mcParticle()) {
          sigmaId = findSigmaPlusMotherOfPhoton(posTrack.template mcParticle_as<aod::McParticles>(), negTrack.template mcParticle_as<aod::McParticles>());
        }
      }
      bool isSignal = sigmaId >= 0;

      auto fillPhotonStep = [&](int step) {
        fillPhotonHist(HIST("hSelectionCounter"), isSignal, step);
      };
      fillPhotonStep(0);

      if (negTrack.eta() < photonDauEtaMin || negTrack.eta() > photonDauEtaMax) {
        continue;
      }
      fillPhotonStep(1);

      if (posTrack.eta() < photonDauEtaMin || posTrack.eta() > photonDauEtaMax) {
        continue;
      }
      fillPhotonStep(2);

      fillPhotonHist(HIST("h2TPCNClsPosVsPt"), isSignal, posTrack.pt(), posTrack.tpcNClsFound());
      fillPhotonHist(HIST("h2TPCNClsNegVsPt"), isSignal, negTrack.pt(), negTrack.tpcNClsFound());
      if (posTrack.tpcNClsFound() < photonDauMinTpcNCls || negTrack.tpcNClsFound() < photonDauMinTpcNCls) {
        continue;
      }
      fillPhotonStep(3);

      fillPhotonHist(HIST("h2TPCNSigmaElPosVsPt"), isSignal, posTrack.pt(), posTrack.tpcNSigmaEl());
      fillPhotonHist(HIST("h2TPCNSigmaElNegVsPt"), isSignal, negTrack.pt(), negTrack.tpcNSigmaEl());
      if (posTrack.tpcNSigmaEl() < photonDauMinTPCNSigmaEl || posTrack.tpcNSigmaEl() > photonDauMaxTPCNSigmaEl ||
          negTrack.tpcNSigmaEl() < photonDauMinTPCNSigmaEl || negTrack.tpcNSigmaEl() > photonDauMaxTPCNSigmaEl) {
        continue;
      }
      if (posTrack.pt() > photonDauMaxPt || negTrack.pt() > photonDauMaxPt) {
        continue;
      }
      fillPhotonStep(4);

      fillPhotonHist(HIST("hDcaV0Daughters"), isSignal, v0.dcaV0daughters());
      if (v0.dcaV0daughters() > photonMaxDCAV0Dau) {
        continue;
      }
      fillPhotonStep(5);

      std::array<float, 3> secVtx{v0.x(), v0.y(), v0.z()};
      float radius = v0.v0radius();
      if (radius < photonMinRadius || radius > photonMaxRadius) {
        continue;
      }
      fillPhotonStep(6);

      std::array<float, 3> pNeg{v0.pxneg(), v0.pyneg(), v0.pzneg()};
      std::array<float, 3> pPos{v0.pxpos(), v0.pypos(), v0.pzpos()};

      float photonOpeningAngle = o2::aod::pwgem::dilepton::utils::pairutil::getOpeningAngle(pPos[0], pPos[1], pPos[2], pNeg[0], pNeg[1], pNeg[2]);
      fillPhotonHist(HIST("hOpeningAngle"), isSignal, photonOpeningAngle);
      if (photonOpeningAngle > photonMaxOpeningAngle) {
        continue;
      }
      fillPhotonStep(7);

      // difference of the daughter tracks' polar angles
      float photonDeltaTheta = 2.f * std::atan(std::exp(-posTrack.eta())) - 2.f * std::atan(std::exp(-negTrack.eta()));
      fillPhotonHist(HIST("hDeltaTheta"), isSignal, photonDeltaTheta);
      if (std::abs(photonDeltaTheta) > photonMaxDeltaTheta) {
        continue;
      }
      fillPhotonStep(8);

      std::array<float, 3> pGamma{pNeg[0] + pPos[0], pNeg[1] + pPos[1], pNeg[2] + pPos[2]};
      float cosPA = RecoDecay::cpa(pv, secVtx, pGamma);
      if (cosPA < photonMinV0cospa) {
        continue;
      }
      fillPhotonStep(9);

      float photonY = RecoDecay::y(pGamma, o2::constants::physics::MassGamma);
      if (photonY < photonMinRapidity || photonY > photonMaxRapidity) {
        continue;
      }
      fillPhotonStep(10);

      float qtarm = v0.qtarm();
      float alpha = v0.alpha();
      if (qtarm > photonMaxQt) {
        continue;
      }
      fillPhotonStep(11);

      if (std::abs(alpha) > photonMaxAlpha) {
        continue;
      }
      fillPhotonStep(12);

      float mGamma = v0.mGamma();
      if (mGamma > photonMaxMass) {
        continue;
      }
      fillPhotonStep(13);
      fillSigmaCounter(sigmaId, 4);

      fillPhotonHist(HIST("hMass"), isSignal, mGamma);
      fillPhotonHist(HIST("hPt"), isSignal, std::hypot(pGamma[0], pGamma[1]));
      fillPhotonHist(HIST("hRadius"), isSignal, radius);
      fillPhotonHist(HIST("h2ArmenterosPodolanski"), isSignal, alpha, qtarm);
      fillPhotonHist(HIST("h2ConvPointXY"), isSignal, secVtx[0], secVtx[1]);

      photons.push_back({secVtx[0], secVtx[1], secVtx[2], pGamma[0], pGamma[1], pGamma[2], mGamma, alpha, qtarm, radius, v0.dcaV0daughters(), cosPA, negTrack, posTrack});
    }

    return photons;
  }

  static uint64_t photonId(int64_t posTrackId, int64_t negTrackId)
  {
    return (static_cast<uint64_t>(posTrackId) << 32) | static_cast<uint32_t>(negTrackId);
  }

  // psi_pair of the photon daughters (similar to PsiPair in PWGEM)
  template <typename TTrack>
  float photonPsiPair(const TTrack& posTrack, const TTrack& negTrack, float convRadius)
  {
    for (const float& offsetR : {60.f, 30.f, 10.f}) {
      auto posTrackPar = getTrackParCov(posTrack);
      auto negTrackPar = getTrackParCov(negTrack);
      posTrackPar.setPID(o2::track::PID::Electron);
      negTrackPar.setPID(o2::track::PID::Electron);
      if (!o2::base::Propagator::Instance()->propagateToR(posTrackPar, convRadius + offsetR) || !o2::base::Propagator::Instance()->propagateToR(negTrackPar, convRadius + offsetR)) {
        continue;
      }
      std::array<float, 3> pPos{};
      std::array<float, 3> pNeg{};
      posTrackPar.getPxPyPzGlo(pPos);
      negTrackPar.getPxPyPzGlo(pNeg);
      return o2::aod::pwgem::dilepton::utils::pairutil::getPsiPair(pPos[0], pPos[1], pPos[2], pNeg[0], pNeg[1], pNeg[2]);
    }
    return 999.f;
  }

  template <typename TCollisions>
  void markGoodCollisions(const TCollisions& collisions)
  {
    mGoodCollision.assign(collisions.size(), false);
    for (const auto& collision : collisions) {
      if (std::abs(collision.posZ()) > cutZVertex || !collision.sel8()) {
        continue;
      }
      mGoodCollision[collision.globalIndex()] = true;
    }
  }

  // photon daughter selection of the self-pairing (TPC histograms filled before cuts)
  template <typename TTrack>
  bool selectPhotonDaughter(const TTrack& track, bool isSignalLeg)
  {
    if (track.eta() < photonDauEtaMin || track.eta() > photonDauEtaMax) {
      return false;
    }
    if (track.sign() > 0) {
      fillPhotonHist(HIST("h2TPCNClsPosVsPt"), isSignalLeg, track.pt(), track.tpcNClsFound());
    } else {
      fillPhotonHist(HIST("h2TPCNClsNegVsPt"), isSignalLeg, track.pt(), track.tpcNClsFound());
    }
    if (track.tpcNClsFound() < photonDauMinTpcNCls) {
      return false;
    }
    if (track.sign() > 0) {
      fillPhotonHist(HIST("h2TPCNSigmaElPosVsPt"), isSignalLeg, track.pt(), track.tpcNSigmaEl());
    } else {
      fillPhotonHist(HIST("h2TPCNSigmaElNegVsPt"), isSignalLeg, track.pt(), track.tpcNSigmaEl());
    }
    if (track.tpcNSigmaEl() < photonDauMinTPCNSigmaEl || track.tpcNSigmaEl() > photonDauMaxTPCNSigmaEl) {
      return false;
    }
    return track.pt() <= photonDauMaxPt;
  }

  // collisions of a track with only TPC time (all inside time range [t - backward, t + forward] plus a margin)
  template <typename TTrack>
  TrackCand tpcOnlyTrackCand(const TTrack& track, const std::vector<double>& collisionBcNS, const std::vector<double>& collisionTimeNS)
  {
    o2::aod::track::extensions::TPCTimeErrEncoding timeEncoding{};
    timeEncoding.encoding.timeErr = track.trackTimeRes();
    double trackTimeNS = collisionBcNS[track.collisionId()] + track.trackTime(); // trackTime() is relative to the BC of its collision
    double timeMin = trackTimeNS - timeEncoding.getDeltaTBwd() - photonPoolTimeMarginNS;
    double timeMax = trackTimeNS + timeEncoding.getDeltaTFwd() + photonPoolTimeMarginNS;
    int firstCollIdx = track.collisionId();
    int lastCollIdx = firstCollIdx;
    for (int collIdx = 0; collIdx < static_cast<int>(collisionTimeNS.size()); ++collIdx) {
      if (collisionTimeNS[collIdx] >= timeMin && collisionTimeNS[collIdx] <= timeMax) {
        firstCollIdx = std::min(firstCollIdx, collIdx);
        lastCollIdx = std::max(lastCollIdx, collIdx);
      }
    }
    histos.fill(HIST("Photon/Inclusive/hTpcTimeRangeNColl"), lastCollIdx - firstCollIdx + 1);
    return TrackCand{.Idxtr = static_cast<int>(track.globalIndex()), .collBracket = {firstCollIdx, lastCollIdx}};
  }

  // electron-positron pairs with |theta+ - theta-| <= photonMaxDeltaTheta that share at least one collision (theta should be similar)
  std::vector<SVCand> findPairsByPolarAngle(const std::vector<float>& thetaOfTrack, int nCollisions)
  {
    // electrons of each collision sorted by polar angle, so a positron only looks at those inside its theta window
    std::vector<std::vector<std::pair<float, int>>> electronsByCollision(nCollisions);
    for (int iElectron = 0; iElectron < static_cast<int>(mElectronPool.size()); ++iElectron) {
      const auto& electron = mElectronPool[iElectron];
      for (int collIdx = electron.collBracket.getMin(); collIdx <= electron.collBracket.getMax(); ++collIdx) {
        electronsByCollision[collIdx].push_back({thetaOfTrack[electron.Idxtr], iElectron});
      }
    }
    for (size_t collIdx = 0; collIdx < electronsByCollision.size(); ++collIdx) {
      std::sort(electronsByCollision[collIdx].begin(), electronsByCollision[collIdx].end());
    }

    std::vector<SVCand> pairs;
    std::vector<int> lastPositronPaired(mElectronPool.size(), -1); // an electron in several collisions is paired once per positron
    for (int iPositron = 0; iPositron < static_cast<int>(mPositronPool.size()); ++iPositron) {
      const auto& positron = mPositronPool[iPositron];
      float thetaMin = thetaOfTrack[positron.Idxtr] - photonMaxDeltaTheta;
      float thetaMax = thetaOfTrack[positron.Idxtr] + photonMaxDeltaTheta;
      for (int collIdx = positron.collBracket.getMin(); collIdx <= positron.collBracket.getMax(); ++collIdx) {
        const auto& electronsOfCollision = electronsByCollision[collIdx];
        auto electronIt = std::lower_bound(electronsOfCollision.begin(), electronsOfCollision.end(), std::make_pair(thetaMin, -1));
        for (; electronIt != electronsOfCollision.end() && electronIt->first <= thetaMax; ++electronIt) {
          int iElectron = electronIt->second;
          if (lastPositronPaired[iElectron] == iPositron) {
            continue;
          }
          lastPositronPaired[iElectron] = iPositron;
          const auto& electron = mElectronPool[iElectron];
          pairs.push_back(SVCand{.tr0Idx = electron.Idxtr, .tr1Idx = positron.Idxtr, .collBracket = electron.collBracket.getOverlap(positron.collBracket)});
        }
      }
    }
    return pairs;
  }

  // cuts before the fit (similar polar angles, and the circles of the two tracks in xy touching at the conversion point)
  template <typename TTrack>
  bool passesPreFitCuts(const TTrack& posTrack, const TTrack& negTrack, float deltaTheta, bool isSignal)
  {
    fillPhotonHist(HIST("hDeltaTheta"), isSignal, deltaTheta);
    if (std::abs(deltaTheta) > photonMaxDeltaTheta) {
      return false;
    }
    float sna = 0.f;
    float csa = 0.f;
    o2::math_utils::CircleXYf_t posCircle;
    o2::math_utils::CircleXYf_t negCircle;
    getTrackParCov(posTrack).getCircleParams(mBz, posCircle, sna, csa);
    getTrackParCov(negTrack).getCircleParams(mBz, negCircle, sna, csa);
    float circleTouchDist = std::abs(std::hypot(posCircle.xC - negCircle.xC, posCircle.yC - negCircle.yC) - (posCircle.rC + negCircle.rC));
    return circleTouchDist <= photonMaxCircleTouchDist;
  }

  static constexpr int LastDaughterCutStep = 4;

  // photon vertex fitted in one collision
  struct PhotonFit {
    float dca = 0.f;
    std::array<float, 3> secVtx{};
    std::array<float, 3> pPos{};
    std::array<float, 3> pNeg{};
    float posTrackZ = 0.f; // z of the positive daughter at the fit, after moving it to the collision's time
  };

  // moves a track with only a TPC time to the time of the collision
  template <typename TCollisions, typename TCollision, typename TTrack>
  bool moveToCollision(const TCollision& collision, const TTrack& track, o2::track::TrackParCov& trackParCov)
  {
    if (!(track.flags() & o2::aod::track::TrackTimeAsym)) {
      return true;
    }
    return mVDriftMgr.moveTPCTrack<aod::BCsWithTimestamps, TCollisions>(collision, track, trackParCov);
  }

  // vertex fit of the photon daughters, in the collinear mode if a daughter has no ITS
  bool fitPhoton(o2::track::TrackParCov posTrackParCov, o2::track::TrackParCov negTrackParCov, bool collinear, PhotonFit& fit)
  {
    posTrackParCov.setPID(o2::track::PID::Electron);
    negTrackParCov.setPID(o2::track::PID::Electron);

    // photon fit settings, reset afterwards since the fitter is shared with the Sigma+ vertex fit
    fitter.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrLUT);
    fitter.setMaxDXYIni(photonMaxDXYIni);
    fitter.setCollinear(collinear);
    int nCand = 0;
    try {
      nCand = fitter.process(posTrackParCov, negTrackParCov);
    } catch (...) {
      nCand = 0;
    }
    fitter.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrNONE);
    fitter.setMaxDXYIni(1e9);
    fitter.setCollinear(false);
    if (nCand == 0 || !fitter.propagateTracksToVertex()) {
      return false;
    }

    fit.dca = std::sqrt(fitter.getChi2AtPCACandidate());
    fit.secVtx = fitter.getPCACandidatePos();
    fitter.getTrack(0).getPxPyPzGlo(fit.pPos);
    fitter.getTrack(1).getPxPyPzGlo(fit.pNeg);
    fit.posTrackZ = posTrackParCov.getZ();
    return true;
  }

  // a photon is kept in every compatible collision where it passes the cuts, and the proton of that collision then decides
  template <bool IsMC, typename TCollisions, typename TTracks>
  std::vector<std::vector<PhotonCand<typename TTracks::iterator>>> findPhotonsSelfPaired(const TCollisions& collisions, const TTracks& tracks, aod::AmbiguousTracks const& ambiguousTracks, aod::BCsWithTimestamps const& bcs)
  {
    std::vector<std::vector<PhotonCand<typename TTracks::iterator>>> photonsByCollision(collisions.size());

    svPhotonPoolCreator.clearPools();
    svPhotonPoolCreator.fillBC2Coll(collisions, bcs);
    mElectronPool.clear();
    mPositronPool.clear();

    // BC start and time of every collision, in ns relative to the BC of the first collision
    std::vector<double> collisionBcNS(collisions.size());
    std::vector<double> collisionTimeNS(collisions.size());
    uint64_t firstGlobalBC = collisions.begin().template bc_as<aod::BCsWithTimestamps>().globalBC();
    for (const auto& collision : collisions) {
      int64_t bcDiff = static_cast<int64_t>(collision.template bc_as<aod::BCsWithTimestamps>().globalBC()) - static_cast<int64_t>(firstGlobalBC);
      collisionBcNS[collision.globalIndex()] = bcDiff * o2::constants::lhc::LHCBunchSpacingNS;
      collisionTimeNS[collision.globalIndex()] = collisionBcNS[collision.globalIndex()] + collision.collisionTime();
    }

    // select the daughter tracks and find the collisions each is compatible with
    std::vector<float> thetaOfTrack(tracks.size(), 0.f);
    std::vector<int> mcPhotonOfTrack(tracks.size(), -1); // MC: Sigma+ photon the track comes from
    std::vector<int> mcSigmaOfTrack(tracks.size(), -1);  // MC: Sigma+ the track comes from
    for (const auto& track : tracks) {
      int mcPhoton = -1;
      int mcSigma = -1;
      if constexpr (IsMC) {
        if (track.has_mcParticle()) {
          mcSigma = findSigmaPlusAncestorOfPhotonDaughter(track.template mcParticle_as<aod::McParticles>(), mcPhoton);
        }
      }
      if (!selectPhotonDaughter(track, mcSigma >= 0)) {
        continue;
      }
      thetaOfTrack[track.globalIndex()] = 2.f * std::atan(std::exp(-track.eta()));
      mcPhotonOfTrack[track.globalIndex()] = mcPhoton;
      mcSigmaOfTrack[track.globalIndex()] = mcSigma;

      if (track.flags() & o2::aod::track::TrackTimeAsym) {
        // track with only a TPC time
        if (!track.has_collision()) {
          continue;
        }
        if (track.sign() < 0) {
          mElectronPool.push_back(tpcOnlyTrackCand(track, collisionBcNS, collisionTimeNS));
        } else {
          mPositronPool.push_back(tpcOnlyTrackCand(track, collisionBcNS, collisionTimeNS));
        }
      } else {
        // track with an ITS time (matched to collisions by svPoolCreator)
        svPhotonPoolCreator.appendTrackCand(track, collisions, track.sign() < 0 ? PDG_t::kElectron : PDG_t::kPositron, ambiguousTracks, bcs);
      }
    }

    // svPoolCreator pools are ordered dau0 pos, dau0 neg, dau1 pos, dau1 neg: electrons are in 1, positrons are in 2
    auto timeMatchedPools = svPhotonPoolCreator.getTrackCandPool();
    mElectronPool.insert(mElectronPool.end(), timeMatchedPools[1].begin(), timeMatchedPools[1].end());
    mPositronPool.insert(mPositronPool.end(), timeMatchedPools[2].begin(), timeMatchedPools[2].end());

    auto pairs = findPairsByPolarAngle(thetaOfTrack, collisions.size());

    // fit and select each pair in every good collision it is compatible with
    uint64_t nTruePairs = 0;
    for (const auto& pair : pairs) {
      int electronIdx = pair.tr0Idx;
      int positronIdx = pair.tr1Idx;
      int sigmaId = -1;
      if constexpr (IsMC) {
        if (mcPhotonOfTrack[electronIdx] >= 0 && mcPhotonOfTrack[electronIdx] == mcPhotonOfTrack[positronIdx]) {
          sigmaId = mcSigmaOfTrack[electronIdx];
          nTruePairs++;
        }
      }
      bool isSignal = sigmaId >= 0;

      auto negTrack = tracks.rawIteratorAt(electronIdx);
      auto posTrack = tracks.rawIteratorAt(positronIdx);
      if (!passesPreFitCuts(posTrack, negTrack, thetaOfTrack[positronIdx] - thetaOfTrack[electronIdx], isSignal)) {
        continue;
      }

      // both daughters with only a TPC time on the same TPC side: moving to another collision shifts both by the same z,
      // so the pair is fitted once and the vertex is only shifted in z for the other collisions
      bool posTpcTime = posTrack.flags() & o2::aod::track::TrackTimeAsym;
      bool negTpcTime = negTrack.flags() & o2::aod::track::TrackTimeAsym;
      bool sameShiftForAllCollisions = posTpcTime && negTpcTime && ((posTrack.tgl() > 0.f) == (negTrack.tgl() > 0.f));
      bool collinearFit = !posTrack.hasITS() || !negTrack.hasITS();
      std::optional<PhotonFit> referenceFit;

      int maxStepReached = LastDaughterCutStep; // furthest selection step reached in any collision (for selection counter histograms)
      auto passStep = [&](int step) { maxStepReached = std::max(maxStepReached, step); };
      bool anyFitOk = false;
      float smallestDca = 1e9f;
      bool openingAngleFilled = false;
      int nCopies = 0;

      for (int collIdx = pair.collBracket.getMin(); collIdx <= pair.collBracket.getMax(); ++collIdx) {
        if (!mGoodCollision[collIdx]) {
          continue;
        }
        auto collision = collisions.rawIteratorAt(collIdx);

        PhotonFit fit;
        if (sameShiftForAllCollisions && referenceFit) {
          auto posTrackParCov = getTrackParCov(posTrack);
          if (!moveToCollision<TCollisions>(collision, posTrack, posTrackParCov)) {
            continue;
          }
          fit = *referenceFit;
          fit.secVtx[2] += posTrackParCov.getZ() - referenceFit->posTrackZ;
        } else {
          auto posTrackParCov = getTrackParCov(posTrack);
          auto negTrackParCov = getTrackParCov(negTrack);
          if (!moveToCollision<TCollisions>(collision, posTrack, posTrackParCov) || !moveToCollision<TCollisions>(collision, negTrack, negTrackParCov)) {
            continue;
          }
          bool fitOk = fitPhoton(posTrackParCov, negTrackParCov, collinearFit, fit);
          if (fitOk) {
            anyFitOk = true;
            smallestDca = std::min(smallestDca, fit.dca);
          }
          if (!fitOk || fit.dca > photonMaxDCAV0Dau) {
            if (sameShiftForAllCollisions) {
              break; // same result in every collision
            }
            continue;
          }
          if (sameShiftForAllCollisions) {
            referenceFit = fit;
          }
        }
        passStep(5);

        float radius = std::hypot(fit.secVtx[0], fit.secVtx[1]);
        if (radius < photonMinRadius || radius > photonMaxRadius) {
          continue;
        }
        passStep(6);

        float photonOpeningAngle = o2::aod::pwgem::dilepton::utils::pairutil::getOpeningAngle(fit.pPos[0], fit.pPos[1], fit.pPos[2], fit.pNeg[0], fit.pNeg[1], fit.pNeg[2]);
        if (!openingAngleFilled) { // once per pair
          openingAngleFilled = true;
          fillPhotonHist(HIST("hOpeningAngle"), isSignal, photonOpeningAngle);
        }
        if (photonOpeningAngle > photonMaxOpeningAngle) {
          continue;
        }
        passStep(8); // opening angle, and delta theta that is cut before the fit

        std::array<float, 3> pGamma{fit.pNeg[0] + fit.pPos[0], fit.pNeg[1] + fit.pPos[1], fit.pNeg[2] + fit.pPos[2]};
        float cosPA = RecoDecay::cpa(std::array{collision.posX(), collision.posY(), collision.posZ()}, fit.secVtx, pGamma);
        if (cosPA < photonMinV0cospa) {
          continue;
        }
        passStep(9);

        float photonY = RecoDecay::y(pGamma, o2::constants::physics::MassGamma);
        if (photonY < photonMinRapidity || photonY > photonMaxRapidity) {
          continue;
        }
        passStep(10);

        float qtarm = v0Qt(fit.pPos[0], fit.pPos[1], fit.pPos[2], fit.pNeg[0], fit.pNeg[1], fit.pNeg[2]);
        float alpha = v0Alpha(fit.pPos[0], fit.pPos[1], fit.pPos[2], fit.pNeg[0], fit.pNeg[1], fit.pNeg[2]);
        if (qtarm > photonMaxQt) {
          continue;
        }
        passStep(11);

        if (std::abs(alpha) > photonMaxAlpha) {
          continue;
        }
        passStep(12);

        float mGamma = RecoDecay::m(std::array{fit.pPos, fit.pNeg}, std::array{o2::constants::physics::MassElectron, o2::constants::physics::MassElectron});
        if (mGamma > photonMaxMass) {
          continue;
        }
        passStep(13);

        fillPhotonHist(HIST("hMass"), isSignal, mGamma);
        fillPhotonHist(HIST("hPt"), isSignal, std::hypot(pGamma[0], pGamma[1]));
        fillPhotonHist(HIST("hRadius"), isSignal, radius);
        fillPhotonHist(HIST("h2ArmenterosPodolanski"), isSignal, alpha, qtarm);
        fillPhotonHist(HIST("h2ConvPointXY"), isSignal, fit.secVtx[0], fit.secVtx[1]);

        nCopies++;
        photonsByCollision[collIdx].push_back({fit.secVtx[0], fit.secVtx[1], fit.secVtx[2], pGamma[0], pGamma[1], pGamma[2], mGamma, alpha, qtarm, radius, fit.dca, cosPA, negTrack, posTrack});
      }

      if (anyFitOk) {
        fillPhotonHist(HIST("hDcaV0Daughters"), isSignal, smallestDca);
      }
      if (nCopies > 0) {
        fillPhotonHist(HIST("hCollisionsPerPhoton"), isSignal, nCopies);
        fillSigmaCounter(sigmaId, 4);
      }
      for (int step = 5; step <= maxStepReached; ++step) {
        fillPhotonHist(HIST("hSelectionCounter"), isSignal, step);
      }
    }

    // steps 0-4 are the daughter-track cuts applied before the pairing, so they count all pairs
    for (int step = 0; step <= LastDaughterCutStep; ++step) {
      histos.fill(HIST("Photon/Inclusive/hSelectionCounter"), step, static_cast<double>(pairs.size()));
      if constexpr (IsMC) {
        histos.fill(HIST("Photon/True/hSelectionCounter"), step, static_cast<double>(nTruePairs));
      }
    }

    return photonsByCollision;
  }

  // proton candidate selection
  template <bool IsMC, typename TTrack>
  bool selectProton(const TTrack& track)
  {
    bool isSignal = false;
    if constexpr (IsMC) {
      if (track.has_mcParticle()) {
        isSignal = findSigmaPlusMotherOfProton(track.template mcParticle_as<aod::McParticles>()) >= 0;
      }
    }
    auto fillProtonStep = [&](int step) {
      fillProtonHist(HIST("hSelectionCounter"), isSignal, step);
    };
    fillProtonStep(0);

    if (track.pt() < protonMinPt) {
      return false;
    }
    fillProtonStep(1);

    if (std::abs(track.eta()) > protonMaxEta) {
      return false;
    }
    fillProtonStep(2);

    fillProtonHist(HIST("h2TPCNClsVsPt"), isSignal, track.pt(), track.tpcNClsFound());
    if (track.tpcNClsFound() < protonMinTpcNCls) {
      return false;
    }
    fillProtonStep(3);

    if (std::abs(track.tpcNSigmaPr()) > protonMaxTPCNSigma) {
      return false;
    }
    fillProtonStep(4);

    if (track.pt() > protonPtMinRequireTOF) {
      if (track.hasTOF()) {
        if (std::abs(track.tofNSigmaPr()) > protonMaxTOFNSigma) {
          return false;
        }
      } else if (protonRequireTofHit) {
        return false;
      }
    }
    fillProtonStep(5);

    fillProtonHist(HIST("h2DcaToPVVsPt"), isSignal, track.pt(), std::abs(track.dcaXY()));
    if (std::abs(track.dcaXY()) < protonMinDcaToPV || std::abs(track.dcaXY()) > protonMaxDcaToPV) {
      return false;
    }
    fillProtonStep(6);

    fillProtonHist(HIST("hPt"), isSignal, track.pt());
    fillProtonHist(HIST("h2TPCNSigmaVsPt"), isSignal, track.pt(), track.tpcNSigmaPr());
    if (track.hasTOF()) {
      fillProtonHist(HIST("h2TOFNSigmaVsPt"), isSignal, track.pt(), track.tofNSigmaPr());
    }

    return true;
  }

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

  // MC: for an e+ or e- from Sigma+ -> pi0 -> photon -> e+e-, returns the Sigma+ index and sets photonIndex (-1 otherwise)
  template <typename TMcPart>
  int findSigmaPlusAncestorOfPhotonDaughter(const TMcPart& mcDaughter, int& photonIndex)
  {
    auto const& mothers = mcDaughter.template mothers_as<aod::McParticles>();
    if (mothers.empty() || mothers.front().pdgCode() != PDG_t::kGamma) {
      return -1;
    }
    auto const& pi0Mothers = mothers.front().template mothers_as<aod::McParticles>();
    if (pi0Mothers.empty() || std::abs(pi0Mothers.front().pdgCode()) != PDG_t::kPi0) {
      return -1;
    }
    auto const& sigmaMothers = pi0Mothers.front().template mothers_as<aod::McParticles>();
    if (sigmaMothers.empty() || std::abs(sigmaMothers.front().pdgCode()) != PDG_t::kSigmaPlus) {
      return -1;
    }
    photonIndex = mothers.front().globalIndex();
    return sigmaMothers.front().globalIndex();
  }

  template <typename TMcPart>
  int findSigmaPlusMotherOfPhoton(const TMcPart& mcPos, const TMcPart& mcNeg)
  {
    int posPhoton = -1;
    int negPhoton = -1;
    int sigmaId = findSigmaPlusAncestorOfPhotonDaughter(mcPos, posPhoton);
    if (sigmaId < 0 || findSigmaPlusAncestorOfPhotonDaughter(mcNeg, negPhoton) < 0 || posPhoton != negPhoton) {
      return -1;
    }
    return sigmaId;
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
  template <bool IsMC, typename TTrack, typename TPhotonTrack>
  void buildSigmaPlusCandidate(const TTrack& protonTrack, const PhotonCand<TPhotonTrack>& photon, const std::array<float, 3>& pv)
  {
    auto posTrack = photon.posTrack;
    auto negTrack = photon.negTrack;
    int photonLegsWithoutIts = (posTrack.hasITS() ? 0 : 1) + (negTrack.hasITS() ? 0 : 1);

    bool isSignal = false;
    bool protonIsSignal = false;
    bool photonIsSignal = false;
    bool collisionIdCheck = false; // MC only: true if the proton's MC collision matches the reconstructed collision
    std::array<float, 3> mcTrueVtx{};
    std::array<float, 3> mcTrueMomProton{};
    std::array<float, 3> mcTrueMomGamma{};
    std::array<float, 3> mcTrueMomSigmaPlus{};
    float decayRadiusMC = -999.f;
    float massMC = -999.f;
    int matchedSigmaId = -1;
    if constexpr (IsMC) {
      if (protonTrack.has_mcParticle()) {
        auto mcProton = protonTrack.template mcParticle_as<aod::McParticles>();

        auto protonCollision = protonTrack.template collision_as<CollisionsFullMC>();
        if (protonCollision.has_mcCollision()) {
          collisionIdCheck = protonCollision.mcCollision().globalIndex() == mcProton.mcCollisionId();
        }

        int protonSigmaIdx = findSigmaPlusMotherOfProton(mcProton);
        protonIsSignal = protonSigmaIdx >= 0;

        if (posTrack.has_mcParticle() && negTrack.has_mcParticle()) {
          auto mcPos = posTrack.template mcParticle_as<aod::McParticles>();
          auto mcNeg = negTrack.template mcParticle_as<aod::McParticles>();

          auto const& posMothers = mcPos.template mothers_as<aod::McParticles>();
          if (!posMothers.empty()) {
            auto mcGamma = posMothers.front();
            mcTrueMomGamma = {mcGamma.px(), mcGamma.py(), mcGamma.pz()};
          }

          int photonSigmaIdx = findSigmaPlusMotherOfPhoton(mcPos, mcNeg);
          photonIsSignal = photonSigmaIdx >= 0;
          isSignal = protonIsSignal && photonIsSignal && photonSigmaIdx == protonSigmaIdx;
        }

        if (isSignal) {
          auto mcSigmaPlusMother = mcProton.template mothers_as<aod::McParticles>().front();
          mcTrueVtx = {mcProton.vx(), mcProton.vy(), mcProton.vz()};
          mcTrueMomProton = {mcProton.px(), mcProton.py(), mcProton.pz()};
          mcTrueMomSigmaPlus = {mcSigmaPlusMother.px(), mcSigmaPlusMother.py(), mcSigmaPlusMother.pz()};
          matchedSigmaId = mcSigmaPlusMother.globalIndex();
          decayRadiusMC = std::hypot(mcTrueVtx[0] - mcSigmaPlusMother.vx(), mcTrueVtx[1] - mcSigmaPlusMother.vy());
          massMC = std::sqrt(mcSigmaPlusMother.e() * mcSigmaPlusMother.e() - mcSigmaPlusMother.p() * mcSigmaPlusMother.p());
        }
      }
    }

    auto fillCandStep = [&](int step) {
      fillCandHist(HIST("hSelectionCounter"), isSignal, step);
      fillCandHist(HIST("h2SelectionCounterVsLegType"), isSignal, step, photonLegsWithoutIts);
    };
    fillCandStep(0); // all pairs

    // reject the pair if the proton track is one of the photon's own daughter tracks
    if (protonTrack.globalIndex() == posTrack.globalIndex() || protonTrack.globalIndex() == negTrack.globalIndex()) {
      return;
    }
    fillCandStep(1); // autocorrelation

    auto protonTrackParCov = getTrackParCov(protonTrack);
    std::array<float, 3> protonOrigPos{};
    protonTrackParCov.getXYZGlo(protonOrigPos);

    std::array<float, 21> zeroCov{};
    auto photonTrackParCov = o2::track::TrackParCov({photon.x, photon.y, photon.z}, {photon.px, photon.py, photon.pz}, zeroCov, 0, true);
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
    fillCandStep(2); // vertex fit
    std::array<float, 3> protonAtPca{};
    std::array<float, 3> photonAtPca{};
    fitter.getTrack(0).getXYZGlo(protonAtPca);
    fitter.getTrack(1).getXYZGlo(photonAtPca);

    float dcaProtonGamma = std::sqrt(fitter.getChi2AtPCACandidate());
    fillCandHist(HIST("hDcaProtonGamma"), isSignal, dcaProtonGamma);
    if (dcaProtonGamma > candMaxDcaProtonGamma) {
      return;
    }
    fillCandStep(3); // DCA(p,gamma)

    // decay vertex between the proton's and the photon's point of closest approach (relative weight given by candVertexProtonWeight)
    std::array<float, 3> secVtx{};
    for (size_t i = 0; i < secVtx.size(); ++i) {
      secVtx[i] = candVertexProtonWeight * protonAtPca[i] + (1.f - candVertexProtonWeight) * photonAtPca[i];
    }
    float radius = std::hypot(secVtx[0], secVtx[1]);
    fillCandHist(HIST("hRadius"), isSignal, radius);
    if (radius < candMinRadius || radius > candMaxRadius) {
      return;
    }
    fillCandStep(4); // radius

    std::array<float, 3> flightVec{secVtx[0] - pv[0], secVtx[1] - pv[1], secVtx[2] - pv[2]};
    std::array<float, 3> nHatOrig = normalize3(flightVec);
    float flightDistance = std::sqrt(dot3(flightVec, flightVec));

    fillCandHist(HIST("hFlightDistance"), isSignal, flightDistance);
    if (flightDistance < candMinFlightDistance || flightDistance > candMaxFlightDistance) {
      return;
    }
    fillCandStep(5); // flight distance

    std::array<float, 3> pProton{};
    std::array<float, 3> pGamma1{};
    fitter.getTrack(0).getPxPyPzGlo(pProton);
    fitter.getTrack(1).getPxPyPzGlo(pGamma1);

    float e1 = std::sqrt(dot3(pGamma1, pGamma1));
    float massPi0 = o2::constants::physics::MassPionNeutral;
    constexpr float CoefADegenerateThreshold = 1e-6f;

    std::array<float, 3> nHat{};
    std::array<float, 3> eOut{};
    std::array<float, 3> eIn{};
    float coefA = 0.f, coefB = 0.f;
    float pGamma2In = 0.f, pGamma2Out = 0.f, tPerp2 = 0.f;

    auto solveForFlightDir = [&](const std::array<float, 3>& nHatUse) {
      nHat = nHatUse;
      eOut = normalize3(cross3(pProton, nHat)); // normal to the decay plane
      eIn = cross3(eOut, nHat);                 // in the decay plane, transverse to n

      float pProtonIn = dot3(pProton, eIn);
      float pGamma1Long = dot3(pGamma1, nHat);
      float pGamma1In = dot3(pGamma1, eIn);
      float pGamma1Out = dot3(pGamma1, eOut);

      pGamma2In = -(pProtonIn + pGamma1In);
      pGamma2Out = -pGamma1Out;
      tPerp2 = pGamma2In * pGamma2In + pGamma2Out * pGamma2Out;

      // m_pi0^2 = 2*(|p_gamma1||p_gamma2| - p_gamma1.p_gamma2) as coefA*x^2 + coefB*x + coefC = 0
      float coefK = massPi0 * massPi0 + 2.f * (pGamma1In * pGamma2In + pGamma1Out * pGamma2Out);
      float coefL = 2.f * pGamma1Long;
      coefA = 4.f * e1 * e1 - coefL * coefL;
      if (std::abs(coefA) < CoefADegenerateThreshold) {
        return -1.f;
      }
      coefB = -2.f * coefK * coefL;
      float coefC = 4.f * e1 * e1 * tPerp2 - coefK * coefK;
      return coefB * coefB - 4.f * coefA * coefC;
    };

    // if the measured flight direction gives no real solution, flight directions on rings of increasing
    // tilt around it are tried (the first ring with a real solution is taken, and on it the direction with the largest discriminant)
    float tiltAngle = 0.f;
    float discriminant = solveForFlightDir(nHatOrig);
    if (discriminant < 0.f) {
      constexpr float MaxAbsNzForZHelperAxis = 0.9f;
      std::array<float, 3> helperAxis = std::abs(nHatOrig[2]) < MaxAbsNzForZHelperAxis ? std::array<float, 3>{0.f, 0.f, 1.f} : std::array<float, 3>{1.f, 0.f, 0.f};
      std::array<float, 3> tiltAxisU = normalize3(cross3(nHatOrig, helperAxis));
      std::array<float, 3> tiltAxisV = cross3(nHatOrig, tiltAxisU);
      int nRings = static_cast<int>(std::lround(candMaxTilt / candTiltStep));
      for (int ring = 1; ring <= nRings && discriminant < 0.f; ++ring) {
        float tilt = ring * candTiltStep;
        float bestRingDisc = -1.f;
        std::array<float, 3> bestRingDir{};
        for (int iAzimuth = 0; iAzimuth < candTiltNAzimuth; ++iAzimuth) {
          float azimuth = o2::constants::math::TwoPI * iAzimuth / candTiltNAzimuth;
          std::array<float, 3> tiltDir{};
          for (size_t i = 0; i < tiltDir.size(); ++i) {
            tiltDir[i] = std::cos(tilt) * nHatOrig[i] + std::sin(tilt) * (std::cos(azimuth) * tiltAxisU[i] + std::sin(azimuth) * tiltAxisV[i]);
          }
          float disc = solveForFlightDir(tiltDir);
          if (disc >= 0.f && disc > bestRingDisc) {
            bestRingDisc = disc;
            bestRingDir = tiltDir;
          }
        }
        if (bestRingDisc >= 0.f) {
          discriminant = solveForFlightDir(bestRingDir);
          tiltAngle = tilt;
        }
      }
    }

    float tiltForHist = discriminant >= 0.f ? tiltAngle : 0.4999f; // no real solution
    fillCandHist(HIST("hTiltAngle"), isSignal, tiltForHist);
    if (discriminant < 0.f) {
      return;
    }
    fillCandStep(6); // real root

    // of the two roots, keep the one giving the mass closest to the nominal Sigma+ mass
    float sqrtDisc = std::sqrt(discriminant);
    std::array<float, 2> roots{(-coefB + sqrtDisc) / (2.f * coefA), (-coefB - sqrtDisc) / (2.f * coefA)};

    float rootCenter = -coefB / (2.f * coefA);
    fillCandHist(HIST("hRootCenter"), isSignal, rootCenter);
    if (candRejectNegRootCenter && rootCenter < 0.f) {
      return;
    }
    if (rootCenter > candMaxRootCenter) {
      return;
    }

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
    fillCandStep(7); // valid root

    std::array<float, 3> pSigma{pProton[0] + pGamma1[0] + bestMomGamma2[0], pProton[1] + pGamma1[1] + bestMomGamma2[1], pProton[2] + pGamma1[2] + bestMomGamma2[2]};
    float ptSigma = std::hypot(pSigma[0], pSigma[1]);

    o2::track::TrackPar sigmaTrackPar({secVtx[0], secVtx[1], secVtx[2]}, {pSigma[0], pSigma[1], pSigma[2]}, protonTrack.sign(), true);
    std::array<float, 2> dcaSigmaToPv{};
    o2::base::Propagator::Instance()->propagateToDCA(o2::math_utils::Point3D<float>{pv[0], pv[1], pv[2]}, sigmaTrackPar, mBz, 2.f,
                                                     o2::base::Propagator::MatCorrType::USEMatCorrNONE, &dcaSigmaToPv);
    float candDcaToPV = std::abs(dcaSigmaToPv[0]);
    fillCandHist(HIST("hDcaToPV"), isSignal, candDcaToPV);
    if (candDcaToPV > candMaxDcaToPV) {
      return;
    }
    fillCandStep(8); // DCA to PV

    // angle between the proton momentum and the flight direction, both rotated back to undo the bending
    // of the proton between its reference point and the decay vertex
    std::array<float, 3> protonPath{secVtx[0] - protonOrigPos[0], secVtx[1] - protonOrigPos[1], secVtx[2] - protonOrigPos[2]};
    float protonPathLength = std::sqrt(dot3(protonPath, protonPath));
    float svRadiusFromPv = std::hypot(secVtx[0] - pv[0], secVtx[1] - pv[1]);
    float protonOrigRadiusFromPv = std::hypot(protonOrigPos[0] - pv[0], protonOrigPos[1] - pv[1]);
    float propDir = (svRadiusFromPv < protonOrigRadiusFromPv) ? -1.f : 1.f;
    float qProton = (protonTrack.sign() < 0) ? -1.f : 1.f;
    float protonP = std::sqrt(dot3(pProton, pProton));
    float rCurve = protonP * 1000.f / (0.2998f * std::abs(mBz));
    float bzSign = (mBz < 0) ? -1.f : 1.f;
    float alphaRot = -propDir * qProton * bzSign * o2::constants::math::PI;
    if (protonPathLength / (2.f * rCurve) < 1.f) {
      alphaRot = -2.f * propDir * qProton * bzSign * std::asin(protonPathLength / (2.f * rCurve));
    }
    std::array<float, 3> pProtonRot{
      pProton[0] * std::cos(alphaRot) - pProton[1] * std::sin(alphaRot),
      pProton[0] * std::sin(alphaRot) + pProton[1] * std::cos(alphaRot),
      pProton[2]};
    float alphaVtxRot = -propDir * qProton * bzSign * o2::constants::math::PI * 0.5f;
    if (flightDistance / (2.f * rCurve) < 1.f) {
      alphaVtxRot = -propDir * qProton * bzSign * std::asin(flightDistance / (2.f * rCurve));
    }
    std::array<float, 3> sigmaVertexVecRot{
      flightVec[0] * std::cos(alphaVtxRot) - flightVec[1] * std::sin(alphaVtxRot),
      flightVec[0] * std::sin(alphaVtxRot) + flightVec[1] * std::cos(alphaVtxRot),
      flightVec[2]};
    float antiSigmaPointingAngle = std::acos(std::clamp(dot3(pProtonRot, sigmaVertexVecRot) / std::sqrt(dot3(pProtonRot, pProtonRot) * dot3(sigmaVertexVecRot, sigmaVertexVecRot)), -1.f, 1.f));
    fillCandHist(HIST("hAntiSigmaPointingAngle"), isSignal, antiSigmaPointingAngle);
    if (antiSigmaPointingAngle < candMinAntiSigmaPointingAngle || antiSigmaPointingAngle > candMaxAntiSigmaPointingAngle) {
      return;
    }

    fillCandHist(HIST("hMassSigmaPlus"), isSignal, bestMass);
    fillCandHist(HIST("h2MassVsPt"), isSignal, ptSigma, bestMass);
    fillCandHist(HIST("h2MassVsTilt"), isSignal, tiltAngle, bestMass);
    if (bestMass > candMaxSigmaMass) {
      return;
    }
    fillCandStep(9); // mass

    float candRapidity = RecoDecay::y(pSigma, o2::constants::physics::MassSigmaPlus);
    fillCandHist(HIST("hRapidity"), isSignal, candRapidity);
    if (std::abs(candRapidity) > candMaxRapidity) {
      return;
    }
    fillCandStep(10); // rapidity

    // photon opening angle from the daughters' track momenta
    float photonOpeningAngle = o2::aod::pwgem::dilepton::utils::pairutil::getOpeningAngle(posTrack.px(), posTrack.py(), posTrack.pz(), negTrack.px(), negTrack.py(), negTrack.pz());
    fillCandHist(HIST("hPhotonOpeningAngle"), isSignal, photonOpeningAngle);
    if (photonOpeningAngle > candMaxPhotonOpeningAngle) {
      return;
    }
    fillCandStep(11); // photon opening angle

    // angle between the photon momentum and the line from the decay vertex to the conversion point
    std::array<float, 3> decVtxToConv{photon.x - secVtx[0], photon.y - secVtx[1], photon.z - secVtx[2]};
    float photonPointingAngle = std::acos(std::clamp(dot3(pGamma1, decVtxToConv) / std::sqrt(dot3(pGamma1, pGamma1) * dot3(decVtxToConv, decVtxToConv)), -1.f, 1.f));
    fillCandHist(HIST("hPhotonPointingAngle"), isSignal, photonPointingAngle);
    if (photonPointingAngle > candMaxPhotonPointingAngle) {
      return;
    }
    fillCandStep(12); // photon pointing angle

    std::array<float, 3> photonDir = normalize3({photon.px, photon.py, photon.pz});
    std::array<float, 3> pvToConv{pv[0] - photon.x, pv[1] - photon.y, pv[2] - photon.z};
    std::array<float, 3> pvToConvCrossDir = cross3(pvToConv, photonDir);
    float photonDcaToPV = std::sqrt(dot3(pvToConvCrossDir, pvToConvCrossDir));
    fillCandHist(HIST("hPhotonDcaToPV"), isSignal, photonDcaToPV);
    if (photonDcaToPV > candMaxPhotonDcaToPV) {
      return;
    }
    fillCandStep(13); // photon DCA to PV

    float photonPt2 = photon.px * photon.px + photon.py * photon.py;
    float tAtDcaXY = photonPt2 > 0.f ? ((pv[0] - photon.x) * photon.px + (pv[1] - photon.y) * photon.py) / photonPt2 : 0.f;
    float photonDcaXYToPV = std::hypot(photon.x + tAtDcaXY * photon.px - pv[0], photon.y + tAtDcaXY * photon.py - pv[1]);
    float photonDcaZToPV = photon.z + tAtDcaXY * photon.pz - pv[2];
    float psiPair = photonPsiPair(posTrack, negTrack, photon.radius);

    SigmaPlusCandidate cand;
    cand.photonId = photonId(posTrack.globalIndex(), negTrack.globalIndex());
    cand.isSignal = isSignal;
    cand.photonLegsWithoutIts = photonLegsWithoutIts;
    cand.matchedSigmaId = matchedSigmaId;
    cand.decVtx = secVtx;
    cand.radius = radius;
    cand.flightDistance = flightDistance;
    cand.dcaProtonGamma = dcaProtonGamma;
    cand.dcaToPV = candDcaToPV;
    cand.tiltAngle = tiltAngle;
    cand.rootCenter = rootCenter;
    cand.antiSigmaPointingAngle = antiSigmaPointingAngle;
    cand.pProton = pProton;
    cand.pGamma1 = pGamma1;
    cand.pGamma2 = bestMomGamma2;
    cand.protonTpcNSigma = protonTrack.tpcNSigmaPr();
    cand.protonTofNSigma = protonTrack.tofNSigmaPr();
    cand.protonSign = protonTrack.sign();
    cand.protonItsNCls = protonTrack.itsNCls();
    cand.protonTpcNCls = protonTrack.tpcNClsFound();
    cand.protonDcaXY = protonTrack.dcaXY();
    cand.protonDcaZ = protonTrack.dcaZ();
    cand.photonMass = photon.mGamma;
    cand.photonAlpha = photon.alpha;
    cand.photonQt = photon.qtarm;
    cand.photonRadius = photon.radius;
    cand.photonOpeningAngle = photonOpeningAngle;
    cand.photonPointingAngle = photonPointingAngle;
    cand.photonDcaDau = photon.dcaDau;
    cand.photonCosPAToPV = photon.cosPAToPV;
    cand.photonDcaXYToPV = photonDcaXYToPV;
    cand.photonDcaZToPV = photonDcaZToPV;
    cand.photonPsiPair = psiPair;
    cand.posTpcNSigmaEl = posTrack.tpcNSigmaEl();
    cand.negTpcNSigmaEl = negTrack.tpcNSigmaEl();
    cand.posItsNCls = posTrack.itsNCls();
    cand.negItsNCls = negTrack.itsNCls();
    cand.posTpcNCls = posTrack.tpcNClsFound();
    cand.negTpcNCls = negTrack.tpcNClsFound();
    cand.posTpcNClsFindable = posTrack.tpcNClsFindable();
    cand.negTpcNClsFindable = negTrack.tpcNClsFindable();
    cand.posTpcChi2NCl = posTrack.tpcChi2NCl();
    cand.negTpcChi2NCl = negTrack.tpcChi2NCl();
    cand.posDcaXY = posTrack.dcaXY();
    cand.posDcaZ = posTrack.dcaZ();
    cand.negDcaXY = negTrack.dcaXY();
    cand.negDcaZ = negTrack.dcaZ();
    cand.collisionIdCheck = collisionIdCheck;
    cand.protonIsSignal = protonIsSignal;
    cand.photonIsSignal = photonIsSignal;
    cand.decVtxMC = mcTrueVtx;
    cand.pProtonMC = mcTrueMomProton;
    cand.pGammaMC = mcTrueMomGamma;
    cand.pSigmaPlusMC = mcTrueMomSigmaPlus;
    cand.decayRadiusMC = decayRadiusMC;
    cand.massMC = massMC;
    mCandidatesOfTimeframe.push_back(cand);
  }

  template <bool IsMC>
  void fillCandidateTables(const SigmaPlusCandidate& cand)
  {
    if constexpr (IsMC) {
      if (fillSlimTables) {
        slimSigmaPlusCandsMC(cand.radius,
                             cand.dcaToPV, cand.tiltAngle, cand.dcaProtonGamma,
                             cand.protonSign,
                             cand.protonDcaXY, cand.protonDcaZ,
                             cand.pProton[0], cand.pProton[1], cand.pProton[2],
                             cand.pGamma1[0], cand.pGamma1[1], cand.pGamma1[2],
                             cand.pGamma2[0], cand.pGamma2[1], cand.pGamma2[2],
                             cand.protonTpcNSigma, cand.protonTofNSigma,
                             cand.posTpcNSigmaEl, cand.negTpcNSigmaEl,
                             cand.photonMass,
                             cand.photonDcaDau, cand.photonCosPAToPV, cand.photonDcaXYToPV, cand.photonDcaZToPV, cand.photonPsiPair,
                             cand.posDcaXY, cand.posDcaZ, cand.negDcaXY, cand.negDcaZ,
                             cand.posTpcNClsFindable, cand.negTpcNClsFindable, cand.posTpcChi2NCl, cand.negTpcChi2NCl,
                             cand.collisionIdCheck,
                             cand.isSignal, cand.protonIsSignal, cand.photonIsSignal,
                             cand.decayRadiusMC, cand.massMC,
                             cand.pSigmaPlusMC[0], cand.pSigmaPlusMC[1], cand.pSigmaPlusMC[2]);
        return;
      }
      sigmaPlusCandsMC(cand.decVtx[0], cand.decVtx[1], cand.decVtx[2],
                       cand.flightDistance, cand.dcaProtonGamma,
                       cand.pProton[0], cand.pProton[1], cand.pProton[2],
                       cand.pGamma1[0], cand.pGamma1[1], cand.pGamma1[2],
                       cand.pGamma2[0], cand.pGamma2[1], cand.pGamma2[2],
                       cand.protonTpcNSigma, cand.protonTofNSigma,
                       cand.posTpcNSigmaEl, cand.negTpcNSigmaEl,
                       cand.photonMass, cand.photonAlpha, cand.photonQt, cand.photonRadius,
                       cand.photonOpeningAngle, cand.photonPointingAngle,
                       cand.rootCenter, cand.antiSigmaPointingAngle, cand.dcaToPV, cand.tiltAngle,
                       cand.protonSign,
                       cand.protonItsNCls, cand.protonTpcNCls, cand.protonDcaXY, cand.protonDcaZ,
                       cand.posItsNCls, cand.posTpcNCls, cand.negItsNCls, cand.negTpcNCls,
                       cand.photonDcaDau, cand.photonCosPAToPV, cand.photonDcaXYToPV, cand.photonDcaZToPV, cand.photonPsiPair,
                       cand.posDcaXY, cand.posDcaZ, cand.negDcaXY, cand.negDcaZ,
                       cand.posTpcNClsFindable, cand.negTpcNClsFindable, cand.posTpcChi2NCl, cand.negTpcChi2NCl,
                       cand.collisionIdCheck,
                       cand.isSignal, cand.protonIsSignal, cand.photonIsSignal,
                       cand.decVtxMC[0], cand.decVtxMC[1], cand.decVtxMC[2],
                       cand.pProtonMC[0], cand.pProtonMC[1], cand.pProtonMC[2],
                       cand.pGammaMC[0], cand.pGammaMC[1], cand.pGammaMC[2],
                       cand.pSigmaPlusMC[0], cand.pSigmaPlusMC[1], cand.pSigmaPlusMC[2],
                       cand.decayRadiusMC, cand.massMC);
    } else {
      if (fillSlimTables) {
        slimSigmaPlusCands(cand.radius,
                           cand.dcaToPV, cand.tiltAngle, cand.dcaProtonGamma,
                           cand.protonSign,
                           cand.protonDcaXY, cand.protonDcaZ,
                           cand.pProton[0], cand.pProton[1], cand.pProton[2],
                           cand.pGamma1[0], cand.pGamma1[1], cand.pGamma1[2],
                           cand.pGamma2[0], cand.pGamma2[1], cand.pGamma2[2],
                           cand.protonTpcNSigma, cand.protonTofNSigma,
                           cand.posTpcNSigmaEl, cand.negTpcNSigmaEl,
                           cand.photonMass,
                           cand.photonDcaDau, cand.photonCosPAToPV, cand.photonDcaXYToPV, cand.photonDcaZToPV, cand.photonPsiPair,
                           cand.posDcaXY, cand.posDcaZ, cand.negDcaXY, cand.negDcaZ,
                           cand.posTpcNClsFindable, cand.negTpcNClsFindable, cand.posTpcChi2NCl, cand.negTpcChi2NCl);
        return;
      }
      sigmaPlusCands(cand.decVtx[0], cand.decVtx[1], cand.decVtx[2],
                     cand.flightDistance, cand.dcaProtonGamma,
                     cand.pProton[0], cand.pProton[1], cand.pProton[2],
                     cand.pGamma1[0], cand.pGamma1[1], cand.pGamma1[2],
                     cand.pGamma2[0], cand.pGamma2[1], cand.pGamma2[2],
                     cand.protonTpcNSigma, cand.protonTofNSigma,
                     cand.posTpcNSigmaEl, cand.negTpcNSigmaEl,
                     cand.photonMass, cand.photonAlpha, cand.photonQt, cand.photonRadius,
                     cand.photonOpeningAngle, cand.photonPointingAngle,
                     cand.rootCenter, cand.antiSigmaPointingAngle, cand.dcaToPV, cand.tiltAngle,
                     cand.protonSign,
                     cand.protonItsNCls, cand.protonTpcNCls, cand.protonDcaXY, cand.protonDcaZ,
                     cand.posItsNCls, cand.posTpcNCls, cand.negItsNCls, cand.negTpcNCls,
                     cand.photonDcaDau, cand.photonCosPAToPV, cand.photonDcaXYToPV, cand.photonDcaZToPV, cand.photonPsiPair,
                     cand.posDcaXY, cand.posDcaZ, cand.negDcaXY, cand.negDcaZ,
                     cand.posTpcNClsFindable, cand.negTpcNClsFindable, cand.posTpcChi2NCl, cand.negTpcChi2NCl);
    }
  }

  // The same photon can be in candidates of several collisions (with candDeduplicatePhotons only the one with the smallest DCA to PV is written)
  template <bool IsMC>
  void deduplicateAndWriteCandidates()
  {
    std::unordered_map<uint64_t, std::vector<int>> candsByPhoton;
    for (int iCand = 0; iCand < static_cast<int>(mCandidatesOfTimeframe.size()); ++iCand) {
      candsByPhoton[mCandidatesOfTimeframe[iCand].photonId].push_back(iCand);
    }
    std::vector<bool> keep(mCandidatesOfTimeframe.size(), !candDeduplicatePhotons);
    for (const auto& [photon, cands] : candsByPhoton) {
      histos.fill(HIST("Candidate/Inclusive/hCandidatesPerPhoton"), cands.size());
      if (candDeduplicatePhotons) {
        int best = cands[0];
        for (const int& iCand : cands) {
          if (mCandidatesOfTimeframe[iCand].dcaToPV < mCandidatesOfTimeframe[best].dcaToPV) {
            best = iCand;
          }
        }
        keep[best] = true;
      }
    }

    for (int iCand = 0; iCand < static_cast<int>(mCandidatesOfTimeframe.size()); ++iCand) {
      if (!keep[iCand]) {
        continue;
      }
      const auto& cand = mCandidatesOfTimeframe[iCand];
      fillCandHist(HIST("hSelectionCounter"), cand.isSignal, 14);
      fillCandHist(HIST("h2SelectionCounterVsLegType"), cand.isSignal, 14, cand.photonLegsWithoutIts);
      if constexpr (IsMC) {
        if (cand.isSignal) {
          fillSigmaCounter(cand.matchedSigmaId, 5);
          auto genIt = mGenSigmaWritten.find(cand.matchedSigmaId);
          if (genIt != mGenSigmaWritten.end()) {
            genIt->second = true;
          }
        }
      }
      fillCandidateTables<IsMC>(cand);
    }
    mCandidatesOfTimeframe.clear();
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    auto* grpmag = ccdb->getForRun<o2::parameters::GRPMagField>(grpmagPath, mRunNumber);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    mBz = grpmag->getNominalL3Field();
    fitter.setBz(mBz);
    if (!mLut) {
      mLut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
    }
    o2::base::Propagator::Instance()->setMatLUT(mLut);
    LOG(info) << "Task initialized for run " << mRunNumber << " with magnetic field " << mBz << " kZG";
  }

  void processData(CollisionsFull const& collisions, aod::V0Datas const& v0s, TracksFull const& tracks, aod::AmbiguousTracks const& ambiguousTracks, aod::BCsWithTimestamps const& bcs)
  {
    markGoodCollisions(collisions);
    std::vector<std::vector<PhotonCand<TracksFull::iterator>>> photonsByCollision;
    if (useCustomVertexer && collisions.size() > 0) {
      auto firstBC = collisions.begin().bc_as<aod::BCsWithTimestamps>();
      initCCDB(firstBC);
      mVDriftMgr.update(firstBC.timestamp());
      photonsByCollision = findPhotonsSelfPaired<false>(collisions, tracks, ambiguousTracks, bcs);
    }

    for (const auto& collision : collisions) {
      if (std::abs(collision.posZ()) > cutZVertex || !collision.sel8()) {
        continue;
      }
      initCCDB(collision.bc_as<aod::BCsWithTimestamps>());
      histos.fill(HIST("hVertexZ"), collision.posZ());
      std::array<float, 3> pv{collision.posX(), collision.posY(), collision.posZ()};

      auto tracksThisCollision = tracks.sliceBy(tracksPerCollision, collision.globalIndex());

      std::vector<PhotonCand<TracksFull::iterator>> acceptedPhotons;
      if (useCustomVertexer) {
        acceptedPhotons = photonsByCollision[collision.globalIndex()];
      } else {
        auto v0sThisCollision = v0s.sliceBy(v0PerCollision, collision.globalIndex());
        acceptedPhotons = findPhotonsFromV0s<false>(v0sThisCollision, tracksThisCollision, pv);
      }

      std::vector<TracksFull::iterator> acceptedProtons;
      for (const auto& track : tracksThisCollision) {
        if (selectProton<false>(track)) {
          acceptedProtons.push_back(track);
        }
      }

      for (const auto& photon : acceptedPhotons) {
        for (const auto& proton : acceptedProtons) {
          buildSigmaPlusCandidate<false>(proton, photon, pv);
        }
      }
    }
    deduplicateAndWriteCandidates<false>();
  }
  PROCESS_SWITCH(Sigmaplusbuilder, processData, "Process data", true);

  void processMc(CollisionsFullMC const& collisions, aod::V0Datas const& v0s, TracksFullMC const& tracks, aod::AmbiguousTracks const& ambiguousTracks, aod::BCsWithTimestamps const& bcs, aod::McParticles const& mcParticles, aod::McCollisions const&)
  {
    // generated Sigma+ -> p pi0 in acceptance
    mGenSigmaWritten.clear();
    for (const auto& mcPart : mcParticles) {
      if (!isSigmaPlusToProtonPi0(mcPart) || std::abs(mcPart.y()) > cutRapMotherMC || mcPart.pt() < cutPtGenMC) {
        continue;
      }
      mGenSigmaWritten[mcPart.globalIndex()] = false;
      histos.fill(HIST("MC/hGenSigmaPlusPt"), mcPart.pt());
      histos.fill(HIST("MC/hSigmaPlusCounter"), 0);
    }

    constexpr int ElectronLeg = 1;
    constexpr int PositronLeg = 2;
    std::unordered_map<int, std::pair<int, int>> legsOfPhoton;
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) {
        continue;
      }
      auto mcParticle = track.template mcParticle_as<aod::McParticles>();
      if (std::abs(mcParticle.pdgCode()) == PDG_t::kProton) {
        fillSigmaCounter(findSigmaPlusMotherOfProton(mcParticle), 1);
      } else if (std::abs(mcParticle.pdgCode()) == PDG_t::kElectron) {
        int photonIndex = -1;
        std::array<float, 3> conversionVertex{};
        int sigmaId = findSigmaPlusMotherOfConversionElectron(mcParticle, photonIndex, conversionVertex);
        if (sigmaId >= 0) {
          auto& legs = legsOfPhoton[photonIndex];
          legs.first = sigmaId;
          legs.second |= (mcParticle.pdgCode() == PDG_t::kElectron ? ElectronLeg : PositronLeg);
        }
      }
    }
    for (const auto& [photonIndex, legs] : legsOfPhoton) {
      if (legs.second == (ElectronLeg | PositronLeg)) {
        fillSigmaCounter(legs.first, 3);
      }
    }

    markGoodCollisions(collisions);
    std::vector<std::vector<PhotonCand<TracksFullMC::iterator>>> photonsByCollision;
    if (useCustomVertexer && collisions.size() > 0) {
      auto firstBC = collisions.begin().bc_as<aod::BCsWithTimestamps>();
      initCCDB(firstBC);
      mVDriftMgr.update(firstBC.timestamp());
      photonsByCollision = findPhotonsSelfPaired<true>(collisions, tracks, ambiguousTracks, bcs);
    }

    for (const auto& collision : collisions) {
      if (std::abs(collision.posZ()) > cutZVertex || !collision.sel8()) {
        continue;
      }
      initCCDB(collision.bc_as<aod::BCsWithTimestamps>());
      histos.fill(HIST("hVertexZ"), collision.posZ());
      std::array<float, 3> pv{collision.posX(), collision.posY(), collision.posZ()};

      auto tracksThisCollision = tracks.sliceBy(tracksPerCollisionMC, collision.globalIndex());

      std::vector<PhotonCand<TracksFullMC::iterator>> acceptedPhotons;
      if (useCustomVertexer) {
        acceptedPhotons = photonsByCollision[collision.globalIndex()];
      } else {
        auto v0sThisCollision = v0s.sliceBy(v0PerCollision, collision.globalIndex());
        acceptedPhotons = findPhotonsFromV0s<true>(v0sThisCollision, tracksThisCollision, pv);
      }

      std::vector<TracksFullMC::iterator> acceptedProtons;
      for (const auto& track : tracksThisCollision) {
        if (selectProton<true>(track)) {
          acceptedProtons.push_back(track);
          if (track.has_mcParticle()) {
            fillSigmaCounter(findSigmaPlusMotherOfProton(track.template mcParticle_as<aod::McParticles>()), 2);
          }
        }
      }

      for (const auto& photon : acceptedPhotons) {
        for (const auto& proton : acceptedProtons) {
          buildSigmaPlusCandidate<true>(proton, photon, pv);
        }
      }
    }
    deduplicateAndWriteCandidates<true>();

    // a table row with the true values for each generated Sigma+ in acceptance without a written candidate
    for (const auto& mcPart : mcParticles) {
      auto genIt = mGenSigmaWritten.find(mcPart.globalIndex());
      if (genIt == mGenSigmaWritten.end() || genIt->second) {
        continue;
      }

      int pdgProton = mcPart.pdgCode() > 0 ? PDG_t::kProton : PDG_t::kProtonBar;
      std::array<float, 3> genDecVtx{-999.f, -999.f, -999.f};
      std::array<float, 3> genMomProton{-999.f, -999.f, -999.f};
      for (const auto& daughter : mcPart.template daughters_as<aod::McParticles>()) {
        if (daughter.pdgCode() == pdgProton) {
          genDecVtx = {daughter.vx(), daughter.vy(), daughter.vz()};
          genMomProton = {daughter.px(), daughter.py(), daughter.pz()};
          break;
        }
      }
      float genDecayRadiusMC = std::hypot(genDecVtx[0] - mcPart.vx(), genDecVtx[1] - mcPart.vy());
      float genMassMC = std::sqrt(mcPart.e() * mcPart.e() - mcPart.p() * mcPart.p());

      if (fillSlimTables) {
        slimSigmaPlusCandsMC(-999.f,
                             -999.f, -999.f, -999.f,
                             0,
                             -999.f, -999.f,
                             -999.f, -999.f, -999.f,
                             -999.f, -999.f, -999.f,
                             -999.f, -999.f, -999.f,
                             -999.f, -999.f,
                             -999.f, -999.f,
                             -999.f,
                             -999.f, -999.f, -999.f, -999.f, -999.f,
                             -999.f, -999.f, -999.f, -999.f,
                             0, 0, -999.f, -999.f,
                             false,
                             false, false, false,
                             genDecayRadiusMC, genMassMC,
                             mcPart.px(), mcPart.py(), mcPart.pz());
        continue;
      }

      sigmaPlusCandsMC(-999.f, -999.f, -999.f,
                       -999.f, -999.f,
                       -999.f, -999.f, -999.f,
                       -999.f, -999.f, -999.f,
                       -999.f, -999.f, -999.f,
                       -999.f, -999.f,
                       -999.f, -999.f,
                       -999.f, -999.f, -999.f, -999.f,
                       -999.f, -999.f,
                       -999.f, -999.f, -999.f, -999.f,
                       0,
                       0, -999, -999.f, -999.f,
                       0, -999, 0, -999,
                       -999.f, -999.f, -999.f, -999.f, -999.f,
                       -999.f, -999.f, -999.f, -999.f,
                       0, 0, -999.f, -999.f,
                       false,
                       false, false, false,
                       genDecVtx[0], genDecVtx[1], genDecVtx[2],
                       genMomProton[0], genMomProton[1], genMomProton[2],
                       -999.f, -999.f, -999.f,
                       mcPart.px(), mcPart.py(), mcPart.pz(),
                       genDecayRadiusMC, genMassMC);
    }
  }
  PROCESS_SWITCH(Sigmaplusbuilder, processMc, "Process MC", false);

  void processFindable(CollisionsFullMC const& collisions, aod::V0Datas const& v0s, TracksFullMC const& tracks, aod::McParticles const&, aod::BCsWithTimestamps const&)
  {
    mGenSigmaWritten.clear();
    constexpr int MinDauTpcCls = 90;
    for (const auto& collision : collisions) {
      if (std::abs(collision.posZ()) > cutZVertex || !collision.sel8()) {
        continue;
      }
      initCCDB(collision.bc_as<aod::BCsWithTimestamps>());
      auto tracksThisCollision = tracks.sliceBy(tracksPerCollisionMC, collision.globalIndex());
      auto v0sThisCollision = v0s.sliceBy(v0PerCollision, collision.globalIndex());
      std::array<float, 3> pv{collision.posX(), collision.posY(), collision.posZ()};
      auto acceptedPhotons = findPhotonsFromV0s<true>(v0sThisCollision, tracksThisCollision, pv);

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

          if (electronTrack.tpcNClsFound() < MinDauTpcCls || electronTrack.sign() > 0) {
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

            if (positronTrack.tpcNClsFound() < MinDauTpcCls || positronTrack.sign() < 0) {
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
            constexpr float ElectronMass = o2::constants::physics::MassElectron;
            float electronP2 = electronTrack.px() * electronTrack.px() + electronTrack.py() * electronTrack.py() + electronTrack.pz() * electronTrack.pz();
            float positronP2 = positronTrack.px() * positronTrack.px() + positronTrack.py() * positronTrack.py() + positronTrack.pz() * positronTrack.pz();
            float pairEnergy = std::sqrt(electronP2 + ElectronMass * ElectronMass) + std::sqrt(positronP2 + ElectronMass * ElectronMass);
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

            histos.fill(HIST("Findable/hPhotonSearchPresence"), 0);
            for (const auto& photon : acceptedPhotons) {
              bool sameChargeMatched = photon.posTrack.globalIndex() == positronTrack.globalIndex() && photon.negTrack.globalIndex() == electronTrack.globalIndex();
              if (sameChargeMatched) {
                histos.fill(HIST("Findable/hPhotonSearchPresence"), 1);
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
