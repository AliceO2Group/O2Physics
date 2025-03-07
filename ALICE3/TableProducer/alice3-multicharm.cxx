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
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//   Decay finder task for ALICE 3
//  *+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
//
//    Uses specific ALICE 3 PID and performance for studying
//    HF decays. Work in progress: use at your own risk!
//

#include <cmath>
#include <array>
#include <cstdlib>
#include <map>
#include <iterator>
#include <utility>

#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "ALICE3/DataModel/OTFTOF.h"
#include "ALICE3/DataModel/RICH.h"
#include "ALICE3/DataModel/A3DecayFinderTables.h"
#include "ALICE3/DataModel/OTFStrangeness.h"
#include "ALICE3/DataModel/OTFMulticharm.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "DetectorsVertexing/PVertexer.h"
#include "DetectorsVertexing/PVertexerHelpers.h"
#include "CommonConstants/PhysicsConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// simple checkers
// #define biton(var, nbit) ((var) |= (static_cast<uint32_t>(1) << (nbit)))
#define bitoff(var, nbit) ((var) &= ~(static_cast<uint32_t>(1) << (nbit))) //((a) &= ~(1ULL<<(b)))
#define bitcheck(var, nbit) ((var) & (static_cast<uint32_t>(1) << (nbit)))

using FullTracksExt = soa::Join<aod::Tracks, aod::TracksCov>;

// For MC association in pre-selection
using labeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
using tofTracks = soa::Join<aod::Tracks, aod::UpgradeTofs>;
using richTracks = soa::Join<aod::Tracks, aod::RICHs>;
using alice3tracks = soa::Join<aod::Tracks, aod::TracksCov, aod::Alice3DecayMaps, aod::McTrackLabels, aod::TracksDCA, aod::TracksExtraA3>;

struct alice3multicharm {
  SliceCache cache;

  Produces<aod::MCharmIndices> multiCharmIdx;
  Produces<aod::MCharmCores> multiCharmCore;

  // Operation and minimisation criteria
  Configurable<bool> fillDerivedTable{"fillDerivedTable", false, "fill MCharm[] tables (careful: memory)"};
  Configurable<float> magneticField{"magneticField", 20.0f, "Magnetic field (in kilogauss)"};
  Configurable<bool> doDCAplots{"doDCAplots", true, "do daughter prong DCA plots for D mesons"};
  Configurable<bool> mcSameMotherCheck{"mcSameMotherCheck", true, "check if tracks come from the same MC mother"};
  Configurable<float> dcaXiCDaughtersSelection{"dcaXiCDaughtersSelection", 0.002f, "DCA between XiC daughters (cm)"};
  Configurable<float> dcaXiCCDaughtersSelection{"dcaXiCCDaughtersSelection", 0.002f, "DCA between XiCC daughters (cm)"};

  Configurable<float> piFromXiC_dcaXYconstant{"piFromXiC_dcaXYconstant", 0.001f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> piFromXiC_dcaXYpTdep{"piFromXiC_dcaXYpTdep", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> piFromXiCC_dcaXYconstant{"piFromXiCC_dcaXYconstant", 0.001f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> piFromXiCC_dcaXYpTdep{"piFromXiCC_dcaXYpTdep", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiFromXiC_dcaXYconstant{"xiFromXiC_dcaXYconstant", 0.001f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiFromXiC_dcaXYpTdep{"xiFromXiC_dcaXYpTdep", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};

  Configurable<float> minPiCPt{"minPiCPt", 0.15, "Minimum pT for XiC pions"};
  Configurable<float> minPiCCPt{"minPiCCPt", 0.3, "Minimum pT for XiCC pions"};
  Configurable<float> minMultiplicity{"minMultiplicity", 0, "Minimum multiplicity"};

  Configurable<float> minXiCRadius{"minXiCRadius", 0.001, "Minimum R2D for XiC decay (cm)"};
  Configurable<float> massWindowXi{"massWindowXi", 0.015, "Mass window around Xi peak"};
  Configurable<float> massWindowXiC{"massWindowXiC", 0.015, "Mass window around XiC peak"};

  ConfigurableAxis axisEta{"axisEta", {80, -4.0f, +4.0f}, "#eta"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  ConfigurableAxis axisDCA{"axisDCA", {200, -100, 100}, "DCA (#mum)"};

  ConfigurableAxis axisXiMass{"axisXiMass", {200, 1.221f, 1.421f}, "Xi Inv Mass (GeV/c^{2})"};
  ConfigurableAxis axisXiCMass{"axisXiCMass", {200, 2.368f, 2.568f}, "XiC Inv Mass (GeV/c^{2})"};
  ConfigurableAxis axisXiCCMass{"axisXiCCMass", {200, 3.521f, 3.721f}, "XiCC Inv Mass (GeV/c^{2})"};

  ConfigurableAxis axisDCAXiCDaughters{"axisDCAXiCDaughters", {200, 0, 100}, "DCA (mum)"};
  ConfigurableAxis axisDCAXiCCDaughters{"axisDCAXiCCDaughters", {200, 0, 100}, "DCA (mum)"};

  ConfigurableAxis axisNConsidered{"axisNConsidered", {200, -0.5f, 199.5f}, "Number of considered track combinations"};

  o2::vertexing::DCAFitterN<2> fitter;
  o2::vertexing::DCAFitterN<3> fitter3;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Partition<aod::McParticles> trueXi = aod::mcparticle::pdgCode == 3312;
  Partition<aod::McParticles> trueXiC = aod::mcparticle::pdgCode == 4232;
  Partition<aod::McParticles> trueXiCC = aod::mcparticle::pdgCode == 4422;

  // filter expressions for D mesons
  static constexpr uint32_t trackSelectionPiFromXiC = 1 << kInnerTOFPion | 1 << kOuterTOFPion | 1 << kRICHPion | 1 << kTruePiFromXiC;
  static constexpr uint32_t trackSelectionPiFromXiCC = 1 << kInnerTOFPion | 1 << kOuterTOFPion | 1 << kRICHPion | 1 << kTruePiFromXiCC;

  // partitions for Xi daughters
  Partition<alice3tracks> tracksPiFromXiC =
    ((aod::a3DecayMap::decayMap & trackSelectionPiFromXiC) == trackSelectionPiFromXiC) && aod::track::signed1Pt > 0.0f && 1.0f / nabs(aod::track::signed1Pt) > minPiCPt&& nabs(aod::track::dcaXY) > piFromXiC_dcaXYconstant + piFromXiC_dcaXYpTdep* nabs(aod::track::signed1Pt);
  Partition<alice3tracks> tracksPiFromXiCC =
    ((aod::a3DecayMap::decayMap & trackSelectionPiFromXiCC) == trackSelectionPiFromXiCC) && aod::track::signed1Pt > 0.0f && 1.0f / nabs(aod::track::signed1Pt) > minPiCCPt&& nabs(aod::track::dcaXY) > piFromXiCC_dcaXYconstant + piFromXiCC_dcaXYpTdep* nabs(aod::track::signed1Pt);

  // Helper struct to pass candidate information
  struct {
    // decay properties
    float dca;
    float mass;
    float pt;
    float eta;
    std::array<float, 3> xyz;
    std::array<float, 3> prong0mom;
    std::array<float, 3> prong1mom;
    std::array<float, 3> prong2mom;
    std::array<float, 21> parentTrackCovMatrix;
  } thisXiCcandidate;

  struct {
    float dca;
    float mass;
    float pt;
    float eta;
    std::array<float, 3> xyz;
    std::array<float, 3> prong0mom;
    std::array<float, 3> prong1mom;
    std::array<float, 21> parentTrackCovMatrix;

    float etaPiCC;

    // charm daughters
    int nSiliconHitsPiCC;
    int nTPCHitsPiCC;
  } thisXiCCcandidate;

  template <typename TTrackType>
  bool buildDecayCandidateTwoBody(TTrackType const& t0, TTrackType const& t1, float mass0, float mass1)
  {
    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}
    // Move close to minima
    int nCand = 0;
    try {
      nCand = fitter.process(t0, t1);
    } catch (...) {
      return false;
    }
    if (nCand == 0) {
      return false;
    }
    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}

    o2::track::TrackParCov t0new = fitter.getTrack(0);
    o2::track::TrackParCov t1new = fitter.getTrack(1);
    t0new.getPxPyPzGlo(thisXiCCcandidate.prong0mom);
    t1new.getPxPyPzGlo(thisXiCCcandidate.prong1mom);

    // get decay vertex coordinates
    const auto& vtx = fitter.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      thisXiCCcandidate.xyz[i] = vtx[i];
    }

    // compute cov mat
    for (int ii = 0; ii < 21; ii++)
      thisXiCCcandidate.parentTrackCovMatrix[ii] = 0.0f;

    std::array<float, 21> covA = {0};
    std::array<float, 21> covB = {0};
    fitter.getTrack(0).getCovXYZPxPyPzGlo(covA);
    fitter.getTrack(1).getCovXYZPxPyPzGlo(covB);

    const int momInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
    for (int i = 0; i < 6; i++) {
      int j = momInd[i];
      thisXiCCcandidate.parentTrackCovMatrix[j] = covA[j] + covB[j];
    }

    auto covVtx = fitter.calcPCACovMatrix();
    thisXiCCcandidate.parentTrackCovMatrix[0] = covVtx(0, 0);
    thisXiCCcandidate.parentTrackCovMatrix[1] = covVtx(1, 0);
    thisXiCCcandidate.parentTrackCovMatrix[2] = covVtx(1, 1);
    thisXiCCcandidate.parentTrackCovMatrix[3] = covVtx(2, 0);
    thisXiCCcandidate.parentTrackCovMatrix[4] = covVtx(2, 1);
    thisXiCCcandidate.parentTrackCovMatrix[5] = covVtx(2, 2);

    // set relevant values
    thisXiCCcandidate.dca = TMath::Sqrt(fitter.getChi2AtPCACandidate());
    if (thisXiCCcandidate.dca > dcaXiCCDaughtersSelection) {
      return false;
    }

    thisXiCCcandidate.mass = RecoDecay::m(array{array{thisXiCCcandidate.prong0mom[0], thisXiCCcandidate.prong0mom[1], thisXiCCcandidate.prong0mom[2]}, array{thisXiCCcandidate.prong1mom[0], thisXiCCcandidate.prong1mom[1], thisXiCCcandidate.prong1mom[2]}}, array{mass0, mass1});
    thisXiCCcandidate.pt = std::hypot(thisXiCCcandidate.prong0mom[0] + thisXiCCcandidate.prong1mom[0], thisXiCCcandidate.prong0mom[1] + thisXiCCcandidate.prong1mom[1]);
    thisXiCCcandidate.eta = RecoDecay::eta(array{thisXiCCcandidate.prong0mom[0] + thisXiCCcandidate.prong1mom[0], thisXiCCcandidate.prong0mom[1] + thisXiCCcandidate.prong1mom[1], thisXiCCcandidate.prong0mom[2] + thisXiCCcandidate.prong1mom[2]});
    return true;
  }

  template <typename TTrackType1, typename TTrackType2, typename TTrackType3>
  bool buildDecayCandidateThreeBody(TTrackType1 const& prong0, TTrackType2 const& prong1, TTrackType3 const& prong2, float p0mass, float p1mass, float p2mass)
  {
    o2::track::TrackParCov t0 = getTrackParCov(prong0);
    o2::track::TrackParCov t1 = getTrackParCov(prong1);
    o2::track::TrackParCov t2 = getTrackParCov(prong2);

    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}
    // Move close to minima
    int nCand = 0;
    try {
      nCand = fitter3.process(t0, t1, t2);
    } catch (...) {
      return false;
    }
    if (nCand == 0) {
      return false;
    }
    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}

    t0 = fitter3.getTrack(0);
    t1 = fitter3.getTrack(1);
    t2 = fitter3.getTrack(2);
    t0.getPxPyPzGlo(thisXiCcandidate.prong0mom);
    t1.getPxPyPzGlo(thisXiCcandidate.prong1mom);
    t2.getPxPyPzGlo(thisXiCcandidate.prong2mom);

    // get decay vertex coordinates
    const auto& vtx = fitter3.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      thisXiCcandidate.xyz[i] = vtx[i];
    }

    // compute cov mat
    for (int ii = 0; ii < 21; ii++)
      thisXiCcandidate.parentTrackCovMatrix[ii] = 0.0f;

    std::array<float, 21> covA = {0};
    std::array<float, 21> covB = {0};
    std::array<float, 21> covC = {0};
    fitter3.getTrack(0).getCovXYZPxPyPzGlo(covA);
    fitter3.getTrack(1).getCovXYZPxPyPzGlo(covB);
    fitter3.getTrack(2).getCovXYZPxPyPzGlo(covC);

    const int momInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
    for (int i = 0; i < 6; i++) {
      int j = momInd[i];
      thisXiCcandidate.parentTrackCovMatrix[j] = covA[j] + covB[j] + covC[j];
    }

    auto covVtx = fitter3.calcPCACovMatrix();
    thisXiCcandidate.parentTrackCovMatrix[0] = covVtx(0, 0);
    thisXiCcandidate.parentTrackCovMatrix[1] = covVtx(1, 0);
    thisXiCcandidate.parentTrackCovMatrix[2] = covVtx(1, 1);
    thisXiCcandidate.parentTrackCovMatrix[3] = covVtx(2, 0);
    thisXiCcandidate.parentTrackCovMatrix[4] = covVtx(2, 1);
    thisXiCcandidate.parentTrackCovMatrix[5] = covVtx(2, 2);

    // set relevant values
    thisXiCcandidate.dca = TMath::Sqrt(fitter3.getChi2AtPCACandidate());
    if (thisXiCcandidate.dca > dcaXiCDaughtersSelection) {
      return false;
    }
    thisXiCcandidate.mass = RecoDecay::m(array{array{thisXiCcandidate.prong0mom[0], thisXiCcandidate.prong0mom[1], thisXiCcandidate.prong0mom[2]}, array{thisXiCcandidate.prong1mom[0], thisXiCcandidate.prong1mom[1], thisXiCcandidate.prong1mom[2]}, array{thisXiCcandidate.prong2mom[0], thisXiCcandidate.prong2mom[1], thisXiCcandidate.prong2mom[2]}}, array{p0mass, p1mass, p2mass});
    thisXiCcandidate.pt = std::hypot(thisXiCcandidate.prong0mom[0] + thisXiCcandidate.prong1mom[0] + thisXiCcandidate.prong2mom[0], thisXiCcandidate.prong0mom[1] + thisXiCcandidate.prong1mom[1] + thisXiCcandidate.prong2mom[1]);
    thisXiCcandidate.eta = RecoDecay::eta(array{thisXiCcandidate.prong0mom[0] + thisXiCcandidate.prong1mom[0] + thisXiCcandidate.prong2mom[0], thisXiCcandidate.prong0mom[1] + thisXiCcandidate.prong1mom[1] + thisXiCcandidate.prong2mom[1], thisXiCcandidate.prong0mom[2] + thisXiCcandidate.prong1mom[2] + thisXiCcandidate.prong2mom[2]});
    return true;
  }

  /// function to check if tracks have the same mother in MC
  template <typename TTrackType1, typename TTrackType2>
  bool checkSameMother(TTrackType1 const& track1, TTrackType2 const& track2)
  {
    bool returnValue = false;
    // Association check
    // There might be smarter ways of doing this in the future
    if (track1.has_mcParticle() && track2.has_mcParticle()) {
      auto mcParticle1 = track1.template mcParticle_as<aod::McParticles>();
      auto mcParticle2 = track2.template mcParticle_as<aod::McParticles>();
      if (mcParticle1.has_mothers() && mcParticle2.has_mothers()) {
        for (auto& mcParticleMother1 : mcParticle1.template mothers_as<aod::McParticles>()) {
          for (auto& mcParticleMother2 : mcParticle2.template mothers_as<aod::McParticles>()) {
            if (mcParticleMother1.globalIndex() == mcParticleMother2.globalIndex()) {
              returnValue = true;
            }
          }
        }
      }
    } // end association check
    return returnValue;
  }

  // Association check for the XiCC pion
  template <typename TTrackType1, typename TTrackType2>
  bool checkSameMotherExtra(TTrackType1 const& track1, TTrackType2 const& track2)
  {
    bool returnValue = false;
    // This might perhaps be a bit excessive
    // Could be joined with `checkSameMother` but leaving as is for now
    if (track1.has_mcParticle() && track2.has_mcParticle()) {
      auto mcParticle1 = track1.template mcParticle_as<aod::McParticles>();
      auto mcParticle2 = track2.template mcParticle_as<aod::McParticles>();
      if (mcParticle1.has_mothers() && mcParticle2.has_mothers()) {
        for (auto& mcParticleMother1 : mcParticle1.template mothers_as<aod::McParticles>()) {
          if (mcParticleMother1.has_mothers()) {
            for (auto& mcParticleGrandMother1 : mcParticleMother1.template mothers_as<aod::McParticles>()) {
              for (auto& mcParticleMother2 : mcParticle2.template mothers_as<aod::McParticles>()) {
                if (mcParticleGrandMother1.globalIndex() == mcParticleMother2.globalIndex()) {
                  returnValue = true;
                }
              }
            }
          }
        }
      }
    } // end association check
    return returnValue;
  }

  void init(InitContext&)
  {
    // initialize O2 2-prong fitter (only once)
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);
    fitter.setWeightedFinalPCA(false);
    fitter.setBz(magneticField);
    fitter.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrNONE);

    fitter3.setPropagateToPCA(true);
    fitter3.setMaxR(200.);
    fitter3.setMinParamChange(1e-3);
    fitter3.setMinRelChi2Change(0.9);
    fitter3.setMaxDZIni(1e9);
    fitter3.setMaxChi2(1e9);
    fitter3.setUseAbsDCA(true);
    fitter3.setWeightedFinalPCA(false);
    fitter3.setBz(magneticField);
    fitter3.setMatCorrType(o2::base::Propagator::MatCorrType::USEMatCorrNONE);

    // This histogram bookkeeps the attempts at DCA minimization and their eventual
    // failure rates.
    // --- 0: attempt XiC, 1: success XiC
    // --- 2: attempt XiCC, 3: success XiCC
    histos.add("hCharmBuilding", "hCharmBuilding", kTH1D, {{10, -0.5, 9.5f}});

    histos.add("h2dGenXi", "h2dGenXi", kTH2D, {axisPt, axisEta});
    histos.add("h2dGenXiC", "h2dGenXiC", kTH2D, {axisPt, axisEta});
    histos.add("h2dGenXiCC", "h2dGenXiCC", kTH2D, {axisPt, axisEta});

    histos.add("hMassXi", "hMassXi", kTH1D, {axisXiMass});
    histos.add("hMassXiC", "hMassXiC", kTH1D, {axisXiCMass});
    histos.add("hMassXiCC", "hMassXiCC", kTH1D, {axisXiCCMass});

    histos.add("hEtaXiCC", "hEtaXiCC", kTH1D, {axisEta});
    histos.add("hPtXiCC", "hPtXiCC", kTH1D, {axisPt});
    histos.add("hMcPtXiCC", "hMcPtXiCC", kTH1D, {axisPt});
    histos.add("h3dMassXiCC", "h3dMassXiCC", kTH3D, {axisPt, axisEta, axisXiCCMass});

    histos.add("hDCAXiCDaughters", "hDCAXiCDaughters", kTH1D, {axisDCAXiCDaughters});
    histos.add("hDCAXiCCDaughters", "hDCAXiCCDaughters", kTH1D, {axisDCAXiCCDaughters});

    // These histograms bookkeep the exact number of combinations attempted
    // CombinationsXiC: triplets Xi-pi-pi considered per Xi
    // CombinationsXiCC: doublets XiC-pi considered per XiC
    histos.add("hCombinationsXiC", "hCombinationsXiC", kTH1D, {axisNConsidered});
    histos.add("hCombinationsXiCC", "hCombinationsXiCC", kTH1D, {axisNConsidered});
    histos.add("hNCollisions", "hNCollisions", kTH1D, {{2, 0.5, 2.5}});
    if (doDCAplots) {
      histos.add("h2dDCAxyVsPtXiFromXiC", "h2dDCAxyVsPtXiFromXiC", kTH2D, {axisPt, axisDCA});
      histos.add("h2dDCAxyVsPtPiFromXiC", "h2dDCAxyVsPtPiFromXiC", kTH2D, {axisPt, axisDCA});
      histos.add("h2dDCAxyVsPtPiFromXiCC", "h2dDCAxyVsPtPiFromXiCC", kTH2D, {axisPt, axisDCA});
    }
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processGenerated(aod::McParticles const&)
  {
    for (auto const& mcParticle : trueXi)
      histos.fill(HIST("h2dGenXi"), mcParticle.pt(), mcParticle.eta());
    for (auto const& mcParticle : trueXiC)
      histos.fill(HIST("h2dGenXiC"), mcParticle.pt(), mcParticle.eta());
    for (auto const& mcParticle : trueXiCC) {
      histos.fill(HIST("h2dGenXiCC"), mcParticle.pt(), mcParticle.eta());
      histos.fill(HIST("hMcPtXiCC"), mcParticle.pt());
    }
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processFindXiCC(aod::Collision const& collision, alice3tracks const& tracks, aod::McParticles const&, aod::UpgradeCascades const& cascades)
  {
    histos.fill(HIST("hNCollisions"), 1);
    if (tracks.size() < minMultiplicity)
      return;

    histos.fill(HIST("hNCollisions"), 2);

    // group with this collision
    // n.b. cascades do not need to be grouped, being used directly in iterator-grouping
    auto tracksPiFromXiCgrouped = tracksPiFromXiC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracksPiFromXiCCgrouped = tracksPiFromXiCC->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (doDCAplots) {
      for (auto const& cascade : cascades) {
        if (cascade.has_cascadeTrack()) {
          auto track = cascade.cascadeTrack_as<alice3tracks>(); // de-reference cascade track
          histos.fill(HIST("h2dDCAxyVsPtXiFromXiC"), track.pt(), track.dcaXY() * 1e+4);
        } else {
          LOGF(info, "Damn, something is wrong");
        }
      }
      for (auto const& track : tracks) {
        if (bitcheck(track.decayMap(), kTruePiFromXiC))
          histos.fill(HIST("h2dDCAxyVsPtPiFromXiC"), track.pt(), track.dcaXY() * 1e+4);
        if (bitcheck(track.decayMap(), kTruePiFromXiCC))
          histos.fill(HIST("h2dDCAxyVsPtPiFromXiCC"), track.pt(), track.dcaXY() * 1e+4);
      }
    }

    for (auto const& xiCand : cascades) {
      histos.fill(HIST("hMassXi"), xiCand.mXi());

      if (std::abs(xiCand.mXi() - o2::constants::physics::MassXiMinus) > massWindowXi)
        continue; // out of mass region

      uint32_t nCombinationsC = 0;
      auto xi = xiCand.cascadeTrack_as<alice3tracks>();    // de-reference cascade track
      auto piFromXi = xiCand.bachTrack_as<alice3tracks>(); // de-reference bach track
      auto piFromLa = xiCand.negTrack_as<alice3tracks>();  // de-reference neg track
      auto prFromLa = xiCand.posTrack_as<alice3tracks>();  // de-reference pos track

      if (!bitcheck(xi.decayMap(), kTrueXiFromXiC))
        continue;

      for (auto const& pi1c : tracksPiFromXiCgrouped) {
        if (mcSameMotherCheck && !checkSameMother(xi, pi1c))
          continue;
        if (xiCand.posTrackId() == pi1c.globalIndex() || xiCand.negTrackId() == pi1c.globalIndex() || xiCand.bachTrackId() == pi1c.globalIndex())
          continue; // avoid using any track that was already used
        if (pi1c.pt() < minPiCPt)
          continue;

        // second pion from XiC decay for starts here
        for (auto const& pi2c : tracksPiFromXiCgrouped) {

          if (mcSameMotherCheck && !checkSameMother(xi, pi2c))
            continue; // keep only if same mother
          if (pi1c.globalIndex() >= pi2c.globalIndex())
            continue; // avoid same-mother, avoid double-counting
          if (xiCand.posTrackId() == pi2c.globalIndex() || xiCand.negTrackId() == pi2c.globalIndex() || xiCand.bachTrackId() == pi2c.globalIndex())
            continue; // avoid using any track that was already used
          if (pi2c.pt() < minPiCPt)
            continue;

          // if I am here, it means this is a triplet to be considered for XiC vertexing.
          // will now attempt to build a three-body decay candidate with these three track rows.

          nCombinationsC++;
          histos.fill(HIST("hCharmBuilding"), 0.0f);
          if (!buildDecayCandidateThreeBody(xi, pi1c, pi2c, o2::constants::physics::MassXiMinus, o2::constants::physics::MassPionCharged, o2::constants::physics::MassPionCharged))
            continue; // failed at building candidate

          if (std::abs(thisXiCcandidate.mass - o2::constants::physics::MassXiCPlus) > massWindowXiC)
            continue; // out of mass region
          histos.fill(HIST("hCharmBuilding"), 1.0f);

          const std::array<float, 3> momentumC = {
            thisXiCcandidate.prong0mom[0] + thisXiCcandidate.prong1mom[0] + thisXiCcandidate.prong2mom[0],
            thisXiCcandidate.prong0mom[1] + thisXiCcandidate.prong1mom[1] + thisXiCcandidate.prong2mom[1],
            thisXiCcandidate.prong0mom[2] + thisXiCcandidate.prong1mom[2] + thisXiCcandidate.prong2mom[2]};

          o2::track::TrackParCov xicTrack(thisXiCcandidate.xyz, momentumC, thisXiCcandidate.parentTrackCovMatrix, +1);

          if (std::hypot(thisXiCcandidate.xyz[0], thisXiCcandidate.xyz[1]) < minXiCRadius)
            continue; // do not take if radius too small, likely a primary combination

          o2::dataformats::DCA dcaInfo;
          float xicdcaXY = 1e+10;
          o2::track::TrackParCov xicTrackCopy(xicTrack); // paranoia

          o2::vertexing::PVertex primaryVertex;
          primaryVertex.setXYZ(collision.posX(), collision.posY(), collision.posZ());

          if (xicTrackCopy.propagateToDCA(primaryVertex, magneticField, &dcaInfo)) {
            xicdcaXY = dcaInfo.getY();
          }

          histos.fill(HIST("hMassXiC"), thisXiCcandidate.mass);
          histos.fill(HIST("hDCAXiCDaughters"), thisXiCcandidate.dca * 1e+4);

          // attempt XiCC finding
          uint32_t nCombinationsCC = 0;
          for (auto const& picc : tracksPiFromXiCCgrouped) {

            if (xiCand.posTrackId() == picc.globalIndex() || xiCand.negTrackId() == picc.globalIndex() || xiCand.bachTrackId() == picc.globalIndex())
              continue; // avoid using any track that was already used
            if (picc.pt() < minPiCCPt)
              continue;
            if (mcSameMotherCheck && !checkSameMotherExtra(xi, picc))
              continue;
            o2::track::TrackParCov piccTrack = getTrackParCov(picc);
            nCombinationsCC++;
            histos.fill(HIST("hCharmBuilding"), 2.0f);
            if (!buildDecayCandidateTwoBody(xicTrack, piccTrack, o2::constants::physics::MassXiCPlus, o2::constants::physics::MassPionCharged))
              continue; // failed at building candidate

            histos.fill(HIST("hCharmBuilding"), 3.0f);
            histos.fill(HIST("hMassXiCC"), thisXiCCcandidate.mass);
            histos.fill(HIST("hPtXiCC"), thisXiCCcandidate.pt);
            histos.fill(HIST("hEtaXiCC"), thisXiCCcandidate.eta);
            histos.fill(HIST("h3dMassXiCC"), thisXiCCcandidate.pt, thisXiCCcandidate.eta, thisXiCCcandidate.mass);
            histos.fill(HIST("hDCAXiCCDaughters"), thisXiCCcandidate.dca * 1e+4);

            const std::array<float, 3> momentumCC = {
              thisXiCCcandidate.prong0mom[0] + thisXiCCcandidate.prong1mom[0],
              thisXiCCcandidate.prong0mom[1] + thisXiCCcandidate.prong1mom[1],
              thisXiCCcandidate.prong0mom[2] + thisXiCCcandidate.prong1mom[2]};

            o2::track::TrackParCov xiccTrack(thisXiCCcandidate.xyz, momentumCC, thisXiCCcandidate.parentTrackCovMatrix, +2);

            float xiccdcaXY = 1e+10;
            if (xiccTrack.propagateToDCA(primaryVertex, magneticField, &dcaInfo)) {
              xiccdcaXY = dcaInfo.getY();
            }

            // produce multi-charm table for posterior analysis
            if (fillDerivedTable) {
              multiCharmCore(
                thisXiCcandidate.dca, thisXiCCcandidate.dca,
                thisXiCcandidate.mass, thisXiCCcandidate.mass,
                thisXiCCcandidate.pt, thisXiCCcandidate.eta,
                xi.nSiliconHits(), piFromXi.nSiliconHits(),
                piFromLa.nSiliconHits(), prFromLa.nSiliconHits(),
                pi1c.nSiliconHits(), pi2c.nSiliconHits(), picc.nSiliconHits(),
                piFromXi.nTPCHits(), piFromLa.nTPCHits(), prFromLa.nTPCHits(),
                pi1c.nTPCHits(), pi2c.nTPCHits(), picc.nTPCHits(),
                xi.dcaXY(), xicdcaXY, xiccdcaXY,
                piFromXi.dcaXY(), piFromLa.dcaXY(), prFromLa.dcaXY(),
                pi1c.dcaXY(), pi2c.dcaXY(), picc.dcaXY());
            }
          }
          histos.fill(HIST("hCombinationsXiCC"), nCombinationsCC);
        }
      }
      histos.fill(HIST("hCombinationsXiC"), nCombinationsC);
    }
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  //*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*
  PROCESS_SWITCH(alice3multicharm, processGenerated, "fill MC-only histograms", true);
  PROCESS_SWITCH(alice3multicharm, processFindXiCC, "find XiCC baryons", true);
  //*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<alice3multicharm>(cfgc)};
}
