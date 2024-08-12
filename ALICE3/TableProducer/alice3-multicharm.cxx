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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

// simple checkers
// #define biton(var, nbit) ((var) |= (static_cast<uint32_t>(1) << (nbit)))
#define bitoff(var, nbit) ((var) &= ~(static_cast<uint32_t>(1) << (nbit))) //((a) &= ~(1ULL<<(b)))
// #define bitcheck(var, nbit) ((var) & (static_cast<uint32_t>(1) << (nbit)))

using FullTracksExt = soa::Join<aod::Tracks, aod::TracksCov>;

// For MC association in pre-selection
using labeledTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
using tofTracks = soa::Join<aod::Tracks, aod::UpgradeTofs>;
using richTracks = soa::Join<aod::Tracks, aod::RICHs>;
using alice3tracks = soa::Join<aod::Tracks, aod::TracksCov, aod::Alice3DecayMaps, aod::McTrackLabels, aod::TracksDCA>;

struct alice3multicharm {
  SliceCache cache;

  // Operation and minimisation criteria
  Configurable<float> magneticField{"magneticField", 20.0f, "Magnetic field (in kilogauss)"};
  Configurable<bool> doDCAplots{"doDCAplots", true, "do daughter prong DCA plots for D mesons"};
  Configurable<bool> mcSameMotherCheck{"mcSameMotherCheck", true, "check if tracks come from the same MC mother"};
  Configurable<float> dcaXiCDaughtersSelection{"dcaXiCDaughtersSelection", 1000.0f, "DCA between XiC daughters (cm)"};
  Configurable<float> dcaXiCCDaughtersSelection{"dcaXiCCDaughtersSelection", 1000.0f, "DCA between XiCC daughters (cm)"};

  Configurable<float> piFromXiC_dcaXYconstant{"piFromXiC_dcaXYconstant", -1.0f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> piFromXiC_dcaXYpTdep{"piFromXiC_dcaXYpTdep", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> piFromXiCC_dcaXYconstant{"piFromXiCC_dcaXYconstant", -1.0f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> piFromXiCC_dcaXYpTdep{"piFromXiCC_dcaXYpTdep", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiFromXiC_dcaXYconstant{"xiFromXiC_dcaXYconstant", -1.0f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> xiFromXiC_dcaXYpTdep{"xiFromXiC_dcaXYpTdep", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};

  ConfigurableAxis axisEta{"axisEta", {8, -4.0f, +4.0f}, "#eta"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  ConfigurableAxis axisDCA{"axisDCA", {200, -100, 100}, "DCA (#mum)"};

  ConfigurableAxis axisXiMass{"axisXiMass", {200, 1.221f, 1.421f}, "Xi Inv Mass (GeV/c^{2})"};
  ConfigurableAxis axisXiCMass{"axisXiCMass", {200, 2.368f, 2.568f}, "XiC Inv Mass (GeV/c^{2})"};
  ConfigurableAxis axisXiCCMass{"axisXiCCMass", {200, 3.521f, 3.721f}, "XiCC Inv Mass (GeV/c^{2})"};

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
    ((aod::a3DecayMap::decayMap & trackSelectionPiFromXiC) == trackSelectionPiFromXiC) && aod::track::signed1Pt > 0.0f && nabs(aod::track::dcaXY) > piFromXiC_dcaXYconstant + piFromXiC_dcaXYpTdep* nabs(aod::track::signed1Pt);
  Partition<alice3tracks> tracksPiFromXiCC =
    ((aod::a3DecayMap::decayMap & trackSelectionPiFromXiCC) == trackSelectionPiFromXiCC) && aod::track::signed1Pt > 0.0f && nabs(aod::track::dcaXY) > piFromXiCC_dcaXYconstant + piFromXiCC_dcaXYpTdep* nabs(aod::track::signed1Pt);

  // Helper struct to pass candidate information
  struct {
    float dca;
    float mass;
    float pt;
    float eta;
    std::array<float, 3> xyz;
    std::array<float, 3> prong0mom;
    std::array<float, 3> prong1mom;
    std::array<float, 3> prong2mom;
    std::array<float, 21> parentTrackCovMatrix;
  } thisCandidate;

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
    t0new.getPxPyPzGlo(thisCandidate.prong0mom);
    t1new.getPxPyPzGlo(thisCandidate.prong1mom);

    // get decay vertex coordinates
    const auto& vtx = fitter.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      thisCandidate.xyz[i] = vtx[i];
    }

    // compute cov mat
    for (int ii = 0; ii < 21; ii++)
      thisCandidate.parentTrackCovMatrix[ii] = 0.0f;

    std::array<float, 21> covA = {0};
    std::array<float, 21> covB = {0};
    fitter.getTrack(0).getCovXYZPxPyPzGlo(covA);
    fitter.getTrack(1).getCovXYZPxPyPzGlo(covB);

    const int momInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
    for (int i = 0; i < 6; i++) {
      int j = momInd[i];
      thisCandidate.parentTrackCovMatrix[j] = covA[j] + covB[j];
    }

    auto covVtx = fitter.calcPCACovMatrix();
    thisCandidate.parentTrackCovMatrix[0] = covVtx(0, 0);
    thisCandidate.parentTrackCovMatrix[1] = covVtx(1, 0);
    thisCandidate.parentTrackCovMatrix[2] = covVtx(1, 1);
    thisCandidate.parentTrackCovMatrix[3] = covVtx(2, 0);
    thisCandidate.parentTrackCovMatrix[4] = covVtx(2, 1);
    thisCandidate.parentTrackCovMatrix[5] = covVtx(2, 2);

    // set relevant values
    thisCandidate.dca = TMath::Sqrt(fitter.getChi2AtPCACandidate());
    thisCandidate.mass = RecoDecay::m(array{array{thisCandidate.prong0mom[0], thisCandidate.prong0mom[1], thisCandidate.prong0mom[2]}, array{thisCandidate.prong1mom[0], thisCandidate.prong1mom[1], thisCandidate.prong1mom[2]}}, array{mass0, mass1});
    thisCandidate.pt = std::hypot(thisCandidate.prong0mom[0] + thisCandidate.prong1mom[0], thisCandidate.prong0mom[1] + thisCandidate.prong1mom[1]);
    thisCandidate.eta = RecoDecay::eta(array{thisCandidate.prong0mom[0] + thisCandidate.prong1mom[0], thisCandidate.prong0mom[1] + thisCandidate.prong1mom[1], thisCandidate.prong0mom[2] + thisCandidate.prong1mom[2]});
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
    t0.getPxPyPzGlo(thisCandidate.prong0mom);
    t1.getPxPyPzGlo(thisCandidate.prong1mom);
    t2.getPxPyPzGlo(thisCandidate.prong2mom);

    // get decay vertex coordinates
    const auto& vtx = fitter3.getPCACandidate();
    for (int i = 0; i < 3; i++) {
      thisCandidate.xyz[i] = vtx[i];
    }

    // compute cov mat
    for (int ii = 0; ii < 21; ii++)
      thisCandidate.parentTrackCovMatrix[ii] = 0.0f;

    std::array<float, 21> covA = {0};
    std::array<float, 21> covB = {0};
    std::array<float, 21> covC = {0};
    fitter3.getTrack(0).getCovXYZPxPyPzGlo(covA);
    fitter3.getTrack(1).getCovXYZPxPyPzGlo(covB);
    fitter3.getTrack(2).getCovXYZPxPyPzGlo(covC);

    const int momInd[6] = {9, 13, 14, 18, 19, 20}; // cov matrix elements for momentum component
    for (int i = 0; i < 6; i++) {
      int j = momInd[i];
      thisCandidate.parentTrackCovMatrix[j] = covA[j] + covB[j] + covC[j];
    }

    auto covVtx = fitter3.calcPCACovMatrix();
    thisCandidate.parentTrackCovMatrix[0] = covVtx(0, 0);
    thisCandidate.parentTrackCovMatrix[1] = covVtx(1, 0);
    thisCandidate.parentTrackCovMatrix[2] = covVtx(1, 1);
    thisCandidate.parentTrackCovMatrix[3] = covVtx(2, 0);
    thisCandidate.parentTrackCovMatrix[4] = covVtx(2, 1);
    thisCandidate.parentTrackCovMatrix[5] = covVtx(2, 2);

    // set relevant values
    thisCandidate.dca = TMath::Sqrt(fitter3.getChi2AtPCACandidate());
    thisCandidate.mass = RecoDecay::m(array{array{thisCandidate.prong0mom[0], thisCandidate.prong0mom[1], thisCandidate.prong0mom[2]}, array{thisCandidate.prong1mom[0], thisCandidate.prong1mom[1], thisCandidate.prong1mom[2]}, array{thisCandidate.prong2mom[0], thisCandidate.prong2mom[1], thisCandidate.prong2mom[2]}}, array{p0mass, p1mass, p2mass});
    thisCandidate.pt = std::hypot(thisCandidate.prong0mom[0] + thisCandidate.prong1mom[0] + thisCandidate.prong2mom[0], thisCandidate.prong0mom[1] + thisCandidate.prong1mom[1] + thisCandidate.prong2mom[1]);
    thisCandidate.eta = RecoDecay::eta(array{thisCandidate.prong0mom[0] + thisCandidate.prong1mom[0] + thisCandidate.prong2mom[0], thisCandidate.prong0mom[1] + thisCandidate.prong1mom[1] + thisCandidate.prong2mom[1], thisCandidate.prong0mom[2] + thisCandidate.prong1mom[2] + thisCandidate.prong2mom[2]});
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

    histos.add("h2dGenXi", "h2dGenXi", kTH2F, {axisPt, axisEta});
    histos.add("h2dGenXiC", "h2dGenXiC", kTH2F, {axisPt, axisEta});
    histos.add("h2dGenXiCC", "h2dGenXiCC", kTH2F, {axisPt, axisEta});

    histos.add("hMassXi", "hMassXi", kTH1F, {axisXiMass});
    histos.add("hMassXiC", "hMassXiC", kTH1F, {axisXiCMass});
    histos.add("hMassXiCC", "hMassXiCC", kTH1F, {axisXiCCMass});

    histos.add("hCombinationsXiC", "hCombinationsXiC", kTH1F, {axisNConsidered});
    histos.add("hCombinationsXiCC", "hCombinationsXiCC", kTH1F, {axisNConsidered});

    if (doDCAplots) {
      histos.add("h2dDCAxyVsPtXiFromXiC", "h2dDCAxyVsPtXiFromXiC", kTH2F, {axisPt, axisDCA});
      histos.add("h2dDCAxyVsPtPiFromXiC", "h2dDCAxyVsPtPiFromXiC", kTH2F, {axisPt, axisDCA});
      histos.add("h2dDCAxyVsPtPiFromXiCC", "h2dDCAxyVsPtPiFromXiCC", kTH2F, {axisPt, axisDCA});
    }
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processGenerated(aod::McParticles const&)
  {
    for (auto const& mcParticle : trueXi)
      histos.fill(HIST("h2dGenXi"), mcParticle.pt(), mcParticle.eta());
    for (auto const& mcParticle : trueXiC)
      histos.fill(HIST("h2dGenXiC"), mcParticle.pt(), mcParticle.eta());
    for (auto const& mcParticle : trueXiCC)
      histos.fill(HIST("h2dGenXiCC"), mcParticle.pt(), mcParticle.eta());
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processFindXiCC(aod::Collision const& collision, alice3tracks const&, aod::McParticles const&, aod::UpgradeCascades const& cascades)
  {
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
      for (auto const& track : tracksPiFromXiCgrouped)
        histos.fill(HIST("h2dDCAxyVsPtPiFromXiC"), track.pt(), track.dcaXY() * 1e+4);
      for (auto const& track : tracksPiFromXiCCgrouped)
        histos.fill(HIST("h2dDCAxyVsPtPiFromXiCC"), track.pt(), track.dcaXY() * 1e+4);
    }

    for (auto const& xiCand : cascades) {
      histos.fill(HIST("hMassXi"), xiCand.mXi());
      auto xi = xiCand.cascadeTrack_as<alice3tracks>(); // de-reference cascade track
      uint32_t nCombinationsC = 0;
      for (auto const& pi1c : tracksPiFromXiCgrouped) {
        if (mcSameMotherCheck && !checkSameMother(xi, pi1c))
          continue;
        for (auto const& pi2c : tracksPiFromXiCgrouped) {
          if (mcSameMotherCheck && !checkSameMother(xi, pi2c))
            continue; // keep only if same mother
          if (pi1c.globalIndex() >= pi2c.globalIndex())
            continue; // avoid same-mother, avoid double-counting

          // if I am here, it means this is a triplet to be considered for XiC vertexing.
          // will now attempt to build a three-body decay candidate with these three track rows.

          nCombinationsC++;
          if (!buildDecayCandidateThreeBody(xi, pi1c, pi2c, 1.32171, 0.139570, 0.139570))
            continue; // failed at building candidate

          const std::array<float, 3> momentumC = {
            thisCandidate.prong0mom[0] + thisCandidate.prong1mom[0] + thisCandidate.prong2mom[0],
            thisCandidate.prong0mom[1] + thisCandidate.prong1mom[1] + thisCandidate.prong2mom[1],
            thisCandidate.prong0mom[2] + thisCandidate.prong1mom[2] + thisCandidate.prong2mom[2]};

          o2::track::TrackParCov xicTrack(thisCandidate.xyz, momentumC, thisCandidate.parentTrackCovMatrix, +1);

          histos.fill(HIST("hMassXiC"), thisCandidate.mass);

          // attempt XiCC finding
          uint32_t nCombinationsCC = 0;
          for (auto const& picc : tracksPiFromXiCCgrouped) {
            // to-do: check same mother here

            nCombinationsCC++;
            o2::track::TrackParCov piccTrack = getTrackParCov(picc);
            if (!buildDecayCandidateTwoBody(xicTrack, piccTrack, 2.46793, 0.139570))
              continue; // failed at building candidate

            histos.fill(HIST("hMassXiCC"), thisCandidate.mass);
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
