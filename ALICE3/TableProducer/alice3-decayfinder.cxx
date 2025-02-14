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

struct alice3decayFinder {
  SliceCache cache;

  // Operation and minimisation criteria
  Configurable<float> magneticField{"magneticField", 20.0f, "Magnetic field (in kilogauss)"};
  Configurable<bool> doDCAplotsD{"doDCAplotsD", true, "do daughter prong DCA plots for D mesons"};
  Configurable<bool> doDCAplotsLc{"doDCAplotsLc", true, "do daughter prong DCA plots for Lc baryons"};
  Configurable<bool> mcSameMotherCheck{"mcSameMotherCheck", true, "check if tracks come from the same MC mother"};
  Configurable<float> dcaDaughtersSelection{"dcaDaughtersSelection", 1000.0f, "DCA between daughters (cm)"};

  Configurable<float> piFromD_dcaXYconstant{"piFromD_dcaXYconstant", -1.0f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> piFromD_dcaXYpTdep{"piFromD_dcaXYpTdep", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> kaFromD_dcaXYconstant{"kaFromD_dcaXYconstant", -1.0f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> kaFromD_dcaXYpTdep{"kaFromD_dcaXYpTdep", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};

  Configurable<float> DCosPA{"DCosPA", 0.99, " Cos of pointing angle: pt < 3 GeV"};
  Configurable<float> DCosPAHighPt{"DCosPAHighPt", 0.995, " Cos of pointing angle: 3 GeV < pt"};
  Configurable<float> DCosPAxy{"DCosPAxy", 0.99, " Cos of pointing angle xy: pt < 3 GeV"};
  Configurable<float> DCosPAxyHighPt{"DCosPAxyHighPt", 0.995, " Cos of pointing angle xy: 3 GeV < pt"};
  Configurable<float> DCosThetaStarLowPt{"DCosThetaStarLowPt", 0.8, "Cos theta; pt < 9"};
  Configurable<float> DCosThetaStarHighPt{"DCosThetaStarHighPt", 0.9, "Cos theta; 9 < pt < 16"};
  Configurable<float> DCosThetaStarVHighPt{"DCosThetaStarVHighPt", 1.0, "Cos theta; 16 < pt"};
  Configurable<float> DDauDecayLength{"DDauDecayLength", 3, "|Normalized dau decay length| > [0]"};

  Configurable<float> piFromLc_dcaXYconstant{"piFromLc_dcaXYconstant", -1.0f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> piFromLc_dcaXYpTdep{"piFromLc_dcaXYpTdep", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> kaFromLc_dcaXYconstant{"kaFromLc_dcaXYconstant", -1.0f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> kaFromLc_dcaXYpTdep{"kaFromLc_dcaXYpTdep", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> prFromLc_dcaXYconstant{"prFromLc_dcaXYconstant", -1.0f, "[0] in |DCAxy| > [0]+[1]/pT"};
  Configurable<float> prFromLc_dcaXYpTdep{"prFromLc_dcaXYpTdep", 0.0, "[1] in |DCAxy| > [0]+[1]/pT"};

  ConfigurableAxis axisEta{"axisEta", {8, -4.0f, +4.0f}, "#eta"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  ConfigurableAxis axisDCA{"axisDCA", {200, -100, 100}, "DCA (#mum)"};
  ConfigurableAxis axisDCADaughters{"axisDCADaughters", {200, 0, 100}, "DCA (#mum)"};
  ConfigurableAxis axisDMass{"axisDMass", {200, 1.765f, 1.965f}, "D Inv Mass (GeV/c^{2})"};
  ConfigurableAxis axisLcMass{"axisLcMass", {200, 2.186f, 2.386f}, "#Lambda_{c} Inv Mass (GeV/c^{2})"};

  o2::vertexing::DCAFitterN<2> fitter;
  o2::vertexing::DCAFitterN<3> fitter3;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Partition<aod::McParticles> trueD = aod::mcparticle::pdgCode == 421;
  Partition<aod::McParticles> trueDbar = aod::mcparticle::pdgCode == -421;
  Partition<aod::McParticles> trueLc = aod::mcparticle::pdgCode == 4122;
  Partition<aod::McParticles> trueLcbar = aod::mcparticle::pdgCode == -4122;

  // filter expressions for D mesons
  static constexpr uint32_t trackSelectionPiPlusFromD = 1 << kInnerTOFPion | 1 << kOuterTOFPion | 1 << kRICHPion | 1 << kTruePiPlusFromD;
  static constexpr uint32_t trackSelectionPiMinusFromD = 1 << kInnerTOFPion | 1 << kOuterTOFPion | 1 << kRICHPion | 1 << kTruePiMinusFromD;
  static constexpr uint32_t trackSelectionKaPlusFromD = 1 << kInnerTOFKaon | 1 << kOuterTOFKaon | 1 << kRICHKaon | 1 << kTrueKaPlusFromD;
  static constexpr uint32_t trackSelectionKaMinusFromD = 1 << kInnerTOFKaon | 1 << kOuterTOFKaon | 1 << kRICHKaon | 1 << kTrueKaMinusFromD;

  // filter expressions for Lambdac baryons
  static constexpr uint32_t trackSelectionPiPlusFromLc = 1 << kInnerTOFPion | 1 << kOuterTOFPion | 1 << kRICHPion | 1 << kTruePiPlusFromLc;
  static constexpr uint32_t trackSelectionKaPlusFromLc = 1 << kInnerTOFKaon | 1 << kOuterTOFKaon | 1 << kRICHKaon | 1 << kTrueKaPlusFromLc;
  static constexpr uint32_t trackSelectionPrPlusFromLc = 1 << kInnerTOFProton | 1 << kOuterTOFProton | 1 << kRICHProton | 1 << kTruePrPlusFromLc;
  static constexpr uint32_t trackSelectionPiMinusFromLc = 1 << kInnerTOFPion | 1 << kOuterTOFPion | 1 << kRICHPion | 1 << kTruePiMinusFromLc;
  static constexpr uint32_t trackSelectionKaMinusFromLc = 1 << kInnerTOFKaon | 1 << kOuterTOFKaon | 1 << kRICHKaon | 1 << kTrueKaMinusFromLc;
  static constexpr uint32_t trackSelectionPrMinusFromLc = 1 << kInnerTOFProton | 1 << kOuterTOFProton | 1 << kRICHProton | 1 << kTruePrMinusFromLc;

  // partitions for D mesons
  Partition<alice3tracks> tracksPiPlusFromD =
    ((aod::a3DecayMap::decayMap & trackSelectionPiPlusFromD) == trackSelectionPiPlusFromD) &&
    aod::track::signed1Pt > 0.0f &&
    nabs(aod::track::dcaXY) > piFromD_dcaXYconstant + piFromD_dcaXYpTdep* nabs(aod::track::signed1Pt);
  Partition<alice3tracks> tracksPiMinusFromD =
    ((aod::a3DecayMap::decayMap & trackSelectionPiMinusFromD) == trackSelectionPiMinusFromD) && aod::track::signed1Pt < 0.0f && nabs(aod::track::dcaXY) > piFromD_dcaXYconstant + piFromD_dcaXYpTdep* nabs(aod::track::signed1Pt);
  Partition<alice3tracks> tracksKaPlusFromD =
    ((aod::a3DecayMap::decayMap & trackSelectionKaPlusFromD) == trackSelectionKaPlusFromD) && aod::track::signed1Pt > 0.0f && nabs(aod::track::dcaXY) > kaFromD_dcaXYconstant + kaFromD_dcaXYpTdep* nabs(aod::track::signed1Pt);
  Partition<alice3tracks> tracksKaMinusFromD =
    ((aod::a3DecayMap::decayMap & trackSelectionKaMinusFromD) == trackSelectionKaMinusFromD) && aod::track::signed1Pt < 0.0f && nabs(aod::track::dcaXY) > kaFromD_dcaXYconstant + kaFromD_dcaXYpTdep* nabs(aod::track::signed1Pt);

  // partitions for Lc baryons
  Partition<alice3tracks> tracksPiPlusFromLc =
    ((aod::a3DecayMap::decayMap & trackSelectionPiPlusFromLc) == trackSelectionPiPlusFromLc) && aod::track::signed1Pt > 0.0f && nabs(aod::track::dcaXY) > piFromLc_dcaXYconstant + piFromLc_dcaXYpTdep* nabs(aod::track::signed1Pt);
  Partition<alice3tracks> tracksKaPlusFromLc =
    ((aod::a3DecayMap::decayMap & trackSelectionKaPlusFromLc) == trackSelectionKaPlusFromLc) && aod::track::signed1Pt > 0.0f && nabs(aod::track::dcaXY) > kaFromLc_dcaXYconstant + kaFromLc_dcaXYpTdep* nabs(aod::track::signed1Pt);
  Partition<alice3tracks> tracksPrPlusFromLc =
    ((aod::a3DecayMap::decayMap & trackSelectionPrPlusFromLc) == trackSelectionPrPlusFromLc) && aod::track::signed1Pt > 0.0f && nabs(aod::track::dcaXY) > prFromLc_dcaXYconstant + prFromLc_dcaXYpTdep* nabs(aod::track::signed1Pt);
  // partitions for Lc baryons
  Partition<alice3tracks> tracksPiMinusFromLc =
    ((aod::a3DecayMap::decayMap & trackSelectionPiMinusFromLc) == trackSelectionPiMinusFromLc) && aod::track::signed1Pt < 0.0f && nabs(aod::track::dcaXY) > piFromLc_dcaXYconstant + piFromLc_dcaXYpTdep* nabs(aod::track::signed1Pt);
  Partition<alice3tracks> tracksKaMinusFromLc =
    ((aod::a3DecayMap::decayMap & trackSelectionKaMinusFromLc) == trackSelectionKaMinusFromLc) && aod::track::signed1Pt < 0.0f && nabs(aod::track::dcaXY) > kaFromLc_dcaXYconstant + kaFromLc_dcaXYpTdep* nabs(aod::track::signed1Pt);
  Partition<alice3tracks> tracksPrMinusFromLc =
    ((aod::a3DecayMap::decayMap & trackSelectionPrMinusFromLc) == trackSelectionPrMinusFromLc) && aod::track::signed1Pt < 0.0f && nabs(aod::track::dcaXY) > prFromLc_dcaXYconstant + prFromLc_dcaXYpTdep* nabs(aod::track::signed1Pt);

  // Helper struct to pass candidate information
  struct {
    float dcaDau;
    float mass;
    std::array<float, 3> posSV;
    std::array<float, 3> P;
    float pt;
    float eta;
    float cosPA;
    float cosPAxy;
    float cosThetaStar;
    float normalizedDecayLength;
  } dmeson;

  struct {
    float dcaDau;
    float mass;
    float pt;
    float eta;
  } lcbaryon;

  template <typename TTrackType>
  bool buildDecayCandidateTwoBody(TTrackType const& posTrackRow, TTrackType const& negTrackRow, float posMass, float negMass)
  {
    o2::track::TrackParCov posTrack = getTrackParCov(posTrackRow);
    o2::track::TrackParCov negTrack = getTrackParCov(negTrackRow);

    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}
    // Move close to minima
    int nCand = 0;
    try {
      nCand = fitter.process(posTrack, negTrack);
    } catch (...) {
      return false;
    }
    if (nCand == 0) {
      return false;
    }
    //}-{}-{}-{}-{}-{}-{}-{}-{}-{}

    posTrack = fitter.getTrack(0);
    negTrack = fitter.getTrack(1);
    std::array<float, 3> posP;
    std::array<float, 3> negP;
    posTrack.getPxPyPzGlo(posP);
    negTrack.getPxPyPzGlo(negP);
    dmeson.dcaDau = TMath::Sqrt(fitter.getChi2AtPCACandidate());

    // return mass
    dmeson.mass = RecoDecay::m(array{array{posP[0], posP[1], posP[2]}, array{negP[0], negP[1], negP[2]}}, array{posMass, negMass});
    dmeson.pt = std::hypot(posP[0] + negP[0], posP[1] + negP[1]);
    dmeson.eta = RecoDecay::eta(array{posP[0] + negP[0], posP[1] + negP[1], posP[2] + negP[2]});
    const auto posSV = fitter.getPCACandidate();
    dmeson.posSV[0] = posSV[0];
    dmeson.posSV[1] = posSV[1];
    dmeson.posSV[2] = posSV[2];
    o2::track::TrackParCov parentTrack = fitter.createParentTrackParCov();
    parentTrack.getPxPyPzGlo(dmeson.P);
    dmeson.cosThetaStar = RecoDecay::cosThetaStar(std::array{std::array{posP[0], posP[1], posP[2]}, std::array{negP[0], negP[1], negP[2]}}, std::array{posMass, negMass}, dmeson.mass, 0);
    return true;
  }

  template <typename TTrackType>
  bool buildDecayCandidateThreeBody(TTrackType const& prong0, TTrackType const& prong1, TTrackType const& prong2, float p0mass, float p1mass, float p2mass)
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
    std::array<float, 3> P0;
    std::array<float, 3> P1;
    std::array<float, 3> P2;
    t0.getPxPyPzGlo(P0);
    t1.getPxPyPzGlo(P1);
    t2.getPxPyPzGlo(P2);

    lcbaryon.dcaDau = TMath::Sqrt(fitter3.getChi2AtPCACandidate());
    if (lcbaryon.dcaDau > dcaDaughtersSelection)
      return false;

    // return mass
    lcbaryon.mass = RecoDecay::m(array{array{P0[0], P0[1], P0[2]}, array{P1[0], P1[1], P1[2]}, array{P2[0], P2[1], P2[2]}}, array{p0mass, p1mass, p2mass});
    lcbaryon.pt = std::hypot(P0[0] + P1[0] + P2[0], P0[1] + P1[1] + P2[1]);
    lcbaryon.eta = RecoDecay::eta(array{P0[0] + P1[0] + P2[0], P0[1] + P1[1] + P2[1], P0[2] + P1[2] + P2[2]});
    return true;
  }

  /// function to check if tracks have the same mother in MC
  template <typename TTrackType>
  bool checkSameMother(TTrackType const& track1, TTrackType const& track2)
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

    if (doprocessFindDmesons) {
      histos.add("h2dGenD", "h2dGenD", kTH2F, {axisPt, axisEta});
      histos.add("h2dGenDbar", "h2dGenDbar", kTH2F, {axisPt, axisEta});
      histos.add("h3dRecD", "h2dRecD", kTH3F, {axisPt, axisEta, axisDMass});
      histos.add("h3dRecDbar", "h2dRecDbar", kTH3F, {axisPt, axisEta, axisDMass});

      histos.add("hMassD", "hMassD", kTH1F, {axisDMass});
      histos.add("hMassDbar", "hMassDbar", kTH1F, {axisDMass});

      histos.add("hDCosPA", "hDCosPA", kTH1F, {{200, 0, 1}});
      histos.add("hDCosPAxy", "hDCosPAxy", kTH1F, {{200, 0, 1}});
      histos.add("hDCosThetaStar", "hDCosThetaStar", kTH1F, {{200, -1, 1}});
      histos.add("hDDauDecayLength", "hDDauDecayLength", kTH1F, {{100, 0, 10}});

      if (doDCAplotsD) {
        histos.add("hDCADDaughters", "hDCADDaughters", kTH1D, {axisDCADaughters});
        histos.add("hDCADbarDaughters", "hDCADbarDaughters", kTH1D, {axisDCA});
        histos.add("h2dDCAxyVsPtPiPlusFromD", "h2dDCAxyVsPtPiPlusFromD", kTH2F, {axisPt, axisDCA});
        histos.add("h2dDCAxyVsPtPiMinusFromD", "h2dDCAxyVsPtPiMinusFromD", kTH2F, {axisPt, axisDCA});
        histos.add("h2dDCAxyVsPtKaPlusFromD", "h2dDCAxyVsPtKaPlusFromD", kTH2F, {axisPt, axisDCA});
        histos.add("h2dDCAxyVsPtKaMinusFromD", "h2dDCAxyVsPtKaMinusFromD", kTH2F, {axisPt, axisDCA});
      }
    }
    if (doprocessFindLcBaryons) {
      histos.add("h2dGenLc", "h2dGenLc", kTH2F, {axisPt, axisEta});
      histos.add("h2dGenLcbar", "h2dGenLcbar", kTH2F, {axisPt, axisEta});
      histos.add("h3dRecLc", "h2dRecLc", kTH3F, {axisPt, axisEta, axisLcMass});
      histos.add("h3dRecLcbar", "h2dRecLcbar", kTH3F, {axisPt, axisEta, axisLcMass});

      histos.add("hMassLc", "hMassLc", kTH1F, {axisLcMass});
      histos.add("hMassLcbar", "hMassLcbar", kTH1F, {axisLcMass});

      if (doDCAplotsD) {
        histos.add("hDCALcDaughters", "hDCALcDaughters", kTH1D, {axisDCADaughters});
        histos.add("hDCALcbarDaughters", "hDCALcbarDaughters", kTH1D, {axisDCA});
        histos.add("h2dDCAxyVsPtPiPlusFromLc", "h2dDCAxyVsPtPiPlusFromLc", kTH2F, {axisPt, axisDCA});
        histos.add("h2dDCAxyVsPtPiMinusFromLc", "h2dDCAxyVsPtPiMinusFromLc", kTH2F, {axisPt, axisDCA});
        histos.add("h2dDCAxyVsPtKaPlusFromLc", "h2dDCAxyVsPtKaPlusFromLc", kTH2F, {axisPt, axisDCA});
        histos.add("h2dDCAxyVsPtKaMinusFromLc", "h2dDCAxyVsPtKaMinusFromLc", kTH2F, {axisPt, axisDCA});
        histos.add("h2dDCAxyVsPtPrPlusFromLc", "h2dDCAxyVsPtPrPlusFromLc", kTH2F, {axisPt, axisDCA});
        histos.add("h2dDCAxyVsPtPrMinusFromLc", "h2dDCAxyVsPtPrMinusFromLc", kTH2F, {axisPt, axisDCA});
      }
    }
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processGenerated(aod::McParticles const&)
  {
    // no grouping for MC particles -> as intended
    if (doprocessFindDmesons) {
      for (auto const& mcParticle : trueD)
        histos.fill(HIST("h2dGenD"), mcParticle.pt(), mcParticle.eta());
      for (auto const& mcParticle : trueDbar)
        histos.fill(HIST("h2dGenDbar"), mcParticle.pt(), mcParticle.eta());
    }
    if (doprocessFindLcBaryons) {
      for (auto const& mcParticle : trueLc)
        histos.fill(HIST("h2dGenLc"), mcParticle.pt(), mcParticle.eta());
      for (auto const& mcParticle : trueLcbar)
        histos.fill(HIST("h2dGenLcbar"), mcParticle.pt(), mcParticle.eta());
    }
  }

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processFindDmesons(aod::Collision const& collision, alice3tracks const&, aod::McParticles const&)
  {
    // group with this collision
    auto tracksPiPlusFromDgrouped = tracksPiPlusFromD->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracksKaMinusFromDgrouped = tracksKaMinusFromD->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracksKaPlusFromDgrouped = tracksKaPlusFromD->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracksPiMinusFromDgrouped = tracksPiMinusFromD->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (doDCAplotsD) {
      for (auto const& track : tracksPiPlusFromDgrouped)
        histos.fill(HIST("h2dDCAxyVsPtPiPlusFromD"), track.pt(), track.dcaXY() * 1e+4);
      for (auto const& track : tracksPiMinusFromDgrouped)
        histos.fill(HIST("h2dDCAxyVsPtPiMinusFromD"), track.pt(), track.dcaXY() * 1e+4);
      for (auto const& track : tracksKaPlusFromDgrouped)
        histos.fill(HIST("h2dDCAxyVsPtKaPlusFromD"), track.pt(), track.dcaXY() * 1e+4);
      for (auto const& track : tracksKaMinusFromDgrouped)
        histos.fill(HIST("h2dDCAxyVsPtKaMinusFromD"), track.pt(), track.dcaXY() * 1e+4);
    }

    // D mesons
    for (auto const& posTrackRow : tracksPiPlusFromDgrouped) {
      for (auto const& negTrackRow : tracksKaMinusFromDgrouped) {
        if (mcSameMotherCheck && !checkSameMother(posTrackRow, negTrackRow))
          continue;
        if (!buildDecayCandidateTwoBody(posTrackRow, negTrackRow, o2::constants::physics::MassPionCharged, o2::constants::physics::MassKaonCharged))
          continue;

        dmeson.cosPA = RecoDecay::cpa(std::array{collision.posX(), collision.posY(), collision.posZ()}, std::array{dmeson.posSV[0], dmeson.posSV[1], dmeson.posSV[2]}, std::array{dmeson.P[0], dmeson.P[1], dmeson.P[2]});
        dmeson.cosPAxy = RecoDecay::cpaXY(std::array{collision.posX(), collision.posY(), collision.posZ()}, std::array{dmeson.posSV[0], dmeson.posSV[1], dmeson.posSV[2]}, std::array{dmeson.P[0], dmeson.P[1], dmeson.P[2]});

        const float dmesonCtau = 0.012301;
        dmeson.normalizedDecayLength = ((dmeson.mass * std::fabs(std::hypot(collision.posX(), collision.posY(), collision.posZ()) - std::hypot(dmeson.posSV[0], dmeson.posSV[1], dmeson.posSV[2]))) / std::hypot(dmeson.P[0], dmeson.P[1], dmeson.P[2])) / dmesonCtau;

        histos.fill(HIST("hDCosPA"), dmeson.cosPA);
        histos.fill(HIST("hDCosPAxy"), dmeson.cosPAxy);
        histos.fill(HIST("hDCosThetaStar"), dmeson.cosThetaStar);
        histos.fill(HIST("hDDauDecayLength"), dmeson.normalizedDecayLength);

        if (dmeson.dcaDau > dcaDaughtersSelection)
          continue;

        if (dmeson.pt <= 3 && dmeson.cosPA < DCosPA)
          continue;
        else if (dmeson.pt > 3 && dmeson.cosPA < DCosPAHighPt)
          continue;

        if (dmeson.pt <= 3 && dmeson.cosPAxy < DCosPAxy)
          continue;
        else if (dmeson.pt > 3 && dmeson.cosPAxy < DCosPAxyHighPt)
          continue;

        if (dmeson.pt <= 9 && std::fabs(dmeson.cosThetaStar) > DCosThetaStarLowPt)
          continue;
        else if (dmeson.pt <= 16 && std::fabs(dmeson.cosThetaStar) > DCosThetaStarHighPt)
          continue;
        else if (dmeson.pt > 16 && std::fabs(dmeson.cosThetaStar) > DCosThetaStarVHighPt)
          continue;

        if (dmeson.normalizedDecayLength > DDauDecayLength)
          continue;

        histos.fill(HIST("hDCADDaughters"), dmeson.dcaDau * 1e+4);
        histos.fill(HIST("hMassD"), dmeson.mass);
        histos.fill(HIST("h3dRecD"), dmeson.pt, dmeson.eta, dmeson.mass);
      }
    }
    // D mesons
    for (auto const& posTrackRow : tracksKaPlusFromDgrouped) {
      for (auto const& negTrackRow : tracksPiMinusFromDgrouped) {
        if (mcSameMotherCheck && !checkSameMother(posTrackRow, negTrackRow))
          continue;
        if (!buildDecayCandidateTwoBody(posTrackRow, negTrackRow, o2::constants::physics::MassKaonCharged, o2::constants::physics::MassPionCharged))
          continue;

        dmeson.cosPA = RecoDecay::cpa(std::array{collision.posX(), collision.posY(), collision.posZ()}, std::array{dmeson.posSV[0], dmeson.posSV[1], dmeson.posSV[2]}, std::array{dmeson.P[0], dmeson.P[1], dmeson.P[2]});
        dmeson.cosPAxy = RecoDecay::cpaXY(std::array{collision.posX(), collision.posY(), collision.posZ()}, std::array{dmeson.posSV[0], dmeson.posSV[1], dmeson.posSV[2]}, std::array{dmeson.P[0], dmeson.P[1], dmeson.P[2]});

        const float dmesonCtau = 0.012301;
        dmeson.normalizedDecayLength = ((dmeson.mass * std::fabs(std::hypot(collision.posX(), collision.posY(), collision.posZ()) - std::hypot(dmeson.posSV[0], dmeson.posSV[1], dmeson.posSV[2]))) / std::hypot(dmeson.P[0], dmeson.P[1], dmeson.P[2])) / dmesonCtau;

        histos.fill(HIST("hDCosPA"), dmeson.cosPA);
        histos.fill(HIST("hDCosPAxy"), dmeson.cosPAxy);
        histos.fill(HIST("hDCosThetaStar"), dmeson.cosThetaStar);
        histos.fill(HIST("hDDauDecayLength"), dmeson.normalizedDecayLength);

        if (dmeson.dcaDau > dcaDaughtersSelection)
          continue;

        if (dmeson.pt <= 3 && dmeson.cosPA < DCosPA)
          continue;
        else if (dmeson.pt > 3 && dmeson.cosPA < DCosPAHighPt)
          continue;

        if (dmeson.pt <= 3 && dmeson.cosPAxy < DCosPAxy)
          continue;
        else if (dmeson.pt > 3 && dmeson.cosPAxy < DCosPAxyHighPt)
          continue;

        if (dmeson.pt <= 9 && std::fabs(dmeson.cosThetaStar) > DCosThetaStarLowPt)
          continue;
        else if (dmeson.pt <= 16 && std::fabs(dmeson.cosThetaStar) > DCosThetaStarHighPt)
          continue;
        else if (dmeson.pt > 16 && std::fabs(dmeson.cosThetaStar) > DCosThetaStarVHighPt)
          continue;

        if (dmeson.normalizedDecayLength > DDauDecayLength)
          continue;

        histos.fill(HIST("hDCADbarDaughters"), dmeson.dcaDau * 1e+4);
        histos.fill(HIST("hMassDbar"), dmeson.mass);
        histos.fill(HIST("h3dRecDbar"), dmeson.pt, dmeson.eta, dmeson.mass);
      }
    }
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*
  void processFindLcBaryons(aod::Collision const& collision, alice3tracks const&, aod::McParticles const&)
  {
    // group with this collision
    auto tracksPiPlusFromLcgrouped = tracksPiPlusFromLc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracksKaPlusFromLcgrouped = tracksKaPlusFromLc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracksPrPlusFromLcgrouped = tracksPrPlusFromLc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    auto tracksPiMinusFromLcgrouped = tracksPiMinusFromLc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracksKaMinusFromLcgrouped = tracksKaMinusFromLc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto tracksPrMinusFromLcgrouped = tracksPrMinusFromLc->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

    if (doDCAplotsLc) {
      for (auto const& track : tracksPiPlusFromLcgrouped)
        histos.fill(HIST("h2dDCAxyVsPtPiPlusFromLc"), track.pt(), track.dcaXY() * 1e+4);
      for (auto const& track : tracksPiMinusFromLcgrouped)
        histos.fill(HIST("h2dDCAxyVsPtPiMinusFromLc"), track.pt(), track.dcaXY() * 1e+4);
      for (auto const& track : tracksKaPlusFromLcgrouped)
        histos.fill(HIST("h2dDCAxyVsPtKaPlusFromLc"), track.pt(), track.dcaXY() * 1e+4);
      for (auto const& track : tracksKaMinusFromLcgrouped)
        histos.fill(HIST("h2dDCAxyVsPtKaMinusFromLc"), track.pt(), track.dcaXY() * 1e+4);
      for (auto const& track : tracksPrPlusFromLcgrouped)
        histos.fill(HIST("h2dDCAxyVsPtPrPlusFromLc"), track.pt(), track.dcaXY() * 1e+4);
      for (auto const& track : tracksPrMinusFromLcgrouped)
        histos.fill(HIST("h2dDCAxyVsPtPrMinusFromLc"), track.pt(), track.dcaXY() * 1e+4);
    }

    // Lc+ baryons +4122 -> +2212 -321 +211
    for (auto const& proton : tracksPrPlusFromLcgrouped) {
      for (auto const& pion : tracksPiPlusFromLcgrouped) {
        if (pion.globalIndex() == proton.globalIndex())
          continue; // avoid self
        for (auto const& kaon : tracksKaMinusFromLcgrouped) {
          if (mcSameMotherCheck && (!checkSameMother(proton, kaon) || !checkSameMother(proton, pion)))
            continue;
          if (!buildDecayCandidateThreeBody(proton, kaon, pion, o2::constants::physics::MassProton, o2::constants::physics::MassKaonCharged, o2::constants::physics::MassPionCharged))
            continue;
          histos.fill(HIST("hDCALcDaughters"), lcbaryon.dcaDau * 1e+4);
          histos.fill(HIST("hMassLc"), lcbaryon.mass);
          histos.fill(HIST("h3dRecLc"), lcbaryon.pt, lcbaryon.eta, lcbaryon.mass);
        }
      }
    }
    // Lc- baryons -4122 -> -2212 +321 -211
    for (auto const& proton : tracksPrMinusFromLcgrouped) {
      for (auto const& pion : tracksPiMinusFromLcgrouped) {
        if (pion.globalIndex() == proton.globalIndex())
          continue; // avoid self
        for (auto const& kaon : tracksKaPlusFromLcgrouped) {
          if (mcSameMotherCheck && (!checkSameMother(proton, kaon) || !checkSameMother(proton, pion)))
            continue;
          if (!buildDecayCandidateThreeBody(proton, kaon, pion, o2::constants::physics::MassProton, o2::constants::physics::MassKaonCharged, o2::constants::physics::MassPionCharged))
            continue;
          histos.fill(HIST("hDCALcbarDaughters"), lcbaryon.dcaDau * 1e+4);
          histos.fill(HIST("hMassLcbar"), lcbaryon.mass);
          histos.fill(HIST("h3dRecLcbar"), lcbaryon.pt, lcbaryon.eta, lcbaryon.mass);
        }
      }
    }
  }
  //*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*+-+*

  //*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*
  PROCESS_SWITCH(alice3decayFinder, processGenerated, "fill MC-only histograms", true);
  PROCESS_SWITCH(alice3decayFinder, processFindDmesons, "find D mesons", true);
  PROCESS_SWITCH(alice3decayFinder, processFindLcBaryons, "find Lc Baryons", true);
  //*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*>-~-<*
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<alice3decayFinder>(cfgc)};
}
