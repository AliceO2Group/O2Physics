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
/// \file hypKfTreeCreator.cxx
/// \brief Creates flat tree for ML analysis
/// \author Janik Ditzel <jditzel@cern.ch> and Michael Hartung <mhartung@cern.ch>

#include "PWGLF/DataModel/LFHypernucleiKfTables.h"

#include "Common/Core/RecoDecay.h"

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include <vector>

using namespace o2;
using namespace o2::framework;
typedef std::array<float, 3> arr3;

namespace
{
enum Decays { kTwoBody = 2,
              kThreeBody = 3 };

struct TrackProperties {
  TrackProperties() : x(0), y(0), z(0), px(0), py(0), pz(0), tpcNcls(0), itsNcls(0), tpcChi2(0), itsChi2(0), itsMeanClsSizeL(0), rigidity(0), tpcSignal(0), tpcNsigma(0), tpcNsigmaNhp(0), tpcNsigmaNlp(0), tofMass(0), dcaXY(0), dcaZ(0), isPvContributor(0), subMass(0) {}
  float x, y, z, px, py, pz;
  uint8_t tpcNcls, itsNcls;
  float tpcChi2, itsChi2, itsMeanClsSizeL;
  float rigidity, tpcSignal, tpcNsigma, tpcNsigmaNhp, tpcNsigmaNlp;
  float tofMass, dcaXY, dcaZ;
  bool isPvContributor;
  float subMass;
};

struct HyperNucleus {
  HyperNucleus() : pdgCode(0), isReconstructed(0), globalIndex(0), species(0), speciesMC(0), isMatter(0), passedEvSel(0), isMatterMC(0), passedEvSelMC(0), isPhysicalPrimary(0), collisionMcTrue(0), mass(0), y(0), pt(0), ct(0), yGen(0), ptGen(0), ctGen(0), cpaPvGen(0), cpaPv(0), cpaSv(0), maxDcaTracks(0), maxDcaTracksSv(0), dcaToPvXY(0), dcaToPvZ(0), dcaToVtxXY(0), dcaToVtxZ(0), devToPvXY(0), chi2(0), pvx(0), pvy(0), pvz(0), svx(0), svy(0), svz(0), px(0), py(0), pz(0), pvxGen(0), pvyGen(0), pvzGen(0), svxGen(0), svyGen(0), svzGen(0), pxGen(0), pyGen(0), pzGen(0), nSingleDaughters(0), mcTrue(0), mcTrueVtx(0), mcPhysicalPrimary(0) {}
  int pdgCode, isReconstructed, globalIndex;
  uint8_t species, speciesMC;
  bool isMatter, passedEvSel, isMatterMC, passedEvSelMC, isPhysicalPrimary, collisionMcTrue;
  float mass, y, pt, ct, yGen, ptGen, ctGen, cpaPvGen, cpaPv, cpaSv, maxDcaTracks, maxDcaTracksSv;
  float dcaToPvXY, dcaToPvZ, dcaToVtxXY, dcaToVtxZ, devToPvXY, chi2;
  float pvx, pvy, pvz, svx, svy, svz, px, py, pz;
  float pvxGen, pvyGen, pvzGen, svxGen, svyGen, svzGen, pxGen, pyGen, pzGen;
  int nSingleDaughters, cent, occu, runNumber;
  bool mcTrue, mcTrueVtx, mcPhysicalPrimary;
  std::vector<TrackProperties> daughterTracks;
  std::vector<float> subDaughterMassVec;
};
} // namespace
namespace o2::aod
{
namespace hypkftree
{
DECLARE_SOA_COLUMN(Y, y, float);
DECLARE_SOA_COLUMN(Ct, ct, float);
DECLARE_SOA_COLUMN(YGen, yGen, float);
DECLARE_SOA_COLUMN(PtGen, ptGen, float);
DECLARE_SOA_COLUMN(CtGen, ctGen, float);
DECLARE_SOA_COLUMN(CpaPvGen, cpaPvGen, float);
DECLARE_SOA_COLUMN(DcaTracks, dcaTracks, float);
DECLARE_SOA_COLUMN(DcaTrackSv, dcaTrackSv, float);
DECLARE_SOA_COLUMN(CosPa, cosPa, double);
DECLARE_SOA_COLUMN(McTrue, mcTrue, bool);
DECLARE_SOA_COLUMN(Pvx, pvx, float);
DECLARE_SOA_COLUMN(Pvy, pvy, float);
DECLARE_SOA_COLUMN(Pvz, pvz, float);
DECLARE_SOA_COLUMN(Tvx, tvx, float);
DECLARE_SOA_COLUMN(Tvy, tvy, float);
DECLARE_SOA_COLUMN(Tvz, tvz, float);
DECLARE_SOA_COLUMN(PxGen, pxGen, float);
DECLARE_SOA_COLUMN(PyGen, pyGen, float);
DECLARE_SOA_COLUMN(PzGen, pzGen, float);
DECLARE_SOA_COLUMN(PvxGen, pvxGen, float);
DECLARE_SOA_COLUMN(PvyGen, pvyGen, float);
DECLARE_SOA_COLUMN(PvzGen, pvzGen, float);
DECLARE_SOA_COLUMN(SvxGen, svxGen, float);
DECLARE_SOA_COLUMN(SvyGen, svyGen, float);
DECLARE_SOA_COLUMN(SvzGen, svzGen, float);
DECLARE_SOA_COLUMN(TvxGen, tvxGen, float);
DECLARE_SOA_COLUMN(TvyGen, tvyGen, float);
DECLARE_SOA_COLUMN(TvzGen, tvzGen, float);
DECLARE_SOA_COLUMN(Centrality, centrality, int);
DECLARE_SOA_COLUMN(Occupancy, occupancy, int);
DECLARE_SOA_COLUMN(RunNumber, runNumber, int);
DECLARE_SOA_COLUMN(PassedEvSelMC, passedEvSelMC, bool);
DECLARE_SOA_COLUMN(SpeciesMC, speciesMC, int8_t); //!
DECLARE_SOA_COLUMN(IsMatter, isMatter, bool);
DECLARE_SOA_COLUMN(IsMatterGen, isMatterGen, bool);
DECLARE_SOA_COLUMN(IsReconstructed, isReconstructed, int);
DECLARE_SOA_COLUMN(CollMcTrue, collMcTrue, bool);
DECLARE_SOA_COLUMN(D1X, d1X, float);
DECLARE_SOA_COLUMN(D1Y, d1Y, float);
DECLARE_SOA_COLUMN(D1Z, d1Z, float);
DECLARE_SOA_COLUMN(D1Px, d1Px, float);
DECLARE_SOA_COLUMN(D1Py, d1Py, float);
DECLARE_SOA_COLUMN(D1Pz, d1Pz, float);
DECLARE_SOA_COLUMN(D1TPCnCls, d1TPCnCls, uint8_t);
DECLARE_SOA_COLUMN(D1TPCchi2, d1TPCchi2, float);
DECLARE_SOA_COLUMN(D1ITSnCls, d1ITSnCls, uint8_t);
DECLARE_SOA_COLUMN(D1ITSchi2, d1ITSchi2, float);
DECLARE_SOA_COLUMN(D1ITSmeanClsSizeL, d1ITSmeanClsSizeL, float);
DECLARE_SOA_COLUMN(D1Rigidity, d1Rigidity, float);
DECLARE_SOA_COLUMN(D1TPCsignal, d1TPCsignal, float);
DECLARE_SOA_COLUMN(D1TPCnSigma, d1TPCnSigma, float);
DECLARE_SOA_COLUMN(D1TPCnSigmaNhp, d1TPCnSigmaNhp, float);
DECLARE_SOA_COLUMN(D1TPCnSigmaNlp, d1TPCnSigmaNlp, float);
DECLARE_SOA_COLUMN(D1TOFmass, d1TOFmass, float);
DECLARE_SOA_COLUMN(D1DcaXY, d1DcaXY, float);
DECLARE_SOA_COLUMN(D1DcaZ, d1DcaZ, float);
DECLARE_SOA_COLUMN(D1IsPvContributor, d1IsPvContributor, bool);
DECLARE_SOA_COLUMN(D2X, d2X, float);
DECLARE_SOA_COLUMN(D2Y, d2Y, float);
DECLARE_SOA_COLUMN(D2Z, d2Z, float);
DECLARE_SOA_COLUMN(D2Px, d2Px, float);
DECLARE_SOA_COLUMN(D2Py, d2Py, float);
DECLARE_SOA_COLUMN(D2Pz, d2Pz, float);
DECLARE_SOA_COLUMN(D2TPCnCls, d2TPCnCls, uint8_t);
DECLARE_SOA_COLUMN(D2TPCchi2, d2TPCchi2, float);
DECLARE_SOA_COLUMN(D2ITSnCls, d2ITSnCls, uint8_t);
DECLARE_SOA_COLUMN(D2ITSchi2, d2ITSchi2, float);
DECLARE_SOA_COLUMN(D2ITSmeanClsSizeL, d2ITSmeanClsSizeL, float);
DECLARE_SOA_COLUMN(D2Rigidity, d2Rigidity, float);
DECLARE_SOA_COLUMN(D2TPCsignal, d2TPCsignal, float);
DECLARE_SOA_COLUMN(D2TPCnSigma, d2TPCnSigma, float);
DECLARE_SOA_COLUMN(D2TPCnSigmaNhp, d2TPCnSigmaNhp, float);
DECLARE_SOA_COLUMN(D2TPCnSigmaNlp, d2TPCnSigmaNlp, float);
DECLARE_SOA_COLUMN(D2TOFmass, d2TOFmass, float);
DECLARE_SOA_COLUMN(D2DcaXY, d2DcaXY, float);
DECLARE_SOA_COLUMN(D2DcaZ, d2DcaZ, float);
DECLARE_SOA_COLUMN(D2IsPvContributor, d2IsPvContributor, bool);
DECLARE_SOA_COLUMN(D3X, d3X, float);
DECLARE_SOA_COLUMN(D3Y, d3Y, float);
DECLARE_SOA_COLUMN(D3Z, d3Z, float);
DECLARE_SOA_COLUMN(D3Px, d3Px, float);
DECLARE_SOA_COLUMN(D3Py, d3Py, float);
DECLARE_SOA_COLUMN(D3Pz, d3Pz, float);
DECLARE_SOA_COLUMN(D3TPCnCls, d3TPCnCls, uint8_t);
DECLARE_SOA_COLUMN(D3TPCchi2, d3TPCchi2, float);
DECLARE_SOA_COLUMN(D3ITSnCls, d3ITSnCls, uint8_t);
DECLARE_SOA_COLUMN(D3ITSchi2, d3ITSchi2, float);
DECLARE_SOA_COLUMN(D3ITSmeanClsSizeL, d3ITSmeanClsSizeL, float);
DECLARE_SOA_COLUMN(D3Rigidity, d3Rigidity, float);
DECLARE_SOA_COLUMN(D3TPCsignal, d3TPCsignal, float);
DECLARE_SOA_COLUMN(D3TPCnSigma, d3TPCnSigma, float);
DECLARE_SOA_COLUMN(D3TPCnSigmaNhp, d3TPCnSigmaNhp, float);
DECLARE_SOA_COLUMN(D3TPCnSigmaNlp, d3TPCnSigmaNlp, float);
DECLARE_SOA_COLUMN(D3TOFmass, d3TOFmass, float);
DECLARE_SOA_COLUMN(D3DcaXY, d3DcaXY, float);
DECLARE_SOA_COLUMN(D3DcaZ, d3DcaZ, float);
DECLARE_SOA_COLUMN(D1d2Mass, d1d2Mass, float);
DECLARE_SOA_COLUMN(D1d3Mass, d1d3Mass, float);
DECLARE_SOA_COLUMN(D2d3Mass, d2d3Mass, float);
DECLARE_SOA_COLUMN(D3IsPvContributor, d3IsPvContributor, bool);
} // namespace hypkftree

#define HYPKFGENBASE hypkftree::SpeciesMC, mcparticle::PdgCode, hypkftree::IsMatterGen, hypkftree::IsReconstructed, hykfmc::IsPhysicalPrimary, hypkftree::PassedEvSelMC, hypkftree::YGen, hypkftree::PtGen, hypkftree::CtGen

#define HYPKFGENEXT hypkftree::CpaPvGen, hypkftree::PxGen, hypkftree::PyGen, hypkftree::PzGen, hypkftree::PvxGen, hypkftree::PvyGen, hypkftree::PvzGen, hypkftree::SvxGen, hypkftree::SvyGen, hypkftree::SvzGen

#define HYPKFHYPNUC hykfmc::Species, hypkftree::IsMatter, hypkftree::Centrality, hypkftree::Occupancy, hypkftree::RunNumber, hykfmccoll::PassedEvSel, hykfhyp::Mass, hypkftree::Y, track::Pt, hypkftree::Ct, hypkftree::CosPa, hypkftree::DcaTracks, hypkftree::DcaTrackSv, hykfhyp::DcaToPvXY, hykfhyp::DcaToPvZ, hykfhyp::DevToPvXY, hykfhyp::Chi2, hypkftree::Pvx, hypkftree::Pvy, hypkftree::Pvz, hykfmc::Svx, hykfmc::Svy, hykfmc::Svz, hykfhyp::Px, hykfhyp::Py, hykfhyp::Pz, hypkftree::CollMcTrue

#define HYPKFHYPNUCMC hypkftree::McTrue, hykfmc::IsPhysicalPrimary

#define HYPKFD1 hypkftree::D1X, hypkftree::D1Y, hypkftree::D1Z, hypkftree::D1Px, hypkftree::D1Py, hypkftree::D1Pz, hypkftree::D1TPCnCls, hypkftree::D1TPCchi2, hypkftree::D1ITSnCls, hypkftree::D1ITSchi2, hypkftree::D1ITSmeanClsSizeL, hypkftree::D1Rigidity, hypkftree::D1TPCsignal, hypkftree::D1TPCnSigma, hypkftree::D1TPCnSigmaNhp, hypkftree::D1TPCnSigmaNlp, hypkftree::D1TOFmass, hypkftree::D1DcaXY, hypkftree::D1DcaZ, hypkftree::D1IsPvContributor

#define HYPKFD2 hypkftree::D2X, hypkftree::D2Y, hypkftree::D2Z, hypkftree::D2Px, hypkftree::D2Py, hypkftree::D2Pz, hypkftree::D2TPCnCls, hypkftree::D2TPCchi2, hypkftree::D2ITSnCls, hypkftree::D2ITSchi2, hypkftree::D2ITSmeanClsSizeL, hypkftree::D2Rigidity, hypkftree::D2TPCsignal, hypkftree::D2TPCnSigma, hypkftree::D2TPCnSigmaNhp, hypkftree::D2TPCnSigmaNlp, hypkftree::D2TOFmass, hypkftree::D2DcaXY, hypkftree::D2DcaZ, hypkftree::D2IsPvContributor

#define HYPKFD3 hypkftree::D3X, hypkftree::D3Y, hypkftree::D3Z, hypkftree::D3Px, hypkftree::D3Py, hypkftree::D3Pz, hypkftree::D3TPCnCls, hypkftree::D3TPCchi2, hypkftree::D3ITSnCls, hypkftree::D3ITSchi2, hypkftree::D3ITSmeanClsSizeL, hypkftree::D3Rigidity, hypkftree::D3TPCsignal, hypkftree::D3TPCnSigma, hypkftree::D3TPCnSigmaNhp, hypkftree::D3TPCnSigmaNlp, hypkftree::D3TOFmass, hypkftree::D3DcaXY, hypkftree::D3DcaZ, hypkftree::D3IsPvContributor

#define HYPKFSDMASS hypkftree::D1d2Mass, hypkftree::D1d3Mass, hypkftree::D2d3Mass

DECLARE_SOA_TABLE(HypKfGens, "AOD", "HYPKFGEN", HYPKFGENBASE);
using HypKfGen = HypKfGens::iterator;

DECLARE_SOA_TABLE(HypKfSingleTwoBodyCandidates, "AOD", "HYPKFCAND2", HYPKFHYPNUC, HYPKFHYPNUCMC, HYPKFD1, HYPKFD2);
using HypKfSingleTwoBodyCandidate = HypKfSingleTwoBodyCandidates::iterator;

DECLARE_SOA_TABLE(HypKfMcSingleTwoBodyCandidates, "AOD", "HYPKFMCCAND2", HYPKFGENBASE, HYPKFGENEXT, HYPKFHYPNUC, HYPKFD1, HYPKFD2);
using HypKfMcSingleTwoBodyCandidate = HypKfMcSingleTwoBodyCandidates::iterator;

DECLARE_SOA_TABLE(HypKfSingleThreeBodyCandidates, "AOD", "HYPKFCAND3", HYPKFHYPNUC, HYPKFHYPNUCMC, HYPKFD1, HYPKFD2, HYPKFD3, HYPKFSDMASS);
using HypKfSingleThreeBodyCandidate = HypKfSingleThreeBodyCandidates::iterator;

DECLARE_SOA_TABLE(HypKfMcSingleThreeBodyCandidates, "AOD", "HYPKFMCCAND3", HYPKFGENBASE, HYPKFGENEXT, HYPKFHYPNUC, HYPKFD1, HYPKFD2, HYPKFD3, HYPKFSDMASS);
using HypKfMcSingleThreeBodyCandidate = HypKfMcSingleThreeBodyCandidates::iterator;

} // namespace o2::aod

struct HypKfTreeCreator {

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Produces<aod::HypKfGens> outputMcGenTable;
  Produces<aod::HypKfSingleTwoBodyCandidates> outputTableTwo;
  Produces<aod::HypKfMcSingleTwoBodyCandidates> outputTableMcTwo;
  Produces<aod::HypKfSingleThreeBodyCandidates> outputTableThree;
  Produces<aod::HypKfMcSingleThreeBodyCandidates> outputTableMcThree;
  PresliceUnsorted<aod::HypKfHypNucs> perMcParticle = aod::hykfhyp::hypKfMcPartId;

  Configurable<int> cfgSpecies{"cfgSpecies", 0, "Select species"};
  Configurable<int> cfgNprimDaughters{"cfgNprimDaughters", 0, "Number of primary daughters"};
  Configurable<bool> cfgMCGenerated{"cfgMCGenerated", false, "create MC generated tree"};
  Configurable<bool> cfgMCReconstructed{"cfgMCReconstructed", false, "create MC reconstructed tree"};
  Configurable<bool> cfgMCCombined{"cfgMCCombined", false, "create MC tree containig generated and reconstructed"};

  bool isMC;
  //___________________________________________________________________________________________________________________________________________________________

  void init(InitContext const&)
  {
    isMC = false;
  }

  //___________________________________________________________________________________________________________________________________________________________
  void fillTable(HyperNucleus& cand)
  {
    if (isMC && cfgMCGenerated)
      outputMcGenTable(
        cand.speciesMC, cand.pdgCode, cand.isMatterMC, cand.isReconstructed, cand.isPhysicalPrimary, cand.passedEvSelMC, cand.yGen, cand.ptGen, cand.ctGen);

    if (!cand.isReconstructed) {
      cand.daughterTracks.resize(4);
      cand.subDaughterMassVec.resize(4);
    }

    if (cfgNprimDaughters == Decays::kTwoBody) {
      const auto& d1 = cand.daughterTracks.at(0);
      const auto& d2 = cand.daughterTracks.at(1);
      if (!isMC || (isMC && cfgMCReconstructed && cand.isReconstructed))
        outputTableTwo(
          cand.species, cand.isMatter, cand.cent, cand.occu, cand.runNumber, cand.passedEvSel, cand.mass, cand.y, cand.pt, cand.ct, cand.cpaPv, cand.maxDcaTracks, cand.maxDcaTracksSv,
          cand.dcaToPvXY, cand.dcaToPvZ, cand.devToPvXY, cand.chi2, cand.pvx, cand.pvy, cand.pvz, cand.svx, cand.svy, cand.svz, cand.px, cand.py, cand.pz, cand.collisionMcTrue,
          cand.mcTrue, cand.mcPhysicalPrimary, d1.x, d1.y, d1.z, d1.px, d1.py, d1.pz, d1.tpcNcls, d1.tpcChi2, d1.itsNcls, d1.itsChi2, d1.itsMeanClsSizeL,
          d1.rigidity, d1.tpcSignal, d1.tpcNsigma, d1.tpcNsigmaNhp, d1.tpcNsigmaNlp, d1.tofMass, d1.dcaXY, d1.dcaZ, d1.isPvContributor,
          d2.x, d2.y, d2.z, d2.px, d2.py, d2.pz, d2.tpcNcls, d2.tpcChi2, d2.itsNcls, d2.itsChi2, d2.itsMeanClsSizeL,
          d2.rigidity, d2.tpcSignal, d2.tpcNsigma, d2.tpcNsigmaNhp, d2.tpcNsigmaNlp, d2.tofMass, d2.dcaXY, d2.dcaZ, d2.isPvContributor);
      if (isMC && cfgMCCombined)
        outputTableMcTwo(
          cand.speciesMC, cand.pdgCode, cand.isMatterMC, cand.isReconstructed, cand.isPhysicalPrimary, cand.passedEvSelMC, cand.yGen, cand.ptGen, cand.ctGen,
          cand.cpaPvGen, cand.pxGen, cand.pyGen, cand.pzGen, cand.pvxGen, cand.pvyGen, cand.pvzGen, cand.svxGen, cand.svyGen, cand.svzGen,
          cand.species, cand.isMatter, cand.cent, cand.occu, cand.runNumber, cand.passedEvSel, cand.mass, cand.y, cand.pt, cand.ct, cand.cpaPv, cand.maxDcaTracks,
          cand.maxDcaTracksSv, cand.dcaToPvXY, cand.dcaToPvZ, cand.devToPvXY,
          cand.chi2, cand.pvx, cand.pvy, cand.pvz, cand.svx, cand.svy, cand.svz, cand.px, cand.py, cand.pz, cand.collisionMcTrue,
          d1.x, d1.y, d1.z, d1.px, d1.py, d1.pz, d1.tpcNcls, d1.tpcChi2, d1.itsNcls, d1.itsChi2, d1.itsMeanClsSizeL,
          d1.rigidity, d1.tpcSignal, d1.tpcNsigma, d1.tpcNsigmaNhp, d1.tpcNsigmaNlp, d1.tofMass, d1.dcaXY, d1.dcaZ, d1.isPvContributor,
          d2.x, d2.y, d2.z, d2.px, d2.py, d2.pz, d2.tpcNcls, d2.tpcChi2, d2.itsNcls, d2.itsChi2, d2.itsMeanClsSizeL,
          d2.rigidity, d2.tpcSignal, d2.tpcNsigma, d2.tpcNsigmaNhp, d2.tpcNsigmaNlp, d2.tofMass, d2.dcaXY, d2.dcaZ, d2.isPvContributor);
    }
    if (cfgNprimDaughters == Decays::kThreeBody) {
      const auto& d1 = cand.daughterTracks.at(0);
      const auto& d2 = cand.daughterTracks.at(1);
      const auto& d3 = cand.daughterTracks.at(2);
      if (!isMC || (isMC && cfgMCReconstructed && cand.isReconstructed))
        outputTableThree(
          cand.species, cand.isMatter, cand.cent, cand.occu, cand.runNumber, cand.passedEvSel, cand.mass, cand.y, cand.pt, cand.ct, cand.cpaPv, cand.maxDcaTracks, cand.maxDcaTracksSv,
          cand.dcaToPvXY, cand.dcaToPvZ, cand.devToPvXY, cand.chi2, cand.pvx, cand.pvy, cand.pvz, cand.svx, cand.svy, cand.svz, cand.px, cand.py, cand.pz, cand.collisionMcTrue,
          cand.mcTrue, cand.mcPhysicalPrimary,
          d1.x, d1.y, d1.z, d1.px, d1.py, d1.pz, d1.tpcNcls, d1.tpcChi2, d1.itsNcls, d1.itsChi2, d1.itsMeanClsSizeL,
          d1.rigidity, d1.tpcSignal, d1.tpcNsigma, d1.tpcNsigmaNhp, d1.tpcNsigmaNlp, d1.tofMass, d1.dcaXY, d1.dcaZ, d1.isPvContributor,
          d2.x, d2.y, d2.z, d2.px, d2.py, d2.pz, d2.tpcNcls, d2.tpcChi2, d2.itsNcls, d2.itsChi2, d2.itsMeanClsSizeL,
          d2.rigidity, d2.tpcSignal, d2.tpcNsigma, d2.tpcNsigmaNhp, d2.tpcNsigmaNlp, d2.tofMass, d2.dcaXY, d2.dcaZ, d2.isPvContributor,
          d3.x, d3.y, d3.z, d3.px, d3.py, d3.pz, d3.tpcNcls, d3.tpcChi2, d3.itsNcls, d3.itsChi2, d3.itsMeanClsSizeL,
          d3.rigidity, d3.tpcSignal, d3.tpcNsigma, d3.tpcNsigmaNhp, d3.tpcNsigmaNlp, d3.tofMass, d3.dcaXY, d3.dcaZ, d3.isPvContributor,
          d1.subMass, d2.subMass, d3.subMass);
      if (isMC && cfgMCCombined)
        outputTableMcThree(
          cand.speciesMC, cand.pdgCode, cand.isMatterMC, cand.isReconstructed, cand.isPhysicalPrimary, cand.passedEvSelMC, cand.yGen, cand.ptGen, cand.ctGen,
          cand.cpaPvGen, cand.pxGen, cand.pyGen, cand.pzGen, cand.pvxGen, cand.pvyGen, cand.pvzGen, cand.svxGen, cand.svyGen, cand.svzGen,
          cand.species, cand.isMatter, cand.cent, cand.occu, cand.runNumber, cand.passedEvSel, cand.mass, cand.y, cand.pt, cand.ct, cand.cpaPv, cand.maxDcaTracks,
          cand.maxDcaTracksSv, cand.dcaToPvXY, cand.dcaToPvZ, cand.devToPvXY, cand.chi2, cand.pvx, cand.pvy, cand.pvz, cand.svx, cand.svy, cand.svz, cand.px, cand.py,
          cand.pz, cand.collisionMcTrue,
          d1.x, d1.y, d1.z, d1.px, d1.py, d1.pz, d1.tpcNcls, d1.tpcChi2, d1.itsNcls, d1.itsChi2, d1.itsMeanClsSizeL,
          d1.rigidity, d1.tpcSignal, d1.tpcNsigma, d1.tpcNsigmaNhp, d1.tpcNsigmaNlp, d1.tofMass, d1.dcaXY, d1.dcaZ, d1.isPvContributor,
          d2.x, d2.y, d2.z, d2.px, d2.py, d2.pz, d2.tpcNcls, d2.tpcChi2, d2.itsNcls, d2.itsChi2, d2.itsMeanClsSizeL,
          d2.rigidity, d2.tpcSignal, d2.tpcNsigma, d2.tpcNsigmaNhp, d2.tpcNsigmaNlp, d2.tofMass, d2.dcaXY, d2.dcaZ, d2.isPvContributor,
          d3.x, d3.y, d3.z, d3.px, d3.py, d3.pz, d3.tpcNcls, d3.tpcChi2, d3.itsNcls, d3.itsChi2, d3.itsMeanClsSizeL,
          d3.rigidity, d3.tpcSignal, d3.tpcNsigma, d3.tpcNsigmaNhp, d3.tpcNsigmaNlp, d3.tofMass, d3.dcaXY, d3.dcaZ, d3.isPvContributor,
          d1.subMass, d2.subMass, d3.subMass);
    }
  }
  //___________________________________________________________________________________________________________________________________________________________
  void fillCandidate(HyperNucleus& cand, aod::HypKfHypNuc const& hypNuc, aod::HypKfHypNucs const&, aod::HypKfColls const&, aod::HypKfTracks const&, aod::HypKfDaughtAdds const&, aod::HypKfSubDs const&)
  {
    cand.daughterTracks.clear();
    cand.subDaughterMassVec.clear();
    auto coll = hypNuc.hypKfColl();
    auto addOns = hypNuc.hypKfDaughtAdd_as<aod::HypKfDaughtAdds>();
    auto posVec = posVector(addOns);
    cand.species = std::abs(hypNuc.species());
    cand.isMatter = hypNuc.isMatter();
    cand.mcTrue = hypNuc.mcTrue();
    cand.cent = coll.centFT0C();
    cand.occu = coll.occupancy();
    cand.runNumber = coll.runNumber();
    cand.passedEvSel = coll.passedEvSel();
    cand.mass = hypNuc.mass();
    cand.y = hypNuc.y();
    cand.pt = hypNuc.pt();
    cand.maxDcaTracks = maxValue(dcaTracksAll(posVec, "XY"));
    cand.maxDcaTracksSv = maxValue(dcaTrackSvAll(posVec, hypNuc, "XY"));
    cand.dcaToPvXY = hypNuc.dcaToPvXY();
    cand.dcaToPvZ = hypNuc.dcaToPvZ();
    cand.dcaToVtxXY = hypNuc.dcaToVtxXY();
    cand.dcaToVtxZ = hypNuc.dcaToVtxZ();
    cand.devToPvXY = hypNuc.devToPvXY();
    cand.chi2 = hypNuc.chi2();
    cand.pvx = coll.posX();
    cand.pvy = coll.posY();
    cand.pvz = coll.posZ();
    cand.svx = hypNuc.svx();
    cand.svy = hypNuc.svy();
    cand.svz = hypNuc.svz();
    cand.px = hypNuc.px();
    cand.py = hypNuc.py();
    cand.pz = hypNuc.pz();

    auto daughterTracks = hypNuc.hypKfTrack_as<aod::HypKfTracks>();
    for (const auto& track : daughterTracks) {
      TrackProperties daughter;
      daughter.tpcNcls = track.tpcNcluster();
      daughter.itsNcls = track.itsNcluster();
      daughter.tpcChi2 = track.tpcChi2NCl();
      daughter.itsChi2 = track.itsChi2NCl();
      daughter.itsMeanClsSizeL = track.itsMeanClsSize() * track.lambda();
      daughter.rigidity = track.rigidity();
      daughter.tpcSignal = track.tpcSignal();
      daughter.tpcNsigma = track.tpcNsigma();
      daughter.tpcNsigmaNhp = track.tpcNsigmaNhp();
      daughter.tpcNsigmaNlp = track.tpcNsigmaNlp();
      daughter.tofMass = track.tofMass();
      daughter.dcaXY = track.dcaXY();
      daughter.dcaZ = track.dcaZ();
      daughter.isPvContributor = track.isPVContributor();
      cand.daughterTracks.push_back(daughter);
    }
    int trackCount = 0;
    for (const auto& addOn : addOns) {
      cand.daughterTracks.at(trackCount).x = addOn.x();
      cand.daughterTracks.at(trackCount).y = addOn.y();
      cand.daughterTracks.at(trackCount).z = addOn.z();
      cand.daughterTracks.at(trackCount).px = addOn.px();
      cand.daughterTracks.at(trackCount).py = addOn.py();
      cand.daughterTracks.at(trackCount).pz = addOn.pz();
      trackCount++;
    }

    cand.nSingleDaughters = trackCount;
    if (cand.nSingleDaughters < Decays::kThreeBody)
      return;

    trackCount = 0;
    auto subDaughters = hypNuc.hypKfSubD_as<aod::HypKfSubDs>();
    for (const auto& subDaughter : subDaughters) {
      cand.daughterTracks.at(trackCount++).subMass = subDaughter.subMass();
    }
  }
  //___________________________________________________________________________________________________________________________________________________________

  void processData(aod::HypKfHypNucs const& hypNucs, aod::HypKfColls const& hypKfColls, aod::HypKfTracks const& hypKfTrks, aod::HypKfDaughtAdds const& hypKfDAdd, aod::HypKfSubDs const& hypKfDSub)
  {
    for (const auto& hypNuc : hypNucs) {
      if (cfgSpecies && std::abs(hypNuc.species()) != cfgSpecies)
        continue;
      HyperNucleus candidate;
      fillCandidate(candidate, hypNuc, hypNucs, hypKfColls, hypKfTrks, hypKfDAdd, hypKfDSub);
      fillTable(candidate);
    }
  }
  PROCESS_SWITCH(HypKfTreeCreator, processData, "data tree", false);

  //___________________________________________________________________________________________________________________________________________________________

  void processMC(aod::HypKfMcParts const& mcHypNucs, aod::HypKfHypNucs const& hypNucs, aod::HypKfMcColls const&, aod::HypKfColls const& hypKfColls, aod::HypKfTracks const& hypKfTrks, aod::HypKfDaughtAdds const& hypKfDAdd, aod::HypKfSubDs const& hypKfDSub)
  {
    isMC = true;
    for (const auto& mcHypNuc : mcHypNucs) {
      if (cfgSpecies && std::abs(mcHypNuc.species()) != cfgSpecies)
        continue;
      auto mcColl = mcHypNuc.hypKfMcColl();
      const auto mcParticleIdx = mcHypNuc.globalIndex();
      auto hypNucsByMc = hypNucs.sliceBy(perMcParticle, mcParticleIdx);
      HyperNucleus candidate, hypNucDaughter;
      candidate.speciesMC = mcHypNuc.species();
      candidate.pdgCode = mcHypNuc.pdgCode();
      candidate.isMatterMC = mcHypNuc.isMatter();
      candidate.isPhysicalPrimary = mcHypNuc.isPhysicalPrimary();
      candidate.mcPhysicalPrimary = mcHypNuc.isPhysicalPrimary();
      candidate.passedEvSelMC = mcColl.passedEvSel();
      candidate.yGen = mcHypNuc.y();
      candidate.ptGen = mcHypNuc.pt();
      candidate.ctGen = ct(mcColl, mcHypNuc);
      candidate.isReconstructed = 0;
      candidate.cpaPvGen = cpa(mcColl, mcHypNuc);
      candidate.pxGen = mcHypNuc.px();
      candidate.pyGen = mcHypNuc.py();
      candidate.pzGen = mcHypNuc.pz();
      candidate.pvxGen = mcColl.posX();
      candidate.pvyGen = mcColl.posY();
      candidate.pvzGen = mcColl.posZ();
      candidate.svxGen = mcHypNuc.svx();
      candidate.svyGen = mcHypNuc.svy();
      candidate.svzGen = mcHypNuc.svz();
      for (const auto& hypNuc : hypNucsByMc) {
        auto coll = hypNuc.hypKfColl();
        if (coll.hypKfMcCollId() == mcHypNuc.hypKfMcCollId()) {
          candidate.collisionMcTrue = true;
        }
        candidate.isReconstructed++;
        fillCandidate(candidate, hypNucs.rawIteratorAt(hypNuc.globalIndex()), hypNucs, hypKfColls, hypKfTrks, hypKfDAdd, hypKfDSub);
      }
      fillTable(candidate);
    }
  }
  PROCESS_SWITCH(HypKfTreeCreator, processMC, "MC tree", false);

  //___________________________________________________________________________________________________________________________________________________________
  std::vector<float> dcaTracksAll(std::vector<arr3>& posVec, TString opt = "")
  {
    std::vector<float> vec;
    int n = posVec.size();
    for (int i = 0; i < (n - 1); i++) {
      for (int j = (i + 1); j < n; j++) {
        vec.push_back(dcaTracks(posVec, i, j, opt));
      }
    }
    return vec;
  }
  template <class T>
  std::vector<float> dcaTrackSvAll(std::vector<arr3>& posVec, T const& hypNuc, TString opt = "")
  {
    std::vector<float> vec;
    for (size_t i = 0; i < posVec.size(); i++) {
      vec.push_back(dcaTrackSv(posVec, i, hypNuc, opt));
    }
    return vec;
  }

  float maxValue(std::vector<float> vec)
  {
    return *max_element(vec.begin(), vec.end());
  }
  float meanValue(std::vector<float> vec)
  {
    float sum = 0;
    for (const auto& value : vec)
      sum += value;
    return sum / vec.size();
  }
  float mean2Value(std::vector<float> vec)
  {
    float sum = 0;
    for (const auto& value : vec)
      sum += (value * value);
    return std::sqrt(sum / vec.size());
  }

  float dcaTracks(std::vector<arr3> v, int track1, int track2, TString opt = "XY")
  {
    if (opt == "XY")
      return RecoDecay::distanceXY(v.at(track1), v.at(track2));
    else if (opt == "Z")
      return std::abs(v.at(track1).at(2) - v.at(track2).at(2));
    else
      return RecoDecay::distance(v.at(track1), v.at(track2));
  }
  template <class T>
  float dcaTrackSv(std::vector<arr3>& v, int track, T const& hypNuc, TString opt = "")
  {
    if (opt == "XY")
      return RecoDecay::distanceXY(v.at(track), decayVtx(hypNuc));
    else if (opt == "Z")
      return std::abs(v.at(track).at(2) - decayVtx(hypNuc).at(2));
    else
      return RecoDecay::distance(v.at(track), decayVtx(hypNuc));
  }
  template <class T>
  std::vector<arr3> posVector(T const& addons)
  {
    std::vector<arr3> v;
    for (const auto& pos : addons) {
      v.push_back(std::array{pos.x(), pos.y(), pos.z()});
    }
    return v;
  }
  template <class T>
  arr3 primVtx(T const& coll)
  {
    return std::array{coll.posX(), coll.posY(), coll.posZ()};
  }
  template <class T>
  arr3 decayVtx(T const& hypNuc)
  {
    return std::array{hypNuc.svx(), hypNuc.svy(), hypNuc.svz()};
  }
  template <class T>
  arr3 momenta(T const& hypNuc)
  {
    return std::array{hypNuc.px(), hypNuc.py(), hypNuc.pz()};
  }
  template <class TColl, class TPart>
  float decayLength(TColl const& coll, TPart const& hypNuc)
  {
    return RecoDecay::distance(primVtx(coll), decayVtx(hypNuc));
  }
  template <class TColl, class TPart>
  float ct(TColl const& coll, TPart const& hypNuc)
  {
    return RecoDecay::ct(momenta(hypNuc), decayLength(coll, hypNuc), hypNuc.mass());
  }
  template <class TColl, class TPart>
  double cpa(TColl const& coll, TPart const& hypNuc)
  {
    return RecoDecay::cpa(primVtx(coll), decayVtx(hypNuc), momenta(hypNuc));
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HypKfTreeCreator>(cfgc)};
}
