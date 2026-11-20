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
/// \file twoFwdTracksEventTableProducer.cxx
/// \brief Produces derived table with two forward tracks per event from UD tables created by upcCandProducerGlobalMuon
///
/// \author Roman Lavicka <roman.lavicka@cern.ch>, Austrian Academy of Sciences & SMI
/// \since  02.10.2026
//

#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"
#include "PWGUD/DataModel/TwoFwdTracksEventTables.h"
#include "PWGUD/DataModel/UDTables.h"

#include <CommonConstants/PhysicsConstants.h>
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

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TH1.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::constants::physics;

namespace
{
constexpr int NumTracks = 2;             // the table stores exactly two forward tracks per event
constexpr int DefaultInt = -999;         // value of integer columns without information
constexpr float DefaultFloat = -999.f;   // value of float columns without information
constexpr uint32_t DefaultMIDBoards = 0; // value of MID-board columns without information

enum EventSelection {
  EvSelAll = 0,
  EvSelFV0Veto,
  EvSelNumContrib,
  EvSelTwoTracks,
  EvSelInvariantMass,
  EvSelSystemPt,
  EvSelOneTrackMomentum,
  EvSelBothTracksMomentum,
  NEvSels
};

enum FwdTrackSelection {
  TrkSelAll = 0,
  TrkSelType,
  TrkSelPt,
  TrkSelEta,
  TrkSelRabs,
  TrkSelPDca,
  TrkSelChi2,
  TrkSelChi2MatchMCHMFT,
  NTrkSels
};

enum TruthTrouble {
  TroubleTooManyMothers = 1,
  TroubleTooManyDaughters,
  TroubleNoTrack,
  TroubleMoreTracks,
  TroubleDifferentCandidates,
  NTroubles
};

// reco-level information of one row, shared by the measured and the simulated table
struct RecoInfo {
  int32_t runNumber{DefaultInt};
  uint64_t bc{0};
  int numContrib{DefaultInt};
  float posX{DefaultFloat};
  float posY{DefaultFloat};
  float posZ{DefaultFloat};
  float totalFV0AmplitudeA{DefaultFloat};
  std::vector<float> amplitudesV0A;
  std::vector<int8_t> ampRelBCsV0A;
  float energyCommonZNA{DefaultFloat};
  float energyCommonZNC{DefaultFloat};
  float timeZNA{DefaultFloat};
  float timeZNC{DefaultFloat};
  float px[NumTracks] = {DefaultFloat, DefaultFloat};
  float py[NumTracks] = {DefaultFloat, DefaultFloat};
  float pz[NumTracks] = {DefaultFloat, DefaultFloat};
  int sign[NumTracks] = {DefaultInt, DefaultInt};
  float time[NumTracks] = {DefaultFloat, DefaultFloat};
  float timeRes[NumTracks] = {DefaultFloat, DefaultFloat};
  int type[NumTracks] = {DefaultInt, DefaultInt};
  int nClusters[NumTracks] = {DefaultInt, DefaultInt};
  float pDca[NumTracks] = {DefaultFloat, DefaultFloat};
  float rAbs[NumTracks] = {DefaultFloat, DefaultFloat};
  float chi2[NumTracks] = {DefaultFloat, DefaultFloat};
  float chi2MatchMCHMID[NumTracks] = {DefaultFloat, DefaultFloat};
  float chi2MatchMCHMFT[NumTracks] = {DefaultFloat, DefaultFloat};
  int mchBitMap[NumTracks] = {DefaultInt, DefaultInt};
  int midBitMap[NumTracks] = {DefaultInt, DefaultInt};
  uint32_t midBoards[NumTracks] = {DefaultMIDBoards, DefaultMIDBoards};
};

// truth-level information of one row of the simulated table
struct TruthInfo {
  int channel{-1};
  bool hasRecoColl{false};
  float motherPx[NumTracks] = {DefaultFloat, DefaultFloat};
  float motherPy[NumTracks] = {DefaultFloat, DefaultFloat};
  float motherPz[NumTracks] = {DefaultFloat, DefaultFloat};
  float daugPx[NumTracks] = {DefaultFloat, DefaultFloat};
  float daugPy[NumTracks] = {DefaultFloat, DefaultFloat};
  float daugPz[NumTracks] = {DefaultFloat, DefaultFloat};
  int daugPdgCode[NumTracks] = {DefaultInt, DefaultInt};
  bool problem{false};
};
} // namespace

struct TwoFwdTracksEventTableProducer {
  Produces<o2::aod::TwoFwdTracks> twoFwdTracks;
  Produces<o2::aod::TrueTwoFwdTracks> trueTwoFwdTracks;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // declare configurables
  Configurable<bool> verboseInfo{"verboseInfo", false, {"Print general info to terminal; default it false."}};

  struct : ConfigurableGroup {
    Configurable<bool> useFV0Veto{"useFV0Veto", true, {"Reject candidates with FV0A activity"}};
    Configurable<float> cutMaxFV0Amp{"cutMaxFV0Amp", 100.f, {"Maximum FV0A amplitude allowed"}};
    Configurable<int> cutFV0RelBcRange{"cutFV0RelBcRange", 0, {"FV0A veto applied in BCs [-range, +range] around the candidate BC; 0 = candidate BC only"}};
    Configurable<bool> useNumContrib{"useNumContrib", false, {"Require given number of forward tracks in the candidate (before track selection)"}};
    Configurable<int> cutNumContrib{"cutNumContrib", 2, {"Number of forward tracks required in the candidate"}};
  } cutSample;

  struct : ConfigurableGroup {
    Configurable<bool> acceptGlobalMuonTrack{"acceptGlobalMuonTrack", true, {"Accept MCH-MID-MFT tracks"}};
    Configurable<bool> acceptGlobalForwardTrack{"acceptGlobalForwardTrack", true, {"Accept MCH-MFT tracks"}};
    Configurable<bool> acceptMuonStandaloneTrack{"acceptMuonStandaloneTrack", true, {"Accept MCH-MID tracks"}};
    Configurable<bool> acceptMCHStandaloneTrack{"acceptMCHStandaloneTrack", true, {"Accept MCH-only tracks"}};
    Configurable<bool> applyFwdTrackSelection{"applyFwdTrackSelection", true, {"Apply kinematic and quality cuts defined below on forward tracks"}};
    Configurable<float> cutMinPt{"cutMinPt", 0.f, {"Forward track cut"}};
    Configurable<float> cutMaxPt{"cutMaxPt", 1e10f, {"Forward track cut"}};
    Configurable<float> cutMinEta{"cutMinEta", -4.f, {"Forward track cut"}};
    Configurable<float> cutMaxEta{"cutMaxEta", -2.5f, {"Forward track cut"}};
    Configurable<float> cutMinRabs{"cutMinRabs", 17.6f, {"Forward track cut on radius at the absorber end (cm)"}};
    Configurable<float> cutMaxRabs{"cutMaxRabs", 89.5f, {"Forward track cut on radius at the absorber end (cm)"}};
    Configurable<float> cutMaxPDcaLowRabs{"cutMaxPDcaLowRabs", 594.f, {"Forward track pDCA cut for rAbs below 26.5 cm"}};
    Configurable<float> cutMaxPDcaHighRabs{"cutMaxPDcaHighRabs", 324.f, {"Forward track pDCA cut for rAbs above 26.5 cm"}};
    Configurable<float> cutMaxChi2{"cutMaxChi2", 1e10f, {"Forward track cut on track chi2"}};
    Configurable<float> cutMaxChi2MatchMCHMFT{"cutMaxChi2MatchMCHMFT", 1e10f, {"Forward track cut on MCH-MFT matching chi2 (MFT tracks only)"}};
  } cutFwdTrack;

  struct : ConfigurableGroup {
    Configurable<bool> preselUseMinMomentumOnBothTracks{"preselUseMinMomentumOnBothTracks", false, {"Both tracks are required to fill requirement on minimum momentum."}};
    Configurable<float> preselMinTrackMomentum{"preselMinTrackMomentum", 0.1, {"Requirement on minimum momentum of the track."}};
    Configurable<float> preselSystemPtCut{"preselSystemPtCut", 2.0, {"By default, cut on maximum system pT."}};
    Configurable<bool> preselUseOppositeSystemPtCut{"preselUseOppositeSystemPtCut", false, {"Negates the system pT cut (cut on minimum system pT)."}};
    Configurable<float> preselMinInvariantMass{"preselMinInvariantMass", 2.0, {"Requirement on minimum system invariant mass."}};
    Configurable<float> preselMaxInvariantMass{"preselMaxInvariantMass", 5.0, {"Requirement on maximum system invariant mass."}};
  } cutPreselect;

  using CandidatesFwd = soa::Join<aod::UDCollisions, aod::UDCollisionsSelsFwd>;
  using CandidateFwd = CandidatesFwd::iterator;
  using FwdTracks = soa::Join<aod::UDFwdTracks, aod::UDFwdTracksExtra>;
  using FwdMcTracks = soa::Join<aod::UDFwdTracks, aod::UDFwdTracksExtra, aod::UDMcFwdTrackLabels>;

  // init
  void init(InitContext&)
  {
    if (verboseInfo)
      printMediumMessage("INIT METHOD");

    histos.add("Reco/hSelections", "Effect of selections;;Number of events (-)", HistType::kTH1D, {{NEvSels, -0.5, static_cast<double>(NEvSels) - 0.5}});
    auto hSel = histos.get<TH1>(HIST("Reco/hSelections"));
    hSel->GetXaxis()->SetBinLabel(EvSelAll + 1, "All");
    hSel->GetXaxis()->SetBinLabel(EvSelFV0Veto + 1, "FV0A veto");
    hSel->GetXaxis()->SetBinLabel(EvSelNumContrib + 1, "N contrib.");
    hSel->GetXaxis()->SetBinLabel(EvSelTwoTracks + 1, "2 sel. tracks");
    hSel->GetXaxis()->SetBinLabel(EvSelInvariantMass + 1, "Inv. mass");
    hSel->GetXaxis()->SetBinLabel(EvSelSystemPt + 1, "System pT");
    hSel->GetXaxis()->SetBinLabel(EvSelOneTrackMomentum + 1, "One trk p");
    hSel->GetXaxis()->SetBinLabel(EvSelBothTracksMomentum + 1, "Both trks p");

    histos.add("Reco/hFwdTrackSelections", "Effect of forward track selections;;Number of tracks (-)", HistType::kTH1D, {{NTrkSels, -0.5, static_cast<double>(NTrkSels) - 0.5}});
    auto hTrkSel = histos.get<TH1>(HIST("Reco/hFwdTrackSelections"));
    hTrkSel->GetXaxis()->SetBinLabel(TrkSelAll + 1, "All");
    hTrkSel->GetXaxis()->SetBinLabel(TrkSelType + 1, "Track type");
    hTrkSel->GetXaxis()->SetBinLabel(TrkSelPt + 1, "pT");
    hTrkSel->GetXaxis()->SetBinLabel(TrkSelEta + 1, "#eta");
    hTrkSel->GetXaxis()->SetBinLabel(TrkSelRabs + 1, "R_{abs}");
    hTrkSel->GetXaxis()->SetBinLabel(TrkSelPDca + 1, "pDCA");
    hTrkSel->GetXaxis()->SetBinLabel(TrkSelChi2 + 1, "#chi^{2}");
    hTrkSel->GetXaxis()->SetBinLabel(TrkSelChi2MatchMCHMFT + 1, "#chi^{2}_{MCH-MFT}");

    histos.add("Reco/hNanalyzedPerRun", "N analyzed events per run;Run number (-);Number of analyzed events (-)", HistType::kTH1D, {{1, 0., 1.}});
    histos.add("Reco/hNselectedPerRun", "N selected events per run;Run number (-);Number of selected events (-)", HistType::kTH1D, {{1, 0., 1.}});

    histos.add("Truth/hTroubles", "Counter of unwanted issues;;Number of troubles (-)", HistType::kTH1D, {{NTroubles - 1, 0.5, static_cast<double>(NTroubles) - 0.5}});
    auto hTroubles = histos.get<TH1>(HIST("Truth/hTroubles"));
    hTroubles->GetXaxis()->SetBinLabel(TroubleTooManyMothers + 1, "> 2 mothers");
    hTroubles->GetXaxis()->SetBinLabel(TroubleTooManyDaughters + 1, "> 2 charged daughters");
    hTroubles->GetXaxis()->SetBinLabel(TroubleNoTrack + 1, "Daughter without track");
    hTroubles->GetXaxis()->SetBinLabel(TroubleMoreTracks + 1, "Daughter in > 1 cand.");
    hTroubles->GetXaxis()->SetBinLabel(TroubleDifferentCandidates + 1, "Daughters in diff. cand.");
  } // end init

  bool isAcceptedTrackType(uint8_t trackType)
  {
    using o2::aod::fwdtrack::ForwardTrackTypeEnum;
    switch (trackType) {
      case ForwardTrackTypeEnum::GlobalMuonTrack:
        return cutFwdTrack.acceptGlobalMuonTrack;
      case ForwardTrackTypeEnum::GlobalForwardTrack:
        return cutFwdTrack.acceptGlobalForwardTrack;
      case ForwardTrackTypeEnum::MuonStandaloneTrack:
        return cutFwdTrack.acceptMuonStandaloneTrack;
      case ForwardTrackTypeEnum::MCHStandaloneTrack:
        return cutFwdTrack.acceptMCHStandaloneTrack;
      default:
        return false;
    }
  }

  static bool hasMFT(uint8_t trackType)
  {
    using o2::aod::fwdtrack::ForwardTrackTypeEnum;
    return trackType == ForwardTrackTypeEnum::GlobalMuonTrack || trackType == ForwardTrackTypeEnum::GlobalForwardTrack;
  }

  template <typename T>
  static ROOT::Math::PxPyPzMVector muonVector(T const& track)
  {
    return ROOT::Math::PxPyPzMVector(track.px(), track.py(), track.pz(), MassMuon);
  }

  template <typename T>
  bool isSelectedFwdTrack(T const& track, bool fillHistos = true)
  {
    auto countPassed = [&](FwdTrackSelection step) {
      if (fillHistos)
        histos.fill(HIST("Reco/hFwdTrackSelections"), step);
    };
    countPassed(TrkSelAll);

    if (!isAcceptedTrackType(track.trackType()))
      return false;
    countPassed(TrkSelType);

    if (!cutFwdTrack.applyFwdTrackSelection)
      return true;

    const auto vec = muonVector(track);
    if (vec.Pt() < cutFwdTrack.cutMinPt || vec.Pt() > cutFwdTrack.cutMaxPt)
      return false;
    countPassed(TrkSelPt);

    if (vec.Eta() < cutFwdTrack.cutMinEta || vec.Eta() > cutFwdTrack.cutMaxEta)
      return false;
    countPassed(TrkSelEta);

    if (track.rAtAbsorberEnd() < cutFwdTrack.cutMinRabs || track.rAtAbsorberEnd() > cutFwdTrack.cutMaxRabs)
      return false;
    countPassed(TrkSelRabs);

    float maxPDca = track.rAtAbsorberEnd() < upchelpers::AbsorberMid ? cutFwdTrack.cutMaxPDcaLowRabs : cutFwdTrack.cutMaxPDcaHighRabs;
    if (track.pDca() > maxPDca)
      return false;
    countPassed(TrkSelPDca);

    if (track.chi2() > cutFwdTrack.cutMaxChi2)
      return false;
    countPassed(TrkSelChi2);

    if (hasMFT(track.trackType()) && track.chi2MatchMCHMFT() > cutFwdTrack.cutMaxChi2MatchMCHMFT)
      return false;
    countPassed(TrkSelChi2MatchMCHMFT);

    return true;
  }

  // largest FV0A amplitude within [-maxRelBc, +maxRelBc] around the candidate BC
  template <typename C>
  static float getMaxFV0Amplitude(C const& cand, int maxRelBc)
  {
    const auto& amps = cand.amplitudesV0A();
    const auto& relBCs = cand.ampRelBCsV0A();
    float maxAmp = 0.f;
    for (std::size_t i = 0; i < amps.size(); ++i) {
      if (std::abs(relBCs[i]) <= maxRelBc)
        maxAmp = std::max(maxAmp, amps[i]);
    }
    return maxAmp;
  }

  static float finiteOrDefault(float value)
  {
    return std::isfinite(value) ? value : DefaultFloat;
  }

  template <typename C>
  static void fillEventInfo(RecoInfo& info, C const& cand)
  {
    info.runNumber = cand.runNumber();
    info.bc = cand.globalBC();
    info.numContrib = cand.numContrib();
    info.posX = cand.posX();
    info.posY = cand.posY();
    info.posZ = cand.posZ();
    info.totalFV0AmplitudeA = getMaxFV0Amplitude(cand, 0);
    const auto& amps = cand.amplitudesV0A();
    const auto& relBCs = cand.ampRelBCsV0A();
    info.amplitudesV0A.assign(amps.begin(), amps.end());
    info.ampRelBCsV0A.assign(relBCs.begin(), relBCs.end());
  }

  // the ZDC table is filled by the producer only for candidates with ZDC activity, i.e. it holds at most one row per candidate
  static void fillZdcInfo(RecoInfo& info, aod::UDZdcsReduced const& zdcs)
  {
    for (const auto& zdc : zdcs) {
      info.energyCommonZNA = finiteOrDefault(zdc.energyCommonZNA());
      info.energyCommonZNC = finiteOrDefault(zdc.energyCommonZNC());
      info.timeZNA = finiteOrDefault(zdc.timeZNA());
      info.timeZNC = finiteOrDefault(zdc.timeZNC());
    }
  }

  template <typename T>
  static void fillTrackInfo(RecoInfo& info, int iTrk, T const& track)
  {
    info.px[iTrk] = track.px();
    info.py[iTrk] = track.py();
    info.pz[iTrk] = track.pz();
    info.sign[iTrk] = track.sign();
    info.time[iTrk] = track.trackTime();
    info.timeRes[iTrk] = track.trackTimeRes();
    info.type[iTrk] = track.trackType();
    info.nClusters[iTrk] = track.nClusters();
    info.pDca[iTrk] = track.pDca();
    info.rAbs[iTrk] = track.rAtAbsorberEnd();
    info.chi2[iTrk] = track.chi2();
    info.chi2MatchMCHMID[iTrk] = track.chi2MatchMCHMID();
    info.chi2MatchMCHMFT[iTrk] = track.chi2MatchMCHMFT();
    info.mchBitMap[iTrk] = track.mchBitMap();
    info.midBitMap[iTrk] = track.midBitMap();
    info.midBoards[iTrk] = track.midBoards();
  }

  // writes reco-level columns, followed by optional truth-level columns of the simulated table
  template <typename TCursor, typename... TTruth>
  static void writeRow(TCursor& cursor, RecoInfo& info, TTruth&&... truth)
  {
    cursor(info.runNumber, info.bc, info.numContrib, info.posX, info.posY, info.posZ,
           info.totalFV0AmplitudeA, info.amplitudesV0A, info.ampRelBCsV0A,
           info.energyCommonZNA, info.energyCommonZNC, info.timeZNA, info.timeZNC,
           info.px, info.py, info.pz, info.sign, info.time, info.timeRes,
           info.type, info.nClusters, info.pDca, info.rAbs, info.chi2, info.chi2MatchMCHMID, info.chi2MatchMCHMFT,
           info.mchBitMap, info.midBitMap, info.midBoards[0], info.midBoards[1],
           std::forward<TTruth>(truth)...);
  }

  static int getTrueChannel(const int pdgCodes[NumTracks])
  {
    const int part1 = enumMyParticle(pdgCodes[0]);
    const int part2 = enumMyParticle(pdgCodes[1]);
    if (part1 == P_ELECTRON && part2 == P_ELECTRON)
      return CH_EE;
    if (part1 == P_MUON && part2 == P_MUON)
      return CH_MUMU;
    if ((part1 == P_ELECTRON && (part2 == P_PION || part2 == P_MUON)) || (part2 == P_ELECTRON && (part1 == P_PION || part1 == P_MUON)))
      return CH_EMUPI;
    return -1;
  }

  void processDataFwd(CandidateFwd const& cand,
                      FwdTracks const& tracks,
                      aod::UDZdcsReduced const& zdcs)
  {
    const char* srun = Form("%d", cand.runNumber());
    histos.get<TH1>(HIST("Reco/hNanalyzedPerRun"))->Fill(srun, 1);

    histos.fill(HIST("Reco/hSelections"), EvSelAll);

    if (cutSample.useFV0Veto && getMaxFV0Amplitude(cand, cutSample.cutFV0RelBcRange) > cutSample.cutMaxFV0Amp)
      return;
    histos.fill(HIST("Reco/hSelections"), EvSelFV0Veto);

    if (cutSample.useNumContrib && cand.numContrib() != cutSample.cutNumContrib)
      return;
    histos.fill(HIST("Reco/hSelections"), EvSelNumContrib);

    RecoInfo info;
    std::array<ROOT::Math::PxPyPzMVector, NumTracks> daug;
    int nSelected = 0;
    for (const auto& track : tracks) {
      if (!isSelectedFwdTrack(track))
        continue;
      if (nSelected < NumTracks) {
        fillTrackInfo(info, nSelected, track);
        daug[nSelected] = muonVector(track);
      }
      nSelected++;
    }

    // Critical selection, without it the rest of the process function will fail
    if (nSelected != NumTracks)
      return;
    histos.fill(HIST("Reco/hSelections"), EvSelTwoTracks);

    const auto mother = daug[0] + daug[1];

    // Apply system selections
    // invariant mass
    if (mother.M() < cutPreselect.preselMinInvariantMass || mother.M() > cutPreselect.preselMaxInvariantMass)
      return;
    histos.fill(HIST("Reco/hSelections"), EvSelInvariantMass);

    // system pt
    if (cutPreselect.preselUseOppositeSystemPtCut ? mother.Pt() < cutPreselect.preselSystemPtCut : mother.Pt() > cutPreselect.preselSystemPtCut)
      return;
    histos.fill(HIST("Reco/hSelections"), EvSelSystemPt);

    // one track momentum
    if (daug[0].P() < cutPreselect.preselMinTrackMomentum && daug[1].P() < cutPreselect.preselMinTrackMomentum)
      return;
    histos.fill(HIST("Reco/hSelections"), EvSelOneTrackMomentum);

    // both tracks momentum
    if (cutPreselect.preselUseMinMomentumOnBothTracks && (daug[0].P() < cutPreselect.preselMinTrackMomentum || daug[1].P() < cutPreselect.preselMinTrackMomentum))
      return;
    histos.fill(HIST("Reco/hSelections"), EvSelBothTracksMomentum);

    histos.get<TH1>(HIST("Reco/hNselectedPerRun"))->Fill(srun, 1);

    fillEventInfo(info, cand);
    fillZdcInfo(info, zdcs);
    writeRow(twoFwdTracks, info);
  }
  PROCESS_SWITCH(TwoFwdTracksEventTableProducer, processDataFwd, "Iterate UD tables with measured data created by upcCandProducerGlobalMuon.", true);

  PresliceUnsorted<aod::UDMcParticles> partPerMcCollision = aod::udmcparticle::udMcCollisionId;
  PresliceUnsorted<FwdMcTracks> trackPerMcParticle = aod::udmcfwdtracklabel::udMcParticleId;

  void processMonteCarlo(aod::UDMcCollisions const& mcCollisions,
                         aod::UDMcParticles const& parts,
                         CandidatesFwd const& candidates,
                         FwdMcTracks const& tracks)
  {
    // start loop over generated collisions
    for (const auto& mcColl : mcCollisions) {
      RecoInfo info;
      TruthInfo truth;
      std::vector<int64_t> daughterIds;

      // store a charged particle as truth daughter; returns false if there are already too many of them
      auto addTrueDaughter = [&](auto const& daughter) {
        // check if it is the charged particle (= no pi0, photon or neutrino)
        if (enumMyParticle(daughter.pdgCode()) == -1)
          return true;
        // check we do not have more than 2 charged daughters in total
        if (static_cast<int>(daughterIds.size()) >= NumTracks) {
          if (verboseInfo)
            printLargeMessage("Truth collision has more than 2 total charged daughters. Breaking the daughter loop.");
          histos.fill(HIST("Truth/hTroubles"), TroubleTooManyDaughters);
          truth.problem = true;
          return false;
        }
        const auto iDaug = daughterIds.size();
        truth.daugPx[iDaug] = daughter.px();
        truth.daugPy[iDaug] = daughter.py();
        truth.daugPz[iDaug] = daughter.pz();
        truth.daugPdgCode[iDaug] = daughter.pdgCode();
        daughterIds.push_back(daughter.globalIndex());
        return true;
      };

      // get particles associated to generated collision
      const auto& partsFromMcColl = parts.sliceBy(partPerMcCollision, mcColl.globalIndex());
      int countMothers = 0;
      for (const auto& particle : partsFromMcColl) {
        // select only mothers with checking if particle has no mother
        if (particle.has_mothers())
          continue;
        countMothers++;
        // check the generated collision does not have more than 2 mothers
        if (countMothers > NumTracks) {
          if (verboseInfo)
            printLargeMessage("Truth collision has more than 2 no mother particles. Breaking the particle loop.");
          histos.fill(HIST("Truth/hTroubles"), TroubleTooManyMothers);
          truth.problem = true;
          break;
        }
        // fill info for each mother
        truth.motherPx[countMothers - 1] = particle.px();
        truth.motherPy[countMothers - 1] = particle.py();
        truth.motherPz[countMothers - 1] = particle.pz();

        if (particle.has_daughters()) {
          for (const auto& daughter : particle.daughters_as<aod::UDMcParticles>()) {
            if (!addTrueDaughter(daughter))
              break;
          }
        } else {
          // motherless final-state particle, e.g. muon from gamma gamma -> mu mu, is its own daughter
          addTrueDaughter(particle);
        }
        if (truth.problem)
          break;
      } // particles

      // get the reconstructed tracks of each daughter (how well the daughter was reconstructed).
      // One muon can be stored in several candidates (as MCH-MID-MFT track in an anchor candidate and as MCH-MID track
      // in an MCH-MID candidate), so collect per daughter the selected track in each candidate, preferring the lowest track type.
      std::array<std::map<int32_t, int64_t>, NumTracks> trackIdPerCand; // candidate id -> track id, per daughter
      int nRecoDaughters = 0;
      for (std::size_t iDaug = 0; iDaug < daughterIds.size(); ++iDaug) {
        auto& candTracks = trackIdPerCand[iDaug];
        for (const auto& trk : tracks.sliceBy(trackPerMcParticle, daughterIds[iDaug])) {
          if (!isSelectedFwdTrack(trk, false))
            continue;
          auto it = candTracks.find(trk.udCollisionId());
          if (it == candTracks.end() || trk.trackType() < tracks.iteratorAt(it->second).trackType())
            candTracks[trk.udCollisionId()] = trk.globalIndex();
        }
        if (candTracks.empty()) {
          if (verboseInfo)
            printLargeMessage("Daughter has no associated track. Skipping this daughter.");
          histos.fill(HIST("Truth/hTroubles"), TroubleNoTrack);
          truth.problem = true;
          continue;
        }
        if (candTracks.size() > 1)
          histos.fill(HIST("Truth/hTroubles"), TroubleMoreTracks);
        nRecoDaughters++;
      } // daughters

      // choose the candidate holding most daughters; on a tie the one with the lowest summed track type (i.e. with MFT)
      std::map<int32_t, std::pair<int, int>> candScore; // candidate id -> (number of daughters, summed track type)
      for (const auto& candTracks : trackIdPerCand) {
        for (const auto& [cand, trkId] : candTracks) {
          candScore[cand].first++;
          candScore[cand].second += tracks.iteratorAt(trkId).trackType();
        }
      }
      int32_t candId = -1;
      std::pair<int, int> bestScore{0, 0};
      for (const auto& [cand, score] : candScore) {
        if (candId < 0 || score.first > bestScore.first || (score.first == bestScore.first && score.second < bestScore.second)) {
          candId = cand;
          bestScore = score;
        }
      }

      // fill info for reconstructed collision and tracks of the chosen candidate
      if (candId >= 0) {
        if (bestScore.first < nRecoDaughters) {
          if (verboseInfo)
            printLargeMessage("Daughters are reconstructed in different candidates.");
          histos.fill(HIST("Truth/hTroubles"), TroubleDifferentCandidates);
          truth.problem = true;
        }
        truth.hasRecoColl = true;
        fillEventInfo(info, candidates.iteratorAt(candId));
        for (int iDaug = 0; iDaug < NumTracks; ++iDaug) {
          auto it = trackIdPerCand[iDaug].find(candId);
          if (it != trackIdPerCand[iDaug].end())
            fillTrackInfo(info, iDaug, tracks.iteratorAt(it->second));
        }
      }

      truth.channel = getTrueChannel(truth.daugPdgCode);

      writeRow(trueTwoFwdTracks, info, // no ZDC info in MC
               truth.channel, truth.hasRecoColl, mcColl.posX(), mcColl.posY(), mcColl.posZ(),
               truth.motherPx, truth.motherPy, truth.motherPz, truth.daugPx, truth.daugPy, truth.daugPz, truth.daugPdgCode, truth.problem);
    } // mccollisions
  }
  PROCESS_SWITCH(TwoFwdTracksEventTableProducer, processMonteCarlo, "Iterate UD tables with simulated data created by upcCandProducerGlobalMuon.", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<TwoFwdTracksEventTableProducer>(cfgc)};
}
