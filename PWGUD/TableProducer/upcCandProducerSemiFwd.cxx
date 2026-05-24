// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file   upcCandProducerSemiFwd.cxx
/// \brief  UPC candidate producer combining forward MCH-MID muons with central-barrel tracks.
///         Anchor: MCH-MID forward track (precise MID timing). Companions: barrel tracks
///         whose BC falls within fBcWindowBarrel of the anchor BC. One UD candidate per anchor BC.
/// \author Roman Lavička, lavicka.roman@gmail.com
/// \since  10.05.2026

#include "PWGUD/Core/UPCCutparHolder.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/LHCConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
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
#include <GlobalTracking/MatchGlobalFwd.h>
#include <MCHTracking/TrackExtrap.h>
#include <ReconstructionDataFormats/GlobalFwdTrack.h>

#include <Math/MatrixRepresentationsStatic.h>
#include <Math/SMatrix.h>
#include <TH1.h>
#include <TMath.h>

#include <cstdint>
#include <cstdlib>
#include <map>
#include <numeric>
#include <string>
#include <vector>

using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcCandProducerSemiFwd {
  bool fDoMC{false};

  std::map<int32_t, int32_t> fNewPartIDs;

  Produces<o2::aod::UDMcCollisions> udMCCollisions;
  Produces<o2::aod::UDMcParticles> udMCParticles;
  Produces<o2::aod::UDMcFwdTrackLabels> udFwdTrackLabels;
  Produces<o2::aod::UDMcTrackLabels> udTrackLabels;

  Produces<o2::aod::UDFwdTracks> udFwdTracks;
  Produces<o2::aod::UDFwdTracksExtra> udFwdTracksExtra;
  Produces<o2::aod::UDFwdIndices> udFwdIndices;

  Produces<o2::aod::UDTracks> udTracks;
  Produces<o2::aod::UDTracksExtra> udTracksExtra;
  Produces<o2::aod::UDTracksDCA> udTracksDCA;
  Produces<o2::aod::UDTracksFlags> udTracksFlags;
  Produces<o2::aod::UDTracksPID> udTracksPID;

  Produces<o2::aod::UDCollisions> eventCandidates;
  Produces<o2::aod::UDCollisionsSelsFwd> eventCandidatesSelsFwd;
  Produces<o2::aod::UDZdcsReduced> udZdcsReduced;

  Configurable<int> fSignalGenID{"fSignalGenID", 1, "Signal generator ID"};

  UPCCutparHolder fUpcCuts = UPCCutparHolder();
  MutableConfigurable<UPCCutparHolder> fUpcCutsConf{"fUpcCutsConf", {}, "UPC event cuts"};
  Configurable<uint64_t> fBcWindowFITAmps{"fBcWindowFITAmps", 20, "BC range for T0A/V0A amplitudes array [-range, +(range-1)]"};
  Configurable<int> fBcWindowBarrel{"fBcWindowBarrel", 20, "Time window for barrel-track to MCH-MID anchor matching"};
  Configurable<float> fMaxFV0Amp{"fMaxFV0Amp", 100.f, "Max FV0 amplitude in the same BC"};

  using ForwardTracks = o2::soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov>;

  using BarrelTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TracksDCA,
                                     o2::aod::pidTPCFullEl, o2::aod::pidTPCFullMu, o2::aod::pidTPCFullPi, o2::aod::pidTPCFullKa, o2::aod::pidTPCFullPr,
                                     o2::aod::TOFSignal, o2::aod::pidTOFbeta,
                                     o2::aod::pidTOFFullEl, o2::aod::pidTOFFullMu, o2::aod::pidTOFFullPi, o2::aod::pidTOFFullKa, o2::aod::pidTOFFullPr>;

  HistogramRegistry histRegistry{"HistRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  int fRun{0};
  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::ccdb::CcdbApi fCCDBApi;
  o2::globaltracking::MatchGlobalFwd fMatching;

  void init(InitContext&)
  {
    fUpcCuts = (UPCCutparHolder)fUpcCutsConf;
    if (fUpcCuts.getUseFwdCuts()) {
      fCCDB->setURL("http://alice-ccdb.cern.ch");
      fCCDB->setCaching(true);
      fCCDB->setLocalObjectValidityChecking();
      fCCDBApi.init("http://alice-ccdb.cern.ch");
    }
    const AxisSpec axisSelFwd{upchelpers::kNFwdSels, 0., static_cast<double>(upchelpers::kNFwdSels), ""};
    histRegistry.add("MuonsSelCounter", "", kTH1F, {axisSelFwd});
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kFwdSelAll + 1, "All");
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kFwdSelPt + 1, "Pt");
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kFwdSelEta + 1, "Eta");
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kFwdSelRabs + 1, "Rabs");
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kFwdSelpDCA + 1, "pDCA");
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kFwdSelChi2 + 1, "Chi2");

    const AxisSpec axisSelBar{upchelpers::kNBarrelSels, 0., static_cast<double>(upchelpers::kNBarrelSels), ""};
    histRegistry.add("BarrelsSelCounter", "", kTH1F, {axisSelBar});
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelAll + 1, "All");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelHasTOF + 1, "HasTOF");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelPt + 1, "Pt");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelEta + 1, "Eta");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelITSNCls + 1, "ITSNCls");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelITSChi2 + 1, "ITSChi2");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelTPCNCls + 1, "TPCNCls");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelTPCChi2 + 1, "TPCChi2");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelDCAXY + 1, "DCAXY");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(upchelpers::kBarrelSelDCAZ + 1, "DCAZ");
  }

  bool fwdCut(const o2::dataformats::GlobalFwdTrack& pft, const ForwardTracks::iterator& fwdTrack)
  {
    histRegistry.fill(HIST("MuonsSelCounter"), upchelpers::kFwdSelAll, 1);
    auto pt = pft.getPt();
    auto eta = pft.getEta();
    auto pdca = fwdTrack.pDca();
    auto rabs = fwdTrack.rAtAbsorberEnd();
    auto chi2 = fwdTrack.chi2();
    bool passPt = pt > fUpcCuts.getFwdPtLow() && pt < fUpcCuts.getFwdPtHigh();
    bool passEta = eta > fUpcCuts.getFwdEtaLow() && eta < fUpcCuts.getFwdEtaHigh();
    bool passRabs = rabs > fUpcCuts.getMuonRAtAbsorberEndLow() && rabs < fUpcCuts.getMuonRAtAbsorberEndHigh();
    bool passPDca = rabs < upchelpers::AbsorberMid ? pdca < fUpcCuts.getMuonPDcaHighFirst() : pdca < fUpcCuts.getMuonPDcaHighSecond();
    bool passChi2 = chi2 > fUpcCuts.getFwdChi2Low() && chi2 < fUpcCuts.getFwdChi2High();
    if (passPt)
      histRegistry.fill(HIST("MuonsSelCounter"), upchelpers::kFwdSelPt, 1);
    if (passEta)
      histRegistry.fill(HIST("MuonsSelCounter"), upchelpers::kFwdSelEta, 1);
    if (passRabs)
      histRegistry.fill(HIST("MuonsSelCounter"), upchelpers::kFwdSelRabs, 1);
    if (passPDca)
      histRegistry.fill(HIST("MuonsSelCounter"), upchelpers::kFwdSelpDCA, 1);
    if (passChi2)
      histRegistry.fill(HIST("MuonsSelCounter"), upchelpers::kFwdSelChi2, 1);
    return passPt && passEta && passRabs && passPDca && passChi2;
  }

  bool barCut(const BarrelTracks::iterator& track)
  {
    histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelAll, 1);
    if (!fUpcCuts.getUseBarCuts())
      return true;
    if (fUpcCuts.getAmbigSwitch() == 1 && !track.isPVContributor())
      return false;
    if (fUpcCuts.getRequireTOF() && !track.hasTOF())
      return false;
    bool passHasTOF = !fUpcCuts.getRequireTOF() || track.hasTOF();
    bool passPt = track.pt() > fUpcCuts.getBarPtLow() && track.pt() < fUpcCuts.getBarPtHigh();
    bool passEta = track.eta() > fUpcCuts.getBarEtaLow() && track.eta() < fUpcCuts.getBarEtaHigh();
    bool passITSNCls = track.itsNCls() >= static_cast<uint8_t>(fUpcCuts.getITSNClusLow()) &&
                       track.itsNCls() <= static_cast<uint8_t>(fUpcCuts.getITSNClusHigh());
    bool passITSChi2 = track.itsChi2NCl() > fUpcCuts.getITSChi2Low() && track.itsChi2NCl() < fUpcCuts.getITSChi2High();
    bool passTPCNCls = track.tpcNClsFound() >= static_cast<int16_t>(fUpcCuts.getTPCNClsLow()) &&
                       track.tpcNClsFound() <= static_cast<int16_t>(fUpcCuts.getTPCNClsHigh());
    bool passTPCChi2 = track.tpcChi2NCl() > fUpcCuts.getTPCChi2Low() && track.tpcChi2NCl() < fUpcCuts.getTPCChi2High();
    bool passDcaZ = track.dcaZ() > fUpcCuts.getDcaZLow() && track.dcaZ() < fUpcCuts.getDcaZHigh();
    bool passDcaXY = true;
    if (fUpcCuts.getCheckMaxDcaXY()) {
      constexpr float DcaXYConst = 0.0105f;
      constexpr float DcaXYSlope = 0.0350f;
      constexpr float DcaXYExp = 1.1f;
      float maxDCA = DcaXYConst + DcaXYSlope / std::pow(track.pt(), DcaXYExp);
      passDcaXY = track.dcaXY() < maxDCA;
    }
    if (passHasTOF)
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelHasTOF, 1);
    if (passPt)
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelPt, 1);
    if (passEta)
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelEta, 1);
    if (passITSNCls)
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelITSNCls, 1);
    if (passITSChi2)
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelITSChi2, 1);
    if (passTPCNCls)
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelTPCNCls, 1);
    if (passTPCChi2)
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelTPCChi2, 1);
    if (passDcaXY)
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelDCAXY, 1);
    if (passDcaZ)
      histRegistry.fill(HIST("BarrelsSelCounter"), upchelpers::kBarrelSelDCAZ, 1);
    return passPt && passEta && passITSNCls && passITSChi2 && passTPCNCls && passTPCChi2 && passDcaXY && passDcaZ;
  }

  void skimMCInfo(o2::aod::McCollisions const& mcCollisions, o2::aod::McParticles const& mcParticles)
  {
    std::vector<int32_t> newEventIDs(mcCollisions.size(), -1);

    int32_t newPartID = 0;
    int32_t newEventID = 0;
    int32_t nMCParticles = mcParticles.size();
    for (int32_t mcPartID = 0; mcPartID < nMCParticles; mcPartID++) {
      const auto& mcPart = mcParticles.iteratorAt(mcPartID);
      if (!mcPart.has_mcCollision())
        continue;
      int32_t mcEventID = mcPart.mcCollisionId();
      const auto& mcEvent = mcCollisions.iteratorAt(mcEventID);
      bool isSignal = mcEvent.getGeneratorId() == fSignalGenID;
      if (!isSignal) {
        continue;
      }
      fNewPartIDs[mcPartID] = newPartID;
      newPartID++;
      if (newEventIDs[mcEventID] == -1) {
        newEventIDs[mcEventID] = newEventID;
        newEventID++;
      }
    }

    std::vector<int32_t> newMotherIDs{};

    for (const auto& item : fNewPartIDs) {
      int32_t mcPartID = item.first;
      const auto& mcPart = mcParticles.iteratorAt(mcPartID);
      int32_t mcEventID = mcPart.mcCollisionId();
      int32_t mappedEventID = newEventIDs[mcEventID];
      if (mcPart.has_mothers()) {
        const auto& motherIDs = mcPart.mothersIds();
        for (const auto& motherID : motherIDs) {
          if (motherID >= nMCParticles) {
            continue;
          }
          auto it = fNewPartIDs.find(motherID);
          if (it != fNewPartIDs.end()) {
            newMotherIDs.push_back(it->second);
          }
        }
      }
      int32_t newDaughterIDs[2] = {-1, -1};
      if (mcPart.has_daughters()) {
        const auto& daughterIDs = mcPart.daughtersIds();
        int32_t firstDaughter = daughterIDs.front();
        int32_t lastDaughter = daughterIDs.back();
        if (firstDaughter >= nMCParticles || lastDaughter >= nMCParticles) {
          continue;
        }
        auto itFirst = fNewPartIDs.find(firstDaughter);
        auto itLast = fNewPartIDs.find(lastDaughter);
        if (itFirst != fNewPartIDs.end() && itLast != fNewPartIDs.end()) {
          newDaughterIDs[0] = fNewPartIDs.at(daughterIDs.front());
          newDaughterIDs[1] = fNewPartIDs.at(daughterIDs.back());
        }
      }
      udMCParticles(mappedEventID, mcPart.pdgCode(), mcPart.getHepMCStatusCode(), mcPart.flags(), newMotherIDs, newDaughterIDs,
                    mcPart.weight(), mcPart.px(), mcPart.py(), mcPart.pz(), mcPart.e());
      newMotherIDs.clear();
    }

    for (int32_t i = 0; i < mcCollisions.size(); i++) {
      if (newEventIDs[i] == -1) {
        continue;
      }
      const auto& mcEvent = mcCollisions.iteratorAt(i);
      const auto& bc = mcEvent.bc();
      udMCCollisions(bc.globalBC(), mcEvent.generatorsID(), mcEvent.posX(), mcEvent.posY(), mcEvent.posZ(),
                     mcEvent.t(), mcEvent.weight(), mcEvent.impactParameter());
    }

    newEventIDs.clear();
  }

  int64_t bcDiff(uint64_t bc1, uint64_t bc2)
  {
    return bc1 > bc2 ? bc1 - bc2 : bc2 - bc1;
  }

  template <typename T>
  T::iterator getStartForScroll(uint64_t inGbc, T& gbcMap)
  {
    auto it1 = gbcMap.lower_bound(inGbc);
    typename T::iterator it;
    if (it1 != gbcMap.end()) {
      auto it2 = it1;
      uint64_t bc1 = it1->first;
      if (it2 != gbcMap.begin())
        --it2;
      uint64_t bc2 = it2->first;
      uint64_t dbc1 = bcDiff(bc1, inGbc);
      uint64_t dbc2 = bcDiff(bc2, inGbc);
      it = (dbc1 <= dbc2) ? it1 : it2;
    } else {
      it = it1;
      --it;
    }
    return it;
  }

  template <typename T, typename F>
  void scrollBackForth(uint64_t inGbc, uint64_t maxDbc, T& gbcMap, F&& func)
  {
    auto it = getStartForScroll(inGbc, gbcMap);
    uint64_t gbc = it->first;
    uint64_t dbc = bcDiff(inGbc, gbc);

    int count = 0;
    while (dbc <= maxDbc) {
      func(it, gbc);
      count++;
      if (it == gbcMap.begin())
        break;
      --it;
      gbc = it->first;
      dbc = bcDiff(inGbc, gbc);
    }

    std::advance(it, count + 1);

    if (it == gbcMap.end())
      return;

    gbc = it->first;
    dbc = bcDiff(inGbc, gbc);

    while (dbc <= maxDbc) {
      func(it, gbc);
      ++it;
      if (it == gbcMap.end())
        break;
      gbc = it->first;
      dbc = bcDiff(inGbc, gbc);
    }
  }

  void getBarrelTrackIds(uint64_t inGbc, std::map<uint64_t, std::vector<int64_t>>& tracksPerBC,
                         uint64_t maxDbc, std::map<int64_t, uint64_t>& outTrkIds)
  {
    auto fillIds = [&outTrkIds](std::map<uint64_t, std::vector<int64_t>>::iterator& inIt, uint64_t gbc) {
      std::vector<int64_t>& ids = inIt->second;
      for (const auto& id : ids)
        outTrkIds[id] = gbc;
    };
    scrollBackForth(inGbc, maxDbc, tracksPerBC, fillIds);
  }

  void getFV0Amplitudes(uint64_t inGbc, o2::aod::FV0As const& fv0s, uint64_t maxDbc,
                        std::map<uint64_t, int64_t>& mapBcs, std::vector<float>& amps, std::vector<int8_t>& relBcs)
  {
    auto fillAmps = [this, &fv0s, &amps, &relBcs, inGbc](std::map<uint64_t, int64_t>::iterator& inIt, uint64_t gbc) {
      int64_t fv0Id = inIt->second;
      const auto& fv0 = fv0s.iteratorAt(fv0Id);
      const auto& amplitudes = fv0.amplitude();
      float totalAmp = std::accumulate(amplitudes.begin(), amplitudes.end(), 0.f);
      if (totalAmp > 0.f) {
        amps.push_back(totalAmp);
        auto relBc = static_cast<int8_t>(bcDiff(gbc, inGbc));
        if (gbc < inGbc)
          relBc *= -1;
        relBcs.push_back(relBc);
      }
    };
    scrollBackForth(inGbc, maxDbc, mapBcs, fillAmps);
  }

  auto propagateToZero(ForwardTracks::iterator const& muon)
  {
    using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
    using SMatrix5 = ROOT::Math::SVector<double, 5>;
    SMatrix5 tpars(muon.x(), muon.y(), muon.phi(), muon.tgl(), muon.signed1Pt());
    std::vector<double> v1{muon.cXX(), muon.cXY(), muon.cYY(), muon.cPhiX(), muon.cPhiY(),
                           muon.cPhiPhi(), muon.cTglX(), muon.cTglY(), muon.cTglPhi(), muon.cTglTgl(),
                           muon.c1PtX(), muon.c1PtY(), muon.c1PtPhi(), muon.c1PtTgl(), muon.c1Pt21Pt2()};
    SMatrix55 tcovs(v1.begin(), v1.end());
    o2::dataformats::GlobalFwdTrack propmuon;
    o2::dataformats::GlobalFwdTrack track;
    track.setParameters(tpars);
    track.setZ(muon.z());
    track.setCovariances(tcovs);
    auto mchTrack = fMatching.FwdtoMCH(track);
    o2::mch::TrackExtrap::extrapToVertex(mchTrack, 0., 0., 0., 0., 0.);
    auto proptrack = fMatching.MCHtoFwd(mchTrack);
    propmuon.setParameters(proptrack.getParameters());
    propmuon.setZ(proptrack.getZ());
    propmuon.setCovariances(proptrack.getCovariances());
    return propmuon;
  }

  bool addToFwdTable(int64_t candId, int64_t trackId, uint64_t gbc, float trackTime, ForwardTracks const& fwdTracks,
                     const o2::aod::McFwdTrackLabels* mcFwdTrackLabels)
  {
    constexpr float NoMchMftMatch = -1.f;
    const auto& track = fwdTracks.iteratorAt(trackId);
    float px;
    float py;
    float pz;
    int sign;
    if (fUpcCuts.getUseFwdCuts()) {
      auto pft = propagateToZero(track);
      bool pass = fwdCut(pft, track);
      if (!pass)
        return false;
      px = pft.getPx();
      py = pft.getPy();
      pz = pft.getPz();
      sign = (pft.getInvQPt() > 0) ? 1 : -1;
    } else {
      px = track.px();
      py = track.py();
      pz = track.pz();
      sign = track.sign();
    }
    udFwdTracks(candId, px, py, pz, sign, gbc, trackTime, track.trackTimeRes());
    udFwdTracksExtra(track.trackType(), track.nClusters(), track.pDca(), track.rAtAbsorberEnd(), track.chi2(),
                     track.chi2MatchMCHMID(), NoMchMftMatch, track.mchBitMap(), track.midBitMap(), track.midBoards());
    udFwdIndices(candId, trackId, track.matchMCHTrackId(), -1);
    if (fDoMC) {
      const auto& label = mcFwdTrackLabels->iteratorAt(trackId);
      uint16_t mcMask = label.mcMask();
      auto it = fNewPartIDs.find(label.mcParticleId());
      int32_t newPartID = it != fNewPartIDs.end() ? it->second : -1;
      udFwdTrackLabels(newPartID, mcMask);
    }
    return true;
  }

  bool addToBarrelTable(int64_t candId, int64_t trackId, uint64_t gbc, float trackTime, BarrelTracks const& barrelTracks,
                        const o2::aod::McTrackLabels* mcTrackLabels)
  {
    const auto& track = barrelTracks.iteratorAt(trackId);
    if (!barCut(track))
      return false;
    int32_t colId = track.collisionId() >= 0 ? track.collisionId() : -1;
    udTracks(candId, track.px(), track.py(), track.pz(), track.sign(), gbc, trackTime, track.trackTimeRes());
    udTracksExtra(track.tpcInnerParam(), track.itsClusterSizes(), track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(),
                  track.tpcNClsFindableMinusCrossedRows(), track.tpcNClsShared(), track.trdPattern(), track.itsChi2NCl(),
                  track.tpcChi2NCl(), track.trdChi2(), track.tofChi2(), track.tpcSignal(), track.tofSignal(), track.trdSignal(),
                  track.length(), track.tofExpMom(), track.detectorMap());
    udTracksDCA(track.dcaZ(), track.dcaXY());
    udTracksFlags(colId, track.isPVContributor());
    udTracksPID(track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                track.beta(), track.betaerror(),
                track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr());
    if (fDoMC) {
      const auto& label = mcTrackLabels->iteratorAt(trackId);
      uint16_t mcMask = label.mcMask();
      auto it = fNewPartIDs.find(label.mcParticleId());
      int32_t newPartID = it != fNewPartIDs.end() ? it->second : -1;
      udTrackLabels(newPartID, mcMask);
    }
    return true;
  }

  void createCandidates(ForwardTracks const& fwdTracks,
                        o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                        BarrelTracks const& barrelTracks,
                        o2::aod::AmbiguousTracks const& ambBarrelTracks,
                        o2::aod::BCs const& bcs,
                        o2::aod::Collisions const& collisions,
                        o2::aod::FV0As const& fv0s,
                        o2::aod::Zdcs const& zdcs,
                        const o2::aod::McFwdTrackLabels* mcFwdTrackLabels,
                        const o2::aod::McTrackLabels* mcBarrelTrackLabels)
  {
    using o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack;

    int32_t runNumber = bcs.iteratorAt(0).runNumber();
    if (fUpcCuts.getUseFwdCuts()) {
      if (runNumber != fRun) {
        fRun = runNumber;
        std::map<std::string, std::string> metadata;
        auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(fCCDBApi, fRun);
        auto ts = soreor.second;
        auto grpmag = fCCDBApi.retrieveFromTFileAny<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", metadata, ts);
        o2::base::Propagator::initFieldFromGRP(grpmag);
        if (!o2::base::GeometryManager::isGeometryLoaded())
          fCCDB->get<TGeoManager>("GLO/Config/GeometryAligned");
        o2::mch::TrackExtrap::setField();
      }
    }

    auto nBcs = bcs.size();
    std::vector<uint64_t> vGlobalBCs(nBcs, 0);
    for (const auto& bc : bcs) {
      vGlobalBCs[bc.globalIndex()] = bc.globalBC();
    }

    auto nCols = collisions.size();
    std::vector<int64_t> vColIndexBCs(nCols, -1);
    for (const auto& col : collisions) {
      vColIndexBCs[col.globalIndex()] = col.bcId();
    }

    std::map<uint64_t, int64_t> mapGlobalBcWithV0A{};
    constexpr float FV0ValidTime = 15.f;
    for (const auto& fv0 : fv0s) {
      if (std::abs(fv0.time()) > FV0ValidTime)
        continue;
      uint64_t globalBC = vGlobalBCs[fv0.bcId()];
      mapGlobalBcWithV0A[globalBC] = fv0.globalIndex();
    }
    auto nFV0s = mapGlobalBcWithV0A.size();

    std::map<uint64_t, int64_t> mapGlobalBcWithZdc{};
    constexpr float ZDCValidTime = 2.f;
    for (const auto& zdc : zdcs) {
      if (std::abs(zdc.timeZNA()) > ZDCValidTime && std::abs(zdc.timeZNC()) > ZDCValidTime)
        continue;
      uint64_t globalBC = vGlobalBCs[zdc.bcId()];
      mapGlobalBcWithZdc[globalBC] = zdc.globalIndex();
    }
    auto nZdcs = mapGlobalBcWithZdc.size();

    auto nFwdTracks = fwdTracks.size();
    auto nAmbFwdTracks = ambFwdTracks.size();
    std::vector<int64_t> vAmbFwdTrackIndex(nFwdTracks, -1);
    std::vector<int64_t> vAmbFwdTrackIndexBCs(nAmbFwdTracks, -1);
    for (const auto& ambTr : ambFwdTracks) {
      vAmbFwdTrackIndex[ambTr.fwdtrackId()] = ambTr.globalIndex();
      vAmbFwdTrackIndexBCs[ambTr.globalIndex()] = ambTr.bcIds()[0];
    }

    auto nBarrelTracks = barrelTracks.size();
    auto nAmbBarrelTracks = ambBarrelTracks.size();
    std::vector<int64_t> vAmbBarTrackIndex(nBarrelTracks, -1);
    std::vector<int64_t> vAmbBarTrackIndexBCs(nAmbBarrelTracks, -1);
    for (const auto& ambTr : ambBarrelTracks) {
      vAmbBarTrackIndex[ambTr.trackId()] = ambTr.globalIndex();
      vAmbBarTrackIndexBCs[ambTr.globalIndex()] = ambTr.bcIds()[0];
    }

    // anchor map: MCH-MID forward muons keyed by their global BC
    std::map<uint64_t, std::vector<int64_t>> mapGlobalBcsWithMCHMIDTrackIds;
    for (const auto& fwdTrack : fwdTracks) {
      if (fwdTrack.trackType() != MuonStandaloneTrack)
        continue;
      auto trackId = fwdTrack.globalIndex();
      int64_t indexBC = vAmbFwdTrackIndex[trackId] < 0 ? vColIndexBCs[fwdTrack.collisionId()] : vAmbFwdTrackIndexBCs[vAmbFwdTrackIndex[trackId]];
      auto globalBC = vGlobalBCs[indexBC] + TMath::FloorNint(fwdTrack.trackTime() / o2::constants::lhc::LHCBunchSpacingNS + 1.);
      mapGlobalBcsWithMCHMIDTrackIds[globalBC].push_back(trackId);
    }

    // companion map: barrel tracks keyed by their global BC
    std::map<uint64_t, std::vector<int64_t>> mapGlobalBcsWithBarrelTrackIds;
    for (const auto& barTrack : barrelTracks) {
      auto trackId = barTrack.globalIndex();
      int64_t bcId = vAmbBarTrackIndex[trackId] < 0
                       ? (barTrack.has_collision() ? vColIndexBCs[barTrack.collisionId()] : -1)
                       : vAmbBarTrackIndexBCs[vAmbBarTrackIndex[trackId]];
      if (bcId < 0)
        continue;
      auto globalBC = vGlobalBCs[bcId] + TMath::FloorNint(barTrack.trackTime() / o2::constants::lhc::LHCBunchSpacingNS + 1.);
      mapGlobalBcsWithBarrelTrackIds[globalBC].push_back(trackId);
    }

    int32_t candId = 0;
    for (const auto& gbcMuids : mapGlobalBcsWithMCHMIDTrackIds) {
      uint64_t globalBcMid = gbcMuids.first;

      // veto on FV0 in the same BC as the anchor
      auto itFv0Id = mapGlobalBcWithV0A.find(globalBcMid);
      if (itFv0Id != mapGlobalBcWithV0A.end()) {
        auto fv0Id = itFv0Id->second;
        const auto& fv0 = fv0s.iteratorAt(fv0Id);
        float fv0Amp = 0.f;
        for (const auto& amp : fv0.amplitude())
          fv0Amp += amp;
        if (fv0Amp > fMaxFV0Amp)
          continue;
      }

      uint16_t numContrib = 0;
      auto& vMuonIds = gbcMuids.second;
      for (const auto& imuon : vMuonIds) {
        if (!addToFwdTable(candId, imuon, globalBcMid, 0., fwdTracks, mcFwdTrackLabels))
          continue;
        numContrib++;
      }
      if (numContrib < 1) // anchor produced no surviving fwd tracks
        continue;

      // gather barrel tracks within ±fBcWindowBarrel of the anchor
      std::map<int64_t, uint64_t> mapBarIdBc{};
      getBarrelTrackIds(globalBcMid, mapGlobalBcsWithBarrelTrackIds, fBcWindowBarrel, mapBarIdBc);
      for (const auto& [ibar, gbc] : mapBarIdBc) {
        float trackTime = (static_cast<int64_t>(gbc) - static_cast<int64_t>(globalBcMid)) * o2::constants::lhc::LHCBunchSpacingNS;
        if (!addToBarrelTable(candId, ibar, gbc, trackTime, barrelTracks, mcBarrelTrackLabels))
          continue;
        numContrib++;
      }

      eventCandidates(globalBcMid, runNumber, 0., 0., 0., 0, numContrib, 0, 0);
      std::vector<float> amplitudesV0A{};
      std::vector<int8_t> relBCsV0A{};
      std::vector<float> amplitudesT0A{};
      std::vector<int8_t> relBCsT0A{};
      if (nFV0s > 0) {
        getFV0Amplitudes(globalBcMid, fv0s, fBcWindowFITAmps, mapGlobalBcWithV0A, amplitudesV0A, relBCsV0A);
      }
      eventCandidatesSelsFwd(0., 0., amplitudesT0A, relBCsT0A, amplitudesV0A, relBCsV0A);
      if (nZdcs > 0) {
        auto itZDC = mapGlobalBcWithZdc.find(globalBcMid);
        if (itZDC != mapGlobalBcWithZdc.end()) {
          const auto& zdc = zdcs.iteratorAt(itZDC->second);
          float timeZNA = zdc.timeZNA();
          float timeZNC = zdc.timeZNC();
          float eComZNA = zdc.energyCommonZNA();
          float eComZNC = zdc.energyCommonZNC();
          udZdcsReduced(candId, timeZNA, timeZNC, eComZNA, eComZNC);
        }
      }
      candId++;
    }

    vGlobalBCs.clear();
    vColIndexBCs.clear();
    mapGlobalBcWithV0A.clear();
    mapGlobalBcWithZdc.clear();
    vAmbFwdTrackIndex.clear();
    vAmbFwdTrackIndexBCs.clear();
    vAmbBarTrackIndex.clear();
    vAmbBarTrackIndexBCs.clear();
    mapGlobalBcsWithMCHMIDTrackIds.clear();
    mapGlobalBcsWithBarrelTrackIds.clear();
  }

  void processSemiFwd(ForwardTracks const& fwdTracks,
                      o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                      BarrelTracks const& barrelTracks,
                      o2::aod::AmbiguousTracks const& ambBarrelTracks,
                      o2::aod::BCs const& bcs,
                      o2::aod::Collisions const& collisions,
                      o2::aod::FV0As const& fv0s,
                      o2::aod::Zdcs const& zdcs)
  {
    fDoMC = false;
    createCandidates(fwdTracks, ambFwdTracks, barrelTracks, ambBarrelTracks, bcs, collisions, fv0s, zdcs,
                     (o2::aod::McFwdTrackLabels*)nullptr, (o2::aod::McTrackLabels*)nullptr);
  }

  void processSemiFwdMC(ForwardTracks const& fwdTracks,
                        o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                        BarrelTracks const& barrelTracks,
                        o2::aod::AmbiguousTracks const& ambBarrelTracks,
                        o2::aod::BCs const& bcs,
                        o2::aod::Collisions const& collisions,
                        o2::aod::FV0As const& fv0s,
                        o2::aod::Zdcs const& zdcs,
                        o2::aod::McCollisions const& mcCollisions,
                        o2::aod::McParticles const& mcParticles,
                        o2::aod::McFwdTrackLabels const& mcFwdTrackLabels,
                        o2::aod::McTrackLabels const& mcBarrelTrackLabels)
  {
    fDoMC = true;
    skimMCInfo(mcCollisions, mcParticles);
    createCandidates(fwdTracks, ambFwdTracks, barrelTracks, ambBarrelTracks, bcs, collisions, fv0s, zdcs,
                     &mcFwdTrackLabels, &mcBarrelTrackLabels);
    fNewPartIDs.clear();
  }

  PROCESS_SWITCH(UpcCandProducerSemiFwd, processSemiFwd, "Produce candidates combining forward MCH-MID and central-barrel tracks", true);
  PROCESS_SWITCH(UpcCandProducerSemiFwd, processSemiFwdMC, "Produce semi-forward candidates with MC information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UpcCandProducerSemiFwd>(cfgc)};
}
