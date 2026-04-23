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

/// \file   upcCandProducerGlobalMuon.cxx
/// \brief  UPC candidate producer for forward muons with MFT support - UNDER DEVELOPMENT
/// \author Roman Lavicka, roman.lavicka@cern.ch
/// \since  11.02.2026

#include "PWGUD/Core/UPCCutparHolder.h"
#include "PWGUD/Core/UPCHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

#include "Common/Core/fwdtrackUtilities.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/FwdDCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Field/MagneticField.h>
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
#include <ReconstructionDataFormats/TrackFwd.h>

#include <Math/MatrixRepresentationsStatic.h>
#include <Math/SMatrix.h>
#include <TGeoGlobalMagField.h>
#include <TH1.h>
#include <TMath.h>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcCandProducerGlobalMuon {
  bool fDoMC{false};

  std::map<int32_t, int32_t> fNewPartIDs;
  std::map<uint32_t, bool> fBestMuonMatch;

  Produces<o2::aod::UDMcCollisions> udMCCollisions;
  Produces<o2::aod::UDMcParticles> udMCParticles;
  Produces<o2::aod::UDMcFwdTrackLabels> udFwdTrackLabels;

  Produces<o2::aod::UDFwdTracks> udFwdTracks;
  Produces<o2::aod::UDFwdTracksExtra> udFwdTracksExtra;
  Produces<o2::aod::UDFwdIndices> udFwdIndices;
  Produces<o2::aod::UDFwdTracksCls> udFwdTrkClusters; // Added for MFT clusters
  Produces<o2::aod::UDCollisions> eventCandidates;
  Produces<o2::aod::UDCollisionsSelsFwd> eventCandidatesSelsFwd;
  Produces<o2::aod::UDZdcsReduced> udZdcsReduced;

  Configurable<int> fSignalGenID{"fSignalGenID", 1, "Signal generator ID"};

  UPCCutparHolder fUpcCuts = UPCCutparHolder();
  MutableConfigurable<UPCCutparHolder> fUpcCutsConf{"fUpcCutsConf", {}, "UPC event cuts"};
  Configurable<uint64_t> fBcWindowFITAmps{"fBcWindowFITAmps", 20, "BC range for T0A/V0A amplitudes array [-range, +(range-1)]"};
  Configurable<int> fBcWindowMCH{"fBcWindowMCH", 20, "Time window for MCH-MID to MCH-only matching for Muon candidates"};
  Configurable<float> fMaxFV0Amp{"fMaxFV0Amp", 100.f, "Max FV0 amplitude in the same BC"};

  // NEW: MFT/Global track support configurables
  Configurable<bool> fEnableMFT{"fEnableMFT", true, "Enable MFT/global track processing"};
  Configurable<bool> fSaveMFTClusters{"fSaveMFTClusters", true, "Save MFT cluster information"};

  // FwdDCAFitter configurables
  Configurable<bool> fPropagateToPCA{"fPropagateToPCA", true, "Propagate tracks to PCA"};
  Configurable<float> fMaxR{"fMaxR", 200.f, "Maximum radius for FwdDCAFitter (cm)"};
  Configurable<float> fMinParamChange{"fMinParamChange", 1.0e-3, "Minimum parameter change for FwdDCAFitter convergence"};
  Configurable<float> fMinRelChi2Change{"fMinRelChi2Change", 0.9, "Minimum relative chi2 change for FwdDCAFitter convergence"};
  Configurable<bool> fUseAbsDCA{"fUseAbsDCA", true, "Use absolute DCA in FwdDCAFitter"};
  Configurable<float> fMaxChi2MatchMCHMFT{"fMaxChi2MatchMCHMFT", 4.f, "Maximum chi2 for MCH-MFT matching quality filter"};
  Configurable<int> fBcWindowMCHMFT{"fBcWindowMCHMFT", 20, "BC window for searching MCH-MFT tracks around MCH-MID-MFT anchors"};
  Configurable<bool> fKeepBestMuonMatch{"fKeepBestMuonMatch", true, "Keep only the best MCH-MFT match per MCH track (lowest chi2)"};

  using ForwardTracks = o2::soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov>;

  HistogramRegistry histRegistry{"HistRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  int fRun{0};
  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::ccdb::CcdbApi fCCDBApi;
  o2::globaltracking::MatchGlobalFwd fMatching;

  // FwdDCAFitter member
  o2::vertexing::FwdDCAFitterN<2> fFwdFitter;

  // Named constants (avoid magic numbers in expressions)
  static constexpr double kBcTimeRoundingOffset = 1.;     // Offset used when rounding trackTime to BC units
  static constexpr uint16_t kMinTracksForPair = 2;        // Minimum tracks required to compute a pair invariant mass
  static constexpr uint16_t kMinTracksForCandidate = 1;   // Minimum contributors required to save a candidate
  static constexpr int kUpperBoundaryToTrackTypeEnum = 2; // Make sure you use MFT tracks

  void init(InitContext&)
  {
    fUpcCuts = (UPCCutparHolder)fUpcCutsConf;
    if (fUpcCuts.getUseFwdCuts() || fEnableMFT) {
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

    // NEW: Add histograms for global track monitoring
    const AxisSpec axisTrackType{5, -0.5, 4.5, "Track Type"};
    histRegistry.add("hTrackTypes", "Track type distribution", kTH1F, {axisTrackType});
    histRegistry.get<TH1>(HIST("hTrackTypes"))->GetXaxis()->SetBinLabel(1, "MuonStandalone");
    histRegistry.get<TH1>(HIST("hTrackTypes"))->GetXaxis()->SetBinLabel(2, "MCHStandalone");
    histRegistry.get<TH1>(HIST("hTrackTypes"))->GetXaxis()->SetBinLabel(3, "GlobalMuon");
    histRegistry.get<TH1>(HIST("hTrackTypes"))->GetXaxis()->SetBinLabel(4, "GlobalFwd");

    const AxisSpec axisEta{100, -4.0, -2.0, "#eta"};
    histRegistry.add("hEtaGlobal", "Global track eta", kTH1F, {axisEta});

    const AxisSpec axisChi2PCA{200, 0., 100., "#chi^{2}_{PCA}"};
    histRegistry.add("hChi2PCA", "Chi2 at PCA from FwdDCAFitter", kTH1F, {axisChi2PCA});

    const AxisSpec axisChi2Match{200, 0., 100., "#chi^{2}_{MCH-MFT}"};
    histRegistry.add("hChi2MatchMCHMFT", "Chi2 of MCH-MFT matching (before cut)", kTH1F, {axisChi2Match});

    const AxisSpec axisMass{500, 0., 10., "m_{inv} (GeV/c^{2})"};
    histRegistry.add("hMassGlobalMuon", "Invariant mass from MCH-MID-MFT tracks only", kTH1F, {axisMass});
    histRegistry.add("hMassGlobalMuonWithMCHMFT", "Invariant mass from MCH-MID-MFT + MCH-MFT tracks", kTH1F, {axisMass});
  }

  bool cut(const o2::dataformats::GlobalFwdTrack& pft, const ForwardTracks::iterator& fwdTrack)
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

  void skimMCInfo(o2::aod::McCollisions const& mcCollisions, o2::aod::McParticles const& mcParticles)
  {
    std::vector<int32_t> newEventIDs(mcCollisions.size(), -1);

    int32_t newPartID = 0;
    int32_t newEventID = 0;
    int32_t nMCParticles = mcParticles.size();
    // loop over MC particles to select only the ones from signal events
    // and calculate new MC table IDs
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

    // storing MC particles
    for (const auto& item : fNewPartIDs) {
      int32_t mcPartID = item.first;
      const auto& mcPart = mcParticles.iteratorAt(mcPartID);
      int32_t mcEventID = mcPart.mcCollisionId();
      int32_t newEventID = newEventIDs[mcEventID];
      // collecting new mother IDs
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
      // collecting new daughter IDs
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
      udMCParticles(newEventID, mcPart.pdgCode(), mcPart.getHepMCStatusCode(), mcPart.flags(), newMotherIDs, newDaughterIDs,
                    mcPart.weight(), mcPart.px(), mcPart.py(), mcPart.pz(), mcPart.e());
      newMotherIDs.clear();
    }

    // storing MC events
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

  // find starting point for scrolling in some BC map -- closest to the input gbc
  template <typename T>
  T::iterator getStartForScroll(uint64_t inGbc, T& gbcMap)
  {
    auto it1 = gbcMap.lower_bound(inGbc);
    typename T::iterator it;
    if (it1 != gbcMap.end()) { // found lower bound
      auto it2 = it1;
      uint64_t bc1 = it1->first;
      if (it2 != gbcMap.begin())
        --it2;
      uint64_t bc2 = it2->first;
      uint64_t dbc1 = bcDiff(bc1, inGbc);
      uint64_t dbc2 = bcDiff(bc2, inGbc);
      it = (dbc1 <= dbc2) ? it1 : it2;
    } else { // ended up in the end
      it = it1;
      --it;
    }
    return it;
  }

  // scroll over gbcMap and do some operation on container
  template <typename T, typename F>
  void scrollBackForth(uint64_t inGbc, uint64_t maxDbc, T& gbcMap, F&& func)
  {
    auto it = getStartForScroll(inGbc, gbcMap);
    uint64_t gbc = it->first;
    uint64_t dbc = bcDiff(inGbc, gbc);

    // start scrolling backward
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

    std::advance(it, count + 1); // move back to the starting point + 1

    if (it == gbcMap.end()) // ended up in the end of map
      return;

    gbc = it->first;
    dbc = bcDiff(inGbc, gbc);

    // start scrolling forward
    while (dbc <= maxDbc) {
      func(it, gbc);
      ++it;
      if (it == gbcMap.end())
        break;
      gbc = it->first;
      dbc = bcDiff(inGbc, gbc);
    }
  }

  void getMchTrackIds(uint64_t inGbc, std::map<uint64_t, std::vector<int64_t>>& mchTracksPerBC,
                      uint64_t maxDbc, std::map<int64_t, uint64_t>& outMchTrkIds)
  {
    auto fillMchIds = [&outMchTrkIds](std::map<uint64_t, std::vector<int64_t>>::iterator& inIt, uint64_t gbc) {
      std::vector<int64_t>& ids = inIt->second;
      for (const auto& id : ids)
        outMchTrkIds[id] = gbc;
    };
    scrollBackForth(inGbc, maxDbc, mchTracksPerBC, fillMchIds);
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

  bool addToFwdTable(int64_t candId, int64_t trackId, uint64_t gbc, float trackTime,
                     ForwardTracks const& fwdTracks, o2::aod::MFTTracks const& mftTracks,
                     const o2::aod::McFwdTrackLabels* mcFwdTrackLabels)
  {
    const auto& track = fwdTracks.iteratorAt(trackId);
    float px, py, pz;
    int sign;

    // Fill track type histogram
    histRegistry.fill(HIST("hTrackTypes"), track.trackType());
    if (track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack ||
        track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalForwardTrack) {
      histRegistry.fill(HIST("hEtaGlobal"), track.eta());
    }

    if (fUpcCuts.getUseFwdCuts()) {
      bool isGlobal = track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack ||
                      track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalForwardTrack;
      o2::dataformats::GlobalFwdTrack pft;
      if (isGlobal && fEnableMFT) {
        // Refit global muon: propagate MCH component to vertex, combine with MFT spatial info
        auto mchTrack = track.matchMCHTrack_as<ForwardTracks>();
        auto propMuon = propagateToZero(mchTrack);
        auto mfttrack = track.matchMFTTrack_as<o2::aod::MFTTracks>();
        using SMatrix5 = ROOT::Math::SVector<double, 5>;
        using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
        SMatrix5 tpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
        SMatrix55 tcovs{};
        o2::track::TrackParCovFwd mft{mfttrack.z(), tpars, tcovs, mfttrack.chi2()};
        pft = o2::aod::fwdtrackutils::refitGlobalMuonCov(propMuon, mft);
      } else {
        pft = propagateToZero(track);
      }
      bool pass = cut(pft, track);
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

    // NEW: Enhanced extra info for global tracks
    float mchmftChi2 = -1.;
    if (track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack ||
        track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalForwardTrack) {
      mchmftChi2 = track.chi2MatchMCHMFT();
    }

    udFwdTracksExtra(track.trackType(), track.nClusters(), track.pDca(), track.rAtAbsorberEnd(),
                     track.chi2(), track.chi2MatchMCHMID(), mchmftChi2,
                     track.mchBitMap(), track.midBitMap(), track.midBoards());

    // NEW: Store MFT index for global tracks
    int64_t mftIndex = -1;
    if (fEnableMFT && (track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack ||
                       track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalForwardTrack)) {
      mftIndex = track.matchMFTTrackId();
    }

    udFwdIndices(candId, trackId, track.matchMCHTrackId(), mftIndex);

    if (fDoMC) {
      const auto& label = mcFwdTrackLabels->iteratorAt(trackId);
      uint16_t mcMask = label.mcMask();
      auto it = fNewPartIDs.find(label.mcParticleId());
      int32_t newPartID = it != fNewPartIDs.end() ? it->second : -1;
      udFwdTrackLabels(newPartID, mcMask);
    }
    return true;
  }

  // NEW: Function to fill MFT cluster information
  void fillFwdClusters(const std::vector<int64_t>& trackIds, o2::aod::FwdTrkCls const& fwdTrkCls)
  {
    if (!fSaveMFTClusters)
      return;

    std::map<int64_t, std::vector<int64_t>> clustersPerTrack;
    for (const auto& cls : fwdTrkCls) {
      clustersPerTrack[cls.fwdtrackId()].push_back(cls.globalIndex());
    }

    int newId = 0;
    for (const auto& trackId : trackIds) {
      auto it = clustersPerTrack.find(trackId);
      if (it != clustersPerTrack.end()) {
        const auto& clusters = it->second;
        for (const auto& clsId : clusters) {
          const auto& clsInfo = fwdTrkCls.iteratorAt(clsId);
          udFwdTrkClusters(newId, clsInfo.x(), clsInfo.y(), clsInfo.z(), clsInfo.clInfo());
        }
      }
      newId++;
    }
  }

  // Convert forward track to TrackParCovFwd for the FwdDCAFitter
  o2::track::TrackParCovFwd fwdToTrackPar(ForwardTracks::iterator const& track)
  {
    using SMatrix5 = ROOT::Math::SVector<double, 5>;
    using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
    SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
    std::vector<double> v1{track.cXX(), track.cXY(), track.cYY(), track.cPhiX(), track.cPhiY(),
                           track.cPhiPhi(), track.cTglX(), track.cTglY(), track.cTglPhi(), track.cTglTgl(),
                           track.c1PtX(), track.c1PtY(), track.c1PtPhi(), track.c1PtTgl(), track.c1Pt21Pt2()};
    SMatrix55 tcovs(v1.begin(), v1.end());
    return o2::track::TrackParCovFwd{track.z(), tpars, tcovs, track.chi2()};
  }

  // Select the best MCH-MFT match per MCH track based on lowest chi2MatchMCHMFT.
  // Multiple global tracks can share the same MCH track with different MFT matches;
  // this function keeps only the best one to reduce combinatorial background.
  void selectBestMuonMatches(ForwardTracks const& fwdTracks)
  {
    fBestMuonMatch.clear();
    std::unordered_map<int, std::pair<float, int>> mCandidates;
    for (const auto& muon : fwdTracks) {
      if (static_cast<int>(muon.trackType()) < kUpperBoundaryToTrackTypeEnum) {
        auto muonID = muon.matchMCHTrackId();
        auto chi2 = muon.chi2MatchMCHMFT();
        if (mCandidates.find(muonID) == mCandidates.end()) {
          mCandidates[muonID] = {chi2, muon.globalIndex()};
        } else {
          if (chi2 < mCandidates[muonID].first) {
            mCandidates[muonID] = {chi2, muon.globalIndex()};
          }
        }
      }
    }
    for (const auto& pairCand : mCandidates) {
      fBestMuonMatch[pairCand.second.second] = true;
    }
  }

  void createCandidates(ForwardTracks const& fwdTracks,
                        o2::aod::MFTTracks const& mftTracks,
                        o2::aod::FwdTrkCls const& fwdTrkCls,
                        o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                        o2::aod::BCs const& bcs,
                        o2::aod::Collisions const& collisions,
                        o2::aod::FV0As const& fv0s,
                        o2::aod::Zdcs const& zdcs,
                        const o2::aod::McFwdTrackLabels* mcFwdTrackLabels)
  {
    using o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalForwardTrack;
    using o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack;
    using o2::aod::fwdtrack::ForwardTrackTypeEnum::MCHStandaloneTrack;
    using o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack;

    int32_t runNumber = bcs.iteratorAt(0).runNumber();
    if (fUpcCuts.getUseFwdCuts() || fEnableMFT) {
      if (runNumber != fRun) {
        fRun = runNumber;
        std::map<std::string, std::string> metadata;
        auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(fCCDBApi, fRun);
        auto ts = soreor.second;
        auto grpmag = fCCDBApi.retrieveFromTFileAny<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", metadata, ts);
        o2::base::Propagator::initFieldFromGRP(grpmag);
        if (!o2::base::GeometryManager::isGeometryLoaded())
          fCCDB->get<TGeoManager>("GLO/Config/GeometryAligned");
        if (fUpcCuts.getUseFwdCuts()) {
          o2::mch::TrackExtrap::setField();
        }
        // Initialize FwdDCAFitter with magnetic field
        if (fEnableMFT) {
          o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
          static constexpr double kCenterMFT[3] = {0, 0, -61.4};
          float bz = field->getBz(kCenterMFT);
          LOG(info) << "FwdDCAFitter magnetic field at MFT center: bZ = " << bz;
          fFwdFitter.setBz(bz);
          fFwdFitter.setPropagateToPCA(fPropagateToPCA);
          fFwdFitter.setMaxR(fMaxR);
          fFwdFitter.setMinParamChange(fMinParamChange);
          fFwdFitter.setMinRelChi2Change(fMinRelChi2Change);
          fFwdFitter.setUseAbsDCA(fUseAbsDCA);
        }
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

    // Select best MCH-MFT match per MCH track before sorting into BC maps
    if (fKeepBestMuonMatch) {
      selectBestMuonMatches(fwdTracks);
    }

    std::map<uint64_t, std::vector<int64_t>> mapGlobalBcsWithMCHMIDTrackIds;
    std::map<uint64_t, std::vector<int64_t>> mapGlobalBcsWithMCHTrackIds;
    std::map<uint64_t, std::vector<int64_t>> mapGlobalBcsWithGlobalMuonTrackIds; // MCH-MID-MFT (good timing from MID)
    std::map<uint64_t, std::vector<int64_t>> mapGlobalBcsWithMCHMFTTrackIds;     // MCH-MFT only (poor timing)

    for (const auto& fwdTrack : fwdTracks) {
      auto trackType = fwdTrack.trackType();

      // Skip if not a relevant track type
      if (trackType != MCHStandaloneTrack &&
          trackType != MuonStandaloneTrack &&
          trackType != GlobalMuonTrack &&
          trackType != GlobalForwardTrack)
        continue;

      // For global tracks, skip if not the best match for this MCH track
      if (fKeepBestMuonMatch && static_cast<int>(trackType) < kUpperBoundaryToTrackTypeEnum) {
        if (fBestMuonMatch.find(fwdTrack.globalIndex()) == fBestMuonMatch.end()) {
          continue;
        }
      }

      auto trackId = fwdTrack.globalIndex();
      int64_t indexBC = vAmbFwdTrackIndex[trackId] < 0 ? vColIndexBCs[fwdTrack.collisionId()] : vAmbFwdTrackIndexBCs[vAmbFwdTrackIndex[trackId]];
      auto globalBC = vGlobalBCs[indexBC] + TMath::FloorNint(fwdTrack.trackTime() / o2::constants::lhc::LHCBunchSpacingNS + kBcTimeRoundingOffset);

      if (trackType == MuonStandaloneTrack) { // MCH-MID
        mapGlobalBcsWithMCHMIDTrackIds[globalBC].push_back(trackId);
      } else if (trackType == MCHStandaloneTrack) { // MCH-only
        mapGlobalBcsWithMCHTrackIds[globalBC].push_back(trackId);
      } else if (trackType == GlobalMuonTrack) { // MCH-MID-MFT: good timing, used as anchor
        histRegistry.fill(HIST("hChi2MatchMCHMFT"), fwdTrack.chi2MatchMCHMFT());
        if (fwdTrack.chi2MatchMCHMFT() > 0 && fwdTrack.chi2MatchMCHMFT() < fMaxChi2MatchMCHMFT) {
          mapGlobalBcsWithGlobalMuonTrackIds[globalBC].push_back(trackId);
        }
      } else if (trackType == GlobalForwardTrack) { // MCH-MFT: poor timing, matched to anchors
        histRegistry.fill(HIST("hChi2MatchMCHMFT"), fwdTrack.chi2MatchMCHMFT());
        if (fwdTrack.chi2MatchMCHMFT() > 0 && fwdTrack.chi2MatchMCHMFT() < fMaxChi2MatchMCHMFT) {
          mapGlobalBcsWithMCHMFTTrackIds[globalBC].push_back(trackId);
        }
      }
    }

    std::vector<int64_t> selTrackIds{}; // For cluster saving

    int32_t candId = 0;

    // Process global tracks: MCH-MID-MFT anchors + MCH-MFT in BC window
    // MCH-MID-MFT tracks have good timing (from MID) and serve as anchors.
    // MCH-MFT tracks have poor timing and are searched in a BC window around anchors.
    if (!mapGlobalBcsWithGlobalMuonTrackIds.empty()) {
      for (const auto& gbc_anchorids : mapGlobalBcsWithGlobalMuonTrackIds) {
        uint64_t globalBcAnchor = gbc_anchorids.first;
        auto itFv0Id = mapGlobalBcWithV0A.find(globalBcAnchor);
        if (itFv0Id != mapGlobalBcWithV0A.end()) {
          auto fv0Id = itFv0Id->second;
          const auto& fv0 = fv0s.iteratorAt(fv0Id);
          float fv0Amp = 0.f;
          for (const auto& amp : fv0.amplitude())
            fv0Amp += amp;
          if (fv0Amp > fMaxFV0Amp)
            continue;
        }

        auto& vAnchorIds = gbc_anchorids.second; // MCH-MID-MFT tracks at this BC

        // Search MCH-MFT tracks in BC window around anchor (analogous to MCH-only matching)
        std::map<int64_t, uint64_t> mapMchMftIdBc{};
        getMchTrackIds(globalBcAnchor, mapGlobalBcsWithMCHMFTTrackIds, fBcWindowMCHMFT, mapMchMftIdBc);

        // Collect all track IDs (anchors + matched MCH-MFT)
        std::vector<int64_t> allTrackIds;
        allTrackIds.reserve(vAnchorIds.size() + mapMchMftIdBc.size());
        for (const auto& id : vAnchorIds)
          allTrackIds.push_back(id);
        for (const auto& [id, gbc] : mapMchMftIdBc)
          allTrackIds.push_back(id);

        // Step 1: Find secondary vertex using FwdDCAFitter on the first track pair
        float bestVtxX = 0., bestVtxY = 0., bestVtxZ = 0.;

        if (allTrackIds.size() >= kMinTracksForPair) {
          const auto& trk1 = fwdTracks.iteratorAt(allTrackIds[0]);
          const auto& trk2 = fwdTracks.iteratorAt(allTrackIds[1]);
          auto pars1 = fwdToTrackPar(trk1);
          auto pars2 = fwdToTrackPar(trk2);
          int procCode = fFwdFitter.process(pars1, pars2);
          if (procCode != 0) {
            auto secondaryVertex = fFwdFitter.getPCACandidate();
            auto chi2PCA = fFwdFitter.getChi2AtPCACandidate();
            bestVtxX = secondaryVertex[0];
            bestVtxY = secondaryVertex[1];
            bestVtxZ = secondaryVertex[2];
            histRegistry.fill(HIST("hChi2PCA"), chi2PCA);
          }
        }

        // Step 2: Write anchor tracks (MCH-MID-MFT)
        constexpr double kMuonMass = o2::constants::physics::MassMuon;
        uint16_t numContrib = 0;
        double sumPx = 0., sumPy = 0., sumPz = 0., sumE = 0.;
        for (const auto& ianchor : vAnchorIds) {
          if (!addToFwdTable(candId, ianchor, globalBcAnchor, 0., fwdTracks, mftTracks, mcFwdTrackLabels))
            continue;
          const auto& trk = fwdTracks.iteratorAt(ianchor);
          double p2 = trk.px() * trk.px() + trk.py() * trk.py() + trk.pz() * trk.pz();
          sumPx += trk.px();
          sumPy += trk.py();
          sumPz += trk.pz();
          sumE += std::sqrt(p2 + kMuonMass * kMuonMass);
          numContrib++;
          selTrackIds.push_back(ianchor);
        }

        // Fill invariant mass from MCH-MID-MFT anchors only
        uint16_t numContribAnch = numContrib;
        if (numContribAnch >= kMinTracksForPair) {
          double mass2 = sumE * sumE - sumPx * sumPx - sumPy * sumPy - sumPz * sumPz;
          histRegistry.fill(HIST("hMassGlobalMuon"), mass2 > 0. ? std::sqrt(mass2) : 0.);
        }

        // Step 3: Write matched MCH-MFT tracks with adjusted track time
        for (const auto& [imchMft, gbc] : mapMchMftIdBc) {
          if (!addToFwdTable(candId, imchMft, gbc, (gbc - globalBcAnchor) * o2::constants::lhc::LHCBunchSpacingNS, fwdTracks, mftTracks, mcFwdTrackLabels))
            continue;
          const auto& trk = fwdTracks.iteratorAt(imchMft);
          double p2 = trk.px() * trk.px() + trk.py() * trk.py() + trk.pz() * trk.pz();
          sumPx += trk.px();
          sumPy += trk.py();
          sumPz += trk.pz();
          sumE += std::sqrt(p2 + kMuonMass * kMuonMass);
          numContrib++;
          selTrackIds.push_back(imchMft);
        }

        // Fill invariant mass including MCH-MFT tracks (only if MCH-MFT tracks were added)
        if (numContrib > numContribAnch && numContrib >= kMinTracksForPair) {
          double mass2 = sumE * sumE - sumPx * sumPx - sumPy * sumPy - sumPz * sumPz;
          histRegistry.fill(HIST("hMassGlobalMuonWithMCHMFT"), mass2 > 0. ? std::sqrt(mass2) : 0.);
        }

        if (numContrib < kMinTracksForCandidate)
          continue;

        eventCandidates(globalBcAnchor, runNumber, bestVtxX, bestVtxY, bestVtxZ, 0, numContrib, 0, 0);
        std::vector<float> amplitudesV0A{};
        std::vector<int8_t> relBCsV0A{};
        std::vector<float> amplitudesT0A{};
        std::vector<int8_t> relBCsT0A{};
        if (nFV0s > 0) {
          getFV0Amplitudes(globalBcAnchor, fv0s, fBcWindowFITAmps, mapGlobalBcWithV0A, amplitudesV0A, relBCsV0A);
        }
        eventCandidatesSelsFwd(0., 0., amplitudesT0A, relBCsT0A, amplitudesV0A, relBCsV0A);
        if (nZdcs > 0) {
          auto itZDC = mapGlobalBcWithZdc.find(globalBcAnchor);
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
    }

    // Process MCH-MID tracks (original logic)
    for (const auto& gbc_muids : mapGlobalBcsWithMCHMIDTrackIds) {
      uint64_t globalBcMid = gbc_muids.first;
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
      auto& vMuonIds = gbc_muids.second;
      // writing MCH-MID tracks
      for (const auto& imuon : vMuonIds) {
        if (!addToFwdTable(candId, imuon, globalBcMid, 0., fwdTracks, mftTracks, mcFwdTrackLabels))
          continue;
        numContrib++;
        selTrackIds.push_back(imuon);
      }
      if (numContrib < kMinTracksForCandidate) // didn't save any MCH-MID tracks
        continue;
      std::map<int64_t, uint64_t> mapMchIdBc{};
      getMchTrackIds(globalBcMid, mapGlobalBcsWithMCHTrackIds, fBcWindowMCH, mapMchIdBc);
      // writing MCH-only tracks
      for (const auto& [imch, gbc] : mapMchIdBc) {
        if (!addToFwdTable(candId, imch, gbc, (gbc - globalBcMid) * o2::constants::lhc::LHCBunchSpacingNS, fwdTracks, mftTracks, mcFwdTrackLabels))
          continue;
        numContrib++;
        selTrackIds.push_back(imch);
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

    // NEW: Fill MFT cluster information
    if (fEnableMFT && fSaveMFTClusters && !selTrackIds.empty()) {
      fillFwdClusters(selTrackIds, fwdTrkCls);
    }

    vGlobalBCs.clear();
    vColIndexBCs.clear();
    mapGlobalBcWithV0A.clear();
    mapGlobalBcWithZdc.clear();
    vAmbFwdTrackIndex.clear();
    vAmbFwdTrackIndexBCs.clear();
    mapGlobalBcsWithMCHMIDTrackIds.clear();
    mapGlobalBcsWithMCHTrackIds.clear();
    mapGlobalBcsWithGlobalMuonTrackIds.clear();
    mapGlobalBcsWithMCHMFTTrackIds.clear();
    selTrackIds.clear();
  }

  void processFwd(ForwardTracks const& fwdTracks,
                  o2::aod::MFTTracks const& mftTracks,
                  o2::aod::FwdTrkCls const& fwdTrkCls,
                  o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                  o2::aod::BCs const& bcs,
                  o2::aod::Collisions const& collisions,
                  o2::aod::FV0As const& fv0s,
                  o2::aod::Zdcs const& zdcs)
  {
    fDoMC = false;
    createCandidates(fwdTracks, mftTracks, fwdTrkCls, ambFwdTracks, bcs, collisions, fv0s, zdcs, (o2::aod::McFwdTrackLabels*)nullptr);
  }

  void processFwdMC(ForwardTracks const& fwdTracks,
                    o2::aod::MFTTracks const& mftTracks,
                    o2::aod::FwdTrkCls const& fwdTrkCls,
                    o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                    o2::aod::BCs const& bcs,
                    o2::aod::Collisions const& collisions,
                    o2::aod::FV0As const& fv0s,
                    o2::aod::Zdcs const& zdcs,
                    o2::aod::McCollisions const& mcCollisions,
                    o2::aod::McParticles const& mcParticles,
                    o2::aod::McFwdTrackLabels const& mcFwdTrackLabels)
  {
    fDoMC = true;
    skimMCInfo(mcCollisions, mcParticles);
    createCandidates(fwdTracks, mftTracks, fwdTrkCls, ambFwdTracks, bcs, collisions, fv0s, zdcs, &mcFwdTrackLabels);
    fNewPartIDs.clear();
  }

  PROCESS_SWITCH(UpcCandProducerGlobalMuon, processFwd, "Produce candidates in forward region", true);
  PROCESS_SWITCH(UpcCandProducerGlobalMuon, processFwdMC, "Produce candidates in forward region with MC information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UpcCandProducerGlobalMuon>(cfgc)};
}
