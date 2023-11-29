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

/// \brief write relevant information for dalitz ee analysis to an AO2D.root file. This file is then the only necessary input to perform pcm analysis.
/// \author daiki.sekihata@cern.ch

#include "Math/Vector4D.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

using MyTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TracksCov,
                           aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                           aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFbeta>;
using MyTracksMC = soa::Join<MyTracks, aod::McTrackLabels>;

struct skimmerPrimaryMuon {
  enum class EM_EEPairType : int {
    kULS = 0,
    kLSpp = +1,
    kLSnn = -1,
  };

  SliceCache cache;
  Preslice<aod::Tracks> perCol = o2::aod::track::collisionId;
  Produces<aod::EMPrimaryMuons> emprimarymuons;
  Produces<aod::EMPrimaryMuonsPrefilterBit> muon_pfb;

  // Configurables
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min. crossed rows"};
  Configurable<float> min_tpc_cr_findable_ratio{"min_tpc_cr_findable_ratio", 0.8, "min. TPC Ncr/Nf ratio"};
  Configurable<int> minitsncls{"minitsncls", 4, "min. number of ITS clusters"};
  Configurable<float> maxchi2tpc{"maxchi2tpc", 5.0, "max. chi2/NclsTPC"};
  Configurable<float> maxchi2its{"maxchi2its", 6.0, "max. chi2/NclsITS"};
  Configurable<float> minpt{"minpt", 0.05, "min pt for track"};
  Configurable<float> maxpt{"maxpt", 1.0, "max pt for track"};
  Configurable<float> maxeta{"maxeta", 0.9, "eta acceptance"};
  Configurable<float> dca_xy_max{"dca_xy_max", 1.0f, "max DCAxy in cm"};
  Configurable<float> dca_z_max{"dca_z_max", 1.0f, "max DCAz in cm"};
  Configurable<float> min_tpcdEdx{"min_tpcdEdx", 30.0, "min TPC dE/dx"};
  Configurable<float> max_tpcdEdx{"max_tpcdEdx", 1e+10, "max TPC dE/dx"};
  Configurable<float> maxPin{"maxPin", 1.0, "max pin for PID"};
  Configurable<float> maxPin_TPC{"maxPin_TPC", 0.2, "max pin for TPC pid only"};
  Configurable<float> maxTPCNsigmaMu_lowPin{"maxTPCNsigmaMu_lowPin", +4.0, "max. TPC n sigma for muon inclusion at low pin"};
  Configurable<float> maxTPCNsigmaMu_highPin{"maxTPCNsigmaMu_highPin", +4.0, "max. TPC n sigma for muon inclusion at high pin"};
  Configurable<float> maxTOFNsigmaMu_highPin{"maxTOFNsigmaMu_highPin", +4.0, "max. TOF n sigma for muon inclusion at high pin"};
  Configurable<float> maxTPCNsigmaEl{"maxTPCNsigmaEl", 1.0, "max. TPC n sigma for electron exclusion"};
  Configurable<float> maxTPCNsigmaPi_lowPin{"maxTPCNsigmaPi_lowPin", -2.0, "max. TPC n sigma for pion exclusion"};
  Configurable<float> maxTOFNsigmaPi{"maxTOFNsigmaPi", -2.0, "max. TPC n sigma for electron exclusion"};
  Configurable<float> maxMmumu{"maxMmumu", 1.1, "max. mee to store ee pairs"};
  Configurable<bool> storeLS{"storeLS", false, "flag to store LS pairs"};

  HistogramRegistry fRegistry{
    "fRegistry",
    {
      {"hNpairs", "hNpairs;pair type;Number of Pairs", {HistType::kTH1F, {{3, -1.5f, +1.5f}}}},
      {"Track/hTPCdEdx_Pin_before", "TPC dE/dx vs. p_{in};p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 1.f}, {200, 0.f, 200.f}}}},
      {"Track/hTPCdEdx_Pin_after", "TPC dE/dx vs. p_{in};p_{in} (GeV/c);TPC dE/dx", {HistType::kTH2F, {{1000, 0.f, 1.f}, {200, 0.f, 200.f}}}},
      {"Track/hTOFbeta_Pin_before", "TOF beta vs. p_{in};p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{1000, 0.f, 1.f}, {240, 0.6f, 1.2f}}}},
      {"Track/hTOFbeta_Pin_after", "TOF beta vs. p_{in};p_{in} (GeV/c);TOF #beta", {HistType::kTH2F, {{1000, 0.f, 1.f}, {240, 0.6f, 1.2f}}}},
      {"Track/hTPCNsigmaEl_before", "TPC n sigma e vs. p_{in};p_{in} (GeV/c);n #sigma_{e}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 1.f}, {100, -5.f, +5.f}}}},
      {"Track/hTOFNsigmaEl_before", "TOF n sigma e vs. p_{in};p_{in} (GeV/c);n #sigma_{e}^{TOF}", {HistType::kTH2F, {{1000, 0.f, 1.f}, {100, -5.f, +5.f}}}},
      {"Track/hTPCNsigmaMu_before", "TPC n sigma #mu vs. p_{in};p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 1.f}, {100, -5.f, +5.f}}}},
      {"Track/hTOFNsigmaMu_before", "TOF n sigma #mu vs. p_{in};p_{in} (GeV/c);n #sigma_{#mu}^{TOF}", {HistType::kTH2F, {{1000, 0.f, 1.f}, {100, -5.f, +5.f}}}},
      {"Track/hTPCNsigmaPi_before", "TPC n sigma #pi vs. p_{in};p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 1.f}, {100, -5.f, +5.f}}}},
      {"Track/hTOFNsigmaPi_before", "TOF n sigma #pi vs. p_{in};p_{in} (GeV/c);n #sigma_{#pi}^{TOF}", {HistType::kTH2F, {{1000, 0.f, 1.f}, {100, -5.f, +5.f}}}},
      {"Track/hTPCNsigmaEl_after", "TPC n sigma e vs. p_{in};p_{in} (GeV/c);n #sigma_{e}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 1.f}, {100, -5.f, +5.f}}}},
      {"Track/hTOFNsigmaEl_after", "TOF n sigma e vs. p_{in};p_{in} (GeV/c);n #sigma_{e}^{TOF}", {HistType::kTH2F, {{1000, 0.f, 1.f}, {100, -5.f, +5.f}}}},
      {"Track/hTPCNsigmaMu_after", "TPC n sigma #mu vs. p_{in};p_{in} (GeV/c);n #sigma_{#mu}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 1.f}, {100, -5.f, +5.f}}}},
      {"Track/hTOFNsigmaMu_after", "TOF n sigma #mu vs. p_{in};p_{in} (GeV/c);n #sigma_{#mu}^{TOF}", {HistType::kTH2F, {{1000, 0.f, 1.f}, {100, -5.f, +5.f}}}},
      {"Track/hTPCNsigmaPi_after", "TPC n sigma #pi vs. p_{in};p_{in} (GeV/c);n #sigma_{#pi}^{TPC}", {HistType::kTH2F, {{1000, 0.f, 1.f}, {100, -5.f, +5.f}}}},
      {"Track/hTOFNsigmaPi_after", "TOF n sigma #pi vs. p_{in};p_{in} (GeV/c);n #sigma_{#pi}^{TOF}", {HistType::kTH2F, {{1000, 0.f, 1.f}, {100, -5.f, +5.f}}}},
      {"Pair/hMmumuPtmumu", "ULS m_{#mu#mu} vs. p_{T,#mu#mu};m_{#mu#mu} (GeV/c^{2});p_{T,#mu#mu} (GeV/c)", {HistType::kTH2F, {{100, 0.2f, 1.2f}, {120, 0.f, 1.2f}}}},
      // for MC primary muon
      {"MC/Primary/hNclsITS", "Ncls ITS;N_{cls}^{ITS}", {HistType::kTH1F, {{8, -0.5f, +7.5f}}}},
      {"MC/Primary/hChi2ITS", "chi2 ITS;#chi^{2}/N_{cls}^{ITS}", {HistType::kTH1F, {{100, 0.f, 10.f}}}},
      {"MC/Primary/hNclsTPC", "Ncls TPC;N_{cls}^{TPC}", {HistType::kTH1F, {{161, -0.5f, +160.5f}}}},
      {"MC/Primary/hNcrTPC", "Ncr TPC;N_{cr}^{TPC}", {HistType::kTH1F, {{161, -0.5f, +160.5f}}}},
      {"MC/Primary/hNfTPC", "Nfindable TPC;N_{f}^{TPC}", {HistType::kTH1F, {{161, -0.5f, +160.5f}}}},
      {"MC/Primary/hChi2TPC", "chi2 TPC;#chi^{2}/N_{cls}^{TPC}", {HistType::kTH1F, {{100, 0.f, 10.f}}}},
      {"MC/Primary/hDCAxy_NclsTPC", "DCA xy vs. N_{cls}^{TPC};DCA_{xy} (cm);N_{cls}^{TPC}", {HistType::kTH2F, {{200, -1.f, +1.f}, {161, -0.5f, +160.5f}}}},
      {"MC/Primary/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", {HistType::kTH2F, {{200, -1.f, +1.f}, {200, -1.f, +1.f}}}},
      {"MC/Primary/hDCA3D", "DCA;DCA_{3D} (cm)", {HistType::kTH1F, {{100, 0.f, +1.f}}}},
      {"MC/Primary/hDCAxy_resolution", "DCA xy resolution;DCA_{xy} resolution (cm)", {HistType::kTH1F, {{100, 0.f, +0.1f}}}},
      {"MC/Primary/hDCAz_resolution", "DCA z resolution;DCA_{z} resolution (cm)", {HistType::kTH1F, {{100, 0.f, +0.1f}}}},
      {"MC/Primary/hDCAxyz_sigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", {HistType::kTH2F, {{200, -10.f, +10.f}, {200, -10.f, +10.f}}}},
      {"MC/Primary/hDCA3D_sigma", "DCA;DCA_{3D} (#sigma)", {HistType::kTH1F, {{100, 0.f, +10.f}}}},
      // for MC seconday muon
      {"MC/Secondary/hNclsITS", "Ncls ITS;N_{cls}^{ITS}", {HistType::kTH1F, {{8, -0.5f, +7.5f}}}},
      {"MC/Secondary/hChi2ITS", "chi2 ITS;#chi^{2}/N_{cls}^{ITS}", {HistType::kTH1F, {{100, 0.f, 10.f}}}},
      {"MC/Secondary/hNclsTPC", "Ncls TPC;N_{cls}^{TPC}", {HistType::kTH1F, {{161, -0.5f, +160.5f}}}},
      {"MC/Secondary/hNcrTPC", "Ncr TPC;N_{cr}^{TPC}", {HistType::kTH1F, {{161, -0.5f, +160.5f}}}},
      {"MC/Secondary/hNfTPC", "Nfindable TPC;N_{f}^{TPC}", {HistType::kTH1F, {{161, -0.5f, +160.5f}}}},
      {"MC/Secondary/hChi2TPC", "chi2 TPC;#chi^{2}/N_{cls}^{TPC}", {HistType::kTH1F, {{100, 0.f, 10.f}}}},
      {"MC/Secondary/hDCAxy_NclsTPC", "DCA xy vs. N_{cls}^{TPC};DCA_{xy} (cm);N_{cls}^{TPC}", {HistType::kTH2F, {{200, -1.f, +1.f}, {161, -0.5f, +160.5f}}}},
      {"MC/Secondary/hDCAxyz", "DCA xy vs. z;DCA_{xy} (cm);DCA_{z} (cm)", {HistType::kTH2F, {{200, -1.f, +1.f}, {200, -1.f, +1.f}}}},
      {"MC/Secondary/hDCA3D", "DCA;DCA_{3D} (cm)", {HistType::kTH1F, {{100, 0.f, +1.f}}}},
      {"MC/Secondary/hDCAxy_resolution", "DCA xy resolution;DCA_{xy} resolution (cm)", {HistType::kTH1F, {{100, 0.f, +0.1f}}}},
      {"MC/Secondary/hDCAz_resolution", "DCA z resolution;DCA_{z} resolution (cm)", {HistType::kTH1F, {{100, 0.f, +0.1f}}}},
      {"MC/Secondary/hDCAxyz_sigma", "DCA xy vs. z;DCA_{xy} (#sigma);DCA_{z} (#sigma)", {HistType::kTH2F, {{200, -10.f, +10.f}, {200, -10.f, +10.f}}}},
      {"MC/Secondary/hDCA3D_sigma", "DCA;DCA_{3D} (#sigma)", {HistType::kTH1F, {{100, 0.f, +10.f}}}},
    },
  };

  std::pair<int8_t, std::set<uint8_t>> itsRequirement = {1, {0, 1, 2}}; // any hits on 3 ITS ib layers.

  void init(InitContext const&) {}

  template <bool isMC, typename TTrack>
  bool checkTrack(TTrack const& track)
  {
    if constexpr (isMC) {
      if (!track.has_mcParticle()) {
        return false;
      }
    }

    if (!track.hasITS() || !track.hasTPC()) {
      return false;
    }

    // if(abs(track.dcaXY()) > dca_xy_max || abs(track.dcaZ()) > dca_z_max){
    //   return false;
    // }

    if (track.itsNCls() < minitsncls) {
      return false;
    }

    auto hits = std::count_if(itsRequirement.second.begin(), itsRequirement.second.end(), [&](auto&& requiredLayer) { return track.itsClusterMap() & (1 << requiredLayer); });
    if (hits < itsRequirement.first) {
      return false;
    }

    if (track.tpcNClsCrossedRows() < mincrossedrows) {
      return false;
    }

    if (track.tpcCrossedRowsOverFindableCls() < min_tpc_cr_findable_ratio) {
      return false;
    }

    return true;
  }

  template <typename TTrack>
  bool isMuon(TTrack const& track)
  {
    if (isInElectronBand(track)) { // reject electron first.
      return false;
    }

    if (track.hasTOF()) {
      return abs(track.tpcNSigmaMu()) < maxTPCNsigmaMu_highPin && abs(track.tofNSigmaMu()) < maxTOFNsigmaMu_highPin && track.tofNSigmaPi() < maxTOFNsigmaPi;
    } else if (track.tpcInnerParam() < maxPin_TPC) {
      return abs(track.tpcNSigmaMu()) < maxTPCNsigmaMu_lowPin && track.tpcNSigmaPi() < maxTPCNsigmaPi_lowPin;
    } else { // muon at high momentum cannot be identified without TOF.
      return false;
    }
  }

  template <typename TTrack>
  bool isInElectronBand(TTrack const& track)
  {
    return abs(track.tpcNSigmaEl()) < maxTPCNsigmaEl;
  }

  template <bool isMC, typename TTracks>
  void fillTrackHistogram(TTracks const& tracks)
  {
    for (auto& track : tracks) {
      if (!checkTrack<isMC>(track)) {
        continue;
      }

      fRegistry.fill(HIST("Track/hTPCdEdx_Pin_before"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/hTOFbeta_Pin_before"), track.tpcInnerParam(), track.beta());
      fRegistry.fill(HIST("Track/hTPCNsigmaMu_before"), track.tpcInnerParam(), track.tpcNSigmaMu());
      fRegistry.fill(HIST("Track/hTOFNsigmaMu_before"), track.tpcInnerParam(), track.tofNSigmaMu());
      fRegistry.fill(HIST("Track/hTPCNsigmaPi_before"), track.tpcInnerParam(), track.tpcNSigmaPi());
      fRegistry.fill(HIST("Track/hTOFNsigmaPi_before"), track.tpcInnerParam(), track.tofNSigmaPi());
      fRegistry.fill(HIST("Track/hTPCNsigmaEl_before"), track.tpcInnerParam(), track.tpcNSigmaEl());
      fRegistry.fill(HIST("Track/hTOFNsigmaEl_before"), track.tpcInnerParam(), track.tofNSigmaEl());

      if constexpr (isMC) {
        auto mctrack = track.template mcParticle_as<aod::McParticles>();
        if (abs(mctrack.pdgCode()) == 13 && mctrack.has_mothers()) {
          auto mp = mctrack.template mothers_first_as<aod::McParticles>();
          if (abs(mp.pdgCode()) == 221 || abs(mp.pdgCode()) == 223 || abs(mp.pdgCode()) == 333 || abs(mp.pdgCode()) == 113 || abs(mp.pdgCode()) == 331) {
            fRegistry.fill(HIST("MC/Primary/hNclsITS"), track.itsNCls());
            fRegistry.fill(HIST("MC/Primary/hChi2ITS"), track.itsChi2NCl());
            fRegistry.fill(HIST("MC/Primary/hNcrTPC"), track.tpcNClsCrossedRows());
            fRegistry.fill(HIST("MC/Primary/hNclsTPC"), track.tpcNClsFound());
            fRegistry.fill(HIST("MC/Primary/hNfTPC"), track.tpcNClsFindable());
            fRegistry.fill(HIST("MC/Primary/hChi2TPC"), track.tpcChi2NCl());
            fRegistry.fill(HIST("MC/Primary/hDCAxy_NclsTPC"), track.dcaXY(), track.tpcNClsFound());
            fRegistry.fill(HIST("MC/Primary/hDCAxyz"), track.dcaXY(), track.dcaZ());
            fRegistry.fill(HIST("MC/Primary/hDCA3D"), sqrt(pow(track.dcaXY(), 2) + pow(track.dcaZ(), 2)));
            fRegistry.fill(HIST("MC/Primary/hDCAxy_resolution"), sqrt(track.cYY()));
            fRegistry.fill(HIST("MC/Primary/hDCAz_resolution"), sqrt(track.cZZ()));
            fRegistry.fill(HIST("MC/Primary/hDCAxyz_sigma"), track.dcaXY() / sqrt(track.cYY()), track.dcaZ() / sqrt(track.cZZ()));
            fRegistry.fill(HIST("MC/Primary/hDCA3D_sigma"), sqrt(pow(track.dcaXY() / sqrt(track.cYY()), 2) + pow(track.dcaZ() / sqrt(track.cZZ()), 2)));
          } else if (abs(mp.pdgCode()) == 211 || abs(mp.pdgCode()) == 321) {
            fRegistry.fill(HIST("MC/Secondary/hNclsITS"), track.itsNCls());
            fRegistry.fill(HIST("MC/Secondary/hChi2ITS"), track.itsChi2NCl());
            fRegistry.fill(HIST("MC/Secondary/hNcrTPC"), track.tpcNClsCrossedRows());
            fRegistry.fill(HIST("MC/Secondary/hNclsTPC"), track.tpcNClsFound());
            fRegistry.fill(HIST("MC/Secondary/hNfTPC"), track.tpcNClsFindable());
            fRegistry.fill(HIST("MC/Secondary/hChi2TPC"), track.tpcChi2NCl());
            fRegistry.fill(HIST("MC/Secondary/hDCAxy_NclsTPC"), track.dcaXY(), track.tpcNClsFound());
            fRegistry.fill(HIST("MC/Secondary/hDCAxyz"), track.dcaXY(), track.dcaZ());
            fRegistry.fill(HIST("MC/Secondary/hDCA3D"), sqrt(pow(track.dcaXY(), 2) + pow(track.dcaZ(), 2)));
            fRegistry.fill(HIST("MC/Secondary/hDCAxy_resolution"), sqrt(track.cYY()));
            fRegistry.fill(HIST("MC/Secondary/hDCAz_resolution"), sqrt(track.cZZ()));
            fRegistry.fill(HIST("MC/Secondary/hDCAxyz_sigma"), track.dcaXY() / sqrt(track.cYY()), track.dcaZ() / sqrt(track.cZZ()));
            fRegistry.fill(HIST("MC/Secondary/hDCA3D_sigma"), sqrt(pow(track.dcaXY() / sqrt(track.cYY()), 2) + pow(track.dcaZ() / sqrt(track.cZZ()), 2)));
          }
        }
      }
    }
  }

  template <typename TTrack>
  void fillTrackTable(TTrack const& track)
  {
    // bool is_stored = std::find(stored_trackIds.begin(), stored_trackIds.end(), track.globalIndex()) != stored_trackIds.end();
    if (std::find(stored_trackIds.begin(), stored_trackIds.end(), track.globalIndex()) == stored_trackIds.end()) {
      emprimarymuons(track.collisionId(), track.globalIndex(), track.sign(),
                     track.pt(), track.eta(), track.phi(), track.dcaXY(), track.dcaZ(),
                     track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                     track.tpcChi2NCl(), track.tpcInnerParam(),
                     track.tpcSignal(), track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                     track.beta(), track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr(),
                     track.itsClusterMap(), track.itsChi2NCl(), track.detectorMap(), track.signed1Pt(), track.cYY(), track.cZZ());
      fRegistry.fill(HIST("Track/hTPCdEdx_Pin_after"), track.tpcInnerParam(), track.tpcSignal());
      fRegistry.fill(HIST("Track/hTOFbeta_Pin_after"), track.tpcInnerParam(), track.beta());
      fRegistry.fill(HIST("Track/hTPCNsigmaMu_after"), track.tpcInnerParam(), track.tpcNSigmaMu());
      fRegistry.fill(HIST("Track/hTOFNsigmaMu_after"), track.tpcInnerParam(), track.tofNSigmaMu());
      fRegistry.fill(HIST("Track/hTPCNsigmaPi_after"), track.tpcInnerParam(), track.tpcNSigmaPi());
      fRegistry.fill(HIST("Track/hTOFNsigmaPi_after"), track.tpcInnerParam(), track.tofNSigmaPi());
      fRegistry.fill(HIST("Track/hTPCNsigmaEl_after"), track.tpcInnerParam(), track.tpcNSigmaEl());
      fRegistry.fill(HIST("Track/hTOFNsigmaEl_after"), track.tpcInnerParam(), track.tofNSigmaEl());
      stored_trackIds.emplace_back(track.globalIndex());
      muon_pfb(0);
    }
  }

  template <bool isMC, EM_EEPairType pairtype, typename TCollision, typename TTracks1, typename TTracks2>
  void fillPairInfo(TCollision const& collision, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    if constexpr (pairtype == EM_EEPairType::kULS) { // ULS
      for (auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(tracks1, tracks2))) {
        if (!checkTrack<isMC>(t1) || !checkTrack<isMC>(t2)) {
          continue;
        }
        if (!isMuon(t1) || !isMuon(t2)) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        fRegistry.fill(HIST("Pair/hMmumuPtmumu"), v12.M(), v12.Pt());

        if (v12.M() > maxMmumu) { // don't store
          continue;
        }
        fRegistry.fill(HIST("hNpairs"), static_cast<int>(pairtype));
        fillTrackTable(t1);
        fillTrackTable(t2);
      }      // end of pairing loop
    } else { // LS
      for (auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(tracks1, tracks2))) {
        if (!checkTrack<isMC>(t1) || !checkTrack<isMC>(t2)) {
          continue;
        }
        if (!isMuon(t1) || !isMuon(t2)) {
          continue;
        }

        ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), o2::constants::physics::MassMuon);
        ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;
        if (v12.M() > maxMmumu) { // don't store
          continue;
        }
        fRegistry.fill(HIST("hNpairs"), static_cast<int>(pairtype));
        fillTrackTable(t1);
        fillTrackTable(t2);
      } // end of pairing loop
    }
  }

  // ============================ FUNCTION DEFINITIONS ====================================================
  std::vector<uint64_t> stored_trackIds;

  Filter trackFilter = minpt < o2::aod::track::pt && o2::aod::track::pt < maxpt && nabs(o2::aod::track::eta) < maxeta && nabs(o2::aod::track::dcaXY) < dca_xy_max && nabs(o2::aod::track::dcaZ) < dca_z_max && o2::aod::track::tpcChi2NCl < maxchi2tpc && o2::aod::track::itsChi2NCl < maxchi2its && min_tpcdEdx < o2::aod::track::tpcSignal && o2::aod::track::tpcSignal < max_tpcdEdx && o2::aod::track::tpcInnerParam < maxPin;

  using MyFilteredTracks = soa::Filtered<MyTracks>;
  Partition<MyFilteredTracks> posTracks = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracks> negTracks = o2::aod::track::signed1Pt < 0.f;
  void processRec(aod::Collisions const& collisions, aod::BCsWithTimestamps const&, MyFilteredTracks const& tracks)
  {
    stored_trackIds.reserve(tracks.size());
    for (auto& collision : collisions) {
      auto posTracks_per_coll = posTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracks->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      fillTrackHistogram<false>(posTracks_per_coll);
      fillTrackHistogram<false>(negTracks_per_coll);

      fillPairInfo<false, EM_EEPairType::kULS>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
      if (storeLS) {
        fillPairInfo<false, EM_EEPairType::kLSpp>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
        fillPairInfo<false, EM_EEPairType::kLSnn>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
      }
    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processRec, "process reconstructed info only", true);

  using MyFilteredTracksMC = soa::Filtered<MyTracksMC>;
  Partition<MyFilteredTracksMC> posTracksMC = o2::aod::track::signed1Pt > 0.f;
  Partition<MyFilteredTracksMC> negTracksMC = o2::aod::track::signed1Pt < 0.f;
  void processMC(soa::Join<aod::McCollisionLabels, aod::Collisions> const& collisions, aod::McCollisions const&, aod::BCsWithTimestamps const&, MyFilteredTracksMC const& tracks, aod::McParticles const&)
  {
    stored_trackIds.reserve(tracks.size());
    for (auto& collision : collisions) {
      if (!collision.has_mcCollision()) {
        continue;
      }

      auto posTracks_per_coll = posTracksMC->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracks_per_coll = negTracksMC->sliceByCached(o2::aod::track::collisionId, collision.globalIndex(), cache);
      fillTrackHistogram<true>(posTracks_per_coll);
      fillTrackHistogram<true>(negTracks_per_coll);

      fillPairInfo<true, EM_EEPairType::kULS>(collision, posTracks_per_coll, negTracks_per_coll); // ULS
      if (storeLS) {
        fillPairInfo<true, EM_EEPairType::kLSpp>(collision, posTracks_per_coll, posTracks_per_coll); // LS++
        fillPairInfo<true, EM_EEPairType::kLSnn>(collision, negTracks_per_coll, negTracks_per_coll); // LS--
      }
    } // end of collision loop

    stored_trackIds.clear();
    stored_trackIds.shrink_to_fit();
  }
  PROCESS_SWITCH(skimmerPrimaryMuon, processMC, "process reconstructed and MC info ", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<skimmerPrimaryMuon>(cfgc, TaskName{"skimmer-primary-muon"})};
}
