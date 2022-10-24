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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"
#include "CommonConstants/LHCConstants.h"
#include "PWGUD/DataModel/UDTables.h"

using namespace o2::framework;
using namespace o2::framework::expressions;

struct UpcCandProducer {
  bool fDoMC{false};
  bool fDoSemiFwd{false};

  std::map<int32_t, int32_t> fNewPartIDs;

  Produces<o2::aod::UDMcCollisions> udMCCollisions;
  Produces<o2::aod::UDMcParticles> udMCParticles;

  Produces<o2::aod::UDFwdTracks> udFwdTracks;
  Produces<o2::aod::UDFwdTracksExtra> udFwdTracksExtra;
  Produces<o2::aod::UDMcFwdTrackLabels> udFwdTrackLabels;

  Produces<o2::aod::UDTracks> udTracks;
  Produces<o2::aod::UDTracksExtra> udTracksExtra;
  Produces<o2::aod::UDTracksDCA> udTracksDCA;
  Produces<o2::aod::UDTracksPID> udTracksPID;
  Produces<o2::aod::UDMcTrackLabels> udTrackLabels;

  Produces<o2::aod::UDCollisions> eventCandidates;
  Produces<o2::aod::UDCollisionsSels> eventCandidatesSels;

  struct FITInfo {
    float ampFT0A = -1;
    float ampFT0C = -1;
    float timeFT0A = -999.;
    float timeFT0C = -999.;
    uint8_t triggerMaskFT0 = 0;
    float ampFDDA = -1;
    float ampFDDC = -1;
    float timeFDDA = -999.;
    float timeFDDC = -999.;
    uint8_t triggerMaskFDD = 0;
    float ampFV0A = -1;
    float timeFV0A = -999.;
    uint8_t triggerMaskFV0A = 0;
    int32_t BBFT0Apf = 0;
    int32_t BGFT0Apf = 0;
    int32_t BBFT0Cpf = 0;
    int32_t BGFT0Cpf = 0;
    int32_t BBFV0Apf = 0;
    int32_t BGFV0Apf = 0;
    int32_t BBFDDApf = 0;
    int32_t BGFDDApf = 0;
    int32_t BBFDDCpf = 0;
    int32_t BGFDDCpf = 0;
  };

  // skimmer flags
  // choose a source of signal MC events
  Configurable<int> fSignalGenID{"signalGenID", 1, "Signal generator ID"};

  // cuts for forward tracks
  Configurable<int> fUseFwdCuts{"useFwdCuts", 1, "Use cuts for forward tracks"};
  Configurable<int> fTrackType{"trackType", 3, "Filter by Fwd. track type: -1 -> no filter, 0 -> MFT-MCH-MID, 2 -> MFT-MCH, 3 -> MCH-MID. See ForwardTrackTypeEnum"};
  // basic
  Configurable<float> fFwdPtLow{"fwdPtLow", 0.5, "Minimal Pt for forward tracks"};
  Configurable<float> fFwdPtHigh{"fwdPtHigh", 4., "Maximal Pt for forward tracks"};
  Configurable<float> fFwdEtaLow{"fwdEtaLow", -4.0, "Maximal Eta for forward tracks"};
  Configurable<float> fFwdEtaHigh{"fwdEtaHigh", -2.5, "Maximal Eta for forward tracks"};
  // quality
  Configurable<float> fMuonRAtAbsorberEndLow{"muonRAtAbsorberEndLow", 17.6, "Minimal muon R at absorber end"};
  Configurable<float> fMuonRAtAbsorberEndHigh{"muonRAtAbsorberEndHigh", 89.5, "Maximal muon R at absorber end"};
  Configurable<float> fMuonPDcaHighFirst{"fMuonPDcaHighFirst", 594.0, "Primary PDCA cut: Maximal value for R < 26.5"};
  Configurable<float> fMuonPDcaHighSecond{"fMuonPDcaHighSecond", 324.0, "Additional PDCA cut: Maximal value for R >= 26.5"};
  Configurable<float> fFwdChi2Low{"fwdChi2Low", 0.0, "Minimal Chi2 for forward tracks"};
  Configurable<float> fFwdChi2High{"fwdChi2High", 10000.0, "Maximal Chi2 for forward tracks"};

  // cuts for central-barrel tracks
  Configurable<int> fUseBarCuts{"useBarCuts", 1, "Use cuts for barrel tracks"};
  // basic
  Configurable<float> fBarPtLow{"barPtLow", 0., "Minimal Pt for barrel tracks"};
  Configurable<float> fBarPtHigh{"barPtHigh", 1000., "Maximal Pt for barrel tracks"};
  Configurable<float> fBarEtaLow{"barEtaLow", -0.9, "Maximal Eta for barrel tracks"};
  Configurable<float> fBarEtaHigh{"barEtaHigh", 0.9, "Maximal Eta for barrel tracks"};
  // quality: ITS
  Configurable<int> fITSNClusLow{"ITSNClusLow", 4, "Minimal number of ITS clusters"};
  Configurable<int> fITSNClusHigh{"ITSNClusHigh", 9, "Maximal number of ITS clusters"};
  Configurable<float> fITSChi2Low{"ITSChi2Low", 0., "Minimal Chi2 in ITS per cluster"};
  Configurable<float> fITSChi2High{"ITSChi2High", 5., "Maximal Chi2 in ITS per cluster"};
  // quality: TPC
  Configurable<int> fTPCNClusCRLow{"TPCNClusCRLow", 70, "Minimal number of TPC clusters (crossed rows)"};
  Configurable<int> fTPCNClusCRHigh{"TPCNClusCRHigh", 161, "Maximal number of TPC clusters (crossed rows)"};
  Configurable<float> fTPCChi2Low{"TPCChi2Low", 0., "Minimal Chi2 in TPC per cluster"};
  Configurable<float> fTPCChi2High{"TPCChi2High", 4., "Maximal Chi2 in TPC per cluster"};
  // quality: DCA
  Configurable<int> fCheckMaxDcaXY{"checkMaxDCACut", 1, "Apply cut on maximal DCA_xy"};
  Configurable<float> fDcaZLow{"dcaZLow", -3., "Minimal DCA_z for barrel tracks"};
  Configurable<float> fDcaZHigh{"dcaZHigh", 3., "Maximal DCA_z for barrel tracks"};
  // quality: TOF
  Configurable<int> fRequireTOF{"requireTOF", 0, "Require all tracks to have TOF matches"};
  // tracks from collisions: consider only tracks from collisions with N tracks less or equal than fMaxNContrib
  Configurable<int> fMaxNContrib{"maxNContrib", 2, "Central barrel: consider tracks from collisions with N contributors <= maxNContrib"};
  Configurable<int> fAmbigSwitch{"ambigSwitch", 0, "Central barrel: 0 -- loop over all tracks, 1 -- loop only over tracks with vertices"};

  // candidate producer flags
  Configurable<int> fCheckTPCPID{"checkTPCPID", 0, "Check TPC PID. Useful for central selection -- see `tpcPIDSwitch` option"};
  Configurable<int> fTPCPIDSwitch{"tpcPIDSwitch", 0, "PID switch: 0 -- two muons/pions, 1 -- two electrons, 2 -- electron + muon/pion"};
  Configurable<int> fFilterFT0{"filterFT0", 0, "Filter candidates by FT0 signals at the same BC"};
  Configurable<int> fNFwdProngs{"nFwdProngs", 2, "Matched forward tracks per candidate"};
  Configurable<int> fNBarProngs{"nBarProngs", 0, "Matched barrel tracks per candidate"};

  // QA histograms to check for tracks after cuts
  HistogramRegistry histRegistry{"HistRegistry", {}, OutputObjHandlingPolicy::AnalysisObject};

  enum FwdSels {
    kFwdSelAll = 0,
    kFwdSelPt,
    kFwdSelEta,
    kFwdSelRabs,
    kFwdSelpDCA,
    kFwdSelChi2,
    kNFwdSels
  };

  enum BarrelSels {
    kBarrelSelAll = 0,
    kBarrelSelHasTOF,
    kBarrelSelPt,
    kBarrelSelEta,
    kBarrelSelITSNCls,
    kBarrelSelITSChi2,
    kBarrelSelTPCNCls,
    kBarrelSelTPCChi2,
    kBarrelSelDCAXY,
    kBarrelSelDCAZ,
    kNBarrelSels
  };

  using BCsWithBcSels = o2::soa::Join<o2::aod::BCs, o2::aod::BcSels>;

  using ForwardTracks = o2::soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov>;

  using BarrelTracks = o2::soa::Join<o2::aod::Tracks, o2::aod::TracksExtra, o2::aod::TracksDCA,
                                     o2::aod::pidTPCFullEl, o2::aod::pidTPCFullMu, o2::aod::pidTPCFullPi, o2::aod::pidTPCFullKa, o2::aod::pidTPCFullPr,
                                     o2::aod::TOFSignal, o2::aod::pidTOFFullEl, o2::aod::pidTOFFullMu, o2::aod::pidTOFFullPi, o2::aod::pidTOFFullKa, o2::aod::pidTOFFullPr>;

  void init(InitContext&)
  {
    const AxisSpec axisSelFwd{kNFwdSels, 0., double(kNFwdSels), ""};
    histRegistry.add("MuonsSelCounter", "", kTH1F, {axisSelFwd});
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(kFwdSelAll + 1, "All");
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(kFwdSelPt + 1, "Pt");
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(kFwdSelEta + 1, "Eta");
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(kFwdSelRabs + 1, "Rabs");
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(kFwdSelpDCA + 1, "pDCA");
    histRegistry.get<TH1>(HIST("MuonsSelCounter"))->GetXaxis()->SetBinLabel(kFwdSelChi2 + 1, "Chi2");

    const AxisSpec axisSelBar{kNBarrelSels, 0., double(kNBarrelSels), ""};
    histRegistry.add("BarrelsSelCounter", "", kTH1F, {axisSelBar});
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(kBarrelSelAll + 1, "All");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(kBarrelSelHasTOF + 1, "HasTOF");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(kBarrelSelPt + 1, "Pt");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(kBarrelSelEta + 1, "Eta");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(kBarrelSelITSNCls + 1, "ITSNCls");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(kBarrelSelITSChi2 + 1, "ITSChi2");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(kBarrelSelTPCNCls + 1, "TPCNCls");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(kBarrelSelTPCChi2 + 1, "TPCChi2");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(kBarrelSelDCAXY + 1, "DCAXY");
    histRegistry.get<TH1>(HIST("BarrelsSelCounter"))->GetXaxis()->SetBinLabel(kBarrelSelDCAZ + 1, "DCAZ");
  }

  bool applyFwdCuts(const ForwardTracks::iterator& track)
  {
    histRegistry.fill(HIST("MuonsSelCounter"), kFwdSelAll, 1);
    // using any cuts at all?
    if (!fUseFwdCuts) {
      return true;
    }
    // check Pt cuts
    float pt = track.pt();
    bool checkPt = pt > fFwdPtLow && pt < fFwdPtHigh;
    if (checkPt) {
      histRegistry.fill(HIST("MuonsSelCounter"), kFwdSelPt, 1);
    }
    // check pseudorapidity cuts
    float eta = track.eta();
    bool checkEta = eta > fFwdEtaLow && eta < fFwdEtaHigh;
    if (checkEta) {
      histRegistry.fill(HIST("MuonsSelCounter"), kFwdSelEta, 1);
    }
    // check muon R
    float r = track.rAtAbsorberEnd();
    bool checkR = r > fMuonRAtAbsorberEndLow && r < fMuonRAtAbsorberEndHigh;
    if (checkR) {
      histRegistry.fill(HIST("MuonsSelCounter"), kFwdSelRabs, 1);
    }
    // check pDCA
    float pDCA = track.pDca();
    bool checkPDCA = r < 26.5 ? pDCA < fMuonPDcaHighFirst : pDCA < fMuonPDcaHighSecond;
    if (checkPDCA) {
      histRegistry.fill(HIST("MuonsSelCounter"), kFwdSelpDCA, 1);
    }
    // check Chi2
    float chi2 = track.chi2();
    bool checkChi2 = chi2 > fFwdChi2Low && chi2 < fFwdChi2High;
    if (checkChi2) {
      histRegistry.fill(HIST("MuonsSelCounter"), kFwdSelChi2, 1);
    }
    bool pass = checkPt && checkEta && checkR && checkPDCA && checkChi2;
    return pass;
  }

  bool applyBarCuts(const BarrelTracks::iterator& track)
  {
    histRegistry.fill(HIST("BarrelsSelCounter"), kBarrelSelAll, 1);
    // using any cuts at all?
    if (!fUseBarCuts) {
      return true;
    }
    // require TOF match
    bool hasTOF = true;
    if (fRequireTOF) {
      hasTOF = track.hasTOF();
      if (hasTOF) {
        histRegistry.fill(HIST("BarrelsSelCounter"), kBarrelSelHasTOF, 1);
      }
    }
    // check Pt cuts
    float pt = track.pt();
    bool checkPt = pt > fBarPtLow && pt < fBarPtHigh;
    if (checkPt) {
      histRegistry.fill(HIST("BarrelsSelCounter"), kBarrelSelPt, 1);
    }
    // check pseudorapidity cuts
    float eta = track.eta();
    bool checkEta = eta > fBarEtaLow && eta < fBarEtaHigh;
    if (checkEta) {
      histRegistry.fill(HIST("BarrelsSelCounter"), kBarrelSelEta, 1);
    }
    // check ITS cuts
    bool checkITSNClus = track.itsNCls() >= static_cast<uint8_t>(fITSNClusLow) && track.itsNCls() <= static_cast<uint8_t>(fITSNClusHigh);
    if (checkITSNClus) {
      histRegistry.fill(HIST("BarrelsSelCounter"), kBarrelSelITSNCls, 1);
    }
    bool checkChi2ITS = track.itsChi2NCl() > fITSChi2Low && track.itsChi2NCl() < fITSChi2High;
    if (checkChi2ITS) {
      histRegistry.fill(HIST("BarrelsSelCounter"), kBarrelSelITSChi2, 1);
    }
    // check TPC cuts
    int checkTPCNClus = track.tpcNClsCrossedRows() > static_cast<int16_t>(fTPCNClusCRLow) && track.tpcNClsCrossedRows() < static_cast<int16_t>(fTPCNClusCRHigh);
    if (checkTPCNClus) {
      histRegistry.fill(HIST("BarrelsSelCounter"), kBarrelSelTPCNCls, 1);
    }
    bool checkChi2TPC = track.tpcChi2NCl() > fTPCChi2Low && track.tpcChi2NCl() < fTPCChi2High;
    if (checkChi2TPC) {
      histRegistry.fill(HIST("BarrelsSelCounter"), kBarrelSelTPCChi2, 1);
    }
    // check DCA
    bool checkMaxDcaXY = true;
    if (fCheckMaxDcaXY) {
      float dca = track.dcaXY();
      float maxDCA = 0.0105f + 0.0350f / pow(pt, 1.1f);
      checkMaxDcaXY = dca < maxDCA;
      if (checkMaxDcaXY) {
        histRegistry.fill(HIST("BarrelsSelCounter"), kBarrelSelDCAXY, 1);
      }
    }
    bool checkDCAZ = track.dcaZ() > fDcaZLow && track.dcaZ() < fDcaZHigh;
    if (checkDCAZ) {
      histRegistry.fill(HIST("BarrelsSelCounter"), kBarrelSelDCAZ, 1);
    }
    bool pass = hasTOF && checkPt && checkEta && checkITSNClus && checkChi2ITS &&
                checkTPCNClus && checkChi2TPC && checkMaxDcaXY && checkDCAZ;
    return pass;
  }

  void skimMCInfo(o2::aod::McCollisions const& mcCollisions,
                  o2::aod::McParticles const& mcParticles,
                  BCsWithBcSels const& bcs)
  {
    std::vector<int32_t> newEventIDs(mcCollisions.size(), -1);

    int32_t newPartID = 0;
    int32_t newEventID = 0;
    int32_t nMCParticles = mcParticles.size();
    // loop over MC particles to select only the ones from signal events
    // and calculate new MC table IDs
    for (int32_t mcPartID = 0; mcPartID < nMCParticles; mcPartID++) {
      const auto& mcPart = mcParticles.iteratorAt(mcPartID);
      int32_t mcEventID = mcPart.mcCollisionId();
      const auto& mcEvent = mcCollisions.iteratorAt(mcEventID);
      bool isSignal = mcEvent.generatorsID() == fSignalGenID;
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

    // storing MC particles
    for (const auto& item : fNewPartIDs) {
      int32_t mcPartID = item.first;
      const auto& mcPart = mcParticles.iteratorAt(mcPartID);
      int32_t mcEventID = mcPart.mcCollisionId();
      int32_t newEventID = newEventIDs[mcEventID];
      // collecting new mother IDs
      const auto& motherIDs = mcPart.mothersIds();
      std::vector<int32_t> newMotherIDs;
      newMotherIDs.reserve(motherIDs.size());
      for (auto motherID : motherIDs) {
        if (motherID >= nMCParticles) {
          continue;
        }
        auto it = fNewPartIDs.find(motherID);
        if (it != fNewPartIDs.end()) {
          newMotherIDs.push_back(it->second);
        }
      }
      // collecting new daughter IDs
      const auto& daughterIDs = mcPart.daughtersIds();
      int32_t newDaughterIDs[2] = {-1, -1};
      if (daughterIDs.size() > 0) {
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
      udMCParticles(newEventID, mcPart.pdgCode(), mcPart.statusCode(), mcPart.flags(), newMotherIDs, newDaughterIDs,
                    mcPart.weight(), mcPart.px(), mcPart.py(), mcPart.pz(), mcPart.e());
    }

    // storing MC events
    for (int32_t i = 0; i < mcCollisions.size(); i++) {
      if (newEventIDs[i] == -1) {
        continue;
      }
      const auto& mcEvent = mcCollisions.iteratorAt(i);
      auto bc = mcEvent.bc_as<BCsWithBcSels>();
      udMCCollisions(bc.globalBC(), mcEvent.generatorsID(), mcEvent.posX(), mcEvent.posY(), mcEvent.posZ(),
                     mcEvent.t(), mcEvent.weight(), mcEvent.impactParameter());
    }

    newEventIDs.clear();
  }

  // trackSwitch: 0 -> forward, 1 -> barrel
  template <int32_t trackSwitch, typename TTracks>
  void filterTracks(TTracks* tracks,
                    o2::aod::Collisions const& collisions,
                    std::unordered_map<int32_t, int32_t>& ambTrIDs,
                    std::vector<int32_t>& filteredTrackIDs)
  {
    for (const auto& tr : *tracks) {
      // skip if track doesn't pass selection
      bool pass = false;
      if constexpr (trackSwitch == 0) {
        if (fTrackType != -1) {
          auto trType = tr.trackType();
          if (trType != fTrackType) {
            continue;
          }
        }
        pass = applyFwdCuts(tr);
      }
      if constexpr (trackSwitch == 1) {
        pass = applyBarCuts(tr);
      }
      if (pass) {
        auto ambIter = ambTrIDs.find(tr.globalIndex());
        if (ambIter == ambTrIDs.end()) {
          int32_t colId = tr.collisionId();
          const auto& col = collisions.iteratorAt(colId);
          if (col.numContrib() > fMaxNContrib) { // skip if multiplicity is too high
            continue;
          }
        } else if (fAmbigSwitch == 1) { // skip ambiguous tracks if needed
          continue;
        }
        filteredTrackIDs.push_back(tr.globalIndex());
      }
    }
  }

  template <typename TTrack, typename TAmbTracks>
  uint64_t getTrackBC(TTrack track,
                      TAmbTracks* ambTracks,
                      std::unordered_map<int32_t, int32_t>& ambTrIDs,
                      o2::aod::Collisions const& collisions,
                      BCsWithBcSels const& bcs)
  {
    uint64_t trackBC = 0;
    auto ambIter = ambTrIDs.find(track.globalIndex());
    if (ambIter == ambTrIDs.end()) {
      const auto& col = track.collision();
      trackBC = col.template bc_as<BCsWithBcSels>().globalBC();
    } else if (fAmbigSwitch != 1) {
      const auto& ambTr = ambTracks->iteratorAt(ambIter->second);
      const auto& bcSlice = ambTr.bc();
      auto first = bcSlice.begin();
      trackBC = first.globalBC();
    }
    int64_t tint = std::round(track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS);
    uint64_t bc = trackBC + tint;
    return bc;
  }

  template <typename TFwdTracks, typename TMcFwdTrackLabels>
  void fillFwdTracks(TFwdTracks* tracks,
                     std::vector<int32_t> const& trackIDs,
                     int32_t candID,
                     uint64_t bc,
                     TMcFwdTrackLabels* mcTrackLabels)
  {
    for (int32_t trackID : trackIDs) {
      const auto& track = tracks->iteratorAt(trackID);
      double trTime = track.trackTime() - std::round(track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS) * o2::constants::lhc::LHCBunchSpacingNS;
      udFwdTracks(candID, track.px(), track.py(), track.pz(), track.sign(), bc, trTime, track.trackTimeRes());
      udFwdTracksExtra(track.nClusters(), track.pDca(), track.rAtAbsorberEnd(), track.chi2(), track.chi2MatchMCHMID(),
                       track.mchBitMap(), track.midBitMap(), track.midBoards());
      // fill MC labels and masks if needed
      if (fDoMC) {
        const auto& label = mcTrackLabels->iteratorAt(trackID);
        uint16_t mcMask = label.mcMask();
        auto it = fNewPartIDs.find(label.mcParticleId());
        // signal tracks should always have an MC particle
        // background tracks have label == -1
        int32_t newPartID = it != fNewPartIDs.end() ? it->second : -1;
        udFwdTrackLabels(newPartID, mcMask);
      }
    }
  }

  template <typename TBarrelTracks, typename TMcTrackLabels>
  void fillBarrelTracks(TBarrelTracks* tracks,
                        std::vector<int32_t> const& trackIDs,
                        int32_t candID,
                        uint64_t bc,
                        TMcTrackLabels* mcTrackLabels)
  {
    for (int32_t trackID : trackIDs) {
      const auto& track = tracks->iteratorAt(trackID);
      double trTime = track.trackTime() - std::round(track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS) * o2::constants::lhc::LHCBunchSpacingNS;
      udTracks(candID, track.px(), track.py(), track.pz(), track.sign(), bc, trTime, track.trackTimeRes());
      udTracksExtra(track.itsClusterMap(), track.tpcNClsFindable(), track.tpcNClsFindableMinusFound(), track.tpcNClsFindableMinusCrossedRows(),
                    track.tpcNClsShared(), track.trdPattern(), track.itsChi2NCl(), track.tpcChi2NCl(), track.trdChi2(), track.tofChi2(),
                    track.tpcSignal(), track.tofSignal(), track.trdSignal(), track.length(), track.tofExpMom(), track.detectorMap());
      udTracksPID(track.tpcNSigmaEl(), track.tpcNSigmaMu(), track.tpcNSigmaPi(), track.tpcNSigmaKa(), track.tpcNSigmaPr(),
                  track.tofNSigmaEl(), track.tofNSigmaMu(), track.tofNSigmaPi(), track.tofNSigmaKa(), track.tofNSigmaPr());
      udTracksDCA(track.dcaZ(), track.dcaXY());
      // fill MC labels and masks if needed
      if (fDoMC) {
        const auto& label = mcTrackLabels->iteratorAt(trackID);
        uint16_t mcMask = label.mcMask();
        int32_t mcPartID = label.mcParticleId();
        auto it = fNewPartIDs.find(mcPartID);
        // signal tracks should always have an MC particle
        // background tracks have label == -1
        int32_t newPartID = it != fNewPartIDs.end() ? it->second : -1;
        udTrackLabels(newPartID, mcMask);
      }
    }
  }

  // filter tracks in the central barrel using TPC PID
  // return `true` if a candidate passes "number of tracks" requirement from `nBarProngs`, `false` otherwise
  // if candidate passes, `filteredTrackIDs` contain track IDs of passed tracks
  template <typename TBarrelTracks>
  bool checkTPCPID(TBarrelTracks* tracks, std::vector<int32_t> const& trackIDs, std::vector<int32_t>& filteredTrackIDs)
  {
    int32_t nIDs = trackIDs.size();
    std::vector<bool> pidFlagsEl(nIDs, false);
    std::vector<bool> pidFlagsPi(nIDs, false);
    int32_t nEl = 0;
    int32_t nPi = 0;

    for (int32_t itr = 0; itr < nIDs; itr++) {
      const auto& tr = tracks->iteratorAt(trackIDs[itr]);
      float tpcNSigmaEl = tr.tpcNSigmaEl();
      float tpcNSigmaPi = tr.tpcNSigmaPi();

      bool isEl = std::abs(tpcNSigmaEl) < 3.0f && std::abs(tpcNSigmaEl) < std::abs(tpcNSigmaPi);
      bool isPi = std::abs(tpcNSigmaPi) < 3.0f && std::abs(tpcNSigmaEl) > std::abs(tpcNSigmaPi);

      pidFlagsEl[itr] = isEl;
      if (isEl) {
        nEl++;
      }

      pidFlagsPi[itr] = isPi;
      if (isPi) {
        nPi++;
      }
    }

    bool pass = false;

    // two muons/pions
    if (fTPCPIDSwitch == 0 && nPi == 2 && nEl == 0) {
      for (int32_t itr = 0; itr < nIDs; itr++) {
        if (pidFlagsPi[itr]) {
          filteredTrackIDs.push_back(trackIDs[itr]);
        }
      }
      pass = true;
    }

    // two electrons
    if (fTPCPIDSwitch == 1 && nEl == 2) {
      for (int32_t itr = 0; itr < nIDs; itr++) {
        if (pidFlagsEl[itr]) {
          filteredTrackIDs.push_back(trackIDs[itr]);
        }
      }
      pass = true;
    }

    // electron + muon/pion
    if (fTPCPIDSwitch == 2 && nEl == 1 && nPi == 1) {
      for (int32_t itr = 0; itr < nIDs; itr++) {
        if (pidFlagsEl[itr]) {
          filteredTrackIDs.push_back(trackIDs[itr]);
        }
        if (pidFlagsPi[itr]) {
          filteredTrackIDs.push_back(trackIDs[itr]);
        }
      }
      pass = true;
    }

    return pass;
  }

  bool checkFT0(FITInfo& info, bool isCentral)
  {
    const uint64_t presBitNum = 16;
    const float ft0DummyTime = 32.767f;
    const float ft0DefaultTime = -999;
    const float fT0CBBlower = -1.0; // ns
    const float fT0CBBupper = 1.0;  // ns
    bool hasNoFT0 = false;
    if (isCentral) {
      bool isBB = TESTBIT(info.BBFT0Apf, presBitNum) || TESTBIT(info.BBFT0Cpf, presBitNum);
      bool isBG = TESTBIT(info.BGFT0Apf, presBitNum) || TESTBIT(info.BGFT0Cpf, presBitNum);
      hasNoFT0 = !isBB && !isBG;
    } else {
      bool checkA = std::abs(info.timeFT0A - ft0DummyTime) < 1e-3 || std::abs(info.timeFT0A - ft0DefaultTime) < 1e-3; // dummy or default time
      bool checkC = info.timeFT0C > fT0CBBlower && info.timeFT0C < fT0CBBupper;
      hasNoFT0 = checkA && checkC;
    }
    return hasNoFT0;
  }

  template <typename TBCs>
  void processFITInfo(TBCs const& bcs,
                      o2::aod::FT0s const& ft0s,
                      o2::aod::FDDs const& fdds,
                      o2::aod::FV0As const& fv0as,
                      std::map<uint64_t, std::pair<std::vector<int32_t>, std::vector<int32_t>>>& bcsMatchedTrIds,
                      std::unordered_map<uint64_t, FITInfo>& bcsWithFIT)
  {
    const uint64_t presBitNum = 16; // bit for present BC

    std::unordered_set<uint64_t> addedBCs;
    std::vector<std::pair<uint64_t, FITInfo>> vbcsWithFIT;

    // gather FIT information for each BC
    // even if there is no FIT info at all -> store dummy "info{}"
    for (const auto& bc : bcs) {
      int32_t ft0Id = bc.foundFT0Id();
      FITInfo info{};
      if (ft0Id != -1) {
        const auto& ft0 = ft0s.iteratorAt(ft0Id);
        info.timeFT0A = ft0.timeA();
        info.timeFT0C = ft0.timeC();
        const auto& ampsA = ft0.amplitudeA();
        const auto& ampsC = ft0.amplitudeC();
        info.ampFT0A = 0.;
        for (auto amp : ampsA) {
          info.ampFT0A += amp;
        }
        info.ampFT0C = 0.;
        for (auto amp : ampsC) {
          info.ampFT0C += amp;
        }
        info.triggerMaskFT0 = ft0.triggerMask();
        if (!bc.selection()[evsel::kNoBGT0A])
          SETBIT(info.BGFT0Apf, presBitNum);
        if (!bc.selection()[evsel::kNoBGT0C])
          SETBIT(info.BGFT0Cpf, presBitNum);
        if (bc.selection()[evsel::kIsBBT0A])
          SETBIT(info.BBFT0Apf, presBitNum);
        if (bc.selection()[evsel::kIsBBT0C])
          SETBIT(info.BBFT0Cpf, presBitNum);
      }
      int32_t fddId = bc.foundFDDId();
      if (fddId != -1) {
        const auto& fdd = fdds.iteratorAt(fddId);
        info.timeFDDA = fdd.timeA();
        info.timeFDDC = fdd.timeC();
        const auto& ampsA = fdd.chargeA();
        const auto& ampsC = fdd.chargeC();
        info.ampFDDA = 0.;
        for (auto amp : ampsA) {
          info.ampFDDA += amp;
        }
        info.ampFDDC = 0.;
        for (auto amp : ampsC) {
          info.ampFDDC += amp;
        }
        info.triggerMaskFDD = fdd.triggerMask();
        if (!bc.selection()[evsel::kNoBGFDA])
          SETBIT(info.BGFDDApf, presBitNum);
        if (!bc.selection()[evsel::kNoBGFDC])
          SETBIT(info.BGFDDCpf, presBitNum);
        if (bc.selection()[evsel::kIsBBFDA])
          SETBIT(info.BBFDDApf, presBitNum);
        if (bc.selection()[evsel::kIsBBFDC])
          SETBIT(info.BBFDDCpf, presBitNum);
      }
      int32_t fv0aId = bc.foundFV0Id();
      if (fv0aId != -1) {
        const auto& fv0a = fv0as.iteratorAt(fv0aId);
        info.timeFV0A = fv0a.time();
        const auto& amps = fv0a.amplitude();
        info.ampFV0A = 0.;
        for (auto amp : amps) {
          info.ampFV0A += amp;
        }
        info.triggerMaskFV0A = fv0a.triggerMask();
        if (!bc.selection()[evsel::kNoBGV0A])
          SETBIT(info.BGFV0Apf, presBitNum);
        if (bc.selection()[evsel::kIsBBV0A])
          SETBIT(info.BBFV0Apf, presBitNum);
      }
      vbcsWithFIT.push_back({bc.globalBC(), info});
      addedBCs.insert(bc.globalBC());
    }

    // add BCs from event candidates if not yet there
    for (const auto& item : bcsMatchedTrIds) {
      uint64_t bc = item.first;
      if (addedBCs.find(bc) == addedBCs.end()) {
        addedBCs.insert(bc);
        FITInfo info{};
        vbcsWithFIT.push_back({bc, info});
      }
    }

    addedBCs.clear();

    std::sort(vbcsWithFIT.begin(), vbcsWithFIT.end(), [](const auto& left, const auto& right) { return left.first < right.first; });

    int32_t nBCsWithFIT = vbcsWithFIT.size();
    for (int32_t ibc = 0; ibc < nBCsWithFIT; ++ibc) {
      uint64_t bc = vbcsWithFIT[ibc].first;
      auto& info = vbcsWithFIT[ibc].second;
      // scrolling back on the original vector
      for (int32_t jbc = ibc - 1; jbc >= 0; --jbc) {
        uint64_t pastbc = vbcsWithFIT[jbc].first;
        uint64_t deltaBC = bc - pastbc;
        if (deltaBC > 16) // range [1, 16]
          break;
        const auto& pastInfo = vbcsWithFIT[jbc].second;
        uint64_t pastBitNum = deltaBC - 1; // bits range in [0, 15]; bit == 16 -> present
        if (TESTBIT(pastInfo.BGFT0Apf, presBitNum))
          SETBIT(info.BGFT0Apf, pastBitNum);
        if (TESTBIT(pastInfo.BGFT0Cpf, presBitNum))
          SETBIT(info.BGFT0Cpf, pastBitNum);
        if (TESTBIT(pastInfo.BBFT0Apf, presBitNum))
          SETBIT(info.BBFT0Apf, pastBitNum);
        if (TESTBIT(pastInfo.BBFT0Cpf, presBitNum))
          SETBIT(info.BBFT0Cpf, pastBitNum);
        if (TESTBIT(pastInfo.BGFDDApf, presBitNum))
          SETBIT(info.BGFDDApf, pastBitNum);
        if (TESTBIT(pastInfo.BGFDDCpf, presBitNum))
          SETBIT(info.BGFDDCpf, pastBitNum);
        if (TESTBIT(pastInfo.BBFDDApf, presBitNum))
          SETBIT(info.BBFDDApf, pastBitNum);
        if (TESTBIT(pastInfo.BBFDDCpf, presBitNum))
          SETBIT(info.BBFDDCpf, pastBitNum);
        if (TESTBIT(pastInfo.BGFV0Apf, presBitNum))
          SETBIT(info.BGFV0Apf, pastBitNum);
        if (TESTBIT(pastInfo.BBFV0Apf, presBitNum))
          SETBIT(info.BBFV0Apf, pastBitNum);
      }
      // scrolling forward on the original vector
      for (int32_t jbc = ibc + 1; jbc <= nBCsWithFIT; ++jbc) {
        uint64_t futbc = vbcsWithFIT[jbc].first;
        uint64_t deltaBC = futbc - bc;
        if (deltaBC > 15) // [1, 15]
          break;
        const auto& futInfo = vbcsWithFIT[jbc].second;
        uint64_t futBitNum = deltaBC + 16; // bits range in [17, 31]
        if (TESTBIT(futInfo.BGFT0Apf, presBitNum))
          SETBIT(info.BGFT0Apf, futBitNum);
        if (TESTBIT(futInfo.BGFT0Cpf, presBitNum))
          SETBIT(info.BGFT0Cpf, futBitNum);
        if (TESTBIT(futInfo.BBFT0Apf, presBitNum))
          SETBIT(info.BBFT0Apf, futBitNum);
        if (TESTBIT(futInfo.BBFT0Cpf, presBitNum))
          SETBIT(info.BBFT0Cpf, futBitNum);
        if (TESTBIT(futInfo.BGFDDApf, presBitNum))
          SETBIT(info.BGFDDApf, futBitNum);
        if (TESTBIT(futInfo.BGFDDCpf, presBitNum))
          SETBIT(info.BGFDDCpf, futBitNum);
        if (TESTBIT(futInfo.BBFDDApf, presBitNum))
          SETBIT(info.BBFDDApf, futBitNum);
        if (TESTBIT(futInfo.BBFDDCpf, presBitNum))
          SETBIT(info.BBFDDCpf, futBitNum);
        if (TESTBIT(futInfo.BGFV0Apf, presBitNum))
          SETBIT(info.BGFV0Apf, futBitNum);
        if (TESTBIT(futInfo.BBFV0Apf, presBitNum))
          SETBIT(info.BBFV0Apf, futBitNum);
      }
      // add current pair to map
      bcsWithFIT.insert({bc, info});
    }

    vbcsWithFIT.clear();
  }

  template <typename TFwdTracks, typename TAmbiguousFwdTracks,
            typename TBarrelTracks, typename TAmbiguousTracks,
            typename TBCs,
            typename TMcFwdTrackLabels, typename TMcTrackLabels>
  void createCandidates(TFwdTracks* fwdTracks,
                        TBarrelTracks* barrelTracks,
                        TAmbiguousFwdTracks* ambFwdTracks,
                        TAmbiguousTracks* ambBarrelTracks,
                        TBCs const& bcs,
                        o2::aod::Collisions const& collisions,
                        o2::aod::FT0s const& ft0s,
                        o2::aod::FDDs const& fdds,
                        o2::aod::FV0As const& fv0as,
                        TMcFwdTrackLabels* mcFwdTrackLabels,
                        TMcTrackLabels* mcBarrelTrackLabels)
  {
    // pairs of global BCs and vectors of matched track IDs:
    // global BC <-> <vector of fwd. trackIDs, vector of barrel trackIDs>
    std::map<uint64_t, std::pair<std::vector<int32_t>, std::vector<int32_t>>> bcsMatchedTrIds;

    // track IDs after cuts
    std::vector<int32_t> filteredFwdTracks;
    std::vector<int32_t> filteredBarrelTracks;

    // trackID -> index in amb. track table
    std::unordered_map<int32_t, int32_t> ambFwdTrIds;
    std::unordered_map<int32_t, int32_t> ambBarrelTrIds;

    // forward matching
    if (fwdTracks != nullptr) {
      if (fAmbigSwitch != 1) {
        for (const auto& ambTr : *ambFwdTracks) {
          auto trId = ambTr.fwdtrackId();
          ambFwdTrIds[trId] = ambTr.globalIndex();
        }
      }
      filterTracks<0>(fwdTracks, collisions, ambFwdTrIds, filteredFwdTracks);
      for (int32_t fwdTrID : filteredFwdTracks) {
        const auto& fwdTr = fwdTracks->iteratorAt(fwdTrID);
        uint64_t bc = getTrackBC(fwdTr, ambFwdTracks, ambFwdTrIds, collisions, bcs);
        // search for BC:
        //  if found -> store track ID to vector of matched tracks
        //  else make a new vector of matched tracks and store track ID
        auto it = bcsMatchedTrIds.find(bc);
        if (it != bcsMatchedTrIds.end()) {
          it->second.first.emplace_back(fwdTr.globalIndex());
        } else {
          bcsMatchedTrIds[bc] = std::make_pair(std::vector<int32_t>(1, fwdTrID), std::vector<int32_t>());
        }
      }
    }

    // central barrel tracks
    if (barrelTracks != nullptr) {
      if (fAmbigSwitch != 1) {
        for (const auto& ambTr : *ambBarrelTracks) {
          auto trId = ambTr.trackId();
          ambBarrelTrIds[trId] = ambTr.globalIndex();
        }
      }
      filterTracks<1>(barrelTracks, collisions, ambBarrelTrIds, filteredBarrelTracks);
      for (int32_t barTrID : filteredBarrelTracks) {
        const auto& barTr = barrelTracks->iteratorAt(barTrID);
        uint64_t bc = getTrackBC(barTr, ambBarrelTracks, ambBarrelTrIds, collisions, bcs);
        // search for BC:
        //  if found -> store track ID to vector of matched tracks
        //  else make a new vector of matched tracks and store track ID
        auto it = bcsMatchedTrIds.find(bc);
        if (it != bcsMatchedTrIds.end()) {
          it->second.second.emplace_back(barTr.globalIndex());
        } else if (!fDoSemiFwd) { // tag central-barrel tracks to forward tracks in semiforward case
          bcsMatchedTrIds[bc] = std::make_pair(std::vector<int32_t>(), std::vector<int32_t>(1, barTrID));
        }
      }
    }

    std::unordered_map<uint64_t, FITInfo> bcsWithFIT;
    processFITInfo(bcs, ft0s, fdds, fv0as, bcsMatchedTrIds, bcsWithFIT);

    // todo: calculate position of UD collision?
    float dummyX = 0.;
    float dummyY = 0.;
    float dummyZ = 0.;

    int32_t runNumber = bcs.iteratorAt(0).runNumber();

    // storing n-prong matches
    int32_t candID = 0;
    for (const auto& item : bcsMatchedTrIds) {
      uint64_t bc = item.first;
      std::vector<int32_t> fwdTrackIDs = item.second.first;
      std::vector<int32_t> barrelTrackIDs = item.second.second;
      int32_t nFwdTracks = fwdTrackIDs.size();
      int32_t nBarTracks = barrelTrackIDs.size();
      // check number of tracks in a candidate
      bool checkForward = nFwdTracks == fNFwdProngs;
      bool checkCentral = nBarTracks == fNBarProngs;
      // check TPC PID if needed
      std::vector<int32_t> filteredTrackIDs;
      if (fCheckTPCPID) {
        checkCentral = checkTPCPID(barrelTracks, barrelTrackIDs, filteredTrackIDs);
        if (checkCentral) {
          barrelTrackIDs.swap(filteredTrackIDs);
          nBarTracks = barrelTrackIDs.size();
        }
      }
      if (!checkForward || !checkCentral) {
        continue;
      }
      int8_t netCharge = 0;
      float RgtrwTOF = 0.;
      uint16_t numContrib = nFwdTracks + nBarTracks;
      for (auto id : barrelTrackIDs) {
        const auto& tr = barrelTracks->iteratorAt(id);
        netCharge += tr.sign();
        if (tr.hasTOF()) {
          RgtrwTOF++;
        }
      }
      RgtrwTOF = nBarTracks != 0 ? RgtrwTOF / (float)nBarTracks : 0.;
      if (RgtrwTOF == 0 && fNBarProngs != 0) { // require at least 1 TOF track in central and semiforward cases
        continue;
      }
      for (auto id : fwdTrackIDs) {
        const auto& tr = fwdTracks->iteratorAt(id);
        netCharge += tr.sign();
      }
      // fetching FT0, FDD, FV0 information
      // if there is no relevant signal, dummy info will be used
      FITInfo fitInfo{};
      auto it = bcsWithFIT.find(bc);
      if (it != bcsWithFIT.end()) {
        fitInfo = it->second;
        if (fFilterFT0) { // check FT0 signal in the same BC
          bool hasNoFT0 = checkFT0(fitInfo, barrelTracks != nullptr && fwdTracks == nullptr);
          if (!hasNoFT0)
            continue;
        }
      }
      // store used tracks
      fillFwdTracks(fwdTracks, fwdTrackIDs, candID, bc, mcFwdTrackLabels);
      fillBarrelTracks(barrelTracks, barrelTrackIDs, candID, bc, mcBarrelTrackLabels);
      eventCandidates(bc, runNumber, dummyX, dummyY, dummyZ, numContrib, netCharge, RgtrwTOF);
      eventCandidatesSels(fitInfo.ampFT0A, fitInfo.ampFT0C, fitInfo.timeFT0A, fitInfo.timeFT0C, fitInfo.triggerMaskFT0,
                          fitInfo.ampFDDA, fitInfo.ampFDDC, fitInfo.timeFDDA, fitInfo.timeFDDC, fitInfo.triggerMaskFDD,
                          fitInfo.ampFV0A, fitInfo.timeFV0A, fitInfo.triggerMaskFV0A,
                          fitInfo.BBFT0Apf, fitInfo.BBFT0Cpf, fitInfo.BGFT0Apf, fitInfo.BGFT0Cpf,
                          fitInfo.BBFV0Apf, fitInfo.BGFV0Apf,
                          fitInfo.BBFDDApf, fitInfo.BBFDDCpf, fitInfo.BGFDDApf, fitInfo.BGFDDCpf);
      candID++;
    }

    ambFwdTrIds.clear();
    ambBarrelTrIds.clear();
    filteredFwdTracks.clear();
    filteredBarrelTracks.clear();
    bcsMatchedTrIds.clear();
  }

  // data processors
  // _________________________________________________

  // create candidates for forward region
  void processFwd(ForwardTracks const& fwdTracks,
                  o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                  BCsWithBcSels const& bcs,
                  o2::aod::Collisions const& collisions,
                  o2::aod::FT0s const& ft0s,
                  o2::aod::FDDs const& fdds,
                  o2::aod::FV0As const& fv0as)
  {
    fDoMC = false;
    fDoSemiFwd = false;
    createCandidates(&fwdTracks, (BarrelTracks*)nullptr,
                     &ambFwdTracks, (o2::aod::AmbiguousTracks*)nullptr,
                     bcs, collisions,
                     ft0s, fdds, fv0as,
                     (o2::aod::McFwdTrackLabels*)nullptr, (o2::aod::McTrackLabels*)nullptr);
  }

  // create candidates for semiforward region
  void processSemiFwd(ForwardTracks const& fwdTracks,
                      o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                      BarrelTracks const& barrelTracks,
                      o2::aod::AmbiguousTracks const& ambTracks,
                      BCsWithBcSels const& bcs,
                      o2::aod::Collisions const& collisions,
                      o2::aod::FT0s const& ft0s,
                      o2::aod::FDDs const& fdds,
                      o2::aod::FV0As const& fv0as)
  {
    fDoMC = false;
    fDoSemiFwd = true;
    createCandidates(&fwdTracks, &barrelTracks,
                     &ambFwdTracks, &ambTracks,
                     bcs, collisions,
                     ft0s, fdds, fv0as,
                     (o2::aod::McFwdTrackLabels*)nullptr, (o2::aod::McTrackLabels*)nullptr);
  }

  // create candidates for central region
  void processCentral(BarrelTracks const& barrelTracks,
                      o2::aod::AmbiguousTracks const& ambBarrelTracks,
                      BCsWithBcSels const& bcs,
                      o2::aod::Collisions const& collisions,
                      o2::aod::FT0s const& ft0s,
                      o2::aod::FDDs const& fdds,
                      o2::aod::FV0As const& fv0as)
  {
    fDoMC = false;
    fDoSemiFwd = false;
    createCandidates((ForwardTracks*)nullptr, &barrelTracks,
                     (o2::aod::AmbiguousFwdTracks*)nullptr, &ambBarrelTracks,
                     bcs, collisions,
                     ft0s, fdds, fv0as,
                     (o2::aod::McFwdTrackLabels*)nullptr, (o2::aod::McTrackLabels*)nullptr);
  }

  // MC processors
  // _________________________________________________

  // create candidates for forward region
  void processFwdMC(ForwardTracks const& fwdTracks,
                    o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                    BCsWithBcSels const& bcs,
                    o2::aod::Collisions const& collisions,
                    o2::aod::FT0s const& ft0s,
                    o2::aod::FDDs const& fdds,
                    o2::aod::FV0As const& fv0as,
                    o2::aod::McCollisions const& mcCollisions, o2::aod::McParticles const& mcParticles,
                    o2::aod::McFwdTrackLabels const& mcFwdTrackLabels)
  {
    fDoMC = true;
    fDoSemiFwd = false;
    skimMCInfo(mcCollisions, mcParticles, bcs);
    createCandidates(&fwdTracks, (BarrelTracks*)nullptr,
                     &ambFwdTracks, (o2::aod::AmbiguousTracks*)nullptr,
                     bcs, collisions,
                     ft0s, fdds, fv0as,
                     &mcFwdTrackLabels, (o2::aod::McTrackLabels*)nullptr);
    fNewPartIDs.clear();
  }

  // create candidates for semiforward region
  void processSemiFwdMC(ForwardTracks const& fwdTracks,
                        o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                        BarrelTracks const& barrelTracks,
                        o2::aod::AmbiguousTracks const& ambTracks,
                        BCsWithBcSels const& bcs,
                        o2::aod::Collisions const& collisions,
                        o2::aod::FT0s const& ft0s,
                        o2::aod::FDDs const& fdds,
                        o2::aod::FV0As const& fv0as,
                        o2::aod::McCollisions const& mcCollisions, o2::aod::McParticles const& mcParticles,
                        o2::aod::McFwdTrackLabels const& mcFwdTrackLabels, o2::aod::McTrackLabels const& mcBarrelTrackLabels)
  {
    fDoMC = true;
    fDoSemiFwd = true;
    skimMCInfo(mcCollisions, mcParticles, bcs);
    createCandidates(&fwdTracks, &barrelTracks,
                     &ambFwdTracks, &ambTracks,
                     bcs, collisions,
                     ft0s, fdds, fv0as,
                     &mcFwdTrackLabels, &mcBarrelTrackLabels);
    fNewPartIDs.clear();
  }

  // create candidates for central region
  void processCentralMC(BarrelTracks const& barrelTracks,
                        o2::aod::AmbiguousTracks const& ambBarrelTracks,
                        BCsWithBcSels const& bcs,
                        o2::aod::Collisions const& collisions,
                        o2::aod::FT0s const& ft0s,
                        o2::aod::FDDs const& fdds,
                        o2::aod::FV0As const& fv0as,
                        o2::aod::McCollisions const& mcCollisions, o2::aod::McParticles const& mcParticles,
                        o2::aod::McTrackLabels const& mcBarrelTrackLabels)
  {
    fDoMC = true;
    fDoSemiFwd = false;
    skimMCInfo(mcCollisions, mcParticles, bcs);
    createCandidates((ForwardTracks*)nullptr, &barrelTracks,
                     (o2::aod::AmbiguousFwdTracks*)nullptr, &ambBarrelTracks,
                     bcs, collisions,
                     ft0s, fdds, fv0as,
                     (o2::aod::McFwdTrackLabels*)nullptr, &mcBarrelTrackLabels);
    fNewPartIDs.clear();
  }

  PROCESS_SWITCH(UpcCandProducer, processFwd, "Produce candidates for forward rapidities", false);
  PROCESS_SWITCH(UpcCandProducer, processSemiFwd, "Produce candidates in semiforward region", false);
  PROCESS_SWITCH(UpcCandProducer, processCentral, "Produce candidates in central region", false);
  PROCESS_SWITCH(UpcCandProducer, processFwdMC, "Produce candidates for forward rapidities with MC information", false);
  PROCESS_SWITCH(UpcCandProducer, processSemiFwdMC, "Produce candidates in semiforward region with MC information", false);
  PROCESS_SWITCH(UpcCandProducer, processCentralMC, "Produce candidates in central region with MC information", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcCandProducer>(cfgc)};
}
