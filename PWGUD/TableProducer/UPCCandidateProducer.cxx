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

  // todo: get parameters from CCDB
  EventSelectionParams par;

  // helper struct
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
    // selection flags
    bool isBBFT0A = false;
    bool isBGFT0A = false;
    bool isBBFT0C = false;
    bool isBGFT0C = false;
    bool isBBFV0A = false;
    bool isBGFV0A = false;
    bool isBBFDDA = false;
    bool isBGFDDA = false;
    bool isBBFDDC = false;
    bool isBGFDDC = false;
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

    // use "default" parameters
    par.fV0ABBlower = -3.0;  // ns
    par.fV0ABBupper = +2.0;  // ns
    par.fV0ABGlower = 2.0;   // ns
    par.fV0ABGupper = 5.0;   // ns
    par.fFDABBlower = -3.0;  // ns
    par.fFDABBupper = +3.0;  // ns
    par.fFDABGlower = 10.0;  // ns
    par.fFDABGupper = 13.0;  // ns
    par.fFDCBBlower = -3.0;  // ns
    par.fFDCBBupper = +3.0;  // ns
    par.fFDCBGlower = -10.0; // ns
    par.fFDCBGupper = -3.0;  // ns
    par.fT0ABBlower = -1.0;  // ns
    par.fT0ABBupper = +1.0;  // ns
    par.fT0CBBlower = -1.0;  // ns
    par.fT0CBBupper = +1.0;  // ns
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

  template <typename TMcCollisions, typename TMcParticles, typename TBCs>
  void skimMCInfo(TMcCollisions const& mcCollisions,
                  TMcParticles const& mcParticles,
                  TBCs const& bcs)
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
      udMCCollisions(mcEvent.bc().globalBC(), mcEvent.generatorsID(), mcEvent.posX(), mcEvent.posY(), mcEvent.posZ(),
                     mcEvent.t(), mcEvent.weight(), mcEvent.impactParameter());
    }

    newEventIDs.clear();
  }

  // trackSwitch: 0 -> forward, 1 -> barrel
  template <int32_t trackSwitch, typename TTracks>
  void filterTracks(TTracks* tracks,
                    o2::aod::Collisions const& collisions,
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
        int32_t colId = tr.collisionId();
        if (colId >= 0) {
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
                      o2::aod::BCs const& bcs)
  {
    uint64_t trackBC = 0;
    int32_t colId = track.collisionId();
    if (colId >= 0) {
      const auto& col = track.collision();
      trackBC = col.bc().globalBC();
    } else if (fAmbigSwitch != 1) {
      const auto& ambTr = ambTracks->iteratorAt(ambTrIDs.at(track.globalIndex()));
      const auto& bcSlice = ambTr.bc();
      auto first = bcSlice.begin();
      trackBC = first.globalBC();
    }
    uint64_t tint = std::round(track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS);
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

  void processFITInfo(uint64_t bc,
                      std::map<uint64_t, int32_t>& bcsWithFT0,
                      std::map<uint64_t, int32_t>& bcsWithFDD,
                      std::map<uint64_t, int32_t>& bcsWithFV0A,
                      o2::aod::FT0s const& ft0s,
                      o2::aod::FDDs const& fdds,
                      o2::aod::FV0As const& fv0as,
                      FITInfo& fitInfo)
  {
    float timeV0ABG = -999.f;
    float timeT0ABG = -999.f;
    float timeT0CBG = -999.f;
    float timeFDABG = -999.f;
    float timeFDCBG = -999.f;

    // check FIT info in the same BC
    auto it = bcsWithFT0.find(bc);
    if (it != bcsWithFT0.end()) {
      const auto& ft0 = ft0s.iteratorAt(it->second);
      fitInfo.timeFT0A = ft0.timeA();
      fitInfo.timeFT0C = ft0.timeC();
      const auto& ampsA = ft0.amplitudeA();
      const auto& ampsC = ft0.amplitudeC();
      fitInfo.ampFT0A = 0.;
      for (auto amp : ampsA) {
        fitInfo.ampFT0A += amp;
      }
      fitInfo.ampFT0C = 0.;
      for (auto amp : ampsC) {
        fitInfo.ampFT0C += amp;
      }
      fitInfo.triggerMaskFT0 = ft0.triggerMask();
    }

    it = bcsWithFDD.find(bc);
    if (it != bcsWithFDD.end()) {
      const auto& fdd = fdds.iteratorAt(it->second);
      fitInfo.timeFDDA = fdd.timeA();
      fitInfo.timeFDDC = fdd.timeC();
      const auto& ampsA = fdd.chargeA();
      const auto& ampsC = fdd.chargeC();
      fitInfo.ampFDDA = 0.;
      for (auto amp : ampsA) {
        fitInfo.ampFDDA += amp;
      }
      fitInfo.ampFDDC = 0.;
      for (auto amp : ampsC) {
        fitInfo.ampFDDC += amp;
      }
      fitInfo.triggerMaskFDD = fdd.triggerMask();
    }

    it = bcsWithFV0A.find(bc);
    if (it != bcsWithFV0A.end()) {
      const auto& fv0a = fv0as.iteratorAt(it->second);
      fitInfo.timeFV0A = fv0a.time();
      const auto& amps = fv0a.amplitude();
      fitInfo.ampFV0A = 0.;
      for (auto amp : amps) {
        fitInfo.ampFV0A += amp;
      }
      fitInfo.triggerMaskFV0A = fv0a.triggerMask();
    }

    // check beam-gas
    it = bcsWithFT0.find(bc - 1);
    if (it != bcsWithFT0.end()) {
      const auto& ft0 = ft0s.iteratorAt(it->second);
      timeT0ABG = ft0.timeA();
      timeT0CBG = ft0.timeC();
    }

    it = bcsWithFDD.find(bc - 5);
    if (it != bcsWithFDD.end()) {
      const auto& ft0 = fdds.iteratorAt(it->second);
      timeFDABG = ft0.timeA();
      timeFDCBG = ft0.timeC();
    }

    it = bcsWithFV0A.find(bc - 1);
    if (it != bcsWithFV0A.end()) {
      const auto& fv0a = fv0as.iteratorAt(it->second);
      timeV0ABG = fv0a.time();
    }

    // beam-gas flags
    fitInfo.isBGFV0A = timeV0ABG > par.fV0ABGlower && timeV0ABG < par.fV0ABGupper;
    fitInfo.isBGFDDA = timeFDABG > par.fFDABGlower && timeFDABG < par.fFDABGupper;
    fitInfo.isBGFDDC = timeFDCBG > par.fFDCBGlower && timeFDCBG < par.fFDCBGupper;
    fitInfo.isBGFT0A = timeT0ABG > par.fT0ABGlower && timeT0ABG < par.fT0ABGupper;
    fitInfo.isBGFT0C = timeT0CBG > par.fT0CBGlower && timeT0CBG < par.fT0CBGupper;

    // beam-beam flags
    fitInfo.isBBFT0A = fitInfo.timeFT0A > par.fT0ABBlower && fitInfo.timeFT0A < par.fT0ABBupper;
    fitInfo.isBBFT0C = fitInfo.timeFT0C > par.fT0CBBlower && fitInfo.timeFT0C < par.fT0CBBupper;
    fitInfo.isBBFV0A = fitInfo.timeFV0A > par.fV0ABBlower && fitInfo.timeFV0A < par.fV0ABBupper;
    fitInfo.isBBFDDA = fitInfo.timeFDDA > par.fFDABBlower && fitInfo.timeFDDA < par.fFDABBupper;
    fitInfo.isBBFDDC = fitInfo.timeFDDC > par.fFDCBBlower && fitInfo.timeFDDC < par.fFDCBBupper;
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
    std::map<uint64_t, int32_t> BCsWithFT0;
    // collect BCs with FT0 signals
    for (const auto& ft0 : ft0s) {
      uint64_t bc = ft0.bc().globalBC();
      BCsWithFT0[bc] = ft0.globalIndex();
    }

    std::map<uint64_t, int32_t> BCsWithFDD;
    // collect BCs with FDD signals
    for (const auto& fdd : fdds) {
      uint64_t bc = fdd.bc().globalBC();
      BCsWithFDD[bc] = fdd.globalIndex();
    }

    std::map<uint64_t, int32_t> BCsWithFV0A;
    // collect BCs with FV0A signals
    for (const auto& fv0a : fv0as) {
      uint64_t bc = fv0a.bc().globalBC();
      BCsWithFV0A[bc] = fv0a.globalIndex();
    }

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
      filterTracks<0>(fwdTracks, collisions, filteredFwdTracks);
      if (fAmbigSwitch != 1) {
        for (const auto& ambTr : *ambFwdTracks) {
          auto trId = ambTr.fwdtrackId();
          ambFwdTrIds[trId] = ambTr.globalIndex();
        }
      }
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
      filterTracks<1>(barrelTracks, collisions, filteredBarrelTracks);
      if (fAmbigSwitch != 1) {
        for (const auto& ambTr : *ambBarrelTracks) {
          auto trId = ambTr.trackId();
          ambBarrelTrIds[trId] = ambTr.globalIndex();
        }
      }
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
      // store used tracks
      fillFwdTracks(fwdTracks, fwdTrackIDs, candID, bc, mcFwdTrackLabels);
      fillBarrelTracks(barrelTracks, barrelTrackIDs, candID, bc, mcBarrelTrackLabels);
      // fetching FT0, FDD, FV0 information
      // if there is no relevant signal, dummy info will be used
      FITInfo fitInfo;
      processFITInfo(bc, BCsWithFT0, BCsWithFDD, BCsWithFV0A, ft0s, fdds, fv0as, fitInfo);
      eventCandidates(bc, runNumber, dummyX, dummyY, dummyZ, numContrib, netCharge, RgtrwTOF);
      eventCandidatesSels(fitInfo.ampFT0A, fitInfo.ampFT0C, fitInfo.timeFT0A, fitInfo.timeFT0C, fitInfo.triggerMaskFT0,
                          fitInfo.ampFDDA, fitInfo.ampFDDC, fitInfo.timeFDDA, fitInfo.timeFDDC, fitInfo.triggerMaskFDD,
                          fitInfo.ampFV0A, fitInfo.timeFV0A, fitInfo.triggerMaskFV0A,
                          fitInfo.isBBFT0A, fitInfo.isBBFT0C, fitInfo.isBGFT0A, fitInfo.isBGFT0C,
                          fitInfo.isBBFV0A, fitInfo.isBGFV0A,
                          fitInfo.isBBFDDA, fitInfo.isBBFDDC, fitInfo.isBGFDDA, fitInfo.isBGFDDC);
      candID++;
    }

    ambFwdTrIds.clear();
    ambBarrelTrIds.clear();
    filteredFwdTracks.clear();
    filteredBarrelTracks.clear();
    bcsMatchedTrIds.clear();
    BCsWithFT0.clear();
    BCsWithFDD.clear();
    BCsWithFV0A.clear();
  }

  // data processors
  // _________________________________________________

  // create candidates for forward region
  void processFwd(ForwardTracks const& fwdTracks,
                  o2::aod::AmbiguousFwdTracks const& ambFwdTracks,
                  o2::aod::BCs const& bcs,
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
                      o2::aod::BCs const& bcs,
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
                      o2::aod::BCs const& bcs,
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
                    o2::aod::BCs const& bcs,
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
                        o2::aod::BCs const& bcs,
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
                        o2::aod::BCs const& bcs,
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
