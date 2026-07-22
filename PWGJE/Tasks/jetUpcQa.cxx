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

#include "PWGJE/Tasks/UPCJetQATables.h"
#include "PWGUD/Core/SGCutParHolder.h"
#include "PWGUD/Core/SGSelector.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>

#include <TPDGCode.h>

#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

#include <cmath>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Stand-alone QA task for validating UPC event tagging for charged-jet studies.
/// The tagging and QA developed here are intended to guide a future integration
/// of the required UPC information into the PWGJE derived-data workflow.
struct JetUpcQa {
  static constexpr int RequiredSingleGapA = 1;
  static constexpr int RequiredSingleGapC = 2;
  static constexpr int RequiredDoubleGap = 3;

  // Output tables
  Produces<aod::UpcJetEvents> upcJetEventTable;
  Produces<aod::UpcJets> upcJetTable;
  Produces<aod::UpcJetTracks> upcJetTrackTable;

  Produces<aod::UpcJetEventsMCD> upcJetEventTableMCD;
  Produces<aod::UpcJetsMCD> upcJetTableMCD;
  Produces<aod::UpcJetTracksMCD> upcJetTrackTableMCD;

  Produces<aod::UpcJetEventsMCP> upcJetEventTableMCP;
  Produces<aod::UpcJetsMCP> upcJetTableMCP;
  Produces<aod::UpcJetTracksMCP> upcJetTrackTableMCP;

  Produces<aod::UpcTracks> upcTrackTable;
  Produces<aod::UpcTracksMCD> upcTrackTableMCD;
  Produces<aod::UpcTracksMCP> upcTrackTableMCP;

  // Type aliases
  using BCsWithSel = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
  using CollisionsData = soa::Join<aod::Collisions, aod::EvSels,
                                   aod::CentFT0Ms, aod::FT0Mults, aod::FV0Mults>;
  using TracksData = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
  using CollisionsMC = soa::Join<aod::Collisions, aod::EvSels,
                                 aod::CentFT0Ms, aod::FT0Mults, aod::FV0Mults>;
  using TracksMC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;

  // Preslices
  // Preslice<TracksData> tracksPerCollision = aod::track::collisionId;
  // Preslice<TracksMC> tracksPerCollisionMC = aod::track::collisionId;

  // SGSelector buffers
  std::vector<float> amplitudesFV0;
  std::vector<float> amplitudesFT0A;
  std::vector<float> amplitudesFT0C;
  std::vector<float> amplitudesFDDA;
  std::vector<float> amplitudesFDDC;

  // Configurables
  Configurable<float> cfgFT0AMax{"cfgFT0AMax", 100.f, "FT0A amplitude threshold"};
  Configurable<float> cfgFT0CMax{"cfgFT0CMax", 50.f, "FT0C amplitude threshold"};
  Configurable<float> cfgFV0AMax{"cfgFV0AMax", 100.f, "FV0A amplitude threshold"};
  Configurable<float> cfgFDDAMax{"cfgFDDAMax", -1.f, "FDDA amplitude (-1 = ignore)"};
  Configurable<float> cfgFDDCMax{"cfgFDDCMax", -1.f, "FDDC amplitude (-1 = ignore)"};
  Configurable<float> cfgFITTime{"cfgFITTime", 4.f, "FIT time window (ns)"};
  Configurable<int> cfgMinNBCs{"cfgMinNBCs", 7, "Minimum number of BCs"};
  Configurable<int> cfgNDtcoll{"cfgNDtcoll", 1, "Time resolution multiplier for compatible BCs"};
  Configurable<int> cfgMinNTracks{"cfgMinNTracks", 0, "Min tracks"};
  Configurable<int> cfgMaxNTracks{"cfgMaxNTracks", 100, "Max tracks"};

  Configurable<float> cfgTrackPtMin{"cfgTrackPtMin", 0.15f, "Min track pT (GeV/c)"};
  Configurable<float> cfgTrackEtaMax{"cfgTrackEtaMax", 0.9f, "Max track |eta|"};
  Configurable<std::string> cfgTrackSelection{"cfgTrackSelection", "kinematicsOnly", "Reconstructed-track selection: kinematicsOnly, globalTracks, or udDefault"};
  Configurable<float> cfgTrackDCAZMax{"cfgTrackDCAZMax", 2.f, "Maximum |DCAz| for udDefault track selection (cm)"};
  Configurable<float> cfgTrackDCAXYMax{"cfgTrackDCAXYMax", 0.f, "Maximum |DCAxy| for udDefault (cm); 0 uses the pT-dependent UD cut"};
  Configurable<float> cfgTrackTPCChi2NClMax{"cfgTrackTPCChi2NClMax", 4.f, "Maximum TPC chi2 per cluster for udDefault"};
  Configurable<int> cfgTrackTPCNClsFindableMin{"cfgTrackTPCNClsFindableMin", 70, "Minimum findable TPC clusters for udDefault"};
  Configurable<float> cfgTrackITSChi2NClMax{"cfgTrackITSChi2NClMax", 36.f, "Maximum ITS chi2 per cluster for udDefault"};
  Configurable<float> cfgJetR{"cfgJetR", 0.4f, "Jet radius"};
  Configurable<float> cfgJetPtMin{"cfgJetPtMin", 5.0f, "Min jet pT (GeV/c)"};
  Configurable<float> cfgJetEtaMax{"cfgJetEtaMax", 0.5f, "Max jet |eta|"};
  Configurable<float> cfgZDCEnergyCut{"cfgZDCEnergyCut", 10.f, "ZDC common energy threshold for Xn tagging"};
  Configurable<float> cfgZDCTimeCut{"cfgZDCTimeCut", 2.f, "ZDC time window for Xn tagging (ns)"};
  Configurable<bool> cfgUseGapTagging{"cfgUseGapTagging", true, "Require SGSelector gap tag for event selection"};
  Configurable<bool> cfgUseZDCTagging{"cfgUseZDCTagging", false, "Require ZDC neutron tag for event selection"};
  Configurable<int> cfgRequiredGapSide{"cfgRequiredGapSide", 0, "Required gap side when gap tagging is enabled (0 = any selected gap, 1 = SingleGapA, 2 = SingleGapC, 3 = DoubleGap)"};
  Configurable<int> cfgRequiredNeutronClass{"cfgRequiredNeutronClass", 0, "Required neutron class when ZDC tagging is enabled (0 = any valid class, 1 = XnXn, 2 = Xn0n, 3 = 0nXn, 4 = 0n0n)"};

  struct ZDCTagInfo {
    float energyCommonZNA = -999.f;
    float energyCommonZNC = -999.f;
    float timeZNA = -999.f;
    float timeZNC = -999.f;
    int neutronClass = 0;
  };

  struct GapSelectionInfo {
    int gapSide = o2::aod::sgselector::NoGap;
    ZDCTagInfo zdc;
  };

  SGSelector sgSelector;
  SGCutParHolder sgCuts;

  HistogramRegistry registry{"registry"};
  enum class TrackSelectionMode {
    KinematicsOnly,
    GlobalTracks,
    UDDefault
  };
  TrackSelectionMode trackSelectionMode = TrackSelectionMode::KinematicsOnly;

  // Init
  void init(InitContext&)
  {
    if (cfgTrackSelection.value == "kinematicsOnly") {
      trackSelectionMode = TrackSelectionMode::KinematicsOnly;
    } else if (cfgTrackSelection.value == "globalTracks") {
      trackSelectionMode = TrackSelectionMode::GlobalTracks;
    } else if (cfgTrackSelection.value == "udDefault") {
      trackSelectionMode = TrackSelectionMode::UDDefault;
    } else {
      LOGP(fatal, "Unknown cfgTrackSelection '{}'. Choose kinematicsOnly, globalTracks, or udDefault.", cfgTrackSelection.value);
    }

    sgCuts.SetNDtcoll(cfgNDtcoll);
    sgCuts.SetMinNBCs(cfgMinNBCs);
    sgCuts.SetNTracks(cfgMinNTracks, cfgMaxNTracks);
    sgCuts.SetMaxFITtime(cfgFITTime);
    sgCuts.SetFITAmpLimits({cfgFV0AMax, cfgFT0AMax, cfgFT0CMax, cfgFDDAMax, cfgFDDCMax});

    // Event-level QA only
    registry.add("hGapSide", "Gap Side (Data);Gap Side;Counts", {HistType::kTH1F, {{6, -1.5, 4.5}}});
    registry.add("hGapSideMCD", "Gap Side (MCD);Gap Side;Counts", {HistType::kTH1F, {{6, -1.5, 4.5}}});

    registry.add("hPosZ", "Vertex Z (Data);z (cm);Counts", {HistType::kTH1F, {{100, -15., 15.}}});
    registry.add("hPosZMCD", "Vertex Z (MCD);z (cm);Counts", {HistType::kTH1F, {{100, -15., 15.}}});
    registry.add("hPosZMCP", "Vertex Z (MCP);z (cm);Counts", {HistType::kTH1F, {{100, -15., 15.}}});

    registry.add("hCentralityFT0M", "Centrality FT0M (Data);Cent;Counts", {HistType::kTH1F, {{100, 0., 100.}}});
    registry.add("hCentMCD", "Centrality FT0M (MCD);Cent;Counts", {HistType::kTH1F, {{100, 0., 100.}}});

    registry.add("hMCDFoundBC", "MCD has_foundBC;0=false 1=true;Counts", {HistType::kTH1F, {{2, -0.5, 1.5}}});
    registry.add("hMCPNCharged", "MCP N charged primaries;N;Counts", {HistType::kTH1F, {{200, -0.5, 199.5}}});

    registry.add("hZNAEnergy", "ZNA common energy (Data);E;Counts", {HistType::kTH1F, {{400, -10., 390.}}});
    registry.add("hZNCEnergy", "ZNC common energy (Data);E;Counts", {HistType::kTH1F, {{400, -10., 390.}}});
    registry.add("hZNATime", "ZNA time (Data);t (ns);Counts", {HistType::kTH1F, {{200, -20., 20.}}});
    registry.add("hZNCTime", "ZNC time (Data);t (ns);Counts", {HistType::kTH1F, {{200, -20., 20.}}});
    registry.add("hNeutronClass", "Neutron class (Data);class;Counts", {HistType::kTH1F, {{5, -0.5, 4.5}}});
    registry.add("hZNAEnergyMCD", "ZNA common energy (MCD);E;Counts", {HistType::kTH1F, {{400, -10., 390.}}});
    registry.add("hZNCEnergyMCD", "ZNC common energy (MCD);E;Counts", {HistType::kTH1F, {{400, -10., 390.}}});
    registry.add("hZNATimeMCD", "ZNA time (MCD);t (ns);Counts", {HistType::kTH1F, {{200, -20., 20.}}});
    registry.add("hZNCTimeMCD", "ZNC time (MCD);t (ns);Counts", {HistType::kTH1F, {{200, -20., 20.}}});
    registry.add("hNeutronClassMCD", "Neutron class (MCD);class;Counts", {HistType::kTH1F, {{5, -0.5, 4.5}}});
    registry.add("hSelectionBits", "Selection bits (Data);0 all 1 gap 2 zdc 3 accepted;Counts", {HistType::kTH1F, {{4, -0.5, 3.5}}});
    registry.add("hSelectionBitsMCD", "Selection bits (MCD);0 all 1 gap 2 zdc 3 accepted;Counts", {HistType::kTH1F, {{4, -0.5, 3.5}}});
  }

  template <typename TTrack>
  bool isSelectedTrack(TTrack const& track) const
  {
    if (std::abs(track.eta()) > cfgTrackEtaMax.value || track.pt() < cfgTrackPtMin.value) {
      return false;
    }
    if (trackSelectionMode == TrackSelectionMode::KinematicsOnly) {
      return true;
    }
    if (trackSelectionMode == TrackSelectionMode::GlobalTracks) {
      return track.isGlobalTrack();
    }
    if (!track.isPVContributor()) {
      return false;
    }
    if (std::abs(track.dcaZ()) > cfgTrackDCAZMax.value) {
      return false;
    }
    float maxDCAxy = cfgTrackDCAXYMax.value;
    if (maxDCAxy <= 0.f) {
      maxDCAxy = 0.0105f + 0.035f / std::pow(track.pt(), 1.1f);
    }
    if (std::abs(track.dcaXY()) > maxDCAxy) {
      return false;
    }
    if (track.tpcChi2NCl() > cfgTrackTPCChi2NClMax.value || track.tpcNClsFindable() < cfgTrackTPCNClsFindableMin.value) {
      return false;
    }
    if (track.itsChi2NCl() > cfgTrackITSChi2NClMax.value) {
      return false;
    }
    return true;
  }

  // Helper: event tagging and selection
  bool isSelectedGap(int gapSide) const
  {
    if (cfgRequiredGapSide.value == RequiredSingleGapA) {
      return gapSide == o2::aod::sgselector::SingleGapA;
    }
    if (cfgRequiredGapSide.value == RequiredSingleGapC) {
      return gapSide == o2::aod::sgselector::SingleGapC;
    }
    if (cfgRequiredGapSide.value == RequiredDoubleGap) {
      return gapSide == o2::aod::sgselector::DoubleGap;
    }
    return gapSide == o2::aod::sgselector::SingleGapA ||
           gapSide == o2::aod::sgselector::SingleGapC ||
           gapSide == o2::aod::sgselector::DoubleGap;
  }

  int getNeutronClass(ZDCTagInfo const& zdc)
  {
    if (zdc.energyCommonZNA < 0.f || zdc.energyCommonZNC < 0.f) {
      return 0;
    }

    bool aXn = zdc.energyCommonZNA > cfgZDCEnergyCut && std::abs(zdc.timeZNA) < cfgZDCTimeCut;
    bool cXn = zdc.energyCommonZNC > cfgZDCEnergyCut && std::abs(zdc.timeZNC) < cfgZDCTimeCut;
    bool a0n = zdc.energyCommonZNA <= cfgZDCEnergyCut;
    bool c0n = zdc.energyCommonZNC <= cfgZDCEnergyCut;

    if (aXn && cXn) {
      return 1; // XnXn
    }
    if (aXn && c0n) {
      return 2; // Xn0n
    }
    if (a0n && cXn) {
      return 3; // 0nXn
    }
    if (a0n && c0n) {
      return 4; // 0n0n
    }
    return 0;
  }

  bool isSelectedZDC(ZDCTagInfo const& zdc) const
  {
    if (cfgRequiredNeutronClass.value > 0) {
      return zdc.neutronClass == cfgRequiredNeutronClass.value;
    }
    return zdc.neutronClass > 0;
  }

  bool isAcceptedEvent(bool isGapTagged, bool isZDCTagged) const
  {
    return (!cfgUseGapTagging.value || isGapTagged) &&
           (!cfgUseZDCTagging.value || isZDCTagged);
  }

  void fillSelectionBits(bool isMCD, bool isGapTagged, bool isZDCTagged, bool isAccepted)
  {
    if (isMCD) {
      registry.fill(HIST("hSelectionBitsMCD"), 0);
      if (isGapTagged) {
        registry.fill(HIST("hSelectionBitsMCD"), 1);
      }
      if (isZDCTagged) {
        registry.fill(HIST("hSelectionBitsMCD"), 2);
      }
      if (isAccepted) {
        registry.fill(HIST("hSelectionBitsMCD"), 3);
      }
      return;
    }

    registry.fill(HIST("hSelectionBits"), 0);
    if (isGapTagged) {
      registry.fill(HIST("hSelectionBits"), 1);
    }
    if (isZDCTagged) {
      registry.fill(HIST("hSelectionBits"), 2);
    }
    if (isAccepted) {
      registry.fill(HIST("hSelectionBits"), 3);
    }
  }

  template <typename TBC>
  ZDCTagInfo getZDCTagInfo(TBC const& bc)
  {
    ZDCTagInfo zdcInfo;
    if (!bc.has_zdc()) {
      return zdcInfo;
    }

    auto const& zdc = bc.zdc();
    zdcInfo.energyCommonZNA = zdc.energyCommonZNA();
    zdcInfo.energyCommonZNC = zdc.energyCommonZNC();
    zdcInfo.timeZNA = zdc.timeZNA();
    zdcInfo.timeZNC = zdc.timeZNC();
    zdcInfo.neutronClass = getNeutronClass(zdcInfo);
    return zdcInfo;
  }

  template <typename TCollision, typename TBCs>
  GapSelectionInfo getGapSide(TCollision const& collision, TBCs const& bcs)
  {
    amplitudesFV0.clear();
    amplitudesFT0A.clear();
    amplitudesFT0C.clear();
    amplitudesFDDA.clear();
    amplitudesFDDC.clear();

    GapSelectionInfo selection;

    if (!collision.has_foundBC()) {
      return selection;
    }

    auto const bc = collision.template foundBC_as<BCsWithSel>();
    auto const bcRange = udhelpers::compatibleBCs(collision, sgCuts.NDtcoll(), bcs, sgCuts.minNBCs());

    auto const result = sgSelector.IsSelected(sgCuts, collision, bcRange, bc,
                                              &amplitudesFV0, &amplitudesFT0A, &amplitudesFT0C,
                                              &amplitudesFDDA, &amplitudesFDDC);

    selection.gapSide = result.value;
    if (result.bc) {
      selection.zdc = getZDCTagInfo(*result.bc);
    }

    if (!isSelectedGap(selection.gapSide)) {
      amplitudesFV0.clear();
      amplitudesFT0A.clear();
      amplitudesFT0C.clear();
      amplitudesFDDA.clear();
      amplitudesFDDC.clear();
    }

    return selection;
  }

  // Helper: jet finding for detector-level tracks
  template <typename TTracks, typename TJetTable, typename TTrkTable>
  void runJetFindingDet(TTracks const& tracks,
                        int evtIdx,
                        TJetTable& jetTable,
                        TTrkTable& trkTable)
  {
    std::vector<fastjet::PseudoJet> fjInput;
    fjInput.reserve(tracks.size());

    for (auto const& t : tracks) {
      if (!isSelectedTrack(t)) {
        continue;
      }

      fastjet::PseudoJet p;
      p.reset_PtYPhiM(t.pt(), t.eta(), t.phi(), o2::constants::physics::MassPiPlus);
      fjInput.push_back(p);
    }

    if (fjInput.empty()) {
      return;
    }

    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, cfgJetR);
    fastjet::AreaDefinition areaDef(fastjet::active_area_explicit_ghosts,
                                    fastjet::GhostedAreaSpec(cfgTrackEtaMax + cfgJetR));
    fastjet::ClusterSequenceArea cs(fjInput, jetDef, areaDef);
    auto jets = fastjet::sorted_by_pt(cs.inclusive_jets(cfgJetPtMin));

    for (auto const& jet : jets) {
      if (std::abs(jet.eta()) > cfgJetEtaMax) {
        continue;
      }

      auto constituents = jet.constituents();
      jetTable(evtIdx, jet.pt(), jet.eta(), jet.phi_std(),
               jet.area(), static_cast<int>(constituents.size()));
      int jetIdx = jetTable.lastIndex();

      for (auto const& c : constituents) {
        trkTable(jetIdx, c.pt(), c.eta(), c.phi_std());
      }
    }
  }

  // Helper: jet finding for particle-level MC
  template <typename TJetTable, typename TTrkTable>
  void runJetFindingMCP(aod::McParticles const& particles,
                        int evtIdx,
                        TJetTable& jetTable,
                        TTrkTable& trkTable)
  {
    std::vector<fastjet::PseudoJet> fjInput;
    int nCharged = 0;

    for (auto const& p : particles) {
      if (!p.isPhysicalPrimary()) {
        continue;
      }
      if (std::abs(p.eta()) > cfgTrackEtaMax) {
        continue;
      }
      if (p.pt() < cfgTrackPtMin) {
        continue;
      }

      auto pdg = std::abs(p.pdgCode());
      bool charged = (pdg == kPiPlus || pdg == kKPlus || pdg == kProton || pdg == kElectron || pdg == kMuonMinus);
      if (!charged) {
        continue;
      }

      fastjet::PseudoJet fpp;
      fpp.reset_PtYPhiM(p.pt(), p.eta(), p.phi(), o2::constants::physics::MassPiPlus);
      fjInput.push_back(fpp);
      ++nCharged;
    }

    registry.fill(HIST("hMCPNCharged"), nCharged);

    if (fjInput.empty()) {
      return;
    }

    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, cfgJetR);
    fastjet::AreaDefinition areaDef(fastjet::active_area_explicit_ghosts,
                                    fastjet::GhostedAreaSpec(cfgTrackEtaMax + cfgJetR));
    fastjet::ClusterSequenceArea cs(fjInput, jetDef, areaDef);
    auto jets = fastjet::sorted_by_pt(cs.inclusive_jets(cfgJetPtMin));

    for (auto const& jet : jets) {
      if (std::abs(jet.eta()) > cfgJetEtaMax) {
        continue;
      }

      auto constituents = jet.constituents();
      jetTable(evtIdx, jet.pt(), jet.eta(), jet.phi_std(),
               jet.area(), static_cast<int>(constituents.size()));
      int jetIdx = jetTable.lastIndex();

      for (auto const& c : constituents) {
        trkTable(jetIdx, c.pt(), c.eta(), c.phi_std());
      }
    }
  }

  template <typename TTracks, typename TTrackTable>
  void fillEventTracksDet(TTracks const& tracks, int evtIdx, TTrackTable& trackTable)
  {
    int nAll = 0;
    int nFilled = 0;
    for (auto const& t : tracks) {
      ++nAll;
      if (!isSelectedTrack(t))
        continue;
      trackTable(evtIdx, t.pt(), t.eta(), t.phi());
      ++nFilled;
    }
    LOGP(info, "fillEventTracksDet: evtIdx={} nAll={} nFilled={}", evtIdx, nAll, nFilled);
  }

  template <typename TTrackTable>
  void fillEventTracksMCP(aod::McParticles const& particles, int evtIdx, TTrackTable& trackTable)
  {
    for (auto const& p : particles) {
      if (!p.isPhysicalPrimary())
        continue;
      if (std::abs(p.eta()) > cfgTrackEtaMax)
        continue;
      if (p.pt() < cfgTrackPtMin)
        continue;

      auto pdg = std::abs(p.pdgCode());
      bool charged = (pdg == kPiPlus || pdg == kKPlus || pdg == kProton ||
                      pdg == kElectron || pdg == kMuonMinus);
      if (!charged)
        continue;

      trackTable(evtIdx, p.pt(), p.eta(), p.phi());
    }
  }

  // Process: Data
  void processData(CollisionsData::iterator const& collision,
                   BCsWithSel const& bcs,
                   aod::Zdcs const&,
                   aod::FT0s const&, aod::FV0As const&, aod::FDDs const&,
                   TracksData const& tracks)
  {
    auto gapInfo = getGapSide(collision, bcs);
    int gapSide = gapInfo.gapSide;
    bool isGapTagged = isSelectedGap(gapSide);
    bool isZDCTagged = isSelectedZDC(gapInfo.zdc);
    bool isAccepted = isAcceptedEvent(isGapTagged, isZDCTagged);
    registry.fill(HIST("hGapSide"), gapSide);
    registry.fill(HIST("hZNAEnergy"), gapInfo.zdc.energyCommonZNA);
    registry.fill(HIST("hZNCEnergy"), gapInfo.zdc.energyCommonZNC);
    registry.fill(HIST("hZNATime"), gapInfo.zdc.timeZNA);
    registry.fill(HIST("hZNCTime"), gapInfo.zdc.timeZNC);
    registry.fill(HIST("hNeutronClass"), gapInfo.zdc.neutronClass);
    fillSelectionBits(false, isGapTagged, isZDCTagged, isAccepted);

    if (!isAccepted) {
      return;
    }

    registry.fill(HIST("hPosZ"), collision.posZ());
    registry.fill(HIST("hCentralityFT0M"), collision.centFT0M());

    upcJetEventTable(gapSide,
                     collision.posZ(),
                     collision.centFT0M(),
                     gapInfo.zdc.energyCommonZNA,
                     gapInfo.zdc.energyCommonZNC,
                     gapInfo.zdc.timeZNA,
                     gapInfo.zdc.timeZNC,
                     gapInfo.zdc.neutronClass,
                     isGapTagged,
                     isZDCTagged);
    int evtIdx = upcJetEventTable.lastIndex();

    std::vector<fastjet::PseudoJet> fjInput;
    fjInput.reserve(256);

    for (auto const& t : tracks) {
      if (t.collisionId() != collision.globalIndex()) {
        continue;
      }
      if (!isSelectedTrack(t)) {
        continue;
      }

      // event-level track table
      upcTrackTable(evtIdx, t.pt(), t.eta(), t.phi());

      // jet finding input
      fastjet::PseudoJet p;
      p.reset_PtYPhiM(t.pt(), t.eta(), t.phi(), o2::constants::physics::MassPiPlus);
      fjInput.push_back(p);
    }

    if (fjInput.empty()) {
      return;
    }

    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, cfgJetR);
    fastjet::AreaDefinition areaDef(
      fastjet::active_area_explicit_ghosts,
      fastjet::GhostedAreaSpec(cfgTrackEtaMax + cfgJetR));

    fastjet::ClusterSequenceArea cs(fjInput, jetDef, areaDef);
    auto jets = fastjet::sorted_by_pt(cs.inclusive_jets(cfgJetPtMin));

    for (auto const& jet : jets) {
      if (std::abs(jet.eta()) > cfgJetEtaMax) {
        continue;
      }

      auto constituents = jet.constituents();

      upcJetTable(evtIdx,
                  jet.pt(),
                  jet.eta(),
                  jet.phi_std(),
                  jet.area(),
                  static_cast<int>(constituents.size()));

      int jetIdx = upcJetTable.lastIndex();

      for (auto const& c : constituents) {
        upcJetTrackTable(jetIdx, evtIdx, c.pt(), c.eta(), c.phi_std());
      }
    }
  }
  PROCESS_SWITCH(JetUpcQa, processData, "Process real data", false);
  // Process: MCD
  void processMCD(CollisionsMC::iterator const& collision,
                  BCsWithSel const& bcs,
                  aod::Zdcs const&,
                  aod::FT0s const&, aod::FV0As const&, aod::FDDs const&,
                  TracksMC const& tracks)
  {
    registry.fill(HIST("hMCDFoundBC"), collision.has_foundBC() ? 1 : 0);

    auto gapInfo = getGapSide(collision, bcs);
    int gapSide = gapInfo.gapSide;
    bool isGapTagged = isSelectedGap(gapSide);
    bool isZDCTagged = isSelectedZDC(gapInfo.zdc);
    bool isAccepted = isAcceptedEvent(isGapTagged, isZDCTagged);
    registry.fill(HIST("hGapSideMCD"), gapSide);
    registry.fill(HIST("hZNAEnergyMCD"), gapInfo.zdc.energyCommonZNA);
    registry.fill(HIST("hZNCEnergyMCD"), gapInfo.zdc.energyCommonZNC);
    registry.fill(HIST("hZNATimeMCD"), gapInfo.zdc.timeZNA);
    registry.fill(HIST("hZNCTimeMCD"), gapInfo.zdc.timeZNC);
    registry.fill(HIST("hNeutronClassMCD"), gapInfo.zdc.neutronClass);
    fillSelectionBits(true, isGapTagged, isZDCTagged, isAccepted);

    if (!isAccepted) {
      return;
    }

    registry.fill(HIST("hPosZMCD"), collision.posZ());
    registry.fill(HIST("hCentMCD"), collision.centFT0M());

    upcJetEventTableMCD(gapSide,
                        collision.posZ(),
                        collision.centFT0M(),
                        gapInfo.zdc.energyCommonZNA,
                        gapInfo.zdc.energyCommonZNC,
                        gapInfo.zdc.timeZNA,
                        gapInfo.zdc.timeZNC,
                        gapInfo.zdc.neutronClass,
                        isGapTagged,
                        isZDCTagged);
    int evtIdx = upcJetEventTableMCD.lastIndex();

    std::vector<fastjet::PseudoJet> fjInput;
    fjInput.reserve(256);

    for (auto const& t : tracks) {
      if (t.collisionId() != collision.globalIndex()) {
        continue;
      }
      if (!isSelectedTrack(t)) {
        continue;
      }

      upcTrackTableMCD(evtIdx, t.pt(), t.eta(), t.phi());

      fastjet::PseudoJet p;
      p.reset_PtYPhiM(t.pt(), t.eta(), t.phi(), o2::constants::physics::MassPiPlus);
      fjInput.push_back(p);
    }

    if (fjInput.empty()) {
      return;
    }

    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, cfgJetR);
    fastjet::AreaDefinition areaDef(
      fastjet::active_area_explicit_ghosts,
      fastjet::GhostedAreaSpec(cfgTrackEtaMax + cfgJetR));

    fastjet::ClusterSequenceArea cs(fjInput, jetDef, areaDef);
    auto jets = fastjet::sorted_by_pt(cs.inclusive_jets(cfgJetPtMin));

    for (auto const& jet : jets) {
      if (std::abs(jet.eta()) > cfgJetEtaMax) {
        continue;
      }

      auto constituents = jet.constituents();

      upcJetTableMCD(evtIdx,
                     jet.pt(),
                     jet.eta(),
                     jet.phi_std(),
                     jet.area(),
                     static_cast<int>(constituents.size()));

      int jetIdx = upcJetTableMCD.lastIndex();

      for (auto const& c : constituents) {
        upcJetTrackTableMCD(jetIdx, c.pt(), c.eta(), c.phi_std());
      }
    }
  }
  PROCESS_SWITCH(JetUpcQa, processMCD, "Process MC detector level", true);

  // Process: MCP
  void processMCP(aod::McCollision const& mcCollision,
                  aod::McParticles const& particles)
  {
    registry.fill(HIST("hPosZMCP"), mcCollision.posZ());

    upcJetEventTableMCP(mcCollision.posZ());
    int evtIdx = upcJetEventTableMCP.lastIndex();
    fillEventTracksMCP(particles, evtIdx, upcTrackTableMCP);
    runJetFindingMCP(particles, evtIdx, upcJetTableMCP, upcJetTrackTableMCP);
  }
  PROCESS_SWITCH(JetUpcQa, processMCP, "Process MC particle level", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<JetUpcQa>(cfgc)};
}
