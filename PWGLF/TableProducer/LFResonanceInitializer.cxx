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

/// \file LFResonanceInitializer.cxx
/// \brief Initializes variables for the resonance candidate producers
///
///
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>

#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/Utils/collisionCuts.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

/// Initializer for the resonance candidate producers
struct reso2initializer {
  float cXiMass = RecoDecay::getMassPDG(3312);

  Produces<aod::ResoCollisions> resoCollisions;
  Produces<aod::ResoTracks> reso2trks;
  Produces<aod::ResoV0s> reso2v0s;
  Produces<aod::ResoCascades> reso2cascades;
  Produces<aod::ResoMCTracks> reso2mctracks;
  Produces<aod::ResoMCParents> reso2mcparents;
  Produces<aod::ResoMCV0s> reso2mcv0s;
  Produces<aod::ResoMCCascades> reso2mccascades;

  // Configurables
  Configurable<bool> ConfIsRun3{"ConfIsRun3", false, "Running on Pilot beam"}; // Choose if running on converted data or pilot beam

  /// Event cuts
  o2::analysis::CollisonCuts colCuts;
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", true, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", kINT7, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", false, "Evt sel: check for offline selection"};

  // Pre-selection cuts
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> pidnSigmaPreSelectionCut{"pidnSigmaPreSelectionCut", 5.0f, "TPC and TOF PID cut (loose, improve performance)"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};

  /// DCA Selections for V0
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.05, "Track DCAr cut to PV Maximum"};
  Configurable<double> cMinV0PosDCArToPVcut{"cMinV0PosDCArToPVcut", 0.05f, "V0 Positive Track DCAr cut to PV Minimum"}; // Pre-selection
  Configurable<double> cMinV0NegDCArToPVcut{"cMinV0NegDCArToPVcut", 0.05f, "V0 Negative Track DCAr cut to PV Minimum"}; // Pre-selection
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};

  Configurable<double> cMinV0Radius{"cMinV0Radius", 5.0, "Minimum V0 radius from PV"};
  Configurable<double> cMaxV0Radius{"cMaxV0Radius", 200.0, "Maximum V0 radius from PV"};
  Configurable<double> cMinV0CosPA{"cMinV0CosPA", 0.995, "Minimum V0 CosPA to PV"};

  /// DCA Selections for Cascades
  Configurable<int> mincrossedrows_cascbach{"mincrossedrows_cascbach", 70, "min crossed rows for bachelor track from cascade"};
  Configurable<double> cMinCascBachDCArToPVcut{"cMinCascBachDCArToPVcut", 0.05f, "Cascade Bachelor Track DCAr cut to PV Minimum"};  // Pre-selection
  Configurable<double> cMaxCascBachDCArToPVcut{"cMaxCascBachDCArToPVcut", 999.0f, "Cascade Bachelor Track DCAr cut to PV Maximum"}; // Pre-selection
  Configurable<double> cMaxCascDCAV0Daughters{"cMaxCascDCAV0Daughters", 1.6, "Cascade DCA between V0 daughters Maximum"};
  Configurable<double> cMaxCascDCACascDaughters{"cMaxCascDCACascDaughters", 1.6, "Cascade DCA between Casc daughters Maximum"};
  Configurable<double> cMinCascCosPA{"cMinCascCosPA", 0.97, "Minimum Cascade CosPA to PV"};
  Configurable<double> cMinCascV0CosPA{"cMinCascV0CosPA", 0.97, "Minimum Cascade V0 CosPA to PV"};
  Configurable<double> cMaxCascV0Radius{"cMaxCascV0Radius", 200.0, "Maximum Cascade V0 radius from PV"};
  Configurable<double> cMinCascV0Radius{"cMinCascV0Radius", 0.0, "Minimum Cascade V0 radius from PV"};
  Configurable<double> cMaxCascRadius{"cMaxCascRadius", 200.0, "Maximum Cascade radius from PV"};
  Configurable<double> cMinCascRadius{"cMinCascRadius", 0.0, "Minimum Cascade radius from PV"};
  Configurable<double> cCascMassResol{"cCascMassResol", 999, "Cascade mass resolution"};

  HistogramRegistry qaRegistry{"QAHistos", {
                                             {"hGoodTrackIndices", "hGoodTrackIndices", {HistType::kTH1F, {{4, 0.0f, 4.0f}}}},
                                             {"hGoodMCTrackIndices", "hGoodMCTrackIndices", {HistType::kTH1F, {{4, 0.0f, 4.0f}}}},
                                             {"hGoodV0Indices", "hGoodV0Indices", {HistType::kTH1F, {{5, 0.0f, 5.0f}}}},
                                             {"hGoodMCV0Indices", "hGoodMCV0Indices", {HistType::kTH1F, {{5, 0.0f, 5.0f}}}},
                                             {"hGoodCascIndices", "hGoodCascIndices", {HistType::kTH1F, {{8, 0.0f, 8.0f}}}},
                                             {"hGoodMCCascIndices", "hGoodMCCascIndices", {HistType::kTH1F, {{8, 0.0f, 8.0f}}}},
                                           },
                               OutputObjHandlingPolicy::AnalysisObject};

  // Pre-filters for efficient process
  // Filter tofPIDFilter = aod::track::tofExpMom < 0.f || ((aod::track::tofExpMom > 0.f) && ((nabs(aod::pidtof::tofNSigmaPi) < pidnSigmaPreSelectionCut) || (nabs(aod::pidtof::tofNSigmaKa) < pidnSigmaPreSelectionCut) || (nabs(aod::pidtof::tofNSigmaPr) < pidnSigmaPreSelectionCut))); // TOF
  Filter tpcPIDFilter = nabs(aod::pidtpc::tpcNSigmaPi) < pidnSigmaPreSelectionCut || nabs(aod::pidtpc::tpcNSigmaKa) < pidnSigmaPreSelectionCut || nabs(aod::pidtpc::tpcNSigmaPr) < pidnSigmaPreSelectionCut; // TPC
  Filter trackFilter = nabs(aod::track::eta) < cfgCutEta;                                                                                                                                                    // Eta cut
  Filter trackCutFilter = requireGlobalTrackInFilter();                                                                                                                                                      // Global track cuts
  Filter collisionFilter = nabs(aod::collision::posZ) < ConfEvtZvtx;

  // MC Resonance parent filter
  Partition<aod::McParticles> selectedMCParticles = (nabs(aod::mcparticle::pdgCode) == 313)        // K*
                                                    || (nabs(aod::mcparticle::pdgCode) == 323)     // K*pm
                                                    || (nabs(aod::mcparticle::pdgCode) == 333)     // phi
                                                    || (nabs(aod::mcparticle::pdgCode) == 9010221) // f_0(980)
                                                    || (nabs(aod::mcparticle::pdgCode) == 10221)   // f_0(1370)
                                                    || (nabs(aod::mcparticle::pdgCode) == 9030221) // f_0(1500)
                                                    || (nabs(aod::mcparticle::pdgCode) == 10331)   // f_0(1710)
                                                    || (nabs(aod::mcparticle::pdgCode) == 113)     // rho(770)
                                                    || (nabs(aod::mcparticle::pdgCode) == 213)     // rho(770)pm
                                                    || (nabs(aod::mcparticle::pdgCode) == 3224)    // Sigma(1385)+
                                                    || (nabs(aod::mcparticle::pdgCode) == 3124)    // Sigma(1385)-
                                                    || (nabs(aod::mcparticle::pdgCode) == 3324)    // Xi(1530)0
                                                    || (nabs(aod::mcparticle::pdgCode) == 123314)  // Xi(1820)0
                                                    || (nabs(aod::mcparticle::pdgCode) == 123324); // Xi(1820)-0

  using ResoEvents = soa::Join<aod::Collisions, aod::EvSels, aod::Mults>;
  using ResoEventsMC = soa::Join<ResoEvents, aod::McCollisionLabels>;
  using ResoTracks = aod::Reso2TracksPIDExt;
  using ResoTracksMC = soa::Join<ResoTracks, aod::McTrackLabels>;
  using ResoV0s = aod::V0Datas;
  using ResoV0sMC = soa::Join<ResoV0s, aod::McV0Labels>;
  using ResoCascades = aod::CascDatas;
  using ResoCascadesMC = soa::Join<ResoCascades, aod::McCascLabels>;

  template <bool isMC, typename CollisionType, typename TrackType>
  bool IsTrackSelected(CollisionType const& collision, TrackType const& track)
  {
    // Track selection
    qaRegistry.fill(HIST("hGoodTrackIndices"), 0.5);
    // MC case can be handled here
    if constexpr (isMC) {
      // MC check
      qaRegistry.fill(HIST("hGoodMCTrackIndices"), 0.5);
    }
    qaRegistry.fill(HIST("hGoodTrackIndices"), 1.5);
    return true;
  }

  template <bool isMC, typename CollisionType, typename V0Type, typename TrackType>
  bool IsV0Selected(CollisionType const& collision, V0Type const& v0, TrackType const& track)
  {
    // V0 selection
    qaRegistry.fill(HIST("hGoodV0Indices"), 0.5);

    auto postrack = v0.template posTrack_as<TrackType>();
    auto negtrack = v0.template negTrack_as<TrackType>();

    if (postrack.tpcNClsCrossedRows() < mincrossedrows)
      return false;
    if (negtrack.tpcNClsCrossedRows() < mincrossedrows)
      return false;
    qaRegistry.fill(HIST("hGoodV0Indices"), 1.5);

    if (fabs(postrack.dcaXY()) < cMinV0PosDCArToPVcut)
      return false;
    if (fabs(negtrack.dcaXY()) < cMinV0NegDCArToPVcut)
      return false;
    qaRegistry.fill(HIST("hGoodV0Indices"), 2.5);

    if ((v0.v0radius() > cMaxV0Radius) || (v0.v0radius() < cMinV0Radius))
      return false;
    qaRegistry.fill(HIST("hGoodV0Indices"), 3.5);
    if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cMinV0CosPA)
      return false;
    qaRegistry.fill(HIST("hGoodV0Indices"), 4.5);

    // MC case can be handled here
    if constexpr (isMC) {
      // MC check
      qaRegistry.fill(HIST("hGoodMCV0Indices"), 0.5);
    }
    return true;
  }

  template <bool isMC, typename CollisionType, typename CascType, typename TrackType>
  bool IsCascSelected(CollisionType const& collision, CascType const& casc, TrackType const& track)
  {
    // V0 selection
    qaRegistry.fill(HIST("hGoodCascIndices"), 0.5);

    auto trackBach = casc.template bachelor_as<TrackType>();
    // auto trackPos = casc.template posTrack_as<TrackType>();
    // auto trackNeg = casc.template negTrack_as<TrackType>();

    // track cuts
    if (trackBach.tpcNClsCrossedRows() < mincrossedrows_cascbach)
      return false;
    qaRegistry.fill(HIST("hGoodCascIndices"), 1.5);

    if (fabs(trackBach.dcaXY()) < cMinCascBachDCArToPVcut)
      return false;
    if (fabs(trackBach.dcaXY()) > cMaxCascBachDCArToPVcut)
      return false;
    qaRegistry.fill(HIST("hGoodCascIndices"), 2.5);

    // DCA daugthers
    if (casc.dcaV0daughters() > cMaxCascDCAV0Daughters)
      return false;
    if (casc.dcacascdaughters() > cMaxCascDCACascDaughters)
      return false;
    qaRegistry.fill(HIST("hGoodCascIndices"), 3.5);

    // CPA cuts
    if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cMinCascCosPA)
      return false;
    if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cMinCascV0CosPA)
      return false;
    qaRegistry.fill(HIST("hGoodCascIndices"), 4.5);

    // V0 radius
    auto v0radius = casc.v0radius();
    if ((v0radius > cMaxCascV0Radius) || (v0radius < cMinCascV0Radius))
      return false;
    qaRegistry.fill(HIST("hGoodCascIndices"), 5.5);

    // Casc radius
    auto cascradius = casc.cascradius();
    if ((cascradius > cMaxCascRadius) || (cascradius < cMinCascRadius))
      return false;
    qaRegistry.fill(HIST("hGoodCascIndices"), 6.5);

    // Casc mass
    auto cascMass = casc.mXi();
    if (abs(cascMass - cXiMass) > cCascMassResol)
      return false;
    qaRegistry.fill(HIST("hGoodCascIndices"), 7.5);

    // MC case can be handled here
    if constexpr (isMC) {
      // MC check
      qaRegistry.fill(HIST("hGoodMCCascIndices"), 0.5);
    }
    return true;
  }

  // Filter for all tracks
  template <bool isMC, typename TrackType, typename CollisionType>
  void fillTracks(CollisionType const& collision, TrackType const& tracks)
  {
    // Loop over tracks
    for (auto& track : tracks) {
      if (!IsTrackSelected<isMC>(collision, track))
        continue;
      // Add PID selection criteria here
      uint8_t tpcPIDselections = 0;
      uint8_t tofPIDselections = 0;
      // TPC PID
      if (std::abs(track.tpcNSigmaPi()) < pidnSigmaPreSelectionCut)
        tpcPIDselections |= aod::resodaughter::PDGtype::kPion;
      if (std::abs(track.tpcNSigmaKa()) < pidnSigmaPreSelectionCut)
        tpcPIDselections |= aod::resodaughter::PDGtype::kKaon;
      if (std::abs(track.tpcNSigmaPr()) < pidnSigmaPreSelectionCut)
        tpcPIDselections |= aod::resodaughter::PDGtype::kProton;
      // TOF PID
      if (track.hasTOF()) {
        tofPIDselections |= aod::resodaughter::PDGtype::kHasTOF;
        if (std::abs(track.tofNSigmaPi()) < pidnSigmaPreSelectionCut)
          tofPIDselections |= aod::resodaughter::PDGtype::kPion;
        if (std::abs(track.tofNSigmaKa()) < pidnSigmaPreSelectionCut)
          tofPIDselections |= aod::resodaughter::PDGtype::kKaon;
        if (std::abs(track.tofNSigmaPr()) < pidnSigmaPreSelectionCut)
          tofPIDselections |= aod::resodaughter::PDGtype::kProton;
      }
      reso2trks(resoCollisions.lastIndex(),
                track.pt(),
                track.px(),
                track.py(),
                track.pz(),
                track.eta(),
                track.phi(),
                track.sign(),
                (uint8_t)track.tpcNClsCrossedRows(),
                track.dcaXY(),
                track.dcaZ(),
                track.x(),
                track.alpha(),
                tpcPIDselections,
                tofPIDselections,
                track.tpcNSigmaPi(),
                track.tpcNSigmaKa(),
                track.tpcNSigmaPr(),
                track.tofNSigmaPi(),
                track.tofNSigmaKa(),
                track.tofNSigmaPr());
      if constexpr (isMC) {
        fillMCTrack(track);
      }
    }
  }

  // Filter for all V0s
  template <bool isMC, typename CollisionType, typename V0Type, typename TrackType>
  void fillV0s(CollisionType const& collision, V0Type const& v0s, TrackType const& tracks)
  {
    int childIDs[2] = {0, 0}; // these IDs are necessary to keep track of the children
    for (auto& v0 : v0s) {
      if (!IsV0Selected<isMC>(collision, v0, tracks))
        continue;
      childIDs[0] = v0.posTrackId();
      childIDs[1] = v0.negTrackId();
      reso2v0s(resoCollisions.lastIndex(),
               v0.pt(),
               v0.px(),
               v0.py(),
               v0.pz(),
               v0.eta(),
               v0.phi(),
               childIDs,
               v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
               v0.dcaV0daughters(), v0.mLambda(), v0.mAntiLambda(),
               v0.v0radius(), v0.x(), v0.y(), v0.z());
      if constexpr (isMC) {
        fillMCV0(v0);
      }
    }
  }

  // Filter for all Cascades
  template <bool isMC, typename CollisionType, typename CascType, typename TrackType>
  void fillCascades(CollisionType const& collision, CascType const& cascades, TrackType const& tracks)
  {
    int childIDs[2] = {0, 0}; // these IDs are necessary to keep track of the children
    for (auto& casc : cascades) {
      if (!IsCascSelected<isMC>(collision, casc, tracks))
        continue;
      childIDs[0] = casc.v0Id();
      childIDs[1] = casc.bachelorId();
      reso2cascades(resoCollisions.lastIndex(),
                    casc.pt(),
                    casc.px(),
                    casc.py(),
                    casc.pz(),
                    casc.eta(),
                    casc.phi(),
                    childIDs,
                    casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
                    casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()),
                    casc.dcaV0daughters(), casc.dcacascdaughters(), casc.mXi(),
                    casc.v0radius(), casc.cascradius(), casc.x(), casc.y(), casc.z());
      if constexpr (isMC) {
        fillMCCascade(casc);
      }
    }
  }

  template <typename TrackType>
  void fillMCTrack(TrackType const& track)
  {
    // ------ Temporal lambda function to prevent error in build
    auto getMothersIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lMothersIndeces{};
      for (auto& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
        lMothersIndeces.push_back(lMother.globalIndex());
      }
      return lMothersIndeces;
    };
    auto getMothersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lMothersPDGs{};
      for (auto& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother pdgcode lMother: %d", lMother.pdgCode());
        lMothersPDGs.push_back(lMother.pdgCode());
      }
      return lMothersPDGs;
    };
    // ------
    std::vector<int> mothers = {-1, -1};
    std::vector<int> motherPDGs = {-1, -1};
    if (track.has_mcParticle()) {
      //
      // Get the MC particle
      const auto& particle = track.mcParticle();
      if (particle.has_mothers()) {
        mothers = getMothersIndeces(particle);
        motherPDGs = getMothersPDGCodes(particle);
      }
      while (mothers.size() > 2) {
        mothers.pop_back();
        motherPDGs.pop_back();
      }
      reso2mctracks(particle.pdgCode(),
                    mothers[0],
                    motherPDGs[0],
                    particle.isPhysicalPrimary(),
                    particle.producedByGenerator());
    } else {
      // No MC particle associated
      reso2mctracks(0,
                    mothers[0],
                    motherPDGs[0],
                    0,
                    0);
    }
  }
  // Additonoal information for MC V0s
  template <typename V0Type>
  void fillMCV0(V0Type const& v0)
  {
    // ------ Temporal lambda function to prevent error in build
    auto getMothersIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lMothersIndeces{};
      for (auto& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
        lMothersIndeces.push_back(lMother.globalIndex());
      }
      return lMothersIndeces;
    };
    auto getMothersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lMothersPDGs{};
      for (auto& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother pdgcode lMother: %d", lMother.pdgCode());
        lMothersPDGs.push_back(lMother.pdgCode());
      }
      return lMothersPDGs;
    };
    auto getDaughtersIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lDaughtersIndeces{};
      for (auto& lDaughter : theMcParticle.template daughters_as<aod::McParticles>()) {
        LOGF(debug, "   daughter index lDaughter: %d", lDaughter.globalIndex());
        if (lDaughter.globalIndex() != 0) {
          lDaughtersIndeces.push_back(lDaughter.globalIndex());
        }
      }
      return lDaughtersIndeces;
    };
    auto getDaughtersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lDaughtersPDGs{};
      for (auto& lDaughter : theMcParticle.template daughters_as<aod::McParticles>()) {
        LOGF(debug, "   daughter pdgcode lDaughter: %d", lDaughter.pdgCode());
        if (lDaughter.globalIndex() != 0) {
          lDaughtersPDGs.push_back(lDaughter.pdgCode());
        }
      }
      return lDaughtersPDGs;
    };
    // ------
    std::vector<int> mothers = {-1, -1};
    std::vector<int> motherPDGs = {-1, -1};
    std::vector<int> daughters = {-1, -1};
    std::vector<int> daughterPDGs = {-1, -1};
    if (v0.has_mcParticle()) {
      auto v0mc = v0.mcParticle();
      if (v0mc.has_mothers()) {
        mothers = getMothersIndeces(v0mc);
        motherPDGs = getMothersPDGCodes(v0mc);
      }
      while (mothers.size() > 2) {
        mothers.pop_back();
        motherPDGs.pop_back();
      }
      if (v0mc.has_daughters()) {
        daughters = getDaughtersIndeces(v0mc);
        daughterPDGs = getDaughtersPDGCodes(v0mc);
      }
      while (daughters.size() > 2) {
        LOGF(info, "daughters.size() is larger than 2");
        daughters.pop_back();
        daughterPDGs.pop_back();
      }
      reso2mcv0s(v0mc.pdgCode(),
                 mothers[0],
                 motherPDGs[0],
                 daughters[0],
                 daughters[1],
                 daughterPDGs[0],
                 daughterPDGs[1],
                 v0mc.isPhysicalPrimary(),
                 v0mc.producedByGenerator());
    } else {
      reso2mcv0s(0,
                 mothers[0],
                 motherPDGs[0],
                 daughters[0],
                 daughters[1],
                 daughterPDGs[0],
                 daughterPDGs[1],
                 0,
                 0);
    }
  }
  // Additonoal information for MC Cascades
  template <typename CascType>
  void fillMCCascade(CascType const& casc)
  {
    // ------ Temporal lambda function to prevent error in build
    auto getMothersIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lMothersIndeces{};
      for (auto& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
        lMothersIndeces.push_back(lMother.globalIndex());
      }
      return lMothersIndeces;
    };
    auto getMothersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lMothersPDGs{};
      for (auto& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother pdgcode lMother: %d", lMother.pdgCode());
        lMothersPDGs.push_back(lMother.pdgCode());
      }
      return lMothersPDGs;
    };
    auto getDaughtersIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lDaughtersIndeces{};
      for (auto& lDaughter : theMcParticle.template daughters_as<aod::McParticles>()) {
        LOGF(debug, "   daughter index lDaughter: %d", lDaughter.globalIndex());
        if (lDaughter.globalIndex() != 0) {
          lDaughtersIndeces.push_back(lDaughter.globalIndex());
        }
      }
      return lDaughtersIndeces;
    };
    auto getDaughtersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lDaughtersPDGs{};
      for (auto& lDaughter : theMcParticle.template daughters_as<aod::McParticles>()) {
        LOGF(debug, "   daughter pdgcode lDaughter: %d", lDaughter.pdgCode());
        if (lDaughter.globalIndex() != 0) {
          lDaughtersPDGs.push_back(lDaughter.pdgCode());
        }
      }
      return lDaughtersPDGs;
    };
    // ------
    std::vector<int> mothers = {-1, -1};
    std::vector<int> motherPDGs = {-1, -1};
    std::vector<int> daughters = {-1, -1};
    std::vector<int> daughterPDGs = {-1, -1};
    if (casc.has_mcParticle()) {
      auto cascmc = casc.mcParticle();
      if (cascmc.has_mothers()) {
        mothers = getMothersIndeces(cascmc);
        motherPDGs = getMothersPDGCodes(cascmc);
      }
      while (mothers.size() > 2) {
        mothers.pop_back();
        motherPDGs.pop_back();
      }
      if (cascmc.has_daughters()) {
        daughters = getDaughtersIndeces(cascmc);
        daughterPDGs = getDaughtersPDGCodes(cascmc);
      }
      while (daughters.size() > 2) {
        LOGF(info, "daughters.size() is larger than 2");
        daughters.pop_back();
        daughterPDGs.pop_back();
      }
      reso2mccascades(cascmc.pdgCode(),
                      mothers[0],
                      motherPDGs[0],
                      daughters[0],
                      daughters[1],
                      daughterPDGs[0],
                      daughterPDGs[1],
                      cascmc.isPhysicalPrimary(),
                      cascmc.producedByGenerator());
    } else {
      reso2mccascades(0,
                      mothers[0],
                      motherPDGs[0],
                      daughters[0],
                      daughters[1],
                      daughterPDGs[0],
                      daughterPDGs[1],
                      0,
                      0);
    }
  }
  // Additonoal information for MC Cascades
  template <typename SelectedMCPartType, typename TotalMCParts>
  void fillMCParticles(SelectedMCPartType const& mcParts, TotalMCParts const& mcParticles)
  {
    for (auto& mcPart : mcParts) {
      std::vector<int> daughterPDGs;
      if (mcPart.has_daughters()) {
        auto daughter01 = mcParticles.rawIteratorAt(mcPart.daughtersIds()[0] - mcParticles.offset());
        auto daughter02 = mcParticles.rawIteratorAt(mcPart.daughtersIds()[1] - mcParticles.offset());
        daughterPDGs = {daughter01.pdgCode(), daughter02.pdgCode()};
      } else {
        daughterPDGs = {-1, -1};
      }
      reso2mcparents(resoCollisions.lastIndex(),
                     mcPart.globalIndex(),
                     mcPart.pdgCode(),
                     daughterPDGs[0], daughterPDGs[1],
                     mcPart.isPhysicalPrimary(),
                     mcPart.producedByGenerator(),
                     mcPart.pt(),
                     mcPart.px(),
                     mcPart.py(),
                     mcPart.pz(),
                     mcPart.eta(),
                     mcPart.phi(),
                     mcPart.y());
      daughterPDGs.clear();
    }
  }

  void init(InitContext&)
  {
    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtTriggerSel, ConfEvtOfflineCheck, ConfIsRun3);
    colCuts.init(&qaRegistry);
  }

  void processTrackData(soa::Filtered<ResoEvents>::iterator const& collision,
                        soa::Filtered<ResoTracks> const& tracks,
                        aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    // Default event selection
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    if (ConfIsRun3) {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFT0M(), collision.multTPC(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    } else {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFV0M(), collision.multTPC(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    }

    fillTracks<false>(collision, tracks);
  }
  PROCESS_SWITCH(reso2initializer, processTrackData, "Process for data", true);

  void processTrackV0Data(soa::Filtered<ResoEvents>::iterator const& collision,
                          soa::Filtered<ResoTracks> const& tracks,
                          ResoV0s const& V0s,
                          aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    // Default event selection
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    if (ConfIsRun3) {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFT0M(), collision.multTPC(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    } else {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFV0M(), collision.multTPC(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    }

    fillTracks<false>(collision, tracks);
    fillV0s<false>(collision, V0s, tracks);
  }
  PROCESS_SWITCH(reso2initializer, processTrackV0Data, "Process for data", true);

  void processTrackV0CascData(soa::Filtered<ResoEvents>::iterator const& collision,
                              soa::Filtered<ResoTracks> const& tracks,
                              ResoV0s const& V0s,
                              ResoCascades const& Cascades,
                              aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    // Default event selection
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    if (ConfIsRun3) {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFT0M(), collision.multTPC(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    } else {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFV0M(), collision.multTPC(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    }

    fillTracks<false>(collision, tracks);
    fillV0s<false>(collision, V0s, tracks);
    fillCascades<false>(collision, Cascades, tracks);
  }
  PROCESS_SWITCH(reso2initializer, processTrackV0CascData, "Process for data", true);

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  void processTrackMC(soa::Filtered<soa::Join<ResoEvents, aod::McCollisionLabels>>::iterator const& collision,
                      aod::McCollisions const& mcCols, soa::Filtered<ResoTracksMC> const& tracks,
                      aod::McParticles const& mcParticles, aod::BCsWithTimestamps const& bcs)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    if (ConfIsRun3) {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFT0M(), collision.multTPC(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    } else {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFV0M(), collision.multTPC(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    }

    // Loop over tracks
    fillTracks<true>(collision, tracks);

    // Loop over all MC particles
    auto mcParts = selectedMCParticles->sliceBy(perMcCollision, collision.mcCollision().globalIndex());
    fillMCParticles(mcParts, mcParticles);
  }
  PROCESS_SWITCH(reso2initializer, processTrackMC, "Process for MC", false);

  void processTrackV0MC(soa::Filtered<soa::Join<ResoEvents, aod::McCollisionLabels>>::iterator const& collision,
                        aod::McCollisions const& mcCols, soa::Filtered<ResoTracksMC> const& tracks,
                        ResoV0sMC const& V0s,
                        aod::McParticles const& mcParticles, aod::BCsWithTimestamps const& bcs)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    if (ConfIsRun3) {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFT0M(), collision.multTPC(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    } else {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFV0M(), collision.multTPC(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    }

    // Loop over tracks
    fillTracks<true>(collision, tracks);
    fillV0s<true>(collision, V0s, tracks);

    // Loop over all MC particles
    auto mcParts = selectedMCParticles->sliceBy(perMcCollision, collision.mcCollision().globalIndex());
    fillMCParticles(mcParts, mcParticles);
  }
  PROCESS_SWITCH(reso2initializer, processTrackV0MC, "Process for MC", false);

  void processTrackV0CascMC(soa::Filtered<soa::Join<ResoEvents, aod::McCollisionLabels>>::iterator const& collision,
                            aod::McCollisions const& mcCols, soa::Filtered<ResoTracksMC> const& tracks,
                            ResoV0sMC const& V0s,
                            ResoCascadesMC const& Cascades,
                            aod::McParticles const& mcParticles, aod::BCsWithTimestamps const& bcs)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    if (ConfIsRun3) {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFT0M(), collision.multTPC(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    } else {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFV0M(), collision.multTPC(), colCuts.computeSphericity(collision, tracks), bc.timestamp());
    }

    // Loop over tracks
    fillTracks<true>(collision, tracks);
    fillV0s<true>(collision, V0s, tracks);
    fillV0s<true>(collision, V0s, tracks);
    fillCascades<true>(collision, Cascades, tracks);

    // Loop over all MC particles
    auto mcParts = selectedMCParticles->sliceBy(perMcCollision, collision.mcCollision().globalIndex());
    fillMCParticles(mcParts, mcParticles);
  }
  PROCESS_SWITCH(reso2initializer, processTrackV0CascMC, "Process for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<reso2initializer>(cfgc, TaskName{"lf-reso2initializer"}),
  };
}
