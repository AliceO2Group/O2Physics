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
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/ASoAHelpers.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/Utils/collisionCuts.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

/// Initializer for the resonance candidate producers
struct reso2initializer {
  enum {
    kECbegin = 0,
    kINEL = 1,
    kINEL10,
    kINELg0,
    kINELg010,
    kTrig,
    kINELg0Trig,
    kINELg010Trig,
    kECend,
  };
  SliceCache cache;
  float cXiMass;
  int mRunNumber;
  int multEstimator;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  Produces<aod::ResoCollisions> resoCollisions;
  Produces<aod::ResoTracks> reso2trks;
  Produces<aod::ResoV0s> reso2v0s;
  Produces<aod::ResoCascades> reso2cascades;
  Produces<aod::ResoMCTracks> reso2mctracks;
  Produces<aod::ResoMCParents> reso2mcparents;
  Produces<aod::ResoMCV0s> reso2mcv0s;
  Produces<aod::ResoMCCascades> reso2mccascades;

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  Configurable<bool> cfgFatalWhenNull{"cfgFatalWhenNull", true, "Fatal when null on ccdb access"};

  // Configurables
  Configurable<bool> ConfIsRun3{"ConfIsRun3", true, "Running on Pilot beam"}; // Choose if running on converted data or pilot beam
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<bool> ConfFillQA{"ConfFillQA", false, "Fill QA histograms"};
  Configurable<bool> ConfBypassCCDB{"ConfBypassCCDB", false, "Bypass loading CCDB part to save CPU time and memory"}; // will be affected to b_z value.

  // Track filter from tpcSkimsTableCreator
  Configurable<int> trackSelection{"trackSelection", 0, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<int> trackSphDef{"trackSphDef", 0, "Spherocity Definition: |pT| = 1 -> 0, otherwise -> 1"};
  Configurable<int> trackSphMin{"trackSphMin", 10, "Number of tracks for Spherocity Calculation"};

  // EventCorrection for MC
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 0.01, 0.1, 1.0, 5.0, 10., 15., 20., 30., 40., 50., 70., 100.0, 105.}, "Binning of the centrality axis"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -20, -15, -10, -7, -5, -3, -2, -1, 0, 1, 2, 3, 5, 7, 10, 15, 20}, "Mixing bins - z-vertex"};

  /// Event cuts
  o2::analysis::CollisonCuts colCuts;
  Configurable<float> ConfEvtZvtx{"ConfEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<bool> ConfEvtTriggerCheck{"ConfEvtTriggerCheck", false, "Evt sel: check for trigger"};
  Configurable<int> ConfEvtTriggerSel{"ConfEvtTriggerSel", 8, "Evt sel: trigger"};
  Configurable<bool> ConfEvtOfflineCheck{"ConfEvtOfflineCheck", true, "Evt sel: check for offline selection"};

  Configurable<std::string> cfgMultName{"cfgMultName", "FT0M", "The name of multiplicity estimator"};

  // Pre-selection cuts
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<float> pidnSigmaPreSelectionCut{"pidnSigmaPreSelectionCut", 5.0f, "TPC and TOF PID cut (loose, improve performance)"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};

  /// DCA Selections for V0
  // DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 2.0, "Track DCAr cut to PV Maximum"};
  Configurable<double> cMinV0PosDCArToPVcut{"cMinV0PosDCArToPVcut", 0.05f, "V0 Positive Track DCAr cut to PV Minimum"}; // Pre-selection
  Configurable<double> cMinV0NegDCArToPVcut{"cMinV0NegDCArToPVcut", 0.05f, "V0 Negative Track DCAr cut to PV Minimum"}; // Pre-selection
  // DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};

  Configurable<double> cMinV0Radius{"cMinV0Radius", 0.0, "Minimum V0 radius from PV"};
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

  HistogramRegistry qaRegistry{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // Pre-filters for efficient process
  // Filter tofPIDFilter = aod::track::tofExpMom < 0.f || ((aod::track::tofExpMom > 0.f) && ((nabs(aod::pidtof::tofNSigmaPi) < pidnSigmaPreSelectionCut) || (nabs(aod::pidtof::tofNSigmaKa) < pidnSigmaPreSelectionCut) || (nabs(aod::pidtof::tofNSigmaPr) < pidnSigmaPreSelectionCut))); // TOF
  Filter trackFilter = (trackSelection.node() == 0) || // from tpcSkimsTableCreator
                       ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                       ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                       ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                       ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                       ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));
  Filter tpcPIDFilter = nabs(aod::pidtpc::tpcNSigmaPi) < pidnSigmaPreSelectionCut || nabs(aod::pidtpc::tpcNSigmaKa) < pidnSigmaPreSelectionCut || nabs(aod::pidtpc::tpcNSigmaPr) < pidnSigmaPreSelectionCut; // TPC
  Filter trackEtaFilter = nabs(aod::track::eta) < cfgCutEta;                                                                                                                                                 // Eta cut
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
                                                    || (nabs(aod::mcparticle::pdgCode) == 3124)    // Lambda(1520)
                                                    || (nabs(aod::mcparticle::pdgCode) == 3324)    // Xi(1530)0
                                                    || (nabs(aod::mcparticle::pdgCode) == 10323)   // K1(1270)+
                                                    || (nabs(aod::mcparticle::pdgCode) == 123314)  // Xi(1820)0
                                                    || (nabs(aod::mcparticle::pdgCode) == 123324); // Xi(1820)-0

  using ResoEvents = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As>;
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
    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodTrackIndices"), 0.5);
    // MC case can be handled here
    if constexpr (isMC) {
      // MC check
      if (ConfFillQA)
        qaRegistry.fill(HIST("hGoodMCTrackIndices"), 0.5);
    }
    // DCAxy cut
    if (fabs(track.dcaXY()) > cMaxDCArToPVcut)
      return false;
    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodTrackIndices"), 1.5);
    // DCAz cut
    if (fabs(track.dcaZ()) > cMaxDCAzToPVcut || fabs(track.dcaZ()) < cMinDCAzToPVcut)
      return false;
    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodTrackIndices"), 2.5);

    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodTrackIndices"), 7.5);
    return true;
  }

  template <bool isMC, typename CollisionType, typename V0Type, typename TrackType>
  bool IsV0Selected(CollisionType const& collision, V0Type const& v0, TrackType const& track)
  {
    // V0 selection
    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodV0Indices"), 0.5);

    auto postrack = v0.template posTrack_as<TrackType>();
    auto negtrack = v0.template negTrack_as<TrackType>();

    if (postrack.tpcNClsCrossedRows() < mincrossedrows)
      return false;
    if (negtrack.tpcNClsCrossedRows() < mincrossedrows)
      return false;
    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodV0Indices"), 1.5);

    if (fabs(postrack.dcaXY()) < cMinV0PosDCArToPVcut)
      return false;
    if (fabs(negtrack.dcaXY()) < cMinV0NegDCArToPVcut)
      return false;
    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodV0Indices"), 2.5);

    if ((v0.v0radius() > cMaxV0Radius) || (v0.v0radius() < cMinV0Radius))
      return false;
    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodV0Indices"), 3.5);
    if (v0.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cMinV0CosPA)
      return false;
    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodV0Indices"), 4.5);

    // MC case can be handled here
    if constexpr (isMC) {
      // MC check
      if (ConfFillQA)
        qaRegistry.fill(HIST("hGoodMCV0Indices"), 0.5);
    }
    return true;
  }

  template <bool isMC, typename CollisionType, typename CascType, typename TrackType>
  bool IsCascSelected(CollisionType const& collision, CascType const& casc, TrackType const& track)
  {
    // V0 selection
    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodCascIndices"), 0.5);

    auto trackBach = casc.template bachelor_as<TrackType>();
    // auto trackPos = casc.template posTrack_as<TrackType>();
    // auto trackNeg = casc.template negTrack_as<TrackType>();

    // track cuts
    if (trackBach.tpcNClsCrossedRows() < mincrossedrows_cascbach)
      return false;
    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodCascIndices"), 1.5);

    if (fabs(trackBach.dcaXY()) < cMinCascBachDCArToPVcut)
      return false;
    if (fabs(trackBach.dcaXY()) > cMaxCascBachDCArToPVcut)
      return false;
    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodCascIndices"), 2.5);

    // DCA daugthers
    if (casc.dcaV0daughters() > cMaxCascDCAV0Daughters)
      return false;
    if (casc.dcacascdaughters() > cMaxCascDCACascDaughters)
      return false;
    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodCascIndices"), 3.5);

    // CPA cuts
    if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cMinCascCosPA)
      return false;
    if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cMinCascV0CosPA)
      return false;
    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodCascIndices"), 4.5);

    // V0 radius
    auto v0radius = casc.v0radius();
    if ((v0radius > cMaxCascV0Radius) || (v0radius < cMinCascV0Radius))
      return false;
    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodCascIndices"), 5.5);

    // Casc radius
    auto cascradius = casc.cascradius();
    if ((cascradius > cMaxCascRadius) || (cascradius < cMinCascRadius))
      return false;
    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodCascIndices"), 6.5);

    // Casc mass
    auto cascMass = casc.mXi();
    if (abs(cascMass - cXiMass) > cCascMassResol)
      return false;
    if (ConfFillQA)
      qaRegistry.fill(HIST("hGoodCascIndices"), 7.5);

    // MC case can be handled here
    if constexpr (isMC) {
      // MC check
      if (ConfFillQA)
        qaRegistry.fill(HIST("hGoodMCCascIndices"), 0.5);
    }
    return true;
  }

  // Centralicity estimator selection
  template <typename ResoColl>
  float CentEst(ResoColl ResoEvents)
  {
    float returnValue = -999.0;
    switch (multEstimator) {
      case 0:
        returnValue = ResoEvents.centFT0M();
      case 1:
        returnValue = ResoEvents.centFT0C();
      case 2:
        returnValue = ResoEvents.centFT0A();
      case 99:
        returnValue = ResoEvents.centFV0A();
      default:
        returnValue = ResoEvents.centFT0M();
    }
    return returnValue;
  }

  // Multiplicity estimator selection
  template <typename ResoColl>
  float MultEst(ResoColl ResoEvents)
  {
    float returnValue = -999.0;
    switch (multEstimator) {
      case 0:
        returnValue = ResoEvents.multFT0M();
      case 1:
        returnValue = ResoEvents.multFT0C();
      case 2:
        returnValue = ResoEvents.multFT0A();
      case 99:
        returnValue = ResoEvents.multFV0A();
      default:
        returnValue = ResoEvents.multFT0M();
    }
    return returnValue;
  }

  /// Compute the spherocity of an event
  /// Important here is that the filter on tracks does not interfere here!
  /// In Run 2 we used here global tracks within |eta| < 0.8
  /// \tparam T type of the tracks
  /// \param tracks All tracks
  /// \return value of the spherocity of the event
  template <typename T>
  float ComputeSpherocity(T const& tracks, int nTracksMin, int spdef)
  {
    // if number of tracks is not enough for spherocity estimation.
    int ntrks = tracks.size();
    if (ntrks < nTracksMin)
      return -99.;

    // start computing spherocity

    float ptSum = 0.;
    for (auto const& track : tracks) {
      if (ConfFillQA) {
        qaRegistry.fill(HIST("Phi"), track.phi());
      }
      if (spdef == 0) {
        ptSum += 1.;
      } else {
        ptSum += track.pt();
      }
    }

    float tempSph = 1.;
    for (int i = 0; i < 360 / 0.1; ++i) {
      float sum = 0., pt = 0.;
      float phiparm = (TMath::Pi() * i * 0.1) / 180.;
      float nx = TMath::Cos(phiparm);
      float ny = TMath::Sin(phiparm);
      for (auto const& trk : tracks) {
        pt = trk.pt();
        if (spdef == 0) {
          pt = 1.;
        }
        float phi = trk.phi();
        float px = pt * TMath::Cos(phi);
        float py = pt * TMath::Sin(phi);
        // sum += pt * abs(sin(phiparm - phi));
        sum += TMath::Abs(px * ny - py * nx);
      }
      float sph = TMath::Power((sum / ptSum), 2);
      if (sph < tempSph)
        tempSph = sph;
    }

    return TMath::Power(TMath::Pi() / 2., 2) * tempSph;
  }

  // Filter for all tracks
  template <bool isMC, typename TrackType, typename CollisionType>
  void fillTracks(CollisionType const& collision, TrackType const& tracks)
  {
    // Loop over tracks
    for (auto& track : tracks) {
      if (!IsTrackSelected<isMC>(collision, track))
        continue;
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
                track.hasTOF(),
                track.tpcNSigmaPi(),
                track.tpcNSigmaKa(),
                track.tpcNSigmaPr(),
                track.tofNSigmaPi(),
                track.tofNSigmaKa(),
                track.tofNSigmaPr(),
                track.tpcSignal(),
                track.passedITSRefit(),
                track.passedTPCRefit(),
                track.isGlobalTrackWoDCA(),
                track.isPrimaryTrack(),
                track.isPVContributor(),
                track.tpcCrossedRowsOverFindableCls(),
                track.itsChi2NCl(),
                track.tpcChi2NCl());
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
    auto getSiblingsIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lSiblingsIndeces{};
      for (auto& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
        for (auto& lDaughter : lMother.template daughters_as<aod::McParticles>()) {
          LOGF(debug, "   daughter index lDaughter: %d", lDaughter.globalIndex());
          if (lDaughter.globalIndex() != 0 && lDaughter.globalIndex() != theMcParticle.globalIndex()) {
            lSiblingsIndeces.push_back(lDaughter.globalIndex());
          }
        }
      }
      return lSiblingsIndeces;
    };
    // ------
    std::vector<int> mothers = {-1, -1};
    std::vector<int> motherPDGs = {-1, -1};
    int siblings[2] = {0, 0};
    std::vector<int> siblings_temp = {-1, -1};
    if (track.has_mcParticle()) {
      //
      // Get the MC particle
      const auto& particle = track.mcParticle();
      if (particle.has_mothers()) {
        mothers = getMothersIndeces(particle);
        motherPDGs = getMothersPDGCodes(particle);
        siblings_temp = getSiblingsIndeces(particle);
      }
      while (mothers.size() > 2) {
        mothers.pop_back();
        motherPDGs.pop_back();
      }
      if (siblings_temp.size() > 0)
        siblings[0] = siblings_temp[0];
      if (siblings_temp.size() > 1)
        siblings[1] = siblings_temp[1];
      reso2mctracks(particle.pdgCode(),
                    mothers[0],
                    motherPDGs[0],
                    siblings,
                    particle.isPhysicalPrimary(),
                    particle.producedByGenerator());
    } else {
      // No MC particle associated
      reso2mctracks(0,
                    mothers[0],
                    motherPDGs[0],
                    siblings,
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
    cXiMass = pdg->GetParticle(3312)->Mass();
    mRunNumber = 0;
    d_bz = 0;
    // Multiplicity estimator selection (0: FT0M, 1: FT0C, 2: FT0A, 99: FV0A)
    if (cfgMultName.value == "FT0M") {
      multEstimator = 0;
    } else if (cfgMultName.value == "FT0C") {
      multEstimator = 1;
    } else if (cfgMultName.value == "FT0A") {
      multEstimator = 2;
    } else if (cfgMultName.value == "FV0A") {
      multEstimator = 99;
    } else {
      multEstimator = 0;
    }

    colCuts.setCuts(ConfEvtZvtx, ConfEvtTriggerCheck, ConfEvtTriggerSel, ConfEvtOfflineCheck, ConfIsRun3);
    colCuts.init(&qaRegistry);
    if (!ConfBypassCCDB) {
      ccdb->setURL(ccdburl.value);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(cfgFatalWhenNull);
      uint64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
      ccdb->setCreatedNotAfter(now); // TODO must become global parameter from the train creation time
    }

    // QA histograms
    AxisSpec idxAxis = {8, 0, 8, "Index"};
    if (ConfFillQA) {
      qaRegistry.add("hGoodTrackIndices", "hGoodTrackIndices", kTH1F, {idxAxis});
      qaRegistry.add("hGoodMCTrackIndices", "hGoodMCTrackIndices", kTH1F, {idxAxis});
      qaRegistry.add("hGoodV0Indices", "hGoodV0Indices", kTH1F, {idxAxis});
      qaRegistry.add("hGoodMCV0Indices", "hGoodMCV0Indices", kTH1F, {idxAxis});
      qaRegistry.add("hGoodCascIndices", "hGoodCascIndices", kTH1F, {idxAxis});
      qaRegistry.add("hGoodMCCascIndices", "hGoodMCCascIndices", kTH1F, {idxAxis});
      qaRegistry.add("Phi", "#phi distribution", kTH1F, {{65, -0.1, 6.4}});
    }
    // MC histograms
    if (doprocessMCGenCount) {
      AxisSpec EvtClassAxis = {kECend - 1, kECbegin + 0.5, kECend - 0.5, "", "event class"};
      AxisSpec ZAxis = {CfgVtxBins, "zaxis"};
      AxisSpec CentAxis = {binsCent, "centrality"};
      qaRegistry.add("Event/totalEventGenMC", "totalEventGenMC", {HistType::kTHnSparseF, {EvtClassAxis}});
      qaRegistry.add("Event/hgenzvtx", "evntclass; zvtex", {HistType::kTHnSparseF, {EvtClassAxis, ZAxis, CentAxis}});
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc) // Simple copy from LambdaKzeroFinder.cxx
  {
    if (ConfBypassCCDB)
      return;
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      ;
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
    // Set magnetic field value once known
    LOGF(info, "Bz set to %f for run: ", d_bz, mRunNumber);
  }

  void processTrackData(soa::Filtered<ResoEvents>::iterator const& collision,
                        soa::Filtered<ResoTracks> const& tracks,
                        aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    initCCDB(bc);
    // Default event selection
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    if (ConfIsRun3) {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), CentEst(collision), MultEst(collision), ComputeSpherocity(tracks, trackSphMin, trackSphDef), d_bz, bc.timestamp());
    } else {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFV0M(), MultEst(collision), ComputeSpherocity(tracks, trackSphMin, trackSphDef), d_bz, bc.timestamp());
    }

    fillTracks<false>(collision, tracks);
  }
  PROCESS_SWITCH(reso2initializer, processTrackData, "Process for data", true);

  PresliceUnsorted<ResoEventsMC> perMcCol = aod::mccollisionlabel::mcCollisionId;
  void processMCGenCount(aod::McCollisions const& mccollisions, aod::McParticles const& mcParticles, ResoEventsMC const& mcCols)
  {
    // Mainly referenced from dndeta_hi task in PWG-MM PAG-Multiplicity (Beomkyu Kim)
    for (auto& mcCollision : mccollisions) { // Gen MC Event loop
      auto particles = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      std::vector<Bool_t> bevtc(kECend, false);
      bevtc[kINEL] = true;
      auto posZ = mcCollision.posZ();
      if (std::abs(posZ) < 10)
        bevtc[kINEL10] = true;
      for (auto& particle : particles) {
        if (!particle.isPhysicalPrimary())
          continue;
        auto kp = pdg->GetParticle(particle.pdgCode());
        if (kp != nullptr) {
          if (std::abs(kp->Charge()) >= 3) {    // 3 quarks
            if (std::abs(particle.eta()) < 1) { // INEL>0 definition
              bevtc[kINELg0] = true;
              break;
            }
          }
        }
      }
      if (bevtc[kINELg0] && bevtc[kINEL10])
        bevtc[kINELg010] = true;

      for (auto ievtc = 1u; ievtc < kECend; ievtc++) {
        if (bevtc[ievtc])
          qaRegistry.fill(HIST("Event/totalEventGenMC"), Double_t(ievtc));
      }

      auto collisionsample = mcCols.sliceBy(perMcCol, mcCollision.globalIndex());
      float cent = -1.0;
      if (collisionsample.size() != 1) { // Prevent no reconstructed collision case
        cent = -1;
      } else {
        for (auto& collision : collisionsample) {
          cent = CentEst(collision);
          if (collision.sel8())
            bevtc[kTrig] = true;
        }
      }
      if (bevtc[kTrig] && bevtc[kINELg0])
        bevtc[kINELg0Trig] = true;
      if (bevtc[kINELg0Trig] && bevtc[kINEL10])
        bevtc[kINELg010Trig] = true;
      for (auto ievtc = 1u; ievtc < kECend; ievtc++) {
        if (bevtc[ievtc])
          qaRegistry.fill(HIST("Event/hgenzvtx"), Double_t(ievtc), posZ, cent);
      }
    }
  }
  PROCESS_SWITCH(reso2initializer, processMCGenCount, "Process for MC", false);

  void processTrackV0Data(soa::Filtered<ResoEvents>::iterator const& collision,
                          soa::Filtered<ResoTracks> const& tracks,
                          ResoV0s const& V0s,
                          aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    initCCDB(bc);
    // Default event selection
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    if (ConfIsRun3) {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), CentEst(collision), MultEst(collision), ComputeSpherocity(tracks, trackSphMin, trackSphDef), d_bz, bc.timestamp());
    } else {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFV0M(), MultEst(collision), ComputeSpherocity(tracks, trackSphMin, trackSphDef), d_bz, bc.timestamp());
    }

    fillTracks<false>(collision, tracks);
    fillV0s<false>(collision, V0s, tracks);
  }
  PROCESS_SWITCH(reso2initializer, processTrackV0Data, "Process for data", false);

  void processTrackV0CascData(soa::Filtered<ResoEvents>::iterator const& collision,
                              soa::Filtered<ResoTracks> const& tracks,
                              ResoV0s const& V0s,
                              ResoCascades const& Cascades,
                              aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    initCCDB(bc);
    // Default event selection
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    if (ConfIsRun3) {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), CentEst(collision), MultEst(collision), ComputeSpherocity(tracks, trackSphMin, trackSphDef), d_bz, bc.timestamp());
    } else {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFV0M(), MultEst(collision), ComputeSpherocity(tracks, trackSphMin, trackSphDef), d_bz, bc.timestamp());
    }

    fillTracks<false>(collision, tracks);
    fillV0s<false>(collision, V0s, tracks);
    fillCascades<false>(collision, Cascades, tracks);
  }
  PROCESS_SWITCH(reso2initializer, processTrackV0CascData, "Process for data", false);

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  void processTrackMC(soa::Filtered<soa::Join<ResoEvents, aod::McCollisionLabels>>::iterator const& collision,
                      aod::McCollisions const& mcCols, soa::Filtered<ResoTracksMC> const& tracks,
                      aod::McParticles const& mcParticles, aod::BCsWithTimestamps const& bcs)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    initCCDB(bc);
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    if (ConfIsRun3) {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), CentEst(collision), MultEst(collision), ComputeSpherocity(tracks, trackSphMin, trackSphDef), d_bz, bc.timestamp());
    } else {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFV0M(), MultEst(collision), ComputeSpherocity(tracks, trackSphMin, trackSphDef), d_bz, bc.timestamp());
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
    initCCDB(bc);
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    if (ConfIsRun3) {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), CentEst(collision), MultEst(collision), ComputeSpherocity(tracks, trackSphMin, trackSphDef), d_bz, bc.timestamp());
    } else {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFV0M(), MultEst(collision), ComputeSpherocity(tracks, trackSphMin, trackSphDef), d_bz, bc.timestamp());
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
    initCCDB(bc);
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    if (ConfIsRun3) {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), CentEst(collision), MultEst(collision), ComputeSpherocity(tracks, trackSphMin, trackSphDef), d_bz, bc.timestamp());
    } else {
      resoCollisions(collision.posX(), collision.posY(), collision.posZ(), collision.multFV0M(), MultEst(collision), ComputeSpherocity(tracks, trackSphMin, trackSphDef), d_bz, bc.timestamp());
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
