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
/// \file resonanceInitializer.cxx
/// \brief Initializes variables for the resonance candidate producers
/// \author Bong-Hwi Lim <bong-hwi.lim@cern.ch>
///

#include <string>
#include <vector>
#include "Common/DataModel/PIDResponse.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/Core/EventPlaneHelper.h"
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
#include "CommonConstants/PhysicsConstants.h"
#include "CommonConstants/MathConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;
using namespace o2::constants::math;

/// Initializer for the resonance candidate producers
struct ResonanceInitializer {
  SliceCache cache;
  int mRunNumber;
  int multEstimator;
  float dBz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  Produces<aod::ResoCollisions> resoCollisions;
  Produces<aod::ResoMCCollisions> resoMCCollisions;
  Produces<aod::ResoTracks> reso2trks;
  Produces<aod::ResoV0s> reso2v0s;
  Produces<aod::ResoCascades> reso2cascades;
  Produces<aod::ResoMCTracks> reso2mctracks;
  Produces<aod::ResoMCParents> reso2mcparents;
  Produces<aod::ResoMCV0s> reso2mcv0s;
  Produces<aod::ResoMCCascades> reso2mccascades;

  // CCDB options
  Configurable<std::string> ccdbURL{"ccdbURL", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  Configurable<bool> cfgFatalWhenNull{"cfgFatalWhenNull", true, "Fatal when null on ccdb access"};

  // Configurables
  Configurable<double> dBzInput{"dBzInput", -999, "bz field, -999 is automatic"};
  Configurable<bool> cfgFillQA{"cfgFillQA", false, "Fill QA histograms"};
  Configurable<bool> cfgBypassCCDB{"cfgBypassCCDB", true, "Bypass loading CCDB part to save CPU time and memory"}; // will be affected to b_z value.

  // Track filter from tpcSkimsTableCreator
  Configurable<int> trackSelection{"trackSelection", 0, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};
  Configurable<int> trackSphDef{"trackSphDef", 0, "Spherocity Definition: |pT| = 1 -> 0, otherwise -> 1"};
  Configurable<int> trackSphMin{"trackSphMin", 10, "Number of tracks for Spherocity Calculation"};

  // EventCorrection for MC
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0., 0.01, 0.1, 1.0, 5.0, 10., 15., 20., 30., 40., 50., 60., 70., 80., 90., 100.0, 105.}, "Binning of the centrality axis"};

  /// Event cuts
  o2::analysis::CollisonCuts colCuts;
  Configurable<float> cfgEvtZvtx{"cfgEvtZvtx", 10.f, "Evt sel: Max. z-Vertex (cm)"};
  Configurable<int> cfgEvtOccupancyInTimeRange{"cfgEvtOccupancyInTimeRange", -1, "Evt sel: maximum track occupancy"};
  Configurable<bool> cfgEvtTriggerCheck{"cfgEvtTriggerCheck", false, "Evt sel: check for trigger"};
  Configurable<bool> cfgEvtOfflineCheck{"cfgEvtOfflineCheck", true, "Evt sel: check for offline selection"};
  Configurable<bool> cfgEvtTriggerTVXSel{"cfgEvtTriggerTVXSel", false, "Evt sel: triggerTVX selection (MB)"};
  Configurable<bool> cfgEvtTFBorderCut{"cfgEvtTFBorderCut", false, "Evt sel: apply TF border cut"};
  Configurable<bool> cfgEvtUseITSTPCvertex{"cfgEvtUseITSTPCvertex", false, "Evt sel: use at lease on ITS-TPC track for vertexing"};
  Configurable<bool> cfgEvtCollInTimeRangeNarrow{"cfgEvtCollInTimeRangeNarrow", false, "Evt sel: apply NoCollInTimeRangeNarrow"};
  Configurable<bool> cfgEvtZvertexTimedifference{"cfgEvtZvertexTimedifference", false, "Evt sel: apply Z-vertex time difference"};
  Configurable<bool> cfgEvtPileupRejection{"cfgEvtPileupRejection", false, "Evt sel: apply pileup rejection"};
  Configurable<bool> cfgEvtNoITSROBorderCut{"cfgEvtNoITSROBorderCut", false, "Evt sel: apply NoITSRO border cut"};
  Configurable<bool> cfgEvtRun2AliEventCuts{"cfgEvtRun2AliEventCuts", true, "Evt sel: apply Run2 AliEventCuts"};
  Configurable<bool> cfgEvtRun2INELgtZERO{"cfgEvtRun2INELgtZERO", false, "Evt sel: apply Run2 INELgtZERO"};

  Configurable<std::string> cfgMultName{"cfgMultName", "FT0M", "The name of multiplicity estimator"};

  // Qvector configuration
  Configurable<bool> cfgBypassQvec{"cfgBypassQvec", true, "Bypass for qvector task"};
  Configurable<int> cfgEvtPl{"cfgEvtPl", 40500, "Configuration of three subsystems for the event plane and its resolution, 10000*RefA + 100*RefB + S, where FT0C:0, FT0A:1, FT0M:2, FV0A:3, BPos:5, BNeg:6"};

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
  Configurable<int> cfgMinCrossedRowsCascBach{"cfgMinCrossedRowsCascBach", 70, "min crossed rows for bachelor track from cascade"};
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
  // Filter tpcPIDFilter = nabs(aod::pidtpc::tpcNSigmaPi) < pidnSigmaPreSelectionCut || nabs(aod::pidtpc::tpcNSigmaKa) < pidnSigmaPreSelectionCut || nabs(aod::pidtpc::tpcNSigmaPr) < pidnSigmaPreSelectionCut; // TPC
  Filter trackFilter = (trackSelection.node() == 0) || // from tpcSkimsTableCreator
                       ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                       ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                       ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                       ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                       ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));
  Filter trackEtaFilter = nabs(aod::track::eta) < cfgCutEta; // Eta cut

  EventPlaneHelper helperEP;

  int evtPlRefAId = static_cast<int>(cfgEvtPl / 10000);
  int evtPlRefBId = static_cast<int>((cfgEvtPl - evtPlRefAId * 10000) / 100);
  int evtPlDetId = cfgEvtPl - evtPlRefAId * 10000 - evtPlRefBId * 100;

  // MC Resonance parent filter
  Partition<aod::McParticles> selectedMCParticles = (nabs(aod::mcparticle::pdgCode) == 313)        // K*
                                                    || (nabs(aod::mcparticle::pdgCode) == 323)     // K*pm
                                                    || (nabs(aod::mcparticle::pdgCode) == 333)     // phi
                                                    || (nabs(aod::mcparticle::pdgCode) == 9010221) // f_0(980)
                                                    || (nabs(aod::mcparticle::pdgCode) == 10221)   // f_0(1370)
                                                    || (nabs(aod::mcparticle::pdgCode) == 9030221) // f_0(1500)
                                                    || (nabs(aod::mcparticle::pdgCode) == 10331)   // f_0(1710)
                                                    || (nabs(aod::mcparticle::pdgCode) == 20223)   // f_1(1285)
                                                    || (nabs(aod::mcparticle::pdgCode) == 20333)   // f_1(1420)
                                                    || (nabs(aod::mcparticle::pdgCode) == 335)     // f_1(1525)
                                                    || (nabs(aod::mcparticle::pdgCode) == 113)     // rho(770)
                                                    || (nabs(aod::mcparticle::pdgCode) == 213)     // rho(770)pm
                                                    || (nabs(aod::mcparticle::pdgCode) == 3224)    // Sigma(1385)+
                                                    || (nabs(aod::mcparticle::pdgCode) == 102134)  // Lambda(1520)
                                                    || (nabs(aod::mcparticle::pdgCode) == 3324)    // Xi(1530)0
                                                    || (nabs(aod::mcparticle::pdgCode) == 10323)   // K1(1270)+
                                                    || (nabs(aod::mcparticle::pdgCode) == 123314)  // Xi(1820)0
                                                    || (nabs(aod::mcparticle::pdgCode) == 123324); // Xi(1820)-0

  using ResoEvents = soa::Join<aod::Collisions, aod::EvSels, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>;
  using ResoRun2Events = soa::Join<aod::Collisions, aod::EvSels, aod::CentRun2V0Ms>;
  using ResoEventsMC = soa::Join<ResoEvents, aod::McCollisionLabels>;
  using ResoRun2EventsMC = soa::Join<ResoEvents, aod::McCollisionLabels>;
  using ResoTracks = aod::Reso2TracksPIDExt;
  using ResoTracksMC = soa::Join<ResoTracks, aod::McTrackLabels>;
  using ResoV0s = aod::V0Datas;
  using ResoV0sMC = soa::Join<ResoV0s, aod::McV0Labels>;
  using ResoCascades = aod::CascDatas;
  using ResoCascadesMC = soa::Join<ResoCascades, aod::McCascLabels>;
  using BCsWithRun2Info = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps>;

  template <bool isMC, typename CollisionType, typename TrackType>
  bool isTrackSelected(CollisionType const&, TrackType const& track)
  {
    // Track selection
    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodTrackIndices"), 0.5);
    // MC case can be handled here
    if constexpr (isMC) {
      // MC check
      if (cfgFillQA)
        qaRegistry.fill(HIST("hGoodMCTrackIndices"), 0.5);
    }
    // DCAxy cut
    if (std::fabs(track.dcaXY()) > cMaxDCArToPVcut)
      return false;
    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodTrackIndices"), 1.5);
    // DCAz cut
    if (std::fabs(track.dcaZ()) > cMaxDCAzToPVcut || std::fabs(track.dcaZ()) < cMinDCAzToPVcut)
      return false;
    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodTrackIndices"), 2.5);

    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodTrackIndices"), 7.5);
    return true;
  }

  template <bool isMC, typename CollisionType, typename V0Type, typename TrackType>
  bool isV0Selected(CollisionType const&, V0Type const& v0, TrackType const&)
  {
    // V0 selection
    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodV0Indices"), 0.5);

    auto postrack = v0.template posTrack_as<TrackType>();
    auto negtrack = v0.template negTrack_as<TrackType>();

    if (postrack.tpcNClsCrossedRows() < mincrossedrows)
      return false;
    if (negtrack.tpcNClsCrossedRows() < mincrossedrows)
      return false;
    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodV0Indices"), 1.5);

    if (std::fabs(postrack.dcaXY()) < cMinV0PosDCArToPVcut)
      return false;
    if (std::fabs(negtrack.dcaXY()) < cMinV0NegDCArToPVcut)
      return false;
    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodV0Indices"), 2.5);

    if ((v0.v0radius() > cMaxV0Radius) || (v0.v0radius() < cMinV0Radius))
      return false;
    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodV0Indices"), 3.5);
    if (v0.v0cosPA() < cMinV0CosPA)
      return false;
    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodV0Indices"), 4.5);

    // MC case can be handled here
    if constexpr (isMC) {
      // MC check
      if (cfgFillQA)
        qaRegistry.fill(HIST("hGoodMCV0Indices"), 0.5);
    }
    return true;
  }

  template <bool isMC, typename CollisionType, typename CascType, typename TrackType>
  bool isCascSelected(CollisionType const& collision, CascType const& casc, TrackType const&)
  {
    // V0 selection
    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodCascIndices"), 0.5);

    auto trackBach = casc.template bachelor_as<TrackType>();
    // auto trackPos = casc.template posTrack_as<TrackType>();
    // auto trackNeg = casc.template negTrack_as<TrackType>();

    // track cuts
    if (trackBach.tpcNClsCrossedRows() < cfgMinCrossedRowsCascBach)
      return false;
    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodCascIndices"), 1.5);

    if (std::fabs(trackBach.dcaXY()) < cMinCascBachDCArToPVcut)
      return false;
    if (std::fabs(trackBach.dcaXY()) > cMaxCascBachDCArToPVcut)
      return false;
    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodCascIndices"), 2.5);

    // DCA daugthers
    if (casc.dcaV0daughters() > cMaxCascDCAV0Daughters)
      return false;
    if (casc.dcacascdaughters() > cMaxCascDCACascDaughters)
      return false;
    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodCascIndices"), 3.5);

    // CPA cuts
    if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cMinCascCosPA)
      return false;
    if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cMinCascV0CosPA)
      return false;
    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodCascIndices"), 4.5);

    // V0 radius
    auto v0radius = casc.v0radius();
    if ((v0radius > cMaxCascV0Radius) || (v0radius < cMinCascV0Radius))
      return false;
    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodCascIndices"), 5.5);

    // Casc radius
    auto cascradius = casc.cascradius();
    if ((cascradius > cMaxCascRadius) || (cascradius < cMinCascRadius))
      return false;
    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodCascIndices"), 6.5);

    // Casc mass
    auto cascMass = casc.mXi();
    if (std::abs(cascMass - MassXiMinus) > cCascMassResol)
      return false;
    if (cfgFillQA)
      qaRegistry.fill(HIST("hGoodCascIndices"), 7.5);

    // MC case can be handled here
    if constexpr (isMC) {
      // MC check
      if (cfgFillQA)
        qaRegistry.fill(HIST("hGoodMCCascIndices"), 0.5);
    }
    return true;
  }

  // Check if the collision is INEL>0
  template <typename MCColl, typename MCPart>
  bool isTrueINEL0(MCColl const& /*mccoll*/, MCPart const& mcparts)
  {
    for (auto const& mcparticle : mcparts) {
      if (!mcparticle.isPhysicalPrimary())
        continue;
      auto p = pdg->GetParticle(mcparticle.pdgCode());
      if (p != nullptr) {
        if (std::abs(p->Charge()) >= 3) {
          if (std::abs(mcparticle.eta()) < 1)
            return true;
        }
      }
    }
    return false;
  }

  // Centralicity estimator selection
  template <typename ResoColl>
  float centEst(ResoColl ResoEvents)
  {
    float returnValue = -999.0;
    switch (multEstimator) {
      case 0:
        returnValue = ResoEvents.centFT0M();
        break;
      case 1:
        returnValue = ResoEvents.centFT0C();
        break;
      case 2:
        returnValue = ResoEvents.centFT0A();
        break;
      case 99:
        returnValue = ResoEvents.centFV0A();
        break;
      default:
        returnValue = ResoEvents.centFT0M();
        break;
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
  float computeSpherocity(T const& tracks, int nTracksMin, int spdef)
  {
    // if number of tracks is not enough for spherocity estimation.
    int ntrks = tracks.size();
    if (ntrks < nTracksMin)
      return -99.;

    // start computing spherocity

    float ptSum = 0.;
    for (auto const& track : tracks) {
      if (cfgFillQA) {
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
      float phiparm = (PI * i * 0.1) / 180.;
      float nx = std::cos(phiparm);
      float ny = std::sin(phiparm);
      for (auto const& trk : tracks) {
        pt = trk.pt();
        if (spdef == 0) {
          pt = 1.;
        }
        float phi = trk.phi();
        float px = pt * std::cos(phi);
        float py = pt * std::sin(phi);
        // sum += pt * abs(sin(phiparm - phi));
        sum += std::abs(px * ny - py * nx);
      }
      float sph = std::pow((sum / ptSum), 2);
      if (sph < tempSph)
        tempSph = sph;
    }

    return std::pow(PIHalf, 2) * tempSph;
  }

  template <typename ResoColl>
  float getEvtPl(ResoColl ResoEvents)
  {
    float returnValue = -999.0;
    if (ResoEvents.qvecAmp()[evtPlDetId] > 1e-8)
      returnValue = helperEP.GetEventPlane(ResoEvents.qvecRe()[evtPlDetId * 4 + 3], ResoEvents.qvecIm()[evtPlDetId * 4 + 3], 2);
    return returnValue;
  }

  template <typename ResoColl>
  float getEvtPlRes(ResoColl ResoEvents, int a, int b)
  {
    float returnValue = -999.0;
    if (ResoEvents.qvecAmp()[a] < 1e-8 || ResoEvents.qvecAmp()[b] < 1e-8)
      return returnValue;
    returnValue = helperEP.GetResolution(helperEP.GetEventPlane(ResoEvents.qvecRe()[a * 4 + 3], ResoEvents.qvecIm()[a * 4 + 3], 2), helperEP.GetEventPlane(ResoEvents.qvecRe()[b * 4 + 3], ResoEvents.qvecIm()[b * 4 + 3], 2), 2);
    return returnValue;
  }

  // Filter for all tracks
  template <bool isMC, typename TrackType, typename CollisionType>
  void fillTracks(CollisionType const& collision, TrackType const& tracks)
  {
    // Loop over tracks
    for (auto const& track : tracks) {
      if (!isTrackSelected<isMC>(collision, track))
        continue;
      reso2trks(resoCollisions.lastIndex(),
                track.globalIndex(),
                track.pt(),
                track.px(),
                track.py(),
                track.pz(),
                track.eta(),
                track.phi(),
                track.sign(),
                (uint8_t)track.tpcNClsCrossedRows(),
                (uint8_t)track.tpcNClsFound(),
                (uint8_t)track.itsNCls(),
                track.dcaXY(),
                track.dcaZ(),
                track.x(),
                track.alpha(),
                track.hasITS(),
                track.hasTPC(),
                track.hasTOF(),
                track.tpcNSigmaPi(),
                track.tpcNSigmaKa(),
                track.tpcNSigmaPr(),
                track.tpcNSigmaEl(),
                track.tofNSigmaPi(),
                track.tofNSigmaKa(),
                track.tofNSigmaPr(),
                track.tofNSigmaEl(),
                track.tpcSignal(),
                track.passedITSRefit(),
                track.passedTPCRefit(),
                track.isGlobalTrackWoDCA(),
                track.isGlobalTrack(),
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
    for (auto const& v0 : v0s) {
      if (!isV0Selected<isMC>(collision, v0, tracks))
        continue;
      childIDs[0] = v0.posTrackId();
      childIDs[1] = v0.negTrackId();
      reso2v0s(resoCollisions.lastIndex(),
               v0.globalIndex(),
               v0.pt(),
               v0.px(),
               v0.py(),
               v0.pz(),
               v0.eta(),
               v0.phi(),
               childIDs,
               v0.template posTrack_as<TrackType>().tpcNSigmaPi(),
               v0.template posTrack_as<TrackType>().tpcNSigmaKa(),
               v0.template posTrack_as<TrackType>().tpcNSigmaPr(),
               v0.template negTrack_as<TrackType>().tpcNSigmaPi(),
               v0.template negTrack_as<TrackType>().tpcNSigmaKa(),
               v0.template negTrack_as<TrackType>().tpcNSigmaPr(),
               v0.template negTrack_as<TrackType>().tofNSigmaPi(),
               v0.template negTrack_as<TrackType>().tofNSigmaKa(),
               v0.template negTrack_as<TrackType>().tofNSigmaPr(),
               v0.template posTrack_as<TrackType>().tofNSigmaPi(),
               v0.template posTrack_as<TrackType>().tofNSigmaKa(),
               v0.template posTrack_as<TrackType>().tofNSigmaPr(),
               v0.v0cosPA(),
               v0.dcaV0daughters(),
               v0.dcapostopv(),
               v0.dcanegtopv(),
               v0.dcav0topv(),
               v0.mLambda(),
               v0.mAntiLambda(),
               v0.mK0Short(),
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
    int childIDs[3] = {0, 0, 0}; // these IDs are necessary to keep track of the children
    for (auto const& casc : cascades) {
      if (!isCascSelected<isMC>(collision, casc, tracks))
        continue;
      childIDs[0] = casc.posTrackId();
      childIDs[1] = casc.negTrackId();
      childIDs[2] = casc.bachelorId();
      reso2cascades(resoCollisions.lastIndex(),
                    casc.globalIndex(),
                    casc.pt(),
                    casc.px(),
                    casc.py(),
                    casc.pz(),
                    casc.eta(),
                    casc.phi(),
                    childIDs,
                    casc.template posTrack_as<TrackType>().tpcNSigmaPi(),
                    casc.template posTrack_as<TrackType>().tpcNSigmaKa(),
                    casc.template posTrack_as<TrackType>().tpcNSigmaPr(),
                    casc.template negTrack_as<TrackType>().tpcNSigmaPi(),
                    casc.template negTrack_as<TrackType>().tpcNSigmaKa(),
                    casc.template negTrack_as<TrackType>().tpcNSigmaPr(),
                    casc.template bachelor_as<TrackType>().tpcNSigmaPi(),
                    casc.template bachelor_as<TrackType>().tpcNSigmaKa(),
                    casc.template bachelor_as<TrackType>().tpcNSigmaPr(),
                    casc.template posTrack_as<TrackType>().tofNSigmaPi(),
                    casc.template posTrack_as<TrackType>().tofNSigmaKa(),
                    casc.template posTrack_as<TrackType>().tofNSigmaPr(),
                    casc.template negTrack_as<TrackType>().tofNSigmaPi(),
                    casc.template negTrack_as<TrackType>().tofNSigmaKa(),
                    casc.template negTrack_as<TrackType>().tofNSigmaPr(),
                    casc.template bachelor_as<TrackType>().tofNSigmaPi(),
                    casc.template bachelor_as<TrackType>().tofNSigmaKa(),
                    casc.template bachelor_as<TrackType>().tofNSigmaPr(),
                    casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()),
                    casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()),
                    casc.dcaV0daughters(),
                    casc.dcacascdaughters(),
                    casc.dcapostopv(),
                    casc.dcanegtopv(),
                    casc.dcabachtopv(),
                    casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()),
                    casc.dcaXYCascToPV(),
                    casc.dcaZCascToPV(),
                    casc.sign(),
                    casc.mLambda(),
                    casc.mXi(),
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
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
        lMothersIndeces.push_back(lMother.globalIndex());
      }
      return lMothersIndeces;
    };
    auto getMothersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lMothersPDGs{};
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother pdgcode lMother: %d", lMother.pdgCode());
        lMothersPDGs.push_back(lMother.pdgCode());
      }
      return lMothersPDGs;
    };
    auto getSiblingsIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lSiblingsIndeces{};
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
        for (auto const& lDaughter : lMother.template daughters_as<aod::McParticles>()) {
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
    std::vector<int> siblingsTemp = {-1, -1};
    if (track.has_mcParticle()) {
      //
      // Get the MC particle
      const auto& particle = track.mcParticle();
      if (particle.has_mothers()) {
        mothers = getMothersIndeces(particle);
        motherPDGs = getMothersPDGCodes(particle);
        siblingsTemp = getSiblingsIndeces(particle);
      }
      while (mothers.size() > 2) {
        mothers.pop_back();
        motherPDGs.pop_back();
      }
      if (siblingsTemp.size() > 0)
        siblings[0] = siblingsTemp[0];
      if (siblingsTemp.size() > 1)
        siblings[1] = siblingsTemp[1];
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
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
        lMothersIndeces.push_back(lMother.globalIndex());
      }
      return lMothersIndeces;
    };
    auto getMothersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lMothersPDGs{};
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother pdgcode lMother: %d", lMother.pdgCode());
        lMothersPDGs.push_back(lMother.pdgCode());
      }
      return lMothersPDGs;
    };
    auto getDaughtersIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lDaughtersIndeces{};
      for (auto const& lDaughter : theMcParticle.template daughters_as<aod::McParticles>()) {
        LOGF(debug, "   daughter index lDaughter: %d", lDaughter.globalIndex());
        if (lDaughter.globalIndex() != 0) {
          lDaughtersIndeces.push_back(lDaughter.globalIndex());
        }
      }
      return lDaughtersIndeces;
    };
    auto getDaughtersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lDaughtersPDGs{};
      for (auto const& lDaughter : theMcParticle.template daughters_as<aod::McParticles>()) {
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
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother index lMother: %d", lMother.globalIndex());
        lMothersIndeces.push_back(lMother.globalIndex());
      }
      return lMothersIndeces;
    };
    auto getMothersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lMothersPDGs{};
      for (auto const& lMother : theMcParticle.template mothers_as<aod::McParticles>()) {
        LOGF(debug, "   mother pdgcode lMother: %d", lMother.pdgCode());
        lMothersPDGs.push_back(lMother.pdgCode());
      }
      return lMothersPDGs;
    };
    auto getDaughtersIndeces = [&](auto const& theMcParticle) {
      std::vector<int> lDaughtersIndeces{};
      for (auto const& lDaughter : theMcParticle.template daughters_as<aod::McParticles>()) {
        LOGF(debug, "   daughter index lDaughter: %d", lDaughter.globalIndex());
        if (lDaughter.globalIndex() != 0) {
          lDaughtersIndeces.push_back(lDaughter.globalIndex());
        }
      }
      return lDaughtersIndeces;
    };
    auto getDaughtersPDGCodes = [&](auto const& theMcParticle) {
      std::vector<int> lDaughtersPDGs{};
      for (auto const& lDaughter : theMcParticle.template daughters_as<aod::McParticles>()) {
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
    for (auto const& mcPart : mcParts) {
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
                     mcPart.y(),
                     mcPart.e(),
                     mcPart.statusCode());
      daughterPDGs.clear();
    }
  }

  template <bool isRun2, typename MCCol, typename MCPart>
  void fillMCCollision(MCCol const& mccol, MCPart const& mcparts, float impactpar = -999.0)
  {
    auto centrality = 0.0;
    if constexpr (!isRun2)
      centrality = centEst(mccol);
    else
      centrality = mccol.centRun2V0M();
    bool inVtx10 = (std::abs(mccol.mcCollision().posZ()) > 10.) ? false : true;
    bool isTrueINELgt0 = isTrueINEL0(mccol, mcparts);
    bool isTriggerTVX = mccol.selection_bit(aod::evsel::kIsTriggerTVX);
    bool isSel8 = mccol.sel8();
    bool isSelected = colCuts.isSelected(mccol);
    resoMCCollisions(inVtx10, isTrueINELgt0, isTriggerTVX, isSel8, isSelected, impactpar);

    // QA for Trigger efficiency
    qaRegistry.fill(HIST("Event/hMCEventIndices"), centrality, aod::resocollision::kINEL);
    if (inVtx10)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), centrality, aod::resocollision::kINEL10);
    if (isTrueINELgt0)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), centrality, aod::resocollision::kINELg0);
    if (inVtx10 && isTrueINELgt0)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), centrality, aod::resocollision::kINELg010);

    // TVX MB trigger
    if (isTriggerTVX)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), centrality, aod::resocollision::kTrig);
    if (isTriggerTVX && inVtx10)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), centrality, aod::resocollision::kTrig10);
    if (isTriggerTVX && isTrueINELgt0)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), centrality, aod::resocollision::kTrigINELg0);
    if (isTriggerTVX && isTrueINELgt0 && inVtx10)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), centrality, aod::resocollision::kTrigINELg010);

    // Sel8 event selection
    if (isSel8)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), centrality, aod::resocollision::kSel8);
    if (isSel8 && inVtx10)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), centrality, aod::resocollision::kSel810);
    if (isSel8 && isTrueINELgt0)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), centrality, aod::resocollision::kSel8INELg0);
    if (isSel8 && isTrueINELgt0 && inVtx10)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), centrality, aod::resocollision::kSel8INELg010);

    // CollisionCuts selection
    if (isSelected)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), centrality, aod::resocollision::kAllCuts);
    if (isSelected && inVtx10)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), centrality, aod::resocollision::kAllCuts10);
    if (isSelected && isTrueINELgt0)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), centrality, aod::resocollision::kAllCutsINELg0);
    if (isSelected && isTrueINELgt0 && inVtx10)
      qaRegistry.fill(HIST("Event/hMCEventIndices"), centrality, aod::resocollision::kAllCutsINELg010);
  }

  void init(InitContext&)
  {
    mRunNumber = 0;
    dBz = 0;
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

    // Check if we are running multiple processes at the same time
    int nProcesses = doprocessTrackData + doprocessTrackV0Data + doprocessTrackV0CascData +
                     doprocessTrackMC + doprocessTrackV0MC + doprocessTrackV0CascMC + doprocessTrackEPData +
                     doprocessTrackDataRun2 + doprocessTrackV0DataRun2 + doprocessTrackV0CascDataRun2 +
                     doprocessTrackMCRun2 + doprocessTrackV0MCRun2 + doprocessTrackV0CascMCRun2;

    if (nProcesses > 1) {
      LOG(fatal) << "Multiple processes are not supported at the same time";
    }

    // Case selector based on the process.
    if (doprocessTrackDataRun2 || doprocessTrackV0DataRun2 || doprocessTrackV0CascDataRun2 || doprocessTrackMCRun2 || doprocessTrackV0MCRun2 || doprocessTrackV0CascMCRun2) {
      colCuts.setCuts(cfgEvtZvtx, cfgEvtTriggerCheck, cfgEvtOfflineCheck, false);
    } else if (doprocessTrackData || doprocessTrackV0Data || doprocessTrackV0CascData || doprocessTrackMC || doprocessTrackV0MC || doprocessTrackV0CascMC || doprocessTrackEPData) {
      colCuts.setCuts(cfgEvtZvtx, cfgEvtTriggerCheck, cfgEvtOfflineCheck, true, false, cfgEvtOccupancyInTimeRange);
    }
    colCuts.init(&qaRegistry);
    colCuts.setTriggerTVX(cfgEvtTriggerTVXSel);
    colCuts.setApplyTFBorderCut(cfgEvtTFBorderCut);
    colCuts.setApplyITSTPCvertex(cfgEvtUseITSTPCvertex);
    colCuts.setApplyCollInTimeRangeNarrow(cfgEvtCollInTimeRangeNarrow);
    colCuts.setApplyZvertexTimedifference(cfgEvtZvertexTimedifference);
    colCuts.setApplyPileupRejection(cfgEvtPileupRejection);
    colCuts.setApplyNoITSROBorderCut(cfgEvtNoITSROBorderCut);
    colCuts.setApplyRun2AliEventCuts(cfgEvtRun2AliEventCuts);
    colCuts.setApplyRun2INELgtZERO(cfgEvtRun2INELgtZERO);
    if (!cfgBypassCCDB) {
      ccdb->setURL(ccdbURL.value);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(cfgFatalWhenNull);
      uint64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
      ccdb->setCreatedNotAfter(now); // TODO must become global parameter from the train creation time
    }

    // QA histograms
    if (doprocessTrackMCRun2 || doprocessTrackV0MCRun2 || doprocessTrackV0CascMCRun2 || doprocessTrackMC || doprocessTrackV0MC || doprocessTrackV0CascMC) {
      AxisSpec centAxis = {binsCent, "Centrality (%)"};
      AxisSpec idxMCAxis = {26, -0.5, 25.5, "Index"};
      qaRegistry.add("Event/hMCEventIndices", "hMCEventIndices", kTH2D, {centAxis, idxMCAxis});
    }
    AxisSpec idxAxis = {8, 0, 8, "Index"};
    if (cfgFillQA) {
      qaRegistry.add("hGoodTrackIndices", "hGoodTrackIndices", kTH1F, {idxAxis});
      qaRegistry.add("hGoodMCTrackIndices", "hGoodMCTrackIndices", kTH1F, {idxAxis});
      qaRegistry.add("hGoodV0Indices", "hGoodV0Indices", kTH1F, {idxAxis});
      qaRegistry.add("hGoodMCV0Indices", "hGoodMCV0Indices", kTH1F, {idxAxis});
      qaRegistry.add("hGoodCascIndices", "hGoodCascIndices", kTH1F, {idxAxis});
      qaRegistry.add("hGoodMCCascIndices", "hGoodMCCascIndices", kTH1F, {idxAxis});
      qaRegistry.add("Phi", "#phi distribution", kTH1F, {{65, -0.1, 6.4}});
    }
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc) // Simple copy from LambdaKzeroFinder.cxx
  {
    if (cfgBypassCCDB)
      return;
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (dBzInput > -990) {
      dBz = dBzInput;
      ;
      o2::parameters::GRPMagField grpmag;
      if (std::fabs(dBz) > 1e-5) {
        grpmag.setL3Current(30000.f / (dBz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3GRPTimestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3GRPTimestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      dBz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3GRPTimestamp << " with magnetic field of " << dBz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3GRPTimestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3GRPTimestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      dBz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3GRPTimestamp << " with magnetic field of " << dBz << " kZG";
    }
    mRunNumber = bc.runNumber();
    // Set magnetic field value once known
    LOGF(info, "Bz set to %f for run: ", dBz, mRunNumber);
  }

  void processDummy(aod::Collisions const& /*collisions*/)
  {
  }
  PROCESS_SWITCH(ResonanceInitializer, processDummy, "Process for dummy", true);

  void processTrackData(ResoEvents::iterator const& collision,
                        soa::Filtered<ResoTracks> const& tracks,
                        aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    initCCDB(bc);
    // Default event selection
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    resoCollisions(collision.globalIndex(), 0, collision.posX(), collision.posY(), collision.posZ(), centEst(collision), computeSpherocity(tracks, trackSphMin, trackSphDef), 0., 0., 0., 0., dBz, bc.timestamp(), collision.trackOccupancyInTimeRange());

    fillTracks<false>(collision, tracks);
  }
  PROCESS_SWITCH(ResonanceInitializer, processTrackData, "Process for data", false);

  void processTrackDataRun2(ResoRun2Events::iterator const& collision,
                            soa::Filtered<ResoTracks> const& tracks,
                            BCsWithRun2Info const&)
  {
    auto bc = collision.bc_as<BCsWithRun2Info>();
    // Default event selection
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQARun2(collision);

    resoCollisions(collision.globalIndex(), 0, collision.posX(), collision.posY(), collision.posZ(), collision.centRun2V0M(), computeSpherocity(tracks, trackSphMin, trackSphDef), 0., 0., 0., 0., dBz, bc.timestamp(), collision.trackOccupancyInTimeRange());

    fillTracks<false>(collision, tracks);
  }
  PROCESS_SWITCH(ResonanceInitializer, processTrackDataRun2, "Process for data", false);

  void processTrackEPData(soa::Join<ResoEvents, aod::Qvectors>::iterator const& collision,
                          soa::Filtered<ResoTracks> const& tracks,
                          aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    initCCDB(bc);
    // Default event selection
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQA(collision);

    resoCollisions(collision.globalIndex(), 0, collision.posX(), collision.posY(), collision.posZ(), centEst(collision), computeSpherocity(tracks, trackSphMin, trackSphDef), getEvtPl(collision), getEvtPlRes(collision, evtPlDetId, evtPlRefAId), getEvtPlRes(collision, evtPlDetId, evtPlRefBId), getEvtPlRes(collision, evtPlRefAId, evtPlRefBId), dBz, bc.timestamp(), collision.trackOccupancyInTimeRange());

    fillTracks<false>(collision, tracks);
  }
  PROCESS_SWITCH(ResonanceInitializer, processTrackEPData, "Process for data and ep ana", false);

  void processTrackV0Data(ResoEvents::iterator const& collision,
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

    resoCollisions(collision.globalIndex(), 0, collision.posX(), collision.posY(), collision.posZ(), centEst(collision), computeSpherocity(tracks, trackSphMin, trackSphDef), 0., 0., 0., 0., dBz, bc.timestamp(), collision.trackOccupancyInTimeRange());

    fillTracks<false>(collision, tracks);
    fillV0s<false>(collision, V0s, tracks);
  }
  PROCESS_SWITCH(ResonanceInitializer, processTrackV0Data, "Process for data", false);

  void processTrackV0DataRun2(ResoRun2Events::iterator const& collision,
                              soa::Filtered<ResoTracks> const& tracks,
                              ResoV0s const& V0s,
                              BCsWithRun2Info const&)
  {
    auto bc = collision.bc_as<BCsWithRun2Info>();
    // Default event selection
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQARun2(collision);

    resoCollisions(collision.globalIndex(), 0, collision.posX(), collision.posY(), collision.posZ(), collision.centRun2V0M(), computeSpherocity(tracks, trackSphMin, trackSphDef), 0., 0., 0., 0., dBz, bc.timestamp(), collision.trackOccupancyInTimeRange());

    fillTracks<false>(collision, tracks);
    fillV0s<false>(collision, V0s, tracks);
  }
  PROCESS_SWITCH(ResonanceInitializer, processTrackV0DataRun2, "Process for data", false);

  void processTrackV0CascData(ResoEvents::iterator const& collision,
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

    resoCollisions(collision.globalIndex(), 0, collision.posX(), collision.posY(), collision.posZ(), centEst(collision), computeSpherocity(tracks, trackSphMin, trackSphDef), 0., 0., 0., 0., dBz, bc.timestamp(), collision.trackOccupancyInTimeRange());

    fillTracks<false>(collision, tracks);
    fillV0s<false>(collision, V0s, tracks);
    fillCascades<false>(collision, Cascades, tracks);
  }
  PROCESS_SWITCH(ResonanceInitializer, processTrackV0CascData, "Process for data", false);

  void processTrackV0CascDataRun2(ResoRun2Events::iterator const& collision,
                                  soa::Filtered<ResoTracks> const& tracks,
                                  ResoV0s const& V0s,
                                  ResoCascades const& Cascades,
                                  BCsWithRun2Info const&)
  {
    auto bc = collision.bc_as<BCsWithRun2Info>();
    // Default event selection
    if (!colCuts.isSelected(collision))
      return;
    colCuts.fillQARun2(collision);

    resoCollisions(collision.globalIndex(), 0, collision.posX(), collision.posY(), collision.posZ(), collision.centRun2V0M(), computeSpherocity(tracks, trackSphMin, trackSphDef), 0., 0., 0., 0., dBz, bc.timestamp(), collision.trackOccupancyInTimeRange());

    fillTracks<false>(collision, tracks);
    fillV0s<false>(collision, V0s, tracks);
    fillCascades<false>(collision, Cascades, tracks);
  }
  PROCESS_SWITCH(ResonanceInitializer, processTrackV0CascDataRun2, "Process for data", false);

  Preslice<aod::McParticles> perMcCollision = aod::mcparticle::mcCollisionId;
  void processTrackMC(soa::Join<ResoEvents, aod::McCollisionLabels>::iterator const& collision,
                      aod::McCollisions const&, soa::Filtered<ResoTracksMC> const& tracks,
                      aod::McParticles const& mcParticles, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    initCCDB(bc);
    colCuts.fillQA(collision);

    resoCollisions(collision.globalIndex(), 0, collision.posX(), collision.posY(), collision.posZ(), centEst(collision), computeSpherocity(tracks, trackSphMin, trackSphDef), 0., 0., 0., 0., dBz, bc.timestamp(), collision.trackOccupancyInTimeRange());

    auto mccollision = collision.mcCollision_as<aod::McCollisions>();
    float impactpar = mccollision.impactParameter();
    fillMCCollision<false>(collision, mcParticles, impactpar);

    // Loop over tracks
    fillTracks<true>(collision, tracks);

    // Loop over all MC particles
    auto mcParts = selectedMCParticles->sliceBy(perMcCollision, collision.mcCollision().globalIndex());
    fillMCParticles(mcParts, mcParticles);
  }
  PROCESS_SWITCH(ResonanceInitializer, processTrackMC, "Process for MC", false);

  void processTrackEPMC(soa::Join<ResoEvents, aod::Qvectors, aod::McCollisionLabels>::iterator const& collision,
                        aod::McCollisions const&, soa::Filtered<ResoTracksMC> const& tracks,
                        aod::McParticles const& mcParticles, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    initCCDB(bc);
    colCuts.fillQA(collision);

    resoCollisions(collision.globalIndex(), 0, collision.posX(), collision.posY(), collision.posZ(), centEst(collision), computeSpherocity(tracks, trackSphMin, trackSphDef), getEvtPl(collision), getEvtPlRes(collision, evtPlDetId, evtPlRefAId), getEvtPlRes(collision, evtPlDetId, evtPlRefBId), getEvtPlRes(collision, evtPlRefAId, evtPlRefBId), dBz, bc.timestamp(), collision.trackOccupancyInTimeRange());
    fillMCCollision<false>(collision, mcParticles);

    // Loop over tracks
    fillTracks<false>(collision, tracks);
    // Loop over all MC particles
    auto mcParts = selectedMCParticles->sliceBy(perMcCollision, collision.mcCollision().globalIndex());
    fillMCParticles(mcParts, mcParticles);
  }
  PROCESS_SWITCH(ResonanceInitializer, processTrackEPMC, "Process for MC and ep ana", false);

  Preslice<aod::McParticles> perMcCollisionRun2 = aod::mcparticle::mcCollisionId;
  void processTrackMCRun2(soa::Join<ResoRun2Events, aod::McCollisionLabels>::iterator const& collision,
                          aod::McCollisions const&, soa::Filtered<ResoTracksMC> const& tracks,
                          aod::McParticles const& mcParticles, BCsWithRun2Info const&)
  {
    auto bc = collision.bc_as<BCsWithRun2Info>();
    colCuts.fillQARun2(collision);

    resoCollisions(collision.globalIndex(), 0, collision.posX(), collision.posY(), collision.posZ(), collision.centRun2V0M(), computeSpherocity(tracks, trackSphMin, trackSphDef), 0., 0., 0., 0., dBz, bc.timestamp(), collision.trackOccupancyInTimeRange());
    fillMCCollision<true>(collision, mcParticles);

    // Loop over tracks
    fillTracks<true>(collision, tracks);

    // Loop over all MC particles
    auto mcParts = selectedMCParticles->sliceBy(perMcCollisionRun2, collision.mcCollision().globalIndex());
    fillMCParticles(mcParts, mcParticles);
  }
  PROCESS_SWITCH(ResonanceInitializer, processTrackMCRun2, "Process for MC", false);

  void processTrackV0MC(soa::Join<ResoEvents, aod::McCollisionLabels>::iterator const& collision,
                        aod::McCollisions const&, soa::Filtered<ResoTracksMC> const& tracks,
                        ResoV0sMC const& V0s,
                        aod::McParticles const& mcParticles, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    initCCDB(bc);
    colCuts.fillQA(collision);

    resoCollisions(collision.globalIndex(), 0, collision.posX(), collision.posY(), collision.posZ(), centEst(collision), computeSpherocity(tracks, trackSphMin, trackSphDef), 0., 0., 0., 0., dBz, bc.timestamp(), collision.trackOccupancyInTimeRange());
    fillMCCollision<false>(collision, mcParticles);

    // Loop over tracks
    fillTracks<true>(collision, tracks);
    fillV0s<true>(collision, V0s, tracks);

    // Loop over all MC particles
    auto mcParts = selectedMCParticles->sliceBy(perMcCollision, collision.mcCollision().globalIndex());
    fillMCParticles(mcParts, mcParticles);
  }
  PROCESS_SWITCH(ResonanceInitializer, processTrackV0MC, "Process for MC", false);

  void processTrackV0MCRun2(soa::Join<ResoRun2Events, aod::McCollisionLabels>::iterator const& collision,
                            aod::McCollisions const&, soa::Filtered<ResoTracksMC> const& tracks,
                            ResoV0sMC const& V0s,
                            aod::McParticles const& mcParticles, BCsWithRun2Info const&)
  {
    auto bc = collision.bc_as<BCsWithRun2Info>();
    colCuts.fillQARun2(collision);

    resoCollisions(collision.globalIndex(), 0, collision.posX(), collision.posY(), collision.posZ(), collision.centRun2V0M(), computeSpherocity(tracks, trackSphMin, trackSphDef), 0., 0., 0., 0., dBz, bc.timestamp(), collision.trackOccupancyInTimeRange());
    fillMCCollision<true>(collision, mcParticles);

    // Loop over tracks
    fillTracks<true>(collision, tracks);
    fillV0s<true>(collision, V0s, tracks);

    // Loop over all MC particles
    auto mcParts = selectedMCParticles->sliceBy(perMcCollision, collision.mcCollision().globalIndex());
    fillMCParticles(mcParts, mcParticles);
  }
  PROCESS_SWITCH(ResonanceInitializer, processTrackV0MCRun2, "Process for MC", false);

  void processTrackV0CascMC(soa::Join<ResoEvents, aod::McCollisionLabels>::iterator const& collision,
                            aod::McCollisions const&, soa::Filtered<ResoTracksMC> const& tracks,
                            ResoV0sMC const& V0s,
                            ResoCascadesMC const& Cascades,
                            aod::McParticles const& mcParticles, aod::BCsWithTimestamps const&)
  {
    auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    initCCDB(bc);
    colCuts.fillQA(collision);

    resoCollisions(collision.globalIndex(), 0, collision.posX(), collision.posY(), collision.posZ(), centEst(collision), computeSpherocity(tracks, trackSphMin, trackSphDef), 0., 0., 0., 0., dBz, bc.timestamp(), collision.trackOccupancyInTimeRange());
    fillMCCollision<false>(collision, mcParticles);

    // Loop over tracks
    fillTracks<true>(collision, tracks);
    fillV0s<true>(collision, V0s, tracks);
    fillCascades<true>(collision, Cascades, tracks);

    // Loop over all MC particles
    auto mcParts = selectedMCParticles->sliceBy(perMcCollision, collision.mcCollision().globalIndex());
    fillMCParticles(mcParts, mcParticles);
  }
  PROCESS_SWITCH(ResonanceInitializer, processTrackV0CascMC, "Process for MC", false);

  void processTrackV0CascMCRun2(soa::Join<ResoRun2Events, aod::McCollisionLabels>::iterator const& collision,
                                aod::McCollisions const&, soa::Filtered<ResoTracksMC> const& tracks,
                                ResoV0sMC const& V0s,
                                ResoCascadesMC const& Cascades,
                                aod::McParticles const& mcParticles, BCsWithRun2Info const&)
  {
    auto bc = collision.bc_as<BCsWithRun2Info>();
    colCuts.fillQARun2(collision);

    resoCollisions(collision.globalIndex(), 0, collision.posX(), collision.posY(), collision.posZ(), collision.centRun2V0M(), computeSpherocity(tracks, trackSphMin, trackSphDef), 0., 0., 0., 0., dBz, bc.timestamp(), collision.trackOccupancyInTimeRange());
    fillMCCollision<true>(collision, mcParticles);

    // Loop over tracks
    fillTracks<true>(collision, tracks);
    fillV0s<true>(collision, V0s, tracks);
    fillCascades<true>(collision, Cascades, tracks);

    // Loop over all MC particles
    auto mcParts = selectedMCParticles->sliceBy(perMcCollision, collision.mcCollision().globalIndex());
    fillMCParticles(mcParts, mcParticles);
  }
  PROCESS_SWITCH(ResonanceInitializer, processTrackV0CascMCRun2, "Process for MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ResonanceInitializer>(cfgc),
  };
}
