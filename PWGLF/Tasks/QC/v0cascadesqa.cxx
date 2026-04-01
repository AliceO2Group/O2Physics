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
/// \file v0cascadesqa.cxx
/// \brief QA task for V0s and Cascades
///
///
/// In case of questions please write to:
/// \author Aimeric Landou (aimeric.landou@cern.ch)
/// \author Chiara De Martin (chiara.de.martin@cern.ch)
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \author Romain Schotter (romain.schotter@cern.ch)

#include <cmath>
// #include <cstdlib>
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/CCDB/ctpRateFetcher.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2::aod::rctsel;

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

enum ParticleType : uint8_t { kPhoton = 0,
                              kK0s,
                              kLambda,
                              kAntiLambda,
                              kXiM,
                              kXiP,
                              kOmegaM,
                              kOmegaP };

// using DaughterTracks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::pidTOFPi, aod::pidTPCPi, aod::pidTOFPr, aod::pidTPCPr, aod::pidTOFKa, aod::pidTPCKa>;
using DaughterTracks = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCKa>;

struct v0cascadesQA {

  HistogramRegistry histos_event{"histos_event", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histos_V0{"histos_V0", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histos_Casc{"histos_Casc", {}, OutputObjHandlingPolicy::AnalysisObject};

  // configurable event properties
  bool isMC = false;
  Configurable<bool> doextraanalysis{"doextraanalysis", false, "Add extra histograms"};
  Configurable<bool> doTreatPiToMuon{"doTreatPiToMuon", true, "Take pi decay into muon into account in MC"};
  Configurable<std::string> irSource{"irSource", "", "Estimator of the interaction rate (Recommended: pp --> T0VTX, Pb-Pb --> ZNC hadronic)"};

  struct : ConfigurableGroup {
    std::string prefix = "eventSelections"; // JSON group name
    Configurable<bool> requireSel8{"requireSel8", true, "require sel8 event selection"};
    Configurable<bool> requireTriggerTVX{"requireTriggerTVX", true, "require FT0 vertex (acceptable FT0C-FT0A time difference) at trigger level"};
    Configurable<bool> rejectITSROFBorder{"rejectITSROFBorder", true, "reject events at ITS ROF border (Run 3 only)"};
    Configurable<bool> rejectTFBorder{"rejectTFBorder", true, "reject events at TF border (Run 3 only)"};
    Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "require events with at least one ITS-TPC track (Run 3 only)"};
    Configurable<bool> requireIsGoodZvtxFT0VsPV{"requireIsGoodZvtxFT0VsPV", false, "require events with PV position along z consistent (within 1 cm) between PV reconstructed using tracks and PV using FT0 A-C time difference (Run 3 only)"};
    Configurable<bool> requireIsVertexTOFmatched{"requireIsVertexTOFmatched", false, "require events with at least one of vertex contributors matched to TOF (Run 3 only)"};
    Configurable<bool> requireIsVertexTRDmatched{"requireIsVertexTRDmatched", false, "require events with at least one of vertex contributors matched to TRD (Run 3 only)"};
    Configurable<bool> rejectSameBunchPileup{"rejectSameBunchPileup", false, "reject collisions in case of pileup with another collision in the same foundBC (Run 3 only)"};
    Configurable<bool> requireNoCollInTimeRangeStd{"requireNoCollInTimeRangeStd", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds (Run 3 only)"};
    Configurable<bool> requireNoCollInTimeRangeStrict{"requireNoCollInTimeRangeStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 10 microseconds (Run 3 only)"};
    Configurable<bool> requireNoCollInTimeRangeNarrow{"requireNoCollInTimeRangeNarrow", false, "reject collisions corrupted by the cannibalism, with other collisions within +/- 2 microseconds (Run 3 only)"};
    Configurable<bool> requireNoCollInROFStd{"requireNoCollInROFStd", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF with mult. above a certain threshold (Run 3 only)"};
    Configurable<bool> requireNoCollInROFStrict{"requireNoCollInROFStrict", false, "reject collisions corrupted by the cannibalism, with other collisions within the same ITS ROF (Run 3 only)"};
    Configurable<bool> requireINEL0{"requireINEL0", false, "require INEL>0 event selection"};
    Configurable<bool> requireINEL1{"requireINEL1", false, "require INEL>1 event selection"};

    Configurable<float> maxZVtxPosition{"maxZVtxPosition", 10., "max Z vtx position"};

    Configurable<bool> useEvtSelInDenomEff{"useEvtSelInDenomEff", false, "Consider event selections in the recoed <-> gen collision association for the denominator (or numerator) of the acc. x eff. (or signal loss)?"};
    Configurable<bool> applyZVtxSelOnMCPV{"applyZVtxSelOnMCPV", false, "Apply Z-vtx cut on the PV of the generated collision?"};
    Configurable<bool> useFT0CbasedOccupancy{"useFT0CbasedOccupancy", false, "Use sum of FT0-C amplitudes for estimating occupancy? (if not, use track-based definition)"};
    // fast check on occupancy
    Configurable<float> minOccupancy{"minOccupancy", -1, "minimum occupancy from neighbouring collisions"};
    Configurable<float> maxOccupancy{"maxOccupancy", -1, "maximum occupancy from neighbouring collisions"};
    // fast check on interaction rate
    Configurable<float> minIR{"minIR", -1, "minimum IR collisions"};
    Configurable<float> maxIR{"maxIR", -1, "maximum IR collisions"};
  } eventSelections;

  struct : ConfigurableGroup {
    std::string prefix = "rctConfigurations"; // JSON group name
    Configurable<std::string> cfgRCTLabel{"cfgRCTLabel", "", "Which detector condition requirements? (CBT, CBT_hadronPID, CBT_electronPID, CBT_calo, CBT_muon, CBT_muon_glo)"};
    Configurable<bool> cfgCheckZDC{"cfgCheckZDC", false, "Include ZDC flags in the bit selection (for Pb-Pb only)"};
    Configurable<bool> cfgTreatLimitedAcceptanceAsBad{"cfgTreatLimitedAcceptanceAsBad", false, "reject all events where the detectors relevant for the specified Runlist are flagged as LimitedAcceptance"};
  } rctConfigurations;

  RCTFlagsChecker rctFlagsChecker{rctConfigurations.cfgRCTLabel.value};

  static constexpr float DefaultLifetimeCuts[1][2] = {{30., 20.}};

  struct : ConfigurableGroup {
    std::string prefix = "v0Selections"; // JSON group name
    Configurable<int> v0TypeSelection{"v0TypeSelection", 1, "select on a certain V0 type (leave negative if no selection desired)"};

    // Selection criteria: acceptance
    Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity"};
    Configurable<float> daughterEtaCut{"daughterEtaCut", 0.8, "max eta for daughters"};

    // Standard 5 topological criteria
    Configurable<float> v0cospa{"v0cospa", 0.97, "min V0 CosPA"};
    Configurable<float> dcav0dau{"dcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> dcanegtopv{"dcanegtopv", .05, "min DCA Neg To PV (cm)"};
    Configurable<float> dcapostopv{"dcapostopv", .05, "min DCA Pos To PV (cm)"};
    Configurable<float> v0radius{"v0radius", 1.2, "minimum V0 radius (cm)"};
    Configurable<float> v0radiusMax{"v0radiusMax", 1E5, "maximum V0 radius (cm)"};
    Configurable<LabeledArray<float>> lifetimecut{"lifetimecut", {DefaultLifetimeCuts[0], 2, {"lifetimecutLambda", "lifetimecutK0s"}}, "lifetimecut"};

    // invariant mass selection
    Configurable<float> compMassRejection{"compMassRejection", -1, "Competing mass rejection (GeV/#it{c}^{2})"};

    // Additional selection on the AP plot (exclusive for K0Short)
    // original equation: lArmPt*5>TMath::Abs(lArmAlpha)
    Configurable<float> armPodCut{"armPodCut", 5.0f, "pT * (cut) > |alpha|, AP cut. Negative: no cut"};

    // Track quality
    Configurable<int> minTPCrows{"minTPCrows", 70, "minimum TPC crossed rows"};
    Configurable<int> minITSclusters{"minITSclusters", -1, "minimum ITS clusters"};
    Configurable<float> maxFractionTPCSharedClusters{"maxFractionTPCSharedClusters", 1e+09, "maximum fraction of TPC shared clusters"};
    Configurable<float> maxITSchi2PerNcls{"maxITSchi2PerNcls", 1e+09, "maximum ITS chi2 per clusters"};
    Configurable<float> maxTPCchi2PerNcls{"maxTPCchi2PerNcls", 1e+09, "maximum TPC chi2 per clusters"};
    Configurable<bool> skipTPConly{"skipTPConly", false, "skip V0s comprised of at least one TPC only prong"};
    Configurable<bool> requirePosITSonly{"requirePosITSonly", false, "require that positive track is ITSonly (overrides TPC quality)"};
    Configurable<bool> requireNegITSonly{"requireNegITSonly", false, "require that negative track is ITSonly (overrides TPC quality)"};
    Configurable<bool> rejectPosITSafterburner{"rejectPosITSafterburner", false, "reject positive track formed out of afterburner ITS tracks"};
    Configurable<bool> rejectNegITSafterburner{"rejectNegITSafterburner", false, "reject negative track formed out of afterburner ITS tracks"};
    Configurable<bool> requirePosITSafterburnerOnly{"requirePosITSafterburnerOnly", false, "require positive track formed out of afterburner ITS tracks"};
    Configurable<bool> requireNegITSafterburnerOnly{"requireNegITSafterburnerOnly", false, "require negative track formed out of afterburner ITS tracks"};
    Configurable<bool> requireAtLeastOneHasTOF{"requireAtLeastOneHasTOF", false, "require that at least one of daughter tracks has an associated TOF signal"};
    Configurable<bool> requirePosHasTOF{"requirePosHasTOF", false, "require that positive track has an associated TOF signal. On false, TOF requirement (if any) is only applied IF the track has an associated TOF signal."};
    Configurable<bool> requireNegHasTOF{"requireNegHasTOF", false, "require that negative track has an associated TOF signal. On false, TOF requirement (if any) is only applied IF the track has an associated TOF signal."};
    Configurable<bool> requirePosHasTRD{"requirePosHasTRD", false, "require that positive track is formed out of TRD track."};
    Configurable<bool> requireNegHasTRD{"requireNegHasTRD", false, "require that negative track is formed out of TRD track."};
    Configurable<int> minTRDclusters{"minTRDclusters", -1, "minimum TRD clusters IF track is formed out of a TRD track"};

    // PID (TPC/TOF)
    Configurable<float> tpcPidNsigmaCutLaPr{"tpcPidNsigmaCutLaPr", 5, "tpcPidNsigmaCutLaPr"};
    Configurable<float> tpcPidNsigmaCutLaPi{"tpcPidNsigmaCutLaPi", 5, "tpcPidNsigmaCutLaPi"};
    Configurable<float> tpcPidNsigmaCutK0Pi{"tpcPidNsigmaCutK0Pi", 5, "tpcPidNsigmaCutK0Pi"};

    // PID (TOF)
    Configurable<float> tofPidNsigmaCutLaPr{"tofPidNsigmaCutLaPr", 1e+6, "tofPidNsigmaCutLaPr"};
    Configurable<float> tofPidNsigmaCutLaPi{"tofPidNsigmaCutLaPi", 1e+6, "tofPidNsigmaCutLaPi"};
    Configurable<float> tofPidNsigmaCutK0Pi{"tofPidNsigmaCutK0Pi", 1e+6, "tofPidNsigmaCutK0Pi"};
  } v0Selections;

  struct : ConfigurableGroup {
    std::string prefix = "cascadeSelections"; // JSON group name
    // Selection criteria: acceptance
    Configurable<float> rapidityCut{"rapidityCut", 0.5, "rapidity"};
    Configurable<float> daughterEtaCut{"daughterEtaCut", 0.8, "max eta for daughters"};

    // Standard 6 topological criteria on V0
    Configurable<float> v0cospa{"v0cospa", 0.97, "min V0 CosPA"};
    Configurable<float> dcav0dau{"dcav0dau", 1.0, "max DCA V0 Daughters (cm)"};
    Configurable<float> dcav0topv{"dcav0topv", .05, "min DCA V0 to PV (cm)"};
    Configurable<float> dcapiontopv{"dcapiontopv", .05, "min DCA Pion To PV (cm)"};
    Configurable<float> dcaprotontopv{"dcaprotontopv", .05, "min DCA Proton To PV (cm)"};
    Configurable<float> v0radius{"v0radius", 1.2, "minimum V0 radius (cm)"};
    Configurable<float> v0radiusMax{"v0radiusMax", 1E5, "maximum V0 radius (cm)"};

    // Standard 6 topological criteria on cascades
    Configurable<float> casccospa{"casccospa", 0.97, "min Cascade CosPA"};
    Configurable<float> dcacascdau{"dcacascdau", 1.0, "max DCA Cascade Daughters (cm)"};
    Configurable<float> dcaxybachbaryontopv{"dcaxybachbaryontopv", -1, "DCAxy Bachelor-Baryon to PV (cm)"};
    Configurable<float> bachbaryoncospa{"bachbaryoncospa", -1, "Bachelor-Baryon CosPA"};
    Configurable<float> dcabachtopv{"dcabachtopv", .05, "min DCA Bachelor To PV (cm)"};
    Configurable<float> cascradius{"cascradius", 0.5, "minimum Cascade radius (cm)"};
    Configurable<float> cascradiusMax{"cascradiusMax", 1E5, "maximum Cascade radius (cm)"};
    Configurable<float> cascProperLifeTime{"cascProperLifeTime", 3, "maximum lifetime (ctau)"};

    // invariant mass selection
    Configurable<float> v0MassWindow{"v0MassWindow", 0.008, "#Lambda mass (GeV/#it{c}^{2})"};
    Configurable<float> compMassRejection{"compMassRejection", 0.008, "Competing mass rejection (GeV/#it{c}^{2})"};

    // Track quality
    Configurable<int> minTPCrows{"minTPCrows", 70, "minimum TPC crossed rows"};
    Configurable<int> minITSclusters{"minITSclusters", -1, "minimum ITS clusters"};
    Configurable<bool> skipTPConly{"skipTPConly", false, "skip V0s comprised of at least one TPC only prong"};
    Configurable<bool> requireBachITSonly{"requireBachITSonly", false, "require that bachelor track is ITSonly (overrides TPC quality)"};
    Configurable<bool> requirePosITSonly{"requirePosITSonly", false, "require that positive track is ITSonly (overrides TPC quality)"};
    Configurable<bool> requireNegITSonly{"requireNegITSonly", false, "require that negative track is ITSonly (overrides TPC quality)"};
    Configurable<bool> requireAtLeastOneHasTOF{"requireAtLeastOneHasTOF", false, "require that at least one of daughter tracks has an associated TOF signal"};
    Configurable<bool> requireBachHasTOF{"requireBachHasTOF", false, "require that bachelor track has an associated TOF signal. On false, TOF requirement (if any) is only applied IF the track has an associated TOF signal."};
    Configurable<bool> requirePosHasTOF{"requirePosHasTOF", false, "require that positive track has an associated TOF signal. On false, TOF requirement (if any) is only applied IF the track has an associated TOF signal."};
    Configurable<bool> requireNegHasTOF{"requireNegHasTOF", false, "require that negative track has an associated TOF signal. On false, TOF requirement (if any) is only applied IF the track has an associated TOF signal."};
    Configurable<bool> requireBachHasTRD{"requireBachHasTRD", false, "require that bachelor track is formed out of TRD track."};
    Configurable<bool> requirePosHasTRD{"requirePosHasTRD", false, "require that positive track is formed out of TRD track."};
    Configurable<bool> requireNegHasTRD{"requireNegHasTRD", false, "require that negative track is formed out of TRD track."};
    Configurable<int> minTRDclusters{"minTRDclusters", -1, "minimum TRD clusters IF track is formed out of a TRD track"};

    // PID (TPC/TOF)
    Configurable<float> tpcPidNsigmaCut{"tpcPidNsigmaCut", 5, "tpcPidNsigmaCut"};
    Configurable<float> tofPidNsigmaCutLaPr{"tofPidNsigmaCutLaPr", 1e+6, "tofPidNsigmaCutLaPr"};
    Configurable<float> tofPidNsigmaCutLaPi{"tofPidNsigmaCutLaPi", 1e+6, "tofPidNsigmaCutLaPi"};
    Configurable<float> tofPidNsigmaCutXiPi{"tofPidNsigmaCutXiPi", 1e+6, "tofPidNsigmaCutXiPi"};
    Configurable<float> tofPidNsigmaCutOmKa{"tofPidNsigmaCutOmKa", 1e+6, "tofPidNsigmaCutOmKa"};
  } cascSelections;

  // CCDB options
  struct : ConfigurableGroup {
    std::string prefix = "ccdbConfigurations"; // JSON group name
    Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
    Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
    Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
    Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
    Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
    Configurable<std::string> mVtxPath{"mVtxPath", "GLO/Calib/MeanVertex", "Path of the mean vertex file"};

    // manual
    Configurable<bool> useCustomMagField{"useCustomMagField", false, "Use custom magnetic field value"};
    Configurable<float> customMagField{"customMagField", 5.0f, "Manually set magnetic field"};
  } ccdbConfigurations;

  o2::ccdb::CcdbApi ccdbApi;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  ctpRateFetcher rateFetcher;

  struct : ConfigurableGroup {
    std::string prefix = "axisConfigurations"; // JSON group name
    // configurable binning of histograms
    ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis axisPtCoarse{"axisPtCoarse", {50, 0.0f, 10.0f}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis axisV0CosPA{"axisV0CosPA", {100, 0.95f, 1.0f}, "V0 cos(PA)"};
    ConfigurableAxis axisV0PA{"axisV0PA", {1000, 0.f, 1.0f}, "Pointing angle (rad)"};
    ConfigurableAxis axisV0Radius{"axisV0Radius", {100, 0.0f, 50.0f}, "V0 radius (cm)"};
    ConfigurableAxis axisV0RadiusCoarse{"axisV0RadiusCoarse", {30, 0.0f, 30.0f}, "V0 radius (cm)"};
    ConfigurableAxis axisV0DecayLength{"axisV0DecayLength", {100, 0.0f, 80.0f}, "Decay length (cm)"};
    ConfigurableAxis axisV0DCANegToPV{"axisV0DCANegToPV", {100, -1.0f, 1.0f}, "DCA neg. to PV (cm)"};
    ConfigurableAxis axisV0DCAPosToPV{"axisV0DCAPosToPV", {100, -1.0f, 1.0f}, "DCA pos. to PV (cm)"};
    ConfigurableAxis axisV0DCAV0Dau{"axisV0DCAV0Dau", {100, 0.0f, 2.f}, "DCA between V0 daughters (cm)"};
    ConfigurableAxis axisLifetimeK0s{"axisLifetimeK0s", {100, 0.0f, 40.0f}, "K^{0}_{S} lifetime (cm)"};
    ConfigurableAxis axisLifetimeLambda{"axisLifetimeLambda", {100, 0.0f, 80.0f}, "#Lambda lifetime (cm)"};
    ConfigurableAxis axisDecayLengthK0s{"axisDecayLengthK0s", {100, 0.0f, 40.0f}, "K^{0}_{S} decay length (cm)"};
    ConfigurableAxis axisDecayLengthLambda{"axisDecayLengthLambda", {100, 0.0f, 80.0f}, "#Lambda decay length (cm)"};
    ConfigurableAxis axisV0DCAV0ToPVK0s{"axisV0DCAV0ToPVK0s", {250, 0.0f, 0.25f}, "DCA V0 to PV (cm)"};
    ConfigurableAxis axisV0DCAV0ToPVLambda{"axisV0DCAV0ToPVLambda", {250, 0.0f, 0.25f}, "DCA V0 to PV (cm)"};
    ConfigurableAxis axisInvMassK0s{"axisInvMassK0s", {200, 0.4f, 0.6f}, "Inv. mass (GeV/#it{c}%{2})"};
    ConfigurableAxis axisInvMassLambda{"axisInvMassLambda", {200, 1.07f, 1.17f}, "Inv. mass (GeV/#it{c}%{2})"};
    ConfigurableAxis axisTPCPIDPion{"axisTPCPIDPion", {100, -10.f, 10.f}, "TPC PID Pion N_{#sigma}"};
    ConfigurableAxis axisTPCPIDProton{"axisTPCPIDProton", {100, -10.f, 10.f}, "TPC PID Proton N_{#sigma}"};
    ConfigurableAxis axisEtaFlag{"axisEtaFlag", {3, -1.5f, 1.5f}, "Pseudorapidity"};
    ConfigurableAxis axisEta{"axisEta", {100, -1.0f, 1.0f}, "Pseudorapidity"};
    ConfigurableAxis axisPhi{"axisPhi", {90, 0.0f, constants::math::TwoPI}, "Azimuthal angle (#rad)"};
    ConfigurableAxis axisITSMapDaughters{"axisITSMapDaughters", {8, -0.5f, 7.5f}, ""};

    ConfigurableAxis axisPtCasc{"axisPtCasc", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "#it{p}_{T} (GeV/#it{c})"};
    ConfigurableAxis axisCascDecayLength{"axisCascDecayLength", {200, 0.0f, 20.0f}, "Decay length (cm)"};
    ConfigurableAxis axisCascRadius{"axisCascRadius", {100, 0.0f, 30.0f}, "Decay length (cm)"};
    ConfigurableAxis axisCascCosPA{"axisCascCosPA", {100, 0.95f, 1.0f}, "Casc cos(PA)"};
    ConfigurableAxis axisCascRapidity{"axisCascRapidity", {100, -1.0f, 1.0f}, "Rapidity"};
    ConfigurableAxis axisCascLifetimeXi{"axisLifetimeXi", {100, 0.0f, 40.0f}, "#Xi lifetime (cm)"};
    ConfigurableAxis axisCascLifetimeOmega{"axisLifetimeOmega", {100, 0.0f, 30.0f}, "#Omega lifetime (cm)"};
    ConfigurableAxis axisCascDCAV0Dau{"axisCascDCAV0Dau", {100, 0.0f, 2.f}, "DCA between casc. daughters (cm)"};
    ConfigurableAxis axisCascDCABachToPV{"axisCascDCABachToPV", {100, -1.0f, 1.0f}, "DCA bach. to PV (cm)"};
    ConfigurableAxis axisInvMassXi{"axisInvMassXi", {80, 1.28f, 1.36f}, "Inv. mass (GeV/#it{c}%{2})"};
    ConfigurableAxis axisInvMassOmega{"axisInvMassOmega", {80, 1.63f, 1.71f}, "Inv. mass (GeV/#it{c}%{2})"};

  } axisConfigurations;

  int dauEtaFlag = 0;

  void init(InitContext const&)
  {
    if (doprocessGenerated || doprocessMonteCarlo) {
      isMC = true;
    }
    // setting CCDB service
    ccdb->setURL(ccdbConfigurations.ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setFatalWhenNull(false);

    // Event Counters
    histos_event.add("hEventCounter", "hEventCounter", {HistType::kTH1F, {{2, 0.0f, 2.0f}}});
    histos_event.add("hEventCounterMC", "hEventCounterMC", {HistType::kTH1F, {{2, 0.0f, 2.0f}}});
    histos_event.add("hEventSelection", "hEventSelection", kTH1D, {{24, -0.5f, +23.5f}});
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(1, "All collisions");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(2, "sel8 cut");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(3, "kIsTriggerTVX");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(4, "kNoITSROFrameBorder");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(5, "kNoTimeFrameBorder");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(6, "posZ cut");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(7, "kIsVertexITSTPC");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(8, "kIsGoodZvtxFT0vsPV");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(9, "kIsVertexTOFmatched");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(10, "kIsVertexTRDmatched");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(11, "kNoSameBunchPileup");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(12, "kNoCollInTimeRangeStd");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(13, "kNoCollInTimeRangeStrict");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(14, "kNoCollInTimeRangeNarrow");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(15, "kNoCollInRofStd");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(16, "kNoCollInRofStrict");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(17, "Below min occup.");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(18, "Above max occup.");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(19, "Below min IR");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(20, "Above max IR");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(21, "INEL>0");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(22, "INEL>1");
    histos_event.get<TH1>(HIST("hEventSelection"))->GetXaxis()->SetBinLabel(23, "RCT flags");

    histos_V0.add("CosPA", "CosPA", kTH1F, {axisConfigurations.axisV0CosPA});
    histos_V0.add("Radius", "Radius", kTH1D, {axisConfigurations.axisV0Radius});
    histos_V0.add("DecayLength", "DecayLength", kTH1F, {axisConfigurations.axisV0DecayLength});
    histos_V0.add("DCANegToPV", "DCANegToPV", kTH1F, {axisConfigurations.axisV0DCANegToPV});
    histos_V0.add("DCAPosToPV", "DCAPosToPV", kTH1F, {axisConfigurations.axisV0DCAPosToPV});
    histos_V0.add("DCAV0Daughters", "DCAV0Daughters", kTH1F, {axisConfigurations.axisV0DCAV0Dau});
    histos_V0.add("CosPAK0s", "CosPAK0s", kTH1F, {axisConfigurations.axisV0CosPA});
    histos_V0.add("CosPALambda", "CosPALambda", kTH1F, {axisConfigurations.axisV0CosPA});
    histos_V0.add("CosPAAntiLambda", "CosPAAntiLambda", kTH1F, {axisConfigurations.axisV0CosPA});
    histos_V0.add("RadiusK0s", "RadiusK0s", kTH1D, {axisConfigurations.axisV0Radius});
    histos_V0.add("RadiusLambda", "RadiusLambda", kTH1D, {axisConfigurations.axisV0Radius});
    histos_V0.add("RadiusAntiLambda", "RadiusAntiLambda", kTH1D, {axisConfigurations.axisV0Radius});
    histos_V0.add("DCANegToPVK0s", "DCANegToPVK0s", kTH1F, {axisConfigurations.axisV0DCANegToPV});
    histos_V0.add("DCANegToPVLambda", "DCANegToPVLambda", kTH1F, {axisConfigurations.axisV0DCANegToPV});
    histos_V0.add("DCANegToPVAntiLambda", "DCANegToPVAntiLambda", kTH1F, {axisConfigurations.axisV0DCANegToPV});
    histos_V0.add("DCAPosToPVK0s", "DCAPosToPVK0s", kTH1F, {axisConfigurations.axisV0DCAPosToPV});
    histos_V0.add("DCAPosToPVLambda", "DCAPosToPVLambda", kTH1F, {axisConfigurations.axisV0DCAPosToPV});
    histos_V0.add("DCAPosToPVAntiLambda", "DCAPosToPVAntiLambda", kTH1F, {axisConfigurations.axisV0DCAPosToPV});
    histos_V0.add("DCAV0DaughtersK0s", "DCAV0DaughtersK0s", kTH1F, {axisConfigurations.axisV0DCAV0Dau});
    histos_V0.add("DCAV0DaughtersLambda", "DCAV0DaughtersLambda", kTH1F, {axisConfigurations.axisV0DCAV0Dau});
    histos_V0.add("DCAV0DaughtersAntiLambda", "DCAV0DaughtersAntiLambda", kTH1F, {axisConfigurations.axisV0DCAV0Dau});
    histos_V0.add("LifetimeK0s", "LifetimeK0s", kTH1F, {axisConfigurations.axisLifetimeK0s});
    histos_V0.add("LifetimeLambda", "LifetimeLambda", kTH1F, {axisConfigurations.axisLifetimeLambda});
    histos_V0.add("LifetimeAntiLambda", "LifetimeAntiLambda", kTH1F, {axisConfigurations.axisLifetimeLambda});
    histos_V0.add("DecayLengthK0s", "DecayLengthK0s", kTH1F, {axisConfigurations.axisDecayLengthK0s});
    histos_V0.add("DecayLengthLambda", "DecayLengthLambda", kTH1F, {axisConfigurations.axisDecayLengthLambda});
    histos_V0.add("DecayLengthAntiLambda", "DecayLengthAntiLambda", kTH1F, {axisConfigurations.axisDecayLengthLambda});
    histos_V0.add("DCAV0ToPVK0s", "DCAV0ToPVK0s", kTH1F, {axisConfigurations.axisV0DCAV0ToPVK0s});
    histos_V0.add("DCAV0ToPVLambda", "DCAV0ToPVLambda", kTH1F, {axisConfigurations.axisV0DCAV0ToPVLambda});
    histos_V0.add("DCAV0ToPVAntiLambda", "DCAV0ToPVAntiLambda", kTH1F, {axisConfigurations.axisV0DCAV0ToPVLambda});
    histos_V0.add("InvMassK0s", "InvMassK0s", kTH3F, {axisConfigurations.axisPt, axisConfigurations.axisInvMassK0s, axisConfigurations.axisEtaFlag});
    histos_V0.add("InvMassLambda", "InvMassLambda", kTH3F, {axisConfigurations.axisPt, axisConfigurations.axisInvMassLambda, axisConfigurations.axisEtaFlag});
    histos_V0.add("InvMassAntiLambda", "InvMassAntiLambda", kTH3F, {axisConfigurations.axisPt, axisConfigurations.axisInvMassLambda, axisConfigurations.axisEtaFlag});
    histos_V0.add("TPCPIDPosPionFromK0s", "TPCPIDPosPionFromK0s", kTH2F, {axisConfigurations.axisPt, axisConfigurations.axisTPCPIDPion});
    histos_V0.add("TPCPIDNegPionFromK0s", "TPCPIDNegPionFromK0s", kTH2F, {axisConfigurations.axisPt, axisConfigurations.axisTPCPIDProton});
    histos_V0.add("TPCPIDPionFromLambda", "TPCPIDPionFromLambda", kTH2F, {axisConfigurations.axisPt, axisConfigurations.axisTPCPIDPion});
    histos_V0.add("TPCPIDProtonFromLambda", "TPCPIDProtonFromLambda", kTH2F, {axisConfigurations.axisPt, axisConfigurations.axisTPCPIDProton});
    histos_V0.add("TPCPIDPionFromAntiLambda", "TPCPIDPionFromAntiLambda", kTH2F, {axisConfigurations.axisPt, axisConfigurations.axisTPCPIDPion});
    histos_V0.add("TPCPIDProtonFromAntiLambda", "TPCPIDProtonFromAntiLambda", kTH2F, {axisConfigurations.axisPt, axisConfigurations.axisTPCPIDProton});
    if (doextraanalysis) {
      histos_V0.add("InvMassK0s_Radius", "InvMassK0s_Radius", kTH2F, {axisConfigurations.axisV0Radius, axisConfigurations.axisInvMassK0s});
      histos_V0.add("InvMassLambda_Radius", "InvMassLambda_Radius", kTH2F, {axisConfigurations.axisV0Radius, axisConfigurations.axisInvMassLambda});
      histos_V0.add("InvMassAntiLambda_Radius", "InvMassAntiLambda_Radius", kTH2F, {axisConfigurations.axisV0Radius, axisConfigurations.axisInvMassLambda});
      histos_V0.add("InvMassK0s_EtaDaughters", "InvMassK0s_EtaDaughters", kTH3F, {axisConfigurations.axisEta, axisConfigurations.axisEta, axisConfigurations.axisInvMassK0s});
      histos_V0.add("InvMassLambda_EtaDaughters", "InvMassLambda_EtaDaughters", kTH3F, {axisConfigurations.axisEta, axisConfigurations.axisEta, axisConfigurations.axisInvMassLambda});
      histos_V0.add("InvMassAntiLambda_EtaDaughters", "InvMassAntiLambda_EtaDaughters", kTH3F, {axisConfigurations.axisEta, axisConfigurations.axisEta, axisConfigurations.axisInvMassLambda});
      histos_V0.add("InvMassK0s_Lifetime", "InvMassK0s_Lifetime", kTH2F, {axisConfigurations.axisLifetimeK0s, axisConfigurations.axisInvMassK0s});
      histos_V0.add("InvMassLambda_Lifetime", "InvMassLambda_Lifetime", kTH2F, {axisConfigurations.axisLifetimeLambda, axisConfigurations.axisInvMassLambda});
      histos_V0.add("InvMassAntiLambda_Lifetime", "InvMassAntiLambda_Lifetime", kTH2F, {axisConfigurations.axisLifetimeLambda, axisConfigurations.axisInvMassLambda});
      histos_V0.add("InvMassK0s_PhiDaughters", "InvMassK0s_PhiDaughters", kTH3F, {axisConfigurations.axisPhi, axisConfigurations.axisPhi, axisConfigurations.axisInvMassK0s});
      histos_V0.add("InvMassLambda_PhiDaughters", "InvMassLambda_PhiDaughters", kTH3F, {axisConfigurations.axisPhi, axisConfigurations.axisPhi, axisConfigurations.axisInvMassLambda});
      histos_V0.add("InvMassAntiLambda_PhiDaughters", "InvMassAntiLambda_PhiDaughters", kTH3F, {axisConfigurations.axisPhi, axisConfigurations.axisPhi, axisConfigurations.axisInvMassLambda});
      histos_V0.add("InvMassK0s_ITSMapDaughters", "InvMassK0s_ITSMapDaughters", kTH3F, {axisConfigurations.axisITSMapDaughters, axisConfigurations.axisITSMapDaughters, axisConfigurations.axisInvMassK0s});
      histos_V0.add("InvMassLambda_ITSMapDaughters", "InvMassLambda_ITSMapDaughters", kTH3F, {axisConfigurations.axisITSMapDaughters, axisConfigurations.axisITSMapDaughters, axisConfigurations.axisInvMassLambda});
      histos_V0.add("InvMassAntiLambda_ITSMapDaughters", "InvMassAntiLambda_ITSMapDaughters", kTH3F, {axisConfigurations.axisITSMapDaughters, axisConfigurations.axisITSMapDaughters, axisConfigurations.axisInvMassLambda});
      histos_V0.add("InvMassK0s_PtRadius", "InvMassK0s_PtRadius", kTH3F, {axisConfigurations.axisPtCoarse, axisConfigurations.axisV0RadiusCoarse, axisConfigurations.axisInvMassK0s});
      histos_V0.add("InvMassLambda_PtRadius", "InvMassLambda_PtRadius", kTH3F, {axisConfigurations.axisPtCoarse, axisConfigurations.axisV0RadiusCoarse, axisConfigurations.axisInvMassLambda});
      histos_V0.add("InvMassAntiLambda_PtRadius", "InvMassAntiLambda_PtRadius", kTH3F, {axisConfigurations.axisPtCoarse, axisConfigurations.axisV0RadiusCoarse, axisConfigurations.axisInvMassLambda});
      histos_V0.add("InvMassK0sVsPtVsPA", "InvMassK0sVsPtVsPA", kTH3F, {axisConfigurations.axisPtCoarse, axisConfigurations.axisV0PA, axisConfigurations.axisInvMassK0s});
      histos_V0.add("InvMassLambdaVsPtVsPA", "InvMassLambdaVsPtVsPA", kTH3F, {axisConfigurations.axisPtCoarse, axisConfigurations.axisV0PA, axisConfigurations.axisInvMassLambda});
      histos_V0.add("InvMassAntiLambdaVsPtVsPA", "InvMassAntiLambdaVsPtVsPA", kTH3F, {axisConfigurations.axisPtCoarse, axisConfigurations.axisV0PA, axisConfigurations.axisInvMassLambda});
    }

    histos_Casc.add("QA_CascadeCandidates", "QA_CascadeCandidates", {HistType::kTH1F, {{10, 0.f, 10.f}}});
    histos_Casc.add("CascCosPA", "CascCosPA", kTH2D, {axisConfigurations.axisV0CosPA, {2, -2, 2}});
    histos_Casc.add("V0CosPA", "V0CosPA", kTH2D, {axisConfigurations.axisV0CosPA, {2, -2, 2}});
    histos_Casc.add("V0CosPAToXi", "V0CosPAToXi", kTH2D, {axisConfigurations.axisV0CosPA, {2, -2, 2}});
    histos_Casc.add("CascDecayLength", "CascDecayLength", kTH2D, {axisConfigurations.axisCascDecayLength, {2, -2, 2}});
    histos_Casc.add("CascRadius", "CascRadius", kTH2F, {axisConfigurations.axisCascRadius, {2, -2, 2}});
    histos_Casc.add("V0Radius", "V0Radius", kTH2F, {axisConfigurations.axisV0Radius, {2, -2, 2}});
    histos_Casc.add("CascRapidityXi", "CascRapidityXi", kTH2F, {axisConfigurations.axisCascRapidity, {2, -2, 2}});
    histos_Casc.add("CascRapidityOmega", "CascRapidityOmega", kTH2F, {axisConfigurations.axisCascRapidity, {2, -2, 2}});
    histos_Casc.add("CascLifetimeXi", "CascLifetimeXi", kTH2F, {axisConfigurations.axisCascLifetimeXi, {2, -2, 2}});
    histos_Casc.add("CascLifetimeOmega", "CascLifetimeOmega", kTH2F, {axisConfigurations.axisCascLifetimeOmega, {2, -2, 2}});
    histos_Casc.add("V0Lifetime", "V0Lifetime", kTH2F, {axisConfigurations.axisLifetimeLambda, {2, -2, 2}});
    histos_Casc.add("CascPt", "CascPt", kTH2F, {axisConfigurations.axisPtCasc, {2, -2, 2}});
    histos_Casc.add("DcaV0Daughters", "DcaV0Daughters", kTH2F, {axisConfigurations.axisV0DCAV0Dau, {2, -2, 2}});
    histos_Casc.add("DcaCascDaughters", "DcaCascDaughters", kTH2F, {axisConfigurations.axisCascDCAV0Dau, {2, -2, 2}});
    histos_Casc.add("DcaV0ToPV", "DcaV0ToPV", kTH2F, {axisConfigurations.axisV0DCAV0ToPVLambda, {2, -2, 2}});
    histos_Casc.add("DcaBachToPV", "DcaBachToPV", kTH2F, {axisConfigurations.axisCascDCABachToPV, {2, -2, 2}});
    histos_Casc.add("DcaPosToPV", "DcaPosToPV", kTH2F, {axisConfigurations.axisV0DCAPosToPV, {2, -2, 2}});
    histos_Casc.add("DcaNegToPV", "DcaNegToPV", kTH2F, {axisConfigurations.axisV0DCANegToPV, {2, -2, 2}});
    histos_Casc.add("InvMassLambda", "InvMassLambda", kTH2F, {axisConfigurations.axisInvMassLambda, {2, -2, 2}});
    histos_Casc.add("InvMassXiPlus", "InvMassXiPlus", kTH3F, {axisConfigurations.axisPtCasc, axisConfigurations.axisInvMassXi, {2, -1.f, 1.f}});
    histos_Casc.add("InvMassXiMinus", "InvMassXiMinus", kTH3F, {axisConfigurations.axisPtCasc, axisConfigurations.axisInvMassXi, {2, -1.f, 1.f}});
    histos_Casc.add("InvMassXiPlus_Radius", "InvMassXiPlus_Radius", kTH2F, {axisConfigurations.axisCascRadius, axisConfigurations.axisInvMassXi});
    histos_Casc.add("InvMassXiMinus_Radius", "InvMassXiMinus_Radius", kTH2F, {axisConfigurations.axisCascRadius, axisConfigurations.axisInvMassXi});
    histos_Casc.add("InvMassOmegaPlus", "InvMassOmegaPlus", kTH3F, {axisConfigurations.axisPtCasc, axisConfigurations.axisInvMassOmega, {2, -1.f, 1.f}});
    histos_Casc.add("InvMassOmegaMinus", "InvMassOmegaMinus", kTH3F, {axisConfigurations.axisPtCasc, axisConfigurations.axisInvMassOmega, {2, -1.f, 1.f}});

    if (isMC) {
      histos_event.add("GeneratedV0s", "GeneratedV0s", kTH3D, {{3, 0.0f, 3.0f}, axisConfigurations.axisPt, axisConfigurations.axisV0Radius});
      histos_event.get<TH1>(HIST("GeneratedV0s"))->GetXaxis()->SetBinLabel(1, "K^{0}_{S}");
      histos_event.get<TH1>(HIST("GeneratedV0s"))->GetXaxis()->SetBinLabel(2, "#Lambda");
      histos_event.get<TH1>(HIST("GeneratedV0s"))->GetXaxis()->SetBinLabel(3, "#bar{#Lambda}");
      histos_event.add("GeneratedCascades", "GeneratedCascades", kTH3D, {{4, 0.0f, 4.0f}, axisConfigurations.axisPtCasc, axisConfigurations.axisCascRadius});
      histos_event.get<TH1>(HIST("GeneratedCascades"))->GetXaxis()->SetBinLabel(1, "#Xi^{#minus}");
      histos_event.get<TH1>(HIST("GeneratedCascades"))->GetXaxis()->SetBinLabel(2, "#bar{#Xi}^{+}");
      histos_event.get<TH1>(HIST("GeneratedCascades"))->GetXaxis()->SetBinLabel(3, "#Omega^{#minus}");
      histos_event.get<TH1>(HIST("GeneratedCascades"))->GetXaxis()->SetBinLabel(4, "#bar{#Omega}^{+}");

      histos_V0.add("InvMassK0sTrue", "InvMassK0sTrue", {HistType::kTH3F, {{100, 0.0f, 10.0f}, {100, 0.f, 50.f}, {200, 0.4f, 0.6f}}});
      histos_V0.add("InvMassLambdaTrue", "InvMassLambdaTrue", {HistType::kTH3F, {{100, 0.0f, 10.0f}, {100, 0.f, 50.f}, {200, 1.07f, 1.17f}}});
      histos_V0.add("InvMassAntiLambdaTrue", "InvMassAntiLambdaTrue", {HistType::kTH3F, {{100, 0.0f, 10.0f}, {100, 0.f, 50.f}, {200, 1.07f, 1.17f}}});

      histos_Casc.add("InvMassXiPlusTrue", "InvMassXiPlusTrue", {HistType::kTH3F, {{100, 0.f, 10.f}, {100, 0.f, 50.f}, {80, 1.28f, 1.36f}}});
      histos_Casc.add("InvMassXiMinusTrue", "InvMassXiMinusTrue", {HistType::kTH3F, {{100, 0.f, 10.f}, {100, 0.f, 50.f}, {80, 1.28f, 1.36f}}});
      histos_Casc.add("InvMassOmegaPlusTrue", "InvMassOmegaPlusTrue", {HistType::kTH3F, {{100, 0.f, 10.f}, {100, 0.f, 50.f}, {80, 1.63f, 1.71f}}});
      histos_Casc.add("InvMassOmegaMinusTrue", "InvMassOmegaMinusTrue", {HistType::kTH3F, {{100, 0.f, 10.f}, {100, 0.f, 50.f}, {80, 1.63f, 1.71f}}});
    }

    // Initialise the RCTFlagsChecker
    rctFlagsChecker.init(rctConfigurations.cfgRCTLabel.value, rctConfigurations.cfgCheckZDC, rctConfigurations.cfgTreatLimitedAcceptanceAsBad);
  }

  template <typename TCollision>
  bool isEventAccepted(TCollision collision, bool fillHists)
  // check whether the collision passes our collision selections
  {
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 0. /* all collisions */);
    }

    if (eventSelections.requireSel8 && !collision.sel8()) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 1 /* sel8 collisions */);
    }

    if (eventSelections.requireTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 2 /* FT0 vertex (acceptable FT0C-FT0A time difference) collisions */);
    }

    if (eventSelections.rejectITSROFBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 3 /* Not at ITS ROF border */);
    }

    if (eventSelections.rejectTFBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 4 /* Not at TF border */);
    }

    if (std::abs(collision.posZ()) > eventSelections.maxZVtxPosition) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 5 /* vertex-Z selected */);
    }

    if (eventSelections.requireIsVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 6 /* Contains at least one ITS-TPC track */);
    }

    if (eventSelections.requireIsGoodZvtxFT0VsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 7 /* PV position consistency check */);
    }

    if (eventSelections.requireIsVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 8 /* PV with at least one contributor matched with TOF */);
    }

    if (eventSelections.requireIsVertexTRDmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTRDmatched)) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 9 /* PV with at least one contributor matched with TRD */);
    }

    if (eventSelections.rejectSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 10 /* Not at same bunch pile-up */);
    }

    if (eventSelections.requireNoCollInTimeRangeStd && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard)) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 11 /* No other collision within +/- 2 microseconds or mult above a certain threshold in -4 - -2 microseconds*/);
    }

    if (eventSelections.requireNoCollInTimeRangeStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStrict)) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 12 /* No other collision within +/- 10 microseconds */);
    }

    if (eventSelections.requireNoCollInTimeRangeNarrow && !collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeNarrow)) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 13 /* No other collision within +/- 2 microseconds */);
    }

    if (eventSelections.requireNoCollInROFStd && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStandard)) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 14 /* No other collision within the same ITS ROF with mult. above a certain threshold */);
    }

    if (eventSelections.requireNoCollInROFStrict && !collision.selection_bit(o2::aod::evsel::kNoCollInRofStrict)) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 15 /* No other collision within the same ITS ROF */);
    }

    float collisionOccupancy = eventSelections.useFT0CbasedOccupancy ? collision.ft0cOccupancyInTimeRange() : collision.trackOccupancyInTimeRange();
    if (eventSelections.minOccupancy >= 0 && collisionOccupancy < eventSelections.minOccupancy) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 16 /* Below min occupancy */);
    }

    if (eventSelections.maxOccupancy >= 0 && collisionOccupancy > eventSelections.maxOccupancy) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 17 /* Above max occupancy */);
    }

    // Fetch interaction rate only if required (in order to limit ccdb calls)
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    double interactionRate = ((eventSelections.minIR >= 0 || eventSelections.maxIR >= 0) && !irSource.value.empty()) ? rateFetcher.fetch(ccdb.service, bc.timestamp(), bc.runNumber(), irSource) * 1.e-3 : -1;
    if (eventSelections.minIR >= 0 && interactionRate < eventSelections.minIR) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 18 /* Below min IR */);
    }

    if (eventSelections.maxIR >= 0 && interactionRate > eventSelections.maxIR) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 19 /* Above max IR */);
    }

    if (eventSelections.requireINEL0 && !collision.isInelGt0()) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 20 /* INEL > 0 */);
    }

    if (eventSelections.requireINEL1 && !collision.isInelGt1()) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 21 /* INEL > 1 */);
    }

    if (!rctConfigurations.cfgRCTLabel.value.empty() && !rctFlagsChecker(collision)) {
      return false;
    }
    if (fillHists) {
      histos_event.fill(HIST("hEventSelection"), 22 /* Pass CBT condition */);
    }

    return true;
  }

  template <typename TV0, typename TCollision>
  bool isV0Accepted(TV0 v0, TCollision collision, float rapidity, int v0Type)
  // precalculate this information so that a check is one mask operation, not many
  {
    // Base topological variables
    if (v0.v0radius() < v0Selections.v0radius)
      return false;
    if (v0.v0radius() > v0Selections.v0radiusMax)
      return false;
    if (std::abs(v0.dcapostopv()) < v0Selections.dcapostopv)
      return false;
    if (std::abs(v0.dcanegtopv()) < v0Selections.dcanegtopv)
      return false;
    if (v0.v0cosPA() < v0Selections.v0cospa)
      return false;
    if (v0.dcaV0daughters() > v0Selections.dcav0dau)
      return false;

    // proper lifetime
    if ((v0Type == kLambda || v0Type == kAntiLambda) && v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0 > v0Selections.lifetimecut->get("lifetimecutLambda"))
      return false;
    if (v0Type == kK0s && v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short > v0Selections.lifetimecut->get("lifetimecutK0s"))
      return false;

    // armenteros (for K0s only)
    if (v0Type == kK0s && v0Selections.armPodCut > 1e-4 && v0.qtarm() * v0Selections.armPodCut < std::abs(v0.alpha()))
      return false;

    // rapidity
    if (std::abs(rapidity) > v0Selections.rapidityCut)
      return false;

    // competing mass rejection
    if ((v0Type == kLambda || v0Type == kAntiLambda) && std::fabs(v0.mK0Short() - o2::constants::physics::MassK0Short) < v0Selections.compMassRejection)
      return false;
    if (v0Type == kK0s && std::fabs(v0.mLambda() - o2::constants::physics::MassLambda0) < v0Selections.compMassRejection)
      return false;

    auto posTrackExtra = v0.template posTrack_as<DaughterTracks>();
    auto negTrackExtra = v0.template negTrack_as<DaughterTracks>();

    // ITS quality flags
    bool posIsFromAfterburner = posTrackExtra.isITSAfterburner();
    bool negIsFromAfterburner = negTrackExtra.isITSAfterburner();

    // reject afterburner track or not
    if (v0Selections.rejectPosITSafterburner && posIsFromAfterburner)
      return false;
    if (v0Selections.rejectNegITSafterburner && negIsFromAfterburner)
      return false;

    // keep afterburner track or not
    if (v0Selections.requirePosITSafterburnerOnly && !posIsFromAfterburner)
      return false;
    if (v0Selections.requireNegITSafterburnerOnly && !negIsFromAfterburner)
      return false;

    // check minium ITS clusters
    if (posTrackExtra.itsNCls() < v0Selections.minITSclusters)
      return false;
    if (negTrackExtra.itsNCls() < v0Selections.minITSclusters)
      return false;

    // check maximum ITS chi2 per clusters
    if (posTrackExtra.itsChi2NCl() > v0Selections.maxITSchi2PerNcls)
      return false;
    if (negTrackExtra.itsChi2NCl() > v0Selections.maxITSchi2PerNcls)
      return false;

    // ITS only tag
    if (v0Selections.requirePosITSonly) {
      if (posTrackExtra.tpcNClsCrossedRows() > 0)
        return false;
    }
    if (v0Selections.requireNegITSonly) {
      if (negTrackExtra.tpcNClsCrossedRows() > 0)
        return false;
    }

    // check minimum TPC crossed rows
    if (posTrackExtra.tpcNClsCrossedRows() < v0Selections.minTPCrows)
      return false;
    if (negTrackExtra.tpcNClsCrossedRows() < v0Selections.minTPCrows)
      return false;

    // check maximum TPC chi2 per clusters
    if (posTrackExtra.tpcChi2NCl() > v0Selections.maxTPCchi2PerNcls)
      return false;
    if (negTrackExtra.tpcChi2NCl() > v0Selections.maxTPCchi2PerNcls)
      return false;

    // check the maximum fraction of allowed shared TPC
    if (posTrackExtra.tpcFractionSharedCls() > v0Selections.maxFractionTPCSharedClusters)
      return false;
    if (negTrackExtra.tpcFractionSharedCls() > v0Selections.maxFractionTPCSharedClusters)
      return false;

    // TPC only tag
    if (v0Selections.skipTPConly) {
      if (posTrackExtra.hasTPC())
        return false;
      if (negTrackExtra.hasTPC())
        return false;
    }

    // TPC PID
    if (v0Type == kK0s) {
      if (std::fabs(posTrackExtra.tpcNSigmaPi()) > v0Selections.tpcPidNsigmaCutK0Pi)
        return false;
      if (std::fabs(negTrackExtra.tpcNSigmaPi()) > v0Selections.tpcPidNsigmaCutK0Pi)
        return false;
    }
    if (v0Type == kLambda) {
      if (std::fabs(posTrackExtra.tpcNSigmaPr()) > v0Selections.tpcPidNsigmaCutLaPr)
        return false;
      if (std::fabs(negTrackExtra.tpcNSigmaPi()) > v0Selections.tpcPidNsigmaCutLaPi)
        return false;
    }
    if (v0Type == kAntiLambda) {
      if (std::fabs(posTrackExtra.tpcNSigmaPi()) > v0Selections.tpcPidNsigmaCutLaPi)
        return false;
      if (std::fabs(negTrackExtra.tpcNSigmaPr()) > v0Selections.tpcPidNsigmaCutLaPr)
        return false;
    }

    // TOF Requirement checks
    if (v0Selections.requirePosHasTOF && !v0.positiveHasTOF()) {
      return false;
    }
    if (v0Selections.requireNegHasTOF && !v0.negativeHasTOF()) {
      return false;
    }

    if (v0Selections.requireAtLeastOneHasTOF && !v0.positiveHasTOF() && !v0.negativeHasTOF()) {
      return false;
    }

    // TOF Nsigma
    if (v0Type == kK0s) {
      if (v0Selections.tofPidNsigmaCutK0Pi < 1e+5 && v0.positiveHasTOF() && std::fabs(v0.tofNSigmaK0PiPlus()) > v0Selections.tofPidNsigmaCutK0Pi) {
        return false;
      }
      if (v0Selections.tofPidNsigmaCutK0Pi < 1e+5 && v0.negativeHasTOF() && std::fabs(v0.tofNSigmaK0PiMinus()) > v0Selections.tofPidNsigmaCutK0Pi) {
        return false;
      }
    }
    if (v0Type == kLambda) {
      if (v0Selections.tofPidNsigmaCutLaPr < 1e+5 && v0.positiveHasTOF() && std::fabs(v0.tofNSigmaLaPr()) > v0Selections.tofPidNsigmaCutLaPr) {
        return false;
      }
      if (v0Selections.tofPidNsigmaCutLaPi < 1e+5 && v0.negativeHasTOF() && std::fabs(v0.tofNSigmaLaPi()) > v0Selections.tofPidNsigmaCutLaPi) {
        return false;
      }
    }
    if (v0Type == kAntiLambda) {
      if (v0Selections.tofPidNsigmaCutLaPi < 1e+5 && v0.positiveHasTOF() && std::fabs(v0.tofNSigmaLaPi()) > v0Selections.tofPidNsigmaCutLaPi) {
        return false;
      }
      if (v0Selections.tofPidNsigmaCutLaPr < 1e+5 && v0.negativeHasTOF() && std::fabs(v0.tofNSigmaLaPr()) > v0Selections.tofPidNsigmaCutLaPr) {
        return false;
      }
    }

    // TRD Requirement checks
    if (v0Selections.requirePosHasTRD && !posTrackExtra.hasTRD()) {
      return false;
    }
    if (v0Selections.requireNegHasTRD && !negTrackExtra.hasTRD()) {
      return false;
    }

    int posTRDhits = 0, negTRDhits = 0;
    for (unsigned int i = 0; i <= 5; i++) {
      if (posTrackExtra.trdPattern() & (1 << i)) {
        posTRDhits++;
      }
      if (negTrackExtra.trdPattern() & (1 << i)) {
        negTRDhits++;
      }
    }

    if (posTrackExtra.hasTRD() && posTRDhits < v0Selections.minTRDclusters) {
      return false;
    }
    if (negTrackExtra.hasTRD() && negTRDhits < v0Selections.minTRDclusters) {
      return false;
    }

    return true;
  }

  template <typename TCascade, typename TCollision>
  bool isCascadeSelected(TCascade casc, TCollision collision, float rapidity, int cascType)
  // precalculate this information so that a check is one mask operation, not many
  {
    //
    // Base topological variables
    //

    // v0 radius min/max selections
    if (casc.v0radius() < cascSelections.v0radius)
      return false;
    if (casc.v0radius() > cascSelections.v0radiusMax)
      return false;
    // DCA proton and pion to PV for Lambda and AntiLambda decay hypotheses
    if (casc.sign() < 0) { // Xi- or Omega- --> positive/negative daughter = proton/pion
      if (std::fabs(casc.dcapostopv()) < cascSelections.dcaprotontopv)
        return false;
      if (std::fabs(casc.dcanegtopv()) < cascSelections.dcapiontopv)
        return false;
    } else { // Xi+ or Omega+ --> positive/negative daughter = pion/proton
      if (std::fabs(casc.dcapostopv()) < cascSelections.dcapiontopv)
        return false;
      if (std::fabs(casc.dcanegtopv()) < cascSelections.dcaprotontopv)
        return false;
    }
    // V0 cosine of pointing angle
    if (casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()) < cascSelections.v0cospa)
      return false;
    // DCA between v0 daughters
    if (casc.dcaV0daughters() > cascSelections.dcav0dau)
      return false;
    // DCA V0 to prim vtx
    if (casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()) < cascSelections.dcav0topv)
      return false;

    // casc radius min/max selections
    if (casc.cascradius() < cascSelections.cascradius)
      return false;
    if (casc.cascradius() > cascSelections.cascradiusMax)
      return false;
    // DCA bachelor selection
    if (std::fabs(casc.dcabachtopv()) < cascSelections.dcabachtopv)
      return false;
    // Bachelor-baryon cosPA selection
    if (casc.bachBaryonCosPA() < cascSelections.bachbaryoncospa)
      return false;
    // DCA bachelor-baryon selection
    if (std::fabs(casc.bachBaryonDCAxyToPV()) < cascSelections.dcaxybachbaryontopv)
      return false;
    // casc cosine of pointing angle
    if (casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()) < cascSelections.casccospa)
      return false;
    // DCA between casc daughters
    if (casc.dcacascdaughters() > cascSelections.dcacascdau)
      return false;

    if (casc.sign() < 0 && (cascType == kXiP || cascType == kOmegaP))
      return false;
    if (casc.sign() > 0 && (cascType == kXiM || cascType == kOmegaM))
      return false;

    //
    // proper lifetime
    float distOverTotMom = std::sqrt(std::pow(casc.x() - collision.posX(), 2) + std::pow(casc.y() - collision.posY(), 2) + std::pow(casc.z() - collision.posZ(), 2)) / (casc.p() + 1E-10);
    if ((cascType == kXiM || cascType == kXiP) && distOverTotMom * o2::constants::physics::MassXiMinus > cascSelections.cascProperLifeTime)
      return false;
    if ((cascType == kOmegaM || cascType == kOmegaP) && distOverTotMom * o2::constants::physics::MassOmegaMinus > cascSelections.cascProperLifeTime)
      return false;

    // rapidity
    if (std::fabs(rapidity) > cascSelections.rapidityCut)
      return false;

    //
    // invariant mass window
    //
    if (std::fabs(casc.mLambda() - o2::constants::physics::MassLambda0) > cascSelections.v0MassWindow)
      return false;

    //
    // competing mass rejection
    //
    if ((cascType == kXiM || cascType == kXiP) && std::fabs(casc.mOmega() - o2::constants::physics::MassOmegaMinus) < cascSelections.compMassRejection)
      return false;
    if ((cascType == kOmegaM || cascType == kOmegaP) && std::fabs(casc.mXi() - o2::constants::physics::MassXiMinus) < cascSelections.compMassRejection)
      return false;

    auto bachTrackExtra = casc.template bachelor_as<DaughterTracks>();
    auto posTrackExtra = casc.template posTrack_as<DaughterTracks>();
    auto negTrackExtra = casc.template negTrack_as<DaughterTracks>();

    // ITS quality flags
    if (bachTrackExtra.itsNCls() < cascSelections.minITSclusters)
      return false;
    if (posTrackExtra.itsNCls() < cascSelections.minITSclusters)
      return false;
    if (negTrackExtra.itsNCls() < cascSelections.minITSclusters)
      return false;

    // TPC quality flags
    if (bachTrackExtra.tpcNClsCrossedRows() < cascSelections.minTPCrows)
      return false;
    if (posTrackExtra.tpcNClsCrossedRows() < cascSelections.minTPCrows)
      return false;
    if (negTrackExtra.tpcNClsCrossedRows() < cascSelections.minTPCrows)
      return false;

    // TPC PID
    if ((cascType == kXiM || cascType == kXiP) && std::fabs(bachTrackExtra.tpcNSigmaPi()) > cascSelections.tpcPidNsigmaCut)
      return false;
    if ((cascType == kOmegaM || cascType == kOmegaP) && std::fabs(bachTrackExtra.tpcNSigmaKa()) > cascSelections.tpcPidNsigmaCut)
      return false;
    if (casc.sign() < 0) { // Xi- or Omega- --> positive/negative daughter = proton/pion
      if (std::fabs(posTrackExtra.tpcNSigmaPr()) > cascSelections.tpcPidNsigmaCut)
        return false;
      if (std::fabs(negTrackExtra.tpcNSigmaPi()) > cascSelections.tpcPidNsigmaCut)
        return false;
    } else { // Xi+ or Omega+ --> positive/negative daughter = pion/proton
      if (std::fabs(posTrackExtra.tpcNSigmaPi()) > cascSelections.tpcPidNsigmaCut)
        return false;
      if (std::fabs(negTrackExtra.tpcNSigmaPr()) > cascSelections.tpcPidNsigmaCut)
        return false;
    }

    // TOF Requirement checks
    if (cascSelections.requireBachHasTOF && !casc.bachelorHasTOF()) {
      return false;
    }
    if (cascSelections.requirePosHasTOF && !casc.positiveHasTOF()) {
      return false;
    }
    if (cascSelections.requireNegHasTOF && !casc.negativeHasTOF()) {
      return false;
    }

    if (cascSelections.requireAtLeastOneHasTOF && !casc.bachelorHasTOF() && !casc.positiveHasTOF() && !casc.negativeHasTOF()) {
      return false;
    }

    //
    // TOF PID in NSigma
    // Bachelor track
    if (casc.bachelorHasTOF()) {
      if ((cascType == kXiM || cascType == kXiP) && std::fabs(casc.tofNSigmaXiPi()) > cascSelections.tofPidNsigmaCutXiPi)
        return false;
      if ((cascType == kOmegaM || cascType == kOmegaP) && std::fabs(casc.tofNSigmaOmKa()) > cascSelections.tofPidNsigmaCutOmKa)
        return false;
    }
    // Positive track
    if (casc.positiveHasTOF()) {
      if (casc.sign() < 0) { // Xi- or Omega- --> positive daughter = proton
        if (cascType == kXiM && std::fabs(casc.tofNSigmaXiLaPr()) > cascSelections.tofPidNsigmaCutLaPr)
          return false;
        if (cascType == kOmegaM && std::fabs(casc.tofNSigmaOmLaPr()) > cascSelections.tofPidNsigmaCutLaPr)
          return false;
      } else { // Xi+ or Omega+ --> positive daughter = pion
        if (cascType == kXiP && std::fabs(casc.tofNSigmaXiLaPi()) > cascSelections.tofPidNsigmaCutLaPi)
          return false;
        if (cascType == kOmegaP && std::fabs(casc.tofNSigmaOmLaPi()) > cascSelections.tofPidNsigmaCutLaPi)
          return false;
      }
    }
    // Negative track
    if (casc.negativeHasTOF()) {
      if (casc.sign() < 0) { // Xi- or Omega- --> negative daughter = pion
        if (cascType == kXiM && std::fabs(casc.tofNSigmaXiLaPr()) > cascSelections.tofPidNsigmaCutLaPi)
          return false;
        if (cascType == kOmegaM && std::fabs(casc.tofNSigmaOmLaPr()) > cascSelections.tofPidNsigmaCutLaPi)
          return false;
      } else { // Xi+ or Omega+ --> negative daughter = proton
        if (cascType == kXiP && std::fabs(casc.tofNSigmaXiLaPi()) > cascSelections.tofPidNsigmaCutLaPr)
          return false;
        if (cascType == kOmegaP && std::fabs(casc.tofNSigmaOmLaPi()) > cascSelections.tofPidNsigmaCutLaPr)
          return false;
      }
    }

    // TRD Requirement checks
    if (cascSelections.requireBachHasTRD && !bachTrackExtra.hasTRD()) {
      return false;
    }
    if (cascSelections.requirePosHasTRD && !posTrackExtra.hasTRD()) {
      return false;
    }
    if (cascSelections.requireNegHasTRD && !negTrackExtra.hasTRD()) {
      return false;
    }

    int bachTRDhits = 0, posTRDhits = 0, negTRDhits = 0;
    for (unsigned int i = 0; i <= 5; i++) {
      if (bachTrackExtra.trdPattern() & (1 << i)) {
        bachTRDhits++;
      }
      if (posTrackExtra.trdPattern() & (1 << i)) {
        posTRDhits++;
      }
      if (negTrackExtra.trdPattern() & (1 << i)) {
        negTRDhits++;
      }
    }

    if (bachTrackExtra.hasTRD() && bachTRDhits < cascSelections.minTRDclusters) {
      return false;
    }
    if (posTrackExtra.hasTRD() && posTRDhits < cascSelections.minTRDclusters) {
      return false;
    }
    if (negTrackExtra.hasTRD() && negTRDhits < cascSelections.minTRDclusters) {
      return false;
    }

    return true;
  }

  template <typename TV0>
  bool checkV0MCAssociation(TV0 v0, int v0Type)
  // precalculate this information so that a check is one mask operation, not many
  {
    if (!v0.isPhysicalPrimary())
      return false;

    bool isPositiveProton = v0.pdgCodePositive() == PDG_t::kProton;
    bool isPositivePion = v0.pdgCodePositive() == PDG_t::kPiPlus || (doTreatPiToMuon && v0.pdgCodePositive() == PDG_t::kMuonPlus);
    bool isNegativeProton = v0.pdgCodeNegative() == PDG_t::kProtonBar;
    bool isNegativePion = v0.pdgCodeNegative() == PDG_t::kPiMinus || (doTreatPiToMuon && v0.pdgCodeNegative() == PDG_t::kMuonMinus);

    if (v0Type == kK0s && v0.pdgCode() == PDG_t::kK0Short && isPositivePion && isNegativePion) {
      return true;
    }
    if (v0Type == kLambda && v0.pdgCode() == PDG_t::kLambda0 && isPositiveProton && isNegativePion) {
      return true;
    }
    if (v0Type == kAntiLambda && v0.pdgCode() == PDG_t::kLambda0Bar && isPositivePion && isNegativeProton) {
      return true;
    }
    return false;
  }

  template <typename TCascade>
  bool checkCascadeMCAssociation(TCascade casc, int cascType)
  // precalculate this information so that a check is one mask operation, not many
  {
    if (!casc.isPhysicalPrimary())
      return false;

    bool isBachelorPionPlus = casc.pdgCodeBachelor() == PDG_t::kPiPlus || (doTreatPiToMuon && casc.pdgCodeBachelor() == PDG_t::kMuonPlus);
    bool isBachelorKaonPlus = casc.pdgCodeBachelor() == PDG_t::kKPlus;
    bool isBachelorPionMinus = casc.pdgCodeBachelor() == PDG_t::kPiMinus || (doTreatPiToMuon && casc.pdgCodeBachelor() == PDG_t::kMuonMinus);
    bool isBachelorKaonMinus = casc.pdgCodeBachelor() == PDG_t::kKMinus;
    bool isPositiveProton = casc.pdgCodePositive() == PDG_t::kProton;
    bool isPositivePion = casc.pdgCodePositive() == PDG_t::kPiPlus || (doTreatPiToMuon && casc.pdgCodePositive() == PDG_t::kMuonPlus);
    bool isNegativeProton = casc.pdgCodeNegative() == PDG_t::kProtonBar;
    bool isNegativePion = casc.pdgCodeNegative() == PDG_t::kPiMinus || (doTreatPiToMuon && casc.pdgCodeNegative() == PDG_t::kMuonMinus);

    if (cascType == kXiM && casc.pdgCode() != PDG_t::kXiMinus || isPositiveProton || isNegativePion || isBachelorPionMinus) {
      return true;
    }
    if (cascType == kXiP && casc.pdgCode() != PDG_t::kXiPlusBar || isPositivePion || isNegativeProton || isBachelorPionPlus) {
      return true;
    }
    if (cascType == kOmegaM && casc.pdgCode() != PDG_t::kOmegaMinus || isPositiveProton || isNegativePion || isBachelorKaonMinus) {
      return true;
    }
    if (cascType == kOmegaP && casc.pdgCode() != PDG_t::kOmegaPlusBar || isPositivePion || isNegativeProton || isBachelorKaonPlus) {
      return true;
    }
    return false;
  }

  ////////////////////////////////////////////
  /////////// QA - Reconstructed /////////////
  ////////////////////////////////////////////

  void processReconstructed(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults>::iterator const& collision, soa::Join<aod::V0Datas, aod::V0TOFPIDs, aod::V0TOFNSigmas> const& fullV0s, soa::Join<aod::CascDatas, aod::CascTOFPIDs, aod::CascTOFNSigmas> const& fullCascades, DaughterTracks const&, aod::BCsWithTimestamps const&)
  {
    histos_event.fill(HIST("hEventCounter"), 0.5);
    if (!isEventAccepted(collision, true)) {
      return;
    }
    histos_event.fill(HIST("hEventCounter"), 1.5);

    for (auto const& v0 : fullV0s) {
      if (std::abs(v0.negativeeta()) > v0Selections.daughterEtaCut ||
          std::abs(v0.positiveeta()) > v0Selections.daughterEtaCut)
        continue; // remove acceptance that's badly reproduced by MC / superfluous in future

      if (v0Selections.v0TypeSelection > -1 && v0.v0Type() != v0Selections.v0TypeSelection)
        continue; // skip V0s that are not standard

      // fillV0s<DaughterTracks>(v0, collision.posX(), collision.posY(), collision.posZ());

      auto posdau = v0.posTrack_as<DaughterTracks>();
      auto negdau = v0.negTrack_as<DaughterTracks>();

      if (v0.negativeeta() < 0. && v0.positiveeta() < 0.) {
        dauEtaFlag = -1;
      } else if (v0.negativeeta() >= 0. && v0.positiveeta() >= 0.) {
        dauEtaFlag = 1;
      } else {
        dauEtaFlag = 0;
      }

      histos_V0.fill(HIST("CosPA"), v0.v0cosPA());
      histos_V0.fill(HIST("Radius"), v0.v0radius());
      histos_V0.fill(HIST("DCANegToPV"), v0.dcanegtopv());
      histos_V0.fill(HIST("DCAPosToPV"), v0.dcapostopv());
      histos_V0.fill(HIST("DCAV0Daughters"), v0.dcaV0daughters());

      float decayLength = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * RecoDecay::sqrtSumOfSquares(v0.px(), v0.py(), v0.pz());
      histos_V0.fill(HIST("DecayLength"), decayLength);

      float CtauLambda = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassLambda0;
      float CtauK0s = v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * o2::constants::physics::MassK0Short;

      // K0Short
      if (isV0Accepted(v0, collision, v0.rapidity(0), kK0s)) {
        histos_V0.fill(HIST("CosPAK0s"), v0.v0cosPA());
        histos_V0.fill(HIST("RadiusK0s"), v0.v0radius());
        histos_V0.fill(HIST("DCANegToPVK0s"), v0.dcanegtopv());
        histos_V0.fill(HIST("DCAPosToPVK0s"), v0.dcapostopv());
        histos_V0.fill(HIST("DCAV0DaughtersK0s"), v0.dcaV0daughters());
        histos_V0.fill(HIST("DecayLengthK0s"), decayLength);
        histos_V0.fill(HIST("LifetimeK0s"), CtauK0s);
        histos_V0.fill(HIST("InvMassK0s"), v0.pt(), v0.mK0Short(), dauEtaFlag);
        histos_V0.fill(HIST("DCAV0ToPVK0s"), v0.dcav0topv());
        histos_V0.fill(HIST("TPCPIDPosPionFromK0s"), v0.pt(), posdau.tpcNSigmaPi());
        histos_V0.fill(HIST("TPCPIDNegPionFromK0s"), v0.pt(), negdau.tpcNSigmaPi());
        if (doextraanalysis) {
          histos_V0.fill(HIST("InvMassK0s_Radius"), v0.v0radius(), v0.mK0Short());
          histos_V0.fill(HIST("InvMassK0s_PtRadius"), v0.pt(), v0.v0radius(), v0.mK0Short());
          histos_V0.fill(HIST("InvMassK0s_Lifetime"), CtauK0s, v0.mK0Short());
          histos_V0.fill(HIST("InvMassK0s_EtaDaughters"), posdau.eta(), negdau.eta(), v0.mK0Short());
          histos_V0.fill(HIST("InvMassK0s_PhiDaughters"), posdau.phi(), negdau.phi(), v0.mK0Short());
          histos_V0.fill(HIST("InvMassK0s_ITSMapDaughters"), posdau.itsNCls(), negdau.itsNCls(), v0.mK0Short());
          histos_V0.fill(HIST("InvMassK0sVsPtVsPA"), v0.pt(), std::acos(v0.v0cosPA()), v0.mK0Short());
        }
      }

      // Lambda
      if (isV0Accepted(v0, collision, v0.rapidity(1), kLambda)) {
        histos_V0.fill(HIST("CosPALambda"), v0.v0cosPA());
        histos_V0.fill(HIST("RadiusLambda"), v0.v0radius());
        histos_V0.fill(HIST("DCANegToPVLambda"), v0.dcanegtopv());
        histos_V0.fill(HIST("DCAPosToPVLambda"), v0.dcapostopv());
        histos_V0.fill(HIST("DCAV0DaughtersLambda"), v0.dcaV0daughters());
        histos_V0.fill(HIST("DecayLengthLambda"), decayLength);
        histos_V0.fill(HIST("LifetimeLambda"), CtauLambda);
        histos_V0.fill(HIST("InvMassLambda"), v0.pt(), v0.mLambda(), dauEtaFlag);
        histos_V0.fill(HIST("DCAV0ToPVLambda"), v0.dcav0topv());
        histos_V0.fill(HIST("TPCPIDPionFromLambda"), v0.pt(), negdau.tpcNSigmaPi());
        histos_V0.fill(HIST("TPCPIDProtonFromLambda"), v0.pt(), posdau.tpcNSigmaPr());
        if (doextraanalysis) {
          histos_V0.fill(HIST("InvMassLambda_Radius"), v0.v0radius(), v0.mLambda());
          histos_V0.fill(HIST("InvMassLambda_PtRadius"), v0.pt(), v0.v0radius(), v0.mLambda());
          histos_V0.fill(HIST("InvMassLambda_Lifetime"), CtauLambda, v0.mLambda());
          histos_V0.fill(HIST("InvMassLambda_EtaDaughters"), posdau.eta(), negdau.eta(), v0.mLambda());
          histos_V0.fill(HIST("InvMassLambda_PhiDaughters"), posdau.phi(), negdau.phi(), v0.mLambda());
          histos_V0.fill(HIST("InvMassLambda_ITSMapDaughters"), posdau.itsNCls(), negdau.itsNCls(), v0.mLambda());
          histos_V0.fill(HIST("InvMassLambdaVsPtVsPA"), v0.pt(), std::acos(v0.v0cosPA()), v0.mLambda());
        }
      }

      // AntiLambda
      if (isV0Accepted(v0, collision, v0.rapidity(2), kAntiLambda)) {
        histos_V0.fill(HIST("CosPAAntiLambda"), v0.v0cosPA());
        histos_V0.fill(HIST("RadiusAntiLambda"), v0.v0radius());
        histos_V0.fill(HIST("DCANegToPVAntiLambda"), v0.dcanegtopv());
        histos_V0.fill(HIST("DCAPosToPVAntiLambda"), v0.dcapostopv());
        histos_V0.fill(HIST("DCAV0DaughtersAntiLambda"), v0.dcaV0daughters());
        histos_V0.fill(HIST("DecayLengthAntiLambda"), decayLength);
        histos_V0.fill(HIST("LifetimeAntiLambda"), CtauLambda);
        histos_V0.fill(HIST("InvMassAntiLambda"), v0.pt(), v0.mAntiLambda(), dauEtaFlag);
        histos_V0.fill(HIST("DCAV0ToPVAntiLambda"), v0.dcav0topv());
        histos_V0.fill(HIST("TPCPIDPionFromAntiLambda"), v0.pt(), posdau.tpcNSigmaPi());
        histos_V0.fill(HIST("TPCPIDProtonFromAntiLambda"), v0.pt(), negdau.tpcNSigmaPr());
        if (doextraanalysis) {
          histos_V0.fill(HIST("InvMassAntiLambda_Radius"), v0.v0radius(), v0.mAntiLambda());
          histos_V0.fill(HIST("InvMassAntiLambda_PtRadius"), v0.pt(), v0.v0radius(), v0.mAntiLambda());
          histos_V0.fill(HIST("InvMassAntiLambda_Lifetime"), CtauLambda, v0.mAntiLambda());
          histos_V0.fill(HIST("InvMassAntiLambda_EtaDaughters"), posdau.eta(), negdau.eta(), v0.mAntiLambda());
          histos_V0.fill(HIST("InvMassAntiLambda_PhiDaughters"), posdau.phi(), negdau.phi(), v0.mAntiLambda());
          histos_V0.fill(HIST("InvMassAntiLambda_ITSMapDaughters"), posdau.itsNCls(), negdau.itsNCls(), v0.mAntiLambda());
          histos_V0.fill(HIST("InvMassAntiLambdaVsPtVsPA"), v0.pt(), std::acos(v0.v0cosPA()), v0.mAntiLambda());
        }
      }
    }

    for (auto const& casc : fullCascades) {
      if (std::abs(casc.negativeeta()) > cascSelections.daughterEtaCut ||
          std::abs(casc.positiveeta()) > cascSelections.daughterEtaCut ||
          std::abs(casc.bacheloreta()) > cascSelections.daughterEtaCut)
        continue; // remove acceptance that's badly reproduced by MC / superfluous in future

      histos_Casc.fill(HIST("CascCosPA"), casc.casccosPA(collision.posX(), collision.posY(), collision.posZ()), casc.sign());
      histos_Casc.fill(HIST("V0CosPA"), casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ()), casc.sign());

      double v0cospatoxi = RecoDecay::cpa(std::array{casc.x(), casc.y(), casc.z()}, array{casc.xlambda(), casc.ylambda(), casc.zlambda()}, std::array{casc.pxpos() + casc.pxneg(), casc.pypos() + casc.pyneg(), casc.pzpos() + casc.pzneg()});

      histos_Casc.fill(HIST("V0CosPAToXi"), v0cospatoxi, casc.sign());
      histos_Casc.fill(HIST("CascRadius"), casc.cascradius(), casc.sign());
      histos_Casc.fill(HIST("V0Radius"), casc.v0radius(), casc.sign());
      histos_Casc.fill(HIST("CascRapidityXi"), casc.yXi(), casc.sign());
      histos_Casc.fill(HIST("CascRapidityOmega"), casc.yOmega(), casc.sign());

      float cascDecayLength = std::sqrt(std::pow(casc.x() - collision.posX(), 2) + std::pow(casc.y() - collision.posY(), 2) + std::pow(casc.z() - collision.posZ(), 2));
      histos_Casc.fill(HIST("CascDecayLength"), cascDecayLength, casc.sign());

      float cascTotalMomentum = RecoDecay::sqrtSumOfSquares(casc.px(), casc.py(), casc.pz());
      float CtauXi = cascDecayLength / (cascTotalMomentum + 1E-10) * o2::constants::physics::MassXi0; // see O2Physics/Common/Core/MC.h for codes and names accepted
      float CtauOmega = cascDecayLength / (cascTotalMomentum + 1E-10) * o2::constants::physics::MassOmegaMinus;

      float v0TotalMomentum = RecoDecay::sqrtSumOfSquares(casc.pxpos() + casc.pxneg(), casc.pypos() + casc.pyneg(), casc.pzpos() + casc.pzneg());
      float v0DecayLength = std::sqrt(std::pow(casc.xlambda() - casc.x(), 2) + std::pow(casc.ylambda() - casc.y(), 2) + std::pow(casc.zlambda() - casc.z(), 2));
      float CtauV0 = v0DecayLength / (v0TotalMomentum + 1E-10) * o2::constants::physics::MassLambda0;

      histos_Casc.fill(HIST("CascLifetimeXi"), CtauXi, casc.sign());
      histos_Casc.fill(HIST("CascLifetimeOmega"), CtauOmega, casc.sign());
      histos_Casc.fill(HIST("V0Lifetime"), CtauV0, casc.sign());
      histos_Casc.fill(HIST("CascPt"), casc.pt(), casc.sign());
      histos_Casc.fill(HIST("DcaV0Daughters"), casc.dcaV0daughters(), casc.sign());
      histos_Casc.fill(HIST("DcaCascDaughters"), casc.dcacascdaughters(), casc.sign());
      histos_Casc.fill(HIST("DcaV0ToPV"), casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()), casc.sign());
      histos_Casc.fill(HIST("DcaBachToPV"), casc.dcabachtopv(), casc.sign());
      histos_Casc.fill(HIST("DcaPosToPV"), casc.dcapostopv(), casc.sign());
      histos_Casc.fill(HIST("DcaNegToPV"), casc.dcanegtopv(), casc.sign());
      histos_Casc.fill(HIST("InvMassLambda"), casc.mLambda(), casc.sign());

      if (isCascadeSelected(casc, collision, casc.rapidity(0), kXiM)) {
        histos_Casc.fill(HIST("InvMassXiMinus"), casc.pt(), casc.mXi(), casc.eta());
        histos_Casc.fill(HIST("InvMassXiMinus_Radius"), casc.cascradius(), casc.mXi());
      }
      if (isCascadeSelected(casc, collision, casc.rapidity(0), kXiP)) {
        histos_Casc.fill(HIST("InvMassXiPlus"), casc.pt(), casc.mXi(), casc.eta());
        histos_Casc.fill(HIST("InvMassXiPlus_Radius"), casc.cascradius(), casc.mXi());
      }
      if (isCascadeSelected(casc, collision, casc.rapidity(2), kOmegaM)) {
        histos_Casc.fill(HIST("InvMassOmegaMinus"), casc.pt(), casc.mOmega(), casc.eta());
      }
      if (isCascadeSelected(casc, collision, casc.rapidity(2), kOmegaP)) {
        histos_Casc.fill(HIST("InvMassOmegaPlus"), casc.pt(), casc.mOmega(), casc.eta());
      }
    }
  }

  ////////////////////////////////
  ////////// QA - MC /////////////
  ////////////////////////////////

  void processMonteCarlo(soa::Join<aod::Collisions, aod::EvSels, aod::PVMults, aod::McCollisionLabels>::iterator const& collision, soa::Join<aod::McCollisions, aod::MultsExtraMC> const&, soa::Join<aod::V0Datas, aod::V0TOFPIDs, aod::V0TOFNSigmas, aod::V0CoreMCLabels> const& fullV0s, soa::Join<aod::V0MCDatas, aod::V0MCCollRefs> const&, soa::Join<aod::CascDatas, aod::CascTOFPIDs, aod::CascTOFNSigmas, aod::CascCoreMCLabels> const& fullCascades, soa::Join<aod::CascMCDatas, aod::CascMCCollRefs> const&, DaughterTracks const&, aod::BCsWithTimestamps const&)
  {
    if (!isEventAccepted(collision, false)) {
      return;
    }
    if (!collision.has_mcCollision())
      return;
    auto mcCollision = collision.mcCollision_as<soa::Join<aod::McCollisions, aod::MultsExtraMC>>();
    // Apply selections on MC collisions
    if (eventSelections.applyZVtxSelOnMCPV && std::abs(mcCollision.posZ()) > eventSelections.maxZVtxPosition) {
      return;
    }
    if (eventSelections.requireINEL0 && mcCollision.multMCNParticlesEta10() < 1) {
      return;
    }
    if (eventSelections.requireINEL1 && mcCollision.multMCNParticlesEta10() < 2) {
      return;
    }

    for (auto const& v0 : fullV0s) {
      if (std::abs(v0.negativeeta()) > v0Selections.daughterEtaCut ||
          std::abs(v0.positiveeta()) > v0Selections.daughterEtaCut)
        continue; // remove acceptance that's badly reproduced by MC / superfluous in future

      if (v0Selections.v0TypeSelection > -1 && v0.v0Type() != v0Selections.v0TypeSelection)
        continue; // skip V0s that are not standard

      if (!v0.has_v0MCCore())
        continue;

      auto v0MC = v0.template v0MCCore_as<soa::Join<aod::V0MCCores, aod::V0MCCollRefs>>();

      // K0Short
      if (isV0Accepted(v0, collision, v0MC.rapidityMC(0), kK0s) && checkV0MCAssociation(v0MC, kK0s)) {
        histos_V0.fill(HIST("InvMassK0sTrue"), v0MC.ptMC(), v0.v0radius(), v0.mK0Short());
      }
      // Lambda
      if (isV0Accepted(v0, collision, v0MC.rapidityMC(1), kLambda) && checkV0MCAssociation(v0MC, kLambda)) {
        histos_V0.fill(HIST("InvMassLambdaTrue"), v0MC.ptMC(), v0.v0radius(), v0.mLambda());
      }
      // AntiLambda
      if (isV0Accepted(v0, collision, v0MC.rapidityMC(2), kAntiLambda) && checkV0MCAssociation(v0MC, kAntiLambda)) {
        histos_V0.fill(HIST("InvMassAntiLambdaTrue"), v0MC.ptMC(), v0.v0radius(), v0.mAntiLambda());
      }
    }

    for (auto const& casc : fullCascades) {
      if (std::abs(casc.negativeeta()) > cascSelections.daughterEtaCut ||
          std::abs(casc.positiveeta()) > cascSelections.daughterEtaCut ||
          std::abs(casc.bacheloreta()) > cascSelections.daughterEtaCut)
        continue; // remove acceptance that's badly reproduced by MC / superfluous in future

      if (!casc.has_cascMCCore())
        continue;

      auto cascMC = casc.template cascMCCore_as<soa::Join<aod::CascMCDatas, aod::CascMCCollRefs>>();

      histos_Casc.fill(HIST("QA_CascadeCandidates"), 0.5);

      if (isCascadeSelected(casc, collision, cascMC.rapidityMC(0), kXiM)) {
        histos_Casc.fill(HIST("QA_CascCandidates"), 1.5);
        if (checkCascadeMCAssociation(cascMC, kXiM)) {
          histos_Casc.fill(HIST("QA_CascadeCandidates"), 2.5);
          histos_Casc.fill(HIST("InvMassXiMinusTrue"), cascMC.ptMC(), casc.cascradius(), casc.mXi());
        }
      }
      if (isCascadeSelected(casc, collision, cascMC.rapidityMC(0), kXiP)) {
        histos_Casc.fill(HIST("QA_CascadeCandidates"), 3.5);
        if (checkCascadeMCAssociation(cascMC, kXiP)) {
          histos_Casc.fill(HIST("QA_CascadeCandidates"), 4.5);
          histos_Casc.fill(HIST("InvMassXiPlusTrue"), cascMC.ptMC(), casc.cascradius(), casc.mXi());
        }
      }
      if (isCascadeSelected(casc, collision, cascMC.rapidityMC(2), kOmegaM)) {
        histos_Casc.fill(HIST("QA_CascCandidates"), 5.5);
        if (checkCascadeMCAssociation(cascMC, kOmegaM)) {
          histos_Casc.fill(HIST("QA_CascadeCandidates"), 6.5);
          histos_Casc.fill(HIST("InvMassOmegaMinusTrue"), cascMC.ptMC(), casc.cascradius(), casc.mOmega());
        }
      }
      if (isCascadeSelected(casc, collision, cascMC.rapidityMC(2), kOmegaP)) {
        histos_Casc.fill(HIST("QA_CascadeCandidates"), 7.5);
        if (checkCascadeMCAssociation(cascMC, kOmegaP)) {
          histos_Casc.fill(HIST("QA_CascadeCandidates"), 8.5);
          histos_Casc.fill(HIST("InvMassOmegaPlusTrue"), cascMC.ptMC(), casc.cascradius(), casc.mOmega());
        }
      }
    }
  }

  ///////////////////////////////////////
  ////////// Collision QA - MC //////////
  ///////////////////////////////////////

  void processGenerated(soa::Join<aod::McCollisions, aod::MultsExtraMC>::iterator const& mcCollision, aod::McParticles const& mcParticles, soa::SmallGroups<o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels, aod::PVMults>> const& collisions)
  {
    // Apply selections on MC collisions
    if (eventSelections.applyZVtxSelOnMCPV && std::abs(mcCollision.posZ()) > eventSelections.maxZVtxPosition) {
      return;
    }
    if (eventSelections.requireINEL0 && mcCollision.multMCNParticlesEta10() < 1) {
      return;
    }
    if (eventSelections.requireINEL1 && mcCollision.multMCNParticlesEta10() < 2) {
      return;
    }

    histos_event.fill(HIST("hEventCounterMC"), 0.5);

    // Check if there is at least one of the reconstructed collisions associated to this MC collision
    // If so, we consider it
    bool atLeastOne = false;
    for (const auto& collision : collisions) {
      if (!isEventAccepted(collision, false)) {
        continue;
      }
      atLeastOne = true;
    }

    if (!atLeastOne) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }

    histos_event.fill(HIST("hEventCounterMC"), 1.5);

    for (auto const& mcparticle : mcParticles) {

      if (!mcparticle.has_daughters()) {
        continue;
      }

      double vx = 0;
      double vy = 0;
      for (auto const& mcparticleDaughter0 : mcparticle.daughters_as<aod::McParticles>()) {
        vx = mcparticleDaughter0.vx();
        vy = mcparticleDaughter0.vy();
        if (vx != 0 && vy != 0)
          break;
      }
      double R_Decay = std::sqrt(vx * vx + vy * vy);

      if (mcparticle.isPhysicalPrimary() && std::abs(mcparticle.y()) < v0Selections.rapidityCut) {
        if (mcparticle.pdgCode() == PDG_t::kK0Short)
          histos_event.fill(HIST("GeneratedV0s"), 0.5, mcparticle.pt(), R_Decay); // K0s
        if (mcparticle.pdgCode() == PDG_t::kLambda0)
          histos_event.fill(HIST("GeneratedV0s"), 1.5, mcparticle.pt(), R_Decay); // Lambda
        if (mcparticle.pdgCode() == PDG_t::kLambda0Bar)
          histos_event.fill(HIST("GeneratedV0s"), 2.5, mcparticle.pt(), R_Decay); // AntiLambda
      }
      if (mcparticle.isPhysicalPrimary() && std::abs(mcparticle.y()) < cascSelections.rapidityCut) {
        if (mcparticle.pdgCode() == PDG_t::kXiMinus)
          histos_event.fill(HIST("GeneratedCascades"), 0.5, mcparticle.pt(), R_Decay); // Xi-
        if (mcparticle.pdgCode() == PDG_t::kXiPlusBar)
          histos_event.fill(HIST("GeneratedCascades"), 1.5, mcparticle.pt(), R_Decay); // Xi+
        if (mcparticle.pdgCode() == PDG_t::kOmegaMinus)
          histos_event.fill(HIST("GeneratedCascades"), 2.5, mcparticle.pt(), R_Decay); // Omega-
        if (mcparticle.pdgCode() == PDG_t::kOmegaPlusBar)
          histos_event.fill(HIST("GeneratedCascades"), 3.5, mcparticle.pt(), R_Decay); // Omega+

        // if (!IsParticleFromOutOfBunchPileupCollision){fill the 1.5, 3.5 etc}   AliPhysics analysis
      }
    }
  }

  PROCESS_SWITCH(v0cascadesQA, processReconstructed, "Process reconstructed event and V0s+cascades in data", true);
  PROCESS_SWITCH(v0cascadesQA, processMonteCarlo, "Process reconstructed event and V0s+cascades in MC", false);
  PROCESS_SWITCH(v0cascadesQA, processGenerated, "Process MC level event and V0s+cascades in MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<v0cascadesQA>(cfgc)};
}
