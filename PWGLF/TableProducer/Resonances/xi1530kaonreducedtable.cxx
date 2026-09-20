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

/// \file xi1530kaonreducedtable.cxx
/// \brief Producer of compact Xi(1530)^0--K reduced tables for femtoscopy.
/// \author Prottay Das, prottay.das@cern.ch

#include "PWGLF/DataModel/LFReducedXi1530KaonTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h>
#include <TH1.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;

namespace
{
constexpr float kTinyMomentum = 1.e-6f;

inline uint8_t toU8(int value)
{
  return static_cast<uint8_t>(std::clamp(value, 0, 255));
}

inline void setTrackSelBit(uint64_t& word, o2::aod::redxistark::XiStarTrackSelBit bit, bool pass)
{
  if (pass) {
    word |= o2::aod::redxistark::xiStarTrackSelMask(bit);
  }
}

inline void setTopoSelBit(uint64_t& word, o2::aod::redxistark::XiStarTopoSelBit bit, bool pass)
{
  if (pass) {
    word |= o2::aod::redxistark::xiStarTopoSelMask(bit);
  }
}

inline void setKaonSelBit(uint8_t& word, o2::aod::redxistarkaon::KaonSelBit bit, bool pass)
{
  if (pass) {
    word |= o2::aod::redxistarkaon::kaonSelMask(bit);
  }
}

} // namespace

struct xi1530kaonreducedtable {
  // CCDB is used only to store the nominal magnetic field in the reduced event
  // row. This keeps the downstream CPR QA self-contained without storing the
  // run number or timestamp.
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdbURL{"ccdbURL", "http://alice-ccdb.cern.ch", "CCDB URL for the magnetic field"};

  // ---------------------------------------------------------------------------
  // Output tables
  // ---------------------------------------------------------------------------
  Produces<aod::RedXiStarKEvents> redEvents;
  Produces<aod::XiStarCandidates> xiStarCandidates;
  Produces<aod::KaonCandidates> kaonCandidates;

  // ---------------------------------------------------------------------------
  // Input table types
  // ---------------------------------------------------------------------------
  using EventCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::Mults, aod::CentFT0Ms>;

  // One joined track type is used both for primary tracks and for the daughters
  // obtained through CascDatas::{posTrack,negTrack,bachelor}_as<TracksFull>().
  using TracksFull = soa::Join<aod::Tracks,
                               aod::TracksExtra,
                               aod::TracksDCA,
                               aod::TrackSelection,
                               aod::pidTPCFullPi,
                               aod::pidTOFFullPi,
                               aod::pidTPCFullKa,
                               aod::pidTOFFullKa,
                               aod::pidTPCFullPr>;

  // ---------------------------------------------------------------------------
  // Event selections -- copied from the Xi(1530) Run-3 analysis philosophy.
  // ---------------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "event";
    Configurable<float> zVtxMax{"zVtxMax", 10.f, "Maximum |z_PV| (cm)"};
    Configurable<bool> requireSel8{"requireSel8", true, "Require collision.sel8()"};
    Configurable<bool> requireINELgt0{"requireINELgt0", true, "Require multNTracksPVeta1 > 0"};
  } eventCuts;

  // ---------------------------------------------------------------------------
  // Prompt / first pion from Xi(1530) -> Xi + pi.
  // ---------------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "piFirst";

    Configurable<float> etaMax{"etaMax", 0.8f, "Maximum |eta|"};
    Configurable<float> ptMin{"ptMin", 0.15f, "Minimum pT (GeV/c)"};
    Configurable<bool> requireGlobalTrack{"requireGlobalTrack", true, "Require isGlobalTrack"};
    Configurable<bool> requirePVContributor{"requirePVContributor", true, "Require isPVContributor"};
    Configurable<int> minITSInnerBarrelClusters{"minITSInnerBarrelClusters", 1, "Minimum ITS inner-barrel clusters"};
    Configurable<int> minITSClusters{"minITSClusters", 3, "Minimum ITS clusters"};

    // DCAxy: configurable; nominal/default form is pT dependent.
    Configurable<bool> dcaXYUsePtDependent{"dcaXYUsePtDependent", true, "Use pT-dependent default DCAxy cut"};
    Configurable<float> dcaXYp0{"dcaXYp0", 0.004f, "DCAxy(pT) p0 (cm)"};
    Configurable<float> dcaXYp1{"dcaXYp1", 0.013f, "DCAxy(pT) p1 (cm GeV/c)"};
    Configurable<float> dcaXYMaxConstant{"dcaXYMaxConstant", 0.1f, "Constant |DCAxy| max if pT-dependent mode is disabled (cm)"};

    // DCAz: independently configurable; nominal/default form is also pT dependent.
    // Two additional alternatives are retained in the track-selection bitmap.
    Configurable<bool> dcaZUsePtDependent{"dcaZUsePtDependent", true, "Use pT-dependent default DCAz cut"};
    Configurable<float> dcaZp0{"dcaZp0", 0.004f, "DCAz(pT) p0 (cm)"};
    Configurable<float> dcaZp1{"dcaZp1", 0.013f, "DCAz(pT) p1 (cm GeV/c)"};
    Configurable<float> dcaZMaxConstant{"dcaZMaxConstant", 0.1f, "Constant default |DCAz| max if pT-dependent mode is disabled (cm)"};
    Configurable<float> dcaZLoose{"dcaZLoose", 1.0f, "Alternative loose |DCAz| maximum (cm)"};
    Configurable<float> dcaZTightVariation{"dcaZTightVariation", 0.1f, "Alternative |DCAz| maximum (cm)"};

    // TPC rows: loose/default/tight from the Xi(1530) systematic table.
    Configurable<int> tpcRowsLoose{"tpcRowsLoose", 70, "Loose minimum crossed rows (strict >)"};
    Configurable<int> tpcRowsDefault{"tpcRowsDefault", 80, "Default minimum crossed rows (strict >)"};
    Configurable<int> tpcRowsTight{"tpcRowsTight", 90, "Tight minimum crossed rows (strict >)"};
  } piFirstCuts;

  // ---------------------------------------------------------------------------
  // PID working points.
  // V0 daughters (Lambda pion/proton): TPC only.
  // Prompt Xi* pion: retain BOTH pure TPC-only and hybrid TPC/TOF choices.
  // Xi bachelor pion: retain BOTH pure TPC-only and hybrid TPC/TOF choices.
  // Hybrid means: if no TOF, TPC only; if TOF exists, use a circular
  // sqrt(nSigmaTPC^2+nSigmaTOF^2) cut.
  // ---------------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "v0Pid";
    Configurable<float> nSigma3{"nSigma3", 3.f, "V0-daughter TPC PID 3-sigma WP"};
    Configurable<float> nSigma4{"nSigma4", 4.f, "V0-daughter TPC PID 4-sigma WP (dedicated Lambda nominal)"};
    Configurable<float> nSigma4p8{"nSigma4p8", 4.8f, "V0 pion TPC PID pp-Xi WP"};
    Configurable<float> nSigma5{"nSigma5", 5.f, "V0-daughter TPC PID 5-sigma WP / pp-Xi proton WP"};
    Configurable<float> nSigma6{"nSigma6", 6.f, "V0-daughter TPC PID 6-sigma loose envelope"};
  } v0PidCuts;

  struct : ConfigurableGroup {
    std::string prefix = "otherTrackPid";
    Configurable<float> nSigma3{"nSigma3", 3.f, "Prompt/bachelor pion PID 3-sigma WP (TPC-only or hybrid)"};
    Configurable<float> nSigma4{"nSigma4", 4.f, "Prompt/bachelor pion PID 4-sigma WP (TPC-only or hybrid)"};
    Configurable<float> nSigma5{"nSigma5", 5.f, "Prompt/bachelor pion PID 5-sigma WP (TPC-only or hybrid)"};
    Configurable<float> nSigma6{"nSigma6", 6.f, "Prompt/bachelor pion PID 6-sigma loose envelope (TPC-only or hybrid)"};
    Configurable<float> bachelorTPC4p8{"bachelorTPC4p8", 4.8f, "TPC-only Xi-bachelor pion pp-Xi WP"};
  } otherTrackPidCuts;

  // Second prompt-pion PID recipe:
  //   pT < ptThreshold  -> TPC-only PID, even if the track has TOF
  //   pT >= ptThreshold -> circular sqrt(TPC^2 + TOF^2) PID
  // By default TOF is required above the threshold.  This can be disabled.
  // TPC and circular limits are independently configurable for each WP.
  struct : ConfigurableGroup {
    std::string prefix = "piFirstThresholdPid";
    Configurable<float> ptThreshold{"ptThreshold", 0.5f, "Prompt-pion pT threshold above which circular TPC+TOF PID is applied (GeV/c)"};
    Configurable<bool> requireTOFAboveThreshold{"requireTOFAboveThreshold", true, "Require valid TOF for prompt pion at/above ptThreshold"};

    Configurable<float> tpcNSigma3{"tpcNSigma3", 3.f, "TPC-only limit for threshold PID 3-sigma WP"};
    Configurable<float> tpcNSigma4{"tpcNSigma4", 4.f, "TPC-only limit for threshold PID 4-sigma WP"};
    Configurable<float> tpcNSigma5{"tpcNSigma5", 5.f, "TPC-only limit for threshold PID 5-sigma WP"};
    Configurable<float> tpcNSigma6{"tpcNSigma6", 6.f, "TPC-only limit for threshold PID 6-sigma loose envelope"};

    Configurable<float> circularNSigma3{"circularNSigma3", 3.f, "Circular TPC+TOF limit for threshold PID 3-sigma WP"};
    Configurable<float> circularNSigma4{"circularNSigma4", 4.f, "Circular TPC+TOF limit for threshold PID 4-sigma WP"};
    Configurable<float> circularNSigma5{"circularNSigma5", 5.f, "Circular TPC+TOF limit for threshold PID 5-sigma WP"};
    Configurable<float> circularNSigma6{"circularNSigma6", 6.f, "Circular TPC+TOF limit for threshold PID 6-sigma loose envelope"};
  } piFirstThresholdPidCuts;

  // ---------------------------------------------------------------------------
  // Cascade daughter quality and fixed cascade acceptance/lifetime cuts.
  // Exact source-defined topology working points are in the next group.
  // ---------------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "cascadeBase";

    // Fixed producer-level acceptance / quality cuts.
    Configurable<float> daughterEtaMax{"daughterEtaMax", 0.8f, "Maximum |eta| of cascade daughters (fixed producer cut)"};
    Configurable<int> daughterTPCrowsMin{"daughterTPCrowsMin", 50, "Broad minimum TPC crossed rows for all cascade daughters (strict >)"};
    Configurable<float> lambdaPionPtMin{"lambdaPionPtMin", 0.2f, "Minimum Lambda-daughter pion pT (GeV/c), fixed producer cut"};
    Configurable<float> lambdaProtonPtMin{"lambdaProtonPtMin", 0.2f, "Minimum Lambda-daughter proton pT (GeV/c), fixed producer cut"};
    Configurable<float> lambdaRapidityMax{"lambdaRapidityMax", 0.5f, "Maximum |y_Lambda|, fixed producer cut"};
    Configurable<float> cascadeEtaMax{"cascadeEtaMax", 0.8f, "Maximum |eta| of Xi candidate (fixed producer cut)"};
    Configurable<float> cascadePtMin{"cascadePtMin", 0.15f, "Minimum Xi-candidate pT (GeV/c), fixed producer cut"};
    Configurable<float> maxV0Radius{"maxV0Radius", 30.f, "Maximum V0/Lambda radius (cm), fixed producer cut"};
    Configurable<float> maxCascRadius{"maxCascRadius", 200.f, "Maximum cascade radius (cm), fixed producer cut"};
  } cascadeBaseCuts;

  // ---------------------------------------------------------------------------
  // Exact topology working points retained in the topology-selection bitmap.
  // The producer applies only the loosest UNION envelope of these values.
  // ---------------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "topology";

    // V0 daughter DCA to PV (minimum, cm).
    Configurable<float> v0PionDcaPV005{"v0PionDcaPV005", 0.05f, "Min pion DCA to PV: 0.05 cm"};
    Configurable<float> v0PionDcaPV006{"v0PionDcaPV006", 0.06f, "Min pion DCA to PV: 0.06 cm (Xi* nominal)"};
    Configurable<float> v0PionDcaPV010{"v0PionDcaPV010", 0.10f, "Min pion DCA to PV: 0.10 cm"};
    Configurable<float> v0PionDcaPV020{"v0PionDcaPV020", 0.20f, "Min pion DCA to PV: 0.20 cm (dedicated Lambda)"};

    Configurable<float> v0ProtonDcaPV005{"v0ProtonDcaPV005", 0.05f, "Min proton DCA to PV: 0.05 cm"};
    Configurable<float> v0ProtonDcaPV006{"v0ProtonDcaPV006", 0.06f, "Min proton DCA to PV: 0.06 cm (Xi* nominal)"};
    Configurable<float> v0ProtonDcaPV007{"v0ProtonDcaPV007", 0.07f, "Min proton DCA to PV: 0.07 cm (dedicated Lambda)"};
    Configurable<float> v0ProtonDcaPV010{"v0ProtonDcaPV010", 0.10f, "Min proton DCA to PV: 0.10 cm"};

    // As requested, keep the Xi(1530)-analysis MINIMUM DCA(V0,PV) convention.
    Configurable<float> v0DcaPV000{"v0DcaPV000", 0.00f, "Min Lambda/V0 DCA to PV: 0.00 cm"};
    Configurable<float> v0DcaPV003{"v0DcaPV003", 0.03f, "Min Lambda/V0 DCA to PV: 0.03 cm"};
    Configurable<float> v0DcaPV010{"v0DcaPV010", 0.10f, "Min Lambda/V0 DCA to PV: 0.10 cm"};

    // Native CascDatas V0-daughter DCA / fitter metric (maximum).
    Configurable<float> v0DcaDaughters1p0{"v0DcaDaughters1p0", 1.0f, "Max V0-daughter DCA/native fitter metric: 1.0"};
    Configurable<float> v0DcaDaughters0p5{"v0DcaDaughters0p5", 0.5f, "Max V0-daughter DCA/native fitter metric: 0.5"};
    Configurable<float> v0DcaDaughters0p1{"v0DcaDaughters0p1", 0.1f, "Max V0-daughter DCA/native fitter metric: 0.1"};

    Configurable<float> v0CosPA097{"v0CosPA097", 0.97f, "Min V0 cosPA: 0.97"};
    Configurable<float> v0CosPA098{"v0CosPA098", 0.98f, "Min V0 cosPA: 0.98 (Xi* nominal)"};
    Configurable<float> v0CosPA09876{"v0CosPA09876", 0.9876f, "Min V0 cosPA: 0.9876 (pp-Xi)"};
    Configurable<float> v0CosPA099{"v0CosPA099", 0.99f, "Min V0 cosPA: 0.99"};
    Configurable<float> v0CosPA0995{"v0CosPA0995", 0.995f, "Min V0 cosPA: 0.995 (dedicated Lambda)"};

    Configurable<float> v0Radius0p9{"v0Radius0p9", 0.9f, "Min V0 radius: 0.9 cm"};
    Configurable<float> v0Radius1p01{"v0Radius1p01", 1.01f, "Min V0 radius: 1.01 cm (Xi* nominal)"};
    Configurable<float> v0Radius1p2{"v0Radius1p2", 1.2f, "Min V0 radius: 1.2 cm (pp-Xi)"};
    Configurable<float> v0Radius2p5{"v0Radius2p5", 2.5f, "Min V0 radius: 2.5 cm"};
    Configurable<float> v0Radius3p0{"v0Radius3p0", 3.0f, "Min V0 radius: 3.0 cm (dedicated Lambda)"};

    Configurable<float> lambdaMass10MeV{"lambdaMass10MeV", 0.010f, "Lambda mass window: 10 MeV/c2"};
    Configurable<float> lambdaMass8MeV{"lambdaMass8MeV", 0.008f, "Lambda mass window: 8 MeV/c2 (Xi* nominal)"};
    Configurable<float> lambdaMass6MeV{"lambdaMass6MeV", 0.006f, "Lambda mass window: 6 MeV/c2"};
    Configurable<float> lambdaMass11p6MeV{"lambdaMass11p6MeV", 0.0116f, "Lambda mass window: 11.6 MeV/c2 (pp-Xi)"};

    Configurable<float> v0LifetimeDefault{"v0LifetimeDefault", 30.f, "Max Lambda proper lifetime: 30 cm"};
    Configurable<float> v0LifetimeVar1{"v0LifetimeVar1", 30.f, "Reserved Lambda lifetime Var1 max; set when defined"};
    Configurable<float> v0LifetimeVar2{"v0LifetimeVar2", 30.f, "Reserved Lambda lifetime Var2 max; set when defined"};

    Configurable<int> v0DaughterRows70{"v0DaughterRows70", 70, "Both Lambda daughters: TPC crossed rows >70"};
    Configurable<int> v0DaughterRowsVar1{"v0DaughterRowsVar1", 70, "Reserved Lambda-daughter crossed-row Var1; set when defined"};
    Configurable<int> v0DaughterRowsVar2{"v0DaughterRowsVar2", 70, "Reserved Lambda-daughter crossed-row Var2; set when defined"};

    // Xi / cascade.
    Configurable<float> cascBachelorDcaPV005{"cascBachelorDcaPV005", 0.05f, "Min bachelor DCA to PV: 0.05 cm (pp-Xi)"};
    Configurable<float> cascBachelorDcaPV006{"cascBachelorDcaPV006", 0.06f, "Min bachelor DCA to PV: 0.06 cm (Xi* nominal)"};
    Configurable<float> cascBachelorDcaPV010{"cascBachelorDcaPV010", 0.10f, "Min bachelor DCA to PV: 0.10 cm"};

    // Native CascDatas cascade-daughter DCA / fitter metric (maximum).
    Configurable<float> cascDcaDaughters1p0{"cascDcaDaughters1p0", 1.0f, "Max cascade-daughter DCA/native fitter metric: 1.0"};
    Configurable<float> cascDcaDaughters0p25{"cascDcaDaughters0p25", 0.25f, "Max cascade-daughter DCA/native fitter metric: 0.25 (Xi* and pp-Xi)"};
    Configurable<float> cascDcaDaughters0p20{"cascDcaDaughters0p20", 0.20f, "Max cascade-daughter DCA/native fitter metric: 0.20"};

    Configurable<float> cascCosPA097{"cascCosPA097", 0.97f, "Min cascade cosPA: 0.97"};
    Configurable<float> cascCosPA098{"cascCosPA098", 0.98f, "Min cascade cosPA: 0.98 (Xi* nominal)"};
    Configurable<float> cascCosPA09947{"cascCosPA09947", 0.9947f, "Min cascade cosPA: 0.9947 (pp-Xi)"};
    Configurable<float> cascCosPA0995{"cascCosPA0995", 0.995f, "Min cascade cosPA: 0.995"};

    Configurable<float> cascRadius0p9{"cascRadius0p9", 0.9f, "Min cascade radius: 0.9 cm"};
    Configurable<float> cascRadius1p0{"cascRadius1p0", 1.0f, "Min cascade radius: 1.0 cm (pp-Xi)"};
    Configurable<float> cascRadius1p01{"cascRadius1p01", 1.01f, "Min cascade radius: 1.01 cm (Xi* nominal)"};
    Configurable<float> cascRadius1p3{"cascRadius1p3", 1.3f, "Min cascade radius: 1.3 cm"};

    Configurable<float> xiMass10MeV{"xiMass10MeV", 0.010f, "Xi mass window: 10 MeV/c2"};
    Configurable<float> xiMass8MeV{"xiMass8MeV", 0.008f, "Xi mass window: 8 MeV/c2 (Xi* nominal)"};
    Configurable<float> xiMass6MeV{"xiMass6MeV", 0.006f, "Xi mass window: 6 MeV/c2"};

    // Xi proper-lifetime working points, mL/p in cm.
    // The default 22.6 cm corresponds approximately to 4.6*c*tau_Xi.
    Configurable<float> xiLifetimeDefault{"xiLifetimeDefault", 22.6f, "Default Xi proper-lifetime maximum mL/p (cm), ~4.6*c*tau_Xi"};
    Configurable<float> xiLifetimeVar1{"xiLifetimeVar1", 15.0f, "Xi proper-lifetime variation 1 maximum mL/p (cm)"};
    Configurable<float> xiLifetimeVar2{"xiLifetimeVar2", 12.0f, "Xi proper-lifetime variation 2 maximum mL/p (cm)"};

    Configurable<float> bachBaryonDCAxyPP{"bachBaryonDCAxyPP", 0.020f, "pp-Xi min |DCAxy| between bachelor and baryon (cm)"};
    Configurable<float> xiRapidityPPMax{"xiRapidityPPMax", 0.5f, "pp-Xi maximum |y_Xi|"};
  } topoCuts;

  // ---------------------------------------------------------------------------
  // Xi(1530) candidate retention.
  // 1.45--1.70 GeV/c2 deliberately contains the existing LS-normalisation
  // sidebands; in particular the right sideband extends to about 1.678 GeV/c2.
  // There is NO Xi(1530) pT cut here.
  // ---------------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "xiStar";
    Configurable<float> massMin{"massMin", 1.45f, "Minimum retained M(Xi pi) (GeV/c2)"};
    Configurable<float> massMax{"massMax", 1.70f, "Maximum retained M(Xi pi) (GeV/c2)"};
    Configurable<float> rapidityMax{"rapidityMax", 0.5f, "Maximum |y_XiStar|"};
  } xiStarCuts;

  // ---------------------------------------------------------------------------
  // External primary K.
  // Raw TPC/TOF PID is retained downstream.  DCA is NOT stored as floats;
  // instead, default/variation decisions are packed into a uint8_t bitmap.
  // ---------------------------------------------------------------------------
  struct : ConfigurableGroup {
    std::string prefix = "kaon";
    Configurable<float> etaMax{"etaMax", 0.8f, "Maximum |eta_K|"};
    Configurable<float> ptMin{"ptMin", 0.15f, "Minimum pT_K (GeV/c)"};
    Configurable<bool> requireGlobalTrack{"requireGlobalTrack", true, "Require isGlobalTrack (same quality baseline as prompt pion)"};
    Configurable<bool> requirePVContributor{"requirePVContributor", true, "Require isPVContributor"};
    Configurable<float> tpcNSigmaEnvelope{"tpcNSigmaEnvelope", 6.f, "Loose |TPC nSigma_K| producer envelope"};

    // Default K DCA working points are pT dependent.  Optional constant Var1/Var2
    // values are disabled by default (<=0).  If enabled, the producer retains the
    // union of the default and enabled alternatives, so looser variations remain
    // available downstream without storing raw DCA values.
    Configurable<bool> dcaXYUsePtDependent{"dcaXYUsePtDependent", true, "Use pT-dependent default DCAxy cut"};
    Configurable<float> dcaXYp0{"dcaXYp0", 0.004f, "K DCAxy(pT) p0 (cm)"};
    Configurable<float> dcaXYp1{"dcaXYp1", 0.013f, "K DCAxy(pT) p1 (cm GeV/c)"};
    Configurable<float> dcaXYMaxConstant{"dcaXYMaxConstant", 0.1f, "Constant default |DCAxy| max if pT-dependent mode is disabled (cm)"};
    Configurable<float> dcaXYVar1Max{"dcaXYVar1Max", -1.f, "Optional K |DCAxy| Var1 max; <=0 disables (cm)"};
    Configurable<float> dcaXYVar2Max{"dcaXYVar2Max", -1.f, "Optional K |DCAxy| Var2 max; <=0 disables (cm)"};

    Configurable<bool> dcaZUsePtDependent{"dcaZUsePtDependent", true, "Use pT-dependent default DCAz cut"};
    Configurable<float> dcaZp0{"dcaZp0", 0.004f, "K DCAz(pT) p0 (cm)"};
    Configurable<float> dcaZp1{"dcaZp1", 0.013f, "K DCAz(pT) p1 (cm GeV/c)"};
    Configurable<float> dcaZMaxConstant{"dcaZMaxConstant", 0.1f, "Constant default |DCAz| max if pT-dependent mode is disabled (cm)"};
    Configurable<float> dcaZVar1Max{"dcaZVar1Max", -1.f, "Optional K |DCAz| Var1 max; <=0 disables (cm)"};
    Configurable<float> dcaZVar2Max{"dcaZVar2Max", -1.f, "Optional K |DCAz| Var2 max; <=0 disables (cm)"};
  } kaonCuts;

  // ---------------------------------------------------------------------------
  // Minimal QA.  Histograms do not contribute to the reduced AOD row size.
  // ---------------------------------------------------------------------------
  HistogramRegistry qa{"qa", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbURL.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    qa.add("hEventCutFlow", "Event cut flow;;events", HistType::kTH1F, {{6, 0.5, 6.5}});
    qa.add("hFT0MPercentile", "FT0M multiplicity percentile of written events;FT0M percentile (%);events", HistType::kTH1F, {{110, 0., 110.}});
    qa.add("hLambdaMass", "Selected loose-envelope Lambda mass;M_{p#pi} (GeV/c^{2});counts", HistType::kTH1F, {{240, 1.09, 1.14}});
    qa.add("hXiMass", "Selected loose-envelope Xi mass;M_{#Lambda#pi} (GeV/c^{2});counts", HistType::kTH1F, {{240, 1.29, 1.35}});
    qa.add("hXiStarMassUS", "US #Xi#pi candidates;M_{#Xi#pi} (GeV/c^{2});counts", HistType::kTH1F, {{500, 1.45, 1.70}});
    qa.add("hXiStarMassLS", "LS #Xi#pi candidates;M_{#Xi#pi} (GeV/c^{2});counts", HistType::kTH1F, {{500, 1.45, 1.70}});
    qa.add("hXiStarMassVsChannel", "#Xi#pi mass vs channel;M_{#Xi#pi} (GeV/c^{2});channel", HistType::kTH2F, {{500, 1.45, 1.70}, {5, -2.5, 2.5}});
    qa.add("hKaonTPCnSigmaVsP", "External K TPC PID;p (GeV/c);n#sigma_{TPC}^{K}", HistType::kTH2F, {{200, 0., 10.}, {160, -8., 8.}});
    qa.add("hKaonTOFnSigmaVsP", "External K TOF PID;p (GeV/c);n#sigma_{TOF}^{K}", HistType::kTH2F, {{200, 0., 10.}, {160, -8., 8.}});
    qa.add("hPiFirstTPCnSigmaVsPt", "Prompt pion TPC PID;p_{T} (GeV/c);n#sigma_{TPC}^{#pi}", HistType::kTH2F, {{200, 0., 10.}, {160, -8., 8.}});
    qa.add("hNObjectsPerWrittenEvent", "Objects in written events;N;type (0=#Xi#pi,1=K)", HistType::kTH2F, {{100, -0.5, 99.5}, {2, -0.5, 1.5}});

    auto h = qa.get<TH1>(HIST("hEventCutFlow"));
    h->GetXaxis()->SetBinLabel(1, "all");
    h->GetXaxis()->SetBinLabel(2, "|z|<10");
    h->GetXaxis()->SetBinLabel(3, "sel8");
    h->GetXaxis()->SetBinLabel(4, "INEL>0");
    h->GetXaxis()->SetBinLabel(5, "XiPi present");
    h->GetXaxis()->SetBinLabel(6, "written");
  }

  // ---------------------------------------------------------------------------
  // Small PODs used only inside one collision.  They keep us from writing any
  // reduced row until we know the event is useful.
  // ---------------------------------------------------------------------------
  struct PromptPionTmp {
    float px{0.f}, py{0.f}, pz{0.f};
    int8_t charge{0};
    int64_t trackIndex{-1};
    uint64_t trackBits{0};
  };

  struct XiTmp {
    float px{0.f}, py{0.f}, pz{0.f};
    float mass{0.f};
    int8_t sign{0};
    int64_t bachelorIndex{-1};
    int64_t posIndex{-1};
    int64_t negIndex{-1};
    uint64_t trackBits{0};
    uint64_t topologyBits{0};
  };

  struct XiStarTmp {
    int8_t channel{0};
    float px{0.f}, py{0.f}, pz{0.f};
    float mass{0.f};
    uint64_t trackBits{0};
    uint64_t topologyBits{0};
    float pionPt{0.f}, pionEta{0.f}, pionPhi{0.f};
    int64_t pionIndex{-1};
    int64_t bachelorIndex{-1};
    int64_t posIndex{-1};
    int64_t negIndex{-1};
  };

  struct KaonTmp {
    int8_t charge{0};
    float px{0.f}, py{0.f}, pz{0.f};
    int64_t trackIndex{-1};
    uint8_t bits{0};
    uint8_t tpcNClsCrossedRows{0};
    float tpcNSigmaKa{999.f};
    float tofNSigmaKa{-999.f};
    bool hasTOF{false};
  };

  // ---------------------------------------------------------------------------
  // Event selection
  // ---------------------------------------------------------------------------
  template <typename TCollision>
  bool passEvent(TCollision const& collision)
  {
    qa.fill(HIST("hEventCutFlow"), 1.0);

    if (std::abs(collision.posZ()) >= eventCuts.zVtxMax) {
      return false;
    }
    qa.fill(HIST("hEventCutFlow"), 2.0);

    if (eventCuts.requireSel8 && !collision.sel8()) {
      return false;
    }
    qa.fill(HIST("hEventCutFlow"), 3.0);

    if (eventCuts.requireINELgt0 && collision.multNTracksPVeta1() < 1) {
      return false;
    }
    qa.fill(HIST("hEventCutFlow"), 4.0);
    return true;
  }

  // ---------------------------------------------------------------------------
  // Prompt-pion selection and bitmap
  // ---------------------------------------------------------------------------
  template <typename TTrack>
  bool passPionHybridPID(TTrack const& track, float nSigmaLimit) const
  {
    const float tpc = track.tpcNSigmaPi();
    if (!track.hasTOF()) {
      return std::abs(tpc) < nSigmaLimit;
    }
    const float tof = track.tofNSigmaPi();
    return std::hypot(tpc, tof) < nSigmaLimit;
  }

  // Threshold-based prompt-pion PID recipe.
  // Below threshold: TPC only, irrespective of whether TOF happens to exist.
  // At/above threshold: circular TPC+TOF.  If TOF is absent, the configurable
  // requireTOFAboveThreshold decides whether the track fails or falls back to TPC.
  template <typename TTrack>
  bool passPromptPionThresholdPID(TTrack const& track,
                                  float tpcNSigmaLimit,
                                  float circularNSigmaLimit) const
  {
    const float tpc = track.tpcNSigmaPi();

    if (track.pt() < piFirstThresholdPidCuts.ptThreshold.value) {
      return std::abs(tpc) < tpcNSigmaLimit;
    }

    if (!track.hasTOF()) {
      if (piFirstThresholdPidCuts.requireTOFAboveThreshold.value) {
        return false;
      }
      return std::abs(tpc) < tpcNSigmaLimit;
    }

    const float tof = track.tofNSigmaPi();
    return std::hypot(tpc, tof) < circularNSigmaLimit;
  }

  static float dcaLimit(float pt, bool usePtDependent, float p0, float p1, float constantMax)
  {
    if (usePtDependent) {
      return p0 + p1 / std::max(pt, kTinyMomentum);
    }
    return constantMax;
  }

  float getMagneticField(uint64_t timestamp)
  {
    auto* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("/GLO/Config/GRPMagField", timestamp);
    if (grpo == nullptr) {
      LOGF(fatal, "GRPMagField object not found in CCDB for timestamp %llu",
           static_cast<unsigned long long>(timestamp));
      return 0.f;
    }
    return 0.1f * grpo->getNominalL3Field(); // kG -> T
  }

  template <typename TTrack>
  bool buildPromptPion(TTrack const& track, PromptPionTmp& out)
  {
    if (std::abs(track.eta()) >= piFirstCuts.etaMax) {
      return false;
    }
    if (track.pt() <= piFirstCuts.ptMin) {
      return false;
    }
    if (piFirstCuts.requireGlobalTrack && !track.isGlobalTrack()) {
      return false;
    }
    if (piFirstCuts.requirePVContributor && !track.isPVContributor()) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < piFirstCuts.minITSInnerBarrelClusters) {
      return false;
    }
    if (track.itsNCls() < piFirstCuts.minITSClusters) {
      return false;
    }

    const float dcaXYDefault = dcaLimit(track.pt(), piFirstCuts.dcaXYUsePtDependent.value,
                                        piFirstCuts.dcaXYp0.value, piFirstCuts.dcaXYp1.value,
                                        piFirstCuts.dcaXYMaxConstant.value);
    if (std::abs(track.dcaXY()) >= dcaXYDefault) {
      return false;
    }

    const float dcaZDefault = dcaLimit(track.pt(), piFirstCuts.dcaZUsePtDependent.value,
                                       piFirstCuts.dcaZp0.value, piFirstCuts.dcaZp1.value,
                                       piFirstCuts.dcaZMaxConstant.value);
    const float dcaZEnvelope = std::max({dcaZDefault, piFirstCuts.dcaZLoose.value,
                                         piFirstCuts.dcaZTightVariation.value});
    if (std::abs(track.dcaZ()) >= dcaZEnvelope) {
      return false;
    }

    if (track.tpcNClsCrossedRows() <= piFirstCuts.tpcRowsLoose) {
      return false;
    }

    // Retain the UNION of all prompt-pion PID choices through their loosest
    // (6-sigma) working points.  The explicit TPC-only envelope is essential:
    // a track with valid TOF may fail the circular hybrid cut while still
    // satisfying the desired TPC-only selection.
    const bool passTPCOnlyEnvelope =
      std::abs(track.tpcNSigmaPi()) < otherTrackPidCuts.nSigma6.value;
    const bool passAvailabilityEnvelope =
      passPionHybridPID(track, otherTrackPidCuts.nSigma6.value);
    const bool passThresholdEnvelope =
      passPromptPionThresholdPID(track,
                                 piFirstThresholdPidCuts.tpcNSigma6.value,
                                 piFirstThresholdPidCuts.circularNSigma6.value);
    if (!(passTPCOnlyEnvelope || passAvailabilityEnvelope || passThresholdEnvelope)) {
      return false;
    }

    uint64_t trackBits = 0;

    setTrackSelBit(trackBits, redxistark::kPiFirstPID3Sigma,
                   passPionHybridPID(track, otherTrackPidCuts.nSigma3.value));
    setTrackSelBit(trackBits, redxistark::kPiFirstPID4Sigma,
                   passPionHybridPID(track, otherTrackPidCuts.nSigma4.value));
    setTrackSelBit(trackBits, redxistark::kPiFirstPID5Sigma,
                   passPionHybridPID(track, otherTrackPidCuts.nSigma5.value));
    setTrackSelBit(trackBits, redxistark::kPiFirstPID6Sigma,
                   passPionHybridPID(track, otherTrackPidCuts.nSigma6.value));

    setTrackSelBit(trackBits, redxistark::kPiFirstPtThresholdPID3Sigma,
                   passPromptPionThresholdPID(track,
                                              piFirstThresholdPidCuts.tpcNSigma3.value,
                                              piFirstThresholdPidCuts.circularNSigma3.value));
    setTrackSelBit(trackBits, redxistark::kPiFirstPtThresholdPID4Sigma,
                   passPromptPionThresholdPID(track,
                                              piFirstThresholdPidCuts.tpcNSigma4.value,
                                              piFirstThresholdPidCuts.circularNSigma4.value));
    setTrackSelBit(trackBits, redxistark::kPiFirstPtThresholdPID5Sigma,
                   passPromptPionThresholdPID(track,
                                              piFirstThresholdPidCuts.tpcNSigma5.value,
                                              piFirstThresholdPidCuts.circularNSigma5.value));
    setTrackSelBit(trackBits, redxistark::kPiFirstPtThresholdPID6Sigma,
                   passPromptPionThresholdPID(track,
                                              piFirstThresholdPidCuts.tpcNSigma6.value,
                                              piFirstThresholdPidCuts.circularNSigma6.value));

    // Pure TPC-only prompt-pion PID bits.  These remain valid irrespective of
    // TOF availability and therefore support an all-pion TPC-only analysis.
    const float absPromptTPCPi = std::abs(track.tpcNSigmaPi());
    setTrackSelBit(trackBits, redxistark::kPiFirstTPC3Sigma,
                   absPromptTPCPi < otherTrackPidCuts.nSigma3.value);
    setTrackSelBit(trackBits, redxistark::kPiFirstTPC4Sigma,
                   absPromptTPCPi < otherTrackPidCuts.nSigma4.value);
    setTrackSelBit(trackBits, redxistark::kPiFirstTPC5Sigma,
                   absPromptTPCPi < otherTrackPidCuts.nSigma5.value);
    setTrackSelBit(trackBits, redxistark::kPiFirstTPC6Sigma,
                   absPromptTPCPi < otherTrackPidCuts.nSigma6.value);

    setTrackSelBit(trackBits, redxistark::kPiFirstTPCRows70, track.tpcNClsCrossedRows() > piFirstCuts.tpcRowsLoose);
    setTrackSelBit(trackBits, redxistark::kPiFirstTPCRows80, track.tpcNClsCrossedRows() > piFirstCuts.tpcRowsDefault);
    setTrackSelBit(trackBits, redxistark::kPiFirstTPCRows90, track.tpcNClsCrossedRows() > piFirstCuts.tpcRowsTight);

    setTrackSelBit(trackBits, redxistark::kPiFirstDCAz1cm, std::abs(track.dcaZ()) < piFirstCuts.dcaZLoose.value);
    setTrackSelBit(trackBits, redxistark::kPiFirstDCAzDefault, std::abs(track.dcaZ()) < dcaZDefault);
    setTrackSelBit(trackBits, redxistark::kPiFirstDCAz0p1cm, std::abs(track.dcaZ()) < piFirstCuts.dcaZTightVariation.value);

    out.px = track.px();
    out.py = track.py();
    out.pz = track.pz();
    out.charge = static_cast<int8_t>(track.sign());
    out.trackIndex = track.globalIndex();
    out.trackBits = trackBits;

    qa.fill(HIST("hPiFirstTPCnSigmaVsPt"), track.pt(), track.tpcNSigmaPi());
    return true;
  }

  // ---------------------------------------------------------------------------
  // External primary-K producer envelope
  // ---------------------------------------------------------------------------
  template <typename TTrack>
  bool buildKaon(TTrack const& track, KaonTmp& out)
  {
    if (std::abs(track.eta()) >= kaonCuts.etaMax) {
      return false;
    }
    if (track.pt() <= kaonCuts.ptMin) {
      return false;
    }
    if (kaonCuts.requireGlobalTrack && !track.isGlobalTrack()) {
      return false;
    }
    if (kaonCuts.requirePVContributor && !track.isPVContributor()) {
      return false;
    }
    if (track.itsNClsInnerBarrel() < piFirstCuts.minITSInnerBarrelClusters) {
      return false;
    }
    if (track.itsNCls() < piFirstCuts.minITSClusters) {
      return false;
    }
    if (track.tpcNClsCrossedRows() <= piFirstCuts.tpcRowsLoose) {
      return false;
    }
    if (std::abs(track.tpcNSigmaKa()) >= kaonCuts.tpcNSigmaEnvelope) {
      return false;
    }

    const float absDcaXY = std::abs(track.dcaXY());
    const float absDcaZ = std::abs(track.dcaZ());
    const float dcaXYDefault = dcaLimit(track.pt(), kaonCuts.dcaXYUsePtDependent.value,
                                        kaonCuts.dcaXYp0.value, kaonCuts.dcaXYp1.value,
                                        kaonCuts.dcaXYMaxConstant.value);
    const float dcaZDefault = dcaLimit(track.pt(), kaonCuts.dcaZUsePtDependent.value,
                                       kaonCuts.dcaZp0.value, kaonCuts.dcaZp1.value,
                                       kaonCuts.dcaZMaxConstant.value);

    const bool passXYDefault = absDcaXY < dcaXYDefault;
    const bool passXYVar1 = kaonCuts.dcaXYVar1Max.value > 0.f && absDcaXY < kaonCuts.dcaXYVar1Max.value;
    const bool passXYVar2 = kaonCuts.dcaXYVar2Max.value > 0.f && absDcaXY < kaonCuts.dcaXYVar2Max.value;
    const bool passZDefault = absDcaZ < dcaZDefault;
    const bool passZVar1 = kaonCuts.dcaZVar1Max.value > 0.f && absDcaZ < kaonCuts.dcaZVar1Max.value;
    const bool passZVar2 = kaonCuts.dcaZVar2Max.value > 0.f && absDcaZ < kaonCuts.dcaZVar2Max.value;

    // Retain the union of enabled DCA working points, but do not store raw DCA.
    if (!(passXYDefault || passXYVar1 || passXYVar2) ||
        !(passZDefault || passZVar1 || passZVar2)) {
      return false;
    }

    uint8_t bits = 0;
    setKaonSelBit(bits, redxistarkaon::kKaonDCAxyDefault, passXYDefault);
    setKaonSelBit(bits, redxistarkaon::kKaonDCAxyVar1, passXYVar1);
    setKaonSelBit(bits, redxistarkaon::kKaonDCAxyVar2, passXYVar2);
    setKaonSelBit(bits, redxistarkaon::kKaonDCAzDefault, passZDefault);
    setKaonSelBit(bits, redxistarkaon::kKaonDCAzVar1, passZVar1);
    setKaonSelBit(bits, redxistarkaon::kKaonDCAzVar2, passZVar2);

    out.charge = static_cast<int8_t>(track.sign());
    out.px = track.px();
    out.py = track.py();
    out.pz = track.pz();
    out.trackIndex = track.globalIndex();
    out.bits = bits;
    out.tpcNClsCrossedRows = toU8(track.tpcNClsCrossedRows());
    out.tpcNSigmaKa = track.tpcNSigmaKa();
    out.hasTOF = track.hasTOF();
    out.tofNSigmaKa = out.hasTOF ? track.tofNSigmaKa() : -999.f;

    const float p = std::sqrt(track.px() * track.px() + track.py() * track.py() + track.pz() * track.pz());
    qa.fill(HIST("hKaonTPCnSigmaVsP"), p, track.tpcNSigmaKa());
    if (track.hasTOF()) {
      qa.fill(HIST("hKaonTOFnSigmaVsP"), p, track.tofNSigmaKa());
    }
    return true;
  }

  // ---------------------------------------------------------------------------
  // Xi candidate from a native CascDatas row.
  // ---------------------------------------------------------------------------
  template <typename TCollision, typename TCascade>
  bool buildXi(TCollision const& collision, TCascade const& casc, XiTmp& out)
  {
    if (casc.sign() == 0) {
      return false;
    }
    if (std::abs(casc.eta()) >= cascadeBaseCuts.cascadeEtaMax) {
      return false;
    }
    if (casc.pt() <= cascadeBaseCuts.cascadePtMin) {
      return false;
    }

    const auto pos = casc.template posTrack_as<TracksFull>();
    const auto neg = casc.template negTrack_as<TracksFull>();
    const auto bach = casc.template bachelor_as<TracksFull>();

    // Fixed broad daughter acceptance.  The >50 crossed-row pp-Xi WP is also
    // explicitly stored as a track-selection bit below.
    if (std::abs(pos.eta()) >= cascadeBaseCuts.daughterEtaMax ||
        std::abs(neg.eta()) >= cascadeBaseCuts.daughterEtaMax ||
        std::abs(bach.eta()) >= cascadeBaseCuts.daughterEtaMax) {
      return false;
    }
    if (pos.tpcNClsCrossedRows() <= cascadeBaseCuts.daughterTPCrowsMin ||
        neg.tpcNClsCrossedRows() <= cascadeBaseCuts.daughterTPCrowsMin ||
        bach.tpcNClsCrossedRows() <= cascadeBaseCuts.daughterTPCrowsMin) {
      return false;
    }

    // Daughter roles:
    //   Xi-      : pos = proton, neg = pion, bachelor = pion
    //   anti-Xi+ : pos = pion,   neg = proton, bachelor = pion
    const bool isXiMinus = casc.sign() < 0;
    const auto& lambdaPion = isXiMinus ? neg : pos;
    const auto& lambdaProton = isXiMinus ? pos : neg;

    // Dedicated Lambda fixed daughter kinematics.
    if (lambdaPion.pt() <= cascadeBaseCuts.lambdaPionPtMin ||
        lambdaProton.pt() <= cascadeBaseCuts.lambdaProtonPtMin) {
      return false;
    }

    const float pLambdaX = casc.pxlambda();
    const float pLambdaY = casc.pylambda();
    const float pLambdaZ = casc.pzlambda();
    ROOT::Math::PxPyPzMVector lambda4(pLambdaX, pLambdaY, pLambdaZ, o2::constants::physics::MassLambda0);
    if (std::abs(lambda4.Rapidity()) >= cascadeBaseCuts.lambdaRapidityMax) {
      return false;
    }

    // PID producer envelope:
    //   * keep TPC-only pion/proton decisions up to 6 sigma;
    //   * for the bachelor, TPC-only 6 sigma is sufficient to retain both the
    //     pp-Xi TPC-only selection and all tighter hybrid selections.
    if (std::abs(bach.tpcNSigmaPi()) >= otherTrackPidCuts.nSigma6.value ||
        std::abs(lambdaPion.tpcNSigmaPi()) >= v0PidCuts.nSigma6.value ||
        std::abs(lambdaProton.tpcNSigmaPr()) >= v0PidCuts.nSigma6.value) {
      return false;
    }

    const float pionDcaPV = std::abs(isXiMinus ? casc.dcanegtopv() : casc.dcapostopv());
    const float protonDcaPV = std::abs(isXiMinus ? casc.dcapostopv() : casc.dcanegtopv());
    const float v0DcaPV = std::abs(casc.dcav0topv(collision.posX(), collision.posY(), collision.posZ()));
    const float v0DcaDaughters = std::abs(casc.dcaV0daughters());
    const float v0CosPA = casc.v0cosPA(collision.posX(), collision.posY(), collision.posZ());
    const float v0Radius = casc.v0radius();
    const float lambdaMassDiff = std::abs(casc.mLambda() - o2::constants::physics::MassLambda0);

    // Lambda proper-lifetime proxy: L(PV->V0)*m_Lambda/p_Lambda.
    const float dxLambda = casc.xlambda() - collision.posX();
    const float dyLambda = casc.ylambda() - collision.posY();
    const float dzLambda = casc.zlambda() - collision.posZ();
    const float lambdaDecayLength = std::sqrt(dxLambda * dxLambda + dyLambda * dyLambda + dzLambda * dzLambda);
    const float pLambda = std::sqrt(pLambdaX * pLambdaX + pLambdaY * pLambdaY + pLambdaZ * pLambdaZ);
    const float v0ProperLifetime = lambdaDecayLength * o2::constants::physics::MassLambda0 / (pLambda + kTinyMomentum);

    const float bachelorDcaPV = std::abs(casc.dcabachtopv());

    // IMPORTANT: these are the native CascDatas stored quantities.
    const float cascDcaDaughters = std::abs(casc.dcacascdaughters());

    const float cascCosPA = casc.casccosPA(collision.posX(), collision.posY(), collision.posZ());
    const float cascRadius = casc.cascradius();
    const float xiMassDiff = std::abs(casc.mXi() - o2::constants::physics::MassXiMinus);
    const float bachBaryonDCAxy = std::abs(casc.bachBaryonDCAxyToPV());

    // Xi proper lifetime and rapidity.
    const float dxXi = casc.x() - collision.posX();
    const float dyXi = casc.y() - collision.posY();
    const float dzXi = casc.z() - collision.posZ();
    const float xiDecayLength = std::sqrt(dxXi * dxXi + dyXi * dyXi + dzXi * dzXi);
    const float pXi = std::sqrt(casc.px() * casc.px() + casc.py() * casc.py() + casc.pz() * casc.pz());
    const float xiProperLifetime = xiDecayLength * casc.mXi() / (pXi + kTinyMomentum);
    ROOT::Math::PxPyPzMVector xi4(casc.px(), casc.py(), casc.pz(), casc.mXi());
    const float xiRapidity = xi4.Rapidity();

    // -----------------------------------------------------------------------
    // Producer envelopes: retain the UNION of all source-defined WPs.
    // -----------------------------------------------------------------------
    const float pionDcaPVEnvelope = std::min({topoCuts.v0PionDcaPV005.value,
                                              topoCuts.v0PionDcaPV006.value,
                                              topoCuts.v0PionDcaPV010.value,
                                              topoCuts.v0PionDcaPV020.value});
    const float protonDcaPVEnvelope = std::min({topoCuts.v0ProtonDcaPV005.value,
                                                topoCuts.v0ProtonDcaPV006.value,
                                                topoCuts.v0ProtonDcaPV007.value,
                                                topoCuts.v0ProtonDcaPV010.value});
    const float v0DcaPVEnvelope = std::min({topoCuts.v0DcaPV000.value,
                                            topoCuts.v0DcaPV003.value,
                                            topoCuts.v0DcaPV010.value});
    const float v0DcaDaughtersEnvelope = std::max({topoCuts.v0DcaDaughters1p0.value,
                                                   topoCuts.v0DcaDaughters0p5.value,
                                                   topoCuts.v0DcaDaughters0p1.value});
    const float v0CosPAEnvelope = std::min({topoCuts.v0CosPA097.value,
                                            topoCuts.v0CosPA098.value,
                                            topoCuts.v0CosPA09876.value,
                                            topoCuts.v0CosPA099.value,
                                            topoCuts.v0CosPA0995.value});
    const float v0RadiusEnvelope = std::min({topoCuts.v0Radius0p9.value,
                                             topoCuts.v0Radius1p01.value,
                                             topoCuts.v0Radius1p2.value,
                                             topoCuts.v0Radius2p5.value,
                                             topoCuts.v0Radius3p0.value});
    const float lambdaMassEnvelope = std::max({topoCuts.lambdaMass10MeV.value,
                                               topoCuts.lambdaMass8MeV.value,
                                               topoCuts.lambdaMass6MeV.value,
                                               topoCuts.lambdaMass11p6MeV.value});
    const float v0LifetimeEnvelope = std::max({topoCuts.v0LifetimeDefault.value,
                                               topoCuts.v0LifetimeVar1.value,
                                               topoCuts.v0LifetimeVar2.value});

    const float bachelorDcaPVEnvelope = std::min({topoCuts.cascBachelorDcaPV005.value,
                                                  topoCuts.cascBachelorDcaPV006.value,
                                                  topoCuts.cascBachelorDcaPV010.value});
    const float cascDcaDaughtersEnvelope = std::max({topoCuts.cascDcaDaughters1p0.value,
                                                     topoCuts.cascDcaDaughters0p25.value,
                                                     topoCuts.cascDcaDaughters0p20.value});
    const float cascCosPAEnvelope = std::min({topoCuts.cascCosPA097.value,
                                              topoCuts.cascCosPA098.value,
                                              topoCuts.cascCosPA09947.value,
                                              topoCuts.cascCosPA0995.value});
    const float xiMassEnvelope = std::max({topoCuts.xiMass10MeV.value,
                                           topoCuts.xiMass8MeV.value,
                                           topoCuts.xiMass6MeV.value});
    const float xiLifetimeEnvelope = std::max({topoCuts.xiLifetimeDefault.value,
                                               topoCuts.xiLifetimeVar1.value,
                                               topoCuts.xiLifetimeVar2.value});

    if (pionDcaPV <= pionDcaPVEnvelope ||
        protonDcaPV <= protonDcaPVEnvelope ||
        v0DcaPV < v0DcaPVEnvelope ||
        v0DcaDaughters >= v0DcaDaughtersEnvelope ||
        v0CosPA <= v0CosPAEnvelope ||
        v0Radius <= v0RadiusEnvelope ||
        v0Radius >= cascadeBaseCuts.maxV0Radius ||
        v0ProperLifetime >= v0LifetimeEnvelope ||
        lambdaMassDiff >= lambdaMassEnvelope ||
        bachelorDcaPV <= bachelorDcaPVEnvelope ||
        cascDcaDaughters >= cascDcaDaughtersEnvelope ||
        cascCosPA <= cascCosPAEnvelope ||
        cascRadius >= cascadeBaseCuts.maxCascRadius ||
        xiMassDiff >= xiMassEnvelope ||
        xiProperLifetime >= xiLifetimeEnvelope) {
      return false;
    }

    uint64_t trackBits = 0;
    uint64_t topologyBits = 0;

    // -----------------------------------------------------------------------
    // Track / PID bits
    // -----------------------------------------------------------------------
    setTrackSelBit(trackBits, redxistark::kXiBachelorPiHybridPID3Sigma,
                   passPionHybridPID(bach, otherTrackPidCuts.nSigma3.value));
    setTrackSelBit(trackBits, redxistark::kXiBachelorPiHybridPID4Sigma,
                   passPionHybridPID(bach, otherTrackPidCuts.nSigma4.value));
    setTrackSelBit(trackBits, redxistark::kXiBachelorPiHybridPID5Sigma,
                   passPionHybridPID(bach, otherTrackPidCuts.nSigma5.value));
    setTrackSelBit(trackBits, redxistark::kXiBachelorPiHybridPID6Sigma,
                   passPionHybridPID(bach, otherTrackPidCuts.nSigma6.value));

    const float absBachelorTPCPi = std::abs(bach.tpcNSigmaPi());
    setTrackSelBit(trackBits, redxistark::kXiBachelorPiTPC3Sigma, absBachelorTPCPi < otherTrackPidCuts.nSigma3.value);
    setTrackSelBit(trackBits, redxistark::kXiBachelorPiTPC4Sigma, absBachelorTPCPi < otherTrackPidCuts.nSigma4.value);
    setTrackSelBit(trackBits, redxistark::kXiBachelorPiTPC4p8Sigma, absBachelorTPCPi < otherTrackPidCuts.bachelorTPC4p8.value);
    setTrackSelBit(trackBits, redxistark::kXiBachelorPiTPC5Sigma, absBachelorTPCPi < otherTrackPidCuts.nSigma5.value);
    setTrackSelBit(trackBits, redxistark::kXiBachelorPiTPC6Sigma, absBachelorTPCPi < otherTrackPidCuts.nSigma6.value);

    const float absLambdaPiTPC = std::abs(lambdaPion.tpcNSigmaPi());
    setTrackSelBit(trackBits, redxistark::kLambdaPiTPC3Sigma, absLambdaPiTPC < v0PidCuts.nSigma3.value);
    setTrackSelBit(trackBits, redxistark::kLambdaPiTPC4Sigma, absLambdaPiTPC < v0PidCuts.nSigma4.value);
    setTrackSelBit(trackBits, redxistark::kLambdaPiTPC4p8Sigma, absLambdaPiTPC < v0PidCuts.nSigma4p8.value);
    setTrackSelBit(trackBits, redxistark::kLambdaPiTPC5Sigma, absLambdaPiTPC < v0PidCuts.nSigma5.value);
    setTrackSelBit(trackBits, redxistark::kLambdaPiTPC6Sigma, absLambdaPiTPC < v0PidCuts.nSigma6.value);

    const float absLambdaPrTPC = std::abs(lambdaProton.tpcNSigmaPr());
    setTrackSelBit(trackBits, redxistark::kLambdaPrTPC3Sigma, absLambdaPrTPC < v0PidCuts.nSigma3.value);
    setTrackSelBit(trackBits, redxistark::kLambdaPrTPC4Sigma, absLambdaPrTPC < v0PidCuts.nSigma4.value);
    setTrackSelBit(trackBits, redxistark::kLambdaPrTPC5Sigma, absLambdaPrTPC < v0PidCuts.nSigma5.value);
    setTrackSelBit(trackBits, redxistark::kLambdaPrTPC6Sigma, absLambdaPrTPC < v0PidCuts.nSigma6.value);

    const bool allCascRows50 =
      pos.tpcNClsCrossedRows() > cascadeBaseCuts.daughterTPCrowsMin &&
      neg.tpcNClsCrossedRows() > cascadeBaseCuts.daughterTPCrowsMin &&
      bach.tpcNClsCrossedRows() > cascadeBaseCuts.daughterTPCrowsMin;
    setTrackSelBit(trackBits, redxistark::kAllCascDaughtersTPCRows50, allCascRows50);

    const bool v0Rows70 =
      pos.tpcNClsCrossedRows() > topoCuts.v0DaughterRows70 &&
      neg.tpcNClsCrossedRows() > topoCuts.v0DaughterRows70;
    const bool v0RowsVar1 =
      pos.tpcNClsCrossedRows() > topoCuts.v0DaughterRowsVar1 &&
      neg.tpcNClsCrossedRows() > topoCuts.v0DaughterRowsVar1;
    const bool v0RowsVar2 =
      pos.tpcNClsCrossedRows() > topoCuts.v0DaughterRowsVar2 &&
      neg.tpcNClsCrossedRows() > topoCuts.v0DaughterRowsVar2;
    setTrackSelBit(trackBits, redxistark::kV0DaughtersTPCRows70, v0Rows70);
    setTrackSelBit(trackBits, redxistark::kV0DaughtersTPCRowsVar1, v0RowsVar1);
    setTrackSelBit(trackBits, redxistark::kV0DaughtersTPCRowsVar2, v0RowsVar2);

    // -----------------------------------------------------------------------
    // V0 / Lambda topology bits
    // -----------------------------------------------------------------------
    setTopoSelBit(topologyBits, redxistark::kV0PionDcaPV005, pionDcaPV > topoCuts.v0PionDcaPV005);
    setTopoSelBit(topologyBits, redxistark::kV0PionDcaPV006, pionDcaPV > topoCuts.v0PionDcaPV006);
    setTopoSelBit(topologyBits, redxistark::kV0PionDcaPV010, pionDcaPV > topoCuts.v0PionDcaPV010);
    setTopoSelBit(topologyBits, redxistark::kV0PionDcaPV020, pionDcaPV > topoCuts.v0PionDcaPV020);

    setTopoSelBit(topologyBits, redxistark::kV0ProtonDcaPV005, protonDcaPV > topoCuts.v0ProtonDcaPV005);
    setTopoSelBit(topologyBits, redxistark::kV0ProtonDcaPV006, protonDcaPV > topoCuts.v0ProtonDcaPV006);
    setTopoSelBit(topologyBits, redxistark::kV0ProtonDcaPV007, protonDcaPV > topoCuts.v0ProtonDcaPV007);
    setTopoSelBit(topologyBits, redxistark::kV0ProtonDcaPV010, protonDcaPV > topoCuts.v0ProtonDcaPV010);

    setTopoSelBit(topologyBits, redxistark::kV0DcaPV000, v0DcaPV > topoCuts.v0DcaPV000);
    setTopoSelBit(topologyBits, redxistark::kV0DcaPV003, v0DcaPV > topoCuts.v0DcaPV003);
    setTopoSelBit(topologyBits, redxistark::kV0DcaPV010, v0DcaPV > topoCuts.v0DcaPV010);

    setTopoSelBit(topologyBits, redxistark::kV0DcaDaughters1p0, v0DcaDaughters < topoCuts.v0DcaDaughters1p0);
    setTopoSelBit(topologyBits, redxistark::kV0DcaDaughters0p5, v0DcaDaughters < topoCuts.v0DcaDaughters0p5);
    setTopoSelBit(topologyBits, redxistark::kV0DcaDaughters0p1, v0DcaDaughters < topoCuts.v0DcaDaughters0p1);

    setTopoSelBit(topologyBits, redxistark::kV0CosPA097, v0CosPA > topoCuts.v0CosPA097);
    setTopoSelBit(topologyBits, redxistark::kV0CosPA098, v0CosPA > topoCuts.v0CosPA098);
    setTopoSelBit(topologyBits, redxistark::kV0CosPA09876, v0CosPA > topoCuts.v0CosPA09876);
    setTopoSelBit(topologyBits, redxistark::kV0CosPA099, v0CosPA > topoCuts.v0CosPA099);
    setTopoSelBit(topologyBits, redxistark::kV0CosPA0995, v0CosPA > topoCuts.v0CosPA0995);

    setTopoSelBit(topologyBits, redxistark::kV0Radius0p9, v0Radius > topoCuts.v0Radius0p9);
    setTopoSelBit(topologyBits, redxistark::kV0Radius1p01, v0Radius > topoCuts.v0Radius1p01);
    setTopoSelBit(topologyBits, redxistark::kV0Radius1p2, v0Radius > topoCuts.v0Radius1p2);
    setTopoSelBit(topologyBits, redxistark::kV0Radius2p5, v0Radius > topoCuts.v0Radius2p5);
    setTopoSelBit(topologyBits, redxistark::kV0Radius3p0, v0Radius > topoCuts.v0Radius3p0);

    setTopoSelBit(topologyBits, redxistark::kLambdaMass10MeV, lambdaMassDiff < topoCuts.lambdaMass10MeV);
    setTopoSelBit(topologyBits, redxistark::kLambdaMass8MeV, lambdaMassDiff < topoCuts.lambdaMass8MeV);
    setTopoSelBit(topologyBits, redxistark::kLambdaMass6MeV, lambdaMassDiff < topoCuts.lambdaMass6MeV);
    setTopoSelBit(topologyBits, redxistark::kLambdaMass11p6MeV, lambdaMassDiff < topoCuts.lambdaMass11p6MeV);

    setTopoSelBit(topologyBits, redxistark::kV0LifetimeDefault, v0ProperLifetime < topoCuts.v0LifetimeDefault);
    setTopoSelBit(topologyBits, redxistark::kV0LifetimeVar1, v0ProperLifetime < topoCuts.v0LifetimeVar1);
    setTopoSelBit(topologyBits, redxistark::kV0LifetimeVar2, v0ProperLifetime < topoCuts.v0LifetimeVar2);

    // -----------------------------------------------------------------------
    // Xi / cascade topology bits
    // -----------------------------------------------------------------------
    setTopoSelBit(topologyBits, redxistark::kCascBachelorDcaPV005, bachelorDcaPV > topoCuts.cascBachelorDcaPV005);
    setTopoSelBit(topologyBits, redxistark::kCascBachelorDcaPV006, bachelorDcaPV > topoCuts.cascBachelorDcaPV006);
    setTopoSelBit(topologyBits, redxistark::kCascBachelorDcaPV010, bachelorDcaPV > topoCuts.cascBachelorDcaPV010);

    setTopoSelBit(topologyBits, redxistark::kCascDcaDaughters1p0, cascDcaDaughters < topoCuts.cascDcaDaughters1p0);
    setTopoSelBit(topologyBits, redxistark::kCascDcaDaughters0p25, cascDcaDaughters < topoCuts.cascDcaDaughters0p25);
    setTopoSelBit(topologyBits, redxistark::kCascDcaDaughters0p20, cascDcaDaughters < topoCuts.cascDcaDaughters0p20);

    setTopoSelBit(topologyBits, redxistark::kCascCosPA097, cascCosPA > topoCuts.cascCosPA097);
    setTopoSelBit(topologyBits, redxistark::kCascCosPA098, cascCosPA > topoCuts.cascCosPA098);
    setTopoSelBit(topologyBits, redxistark::kCascCosPA09947, cascCosPA > topoCuts.cascCosPA09947);
    setTopoSelBit(topologyBits, redxistark::kCascCosPA0995, cascCosPA > topoCuts.cascCosPA0995);

    setTopoSelBit(topologyBits, redxistark::kCascRadius0p9, cascRadius > topoCuts.cascRadius0p9);
    setTopoSelBit(topologyBits, redxistark::kCascRadius1p0, cascRadius > topoCuts.cascRadius1p0);
    setTopoSelBit(topologyBits, redxistark::kCascRadius1p01, cascRadius > topoCuts.cascRadius1p01);
    setTopoSelBit(topologyBits, redxistark::kCascRadius1p3, cascRadius > topoCuts.cascRadius1p3);

    setTopoSelBit(topologyBits, redxistark::kXiMass10MeV, xiMassDiff < topoCuts.xiMass10MeV);
    setTopoSelBit(topologyBits, redxistark::kXiMass8MeV, xiMassDiff < topoCuts.xiMass8MeV);
    setTopoSelBit(topologyBits, redxistark::kXiMass6MeV, xiMassDiff < topoCuts.xiMass6MeV);

    setTopoSelBit(topologyBits, redxistark::kXiLifetimeDefault, xiProperLifetime < topoCuts.xiLifetimeDefault.value);
    setTopoSelBit(topologyBits, redxistark::kXiLifetimeVar1, xiProperLifetime < topoCuts.xiLifetimeVar1.value);
    setTopoSelBit(topologyBits, redxistark::kXiLifetimeVar2, xiProperLifetime < topoCuts.xiLifetimeVar2.value);

    setTopoSelBit(topologyBits, redxistark::kBachBaryonDCAxy0020, bachBaryonDCAxy > topoCuts.bachBaryonDCAxyPP);
    setTopoSelBit(topologyBits, redxistark::kXiRapidity05, std::abs(xiRapidity) < topoCuts.xiRapidityPPMax);

    out.px = casc.px();
    out.py = casc.py();
    out.pz = casc.pz();
    out.mass = casc.mXi();
    out.sign = static_cast<int8_t>(casc.sign());
    out.bachelorIndex = bach.globalIndex();
    out.posIndex = pos.globalIndex();
    out.negIndex = neg.globalIndex();
    out.trackBits = trackBits;
    out.topologyBits = topologyBits;

    qa.fill(HIST("hLambdaMass"), casc.mLambda());
    qa.fill(HIST("hXiMass"), casc.mXi());
    return true;
  }

  // ---------------------------------------------------------------------------
  // Xi + prompt-pion -> Xi(1530) candidate.
  // ---------------------------------------------------------------------------
  bool buildXiStar(XiTmp const& xi, PromptPionTmp const& pion, XiStarTmp& out)
  {
    // Never allow the prompt pion to reuse a cascade daughter track.
    if (pion.trackIndex == xi.bachelorIndex ||
        pion.trackIndex == xi.posIndex ||
        pion.trackIndex == xi.negIndex) {
      return false;
    }

    ROOT::Math::PxPyPzMVector xi4(xi.px, xi.py, xi.pz, xi.mass);
    ROOT::Math::PxPyPzMVector pion4(pion.px, pion.py, pion.pz, o2::constants::physics::MassPionCharged);
    const auto mother = xi4 + pion4;

    if (std::abs(mother.Rapidity()) >= xiStarCuts.rapidityMax) {
      return false;
    }
    if (mother.M() <= xiStarCuts.massMin || mother.M() >= xiStarCuts.massMax) {
      return false;
    }

    int8_t channel = 0;
    if (xi.sign < 0) { // Xi-
      channel = pion.charge > 0 ? redxistark::kXiStarUS : redxistark::kXiStarLS;
    } else { // anti-Xi+
      channel = pion.charge < 0 ? redxistark::kAntiXiStarUS : redxistark::kAntiXiStarLS;
    }

    out.channel = channel;
    out.px = mother.Px();
    out.py = mother.Py();
    out.pz = mother.Pz();
    out.mass = mother.M();
    out.trackBits = xi.trackBits | pion.trackBits;
    out.topologyBits = xi.topologyBits;
    out.pionPt = std::hypot(pion.px, pion.py);
    out.pionEta = out.pionPt > 0.f ? std::asinh(pion.pz / out.pionPt) : 0.f;
    out.pionPhi = std::atan2(pion.py, pion.px);
    out.pionIndex = pion.trackIndex;
    out.bachelorIndex = xi.bachelorIndex;
    out.posIndex = xi.posIndex;
    out.negIndex = xi.negIndex;

    if (std::abs(channel) == 1) {
      qa.fill(HIST("hXiStarMassUS"), out.mass);
    } else {
      qa.fill(HIST("hXiStarMassLS"), out.mass);
    }
    qa.fill(HIST("hXiStarMassVsChannel"), out.mass, static_cast<float>(channel));
    return true;
  }

  // ---------------------------------------------------------------------------
  // Main producer
  // ---------------------------------------------------------------------------
  void processData(EventCandidates::iterator const& collision,
                   aod::BCsWithTimestamps const&,
                   TracksFull const& tracks,
                   aod::CascDatas const& cascades)
  {
    if (!passEvent(collision)) {
      return;
    }

    const float ft0MPercentile = collision.centFT0M();
    qa.fill(HIST("hFT0MPercentile"), ft0MPercentile);

    std::vector<PromptPionTmp> promptPions;
    std::vector<KaonTmp> kaons;
    std::vector<XiTmp> xis;
    std::vector<XiStarTmp> xiStars;

    // CPU protection for Hyperloop: Xi candidates are much rarer than ordinary
    // primary tracks.  Select cascades first and do not scan all primary tracks
    // in events that contain no Xi candidate inside the loose envelope.
    if (cascades.size() == 0) {
      return;
    }
    xis.reserve(cascades.size());
    for (auto const& casc : cascades) {
      XiTmp xi;
      if (buildXi(collision, casc, xi)) {
        xis.push_back(xi);
      }
    }
    if (xis.empty()) {
      return;
    }

    // One pass over primary tracks: build prompt-pion and external-K lists.
    for (auto const& track : tracks) {
      PromptPionTmp pion;
      if (buildPromptPion(track, pion)) {
        promptPions.push_back(pion);
      }

      KaonTmp kaon;
      if (buildKaon(track, kaon)) {
        kaons.push_back(kaon);
      }
    }

    if (promptPions.empty()) {
      return;
    }

    // Build both US and LS Xi-pi candidates.  No k* selection is made here.
    xiStars.reserve(xis.size() * std::min<std::size_t>(promptPions.size(), 8));
    for (auto const& xi : xis) {
      for (auto const& pion : promptPions) {
        XiStarTmp candidate;
        if (buildXiStar(xi, pion, candidate)) {
          xiStars.push_back(candidate);
        }
      }
    }

    if (xiStars.empty()) {
      return;
    }

    qa.fill(HIST("hEventCutFlow"), 5.0);

    const auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    const float bz = getMagneticField(bc.timestamp());

    redEvents(collision.posZ(), collision.numContrib(), ft0MPercentile, bz);
    const auto redEventIndex = redEvents.lastIndex();

    for (auto const& candidate : xiStars) {
      xiStarCandidates(redEventIndex,
                       candidate.channel,
                       candidate.px,
                       candidate.py,
                       candidate.pz,
                       candidate.mass,
                       candidate.trackBits,
                       candidate.topologyBits,
                       candidate.pionPt,
                       candidate.pionEta,
                       candidate.pionPhi,
                       candidate.pionIndex,
                       candidate.bachelorIndex,
                       candidate.posIndex,
                       candidate.negIndex);
    }

    for (auto const& kaon : kaons) {
      kaonCandidates(redEventIndex,
                     kaon.charge,
                     kaon.px,
                     kaon.py,
                     kaon.pz,
                     kaon.trackIndex,
                     kaon.bits,
                     kaon.tpcNClsCrossedRows,
                     kaon.tpcNSigmaKa,
                     kaon.tofNSigmaKa,
                     kaon.hasTOF);
    }

    qa.fill(HIST("hNObjectsPerWrittenEvent"), static_cast<float>(xiStars.size()), 0.f);
    qa.fill(HIST("hNObjectsPerWrittenEvent"), static_cast<float>(kaons.size()), 1.f);
    qa.fill(HIST("hEventCutFlow"), 6.0);
  }

  PROCESS_SWITCH(xi1530kaonreducedtable, processData, "Produce reduced Xi(1530)-K femtoscopy tables", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<xi1530kaonreducedtable>(cfgc)};
}
