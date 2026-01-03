// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file trackHistManager.h
/// \brief histogram manager for track histograms
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTO_CORE_TRACKHISTMANAGER_H_
#define PWGCF_FEMTO_CORE_TRACKHISTMANAGER_H_

#include "PWGCF/Femto/Core/femtoUtils.h"
#include "PWGCF/Femto/Core/histManager.h"
#include "PWGCF/Femto/Core/modes.h"

#include "CommonConstants/MathConstants.h"
#include "Framework/Configurable.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/HistogramSpec.h"

#include "TH1.h"

#include <array>
#include <map>
#include <string>
#include <string_view>
#include <vector>

namespace o2::analysis::femto
{
namespace trackhistmanager
{

// enum for track histograms
enum TrackHist {
  // kinemtics
  kPt,
  kEta,
  kPhi,
  kSign,
  // qa variables
  kPAtPv,
  kPTpc,
  kItsCluster,
  kItsClusterIb,
  kTpcCrossedRows,
  kTpcCluster,
  kTpcClusterOverCrossedRows,
  kTpcClusterShared,
  kTpcClusterFractionShared,
  // 2d qa
  kPtVsEta,
  kPtVsPhi,
  kPhiVsEta,
  kPtVsItsCluster,
  kPtVsTpcCluster,
  kPtVsTpcCrossedRows,
  kPtVsTpcClusterOverCrossedRows,
  kPtVsTpcClusterShared,
  kPtVsTpcClusterFractionShared,
  kTpcClusterVsTpcCrossedRows,
  kTpcClusterVsTpcClusterShared,
  kPtVsDcaxy,
  kPtVsDcaz,
  kPtVsDca,
  // its pid
  kItsSignal,
  kItsElectron,
  kItsPion,
  kItsKaon,
  kItsProton,
  kItsDeuteron,
  kItsTriton,
  kItsHelium,
  // tpc pid
  kTpcSignal,
  kTpcElectron,
  kTpcPion,
  kTpcKaon,
  kTpcProton,
  kTpcDeuteron,
  kTpcTriton,
  kTpcHelium,
  // tof pid
  kTofBeta,
  kTofMass,
  kTofElectron,
  kTofPion,
  kTofKaon,
  kTofProton,
  kTofDeuteron,
  kTofTriton,
  kTofHelium,
  // tpc+its pid
  kTpcitsElectron,
  kTpcitsPion,
  kTpcitsKaon,
  kTpcitsProton,
  kTpcitsDeuteron,
  kTpcitsTriton,
  kTpcitsHelium,
  // tpc+tof pid
  kTpctofElectron,
  kTpctofPion,
  kTpctofKaon,
  kTpctofProton,
  kTpctofDeuteron,
  kTpctofTriton,
  kTpctofHelium,
  // mc
  kOrigin,
  kPdg,
  kPdgMother,
  kPdgPartonicMother,
  kTruePtVsPt,
  kTrueEtaVsEta,
  kTruePhiVsPhi,
  // histograms for fraction estimation of tracks
  kNoMcParticle,
  kPrimary,
  kFromWrongCollision,
  kFromMaterial,
  kMissidentified,
  kSecondary1,
  kSecondary2,
  kSecondary3,
  kSecondaryOther,

  kTrackHistLast
};

constexpr std::size_t MaxSecondary = 3;

template <const char* Prefix>
struct ConfTrackBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};
  o2::framework::ConfigurableAxis pdgCodes{"pdgCodes", {{8001, -4000.5, 4000.5}}, "MC ONLY: PDG codes of selected tracks"};
};

constexpr const char PrefixTrackBinning1[] = "TrackBinning1";
constexpr const char PrefixTrackBinning2[] = "TrackBinning2";
constexpr const char PrefixResonancePosDauBinning[] = "ResonancePosDauBinning";
constexpr const char PrefixResonanceNegDauBinning[] = "ResonanceNegDauBinning";
constexpr const char PrefixV0PosDauBinning[] = "V0PosDauBinning";
constexpr const char PrefixV0NegDauBinning[] = "V0NegDauBinning";
constexpr const char PrefixCascadePosDauBinning[] = "CascadePosDauBinning";
constexpr const char PrefixCascadeNegDauBinning[] = "CascadeNegDauBinning";
constexpr const char PrefixCascadeBachelorBinning[] = "CascadeBachelorBinning";
constexpr const char PrefixKinkChaDauBinning[] = "KinkChaDauBinning";

using ConfTrackBinning1 = ConfTrackBinning<PrefixTrackBinning1>;
using ConfTrackBinning2 = ConfTrackBinning<PrefixTrackBinning2>;
using ConfResonancePosDauBinning = ConfTrackBinning<PrefixResonancePosDauBinning>;
using ConfResonanceNegDauBinning = ConfTrackBinning<PrefixResonanceNegDauBinning>;
using ConfV0PosDauBinning = ConfTrackBinning<PrefixV0PosDauBinning>;
using ConfV0NegDauBinning = ConfTrackBinning<PrefixV0NegDauBinning>;
using ConfCascadePosDauBinning = ConfTrackBinning<PrefixCascadePosDauBinning>;
using ConfCascadeNegDauBinning = ConfTrackBinning<PrefixCascadeNegDauBinning>;
using ConfCascadeBachelorBinning = ConfTrackBinning<PrefixCascadeBachelorBinning>;
using ConfKinkChaDauBinning = ConfTrackBinning<PrefixKinkChaDauBinning>;

template <const char* Prefix>
struct ConfTrackQaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<int> momentumType{"momentumType", static_cast<int>(modes::MomentumType::kPAtPv), "Momentum on x-axis (0->Pt, 1->P at PV, 2->P at TPC inner wall)"};
  o2::framework::Configurable<bool> plot2d{"plot2d", true, "Generate various 2D QA plots"};
  o2::framework::Configurable<bool> plotElectronPid{"plotElectronPid", true, "Generate plots for Electron PID"};
  o2::framework::Configurable<bool> plotPionPid{"plotPionPid", true, "Generate plots for Pion PID"};
  o2::framework::Configurable<bool> plotKaonPid{"plotKaonPid", true, "Generate plots for Kaon PID"};
  o2::framework::Configurable<bool> plotProtonPid{"plotProtonPid", true, "Generate plots for Proton PID"};
  o2::framework::Configurable<bool> plotDeuteronPid{"plotDeuteronPid", true, "Generate plots for Deuteron PID"};
  o2::framework::Configurable<bool> plotTritonPid{"plotTritonPid", true, "Generate plots for Triton PID"};
  o2::framework::Configurable<bool> plotHeliumPid{"plotHeliumPid", true, "Generate plots for Helium PID"};
  o2::framework::Configurable<bool> plotOrigins{"plotOrigins", true, "MC ONLY: Plot pt vs DCAxy vs DCAz for different particle origins"};
  o2::framework::Configurable<std::vector<int>> pdgCodesForMothersOfSecondary{"pdgCodesForMothersOfSecondary", {3122}, "MC ONLY: PDG codes of mothers of secondaries (Max 3 will be considered)"};
  o2::framework::ConfigurableAxis itsCluster{"itsCluster", {{8, -0.5, 7.5}}, "ITS cluster"};
  o2::framework::ConfigurableAxis itsClusterIb{"itsClusterIb", {{4, -0.5, 3.5}}, "ITS cluster in inner barrel"};
  o2::framework::ConfigurableAxis tpcCrossedRows{"tpcCrossedRows", {{161, -0.5, 160.5}}, "TPC cluster"};
  o2::framework::ConfigurableAxis tpcCluster{"tpcCluster", {{161, -0.5, 160.5}}, "TPC cluster"};
  o2::framework::ConfigurableAxis tpcClusterOverCrossedRows{"tpcClusterOverCrossedRows", {{75, 0, 1.5}}, "TPC cluster over TPC crossed rows"};
  o2::framework::ConfigurableAxis tpcClusterShared{"tpcClusterShared", {{161, -0.5, 160.5}}, "TPC cluster shared"};
  o2::framework::ConfigurableAxis tpcClusterFractionShared{"tpcClusterFractionShared", {{60, 0, 1.2}}, "TPC cluster fraction shared"};
  o2::framework::ConfigurableAxis dcaXy{"dcaXy", {{300, -0.3, 0.3}}, "DCA_xy"};
  o2::framework::ConfigurableAxis dcaZ{"dcaZ", {{300, -0.3, 0.3}}, "DCA_Z"};
  o2::framework::ConfigurableAxis dca{"dca", {{300, 0, 0.3}}, "DCA"};
  o2::framework::ConfigurableAxis p{"p", {{300, 0, 6}}, "Momentum axis"};
  o2::framework::ConfigurableAxis itsSignal{"itsSignal", {{150, 0, 15}}, "ITS Signal"};
  o2::framework::ConfigurableAxis itsElectron{"itsElectron", {{300, -3, 3}}, "ITS PID for electron"};
  o2::framework::ConfigurableAxis itsPion{"itsPion", {{300, -3, 3}}, "ITS PID for pion"};
  o2::framework::ConfigurableAxis itsKaon{"itsKaon", {{300, -3, 3}}, "ITS PID for kaon"};
  o2::framework::ConfigurableAxis itsProton{"itsProton", {{300, -3, 3}}, "ITS PID for proton"};
  o2::framework::ConfigurableAxis itsDeuteron{"itsDeuteron", {{300, -3, 3}}, "ITS PID for deuteron"};
  o2::framework::ConfigurableAxis itsTriton{"itsTriton", {{300, -3, 3}}, "ITS PID for triton"};
  o2::framework::ConfigurableAxis itsHelium{"itsHelium", {{300, -3, 3}}, "ITS PID for helium"};
  o2::framework::ConfigurableAxis tpcSignal{"tpcSignal", {{150, 0, 150}}, "TPC Signal"};
  o2::framework::ConfigurableAxis tpcElectron{"tpcElectron", {{300, -3, 3}}, "TPC PID for electron"};
  o2::framework::ConfigurableAxis tpcPion{"tpcPion", {{300, -3, 3}}, "TPC PID for pion"};
  o2::framework::ConfigurableAxis tpcKaon{"tpcKaon", {{300, -3, 3}}, "TPC PID for kaon"};
  o2::framework::ConfigurableAxis tpcProton{"tpcProton", {{300, -3, 3}}, "TPC PID for proton"};
  o2::framework::ConfigurableAxis tpcDeuteron{"tpcDeuteron", {{300, -3, 3}}, "TPC PID for deuteron"};
  o2::framework::ConfigurableAxis tpcTriton{"tpcTriton", {{300, -3, 3}}, "TPC PID for triton"};
  o2::framework::ConfigurableAxis tpcHelium{"tpcHelium", {{300, -3, 3}}, "TPC PID for helium"};
  o2::framework::ConfigurableAxis tofBeta{"tofBeta", {{150, 0, 1.5}}, "TOF Signal"};
  o2::framework::ConfigurableAxis tofMass{"tofMass", {{150, 0, 1.5}}, "TOF Mass"};
  o2::framework::ConfigurableAxis tofElectron{"tofElectron", {{300, -3, 3}}, "TOF PID for electron"};
  o2::framework::ConfigurableAxis tofPion{"tofPion", {{300, -3, 3}}, "TOF PID for pion"};
  o2::framework::ConfigurableAxis tofKaon{"tofKaon", {{300, -3, 3}}, "TOF PID for kaon"};
  o2::framework::ConfigurableAxis tofProton{"tofProton", {{300, -3, 3}}, "TOF PID for proton"};
  o2::framework::ConfigurableAxis tofDeuteron{"tofDeuteron", {{300, -3, 3}}, "TOF PID for deuteron"};
  o2::framework::ConfigurableAxis tofTriton{"tofTriton", {{300, -3, 3}}, "TOF PID for triton"};
  o2::framework::ConfigurableAxis tofHelium{"tofHelium", {{300, -3, 3}}, "TOF PID for helium"};
  o2::framework::ConfigurableAxis tpcitsElectron{"tpcitsElectron", {{300, 0, 3}}, "tpcits PID for electron"};
  o2::framework::ConfigurableAxis tpcitsPion{"tpcitsPion", {{300, 0, 3}}, "TPCITS PID for pion"};
  o2::framework::ConfigurableAxis tpcitsKaon{"tpcitsKaon", {{300, 0, 3}}, "TPCITS PID for kaon"};
  o2::framework::ConfigurableAxis tpcitsProton{"tpcitsProton", {{300, 0, 3}}, "TPCITS PID for proton"};
  o2::framework::ConfigurableAxis tpcitsDeuteron{"tpcitsDeuteron", {{300, 0, 3}}, "TPCITS PID for deuteron"};
  o2::framework::ConfigurableAxis tpcitsTriton{"tpcitsTriton", {{300, 0, 3}}, "TPCITS PID for triton"};
  o2::framework::ConfigurableAxis tpcitsHelium{"tpcitsHelium", {{300, 0, 3}}, "TPCITS PID for helium"};
  o2::framework::ConfigurableAxis tpctofElectron{"tpctofElectron", {{300, 0, 3}}, "TPCTOF PID for electron"};
  o2::framework::ConfigurableAxis tpctofPion{"tpctofPion", {{300, 0, 3}}, "TPCTOF PID for pion"};
  o2::framework::ConfigurableAxis tpctofKaon{"tpctofKaon", {{300, 0, 3}}, "TPCTOF PID for kaon"};
  o2::framework::ConfigurableAxis tpctofProton{"tpctofProton", {{300, 0, 3}}, "TPCTOF PID for proton"};
  o2::framework::ConfigurableAxis tpctofDeuteron{"tpctofDeuteron", {{300, 0, 3}}, "TPCTOF PID for deuteron"};
  o2::framework::ConfigurableAxis tpctofTriton{"tpctofTriton", {{300, 0, 3}}, "TPCTOF PID for triton"};
  o2::framework::ConfigurableAxis tpctofHelium{"tpctofHelium", {{300, 0, 3}}, "TPCTOF PID for helium"};
};

constexpr const char PrefixTrackQaBinning1[] = "TrackQaBinning1";
constexpr const char PrefixTrackQaBinning2[] = "TrackQaBinning2";
constexpr const char PrefixResonancePosDauQaBinning[] = "ResonancePosDauQaBinning";
constexpr const char PrefixResonanceNegDauQaBinning[] = "ResonanceNegDauQaBinning";
constexpr const char PrefixV0PosDauQaBinning[] = "V0PosDauQaBinning";
constexpr const char PrefixV0NegDauQaBinning[] = "V0NegDauQaBinning";
constexpr const char PrefixCascadePosDauQaBinning[] = "CascadePosDauQaBinning";
constexpr const char PrefixCascadeNegDauQaBinning[] = "CascadeNegDauQaBinning";
constexpr const char PrefixCascadeBachelorQaBinning[] = "CascadeBachelorQaBinning";
constexpr const char PrefixKinkChaDauQaBinning[] = "KinkChaDauQaBinning";

using ConfTrackQaBinning1 = ConfTrackQaBinning<PrefixTrackQaBinning1>;
using ConfTrackQaBinning2 = ConfTrackQaBinning<PrefixTrackQaBinning2>;
using ConfResonancePosDauQaBinning = ConfTrackQaBinning<PrefixResonancePosDauQaBinning>;
using ConfResonanceNegDauQaBinning = ConfTrackQaBinning<PrefixResonanceNegDauQaBinning>;
using ConfV0PosDauQaBinning = ConfTrackQaBinning<PrefixV0PosDauQaBinning>;
using ConfV0NegDauQaBinning = ConfTrackQaBinning<PrefixV0NegDauQaBinning>;
using ConfCascadePosDauQaBinning = ConfTrackQaBinning<PrefixCascadePosDauQaBinning>;
using ConfCascadeNegDauQaBinning = ConfTrackQaBinning<PrefixCascadeNegDauQaBinning>;
using ConfCascadeBachelorQaBinning = ConfTrackQaBinning<PrefixCascadeBachelorQaBinning>;
using ConfKinkChaDauQaBinning = ConfTrackQaBinning<PrefixKinkChaDauQaBinning>;

// must be in sync with enum TrackVariables
// the enum gives the correct index in the array
constexpr std::array<histmanager::HistInfo<TrackHist>, kTrackHistLast>
  HistTable = {
    {
      {kPt, o2::framework::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
      {kEta, o2::framework::kTH1F, "hEta", "Pseudorapidity; #eta; Entries"},
      {kPhi, o2::framework::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
      {kSign, o2::framework::kTH1F, "hSign", "Sign of charge ; Sign; Entries"},
      {kPAtPv, o2::framework::kTH1F, "hPAtPv", "Momentum at Primary vertex; p_{vertex}; Entries"},
      {kPTpc, o2::framework::kTH1F, "hPTpc", "Momentum at inner wall of TPC; p_{TPC}; Entries"},
      {kItsCluster, o2::framework::kTH1F, "hItsCluster", "ITS cluster; ITS cluster; Entries"},
      {kItsClusterIb, o2::framework::kTH1F, "hItsClusterIb", "ITS cluster in inner barrel; ITS IB cluster; Entries"},
      {kTpcCrossedRows, o2::framework::kTH1F, "hTpcCrossedRows", "TPC crossed rows; TPC crossed rows; Entries"},
      {kTpcCluster, o2::framework::kTH1F, "hTpcCluster", "TPC cluster found; TPC cluster found; Entries"},
      {kTpcClusterOverCrossedRows, o2::framework::kTH1F, "hTpcClusterOverCrossedRows", "TPC cluster found  over TPC crossed rows; TPC cluster found / Tpc crossed rows; Entries"},
      {kTpcClusterShared, o2::framework::kTH1F, "hTpcClusterShared", "TPC cluster shared; TPC cluster shared ; Entries"},
      {kTpcClusterFractionShared, o2::framework::kTH1F, "hTpcClusterFractionShared", "TPC cluster fraction shared; TPC cluster found / TPC cluster shared ; Entries"},
      {kPtVsEta, o2::framework::kTH2F, "hPtVsEta", "p_{T} vs #eta; p_{T} (GeV/#it{c}) ; #eta"},
      {kPtVsPhi, o2::framework::kTH2F, "hPtVsPhi", "p_{T} vs #varphi; p_{T} (GeV/#it{c}) ; #varphi"},
      {kPhiVsEta, o2::framework::kTH2F, "hPhiVsEta", "#varphi vs #eta; #varphi ; #eta"},
      {kPtVsItsCluster, o2::framework::kTH2F, "hPtVsItsCluster", "p_{T} vs ITS cluster; p_{T} (GeV/#it{c}) ; ITS cluster"},
      {kPtVsTpcCluster, o2::framework::kTH2F, "hPtVsTpcCluster", "p_{T} vs TPC cluster found; p_{T} (GeV/#it{c}) ; TPC cluster found"},
      {kPtVsTpcCrossedRows, o2::framework::kTH2F, "hPtVsTpcCrossedRows", "p_{T} vs TPC crossed rows; p_{T} (GeV/#it{c}) ; TPC crossed rows"},
      {kPtVsTpcClusterOverCrossedRows, o2::framework::kTH2F, "hPtVsTpcClusterOverCrossedRows", "p_{T} vs TPC cluster found over crossed rows; p_{T} (GeV/#it{c}) ; TPC cluster found / TPC crossed rows"},
      {kPtVsTpcClusterShared, o2::framework::kTH2F, "hPtVsTpcClusterShared", "p_{T} vs TPC cluster shared; p_{T} (GeV/#it{c}) ; TPC cluster shared"},
      {kPtVsTpcClusterFractionShared, o2::framework::kTH2F, "hPtVsTpcClusterSharedFraction", "p_{T} vs TPC cluster shared over TPC cluster found; p_{T} (GeV/#it{c}) ; TPC cluster shared / TPC cluster found"},
      {kTpcClusterVsTpcCrossedRows, o2::framework::kTH2F, "hTpcClusterVsTpcCrossedRows", "TPC cluster found vs TPC crossed rows; TPC cluster found; TPC crossed rows"},
      {kTpcClusterVsTpcClusterShared, o2::framework::kTH2F, "hTpcClusterVsTpcClusterShared", "TPC cluster found vs TPC cluster shared; TPC cluster found; TPC cluster shared"},
      {kPtVsDcaxy, o2::framework::kTH2F, "hPtVsDcaxy", "p_{T} vs DCA_{XY}; p_{T} (GeV/#it{c}); DCA_{XY} (cm)"},
      {kPtVsDcaz, o2::framework::kTH2F, "hPtVsDcaz", "p_{T} vs DCA_{Z}; p_{T} (GeV/#it{c}); DCA_{Z} (cm)"},
      {kPtVsDca, o2::framework::kTH2F, "hPtVsDca", "p_{T} vs DCA; p_{T} (GeV/#it{c}); DCA (cm)"},
      {kItsSignal, o2::framework::kTH2F, "hItsSignal", "ITS Signal; p (GeV/#it{c}) ; <ITS Cluster Size> x <cos #lambda>"},
      {kItsElectron, o2::framework::kTH2F, "hItsPidElectron", "ITS PID Electron; p (GeV/#it{c}) ; n#sigma_{TPC,el}"},
      {kItsPion, o2::framework::kTH2F, "hItsPidPion", "ITS PID Pion; p (GeV/#it{c}) ; n#sigma_{ITS,pi}"},
      {kItsKaon, o2::framework::kTH2F, "hItsPidKaon", "ITS PID Kaon; p (GeV/#it{c}) ; n#sigma_{ITS,ka}"},
      {kItsProton, o2::framework::kTH2F, "hItsPidProton", "ITS PID Proton; p (GeV/#it{c}) ; n#sigma_{ITS,pr}"},
      {kItsDeuteron, o2::framework::kTH2F, "hItsPidDeuteron", "ITS PID Deuteron; p (GeV/#it{c}) ; n#sigma_{ITS,de}"},
      {kItsTriton, o2::framework::kTH2F, "hItsPidTriton", "ITS PID Triton; p (GeV/#it{c}) ; n#sigma_{ITS,tr}"},
      {kItsHelium, o2::framework::kTH2F, "hItsPidHelium", "ITS PID Helium; p (GeV/#it{c}) ; n#sigma_{ITS,he}"},
      {kTpcSignal, o2::framework::kTH2F, "hTpcSignal", "TPC Signal; p (GeV/#it{c}) ; TPC Signal"},
      {kTpcElectron, o2::framework::kTH2F, "hTpcPidElectron", "TPC PID Electron; p (GeV/#it{c}) ; n#sigma_{TPC,el}"},
      {kTpcPion, o2::framework::kTH2F, "hTpcPidPion", "TPC PID Pion; p (GeV/#it{c}) ; n#sigma_{TPC,pi}"},
      {kTpcKaon, o2::framework::kTH2F, "hTpcPidKaon", "TPC PID Kaon; p (GeV/#it{c}) ; n#sigma_{TPC,ka}"},
      {kTpcProton, o2::framework::kTH2F, "hTpcPidProton", "TPC PID Proton; p (GeV/#it{c}) ; n#sigma_{TPC,pr}"},
      {kTpcDeuteron, o2::framework::kTH2F, "hTpcPidDeuteron", "TPC PID Deuteron; p (GeV/#it{c}) ; n#sigma_{TPC,de}"},
      {kTpcTriton, o2::framework::kTH2F, "hTpcPidTriton", "TPC PID Triton; p (GeV/#it{c}) ; n#sigma_{TPC,tr}"},
      {kTpcHelium, o2::framework::kTH2F, "hTpcPidHelium", "TPC PID Helium; p (GeV/#it{c}) ; n#sigma_{TPC,he}"},
      {kTofBeta, o2::framework::kTH2F, "hTofBeta", "TOF #beta; p (GeV/#it{c}) ; TOF #beta"},
      {kTofMass, o2::framework::kTH2F, "hTofMass", "TOF mass; p (GeV/#it{c}) ; m_{TOF} (GeV/#it{c}^{2})"},
      {kTofElectron, o2::framework::kTH2F, "hTofPidElectron", "TOF PID Electron; p (GeV/#it{c}) ; n#sigma_{TOF,el}"},
      {kTofPion, o2::framework::kTH2F, "hTofPidPion", "TOF PID Pion; p (GeV/#it{c}) ; n#sigma_{TOF,pi}"},
      {kTofKaon, o2::framework::kTH2F, "hTofPidKaon", "TOF PID Kaon; p (GeV/#it{c}) ; n#sigma_{TOF,ka}"},
      {kTofProton, o2::framework::kTH2F, "hTofPidProton", "TOF PID Proton; p (GeV/#it{c}) ; n#sigma_{TOF,pr}"},
      {kTofDeuteron, o2::framework::kTH2F, "hTofPidDeuteron", "TOF PID Deuteron; p (GeV/#it{c}) ; n#sigma_{TOF,de}"},
      {kTofTriton, o2::framework::kTH2F, "hTofPidTriton", "TOF PID Triton; p (GeV/#it{c}) ; n#sigma_{TOF,tr}"},
      {kTofHelium, o2::framework::kTH2F, "hTofPidHelium", "TOF PID Helium; p (GeV/#it{c}) ; n#sigma_{TOF,he}"},
      {kTpcitsElectron, o2::framework::kTH2F, "hTpcitsPidElectron", "its PID Electron; p (GeV/#it{c}) ; n#sigma_{its,el}"},
      {kTpcitsPion, o2::framework::kTH2F, "hTpcitsPidPion", "TPC+ITS PID Pion; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,pi}^{2}+n#sigma_{its,pi}^{2}}"},
      {kTpcitsKaon, o2::framework::kTH2F, "hTpcitsPidKaon", "TPC+ITS PID Kaon; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,ka}^{2}+n#sigma_{its,ka}^{2}}"},
      {kTpcitsProton, o2::framework::kTH2F, "hTpcitsPidProton", "TPC+ITS PID Proton; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,pr}^{2}+n#sigma_{its,pr}^{2}}"},
      {kTpcitsDeuteron, o2::framework::kTH2F, "hTpcitsPidDeuteron", "TPC+ITS PID Deuteron; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,de}^{2}+n#sigma_{its,de}^{2}}"},
      {kTpcitsTriton, o2::framework::kTH2F, "hTpcitsPidTriton", "TPC+ITS PID Triton; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,tr}^{2}+n#sigma_{its,tr}^{2}}"},
      {kTpcitsHelium, o2::framework::kTH2F, "hTpcitsPidHelium", "TPC+ITS PID Helium; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,he}^{2}+n#sigma_{its,he}^{2}}"},
      {kTpctofElectron, o2::framework::kTH2F, "hTpctofPidElectron", "TOF PID Electron; p (GeV/#it{c}) ; n#sigma_{TOF,el}"},
      {kTpctofPion, o2::framework::kTH2F, "hTpctofPidPion", "TPC+TOF PID Pion; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,pi}^{2}+n#sigma_{TOF,pi}^{2}}"},
      {kTpctofKaon, o2::framework::kTH2F, "hTpctofPidKaon", "TPC+TOF PID Kaon; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,ka}^{2}+n#sigma_{TOF,ka}^{2}}"},
      {kTpctofProton, o2::framework::kTH2F, "hTpctofPidProton", "TPC+TOF PID Proton; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,pr}^{2}+n#sigma_{TOF,pr}^{2}}"},
      {kTpctofDeuteron, o2::framework::kTH2F, "hTpctofPidDeuteron", "TPC+TOF PID Deuteron; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,de}^{2}+n#sigma_{TOF,de}^{2}}"},
      {kTpctofTriton, o2::framework::kTH2F, "hTpctofPidTriton", "TPC+TOF PID Triton; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,tr}^{2}+n#sigma_{TOF,tr}^{2}}"},
      {kTpctofHelium, o2::framework::kTH2F, "hTpctofPidHelium", "TPC+TOF PID Helium; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,he}^{2}+n#sigma_{TOF,he}^{2}}"},
      {kOrigin, o2::framework::kTH1F, "hOrigin", "Status Codes (=Origin); Status Code; Entries"},
      {kPdg, o2::framework::kTH1F, "hPdg", "PDG Codes of selected tracks; PDG Code; Entries"},
      {kPdgMother, o2::framework::kTH1F, "hPdgMother", "PDG Codes of mother of selected tracks; PDG Code; Entries"},
      {kPdgPartonicMother, o2::framework::kTH1F, "hPdgPartonicMother", "PDG Codes of partonic mother selected tracks; PDG Code; Entries"},
      {kTruePtVsPt, o2::framework::kTH2F, "hTruePtVsPt", "True transverse momentum vs transverse momentum; p_{T,True} (GeV/#it{c}); p_{T,True} (GeV/#it{c})"},
      {kTrueEtaVsEta, o2::framework::kTH2F, "hTrueEtaVsEta", "True pseudorapdity vs pseudorapdity; #eta_{True}; #eta"},
      {kTruePhiVsPhi, o2::framework::kTH2F, "hTruePhiVsPhi", "True azimuthal angle vs azimuthal angle; #varphi_{True}; #varphi"},
      {kNoMcParticle, o2::framework::kTHnSparseF, "hNoMcParticle", "Wrongly reconstructed particles; p_{T} (GeV/#it{c}); DCA_{xy} (cm); DCA_{z} (cm)"},
      {kPrimary, o2::framework::kTHnSparseF, "hPrimary", "Primary particles; p_{T} (GeV/#it{c}); DCA_{xy} (cm); DCA_{z} (cm)"},
      {kFromWrongCollision, o2::framework::kTHnSparseF, "hFromWrongCollision", "Particles associated to wrong collision; p_{T} (GeV/#it{c}); DCA_{xy} (cm); DCA_{z} (cm)"},
      {kFromMaterial, o2::framework::kTHnSparseF, "hFromMaterial", "Particles from material; p_{T} (GeV/#it{c}); DCA_{xy} (cm); DCA_{z} (cm)"},
      {kMissidentified, o2::framework::kTHnSparseF, "hMissidentified", "Missidentified particles (fake/wrong PDG code); p_{T} (GeV/#it{c}); DCA_{xy} (cm); DCA_{z} (cm)"},
      {kSecondary1, o2::framework::kTHnSparseF, "hFromSecondary1", "Particles from secondary decay; p_{T} (GeV/#it{c}); DCA_{xy} (cm); DCA_{z} (cm)"},
      {kSecondary2, o2::framework::kTHnSparseF, "hFromSecondary2", "Particles from seconary decay; p_{T} (GeV/#it{c}); DCA_{xy} (cm); DCA_{z} (cm)"},
      {kSecondary3, o2::framework::kTHnSparseF, "hFromSecondary3", "Particles from seconary decay; p_{T} (GeV/#it{c}); DCA_{xy} (cm); DCA_{z} (cm)"},
      {kSecondaryOther, o2::framework::kTHnSparseF, "hFromSecondaryOther", "Particles from every other seconary decay; p_{T} (GeV/#it{c}); DCA_{xy} (cm); DCA_{z} (cm)"},
    }};

#define TRACK_HIST_ANALYSIS_MAP(conf) \
  {kPt, {conf.pt}},                   \
    {kEta, {conf.eta}},               \
    {kPhi, {conf.phi}},               \
    {kSign, {conf.sign}},

#define TRACK_HIST_QA_MAP(confAnalysis, confQa)                                            \
  {kPAtPv, {confQa.p}},                                                                    \
    {kPTpc, {confQa.p}},                                                                   \
    {kItsCluster, {confQa.itsCluster}},                                                    \
    {kItsClusterIb, {confQa.itsClusterIb}},                                                \
    {kPtVsEta, {confAnalysis.pt, confAnalysis.eta}},                                       \
    {kPtVsPhi, {confAnalysis.pt, confAnalysis.phi}},                                       \
    {kPhiVsEta, {confAnalysis.phi, confAnalysis.eta}},                                     \
    {kPtVsItsCluster, {confAnalysis.pt, confQa.itsCluster}},                               \
    {kPtVsTpcCluster, {confAnalysis.pt, confQa.tpcCluster}},                               \
    {kPtVsTpcCrossedRows, {confAnalysis.pt, confQa.tpcCrossedRows}},                       \
    {kPtVsTpcClusterOverCrossedRows, {confAnalysis.pt, confQa.tpcClusterOverCrossedRows}}, \
    {kPtVsTpcClusterShared, {confAnalysis.pt, confQa.tpcClusterShared}},                   \
    {kPtVsTpcClusterFractionShared, {confAnalysis.pt, confQa.tpcClusterFractionShared}},   \
    {kTpcClusterVsTpcCrossedRows, {confQa.tpcCluster, confQa.tpcCrossedRows}},             \
    {kTpcClusterVsTpcClusterShared, {confQa.tpcCluster, confQa.tpcClusterShared}},         \
    {kTpcCrossedRows, {confQa.tpcCrossedRows}},                                            \
    {kTpcCluster, {confQa.tpcCluster}},                                                    \
    {kTpcClusterOverCrossedRows, {confQa.tpcClusterOverCrossedRows}},                      \
    {kTpcClusterShared, {confQa.tpcClusterShared}},                                        \
    {kTpcClusterFractionShared, {confQa.tpcClusterFractionShared}},                        \
    {kPtVsDcaxy, {confAnalysis.pt, confQa.dcaXy}},                                         \
    {kPtVsDcaz, {confAnalysis.pt, confQa.dcaZ}},                                           \
    {kPtVsDca, {confAnalysis.pt, confQa.dca}},                                             \
    {kItsSignal, {confQa.p, confQa.itsSignal}},                                            \
    {kItsElectron, {confQa.p, confQa.itsElectron}},                                        \
    {kItsPion, {confQa.p, confQa.itsPion}},                                                \
    {kItsKaon, {confQa.p, confQa.itsKaon}},                                                \
    {kItsProton, {confQa.p, confQa.itsProton}},                                            \
    {kItsDeuteron, {confQa.p, confQa.itsDeuteron}},                                        \
    {kItsTriton, {confQa.p, confQa.itsTriton}},                                            \
    {kItsHelium, {confQa.p, confQa.itsHelium}},                                            \
    {kTpcSignal, {confQa.p, confQa.tpcSignal}},                                            \
    {kTpcElectron, {confQa.p, confQa.tpcElectron}},                                        \
    {kTpcPion, {confQa.p, confQa.tpcPion}},                                                \
    {kTpcKaon, {confQa.p, confQa.tpcKaon}},                                                \
    {kTpcProton, {confQa.p, confQa.tpcProton}},                                            \
    {kTpcDeuteron, {confQa.p, confQa.tpcDeuteron}},                                        \
    {kTpcTriton, {confQa.p, confQa.tpcTriton}},                                            \
    {kTpcHelium, {confQa.p, confQa.tpcHelium}},                                            \
    {kTofBeta, {confQa.p, confQa.tofBeta}},                                                \
    {kTofMass, {confQa.p, confQa.tofMass}},                                                \
    {kTofElectron, {confQa.p, confQa.tofElectron}},                                        \
    {kTofPion, {confQa.p, confQa.tofPion}},                                                \
    {kTofKaon, {confQa.p, confQa.tofKaon}},                                                \
    {kTofProton, {confQa.p, confQa.tofProton}},                                            \
    {kTofDeuteron, {confQa.p, confQa.tofDeuteron}},                                        \
    {kTofTriton, {confQa.p, confQa.tofTriton}},                                            \
    {kTofHelium, {confQa.p, confQa.tofHelium}},                                            \
    {kTpcitsElectron, {confQa.p, confQa.tpcitsElectron}},                                  \
    {kTpcitsPion, {confQa.p, confQa.tpcitsPion}},                                          \
    {kTpcitsKaon, {confQa.p, confQa.tpcitsKaon}},                                          \
    {kTpcitsProton, {confQa.p, confQa.tpcitsProton}},                                      \
    {kTpcitsDeuteron, {confQa.p, confQa.tpcitsDeuteron}},                                  \
    {kTpcitsTriton, {confQa.p, confQa.tpcitsTriton}},                                      \
    {kTpcitsHelium, {confQa.p, confQa.tpcitsHelium}},                                      \
    {kTpctofElectron, {confQa.p, confQa.tpctofElectron}},                                  \
    {kTpctofPion, {confQa.p, confQa.tpctofPion}},                                          \
    {kTpctofKaon, {confQa.p, confQa.tpctofKaon}},                                          \
    {kTpctofProton, {confQa.p, confQa.tpctofProton}},                                      \
    {kTpctofDeuteron, {confQa.p, confQa.tpctofDeuteron}},                                  \
    {kTpctofTriton, {confQa.p, confQa.tpctofTriton}},                                      \
    {kTpctofHelium, {confQa.p, confQa.tpctofHelium}},

#define TRACK_HIST_MC_MAP(conf)            \
  {kTruePtVsPt, {conf.pt, conf.pt}},       \
    {kTrueEtaVsEta, {conf.eta, conf.eta}}, \
    {kTruePhiVsPhi, {conf.phi, conf.phi}}, \
    {kPdg, {conf.pdgCodes}},               \
    {kPdgMother, {conf.pdgCodes}},         \
    {kPdgPartonicMother, {conf.pdgCodes}},

#define TRACK_HIST_MC_QA_MAP(confAnalysis, confQa)                       \
  {kNoMcParticle, {confAnalysis.pt, confQa.dcaXy, confQa.dcaZ}},         \
    {kPrimary, {confAnalysis.pt, confQa.dcaXy, confQa.dcaZ}},            \
    {kFromWrongCollision, {confAnalysis.pt, confQa.dcaXy, confQa.dcaZ}}, \
    {kFromMaterial, {confAnalysis.pt, confQa.dcaXy, confQa.dcaZ}},       \
    {kMissidentified, {confAnalysis.pt, confQa.dcaXy, confQa.dcaZ}},     \
    {kSecondary1, {confAnalysis.pt, confQa.dcaXy, confQa.dcaZ}},         \
    {kSecondary2, {confAnalysis.pt, confQa.dcaXy, confQa.dcaZ}},         \
    {kSecondary3, {confAnalysis.pt, confQa.dcaXy, confQa.dcaZ}},         \
    {kSecondaryOther, {confAnalysis.pt, confQa.dcaXy, confQa.dcaZ}},

template <typename T>
auto makeTrackHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<TrackHist, std::vector<framework::AxisSpec>>{
    TRACK_HIST_ANALYSIS_MAP(confBinningAnalysis)};
}

template <typename T>
auto makeTrackMcHistSpecMap(T const& confBinningAnalysis)
{
  return std::map<TrackHist, std::vector<framework::AxisSpec>>{
    TRACK_HIST_ANALYSIS_MAP(confBinningAnalysis)
      TRACK_HIST_MC_MAP(confBinningAnalysis)};
};

template <typename T1, typename T2>
auto makeTrackQaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<TrackHist, std::vector<framework::AxisSpec>>{
    TRACK_HIST_ANALYSIS_MAP(confBinningAnalysis)
      TRACK_HIST_QA_MAP(confBinningAnalysis, confBinningQa)};
}

template <typename T1, typename T2>
auto makeTrackMcQaHistSpecMap(T1 const& confBinningAnalysis, T2 const& confBinningQa)
{
  return std::map<TrackHist, std::vector<framework::AxisSpec>>{
    TRACK_HIST_ANALYSIS_MAP(confBinningAnalysis)
      TRACK_HIST_QA_MAP(confBinningAnalysis, confBinningQa)
        TRACK_HIST_MC_MAP(confBinningAnalysis)
          TRACK_HIST_MC_QA_MAP(confBinningAnalysis, confBinningQa)};
};

#undef TRACK_HIST_ANALYSIS_MAP
#undef TRACK_HIST_QA_MAP
#undef TRACK_HIST_MC_MAP
#undef TRACK_HIST_MC_QA_MAP

constexpr char PrefixTrackQa[] = "TrackQA/";
constexpr char PrefixTrack1[] = "Track1/";
constexpr char PrefixTrack2[] = "Track2/";
constexpr char PrefixTrack3[] = "Track3/";

constexpr char PrefixResonancePosDaughter[] = "ResonancePosDau/";
constexpr char PrefixResonanceNegDaughter[] = "ResonanceNegDau/";
constexpr char PrefixResonancePosDaughterQa[] = "ResonancePosDauQa/";
constexpr char PrefixResonanceNegDaughterQa[] = "ResonanceNegDauQa/";

constexpr char PrefixV01PosDaughter[] = "V01PosDau/";
constexpr char PrefixV01NegDaughter[] = "V01NegDau/";
constexpr char PrefixV02PosDaughter[] = "V02PosDau/";
constexpr char PrefixV02NegDaughter[] = "V02NegDau/";
constexpr char PrefixV0PosDaughterQa[] = "V0PosDauQa/";
constexpr char PrefixV0NegDaughterQa[] = "V0NegDauQa/";

constexpr char PrefixCascadePosDaughter[] = "CascadePosDau/";
constexpr char PrefixCascadeNegDaughter[] = "CascadeNegDau/";
constexpr char PrefixCascadeBachelor[] = "CascadeBachelor/";
constexpr char PrefixCascadePosDaughterQa[] = "CascadePosDauQa/";
constexpr char PrefixCascadeNegDaughterQa[] = "CascadeNegDauQa/";
constexpr char PrefixCascadeBachelorQa[] = "CascadeBachelorQa/";

constexpr char PrefixKinkChaDaughter[] = "KinkChaDau/";
constexpr char PrefixKinkChaDaughterQa[] = "KinkChaDauQa/";

constexpr std::string_view AnalysisDir = "Kinematics/";
constexpr std::string_view QaDir = "QA/";
constexpr std::string_view PidDir = "PID/";
constexpr std::string_view McDir = "MC/";

template <const char* prefix>
class TrackHistManager
{
 public:
  TrackHistManager() = default;
  ~TrackHistManager() = default;

  // init for analysis and mc
  template <modes::Mode mode, typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<TrackHist, std::vector<o2::framework::AxisSpec>> const& Specs,
            T const& ConfTrackSelection)
  {
    mHistogramRegistry = registry;
    mAbsCharge = std::abs(ConfTrackSelection.chargeAbs.value);
    mPdgCode = std::abs(ConfTrackSelection.pdgCodeAbs.value) * ConfTrackSelection.chargeSign.value;
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      this->initAnalysis(Specs);
    }
    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      this->initQa(Specs);
    }
    if constexpr (isFlagSet(mode, modes::Mode::kMc)) {
      this->initMc(Specs);
    }
  }

  template <modes::Mode mode>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<TrackHist, std::vector<o2::framework::AxisSpec>> const& Specs,
            int ChargeAbs,
            int ChargeSign,
            int PdgCodeAbs)
  {
    mHistogramRegistry = registry;
    mAbsCharge = std::abs(ChargeAbs);
    mPdgCode = std::abs(PdgCodeAbs) * ChargeSign;
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      this->initAnalysis(Specs);
    }
    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      this->initQa(Specs);
    }
    if constexpr (isFlagSet(mode, modes::Mode::kMc)) {
      this->initMc(Specs);
    }
  }

  // init for analysis and qa and mc
  template <modes::Mode mode, typename T1, typename T2>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<TrackHist, std::vector<o2::framework::AxisSpec>> const& Specs,
            T1 const& ConfTrackSelection,
            T2 const& ConfBinningQa)
  {
    this->template enableOptionalHistograms<T2>(ConfBinningQa);
    this->template init<mode>(registry, Specs, ConfTrackSelection);
  }

  template <modes::Mode mode, typename T>
  void init(o2::framework::HistogramRegistry* registry,
            std::map<TrackHist, std::vector<o2::framework::AxisSpec>> const& Specs,
            int ChargeAbs,
            int ChargeSign,
            int PdgCodeAbs,
            T const& ConfBinningQa)
  {
    this->template enableOptionalHistograms<T>(ConfBinningQa);
    this->template init<mode>(registry, Specs, ChargeAbs, ChargeSign, PdgCodeAbs);
  }

  template <modes::Mode mode, typename T1, typename T2>
  void fill(T1 const& track, T2 const& /*trackTable*/)
  {
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      this->fillAnalysis(track);
    }
    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      this->fillQa(track);
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4, typename T5>
  void fill(T1 const& track, T2 const& /*trackTable*/, T3 const& mcParticles, T4 const& mcMothers, T5 const& mcPartonicMothers)
  {
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      this->fillAnalysis(track);
    }
    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      this->fillQa(track);
    }
    if constexpr (isFlagSet(mode, modes::Mode::kMc)) {
      this->template fillMc<mode>(track, mcParticles, mcMothers, mcPartonicMothers);
    }
  }

 private:
  template <typename T>
  void enableOptionalHistograms(T const& ConfBinningQa)
  {
    mPlot2d = ConfBinningQa.plot2d.value;
    mPlotElectronPid = ConfBinningQa.plotElectronPid.value;
    mPlotPionPid = ConfBinningQa.plotPionPid.value;
    mPlotKaonPid = ConfBinningQa.plotKaonPid.value;
    mPlotProtonPid = ConfBinningQa.plotProtonPid.value;
    mPlotDeuteronPid = ConfBinningQa.plotDeuteronPid.value;
    mPlotTritonPid = ConfBinningQa.plotTritonPid.value;
    mPlotHeliumPid = ConfBinningQa.plotHeliumPid.value;
    mMomentumType = static_cast<modes::MomentumType>(ConfBinningQa.momentumType.value);

    mPlotOrigins = ConfBinningQa.plotOrigins.value;
    mPlotNSecondaries = ConfBinningQa.pdgCodesForMothersOfSecondary.value.size();

    for (std::size_t i = 0; i < MaxSecondary; i++) {
      if (i < ConfBinningQa.pdgCodesForMothersOfSecondary.value.size()) {
        mPdgCodesSecondaryMother.at(i) = std::abs(ConfBinningQa.pdgCodesForMothersOfSecondary.value.at(i));
      } else {
        mPdgCodesSecondaryMother.at(i) = 0;
      }
    }
  }

  void initAnalysis(std::map<TrackHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    std::string analysisDir = std::string(prefix) + std::string(AnalysisDir);
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPt, HistTable), getHistDesc(kPt, HistTable), getHistType(kPt, HistTable), {Specs.at(kPt)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kEta, HistTable), getHistDesc(kEta, HistTable), getHistType(kEta, HistTable), {Specs.at(kEta)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kPhi, HistTable), getHistDesc(kPhi, HistTable), getHistType(kPhi, HistTable), {Specs.at(kPhi)});
    mHistogramRegistry->add(analysisDir + getHistNameV2(kSign, HistTable), getHistDesc(kSign, HistTable), getHistType(kSign, HistTable), {Specs.at(kSign)});
  }

  void initQa(std::map<TrackHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    std::string qaDir = std::string(prefix) + std::string(QaDir);

    mHistogramRegistry->add(qaDir + getHistNameV2(kPAtPv, HistTable), getHistDesc(kPAtPv, HistTable), getHistType(kPAtPv, HistTable), {Specs.at(kPAtPv)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kPTpc, HistTable), getHistDesc(kPTpc, HistTable), getHistType(kPTpc, HistTable), {Specs.at(kPTpc)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kItsCluster, HistTable), getHistDesc(kItsCluster, HistTable), getHistType(kItsCluster, HistTable), {Specs.at(kItsCluster)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kItsClusterIb, HistTable), getHistDesc(kItsClusterIb, HistTable), getHistType(kItsClusterIb, HistTable), {Specs.at(kItsClusterIb)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kTpcCrossedRows, HistTable), getHistDesc(kTpcCrossedRows, HistTable), getHistType(kTpcCrossedRows, HistTable), {Specs.at(kTpcCrossedRows)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kTpcCluster, HistTable), getHistDesc(kTpcCluster, HistTable), getHistType(kTpcCluster, HistTable), {Specs.at(kTpcCluster)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kTpcClusterOverCrossedRows, HistTable), getHistDesc(kTpcClusterOverCrossedRows, HistTable), getHistType(kTpcClusterOverCrossedRows, HistTable), {Specs.at(kTpcClusterOverCrossedRows)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kTpcClusterShared, HistTable), getHistDesc(kTpcClusterShared, HistTable), getHistType(kTpcClusterShared, HistTable), {Specs.at(kTpcClusterShared)});
    mHistogramRegistry->add(qaDir + getHistNameV2(kTpcClusterFractionShared, HistTable), getHistDesc(kTpcClusterFractionShared, HistTable), getHistType(kTpcClusterFractionShared, HistTable), {Specs.at(kTpcClusterFractionShared)});

    // qa 2d
    if (mPlot2d) {
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsEta, HistTable), getHistDesc(kPtVsEta, HistTable), getHistType(kPtVsEta, HistTable), {Specs.at(kPtVsEta)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsPhi, HistTable), getHistDesc(kPtVsPhi, HistTable), getHistType(kPtVsPhi, HistTable), {Specs.at(kPtVsPhi)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPhiVsEta, HistTable), getHistDesc(kPhiVsEta, HistTable), getHistType(kPhiVsEta, HistTable), {Specs.at(kPhiVsEta)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsItsCluster, HistTable), getHistDesc(kPtVsItsCluster, HistTable), getHistType(kPtVsItsCluster, HistTable), {Specs.at(kPtVsItsCluster)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsTpcCluster, HistTable), getHistDesc(kPtVsTpcCluster, HistTable), getHistType(kPtVsTpcCluster, HistTable), {Specs.at(kPtVsTpcCluster)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsTpcCrossedRows, HistTable), getHistDesc(kPtVsTpcCrossedRows, HistTable), getHistType(kPtVsTpcCrossedRows, HistTable), {Specs.at(kPtVsTpcCrossedRows)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsTpcClusterOverCrossedRows, HistTable), getHistDesc(kPtVsTpcClusterOverCrossedRows, HistTable), getHistType(kPtVsTpcClusterOverCrossedRows, HistTable), {Specs.at(kPtVsTpcClusterOverCrossedRows)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsTpcClusterShared, HistTable), getHistDesc(kPtVsTpcClusterShared, HistTable), getHistType(kPtVsTpcClusterShared, HistTable), {Specs.at(kPtVsTpcClusterShared)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsTpcClusterFractionShared, HistTable), getHistDesc(kPtVsTpcClusterFractionShared, HistTable), getHistType(kPtVsTpcClusterFractionShared, HistTable), {Specs.at(kPtVsTpcClusterFractionShared)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kTpcClusterVsTpcCrossedRows, HistTable), getHistDesc(kTpcClusterVsTpcCrossedRows, HistTable), getHistType(kTpcClusterVsTpcCrossedRows, HistTable), {Specs.at(kTpcClusterVsTpcCrossedRows)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kTpcClusterVsTpcClusterShared, HistTable), getHistDesc(kTpcClusterVsTpcClusterShared, HistTable), getHistType(kTpcClusterVsTpcClusterShared, HistTable), {Specs.at(kTpcClusterVsTpcClusterShared)});
      // dca
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsDcaxy, HistTable), getHistDesc(kPtVsDcaxy, HistTable), getHistType(kPtVsDcaxy, HistTable), {Specs.at(kPtVsDcaxy)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsDcaz, HistTable), getHistDesc(kPtVsDcaz, HistTable), getHistType(kPtVsDcaz, HistTable), {Specs.at(kPtVsDcaz)});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsDca, HistTable), getHistDesc(kPtVsDca, HistTable), getHistType(kPtVsDca, HistTable), {Specs.at(kPtVsDca)});
    }

    std::string pidDir = std::string(prefix) + std::string(PidDir);

    mHistogramRegistry->add(pidDir + getHistNameV2(kItsSignal, HistTable), getHistDesc(kItsSignal, HistTable), getHistType(kItsSignal, HistTable), {Specs.at(kItsSignal)});
    mHistogramRegistry->add(pidDir + getHistNameV2(kTpcSignal, HistTable), getHistDesc(kTpcSignal, HistTable), getHistType(kTpcSignal, HistTable), {Specs.at(kTpcSignal)});
    mHistogramRegistry->add(pidDir + getHistNameV2(kTofBeta, HistTable), getHistDesc(kTofBeta, HistTable), getHistType(kTofBeta, HistTable), {Specs.at(kTofBeta)});
    mHistogramRegistry->add(pidDir + getHistNameV2(kTofMass, HistTable), getHistDesc(kTofMass, HistTable), getHistType(kTofMass, HistTable), {Specs.at(kTofMass)});

    if (mPlotElectronPid) {
      mHistogramRegistry->add(pidDir + getHistNameV2(kItsElectron, HistTable), getHistDesc(kItsElectron, HistTable), getHistType(kItsElectron, HistTable), {Specs.at(kItsElectron)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcElectron, HistTable), getHistDesc(kTpcElectron, HistTable), getHistType(kTpcElectron, HistTable), {Specs.at(kTpcElectron)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTofElectron, HistTable), getHistDesc(kTofElectron, HistTable), getHistType(kTofElectron, HistTable), {Specs.at(kTofElectron)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcitsElectron, HistTable), getHistDesc(kTpcitsElectron, HistTable), getHistType(kTpcitsElectron, HistTable), {Specs.at(kTpcitsElectron)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpctofElectron, HistTable), getHistDesc(kTpctofElectron, HistTable), getHistType(kTpctofElectron, HistTable), {Specs.at(kTpctofElectron)});
    }

    if (mPlotPionPid) {
      mHistogramRegistry->add(pidDir + getHistNameV2(kItsPion, HistTable), getHistDesc(kItsPion, HistTable), getHistType(kItsPion, HistTable), {Specs.at(kItsPion)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcPion, HistTable), getHistDesc(kTpcPion, HistTable), getHistType(kTpcPion, HistTable), {Specs.at(kTpcPion)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTofPion, HistTable), getHistDesc(kTofPion, HistTable), getHistType(kTofPion, HistTable), {Specs.at(kTofPion)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcitsPion, HistTable), getHistDesc(kTpcitsPion, HistTable), getHistType(kTpcitsPion, HistTable), {Specs.at(kTpcitsPion)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpctofPion, HistTable), getHistDesc(kTpctofPion, HistTable), getHistType(kTpctofPion, HistTable), {Specs.at(kTpctofPion)});
    }

    if (mPlotKaonPid) {
      mHistogramRegistry->add(pidDir + getHistNameV2(kItsKaon, HistTable), getHistDesc(kItsKaon, HistTable), getHistType(kItsKaon, HistTable), {Specs.at(kItsKaon)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcKaon, HistTable), getHistDesc(kTpcKaon, HistTable), getHistType(kTpcKaon, HistTable), {Specs.at(kTpcKaon)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTofKaon, HistTable), getHistDesc(kTofKaon, HistTable), getHistType(kTofKaon, HistTable), {Specs.at(kTofKaon)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcitsKaon, HistTable), getHistDesc(kTpcitsKaon, HistTable), getHistType(kTpcitsKaon, HistTable), {Specs.at(kTpcitsKaon)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpctofKaon, HistTable), getHistDesc(kTpctofKaon, HistTable), getHistType(kTpctofKaon, HistTable), {Specs.at(kTpctofKaon)});
    }

    if (mPlotProtonPid) {
      mHistogramRegistry->add(pidDir + getHistNameV2(kItsProton, HistTable), getHistDesc(kItsProton, HistTable), getHistType(kItsProton, HistTable), {Specs.at(kItsProton)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcProton, HistTable), getHistDesc(kTpcProton, HistTable), getHistType(kTpcProton, HistTable), {Specs.at(kTpcProton)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTofProton, HistTable), getHistDesc(kTofProton, HistTable), getHistType(kTofProton, HistTable), {Specs.at(kTofProton)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcitsProton, HistTable), getHistDesc(kTpcitsProton, HistTable), getHistType(kTpcitsProton, HistTable), {Specs.at(kTpcitsProton)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpctofProton, HistTable), getHistDesc(kTpctofProton, HistTable), getHistType(kTpctofProton, HistTable), {Specs.at(kTpctofProton)});
    }

    if (mPlotDeuteronPid) {
      mHistogramRegistry->add(pidDir + getHistNameV2(kItsDeuteron, HistTable), getHistDesc(kItsDeuteron, HistTable), getHistType(kItsDeuteron, HistTable), {Specs.at(kItsDeuteron)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcDeuteron, HistTable), getHistDesc(kTpcDeuteron, HistTable), getHistType(kTpcDeuteron, HistTable), {Specs.at(kTpcDeuteron)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTofDeuteron, HistTable), getHistDesc(kTofDeuteron, HistTable), getHistType(kTofDeuteron, HistTable), {Specs.at(kTofDeuteron)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcitsDeuteron, HistTable), getHistDesc(kTpcitsDeuteron, HistTable), getHistType(kTpcitsDeuteron, HistTable), {Specs.at(kTpcitsDeuteron)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpctofDeuteron, HistTable), getHistDesc(kTpctofDeuteron, HistTable), getHistType(kTpctofDeuteron, HistTable), {Specs.at(kTpctofDeuteron)});
    }

    if (mPlotTritonPid) {
      mHistogramRegistry->add(pidDir + getHistNameV2(kItsTriton, HistTable), getHistDesc(kItsTriton, HistTable), getHistType(kItsTriton, HistTable), {Specs.at(kItsTriton)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcTriton, HistTable), getHistDesc(kTpcTriton, HistTable), getHistType(kTpcTriton, HistTable), {Specs.at(kTpcTriton)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTofTriton, HistTable), getHistDesc(kTofTriton, HistTable), getHistType(kTofTriton, HistTable), {Specs.at(kTofTriton)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcitsTriton, HistTable), getHistDesc(kTpcitsTriton, HistTable), getHistType(kTpcitsTriton, HistTable), {Specs.at(kTpcitsTriton)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpctofTriton, HistTable), getHistDesc(kTpctofTriton, HistTable), getHistType(kTpctofTriton, HistTable), {Specs.at(kTpctofTriton)});
    }

    if (mPlotHeliumPid) {
      mHistogramRegistry->add(pidDir + getHistNameV2(kItsHelium, HistTable), getHistDesc(kItsHelium, HistTable), getHistType(kItsHelium, HistTable), {Specs.at(kItsHelium)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcHelium, HistTable), getHistDesc(kTpcHelium, HistTable), getHistType(kTpcHelium, HistTable), {Specs.at(kTpcHelium)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTofHelium, HistTable), getHistDesc(kTofHelium, HistTable), getHistType(kTofHelium, HistTable), {Specs.at(kTofHelium)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcitsHelium, HistTable), getHistDesc(kTpcitsHelium, HistTable), getHistType(kTpcitsHelium, HistTable), {Specs.at(kTpcitsHelium)});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpctofHelium, HistTable), getHistDesc(kTpctofHelium, HistTable), getHistType(kTpctofHelium, HistTable), {Specs.at(kTpctofHelium)});
    }
  }

  void initMc(std::map<TrackHist, std::vector<o2::framework::AxisSpec>> const& Specs)
  {
    std::string mcDir = std::string(prefix) + std::string(McDir);
    mHistogramRegistry->add(mcDir + getHistNameV2(kTruePtVsPt, HistTable), getHistDesc(kTruePtVsPt, HistTable), getHistType(kTruePtVsPt, HistTable), {Specs.at(kTruePtVsPt)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTrueEtaVsEta, HistTable), getHistDesc(kTrueEtaVsEta, HistTable), getHistType(kTrueEtaVsEta, HistTable), {Specs.at(kTrueEtaVsEta)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kTruePhiVsPhi, HistTable), getHistDesc(kTruePhiVsPhi, HistTable), getHistType(kTruePhiVsPhi, HistTable), {Specs.at(kTruePhiVsPhi)});

    // mc origin can be configured here
    const framework::AxisSpec axisOrigin = {static_cast<int>(modes::McOrigin::kMcOriginLast), -0.5, static_cast<double>(modes::McOrigin::kMcOriginLast) - 0.5};
    mHistogramRegistry->add(mcDir + getHistNameV2(kOrigin, HistTable), getHistDesc(kOrigin, HistTable), getHistType(kOrigin, HistTable), {axisOrigin});
    mHistogramRegistry->get<TH1>(HIST(prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kNoMcParticle), modes::mcOriginToString(modes::McOrigin::kNoMcParticle));
    mHistogramRegistry->get<TH1>(HIST(prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kFromWrongCollision), modes::mcOriginToString(modes::McOrigin::kFromWrongCollision));
    mHistogramRegistry->get<TH1>(HIST(prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kPhysicalPrimary), modes::mcOriginToString(modes::McOrigin::kPhysicalPrimary));
    mHistogramRegistry->get<TH1>(HIST(prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kFromSecondaryDecay), modes::mcOriginToString(modes::McOrigin::kFromSecondaryDecay));
    mHistogramRegistry->get<TH1>(HIST(prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kFromMaterial), modes::mcOriginToString(modes::McOrigin::kFromMaterial));
    mHistogramRegistry->get<TH1>(HIST(prefix) + HIST(McDir) + HIST(histmanager::getHistName(kOrigin, HistTable)))->GetXaxis()->SetBinLabel(1 + static_cast<int>(modes::McOrigin::kMissidentified), modes::mcOriginToString(modes::McOrigin::kMissidentified));

    mHistogramRegistry->add(mcDir + getHistNameV2(kPdg, HistTable), getHistDesc(kPdg, HistTable), getHistType(kPdg, HistTable), {Specs.at(kPdg)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdgMother, HistTable), getHistDesc(kPdgMother, HistTable), getHistType(kPdgMother, HistTable), {Specs.at(kPdgMother)});
    mHistogramRegistry->add(mcDir + getHistNameV2(kPdgPartonicMother, HistTable), getHistDesc(kPdgPartonicMother, HistTable), getHistType(kPdgPartonicMother, HistTable), {Specs.at(kPdgPartonicMother)});

    if (mPlotOrigins) {
      mHistogramRegistry->add(mcDir + getHistNameV2(kNoMcParticle, HistTable), getHistDesc(kNoMcParticle, HistTable), getHistType(kNoMcParticle, HistTable), {Specs.at(kNoMcParticle)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kPrimary, HistTable), getHistDesc(kPrimary, HistTable), getHistType(kPrimary, HistTable), {Specs.at(kPrimary)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kFromWrongCollision, HistTable), getHistDesc(kFromWrongCollision, HistTable), getHistType(kFromWrongCollision, HistTable), {Specs.at(kFromWrongCollision)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kFromMaterial, HistTable), getHistDesc(kFromMaterial, HistTable), getHistType(kFromMaterial, HistTable), {Specs.at(kFromMaterial)});
      mHistogramRegistry->add(mcDir + getHistNameV2(kMissidentified, HistTable), getHistDesc(kMissidentified, HistTable), getHistType(kMissidentified, HistTable), {Specs.at(kMissidentified)});

      if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel1) {
        mHistogramRegistry->add(mcDir + getHistNameV2(kSecondary1, HistTable), getHistDesc(kSecondary1, HistTable), getHistType(kSecondary1, HistTable), {Specs.at(kSecondary1)});
      }
      if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel2) {
        mHistogramRegistry->add(mcDir + getHistNameV2(kSecondary2, HistTable), getHistDesc(kSecondary2, HistTable), getHistType(kSecondary2, HistTable), {Specs.at(kSecondary2)});
      }
      if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel3) {
        mHistogramRegistry->add(mcDir + getHistNameV2(kSecondary3, HistTable), getHistDesc(kSecondary3, HistTable), getHistType(kSecondary3, HistTable), {Specs.at(kSecondary3)});
      }
      mHistogramRegistry->add(mcDir + getHistNameV2(kSecondaryOther, HistTable), getHistDesc(kSecondaryOther, HistTable), getHistType(kSecondaryOther, HistTable), {Specs.at(kSecondaryOther)});
    }
  }

  template <typename T>
  void fillAnalysis(T const& track)
  {
    mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt, HistTable)), mAbsCharge * track.pt());
    mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kEta, HistTable)), track.eta());
    mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPhi, HistTable)), track.phi());
    mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kSign, HistTable)), track.sign());
  }

  template <typename T>
  void fillQa(T const& track)
  {
    mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kPAtPv, HistTable)), track.p());
    mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kPTpc, HistTable)), track.tpcInnerParam());
    mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kItsCluster, HistTable)), static_cast<float>(track.itsNCls()));
    mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kItsClusterIb, HistTable)), static_cast<float>(track.itsNClsInnerBarrel()));
    mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kTpcCrossedRows, HistTable)), static_cast<float>(track.tpcNClsCrossedRows()));
    mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kTpcCluster, HistTable)), static_cast<float>(track.tpcNClsFound()));
    mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kTpcClusterOverCrossedRows, HistTable)), static_cast<float>(track.tpcNClsFound()) / static_cast<float>(track.tpcNClsCrossedRows()));
    mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kTpcClusterShared, HistTable)), static_cast<float>(track.tpcNClsShared()));
    mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kTpcClusterFractionShared, HistTable)), track.tpcSharedOverFound());

    if (mPlot2d) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kPtVsEta, HistTable)), mAbsCharge * track.pt(), track.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kPtVsPhi, HistTable)), mAbsCharge * track.pt(), track.phi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kPhiVsEta, HistTable)), track.phi(), track.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kPtVsItsCluster, HistTable)), mAbsCharge * track.pt(), static_cast<float>(track.itsNCls()));
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kPtVsTpcCluster, HistTable)), mAbsCharge * track.pt(), static_cast<float>(track.tpcNClsFound()));
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kPtVsTpcCrossedRows, HistTable)), mAbsCharge * track.pt(), static_cast<float>(track.tpcNClsCrossedRows()));
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kPtVsTpcClusterOverCrossedRows, HistTable)), mAbsCharge * track.pt(), static_cast<float>(track.tpcNClsFound()) / static_cast<float>(track.tpcNClsCrossedRows()));
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kPtVsTpcClusterShared, HistTable)), mAbsCharge * track.pt(), static_cast<float>(track.tpcNClsShared()));
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kPtVsTpcClusterFractionShared, HistTable)), mAbsCharge * track.pt(), static_cast<float>(track.tpcNClsShared()) / static_cast<float>(track.tpcNClsFound()));
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kTpcClusterVsTpcCrossedRows, HistTable)), static_cast<float>(track.tpcNClsFound()), static_cast<float>(track.tpcNClsCrossedRows()));
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kTpcClusterVsTpcClusterShared, HistTable)), static_cast<float>(track.tpcNClsFound()), static_cast<float>(track.tpcNClsShared()));
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kPtVsDcaxy, HistTable)), mAbsCharge * track.pt(), track.dcaXY());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kPtVsDcaz, HistTable)), mAbsCharge * track.pt(), track.dcaZ());
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kPtVsDca, HistTable)), mAbsCharge * track.pt(), track.dca());
    }

    float momentum = 0.f;
    if (mMomentumType == modes::MomentumType::kPt) {
      momentum = mAbsCharge * track.pt();
    } else if (mMomentumType == modes::MomentumType::kPAtPv) {
      momentum = mAbsCharge * track.p();
    } else if (mMomentumType == modes::MomentumType::kPTpc) {
      momentum = track.tpcInnerParam();
    } else {
      LOG(warn) << "Invalid momentum type for PID plots";
      momentum = 0;
    }

    mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kItsSignal, HistTable)), momentum, o2::analysis::femto::utils::itsSignal(track));
    mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcSignal, HistTable)), momentum, track.tpcSignal());
    mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofBeta, HistTable)), momentum, track.tofBeta());
    mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofMass, HistTable)), momentum, track.tofMass());

    if (mPlotElectronPid) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kItsElectron, HistTable)), momentum, track.itsNSigmaEl());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcElectron, HistTable)), momentum, track.tpcNSigmaEl());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofElectron, HistTable)), momentum, track.tofNSigmaEl());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcitsElectron, HistTable)), momentum, track.tpcitsNSigmaEl());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpctofElectron, HistTable)), momentum, track.tpctofNSigmaEl());
    }

    if (mPlotPionPid) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kItsPion, HistTable)), momentum, track.itsNSigmaPi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcPion, HistTable)), momentum, track.tpcNSigmaPi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofPion, HistTable)), momentum, track.tofNSigmaPi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcitsPion, HistTable)), momentum, track.tpcitsNSigmaPi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpctofPion, HistTable)), momentum, track.tpctofNSigmaPi());
    }

    if (mPlotKaonPid) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kItsKaon, HistTable)), momentum, track.itsNSigmaKa());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcKaon, HistTable)), momentum, track.tpcNSigmaKa());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofKaon, HistTable)), momentum, track.tofNSigmaKa());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcitsKaon, HistTable)), momentum, track.tpcitsNSigmaKa());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpctofKaon, HistTable)), momentum, track.tpctofNSigmaKa());
    }

    if (mPlotProtonPid) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kItsProton, HistTable)), momentum, track.itsNSigmaPr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcProton, HistTable)), momentum, track.tpcNSigmaPr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofProton, HistTable)), momentum, track.tofNSigmaPr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcitsProton, HistTable)), momentum, track.tpcitsNSigmaPr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpctofProton, HistTable)), momentum, track.tpctofNSigmaPr());
    }

    if (mPlotDeuteronPid) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kItsDeuteron, HistTable)), momentum, track.itsNSigmaDe());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcDeuteron, HistTable)), momentum, track.tpcNSigmaDe());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofDeuteron, HistTable)), momentum, track.tofNSigmaDe());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcitsDeuteron, HistTable)), momentum, track.tpcitsNSigmaDe());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpctofDeuteron, HistTable)), momentum, track.tpctofNSigmaDe());
    }

    if (mPlotTritonPid) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kItsTriton, HistTable)), momentum, track.itsNSigmaTr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcTriton, HistTable)), momentum, track.tpcNSigmaTr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofTriton, HistTable)), momentum, track.tofNSigmaTr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcitsTriton, HistTable)), momentum, track.tpcitsNSigmaTr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpctofTriton, HistTable)), momentum, track.tpctofNSigmaTr());
    }

    if (mPlotHeliumPid) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kItsHelium, HistTable)), momentum, track.itsNSigmaHe());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcHelium, HistTable)), momentum, track.tpcNSigmaHe());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofHelium, HistTable)), momentum, track.tofNSigmaHe());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcitsHelium, HistTable)), momentum, track.tpcitsNSigmaHe());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpctofHelium, HistTable)), momentum, track.tpctofNSigmaHe());
    }
  }

  template <modes::Mode mode, typename T1, typename T2, typename T3, typename T4>
  void fillMc(T1 const& track, T2 const& /*mcParticles*/, T3 const& /*mcMothers*/, T4 const& /*mcPartonicMothers*/)
  {
    // No MC Particle
    if (!track.has_fMcParticle()) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), 0);
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), static_cast<float>(modes::McOrigin::kNoMcParticle));
      if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
        if (mPlotOrigins) {
          mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kNoMcParticle, HistTable)), track.pt(), track.dcaXY(), track.dcaZ());
        }
      }
      return;
    }

    // Retrieve MC particle
    auto mcParticle = track.template fMcParticle_as<T2>();

    // missidentifed particles are special case
    // whether a particle is missidentfied or not cannot be known by the producer so we check it here
    bool isMissidentified = mcParticle.pdgCode() != mPdgCode;

    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTruePtVsPt, HistTable)), mcParticle.pt(), track.pt());
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTrueEtaVsEta, HistTable)), mcParticle.eta(), track.eta());
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kTruePhiVsPhi, HistTable)), mcParticle.phi(), track.phi());
    if (isMissidentified) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), static_cast<int>(modes::McOrigin::kMissidentified));
    } else {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kOrigin, HistTable)), mcParticle.origin());
    }
    mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPdg, HistTable)), mcParticle.pdgCode());

    // get mother
    if (track.has_fMcMother()) {
      auto mother = track.template fMcMother_as<T3>();
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPdgMother, HistTable)), mother.pdgCode());
    } else {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPdgMother, HistTable)), 0);
    }

    // get partonic mother
    if (track.has_fMcPartMoth()) {
      auto partonicMother = track.template fMcPartMoth_as<T4>();
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPdgPartonicMother, HistTable)), partonicMother.pdgCode());
    } else {
      mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPdgPartonicMother, HistTable)), 0);
    }

    if constexpr (modes::isFlagSet(mode, modes::Mode::kQa)) {
      if (mPlotOrigins) {
        // check first if particle is missidentified
        if (isMissidentified) {
          // if it is, we fill it as such
          mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kMissidentified, HistTable)), track.pt(), track.dcaXY(), track.dcaZ());
        } else {
          // if not, we fill it acccoridng to its origin
          switch (static_cast<modes::McOrigin>(mcParticle.origin())) {
            case modes::McOrigin::kPhysicalPrimary:
              mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kPrimary, HistTable)), track.pt(), track.dcaXY(), track.dcaZ());
              break;
            case modes::McOrigin::kFromWrongCollision:
              mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kFromWrongCollision, HistTable)), track.pt(), track.dcaXY(), track.dcaZ());
              break;
            case modes::McOrigin::kFromMaterial:
              mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kFromMaterial, HistTable)), track.pt(), track.dcaXY(), track.dcaZ());
              break;
            case modes::McOrigin::kFromSecondaryDecay:
              if (track.has_fMcMother()) {
                auto mother = track.template fMcMother_as<T3>();
                int motherPdgCode = std::abs(mother.pdgCode());
                // Switch on PDG of the mother
                if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel1 && motherPdgCode == mPdgCodesSecondaryMother[0]) {
                  mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kSecondary1, HistTable)), track.pt(), track.dcaXY(), track.dcaZ());
                } else if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel2 && motherPdgCode == mPdgCodesSecondaryMother[1]) {
                  mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kSecondary2, HistTable)), track.pt(), track.dcaXY(), track.dcaZ());
                } else if (mPlotNSecondaries >= histmanager::kSecondaryPlotLevel3 && motherPdgCode == mPdgCodesSecondaryMother[2]) {
                  mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kSecondary3, HistTable)), track.pt(), track.dcaXY(), track.dcaZ());
                } else {
                  mHistogramRegistry->fill(HIST(prefix) + HIST(McDir) + HIST(getHistName(kSecondaryOther, HistTable)), track.pt(), track.dcaXY(), track.dcaZ());
                }
              }
              break;
            default:
              LOG(warn) << "Encounted partilce with unknown origin!";
              break;
          }
        }
      }
    }
  }

  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  int mAbsCharge = 1;
  int mPdgCode = 0;
  bool mPlot2d = false;
  bool mPlotElectronPid = false;
  bool mPlotPionPid = false;
  bool mPlotKaonPid = false;
  bool mPlotProtonPid = false;
  bool mPlotDeuteronPid = false;
  bool mPlotTritonPid = false;
  bool mPlotHeliumPid = false;
  bool mPlotOrigins = false;
  int mPlotNSecondaries = 0;
  std::array<int, MaxSecondary> mPdgCodesSecondaryMother = {0};
  modes::MomentumType mMomentumType = modes::MomentumType::kPAtPv;
};
}; // namespace trackhistmanager
// aespace trackhistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_TRACKHISTMANAGER_H_
