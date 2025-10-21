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
  kItsCluster,
  kItsClusterIb,
  kTpcCrossedRows,
  kTpcCluster,
  kTpcClusterOverCrossedRows,
  kTpcClusterShared,
  kTpcClusterFractionShared,
  // kDcaxy,
  // kDcaz,
  // kDca,
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
  kTrackHistLast
};

template <const char* Prefix>
struct ConfTrackBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::ConfigurableAxis pt{"pt", {{600, 0, 6}}, "Pt"};
  o2::framework::ConfigurableAxis eta{"eta", {{300, -1.5, 1.5}}, "Eta"};
  o2::framework::ConfigurableAxis phi{"phi", {{720, 0, 1.f * o2::constants::math::TwoPI}}, "Phi"};
  o2::framework::ConfigurableAxis sign{"sign", {{3, -1.5, 1.5}}, "Sign"};
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
using ConfKinkChaDauBinning = ConfTrackBinning<PrefixKinkChaDauBinning>;
using ConfCascadePosDauBinning = ConfTrackBinning<PrefixCascadePosDauBinning>;
using ConfCascadeNegDauBinning = ConfTrackBinning<PrefixCascadeNegDauBinning>;
using ConfCascadeBachelorBinning = ConfTrackBinning<PrefixCascadeBachelorBinning>;
using ConfKinkChaDauBinning = ConfTrackBinning<PrefixKinkChaDauBinning>;

template <const char* Prefix>
struct ConfTrackQaBinning : o2::framework::ConfigurableGroup {
  std::string prefix = Prefix;
  o2::framework::Configurable<int> momentumType{"momentumType", 0, "Momentum on x-axis (0->Pt, 1->P at PV"};
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
constexpr const char PrefixTrackQaBinning2[] = "TrackQaBinning1";
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
constexpr std::array<histmanager::HistInfo<TrackHist>, kTrackHistLast> HistTable = {
  {{kPt, o2::framework::kTH1F, "hPt", "Transverse Momentum; p_{T} (GeV/#it{c}); Entries"},
   {kEta, o2::framework::kTH1F, "hEta", "Pseudorapdity; #eta; Entries"},
   {kPhi, o2::framework::kTH1F, "hPhi", "Azimuthal angle; #varphi; Entries"},
   {kSign, o2::framework::kTH1F, "hSign", "Sign of charge ; Sign; Entries"},
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
   {kItsElectron, o2::framework::kTH2F, "hItsPidElectron", "TPC PID Electron; p (GeV/#it{c}) ; n#sigma_{TPC,el}"},
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
   {kTpctofHelium, o2::framework::kTH2F, "hTpctofPidHelium", "TPC+TOF PID Helium; p (GeV/#it{c}) ; #sqrt{n#sigma_{TPC,he}^{2}+n#sigma_{TOF,he}^{2}}"}}};

template <typename T>
auto makeTrackHistSpecMap(const T& confBinningAnalysis)
{
  return std::map<TrackHist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kSign, {confBinningAnalysis.sign}}};
};

template <typename T1, typename T2>
auto makeTrackQaHistSpecMap(const T1& confBinningAnalysis, const T2 confiBinningQa)
{
  return std::map<TrackHist, std::vector<framework::AxisSpec>>{
    {kPt, {confBinningAnalysis.pt}},
    {kEta, {confBinningAnalysis.eta}},
    {kPhi, {confBinningAnalysis.phi}},
    {kSign, {confBinningAnalysis.sign}},
    {kItsCluster, {confiBinningQa.itsCluster}},
    {kItsClusterIb, {confiBinningQa.itsClusterIb}},
    {kPtVsEta, {confBinningAnalysis.pt, confBinningAnalysis.eta}},
    {kPtVsPhi, {confBinningAnalysis.pt, confBinningAnalysis.phi}},
    {kPhiVsEta, {confBinningAnalysis.phi, confBinningAnalysis.eta}},
    {kPtVsItsCluster, {confBinningAnalysis.pt, confiBinningQa.itsCluster}},
    {kPtVsTpcCluster, {confBinningAnalysis.pt, confiBinningQa.tpcCluster}},
    {kPtVsTpcCrossedRows, {confBinningAnalysis.pt, confiBinningQa.tpcCrossedRows}},
    {kPtVsTpcClusterOverCrossedRows, {confBinningAnalysis.pt, confiBinningQa.tpcClusterOverCrossedRows}},
    {kPtVsTpcClusterShared, {confBinningAnalysis.pt, confiBinningQa.tpcClusterShared}},
    {kPtVsTpcClusterFractionShared, {confBinningAnalysis.pt, confiBinningQa.tpcClusterFractionShared}},
    {kTpcClusterVsTpcCrossedRows, {confiBinningQa.tpcCluster, confiBinningQa.tpcCrossedRows}},
    {kTpcClusterVsTpcClusterShared, {confiBinningQa.tpcCluster, confiBinningQa.tpcClusterShared}},
    {kTpcCrossedRows, {confiBinningQa.tpcCrossedRows}},
    {kTpcCluster, {confiBinningQa.tpcCluster}},
    {kTpcClusterOverCrossedRows, {confiBinningQa.tpcClusterOverCrossedRows}},
    {kTpcClusterShared, {confiBinningQa.tpcClusterShared}},
    {kTpcClusterFractionShared, {confiBinningQa.tpcClusterFractionShared}},
    {kPtVsDcaxy, {confBinningAnalysis.pt, confiBinningQa.dcaXy}},
    {kPtVsDcaz, {confBinningAnalysis.pt, confiBinningQa.dcaZ}},
    {kPtVsDca, {confBinningAnalysis.pt, confiBinningQa.dca}},
    {kItsSignal, {confiBinningQa.p, confiBinningQa.itsSignal}},
    {kItsElectron, {confiBinningQa.p, confiBinningQa.itsElectron}},
    {kItsPion, {confiBinningQa.p, confiBinningQa.itsPion}},
    {kItsKaon, {confiBinningQa.p, confiBinningQa.itsKaon}},
    {kItsProton, {confiBinningQa.p, confiBinningQa.itsProton}},
    {kItsDeuteron, {confiBinningQa.p, confiBinningQa.itsDeuteron}},
    {kItsTriton, {confiBinningQa.p, confiBinningQa.itsTriton}},
    {kItsHelium, {confiBinningQa.p, confiBinningQa.itsHelium}},
    {kTpcSignal, {confiBinningQa.p, confiBinningQa.tpcSignal}},
    {kTpcElectron, {confiBinningQa.p, confiBinningQa.tpcElectron}},
    {kTpcPion, {confiBinningQa.p, confiBinningQa.tpcPion}},
    {kTpcKaon, {confiBinningQa.p, confiBinningQa.tpcKaon}},
    {kTpcProton, {confiBinningQa.p, confiBinningQa.tpcProton}},
    {kTpcDeuteron, {confiBinningQa.p, confiBinningQa.tpcDeuteron}},
    {kTpcTriton, {confiBinningQa.p, confiBinningQa.tpcTriton}},
    {kTpcHelium, {confiBinningQa.p, confiBinningQa.tpcHelium}},
    {kTofBeta, {confiBinningQa.p, confiBinningQa.tofBeta}},
    {kTofMass, {confiBinningQa.p, confiBinningQa.tofMass}},
    {kTofElectron, {confiBinningQa.p, confiBinningQa.tofElectron}},
    {kTofPion, {confiBinningQa.p, confiBinningQa.tofPion}},
    {kTofKaon, {confiBinningQa.p, confiBinningQa.tofKaon}},
    {kTofProton, {confiBinningQa.p, confiBinningQa.tofProton}},
    {kTofDeuteron, {confiBinningQa.p, confiBinningQa.tofDeuteron}},
    {kTofTriton, {confiBinningQa.p, confiBinningQa.tofTriton}},
    {kTofHelium, {confiBinningQa.p, confiBinningQa.tofHelium}},
    {kTpcitsElectron, {confiBinningQa.p, confiBinningQa.tpcitsElectron}},
    {kTpcitsPion, {confiBinningQa.p, confiBinningQa.tpcitsPion}},
    {kTpcitsKaon, {confiBinningQa.p, confiBinningQa.tpcitsKaon}},
    {kTpcitsProton, {confiBinningQa.p, confiBinningQa.tpcitsProton}},
    {kTpcitsDeuteron, {confiBinningQa.p, confiBinningQa.tpcitsDeuteron}},
    {kTpcitsTriton, {confiBinningQa.p, confiBinningQa.tpcitsTriton}},
    {kTpcitsHelium, {confiBinningQa.p, confiBinningQa.tpcitsHelium}},
    {kTpctofElectron, {confiBinningQa.p, confiBinningQa.tpctofElectron}},
    {kTpctofPion, {confiBinningQa.p, confiBinningQa.tpctofPion}},
    {kTpctofKaon, {confiBinningQa.p, confiBinningQa.tpctofKaon}},
    {kTpctofProton, {confiBinningQa.p, confiBinningQa.tpctofProton}},
    {kTpctofDeuteron, {confiBinningQa.p, confiBinningQa.tpctofDeuteron}},
    {kTpctofTriton, {confiBinningQa.p, confiBinningQa.tpctofTriton}},
    {kTpctofHelium, {confiBinningQa.p, confiBinningQa.tpctofHelium}}};
};

constexpr char PrefixTrackQa[] = "TrackQA/";
constexpr char PrefixTrack1[] = "Track1/";
constexpr char PrefixTrack2[] = "Track2/";
constexpr char PrefixTrack3[] = "Track3/";

constexpr char PrefixResonancePosDaughter[] = "ResonancePosDau/";
constexpr char PrefixResonanceNegDaughter[] = "ResonanceNegDau/";
constexpr char PrefixResonancePosDaughterQa[] = "ResonancePosDauQa/";
constexpr char PrefixResonanceNegDaughterQa[] = "ResonanceNegDauQa/";

constexpr char PrefixV0PosDaughter1[] = "V0PosDau1/";
constexpr char PrefixV0NegDaughter1[] = "V0NegDau1/";
constexpr char PrefixV0PosDaughter2[] = "V0PosDau2/";
constexpr char PrefixV0NegDaughter2[] = "V0NegDau2/";
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

/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
template <const char* prefix, modes::Mode mode>
class TrackHistManager
{
 public:
  /// Destructor
  virtual ~TrackHistManager() = default;

  void init(o2::framework::HistogramRegistry* registry, std::map<TrackHist, std::vector<o2::framework::AxisSpec>> Specs, int absCharge = 1, int momentumTypeForPid = 0)
  {
    mHistogramRegistry = registry;
    mAbsCharge = absCharge; // stored absolute charge of the track to scale the momentum in case of Z!=1 (case only for He3)

    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      std::string analysisDir = std::string(prefix) + std::string(AnalysisDir);

      mHistogramRegistry->add(analysisDir + getHistNameV2(kPt, HistTable), getHistDesc(kPt, HistTable), getHistType(kPt, HistTable), {Specs[kPt]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kEta, HistTable), getHistDesc(kEta, HistTable), getHistType(kEta, HistTable), {Specs[kEta]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kPhi, HistTable), getHistDesc(kPhi, HistTable), getHistType(kPhi, HistTable), {Specs[kPhi]});
      mHistogramRegistry->add(analysisDir + getHistNameV2(kSign, HistTable), getHistDesc(kSign, HistTable), getHistType(kSign, HistTable), {Specs[kSign]});
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {

      std::string qaDir = std::string(prefix) + std::string(QaDir);

      mHistogramRegistry->add(qaDir + getHistNameV2(kItsCluster, HistTable), getHistDesc(kItsCluster, HistTable), getHistType(kItsCluster, HistTable), {Specs[kItsCluster]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kItsClusterIb, HistTable), getHistDesc(kItsClusterIb, HistTable), getHistType(kItsClusterIb, HistTable), {Specs[kItsClusterIb]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kTpcCrossedRows, HistTable), getHistDesc(kTpcCrossedRows, HistTable), getHistType(kTpcCrossedRows, HistTable), {Specs[kTpcCrossedRows]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kTpcCluster, HistTable), getHistDesc(kTpcCluster, HistTable), getHistType(kTpcCluster, HistTable), {Specs[kTpcCluster]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kTpcClusterOverCrossedRows, HistTable), getHistDesc(kTpcClusterOverCrossedRows, HistTable), getHistType(kTpcClusterOverCrossedRows, HistTable), {Specs[kTpcClusterOverCrossedRows]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kTpcClusterShared, HistTable), getHistDesc(kTpcClusterShared, HistTable), getHistType(kTpcClusterShared, HistTable), {Specs[kTpcClusterShared]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kTpcClusterFractionShared, HistTable), getHistDesc(kTpcClusterFractionShared, HistTable), getHistType(kTpcClusterFractionShared, HistTable), {Specs[kTpcClusterFractionShared]});

      // qa 2d
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsEta, HistTable), getHistDesc(kPtVsEta, HistTable), getHistType(kPtVsEta, HistTable), {Specs[kPtVsEta]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsPhi, HistTable), getHistDesc(kPtVsPhi, HistTable), getHistType(kPtVsPhi, HistTable), {Specs[kPtVsPhi]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPhiVsEta, HistTable), getHistDesc(kPhiVsEta, HistTable), getHistType(kPhiVsEta, HistTable), {Specs[kPhiVsEta]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsItsCluster, HistTable), getHistDesc(kPtVsItsCluster, HistTable), getHistType(kPtVsItsCluster, HistTable), {Specs[kPtVsItsCluster]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsTpcCluster, HistTable), getHistDesc(kPtVsTpcCluster, HistTable), getHistType(kPtVsTpcCluster, HistTable), {Specs[kPtVsTpcCluster]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsTpcCrossedRows, HistTable), getHistDesc(kPtVsTpcCrossedRows, HistTable), getHistType(kPtVsTpcCrossedRows, HistTable), {Specs[kPtVsTpcCrossedRows]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsTpcClusterOverCrossedRows, HistTable), getHistDesc(kPtVsTpcClusterOverCrossedRows, HistTable), getHistType(kPtVsTpcClusterOverCrossedRows, HistTable), {Specs[kPtVsTpcClusterOverCrossedRows]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsTpcClusterShared, HistTable), getHistDesc(kPtVsTpcClusterShared, HistTable), getHistType(kPtVsTpcClusterShared, HistTable), {Specs[kPtVsTpcClusterShared]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsTpcClusterFractionShared, HistTable), getHistDesc(kPtVsTpcClusterFractionShared, HistTable), getHistType(kPtVsTpcClusterFractionShared, HistTable), {Specs[kPtVsTpcClusterFractionShared]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kTpcClusterVsTpcCrossedRows, HistTable), getHistDesc(kTpcClusterVsTpcCrossedRows, HistTable), getHistType(kTpcClusterVsTpcCrossedRows, HistTable), {Specs[kTpcClusterVsTpcCrossedRows]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kTpcClusterVsTpcClusterShared, HistTable), getHistDesc(kTpcClusterVsTpcClusterShared, HistTable), getHistType(kTpcClusterVsTpcClusterShared, HistTable), {Specs[kTpcClusterVsTpcClusterShared]});

      // dca
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsDcaxy, HistTable), getHistDesc(kPtVsDcaxy, HistTable), getHistType(kPtVsDcaxy, HistTable), {Specs[kPtVsDcaxy]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsDcaz, HistTable), getHistDesc(kPtVsDcaz, HistTable), getHistType(kPtVsDcaz, HistTable), {Specs[kPtVsDcaz]});
      mHistogramRegistry->add(qaDir + getHistNameV2(kPtVsDca, HistTable), getHistDesc(kPtVsDca, HistTable), getHistType(kPtVsDca, HistTable), {Specs[kPtVsDca]});

      std::string pidDir = std::string(prefix) + std::string(PidDir);

      mMomentumType = static_cast<modes::MomentumType>(momentumTypeForPid);

      mHistogramRegistry->add(pidDir + getHistNameV2(kItsSignal, HistTable), getHistDesc(kItsSignal, HistTable), getHistType(kItsSignal, HistTable), {Specs[kItsSignal]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kItsElectron, HistTable), getHistDesc(kItsElectron, HistTable), getHistType(kItsElectron, HistTable), {Specs[kItsElectron]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kItsPion, HistTable), getHistDesc(kItsPion, HistTable), getHistType(kItsPion, HistTable), {Specs[kItsPion]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kItsKaon, HistTable), getHistDesc(kItsKaon, HistTable), getHistType(kItsKaon, HistTable), {Specs[kItsKaon]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kItsProton, HistTable), getHistDesc(kItsProton, HistTable), getHistType(kItsProton, HistTable), {Specs[kItsProton]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kItsDeuteron, HistTable), getHistDesc(kItsDeuteron, HistTable), getHistType(kItsDeuteron, HistTable), {Specs[kItsDeuteron]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kItsTriton, HistTable), getHistDesc(kItsTriton, HistTable), getHistType(kItsTriton, HistTable), {Specs[kItsTriton]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kItsHelium, HistTable), getHistDesc(kItsHelium, HistTable), getHistType(kItsHelium, HistTable), {Specs[kItsHelium]});

      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcSignal, HistTable), getHistDesc(kTpcSignal, HistTable), getHistType(kTpcSignal, HistTable), {Specs[kTpcSignal]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcElectron, HistTable), getHistDesc(kTpcElectron, HistTable), getHistType(kTpcElectron, HistTable), {Specs[kTpcElectron]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcPion, HistTable), getHistDesc(kTpcPion, HistTable), getHistType(kTpcPion, HistTable), {Specs[kTpcPion]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcKaon, HistTable), getHistDesc(kTpcKaon, HistTable), getHistType(kTpcKaon, HistTable), {Specs[kTpcKaon]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcProton, HistTable), getHistDesc(kTpcProton, HistTable), getHistType(kTpcProton, HistTable), {Specs[kTpcProton]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcDeuteron, HistTable), getHistDesc(kTpcDeuteron, HistTable), getHistType(kTpcDeuteron, HistTable), {Specs[kTpcDeuteron]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcTriton, HistTable), getHistDesc(kTpcTriton, HistTable), getHistType(kTpcTriton, HistTable), {Specs[kTpcTriton]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcHelium, HistTable), getHistDesc(kTpcHelium, HistTable), getHistType(kTpcHelium, HistTable), {Specs[kTpcHelium]});

      mHistogramRegistry->add(pidDir + getHistNameV2(kTofBeta, HistTable), getHistDesc(kTofBeta, HistTable), getHistType(kTofBeta, HistTable), {Specs[kTofBeta]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTofMass, HistTable), getHistDesc(kTofMass, HistTable), getHistType(kTofMass, HistTable), {Specs[kTofMass]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTofElectron, HistTable), getHistDesc(kTofElectron, HistTable), getHistType(kTofElectron, HistTable), {Specs[kTofElectron]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTofPion, HistTable), getHistDesc(kTofPion, HistTable), getHistType(kTofPion, HistTable), {Specs[kTofPion]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTofKaon, HistTable), getHistDesc(kTofKaon, HistTable), getHistType(kTofKaon, HistTable), {Specs[kTofKaon]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTofProton, HistTable), getHistDesc(kTofProton, HistTable), getHistType(kTofProton, HistTable), {Specs[kTofProton]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTofDeuteron, HistTable), getHistDesc(kTofDeuteron, HistTable), getHistType(kTofDeuteron, HistTable), {Specs[kTofDeuteron]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTofTriton, HistTable), getHistDesc(kTofTriton, HistTable), getHistType(kTofTriton, HistTable), {Specs[kTofTriton]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTofHelium, HistTable), getHistDesc(kTofHelium, HistTable), getHistType(kTofHelium, HistTable), {Specs[kTofHelium]});

      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcitsElectron, HistTable), getHistDesc(kTpcitsElectron, HistTable), getHistType(kTpcitsElectron, HistTable), {Specs[kTpcitsElectron]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcitsPion, HistTable), getHistDesc(kTpcitsPion, HistTable), getHistType(kTpcitsPion, HistTable), {Specs[kTpcitsPion]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcitsKaon, HistTable), getHistDesc(kTpcitsKaon, HistTable), getHistType(kTpcitsKaon, HistTable), {Specs[kTpcitsKaon]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcitsProton, HistTable), getHistDesc(kTpcitsProton, HistTable), getHistType(kTpcitsProton, HistTable), {Specs[kTpcitsProton]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcitsDeuteron, HistTable), getHistDesc(kTpcitsDeuteron, HistTable), getHistType(kTpcitsDeuteron, HistTable), {Specs[kTpcitsDeuteron]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcitsTriton, HistTable), getHistDesc(kTpcitsTriton, HistTable), getHistType(kTpcitsTriton, HistTable), {Specs[kTpcitsTriton]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpcitsHelium, HistTable), getHistDesc(kTpcitsHelium, HistTable), getHistType(kTpcitsHelium, HistTable), {Specs[kTpcitsHelium]});

      mHistogramRegistry->add(pidDir + getHistNameV2(kTpctofElectron, HistTable), getHistDesc(kTpctofElectron, HistTable), getHistType(kTpctofElectron, HistTable), {Specs[kTpctofElectron]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpctofPion, HistTable), getHistDesc(kTpctofPion, HistTable), getHistType(kTpctofPion, HistTable), {Specs[kTpctofPion]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpctofKaon, HistTable), getHistDesc(kTpctofKaon, HistTable), getHistType(kTpctofKaon, HistTable), {Specs[kTpctofKaon]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpctofProton, HistTable), getHistDesc(kTpctofProton, HistTable), getHistType(kTpctofProton, HistTable), {Specs[kTpctofProton]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpctofDeuteron, HistTable), getHistDesc(kTpctofDeuteron, HistTable), getHistType(kTpctofDeuteron, HistTable), {Specs[kTpctofDeuteron]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpctofTriton, HistTable), getHistDesc(kTpctofTriton, HistTable), getHistType(kTpctofTriton, HistTable), {Specs[kTpctofTriton]});
      mHistogramRegistry->add(pidDir + getHistNameV2(kTpctofHelium, HistTable), getHistDesc(kTpctofHelium, HistTable), getHistType(kTpctofHelium, HistTable), {Specs[kTpctofHelium]});
    }
  }

  template <typename T1, typename T2>
  void fill(T1 const& track, T2 const& /*trackTable*/)
  {
    if constexpr (isFlagSet(mode, modes::Mode::kAnalysis)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPt, HistTable)), mAbsCharge * track.pt());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kEta, HistTable)), track.eta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kPhi, HistTable)), track.phi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(AnalysisDir) + HIST(getHistName(kSign, HistTable)), track.sign());
    }

    if constexpr (isFlagSet(mode, modes::Mode::kQa)) {
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kItsCluster, HistTable)), static_cast<float>(track.itsNCls()));
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kItsClusterIb, HistTable)), static_cast<float>(track.itsNClsInnerBarrel()));
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kTpcCrossedRows, HistTable)), static_cast<float>(track.tpcNClsCrossedRows()));
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kTpcCluster, HistTable)), static_cast<float>(track.tpcNClsFound()));
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kTpcClusterOverCrossedRows, HistTable)), static_cast<float>(track.tpcNClsFound()) / static_cast<float>(track.tpcNClsCrossedRows()));
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kTpcClusterShared, HistTable)), static_cast<float>(track.tpcNClsShared()));
      mHistogramRegistry->fill(HIST(prefix) + HIST(QaDir) + HIST(getHistName(kTpcClusterFractionShared, HistTable)), track.tpcSharedOverFound());

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

      float momentum = 0.f;
      if (mMomentumType == modes::MomentumType::kPAtPv) {
        momentum = mAbsCharge * track.p();
      } else if (mMomentumType == modes::MomentumType::kPt) {
        momentum = mAbsCharge * track.pt();
      }

      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kItsSignal, HistTable)), momentum, o2::analysis::femto::utils::itsSignal(track));
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kItsElectron, HistTable)), momentum, track.itsNSigmaEl());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kItsPion, HistTable)), momentum, track.itsNSigmaPi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kItsKaon, HistTable)), momentum, track.itsNSigmaKa());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kItsProton, HistTable)), momentum, track.itsNSigmaPr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kItsDeuteron, HistTable)), momentum, track.itsNSigmaDe());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kItsTriton, HistTable)), momentum, track.itsNSigmaTr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kItsHelium, HistTable)), momentum, track.itsNSigmaHe());

      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcSignal, HistTable)), momentum, track.tpcSignal());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcElectron, HistTable)), momentum, track.tpcNSigmaEl());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcPion, HistTable)), momentum, track.tpcNSigmaPi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcKaon, HistTable)), momentum, track.tpcNSigmaKa());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcProton, HistTable)), momentum, track.tpcNSigmaPr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcDeuteron, HistTable)), momentum, track.tpcNSigmaDe());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcTriton, HistTable)), momentum, track.tpcNSigmaTr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcHelium, HistTable)), momentum, track.tpcNSigmaHe());

      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofBeta, HistTable)), momentum, track.tofBeta());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofMass, HistTable)), momentum, track.tofMass());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofElectron, HistTable)), momentum, track.tofNSigmaEl());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofPion, HistTable)), momentum, track.tofNSigmaPi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofKaon, HistTable)), momentum, track.tofNSigmaKa());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofProton, HistTable)), momentum, track.tofNSigmaPr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofDeuteron, HistTable)), momentum, track.tofNSigmaDe());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofTriton, HistTable)), momentum, track.tofNSigmaTr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTofHelium, HistTable)), momentum, track.tofNSigmaHe());

      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcitsElectron, HistTable)), momentum, track.tpcitsNSigmaEl());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcitsPion, HistTable)), momentum, track.tpcitsNSigmaPi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcitsKaon, HistTable)), momentum, track.tpcitsNSigmaKa());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcitsProton, HistTable)), momentum, track.tpcitsNSigmaPr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcitsDeuteron, HistTable)), momentum, track.tpcitsNSigmaDe());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcitsTriton, HistTable)), momentum, track.tpcitsNSigmaTr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpcitsHelium, HistTable)), momentum, track.tpcitsNSigmaHe());

      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpctofElectron, HistTable)), momentum, track.tpctofNSigmaEl());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpctofPion, HistTable)), momentum, track.tpctofNSigmaPi());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpctofKaon, HistTable)), momentum, track.tpctofNSigmaKa());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpctofProton, HistTable)), momentum, track.tpctofNSigmaPr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpctofDeuteron, HistTable)), momentum, track.tpctofNSigmaDe());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpctofTriton, HistTable)), momentum, track.tpctofNSigmaTr());
      mHistogramRegistry->fill(HIST(prefix) + HIST(PidDir) + HIST(getHistName(kTpctofHelium, HistTable)), momentum, track.tpctofNSigmaHe());
    }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry = nullptr;
  float mAbsCharge = 1;
  modes::MomentumType mMomentumType = modes::MomentumType::kPAtPv;
};
}; // namespace trackhistmanager
// aespace trackhistmanager
}; // namespace o2::analysis::femto
#endif // PWGCF_FEMTO_CORE_TRACKHISTMANAGER_H_
