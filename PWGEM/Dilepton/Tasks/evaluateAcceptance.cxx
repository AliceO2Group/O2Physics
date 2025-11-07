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

/// \file evaluateAcceptance.cxx
/// \brief a task to evaluate pair acceptance in MC
/// \author daiki.sekihata@cern.ch

#include "PWGEM/Dilepton/Utils/MCUtilities.h"
#include "PWGEM/Dilepton/Utils/PairUtilities.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TableHelper.h"

// #include "Common/Core/trackUtilities.h"
// #include "Common/DataModel/Centrality.h"
// #include "Common/DataModel/CollisionAssociationTables.h"
// #include "Common/DataModel/EventSelection.h"
// #include "Common/DataModel/Multiplicity.h"
// #include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/LHCConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DataFormatsParameters/GRPLHCIFData.h"

// #include "DataFormatsCalibration/MeanVertexObject.h"
// #include "DataFormatsParameters/GRPMagField.h"
// #include "DataFormatsParameters/GRPObject.h"
// #include "DetectorsBase/GeometryManager.h"
// #include "DetectorsBase/Propagator.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"

#include <algorithm>
#include <array>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::soa;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::aod::pwgem::dilepton::utils::mcutil;
using namespace o2::aod::pwgem::dilepton::utils::pairutil;

struct evaluateAcceptance {
  // Configurables
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  Configurable<int> cfgEventGeneratorType{"cfgEventGeneratorType", -1, "if positive, select event generator type. i.e. gap or signal"};
  Configurable<int> cfgPdgLepton{"cfgPdgLepton", 11, "pdg code 11 or 13"};
  ConfigurableAxis ConfMllBins{"ConfMllBins", {400, 0, 4}, "mll bins"};
  ConfigurableAxis ConfPtllBins{"ConfPtllBins", {100, 0, 10}, "pTll bins"};
  ConfigurableAxis ConfYllBins{"ConfYllBins", {400, -10, +10}, "yll bins"};
  ConfigurableAxis ConfCosThetaBins{"ConfCosThetaBins", {40, -1, +1}, "cos theta bins for polarization"};
  ConfigurableAxis ConfPhiBins{"ConfPhiBins", {72, -M_PI, M_PI}, "phi bins for polarization"};
  ConfigurableAxis ConfQuadMomBins{"ConfQuadMomBins", {150, -0.5, 1}, "quadrupole moment bins for polarization"};
  ConfigurableAxis ConfPtlBins{"ConfPtlBins", {200, 0, 10}, "pTl bins"};
  ConfigurableAxis ConfEtalBins{"ConfEtalBins", {200, -10, 10}, "etal bins"};

  HistogramRegistry fRegistry{"fRegistry"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = 0;

  float beamM1 = o2::constants::physics::MassProton; // mass of beam
  float beamM2 = o2::constants::physics::MassProton; // mass of beam
  float beamE1 = 0.f;                                // beam energy
  float beamE2 = 0.f;                                // beam energy
  float beamP1 = 0.f;                                // beam momentum
  float beamP2 = 0.f;                                // beam momentum

  float leptonM1 = 0.f;
  float leptonM2 = 0.f;
  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    if (cfgPdgLepton.value == 11) {
      leptonM1 = o2::constants::physics::MassElectron;
      leptonM2 = o2::constants::physics::MassElectron;
    } else if (cfgPdgLepton.value == 13) {
      leptonM1 = o2::constants::physics::MassMuon;
      leptonM2 = o2::constants::physics::MassMuon;
    } else {
      LOGF(fatal, "pdg code must be 11 or 13.");
    }

    addHistograms();
  }

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", bc.timestamp());
    int beamZ1 = grplhcif->getBeamZ(o2::constants::lhc::BeamC);
    int beamZ2 = grplhcif->getBeamZ(o2::constants::lhc::BeamA);
    int beamA1 = grplhcif->getBeamA(o2::constants::lhc::BeamC);
    int beamA2 = grplhcif->getBeamA(o2::constants::lhc::BeamA);
    beamE1 = grplhcif->getBeamEnergyPerNucleonInGeV(o2::constants::lhc::BeamC);
    beamE2 = grplhcif->getBeamEnergyPerNucleonInGeV(o2::constants::lhc::BeamA);
    beamM1 = o2::constants::physics::MassProton * beamA1;
    beamM2 = o2::constants::physics::MassProton * beamA2;
    beamP1 = std::sqrt(std::pow(beamE1, 2) - std::pow(beamM1, 2));
    beamP2 = std::sqrt(std::pow(beamE2, 2) - std::pow(beamM2, 2));
    LOGF(info, "beamZ1 = %d, beamZ2 = %d, beamA1 = %d, beamA2 = %d, beamE1 = %f (GeV), beamE2 = %f (GeV), beamM1 = %f (GeV), beamM2 = %f (GeV), beamP1 = %f (GeV), beamP2 = %f (GeV)", beamZ1, beamZ2, beamA1, beamA2, beamE1, beamE2, beamM1, beamM2, beamP1, beamP2);
    mRunNumber = bc.runNumber();
  }

  static constexpr std::string_view pair_sign_types[3] = {"uls/", "lspp/", "lsmm/"};
  static constexpr std::string_view dilepton_source_types[20] = {
    "sm/Pi0/",                // 0
    "sm/Eta/",                // 1
    "sm/EtaPrime/",           // 2
    "sm/Rho/",                // 3
    "sm/Omega/",              // 4
    "sm/Omega2ll/",           // 5
    "sm/Phi/",                // 6
    "sm/Phi2ll/",             // 7
    "sm/PromptJPsi/",         // 8
    "sm/NonPromptJPsi/",      // 9
    "sm/PromptPsi2S/",        // 10
    "sm/NonPromptPsi2S/",     // 11
    "sm/Upsilon1S/",          // 12
    "sm/Upsilon2S/",          // 13
    "sm/Upsilon3S/",          // 14
    "ccbar/c2l_c2l/",         // 15
    "bbbar/b2l_b2l/",         // 16
    "bbbar/b2c2l_b2c2l/",     // 17
    "bbbar/b2c2l_b2l_sameb/", // 18
    "bbbar/b2c2l_b2l_diffb/"  // 19
  }; // unordered_map is better, but cannot be constexpr.

  void addHistograms()
  {
    auto hCollisionCounter = fRegistry.add<TH1>("Event/hCollisionCounter", "collision counter", kTH1D, {{2, -0.5f, 1.5f}}, false);
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "accepted");

    const AxisSpec axisMll{ConfMllBins, "m_{ll} (GeV/c^{2})"};
    const AxisSpec axisPtll{ConfPtllBins, "p_{T,ll} (GeV/c)"};
    const AxisSpec axisYll{ConfYllBins, "y_{ll}"};
    const AxisSpec axisCosThetaCS{ConfCosThetaBins, "cos(#theta^{CS})"};
    const AxisSpec axisPhiCS{ConfPhiBins, "#varphi^{CS} (rad.)"};
    const AxisSpec axisQuadMomCS{ConfQuadMomBins, "#frac{3 cos^{2}(#theta^{CS}) #minus 1}{2}"};
    const AxisSpec axisCosThetaHX{ConfCosThetaBins, "cos(#theta^{HX})"};
    const AxisSpec axisPhiHX{ConfPhiBins, "#varphi^{HX} (rad.)"};
    const AxisSpec axisQuadMomHX{ConfQuadMomBins, "#frac{3 cos^{2}(#theta^{HX}) #minus 1}{2}"};

    const AxisSpec axisPtl1{ConfPtlBins, "p_{T,l1} (GeV/c)"};
    const AxisSpec axisPtl2{ConfPtlBins, "p_{T,l2} (GeV/c)"};
    const AxisSpec axisEtal1{ConfEtalBins, "#eta_{l1}"};
    const AxisSpec axisEtal2{ConfEtalBins, "#eta_{l2}"};

    // for pairs
    fRegistry.add("Generated/sm/Pi0/uls/hs", "gen. dilepton", kTHnSparseD, {axisMll, axisPtll, axisYll, axisCosThetaCS, axisPhiCS, axisQuadMomCS, axisCosThetaHX, axisPhiHX, axisQuadMomHX, axisPtl1, axisPtl2, axisEtal1, axisEtal2}, true);
    fRegistry.addClone("Generated/sm/Pi0/uls/", "Generated/sm/Pi0/lspp/");
    fRegistry.addClone("Generated/sm/Pi0/uls/", "Generated/sm/Pi0/lsmm/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/Eta/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/EtaPrime/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/Rho/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/Omega/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/Omega2ll/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/Phi/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/Phi2ll/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/PromptJPsi/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/NonPromptJPsi/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/PromptPsi2S/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/NonPromptPsi2S/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/Upsilon1S/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/Upsilon2S/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/sm/Upsilon3S/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/ccbar/c2l_c2l/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/bbbar/b2l_b2l/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/bbbar/b2c2l_b2c2l/");
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/bbbar/b2c2l_b2l_sameb/"); // ULS
    fRegistry.addClone("Generated/sm/Pi0/", "Generated/bbbar/b2c2l_b2l_diffb/"); // LS
  }

  template <int pairSignId, typename TLepton, typename TMCParticles>
  void fillGenPairInfo(TLepton const& t1, TLepton const& t2, TMCParticles const& mcParticles)
  {
    if (!t1.isPhysicalPrimary() && !t1.producedByGenerator()) {
      return;
    }
    if (!t2.isPhysicalPrimary() && !t2.producedByGenerator()) {
      return;
    }

    int mother_id = std::max({FindSMULS(t1, t2, mcParticles), FindSMULS(t2, t1, mcParticles), FindSMLSPP(t1, t2, mcParticles), FindSMLSMM(t1, t2, mcParticles)});
    int hfee_type = IsHF(t1, t2, mcParticles);
    if (mother_id < 0 && hfee_type < 0) {
      return;
    }

    ROOT::Math::PtEtaPhiMVector v1(t1.pt(), t1.eta(), t1.phi(), leptonM1);
    ROOT::Math::PtEtaPhiMVector v2(t2.pt(), t2.eta(), t2.phi(), leptonM2);
    ROOT::Math::PtEtaPhiMVector v12 = v1 + v2;

    int sign1 = t1.pdgCode() > 0 ? -1 : 1;
    std::array<float, 4> arrP1 = {t1.px(), t1.py(), t1.pz(), leptonM1};
    std::array<float, 4> arrP2 = {t2.px(), t2.py(), t2.pz(), leptonM2};
    float cosThetaCS = 999, phiCS = 999.f;
    float cosThetaHX = 999, phiHX = 999.f;
    o2::aod::pwgem::dilepton::utils::pairutil::getAngleCS(arrP1, arrP2, beamE1, beamE2, beamP1, beamP2, sign1, cosThetaCS, phiCS);
    o2::aod::pwgem::dilepton::utils::pairutil::getAngleHX(arrP1, arrP2, beamE1, beamE2, beamP1, beamP2, sign1, cosThetaHX, phiHX);
    o2::math_utils::bringToPMPi(phiCS);
    o2::math_utils::bringToPMPi(phiHX);
    float quadmomCS = (3.f * std::pow(cosThetaCS, 2) - 1.f) / 2.f;
    float quadmomHX = (3.f * std::pow(cosThetaHX, 2) - 1.f) / 2.f;

    if (mother_id > -1) {
      auto mcmother = mcParticles.iteratorAt(mother_id);
      int nd = mcmother.daughtersIds()[1] - mcmother.daughtersIds()[0] + 1; // number of daughters
      switch (std::abs(mcmother.pdgCode())) {
        case 111:
          fRegistry.fill(HIST("Generated/sm/Pi0/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          break;
        case 221:
          fRegistry.fill(HIST("Generated/sm/Eta/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          break;
        case 331:
          fRegistry.fill(HIST("Generated/sm/EtaPrime/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          break;
        case 113:
          fRegistry.fill(HIST("Generated/sm/Rho/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          break;
        case 223:
          fRegistry.fill(HIST("Generated/sm/Omega/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          if (nd == 2) {
            fRegistry.fill(HIST("Generated/sm/Omega2ll/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          }
          break;
        case 333:
          fRegistry.fill(HIST("Generated/sm/Phi/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          if (nd == 2) {
            fRegistry.fill(HIST("Generated/sm/Phi2ll/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          }
          break;
        case 443:
          if (IsFromBeauty(mcmother, mcParticles) > 0) {
            fRegistry.fill(HIST("Generated/sm/NonPromptJPsi/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          } else {
            fRegistry.fill(HIST("Generated/sm/PromptJPsi/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          }
          break;
        case 100443:
          if (IsFromBeauty(mcmother, mcParticles) > 0) {
            fRegistry.fill(HIST("Generated/sm/NonPromptPsi2S/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          } else {
            fRegistry.fill(HIST("Generated/sm/PromptPsi2S/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          }
          break;
        case 553:
          fRegistry.fill(HIST("Generated/sm/Upsilon1S/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          break;
        case 100553:
          fRegistry.fill(HIST("Generated/sm/Upsilon2S/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          break;
        case 200553:
          fRegistry.fill(HIST("Generated/sm/Upsilon3S/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          break;
        default:
          break;
      }
    } else if (hfee_type > -1) {
      switch (hfee_type) {
        case static_cast<int>(EM_HFeeType::kCe_Ce):
          fRegistry.fill(HIST("Generated/ccbar/c2l_c2l/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          break;
        case static_cast<int>(EM_HFeeType::kBe_Be):
          fRegistry.fill(HIST("Generated/bbbar/b2l_b2l/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          break;
        case static_cast<int>(EM_HFeeType::kBCe_BCe):
          fRegistry.fill(HIST("Generated/bbbar/b2c2l_b2c2l/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          break;
        case static_cast<int>(EM_HFeeType::kBCe_Be_SameB):
          fRegistry.fill(HIST("Generated/bbbar/b2c2l_b2l_sameb/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          break;
        case static_cast<int>(EM_HFeeType::kBCe_Be_DiffB):
          fRegistry.fill(HIST("Generated/bbbar/b2c2l_b2l_diffb/") + HIST(pair_sign_types[pairSignId]) + HIST("hs"), v12.M(), v12.Pt(), v12.Rapidity(), cosThetaCS, phiCS, quadmomCS, cosThetaHX, phiHX, quadmomHX, t1.pt(), t2.pt(), t1.eta(), t2.eta());
          break;
        default:
          break;
      }
    }
  }

  template <typename TTrack, typename TMCParticles>
  int FindSMULS(TTrack const& t1mc, TTrack const& t2mc, TMCParticles const& mcParticles)
  {
    int arr[] = {
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -cfgPdgLepton, cfgPdgLepton, 111, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -cfgPdgLepton, cfgPdgLepton, 221, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -cfgPdgLepton, cfgPdgLepton, 331, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -cfgPdgLepton, cfgPdgLepton, 113, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -cfgPdgLepton, cfgPdgLepton, 223, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -cfgPdgLepton, cfgPdgLepton, 333, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -cfgPdgLepton, cfgPdgLepton, 443, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -cfgPdgLepton, cfgPdgLepton, 100443, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -cfgPdgLepton, cfgPdgLepton, 553, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -cfgPdgLepton, cfgPdgLepton, 100553, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -cfgPdgLepton, cfgPdgLepton, 200553, mcParticles)};
    int size = sizeof(arr) / sizeof(*arr);
    int max = *std::max_element(arr, arr + size);
    return max;
  }

  template <typename TTrack, typename TMCParticles>
  int FindSMLSPP(TTrack const& t1mc, TTrack const& t2mc, TMCParticles const& mcParticles)
  {
    int arr[] = {
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -cfgPdgLepton, -cfgPdgLepton, 221, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -cfgPdgLepton, -cfgPdgLepton, 331, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -cfgPdgLepton, -cfgPdgLepton, 113, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -cfgPdgLepton, -cfgPdgLepton, 223, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, -cfgPdgLepton, -cfgPdgLepton, 333, mcParticles)};
    int size = sizeof(arr) / sizeof(*arr);
    int max = *std::max_element(arr, arr + size);
    return max;
  }

  template <typename TTrack, typename TMCParticles>
  int FindSMLSMM(TTrack const& t1mc, TTrack const& t2mc, TMCParticles const& mcParticles)
  {
    int arr[] = {
      FindCommonMotherFrom2Prongs(t1mc, t2mc, cfgPdgLepton, cfgPdgLepton, 221, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, cfgPdgLepton, cfgPdgLepton, 331, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, cfgPdgLepton, cfgPdgLepton, 113, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, cfgPdgLepton, cfgPdgLepton, 223, mcParticles),
      FindCommonMotherFrom2Prongs(t1mc, t2mc, cfgPdgLepton, cfgPdgLepton, 333, mcParticles)};
    int size = sizeof(arr) / sizeof(*arr);
    int max = *std::max_element(arr, arr + size);
    return max;
  }

  SliceCache cache;
  Preslice<aod::McParticles> perMCCollision = o2::aod::mcparticle::mcCollisionId;
  Partition<aod::McParticles> posLeptons = o2::aod::mcparticle::pdgCode == -cfgPdgLepton; // l+
  Partition<aod::McParticles> negLeptons = o2::aod::mcparticle::pdgCode == cfgPdgLepton;  // l-

  void process(aod::McCollisions const& mcCollisions, aod::McParticles const& mcParticles, aod::BCsWithTimestamps const&)
  {
    for (const auto& mcCollision : mcCollisions) {
      auto bc = mcCollision.template bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);

      fRegistry.fill(HIST("Event/hCollisionCounter"), 0);
      if (cfgEventGeneratorType >= 0 && mcCollision.getSubGeneratorId() != cfgEventGeneratorType) {
        continue;
      }
      fRegistry.fill(HIST("Event/hCollisionCounter"), 1);

      auto posLeptons_per_coll = posLeptons->sliceByCached(o2::aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      auto negLeptons_per_coll = negLeptons->sliceByCached(o2::aod::mcparticle::mcCollisionId, mcCollision.globalIndex(), cache);
      // LOGF(info, "mcCollision.globalIndex() = %d, posLeptons_per_coll.size() = %d, negLeptons_per_coll.size() = %d", mcCollision.globalIndex(), posLeptons_per_coll.size(), negLeptons_per_coll.size());

      for (const auto& [t1, t2] : combinations(CombinationsFullIndexPolicy(posLeptons_per_coll, negLeptons_per_coll))) { // ULS
        if (!(t1.isPhysicalPrimary() || t1.producedByGenerator()) || !(t2.isPhysicalPrimary() || t2.producedByGenerator())) {
          continue;
        }
        fillGenPairInfo<0>(t1, t2, mcParticles);
      } // end of ULS pairing

      for (const auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(posLeptons_per_coll, posLeptons_per_coll))) { // LS++
        if (!(t1.isPhysicalPrimary() || t1.producedByGenerator()) || !(t2.isPhysicalPrimary() || t2.producedByGenerator())) {
          continue;
        }
        fillGenPairInfo<1>(t1, t2, mcParticles);
      } // end of LS++ pairing

      for (const auto& [t1, t2] : combinations(CombinationsStrictlyUpperIndexPolicy(negLeptons_per_coll, negLeptons_per_coll))) { // LS--
        if (!(t1.isPhysicalPrimary() || t1.producedByGenerator()) || !(t2.isPhysicalPrimary() || t2.producedByGenerator())) {
          continue;
        }
        fillGenPairInfo<2>(t1, t2, mcParticles);
      } // end of LS++ pairing

    } // end of mc collision loop
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<evaluateAcceptance>(cfgc, TaskName{"evaluate-acceptance"})};
}
