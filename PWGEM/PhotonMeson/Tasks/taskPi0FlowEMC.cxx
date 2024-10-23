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

/// \file taskPi0FlowEMC.cxx
/// \brief Analysis task for neutral pion flow with EMCal
///
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "CCDB/BasicCCDBManager.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Qvectors.h"

#include "PWGEM/Dilepton/Utils/EMTrackUtilities.h"
#include "PWGEM/PhotonMeson/Core/EMCPhotonCut.h"
#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/emcalHistoDefinitions.h"
#include "PWGEM/PhotonMeson/Utils/PairUtilities.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"
#include "PWGEM/PhotonMeson/Utils/NMHistograms.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photonmeson::photonpair;
using namespace o2::aod::pwgem::photon;
using namespace o2::aod::pwgem::dilepton::utils;

enum QvecEstimator { FT0M = 0,
                     FT0A = 1,
                     FT0C,
                     TPCPos,
                     TPCNeg,
                     TPCTot };

enum CentralityEstimator { None = 0,
                           CFT0A = 1,
                           CFT0C,
                           CFT0M,
                           NCentralityEstimators
};

struct EMfTaskPi0Flow {
  // configurable for flow
  Configurable<int> harmonic{"harmonic", 2, "harmonic number"};
  Configurable<int> qvecDetector{"qvecDetector", 0, "Detector for Q vector estimation (FT0M: 0, FT0A: 1, FT0C: 2, TPC Pos: 3, TPC Neg: 4, TPC Tot: 5)"};
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3)"};
  Configurable<bool> saveEpResoHisto{"saveEpResoHisto", false, "Flag to save event plane resolution histogram"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // configurable axis
  ConfigurableAxis thnConfigAxisInvMass{"thnConfigAxisInvMass", {400, 0.0, 0.8}, ""};
  ConfigurableAxis thnConfigAxisPt{"thnConfigAxisPt", {100, 0., 20.}, ""};
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {100, 0., 100.}, ""};
  ConfigurableAxis thnConfigAxisCosNPhi{"thnConfigAxisCosNPhi", {100, -1., 1.}, ""};
  ConfigurableAxis thnConfigAxisCosDeltaPhi{"thnConfigAxisCosDeltaPhi", {100, -1., 1.}, ""};
  ConfigurableAxis thnConfigAxisScalarProd{"thnConfigAxisScalarProd", {100, 0., 1.}, ""};

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcut_group";
    Configurable<float> cfgZvtxMax{"cfgZvtxMax", 10.f, "max. Zvtx"};
    Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "require sel8 in event cut"};
    Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND in event cut"};
    Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", false, "require No time frame border in event cut"};
    Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border in event cut"};
    Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"}; // ITS-TPC matched track contributes PV.
    Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
    Configurable<bool> cfgRequireEMCReadoutInMB{"cfgRequireEMCReadoutInMB", true, "require the EMC to be read out in an MB collision (kTVXinEMC)"};
    Configurable<bool> cfgRequireEMCHardwareTriggered{"cfgRequireEMCHardwareTriggered", false, "require the EMC to be hardware triggered (kEMC7 or kDMC7)"};
    Configurable<int> cfgOccupancyMin{"cfgOccupancyMin", -1, "min. occupancy"};
    Configurable<int> cfgOccupancyMax{"cfgOccupancyMax", 1000000000, "max. occupancy"};
    Configurable<bool> onlyKeepWeightedEvents{"onlyKeepWeightedEvents", false, "flag to keep only weighted events (for JJ MCs) and remove all MB events (with weight = 1)"};
  } eventcuts;

  EMCPhotonCut fEMCCut;
  struct : ConfigurableGroup {
    std::string prefix = "emccut_group";
    Configurable<float> minOpenAngle{"minOpenAngle", 0.0202, "apply min opening angle"};
    Configurable<float> EMC_minTime{"EMC_minTime", -25., "Minimum cluster time for EMCal time cut"};
    Configurable<float> EMC_maxTime{"EMC_maxTime", +30., "Maximum cluster time for EMCal time cut"};
    Configurable<float> EMC_minM02{"EMC_minM02", 0.1, "Minimum M02 for EMCal M02 cut"};
    Configurable<float> EMC_maxM02{"EMC_maxM02", 0.7, "Maximum M02 for EMCal M02 cut"};
    Configurable<float> EMC_minE{"EMC_minE", 0.7, "Minimum cluster energy for EMCal energy cut"};
    Configurable<int> EMC_minNCell{"EMC_minNCell", 1, "Minimum number of cells per cluster for EMCal NCell cut"};
    Configurable<std::vector<float>> EMC_TM_Eta{"EMC_TM_Eta", {0.01f, 4.07f, -2.5f}, "|eta| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<std::vector<float>> EMC_TM_Phi{"EMC_TM_Phi", {0.015f, 3.65f, -2.f}, "|phi| <= [0]+(pT+[1])^[2] for EMCal track matching"};
    Configurable<float> EMC_Eoverp{"EMC_Eoverp", 1.75, "Minimum cluster energy over track momentum for EMCal track matching"};
    Configurable<bool> EMC_UseExoticCut{"EMC_UseExoticCut", true, "FLag to use the EMCal exotic cluster cut"};
    Configurable<bool> EMC_UseTM{"EMC_UseTM", false, "flag to use EMCal track matching cut or not"};
  } emccuts;

  using CollsWithQvecs = soa::Join<aod::EMEvents, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>;
  using EMCalPhotons = soa::Join<aod::EMCEMEventIds, aod::SkimEMCClusters>;

  SliceCache cache;
  EventPlaneHelper epHelper;
  Preslice<EMCalPhotons> perCollision_emc = aod::emccluster::emeventId;
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void DefineEMEventCut()
  {
    fEMEventCut = EMPhotonEventCut("fEMEventCut", "fEMEventCut");
    fEMEventCut.SetRequireSel8(eventcuts.cfgRequireSel8);
    fEMEventCut.SetRequireFT0AND(eventcuts.cfgRequireFT0AND);
    fEMEventCut.SetZvtxRange(-eventcuts.cfgZvtxMax, +eventcuts.cfgZvtxMax);
    fEMEventCut.SetRequireNoTFB(eventcuts.cfgRequireNoTFB);
    fEMEventCut.SetRequireNoITSROFB(eventcuts.cfgRequireNoITSROFB);
    fEMEventCut.SetRequireNoSameBunchPileup(eventcuts.cfgRequireNoSameBunchPileup);
    fEMEventCut.SetRequireVertexITSTPC(eventcuts.cfgRequireVertexITSTPC);
    fEMEventCut.SetRequireGoodZvtxFT0vsPV(eventcuts.cfgRequireGoodZvtxFT0vsPV);
    fEMEventCut.SetRequireEMCReadoutInMB(eventcuts.cfgRequireEMCReadoutInMB);
    fEMEventCut.SetRequireEMCHardwareTriggered(eventcuts.cfgRequireEMCHardwareTriggered);
    fEMEventCut.SetOccupancyRange(eventcuts.cfgOccupancyMin, eventcuts.cfgOccupancyMax);
  }

  void DefineEMCCut()
  {
    fEMCCut = EMCPhotonCut("fEMCCut", "fEMCCut");
    const float a = emccuts.EMC_TM_Eta->at(0);
    const float b = emccuts.EMC_TM_Eta->at(1);
    const float c = emccuts.EMC_TM_Eta->at(2);

    const float d = emccuts.EMC_TM_Phi->at(0);
    const float e = emccuts.EMC_TM_Phi->at(1);
    const float f = emccuts.EMC_TM_Phi->at(2);
    LOGF(info, "EMCal track matching parameters : a = %f, b = %f, c = %f, d = %f, e = %f, f = %f", a, b, c, d, e, f);
    fEMCCut.SetTrackMatchingEta([a, b, c](float pT) { return a + pow(pT + b, c); });
    fEMCCut.SetTrackMatchingPhi([d, e, f](float pT) { return d + pow(pT + e, f); });
    fEMCCut.SetMinEoverP(emccuts.EMC_Eoverp);

    fEMCCut.SetMinE(emccuts.EMC_minE);
    fEMCCut.SetMinNCell(emccuts.EMC_minNCell);
    fEMCCut.SetM02Range(emccuts.EMC_minM02, emccuts.EMC_maxM02);
    fEMCCut.SetTimeRange(emccuts.EMC_minTime, emccuts.EMC_maxTime);
    fEMCCut.SetUseExoticCut(emccuts.EMC_UseExoticCut);
  }

  void init(InitContext&)
  {
    if (harmonic != 2 && harmonic != 3) {
      LOG(info) << "Harmonic was set to " << harmonic << " but can only be 2 or 3!";
    }

    DefineEMEventCut();
    DefineEMCCut();
    fEMCCut.SetUseTM(emccuts.EMC_UseTM); // disables TM
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&registry);

    const AxisSpec thnAxisInvMass{thnConfigAxisInvMass, "#it{M}_{#gamma#gamma} (GeV/#it{c}^{2})"};
    const AxisSpec thnAxisPt{thnConfigAxisPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality (%)"};
    const AxisSpec thnAxisCosNPhi{thnConfigAxisCosNPhi, Form("cos(%d#varphi)", harmonic.value)};
    const AxisSpec thnAxisCosDeltaPhi{thnConfigAxisCosDeltaPhi, Form("cos(%d(#varphi - #Psi_{sub}))", harmonic.value)};
    const AxisSpec thnAxisScalarProd{thnConfigAxisScalarProd, "SP"};

    registry.add("hSparsePi0Flow", "THn for SP", HistType::kTHnSparseF, {thnAxisInvMass, thnAxisPt, thnAxisCent, thnAxisCosNPhi, thnAxisCosDeltaPhi, thnAxisScalarProd});
    registry.add("spReso/hSpResoFT0cFT0a", "hSpResoFT0cFT0a; centrality; Q_{FT0c} #bullet Q_{FT0a}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0cTPCpos", "hSpResoFT0cTPCpos; centrality; Q_{FT0c} #bullet Q_{TPCpos}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0cTPCneg", "hSpResoFT0cTPCneg; centrality; Q_{FT0c} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0cTPCtot", "hSpResoFT0cTPCtot; centrality; Q_{FT0c} #bullet Q_{TPCtot}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0aTPCpos", "hSpResoFT0aTPCpos; centrality; Q_{FT0a} #bullet Q_{TPCpos}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0aTPCneg", "hSpResoFT0aTPCneg; centrality; Q_{FT0a} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0aTPCtot", "hSpResoFT0aTPCtot; centrality; Q_{FT0m} #bullet Q_{TPCtot}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0mTPCpos", "hSpResoFT0mTPCpos; centrality; Q_{FT0m} #bullet Q_{TPCpos}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0mTPCneg", "hSpResoFT0mTPCneg; centrality; Q_{FT0m} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoFT0mTPCtot", "hSpResoFT0mTPCtot; centrality; Q_{FT0m} #bullet Q_{TPCtot}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});
    registry.add("spReso/hSpResoTPCposTPCneg", "hSpResoTPCposTPCneg; centrality; Q_{TPCpos} #bullet Q_{TPCneg}", {HistType::kTH2F, {thnAxisCent, thnAxisScalarProd}});

    if (saveEpResoHisto) {
      registry.add("epReso/hEpResoFT0cFT0a", "hEpResoFT0cFT0a; centrality; #Delta#Psi_{sub}", {HistType::kTProfile2D, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0cTPCpos", "hEpResoFT0cTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTProfile2D, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0cTPCneg", "hEpResoFT0cTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTProfile2D, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0cTPCtot", "hEpResoFT0cTPCtot; centrality; #Delta#Psi_{sub}", {HistType::kTProfile2D, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0aTPCpos", "hEpResoFT0aTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTProfile2D, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0aTPCneg", "hEpResoFT0aTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTProfile2D, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0aTPCtot", "hEpResoFT0aTPCtot; centrality; #Delta#Psi_{sub}", {HistType::kTProfile2D, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0mTPCpos", "hEpResoFT0mTPCpos; centrality; #Delta#Psi_{sub}", {HistType::kTProfile2D, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0mTPCneg", "hEpResoFT0mTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTProfile2D, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoFT0mTPCtot", "hEpResoFT0mTPCtot; centrality; #Delta#Psi_{sub}", {HistType::kTProfile2D, {thnAxisCent, thnAxisCosNPhi}});
      registry.add("epReso/hEpResoTPCposTPCneg", "hEpResoTPCposTPCneg; centrality; #Delta#Psi_{sub}", {HistType::kTProfile2D, {thnAxisCent, thnAxisCosNPhi}});
    }

    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }; // end init

  /// Compute the delta psi in the range [0, pi/harmonic]
  /// \param psi1 is the first angle
  /// \param psi2 is the second angle
  float getDeltaPsiInRange(float psi1, float psi2)
  {
    float deltaPsi = psi1 - psi2;
    if (std::abs(deltaPsi) > constants::math::PI / harmonic) {
      if (deltaPsi > 0.)
        deltaPsi -= constants::math::TwoPI / harmonic;
      else
        deltaPsi += constants::math::TwoPI / harmonic;
    }
    return deltaPsi;
  }

  /// Fill THnSparse
  /// \param mass is the invariant mass of the candidate
  /// \param pt is the transverse momentum of the candidate
  /// \param cent is the centrality of the collision
  /// \param cosNPhi is the cosine of the n*phi angle
  /// \param cosDeltaPhi is the cosine of the n*(phi - evtPl) angle
  /// \param sp is the scalar product
  void fillThn(float& mass,
               float& pt,
               float& cent,
               float& cosNPhi,
               float& cosDeltaPhi,
               float& sp)
  {
    registry.fill(HIST("hSparsePi0Flow"), mass, pt, cent, cosNPhi, cosDeltaPhi, sp);
  }

  /// Get the centrality
  /// \param collision is the collision with the centrality information
  float getCentrality(CollsWithQvecs::iterator const& collision)
  {
    float cent = -999.;
    switch (centEstimator) {
      case CentralityEstimator::CFT0M:
        cent = collision.centFT0M();
        break;
      case CentralityEstimator::CFT0A:
        cent = collision.centFT0A();
        break;
      case CentralityEstimator::CFT0C:
        cent = collision.centFT0C();
        break;
      default:
        LOG(warning) << "Centrality estimator not valid. Possible values are T0M, T0A, T0C. Fallback to T0C";
        cent = collision.centFT0C();
        break;
    }
    return cent;
  }

  /// Get the Q vector
  /// \param collision is the collision with the Q vector information
  std::vector<float> getQvec(CollsWithQvecs::iterator const& collision)
  {
    float xQVec = -999.;
    float yQVec = -999.;
    switch (qvecDetector) {
      case QvecEstimator::FT0M:
        if (harmonic == 2) {
          xQVec = collision.q2xft0m();
          yQVec = collision.q2yft0m();
        } else if (harmonic == 3) {
          xQVec = collision.q3xft0m();
          yQVec = collision.q3yft0m();
        }
        break;
      case QvecEstimator::FT0A:
        if (harmonic == 2) {
          xQVec = collision.q2xft0a();
          yQVec = collision.q2yft0a();
        } else if (harmonic == 3) {
          xQVec = collision.q3xft0a();
          yQVec = collision.q3yft0a();
        }
        break;
      case QvecEstimator::FT0C:
        if (harmonic == 2) {
          xQVec = collision.q2xft0c();
          yQVec = collision.q2yft0c();
        } else if (harmonic == 3) {
          xQVec = collision.q3xft0c();
          yQVec = collision.q3yft0c();
        }
        break;
      case QvecEstimator::TPCPos:
        if (harmonic == 2) {
          xQVec = collision.q2xbpos();
          yQVec = collision.q2ybpos();
        } else if (harmonic == 3) {
          xQVec = collision.q3xbpos();
          yQVec = collision.q3ybpos();
        }
        break;
      case QvecEstimator::TPCNeg:
        if (harmonic == 2) {
          xQVec = collision.q2xbneg();
          yQVec = collision.q2ybneg();
        } else if (harmonic == 3) {
          xQVec = collision.q3xbneg();
          yQVec = collision.q3ybneg();
        }
        break;
      case QvecEstimator::TPCTot:
        if (harmonic == 2) {
          xQVec = collision.q2xbtot();
          yQVec = collision.q2ybtot();
        } else if (harmonic == 3) {
          xQVec = collision.q3xbtot();
          yQVec = collision.q3ybtot();
        }
        break;
      default:
        LOG(warning) << "Q vector estimator not valid. Please choose between FT0M, FT0A, FT0C, TPC Pos, TPC Neg. Fallback to FT0M";
        if (harmonic == 2) {
          xQVec = collision.q2xft0m();
          yQVec = collision.q2yft0m();
        } else if (harmonic == 3) {
          xQVec = collision.q3xft0m();
          yQVec = collision.q3yft0m();
        }
        break;
    }
    return {xQVec, yQVec};
  }

  /// Compute the scalar product
  /// \param collision is the collision with the Q vector information and event plane
  /// \param meson are the selected candidates
  template <typename T1>
  void runFlowAnalysis(CollsWithQvecs::iterator const& collision, T1 const& meson)
  {
    std::vector<float> qVecs = getQvec(collision);
    float xQVec = qVecs[0];
    float yQVec = qVecs[1];
    float evtPl = epHelper.GetEventPlane(xQVec, yQVec, harmonic);
    float cent = getCentrality(collision);

    float massCand = 0.;
    float ptCand = meson.pt();
    float phiCand = meson.phi();

    float cosNPhi = std::cos(harmonic * phiCand);
    float sinNPhi = std::sin(harmonic * phiCand);
    float scalprodCand = cosNPhi * xQVec + sinNPhi * yQVec;
    float cosDeltaPhi = std::cos(harmonic * (phiCand - evtPl));

    fillThn(massCand, ptCand, cent, cosNPhi, cosDeltaPhi, scalprodCand);
  }

  // Ds with rectangular cuts
  void processEMCal(CollsWithQvecs const& collisions, EMCalPhotons const& clusters)
  {
    for (auto& collision : collisions) {
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&registry, collision);
      if (!(fEMEventCut.IsSelected(collision))) {
        // no selection on the centrality is applied on purpose to allow for the resolution study in post-processing
        return;
      }
      o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&registry, collision);
      auto photons_per_collision = clusters.sliceBy(perCollision_emc, collision.globalIndex());
      for (auto& [g1, g2] : combinations(CombinationsStrictlyUpperIndexPolicy(photons_per_collision, photons_per_collision))) {
        if (!(fEMCCut.IsSelected<EMCalPhotons::iterator>(g1)) || !(fEMCCut.IsSelected<EMCalPhotons::iterator>(g2))) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v1(g1.pt(), g1.eta(), g1.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector v2(g2.pt(), g2.eta(), g2.phi(), 0.);
        ROOT::Math::PtEtaPhiMVector vPhoton = v1 + v2;
        runFlowAnalysis(collision, vPhoton);
      }
    }
  }
  PROCESS_SWITCH(EMfTaskPi0Flow, processEMCal, "Process EMCal Pi0 candidates", false);

  // Resolution
  void processResolution(CollsWithQvecs::iterator const& collision)
  {
    // we don't need to require EMCal readout for the resolution
    fEMEventCut.SetRequireEMCReadoutInMB(false);
    if (!(fEMEventCut.IsSelected(collision))) {
      // no selection on the centrality is applied on purpose to allow for the resolution study in post-processing
      return;
    }

    float centrality = getCentrality(collision); // centrality not updated in the rejection mask function
    float xQVecFT0a = -999.f;
    float yQVecFT0a = -999.f;
    float xQVecFT0c = -999.f;
    float yQVecFT0c = -999.f;
    float xQVecFT0m = -999.f;
    float yQVecFT0m = -999.f;
    float xQVecBPos = -999.f;
    float yQVecBPos = -999.f;
    float xQVecBNeg = -999.f;
    float yQVecBNeg = -999.f;
    float xQVecBTot = -999.f;
    float yQVecBTot = -999.f;
    if (harmonic == 2) {
      xQVecFT0a = collision.q2xft0a();
      yQVecFT0a = collision.q2yft0a();
      xQVecFT0c = collision.q2xft0c();
      yQVecFT0c = collision.q2yft0c();
      xQVecFT0m = collision.q2xft0m();
      yQVecFT0m = collision.q2yft0m();
      xQVecBPos = collision.q2xbpos();
      yQVecBPos = collision.q2ybpos();
      xQVecBNeg = collision.q2xbneg();
      yQVecBNeg = collision.q2ybneg();
      xQVecBTot = collision.q2xbtot();
      yQVecBTot = collision.q2ybtot();
    } else if (harmonic == 3) {
      xQVecFT0a = collision.q3xft0a();
      yQVecFT0a = collision.q3yft0a();
      xQVecFT0c = collision.q3xft0c();
      yQVecFT0c = collision.q3yft0c();
      xQVecFT0m = collision.q3xft0m();
      yQVecFT0m = collision.q3yft0m();
      xQVecBPos = collision.q3xbpos();
      yQVecBPos = collision.q3ybpos();
      xQVecBNeg = collision.q3xbneg();
      yQVecBNeg = collision.q3ybneg();
      xQVecBTot = collision.q3xbtot();
      yQVecBTot = collision.q3ybtot();
    }
    registry.fill(HIST("spReso/hSpResoFT0cFT0a"), centrality, xQVecFT0c * xQVecFT0a + yQVecFT0c * yQVecFT0a);
    registry.fill(HIST("spReso/hSpResoFT0cTPCpos"), centrality, xQVecFT0c * xQVecBPos + yQVecFT0c * yQVecBPos);
    registry.fill(HIST("spReso/hSpResoFT0cTPCneg"), centrality, xQVecFT0c * xQVecBNeg + yQVecFT0c * yQVecBNeg);
    registry.fill(HIST("spReso/hSpResoFT0cTPCtot"), centrality, xQVecFT0c * xQVecBTot + yQVecFT0c * yQVecBTot);
    registry.fill(HIST("spReso/hSpResoFT0aTPCpos"), centrality, xQVecFT0a * xQVecBPos + yQVecFT0a * yQVecBPos);
    registry.fill(HIST("spReso/hSpResoFT0aTPCneg"), centrality, xQVecFT0a * xQVecBNeg + yQVecFT0a * yQVecBNeg);
    registry.fill(HIST("spReso/hSpResoFT0aTPCtot"), centrality, xQVecFT0a * xQVecBTot + yQVecFT0a * yQVecBTot);
    registry.fill(HIST("spReso/hSpResoFT0mTPCpos"), centrality, xQVecFT0m * xQVecBPos + yQVecFT0m * yQVecBPos);
    registry.fill(HIST("spReso/hSpResoFT0mTPCneg"), centrality, xQVecFT0m * xQVecBNeg + yQVecFT0m * yQVecBNeg);
    registry.fill(HIST("spReso/hSpResoFT0mTPCtot"), centrality, xQVecFT0m * xQVecBTot + yQVecFT0m * yQVecBTot);
    registry.fill(HIST("spReso/hSpResoTPCposTPCneg"), centrality, xQVecBPos * xQVecBNeg + yQVecBPos * yQVecBNeg);

    if (saveEpResoHisto) {
      float epFT0a = epHelper.GetEventPlane(xQVecFT0a, yQVecFT0a, harmonic);
      float epFT0c = epHelper.GetEventPlane(xQVecFT0c, yQVecFT0c, harmonic);
      float epFT0m = epHelper.GetEventPlane(xQVecFT0m, yQVecFT0m, harmonic);
      float epBPoss = epHelper.GetEventPlane(xQVecBPos, yQVecBPos, harmonic);
      float epBNegs = epHelper.GetEventPlane(xQVecBNeg, yQVecBNeg, harmonic);
      float epBTots = epHelper.GetEventPlane(xQVecBTot, yQVecBTot, harmonic);

      registry.fill(HIST("epReso/hEpResoFT0cFT0a"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epFT0a)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFT0cTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0c, epBTots)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFT0aTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0a, epBTots)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFT0mTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFT0m, epBTots)));
      registry.fill(HIST("epReso/hEpResoTPCposTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epBPoss, epBNegs)));
    }
  }
  PROCESS_SWITCH(EMfTaskPi0Flow, processResolution, "Process resolution", false);

}; // End struct EMfTaskPi0Flow

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<EMfTaskPi0Flow>(cfgc)};
}
