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

/// \file taskFlowReso.cxx
/// \brief Analysis task for flow resolution
/// \author M. Hemmer, marvin.hemmer@cern.ch

#include "PWGEM/PhotonMeson/Core/EMPhotonEventCut.h"
#include "PWGEM/PhotonMeson/DataModel/gammaTables.h"
#include "PWGEM/PhotonMeson/Utils/EventHistograms.h"

#include "Common/Core/EventPlaneHelper.h"
#include "Common/Core/RecoDecay.h"

#include <CommonConstants/MathConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <TString.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <utility>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::aod::pwgem::photon;

enum QvecEstimator {
  FT0M = 0,
  FT0A = 1,
  FT0C,
  TPCPos,
  TPCNeg,
  TPCTot,
  FV0A
};

enum CentralityEstimator {
  None = 0,
  CFT0A = 1,
  CFT0C,
  CFT0M,
  NCentralityEstimators
};

enum Harmonics {
  kNone = 0,
  kDirect = 1,
  kElliptic = 2,
  kTriangluar = 3,
  kQuadrangular = 4,
  kPentagonal = 5,
  kHexagonal = 6,
  kHeptagonal = 7,
  kOctagonal = 8
};

struct TaskFlowReso {

  // configurable for flow
  Configurable<int> harmonic{"harmonic", 2, "harmonic number"};
  Configurable<int> qvecDetector{"qvecDetector", 0, "Detector for Q vector estimation (FT0M: 0, FT0A: 1, FT0C: 2, TPC Pos: 3, TPC Neg: 4, TPC Tot: 5, FV0A: 6)"};
  Configurable<int> qvecSubADetector{"qvecSubADetector", 3, "Sub A Detector for Q vector estimation for resolution (FT0M: 0, FT0A: 1, FT0C: 2, TPC Pos: 3, TPC Neg: 4, TPC Tot: 5, FV0A: 6)"};
  Configurable<int> qvecSubBDetector{"qvecSubBDetector", 4, "Sub B Detector for Q vector estimation for resolution (FT0M: 0, FT0A: 1, FT0C: 2, TPC Pos: 3, TPC Neg: 4, TPC Tot: 5, FV0A: 6)"};
  Configurable<int> centEstimator{"centEstimator", 2, "Centrality estimation (FT0A: 1, FT0C: 2, FT0M: 3)"};
  Configurable<bool> saveEpResoHisto{"saveEpResoHisto", true, "Flag to save event plane resolution histogram"};
  Configurable<bool> saveSPResoHist{"saveSPResoHist", true, "Flag to save scalar product resolution histogram"};
  Configurable<float> cfgMaxQVector{"cfgMaxQVector", 999999999.f, "Maximum allowed absolute QVector value."};

  // configurable axis
  ConfigurableAxis thnConfigAxisCent{"thnConfigAxisCent", {20, 0., 100.}, "centrality axis for the current event"};
  ConfigurableAxis thnConfigAxisCosNPhi{"thnConfigAxisCosNPhi", {100, -1., 1.}, "cos(n*phi) axis for the current event"};

  EMPhotonEventCut fEMEventCut;
  struct : ConfigurableGroup {
    std::string prefix = "eventcuts";
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
    Configurable<int> cfgTrackOccupancyMin{"cfgTrackOccupancyMin", -1, "min. track occupancy"};
    Configurable<int> cfgTrackOccupancyMax{"cfgTrackOccupancyMax", 1000000000, "max. track occupancy"};
    Configurable<float> cfgFT0COccupancyMin{"cfgFT0COccupancyMin", -1, "min. FT0C occupancy"};
    Configurable<float> cfgFT0COccupancyMax{"cfgFT0COccupancyMax", 1000000000, "max. FT0C occupancy"};
    Configurable<float> cfgMinCent{"cfgMinCent", 0, "min. centrality (%)"};
    Configurable<float> cfgMaxCent{"cfgMaxCent", 90, "max. centrality (%)"};
    Configurable<bool> onlyKeepWeightedEvents{"onlyKeepWeightedEvents", false, "flag to keep only weighted events (for JJ MCs) and remove all MB events (with weight = 1)"};
  } eventcuts;

  SliceCache cache;
  EventPlaneHelper epHelper;

  using CollsWithQvecs = soa::Join<aod::EMEvents, aod::EMEventsAlias, aod::EMEventsMult, aod::EMEventsCent, aod::EMEventsQvec>;
  using CollWithQvec = CollsWithQvecs::iterator;

  static constexpr std::size_t NQVecEntries = 6;

  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject, false, false};

  void defineEMEventCut()
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
  }

  void init(InitContext&)
  {
    if (harmonic != kElliptic && harmonic != kTriangluar) {
      LOG(info) << "Harmonic was set to " << harmonic << " but can only be 2 or 3!";
    }

    defineEMEventCut();
    o2::aod::pwgem::photonmeson::utils::eventhistogram::addEventHistograms(&registry);

    const AxisSpec thnAxisCent{thnConfigAxisCent, "Centrality (%)"};
    const AxisSpec thnAxisCosNPhi{thnConfigAxisCosNPhi, Form("cos(%d#varphi)", harmonic.value)};
    const AxisSpec thAxisPsi{360 / harmonic.value, -(1. / static_cast<float>(harmonic.value)) * std::numbers::pi_v<float>, (1. / static_cast<float>(harmonic.value)) * std::numbers::pi_v<float>, Form("#Psi_{%d}", harmonic.value)};
    const AxisSpec thAxisCN{8, 0.5, 8.5, "#it{c}_{n}"};
    const AxisSpec thAxisSN{8, 0.5, 8.5, "#it{s}_{n}"};
    const AxisSpec thAxisAzimuth{360, -o2::constants::math::PI, o2::constants::math::PI, "#it{#varphi} (rad)"};

    if (saveSPResoHist.value) {
      registry.add("spReso/hSpResoFT0cFT0a", "hSpResoFT0cFT0a; centrality; Q_{FT0c} #bullet Q_{FT0a}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0cTPCpos", "hSpResoFT0cTPCpos; centrality; Q_{FT0c} #bullet Q_{TPCpos}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0cTPCneg", "hSpResoFT0cTPCneg; centrality; Q_{FT0c} #bullet Q_{TPCneg}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0cTPCtot", "hSpResoFT0cTPCtot; centrality; Q_{FT0c} #bullet Q_{TPCtot}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0aTPCpos", "hSpResoFT0aTPCpos; centrality; Q_{FT0a} #bullet Q_{TPCpos}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0aTPCneg", "hSpResoFT0aTPCneg; centrality; Q_{FT0a} #bullet Q_{TPCneg}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0aTPCtot", "hSpResoFT0aTPCtot; centrality; Q_{FT0m} #bullet Q_{TPCtot}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0mTPCpos", "hSpResoFT0mTPCpos; centrality; Q_{FT0m} #bullet Q_{TPCpos}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0mTPCneg", "hSpResoFT0mTPCneg; centrality; Q_{FT0m} #bullet Q_{TPCneg}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFT0mTPCtot", "hSpResoFT0mTPCtot; centrality; Q_{FT0m} #bullet Q_{TPCtot}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoTPCposTPCneg", "hSpResoTPCposTPCneg; centrality; Q_{TPCpos} #bullet Q_{TPCneg}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFV0aFT0c", "hSpResoFV0aFT0c; centrality; Q_{FV0a} #bullet Q_{FT0c}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFV0aTPCpos", "hSpResoFV0aTPCpos; centrality; Q_{FV0a} #bullet Q_{TPCpos}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFV0aTPCneg", "hSpResoFV0aTPCneg; centrality; Q_{FV0a} #bullet Q_{TPCneg}", HistType::kTProfile, {thnAxisCent});
      registry.add("spReso/hSpResoFV0aTPCtot", "hSpResoFV0aTPCtot; centrality; Q_{FV0a} #bullet Q_{TPCtot}", HistType::kTProfile, {thnAxisCent});
    }

    if (saveEpResoHisto.value) {
      registry.add("hEventPlaneAngleFT0C", "hEventPlaneAngleFT0C", HistType::kTH2D, {thnAxisCent, thAxisPsi});
      registry.add("hEventPlaneAngleFT0A", "hEventPlaneAngleFT0A", HistType::kTH2D, {thnAxisCent, thAxisPsi});
      registry.add("hEventPlaneAngleFT0M", "hEventPlaneAngleFT0M", HistType::kTH2D, {thnAxisCent, thAxisPsi});
      registry.add("hEventPlaneAngleFV0A", "hEventPlaneAngleFV0A", HistType::kTH2D, {thnAxisCent, thAxisPsi});
      registry.add("hEventPlaneAngleTPCpos", "hEventPlaneAngleTPCpos", HistType::kTH2D, {thnAxisCent, thAxisPsi});
      registry.add("hEventPlaneAngleTPCneg", "hEventPlaneAngleTPCneg", HistType::kTH2D, {thnAxisCent, thAxisPsi});
      registry.add("epReso/hEpResoFT0cFT0a", "hEpResoFT0cFT0a; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0cTPCpos", "hEpResoFT0cTPCpos; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0cTPCneg", "hEpResoFT0cTPCneg; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0cTPCtot", "hEpResoFT0cTPCtot; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0aTPCpos", "hEpResoFT0aTPCpos; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0aTPCneg", "hEpResoFT0aTPCneg; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0aTPCtot", "hEpResoFT0aTPCtot; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0mTPCpos", "hEpResoFT0mTPCpos; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0mTPCneg", "hEpResoFT0mTPCneg; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFT0mTPCtot", "hEpResoFT0mTPCtot; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoTPCposTPCneg", "hEpResoTPCposTPCneg; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFV0aFT0c", "hEpResoFV0aFT0c; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFV0aTPCpos", "hEpResoFV0aTPCpos; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFV0aTPCneg", "hEpResoFV0aTPCneg; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpResoFV0aTPCtot", "hEpResoFV0aTPCtot; centrality; #Delta#Psi_{sub}", HistType::kTH2D, {thnAxisCent, thnAxisCosNPhi});
      registry.add("epReso/hEpCosCoefficientsFT0c", "hEpCosCoefficientsFT0c; centrality; c_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisCN});
      registry.add("epReso/hEpSinCoefficientsFT0c", "hEpSinCoefficientsFT0c; centrality; s_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisSN});
      registry.add("epReso/hEpCosCoefficientsFT0a", "hEpCosCoefficientsFT0a; centrality; c_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisCN});
      registry.add("epReso/hEpSinCoefficientsFT0a", "hEpSinCoefficientsFT0a; centrality; s_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisSN});
      registry.add("epReso/hEpCosCoefficientsFV0a", "hEpCosCoefficientsFV0a; centrality; c_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisCN});
      registry.add("epReso/hEpSinCoefficientsFV0a", "hEpSinCoefficientsFV0a; centrality; s_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisSN});
      registry.add("epReso/hEpCosCoefficientsFT0m", "hEpCosCoefficientsFT0m; centrality; c_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisCN});
      registry.add("epReso/hEpSinCoefficientsFT0m", "hEpSinCoefficientsFT0m; centrality; s_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisSN});
      registry.add("epReso/hEpCosCoefficientsTPCpos", "hEpCosCoefficientsTPCpos; centrality; c_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisCN});
      registry.add("epReso/hEpSinCoefficientsTPCpos", "hEpSinCoefficientsTPCpos; centrality; s_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisSN});
      registry.add("epReso/hEpCosCoefficientsTPCneg", "hEpCosCoefficientsTPCneg; centrality; c_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisCN});
      registry.add("epReso/hEpSinCoefficientsTPCneg", "hEpSinCoefficientsTPCneg; centrality; s_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisSN});
      registry.add("epReso/hEpCosCoefficientsTPCTots", "hEpCosCoefficientsTPCTots; centrality; c_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisCN});
      registry.add("epReso/hEpSinCoefficientsTPCTots", "hEpSinCoefficientsTPCTots; centrality; s_{n}", HistType::kTProfile2D, {thnAxisCent, thAxisSN});
      registry.add("QVector/hQVecMeanRVsPhiFT0a", "hQVecMeanRVsPhiFT0a; centrality; #it{#varphi} (rad), <#it{r}> (a.u.)", HistType::kTProfile2D, {thnAxisCent, thAxisAzimuth});
      registry.add("QVector/hQVecMeanRVsPhiFT0c", "hQVecMeanRVsPhiFT0c; centrality; #it{#varphi} (rad), <#it{r}> (a.u.)", HistType::kTProfile2D, {thnAxisCent, thAxisAzimuth});
      registry.add("QVector/hQVecMeanRVsPhiFT0m", "hQVecMeanRVsPhiFT0m; centrality; #it{#varphi} (rad), <#it{r}> (a.u.)", HistType::kTProfile2D, {thnAxisCent, thAxisAzimuth});
      registry.add("QVector/hQVecMeanRVsPhiFV0a", "hQVecMeanRVsPhiFV0a; centrality; #it{#varphi} (rad), <#it{r}> (a.u.)", HistType::kTProfile2D, {thnAxisCent, thAxisAzimuth});
      registry.add("QVector/hQVecMeanRVsPhiTPCpos", "hQVecMeanRVsPhiTPCpos; centrality; #it{#varphi} (rad), <#it{r}> (a.u.)", HistType::kTProfile2D, {thnAxisCent, thAxisAzimuth});
      registry.add("QVector/hQVecMeanRVsPhiTPCneg", "hQVecMeanRVsPhiTPCneg; centrality; #it{#varphi} (rad), <#it{r}> (a.u.)", HistType::kTProfile2D, {thnAxisCent, thAxisAzimuth});
      registry.add("QVector/hQVecMeanRVsPhiTPCTot", "hQVecMeanRVsPhiTPCTot; centrality; #it{#varphi} (rad), <#it{r}> (a.u.)", HistType::kTProfile2D, {thnAxisCent, thAxisAzimuth});
    }
  }; // end init

  /// Compute the delta psi in the range [0, pi/harmonic]
  /// \param psi1 is the first angle
  /// \param psi2 is the second angle
  float getDeltaPsiInRange(float psi1, float psi2)
  {
    float deltaPsi = psi1 - psi2;
    return RecoDecay::constrainAngle(deltaPsi, 0.f, harmonic);
  }

  /// Get the centrality
  /// \param collision is the collision with the centrality information
  template <o2::soa::is_iterator TCollision>
  float getCentrality(TCollision const& collision)
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

  /// Get all used Q vector
  /// \param collision is the collision with the Q vector information
  template <o2::soa::is_iterator TCollision>
  std::array<float, NQVecEntries> getAllQvec(TCollision const& collision)
  {
    // Retrieve the Q vectors using the helper function for each detector
    auto [xQVecMain, yQVecMain] = getQvec(collision, qvecDetector);
    auto [xQVecSubA, yQVecSubA] = getQvec(collision, qvecSubADetector);
    auto [xQVecSubB, yQVecSubB] = getQvec(collision, qvecSubBDetector);

    return {xQVecMain, yQVecMain, xQVecSubA, yQVecSubA, xQVecSubB, yQVecSubB};
  }

  /// Get the Q vector
  /// \param collision is the collision with the Q vector information
  template <o2::soa::is_iterator TCollision>
  std::pair<float, float> getQvec(TCollision const& collision, int detector)
  {
    float xQVec = -999.f;
    float yQVec = -999.f;

    switch (detector) {
      case QvecEstimator::FT0M:
        if (harmonic == kElliptic) {
          xQVec = collision.q2xft0m();
          yQVec = collision.q2yft0m();
        } else if (harmonic == kTriangluar) {
          xQVec = collision.q3xft0m();
          yQVec = collision.q3yft0m();
        }
        break;
      case QvecEstimator::FT0A:
        if (harmonic == kElliptic) {
          xQVec = collision.q2xft0a();
          yQVec = collision.q2yft0a();
        } else if (harmonic == kTriangluar) {
          xQVec = collision.q3xft0a();
          yQVec = collision.q3yft0a();
        }
        break;
      case QvecEstimator::FT0C:
        if (harmonic == kElliptic) {
          xQVec = collision.q2xft0c();
          yQVec = collision.q2yft0c();
        } else if (harmonic == kTriangluar) {
          xQVec = collision.q3xft0c();
          yQVec = collision.q3yft0c();
        }
        break;
      case QvecEstimator::TPCPos:
        if (harmonic == kElliptic) {
          xQVec = collision.q2xbpos();
          yQVec = collision.q2ybpos();
        } else if (harmonic == kTriangluar) {
          xQVec = collision.q3xbpos();
          yQVec = collision.q3ybpos();
        }
        break;
      case QvecEstimator::TPCNeg:
        if (harmonic == kElliptic) {
          xQVec = collision.q2xbneg();
          yQVec = collision.q2ybneg();
        } else if (harmonic == kTriangluar) {
          xQVec = collision.q3xbneg();
          yQVec = collision.q3ybneg();
        }
        break;
      case QvecEstimator::TPCTot:
        if (harmonic == kElliptic) {
          xQVec = collision.q2xbtot();
          yQVec = collision.q2ybtot();
        } else if (harmonic == kTriangluar) {
          xQVec = collision.q3xbtot();
          yQVec = collision.q3ybtot();
        }
        break;
      case QvecEstimator::FV0A:
        if (harmonic == kElliptic) {
          xQVec = collision.q2xfv0a();
          yQVec = collision.q2yfv0a();
        } else if (harmonic == kTriangluar) {
          xQVec = collision.q3xfv0a();
          yQVec = collision.q3yfv0a();
        }
        break;
      default:
        LOG(warning) << "Q vector estimator not valid. Falling back to FT0M";
        if (harmonic == kElliptic) {
          xQVec = collision.q2xft0m();
          yQVec = collision.q2yft0m();
        } else if (harmonic == kTriangluar) {
          xQVec = collision.q3xft0m();
          yQVec = collision.q3yft0m();
        }
        break;
    }
    return {xQVec, yQVec};
  }

  /// Check if the QVector values are within reasonable range
  /// \param collision is the collision with the Q vector information
  bool isQvecGood(std::array<float, NQVecEntries> const& QVecs)
  {
    bool isgood = true;
    for (const auto& QVec : QVecs) {
      if (std::fabs(QVec) > cfgMaxQVector) {
        isgood = false;
        break;
      }
    }
    return isgood;
  }

  // Resolution
  void processResolution(CollWithQvec const& collision)
  {
    o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<0>(&registry, collision);
    if (!(fEMEventCut.IsSelected(collision))) {
      // no selection on the centrality is applied on purpose to allow for the resolution study in post-processing
      return;
    }
    if (!(eventcuts.cfgFT0COccupancyMin <= collision.ft0cOccupancyInTimeRange() && collision.ft0cOccupancyInTimeRange() < eventcuts.cfgFT0COccupancyMax)) {
      return;
    }
    float cent = getCentrality(collision);
    if (cent < eventcuts.cfgMinCent || cent > eventcuts.cfgMaxCent) {
      // event selection
      return;
    }
    if (!isQvecGood(getAllQvec(collision))) {
      // selection based on QVector
      return;
    }
    o2::aod::pwgem::photonmeson::utils::eventhistogram::fillEventInfo<1>(&registry, collision);
    registry.fill(HIST("Event/before/hCollisionCounter"), 12.0); // accepted
    registry.fill(HIST("Event/after/hCollisionCounter"), 12.0);  // accepted

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
    float xQVecFV0a = -999.f;
    float yQVecFV0a = -999.f;
    if (harmonic == kElliptic) {
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
      xQVecFV0a = collision.q2xfv0a();
      yQVecFV0a = collision.q2yfv0a();
    } else if (harmonic == kTriangluar) {
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
      xQVecFV0a = collision.q3xfv0a();
      yQVecFV0a = collision.q3yfv0a();
    }

    if (saveSPResoHist) {
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
      registry.fill(HIST("spReso/hSpResoFV0aFT0c"), centrality, xQVecFV0a * xQVecFT0c + yQVecFV0a * yQVecFT0c);
      registry.fill(HIST("spReso/hSpResoFV0aTPCpos"), centrality, xQVecFV0a * xQVecBPos + yQVecFV0a * yQVecBPos);
      registry.fill(HIST("spReso/hSpResoFV0aTPCneg"), centrality, xQVecFV0a * xQVecBNeg + yQVecFV0a * yQVecBNeg);
      registry.fill(HIST("spReso/hSpResoFV0aTPCtot"), centrality, xQVecFV0a * xQVecBTot + yQVecFV0a * yQVecBTot);
    }

    if (saveEpResoHisto) {
      float epFT0a = epHelper.GetEventPlane(xQVecFT0a, yQVecFT0a, harmonic);
      float epFT0c = epHelper.GetEventPlane(xQVecFT0c, yQVecFT0c, harmonic);
      float epFT0m = epHelper.GetEventPlane(xQVecFT0m, yQVecFT0m, harmonic);
      float epBPoss = epHelper.GetEventPlane(xQVecBPos, yQVecBPos, harmonic);
      float epBNegs = epHelper.GetEventPlane(xQVecBNeg, yQVecBNeg, harmonic);
      float epBTots = epHelper.GetEventPlane(xQVecBTot, yQVecBTot, harmonic);
      float epFV0a = epHelper.GetEventPlane(xQVecFV0a, yQVecFV0a, harmonic);

      registry.fill(HIST("hEventPlaneAngleFT0M"), centrality, epFT0m);
      registry.fill(HIST("hEventPlaneAngleFT0A"), centrality, epFT0c);
      registry.fill(HIST("hEventPlaneAngleFT0C"), centrality, epFT0a);
      registry.fill(HIST("hEventPlaneAngleTPCpos"), centrality, epBPoss);
      registry.fill(HIST("hEventPlaneAngleTPCneg"), centrality, epBNegs);
      registry.fill(HIST("hEventPlaneAngleFV0A"), centrality, epFV0a);

      registry.fill(HIST("QVector/hQVecMeanRVsPhiFT0a"), centrality, std::atan2(yQVecFT0a, xQVecFT0a), std::hypot(xQVecFT0a, yQVecFT0a));
      registry.fill(HIST("QVector/hQVecMeanRVsPhiFT0c"), centrality, std::atan2(yQVecFT0c, xQVecFT0c), std::hypot(xQVecFT0c, yQVecFT0c));
      registry.fill(HIST("QVector/hQVecMeanRVsPhiFT0m"), centrality, std::atan2(yQVecFT0m, xQVecFT0m), std::hypot(xQVecFT0m, yQVecFT0m));
      registry.fill(HIST("QVector/hQVecMeanRVsPhiTPCpos"), centrality, std::atan2(yQVecBPos, xQVecBPos), std::hypot(xQVecBPos, yQVecBPos));
      registry.fill(HIST("QVector/hQVecMeanRVsPhiTPCneg"), centrality, std::atan2(yQVecBNeg, xQVecBNeg), std::hypot(xQVecBNeg, yQVecBNeg));
      registry.fill(HIST("QVector/hQVecMeanRVsPhiTPCTot"), centrality, std::atan2(yQVecBTot, xQVecBTot), std::hypot(xQVecBTot, yQVecBTot));
      registry.fill(HIST("QVector/hQVecMeanRVsPhiFV0a"), centrality, std::atan2(yQVecFV0a, xQVecFV0a), std::hypot(xQVecFV0a, yQVecFV0a));

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
      registry.fill(HIST("epReso/hEpResoFV0aFT0c"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFV0a, epFT0c)));
      registry.fill(HIST("epReso/hEpResoFV0aTPCpos"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFV0a, epBPoss)));
      registry.fill(HIST("epReso/hEpResoFV0aTPCneg"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFV0a, epBNegs)));
      registry.fill(HIST("epReso/hEpResoFV0aTPCtot"), centrality, std::cos(harmonic * getDeltaPsiInRange(epFV0a, epBTots)));
      for (int n = 1; n <= kOctagonal; n++) {
        registry.fill(HIST("epReso/hEpCosCoefficientsFT0c"), centrality, n, std::cos(n * harmonic * epFT0c));
        registry.fill(HIST("epReso/hEpSinCoefficientsFT0c"), centrality, n, std::sin(n * harmonic * epFT0c));
        registry.fill(HIST("epReso/hEpCosCoefficientsFT0a"), centrality, n, std::cos(n * harmonic * epFT0a));
        registry.fill(HIST("epReso/hEpSinCoefficientsFT0a"), centrality, n, std::sin(n * harmonic * epFT0a));
        registry.fill(HIST("epReso/hEpCosCoefficientsFT0m"), centrality, n, std::cos(n * harmonic * epFT0m));
        registry.fill(HIST("epReso/hEpSinCoefficientsFT0m"), centrality, n, std::sin(n * harmonic * epFT0m));
        registry.fill(HIST("epReso/hEpCosCoefficientsTPCpos"), centrality, n, std::cos(n * harmonic * epBPoss));
        registry.fill(HIST("epReso/hEpSinCoefficientsTPCpos"), centrality, n, std::sin(n * harmonic * epBPoss));
        registry.fill(HIST("epReso/hEpCosCoefficientsTPCneg"), centrality, n, std::cos(n * harmonic * epBNegs));
        registry.fill(HIST("epReso/hEpSinCoefficientsTPCneg"), centrality, n, std::sin(n * harmonic * epBNegs));
        registry.fill(HIST("epReso/hEpCosCoefficientsTPCTots"), centrality, n, std::cos(n * harmonic * epBTots));
        registry.fill(HIST("epReso/hEpSinCoefficientsTPCTots"), centrality, n, std::sin(n * harmonic * epBTots));
        registry.fill(HIST("epReso/hEpCosCoefficientsFV0a"), centrality, n, std::cos(n * harmonic * epFV0a));
        registry.fill(HIST("epReso/hEpSinCoefficientsFV0a"), centrality, n, std::sin(n * harmonic * epFV0a));
      }
    }
  }
  PROCESS_SWITCH(TaskFlowReso, processResolution, "Process resolution", true);

}; // End struct TaskFlowReso

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<TaskFlowReso>(cfgc)};
}
