// Copyright 2019-2020 CERN and copyright holders of the ALICE Collaboration.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".

/// \file f0980FlowSP.cxx
/// \brief pi+pi- invariant-mass analysis and f0(980) v2 with scalar-product and event-plane methods

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/Qvectors.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/ASoA.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <Math/Vector4D.h>

#include <array>
#include <cinttypes>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;

namespace
{
using P4 = ROOT::Math::PxPyPzMVector;
constexpr float InvalidQ = -900.f;
constexpr int ProtonMassNumber = 1;
constexpr int OxygenMassNumber = 16;
} // namespace

struct F0980FlowSP {
  SliceCache cache;
  HistogramRegistry registry{"registry"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  enum class CollisionSystem : int {
    PP = 0,
    OO = 1,
    PO = 2,
    Unsupported = 3
  };

  enum class PionPidMode : int {
    TPCWithTOFWhenAvailable = 0,
    TPCOnly = 1,
    TPCTOFRequired = 2,
    MomentumDependentTOF = 4
  };

  enum class EventSelectionStep : int {
    All,
    Sel8,
    NoSameBunchPileup,
    NoTimeFrameBorder,
    NoITSROFrameBorder,
    ITSTPCVertex,
    Centrality,
    QVectorCalibrated,
    QVectorValid,
    Count
  };

  enum class PPEventSelectionStep : int {
    All,
    Sel8,
    NoSameBunchPileup,
    NoTimeFrameBorder,
    NoITSROFrameBorder,
    ITSTPCVertex,
    Inelastic,
    QVectorValid,
    Count
  };

  enum class QVectorStatus {
    IndexUnavailable,
    Uncalibrated,
    Invalid,
    Valid
  };

  struct QVectorData {
    double qAx = 0.;
    double qAy = 0.;
    double qBx = 0.;
    double qBy = 0.;
    double qCx = 0.;
    double qCy = 0.;
    double qAMagnitude = 0.;
    double qBMagnitude = 0.;
    double qCMagnitude = 0.;

    bool hasValidEventPlanes() const
    {
      return qAMagnitude > 0. && qBMagnitude > 0. && qCMagnitude > 0.;
    }
  };

  int cachedRunNumber = -1;
  CollisionSystem cachedCollisionSystem = CollisionSystem::Unsupported;

  // Pair and daughter selection. Defaults follow rho770analysis.cxx.
  Configurable<float> cfgMinPt{"cfgMinPt", 0.15f, "Minimum daughter transverse momentum"};
  Configurable<float> cfgMaxEta{"cfgMaxEta", 0.8f, "Maximum daughter pseudorapidity"};
  Configurable<float> cfgMinRap{"cfgMinRap", -0.5f, "Minimum pair rapidity"};
  Configurable<float> cfgMaxRap{"cfgMaxRap", 0.5f, "Maximum pair rapidity"};
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Require primary track"};
  Configurable<bool> cfgGlobalWoDCA{"cfgGlobalWoDCA", true, "Require global track without DCA cuts"};
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "Require PV contributor"};
  Configurable<float> cfgMaxDCAxy{"cfgMaxDCAxy", 0.1f, "Maximum |DCAxy| (cm)"};
  Configurable<float> cfgMaxDCAz{"cfgMaxDCAz", 0.1f, "Maximum |DCAz| (cm)"};
  Configurable<float> cfgMaxTPCChi2NCl{"cfgMaxTPCChi2NCl", 4.f, "Maximum TPC chi2 per cluster"};
  Configurable<float> cfgMaxITSChi2NCl{"cfgMaxITSChi2NCl", 36.f, "Maximum ITS chi2 per cluster"};
  Configurable<int> cfgMinTPCRows{"cfgMinTPCRows", 70, "Minimum TPC crossed rows"};
  Configurable<int> cfgMinITSClusters{"cfgMinITSClusters", 5, "Minimum ITS clusters"};

  Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "Require sel8"};
  Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", true, "Require kNoSameBunchPileup"};
  Configurable<bool> cfgRequireNoTimeFrameBorder{"cfgRequireNoTimeFrameBorder", true, "Require kNoTimeFrameBorder"};
  Configurable<bool> cfgRequireNoITSROFrameBorder{"cfgRequireNoITSROFrameBorder", true, "Require kNoITSROFrameBorder"};
  Configurable<bool> cfgRequireITSTPCVertex{"cfgRequireITSTPCVertex", true, "Require kIsVertexITSTPC"};
  Configurable<int> cfgCollisionSystem{"cfgCollisionSystem", -1, "Collision system: -1 auto from GRPLHCIF, 0 pp, 1 OO, 2 pO"};

  Configurable<int> cfgPIDMode{"cfgPIDMode", 0, "Pion PID: 0 TPC plus TOF when available, 1 TPC only, 2 require both TPC and TOF, 4 TPC-only below momentum threshold and require TOF above it"};
  Configurable<float> cfgTPCNSigmaPi{"cfgTPCNSigmaPi", 3.f, "TPC pion n-sigma with TOF"};
  Configurable<float> cfgTPCNSigmaPiNoTOF{"cfgTPCNSigmaPiNoTOF", 3.f, "TPC pion n-sigma without TOF"};
  Configurable<float> cfgTOFNSigmaPi{"cfgTOFNSigmaPi", 3.f, "TOF pion n-sigma"};
  Configurable<float> cfgTOFRequiredAboveP{"cfgTOFRequiredAboveP", 0.6f, "For PID mode 4, require TOF at and above this total momentum (GeV/c)"};
  // The vector tables contain one entry per harmonic requested from qVectorsTable.
  Configurable<size_t> cfgHarmonicIndex{"cfgHarmonicIndex", 0, "Index of n=2 in Q-vector vector columns"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0.f, "Minimum FT0C centrality for OO"};
  Configurable<float> cfgCentMax{"cfgCentMax", 100.f, "Maximum FT0C centrality for OO"};

  ConfigurableAxis axisMass{"axisMass", {400, 0.2, 2.2}, "m_{pi pi} (GeV/c^{2})"};
  ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0., 0.5, 1., 1.5, 2., 2.5, 3., 4., 6., 8., 10., 13., 20.}, "pair pT (GeV/c)"};
  ConfigurableAxis axisCent{"axisCent", {VARIABLE_WIDTH, 0., 5., 10., 20., 30., 40., 50., 60., 70., 80., 100.}, "FT0C centrality (%)"};
  ConfigurableAxis axisMult{"axisMult", {VARIABLE_WIDTH, 0., 1., 5., 10., 20., 30., 40., 50., 60., 70., 80., 90., 100., 105.}, "N_{PV contributors}"};
  ConfigurableAxis axisRun{"axisRun", {2001, 563999.5, 566000.5}, "run number"};

  using PPCollisions = soa::Join<aod::Collisions, aod::EvSels,
                                 aod::QvectorFT0CVecs, aod::QvectorFT0AVecs,
                                 aod::QvectorTPCposVecs>;
  using OOCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs,
                                 aod::QvectorFT0CVecs, aod::QvectorFT0AVecs,
                                 aod::QvectorTPCposVecs>;
  using POCollisions = OOCollisions;
  using TracksTPC = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                              aod::TrackSelection, aod::pidTPCPi>;
  using TracksTPCTOF = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA,
                                 aod::TrackSelection, aod::pidTPCPi, aod::pidTOFPi,
                                 aod::pidTOFbeta>;

  template <typename Collision>
  CollisionSystem collisionSystem(Collision const& collision)
  {
    if (cfgCollisionSystem.value >= static_cast<int>(CollisionSystem::PP) && cfgCollisionSystem.value <= static_cast<int>(CollisionSystem::PO)) {
      return static_cast<CollisionSystem>(cfgCollisionSystem.value);
    }

    const auto bc = collision.template bc_as<aod::BCsWithTimestamps>();
    if (bc.runNumber() == cachedRunNumber) {
      return cachedCollisionSystem;
    }

    auto* grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", bc.timestamp());
    if (grplhcif == nullptr) {
      LOGF(fatal, "Cannot load GLO/Config/GRPLHCIF for run %d and timestamp %" PRIu64, bc.runNumber(), bc.timestamp());
    }

    const int beamAA = grplhcif->getBeamA(o2::constants::lhc::BeamA);
    const int beamCA = grplhcif->getBeamA(o2::constants::lhc::BeamC);
    if (beamAA == ProtonMassNumber && beamCA == ProtonMassNumber) {
      cachedCollisionSystem = CollisionSystem::PP;
    } else if (beamAA == OxygenMassNumber && beamCA == OxygenMassNumber) {
      cachedCollisionSystem = CollisionSystem::OO;
    } else if ((beamAA == ProtonMassNumber && beamCA == OxygenMassNumber) ||
               (beamAA == OxygenMassNumber && beamCA == ProtonMassNumber)) {
      cachedCollisionSystem = CollisionSystem::PO;
    } else {
      cachedCollisionSystem = CollisionSystem::Unsupported;
    }
    cachedRunNumber = bc.runNumber();
    LOGF(info, "Automatically identified collision system for run %d: beam A values %d and %d, mode %d", cachedRunNumber, beamAA, beamCA, static_cast<int>(cachedCollisionSystem));
    return cachedCollisionSystem;
  }

  template <typename Collision>
  QVectorStatus getQVectorData(Collision const& collision, QVectorData& qVector) const
  {
    const size_t harmonicIndex = cfgHarmonicIndex.value;
    if (harmonicIndex >= collision.qvecFT0CReVec().size() ||
        harmonicIndex >= collision.qvecFT0CImVec().size() ||
        harmonicIndex >= collision.qvecFT0AReVec().size() ||
        harmonicIndex >= collision.qvecFT0AImVec().size() ||
        harmonicIndex >= collision.qvecTPCposReVec().size() ||
        harmonicIndex >= collision.qvecTPCposImVec().size()) {
      return QVectorStatus::IndexUnavailable;
    }
    if (!collision.isCalibrated()) {
      return QVectorStatus::Uncalibrated;
    }

    qVector.qAx = collision.qvecFT0CReVec()[harmonicIndex];
    qVector.qAy = collision.qvecFT0CImVec()[harmonicIndex];
    qVector.qBx = collision.qvecFT0AReVec()[harmonicIndex];
    qVector.qBy = collision.qvecFT0AImVec()[harmonicIndex];
    qVector.qCx = collision.qvecTPCposReVec()[harmonicIndex];
    qVector.qCy = collision.qvecTPCposImVec()[harmonicIndex];
    if (qVector.qAx < InvalidQ || qVector.qAy < InvalidQ ||
        qVector.qBx < InvalidQ || qVector.qBy < InvalidQ ||
        qVector.qCx < InvalidQ || qVector.qCy < InvalidQ) {
      return QVectorStatus::Invalid;
    }

    qVector.qAMagnitude = std::hypot(qVector.qAx, qVector.qAy);
    qVector.qBMagnitude = std::hypot(qVector.qBx, qVector.qBy);
    qVector.qCMagnitude = std::hypot(qVector.qCx, qVector.qCy);
    return QVectorStatus::Valid;
  }

  void init(InitContext const&)
  {
    // Counts and weighted SP numerators have identical axes. Their ratio is
    // <u_2 q_A*> in each (mass, pT, centrality) bin.
    registry.add("PP/hMassUS", "pp unlike-sign;mass;pT;multiplicity", HistType::kTHnSparseF, {axisMass, axisPt, axisMult});
    registry.add("PP/hMassLSpp", "pp ++;mass;pT;multiplicity", HistType::kTHnSparseF, {axisMass, axisPt, axisMult});
    registry.add("PP/hMassLSmm", "pp --;mass;pT;multiplicity", HistType::kTHnSparseF, {axisMass, axisPt, axisMult});
    registry.add("PP/hSpNumUS", "sum Re(u2 qFT0C*), unlike-sign;mass;pT;multiplicity", HistType::kTHnSparseF, {axisMass, axisPt, axisMult});
    registry.add("PP/hSpNumLSpp", "sum Re(u2 qFT0C*), ++;mass;pT;multiplicity", HistType::kTHnSparseF, {axisMass, axisPt, axisMult});
    registry.add("PP/hSpNumLSmm", "sum Re(u2 qFT0C*), --;mass;pT;multiplicity", HistType::kTHnSparseF, {axisMass, axisPt, axisMult});
    registry.add("PP/hEpNumUS", "sum cos(2(phi-Psi2FT0C)), unlike-sign;mass;pT;multiplicity", HistType::kTHnSparseF, {axisMass, axisPt, axisMult});
    registry.add("PP/hEpNumLSpp", "sum cos(2(phi-Psi2FT0C)), ++;mass;pT;multiplicity", HistType::kTHnSparseF, {axisMass, axisPt, axisMult});
    registry.add("PP/hEpNumLSmm", "sum cos(2(phi-Psi2FT0C)), --;mass;pT;multiplicity", HistType::kTHnSparseF, {axisMass, axisPt, axisMult});
    registry.add("PP/pQAB", "<qFT0C qFT0A*>;multiplicity", HistType::kTProfile, {axisMult});
    registry.add("PP/pQAC", "<qFT0C qTPCpos*>;multiplicity", HistType::kTProfile, {axisMult});
    registry.add("PP/pQBC", "<qFT0A qTPCpos*>;multiplicity", HistType::kTProfile, {axisMult});
    registry.add("PP/pQABIm", "Im(qFT0C qFT0A*);multiplicity", HistType::kTProfile, {axisMult});
    registry.add("PP/pQACIm", "Im(qFT0C qTPCpos*);multiplicity", HistType::kTProfile, {axisMult});
    registry.add("PP/pQBCIm", "Im(qFT0A qTPCpos*);multiplicity", HistType::kTProfile, {axisMult});
    registry.add("PP/pCosPsiAB", "<cos(2(PsiFT0C-PsiFT0A))>;multiplicity", HistType::kTProfile, {axisMult});
    registry.add("PP/pCosPsiAC", "<cos(2(PsiFT0C-PsiTPCpos))>;multiplicity", HistType::kTProfile, {axisMult});
    registry.add("PP/pCosPsiBC", "<cos(2(PsiFT0A-PsiTPCpos))>;multiplicity", HistType::kTProfile, {axisMult});
    registry.add("PP/QA/Event/hSelection", "pp event-selection cut flow;selection step;events", HistType::kTH1F, {{static_cast<int>(PPEventSelectionStep::Count), -0.5, static_cast<double>(PPEventSelectionStep::Count) - 0.5}});
    registry.add("PP/QA/Event/hMultiplicityBefore", "pp events before selection;N_{PV contributors};events", HistType::kTH1F, {{106, -0.5, 105.5}});
    registry.add("PP/QA/Event/hMultiplicityAfter", "pp events after selection;N_{PV contributors};events", HistType::kTH1F, {{106, -0.5, 105.5}});
    registry.add("PP/QA/Event/hVertexZBefore", "pp events before selection;z_{vtx} (cm);events", HistType::kTH1F, {{240, -30., 30.}});
    registry.add("PP/QA/Event/hVertexZAfter", "pp events after selection;z_{vtx} (cm);events", HistType::kTH1F, {{240, -30., 30.}});
    registry.add("PP/QA/Qvector/hFT0C", "corrected q_{2,FT0C};q_{x};q_{y}", HistType::kTH2F, {{200, -1., 1.}, {200, -1., 1.}});
    registry.add("PP/QA/Qvector/hFT0A", "corrected q_{2,FT0A};q_{x};q_{y}", HistType::kTH2F, {{200, -1., 1.}, {200, -1., 1.}});
    registry.add("PP/QA/Qvector/hTPCpos", "corrected q_{2,TPCpos};q_{x};q_{y}", HistType::kTH2F, {{200, -1., 1.}, {200, -1., 1.}});
    registry.add("PP/QA/Qvector/hMagnitudeFT0C", "FT0C q_{2} magnitude;multiplicity;|q_{2}|", HistType::kTH2F, {axisMult, {200, 0., 1.5}});
    registry.add("PP/QA/Qvector/hMagnitudeFT0A", "FT0A q_{2} magnitude;multiplicity;|q_{2}|", HistType::kTH2F, {axisMult, {200, 0., 1.5}});
    registry.add("PP/QA/Qvector/hMagnitudeTPCpos", "TPCpos q_{2} magnitude;multiplicity;|q_{2}|", HistType::kTH2F, {axisMult, {200, 0., 1.5}});
    registry.add("PP/QA/Qvector/hPsi2FT0C", "FT0C event-plane angle;multiplicity;#Psi_{2}", HistType::kTH2F, {axisMult, {180, -1.5708, 1.5708}});
    registry.add("PP/QA/Qvector/hPsi2FT0A", "FT0A event-plane angle;multiplicity;#Psi_{2}", HistType::kTH2F, {axisMult, {180, -1.5708, 1.5708}});
    registry.add("PP/QA/Qvector/hPsi2TPCpos", "TPCpos event-plane angle;multiplicity;#Psi_{2}", HistType::kTH2F, {axisMult, {180, -1.5708, 1.5708}});
    registry.add("PP/QA/Qvector/hDotAB", "Re(q_{FT0C}q_{FT0A}^{*});multiplicity;q_{A}#upointq_{B}", HistType::kTH2F, {axisMult, {240, -1.2, 1.2}});
    registry.add("PP/QA/Qvector/hDotAC", "Re(q_{FT0C}q_{TPCpos}^{*});multiplicity;q_{A}#upointq_{C}", HistType::kTH2F, {axisMult, {240, -1.2, 1.2}});
    registry.add("PP/QA/Qvector/hDotBC", "Re(q_{FT0A}q_{TPCpos}^{*});multiplicity;q_{B}#upointq_{C}", HistType::kTH2F, {axisMult, {240, -1.2, 1.2}});
    registry.add("PP/QA/Qvector/pQxFT0CVsMult", "mean corrected q_{x,FT0C};multiplicity;<#it{q}_{x}>", HistType::kTProfile, {axisMult});
    registry.add("PP/QA/Qvector/pQyFT0CVsMult", "mean corrected q_{y,FT0C};multiplicity;<#it{q}_{y}>", HistType::kTProfile, {axisMult});
    registry.add("PP/QA/Qvector/pQxFT0AVsMult", "mean corrected q_{x,FT0A};multiplicity;<#it{q}_{x}>", HistType::kTProfile, {axisMult});
    registry.add("PP/QA/Qvector/pQyFT0AVsMult", "mean corrected q_{y,FT0A};multiplicity;<#it{q}_{y}>", HistType::kTProfile, {axisMult});
    registry.add("PP/QA/Qvector/pQxTPCposVsMult", "mean corrected q_{x,TPCpos};multiplicity;<#it{q}_{x}>", HistType::kTProfile, {axisMult});
    registry.add("PP/QA/Qvector/pQyTPCposVsMult", "mean corrected q_{y,TPCpos};multiplicity;<#it{q}_{y}>", HistType::kTProfile, {axisMult});
    registry.add("PP/QA/Qvector/pQxFT0CVsRun", "mean corrected q_{x,FT0C};run number;<#it{q}_{x}>", HistType::kTProfile, {axisRun});
    registry.add("PP/QA/Qvector/pQyFT0CVsRun", "mean corrected q_{y,FT0C};run number;<#it{q}_{y}>", HistType::kTProfile, {axisRun});
    registry.add("PP/QA/Qvector/pQxFT0AVsRun", "mean corrected q_{x,FT0A};run number;<#it{q}_{x}>", HistType::kTProfile, {axisRun});
    registry.add("PP/QA/Qvector/pQyFT0AVsRun", "mean corrected q_{y,FT0A};run number;<#it{q}_{y}>", HistType::kTProfile, {axisRun});
    registry.add("PP/QA/Qvector/pQxTPCposVsRun", "mean corrected q_{x,TPCpos};run number;<#it{q}_{x}>", HistType::kTProfile, {axisRun});
    registry.add("PP/QA/Qvector/pQyTPCposVsRun", "mean corrected q_{y,TPCpos};run number;<#it{q}_{y}>", HistType::kTProfile, {axisRun});

    registry.add("OO/hMassUS", "OO unlike-sign;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("OO/hMassLSpp", "OO ++;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("OO/hMassLSmm", "OO --;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("OO/hSpNumUS", "sum Re(u2 qFT0C*), unlike-sign;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("OO/hSpNumLSpp", "sum Re(u2 qFT0C*), ++;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("OO/hSpNumLSmm", "sum Re(u2 qFT0C*), --;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("OO/hEpNumUS", "sum cos(2(phi-Psi2FT0C)), unlike-sign;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("OO/hEpNumLSpp", "sum cos(2(phi-Psi2FT0C)), ++;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("OO/hEpNumLSmm", "sum cos(2(phi-Psi2FT0C)), --;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});

    registry.add("OO/pQAB", "<qFT0C qFT0A*>;centrality", HistType::kTProfile, {axisCent});
    registry.add("OO/pQAC", "<qFT0C qTPCpos*>;centrality", HistType::kTProfile, {axisCent});
    registry.add("OO/pQBC", "<qFT0A qTPCpos*>;centrality", HistType::kTProfile, {axisCent});
    registry.add("OO/pQABIm", "Im(qFT0C qFT0A*);centrality", HistType::kTProfile, {axisCent});
    registry.add("OO/pQACIm", "Im(qFT0C qTPCpos*);centrality", HistType::kTProfile, {axisCent});
    registry.add("OO/pQBCIm", "Im(qFT0A qTPCpos*);centrality", HistType::kTProfile, {axisCent});
    registry.add("OO/pCosPsiAB", "<cos(2(PsiFT0C-PsiFT0A))>;centrality", HistType::kTProfile, {axisCent});
    registry.add("OO/pCosPsiAC", "<cos(2(PsiFT0C-PsiTPCpos))>;centrality", HistType::kTProfile, {axisCent});
    registry.add("OO/pCosPsiBC", "<cos(2(PsiFT0A-PsiTPCpos))>;centrality", HistType::kTProfile, {axisCent});

    registry.add("PO/hMassUS", "pO unlike-sign;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("PO/hMassLSpp", "pO ++;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("PO/hMassLSmm", "pO --;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("PO/hSpNumUS", "sum Re(u2 qFT0C*), unlike-sign;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("PO/hSpNumLSpp", "sum Re(u2 qFT0C*), ++;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("PO/hSpNumLSmm", "sum Re(u2 qFT0C*), --;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("PO/hEpNumUS", "sum cos(2(phi-Psi2FT0C)), unlike-sign;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("PO/hEpNumLSpp", "sum cos(2(phi-Psi2FT0C)), ++;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("PO/hEpNumLSmm", "sum cos(2(phi-Psi2FT0C)), --;mass;pT;centrality", HistType::kTHnSparseF, {axisMass, axisPt, axisCent});
    registry.add("PO/pQAB", "<qFT0C qFT0A*>;centrality", HistType::kTProfile, {axisCent});
    registry.add("PO/pQAC", "<qFT0C qTPCpos*>;centrality", HistType::kTProfile, {axisCent});
    registry.add("PO/pQBC", "<qFT0A qTPCpos*>;centrality", HistType::kTProfile, {axisCent});
    registry.add("PO/pQABIm", "Im(qFT0C qFT0A*);centrality", HistType::kTProfile, {axisCent});
    registry.add("PO/pQACIm", "Im(qFT0C qTPCpos*);centrality", HistType::kTProfile, {axisCent});
    registry.add("PO/pQBCIm", "Im(qFT0A qTPCpos*);centrality", HistType::kTProfile, {axisCent});
    registry.add("PO/pCosPsiAB", "<cos(2(PsiFT0C-PsiFT0A))>;centrality", HistType::kTProfile, {axisCent});
    registry.add("PO/pCosPsiAC", "<cos(2(PsiFT0C-PsiTPCpos))>;centrality", HistType::kTProfile, {axisCent});
    registry.add("PO/pCosPsiBC", "<cos(2(PsiFT0A-PsiTPCpos))>;centrality", HistType::kTProfile, {axisCent});

    registry.add("QA/Event/hSelection", "ion-collision event-selection cut flow;selection step;events", HistType::kTH1F, {{static_cast<int>(EventSelectionStep::Count), -0.5, static_cast<double>(EventSelectionStep::Count) - 0.5}});
    registry.add("QA/Event/hCentrality", "selected ion collisions;FT0C centrality (%);events", HistType::kTH1F, {{100, 0., 100.}});
    registry.add("QA/Event/hCentralityBefore", "ion-collision events before selection;FT0C centrality (%);events", HistType::kTH1F, {{100, 0., 100.}});
    registry.add("QA/Event/hCentralityAfter", "ion-collision events after selection;FT0C centrality (%);events", HistType::kTH1F, {{100, 0., 100.}});
    registry.add("QA/Event/hVertexZBefore", "ion-collision events before selection;z_{vtx} (cm);events", HistType::kTH1F, {{240, -30., 30.}});
    registry.add("QA/Event/hVertexZAfter", "ion-collision events after selection;z_{vtx} (cm);events", HistType::kTH1F, {{240, -30., 30.}});

    registry.add("QA/Qvector/hFT0C", "corrected q_{2,FT0C};q_{x};q_{y}", HistType::kTH2F, {{200, -1., 1.}, {200, -1., 1.}});
    registry.add("QA/Qvector/hFT0A", "corrected q_{2,FT0A};q_{x};q_{y}", HistType::kTH2F, {{200, -1., 1.}, {200, -1., 1.}});
    registry.add("QA/Qvector/hTPCpos", "corrected q_{2,TPCpos};q_{x};q_{y}", HistType::kTH2F, {{200, -1., 1.}, {200, -1., 1.}});
    registry.add("QA/Qvector/hMagnitudeFT0C", "FT0C q_{2} magnitude;FT0C centrality (%);|q_{2}|", HistType::kTH2F, {axisCent, {200, 0., 1.5}});
    registry.add("QA/Qvector/hMagnitudeFT0A", "FT0A q_{2} magnitude;FT0C centrality (%);|q_{2}|", HistType::kTH2F, {axisCent, {200, 0., 1.5}});
    registry.add("QA/Qvector/hMagnitudeTPCpos", "TPCpos q_{2} magnitude;FT0C centrality (%);|q_{2}|", HistType::kTH2F, {axisCent, {200, 0., 1.5}});
    registry.add("QA/Qvector/hPsi2FT0C", "FT0C event-plane angle;FT0C centrality (%);#Psi_{2}", HistType::kTH2F, {axisCent, {180, -1.5708, 1.5708}});
    registry.add("QA/Qvector/hPsi2FT0A", "FT0A event-plane angle;FT0C centrality (%);#Psi_{2}", HistType::kTH2F, {axisCent, {180, -1.5708, 1.5708}});
    registry.add("QA/Qvector/hPsi2TPCpos", "TPCpos event-plane angle;FT0C centrality (%);#Psi_{2}", HistType::kTH2F, {axisCent, {180, -1.5708, 1.5708}});
    registry.add("QA/Qvector/hDotAB", "Re(q_{FT0C}q_{FT0A}^{*});FT0C centrality (%);q_{A}#upointq_{B}", HistType::kTH2F, {axisCent, {240, -1.2, 1.2}});
    registry.add("QA/Qvector/hDotAC", "Re(q_{FT0C}q_{TPCpos}^{*});FT0C centrality (%);q_{A}#upointq_{C}", HistType::kTH2F, {axisCent, {240, -1.2, 1.2}});
    registry.add("QA/Qvector/hDotBC", "Re(q_{FT0A}q_{TPCpos}^{*});FT0C centrality (%);q_{B}#upointq_{C}", HistType::kTH2F, {axisCent, {240, -1.2, 1.2}});
    registry.add("QA/Qvector/pQxFT0CVsCent", "mean corrected q_{x,FT0C};FT0C centrality (%);<#it{q}_{x}>", HistType::kTProfile, {axisCent});
    registry.add("QA/Qvector/pQyFT0CVsCent", "mean corrected q_{y,FT0C};FT0C centrality (%);<#it{q}_{y}>", HistType::kTProfile, {axisCent});
    registry.add("QA/Qvector/pQxFT0AVsCent", "mean corrected q_{x,FT0A};FT0C centrality (%);<#it{q}_{x}>", HistType::kTProfile, {axisCent});
    registry.add("QA/Qvector/pQyFT0AVsCent", "mean corrected q_{y,FT0A};FT0C centrality (%);<#it{q}_{y}>", HistType::kTProfile, {axisCent});
    registry.add("QA/Qvector/pQxTPCposVsCent", "mean corrected q_{x,TPCpos};FT0C centrality (%);<#it{q}_{x}>", HistType::kTProfile, {axisCent});
    registry.add("QA/Qvector/pQyTPCposVsCent", "mean corrected q_{y,TPCpos};FT0C centrality (%);<#it{q}_{y}>", HistType::kTProfile, {axisCent});
    registry.add("QA/Qvector/pQxFT0CVsRun", "mean corrected q_{x,FT0C};run number;<#it{q}_{x}>", HistType::kTProfile, {axisRun});
    registry.add("QA/Qvector/pQyFT0CVsRun", "mean corrected q_{y,FT0C};run number;<#it{q}_{y}>", HistType::kTProfile, {axisRun});
    registry.add("QA/Qvector/pQxFT0AVsRun", "mean corrected q_{x,FT0A};run number;<#it{q}_{x}>", HistType::kTProfile, {axisRun});
    registry.add("QA/Qvector/pQyFT0AVsRun", "mean corrected q_{y,FT0A};run number;<#it{q}_{y}>", HistType::kTProfile, {axisRun});
    registry.add("QA/Qvector/pQxTPCposVsRun", "mean corrected q_{x,TPCpos};run number;<#it{q}_{x}>", HistType::kTProfile, {axisRun});
    registry.add("QA/Qvector/pQyTPCposVsRun", "mean corrected q_{y,TPCpos};run number;<#it{q}_{y}>", HistType::kTProfile, {axisRun});

    auto eventSelection = registry.get<TH1>(HIST("QA/Event/hSelection"));
    constexpr std::array eventSelectionLabels = {
      "all", "sel8", "no same-bunch pileup", "no TF border",
      "no ITS ROF border", "ITS-TPC vertex", "centrality", "Q calibrated", "Q valid"};
    for (size_t i = 0; i < eventSelectionLabels.size(); ++i) {
      eventSelection->GetXaxis()->SetBinLabel(i + 1, eventSelectionLabels[i]);
    }
    auto ppEventSelection = registry.get<TH1>(HIST("PP/QA/Event/hSelection"));
    constexpr std::array ppEventSelectionLabels = {
      "all", "sel8", "no same-bunch pileup", "no TF border",
      "no ITS ROF border", "ITS-TPC vertex", "INEL > 0", "Q valid"};
    for (size_t i = 0; i < ppEventSelectionLabels.size(); ++i) {
      ppEventSelection->GetXaxis()->SetBinLabel(i + 1, ppEventSelectionLabels[i]);
    }
    registry.add("QA/TrackSelected/hPt", "track-selected candidates;p_{T} (GeV/c);counts", HistType::kTH1F, {{200, 0., 20.}});
    registry.add("QA/TrackSelected/hEta", "track-selected candidates;#eta;counts", HistType::kTH1F, {{160, -0.8, 0.8}});
    registry.add("QA/TrackSelected/hDCAxy", "track-selected candidates;DCA_{xy} (cm);counts", HistType::kTH1F, {{200, -0.1, 0.1}});
    registry.add("QA/TrackSelected/hDCAz", "track-selected candidates;DCA_{z} (cm);counts", HistType::kTH1F, {{200, -0.1, 0.1}});
    registry.add("QA/TrackSelected/hTPCCrossedRows", "track-selected candidates;TPC crossed rows;counts", HistType::kTH1F, {{165, -0.5, 164.5}});
    registry.add("QA/TrackSelected/hTPCNSigmaPiVsP", "track-selected candidates before pion PID;p (GeV/c);n#sigma_{TPC}^{#pi}", HistType::kTH2F, {{200, 0., 20.}, {120, -6., 6.}});
    registry.add("QA/TrackSelected/hTPCSignalVsP", "track-selected candidates before pion PID;p (GeV/c);TPC dE/dx", HistType::kTH2F, {{200, 0., 20.}, {500, 0., 500.}});
    registry.add("QA/TrackSelected/hTOFNSigmaPiVsP", "track-selected candidates with TOF before pion PID;p (GeV/c);n#sigma_{TOF}^{#pi}", HistType::kTH2F, {{200, 0., 20.}, {120, -6., 6.}});
    registry.add("QA/TrackSelected/hTOFBetaVsP", "track-selected candidates with TOF before pion PID;p (GeV/c);#beta_{TOF}", HistType::kTH2F, {{200, 0., 20.}, {300, 0., 1.5}});
    registry.add("QA/SelectedPions/hTPCNSigmaPiVsP", "selected pion candidates;p (GeV/c);n#sigma_{TPC}^{#pi}", HistType::kTH2F, {{200, 0., 20.}, {120, -6., 6.}});
    registry.add("QA/SelectedPions/hTPCSignalVsP", "selected pion candidates after pion PID;p (GeV/c);TPC dE/dx", HistType::kTH2F, {{200, 0., 20.}, {500, 0., 500.}});
    registry.add("QA/SelectedPions/hTOFNSigmaPiVsP", "selected pion candidates with TOF after pion PID;p (GeV/c);n#sigma_{TOF}^{#pi}", HistType::kTH2F, {{200, 0., 20.}, {120, -6., 6.}});
    registry.add("QA/SelectedPions/hTOFBetaVsP", "selected pion candidates with TOF after pion PID;p (GeV/c);#beta_{TOF}", HistType::kTH2F, {{200, 0., 20.}, {300, 0., 1.5}});
  }

  template <typename T>
  bool selectTrack(T const& track) const
  {
    if (track.pt() < cfgMinPt.value || std::abs(track.eta()) > cfgMaxEta.value) {
      return false;
    }
    if (std::abs(track.dcaXY()) >= cfgMaxDCAxy.value || std::abs(track.dcaZ()) >= cfgMaxDCAz.value) {
      return false;
    }
    if (cfgPrimaryTrack.value && !track.isPrimaryTrack()) {
      return false;
    }
    if (cfgGlobalWoDCA.value && !track.isGlobalTrackWoDCA()) {
      return false;
    }
    if (cfgPVContributor.value && !track.isPVContributor()) {
      return false;
    }
    return track.tpcChi2NCl() < cfgMaxTPCChi2NCl.value &&
           track.itsChi2NCl() < cfgMaxITSChi2NCl.value &&
           track.tpcNClsCrossedRows() > cfgMinTPCRows.value &&
           track.itsNCls() > cfgMinITSClusters.value;
  }

  template <typename C>
  bool selectEvent(C const& collision) const
  {
    return (!cfgRequireSel8.value || collision.sel8()) &&
           (!cfgRequireNoSameBunchPileup.value || collision.selection_bit(aod::evsel::kNoSameBunchPileup)) &&
           (!cfgRequireNoTimeFrameBorder.value || collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) &&
           (!cfgRequireNoITSROFrameBorder.value || collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) &&
           (!cfgRequireITSTPCVertex.value || collision.selection_bit(aod::evsel::kIsVertexITSTPC));
  }

  template <typename C>
  bool selectOOEventWithQA(C const& collision, float centrality)
  {
    registry.fill(HIST("QA/Event/hSelection"), static_cast<int>(EventSelectionStep::All));
    registry.fill(HIST("QA/Event/hCentralityBefore"), centrality);
    registry.fill(HIST("QA/Event/hVertexZBefore"), collision.posZ());
    if (cfgRequireSel8.value && !collision.sel8()) {
      return false;
    }
    registry.fill(HIST("QA/Event/hSelection"), static_cast<int>(EventSelectionStep::Sel8));
    if (cfgRequireNoSameBunchPileup.value && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    registry.fill(HIST("QA/Event/hSelection"), static_cast<int>(EventSelectionStep::NoSameBunchPileup));
    if (cfgRequireNoTimeFrameBorder.value && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    registry.fill(HIST("QA/Event/hSelection"), static_cast<int>(EventSelectionStep::NoTimeFrameBorder));
    if (cfgRequireNoITSROFrameBorder.value && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    registry.fill(HIST("QA/Event/hSelection"), static_cast<int>(EventSelectionStep::NoITSROFrameBorder));
    if (cfgRequireITSTPCVertex.value && !collision.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    registry.fill(HIST("QA/Event/hSelection"), static_cast<int>(EventSelectionStep::ITSTPCVertex));
    if (centrality < cfgCentMin.value || centrality >= cfgCentMax.value) {
      return false;
    }
    registry.fill(HIST("QA/Event/hSelection"), static_cast<int>(EventSelectionStep::Centrality));
    registry.fill(HIST("QA/Event/hCentrality"), centrality);
    registry.fill(HIST("QA/Event/hCentralityAfter"), centrality);
    registry.fill(HIST("QA/Event/hVertexZAfter"), collision.posZ());
    return true;
  }

  void fillQVectorQA(float centrality, int runNumber,
                     double qAx, double qAy,
                     double qBx, double qBy,
                     double qCx, double qCy)
  {
    registry.fill(HIST("QA/Qvector/hFT0C"), qAx, qAy);
    registry.fill(HIST("QA/Qvector/hFT0A"), qBx, qBy);
    registry.fill(HIST("QA/Qvector/hTPCpos"), qCx, qCy);
    registry.fill(HIST("QA/Qvector/hMagnitudeFT0C"), centrality, std::hypot(qAx, qAy));
    registry.fill(HIST("QA/Qvector/hMagnitudeFT0A"), centrality, std::hypot(qBx, qBy));
    registry.fill(HIST("QA/Qvector/hMagnitudeTPCpos"), centrality, std::hypot(qCx, qCy));
    registry.fill(HIST("QA/Qvector/hPsi2FT0C"), centrality, 0.5 * std::atan2(qAy, qAx));
    registry.fill(HIST("QA/Qvector/hPsi2FT0A"), centrality, 0.5 * std::atan2(qBy, qBx));
    registry.fill(HIST("QA/Qvector/hPsi2TPCpos"), centrality, 0.5 * std::atan2(qCy, qCx));
    registry.fill(HIST("QA/Qvector/hDotAB"), centrality, qAx * qBx + qAy * qBy);
    registry.fill(HIST("QA/Qvector/hDotAC"), centrality, qAx * qCx + qAy * qCy);
    registry.fill(HIST("QA/Qvector/hDotBC"), centrality, qBx * qCx + qBy * qCy);
    registry.fill(HIST("QA/Qvector/pQxFT0CVsCent"), centrality, qAx);
    registry.fill(HIST("QA/Qvector/pQyFT0CVsCent"), centrality, qAy);
    registry.fill(HIST("QA/Qvector/pQxFT0AVsCent"), centrality, qBx);
    registry.fill(HIST("QA/Qvector/pQyFT0AVsCent"), centrality, qBy);
    registry.fill(HIST("QA/Qvector/pQxTPCposVsCent"), centrality, qCx);
    registry.fill(HIST("QA/Qvector/pQyTPCposVsCent"), centrality, qCy);
    registry.fill(HIST("QA/Qvector/pQxFT0CVsRun"), runNumber, qAx);
    registry.fill(HIST("QA/Qvector/pQyFT0CVsRun"), runNumber, qAy);
    registry.fill(HIST("QA/Qvector/pQxFT0AVsRun"), runNumber, qBx);
    registry.fill(HIST("QA/Qvector/pQyFT0AVsRun"), runNumber, qBy);
    registry.fill(HIST("QA/Qvector/pQxTPCposVsRun"), runNumber, qCx);
    registry.fill(HIST("QA/Qvector/pQyTPCposVsRun"), runNumber, qCy);
  }

  void fillPPQVectorQA(float multiplicity, int runNumber,
                       double qAx, double qAy,
                       double qBx, double qBy,
                       double qCx, double qCy)
  {
    registry.fill(HIST("PP/QA/Qvector/hFT0C"), qAx, qAy);
    registry.fill(HIST("PP/QA/Qvector/hFT0A"), qBx, qBy);
    registry.fill(HIST("PP/QA/Qvector/hTPCpos"), qCx, qCy);
    registry.fill(HIST("PP/QA/Qvector/hMagnitudeFT0C"), multiplicity, std::hypot(qAx, qAy));
    registry.fill(HIST("PP/QA/Qvector/hMagnitudeFT0A"), multiplicity, std::hypot(qBx, qBy));
    registry.fill(HIST("PP/QA/Qvector/hMagnitudeTPCpos"), multiplicity, std::hypot(qCx, qCy));
    registry.fill(HIST("PP/QA/Qvector/hPsi2FT0C"), multiplicity, 0.5 * std::atan2(qAy, qAx));
    registry.fill(HIST("PP/QA/Qvector/hPsi2FT0A"), multiplicity, 0.5 * std::atan2(qBy, qBx));
    registry.fill(HIST("PP/QA/Qvector/hPsi2TPCpos"), multiplicity, 0.5 * std::atan2(qCy, qCx));
    registry.fill(HIST("PP/QA/Qvector/hDotAB"), multiplicity, qAx * qBx + qAy * qBy);
    registry.fill(HIST("PP/QA/Qvector/hDotAC"), multiplicity, qAx * qCx + qAy * qCy);
    registry.fill(HIST("PP/QA/Qvector/hDotBC"), multiplicity, qBx * qCx + qBy * qCy);
    registry.fill(HIST("PP/QA/Qvector/pQxFT0CVsMult"), multiplicity, qAx);
    registry.fill(HIST("PP/QA/Qvector/pQyFT0CVsMult"), multiplicity, qAy);
    registry.fill(HIST("PP/QA/Qvector/pQxFT0AVsMult"), multiplicity, qBx);
    registry.fill(HIST("PP/QA/Qvector/pQyFT0AVsMult"), multiplicity, qBy);
    registry.fill(HIST("PP/QA/Qvector/pQxTPCposVsMult"), multiplicity, qCx);
    registry.fill(HIST("PP/QA/Qvector/pQyTPCposVsMult"), multiplicity, qCy);
    registry.fill(HIST("PP/QA/Qvector/pQxFT0CVsRun"), runNumber, qAx);
    registry.fill(HIST("PP/QA/Qvector/pQyFT0CVsRun"), runNumber, qAy);
    registry.fill(HIST("PP/QA/Qvector/pQxFT0AVsRun"), runNumber, qBx);
    registry.fill(HIST("PP/QA/Qvector/pQyFT0AVsRun"), runNumber, qBy);
    registry.fill(HIST("PP/QA/Qvector/pQxTPCposVsRun"), runNumber, qCx);
    registry.fill(HIST("PP/QA/Qvector/pQyTPCposVsRun"), runNumber, qCy);
  }

  template <typename T>
  bool selectPion(T const& track) const
  {
    const float tpc = std::abs(track.tpcNSigmaPi());
    const float tof = std::abs(track.tofNSigmaPi());
    switch (static_cast<PionPidMode>(cfgPIDMode.value)) {
      case PionPidMode::TPCWithTOFWhenAvailable:
        if (track.hasTOF()) {
          return tpc < cfgTPCNSigmaPi.value &&
                 tof < cfgTOFNSigmaPi.value;
        }
        return tpc < cfgTPCNSigmaPiNoTOF.value;
      case PionPidMode::TPCOnly:
        return tpc < cfgTPCNSigmaPiNoTOF.value;
      case PionPidMode::TPCTOFRequired:
        return track.hasTOF() &&
               tpc < cfgTPCNSigmaPi.value &&
               tof < cfgTOFNSigmaPi.value;
      case PionPidMode::MomentumDependentTOF:
        if (track.p() < cfgTOFRequiredAboveP.value) {
          return tpc < cfgTPCNSigmaPiNoTOF.value;
        }
        return track.hasTOF() &&
               tpc < cfgTPCNSigmaPi.value &&
               tof < cfgTOFNSigmaPi.value;
      default:
        return false;
    }
  }

  template <typename T>
  bool selectPionTPCOnly(T const& track) const
  {
    return std::abs(track.tpcNSigmaPi()) < cfgTPCNSigmaPiNoTOF.value;
  }

  template <typename Tracks, typename TrackSelect, typename PionSelect, typename FillPair>
  void loopPairsSelected(Tracks const& tracks, TrackSelect&& trackSelect, PionSelect&& pionSelect, FillPair&& fillPair)
  {
    struct SelectedPion {
      float px = 0.f;
      float py = 0.f;
      float pz = 0.f;
      int8_t sign = 0;
    };

    // Apply daughter-quality and PID selections once per track before building
    // pairs.  In OO events this avoids constructing all unfiltered O(N^2)
    // combinations and repeatedly evaluating the same track cuts.
    std::vector<SelectedPion> pions;
    pions.reserve(tracks.size());
    for (auto const& track : tracks) {
      if (!trackSelect(track)) {
        continue;
      }
      registry.fill(HIST("QA/TrackSelected/hPt"), track.pt());
      registry.fill(HIST("QA/TrackSelected/hEta"), track.eta());
      registry.fill(HIST("QA/TrackSelected/hDCAxy"), track.dcaXY());
      registry.fill(HIST("QA/TrackSelected/hDCAz"), track.dcaZ());
      registry.fill(HIST("QA/TrackSelected/hTPCCrossedRows"), track.tpcNClsCrossedRows());
      registry.fill(HIST("QA/TrackSelected/hTPCNSigmaPiVsP"), track.p(), track.tpcNSigmaPi());
      registry.fill(HIST("QA/TrackSelected/hTPCSignalVsP"), track.p(), track.tpcSignal());
      if constexpr (requires { track.tofNSigmaPi(); }) {
        if (track.hasTOF()) {
          registry.fill(HIST("QA/TrackSelected/hTOFNSigmaPiVsP"), track.p(), track.tofNSigmaPi());
          registry.fill(HIST("QA/TrackSelected/hTOFBetaVsP"), track.p(), track.beta());
        }
      }
      if (!pionSelect(track)) {
        continue;
      }
      registry.fill(HIST("QA/SelectedPions/hTPCNSigmaPiVsP"), track.p(), track.tpcNSigmaPi());
      registry.fill(HIST("QA/SelectedPions/hTPCSignalVsP"), track.p(), track.tpcSignal());
      if constexpr (requires { track.tofNSigmaPi(); }) {
        if (track.hasTOF()) {
          registry.fill(HIST("QA/SelectedPions/hTOFNSigmaPiVsP"), track.p(), track.tofNSigmaPi());
          registry.fill(HIST("QA/SelectedPions/hTOFBetaVsP"), track.p(), track.beta());
        }
      }
      pions.push_back({track.px(), track.py(), track.pz(), static_cast<int8_t>(track.sign())});
    }

    const double mpi = o2::constants::physics::MassPionCharged;
    for (size_t i = 0; i < pions.size(); ++i) {
      for (size_t j = i + 1; j < pions.size(); ++j) {
        auto const& track1 = pions[i];
        auto const& track2 = pions[j];
        P4 pion1(track1.px, track1.py, track1.pz, mpi);
        P4 pion2(track2.px, track2.py, track2.pz, mpi);
        auto pair = pion1 + pion2;
        if (pair.Rapidity() < cfgMinRap.value || pair.Rapidity() > cfgMaxRap.value) {
          continue;
        }
        fillPair(track1.sign, track2.sign, pair);
      }
    }
  }

  template <typename Tracks, typename FillPair>
  void loopPairs(Tracks const& tracks, FillPair&& fillPair)
  {
    loopPairsSelected(tracks, [&](auto const& track) { return selectTrack(track); }, [&](auto const& track) { return selectPion(track); }, std::forward<FillPair>(fillPair));
  }

  void processPP(PPCollisions::iterator const& collision, TracksTPCTOF const& tracks,
                 aod::BCsWithTimestamps const&)
  {
    if (collisionSystem(collision) != CollisionSystem::PP) {
      return;
    }
    const float multiplicity = collision.numContrib();
    registry.fill(HIST("PP/QA/Event/hSelection"), static_cast<int>(PPEventSelectionStep::All));
    registry.fill(HIST("PP/QA/Event/hMultiplicityBefore"), multiplicity);
    registry.fill(HIST("PP/QA/Event/hVertexZBefore"), collision.posZ());
    if (cfgRequireSel8.value && !collision.sel8()) {
      return;
    }
    registry.fill(HIST("PP/QA/Event/hSelection"), static_cast<int>(PPEventSelectionStep::Sel8));
    if (cfgRequireNoSameBunchPileup.value && !collision.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return;
    }
    registry.fill(HIST("PP/QA/Event/hSelection"), static_cast<int>(PPEventSelectionStep::NoSameBunchPileup));
    if (cfgRequireNoTimeFrameBorder.value && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    registry.fill(HIST("PP/QA/Event/hSelection"), static_cast<int>(PPEventSelectionStep::NoTimeFrameBorder));
    if (cfgRequireNoITSROFrameBorder.value && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return;
    }
    registry.fill(HIST("PP/QA/Event/hSelection"), static_cast<int>(PPEventSelectionStep::NoITSROFrameBorder));
    if (cfgRequireITSTPCVertex.value && !collision.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return;
    }
    registry.fill(HIST("PP/QA/Event/hSelection"), static_cast<int>(PPEventSelectionStep::ITSTPCVertex));
    if (collision.numContrib() < 1) {
      return;
    }
    registry.fill(HIST("PP/QA/Event/hSelection"), static_cast<int>(PPEventSelectionStep::Inelastic));
    registry.fill(HIST("PP/QA/Event/hMultiplicityAfter"), multiplicity);
    registry.fill(HIST("PP/QA/Event/hVertexZAfter"), collision.posZ());

    QVectorData qVector;
    if (getQVectorData(collision, qVector) != QVectorStatus::Valid) {
      return;
    }
    registry.fill(HIST("PP/QA/Event/hSelection"), static_cast<int>(PPEventSelectionStep::QVectorValid));

    const auto [qAx, qAy, qBx, qBy, qCx, qCy, qAMagnitude, qBMagnitude, qCMagnitude] = qVector;

    const int runNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    fillPPQVectorQA(multiplicity, runNumber, qAx, qAy, qBx, qBy, qCx, qCy);

    registry.fill(HIST("PP/pQAB"), multiplicity, qAx * qBx + qAy * qBy);
    registry.fill(HIST("PP/pQAC"), multiplicity, qAx * qCx + qAy * qCy);
    registry.fill(HIST("PP/pQBC"), multiplicity, qBx * qCx + qBy * qCy);
    registry.fill(HIST("PP/pQABIm"), multiplicity, qAy * qBx - qAx * qBy);
    registry.fill(HIST("PP/pQACIm"), multiplicity, qAy * qCx - qAx * qCy);
    registry.fill(HIST("PP/pQBCIm"), multiplicity, qBy * qCx - qBx * qCy);
    const bool validEventPlanes = qVector.hasValidEventPlanes();
    if (validEventPlanes) {
      registry.fill(HIST("PP/pCosPsiAB"), multiplicity, (qAx * qBx + qAy * qBy) / (qAMagnitude * qBMagnitude));
      registry.fill(HIST("PP/pCosPsiAC"), multiplicity, (qAx * qCx + qAy * qCy) / (qAMagnitude * qCMagnitude));
      registry.fill(HIST("PP/pCosPsiBC"), multiplicity, (qBx * qCx + qBy * qCy) / (qBMagnitude * qCMagnitude));
    }

    loopPairs(tracks, [&](int8_t sign1, int8_t sign2, P4 const& pair) {
      const double ux = std::cos(2. * pair.Phi());
      const double uy = std::sin(2. * pair.Phi());
      const double sp = ux * qAx + uy * qAy;
      const double ep = validEventPlanes ? sp / qAMagnitude : 0.;
      if (sign1 * sign2 < 0) {
        registry.fill(HIST("PP/hMassUS"), pair.M(), pair.Pt(), multiplicity);
        registry.fill(HIST("PP/hSpNumUS"), pair.M(), pair.Pt(), multiplicity, sp);
        if (validEventPlanes) {
          registry.fill(HIST("PP/hEpNumUS"), pair.M(), pair.Pt(), multiplicity, ep);
        }
      } else if (sign1 > 0 && sign2 > 0) {
        registry.fill(HIST("PP/hMassLSpp"), pair.M(), pair.Pt(), multiplicity);
        registry.fill(HIST("PP/hSpNumLSpp"), pair.M(), pair.Pt(), multiplicity, sp);
        if (validEventPlanes) {
          registry.fill(HIST("PP/hEpNumLSpp"), pair.M(), pair.Pt(), multiplicity, ep);
        }
      } else if (sign1 < 0 && sign2 < 0) {
        registry.fill(HIST("PP/hMassLSmm"), pair.M(), pair.Pt(), multiplicity);
        registry.fill(HIST("PP/hSpNumLSmm"), pair.M(), pair.Pt(), multiplicity, sp);
        if (validEventPlanes) {
          registry.fill(HIST("PP/hEpNumLSmm"), pair.M(), pair.Pt(), multiplicity, ep);
        }
      }
    });
  }
  PROCESS_SWITCH(F0980FlowSP, processPP, "Process pp invariant-mass scalar-product and event-plane analyses", true);

  void processOO(OOCollisions::iterator const& collision, TracksTPCTOF const& tracks,
                 aod::BCsWithTimestamps const&)
  {
    if (collisionSystem(collision) != CollisionSystem::OO) {
      return;
    }
    const float centrality = collision.centFT0C();
    if (!selectOOEventWithQA(collision, centrality)) {
      return;
    }

    QVectorData qVector;
    const auto qVectorStatus = getQVectorData(collision, qVector);
    if (qVectorStatus == QVectorStatus::IndexUnavailable || qVectorStatus == QVectorStatus::Uncalibrated) {
      return;
    }
    registry.fill(HIST("QA/Event/hSelection"), static_cast<int>(EventSelectionStep::QVectorCalibrated));
    if (qVectorStatus != QVectorStatus::Valid) {
      return;
    }
    registry.fill(HIST("QA/Event/hSelection"), static_cast<int>(EventSelectionStep::QVectorValid));

    const auto [qAx, qAy, qBx, qBy, qCx, qCy, qAMagnitude, qBMagnitude, qCMagnitude] = qVector;

    const int runNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    fillQVectorQA(centrality, runNumber, qAx, qAy, qBx, qBy, qCx, qCy);
    registry.fill(HIST("OO/pQAB"), centrality, qAx * qBx + qAy * qBy);
    registry.fill(HIST("OO/pQAC"), centrality, qAx * qCx + qAy * qCy);
    registry.fill(HIST("OO/pQBC"), centrality, qBx * qCx + qBy * qCy);
    registry.fill(HIST("OO/pQABIm"), centrality, qAy * qBx - qAx * qBy);
    registry.fill(HIST("OO/pQACIm"), centrality, qAy * qCx - qAx * qCy);
    registry.fill(HIST("OO/pQBCIm"), centrality, qBy * qCx - qBx * qCy);
    const bool validEventPlanes = qVector.hasValidEventPlanes();
    if (validEventPlanes) {
      registry.fill(HIST("OO/pCosPsiAB"), centrality, (qAx * qBx + qAy * qBy) / (qAMagnitude * qBMagnitude));
      registry.fill(HIST("OO/pCosPsiAC"), centrality, (qAx * qCx + qAy * qCy) / (qAMagnitude * qCMagnitude));
      registry.fill(HIST("OO/pCosPsiBC"), centrality, (qBx * qCx + qBy * qCy) / (qBMagnitude * qCMagnitude));
    }

    loopPairs(tracks, [&](int8_t sign1, int8_t sign2, P4 const& pair) {
      const double ux = std::cos(2. * pair.Phi());
      const double uy = std::sin(2. * pair.Phi());
      const double sp = ux * qAx + uy * qAy;
      const double ep = validEventPlanes ? (ux * qAx + uy * qAy) / qAMagnitude : 0.;
      if (sign1 * sign2 < 0) {
        registry.fill(HIST("OO/hMassUS"), pair.M(), pair.Pt(), centrality);
        registry.fill(HIST("OO/hSpNumUS"), pair.M(), pair.Pt(), centrality, sp);
        if (validEventPlanes) {
          registry.fill(HIST("OO/hEpNumUS"), pair.M(), pair.Pt(), centrality, ep);
        }
      } else if (sign1 > 0 && sign2 > 0) {
        registry.fill(HIST("OO/hMassLSpp"), pair.M(), pair.Pt(), centrality);
        registry.fill(HIST("OO/hSpNumLSpp"), pair.M(), pair.Pt(), centrality, sp);
        if (validEventPlanes) {
          registry.fill(HIST("OO/hEpNumLSpp"), pair.M(), pair.Pt(), centrality, ep);
        }
      } else if (sign1 < 0 && sign2 < 0) {
        registry.fill(HIST("OO/hMassLSmm"), pair.M(), pair.Pt(), centrality);
        registry.fill(HIST("OO/hSpNumLSmm"), pair.M(), pair.Pt(), centrality, sp);
        if (validEventPlanes) {
          registry.fill(HIST("OO/hEpNumLSmm"), pair.M(), pair.Pt(), centrality, ep);
        }
      }
    });
  }
  PROCESS_SWITCH(F0980FlowSP, processOO, "Process OO invariant-mass scalar-product and event-plane analyses", true);

  void processPO(POCollisions::iterator const& collision, TracksTPCTOF const& tracks,
                 aod::BCsWithTimestamps const&)
  {
    if (collisionSystem(collision) != CollisionSystem::PO) {
      return;
    }
    const float centrality = collision.centFT0C();
    if (!selectOOEventWithQA(collision, centrality)) {
      return;
    }

    QVectorData qVector;
    const auto qVectorStatus = getQVectorData(collision, qVector);
    if (qVectorStatus == QVectorStatus::IndexUnavailable || qVectorStatus == QVectorStatus::Uncalibrated) {
      return;
    }
    registry.fill(HIST("QA/Event/hSelection"), static_cast<int>(EventSelectionStep::QVectorCalibrated));
    if (qVectorStatus != QVectorStatus::Valid) {
      return;
    }
    registry.fill(HIST("QA/Event/hSelection"), static_cast<int>(EventSelectionStep::QVectorValid));

    const auto [qAx, qAy, qBx, qBy, qCx, qCy, qAMagnitude, qBMagnitude, qCMagnitude] = qVector;

    const int runNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    fillQVectorQA(centrality, runNumber, qAx, qAy, qBx, qBy, qCx, qCy);
    registry.fill(HIST("PO/pQAB"), centrality, qAx * qBx + qAy * qBy);
    registry.fill(HIST("PO/pQAC"), centrality, qAx * qCx + qAy * qCy);
    registry.fill(HIST("PO/pQBC"), centrality, qBx * qCx + qBy * qCy);
    registry.fill(HIST("PO/pQABIm"), centrality, qAy * qBx - qAx * qBy);
    registry.fill(HIST("PO/pQACIm"), centrality, qAy * qCx - qAx * qCy);
    registry.fill(HIST("PO/pQBCIm"), centrality, qBy * qCx - qBx * qCy);
    const bool validEventPlanes = qVector.hasValidEventPlanes();
    if (validEventPlanes) {
      registry.fill(HIST("PO/pCosPsiAB"), centrality, (qAx * qBx + qAy * qBy) / (qAMagnitude * qBMagnitude));
      registry.fill(HIST("PO/pCosPsiAC"), centrality, (qAx * qCx + qAy * qCy) / (qAMagnitude * qCMagnitude));
      registry.fill(HIST("PO/pCosPsiBC"), centrality, (qBx * qCx + qBy * qCy) / (qBMagnitude * qCMagnitude));
    }

    loopPairs(tracks, [&](int8_t sign1, int8_t sign2, P4 const& pair) {
      const double ux = std::cos(2. * pair.Phi());
      const double uy = std::sin(2. * pair.Phi());
      const double sp = ux * qAx + uy * qAy;
      const double ep = validEventPlanes ? sp / qAMagnitude : 0.;
      if (sign1 * sign2 < 0) {
        registry.fill(HIST("PO/hMassUS"), pair.M(), pair.Pt(), centrality);
        registry.fill(HIST("PO/hSpNumUS"), pair.M(), pair.Pt(), centrality, sp);
        if (validEventPlanes) {
          registry.fill(HIST("PO/hEpNumUS"), pair.M(), pair.Pt(), centrality, ep);
        }
      } else if (sign1 > 0 && sign2 > 0) {
        registry.fill(HIST("PO/hMassLSpp"), pair.M(), pair.Pt(), centrality);
        registry.fill(HIST("PO/hSpNumLSpp"), pair.M(), pair.Pt(), centrality, sp);
        if (validEventPlanes) {
          registry.fill(HIST("PO/hEpNumLSpp"), pair.M(), pair.Pt(), centrality, ep);
        }
      } else if (sign1 < 0 && sign2 < 0) {
        registry.fill(HIST("PO/hMassLSmm"), pair.M(), pair.Pt(), centrality);
        registry.fill(HIST("PO/hSpNumLSmm"), pair.M(), pair.Pt(), centrality, sp);
        if (validEventPlanes) {
          registry.fill(HIST("PO/hEpNumLSmm"), pair.M(), pair.Pt(), centrality, ep);
        }
      }
    });
  }
  PROCESS_SWITCH(F0980FlowSP, processPO, "Process pO invariant-mass scalar-product and event-plane analyses", true);

  void processOOTPCOnly(OOCollisions::iterator const& collision, TracksTPC const& tracks,
                        aod::BCsWithTimestamps const&)
  {
    if (collisionSystem(collision) != CollisionSystem::OO) {
      return;
    }
    const float centrality = collision.centFT0C();
    if (!selectOOEventWithQA(collision, centrality)) {
      return;
    }

    QVectorData qVector;
    const auto qVectorStatus = getQVectorData(collision, qVector);
    if (qVectorStatus == QVectorStatus::IndexUnavailable || qVectorStatus == QVectorStatus::Uncalibrated) {
      return;
    }
    registry.fill(HIST("QA/Event/hSelection"), static_cast<int>(EventSelectionStep::QVectorCalibrated));
    if (qVectorStatus != QVectorStatus::Valid) {
      return;
    }
    registry.fill(HIST("QA/Event/hSelection"), static_cast<int>(EventSelectionStep::QVectorValid));

    const auto [qAx, qAy, qBx, qBy, qCx, qCy, qAMagnitude, qBMagnitude, qCMagnitude] = qVector;

    const int runNumber = collision.bc_as<aod::BCsWithTimestamps>().runNumber();
    fillQVectorQA(centrality, runNumber, qAx, qAy, qBx, qBy, qCx, qCy);
    registry.fill(HIST("OO/pQAB"), centrality, qAx * qBx + qAy * qBy);
    registry.fill(HIST("OO/pQAC"), centrality, qAx * qCx + qAy * qCy);
    registry.fill(HIST("OO/pQBC"), centrality, qBx * qCx + qBy * qCy);
    registry.fill(HIST("OO/pQABIm"), centrality, qAy * qBx - qAx * qBy);
    registry.fill(HIST("OO/pQACIm"), centrality, qAy * qCx - qAx * qCy);
    registry.fill(HIST("OO/pQBCIm"), centrality, qBy * qCx - qBx * qCy);
    const bool validEventPlanes = qVector.hasValidEventPlanes();
    if (validEventPlanes) {
      registry.fill(HIST("OO/pCosPsiAB"), centrality, (qAx * qBx + qAy * qBy) / (qAMagnitude * qBMagnitude));
      registry.fill(HIST("OO/pCosPsiAC"), centrality, (qAx * qCx + qAy * qCy) / (qAMagnitude * qCMagnitude));
      registry.fill(HIST("OO/pCosPsiBC"), centrality, (qBx * qCx + qBy * qCy) / (qBMagnitude * qCMagnitude));
    }

    loopPairsSelected(tracks, [&](auto const& track) { return selectTrack(track); }, [&](auto const& track) { return selectPionTPCOnly(track); }, [&](int8_t sign1, int8_t sign2, P4 const& pair) {
                        const double ux = std::cos(2. * pair.Phi());
                        const double uy = std::sin(2. * pair.Phi());
                        const double sp = ux * qAx + uy * qAy;
                        const double ep = validEventPlanes ? (ux * qAx + uy * qAy) / qAMagnitude : 0.;
                        if (sign1 * sign2 < 0) {
                          registry.fill(HIST("OO/hMassUS"), pair.M(), pair.Pt(), centrality);
                          registry.fill(HIST("OO/hSpNumUS"), pair.M(), pair.Pt(), centrality, sp);
                          if (validEventPlanes) {
                            registry.fill(HIST("OO/hEpNumUS"), pair.M(), pair.Pt(), centrality, ep);
                          }
                        } else if (sign1 > 0 && sign2 > 0) {
                          registry.fill(HIST("OO/hMassLSpp"), pair.M(), pair.Pt(), centrality);
                          registry.fill(HIST("OO/hSpNumLSpp"), pair.M(), pair.Pt(), centrality, sp);
                          if (validEventPlanes) {
                            registry.fill(HIST("OO/hEpNumLSpp"), pair.M(), pair.Pt(), centrality, ep);
                          }
                        } else if (sign1 < 0 && sign2 < 0) {
                          registry.fill(HIST("OO/hMassLSmm"), pair.M(), pair.Pt(), centrality);
                          registry.fill(HIST("OO/hSpNumLSmm"), pair.M(), pair.Pt(), centrality, sp);
                          if (validEventPlanes) {
                            registry.fill(HIST("OO/hEpNumLSmm"), pair.M(), pair.Pt(), centrality, ep);
                          }
                        } });
  }
  PROCESS_SWITCH(F0980FlowSP, processOOTPCOnly, "Process OO with TPC-only pion PID", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<F0980FlowSP>(cfgc)};
}
