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
/// \brief A filter task for strangeness triggers
/// \author Chiara De Martin (chiara.de.martin@cern.ch)
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \since June 1, 2021

#include "../filterTables.h"

#include "PWGLF/DataModel/LFKinkDecayTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/Utils/strangenessBuilderHelper.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "CommonConstants/PhysicsConstants.h"
#include "DCAFitter/DCAFitterN.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/TrackParametrization.h"

#include "TVector3.h"
#include <Math/GenVector/Boost.h>
#include <TLorentzVector.h>

#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

namespace stfilter
{
enum species { Xi = 0,
               Omega = 1 };
constexpr double massSigmaParameters[4][2]{
  {4.9736e-3, 0.006815},
  {-2.39594, -2.257},
  {1.8064e-3, 0.00138},
  {1.03468e-1, 0.1898}};
static const std::vector<std::string> massSigmaParameterNames{"p0", "p1", "p2", "p3"};
static const std::vector<std::string> speciesNames{"Xi", "Omega"};
} // namespace stfilter

float CalculateDCAStraightToPV(float X, float Y, float Z, float Px, float Py, float Pz, float pvX, float pvY, float pvZ)
{
  return std::hypot((pvY - Y) * Pz - (pvZ - Z) * Py, (pvX - X) * Pz - (pvZ - Z) * Px, (pvX - X) * Py - (pvY - Y) * Px) / std::sqrt(Px * Px + Py * Py + Pz * Pz);
}

struct strangenessFilter {

  // Recall the output table
  Produces<aod::StrangenessFilters> strgtable;

  // Define a histograms and registries
  HistogramRegistry QAHistos{"QAHistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry QAHistosTopologicalVariables{"QAHistosTopologicalVariables", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry QAHistosTriggerParticles{"QAHistosTriggerParticles", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry QAHistosStrangenessTracking{"QAHistosStrangenessTracking", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry QAHistosSigma{"QAHistosSigma", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  HistogramRegistry EventsvsMultiplicity{"EventsvsMultiplicity", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  OutputObj<TH1D> hProcessedEvents{TH1D("hProcessedEvents", "Strangeness - event filtered;; Number of events", 21, -1., 20.)};
  OutputObj<TH1F> hCandidate{TH1F("hCandidate", "; Candidate pass selection; Number of events", 30, 0., 30.)};
  OutputObj<TH1F> hEvtvshMinPt{TH1F("hEvtvshMinPt", " Number of h-Omega events with pT_h higher than thrd; min p_{T, trigg} (GeV/c); Number of events", 11, 0., 11.)};

  // Dedicated selection criteria for lambda-lambda
  struct : ConfigurableGroup {
    Configurable<float> cfgv0radiusMin{"cfgv0radiusMin", 1.2, "minimum decay radius"};
    Configurable<float> cfgDCAPosToPVMin{"cfgDCAPosToPVMin", 0.05, "minimum DCA to PV for positive track"};
    Configurable<float> cfgDCANegToPVMin{"cfgDCANegToPVMin", 0.2, "minimum DCA to PV for negative track"};
    Configurable<float> cfgv0CosPA{"cfgv0CosPA", 0.995, "minimum v0 cosine"};
    Configurable<float> cfgDCAV0Dau{"cfgDCAV0Dau", 1.0, "maximum DCA between daughters"};
    Configurable<float> cfgV0PtMin{"cfgV0PtMin", 0, "minimum pT for lambda"};
    Configurable<float> cfgV0RapMin{"cfgV0RapMin", -0.5, "maximum rapidity"};
    Configurable<float> cfgV0RapMax{"cfgV0RapMax", 0.5, "maximum rapidity"};
    Configurable<float> cfgV0LifeTime{"cfgV0LifeTime", 30., "maximum lambda lifetime"};
    Configurable<int16_t> cfgDaughTPCnclsMin{"cfgDaughTPCnclsMin", 70, "minimum fired crossed rows"};
    Configurable<uint8_t> cfgITSNclus{"cfgITSNclus", 1, "minimum its cluster"};
    Configurable<float> cfgRCrossedFindable{"cfgRCrossedFindable", 0.0, "minimum ratio of crossed rows over findable clusters"};
    Configurable<float> cfgDaughPIDCutsTPCPr{"cfgDaughPIDCutsTPCPr", 5, "proton nsigma for TPC"};
    Configurable<float> cfgDaughPIDCutsTPCPi{"cfgDaughPIDCutsTPCPi", 5, "pion nsigma for TPC"};
    Configurable<float> cfgDaughEtaMin{"cfgDaughEtaMin", -0.8, "minimum daughter eta"};
    Configurable<float> cfgDaughEtaMax{"cfgDaughEtaMax", 0.8, "maximum daughter eta"};
    Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.5, "minimum daughter proton pt"};
    Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.5, "minimum daughter pion pt"};
    Configurable<float> cfgLambdaMassWindow{"cfgLambdaMassWindow", 0.01, "window for lambda mass selection"};
    Configurable<float> cfgCompV0Rej{"cfgCompV0Rej", 0.01, "competing V0 rejection"};
    Configurable<float> cfgMinCPAV0V0{"cfgMinCPAV0V0", 0.8, "minimum CPA of v0v0"};
    Configurable<float> cfgMaxRadiusV0V0{"cfgMaxRadiusV0V0", 10.0, "maximum radius of v0v0"};
    Configurable<float> cfgMaxDistanceV0V0{"cfgMaxDistanceV0V0", 5.0, "maximum distance of v0v0"};
    Configurable<float> cfgMaxDCAV0V0{"cfgMaxDCAV0V0", 5.0, "maximum DCA of v0v0"};
  } cfgLLCuts;

  // Selection criteria for cascades
  Configurable<bool> useCascadeMomentumAtPrimVtx{"useCascadeMomentumAtPrimVtx", false, "use cascade momentum at PV"};
  Configurable<bool> doextraQA{"doextraQA", 1, "do extra QA"};
  Configurable<float> cutzvertex{"cutzvertex", 100.0f, "Accepted z-vertex range"};
  Configurable<int> tpcmincrossedrows{"tpcmincrossedrows", 50, "Min number of crossed TPC rows"};
  Configurable<float> v0cospa{"v0cospa", 0.95, "V0 CosPA"};
  Configurable<float> casccospaxi{"casccospaxi", 0.95, "Casc CosPA"};
  Configurable<float> casccospaomega{"casccospaomega", 0.95, "Casc CosPA"};
  Configurable<float> dcav0dau{"dcav0dau", 2.0, "DCA V0 Daughters"};
  Configurable<float> dcacascdau{"dcacascdau", 2.0, "DCA Casc Daughters"};
  Configurable<float> dcamesontopv{"dcamesontopv", 0.05, "DCA Meson To PV"};
  Configurable<float> dcabaryontopv{"dcabaryontopv", 0.05, "DCA Baryon To PV"};
  Configurable<float> dcabachtopv{"dcabachtopv", 0.05, "DCA Bach To PV"};
  Configurable<float> dcav0topv{"dcav0topv", 0.0, "DCA V0 To PV"};
  Configurable<float> v0radius{"v0radius", 1.0, "V0 Radius"};
  Configurable<float> cascradius{"cascradius", 0.6, "cascradius"};
  Configurable<float> rapidity{"rapidity", 2, "rapidity"};
  Configurable<float> eta{"eta", 2, "Eta"};
  Configurable<float> minpt{"minpt", 0.5, "minpt"};
  Configurable<float> etadau{"etadau", 0.9, "EtaDaughters"};
  Configurable<float> masslambdalimit{"masslambdalimit", 0.02, "masslambdalimit"};
  Configurable<float> omegarej{"omegarej", 0.005, "omegarej"};
  Configurable<float> xirej{"xirej", 0.008, "xirej"}; // merge the two rejection variables into one?
  Configurable<float> ximasswindow{"ximasswindow", 0.075, "Xi Mass Window"};
  Configurable<float> omegamasswindow{"omegamasswindow", 0.075, "Omega Mass Window"}; // merge the two windows variables into one?
  Configurable<int> properlifetimefactor{"properlifetimefactor", 5, "Proper Lifetime cut"};
  Configurable<float> lowerradiusXiYN{"lowerradiusXiYN", 24.39, "Cascade lower radius for single Xi trigger"};
  Configurable<float> lowerradiusOmega{"lowerradiusOmega", 19.0, "Omega lower radius for high radius Omega trigger"};
  Configurable<float> upperradiusOmega{"upperradiusOmega", 19.0, "Omega upper radius for low radius Omega trigger"};
  Configurable<float> nsigmatpcpi{"nsigmatpcpi", 6, "N Sigmas TPC pi"};
  Configurable<float> nsigmatpcka{"nsigmatpcka", 6, "N Sigmas TPC ka"};
  Configurable<float> nsigmatpcpr{"nsigmatpcpr", 6, "N Sigmas TPC pr"};
  Configurable<bool> hastof{"hastof", 1, "Has TOF (OOB condition)"};
  Configurable<float> ptthrtof{"ptthrtof", 1.0, "Pt threshold to apply TOF condition"};
  Configurable<bool> sel8{"sel8", 0, "Apply sel8 event selection"};
  Configurable<bool> isTriggerTVX{"isTriggerTVX", 1, "Require TVX"};
  struct : ConfigurableGroup {
    Configurable<int> HMTrgSelectionForOmegaFT0M{"HMTrgSelectionForOmegaFT0M", 1, "0: none, 1: normalised FT0M, 2: FT0M "};
    Configurable<int> LowLimitHMTrgOmegaT0M{"LowLimitHMTrgOmegaT0M", 3100, "T0M"};
    Configurable<int> LowLimitHMTrgOmegaT0MNorm{"LowLimitHMTrgOmegaT0MNorm", 70, "normalised T0M selection [2] of multFiler"};
    Configurable<int> LowLimitHMTrgT0MNorm{"LowLimitHMTrgT0MNorm", 140, "normalised T0M selection [2] of multFiler"};
    Configurable<int> HMTrgSelectionForOmegaTrks{"HMTrgSelectionForOmegaTrks", 1, "0: none, 1: GlobalMult,2: selectTrack"};
    Configurable<int> LowLimitHMTrgOmegaTrkGlob{"LowLimitHMTrgOmegaTrkGlob", 45, "Omega HM GlobalMult"};
    Configurable<int> LowLimitHMTrgOmegaTrkSel{"LowLimitHMTrgOmegaTrkSel", 50, "Omega HM selectTrackHMO"};
    Configurable<int> LowLimitHMTrgTrkGlob{"LowLimitHMTrgTrksGlob", 100, "HM Omega normalisation GlobalMult"};
    Configurable<int> LowLimitHMTrgTrkSel{"LowLimitHMTrgTrkSel", 50, "HM Omega normalisation selectTrackHMO"};
    Configurable<float> hEtaHM{"hEtaHM", 1.0f, "Eta range for particles defining HM events"};
    Configurable<float> hMinPtHM{"hMinPtHM", 0.2f, "Min pt for particles defining HM events"};
  } cfgHMOmegaCuts;
  Configurable<float> avPyT0C{"avPyT0C", 8.83, "nch from pythia T0C"};
  Configurable<float> avPyT0A{"avPyT0A", 8.16, "nch from pythia T0A"};
  Configurable<bool> isTimeFrameBorderCut{"isTimeFrameBorderCut", 1, "Apply timeframe border cut"};
  Configurable<bool> useSigmaBasedMassCutXi{"useSigmaBasedMassCutXi", true, "Mass window based on n*sigma instead of fixed"};
  Configurable<bool> useSigmaBasedMassCutOmega{"useSigmaBasedMassCutOmega", true, "Mass window based on n*sigma instead of fixed"};
  Configurable<float> massWindowOmegaNsigma{"massWindowOmegaNsigma", 6, "Inv. mass window for tracked Omega"};
  Configurable<float> massWindowXiNsigma{"massWindowXiNsigma", 6, "Inv. mass window for tracked Xi"};

  // Selections criteria for tracks
  struct : ConfigurableGroup {
    Configurable<float> hEta{"hEta", 0.9f, "Eta range for trigger particles"};
    Configurable<float> hMinPt{"hMinPt", 1.0f, "Min pt for trigger particles"};
    Configurable<bool> isTrackFilter{"isTrackFilter", true, "Apply track myTrackSelections"};
  } cfgTrackCuts;

  // Settings for strangeness tracking filter
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<bool> propToDCA{"propToDCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};
  Configurable<int> minNoClsTrackedCascade{"minNoClsTrackedCascade", 70, "Minimum number of clusters required for daughters of tracked cascades"};
  Configurable<float> minPtTrackedCascade{"minPtTrackedCascade", 0., "Min. pt for tracked cascades"};
  Configurable<bool> useNsigmaCutTrackedXi{"useNsigmaCutTrackedXi", true, "Mass window based on n*sigma instead of fixed"};
  Configurable<bool> useNsigmaCutTrackedOmega{"useNsigmaCutTrackedOmega", true, "Mass window based on n*sigma instead of fixed"};
  Configurable<float> massWindowTrackedOmegaNsigma{"massWindowTrackedOmegaNsigma", 6, "Inv. mass window for tracked Omega"};
  Configurable<float> massWindowTrackedXiNsigma{"massWindowTrackedXiNsigma", 6, "Inv. mass window for tracked Xi"};
  Configurable<float> massWindowTrackedOmega{"massWindowTrackedOmega", 0.05, "Inv. mass window for tracked Omega"};
  Configurable<float> massWindowXiExclTrackedOmega{"massWindowXiExclTrackedOmega", 0.005, "Inv. mass window for exclusion of Xi for tracked Omega-"};
  Configurable<float> massWindowTrackedXi{"massWindowTrackedXi", 0.05, "Inv. mass window for tracked Xi"};
  Configurable<float> massWindowLambda{"massWindowLambda", 0.05, "Inv. mass window for Lambda (ST)"};
  Configurable<float> maxMatchingChi2TrackedCascade{"maxMatchingChi2TrackedCascade", 2000., "Max matching chi2 for tracked cascades"};
  Configurable<bool> recalculateMasses{"recalculateMasses", true, "Recalculate Xi/Omega masses"};
  Configurable<float> maxNSigmaBachelorTrackedXi{"maxNSigmaBachelorTrackedXi", 4., "Max Nsigma for bachelor of tracked Xi (pi)"};
  Configurable<float> maxNSigmaBachelorTrackedOmega{"maxNSigmaBachelorTrackedOmega", 4., "Max Nsigma for bachelor of tracked Xi (Ka)"};
  Configurable<float> maxNSigmaV0PrTrackedCascade{"maxNSigmaV0PrTrackedCascade", 4., "Max Nsigma for proton from V0 fromtracked Xi"};
  Configurable<float> maxNSigmaV0PiTrackedCascade{"maxNSigmaV0PiTrackedCascade", 4., "Max Nsigma for pion from V0 fromtracked Xi"};
  Configurable<float> minDcaTrackedXi{"minDcaTrackedXi", -1., "Minimum DCA for tracked cascades"};
  Configurable<float> maxCpaTrackedXi{"maxCpaTrackedXi", 1., "Maximum CPA for tracked cascades"};
  Configurable<float> minDcaTrackedOmega{"minDcaTrackedOmega", -1., "Minimum DCA for tracked cascades (ST)"};
  Configurable<float> maxCpaTrackedOmega{"maxCpaTrackedOmega", 1., "Maximum CPA for tracked cascades (ST)"};

  // Settings for sigmaplus filter
  struct : ConfigurableGroup {
    Configurable<float> nsigmatpcSigma{"nsigmatpcSigma", 3.f, "N Sigmas TPC Sigma"};
    Configurable<float> nsigmatofSigma{"nsigmatofSigma", 3.f, "N Sigmas TOF Sigma"};
    Configurable<float> minMassSigma{"minMassSigma", 1.15f, "min mass for Sigma"};
    Configurable<float> maxMassSigma{"maxMassSigma", 1.25f, "max mass for Sigma"};
    Configurable<float> minPtSigma{"minPtSigma", 1.2f, "min Pt for Sigma"};
    Configurable<float> minQtAPSigma{"minQtAPSigma", 0.15f, "min QtAP for Sigma"};
    Configurable<float> maxQtAPSigma{"maxQtAPSigma", 0.20f, "max QtAP for Sigma"};
    Configurable<float> maxDCAtoPVSigma{"maxDCAtoPVSigma", 0.1f, "Max DCA to primary vertex for Sigma candidates (cm)"};
    Configurable<float> minRadiusSigma{"minRadiusSigma", 19.f, "Min radius for Sigma+ decay vertex (cm)"};
    Configurable<float> minCosPASigma{"minCosPASigma", 0.995f, "Min Cosine of pointing angle for Sigma candidates"};
    Configurable<float> minPtProtonTOF{"minPtProtonTOF", 0.75f, "Min Pt for proton to have TOF signal (GeV/c)"};
    Configurable<float> maxKStarSigmaProton{"maxKStarSigmaProton", 0.8f, "Max k* for Sigma-Proton pairs (GeV/c)"};
  } cfgSigma;

  Configurable<LabeledArray<double>> parSigmaMass{
    "parSigmaMass",
    {stfilter::massSigmaParameters[0], 4, 2,
     stfilter::massSigmaParameterNames, stfilter::speciesNames},
    "Mass resolution parameters: [0]*exp([1]*x)+[2]*exp([3]*x)"};

  // helper object
  o2::pwglf::strangenessBuilderHelper mStraHelper;
  o2::vertexing::DCAFitterN<2> mDCAFitter;

  /// CCDB and info/objects to be fetched from it
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = 0;
  float mBz = 0.;
  std::vector<double>* mMeanMultT0C;
  std::vector<double>* mMeanMultT0A;

  bool selectTrack(const auto& track)
  {
    return track.pt() > cfgTrackCuts.hMinPt && std::abs(track.eta()) < cfgTrackCuts.hEta && track.tpcNClsCrossedRows() >= tpcmincrossedrows && track.tpcChi2NCl() <= 4.f && track.itsChi2NCl() <= 36.f && (track.itsClusterMap() & 0x7) != 0;
  }
  bool selectTrackOHM(const auto& track)
  {
    return track.pt() > cfgHMOmegaCuts.hMinPtHM && std::abs(track.eta()) < cfgHMOmegaCuts.hEtaHM && track.tpcNClsCrossedRows() >= tpcmincrossedrows && track.tpcChi2NCl() <= 4.f && track.itsChi2NCl() <= 36.f && (track.itsClusterMap() & 0x7) != 0;
  }
  float getV0V0DCA(TVector3 v01pos, TVector3 v01mom, TVector3 v02pos, TVector3 v02mom)
  {
    TVector3 posdiff = v02pos - v01pos;
    TVector3 cross = v01mom.Cross(v02mom);
    TVector3 dcaVec = (posdiff.Dot(cross) / cross.Mag2()) * cross;
    return dcaVec.Mag();
  }
  float getV0V0CPA(TVector3 v01mom, TVector3 v02mom)
  {
    return v01mom.Dot(v02mom) / (v01mom.Mag() * v02mom.Mag());
  }
  float getV0V0Distance(TVector3 v01pos, TVector3 v02pos)
  {
    TVector3 posdiff = v02pos - v01pos;
    return posdiff.Mag();
  }
  float getV0V0Radius(TVector3 v01pos, TVector3 v01mom, TVector3 v02pos, TVector3 v02mom)
  {
    TVector3 posdiff = v02pos - v01pos;
    v01mom *= 1. / v01mom.Mag();
    v02mom *= 1. / v02mom.Mag();
    float dd = 1. - TMath::Power(v01mom.Dot(v02mom), 2);
    if (dd < 1e-5)
      return 999;
    float tt = posdiff.Dot(v01mom - v01mom.Dot(v02mom) * v02mom) / dd;
    float ss = -posdiff.Dot(v02mom - v01mom.Dot(v02mom) * v01mom) / dd;
    TVector3 radVec = v01pos + v02pos + tt * v01mom + ss * v02mom;
    radVec *= 0.5;
    return radVec.Mag();
  }
  bool isSelectedV0V0(TVector3 v01pos, TVector3 v01mom, TVector3 v02pos, TVector3 v02mom)
  {
    if (getV0V0DCA(v01pos, v01mom, v02pos, v02mom) > cfgLLCuts.cfgMaxDCAV0V0)
      return false;
    if (getV0V0CPA(v01mom, v02mom) < cfgLLCuts.cfgMinCPAV0V0)
      return false;
    if (getV0V0Distance(v01pos, v02pos) > cfgLLCuts.cfgMaxDistanceV0V0)
      return false;
    if (getV0V0Radius(v01pos, v01mom, v02pos, v02mom) > cfgLLCuts.cfgMaxRadiusV0V0)
      return false;

    return true;
  }

  template <typename T>
  bool isSelectedV0Daughter(T const& track)
  {
    if (track.tpcNClsCrossedRows() < cfgLLCuts.cfgDaughTPCnclsMin)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < cfgLLCuts.cfgRCrossedFindable)
      return false;
    if (track.itsNCls() < cfgLLCuts.cfgITSNclus)
      return false;
    if (track.eta() > cfgLLCuts.cfgDaughEtaMax)
      return false;
    if (track.eta() < cfgLLCuts.cfgDaughEtaMin)
      return false;

    return true;
  }
  template <typename T>
  bool isSelectedV0DaughterPID(T const& track, int pid) // pid 0: proton, pid 1: pion
  {
    if (pid == 0 && std::abs(track.tpcNSigmaPr()) > cfgLLCuts.cfgDaughPIDCutsTPCPr)
      return false;
    if (pid == 1 && std::abs(track.tpcNSigmaPi()) > cfgLLCuts.cfgDaughPIDCutsTPCPi)
      return false;
    if (pid == 0 && track.pt() < cfgLLCuts.cfgDaughPrPt)
      return false;
    if (pid == 1 && track.pt() < cfgLLCuts.cfgDaughPiPt)
      return false;

    return true;
  }

  float getAlphaAP(const std::array<float, 3>& momMother, const std::array<float, 3>& momKink)
  {
    std::array<float, 3> momMissing = {momMother[0] - momKink[0], momMother[1] - momKink[1], momMother[2] - momKink[2]};
    float lQlP = std::inner_product(momMother.begin(), momMother.end(), momKink.begin(), 0.f);
    float lQlN = std::inner_product(momMother.begin(), momMother.end(), momMissing.begin(), 0.f);
    return (lQlP - lQlN) / (lQlP + lQlN);
  }

  float getQtAP(const std::array<float, 3>& momMother, const std::array<float, 3>& momKink)
  {
    float dp = std::inner_product(momMother.begin(), momMother.end(), momKink.begin(), 0.f);
    float p2V0 = std::inner_product(momMother.begin(), momMother.end(), momMother.begin(), 0.f);
    float p2A = std::inner_product(momKink.begin(), momKink.end(), momKink.begin(), 0.f);
    return std::sqrt(p2A - dp * dp / p2V0);
  }

  float getCosPA(const std::array<float, 3>& momMother, const std::array<float, 3>& decayVertex, const std::array<float, 3>& primaryVertex)
  {
    std::array<float, 3> decayVec = {decayVertex[0] - primaryVertex[0], decayVertex[1] - primaryVertex[1], decayVertex[2] - primaryVertex[2]};
    float dotProduct = std::inner_product(momMother.begin(), momMother.end(), decayVec.begin(), 0.f);
    float momMotherMag = std::sqrt(std::inner_product(momMother.begin(), momMother.end(), momMother.begin(), 0.f));
    float decayVecMag = std::sqrt(std::inner_product(decayVec.begin(), decayVec.end(), decayVec.begin(), 0.f));
    return dotProduct / (momMotherMag * decayVecMag);
  }

  float getKStar(TLorentzVector const& part1, TLorentzVector const& part2)
  {
    TLorentzVector trackSum, PartOneCMS, PartTwoCMS, trackRelK;
    trackSum = part1 + part2;
    const float beta = trackSum.Beta();
    const float betax = beta * std::cos(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betay = beta * std::sin(trackSum.Phi()) * std::sin(trackSum.Theta());
    const float betaz = beta * std::cos(trackSum.Theta());
    PartOneCMS.SetXYZM(part1.Px(), part1.Py(), part1.Pz(), part1.M());
    PartTwoCMS.SetXYZM(part2.Px(), part2.Py(), part2.Pz(), part2.M());
    const ROOT::Math::Boost boostPRF = ROOT::Math::Boost(-betax, -betay, -betaz);
    PartOneCMS = boostPRF(PartOneCMS);
    PartTwoCMS = boostPRF(PartTwoCMS);
    trackRelK = PartOneCMS - PartTwoCMS;
    return 0.5 * trackRelK.P();
  }

  void init(o2::framework::InitContext&)
  {
    // set V0 parameters in the helper
    mStraHelper.v0selections.minCrossedRows = tpcmincrossedrows;
    if (dcamesontopv <= dcabaryontopv)
      mStraHelper.v0selections.dcanegtopv = dcamesontopv;
    else
      mStraHelper.v0selections.dcanegtopv = dcabaryontopv; // get the minimum one
    if (dcamesontopv <= dcabaryontopv)
      mStraHelper.v0selections.dcapostopv = dcamesontopv;
    else
      mStraHelper.v0selections.dcapostopv = dcabaryontopv; // get the minimum one
    mStraHelper.v0selections.v0cospa = v0cospa;
    mStraHelper.v0selections.dcav0dau = dcav0dau;
    mStraHelper.v0selections.v0radius = v0radius;
    mStraHelper.v0selections.maxDaughterEta = etadau;

    // set cascade parameters in the helper
    mStraHelper.cascadeselections.minCrossedRows = tpcmincrossedrows;
    mStraHelper.cascadeselections.dcabachtopv = dcabachtopv;
    mStraHelper.cascadeselections.cascradius = cascradius;
    if (casccospaxi <= casccospaomega)
      mStraHelper.cascadeselections.casccospa = casccospaxi;
    else
      mStraHelper.cascadeselections.casccospa = casccospaomega; // get the minimum one
    mStraHelper.cascadeselections.dcacascdau = dcacascdau;
    mStraHelper.cascadeselections.lambdaMassWindow = masslambdalimit;
    mStraHelper.cascadeselections.maxDaughterEta = etadau;

    hProcessedEvents->GetXaxis()->SetBinLabel(1, "Events processed");
    hProcessedEvents->GetXaxis()->SetBinLabel(2, "Event selection");
    hProcessedEvents->GetXaxis()->SetBinLabel(3, "Events w/ high-#it{p}_{T} hadron");
    hProcessedEvents->GetXaxis()->SetBinLabel(4, aod::filtering::Omega::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(5, aod::filtering::hadronOmega::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(6, aod::filtering::DoubleXi::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(7, aod::filtering::TripleXi::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(8, aod::filtering::QuadrupleXi::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(9, aod::filtering::SingleXiYN::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(10, aod::filtering::OmegaLargeRadius::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(11, "#Xi");
    hProcessedEvents->GetXaxis()->SetBinLabel(12, aod::filtering::TrackedXi::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(13, aod::filtering::TrackedOmega::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(14, aod::filtering::OmegaHighMult::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(15, aod::filtering::DoubleOmega::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(16, aod::filtering::OmegaXi::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(17, "LL");
    hProcessedEvents->GetXaxis()->SetBinLabel(18, aod::filtering::OmegaHighMultTrk::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(19, aod::filtering::HighMultFT0M::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(20, aod::filtering::HighMultTrk::columnLabel());
    hProcessedEvents->GetXaxis()->SetBinLabel(21, aod::filtering::SigmaProton::columnLabel());

    hCandidate->GetXaxis()->SetBinLabel(1, "All");
    hCandidate->GetXaxis()->SetBinLabel(2, "PassBuilderSel");
    hCandidate->GetXaxis()->SetBinLabel(3, "DCA_meson");
    hCandidate->GetXaxis()->SetBinLabel(4, "DCA_baryon");
    hCandidate->GetXaxis()->SetBinLabel(5, "TPCNsigma_pion");
    hCandidate->GetXaxis()->SetBinLabel(6, "TPCNsigma_proton");
    hCandidate->GetXaxis()->SetBinLabel(7, "Eta_dau");
    hCandidate->GetXaxis()->SetBinLabel(8, "DCABachToPV");
    hCandidate->GetXaxis()->SetBinLabel(9, "V0Radius");
    hCandidate->GetXaxis()->SetBinLabel(10, "CascRadius");
    hCandidate->GetXaxis()->SetBinLabel(11, "V0CosPA");
    hCandidate->GetXaxis()->SetBinLabel(12, "DCAV0Dau");
    hCandidate->GetXaxis()->SetBinLabel(13, "DCACascDau");
    hCandidate->GetXaxis()->SetBinLabel(14, "MassLambdaLimit");
    hCandidate->GetXaxis()->SetBinLabel(15, "Eta");
    hCandidate->GetXaxis()->SetBinLabel(16, "HasTOFOneLeg");
    hCandidate->GetXaxis()->SetBinLabel(17, "CascCosPA");
    hCandidate->GetXaxis()->SetBinLabel(18, "DCAV0ToPV");
    hCandidate->GetXaxis()->SetBinLabel(19, "ProperLifeTime");
    hCandidate->GetXaxis()->SetBinLabel(20, "Rapidity");

    std::vector<double> centBinning = {0., 1., 5., 10., 20., 30., 40., 50., 70., 100.};
    AxisSpec multAxisNTPV = {100, 0.0f, 100.0f, "N. tracks PV estimator"};
    AxisSpec multAxisT0M = {600, 0.0f, 6000.0f, "T0M multiplicity estimator"};
    AxisSpec multAxisT0MNorm = {200, 0.0f, 200.0f, "Normalised T0M multiplicity estimator"};
    AxisSpec multAxisTrack = {150, 0.0f, 150.0f, "Track multiplicity"};
    AxisSpec multAxisV0A = {500, 0.0f, 25000.0f, "V0A multiplicity estimator"};
    AxisSpec ximassAxis = {200, 1.28f, 1.36f};
    AxisSpec omegamassAxis = {200, 1.59f, 1.75f};
    AxisSpec ptAxis = {100, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec pTPCAxis = {100, 0.0f, 10.0f, "#it{p} TPC (GeV/#it{c})"};
    AxisSpec etaAxis = {200, -2.0f, 2.0f, "#eta"};
    AxisSpec phiAxis = {100, -TMath::Pi() / 2, 3. * TMath::Pi() / 2, "#varphi"};
    AxisSpec ptTriggAxis = {150, 0.0f, 15.0f, "#it{p}_{T} (GeV/#it{c})"};

    // general QA histograms
    QAHistos.add("hVtxZ", "Z-Vertex distribution after selection;Z (cm)", HistType::kTH1F, {{100, -50, 50}});
    QAHistos.add("hMassXiBefSelvsPt", "hMassXiBefSelvsPt", HistType::kTH2F, {ximassAxis, ptAxis});
    QAHistos.add("hMassOmegaBefSelvsPt", "hMassOmegaBefSelvsPt", HistType::kTH2F, {omegamassAxis, ptAxis});
    QAHistos.add("hMassXiAfterSelvsPt", "hMassXiAfterSelvsPt", HistType::kTH2F, {ximassAxis, ptAxis});
    QAHistos.add("hMassOmegaAfterSelvsPt", "hMassOmegaAfterSelvsPt", HistType::kTH2F, {omegamassAxis, ptAxis});
    QAHistos.add("hPtXi", "pt distribution of selected Xi candidates", HistType::kTH1F, {ptAxis});
    QAHistos.add("hPtOmega", "pt distribution of selected Omega candidates", HistType::kTH1F, {ptAxis});
    QAHistos.add("hEtaXi", "eta distribution of selected Xi candidates", HistType::kTH1F, {etaAxis});
    QAHistos.add("hEtaOmega", "eta distribution of selected Omega candidates", HistType::kTH1F, {etaAxis});

    // topological variables distributions
    QAHistosTopologicalVariables.add("hCascCosPAXi", "hCascCosPAXi", HistType::kTH1F, {{350, 0.65f, 1.0f}});
    QAHistosTopologicalVariables.add("hV0CosPAXi", "hV0CosPAXi", HistType::kTH1F, {{250, 0.75f, 1.0f}});
    QAHistosTopologicalVariables.add("hCascRadiusXi", "hCascRadiusXi", HistType::kTH1F, {{500, 0.0f, 50.0f}});
    QAHistosTopologicalVariables.add("hV0RadiusXi", "hV0RadiusXi", HistType::kTH1F, {{500, 0.0f, 50.0f}});
    QAHistosTopologicalVariables.add("hDCAV0DaughtersXi", "hDCAV0DaughtersXi", HistType::kTH1F, {{110, 0.0f, 2.2f}});
    QAHistosTopologicalVariables.add("hDCACascDaughtersXi", "hDCACascDaughtersXi", HistType::kTH1F, {{110, 0.0f, 2.2f}});
    QAHistosTopologicalVariables.add("hDCAV0ToPVXi", "hDCAV0ToPVXi", HistType::kTH1F, {{220, 0.0f, 2.2f}});
    QAHistosTopologicalVariables.add("hDCABachToPVXi", "|hDCABachToPVXi|", HistType::kTH1F, {{400, 0.0f, 2.0f}});
    QAHistosTopologicalVariables.add("hDCAPosToPVXi", "|hDCAPosToPVXi|", HistType::kTH1F, {{400, 0.0f, 2.0f}});
    QAHistosTopologicalVariables.add("hDCANegToPVXi", "|hDCANegToPVXi|", HistType::kTH1F, {{400, 0.0f, 2.0f}});
    QAHistosTopologicalVariables.add("hInvMassLambdaXi", "InvMassLambdaXi", HistType::kTH1F, {{200, 1.07f, 1.17f}});
    QAHistosTopologicalVariables.add("hProperLifetimeXi", "Proper Lifetime Xi", HistType::kTH1F, {{50, 0, 50}});
    //
    QAHistosTopologicalVariables.add("hCascCosPAOmega", "hCascCosPAOmega", HistType::kTH1F, {{350, 0.65f, 1.0f}});
    QAHistosTopologicalVariables.add("hV0CosPAOmega", "hV0CosPAOmega", HistType::kTH1F, {{250, 0.75f, 1.0f}});
    QAHistosTopologicalVariables.add("hCascRadiusOmega", "hCascRadiusOmega", HistType::kTH1F, {{500, 0.0f, 50.0f}});
    QAHistosTopologicalVariables.add("hV0RadiusOmega", "hV0RadiusOmega", HistType::kTH1F, {{500, 0.0f, 50.0f}});
    QAHistosTopologicalVariables.add("hDCAV0DaughtersOmega", "hDCAV0DaughtersOmega", HistType::kTH1F, {{110, 0.0f, 2.2f}});
    QAHistosTopologicalVariables.add("hDCACascDaughtersOmega", "hDCACascDaughtersOmega", HistType::kTH1F, {{110, 0.0f, 2.2f}});
    QAHistosTopologicalVariables.add("hDCAV0ToPVOmega", "hDCAV0ToPVOmega", HistType::kTH1F, {{220, 0.0f, 2.2f}});
    QAHistosTopologicalVariables.add("hDCABachToPVOmega", "hDCABachToPVOmega", HistType::kTH1F, {{400, 0.0f, 2.0f}});
    QAHistosTopologicalVariables.add("hDCAPosToPVOmega", "hDCAPosToPVOmega", HistType::kTH1F, {{400, 0.0f, 2.0f}});
    QAHistosTopologicalVariables.add("hDCANegToPVOmega", "hDCANegToPVOmega", HistType::kTH1F, {{400, 0.0f, 2.0f}});
    QAHistosTopologicalVariables.add("hInvMassLambdaOmega", "InvMassLambdaOmega", HistType::kTH1F, {{200, 1.07f, 1.17f}});
    QAHistosTopologicalVariables.add("hProperLifetimeOmega", "Proper Lifetime Omega", HistType::kTH1F, {{50, 0, 50}});
    QAHistosTopologicalVariables.add("hCascRadiusOmegaLargeR", "hCascRadiusOmegaLargeR", HistType::kTH1F, {{500, 0.0f, 50.0f}});
    QAHistosTopologicalVariables.add("hCascRadiusXiYN", "hCascRadiusXiYN", HistType::kTH1F, {{500, 0.0f, 50.0f}});

    // trigger particles QA
    QAHistosTriggerParticles.add("hTriggeredParticlesAllEv", "Distribution of #tracks w/ pt > pt,trigg,min", HistType::kTH1F, {{20, 0.5, 20.5, "Trigger counter"}});
    QAHistosTriggerParticles.add("hTriggeredParticlesSelEv", "Distribution of #tracks w/ pt > pt,trigg,min (after sel)", HistType::kTH1F, {{20, 0.5, 20.5, "Trigger counter"}});
    QAHistosTriggerParticles.add("hPtTriggerAllEv", "hPtTriggerAllEv", HistType::kTH1F, {{300, 0, 30, "Pt of trigger particles"}});
    QAHistosTriggerParticles.add("hPtTriggerSelEv", "hPtTriggerSelEv", HistType::kTH1F, {{300, 0, 30, "Pt of trigger particles after selections"}});
    QAHistosTriggerParticles.add("hEtaTriggerAllEv", "hEtaTriggerAllEv", HistType::kTH2F, {{180, -1.4, 1.4, "Eta of trigger particles"}, {ptTriggAxis}});
    QAHistosTriggerParticles.add("hPhiTriggerAllEv", "hPhiTriggerAllEv", HistType::kTH2F, {{100, 0, 2 * TMath::Pi(), "Phi of trigger particles"}, {ptTriggAxis}});

    EventsvsMultiplicity.add("AllEventsvsMultiplicityFT0M", "T0M distribution of all events", HistType::kTH1F, {multAxisT0M});
    EventsvsMultiplicity.add("AllEventsvsMultiplicityFT0MwOmega", "T0M distribution of events w/ Omega candidate", HistType::kTH1F, {multAxisT0M});
    EventsvsMultiplicity.add("AllEventsvsMultiplicityFT0MNorm", "T0M Normalised of all events", HistType::kTH1F, {multAxisT0MNorm});
    EventsvsMultiplicity.add("AllEventsvsMultiplicityFT0MwOmegaNorm", "T0M distribution of events w/ Omega candidate - Normalised FT0M", HistType::kTH1F, {multAxisT0MNorm});
    EventsvsMultiplicity.add("AllEventsvsMultiplicityFT0MNoFT0", "T0M distribution of events without FT0", HistType::kTH1F, {multAxisT0M});
    EventsvsMultiplicity.add("AllEventsvsMultiplicityTracks", "Track distribution of all events", HistType::kTH1F, {multAxisTrack});
    EventsvsMultiplicity.add("AllEventsvsMultiplicityTrackswOmega", "Track distribution of events w/ Omega candidate", HistType::kTH1F, {multAxisTrack});
    EventsvsMultiplicity.add("AllEventsvsMultiplicityTracksGlob", "MultGlob Track distribution of all events", HistType::kTH1F, {multAxisTrack});
    EventsvsMultiplicity.add("AllEventsvsMultiplicityTracksGlobwOmega", "MultGlob Track distribution of events w/ Omega candidate", HistType::kTH1F, {multAxisTrack});
    EventsvsMultiplicity.add("AllEventsvsMultiplicityFT0MTrackswOmega", "Track distribution of events w/ Omega candidate", HistType::kTH1F, {multAxisTrack});
    EventsvsMultiplicity.add("AllEventsvsMultiplicityFT0MTracksGlobwOmega", "MultGlob Track distribution of events w/ Omega candidate", HistType::kTH1F, {multAxisTrack});
    EventsvsMultiplicity.add("AllEventsvsMultiplicityFT0MTrackswOmega2D", "2D Track vs FT0M normalised distribution of events w/ Omega candidate", HistType::kTH2F, {multAxisTrack, multAxisT0MNorm});
    EventsvsMultiplicity.add("AllEventsvsMultiplicityFT0MTracksGlobwOmega2D", "2D Track vs FT0M normalised distribution of events w/ Omega candidate", HistType::kTH2F, {multAxisTrack, multAxisT0MNorm});
    if (doextraQA) {
      EventsvsMultiplicity.add("AllEventsvsMultiplicityZeqV0A", "ZeqV0A distribution of all events", HistType::kTH1F, {multAxisV0A});
      EventsvsMultiplicity.add("hadEventsvsMultiplicityZeqV0A", "ZeqV0A distribution of events with hight pT hadron", HistType::kTH1F, {multAxisV0A});
      EventsvsMultiplicity.add("hadEventsvsMultiplicityZeqV0AvsPt", "ZeqV0A distribution of events with hight pT hadron", HistType::kTH2F, {{multAxisV0A}, {11, 0, 11}});

      EventsvsMultiplicity.add("AllEventsvsMultiplicityZeqT0M", "ZeqT0M distribution of all events", HistType::kTH1F, {multAxisT0M});
      EventsvsMultiplicity.add("hadEventsvsMultiplicityZeqT0M", "ZeqT0M distribution of events with hight pT hadron", HistType::kTH1F, {multAxisT0M});
      EventsvsMultiplicity.add("hadEventsvsMultiplicityZeqT0MvsPt", "ZeqT0M distribution of events with hight pT hadron", HistType::kTH2F, {{multAxisT0M}, {11, 0, 11}});

      EventsvsMultiplicity.add("AllEventsvsMultiplicityZeqNTracksPV", "ZeqNTracksPV distribution of all events", HistType::kTH1F, {multAxisNTPV});
      EventsvsMultiplicity.add("hadEventsvsMultiplicityZeqNTracksPV", "ZeqNTracksPV distribution of events with hight pT hadron", HistType::kTH1F, {multAxisNTPV});
      EventsvsMultiplicity.add("hadEventsvsMultiplicityZeqNTracksPVvsPt", "ZeqNTracksPV distribution of events with hight pT hadron", HistType::kTH2F, {{multAxisNTPV}, {11, 0, 11}});

      // additional QA histos
      QAHistos.add("hTPCNsigmaXiBachPiPlus", "nsigma TPC distribution bachelor pion+", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaXiV0PiPlus", "nsigma TPC distribution pi+", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaXiV0Proton", "nsigma TPC distribution proton", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaXiBachPiMinus", "nsigma TPC distribution bachelor pion-", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaXiV0PiMinus", "nsigma TPC distribution pi-", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaXiV0AntiProton", "nsigma TPC distribution antiproton", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaOmegaBachKaPlus", "nsigma TPC distribution bachelor kaon+", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaOmegaV0PiPlus", "nsigma TPC distribution pi+", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaOmegaV0Proton", "nsigma TPC distribution proton", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaOmegaBachKaMinus", "nsigma TPC distribution bachelor kaon-", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaOmegaV0PiMinus", "nsigma TPC distribution pi-", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hTPCNsigmaOmegaV0AntiProton", "nsigma TPC distribution antiproton", HistType::kTH2F, {{80, -10, 10}, {pTPCAxis}});
      QAHistos.add("hHasTOFBachKa", "bachelor kaon has TOF", HistType::kTH2F, {{2, 0, 2}, {ptAxis}});
      QAHistos.add("hHasTOFBachPi", "bachelor pi has TOF", HistType::kTH2F, {{2, 0, 2}, {ptAxis}});
      QAHistos.add("hHasTOFPr", "pr dau has TOF", HistType::kTH2F, {{2, 0, 2}, {ptAxis}});
      QAHistos.add("hHasTOFPi", "pi dau has TOF", HistType::kTH2F, {{2, 0, 2}, {ptAxis}});
      QAHistos.add("hRapXi", "Rap Xi", HistType::kTH1F, {{100, -1, 1}});
      QAHistos.add("hRapOmega", "Rap Omega", HistType::kTH1F, {{100, -1, 1}});

      QAHistosStrangenessTracking.add("hStRVsPtTrkCasc", "Tracked cascades;p_{T} (GeV/#it{c});R (cm)", HistType::kTH2D, {{200, 0., 10.}, {200, 0., 50}});
      QAHistosStrangenessTracking.add("hMassOmegaTrkCasc", "Tracked cascades;m_{#Omega} (GeV/#it{c}^{2})", HistType::kTH1D, {{1000, 1., 3.}});
      QAHistosStrangenessTracking.add("hMassXiTrkCasc", "Tracked cascades;m_{#Xi} (GeV/#it{c}^{2})", HistType::kTH1D, {{1000, 1., 3.}});
      QAHistosStrangenessTracking.add("hMassV0TrkCasc", "Tracked cascades;m_{V^{0}} (GeV/#it{c}^{2})", HistType::kTH1D, {{1000, 1., 3.}});
      QAHistosStrangenessTracking.add("hMatchChi2TrkCasc", "Tracked cascades;#chi^{2}", HistType::kTH1D, {{1000, 0., 2000.}});
      QAHistosStrangenessTracking.add("hMatchChi2TrkCascSelectedXi", "Tracked cascades;#chi^{2}", HistType::kTH1D, {{1000, 0., 2000.}});
      QAHistosStrangenessTracking.add("hMatchChi2TrkCascSelectedOmega", "Tracked cascades;#chi^{2}", HistType::kTH1D, {{1000, 0., 2000.}});
      QAHistosStrangenessTracking.add("hMassOmegaVsMatchChi2TrkCasc", "Tracked cascades;m_{#Omega} (GeV/#it{c}^{2});#chi^{2}", HistType::kTH2D, {{1000, 1., 3.}, {1000, 0., 2000.}});
      QAHistosStrangenessTracking.add("hMassXiVsMatchChi2TrkCasc", "Tracked cascades;m_{#Xi} (GeV/#it{c}^{2});#chi^{2}", HistType::kTH2D, {{1000, 1., 3.}, {1000, 0., 2000.}});
      QAHistosStrangenessTracking.add("hMassOmegaVsTopChi2TrkCasc", "Tracked cascades;m_{#Omega} (GeV/#it{c}^{2});#chi^{2}", HistType::kTH2D, {{1000, 1., 3.}, {1000, 0., 2000.}});
      QAHistosStrangenessTracking.add("hMassXiVsTopChi2TrkCasc", "Tracked cascades;m_{#Xi} (GeV/#it{c}^{2});#chi^{2}", HistType::kTH2D, {{1000, 1., 3.}, {1000, 0., 2000.}});
      QAHistosStrangenessTracking.add("hNSigmaTpcPiTrkCascBachelor", "Tracked cascades;N_{#sigma, #pi}", HistType::kTH1D, {{100, -5., 5.}});
      QAHistosStrangenessTracking.add("hNSigmaTpcKaTrkCascBachelor", "Tracked cascades;N_{#sigma, K}", HistType::kTH1D, {{100, -5., 5.}});
      QAHistosStrangenessTracking.add("hNSigmaTpcPrTrkCascV0", "Tracked cascades;N_{#sigma, p}", HistType::kTH1D, {{100, -5., 5.}});
      QAHistosStrangenessTracking.add("hNSigmaTpcPiTrkCascV0", "Tracked cascades;N_{#sigma, #pi}", HistType::kTH1D, {{100, -5., 5.}});
      QAHistosStrangenessTracking.add("hNSigmaTpcPrTrkCascV0SelectedXi", "Tracked cascades;N_{#sigma, p}", HistType::kTH1D, {{100, -5., 5.}});
      QAHistosStrangenessTracking.add("hNSigmaTpcPiTrkCascV0SelectedXi", "Tracked cascades;N_{#sigma, #pi}", HistType::kTH1D, {{100, -5., 5.}});
      QAHistosStrangenessTracking.add("hNSigmaTpcPiTrkCascBachelorSelectedXi", "Tracked cascades;N_{#sigma, #pi}", HistType::kTH1D, {{100, -5., 5.}});
      QAHistosStrangenessTracking.add("hNSigmaTpcPrTrkCascV0SelectedOmega", "Tracked cascades;N_{#sigma, p}", HistType::kTH1D, {{100, -5., 5.}});
      QAHistosStrangenessTracking.add("hNSigmaTpcPiTrkCascV0SelectedOmega", "Tracked cascades;N_{#sigma, #pi}", HistType::kTH1D, {{100, -5., 5.}});
      QAHistosStrangenessTracking.add("hNSigmaTpcKaTrkCascBachelorSelectedOmega", "Tracked cascades;N_{#sigma, K}", HistType::kTH1D, {{100, -5., 5.}});
      QAHistosStrangenessTracking.add("hMassH3LTrkV0", "Tracked V0;m_{H3L} (GeV/#it{c}^{2})", HistType::kTH1D, {{1000, 2.8, 3.8}});
      QAHistosStrangenessTracking.add("hMassH4LTrkV0", "Tracked V0;m_{H4L} (GeV/#it{c}^{2})", HistType::kTH1D, {{1000, 3.8, 4.8}});
      QAHistosStrangenessTracking.add("hMassH3LTrk3body", "Tracked 3body;m_{H3L} (GeV/#it{c}^{2})", HistType::kTH1D, {{200, 0., 10.}});
      QAHistosStrangenessTracking.add("hMassHe4LTrk3body", "Tracked 3body;m_{He4L} (GeV/#it{c}^{2})", HistType::kTH1D, {{200, 0., 10.}});
      QAHistosStrangenessTracking.add("hDcaXY", "DCA;DCA_{xy} (cm)", HistType::kTH1D, {{200, -.5, .5}});
      QAHistosStrangenessTracking.add("hDcaXYSelectedXi", "DCA;DCA_{xy} (cm)", HistType::kTH1D, {{200, -.05, .05}});
      QAHistosStrangenessTracking.add("hDcaXYSelectedOmega", "DCA;DCA_{xy} (cm)", HistType::kTH1D, {{200, -.05, .05}});
      QAHistosStrangenessTracking.add("hDcaXYVsPt", "DCA;p_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2D, {{200, 0., 10.}, {200, -.5, .5}});
      QAHistosStrangenessTracking.add("hDcaXYVsPtSelectedXi", "DCA;p_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2D, {{200, 0., 10.}, {200, -.5, .5}});
      QAHistosStrangenessTracking.add("hDcaXYVsPtSelectedOmega", "DCA;p_{T} (GeV/#it{c});DCA_{xy} (cm)", HistType::kTH2D, {{200, 0., 10.}, {200, -.5, .5}});
      QAHistosStrangenessTracking.add("hDcaZ", "DCA;DCA_{z} (cm)", HistType::kTH1D, {{200, -.5, .5}});
      QAHistosStrangenessTracking.add("hDcaZSelectedXi", "DCA;DCA_{z} (cm)", HistType::kTH1D, {{200, -.05, .05}});
      QAHistosStrangenessTracking.add("hDcaZSelectedOmega", "DCA;DCA_{z} (cm)", HistType::kTH1D, {{200, -.05, .05}});
      QAHistosStrangenessTracking.add("hDcaZVsPt", "DCA;p_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2D, {{200, 0., 10.}, {200, -.5, .5}});
      QAHistosStrangenessTracking.add("hDcaZVsPtSelectedXi", "DCA;p_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2D, {{200, 0., 10.}, {200, -.5, .5}});
      QAHistosStrangenessTracking.add("hDcaZVsPtSelectedOmega", "DCA;p_{T} (GeV/#it{c});DCA_{z} (cm)", HistType::kTH2D, {{200, 0., 10.}, {200, -.5, .5}});
      QAHistosStrangenessTracking.add("hDcaVsPt", "DCA;DCA (cm);p_{T} (GeV/#it{c})", HistType::kTH2D, {{200, 0., .5}, {200, 0., 10.}});
      QAHistosStrangenessTracking.add("hDcaVsR", "DCA;DCA (cm);R (cm)", HistType::kTH2D, {{200, 0., .5}, {200, 0., 10.}});
      QAHistosStrangenessTracking.add("hDecayRadius", "Decay radius;R (cm)", HistType::kTH1D, {{100, 0., 30.}});
      QAHistosStrangenessTracking.add("hDecayRadiusSelectedXi", "Decay radius;R (cm)", HistType::kTH1D, {{100, 0., 30.}});
      QAHistosStrangenessTracking.add("hDecayRadiusSelectedOmega", "Decay radius;R (cm)", HistType::kTH1D, {{100, 0., 30.}});
      QAHistosStrangenessTracking.add("hDecayRadiusVsXiMass", "Decay radius;R (cm);m (GeV/#it{c}^2)", HistType::kTH2D, {{100, 0., 30.}, {1000, 1., 2.}});
      QAHistosStrangenessTracking.add("hDecayRadiusVsXiMassSelected", "Decay radius;R (cm);m (GeV/#it{c}^2)", HistType::kTH2D, {{100, 0., 30.}, {1000, 1., 2.}});
      QAHistosStrangenessTracking.add("hDecayRadiusVsOmegaMass", "Decay radius;R (cm);m (GeV/#it{c}^2)", HistType::kTH2D, {{100, 0., 30.}, {1000, 1., 2.}});
      QAHistosStrangenessTracking.add("hDecayRadiusVsOmegaMassSelected", "Decay radius;R (cm);m (GeV/#it{c}^2)", HistType::kTH2D, {{100, 0., 30.}, {1000, 1., 2.}});
      QAHistosStrangenessTracking.add("hCpa", "cpa;cpa", HistType::kTH1D, {{500, .995, 1.}});
      QAHistosStrangenessTracking.add("hCpaSelectedXi", "cpa;cpa", HistType::kTH1D, {{500, .995, 1.}});
      QAHistosStrangenessTracking.add("hCpaSelectedOmega", "cpa;cpa", HistType::kTH1D, {{500, .995, 1.}});
      QAHistosStrangenessTracking.add("hPtCascCand", "cascades;p_{T} (GeV/#it{c})", HistType::kTH1D, {{200, 0., 10.}});
      QAHistosStrangenessTracking.add("hPtCascTracked", "tracked cascades;p_{T} (GeV/#it{c})", HistType::kTH1D, {{200, 0., 10.}});
      QAHistosStrangenessTracking.add("hPtVsMassTrkXi", "cascades;p_{T} (GeV/#it{c});m (GeV/#it{c}^2)", HistType::kTH2D, {{200, 0., 10.}, {1000, 1.2, 1.7}});
      QAHistosStrangenessTracking.add("hPtVsMassTrkOmega", "cascades;p_{T} (GeV/#it{c});m (GeV/#it{c}^2)", HistType::kTH2D, {{200, 0., 10.}, {1000, 1.6, 2.1}});
      QAHistosStrangenessTracking.add("hPtVsMassTrkXiSelected", "cascades;p_{T} (GeV/#it{c});m (GeV/#it{c}^2)", HistType::kTH2D, {{200, 0., 10.}, {1000, 1.2, 1.7}});
      QAHistosStrangenessTracking.add("hPtVsMassTrkOmegaSelected", "cascades;p_{T} (GeV/#it{c});m (GeV/#it{c}^2)", HistType::kTH2D, {{200, 0., 10.}, {1000, 1.6, 2.1}});

      // Sigma QA histograms
      QAHistosSigma.add("hPtVsMassSigmaPlus", ";p_{T} (GeV/#it{c});m (GeV/#it{c}^2)", HistType::kTH2D, {{20, -5, 5.}, {50, 1.1, 1.3}});
      QAHistosSigma.add("hDecayRadiusSigma", ";R (cm)", HistType::kTH1D, {{100, 0., 40.}});
      QAHistosSigma.add("hPtNSigmaTPCPrPair", ";p_{T} (GeV/#it{c});#sigma_{TPC}", HistType::kTH2D, {{200, -5., 5.}, {200, -5., 5.}});
      QAHistosSigma.add("hPtNSigmaTOFPrPair", ";p_{T} (GeV/#it{c});#sigma_{TOF}", HistType::kTH2D, {{200, -5., 5.}, {200, -5., 5.}});
      QAHistosSigma.add("hKStarSigmaPr", ";k*", HistType::kTH1D, {{200, 0., 2}});
    }
  }

  void initCCDB(int run)
  {
    if (run != mRunNumber) {
      mRunNumber = run;
      o2::parameters::GRPMagField* grpmag = ccdb->getForRun<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", run);
      o2::base::Propagator::initFieldFromGRP(grpmag);
      mBz = static_cast<float>(grpmag->getNominalL3Field());
      if (cfgHMOmegaCuts.HMTrgSelectionForOmegaFT0M == 1) {
        mMeanMultT0C = ccdb->getForRun<std::vector<double>>("Users/e/ekryshen/meanT0C", run);
        mMeanMultT0A = ccdb->getForRun<std::vector<double>>("Users/e/ekryshen/meanT0A", run);
      }
      mDCAFitter.setBz(mBz);
      mDCAFitter.setPropagateToPCA(propToDCA);
      mDCAFitter.setMaxR(maxR);
      mDCAFitter.setMaxDZIni(maxDZIni);
      mDCAFitter.setMinParamChange(minParamChange);
      mDCAFitter.setMinRelChi2Change(minRelChi2Change);
      mDCAFitter.setUseAbsDCA(useAbsDCA);
      mStraHelper.fitter.setBz(mBz);
    }
    if (!mStraHelper.lut) { /// done only once
      ccdb->setURL(ccdbUrl);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setFatalWhenNull(true);
      auto* lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>("GLO/Param/MatLUT"));
      o2::base::Propagator::Instance()->setMatLUT(lut);
      mStraHelper.lut = lut;
    }
  }

  // Tables
  using CollisionCandidates = soa::Join<aod::Collisions, aod::EvSels, aod::MultZeqs, aod::FT0Mults, aod::MultsGlobal>::iterator;
  using TrackCandidates = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::pidTPCLfFullPi, aod::pidTPCLfFullPr, aod::pidTPCLfFullKa, aod::pidTOFFullPr>;

  float getMassWindow(const stfilter::species s, const float pt, const float nsigma = 6)
  {
    const auto sigma = parSigmaMass->get(0u, s) * exp(parSigmaMass->get(1, s) * pt) + parSigmaMass->get(2, s) * exp(parSigmaMass->get(3, s) * pt);
    return nsigma * sigma;
  }

  ////////////////////////////////////////////////////////
  ////////// Strangeness Filter //////////////////////////
  ////////////////////////////////////////////////////////

  void fillTriggerTable(bool keepEvent[])
  {
    strgtable(keepEvent[0], keepEvent[1], keepEvent[2], keepEvent[3], keepEvent[4], keepEvent[5], keepEvent[6], keepEvent[7], keepEvent[8], keepEvent[9], keepEvent[10], keepEvent[11], keepEvent[12], keepEvent[13], keepEvent[14], keepEvent[15], keepEvent[16]);
  }

  void process(CollisionCandidates const& collision, TrackCandidates const& tracks, aod::Cascades const& cascadesBase, aod::AssignedTrackedCascades const& trackedCascades, aod::KinkCands const& kinkCands, aod::AssignedTrackedV0s const& /*trackedV0s*/, aod::AssignedTracked3Bodys const& /*tracked3Bodys*/, aod::V0s const& v0Base, aod::BCs const&, aod::FT0s const& /*ft0s*/)
  {
    // Is event good? [0] = Omega, [1] = high-pT hadron + Omega, [2] = 2Xi, [3] = 3Xi, [4] = 4Xi, [5] single-Xi, [6] Omega with high radius
    // [7] tracked Xi, [8] tracked Omega, [9] Omega + high mult event
    bool keepEvent[17]{}; // explicitly zero-initialised
    std::vector<std::array<int64_t, 2>> v0sFromOmegaID;
    std::vector<std::array<int64_t, 2>> v0sFromXiID;

    initCCDB(collision.bc().runNumber());

    hProcessedEvents->Fill(-0.5);
    if (sel8 && !collision.sel8()) {
      fillTriggerTable(keepEvent);
      return;
    }
    if (isTriggerTVX && !collision.selection_bit(aod::evsel::kIsTriggerTVX)) {
      fillTriggerTable(keepEvent);
      return;
    }
    if (isTimeFrameBorderCut && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      fillTriggerTable(keepEvent);
      return;
    }
    // all processed events after event selection
    hProcessedEvents->Fill(0.5);

    if (std::fabs(collision.posZ()) > cutzvertex) {
      fillTriggerTable(keepEvent);
      return;
    }
    QAHistos.fill(HIST("hVtxZ"), collision.posZ());
    if (doextraQA) {
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityZeqV0A"), collision.multZeqFV0A());
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityZeqT0M"), collision.multZeqFT0A() + collision.multZeqFT0C());
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityZeqNTracksPV"), collision.multZeqNTracksPV());
    }

    Bool_t isHighMultEvent = 0;         // tail
    Bool_t isHighMultEventOmegaCut = 0; // Omega HM cut
    float multFT0MNorm = 0.f;
    Bool_t isHighMultEventTrk = 0;         // tail
    Bool_t isHighMultEventTrkOmegaCut = 0; // Omega HM cut

    float multTrack = 0.f;
    if (cfgHMOmegaCuts.HMTrgSelectionForOmegaFT0M == 1) {
      float meanMultT0C = 0.f;
      float fac_FT0C_ebe = 1.;
      meanMultT0C = (*mMeanMultT0C)[0];
      if (meanMultT0C > 0) {
        fac_FT0C_ebe = avPyT0C / meanMultT0C;
      }
      float meanMultT0A = 0.f;
      meanMultT0A = (*mMeanMultT0A)[0];
      float fac_FT0A_ebe = 1.;
      if (meanMultT0A > 0) {
        fac_FT0A_ebe = avPyT0A / meanMultT0A;
      }
      LOG(debug) << "Mean mults t0:" << fac_FT0A_ebe << " " << fac_FT0C_ebe;
      if (collision.has_foundFT0()) {
        static int ampneg = 0;
        auto ft0 = collision.foundFT0();
        float sumAmpFT0C = 0.f;
        for (std::size_t i_c = 0; i_c < ft0.amplitudeC().size(); i_c++) {
          float amplitude = ft0.amplitudeC()[i_c];
          sumAmpFT0C += amplitude;
        }
        float sumAmpFT0A = 0.f;
        for (std::size_t i_a = 0; i_a < ft0.amplitudeA().size(); i_a++) {
          float amplitude = ft0.amplitudeA()[i_a];
          sumAmpFT0A += amplitude;
        }
        const int nEta5 = 2; // FT0C + FT0A
        float weigthsEta5[nEta5] = {0.0490638, 0.010958415};
        if (sumAmpFT0C >= 0 || sumAmpFT0A >= 0) {
          if (meanMultT0A > 0 && meanMultT0C > 0) {
            multFT0MNorm = sumAmpFT0C * fac_FT0C_ebe + sumAmpFT0A * fac_FT0A_ebe;
          } else {
            multFT0MNorm = sumAmpFT0C * weigthsEta5[0] + sumAmpFT0A * weigthsEta5[1];
          }
          LOG(debug) << "meanMult:" << multFT0MNorm << " multFT0M:" << collision.multFT0M();
          if (sumAmpFT0A < 0 || sumAmpFT0C < 0) {
            // LOG(info) << "ampa: " << sumAmpFT0A << " ampc:" << sumAmpFT0C;
            ampneg++;
          }
          EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityFT0MNorm"), multFT0MNorm);
          if (multFT0MNorm > cfgHMOmegaCuts.LowLimitHMTrgOmegaT0MNorm) {
            isHighMultEventOmegaCut = 1;
            LOG(debug) << "Found FT0 using norm mult";
          }
          if (multFT0MNorm > cfgHMOmegaCuts.LowLimitHMTrgT0MNorm) {
            isHighMultEvent = 1;
          }
        } else {
          LOG(warn) << "Found FT0 but, bith amplitudes are <=0 ";
          EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityFT0MNorm"), 148);
          EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityFT0MNoFT0"), collision.multFT0M());
        }
        if (ampneg) {
          LOG(warn) << "# of negative amplitudes:" << ampneg;
        }
      } else {
        LOG(debug) << "FT0 not Found, using FT0M";
        EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityFT0MNorm"), 149);
        EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityFT0MNoFT0"), collision.multFT0M());
      }
    } else if (cfgHMOmegaCuts.HMTrgSelectionForOmegaFT0M == 2) {
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityFT0M"), collision.multFT0M());
      if (collision.multFT0M() > cfgHMOmegaCuts.LowLimitHMTrgOmegaT0M) {
        isHighMultEventOmegaCut = 1;
      }
    }
    if (cfgHMOmegaCuts.HMTrgSelectionForOmegaTrks == 1) {
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityTracksGlob"), collision.multNTracksGlobal());
      if (collision.multNTracksGlobal() > cfgHMOmegaCuts.LowLimitHMTrgOmegaTrkGlob) {
        isHighMultEventTrkOmegaCut = 1;
      }
      if (collision.multNTracksGlobal() > cfgHMOmegaCuts.LowLimitHMTrgTrkGlob) {
        isHighMultEventTrk = 1;
      }
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityTracks"), collision.multNTracksGlobal());

    } else if (cfgHMOmegaCuts.HMTrgSelectionForOmegaTrks == 2) {
      for (auto& track : tracks) {
        if (selectTrackOHM(track)) {
          multTrack++;
        }
      }
      if (multTrack > cfgHMOmegaCuts.LowLimitHMTrgOmegaTrkSel) {
        isHighMultEventTrkOmegaCut = 1;
      }
      if (multTrack > cfgHMOmegaCuts.LowLimitHMTrgTrkSel) {
        isHighMultEventTrk = 1;
      }
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityTracks"), multTrack);
    }
    // constants
    const float ctauxi = 4.91;     // from PDG
    const float ctauomega = 2.461; // from PDG

    // variables
    float xipos = -1.;
    float xiproperlifetime = -1.;
    float omegaproperlifetime = -1.;
    float xiptotmom = -1.;
    int xicounter = 0;
    int xicounterYN = 0;
    int omegacounter = 0;
    int omegalargeRcounter = 0;
    const std::array<float, 3> pvPos{collision.posX(), collision.posY(), collision.posZ()};
    float pvX = 0.0f, pvY = 0.0f, pvZ = 0.0f;

    // strangeness tracking selection
    const auto primaryVertex = getPrimaryVertex(collision);
    o2::dataformats::DCA impactParameterTrk;

    std::vector<std::tuple<int64_t, int64_t, TVector3, TVector3>> v0sSelTuple;
    for (auto& v00 : v0Base) { // loop over v0 for pre selection
      hCandidate->Fill(0.5);   // All candidates

      if (v00.v0Type() != 1) {
        continue;
      }

      const auto posTrack0 = v00.posTrack_as<TrackCandidates>();
      const auto negTrack0 = v00.negTrack_as<TrackCandidates>();

      if (!isSelectedV0Daughter(posTrack0) || !isSelectedV0Daughter(negTrack0)) {
        continue;
      }

      auto trackParPos0 = getTrackParCov(posTrack0);
      auto trackParNeg0 = getTrackParCov(negTrack0);

      if (!mStraHelper.buildV0Candidate(v00.collisionId(), pvPos[0], pvPos[1], pvPos[2], posTrack0, negTrack0, trackParPos0, trackParNeg0)) {
        continue;
      }

      if (std::hypot(mStraHelper.v0.position[0], mStraHelper.v0.position[1]) < cfgLLCuts.cfgv0radiusMin) {
        continue;
      }
      if (std::fabs(mStraHelper.v0.positiveDCAxy) < cfgLLCuts.cfgDCAPosToPVMin) {
        continue;
      }
      if (std::fabs(mStraHelper.v0.negativeDCAxy) < cfgLLCuts.cfgDCANegToPVMin) {
        continue;
      }
      if (TMath::Cos(mStraHelper.v0.pointingAngle) < cfgLLCuts.cfgv0CosPA) {
        continue;
      }
      if (std::fabs(mStraHelper.v0.daughterDCA) > cfgLLCuts.cfgDCAV0Dau) {
        continue;
      }
      if (std::hypot(mStraHelper.v0.momentum[0], mStraHelper.v0.momentum[1]) < cfgLLCuts.cfgV0PtMin) {
        continue;
      }
      double yLambda = RecoDecay::y(array{mStraHelper.v0.momentum[0], mStraHelper.v0.momentum[1], mStraHelper.v0.momentum[2]}, o2::constants::physics::MassLambda0);
      if (yLambda < cfgLLCuts.cfgV0RapMin) {
        continue;
      }
      if (yLambda > cfgLLCuts.cfgV0RapMax) {
        continue;
      }
      double distovertotmom = std::hypot(mStraHelper.v0.position[0] - collision.posX(), mStraHelper.v0.position[1] - collision.posY(), mStraHelper.v0.position[2] - collision.posZ()) / (std::hypot(mStraHelper.v0.momentum[0], mStraHelper.v0.momentum[1], mStraHelper.v0.momentum[2]) + 1e-13);
      if (distovertotmom * o2::constants::physics::MassLambda0 > cfgLLCuts.cfgV0LifeTime) {
        continue;
      }

      int Tag = 0;
      if (isSelectedV0DaughterPID(posTrack0, 0) && isSelectedV0DaughterPID(negTrack0, 1)) {
        if (cfgLLCuts.cfgLambdaMassWindow > std::fabs(mStraHelper.v0.massLambda - o2::constants::physics::MassLambda0)) {
          if (cfgLLCuts.cfgCompV0Rej < std::fabs(mStraHelper.v0.massK0Short - o2::constants::physics::MassLambda0)) {
            Tag++;
          }
        }
      } // lambda
      if (isSelectedV0DaughterPID(posTrack0, 1) && isSelectedV0DaughterPID(negTrack0, 0)) {
        if (cfgLLCuts.cfgLambdaMassWindow > std::fabs(mStraHelper.v0.massAntiLambda - o2::constants::physics::MassLambda0)) {
          if (cfgLLCuts.cfgCompV0Rej < std::fabs(mStraHelper.v0.massK0Short - o2::constants::physics::MassLambda0)) {
            Tag++;
          }
        }
      } // anti lambda
      if (Tag != 1) { // Select when only one hypothesis is satisfied
        continue;
      }

      TVector3 v0pos(mStraHelper.v0.position[0], mStraHelper.v0.position[1], mStraHelper.v0.position[2]);
      TVector3 v0mom(mStraHelper.v0.momentum[0], mStraHelper.v0.momentum[1], mStraHelper.v0.momentum[2]);

      v0sSelTuple.emplace_back(posTrack0.globalIndex(), negTrack0.globalIndex(), v0pos, v0mom);
    }

    for (size_t i = 0; i < v0sSelTuple.size(); ++i) {
      for (size_t j = i + 1; j < v0sSelTuple.size(); ++j) {
        auto d00 = std::get<0>(v0sSelTuple[i]);
        auto d01 = std::get<1>(v0sSelTuple[i]);
        auto d10 = std::get<0>(v0sSelTuple[j]);
        auto d11 = std::get<1>(v0sSelTuple[j]);
        if (d00 == d10 || d00 == d11 || d01 == d10 || d01 == d11) {
          continue;
        }
        auto v00pos = std::get<2>(v0sSelTuple[i]);
        auto v00mom = std::get<3>(v0sSelTuple[i]);
        auto v01pos = std::get<2>(v0sSelTuple[j]);
        auto v01mom = std::get<3>(v0sSelTuple[j]);
        if (isSelectedV0V0(v00pos, v00mom, v01pos, v01mom)) {
          keepEvent[12] = true;
        }
      }
    }

    for (auto& casc : cascadesBase) { // loop over cascades
      hCandidate->Fill(0.5);          // All candidates

      const auto bachTrack = casc.bachelor_as<TrackCandidates>();
      const auto v0Dau = casc.v0_as<o2::aod::V0s>();
      const auto negTrack = v0Dau.negTrack_as<TrackCandidates>();
      const auto posTrack = v0Dau.posTrack_as<TrackCandidates>();

      if (!mStraHelper.buildCascadeCandidate(casc.collisionId(), pvPos[0], pvPos[1], pvPos[2], posTrack, negTrack, bachTrack, -1, useCascadeMomentumAtPrimVtx, -1)) {
        continue;
      }
      hCandidate->Fill(1.5); // Built and selected candidates in StraBuilder

      bool isXi = false;
      bool isXiYN = false;
      bool isOmega = false;
      bool isOmegalargeR = false;

      // QA
      double massXi = mStraHelper.cascade.massXi;
      double massOmega = mStraHelper.cascade.massOmega;
      double ptCasc = RecoDecay::sqrtSumOfSquares(mStraHelper.cascade.cascadeMomentum[0], mStraHelper.cascade.cascadeMomentum[1]);
      QAHistos.fill(HIST("hMassXiBefSelvsPt"), massXi, ptCasc);
      QAHistos.fill(HIST("hMassOmegaBefSelvsPt"), massOmega, ptCasc);
      // Position
      xipos = std::hypot(mStraHelper.cascade.cascadePosition[0] - collision.posX(), mStraHelper.cascade.cascadePosition[1] - collision.posY(), mStraHelper.cascade.cascadePosition[2] - collision.posZ());
      // Total momentum
      xiptotmom = std::hypot(mStraHelper.cascade.cascadeMomentum[0], mStraHelper.cascade.cascadeMomentum[1], mStraHelper.cascade.cascadeMomentum[2]);
      // Proper lifetime
      xiproperlifetime = o2::constants::physics::MassXiMinus * xipos / (xiptotmom + 1e-13);
      omegaproperlifetime = o2::constants::physics::MassOmegaMinus * xipos / (xiptotmom + 1e-13);
      // Radii
      double Cascv0radius = std::hypot(mStraHelper.cascade.v0Position[0], mStraHelper.cascade.v0Position[1]);
      double Casccascradius = std::hypot(mStraHelper.cascade.cascadePosition[0], mStraHelper.cascade.cascadePosition[1]);
      // Rapidity
      double etaCasc = RecoDecay::eta(std::array{mStraHelper.cascade.cascadeMomentum[0], mStraHelper.cascade.cascadeMomentum[1], mStraHelper.cascade.cascadeMomentum[2]});
      // pointing angle
      double v0DauCPA = RecoDecay::cpa(pvPos, array{mStraHelper.cascade.v0Position[0], mStraHelper.cascade.v0Position[1], mStraHelper.cascade.v0Position[2]}, array{mStraHelper.cascade.positiveMomentum[0] + mStraHelper.cascade.negativeMomentum[0], mStraHelper.cascade.positiveMomentum[1] + mStraHelper.cascade.negativeMomentum[1], mStraHelper.cascade.positiveMomentum[2] + mStraHelper.cascade.negativeMomentum[2]});
      double cascCPA = RecoDecay::cpa(
        pvPos,
        array{mStraHelper.cascade.cascadePosition[0], mStraHelper.cascade.cascadePosition[1], mStraHelper.cascade.cascadePosition[2]},
        array{mStraHelper.cascade.positiveMomentum[0] + mStraHelper.cascade.negativeMomentum[0] + mStraHelper.cascade.bachelorMomentum[0], mStraHelper.cascade.positiveMomentum[1] + mStraHelper.cascade.negativeMomentum[1] + mStraHelper.cascade.bachelorMomentum[1], mStraHelper.cascade.positiveMomentum[2] + mStraHelper.cascade.negativeMomentum[2] + mStraHelper.cascade.bachelorMomentum[2]});
      // dca V0 to PV
      double DCAV0ToPV = CalculateDCAStraightToPV(
        mStraHelper.cascade.v0Position[0], mStraHelper.cascade.v0Position[1], mStraHelper.cascade.v0Position[2],
        mStraHelper.cascade.positiveMomentum[0] + mStraHelper.cascade.negativeMomentum[0],
        mStraHelper.cascade.positiveMomentum[1] + mStraHelper.cascade.negativeMomentum[1],
        mStraHelper.cascade.positiveMomentum[2] + mStraHelper.cascade.negativeMomentum[2],
        pvX, pvY, pvZ);
      // massLambda
      double LambdaMass = 0;
      if (mStraHelper.cascade.charge < 0) {
        LambdaMass = RecoDecay::m(array{array{mStraHelper.cascade.positiveMomentum[0], mStraHelper.cascade.positiveMomentum[1], mStraHelper.cascade.positiveMomentum[2]}, array{mStraHelper.cascade.negativeMomentum[0], mStraHelper.cascade.negativeMomentum[1], mStraHelper.cascade.negativeMomentum[2]}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged});
      } else {
        LambdaMass = RecoDecay::m(array{array{mStraHelper.cascade.positiveMomentum[0], mStraHelper.cascade.positiveMomentum[1], mStraHelper.cascade.positiveMomentum[2]}, array{mStraHelper.cascade.negativeMomentum[0], mStraHelper.cascade.negativeMomentum[1], mStraHelper.cascade.negativeMomentum[2]}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton});
      }

      // rapidity
      double yXi = RecoDecay::y(array{mStraHelper.cascade.bachelorMomentum[0] + mStraHelper.cascade.positiveMomentum[0] + mStraHelper.cascade.negativeMomentum[0], mStraHelper.cascade.bachelorMomentum[1] + mStraHelper.cascade.positiveMomentum[1] + mStraHelper.cascade.negativeMomentum[1], mStraHelper.cascade.bachelorMomentum[2] + mStraHelper.cascade.positiveMomentum[2] + mStraHelper.cascade.negativeMomentum[2]}, o2::constants::physics::MassXiMinus);
      double yOmega = RecoDecay::y(array{mStraHelper.cascade.bachelorMomentum[0] + mStraHelper.cascade.positiveMomentum[0] + mStraHelper.cascade.negativeMomentum[0], mStraHelper.cascade.bachelorMomentum[1] + mStraHelper.cascade.positiveMomentum[1] + mStraHelper.cascade.negativeMomentum[1], mStraHelper.cascade.bachelorMomentum[2] + mStraHelper.cascade.positiveMomentum[2] + mStraHelper.cascade.negativeMomentum[2]}, o2::constants::physics::MassOmegaMinus);

      if (mStraHelper.cascade.charge > 0) {
        if (std::fabs(mStraHelper.cascade.positiveDCAxy) < dcamesontopv) {
          continue;
        }
        hCandidate->Fill(2.5);
        if (std::fabs(mStraHelper.cascade.negativeDCAxy) < dcabaryontopv) {
          continue;
        }
        hCandidate->Fill(3.5);
        if (std::fabs(posTrack.tpcNSigmaPi()) > nsigmatpcpi) {
          continue;
        }
        hCandidate->Fill(4.5);
        if (std::fabs(negTrack.tpcNSigmaPr()) > nsigmatpcpr) {
          continue;
        }
        hCandidate->Fill(5.5);
      } else if (mStraHelper.cascade.charge < 0) {
        if (std::fabs(mStraHelper.cascade.negativeDCAxy) < dcamesontopv) {
          continue;
        }
        hCandidate->Fill(2.5);
        if (std::fabs(mStraHelper.cascade.positiveDCAxy) < dcabaryontopv) {
          continue;
        }
        hCandidate->Fill(3.5);
        if (std::fabs(negTrack.tpcNSigmaPi()) > nsigmatpcpi) {
          continue;
        }
        hCandidate->Fill(4.5);
        if (std::fabs(posTrack.tpcNSigmaPr()) > nsigmatpcpr) {
          continue;
        }
        hCandidate->Fill(5.5);
      }
      hCandidate->Fill(6.5); // OLD: eta dau (selection now applied in strangeness helper)
      hCandidate->Fill(7.5); // OLD: bachtopv (selection now applied in strangeness helper)

      // not striclty needed as selection are applied beforehand - just as QA (no change in number expected)
      if (Cascv0radius < v0radius) {
        continue;
      }
      hCandidate->Fill(8.5);
      if (Casccascradius < cascradius) {
        continue;
      }
      hCandidate->Fill(9.5);
      if (v0DauCPA < v0cospa) {
        continue;
      }
      hCandidate->Fill(10.5);
      if (mStraHelper.cascade.v0DaughterDCA > dcav0dau) {
        continue;
      }
      hCandidate->Fill(11.5);
      if (mStraHelper.cascade.cascadeDaughterDCA > dcacascdau) {
        continue;
      }
      hCandidate->Fill(12.5);
      if (std::fabs(LambdaMass - constants::physics::MassLambda) > masslambdalimit) {
        continue;
      }
      hCandidate->Fill(13.5);
      if (std::fabs(etaCasc) > eta) {
        continue;
      }
      hCandidate->Fill(14.5);
      if (hastof &&
          (!posTrack.hasTOF() && posTrack.pt() > ptthrtof) &&
          (!negTrack.hasTOF() && negTrack.pt() > ptthrtof) &&
          (!bachTrack.hasTOF() && bachTrack.pt() > ptthrtof)) {
        continue;
      }
      hCandidate->Fill(15.5);

      // Fill selections QA for Xi
      if (cascCPA > casccospaxi) {
        hCandidate->Fill(16.5);
        if (cascCPA > dcav0topv) {
          hCandidate->Fill(17.5);
          if (xiproperlifetime < properlifetimefactor * ctauxi) {
            hCandidate->Fill(18.5);
            if (std::fabs(yXi) < rapidity) {
              hCandidate->Fill(19.5);
            }
          }
        }
      }

      const auto deltaMassXi = useSigmaBasedMassCutXi ? getMassWindow(stfilter::species::Xi, ptCasc) : ximasswindow;
      const auto deltaMassOmega = useSigmaBasedMassCutOmega ? getMassWindow(stfilter::species::Omega, ptCasc) : omegamasswindow;

      isXi = (std::fabs(bachTrack.tpcNSigmaPi()) < nsigmatpcpi) &&
             (cascCPA > casccospaxi) &&
             (DCAV0ToPV > dcav0topv) &&
             (std::fabs(massXi - o2::constants::physics::MassXiMinus) < deltaMassXi) &&
             (std::fabs(massOmega - o2::constants::physics::MassOmegaMinus) > omegarej) &&
             (xiproperlifetime < properlifetimefactor * ctauxi) &&
             (std::fabs(yXi) < rapidity);
      isXiYN = (std::fabs(bachTrack.tpcNSigmaPi()) < nsigmatpcpi) &&
               (Casccascradius > lowerradiusXiYN) &&
               (std::fabs(massXi - o2::constants::physics::MassXiMinus) < deltaMassXi) &&
               (std::fabs(massOmega - o2::constants::physics::MassOmegaMinus) > omegarej) &&
               (xiproperlifetime < properlifetimefactor * ctauxi) &&
               (std::fabs(yXi) < rapidity);
      isOmega = (std::fabs(bachTrack.tpcNSigmaKa()) < nsigmatpcka) &&
                (cascCPA > casccospaomega) &&
                (DCAV0ToPV > dcav0topv) &&
                (std::fabs(massOmega - o2::constants::physics::MassOmegaMinus) < deltaMassOmega) &&
                (std::fabs(massXi - o2::constants::physics::MassXiMinus) > xirej) &&
                (Casccascradius < upperradiusOmega) &&
                (omegaproperlifetime < properlifetimefactor * ctauomega) &&
                (std::fabs(yOmega) < rapidity);
      isOmegalargeR = (std::fabs(bachTrack.tpcNSigmaKa()) < nsigmatpcka) &&
                      (cascCPA > casccospaomega) &&
                      (DCAV0ToPV > dcav0topv) &&
                      (Casccascradius > lowerradiusOmega) &&
                      (std::fabs(massOmega - o2::constants::physics::MassOmegaMinus) < deltaMassOmega) &&
                      (std::fabs(massXi - o2::constants::physics::MassXiMinus) > xirej) &&
                      (omegaproperlifetime < properlifetimefactor * ctauomega) &&
                      (std::fabs(yOmega) < rapidity);

      if (isXi) {
        QAHistos.fill(HIST("hMassXiAfterSelvsPt"), massXi, ptCasc);
        QAHistos.fill(HIST("hPtXi"), ptCasc);
        QAHistos.fill(HIST("hEtaXi"), etaCasc);
        QAHistosTopologicalVariables.fill(HIST("hProperLifetimeXi"), xiproperlifetime);
        QAHistosTopologicalVariables.fill(HIST("hCascCosPAXi"), cascCPA);
        QAHistosTopologicalVariables.fill(HIST("hV0CosPAXi"), v0DauCPA);
        QAHistosTopologicalVariables.fill(HIST("hCascRadiusXi"), Casccascradius);
        QAHistosTopologicalVariables.fill(HIST("hV0RadiusXi"), Cascv0radius);
        QAHistosTopologicalVariables.fill(HIST("hDCAV0ToPVXi"), DCAV0ToPV);
        QAHistosTopologicalVariables.fill(HIST("hDCAV0DaughtersXi"), mStraHelper.cascade.v0DaughterDCA);
        QAHistosTopologicalVariables.fill(HIST("hDCACascDaughtersXi"), mStraHelper.cascade.cascadeDaughterDCA);
        QAHistosTopologicalVariables.fill(HIST("hDCABachToPVXi"), std::fabs(mStraHelper.cascade.bachelorDCAxy));
        QAHistosTopologicalVariables.fill(HIST("hDCAPosToPVXi"), std::fabs(mStraHelper.cascade.positiveDCAxy));
        QAHistosTopologicalVariables.fill(HIST("hDCANegToPVXi"), std::fabs(mStraHelper.cascade.negativeDCAxy));
        QAHistosTopologicalVariables.fill(HIST("hInvMassLambdaXi"), LambdaMass);

        if (doextraQA) {
          QAHistos.fill(HIST("hHasTOFBachPi"), bachTrack.hasTOF(), bachTrack.pt());
          // QA PID
          if (mStraHelper.cascade.charge > 0) {
            QAHistos.fill(HIST("hTPCNsigmaXiBachPiPlus"), bachTrack.tpcNSigmaPi(), bachTrack.tpcInnerParam());
            QAHistos.fill(HIST("hTPCNsigmaXiV0PiPlus"), posTrack.tpcNSigmaPi(), posTrack.tpcInnerParam());
            QAHistos.fill(HIST("hTPCNsigmaXiV0AntiProton"), negTrack.tpcNSigmaPr(), negTrack.tpcInnerParam());
            QAHistos.fill(HIST("hHasTOFPi"), posTrack.hasTOF(), posTrack.pt());
            QAHistos.fill(HIST("hHasTOFPr"), negTrack.hasTOF(), negTrack.pt());
          } else {
            QAHistos.fill(HIST("hTPCNsigmaXiBachPiMinus"), bachTrack.tpcNSigmaPi(), bachTrack.tpcInnerParam());
            QAHistos.fill(HIST("hTPCNsigmaXiV0Proton"), posTrack.tpcNSigmaPr(), posTrack.tpcInnerParam());
            QAHistos.fill(HIST("hTPCNsigmaXiV0PiMinus"), negTrack.tpcNSigmaPi(), negTrack.tpcInnerParam());
            QAHistos.fill(HIST("hHasTOFPr"), posTrack.hasTOF(), posTrack.pt());
            QAHistos.fill(HIST("hHasTOFPi"), negTrack.hasTOF(), negTrack.pt());
          }
          QAHistos.fill(HIST("hRapXi"), yXi);
        }

        // Count number of Xi candidates
        xicounter++;
        //        v0sFromXiID.push_back({casc.posTrackId(), casc.negTrackId()});
        v0sFromXiID.push_back({posTrack.globalIndex(), negTrack.globalIndex()});
      }

      if (isXiYN) {
        // Xis for YN interactions
        xicounterYN++;
        QAHistosTopologicalVariables.fill(HIST("hCascRadiusXiYN"), Casccascradius);
      }
      if (isOmega) {
        QAHistos.fill(HIST("hMassOmegaAfterSelvsPt"), massOmega, ptCasc);
        QAHistos.fill(HIST("hPtOmega"), ptCasc);
        QAHistos.fill(HIST("hEtaOmega"), etaCasc);
        QAHistosTopologicalVariables.fill(HIST("hProperLifetimeOmega"), omegaproperlifetime);
        QAHistosTopologicalVariables.fill(HIST("hCascCosPAOmega"), cascCPA);
        QAHistosTopologicalVariables.fill(HIST("hV0CosPAOmega"), v0DauCPA);
        QAHistosTopologicalVariables.fill(HIST("hCascRadiusOmega"), Casccascradius);
        QAHistosTopologicalVariables.fill(HIST("hV0RadiusOmega"), Cascv0radius);
        QAHistosTopologicalVariables.fill(HIST("hDCAV0ToPVOmega"), DCAV0ToPV);
        QAHistosTopologicalVariables.fill(HIST("hDCAV0DaughtersOmega"), mStraHelper.cascade.v0DaughterDCA);
        QAHistosTopologicalVariables.fill(HIST("hDCACascDaughtersOmega"), mStraHelper.cascade.cascadeDaughterDCA);
        QAHistosTopologicalVariables.fill(HIST("hDCABachToPVOmega"), std::fabs(mStraHelper.cascade.bachelorDCAxy));
        QAHistosTopologicalVariables.fill(HIST("hDCAPosToPVOmega"), std::fabs(mStraHelper.cascade.positiveDCAxy));
        QAHistosTopologicalVariables.fill(HIST("hDCANegToPVOmega"), std::fabs(mStraHelper.cascade.negativeDCAxy));
        QAHistosTopologicalVariables.fill(HIST("hInvMassLambdaOmega"), LambdaMass);

        if (doextraQA) {

          // QA PID
          if (mStraHelper.cascade.charge > 0) {
            QAHistos.fill(HIST("hTPCNsigmaOmegaBachKaPlus"), bachTrack.tpcNSigmaKa(), bachTrack.tpcInnerParam());
            QAHistos.fill(HIST("hTPCNsigmaOmegaV0PiPlus"), posTrack.tpcNSigmaPi(), posTrack.tpcInnerParam());
            QAHistos.fill(HIST("hTPCNsigmaOmegaV0AntiProton"), negTrack.tpcNSigmaPr(), negTrack.tpcInnerParam());
            QAHistos.fill(HIST("hHasTOFPi"), posTrack.hasTOF(), posTrack.pt());
            QAHistos.fill(HIST("hHasTOFPr"), negTrack.hasTOF(), negTrack.pt());
          } else {
            QAHistos.fill(HIST("hTPCNsigmaOmegaBachKaMinus"), bachTrack.tpcNSigmaKa(), bachTrack.tpcInnerParam());
            QAHistos.fill(HIST("hTPCNsigmaOmegaV0Proton"), posTrack.tpcNSigmaPr(), posTrack.tpcInnerParam());
            QAHistos.fill(HIST("hTPCNsigmaOmegaV0PiMinus"), negTrack.tpcNSigmaPi(), negTrack.tpcInnerParam());
            QAHistos.fill(HIST("hHasTOFPr"), posTrack.hasTOF(), posTrack.pt());
            QAHistos.fill(HIST("hHasTOFPi"), negTrack.hasTOF(), negTrack.pt());
          }
          QAHistos.fill(HIST("hHasTOFBachKa"), bachTrack.hasTOF(), bachTrack.pt());
          QAHistos.fill(HIST("hRapOmega"), yOmega);
        }

        // Count number of Omega candidates
        omegacounter++;
        v0sFromOmegaID.push_back({posTrack.globalIndex(), negTrack.globalIndex()});
      }

      if (isOmegalargeR) {
        omegalargeRcounter++;
        QAHistosTopologicalVariables.fill(HIST("hCascRadiusOmegaLargeR"), Casccascradius);
      }
    } // end loop over cascades
    // Omega trigger definition
    keepEvent[0] = omegacounter > 0;

    std::array<bool, 11> EvtwhMinPt{false};
    std::array<float, 11> ThrdPt;
    for (int i = 0; i < 11; i++) {
      ThrdPt[i] = static_cast<float>(i);
    }

    // QA tracks
    int triggcounterAllEv = 0;
    for (auto track : tracks) { // start loop over tracks
      if (cfgTrackCuts.isTrackFilter && !selectTrack(track)) {
        continue;
      }
      triggcounterAllEv++;
      QAHistosTriggerParticles.fill(HIST("hPtTriggerAllEv"), track.pt());
      QAHistosTriggerParticles.fill(HIST("hPhiTriggerAllEv"), track.phi(), track.pt());
      QAHistosTriggerParticles.fill(HIST("hEtaTriggerAllEv"), track.eta(), track.pt());
      for (size_t i = 0; i < ThrdPt.size(); i++) {
        EvtwhMinPt[i] = track.pt() > ThrdPt[i];
      }

      // High-pT hadron + Omega trigger definition
      if (omegacounter > 0) {
        keepEvent[1] = true;
        QAHistosTriggerParticles.fill(HIST("hPtTriggerSelEv"), track.pt());
      }
    } // end loop over tracks
    for (int i = 0; i < 11; i++) {
      if (EvtwhMinPt[i]) {
        EventsvsMultiplicity.fill(HIST("hadEventsvsMultiplicityZeqV0AvsPt"), collision.multZeqFV0A(), i + 0.5);
        EventsvsMultiplicity.fill(HIST("hadEventsvsMultiplicityZeqT0MvsPt"), collision.multZeqFT0A() + collision.multZeqFT0C(), i + 0.5);
        EventsvsMultiplicity.fill(HIST("hadEventsvsMultiplicityZeqNTracksPVvsPt"), collision.multZeqNTracksPV(), i + 0.5);
      }
    }
    if (triggcounterAllEv > 0) {
      hProcessedEvents->Fill(1.5);
      if (doextraQA) {
        EventsvsMultiplicity.fill(HIST("hadEventsvsMultiplicityZeqV0A"), collision.multZeqFV0A());
        EventsvsMultiplicity.fill(HIST("hadEventsvsMultiplicityZeqT0M"), collision.multZeqFT0A() + collision.multZeqFT0C());
        EventsvsMultiplicity.fill(HIST("hadEventsvsMultiplicityZeqNTracksPV"), collision.multZeqNTracksPV());
      }
    }
    QAHistosTriggerParticles.fill(HIST("hTriggeredParticlesAllEv"), triggcounterAllEv);
    if (keepEvent[1]) {
      QAHistosTriggerParticles.fill(HIST("hTriggeredParticlesSelEv"), triggcounterAllEv);
      for (size_t i = 0; i < EvtwhMinPt.size(); i++) {
        if (EvtwhMinPt[i]) {
          hEvtvshMinPt->Fill(i + 0.5);
        }
      }
    }

    // Double/triple/quad Xi trigger definition
    if (v0sFromXiID.size() > 0) {
      std::set<std::array<int64_t, 2>> uniqueXis = {v0sFromXiID.begin(), v0sFromXiID.end()};
      if (uniqueXis.size() > 1) {
        keepEvent[2] = true;
      }
      if (uniqueXis.size() > 2) {
        keepEvent[3] = true;
      }
      if (uniqueXis.size() > 3) {
        keepEvent[4] = true;
      }
    }

    // Double Omega trigger definition
    if (v0sFromOmegaID.size() > 0) {
      std::set<std::array<int64_t, 2>> uniqueOmegas = {v0sFromOmegaID.begin(), v0sFromOmegaID.end()};
      if (uniqueOmegas.size() > 1) {
        keepEvent[10] = true;
      }
    }

    // Omega + Xi trigger definition
    if (v0sFromOmegaID.size() > 0 && v0sFromXiID.size() > 0) {
      std::set<std::array<int64_t, 2>> uniqueOmegas = {v0sFromOmegaID.begin(), v0sFromOmegaID.end()};
      std::set<std::array<int64_t, 2>> uniqueXis = {v0sFromXiID.begin(), v0sFromXiID.end()};
      if (uniqueOmegas.size() > 1 || uniqueXis.size() > 1) {
        keepEvent[11] = true;
      } else {
        // keep only if there is at least one non-overlapping v0
        for (auto v0Omega : uniqueOmegas) {
          if (uniqueXis.find(v0Omega) == uniqueXis.end()) {
            keepEvent[11] = true;
            break;
          }
        }
      }
    }

    // Single-Xi (YN) trigger definition
    if (xicounterYN > 0) {
      keepEvent[5] = true;
    }

    // Omega with high radius trigger definition
    if (omegalargeRcounter > 0) {
      keepEvent[6] = true;
    }

    // Omega in high multiplicity events
    if (omegacounter > 0 && isHighMultEventOmegaCut) {
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityFT0MwOmega"), collision.multFT0M());
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityFT0MwOmegaNorm"), multFT0MNorm);
    }
    if (omegacounter > 0 && isHighMultEventTrkOmegaCut) {
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityTrackswOmega"), multTrack);
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityTracksGlobwOmega"), collision.multNTracksGlobal());
    }
    if (omegacounter > 0 && (isHighMultEventOmegaCut || isHighMultEventTrkOmegaCut)) { // to compute "OR" selectivity
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityFT0MTrackswOmega"), multTrack);
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityFT0MTracksGlobwOmega"), collision.multNTracksGlobal());
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityFT0MTrackswOmega2D"), multTrack, multFT0MNorm);
      EventsvsMultiplicity.fill(HIST("AllEventsvsMultiplicityFT0MTracksGlobwOmega2D"), collision.multNTracksGlobal(), multFT0MNorm);
    }
    if (omegacounter > 0 && isHighMultEventOmegaCut) {
      keepEvent[9] = true;
    }
    if (omegacounter > 0 && isHighMultEventTrkOmegaCut) {
      keepEvent[13] = true;
    }
    if (isHighMultEvent) { // Normalisation tail
      keepEvent[14] = true;
    }
    if (isHighMultEventTrk) { // Normalisation tail
      keepEvent[15] = true;
    }
    for (const auto& trackedCascade : trackedCascades) {
      const auto trackCasc = trackedCascade.track_as<TrackCandidates>();
      QAHistosStrangenessTracking.fill(HIST("hPtCascTracked"), trackCasc.pt());
      QAHistosStrangenessTracking.fill(HIST("hStRVsPtTrkCasc"), trackCasc.pt(), RecoDecay::sqrtSumOfSquares(trackCasc.x(), trackCasc.y()));
      QAHistosStrangenessTracking.fill(HIST("hMatchChi2TrkCasc"), trackedCascade.matchingChi2());

      auto trackParCovTrk = getTrackParCov(trackCasc);
      o2::base::Propagator::Instance()->propagateToDCA(primaryVertex, trackParCovTrk, mBz, 2.f, o2::base::Propagator::MatCorrType::USEMatCorrLUT, &impactParameterTrk);

      QAHistosStrangenessTracking.fill(HIST("hDcaXY"), impactParameterTrk.getY());
      QAHistosStrangenessTracking.fill(HIST("hDcaXYVsPt"), trackParCovTrk.getPt(), impactParameterTrk.getY());
      QAHistosStrangenessTracking.fill(HIST("hDcaZ"), impactParameterTrk.getZ());
      QAHistosStrangenessTracking.fill(HIST("hDcaZVsPt"), trackParCovTrk.getPt(), impactParameterTrk.getZ());
      QAHistosStrangenessTracking.fill(HIST("hDcaVsPt"), impactParameterTrk.getY(), trackCasc.pt());
      QAHistosStrangenessTracking.fill(HIST("hDcaVsR"), impactParameterTrk.getY(), RecoDecay::sqrtSumOfSquares(trackCasc.x(), trackCasc.y()));
      const auto decayRadius = RecoDecay::sqrtSumOfSquares(trackedCascade.decayX(), trackedCascade.decayY());
      QAHistosStrangenessTracking.fill(HIST("hDecayRadius"), decayRadius);

      // const auto itsTrack = trackedCascade.itsTrack();
      const auto cascade = trackedCascade.cascade();
      const auto bachelor = cascade.bachelor_as<TrackCandidates>();
      const auto v0 = cascade.v0_as<o2::aod::V0s>();
      const auto negTrack = v0.negTrack_as<TrackCandidates>();
      const auto posTrack = v0.posTrack_as<TrackCandidates>();

      if (!posTrack.hasTPC() || !negTrack.hasTPC() || !bachelor.hasTPC() ||
          posTrack.tpcNClsFindable() < minNoClsTrackedCascade ||
          negTrack.tpcNClsFindable() < minNoClsTrackedCascade ||
          bachelor.tpcNClsFindable() < minNoClsTrackedCascade) {
        continue;
      }

      std::array<double, 2> masses{o2::constants::physics::MassProton, o2::constants::physics::MassPiMinus};
      std::array<std::array<float, 3>, 2> momenta;
      std::array<double, 2> nsigma;
      if (trackCasc.sign() < 0) {
        // Omega-, Xi-
        momenta[0] = {posTrack.px(), posTrack.py(), posTrack.pz()};
        momenta[1] = {negTrack.px(), negTrack.py(), negTrack.pz()};
        nsigma[0] = posTrack.tpcNSigmaPr();
        nsigma[1] = negTrack.tpcNSigmaPi();
      } else {
        // Omega+, Xi+
        momenta[0] = {negTrack.px(), negTrack.py(), negTrack.pz()};
        momenta[1] = {posTrack.px(), posTrack.py(), posTrack.pz()};
        nsigma[0] = negTrack.tpcNSigmaPr();
        nsigma[1] = posTrack.tpcNSigmaPi();
      }

      const auto v0mass = RecoDecay::m(momenta, masses);
      QAHistosStrangenessTracking.fill(HIST("hMassV0TrkCasc"), v0mass);
      QAHistosStrangenessTracking.fill(HIST("hNSigmaTpcPrTrkCascV0"), nsigma[0]);
      QAHistosStrangenessTracking.fill(HIST("hNSigmaTpcPiTrkCascV0"), nsigma[1]);
      QAHistosStrangenessTracking.fill(HIST("hNSigmaTpcPiTrkCascBachelor"), bachelor.tpcNSigmaPi());
      QAHistosStrangenessTracking.fill(HIST("hNSigmaTpcKaTrkCascBachelor"), bachelor.tpcNSigmaKa());

      // track propagation
      o2::track::TrackParCov trackParCovV0;
      o2::track::TrackPar trackParV0;
      o2::track::TrackPar trackParBachelor;
      float cpa = -1;
      if (mDCAFitter.process(getTrackParCov(negTrack), getTrackParCov(posTrack))) {
        trackParCovV0 = mDCAFitter.createParentTrackParCov(0);
        if (mDCAFitter.process(trackParCovV0, getTrackParCov(bachelor))) {
          trackParV0 = mDCAFitter.getTrackParamAtPCA(0);
          trackParBachelor = mDCAFitter.getTrackParamAtPCA(1);
          trackParV0.getPxPyPzGlo(momenta[0]);
          trackParBachelor.getPxPyPzGlo(momenta[1]);
          std::array<float, 3> pVec;
          mDCAFitter.createParentTrackParCov().getPxPyPzGlo(pVec);
          cpa = RecoDecay::cpa(pvPos, mDCAFitter.getPCACandidate(), pVec);
          QAHistosStrangenessTracking.fill(HIST("hCpa"), cpa);
        } else {
          continue;
        }
      } else {
        continue;
      }

      // Omega hypothesis
      masses = {o2::constants::physics::MassLambda0, o2::constants::physics::MassKPlus};
      const auto massOmega = recalculateMasses ? RecoDecay::m(momenta, masses) : trackedCascade.omegaMass();
      QAHistosStrangenessTracking.fill(HIST("hMassOmegaTrkCasc"), massOmega);
      QAHistosStrangenessTracking.fill(HIST("hMassOmegaVsMatchChi2TrkCasc"), massOmega, trackedCascade.matchingChi2());
      QAHistosStrangenessTracking.fill(HIST("hMassOmegaVsTopChi2TrkCasc"), massOmega, trackedCascade.topologyChi2());
      QAHistosStrangenessTracking.fill(HIST("hPtVsMassTrkOmega"), trackCasc.pt(), massOmega);

      // Xi hypothesis
      masses = {o2::constants::physics::MassLambda0, o2::constants::physics::MassPiPlus};
      const auto massXi = recalculateMasses ? RecoDecay::m(momenta, masses) : trackedCascade.xiMass();
      QAHistosStrangenessTracking.fill(HIST("hMassXiTrkCasc"), massXi);
      QAHistosStrangenessTracking.fill(HIST("hPtVsMassTrkXi"), trackCasc.pt(), massXi);
      QAHistosStrangenessTracking.fill(HIST("hMassXiVsMatchChi2TrkCasc"), massXi, trackedCascade.matchingChi2());
      QAHistosStrangenessTracking.fill(HIST("hMassXiVsTopChi2TrkCasc"), massXi, trackedCascade.topologyChi2());

      QAHistosStrangenessTracking.fill(HIST("hDecayRadiusVsXiMass"), decayRadius, massXi);
      QAHistosStrangenessTracking.fill(HIST("hDecayRadiusVsOmegaMass"), decayRadius, massOmega);

      if ((trackCasc.pt() > minPtTrackedCascade) &&
          // (trackedCascade.matchingChi2() < maxMatchingChi2TrackedCascade) &&
          (std::abs(v0mass - o2::constants::physics::MassLambda0) < massWindowLambda) &&
          (std::abs(nsigma[0]) < maxNSigmaV0PrTrackedCascade) &&
          (std::abs(nsigma[1]) < maxNSigmaV0PiTrackedCascade)) {
        // Xi
        const auto deltaMassTrackedXi = useNsigmaCutTrackedXi ? getMassWindow(stfilter::species::Xi, trackCasc.pt(), massWindowTrackedXiNsigma) : massWindowTrackedXi;
        if ((std::abs(massXi - o2::constants::physics::MassXiMinus) < deltaMassTrackedXi) &&
            (std::abs(impactParameterTrk.getY()) >= minDcaTrackedXi) &&
            (cpa <= maxCpaTrackedOmega) &&
            (std::abs(bachelor.tpcNSigmaPi()) < maxNSigmaBachelorTrackedXi)) {
          keepEvent[7] = true;
          QAHistosStrangenessTracking.fill(HIST("hDcaXYSelectedXi"), impactParameterTrk.getY());
          QAHistosStrangenessTracking.fill(HIST("hDcaZSelectedXi"), impactParameterTrk.getZ());
          QAHistosStrangenessTracking.fill(HIST("hPtVsMassTrkXiSelected"), trackCasc.pt(), massXi);
          QAHistosStrangenessTracking.fill(HIST("hNSigmaTpcPrTrkCascV0SelectedXi"), nsigma[0]);
          QAHistosStrangenessTracking.fill(HIST("hNSigmaTpcPiTrkCascV0SelectedXi"), nsigma[1]);
          QAHistosStrangenessTracking.fill(HIST("hNSigmaTpcPiTrkCascBachelorSelectedXi"), bachelor.tpcNSigmaPi());
          QAHistosStrangenessTracking.fill(HIST("hMatchChi2TrkCascSelectedXi"), trackedCascade.matchingChi2());
          QAHistosStrangenessTracking.fill(HIST("hCpaSelectedXi"), cpa);
          QAHistosStrangenessTracking.fill(HIST("hDecayRadiusSelectedXi"), decayRadius);
          QAHistosStrangenessTracking.fill(HIST("hDecayRadiusVsXiMassSelected"), decayRadius, massXi);
          QAHistosStrangenessTracking.fill(HIST("hDcaXYVsPtSelectedXi"), trackParCovTrk.getPt(), impactParameterTrk.getY());
          QAHistosStrangenessTracking.fill(HIST("hDcaZVsPtSelectedXi"), trackParCovTrk.getPt(), impactParameterTrk.getZ());
        }
        // Omega
        const auto deltaMassTrackedOmega = useNsigmaCutTrackedOmega ? getMassWindow(stfilter::species::Omega, trackCasc.pt(), massWindowTrackedOmegaNsigma) : massWindowTrackedOmega;
        if ((std::abs(massOmega - o2::constants::physics::MassOmegaMinus) < deltaMassTrackedOmega) &&
            (std::abs(massXi - o2::constants::physics::MassXiMinus) >= massWindowXiExclTrackedOmega) &&
            (std::abs(impactParameterTrk.getY()) >= minDcaTrackedOmega) &&
            (cpa <= maxCpaTrackedOmega) &&
            (std::abs(bachelor.tpcNSigmaKa()) < maxNSigmaBachelorTrackedOmega)) {
          keepEvent[8] = true;
          QAHistosStrangenessTracking.fill(HIST("hDcaXYSelectedOmega"), impactParameterTrk.getY());
          QAHistosStrangenessTracking.fill(HIST("hDcaZSelectedOmega"), impactParameterTrk.getZ());
          QAHistosStrangenessTracking.fill(HIST("hNSigmaTpcPrTrkCascV0SelectedOmega"), nsigma[0]);
          QAHistosStrangenessTracking.fill(HIST("hNSigmaTpcPiTrkCascV0SelectedOmega"), nsigma[1]);
          QAHistosStrangenessTracking.fill(HIST("hNSigmaTpcKaTrkCascBachelorSelectedOmega"), bachelor.tpcNSigmaKa());
          QAHistosStrangenessTracking.fill(HIST("hPtVsMassTrkOmegaSelected"), trackCasc.pt(), massOmega);
          QAHistosStrangenessTracking.fill(HIST("hMatchChi2TrkCascSelectedOmega"), trackedCascade.matchingChi2());
          QAHistosStrangenessTracking.fill(HIST("hCpaSelectedOmega"), cpa);
          QAHistosStrangenessTracking.fill(HIST("hDecayRadiusSelectedOmega"), decayRadius);
          QAHistosStrangenessTracking.fill(HIST("hDecayRadiusVsOmegaMassSelected"), decayRadius, massOmega);
          QAHistosStrangenessTracking.fill(HIST("hDcaXYVsPtSelectedOmega"), trackParCovTrk.getPt(), impactParameterTrk.getY());
          QAHistosStrangenessTracking.fill(HIST("hDcaZVsPtSelectedOmega"), trackParCovTrk.getPt(), impactParameterTrk.getZ());
        }
      }
    }
    // // Sigma - proton trigger definition
    for (const auto& kinkCand : kinkCands) {
      auto dauTrack = kinkCand.trackDaug_as<TrackCandidates>();
      if (!dauTrack.hasTPC() || !dauTrack.hasTOF()) {
        continue;
      }
      if (std::abs(dauTrack.tpcNSigmaPr()) > cfgSigma.nsigmatpcSigma || std::abs(dauTrack.tofNSigmaPr()) > cfgSigma.nsigmatofSigma) {
        continue;
      }
      if (kinkCand.ptMoth() < cfgSigma.minPtSigma) {
        continue;
      }
      if (kinkCand.mSigmaPlus() < cfgSigma.minMassSigma || kinkCand.mSigmaPlus() > cfgSigma.maxMassSigma) {
        continue;
      }
      QAHistosSigma.fill(HIST("hPtVsMassSigmaPlus"), kinkCand.mothSign() * kinkCand.ptMoth(), kinkCand.mSigmaPlus());
      std::array<float, 3> momMoth = {kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth()};
      std::array<float, 3> momDaug = {kinkCand.pxDaug(), kinkCand.pyDaug(), kinkCand.pzDaug()};
      std::array<float, 3> primaryVtx = {collision.posX(), collision.posY(), collision.posZ()};
      std::array<float, 3> decayVtx = {kinkCand.xDecVtx(), kinkCand.yDecVtx(), kinkCand.zDecVtx()};
      float qtAP = getQtAP(momMoth, momDaug);
      float alphaAP = getAlphaAP(momMoth, momDaug);
      float cosPA = getCosPA(momMoth, decayVtx, primaryVtx);
      if (alphaAP < 0) {
        continue;
      }
      if (qtAP < cfgSigma.minQtAPSigma || qtAP > cfgSigma.maxQtAPSigma) {
        continue;
      }
      if (cosPA < cfgSigma.minCosPASigma) {
        continue;
      }
      if (std::abs(kinkCand.dcaMothPv()) > cfgSigma.maxDCAtoPVSigma) {
        continue;
      }
      float decRad = std::hypot(kinkCand.xDecVtx(), kinkCand.yDecVtx());
      if (decRad < cfgSigma.minRadiusSigma) {
        continue;
      }
      QAHistosSigma.fill(HIST("hDecayRadiusSigma"), decRad);
      // pair a proton
      bool isProtonPaired = false;
      for (auto track : tracks) {
        if (track.globalIndex() == dauTrack.globalIndex()) {
          continue;
        }
        if (std::abs(track.tpcNSigmaPr()) > cfgSigma.nsigmatpcSigma) {
          continue;
        }
        QAHistosSigma.fill(HIST("hPtNSigmaTPCPrPair"), track.sign() * track.pt(), track.tpcNSigmaPr());
        TLorentzVector sigmaVec, protonVec;
        sigmaVec.SetXYZM(kinkCand.pxMoth(), kinkCand.pyMoth(), kinkCand.pzMoth(), o2::constants::physics::MassSigmaPlus);
        protonVec.SetXYZM(track.px(), track.py(), track.pz(), o2::constants::physics::MassProton);
        float kstar = getKStar(sigmaVec, protonVec);
        if (kstar > cfgSigma.maxKStarSigmaProton) {
          continue;
        }
        QAHistosSigma.fill(HIST("hKStarSigmaPr"), kstar);
        if (track.pt() < cfgSigma.minPtProtonTOF) {
          isProtonPaired = true;
          break;
        }
        if (!track.hasTOF()) {
          continue;
        }
        QAHistosSigma.fill(HIST("hPtNSigmaTOFPrPair"), track.sign() * track.pt(), track.tofNSigmaPr());
        if (std::abs(track.tofNSigmaPr()) > cfgSigma.nsigmatofSigma) {
          continue;
        }
        isProtonPaired = true;
        break;
      }
      if (isProtonPaired) {
        keepEvent[16] = true;
      }
    }

    // Fill centrality dependent histos
    if (keepEvent[0]) {
      hProcessedEvents->Fill(2.5);
    }
    if (keepEvent[1]) {
      hProcessedEvents->Fill(3.5);
    }
    if (keepEvent[2]) {
      hProcessedEvents->Fill(4.5);
    }
    if (keepEvent[3]) {
      hProcessedEvents->Fill(5.5);
    }
    if (keepEvent[4]) {
      hProcessedEvents->Fill(6.5);
    }
    if (keepEvent[5]) {
      hProcessedEvents->Fill(7.5);
    }
    if (keepEvent[6]) {
      hProcessedEvents->Fill(8.5);
    }
    if (xicounter > 0) {
      hProcessedEvents->Fill(9.5);
    }
    if (keepEvent[7]) {
      hProcessedEvents->Fill(10.5);
    }
    if (keepEvent[8]) {
      hProcessedEvents->Fill(11.5);
    }
    if (keepEvent[9]) {
      hProcessedEvents->Fill(12.5);
    }
    if (keepEvent[10]) {
      hProcessedEvents->Fill(13.5);
    }
    if (keepEvent[11]) {
      hProcessedEvents->Fill(14.5);
    }
    if (keepEvent[12]) {
      hProcessedEvents->Fill(15.5);
    }
    if (keepEvent[13]) {
      hProcessedEvents->Fill(16.5);
    }
    if (keepEvent[14]) {
      hProcessedEvents->Fill(17.5);
    }
    if (keepEvent[15]) {
      hProcessedEvents->Fill(18.5);
    }
    if (keepEvent[16]) {
      hProcessedEvents->Fill(19.5);
    }
    // Filling the table
    fillTriggerTable(keepEvent);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<strangenessFilter>(cfgc, TaskName{"lf-strangeness-filter"})};
}
