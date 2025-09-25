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

/// \file kaonIsospinFluctuations.cxx
/// \brief Kaon Isospin fluctuations
///
/// \author Rahul Verma (rahul.verma@iitb.ac.in) :: Sadhana Dash (sadhana@phy.iitb.ac.in)

#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Tools/TListHandler.h"

#include "CCDB/BasicCCDBManager.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include <Math/LorentzVector.h>
#include <Math/PxPyPzM4D.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <string>
#include <typeindex>
#include <vector>

// --------------------------------------------------------------------------------------------------
// Prim Vtx Table Creator to check contamination from prim vtx.
// --------------------------------------------------------------------------------------------------

namespace o2::aod
{
using Collision = aod::Collisions::iterator;
namespace primvtxcndt
{
DECLARE_SOA_INDEX_COLUMN_FULL(Collision, collision, int, Collisions, ""); //! primVtxTable to collision
DECLARE_SOA_INDEX_COLUMN_FULL(PosTrack, posTrack, int, Tracks, "_Pos");   //! Positive track
DECLARE_SOA_INDEX_COLUMN_FULL(NegTrack, negTrack, int, Tracks, "_Neg");   //! Negative track

DECLARE_SOA_COLUMN(Pt, pt, float);
DECLARE_SOA_COLUMN(Eta, eta, float);
DECLARE_SOA_COLUMN(Phi, phi, float);

DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);

DECLARE_SOA_COLUMN(MPhi1020, mPhi1020, float);
DECLARE_SOA_COLUMN(MJPsiToEE, mJPsiToEE, float);
DECLARE_SOA_COLUMN(MJPsiToMuMu, mJPsiToMuMu, float);
DECLARE_SOA_COLUMN(MKStar892, mKStar892, float);
DECLARE_SOA_COLUMN(MKStar892Bar, mKStar892Bar, float);
DECLARE_SOA_COLUMN(MRho770, mRho770, float);

DECLARE_SOA_COLUMN(EPhi1020, ePhi1020, float);
DECLARE_SOA_COLUMN(EJPsiToEE, eJPsiToEE, float);
DECLARE_SOA_COLUMN(EJPsiToMuMu, eJPsiToMuMu, float);
DECLARE_SOA_COLUMN(EKStar892, eKStar892, float);
DECLARE_SOA_COLUMN(EKStar892Bar, eKStar892Bar, float);
DECLARE_SOA_COLUMN(ERho770, eRho770, float);

} // namespace primvtxcndt

DECLARE_SOA_TABLE(PrimVtxCndts, "AOD", "PRIMVTXCNDTS", o2::soa::Index<>,
                  o2::aod::primvtxcndt::CollisionId,
                  o2::aod::primvtxcndt::PosTrackId,
                  o2::aod::primvtxcndt::NegTrackId,
                  o2::aod::primvtxcndt::Pt,
                  o2::aod::primvtxcndt::Eta,
                  o2::aod::primvtxcndt::Phi,
                  o2::aod::primvtxcndt::Px,
                  o2::aod::primvtxcndt::Py,
                  o2::aod::primvtxcndt::Pz,
                  o2::aod::primvtxcndt::MPhi1020,
                  o2::aod::primvtxcndt::MJPsiToEE,
                  o2::aod::primvtxcndt::MJPsiToMuMu,
                  o2::aod::primvtxcndt::MKStar892,
                  o2::aod::primvtxcndt::MKStar892Bar,
                  o2::aod::primvtxcndt::MRho770);
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics; // for constants

#define ID_BIT_PI 0 // Identificationi bits for PID checks
#define ID_BIT_KA 1
#define ID_BIT_PR 2
#define ID_BIT_EL 3
#define ID_BIT_MU 4
#define ID_BIT_DE 5

#define BIT_IS_K0S 0
#define BIT_IS_LAMBDA 1
#define BIT_IS_ANTILAMBDA 2
#define BIT_IS_GAMMA 3
#define BIT_IS_PHI_1020 4
#define BIT_IS_JPSI_TO_EE 5
#define BIT_IS_JPSI_TO_MUMU 6
#define BIT_IS_RHO_770 7
#define BIT_IS_KSTAR_892 8
#define BIT_IS_KSTAR_892_BAR 9
#define BIT_IS_LAMBDA_1520 10
#define BIT_IS_DELTA_1232 11

#define BIT_POS_DAU_HAS_SAME_COLL 0
#define BIT_NEG_DAU_HAS_SAME_COLL 1
#define BIT_BOTH_DAU_HAS_SAME_COLL 2

#define BIT_POS_DAU_HAS_SAME_COLL_BC 0
#define BIT_NEG_DAU_HAS_SAME_COLL_BC 1
#define BIT_BOTH_DAU_HAS_SAME_COLL_BC 2

#define BITSET(mask, ithBit) ((mask) |= (1 << (ithBit)))    // avoid name bitset as std::bitset is already there
#define BITCHECK(mask, ithBit) ((mask) & (1 << (ithBit)))   // bit check will return int value of (1<<ithBit), not bool, use BITCHECK != 0 in Analysis
#define OTHERBITS(mask, ithBit) ((mask) & ~(1 << (ithBit))) // Returns new value with ithBit set to zero and all other bits unmodified; does not modify the original mask

#define BOOL_BITCHECK(mask, ithBit) (((mask) & (1 << (ithBit))) != 0)     // returns true or false, bit check will return int value of (1<<ithBit)
#define BOOL_OTHERBITSON(mask, ithBit) (((mask) & ~(1 << (ithBit))) != 0) // returns true or false, checks if any of the other bit execpt ithBit is non zero or not.

// Precomputed powers of 10 up to 1e18 (19 elements) //int64_t has max 19 digits
static const double doublePowersOf10[] = {
  1e0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9,
  1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18};

static const int64_t intPowersOf10[] = {
  1LL, 10LL, 100LL, 1000LL, 10000LL, 100000LL, 1000000LL,
  10000000LL, 100000000LL, 1000000000LL, 10000000000LL,
  100000000000LL, 1000000000000LL, 10000000000000LL,
  100000000000000LL, 1000000000000000LL,
  10000000000000000LL, 100000000000000000LL,
  1000000000000000000LL};

template <typename T>
inline T getCfg(const auto& cfgLabelledArray, int32_t row, int32_t col)
{
  double val = cfgLabelledArray->get(row, col); // Always stored as double

  if constexpr (std::is_same_v<T, bool>) {
    return val != 0.0; // will return anything nonzero as true
  } else if constexpr (std::is_same_v<T, float>) {
    return static_cast<float>(val);
  } else if constexpr (std::is_same_v<T, double>) {
    return val;
  } else if constexpr (std::is_same_v<T, int>) {
    return static_cast<int>(val);
  } else {
    static_assert(
      std::is_same_v<T, bool> || std::is_same_v<T, float> || std::is_same_v<T, double> || std::is_same_v<T, int>,
      "getCfg<T>() only supports T = bool, float, double or int");
  }
}

template <typename T>
double calculateTime(T Start)
{
  auto stopTime = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> duration = stopTime - Start; // precise seconds (double)
  double timeInSeconds = duration.count();
  return timeInSeconds;
}
template <typename T>
void printTime(T Start, std::string String)
{
  LOG(info) << String << calculateTime(Start) << " seconds";
}

void printHistInfo(std::string String, const auto& hist)
{
  int64_t iEntries = hist->GetEntries();
  double mean = hist->GetMean();
  int64_t overflow = hist->GetBinContent(0);
  int64_t underflow = hist->GetBinContent(hist->GetNbinsX() + 1);
  double totalSum = iEntries * mean;
  LOG(info) << String << " :: totalSum = " << totalSum << " :: iEntries = " << iEntries << " :: mean = " << mean << " :: overflow = " << overflow << " :: underflow = " << underflow;
}

enum PrimVertexParicleCutTypeEnum {
  kPrimTrkPhi1020 = 0,
  kPrimTrkJPsiToEE,
  kPrimTrkJPsiToMuMu,
  kPrimTrkKStar892,
  kPrimTrkKStar892Bar,
  kPrimTrkRho770,
  primVtxTrkEnumSize
};

enum PrimVertexParicleCutFlagEnum {
  kPrimMassLow = 0,
  kPrimMassUp,
  kPrimDoRecoCheck,
  kPrimFillPreSel,
  kPrimFillPreSelDau,
  kPrimFillPostSel,
  kPrimFillPostSelDau,
  kPrimCheckPrimVtxContm,
  kPrimCheckV0DecayContm,
  kPrimCheckMPMSigma
};

struct PrimVtxParticleTable {
  Produces<aod::PrimVtxCndts> genPrimVtxCndts;

  HistogramRegistry primVtxQAPlots{"primVtxQAPlots", {}, OutputObjHandlingPolicy::AnalysisObject, true};

  static constexpr double DefaultWideMassCutValues[6][3] = {{1.013, 1.026, 1}, {3.035, 3.155, 1}, {3.035, 3.155, 1}, {0.875, 0.915, 1}, {0.875, 0.915, 1}, {0.650, 0.900, 1}};
  Configurable<LabeledArray<double>> cfgPrimVtxCndtsMassCuts{"cfgPrimVtxCndtsMassCuts", {&DefaultWideMassCutValues[0][0], 6, 3, {"Phi1020", "JPsiToEE", "JPsiToMuMu", "KStar892", "KStar892Bar", "Rho770"}, {"massLow", "massUp", "doRecoCheck"}}, "cut values for primary vertex reconstruction to reduce pair memory storage"};

  struct : ConfigurableGroup {
    Configurable<bool> printDebugMessages{"printDebugMessages", false, "printDebugMessages"};
    Configurable<bool> resetHistograms{"resetHistograms", false, "resetHistograms"};
  } cfgDebug;

  struct : ConfigurableGroup {
    Configurable<float> trackPtLow{"trackPtLow", 0.15, "trackPtLow"};
    Configurable<float> trackEtaUp{"trackEtaUp", 0.8, "trackEtaUp"};
    Configurable<int> trackTPCNClsCrossedRowsLow{"trackTPCNClsCrossedRowsLow", 60, "trackTPCNClsCrossedRowsLow"};
  } analysisCuts;

  void init(InitContext const&)
  {
    // Printing the Stored Registry information
    primVtxQAPlots.add("mPhi1020", "mPhi1020", kTH1F, {{200, 0.99f, 1.07f, "#it{M}_{inv}^{#phi} [GeV/#it{c}^{2}]"}});
    primVtxQAPlots.add("mJPsiToEE", "mJPsiToEE", kTH1F, {{200, 3.00f, 3.20f, "#it{M}_{inv}^{J/#psi (ee)} [GeV/#it{c}^{2}]"}});
    primVtxQAPlots.add("mJPsiToMuMu", "mJPsiToMuMu", kTH1F, {{200, 3.00f, 3.20f, "#it{M}_{inv}^{J/#psi (MuMu)} [GeV/#it{c}^{2}]"}});
    primVtxQAPlots.add("mKStar892", "mKStar892", kTH1F, {{200, 0.860f, 0.920f, "#it{M}_{inv}^{K^{*}(892)} [GeV/#it{c}^{2}]"}});
    primVtxQAPlots.add("mKStar892Bar", "mKStar892Bar", kTH1F, {{200, 0.860f, 0.920f, "#it{M}_{inv}^{K^{*}(892)Bar} [GeV/#it{c}^{2}]"}});
    primVtxQAPlots.add("mRho770", "mRho770", kTH1F, {{200, 0.600, 0.950f, "#it{M}_{inv}^{#rho(770)} [GeV/#it{c}^{2}]"}});
    if (cfgDebug.printDebugMessages) {
      LOG(info) << "Printing Stored Registry Information";
      primVtxQAPlots.print();
    }
  }

  // // Event Filter //Configurable filter to be added in the future.
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);
  // Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutZvertex);

  // // Track Filter
  Filter ptFilter = o2::aod::track::pt > analysisCuts.trackPtLow; //-->ohterFilters to be added later like && (o2::aod::track::pt) < cfgTrackCuts.cfgTrk06PtHigh;
  Filter etaFilter = (nabs(o2::aod::track::eta) < analysisCuts.trackEtaUp);

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels>;
  using MyTracks = soa::Filtered<soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;

  // For manual sliceBy
  SliceCache cache;

  Preslice<MyTracks> tracksPerCollisionPreslice = o2::aod::track::collisionId;

  Partition<MyTracks> posTracks = aod::track::signed1Pt > 0.0f;
  Partition<MyTracks> negTracks = aod::track::signed1Pt < 0.0f;

  std::chrono::high_resolution_clock::time_point start0 = std::chrono::high_resolution_clock::now();
  int64_t dfCount = 0;

  /// Reconstruction of particles from primary vertex
  /// to get contamination sources from primary vertex
  /// check the signal over background ratio of the particles. If it is very bad then remove the particle
  /// phi(1020) -> K+ + K-         :: has clean Peak and good abundance
  /// J/ψ       -> e+/e-, mu+/mu-  :: has clean Peak and may have abundance in high multiplicity collisions
  /// K(892)*   -> K + pi,         :: // Many decay modes
  /// ρ(770)    -> pi+ + pi-       :: // Is this reliable, High background and broad peak is expected
  /// Δ(1232)   -> pπ+             :: // To do :: Add this in future, presently don't know how good the peak is.
  /// Λ(1520)   ->                 :: // To do :: Add this in future
  /// Upsilon for muon checking and muon pion distinction in future.

  void processData(aod::BCsWithTimestamps const&, MyCollisions const& collisions, MyTracks const&)
  {
    dfCount++;
    auto start1 = std::chrono::high_resolution_clock::now();

    float px = -999.0;
    float py = -999.0;
    float pz = -999.0;
    float p = -999.0;
    float pt = -999.0;
    float eta = -999.0;
    float phi = -999.0;

    float mPhi1020 = -999.0;
    float mJPsiToEE = -999.0;
    float mJPsiToMuMu = -999.0;
    float mKStar892 = -999.0;
    float mKStar892Bar = -999.0;
    float mRho770 = -999.0;

    float ePosKa = -999.0;
    float eNegKa = -999.0;
    float ePosPi = -999.0;
    float eNegPi = -999.0;
    float ePosEl = -999.0;
    float eNegEl = -999.0;
    float ePosMu = -999.0;
    float eNegMu = -999.0;

    // phi(1020)
    static const float phiMassLow = getCfg<float>(cfgPrimVtxCndtsMassCuts, kPrimTrkPhi1020, kPrimMassLow);
    static const float phiMassUp = getCfg<float>(cfgPrimVtxCndtsMassCuts, kPrimTrkPhi1020, kPrimMassUp);
    static const bool doPhi1020 = getCfg<bool>(cfgPrimVtxCndtsMassCuts, kPrimTrkPhi1020, kPrimDoRecoCheck);

    // J/psi -> e+e-
    static const float jpsiEEMassLow = getCfg<float>(cfgPrimVtxCndtsMassCuts, kPrimTrkJPsiToEE, kPrimMassLow);
    static const float jpsiEEMassUp = getCfg<float>(cfgPrimVtxCndtsMassCuts, kPrimTrkJPsiToEE, kPrimMassUp);
    static const bool doJPsiToEE = getCfg<bool>(cfgPrimVtxCndtsMassCuts, kPrimTrkJPsiToEE, kPrimDoRecoCheck);

    // J/psi -> mu+mu-
    static const float jpsiMuMuMassLow = getCfg<float>(cfgPrimVtxCndtsMassCuts, kPrimTrkJPsiToMuMu, kPrimMassLow);
    static const float jpsiMuMuMassUp = getCfg<float>(cfgPrimVtxCndtsMassCuts, kPrimTrkJPsiToMuMu, kPrimMassUp);
    static const bool doJPsiToMuMu = getCfg<bool>(cfgPrimVtxCndtsMassCuts, kPrimTrkJPsiToMuMu, kPrimDoRecoCheck);

    // K*(892)
    static const float kstarMassLow = getCfg<float>(cfgPrimVtxCndtsMassCuts, kPrimTrkKStar892, kPrimMassLow);
    static const float kstarMassUp = getCfg<float>(cfgPrimVtxCndtsMassCuts, kPrimTrkKStar892, kPrimMassUp);
    static const bool doKStar892 = getCfg<bool>(cfgPrimVtxCndtsMassCuts, kPrimTrkKStar892, kPrimDoRecoCheck);

    // K*(892)Bar
    static const float kstarBarMassLow = getCfg<float>(cfgPrimVtxCndtsMassCuts, kPrimTrkKStar892Bar, kPrimMassLow);
    static const float kstarBarMassUp = getCfg<float>(cfgPrimVtxCndtsMassCuts, kPrimTrkKStar892Bar, kPrimMassUp);
    static const bool doKStar892Bar = getCfg<bool>(cfgPrimVtxCndtsMassCuts, kPrimTrkKStar892Bar, kPrimDoRecoCheck);

    // ρ(770)
    static const float rhoMassLow = getCfg<float>(cfgPrimVtxCndtsMassCuts, kPrimTrkRho770, kPrimMassLow);
    static const float rhoMassUp = getCfg<float>(cfgPrimVtxCndtsMassCuts, kPrimTrkRho770, kPrimMassUp);
    static const bool doRho770 = getCfg<bool>(cfgPrimVtxCndtsMassCuts, kPrimTrkRho770, kPrimDoRecoCheck);

    bool fillTable = false;

    auto collLoop1Start = std::chrono::high_resolution_clock::now();
    // Background Estimation From Primary Vertex
    for (const auto& collision : collisions) {
      if (!collision.sel8()) {
        continue;
      }
      auto posTracksPerColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
      auto negTracksPerColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);

      for (const auto& posTrack : posTracksPerColl) {
        if (posTrack.tpcNClsCrossedRows() < analysisCuts.trackTPCNClsCrossedRowsLow)
          continue;
        for (const auto& negTrack : negTracksPerColl) {
          if (negTrack.tpcNClsCrossedRows() < analysisCuts.trackTPCNClsCrossedRowsLow)
            continue;
          px = posTrack.px() + negTrack.px();
          py = posTrack.py() + negTrack.pz();
          pz = posTrack.py() + negTrack.pz();
          p = RecoDecay::p(px, py, pz);                           // definition = std::sqrt(px*px+py*py+pz*pz);
          pt = RecoDecay::pt(px, py);                             // definition = std::sqrt(px*px+py*py);
          eta = RecoDecay::eta(std::array<float, 3>{px, py, pz}); // definition = 0.5 * std::log((p + pz) / (p - pz + 1e-12)); // 1e-12 protect divide-by-zero
          phi = RecoDecay::phi(px, py);                           // definition = std::atan2(py, px);

          // Precompute energies for fast speed
          ePosKa = RecoDecay::e(posTrack.px(), posTrack.py(), posTrack.pz(), MassKPlus);
          eNegKa = RecoDecay::e(negTrack.px(), negTrack.py(), negTrack.pz(), MassKPlus);
          ePosPi = RecoDecay::e(posTrack.px(), posTrack.py(), posTrack.pz(), MassPiPlus);
          eNegPi = RecoDecay::e(negTrack.px(), negTrack.py(), negTrack.pz(), MassPiPlus);

          fillTable = false;
          // phi(1020) -> K+ + K-
          if (doPhi1020) {
            mPhi1020 = RecoDecay::m(p, ePosKa + eNegKa);
            if (phiMassLow < mPhi1020 && mPhi1020 < phiMassUp)
              fillTable = true;
          }

          // J/ψ       -> e+/e-
          if (doJPsiToEE) {
            ePosEl = RecoDecay::e(posTrack.px(), posTrack.py(), posTrack.pz(), MassElectron);
            eNegEl = RecoDecay::e(negTrack.px(), negTrack.py(), negTrack.pz(), MassElectron);
            mJPsiToEE = RecoDecay::m(p, ePosEl + eNegEl);
            if (jpsiEEMassLow < mJPsiToEE && mJPsiToEE < jpsiEEMassUp)
              fillTable = true;
          }

          // J/ψ       -> mu+/mu-
          if (doJPsiToMuMu) {
            ePosMu = RecoDecay::e(posTrack.px(), posTrack.py(), posTrack.pz(), MassMuonMinus);
            eNegMu = RecoDecay::e(negTrack.px(), negTrack.py(), negTrack.pz(), MassMuonMinus);
            mJPsiToMuMu = RecoDecay::m(p, ePosMu + eNegMu);
            if (jpsiMuMuMassLow < mJPsiToMuMu && mJPsiToMuMu < jpsiMuMuMassUp)
              fillTable = true;
          }

          // K(892)*   -> K+ + pi-
          if (doKStar892) {
            mKStar892 = RecoDecay::m(p, ePosKa + eNegPi);
            if (kstarMassLow < mKStar892 && mKStar892 < kstarMassUp)
              fillTable = true;
          }

          // K(892)*Bar -> K- + pi+
          if (doKStar892Bar) {
            mKStar892Bar = RecoDecay::m(p, ePosPi + eNegKa);
            if (kstarBarMassLow < mKStar892Bar && mKStar892Bar < kstarBarMassUp)
              fillTable = true;
          }
          // ρ(770)    -> pi+ + pi-
          if (doRho770) {
            mRho770 = RecoDecay::m(p, ePosPi + eNegPi);
            if (rhoMassLow < mRho770 && mRho770 < rhoMassUp)
              fillTable = true;
          }

          // Δ(1232)   -> p + π :: // To do :: Add this in future
          // Λ(1520)   ->       :: // To do :: Add this in future

          primVtxQAPlots.fill(HIST("mPhi1020"), mPhi1020);
          primVtxQAPlots.fill(HIST("mJPsiToEE"), mJPsiToEE);
          primVtxQAPlots.fill(HIST("mJPsiToMuMu"), mJPsiToMuMu);
          primVtxQAPlots.fill(HIST("mKStar892"), mKStar892);
          primVtxQAPlots.fill(HIST("mKStar892Bar"), mKStar892Bar);
          primVtxQAPlots.fill(HIST("mRho770"), mRho770);

          if (fillTable) {
            genPrimVtxCndts(collision.globalIndex(), posTrack.globalIndex(), negTrack.globalIndex(),
                            pt, eta, phi, px, py, pz,
                            mPhi1020, mJPsiToEE, mJPsiToMuMu, mKStar892, mKStar892Bar, mRho770);
          }
        } // negTrackLoop
      } // posTrackLoop
    } // collision loop

    if (cfgDebug.printDebugMessages) {
      printTime(collLoop1Start, Form("DEBUG :: df_%lld :: collLoop1 Time :: ", dfCount));
      printTime(start1, Form("DEBUG :: df_%lld :: DF Reading :: DF Processing Time :: ", dfCount));
      printTime(start0, Form("DEBUG :: df_%lld :: DF Reading :: DF Elapsed Time :: ", dfCount));
    }
  } // Process Function Ends
  PROCESS_SWITCH(PrimVtxParticleTable, processData, "Process for Data", true);
};

// Stable pdg masses and width
// Constant mass and widths to be used //==> To Do :: Find and Fix these to the correct values
constexpr std::array<float, 6> MPrimVtxCndt = {MassPhi, MassJPsi, MassJPsi, MassK0Star892, MassK0Star892, 0.77526f};
constexpr std::array<float, 4> MV0DecayCndt = {MassK0Short, MassLambda0, MassLambda0Bar, MassGamma};

constexpr std::array<float, 6> WPrimVtxCndt = {0.004249f, 0.0000929f, 0.0000929f, 0.0473f, 0.0473f, 0.1491f}; // from acutal pdg widths
constexpr std::array<float, 4> WV0DecayCndt = {0.01f, 0.01f, 0.01f, 0.01f};                                   // actual pdg widths are very small, using 0.01 for now

static const std::string defaultSparseAxis[10][3] = {
  {"centrality", "centrality", "centrality"},
  {"nK0s/kPos/iNTrk", "nK0s/kPos/iNTrk", "nKa/kPos/iNTrk"},
  {"nKa/kPos/iNTrk", "nKaon/kPos/iNTrk", "nKa/kNeg/iNTrk"},
  {"nKa/kNeg/iNTrk", "(nK0s)^{2}/kPos/iNTrk", "(nKaPlus)^{2}/kPos/iNTrk"},
  {"nK0s/kPos/effSum", "(nKaon)^{2}/kPos/iNTrk", "(nKaMinus)^{2}/kPos/iNTrk"},
  {"nKa/kPos/effSum", "(nK0s*nKaon)/kPos/iNTrk", "(nKaPlus*nKaMinus)/kPos/iNTrk"},
  {"nKa/kNeg/effSum", "nTrack/kPos/iNTrk", "nTrack/kPos/iNTrk"},
  {"", "(nKaPlus)^{2}/kPos/iNTrk", "nKaon/kPos/iNTrk"},
  {"", "(nKaMinus)^{2}/kPos/iNTrk", "(nK0s)^{2}/kPos/iNTrk"},
  {"", "(nKaPlus*nKaMinus)/kPos/iNTrk", "(nK0s*nKaon)/kPos/iNTrk"},
};

struct KaonIsospinFluctuations {
  // Hisogram registry:
  HistogramRegistry anlysisTimingInfo{"anlysisTimingInfo", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoEvent{"recoEvent", {}, OutputObjHandlingPolicy::AnalysisObject, true};

  HistogramRegistry recoV0sFull{"recoV0sFull", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoV0sPostK0sCheck{"recoV0sPostK0sCheck", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoV0sPostMassCut{"recoV0sPostMassCut", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoV0sPostK0sSelectionCut{"recoV0sPostK0sSelectionCut", {}, OutputObjHandlingPolicy::AnalysisObject, true};

  HistogramRegistry recoV0sPostLambdaCheck{"recoV0sPostLambdaCheck", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoV0sPostAntiLambdaCheck{"recoV0sPostAntiLambdaCheck", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoV0sPostGammaCheck{"recoV0sPostGammaCheck", {}, OutputObjHandlingPolicy::AnalysisObject, true};

  HistogramRegistry recoPhi1020{"recoPhi1020", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoJPsiToEE{"recoJPsiToEE", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoJPsiToMuMu{"recoJPsiToMuMu", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoKStar892{"recoKStar892", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoKStar892Bar{"recoKStar892Bar", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoRho770{"recoRho770", {}, OutputObjHandlingPolicy::AnalysisObject, true};

  HistogramRegistry recoK0s{"recoK0s", {}, OutputObjHandlingPolicy::AnalysisObject, true};

  HistogramRegistry recoTracks{"recoTracks", {}, OutputObjHandlingPolicy::AnalysisObject, true};

  HistogramRegistry recoAnalysisPi{"recoAnalysisPi", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoAnalysisKa{"recoAnalysisKa", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoAnalysisPr{"recoAnalysisPr", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoAnalysisEl{"recoAnalysisEl", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoAnalysisDe{"recoAnalysisDe", {}, OutputObjHandlingPolicy::AnalysisObject, true};

  HistogramRegistry recoRejectedPi{"recoRejectedPi", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoRejectedKa{"recoRejectedKa", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoRejectedPr{"recoRejectedPr", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoRejectedEl{"recoRejectedEl", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry recoRejectedDe{"recoRejectedDe", {}, OutputObjHandlingPolicy::AnalysisObject, true};

  HistogramRegistry recoAnalysis{"recoAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject, true};
  HistogramRegistry genAnalysis{"genAnalysis", {}, OutputObjHandlingPolicy::AnalysisObject, true};

  OutputObj<TList> fOutput{"TListForAnalysis", OutputObjHandlingPolicy::AnalysisObject, OutputObjSourceType::OutputObjSource};
  TListHandler listHandler;
  std::array<TListHandler, 6> primVtxTLHs;

  // PDG data base
  // Service<o2::framework::O2DatabasePDG> pdgDB;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Configurables
  // TListHandler config
  struct : ConfigurableGroup {
    Configurable<bool> makeNestedTList{"makeNestedTList", false, "makeNestedTList"};
    Configurable<bool> printDebugMessages{"printDebugMessages", true, "printDebugMessages"};
    Configurable<bool> resetHistograms{"resetHistograms", true, "resetHistograms"};
  } cfgDebug;

  // Event Selection
  Configurable<float> cutZvertex{"cutZvertex", 8.0f, "Accepted z-vertex range (cm)"};

  // Histogram Configurables
  struct : ConfigurableGroup {
    Configurable<int> centAxis01Bins{"centAxis01Bins", 1020, "No of bins in centrality axis"};
    Configurable<double> centAxis02XLow{"centAxis02XLow", -1.0, "centAxis02XLow"};
    Configurable<double> centAxis03XUp{"centAxis03XUp", 101.0, "centAxis03XUp"};
    Configurable<int> centAxis04Type{"centAxis04Type", 0, "centAxis04Type"}; // 0 -> centFT0C, 1-> centFT0M, 2->multFT0M, 3->multFT0C
    Configurable<std::string> centAxis05Name{"centAxis05Name", "centFT0C", "centAxis05Name"};
  } cfgCentAxis;

  // Configurable parameters for V0 selection
  Configurable<float> v0settingDcaPosToPV{"v0settingDcaPosToPV", 0.06, "DCA Pos to PV"};
  Configurable<float> v0settingDcaNegToPV{"v0settingDcaNegToPV", 0.06, "DCA Neg to PV"};
  Configurable<float> v0settingDcaV0Dau{"v0settingDcaV0Dau", 1, "DCA V0 Daughters"};
  Configurable<double> v0settingCosPA{"v0settingCosPA", 0.98, "V0 CosPA"};
  Configurable<float> v0settingRadius{"v0settingRadius", 0.5, "v0radius"};

  static constexpr double DefaulV0CutValues[4][15] = {{0.48, 0.515, 0.1, 1.5, 0.5, 0.2, 0.2, 1, 1, 1, 1, 1, 1, 1, 1}, {1.10, 1.12, -0.1, 999.0, 999.0, 999.0, 999.0, 1, 1, 1, 1, 1, 1, 1, 1}, {1.10, 1.12, -0.1, 999.0, 999.0, 999.0, 999.0, 1, 1, 1, 1, 1, 1, 1, 1}, {-0.1, 0.03, -0.1, 999.0, 999.0, 0.01, 0.2, 1, 1, 1, 1, 1, 1, 1, 1}};
  Configurable<LabeledArray<double>> cfgV0ParticleCuts = {"cfgV0ParticleCuts", {&DefaulV0CutValues[0][0], 4, 15, {"K0s", "Lambda", "AntiLambda", "Gamma"}, {"V0MLow", "V0MUp", "V0LowPt", "V0HighPt", "V0Rapitidy", "V0ARMcut1", "V0ARMcut2", "V0DoRecoCheck", "V0FillPreSel", "V0FillPreSelDau", "V0FillPostSel", "V0FillPostSelDau", "V0CheckPrimVtxContm", "V0CheckV0DecayContm", "V0CheckMPMSigma"}}, "cut values for V0 reconstruction"};

  static constexpr double DefaulPrimVtxCutValues[6][10] = {{1.013, 1.026, 1, 1, 1, 1, 1, 1, 1, 1}, {3.035, 3.155, 1, 1, 1, 1, 1, 1, 1, 1}, {3.035, 3.155, 1, 1, 1, 1, 1, 1, 1, 1}, {0.875, 0.915, 1, 1, 1, 1, 1, 1, 1, 1}, {0.875, 0.915, 1, 1, 1, 1, 1, 1, 1, 1}, {0.650, 0.900, 1, 1, 1, 1, 1, 1, 1, 1}};
  Configurable<LabeledArray<double>> cfgPrimVtxParticleCuts{"cfgPrimVtxParticleCuts", {&DefaulPrimVtxCutValues[0][0], 6, 10, {"Phi1020", "JPsiToEE", "JPsiToMuMu", "KStar892", "KStar892Bar", "Rho770"}, {"massLow", "massUp", "doRecoCheck", "fillPreSelQA", "fillPreSelDauQA", "fillPostSelQA", "fillPostSelDauQA", "checkPrimVtxContm", "checkV0DecayContm", "checkMPMSigma"}}, "cut values for primary vertex reconstruction"};

  static constexpr double DefaulTrackCountConfigValues[4][8] = {{1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1}};
  // Track Configurables
  struct : ConfigurableGroup {
    Configurable<int> cfgTrk01TpcNClsCrossedRows{"cfgTrk01TpcNClsCrossedRows", 70, "cfgTrk01TpcNClsCrossedRows"};
    Configurable<float> cfgTrk02dcaXY{"cfgTrk02dcaXY", 0.2, "cfgTrk02dcaXY"};
    Configurable<float> cfgTrk03dcaZ{"cfgTrk03dcaZ", 2.0, "cfgTrk03dcaZ"};
    Configurable<float> cfgTrk04Eta{"cfgTrk04Eta", 0.8, "cfgTrk04Eta"};
    Configurable<float> cfgTrk05PtLow{"cfgTrk05PtLow", 0.15, "cfgTrk05PtLow"};
    Configurable<float> cfgTrk06PtHigh{"cfgTrk06PtHigh", 2.0, "cfgTrk06PtHigh"};
    Configurable<bool> cfgTrk07DoVGselTrackCheck{"cfgTrk07DoVGselTrackCheck", false, "cfgTrk07DoVGselTrackCheck"};
    Configurable<LabeledArray<double>> countSetting{"countSetting", {&DefaulTrackCountConfigValues[0][0], 4, 8, {"Pi", "Ka", "Pr", "El"}, {"countPrompt", "countNonPrompt", "countDauContam", "countRejected", "fillPrompt", "fillNonPrompt", "fillDauContam", "fillRejected"}}, "configurables for identified particle count"}; //
  } cfgTrackCuts;

  // Configurables for particle Identification
  static constexpr double DefaulPIDcheckValues[6][10] = {{0.7, 0, 2.0, 2.0, 4.0, 0, 2.0, 2.0, 4.0, 1}, {0.8, 0, 2.0, 2.0, 4.0, 0, 2.0, 2.0, 4.0, 1}, {0.8, 0, 2.0, 2.0, 4.0, 0, 2.0, 2.0, 4.0, 1}, {0.4, 0, 2.0, 2.0, 4.0, 0, 2.0, 2.0, 4.0, 1}, {0.4, 0, 2.0, 2.0, 4.0, 0, 2.0, 2.0, 4.0, 1}, {0.4, 0, 2.0, 2.0, 4.0, 0, 2.0, 2.0, 4.0, 1}};
  static constexpr double DefaulPidVetoValues[6][4] = {{1, 1, 3.0, 3.0}, {1, 1, 3.0, 3.0}, {1, 1, 3.0, 3.0}, {0, 0, 3.0, 3.0}, {0, 0, 3.0, 3.0}, {0, 0, 3.0, 3.0}};
  struct : ConfigurableGroup {
    Configurable<bool> cfgId01CheckVetoCut{"cfgId01CheckVetoCut", false, "cfgId01CheckVetoCut"};
    Configurable<bool> cfgId02DoElRejection{"cfgId02DoElRejection", true, "cfgId02DoElRejection"};
    Configurable<bool> cfgId03DoDeRejection{"cfgId03DoDeRejection", false, "cfgId03DoDeRejection"};
    Configurable<bool> cfgId04DoPdependentId{"cfgId04DoPdependentId", true, "cfgId04DoPdependentId"};
    Configurable<bool> cfgId05DoTpcInnerParamId{"cfgId05DoTpcInnerParamId", false, "cfgId05DoTpcInnerParamId"};
    Configurable<LabeledArray<double>> pidConfigSetting{"pidConfigSetting", {&DefaulPIDcheckValues[0][0], 6, 10, {"Pi", "Ka", "Pr", "El", "Mu", "De"}, {"ThrPforTOF", "IdCutTypeLowP", "NSigmaTPCLowP", "NSigmaTOFLowP", "NSigmaRadLowP", "IdCutTypeHighP", "NSigmaTPCHighP", "NSigmaTOFHighP", "NSigmaRadHighP", "doVetoOthers"}}, "cut values for particle identification"}; //
    Configurable<LabeledArray<double>> pidVetoSetting{"pidVetoSetting", {&DefaulPidVetoValues[0][0], 6, 4, {"Pi", "Ka", "Pr", "El", "Mu", "De"}, {"doVetoTPC", "doVetoTOF", "vetoTPC", "vetoTOF"}}, "veto cut for particle Identifiation"};                                                                                                                                    //
  } cfgIdCut;

  struct : ConfigurableGroup {

  } cfgVetoIdCut;

  // configurable for process functions to  reduce memory usage
  struct : ConfigurableGroup {
    Configurable<bool> cfgFill01V0TableFull{"cfgFill01V0TableFull", true, "cfgFill01V0TableFull"};
    Configurable<bool> cfgFill02V0TablePostK0sCheck{"cfgFill02V0TablePostK0sCheck", true, "cfgFill02V0TablePostK0sCheck"};
    Configurable<bool> cfgFill03v0TablePostK0sMassCut{"cfgFill03v0TablePostK0sMassCut", true, "cfgFill03v0TablePostK0sMassCut"};
    Configurable<bool> cfgFill04V0TablePostK0sSelectionCut{"cfgFill04V0TablePostK0sSelectionCut", true, "cfgFill04V0TablePostK0sSelectionCut"};

    Configurable<bool> cfgFill05V0TablePostLambdaCheck{"cfgFill05V0TablePostLambdaCheck", true, "cfgFill05V0TablePostLambdaCheck"};
    Configurable<bool> cfgFill06V0TablePostAntiLambdaCheck{"cfgFill06V0TablePostAntiLambdaCheck", true, "cfgFill06V0TablePostAntiLambdaCheck"};
    Configurable<bool> cfgFill07V0TablePostGammaCheck{"cfgFill07V0TablePostGammaCheck", true, "cfgFill07V0TablePostGammaCheck"};

    Configurable<bool> cfgFill05RecoK0sPreSel{"cfgFill05RecoK0sPreSel", true, "cfgFill05RecoK0sPreSel"};
    Configurable<bool> cfgFill06RecoK0sPostSel{"cfgFill06RecoK0sPostSel", true, "cfgFill06RecoK0sPostSel"};

    Configurable<bool> cfgFill07RecoTrackPreSel{"cfgFill07RecoTrackPreSel", true, "cfgFill07RecoTrackPreSel"};
    Configurable<bool> cfgFill08RecoTrackPostSel{"cfgFill08RecoTrackPostSel", true, "cfgFill08RecoTrackPostSel"};

    Configurable<bool> cfgFill09PiQA{"cfgFill09PiQA", true, "cfgFill09PiQA"};
    Configurable<bool> cfgFill10KaQA{"cfgFill10KaQA", true, "cfgFill10KaQA"};
    Configurable<bool> cfgFill11PrQA{"cfgFill11PrQA", true, "cfgFill11PrQA"};
    Configurable<bool> cfgFill12ElQA{"cfgFill12ElQA", true, "cfgFill12ElQA"};
    Configurable<bool> cfgFill13DeQA{"cfgFill13DeQA", true, "cfgFill13DeQA"};
  } cfgFill;

  struct : ConfigurableGroup {
    Configurable<bool> cfgSim01VtxZCheck{"cfgSim01VtxZCheck", 0, "cfgSim01VtxZCheck"};
    Configurable<bool> cfgSim02CountFinalParticles{"cfgSim02CountFinalParticles", 1, "cfgSim02CountFinalParticles"};
    Configurable<bool> cfgSim03CountNonFinalParticles{"cfgSim03CountNonFinalParticles", 0, "cfgSim03CountNonFinalParticles"};
    Configurable<bool> cfgSim04CountPhysicalPrimAndFinalParticles{"cfgSim04CountPhysicalPrimAndFinalParticles", 0, "cfgSim04CountPhysicalPrimAndFinalParticles"};
    Configurable<bool> cfgSim05doFWDPtDependentCheck{"cfgSim05doFWDPtDependentCheck", 1, "cfgSim05doFWDPtDependentCheck"};
    Configurable<float> cfgSim06FWDPtCut{"cfgSim06FWDPtCut", 0.2, "cfgSim06FWDPtCut"};
    Configurable<std::vector<int>> cfgSim07FinalParticleIdList{"cfgSim07FinalParticleIdList", {11, -11, 13, -13, 15, -15, 211, -211, 321, -321, 2212, -2212}, "cfgSim07FinalParticleIdList"};
    Configurable<std::vector<int>> cfgSim08NonFinalParticleIdList{"cfgSim08NonFinalParticleIdList", {11, -11, 13, -13, 15, -15, 211, -211, 321, -321, 2212, -2212}, "cfgSim08NonFinalParticleIdList"};
    Configurable<float> cfgSim09V0CLow{"cfgSim09V0CLow", -3.7, "cfgSim09V0CLow"};
    Configurable<float> cfgSim09V0CUp{"cfgSim09V0CUp", -1.7, "cfgSim09V0CUp"};
    Configurable<float> cfgSim09V0ALow{"cfgSim09V0ALow", 2.8, "cfgSim09V0ALow"};
    Configurable<float> cfgSim09V0AUp{"cfgSim09V0AUp", 5.1, "cfgSim09V0AUp"};
  } cfgSim;

  struct : ConfigurableGroup {
    Configurable<std::string> cfgCCDB01URL{"cfgCCDB01URL", "http://ccdb-test.cern.ch:8080", "cfgCCDB01URL"};
    Configurable<std::string> cfgCCDB02Path{"cfgCCDB02Path", "Users/r/raverma/KaonIsospinFluctuation/correctionWeights/", "cfgCCDB02Path"};
    Configurable<int> cfgCCDB03SOR{"cfgCCDB03SOR", 1, "cfgCCDB03SOR"};
    Configurable<std::string> cfgCCDB04StrHistP{"cfgCCDB04StrHistP", "h01_p", "cfgCCDB04StrHistP"};
    Configurable<std::string> cfgCCDB05StrHistPt{"cfgCCDB05StrHistPt", "h02_pt", "cfgCCDB05StrHistPt"};
    Configurable<std::string> cfgCCDB06StrHistPtEta{"cfgCCDB06StrHistPtEta", "h20_pt_eta", "cfgCCDB06StrHistPtEta"};
  } cfgCCDB;

  struct : ConfigurableGroup {
    Configurable<int> cfgEffCorrType = {"cfgEffCorrType", 2, "cfgEffCorrType"};
    // Add new Configurables related to efficency axis in future
  } cfgEffCorr;

  static constexpr int DefaultSparseHistValues[3][2] = {{4, 1}, {6, 1}, {6, 1}};
  static constexpr double DefaultAxisValues[31][3] = {
    {500, -1.5, 498.5 /*nPi*/},
    {500, -1.5, 498.5 /*nKa*/},
    {500, -1.5, 498.5 /*nPr*/},
    {500, -1.5, 498.5 /*nEl*/},
    {500, -1.5, 498.5 /*nMu*/},
    {500, -1.5, 498.5 /*nDe*/},
    {100, -1.5, 98.5 /*nK0s*/},
    {100, -1.5, 98.5 /*nLambda*/},
    {100, -1.5, 98.5 /*nAntiLambda*/},
    {100, -1.5, 98.5 /*nGamma*/},
    {100, -1.5, 98.5 /*nPhi1020*/},
    {100, -1.5, 98.5 /*nJPsiToEE*/},
    {100, -1.5, 98.5 /*nJPsiToMuMu*/},
    {100, -1.5, 98.5 /*nKStar892*/},
    {100, -1.5, 98.5 /*nKStar892Bar*/},
    {100, -1.5, 98.5 /*nRho770*/},
    {100, -1.5, 98.5 /*nRejectedPi*/},
    {100, -1.5, 98.5 /*nRejectedKa*/},
    {100, -1.5, 98.5 /*nRejectedPr*/},
    {100, -1.5, 98.5 /*nRejectedEl*/},
    {100, -1.5, 98.5 /*nRejectedMu*/},
    {100, -1.5, 98.5 /*nRejectedDe*/},
    {500, -1.5, 498.5 /*nKaon*/},
    {10000, -1.5, 9998.5 /*(nK0s)^{2}*/},
    {250000, -1.5, 249998.5 /*(nKaon)^{2}*/},
    {500, -1.5, 498.5 /*(nK0s*nKaon)*/},
    {250000, -1.5, 249998.5 /*(nKaPlus)^{2}*/},
    {250000, -1.5, 249998.5 /*(nKaMinus)^{2}*/},
    {250000, -1.5, 249998.5 /*(nKaPlus*nKaMinus)*/},
    {1020, -1.0, 101.0 /*centrality*/},
    {2000, -1.5, 1998.5 /*nTrack*/}};

  struct : ConfigurableGroup {
    Configurable<LabeledArray<int>> sparseSetting{"sparseSetting", {&DefaultSparseHistValues[0][0], 3, 2, {"hSparse0", "hSparse1", "hSparse2"}, {"nAxis", "sparseType"}}, "configuration of sparse histogram"};
    Configurable<LabeledArray<std::string>> sparseAxisStr{"sparseAxisStr", {&defaultSparseAxis[0][0], 10, 3, {"Axis0", "Axis1", "Axis2", "Axis3", "Axis4", "Axis5", "Axis6", "Axis7", "Axis8", "Axis9"}, {"sparse1", "sparse2", "sparse3"}}, "configuration of sparse histogramee"};
    // Configurable<std::vector<std::string>> cfgSparse1qqAxis{"cfgSparse1qqAxis", {"axis1", "col1", "axis2", "col2"}, "cfgSparse1qqAxis"};

    Configurable<LabeledArray<double>> axisSetting{"axisSetting", {&DefaultAxisValues[0][0], 31, 3, {"nPi", "nKa", "nPr", "nEl", "nMu", "nDe", "nK0s", "nLambda", "nAntiLambda", "nGamma", "nPhi1020", "nJPsiToEE", "nJPsiToMuMu", "nKStar892", "nKStar892Bar", "nRho770", "nRejectedPi", "nRejectedKa", "nRejectedPr", "nRejectedEl", "nRejectedMu", "nRejectedDe", "nKaon", "(nK0s)^{2}", "(nKaon)^{2}", "(nK0s*nKaon)", "(nKaPlus)^{2}", "(nKaMinus)^{2}", "(nKaPlus*nKaMinus)", "centrality", "nTrack"}, {"nBins", "xLow", "xUp"}}, "particle count axis bin configuration"};
  } cfgAxis;

  std::string getModifiedStr(const std::string& myString)
  {
    size_t pos = myString.rfind('/');
    if (pos != std::string::npos) {
      std::string subString = myString.substr(0, pos); // remove "/" from end of the string
      return subString;
    } else {
      return myString;
    }
  }

  enum V0ParicleTypeEnum {
    kV0TrkK0s = 0,
    kV0TrkLambda,
    kV0TrkAntiLambda,
    kV0TrkGamma,
    v0TrkEnumSize
  };

  enum V0ParicleCutTypeEnum {
    kV0MLow = 0,
    kV0MUp,
    kV0LowPt,
    kV0HighPt,
    kV0Rapitidy,
    kV0ARMcut1,
    kV0ARMcut2,
    kV0DoRecoCheck,
    kV0FillPreSel,
    kV0FillPreSelDau,
    kV0FillPostSel,
    kV0FillPostSelDau,
    kV0CheckPrimVtxContm,
    kV0CheckV0DecayContm,
    kV0CheckMPMSigma
  };

  enum FillFlagEnum {
    kFillPreSel = 0,
    kFillPostSel,
    kFillSimple
  };

  static constexpr std::string_view FillModeDire[] = {
    "PreSel/",
    "PostSel/",
    ""};

  enum IdentifiedParticleCountFlagEnum {
    kCountPrompt = 0,
    kCountNonPrompt,
    kCountDauContam,
    kCountRejected,
    kFillPrompt,
    kFillNonPrompt,
    kFillDauContam,
    kFillRejected
  };

  enum ProcessTypeEnum {
    doDataProcessing = 0,
    doRecoProcessing,
    doPurityProcessing,
    doGenProcessing,
    doSimProcessing,
    processTypeEnumSize
  };

  static constexpr std::string_view ProcessTypeDire[] = {
    "DataProcessing",
    "RecoProcessing",
    "PurityProcessing",
    "GenProcessing",
    "SimProcessing"};

  enum PidEnum {
    kPi = 0, // dont use kPion, kKaon, as these enumeration
    kKa,     // are already defined in $ROOTSYS/root/include/TPDGCode.h
    kPr,
    kEl,
    kMu,
    kDe,
    kV0K0s,
    kV0Lambda,
    kV0AntiLambda,
    kV0Gamma,
    kPrimPhi1020,
    kPrimJPsiToEE,
    kPrimJPsiToMuMu,
    kPrimKStar892,
    kPrimKStar892Bar,
    kPrimRho770,
    kRejectedPi,
    kRejectedKa,
    kRejectedPr,
    kRejectedEl,
    kRejectedMu,
    kRejectedDe,
    kNKaon,
    kNK0sSq,
    kNKaonSq,
    kNK0sNKaonProd,
    kNKaPosSq,
    kNKaNegSq,
    kNKaPosNKaNegProd,
    kCentrality,
    kNTrack,
    kNonePid,
    pidEnumSize
  };

  static constexpr std::string_view PidDire[] = {
    "Pi/",
    "Ka/",
    "Pr/",
    "El/",
    "Mu/",
    "De/",
    "V0K0s/",
    "V0Lambda/",
    "V0AntiLambda/",
    "V0Gamma/",
    "PrimPhi1020/",
    "PrimJPsiToEE/",
    "PrimJPsiToMuMu/",
    "PrimKStar892/",
    "PrimKStar892Bar/",
    "PrimRho770/",
    "RejectedPi/",
    "RejectedKa/",
    "RejectedPr/",
    "RejectedEl/",
    "RejectedMu/",
    "RejectedDe/",
    "nKaon",
    "(nK0s)^{2}",
    "(nKaon)^{2}",
    "(nK0s*nKaon)",
    "(nKaPlus)^{2}",
    "(nKaMinus)^{2}",
    "(nKaPlus*nKaMinus)",
    "centrality",
    "nTrack",
    ""};

  enum CountAxisEnum {
    kNBin = 0,
    kXLow,
    kXUp
  };

  enum CutSettingEnum {
    kThrPforTOF = 0,
    kIdCutTypeLowP,
    kNSigmaTPCLowP,
    kNSigmaTOFLowP,
    kNSigmaRadLowP,
    kIdCutTypeHighP,
    kNSigmaTPCHighP,
    kNSigmaTOFHighP,
    kNSigmaRadHighP,
    kDoVetoOthers
  };

  enum VetoSettingEnum {
    kDoVetoTPC = 0,
    kDoVetoTOF,
    kVetoTPC,
    kVetoTOF
  };

  enum SignModeEnum {
    kPos = 0,
    kNeg,
    kNeu,
    kNoneSign,
    signModeEnumSize
  };

  static constexpr std::string_view SignDire[] = {
    "Pos/",
    "Neg/",
    "Neu/",
    ""};

  enum CountTypeEnum {
    kIntCount = 0,
    kFloatCount,
    kEffWeightSum,
    countTypeEnumSize
  };

  static constexpr std::string_view CountDire[] = {
    "intCount/",
    "floatCount/",
    "effWeightSum/"};

  enum RejectionTagEnum {
    kPassed = 0,
    kFailTpcNClsCrossedRows,
    kFailTrkdcaXY,
    kFailTrkdcaZ,
    kFailGlobalTrack,
    kFailVGSelCheck,
    kFailV0K0ShortDaughter,
    kFailV0LambdaDaughter,
    kFailV0AntiLambdaDaughter,
    kFailPhi1020Daughter,
    kFailJPsiToEEDaughter,
    kFailJPsiToMuMuDaughter,
    kFailKStar892Daughter,
    kFailKStar892BarDaughter,
    kFailkPrimTrkRho770Daughter,
    kRejectionTagEnumSize
  };

  static constexpr std::string_view TrackTagDire[]{
    "kPassed",
    "kFailTpcNClsCrossedRows",
    "kFailTrkdcaXY",
    "kFailTrkdcaZ",
    "kFailGlobalTrack",
    "kFailVGSelCheck",
    "kFailV0K0ShortDaughter",
    "kFailV0LambdaDaughter",
    "kFailV0AntiLambdaDaughter",
    "kFailPhi1020Daughter",
    "kFailJPsiToEEDaughter",
    "kFailJPsiToMuMuDaughter",
    "kFailKStar892Daughter",
    "kFailKStar892BarDaughter",
    "kFailkPrimTrkRho770Daughter",
  };

  enum CentEnum {
    kCentFT0C = 0,
    kCentFT0M,
    kMultFT0M,
    kMultFT0C
  };

  // Function to split a string by '/' and return tokens
  std::vector<std::string> tokenizeBySlash(const std::string& input)
  {
    std::vector<std::string> tokens;
    std::stringstream ss(input);
    std::string token;

    while (std::getline(ss, token, '/')) {
      if (!token.empty()) {
        tokens.push_back(token);
      }
    }

    return tokens;
  }

  void getSparseInfoFromTokens(const auto& cfgList, const auto& nameTokens, auto& cfgAxisNumber, auto& signType, auto& variableType, const int& row, const int& col)
  {
    // nameTokens[0] is the particle type;
    if (nameTokens.size() == 0) {
      cfgAxisNumber = -1;
      signType = -1;
      variableType = -1;
      return;
    }
    for (uint i = 0; i < cfgList->getLabelsRows().size(); i++) {
      if (nameTokens[0] == cfgList->getLabelsRows()[i]) {
        cfgAxisNumber = i;
        break;
      }
      cfgAxisNumber = -1;
    }
    static const int three = 3;
    if (nameTokens.size() <= 1) {
      signType = kPos;
      variableType = -1;
    } else if (nameTokens.size() > three) {
      LOG(fatal) << "DEBUG :: invalid token name :: row = " << row << " :: col = " << col;
    } else {
      // nameToken[1] is the signType
      if (nameTokens[1] == "kPos") {
        signType = kPos;
      } else if (nameTokens[1] == "kNeg") {
        signType = kNeg;
      } else {
        signType = -1;
        LOG(fatal) << "DEBUG :: unidentifed sparse name value :: signType = " << nameTokens[1] << " :: ERROR in row = " << row << " :: col = " << col;
      }

      // nameToken[2] is the variableType
      if (nameTokens[2] == "iNTrk") {
        variableType = kIntCount;
      } else if (nameTokens[2] == "fNTrk") {
        variableType = kFloatCount;
      } else if (nameTokens[2] == "effSum") {
        variableType = kEffWeightSum;
      } else {
        variableType = -1;
        LOG(fatal) << "DEBUG :: unidentifed sparse name value :: variable Type = " << nameTokens[2] << " :: ERROR in row = " << row << " :: col = " << col;
      }
    }
  }

  // Efficieny containing Histograms.
  TH2F* hPtEtaForBinSearch = nullptr;
  std::vector<std::array<TH2F*, kNeg + 1>> hPtEtaForEffCorrection{kPrimRho770 + 1, std::array<TH2F*, kNeg + 1>{}};

  // 3 sparse histos ==> each with 10 axes ==> each axis with 3 parameter (0 - pidType, 1- SignType, 2-VariableType)
  std::array<std::array<std::array<int, 3>, 10>, 3> sparseHistFillerIndices{};

  void init(InitContext const&)
  {
    TList* fOutputList = new TList();
    fOutputList->SetOwner(true);
    fOutput.setObject(fOutputList);
    listHandler(fOutputList, cfgDebug.makeNestedTList);

    primVtxTLHs[0](fOutputList, cfgDebug.makeNestedTList);
    primVtxTLHs[1](fOutputList, cfgDebug.makeNestedTList);
    primVtxTLHs[2](fOutputList, cfgDebug.makeNestedTList);
    primVtxTLHs[3](fOutputList, cfgDebug.makeNestedTList);
    primVtxTLHs[4](fOutputList, cfgDebug.makeNestedTList);
    primVtxTLHs[5](fOutputList, cfgDebug.makeNestedTList);

    auto& mgr = o2::ccdb::BasicCCDBManager::instance();
    mgr.setURL(cfgCCDB.cfgCCDB01URL);
    mgr.setCaching(true);
    auto ccdbObj = mgr.getForTimeStamp<TList>(cfgCCDB.cfgCCDB02Path, cfgCCDB.cfgCCDB03SOR);
    if (!ccdbObj) {
      if (cfgDebug.printDebugMessages)
        LOG(info) << "DEBUG :: CCDB OBJECT NOT FOUND";
    } else {
      if (cfgDebug.printDebugMessages)
        LOG(info) << "DEBUG :: CCDB OBJECT FOUND";
    }

    ccdbObj->Print();

    hPtEtaForBinSearch = reinterpret_cast<TH2F*>(ccdbObj->FindObject("hPtEta"));
    if (cfgDebug.printDebugMessages)
      LOG(info) << "DEBUG :: Obj Name = " << hPtEtaForBinSearch->GetName() << " :: entries = " << hPtEtaForBinSearch->GetEntries();
    std::string name = "";
    for (int i = 0; i <= kPrimRho770; i++) {
      for (int j = 0; j < (kNeg + 1); j++) {
        if (i > kDe) {
          if (j == kNeg) {
            continue;
          }
        }
        name = "hPtEta" + getModifiedStr(static_cast<std::string>(PidDire[i])) + getModifiedStr(static_cast<std::string>(SignDire[j]));
        hPtEtaForEffCorrection[i][j] = reinterpret_cast<TH2F*>(ccdbObj->FindObject(name.c_str()));
        if (cfgDebug.printDebugMessages)
          LOG(info) << "DEBUG :: Obj Name = " << hPtEtaForEffCorrection[i][j]->GetName() << " :: entries = " << hPtEtaForBinSearch->GetEntries();
      }
    }

    mgr.setURL("http://alice-ccdb.cern.ch"); // RESET the URL otherwise the other process functions which contains ccdb lookups will fail

    // Axes
    const AxisSpec axisV0Mass = {1800, -0.1f, 1.7f, "#it{M}_{V0} [GeV/#it{c}^{2}]"};
    const AxisSpec axisK0sMass = {200, 0.40f, 0.60f, "#it{M}_{inv}^{K_{0}^{s}} [GeV/#it{c}^{2}]"};
    const AxisSpec axisLambdaMass = {200, 1.f, 1.2f, "#it{M}_{inv}^{#Lambda} [GeV/#it{c}^{2}]"};
    const AxisSpec axisGammaMass = {120, -0.01f, 0.05f, "#it{M}_{inv}^{#gamma} [GeV/#it{c}^{2}]"};

    const AxisSpec axisPhiMass = {200, 0.99f, 1.07f, "#it{M}_{inv}^{#phi} [GeV/#it{c}^{2}]"};
    const AxisSpec axisJPsiMass = {200, 3.00f, 3.20f, "#it{M}_{inv}^{J/#psi} [GeV/#it{c}^{2}]"};
    const AxisSpec axisKStar892Mass = {200, 0.860f, 0.910f, "#it{M}_{inv}^{K^{*}(892)} [GeV/#it{c}^{2}]"};
    const AxisSpec axisRho770Mass = {200, 0.600, 0.950f, "#it{M}_{inv}^{#rho(770)} [GeV/#it{c}^{2}]"};

    const AxisSpec variableAxis = {{VARIABLE_WIDTH, -10, -5, -3, -1, 1, 3, 5, 10}, "variableAxis"};
    const AxisSpec timingAxisUpto50 = {510, -1, 50, "seconds"};
    const AxisSpec timingAxisUpto10 = {110, -1, 10, "seconds"};
    const AxisSpec timingAxisUpto5 = {300, -1, 5, "seconds"};
    const AxisSpec timingAxisUpto1 = {110, -0.1, 1, "seconds"};

    const AxisSpec axisVertexZ = {30, -15., 15., "vrtx_{Z} [cm]"};

    AxisSpec axisCent = {cfgCentAxis.centAxis01Bins, cfgCentAxis.centAxis02XLow, cfgCentAxis.centAxis03XUp, "centFT0C(percentile)"};

    const AxisSpec axisP = {200, 0.0f, 10.0f, "#it{p} (GeV/#it{c})"};
    const AxisSpec axisPt = {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisTPCInnerParam = {200, 0.0f, 10.0f, "#it{p}_{tpcInnerParam} (GeV/#it{c})"};
    const AxisSpec axisTOFExpMom = {200, 0.0f, 10.0f, "#it{p}_{tofExpMom} (GeV/#it{c})"};

    const AxisSpec axisEta = {100, -5, 5, "#eta"};
    const AxisSpec axisPhi = {110, -1, 10, "#phi (radians)"};
    const AxisSpec axisRapidity = {200, -5, 5, "Rapidity (y)"};
    const AxisSpec axisDcaXY = {100, -5, 5, "dcaXY"};
    const AxisSpec axisDcaZ = {100, -5, 5, "dcaZ"};
    const AxisSpec axisDcaXYwide = {2000, -100, 100, "dcaXY"};
    const AxisSpec axisDcaZwide = {2000, -100, 100, "dcaZ"};
    const AxisSpec axisSign = {10, -5, 5, "track.sign"};

    const AxisSpec axisTPCSignal = {100, -1, 1000, "tpcSignal"};
    const AxisSpec axisTOFBeta = {40, -2.0, 2.0, "tofBeta"};

    const AxisSpec axisTPCSignalFine = {10010, -1, 1000, "tpcSignal"};
    const AxisSpec axisTOFBetaFine = {10010, -1, 1000, "tpcSignal"};

    const AxisSpec axisTPCNSigma = {200, -10.0, 10.0, "n#sigma_{TPC}"};
    const AxisSpec axisTOFNSigma = {200, -10.0, 10.0, "n#sigma_{TOF}"};

    const AxisSpec axisTPCNSigmaPi = {200, -10.0, 10.0, "n#sigma_{TPC}^{Pi}"};
    const AxisSpec axisTOFNSigmaPi = {200, -10.0, 10.0, "n#sigma_{TOF}^{Pi}"};

    const AxisSpec axisTPCNClsCrossedRows = {200, -1.5, 198.5, "tpcNClsCrossedRows"};
    const AxisSpec axisIsPVContributor = {4, -1, 3, "isPVContributor"};
    const AxisSpec axisIsGlobalTrack = {4, -1, 3, "isGobalTrack"};
    const AxisSpec axisIsK0sDau = {4, -1, 3, "isK0sDau"};

    const AxisSpec axisDcapostopv = {1000, -50, 50, "dcapostopv"};
    const AxisSpec axisDcanegtopv = {1000, -50, 50, "dcanegtopv"};
    const AxisSpec axisDcaV0daughters = {2000, -10.0, 10.0, "dcaV0daughters"};
    const AxisSpec axisV0cosPA = {1400, -0.1, 1.3, "v0cosPA"};
    const AxisSpec axisV0radius = {600, -10, 50, "v0radius"};

    const AxisSpec axisIParticleCount1 = {60, -10, 50, "particleCount"};
    const AxisSpec axisIParticleCount2 = {260, -10, 250, "particleCount"};
    const AxisSpec axisIParticleCount3 = {1060, -10, 1050, "particleCount"};

    const AxisSpec axisfParticleCount1 = {600, -10, 50, "particleCount{float}"};
    const AxisSpec axisfParticleCount2 = {2600, -10, 250, "particleCount{float}"};

    const AxisSpec axisEffSum1 = {600, -10, 50, "effWeightSum"};
    const AxisSpec axisEffSum2 = {2600, -10, 250, "effWeightSum"};

    const AxisSpec axisArmenterosAlpha = {100, -1.0, 1.0, "ArmenterosAlpha"};
    const AxisSpec axisArmenterosQt = {150, 0, 0.3, "ArmenterosQt"};

    const AxisSpec axisCentBins = {103, -1.5, 101.5, "centrality"};
    const AxisSpec axisPtBins = {110, -0.5, 10.5, "p_{T} (GeV/c)"};

    HistogramConfigSpec histPDcaXY({HistType::kTH2F, {axisP, axisDcaXY}});
    HistogramConfigSpec histPtDcaXY({HistType::kTH2F, {axisPt, axisDcaXY}});
    HistogramConfigSpec histTpcInnerParamDcaXY({HistType::kTH2F, {axisTPCInnerParam, axisDcaXY}});
    HistogramConfigSpec histTofExpMomDcaXY({HistType::kTH2F, {axisTOFExpMom, axisDcaXY}});

    HistogramConfigSpec histPDcaZ({HistType::kTH2F, {axisP, axisDcaZ}});
    HistogramConfigSpec histPtDcaZ({HistType::kTH2F, {axisPt, axisDcaZ}});
    HistogramConfigSpec histTpcInnerParamDcaZ({HistType::kTH2F, {axisTPCInnerParam, axisDcaZ}});
    HistogramConfigSpec histTofExpMomDcaZ({HistType::kTH2F, {axisTOFExpMom, axisDcaZ}});

    HistogramConfigSpec histPPt({HistType::kTH2F, {axisP, axisPt}});
    HistogramConfigSpec histPTpcInnerParam({HistType::kTH2F, {axisP, axisTPCInnerParam}});
    HistogramConfigSpec histPTofExpMom({HistType::kTH2F, {axisP, axisTOFExpMom}});

    HistogramConfigSpec histPTpcSignal({HistType::kTH2F, {axisP, axisTPCSignal}});
    HistogramConfigSpec histTpcInnerParamTpcSignal({HistType::kTH2F, {axisTPCInnerParam, axisTPCSignal}});
    HistogramConfigSpec histTofExpMomTpcSignal({HistType::kTH2F, {axisTOFExpMom, axisTPCSignal}});

    HistogramConfigSpec histPBeta({HistType::kTH2F, {axisP, axisTOFBeta}});
    HistogramConfigSpec histTpcInnerParamBeta({HistType::kTH2F, {axisTPCInnerParam, axisTOFBeta}});
    HistogramConfigSpec histTofExpMomBeta({HistType::kTH2F, {axisTOFExpMom, axisTOFBeta}});

    HistogramConfigSpec histPTpcNSigma({HistType::kTH2F, {axisP, axisTPCNSigma}});
    HistogramConfigSpec histPtTpcNSigma({HistType::kTH2F, {axisPt, axisTPCNSigma}});
    HistogramConfigSpec histTpcInnerParamTpcNSigma({HistType::kTH2F, {axisTPCInnerParam, axisTPCNSigma}});
    HistogramConfigSpec histTofExpMomTpcNSigma({HistType::kTH2F, {axisTOFExpMom, axisTPCNSigma}});
    HistogramConfigSpec histPTofNSigma({HistType::kTH2F, {axisP, axisTOFNSigma}});
    HistogramConfigSpec histPtTofNSigma({HistType::kTH2F, {axisPt, axisTOFNSigma}});
    HistogramConfigSpec histTpcInnerParamTofNSigma({HistType::kTH2F, {axisTPCInnerParam, axisTOFNSigma}});
    HistogramConfigSpec histTofExpMomTofNSigma({HistType::kTH2F, {axisTOFExpMom, axisTOFNSigma}});
    HistogramConfigSpec histTpcNSigmaTofNSigma({HistType::kTH2F, {axisTPCNSigma, axisTOFNSigma}});

    auto addV0Histos = [axisV0Mass, axisK0sMass, axisLambdaMass, axisGammaMass,
                        axisDcapostopv, axisDcanegtopv, axisDcaV0daughters, axisV0cosPA, axisV0radius,
                        axisP, axisPt, axisEta, axisPhi, axisRapidity, axisArmenterosAlpha, axisArmenterosQt](auto& histReg, const std::string& basePath) {
      histReg.add((basePath + "h01_K0s_Mass").c_str(), "K0s_Mass", {HistType::kTH1F, {axisK0sMass}});
      histReg.add((basePath + "h02_Lambda_Mass").c_str(), "Lambda_Mass", {HistType::kTH1F, {axisLambdaMass}});
      histReg.add((basePath + "h03_AntiLambda_Mass").c_str(), "AntiLambda_Mass", {HistType::kTH1F, {axisLambdaMass}});
      histReg.add((basePath + "h04_V0Gamma_Mass").c_str(), "Gamma_Mass", {HistType::kTH1F, {axisGammaMass}});
      histReg.add((basePath + "h04_v0DaughterCollisionIndexTag").c_str(), "hV0s_K0s_v0DaughterCollisionIndexTag", {HistType::kTH1D, {{26, -1.0, 12.0}}});
      histReg.add((basePath + "h04_v0DauBCTag").c_str(), "hV0s_v0DauBCTag", {HistType::kTH1D, {{26, -1.0, 12.0}}});
      histReg.add((basePath + "h05_V0Tag").c_str(), "V0Tag", {HistType::kTH1F, {{37, -2, 35}}}); // 001 = Kaon, 010 = Lambda, 100 = AnitLambda

      // Topological Cuts
      histReg.add((basePath + "h06_dcapostopv").c_str(), "dcapostopv", kTH1F, {axisDcapostopv});
      histReg.add((basePath + "h07_dcanegtopv").c_str(), "dcanegtopv", kTH1F, {axisDcanegtopv});
      histReg.add((basePath + "h08_dcaV0daughters").c_str(), "dcaV0daughters", kTH1F, {axisDcaV0daughters});
      histReg.add((basePath + "h09_v0cosPA").c_str(), "v0cosPA", kTH1F, {axisV0cosPA});
      histReg.add((basePath + "h10_v0radius").c_str(), "v0radius", kTH1F, {axisV0radius});

      // K0s-FullInformation

      if (basePath == "v0Table/postK0sCheck/" || basePath == "v0Table/postK0sMassCut/" || basePath == "v0Table/postK0sSelectionCut/") {
        histReg.add((basePath + "h11_mass").c_str(), "mass", kTH1F, {axisK0sMass});
      } else if (basePath == "v0Table/postLambdaCheck/" || basePath == "v0Table/postAntiLambdaCheck/") {
        histReg.add((basePath + "h11_mass").c_str(), "mass", kTH1F, {axisLambdaMass});
      } else if (basePath == "v0Table/postGammaCheck/") {
        histReg.add((basePath + "h11_mass").c_str(), "mass", kTH1F, {axisGammaMass});
      } else if (basePath == "v0Table/Full/") {
        histReg.add((basePath + "h11_mass").c_str(), "mass", kTH1F, {axisV0Mass});
      }

      histReg.add((basePath + "h12_p").c_str(), "p", kTH1F, {axisP});
      histReg.add((basePath + "h13_pt").c_str(), "pt", kTH1F, {axisPt});
      histReg.add((basePath + "h14_eta").c_str(), "eta", kTH1F, {axisEta});
      histReg.add((basePath + "h15_phi").c_str(), "phi", kTH1F, {axisPhi});
      histReg.add((basePath + "h16_rapidity0_hypMass").c_str(), "rapidity0_hypMass", kTH1F, {axisRapidity});
      histReg.add((basePath + "h16_rapidity1_invMass").c_str(), "rapidity1_invMass", kTH1F, {axisRapidity});
      histReg.add((basePath + "h17_alpha").c_str(), "alpha", kTH1F, {axisArmenterosAlpha});
      histReg.add((basePath + "h18_qtarm").c_str(), "qtarm", kTH1F, {axisArmenterosQt});
      histReg.add((basePath + "h19_alpha_qtarm").c_str(), "alpha_qtarm", kTH2F, {axisArmenterosAlpha, axisArmenterosQt});
      histReg.add((basePath + "h20_pt_eta").c_str(), "pt_eta", kTH2F, {axisPt, axisEta});
      histReg.add((basePath + "h20_pt_eta_effCorr").c_str(), "pt_eta_effCorr", kTH2F, {axisPt, axisEta});
    };

    auto addTrackQAHistos = [axisP, axisPt, axisTPCInnerParam, axisTOFExpMom, axisEta, axisPhi, axisRapidity,
                             axisIsPVContributor, axisIsGlobalTrack, axisSign, axisDcaXY, axisDcaZ, axisDcaXYwide, axisDcaZwide](auto& histReg, const std::string& basePath) {
      histReg.add((basePath + "h01_p").c_str(), "p", kTH1F, {axisP});
      histReg.add((basePath + "h02_pt").c_str(), "pt", kTH1F, {axisPt});
      histReg.add((basePath + "h03_tpcInnerParam").c_str(), "tpcInnerParam", kTH1F, {axisTPCInnerParam});
      histReg.add((basePath + "h04_tofExpMom").c_str(), "tofExpMom", kTH1F, {axisTOFExpMom});
      histReg.add((basePath + "h05_eta").c_str(), "eta", kTH1F, {axisEta});
      histReg.add((basePath + "h06_phi").c_str(), "phi", kTH1F, {axisPhi});
      histReg.add((basePath + "h07_rapidity").c_str(), "rapidity", kTH1F, {axisRapidity});
      histReg.add((basePath + "h08_isPVContributor").c_str(), "isPVContributor", kTH1F, {axisIsPVContributor});
      histReg.add((basePath + "h09_isGlobalTrack").c_str(), "isGlobalTrack", kTH1F, {axisIsGlobalTrack});
      histReg.add((basePath + "h10_sign").c_str(), "sign", kTH1D, {axisSign});
      histReg.add((basePath + "h11_dcaXY").c_str(), "dcaXY", kTH1F, {axisDcaXY});
      histReg.add((basePath + "h12_dcaZ").c_str(), "dcaZ", kTH1F, {axisDcaZ});

      histReg.add((basePath + "h17_dcaXYwide").c_str(), "dcaXYwide", kTH1F, {axisDcaXYwide});
      histReg.add((basePath + "h18_dcaZwide").c_str(), "dcaZwide", kTH1F, {axisDcaZwide});
      histReg.add((basePath + "h20_pt_eta").c_str(), "pt_eta", kTH2F, {axisPt, axisEta});
      histReg.add((basePath + "h20_pt_eta_effCorr").c_str(), "pt_eta_effCorr", kTH2F, {axisPt, axisEta});
    };

    auto addTrackMomentumQAHistos = [histPPt, histPTpcInnerParam, histPTofExpMom,
                                     histPTpcSignal, histTpcInnerParamTpcSignal, histTofExpMomTpcSignal,
                                     histPBeta, histTpcInnerParamBeta, histTofExpMomBeta](auto& histReg, const std::string& basePath) {
      histReg.add((basePath + "h21_p_pt").c_str(), "p_pt", histPPt);
      histReg.add((basePath + "h22_p_tpcInnerParam").c_str(), "p_tpcInnerParam", histPTpcInnerParam);
      histReg.add((basePath + "h23_p_tofExpMom").c_str(), "p_tofExpMom", histPTofExpMom);
      // tpcSignal
      histReg.add((basePath + "h24_p_tpcSignal").c_str(), "p_tpcSignal", histPTpcSignal);
      histReg.add((basePath + "h25_tpcInnerParam_tpcSignal").c_str(), "tpcInnerParam_tpcSignal", histTpcInnerParamTpcSignal);
      histReg.add((basePath + "h26_tofExpMom_tpcSignal").c_str(), "tofExpMom_tpcSignal", histTofExpMomTpcSignal);
      // tofBeta
      histReg.add((basePath + "h27_p_beta").c_str(), "p_beta", histPBeta);
      histReg.add((basePath + "h28_tpcInnerParam_beta").c_str(), "tpcInnerParam_beta", histTpcInnerParamBeta);
      histReg.add((basePath + "h29_tofExpMom_beta").c_str(), "tofExpMom_beta", histTofExpMomBeta);
    };

    auto addTrackIdQAHistos = [histPTpcNSigma, histPtTpcNSigma, histTpcInnerParamTpcNSigma, histTofExpMomTpcNSigma,
                               histPTofNSigma, histPtTofNSigma, histTpcInnerParamTofNSigma, histTofExpMomTofNSigma, histTpcNSigmaTofNSigma](auto& histReg, const std::string& basePath) {
      // Look at Identification
      histReg.add((basePath + "h31_p_tpcNSigma").c_str(), "p_tpcNSigma", histPTpcNSigma);
      histReg.add((basePath + "h32_pt_tpcNSigma").c_str(), "pt_tpcNSigma", histPtTpcNSigma);
      histReg.add((basePath + "h33_tpcInnerParam_tpcNSigma").c_str(), "tpcInnerParam_tpcNSigma", histTpcInnerParamTpcNSigma);
      histReg.add((basePath + "h34_tofExpMom_tpcNSigma").c_str(), "tofExpMom_tpcNSigma", histTofExpMomTpcNSigma);
      histReg.add((basePath + "h35_p_tofNSigma").c_str(), "p_tofNSigma", histPTofNSigma);
      histReg.add((basePath + "h36_pt_tofNSigma").c_str(), "pt_tofNSigma", histPtTofNSigma);
      histReg.add((basePath + "h37_tpcInnerParam_tofNSigma").c_str(), "tpcInnerParam_tofNSigma", histTpcInnerParamTofNSigma);
      histReg.add((basePath + "h38_tofExpMom_tofNSigma").c_str(), "tofExpMom_tofNSigma", histTofExpMomTofNSigma);
      histReg.add((basePath + "h39_tpcNSigma_tofNSigma").c_str(), "tpcNSigma_tofNSigma", histTpcNSigmaTofNSigma);
    };

    auto addTrackDcaQAHistos = [histPDcaXY, histPtDcaXY, histTpcInnerParamDcaXY, histTofExpMomDcaXY,
                                histPDcaZ, histPtDcaZ, histTpcInnerParamDcaZ, histTofExpMomDcaZ](auto& histReg, const std::string& basePath) {
      // DcaXY
      histReg.add((basePath + "h41_p_dcaXY").c_str(), "p_dcaXY", histPDcaXY);
      histReg.add((basePath + "h42_pt_dcaXY").c_str(), "pt_dcaXY", histPtDcaXY);
      histReg.add((basePath + "h43_tpcInnerParam_dcaXY").c_str(), "tpcInnerParam_dcaXY", histTpcInnerParamDcaXY);
      histReg.add((basePath + "h44_tofExpMom_dcaXY").c_str(), "tofExpMom_dcaXY", histTofExpMomDcaXY);

      // DcaZ
      histReg.add((basePath + "h45_p_dcaZ").c_str(), "p_dcaZ", histPDcaZ);
      histReg.add((basePath + "h46_pt_dcaZ").c_str(), "pt_dcaZ", histPtDcaZ);
      histReg.add((basePath + "h47_tpcInnerParam_dcaZ").c_str(), "tpcInnerParam_dcaZ", histTpcInnerParamDcaZ);
      histReg.add((basePath + "h48_tofExpMom_dcaZ").c_str(), "tofExpMom_dcaZ", histTofExpMomDcaZ);
    };

    auto addAllTrackQAHistos = [addTrackQAHistos, addTrackMomentumQAHistos, addTrackIdQAHistos, addTrackDcaQAHistos](auto& histReg, const std::string& basePath) {
      addTrackQAHistos(histReg, basePath);
      addTrackMomentumQAHistos(histReg, basePath);
      addTrackIdQAHistos(histReg, basePath);
      addTrackDcaQAHistos(histReg, basePath);
    };

    auto addAllV0Histos = [addV0Histos, addAllTrackQAHistos](auto& histReg, const std::string& basePath, const std::string& posPidPath, const std::string& negPidPath, const bool fillMotherQA, const bool fillDauTrackQA) {
      if (fillMotherQA) {
        addV0Histos(histReg, basePath);
      }
      if (fillDauTrackQA) {
        addAllTrackQAHistos(histReg, basePath + posPidPath + "Pos/tpcId/");                                                  // K0s daughter info
        histReg.addClone((basePath + posPidPath + "Pos/tpcId/").c_str(), (basePath + posPidPath + "Pos/tpctofId/").c_str()); // for identification using tof+tpc
        histReg.addClone((basePath + posPidPath + "Pos/").c_str(), (basePath + negPidPath + "Neg/").c_str());                // For both daughters positive and negative
      }
    };

    auto addCountHist = [](auto& histReg, const std::string& basePath, const std::string& histTitle) {
      histReg.add(basePath.c_str(), histTitle.c_str(), {HistType::kTH1D, {{44, -2, 20}}});
    };

    auto addCountLabels = [](const auto& h) {
      h->GetXaxis()->SetBinLabel(h->GetXaxis()->FindBin(int(0)), "Total Repeats / Fake Counts");
      static const int twenty = 20;
      for (int i = 1; i < twenty; i++) {
        h->GetXaxis()->SetBinLabel(h->GetXaxis()->FindBin(i), Form("Counted %d Times", i));
      }
    };

    // Full V0 Table information -- First V0 Loop
    addAllV0Histos(recoV0sFull, "v0Table/Full/", "Pi/", "Pi/", true, true);
    recoV0sFull.addClone("v0Table/Full/Pi/Pos/tpcId/", "v0Table/Full/Pi/Pos/NoId/");
    recoV0sFull.addClone("v0Table/Full/Pi/Pos/NoId/", "v0Table/Full/Pi/Neg/NoId/");
    recoV0sFull.add("v0Table/Full/hTrueV0TagCount", "hTrueV0TagCount", {HistType::kTH1F, {{12, -2, 10}}});

    if (cfgFill.cfgFill02V0TablePostK0sCheck) {
      addAllV0Histos(recoV0sPostK0sCheck, "v0Table/postK0sCheck/", "Pi/", "Pi/", true, true);
    }
    if (cfgFill.cfgFill03v0TablePostK0sMassCut) {
      addAllV0Histos(recoV0sPostMassCut, "v0Table/postK0sMassCut/", "Pi/", "Pi/", true, true);
    }

    if (cfgFill.cfgFill04V0TablePostK0sSelectionCut) {
      addAllV0Histos(recoV0sPostK0sSelectionCut, "v0Table/postK0sSelectionCut/", "Pi/", "Pi/", true, true);
      addCountHist(recoV0sPostK0sSelectionCut, "v0Table/postK0sSelectionCut/nCommonDauOfDifferentK0s", "nCommonDauOfDifferentK0s");
      addCountLabels(recoV0sPostK0sSelectionCut.get<TH1>(HIST("v0Table/postK0sSelectionCut/nCommonDauOfDifferentK0s")));
    }

    if (cfgFill.cfgFill05V0TablePostLambdaCheck) {
      addAllV0Histos(recoV0sPostLambdaCheck, "v0Table/postLambdaCheck/", "Pr/", "Pi/", true, true);
      addCountHist(recoV0sPostLambdaCheck, "v0Table/postLambdaCheck/nCommonDauOfDifferentLambdas", "nCommonDauOfDifferentLambdas");
      addCountLabels(recoV0sPostLambdaCheck.get<TH1>(HIST("v0Table/postLambdaCheck/nCommonDauOfDifferentLambdas")));
    }

    if (cfgFill.cfgFill06V0TablePostAntiLambdaCheck) {
      addAllV0Histos(recoV0sPostAntiLambdaCheck, "v0Table/postAntiLambdaCheck/", "Pi/", "Pr/", true, true);
      addCountHist(recoV0sPostAntiLambdaCheck, "v0Table/postAntiLambdaCheck/nCommonDauOfDifferentAntiLambdas", "nCommonDauOfDifferentAntiLambdas");
      addCountLabels(recoV0sPostAntiLambdaCheck.get<TH1>(HIST("v0Table/postAntiLambdaCheck/nCommonDauOfDifferentAntiLambdas")));
    }

    if (cfgFill.cfgFill07V0TablePostGammaCheck) {
      addAllV0Histos(recoV0sPostGammaCheck, "v0Table/postGammaCheck/", "El/", "El/", true, true);
      addCountHist(recoV0sPostGammaCheck, "v0Table/postGammaCheck/nCommonDauOfDifferentGammas", "nCommonDauOfDifferentGammas");
      addCountLabels(recoV0sPostGammaCheck.get<TH1>(HIST("v0Table/postGammaCheck/nCommonDauOfDifferentGammas")));
    }

    // Event Selection
    recoEvent.add("recoEvent/ProcessType", "ProcessType", {HistType::kTH1D, {{20, -1, 9}}});
    recoEvent.add("recoEvent/PreSel/CollType", "CollType", {HistType::kTH1D, {{20, -1, 9}}});
    recoEvent.add("recoEvent/PostSel/CollType", "CollType", {HistType::kTH1D, {{20, -1, 9}}});

    recoEvent.add("recoEvent/h01_CollisionCount", "CollisionCount", {HistType::kTH1D, {{1, 0, 1}}});
    recoEvent.add("recoEvent/h02_VertexXRec", "VertexXRec", {HistType::kTH1D, {{1000, -0.2, 0.2}}});
    recoEvent.add("recoEvent/h03_VertexYRec", "VertexYRec", {HistType::kTH1D, {{1000, -0.2, 0.2}}});
    recoEvent.add("recoEvent/h04_VertexZRec", "VertexZRec", {HistType::kTH1F, {axisVertexZ}});
    recoEvent.add("recoEvent/h05_Centrality", "Centrality", {HistType::kTH1F, {axisCent}});
    recoEvent.add("recoEvent/h06_V0Size", "V0Size", {HistType::kTH1F, {{60, -10, 50}}});
    recoEvent.add("recoEvent/h07_TracksSize", "TracksSize", {HistType::kTH1F, {axisIParticleCount2}});
    recoEvent.add("recoEvent/h08_nTrack", "nTrack", {HistType::kTH1F, {axisIParticleCount2}});

    auto addParticleCountHistos = [this,
                                   axisIParticleCount1, axisIParticleCount2,
                                   axisfParticleCount1, axisfParticleCount2,
                                   axisEffSum1, axisEffSum2](auto& histReg, const std::string& basePath) {
      std::string histTitle = "";
      std::string signLabel = "";
      std::string prefixLabel = "";
      auto axisOfCounts = axisIParticleCount1;
      for (int k = 0; k < countTypeEnumSize; k++) {
        if (k == kIntCount) {
          prefixLabel = "iN";
        }
        if (k == kFloatCount) {
          prefixLabel = "fN";
        }
        if (k == kEffWeightSum) {
          prefixLabel = "effWeightSum";
        }
        for (int i = 0; i < (kRejectedDe + 1); i++) {
          if (i == kNonePid) {
            continue;
          }
          for (int j = kPos; j <= kNeg; j++) {
            signLabel = SignDire[j].data();
            if (kDe < i && i < kRejectedPi) {
              signLabel = "";
              if (j == 1)
                continue;
            }
            if (i == kPi || i == kRejectedPi) {
              if (k == kIntCount) {
                axisOfCounts = axisIParticleCount2;
              }
              if (k == kFloatCount) {
                axisOfCounts = axisfParticleCount2;
              }
              if (k == kEffWeightSum) {
                axisOfCounts = axisEffSum2;
              }
            } else {
              if (k == kIntCount) {
                axisOfCounts = axisIParticleCount1;
              }
              if (k == kFloatCount) {
                axisOfCounts = axisfParticleCount1;
              }
              if (k == kEffWeightSum) {
                axisOfCounts = axisEffSum1;
              }
            }
            std::ostringstream histName;
            histName << basePath << CountDire[k] << PidDire[i] << signLabel << "hCount";
            histTitle = prefixLabel + getModifiedStr(static_cast<std::string>(PidDire[i])) + getModifiedStr(signLabel);
            histReg.add(histName.str().c_str(), histTitle.c_str(), kTH1F, {axisOfCounts});
          }
        }
      }
    };

    addParticleCountHistos(recoEvent, "recoEvent/");

    for (int i = 0; i < processTypeEnumSize; i++) {
      recoEvent.get<TH1>(HIST("recoEvent/ProcessType"))->GetXaxis()->SetBinLabel(recoEvent.get<TH1>(HIST("recoEvent/ProcessType"))->FindBin(i), ProcessTypeDire[i].data());
    }

    auto addPrimVtxParticleQAHistos = [axisPhiMass, axisJPsiMass, axisKStar892Mass, axisRho770Mass,
                                       axisP, axisPt, axisEta, axisPhi, axisRapidity](auto& histReg, const std::string& basePath) {
      if (basePath == "recoPhi1020/PreSel/" || basePath == "recoPhi1020/PostSel/") {
        histReg.add((basePath + "h11_mass").c_str(), "mass", kTH1F, {axisPhiMass});
      } else if (basePath == "recoJPsiToEE/PreSel/" || basePath == "recoJPsiToEE/PostSel/") {
        histReg.add((basePath + "h11_mass").c_str(), "mass", kTH1F, {axisJPsiMass});
      } else if (basePath == "recoJPsiToMuMu/PreSel/" || basePath == "recoJPsiToMuMu/PostSel/") {
        histReg.add((basePath + "h11_mass").c_str(), "mass", kTH1F, {axisJPsiMass});
      } else if (basePath == "recoKStar892/PreSel/" || basePath == "recoKStar892/PostSel/") {
        histReg.add((basePath + "h11_mass").c_str(), "mass", kTH1F, {axisKStar892Mass});
      } else if (basePath == "recoKStar892Bar/PreSel/" || basePath == "recoKStar892Bar/PostSel/") {
        histReg.add((basePath + "h11_mass").c_str(), "mass", kTH1F, {axisKStar892Mass});
      } else if (basePath == "recoRho770/PreSel/" || basePath == "recoRho770/PostSel/") {
        histReg.add((basePath + "h11_mass").c_str(), "mass", kTH1F, {axisRho770Mass});
      }

      histReg.add((basePath + "h12_p").c_str(), "p", kTH1F, {axisP});
      histReg.add((basePath + "h13_pt").c_str(), "pt", kTH1F, {axisPt});
      histReg.add((basePath + "h14_eta").c_str(), "eta", kTH1F, {axisEta});
      histReg.add((basePath + "h15_phi").c_str(), "phi", kTH1F, {axisPhi});
      histReg.add((basePath + "h16_rapidity1_invMass").c_str(), "rapidity1_invMass", kTH1F, {axisRapidity});
      histReg.add((basePath + "h20_pt_eta").c_str(), "pt_eta", kTH2F, {axisPt, axisEta});
      histReg.add((basePath + "h20_pt_eta_effCorr").c_str(), "pt_eta_effCorr", kTH2F, {axisPt, axisEta});
    };

    auto addAllPrimVtxParticleHistos = [addPrimVtxParticleQAHistos, addAllTrackQAHistos](auto& histReg, const std::string& basePath, const std::string& posPidPath, const std::string& negPidPath, const bool fillMotherQA, const bool fillDauTrackQA) {
      if (fillMotherQA) {
        addPrimVtxParticleQAHistos(histReg, basePath);
      }
      if (fillDauTrackQA) {
        addAllTrackQAHistos(histReg, basePath + posPidPath + "Pos/tpcId/");                                                  // Daughter info
        histReg.addClone((basePath + posPidPath + "Pos/tpcId/").c_str(), (basePath + posPidPath + "Pos/tpctofId/").c_str()); // for identification using tof+tpc
        histReg.addClone((basePath + posPidPath + "Pos/").c_str(), (basePath + negPidPath + "Neg/").c_str());                // For both daughters positive and negative
      }
    };

    // Pre Sel and post Sel Histos related to particle originating from primary vertex
    std::string basePath = "", posPidPath = "", negPidPath = "", commonCountName = "", commonCountTitle = "";
    bool fillMotherQA = false, fillDauTrackQA = false;
    // Pre Sel and post Sel Histos  for Phi1020, J/psi to ee, J/psi to MuMu, K*(892), K*(892)Bar, Rho770
    for (int iPidType = kPrimTrkPhi1020; iPidType <= kPrimTrkRho770; iPidType++) {
      if (iPidType == kPrimTrkPhi1020) {
        basePath = "recoPhi1020/";
        posPidPath = "Ka/";
        negPidPath = "Ka/";
        commonCountName = "recoPhi1020/PreSel/nCommonKaOfDifferentPhi";
        commonCountTitle = "nCommonKaOfDifferentPhi";
      }
      if (iPidType == kPrimTrkJPsiToEE) {
        basePath = "recoJPsiToEE/";
        posPidPath = "El/";
        negPidPath = "El/";
        commonCountName = "recoJPsiToEE/PreSel/nCommonKaOfDifferentJPsiToEE";
        commonCountTitle = "nCommonKaOfDifferentJPsiToEE";
      }
      if (iPidType == kPrimTrkJPsiToMuMu) {
        basePath = "recoJPsiToMuMu/";
        posPidPath = "Mu/";
        negPidPath = "Mu/";
        commonCountName = "recoJPsiToMuMu/PreSel/nCommonKaOfDifferentJPsiToMuMu";
        commonCountTitle = "nCommonKaOfDifferentJPsiToMuMu";
      }
      if (iPidType == kPrimTrkKStar892) {
        basePath = "recoKStar892/";
        posPidPath = "Ka/";
        negPidPath = "Pi/";
        commonCountName = "recoKStar892/PreSel/nCommonKaOfDifferentKStar892";
        commonCountTitle = "nCommonKaOfDifferentKStar892";
      }
      if (iPidType == kPrimTrkKStar892Bar) {
        basePath = "recoKStar892Bar/";
        posPidPath = "Pi/";
        negPidPath = "Ka/";
        commonCountName = "recoKStar892Bar/PreSel/nCommonKaOfDifferentKStar892Bar";
        commonCountTitle = "nCommonKaOfDifferentKStar892Bar";
      }
      if (iPidType == kPrimTrkRho770) {
        basePath = "recoRho770/";
        posPidPath = "Pi/";
        negPidPath = "Pi/";
        commonCountName = "recoRho770/PreSel/nCommonKaOfDifferentRho770";
        commonCountTitle = "nCommonKaOfDifferentRho770";
      }

      if (getCfg<bool>(cfgPrimVtxParticleCuts, iPidType, kPrimDoRecoCheck)) {
        if (getCfg<bool>(cfgPrimVtxParticleCuts, iPidType, kPrimFillPreSel)) {
          fillMotherQA = getCfg<bool>(cfgPrimVtxParticleCuts, iPidType, kPrimFillPreSel);
          fillDauTrackQA = getCfg<bool>(cfgPrimVtxParticleCuts, iPidType, kPrimFillPreSelDau);
          addAllPrimVtxParticleHistos(primVtxTLHs[iPidType], basePath + "PreSel/", posPidPath, negPidPath, fillMotherQA, fillDauTrackQA);
          addCountHist(primVtxTLHs[iPidType], commonCountName, commonCountTitle);
          if (iPidType == kPrimTrkPhi1020) {
            addCountLabels(primVtxTLHs[iPidType].get<TH1>(HIST("recoPhi1020/PreSel/nCommonKaOfDifferentPhi")));
          }
          if (iPidType == kPrimTrkJPsiToEE) {
            addCountLabels(primVtxTLHs[iPidType].get<TH1>(HIST("recoJPsiToEE/PreSel/nCommonKaOfDifferentJPsiToEE")));
          }
          if (iPidType == kPrimTrkJPsiToMuMu) {
            addCountLabels(primVtxTLHs[iPidType].get<TH1>(HIST("recoJPsiToMuMu/PreSel/nCommonKaOfDifferentJPsiToMuMu")));
          }
          if (iPidType == kPrimTrkKStar892) {
            addCountLabels(primVtxTLHs[iPidType].get<TH1>(HIST("recoKStar892/PreSel/nCommonKaOfDifferentKStar892")));
          }
          if (iPidType == kPrimTrkKStar892Bar) {
            addCountLabels(primVtxTLHs[iPidType].get<TH1>(HIST("recoKStar892Bar/PreSel/nCommonKaOfDifferentKStar892Bar")));
          }
          if (iPidType == kPrimTrkRho770) {
            addCountLabels(primVtxTLHs[iPidType].get<TH1>(HIST("recoRho770/PreSel/nCommonKaOfDifferentRho770")));
          }
        }
        if (getCfg<bool>(cfgPrimVtxParticleCuts, iPidType, kPrimFillPostSel)) {
          fillMotherQA = getCfg<bool>(cfgPrimVtxParticleCuts, iPidType, kPrimFillPostSel);
          fillDauTrackQA = getCfg<bool>(cfgPrimVtxParticleCuts, iPidType, kPrimFillPostSelDau);
          addAllPrimVtxParticleHistos(primVtxTLHs[iPidType], basePath + "PostSel/", posPidPath, negPidPath, fillMotherQA, fillDauTrackQA);
        }
      }
    }

    // K0s Table
    if (getCfg<bool>(cfgPrimVtxParticleCuts, kV0TrkK0s, kV0DoRecoCheck)) {
      fillMotherQA = getCfg<bool>(cfgPrimVtxParticleCuts, kV0TrkK0s, kV0FillPreSel);
      fillDauTrackQA = getCfg<bool>(cfgPrimVtxParticleCuts, kV0TrkK0s, kV0FillPreSelDau);
      addAllV0Histos(recoK0s, "recoK0s/PreSel/", "Pi/", "Pi/", fillMotherQA, fillDauTrackQA);
      if (fillDauTrackQA) {
        recoK0s.addClone("recoK0s/PreSel/Pi/Pos/tpcId/", "recoK0s/PreSel/Pi/Pos/NoId/");
      } // for unidentified
      if (fillDauTrackQA) {
        recoK0s.addClone("recoK0s/PreSel/Pi/Pos/NoId/", "recoK0s/PreSel/Pi/Neg/NoId/");
      } // for unidentified

      fillMotherQA = getCfg<bool>(cfgPrimVtxParticleCuts, kV0TrkK0s, kV0FillPostSel);
      fillDauTrackQA = getCfg<bool>(cfgPrimVtxParticleCuts, kV0TrkK0s, kV0FillPostSelDau);
      addAllV0Histos(recoK0s, "recoK0s/PostSel/", "Pi/", "Pi/", fillMotherQA, fillDauTrackQA);
      recoK0s.add("recoK0s/PostSel/mK0s_vs_cent", "mK0s_vs_cent", kTH2F, {axisCent, axisK0sMass});
      recoK0s.add("recoK0s/PostSel/mK0s_vs_cent_effcorr", "mK0s_vs_cent_effcorr", kTH2F, {axisCent, axisK0sMass});
      recoK0s.add("recoK0s/PostSel/mK0s_vs_cent_vs_Pt", "mK0s_vs_cent_vs_Pt", kTH3F, {axisK0sMass, axisCentBins, axisPtBins});
    }

    if (getCfg<bool>(cfgPrimVtxParticleCuts, kV0TrkLambda, kV0DoRecoCheck)) {
      fillMotherQA = getCfg<bool>(cfgPrimVtxParticleCuts, kV0TrkK0s, kV0FillPostSel);
      fillDauTrackQA = getCfg<bool>(cfgPrimVtxParticleCuts, kV0TrkK0s, kV0FillPostSelDau);
      addAllV0Histos(listHandler, "recoLambda/PostSel/", "Pr/", "Pi/", fillMotherQA, fillDauTrackQA);
      listHandler.add("recoLambda/PostSel/mLambda_vs_cent", "mLambda_vs_cent", kTH2F, {axisCent, axisLambdaMass});
      listHandler.add("recoLambda/PostSel/mLambda_vs_cent_effcorr", "mLambda_vs_cent_effcorr", kTH2F, {axisCent, axisLambdaMass});
      listHandler.add("recoLambda/PostSel/mLambda_vs_cent_vs_Pt", "mLambda_vs_cent_vs_Pt", kTH3F, {axisLambdaMass, axisCentBins, axisPtBins});
    }

    if (getCfg<bool>(cfgPrimVtxParticleCuts, kV0TrkAntiLambda, kV0DoRecoCheck)) {
      fillMotherQA = getCfg<bool>(cfgPrimVtxParticleCuts, kV0TrkAntiLambda, kV0FillPostSel);
      fillDauTrackQA = getCfg<bool>(cfgPrimVtxParticleCuts, kV0TrkAntiLambda, kV0FillPostSelDau);
      addAllV0Histos(listHandler, "recoAntiLambda/PostSel/", "Pi/", "Pr/", fillMotherQA, fillDauTrackQA);
      listHandler.add("recoAntiLambda/PostSel/mAntiLambda_vs_cent", "mAntiLambda_vs_cent", kTH2F, {axisCent, axisLambdaMass});
      listHandler.add("recoAntiLambda/PostSel/mAntiLambda_vs_cent_effcorr", "mAntiLambda_vs_cent_effcorr", kTH2F, {axisCent, axisLambdaMass});
      listHandler.add("recoAntiLambda/PostSel/mAntiLambda_vs_cent_vs_Pt", "mAntiLambda_vs_cent_vs_Pt", kTH3F, {axisLambdaMass, axisCentBins, axisPtBins});
    }

    if (getCfg<bool>(cfgPrimVtxParticleCuts, kV0TrkGamma, kV0DoRecoCheck)) {
      fillMotherQA = getCfg<bool>(cfgPrimVtxParticleCuts, kV0TrkGamma, kV0FillPostSel);
      fillDauTrackQA = getCfg<bool>(cfgPrimVtxParticleCuts, kV0TrkGamma, kV0FillPostSelDau);
      addAllV0Histos(listHandler, "recoGamma/PostSel/", "El/", "El/", fillMotherQA, fillDauTrackQA);
      listHandler.add("recoGamma/PostSel/mGamma_vs_cent", "mGamma_vs_cent", kTH2F, {axisCent, axisGammaMass});
      listHandler.add("recoGamma/PostSel/mGamma_vs_cent_effcorr", "mGamma_vs_cent_effcorr", kTH2F, {axisCent, axisGammaMass});
      listHandler.add("recoGamma/PostSel/mGamma_vs_cent_vs_Pt", "mGamma_vs_cent_vs_Pt", kTH3F, {axisGammaMass, axisCentBins, axisPtBins});
    }

    // Tracks reconstruction
    // FullTrack
    addTrackQAHistos(recoTracks, "recoTracks/PreSel/");
    addTrackMomentumQAHistos(recoTracks, "recoTracks/PreSel/");
    addTrackDcaQAHistos(recoTracks, "recoTracks/PreSel/");
    addTrackIdQAHistos(recoTracks, "recoTracks/PreSel/Pi/Pos/NoId/"); // Look at Pion

    recoTracks.addClone("recoTracks/PreSel/Pi/Pos/", "recoTracks/PreSel/Pi/Neg/"); // Look at Pion
    recoTracks.addClone("recoTracks/PreSel/Pi/", "recoTracks/PreSel/Ka/");         // Kaon
    recoTracks.addClone("recoTracks/PreSel/Pi/", "recoTracks/PreSel/Pr/");         // Proton
    recoTracks.addClone("recoTracks/PreSel/Pi/", "recoTracks/PreSel/El/");         // Electron
    recoTracks.addClone("recoTracks/PreSel/Pi/", "recoTracks/PreSel/De/");         // Deuteron
    // Write Code for naming the axis for Identified Particles in case of full track information

    // Selected Track
    recoTracks.addClone("recoTracks/PreSel/", "recoTracks/PostSel/");
    //

    auto addAllAnalysisTrackHistos = [addAllTrackQAHistos](auto& histReg, const std::string& basePath) {
      addAllTrackQAHistos(histReg, basePath + "Pos/tpcId/");
      histReg.addClone(basePath + "Pos/tpcId/", basePath + "Pos/tpctofId/");
      histReg.addClone(basePath + "Pos/", basePath + "Neg/");
    };
    // Analysis

    // Pi QA
    addAllAnalysisTrackHistos(recoAnalysisPi, "recoAnalysis/Prompt/Pi/");
    recoAnalysisPi.addClone("recoAnalysis/Prompt/Pi/", "recoAnalysis/NonPrompt/Pi/");
    recoAnalysisPi.addClone("recoAnalysis/Prompt/Pi/", "recoAnalysis/DauContam/Pi/");
    addAllAnalysisTrackHistos(recoRejectedPi, "recoAnalysis/Rejected/Pi/");

    // //Ka QA
    addAllAnalysisTrackHistos(recoAnalysisKa, "recoAnalysis/Prompt/Ka/");
    recoAnalysisKa.addClone("recoAnalysis/Prompt/Ka/", "recoAnalysis/NonPrompt/Ka/");
    recoAnalysisKa.addClone("recoAnalysis/Prompt/Ka/", "recoAnalysis/DauContam/Ka/");
    addAllAnalysisTrackHistos(recoRejectedKa, "recoAnalysis/Rejected/Ka/");

    // //Pr QA
    addAllAnalysisTrackHistos(recoAnalysisPr, "recoAnalysis/Prompt/Pr/");
    recoAnalysisPr.addClone("recoAnalysis/Prompt/Pr/", "recoAnalysis/NonPrompt/Pr/");
    recoAnalysisPr.addClone("recoAnalysis/Prompt/Pr/", "recoAnalysis/DauContam/Pr/");
    addAllAnalysisTrackHistos(recoRejectedPr, "recoAnalysis/Rejected/Pr/");

    // //El QA
    addAllAnalysisTrackHistos(recoAnalysisEl, "recoAnalysis/Prompt/El/");
    recoAnalysisEl.addClone("recoAnalysis/Prompt/El/", "recoAnalysis/NonPrompt/El/");
    recoAnalysisEl.addClone("recoAnalysis/Prompt/El/", "recoAnalysis/DauContam/El/");
    addAllAnalysisTrackHistos(recoRejectedEl, "recoAnalysis/Rejected/El/");

    // //De QA
    addAllAnalysisTrackHistos(recoAnalysisDe, "recoAnalysis/Prompt/De/");
    recoAnalysisDe.addClone("recoAnalysis/Prompt/De/", "recoAnalysis/NonPrompt/De/");
    recoAnalysisDe.addClone("recoAnalysis/Prompt/De/", "recoAnalysis/DauContam/De/");
    addAllAnalysisTrackHistos(recoRejectedDe, "recoAnalysis/Rejected/De/");

    // Analysis
    recoAnalysis.add("recoAnalysis/SelectedTrack_IdentificationTag", "SelectedTrack_IdentificationTag", kTH1D, {{34, -1.5, 32.5, "trackTAG"}});
    recoAnalysis.add("recoAnalysis/RejectedTrack_RejectionTag", "RejectedTrack_RejectionTag", kTH1D, {{50, -1.5, 23.5, "rejectionTAG"}});
    // by default the track is accepted
    for (int i = 0; i < kRejectionTagEnumSize; i++) {
      recoAnalysis.get<TH1>(HIST("recoAnalysis/RejectedTrack_RejectionTag"))->GetXaxis()->SetBinLabel(recoAnalysis.get<TH1>(HIST("recoAnalysis/RejectedTrack_RejectionTag"))->FindBin(i), TrackTagDire[i].data());
    }

    recoAnalysis.add("recoAnalysis/centFT0A_vs_NTrk", "centFT0A_vs_NTrk", kTH2F, {{1030, -1.5, 101.5, "centFT0A"}, {203, -1.5, 201.5, "nTrk"}});
    recoAnalysis.add("recoAnalysis/centFT0C_vs_NTrk", "centFT0C_vs_NTrk", kTH2F, {{1030, -1.5, 101.5, "centFT0C"}, {203, -1.5, 201.5, "nTrk"}});
    recoAnalysis.add("recoAnalysis/centFT0M_vs_NTrk", "centFT0M_vs_NTrk", kTH2F, {{1030, -1.5, 101.5, "centFT0M"}, {203, -1.5, 201.5, "nTrk"}});

    recoAnalysis.add("recoAnalysis/centFT0A_vs_centFT0M", "centFT0A_vs_centFT0M", kTH2F, {{1030, -1.5, 101.5, "centFT0A"}, {1030, -1.5, 101.5, "centFT0M"}});
    recoAnalysis.add("recoAnalysis/centFT0C_vs_centFT0A", "centFT0C_vs_centFT0A", kTH2F, {{1030, -1.5, 101.5, "centFT0C"}, {1030, -1.5, 101.5, "centFT0A"}});
    recoAnalysis.add("recoAnalysis/centFT0M_vs_centFT0C", "centFT0M_vs_centFT0C", kTH2F, {{1030, -1.5, 101.5, "centFT0M"}, {1030, -1.5, 101.5, "centFT0C"}});
    // Add histograms for centrality vs multiplicity later

    // Declaring Sparse Histograms
    std::vector<AxisSpec> axisList0(10, {1, 0, 1, "dummyAxis"});
    std::vector<AxisSpec> axisList1(10, {1, 0, 1, "dummyAxis"});
    std::vector<AxisSpec> axisList2(10, {1, 0, 1, "dummyAxis"});

    for (auto& matrix : sparseHistFillerIndices) { // o2-linter: disable=const-ref-in-for-loop (element of 3D array not constant)
      for (auto& row : matrix) {                   // o2-linter: disable=const-ref-in-for-loop (element of 3D array not constant)
        for (auto& val : row) {                    // o2-linter: disable=const-ref-in-for-loop (element of 3D array not constant)
          val = -1;
        }
      }
    }

    int iPidMode;
    int signType;
    int variableType;
    int nBins;
    double xLow;
    double xUp;
    std::string axisTitle, signLabel, countTypeLabel;
    std::vector<std::string> nameTokens;
    for (uint iAxis = 0; iAxis < cfgAxis.sparseAxisStr->getLabelsRows().size(); iAxis++) {
      for (uint iSparse = 0; iSparse < cfgAxis.sparseAxisStr->getLabelsCols().size(); iSparse++) {
        nameTokens = tokenizeBySlash(cfgAxis.sparseAxisStr->get(iAxis, iSparse));
        getSparseInfoFromTokens(cfgAxis.axisSetting, nameTokens, iPidMode, signType, variableType, iAxis, iSparse);

        if (iAxis >= static_cast<uint>(cfgAxis.sparseSetting->get(iSparse, "nAxis")))
          continue;

        nBins = getCfg<int>(cfgAxis.axisSetting, iPidMode, kNBin);
        xLow = getCfg<double>(cfgAxis.axisSetting, iPidMode, kXLow);
        xUp = getCfg<double>(cfgAxis.axisSetting, iPidMode, kXUp);

        signLabel = "";
        countTypeLabel = "";

        if (signType == kPos) {
          signLabel = "Pos";
        } else if (signType == kNeg) {
          signLabel = "Neg";
        }

        if (iPidMode > kRejectedDe) {
          signLabel = "";
        }

        if (variableType == kIntCount) {
          countTypeLabel = "i_";
        } else if (variableType == kFloatCount) {
          countTypeLabel = "f_";
        } else if (variableType == kEffWeightSum) {
          countTypeLabel = "ews_";
        }

        axisTitle = countTypeLabel + cfgAxis.axisSetting->getLabelsRows()[iPidMode] + signLabel;
        static const int two = 2;
        if (iSparse == 0) {
          axisList0[iAxis] = AxisSpec(nBins, xLow, xUp, axisTitle);
        }
        if (iSparse == 1) {
          axisList1[iAxis] = AxisSpec(nBins, xLow, xUp, axisTitle);
        }
        if (iSparse == two) {
          axisList2[iAxis] = AxisSpec(nBins, xLow, xUp, axisTitle);
        }

        sparseHistFillerIndices[iSparse][iAxis][0] = iPidMode;
        sparseHistFillerIndices[iSparse][iAxis][1] = signType;
        sparseHistFillerIndices[iSparse][iAxis][2] = variableType;

        if (iPidMode == -1) {
          LOG(fatal) << "DEBUG :: axis Not found in the List :: iAxis = " << iAxis << " :: iSparse = " << iSparse;
        }
        if (cfgDebug.printDebugMessages)
          LOG(info) << "DEBUG :: row = " << iAxis << " :: col = " << iSparse << " :: iPidMode = " << iPidMode << " :: sign = " << signType << " :: variable = " << variableType;
        if (cfgDebug.printDebugMessages)
          LOG(info) << "DEBUG :: axisTitle == " << axisTitle;
      }
      if (cfgDebug.printDebugMessages)
        LOG(info) << "DEBUG :: ";
    }

    auto createSparseHisotgrams = [](auto& histReg, const std::string& basePath, const std::string& title, const auto& axisList, const int& nAxis, const int& sparseType) {
      std::vector<o2::framework::AxisSpec> selectedAxes(axisList.begin(), axisList.begin() + nAxis);
      if (sparseType == 0) {
        histReg.add((basePath + title).c_str(), title.c_str(), kTHnSparseF, selectedAxes);
      }
      if (sparseType == 1) {
        histReg.add((basePath + title).c_str(), title.c_str(), kTHnSparseD, selectedAxes);
      }
    };

    createSparseHisotgrams(recoAnalysis, "recoAnalysis/", "Sparse0", axisList0, cfgAxis.sparseSetting->get("hSparse0", "nAxis"), cfgAxis.sparseSetting->get("hSparse0", "sparseType"));
    createSparseHisotgrams(recoAnalysis, "recoAnalysis/", "Sparse1", axisList1, cfgAxis.sparseSetting->get("hSparse1", "nAxis"), cfgAxis.sparseSetting->get("hSparse1", "sparseType"));
    createSparseHisotgrams(recoAnalysis, "recoAnalysis/", "Sparse2", axisList2, cfgAxis.sparseSetting->get("hSparse2", "nAxis"), cfgAxis.sparseSetting->get("hSparse2", "sparseType"));

    // Analysis Timing inf
    anlysisTimingInfo.add("time00_dfProcessingTime", "time00_dfProcessingTime", kTH1F, {timingAxisUpto50});
    anlysisTimingInfo.add("time01_executeV0loop", "time01_executeV0loop", kTH1F, {timingAxisUpto10});
    anlysisTimingInfo.add("time02_v0SortTime", "time02_v0SortTime", kTH1F, {timingAxisUpto10});
    anlysisTimingInfo.add("time03_collisionLoopTime", "time03_collisionLoopTime", kTH1F, {timingAxisUpto50});
    anlysisTimingInfo.add("time04_manualGrouping", "time04_manualGrouping", kTH1F, {timingAxisUpto1});
    anlysisTimingInfo.add("time05_primVtxLoop1Time", "time05_primVtxLoop1Time", kTH1F, {timingAxisUpto10});
    anlysisTimingInfo.add("time06_primVtxSortTime", "time06_primVtxSortTime", kTH1F, {timingAxisUpto1});
    anlysisTimingInfo.add("time07_primVtxLoop2Time", "time07_primVtxLoop2Time", kTH1F, {timingAxisUpto10});
    anlysisTimingInfo.add("time08_seconV0LoopTime", "time08_seconV0LoopTime", kTH1F, {timingAxisUpto1});
    anlysisTimingInfo.add("time09_trackLoopTime", "time09_trackLoopTime", kTH1F, {timingAxisUpto10});
    anlysisTimingInfo.add("time10_hSparseFillTime", "time10_hSparseFillTime", kTH1F, {timingAxisUpto10});
    anlysisTimingInfo.add("time11_EventInfoFillTime", "time11_EventInfoFillTime", kTH1F, {timingAxisUpto10});

    genAnalysis.add("genAnalysis/V0K0s/h12_p", "p", kTH1F, {axisP});
    genAnalysis.add("genAnalysis/V0K0s/h13_pt", "pt", kTH1F, {axisPt});
    genAnalysis.add("genAnalysis/V0K0s/h14_eta", "eta", kTH1F, {axisEta});
    genAnalysis.add("genAnalysis/V0K0s/h15_phi", "phi", kTH1F, {axisPhi});
    genAnalysis.add("genAnalysis/V0K0s/h16_rapidity", "rapidity", kTH1F, {axisRapidity});
    genAnalysis.add("genAnalysis/V0K0s/h20_pt_eta", "pt_eta", kTH2F, {axisPt, axisEta});
    genAnalysis.addClone("genAnalysis/V0K0s/", "genAnalysis/Pi/");
    genAnalysis.addClone("genAnalysis/V0K0s/", "genAnalysis/Ka/");
    genAnalysis.addClone("genAnalysis/V0K0s/", "genAnalysis/Pr/");
    genAnalysis.addClone("genAnalysis/V0K0s/", "genAnalysis/El/");
    genAnalysis.addClone("genAnalysis/V0K0s/", "genAnalysis/De/");

    // Printing the Stored Registry information
    if (cfgDebug.printDebugMessages) {
      LOG(info) << "Printing Stored Registry Information";
      anlysisTimingInfo.print();
      recoEvent.print();
      recoV0sFull.print();
      recoV0sPostK0sCheck.print();
      recoV0sPostMassCut.print();
      recoV0sPostK0sSelectionCut.print();
      recoV0sPostLambdaCheck.print();
      recoV0sPostAntiLambdaCheck.print();
      recoV0sPostGammaCheck.print();
      recoPhi1020.print();
      recoJPsiToEE.print();
      recoJPsiToMuMu.print();
      recoKStar892.print();
      recoKStar892Bar.print();
      recoRho770.print();
      recoK0s.print();
      recoTracks.print();
      recoAnalysisPi.print();
      recoAnalysisKa.print();
      recoAnalysisPr.print();
      recoAnalysisEl.print();
      recoAnalysisDe.print();
      recoRejectedPi.print();
      recoRejectedKa.print();
      recoRejectedPr.print();
      recoRejectedEl.print();
      recoRejectedDe.print();
      recoAnalysis.print();
      genAnalysis.print();
      listHandler.print();
      for (uint i = 0; i < primVtxTLHs.size(); i++) {
        primVtxTLHs[i].print();
      }
    }
  }

  enum IdentificationType {
    kTPCidentified = 0,
    kTOFidentified,
    kTPCTOFidentified,
    kUnidentified
  };

  enum TpcTofCutType {
    kRectangularCut = 0,
    kCircularCut,
    kEllipsoidalCut
  };

  enum HistRegEnum {
    v0TableFull = 0,
    v0TablePostK0sCheck,
    v0TablePostK0sMassCut,
    v0TablePostK0sSelectionCut,
    v0TablePostLambdaCheck,
    v0TablePostAntiLambdaCheck,
    v0TablePostGammaCheck,
    Phi1020,
    JPsiToEE,
    JPsiToMuMu,
    KStar892,
    KStar892Bar,
    Rho770,
    recoK0sPreSel,
    recoK0sPostSel,
    recoLambdaPostSel,
    recoAntiLambdaPostSel,
    recoGammaPostSel,
    recoTrackPreSel,
    recoTrackPostSel,
    recoAnalysisPrompt,
    recoAnalysisNonPrompt,
    recoAnalysisDauContam,
    recoAnalysisRejected,
    recoAnalysisDir,
    genAnalysisDir,
    eventPreSel,
    eventPostSel
  };

  static constexpr std::string_view HistRegDire[] = {
    "v0Table/Full/",
    "v0Table/postK0sCheck/",
    "v0Table/postK0sMassCut/",
    "v0Table/postK0sSelectionCut/",
    "v0Table/postLambdaCheck/",
    "v0Table/postAntiLambdaCheck/",
    "v0Table/postGammaCheck/",
    "recoPhi1020/",
    "recoJPsiToEE/",
    "recoJPsiToMuMu/",
    "recoKStar892/",
    "recoKStar892Bar/",
    "recoRho770/",
    "recoK0s/PreSel/",
    "recoK0s/PostSel/",
    "recoLambda/PostSel/",
    "recoAntiLambda/PostSel/",
    "recoGamma/PostSel/",
    "recoTracks/PreSel/",
    "recoTracks/PostSel/",
    "recoAnalysis/Prompt/",
    "recoAnalysis/NonPrompt/",
    "recoAnalysis/DauContam/",
    "recoAnalysis/Rejected/",
    "recoAnalysis/",
    "genAnalysis/",
    "recoEvent/PreSel/",
    "recoEvent/PostSel/"};

  enum DetEnum {
    tpcId = 0,
    tofId,
    tpctofId,
    NoId,
    sumId,
    KNoneDet,
    detEnumSize
  };

  static constexpr std::string_view DetDire[] = {
    "tpcId/",
    "tofId/",
    "tpctofId/",
    "NoId/",
    "sumId/",
    ""};

  // vetoRejection for particles //From Victor Luis Gonzalez Sebastian's analysis note for balance functions
  template <typename T>
  bool selTrackForId(const T& track)
  {
    static const float minus3 = -3.0;
    static const float five = 5.0;
    static const float three = 3.0;
    if (minus3 < track.tpcNSigmaEl() && track.tpcNSigmaEl() < five &&
        std::fabs(track.tpcNSigmaPi()) > three &&
        std::fabs(track.tpcNSigmaKa()) > three &&
        std::fabs(track.tpcNSigmaPr()) > three) {
      return false;
    } else {
      return true;
    }
  }

  // If vetoIdOthers = true; it passed all veto checks
  // If vetoIdOthers = false; it failed veto check with some other particle
  template <int pidMode, typename T>
  bool vetoIdOthersTPC(const T& track)
  {
    // Static is only run once, ever.
    static const bool doVetoTPC[6] = {
      getCfg<bool>(cfgIdCut.pidVetoSetting, kPi, kDoVetoTPC),
      getCfg<bool>(cfgIdCut.pidVetoSetting, kKa, kDoVetoTPC),
      getCfg<bool>(cfgIdCut.pidVetoSetting, kPr, kDoVetoTPC),
      getCfg<bool>(cfgIdCut.pidVetoSetting, kEl, kDoVetoTPC),
      getCfg<bool>(cfgIdCut.pidVetoSetting, kMu, kDoVetoTPC),
      getCfg<bool>(cfgIdCut.pidVetoSetting, kDe, kDoVetoTPC)};

    static const float vetoTPC[6] = {
      getCfg<float>(cfgIdCut.pidVetoSetting, kPi, kVetoTPC),
      getCfg<float>(cfgIdCut.pidVetoSetting, kKa, kVetoTPC),
      getCfg<float>(cfgIdCut.pidVetoSetting, kPr, kVetoTPC),
      getCfg<float>(cfgIdCut.pidVetoSetting, kEl, kVetoTPC),
      getCfg<float>(cfgIdCut.pidVetoSetting, kMu, kVetoTPC),
      getCfg<float>(cfgIdCut.pidVetoSetting, kDe, kVetoTPC)};

    if constexpr (pidMode != kPi)
      if (doVetoTPC[kPi] && std::fabs(track.tpcNSigmaPi()) < vetoTPC[kPi])
        return false;
    if constexpr (pidMode != kKa)
      if (doVetoTPC[kKa] && std::fabs(track.tpcNSigmaKa()) < vetoTPC[kKa])
        return false;
    if constexpr (pidMode != kPr)
      if (doVetoTPC[kPr] && std::fabs(track.tpcNSigmaPr()) < vetoTPC[kPr])
        return false;
    if constexpr (pidMode != kEl)
      if (doVetoTPC[kEl] && std::fabs(track.tpcNSigmaEl()) < vetoTPC[kEl])
        return false;
    if constexpr (pidMode != kMu)
      if (doVetoTPC[kMu] && std::fabs(track.tpcNSigmaMu()) < vetoTPC[kMu])
        return false;
    if constexpr (pidMode != kDe)
      if (doVetoTPC[kDe] && std::fabs(track.tpcNSigmaDe()) < vetoTPC[kDe])
        return false;

    return true;
  }

  template <int pidMode, typename T>
  bool vetoIdOthersTOF(const T& track)
  {
    // Only computed once
    static const bool doVetoTOF[6] = {
      getCfg<bool>(cfgIdCut.pidVetoSetting, kPi, kDoVetoTOF),
      getCfg<bool>(cfgIdCut.pidVetoSetting, kKa, kDoVetoTOF),
      getCfg<bool>(cfgIdCut.pidVetoSetting, kPr, kDoVetoTOF),
      getCfg<bool>(cfgIdCut.pidVetoSetting, kEl, kDoVetoTOF),
      getCfg<bool>(cfgIdCut.pidVetoSetting, kMu, kDoVetoTOF),
      getCfg<bool>(cfgIdCut.pidVetoSetting, kDe, kDoVetoTOF)};

    static const float vetoTOF[6] = {
      getCfg<float>(cfgIdCut.pidVetoSetting, kPi, kVetoTOF),
      getCfg<float>(cfgIdCut.pidVetoSetting, kKa, kVetoTOF),
      getCfg<float>(cfgIdCut.pidVetoSetting, kPr, kVetoTOF),
      getCfg<float>(cfgIdCut.pidVetoSetting, kEl, kVetoTOF),
      getCfg<float>(cfgIdCut.pidVetoSetting, kMu, kVetoTOF),
      getCfg<float>(cfgIdCut.pidVetoSetting, kDe, kVetoTOF)};

    if constexpr (pidMode != kPi)
      if (doVetoTOF[kPi] && std::fabs(track.tofNSigmaPi()) < vetoTOF[kPi])
        return false;
    if constexpr (pidMode != kKa)
      if (doVetoTOF[kKa] && std::fabs(track.tofNSigmaKa()) < vetoTOF[kKa])
        return false;
    if constexpr (pidMode != kPr)
      if (doVetoTOF[kPr] && std::fabs(track.tofNSigmaPr()) < vetoTOF[kPr])
        return false;
    if constexpr (pidMode != kEl)
      if (doVetoTOF[kEl] && std::fabs(track.tofNSigmaEl()) < vetoTOF[kEl])
        return false;
    if constexpr (pidMode != kMu)
      if (doVetoTOF[kMu] && std::fabs(track.tofNSigmaMu()) < vetoTOF[kMu])
        return false;
    if constexpr (pidMode != kDe)
      if (doVetoTOF[kDe] && std::fabs(track.tofNSigmaDe()) < vetoTOF[kDe])
        return false;

    return true;
  }

  template <int pidMode, typename T>
  bool vetoIdOthersTPCTOF(const T& track)
  {
    // If either veto fails, reject
    return vetoIdOthersTPC<pidMode>(track) && vetoIdOthersTOF<pidMode>(track);
  }

  template <int pidMode, typename T>
  bool selIdRectangularCut(const T& track, const float& nSigmaTPC, const float& nSigmaTOF)
  {
    if constexpr (pidMode == kPi) {
      return std::fabs(track.tpcNSigmaPi()) < nSigmaTPC && std::fabs(track.tofNSigmaPi()) < nSigmaTOF;
    } else if constexpr (pidMode == kKa) {
      return std::fabs(track.tpcNSigmaKa()) < nSigmaTPC && std::fabs(track.tofNSigmaKa()) < nSigmaTOF;
    } else if constexpr (pidMode == kPr) {
      return std::fabs(track.tpcNSigmaPr()) < nSigmaTPC && std::fabs(track.tofNSigmaPr()) < nSigmaTOF;
    } else if constexpr (pidMode == kEl) {
      return std::fabs(track.tpcNSigmaEl()) < nSigmaTPC && std::fabs(track.tofNSigmaEl()) < nSigmaTOF;
    } else if constexpr (pidMode == kMu) {
      return std::fabs(track.tpcNSigmaMu()) < nSigmaTPC && std::fabs(track.tofNSigmaMu()) < nSigmaTOF;
    } else if constexpr (pidMode == kDe) {
      return std::fabs(track.tpcNSigmaDe()) < nSigmaTPC && std::fabs(track.tofNSigmaDe()) < nSigmaTOF;
    } else {
      return false;
    }
  }

  template <int pidMode, typename T>
  bool selIdEllipsoidalCut(const T& track, const float& nSigmaTPC, const float& nSigmaTOF)
  {
    if constexpr (pidMode == kPi) {
      float tpc = track.tpcNSigmaPi() / nSigmaTPC;
      float tof = track.tofNSigmaPi() / nSigmaTOF;
      return (tpc * tpc + tof * tof) < 1.0f;
    } else if constexpr (pidMode == kKa) {
      float tpc = track.tpcNSigmaKa() / nSigmaTPC;
      float tof = track.tofNSigmaKa() / nSigmaTOF;
      return (tpc * tpc + tof * tof) < 1.0f;
    } else if constexpr (pidMode == kPr) {
      float tpc = track.tpcNSigmaPr() / nSigmaTPC;
      float tof = track.tofNSigmaPr() / nSigmaTOF;
      return (tpc * tpc + tof * tof) < 1.0f;
    } else if constexpr (pidMode == kEl) {
      float tpc = track.tpcNSigmaEl() / nSigmaTPC;
      float tof = track.tofNSigmaEl() / nSigmaTOF;
      return (tpc * tpc + tof * tof) < 1.0f;
    } else if constexpr (pidMode == kMu) {
      float tpc = track.tpcNSigmaMu() / nSigmaTPC;
      float tof = track.tofNSigmaMu() / nSigmaTOF;
      return (tpc * tpc + tof * tof) < 1.0f;
    } else if constexpr (pidMode == kDe) {
      float tpc = track.tpcNSigmaDe() / nSigmaTPC;
      float tof = track.tofNSigmaDe() / nSigmaTOF;
      return (tpc * tpc + tof * tof) < 1.0f;
    } else {
      return false;
    }
  }

  template <int pidMode, typename T>
  constexpr bool selIdCircularCut(const T& track, const float& nSigmaSquaredRad)
  {
    float tpc = 0.f;
    float tof = 0.f;

    if constexpr (pidMode == kPi) {
      tpc = track.tpcNSigmaPi();
      tof = track.tofNSigmaPi();
    } else if constexpr (pidMode == kKa) {
      tpc = track.tpcNSigmaKa();
      tof = track.tofNSigmaKa();
    } else if constexpr (pidMode == kPr) {
      tpc = track.tpcNSigmaPr();
      tof = track.tofNSigmaPr();
    } else if constexpr (pidMode == kEl) {
      tpc = track.tpcNSigmaEl();
      tof = track.tofNSigmaEl();
    } else if constexpr (pidMode == kMu) {
      tpc = track.tpcNSigmaMu();
      tof = track.tofNSigmaMu();
    } else if constexpr (pidMode == kDe) {
      tpc = track.tpcNSigmaDe();
      tof = track.tofNSigmaDe();
    } else {
      return false; // unknown pidMode
    }

    return (tpc * tpc + tof * tof) < nSigmaSquaredRad;
  }

  // Check if TOF is reliable
  template <typename T>
  inline bool checkReliableTOF(const T& track)
  {
    // which check makes the information of TOF relaiable? should track.beta() be checked? e.g.:
    // return track.hasTOF() && track.beta() > 0.0f;
    return track.hasTOF();
  }

  template <int pidMode, typename T>
  bool idTPC(const T& track, const float& nSigmaTPC)
  {
    static const bool doVetoOthers = getCfg<bool>(cfgIdCut.pidConfigSetting, pidMode, kDoVetoOthers);
    if (doVetoOthers) {
      if (!vetoIdOthersTPC<pidMode>(track)) {
        // If vetoIdOthers = true; it passed all veto checks
        // If vetoIdOthers = false; it failed veto check with some other particle
        return false;
      }
    }
    if constexpr (pidMode == kPi) {
      return std::fabs(track.tpcNSigmaPi()) < nSigmaTPC;
    } else if constexpr (pidMode == kKa) {
      return std::fabs(track.tpcNSigmaKa()) < nSigmaTPC;
    } else if constexpr (pidMode == kPr) {
      return std::fabs(track.tpcNSigmaPr()) < nSigmaTPC;
    } else if constexpr (pidMode == kEl) {
      return std::fabs(track.tpcNSigmaEl()) < nSigmaTPC;
    } else if constexpr (pidMode == kMu) {
      return std::fabs(track.tpcNSigmaMu()) < nSigmaTPC;
    } else if constexpr (pidMode == kDe) {
      return std::fabs(track.tpcNSigmaDe()) < nSigmaTPC;
    } else {
      return false;
    }
  }

  template <int pidMode, typename T>
  bool idTPCTOF(const T& track, const int& pidCutType, const float& nSigmaTPC, const float& nSigmaTOF, const float& nSigmaSquaredRad)
  {
    static const bool doVetoOthers = getCfg<bool>(cfgIdCut.pidConfigSetting, pidMode, kDoVetoOthers);
    if (doVetoOthers) {
      if (!vetoIdOthersTPCTOF<pidMode>(track)) {
        // If vetoIdOthers = true; it passed all veto checks
        // If vetoIdOthers = false; it failed veto check with some other particle
        return false;
      }
    }
    switch (pidCutType) {
      case kRectangularCut:
        return selIdRectangularCut<pidMode>(track, nSigmaTPC, nSigmaTOF);
      case kCircularCut:
        return selIdCircularCut<pidMode>(track, nSigmaSquaredRad);
      case kEllipsoidalCut:
        return selIdEllipsoidalCut<pidMode>(track, nSigmaTPC, nSigmaTOF);
      default:
        return false;
    }
  }

  template <int pidMode, typename T>
  bool selPdependent(const T& track, int& IdMethod)
  {
    // Static cache inside function - initialized once on first call
    static const float thrPforTOF = getCfg<float>(cfgIdCut.pidConfigSetting, pidMode, kThrPforTOF);
    static const int idCutTypeLowP = getCfg<int>(cfgIdCut.pidConfigSetting, pidMode, kIdCutTypeLowP);
    static const float nSigmaTPCLowP = getCfg<float>(cfgIdCut.pidConfigSetting, pidMode, kNSigmaTPCLowP);
    static const float nSigmaTOFLowP = getCfg<float>(cfgIdCut.pidConfigSetting, pidMode, kNSigmaTOFLowP);
    static const float nSigmaRadLowP = getCfg<float>(cfgIdCut.pidConfigSetting, pidMode, kNSigmaRadLowP);
    static const int idCutTypeHighP = getCfg<int>(cfgIdCut.pidConfigSetting, pidMode, kIdCutTypeHighP);
    static const float nSigmaTPCHighP = getCfg<float>(cfgIdCut.pidConfigSetting, pidMode, kNSigmaTPCHighP);
    static const float nSigmaTOFHighP = getCfg<float>(cfgIdCut.pidConfigSetting, pidMode, kNSigmaTOFHighP);
    static const float nSigmaRadHighP = getCfg<float>(cfgIdCut.pidConfigSetting, pidMode, kNSigmaRadHighP);

    if (track.p() < thrPforTOF) {
      if (checkReliableTOF(track)) {
        if (idTPCTOF<pidMode>(track, idCutTypeLowP, nSigmaTPCLowP, nSigmaTOFLowP, nSigmaRadLowP)) {
          IdMethod = kTPCTOFidentified;
          return true;
        }
      } else {
        if (idTPC<pidMode>(track, nSigmaTPCLowP)) {
          IdMethod = kTPCidentified;
          return true;
        }
      }
    } else {
      if (checkReliableTOF(track)) {
        if (idTPCTOF<pidMode>(track, idCutTypeHighP, nSigmaTPCHighP, nSigmaTOFHighP, nSigmaRadHighP)) {
          IdMethod = kTPCTOFidentified;
          return true;
        }
      }
    }
    return false;
  }

  //_______________________________Identification Funtions Depending on the tpcInnerParam _______________________________
  // tpc Selections
  template <typename T>
  bool selPionTPCInnerParam(const T& track)
  {
    static const float nSigmaTPCLowP = getCfg<float>(cfgIdCut.pidConfigSetting, kPi, kNSigmaTPCLowP);
    static const float nSigmaTPCHighP = getCfg<float>(cfgIdCut.pidConfigSetting, kPi, kNSigmaTPCHighP);
    constexpr float kTPCInnerLow = 0.05f;
    constexpr float kTPCInnerHigh = 0.70f;

    if (vetoIdOthersTPC<kPi>(track)) {
      if (kTPCInnerLow <= track.tpcInnerParam() && track.tpcInnerParam() < kTPCInnerHigh && std::abs(track.tpcNSigmaPi()) < nSigmaTPCLowP) {
        return true;
      }
      if (kTPCInnerHigh <= track.tpcInnerParam() && std::abs(track.tpcNSigmaPi()) < nSigmaTPCHighP) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selKaonTPCInnerParam(const T& track)
  {
    static const float nSigmaTPCLowP = getCfg<float>(cfgIdCut.pidConfigSetting, kKa, kNSigmaTPCLowP);
    static const float nSigmaTPCHighP = getCfg<float>(cfgIdCut.pidConfigSetting, kKa, kNSigmaTPCHighP);
    constexpr float kTPCInnerLow = 0.05f;
    constexpr float kTPCInnerHigh = 0.70f;
    if (vetoIdOthersTPC<kKa>(track)) {
      if (kTPCInnerLow <= track.tpcInnerParam() && track.tpcInnerParam() < kTPCInnerHigh && std::abs(track.tpcNSigmaKa()) < nSigmaTPCLowP) {
        return true;
      }
      if (kTPCInnerHigh <= track.tpcInnerParam() && std::abs(track.tpcNSigmaKa()) < nSigmaTPCHighP) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selProtonTPCInnerParam(const T& track)
  {
    static const float nSigmaTPCLowP = getCfg<float>(cfgIdCut.pidConfigSetting, kPr, kNSigmaTPCLowP);
    static const float nSigmaTPCHighP = getCfg<float>(cfgIdCut.pidConfigSetting, kPr, kNSigmaTPCHighP);
    constexpr float kTPCInnerLow = 0.05f;
    constexpr float kTPCInnerHigh = 1.60f;

    if (vetoIdOthersTPC<kPr>(track)) {
      if (kTPCInnerLow <= track.tpcInnerParam() && track.tpcInnerParam() < kTPCInnerHigh && std::abs(track.tpcNSigmaPr()) < nSigmaTPCLowP) {
        return true;
      }
      if (kTPCInnerHigh <= track.tpcInnerParam() && std::abs(track.tpcNSigmaPr()) < nSigmaTPCHighP) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selDeuteronTPCInnerParam(const T& track)
  {
    constexpr float kTPCInnerLow = 0.05f;
    constexpr float kTPCInnerHigh = 1.80f;
    constexpr float kNSigmaLowP = 3.0f;
    constexpr float kNSigmaHighP = 2.0f;

    if (vetoIdOthersTPC<kDe>(track)) {
      if (kTPCInnerLow <= track.tpcInnerParam() && track.tpcInnerParam() < kTPCInnerHigh && std::abs(track.tpcNSigmaDe()) < kNSigmaLowP) {
        return true;
      }
      if (kTPCInnerHigh <= track.tpcInnerParam() && std::abs(track.tpcNSigmaDe()) < kNSigmaHighP) {
        return true;
      }
    }
    return false;
  }

  template <typename T>
  bool selElectronTPCInnerParam(const T& track)
  {
    constexpr float kNSigmaMax = 3.0f;
    if (track.tpcNSigmaEl() < kNSigmaMax && track.tpcNSigmaPi() > kNSigmaMax && track.tpcNSigmaKa() > kNSigmaMax && track.tpcNSigmaPr() > kNSigmaMax && track.tpcNSigmaDe() > kNSigmaMax) {
      return true;
    }
    return false;
  }
  //
  //_____________________________________________________TOF selection Functions _______________________________________________________________________
  // TOF Selections
  // Pion
  template <typename T>
  bool selPionTOF(const T& track)
  {
    // Constants to avoid magic numbers and repeated getCfg calls
    static constexpr float ThresholdP = 0.75f;
    static const float nSigmaTPCLowP = getCfg<float>(cfgIdCut.pidConfigSetting, kPi, kNSigmaTPCLowP);
    static const float nSigmaTOFLowP = getCfg<float>(cfgIdCut.pidConfigSetting, kPi, kNSigmaTOFLowP);
    static const float nSigmaTPCHighP = getCfg<float>(cfgIdCut.pidConfigSetting, kPi, kNSigmaTPCHighP);
    static const float nSigmaTOFHighP = getCfg<float>(cfgIdCut.pidConfigSetting, kPi, kNSigmaTOFHighP);
    if (vetoIdOthersTOF<kPi>(track)) {
      if (track.p() <= ThresholdP && std::abs(track.tpcNSigmaPi()) < nSigmaTPCLowP && std::abs(track.tofNSigmaPi()) < nSigmaTOFLowP) {
        return true;
      } else if (ThresholdP < track.p() // after p = 0.75, Pi and Ka lines of nSigma 3.0 will start intersecting
                 && std::abs(track.tpcNSigmaPi()) < nSigmaTPCHighP && std::abs(track.tofNSigmaPi()) < nSigmaTOFHighP) {
        return true;
      }
    }
    return false;
  }

  // Kaon
  template <typename T>
  bool selKaonTOF(const T& track)
  {
    static constexpr float ThresholdPLow = 0.75f; // π-K separation
    static constexpr float ThresholdPUp = 1.30f;  // K-p separation
    static const float nSigmaTPCLowP = getCfg<float>(cfgIdCut.pidConfigSetting, kKa, kNSigmaTPCLowP);
    static const float nSigmaTOFLowP = getCfg<float>(cfgIdCut.pidConfigSetting, kKa, kNSigmaTOFLowP);
    static const float nSigmaTPCHighP = getCfg<float>(cfgIdCut.pidConfigSetting, kKa, kNSigmaTPCHighP);
    static const float nSigmaTOFHighP = getCfg<float>(cfgIdCut.pidConfigSetting, kKa, kNSigmaTOFHighP);

    if (vetoIdOthersTOF<kKa>(track)) {
      if (track.p() <= ThresholdPLow && std::abs(track.tpcNSigmaKa()) < nSigmaTPCLowP && std::abs(track.tofNSigmaKa()) < nSigmaTOFLowP) {
        return true;
      }
      if (ThresholdPLow < track.p() && track.p() <= ThresholdPUp // after 0.75 Pi and Ka lines of nSigma 3.0 will start intersecting
          && std::abs(track.tpcNSigmaKa()) < nSigmaTPCLowP && std::abs(track.tofNSigmaKa()) < nSigmaTOFLowP) {
        return true;
      }
      if (ThresholdPUp < track.p() // after 1.30 Pr and Ka lines of nSigma 3.0 will start intersecting
          && std::abs(track.tpcNSigmaKa()) < nSigmaTPCHighP && std::abs(track.tofNSigmaKa()) < nSigmaTOFHighP) {
        return true;
      }
    }
    return false;
  }

  // Proton
  template <typename T>
  bool selProtonTOF(const T& track)
  {
    // Static config values (fetched once per template instantiation)
    static constexpr float ThresholdPLow = 1.30f; // Kaon-proton separation
    static constexpr float ThresholdPUp = 3.10f;  // Proton-deuteron separation
    static const float nSigmaTPCLowP = getCfg<float>(cfgIdCut.pidConfigSetting, kPr, kNSigmaTPCLowP);
    static const float nSigmaTOFLowP = getCfg<float>(cfgIdCut.pidConfigSetting, kPr, kNSigmaTOFLowP);
    static const float nSigmaTPCHighP = getCfg<float>(cfgIdCut.pidConfigSetting, kPr, kNSigmaTPCHighP);
    static const float nSigmaTOFHighP = getCfg<float>(cfgIdCut.pidConfigSetting, kPr, kNSigmaTOFHighP);

    if (vetoIdOthersTOF<kPr>(track)) {
      if (track.p() <= ThresholdPLow && std::abs(track.tpcNSigmaPr()) < nSigmaTPCLowP && std::abs(track.tofNSigmaPr()) < nSigmaTOFLowP) {
        return true;
      }
      if (ThresholdPLow < track.p() && track.p() <= ThresholdPUp                                            // after 1.30 Pr and Ka lines of nSigma 3.0 will start intersecting
          && std::abs(track.tpcNSigmaPr()) < nSigmaTPCLowP && std::abs(track.tofNSigmaPr()) < nSigmaTOFLowP // Some Deuteron contamination is still coming in p dependent cuts
      ) {
        return true;
      }
      if (ThresholdPUp < track.p() // after 3.10 Pr and De lines of nSigma 3.0 will start intersecting
          && std::abs(track.tpcNSigmaPr()) < nSigmaTPCHighP && std::abs(track.tofNSigmaPr()) < nSigmaTOFHighP) {
        return true;
      }
    }
    return false;
  }

  // Deuteron
  template <typename T>
  bool selDeuteronTOF(const T& track)
  {
    static constexpr float ThresholdP = 3.10f;
    static constexpr float TpcSigmaLowCut = 3.0f;
    static constexpr float TofSigmaLowCut = 3.0f;
    static constexpr float TpcSigmaHighCut = 2.0f;
    static constexpr float TofSigmaHighCut = 2.0f;

    if (vetoIdOthersTOF<kDe>(track)) {
      if (track.p() <= ThresholdP && std::abs(track.tpcNSigmaDe()) < TpcSigmaLowCut && std::abs(track.tofNSigmaDe()) < TofSigmaLowCut) {
        return true;
      }
      if (ThresholdP < track.p() // after 3.10 De and Pr lines of nSigma 3.0 will start intersecting
          && std::abs(track.tpcNSigmaDe()) < TpcSigmaHighCut && std::abs(track.tofNSigmaDe()) < TofSigmaHighCut) {
        return true;
      }
    }
    return false;
  }

  // Electron
  template <typename T>
  bool selElectronTOF(const T& track)
  {
    if (!vetoIdOthersTOF<kEl>(track)) {
      return false;
    }
    static const float kNSigmaRadHighP = getCfg<float>(cfgIdCut.pidConfigSetting, kEl, kNSigmaRadHighP);
    const float sumSq = track.tpcNSigmaEl() * track.tpcNSigmaEl() + track.tofNSigmaEl() * track.tofNSigmaEl();
    return sumSq < kNSigmaRadHighP;
  }

  //______________________________Identification Functions________________________________________________________________
  // Pion
  template <typename T>
  bool selPion(const T& track, int& IdMethod)
  {
    if (cfgIdCut.cfgId04DoPdependentId) {
      return selPdependent<kPi>(track, IdMethod);
    } else if (cfgIdCut.cfgId05DoTpcInnerParamId) {
      if (selPionTPCInnerParam(track)) {
        IdMethod = kTPCidentified;
        return true;
      } else if (track.hasTOF() && track.beta() > 0.0 && selPionTOF(track)) {
        IdMethod = kTOFidentified;
        return true;
      }
      return false;
    }
    return false;
  }

  // Kaon
  template <typename T>
  bool selKaon(const T& track, int& IdMethod)
  {
    if (cfgIdCut.cfgId04DoPdependentId) {
      return selPdependent<kKa>(track, IdMethod);
    } else if (cfgIdCut.cfgId05DoTpcInnerParamId) {
      if (selKaonTPCInnerParam(track)) {
        IdMethod = kTPCidentified;
        return true;
      } else if (track.hasTOF() && track.beta() > 0.0 && selKaonTOF(track)) {
        IdMethod = kTOFidentified;
        return true;
      }
      return false;
    }
    return false;
  }

  // Proton
  template <typename T>
  bool selProton(const T& track, int& IdMethod)
  {
    if (cfgIdCut.cfgId04DoPdependentId) {
      return selPdependent<kPr>(track, IdMethod);
    } else if (cfgIdCut.cfgId05DoTpcInnerParamId) {
      if (selProtonTPCInnerParam(track)) {
        IdMethod = kTPCidentified;
        return true;
      } else if (track.hasTOF() && track.beta() > 0.0 && selProtonTOF(track)) {
        IdMethod = kTOFidentified;
        return true;
      }
      return false;
    }
    return false;
  }

  // Electron
  template <typename T>
  bool selElectron(const T& track, int& IdMethod)
  {
    if (cfgIdCut.cfgId04DoPdependentId) {
      return selPdependent<kEl>(track, IdMethod);
    } else if (cfgIdCut.cfgId05DoTpcInnerParamId) {
      if (selElectronTPCInnerParam(track)) {
        IdMethod = kTPCidentified;
        return true;
      } else if (track.hasTOF() && track.beta() > 0.0 && selElectronTOF(track)) {
        IdMethod = kTOFidentified;
        return true;
      }
      return false;
    }
    return false;
  }

  // Muon
  template <typename T>
  bool selMuon(const T& track, int& IdMethod)
  {
    if (cfgIdCut.cfgId04DoPdependentId) {
      return selPdependent<kMu>(track, IdMethod);
    }
    return false;
  }

  // Deuteron
  template <typename T>
  bool selDeuteron(const T& track, int& IdMethod)
  {
    if (cfgIdCut.cfgId04DoPdependentId) {
      return selPdependent<kPr>(track, IdMethod);
    } else if (cfgIdCut.cfgId05DoTpcInnerParamId) {
      if (selDeuteronTPCInnerParam(track)) {
        IdMethod = kTPCidentified;
        return true;
      } else if (track.hasTOF() && track.beta() > 0.0 && selDeuteronTOF(track)) {
        IdMethod = kTOFidentified;
        return true;
      }
      return false;
    }
    return false;
  }

  template <int particleMode, typename T>
  inline bool selParticle(const T& track, int& idMethod)
  {
    idMethod = kUnidentified;
    if constexpr (particleMode == kPi) {
      return selPion(track, idMethod);
    }
    if constexpr (particleMode == kKa) {
      return selKaon(track, idMethod);
    }
    if constexpr (particleMode == kPr) {
      return selProton(track, idMethod);
    }
    if constexpr (particleMode == kEl) {
      return selElectron(track, idMethod);
    }
    if constexpr (particleMode == kMu) {
      return selMuon(track, idMethod);
    }
    return false;
  }

  //____________________________Other Functions and objects usd in Analysis___________________
  // Brian Kernighan’s Algorithm
  int countBits(uint8_t x)
  {
    int count = 0;
    while (x) {
      x &= (x - 1); // Clear the least significant bit set
      count++;
    }
    return count;
  }

  template <typename ReturnType>
  using GenericCaster = std::function<ReturnType(void*)>;

  template <typename ReturnType>
  std::vector<GenericCaster<ReturnType>> setupCasters(const std::vector<std::type_index>& types)
  {
    std::vector<GenericCaster<ReturnType>> casters;

    for (const auto& t : types) {
      if (t == typeid(int)) {
        casters.emplace_back([](void* ptr) { return static_cast<ReturnType>(*static_cast<int*>(ptr)); });
      } else if (t == typeid(float)) {
        casters.emplace_back([](void* ptr) { return static_cast<ReturnType>(*static_cast<float*>(ptr)); });
      } else if (t == typeid(double)) {
        casters.emplace_back([](void* ptr) { return static_cast<ReturnType>(*static_cast<double*>(ptr)); });
      } else {
        casters.emplace_back([](void*) {
          std::cerr << "Unsupported type\n";
          return ReturnType{};
        });
      }
    }

    return casters;
  }

  template <typename T>
  inline bool isSorted(const T& vec)
  {
    return std::is_sorted(vec.begin(), vec.end(), [](const MyTrackData& a, const MyTrackData& b) { return a.globalIndex < b.globalIndex; });
    // if a.globalIndex == b.globalIndex  then it will give false result, so remove duplicates before checking.  a.globalIndex <= b.globalIndex; was not working as expected.
  }

  template <typename T>
  inline void sortVector(T& vec)
  {
    std::sort(vec.begin(), vec.end(), [](const MyTrackData& a, const MyTrackData& b) { return a.globalIndex < b.globalIndex; });
  }

  template <typename T, typename H>
  void findRepeatedEntries(std::vector<H> ParticleList, T hist)
  {
    for (uint ii = 0; ii < ParticleList.size(); ii++) {
      int nCommonCount = 0; // checking the repeat number of track
      for (uint jj = 0; jj < ParticleList.size(); jj++) {
        if (ParticleList[jj].globalIndex == ParticleList[ii].globalIndex) {
          if (jj < ii) {
            break;
          } // break if it was already counted
          nCommonCount++; // To Calculate no of times the entry was repeated
        }
      }
      hist->Fill(nCommonCount);
    }
  }

  void removeDuplicates(auto& fullDauList)
  {
    auto last = std::unique(fullDauList.begin(), fullDauList.end(), [](const MyTrackData& a, const MyTrackData& b) { return a.globalIndex == b.globalIndex; });
    fullDauList.erase(last, fullDauList.end()); // std::unique only moves duplicates to end of the vector // last is the iterator position from where duplicate entries start
  }

  template <typename T>
  void printList(const T& vec)
  {
    int i = 0;
    for (const auto& v : vec) {
      LOG(info) << "DEBUG :: i = " << i << " :: GI = " << v.globalIndex;
      i++;
    }
  }

  template <typename T, typename U>
  int binarySearchAnyList(const T& ParticleList, const U& key)
  {
    if (ParticleList.empty())
      return -1;
    int low = 0;
    int high = ParticleList.size() - 1;
    while (low <= high) {
      int mid = low + (high - low) / 2;
      if (ParticleList[mid].globalIndex == key) {
        return mid;
      }
      if (ParticleList[mid].globalIndex < key)
        low = mid + 1;
      else
        high = mid - 1;
    }
    return -1; // If we reach here, then element was not present
  }

  template <bool checkTag, typename T>
  void executeSortPairDaughters(auto& posDauList, auto& negDauList, T hist)
  {
    sortVector(posDauList);
    sortVector(negDauList);
    if constexpr (checkTag) {
      if (!isSorted(posDauList)) {
        LOG(warning) << "DEBUG :: ERROR ERROR ERROR :: Usorted even after sorting";
        sortVector(posDauList);
        printList(posDauList);
      }
      if (!isSorted(negDauList)) {
        LOG(warning) << "DEBUG :: ERROR ERROR ERROR :: Usorted even after sorting";
        sortVector(negDauList);
        printList(negDauList);
      }
    }
    findRepeatedEntries(posDauList, hist);
    findRepeatedEntries(negDauList, hist);

    removeDuplicates(posDauList);
    removeDuplicates(negDauList);

    if (!isSorted(posDauList)) {
      LOG(fatal) << "DEBUG :: ERROR ERROR ERROR :: Usorted even after removing duplicates";
      sortVector(posDauList);
      printList(posDauList);
    }
    if (!isSorted(negDauList)) {
      LOG(fatal) << "DEBUG :: ERROR ERROR ERROR :: Usorted even after removing duplicates";
      sortVector(negDauList);
      printList(negDauList);
    }
  }

  template <bool isPrimVtxCndtList, typename T, typename U>
  void contaminatoinFlagCheck(const T& decayCndtDauList, const U& posTrack, const U& negTrack, auto& posTrackMotherFlag, auto& negTrackMotherFlag, auto& mSigma, auto& mpBit)
  {
    posTrackMotherFlag = 0;
    negTrackMotherFlag = 0;
    mSigma = doublePowersOf10[8];
    mpBit = 0;

    uint start = 0;
    uint end = 0;
    if constexpr (isPrimVtxCndtList) {
      start = kPrimTrkPhi1020;
      end = kPrimTrkRho770;
    } else {
      start = kV0TrkK0s;
      end = kV0TrkGamma;
    }

    for (uint i = start; i <= end; i++) {
      // Positive Daughter List
      int idx = binarySearchAnyList(decayCndtDauList[i][kPos], posTrack.globalIndex());
      if (idx != -1) {
        BITSET(posTrackMotherFlag, i);
        if (decayCndtDauList[i][kPos][idx].massSigma < mSigma) {
          mSigma = decayCndtDauList[i][kPos][idx].massSigma;
          mpBit = static_cast<uint8_t>(1 << i);
        }
        if (decayCndtDauList[i][kPos][idx].globalIndex != posTrack.globalIndex()) {
          LOG(fatal) << "DEBUG :: Wrong match triggered";
        }
      }
      // Negative Daughter List
      idx = binarySearchAnyList(decayCndtDauList[i][kNeg], negTrack.globalIndex());
      if (idx != -1) {
        BITSET(negTrackMotherFlag, i);
        if (decayCndtDauList[i][kNeg][idx].massSigma < mSigma) {
          mSigma = decayCndtDauList[i][kNeg][idx].massSigma;
          mpBit = static_cast<uint8_t>(1 << i);
        }
        if (decayCndtDauList[i][kNeg][idx].globalIndex != negTrack.globalIndex()) {
          LOG(fatal) << "DEBUG :: Wrong match triggered";
        }
      }
    }
  }

  //____________
  template <typename T>
  bool checkTrackSelection(const T& track, int& rejectionTag)
  {
    if (track.tpcNClsCrossedRows() < cfgTrackCuts.cfgTrk01TpcNClsCrossedRows) {
      rejectionTag = kFailTpcNClsCrossedRows;
      return false;
    }
    if (std::fabs(track.dcaXY()) > cfgTrackCuts.cfgTrk02dcaXY) {
      rejectionTag = kFailTrkdcaXY;
      return false;
    }
    if (std::fabs(track.dcaZ()) > cfgTrackCuts.cfgTrk03dcaZ) {
      rejectionTag = kFailTrkdcaZ;
      return false;
    }
    if (!track.isGlobalTrack()) {
      rejectionTag = kFailGlobalTrack;
      return false;
    }
    if (cfgTrackCuts.cfgTrk07DoVGselTrackCheck) {
      if (!selTrackForId(track)) {
        rejectionTag = kFailVGSelCheck;
        return false;
      }
    }
    return true;
  }

  template <typename T>
  bool selK0s(const T& v0)
  {
    static const float mLow = getCfg<float>(cfgV0ParticleCuts, kV0TrkK0s, kV0MLow);
    static const float mUp = getCfg<float>(cfgV0ParticleCuts, kV0TrkK0s, kV0MUp);
    static const float ptLow = getCfg<float>(cfgV0ParticleCuts, kV0TrkK0s, kV0LowPt);
    static const float ptHigh = getCfg<float>(cfgV0ParticleCuts, kV0TrkK0s, kV0HighPt);
    static const float yCut = getCfg<float>(cfgV0ParticleCuts, kV0TrkK0s, kV0Rapitidy);
    static const float armFactor = getCfg<float>(cfgV0ParticleCuts, kV0TrkK0s, kV0ARMcut1);

    if (mLow < v0.mK0Short() && v0.mK0Short() < mUp &&
        ptLow < v0.pt() && v0.pt() < ptHigh &&
        std::abs(v0.rapidity(MassK0Short)) < yCut &&
        v0.qtarm() > (armFactor * std::abs(v0.alpha()))) {
      return true;
    }
    return false;
  }

  template <typename T>
  bool selLambda(const T& v0)
  {
    static const float mLow = getCfg<float>(cfgV0ParticleCuts, kV0TrkLambda, kV0MLow);
    static const float mUp = getCfg<float>(cfgV0ParticleCuts, kV0TrkLambda, kV0MUp);

    if (mLow < v0.mLambda() && v0.mLambda() < mUp) {
      return true;
    } else {
      return false;
    }
  }

  template <typename T>
  bool selAntiLambda(const T& v0)
  {
    static const float mLow = getCfg<float>(cfgV0ParticleCuts, kV0TrkAntiLambda, kV0MLow);
    static const float mUp = getCfg<float>(cfgV0ParticleCuts, kV0TrkAntiLambda, kV0MUp);

    if (mLow < v0.mAntiLambda() && v0.mAntiLambda() < mUp) {
      return true;
    } else {
      return false;
    }
  }

  template <typename T>
  bool selGamma(const T& v0)
  {
    static const float mLow = getCfg<float>(cfgV0ParticleCuts, kV0TrkGamma, kV0MLow);
    static const float mUp = getCfg<float>(cfgV0ParticleCuts, kV0TrkGamma, kV0MUp);
    static const float armCut = getCfg<float>(cfgV0ParticleCuts, kV0TrkGamma, kV0ARMcut1);
    static const float alphaCut = getCfg<float>(cfgV0ParticleCuts, kV0TrkGamma, kV0ARMcut2);

    if (v0.mGamma() > mLow && v0.mGamma() < mUp &&
        v0.qtarm() < armCut &&
        std::abs(v0.alpha()) < alphaCut) {
      return true;
    } else {
      return false;
    }
  }

  template <int Mode, int fillMode, int pidMode, int signMode, int detMode, typename H, typename T>
  void fillTrackQAHistos(H histReg, const T& track, const double& particleMass, const float& effWeight)
  {
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h01_p"), track.p());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h02_pt"), track.pt());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h03_tpcInnerParam"), track.tpcInnerParam());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h04_tofExpMom"), track.tofExpMom());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h05_eta"), track.eta());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h06_phi"), track.phi());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h07_rapidity"), track.rapidity(particleMass));
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h08_isPVContributor"), track.isPVContributor());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h09_isGlobalTrack"), track.isGlobalTrack());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h10_sign"), track.sign());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h11_dcaXY"), track.dcaXY());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h12_dcaZ"), track.dcaZ());

    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h17_dcaXYwide"), track.dcaXY());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h18_dcaZwide"), track.dcaZ());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h20_pt_eta"), track.pt(), track.eta());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h20_pt_eta_effCorr"), track.pt(), track.eta(), effWeight);
  }

  template <int Mode, int fillMode, int pidMode, int signMode, int detMode, typename H, typename T>
  void fillTrackMomentumQAHistos(H histReg, const T& track)
  {
    // momemtum
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h21_p_pt"), track.p(), track.pt());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h22_p_tpcInnerParam"), track.p(), track.tpcInnerParam());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h23_p_tofExpMom"), track.p(), track.tofExpMom());

    // tpcSignal
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h24_p_tpcSignal"), track.p(), track.tpcSignal());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h25_tpcInnerParam_tpcSignal"), track.tpcInnerParam(), track.tpcSignal());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h26_tofExpMom_tpcSignal"), track.tofExpMom(), track.tpcSignal());

    // tofBeta
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h27_p_beta"), track.p(), track.beta());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h28_tpcInnerParam_beta"), track.tpcInnerParam(), track.beta());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h29_tofExpMom_beta"), track.tofExpMom(), track.beta());
  }

  template <int Mode, int fillMode, int pidMode, int signMode, int detMode, typename H, typename T>
  void fillTrackIdQAHistos(H& histReg, const T& track)
  {
    float tpcNSigmaVal = -999, tofNSigmaVal = -999;
    switch (pidMode) {
      case kPi:
        if (!cfgFill.cfgFill09PiQA)
          return;
        tpcNSigmaVal = track.tpcNSigmaPi();
        tofNSigmaVal = track.tofNSigmaPi();
        break;
      case kKa:
        if (!cfgFill.cfgFill10KaQA)
          return;
        tpcNSigmaVal = track.tpcNSigmaKa();
        tofNSigmaVal = track.tofNSigmaKa();
        break;
      case kPr:
        if (!cfgFill.cfgFill11PrQA)
          return;
        tpcNSigmaVal = track.tpcNSigmaPr();
        tofNSigmaVal = track.tofNSigmaPr();
        break;
      case kEl:
        if (!cfgFill.cfgFill12ElQA)
          return;
        tpcNSigmaVal = track.tpcNSigmaEl();
        tofNSigmaVal = track.tofNSigmaEl();
        break;
      case kDe:
        if (!cfgFill.cfgFill13DeQA)
          return;
        tpcNSigmaVal = track.tpcNSigmaDe();
        tofNSigmaVal = track.tofNSigmaDe();
        break;
      default:
        tpcNSigmaVal = -999, tofNSigmaVal = -999;
        break;
    }
    // NSigma
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h31_p_tpcNSigma"), track.p(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h32_pt_tpcNSigma"), track.pt(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h33_tpcInnerParam_tpcNSigma"), track.tpcInnerParam(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h34_tofExpMom_tpcNSigma"), track.tofExpMom(), tpcNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h35_p_tofNSigma"), track.p(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h36_pt_tofNSigma"), track.pt(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h37_tpcInnerParam_tofNSigma"), track.tpcInnerParam(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h38_tofExpMom_tofNSigma"), track.tofExpMom(), tofNSigmaVal);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h39_tpcNSigma_tofNSigma"), tpcNSigmaVal, tofNSigmaVal);
  }

  template <int Mode, int fillMode, int pidMode, int signMode, int detMode, typename H, typename T>
  void fillTrackDcaQAHistos(H& histReg, const T& track)
  {
    // DcaXY
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h41_p_dcaXY"), track.p(), track.dcaXY());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h42_pt_dcaXY"), track.pt(), track.dcaXY());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h43_tpcInnerParam_dcaXY"), track.tpcInnerParam(), track.dcaXY());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h44_tofExpMom_dcaXY"), track.tofExpMom(), track.dcaXY());

    // DcaZ
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h45_p_dcaZ"), track.p(), track.dcaZ());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h46_pt_dcaZ"), track.pt(), track.dcaZ());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h47_tpcInnerParam_dcaZ"), track.tpcInnerParam(), track.dcaZ());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST(DetDire[detMode]) + HIST("h48_tofExpMom_dcaZ"), track.tofExpMom(), track.dcaZ());
  }

  template <int Mode, int fillMode, int pidMode, int signMode, int detMode, typename H, typename T>
  void fillDaughterQA(H histReg, const T& track, const double& particleMass, const float& effWeight)
  {
    fillTrackQAHistos<Mode, fillMode, pidMode, signMode, detMode>(histReg, track, particleMass, effWeight);
    fillTrackMomentumQAHistos<Mode, fillMode, pidMode, signMode, detMode>(histReg, track);
    fillTrackIdQAHistos<Mode, fillMode, pidMode, signMode, detMode>(histReg, track);
    fillTrackDcaQAHistos<Mode, fillMode, pidMode, signMode, detMode>(histReg, track);
  }

  template <int Mode, typename H, typename T>
  void fillV0Histos(H histReg, const T& v0, const int& v0Tag, const int& v0DauCollisionIndexTag, const int& v0DauBCTag, const float& effWeight)
  {
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h01_K0s_Mass"), v0.mK0Short());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h02_Lambda_Mass"), v0.mLambda());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h03_AntiLambda_Mass"), v0.mAntiLambda());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h04_V0Gamma_Mass"), v0.mGamma());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h04_v0DaughterCollisionIndexTag"), v0DauCollisionIndexTag);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h04_v0DauBCTag"), v0DauBCTag);
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h05_V0Tag"), v0Tag);

    // Topological Cuts
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h06_dcapostopv"), v0.dcapostopv());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h07_dcanegtopv"), v0.dcanegtopv());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h08_dcaV0daughters"), v0.dcaV0daughters());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h09_v0cosPA"), v0.v0cosPA());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h10_v0radius"), v0.v0radius());

    // K0s-FullInformation
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h12_p"), v0.p());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h13_pt"), v0.pt());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h14_eta"), v0.eta());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h15_phi"), v0.phi());

    if constexpr (Mode == v0TablePostK0sCheck || Mode == v0TablePostK0sMassCut || Mode == v0TablePostK0sSelectionCut) {
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h11_mass"), v0.mK0Short());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h16_rapidity0_hypMass"), v0.rapidity(v0.mK0Short()));
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h16_rapidity1_invMass"), v0.rapidity(MassK0Short));
    } else if constexpr (Mode == v0TablePostLambdaCheck) {
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h11_mass"), v0.mLambda());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h16_rapidity0_hypMass"), v0.rapidity(v0.mLambda()));
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h16_rapidity1_invMass"), v0.rapidity(MassLambda0));
    } else if constexpr (Mode == v0TablePostAntiLambdaCheck) {
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h11_mass"), v0.mAntiLambda());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h16_rapidity0_hypMass"), v0.rapidity(v0.mAntiLambda()));
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h16_rapidity1_invMass"), v0.rapidity(MassLambda0Bar));
    } else if constexpr (Mode == v0TablePostGammaCheck) {
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h11_mass"), v0.mGamma());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h16_rapidity0_hypMass"), v0.rapidity(v0.mGamma()));
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h16_rapidity1_invMass"), v0.rapidity(MassGamma));
    } else if constexpr (Mode == v0TableFull) {
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h11_mass"), v0.mK0Short());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h11_mass"), v0.mLambda());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h11_mass"), v0.mAntiLambda());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h11_mass"), v0.mGamma());
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h16_rapidity0_hypMass"), v0.rapidity(0.0));
      histReg.fill(HIST(HistRegDire[Mode]) + HIST("h16_rapidity1_invMass"), v0.rapidity(0.0));
    }

    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h17_alpha"), v0.alpha());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h18_qtarm"), v0.qtarm());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h19_alpha_qtarm"), v0.alpha(), v0.qtarm());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h20_pt_eta"), v0.pt(), v0.eta());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST("h20_pt_eta_effCorr"), v0.pt(), v0.eta(), effWeight);
  }

  template <int Mode, int fillMode, int posPidMode, int negPidMode, typename H, typename T, typename U>
  void fillV0QA(H histReg, const T& v0, const U& posDaughterTrack, const U& negDaughterTrack, const int& ptEtaBinPosTrk, const int& ptEtaBinNegTrk, const int& v0Tag, const int& v0DauCollisionIndexTag, const int& v0DauBCTag, const int& posDauIdMethod, const int& negDauIdMethod, const float& v0EffWeight, const bool& fillMotherQA, const bool& fillDauTrackQA)
  {
    float posDauMass = -999., negDauMass = -999.;
    if constexpr (Mode == v0TableFull || Mode == v0TablePostK0sCheck || Mode == v0TablePostK0sMassCut ||
                  Mode == v0TablePostK0sSelectionCut || Mode == recoK0sPreSel || Mode == recoK0sPostSel) {
      posDauMass = MassPiPlus;
      negDauMass = MassPiPlus;
    } else if constexpr (Mode == v0TablePostLambdaCheck || Mode == recoLambdaPostSel) {
      posDauMass = MassProton;
      negDauMass = MassPiPlus;
    } else if constexpr (Mode == v0TablePostAntiLambdaCheck || Mode == recoAntiLambdaPostSel) {
      posDauMass = MassPiPlus;
      negDauMass = MassProtonBar;
    } else if constexpr (Mode == v0TablePostGammaCheck || Mode == recoGammaPostSel) {
      posDauMass = MassElectron;
      negDauMass = MassElectron;
    } else {
      static_assert(Mode == v0TableFull || Mode == v0TablePostK0sCheck || Mode == v0TablePostK0sMassCut ||
                      Mode == v0TablePostK0sSelectionCut || Mode == v0TablePostLambdaCheck || Mode == v0TablePostLambdaCheck ||
                      Mode == v0TablePostGammaCheck || Mode == recoK0sPreSel || Mode == recoK0sPostSel ||
                      Mode == recoLambdaPostSel || Mode == recoAntiLambdaPostSel || Mode == recoGammaPostSel,
                    "Unsupported particleMode");
    }

    if (fillMotherQA) {
      fillV0Histos<Mode>(histReg, v0, v0Tag, v0DauCollisionIndexTag, v0DauBCTag, v0EffWeight);
    }
    if (fillDauTrackQA) {
      float posDauEffWeight = hPtEtaForEffCorrection[posPidMode][kPos]->GetBinContent(ptEtaBinPosTrk);
      float negDauEffWeight = hPtEtaForEffCorrection[negPidMode][kNeg]->GetBinContent(ptEtaBinNegTrk);
      if (posDauIdMethod == kTPCidentified) {
        fillDaughterQA<Mode, fillMode, posPidMode, kPos, tpcId>(histReg, posDaughterTrack, posDauMass, posDauEffWeight);
      } else if (posDauIdMethod == kTPCTOFidentified) {
        fillDaughterQA<Mode, fillMode, posPidMode, kPos, tpctofId>(histReg, posDaughterTrack, posDauMass, posDauEffWeight);
      } else if (posDauIdMethod == kUnidentified) {
        fillDaughterQA<Mode, fillMode, posPidMode, kPos, NoId>(histReg, posDaughterTrack, posDauMass, posDauEffWeight);
      }
      // LOG(info)<<"DEBUG :: Point 03";
      if (negDauIdMethod == kTPCidentified) {
        fillDaughterQA<Mode, fillMode, negPidMode, kNeg, tpcId>(histReg, negDaughterTrack, negDauMass, negDauEffWeight);
      } else if (negDauIdMethod == kTPCTOFidentified) {
        fillDaughterQA<Mode, fillMode, negPidMode, kNeg, tpctofId>(histReg, negDaughterTrack, negDauMass, negDauEffWeight);
      } else if (negDauIdMethod == kUnidentified) {
        fillDaughterQA<Mode, fillMode, negPidMode, kNeg, NoId>(histReg, negDaughterTrack, negDauMass, negDauEffWeight);
      }
    }
  }

  template <typename T>
  int findV0Tag(const T& posDaughterTrack, const T& negDaughterTrack, int& posPiIdMethod, int& posPrIdMethod, int& posElIdMethod, int& negPiIdMethod, int& negPrIdMethod, int& negElIdMethod)
  {
    bool posIsPion = false;
    bool posIsProton = false;
    bool posIsElectron = false;
    bool negIsPion = false;
    bool negIsProton = false;
    bool negIsElectron = false;

    int v0TagValue = 0;

    // Check if positive track is pion or proton or Electron
    if (selPion(posDaughterTrack, posPiIdMethod))
      posIsPion = true; // Coming From K0s    -> PiPlus + PiMinus and AntiLambda -> PiPlus + AntiProton
    if (selProton(posDaughterTrack, posPrIdMethod))
      posIsProton = true; // Coming From Lambda -> proton + PiMinus
    if (selElectron(posDaughterTrack, posElIdMethod))
      posIsElectron = true; // Coming From Gamma -> ElPlus + ElMinus
    if (selPion(negDaughterTrack, negPiIdMethod))
      negIsPion = true; // Coming From K0s       -> PiPlus + PiMinus and Lambda -> proton + PiMinus
    if (selProton(negDaughterTrack, negPrIdMethod))
      negIsProton = true; // Coming From AntiLambda -> PiPlus + AntiProton
    if (selElectron(negDaughterTrack, negElIdMethod))
      negIsElectron = true; // Coming From Gamma -> ElPlus + ElMinus

    if (posIsPion && negIsPion) {
      BITSET(v0TagValue, BIT_IS_K0S);
    }
    if (posIsProton && negIsPion) {
      BITSET(v0TagValue, BIT_IS_LAMBDA);
    }
    if (posIsPion && negIsProton) {
      BITSET(v0TagValue, BIT_IS_ANTILAMBDA);
    }
    if (posIsElectron && negIsElectron) {
      BITSET(v0TagValue, BIT_IS_GAMMA);
    }
    return v0TagValue;
  }

  template <typename T, typename U>
  int findCollisionIndexTag(const T& v0, const U& posDaughterTrack, const U& negDaughterTrack)
  {
    int v0daughterCollisionIndexTag = 0;
    if (v0.collisionId() == posDaughterTrack.collisionId()) {
      BITSET(v0daughterCollisionIndexTag, BIT_POS_DAU_HAS_SAME_COLL);
    }
    if (v0.collisionId() == negDaughterTrack.collisionId()) {
      BITSET(v0daughterCollisionIndexTag, BIT_NEG_DAU_HAS_SAME_COLL);
    }
    if (posDaughterTrack.collisionId() == negDaughterTrack.collisionId()) {
      BITSET(v0daughterCollisionIndexTag, BIT_BOTH_DAU_HAS_SAME_COLL);
    }
    return v0daughterCollisionIndexTag;
  }

  template <typename collisionType, typename T, typename U>
  int findBCTag(const T& v0, const U& posDaughterTrack, const U& negDaughterTrack)
  {
    int v0daughterBCTag = 0;
    uint64_t v0GlobalBC = v0.template collision_as<collisionType>().template bc_as<aod::BCsWithTimestamps>().globalBC();
    uint64_t posDauGlobalBC = posDaughterTrack.template collision_as<collisionType>().template bc_as<aod::BCsWithTimestamps>().globalBC();
    uint64_t negDauGlobalBC = negDaughterTrack.template collision_as<collisionType>().template bc_as<aod::BCsWithTimestamps>().globalBC();

    if (v0GlobalBC == posDauGlobalBC) {
      BITSET(v0daughterBCTag, BIT_POS_DAU_HAS_SAME_COLL_BC);
    }
    if (v0GlobalBC == negDauGlobalBC) {
      BITSET(v0daughterBCTag, BIT_NEG_DAU_HAS_SAME_COLL_BC);
    }
    if (posDauGlobalBC == negDauGlobalBC) {
      BITSET(v0daughterBCTag, BIT_BOTH_DAU_HAS_SAME_COLL_BC);
    }
    return v0daughterBCTag;
  }

  template <int analysisType, typename collisionType, typename T, typename U>
  void executeV0loop(const T& posDaughterTrack, const T& negDaughterTrack, const U& v0,
                     auto& idMethodPi, auto& idMethodPr, auto& idMethodEl,
                     int& v0Tag, int& trueV0TagValue,
                     int& v0DauCollisionIndexTag, int& v0DauBCTag, auto& v0CndtDauList,
                     int& v0PtEtaBin, int& posDauPtEtaBin, int& negDauPtEtaBin,
                     float& v0K0sEffWeight, float& v0LambdaEffWeight, float& v0AntiLambdaEffWeight, float& v0GammaEffWeight,
                     const auto& v0DecayTrueMcTag)
  {
    idMethodPi[kPos] = kUnidentified;
    idMethodPr[kPos] = kUnidentified;
    idMethodEl[kPos] = kUnidentified;
    idMethodPi[kNeg] = kUnidentified;
    idMethodPr[kNeg] = kUnidentified;
    idMethodEl[kNeg] = kUnidentified;

    v0Tag = findV0Tag(posDaughterTrack, negDaughterTrack, idMethodPi[kPos], idMethodPr[kPos], idMethodEl[kPos], idMethodPi[kNeg], idMethodPr[kNeg], idMethodEl[kNeg]);
    v0DauCollisionIndexTag = findCollisionIndexTag(v0, posDaughterTrack, negDaughterTrack);
    v0DauBCTag = findBCTag<collisionType>(v0, posDaughterTrack, negDaughterTrack);

    static bool fillMotherQA = true;
    static bool fillDauTrackQA = true;
    trueV0TagValue = 0;
    if (cfgFill.cfgFill01V0TableFull) {
      fillV0QA<v0TableFull, kFillSimple, kPi, kPi>(recoV0sFull, v0, posDaughterTrack, negDaughterTrack, posDauPtEtaBin, negDauPtEtaBin, v0Tag, v0DauCollisionIndexTag, v0DauBCTag, idMethodPi[kPos], idMethodPi[kNeg], 1.0, fillMotherQA, fillDauTrackQA);
    }

    // cut on dynamic columns for v0 particles
    if (v0.v0cosPA() < v0settingCosPA)
      return; // in place of continue;
    if (v0.v0radius() < v0settingRadius)
      return; // in place of continue;

    // K0s Analysis
    if (BOOL_BITCHECK(v0Tag, BIT_IS_K0S) && isTrueMcMatch<analysisType>(v0DecayTrueMcTag, kV0TrkK0s)) {
      v0K0sEffWeight = hPtEtaForEffCorrection[kV0K0s][kPos]->GetBinContent(v0PtEtaBin);

      if (cfgFill.cfgFill02V0TablePostK0sCheck) {
        fillV0QA<v0TablePostK0sCheck, kFillSimple, kPi, kPi>(recoV0sPostK0sCheck, v0, posDaughterTrack, negDaughterTrack, posDauPtEtaBin, negDauPtEtaBin, v0Tag, v0DauCollisionIndexTag, v0DauBCTag, idMethodPi[kPos], idMethodPi[kNeg], v0K0sEffWeight, fillMotherQA, fillDauTrackQA);
      }
      // K0s mass cut
      if (cfgFill.cfgFill03v0TablePostK0sMassCut) {
        static const float k0sMassLow = getCfg<float>(cfgV0ParticleCuts, kV0TrkK0s, kV0MLow);
        static const float k0sMassHigh = getCfg<float>(cfgV0ParticleCuts, kV0TrkK0s, kV0MUp);
        if (k0sMassLow < v0.mK0Short() && v0.mK0Short() < k0sMassHigh) {
          fillV0QA<v0TablePostK0sMassCut, kFillSimple, kPi, kPi>(recoV0sPostMassCut, v0, posDaughterTrack, negDaughterTrack, posDauPtEtaBin, negDauPtEtaBin, v0Tag, v0DauCollisionIndexTag, v0DauBCTag, idMethodPi[kPos], idMethodPi[kNeg], v0K0sEffWeight, fillMotherQA, fillDauTrackQA);
        }
      }

      // Final K0s Selection.
      if (selK0s(v0)) {
        BITSET(trueV0TagValue, BIT_IS_K0S);
        v0CndtDauList[kV0TrkK0s][kPos].emplace_back(posDaughterTrack.globalIndex(), std::abs((v0.mK0Short() - MV0DecayCndt[kV0TrkK0s]) / WV0DecayCndt[kV0TrkK0s]));
        v0CndtDauList[kV0TrkK0s][kNeg].emplace_back(negDaughterTrack.globalIndex(), std::abs((v0.mK0Short() - MV0DecayCndt[kV0TrkK0s]) / WV0DecayCndt[kV0TrkK0s]));
        k0sTagIndexList.push_back(v0.globalIndex());
        if (cfgFill.cfgFill04V0TablePostK0sSelectionCut) {
          fillV0QA<v0TablePostK0sSelectionCut, kFillSimple, kPi, kPi>(recoV0sPostK0sSelectionCut, v0, posDaughterTrack, negDaughterTrack, posDauPtEtaBin, negDauPtEtaBin, v0Tag, v0DauCollisionIndexTag, v0DauBCTag, idMethodPi[kPos], idMethodPi[kNeg], v0K0sEffWeight, fillMotherQA, fillDauTrackQA);
        }
      }
    } // End of K0s block

    // Lambda Analysis
    if (BOOL_BITCHECK(v0Tag, BIT_IS_LAMBDA) && isTrueMcMatch<analysisType>(v0DecayTrueMcTag, kV0TrkLambda)) {
      if (selLambda(v0)) {
        BITSET(trueV0TagValue, BIT_IS_LAMBDA);
        v0LambdaEffWeight = hPtEtaForEffCorrection[kV0Lambda][kPos]->GetBinContent(v0PtEtaBin);

        v0CndtDauList[kV0TrkLambda][kPos].emplace_back(posDaughterTrack.globalIndex(), std::abs((v0.mLambda() - MV0DecayCndt[kV0TrkLambda]) / WV0DecayCndt[kV0TrkLambda]));
        v0CndtDauList[kV0TrkLambda][kNeg].emplace_back(negDaughterTrack.globalIndex(), std::abs((v0.mLambda() - MV0DecayCndt[kV0TrkLambda]) / WV0DecayCndt[kV0TrkLambda]));
        if (cfgFill.cfgFill05V0TablePostLambdaCheck) {
          fillV0QA<v0TablePostLambdaCheck, kFillSimple, kPr, kPi>(recoV0sPostLambdaCheck, v0, posDaughterTrack, negDaughterTrack, posDauPtEtaBin, negDauPtEtaBin, v0Tag, v0DauCollisionIndexTag, v0DauBCTag, idMethodPr[kPos], idMethodPi[kNeg], v0LambdaEffWeight, fillMotherQA, fillDauTrackQA);
        }
      }
    }

    // AntiLambda Analysis
    if (BOOL_BITCHECK(v0Tag, BIT_IS_ANTILAMBDA) && isTrueMcMatch<analysisType>(v0DecayTrueMcTag, kV0TrkAntiLambda)) {
      if (selAntiLambda(v0)) {
        BITSET(trueV0TagValue, BIT_IS_ANTILAMBDA);
        v0AntiLambdaEffWeight = hPtEtaForEffCorrection[kV0AntiLambda][kPos]->GetBinContent(v0PtEtaBin);

        v0CndtDauList[kV0TrkAntiLambda][kPos].emplace_back(posDaughterTrack.globalIndex(), std::abs((v0.mAntiLambda() - MV0DecayCndt[kV0TrkAntiLambda]) / WV0DecayCndt[kV0TrkAntiLambda]));
        v0CndtDauList[kV0TrkAntiLambda][kNeg].emplace_back(negDaughterTrack.globalIndex(), std::abs((v0.mAntiLambda() - MV0DecayCndt[kV0TrkAntiLambda]) / WV0DecayCndt[kV0TrkAntiLambda]));
        if (cfgFill.cfgFill06V0TablePostAntiLambdaCheck) {
          fillV0QA<v0TablePostAntiLambdaCheck, kFillSimple, kPi, kPr>(recoV0sPostAntiLambdaCheck, v0, posDaughterTrack, negDaughterTrack, posDauPtEtaBin, negDauPtEtaBin, v0Tag, v0DauCollisionIndexTag, v0DauBCTag, idMethodPi[kPos], idMethodPr[kNeg], v0AntiLambdaEffWeight, fillMotherQA, fillDauTrackQA);
        }
      }
    }

    // Gamma Analysis
    if (BOOL_BITCHECK(v0Tag, BIT_IS_GAMMA) && isTrueMcMatch<analysisType>(v0DecayTrueMcTag, kV0TrkGamma)) {
      if (selGamma(v0)) {
        BITSET(trueV0TagValue, BIT_IS_GAMMA);
        v0GammaEffWeight = hPtEtaForEffCorrection[kV0Gamma][kPos]->GetBinContent(v0PtEtaBin);

        v0CndtDauList[kV0TrkGamma][kPos].emplace_back(posDaughterTrack.globalIndex(), std::abs((v0.mGamma() - MV0DecayCndt[kV0TrkGamma]) / WV0DecayCndt[kV0TrkGamma]));
        v0CndtDauList[kV0TrkGamma][kNeg].emplace_back(negDaughterTrack.globalIndex(), std::abs((v0.mGamma() - MV0DecayCndt[kV0TrkGamma]) / WV0DecayCndt[kV0TrkGamma]));
        if (cfgFill.cfgFill07V0TablePostGammaCheck) {
          fillV0QA<v0TablePostGammaCheck, kFillSimple, kEl, kEl>(recoV0sPostGammaCheck, v0, posDaughterTrack, negDaughterTrack, posDauPtEtaBin, negDauPtEtaBin, v0Tag, v0DauCollisionIndexTag, v0DauBCTag, idMethodEl[kPos], idMethodEl[kNeg], v0GammaEffWeight, fillMotherQA, fillDauTrackQA);
        }
      }
    }
    recoV0sFull.fill(HIST(HistRegDire[v0TableFull]) + HIST("hTrueV0TagCount"), trueV0TagValue);
  }

  template <int particleMode, int posPidMode, int negPidMode, int fillMode, int32_t row, typename H, typename P, typename T> //, typename L>
  void checkV0Particle(H& histReg, const P& v0, const T& posTrack, const T& negTrack,
                       const int& v0PtEtaBin, const int& posDauPtEtaBin, const int& negDauPtEtaBin,
                       const int& v0Tag, const int& v0DauCollisionIndexTag, const int& v0DauBCTag,
                       const int& idMethodPosDau, const int& idMethodNegDau,
                       auto& iNTrkV0, auto& fNTrkV0, auto& effWeightSum,
                       const auto& centrality,
                       const auto& posTrackPrimVtxMotherFlag, const auto& negTrackPrimVtxMotherFlag, const auto& posTrackV0MotherFlag, const auto& negTrackV0MotherFlag, const int& requiredBit, const auto& mpBit)
  {
    static const bool checkPrimVtxContm = getCfg<bool>(cfgV0ParticleCuts, row, kV0CheckPrimVtxContm);
    static const bool checkV0DecayContm = getCfg<bool>(cfgV0ParticleCuts, row, kV0CheckV0DecayContm);
    static const bool checkMPMSigma = getCfg<bool>(cfgV0ParticleCuts, row, kV0CheckMPMSigma);
    static const bool fillPostSel = getCfg<bool>(cfgV0ParticleCuts, row, kV0FillPostSel);
    static const bool fillPostSelDau = getCfg<bool>(cfgV0ParticleCuts, row, kV0FillPostSelDau);

    float v0EffWeight = 0;
    if constexpr (particleMode == recoK0sPostSel) {
      v0EffWeight = hPtEtaForEffCorrection[kV0K0s][kPos]->GetBinContent(v0PtEtaBin);
    }
    if constexpr (particleMode == recoLambdaPostSel) {
      v0EffWeight = hPtEtaForEffCorrection[kV0Lambda][kPos]->GetBinContent(v0PtEtaBin);
    }
    if constexpr (particleMode == recoAntiLambdaPostSel) {
      v0EffWeight = hPtEtaForEffCorrection[kV0AntiLambda][kPos]->GetBinContent(v0PtEtaBin);
    }
    if constexpr (particleMode == recoGammaPostSel) {
      v0EffWeight = hPtEtaForEffCorrection[kV0Gamma][kPos]->GetBinContent(v0PtEtaBin);
    }

    bool countIt = true;
    bool fillMotherQA = false;
    bool fillDauTrackQA = false;

    if constexpr (fillMode == kFillPostSel) {
      // Since it is selected it must have the that track's bit
      if (!BOOL_BITCHECK(posTrackV0MotherFlag, row)) {
        LOG(error) << "pos Daughter Track tagged in postSel, but not found in preSel";
      }
      if (!BOOL_BITCHECK(negTrackV0MotherFlag, row)) {
        LOG(error) << "neg Daughter Track tagged in postSel, but not found in preSel";
      }

      // check other bits to know if dauhter is contamination from other particles.
      // If contamination is to be checked, then remove particles with contamination.
      if (checkPrimVtxContm) {
        if (BOOL_OTHERBITSON(posTrackPrimVtxMotherFlag, row)) {
          countIt = false;
        }
        if (BOOL_OTHERBITSON(negTrackPrimVtxMotherFlag, row)) {
          countIt = false;
        }
      }
      if (checkV0DecayContm) {
        if (BOOL_OTHERBITSON(posTrackV0MotherFlag, row)) {
          countIt = false;
        }
        if (BOOL_OTHERBITSON(negTrackV0MotherFlag, row)) {
          countIt = false;
        }
        // If it is the most probable mass case then count it.
        if (countIt == false) {
          // LOG(info)<<"DEBUG :: cout check :: rqBit = "<<requiredBit<<" :: row = "<<row<<" : "<<std::bitset<8>(1<<row)<<" :: bool_check "<<BOOL_BITCHECK(mpBit, row)<<" :: "<<getCfg<bool>(cfgPrimVtxParticleCuts, row, kPrimCheckMPMSigma)<<" :: "<<std::bitset<8>(mpBit);
          if (requiredBit == 1 && BOOL_BITCHECK(mpBit, row) && checkMPMSigma) {
            countIt = true;
            // LOG(info)<<"DEBUG :: counted at row = "<<row<<" :: "<<std::bitset<8>(1<<row);
          }
        }
      }

      if (!countIt) {
        iNTrkV0[kNeg] += 1;
        fNTrkV0[kNeg] += v0EffWeight;
        effWeightSum[kNeg] += v0EffWeight;
        return;
      }
    }

    fillV0QA<particleMode, kFillSimple, posPidMode, negPidMode>(histReg, v0, posTrack, negTrack, posDauPtEtaBin, negDauPtEtaBin, v0Tag, v0DauCollisionIndexTag, v0DauBCTag, idMethodPosDau, idMethodNegDau, v0EffWeight, fillMotherQA, fillDauTrackQA);
    iNTrkV0[kPos] += 1;
    fNTrkV0[kPos] += v0EffWeight; // efficiency correction applied here
    effWeightSum[kPos] += v0EffWeight;

    if constexpr (particleMode == recoK0sPostSel) {
      histReg.fill(HIST(HistRegDire[particleMode]) + HIST("mK0s_vs_cent"), centrality, v0.mK0Short());                      // centrality dependent mass
      histReg.fill(HIST(HistRegDire[particleMode]) + HIST("mK0s_vs_cent_effcorr"), centrality, v0.mK0Short(), v0EffWeight); // centrality dependent mass
      histReg.fill(HIST(HistRegDire[particleMode]) + HIST("mK0s_vs_cent_vs_Pt"), v0.mK0Short(), centrality, v0.pt());       // centrality and pt dependent mass
    }
    if constexpr (particleMode == recoLambdaPostSel) {
      histReg.fill(HIST(HistRegDire[particleMode]) + HIST("mLambda_vs_cent"), centrality, v0.mLambda());                      // centrality dependent mass
      histReg.fill(HIST(HistRegDire[particleMode]) + HIST("mLambda_vs_cent_effcorr"), centrality, v0.mLambda(), v0EffWeight); // centrality dependent mass
      histReg.fill(HIST(HistRegDire[particleMode]) + HIST("mLambda_vs_cent_vs_Pt"), v0.mLambda(), centrality, v0.pt());       // centrality and pt dependent mass
    }
    if constexpr (particleMode == recoAntiLambdaPostSel) {
      histReg.fill(HIST(HistRegDire[particleMode]) + HIST("mAntiLambda_vs_cent"), centrality, v0.mAntiLambda());                      // centrality dependent mass
      histReg.fill(HIST(HistRegDire[particleMode]) + HIST("mAntiLambda_vs_cent_effcorr"), centrality, v0.mAntiLambda(), v0EffWeight); // centrality dependent mass
      histReg.fill(HIST(HistRegDire[particleMode]) + HIST("mAntiLambda_vs_cent_vs_Pt"), v0.mAntiLambda(), centrality, v0.pt());       // centrality dependent mass
    }
    if constexpr (particleMode == recoGammaPostSel) {
      histReg.fill(HIST(HistRegDire[particleMode]) + HIST("mGamma_vs_cent"), centrality, v0.mAntiLambda());                      // centrality dependent mass
      histReg.fill(HIST(HistRegDire[particleMode]) + HIST("mGamma_vs_cent_effcorr"), centrality, v0.mAntiLambda(), v0EffWeight); // centrality dependent mass
      histReg.fill(HIST(HistRegDire[particleMode]) + HIST("mGamma_vs_cent_vs_Pt"), v0.mAntiLambda(), centrality, v0.pt());       // centrality dependent mass
    }
  }

  template <int analysisType, typename collisionType, int fillMode, typename U, typename T, typename W>
  void executeV0InCollisionloop(const U& v0, const T& posDaughterTrack, const T& negDaughterTrack,
                                W& idMethodPi, W& idMethodPr, W& idMethodEl,
                                int& v0Tag, int& v0DauCollisionIndexTag, int& v0DauBCTag, auto& iNTrk, auto& fNTrk, auto& effWeightSum, const float& centrality,
                                const auto& posTrackPrimVtxMotherFlag, const auto& negTrackPrimVtxMotherFlag, const auto& posTrackV0MotherFlag, const auto& negTrackV0MotherFlag, const int& requiredBit, const auto& mpBit,
                                const auto& v0DecayTrueMcTag)
  {
    idMethodPi[kPos] = kUnidentified;
    idMethodPr[kPos] = kUnidentified;
    idMethodEl[kPos] = kUnidentified;
    idMethodPi[kNeg] = kUnidentified;
    idMethodPr[kNeg] = kUnidentified;
    idMethodEl[kNeg] = kUnidentified;

    int v0PtEtaBin = hPtEtaForEffCorrection[kV0K0s][kPos]->FindBin(v0.pt(), v0.eta());
    int posDauPtEtaBin = hPtEtaForEffCorrection[kPi][kPos]->FindBin(posDaughterTrack.pt(), posDaughterTrack.eta());
    int negDauPtEtaBin = hPtEtaForEffCorrection[kPi][kNeg]->FindBin(negDaughterTrack.pt(), negDaughterTrack.eta());

    v0Tag = findV0Tag(posDaughterTrack, negDaughterTrack, idMethodPi[kPos], idMethodPr[kPos], idMethodEl[kPos], idMethodPi[kNeg], idMethodPr[kNeg], idMethodEl[kNeg]);
    v0DauCollisionIndexTag = findCollisionIndexTag(v0, posDaughterTrack, negDaughterTrack);
    v0DauBCTag = findBCTag<collisionType>(v0, posDaughterTrack, negDaughterTrack);

    if (cfgFill.cfgFill05RecoK0sPreSel) {
      fillV0QA<recoK0sPreSel, kFillSimple, kPi, kPi>(recoK0s, v0, posDaughterTrack, negDaughterTrack, posDauPtEtaBin, negDauPtEtaBin, v0Tag, v0DauCollisionIndexTag, v0DauBCTag, idMethodPi[kPos], idMethodPi[kNeg], 1.0, true, true);
    }

    // K0s Analysis
    if (BOOL_BITCHECK(v0Tag, BIT_IS_K0S) && selK0s(v0) && isTrueMcMatch<analysisType>(v0DecayTrueMcTag, kV0TrkK0s)) {
      checkV0Particle<recoK0sPostSel, kPi, kPi, fillMode, kV0TrkK0s>(recoK0s,
                                                                     v0, posDaughterTrack, negDaughterTrack,
                                                                     v0PtEtaBin, posDauPtEtaBin, negDauPtEtaBin,
                                                                     v0Tag, v0DauCollisionIndexTag, v0DauBCTag,
                                                                     idMethodPi[kPos], idMethodPi[kNeg],
                                                                     iNTrk[kV0K0s], fNTrk[kV0K0s], effWeightSum[kV0K0s], centrality,
                                                                     posTrackPrimVtxMotherFlag, negTrackPrimVtxMotherFlag, posTrackV0MotherFlag, negTrackV0MotherFlag,
                                                                     requiredBit, mpBit);
    } // End of K0s block

    if (BOOL_BITCHECK(v0Tag, BIT_IS_LAMBDA) && selLambda(v0) && isTrueMcMatch<analysisType>(v0DecayTrueMcTag, kV0TrkLambda)) {
      checkV0Particle<recoLambdaPostSel, kPr, kPi, fillMode, kV0TrkLambda>(listHandler, // recoLambda,
                                                                           v0, posDaughterTrack, negDaughterTrack,
                                                                           v0PtEtaBin, posDauPtEtaBin, negDauPtEtaBin,
                                                                           v0Tag, v0DauCollisionIndexTag, v0DauBCTag,
                                                                           idMethodPr[kPos], idMethodPi[kNeg],
                                                                           iNTrk[kV0Lambda], fNTrk[kV0Lambda], effWeightSum[kV0Lambda], centrality,
                                                                           posTrackPrimVtxMotherFlag, negTrackPrimVtxMotherFlag, posTrackV0MotherFlag, negTrackV0MotherFlag,
                                                                           requiredBit, mpBit);
    } // End of Lambda block

    if (BOOL_BITCHECK(v0Tag, BIT_IS_ANTILAMBDA) && selAntiLambda(v0) && isTrueMcMatch<analysisType>(v0DecayTrueMcTag, kV0TrkAntiLambda)) {
      checkV0Particle<recoAntiLambdaPostSel, kPi, kPr, fillMode, kV0TrkAntiLambda>(listHandler, // recoAntiLambda,
                                                                                   v0, posDaughterTrack, negDaughterTrack,
                                                                                   v0PtEtaBin, posDauPtEtaBin, negDauPtEtaBin,
                                                                                   v0Tag, v0DauCollisionIndexTag, v0DauBCTag,
                                                                                   idMethodPi[kPos], idMethodPr[kNeg],
                                                                                   iNTrk[kV0AntiLambda], fNTrk[kV0AntiLambda], effWeightSum[kV0AntiLambda], centrality,
                                                                                   posTrackPrimVtxMotherFlag, negTrackPrimVtxMotherFlag, posTrackV0MotherFlag, negTrackV0MotherFlag,
                                                                                   requiredBit, mpBit);
    } // End of AntiLambda block

    if (BOOL_BITCHECK(v0Tag, BIT_IS_GAMMA) && selGamma(v0) && isTrueMcMatch<analysisType>(v0DecayTrueMcTag, kV0TrkGamma)) {
      checkV0Particle<recoGammaPostSel, kEl, kEl, fillMode, kV0TrkGamma>(listHandler, // recoGamma,
                                                                         v0, posDaughterTrack, negDaughterTrack,
                                                                         v0PtEtaBin, posDauPtEtaBin, negDauPtEtaBin,
                                                                         v0Tag, v0DauCollisionIndexTag, v0DauBCTag,
                                                                         idMethodEl[kPos], idMethodEl[kNeg],
                                                                         iNTrk[kV0Gamma], fNTrk[kV0Gamma], effWeightSum[kV0Gamma], centrality,
                                                                         posTrackPrimVtxMotherFlag, negTrackPrimVtxMotherFlag, posTrackV0MotherFlag, negTrackV0MotherFlag,
                                                                         requiredBit, mpBit);
    } // End of Gamma block
  }

  template <int Mode, int fillMode, typename H, typename T>
  void fillPrimVtxParticleQAHistos(H histReg, const T& mother, const float& effWeight)
  {
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST("h11_mass"), mother.M());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST("h12_p"), mother.P());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST("h13_pt"), mother.Pt());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST("h14_eta"), mother.Eta());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST("h15_phi"), RecoDecay::constrainAngle(mother.Phi(), 0.f, 1U));
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST("h16_rapidity1_invMass"), mother.Rapidity());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST("h20_pt_eta"), mother.Pt(), mother.Eta());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(FillModeDire[fillMode]) + HIST("h20_pt_eta_effCorr"), mother.Pt(), mother.Eta(), effWeight);
  };

  template <int particleMode, int posPidMode, int negPidMode, int fillMode, typename H, typename P, typename T>
  void fillPrimVtxCndtPair(H& histReg, const P& mother, const T& posTrack, const T& negTrack, const int& ptEtaBinPosTrk, const int& ptEtaBinNegTrk, const int& idMethodPosDau, const int& idMethodNegDau, int& iNTrkPrim, float& fNTrkPrim, float& effWeightSum, const bool& fillMotherQA, const bool& fillDauTrackQA)
  {
    float posDauMass = -999., negDauMass = -999.;
    if constexpr (particleMode == Phi1020) {
      posDauMass = MassKPlus;
      negDauMass = MassKPlus;
    } else if constexpr (particleMode == JPsiToEE) {
      posDauMass = MassElectron;
      negDauMass = MassElectron;
    } else if constexpr (particleMode == JPsiToMuMu) {
      posDauMass = MassMuonMinus;
      negDauMass = MassMuonMinus;
    } else if constexpr (particleMode == KStar892) {
      posDauMass = MassKPlus;
      negDauMass = MassPiPlus;
    } else if constexpr (particleMode == KStar892Bar) {
      posDauMass = MassPiPlus;
      negDauMass = MassKPlus;
    } else if constexpr (particleMode == Rho770) {
      posDauMass = MassPiPlus;
      negDauMass = MassPiPlus;
    } else {
      static_assert(particleMode == Phi1020 || particleMode == JPsiToEE || particleMode == JPsiToMuMu ||
                      particleMode == KStar892 || particleMode == KStar892Bar || particleMode == Rho770,
                    "Unsupported particleMode");
    }

    float effWeight = 1.0; // get efficiency for the PrimVtx particles in future.
    iNTrkPrim += 1;
    fNTrkPrim += effWeight;
    effWeightSum += effWeight;
    if (fillMotherQA) {
      fillPrimVtxParticleQAHistos<particleMode, fillMode>(histReg, mother, effWeight);
    }
    if (fillDauTrackQA) {
      float posDauEffWeight = hPtEtaForEffCorrection[posPidMode][kPos]->GetBinContent(ptEtaBinPosTrk);
      float negDauEffWeight = hPtEtaForEffCorrection[negPidMode][kNeg]->GetBinContent(ptEtaBinNegTrk);

      if (idMethodPosDau == kTPCidentified) {
        fillDaughterQA<particleMode, fillMode, posPidMode, kPos, tpcId>(histReg, posTrack, posDauMass, posDauEffWeight);
      } else if (idMethodPosDau == kTPCTOFidentified) {
        fillDaughterQA<particleMode, fillMode, posPidMode, kPos, tpctofId>(histReg, posTrack, posDauMass, posDauEffWeight);
      } else if (idMethodPosDau == kUnidentified) {
        fillDaughterQA<particleMode, fillMode, posPidMode, kPos, NoId>(histReg, posTrack, posDauMass, posDauEffWeight);
      }

      if (idMethodNegDau == kTPCidentified) {
        fillDaughterQA<particleMode, fillMode, negPidMode, kNeg, tpcId>(histReg, negTrack, negDauMass, negDauEffWeight);
      } else if (idMethodNegDau == kTPCTOFidentified) {
        fillDaughterQA<particleMode, fillMode, negPidMode, kNeg, tpctofId>(histReg, negTrack, negDauMass, negDauEffWeight);
      } else if (idMethodNegDau == kUnidentified) {
        fillDaughterQA<particleMode, fillMode, negPidMode, kNeg, NoId>(histReg, negTrack, negDauMass, negDauEffWeight);
      }
    }
  }

  using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>>;
  template <int particleMode, int posPidMode, int negPidMode, int fillMode, int32_t row, typename H, typename P, typename T, typename L>
  void checkPrimVtxParticle(H& histReg, const P& mother, const T& posTrack, const T& negTrack,
                            L& posDauList, L& negDauList, int& primVtxCndtTag,
                            const int ptEtaBinPosTrk, const int ptEtaBinNegTrk,
                            const int idMethodPosDau, const int idMethodNegDau,
                            auto& iNTrkPrim, auto& fNTrkPrim, auto& effWeightSum,
                            const auto posTrackPrimVtxMotherFlag, const auto negTrackPrimVtxMotherFlag, const auto posTrackV0MotherFlag, const auto negTrackV0MotherFlag, const int requiredBit, const auto mpBit)
  {

    // Cache config values once per template instantiation
    static const float massLow = getCfg<float>(cfgPrimVtxParticleCuts, row, kPrimMassLow);
    static const float massHigh = getCfg<float>(cfgPrimVtxParticleCuts, row, kPrimMassUp);

    static const bool checkPrimVtxContm = getCfg<bool>(cfgPrimVtxParticleCuts, row, kPrimCheckPrimVtxContm);
    static const bool checkV0DecayContm = getCfg<bool>(cfgPrimVtxParticleCuts, row, kPrimCheckV0DecayContm);
    static const bool checkMPMSigma = getCfg<bool>(cfgPrimVtxParticleCuts, row, kPrimCheckMPMSigma);

    static const bool fillPreSel = getCfg<bool>(cfgPrimVtxParticleCuts, row, kPrimFillPreSel);
    static const bool fillPreSelDau = getCfg<bool>(cfgPrimVtxParticleCuts, row, kPrimFillPreSelDau);
    static const bool fillPostSel = getCfg<bool>(cfgPrimVtxParticleCuts, row, kPrimFillPostSel);
    static const bool fillPostSelDau = getCfg<bool>(cfgPrimVtxParticleCuts, row, kPrimFillPostSelDau);

    static LorentzVector vecMother;
    uint16_t flagBit = 0;

    // Mass window check & vector assignment
    if constexpr (particleMode == Phi1020) {
      if (!(massLow < mother.mPhi1020() && mother.mPhi1020() < massHigh))
        return;
      BITSET(primVtxCndtTag, BIT_IS_PHI_1020);
      flagBit = 1 << BIT_IS_PHI_1020;
      vecMother = LorentzVector(mother.px(), mother.py(), mother.pz(), mother.mPhi1020());
    } else if constexpr (particleMode == JPsiToEE) {
      if (!(massLow < mother.mJPsiToEE() && mother.mJPsiToEE() < massHigh))
        return;
      BITSET(primVtxCndtTag, BIT_IS_JPSI_TO_EE);
      flagBit = 1 << BIT_IS_JPSI_TO_EE;
      vecMother = LorentzVector(mother.px(), mother.py(), mother.pz(), mother.mJPsiToEE());
    } else if constexpr (particleMode == JPsiToMuMu) {
      if (!(massLow < mother.mJPsiToMuMu() && mother.mJPsiToMuMu() < massHigh))
        return;
      BITSET(primVtxCndtTag, BIT_IS_JPSI_TO_MUMU);
      flagBit = 1 << BIT_IS_JPSI_TO_MUMU;
      vecMother = LorentzVector(mother.px(), mother.py(), mother.pz(), mother.mJPsiToMuMu());
    } else if constexpr (particleMode == KStar892) {
      if (!(massLow < mother.mKStar892() && mother.mKStar892() < massHigh))
        return;
      BITSET(primVtxCndtTag, BIT_IS_KSTAR_892);
      flagBit = 1 << BIT_IS_KSTAR_892;
      vecMother = LorentzVector(mother.px(), mother.py(), mother.pz(), mother.mKStar892());
    } else if constexpr (particleMode == KStar892Bar) {
      if (!(massLow < mother.mKStar892Bar() && mother.mKStar892Bar() < massHigh))
        return;
      BITSET(primVtxCndtTag, BIT_IS_KSTAR_892_BAR);
      flagBit = 1 << BIT_IS_KSTAR_892_BAR;
      vecMother = LorentzVector(mother.px(), mother.py(), mother.pz(), mother.mKStar892Bar());
    } else if constexpr (particleMode == Rho770) {
      if (!(massLow < mother.mRho770() && mother.mRho770() < massHigh))
        return;
      BITSET(primVtxCndtTag, BIT_IS_RHO_770);
      flagBit = 1 << BIT_IS_RHO_770;
      vecMother = LorentzVector(mother.px(), mother.py(), mother.pz(), mother.mRho770());
    } else {
      static_assert(particleMode == Phi1020 || particleMode == JPsiToEE || particleMode == JPsiToMuMu ||
                      particleMode == KStar892 || particleMode == KStar892Bar || particleMode == Rho770,
                    "Unsupported particleMode");
    }

    bool countIt = true;
    bool fillMotherQA = false;
    bool fillDauTrackQA = false;

    if constexpr (fillMode == kFillPreSel) {
      posDauList.emplace_back(posTrack.globalIndex(), std::abs((vecMother.M() - MPrimVtxCndt[row]) / WPrimVtxCndt[row]));
      negDauList.emplace_back(negTrack.globalIndex(), std::abs((vecMother.M() - MPrimVtxCndt[row]) / WPrimVtxCndt[row]));

      fillMotherQA = fillPreSel;
      fillDauTrackQA = fillPreSelDau;
    } else if constexpr (fillMode == kFillPostSel) {
      // Since it is selected it must have the that track's bit
      if (!BOOL_BITCHECK(posTrackPrimVtxMotherFlag, row)) {
        LOG(error) << "pos Daughter Track tagged in postSel, but not found in preSel";
      }
      if (!BOOL_BITCHECK(negTrackPrimVtxMotherFlag, row)) {
        LOG(error) << "neg Daughter Track tagged in postSel, but not found in preSel";
      }

      // check other bits to know if dauhter is contamination from other particles.
      // If contamination is to be check, remove particles with contamination.
      if (checkPrimVtxContm) {
        if (BOOL_OTHERBITSON(posTrackPrimVtxMotherFlag, row) || BOOL_OTHERBITSON(negTrackPrimVtxMotherFlag, row)) {
          countIt = false;
        }
        // If it is the most probable mass case then count it.
        if (countIt == false) {
          if (requiredBit == 0 && BOOL_BITCHECK(mpBit, row) && checkMPMSigma) {
            countIt = true;
          }
        }
      }
      if (checkV0DecayContm) {
        if (BOOL_OTHERBITSON(posTrackV0MotherFlag, row) || BOOL_OTHERBITSON(negTrackV0MotherFlag, row)) {
          countIt = false;
        }
      }
      if (!countIt) {
        iNTrkPrim[kNeg] += 1;
        fNTrkPrim[kNeg] += 1;
        effWeightSum[kNeg] += 1;
        return;
      }
      fillMotherQA = fillPostSel;
      fillDauTrackQA = fillPostSelDau;
    }
    fillPrimVtxCndtPair<particleMode, posPidMode, negPidMode, fillMode>(histReg,
                                                                        vecMother, posTrack, negTrack, ptEtaBinPosTrk, ptEtaBinNegTrk,
                                                                        idMethodPosDau, idMethodNegDau,
                                                                        iNTrkPrim[kPos], fNTrkPrim[kPos], effWeightSum[kPos], fillMotherQA, fillDauTrackQA);
  }

  template <int analysisType, int fillMode, typename P, typename T, typename D, typename B, typename C, typename L>
  void executePrimVtxCndtInCollisionloop(const P& mother, const T& posTrack, const T& negTrack,
                                         const int& ptEtaBinPosTrk, const int& ptEtaBinNegTrk, int& primVtxCndtTag, L& primVtxCndtDauList, D& idMethodSignTrk,
                                         C& checkTrackId, B& posDauIs, B& negDauIs,
                                         auto& iNTrkPrim, auto& fNTrkPrim, auto& effWeightSum,
                                         const auto posTrackPrimVtxMotherFlag, const auto negTrackPrimVtxMotherFlag, const auto posTrackV0MotherFlag, const auto negTrackV0MotherFlag, const int requiredBit, const auto mpBit,
                                         const bool doPhi1020, const bool doJPsiToEE, const bool doJPsiToMuMu, const bool doKStar892, const bool doKStar892Bar, const bool doRho770,
                                         const auto primVtxCndtMcTag)
  {
    for (int i = kPi; i <= kMu; i++) {
      posDauIs[i] = false;
      negDauIs[i] = false;
    }

    // do Identification of only what is needed
    if (checkTrackId[kPi][kPos]) {
      if (selParticle<kPi>(posTrack, idMethodSignTrk[kPi][kPos]))
        posDauIs[kPi] = true;
    }
    if (checkTrackId[kPi][kNeg]) {
      if (selParticle<kPi>(negTrack, idMethodSignTrk[kPi][kNeg]))
        negDauIs[kPi] = true;
    }
    if (checkTrackId[kKa][kPos]) {
      if (selParticle<kKa>(posTrack, idMethodSignTrk[kKa][kPos]))
        posDauIs[kKa] = true;
    }
    if (checkTrackId[kKa][kNeg]) {
      if (selParticle<kKa>(negTrack, idMethodSignTrk[kKa][kNeg]))
        negDauIs[kKa] = true;
    }
    if (checkTrackId[kEl][kPos]) {
      if (selParticle<kEl>(posTrack, idMethodSignTrk[kEl][kPos]))
        posDauIs[kEl] = true;
    }
    if (checkTrackId[kEl][kNeg]) {
      if (selParticle<kEl>(negTrack, idMethodSignTrk[kEl][kNeg]))
        negDauIs[kEl] = true;
    }
    if (checkTrackId[kMu][kPos]) {
      if (selParticle<kMu>(posTrack, idMethodSignTrk[kMu][kPos]))
        posDauIs[kMu] = true;
    }
    if (checkTrackId[kMu][kNeg]) {
      if (selParticle<kMu>(negTrack, idMethodSignTrk[kMu][kNeg]))
        negDauIs[kMu] = true;
    }

    // phi(1020) -> K+ + K-
    if (doPhi1020 && posDauIs[kKa] && negDauIs[kKa] && isTrueMcMatch<analysisType>(primVtxCndtMcTag, kPrimPhi1020)) {
      checkPrimVtxParticle<Phi1020, kKa, kKa, fillMode, kPrimTrkPhi1020>(primVtxTLHs[kPrimTrkPhi1020],
                                                                         mother, posTrack, negTrack, primVtxCndtDauList[kPrimTrkPhi1020][kPos], primVtxCndtDauList[kPrimTrkPhi1020][kNeg], primVtxCndtTag,
                                                                         ptEtaBinPosTrk, ptEtaBinNegTrk,
                                                                         idMethodSignTrk[kKa][kPos], idMethodSignTrk[kKa][kNeg],
                                                                         iNTrkPrim[kPrimPhi1020], fNTrkPrim[kPrimPhi1020], effWeightSum[kPrimPhi1020],
                                                                         posTrackPrimVtxMotherFlag, negTrackPrimVtxMotherFlag, posTrackV0MotherFlag, negTrackV0MotherFlag,
                                                                         requiredBit, mpBit);
    }

    // J/ψ       -> e+/e-
    if (doJPsiToEE && posDauIs[kEl] && negDauIs[kEl] && isTrueMcMatch<analysisType>(primVtxCndtMcTag, kPrimJPsiToEE)) {
      checkPrimVtxParticle<JPsiToEE, kEl, kEl, fillMode, kPrimTrkJPsiToEE>(primVtxTLHs[kPrimTrkJPsiToEE],
                                                                           mother, posTrack, negTrack, primVtxCndtDauList[kPrimTrkJPsiToEE][kPos], primVtxCndtDauList[kPrimTrkJPsiToEE][kNeg], primVtxCndtTag,
                                                                           ptEtaBinPosTrk, ptEtaBinNegTrk,
                                                                           idMethodSignTrk[kEl][kPos], idMethodSignTrk[kEl][kNeg],
                                                                           iNTrkPrim[kPrimJPsiToEE], fNTrkPrim[kPrimJPsiToEE], effWeightSum[kPrimJPsiToEE],
                                                                           posTrackPrimVtxMotherFlag, negTrackPrimVtxMotherFlag, posTrackV0MotherFlag, negTrackV0MotherFlag,
                                                                           requiredBit, mpBit);
    }

    // J/ψ       -> mu+/mu-
    if (doJPsiToMuMu && posDauIs[kMu] && negDauIs[kMu] && isTrueMcMatch<analysisType>(primVtxCndtMcTag, kPrimJPsiToMuMu)) {
      checkPrimVtxParticle<JPsiToMuMu, kMu, kMu, fillMode, kPrimTrkJPsiToMuMu>(primVtxTLHs[kPrimTrkJPsiToMuMu],
                                                                               mother, posTrack, negTrack, primVtxCndtDauList[kPrimTrkJPsiToMuMu][kPos], primVtxCndtDauList[kPrimTrkJPsiToMuMu][kNeg], primVtxCndtTag,
                                                                               ptEtaBinPosTrk, ptEtaBinNegTrk,
                                                                               idMethodSignTrk[kMu][kPos], idMethodSignTrk[kMu][kNeg],
                                                                               iNTrkPrim[kPrimJPsiToMuMu], fNTrkPrim[kPrimJPsiToMuMu], effWeightSum[kPrimJPsiToMuMu],
                                                                               posTrackPrimVtxMotherFlag, negTrackPrimVtxMotherFlag, posTrackV0MotherFlag, negTrackV0MotherFlag,
                                                                               requiredBit, mpBit);
    }

    // K*(892)^0 → K⁺ + π⁻
    if (doKStar892 && posDauIs[kKa] && negDauIs[kPi] && isTrueMcMatch<analysisType>(primVtxCndtMcTag, kPrimKStar892)) {
      checkPrimVtxParticle<KStar892, kKa, kPi, fillMode, kPrimTrkKStar892>(primVtxTLHs[kPrimTrkKStar892],
                                                                           mother, posTrack, negTrack, primVtxCndtDauList[kPrimTrkKStar892][kPos], primVtxCndtDauList[kPrimTrkKStar892][kNeg], primVtxCndtTag,
                                                                           ptEtaBinPosTrk, ptEtaBinNegTrk,
                                                                           idMethodSignTrk[kKa][kPos], idMethodSignTrk[kPi][kNeg],
                                                                           iNTrkPrim[kPrimKStar892], fNTrkPrim[kPrimKStar892], effWeightSum[kPrimKStar892],
                                                                           posTrackPrimVtxMotherFlag, negTrackPrimVtxMotherFlag, posTrackV0MotherFlag, negTrackV0MotherFlag,
                                                                           requiredBit, mpBit);
    }

    // K*(892)^0-bar → K⁻ + π⁺
    if (doKStar892Bar && posDauIs[kPi] && negDauIs[kKa] && isTrueMcMatch<analysisType>(primVtxCndtMcTag, kPrimKStar892Bar)) {
      checkPrimVtxParticle<KStar892Bar, kPi, kKa, fillMode, kPrimTrkKStar892Bar>(primVtxTLHs[kPrimTrkKStar892Bar],
                                                                                 mother, posTrack, negTrack, primVtxCndtDauList[kPrimTrkKStar892Bar][kPos], primVtxCndtDauList[kPrimTrkKStar892Bar][kNeg], primVtxCndtTag,
                                                                                 ptEtaBinPosTrk, ptEtaBinNegTrk,
                                                                                 idMethodSignTrk[kPi][kPos], idMethodSignTrk[kKa][kNeg],
                                                                                 iNTrkPrim[kPrimKStar892Bar], fNTrkPrim[kPrimKStar892Bar], effWeightSum[kPrimKStar892Bar],
                                                                                 posTrackPrimVtxMotherFlag, negTrackPrimVtxMotherFlag, posTrackV0MotherFlag, negTrackV0MotherFlag,
                                                                                 requiredBit, mpBit);
    }

    // ρ(770)    -> pi+ + pi-
    if (doRho770 && posDauIs[kPi] && negDauIs[kPi] && isTrueMcMatch<analysisType>(primVtxCndtMcTag, kPrimRho770)) {
      checkPrimVtxParticle<Rho770, kPi, kPi, fillMode, kPrimTrkRho770>(primVtxTLHs[kPrimTrkRho770],
                                                                       mother, posTrack, negTrack, primVtxCndtDauList[kPrimTrkRho770][kPos], primVtxCndtDauList[kPrimTrkRho770][kNeg], primVtxCndtTag,
                                                                       ptEtaBinPosTrk, ptEtaBinNegTrk,
                                                                       idMethodSignTrk[kPi][kPos], idMethodSignTrk[kPi][kNeg],
                                                                       iNTrkPrim[kPrimRho770], fNTrkPrim[kPrimRho770], effWeightSum[kPrimRho770],
                                                                       posTrackPrimVtxMotherFlag, negTrackPrimVtxMotherFlag, posTrackV0MotherFlag, negTrackV0MotherFlag,
                                                                       requiredBit, mpBit);
    }
  }

  void getTrackDecayInfoBit(const auto& track, auto& decayDauTagBit, const auto& v0CndtDauList, const auto& primVtxCndtDauList,
                            bool doV0K0s, bool doV0Lambda, bool doV0AntiLambda, bool doV0Gamma,
                            bool doPhi1020, bool doJPsiToEE, bool doJPsiToMuMu, bool doKStar892, bool doKStar892Bar, bool doRho770)
  {
    decayDauTagBit = 0;
    const int chargeIndex = (track.signed1Pt() > 0) ? kPos : (track.signed1Pt() < 0) ? kNeg
                                                                                     : -1;
    if (chargeIndex == -1)
      LOG(fatal) << "DEBUG :: Unsigned track found";

    const auto globalIdx = track.globalIndex();
    // V0 channel checks (kPos/kNeg differ by chargeIndex)
    if (doV0K0s && binarySearchAnyList(v0CndtDauList[kV0K0s][chargeIndex], globalIdx) != -1)
      BITSET(decayDauTagBit, ID_BIT_PI);

    if (doV0Lambda) {
      if (binarySearchAnyList(v0CndtDauList[kV0Lambda][chargeIndex], globalIdx) != -1) {
        BITSET(decayDauTagBit, (chargeIndex == kPos) ? ID_BIT_PR : ID_BIT_PI);
      }
    }

    if (doV0AntiLambda) {
      if (binarySearchAnyList(v0CndtDauList[kV0AntiLambda][chargeIndex], globalIdx) != -1) {
        BITSET(decayDauTagBit, (chargeIndex == kPos) ? ID_BIT_PI : ID_BIT_PR);
      }
    }

    if (doV0Gamma && binarySearchAnyList(v0CndtDauList[kV0Gamma][chargeIndex], globalIdx) != -1)
      BITSET(decayDauTagBit, ID_BIT_EL);

    // PrimVtx candidate decays
    if (doPhi1020 && binarySearchAnyList(primVtxCndtDauList[kPrimTrkPhi1020][chargeIndex], globalIdx) != -1)
      BITSET(decayDauTagBit, ID_BIT_KA);

    if (doJPsiToEE && binarySearchAnyList(primVtxCndtDauList[kPrimTrkJPsiToEE][chargeIndex], globalIdx) != -1)
      BITSET(decayDauTagBit, ID_BIT_EL);

    if (doJPsiToMuMu && binarySearchAnyList(primVtxCndtDauList[kPrimTrkJPsiToMuMu][chargeIndex], globalIdx) != -1)
      BITSET(decayDauTagBit, ID_BIT_MU);

    if (doKStar892 && binarySearchAnyList(primVtxCndtDauList[kPrimTrkKStar892][chargeIndex], globalIdx) != -1)
      BITSET(decayDauTagBit, (chargeIndex == kPos) ? ID_BIT_KA : ID_BIT_PI);

    if (doKStar892Bar && binarySearchAnyList(primVtxCndtDauList[kPrimTrkKStar892Bar][chargeIndex], globalIdx) != -1)
      BITSET(decayDauTagBit, (chargeIndex == kPos) ? ID_BIT_PI : ID_BIT_KA);

    if (doRho770 && binarySearchAnyList(primVtxCndtDauList[kPrimTrkRho770][chargeIndex], globalIdx) != -1)
      BITSET(decayDauTagBit, ID_BIT_PI);
  }

  template <int Mode, int fillMode, int pidMode, int signMode, int detMode, typename H, typename T>
  void fillAllTypeQA(H& histReg, const T& track, const double& particleMass, const float& effWeight)
  {
    fillTrackQAHistos<Mode, fillMode, pidMode, signMode, detMode>(histReg, track, particleMass, effWeight);
    fillTrackMomentumQAHistos<Mode, fillMode, pidMode, signMode, detMode>(histReg, track);
    fillTrackIdQAHistos<Mode, fillMode, pidMode, signMode, detMode>(histReg, track);
    fillTrackDcaQAHistos<Mode, fillMode, pidMode, signMode, detMode>(histReg, track);
  }

  template <int Mode, int fillMode, int pidMode, int signMode, typename H, typename T>
  void fillIdentifiedTrackFullQA(H& histReg, const T& track, const int& idMethod, const double& particleMass, const float& effWeight)
  {
    if (idMethod == kTPCidentified) {
      fillAllTypeQA<Mode, fillMode, pidMode, signMode, tpcId>(histReg, track, particleMass, effWeight);
    } else if (idMethod == kTPCTOFidentified) {
      fillAllTypeQA<Mode, fillMode, pidMode, signMode, tpctofId>(histReg, track, particleMass, effWeight);
    } else if (idMethod == kUnidentified) {
      fillAllTypeQA<Mode, fillMode, pidMode, signMode, tpctofId>(histReg, track, particleMass, effWeight);
    }
  }

  template <int pidMode, int signMode, int idBit, typename H, typename T>
  void fillIdentifiedTrackQASign(H& histReg1, H& histReg2, const T& track, const double& particleMass, const int& idMethod,
                                 int& iNRejected, float& fNRejected, int& iNTrk, float& fNTrk, float& effWeight, float& effWeightSum,
                                 const int& ptEtaBin, const int& decayDauTagBit)
  {
    // Cache config values once per template instantiation
    static const bool countPrompt = getCfg<bool>(cfgTrackCuts.countSetting, pidMode, kCountPrompt);
    static const bool countNonPrompt = getCfg<bool>(cfgTrackCuts.countSetting, pidMode, kCountNonPrompt);
    static const bool countDauContam = getCfg<bool>(cfgTrackCuts.countSetting, pidMode, kCountDauContam);

    static const bool fillPrompt = getCfg<bool>(cfgTrackCuts.countSetting, pidMode, kFillPrompt);
    static const bool fillNonPrompt = getCfg<bool>(cfgTrackCuts.countSetting, pidMode, kFillNonPrompt);
    static const bool fillDauContam = getCfg<bool>(cfgTrackCuts.countSetting, pidMode, kFillDauContam);
    static const bool fillRejected = getCfg<bool>(cfgTrackCuts.countSetting, pidMode, kFillRejected);

    bool isPrompTrk = (decayDauTagBit == 0);
    bool isNonPromptTrk = BOOL_BITCHECK(decayDauTagBit, idBit);
    bool isDauContam = BOOL_OTHERBITSON(decayDauTagBit, idBit);

    bool countIt = (countPrompt && isPrompTrk) || (countNonPrompt && isNonPromptTrk) || (countDauContam && isDauContam);
    effWeight = hPtEtaForEffCorrection[pidMode][signMode]->GetBinContent(ptEtaBin);

    if (countIt) {
      iNTrk += 1;
      fNTrk += effWeight;
      effWeightSum += effWeight;
      if (fillPrompt && isPrompTrk) {
        fillIdentifiedTrackFullQA<recoAnalysisPrompt, kFillSimple, pidMode, signMode>(histReg1, track, idMethod, particleMass, effWeight);
      }
      if (fillNonPrompt && isNonPromptTrk) {
        fillIdentifiedTrackFullQA<recoAnalysisNonPrompt, kFillSimple, pidMode, signMode>(histReg1, track, idMethod, particleMass, effWeight);
      }
      if (fillDauContam && isDauContam) {
        fillIdentifiedTrackFullQA<recoAnalysisDauContam, kFillSimple, pidMode, signMode>(histReg1, track, idMethod, particleMass, effWeight);
      }
    } else {
      iNRejected += 1;
      fNRejected += effWeight;
      if (fillRejected) {
        fillIdentifiedTrackFullQA<recoAnalysisRejected, kFillSimple, pidMode, signMode>(histReg2, track, idMethod, particleMass, effWeight);
      }
    }
  }

  template <int pidMode, int idBit, typename H, typename T>
  void executeIdentifiedTrackPart(const bool& trackIsType, H& histReg1, H& histReg2, const T& track, const double& particleMass, const int& idMethod,
                                  auto& iNRejected, auto& fNRejected, auto& iNTrk, auto& fNTrk, auto& effWeight, auto& effWeightSum,
                                  const int& ptEtaBin, const int& decayDauTagBit)
  {
    if (!trackIsType)
      return;
    if (track.signed1Pt() > 0) {
      fillIdentifiedTrackQASign<pidMode, kPos, idBit>(histReg1, histReg2, track, particleMass, idMethod,
                                                      iNRejected[kPos], fNRejected[kPos], iNTrk[kPos], fNTrk[kPos], effWeight[kPos], effWeightSum[kPos],
                                                      ptEtaBin, decayDauTagBit);
    } else {
      fillIdentifiedTrackQASign<pidMode, kNeg, idBit>(histReg1, histReg2, track, particleMass, idMethod,
                                                      iNRejected[kNeg], fNRejected[kNeg], iNTrk[kNeg], fNTrk[kNeg], effWeight[kNeg], effWeightSum[kNeg],
                                                      ptEtaBin, decayDauTagBit);
    }
  }

  template <typename T>
  void executeTrackQAPart(const T& track, int& rejectionTag, int& nTrack, bool& isAcceptedTrack)
  {
    if (cfgFill.cfgFill07RecoTrackPreSel) {
      fillTrackQA<recoTrackPreSel>(recoTracks, track);
    }
    rejectionTag = kPassed;
    if (!checkTrackSelection(track, rejectionTag)) {
      recoAnalysis.fill(HIST("recoAnalysis/RejectedTrack_RejectionTag"), rejectionTag);
      isAcceptedTrack = false;
      return; // for continue;
    }
    isAcceptedTrack = true;
    if (cfgFill.cfgFill08RecoTrackPostSel) {
      fillTrackQA<recoTrackPostSel>(recoTracks, track);
    }
    nTrack++;
  }

  template <typename T>
  void executeTrackAnalysisPart(const T& track, const int& trackIdTag, const int& ptEtaBin, const int& decayDauTagBit, const auto& idMethod, const auto& trackIs,
                                auto& iNTrkPi, auto& fNTrkPi, auto& effWeightPi, auto& effWeightSumPi, auto& iNTrkRejectedPi, auto& fNTrkRejectedPi,
                                auto& iNTrkKa, auto& fNTrkKa, auto& effWeightKa, auto& effWeightSumKa, auto& iNTrkRejectedKa, auto& fNTrkRejectedKa,
                                auto& iNTrkPr, auto& fNTrkPr, auto& effWeightPr, auto& effWeightSumPr, auto& iNTrkRejectedPr, auto& fNTrkRejectedPr,
                                auto& iNTrkEl, auto& fNTrkEl, auto& effWeightEl, auto& effWeightSumEl, auto& iNTrkRejectedEl, auto& fNTrkRejectedEl,
                                auto& iNTrkDe, auto& fNTrkDe, auto& effWeightDe, auto& effWeightSumDe, auto& iNTrkRejectedDe, auto& fNTrkRejectedDe)
  {

    executeIdentifiedTrackPart<kPi, ID_BIT_PI>(trackIs[kPi], recoAnalysisPi, recoRejectedPi, track, MassPiPlus, idMethod[kPi],
                                               iNTrkRejectedPi, fNTrkRejectedPi, iNTrkPi, fNTrkPi, effWeightPi, effWeightSumPi,
                                               ptEtaBin, decayDauTagBit);

    executeIdentifiedTrackPart<kKa, ID_BIT_KA>(trackIs[kKa], recoAnalysisKa, recoRejectedKa, track, MassKPlus, idMethod[kKa],
                                               iNTrkRejectedKa, fNTrkRejectedKa, iNTrkKa, fNTrkKa, effWeightKa, effWeightSumKa,
                                               ptEtaBin, decayDauTagBit);

    executeIdentifiedTrackPart<kPr, ID_BIT_PR>(trackIs[kPr], recoAnalysisPr, recoRejectedPr, track, MassProton, idMethod[kPr],
                                               iNTrkRejectedPr, fNTrkRejectedPr, iNTrkPr, fNTrkPr, effWeightPr, effWeightSumPr,
                                               ptEtaBin, decayDauTagBit);

    executeIdentifiedTrackPart<kEl, ID_BIT_EL>(trackIs[kEl], recoAnalysisEl, recoRejectedEl, track, MassElectron, idMethod[kEl],
                                               iNTrkRejectedEl, fNTrkRejectedEl, iNTrkEl, fNTrkEl, effWeightEl, effWeightSumEl,
                                               ptEtaBin, decayDauTagBit);

    executeIdentifiedTrackPart<kDe, ID_BIT_DE>(trackIs[kDe], recoAnalysisDe, recoRejectedDe, track, MassDeuteron, idMethod[kDe],
                                               iNTrkRejectedDe, fNTrkRejectedDe, iNTrkDe, fNTrkDe, effWeightDe, effWeightSumDe,
                                               ptEtaBin, decayDauTagBit);

    recoAnalysis.fill(HIST("recoAnalysis/SelectedTrack_IdentificationTag"), trackIdTag);
  }

  template <int Mode, typename H, typename T>
  void fillTrackQA(H& histReg, const T& track)
  {
    fillTrackQAHistos<Mode, kFillSimple, kNonePid, kNoneSign, KNoneDet>(histReg, track, 0, 1); // 0 and 1 are mass and weight
    fillTrackMomentumQAHistos<Mode, kFillSimple, kNonePid, kNoneSign, KNoneDet>(histReg, track);
    fillTrackDcaQAHistos<Mode, kFillSimple, kNonePid, kNoneSign, KNoneDet>(histReg, track);

    if (track.sign() > 0) {
      fillTrackIdQAHistos<Mode, kFillSimple, kPi, kPos, NoId>(histReg, track); // Look at Pion
      fillTrackIdQAHistos<Mode, kFillSimple, kKa, kPos, NoId>(histReg, track); // Look at Kaon
      fillTrackIdQAHistos<Mode, kFillSimple, kPr, kPos, NoId>(histReg, track); // Look at Proton
      fillTrackIdQAHistos<Mode, kFillSimple, kEl, kPos, NoId>(histReg, track); // Look at Electron
      fillTrackIdQAHistos<Mode, kFillSimple, kDe, kPos, NoId>(histReg, track); // Look at Deuteron
    } else {
      fillTrackIdQAHistos<Mode, kFillSimple, kPi, kNeg, NoId>(histReg, track); // Look at Pion
      fillTrackIdQAHistos<Mode, kFillSimple, kKa, kNeg, NoId>(histReg, track); // Look at Kaon
      fillTrackIdQAHistos<Mode, kFillSimple, kPr, kNeg, NoId>(histReg, track); // Look at Proton
      fillTrackIdQAHistos<Mode, kFillSimple, kEl, kNeg, NoId>(histReg, track); // Look at Electron
      fillTrackIdQAHistos<Mode, kFillSimple, kDe, kNeg, NoId>(histReg, track); // Look at Deuteron
    }
  }

  template <int pidMode, int signMode>
  void fillTrackCountQA(const int& iNTrk, const float& fNTrk, const float& effWeightSum)
  {
    recoEvent.fill(HIST("recoEvent/intCount/") + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST("hCount"), iNTrk);
    recoEvent.fill(HIST("recoEvent/floatCount/") + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST("hCount"), fNTrk);
    recoEvent.fill(HIST("recoEvent/effWeightSum/") + HIST(PidDire[pidMode]) + HIST(SignDire[signMode]) + HIST("hCount"), effWeightSum);
  }

  template <int dirMode, typename T>
  void fillCollQA(const T& col)
  {
    if (col.sel8())
      recoEvent.fill(HIST(HistRegDire[dirMode]) + HIST("CollType"), 0);
    if (col.selection_bit(o2::aod::evsel::kNoSameBunchPileup))
      recoEvent.fill(HIST(HistRegDire[dirMode]) + HIST("CollType"), 1);
    if (col.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV))
      recoEvent.fill(HIST(HistRegDire[dirMode]) + HIST("CollType"), 2);
    if (col.selection_bit(o2::aod::evsel::kIsVertexITSTPC))
      recoEvent.fill(HIST(HistRegDire[dirMode]) + HIST("CollType"), 3);
  }

  template <typename C, typename T>
  void executeEventInfoPart(const C& collision, const float& centrality, const int v0TableSize, const T& tracksTableSize,
                            const int& nTrack, const auto& iNTrk, const auto& fNTrk, const auto& effWeightSum)
  {
    // Collisions QA
    recoEvent.fill(HIST("recoEvent/h01_CollisionCount"), 0.5);
    recoEvent.fill(HIST("recoEvent/h02_VertexXRec"), collision.posX());
    recoEvent.fill(HIST("recoEvent/h03_VertexYRec"), collision.posY());
    recoEvent.fill(HIST("recoEvent/h04_VertexZRec"), collision.posZ());
    recoEvent.fill(HIST("recoEvent/h05_Centrality"), centrality);
    recoEvent.fill(HIST("recoEvent/h06_V0Size"), v0TableSize);
    recoEvent.fill(HIST("recoEvent/h07_TracksSize"), tracksTableSize);
    recoEvent.fill(HIST("recoEvent/h08_nTrack"), nTrack);

    fillTrackCountQA<kPi, kPos>(iNTrk[kPi][kPos], fNTrk[kPi][kPos], effWeightSum[kPi][kPos]);
    fillTrackCountQA<kKa, kPos>(iNTrk[kKa][kPos], fNTrk[kKa][kPos], effWeightSum[kKa][kPos]);
    fillTrackCountQA<kPr, kPos>(iNTrk[kPr][kPos], fNTrk[kPr][kPos], effWeightSum[kPr][kPos]);
    fillTrackCountQA<kEl, kPos>(iNTrk[kEl][kPos], fNTrk[kEl][kPos], effWeightSum[kEl][kPos]);
    fillTrackCountQA<kMu, kPos>(iNTrk[kMu][kPos], fNTrk[kMu][kPos], effWeightSum[kMu][kPos]);
    fillTrackCountQA<kDe, kPos>(iNTrk[kDe][kPos], fNTrk[kDe][kPos], effWeightSum[kDe][kPos]);

    fillTrackCountQA<kPi, kNeg>(iNTrk[kPi][kNeg], fNTrk[kPi][kNeg], effWeightSum[kPi][kNeg]);
    fillTrackCountQA<kKa, kNeg>(iNTrk[kKa][kNeg], fNTrk[kKa][kNeg], effWeightSum[kKa][kNeg]);
    fillTrackCountQA<kPr, kNeg>(iNTrk[kPr][kNeg], fNTrk[kPr][kNeg], effWeightSum[kPr][kNeg]);
    fillTrackCountQA<kEl, kNeg>(iNTrk[kEl][kNeg], fNTrk[kEl][kNeg], effWeightSum[kEl][kNeg]);
    fillTrackCountQA<kMu, kNeg>(iNTrk[kMu][kNeg], fNTrk[kMu][kNeg], effWeightSum[kMu][kNeg]);
    fillTrackCountQA<kDe, kNeg>(iNTrk[kDe][kNeg], fNTrk[kDe][kNeg], effWeightSum[kDe][kNeg]);

    fillTrackCountQA<kV0K0s, kNoneSign>(iNTrk[kV0K0s][kPos], fNTrk[kV0K0s][kPos], effWeightSum[kV0K0s][kPos]);
    fillTrackCountQA<kV0Lambda, kNoneSign>(iNTrk[kV0Lambda][kPos], fNTrk[kV0Lambda][kPos], effWeightSum[kV0Lambda][kPos]);
    fillTrackCountQA<kV0AntiLambda, kNoneSign>(iNTrk[kV0AntiLambda][kPos], fNTrk[kV0AntiLambda][kPos], effWeightSum[kV0AntiLambda][kPos]);
    fillTrackCountQA<kV0Gamma, kNoneSign>(iNTrk[kV0Gamma][kPos], fNTrk[kV0Gamma][kPos], effWeightSum[kV0Gamma][kPos]);

    fillTrackCountQA<kPrimPhi1020, kNoneSign>(iNTrk[kPrimPhi1020][kPos], fNTrk[kPrimPhi1020][kPos], effWeightSum[kPrimPhi1020][kPos]);
    fillTrackCountQA<kPrimJPsiToEE, kNoneSign>(iNTrk[kPrimJPsiToEE][kPos], fNTrk[kPrimJPsiToEE][kPos], effWeightSum[kPrimJPsiToEE][kPos]);
    fillTrackCountQA<kPrimJPsiToMuMu, kNoneSign>(iNTrk[kPrimJPsiToMuMu][kPos], fNTrk[kPrimJPsiToMuMu][kPos], effWeightSum[kPrimJPsiToMuMu][kPos]);
    fillTrackCountQA<kPrimKStar892, kNoneSign>(iNTrk[kPrimKStar892][kPos], fNTrk[kPrimKStar892][kPos], effWeightSum[kPrimKStar892][kPos]);
    fillTrackCountQA<kPrimKStar892Bar, kNoneSign>(iNTrk[kPrimKStar892Bar][kPos], fNTrk[kPrimKStar892Bar][kPos], effWeightSum[kPrimKStar892Bar][kPos]);
    fillTrackCountQA<kPrimRho770, kNoneSign>(iNTrk[kPrimRho770][kPos], fNTrk[kPrimRho770][kPos], effWeightSum[kPrimRho770][kPos]);

    fillTrackCountQA<kRejectedPi, kPos>(iNTrk[kRejectedPi][kPos], fNTrk[kRejectedPi][kPos], -1);
    fillTrackCountQA<kRejectedKa, kPos>(iNTrk[kRejectedKa][kPos], fNTrk[kRejectedKa][kPos], -1);
    fillTrackCountQA<kRejectedPr, kPos>(iNTrk[kRejectedPr][kPos], fNTrk[kRejectedPr][kPos], -1);
    fillTrackCountQA<kRejectedEl, kPos>(iNTrk[kRejectedEl][kPos], fNTrk[kRejectedEl][kPos], -1);
    fillTrackCountQA<kRejectedMu, kPos>(iNTrk[kRejectedMu][kPos], fNTrk[kRejectedMu][kPos], -1);
    fillTrackCountQA<kRejectedDe, kPos>(iNTrk[kRejectedDe][kPos], fNTrk[kRejectedDe][kPos], -1);

    fillTrackCountQA<kRejectedPi, kNeg>(iNTrk[kRejectedPi][kNeg], fNTrk[kRejectedPi][kNeg], -1);
    fillTrackCountQA<kRejectedKa, kNeg>(iNTrk[kRejectedKa][kNeg], fNTrk[kRejectedKa][kNeg], -1);
    fillTrackCountQA<kRejectedPr, kNeg>(iNTrk[kRejectedPr][kNeg], fNTrk[kRejectedPr][kNeg], -1);
    fillTrackCountQA<kRejectedEl, kNeg>(iNTrk[kRejectedEl][kNeg], fNTrk[kRejectedEl][kNeg], -1);
    fillTrackCountQA<kRejectedMu, kPos>(iNTrk[kRejectedMu][kPos], fNTrk[kRejectedMu][kPos], -1);
    fillTrackCountQA<kRejectedDe, kNeg>(iNTrk[kRejectedDe][kNeg], fNTrk[kRejectedDe][kNeg], -1);
  }

  template <typename T, typename U, typename V, typename W>
  inline void fillAllSparses(
    const T& hists,
    const U& varPtrs,
    const V& casters,
    W& buffers)
  {
    // Loop over the sparse histograms
    for (uint iHist = 0; iHist < hists.size(); iHist++) {
      const int nDim = static_cast<int>(varPtrs[iHist].size());
      for (int i = 0; i < nDim; ++i) {
        buffers[iHist][i] = casters[iHist][i](varPtrs[iHist][i]);
      }
      hists[iHist]->Fill(buffers[iHist].data());
    }
  }

  //________________MC__Analysis related funtions_________________________
  template <int analysisType>
  inline bool isTrueMcMatch(uint8_t mcTag, int idBit)
  {
    if constexpr (analysisType == doPurityProcessing) {
      return BOOL_BITCHECK(mcTag, idBit);
    } else {
      return true; // Always allow if not doing purity processing
    }
  }

  template <int Mode, int pidMode, typename H, typename T>
  void fillGenTrackQA(H& histReg, const T& mcTrack)
  {
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST("h12_p"), mcTrack.p());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST("h13_pt"), mcTrack.pt());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST("h14_eta"), mcTrack.eta());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST("h15_phi"), mcTrack.phi());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST("h16_rapidity"), mcTrack.y());
    histReg.fill(HIST(HistRegDire[Mode]) + HIST(PidDire[pidMode]) + HIST("h20_pt_eta"), mcTrack.pt(), mcTrack.eta());
  }

  void getV0McBitTag(auto& v0DecayTrueMcTag, const auto& v0mcparticle, const auto& posDauMcPart, const auto& negDauMcPart)
  {
    v0DecayTrueMcTag = 0;
    const int v0Code = v0mcparticle.pdgCode();
    const int posCode = posDauMcPart.pdgCode();
    const int negCode = negDauMcPart.pdgCode();

    if (v0Code == kK0Short && posCode == kPiPlus && negCode == kPiMinus) {
      BITSET(v0DecayTrueMcTag, kV0TrkK0s);
    } else if (v0Code == kLambda0 && posCode == kProton && negCode == kPiMinus) {
      BITSET(v0DecayTrueMcTag, kV0TrkLambda);
    } else if (v0Code == kLambda0Bar && posCode == kPiPlus && negCode == kProtonBar) {
      BITSET(v0DecayTrueMcTag, kV0TrkAntiLambda);
    } else if (v0Code == kGamma && posCode == kPositron && negCode == kElectron) {
      BITSET(v0DecayTrueMcTag, kV0TrkGamma);
    }
  }

  template <bool Mode>
  void getPvcMcBitTag(auto& primVtxCndtMcTag, const auto& pvcMcParticle, const auto& posDauMcPart, const auto& negDauMcPart)
  {
    primVtxCndtMcTag = 0;

    const int posCode = posDauMcPart.pdgCode();
    const int negCode = negDauMcPart.pdgCode();

    if constexpr (Mode) {
      const int pvcCode = pvcMcParticle.pdgCode();
      if (pvcCode == kPhi && posCode == kKPlus && negCode == kKMinus) {
        BITSET(primVtxCndtMcTag, kPrimPhi1020);
        return;
      }
      if (pvcCode == kJPsi && posCode == kPositron && negCode == kElectron) {
        BITSET(primVtxCndtMcTag, kPrimJPsiToEE);
        return;
      }
      if (pvcCode == kJPsi && posCode == kMuonPlus && negCode == kMuonMinus) {
        BITSET(primVtxCndtMcTag, kPrimJPsiToMuMu);
        return;
      }
      if (pvcCode == kK0Star892 && posCode == kKPlus && negCode == kPiMinus) {
        BITSET(primVtxCndtMcTag, kPrimKStar892);
        return;
      }
      if (pvcCode == -kK0Star892 && posCode == kPiPlus && negCode == kKMinus) {
        BITSET(primVtxCndtMcTag, kPrimKStar892Bar);
        return;
      }
      if (pvcCode == kRho770_0 && posCode == kPiPlus && negCode == kPiMinus) {
        BITSET(primVtxCndtMcTag, kPrimRho770);
        return;
      }
    } else {
      if (posCode == kKPlus && negCode == kKMinus) {
        BITSET(primVtxCndtMcTag, kPrimPhi1020);
        return;
      }
      if (posCode == kPositron && negCode == kElectron) {
        BITSET(primVtxCndtMcTag, kPrimJPsiToEE);
        return;
      }
      if (posCode == kMuonPlus && negCode == kMuonMinus) {
        BITSET(primVtxCndtMcTag, kPrimJPsiToMuMu);
        return;
      }
      if (posCode == kKPlus && negCode == kPiMinus) {
        BITSET(primVtxCndtMcTag, kPrimKStar892);
        return;
      }
      if (posCode == kPiPlus && negCode == kKMinus) {
        BITSET(primVtxCndtMcTag, kPrimKStar892Bar);
        return;
      }
      if (posCode == kPiPlus && negCode == kPiMinus) {
        BITSET(primVtxCndtMcTag, kPrimRho770);
        return;
      }
    }
  }

  void getTrackMcBitTag(auto& trackTrueMcTag, const auto& mcPart)
  {
    trackTrueMcTag = 0;
    const int absPdg = std::abs(mcPart.pdgCode());

    // Use if-else ladder for faster branch prediction
    if (absPdg == kPiPlus) {
      BITSET(trackTrueMcTag, ID_BIT_PI);
    } else if (absPdg == kKPlus) {
      BITSET(trackTrueMcTag, ID_BIT_KA);
    } else if (absPdg == kProton) {
      BITSET(trackTrueMcTag, ID_BIT_PR);
    } else if (absPdg == kElectron) {
      BITSET(trackTrueMcTag, ID_BIT_EL);
    } else if (absPdg == kMuonMinus) {
      BITSET(trackTrueMcTag, ID_BIT_MU);
    } else if (absPdg == kDeuteron) {
      BITSET(trackTrueMcTag, ID_BIT_DE);
    }
  }

  template <typename T>
  void getV0MCount(const T& mcTrack, auto& multV0M)
  {
    if ((cfgSim.cfgSim09V0CLow < mcTrack.eta() && mcTrack.eta() < cfgSim.cfgSim09V0CUp) || (cfgSim.cfgSim09V0ALow < mcTrack.eta() && mcTrack.eta() < cfgSim.cfgSim09V0AUp)) {
      if (cfgSim.cfgSim05doFWDPtDependentCheck) {
        if (mcTrack.pt() > cfgSim.cfgSim06FWDPtCut) {
          multV0M++; // V0C: at -3.7 < η < -1.7 (backward direction).
                     // V0A: at 2.8 < η < 5.1 (forward direction).
        }
      } else {
        multV0M++;
      }
    }
  }

  // Event Filter
  Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutZvertex);

  // Track Filter
  Filter ptFilter = (o2::aod::track::pt) > cfgTrackCuts.cfgTrk05PtLow && (o2::aod::track::pt) < cfgTrackCuts.cfgTrk06PtHigh;
  Filter etaFilter = (nabs(o2::aod::track::eta) < cfgTrackCuts.cfgTrk04Eta);

  // Filters on V0s
  Filter preFilterv0 = (nabs(aod::v0data::dcapostopv) > v0settingDcaPosToPV &&
                        nabs(aod::v0data::dcanegtopv) > v0settingDcaNegToPV &&
                        aod::v0data::dcaV0daughters < v0settingDcaV0Dau);

  using MyCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels,
                                               aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::Mults>>;

  using MyTracks = soa::Filtered<soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                           aod::TOFSignal, aod::pidTOFbeta, aod::pidTOFmass,
                                           aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullDe,
                                           aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullDe>>;

  using MyV0s = soa::Filtered<aod::V0Datas>;
  using MyPrimVtxCndts = aod::PrimVtxCndts;

  // For manual sliceBy
  SliceCache cache;

  Preslice<MyTracks> tracksPerCollisionPreslice = o2::aod::track::collisionId;
  Preslice<MyV0s> v0sPerCollisionPreslice = o2::aod::track::collisionId;
  Preslice<MyPrimVtxCndts> primVtxCndtsPreslice = o2::aod::primvtxcndt::collisionId;

  Partition<MyTracks> posTracks = aod::track::signed1Pt > 0.0f;
  Partition<MyTracks> negTracks = aod::track::signed1Pt < 0.0f;

  using MyCollisionsWithMcLabels = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels,
                                                           aod::CentFT0Ms, aod::CentFT0Cs, aod::Mults, aod::McCollisionLabels>>;

  using MyTracksWithMcLabels = soa::Filtered<soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                                                       aod::TOFSignal, aod::pidTOFbeta, // aod::pidTOFmass,
                                                       aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullDe,
                                                       aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullDe, aod::McTrackLabels>>;

  using MyV0sWithMcLabels = soa::Filtered<soa::Join<aod::V0Datas, aod::McV0Labels>>;

  Partition<MyTracksWithMcLabels> posTracksWithMcLabels = aod::track::signed1Pt > 0.0f;
  Partition<MyTracksWithMcLabels> negTracksWithMcLabels = aod::track::signed1Pt < 0.0f;

  Preslice<aod::McParticles> mcTracksPerMcCollisionPreslice = o2::aod::mcparticle::mcCollisionId;
  using MyMcCollisions = aod::McCollisions;

  // Declaring vectors outside the process to avoid slight overhead for stack allocation and deallocation during each iteration.
  struct MyTrackData {
    int64_t globalIndex;
    float massSigma;
    explicit MyTrackData(int64_t gi, float mS) : globalIndex(gi), massSigma(mS) {}
  };

  std::vector<int64_t> k0sTagIndexList;
  std::array<std::array<std::vector<MyTrackData>, 2>, v0TrkEnumSize> v0CndtDauList;
  std::array<std::array<std::vector<MyTrackData>, 2>, primVtxTrkEnumSize> primVtxCndtDauList;

  int64_t dfCount = 0;
  std::chrono::high_resolution_clock::time_point start0 = std::chrono::high_resolution_clock::now();

  template <int analysisType, typename C, typename V, typename T, typename P>
  void executeAnalysis(const C& collisions, const V& V0s, const T& tracks, const P& primVtxCndts)
  {
    // Time related variables
    dfCount++;
    auto start1 = std::chrono::high_resolution_clock::now();
    auto v0Start = std::chrono::high_resolution_clock::now();
    auto v0SortStart = std::chrono::high_resolution_clock::now();
    auto collLoopStart = std::chrono::high_resolution_clock::now();
    auto manualGroupingStart = std::chrono::high_resolution_clock::now();
    auto primVtxLoop1Start = std::chrono::high_resolution_clock::now();
    auto primVtxSortStart = std::chrono::high_resolution_clock::now();
    auto primVtxLoop2Start = std::chrono::high_resolution_clock::now();
    auto secondV0LoopStart = std::chrono::high_resolution_clock::now();
    auto trackLoopStart = std::chrono::high_resolution_clock::now();
    auto sparseFillStart = std::chrono::high_resolution_clock::now();
    auto fillEventInfoStart = std::chrono::high_resolution_clock::now();

    const int pidVecSize = kNKaPosNKaNegProd + 1; // kPrimRho770 is last enum value of different particles. and kNKaPosNKaNegProd is of last count types
    // Variables for  collision loop and track counting
    float centrality = 0;
    int nTrack;
    std::array<std::array<int, 2>, pidVecSize> iNTrk;
    std::array<std::array<float, 2>, pidVecSize> fNTrk;
    std::array<std::array<float, 2>, pidVecSize> effWeightSum;

    // Variables for sparse histogram filling
    // Get Size of all Three sparse histograms
    std::array<int, 3> hSparseNAxis;
    hSparseNAxis[0] = getCfg<int>(cfgAxis.sparseSetting, 0, 0);
    hSparseNAxis[1] = getCfg<int>(cfgAxis.sparseSetting, 1, 0);
    hSparseNAxis[2] = getCfg<int>(cfgAxis.sparseSetting, 2, 0);

    // Creating pointers to directly fill desired value in histogram to avoid any looking for variables.
    // Variable pointers for values to be filled in sparse histos, values can be int or float so need variable pointer
    std::array<std::vector<void*>, 3> hSparseVarPtrs{{std::vector<void*>(hSparseNAxis[0], nullptr),
                                                      std::vector<void*>(hSparseNAxis[1], nullptr),
                                                      std::vector<void*>(hSparseNAxis[2], nullptr)}};

    // typeid of pointers to recast them while filling.
    std::array<std::vector<std::type_index>, 3> hSparseVarTypes{{std::vector<std::type_index>(hSparseNAxis[0], typeid(void)),
                                                                 std::vector<std::type_index>(hSparseNAxis[1], typeid(void)),
                                                                 std::vector<std::type_index>(hSparseNAxis[2], typeid(void))}};

    // Obtaining pointers to the variables and their typeid
    int iPidMode, signType, variableType;
    int nSparse = 3;
    for (int iSparse = 0; iSparse < nSparse; iSparse++) {           // Looping over sparse histos
      for (int iAxis = 0; iAxis < hSparseNAxis[iSparse]; iAxis++) { // Looping over the available axis.
        // For each axis
        iPidMode = sparseHistFillerIndices[iSparse][iAxis][0];
        signType = sparseHistFillerIndices[iSparse][iAxis][1];
        variableType = sparseHistFillerIndices[iSparse][iAxis][2];

        if (signType < 0 || signType > 1) {
          LOG(fatal) << "DEBUG :: signType is wrong, stopping the process :: signType == " << signType;
        }
        // Get the pointers
        if (iPidMode == kCentrality) {
          hSparseVarPtrs[iSparse][iAxis] = &centrality;
          hSparseVarTypes[iSparse][iAxis] = typeid(centrality);
        } else if (iPidMode == kNTrack) {
          hSparseVarPtrs[iSparse][iAxis] = &nTrack;
          hSparseVarTypes[iSparse][iAxis] = typeid(nTrack);
        } else {
          if constexpr (analysisType == doDataProcessing || analysisType == doRecoProcessing || analysisType == doPurityProcessing) {
            if (variableType == kIntCount) {
              hSparseVarPtrs[iSparse][iAxis] = &iNTrk[iPidMode][signType];
              hSparseVarTypes[iSparse][iAxis] = typeid(iNTrk[iPidMode][signType]);
            } else if (variableType == kFloatCount) {
              hSparseVarPtrs[iSparse][iAxis] = &fNTrk[iPidMode][signType];
              hSparseVarTypes[iSparse][iAxis] = typeid(fNTrk[iPidMode][signType]);
            } else if (variableType == kEffWeightSum) {
              hSparseVarPtrs[iSparse][iAxis] = &effWeightSum[iPidMode][signType];
              hSparseVarTypes[iSparse][iAxis] = typeid(effWeightSum[iPidMode][signType]);
            } else {
              LOG(fatal) << "DEBUG :: Varaible type mismatch is there :: iSparse = " << iSparse << " :: iAxis = " << iAxis << " :: variableType = " << variableType;
            }
          } else if constexpr (analysisType == doGenProcessing || analysisType == doSimProcessing) {
            if (variableType == kIntCount || variableType == kFloatCount || variableType == kEffWeightSum) {
              hSparseVarPtrs[iSparse][iAxis] = &iNTrk[iPidMode][signType];
              hSparseVarTypes[iSparse][iAxis] = typeid(iNTrk[iPidMode][signType]);
            } else {
              LOG(fatal) << "DEBUG :: Varaibel type mismatch is there :: iSparse = " << iSparse << " :: iAxis = " << iAxis << " :: variableType = " << variableType;
            }
          }
        }
      } // End of Axis Loop.
    } // End of Sparse histogram loop

    // type caster lambdas for fast addressing of pointers and sparse filling without any overhead
    std::array<std::vector<GenericCaster<double>>, 3> doubleCasters{{setupCasters<double>(hSparseVarTypes[0]),
                                                                     setupCasters<double>(hSparseVarTypes[1]),
                                                                     setupCasters<double>(hSparseVarTypes[2])}};

    // create buffers for storing the values and filling them
    // std::shared_ptr<THnSparse>
    std::array<std::shared_ptr<THnSparse>, 3> hSparses;
    std::array<std::array<double, 8>, 3> buffers;
    // preallocated stack buffers for obtaining variables // To Use hSparse->Fill(buffer.data());

    hSparses[0] = recoAnalysis.get<THnSparse>(HIST("recoAnalysis/Sparse0"));
    hSparses[1] = recoAnalysis.get<THnSparse>(HIST("recoAnalysis/Sparse1"));
    hSparses[2] = recoAnalysis.get<THnSparse>(HIST("recoAnalysis/Sparse2"));
    // End of sparse histogram variables.

    if constexpr (analysisType == doDataProcessing || analysisType == doRecoProcessing || analysisType == doPurityProcessing) {

      // Declaring variables outside the loop to avoid slight overhead for stack allocation and deallocation during each iteration.
      // Variables and Efficiency Weigths for V0s, Prim Vtx particles and  Tracks
      std::array<std::array<int, 2>, 6> idMethodSignTrk; // 2 for pos and neg, 6 for Pi, Ka, Pr, El, Mu, De
      std::array<std::array<float, 2>, pidVecSize> effWeight;

      // Varaibles for V0 Loop
      int ptEtaBinV0 = -1, ptEtaBinPosDau = -1, ptEtaBinNegDau = -1;
      int v0Tag = 0, trueV0TagValue = 0;
      int v0DauCollisionIndexTag = 0, v0DauBCTag = 0;

      // Variables for track loop
      int rejectionTag = 0;
      bool isAcceptedTrack = true;
      int nTrack = 0;
      int trackIdTag = 0;
      int ptEtaBinTrk = -1;

      std::array<bool, 6> trackIs;
      std::array<int, 6> idMethodTrk;

      int decayDauTagBit = 0;

      const bool doV0K0s = getCfg<bool>(cfgPrimVtxParticleCuts, kV0K0s, kV0DoRecoCheck);
      const bool doV0Lambda = getCfg<bool>(cfgPrimVtxParticleCuts, kV0Lambda, kV0DoRecoCheck);
      const bool doV0AntiLambda = getCfg<bool>(cfgPrimVtxParticleCuts, kV0AntiLambda, kV0DoRecoCheck);
      const bool doV0Gamma = getCfg<bool>(cfgPrimVtxParticleCuts, kV0Gamma, kV0DoRecoCheck);

      // Background Estimation From Primary Vertex
      const bool doPhi1020 = getCfg<bool>(cfgPrimVtxParticleCuts, kPrimTrkPhi1020, kPrimDoRecoCheck);
      const bool doJPsiToEE = getCfg<bool>(cfgPrimVtxParticleCuts, kPrimTrkJPsiToEE, kPrimDoRecoCheck);
      const bool doJPsiToMuMu = getCfg<bool>(cfgPrimVtxParticleCuts, kPrimTrkJPsiToMuMu, kPrimDoRecoCheck);
      const bool doKStar892 = getCfg<bool>(cfgPrimVtxParticleCuts, kPrimTrkKStar892, kPrimDoRecoCheck);
      const bool doKStar892Bar = getCfg<bool>(cfgPrimVtxParticleCuts, kPrimTrkKStar892Bar, kPrimDoRecoCheck);
      const bool doRho770 = getCfg<bool>(cfgPrimVtxParticleCuts, kPrimTrkRho770, kPrimDoRecoCheck);

      bool checkTrackId[6][2];
      checkTrackId[kPi][kPos] = doRho770 || doKStar892Bar;
      checkTrackId[kPi][kNeg] = doRho770 || doKStar892;
      checkTrackId[kKa][kPos] = doPhi1020 || doKStar892;
      checkTrackId[kKa][kNeg] = doPhi1020 || doKStar892Bar;
      checkTrackId[kEl][kPos] = doJPsiToEE;
      checkTrackId[kEl][kNeg] = doJPsiToEE;
      checkTrackId[kMu][kPos] = doJPsiToMuMu;
      checkTrackId[kMu][kNeg] = doJPsiToMuMu;

      bool posDauIs[6];
      bool negDauIs[6];

      uint8_t posTrackPrimVtxMotherFlag = 0;
      uint8_t negTrackPrimVtxMotherFlag = 0;

      uint8_t posTrackV0MotherFlag = 0;
      uint8_t negTrackV0MotherFlag = 0;

      int iTotalPrimVtxCndtTrk = 0;
      float fTotalPrimVtxCndtTrk = 0;

      int count1;
      int count2;
      int count4;
      int count5;
      int count6;

      std::array<float, 2> massSigma = {static_cast<float>(doublePowersOf10[8]), static_cast<float>(doublePowersOf10[8])};
      std::array<uint8_t, 2> mpBit = {0, 0};
      int requiredBit = -1;

      uint8_t v0DecayTrueMcTag = 0xFF; // will turn all bits on
      uint8_t primVtxCndtMcTag = 0xFF;
      uint8_t trackTrueMcTag = 0xFF;

      // Processing starts from here.
      k0sTagIndexList.clear();
      for (int i = kV0TrkK0s; i <= kV0TrkGamma; i++) {
        v0CndtDauList[i][kPos].clear();
        v0CndtDauList[i][kNeg].clear();
      }

      v0Start = std::chrono::high_resolution_clock::now();
      for (const auto& v0 : V0s) {
        const auto& posDaughterTrack = v0.template posTrack_as<T>();
        const auto& negDaughterTrack = v0.template negTrack_as<T>();

        v0DecayTrueMcTag = 0xFF; // By Default all bits are on, if Mc process is running the Tag will be recalculated.
        //__________________________Reco Level _______________________________________________________________
        if constexpr (analysisType == doRecoProcessing || analysisType == doPurityProcessing) {

          if (!v0.has_mcParticle() || !posDaughterTrack.has_mcParticle() || !negDaughterTrack.has_mcParticle()) {
            continue; // Reco Level Check + Truth Level Check
          }

          auto v0mcparticle = v0.mcParticle(); // Reco Level Check + Truth Level Check
          if (!v0mcparticle.isPhysicalPrimary())
            continue;

          //___________________________Truth Level_________________________________________________________
          if constexpr (analysisType == doPurityProcessing) {
            auto posDauMcPart = posDaughterTrack.mcParticle();
            auto negDauMcPart = negDaughterTrack.mcParticle();
            getV0McBitTag(v0DecayTrueMcTag, v0mcparticle, posDauMcPart, negDauMcPart);
            if (countBits(v0DecayTrueMcTag) > 1) {
              LOG(fatal) << "DEBUG :: countBits(v0DecayTrueMcTag) > 1  :: v0DecayTrueMcTag = " << std::bitset<8>(v0DecayTrueMcTag);
            }
          } //____Truth_Level_Over______
        } //____Reco_Level_Over______

        ptEtaBinV0 = hPtEtaForEffCorrection[kV0K0s][kPos]->FindBin(v0.pt(), v0.eta());
        ptEtaBinPosDau = hPtEtaForEffCorrection[kPi][kPos]->FindBin(posDaughterTrack.pt(), posDaughterTrack.eta());
        ptEtaBinNegDau = hPtEtaForEffCorrection[kPi][kNeg]->FindBin(negDaughterTrack.pt(), negDaughterTrack.eta());

        executeV0loop<analysisType, C>(posDaughterTrack, negDaughterTrack, v0,
                                       idMethodSignTrk[kPi], idMethodSignTrk[kPr], idMethodSignTrk[kEl],
                                       v0Tag, trueV0TagValue, v0DauCollisionIndexTag, v0DauBCTag, v0CndtDauList,
                                       ptEtaBinV0, ptEtaBinPosDau, ptEtaBinNegDau,
                                       effWeight[kV0K0s][kPos], effWeight[kV0Lambda][kPos], effWeight[kV0AntiLambda][kPos], effWeight[kV0Gamma][kPos],
                                       v0DecayTrueMcTag);
      } // End of V0s Loop

      anlysisTimingInfo.fill(HIST("time01_executeV0loop"), calculateTime(v0Start));

      v0SortStart = std::chrono::high_resolution_clock::now();
      executeSortPairDaughters<false>(v0CndtDauList[kV0TrkK0s][kPos], v0CndtDauList[kV0TrkK0s][kNeg], recoV0sPostK0sSelectionCut.get<TH1>(HIST(HistRegDire[v0TablePostK0sSelectionCut]) + HIST("nCommonDauOfDifferentK0s")));
      executeSortPairDaughters<false>(v0CndtDauList[kV0TrkLambda][kPos], v0CndtDauList[kV0TrkLambda][kNeg], recoV0sPostLambdaCheck.get<TH1>(HIST(HistRegDire[v0TablePostLambdaCheck]) + HIST("nCommonDauOfDifferentLambdas")));
      executeSortPairDaughters<false>(v0CndtDauList[kV0TrkAntiLambda][kPos], v0CndtDauList[kV0TrkAntiLambda][kNeg], recoV0sPostAntiLambdaCheck.get<TH1>(HIST(HistRegDire[v0TablePostAntiLambdaCheck]) + HIST("nCommonDauOfDifferentAntiLambdas")));
      executeSortPairDaughters<false>(v0CndtDauList[kV0TrkGamma][kPos], v0CndtDauList[kV0TrkGamma][kNeg], recoV0sPostGammaCheck.get<TH1>(HIST(HistRegDire[v0TablePostGammaCheck]) + HIST("nCommonDauOfDifferentGammas")));
      anlysisTimingInfo.fill(HIST("time02_v0SortTime"), calculateTime(v0SortStart));

      collLoopStart = std::chrono::high_resolution_clock::now();
      for (const auto& collision : collisions) {

        //__________________Reco Level__________________________________________________
        if constexpr (analysisType == doRecoProcessing || analysisType == doPurityProcessing) {
          if (!collision.has_mcCollision()) {
            LOG(warning) << "No MC collision for this collision, skip...";
            continue;
          }
        }

        fillCollQA<eventPreSel>(collision);

        centrality = collision.centFT0C();
        if (cfgCentAxis.centAxis04Type == kCentFT0M) {
          centrality = collision.centFT0M();
        } else if (cfgCentAxis.centAxis04Type == kMultFT0M) {
          centrality = collision.multFT0M();
        } else if (cfgCentAxis.centAxis04Type == kMultFT0C) {
          centrality = collision.multFT0C();
        }

        manualGroupingStart = std::chrono::high_resolution_clock::now();
        // group tracks, v0s manually
        const uint64_t collIdx = collision.globalIndex();
        const auto tracksTablePerColl = tracks.sliceBy(tracksPerCollisionPreslice, collIdx);
        const auto v0sTablePerColl = V0s.sliceBy(v0sPerCollisionPreslice, collIdx);
        const auto primVtxCndtsPerColl = primVtxCndts.sliceBy(primVtxCndtsPreslice, collIdx);
        anlysisTimingInfo.fill(HIST("time04_manualGrouping"), calculateTime(manualGroupingStart));

        primVtxLoop1Start = std::chrono::high_resolution_clock::now();

        for (int i = kPrimPhi1020; i <= kPrimRho770; i++) {
          iNTrk[i][kPos] = 0;
          fNTrk[i][kPos] = 0;
          effWeightSum[i][kPos] = 0;
          iNTrk[i][kNeg] = 0;
          fNTrk[i][kNeg] = 0;
          effWeightSum[i][kNeg] = 0;
        }

        for (uint i = 0; i < primVtxTrkEnumSize; i++) {
          primVtxCndtDauList[i][kPos].clear();
          primVtxCndtDauList[i][kNeg].clear();
        }

        // First Prim Vtx Loop
        for (const auto& primVtxCndt : primVtxCndtsPerColl) {
          const auto& posTrack = primVtxCndt.template posTrack_as<T>();
          const auto& negTrack = primVtxCndt.template negTrack_as<T>();
          if (!checkTrackSelection(posTrack, rejectionTag)) {
            continue;
          }
          if (!checkTrackSelection(negTrack, rejectionTag)) {
            continue;
          }

          primVtxCndtMcTag = 0xFF; // By Default all bits are on, if Mc process is running the Tag will be recalculated.
          //__________________________Reco Level _______________________________________________________________
          if constexpr (analysisType == doRecoProcessing || analysisType == doPurityProcessing) {

            if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle()) {
              continue; // Reco Level Check // Fix in Future :: !primVtxCndt.has_mcParticle() is not implement for now.
            }

            auto pvcMcParticle = nullptr; // To correct this in future
            //-->Fix This ::  auto pvcMcParticle = v0.mcParticle(); //Reco Level Check
            //-->Fix This ::  if (!pvcMcParticle.isPhysicalPrimary())
            //-->Fix This ::    continue;
            //

            //___________________________Truth Level_________________________________________________________
            if constexpr (analysisType == doPurityProcessing) {
              auto posDauMcPart = posTrack.mcParticle();
              auto negDauMcPart = negTrack.mcParticle();
              getPvcMcBitTag<false>(primVtxCndtMcTag, pvcMcParticle, posDauMcPart, negDauMcPart);
              if (countBits(primVtxCndtMcTag) > 1) {
                LOG(fatal) << "DEBUG :: countBits(primVtxCndtMcTag) > 1  :: primVtxCndtMcTag = " << std::bitset<8>(primVtxCndtMcTag);
              }
            } //____Truth_Level_Over______
          } //____Reco_Level_Over______

          int ptEtaBinPosTrk = hPtEtaForBinSearch->FindBin(posTrack.pt(), posTrack.eta()); // Find Track Bin for efficiency correction
          int ptEtaBinNegTrk = hPtEtaForBinSearch->FindBin(negTrack.pt(), negTrack.eta()); // Find Track Bin for efficiency correction
          int primVtxCndtTag = 0;

          executePrimVtxCndtInCollisionloop<analysisType, kFillPreSel>(primVtxCndt, posTrack, negTrack,
                                                                       ptEtaBinPosTrk, ptEtaBinNegTrk, primVtxCndtTag, primVtxCndtDauList,
                                                                       idMethodSignTrk, checkTrackId, posDauIs, negDauIs,
                                                                       iNTrk, fNTrk, effWeightSum,
                                                                       posTrackPrimVtxMotherFlag, negTrackPrimVtxMotherFlag, posTrackV0MotherFlag, negTrackV0MotherFlag, 0, 0,
                                                                       doPhi1020, doJPsiToEE, doJPsiToMuMu, doKStar892, doKStar892Bar, doRho770,
                                                                       primVtxCndtMcTag);
        }
        anlysisTimingInfo.fill(HIST("time05_primVtxLoop1Time"), calculateTime(primVtxLoop1Start));

        // If all calculated particles are zero, skip unnecessary code and jump to the revelant part;
        iTotalPrimVtxCndtTrk = iNTrk[kPrimPhi1020][kPos] + iNTrk[kPrimJPsiToEE][kPos] + iNTrk[kPrimJPsiToMuMu][kPos] + iNTrk[kPrimKStar892][kPos] + iNTrk[kPrimKStar892Bar][kPos] + iNTrk[kPrimRho770][kPos];
        fTotalPrimVtxCndtTrk = fNTrk[kPrimPhi1020][kPos] + fNTrk[kPrimJPsiToEE][kPos] + fNTrk[kPrimJPsiToMuMu][kPos] + fNTrk[kPrimKStar892][kPos] + fNTrk[kPrimKStar892Bar][kPos] + fNTrk[kPrimRho770][kPos];
        if (iTotalPrimVtxCndtTrk == 0 && fTotalPrimVtxCndtTrk == 0) {
          goto skipSecondPrimVtxCndtLoop;
        }

        primVtxSortStart = std::chrono::high_resolution_clock::now();
        executeSortPairDaughters<false>(primVtxCndtDauList[kPrimTrkPhi1020][kPos], primVtxCndtDauList[kPrimTrkPhi1020][kNeg], primVtxTLHs[0].get<TH1>(HIST(HistRegDire[Phi1020]) + HIST(FillModeDire[kFillPreSel]) + HIST("nCommonKaOfDifferentPhi")));
        executeSortPairDaughters<false>(primVtxCndtDauList[kPrimTrkJPsiToEE][kPos], primVtxCndtDauList[kPrimTrkJPsiToEE][kNeg], primVtxTLHs[1].get<TH1>(HIST(HistRegDire[JPsiToEE]) + HIST(FillModeDire[kFillPreSel]) + HIST("nCommonKaOfDifferentJPsiToEE")));
        // To be updated in future : executeSortPairDaughters<false>(primVtxCndtDauList[kPrimTrkJPsiToMuMu ][kPos] , primVtxCndtDauList[kPrimTrkJPsiToMuMu ][kNeg], recoJPsiToMuMu.get<TH1>(HIST(HistRegDire[JPsiToMuMu])+ HIST(FillModeDire[kFillPreSel]) + HIST("nCommonKaOfDifferentJPsiToMuMu")));
        executeSortPairDaughters<false>(primVtxCndtDauList[kPrimTrkKStar892][kPos], primVtxCndtDauList[kPrimTrkKStar892][kNeg], primVtxTLHs[3].get<TH1>(HIST(HistRegDire[KStar892]) + HIST(FillModeDire[kFillPreSel]) + HIST("nCommonKaOfDifferentKStar892")));
        executeSortPairDaughters<false>(primVtxCndtDauList[kPrimTrkKStar892Bar][kPos], primVtxCndtDauList[kPrimTrkKStar892Bar][kNeg], primVtxTLHs[4].get<TH1>(HIST(HistRegDire[KStar892Bar]) + HIST(FillModeDire[kFillPreSel]) + HIST("nCommonKaOfDifferentKStar892Bar")));
        executeSortPairDaughters<false>(primVtxCndtDauList[kPrimTrkRho770][kPos], primVtxCndtDauList[kPrimTrkRho770][kNeg], primVtxTLHs[5].get<TH1>(HIST(HistRegDire[Rho770]) + HIST(FillModeDire[kFillPreSel]) + HIST("nCommonKaOfDifferentRho770")));

        count1 = iNTrk[kPrimPhi1020][kPos];
        count2 = iNTrk[kPrimJPsiToEE][kPos];
        count4 = iNTrk[kPrimKStar892][kPos];
        count5 = iNTrk[kPrimKStar892Bar][kPos];
        count6 = iNTrk[kPrimRho770][kPos];

        if (iNTrk[kPrimPhi1020][kNeg] != 0) {
          LOG(fatal) << "DEBUG :: fatal Error :: something wrong";
        }
        if (iNTrk[kPrimJPsiToEE][kNeg] != 0) {
          LOG(fatal) << "DEBUG :: fatal Error :: something wrong";
        }
        if (iNTrk[kPrimKStar892][kNeg] != 0) {
          LOG(fatal) << "DEBUG :: fatal Error :: something wrong";
        }
        if (iNTrk[kPrimKStar892Bar][kNeg] != 0) {
          LOG(fatal) << "DEBUG :: fatal Error :: something wrong";
        }
        if (iNTrk[kPrimRho770][kNeg] != 0) {
          LOG(fatal) << "DEBUG :: fatal Error :: something wrong";
        }

        anlysisTimingInfo.fill(HIST("time06_primVtxSortTime"), calculateTime(primVtxSortStart));
        primVtxLoop2Start = std::chrono::high_resolution_clock::now();
        // Second loop for proper counting i.e contamination removed counting.
        for (int i = kPrimPhi1020; i <= kPrimRho770; i++) {
          iNTrk[i][kPos] = 0;
          fNTrk[i][kPos] = 0;
          effWeightSum[i][kPos] = 0;
          iNTrk[i][kNeg] = 0;
          fNTrk[i][kNeg] = 0;
          effWeightSum[i][kNeg] = 0;
        }

        // Second Prim Vtx Loop
        for (const auto& primVtxCndt : primVtxCndtsPerColl) {
          const auto& posTrack = primVtxCndt.template posTrack_as<T>();
          const auto& negTrack = primVtxCndt.template negTrack_as<T>();
          if (!checkTrackSelection(posTrack, rejectionTag)) {
            continue;
          }
          if (!checkTrackSelection(negTrack, rejectionTag)) {
            continue;
          }

          primVtxCndtMcTag = 0xFF; // By Default all bits are on, if Mc process is running the Tag will be recalculated.
          //__________________________Reco Level _______________________________________________________________
          if constexpr (analysisType == doRecoProcessing || analysisType == doPurityProcessing) {

            if (!posTrack.has_mcParticle() || !negTrack.has_mcParticle()) {
              continue; // Reco Level Check // Fix in Future :: !primVtxCndt.has_mcParticle() is not implement for now.
            }

            auto pvcMcParticle = nullptr; // To correct this in future
            //-->Fix This ::  auto pvcMcParticle = v0.mcParticle(); //Reco Level Check
            //-->Fix This ::  if (!pvcMcParticle.isPhysicalPrimary())
            //-->Fix This ::    continue;
            //

            //___________________________Truth Level_________________________________________________________
            if constexpr (analysisType == doPurityProcessing) {
              auto posDauMcPart = posTrack.mcParticle();
              auto negDauMcPart = negTrack.mcParticle();
              getPvcMcBitTag<false>(primVtxCndtMcTag, pvcMcParticle, posDauMcPart, negDauMcPart);
              if (countBits(primVtxCndtMcTag) > 1) {
                LOG(fatal) << "DEBUG :: countBits(primVtxCndtMcTag) > 1  :: primVtxCndtMcTag = " << std::bitset<8>(primVtxCndtMcTag);
              }
            } //____Truth_Level_Over______
          } //____Reco_Level_Over______

          int ptEtaBinPosTrk = hPtEtaForBinSearch->FindBin(posTrack.pt(), posTrack.eta()); // Find Track Bin for efficiency correction
          int ptEtaBinNegTrk = hPtEtaForBinSearch->FindBin(negTrack.pt(), negTrack.eta()); // Find Track Bin for efficiency correction
          int primVtxCndtTag = 0;

          contaminatoinFlagCheck<true>(primVtxCndtDauList, posTrack, negTrack, posTrackPrimVtxMotherFlag, negTrackPrimVtxMotherFlag, massSigma[0], mpBit[0]);
          contaminatoinFlagCheck<false>(v0CndtDauList, posTrack, negTrack, posTrackV0MotherFlag, negTrackV0MotherFlag, massSigma[1], mpBit[1]);

          requiredBit = (massSigma[0] < massSigma[1]) ? 0 : 1; // 0 for primVtxCndt, 1 for v0DecayCndt;
          // If requiredBit is 1, daughter track likely belongs to a V0, not this primCandt

          executePrimVtxCndtInCollisionloop<analysisType, kFillPostSel>(primVtxCndt, posTrack, negTrack,
                                                                        ptEtaBinPosTrk, ptEtaBinNegTrk, primVtxCndtTag, primVtxCndtDauList,
                                                                        idMethodSignTrk, checkTrackId, posDauIs, negDauIs,
                                                                        iNTrk, fNTrk, effWeightSum,
                                                                        posTrackPrimVtxMotherFlag, negTrackPrimVtxMotherFlag, posTrackV0MotherFlag, negTrackV0MotherFlag, requiredBit, mpBit[requiredBit],
                                                                        doPhi1020, doJPsiToEE, doJPsiToMuMu, doKStar892, doKStar892Bar, doRho770, primVtxCndtMcTag);
        } // PrimVtx Particle counting is over;
        anlysisTimingInfo.fill(HIST("time07_primVtxLoop2Time"), calculateTime(primVtxLoop2Start));

        if (count1 != iNTrk[kPrimPhi1020][kNeg] + iNTrk[kPrimPhi1020][kPos]) {
          LOG(fatal) << "DEBUG :: something failed in second prim vtx candidate loop";
        }
        if (count2 != iNTrk[kPrimJPsiToEE][kNeg] + iNTrk[kPrimJPsiToEE][kPos]) {
          LOG(fatal) << "DEBUG :: something failed in second prim vtx candidate loop";
        }
        if (count4 != iNTrk[kPrimKStar892][kNeg] + iNTrk[kPrimKStar892][kPos]) {
          LOG(fatal) << "DEBUG :: something failed in second prim vtx candidate loop";
        }
        if (count5 != iNTrk[kPrimKStar892Bar][kNeg] + iNTrk[kPrimKStar892Bar][kPos]) {
          LOG(fatal) << "DEBUG :: something failed in second prim vtx candidate loop";
        }
        if (count6 != iNTrk[kPrimRho770][kNeg] + iNTrk[kPrimRho770][kPos]) {
          LOG(fatal) << "DEBUG :: something failed in second prim vtx candidate loop";
        }

      skipSecondPrimVtxCndtLoop:

        secondV0LoopStart = std::chrono::high_resolution_clock::now();
        // Intialise V0 related variables to zero
        for (int i = kV0K0s; i <= kV0Gamma; i++) {
          iNTrk[i][kPos] = 0;
          fNTrk[i][kPos] = 0;
          effWeightSum[i][kPos] = 0;
          iNTrk[i][kNeg] = 0;
          fNTrk[i][kNeg] = 0;
          effWeightSum[i][kNeg] = 0;
        }

        for (const auto& v0 : v0sTablePerColl) {
          const auto& posDaughterTrack = v0.template posTrack_as<T>();
          const auto& negDaughterTrack = v0.template negTrack_as<T>();

          v0DecayTrueMcTag = 0xFF; // By Default all bits are on, if Mc process is running the Tag will be recalculated.
          //__________________________Reco Level _______________________________________________________________
          if constexpr (analysisType == doRecoProcessing || analysisType == doPurityProcessing) {

            if (!v0.has_mcParticle() || !posDaughterTrack.has_mcParticle() || !negDaughterTrack.has_mcParticle()) {
              continue; // Reco Level Check
            }

            auto v0mcparticle = v0.mcParticle(); // Reco Level Check
            if (!v0mcparticle.isPhysicalPrimary())
              continue;

            //___________________________Truth Level_________________________________________________________
            if constexpr (analysisType == doPurityProcessing) {
              auto posDauMcPart = posDaughterTrack.mcParticle();
              auto negDauMcPart = negDaughterTrack.mcParticle();
              getV0McBitTag(v0DecayTrueMcTag, v0mcparticle, posDauMcPart, negDauMcPart);
              if (countBits(v0DecayTrueMcTag) > 1) {
                LOG(fatal) << "DEBUG :: countBits(v0DecayTrueMcTag) > 1  :: v0DecayTrueMcTag = " << std::bitset<8>(v0DecayTrueMcTag);
              }
            } //____Truth_Level_Over______
          } //____Reco_Level_Over______

          contaminatoinFlagCheck<true>(primVtxCndtDauList, posDaughterTrack, negDaughterTrack, posTrackPrimVtxMotherFlag, negTrackPrimVtxMotherFlag, massSigma[0], mpBit[0]);
          contaminatoinFlagCheck<false>(v0CndtDauList, posDaughterTrack, negDaughterTrack, posTrackV0MotherFlag, negTrackV0MotherFlag, massSigma[1], mpBit[1]);

          requiredBit = (massSigma[0] < massSigma[1]) ? 0 : 1; // 0 for primVtxCndt, 1 for v0DecayCndt;
          // If requiredBit is 0, daughter track likely belongs to a primVtxCndt, not this V0

          executeV0InCollisionloop<analysisType, C, kFillPostSel>(v0, posDaughterTrack, negDaughterTrack,
                                                                  idMethodSignTrk[kPi], idMethodSignTrk[kPr], idMethodSignTrk[kEl],
                                                                  v0Tag, v0DauCollisionIndexTag, v0DauBCTag, iNTrk, fNTrk, effWeightSum, centrality,
                                                                  posTrackPrimVtxMotherFlag, negTrackPrimVtxMotherFlag, posTrackV0MotherFlag, negTrackV0MotherFlag, requiredBit, mpBit[requiredBit],
                                                                  v0DecayTrueMcTag);
        } // End of V0s Loop
        anlysisTimingInfo.fill(HIST("time08_seconV0LoopTime"), calculateTime(secondV0LoopStart));

        trackLoopStart = std::chrono::high_resolution_clock::now();

        // Track count related variables
        for (int iPid = kPi; iPid <= kDe; iPid++) {
          iNTrk[iPid][kPos] = 0;
          iNTrk[iPid][kNeg] = 0;
          fNTrk[iPid][kPos] = 0;
          fNTrk[iPid][kNeg] = 0;
          effWeightSum[iPid][kPos] = 0;
          effWeightSum[iPid][kNeg] = 0;
        }

        // Rejected Track count
        for (int iPid = kRejectedPi; iPid <= kRejectedDe; iPid++) {
          iNTrk[iPid][kPos] = 0;
          iNTrk[iPid][kNeg] = 0;
          fNTrk[iPid][kPos] = 0;
          fNTrk[iPid][kNeg] = 0;
        }

        nTrack = 0;
        for (const auto& track : tracksTablePerColl) {

          trackTrueMcTag = 0xFF;
          //_________________________Reco___Level___________________________________________
          if constexpr (analysisType == doRecoProcessing || analysisType == doPurityProcessing) {
            if (!track.has_mcParticle()) {
              LOG(warning) << "No MC Particle for this track, skip...";
              continue;
            }

            auto mcPart = track.mcParticle();
            if (!mcPart.isPhysicalPrimary()) {
              continue;
            }
            if constexpr (analysisType == doPurityProcessing) {
              getTrackMcBitTag(trackTrueMcTag, mcPart);
            }
          }

          isAcceptedTrack = true; // by default track is accepted

          executeTrackQAPart(track, rejectionTag, nTrack, isAcceptedTrack);
          if (!isAcceptedTrack) {
            continue;
          }

          // Do Proper Track Identification
          for (int iPid = kPi; iPid <= kDe; iPid++) {
            trackIs[iPid] = false;
            idMethodTrk[iPid] = kUnidentified;
          }

          trackIdTag = 0;

          ptEtaBinTrk = hPtEtaForBinSearch->FindBin(track.pt(), track.eta()); // Find Track Bin for efficiency correction

          if (selPion(track, idMethodTrk[kPi]) && isTrueMcMatch<analysisType>(trackTrueMcTag, ID_BIT_PI)) {
            trackIs[kPi] = true;
            BITSET(trackIdTag, ID_BIT_PI);
          }
          if (selKaon(track, idMethodTrk[kKa]) && isTrueMcMatch<analysisType>(trackTrueMcTag, ID_BIT_KA)) {
            trackIs[kKa] = true;
            BITSET(trackIdTag, ID_BIT_KA);
          }
          if (selProton(track, idMethodTrk[kPr]) && isTrueMcMatch<analysisType>(trackTrueMcTag, ID_BIT_PR)) {
            trackIs[kPr] = true;
            BITSET(trackIdTag, ID_BIT_PR);
          }
          if (selElectron(track, idMethodTrk[kEl]) && isTrueMcMatch<analysisType>(trackTrueMcTag, ID_BIT_EL)) {
            trackIs[kEl] = true;
            BITSET(trackIdTag, ID_BIT_EL);
          }
          if (selMuon(track, idMethodTrk[kMu]) && isTrueMcMatch<analysisType>(trackTrueMcTag, ID_BIT_MU)) {
            trackIs[kMu] = true;
            BITSET(trackIdTag, ID_BIT_MU);
          }
          if (selDeuteron(track, idMethodTrk[kDe]) && isTrueMcMatch<analysisType>(trackTrueMcTag, ID_BIT_DE)) {
            trackIs[kDe] = true;
            BITSET(trackIdTag, ID_BIT_DE);
          }

          getTrackDecayInfoBit(track, decayDauTagBit, v0CndtDauList, primVtxCndtDauList,
                               doV0K0s, doV0Lambda, doV0AntiLambda, doV0Gamma,
                               doPhi1020, doJPsiToEE, doJPsiToMuMu, doKStar892, doKStar892Bar, doRho770);

          executeTrackAnalysisPart(track, trackIdTag, ptEtaBinTrk, decayDauTagBit, idMethodTrk, trackIs,
                                   iNTrk[kPi], fNTrk[kPi], effWeight[kPi], effWeightSum[kPi], iNTrk[kRejectedPi], fNTrk[kRejectedPi],
                                   iNTrk[kKa], fNTrk[kKa], effWeight[kKa], effWeightSum[kKa], iNTrk[kRejectedKa], fNTrk[kRejectedKa],
                                   iNTrk[kPr], fNTrk[kPr], effWeight[kPr], effWeightSum[kPr], iNTrk[kRejectedPr], fNTrk[kRejectedPr],
                                   iNTrk[kEl], fNTrk[kEl], effWeight[kEl], effWeightSum[kEl], iNTrk[kRejectedEl], fNTrk[kRejectedEl],
                                   iNTrk[kDe], fNTrk[kDe], effWeight[kDe], effWeightSum[kDe], iNTrk[kRejectedDe], fNTrk[kRejectedDe]);
        } // track loop ends
        anlysisTimingInfo.fill(HIST("time09_trackLoopTime"), calculateTime(trackLoopStart));

        iNTrk[kNKaon][kPos] = iNTrk[kKa][kPos] + iNTrk[kKa][kNeg];
        iNTrk[kNK0sSq][kPos] = iNTrk[kV0K0s][kPos] * iNTrk[kV0K0s][kPos];
        iNTrk[kNKaonSq][kPos] = iNTrk[kNKaon][kPos] * iNTrk[kNKaon][kPos];
        iNTrk[kNK0sNKaonProd][kPos] = iNTrk[kV0K0s][kPos] * iNTrk[kNKaon][kPos];
        iNTrk[kNKaPosSq][kPos] = iNTrk[kKa][kPos] * iNTrk[kKa][kPos];
        iNTrk[kNKaNegSq][kPos] = iNTrk[kKa][kNeg] * iNTrk[kKa][kNeg];
        iNTrk[kNKaPosNKaNegProd][kPos] = iNTrk[kKa][kPos] * iNTrk[kKa][kNeg];

        fNTrk[kNKaon][kPos] = fNTrk[kKa][kPos] + fNTrk[kKa][kNeg];
        fNTrk[kNK0sSq][kPos] = fNTrk[kV0K0s][kPos] * fNTrk[kV0K0s][kPos];
        fNTrk[kNKaonSq][kPos] = fNTrk[kNKaon][kPos] * fNTrk[kNKaon][kPos];
        fNTrk[kNK0sNKaonProd][kPos] = fNTrk[kV0K0s][kPos] * fNTrk[kNKaon][kPos];
        fNTrk[kNKaPosSq][kPos] = fNTrk[kKa][kPos] * fNTrk[kKa][kPos];
        fNTrk[kNKaNegSq][kPos] = fNTrk[kKa][kNeg] * fNTrk[kKa][kNeg];
        fNTrk[kNKaPosNKaNegProd][kPos] = fNTrk[kKa][kPos] * fNTrk[kKa][kNeg];

        sparseFillStart = std::chrono::high_resolution_clock::now();
        fillAllSparses(hSparses, hSparseVarPtrs, doubleCasters, buffers);
        anlysisTimingInfo.fill(HIST("time10_hSparseFillTime"), calculateTime(sparseFillStart));

        fillEventInfoStart = std::chrono::high_resolution_clock::now();
        recoAnalysis.fill(HIST("recoAnalysis/centFT0C_vs_NTrk"), collision.centFT0C(), nTrack);

        if constexpr (analysisType == doDataProcessing) {
          recoAnalysis.fill(HIST("recoAnalysis/centFT0A_vs_NTrk"), collision.centFT0A(), nTrack);
          recoAnalysis.fill(HIST("recoAnalysis/centFT0M_vs_NTrk"), collision.centFT0M(), nTrack);

          recoAnalysis.fill(HIST("recoAnalysis/centFT0A_vs_centFT0M"), collision.centFT0A(), collision.centFT0M());
          recoAnalysis.fill(HIST("recoAnalysis/centFT0C_vs_centFT0A"), collision.centFT0C(), collision.centFT0A());
          recoAnalysis.fill(HIST("recoAnalysis/centFT0M_vs_centFT0C"), collision.centFT0M(), collision.centFT0C());
        }

        executeEventInfoPart(collision, centrality, v0sTablePerColl.size(), tracksTablePerColl.size(),
                             nTrack, iNTrk, fNTrk, effWeightSum);

        anlysisTimingInfo.fill(HIST("time11_EventInfoFillTime"), calculateTime(fillEventInfoStart));
        fillCollQA<eventPostSel>(collision);
      } // collision loop ends
      anlysisTimingInfo.fill(HIST("time03_collisionLoopTime"), calculateTime(collLoopStart));
      anlysisTimingInfo.fill(HIST("time00_dfProcessingTime"), calculateTime(start1));
    } else if constexpr (analysisType == doGenProcessing || analysisType == doSimProcessing) {

      auto finalParticleIdList = (std::vector<int>)cfgSim.cfgSim07FinalParticleIdList;
      auto nonFinalParticleIdList = (std::vector<int>)cfgSim.cfgSim08NonFinalParticleIdList;
      if constexpr (analysisType == doSimProcessing) {
        std::sort(finalParticleIdList.begin(), finalParticleIdList.end());
        std::sort(nonFinalParticleIdList.begin(), nonFinalParticleIdList.end());
      }

      MyMcCollisions::iterator mcColl;
      const auto& mcParticles = tracks; // In Gen and Sim, tracks are actually mcParticles from mcParticle table
      for (const auto& collision : collisions) {

        if constexpr (analysisType == doGenProcessing) {
          // In Gen processing the collision not mcCollision, it is aod::collision joined with mc collisoin labels.
          if (!collision.has_mcCollision()) {
            continue;
          }
          centrality = -1;
          mcColl = collision.mcCollision();

          centrality = collision.centFT0C();
          if (cfgCentAxis.centAxis04Type == kCentFT0M) {
            centrality = collision.centFT0M();
          } else if (cfgCentAxis.centAxis04Type == kMultFT0M) {
            centrality = collision.multFT0M();
          } else if (cfgCentAxis.centAxis04Type == kMultFT0C) {
            centrality = collision.multFT0C();
          }
        } else if constexpr (analysisType == doSimProcessing) {
          // In sim processing the collision is itself montecarlo collision;
          mcColl = collision;
          if (cfgSim.cfgSim01VtxZCheck) {
            if (std::abs(mcColl.posZ()) >= cutZvertex) {
              continue;
            }
          }
        }

        // group over mcParticles
        const auto mcTracksTablePerMcColl = mcParticles.sliceBy(mcTracksPerMcCollisionPreslice, mcColl.globalIndex());

        for (int iPid = 0; iPid < pidVecSize; iPid++) {
          iNTrk[iPid][kPos] = 0;
          iNTrk[iPid][kNeg] = 0;
          fNTrk[iPid][kPos] = 0;
          fNTrk[iPid][kNeg] = 0;
          effWeightSum[iPid][kPos] = 0;
          effWeightSum[iPid][kNeg] = 0;
        }

        int multV0M = 0;
        for (const auto& mcTrack : mcTracksTablePerMcColl) {

          //______________Simulation Part__________________________________________________
          if constexpr (analysisType == doSimProcessing) {
            if (cfgSim.cfgSim02CountFinalParticles) {
              if (!mcTrack.has_daughters() && std::binary_search(finalParticleIdList.begin(), finalParticleIdList.end(), mcTrack.pdgCode())) {
                getV0MCount(mcTrack, multV0M);
              }
            }
            if (cfgSim.cfgSim03CountNonFinalParticles) {
              if (mcTrack.has_daughters() && std::binary_search(nonFinalParticleIdList.begin(), nonFinalParticleIdList.end(), mcTrack.pdgCode())) {
                getV0MCount(mcTrack, multV0M);
              }
            }

            if (cfgSim.cfgSim04CountPhysicalPrimAndFinalParticles) {
              if (!mcTrack.has_daughters() && mcTrack.isPhysicalPrimary()) {
                if (!(std::abs(mcTrack.pdgCode()) == kNuE || std::abs(mcTrack.pdgCode()) == kNuMu || std::abs(mcTrack.pdgCode()) == kNuTau)) {
                  // Removed invisible neutrinos;
                  getV0MCount(mcTrack, multV0M);
                }
              }
            }

            centrality = multV0M;
          } //______________Simulation Related Part is over__________________________________________________

          if (!mcTrack.isPhysicalPrimary()) {
            continue;
          }

          if (mcTrack.pdgCode() == kK0Short &&
              getCfg<float>(cfgV0ParticleCuts, kV0TrkK0s, kV0LowPt) < mcTrack.pt() && mcTrack.pt() < getCfg<float>(cfgV0ParticleCuts, kV0TrkK0s, kV0HighPt) &&
              std::abs(mcTrack.y()) < getCfg<float>(cfgV0ParticleCuts, kV0TrkK0s, kV0Rapitidy)) {
            iNTrk[kV0K0s][kPos]++;
            fillGenTrackQA<genAnalysisDir, kV0K0s>(genAnalysis, mcTrack);
          }

          if (mcTrack.pt() <= cfgTrackCuts.cfgTrk05PtLow || mcTrack.pt() >= cfgTrackCuts.cfgTrk06PtHigh || std::abs(mcTrack.eta()) >= cfgTrackCuts.cfgTrk04Eta) {
            continue; // For track loop
          }

          if (mcTrack.pdgCode() == kPiPlus) {
            fillGenTrackQA<genAnalysisDir, kPi>(genAnalysis, mcTrack);
            iNTrk[kPi][kPos]++;
          } else if (mcTrack.pdgCode() == kPiMinus) {
            fillGenTrackQA<genAnalysisDir, kPi>(genAnalysis, mcTrack);
            iNTrk[kPi][kNeg]++;
          } else if (mcTrack.pdgCode() == kKPlus) {
            fillGenTrackQA<genAnalysisDir, kKa>(genAnalysis, mcTrack);
            iNTrk[kKa][kPos]++;
          } else if (mcTrack.pdgCode() == kKMinus) {
            fillGenTrackQA<genAnalysisDir, kKa>(genAnalysis, mcTrack);
            iNTrk[kKa][kNeg]++;
          } else if (mcTrack.pdgCode() == kProton) {
            fillGenTrackQA<genAnalysisDir, kPr>(genAnalysis, mcTrack);
            iNTrk[kPr][kPos]++;
          } else if (mcTrack.pdgCode() == kProtonBar) {
            fillGenTrackQA<genAnalysisDir, kPr>(genAnalysis, mcTrack);
            iNTrk[kPr][kNeg]++;
          } else if (mcTrack.pdgCode() == kElectron) {
            fillGenTrackQA<genAnalysisDir, kEl>(genAnalysis, mcTrack);
            iNTrk[kEl][kNeg]++;
          } else if (mcTrack.pdgCode() == kPositron) {
            fillGenTrackQA<genAnalysisDir, kEl>(genAnalysis, mcTrack);
            iNTrk[kEl][kPos]++;
          } else if (mcTrack.pdgCode() == kDeuteron) {
            fillGenTrackQA<genAnalysisDir, kDe>(genAnalysis, mcTrack);
            iNTrk[kDe][kPos]++;
          } else if (mcTrack.pdgCode() == -kDeuteron) {
            fillGenTrackQA<genAnalysisDir, kDe>(genAnalysis, mcTrack);
            iNTrk[kDe][kNeg]++;
          }
          nTrack++;
        } // mcTrack loop is over

        iNTrk[kNKaon][kPos] = iNTrk[kKa][kPos] + iNTrk[kKa][kNeg];
        iNTrk[kNK0sSq][kPos] = iNTrk[kV0K0s][kPos] * iNTrk[kV0K0s][kPos];
        iNTrk[kNKaonSq][kPos] = iNTrk[kNKaon][kPos] * iNTrk[kNKaon][kPos];
        iNTrk[kNK0sNKaonProd][kPos] = iNTrk[kV0K0s][kPos] * iNTrk[kNKaon][kPos];
        iNTrk[kNKaPosSq][kPos] = iNTrk[kKa][kPos] * iNTrk[kKa][kPos];
        iNTrk[kNKaNegSq][kPos] = iNTrk[kKa][kNeg] * iNTrk[kKa][kNeg];
        iNTrk[kNKaPosNKaNegProd][kPos] = iNTrk[kKa][kPos] * iNTrk[kKa][kNeg];

        sparseFillStart = std::chrono::high_resolution_clock::now();
        fillAllSparses(hSparses, hSparseVarPtrs, doubleCasters, buffers);
        anlysisTimingInfo.fill(HIST("time10_hSparseFillTime"), calculateTime(sparseFillStart));

        executeEventInfoPart(mcColl, centrality, 0, 0,
                             nTrack, iNTrk, fNTrk, effWeightSum);
      } // collision loop is over
      anlysisTimingInfo.fill(HIST("time03_collisionLoopTime"), calculateTime(collLoopStart));
      anlysisTimingInfo.fill(HIST("time00_dfProcessingTime"), calculateTime(start1));
    }

    // Debug messages for checking time consuming parts
    if (cfgDebug.printDebugMessages) {
      printHistInfo(Form("DEBUG :: df_%lld :: time00_dfProcessingTime  = ", dfCount), anlysisTimingInfo.get<TH1>(HIST("time00_dfProcessingTime")));
      printHistInfo(Form("DEBUG :: df_%lld :: time01_executeV0loop     = ", dfCount), anlysisTimingInfo.get<TH1>(HIST("time01_executeV0loop")));
      printHistInfo(Form("DEBUG :: df_%lld :: time02_v0SortTime        = ", dfCount), anlysisTimingInfo.get<TH1>(HIST("time02_v0SortTime")));
      printHistInfo(Form("DEBUG :: df_%lld :: time03_collisionLoopTime = ", dfCount), anlysisTimingInfo.get<TH1>(HIST("time03_collisionLoopTime")));
      printHistInfo(Form("DEBUG :: df_%lld :: time04_manualGrouping    = ", dfCount), anlysisTimingInfo.get<TH1>(HIST("time04_manualGrouping")));
      printHistInfo(Form("DEBUG :: df_%lld :: time05_primVtxLoop1Time  = ", dfCount), anlysisTimingInfo.get<TH1>(HIST("time05_primVtxLoop1Time")));
      printHistInfo(Form("DEBUG :: df_%lld :: time06_primVtxSortTime   = ", dfCount), anlysisTimingInfo.get<TH1>(HIST("time06_primVtxSortTime")));
      printHistInfo(Form("DEBUG :: df_%lld :: time07_primVtxLoop2Time  = ", dfCount), anlysisTimingInfo.get<TH1>(HIST("time07_primVtxLoop2Time")));
      printHistInfo(Form("DEBUG :: df_%lld :: time08_seconV0LoopTime   = ", dfCount), anlysisTimingInfo.get<TH1>(HIST("time08_seconV0LoopTime")));
      printHistInfo(Form("DEBUG :: df_%lld :: time09_trackLoopTime     = ", dfCount), anlysisTimingInfo.get<TH1>(HIST("time09_trackLoopTime")));
      printHistInfo(Form("DEBUG :: df_%lld :: time10_hSparseFillTime   = ", dfCount), anlysisTimingInfo.get<TH1>(HIST("time10_hSparseFillTime")));
      printHistInfo(Form("DEBUG :: df_%lld :: time11_EventInfoFillTime = ", dfCount), anlysisTimingInfo.get<TH1>(HIST("time11_EventInfoFillTime")));
      LOG(info) << "DEBUG ::";
    }
    if (cfgDebug.resetHistograms) {
      anlysisTimingInfo.get<TH1>(HIST("time00_dfProcessingTime"))->Reset();
      anlysisTimingInfo.get<TH1>(HIST("time01_executeV0loop"))->Reset();
      anlysisTimingInfo.get<TH1>(HIST("time02_v0SortTime"))->Reset();
      anlysisTimingInfo.get<TH1>(HIST("time03_collisionLoopTime"))->Reset();
      anlysisTimingInfo.get<TH1>(HIST("time04_manualGrouping"))->Reset();
      anlysisTimingInfo.get<TH1>(HIST("time05_primVtxLoop1Time"))->Reset();
      anlysisTimingInfo.get<TH1>(HIST("time06_primVtxSortTime"))->Reset();
      anlysisTimingInfo.get<TH1>(HIST("time07_primVtxLoop2Time"))->Reset();
      anlysisTimingInfo.get<TH1>(HIST("time08_seconV0LoopTime"))->Reset();
      anlysisTimingInfo.get<TH1>(HIST("time09_trackLoopTime"))->Reset();
      anlysisTimingInfo.get<TH1>(HIST("time10_hSparseFillTime"))->Reset();
      anlysisTimingInfo.get<TH1>(HIST("time11_EventInfoFillTime"))->Reset();
    }
  } // end of execute analyis block

  //____________________________________Process Funtion For Analysis Starts Here____________________________________//

  void processData(aod::BCsWithTimestamps const&, MyCollisions const& collisions, MyV0s const& V0s, MyTracks const& tracks, MyPrimVtxCndts const& primVtxCndts)
  {
    recoEvent.fill(HIST("recoEvent/ProcessType"), doDataProcessing);
    executeAnalysis<doDataProcessing>(collisions, V0s, tracks, primVtxCndts);
  }
  PROCESS_SWITCH(KaonIsospinFluctuations, processData, "Process for Data", true);

  void processReco(aod::BCsWithTimestamps const&, MyCollisionsWithMcLabels const& collisions, MyV0sWithMcLabels const& V0s, MyTracksWithMcLabels const& tracks, MyPrimVtxCndts const& primVtxCndts, aod::McParticles const&)
  {
    recoEvent.fill(HIST("recoEvent/ProcessType"), doRecoProcessing);
    executeAnalysis<doRecoProcessing>(collisions, V0s, tracks, primVtxCndts);
  }
  PROCESS_SWITCH(KaonIsospinFluctuations, processReco, "Process for Reco", false);

  void processPurity(aod::BCsWithTimestamps const&, MyCollisionsWithMcLabels const& collisions, MyV0sWithMcLabels const& V0s, MyTracksWithMcLabels const& tracks, MyPrimVtxCndts const& primVtxCndts, aod::McParticles const&)
  {
    recoEvent.fill(HIST("recoEvent/ProcessType"), doPurityProcessing);
    executeAnalysis<doPurityProcessing>(collisions, V0s, tracks, primVtxCndts);
  }
  PROCESS_SWITCH(KaonIsospinFluctuations, processPurity, "Process for Purity", false);

  void processGen(MyMcCollisions const&, MyCollisionsWithMcLabels const& collisions, aod::McParticles const& mcParticles)
  {
    recoEvent.fill(HIST("recoEvent/ProcessType"), doGenProcessing);
    executeAnalysis<doGenProcessing>(collisions, nullptr, mcParticles, nullptr);
  }
  PROCESS_SWITCH(KaonIsospinFluctuations, processGen, "Process for Gen", false);

  void processSim(MyMcCollisions const& mcCollisions, aod::McParticles const& mcParticles)
  {
    recoEvent.fill(HIST("recoEvent/ProcessType"), doSimProcessing);
    executeAnalysis<doSimProcessing>(mcCollisions, nullptr, mcParticles, nullptr);
  }
  PROCESS_SWITCH(KaonIsospinFluctuations, processSim, "Process for Sim", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<PrimVtxParticleTable>(cfgc),
    adaptAnalysisTask<KaonIsospinFluctuations>(cfgc)};
}
