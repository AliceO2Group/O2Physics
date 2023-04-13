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

/// \file task-mini.cxx
/// \brief Mini version of the HF analysis chain
///
/// \author Vít Kučera <vit.kucera@cern.ch>, CERN

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "DetectorsVertexing/DCAFitterN.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelectorPID.h"
#include "PWGHF/Core/HFSelectorCuts.h"

#include <algorithm>

using namespace o2::analysis;

// Skimming =====================================================================

namespace o2::aod
{
namespace hf_seltrack
{
// Track selection columns
DECLARE_SOA_COLUMN(IsSelProng, isSelProng, int); //!
DECLARE_SOA_COLUMN(PxProng, pxProng, float);     //!
DECLARE_SOA_COLUMN(PyProng, pyProng, float);     //!
DECLARE_SOA_COLUMN(PzProng, pzProng, float);     //!
} // namespace hf_seltrack

// Track selection table
DECLARE_SOA_TABLE(HFSelTrack, "AOD", "HFSELTRACK", //!
                  hf_seltrack::IsSelProng,
                  hf_seltrack::PxProng,
                  hf_seltrack::PyProng,
                  hf_seltrack::PzProng);

using BigTracks = soa::Join<Tracks, TracksCov, TracksExtra, HFSelTrack>;
using BigTracksExtended = soa::Join<BigTracks, aod::TracksExtended>;
using BigTracksPID = soa::Join<BigTracks,
                               aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                               aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
using BigTracksPIDExtended = soa::Join<BigTracksPID, aod::TracksExtended>;

namespace hf_track_index
{
// Track index skim columns
DECLARE_SOA_INDEX_COLUMN_FULL(Index0, index0, int, Tracks, "_0"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(Index1, index1, int, Tracks, "_1"); //!
} // namespace hf_track_index

// Track index skim table
DECLARE_SOA_TABLE(HfTrackIndexProng2, "AOD", "HFTRACKIDXP2", //!
                  hf_track_index::Index0Id,
                  hf_track_index::Index1Id);

} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

static const double massPi = RecoDecay::getMassPDG(kPiPlus);
static const double massK = RecoDecay::getMassPDG(kKPlus);
static const auto arrMassPiK = std::array{massPi, massK};
static const auto arrMassKPi = std::array{massK, massPi};

using TracksAll = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksExtended>;

/// Track selection
struct HfTagSelTracks {
  Produces<aod::HFSelTrack> rowSelectedTrack;

  // 2-prong cuts
  Configurable<double> pTTrackMin{"pTTrackMin", -1., "min. track pT for 2 prong candidate"};
  Configurable<double> etaTrackMax{"etaTrackMax", 4., "max. pseudorapidity for 2 prong candidate"};
  Configurable<double> dcaTrackMin{"dcaTrackMin", 0.0025, "min. DCA for 2 prong candidate"};

  HistogramRegistry registry{
    "registry",
    {{"hPtNoCuts", "all tracks;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     // 2-prong histograms
     {"hPtCuts2Prong", "tracks selected for 2-prong vertexing;#it{p}_{T}^{track} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hDCAToPrimXYVsPtCuts2Prong", "tracks selected for 2-prong vertexing;#it{p}_{T}^{track} (GeV/#it{c});DCAxy to prim. vtx. (cm);entries", {HistType::kTH2F, {{100, 0., 10.}, {400, -2., 2.}}}},
     {"hEtaCuts2Prong", "tracks selected for 2-prong vertexing;#it{#eta};entries", {HistType::kTH1F, {{static_cast<int>(1.2 * etaTrackMax * 100), -1.2 * etaTrackMax, 1.2 * etaTrackMax}}}}}};

  void process(aod::Collisions const& collisions,
               TracksAll const& tracks)
  {
    for (auto& track : tracks) {

      int statusProng = 1;

      auto trackPt = track.pt();
      auto trackEta = track.eta();

      registry.fill(HIST("hPtNoCuts"), trackPt);

      // pT cut
      if (trackPt < pTTrackMin) {
        statusProng = 0;
      }

      // eta cut
      if (statusProng > 0 && std::abs(trackEta) > etaTrackMax) {
        statusProng = 0;
      }

      // DCA cut
      auto dcaXY = track.dcaXY();
      if (statusProng > 0 && std::abs(dcaXY) < dcaTrackMin) {
        statusProng = 0;
      }

      // fill histograms
      if (statusProng > 0) {
        registry.fill(HIST("hPtCuts2Prong"), trackPt);
        registry.fill(HIST("hEtaCuts2Prong"), trackEta);
        registry.fill(HIST("hDCAToPrimXYVsPtCuts2Prong"), trackPt, dcaXY);
      }

      // fill table row
      rowSelectedTrack(statusProng, track.px(), track.py(), track.pz());
    }
  }
};

/// Track index skim creator
/// Pre-selection of 2-prong secondary vertices
struct HfTrackIndexSkimsCreator {
  Produces<aod::HfTrackIndexProng2> rowTrackIndexProng2;

  // vertexing parameters
  Configurable<double> magneticField{"magneticField", 5., "magnetic field [kG]"};
  Configurable<bool> propToDCA{"propToDCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};

  HistogramRegistry registry{
    "registry",
    {// 2-prong histograms
     {"hVtx2ProngX", "2-prong candidates;#it{x}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}}},
     {"hVtx2ProngY", "2-prong candidates;#it{y}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}}},
     {"hVtx2ProngZ", "2-prong candidates;#it{z}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -20., 20.}}}},
     {"hmassD0ToPiK", "D^{0} candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}}}}};

  Filter filterSelectTracks = aod::hf_seltrack::isSelProng > 0;

  using SelectedTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra, aod::TracksExtended, aod::HFSelTrack>>;

  void process(
    aod::Collision const& collision,
    aod::BCs const& bcs,
    SelectedTracks const& tracks)
  {
    // 2-prong vertex fitter
    o2::vertexing::DCAFitterN<2> df2;
    df2.setBz(magneticField);
    df2.setPropagateToPCA(propToDCA);
    df2.setMaxR(maxR);
    df2.setMaxDZIni(maxDZIni);
    df2.setMinParamChange(minParamChange);
    df2.setMinRelChi2Change(minRelChi2Change);
    df2.setUseAbsDCA(useAbsDCA);

    // first loop over positive tracks
    for (auto& trackPos1 : tracks) {
      if (trackPos1.signed1Pt() < 0) {
        continue;
      }
      if (!trackPos1.isSelProng()) {
        continue;
      }

      auto trackParVarPos1 = getTrackParCov(trackPos1);

      // first loop over negative tracks
      for (auto& trackNeg1 : tracks) {
        if (trackNeg1.signed1Pt() > 0) {
          continue;
        }
        if (!trackNeg1.isSelProng()) {
          continue;
        }

        auto trackParVarNeg1 = getTrackParCov(trackNeg1);

        // 2-prong preselections

        // auto arrMom = std::array{
        //   std::array{trackPos1.pxProng(), trackPos1.pyProng(), trackPos1.pzProng()},
        //   std::array{trackNeg1.pxProng(), trackNeg1.pyProng(), trackNeg1.pzProng()}};

        // auto pT = RecoDecay::Pt(arrMom[0], arrMom[1]);
        // auto massPiK = RecoDecay::M2(arrMom, arrMassPiK);

        // secondary vertex reconstruction and further 2-prong selections
        if (df2.process(trackParVarPos1, trackParVarNeg1) > 0) {
          // auto primaryVertex = std::array{collision.posX(), collision.posY(), collision.posZ()};
          //  get secondary vertex
          const auto& secondaryVertex = df2.getPCACandidate();
          // get track momenta
          array<float, 3> pvec0;
          array<float, 3> pvec1;
          df2.getTrack(0).getPxPyPzGlo(pvec0);
          df2.getTrack(1).getPxPyPzGlo(pvec1);

          // auto pVecCand = RecoDecay::PVec(pvec0, pvec1);
          // auto pTCand = RecoDecay::Pt(pVecCand);
          //  2-prong selections after secondary vertex
          // auto cpa = RecoDecay::CPA(primaryVertex, secondaryVertex, pVecCand);
          std::array<std::array<float, 3>, 2> arrMom = {pvec0, pvec1};

          // fill table row
          rowTrackIndexProng2(trackPos1.globalIndex(),
                              trackNeg1.globalIndex());

          // fill histograms
          registry.fill(HIST("hVtx2ProngX"), secondaryVertex[0]);
          registry.fill(HIST("hVtx2ProngY"), secondaryVertex[1]);
          registry.fill(HIST("hVtx2ProngZ"), secondaryVertex[2]);
          auto mass2Prong = RecoDecay::M(arrMom, arrMassPiK);
          registry.fill(HIST("hmassD0ToPiK"), mass2Prong);
        }
      }
    }
  }
};

// Candidate creation =====================================================================

namespace o2::aod
{
// 2-prong decay properties
namespace hf_cand_prong2
{
// Candidate columns
// collision properties
DECLARE_SOA_INDEX_COLUMN(Collision, collision); //!
// secondary vertex
DECLARE_SOA_COLUMN(XSecondaryVertex, xSecondaryVertex, float); //!
DECLARE_SOA_COLUMN(YSecondaryVertex, ySecondaryVertex, float); //!
DECLARE_SOA_COLUMN(ZSecondaryVertex, zSecondaryVertex, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(RSecondaryVertex, rSecondaryVertex, //!
                           [](float xVtxS, float yVtxS) -> float { return RecoDecay::sqrtSumOfSquares(xVtxS, yVtxS); });
// prong properties
DECLARE_SOA_COLUMN(PxProng0, pxProng0, float); //!
DECLARE_SOA_COLUMN(PyProng0, pyProng0, float); //!
DECLARE_SOA_COLUMN(PzProng0, pzProng0, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(PtProng0, ptProng0, //!
                           [](float px, float py) -> float { return RecoDecay::Pt(px, py); });
DECLARE_SOA_COLUMN(PxProng1, pxProng1, float); //!
DECLARE_SOA_COLUMN(PyProng1, pyProng1, float); //!
DECLARE_SOA_COLUMN(PzProng1, pzProng1, float); //!
DECLARE_SOA_DYNAMIC_COLUMN(PtProng1, ptProng1, //!
                           [](float px, float py) -> float { return RecoDecay::Pt(px, py); });
// candidate properties
DECLARE_SOA_DYNAMIC_COLUMN(DecayLength, decayLength, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS) -> float { return RecoDecay::distance(array{xVtxP, yVtxP, zVtxP}, array{xVtxS, yVtxS, zVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //!
                           [](float px, float py) -> float { return RecoDecay::Pt(px, py); });
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //!
                              float, 1.f * pxProng0 + 1.f * pxProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //!
                              float, 1.f * pyProng0 + 1.f * pyProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //!
                              float, 1.f * pzProng0 + 1.f * pzProng1);
DECLARE_SOA_DYNAMIC_COLUMN(M, m, //!
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const array<double, 2>& m) -> float { return RecoDecay::M(array{array{px0, py0, pz0}, array{px1, py1, pz1}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(CPA, cpa, //!
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::CPA(array{xVtxP, yVtxP, zVtxP}, array{xVtxS, yVtxS, zVtxS}, array{px, py, pz}); });

template <typename T>
auto InvMassD0(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(kPiPlus), RecoDecay::getMassPDG(kKPlus)});
}

template <typename T>
auto InvMassD0bar(const T& candidate)
{
  return candidate.m(array{RecoDecay::getMassPDG(kKPlus), RecoDecay::getMassPDG(kPiPlus)});
}
} // namespace hf_cand_prong2

// Candidate table
DECLARE_SOA_TABLE(HfCandProng2Base, "AOD", "HFCANDP2BASE", //!
                  hf_cand_prong2::CollisionId,
                  collision::PosX, collision::PosY, collision::PosZ,
                  hf_cand_prong2::XSecondaryVertex, hf_cand_prong2::YSecondaryVertex, hf_cand_prong2::ZSecondaryVertex,
                  /* dynamic columns */ hf_cand_prong2::RSecondaryVertex<hf_cand_prong2::XSecondaryVertex, hf_cand_prong2::YSecondaryVertex>,
                  hf_cand_prong2::DecayLength<collision::PosX, collision::PosY, collision::PosZ, hf_cand_prong2::XSecondaryVertex, hf_cand_prong2::YSecondaryVertex, hf_cand_prong2::ZSecondaryVertex>,
                  /* prong 0 */ hf_cand_prong2::PtProng0<hf_cand_prong2::PxProng0, hf_cand_prong2::PyProng0>,
                  hf_cand_prong2::PxProng0, hf_cand_prong2::PyProng0, hf_cand_prong2::PzProng0,
                  /* prong 1 */ hf_cand_prong2::PtProng1<hf_cand_prong2::PxProng1, hf_cand_prong2::PyProng1>,
                  hf_cand_prong2::PxProng1, hf_cand_prong2::PyProng1, hf_cand_prong2::PzProng1,
                  hf_track_index::Index0Id, hf_track_index::Index1Id,
                  /* dynamic columns */
                  hf_cand_prong2::M<hf_cand_prong2::PxProng0, hf_cand_prong2::PyProng0, hf_cand_prong2::PzProng0, hf_cand_prong2::PxProng1, hf_cand_prong2::PyProng1, hf_cand_prong2::PzProng1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand_prong2::CPA<collision::PosX, collision::PosY, collision::PosZ, hf_cand_prong2::XSecondaryVertex, hf_cand_prong2::YSecondaryVertex, hf_cand_prong2::ZSecondaryVertex, hf_cand_prong2::Px, hf_cand_prong2::Py, hf_cand_prong2::Pz>,
                  hf_cand_prong2::Pt<hf_cand_prong2::Px, hf_cand_prong2::Py>);

// extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandProng2Ext, HfCandProng2Base, "HFCANDP2EXT", //!
                                hf_cand_prong2::Px, hf_cand_prong2::Py, hf_cand_prong2::Pz);

using HfCandProng2 = HfCandProng2Ext;

} // namespace o2::aod

/// Candidate creator
/// Reconstruction of heavy-flavour 2-prong decay candidates
struct HfCandidateCreator2Prong {
  Produces<aod::HfCandProng2Base> rowCandidateBase;

  Configurable<double> magneticField{"magneticField", 5., "magnetic field [kG]"};
  Configurable<bool> propToDCA{"propToDCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};

  OutputObj<TH1F> hMass{TH1F("hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", 500, 0., 5.)};

  double massPiK{0.};
  double massKPi{0.};

  void process(aod::Collisions const& collisions,
               aod::HfTrackIndexProng2 const& rowsTrackIndexProng2,
               aod::BigTracks const& tracks)
  {
    // 2-prong vertex fitter
    o2::vertexing::DCAFitterN<2> df;
    df.setBz(magneticField);
    df.setPropagateToPCA(propToDCA);
    df.setMaxR(maxR);
    df.setMaxDZIni(maxDZIni);
    df.setMinParamChange(minParamChange);
    df.setMinRelChi2Change(minRelChi2Change);
    df.setUseAbsDCA(useAbsDCA);

    // loop over pairs of track indices
    for (const auto& rowTrackIndexProng2 : rowsTrackIndexProng2) {
      auto track0 = rowTrackIndexProng2.index0_as<aod::BigTracks>();
      auto track1 = rowTrackIndexProng2.index1_as<aod::BigTracks>();
      auto trackParVarPos1 = getTrackParCov(track0);
      auto trackParVarNeg1 = getTrackParCov(track1);
      auto collision = track0.collision();

      // reconstruct the 2-prong secondary vertex
      if (df.process(trackParVarPos1, trackParVarNeg1) == 0) {
        continue;
      }
      const auto& secondaryVertex = df.getPCACandidate();
      auto trackParVar0 = df.getTrack(0);
      auto trackParVar1 = df.getTrack(1);

      // get track momenta
      array<float, 3> pvec0;
      array<float, 3> pvec1;
      trackParVar0.getPxPyPzGlo(pvec0);
      trackParVar1.getPxPyPzGlo(pvec1);

      // fill candidate table rows
      rowCandidateBase(collision.globalIndex(),
                       collision.posX(), collision.posY(), collision.posZ(),
                       secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                       pvec0[0], pvec0[1], pvec0[2],
                       pvec1[0], pvec1[1], pvec1[2],
                       rowTrackIndexProng2.index0Id(), rowTrackIndexProng2.index1Id());

      // fill histograms
      // calculate invariant masses
      auto arrayMomenta = std::array{pvec0, pvec1};
      massPiK = RecoDecay::M(arrayMomenta, arrMassPiK);
      massKPi = RecoDecay::M(arrayMomenta, arrMassKPi);
      hMass->Fill(massPiK);
      // hMass->Fill(massKPi);
    }
  }
};

/// Helper extension task
/// Extends the base table with expression columns.
struct HfCandidateCreator2ProngExpressions {
  Spawns<aod::HfCandProng2Ext> rowCandidateProng2;
  void init(InitContext const&) {}
};

// Candidate selection =====================================================================

namespace o2::aod
{
namespace hf_selcandidate_d0
{
// Candidate selection columns
DECLARE_SOA_COLUMN(IsSelD0, isSelD0, int);       //!
DECLARE_SOA_COLUMN(IsSelD0bar, isSelD0bar, int); //!
} // namespace hf_selcandidate_d0

// Candidate selection table
DECLARE_SOA_TABLE(HfSelCandidateD0, "AOD", "HFSELCANDD0", //!
                  hf_selcandidate_d0::IsSelD0,
                  hf_selcandidate_d0::IsSelD0bar);
} // namespace o2::aod

/// D0 candidate selector
struct HfCandidateSelectorD0 {
  Produces<aod::HfSelCandidateD0> hfSelD0Candidate;

  Configurable<double> pTCandMin{"pTCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> pTCandMax{"pTCandMax", 50., "Upper bound of candidate pT"};
  // TPC
  Configurable<double> pTPidTpcMin{"pTPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> pTPidTpcMax{"pTPidTpcMax", 5., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpc{"nSigmaTpc", 3., "Nsigma cut on TPC only"};
  // topological cuts
  Configurable<double> cpaMin{"cpaMin", 0.98, "Min. cosine of pointing angle"};
  Configurable<double> massWindow{"massWindow", 0.4, "Half-width of the invariant-mass window"};

  /// Conjugate-independent topological cuts
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T>
  bool selectionTopol(const T& candidate)
  {
    // check that the candidate pT is within the analysis range
    if (candidate.pt() < pTCandMin || candidate.pt() >= pTCandMax) {
      return false;
    }
    // cosine of pointing angle
    if (candidate.cpa() < cpaMin) {
      return false;
    }
    return true;
  }

  /// Conjugate-dependent topological cuts
  /// \param candidate is candidate
  /// \param trackPion is the track with the pion hypothesis
  /// \param trackKaon is the track with the kaon hypothesis
  /// \note trackPion = positive and trackKaon = negative for D0 selection and inverse for D0bar
  /// \return true if candidate passes all cuts for the given Conjugate
  template <typename T1, typename T2>
  bool selectionTopolConjugate(const T1& candidate, const T2& trackPion, const T2& trackKaon)
  {
    // invariant-mass cut
    if (trackPion.sign() > 0) {
      if (std::abs(InvMassD0(candidate) - RecoDecay::getMassPDG(pdg::Code::kD0)) > massWindow) {
        return false;
      }
    } else {
      if (std::abs(InvMassD0bar(candidate) - RecoDecay::getMassPDG(pdg::Code::kD0)) > massWindow) {
        return false;
      }
    }
    return true;
  }

  void process(aod::HfCandProng2 const& candidates, aod::BigTracksPIDExtended const&)
  {
    TrackSelectorPID selectorPion(kPiPlus);
    selectorPion.setRangePtTPC(pTPidTpcMin, pTPidTpcMax);
    selectorPion.setRangeNSigmaTPC(-nSigmaTpc, nSigmaTpc);

    TrackSelectorPID selectorKaon(selectorPion);
    selectorKaon.setPDG(kKPlus);

    // looping over 2-prong candidates
    for (auto& candidate : candidates) {

      // final selection flag: 0 - rejected, 1 - accepted
      int statusD0 = 0;
      int statusD0bar = 0;

      auto trackPos = candidate.index0_as<aod::BigTracksPIDExtended>(); // positive daughter
      auto trackNeg = candidate.index1_as<aod::BigTracksPIDExtended>(); // negative daughter

      // conjugate-independent topological selection
      if (!selectionTopol(candidate)) {
        hfSelD0Candidate(statusD0, statusD0bar);
        continue;
      }

      // conjugate-dependent topological selection for D0
      bool topolD0 = selectionTopolConjugate(candidate, trackPos, trackNeg);
      // conjugate-dependent topological selection for D0bar
      bool topolD0bar = selectionTopolConjugate(candidate, trackNeg, trackPos);

      if (!topolD0 && !topolD0bar) {
        hfSelD0Candidate(statusD0, statusD0bar);
        continue;
      }

      // track-level PID selection
      int pidTrackPosKaon = selectorKaon.getStatusTrackPIDAll(trackPos);
      int pidTrackPosPion = selectorPion.getStatusTrackPIDAll(trackPos);
      int pidTrackNegKaon = selectorKaon.getStatusTrackPIDAll(trackNeg);
      int pidTrackNegPion = selectorPion.getStatusTrackPIDAll(trackNeg);

      int pidD0 = -1;
      int pidD0bar = -1;

      if (pidTrackPosPion == TrackSelectorPID::Status::PIDAccepted &&
          pidTrackNegKaon == TrackSelectorPID::Status::PIDAccepted) {
        pidD0 = 1; // accept D0
      } else if (pidTrackPosPion == TrackSelectorPID::Status::PIDRejected ||
                 pidTrackNegKaon == TrackSelectorPID::Status::PIDRejected ||
                 pidTrackNegPion == TrackSelectorPID::Status::PIDAccepted ||
                 pidTrackPosKaon == TrackSelectorPID::Status::PIDAccepted) {
        pidD0 = 0; // exclude D0
      }

      if (pidTrackNegPion == TrackSelectorPID::Status::PIDAccepted &&
          pidTrackPosKaon == TrackSelectorPID::Status::PIDAccepted) {
        pidD0bar = 1; // accept D0bar
      } else if (pidTrackPosPion == TrackSelectorPID::Status::PIDAccepted ||
                 pidTrackNegKaon == TrackSelectorPID::Status::PIDAccepted ||
                 pidTrackNegPion == TrackSelectorPID::Status::PIDRejected ||
                 pidTrackPosKaon == TrackSelectorPID::Status::PIDRejected) {
        pidD0bar = 0; // exclude D0bar
      }

      if (pidD0 == 0 && pidD0bar == 0) {
        hfSelD0Candidate(statusD0, statusD0bar);
        continue;
      }

      if ((pidD0 == -1 || pidD0 == 1) && topolD0) {
        statusD0 = 1; // identified as D0
      }
      if ((pidD0bar == -1 || pidD0bar == 1) && topolD0bar) {
        statusD0bar = 1; // identified as D0bar
      }

      hfSelD0Candidate(statusD0, statusD0bar);
    }
  }
};

// Analysis task =====================================================================

/// D0 analysis task
struct HfTaskD0 {
  HistogramRegistry registry{
    "registry",
    {}};

  Configurable<int> flagSelCandD0{"flagSelCandD0", 1, "Selection flag for D0"};
  Configurable<int> flagSelCandD0bar{"flagSelCandD0bar", 1, "Selection flag for D0 bar"};

  void init(o2::framework::InitContext&)
  {
    const TString strTitle = "D^{0} candidates";
    const TString strPT = "#it{p}_{T} (GeV/#it{c})";
    const TString strEntries = "entries";
    registry.add("hPTCand", strTitle + ";" + strPT + ";" + strEntries, {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hMass", strTitle + ";" + "inv. mass (#pi K) (GeV/#it{c}^{2})" + ";" + strEntries, {HistType::kTH1F, {{500, 0., 5.}}});
    registry.add("hCPA", strTitle + ";" + "cosine of pointing angle" + ";" + strPT + ";" + strEntries, {HistType::kTH2F, {{110, -1.1, 1.1}, {100, 0., 10.}}});
  }

  Partition<soa::Join<aod::HfCandProng2, aod::HfSelCandidateD0>> selectedD0Candidates = aod::hf_selcandidate_d0::isSelD0 >= flagSelCandD0 || aod::hf_selcandidate_d0::isSelD0bar >= flagSelCandD0bar;

  void process(soa::Join<aod::HfCandProng2, aod::HfSelCandidateD0>& candidates)
  {
    for (auto& candidate : selectedD0Candidates) {
      if (candidate.isSelD0() >= flagSelCandD0) {
        registry.fill(HIST("hMass"), InvMassD0(candidate));
      }
      if (candidate.isSelD0bar() >= flagSelCandD0bar) {
        registry.fill(HIST("hMass"), InvMassD0bar(candidate));
      }
      registry.fill(HIST("hPTCand"), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
    }
  }
};

// Add all tasks in the workflow specification.
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTagSelTracks>(cfgc),
    adaptAnalysisTask<HfTrackIndexSkimsCreator>(cfgc),
    adaptAnalysisTask<HfCandidateCreator2Prong>(cfgc),
    adaptAnalysisTask<HfCandidateCreator2ProngExpressions>(cfgc),
    adaptAnalysisTask<HfCandidateSelectorD0>(cfgc),
    adaptAnalysisTask<HfTaskD0>(cfgc)};
}
