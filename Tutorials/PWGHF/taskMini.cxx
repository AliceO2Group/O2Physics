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
/// \author Vít Kučera <vit.kucera@cern.ch>, Inha University

#include <algorithm>

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelectorPID.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "DCAFitter/DCAFitterN.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"

#include "PWGHF/Core/HfHelper.h"
#include "PWGHF/Core/SelectorCuts.h"

using namespace o2;
using namespace o2::analysis;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Track selection =====================================================================

namespace o2::aod
{
namespace hf_seltrack
{
// Track selection columns
DECLARE_SOA_COLUMN(IsSelProng, isSelProng, bool); //! prong selection flag
} // namespace hf_seltrack

// Track selection table
DECLARE_SOA_TABLE(HfSelTrack, "AOD", "HFSELTRACK", //! track selection table
                  hf_seltrack::IsSelProng);
} // namespace o2::aod

/// Track selection
struct HfTagSelTracks {
  Produces<aod::HfSelTrack> rowSelectedTrack;

  // 2-prong cuts
  Configurable<double> ptTrackMin{"ptTrackMin", -1., "min. track pT for 2 prong candidate"};
  Configurable<double> etaTrackMax{"etaTrackMax", 4., "max. pseudorapidity for 2 prong candidate"};
  Configurable<double> dcaTrackMin{"dcaTrackMin", 0.0025, "min. DCA for 2 prong candidate"};

  using TracksWithDca = soa::Join<aod::Tracks, aod::TracksDCA>;

  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    const TString strTitle = "D^{0} candidates";
    const TString strPt = "#it{p}_{T}^{track} (GeV/#it{c})";
    const TString strEntries = "entries";
    registry.add("hPtNoCuts", "all tracks;" + strPt + ";" + strEntries, {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hPtCuts2Prong", "tracks selected for 2-prong vertexing;" + strPt + ";" + strEntries, {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hPtVsDcaXYToPvCuts2Prong", "tracks selected for 2-prong vertexing;" + strPt + ";" + "DCAxy to prim. vtx. (cm)" + ";" + strEntries, {HistType::kTH2F, {{100, 0., 10.}, {400, -2., 2.}}});
    registry.add("hEtaCuts2Prong", "tracks selected for 2-prong vertexing;#it{#eta};" + strEntries, {HistType::kTH1F, {{static_cast<int>(1.2 * etaTrackMax * 100), -1.2 * etaTrackMax, 1.2 * etaTrackMax}}});
  }

  void process(TracksWithDca const& tracks)
  {
    for (const auto& track : tracks) {
      bool statusProng = true;

      auto ptTrack = track.pt();
      registry.fill(HIST("hPtNoCuts"), ptTrack);

      // pT cut
      if (ptTrack < ptTrackMin) {
        statusProng = false;
      }

      // eta cut
      auto etaTrack = track.eta();
      if (statusProng && std::abs(etaTrack) > etaTrackMax) {
        statusProng = false;
      }

      // DCA cut
      auto dcaXY = track.dcaXY();
      if (statusProng && std::abs(dcaXY) < dcaTrackMin) {
        statusProng = false;
      }

      // fill histograms
      if (statusProng) {
        registry.fill(HIST("hPtCuts2Prong"), ptTrack);
        registry.fill(HIST("hEtaCuts2Prong"), etaTrack);
        registry.fill(HIST("hPtVsDcaXYToPvCuts2Prong"), ptTrack, dcaXY);
      }

      // fill table row
      rowSelectedTrack(statusProng);
    }
  }
};

// Track index skimming =====================================================================

namespace o2::aod
{
namespace hf_track_index
{
// Track index skim columns
DECLARE_SOA_INDEX_COLUMN_FULL(Prong0, prong0, int, Tracks, "_0"); //! prong 0
DECLARE_SOA_INDEX_COLUMN_FULL(Prong1, prong1, int, Tracks, "_1"); //! prong 1
} // namespace hf_track_index

// Track index skim table
DECLARE_SOA_TABLE(HfTrackIndexProng2, "AOD", "HFTRACKIDXP2", //! table with prongs indices
                  hf_track_index::Prong0Id,
                  hf_track_index::Prong1Id);
} // namespace o2::aod

/// Track index skim creator
/// Pre-selection of 2-prong secondary vertices
struct HfTrackIndexSkimCreator {
  Produces<aod::HfTrackIndexProng2> rowTrackIndexProng2;

  // vertexing parameters
  Configurable<double> magneticField{"magneticField", 5., "magnetic field [kG]"};
  Configurable<bool> propToDCA{"propToDCA", true, "create tracks version propagated to PCA"};
  Configurable<bool> useAbsDCA{"useAbsDCA", true, "Minimise abs. distance rather than chi2"};
  Configurable<double> maxR{"maxR", 200., "reject PCA's above this radius"};
  Configurable<double> maxDZIni{"maxDZIni", 4., "reject (if>0) PCA candidate if tracks DZ exceeds threshold"};
  Configurable<double> minParamChange{"minParamChange", 1.e-3, "stop iterations if largest change of any X is smaller than this"};
  Configurable<double> minRelChi2Change{"minRelChi2Change", 0.9, "stop iterations if chi2/chi2old > this"};

  HfHelper hfHelper;

  using SelectedTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksCov, aod::HfSelTrack>>;

  Filter filterSelectTracks = aod::hf_seltrack::isSelProng == true;

  HistogramRegistry registry{
    "registry",
    {// 2-prong histograms
     {"hVtx2ProngX", "2-prong candidates;#it{x}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}}},
     {"hVtx2ProngY", "2-prong candidates;#it{y}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -2., 2.}}}},
     {"hVtx2ProngZ", "2-prong candidates;#it{z}_{sec. vtx.} (cm);entries", {HistType::kTH1F, {{1000, -20., 20.}}}},
     {"hMassD0ToPiK", "D^{0} candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", {HistType::kTH1F, {{500, 0., 5.}}}}}};

  void init(InitContext&)
  {
  }

  void process(aod::Collision const&,
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

    // loop over positive tracks
    for (const auto& trackPos1 : tracks) {
      if (trackPos1.signed1Pt() < 0) {
        continue;
      }
      auto trackParVarPos1 = getTrackParCov(trackPos1);

      // loop over negative tracks
      for (const auto& trackNeg1 : tracks) {
        if (trackNeg1.signed1Pt() > 0) {
          continue;
        }
        auto trackParVarNeg1 = getTrackParCov(trackNeg1);

        // secondary vertex reconstruction and further 2-prong selections
        if (df2.process(trackParVarPos1, trackParVarNeg1) == 0) {
          continue;
        }
        //  get secondary vertex
        const auto& secondaryVertex = df2.getPCACandidate();
        // get track momenta
        std::array<float, 3> pVec0;
        std::array<float, 3> pVec1;
        df2.getTrack(0).getPxPyPzGlo(pVec0);
        df2.getTrack(1).getPxPyPzGlo(pVec1);

        // fill table row
        rowTrackIndexProng2(trackPos1.globalIndex(),
                            trackNeg1.globalIndex());

        // fill histograms
        registry.fill(HIST("hVtx2ProngX"), secondaryVertex[0]);
        registry.fill(HIST("hVtx2ProngY"), secondaryVertex[1]);
        registry.fill(HIST("hVtx2ProngZ"), secondaryVertex[2]);
        std::array<std::array<float, 3>, 2> arrMom = {pVec0, pVec1};
        auto mass2Prong = RecoDecay::m(arrMom, std::array{o2::analysis::pdg::MassPiPlus, o2::analysis::pdg::MassKPlus});
        registry.fill(HIST("hMassD0ToPiK"), mass2Prong);
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
DECLARE_SOA_INDEX_COLUMN(Collision, collision); //! collisions
// secondary vertex
DECLARE_SOA_COLUMN(XSecondaryVertex, xSecondaryVertex, float); //! x coordinate of the secondary vertex
DECLARE_SOA_COLUMN(YSecondaryVertex, ySecondaryVertex, float); //! y coordinate of the secondary vertex
DECLARE_SOA_COLUMN(ZSecondaryVertex, zSecondaryVertex, float); //! z coordinate of the secondary vertex
DECLARE_SOA_DYNAMIC_COLUMN(RSecondaryVertex, rSecondaryVertex, //! radius of the secondary vertex
                           [](float xVtxS, float yVtxS) -> float { return RecoDecay::sqrtSumOfSquares(xVtxS, yVtxS); });
// prong properties
DECLARE_SOA_COLUMN(PxProng0, pxProng0, float); //! px of prong 0
DECLARE_SOA_COLUMN(PyProng0, pyProng0, float); //! py of prong 0
DECLARE_SOA_COLUMN(PzProng0, pzProng0, float); //! pz of prong 0
DECLARE_SOA_DYNAMIC_COLUMN(PtProng0, ptProng0, //! pt of prong 0
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_COLUMN(PxProng1, pxProng1, float); //! px of prong 1
DECLARE_SOA_COLUMN(PyProng1, pyProng1, float); //! py of prong 1
DECLARE_SOA_COLUMN(PzProng1, pzProng1, float); //! pz of prong 1
DECLARE_SOA_DYNAMIC_COLUMN(PtProng1, ptProng1, //! pt of prong 1
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
// candidate properties
DECLARE_SOA_DYNAMIC_COLUMN(DecayLength, decayLength, //! decay length of candidate
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS) -> float { return RecoDecay::distance(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}); });
DECLARE_SOA_DYNAMIC_COLUMN(Pt, pt, //! pt of candidate
                           [](float px, float py) -> float { return RecoDecay::pt(px, py); });
DECLARE_SOA_EXPRESSION_COLUMN(Px, px, //! px of candidate
                              float, 1.f * pxProng0 + 1.f * pxProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Py, py, //! py of candidate
                              float, 1.f * pyProng0 + 1.f * pyProng1);
DECLARE_SOA_EXPRESSION_COLUMN(Pz, pz, //! pz of candidate
                              float, 1.f * pzProng0 + 1.f * pzProng1);
DECLARE_SOA_DYNAMIC_COLUMN(M, m, //! invariant mass of candidate
                           [](float px0, float py0, float pz0, float px1, float py1, float pz1, const std::array<double, 2>& m) -> float { return RecoDecay::m(std::array{std::array{px0, py0, pz0}, std::array{px1, py1, pz1}}, m); });
DECLARE_SOA_DYNAMIC_COLUMN(CPA, cpa, //! cosine of pointing angle of candidate
                           [](float xVtxP, float yVtxP, float zVtxP, float xVtxS, float yVtxS, float zVtxS, float px, float py, float pz) -> float { return RecoDecay::cpa(std::array{xVtxP, yVtxP, zVtxP}, std::array{xVtxS, yVtxS, zVtxS}, std::array{px, py, pz}); });
} // namespace hf_cand_prong2

// Candidate table
DECLARE_SOA_TABLE(HfCandProng2Base, "AOD", "HFCANDP2BASE", //! 2-prong candidate table
                  hf_cand_prong2::CollisionId,
                  collision::PosX, collision::PosY, collision::PosZ,
                  hf_cand_prong2::XSecondaryVertex, hf_cand_prong2::YSecondaryVertex, hf_cand_prong2::ZSecondaryVertex,
                  /* dynamic columns */ hf_cand_prong2::RSecondaryVertex<hf_cand_prong2::XSecondaryVertex, hf_cand_prong2::YSecondaryVertex>,
                  hf_cand_prong2::DecayLength<collision::PosX, collision::PosY, collision::PosZ, hf_cand_prong2::XSecondaryVertex, hf_cand_prong2::YSecondaryVertex, hf_cand_prong2::ZSecondaryVertex>,
                  /* prong 0 */ hf_cand_prong2::PtProng0<hf_cand_prong2::PxProng0, hf_cand_prong2::PyProng0>,
                  hf_cand_prong2::PxProng0, hf_cand_prong2::PyProng0, hf_cand_prong2::PzProng0,
                  /* prong 1 */ hf_cand_prong2::PtProng1<hf_cand_prong2::PxProng1, hf_cand_prong2::PyProng1>,
                  hf_cand_prong2::PxProng1, hf_cand_prong2::PyProng1, hf_cand_prong2::PzProng1,
                  hf_track_index::Prong0Id, hf_track_index::Prong1Id,
                  /* dynamic columns */
                  hf_cand_prong2::M<hf_cand_prong2::PxProng0, hf_cand_prong2::PyProng0, hf_cand_prong2::PzProng0, hf_cand_prong2::PxProng1, hf_cand_prong2::PyProng1, hf_cand_prong2::PzProng1>,
                  /* dynamic columns that use candidate momentum components */
                  hf_cand_prong2::CPA<collision::PosX, collision::PosY, collision::PosZ, hf_cand_prong2::XSecondaryVertex, hf_cand_prong2::YSecondaryVertex, hf_cand_prong2::ZSecondaryVertex, hf_cand_prong2::Px, hf_cand_prong2::Py, hf_cand_prong2::Pz>,
                  hf_cand_prong2::Pt<hf_cand_prong2::Px, hf_cand_prong2::Py>);

// Extended table with expression columns that can be used as arguments of dynamic columns
DECLARE_SOA_EXTENDED_TABLE_USER(HfCandProng2Ext, HfCandProng2Base, "HFCANDP2EXT", //! extension table for the 2-prong candidate table
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

  HfHelper hfHelper;

  double massPiK{0.};
  double massKPi{0.};

  using TracksWithCov = soa::Join<Tracks, TracksCov>;

  OutputObj<TH1F> hMass{TH1F("hMass", "2-prong candidates;inv. mass (#pi K) (GeV/#it{c}^{2});entries", 500, 0., 5.)};

  void process(aod::Collisions const&,
               aod::HfTrackIndexProng2 const& rowsTrackIndexProng2,
               TracksWithCov const&)
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
      auto track0 = rowTrackIndexProng2.prong0_as<TracksWithCov>();
      auto track1 = rowTrackIndexProng2.prong1_as<TracksWithCov>();
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
      std::array<float, 3> pVec0;
      std::array<float, 3> pVec1;
      trackParVar0.getPxPyPzGlo(pVec0);
      trackParVar1.getPxPyPzGlo(pVec1);

      // fill candidate table rows
      rowCandidateBase(collision.globalIndex(),
                       collision.posX(), collision.posY(), collision.posZ(),
                       secondaryVertex[0], secondaryVertex[1], secondaryVertex[2],
                       pVec0[0], pVec0[1], pVec0[2],
                       pVec1[0], pVec1[1], pVec1[2],
                       rowTrackIndexProng2.prong0Id(), rowTrackIndexProng2.prong1Id());

      // fill histograms
      // calculate invariant masses
      auto arrayMomenta = std::array{pVec0, pVec1};
      massPiK = RecoDecay::m(arrayMomenta, std::array{o2::analysis::pdg::MassPiPlus, o2::analysis::pdg::MassKPlus});
      massKPi = RecoDecay::m(arrayMomenta, std::array{o2::analysis::pdg::MassKPlus, o2::analysis::pdg::MassPiPlus});
      hMass->Fill(massPiK);
      // hMass->Fill(massKPi);
    }
  }
};

/// Helper extension task
/// Extends the base table with expression columns (see the HfCandProng2Ext table).
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
DECLARE_SOA_COLUMN(IsSelD0, isSelD0, int);       //! selection flag for D0
DECLARE_SOA_COLUMN(IsSelD0bar, isSelD0bar, int); //! selection flag for D0 bar
} // namespace hf_selcandidate_d0

// Candidate selection table
DECLARE_SOA_TABLE(HfSelCandidateD0, "AOD", "HFSELCANDD0", //! table with D0 selection flags
                  hf_selcandidate_d0::IsSelD0,
                  hf_selcandidate_d0::IsSelD0bar);
} // namespace o2::aod

/// D0 candidate selector
struct HfCandidateSelectorD0 {
  Produces<aod::HfSelCandidateD0> hfSelD0Candidate;

  Configurable<double> ptCandMin{"ptCandMin", 0., "Lower bound of candidate pT"};
  Configurable<double> ptCandMax{"ptCandMax", 50., "Upper bound of candidate pT"};
  // TPC
  Configurable<double> ptPidTpcMin{"ptPidTpcMin", 0.15, "Lower bound of track pT for TPC PID"};
  Configurable<double> ptPidTpcMax{"ptPidTpcMax", 5., "Upper bound of track pT for TPC PID"};
  Configurable<double> nSigmaTpc{"nSigmaTpc", 3., "Nsigma cut on TPC only"};
  // topological cuts
  Configurable<double> cpaMin{"cpaMin", 0.98, "Min. cosine of pointing angle"};
  Configurable<double> massWindow{"massWindow", 0.4, "Half-width of the invariant-mass window"};

  HfHelper hfHelper;
  TrackSelectorPi selectorPion;
  TrackSelectorKa selectorKaon;

  using TracksWithPid = soa::Join<Tracks,
                                  aod::pidTPCFullPi, aod::pidTPCFullKa,
                                  aod::pidTOFFullPi, aod::pidTOFFullKa>;

  void init(InitContext const&)
  {
    selectorPion.setRangePtTpc(ptPidTpcMin, ptPidTpcMax);
    selectorPion.setRangeNSigmaTpc(-nSigmaTpc, nSigmaTpc);
    selectorKaon = selectorPion;
  }

  /// Conjugate-independent topological cuts
  /// \param candidate is candidate
  /// \return true if candidate passes all cuts
  template <typename T>
  bool selectionTopol(const T& candidate)
  {
    // check that the candidate pT is within the analysis range
    if (candidate.pt() < ptCandMin || candidate.pt() >= ptCandMax) {
      return false;
    }
    // cosine of pointing angle
    if (candidate.cpa() < cpaMin) {
      return false;
    }
    return true;
  }

  /// Conjugate-dependent topological cuts
  /// \param candidate candidate
  /// \param trackPion the track with the pion hypothesis
  /// \param trackKaon the track with the kaon hypothesis
  /// \note trackPion = positive and trackKaon = negative for D0 selection and inverse for D0bar
  /// \return true if candidate passes all cuts for the given conjugate
  template <typename T1, typename T2>
  bool selectionTopolConjugate(const T1& candidate, const T2& trackPion, const T2& trackKaon)
  {
    // invariant-mass cut
    if (trackPion.sign() > 0) {
      if (std::abs(hfHelper.invMassD0ToPiK(candidate) - o2::analysis::pdg::MassD0) > massWindow) {
        return false;
      }
    } else {
      if (std::abs(hfHelper.invMassD0barToKPi(candidate) - o2::analysis::pdg::MassD0) > massWindow) {
        return false;
      }
    }
    return true;
  }

  void process(aod::HfCandProng2 const& candidates,
               TracksWithPid const&)
  {
    // looping over 2-prong candidates
    for (const auto& candidate : candidates) {

      // final selection flag: 0 - rejected, 1 - accepted
      int statusD0 = 0;
      int statusD0bar = 0;

      auto trackPos = candidate.prong0_as<TracksWithPid>(); // positive daughter
      auto trackNeg = candidate.prong1_as<TracksWithPid>(); // negative daughter

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
      int pidTrackPosKaon = selectorKaon.statusTpcOrTof(trackPos);
      int pidTrackPosPion = selectorPion.statusTpcOrTof(trackPos);
      int pidTrackNegKaon = selectorKaon.statusTpcOrTof(trackNeg);
      int pidTrackNegPion = selectorPion.statusTpcOrTof(trackNeg);

      int pidD0 = -1;
      int pidD0bar = -1;

      if (pidTrackPosPion == TrackSelectorPID::Accepted &&
          pidTrackNegKaon == TrackSelectorPID::Accepted) {
        pidD0 = 1; // accept D0
      } else if (pidTrackPosPion == TrackSelectorPID::Rejected ||
                 pidTrackNegKaon == TrackSelectorPID::Rejected ||
                 pidTrackNegPion == TrackSelectorPID::Accepted ||
                 pidTrackPosKaon == TrackSelectorPID::Accepted) {
        pidD0 = 0; // exclude D0
      }

      if (pidTrackNegPion == TrackSelectorPID::Accepted &&
          pidTrackPosKaon == TrackSelectorPID::Accepted) {
        pidD0bar = 1; // accept D0bar
      } else if (pidTrackPosPion == TrackSelectorPID::Accepted ||
                 pidTrackNegKaon == TrackSelectorPID::Accepted ||
                 pidTrackNegPion == TrackSelectorPID::Rejected ||
                 pidTrackPosKaon == TrackSelectorPID::Rejected) {
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
  Configurable<int> selectionFlagD0{"selectionFlagD0", 1, "Selection flag for D0"};
  Configurable<int> selectionFlagD0bar{"selectionFlagD0bar", 1, "Selection flag for D0 bar"};

  HfHelper hfHelper;

  Partition<soa::Join<aod::HfCandProng2, aod::HfSelCandidateD0>> selectedD0Candidates = aod::hf_selcandidate_d0::isSelD0 >= selectionFlagD0 || aod::hf_selcandidate_d0::isSelD0bar >= selectionFlagD0bar;

  HistogramRegistry registry{
    "registry",
    {}};

  void init(InitContext&)
  {
    const TString strTitle = "D^{0} candidates";
    const TString strPt = "#it{p}_{T} (GeV/#it{c})";
    const TString strEntries = "entries";
    registry.add("hPtCand", strTitle + ";" + strPt + ";" + strEntries, {HistType::kTH1F, {{100, 0., 10.}}});
    registry.add("hMass", strTitle + ";" + "inv. mass (#pi K) (GeV/#it{c}^{2})" + ";" + strEntries, {HistType::kTH1F, {{500, 0., 5.}}});
    registry.add("hCpaVsPtCand", strTitle + ";" + "cosine of pointing angle" + ";" + strPt + ";" + strEntries, {HistType::kTH2F, {{110, -1.1, 1.1}, {100, 0., 10.}}});
  }

  void process(soa::Join<aod::HfCandProng2, aod::HfSelCandidateD0> const& candidates)
  {
    for (const auto& candidate : selectedD0Candidates) {
      if (candidate.isSelD0() >= selectionFlagD0) {
        registry.fill(HIST("hMass"), hfHelper.invMassD0ToPiK(candidate));
      }
      if (candidate.isSelD0bar() >= selectionFlagD0bar) {
        registry.fill(HIST("hMass"), hfHelper.invMassD0barToKPi(candidate));
      }
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hCpaVsPtCand"), candidate.cpa(), candidate.pt());
    }
  }
};

// Add all tasks in the workflow specification.
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HfTagSelTracks>(cfgc),
    adaptAnalysisTask<HfTrackIndexSkimCreator>(cfgc),
    adaptAnalysisTask<HfCandidateCreator2Prong>(cfgc),
    adaptAnalysisTask<HfCandidateCreator2ProngExpressions>(cfgc),
    adaptAnalysisTask<HfCandidateSelectorD0>(cfgc),
    adaptAnalysisTask<HfTaskD0>(cfgc)};
}
