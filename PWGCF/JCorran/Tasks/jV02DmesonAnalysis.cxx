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

/// \file jV02DmesonAnalysis.cxx
/// \brief D0 v02(pT) on the JCorran derived tables
/// \author Dong Jo Kim (djkim@jyu.fi)
///
/// v02(pT) = <X Y Z(pT)>_c / (<X Y>_c <Z(pT)>).
/// X and Y are the corrected Q2 of two forward detectors (default FT0A x FT0C)
/// from the combined Qvectors table, with the same indexing as jEPDzeroFlowAnalysis.
/// Z_D(pT) = N_D(pT) / N_ch. N_D comes from CF2ProngTracks (D0 and D0bar).
/// N_ch comes from CFTracks. Both derived tables are indexed to CFCollisions,
/// so the original collision is matched through CFCollRefs (one row per kept
/// collision, written by filter-correlations with cfgTransientTables).
/// The D azimuth is not part of the Q-vector.
/// Do not divide by an extra resolution: it cancels in the ratio.
///
/// Definition: arXiv:2506.18690, Appendix A.

#include "PWGCF/DataModel/CorrelationsDerived.h"

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Qvectors.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <TAxis.h>
#include <TProfile2D.h>

#include <cmath>
#include <cstdint>
#include <random>
#include <string>
#include <unordered_set>
#include <vector>

using namespace o2;
using namespace o2::framework;

#define V02_FILL_QA(NAME, STEP)                            \
  for (int qaBin = 1; qaBin <= (STEP); ++qaBin) {          \
    histos.fill(HIST(NAME), static_cast<float>(qaBin));   \
  }

struct JV02DmesonAnalysis {
  enum Detector {
    kFT0C = 0,
    kFT0A,
    kFT0M,
    kFV0A,
    kTPCpos,
    kTPCneg,
    kTPCall
  };
  static constexpr int NDetectors = 7;
  static constexpr int NCorrections = 4;
  static constexpr int FullyCorrected = 3;
  static constexpr int StepAccepted = 5;
  static constexpr float InvalidQvector = 900.f;

  using MyCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::CentFT0Cs, aod::Qvectors>;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Configurable<std::string> cfgCentEst{"cfgCentEst", "FT0C", "Centrality estimator: FT0C or FT0M"};
  Configurable<std::string> cfgRefAName{"cfgRefAName", "FT0A", "Reference detector A: FT0A, FT0C, FT0M, FV0A, TPCpos, TPCneg, TPCall"};
  Configurable<std::string> cfgRefBName{"cfgRefBName", "FT0C", "Reference detector B"};
  Configurable<int> cfgQvecMode{"cfgQvecMode", 0, "Harmonic slot in qvector-table cfgnMods. 0 is v2 when cfgnMods starts at 2"};
  Configurable<int> cfgNSubsample{"cfgNSubsample", 10, "Number of bootstrap subsamples"};
  Configurable<int> cfgSeed{"cfgSeed", 1, "Seed for the subsample index"};
  Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 999999, "Maximum track occupancy"};
  Configurable<int> cfgMinOccupancy{"cfgMinOccupancy", 0, "Minimum track occupancy"};
  Configurable<float> cfgVtxZ{"cfgVtxZ", 10.f, "Maximum |z_vtx| (cm)"};
  Configurable<float> cfgCentMin{"cfgCentMin", 0.f, "Minimum centrality"};
  Configurable<float> cfgCentMax{"cfgCentMax", 80.f, "Maximum centrality"};
  Configurable<float> cfgEtaMax{"cfgEtaMax", 0.8f, "Maximum |eta| of CF tracks used for N_ch"};
  Configurable<float> cfgPtMinCh{"cfgPtMinCh", 0.2f, "Minimum charged-track pT (GeV/c)"};
  Configurable<float> cfgPtMaxCh{"cfgPtMaxCh", 3.f, "Maximum charged-track pT (GeV/c)"};
  Configurable<float> cfgPtMinD{"cfgPtMinD", 1.f, "Minimum D0 pT (GeV/c)"};
  Configurable<float> cfgPtMaxD{"cfgPtMaxD", 24.f, "Maximum D0 pT (GeV/c)"};
  Configurable<float> cfgMinNch{"cfgMinNch", 10.f, "Minimum N_ch"};
  Configurable<float> cfgMinSumAmpl{"cfgMinSumAmpl", 1e-4f, "Minimum Q-vector amplitude"};
  Configurable<bool> cfgRequireSel8{"cfgRequireSel8", true, "Require sel8"};
  Configurable<bool> cfgCutNoSameBunchPileup{"cfgCutNoSameBunchPileup", true, "Require kNoSameBunchPileup"};
  Configurable<bool> cfgCutGoodZvtxFT0vsPV{"cfgCutGoodZvtxFT0vsPV", true, "Require kIsGoodZvtxFT0vsPV"};
  Configurable<bool> cfgCutNoTimeFrameBorder{"cfgCutNoTimeFrameBorder", true, "Require kNoTimeFrameBorder"};
  Configurable<bool> cfgCutNoITSROFrameBorder{"cfgCutNoITSROFrameBorder", true, "Require kNoITSROFrameBorder"};
  Configurable<bool> cfgRejectDaughters{"cfgRejectDaughters", false, "Drop D0 daughter CF tracks from N_ch"};

  ConfigurableAxis cfgAxisCent{"cfgAxisCent", {8, 0., 80.}, "Centrality axis"};
  ConfigurableAxis cfgAxisPtD{"cfgAxisPtD", {VARIABLE_WIDTH, 1.f, 2.f, 3.f, 4.f, 5.f, 6.f, 8.f, 12.f, 16.f, 24.f}, "D0 pT axis"};
  ConfigurableAxis cfgAxisPtHad{"cfgAxisPtHad", {VARIABLE_WIDTH, 0.2f, 0.4f, 0.6f, 0.8f, 1.f, 1.2f, 1.5f, 2.f, 2.5f, 3.f}, "Charged-hadron pT axis"};
  ConfigurableAxis cfgAxisMass{"cfgAxisMass", {70, 1.70, 2.05}, "D0 invariant-mass axis"};

  int mRefA = kFT0A;
  int mRefB = kFT0C;
  int mNSubsample = 1;
  Preslice<aod::CFTracks> mPerCFTrack = aod::cftrack::cfCollisionId;
  Preslice<aod::CF2ProngTracks> mPerCFCand = aod::cftrack::cfCollisionId;
  const TAxis* mHadPtAxis = nullptr;
  std::vector<float> mHadCounts;
  std::mt19937 mRng{1};
  std::uniform_int_distribution<int> mSampleDist;

  int detectorId(const std::string& name) const
  {
    if (name == "FT0C") {
      return kFT0C;
    }
    if (name == "FT0A") {
      return kFT0A;
    }
    if (name == "FT0M") {
      return kFT0M;
    }
    if (name == "FV0A") {
      return kFV0A;
    }
    if (name == "TPCpos") {
      return kTPCpos;
    }
    if (name == "TPCneg") {
      return kTPCneg;
    }
    if (name == "TPCall") {
      return kTPCall;
    }
    return -1;
  }

  void init(InitContext const&)
  {
    mRefA = detectorId(cfgRefAName.value);
    mRefB = detectorId(cfgRefBName.value);
    if (mRefA < 0 || mRefB < 0 || mRefA == mRefB) {
      LOGF(fatal, "cfgRefAName and cfgRefBName must be two different detectors");
    }
    mNSubsample = cfgNSubsample > 0 ? cfgNSubsample.value : 1;
    mRng.seed(static_cast<std::mt19937::result_type>(cfgSeed));
    mSampleDist = std::uniform_int_distribution<int>(0, mNSubsample - 1);

    const AxisSpec subsampleAxis{mNSubsample, 0., static_cast<double>(mNSubsample), "subsample"};
    const AxisSpec qaAxis{6, 0.5, 6.5, "1 seen; 2 sel8; 3 vertex; 4 quality; 5 centrality; 6 accepted"};

    if (doprocessReference) {
      histos.add("hadron/hEventQa", "1 seen; 2 sel8; 3 vertex; 4 quality; 5 centrality; 6 accepted", {HistType::kTH1F, {qaAxis}});
      histos.add("hadron/hEvents", "Selected events;centrality (%);events", {HistType::kTH1F, {cfgAxisCent}});
      histos.add("hadron/Prof_XY", "Bulk #LTXY#GT;centrality (%);#LTXY#GT", {HistType::kTProfile, {cfgAxisCent}});
      histos.add("hadron/Prof_Nch", "N_{ch};centrality (%);#LTN_{ch}#GT", {HistType::kTProfile, {cfgAxisCent}});
      histos.add("hadron/Prof_Z", "Z = n(p_{T})/N_{ch};centrality (%);p_{T} (GeV/c);#LTZ#GT", {HistType::kTProfile2D, {cfgAxisCent, cfgAxisPtHad}});
      histos.add("hadron/Prof_XYZ", "XY Z;centrality (%);p_{T} (GeV/c);#LTXYZ#GT", {HistType::kTProfile2D, {cfgAxisCent, cfgAxisPtHad}});
      histos.add("hadron/Prof_XY_sub", "Bulk #LTXY#GT per subsample;subsample;centrality (%);#LTXY#GT", {HistType::kTProfile2D, {subsampleAxis, cfgAxisCent}});
      histos.add("hadron/Prof_Z_sub", "Z per subsample;subsample;centrality (%);p_{T} (GeV/c);#LTZ#GT", {HistType::kTProfile3D, {subsampleAxis, cfgAxisCent, cfgAxisPtHad}});
      histos.add("hadron/Prof_XYZ_sub", "XYZ per subsample;subsample;centrality (%);p_{T} (GeV/c);#LTXYZ#GT", {HistType::kTProfile3D, {subsampleAxis, cfgAxisCent, cfgAxisPtHad}});
      mHadPtAxis = histos.get<TProfile2D>(HIST("hadron/Prof_Z"))->GetYaxis();
      mHadCounts.assign(mHadPtAxis->GetNbins(), 0.f);
    }
    if (doprocessD0) {
      histos.add("D0/hEventQa", "1 seen; 2 sel8; 3 vertex; 4 quality; 5 centrality; 6 accepted", {HistType::kTH1F, {qaAxis}});
      histos.add("D0/hEvents", "Selected events;centrality (%);events", {HistType::kTH1F, {cfgAxisCent}});
      histos.add("D0/Prof_XY", "Bulk #LTXY#GT;centrality (%);#LTXY#GT", {HistType::kTProfile, {cfgAxisCent}});
      histos.add("D0/Prof_Nch", "N_{ch};centrality (%);#LTN_{ch}#GT", {HistType::kTProfile, {cfgAxisCent}});
      histos.add("D0/hMassRaw", "Weight 1/N_{ch};centrality (%);p_{T} (GeV/c);m (GeV/c^{2})", {HistType::kTH3F, {cfgAxisCent, cfgAxisPtD, cfgAxisMass}}, true);
      histos.add("D0/hMassCorr", "Weight XY/N_{ch};centrality (%);p_{T} (GeV/c);m (GeV/c^{2})", {HistType::kTH3F, {cfgAxisCent, cfgAxisPtD, cfgAxisMass}}, true);
      histos.add("D0/hMassRawSub", "Weight 1/N_{ch};subsample;centrality (%);p_{T} (GeV/c);m (GeV/c^{2})", {HistType::kTHnSparseF, {subsampleAxis, cfgAxisCent, cfgAxisPtD, cfgAxisMass}}, true);
      histos.add("D0/hMassCorrSub", "Weight XY/N_{ch};subsample;centrality (%);p_{T} (GeV/c);m (GeV/c^{2})", {HistType::kTHnSparseF, {subsampleAxis, cfgAxisCent, cfgAxisPtD, cfgAxisMass}}, true);
    }
  }

  template <typename TCollision>
  float centralityOf(const TCollision& collision)
  {
    if (cfgCentEst.value == "FT0M") {
      return collision.centFT0M();
    }
    return collision.centFT0C();
  }

  template <typename TCollision>
  int passedEventStep(const TCollision& collision, float& centrality)
  {
    if (cfgRequireSel8 && !collision.sel8()) {
      return 1;
    }
    if (std::abs(collision.posZ()) > cfgVtxZ) {
      return 2;
    }
    const bool rejected =
      (cfgCutNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) ||
      (cfgCutGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) ||
      (cfgCutNoTimeFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) ||
      (cfgCutNoITSROFrameBorder && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) ||
      collision.trackOccupancyInTimeRange() > cfgMaxOccupancy ||
      collision.trackOccupancyInTimeRange() < cfgMinOccupancy;
    if (rejected) {
      return 3;
    }
    centrality = centralityOf(collision);
    if (centrality < cfgCentMin || centrality > cfgCentMax) {
      return 4;
    }
    return 5;
  }

  template <typename TCollision>
  bool referenceProduct(const TCollision& collision, float& xy)
  {
    const auto& re = collision.qvecRe();
    const auto& im = collision.qvecIm();
    const auto& amp = collision.qvecAmp();
    const int idxA = cfgQvecMode * NDetectors * NCorrections + mRefA * NCorrections + FullyCorrected;
    const int idxB = cfgQvecMode * NDetectors * NCorrections + mRefB * NCorrections + FullyCorrected;
    const int ampA = cfgQvecMode * NDetectors + mRefA;
    const int ampB = cfgQvecMode * NDetectors + mRefB;
    if (idxA >= static_cast<int>(re.size()) || idxB >= static_cast<int>(im.size()) || ampA >= static_cast<int>(amp.size()) || ampB >= static_cast<int>(amp.size())) {
      return false;
    }
    if (!(amp[ampA] > cfgMinSumAmpl) || !(amp[ampB] > cfgMinSumAmpl)) {
      return false;
    }
    if (re[idxA] > InvalidQvector || re[idxB] > InvalidQvector) {
      return false;
    }
    xy = re[idxA] * re[idxB] + im[idxA] * im[idxB];
    return std::isfinite(xy);
  }

  template <typename TTracks>
  float countCharged(const TTracks& tracks, const std::unordered_set<int64_t>& daughters, bool differential)
  {
    float nch = 0.f;
    if (differential) {
      std::fill(mHadCounts.begin(), mHadCounts.end(), 0.f);
    }
    for (const auto& track : tracks) {
      const float pt = track.pt();
      if (track.sign() == 0 || pt < cfgPtMinCh || pt >= cfgPtMaxCh || std::abs(track.eta()) > cfgEtaMax) {
        continue;
      }
      if (cfgRejectDaughters && daughters.contains(static_cast<int64_t>(track.globalIndex()))) {
        continue;
      }
      nch += 1.f;
      if (differential && mHadPtAxis != nullptr) {
        const int bin = mHadPtAxis->FindBin(pt);
        if (bin >= 1 && bin <= mHadPtAxis->GetNbins()) {
          mHadCounts[bin - 1] += 1.f;
        }
      }
    }
    return nch;
  }

  static bool isD0Decay(uint8_t decay)
  {
    return decay == aod::cf2prongtrack::D0ToPiK || decay == aod::cf2prongtrack::D0barToKPi || decay == aod::cf2prongtrack::D0barToKPiExclusive;
  }

  template <typename TRefs>
  static int cfCollisionIndex(const TRefs& refs)
  {
    if (refs.size() != 1) {
      return -1;
    }
    return static_cast<int>(refs.begin().globalIndex());
  }

  void processReference(MyCollisions::iterator const& collision, aod::BCs const&, aod::CFCollisions const&, aod::CFCollRefs const& refs, aod::CFTracks const& tracks)
  {
    float cent = -1.f;
    float xy = 0.f;
    const int step = passedEventStep(collision, cent);
    V02_FILL_QA("hadron/hEventQa", step);
    const int cfId = cfCollisionIndex(refs);
    if (step < StepAccepted || cfId < 0 || !referenceProduct(collision, xy)) {
      return;
    }
    const std::unordered_set<int64_t> noDaughters;
    const float nch = countCharged(tracks.sliceBy(mPerCFTrack, cfId), noDaughters, true);
    if (nch < cfgMinNch) {
      return;
    }
    histos.fill(HIST("hadron/hEventQa"), 6.f);
    const int sample = mSampleDist(mRng);
    histos.fill(HIST("hadron/hEvents"), cent);
    histos.fill(HIST("hadron/Prof_XY"), cent, xy);
    histos.fill(HIST("hadron/Prof_Nch"), cent, nch);
    histos.fill(HIST("hadron/Prof_XY_sub"), static_cast<float>(sample), cent, xy);
    const int nBins = mHadPtAxis->GetNbins();
    for (int ibin = 1; ibin <= nBins; ++ibin) {
      const float pt = static_cast<float>(mHadPtAxis->GetBinCenter(ibin));
      const float z = mHadCounts[ibin - 1] / nch;
      histos.fill(HIST("hadron/Prof_Z"), cent, pt, z);
      histos.fill(HIST("hadron/Prof_XYZ"), cent, pt, xy * z);
      histos.fill(HIST("hadron/Prof_Z_sub"), static_cast<float>(sample), cent, pt, z);
      histos.fill(HIST("hadron/Prof_XYZ_sub"), static_cast<float>(sample), cent, pt, xy * z);
    }
  }
  PROCESS_SWITCH(JV02DmesonAnalysis, processReference, "Charged-hadron v02(pT) closure", true);

  void processD0(MyCollisions::iterator const& collision, aod::BCs const&, aod::CFCollisions const&, aod::CFCollRefs const& refs, aod::CF2ProngTracks const& candidates, aod::CFTracks const& tracks)
  {
    float cent = -1.f;
    float xy = 0.f;
    const int step = passedEventStep(collision, cent);
    V02_FILL_QA("D0/hEventQa", step);
    const int cfId = cfCollisionIndex(refs);
    if (step < StepAccepted || cfId < 0 || !referenceProduct(collision, xy)) {
      return;
    }
    auto eventCandidates = candidates.sliceBy(mPerCFCand, cfId);
    auto eventTracks = tracks.sliceBy(mPerCFTrack, cfId);

    std::unordered_set<int64_t> daughters;
    if (cfgRejectDaughters) {
      for (const auto& cand : eventCandidates) {
        if (!isD0Decay(cand.decay())) {
          continue;
        }
        daughters.insert(static_cast<int64_t>(cand.cfTrackProng0Id()));
        daughters.insert(static_cast<int64_t>(cand.cfTrackProng1Id()));
      }
    }
    const float nch = countCharged(eventTracks, daughters, false);
    if (nch < cfgMinNch) {
      return;
    }
    histos.fill(HIST("D0/hEventQa"), 6.f);
    const int sample = mSampleDist(mRng);
    const float sampleCoord = static_cast<float>(sample);
    histos.fill(HIST("D0/hEvents"), cent);
    histos.fill(HIST("D0/Prof_XY"), cent, xy);
    histos.fill(HIST("D0/Prof_Nch"), cent, nch);

    const float wRaw = 1.f / nch;
    const float wCorr = xy / nch;
    for (const auto& cand : eventCandidates) {
      if (!isD0Decay(cand.decay())) {
        continue;
      }
      const float pt = cand.pt();
      if (pt < cfgPtMinD || pt >= cfgPtMaxD) {
        continue;
      }
      const float mass = cand.invMass();
      histos.fill(HIST("D0/hMassRaw"), cent, pt, mass, wRaw);
      histos.fill(HIST("D0/hMassCorr"), cent, pt, mass, wCorr);
      histos.fill(HIST("D0/hMassRawSub"), sampleCoord, cent, pt, mass, wRaw);
      histos.fill(HIST("D0/hMassCorrSub"), sampleCoord, cent, pt, mass, wCorr);
    }
  }
  PROCESS_SWITCH(JV02DmesonAnalysis, processD0, "D0 and D0bar v02(pT) from CF2ProngTracks", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<JV02DmesonAnalysis>(cfgc)};
}
