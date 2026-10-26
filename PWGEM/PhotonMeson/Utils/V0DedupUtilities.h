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

/// \file V0DedupUtilities.h
/// \brief helpers of the photon conversion builder: conflict-group solver of the deduplication and the geometry of the photon finder
/// \author Stefanie Mrozinski, stefanie.mrozinski@cern.ch

#ifndef PWGEM_PHOTONMESON_UTILS_V0DEDUPUTILITIES_H_
#define PWGEM_PHOTONMESON_UTILS_V0DEDUPUTILITIES_H_

#include "PWGEM/PhotonMeson/Utils/PCMUtilities.h"

#include <ReconstructionDataFormats/HelixHelper.h>
#include <ReconstructionDataFormats/Track.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

namespace o2::pwgem::photonmeson
{

//_______________________________________________________________________
/// \brief q_inv of two massless particles (photons)
/// \param px1 px of the first photon (GeV/c), same for py1, pz1
/// \param px2 px of the second photon (GeV/c), same for py2, pz2
/// \return sqrt(2 (|p1| |p2| - p1 . p2)) (GeV/c)
inline float qInv(const float px1, const float py1, const float pz1, const float px2, const float py2, const float pz2)
{
  return std::sqrt(std::max(0.f, 2.f * (std::hypot(px1, py1, pz1) * std::hypot(px2, py2, pz2) - px1 * px2 - py1 * py2 - pz1 * pz2)));
}

/// \brief result of solveConflictGroup
struct ConflictGroupSolution {
  uint32_t bestMask = 0;
  int bestCount = 0;
  float bestScore = 0.f;
  std::vector<int> bestWithCount;   // per candidate: size of the best subset that contains it
  std::vector<float> bestWithScore; // per candidate: score sum of that subset
};

//_______________________________________________________________________
/// \brief best leg-disjoint subset of a conflict group of the deduplication (mode 2), exact enumeration of all subsets
/// \param legMask per candidate: bit mask of its legs (group-local bits); at most 16 candidates
/// \param score per candidate: dedup score, lower = better
/// \param s0 objective: <= 0 maximum number of candidates, then minimum score sum; > 0 maximum sum of max(s0 - score, 1e-3),
/// ties by the score sum
/// \param memberBest also find, per candidate, the best subset that contains it (MC diagnostics)
/// \return mask of the best subset, its size and score sum, and the per-candidate results
inline ConflictGroupSolution solveConflictGroup(std::vector<uint64_t> const& legMask, std::vector<float> const& score, const float s0, const bool memberBest)
{
  const int nMembers = static_cast<int>(legMask.size());
  const bool weighted = s0 > 0.f;
  constexpr float MinWeight = 1e-3f;
  constexpr float WeightTolerance = 1e-6f;
  ConflictGroupSolution sol;
  sol.bestWithCount.assign(nMembers, -1);
  sol.bestWithScore.assign(nMembers, 0.f);
  float bestWeight = -1.f;
  for (uint32_t mask = 1; mask < (1u << nMembers); ++mask) {
    uint64_t legs = 0;
    float scoreSum = 0.f;
    float weightSum = 0.f;
    int count = 0;
    bool disjoint = true;
    for (int k = 0; k < nMembers; ++k) {
      if ((mask & (1u << k)) == 0) {
        continue;
      }
      if ((legs & legMask[k]) != 0) {
        disjoint = false;
        break;
      }
      legs |= legMask[k];
      scoreSum += score[k];
      weightSum += std::max(s0 - score[k], MinWeight);
      ++count;
    }
    if (!disjoint) {
      continue;
    }
    const bool better = weighted ? (weightSum > bestWeight + WeightTolerance || (std::fabs(weightSum - bestWeight) <= WeightTolerance && scoreSum < sol.bestScore))
                                 : (count > sol.bestCount || (count == sol.bestCount && scoreSum < sol.bestScore));
    if (better) {
      sol.bestCount = count;
      sol.bestScore = scoreSum;
      bestWeight = weightSum;
      sol.bestMask = mask;
    }
    if (memberBest) {
      for (int k = 0; k < nMembers; ++k) {
        if ((mask & (1u << k)) == 0) {
          continue;
        }
        if (count > sol.bestWithCount[k] || (count == sol.bestWithCount[k] && scoreSum < sol.bestWithScore[k])) {
          sol.bestWithCount[k] = count;
          sol.bestWithScore[k] = scoreSum;
        }
      }
    }
  }
  return sol;
}

/// \brief steps of the photon finder (pairing of free legs at the touching point of their helix circles); evalFinderPair returns the last step passed
enum FinderStep {
  kFindTglWindow = 0,
  kFindTime,
  kFindGap,
  kFindRTouch,
  kFindPointing,
  kFindAlongMotion,
  kFindZMatch,
  kFindCollision,
  kNFindStep
};

/// \brief collision for the photon finder: time (ns, relative to the first BC of the DF), index, vertex (cm)
struct FinderCollision {
  double t = 0.;
  int64_t id = -1;
  float x = 0.f, y = 0.f, z = 0.f;
};

/// \brief leg for the photon finder: track parameters, helix circle, TPC drift side and the time window (ns) of the track
struct FinderLeg {
  int64_t id = -1;
  float tgl = 0.f;
  bool noITS = false;
  bool movable = false;
  int side = 0;
  double t = 0., tMin = 0., tMax = 0.;
  o2::track::TrackPar par;
  o2::track::TrackAuxPar circle;
};

/// \brief cuts of the photon finder (cm) and the magnetic field (kG)
struct FinderCuts {
  float maxDeltaTgl = 0.f;
  float maxCircleGap = 0.f;
  float minRTPConly = 0.f;
  float maxRReqITS = 0.f;
  float rMargin = 0.f;
  float maxR = 0.f;
  float minCosXY = 0.f;
  float maxDeltaZ = 0.f;
  float maxDz = 0.f;
  float bz = 0.f;
};

//_______________________________________________________________________
/// \brief geometric test of one leg pair of the photon finder: tgl, time windows, circle gap, touching radius, xy pointing to the
/// beam, direction of motion, z match at the touching point (incl. the drift of TPC-only legs), collision by z
/// \param p positive leg
/// \param e negative leg
/// \param colls collisions of the DF, sorted by time
/// \param cuts cuts of the finder
/// \param xB x of the beam position (cm)
/// \param yB y of the beam position (cm)
/// \param vDrift TPC drift velocity (cm/ns)
/// \param colId output: index of the matched collision, -1 if none
/// \return last FinderStep passed, -1 if the tgl window fails; kFindCollision if a collision matches
inline int evalFinderPair(FinderLeg const& p, FinderLeg const& e, std::vector<FinderCollision> const& colls, FinderCuts const& cuts, const float xB, const float yB, const double vDrift, int64_t& colId)
{
  colId = -1;
  if (std::fabs(p.tgl - e.tgl) > cuts.maxDeltaTgl) {
    return -1;
  }
  const double tLo = std::max(p.tMin, e.tMin);
  const double tHi = std::min(p.tMax, e.tMax);
  if (tLo > tHi) {
    return kFindTglWindow;
  }
  const auto& cp = p.circle;
  const auto& ce = e.circle;
  const float dxC = ce.xC - cp.xC;
  const float dyC = ce.yC - cp.yC;
  const float dC = std::hypot(dxC, dyC);
  if (std::fabs(circleGap(cp, ce)) > cuts.maxCircleGap) {
    return kFindTime;
  }
  const float xT = (cp.xC * ce.rC + ce.xC * cp.rC) / (cp.rC + ce.rC);
  const float yT = (cp.yC * ce.rC + ce.yC * cp.rC) / (cp.rC + ce.rC);
  const float rTouch = std::hypot(xT, yT);
  const float rMin = (p.noITS && e.noITS ? cuts.minRTPConly : cuts.maxRReqITS) - cuts.rMargin;
  if (rTouch > cuts.maxR || rTouch < rMin) {
    return kFindGap;
  }
  const float dB = std::hypot(xT - xB, yT - yB);
  if (!(dC > 0.f) || !(dB > 0.f) || std::fabs((xT - xB) * (-dyC / dC) + (yT - yB) * (dxC / dC)) / dB < cuts.minCosXY) {
    return kFindRTouch;
  }
  auto nearest = [&](o2::track::TrackAuxPar const& c, float& x, float& y) {
    const float d = std::max(std::hypot(xT - c.xC, yT - c.yC), 1e-6f);
    x = c.xC + c.rC * (xT - c.xC) / d;
    y = c.yC + c.rC * (yT - c.yC) / d;
  };
  float xP = 0.f, yP = 0.f, xE = 0.f, yE = 0.f;
  nearest(cp, xP, yP);
  nearest(ce, xE, yE);
  const auto hp = helixPointAt(p.par, xP, yP, xB, yB, cuts.bz);
  const auto he = helixPointAt(e.par, xE, yE, xB, yB, cuts.bz);
  if (!hp.ok || !he.ok || hp.inward || he.inward) {
    return kFindPointing;
  }
  const double b = vDrift * (p.side - e.side);
  const double a = (hp.z - p.side * vDrift * p.t) - (he.z - e.side * vDrift * e.t);
  const double tBest = b == 0. ? tLo : std::clamp(-a / b, tLo, tHi);
  if (std::fabs(a + b * tBest) > cuts.maxDeltaZ) {
    return kFindAlongMotion;
  }
  const float ptG = std::hypot(hp.p[0] + he.p[0], hp.p[1] + he.p[1]);
  const float tglG = ptG > 0.f ? (hp.p[2] + he.p[2]) / ptG : 0.f;
  float best = std::numeric_limits<float>::max();
  for (auto it = std::lower_bound(colls.begin(), colls.end(), tLo, [](const FinderCollision& c, double t) { return c.t < t; }); it != colls.end() && it->t <= tHi; ++it) {
    const auto zPos = static_cast<float>(hp.z + p.side * vDrift * (it->t - p.t));
    const auto zEle = static_cast<float>(he.z + e.side * vDrift * (it->t - e.t));
    const float deltaZ = zPos - zEle;
    const float dz = (zPos * ce.rC + zEle * cp.rC) / (cp.rC + ce.rC) - std::hypot(xT - it->x, yT - it->y) * tglG - it->z;
    const float dist = std::hypot(dz, deltaZ);
    if (dist < best && std::fabs(dz) < cuts.maxDz && std::fabs(deltaZ) < cuts.maxDeltaZ) {
      best = dist;
      colId = it->id;
    }
  }
  return colId >= 0 ? kFindCollision : kFindZMatch;
}

/// \brief input features of the finder purity score (naive-Bayes log-likelihood ratio, finderScore_v1, PbPb MC LHC26a5a)
enum FinderScoreFeature {
  kScoreCosPA = 0, // 3D cosine of the pointing angle
  kScorePCA,       // distance of the legs at the conversion point (cm)
  kScoreMee,       // e+e- mass at the conversion point (GeV/c^2)
  kScoreGap,       // helix-circle gap of the legs (cm)
  kScoreDTgl,      // |tgl+ - tgl-|
  kScorePsiPair,
  kScoreMaxNSigmaEl, // max(|n sigma_e|) of the legs in the TPC
  kScoreMinNSigmaPi, // min(|n sigma_pi|) of the legs in the TPC
  kNScoreFeatures
};

/// \brief log-likelihood-ratio table of one feature: nBins uniform bins in [lo, hi), llr[0] underflow, llr[nBins + 1] overflow
/// (same bin lookup as TH1::FindBin)
struct FinderScoreTable {
  int nBins;
  double lo;
  double hi;
  const double* llr;
};

/// \brief log(p_true / p_fake) per bin of each feature (finderScore_v1)
inline constexpr std::array<double, 52> FinderLlrCosPA = {-0.42326231102291223, -3.0746401250262112, -2.902679468526141, -2.7737498628462491, -2.8653323764415854, -2.9977154358984874, -2.879287982624958, -2.9449559801189884, -2.8305718792936898, -2.7610493575753212, -2.6743115383729115, -2.739759837578176, -2.7445126982363175, -2.6010956990444813, -2.5377329273302194, -2.4413963521466235, -2.5007549142569676, -2.5484647707016048, -2.3339023465737707, -2.3706555897531327, -2.3350661347670365, -2.3987685637953868, -2.2562477335729185, -2.2430310349096305, -2.1429211398936303, -2.1233771903398493, -2.0022839108013391, -1.9558000169723206, -1.8729479034124974, -1.9065302845458256, -1.8542708039089619, -1.7230031706190478, -1.6490950919627059, -1.633800343346663, -1.5180730364736743, -1.389137677593709, -1.3181114339137541, -1.2615432549457208, -1.1021895055813662, -0.97913284177589832, -0.97669632321874755, -0.81070582095170196, -0.70030094035730539, -0.55130110346066707, -0.34392181466201621, -0.16444806113591265, 0.025741391801742678, 0.22392951134446418, 0.48449501950084528, 0.74029958080748082, 1.3544181524605163, 2.1416870464386246};
inline constexpr std::array<double, 32> FinderLlrPca = {-0.42328563731423596, 0.74173300442294698, 0.60064005460866021, 0.48774214670698213, 0.43172212684872524, 0.34949695425445604, 0.2876013953788461, 0.13114334603545452, -0.013681195276045386, -0.10436433242742214, -0.27099684361934573, -0.38200720615175005, -0.47403001970996045, -0.57728096929145511, -0.65985610396473016, -0.75403841222733481, -0.76619973024014831, -0.76182156236311549, -0.84274340404963699, -0.93501651159467125, -0.99602554213200045, -1.029448537116822, -1.1240438482723172, -1.124554873172243, -1.2439733445571921, -1.2077440689826793, -1.3081144381420802, -1.3626265850153199, -1.3892300573742717, -1.3978783977548588, -1.4463232066948644, -0.42328563731423596};
inline constexpr std::array<double, 27> FinderLlrMee = {-0.42329146929478079, 0.82940375143002754, 0.6225570701387575, 0.30541176737301567, -0.042962035242895512, -0.43745174446394713, -0.78901712806679281, -1.1834130413160895, -1.5486037315484942, -1.8446313857300292, -2.2387380284795446, -2.1811493868471543, -0.42329146929478079, -0.42329146929478079, -0.42329146929478079, -0.42329146929478079, -0.42329146929478079, -0.42329146929478079, -0.42329146929478079, -0.42329146929478079, -0.42329146929478079, -0.42329146929478079, -0.42329146929478079, -0.42329146929478079, -0.42329146929478079, -0.42329146929478079, -0.42329146929478079};
inline constexpr std::array<double, 22> FinderLlrGap = {-0.42329730143842487, -3.0127420263468552, -2.6803080781356927, -2.2310834132674926, -1.852433867664681, -1.4426351719772665, -1.1194240909923954, -0.72769983474283484, -0.24592785900768224, 0.40475570320106047, 1.4082550958014817, 1.1535415604659081, 0.43945228763843225, 0.14326171246634631, -0.19523294768817565, -0.39258756076294099, -0.41539680606824775, -0.42329730143842487, -0.42329730143842487, -0.42329730143842487, -0.42329730143842487, -0.42329730143842487};
inline constexpr std::array<double, 17> FinderLlrDTgl = {-0.42330313374517559, 0.84513809550928687, 0.46735826574627454, 0.13082343003051436, -0.12182979016333238, -0.33846868362980387, -0.52696466906551476, -0.72450796222957359, -0.89155673129371438, -1.0310133083500654, -1.2327015991375558, -1.3726345840013925, -1.4519781151564646, -1.7004785855213946, -1.7559986132128569, -1.818203125662764, -0.42330313374517559};
inline constexpr std::array<double, 18> FinderLlrPsiPair = {-0.42330196727077668, 1.5226081817845367, 0.42399589311642699, 1.6969615689293145, 1.6135799599902634, 0.44573587975283285, 0.50095693425255516, 0.43535965176674207, 0.0028547913604746104, -0.0059523726814127804, 0.61779960880440121, 0.82263503615419142, 0.81841116503800626, 0.99808371366038406, 0.57522686284035063, 0.80047346435133904, 2.521137011895664, 2.2392858597546761};
inline constexpr std::array<double, 10> FinderLlrMaxNSigmaEl = {-0.42331129924865812, 0.58307082727859216, 0.52708190318359982, 0.41704446150850366, 0.16536503067868297, -0.38368187287820599, -1.4189814933894278, -3.0185512358113016, -5.190101150118414, -0.42331129924865812};
inline constexpr std::array<double, 12> FinderLlrMinNSigmaPi = {-0.42330896621503916, -0.9609371208084968, -0.46431089172823747, -0.06050992822674324, 0.16077238768804142, 0.3920347340880036, 0.64255032840474002, 0.82818496003043784, 0.96447986209716363, 1.0700938177484094, 1.1082284119787249, 1.017209914244056};

/// \brief binning and LLR table of each feature, in the order of FinderScoreFeature
inline constexpr std::array<FinderScoreTable, kNScoreFeatures> FinderScoreTables = {{
  {50, 0.99000000953674316, 1, FinderLlrCosPA.data()},
  {30, 0, 3, FinderLlrPca.data()},
  {25, 0, 0.05000000074505806, FinderLlrMee.data()},
  {20, -5, 5, FinderLlrGap.data()},
  {15, 0, 0.15000000596046448, FinderLlrDTgl.data()},
  {16, -1.6000000238418579, 1.6000000238418579, FinderLlrPsiPair.data()},
  {8, 0, 4, FinderLlrMaxNSigmaEl.data()},
  {10, 0, 10, FinderLlrMinNSigmaPi.data()},
}};

//_______________________________________________________________________
/// \brief purity score of a finder candidate: log(p_true / p_fake) summed over the 8 features (naive Bayes)
/// \param x feature values in the order of FinderScoreFeature
/// \return score; > 0: more likely a true photon
inline float finderPurityScore(std::array<float, kNScoreFeatures> const& x)
{
  float score = 0.f;
  for (int i = 0; i < kNScoreFeatures; ++i) {
    const auto& t = FinderScoreTables[i];
    const double v = x[i];
    int bin = 0;
    if (v < t.lo) {
      bin = 0;
    } else if (!(v < t.hi)) {
      bin = t.nBins + 1;
    } else {
      bin = 1 + static_cast<int>(t.nBins * (v - t.lo) / (t.hi - t.lo));
    }
    score += t.llr[bin];
  }
  return score;
}

} // namespace o2::pwgem::photonmeson

#endif // PWGEM_PHOTONMESON_UTILS_V0DEDUPUTILITIES_H_
