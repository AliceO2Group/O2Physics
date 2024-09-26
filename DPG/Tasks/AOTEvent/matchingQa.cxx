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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsFT0/Digit.h"
#include "MetadataHelper.h"

using namespace o2;
using namespace o2::framework;

using BCsRun3 = soa::Join<aod::BCs, aod::Run3MatchedToBCSparse>;
using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra>;
using FullTracksIUwithLabels = soa::Join<aod::TracksIU, aod::TracksExtra, aod::McTrackLabels>;
float bcNS = o2::constants::lhc::LHCBunchSpacingNS;
int32_t nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;

MetadataHelper metadataInfo; // Metadata helper

struct MatchingQaTask {
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  Preslice<FullTracksIU> perCollision = aod::track::collisionId;
  Configurable<int> customOrbitOffset{"customOrbitOffset", 0, "customOrbitOffset for MC"};
  Configurable<bool> isLowFlux{"isLowFlux", 0, "1 - low flux (pp, pPb), 0 - high flux (PbPb)"};
  Configurable<bool> useTimeDiff{"useTimeDiff", 1, "use time difference for selection"};
  Configurable<bool> useVtxDiff{"useVtxDiff", 1, "use vertex difference for selection"};
  Configurable<bool> removeTOFmatches{"removeTOFmatches", 1, "remove TVX bcs matched to collisions with TOF tracks"};
  Configurable<bool> removeColsWithAmbiguousTOF{"removeColsWithAmbiguousTOF", 0, "remove collisions with ambiguous TOF signals"};
  int lastRun = -1;
  int64_t bcSOR = -1;     // global bc of the start of the first orbit
  int32_t nBCsPerTF = -1; // duration of TF in bcs, should be 128*3564 or 32*3564
  int32_t offsetITSROF = 64;
  int32_t nBCsPerITSROF = 198;
  std::vector<int> vFoundBCindex;
  std::vector<int> vNumITStracks;
  std::vector<int> vNumTOFtracks;
  std::vector<int> vNumTRDtracks;
  std::vector<int> vNumTPCtracks;

  bool isGoodBC(int64_t globalBC, bool fillHistos = 0)
  {
    // kNoTimeFrameBorder
    int64_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;
    if (fillHistos)
      histos.fill(HIST("hBcInTFall"), bcInTF);
    if (bcInTF < 300 || bcInTF > nBCsPerTF - 4000)
      return 0;
    if (fillHistos)
      histos.fill(HIST("hBcInTFcut"), bcInTF);
    // kNoITSROFrameBorder
    uint16_t bcInITSROF = (globalBC + nBCsPerOrbit - offsetITSROF) % nBCsPerITSROF;
    if (fillHistos)
      histos.fill(HIST("hBcInITSROFall"), bcInITSROF);
    if (bcInITSROF < 10 || bcInITSROF > nBCsPerITSROF - 20)
      return 0;
    if (fillHistos)
      histos.fill(HIST("hBcInITSROFcut"), bcInITSROF);
    return 1;
  }

  void init(InitContext&)
  {
    if (metadataInfo.isFullyDefined()) {
      if (!metadataInfo.isMC()) {
        doprocessMC.value = false;
      }
    }

    const AxisSpec axisNcontrib{isLowFlux ? 150 : 8000, 0., isLowFlux ? 150. : 8000., "n contributors"};
    const AxisSpec axisColTimeRes{1500, 0., 1500., "collision time resolution (ns)"};
    const AxisSpec axisFraction{1000, 0., 1., ""};
    const AxisSpec axisBcDiff{800, -400., 400., "bc diff"};
    const AxisSpec axisBcs{nBCsPerOrbit, 0., static_cast<float>(nBCsPerOrbit), "bc"};
    const AxisSpec axisMultT0C{200, 0., isLowFlux ? 1600. : 60000., "Rec. mult. T0C"};
    const AxisSpec axisZvtxDiff{200, -20., 20., "Zvtx difference, cm"};

    histos.add("hRecMultT0C", "", kTH1D, {axisMultT0C});

    histos.add("hRecMultT0CvsNcontrib", "", kTH2D, {axisMultT0C, axisNcontrib});
    histos.add("hRecMultT0CvsNcontribTPC", "", kTH2D, {axisMultT0C, axisNcontrib});
    histos.add("hRecMultT0CvsNcontribTOF", "", kTH2D, {axisMultT0C, axisNcontrib});
    histos.add("hRecMultT0CvsNcontribTRD", "", kTH2D, {axisMultT0C, axisNcontrib});

    histos.add("hBCsITS", "", kTH1F, {axisBcs});
    histos.add("hNcontribCandidates", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribSigma", "", kTH2F, {axisNcontrib, axisColTimeRes});

    histos.add("hNcontribAll", "", kTH1F, {axisNcontrib});

    histos.add("hNcontribCol", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribColTOF", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribColTRD", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribColTPC", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribColITS", "", kTH1F, {axisNcontrib});

    histos.add("hNcontribAcc", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribAccTOF", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribAccTRD", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribAccTPC", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribAccITS", "", kTH1F, {axisNcontrib});

    histos.add("hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("hColBcDiffVsNcontribTOF", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("hColBcDiffVsNcontribTRD", "", kTH2F, {axisNcontrib, axisBcDiff});
    histos.add("hColBcDiffVsNcontribTPC", "", kTH2F, {axisNcontrib, axisBcDiff});

    histos.add("hZvtxDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisZvtxDiff});
    histos.add("hZvtxDiffVsNcontribTOF", "", kTH2F, {axisNcontrib, axisZvtxDiff});
    histos.add("hZvtxDiffVsNcontribTPC", "", kTH2F, {axisNcontrib, axisZvtxDiff});
    histos.add("hZvtxDiffVsNcontribTRD", "", kTH2F, {axisNcontrib, axisZvtxDiff});

    histos.add("hNcontribUnambiguous", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribUnambiguousTOF", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribUnambiguousTRD", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribUnambiguousTPC", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribUnambiguousITS", "", kTH1F, {axisNcontrib});

    histos.add("hNcontribMis", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribMisTOF", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribMisTRD", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribMisTPC", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribMisITS", "", kTH1F, {axisNcontrib});

    histos.add("hNcontribColMostlyOk", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribColMostlyOkTOF", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribColMostlyOkTPC", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribColMostlyOkTRD", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribColMostlyOkITS", "", kTH1F, {axisNcontrib});

    histos.add("hNcontribColMostlyOkMis", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribColMostlyOkMisTOF", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribColMostlyOkMisTPC", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribColMostlyOkMisTRD", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribColMostlyOkMisITS", "", kTH1F, {axisNcontrib});

    histos.add("hNcontribAllContribAll", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribAllContribWrong", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribAllFractionWrong", "", kTH2F, {axisNcontrib, axisFraction});
    histos.add("hNcontribTvxMostlyOk", "", kTH1F, {axisNcontrib});
  }

  int32_t findClosest(int64_t globalBC, std::map<int64_t, int32_t>& bcs)
  {
    auto it = bcs.lower_bound(globalBC);
    int64_t bc1 = it->first;
    int32_t index1 = it->second;
    if (it != bcs.begin())
      --it;
    int64_t bc2 = it->first;
    int32_t index2 = it->second;
    int64_t dbc1 = std::abs(bc1 - globalBC);
    int64_t dbc2 = std::abs(bc2 - globalBC);
    return (dbc1 <= dbc2) ? index1 : index2;
  }

  void process(aod::Collisions const& cols, FullTracksIU const& tracks, BCsRun3 const& bcs, aod::FT0s const& ft0s)
  {
    int run = bcs.iteratorAt(0).runNumber();
    if (run != lastRun) {
      lastRun = run;
      auto runDuration = ccdb->getRunDuration(run, true);
      int64_t tsSOR = runDuration.first;
      auto ctpx = ccdb->getForTimeStamp<std::vector<Long64_t>>("CTP/Calib/OrbitReset", tsSOR);
      int64_t tsOrbitReset = (*ctpx)[0];
      uint32_t nOrbitsPerTF = run < 534133 ? 128 : 32;
      int64_t orbitSOR = (tsSOR * 1000 - tsOrbitReset) / o2::constants::lhc::LHCOrbitMUS;
      orbitSOR = orbitSOR / nOrbitsPerTF * nOrbitsPerTF;
      bcSOR = orbitSOR * nBCsPerOrbit + customOrbitOffset * nBCsPerOrbit;
      nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;
      nBCsPerITSROF = (run >= 543437 && run <= 545367) ? 594 : 198;
      const AxisSpec axisBcsInTF{nBCsPerTF, 0., static_cast<float>(nBCsPerTF), "bc"};
      const AxisSpec axisBcsInITSROF{nBCsPerITSROF, 0., static_cast<float>(nBCsPerITSROF), "bc"};
      histos.add("hBcInTFall", "", kTH1F, {axisBcsInTF});
      histos.add("hBcInTFcut", "", kTH1F, {axisBcsInTF});
      histos.add("hBcInITSROFall", "", kTH1F, {axisBcsInITSROF});
      histos.add("hBcInITSROFcut", "", kTH1F, {axisBcsInITSROF});
    }

    int nCols = cols.size();
    vFoundBCindex.resize(nCols);
    vNumITStracks.resize(nCols);
    vNumTOFtracks.resize(nCols);
    vNumTRDtracks.resize(nCols);
    vNumTPCtracks.resize(nCols);
    std::fill(vFoundBCindex.begin(), vFoundBCindex.end(), -1);
    std::fill(vNumITStracks.begin(), vNumITStracks.end(), 0);
    std::fill(vNumTOFtracks.begin(), vNumTOFtracks.end(), 0);
    std::fill(vNumTRDtracks.begin(), vNumTRDtracks.end(), 0);
    std::fill(vNumTPCtracks.begin(), vNumTPCtracks.end(), 0);

    std::vector<float> vTOFtracksSumWeightedTimes(cols.size(), 0);
    std::vector<float> vTRDtracksSumWeightedTimes(cols.size(), 0);
    std::vector<float> vTPCtracksSumWeightedTimes(cols.size(), 0);
    std::vector<float> vITStracksSumWeightedTimes(cols.size(), 0);
    std::vector<float> vTOFtracksSumWeights(cols.size(), 0);
    std::vector<float> vTRDtracksSumWeights(cols.size(), 0);
    std::vector<float> vTPCtracksSumWeights(cols.size(), 0);
    std::vector<float> vITStracksSumWeights(cols.size(), 0);
    std::vector<float> vMinTimeTOFtracks(cols.size(), 10000);
    std::vector<float> vMaxTimeTOFtracks(cols.size(), -10000);
    std::vector<float> vWeightedSigma(cols.size(), 0);
    std::map<int64_t, int32_t> mapGlobalBcWithTVX;
    std::map<int64_t, float> mapGlobalBcVtxZ;
    std::map<int64_t, float> mapGlobalBcMultT0C;
    std::map<int64_t, float> mapGlobalBcVtxZ2;

    int nBCs = bcs.size();
    std::vector<uint64_t> vGlobalBCs(nBCs, 0);

    for (auto& bc : bcs) {
      vGlobalBCs[bc.globalIndex()] = bc.globalBC();
    }

    for (auto& ft0 : ft0s) {
      auto bc = ft0.bc_as<BCsRun3>();
      int64_t globalBC = bc.globalBC();
      if (ft0.triggerMask() & BIT(o2::ft0::Triggers::bitVertex)) {
        mapGlobalBcWithTVX[globalBC] = bc.globalIndex();
        mapGlobalBcVtxZ[globalBC] = ft0.posZ();
        mapGlobalBcVtxZ2[globalBC] = ft0.posZ();
        mapGlobalBcMultT0C[globalBC] = ft0.sumAmpC();
      }
    }

    for (auto& track : tracks) {
      // DataFormats/Detectors/GlobalTracking/include/DataFormatsGlobalTracking/RecoContainerCreateTracksVariadic.h
      // Time for different track types:
      // ITS-TPC-TRD-TOF: time from TOF +/- 10 ns
      // ITS-TPC-TRD: time from TRD +/- 5 ns
      // ITS-TPC: time from ITS-TPC matching

      // ITS and colId requirements are redundant for contributors
      int32_t colId = track.collisionId();

      if (!track.isPVContributor() || colId < 0 || !track.hasITS())
        continue;

      float trackTime = track.trackTime();
      float trackTimeRes = track.trackTimeRes();
      float w = 1. / (trackTimeRes * trackTimeRes);
      if (track.hasTOF()) {
        vNumTOFtracks[colId]++;
        vTOFtracksSumWeightedTimes[colId] += trackTime * w;
        vTOFtracksSumWeights[colId] += w;
        if (vMinTimeTOFtracks[colId] > trackTime)
          vMinTimeTOFtracks[colId] = trackTime;
        if (vMaxTimeTOFtracks[colId] < trackTime)
          vMaxTimeTOFtracks[colId] = trackTime;
      } else if (track.hasTRD()) {
        vNumTRDtracks[colId]++;
        vTRDtracksSumWeightedTimes[colId] += trackTime * w;
        vTRDtracksSumWeights[colId] += w;
      } else if (track.hasTPC()) {
        vNumTPCtracks[colId]++;
        vTPCtracksSumWeightedTimes[colId] += trackTime * w;
        vTPCtracksSumWeights[colId] += w;
      } else {
        vNumITStracks[colId]++;
        vITStracksSumWeightedTimes[colId] += trackTime * w;
        vITStracksSumWeights[colId] += w;
      }
    }

    for (auto& col : cols) {
      int32_t colId = col.globalIndex();

      if (vNumTOFtracks[colId] == 0)
        continue;
      auto bc = col.bc_as<BCsRun3>();
      int64_t globalBC = bc.globalBC();
      // todo: bypass ambiguous collisions with TOF tracks pointing to different bcs
      float weightedTime = vTOFtracksSumWeightedTimes[colId] / vTOFtracksSumWeights[colId];
      int64_t tofGlobalBC = globalBC + TMath::Nint(weightedTime / bcNS);
      int32_t foundBC = findClosest(tofGlobalBC, mapGlobalBcWithTVX);
      // todo: check what to do if foundBC is too far from tofGlobalBC
      if (fabs(bcs.iteratorAt(foundBC).globalBC() - tofGlobalBC) > 0) {
        foundBC = -1;
      }
      vFoundBCindex[colId] = foundBC;
      if (removeTOFmatches)
        mapGlobalBcVtxZ.erase(globalBC);
    }

    // second loop to match collisions with poor time resolution
    for (auto& col : cols) {
      int32_t colId = col.globalIndex();
      if (vNumTOFtracks[colId] > 0)
        continue;
      if (vNumTRDtracks[colId] > 0)
        continue;
      if (vNumTPCtracks[colId] == 0)
        continue;
      auto bc = col.bc_as<BCsRun3>();
      int64_t globalBC = bc.globalBC();
      float weightedTime = vTPCtracksSumWeightedTimes[colId] / vTPCtracksSumWeights[colId];
      float weightedSigma = sqrt(1. / vTPCtracksSumWeights[colId]);
      int64_t deltaBC = std::ceil(weightedSigma / bcNS * 3);
      int64_t tpcGlobalBC = globalBC + TMath::Nint(weightedTime / bcNS);

      // if (deltaBC>20) deltaBC = 20;

      int64_t minBC = tpcGlobalBC - deltaBC;
      int64_t maxBC = tpcGlobalBC + deltaBC;

      if (vNumITStracks[colId] > 0) {
        float weightedTimeITS = vITStracksSumWeightedTimes[colId] / vITStracksSumWeights[colId];
        int64_t itsGlobalBC = globalBC + TMath::Nint(weightedTimeITS / bcNS);
        int64_t minROF = (itsGlobalBC - offsetITSROF) / nBCsPerITSROF * nBCsPerITSROF + offsetITSROF;
        int64_t maxROF = minROF + nBCsPerITSROF;
        histos.fill(HIST("hBCsITS"), minROF % 3564);
        minBC = minBC < minROF ? minROF : minBC;
        maxBC = maxBC > maxROF ? maxROF : maxBC;
        // LOGP(info,"{} {} {} {}",minBC, maxBC, minROF, maxROF);
      }
      if (minBC > maxBC)
        continue;

      int32_t nContrib = col.numContrib();
      float zVtxCol = col.posZ();
      float zVtxSigma = 2.7 * pow(nContrib, -0.466) + 0.024;
      zVtxSigma += 1.0; // additional uncertainty due to imperfectections of FT0 time calibration

      // QA
      vWeightedSigma[colId] = weightedSigma;

      // todo: check upper bound
      auto itMin = mapGlobalBcVtxZ.lower_bound(minBC);
      auto itMax = mapGlobalBcVtxZ.upper_bound(maxBC);

      float bestChi2 = 1e+10;
      int64_t globalBcBest = 0;

      int nCandidates = 0;
      for (std::map<int64_t, float>::iterator it = itMin; it != itMax; ++it) {
        float zVtxFT0 = it->second;
        float zVtxDiff = zVtxFT0 - zVtxCol;
        float timeDiff = bcNS * (tpcGlobalBC - globalBC);
        float chi2 = 0;
        chi2 += useVtxDiff ? pow(zVtxDiff / zVtxSigma, 2) : 0.;
        chi2 += useTimeDiff ? pow(timeDiff / weightedSigma, 2) : 0.;

        if (chi2 < bestChi2) {
          bestChi2 = chi2;
          globalBcBest = it->first;
        }
        nCandidates++;
      }
      if (nCandidates > 100)
        LOGP(info, "{} {}", minBC, maxBC);

      histos.fill(HIST("hNcontribCandidates"), nContrib, nCandidates);

      if (globalBcBest != 0)
        vFoundBCindex[colId] = mapGlobalBcWithTVX[globalBcBest];

    } // second loop

    // QA
    for (auto& ft0 : ft0s) {
      histos.fill(HIST("hRecMultT0C"), ft0.sumAmpC());
    }

    for (auto& col : cols) {
      int64_t globalBC = col.bc_as<BCsRun3>().globalBC();
      if (!isGoodBC(globalBC, 1))
        continue;

      int32_t colId = col.globalIndex();
      int32_t nContrib = col.numContrib();
      // float timeRes = col.collisionTimeRes();
      int32_t foundBC = vFoundBCindex[colId];
      bool isGoodTOF = vMaxTimeTOFtracks[colId] - vMinTimeTOFtracks[colId] < 50;
      bool isFoundTVX = foundBC >= 0;

      float zVtxDiff = 1e+10;
      float multT0C = 0;

      if (foundBC >= 0 && foundBC < bcs.size()) {
        auto bc = bcs.iteratorAt(foundBC);
        // LOGP(info,"{}",bc.has_ft0());
        if (bc.has_ft0()) {
          zVtxDiff = bc.ft0().posZ() - col.posZ();
          multT0C = bc.ft0().sumAmpC();
        }
      }

      histos.fill(HIST("hNcontribAll"), nContrib);
      if (removeColsWithAmbiguousTOF && !isGoodTOF) {
        continue;
      }

      histos.fill(HIST("hNcontribCol"), nContrib);
      if (isFoundTVX) {
        histos.fill(HIST("hNcontribAcc"), nContrib);
        histos.fill(HIST("hRecMultT0CvsNcontrib"), multT0C, nContrib);
        histos.fill(HIST("hZvtxDiffVsNcontrib"), nContrib, zVtxDiff);
      }

      // search for nearest ft0a&ft0c entry
      int32_t indexClosestTVX = findClosest(globalBC, mapGlobalBcWithTVX);
      int bcDiff = static_cast<int>(globalBC - vGlobalBCs[indexClosestTVX]);
      histos.fill(HIST("hColBcDiffVsNcontrib"), nContrib, bcDiff);

      if (vNumTOFtracks[colId] > 0) {
        histos.fill(HIST("hNcontribColTOF"), nContrib);
        if (isFoundTVX) {
          histos.fill(HIST("hNcontribAccTOF"), nContrib);
          histos.fill(HIST("hRecMultT0CvsNcontribTOF"), multT0C, nContrib);
          histos.fill(HIST("hZvtxDiffVsNcontribTOF"), nContrib, zVtxDiff);
        }
        histos.fill(HIST("hColBcDiffVsNcontribTOF"), nContrib, bcDiff);
      } else if (vNumTRDtracks[colId] > 0) {
        histos.fill(HIST("hNcontribColTRD"), nContrib);
        if (isFoundTVX) {
          histos.fill(HIST("hNcontribAccTRD"), nContrib);
          histos.fill(HIST("hRecMultT0CvsNcontribTRD"), multT0C, nContrib);
          histos.fill(HIST("hZvtxDiffVsNcontribTRD"), nContrib, zVtxDiff);
        }
        histos.fill(HIST("hColBcDiffVsNcontribTRD"), nContrib, bcDiff);
      } else if (vNumTPCtracks[colId] > 0) {
        histos.fill(HIST("hNcontribSigma"), nContrib, vWeightedSigma[colId]);
        histos.fill(HIST("hNcontribColTPC"), nContrib);
        histos.fill(HIST("hColBcDiffVsNcontribTPC"), nContrib, bcDiff);
        if (isFoundTVX) {
          histos.fill(HIST("hNcontribAccTPC"), nContrib);
          histos.fill(HIST("hRecMultT0CvsNcontribTPC"), multT0C, nContrib);
          histos.fill(HIST("hZvtxDiffVsNcontribTPC"), nContrib, zVtxDiff);
        }
      } else if (vNumITStracks[colId] > 0) {
        histos.fill(HIST("hNcontribColITS"), nContrib);
        if (isFoundTVX)
          histos.fill(HIST("hNcontribAccITS"), nContrib);
      }
    }
  }

  void processMC(
    aod::McCollisions const& mcCols,
    soa::Join<aod::Collisions, aod::McCollisionLabels> const& cols,
    FullTracksIUwithLabels const& tracks,
    BCsRun3 const& /*bcs*/,
    aod::FT0s const& /*ft0s*/,
    aod::McParticles const& mcParts)
  {

    std::vector<int> vLabel(cols.size(), -1);
    std::vector<bool> vIsAmbiguousLabel(cols.size(), 0);
    std::vector<int> vNumWrongContributors(cols.size(), 0);

    for (auto& track : tracks) {
      if (!track.isPVContributor())
        continue;
      int32_t colId = track.collisionId();
      auto col = cols.iteratorAt(colId);
      int32_t mcColIdFromCollision = col.mcCollisionId();
      if (mcColIdFromCollision < 0)
        continue;
      int mcId = track.mcParticleId();
      if (mcId < 0 || mcId >= mcParts.size())
        continue;
      auto mcPart = mcParts.iteratorAt(mcId);
      int32_t mcColId = mcPart.mcCollisionId();
      if (mcColId < 0)
        continue;
      if (mcColIdFromCollision != mcColId)
        vNumWrongContributors[colId]++;

      if (vLabel[colId] != -1 && vLabel[colId] != mcColId) {
        vIsAmbiguousLabel[colId] = 1;
      }
      vLabel[colId] = mcColId;
    }

    for (auto& col : cols) {
      int64_t globalBC = col.bc_as<BCsRun3>().globalBC();
      if (!isGoodBC(globalBC))
        continue;

      int32_t colId = col.globalIndex();
      int32_t nContrib = col.numContrib();

      histos.fill(HIST("hNcontribAllContribAll"), nContrib, nContrib);
      histos.fill(HIST("hNcontribAllContribWrong"), nContrib, vNumWrongContributors[colId]);
      histos.fill(HIST("hNcontribAllFractionWrong"), nContrib, static_cast<float>(vNumWrongContributors[colId]) / nContrib);

      if (static_cast<float>(vNumWrongContributors[colId]) / nContrib > 0.1)
        continue;

      histos.fill(HIST("hNcontribColMostlyOk"), nContrib);
      if (vNumTOFtracks[colId] > 0) {
        histos.fill(HIST("hNcontribColMostlyOkTOF"), nContrib);
      } else if (vNumTPCtracks[colId] > 0) {
        histos.fill(HIST("hNcontribColMostlyOkTPC"), nContrib);
      } else if (vNumTRDtracks[colId] > 0) {
        histos.fill(HIST("hNcontribColMostlyOkTRD"), nContrib);
      } else if (vNumITStracks[colId] > 0) {
        histos.fill(HIST("hNcontribColMostlyOkITS"), nContrib);
      }

      int32_t foundBC = vFoundBCindex[colId];
      int32_t mcColId = vLabel[colId];
      auto mcCol = mcCols.iteratorAt(mcColId);
      auto mcBC = mcCol.bc_as<BCsRun3>();
      // int64_t mcGlobalBC = mcBC.globalBC();
      bool isMcTVX = mcBC.has_ft0() ? mcBC.ft0().triggerMask() & BIT(o2::ft0::Triggers::bitVertex) : 0;
      if (isMcTVX)
        histos.fill(HIST("hNcontribTvxMostlyOk"), nContrib);

      if (foundBC >= 0 && foundBC != mcBC.globalIndex()) {
        // Analyse mismatches
        histos.fill(HIST("hNcontribColMostlyOkMis"), nContrib);
        if (vNumTOFtracks[colId] > 0) {
          histos.fill(HIST("hNcontribColMostlyOkMisTOF"), nContrib);
        } else if (vNumTPCtracks[colId] > 0) {
          histos.fill(HIST("hNcontribColMostlyOkMisTPC"), nContrib);
        } else if (vNumTRDtracks[colId] > 0) {
          histos.fill(HIST("hNcontribColMostlyOkMisTRD"), nContrib);
        } else if (vNumITStracks[colId] > 0) {
          histos.fill(HIST("hNcontribColMostlyOkMisITS"), nContrib);
        }
      }

      if (vIsAmbiguousLabel[colId])
        continue;

      histos.fill(HIST("hNcontribUnambiguous"), nContrib);
      if (vNumTOFtracks[colId] > 0) {
        histos.fill(HIST("hNcontribUnambiguousTOF"), nContrib);
      } else if (vNumTPCtracks[colId] > 0) {
        histos.fill(HIST("hNcontribUnambiguousTPC"), nContrib);
      } else if (vNumTRDtracks[colId] > 0) {
        histos.fill(HIST("hNcontribUnambiguousTRD"), nContrib);
      } else if (vNumITStracks[colId] > 0) {
        histos.fill(HIST("hNcontribUnambiguousITS"), nContrib);
      }

      if (foundBC >= 0 && foundBC != mcBC.globalIndex()) {
        // Analyse mismatches
        histos.fill(HIST("hNcontribMis"), nContrib);
        if (vNumTOFtracks[colId] > 0) {
          histos.fill(HIST("hNcontribMisTOF"), nContrib);
        } else if (vNumTPCtracks[colId] > 0) {
          histos.fill(HIST("hNcontribMisTPC"), nContrib);
        } else if (vNumTRDtracks[colId] > 0) {
          histos.fill(HIST("hNcontribMisTRD"), nContrib);
        } else if (vNumITStracks[colId] > 0) {
          histos.fill(HIST("hNcontribMisITS"), nContrib);
        }
      }
    }
  }

  PROCESS_SWITCH(MatchingQaTask, processMC, "", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Parse the metadata
  metadataInfo.initMetadata(cfgc);

  return WorkflowSpec{
    adaptAnalysisTask<MatchingQaTask>(cfgc)};
}
