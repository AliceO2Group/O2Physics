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
/// \file upcPileUp.cxx
/// \brief Pile-up probability in FT0 for UPC analyses, per colliding bunch and time frame
///
/// For each colliding bunch crossing it counts
///   - denominator A: crossings with no activity on the C side (FT0C, FDDC, ZNC)
///   - numerator A  : of those, the ones with a signal in FT0A above threshold
/// and the mirror for the C side (A side empty: FT0A, FV0A, FDDA, ZNA).
/// The ratio numerator / denominator is the pile-up probability.
///
/// The BCs table only stores crossings with some activity. Completely empty
/// crossings are not in it, but they belong in the denominator. The task
/// therefore computes the number of analysed crossings of every colliding
/// bunch in every time frame, applying the same time-frame and ITS ROF border
/// cuts as the event selection, and the denominator is that number minus the
/// stored crossings with activity on the other side.
///
/// \author David Arevalos <david.arevalos.bueckert@cern.ch>

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonDataFormat/BunchFilling.h>
#include <DataFormatsITSMFT/DPLAlpideParam.h>
#include <DataFormatsParameters/AggregatedRunInfo.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>
#include <TString.h>

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

using BCsTSsSels = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

struct UpcPileUp {

  // services
  Service<o2::ccdb::BasicCCDBManager> ccdb{};

  // histograms, booked per run and kept in the maps by name
  HistogramRegistry bcTH1Registry{"bcTH1Registry", {}};
  std::map<std::string, std::shared_ptr<TH1>> bcTH1Pointers;
  HistogramRegistry bcTH2Registry{"bcTH2Registry", {}};
  std::map<std::string, std::shared_ptr<TH2>> bcTH2Pointers;

  // histograms of the current run, cached so that the loops do no name lookup
  std::shared_ptr<TH1> hBcSel, hTF, hBorderCheck;
  std::shared_ptr<TH2> hCross, hVtx, hVtxCe, hDenA, hDenC;
  std::vector<std::shared_ptr<TH2>> hNumA, hNumC, hNumAZn, hNumCZn;

  // filling scheme
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternB;
  std::vector<int> bcbIdx;  // index of each colliding bunch, -1 otherwise
  std::vector<int> bcbList; // bcID within the orbit of each colliding bunch
  int nbcB = 0;

  // ITS ROF info
  int rofLength = -1; // number of bcs per ROF
  int rofShift = -1;  // bc shift of the ITS ROF

  // border margins, as in the event selection
  int tfStartMargin = 0;
  int tfEndMargin = 0;
  int rofStartMargin = 0;
  int rofEndMargin = 0;

  // run info
  int runNumberBc = 0; // run number used to process BCs
  int64_t tsMid = 0;   // middle of the run, to access run-wise CCDB entries
  int64_t orbitsPerTF = 0;
  int64_t bcSOR = 0;         // first bc of the first orbit
  int64_t nBCsPerTF = 0;     // duration of TF in bcs
  int64_t nTF = 0;           // number of time frames in run
  std::set<int64_t> seenTFs; // time frames already counted in this run

  // thresholds used, filled in init from the configurable
  std::vector<float> thresholds;
  int nThresholds = 0;

  // analysed crossings of each colliding bunch in the current time frame
  std::vector<double> nCrossTF;

  // histogram binning constants
  static constexpr int NBinsBcSel = 4;
  static constexpr int NBinsBorderCheck = 2;

  // FDD channel pairing, as in upcVmRof
  static constexpr int NFddPairs = 4;
  static constexpr int FddPairOffset = 4;

  // configurables
  Configurable<int> tfPerBin{"tfPerBin", 10000, "timeframes per bin, 1e4 means some 28 s"};
  Configurable<int> ft0VtxIdx{"ft0VtxIdx", 2, "index of the FT0VTX input in the CTP input mask"};
  Configurable<int> ft0CeIdx{"ft0CeIdx", 4, "index of the FT0CE input in the CTP input mask"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  // signal and veto definition
  Configurable<std::vector<float>> ampThresholds{"ampThresholds", {100.0f, 150.0f, 200.0f}, "amplitude thresholds defining a signal in FT0"};
  Configurable<float> maxAbsTimeFT0{"maxAbsTimeFT0", 4.0, "max |time| in a FT0 side (ns)"};
  Configurable<float> maxAbsTimeFV0{"maxAbsTimeFV0", 4.0, "max |time| in FV0A (ns)"};
  Configurable<float> maxAbsTimeFDD{"maxAbsTimeFDD", 4.0, "max |time| in a FDD side (ns)"};
  Configurable<float> maxAbsTimeZN{"maxAbsTimeZN", 2.0, "max |time| in a ZN (ns)"};
  Configurable<bool> fillZnVariation{"fillZnVariation", true, "also fill the numerators asking in addition for a signal in ZNA or ZNC"};

  // border margins, must match the ones used by the event selection; -1 takes them from CCDB
  Configurable<int> confTFStartBorderMargin{"confTFStartBorderMargin", -1, "bcs cut at the start of the TF, -1 from CCDB"};
  Configurable<int> confTFEndBorderMargin{"confTFEndBorderMargin", -1, "bcs cut at the end of the TF, -1 from CCDB"};
  Configurable<int> confROFStartBorderMargin{"confROFStartBorderMargin", -1, "bcs cut at the start of the ITS ROF, -1 from CCDB"};
  Configurable<int> confROFEndBorderMargin{"confROFEndBorderMargin", -1, "bcs cut at the end of the ITS ROF, -1 from CCDB"};

  //--------------------------------------------------------------------------------
  // general run information
  void getRunInfo(int run)
  {
    auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(ccdb->instance(), run);
    tsMid = runInfo.sor / 2 + runInfo.eor / 2;
    orbitsPerTF = runInfo.orbitsPerTF;
    bcSOR = runInfo.orbitSOR * o2::constants::lhc::LHCMaxBunches;
    nBCsPerTF = orbitsPerTF * o2::constants::lhc::LHCMaxBunches;
    nTF = (runInfo.orbitEOR - runInfo.orbitSOR + orbitsPerTF - 1) / orbitsPerTF; // rounded up
  } // end getRunInfo()

  //--------------------------------------------------------------------------------
  // filling scheme
  void getFillingScheme()
  {
    auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", tsMid);
    auto beamPatternA = grplhcif->getBunchFilling().getBeamPattern(0);
    auto beamPatternC = grplhcif->getBunchFilling().getBeamPattern(1);
    bcPatternB = beamPatternA & beamPatternC;
    bcbIdx.assign(o2::constants::lhc::LHCMaxBunches, -1);
    bcbList.clear();
    nbcB = 0;
    for (int i = 0; i < o2::constants::lhc::LHCMaxBunches; i++) {
      if (bcPatternB.test(i)) {
        bcbIdx[i] = nbcB;
        bcbList.push_back(i);
        nbcB++;
      }
    }
    nCrossTF.assign(nbcB, 0.);
  } // end getFillingScheme()

  //--------------------------------------------------------------------------------
  // ITS ROF info and border margins, taken as in the event selection
  void getBorderInfo()
  {
    auto alppar = ccdb->getForTimeStamp<o2::itsmft::DPLAlpideParam<0>>("ITS/Config/AlpideParam", tsMid);
    rofShift = alppar->roFrameBiasInBC;
    rofLength = alppar->roFrameLengthInBC;

    auto par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", tsMid);
    tfStartMargin = confTFStartBorderMargin < 0 ? par->fTimeFrameStartBorderMargin : confTFStartBorderMargin.value;
    tfEndMargin = confTFEndBorderMargin < 0 ? par->fTimeFrameEndBorderMargin : confTFEndBorderMargin.value;
    rofStartMargin = confROFStartBorderMargin < 0 ? par->fITSROFrameStartBorderMargin : confROFStartBorderMargin.value;
    rofEndMargin = confROFEndBorderMargin < 0 ? par->fITSROFrameEndBorderMargin : confROFEndBorderMargin.value;
    LOGP(info, "border margins: TF {} {}, ITS ROF {} {}", tfStartMargin, tfEndMargin, rofStartMargin, rofEndMargin);
  } // end getBorderInfo()

  //--------------------------------------------------------------------------------
  int64_t getBcWithinOrbit(int64_t globalBC)
  {
    return (globalBC % o2::constants::lhc::LHCMaxBunches);
  }

  //--------------------------------------------------------------------------------
  int64_t getTimeFrame(int64_t globalBC)
  {
    return (globalBC - bcSOR) / nBCsPerTF;
  }

  //--------------------------------------------------------------------------------
  // same conditions as kNoTimeFrameBorder and kNoITSROFrameBorder in the event selection
  bool passesBorders(int64_t globalBC)
  {
    int64_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;
    if (!(bcInTF > tfStartMargin && bcInTF < nBCsPerTF - tfEndMargin)) {
      return false;
    }
    int64_t bcInRof = (globalBC + o2::constants::lhc::LHCMaxBunches - rofShift) % rofLength;
    return (bcInRof > rofStartMargin && bcInRof < rofLength - rofEndMargin);
  }

  //--------------------------------------------------------------------------------
  // selection flags of a stored bc
  bool checkBcFlags(const auto& bc)
  {
    hBcSel->Fill(0);
    if (!bc.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    hBcSel->Fill(1);
    if (!bc.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    hBcSel->Fill(2);
    return true;
  } // end checkBcFlags()

  //--------------------------------------------------------------------------------
  // analysed crossings of every colliding bunch in a time frame: all of them,
  // stored or not, that pass the same border cuts as the stored ones.
  // They go into the crossings histogram and, as a start, into both
  // denominators; stored crossings with activity on the other side are removed
  // later, one by one.
  void fillCrossings(int64_t tf)
  {
    const int64_t tfStartBC = bcSOR + tf * nBCsPerTF;
    std::fill(nCrossTF.begin(), nCrossTF.end(), 0.);
    for (int64_t orbit = 0; orbit < orbitsPerTF; orbit++) {
      for (const auto& bcInOrbit : bcbList) {
        int64_t globalBC = tfStartBC + orbit * o2::constants::lhc::LHCMaxBunches + bcInOrbit;
        if (passesBorders(globalBC)) {
          nCrossTF[bcbIdx[bcInOrbit]] += 1.;
        }
      }
    }
    for (int ib = 0; ib < nbcB; ib++) {
      if (nCrossTF[ib] > 0.) {
        hCross->Fill(tf, ib, nCrossTF[ib]);
        hDenA->Fill(tf, ib, nCrossTF[ib]);
        hDenC->Fill(tf, ib, nCrossTF[ib]);
      }
    }
  } // end fillCrossings()

  //--------------------------------------------------------------------------------
  // FDD amplitude of a side, using the channel pairing of upcVmRof
  // channelPairs = {{0, 4}, {1, 5}, {2, 6}, {3, 7}}
  float sumFddAmplitude(const auto& charges)
  {
    float amp = 0;
    for (int i = 0; i < NFddPairs; i++) {
      if (charges[i] > 0 && charges[i + FddPairOffset] > 0) {
        amp += 0.5 * (charges[i] + charges[i + FddPairOffset]);
      }
    }
    return amp;
  } // end sumFddAmplitude()

  //--------------------------------------------------------------------------------
  // book the histograms of a run
  void addBcHistos(int run)
  {
    if (bcTH1Pointers[Form("bc/%d/bcSel_H", run)]) {
      return;
    }

    const int nBunches = o2::constants::lhc::LHCMaxBunches;
    const double maxBunch = static_cast<double>(nBunches) - 0.5;
    const int nBinsTF = static_cast<int>(nTF / tfPerBin) + 1;
    const double maxTF = static_cast<double>(nBinsTF * tfPerBin - 1) - 0.5;
    const AxisSpec axisTF{nBinsTF, -0.5, maxTF, "TF"};
    const AxisSpec axisBcb{nbcB, -0.5, nbcB - 0.5, "bc-B idx"};

    // filling scheme and bookkeeping
    bcTH1Pointers[Form("bc/%d/bcPatternB_H", run)] = bcTH1Registry.add<TH1>(Form("bc/%d/bcPatternB_H", run), "colliding bunches; bcID;",
                                                                            {HistType::kTH1D, {{nBunches, -0.5, maxBunch}}});
    bcTH1Pointers[Form("bc/%d/bcSel_H", run)] = bcTH1Registry.add<TH1>(Form("bc/%d/bcSel_H", run), "bc selection: all, no TF border, no ROF border, colliding; ; Counter",
                                                                       {HistType::kTH1D, {{NBinsBcSel, -0.5, NBinsBcSel - 0.5}}});
    bcTH1Pointers[Form("bc/%d/tf_H", run)] = bcTH1Registry.add<TH1>(Form("bc/%d/tf_H", run), "analysed time frames; TF; Counts",
                                                                    {HistType::kTH1D, {axisTF}});
    // stored bcs where the border cut computed here agrees / disagrees with the event selection
    bcTH1Pointers[Form("bc/%d/borderCheck_H", run)] = bcTH1Registry.add<TH1>(Form("bc/%d/borderCheck_H", run), "border cut vs event selection: agree, disagree; ; Counter",
                                                                             {HistType::kTH1D, {{NBinsBorderCheck, -0.5, NBinsBorderCheck - 0.5}}});

    // thresholds, stored as bin labels so that they survive merging
    bcTH1Pointers[Form("bc/%d/thresholds_H", run)] = bcTH1Registry.add<TH1>(Form("bc/%d/thresholds_H", run), "amplitude thresholds; threshold index; ",
                                                                            {HistType::kTH1D, {{nThresholds, -0.5, nThresholds - 0.5}}});
    for (int i = 0; i < nThresholds; i++) {
      bcTH1Pointers[Form("bc/%d/thresholds_H", run)]->GetXaxis()->SetBinLabel(i + 1, Form("%g", thresholds[i]));
      bcTH1Pointers[Form("bc/%d/thresholds_H", run)]->SetBinContent(i + 1, 1);
    }
    for (const auto& bcInOrbit : bcbList) {
      bcTH1Pointers[Form("bc/%d/bcPatternB_H", run)]->Fill(bcInOrbit);
    }

    // analysed crossings, stored or not
    bcTH2Pointers[Form("bc/%d/nCross_bcb_H", run)] = bcTH2Registry.add<TH2>(Form("bc/%d/nCross_bcb_H", run), "analysed crossings; TF; bc-B idx; Counter",
                                                                            {HistType::kTH2D, {axisTF, axisBcb}});

    // triggers from the CTP input mask, to compute mu
    bcTH2Pointers[Form("bc/%d/ft0Vtx_bcb_H", run)] = bcTH2Registry.add<TH2>(Form("bc/%d/ft0Vtx_bcb_H", run), "ft0Vtx triggers; TF; bc-B idx; Counter",
                                                                            {HistType::kTH2D, {axisTF, axisBcb}});
    bcTH2Pointers[Form("bc/%d/ft0VtxCe_bcb_H", run)] = bcTH2Registry.add<TH2>(Form("bc/%d/ft0VtxCe_bcb_H", run), "ft0VtxCe triggers; TF; bc-B idx; Counter",
                                                                              {HistType::kTH2D, {axisTF, axisBcb}});

    // pile-up probability. The denominators are filled with weights, so their
    // bin errors are not meaningful: compute the uncertainty of the ratio
    // binomially from the contents.
    bcTH2Pointers[Form("bc/%d/puDenA_bcb_H", run)] = bcTH2Registry.add<TH2>(Form("bc/%d/puDenA_bcb_H", run), "crossings empty on the C side; TF; bc-B idx; Counter",
                                                                            {HistType::kTH2D, {axisTF, axisBcb}});
    bcTH2Pointers[Form("bc/%d/puDenC_bcb_H", run)] = bcTH2Registry.add<TH2>(Form("bc/%d/puDenC_bcb_H", run), "crossings empty on the A side; TF; bc-B idx; Counter",
                                                                            {HistType::kTH2D, {axisTF, axisBcb}});
    for (int i = 0; i < nThresholds; i++) {
      bcTH2Pointers[Form("bc/%d/puNumA_t%d_bcb_H", run, i)] = bcTH2Registry.add<TH2>(Form("bc/%d/puNumA_t%d_bcb_H", run, i), Form("FT0A amplitude > %g, C side empty; TF; bc-B idx; Counter", thresholds[i]),
                                                                                     {HistType::kTH2D, {axisTF, axisBcb}});
      bcTH2Pointers[Form("bc/%d/puNumC_t%d_bcb_H", run, i)] = bcTH2Registry.add<TH2>(Form("bc/%d/puNumC_t%d_bcb_H", run, i), Form("FT0C amplitude > %g, A side empty; TF; bc-B idx; Counter", thresholds[i]),
                                                                                     {HistType::kTH2D, {axisTF, axisBcb}});
      // variation asking in addition for a signal in ZNA or ZNC. Since the ZN of
      // the other side is part of the empty condition, this is in practice the
      // ZN of the same side.
      if (fillZnVariation) {
        bcTH2Pointers[Form("bc/%d/puNumAZn_t%d_bcb_H", run, i)] = bcTH2Registry.add<TH2>(Form("bc/%d/puNumAZn_t%d_bcb_H", run, i), Form("FT0A amplitude > %g and ZN, C side empty; TF; bc-B idx; Counter", thresholds[i]),
                                                                                         {HistType::kTH2D, {axisTF, axisBcb}});
        bcTH2Pointers[Form("bc/%d/puNumCZn_t%d_bcb_H", run, i)] = bcTH2Registry.add<TH2>(Form("bc/%d/puNumCZn_t%d_bcb_H", run, i), Form("FT0C amplitude > %g and ZN, A side empty; TF; bc-B idx; Counter", thresholds[i]),
                                                                                         {HistType::kTH2D, {axisTF, axisBcb}});
      }
    }
  } // end addBcHistos()

  //--------------------------------------------------------------------------------
  // cache the histograms of a run
  void setHistoPointers(int run)
  {
    hBcSel = bcTH1Pointers[Form("bc/%d/bcSel_H", run)];
    hTF = bcTH1Pointers[Form("bc/%d/tf_H", run)];
    hBorderCheck = bcTH1Pointers[Form("bc/%d/borderCheck_H", run)];
    hCross = bcTH2Pointers[Form("bc/%d/nCross_bcb_H", run)];
    hVtx = bcTH2Pointers[Form("bc/%d/ft0Vtx_bcb_H", run)];
    hVtxCe = bcTH2Pointers[Form("bc/%d/ft0VtxCe_bcb_H", run)];
    hDenA = bcTH2Pointers[Form("bc/%d/puDenA_bcb_H", run)];
    hDenC = bcTH2Pointers[Form("bc/%d/puDenC_bcb_H", run)];
    hNumA.assign(nThresholds, nullptr);
    hNumC.assign(nThresholds, nullptr);
    hNumAZn.assign(nThresholds, nullptr);
    hNumCZn.assign(nThresholds, nullptr);
    for (int i = 0; i < nThresholds; i++) {
      hNumA[i] = bcTH2Pointers[Form("bc/%d/puNumA_t%d_bcb_H", run, i)];
      hNumC[i] = bcTH2Pointers[Form("bc/%d/puNumC_t%d_bcb_H", run, i)];
      if (fillZnVariation) {
        hNumAZn[i] = bcTH2Pointers[Form("bc/%d/puNumAZn_t%d_bcb_H", run, i)];
        hNumCZn[i] = bcTH2Pointers[Form("bc/%d/puNumCZn_t%d_bcb_H", run, i)];
      }
    }
  } // end setHistoPointers()

  //--------------------------------------------------------------------------------
  void init(InitContext&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    thresholds = static_cast<std::vector<float>>(ampThresholds);
    nThresholds = static_cast<int>(thresholds.size());
  } // end init()

  //--------------------------------------------------------------------------------
  void processBCs(BCsTSsSels const& bcs,
                  aod::FV0As const&, aod::FT0s const&, aod::FDDs const&,
                  aod::Zdcs const&)
  {
    if (bcs.size() == 0) {
      return;
    }

    // new run
    const auto& bcAt0 = bcs.iteratorAt(0);
    if (runNumberBc != bcAt0.runNumber()) {
      runNumberBc = bcAt0.runNumber();
      getRunInfo(runNumberBc);
      getFillingScheme();
      getBorderInfo();
      addBcHistos(runNumberBc);
      setHistoPointers(runNumberBc);
      seenTFs.clear();
    }

    for (const auto& bc : bcs) {
      int64_t thisBC = getBcWithinOrbit(bc.globalBC());
      int64_t thisTF = getTimeFrame(bc.globalBC());

      // first time this time frame is seen: count it and add all its analysed
      // crossings. The set protects against a time frame split over two
      // dataframes being counted twice.
      if (seenTFs.insert(thisTF).second) {
        hTF->Fill(thisTF);
        fillCrossings(thisTF);
      }

      // check that the border cut computed here agrees with the event selection
      bool passesFlags = bc.selection_bit(aod::evsel::kNoTimeFrameBorder) && bc.selection_bit(aod::evsel::kNoITSROFrameBorder);
      hBorderCheck->Fill(passesFlags == passesBorders(bc.globalBC()) ? 0 : 1);

      if (!checkBcFlags(bc)) {
        continue;
      }
      if (!bcPatternB.test(thisBC)) {
        continue;
      }
      hBcSel->Fill(3);
      int thisBcb = bcbIdx[thisBC];

      // triggers from the CTP input mask
      std::bitset<64> mask = bc.inputMask();
      if (mask[ft0VtxIdx]) {
        hVtx->Fill(thisTF, thisBcb);
        if (mask[ft0CeIdx]) {
          hVtxCe->Fill(thisTF, thisBcb);
        }
      }

      //------------------------------------------------------------------------------
      // detector information; amplitudes are only meaningful when in time

      // FT0
      bool ft0aInTime = false;
      bool ft0cInTime = false;
      float ampFT0A = 0;
      float ampFT0C = 0;
      if (bc.has_foundFT0()) {
        if (bc.foundFT0().isValidTimeA() && std::abs(bc.foundFT0().timeA()) < maxAbsTimeFT0) {
          ft0aInTime = true;
          ampFT0A = bc.foundFT0().sumAmpA();
        }
        if (bc.foundFT0().isValidTimeC() && std::abs(bc.foundFT0().timeC()) < maxAbsTimeFT0) {
          ft0cInTime = true;
          ampFT0C = bc.foundFT0().sumAmpC();
        }
      }

      // FV0A
      bool fv0aSignal = false;
      if (bc.has_foundFV0() && std::abs(bc.foundFV0().time()) < maxAbsTimeFV0) {
        auto amps = bc.foundFV0().amplitude();
        fv0aSignal = (std::accumulate(amps.begin(), amps.end(), 0.f) > 0);
      }

      // FDD
      bool fddaSignal = false;
      bool fddcSignal = false;
      if (bc.has_foundFDD()) {
        if (std::abs(bc.foundFDD().timeA()) < maxAbsTimeFDD) {
          fddaSignal = (sumFddAmplitude(bc.foundFDD().chargeA()) > 0);
        }
        if (std::abs(bc.foundFDD().timeC()) < maxAbsTimeFDD) {
          fddcSignal = (sumFddAmplitude(bc.foundFDD().chargeC()) > 0);
        }
      }

      // ZN: an entry in time, no requirement on the energy
      bool znaSignal = false;
      bool zncSignal = false;
      if (bc.has_zdc()) {
        float tZNA = bc.zdc().timeZNA();
        float tZNC = bc.zdc().timeZNC();
        znaSignal = std::isfinite(tZNA) && std::abs(tZNA) < maxAbsTimeZN;
        zncSignal = std::isfinite(tZNC) && std::abs(tZNC) < maxAbsTimeZN;
      }
      bool znSignal = znaSignal || zncSignal;

      // a side is empty if no detector has a signal in time, at any amplitude above zero
      bool ft0aSignal = ft0aInTime && (ampFT0A > 0);
      bool ft0cSignal = ft0cInTime && (ampFT0C > 0);
      bool noSignalC = !ft0cSignal && !fddcSignal && !zncSignal;
      bool noSignalA = !ft0aSignal && !fv0aSignal && !fddaSignal && !znaSignal;

      //------------------------------------------------------------------------------
      // pile-up counters. The denominators already hold all analysed crossings
      // of this bunch; a stored crossing with activity on the other side is
      // removed from them.
      if (noSignalC) {
        for (int i = 0; i < nThresholds; i++) {
          if (ft0aInTime && ampFT0A > thresholds[i]) {
            hNumA[i]->Fill(thisTF, thisBcb);
            if (fillZnVariation && znSignal) {
              hNumAZn[i]->Fill(thisTF, thisBcb);
            }
          }
        }
      } else {
        hDenA->Fill(thisTF, thisBcb, -1.);
      }
      if (noSignalA) {
        for (int i = 0; i < nThresholds; i++) {
          if (ft0cInTime && ampFT0C > thresholds[i]) {
            hNumC[i]->Fill(thisTF, thisBcb);
            if (fillZnVariation && znSignal) {
              hNumCZn[i]->Fill(thisTF, thisBcb);
            }
          }
        }
      } else {
        hDenC->Fill(thisTF, thisBcb, -1.);
      }
    } // loop over bcs
  } // end processBCs()
  PROCESS_SWITCH(UpcPileUp, processBCs, "pile-up probability in FT0", true);

}; // end of struct UpcPileUp

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcPileUp>(cfgc)};
}
