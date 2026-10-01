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
///
/// \brief Read output table from ZDC light ion task
/// \author chiara.oppedisano@cern.ch
//

#include "Common/DataModel/ZDCLightIons.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TH2.h>

#include <cstdio>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

#define CHECK_BIT(var, pos) (((var) >> (pos)) & 1)

struct ZDCLIAnalysis {

  // Configurable
  Configurable<bool> selectBC{"selectBC", 0, "Select BC events"};
  Configurable<bool> selectOnlyB{"selectOnlyB", 0, "Select BC with A && C"};
  Configurable<bool> selectColl{"selectColl", 0, "Select ALICE collision events"};
  //
  Configurable<uint64_t> tStampOffset{"tStampOffset", 0, "offset value for timestamp"};
  Configurable<int> nBinstStamp{"nBinstStamp", 1000, "no. bins in histo vs. timestamp"};
  Configurable<float> tStampMax{"tStampMax", 1000, ", maximum value for timestamp"};
  //
  Configurable<bool> tdcCut{"tdcCut", false, "Flag for TDC cut"};
  Configurable<float> tdcZNmincut{"tdcZNmincut", -1.5, "Min. ZN TDC cut value"};
  Configurable<float> tdcZNmaxcut{"tdcZNmaxcut", 1.5, "Max. ZN TDC cut value"};
  Configurable<float> tdcZPmincut{"tdcZPmincut", -1.5, "Min. ZP TDC cut value"};
  Configurable<float> tdcZPmaxcut{"tdcZPmaxcut", 1.5, "Max. ZP TDC cut value"};
  //
  Configurable<int> nBinsADC{"nBinsADC", 1000, "n bins 4 ZDC ADCs"};
  Configurable<int> nBinsAmpZN{"nBinsAmpZN", 1025, "n bins 4 ZN amplitudes"};
  Configurable<int> nBinsAmpZP{"nBinsAmpZP", 1025, "n bins 4 ZP amplitudes"};
  Configurable<int> nBinsTDC{"nBinsTDC", 480, "n bins 4 TDCs"};
  Configurable<int> nBinsFit{"nBinsFit", 1000, "n bins 4 FIT"};
  Configurable<float> MaxZN{"MaxZN", 4099.5, "Max 4 ZN histos"};
  Configurable<float> MaxZP{"MaxZP", 3099.5, "Max 4 ZP histos"};
  Configurable<float> MaxZEM{"MaxZEM", 3099.5, "Max 4 ZEM histos"};
  //
  Configurable<float> MaxMultFV0{"MaxMultFV0", 3000, "Max 4 FV0 histos"};
  Configurable<float> MaxMultFT0{"MaxMultFT0", 3000, "Max 4 FT0 histos"};
  //
  Configurable<float> enCalibZNA{"enCalibZNA", 1.0, "Energy calibration ZNA"};
  Configurable<float> enCalibZNC{"enCalibZNC", 1.0, "Energy calibration ZNC"};
  Configurable<float> enCalibZPA{"enCalibZPA", 1.0, "Energy calibration ZPA"};
  Configurable<float> enCalibZPC{"enCalibZPC", 1.0, "Energy calibration ZPC"};
  //
  Configurable<bool> applyZDCcut{"applyZDCcut", false, "Apply ZDC cut for light ion analysis"};
  Configurable<float> zdcaCutLow{"zdcaCutLow", 12., "ZDCA cut for light ion analysis"};
  Configurable<float> zdccCutLow{"zdccCutLow", 12., "ZDCC cut for light ion analysis"};
  //
  Configurable<bool> selectZvtx{"selectZvtx", true, "Activate Z vertex selection"};
  Configurable<bool> sel8{"sel8", true, "Activate sel8 selection"};
  Configurable<bool> triggerTVX{"triggerTVX", true, "Activate trigger TVX selection"};
  Configurable<bool> doOccupancySel{"doOccupancySel", false, "Activate occupancy selection"};
  Configurable<bool> noSameBunchPileupCut{"noSameBunchPileupCut", true, "Activate no same bunch pileup selection"};
  Configurable<bool> isGoodZvtxFT0vsPV{"isGoodZvtxFT0vsPV", true, "Activate is good Z vertex FT0 vs PV selection"};
  Configurable<bool> noCollInTimeRangeStandard{"noCollInTimeRangeStandard", true, "Activate no collision in time range standard selection"};
  Configurable<bool> noTimeFrameBorder{"noTimeFrameBorder", true, "Activate no time frame border selection"};
  Configurable<bool> noITSROFFrameBorder{"noITSROFFrameBorder", true, "Activate no ITS ROF frame border selection"};
  Configurable<bool> isGoodITSLayersAll{"isGoodITSLayersAll", false, "Activate is good ITS layers all selection"};
  //
  HistogramRegistry registry{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    registry.add("hBCmask", "mask; BC; counts", {HistType::kTH1F, {{5, 0., 5.}}});
    registry.add("hcounts", "counts; selections; counts", {HistType::kTH1F, {{14, 0., 14.}}});
    registry.add("hzvertex", "z vertex; z_vertex (cm); Entries", {HistType::kTH1F, {{200, -20., 20.}}});
    //
    registry.add("hZNApmc", "ZNA pmc; ZNA amplitude; Entries", {HistType::kTH1F, {{nBinsAmpZN, -0.5, MaxZN}}});
    registry.add("hZPApmc", "ZPA pmc; ZPA amplitude; Entries", {HistType::kTH1F, {{nBinsAmpZP, -0.5, MaxZP}}});
    registry.add("hZNCpmc", "ZNC pmc; ZNC amplitude; Entries", {HistType::kTH1F, {{nBinsAmpZN, -0.5, MaxZN}}});
    registry.add("hZPCpmc", "ZPC pmc; ZPC amplitude; Entries", {HistType::kTH1F, {{nBinsAmpZP, -0.5, MaxZP}}});
    registry.add("hZEM", "ZEM; ZEM1+ZEM2 amplitude; Entries", {HistType::kTH1F, {{nBinsAmpZP, -0.5, MaxZEM}}});
    registry.add("hZDCA", "ZNA+ZPA; ZNA+ZPA; Entries", {HistType::kTH1F, {{nBinsAmpZN + nBinsAmpZP, -0.5, MaxZN + MaxZP}}});
    registry.add("hZDCC", "ZNC+ZPC; ZNC+ZPC; Entries", {HistType::kTH1F, {{nBinsAmpZN + nBinsAmpZP, -0.5, MaxZN + MaxZP}}});
    //
    registry.add("hZDCApmcwZDCCcut", "ZDCA w. ZDCC cut; ZDCA amplitude; Entries", {HistType::kTH1F, {{nBinsAmpZN + nBinsAmpZP, -0.5, MaxZN + MaxZP}}});
    registry.add("hZDCCpmcwZDCAcut", "ZDCC w. ZDCA cut; ZDCC amplitude; Entries", {HistType::kTH1F, {{nBinsAmpZN + nBinsAmpZP, -0.5, MaxZN + MaxZP}}});
    //
    registry.add("hZNAtdc", "ZNA tdc; ZNA tdc; Entries", {HistType::kTH1F, {{520, -14., 12.}}});
    registry.add("hZPAtdc", "ZPA tdc; ZPA tdc; Entries", {HistType::kTH1F, {{520, -14., 12.}}});
    registry.add("hZNCtdc", "ZNC tdc; ZNC tdc; Entries", {HistType::kTH1F, {{520, -14., 12.}}});
    registry.add("hZPCtdc", "ZPC tdc; ZPC tdc; Entries", {HistType::kTH1F, {{520, -14., 12.}}});
    //
    registry.add("hZDCCvsA", "ZDC side C vs. side A; ZDCA; ZDCC", {HistType::kTH2F, {{{nBinsAmpZN + nBinsAmpZP, -0.5, MaxZN + MaxZP}, {nBinsAmpZN + nBinsAmpZP, -0.5, MaxZN + MaxZP}}}});
    //
    registry.add("hZNAamplvsADC", "ZNA amplitude vs. ADC; ZNA ADC; ZNA amplitude", {HistType::kTH2F, {{{nBinsAmpZN, -0.5, 3. * MaxZN}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZNCamplvsADC", "ZNC amplitude vs. ADC; ZNC ADC; ZNC amplitude", {HistType::kTH2F, {{{nBinsAmpZN, -0.5, 3. * MaxZN}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZPAamplvsADC", "ZPA amplitude vs. ADC; ZPA ADC; ZPA amplitude", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, 3. * MaxZP}, {nBinsAmpZP, -0.5, MaxZP}}}});
    registry.add("hZPCamplvsADC", "ZPC amplitude vs. ADC; ZPC ADC; ZPC amplitude", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, 3. * MaxZP}, {nBinsAmpZP, -0.5, MaxZP}}}});
    //
    registry.add("hZNvsZEM", "ZN vs ZEM; ZEM; ZNA+ZNC", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, MaxZEM}, {nBinsAmpZN, -0.5, 2. * MaxZN}}}});
    registry.add("hZPvsZEM", "ZP vs ZEM; ZEM; ZPA+ZPC", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, MaxZEM}, {nBinsAmpZP, -0.5, 2. * MaxZP}}}});
    registry.add("hZNAvsZEM", "ZNA vs ZEM; ZEM; ZNA", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, MaxZEM}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZNCvsZEM", "ZNC vs ZEM; ZEM; ZNC", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, MaxZEM}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZPAvsZEM", "ZPA vs ZEM; ZEM; ZPA", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, MaxZEM}, {nBinsAmpZP, -0.5, MaxZP}}}});
    registry.add("hZPCvsZEM", "ZPC vs ZEM; ZEM; ZPC", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, MaxZEM}, {nBinsAmpZP, -0.5, MaxZP}}}});
    //
    registry.add("hZNAvsZNC", "ZNA vs ZNC; ZNC; ZNA", {HistType::kTH2F, {{{nBinsAmpZN, -0.5, MaxZN}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZPAvsZPC", "ZPA vs ZPC; ZPC; ZPA", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, MaxZP}, {nBinsAmpZP, -0.5, MaxZP}}}});
    registry.add("hZNAvsZPA", "ZNA vs ZPA; ZPA; ZNA", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, MaxZP}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZNCvsZPC", "ZNC vs ZPC; ZPC; ZNC", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, MaxZP}, {nBinsAmpZN, -0.5, MaxZN}}}});
    //
    registry.add("hZNCcvsZNCsum", "ZNC PMC vs PMsum; ZNCC ADC; ZNC sum", {HistType::kTH2F, {{{nBinsADC, -0.5, 3. * MaxZN}, {nBinsADC, -0.5, 3. * MaxZN}}}});
    registry.add("hZNAcvsZNAsum", "ZNA PMC vs PMsum; ZNAC ADC; ZNA sum", {HistType::kTH2F, {{{nBinsADC, -0.5, 3. * MaxZN}, {nBinsADC, -0.5, 3. * MaxZN}}}});
    //
    registry.add("hZNCvstdc", "ZNC vs tdc; ZNC TDC (ns); ZNC amplitude", {HistType::kTH2F, {{{480, -13.5, 11.45}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZNAvstdc", "ZNA vs tdc; ZNA TDC (ns); ZNA amplitude", {HistType::kTH2F, {{{480, -13.5, 11.45}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZPCvstdc", "ZPC vs tdc; ZPC TDC (ns); ZPC amplitude", {HistType::kTH2F, {{{480, -13.5, 11.45}, {nBinsAmpZP, -0.5, MaxZP}}}});
    registry.add("hZPAvstdc", "ZPA vs tdc; ZPA TDC (ns); ZPA amplitude", {HistType::kTH2F, {{{480, -13.5, 11.45}, {nBinsAmpZP, -0.5, MaxZP}}}});
    //
    registry.add("hZNvsV0A", "ZN vs V0A", {HistType::kTH2F, {{{nBinsFit, 0., MaxMultFV0}, {nBinsAmpZN, -0.5, 2. * MaxZN}}}});
    registry.add("hZNAvsFT0A", "ZNA vs FT0A", {HistType::kTH2F, {{{nBinsFit, 0., MaxMultFT0}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZNCvsFT0C", "ZNC vs FT0C", {HistType::kTH2F, {{{nBinsFit, 0., MaxMultFT0}, {nBinsAmpZN, -0.5, MaxZN}}}});
    //
    registry.add("hZNAvscentrFT0A", "ZNA vs centrality FT0A", {HistType::kTH2F, {{{100, 0., 100.}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZNAvscentrFT0C", "ZNA vs centrality FT0C", {HistType::kTH2F, {{{100, 0., 100.}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZNAvscentrFT0M", "ZNA vs centrality FT0M", {HistType::kTH2F, {{{100, 0., 100.}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZPAvscentrFT0A", "ZPA vs centrality FT0A", {HistType::kTH2F, {{{100, 0., 100.}, {nBinsAmpZP, -0.5, MaxZP}}}});
    registry.add("hZPAvscentrFT0C", "ZPA vs centrality FT0C", {HistType::kTH2F, {{{100, 0., 100.}, {nBinsAmpZP, -0.5, MaxZP}}}});
    registry.add("hZPAvscentrFT0M", "ZPA vs centrality FT0M", {HistType::kTH2F, {{{100, 0., 100.}, {nBinsAmpZP, -0.5, MaxZP}}}});
    registry.add("hZNCvscentrFT0A", "ZNC vs centrality FT0A", {HistType::kTH2F, {{{100, 0., 100.}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZNCvscentrFT0C", "ZNC vs centrality FT0C", {HistType::kTH2F, {{{100, 0., 100.}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZNCvscentrFT0M", "ZNC vs centrality FT0M", {HistType::kTH2F, {{{100, 0., 100.}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZPCvscentrFT0A", "ZPC vs centrality FT0A", {HistType::kTH2F, {{{100, 0., 100.}, {nBinsAmpZP, -0.5, MaxZP}}}});
    registry.add("hZPCvscentrFT0C", "ZPC vs centrality FT0C", {HistType::kTH2F, {{{100, 0., 100.}, {nBinsAmpZP, -0.5, MaxZP}}}});
    registry.add("hZPCvscentrFT0M", "ZPC vs centrality FT0M", {HistType::kTH2F, {{{100, 0., 100.}, {nBinsAmpZP, -0.5, MaxZP}}}});
    //
    registry.add("hNnZNAvscentrFT0C", "N_{neutrons} in ZNA vs FT0C", {HistType::kTH2F, {{{100, 0., 100.}, {12, -0.5, 11.5}}}});
    registry.add("hNnZNCvscentrFT0C", "N_{neutrons} in ZNC vs FT0C", {HistType::kTH2F, {{{100, 0., 100.}, {12, -0.5, 11.5}}}});
    registry.add("hNpZPAvscentrFT0C", "N_{protons} in ZPA vs FT0C", {HistType::kTH2F, {{{100, 0., 100.}, {12, -0.5, 11.5}}}});
    registry.add("hNpZPCvscentrFT0C", "N_{protons} in ZPC vs FT0C", {HistType::kTH2F, {{{100, 0., 100.}, {12, -0.5, 11.5}}}});
    registry.add("hNnZNAvscentrFT0M", "N_{neutrons} in ZNA vs FT0M", {HistType::kTH2F, {{{100, 0., 100.}, {12, -0.5, 11.5}}}});
    registry.add("hNnZNCvscentrFT0M", "N_{neutrons} in ZNC vs FT0M", {HistType::kTH2F, {{{100, 0., 100.}, {12, -0.5, 11.5}}}});
    registry.add("hNpZPAvscentrFT0M", "N_{protons} in ZPA vs FT0M", {HistType::kTH2F, {{{100, 0., 100.}, {12, -0.5, 11.5}}}});
    registry.add("hNpZPCvscentrFT0M", "N_{protons} in ZPC vs FT0M", {HistType::kTH2F, {{{100, 0., 100.}, {12, -0.5, 11.5}}}});
    //
    registry.add("hNpvsNnZNA", "N_{protons} vs N_{neutrons} in ZNA", {HistType::kTH2F, {{{12, -0.5, 11.5}, {12, -0.5, 11.5}}}});
    registry.add("hNpvsNnZNC", "N_{protons} vs N_{neutrons} in ZNC", {HistType::kTH2F, {{{12, -0.5, 11.5}, {12, -0.5, 11.5}}}});
    //
    registry.add("hZNAvstimestamp", "ZNA vs timestamp", {HistType::kTH2F, {{{nBinstStamp, 0., tStampMax}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZNCvstimestamp", "ZNC vs timestamp", {HistType::kTH2F, {{{nBinstStamp, 0., tStampMax}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZPAvstimestamp", "ZPA vs timestamp", {HistType::kTH2F, {{{nBinstStamp, 0., tStampMax}, {nBinsAmpZP, -0.5, MaxZP}}}});
    registry.add("hZPCvstimestamp", "ZPC vs timestamp", {HistType::kTH2F, {{{nBinstStamp, 0., tStampMax}, {nBinsAmpZP, -0.5, MaxZP}}}});
  }

  void process(aod::ZDCLightIons const& zdclightions)
  {
    for (auto const& zdc : zdclightions) {
      auto tdczna = zdc.znaTdc();
      auto tdcznc = zdc.zncTdc();
      auto tdczpa = zdc.zpaTdc();
      auto tdczpc = zdc.zpcTdc();
      auto tdczem1 = zdc.zem1Tdc();
      auto tdczem2 = zdc.zem2Tdc();
      auto zna = zdc.znaAmpl();
      auto znaADC = zdc.znaPmc();
      auto znapm1 = zdc.znaPm1();
      auto znapm2 = zdc.znaPm2();
      auto znapm3 = zdc.znaPm3();
      auto znapm4 = zdc.znaPm4();
      auto znc = zdc.zncAmpl();
      auto zncADC = zdc.zncPmc();
      auto zncpm1 = zdc.zncPm1();
      auto zncpm2 = zdc.zncPm2();
      auto zncpm3 = zdc.zncPm3();
      auto zncpm4 = zdc.zncPm4();
      auto zpa = zdc.zpaAmpl();
      auto zpaADC = zdc.zpaPmc();
      auto zpc = zdc.zpcAmpl();
      auto zpcADC = zdc.zpcPmc();
      auto zem1 = zdc.zem1Ampl();
      auto zem2 = zdc.zem2Ampl();
      auto multFT0A = zdc.multFt0a();
      auto multFT0C = zdc.multFt0c();
      auto multV0A = zdc.multV0a();
      auto zvtx = zdc.vertexZ();
      auto centrFT0C = zdc.centralityFt0c();
      auto centrFT0A = zdc.centralityFt0a();
      auto centrFT0M = zdc.centralityFt0m();
      auto timestamp = zdc.timestamp();
      auto selectionBits = zdc.selectionBits();
      auto bcMask = zdc.bcMask();

      bool isZNAtdc = false;
      bool isZNCtdc = false;
      bool isZPAtdc = false;
      bool isZPCtdc = false;
      bool isZEMtdc = false;

      if (tdczna > -99999.0) {
        isZNAtdc = true;
      }
      if (tdcznc > -99999.0) {
        isZNCtdc = true;
      }
      if (tdczpa > -99999.0) {
        isZPAtdc = true;
      }
      if (tdczpc > -99999.0) {
        isZPCtdc = true;
      }
      if (tdczem1 > -99999.0 || tdczem2 > -99999.0) {
        isZEMtdc = true;
      }

      if (tdcCut) { // TDC cuts applied
        if ((tdczna < tdcZNmincut) || (tdczna > tdcZNmaxcut))
          isZNAtdc = false;
        if ((tdcznc < tdcZNmincut) || (tdcznc > tdcZNmaxcut))
          isZNCtdc = false;
        if ((tdczpa < tdcZPmincut) || (tdczpa > tdcZPmaxcut))
          isZPAtdc = false;
        if ((tdczpc < tdcZPmincut) || (tdczpc > tdcZPmaxcut))
          isZPCtdc = false;
      }

      bool eventSelected = false;

      // for BC events -------
      if (selectBC) {
        // printf("BCmask: %x \n", bcMask);

        registry.get<TH1>(HIST("hBCmask"))->Fill(0., 1.);
        auto isB = CHECK_BIT(bcMask, 0);
        if (isB)
          registry.get<TH1>(HIST("hBCmask"))->Fill(1., 1.);
        if (CHECK_BIT(bcMask, 1))
          registry.get<TH1>(HIST("hBCmask"))->Fill(2., 1.);
        if (CHECK_BIT(bcMask, 2))
          registry.get<TH1>(HIST("hBCmask"))->Fill(3., 1.);
        if (CHECK_BIT(bcMask, 3))
          registry.get<TH1>(HIST("hBCmask"))->Fill(4., 1.);
        //
        if (selectOnlyB && isB)
          eventSelected = true;
        else if (!selectOnlyB)
          eventSelected = true;
      }

      // for collision events -------
      if (selectColl) {
        bool zvtxSel = false;
        if (selectZvtx && CHECK_BIT(selectionBits, 0))
          zvtxSel = true;
        else if (!selectZvtx)
          zvtxSel = true;
        //
        bool ottoSel = false;
        if (sel8 && CHECK_BIT(selectionBits, 1))
          ottoSel = true;
        else if (!sel8)
          ottoSel = true;
        //
        bool isdoOccupancySel = false;
        if (doOccupancySel && CHECK_BIT(selectionBits, 2))
          isdoOccupancySel = true;
        else if (!doOccupancySel)
          isdoOccupancySel = true;
        //
        bool isnoSameBunchPileupCut = false;
        if (noSameBunchPileupCut && CHECK_BIT(selectionBits, 3))
          isnoSameBunchPileupCut = true;
        else if (!noSameBunchPileupCut)
          isnoSameBunchPileupCut = true;
        //
        bool isGoodZvtxFT0vsPVsel = false;
        if (isGoodZvtxFT0vsPV && CHECK_BIT(selectionBits, 4))
          isGoodZvtxFT0vsPVsel = true;
        else if (!isGoodZvtxFT0vsPV)
          isGoodZvtxFT0vsPVsel = true;
        //
        bool isnoCollInTimeRangeStandard = false;
        if (noCollInTimeRangeStandard && CHECK_BIT(selectionBits, 5))
          isnoCollInTimeRangeStandard = true;
        else if (!noCollInTimeRangeStandard)
          isnoCollInTimeRangeStandard = true;
        //
        bool isnoTimeFrameBorder = false;
        if (noTimeFrameBorder && CHECK_BIT(selectionBits, 6))
          isnoTimeFrameBorder = true;
        else if (!noTimeFrameBorder)
          isnoTimeFrameBorder = true;
        //
        bool isnoITSROFFrameBorder = false;
        if (noITSROFFrameBorder && CHECK_BIT(selectionBits, 7))
          isnoITSROFFrameBorder = true;
        else if (!noITSROFFrameBorder)
          isnoITSROFFrameBorder = true;
        //
        bool isGoodITSLayersAllsel = false;
        if (isGoodITSLayersAll && CHECK_BIT(selectionBits, 8))
          isGoodITSLayersAllsel = true;
        else if (!isGoodITSLayersAll)
          isGoodITSLayersAllsel = true;
        //
        bool istriggerTVX = false;
        if (triggerTVX && CHECK_BIT(selectionBits, 9))
          istriggerTVX = true;
        else if (!triggerTVX)
          istriggerTVX = true;

        if (zvtxSel && ottoSel && istriggerTVX && isdoOccupancySel && isnoSameBunchPileupCut && isGoodZvtxFT0vsPVsel && isnoCollInTimeRangeStandard && isnoTimeFrameBorder && isnoITSROFFrameBorder && isGoodITSLayersAllsel)
          eventSelected = true;
        // if (zvtxSel && ottoSel && isnoSameBunchPileupCut) eventSelected = true;

        if (eventSelected) {
          registry.get<TH1>(HIST("hcounts"))->Fill(0., 1.);
          if (isZNAtdc)
            registry.get<TH1>(HIST("hcounts"))->Fill(1., 1.);
          if (isZPAtdc)
            registry.get<TH1>(HIST("hcounts"))->Fill(2., 1.);
          if (isZNCtdc)
            registry.get<TH1>(HIST("hcounts"))->Fill(3., 1.);
          if (isZPCtdc)
            registry.get<TH1>(HIST("hcounts"))->Fill(4., 1.);
          if (isZNAtdc || isZNCtdc)
            registry.get<TH1>(HIST("hcounts"))->Fill(5., 1.);
          if (isZNAtdc || isZPAtdc)
            registry.get<TH1>(HIST("hcounts"))->Fill(6., 1.);
          if (isZNCtdc || isZPCtdc)
            registry.get<TH1>(HIST("hcounts"))->Fill(7., 1.);
          if (isZNAtdc || isZPAtdc || isZNCtdc || isZPCtdc)
            registry.get<TH1>(HIST("hcounts"))->Fill(8., 1.);
          if (isZNAtdc && isZNCtdc)
            registry.get<TH1>(HIST("hcounts"))->Fill(9., 1.);
          if (isZNAtdc && isZPAtdc)
            registry.get<TH1>(HIST("hcounts"))->Fill(10., 1.);
          if (isZNCtdc && isZPCtdc)
            registry.get<TH1>(HIST("hcounts"))->Fill(11., 1.);
          if (isZEMtdc)
            registry.get<TH1>(HIST("hcounts"))->Fill(12., 1.);
        }
      }

      if (eventSelected) {

        if (enCalibZNA > 0.) {
          zna *= enCalibZNA;
          znaADC *= enCalibZNA;
          znapm1 *= enCalibZNA;
          znapm2 *= enCalibZNA;
          znapm3 *= enCalibZNA;
          znapm4 *= enCalibZNA;
        }
        if (enCalibZNC > 0.) {
          znc *= enCalibZNC;
          zncADC *= enCalibZNC;
          zncpm1 *= enCalibZNC;
          zncpm2 *= enCalibZNC;
          zncpm3 *= enCalibZNC;
          zncpm4 *= enCalibZNC;
        }
        if (enCalibZPA > 0.) {
          zpa *= enCalibZPA;
          zpaADC *= enCalibZPA;
        }
        if (enCalibZPC > 0.) {
          zpc *= enCalibZPC;
          zpcADC *= enCalibZPC;
        }

        registry.get<TH1>(HIST("hzvertex"))->Fill(zvtx);

        if (!applyZDCcut) {

          if (isZNAtdc)
            registry.get<TH1>(HIST("hZNApmc"))->Fill(zna);
          if (isZNCtdc)
            registry.get<TH1>(HIST("hZNCpmc"))->Fill(znc);
          if (isZPAtdc)
            registry.get<TH1>(HIST("hZPApmc"))->Fill(zpa);
          if (isZPCtdc)
            registry.get<TH1>(HIST("hZPCpmc"))->Fill(zpc);
          //
          if (isZNAtdc || isZPAtdc)
            registry.get<TH1>(HIST("hZDCA"))->Fill(zna + zpa);
          if (isZNCtdc || isZPCtdc)
            registry.get<TH1>(HIST("hZDCC"))->Fill(znc + zpc);
          //
          if (isZNAtdc)
            registry.get<TH2>(HIST("hZNAamplvsADC"))->Fill(znaADC, zna);
          if (isZNCtdc)
            registry.get<TH2>(HIST("hZNCamplvsADC"))->Fill(zncADC, znc);
          if (isZPAtdc)
            registry.get<TH2>(HIST("hZPAamplvsADC"))->Fill(zpaADC, zpa);
          if (isZPCtdc)
            registry.get<TH2>(HIST("hZPCamplvsADC"))->Fill(zpcADC, zpc);
          //
          if (isZNAtdc || isZNCtdc)
            registry.get<TH2>(HIST("hZNAvsZNC"))->Fill(znc, zna);
          if (isZPAtdc || isZPCtdc)
            registry.get<TH2>(HIST("hZPAvsZPC"))->Fill(zpc, zpa);
          if (isZNAtdc || isZPAtdc)
            registry.get<TH2>(HIST("hZNAvsZPA"))->Fill(zpa, zna);
          if (isZNCtdc || isZPCtdc)
            registry.get<TH2>(HIST("hZNCvsZPC"))->Fill(zpc, znc);
          //
          if (isZNAtdc)
            registry.get<TH2>(HIST("hZNAvstdc"))->Fill(tdczna, zna);
          if (isZNCtdc)
            registry.get<TH2>(HIST("hZNCvstdc"))->Fill(tdcznc, znc);
          if (isZPAtdc)
            registry.get<TH2>(HIST("hZPAvstdc"))->Fill(tdczpa, zpa);
          if (isZPCtdc)
            registry.get<TH2>(HIST("hZPCvstdc"))->Fill(tdczpc, zpc);
          //
          if (isZNAtdc)
            registry.get<TH2>(HIST("hZNAcvsZNAsum"))->Fill(0.25 * (znapm1 + znapm2 + znapm3 + znapm4), zna);
          if (isZNCtdc)
            registry.get<TH2>(HIST("hZNCcvsZNCsum"))->Fill(0.25 * (zncpm1 + zncpm2 + zncpm3 + zncpm4), znc);
          //
          if (isZNAtdc || isZNCtdc)
            registry.get<TH2>(HIST("hZNvsV0A"))->Fill(multV0A / 100., zna + znc);
          if (isZNAtdc)
            registry.get<TH2>(HIST("hZNAvsFT0A"))->Fill((multFT0A) / 100., zna);
          if (isZNCtdc)
            registry.get<TH2>(HIST("hZNCvsFT0C"))->Fill((multFT0C) / 100., znc);
        } else {
          bool isZNChigh = false;
          if (znc >= zdccCutLow)
            isZNChigh = true;
          bool isZNAhigh = false;
          if (zna >= zdcaCutLow)
            isZNAhigh = true;
          bool isZPAhigh = false;
          if (zpa >= zdcaCutLow)
            isZPAhigh = true;
          bool isZPChigh = false;
          if (zpc >= zdccCutLow)
            isZPChigh = true;
          //
          if (isZNAtdc && isZNChigh)
            registry.get<TH1>(HIST("hZNApmc"))->Fill(zna);
          if (isZNCtdc && isZNAhigh)
            registry.get<TH1>(HIST("hZNCpmc"))->Fill(znc);
          if (isZPAtdc && isZPChigh)
            registry.get<TH1>(HIST("hZPApmc"))->Fill(zpa);
          if (isZPCtdc && isZPAhigh)
            registry.get<TH1>(HIST("hZPCpmc"))->Fill(zpc);
          //
          if ((isZNAtdc || isZPAtdc) && (isZNChigh && isZPChigh))
            registry.get<TH1>(HIST("hZDCA"))->Fill(zna + zpa);
          if ((isZNCtdc || isZPCtdc) && (isZNAhigh && isZPAhigh))
            registry.get<TH1>(HIST("hZDCC"))->Fill(znc + zpc);
          //
          if (isZNAtdc && isZNChigh)
            registry.get<TH2>(HIST("hZNAamplvsADC"))->Fill(znaADC, zna);
          if (isZNCtdc && isZNAhigh)
            registry.get<TH2>(HIST("hZNCamplvsADC"))->Fill(zncADC, znc);
          if (isZPAtdc && isZPChigh)
            registry.get<TH2>(HIST("hZPAamplvsADC"))->Fill(zpaADC, zpa);
          if (isZPCtdc && isZPAhigh)
            registry.get<TH2>(HIST("hZPCamplvsADC"))->Fill(zpcADC, zpc);
          //
          if (isZNAtdc && isZNChigh)
            registry.get<TH2>(HIST("hZNAvstdc"))->Fill(tdczna, zna);
          if (isZNCtdc && isZNAhigh)
            registry.get<TH2>(HIST("hZNCvstdc"))->Fill(tdcznc, znc);
          if (isZPAtdc && isZPChigh)
            registry.get<TH2>(HIST("hZPAvstdc"))->Fill(tdczpa, zpa);
          if (isZPCtdc && isZPAhigh)
            registry.get<TH2>(HIST("hZPCvstdc"))->Fill(tdczpc, zpc);
          //
          if (isZNAtdc && isZNChigh)
            registry.get<TH2>(HIST("hZNAcvsZNAsum"))->Fill(0.25 * (znapm1 + znapm2 + znapm3 + znapm4), zna);
          if (isZNCtdc && isZNAhigh)
            registry.get<TH2>(HIST("hZNCcvsZNCsum"))->Fill(0.25 * (zncpm1 + zncpm2 + zncpm3 + zncpm4), znc);
          //
          if (isZNAtdc && isZNChigh)
            registry.get<TH2>(HIST("hZNAvsFT0A"))->Fill((multFT0A) / 100., zna);
          if (isZNCtdc && isZNAhigh)
            registry.get<TH2>(HIST("hZNCvsFT0C"))->Fill((multFT0C) / 100., znc);
        }

        // Timestamp
        /*if (tStampOffset > timestamp) {
          printf("\n\n #################  OFFSET timestamp too large!!!!!!!!!!!!!!!!!!!!!!!!!! >  timestamp %lu \n\n", timestamp);
          return;
        }*/
        // float tsh = (timestamp / 1000.) - (tStampOffset / 1000.); // in hours
        /*if (timestamp > tStampMax) {
          printf("\n\n MAXIMUM timestamp too small!!!!!!!!!!!!!!!!!!!!!!!!!! > timestamp-offset %f \n\n", timestamp);
          return;
        }*/

        if (!applyZDCcut) {
          if (isZNAtdc)
            registry.get<TH2>(HIST("hZNAvstimestamp"))->Fill(timestamp, zna);
          if (isZNCtdc)
            registry.get<TH2>(HIST("hZNCvstimestamp"))->Fill(timestamp, znc);
          if (isZPAtdc)
            registry.get<TH2>(HIST("hZPAvstimestamp"))->Fill(timestamp, zpa);
          if (isZPCtdc)
            registry.get<TH2>(HIST("hZPCvstimestamp"))->Fill(timestamp, zpc);
        } else {
          if (isZNAtdc && znc >= zdccCutLow) {
            registry.get<TH2>(HIST("hZNAvstimestamp"))->Fill(timestamp, zna);
          }
          if (isZPAtdc && zpc >= zdccCutLow) {
            registry.get<TH2>(HIST("hZPAvstimestamp"))->Fill(timestamp, zpa);
          }
          if ((isZNAtdc || isZPAtdc) && (znc >= zdccCutLow && zpc >= zdccCutLow)) {
            registry.get<TH1>(HIST("hZDCApmcwZDCCcut"))->Fill(zna + zpa);
          }
          if (isZNCtdc && zna >= zdcaCutLow) {
            registry.get<TH2>(HIST("hZNCvstimestamp"))->Fill(timestamp, znc);
          }
          if (isZPCtdc && zpa >= zdcaCutLow) {
            registry.get<TH2>(HIST("hZPCvstimestamp"))->Fill(timestamp, zpc);
          }
          if ((isZNCtdc || isZPCtdc) && (zna >= zdcaCutLow && zpa >= zdcaCutLow)) {
            registry.get<TH1>(HIST("hZDCCpmcwZDCAcut"))->Fill(znc + zpc);
          }
        }
        // TDCs
        if (isZNAtdc)
          registry.get<TH1>(HIST("hZNAtdc"))->Fill(tdczna);
        if (isZNCtdc)
          registry.get<TH1>(HIST("hZNCtdc"))->Fill(tdcznc);
        if (isZPAtdc)
          registry.get<TH1>(HIST("hZPAtdc"))->Fill(tdczpa);
        if (isZPCtdc)
          registry.get<TH1>(HIST("hZPCtdc"))->Fill(tdczpc);
        //
        if (isZEMtdc) {
          registry.get<TH1>(HIST("hZEM"))->Fill(zem1 + zem2);
          registry.get<TH2>(HIST("hZNAvsZEM"))->Fill(zem1 + zem2, zna);
          registry.get<TH2>(HIST("hZNCvsZEM"))->Fill(zem1 + zem2, znc);
          registry.get<TH2>(HIST("hZPAvsZEM"))->Fill(zem1 + zem2, zpa);
          registry.get<TH2>(HIST("hZPCvsZEM"))->Fill(zem1 + zem2, zpc);
        }
        if (isZNAtdc || isZNCtdc)
          registry.get<TH2>(HIST("hZNvsZEM"))->Fill(zem1 + zem2, zna + znc);
        if (isZPAtdc || isZPCtdc)
          registry.get<TH2>(HIST("hZPvsZEM"))->Fill(zem1 + zem2, zpa + zpc);
        //
        if (isZNAtdc || isZNCtdc || isZPAtdc || isZPCtdc)
          registry.get<TH2>(HIST("hZDCCvsA"))->Fill(zna + zpa, znc + zpc);
        //
        if (centrFT0C > -1. && centrFT0C < 101.) {
          if (isZNAtdc)
            registry.get<TH2>(HIST("hZNAvscentrFT0C"))->Fill(centrFT0C, zna);
          if (isZPAtdc)
            registry.get<TH2>(HIST("hZPAvscentrFT0C"))->Fill(centrFT0C, zpa);
          if (isZNCtdc)
            registry.get<TH2>(HIST("hZNCvscentrFT0C"))->Fill(centrFT0C, znc);
          if (isZPCtdc)
            registry.get<TH2>(HIST("hZPCvscentrFT0C"))->Fill(centrFT0C, zpc);
          registry.get<TH2>(HIST("hNnZNAvscentrFT0C"))->Fill(centrFT0C, zna / 2.680);
          registry.get<TH2>(HIST("hNnZNCvscentrFT0C"))->Fill(centrFT0C, znc / 2.680);
          registry.get<TH2>(HIST("hNpZPAvscentrFT0C"))->Fill(centrFT0C, zpa / 2.680);
          registry.get<TH2>(HIST("hNpZPCvscentrFT0C"))->Fill(centrFT0C, zpc / 2.680);
        }
        if (centrFT0A > -1. && centrFT0A < 101.) {
          if (isZNAtdc)
            registry.get<TH2>(HIST("hZNAvscentrFT0A"))->Fill(centrFT0A, zna);
          if (isZPAtdc)
            registry.get<TH2>(HIST("hZPAvscentrFT0A"))->Fill(centrFT0A, zpa);
          if (isZNCtdc)
            registry.get<TH2>(HIST("hZNCvscentrFT0A"))->Fill(centrFT0A, znc);
          if (isZPCtdc)
            registry.get<TH2>(HIST("hZPCvscentrFT0A"))->Fill(centrFT0A, zpc);
        }
        if (centrFT0M > -1. && centrFT0M < 101.) {
          if (isZNAtdc)
            registry.get<TH2>(HIST("hZNAvscentrFT0M"))->Fill(centrFT0M, zna);
          if (isZPAtdc)
            registry.get<TH2>(HIST("hZPAvscentrFT0M"))->Fill(centrFT0M, zpa);
          if (isZNCtdc)
            registry.get<TH2>(HIST("hZNCvscentrFT0M"))->Fill(centrFT0M, znc);
          if (isZPCtdc)
            registry.get<TH2>(HIST("hZPCvscentrFT0M"))->Fill(centrFT0M, zpc);
          registry.get<TH2>(HIST("hNnZNAvscentrFT0M"))->Fill(centrFT0M, zna / 2.680);
          registry.get<TH2>(HIST("hNnZNCvscentrFT0M"))->Fill(centrFT0M, znc / 2.680);
          registry.get<TH2>(HIST("hNpZPAvscentrFT0M"))->Fill(centrFT0M, zpa / 2.680);
          registry.get<TH2>(HIST("hNpZPCvscentrFT0M"))->Fill(centrFT0M, zpc / 2.680);
          registry.get<TH2>(HIST("hNpvsNnZNA"))->Fill(zna / 2.680, zpa / 2.680);
          registry.get<TH2>(HIST("hNpvsNnZNC"))->Fill(znc / 2.680, zpc / 2.680);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ZDCLIAnalysis>(cfgc) //
  };
}
