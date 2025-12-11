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
#include <Framework/HistogramRegistry.h>
#include <Framework/runDataProcessing.h>

#include <TH1F.h>
#include <TH2F.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

struct ZDCLIAnalysis {

  // Configurable number of bins
  Configurable<bool> useZvtx{"useZvtx", false, "If true uses Z_vertex"};
  Configurable<float> zVval{"zVval", 10., "Z_vertex cut value"};
  Configurable<uint64_t> tStampOffset{"tStampOffset", 0, "offset value for timestamp"};
  Configurable<int> nBinstStamp{"nBinstStamp", 1000, "no. bins in histo vs. timestamp"};
  Configurable<float> tStampMax{"tStampMax", 1000, ",maximum value for timestamp"};
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
  HistogramRegistry registry{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext const&)
  {
    registry.add("hZNApmc", "ZNApmc; ZNA amplitude; Entries", {HistType::kTH1F, {{nBinsAmpZN, -0.5, MaxZN}}});
    registry.add("hZPApmc", "ZPApmc; ZPA amplitude; Entries", {HistType::kTH1F, {{nBinsAmpZP, -0.5, MaxZP}}});
    registry.add("hZNCpmc", "ZNCpmc; ZNC amplitude; Entries", {HistType::kTH1F, {{nBinsAmpZN, -0.5, MaxZN}}});
    registry.add("hZPCpmc", "ZPCpmc; ZPC amplitude; Entries", {HistType::kTH1F, {{nBinsAmpZP, -0.5, MaxZP}}});
    registry.add("hZEM", "ZEM; ZEM1+ZEM2 amplitude; Entries", {HistType::kTH1F, {{nBinsAmpZP, -0.5, MaxZEM}}});
    registry.add("hZNAamplvsADC", "ZNA amplitude vs. ADC; ZNA ADC; ZNA amplitude", {HistType::kTH2F, {{{nBinsAmpZN, -0.5, 3. * MaxZN}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZNCamplvsADC", "ZNC amplitude vs. ADC; ZNC ADC; ZNC amplitude", {HistType::kTH2F, {{{nBinsAmpZN, -0.5, 3. * MaxZN}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZPAamplvsADC", "ZPA amplitude vs. ADC; ZPA ADC; ZPA amplitude", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, 3. * MaxZP}, {nBinsAmpZP, -0.5, MaxZP}}}});
    registry.add("hZPCamplvsADC", "ZPC amplitude vs. ADC; ZPC ADC; ZPC amplitude", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, 3. * MaxZP}, {nBinsAmpZP, -0.5, MaxZP}}}});
    registry.add("hZNvsZEM", "ZN vs ZEM; ZEM; ZNA+ZNC", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, MaxZEM}, {nBinsAmpZN, -0.5, 2. * MaxZN}}}});
    registry.add("hZNAvsZNC", "ZNA vs ZNC; ZNC; ZNA", {HistType::kTH2F, {{{nBinsAmpZN, -0.5, MaxZN}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZPAvsZPC", "ZPA vs ZPC; ZPC; ZPA", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, MaxZP}, {nBinsAmpZP, -0.5, MaxZP}}}});
    registry.add("hZNAvsZPA", "ZNA vs ZPA; ZPA; ZNA", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, MaxZP}, {nBinsAmpZN, -0.5, MaxZN}}}});
    registry.add("hZNCvsZPC", "ZNC vs ZPC; ZPC; ZNC", {HistType::kTH2F, {{{nBinsAmpZP, -0.5, MaxZP}, {nBinsAmpZN, -0.5, MaxZN}}}});
    //
    registry.add("hZNCcvsZNCsum", "ZNC PMC vs PMsum; ZNCC ADC; ZNCsum", {HistType::kTH2F, {{{nBinsADC, -0.5, 3. * MaxZN}, {nBinsADC, -0.5, 3. * MaxZN}}}});
    registry.add("hZNAcvsZNAsum", "ZNA PMC vs PMsum; ZNAsum", {HistType::kTH2F, {{{nBinsADC, -0.5, 3. * MaxZN}, {nBinsADC, -0.5, 3. * MaxZN}}}});
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
      // auto selectionBits = zdc.selectionBits();

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

      registry.get<TH1>(HIST("hZNApmc"))->Fill(zna);
      registry.get<TH1>(HIST("hZNCpmc"))->Fill(znc);
      registry.get<TH1>(HIST("hZPApmc"))->Fill(zpa);
      registry.get<TH1>(HIST("hZPCpmc"))->Fill(zpc);
      registry.get<TH1>(HIST("hZEM"))->Fill(zem1 + zem2);
      //
      registry.get<TH2>(HIST("hZNAamplvsADC"))->Fill(znaADC, zna);
      registry.get<TH2>(HIST("hZNCamplvsADC"))->Fill(zncADC, znc);
      registry.get<TH2>(HIST("hZPAamplvsADC"))->Fill(zpaADC, zpa);
      registry.get<TH2>(HIST("hZPCamplvsADC"))->Fill(zpcADC, zpc);
      //
      registry.get<TH2>(HIST("hZNvsZEM"))->Fill(zem1 + zem2, zna + znc);
      registry.get<TH2>(HIST("hZNAvsZNC"))->Fill(znc, zna);
      registry.get<TH2>(HIST("hZPAvsZPC"))->Fill(zpc, zpa);
      registry.get<TH2>(HIST("hZNAvsZPA"))->Fill(zpa, zna);
      registry.get<TH2>(HIST("hZNCvsZPC"))->Fill(zpc, znc);
      //
      registry.get<TH2>(HIST("hZNAvstdc"))->Fill(tdczna, zna);
      registry.get<TH2>(HIST("hZNCvstdc"))->Fill(tdcznc, znc);
      registry.get<TH2>(HIST("hZPAvstdc"))->Fill(tdczpa, zpa);
      registry.get<TH2>(HIST("hZPCvstdc"))->Fill(tdczpc, zpc);
      //
      registry.get<TH2>(HIST("hZNAcvsZNAsum"))->Fill(0.25 * (znapm1 + znapm2 + znapm3 + znapm4), zna);
      registry.get<TH2>(HIST("hZNCcvsZNCsum"))->Fill(0.25 * (zncpm1 + zncpm2 + zncpm3 + zncpm4), znc);
      //
      registry.get<TH2>(HIST("hZNvsV0A"))->Fill(multV0A / 100., zna + znc);
      registry.get<TH2>(HIST("hZNAvsFT0A"))->Fill((multFT0A) / 100., zna);
      registry.get<TH2>(HIST("hZNCvsFT0C"))->Fill((multFT0C) / 100., znc);
      //
      registry.get<TH2>(HIST("hZNAvscentrFT0A"))->Fill(centrFT0A, zna);
      registry.get<TH2>(HIST("hZNAvscentrFT0C"))->Fill(centrFT0C, zna);
      registry.get<TH2>(HIST("hZNAvscentrFT0M"))->Fill(centrFT0M, zna);
      registry.get<TH2>(HIST("hZPAvscentrFT0A"))->Fill(centrFT0A, zpa);
      registry.get<TH2>(HIST("hZPAvscentrFT0C"))->Fill(centrFT0C, zpa);
      registry.get<TH2>(HIST("hZPAvscentrFT0M"))->Fill(centrFT0M, zpa);
      registry.get<TH2>(HIST("hZNCvscentrFT0A"))->Fill(centrFT0A, znc);
      registry.get<TH2>(HIST("hZNCvscentrFT0C"))->Fill(centrFT0C, znc);
      registry.get<TH2>(HIST("hZNCvscentrFT0M"))->Fill(centrFT0M, znc);
      registry.get<TH2>(HIST("hZPCvscentrFT0A"))->Fill(centrFT0A, zpc);
      registry.get<TH2>(HIST("hZPCvscentrFT0C"))->Fill(centrFT0C, zpc);
      registry.get<TH2>(HIST("hZPCvscentrFT0M"))->Fill(centrFT0M, zpc);
      //
      /*if (tStampOffset > timestamp) {
        printf("\n\n #################  OFFSET timestamp too large!!!!!!!!!!!!!!!!!!!!!!!!!! >  timestamp %llu \n\n", timestamp);
        return;
      }*/
      float tsh = (timestamp / 1000.) - (tStampOffset / 1000.); // in hours
      /*if (tsh > tStampMax) {
        printf("\n\n MAXIMUM timestamp too small!!!!!!!!!!!!!!!!!!!!!!!!!! > timestamp-offset %f \n\n", tsh);
        return;
      }*/
      registry.get<TH2>(HIST("hZNAvstimestamp"))->Fill(tsh, zna);
      registry.get<TH2>(HIST("hZNCvstimestamp"))->Fill(tsh, znc);
      registry.get<TH2>(HIST("hZPAvstimestamp"))->Fill(tsh, zpa);
      registry.get<TH2>(HIST("hZPCvstimestamp"))->Fill(tsh, zpc);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ZDCLIAnalysis>(cfgc) //
  };
}
