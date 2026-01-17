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

/// \file eventSelectionQa.cxx
/// \brief Event selection QA task
///
/// \author Evgeny Kryshen <evgeny.kryshen@cern.ch> and Igor Altsybeev <Igor.Altsybeev@cern.ch>

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "Common/DataModel/EventSelection.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonDataFormat/BunchFilling.h>
#include <CommonDataFormat/TimeStamp.h>
#include <DataFormatsITSMFT/TimeDeadMap.h>
#include <DataFormatsParameters/AggregatedRunInfo.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ITSMFTBase/DPLAlpideParam.h>
#include <ITSMFTReconstruction/ChipMappingITS.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <TH1.h>
#include <TMath.h>
#include <TString.h>

#include <sys/types.h>

#include <bitset>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

using namespace o2::framework;
using namespace o2;
using namespace o2::aod::evsel;

using BCsRun2 = soa::Join<aod::BCs, aod::Run2BCInfos, aod::Timestamps, aod::BcSels, aod::Run2MatchedToBCSparse>;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using ColEvSels = soa::Join<aod::Collisions, aod::EvSels>;
using FullTracksIU = soa::Join<aod::TracksIU, aod::TracksExtra>;
using FullTracksIUwithLabels = soa::Join<aod::TracksIU, aod::TracksExtra, aod::McTrackLabels>;

struct EventSelectionQaTask {
  Configurable<bool> isMC{"isMC", 0, "0 - data, 1 - MC"};
  Configurable<int32_t> nGlobalBCs{"nGlobalBCs", 100000, "number of global bcs for detailed monitoring"};
  Configurable<bool> isLowFlux{"isLowFlux", 1, "1 - low flux (pp, pPb), 0 - high flux (PbPb)"};
  Configurable<bool> fillITSdeadStaveHists{"fillITSdeadStaveHists", 0, "0 - no, 1 - yes"};
  Configurable<bool> fillTPCnClsVsOccupancyHists{"fillTPCnClsVsOccupancyHists", 0, "0 - no, 1 - yes"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  static const int32_t nBCsPerOrbit = o2::constants::lhc::LHCMaxBunches;
  int32_t lastRun = -1;
  int64_t nOrbits = 1;                             // number of orbits, setting 1 for unanchored MC
  int64_t orbitSOR = 0;                            // first orbit, setting 0 for unanchored MC
  int64_t bcSOR = 0;                               // global bc of the start of the first orbit, setting 0 for unanchored MC
  int32_t nOrbitsPerTF = 128;                      // 128 in 2022, 32 in 2023, setting 128 for unanchored MC
  int64_t nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit; // duration of TF in bcs
  int rofOffset = -1;                              // ITS ROF offset, in bc
  int rofLength = -1;                              // ITS ROF length, in bc

  std::bitset<nBCsPerOrbit> bcPatternA;
  std::bitset<nBCsPerOrbit> bcPatternC;
  std::bitset<nBCsPerOrbit> bcPatternB;
  SliceCache cache;
  Partition<aod::Tracks> tracklets = (aod::track::trackType == static_cast<uint8_t>(o2::aod::track::TrackTypeEnum::Run2Tracklet));

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

  void init(InitContext&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    const AxisSpec axisMultV0M{1000, 0., isLowFlux ? 40000. : 40000., "V0M multiplicity"};
    const AxisSpec axisMultV0A{1000, 0., isLowFlux ? 40000. : 200000., "V0A multiplicity"};
    const AxisSpec axisMultV0C{1000, 0., isLowFlux ? 30000. : 30000., "V0C multiplicity"};
    const AxisSpec axisMultT0A{1000, 0., isLowFlux ? 10000. : 200000., "T0A multiplicity"};
    const AxisSpec axisMultT0C{1000, 0., isLowFlux ? 2000. : 70000., "T0C multiplicity"};
    const AxisSpec axisMultT0M{1000, 0., isLowFlux ? 12000. : 270000., "T0M multiplicity"};
    const AxisSpec axisMultFDA{1000, 0., isLowFlux ? 50000. : 40000., "FDA multiplicity"};
    const AxisSpec axisMultFDC{1000, 0., isLowFlux ? 50000. : 40000., "FDC multiplicity"};
    const AxisSpec axisMultZNA{1000, 0., isLowFlux ? 1000. : 400., "ZNA multiplicity"};
    const AxisSpec axisMultZNC{1000, 0., isLowFlux ? 1000. : 400., "ZNC multiplicity"};
    const AxisSpec axisNtracklets{200, 0., isLowFlux ? 200. : 6000., "n tracklets"};
    const AxisSpec axisNclusters{200, 0., isLowFlux ? 1000. : 20000., "n clusters"};
    const AxisSpec axisMultOnlineV0M{400, 0., isLowFlux ? 8000. : 40000., "Online V0M"};
    const AxisSpec axisMultOnlineFOR{300, 0., isLowFlux ? 300. : 1200., "Online FOR"};
    const AxisSpec axisMultOflineFOR{300, 0., isLowFlux ? 300. : 1200., "Ofline FOR"};

    const AxisSpec axisTime{700, -35., 35., ""};
    const AxisSpec axisTimeDif{100, -10., 10., ""};
    const AxisSpec axisTimeSum{100, -10., 10., ""};
    const AxisSpec axisGlobalBCs{nGlobalBCs, 0., static_cast<double>(nGlobalBCs), ""};
    const AxisSpec axisBCs{nBCsPerOrbit, 0., static_cast<double>(nBCsPerOrbit), ""};
    const AxisSpec axisNcontrib{200, 0., isLowFlux ? 200. : 8000., "n contributors"};
    const AxisSpec axisEta{100, -1., 1., "track #eta"};
    const AxisSpec axisColTimeRes{1500, 0., 1500., "collision time resolution (ns)"};
    const AxisSpec axisBcDif{600, -300., 300., "collision bc difference"};
    const AxisSpec axisAliases{kNaliases, 0., static_cast<double>(kNaliases), ""};
    const AxisSpec axisSelections{kNsel, 0., static_cast<double>(kNsel), ""};
    const AxisSpec axisVtxZ{500, -25., 25., ""};
    const AxisSpec axisVtxXY{500, -1., 1., ""};

    histos.add("hTimeV0Aall", "All bcs;V0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeV0Call", "All bcs;V0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNAall", "All bcs;ZNA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNCall", "All bcs;ZNC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Aall", "All bcs;T0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Call", "All bcs;T0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDAall", "All bcs;FDA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDCall", "All bcs;FDC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZACall", "All bcs; ZNC-ZNA time (ns); ZNC+ZNA time (ns)", kTH2F, {axisTimeDif, axisTimeSum});
    histos.add("hTimeV0Abga", "BeamA-only bcs;V0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeV0Cbga", "BeamA-only bcs;V0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNAbga", "BeamA-only bcs;ZNA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNCbga", "BeamA-only bcs;ZNC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Abga", "BeamA-only bcs;T0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Cbga", "BeamA-only bcs;T0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDAbga", "BeamA-only bcs;FDA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDCbga", "BeamA-only bcs;FDC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeV0Abgc", "BeamC-only bcs;V0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeV0Cbgc", "BeamC-only bcs;V0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNAbgc", "BeamC-only bcs;ZNA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNCbgc", "BeamC-only bcs;ZNC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Abgc", "BeamC-only bcs;T0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Cbgc", "BeamC-only bcs;T0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDAbgc", "BeamC-only bcs;FDA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDCbgc", "BeamC-only bcs;FDC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeV0Aref", "Reference bcs;V0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeV0Cref", "Reference bcs;V0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNAref", "Reference bcs;ZNA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNCref", "Reference bcs;ZNC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Aref", "Reference bcs;T0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Cref", "Reference bcs;T0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDAref", "Reference bcs;FDA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDCref", "Reference bcs;FDC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZACref", "Reference bcs; ZNC-ZNA time (ns); ZNC+ZNA time (ns)", kTH2F, {axisTimeDif, axisTimeSum});
    histos.add("hTimeV0Acol", "All events;V0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeV0Ccol", "All events;V0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNAcol", "All events;ZNA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNCcol", "All events;ZNC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Acol", "All events;T0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Ccol", "All events;T0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDAcol", "All events;FDA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDCcol", "All events;FDC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZACcol", "All events; ZNC-ZNA time (ns); ZNC+ZNA time (ns)", kTH2F, {axisTimeDif, axisTimeSum});
    histos.add("hTimeV0Aacc", "Accepted events;V0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeV0Cacc", "Accepted events;V0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNAacc", "Accepted events;ZNA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZNCacc", "Accepted events;ZNC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Aacc", "Accepted events;T0A time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeT0Cacc", "Accepted events;T0C time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDAacc", "Accepted events;FDA time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeFDCacc", "Accepted events;FDC time (ns);Entries", kTH1F, {axisTime});
    histos.add("hTimeZACacc", "Accepted events; ZNC-ZNA time (ns); ZNC+ZNA time (ns)", kTH2F, {axisTimeDif, axisTimeSum});
    histos.add("hSPDClsVsTklCol", "All events", kTH2F, {axisNtracklets, axisNclusters});
    histos.add("hV0C012vsTklCol", "All events;n tracklets;V0C012 multiplicity", kTH2F, {axisNtracklets, axisMultV0C});
    histos.add("hV0MOnVsOfCol", "All events", kTH2F, {axisMultV0M, axisMultOnlineV0M});
    histos.add("hSPDOnVsOfCol", "All events", kTH2F, {axisMultOflineFOR, axisMultOnlineFOR});
    histos.add("hV0C3vs012Col", "All events;V0C012 multiplicity;V0C3 multiplicity", kTH2F, {axisMultV0C, axisMultV0C});
    histos.add("hSPDClsVsTklAcc", "Accepted events", kTH2F, {axisNtracklets, axisNclusters});
    histos.add("hV0C012vsTklAcc", "Accepted events;n tracklets;V0C012 multiplicity", kTH2F, {axisNtracklets, axisMultV0C});
    histos.add("hV0MOnVsOfAcc", "Accepted events", kTH2F, {axisMultV0M, axisMultOnlineV0M});
    histos.add("hSPDOnVsOfAcc", "Accepted events", kTH2F, {axisMultOflineFOR, axisMultOnlineFOR});
    histos.add("hV0C3vs012Acc", "Accepted events;V0C012 multiplicity;V0C3 multiplicity", kTH2F, {axisMultV0C, axisMultV0C});

    histos.add("hColCounterAll", "", kTH1F, {axisAliases});
    histos.add("hColCounterAcc", "", kTH1F, {axisAliases});
    histos.add("hBcCounterAll", "", kTH1F, {axisAliases});
    histos.add("hSelCounter", "", kTH1F, {axisSelections});
    histos.add("hSelMask", "", kTH1F, {axisSelections});

    histos.add("hGlobalBcAll", "", kTH1F, {axisGlobalBCs});
    histos.add("hGlobalBcCol", "", kTH1F, {axisGlobalBCs});
    histos.add("hGlobalBcFT0", "", kTH1F, {axisGlobalBCs});
    histos.add("hGlobalBcFV0", "", kTH1F, {axisGlobalBCs});
    histos.add("hGlobalBcFDD", "", kTH1F, {axisGlobalBCs});
    histos.add("hGlobalBcZDC", "", kTH1F, {axisGlobalBCs});

    histos.add("hBcA", "", kTH1F, {axisBCs});
    histos.add("hBcC", "", kTH1F, {axisBCs});
    histos.add("hBcB", "", kTH1F, {axisBCs});
    histos.add("hBcAll", "", kTH1F, {axisBCs});
    histos.add("hBcCol", "", kTH1F, {axisBCs});
    histos.add("hBcTVX", "", kTH1F, {axisBCs});
    histos.add("hBcFT0", "", kTH1F, {axisBCs});
    histos.add("hBcFV0", "", kTH1F, {axisBCs});
    histos.add("hBcFDD", "", kTH1F, {axisBCs});
    histos.add("hBcZDC", "", kTH1F, {axisBCs});
    histos.add("hBcColTOF", "", kTH1F, {axisBCs});
    histos.add("hBcColTRD", "", kTH1F, {axisBCs});
    histos.add("hBcTrackTOF", "", kTH1F, {axisBCs});
    histos.add("hBcTrackTRD", "", kTH1F, {axisBCs});

    histos.add("hMultV0Aall", "All bcs", kTH1F, {axisMultV0A});
    histos.add("hMultV0Call", "All bcs", kTH1F, {axisMultV0C});
    histos.add("hMultZNAall", "All bcs", kTH1F, {axisMultZNA});
    histos.add("hMultZNCall", "All bcs", kTH1F, {axisMultZNC});
    histos.add("hMultT0Aall", "All bcs", kTH1F, {axisMultT0A});
    histos.add("hMultT0Call", "All bcs", kTH1F, {axisMultT0C});
    histos.add("hMultFDAall", "All bcs", kTH1F, {axisMultFDA});
    histos.add("hMultFDCall", "All bcs", kTH1F, {axisMultFDC});
    histos.add("hMultV0Aref", "Reference bcs", kTH1F, {axisMultV0A});
    histos.add("hMultV0Cref", "Reference bcs", kTH1F, {axisMultV0C});
    histos.add("hMultZNAref", "Reference bcs", kTH1F, {axisMultZNA});
    histos.add("hMultZNCref", "Reference bcs", kTH1F, {axisMultZNC});
    histos.add("hMultT0Aref", "Reference bcs", kTH1F, {axisMultT0A});
    histos.add("hMultT0Cref", "Reference bcs", kTH1F, {axisMultT0C});
    histos.add("hMultFDAref", "Reference bcs", kTH1F, {axisMultFDA});
    histos.add("hMultFDCref", "Reference bcs", kTH1F, {axisMultFDC});
    histos.add("hMultV0Mcol", "All events", kTH1F, {axisMultV0M});
    histos.add("hMultV0Acol", "All events", kTH1F, {axisMultV0A});
    histos.add("hMultV0Ccol", "All events", kTH1F, {axisMultV0C});
    histos.add("hMultZNAcol", "All events", kTH1F, {axisMultZNA});
    histos.add("hMultZNCcol", "All events", kTH1F, {axisMultZNC});
    histos.add("hMultT0Acol", "All events", kTH1F, {axisMultT0A});
    histos.add("hMultT0Ccol", "All events", kTH1F, {axisMultT0C});
    histos.add("hMultFDAcol", "All events", kTH1F, {axisMultFDA});
    histos.add("hMultFDCcol", "All events", kTH1F, {axisMultFDC});
    histos.add("hMultV0Macc", "Accepted events", kTH1F, {axisMultV0M});
    histos.add("hMultV0Aacc", "Accepted events", kTH1F, {axisMultV0A});
    histos.add("hMultV0Cacc", "Accepted events", kTH1F, {axisMultV0C});
    histos.add("hMultZNAacc", "Accepted events", kTH1F, {axisMultZNA});
    histos.add("hMultZNCacc", "Accepted events", kTH1F, {axisMultZNC});
    histos.add("hMultT0Aacc", "Accepted events", kTH1F, {axisMultT0A});
    histos.add("hMultT0Cacc", "Accepted events", kTH1F, {axisMultT0C});
    histos.add("hMultFDAacc", "Accepted events", kTH1F, {axisMultFDA});
    histos.add("hMultFDCacc", "Accepted events", kTH1F, {axisMultFDC});

    histos.add("hMultT0Abga", "A-side beam-gas events", kTH1F, {axisMultT0A});
    histos.add("hMultT0Abgc", "C-side beam-gas events", kTH1F, {axisMultT0A});
    histos.add("hMultT0Cbga", "A-side beam-gas events", kTH1F, {axisMultT0C});
    histos.add("hMultT0Cbgc", "C-side beam-gas events", kTH1F, {axisMultT0C});

    histos.add("hMultT0Mall", "BCs with collisions", kTH1F, {axisMultT0M});
    histos.add("hMultT0Mref", "", kTH1F, {axisMultT0M});
    histos.add("hMultT0Mtvx", "", kTH1F, {axisMultT0M});
    histos.add("hMultT0Mzac", "", kTH1F, {axisMultT0M});
    histos.add("hMultT0Mpup", "BCs with pileup", kTH1F, {axisMultT0M});

    histos.add("hMultT0Atvx", "", kTH1F, {axisMultT0A});
    histos.add("hMultT0Ctvx", "", kTH1F, {axisMultT0C});
    histos.add("hMultT0Azac", "", kTH1F, {axisMultT0A});
    histos.add("hMultT0Czac", "", kTH1F, {axisMultT0C});

    histos.add("hColTimeResVsNcontrib", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("hColTimeResVsNcontribITSonly", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("hColTimeResVsNcontribWithTOF", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("hColTimeResVsNcontribWithTRD", "", kTH2F, {axisNcontrib, axisColTimeRes});
    histos.add("hColBcDiffVsNcontrib", "", kTH2F, {axisNcontrib, axisBcDif});
    histos.add("hColBcDiffVsNcontribITSonly", "", kTH2F, {axisNcontrib, axisBcDif});
    histos.add("hColBcDiffVsNcontribWithTOF", "", kTH2F, {axisNcontrib, axisBcDif});
    histos.add("hColBcDiffVsNcontribWithTRD", "", kTH2F, {axisNcontrib, axisBcDif});

    histos.add("hITStrackBcDiff", "", kTH1F, {axisBcDif});
    histos.add("hTrackBcDiffVsEta", "", kTH2F, {axisEta, axisBcDif});
    histos.add("hTrackBcDiffVsEtaAll", "", kTH2F, {axisEta, axisBcDif});

    histos.add("hNcontribCol", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribAcc", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribMis", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribColTOF", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribColTRD", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribAccTOF", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribAccTRD", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribMisTOF", "", kTH1F, {axisNcontrib});

    histos.add("hMultT0MVsNcontribTVX", "", kTH2F, {axisMultT0M, axisNcontrib});          // before ITS RO Frame border cut
    histos.add("hMultT0MVsNcontribTVXTFcuts", "", kTH2F, {axisMultT0M, axisNcontrib});    // before ITS RO Frame border cut
    histos.add("hMultT0MVsNcontribTVXROFcuts", "", kTH2F, {axisMultT0M, axisNcontrib});   // after ITS RO Frame border cut
    histos.add("hMultT0MVsNcontribTVXTFROFcuts", "", kTH2F, {axisMultT0M, axisNcontrib}); // after ITS RO Frame border cut
    // histos.add("hMultT0MVsNcontribAcc", "", kTH2F, {axisMultT0M, axisNcontrib});          // before ITS RO Frame border cut

    histos.add("hMultV0AVsNcontribTVX", "", kTH2F, {axisMultV0A, axisNcontrib});            // before ITS RO Frame border cut
    histos.add("hMultV0AVsNcontribTVXTFcuts", "", kTH2F, {axisMultV0A, axisNcontrib});      // before ITS RO Frame border cut
    histos.add("hMultV0AVsNcontribTVXROFcuts", "", kTH2F, {axisMultV0A, axisNcontrib});     // before ITS RO Frame border cut
    histos.add("hMultV0AVsNcontribTVXTFROFcuts", "", kTH2F, {axisMultV0A, axisNcontrib});   // after ITS RO Frame border cut
    histos.add("hMultV0AVsNcontribIsVertexITSTPC", "", kTH2F, {axisMultV0A, axisNcontrib}); // after good vertex cut
    histos.add("hMultV0AVsNcontribGood", "", kTH2F, {axisMultV0A, axisNcontrib});           // after pileup check

    // histos.add("hFoundBcForMultV0AVsNcontribAcc", "", kTH1F, {axisBCs});      // bc distribution for V0A-vs-Ncontrib accepted
    histos.add("hFoundBcForMultV0AVsNcontribOutliers", "", kTH1F, {axisBCs}); // bc distribution for V0A-vs-Ncontrib outliers
    histos.add("hFoundBcAfterROFborderCut", "", kTH1F, {axisBCs});            // bc distribution for V0A-vs-Ncontrib after ITS-ROF border cut

    histos.add("hVtxFT0VsVtxCol", "", kTH2F, {axisVtxZ, axisVtxZ});                // FT0-vertex vs z-vertex from collisions
    histos.add("hVtxFT0MinusVtxCol", "", kTH1F, {axisVtxZ});                       // FT0-vertex minus z-vertex from collisions
    histos.add("hVtxFT0MinusVtxColVsMultT0M", "", kTH2F, {axisVtxZ, axisMultT0M}); // FT0-vertex minus z-vertex from collisions vs multiplicity

    histos.add("hFoundBc", "", kTH1F, {axisBCs});            // distribution of found bcs (for ITS ROF studies)
    histos.add("hFoundBcTOF", "", kTH1F, {axisBCs});         // distribution of found bcs (TOF-matched vertex)
    histos.add("hFoundBcNcontrib", "", kTH1F, {axisBCs});    // accumulated distribution of n contributors vs found bc (for ITS ROF studies)
    histos.add("hFoundBcNcontribTOF", "", kTH1F, {axisBCs}); // accumulated distribution of n contributors vs found bc (TOF-matched vertex)

    // MC histograms
    histos.add("hGlobalBcColMC", "", kTH1F, {axisGlobalBCs});
    histos.add("hBcColMC", "", kTH1F, {axisBCs});
    histos.add("hVertexXMC", "", kTH1F, {axisVtxXY});
    histos.add("hVertexYMC", "", kTH1F, {axisVtxXY});
    histos.add("hVertexZMC", "", kTH1F, {axisVtxZ});
    histos.add("hNcontribColFromMC", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribAccFromMC", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribMisFromMC", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribColFromData", "", kTH1F, {axisNcontrib});
    histos.add("hNcontribAccFromData", "", kTH1F, {axisNcontrib});

    for (int i = 0; i < kNsel; i++) {
      histos.get<TH1>(HIST("hSelCounter"))->GetXaxis()->SetBinLabel(i + 1, selectionLabels[i]);
      histos.get<TH1>(HIST("hSelMask"))->GetXaxis()->SetBinLabel(i + 1, selectionLabels[i]);
    }
    for (int i = 0; i < kNaliases; i++) {
      histos.get<TH1>(HIST("hColCounterAll"))->GetXaxis()->SetBinLabel(i + 1, aliasLabels[i].data());
      histos.get<TH1>(HIST("hColCounterAcc"))->GetXaxis()->SetBinLabel(i + 1, aliasLabels[i].data());
      histos.get<TH1>(HIST("hBcCounterAll"))->GetXaxis()->SetBinLabel(i + 1, aliasLabels[i].data());
    }

    // ROF border QA
    histos.add("ITSROFborderQA/hFoundBC_kTVX_counter_ITSTPCtracks", "", kTH1D, {axisBCs});
    histos.add("ITSROFborderQA/hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks", "", kTH1D, {axisBCs});

    // occupancy QA
    if (!isLowFlux) {
      histos.add("occupancyQA/hOccupancyByTracks", "", kTH1D, {{15002, -1.5, 15000.5}});
      histos.add("occupancyQA/hOccupancyByFT0C", "", kTH1D, {{15002, -20, 150000}});
      histos.add("occupancyQA/hOccupancyByFT0CvsByTracks", "", kTH2D, {{150, 0, 15000}, {150, 0, 150000}});

      // 3D histograms: nGlobalTracks with cls567 as y-axis, V0A as x-axis:
      const AxisSpec axisNtracksPV{200, -0.5, 5000 - 0.5, "n ITS PV tracks"};
      const AxisSpec axisNtracksPVTPC{160, -0.5, 4000 - 0.5, "n ITS-TPC PV tracks"};
      const AxisSpec axisNtracksTPConly{160, -0.5, 8000 - 0.5, "n TPC-only tracks"};
      const AxisSpec axisMultV0AForOccup{20, 0., static_cast<float>(200000), "mult V0A"};
      const AxisSpec axisOccupancyTracks{150, 0., 15000, "occupancy (n ITS tracks weighted)"};
      histos.add("occupancyQA/hNumTracksPV_vs_V0A_vs_occupancy", "", kTH3F, {axisMultV0AForOccup, axisNtracksPV, axisOccupancyTracks});
      histos.add("occupancyQA/hNumTracksPVTPC_vs_V0A_vs_occupancy", "", kTH3F, {axisMultV0AForOccup, axisNtracksPVTPC, axisOccupancyTracks});
      histos.add("occupancyQA/hNumTracksPVTPCLooseCuts_vs_V0A_vs_occupancy", "", kTH3F, {axisMultV0AForOccup, axisNtracksPVTPC, axisOccupancyTracks});
      histos.add("occupancyQA/hNumTracksITS_vs_V0A_vs_occupancy", "", kTH3F, {axisMultV0AForOccup, axisNtracksPV, axisOccupancyTracks});
      histos.add("occupancyQA/hNumTracksITSTPC_vs_V0A_vs_occupancy", "", kTH3F, {axisMultV0AForOccup, axisNtracksPVTPC, axisOccupancyTracks});
      histos.add("occupancyQA/hNumTracksPV_vs_V0A_vs_occupancy_NarrowDeltaTimeCut", "", kTH3F, {axisMultV0AForOccup, axisNtracksPV, axisOccupancyTracks});
      histos.add("occupancyQA/hNumTracksPVTPC_vs_V0A_vs_occupancy_NarrowDeltaTimeCut", "", kTH3F, {axisMultV0AForOccup, axisNtracksPVTPC, axisOccupancyTracks});
      histos.add("occupancyQA/hNumTracksPV_vs_V0A_vs_occupancy_StandardDeltaTimeCut", "", kTH3F, {axisMultV0AForOccup, axisNtracksPV, axisOccupancyTracks});
      histos.add("occupancyQA/hNumTracksPVTPC_vs_V0A_vs_occupancy_StandardDeltaTimeCut", "", kTH3F, {axisMultV0AForOccup, axisNtracksPVTPC, axisOccupancyTracks});
      histos.add("occupancyQA/hNumTracksPV_vs_V0A_vs_occupancy_GoodITSLayersAllCut", "", kTH3F, {axisMultV0AForOccup, axisNtracksPV, axisOccupancyTracks});
      histos.add("occupancyQA/hNumTracksPVTPC_vs_V0A_vs_occupancy_GoodITSLayersAllCut", "", kTH3F, {axisMultV0AForOccup, axisNtracksPVTPC, axisOccupancyTracks});
      // requested by TPC experts: nTPConly tracks vs occupancy
      histos.add("occupancyQA/hNumTracksTPConly_vs_V0A_vs_occupancy", "", kTH3F, {axisMultV0AForOccup, axisNtracksTPConly, axisOccupancyTracks});
      histos.add("occupancyQA/hNumTracksTPConlyNoITS_vs_V0A_vs_occupancy", "", kTH3F, {axisMultV0AForOccup, axisNtracksTPConly, axisOccupancyTracks});
      // request from experts to add track properties vs occupancy, to compare data vs MC
      if (fillTPCnClsVsOccupancyHists) {
        const AxisSpec axisOccupancyForTrackQA{60, 0., 15000, "occupancy (n ITS tracks weighted)"};
        const AxisSpec axisNTPCcls{150, 0, 150, "n TPC clusters"};
        histos.add("occupancyQA/tpcNClsFound_vs_V0A_vs_occupancy", "", kTH3F, {axisMultV0AForOccup, axisNTPCcls, axisOccupancyForTrackQA});
        histos.add("occupancyQA/tpcNClsFindable_vs_V0A_vs_occupancy", "", kTH3F, {axisMultV0AForOccup, axisNTPCcls, axisOccupancyForTrackQA});
        histos.add("occupancyQA/tpcNClsShared_vs_V0A_vs_occupancy", "", kTH3F, {axisMultV0AForOccup, axisNTPCcls, axisOccupancyForTrackQA});
        histos.add("occupancyQA/tpcNCrossedRows_vs_V0A_vs_occupancy", "", kTH3F, {axisMultV0AForOccup, axisNTPCcls, axisOccupancyForTrackQA});
        const AxisSpec axisChi2TPC{150, 0, 15, "chi2Ncl TPC"};
        histos.add("occupancyQA/tpcChi2_vs_V0A_vs_occupancy", "", kTH3F, {axisMultV0AForOccup, axisChi2TPC, axisOccupancyForTrackQA});
      }

      // ITS in-ROF occupancy
      histos.add("occupancyQA/hITSTracks_ev1_vs_ev2_2coll_in_ROF", ";nITStracks event #1;nITStracks event #2", kTH2D, {{200, 0., 6000}, {200, 0., 6000}});
      histos.add("occupancyQA/hITSTracks_ev1_vs_ev2_2coll_in_ROF_UPC", ";nITStracks event #1;nITStracks event #2", kTH2D, {{41, -0.5, 40.5}, {41, -0.5, 40.5}});
      histos.add("occupancyQA/hITSTracks_ev1_vs_ev2_2coll_in_ROF_nonUPC", ";nITStracks event #1;nITStracks event #2", kTH2D, {{200, 0., 6000}, {200, 0., 6000}});

      histos.add("occupancyQA/dEdx_vs_centr_vs_occup_narrow_p_win", "dE/dx", kTH3F, {{20, 0, 4000, "n PV tracks"}, {60, 0, 15000, "occupancy"}, {800, 0.0, 800.0, "dE/dx (a. u.)"}});
    }
  }

  void processRun2(
    ColEvSels const& cols,
    BCsRun2 const& bcs,
    aod::Zdcs const&,
    aod::FV0As const&,
    aod::FV0Cs const&,
    aod::FT0s const&,
    aod::FDDs const&)
  {
    bool isINT1period = 0;

    int run = bcs.iteratorAt(0).runNumber();
    if (run != lastRun) {
      lastRun = run;
      auto firstBC = bcs.iteratorAt(0);
      EventSelectionParams* par = ccdb->getForTimeStamp<EventSelectionParams>("EventSelection/EventSelectionParams", firstBC.timestamp());
      bool* applySelection = par->getSelection(0);
      for (int i = 0; i < kNsel; i++) {
        histos.get<TH1>(HIST("hSelMask"))->SetBinContent(i + 1, applySelection[i]);
      }
      isINT1period = run <= 136377 || (run >= 144871 && run <= 159582);
    }

    // bc-based event selection qa
    for (const auto& bc : bcs) {
      for (int iAlias = 0; iAlias < kNaliases; iAlias++) {
        histos.fill(HIST("hBcCounterAll"), iAlias, bc.alias_bit(iAlias));
      }
    }

    // collision-based event selection qa
    for (const auto& col : cols) {
      bool sel1 = col.selection_bit(kIsINT1) && col.selection_bit(kNoBGV0A) && col.selection_bit(kNoBGV0C) && col.selection_bit(kNoTPCLaserWarmUp) && col.selection_bit(kNoTPCHVdip);

      for (int iAlias = 0; iAlias < kNaliases; iAlias++) {
        if (!col.alias_bit(iAlias)) {
          continue;
        }
        histos.fill(HIST("hColCounterAll"), iAlias, 1);
        if ((!isINT1period && col.sel7()) || (isINT1period && sel1)) {
          histos.fill(HIST("hColCounterAcc"), iAlias, 1);
        }
      }

      bool mb = isMC;
      mb |= !isINT1period && col.alias_bit(kINT7);
      mb |= isINT1period && col.alias_bit(kINT1);
      // further checks just on minimum bias triggers
      if (!mb) {
        continue;
      }
      for (int i = 0; i < kNsel; i++) {
        histos.fill(HIST("hSelCounter"), i, col.selection_bit(i));
      }

      const auto& bc = col.bc_as<BCsRun2>();
      uint64_t globalBC = bc.globalBC();
      // uint64_t orbit = globalBC / nBCsPerOrbit;
      int localBC = globalBC % nBCsPerOrbit;
      histos.fill(HIST("hGlobalBcAll"), globalBC - bcSOR);
      // histos.fill(HIST("hOrbitAll"), orbit - orbitSOR);
      histos.fill(HIST("hBcAll"), localBC);
      if (col.selection_bit(kIsBBV0A) || col.selection_bit(kIsBBV0C)) {
        histos.fill(HIST("hGlobalBcFV0"), globalBC - bcSOR);
        // histos.fill(HIST("hOrbitFV0"), orbit - orbitSOR);
        histos.fill(HIST("hBcFV0"), localBC);
      }
      if (col.selection_bit(kIsBBT0A) || col.selection_bit(kIsBBT0C)) {
        histos.fill(HIST("hGlobalBcFT0"), globalBC - bcSOR);
        // histos.fill(HIST("hOrbitFT0"), orbit - orbitSOR);
        histos.fill(HIST("hBcFT0"), localBC);
      }
      if (col.selection_bit(kIsBBFDA) || col.selection_bit(kIsBBFDC)) {
        histos.fill(HIST("hGlobalBcFDD"), globalBC - bcSOR);
        // histos.fill(HIST("hOrbitFDD"), orbit - orbitSOR);
        histos.fill(HIST("hBcFDD"), localBC);
      }

      // Calculate V0 multiplicity per ring
      float multRingV0A[5] = {0.};
      float multRingV0C[4] = {0.};
      float multV0A = 0;
      float multV0C = 0;
      if (bc.has_fv0a()) {
        for (unsigned int i = 0; i < bc.fv0a().amplitude().size(); ++i) {
          int ring = bc.fv0a().channel()[i] / 8;
          multRingV0A[ring] += bc.fv0a().amplitude()[i];
          multV0A += bc.fv0a().amplitude()[i];
        }
      }

      if (bc.has_fv0c()) {
        for (unsigned int i = 0; i < bc.fv0c().amplitude().size(); ++i) {
          int ring = bc.fv0c().channel()[i] / 8;
          multRingV0C[ring] += bc.fv0c().amplitude()[i];
          multV0C += bc.fv0c().amplitude()[i];
        }
      }

      float timeZNA = bc.has_zdc() ? bc.zdc().timeZNA() : -999.f;
      float timeZNC = bc.has_zdc() ? bc.zdc().timeZNC() : -999.f;
      float timeV0A = bc.has_fv0a() ? bc.fv0a().time() : -999.f;
      float timeV0C = bc.has_fv0c() ? bc.fv0c().time() : -999.f;
      float timeT0A = bc.has_ft0() ? bc.ft0().timeA() : -999.f;
      float timeT0C = bc.has_ft0() ? bc.ft0().timeC() : -999.f;
      float timeFDA = bc.has_fdd() ? bc.fdd().timeA() : -999.f;
      float timeFDC = bc.has_fdd() ? bc.fdd().timeC() : -999.f;
      float znSum = timeZNA + timeZNC;
      float znDif = timeZNA - timeZNC;
      float ofSPD = bc.spdFiredChipsL0() + bc.spdFiredChipsL1();
      float onSPD = bc.spdFiredFastOrL0() + bc.spdFiredFastOrL1();
      float multV0M = multV0A + multV0C;
      float multRingV0C3 = multRingV0C[3];
      float multRingV0C012 = multV0C - multRingV0C3;
      float onV0M = bc.v0TriggerChargeA() + bc.v0TriggerChargeC();
      float ofV0M = multV0A + multV0C - multRingV0A[0];
      int spdClusters = bc.spdClustersL0() + bc.spdClustersL1();

      auto trackletsGrouped = tracklets->sliceByCached(aod::track::collisionId, col.globalIndex(), cache);
      int nTracklets = trackletsGrouped.size();

      float multFDA = 0;
      float multFDC = 0;
      float multT0A = bc.has_ft0() ? bc.ft0().sumAmpA() : -999.f;
      float multT0C = bc.has_ft0() ? bc.ft0().sumAmpC() : -999.f;

      if (bc.has_fdd()) {
        auto fdd = bc.fdd();
        for (const auto& amplitude : fdd.chargeA()) {
          multFDA += amplitude;
        }
        for (const auto& amplitude : fdd.chargeC()) {
          multFDC += amplitude;
        }
      }
      float multZNA = bc.has_zdc() ? bc.zdc().energyCommonZNA() : 0;
      float multZNC = bc.has_zdc() ? bc.zdc().energyCommonZNC() : 0;

      histos.fill(HIST("hMultV0Mcol"), multV0M);
      histos.fill(HIST("hMultV0Acol"), multV0A);
      histos.fill(HIST("hMultV0Ccol"), multV0C);
      histos.fill(HIST("hMultZNAcol"), multZNA);
      histos.fill(HIST("hMultZNCcol"), multZNC);
      histos.fill(HIST("hMultT0Acol"), multT0A);
      histos.fill(HIST("hMultT0Ccol"), multT0C);
      histos.fill(HIST("hMultFDAcol"), multFDA);
      histos.fill(HIST("hMultFDCcol"), multFDC);

      histos.fill(HIST("hTimeV0Acol"), timeV0A);
      histos.fill(HIST("hTimeV0Ccol"), timeV0C);
      histos.fill(HIST("hTimeZNAcol"), timeZNA);
      histos.fill(HIST("hTimeZNCcol"), timeZNC);
      histos.fill(HIST("hTimeT0Acol"), timeT0A);
      histos.fill(HIST("hTimeT0Ccol"), timeT0C);
      histos.fill(HIST("hTimeFDAcol"), timeFDA);
      histos.fill(HIST("hTimeFDCcol"), timeFDC);
      histos.fill(HIST("hTimeZACcol"), znDif, znSum);
      histos.fill(HIST("hSPDClsVsTklCol"), nTracklets, spdClusters);
      histos.fill(HIST("hSPDOnVsOfCol"), ofSPD, onSPD);
      histos.fill(HIST("hV0MOnVsOfCol"), ofV0M, onV0M);
      histos.fill(HIST("hV0C3vs012Col"), multRingV0C012, multRingV0C3);
      histos.fill(HIST("hV0C012vsTklCol"), nTracklets, multRingV0C012);

      // filling plots for accepted events
      bool accepted = 0;
      accepted |= !isINT1period & col.sel7();
      accepted |= isINT1period & sel1;
      if (!accepted) {
        continue;
      }

      histos.fill(HIST("hMultV0Macc"), multV0M);
      histos.fill(HIST("hMultV0Aacc"), multV0A);
      histos.fill(HIST("hMultV0Cacc"), multV0C);
      histos.fill(HIST("hMultZNAacc"), multZNA);
      histos.fill(HIST("hMultZNCacc"), multZNC);
      histos.fill(HIST("hMultT0Aacc"), multT0A);
      histos.fill(HIST("hMultT0Cacc"), multT0C);
      histos.fill(HIST("hMultFDAacc"), multFDA);
      histos.fill(HIST("hMultFDCacc"), multFDC);

      histos.fill(HIST("hTimeV0Aacc"), timeV0A);
      histos.fill(HIST("hTimeV0Cacc"), timeV0C);
      histos.fill(HIST("hTimeZNAacc"), timeZNA);
      histos.fill(HIST("hTimeZNCacc"), timeZNC);
      histos.fill(HIST("hTimeT0Aacc"), timeT0A);
      histos.fill(HIST("hTimeT0Cacc"), timeT0C);
      histos.fill(HIST("hTimeFDAacc"), timeFDA);
      histos.fill(HIST("hTimeFDCacc"), timeFDC);
      histos.fill(HIST("hTimeZACacc"), znDif, znSum);
      histos.fill(HIST("hSPDClsVsTklAcc"), nTracklets, spdClusters);
      histos.fill(HIST("hSPDOnVsOfAcc"), ofSPD, onSPD);
      histos.fill(HIST("hV0MOnVsOfAcc"), ofV0M, onV0M);
      histos.fill(HIST("hV0C3vs012Acc"), multRingV0C012, multRingV0C3);
      histos.fill(HIST("hV0C012vsTklAcc"), nTracklets, multRingV0C012);
    }
  }
  PROCESS_SWITCH(EventSelectionQaTask, processRun2, "Process Run2 event selection QA", true);

  Preslice<FullTracksIU> perCollision = aod::track::collisionId;
  // Preslice<ColEvSels> perFoundBC = aod::evsel::foundBCId;

  void processRun3(
    ColEvSels const& cols,
    FullTracksIU const& tracks,
    aod::AmbiguousTracks const& ambTracks,
    BCsRun3 const& bcs,
    aod::Zdcs const&,
    aod::FV0As const&,
    aod::FT0s const&,
    aod::FDDs const&)
  {
    int run = bcs.iteratorAt(0).runNumber();

    if (run != lastRun) {
      lastRun = run;
      int64_t tsSOR = 0; // dummy start-of-run timestamp for unanchored MC
      int64_t tsEOR = 1; // dummy end-of-run timestamp for unanchored MC
      if (run >= 500000) {
        auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(o2::ccdb::BasicCCDBManager::instance(), run);
        // first bc of the first orbit
        bcSOR = runInfo.orbitSOR * nBCsPerOrbit;
        // number of orbits per TF
        nOrbitsPerTF = runInfo.orbitsPerTF;
        // duration of TF in bcs
        nBCsPerTF = nOrbitsPerTF * nBCsPerOrbit;
        // first orbit
        orbitSOR = runInfo.orbitSOR;
        // total number of orbits
        nOrbits = runInfo.orbitEOR - runInfo.orbitSOR;
        // start-of-run timestamp
        tsSOR = runInfo.sor;
        // end-of-run timestamp
        tsEOR = runInfo.eor;

        // extract ITS ROF parameters
        int64_t ts = bcs.iteratorAt(0).timestamp();
        auto alppar = ccdb->getForTimeStamp<o2::itsmft::DPLAlpideParam<0>>("ITS/Config/AlpideParam", ts);
        rofOffset = alppar->roFrameBiasInBC;
        rofLength = alppar->roFrameLengthInBC;
        LOGP(info, "rofOffset={} rofLength={}", rofOffset, rofLength);
        LOGP(info, "nOrbitsPerTF={} nBCsPerTF={}", nOrbitsPerTF, nBCsPerTF);

        // bc patterns
        auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", (tsSOR + tsEOR) / 2);
        auto beamPatternA = grplhcif->getBunchFilling().getBeamPattern(0);
        auto beamPatternC = grplhcif->getBunchFilling().getBeamPattern(1);
        bcPatternA = beamPatternA & ~beamPatternC;
        bcPatternC = ~beamPatternA & beamPatternC;
        bcPatternB = beamPatternA & beamPatternC;

        // fill once
        for (int i = 0; i < nBCsPerOrbit; i++) {
          histos.fill(HIST("hBcA"), i, bcPatternA[i] ? 1. : 0.);
          histos.fill(HIST("hBcB"), i, bcPatternB[i] ? 1. : 0.);
          histos.fill(HIST("hBcC"), i, bcPatternC[i] ? 1. : 0.);
        }

        // fill ITS dead maps
        if (fillITSdeadStaveHists) {
          o2::itsmft::TimeDeadMap* itsDeadMap = ccdb->getForTimeStamp<o2::itsmft::TimeDeadMap>("ITS/Calib/TimeDeadMap", (tsSOR + tsEOR) / 2);
          auto itsDeadMapOrbits = itsDeadMap->getEvolvingMapKeys(); // roughly every second, ~350 TFs = 350x32 orbits
          if (itsDeadMapOrbits.size() > 0) {
            std::vector<double> itsDeadMapOrbitsDouble(itsDeadMapOrbits.begin(), itsDeadMapOrbits.end());
            const AxisSpec axisItsDeadMapOrbits{itsDeadMapOrbitsDouble};

            for (int l = 0; l < o2::itsmft::ChipMappingITS::NLayers; l++) {
              int nChips = o2::itsmft::ChipMappingITS::getNChipsOnLayer(l);
              double idFirstChip = o2::itsmft::ChipMappingITS::getFirstChipsOnLayer(l);
              // int nStaves = o2::itsmft::ChipMappingITS::getNStavesOnLr(l);
              // double idFirstStave = o2::itsmft::ChipMappingITS::getFirstStavesOnLr(l);
              histos.add(Form("hDeadChipsVsOrbitL%d", l), Form(";orbit; chip; Layer %d", l), kTH2C, {axisItsDeadMapOrbits, {nChips, idFirstChip, idFirstChip + nChips}});
              histos.add(Form("hNumberOfInactiveChipsVsOrbitL%d", l), Form(";orbit; Layer %d", l), kTH1I, {axisItsDeadMapOrbits});
            }

            std::vector<uint16_t> vClosest;
            std::bitset<o2::itsmft::ChipMappingITS::getNChips()> alwaysDeadChips;
            std::bitset<o2::itsmft::ChipMappingITS::getNChips()> deadChips;
            alwaysDeadChips.set();
            for (const auto& orbit : itsDeadMapOrbits) {
              itsDeadMap->getMapAtOrbit(orbit, vClosest);
              deadChips.reset();
              for (size_t iel = 0; iel < vClosest.size(); iel++) {
                uint16_t w1 = vClosest[iel];
                bool isLastInSequence = (w1 & 0x8000) == 0;
                uint16_t w2 = isLastInSequence ? w1 + 1 : vClosest[iel + 1];
                uint16_t chipId1 = w1 & 0x7FFF;
                uint16_t chipId2 = w2 & 0x7FFF;
                // dead chips are stored as ranges
                // vClosest contains first and last chip ids in the range
                // last chip id in the range is marked with 0x8000 bit set to 1
                for (int chipId = chipId1; chipId < chipId2; chipId++) {
                  histos.fill(HIST("hDeadChipsVsOrbitL0"), orbit, chipId, 1);
                  histos.fill(HIST("hDeadChipsVsOrbitL1"), orbit, chipId, 1);
                  histos.fill(HIST("hDeadChipsVsOrbitL2"), orbit, chipId, 1);
                  histos.fill(HIST("hDeadChipsVsOrbitL3"), orbit, chipId, 1);
                  histos.fill(HIST("hDeadChipsVsOrbitL4"), orbit, chipId, 1);
                  histos.fill(HIST("hDeadChipsVsOrbitL5"), orbit, chipId, 1);
                  histos.fill(HIST("hDeadChipsVsOrbitL6"), orbit, chipId, 1);
                  deadChips.set(chipId);
                }
              }
              alwaysDeadChips &= deadChips; // chips active in the current orbit are set to 0
            }
            // std::cout << alwaysDeadChips << std::endl;

            // filling histograms with number of inactive chips per layer vs orbit (ignoring always inactive)
            for (const auto& orbit : itsDeadMapOrbits) {
              itsDeadMap->getMapAtOrbit(orbit, vClosest);
              std::vector<int16_t> nInactiveChips(o2::itsmft::ChipMappingITS::NLayers, 0);
              for (size_t iel = 0; iel < vClosest.size(); iel++) {
                uint16_t w1 = vClosest[iel];
                bool isLastInSequence = (w1 & 0x8000) == 0;
                uint16_t w2 = isLastInSequence ? w1 + 1 : vClosest[iel + 1];
                uint16_t chipId1 = w1 & 0x7FFF;
                uint16_t chipId2 = w2 & 0x7FFF;
                for (int chipId = chipId1; chipId < chipId2; chipId++) {
                  if (alwaysDeadChips[chipId]) // skip always inactive chips
                    continue;
                  int32_t layer = o2::itsmft::ChipMappingITS::getLayer(chipId);
                  nInactiveChips[layer]++;
                }
              }
              histos.fill(HIST("hNumberOfInactiveChipsVsOrbitL0"), orbit, nInactiveChips[0]);
              histos.fill(HIST("hNumberOfInactiveChipsVsOrbitL1"), orbit, nInactiveChips[1]);
              histos.fill(HIST("hNumberOfInactiveChipsVsOrbitL2"), orbit, nInactiveChips[2]);
              histos.fill(HIST("hNumberOfInactiveChipsVsOrbitL3"), orbit, nInactiveChips[3]);
              histos.fill(HIST("hNumberOfInactiveChipsVsOrbitL4"), orbit, nInactiveChips[4]);
              histos.fill(HIST("hNumberOfInactiveChipsVsOrbitL5"), orbit, nInactiveChips[5]);
              histos.fill(HIST("hNumberOfInactiveChipsVsOrbitL6"), orbit, nInactiveChips[6]);
            }
          }
        } // end of fill ITS dead maps
      } // run >= 500000

      // create orbit-axis histograms on the fly with binning based on info from GRP if GRP is available
      // otherwise default orbitSOR and nOrbits will be used
      const AxisSpec axisOrbits{static_cast<int>(nOrbits / nOrbitsPerTF), 0., static_cast<double>(nOrbits), ""};
      histos.add("hOrbitAll", "", kTH1F, {axisOrbits});
      histos.add("hOrbitCol", "", kTH1F, {axisOrbits});
      histos.add("hOrbitAcc", "", kTH1F, {axisOrbits});
      histos.add("hOrbitTVX", "", kTH1F, {axisOrbits});
      histos.add("hOrbitFT0", "", kTH1F, {axisOrbits});
      histos.add("hOrbitFV0", "", kTH1F, {axisOrbits});
      histos.add("hOrbitFDD", "", kTH1F, {axisOrbits});
      histos.add("hOrbitZDC", "", kTH1F, {axisOrbits});
      histos.add("hOrbitColMC", "", kTH1F, {axisOrbits});

      const AxisSpec axisBCinTF{static_cast<int>(nBCsPerTF), 0, static_cast<double>(nBCsPerTF), "bc in TF"};
      histos.add("hNcontribVsBcInTF", ";bc in TF; n vertex contributors", kTH1F, {axisBCinTF});
      histos.add("hNcontribAfterCutsVsBcInTF", ";bc in TF; n vertex contributors", kTH1F, {axisBCinTF});
      histos.add("hNcolMCVsBcInTF", ";bc in TF; n MC collisions", kTH1F, {axisBCinTF});
      histos.add("hNcolVsBcInTF", ";bc in TF; n collisions", kTH1F, {axisBCinTF});
      histos.add("hNcolVsBcInTFafterTFborderCut", ";bc in TF; n collisions", kTH1F, {axisBCinTF});
      histos.add("hNtvxVsBcInTF", ";bc in TF; n TVX triggers", kTH1F, {axisBCinTF});

      double minSec = floor(tsSOR / 1000.);
      double maxSec = ceil(tsEOR / 1000.);
      const AxisSpec axisSeconds{maxSec - minSec < 1000 ? static_cast<int>(maxSec - minSec) : 1000, minSec, maxSec, "seconds"};
      const AxisSpec axisBcDif{600, -300., 300., "bc difference"};
      histos.add("hSecondsTVXvsBcDif", "", kTH2F, {axisSeconds, axisBcDif});
      histos.add("hSecondsTVXvsBcDifAll", "", kTH2F, {axisSeconds, axisBcDif});
    }

    // background studies
    for (const auto& bc : bcs) {
      // make sure previous bcs are empty to clean-up other activity
      uint64_t globalBC = bc.globalBC();
      int deltaIndex = 0;  // backward move counts
      int deltaBC = 0;     // current difference wrt globalBC
      int maxDeltaBC = 10; // maximum difference
      bool pastActivityFT0 = 0;
      bool pastActivityFDD = 0;
      bool pastActivityFV0 = 0;
      while (deltaBC < maxDeltaBC) {
        if (bc.globalIndex() - deltaIndex < 0) {
          break;
        }
        deltaIndex++;
        const auto& bcPast = bcs.iteratorAt(bc.globalIndex() - deltaIndex);
        deltaBC = globalBC - bcPast.globalBC();
        if (deltaBC < maxDeltaBC) {
          pastActivityFT0 |= bcPast.has_ft0();
          pastActivityFV0 |= bcPast.has_fv0a();
          pastActivityFDD |= bcPast.has_fdd();
        }
      }

      bool pastActivity = pastActivityFT0 | pastActivityFV0 | pastActivityFDD;

      int localBC = bc.globalBC() % nBCsPerOrbit;
      float timeV0A = bc.has_fv0a() ? bc.fv0a().time() : -999.f;
      float timeT0A = bc.has_ft0() ? bc.ft0().timeA() : -999.f;
      float timeT0C = bc.has_ft0() ? bc.ft0().timeC() : -999.f;
      float timeFDA = bc.has_fdd() ? bc.fdd().timeA() : -999.f;
      float timeFDC = bc.has_fdd() ? bc.fdd().timeC() : -999.f;
      if (bcPatternA[(localBC + 5) % nBCsPerOrbit] && !pastActivity && !bc.has_ft0()) {
        histos.fill(HIST("hTimeFDAbga"), timeFDA);
        histos.fill(HIST("hTimeFDCbga"), timeFDC);
      }
      if (bcPatternC[(localBC + 5) % nBCsPerOrbit] && !pastActivity && !bc.has_ft0()) {
        histos.fill(HIST("hTimeFDAbgc"), timeFDA);
        histos.fill(HIST("hTimeFDCbgc"), timeFDC);
      }
      if (bcPatternA[(localBC + 1) % nBCsPerOrbit] && !pastActivity && !bc.has_ft0()) {
        histos.fill(HIST("hTimeT0Abga"), timeT0A);
        histos.fill(HIST("hTimeT0Cbga"), timeT0C);
        histos.fill(HIST("hTimeV0Abga"), timeV0A);
      }
      if (bcPatternC[(localBC + 1) % nBCsPerOrbit] && !pastActivity && !bc.has_ft0()) {
        histos.fill(HIST("hTimeT0Abgc"), timeT0A);
        histos.fill(HIST("hTimeT0Cbgc"), timeT0C);
      }
    }

    // vectors of TVX flags used for past-future studies
    int nBCs = bcs.size();
    std::vector<bool> vIsTVX(nBCs, 0);
    std::vector<uint64_t> vGlobalBCs(nBCs, 0);

    // bc-based event selection qa
    for (const auto& bc : bcs) {
      if (!bc.has_ft0())
        continue;
      float multT0A = bc.ft0().sumAmpA();
      float multT0C = bc.ft0().sumAmpC();
      histos.fill(HIST("hMultT0Mref"), multT0A + multT0C);
      if (!bc.selection_bit(kIsTriggerTVX))
        continue;
      histos.fill(HIST("hMultT0Mtvx"), multT0A + multT0C);
      histos.fill(HIST("hMultT0Atvx"), multT0A);
      histos.fill(HIST("hMultT0Ctvx"), multT0C);
      if (!bc.selection_bit(kIsBBZAC))
        continue;
      histos.fill(HIST("hMultT0Mzac"), multT0A + multT0C);
      histos.fill(HIST("hMultT0Azac"), multT0A);
      histos.fill(HIST("hMultT0Czac"), multT0C);
    }

    // bc-based event selection qa
    for (const auto& bc : bcs) {
      for (int iAlias = 0; iAlias < kNaliases; iAlias++) {
        histos.fill(HIST("hBcCounterAll"), iAlias, bc.alias_bit(iAlias));
      }
      uint64_t globalBC = bc.globalBC();
      uint64_t orbit = globalBC / nBCsPerOrbit;
      int localBC = globalBC % nBCsPerOrbit;
      float timeZNA = bc.has_zdc() ? bc.zdc().timeZNA() : -999.f;
      float timeZNC = bc.has_zdc() ? bc.zdc().timeZNC() : -999.f;
      float timeV0A = bc.has_fv0a() ? bc.fv0a().time() : -999.f;
      float timeT0A = bc.has_ft0() ? bc.ft0().timeA() : -999.f;
      float timeT0C = bc.has_ft0() ? bc.ft0().timeC() : -999.f;
      float timeFDA = bc.has_fdd() ? bc.fdd().timeA() : -999.f;
      float timeFDC = bc.has_fdd() ? bc.fdd().timeC() : -999.f;
      histos.fill(HIST("hTimeV0Aall"), timeV0A);
      histos.fill(HIST("hTimeZNAall"), timeZNA);
      histos.fill(HIST("hTimeZNCall"), timeZNC);
      histos.fill(HIST("hTimeT0Aall"), timeT0A);
      histos.fill(HIST("hTimeT0Call"), timeT0C);
      histos.fill(HIST("hTimeFDAall"), timeFDA);
      histos.fill(HIST("hTimeFDCall"), timeFDC);
      if (bcPatternB[localBC]) {
        histos.fill(HIST("hTimeV0Aref"), timeV0A);
        histos.fill(HIST("hTimeZNAref"), timeZNA);
        histos.fill(HIST("hTimeZNCref"), timeZNC);
        histos.fill(HIST("hTimeT0Aref"), timeT0A);
        histos.fill(HIST("hTimeT0Cref"), timeT0C);
        histos.fill(HIST("hTimeFDAref"), timeFDA);
        histos.fill(HIST("hTimeFDCref"), timeFDC);
      }

      histos.fill(HIST("hGlobalBcAll"), globalBC - bcSOR);
      histos.fill(HIST("hOrbitAll"), orbit - orbitSOR);
      histos.fill(HIST("hBcAll"), localBC);

      if (bc.selection_bit(kIsTriggerTVX)) {
        histos.fill(HIST("hOrbitTVX"), orbit - orbitSOR);
        histos.fill(HIST("hBcTVX"), localBC);
      }

      // FV0
      if (bc.has_fv0a()) {
        histos.fill(HIST("hGlobalBcFV0"), globalBC - bcSOR);
        histos.fill(HIST("hOrbitFV0"), orbit - orbitSOR);
        histos.fill(HIST("hBcFV0"), localBC);
        float multV0A = 0;
        for (const auto& amplitude : bc.fv0a().amplitude()) {
          multV0A += amplitude;
        }
        histos.fill(HIST("hMultV0Aall"), multV0A);
        if (bcPatternB[localBC]) {
          histos.fill(HIST("hMultV0Aref"), multV0A);
        }
      }

      // FT0
      if (bc.has_ft0()) {
        histos.fill(HIST("hGlobalBcFT0"), globalBC - bcSOR);
        histos.fill(HIST("hOrbitFT0"), orbit - orbitSOR);
        histos.fill(HIST("hBcFT0"), localBC);
        float multT0A = bc.ft0().sumAmpA();
        float multT0C = bc.ft0().sumAmpC();
        histos.fill(HIST("hMultT0Aall"), multT0A);
        histos.fill(HIST("hMultT0Call"), multT0C);
        if (bcPatternB[localBC]) {
          histos.fill(HIST("hMultT0Aref"), multT0A);
          histos.fill(HIST("hMultT0Cref"), multT0C);
        }
        if (bc.selection_bit(kIsTriggerTVX)) {
          int64_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;
          histos.fill(HIST("hNtvxVsBcInTF"), bcInTF);
        }
        if (!bc.selection_bit(kNoBGFDA) && bc.selection_bit(kIsTriggerTVX)) {
          histos.fill(HIST("hMultT0Abga"), multT0A);
          histos.fill(HIST("hMultT0Cbga"), multT0C);
        }
        if (!bc.selection_bit(kNoBGFDC) && bc.selection_bit(kIsTriggerTVX)) {
          histos.fill(HIST("hMultT0Abgc"), multT0A);
          histos.fill(HIST("hMultT0Cbgc"), multT0C);
        }
      }

      // FDD
      if (bc.has_fdd()) {
        histos.fill(HIST("hGlobalBcFDD"), globalBC - bcSOR);
        histos.fill(HIST("hOrbitFDD"), orbit - orbitSOR);
        histos.fill(HIST("hBcFDD"), localBC);

        auto fdd = bc.fdd();
        float multFDA = 0;
        for (const auto& amplitude : fdd.chargeA()) {
          multFDA += amplitude;
        }
        float multFDC = 0;
        for (const auto& amplitude : fdd.chargeC()) {
          multFDC += amplitude;
        }
        histos.fill(HIST("hMultFDAall"), multFDA);
        histos.fill(HIST("hMultFDCall"), multFDC);
        if (bcPatternB[localBC]) {
          histos.fill(HIST("hMultFDAref"), multFDA);
          histos.fill(HIST("hMultFDCref"), multFDC);
        }
      }

      // ZDC
      if (bc.has_zdc()) {
        histos.fill(HIST("hGlobalBcZDC"), globalBC - bcSOR);
        histos.fill(HIST("hOrbitZDC"), orbit - orbitSOR);
        histos.fill(HIST("hBcZDC"), localBC);
        float multZNA = bc.zdc().energyCommonZNA();
        float multZNC = bc.zdc().energyCommonZNC();
        histos.fill(HIST("hMultZNAall"), multZNA);
        histos.fill(HIST("hMultZNCall"), multZNC);
        if (bcPatternB[localBC]) {
          histos.fill(HIST("hMultZNAref"), multZNA);
          histos.fill(HIST("hMultZNCref"), multZNC);
        }
      }

      // fill TVX flags for past-future searches
      int indexBc = bc.globalIndex();
      vIsTVX[indexBc] = bc.selection_bit(kIsTriggerTVX);
      vGlobalBCs[indexBc] = globalBC;
    }

    // map for pileup checks
    std::vector<int> vCollisionsPerBc(bcs.size(), 0);
    for (const auto& col : cols) {
      if (col.foundBCId() < 0 || col.foundBCId() >= bcs.size())
        continue;
      vCollisionsPerBc[col.foundBCId()]++;
    }

    // build map from track index to ambiguous track index
    std::unordered_map<int32_t, int32_t> mapAmbTrIds;
    for (const auto& ambTrack : ambTracks) {
      mapAmbTrIds[ambTrack.trackId()] = ambTrack.globalIndex();
    }

    // create maps from globalBC to bc index for TVX or FT0-OR fired bcs
    // to be used for closest TVX (FT0-OR) searches
    std::map<int64_t, int32_t> mapGlobalBcWithTVX;
    std::map<int64_t, int32_t> mapGlobalBcWithTOR;
    for (const auto& bc : bcs) {
      int64_t globalBC = bc.globalBC();
      // skip non-colliding bcs for data and anchored runs
      if (run >= 500000 && bcPatternB[globalBC % nBCsPerOrbit] == 0) {
        continue;
      }
      if (bc.selection_bit(kIsBBT0A) || bc.selection_bit(kIsBBT0C)) {
        mapGlobalBcWithTOR[globalBC] = bc.globalIndex();
      }
      if (bc.selection_bit(kIsTriggerTVX)) {
        mapGlobalBcWithTVX[globalBC] = bc.globalIndex();
      }
    }

    // Fill track bc distributions (all tracks including ambiguous)
    for (const auto& track : tracks) {
      auto mapAmbTrIdsIt = mapAmbTrIds.find(track.globalIndex());
      int ambTrId = mapAmbTrIdsIt == mapAmbTrIds.end() ? -1 : mapAmbTrIdsIt->second;

      // special check to avoid crashes (in particular, on some MC Pb-Pb datasets)
      // (related to shifts in ambiguous tracks association to bc slices (off by 1) - see https://mattermost.web.cern.ch/alice/pl/g9yaaf3tn3g4pgn7c1yex9copy
      if (ambTrId >= 0 && (ambTracks.iteratorAt(ambTrId).bcIds()[0] >= bcs.size()))
        continue;

      int indexBc = ambTrId < 0 ? track.collision_as<ColEvSels>().bc_as<BCsRun3>().globalIndex() : ambTracks.iteratorAt(ambTrId).bc_as<BCsRun3>().begin().globalIndex();
      auto bc = bcs.iteratorAt(indexBc);
      int64_t globalBC = bc.globalBC() + floor(track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS);

      int32_t indexClosestTVX = findClosest(globalBC, mapGlobalBcWithTVX);
      int bcDiff = static_cast<int>(globalBC - vGlobalBCs[indexClosestTVX]);
      if (track.hasTOF() || track.hasTRD() || !track.hasITS() || !track.hasTPC() || track.pt() < 1)
        continue;
      histos.fill(HIST("hTrackBcDiffVsEtaAll"), track.eta(), bcDiff);
      if (track.eta() < -0.2 || track.eta() > 0.2)
        continue;
      histos.fill(HIST("hSecondsTVXvsBcDifAll"), bc.timestamp() / 1000., bcDiff);
    }

    // collision-based event selection qa
    std::vector<int64_t> vFoundGlobalBC(cols.size(), 0);   // global BCs for collisions
    std::vector<float> vCollVz(cols.size(), 0);            // vector with vZ positions for each collision
    std::vector<bool> vIsSel8(cols.size(), 0);             // vector with sel8 decisions
    std::vector<int> vTracksITS567perColl(cols.size(), 0); // counter of tracks per collision for occupancy studies
    for (const auto& col : cols) {
      for (int iAlias = 0; iAlias < kNaliases; iAlias++) {
        if (!col.alias_bit(iAlias)) {
          continue;
        }
        histos.fill(HIST("hColCounterAll"), iAlias, 1);
        if (!col.sel8()) {
          continue;
        }
        histos.fill(HIST("hColCounterAcc"), iAlias, 1);
      }

      for (int i = 0; i < kNsel; i++) {
        histos.fill(HIST("hSelCounter"), i, col.selection_bit(i));
      }

      auto bc = col.bc_as<BCsRun3>();
      uint64_t globalBC = bc.globalBC();
      uint64_t orbit = globalBC / nBCsPerOrbit;
      int localBC = globalBC % nBCsPerOrbit;
      histos.fill(HIST("hGlobalBcCol"), globalBC - bcSOR);
      histos.fill(HIST("hOrbitCol"), orbit - orbitSOR);
      histos.fill(HIST("hBcCol"), localBC);
      if (col.sel8()) {
        histos.fill(HIST("hOrbitAcc"), orbit - orbitSOR);
      }

      int32_t colIndex = col.globalIndex();
      vFoundGlobalBC[colIndex] = globalBC;
      vCollVz[colIndex] = col.posZ();
      vIsSel8[colIndex] = col.sel8();

      // search for nearest ft0a&ft0c entry
      int32_t indexClosestTVX = findClosest(globalBC, mapGlobalBcWithTVX);
      int bcDiff = static_cast<int>(globalBC - vGlobalBCs[indexClosestTVX]);

      int nContributors = col.numContrib();
      float timeRes = col.collisionTimeRes();
      int64_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;
      histos.fill(HIST("hNcontribCol"), nContributors);
      histos.fill(HIST("hNcontribVsBcInTF"), bcInTF, nContributors);
      histos.fill(HIST("hNcolVsBcInTF"), bcInTF);
      if (col.selection_bit(kNoTimeFrameBorder))
        histos.fill(HIST("hNcolVsBcInTFafterTFborderCut"), bcInTF);
      histos.fill(HIST("hColBcDiffVsNcontrib"), nContributors, bcDiff);
      histos.fill(HIST("hColTimeResVsNcontrib"), nContributors, timeRes);
      if (!col.selection_bit(kIsVertexITSTPC)) {
        histos.fill(HIST("hColBcDiffVsNcontribITSonly"), nContributors, bcDiff);
        histos.fill(HIST("hColTimeResVsNcontribITSonly"), nContributors, timeRes);
      }
      if (col.selection_bit(kIsVertexTOFmatched)) {
        histos.fill(HIST("hColBcDiffVsNcontribWithTOF"), nContributors, bcDiff);
        histos.fill(HIST("hColTimeResVsNcontribWithTOF"), nContributors, timeRes);
        histos.fill(HIST("hNcontribColTOF"), nContributors);
        histos.fill(HIST("hBcColTOF"), localBC);
        if (col.sel8()) {
          histos.fill(HIST("hNcontribAccTOF"), nContributors);
        }
      }
      if (col.selection_bit(kIsVertexTRDmatched)) {
        histos.fill(HIST("hColBcDiffVsNcontribWithTRD"), nContributors, bcDiff);
        histos.fill(HIST("hColTimeResVsNcontribWithTRD"), nContributors, timeRes);
        histos.fill(HIST("hNcontribColTRD"), nContributors);
        histos.fill(HIST("hBcColTRD"), localBC);
        if (col.sel8()) {
          histos.fill(HIST("hNcontribAccTRD"), nContributors);
        }
      }

      const auto& foundBC = col.foundBC_as<BCsRun3>();

      float timeZNA = foundBC.has_zdc() ? foundBC.zdc().timeZNA() : -999.f;
      float timeZNC = foundBC.has_zdc() ? foundBC.zdc().timeZNC() : -999.f;
      float timeV0A = foundBC.has_fv0a() ? foundBC.fv0a().time() : -999.f;
      float timeT0A = foundBC.has_ft0() ? foundBC.ft0().timeA() : -999.f;
      float timeT0C = foundBC.has_ft0() ? foundBC.ft0().timeC() : -999.f;
      float timeFDA = foundBC.has_fdd() ? foundBC.fdd().timeA() : -999.f;
      float timeFDC = foundBC.has_fdd() ? foundBC.fdd().timeC() : -999.f;
      float znSum = timeZNA + timeZNC;
      float znDif = timeZNA - timeZNC;

      histos.fill(HIST("hTimeV0Acol"), timeV0A);
      histos.fill(HIST("hTimeZNAcol"), timeZNA);
      histos.fill(HIST("hTimeZNCcol"), timeZNC);
      histos.fill(HIST("hTimeT0Acol"), timeT0A);
      histos.fill(HIST("hTimeT0Ccol"), timeT0C);
      histos.fill(HIST("hTimeFDAcol"), timeFDA);
      histos.fill(HIST("hTimeFDCcol"), timeFDC);
      histos.fill(HIST("hTimeZACcol"), znDif, znSum);

      // FT0
      float multT0A = foundBC.has_ft0() ? foundBC.ft0().sumAmpA() : -999.f;
      float multT0C = foundBC.has_ft0() ? foundBC.ft0().sumAmpC() : -999.f;

      // FV0
      float multV0A = 0;
      if (foundBC.has_fv0a()) {
        for (const auto& amplitude : foundBC.fv0a().amplitude()) {
          multV0A += amplitude;
        }
      }
      // FDD
      float multFDA = 0;
      float multFDC = 0;
      if (foundBC.has_fdd()) {
        auto fdd = foundBC.fdd();
        for (const auto& amplitude : fdd.chargeA()) {
          multFDA += amplitude;
        }
        for (const auto& amplitude : fdd.chargeC()) {
          multFDC += amplitude;
        }
      }

      // ZDC
      float multZNA = foundBC.has_zdc() ? foundBC.zdc().energyCommonZNA() : -999.f;
      float multZNC = foundBC.has_zdc() ? foundBC.zdc().energyCommonZNC() : -999.f;

      histos.fill(HIST("hMultT0Acol"), multT0A);
      histos.fill(HIST("hMultT0Ccol"), multT0C);
      histos.fill(HIST("hMultV0Acol"), multV0A);
      histos.fill(HIST("hMultFDAcol"), multFDA);
      histos.fill(HIST("hMultFDCcol"), multFDC);
      histos.fill(HIST("hMultZNAcol"), multZNA);
      histos.fill(HIST("hMultZNCcol"), multZNC);

      // count tracks of different types
      auto tracksGrouped = tracks.sliceBy(perCollision, colIndex);
      int nPV = 0;
      int nTPConly = 0;
      // int nTPConlyWithDeDxCut = 0;
      int nTPConlyNoITS = 0;
      int nContributorsAfterEtaTPCCuts = 0;
      int nContributorsAfterEtaTPCLooseCuts = 0;

      int nTracksITS = 0;
      int nTracksITSTPC = 0;

      bool isTVX = col.selection_bit(kIsTriggerTVX);

      int occupancyByTracks = col.trackOccupancyInTimeRange();
      float occupancyByFT0C = col.ft0cOccupancyInTimeRange();

      for (const auto& track : tracksGrouped) {
        int trackBcDiff = bcDiff + track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS;

        if (track.hasTPC() && std::fabs(track.eta()) < 0.8 && track.pt() > 0.2 && track.tpcNClsFound() > 50 && track.tpcNClsCrossedRows() > 50 && track.tpcChi2NCl() < 4) {
          nTPConly++;
          // if (track.tpcSignal() > 20)
          // nTPConlyWithDeDxCut++;
          if (!track.hasITS())
            nTPConlyNoITS++;
        }

        if (std::fabs(track.eta()) < 0.8 && track.pt() > 0.2) {
          if (track.hasITS()) {
            nTracksITS++;
            if (track.hasTPC())
              nTracksITSTPC++;
          }
        }

        if (!track.isPVContributor())
          continue;

        if (track.itsNCls() >= 5)
          vTracksITS567perColl[colIndex]++;

        // high-quality contributors for ROF border QA and occupancy study
        if (std::fabs(track.eta()) < 0.8 && track.pt() > 0.2 && track.itsNCls() >= 5) {
          nPV++;
          if (track.hasTPC()) {
            nContributorsAfterEtaTPCLooseCuts++;

            if (!isLowFlux && fillTPCnClsVsOccupancyHists && col.sel8() && col.selection_bit(kNoSameBunchPileup) && fabs(col.posZ()) < 10 && occupancyByTracks >= 0) {
              histos.fill(HIST("occupancyQA/tpcNClsFound_vs_V0A_vs_occupancy"), multV0A, track.tpcNClsFound(), occupancyByTracks);
              histos.fill(HIST("occupancyQA/tpcNClsFindable_vs_V0A_vs_occupancy"), multV0A, track.tpcNClsFindable(), occupancyByTracks);
              histos.fill(HIST("occupancyQA/tpcNClsShared_vs_V0A_vs_occupancy"), multV0A, track.tpcNClsShared(), occupancyByTracks);
              histos.fill(HIST("occupancyQA/tpcChi2_vs_V0A_vs_occupancy"), multV0A, track.tpcChi2NCl(), occupancyByTracks);
              int tpcNClsFindableMinusCrossedRowsCorrected = track.tpcNClsFindableMinusCrossedRows();
              // correct for a buggy behaviour due to int8 and uint8 difference:
              if (tpcNClsFindableMinusCrossedRowsCorrected < -70)
                tpcNClsFindableMinusCrossedRowsCorrected += 256;
              histos.fill(HIST("occupancyQA/tpcNCrossedRows_vs_V0A_vs_occupancy"), multV0A, track.tpcNClsFindable() - tpcNClsFindableMinusCrossedRowsCorrected, occupancyByTracks);
            }
          } // end of hasTPC
          if (col.sel8() && fabs(col.posZ()) < 10 && track.tpcNClsFound() > 50 && track.tpcNClsCrossedRows() > 80 && track.itsChi2NCl() < 36 && track.tpcChi2NCl() < 4) {
            nContributorsAfterEtaTPCCuts++;
            // ROF border QA
            histos.fill(HIST("ITSROFborderQA/hFoundBC_kTVX_nITSlayers_for_ITSTPCtracks"), localBC, track.itsNCls());
            histos.fill(HIST("ITSROFborderQA/hFoundBC_kTVX_counter_ITSTPCtracks"), localBC);
          }
        }
        if (!track.hasTPC())
          histos.fill(HIST("hITStrackBcDiff"), trackBcDiff);
        if (track.hasTOF()) {
          histos.fill(HIST("hBcTrackTOF"), (globalBC + TMath::FloorNint(track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS)) % nBCsPerOrbit);
        } else if (track.hasTRD()) {
          histos.fill(HIST("hBcTrackTRD"), (globalBC + TMath::Nint(track.trackTime() / o2::constants::lhc::LHCBunchSpacingNS)) % nBCsPerOrbit);
        }
        if (track.hasTOF() || track.hasTRD() || !track.hasITS() || !track.hasTPC() || track.pt() < 1)
          continue;
        histos.fill(HIST("hTrackBcDiffVsEta"), track.eta(), trackBcDiff);
        if (track.eta() < -0.2 || track.eta() > 0.2)
          continue;
        histos.fill(HIST("hSecondsTVXvsBcDif"), bc.timestamp() / 1000., trackBcDiff);
      } // end of track loop

      histos.fill(HIST("hNcontribAfterCutsVsBcInTF"), bcInTF, nContributorsAfterEtaTPCCuts);

      if (!isLowFlux && col.sel8() && col.selection_bit(kNoSameBunchPileup) && fabs(col.posZ()) < 10) {
        histos.fill(HIST("occupancyQA/hOccupancyByTracks"), occupancyByTracks);
        histos.fill(HIST("occupancyQA/hOccupancyByFT0C"), occupancyByFT0C);
        if (occupancyByTracks >= 0) {
          histos.fill(HIST("occupancyQA/hOccupancyByFT0CvsByTracks"), occupancyByTracks, occupancyByFT0C);
          histos.fill(HIST("occupancyQA/hNumTracksPV_vs_V0A_vs_occupancy"), multV0A, nPV, occupancyByTracks);
          histos.fill(HIST("occupancyQA/hNumTracksPVTPC_vs_V0A_vs_occupancy"), multV0A, nContributorsAfterEtaTPCCuts, occupancyByTracks);
          histos.fill(HIST("occupancyQA/hNumTracksPVTPCLooseCuts_vs_V0A_vs_occupancy"), multV0A, nContributorsAfterEtaTPCLooseCuts, occupancyByTracks);
          histos.fill(HIST("occupancyQA/hNumTracksITS_vs_V0A_vs_occupancy"), multV0A, nTracksITS, occupancyByTracks);
          histos.fill(HIST("occupancyQA/hNumTracksITSTPC_vs_V0A_vs_occupancy"), multV0A, nTracksITSTPC, occupancyByTracks);
          if (col.selection_bit(kNoCollInTimeRangeNarrow)) {
            histos.fill(HIST("occupancyQA/hNumTracksPV_vs_V0A_vs_occupancy_NarrowDeltaTimeCut"), multV0A, nPV, occupancyByTracks);
            histos.fill(HIST("occupancyQA/hNumTracksPVTPC_vs_V0A_vs_occupancy_NarrowDeltaTimeCut"), multV0A, nContributorsAfterEtaTPCCuts, occupancyByTracks);
          }
          if (col.selection_bit(kNoCollInTimeRangeStandard)) {
            histos.fill(HIST("occupancyQA/hNumTracksPV_vs_V0A_vs_occupancy_StandardDeltaTimeCut"), multV0A, nPV, occupancyByTracks);
            histos.fill(HIST("occupancyQA/hNumTracksPVTPC_vs_V0A_vs_occupancy_StandardDeltaTimeCut"), multV0A, nContributorsAfterEtaTPCCuts, occupancyByTracks);
          }
          if (col.selection_bit(kIsGoodITSLayersAll)) {
            histos.fill(HIST("occupancyQA/hNumTracksPV_vs_V0A_vs_occupancy_GoodITSLayersAllCut"), multV0A, nPV, occupancyByTracks);
            histos.fill(HIST("occupancyQA/hNumTracksPVTPC_vs_V0A_vs_occupancy_GoodITSLayersAllCut"), multV0A, nContributorsAfterEtaTPCCuts, occupancyByTracks);
          }
          histos.fill(HIST("occupancyQA/hNumTracksTPConly_vs_V0A_vs_occupancy"), multV0A, nTPConly, occupancyByTracks);
          histos.fill(HIST("occupancyQA/hNumTracksTPConlyNoITS_vs_V0A_vs_occupancy"), multV0A, nTPConlyNoITS, occupancyByTracks);

          // dE/dx QA for a narrow pT bin
          for (const auto& track : tracksGrouped) {
            if (!track.isPVContributor())
              continue;
            if (std::fabs(track.eta()) < 0.8 && track.pt() > 0.2 && track.itsNCls() >= 5) {
              float signedP = track.sign() * track.tpcInnerParam();
              if (std::fabs(signedP) > 0.38 && std::fabs(signedP) < 0.4 && track.tpcNClsFound() > 50 && track.tpcNClsCrossedRows() > 80 && track.itsChi2NCl() < 36 && track.tpcChi2NCl() < 4) {
                float dEdx = track.tpcSignal();
                histos.fill(HIST("occupancyQA/dEdx_vs_centr_vs_occup_narrow_p_win"), nPV, occupancyByTracks, dEdx);
              }
            }
          }
        }
      }

      // filling plots for events passing basic TVX selection
      if (!isTVX) {
        continue;
      }
      histos.fill(HIST("hMultT0MVsNcontribTVX"), multT0A + multT0C, nContributors);
      histos.fill(HIST("hMultV0AVsNcontribTVX"), multV0A, nContributors);

      // z-vertex from FT0 vs PV
      if (foundBC.has_ft0()) {
        histos.fill(HIST("hVtxFT0VsVtxCol"), foundBC.ft0().posZ(), col.posZ());
        histos.fill(HIST("hVtxFT0MinusVtxCol"), foundBC.ft0().posZ() - col.posZ());
        histos.fill(HIST("hVtxFT0MinusVtxColVsMultT0M"), foundBC.ft0().posZ() - col.posZ(), multT0A + multT0C);
      }

      int foundLocalBC = foundBC.globalBC() % nBCsPerOrbit;

      if (col.selection_bit(kNoITSROFrameBorder)) {
        histos.fill(HIST("hMultT0MVsNcontribTVXROFcuts"), multT0A + multT0C, nContributors);
        histos.fill(HIST("hMultV0AVsNcontribTVXROFcuts"), multV0A, nContributors);
      }

      if (col.selection_bit(kNoTimeFrameBorder)) {
        histos.fill(HIST("hMultT0MVsNcontribTVXTFcuts"), multT0A + multT0C, nContributors);
        histos.fill(HIST("hMultV0AVsNcontribTVXTFcuts"), multV0A, nContributors);

        // histos.fill(HIST("hFoundBcForMultV0AVsNcontribAcc"), foundLocalBC);
        histos.fill(HIST("hFoundBc"), foundLocalBC);
        histos.fill(HIST("hFoundBcNcontrib"), foundLocalBC, nContributors);
        if (col.selection_bit(kIsVertexTOFmatched)) {
          histos.fill(HIST("hFoundBcTOF"), foundLocalBC);
          histos.fill(HIST("hFoundBcNcontribTOF"), foundLocalBC, nContributors);
        }
        if (nContributors < 0.043 * multV0A - 860) {
          histos.fill(HIST("hFoundBcForMultV0AVsNcontribOutliers"), foundLocalBC);
        }
        if (col.selection_bit(kNoITSROFrameBorder)) {
          histos.fill(HIST("hMultT0MVsNcontribTVXTFROFcuts"), multT0A + multT0C, nContributors);
          histos.fill(HIST("hMultV0AVsNcontribTVXTFROFcuts"), multV0A, nContributors);

          histos.fill(HIST("hFoundBcAfterROFborderCut"), foundLocalBC);
        }
      }

      // filling plots for accepted events
      if (!col.sel8()) {
        continue;
      }

      if (col.selection_bit(kIsVertexITSTPC)) {
        histos.fill(HIST("hMultV0AVsNcontribIsVertexITSTPC"), multV0A, nContributors);
        if (col.selection_bit(kNoSameBunchPileup)) {
          histos.fill(HIST("hMultV0AVsNcontribGood"), multV0A, nContributors);
        }
      }

      if (!col.selection_bit(kNoSameBunchPileup)) {
        histos.fill(HIST("hMultT0Mpup"), multT0A + multT0C);
      }

      // histos.fill(HIST("hMultT0MVsNcontribAcc"), multT0A + multT0C, nContributors);
      histos.fill(HIST("hTimeV0Aacc"), timeV0A);
      histos.fill(HIST("hTimeZNAacc"), timeZNA);
      histos.fill(HIST("hTimeZNCacc"), timeZNC);
      histos.fill(HIST("hTimeT0Aacc"), timeT0A);
      histos.fill(HIST("hTimeT0Cacc"), timeT0C);
      histos.fill(HIST("hTimeFDAacc"), timeFDA);
      histos.fill(HIST("hTimeFDCacc"), timeFDC);
      histos.fill(HIST("hTimeZACacc"), znDif, znSum);
      histos.fill(HIST("hMultT0Aacc"), multT0A);
      histos.fill(HIST("hMultT0Cacc"), multT0C);
      histos.fill(HIST("hMultV0Aacc"), multV0A);
      histos.fill(HIST("hMultFDAacc"), multFDA);
      histos.fill(HIST("hMultFDCacc"), multFDC);
      histos.fill(HIST("hMultZNAacc"), multZNA);
      histos.fill(HIST("hMultZNCacc"), multZNC);
      histos.fill(HIST("hNcontribAcc"), nContributors);
    } // collisions

    // ### in-ROF occupancy QA
    if (!isLowFlux) {
      std::vector<std::vector<int>> vCollsInSameITSROF;
      // save indices of collisions in same ROF
      for (const auto& col : cols) {
        int32_t colIndex = col.globalIndex();
        int64_t foundGlobalBC = vFoundGlobalBC[colIndex];
        int64_t tfId = (foundGlobalBC - bcSOR) / nBCsPerTF;
        int64_t rofId = (foundGlobalBC + 3564 - rofOffset) / rofLength;
        std::vector<int> vAssocToSameROF;
        // find all collisions in the same ROF before a given collision
        int32_t minColIndex = colIndex - 1;
        while (minColIndex >= 0) {
          int64_t thisBC = vFoundGlobalBC[minColIndex];
          // check if this is still the same TF
          int64_t thisTFid = (thisBC - bcSOR) / nBCsPerTF;
          if (thisTFid != tfId)
            break;
          int64_t thisRofId = (thisBC + 3564 - rofOffset) / rofLength;

          // check if we are within the same ROF
          if (thisRofId != rofId)
            break;
          vAssocToSameROF.push_back(minColIndex);
          minColIndex--;
        }
        // find all collisions in the same ROF after the current one
        int32_t maxColIndex = colIndex + 1;
        while (maxColIndex < cols.size()) {
          int64_t thisBC = vFoundGlobalBC[maxColIndex];
          int64_t thisTFid = (thisBC - bcSOR) / nBCsPerTF;
          if (thisTFid != tfId)
            break;
          int64_t thisRofId = (thisBC + 3564 - rofOffset) / rofLength;
          if (thisRofId != rofId)
            break;
          vAssocToSameROF.push_back(maxColIndex);
          maxColIndex++;
        }
        vCollsInSameITSROF.push_back(vAssocToSameROF);
      } // end of in-ROF occupancy 1st loop

      // nTrack correlations in ROFs with 2 collisions inside
      for (const auto& col : cols) {
        int32_t colIndex = col.globalIndex();
        if (!col.sel8() || !col.selection_bit(kNoSameBunchPileup))
          continue;
        if (vCollsInSameITSROF[colIndex].size() != 1) // analyse only cases with 2 collisions in the same ROF
          continue;
        float vZ = col.posZ();
        float nPV = vTracksITS567perColl[colIndex];

        ushort flags = col.flags();
        bool isVertexUPC = flags & dataformats::Vertex<o2::dataformats::TimeStamp<int>>::Flags::UPCMode; // is vertex with UPC settings

        // the second collision in ROF
        std::vector<int> vAssocToSameROF = vCollsInSameITSROF[colIndex];
        int thisColIndex = vAssocToSameROF[0];
        float vZassoc = vCollVz[thisColIndex];               // vZ of the second collision in the same ROF
        float nPVassoc = vTracksITS567perColl[thisColIndex]; // n PV tracks of the second collision in the same ROF
        if (std::fabs(vZ) < 10 && std::fabs(vZassoc) < 10 && thisColIndex > colIndex && vIsSel8[thisColIndex]) {
          histos.fill(HIST("occupancyQA/hITSTracks_ev1_vs_ev2_2coll_in_ROF"), nPV, nPVassoc);
          if (isVertexUPC)
            histos.fill(HIST("occupancyQA/hITSTracks_ev1_vs_ev2_2coll_in_ROF_UPC"), nPV, nPVassoc);
          else
            histos.fill(HIST("occupancyQA/hITSTracks_ev1_vs_ev2_2coll_in_ROF_nonUPC"), nPV, nPVassoc);
        }
      }
    } // end of in-ROF occupancy QA

    // TVX efficiency after TF and ITS ROF border cuts
    for (const auto& col : cols) {
      if (!col.selection_bit(kNoTimeFrameBorder) || !col.selection_bit(kNoITSROFrameBorder))
        continue;

      uint32_t nContrib = col.numContrib();
      histos.fill(HIST("hNcontribColFromData"), nContrib);
      if (!col.selection_bit(kIsTriggerTVX))
        continue;

      histos.fill(HIST("hNcontribAccFromData"), nContrib);
    }
  }
  PROCESS_SWITCH(EventSelectionQaTask, processRun3, "Process Run3 event selection QA", false);

  Partition<FullTracksIUwithLabels> pvTracks = ((aod::track::flags & static_cast<uint32_t>(o2::aod::track::PVContributor)) == static_cast<uint32_t>(o2::aod::track::PVContributor));
  void processMCRun3(aod::McCollisions const& mcCols, soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels> const& cols, FullTracksIUwithLabels const&, BCsRun3 const&, aod::FT0s const&, aod::McParticles const& mcParts)
  {
    for (const auto& mcCol : mcCols) {
      auto bc = mcCol.bc_as<BCsRun3>();
      uint64_t globalBC = bc.globalBC();
      uint64_t orbit = globalBC / nBCsPerOrbit;
      int localBC = globalBC % nBCsPerOrbit;
      int64_t bcInTF = (globalBC - bcSOR) % nBCsPerTF;
      histos.fill(HIST("hGlobalBcColMC"), globalBC - bcSOR);
      histos.fill(HIST("hOrbitColMC"), orbit - orbitSOR);
      histos.fill(HIST("hBcColMC"), localBC);
      histos.fill(HIST("hVertexXMC"), mcCol.posX());
      histos.fill(HIST("hVertexYMC"), mcCol.posY());
      histos.fill(HIST("hVertexZMC"), mcCol.posZ());
      histos.fill(HIST("hNcolMCVsBcInTF"), bcInTF);
    }

    for (const auto& col : cols) {
      int32_t mcColIdFromCollision = col.mcCollisionId();
      // check if collision is built from tracks originating from different MC collisions
      bool isCollisionAmbiguous = 0;
      const auto& colPvTracks = pvTracks.sliceByCached(aod::track::collisionId, col.globalIndex(), cache);
      for (const auto& track : colPvTracks) {
        int32_t mcPartId = track.mcParticleId();
        int32_t mcColId = mcPartId >= 0 ? mcParts.iteratorAt(mcPartId).mcCollisionId() : -1;
        if (mcColId < 0 || mcColIdFromCollision != mcColId) {
          isCollisionAmbiguous = 1;
          break;
        }
      }

      // skip ambiguous collisions
      if (isCollisionAmbiguous)
        continue;

      // skip collisions at the borders of TF and ITS ROF
      if (!col.selection_bit(kNoTimeFrameBorder) || !col.selection_bit(kNoITSROFrameBorder))
        continue;

      uint32_t nContrib = col.numContrib();
      histos.fill(HIST("hNcontribColFromMC"), nContrib);
      if (!col.selection_bit(kIsTriggerTVX))
        continue;

      histos.fill(HIST("hNcontribAccFromMC"), nContrib);

      int64_t rcBC = col.foundBC_as<BCsRun3>().globalBC();
      int64_t mcBC = col.mcCollision().bc_as<BCsRun3>().globalBC();

      if (mcBC != rcBC) {
        histos.fill(HIST("hNcontribMisFromMC"), nContrib);
      }
    }
  }
  PROCESS_SWITCH(EventSelectionQaTask, processMCRun3, "Process Run3 MC event selection QA", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<EventSelectionQaTask>(cfgc)};
}
