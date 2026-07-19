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
/// \file upcVmRof.cxx
/// \brief analysis of UPC vector meson production ROF by ROF
///
/// \author Guillermo Contreras (jesus.guillermo.contreras.nuno@cern.ch), Czech Technical University in Prague

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/LHCConstants.h>
#include <CommonDataFormat/BunchFilling.h>
#include <CommonDataFormat/TimeStamp.h>
#include <DataFormatsParameters/AggregatedRunInfo.h>
#include <DataFormatsParameters/GRPLHCIFData.h>
#include <DataFormatsITSMFT/DPLAlpideParam.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <TH1.h>
#include <TH2.h>
#include <TString.h>

#include <bitset>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using BCsTSsSels = soa::Join<aod::BCsWithTimestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;
using ColSels = soa::Join<aod::Collisions, aod::EvSels>;
using TRKs = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection,
                       aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                       aod::pidTPCFullDe, aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl>;

using ColSel = ColSels::iterator;

namespace o2::aod
{
namespace datarows
{
// collision info
DECLARE_SOA_COLUMN(RunNumber, runNumber, int32_t);
DECLARE_SOA_COLUMN(PosX, posX, float);
DECLARE_SOA_COLUMN(PosY, posY, float);
DECLARE_SOA_COLUMN(PosZ, posZ, float);
DECLARE_SOA_COLUMN(Chi2, chi2, float);
DECLARE_SOA_COLUMN(LocalBC, localBC, int);
DECLARE_SOA_COLUMN(LocalTF, localTF, int);
DECLARE_SOA_COLUMN(LocalROF, localROF, int);
DECLARE_SOA_COLUMN(UpcFlag, upcFlag, int);

// FIT info
DECLARE_SOA_COLUMN(AmplitudeFT0A, amplitudeFT0A, float);
DECLARE_SOA_COLUMN(AmplitudeFT0C, amplitudeFT0C, float);
DECLARE_SOA_COLUMN(AmplitudeFV0A, amplitudeFV0A, float);
DECLARE_SOA_COLUMN(AmplitudeFDDA, amplitudeFDDA, float);
DECLARE_SOA_COLUMN(AmplitudeFDDC, amplitudeFDDC, float);
DECLARE_SOA_COLUMN(TimeFT0A, timeFT0A, float);
DECLARE_SOA_COLUMN(TimeFT0C, timeFT0C, float);
DECLARE_SOA_COLUMN(TimeFV0A, timeFV0A, float);
DECLARE_SOA_COLUMN(TimeFDDA, timeFDDA, float);
DECLARE_SOA_COLUMN(TimeFDDC, timeFDDC, float);
DECLARE_SOA_COLUMN(ChannelsFT0A, channelsFT0A, int);
DECLARE_SOA_COLUMN(ChannelsFT0C, channelsFT0C, int);
DECLARE_SOA_COLUMN(ChannelsFV0A, channelsFV0A, int);
DECLARE_SOA_COLUMN(ChannelsFDDA, channelsFDDA, int);
DECLARE_SOA_COLUMN(ChannelsFDDC, channelsFDDC, int);

// ZDC info
DECLARE_SOA_COLUMN(EnergyCommonZNA, energyCommonZNA, float);
DECLARE_SOA_COLUMN(EnergyCommonZNC, energyCommonZNC, float);
DECLARE_SOA_COLUMN(TimeZNA, timeZNA, float);
DECLARE_SOA_COLUMN(TimeZNC, timeZNC, float);

// track info
DECLARE_SOA_COLUMN(Pt1, pt1, float);
DECLARE_SOA_COLUMN(Eta1, eta1, float);
DECLARE_SOA_COLUMN(Phi1, phi1, float);
DECLARE_SOA_COLUMN(Q1, q1, int);
DECLARE_SOA_COLUMN(PidPion1, pidPion1, float);
DECLARE_SOA_COLUMN(PidElectron1, pidElectron1, float);
DECLARE_SOA_COLUMN(PidKaon1, pidKaon1, float);
DECLARE_SOA_COLUMN(PidProton1, pidProton1, float);
DECLARE_SOA_COLUMN(Pt2, pt2, float);
DECLARE_SOA_COLUMN(Eta2, eta2, float);
DECLARE_SOA_COLUMN(Phi2, phi2, float);
DECLARE_SOA_COLUMN(Q2, q2, int);
DECLARE_SOA_COLUMN(PidPion2, pidPion2, float);
DECLARE_SOA_COLUMN(PidElectron2, pidElectron2, float);
DECLARE_SOA_COLUMN(PidKaon2, pidKaon2, float);
DECLARE_SOA_COLUMN(PidProton2, pidProton2, float);
DECLARE_SOA_COLUMN(Pt3, pt3, float);
DECLARE_SOA_COLUMN(Eta3, eta3, float);
DECLARE_SOA_COLUMN(Phi3, phi3, float);
DECLARE_SOA_COLUMN(Q3, q3, int);
DECLARE_SOA_COLUMN(PidPion3, pidPion3, float);
DECLARE_SOA_COLUMN(PidElectron3, pidElectron3, float);
DECLARE_SOA_COLUMN(PidKaon3, pidKaon3, float);
DECLARE_SOA_COLUMN(PidProton3, pidProton3, float);
DECLARE_SOA_COLUMN(Pt4, pt4, float);
DECLARE_SOA_COLUMN(Eta4, eta4, float);
DECLARE_SOA_COLUMN(Phi4, phi4, float);
DECLARE_SOA_COLUMN(Q4, q4, int);
DECLARE_SOA_COLUMN(PidPion4, pidPion4, float);
DECLARE_SOA_COLUMN(PidElectron4, pidElectron4, float);
DECLARE_SOA_COLUMN(PidKaon4, pidKaon4, float);
DECLARE_SOA_COLUMN(PidProton4, pidProton4, float);
} // namespace datarows

DECLARE_SOA_TABLE(TwoTrkTable, "AOD", "TWOTRKTABLE",
                  datarows::RunNumber, datarows::PosX, datarows::PosY, datarows::PosZ, datarows::Chi2,
                  datarows::LocalBC, datarows::LocalTF, datarows::LocalROF, datarows::UpcFlag,
                  datarows::AmplitudeFT0A, datarows::AmplitudeFT0C, datarows::AmplitudeFV0A, datarows::AmplitudeFDDA, datarows::AmplitudeFDDC,
                  datarows::TimeFT0A, datarows::TimeFT0C, datarows::TimeFV0A, datarows::TimeFDDA, datarows::TimeFDDC,
                  datarows::ChannelsFT0A, datarows::ChannelsFT0C, datarows::ChannelsFV0A, datarows::ChannelsFDDA, datarows::ChannelsFDDC,
                  datarows::EnergyCommonZNA, datarows::EnergyCommonZNC, datarows::TimeZNA, datarows::TimeZNC,
                  datarows::Pt1, datarows::Eta1, datarows::Phi1, datarows::Q1, datarows::PidPion1, datarows::PidElectron1, datarows::PidKaon1, datarows::PidProton1,
                  datarows::Pt2, datarows::Eta2, datarows::Phi2, datarows::Q2, datarows::PidPion2, datarows::PidElectron2, datarows::PidKaon2, datarows::PidProton2);
DECLARE_SOA_TABLE(FourTrkTable, "AOD", "FOURTRKTABLE",
                  datarows::RunNumber, datarows::PosX, datarows::PosY, datarows::PosZ, datarows::Chi2,
                  datarows::LocalBC, datarows::LocalTF, datarows::LocalROF, datarows::UpcFlag,
                  datarows::AmplitudeFT0A, datarows::AmplitudeFT0C, datarows::AmplitudeFV0A, datarows::AmplitudeFDDA, datarows::AmplitudeFDDC,
                  datarows::TimeFT0A, datarows::TimeFT0C, datarows::TimeFV0A, datarows::TimeFDDA, datarows::TimeFDDC,
                  datarows::ChannelsFT0A, datarows::ChannelsFT0C, datarows::ChannelsFV0A, datarows::ChannelsFDDA, datarows::ChannelsFDDC,
                  datarows::EnergyCommonZNA, datarows::EnergyCommonZNC, datarows::TimeZNA, datarows::TimeZNC,
                  datarows::Pt1, datarows::Eta1, datarows::Phi1, datarows::Q1, datarows::PidPion1, datarows::PidElectron1, datarows::PidKaon1, datarows::PidProton1,
                  datarows::Pt2, datarows::Eta2, datarows::Phi2, datarows::Q2, datarows::PidPion2, datarows::PidElectron2, datarows::PidKaon2, datarows::PidProton2,
                  datarows::Pt3, datarows::Eta3, datarows::Phi3, datarows::Q3, datarows::PidPion3, datarows::PidElectron3, datarows::PidKaon3, datarows::PidProton3,
                  datarows::Pt4, datarows::Eta4, datarows::Phi4, datarows::Q4, datarows::PidPion4, datarows::PidElectron4, datarows::PidKaon4, datarows::PidProton4);
} // namespace o2::aod

struct UpcVmRof {

  // output
  Produces<o2::aod::TwoTrkTable> twoTrkTable;
  Produces<o2::aod::FourTrkTable> fourTrkTable;

  // services
  Service<o2::ccdb::BasicCCDBManager> ccdb{}; // access to database

  // histograms
  HistogramRegistry bcTH1Registry{"bcTH1Registry", {}};
  std::map<std::string, std::shared_ptr<TH1>> bcTH1Pointers;
  HistogramRegistry bcTH2Registry{"bcTH2Registry", {}};
  std::map<std::string, std::shared_ptr<TH2>> bcTH2Pointers;
  HistogramRegistry colTH1Registry{"colTH1Registry", {}};
  std::map<std::string, std::shared_ptr<TH1>> colTH1Pointers;

  // variables to store filling scheme info
  std::bitset<o2::constants::lhc::LHCMaxBunches> beamPatternA;
  std::bitset<o2::constants::lhc::LHCMaxBunches> beamPatternC;
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternA;
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternB;
  std::bitset<o2::constants::lhc::LHCMaxBunches> bcPatternC;
  std::vector<int> bcbIdx;
  int nbcB = 0;

  // variables to store ITS ROF info
  int rofPerOrbit = -1;    // number of rofs per orbit
  int rofLength = -1;    // number of bcs per ROF
  int rofShift = -1;    // bc shift of ITS.

  // variables to store run info
  int runNumberBc = 0;     // run number used to process BCs
  int runNumberCol = 0;    // run number used to process collisions
  int64_t sor = 0;         // best known timestamp for the start of run
  int64_t orbitsPerTF = 0; // number of orbits per TF
  int64_t bcSOR = 0;       // first bc of the first orbit
  int64_t nBCsPerTF = 0;   // duration of TF in bcs
  int64_t currentTF = -1;  // current time frame being looked at
  int64_t nTF = 0;         // number of time frames in run

  // constant related to trigger mask indices
  // https://github.com/AliceO2Group/AliceO2/blob/6c0251c35e5cbf6028d2bd6f7e30f9a4fd348e38/DataFormats/Detectors/CTP/src/Configuration.cxx#L1123
  // PbPb triggers: 1ZNC, FV0CH+FT0VTX, OO: 1ZNC, FT0CE+FT0VTX
  static constexpr int Ft0VtxIdx = 2;
  static constexpr int Ft0CeIdx = 4;

  // number of tracks for the selected event topologies
  static constexpr int NTrksTwoBody = 2;
  static constexpr int NTrksFourBody = 4;

  // reconstruction modes
  static constexpr int stdReco = 0;
  static constexpr int upcReco = 1;
  
  // information for selection collisions
  Configurable<float> maxAbsPosZ{"maxAbsPosZ", 10.0, "max |Z| position of vtx"};
  Configurable<float> maxAbsTimeFT0{"maxAbsTimeFT0", 4.0, "max |time| in a FT0 side"};
  Configurable<float> maxAmpFT0{"maxAmpFT0", 150.0, "max amplitude for signals in a FT0 side"};
  Configurable<float> maxTrkTpcChi2{"maxTrkTpcChi2", 4.0, "max chi2 of TPC track"};
  Configurable<float> maxTrkItsChi2{"maxTrkItsChi2", 36.0, "max chi2 of ITS track"};
  Configurable<float> minTrkTpcClusters{"minTrkTpcClusters", 70.0, "minimum number of TPC clusters associated to the track"};
  Configurable<float> maxTrkDcaZ{"maxTrkDcaZ", 2.0, "max DCA in z of track to vtx (cm)"};
  Configurable<int> tfPerBin{"tfPerBin", 10000, "timeframes per bin 1e4 means some 28 s"};



  //--------------------------------------------------------------------------------
  // get ITS ROF info
  // code from https://github.com/AliceO2Group/O2Physics/blob/master/Common/Tools/EventSelectionModule.h#L779-L780
  void getRofInfo()
  {
    auto alppar = ccdb->getForTimeStamp<o2::itsmft::DPLAlpideParam<0>>("ITS/Config/AlpideParam", sor);
    rofShift = alppar->roFrameBiasInBC;
    rofLength = alppar->roFrameLengthInBC;
    rofPerOrbit = static_cast<int>(o2::constants::lhc::LHCMaxBunches / rofLength);
  }

  //--------------------------------------------------------------------------------
  // get filling scheme
  void getFillingScheme()
  {
    // get the info
    auto grplhcif = ccdb->getForTimeStamp<o2::parameters::GRPLHCIFData>("GLO/Config/GRPLHCIF", sor);
    beamPatternA = grplhcif->getBunchFilling().getBeamPattern(0);
    beamPatternC = grplhcif->getBunchFilling().getBeamPattern(1);
    bcPatternA = beamPatternA & ~beamPatternC;
    bcPatternB = beamPatternA & beamPatternC;
    bcPatternC = ~beamPatternA & beamPatternC;
    // define internal indices to access the correct bin in histogram
    bcbIdx.clear();
    nbcB = 0;
    for (int i = 0; i < o2::constants::lhc::LHCMaxBunches; i++) {
      bcbIdx.push_back(-1);
      if (bcPatternB.test(i)) {
        bcbIdx[i] = nbcB;
        nbcB++;
      }
    }
  } // end getFillingScheme()

  //--------------------------------------------------------------------------------
  // store bc patterns
  void fillBcPatternHistos(int run)
  {
    for (int i = 0; i < o2::constants::lhc::LHCMaxBunches; i++) {
      if (bcPatternA.test(i)) {
        bcTH1Pointers[Form("bc/%d/bcPatternA_H", run)]->Fill(i);
      }
      if (bcPatternB.test(i)) {
        bcTH1Pointers[Form("bc/%d/bcPatternB_H", run)]->Fill(i);
      }
      if (bcPatternC.test(i)) {
        bcTH1Pointers[Form("bc/%d/bcPatternC_H", run)]->Fill(i);
      }
    }
  } // end fillBcPatternHistos

  //--------------------------------------------------------------------------------
  // get general run information
  void getRunInfo(int run)
  {
    auto runInfo = o2::parameters::AggregatedRunInfo::buildAggregatedRunInfo(ccdb->instance(), run);
    sor = runInfo.sor; // in ms
    auto orbitSOR = runInfo.orbitSOR;
    auto orbitEOR = runInfo.orbitEOR;
    orbitsPerTF = runInfo.orbitsPerTF;
    bcSOR = orbitSOR * o2::constants::lhc::LHCMaxBunches;        // first bc of the first orbit
    nBCsPerTF = orbitsPerTF * o2::constants::lhc::LHCMaxBunches; // duration of TF in bcs
    nTF = std::ceil((orbitEOR - orbitSOR) / orbitsPerTF);
  } // end getRunInfo()

  //--------------------------------------------------------------------------------
  // compute  Bc within the orbit
  int64_t getBcWithinOrbit(int64_t globalBC)
  {
    return (globalBC % o2::constants::lhc::LHCMaxBunches);
  }

  //--------------------------------------------------------------------------------
  // compute TF for this BC
  int64_t getTimeFrame(int64_t globalBC)
  {
    return (globalBC - bcSOR) / nBCsPerTF;
  }

  //--------------------------------------------------------------------------------
  // compute ROF for this BC
  int getRof(int64_t thisBC)
  {
    int64_t bctmp = thisBC - rofShift;
    if (bctmp < 0) {
      bctmp = o2::constants::lhc::LHCMaxBunches - rofShift - 1;
    }
    return static_cast<int>(bctmp / rofLength);
  }

  //--------------------------------------------------------------------------------
  // check flags for a bc
  bool checkBcFlags(const auto& bc, int run)
  {
    bcTH1Pointers[Form("bc/%d/bcSel_H", run)]->Fill(0);
    if (!bc.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    bcTH1Pointers[Form("bc/%d/bcSel_H", run)]->Fill(1);
    if (!bc.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    bcTH1Pointers[Form("bc/%d/bcSel_H", run)]->Fill(2);
    return true;
  } // end checkBcFlags

  //--------------------------------------------------------------------------------
  // check flags for a colllision
  bool checkColFlags(const auto& col, int run)
  {
    colTH1Pointers[Form("col/%d/colSel_H", run)]->Fill(0);
    if (!col.has_foundBC()) {
      return false;
    }
    colTH1Pointers[Form("col/%d/colSel_H", run)]->Fill(1);
    if (!col.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    colTH1Pointers[Form("col/%d/colSel_H", run)]->Fill(2);
    if (!col.selection_bit(aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    colTH1Pointers[Form("col/%d/colSel_H", run)]->Fill(3);
    if (!col.selection_bit(aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    colTH1Pointers[Form("col/%d/colSel_H", run)]->Fill(4);
    if (!col.selection_bit(aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    colTH1Pointers[Form("col/%d/colSel_H", run)]->Fill(5);
    // missing: add rct flags
    return true;
  } // end checkColFlags

  //--------------------------------------------------------------------------------
  // add histos used in processBCs
  void addBcHistos(int run)
  {
    // check if histos already created
    if (bcTH1Pointers[Form("bc/%d/bcPatternA_H", run)]) {
      return;
    }

    // filling scheme
    bcTH1Pointers[Form("bc/%d/bcPatternA_H", run)] = bcTH1Registry.add<TH1>(Form("bc/%d/bcPatternA_H", run), "Pattern of bc-A; bcID;",
                                                                            {HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, static_cast<double>(o2::constants::lhc::LHCMaxBunches) - 0.5}}});
    bcTH1Pointers[Form("bc/%d/bcPatternB_H", run)] = bcTH1Registry.add<TH1>(Form("bc/%d/bcPatternB_H", run), "Pattern of bc-B; bcID;",
                                                                            {HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, static_cast<double>(o2::constants::lhc::LHCMaxBunches) - 0.5}}});
    bcTH1Pointers[Form("bc/%d/bcPatternC_H", run)] = bcTH1Registry.add<TH1>(Form("bc/%d/bcPatternC_H", run), "Pattern of bc-C; bcID;",
                                                                            {HistType::kTH1D, {{o2::constants::lhc::LHCMaxBunches, -0.5, static_cast<double>(o2::constants::lhc::LHCMaxBunches) - 0.5}}});

    // bc sel and tf info
    int nBinsTF = static_cast<int>(nTF / tfPerBin) + 1;
    int lastTFinHisto = (nBinsTF * tfPerBin) - 1; // first TF is zero
    bcTH1Pointers[Form("bc/%d/bcSel_H", run)] = bcTH1Registry.add<TH1>(Form("bc/%d/bcSel_H", run), "bc selection counter; selID; Counter",
                                                                       {HistType::kTH1D, {{4, -0.5, 3.5}}});
    bcTH1Pointers[Form("bc/%d/tf_H", run)] = bcTH1Registry.add<TH1>(Form("bc/%d/tf_H", run), "analysed time frames;TF;Counts",
                                                                    {HistType::kTH1D, {{nBinsTF, -0.5, static_cast<double>(lastTFinHisto) - 0.5}}});
    // trigger info per rof
    bcTH2Pointers[Form("bc/%d/ft0Vtx_H", run)] = bcTH2Registry.add<TH2>(Form("bc/%d/ft0Vtx_H", run), "ft0Vtx triggers; TF; ROF; Counter",
                                                                        {HistType::kTH2F, {{nBinsTF, -0.5, static_cast<double>(lastTFinHisto) - 0.5}, {rofPerOrbit, -0.5, rofPerOrbit - 0.5}}});
    bcTH2Pointers[Form("bc/%d/ft0VtxCe_H", run)] = bcTH2Registry.add<TH2>(Form("bc/%d/ft0VtxCe_H", run), "ft0VtxCe triggers; TF; ROF; Counter",
                                                                          {HistType::kTH2F, {{nBinsTF, -0.5, static_cast<double>(lastTFinHisto) - 0.5}, {rofPerOrbit, -0.5, rofPerOrbit - 0.5}}});
    // trigger info per bcb
    bcTH2Pointers[Form("bc/%d/ft0Vtx_bcb_H", run)] = bcTH2Registry.add<TH2>(Form("bc/%d/ft0Vtx_bcb_H", run), "ft0Vtx triggers; TF; bc-B idx; Counter",
                                                                            {HistType::kTH2F, {{nBinsTF, -0.5, static_cast<double>(lastTFinHisto) - 0.5}, {nbcB, -0.5, nbcB - 0.5}}});
    bcTH2Pointers[Form("bc/%d/ft0VtxCe_bcb_H", run)] = bcTH2Registry.add<TH2>(Form("bc/%d/ft0VtxCe_bcb_H", run), "ft0Vtx triggers; TF; bc-B idx; Counter",
                                                                              {HistType::kTH2F, {{nBinsTF, -0.5, static_cast<double>(lastTFinHisto) - 0.5}, {nbcB, -0.5, nbcB - 0.5}}});
  } // addBcHistos

  //--------------------------------------------------------------------------------
  // add histos used in processBCs
  void addColHistos(int run)
  {
    // check if histos already created
    if (colTH1Pointers[Form("col/%d/colSel_H", run)]) {
      return;
    }

    // counters
    colTH1Pointers[Form("col/%d/colSel_H", run)] = colTH1Registry.add<TH1>(Form("col/%d/colSel_H", run),
                                                                           "collision selection counter; selID; Counter",
                                                                           {HistType::kTH1D, {{25, -0.5, 24.5}}});
    colTH1Pointers[Form("col/%d/trkSel_H", run)] = colTH1Registry.add<TH1>(Form("col/%d/trkSel_H", run),
                                                                           "track selection counter; selID; Counter",
                                                                           {HistType::kTH1D, {{10, -0.5, 9.5}}});
    // track properties
    colTH1Pointers[Form("col/%d/tpcChi2_H", run)] = colTH1Registry.add<TH1>(Form("col/%d/tpcChi2_H", run),
                                                                            "track #chi^2 in TPC; #chi^2; Entries",
                                                                            {HistType::kTH1D, {{100, 0, 5}}});
    colTH1Pointers[Form("col/%d/itsChi2_H", run)] = colTH1Registry.add<TH1>(Form("col/%d/itsChi2_H", run),
                                                                            "track #chi^2 in ITS; #chi^2; Entries",
                                                                            {HistType::kTH1D, {{120, 0, 40}}});
    colTH1Pointers[Form("col/%d/tpcNcls_H", run)] = colTH1Registry.add<TH1>(Form("col/%d/tpcNcls_H", run),
                                                                            "track clusters in TPC; n-clusters; Entries",
                                                                            {HistType::kTH1D, {{160, 0, 160}}});
    colTH1Pointers[Form("col/%d/tpcXoF_H", run)] = colTH1Registry.add<TH1>(Form("col/%d/tpcXoF_H", run),
                                                                           "track crossed rows over findable TPC clusters; XoF; Entries",
                                                                           {HistType::kTH1D, {{150, 0.5, 2.0}}});
    colTH1Pointers[Form("col/%d/trkDcaZ_H", run)] = colTH1Registry.add<TH1>(Form("col/%d/trkDcaZ_H", run),
                                                                            "track z-DCA; z-dca (cm); Entries",
                                                                            {HistType::kTH1D, {{100, -0.5, 0.5}}});
    colTH1Pointers[Form("col/%d/trkDcaXY_H", run)] = colTH1Registry.add<TH1>(Form("col/%d/trkDcaXY_H", run),
                                                                             "track xy-DCA; xy-dca (cm); Entries",
                                                                             {HistType::kTH1D, {{100, -0.1, 0.1}}});
    // selected events per TF
    colTH1Pointers[Form("col/%d/twoTrkTF_H", run)] = colTH1Registry.add<TH1>(Form("col/%d/twoTrkTF_H", run), "Number of selected 2-trk events per;TF;Entries",
                                                                             {HistType::kTH1D, {{static_cast<int>(nTF / tfPerBin), -0.5, static_cast<double>(nTF) - 0.5}}});
    colTH1Pointers[Form("col/%d/fourTrkTF_H", run)] = colTH1Registry.add<TH1>(Form("col/%d/fourTrkTF_H", run), "Number of selected 4-trk events per;TF;Entries",
                                                                              {HistType::kTH1D, {{static_cast<int>(nTF / tfPerBin), -0.5, static_cast<double>(nTF) - 0.5}}});

  } // addColHistos

  //--------------------------------------------------------------------------------
  // initialization
  void init(InitContext&)
  {
    // set access to db
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
  } // end init()

  //--------------------------------------------------------------------------------
  // process BCs to get trigger information
  void processBCs(BCsTSsSels const& bcs)
  {
    // check if new run
    // -- assuming that each call to processBCs involves only one run
    const auto& bcAt0 = bcs.iteratorAt(0);
    if (runNumberBc != bcAt0.runNumber()) { // new run
      runNumberBc = bcAt0.runNumber();
      getRunInfo(runNumberBc);
      getFillingScheme();
      getRofInfo();
      addBcHistos(runNumberBc);
      fillBcPatternHistos(runNumberBc);
    }

    //--------------------------------------------------------------------------------
    for (const auto& bc : bcs) {
      // get info for this bc
      int64_t thisBC = getBcWithinOrbit(bc.globalBC());
      int64_t thisTF = getTimeFrame(bc.globalBC());
      int thisROF = getRof(thisBC);

      // count TF seen
      if (thisTF != currentTF) {
        currentTF = thisTF;
        bcTH1Pointers[Form("bc/%d/tf_H", runNumberBc)]->Fill(thisTF);
      }

      // check that the bc pass the selection
      if (!checkBcFlags(bc, runNumberBc)) {
        continue;
      }

      // consider only b-bcs
      if (!bcPatternB.test(thisBC)) {
        continue;
      }
      bcTH1Pointers[Form("bc/%d/bcSel_H", runNumberBc)]->Fill(3);

      // get triggers
      std::bitset<64> mask = bc.inputMask();
      bool ft0vtxTrg = mask[Ft0VtxIdx];
      bool ft0ceTrg = mask[Ft0CeIdx];
      if (ft0vtxTrg) {
        bcTH2Pointers[Form("bc/%d/ft0Vtx_H", runNumberBc)]->Fill(thisTF, thisROF);
        bcTH2Pointers[Form("bc/%d/ft0Vtx_bcb_H", runNumberBc)]->Fill(thisTF, bcbIdx[thisBC]);
        if (ft0ceTrg) {
          bcTH2Pointers[Form("bc/%d/ft0VtxCe_H", runNumberBc)]->Fill(thisTF, thisROF);
          bcTH2Pointers[Form("bc/%d/ft0VtxCe_bcb_H", runNumberBc)]->Fill(thisTF, bcbIdx[thisBC]);
        }
      }
    } // loop over bcs

    // reset current time frame
    currentTF = -1;
  } // end processBCs
  PROCESS_SWITCH(UpcVmRof, processBCs, "get BCs and trigger information", true);

  //--------------------------------------------------------------------------------
  // get collision information
  void processCols(ColSel const& col, BCsTSsSels const&, TRKs const& tracks,
                   aod::FV0As const&, aod::FT0s const&, aod::FDDs const&,
                   aod::Zdcs const&)
  {
    // get info for this bc
    auto bc = col.template foundBC_as<BCsTSsSels>();
    if (runNumberCol != bc.runNumber()) { // new run
      runNumberCol = bc.runNumber();
      getRunInfo(runNumberCol);
      getFillingScheme();
      addColHistos(runNumberCol);
    }
    int64_t thisBC = getBcWithinOrbit(bc.globalBC());
    int64_t thisTF = getTimeFrame(bc.globalBC());
    int64_t thisROF = getRof(thisBC);

    // select collision
    if (!checkColFlags(col, runNumberCol)) {
      return;
    }

    // accept only -B bcs
    if (!bcPatternB.test(thisBC)) {
      return;
    }
    colTH1Pointers[Form("col/%d/colSel_H", runNumberCol)]->Fill(10);

    // select on zVtx
    if (std::abs(col.posZ()) > maxAbsPosZ) {
      return;
    }
    colTH1Pointers[Form("col/%d/colSel_H", runNumberCol)]->Fill(11);

    // select number of contributors
    bool isTwoContributors = (col.numContrib() == NTrksTwoBody);
    bool isFourContributors = (col.numContrib() == NTrksFourBody);
    if ( !isTwoContributors && !isFourContributors) {
      return;
    }
    if (isTwoContributors) {
      colTH1Pointers[Form("col/%d/colSel_H", runNumberCol)]->Fill(16);
    }
    if (isFourContributors) {
      colTH1Pointers[Form("col/%d/colSel_H", runNumberCol)]->Fill(18);
    }

    // select tracks
    std::vector<decltype(tracks.begin())> selTrks;
    colTH1Pointers[Form("col/%d/trkSel_H", runNumberCol)]->Fill(0);
    for (const auto& track : tracks) {
      if (!track.isPVContributor()) {
        continue;
      }
      colTH1Pointers[Form("col/%d/trkSel_H", runNumberCol)]->Fill(1);
      if (!track.hasITS()) {
        continue;
      }
      colTH1Pointers[Form("col/%d/trkSel_H", runNumberCol)]->Fill(2);
      if (!track.hasTPC()) {
        continue;
      }
      colTH1Pointers[Form("col/%d/trkSel_H", runNumberCol)]->Fill(3);
      colTH1Pointers[Form("col/%d/tpcChi2_H", runNumberCol)]->Fill(track.tpcChi2NCl());
      colTH1Pointers[Form("col/%d/itsChi2_H", runNumberCol)]->Fill(track.itsChi2NCl());
      colTH1Pointers[Form("col/%d/tpcNcls_H", runNumberCol)]->Fill(track.tpcNClsFindable());
      colTH1Pointers[Form("col/%d/tpcXoF_H", runNumberCol)]->Fill(track.tpcCrossedRowsOverFindableCls());
      colTH1Pointers[Form("col/%d/trkDcaZ_H", runNumberCol)]->Fill(track.dcaZ());
      colTH1Pointers[Form("col/%d/trkDcaXY_H", runNumberCol)]->Fill(track.dcaXY());
      if (track.tpcChi2NCl() > maxTrkTpcChi2) {
        continue;
      }
      colTH1Pointers[Form("col/%d/trkSel_H", runNumberCol)]->Fill(4);
      if (track.itsChi2NCl() > maxTrkItsChi2) {
        continue;
      }
      colTH1Pointers[Form("col/%d/trkSel_H", runNumberCol)]->Fill(5);
      if (track.tpcNClsFindable() < minTrkTpcClusters) {
        continue;
      }
      colTH1Pointers[Form("col/%d/trkSel_H", runNumberCol)]->Fill(6);
      if (std::abs(track.dcaZ()) > maxTrkDcaZ) {
        continue;
      }
      colTH1Pointers[Form("col/%d/trkSel_H", runNumberCol)]->Fill(7);
      float maxTrkDcaXY = 0.0105 + 0.035 / std::pow(track.pt(), 1.1);
      if (std::abs(track.dcaXY()) > maxTrkDcaXY) {
        continue;
      }
      colTH1Pointers[Form("col/%d/trkSel_H", runNumberCol)]->Fill(8);
      selTrks.push_back(track);
    }
    bool isTwoBody = isTwoContributors && (selTrks.size() == NTrksTwoBody);
    bool isFourBody = isFourContributors && (selTrks.size() == NTrksFourBody);    
    if (!isTwoBody && !isFourBody) {
      return;
    }
    
    //  selected events
    if (isTwoBody) {
      colTH1Pointers[Form("col/%d/colSel_H", runNumberCol)]->Fill(17);
    }
    if (isFourBody) {
      colTH1Pointers[Form("col/%d/colSel_H", runNumberCol)]->Fill(19);
    }

    // FT0 selection
    float aFT0A = 0;
    float aFT0C = 0;
    float tFT0A = 33; // default time to mark events without FT0 info
    float tFT0C = 33; // default time to mark events without FT0 info
    int nFT0A = 0;
    int nFT0C = 0;
    if (bc.has_foundFT0()) {
      // a side
      if (bc.foundFT0().isValidTimeA()) { // valid time
        tFT0A = bc.foundFT0().timeA();
        if (std::abs(tFT0A) > maxAbsTimeFT0)
          return;
        colTH1Pointers[Form("col/%d/colSel_H", runNumberCol)]->Fill(12);
        aFT0A = bc.foundFT0().sumAmpA();
        if (aFT0A > maxAmpFT0)
          return;
        colTH1Pointers[Form("col/%d/colSel_H", runNumberCol)]->Fill(13);
        nFT0A = (bc.foundFT0().amplitudeA()).size();
      } // a side
      // c side
      if (bc.foundFT0().isValidTimeC()) { // valid time
        tFT0C = bc.foundFT0().timeC();
        if (std::abs(tFT0C) > maxAbsTimeFT0)
          return;
        colTH1Pointers[Form("col/%d/colSel_H", runNumberCol)]->Fill(14);
        aFT0C = bc.foundFT0().sumAmpC();
        if (aFT0C > maxAmpFT0)
          return;
        colTH1Pointers[Form("col/%d/colSel_H", runNumberCol)]->Fill(15);
        nFT0C = (bc.foundFT0().amplitudeC()).size();
      } // c side
    } // FT0 selection

    // final number of selected events
    if (isTwoBody) {
      colTH1Pointers[Form("col/%d/colSel_H", runNumberCol)]->Fill(20);
    }
    if (isFourBody) {
      colTH1Pointers[Form("col/%d/colSel_H", runNumberCol)]->Fill(21);
    }

    // get FV0A info
    float aFV0A = 0;
    float tFV0A = 33; // default time to mark events without FV0 info
    int nFV0A = 0;
    if (bc.has_foundFV0()) {
      tFV0A = bc.foundFV0().time();
      auto v = bc.foundFV0().amplitude();
      aFV0A = std::accumulate(v.begin(), v.end(), 0.f);
      nFV0A = v.size();
    } // FV0A info

    // get FDD info
    float aFDDA = 0;
    float tFDDA = 33; // default time to mark events without FDD info
    int nFDDA = 0;
    float aFDDC = 0;
    float tFDDC = 33; // default time to mark events without FDD info
    int nFDDC = 0;
    if (bc.has_foundFDD()) {
      tFDDA = bc.foundFDD().timeA();
      auto vA = bc.foundFDD().chargeA();
      // channelPairs = {{0, 4}, {1, 5}, {2, 6}, {3, 7}};
      if (vA[0] > 0 && vA[4] > 0) {
        aFDDA += 0.5 * (vA[0] + vA[4]);
        nFDDA++;
      }
      if (vA[1] > 0 && vA[5] > 0) {
        aFDDA += 0.5 * (vA[1] + vA[5]);
        nFDDA++;
      }
      if (vA[2] > 0 && vA[6] > 0) {
        aFDDA += 0.5 * (vA[2] + vA[6]);
        nFDDA++;
      }
      if (vA[3] > 0 && vA[7] > 0) {
        aFDDA += 0.5 * (vA[3] + vA[7]);
        nFDDA++;
      }
      tFDDC = bc.foundFDD().timeC();
      auto vC = bc.foundFDD().chargeC();
      // channelPairs = {{0, 4}, {1, 5}, {2, 6}, {3, 7}};
      if (vC[0] > 0 && vC[4] > 0) {
        aFDDC += 0.5 * (vC[0] + vC[4]);
        nFDDC++;
      }
      if (vC[1] > 0 && vC[5] > 0) {
        aFDDC += 0.5 * (vC[1] + vC[5]);
        nFDDC++;
      }
      if (vC[2] > 0 && vC[6] > 0) {
        aFDDC += 0.5 * (vC[2] + vC[6]);
        nFDDC++;
      }
      if (vC[3] > 0 && vC[7] > 0) {
        aFDDC += 0.5 * (vC[3] + vC[7]);
        nFDDC++;
      }
    } // FDD info

    // get ZDC info
    float tZNA = -999; // default time to mark events without ZN info
    float tZNC = -999; // default time to mark events without ZN info
    float eZNA = -999;
    float eZNC = -999;
    if (bc.has_zdc()) {
      tZNA = (bc.zdc()).timeZNA();
      tZNC = (bc.zdc()).timeZNC();
      eZNA = (bc.zdc()).energyCommonZNA();
      eZNC = (bc.zdc()).energyCommonZNC();
      if (!std::isfinite(tZNA))
        tZNA = -999;
      if (!std::isfinite(tZNC))
        tZNC = -999;
      if (!std::isfinite(eZNA))
        eZNA = -999;
      if (!std::isfinite(eZNC))
        eZNC = -999;
    } // ZDC info

    // fill output table
    //  int recoFlag = TESTBIT(col.flags(), dataformats::Vertex<o2::dataformats::TimeStamp<int>>::Flags::UPCMode) ? upcReco : stdReco;
    const int recoFlag= (col.flags() & dataformats::Vertex<o2::dataformats::TimeStamp<int>>::Flags::UPCMode) ? upcReco : stdReco;
    if (isTwoBody) {
      colTH1Pointers[Form("col/%d/twoTrkTF_H", runNumberCol)]->Fill(thisTF);
      twoTrkTable(runNumberCol, col.posX(), col.posY(), col.posZ(), col.chi2(), thisBC, thisTF, thisROF, recoFlag,
                  aFT0A, aFT0C, aFV0A, aFDDA, aFDDC, tFT0A, tFT0C, tFV0A, tFDDA, tFDDC, nFT0A, nFT0C, nFV0A, nFDDA, nFDDC,
                  eZNA, eZNC, tZNA, tZNC,
                  selTrks[0].pt(), selTrks[0].eta(), selTrks[0].phi(), selTrks[0].sign(),
                  selTrks[0].tpcNSigmaPi(), selTrks[0].tpcNSigmaEl(), selTrks[0].tpcNSigmaKa(), selTrks[0].tpcNSigmaPr(),
                  selTrks[1].pt(), selTrks[1].eta(), selTrks[1].phi(), selTrks[1].sign(),
                  selTrks[1].tpcNSigmaPi(), selTrks[1].tpcNSigmaEl(), selTrks[1].tpcNSigmaKa(), selTrks[1].tpcNSigmaPr());
    }
    if (isFourBody) {
      colTH1Pointers[Form("col/%d/fourTrkTF_H", runNumberCol)]->Fill(thisTF);
      fourTrkTable(runNumberCol, col.posX(), col.posY(), col.posZ(), col.chi2(), thisBC, thisTF, thisROF, recoFlag,
                   aFT0A, aFT0C, aFV0A, aFDDA, aFDDC, tFT0A, tFT0C, tFV0A, tFDDA, tFDDC, nFT0A, nFT0C, nFV0A, nFDDA, nFDDC,
                   eZNA, eZNC, tZNA, tZNC,
                   selTrks[0].pt(), selTrks[0].eta(), selTrks[0].phi(), selTrks[0].sign(),
                   selTrks[0].tpcNSigmaPi(), selTrks[0].tpcNSigmaEl(), selTrks[0].tpcNSigmaKa(), selTrks[0].tpcNSigmaPr(),
                   selTrks[1].pt(), selTrks[1].eta(), selTrks[1].phi(), selTrks[1].sign(),
                   selTrks[1].tpcNSigmaPi(), selTrks[1].tpcNSigmaEl(), selTrks[1].tpcNSigmaKa(), selTrks[1].tpcNSigmaPr(),
                   selTrks[2].pt(), selTrks[2].eta(), selTrks[2].phi(), selTrks[2].sign(),
                   selTrks[2].tpcNSigmaPi(), selTrks[2].tpcNSigmaEl(), selTrks[2].tpcNSigmaKa(), selTrks[2].tpcNSigmaPr(),
                   selTrks[3].pt(), selTrks[3].eta(), selTrks[3].phi(), selTrks[3].sign(),
                   selTrks[3].tpcNSigmaPi(), selTrks[3].tpcNSigmaEl(), selTrks[3].tpcNSigmaKa(), selTrks[3].tpcNSigmaPr());
    }

  } // end processCol
  PROCESS_SWITCH(UpcVmRof, processCols, "get collisions and track information", true);

}; // end of struct UpcVmRof

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcVmRof>(cfgc)};
}
