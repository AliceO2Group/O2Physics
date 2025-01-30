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
/// \file upcTauCentralBarrelRL.cxx
/// \brief Personal task to analyze tau events from UPC collisions
///
/// \author Roman Lavicka <roman.lavicka@cern.ch>, Austrian Academy of Sciences & SMI
/// \since  12.07.2022
//

// C++ headers
#include <set>
#include <utility>
#include <algorithm>
#include <vector>

// O2 headers
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

// O2Physics headers
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGUD/Core/UPCTauCentralBarrelHelperRL.h"
#include "PWGUD/DataModel/UDTables.h"
#include "PWGUD/Core/SGSelector.h"

// ROOT headers
#include "TLorentzVector.h"
#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;

struct UpcTauRl {

  // Global varialbes
  bool isMC = false;
  Service<o2::framework::O2DatabasePDG> pdg;
  SGSelector sgSelector;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  // declare configurables
  Configurable<bool> verboseInfo{"verboseInfo", false, {"Print general info to terminal; default it false."}};
  Configurable<bool> doMainHistos{"doMainHistos", false, {"Fill main histos"}};
  Configurable<bool> doPIDhistos{"doPIDhistos", false, {"Fill PID histos"}};
  Configurable<bool> doTruthHistos{"doTruthHistos", false, {"Do histograms specific for generated events/particles"}};
  Configurable<bool> doMCtrueElectronCheck{"doMCtrueElectronCheck", false, {"Check if track hypothesis corresponds to MC truth. If no, it cuts."}};
  Configurable<bool> oppositeMCtrueElectronCheck{"oppositeMCtrueElectronCheck", false, {"While doMCtrueElectronCheck is true, check if track hypothesis corresponds to MC truth. If yes, it cuts."}};
  Configurable<bool> doTwoTracks{"doTwoTracks", false, {"Define histos for two tracks and allow to fill them"}};
  Configurable<bool> doFourTracks{"doFourTracks", false, {"Define histos for four tracks and allow to fill them"}};
  Configurable<bool> doSixTracks{"doSixTracks", false, {"Define histos for six tracks and allow to fill them"}};

  struct : ConfigurableGroup {
    Configurable<int> whichGapSide{"whichGapSide", 2, {"0 for side A, 1 for side C, 2 for both sides"}};
    Configurable<bool> useTrueGap{"useTrueGap", true, {"Calculate gapSide for a given FV0/FT0/ZDC thresholds"}};
    Configurable<float> cutTrueGapSideFV0{"cutTrueGapSideFV0", -1, "FV0A threshold for SG selector"};
    Configurable<float> cutTrueGapSideFT0A{"cutTrueGapSideFT0A", 150., "FT0A threshold for SG selector"};
    Configurable<float> cutTrueGapSideFT0C{"cutTrueGapSideFT0C", 50., "FT0C threshold for SG selector"};
    Configurable<float> cutTrueGapSideZDC{"cutTrueGapSideZDC", 0., "ZDC threshold for SG selector. 0 is <1n, 4.2 is <2n, 6.7 is <3n, 9.5 is <4n, 12.5 is <5n"};
    Configurable<float> cutFITtime{"cutFITtime", 40., "Maximum FIT time allowed. Default is 40ns"};
    Configurable<bool> applyAcceptanceSelection{"applyAcceptanceSelection", false, {"Select events in ALICE CB acceptance set with cutTrackEta"}};
    Configurable<float> cutTrackEta{"cutTrackEta", 0.9, "Cut on central barrel track eta in absolute values."};
  } cutSample;

  struct : ConfigurableGroup {
    Configurable<bool> applyGlobalTrackSelection{"applyGlobalTrackSelection", false, {"Applies cut on here defined global tracks"}};
    Configurable<float> cutMinPt{"cutMinPt", 0.1f, {"Global track cut"}};
    Configurable<float> cutMaxPt{"cutMaxPt", 1e10f, {"Global track cut"}};
    Configurable<float> cutMinEta{"cutMinEta", -0.8f, {"Global track cut"}};
    Configurable<float> cutMaxEta{"cutMaxEta", 0.8f, {"Global track cut"}};
    Configurable<float> cutMaxDCAz{"cutMaxDCAz", 2.f, {"Global track cut"}};
    Configurable<float> cutMaxDCAxy{"cutMaxDCAxy", 1e10f, {"Global track cut"}};
    Configurable<bool> applyPtDependentDCAxy{"applyPtDependentDCAxy", false, {"Global track cut"}};
    Configurable<bool> cutHasITS{"cutHasITS", true, {"Global track cut"}};
    Configurable<int> cutMinITSnCls{"cutMinITSnCls", 1, {"Global track cut"}};
    Configurable<float> cutMaxITSchi2{"cutMaxITSchi2", 36.f, {"Global track cut"}};
    Configurable<int> cutITShitsRule{"cutITShitsRule", 0, {"Global track cut"}};
    Configurable<bool> cutHasTPC{"cutHasTPC", true, {"Global track cut"}};
    Configurable<int> cutMinTPCnCls{"cutMinTPCnCls", 1, {"Global track cut"}};
    Configurable<int> cutMinTPCnClsXrows{"cutMinTPCnClsXrows", 70, {"Global track cut"}};
    Configurable<float> cutMinTPCnClsXrowsOverNcls{"cutMinTPCnClsXrowsOverNcls", 0.8f, {"Global track cut"}};
    Configurable<float> cutMaxTPCchi2{"cutMaxTPCchi2", 4.f, {"Global track cut"}};
  } cutGlobalTrack;

  struct : ConfigurableGroup {
    Configurable<bool> useThresholdsPID{"useThresholdsPID", false, {"Switch off smaller-sigma-wins pidZ."}};
    Configurable<bool> applyTauEventSelection{"applyTauEventSelection", true, {"Select tau event."}};
    Configurable<bool> cutOppositeCharge{"cutOppositeCharge", true, {"Tracks have opposite charge."}};
    Configurable<float> cutMaxAcoplanarity{"cutMaxAcoplanarity", 4 * o2::constants::math::PI / 5, {"Opening angle of the tracks. What is more goes away."}};
    Configurable<float> cutMinAcoplanarity{"cutMinAcoplanarity", 2 * o2::constants::math::PI / 5, {"Opening angle of the tracks. What is less goes away."}};
    Configurable<bool> cutElectronHasTOF{"cutElectronHasTOF", true, {"Electron is required to hit TOF."}};
    Configurable<bool> cutGoodElectron{"cutGoodElectron", true, {"Select good electron."}};
    Configurable<bool> cutOutRho{"cutOutRho", false, {"Cut out rho mass under two tracks are pions hypothesis"}};
    Configurable<bool> cutOnRho{"cutOnRho", false, {"Cut on rho mass under two tracks are pions hypothesis"}};
    Configurable<float> cutMinRhoMass{"cutMinRhoMass", 0.6, {"Lower limit on the rho mass region for cut"}};
    Configurable<float> cutMaxRhoMass{"cutMaxRhoMass", 0.95, {"Higher limit on the rho mass region for cut"}};
    Configurable<float> cutMinElectronNsigmaEl{"cutMinElectronNsigmaEl", 2.0, {"Good el hypo in. Upper n sigma cut on el hypo of selected electron. What is more goes away."}};
    Configurable<float> cutMaxElectronNsigmaEl{"cutMaxElectronNsigmaEl", -1.0, {"Good el hypo in. Lower n sigma cut on el hypo of selected electron. What is less goes away."}};
    Configurable<float> cutMinElectronNsigmaPi{"cutMinElectronNsigmaPi", -4.0, {"Good pi hypo out. Lower n sigma cut on pi hypo of selected electron. What is more till upper cut goes away."}};
    Configurable<float> cutMaxElectronNsigmaPi{"cutMaxElectronNsigmaPi", 4.0, {"Good pi hypo out. Upper n sigma cut on pi hypo of selected electron. What is less till lower cut goes away."}};
    Configurable<float> cutMinElectronNsigmaKa{"cutMinElectronNsigmaKa", -4.0, {"Good Ka hypo out. Lower n sigma cut on Ka hypo of selected electron. What is more till upper cut goes away."}};
    Configurable<float> cutMaxElectronNsigmaKa{"cutMaxElectronNsigmaKa", 4.0, {"Good Ka hypo out. Upper n sigma cut on Ka hypo of selected electron. What is less till lower cut goes away."}};
    Configurable<float> cutMinElectronNsigmaPr{"cutMinElectronNsigmaPr", -4.0, {"Good Pr hypo out. Lower n sigma cut on Pr hypo of selected electron. What is more till upper cut goes away."}};
    Configurable<float> cutMaxElectronNsigmaPr{"cutMaxElectronNsigmaPr", 4.0, {"Good Pr hypo out. Upper n sigma cut on Pr hypo of selected electron. What is less till lower cut goes away."}};
    Configurable<float> cutMinElectronTofNsigmaKa{"cutMinElectronTofNsigmaKa", -4.0, {"Good Ka TOF hypo out. Lower n sigma cut on Ka TOF hypo of selected electron. What is more till upper cut goes away."}};
    Configurable<float> cutMaxElectronTofNsigmaKa{"cutMaxElectronTofNsigmaKa", 4.0, {"Good Ka TOF hypo out. Upper n sigma cut on Ka TOF hypo of selected electron. What is less till lower cut goes away."}};
    Configurable<bool> cutPionHasTOF{"cutPionHasTOF", true, {"Pion is required to hit TOF."}};
    Configurable<bool> cutGoodMupion{"cutGoodMupion", true, {"Select good muon/pion."}};
    Configurable<float> cutMinPionNsigmaPi{"cutMinPionNsigmaPi", 4.0, {"Good pi hypo in. Upper n sigma cut on pi hypo of selected electron. What is more goes away."}};
    Configurable<float> cutMaxPionNsigmaPi{"cutMaxPionNsigmaPi", -4.0, {"Good pi hypo in. Lower n sigma cut on pi hypo of selected electron. What is less goes away."}};
    Configurable<float> cutMinPionNsigmaKa{"cutMinPionNsigmaKa", -4.0, {"Good Ka hypo out. Lower n sigma cut on Ka hypo of selected electron. What is more till upper cut goes away."}};
    Configurable<float> cutMaxPionNsigmaKa{"cutMaxPionNsigmaKa", 4.0, {"Good Ka hypo out. Upper n sigma cut on Ka hypo of selected electron. What is less till lower cut goes away."}};
    Configurable<float> cutElectronPt{"cutElectronPt", 0.9, {"Pt, where PiKaon invariant mass histos will split."}};
  } cutTauEvent;

  struct : ConfigurableGroup {
    Configurable<bool> usePIDwTOF{"usePIDwTOF", false, {"Determine whether also TOF should be used in testPIDhypothesis"}};
    Configurable<bool> useScutTOFinTPC{"useScutTOFinTPC", true, {"Determine whether cut on TOF n sigma should be used after TPC-based decision in testPIDhypothesis"}};
    Configurable<float> cutSiTPC{"cutSiTPC", 35.f, {"n sigma TPC cut on all particles in absolut values for testPIDhypothesis"}};
    Configurable<float> cutSiTOF{"cutSiTOF", 35.f, {"n sigma TOF cut on all particles in absolut values for testPIDhypothesis"}};
  } cutPID;

  struct : ConfigurableGroup {
    ConfigurableAxis zzAxisNtracks{"zzAxisNtracks", {30, -0.5, 29.5}, "Number of tracks in collision"};
    ConfigurableAxis zzAxisNparticles{"zzAxisNparticles", {100, -0.5, 99.5}, "Number of particles in collision"};
    ConfigurableAxis zzAxisZvtx{"zzAxisZvtx", {40, -20., 20.}, "Z-vertex position (cm)"};
    ConfigurableAxis zzAxisInvMass{"zzAxisInvMass", {400, 1., 5.}, "Invariant mass (GeV/c^{2})"};
    ConfigurableAxis zzAxisInvMassWide{"zzAxisInvMassWide", {1000, 0., 10.}, "Invariant mass (GeV/c^{2}), wider range"};
    ConfigurableAxis zzAxisMom{"zzAxisMom", {400, 0., 2.}, "Momentum (GeV/c)"};
    ConfigurableAxis zzAxisMomWide{"zzAxisMomWide", {1000, 0., 10.}, "Momentum (GeV/c), wider range"};
    ConfigurableAxis zzAxisMomSigned{"zzAxisMomSigned", {800, -2., 2.}, "Signed momentum (GeV/c)"};
    ConfigurableAxis zzAxisPt{"zzAxisPt", {400, 0., 2.}, "Transversal momentum (GeV/c)"};
    ConfigurableAxis zzAxisPhi{"zzAxisPhi", {64, -o2::constants::math::TwoPI, o2::constants::math::TwoPI}, "Azimuthal angle (a.y.)"};
    ConfigurableAxis zzAxisModPhi{"zzAxisModPhi", {400, 0., .4}, "Track fmod(#phi,#pi/9)"};
    ConfigurableAxis zzAxisEta{"zzAxisEta", {50, -1.2, 1.2}, "Pseudorapidity (a.u.)"};
    ConfigurableAxis zzAxisRap{"zzAxisRap", {50, -1.2, 1.2}, "Rapidity (a.u.)"};
    ConfigurableAxis zzAxisFraction{"zzAxisFraction", {500, 0., 1.}, "Fraction (-)"};
    ConfigurableAxis zzAxisMirrorFraction{"zzAxisMirrorFraction", {500, 0., 1.}, "Fraction (-)"};
    ConfigurableAxis zzAxisAcoplanarity{"zzAxisAcoplanarity", {32, 0.0, o2::constants::math::PI}, "Acoplanarity (rad)"};
    ConfigurableAxis zzAxisCollinearity{"zzAxisCollinearity", {200, 0, 20}, "Collinearity (-)"};
    ConfigurableAxis zzAxisTPCdEdx{"zzAxisTPCdEdx", {2000, 0., 200.}, "TPC dE/dx (a.u.)"};
    ConfigurableAxis zzAxisTOFsignal{"zzAxisTOFsignal", {2500, -10000., 40000.}, "TOF signal (a.u.)"};
    ConfigurableAxis zzAxisNsigma{"zzAxisNsigma", {200, -10., 10.}, "n sigma"};
    ConfigurableAxis zzAxisDCA{"zzAxisDCA", {100, -0.5, 0.5}, "DCA (cm)"};
    ConfigurableAxis zzAxisAvgITSclsSizes{"zzAxisAvgITSclsSizes", {500, 0., 10.}, "ITS average cluster size"};
    ConfigurableAxis zzAxisITSnCls{"zzAxisITSnCls", {8, -0.5, 7.5}, "ITS n clusters"};
    ConfigurableAxis zzAxisITSchi2{"zzAxisITSchi2", {100, 0, 50}, "UTS chi2"};
    ConfigurableAxis zzAxisTPCnCls{"zzAxisTPCnCls", {165, -0.5, 164.5}, "TPC n clusters"};
    ConfigurableAxis zzAxisTPCxRwsFrac{"zzAxisTPCxRwsFrac", {200, 0.0, 2.0}, "TPC fraction of crossed raws"};
    ConfigurableAxis zzAxisTPCchi2{"zzAxisTPCchi2", {100, 0, 10}, "TPC chi2"};
    ConfigurableAxis zzAxisFITtime{"zzAxisFITtime", {201, -40.5, 40.5}, "FIT time in ns"};
    ConfigurableAxis zzAxisFITamplitude{"zzAxisFITamplitude", {1000, 0., 1000.}, "FIT amplitude"};

    AxisSpec zzAxisChannels{CH_ENUM_COUNTER, -0.5, +CH_ENUM_COUNTER - 0.5, "Channels (-)"};
  } confAxis;

  using FullUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags>;
  using FullUDCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels>::iterator;
  using FullSGUDCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::SGCollisions, aod::UDZdcsReduced>::iterator;
  using FullMCUDTracks = soa::Join<aod::UDTracks, aod::UDTracksExtra, aod::UDTracksDCA, aod::UDTracksPID, aod::UDTracksFlags, aod::UDMcTrackLabels>;
  using FullMCUDCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::UDMcCollsLabels>::iterator;
  using FullMCSGUDCollision = soa::Join<aod::UDCollisions, aod::UDCollisionsSels, aod::SGCollisions, aod::UDMcCollsLabels>::iterator;

  // init
  void init(InitContext&)
  {
    if (verboseInfo)
      printMediumMessage("INIT METHOD");

    mySetITShitsRule(cutGlobalTrack.cutITShitsRule);

    if (doMainHistos) {
      histos.add("Events/hCountCollisions", ";;Number of  analysed collision (-)", HistType::kTH1D, {{1, 0.5, 1.5}});
      histos.add("Events/UDtableGapSide", ";GapSide value from UD table (-);Number of events (-)", HistType::kTH1D, {{4, -1.5, 2.5}});
      histos.add("Events/TrueGapSideDiffToTableValue", ";Difference trueGapSide from SGselector and gapSide from UD table (-);Number of events (-)", HistType::kTH1D, {{7, -3.5, 3.5}});
      histos.add("Events/hNreconstructedTracks", ";Number of tracks in a collision (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisNtracks});
      histos.add("Events/hNreconstructedPVGT", ";Number of good track particles from primary vertex in a collision (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisNtracks});
      histos.add("Events/hNreconstructedNotPVGT", ";Number of good track particles from NOT primary vertex in a collision (-);Number of events (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisNtracks});
      histos.add("Events/hNreconstructedPVGTelectrons", ";Number of good track identified electrons from primary vertex in a collision (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisNtracks});
      histos.add("Events/hNreconstructedPVGTmuons", ";Number of good track identified muons from primary vertex in a collision (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisNtracks});
      histos.add("Events/hNreconstructedPVGTpions", ";Number of good track identified pions from primary vertex in a collision (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisNtracks});
      histos.add("Events/hNreconstructedPVGTothers", ";Number of good track NOT identified electron/muon/pion particles from primary vertex in a collision (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisNtracks});
      histos.add("Events/hChannels", ";Channels (-);Number of events (-)", HistType::kTH1D, {{confAxis.zzAxisChannels}});
      histos.add("Events/FIT/hAmplitudeFT0A", ";Amplitude (-);Number of events (-)", HistType::kTH1F, {{confAxis.zzAxisFITamplitude}});
      histos.add("Events/FIT/hAmplitudeFT0C", ";Amplitude (-);Number of events (-)", HistType::kTH1F, {{confAxis.zzAxisFITamplitude}});
      histos.add("Events/FIT/hAmplitudeFDDA", ";Amplitude (-);Number of events (-)", HistType::kTH1F, {{confAxis.zzAxisFITamplitude}});
      histos.add("Events/FIT/hAmplitudeFDDC", ";Amplitude (-);Number of events (-)", HistType::kTH1F, {{confAxis.zzAxisFITamplitude}});
      histos.add("Events/FIT/hAmplitudeFV0A", ";Amplitude (-);Number of events (-)", HistType::kTH1F, {{confAxis.zzAxisFITamplitude}});
      histos.add("Events/FIT/hTimeFT0A", ";Time (ns);Number of events (-)", HistType::kTH1F, {{confAxis.zzAxisFITtime}});
      histos.add("Events/FIT/hTimeFT0C", ";Time (ns);Number of events (-)", HistType::kTH1F, {{confAxis.zzAxisFITtime}});
      histos.add("Events/FIT/hTimeFDDA", ";Time (ns);Number of events (-)", HistType::kTH1F, {{confAxis.zzAxisFITtime}});
      histos.add("Events/FIT/hTimeFDDC", ";Time (ns);Number of events (-)", HistType::kTH1F, {{confAxis.zzAxisFITtime}});
      histos.add("Events/FIT/hTimeFV0A", ";Time (ns);Number of events (-)", HistType::kTH1F, {{confAxis.zzAxisFITtime}});
      histos.add("Events/FIT/hTimeFT0AvsFT0C", ";FT0A time (ns);FT0C time (ns)", HistType::kTH2F, {{confAxis.zzAxisFITtime}, {confAxis.zzAxisFITtime}});
      histos.add("Events/FIT/hTimeFT0CvsFDDA", ";FT0C time (ns);FDDA time (ns)", HistType::kTH2F, {{confAxis.zzAxisFITtime}, {confAxis.zzAxisFITtime}});
      histos.add("Events/FIT/hTimeFDDAvsFDDC", ";FDDA time (ns);FDDC time (ns)", HistType::kTH2F, {{confAxis.zzAxisFITtime}, {confAxis.zzAxisFITtime}});
      histos.add("Events/FIT/hTimeFDDCvsFV0A", ";FDDC time (ns);FV0A time (ns)", HistType::kTH2F, {{confAxis.zzAxisFITtime}, {confAxis.zzAxisFITtime}});
      histos.add("Events/FIT/hTimeFV0AvsFT0A", ";FV0A time (ns);FT0A time (ns)", HistType::kTH2F, {{confAxis.zzAxisFITtime}, {confAxis.zzAxisFITtime}});

      histos.add("Tracks/raw/hTrackZ", ";Track z-vertex (cm);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisZvtx});
      histos.add("Tracks/raw/hTrackP", ";Track #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("Tracks/raw/hTrackPt", ";Track #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("Tracks/raw/hTrackPhi", ";Track #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("Tracks/raw/hTrackPtvsModPhi", ";Track #it{p_{T}} (GeV/c);Track fmod(#phi,#pi/9)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisModPhi});
      histos.add("Tracks/raw/hTrackPtvsModPhiTOF", ";Track #it{p_{T}} (GeV/c);Track fmod(#phi,#pi/9)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisModPhi});
      histos.add("Tracks/raw/hTrackEta", ";Track #eta (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisEta});
      histos.add("Tracks/raw/hTrackDcaXY", ";Track DCA_{XY} (cm);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisDCA});
      histos.add("Tracks/raw/hTrackPtvsDcaXY", ";Track #it{p_{T}} (GeV/c);Track DCA_{XY} (cm)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisDCA});
      histos.add("Tracks/raw/hTrackDcaZ", ";Track DCA_{Z} (cm);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisDCA});
      histos.add("Tracks/raw/ITS/itsNCls", "number of found ITS clusters;# clusters ITS", kTH1D, {confAxis.zzAxisITSnCls});
      histos.add("Tracks/raw/ITS/itsChi2NCl", "chi2 per ITS cluster;chi2 / cluster ITS", kTH1D, {confAxis.zzAxisITSchi2});
      histos.add("Tracks/raw/TPC/tpcNClsFindable", "number of findable TPC clusters;# findable clusters TPC", kTH1D, {confAxis.zzAxisTPCnCls});
      histos.add("Tracks/raw/TPC/tpcNClsFound", "number of found TPC clusters;# clusters TPC", kTH1D, {confAxis.zzAxisTPCnCls});
      histos.add("Tracks/raw/TPC/tpcCrossedRows", "number of crossed TPC rows;# crossed rows TPC", kTH1D, {confAxis.zzAxisTPCnCls});
      histos.add("Tracks/raw/TPC/tpcCrossedRowsOverFindableCls", "crossed TPC rows over findable clusters;crossed rows / findable clusters TPC", kTH1D, {confAxis.zzAxisTPCxRwsFrac});
      histos.add("Tracks/raw/TPC/tpcChi2NCl", "chi2 per cluster in TPC;chi2 / cluster TPC", kTH1D, {confAxis.zzAxisTPCchi2});
      histos.add("Tracks/GoodTrack/hTrackZ", ";Track z-vertex (cm);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisZvtx});
      histos.add("Tracks/GoodTrack/hTrackP", ";Track #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("Tracks/GoodTrack/hTrackPt", ";Track #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("Tracks/GoodTrack/hTrackPhi", ";Track #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("Tracks/GoodTrack/hTrackPtvsModPhi", ";Track #it{p_{T}} (GeV/c);Track fmod(#phi,#pi/9)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisModPhi});
      histos.add("Tracks/GoodTrack/hTrackPtvsModPhiTOF", ";Track #it{p_{T}} (GeV/c);Track fmod(#phi,#pi/9)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisModPhi});
      histos.add("Tracks/GoodTrack/hTrackEta", ";Track #eta (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisEta});
      histos.add("Tracks/GoodTrack/hTrackDcaXY", ";Track DCA_{XY} (cm);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisDCA});
      histos.add("Tracks/GoodTrack/hTrackPtvsDcaXY", ";Track #it{p_{T}} (GeV/c);Track DCA_{XY} (cm)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisDCA});
      histos.add("Tracks/GoodTrack/hTrackDcaZ", ";Track DCA_{Z} (cm);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisDCA});
      histos.add("Tracks/GoodTrack/ITS/itsNCls", "number of found ITS clusters;# clusters ITS", kTH1D, {confAxis.zzAxisITSnCls});
      histos.add("Tracks/GoodTrack/ITS/itsChi2NCl", "chi2 per ITS cluster;chi2 / cluster ITS", kTH1D, {confAxis.zzAxisITSchi2});
      histos.add("Tracks/GoodTrack/TPC/tpcNClsFindable", "number of findable TPC clusters;# findable clusters TPC", kTH1D, {confAxis.zzAxisTPCnCls});
      histos.add("Tracks/GoodTrack/TPC/tpcNClsFound", "number of found TPC clusters;# clusters TPC", kTH1D, {confAxis.zzAxisTPCnCls});
      histos.add("Tracks/GoodTrack/TPC/tpcCrossedRows", "number of crossed TPC rows;# crossed rows TPC", kTH1D, {confAxis.zzAxisTPCnCls});
      histos.add("Tracks/GoodTrack/TPC/tpcCrossedRowsOverFindableCls", "crossed TPC rows over findable clusters;crossed rows / findable clusters TPC", kTH1D, {confAxis.zzAxisTPCxRwsFrac});
      histos.add("Tracks/GoodTrack/TPC/tpcChi2NCl", "chi2 per cluster in TPC;chi2 / cluster TPC", kTH1D, {confAxis.zzAxisTPCchi2});
    }

    if (doPIDhistos) {
      histos.add("Tracks/raw/PID/hTPCsignalVsZ", "All tracks;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisZvtx, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/raw/PID/hTPCsignalVsP", "All tracks;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/raw/PID/hTPCsignalVsPt", "All tracks;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/raw/PID/hTPCsignalVsEta", "All tracks;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisEta, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/raw/PID/hTPCsignalVsPhi", "All tracks;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/raw/PID/hTOFsignalVsP", "All tracks;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("Tracks/raw/PID/hTPCnSigmaElVsP", ";Track #it{p} (GeV/c);n#sigma^{e}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("Tracks/raw/PID/hTPCnSigmaMuVsP", ";Track #it{p} (GeV/c);n#sigma^{#mu}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("Tracks/raw/PID/hTPCnSigmaPiVsP", ";Track #it{p} (GeV/c);n#sigma^{#pi}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("Tracks/raw/PID/hTPCnSigmaKaVsP", ";Track #it{p} (GeV/c);n#sigma^{K}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("Tracks/raw/PID/hTPCnSigmaPrVsP", ";Track #it{p} (GeV/c);n#sigma^{p}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("Tracks/raw/PID/PosCharge/hTPCsignalVsZ", "Positively charged tracks;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisZvtx, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/raw/PID/PosCharge/hTPCsignalVsP", "Positively charged tracks;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/raw/PID/PosCharge/hTPCsignalVsPt", "Positively charged tracks;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/raw/PID/PosCharge/hTPCsignalVsEta", "Positively charged tracks;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisEta, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/raw/PID/PosCharge/hTPCsignalVsPhi", "Positively charged tracks;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/raw/PID/PosCharge/hTOFsignalVsP", "Positively charged tracks;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("Tracks/raw/PID/NegCharge/hTPCsignalVsZ", "Negatively charged tracks;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisZvtx, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/raw/PID/NegCharge/hTPCsignalVsP", "Negatively charged tracks;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/raw/PID/NegCharge/hTPCsignalVsPt", "Negatively charged tracks;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/raw/PID/NegCharge/hTPCsignalVsEta", "Negatively charged tracks;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisEta, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/raw/PID/NegCharge/hTPCsignalVsPhi", "Negatively charged tracks;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/raw/PID/NegCharge/hTOFsignalVsP", "Negatively charged tracks;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("Tracks/GoodTrack/PID/hTPCsignalVsZ", "All good tracks;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisZvtx, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/hTPCsignalVsP", "All good tracks;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/hTPCsignalVsPt", "All good tracks;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/hTPCsignalVsEta", "All good tracks;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisEta, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/hTPCsignalVsPhi", "All good tracks;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/hTPCnSigmaElVsP", ";Track #it{p} (GeV/c);n#sigma^{e}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/hTPCnSigmaMuVsP", ";Track #it{p} (GeV/c);n#sigma^{#mu}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/hTPCnSigmaPiVsP", ";Track #it{p} (GeV/c);n#sigma^{#pi}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/hTPCnSigmaKaVsP", ";Track #it{p} (GeV/c);n#sigma^{K}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/hTPCnSigmaPrVsP", ";Track #it{p} (GeV/c);n#sigma^{p}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/hTOFsignalVsP", "All good tracks;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsZ", "Positively charged good tracks;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisZvtx, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsP", "Positively charged good tracks;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsPt", "Positively charged good tracks;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsEta", "Positively charged good tracks;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisEta, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsPhi", "Positively charged good tracks;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/PosCharge/hTOFsignalVsP", "Positively charged good tracks;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsZ", "Negatively charged good tracks;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisZvtx, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsP", "Negatively charged good tracks;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsPt", "Negatively charged good tracks;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsEta", "Negatively charged good tracks;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisEta, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsPhi", "Negatively charged good tracks;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/NegCharge/hTOFsignalVsP", "Negatively charged good tracks;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("Tracks/GoodTrack/PID/Electron/hTPCsignalVsZ", "Identified electrons;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisZvtx, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Electron/hTPCsignalVsP", "Identified electrons;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Electron/hTPCsignalVsPt", "Identified electrons;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Electron/hTPCsignalVsEta", "Identified electrons;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisEta, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Electron/hTPCsignalVsPhi", "Identified electrons;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Electron/hTOFsignalVsP", "Identified electrons;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("Tracks/GoodTrack/PID/Electron/hTPCnSigmaVsP", "Identified electrons;Track #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Electron/hTOFnSigmaVsP", "Identified electrons;Track #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsMu", "Identified electrons;n#sigma^{#it{e}}_{TPC} (arb. units);n#sigma^{#mu}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsMu", "Identified electrons;n#sigma^{#it{e}}_{TOF} (arb. units);n#sigma^{#mu}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsPi", "Identified electrons;n#sigma^{#it{e}}_{TPC} (arb. units);n#sigma^{#pi}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsPi", "Identified electrons;n#sigma^{#it{e}}_{TOF} (arb. units);n#sigma^{#pi}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsKa", "Identified electrons;n#sigma^{#it{e}}_{TPC} (arb. units);n#sigma^{#it{K}}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsKa", "Identified electrons;n#sigma^{#it{e}}_{TOF} (arb. units);n#sigma^{#it{K}}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsPr", "Identified electrons;n#sigma^{#it{e}}_{TPC} (arb. units);n#sigma^{p}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsPr", "Identified electrons;n#sigma^{#it{e}}_{TOF} (arb. units);n#sigma^{p}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Muon/hTPCsignalVsZ", "Identified muons;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisZvtx, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Muon/hTPCsignalVsP", "Identified muons;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Muon/hTPCsignalVsPt", "Identified muons;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Muon/hTPCsignalVsEta", "Identified muons;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisEta, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Muon/hTPCsignalVsPhi", "Identified muons;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Muon/hTOFsignalVsP", "Identified muons;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("Tracks/GoodTrack/PID/Muon/hTPCnSigmaVsP", "Identified muons;Track #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Muon/hTOFnSigmaVsP", "Identified muons;Track #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Muon/hTPCnSigmaMuVsEl", "Identified muons;n#sigma^{#mu}_{TPC} (arb. units);n#sigma^{#it{e}}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Muon/hTOFnSigmaMuVsEl", "Identified muons;n#sigma^{#mu}_{TOF} (arb. units);n#sigma^{#it{e}}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Muon/hTPCnSigmaMuVsPi", "Identified muons;n#sigma^{#mu}_{TPC} (arb. units);n#sigma^{#pi}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Muon/hTOFnSigmaMuVsPi", "Identified muons;n#sigma^{#mu}_{TOF} (arb. units);n#sigma^{#pi}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Muon/hTPCnSigmaMuVsKa", "Identified muons;n#sigma^{#mu}_{TPC} (arb. units);n#sigma^{#it{K}}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Muon/hTOFnSigmaMuVsKa", "Identified muons;n#sigma^{#mu}_{TOF} (arb. units);n#sigma^{#it{K}}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Muon/hTPCnSigmaMuVsPr", "Identified muons;n#sigma^{#mu}_{TPC} (arb. units);n#sigma^{p}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Muon/hTOFnSigmaMuVsPr", "Identified muons;n#sigma^{#mu}_{TOF} (arb. units);n#sigma^{p}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Pion/hTPCsignalVsZ", "Identified pions;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisZvtx, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Pion/hTPCsignalVsP", "Identified pions;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Pion/hTPCsignalVsPt", "Identified pions;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Pion/hTPCsignalVsEta", "Identified pions;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisEta, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Pion/hTPCsignalVsPhi", "Identified pions;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Pion/hTOFsignalVsP", "Identified pions;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("Tracks/GoodTrack/PID/Pion/hTPCnSigmaVsP", "Identified pions;Track #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Pion/hTOFnSigmaVsP", "Identified pions;Track #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsEl", "Identified pions;n#sigma^{#pi}_{TPC} (arb. units);n#sigma^{#it{e}}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsEl", "Identified pions;n#sigma^{#pi}_{TOF} (arb. units);n#sigma^{#it{e}}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsMu", "Identified pions;n#sigma^{#pi}_{TPC} (arb. units);n#sigma^{#mu}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsMu", "Identified pions;n#sigma^{#pi}_{TOF} (arb. units);n#sigma^{#mu}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsKa", "Identified pions;n#sigma^{#pi}_{TPC} (arb. units);n#sigma^{#it{K}}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsKa", "Identified pions;n#sigma^{#pi}_{TOF} (arb. units);n#sigma^{#it{K}}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsPr", "Identified pions;n#sigma^{#pi}_{TPC} (arb. units);n#sigma^{p}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsPr", "Identified pions;n#sigma^{#pi}_{TOF} (arb. units);n#sigma^{p}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("Tracks/GoodTrack/PID/Others/hTPCsignalVsZ", "Identified NOT electron/Muon/Pion;Track z-vertex (cm);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisZvtx, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Others/hTPCsignalVsP", "Identified NOT electron/Muon/Pion;Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Others/hTPCsignalVsPt", "Identified NOT electron/Muon/Pion;Track #it{p_{#rm T}} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Others/hTPCsignalVsEta", "Identified NOT electron/Muon/Pion;Track #eta (-);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisEta, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Others/hTPCsignalVsPhi", "Identified NOT electron/Muon/Pion;Track #phi (rad);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisTPCdEdx});
      histos.add("Tracks/GoodTrack/PID/Others/hTOFsignalVsP", "Identified NOT electron/Muon/Pion;Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
    }

    if (doTwoTracks) {
      histos.add("EventTwoTracks/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventTwoTracks/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/hInvariantMassWideNoMothers", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/hInvariantMassWideAllPionMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/hInvariantMassWideAllPionMassPtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/hInvariantMassWideAllPionMassTOF", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisAcoplanarity});
      histos.add("EventTwoTracks/hCollinearity", ";#DeltaR (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisCollinearity});
      histos.add("EventTwoTracks/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("EventTwoTracks/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventTwoTracks/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisRap});
      histos.add("EventTwoTracks/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/hDaughtersEnergyAsymmetry", ";(E_{electron} - E_{#mu/#pi}) / (E_{electron} + E_{#mu/#pi});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMirrorFraction});
      histos.add("EventTwoTracks/hDaughtersMomentaAsymmetry", ";(#it{p}_{electron} - #it{p}_{#mu/#pi}) / (#it{p}_{electron} + #it{p}_{#mu/#pi});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMirrorFraction});
      histos.add("EventTwoTracks/hDaughtersPtAsymmetry", ";(#it{p_{T}}_{electron} - #it{p_{T}}_{#mu/#pi}) / (#it{p_{T}}_{electron} + #it{p_{T}}_{#mu/#pi});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMirrorFraction});
      histos.add("EventTwoTracks/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisMom});
      histos.add("EventTwoTracks/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMomWide, confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {confAxis.zzAxisRap, confAxis.zzAxisRap});
      histos.add("EventTwoTracks/hDaughtersEnergyFractions", ";E_{daughter 1} / E_{tot} (-);E_{daughter 1} / E_{tot} (-)", HistType::kTH2D, {confAxis.zzAxisFraction, confAxis.zzAxisFraction});
      histos.add("EventTwoTracks/hDaughtersPvsITSclusterSize", ";Average ITS cluster size;Daughter #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisAvgITSclsSizes, confAxis.zzAxisMomSigned});
      histos.add("EventTwoTracks/hDaughtersPvsITSclusterSizeXcos", ";Average ITS cluster size x cos(#lambda);Daughter #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisAvgITSclsSizes, confAxis.zzAxisMomSigned});
      histos.add("EventTwoTracks/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("EventTwoTracks/PID/hTOFsignalVsP", ";Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("EventTwoTracks/PID/hTPCnSigmaElVsP", ";Track #it{p} (GeV/c);n#sigma^{e}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/PID/hTPCnSigmaMuVsP", ";Track #it{p} (GeV/c);n#sigma^{#mu}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/PID/hTPCnSigmaPiVsP", ";Track #it{p} (GeV/c);n#sigma^{#pi}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/PID/hTPCnSigmaKaVsP", ";Track #it{p} (GeV/c);n#sigma^{K}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/PID/hTPCnSigmaPrVsP", ";Track #it{p} (GeV/c);n#sigma^{p}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/PID/NoPID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("EventTwoTracks/PID/NoPID/hTPCnSigmaElVsP", ";Track #it{p} (GeV/c);n#sigma^{e}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/PID/NoPID/hTPCnSigmaMuVsP", ";Track #it{p} (GeV/c);n#sigma^{#mu}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/PID/NoPID/hTPCnSigmaPiVsP", ";Track #it{p} (GeV/c);n#sigma^{#pi}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/PID/NoPID/hTPCnSigmaKaVsP", ";Track #it{p} (GeV/c);n#sigma^{K}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/PID/NoPID/hTPCnSigmaPrVsP", ";Track #it{p} (GeV/c);n#sigma^{p}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});

      histos.add("EventTwoTracks/TwoElectrons/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventTwoTracks/TwoElectrons/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/TwoElectrons/hInvariantMassWidePtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/TwoElectrons/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisAcoplanarity});
      histos.add("EventTwoTracks/TwoElectrons/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("EventTwoTracks/TwoElectrons/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/TwoElectrons/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventTwoTracks/TwoElectrons/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/TwoElectrons/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisRap});
      histos.add("EventTwoTracks/TwoElectrons/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisMom});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMomWide, confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhi", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisModPhi});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiTOF", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisModPhi});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiPtCut", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisModPhi});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiPtCutTOF", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisModPhi});
      histos.add("EventTwoTracks/TwoElectrons/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {confAxis.zzAxisRap, confAxis.zzAxisRap});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingP", ";Leading #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingPwide", ";Leading #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingPt", ";Leading #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingPhi", ";Leading #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingRapidity", ";Leading #it{y} (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisRap});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingPvsOtherP", ";Leading #it{p} (GeV/c); Other #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisMom});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingPwideVsOtherPwide", ";Leading #it{p} (GeV/c); Other #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMomWide, confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingPtVsOtherPt", ";Leading #it{p_{T}} (GeV/c); Other #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingPhiVsOtherPhi", ";Leading #phi (rad); Other #phi (rad)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/TwoElectrons/hLeadingRapVsOtherRap", ";Leading #it{y} (-); Other #it{y} (-)", HistType::kTH2D, {confAxis.zzAxisRap, confAxis.zzAxisRap});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsLP", ";Leading #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsOP", ";Other #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsP", ";Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsLP", ";Leading #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsOP", ";Other #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsP", ";Track #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsLP", ";Leading #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsOP", ";Other #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsP", ";Track #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsLP", ";Leading #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsOP", ";Other #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});

      histos.add("EventTwoTracks/TwoMuons/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventTwoTracks/TwoMuons/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/TwoMuons/hInvariantMassWidePtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/TwoMuons/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisAcoplanarity});
      histos.add("EventTwoTracks/TwoMuons/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("EventTwoTracks/TwoMuons/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/TwoMuons/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventTwoTracks/TwoMuons/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/TwoMuons/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisRap});
      histos.add("EventTwoTracks/TwoMuons/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisMom});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMomWide, confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhi", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisModPhi});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiTOF", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisModPhi});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiPtCut", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisModPhi});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiPtCutTOF", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisModPhi});
      histos.add("EventTwoTracks/TwoMuons/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {confAxis.zzAxisRap, confAxis.zzAxisRap});
      histos.add("EventTwoTracks/TwoMuons/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});

      histos.add("EventTwoTracks/TwoPions/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventTwoTracks/TwoPions/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/TwoPions/hInvariantMassWidePtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/TwoPions/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisAcoplanarity});
      histos.add("EventTwoTracks/TwoPions/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("EventTwoTracks/TwoPions/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/TwoPions/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventTwoTracks/TwoPions/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/TwoPions/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisRap});
      histos.add("EventTwoTracks/TwoPions/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisMom});
      histos.add("EventTwoTracks/TwoPions/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMomWide, confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/TwoPions/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/TwoPions/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/TwoPions/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/TwoPions/hDaughtersPtvsModPhi", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisModPhi});
      histos.add("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiTOF", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisModPhi});
      histos.add("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiPtCut", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisModPhi});
      histos.add("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiPtCutTOF", ";Daughter #it{p_{T}} (GeV/c);Daughter fmod(#phi,#pi/9)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisModPhi});
      histos.add("EventTwoTracks/TwoPions/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {confAxis.zzAxisRap, confAxis.zzAxisRap});
      histos.add("EventTwoTracks/TwoPions/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});

      histos.add("EventTwoTracks/ElectronMuon/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventTwoTracks/ElectronMuon/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/ElectronMuon/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisAcoplanarity});
      histos.add("EventTwoTracks/ElectronMuon/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("EventTwoTracks/ElectronMuon/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/ElectronMuon/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronMuon/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/ElectronMuon/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisRap});
      histos.add("EventTwoTracks/ElectronMuon/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronMuon/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisMom});
      histos.add("EventTwoTracks/ElectronMuon/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMomWide, confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/ElectronMuon/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronMuon/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/ElectronMuon/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {confAxis.zzAxisRap, confAxis.zzAxisRap});
      histos.add("EventTwoTracks/ElectronMuon/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});

      histos.add("EventTwoTracks/ElectronPion/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventTwoTracks/ElectronPion/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/ElectronPion/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisAcoplanarity});
      histos.add("EventTwoTracks/ElectronPion/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("EventTwoTracks/ElectronPion/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/ElectronPion/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronPion/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/ElectronPion/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisRap});
      histos.add("EventTwoTracks/ElectronPion/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronPion/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisMom});
      histos.add("EventTwoTracks/ElectronPion/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMomWide, confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/ElectronPion/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronPion/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/ElectronPion/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {confAxis.zzAxisRap, confAxis.zzAxisRap});
      histos.add("EventTwoTracks/ElectronPion/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});

      histos.add("EventTwoTracks/MuonPion/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventTwoTracks/MuonPion/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/MuonPion/hInvariantMassWidePtCut", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/MuonPion/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisAcoplanarity});
      histos.add("EventTwoTracks/MuonPion/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("EventTwoTracks/MuonPion/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/MuonPion/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventTwoTracks/MuonPion/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/MuonPion/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisRap});
      histos.add("EventTwoTracks/MuonPion/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/MuonPion/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisMom});
      histos.add("EventTwoTracks/MuonPion/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMomWide, confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/MuonPion/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/MuonPion/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/MuonPion/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {confAxis.zzAxisRap, confAxis.zzAxisRap});
      histos.add("EventTwoTracks/MuonPion/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});

      histos.add("EventTwoTracks/ElectronMuPi/hNeventsPtCuts", ";Selection (-);Number of events (-)", HistType::kTH1D, {{20, -0.5, 19.5}});
      histos.add("EventTwoTracks/ElectronMuPi/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventTwoTracks/ElectronMuPi/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/ElectronMuPi/PionsSelection/hInvariantMass", ";Invariant mass (#pi^{+}#pi^{-}) (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventTwoTracks/ElectronMuPi/PionsSelection/hInvariantMassWide", ";Invariant mass (#pi^{+}#pi^{-}) (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/ElectronMuPi/KstarSelection/hInvariantMass", ";Invariant mass (K#pi) (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventTwoTracks/ElectronMuPi/KstarSelection/hInvariantMassWide", ";Invariant mass (K#pi) (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/ElectronMuPi/KstarSelection/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronMuPi/KstarSelection/hMotherMassVsElectronP", ";Invariant mass (K#pi) (GeV/c^{2});Electron #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/ElectronMuPi/KstarSelection/hIMKvsIMe", ";#pi+K invariant mass (GeV/c^{2});#pi+e invariant mass (GeV/c^{2})", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutLow/hInvariantMass", ";Invariant mass (K#pi) (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutLow/hInvariantMassWide", ";Invariant mass (K#pi) (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutLow/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutLow/hIMKvsIMe", ";#pi+K invariant mass (GeV/c^{2});#pi+e invariant mass (GeV/c^{2})", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutLow/hElectronPwideVsOtherPwide", ";Electron #it{p} (GeV/c); #mu/#pi #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMomWide, confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutHigh/hInvariantMass", ";Invariant mass (K#pi) (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutHigh/hInvariantMassWide", ";Invariant mass (K#pi) (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutHigh/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutHigh/hIMKvsIMe", ";#pi+K invariant mass (GeV/c^{2});#pi+e invariant mass (GeV/c^{2})", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutHigh/hElectronPwideVsOtherPwide", ";Electron #it{p} (GeV/c); #mu/#pi #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMomWide, confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/ElectronMuPi/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisAcoplanarity});
      histos.add("EventTwoTracks/ElectronMuPi/hCollinearity", ";#DeltaR (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisCollinearity});
      histos.add("EventTwoTracks/ElectronMuPi/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("EventTwoTracks/ElectronMuPi/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/ElectronMuPi/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronMuPi/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/ElectronMuPi/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisRap});
      histos.add("EventTwoTracks/ElectronMuPi/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronMuPi/hMotherMassVsElectronP", ";Invariant mass (GeV/c^{2});Electron #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/ElectronMuPi/hMotherMassVsAcoplanarity", ";Invariant mass (GeV/c^{2});#Delta#phi (rad)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisAcoplanarity});
      histos.add("EventTwoTracks/ElectronMuPi/hElectronPt", ";Electron #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronMuPi/hElectronPtWide", ";Electron #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/ElectronMuPi/hDaughtersEnergyAsymmetry", ";(E_{electron} - E_{#mu/#pi}) / (E_{electron} + E_{#mu/#pi});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMirrorFraction});
      histos.add("EventTwoTracks/ElectronMuPi/hDaughtersMomentaAsymmetry", ";(#it{p}_{electron} - #it{p}_{#mu/#pi}) / (#it{p}_{electron} + #it{p}_{#mu/#pi});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMirrorFraction});
      histos.add("EventTwoTracks/ElectronMuPi/hDaughtersPtAsymmetry", ";(#it{p_{T}}_{electron} - #it{p_{T}}_{#mu/#pi}) / (#it{p_{T}}_{electron} + #it{p_{T}}_{#mu/#pi});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMirrorFraction});
      histos.add("EventTwoTracks/ElectronMuPi/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisMom});
      histos.add("EventTwoTracks/ElectronMuPi/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMomWide, confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/ElectronMuPi/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronMuPi/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/ElectronMuPi/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {confAxis.zzAxisRap, confAxis.zzAxisRap});
      histos.add("EventTwoTracks/ElectronMuPi/hDaughtersEnergyFractions", ";E_{electron} / E_{tot} (-);E_{#mu/#pi} / E_{tot} (-)", HistType::kTH2D, {confAxis.zzAxisFraction, confAxis.zzAxisFraction});
      histos.add("EventTwoTracks/ElectronMuPi/hElectronPvsOtherP", ";Electron #it{p} (GeV/c); #mu/#pi #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisMom});
      histos.add("EventTwoTracks/ElectronMuPi/hElectronPwideVsOtherPwide", ";Electron #it{p} (GeV/c); #mu/#pi #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMomWide, confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/ElectronMuPi/hElectronPvsAcoplanarity", ";Electron #it{p} (GeV/c); #Delta#phi (rad)", HistType::kTH2D, {confAxis.zzAxisMomWide, confAxis.zzAxisAcoplanarity});
      histos.add("EventTwoTracks/ElectronMuPi/hOtherPvsAcoplanarity", ";#mu/#pi #it{p} (GeV/c); #Delta#phi (rad)", HistType::kTH2D, {confAxis.zzAxisMomWide, confAxis.zzAxisAcoplanarity});
      histos.add("EventTwoTracks/ElectronMuPi/hElectronPtVsOtherPt", ";Electron #it{p_{T}} (GeV/c); #mu/#pi #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronMuPi/hElectronPhiVsOtherPhi", ";Electron #phi (rad); #mu/#pi #phi (rad)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/ElectronMuPi/hElectronRapVsOtherRap", ";Electron #it{y} (-); #mu/#pi #it{y} (-)", HistType::kTH2D, {confAxis.zzAxisRap, confAxis.zzAxisRap});

      histos.add("EventTwoTracks/ElectronMuPi/PID/mcTruth/nSigmaTPC1", "Paul's way;True electron #it{p} (GeV/c);n#sigma^{e}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/mcTruth/nSigmaTPC2", "Paul's way;True not-electron #it{p} (GeV/c);n#sigma^{e}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});

      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCsignalVsEPofE", ";Electron #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCsignalVsOPofO", ";#mu/#pi #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsP", ";Track #it{p} (GeV/c);n#sigma^{e}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsEPofE", ";Electron #it{p} (GeV/c);n#sigma^{e}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsPPofE", ";Electron #it{p} (GeV/c);n#sigma^{#pi}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaEvsnSigmaPofE", ";Electron n#sigma^{e}_{TPC} (arb. units);Electron n#sigma^{#pi}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsEPofO", ";Non-electron #it{p} (GeV/c);n#sigma^{e}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsMPofO", ";Non-electron #it{p} (GeV/c);n#sigma^{#mu}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsPPofO", ";Non-electron #it{p} (GeV/c);n#sigma^{#pi}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaEvsnSigmaPofO", ";Non-electron n#sigma^{e}_{TPC} (arb. units);Non-electron n#sigma^{#pi}_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsP", ";Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsP", ";Track #it{p} (GeV/c);n#sigma^{e}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsEPofE", ";Electron #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsEPofE", ";Electron #it{p} (GeV/c);n#sigma^{e}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsPPofE", ";Electron #it{p} (GeV/c);n#sigma^{#pi}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaEvsnSigmaPofE", ";Electron n#sigma^{e}_{TOF} (arb. units);Electron n#sigma^{#pi}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsOPofO", ";Not-electron #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsEPofO", ";Not-electron #it{p} (GeV/c);n#sigma^{e}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsMPofO", ";Not-electron #it{p} (GeV/c);n#sigma^{#mu}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsPPofO", ";Not-electron #it{p} (GeV/c);n#sigma^{#pi}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaEvsnSigmaPofO", ";Not-electron n#sigma^{e}_{TOF} (arb. units);Not-electron n#sigma^{#pi}_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisNsigma, confAxis.zzAxisNsigma});

      histos.add("EventTwoTracks/ElectronOther/hNeventsPtCuts", ";Selection (-);Number of events (-)", HistType::kTH1D, {{20, -0.5, 19.5}});
      histos.add("EventTwoTracks/ElectronOther/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventTwoTracks/ElectronOther/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventTwoTracks/ElectronOther/hAcoplanarity", ";#Delta#phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisAcoplanarity});
      histos.add("EventTwoTracks/ElectronOther/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("EventTwoTracks/ElectronOther/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/ElectronOther/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronOther/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/ElectronOther/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisRap});
      histos.add("EventTwoTracks/ElectronOther/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronOther/hElectronPt", ";Electron #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronOther/hElectronPtWide", ";Electron #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/ElectronOther/hDaughtersP", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisMom});
      histos.add("EventTwoTracks/ElectronOther/hDaughtersPwide", ";Daughter 1 #it{p} (GeV/c);Daughter 2 #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMomWide, confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/ElectronOther/hDaughtersPt", ";Daughter 1 #it{p_{T}} (GeV/c);Daughter 2 #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronOther/hDaughtersPhi", ";Daughter 1 #phi (rad);Daughter 2 #phi (rad)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/ElectronOther/hDaughtersRapidity", ";Daughter 1 #it{y} (-);Daughter 2 #it{y} (-)", HistType::kTH2D, {confAxis.zzAxisRap, confAxis.zzAxisRap});
      histos.add("EventTwoTracks/ElectronOther/hElectronPvsOtherP", ";Electron #it{p} (GeV/c); Other #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisMom});
      histos.add("EventTwoTracks/ElectronOther/hElectronPwideVsOtherPwide", ";Electron #it{p} (GeV/c); Other #it{p} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisMomWide, confAxis.zzAxisMomWide});
      histos.add("EventTwoTracks/ElectronOther/hElectronPtVsOtherPt", ";Electron #it{p_{T}} (GeV/c); Other #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisPt, confAxis.zzAxisPt});
      histos.add("EventTwoTracks/ElectronOther/hElectronPhiVsOtherPhi", ";Electron #phi (rad); Other #phi (rad)", HistType::kTH2D, {confAxis.zzAxisPhi, confAxis.zzAxisPhi});
      histos.add("EventTwoTracks/ElectronOther/hElectronRapVsOtherRap", ";Electron #it{y} (-); Other #it{y} (-)", HistType::kTH2D, {confAxis.zzAxisRap, confAxis.zzAxisRap});
      histos.add("EventTwoTracks/ElectronOther/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("EventTwoTracks/ElectronOther/PID/hTPCsignalVsEP", ";Electron #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("EventTwoTracks/ElectronOther/PID/hTPCsignalVsOP", ";#it{e}/#mu/#pi #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
      histos.add("EventTwoTracks/ElectronOther/PID/hTOFsignalVsP", ";Track #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("EventTwoTracks/ElectronOther/PID/hTOFsignalVsEP", ";Electron #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("EventTwoTracks/ElectronOther/PID/hTOFsignalVsOP", ";Other #it{p} (GeV/c);TOF signal (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTOFsignal});
      histos.add("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsP", ";Track #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsEP", ";Electron #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsMP", ";Muon #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsPP", ";Pion #it{p} (GeV/c);n#sigma_{TPC} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsP", ";Track #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsEP", ";Electron #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsMP", ";Muon #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
      histos.add("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsPP", ";Pion #it{p} (GeV/c);n#sigma_{TOF} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisNsigma});
    }

    if (doFourTracks) {
      histos.add("EventFourTracks/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventFourTracks/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventFourTracks/hInvariantMassWideNoMothers", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventFourTracks/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("EventFourTracks/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventFourTracks/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventFourTracks/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("EventFourTracks/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisRap});
      histos.add("EventFourTracks/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventFourTracks/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});

      histos.add("EventFourTracks/WithElectron/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventFourTracks/WithElectron/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventFourTracks/WithElectron/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("EventFourTracks/WithElectron/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventFourTracks/WithElectron/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventFourTracks/WithElectron/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("EventFourTracks/WithElectron/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisRap});
      histos.add("EventFourTracks/WithElectron/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventFourTracks/WithElectron/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});

      histos.add("EventFourTracks/WithMuon/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventFourTracks/WithMuon/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventFourTracks/WithMuon/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("EventFourTracks/WithMuon/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventFourTracks/WithMuon/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventFourTracks/WithMuon/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("EventFourTracks/WithMuon/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisRap});
      histos.add("EventFourTracks/WithMuon/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventFourTracks/WithMuon/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});

      histos.add("EventFourTracks/WithPion/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventFourTracks/WithPion/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventFourTracks/WithPion/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("EventFourTracks/WithPion/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventFourTracks/WithPion/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventFourTracks/WithPion/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("EventFourTracks/WithPion/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisRap});
      histos.add("EventFourTracks/WithPion/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventFourTracks/WithPion/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
    }

    if (doSixTracks) {
      histos.add("EventSixTracks/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventSixTracks/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventSixTracks/hInvariantMassWideNoMothers", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventSixTracks/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("EventSixTracks/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventSixTracks/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventSixTracks/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("EventSixTracks/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisRap});
      histos.add("EventSixTracks/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventSixTracks/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});

      histos.add("EventSixTracks/SixPions/hInvariantMass", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMass});
      histos.add("EventSixTracks/SixPions/hInvariantMassWide", ";Invariant mass (GeV/c^{2});Number of events (-)", HistType::kTH1D, {confAxis.zzAxisInvMassWide});
      histos.add("EventSixTracks/SixPions/hMotherP", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("EventSixTracks/SixPions/hMotherPwide", ";Mother #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMomWide});
      histos.add("EventSixTracks/SixPions/hMotherPt", ";Mother #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("EventSixTracks/SixPions/hMotherPhi", ";Mother #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("EventSixTracks/SixPions/hMotherRapidity", ";Mother #it{y} (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisRap});
      histos.add("EventSixTracks/SixPions/hMotherMassVsPt", ";Invariant mass (GeV/c^{2});Mother #it{p_{T}} (GeV/c)", HistType::kTH2D, {confAxis.zzAxisInvMassWide, confAxis.zzAxisPt});
      histos.add("EventSixTracks/SixPions/PID/hTPCsignalVsP", ";Track #it{p} (GeV/c);TPC d#it{E}/d#it{x} (arb. units)", HistType::kTH2D, {confAxis.zzAxisMom, confAxis.zzAxisTPCdEdx});
    }

    if (doTruthHistos) {
      histos.add("Events/Truth/hCountCollisions", ";;Number of generated collision (-)", HistType::kTH1D, {{1, 0.5, 1.5}});
      histos.add("Events/Truth/hChannels", ";Channels (-);Number of events (-)", HistType::kTH1D, {{confAxis.zzAxisChannels}});
      histos.add("Events/Truth/hPDGcodesAll", ";PDG codes of all particles (-);Number of events (-)", HistType::kTH1D, {{2001, -1000, 1000}});
      histos.add("Events/Truth/hPDGcodesNoMother", ";PDG codes of particles without mother (-);Number of events (-)", HistType::kTH1D, {{2001, -1000, 1000}});
      histos.add("Events/Truth/hPDGcodesTauDaughters", ";PDG codes of daughters of particles without mother (-);Number of events (-)", HistType::kTH1D, {{2001, -1000, 1000}});
      histos.add("Events/Truth/hNparticles", ";Number of particles in a collision (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisNparticles});
      histos.add("Events/Truth/hNtauDaughters", ";Number of daughters of no-mother particle in a collision (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisNparticles});
      histos.add("Events/Truth/hNelectrons", ";Number of electrons in a collision (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisNparticles});
      histos.add("Events/Truth/hNmuons", ";Number of muons in a collision (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisNparticles});
      histos.add("Events/Truth/hNpions", ";Number of pions in a collision (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisNparticles});
      histos.add("Events/Truth/hNphysPartVsNwoutMotherParts", ";Number of physical primary particles (-);Number of particles without mother(-)", HistType::kTH2D, {confAxis.zzAxisNparticles, confAxis.zzAxisNparticles});
      histos.add("Tracks/Truth/hTauP", ";Tau #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("Tracks/Truth/hTauPt", ";Tau #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("Tracks/Truth/hTauPhi", ";Tau #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("Tracks/Truth/hTauEta", ";Tau #eta (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisEta});
      histos.add("Tracks/Truth/hElectronP", ";Electron #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("Tracks/Truth/hElectronPt", ";Electron #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("Tracks/Truth/hElectronPhi", ";Electron #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("Tracks/Truth/hElectronEta", ";Electron #eta (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisEta});
      histos.add("Tracks/Truth/hMuonP", ";Muon #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("Tracks/Truth/hMuonPt", ";Muon #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("Tracks/Truth/hMuonPhi", ";Muon #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("Tracks/Truth/hMuonEta", ";Muon #eta (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisEta});
      histos.add("Tracks/Truth/hPionP", ";Pion #it{p} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisMom});
      histos.add("Tracks/Truth/hPionPt", ";Pion #it{p_{T}} (GeV/c);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPt});
      histos.add("Tracks/Truth/hPionPhi", ";Pion #phi (rad);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisPhi});
      histos.add("Tracks/Truth/hPionEta", ";Pion #eta (-);Number of events (-)", HistType::kTH1D, {confAxis.zzAxisEta});
    }

  } // end init

  // run (always called before process :( )
  void run(ProcessingContext&)
  {

    if (verboseInfo)
      printLargeMessage("RUN METHOD");

  } // end run

  std::vector<std::pair<int8_t, std::set<uint8_t>>> cutMyRequiredITSHits{};

  void mySetRequireHitsInITSLayers(int8_t minNRequiredHits, std::set<uint8_t> requiredLayers)
  {
    // layer 0 corresponds to the the innermost ITS layer
    cutMyRequiredITSHits.push_back(std::make_pair(minNRequiredHits, requiredLayers));
  }

  void mySetITShitsRule(int matching)
  {
    switch (matching) {
      case 0: // Run3ITSibAny
        mySetRequireHitsInITSLayers(1, {0, 1, 2});
        break;
      case 1: // Run3ITSibTwo
        mySetRequireHitsInITSLayers(2, {0, 1, 2});
        break;
      case 2: // Run3ITSallAny
        mySetRequireHitsInITSLayers(1, {0, 1, 2, 3, 4, 5, 6});
        break;
      case 3: // Run3ITSall7Layers
        mySetRequireHitsInITSLayers(7, {0, 1, 2, 3, 4, 5, 6});
        break;
      default:
        LOG(fatal) << "You chose wrong ITS matching";
        break;
    }
  }

  bool isFulfillsITSHitRequirementsReinstatement(uint8_t itsClusterMap) const
  {
    constexpr uint8_t bit = 1;
    for (const auto& kITSrequirement : cutMyRequiredITSHits) {
      auto hits = std::count_if(kITSrequirement.second.begin(), kITSrequirement.second.end(), [&](auto&& requiredLayer) { return itsClusterMap & (bit << requiredLayer); });
      if ((kITSrequirement.first == -1) && (hits > 0)) {
        return false; // no hits were required in specified layers
      } else if (hits < kITSrequirement.first) {
        return false; // not enough hits found in specified layers
      }
    }
    return true;
  }

  template <typename T>
  bool isGlobalTrackReinstatement(T const& track)
  {
    // kInAcceptance copy
    if (track.pt() < cutGlobalTrack.cutMinPt || track.pt() > cutGlobalTrack.cutMaxPt)
      return false;
    if (eta(track.px(), track.py(), track.pz()) < cutGlobalTrack.cutMinEta || eta(track.px(), track.py(), track.pz()) > cutGlobalTrack.cutMaxEta)
      return false;
    // kPrimaryTracks
    // GoldenChi2 cut is only for Run 2
    if (std::abs(track.dcaZ()) > cutGlobalTrack.cutMaxDCAz)
      return false;
    if (cutGlobalTrack.applyPtDependentDCAxy) {
      float maxDCA = 0.0182f + 0.0350f / std::pow(track.pt(), 1.01f);
      if (std::abs(track.dcaXY()) > maxDCA)
        return false;
    } else {
      if (std::abs(track.dcaXY()) > cutGlobalTrack.cutMaxDCAxy)
        return false;
    }
    // kQualityTrack
    // TrackType is always 1 as per definition of processed Run3 AO2Ds
    // ITS
    if (cutGlobalTrack.cutHasITS && !track.hasITS())
      return false; // ITS refit
    if (track.itsNCls() < cutGlobalTrack.cutMinITSnCls)
      return false;
    if (track.itsChi2NCl() > cutGlobalTrack.cutMaxITSchi2)
      return false;
    if (!isFulfillsITSHitRequirementsReinstatement(track.itsClusterMap()))
      return false;
    //  TPC
    if (cutGlobalTrack.cutHasTPC && !track.hasTPC())
      return false; // TPC refit
    if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < cutGlobalTrack.cutMinTPCnCls)
      return false; // tpcNClsFound()
    if (track.tpcNClsCrossedRows() < cutGlobalTrack.cutMinTPCnClsXrows)
      return false;
    if ((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())) < cutGlobalTrack.cutMinTPCnClsXrowsOverNcls)
      return false;
    if (track.tpcChi2NCl() > cutGlobalTrack.cutMaxTPCchi2)
      return false; // TPC chi2

    return true;
  }

  template <typename T>
  bool reinstallRun2JpsiTrackSelection(T const& track)
  {
    // kInAcceptance copy
    if (eta(track.px(), track.py(), track.pz()) < -0.8 || eta(track.px(), track.py(), track.pz()) > 0.8)
      return false;
    // kPrimaryTracks
    if (std::abs(track.dcaZ()) > 2.0)
      return false;
    float maxDCA = 0.0105f + 0.0350f / std::pow(track.pt(), 1.1f);
    if (std::abs(track.dcaXY()) > maxDCA)
      return false;
    // kQualityTrack
    // ITS
    if (!track.hasITS())
      return false; // ITS refit
    if (track.itsChi2NCl() > 36.)
      return false;
    if (!isFulfillsITSHitRequirementsReinstatement(track.itsClusterMap()))
      return false;
    // TPC
    if (!track.hasTPC())
      return false; // TPC refit
    if ((track.tpcNClsFindable() - track.tpcNClsFindableMinusFound()) < 70)
      return false; // tpcNClsFound()
    if ((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())) < 0.8)
      return false;
    if (track.tpcChi2NCl() > 4.)
      return false;
    // TOF
    if (!track.hasTOF())
      return false;

    return true;
  }

  template <typename C, typename T>
  bool reinstallRun2JpsiEventSelection(C const& collision, T const& trk1, T const& trk2, float rapMother, float aco)
  {
    // tracks
    if (!reinstallRun2JpsiTrackSelection(trk1))
      return false;
    if (!reinstallRun2JpsiTrackSelection(trk2))
      return false;
    if (trk1.sign() * trk2.sign() > 0)
      return false; // opposite sign
    if ((trk1.tpcNSigmaMu() * trk1.tpcNSigmaMu() + trk2.tpcNSigmaMu() * trk2.tpcNSigmaMu()) >
        (trk1.tpcNSigmaEl() * trk1.tpcNSigmaEl() + trk2.tpcNSigmaEl() * trk2.tpcNSigmaEl()))
      return false; // definitely muon than electron
    // event
    if (collision.posZ() > 15.)
      return false;
    if (rapMother < -0.8 || rapMother > 0.8)
      return false;
    if (aco > 4 * o2::constants::math::PI / 5) // max opening angle 144 degrees (I hope, check)
      return false;

    return true;
  }

  float getPhiModN(float phimodn, int sign, int fieldpolarity)
  {
    if (fieldpolarity < 0) // for negative polarity field
      phimodn = o2::constants::math::TwoPI - phimodn;
    if (sign < 0) // for negative charge
      phimodn = o2::constants::math::TwoPI - phimodn;
    phimodn += o2::constants::math::PI / 18.0; // to center gap in the middle
    return std::fmod(phimodn, o2::constants::math::PI / 9.0);
  }

  template <typename C>
  bool isGoodFITtime(C const& coll, float maxFITtime)
  {

    // FTOA
    if ((std::abs(coll.timeFT0A()) > maxFITtime) && coll.timeFT0A() > -998.)
      return false;

    // FTOC
    if ((std::abs(coll.timeFT0C()) > maxFITtime) && coll.timeFT0A() > -998.)
      return false;

    return true;
  }

  template <typename T>
  bool isElectronCandidate(T const& electronCandidate)
  // Loose criterium to find electron-like particle
  // Requiring TOF to avoid double-counting pions/electrons
  {
    if (electronCandidate.tpcNSigmaEl() < -2.0 || electronCandidate.tpcNSigmaEl() > 4.0)
      return false;
    if (!electronCandidate.hasTOF())
      return false;
    return true;
  }

  template <typename T>
  bool isMuPionCandidate(T const& muPionCandidate)
  // Loose criterium to find muon/pion-like particle
  {
    if (muPionCandidate.tpcNSigmaMu() < -5.0 || muPionCandidate.tpcNSigmaMu() > 5.0)
      return false;
    if (muPionCandidate.tpcNSigmaPi() < -5.0 || muPionCandidate.tpcNSigmaPi() > 5.0)
      return false;
    return true;
  }

  template <typename T>
  bool selectedGoodElectron(T const& electronCandidate)
  {
    if (electronCandidate.tpcNSigmaEl() < cutTauEvent.cutMaxElectronNsigmaEl || electronCandidate.tpcNSigmaEl() > cutTauEvent.cutMinElectronNsigmaEl)
      return false;
    if (electronCandidate.tpcNSigmaPi() > cutTauEvent.cutMinElectronNsigmaPi && electronCandidate.tpcNSigmaPi() < cutTauEvent.cutMaxElectronNsigmaPi)
      return false;
    if (cutTauEvent.cutElectronHasTOF && !electronCandidate.hasTOF())
      return false;
    if (electronCandidate.hasTOF()) {
      if (electronCandidate.tofNSigmaPr() > cutTauEvent.cutMinElectronNsigmaPr && electronCandidate.tofNSigmaPr() < cutTauEvent.cutMaxElectronNsigmaPr)
        return false;
      if (momentum(electronCandidate.px(), electronCandidate.py(), electronCandidate.pz()) < 1.0) {
        if (electronCandidate.tofNSigmaKa() > cutTauEvent.cutMinElectronTofNsigmaKa && electronCandidate.tofNSigmaKa() < cutTauEvent.cutMaxElectronTofNsigmaKa)
          return false;
      } else {
        if (electronCandidate.tpcNSigmaKa() > cutTauEvent.cutMinElectronNsigmaKa && electronCandidate.tpcNSigmaKa() < cutTauEvent.cutMaxElectronNsigmaKa)
          return false;
      }
    }
    return true;
  }

  template <typename T>
  bool selectedGoodPion(T const& pionCandidate)
  {
    if (pionCandidate.tpcNSigmaPi() < cutTauEvent.cutMaxPionNsigmaPi || pionCandidate.tpcNSigmaPi() > cutTauEvent.cutMinPionNsigmaPi)
      return false;
    if (cutTauEvent.cutPionHasTOF && !pionCandidate.hasTOF())
      return false;
    if (pionCandidate.hasTOF()) {
      if (pionCandidate.tofNSigmaPr() > cutTauEvent.cutMinElectronNsigmaPr && pionCandidate.tofNSigmaPr() < cutTauEvent.cutMaxElectronNsigmaPr)
        return false;
      if (momentum(pionCandidate.px(), pionCandidate.py(), pionCandidate.pz()) < 1.0) {
        if (pionCandidate.tofNSigmaKa() > cutTauEvent.cutMinElectronTofNsigmaKa && pionCandidate.tofNSigmaKa() < cutTauEvent.cutMaxElectronTofNsigmaKa)
          return false;
      }
    }
    return true;
  }

  template <typename T>
  bool selectedTauEvent(T const& trkDaug1, T const& trkDaug2)
  {
    TLorentzVector mother, daug[2], motherOfPions, pion[2];
    daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
    daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
    if (cutTauEvent.useThresholdsPID) {
      if (isElectronCandidate(trkDaug1))
        daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(MassElectron, trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      if (isElectronCandidate(trkDaug2))
        daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(MassElectron, trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      if (isMuPionCandidate(trkDaug1))
        daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(MassPionCharged, trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      if (isMuPionCandidate(trkDaug2))
        daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(MassPionCharged, trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
    }
    mother = daug[0] + daug[1];
    pion[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(MassPionCharged, trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
    pion[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(MassPionCharged, trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
    motherOfPions = pion[0] + pion[1];

    int enumTrk1 = (cutTauEvent.useThresholdsPID ? (isElectronCandidate(trkDaug1) ? P_ELECTRON : P_PION) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)));
    if (cutTauEvent.cutOppositeCharge && (trkDaug1.sign() * trkDaug2.sign() > 0))
      return false;
    if (calculateAcoplanarity(daug[0].Phi(), daug[1].Phi()) > cutTauEvent.cutMaxAcoplanarity)
      return false;
    if (calculateAcoplanarity(daug[0].Phi(), daug[1].Phi()) < cutTauEvent.cutMinAcoplanarity)
      return false;
    bool goodElectron = (enumTrk1 == P_ELECTRON) ? selectedGoodElectron(trkDaug1) : selectedGoodElectron(trkDaug2);
    if (cutTauEvent.cutGoodElectron && !goodElectron)
      return false;
    bool goodMupion = ((enumTrk1 == P_MUON) || (enumTrk1 == P_PION)) ? selectedGoodPion(trkDaug1) : selectedGoodPion(trkDaug2);
    if (cutTauEvent.cutGoodMupion && !goodMupion)
      return false;
    if (cutTauEvent.cutOutRho && (motherOfPions.M() > cutTauEvent.cutMinRhoMass && motherOfPions.M() < cutTauEvent.cutMaxRhoMass))
      return false;
    if (cutTauEvent.cutOnRho && (motherOfPions.M() > cutTauEvent.cutMaxRhoMass || motherOfPions.M() < cutTauEvent.cutMinRhoMass))
      return false;
    return true;
  }

  template <typename Ts>
  void fillHistograms(Ts const& reconstructedBarrelTracks)
  {

    histos.get<TH1>(HIST("Events/hCountCollisions"))->Fill(1);
    histos.get<TH1>(HIST("Events/hNreconstructedTracks"))->Fill(reconstructedBarrelTracks.size());

    // Loop over tracks without selections
    for (const auto& track : reconstructedBarrelTracks) {
      float trkPx = track.px();
      float trkPy = track.py();
      float trkPz = track.pz();
      //          histos.get<TH1>(HIST("Tracks/raw/hTrackZ"))->Fill(track.z());
      histos.get<TH1>(HIST("Tracks/raw/hTrackP"))->Fill(momentum(trkPx, trkPy, trkPz));
      histos.get<TH1>(HIST("Tracks/raw/hTrackPt"))->Fill(track.pt());
      histos.get<TH1>(HIST("Tracks/raw/hTrackPhi"))->Fill(phi(trkPx, trkPy));
      histos.get<TH1>(HIST("Tracks/raw/hTrackEta"))->Fill(eta(trkPx, trkPy, trkPz));
      histos.get<TH2>(HIST("Tracks/raw/hTrackPtvsModPhi"))->Fill(track.pt(), std::fmod(phi(trkPx, trkPy), o2::constants::math::PI / 9));
      if (track.hasTOF())
        histos.get<TH2>(HIST("Tracks/raw/hTrackPtvsModPhiTOF"))->Fill(track.pt(), std::fmod(phi(trkPx, trkPy), o2::constants::math::PI / 9));
      histos.get<TH1>(HIST("Tracks/raw/hTrackDcaXY"))->Fill(track.dcaXY());
      histos.get<TH2>(HIST("Tracks/raw/hTrackPtvsDcaXY"))->Fill(track.pt(), track.dcaXY());
      histos.get<TH1>(HIST("Tracks/raw/hTrackDcaZ"))->Fill(track.dcaZ());
      histos.get<TH1>(HIST("Tracks/raw/ITS/itsNCls"))->Fill(track.itsNCls());
      histos.get<TH1>(HIST("Tracks/raw/ITS/itsChi2NCl"))->Fill(track.itsChi2NCl());
      histos.get<TH1>(HIST("Tracks/raw/TPC/tpcNClsFindable"))->Fill(track.tpcNClsFindable());
      histos.get<TH1>(HIST("Tracks/raw/TPC/tpcNClsFound"))->Fill(track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
      histos.get<TH1>(HIST("Tracks/raw/TPC/tpcCrossedRows"))->Fill(track.tpcNClsCrossedRows());
      histos.get<TH1>(HIST("Tracks/raw/TPC/tpcCrossedRowsOverFindableCls"))->Fill((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())));
      histos.get<TH1>(HIST("Tracks/raw/TPC/tpcChi2NCl"))->Fill(track.tpcChi2NCl());
    } // Loop over tracks without selections

    int countPVGT = 0;
    int countPVGTselected = 0;
    int countPVGTelectrons = 0;
    int countPVGTmuons = 0;
    int countPVGTpions = 0;
    int countPVGTothers = 0;
    int countTOFtracks = 0;
    int countPVGTelmupiAlt = 0;
    int countPVGTelectronsAlt = 0;
    int countPVGTmupionsAlt = 0;
    std::vector<int> vecPVidx;
    std::vector<int> vecPVnewPIDidx;
    // Loop over tracks with selections
    for (const auto& track : reconstructedBarrelTracks) {
      if (!track.isPVContributor())
        continue;
      if (cutGlobalTrack.applyGlobalTrackSelection) {
        if (isGlobalTrackReinstatement(track) != 1)
          continue;
      }
      countPVGT++;
      float trkPx = track.px();
      float trkPy = track.py();
      float trkPz = track.pz();
      //          histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackZ"))->Fill(track.z());
      histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackP"))->Fill(momentum(trkPx, trkPy, trkPz));
      histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackPt"))->Fill(track.pt());
      histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackPhi"))->Fill(phi(trkPx, trkPy));
      histos.get<TH2>(HIST("Tracks/GoodTrack/hTrackPtvsModPhi"))->Fill(track.pt(), std::fmod(phi(trkPx, trkPy), o2::constants::math::PI / 9));
      if (track.hasTOF())
        histos.get<TH2>(HIST("Tracks/GoodTrack/hTrackPtvsModPhiTOF"))->Fill(track.pt(), std::fmod(phi(trkPx, trkPy), o2::constants::math::PI / 9));
      histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackEta"))->Fill(eta(trkPx, trkPy, trkPz));
      histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackDcaXY"))->Fill(track.dcaXY());
      histos.get<TH2>(HIST("Tracks/GoodTrack/hTrackPtvsDcaXY"))->Fill(track.pt(), track.dcaXY());
      histos.get<TH1>(HIST("Tracks/GoodTrack/hTrackDcaZ"))->Fill(track.dcaZ());
      histos.get<TH1>(HIST("Tracks/GoodTrack/ITS/itsNCls"))->Fill(track.itsNCls());
      histos.get<TH1>(HIST("Tracks/GoodTrack/ITS/itsChi2NCl"))->Fill(track.itsChi2NCl());
      histos.get<TH1>(HIST("Tracks/GoodTrack/TPC/tpcNClsFindable"))->Fill(track.tpcNClsFindable());
      histos.get<TH1>(HIST("Tracks/GoodTrack/TPC/tpcNClsFound"))->Fill(track.tpcNClsFindable() - track.tpcNClsFindableMinusFound());
      histos.get<TH1>(HIST("Tracks/GoodTrack/TPC/tpcCrossedRows"))->Fill(track.tpcNClsCrossedRows());
      histos.get<TH1>(HIST("Tracks/GoodTrack/TPC/tpcCrossedRowsOverFindableCls"))->Fill((static_cast<float>(track.tpcNClsCrossedRows()) / static_cast<float>(track.tpcNClsFindable())));
      histos.get<TH1>(HIST("Tracks/GoodTrack/TPC/tpcChi2NCl"))->Fill(track.tpcChi2NCl());
      if (track.hasTOF())
        countTOFtracks++;
      int hypothesisID = testPIDhypothesis(track, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC);
      if (hypothesisID == P_ELECTRON || hypothesisID == P_MUON || hypothesisID == P_PION) {
        countPVGTselected++;
        vecPVidx.push_back(track.index());
        if (hypothesisID == P_ELECTRON) {
          countPVGTelectrons++;
        } else if (hypothesisID == P_MUON) {
          countPVGTmuons++;
        } else {
          countPVGTpions++;
        }
      } else {
        countPVGTothers++;
      }
      // alternative selection
      if (isElectronCandidate(track)) {
        countPVGTelmupiAlt++;
        countPVGTelectronsAlt++;
        vecPVnewPIDidx.push_back(track.index());
      }
      if (isMuPionCandidate(track)) {
        countPVGTelmupiAlt++;
        countPVGTmupionsAlt++;
        vecPVnewPIDidx.push_back(track.index());
      }
    } // Loop over tracks with selections

    histos.get<TH1>(HIST("Events/hNreconstructedPVGT"))->Fill(countPVGT);
    histos.get<TH1>(HIST("Events/hNreconstructedNotPVGT"))->Fill(reconstructedBarrelTracks.size() - countPVGT);
    histos.get<TH1>(HIST("Events/hNreconstructedPVGTelectrons"))->Fill(countPVGTelectrons);
    histos.get<TH1>(HIST("Events/hNreconstructedPVGTmuons"))->Fill(countPVGTmuons);
    histos.get<TH1>(HIST("Events/hNreconstructedPVGTpions"))->Fill(countPVGTpions);
    histos.get<TH1>(HIST("Events/hNreconstructedPVGTothers"))->Fill(countPVGTothers);

    bool isTwoSelectedTracks = (cutTauEvent.useThresholdsPID ? countPVGTelmupiAlt == 2 : countPVGTselected == 2);
    bool isElMuPion = (cutTauEvent.useThresholdsPID ? (countPVGTelectronsAlt == 1 && countPVGTmupionsAlt == 1) : ((countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)));
    if (isTwoSelectedTracks && doTwoTracks) {
      TLorentzVector mother, daug[2], motherOfPions, pion[2], motherOfMuons, muon[2];

      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(cutTauEvent.useThresholdsPID ? vecPVnewPIDidx[0] : vecPVidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(cutTauEvent.useThresholdsPID ? vecPVnewPIDidx[1] : vecPVidx[1]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      if (cutTauEvent.useThresholdsPID) {
        if (isElectronCandidate(trkDaug1))
          daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(MassElectron, trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
        if (isElectronCandidate(trkDaug2))
          daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(MassElectron, trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
        if (isMuPionCandidate(trkDaug1))
          daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(MassPionCharged, trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
        if (isMuPionCandidate(trkDaug2))
          daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(MassPionCharged, trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      }
      mother = daug[0] + daug[1];
      pion[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(MassPionCharged, trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      pion[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(MassPionCharged, trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      motherOfPions = pion[0] + pion[1];
      muon[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(MassMuon, trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      muon[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(MassMuon, trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      motherOfMuons = muon[0] + muon[1];
      const auto acoplanarity = calculateAcoplanarity(daug[0].Phi(), daug[1].Phi());
      const auto collinearity = calculateCollinearity(daug[0].Eta(), daug[1].Eta(), daug[0].Phi(), daug[1].Phi());
      if (cutTauEvent.applyTauEventSelection && !selectedTauEvent(trkDaug1, trkDaug2)) {
        return;
      }

      histos.get<TH1>(HIST("EventTwoTracks/hInvariantMass"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventTwoTracks/hInvariantMassWide"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventTwoTracks/hInvariantMassWideAllPionMass"))->Fill(motherOfPions.M());
      histos.get<TH1>(HIST("EventTwoTracks/hAcoplanarity"))->Fill(acoplanarity);
      histos.get<TH1>(HIST("EventTwoTracks/hCollinearity"))->Fill(collinearity);
      histos.get<TH1>(HIST("EventTwoTracks/hMotherP"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventTwoTracks/hMotherPwide"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventTwoTracks/hMotherPt"))->Fill(mother.Pt());
      histos.get<TH1>(HIST("EventTwoTracks/hMotherPhi"))->Fill(mother.Phi());
      histos.get<TH1>(HIST("EventTwoTracks/hMotherRapidity"))->Fill(mother.Rapidity());
      histos.get<TH2>(HIST("EventTwoTracks/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
      histos.get<TH1>(HIST("EventTwoTracks/hDaughtersEnergyAsymmetry"))->Fill((daug[0].E() - daug[1].E()) / (daug[0].E() + daug[1].E()));
      histos.get<TH1>(HIST("EventTwoTracks/hDaughtersMomentaAsymmetry"))->Fill((daug[0].P() - daug[1].P()) / (daug[0].P() + daug[1].P()));
      histos.get<TH1>(HIST("EventTwoTracks/hDaughtersPtAsymmetry"))->Fill((daug[0].Pt() - daug[1].Pt()) / (daug[0].Pt() + daug[1].Pt()));
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersEnergyFractions"))->Fill(daug[0].E() / (daug[0].E() + daug[1].E()), daug[1].E() / (daug[0].E() + daug[1].E()));
      if (motherOfPions.Pt() < 0.2) {
        histos.get<TH1>(HIST("EventTwoTracks/hInvariantMassWideAllPionMassPtCut"))->Fill(motherOfPions.M());
      }
      if (countTOFtracks == 2) {
        histos.get<TH1>(HIST("EventTwoTracks/hInvariantMassWideAllPionMassTOF"))->Fill(motherOfPions.M());
      }
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPvsITSclusterSize"))->Fill(getAvgITSClSize(trkDaug1), trkDaug1.sign() * daug[0].P());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPvsITSclusterSize"))->Fill(getAvgITSClSize(trkDaug2), trkDaug2.sign() * daug[1].P());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPvsITSclusterSizeXcos"))->Fill(getAvgITSClSize(trkDaug1) * getCosLambda(trkDaug1), trkDaug1.sign() * daug[0].P());
      histos.get<TH2>(HIST("EventTwoTracks/hDaughtersPvsITSclusterSizeXcos"))->Fill(getAvgITSClSize(trkDaug2) * getCosLambda(trkDaug2), trkDaug2.sign() * daug[1].P());

      // ee, mm, em, pp, ep, mp, pppp, eppp, mppp, pppppp
      if (countPVGTelectrons == 2) {
        histos.get<TH1>(HIST("Events/hChannels"))->Fill(CH_EE);
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhi"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhi"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hLeadingP"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].P() : daug[1].P()));
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hLeadingPwide"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].P() : daug[1].P()));
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hLeadingPt"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Pt() : daug[1].Pt()));
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hLeadingPhi"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Phi() : daug[1].Phi()));
        histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hLeadingRapidity"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Rapidity() : daug[1].Rapidity()));
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hLeadingPvsOtherP"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].P() : daug[1].P()), ((daug[0].P() > daug[1].P()) ? daug[1].P() : daug[0].P()));
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hLeadingPwideVsOtherPwide"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].P() : daug[1].P()), ((daug[0].P() > daug[1].P()) ? daug[1].P() : daug[0].P()));
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hLeadingPtVsOtherPt"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Pt() : daug[1].Pt()), ((daug[0].P() > daug[1].P()) ? daug[1].Pt() : daug[0].Pt()));
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hLeadingPhiVsOtherPhi"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Phi() : daug[1].Phi()), ((daug[0].P() > daug[1].P()) ? daug[1].Phi() : daug[0].Phi()));
        histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hLeadingRapVsOtherRap"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Rapidity() : daug[1].Rapidity()), ((daug[0].P() > daug[1].P()) ? daug[1].Rapidity() : daug[0].Rapidity()));
        if (mother.Pt() < 0.2) {
          histos.get<TH1>(HIST("EventTwoTracks/TwoElectrons/hInvariantMassWidePtCut"))->Fill(mother.M());
          histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiPtCut"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
          histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiPtCut"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
          if (countTOFtracks == 2) {
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiPtCutTOF"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiPtCutTOF"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
          }
        }
        if (countTOFtracks == 2) {
          histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiTOF"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
          histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/hDaughtersPtvsModPhiTOF"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
        }
      }
      if (countPVGTmuons == 2) {
        histos.get<TH1>(HIST("Events/hChannels"))->Fill(CH_MUMU);
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhi"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
        histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhi"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
        if (mother.Pt() < 0.2) {
          histos.get<TH1>(HIST("EventTwoTracks/TwoMuons/hInvariantMassWidePtCut"))->Fill(mother.M());
          histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiPtCut"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
          histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiPtCut"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
          if (countTOFtracks == 2) {
            histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiPtCutTOF"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
            histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiPtCutTOF"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
          }
        }
        if (countTOFtracks == 2) {
          histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiTOF"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
          histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/hDaughtersPtvsModPhiTOF"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
        }
      }
      if (countPVGTelectrons == 1 && countPVGTmuons == 1) {
        histos.get<TH1>(HIST("Events/hChannels"))->Fill(CH_EMU);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuon/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
      }
      if (countPVGTpions == 2) {
        histos.get<TH1>(HIST("Events/hChannels"))->Fill(CH_PIPI);
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPtvsModPhi"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
        histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPtvsModPhi"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
        if (mother.Pt() < 0.2) {
          histos.get<TH1>(HIST("EventTwoTracks/TwoPions/hInvariantMassWidePtCut"))->Fill(mother.M());
          histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiPtCut"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
          histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiPtCut"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
          if (countTOFtracks == 2) {
            histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiPtCutTOF"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
            histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiPtCutTOF"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
          }
        }
        if (countTOFtracks == 2) {
          histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiTOF"))->Fill(daug[0].P(), getPhiModN(daug[0].Phi(), trkDaug1.sign(), 1));
          histos.get<TH2>(HIST("EventTwoTracks/TwoPions/hDaughtersPtvsModPhiTOF"))->Fill(daug[1].P(), getPhiModN(daug[1].Phi(), trkDaug2.sign(), 1));
        }
      }
      if (countPVGTelectrons == 1 && countPVGTpions == 1) {
        histos.get<TH1>(HIST("Events/hChannels"))->Fill(CH_EPI);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronPion/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
      }
      if (countPVGTpions == 1 && countPVGTmuons == 1) {
        histos.get<TH1>(HIST("Events/hChannels"))->Fill(CH_MUPI);
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/MuonPion/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        if (mother.Pt() < 0.2) {
          histos.get<TH1>(HIST("EventTwoTracks/MuonPion/hInvariantMassWidePtCut"))->Fill(mother.M());
        }
      }
      if (isElMuPion) {
        double electronPt = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt();
        double electronP = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[0].P() : daug[1].P();
        double electronE = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[0].E() : daug[1].E();
        double mupionPt = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[1].Pt() : daug[0].Pt();
        double mupionP = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[1].P() : daug[0].P();
        double mupionE = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[1].E() : daug[0].E();

        TLorentzVector motherOfPiKaon, kaon;
        if (isElectronCandidate(trkDaug1)) {
          kaon.SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(MassKaonCharged, trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
          motherOfPiKaon = kaon + daug[1];
        } else {
          kaon.SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(MassKaonCharged, trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
          motherOfPiKaon = daug[0] + kaon;
        }
        histos.get<TH1>(HIST("Events/hChannels"))->Fill(CH_EMUPI);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/PionsSelection/hInvariantMass"))->Fill(motherOfPions.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/PionsSelection/hInvariantMassWide"))->Fill(motherOfPions.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/KstarSelection/hInvariantMass"))->Fill(motherOfPiKaon.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/KstarSelection/hInvariantMassWide"))->Fill(motherOfPiKaon.M());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/KstarSelection/hMotherMassVsPt"))->Fill(motherOfPiKaon.M(), motherOfPiKaon.Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/KstarSelection/hMotherMassVsElectronP"))->Fill(motherOfPiKaon.M(), electronP);
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/KstarSelection/hIMKvsIMe"))->Fill(motherOfPiKaon.M(), mother.M());
        if (electronPt < cutTauEvent.cutElectronPt) {
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutLow/hInvariantMass"))->Fill(motherOfPiKaon.M());
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutLow/hInvariantMassWide"))->Fill(motherOfPiKaon.M());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutLow/hMotherMassVsPt"))->Fill(motherOfPiKaon.M(), motherOfPiKaon.Pt());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutLow/hIMKvsIMe"))->Fill(motherOfPiKaon.M(), mother.M());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutLow/hElectronPwideVsOtherPwide"))->Fill(electronP, mupionP);
        } else {
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutHigh/hInvariantMass"))->Fill(motherOfPiKaon.M());
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutHigh/hInvariantMassWide"))->Fill(motherOfPiKaon.M());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutHigh/hMotherMassVsPt"))->Fill(motherOfPiKaon.M(), motherOfPiKaon.Pt());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutHigh/hIMKvsIMe"))->Fill(motherOfPiKaon.M(), mother.M());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/KstarSelection/eMomentaCutHigh/hElectronPwideVsOtherPwide"))->Fill(electronP, mupionP);
        }
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hCollinearity"))->Fill(collinearity);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hMotherMassVsElectronP"))->Fill(mother.M(), electronP);
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hMotherMassVsAcoplanarity"))->Fill(mother.M(), acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hElectronPt"))->Fill(electronPt);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hElectronPtWide"))->Fill(electronPt);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersEnergyAsymmetry"))->Fill((daug[0].E() - daug[1].E()) / (daug[0].E() + daug[1].E()));
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersMomentaAsymmetry"))->Fill((daug[0].P() - daug[1].P()) / (daug[0].P() + daug[1].P()));
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersPtAsymmetry"))->Fill((daug[0].Pt() - daug[1].Pt()) / (daug[0].Pt() + daug[1].Pt()));
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hDaughtersEnergyFractions"))->Fill(electronE / (electronE + mupionE), mupionE / (electronE + mupionE));
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hElectronPvsOtherP"))->Fill(electronP, mupionP);
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hElectronPwideVsOtherPwide"))->Fill(electronP, mupionP);
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hElectronPvsAcoplanarity"))->Fill(electronP, acoplanarity);
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hOtherPvsAcoplanarity"))->Fill(mupionP, acoplanarity);
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hElectronPtVsOtherPt"))->Fill(electronPt, mupionPt);
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hElectronPhiVsOtherPhi"))->Fill(isElectronCandidate(trkDaug1) ? daug[0].Phi() : daug[1].Phi(), isElectronCandidate(trkDaug1) ? daug[1].Phi() : daug[0].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/hElectronRapVsOtherRap"))->Fill(isElectronCandidate(trkDaug1) ? daug[0].Rapidity() : daug[1].Rapidity(), isElectronCandidate(trkDaug1) ? daug[1].Rapidity() : daug[0].Rapidity());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(0);
        if (mother.Pt() < 9.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(1);
        if (mother.Pt() < 8.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(2);
        if (mother.Pt() < 7.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(3);
        if (mother.Pt() < 6.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(4);
        if (mother.Pt() < 5.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(5);
        if (mother.Pt() < 4.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(6);
        if (mother.Pt() < 3.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(7);
        if (mother.Pt() < 2.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(8);
        if (mother.Pt() < 1.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(9);
        if (electronPt > 0.1 && electronPt < 1.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(10);
        if (electronPt > 1. && electronPt < 2.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(11);
        if (electronPt > 2. && electronPt < 100.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronMuPi/hNeventsPtCuts"))->Fill(12);
      }
      if ((countPVGTelectrons == 2) || (countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)) {
        double electronPt = (enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt();
        if (countPVGTelectrons == 2)
          electronPt = (daug[0].Pt() > daug[1].Pt()) ? daug[0].Pt() : daug[1].Pt();
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hAcoplanarity"))->Fill(acoplanarity);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hElectronPt"))->Fill(electronPt);
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hElectronPtWide"))->Fill(electronPt);
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hDaughtersP"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hDaughtersPwide"))->Fill(daug[0].P(), daug[1].P());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hDaughtersPt"))->Fill(daug[0].Pt(), daug[1].Pt());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hDaughtersPhi"))->Fill(daug[0].Phi(), daug[1].Phi());
        histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hDaughtersRapidity"))->Fill(daug[0].Rapidity(), daug[1].Rapidity());
        if (countPVGTelectrons == 2) {
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronPvsOtherP"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].P() : daug[1].P()), ((daug[0].P() > daug[1].P()) ? daug[1].P() : daug[0].P()));
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronPwideVsOtherPwide"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].P() : daug[1].P()), ((daug[0].P() > daug[1].P()) ? daug[1].P() : daug[0].P()));
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronPtVsOtherPt"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Pt() : daug[1].Pt()), ((daug[0].P() > daug[1].P()) ? daug[1].Pt() : daug[0].Pt()));
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronPhiVsOtherPhi"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Phi() : daug[1].Phi()), ((daug[0].P() > daug[1].P()) ? daug[1].Phi() : daug[0].Phi()));
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronRapVsOtherRap"))->Fill(((daug[0].P() > daug[1].P()) ? daug[0].Rapidity() : daug[1].Rapidity()), ((daug[0].P() > daug[1].P()) ? daug[1].Rapidity() : daug[0].Rapidity()));
        } else {
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronPvsOtherP"))->Fill((enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[0].P() : daug[1].P(), (enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[1].P() : daug[0].P());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronPwideVsOtherPwide"))->Fill((enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[0].P() : daug[1].P(), (enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[1].P() : daug[0].P());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronPtVsOtherPt"))->Fill((enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt(), (enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[1].Pt() : daug[0].Pt());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronPhiVsOtherPhi"))->Fill((enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[0].Phi() : daug[1].Phi(), (enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[1].Phi() : daug[0].Phi());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/hElectronRapVsOtherRap"))->Fill((enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[0].Rapidity() : daug[1].Rapidity(), (enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[1].Rapidity() : daug[0].Rapidity());
        }
        histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(0);
        if (mother.Pt() < 9.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(1);
        if (mother.Pt() < 8.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(2);
        if (mother.Pt() < 7.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(3);
        if (mother.Pt() < 6.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(4);
        if (mother.Pt() < 5.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(5);
        if (mother.Pt() < 4.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(6);
        if (mother.Pt() < 3.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(7);
        if (mother.Pt() < 2.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(8);
        if (mother.Pt() < 1.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(9);
        if (electronPt > 0.1 && electronPt < 1.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(10);
        if (electronPt > 1. && electronPt < 2.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(11);
        if (electronPt > 2. && electronPt < 100.)
          histos.get<TH1>(HIST("EventTwoTracks/ElectronOther/hNeventsPtCuts"))->Fill(12);
      }
    } else if (countPVGTselected == 4 && doFourTracks) {

      TLorentzVector mother, daug[4];
      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(vecPVidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(vecPVidx[1]);
      const auto& trkDaug3 = reconstructedBarrelTracks.iteratorAt(vecPVidx[2]);
      const auto& trkDaug4 = reconstructedBarrelTracks.iteratorAt(vecPVidx[3]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      daug[2].SetPxPyPzE(trkDaug3.px(), trkDaug3.py(), trkDaug3.pz(), energy(pdg->Mass(trackPDG(trkDaug3, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug3.px(), trkDaug3.py(), trkDaug3.pz()));
      daug[3].SetPxPyPzE(trkDaug4.px(), trkDaug4.py(), trkDaug4.pz(), energy(pdg->Mass(trackPDG(trkDaug4, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug4.px(), trkDaug4.py(), trkDaug4.pz()));
      mother = daug[0] + daug[1] + daug[2] + daug[3];

      histos.get<TH1>(HIST("EventFourTracks/hInvariantMass"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventFourTracks/hInvariantMassWide"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventFourTracks/hMotherP"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventFourTracks/hMotherPwide"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventFourTracks/hMotherPt"))->Fill(mother.Pt());
      histos.get<TH1>(HIST("EventFourTracks/hMotherPhi"))->Fill(mother.Phi());
      histos.get<TH1>(HIST("EventFourTracks/hMotherRapidity"))->Fill(mother.Rapidity());
      histos.get<TH2>(HIST("EventFourTracks/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());

      // ee, mm, em, pp, ep, mp, pppp, eppp, mppp, pppppp
      if (countPVGTpions == 4) {
        histos.get<TH1>(HIST("Events/hChannels"))->Fill(CH_FOURPI);
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventFourTracks/WithPion/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventFourTracks/WithPion/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
      }
      if (countPVGTelectrons == 1 && countPVGTpions == 3) {
        histos.get<TH1>(HIST("Events/hChannels"))->Fill(CH_ETHREEPI);
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventFourTracks/WithElectron/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventFourTracks/WithElectron/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
      }
      if (countPVGTpions == 3 && countPVGTmuons == 1) {
        histos.get<TH1>(HIST("Events/hChannels"))->Fill(CH_MUTHREEPI);
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventFourTracks/WithMuon/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventFourTracks/WithMuon/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
      }
    } else if (countPVGTselected == 6 && doSixTracks) {
      TLorentzVector mother, daug[6];
      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(vecPVidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(vecPVidx[1]);
      const auto& trkDaug3 = reconstructedBarrelTracks.iteratorAt(vecPVidx[2]);
      const auto& trkDaug4 = reconstructedBarrelTracks.iteratorAt(vecPVidx[3]);
      const auto& trkDaug5 = reconstructedBarrelTracks.iteratorAt(vecPVidx[4]);
      const auto& trkDaug6 = reconstructedBarrelTracks.iteratorAt(vecPVidx[5]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      daug[2].SetPxPyPzE(trkDaug3.px(), trkDaug3.py(), trkDaug3.pz(), energy(pdg->Mass(trackPDG(trkDaug3, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug3.px(), trkDaug3.py(), trkDaug3.pz()));
      daug[3].SetPxPyPzE(trkDaug4.px(), trkDaug4.py(), trkDaug4.pz(), energy(pdg->Mass(trackPDG(trkDaug4, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug4.px(), trkDaug4.py(), trkDaug4.pz()));
      daug[4].SetPxPyPzE(trkDaug5.px(), trkDaug5.py(), trkDaug5.pz(), energy(pdg->Mass(trackPDG(trkDaug5, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug5.px(), trkDaug5.py(), trkDaug5.pz()));
      daug[5].SetPxPyPzE(trkDaug6.px(), trkDaug6.py(), trkDaug6.pz(), energy(pdg->Mass(trackPDG(trkDaug6, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug6.px(), trkDaug6.py(), trkDaug6.pz()));
      mother = daug[0] + daug[1] + daug[2] + daug[3] + daug[4] + daug[5];

      histos.get<TH1>(HIST("EventSixTracks/hInvariantMass"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventSixTracks/hInvariantMassWide"))->Fill(mother.M());
      histos.get<TH1>(HIST("EventSixTracks/hMotherP"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventSixTracks/hMotherPwide"))->Fill(mother.P());
      histos.get<TH1>(HIST("EventSixTracks/hMotherPt"))->Fill(mother.Pt());
      histos.get<TH1>(HIST("EventSixTracks/hMotherPhi"))->Fill(mother.Phi());
      histos.get<TH1>(HIST("EventSixTracks/hMotherRapidity"))->Fill(mother.Rapidity());
      histos.get<TH2>(HIST("EventSixTracks/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());

      // ee, mm, em, pp, ep, mp, pppp, eppp, mppp, pppppp
      if (countPVGTpions == 6) {
        histos.get<TH1>(HIST("Events/hChannels"))->Fill(CH_SIXPI);
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hInvariantMass"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hInvariantMassWide"))->Fill(mother.M());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hMotherP"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hMotherPwide"))->Fill(mother.P());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hMotherPt"))->Fill(mother.Pt());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hMotherPhi"))->Fill(mother.Phi());
        histos.get<TH1>(HIST("EventSixTracks/SixPions/hMotherRapidity"))->Fill(mother.Rapidity());
        histos.get<TH2>(HIST("EventSixTracks/SixPions/hMotherMassVsPt"))->Fill(mother.M(), mother.Pt());
      }
    } else {
      printDebugMessage("Other particles");
    }

  } // end fillHistograms

  template <typename C, typename Ts>
  void fillPIDhistograms(C const& /*reconstructedCollision*/, Ts const& reconstructedBarrelTracks)
  {

    // Loop over tracks without selections
    for (const auto& track : reconstructedBarrelTracks) {
      float trkPx = track.px();
      float trkPy = track.py();
      float trkPz = track.pz();
      if (track.hasTPC()) {
        //          histos.get<TH2>(HIST("Tracks/raw/PID/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/raw/PID/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/raw/PID/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/raw/PID/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/raw/PID/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/raw/PID/hTPCnSigmaElVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaEl());
        histos.get<TH2>(HIST("Tracks/raw/PID/hTPCnSigmaMuVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaMu());
        histos.get<TH2>(HIST("Tracks/raw/PID/hTPCnSigmaPiVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaPi());
        histos.get<TH2>(HIST("Tracks/raw/PID/hTPCnSigmaKaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaKa());
        histos.get<TH2>(HIST("Tracks/raw/PID/hTPCnSigmaPrVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaPr());
        if (track.hasTOF())
          histos.get<TH2>(HIST("Tracks/raw/PID/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
        if (track.sign() == 1) {
          //                  histos.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
          if (track.hasTOF())
            histos.get<TH2>(HIST("Tracks/raw/PID/PosCharge/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
        } else if (track.sign() == -1) {
          //                  histos.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTPCsignalVsZ"))->Fill(track.z(),track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
          if (track.hasTOF())
            histos.get<TH2>(HIST("Tracks/raw/PID/NegCharge/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
        } else {
          printMediumMessage("Track has no charge");
        }
      }
    } // Loop over tracks without selections

    int countPVGT = 0;
    int countPVGTselected = 0;
    int countPVGTelectrons = 0;
    int countPVGTmuons = 0;
    int countPVGTpions = 0;
    int countPVGTelmupiAlt = 0;
    int countPVGTelectronsAlt = 0;
    int countPVGTmupionsAlt = 0;
    std::vector<int> vecPVidx;
    std::vector<int> vecPVnoPIDidx;
    std::vector<int> vecPVnewPIDidx;
    // Loop over tracks with selections
    for (const auto& track : reconstructedBarrelTracks) {
      if (!track.isPVContributor())
        continue;
      if (cutGlobalTrack.applyGlobalTrackSelection && !isGlobalTrackReinstatement(track)) {
        continue;
      }
      countPVGT++;
      vecPVnoPIDidx.push_back(track.index());
      float trkPx = track.px();
      float trkPy = track.py();
      float trkPz = track.pz();
      if (track.hasTPC()) {
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCnSigmaElVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaEl());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCnSigmaMuVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaMu());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCnSigmaPiVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaPi());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCnSigmaKaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaKa());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/hTPCnSigmaPrVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaPr());
        if (track.hasTOF())
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
        if (track.sign() == 1) {
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/PosCharge/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
          if (track.hasTOF())
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/PosCharge/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
        } else if (track.sign() == -1) {
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/NegCharge/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
          if (track.hasTOF())
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/NegCharge/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
        } else {
          printMediumMessage("Track has no charge");
        }
      }
      int hypothesisID = testPIDhypothesis(track, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC);
      if (hypothesisID == P_ELECTRON || hypothesisID == P_MUON || hypothesisID == P_PION) {
        countPVGTselected++;
        vecPVidx.push_back(track.index());
        if (hypothesisID == P_ELECTRON) {
          countPVGTelectrons++;
          if (!cutTauEvent.useThresholdsPID) {
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCnSigmaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaEl());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsMu"))->Fill(track.tpcNSigmaEl(), track.tpcNSigmaMu());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsPi"))->Fill(track.tpcNSigmaEl(), track.tpcNSigmaPi());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsKa"))->Fill(track.tpcNSigmaEl(), track.tpcNSigmaKa());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsPr"))->Fill(track.tpcNSigmaEl(), track.tpcNSigmaPr());
            if (track.hasTOF()) {
              histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
              histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFnSigmaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofNSigmaEl());
              histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsMu"))->Fill(track.tofNSigmaEl(), track.tofNSigmaMu());
              histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsPi"))->Fill(track.tofNSigmaEl(), track.tofNSigmaPi());
              histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsKa"))->Fill(track.tofNSigmaEl(), track.tofNSigmaKa());
              histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsPr"))->Fill(track.tofNSigmaEl(), track.tofNSigmaPr());
            }
          }
        } else if (hypothesisID == P_MUON) {
          countPVGTmuons++;
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCnSigmaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaMu());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCnSigmaMuVsEl"))->Fill(track.tpcNSigmaMu(), track.tpcNSigmaEl());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCnSigmaMuVsPi"))->Fill(track.tpcNSigmaMu(), track.tpcNSigmaPi());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCnSigmaMuVsKa"))->Fill(track.tpcNSigmaMu(), track.tpcNSigmaKa());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTPCnSigmaMuVsPr"))->Fill(track.tpcNSigmaMu(), track.tpcNSigmaPr());
          if (track.hasTOF()) {
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTOFnSigmaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofNSigmaMu());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTOFnSigmaMuVsEl"))->Fill(track.tofNSigmaMu(), track.tofNSigmaEl());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTOFnSigmaMuVsPi"))->Fill(track.tofNSigmaMu(), track.tofNSigmaPi());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTOFnSigmaMuVsKa"))->Fill(track.tofNSigmaMu(), track.tofNSigmaKa());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Muon/hTOFnSigmaMuVsPr"))->Fill(track.tofNSigmaMu(), track.tofNSigmaPr());
          }
        } else {
          countPVGTpions++;
          if (!cutTauEvent.useThresholdsPID) {
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCnSigmaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaPi());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsEl"))->Fill(track.tpcNSigmaPi(), track.tpcNSigmaEl());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsMu"))->Fill(track.tpcNSigmaPi(), track.tpcNSigmaMu());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsKa"))->Fill(track.tpcNSigmaPi(), track.tpcNSigmaKa());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsPr"))->Fill(track.tpcNSigmaPi(), track.tpcNSigmaPr());
            if (track.hasTOF()) {
              histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
              histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFnSigmaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofNSigmaPi());
              histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsEl"))->Fill(track.tofNSigmaPi(), track.tofNSigmaEl());
              histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsMu"))->Fill(track.tofNSigmaPi(), track.tofNSigmaMu());
              histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsKa"))->Fill(track.tofNSigmaPi(), track.tofNSigmaKa());
              histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsPr"))->Fill(track.tofNSigmaPi(), track.tofNSigmaPr());
            }
          }
        }
      } else {
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Others/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Others/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Others/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
        histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Others/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
        if (track.hasTOF()) {
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Others/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
        }
      }
      // alternative selection
      if (isElectronCandidate(track)) {
        countPVGTelmupiAlt++;
        countPVGTelectronsAlt++;
        vecPVnewPIDidx.push_back(track.index());
        if (cutTauEvent.useThresholdsPID) {
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCnSigmaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaEl());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsMu"))->Fill(track.tpcNSigmaEl(), track.tpcNSigmaMu());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsPi"))->Fill(track.tpcNSigmaEl(), track.tpcNSigmaPi());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsKa"))->Fill(track.tpcNSigmaEl(), track.tpcNSigmaKa());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTPCnSigmaElVsPr"))->Fill(track.tpcNSigmaEl(), track.tpcNSigmaPr());
          if (track.hasTOF()) {
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFnSigmaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofNSigmaEl());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsMu"))->Fill(track.tofNSigmaEl(), track.tofNSigmaMu());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsPi"))->Fill(track.tofNSigmaEl(), track.tofNSigmaPi());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsKa"))->Fill(track.tofNSigmaEl(), track.tofNSigmaKa());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Electron/hTOFnSigmaElVsPr"))->Fill(track.tofNSigmaEl(), track.tofNSigmaPr());
          }
        }
      }
      if (isMuPionCandidate(track)) {
        countPVGTelmupiAlt++;
        countPVGTmupionsAlt++;
        vecPVnewPIDidx.push_back(track.index());
        if (cutTauEvent.useThresholdsPID) {
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsPt"))->Fill(track.pt(), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsEta"))->Fill(eta(trkPx, trkPy, trkPz), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCsignalVsPhi"))->Fill(phi(trkPx, trkPy), track.tpcSignal());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCnSigmaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tpcNSigmaPi());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsEl"))->Fill(track.tpcNSigmaPi(), track.tpcNSigmaEl());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsMu"))->Fill(track.tpcNSigmaPi(), track.tpcNSigmaMu());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsKa"))->Fill(track.tpcNSigmaPi(), track.tpcNSigmaKa());
          histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTPCnSigmaPiVsPr"))->Fill(track.tpcNSigmaPi(), track.tpcNSigmaPr());
          if (track.hasTOF()) {
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFsignalVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofSignal());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFnSigmaVsP"))->Fill(momentum(trkPx, trkPy, trkPz), track.tofNSigmaPi());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsEl"))->Fill(track.tofNSigmaPi(), track.tofNSigmaEl());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsMu"))->Fill(track.tofNSigmaPi(), track.tofNSigmaMu());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsKa"))->Fill(track.tofNSigmaPi(), track.tofNSigmaKa());
            histos.get<TH2>(HIST("Tracks/GoodTrack/PID/Pion/hTOFnSigmaPiVsPr"))->Fill(track.tofNSigmaPi(), track.tofNSigmaPr());
          }
        }
      }

    } // Loop over tracks with selections

    if (countPVGT == 2 && doTwoTracks) {
      TLorentzVector daug[2], pion[2], muon[2];
      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(vecPVnoPIDidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(vecPVnoPIDidx[1]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      if (cutTauEvent.applyTauEventSelection && !selectedTauEvent(trkDaug1, trkDaug2)) {
        return;
      }
      if (cutTauEvent.useThresholdsPID) {
        if (isElectronCandidate(trkDaug1))
          daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(MassElectron, trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
        if (isElectronCandidate(trkDaug2))
          daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(MassElectron, trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
        if (isMuPionCandidate(trkDaug1))
          daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(MassPionCharged, trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
        if (isMuPionCandidate(trkDaug2))
          daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(MassPionCharged, trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      }

      if (trkDaug1.hasTPC()) {
        histos.get<TH2>(HIST("EventTwoTracks/PID/NoPID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        histos.get<TH2>(HIST("EventTwoTracks/PID/NoPID/hTPCnSigmaElVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaEl());
        histos.get<TH2>(HIST("EventTwoTracks/PID/NoPID/hTPCnSigmaMuVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaMu());
        histos.get<TH2>(HIST("EventTwoTracks/PID/NoPID/hTPCnSigmaPiVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaPi());
        histos.get<TH2>(HIST("EventTwoTracks/PID/NoPID/hTPCnSigmaKaVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaKa());
        histos.get<TH2>(HIST("EventTwoTracks/PID/NoPID/hTPCnSigmaPrVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaPr());
      }
      if (trkDaug2.hasTPC()) {
        histos.get<TH2>(HIST("EventTwoTracks/PID/NoPID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        histos.get<TH2>(HIST("EventTwoTracks/PID/NoPID/hTPCnSigmaElVsP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaEl());
        histos.get<TH2>(HIST("EventTwoTracks/PID/NoPID/hTPCnSigmaMuVsP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaMu());
        histos.get<TH2>(HIST("EventTwoTracks/PID/NoPID/hTPCnSigmaPiVsP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaPi());
        histos.get<TH2>(HIST("EventTwoTracks/PID/NoPID/hTPCnSigmaKaVsP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaKa());
        histos.get<TH2>(HIST("EventTwoTracks/PID/NoPID/hTPCnSigmaPrVsP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaPr());
      }
    }

    bool isTwoSelectedTracks = (cutTauEvent.useThresholdsPID ? countPVGTelmupiAlt == 2 : countPVGTselected == 2);
    bool isElMuPion = (cutTauEvent.useThresholdsPID ? (countPVGTelectronsAlt == 1 && countPVGTmupionsAlt == 1) : ((countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)));
    if (isTwoSelectedTracks && doTwoTracks) {
      TLorentzVector daug[2], pion[2], muon[2];
      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(cutTauEvent.useThresholdsPID ? vecPVnewPIDidx[0] : vecPVidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(cutTauEvent.useThresholdsPID ? vecPVnewPIDidx[1] : vecPVidx[1]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      pion[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(MassPionCharged, trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      pion[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(MassPionCharged, trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      muon[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(MassMuon, trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      muon[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(MassMuon, trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      if (cutTauEvent.applyTauEventSelection && !selectedTauEvent(trkDaug1, trkDaug2)) {
        return;
      }
      if (cutTauEvent.useThresholdsPID) {
        if (isElectronCandidate(trkDaug1))
          daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(MassElectron, trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
        if (isElectronCandidate(trkDaug2))
          daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(MassElectron, trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
        if (isMuPionCandidate(trkDaug1))
          daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(MassPionCharged, trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
        if (isMuPionCandidate(trkDaug2))
          daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(MassPionCharged, trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      }

      if (trkDaug1.hasTPC()) {
        histos.get<TH2>(HIST("EventTwoTracks/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        histos.get<TH2>(HIST("EventTwoTracks/PID/hTPCnSigmaElVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaEl());
        histos.get<TH2>(HIST("EventTwoTracks/PID/hTPCnSigmaMuVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaMu());
        histos.get<TH2>(HIST("EventTwoTracks/PID/hTPCnSigmaPiVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaPi());
        histos.get<TH2>(HIST("EventTwoTracks/PID/hTPCnSigmaKaVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaKa());
        histos.get<TH2>(HIST("EventTwoTracks/PID/hTPCnSigmaPrVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaPr());
        if (trkDaug1.hasTOF()) {
          histos.get<TH2>(HIST("EventTwoTracks/PID/hTOFsignalVsP"))->Fill(daug[0].P(), trkDaug1.tofSignal());
        }
        if (countPVGTelectrons == 2) {
          histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
          histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaEl());
          if (trkDaug1.hasTOF()) {
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsP"))->Fill(daug[0].P(), trkDaug1.tofSignal());
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsP"))->Fill(daug[0].P(), trkDaug1.tofNSigmaEl());
          }
        }
        if (countPVGTmuons == 2)
          histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTpions == 2)
          histos.get<TH2>(HIST("EventTwoTracks/TwoPions/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTmuons == 1)
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 1)
          histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTpions == 1 && countPVGTmuons == 1)
          histos.get<TH2>(HIST("EventTwoTracks/MuonPion/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (isElMuPion) {
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaEl());
          if (trkDaug1.hasTOF()) {
            histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsP"))->Fill(daug[0].P(), trkDaug1.tofSignal());
            histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsP"))->Fill(daug[0].P(), trkDaug1.tofNSigmaEl());
          }
        }
        if ((countPVGTelectrons == 2) || (countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)) {
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaEl());
          if (trkDaug1.hasTOF()) {
            histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFsignalVsP"))->Fill(daug[0].P(), trkDaug1.tofSignal());
            histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsP"))->Fill(daug[0].P(), trkDaug1.tofNSigmaEl());
          }
        }
      }
      if (trkDaug2.hasTPC()) {
        histos.get<TH2>(HIST("EventTwoTracks/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        histos.get<TH2>(HIST("EventTwoTracks/PID/hTPCnSigmaElVsP"))->Fill(daug[1].P(), trkDaug1.tpcNSigmaEl());
        histos.get<TH2>(HIST("EventTwoTracks/PID/hTPCnSigmaMuVsP"))->Fill(daug[1].P(), trkDaug1.tpcNSigmaMu());
        histos.get<TH2>(HIST("EventTwoTracks/PID/hTPCnSigmaPiVsP"))->Fill(daug[1].P(), trkDaug1.tpcNSigmaPi());
        histos.get<TH2>(HIST("EventTwoTracks/PID/hTPCnSigmaKaVsP"))->Fill(daug[1].P(), trkDaug1.tpcNSigmaKa());
        histos.get<TH2>(HIST("EventTwoTracks/PID/hTPCnSigmaPrVsP"))->Fill(daug[1].P(), trkDaug1.tpcNSigmaPr());
        if (trkDaug2.hasTOF()) {
          histos.get<TH2>(HIST("EventTwoTracks/PID/hTOFsignalVsP"))->Fill(daug[1].P(), trkDaug2.tofSignal());
        }
        if (countPVGTelectrons == 2) {
          histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
          histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaEl());
          if (trkDaug2.hasTOF()) {
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsP"))->Fill(daug[1].P(), trkDaug2.tofSignal());
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsP"))->Fill(daug[1].P(), trkDaug2.tofNSigmaEl());
          }
        }
        if (countPVGTmuons == 2)
          histos.get<TH2>(HIST("EventTwoTracks/TwoMuons/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTpions == 2)
          histos.get<TH2>(HIST("EventTwoTracks/TwoPions/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTmuons == 1)
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuon/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 1)
          histos.get<TH2>(HIST("EventTwoTracks/ElectronPion/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTpions == 1 && countPVGTmuons == 1)
          histos.get<TH2>(HIST("EventTwoTracks/MuonPion/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (isElMuPion) {
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaEl());
          if (trkDaug2.hasTOF()) {
            histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsP"))->Fill(daug[1].P(), trkDaug2.tofSignal());
            histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsP"))->Fill(daug[1].P(), trkDaug2.tofNSigmaEl());
          }
        }
        if ((countPVGTelectrons == 2) || (countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)) {
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaEl());
          if (trkDaug2.hasTOF()) {
            histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFsignalVsP"))->Fill(daug[1].P(), trkDaug2.tofSignal());
            histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsP"))->Fill(daug[1].P(), trkDaug2.tofNSigmaEl());
          }
        }
      }
      if (trkDaug1.hasTPC() && trkDaug2.hasTPC()) {
        if (countPVGTelectrons == 2) {
          if (daug[0].P() > daug[1].P()) {
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsLP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsOP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsLP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaEl());
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsOP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaEl());
            if (trkDaug1.hasTOF()) {
              histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsLP"))->Fill(daug[0].P(), trkDaug1.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsLP"))->Fill(daug[0].P(), trkDaug1.tofNSigmaEl());
            }
            if (trkDaug2.hasTOF()) {
              histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsOP"))->Fill(daug[1].P(), trkDaug2.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsOP"))->Fill(daug[1].P(), trkDaug2.tofNSigmaEl());
            }
          } else {
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsOP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCsignalVsLP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsOP"))->Fill(daug[0].P(), trkDaug1.tpcNSigmaEl());
            histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTPCnSigmaVsLP"))->Fill(daug[1].P(), trkDaug2.tpcNSigmaEl());
            if (trkDaug1.hasTOF()) {
              histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsOP"))->Fill(daug[0].P(), trkDaug1.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsOP"))->Fill(daug[0].P(), trkDaug1.tofNSigmaEl());
            }
            if (trkDaug2.hasTOF()) {
              histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFsignalVsLP"))->Fill(daug[1].P(), trkDaug2.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/TwoElectrons/PID/hTOFnSigmaVsLP"))->Fill(daug[1].P(), trkDaug2.tofNSigmaEl());
            }
          }
        }
        if (!isMC && isElMuPion) {
          double electronPt = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt();
          double electronPID = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcSignal() : trkDaug2.tpcSignal();
          double electronNsigmaEl = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcNSigmaEl() : trkDaug2.tpcNSigmaEl();
          double electronNsigmaPi = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcNSigmaPi() : trkDaug2.tpcNSigmaPi();
          double otherPt = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug2) : enumMyParticle(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt();
          double otherPID = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug2) : enumMyParticle(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcSignal() : trkDaug2.tpcSignal();
          double otherNsigmaMu = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug2) : enumMyParticle(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcNSigmaMu() : trkDaug2.tpcNSigmaMu();
          double otherNsigmaPi = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug2) : enumMyParticle(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcNSigmaPi() : trkDaug2.tpcNSigmaPi();
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCsignalVsEPofE"))->Fill(electronPt, electronPID);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCsignalVsOPofO"))->Fill(otherPt, otherPID);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsEPofE"))->Fill(electronPt, electronNsigmaEl);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsMPofO"))->Fill(otherPt, otherNsigmaMu);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsPPofO"))->Fill(otherPt, otherNsigmaPi);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaEvsnSigmaPofE"))->Fill(electronNsigmaEl, electronNsigmaPi);
          if (trkDaug1.hasTOF()) {
            if (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsEPofE"))->Fill(electronPt, trkDaug1.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsEPofE"))->Fill(electronPt, trkDaug1.tofNSigmaEl());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaEvsnSigmaPofE"))->Fill(trkDaug1.tofNSigmaEl(), trkDaug1.tofNSigmaPi());
            } else {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsOPofO"))->Fill(otherPt, trkDaug1.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsMPofO"))->Fill(otherPt, trkDaug1.tofNSigmaMu());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsPPofO"))->Fill(otherPt, trkDaug1.tofNSigmaPi());
            }
          }
          if (trkDaug2.hasTOF()) {
            if (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug2) : enumMyParticle(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsEPofE"))->Fill(electronPt, trkDaug2.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsEPofE"))->Fill(electronPt, trkDaug2.tofNSigmaEl());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaEvsnSigmaPofE"))->Fill(trkDaug2.tofNSigmaEl(), trkDaug2.tofNSigmaPi());
            } else {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsOPofO"))->Fill(otherPt, trkDaug2.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsMPofO"))->Fill(otherPt, trkDaug2.tofNSigmaMu());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsPPofO"))->Fill(otherPt, trkDaug2.tofNSigmaPi());
            }
          }
        }
        if ((countPVGTelectrons == 2) || (countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)) {
          double electronPt = (enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt();
          double electronPID = (enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcSignal() : trkDaug2.tpcSignal();
          double electronNsigmaEl = (enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcNSigmaEl() : trkDaug2.tpcNSigmaEl();
          double otherPt = (enumMyParticle(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt();
          double otherPID = (enumMyParticle(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcSignal() : trkDaug2.tpcSignal();
          double otherNsigmaMu = (enumMyParticle(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcNSigmaMu() : trkDaug2.tpcNSigmaMu();
          double otherNsigmaPi = (enumMyParticle(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcNSigmaPi() : trkDaug2.tpcNSigmaPi();
          if (countPVGTelectrons == 2) {
            electronPt = (daug[0].Pt() > daug[1].Pt()) ? daug[0].Pt() : daug[1].Pt();
            electronPID = (daug[0].Pt() > daug[1].Pt()) ? trkDaug1.tpcSignal() : trkDaug2.tpcSignal();
            electronNsigmaEl = (daug[0].Pt() > daug[1].Pt()) ? trkDaug1.tpcNSigmaEl() : trkDaug2.tpcNSigmaEl();
            otherPt = (daug[0].Pt() > daug[1].Pt()) ? daug[1].Pt() : daug[0].Pt();
            otherPID = (daug[0].Pt() > daug[1].Pt()) ? trkDaug2.tpcSignal() : trkDaug1.tpcSignal();
            otherNsigmaMu = (daug[0].Pt() > daug[1].Pt()) ? trkDaug2.tpcNSigmaMu() : trkDaug1.tpcNSigmaMu();
            otherNsigmaPi = (daug[0].Pt() > daug[1].Pt()) ? trkDaug2.tpcNSigmaPi() : trkDaug1.tpcNSigmaPi();
          }
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCsignalVsEP"))->Fill(electronPt, electronPID);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCsignalVsOP"))->Fill(otherPt, otherPID);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsEP"))->Fill(electronPt, electronNsigmaEl);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsMP"))->Fill(otherPt, otherNsigmaMu);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTPCnSigmaVsPP"))->Fill(otherPt, otherNsigmaPi);
          if (trkDaug1.hasTOF()) {
            if (enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFsignalVsEP"))->Fill(electronPt, trkDaug1.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsEP"))->Fill(electronPt, trkDaug1.tofNSigmaEl());
            } else {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFsignalVsOP"))->Fill(otherPt, trkDaug1.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsMP"))->Fill(otherPt, trkDaug1.tofNSigmaMu());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsPP"))->Fill(otherPt, trkDaug1.tofNSigmaPi());
            }
          }
          if (trkDaug2.hasTOF()) {
            if (enumMyParticle(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFsignalVsEP"))->Fill(electronPt, trkDaug2.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsEP"))->Fill(electronPt, trkDaug2.tofNSigmaEl());
            } else {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFsignalVsOP"))->Fill(otherPt, trkDaug2.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsMP"))->Fill(otherPt, trkDaug2.tofNSigmaMu());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronOther/PID/hTOFnSigmaVsPP"))->Fill(otherPt, trkDaug2.tofNSigmaPi());
            }
          }
        }
      }
    } else if (countPVGTselected == 4 && doFourTracks) {

      TLorentzVector daug[4];
      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(vecPVidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(vecPVidx[1]);
      const auto& trkDaug3 = reconstructedBarrelTracks.iteratorAt(vecPVidx[2]);
      const auto& trkDaug4 = reconstructedBarrelTracks.iteratorAt(vecPVidx[3]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      daug[2].SetPxPyPzE(trkDaug3.px(), trkDaug3.py(), trkDaug3.pz(), energy(pdg->Mass(trackPDG(trkDaug3, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug3.px(), trkDaug3.py(), trkDaug3.pz()));
      daug[3].SetPxPyPzE(trkDaug4.px(), trkDaug4.py(), trkDaug4.pz(), energy(pdg->Mass(trackPDG(trkDaug4, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug4.px(), trkDaug4.py(), trkDaug4.pz()));

      if (trkDaug1.hasTPC()) {
        histos.get<TH2>(HIST("EventFourTracks/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTpions == 4)
          histos.get<TH2>(HIST("EventFourTracks/WithPion/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 3)
          histos.get<TH2>(HIST("EventFourTracks/WithElectron/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTpions == 3 && countPVGTmuons == 1)
          histos.get<TH2>(HIST("EventFourTracks/WithMuon/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
      }
      if (trkDaug2.hasTPC()) {
        histos.get<TH2>(HIST("EventFourTracks/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTpions == 4)
          histos.get<TH2>(HIST("EventFourTracks/WithPion/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 3)
          histos.get<TH2>(HIST("EventFourTracks/WithElectron/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTpions == 3 && countPVGTmuons == 1)
          histos.get<TH2>(HIST("EventFourTracks/WithMuon/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
      }
      if (trkDaug3.hasTPC()) {
        histos.get<TH2>(HIST("EventFourTracks/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
        if (countPVGTpions == 4)
          histos.get<TH2>(HIST("EventFourTracks/WithPion/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 3)
          histos.get<TH2>(HIST("EventFourTracks/WithElectron/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
        if (countPVGTpions == 3 && countPVGTmuons == 1)
          histos.get<TH2>(HIST("EventFourTracks/WithMuon/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
      }
      if (trkDaug4.hasTPC()) {
        histos.get<TH2>(HIST("EventFourTracks/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
        if (countPVGTpions == 4)
          histos.get<TH2>(HIST("EventFourTracks/WithPion/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
        if (countPVGTelectrons == 1 && countPVGTpions == 3)
          histos.get<TH2>(HIST("EventFourTracks/WithElectron/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
        if (countPVGTpions == 3 && countPVGTmuons == 1)
          histos.get<TH2>(HIST("EventFourTracks/WithMuon/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
      }
    } else if (countPVGTselected == 6 && doSixTracks) {
      TLorentzVector daug[6];
      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(vecPVidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(vecPVidx[1]);
      const auto& trkDaug3 = reconstructedBarrelTracks.iteratorAt(vecPVidx[2]);
      const auto& trkDaug4 = reconstructedBarrelTracks.iteratorAt(vecPVidx[3]);
      const auto& trkDaug5 = reconstructedBarrelTracks.iteratorAt(vecPVidx[4]);
      const auto& trkDaug6 = reconstructedBarrelTracks.iteratorAt(vecPVidx[5]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      daug[2].SetPxPyPzE(trkDaug3.px(), trkDaug3.py(), trkDaug3.pz(), energy(pdg->Mass(trackPDG(trkDaug3, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug3.px(), trkDaug3.py(), trkDaug3.pz()));
      daug[3].SetPxPyPzE(trkDaug4.px(), trkDaug4.py(), trkDaug4.pz(), energy(pdg->Mass(trackPDG(trkDaug4, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug4.px(), trkDaug4.py(), trkDaug4.pz()));
      daug[4].SetPxPyPzE(trkDaug5.px(), trkDaug5.py(), trkDaug5.pz(), energy(pdg->Mass(trackPDG(trkDaug5, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug5.px(), trkDaug5.py(), trkDaug5.pz()));
      daug[5].SetPxPyPzE(trkDaug6.px(), trkDaug6.py(), trkDaug6.pz(), energy(pdg->Mass(trackPDG(trkDaug6, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug6.px(), trkDaug6.py(), trkDaug6.pz()));

      if (trkDaug1.hasTPC()) {
        histos.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
        if (countPVGTpions == 6)
          histos.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[0].P(), trkDaug1.tpcSignal());
      }
      if (trkDaug2.hasTPC()) {
        histos.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
        if (countPVGTpions == 6)
          histos.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[1].P(), trkDaug2.tpcSignal());
      }
      if (trkDaug3.hasTPC()) {
        histos.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
        if (countPVGTpions == 6)
          histos.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[2].P(), trkDaug3.tpcSignal());
      }
      if (trkDaug4.hasTPC()) {
        histos.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
        if (countPVGTpions == 6)
          histos.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[3].P(), trkDaug4.tpcSignal());
      }
      if (trkDaug5.hasTPC()) {
        histos.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[4].P(), trkDaug5.tpcSignal());
        if (countPVGTpions == 6)
          histos.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[4].P(), trkDaug5.tpcSignal());
      }
      if (trkDaug6.hasTPC()) {
        histos.get<TH2>(HIST("EventSixTracks/PID/hTPCsignalVsP"))->Fill(daug[5].P(), trkDaug6.tpcSignal());
        if (countPVGTpions == 6)
          histos.get<TH2>(HIST("EventSixTracks/SixPions/PID/hTPCsignalVsP"))->Fill(daug[5].P(), trkDaug6.tpcSignal());
      }
    } else {
      printDebugMessage("Other particles");
    }

  } // end fillPIDhistograms

  void fillMCPIDhistograms(FullMCUDTracks const& reconstructedBarrelTracks)
  {

    int countPVGTselected = 0;
    int countPVGTelectrons = 0;
    int countPVGTmuons = 0;
    int countPVGTpions = 0;
    int countPVGTelmupiAlt = 0;
    int countPVGTelectronsAlt = 0;
    int countPVGTmupionsAlt = 0;
    std::vector<int> vecPVidx;
    std::vector<int> vecPVnewPIDidx;
    // Loop over tracks with selections
    for (const auto& track : reconstructedBarrelTracks) {
      if (!track.isPVContributor())
        continue;
      if (cutGlobalTrack.applyGlobalTrackSelection) {
        if (isGlobalTrackReinstatement(track) != 1)
          continue;
      }
      int hypothesisID = testPIDhypothesis(track, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC);
      if (hypothesisID == P_ELECTRON || hypothesisID == P_MUON || hypothesisID == P_PION) {
        countPVGTselected++;
        vecPVidx.push_back(track.index());
        if (hypothesisID == P_ELECTRON) {
          countPVGTelectrons++;
        } else if (hypothesisID == P_MUON) {
          countPVGTmuons++;
        } else {
          countPVGTpions++;
        }
      }
      // alternative selection
      if (isElectronCandidate(track)) {
        countPVGTelmupiAlt++;
        countPVGTelectronsAlt++;
        vecPVnewPIDidx.push_back(track.index());
      }
      if (isMuPionCandidate(track)) {
        countPVGTelmupiAlt++;
        countPVGTmupionsAlt++;
        vecPVnewPIDidx.push_back(track.index());
      }

    } // Loop over tracks with selections

    bool isTwoSelectedTracks = (cutTauEvent.useThresholdsPID ? countPVGTelmupiAlt == 2 : countPVGTselected == 2);
    bool isElMuPion = (cutTauEvent.useThresholdsPID ? (countPVGTelectronsAlt == 1 && countPVGTmupionsAlt == 1) : ((countPVGTelectrons == 1 && countPVGTmuons == 1) || (countPVGTelectrons == 1 && countPVGTpions == 1)));
    if (isTwoSelectedTracks && doTwoTracks) {
      TLorentzVector daug[2];
      const auto& trkDaug1 = reconstructedBarrelTracks.iteratorAt(cutTauEvent.useThresholdsPID ? vecPVnewPIDidx[0] : vecPVidx[0]);
      const auto& trkDaug2 = reconstructedBarrelTracks.iteratorAt(cutTauEvent.useThresholdsPID ? vecPVnewPIDidx[1] : vecPVidx[1]);
      daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(pdg->Mass(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
      daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(pdg->Mass(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)), trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      if (cutTauEvent.applyTauEventSelection && !selectedTauEvent(trkDaug1, trkDaug2)) {
        return;
      }
      if (cutTauEvent.useThresholdsPID) {
        if (isElectronCandidate(trkDaug1))
          daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(MassElectron, trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
        if (isElectronCandidate(trkDaug2))
          daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(MassElectron, trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
        if (isMuPionCandidate(trkDaug1))
          daug[0].SetPxPyPzE(trkDaug1.px(), trkDaug1.py(), trkDaug1.pz(), energy(MassPionCharged, trkDaug1.px(), trkDaug1.py(), trkDaug1.pz()));
        if (isMuPionCandidate(trkDaug2))
          daug[1].SetPxPyPzE(trkDaug2.px(), trkDaug2.py(), trkDaug2.pz(), energy(MassPionCharged, trkDaug2.px(), trkDaug2.py(), trkDaug2.pz()));
      }
      if (trkDaug1.hasTPC() && trkDaug2.hasTPC()) {

        if (isMC && isElMuPion) {
          int pid = 0;
          if (trkDaug1.has_udMcParticle()) {
            const auto& part = trkDaug1.udMcParticle();
            pid = std::abs(part.pdgCode());
            if (pid == 11) {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/mcTruth/nSigmaTPC1"))->Fill(trkDaug1.pt(), trkDaug1.tpcNSigmaEl(), 1.);
            } else {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/mcTruth/nSigmaTPC2"))->Fill(trkDaug1.pt(), trkDaug1.tpcNSigmaEl(), 1.);
            }
          }
          if (trkDaug2.has_udMcParticle()) {
            const auto& part = trkDaug2.udMcParticle();
            pid = std::abs(part.pdgCode());
            if (pid == 11) {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/mcTruth/nSigmaTPC1"))->Fill(trkDaug2.pt(), trkDaug2.tpcNSigmaEl(), 1.);
            } else {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/mcTruth/nSigmaTPC2"))->Fill(trkDaug2.pt(), trkDaug2.tpcNSigmaEl(), 1.);
            }
          }
          bool isNotTrueElectron = false;
          if (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) {
            if (trkDaug1.has_udMcParticle()) {
              const auto& particle = trkDaug1.udMcParticle();
              if (enumMyParticle(particle.pdgCode()) != P_ELECTRON)
                isNotTrueElectron = true;
            }
          } else {
            if (trkDaug2.has_udMcParticle()) {
              const auto& particle = trkDaug2.udMcParticle();
              if (enumMyParticle(particle.pdgCode()) != P_ELECTRON)
                isNotTrueElectron = true;
            }
          }
          if (oppositeMCtrueElectronCheck) {
            if (doMCtrueElectronCheck && !isNotTrueElectron)
              return;
          } else {
            if (doMCtrueElectronCheck && isNotTrueElectron)
              return;
          }

          double electronPt = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt();
          double electronPID = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcSignal() : trkDaug2.tpcSignal();
          double electronNsigmaEl = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcNSigmaEl() : trkDaug2.tpcNSigmaEl();
          double electronNsigmaPi = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcNSigmaPi() : trkDaug2.tpcNSigmaPi();
          double otherPt = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug2) : enumMyParticle(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? daug[0].Pt() : daug[1].Pt();
          double otherPID = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug2) : enumMyParticle(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcSignal() : trkDaug2.tpcSignal();
          double otherNsigmaEl = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug2) : enumMyParticle(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcNSigmaEl() : trkDaug2.tpcNSigmaEl();
          double otherNsigmaMu = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug2) : enumMyParticle(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcNSigmaMu() : trkDaug2.tpcNSigmaMu();
          double otherNsigmaPi = (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug2) : enumMyParticle(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) ? trkDaug1.tpcNSigmaPi() : trkDaug2.tpcNSigmaPi();
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCsignalVsEPofE"))->Fill(electronPt, electronPID);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsEPofE"))->Fill(electronPt, electronNsigmaEl);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsPPofE"))->Fill(electronPt, electronNsigmaPi);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaEvsnSigmaPofE"))->Fill(electronNsigmaEl, electronNsigmaPi);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCsignalVsOPofO"))->Fill(otherPt, otherPID);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsEPofO"))->Fill(otherPt, otherNsigmaEl);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsMPofO"))->Fill(otherPt, otherNsigmaMu);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaVsPPofO"))->Fill(otherPt, otherNsigmaPi);
          histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTPCnSigmaEvsnSigmaPofO"))->Fill(otherNsigmaEl, otherNsigmaPi);
          if (trkDaug1.hasTOF()) {
            if (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug1) : enumMyParticle(trackPDG(trkDaug1, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsEPofE"))->Fill(electronPt, trkDaug1.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsEPofE"))->Fill(electronPt, trkDaug1.tofNSigmaEl());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsPPofE"))->Fill(electronPt, trkDaug1.tofNSigmaPi());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaEvsnSigmaPofE"))->Fill(trkDaug1.tofNSigmaEl(), trkDaug1.tofNSigmaPi());
            } else {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsOPofO"))->Fill(otherPt, trkDaug1.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsEPofO"))->Fill(otherPt, trkDaug1.tofNSigmaEl());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsMPofO"))->Fill(otherPt, trkDaug1.tofNSigmaMu());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsPPofO"))->Fill(otherPt, trkDaug1.tofNSigmaPi());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaEvsnSigmaPofO"))->Fill(trkDaug1.tofNSigmaEl(), trkDaug1.tofNSigmaPi());
            }
          }
          if (trkDaug2.hasTOF()) {
            if (cutTauEvent.useThresholdsPID ? isElectronCandidate(trkDaug2) : enumMyParticle(trackPDG(trkDaug2, cutPID.cutSiTPC, cutPID.cutSiTOF, cutPID.usePIDwTOF, cutPID.useScutTOFinTPC)) == P_ELECTRON) {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsEPofE"))->Fill(electronPt, trkDaug2.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsEPofE"))->Fill(electronPt, trkDaug2.tofNSigmaEl());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsPPofE"))->Fill(electronPt, trkDaug2.tofNSigmaPi());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaEvsnSigmaPofE"))->Fill(trkDaug2.tofNSigmaEl(), trkDaug2.tofNSigmaPi());
            } else {
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFsignalVsOPofO"))->Fill(otherPt, trkDaug2.tofSignal());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsEPofO"))->Fill(otherPt, trkDaug2.tofNSigmaEl());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsMPofO"))->Fill(otherPt, trkDaug2.tofNSigmaMu());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaVsPPofO"))->Fill(otherPt, trkDaug2.tofNSigmaPi());
              histos.get<TH2>(HIST("EventTwoTracks/ElectronMuPi/PID/hTOFnSigmaEvsnSigmaPofO"))->Fill(trkDaug2.tofNSigmaEl(), trkDaug2.tofNSigmaPi());
            }
          }
        }
      }
    }

  } // end fillMCPIDhistograms

  template <typename C>
  void fillFIThistograms(C const& reconstructedCollision)
  {

    histos.get<TH1>(HIST("Events/FIT/hAmplitudeFT0A"))->Fill(reconstructedCollision.totalFT0AmplitudeA());
    histos.get<TH1>(HIST("Events/FIT/hAmplitudeFT0C"))->Fill(reconstructedCollision.totalFT0AmplitudeC());
    histos.get<TH1>(HIST("Events/FIT/hAmplitudeFDDA"))->Fill(reconstructedCollision.totalFDDAmplitudeA());
    histos.get<TH1>(HIST("Events/FIT/hAmplitudeFDDC"))->Fill(reconstructedCollision.totalFDDAmplitudeC());
    histos.get<TH1>(HIST("Events/FIT/hAmplitudeFV0A"))->Fill(reconstructedCollision.totalFV0AmplitudeA());

    histos.get<TH1>(HIST("Events/FIT/hTimeFT0A"))->Fill(reconstructedCollision.timeFT0A());
    histos.get<TH1>(HIST("Events/FIT/hTimeFT0C"))->Fill(reconstructedCollision.timeFT0C());
    histos.get<TH1>(HIST("Events/FIT/hTimeFDDA"))->Fill(reconstructedCollision.timeFDDA());
    histos.get<TH1>(HIST("Events/FIT/hTimeFDDC"))->Fill(reconstructedCollision.timeFDDC());
    histos.get<TH1>(HIST("Events/FIT/hTimeFV0A"))->Fill(reconstructedCollision.timeFV0A());

    histos.get<TH2>(HIST("Events/FIT/hTimeFT0AvsFT0C"))->Fill(reconstructedCollision.timeFT0A(), reconstructedCollision.timeFT0C());
    histos.get<TH2>(HIST("Events/FIT/hTimeFT0CvsFDDA"))->Fill(reconstructedCollision.timeFT0C(), reconstructedCollision.timeFDDA());
    histos.get<TH2>(HIST("Events/FIT/hTimeFDDAvsFDDC"))->Fill(reconstructedCollision.timeFDDA(), reconstructedCollision.timeFDDC());
    histos.get<TH2>(HIST("Events/FIT/hTimeFDDCvsFV0A"))->Fill(reconstructedCollision.timeFDDC(), reconstructedCollision.timeFV0A());
    histos.get<TH2>(HIST("Events/FIT/hTimeFV0AvsFT0A"))->Fill(reconstructedCollision.timeFV0A(), reconstructedCollision.timeFT0A());
  }

  void fillTruthHistograms(aod::UDMcParticles const& particles)
  {
    histos.get<TH1>(HIST("Events/Truth/hCountCollisions"))->Fill(1);
    histos.get<TH1>(HIST("Events/Truth/hNparticles"))->Fill(particles.size());
    histos.get<TH2>(HIST("Events/Truth/hNphysPartVsNwoutMotherParts"))->Fill(countPhysicalPrimary(particles), countParticlesWithoutMother(particles));

    int countElectrons = 0;
    int countMuons = 0;
    int countPions = 0;

    for (const auto& particle : particles) {
      histos.get<TH1>(HIST("Events/Truth/hPDGcodesAll"))->Fill(particle.pdgCode());
      //        if (!particle.isPhysicalPrimary()) continue;
      if (particle.has_mothers())
        continue;
      histos.get<TH1>(HIST("Events/Truth/hPDGcodesNoMother"))->Fill(particle.pdgCode());
      histos.get<TH1>(HIST("Tracks/Truth/hTauPt"))->Fill(pt(particle.px(), particle.py()));
      histos.get<TH1>(HIST("Tracks/Truth/hTauP"))->Fill(momentum(particle.px(), particle.py(), particle.pz()));
      histos.get<TH1>(HIST("Tracks/Truth/hTauPhi"))->Fill(phi(particle.px(), particle.py()));
      histos.get<TH1>(HIST("Tracks/Truth/hTauEta"))->Fill(eta(particle.px(), particle.py(), particle.pz()));
      const auto& daughters = particle.daughters_as<aod::UDMcParticles>();
      histos.get<TH1>(HIST("Events/Truth/hNtauDaughters"))->Fill(daughters.size());
      for (const auto& daughter : daughters) {
        histos.get<TH1>(HIST("Events/Truth/hPDGcodesTauDaughters"))->Fill(daughter.pdgCode());
        if (enumMyParticle(daughter.pdgCode()) == P_ELECTRON) {
          countElectrons++;
          histos.get<TH1>(HIST("Tracks/Truth/hElectronPt"))->Fill(pt(daughter.px(), daughter.py()));
          histos.get<TH1>(HIST("Tracks/Truth/hElectronP"))->Fill(momentum(daughter.px(), daughter.py(), daughter.pz()));
          histos.get<TH1>(HIST("Tracks/Truth/hElectronPhi"))->Fill(phi(daughter.px(), daughter.py()));
          histos.get<TH1>(HIST("Tracks/Truth/hElectronEta"))->Fill(eta(daughter.px(), daughter.py(), daughter.pz()));
        }
        if (enumMyParticle(daughter.pdgCode()) == P_MUON) {
          countMuons++;
          histos.get<TH1>(HIST("Tracks/Truth/hMuonPt"))->Fill(pt(daughter.px(), daughter.py()));
          histos.get<TH1>(HIST("Tracks/Truth/hMuonP"))->Fill(momentum(daughter.px(), daughter.py(), daughter.pz()));
          histos.get<TH1>(HIST("Tracks/Truth/hMuonPhi"))->Fill(phi(daughter.px(), daughter.py()));
          histos.get<TH1>(HIST("Tracks/Truth/hMuonEta"))->Fill(eta(daughter.px(), daughter.py(), daughter.pz()));
        }
        if (enumMyParticle(daughter.pdgCode()) == P_PION) {
          countPions++;
          histos.get<TH1>(HIST("Tracks/Truth/hPionPt"))->Fill(pt(daughter.px(), daughter.py()));
          histos.get<TH1>(HIST("Tracks/Truth/hPionP"))->Fill(momentum(daughter.px(), daughter.py(), daughter.pz()));
          histos.get<TH1>(HIST("Tracks/Truth/hPionPhi"))->Fill(phi(daughter.px(), daughter.py()));
          histos.get<TH1>(HIST("Tracks/Truth/hPionEta"))->Fill(eta(daughter.px(), daughter.py(), daughter.pz()));
        }
      }
    }

    histos.get<TH1>(HIST("Events/Truth/hNelectrons"))->Fill(countElectrons);
    histos.get<TH1>(HIST("Events/Truth/hNmuons"))->Fill(countMuons);
    histos.get<TH1>(HIST("Events/Truth/hNpions"))->Fill(countPions);

    if (countElectrons == 2 && countMuons == 0 && countPions == 0)
      histos.get<TH1>(HIST("Events/Truth/hChannels"))->Fill(CH_EE);
    if (countElectrons == 1 && countMuons == 1 && countPions == 0)
      histos.get<TH1>(HIST("Events/Truth/hChannels"))->Fill(CH_EMU);
    if (countElectrons == 1 && countMuons == 0 && countPions == 1)
      histos.get<TH1>(HIST("Events/Truth/hChannels"))->Fill(CH_EPI);
    if ((countElectrons == 1 && countMuons == 1 && countPions == 0) || (countElectrons == 1 && countMuons == 0 && countPions == 1))
      histos.get<TH1>(HIST("Events/Truth/hChannels"))->Fill(CH_EMUPI);
    if (countElectrons == 0 && countMuons == 2 && countPions == 0)
      histos.get<TH1>(HIST("Events/Truth/hChannels"))->Fill(CH_MUMU);
    if (countElectrons == 0 && countMuons == 1 && countPions == 1)
      histos.get<TH1>(HIST("Events/Truth/hChannels"))->Fill(CH_MUPI);
    if (countElectrons == 0 && countMuons == 0 && countPions == 2)
      histos.get<TH1>(HIST("Events/Truth/hChannels"))->Fill(CH_PIPI);
    if (countElectrons == 0 && countMuons == 0 && countPions == 4)
      histos.get<TH1>(HIST("Events/Truth/hChannels"))->Fill(CH_FOURPI);
    if (countElectrons == 1 && countMuons == 0 && countPions == 3)
      histos.get<TH1>(HIST("Events/Truth/hChannels"))->Fill(CH_ETHREEPI);
    if (countElectrons == 0 && countMuons == 1 && countPions == 3)
      histos.get<TH1>(HIST("Events/Truth/hChannels"))->Fill(CH_MUTHREEPI);
    if (countElectrons == 0 && countMuons == 0 && countPions == 6)
      histos.get<TH1>(HIST("Events/Truth/hChannels"))->Fill(CH_SIXPI);
  }

  void processDataDG(FullUDCollision const& reconstructedCollision,
                     FullUDTracks const& reconstructedBarrelTracks)
  {

    if (!isGoodFITtime(reconstructedCollision, cutSample.cutFITtime))
      return;

    if (doMainHistos) {
      fillHistograms(reconstructedBarrelTracks);
      fillFIThistograms(reconstructedCollision);
    }
    if (doPIDhistos)
      fillPIDhistograms(reconstructedCollision, reconstructedBarrelTracks);

  } // end processDataDG

  void processDataSG(FullSGUDCollision const& reconstructedCollision,
                     FullUDTracks const& reconstructedBarrelTracks)
  {

    int gapSide = reconstructedCollision.gapSide();
    int trueGapSide = sgSelector.trueGap(reconstructedCollision, cutSample.cutTrueGapSideFV0, cutSample.cutTrueGapSideFT0A, cutSample.cutTrueGapSideFT0C, cutSample.cutTrueGapSideZDC);
    histos.fill(HIST("Events/UDtableGapSide"), gapSide);
    histos.fill(HIST("Events/TrueGapSideDiffToTableValue"), gapSide - trueGapSide);
    if (cutSample.useTrueGap)
      gapSide = trueGapSide;

    if (gapSide != cutSample.whichGapSide)
      return;

    if (!isGoodFITtime(reconstructedCollision, cutSample.cutFITtime))
      return;

    if (doMainHistos) {
      fillHistograms(reconstructedBarrelTracks);
      fillFIThistograms(reconstructedCollision);
    }
    if (doPIDhistos)
      fillPIDhistograms(reconstructedCollision, reconstructedBarrelTracks);

  } // end processDataSG

  void processMCrecDG(FullMCUDCollision const& reconstructedCollision,
                      FullMCUDTracks const& reconstructedBarrelTracks,
                      aod::UDMcParticles const&)
  {
    isMC = true;

    if (!isGoodFITtime(reconstructedCollision, cutSample.cutFITtime))
      return;

    if (cutSample.applyAcceptanceSelection) {
      for (const auto& track : reconstructedBarrelTracks) {
        if (!track.isPVContributor())
          continue;
        if (std::abs(eta(track.px(), track.py(), track.py())) > cutSample.cutTrackEta)
          return;
      }
    }

    if (doMainHistos) {
      fillHistograms(reconstructedBarrelTracks);
      fillFIThistograms(reconstructedCollision);
    }

    if (doPIDhistos) {
      fillPIDhistograms(reconstructedCollision, reconstructedBarrelTracks);
      fillMCPIDhistograms(reconstructedBarrelTracks);
    }

  } // end processMCrecDG

  void processMCrecSG(FullMCSGUDCollision const& reconstructedCollision,
                      FullMCUDTracks const& reconstructedBarrelTracks,
                      aod::UDMcParticles const&)
  {
    isMC = true;

    int gapSide = reconstructedCollision.gapSide();
    histos.fill(HIST("Events/UDtableGapSide"), gapSide);

    if (gapSide != cutSample.whichGapSide)
      return;

    if (!isGoodFITtime(reconstructedCollision, cutSample.cutFITtime))
      return;

    if (cutSample.applyAcceptanceSelection) {
      for (const auto& track : reconstructedBarrelTracks) {
        if (!track.isPVContributor())
          continue;
        if (std::abs(eta(track.px(), track.py(), track.py())) > cutSample.cutTrackEta)
          return;
      }
    }

    if (doMainHistos) {
      fillHistograms(reconstructedBarrelTracks);
      fillFIThistograms(reconstructedCollision);
    }

    if (doPIDhistos) {
      fillPIDhistograms(reconstructedCollision, reconstructedBarrelTracks);
      fillMCPIDhistograms(reconstructedBarrelTracks);
    }

  } // end processMCrecDG

  void processMCgen(aod::UDMcCollision const& /*generatedCollision*/,
                    aod::UDMcParticles const& particles)
  {
    isMC = true;

    if (cutSample.applyAcceptanceSelection) {
      for (const auto& particle : particles) {
        if (particle.has_mothers())
          continue;
        //        printLargeMessage(Form("GENE: eta %.3f cut %.2f",std::abs(eta(particle.px(), particle.py(), particle.py())),static_cast<double>(cutTrackEta)));
        if (std::abs(eta(particle.px(), particle.py(), particle.py())) > cutSample.cutTrackEta)
          return;
      }
    }

    if (doTruthHistos) {
      fillTruthHistograms(particles);
    }

  } // end processMCgenDG

  void processTestMC(FullMCUDCollision const& /*reconstructedCollision*/,
                     FullMCUDTracks const& /*reconstructedBarrelTracks*/,
                     aod::UDMcCollisions const&,
                     aod::UDMcParticles const&)
  {
    // if (reconstructedCollision.has_udMcCollision()) {
    //   const auto& generatedCollision = reconstructedCollision.udMcCollision();
    //   printDebugMessage(Form("%lli udMcCollision found", generatedCollision.size())); // FIXME: Type of size() is not invariant.
    // }

    // const auto& track = reconstructedBarrelTracks.iteratorAt(0);
    // if (track.size() && track.has_udMcParticle()) {
    //   const auto& particle = track.udMcParticle();
    //   printDebugMessage(Form("%lli udMcParticle found", particle.size())); // FIXME: Type of size() is not invariant.
    // }

  } // end processTestMC

  PROCESS_SWITCH(UpcTauRl, processDataDG, "Iterate UD tables with measured data created by DG-Candidate-Producer.", false);
  PROCESS_SWITCH(UpcTauRl, processDataSG, "Iterate UD tables with measured data created by SG-Candidate-Producer.", false);
  PROCESS_SWITCH(UpcTauRl, processMCrecDG, "Iterate Monte Carlo UD tables with reconstructed data created by DG-Candidate-Producer. Similar to processDataDG but uses association to truth level.", false);
  PROCESS_SWITCH(UpcTauRl, processMCrecSG, "Iterate Monte Carlo UD tables with reconstructed data created by SG-Candidate-Producer. Similar to processDataSG but uses association to truth level and trueGap is not available.", false);
  PROCESS_SWITCH(UpcTauRl, processMCgen, "Iterate Monte Carlo UD tables with truth data.", false);
  PROCESS_SWITCH(UpcTauRl, processTestMC, "Simple test of indices in MC sample.", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcTauRl>(cfgc)};
}
