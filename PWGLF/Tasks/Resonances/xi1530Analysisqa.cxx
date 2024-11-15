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

/// \file xi1530analysisqa.cxx
/// \brief Reconstruction of track-track decay resonance candidates
///
/// \author Min-jae Kim <minjae.kim@cern.ch>, Bong-Hwi Lim <bong-hwi.lim@cern.ch>
#include <TLorentzVector.h>
#include "TF1.h"
#include "TRandom3.h"

#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Framework/AnalysisTask.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/RecoDecay.h"
#include "Framework/O2DatabasePDGPlugin.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;
Service<o2::framework::O2DatabasePDG> pdgDB;

enum {
  kData = 0,
  kLS,
  kMixing,
  kMCReco,
  kMCTrue,
  kMCTruePS,
  kINEL10,
  kINELg010,
  kAllType
};

struct xi1530analysisqa {

  // Basic set-up //
  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  using ResoMCCols = soa::Join<aod::ResoCollisions, aod::ResoMCCollisions>;

  Configurable<float> cMassXiminus{"cMassXiminus", 1.31486, "Mass of Xi baryon"};

  // associated with histograms
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
  ConfigurableAxis binsPtQA{"binsPtQA", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0}, "Binning of the pT axis"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0}, "Binning of the centrality axis"};

  Configurable<float> cInvMassStart{"cInvMassStart", 1.2, "Invariant mass start"};
  Configurable<float> cInvMassEnd{"cInvMassEnd", 2.1, "Invariant mass end"};
  Configurable<int> cInvMassBins{"cInvMassBins", 900, "Invariant mass binning"};

  Configurable<int> cPIDBins{"cPIDBins", 65, "PID binning"};
  Configurable<float> cPIDQALimit{"cPIDQALimit", 6.5, "PID QA limit"};
  Configurable<int> cDCABins{"cDCABins", 150, "DCA binning"};

  Configurable<bool> invmass1D{"invmass1D", true, "Invariant mass 1D"};
  Configurable<bool> study_antiparticle{"study_antiparticle", true, "Study anti-particles separately"};
  Configurable<bool> PIDplots{"PIDplots", true, "Make TPC and TOF PID plots"};

  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of events to mix"};
  ConfigurableAxis CfgVtxBins{"CfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis CfgMultBins{"CfgMultBins", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Mixing bins - z-vertex"};

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  // Track selections (Execpt DCA selelctions) //

  // Primary track selections
  Configurable<double> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};
  Configurable<double> cMaxetacut{"cMaxetacut", 0.8, "Track maximum eta cut"};

  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", true, "Global track selection without DCA"};
  Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", true, "Global track selection"};
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};
  Configurable<bool> additionalQAplots{"additionalQAplots", true, "Additional QA plots"};
  Configurable<bool> additionalQAeventPlots{"additionalQAeventPlots", true, "Additional QA event plots"};
  Configurable<bool> additionalMEPlots{"additionalMEPlots", true, "Additional Mixed event plots"};

  Configurable<bool> tof_at_high_pt{"tof_at_high_pt", true, "Use TOF at high pT"};

  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"}; // Minmimum
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 0, "Number of TPC cluster"}; // Minmimum

  Configurable<float> cfgRatioTPCRowsOverFindableCls{"cfgRatioTPCRowsOverFindableCls", 0.0f, "TPC Crossed Rows to Findable Clusters"}; // Minmimum

  Configurable<float> cfgITSChi2NCl{"cfgITSChi2NCl", 999.0, "ITS Chi2/NCl"}; // Maximum
  Configurable<float> cfgTPCChi2NCl{"cfgTPCChi2NCl", 999.0, "TPC Chi2/NCl"}; // Maximum

  Configurable<bool> cfgUseTPCRefit{"cfgUseTPCRefit", true, "Require TPC Refit"};
  Configurable<bool> cfgUseITSRefit{"cfgUseITSRefit", true, "Require ITS Refit"};

  Configurable<bool> cfgHasITS{"cfgHasITS", true, "Require ITS"};
  Configurable<bool> cfgHasTPC{"cfgHasTPC", true, "Require TPC"};
  Configurable<bool> cfgHasTOF{"cfgHasTOF", true, "Require TOF"};

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  // DCA selections //

  // Primary track DCAr to PV
  Configurable<double> cMaxDCArToPVcut{"cMaxDCArToPVcut", 0.5, "Track DCAr cut to PV Maximum"};

  // Primary track DCAz to PV
  Configurable<double> cMaxDCAzToPVcut{"cMaxDCAzToPVcut", 2.0, "Track DCAz cut to PV Maximum"};
  Configurable<double> cMinDCAzToPVcut{"cMinDCAzToPVcut", 0.0, "Track DCAz cut to PV Minimum"};

  // Topological selections for V0
  Configurable<double> cDCALambdaDaugtherscut{"cDCALambdaDaugtherscut", 1.4, "Lambda dauthers DCA cut"};
  Configurable<double> cDCALambdaToPVcut{"cDCALambdaToPVcut", 0.07, "Lambda DCA cut to PV"};
  Configurable<double> cDCAPionToPVcut{"cMinDCApion", 0.05, "pion DCA cut to PV"};
  Configurable<double> cDCAProtonToPVcut{"cMinDCAproton", 0.05, "proton DCA cut to PV"};
  Configurable<double> cCosV0cut{"cCosV0cut", 0.97, "Cosine Pointing angle for V0"};
  Configurable<double> cMaxV0radiuscut{"cMaxV0radiuscut", 100., "V0 radius cut Maximum"};
  Configurable<double> cMinV0radiuscut{"cMinV0radiuscut", 0.2, "V0 radius cut Minimum"};
  // Configurable<double> cMasswindowV0cut{"cV0Masswindowcut", 0.007, "V0 Mass window cut"}; // How to ?

  // Topological selections for Cascade
  Configurable<double> cDCABachlorToPVcut{"cDCABachlorToPVcut", 0.015, "Bachelor DCA cut to PV"};
  Configurable<double> cDCAXiDaugtherscut{"cDCAXiDaugtherscut", 1.6, "Xi- DCA cut to PV"};
  Configurable<double> cCosPACasc{"cCosPACasc", 0.97, "Cosine Pointing angle for Cascade"};
  Configurable<double> cMaxCascradiuscut{"cMaxCascradiuscut", 100., "Cascade radius cut Maximum"};
  Configurable<double> cMinCascradiuscut{"cMinCascradiuscut", 0.2, "Cascade radius cut Minimum"};
  Configurable<double> cMasswindowCasccut{"cMasswindowCasccut", 0.007, "Cascade Mass window cut"};

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  // PID Selections//

  // PID Selections for Pion First
  Configurable<double> cMaxTPCnSigmaPionFirst{"cMaxTPCnSigmaPionFirst", 3.0, "TPC nSigma cut for Pion First"};
  Configurable<double> cMaxTOFnSigmaPionFirst{"cMaxTOFnSigmaPionFirst", 3.0, "TOF nSigma cut for Pion First"};

  Configurable<double> nsigmaCutCombinedPionFirst{"nsigmaCutCombinedPionFirst", -999, "Combined nSigma cut for Pion First"};

  Configurable<bool> cUseOnlyTOFTrackPionFirst{"cUseOnlyTOFTrackPionFirst", true, "Use only TOF track for PID selection Pion First"};
  Configurable<bool> cByPassTOFPionFirst{"cByPassTOFPionFirst", true, "By pass TOF Pion First PID selection"};

  // // PID Selections for Pion Bachelor
  // Configurable<double> cMaxTPCnSigmaPionBachelor{"cMaxTPCnSigmaPionBachelor", 3.0, "TPC nSigma cut for Pion Bachelor"};
  // Configurable<double> cMaxTOFnSigmaPionBachelor{"cMaxTOFnSigmaPionBachelor", 3.0, "TOF nSigma cut for Pion Bachelor"};

  // Configurable<double> nsigmaCutCombinedPionBachelor{"nsigmaCutCombinedPionBachelor", -999, "Combined nSigma cut for Pion Bachelor"};

  // Configurable<bool> cUseOnlyTOFTrackPionBachelor{"cUseOnlyTOFTrackPionBachelor", false, "Use only TOF track for PID selection Pion Bachelor"};

  // // PID Selections for Pion
  // Configurable<double> cMaxTPCnSigmaPion{"cMaxTPCnSigmaPion", 3.0, "TPC nSigma cut for Pion"};
  // Configurable<double> cMaxTOFnSigmaPion{"cMaxTOFnSigmaPion", 3.0, "TOF nSigma cut for Pion"};

  // Configurable<double> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", -999, "Combined nSigma cut for Pion"};

  // Configurable<bool> cUseOnlyTOFTrackPion{"cUseOnlyTOFTrackPion", false, "Use only TOF track for PID selection Pion"};

  // // PID Selections for Proton
  // Configurable<double> cMaxTPCnSigmaProton{"cMaxTPCnSigmaProton", 3.0, "TPC nSigma cut for Proton"};
  // Configurable<double> cMaxTOFnSigmaProton{"cMaxTOFnSigmaProton", 3.0, "TOF nSigma cut for Proton"};

  // Configurable<double> nsigmaCutCombinedProton{"nsigmaCutCombinedProton", -999, "Combined nSigma cut for Proton"};

  // Configurable<bool> cUseOnlyTOFTrackProton{"cUseOnlyTOFTrackProton", false, "Use only TOF track for PID selection Proton"};

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  // MC Event selection //
  Configurable<float> cZvertCutMC{"cZvertCutMC", 10.0, "MC Z-vertex cut"};

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  // Cuts on mother particle
  Configurable<bool> cfgCutsOnMother{"cfgCutsOnMother", true, "Enamble additional cuts on mother"};
  Configurable<double> cMaxPtMotherCut{"cMaxPtMotherCut", 15.0, "Maximum pt of mother cut"};
  Configurable<double> cMaxMinvMotherCut{"cMaxMinvMotherCut", 1.5, "Maximum Minv of mother cut"};
  Configurable<int> cetaphiBins{"cetaphiBins", 400, "number of eta and phi bins"};
  TRandom* rn = new TRandom();

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  void init(o2::framework::InitContext&)
  {
    AxisSpec centAxis = {binsCent, "FT0M (%)"};
    AxisSpec dcaxyAxis = {cDCABins, 0.0, 3.0, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {cDCABins, 0.0, 3.0, "DCA_{#it{z}} (cm)"};
    AxisSpec mcLabelAxis = {5, -0.5, 4.5, "MC Label"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {cInvMassBins, cInvMassStart, cInvMassEnd, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec pidQAAxis = {cPIDBins, -cPIDQALimit, cPIDQALimit};
    AxisSpec FlagAxis = {9, 0, 9, "Flags"};

    if (additionalQAeventPlots) {
      // Test on Mixed event
      histos.add("TestME/hCollisionIndexSameE", "coll index sameE", HistType::kTH1F, {{500, 0.0f, 500.0f}});
      histos.add("TestME/hCollisionIndexMixedE", "coll index mixedE", HistType::kTH1F, {{500, 0.0f, 500.0f}});

      histos.add("TestME/hnTrksSameE", "n tracks per event SameE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});
      histos.add("TestME/hnTrksMixedE", "n tracks per event MixedE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});

      histos.add("TestME/hPairsCounterSameE", "tot n pairs sameE", HistType::kTH1F, {{1, 0.5f, 1.5f}});
      histos.add("TestME/hPairsCounterMixedE", "tot n pairs mixedE", HistType::kTH1F, {{1, 0.5f, 1.5f}});

      // event histograms
      histos.add("QAevent/hEvtCounterSameE", "Number of analyzed Same Events", HistType::kTH1F, {{1, 0.5, 1.5}});
      histos.add("QAevent/hVertexZSameE", "Collision Vertex Z position", HistType::kTH1F, {{100, -15., 15.}});
      histos.add("QAevent/hMultiplicityPercentSameE", "Multiplicity percentile of collision", HistType::kTH1F, {{120, 0.0f, 120.0f}});

      histos.add("QAevent/hEvtCounterMixedE", "Number of analyzed Mixed Events", HistType::kTH1F, {{1, 0.5, 1.5}});
      histos.add("QAevent/hVertexZMixedE", "Collision Vertex Z position", HistType::kTH1F, {{100, -15., 15.}});
      histos.add("QAevent/hMultiplicityPercentMixedE", "Multiplicity percentile of collision", HistType::kTH1F, {{120, 0.0f, 120.0f}});
    }

    if (invmass1D) {
      histos.add("Xi1530invmassDS", "Invariant mass of Xi(1530)0 differnt sign", kTH1F, {invMassAxis});
      histos.add("Xi1530invmassLS", "Invariant mass of Xi(1530)0 like sign", kTH1F, {invMassAxis});
      histos.add("Xi1530invmassME", "Invariant mass of Xi(1530)0 mixed event", kTH1F, {invMassAxis});

      if (study_antiparticle) {
        histos.add("Xi1530invmassDSAnti", "Invariant mass of Anti-Xi(1530)0 differnt sign", kTH1F, {invMassAxis});
        histos.add("Xi1530invmassLSAnti", "Invariant mass of Anti-Xi(1530)0 like sign", kTH1F, {invMassAxis});
      }
    }

    if (additionalMEPlots) {
      histos.add("Xi1530invmassME_DS", "Invariant mass of Xi(1530)0 mixed event DS", kTH1F, {invMassAxis});
      histos.add("Xi1530invmassME_DSAnti", "Invariant mass of Xi(1530)0 mixed event DSAnti", kTH1F, {invMassAxis});
    }

    if (additionalQAplots) {
      // TPC ncluster distirbutions
      histos.add("TPCncluster/TPCnclusterpifirst", "TPC ncluster distribution", kTH1F, {{160, 0, 160, "TPC nCluster"}});
      // histos.add("TPCncluster/TPCnclusterpr", "TPC ncluster distribution", kTH1F, {{160, 0, 160, "TPC nCluster"}}); // can't use TPC info. for cascades!!
      histos.add("TPCncluster/TPCnclusterPhipifirst", "TPC ncluster vs phi", kTH2F, {{160, 0, 160, "TPC nCluster"}, {63, 0, 6.28, "#phi"}});
      // histos.add("TPCncluster/TPCnclusterPhipr", "TPC ncluster vs phi", kTH2F, {{160, 0, 160, "TPC nCluster"}, {63, 0, 6.28, "#phi"}}); // can't use TPC info. for cascades!!

      // Multiplicity correlation calibrations
      histos.add("MultCalib/centglopi_before", "Centrality vs Global-Tracks", kTH2F, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Global Tracks"}});
      histos.add("MultCalib/GloPVpi_before", "Global tracks vs PV tracks", kTH2F, {{500, 0, 5000, "Global tracks"}, {500, 0, 5000, "PV tracks"}});

      histos.add("MultCalib/centglopi_after", "Centrality vs Global-Tracks", kTH2F, {{110, 0, 110, "Centrality"}, {500, 0, 5000, "Global Tracks"}});
      histos.add("MultCalib/GloPVpi_after", "Global tracks vs PV tracks", kTH2F, {{500, 0, 5000, "Global tracks"}, {500, 0, 5000, "PV tracks"}});
    }

    // DCA QA to candidates for first pion and Xi-
    histos.add("QAbefore/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAbefore/trkDCAxy_Xi", "DCAxy distribution of Xi- track candidates", HistType::kTH1F, {dcaxyAxis});

    histos.add("QAbefore/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAbefore/trkDCAz_Xi", "DCAz distribution of Xi- track candidates", HistType::kTH1F, {dcazAxis});

    histos.add("QAafter/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxis});
    histos.add("QAafter/trkDCAxy_Xi", "DCAxy distribution of Xi- track candidates", HistType::kTH1F, {dcaxyAxis});

    histos.add("QAafter/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxis});
    histos.add("QAafter/trkDCAz_Xi", "DCAz distribution of Xi- track candidates", HistType::kTH1F, {dcazAxis});

    // pT QA to candidates for first pion, Xi
    histos.add("QAbefore/trkpT_pi", "pT distribution of pion track candidates", kTH1F, {ptAxis});
    histos.add("QAbefore/trkpT_Xi", "pT distribution of Xi- track candidates", kTH1F, {ptAxis});

    histos.add("QAafter/trkpT_pi", "pT distribution of pion track candidates", kTH1F, {ptAxis});
    histos.add("QAafter/trkpT_Xi", "pT distribution of Xi- track candidates", kTH1F, {ptAxis});

    // Can't use PID for cascades !!
    if (PIDplots) {
      histos.add("QAbefore/TOF_TPC_Map_pi_first_all", "TOF + TPC Combined PID for Pion_{First};#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QAbefore/TOF_Nsigma_pi_first_all", "TOF NSigma for Pion_{First};#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAbefore/TPC_Nsigma_pi_first_all", "TPC NSigma for Pion_{First};#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

      // histos.add("QAbefore/TOF_TPC_Map_pr_all", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{proton};#sigma_{TPC}^{proton}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      // histos.add("QAbefore/TOF_Nsigma_pr_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{proton};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      // histos.add("QAbefore/TPC_Nsigma_pr_all", "TPC NSigma for proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{proton};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

      // PID QA after cuts
      histos.add("QAafter/TOF_TPC_Map_pi_first_all", "TOF + TPC Combined PID for Pion_{First};#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QAafter/TOF_Nsigma_pi_first_all", "TOF NSigma for Pion_{First};#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAafter/TPC_Nsigma_pi_first_all", "TPC NSigma for Pion_{First};#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

      //     histos.add("QAafter/TOF_TPC_Map_pr_all", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{proton};#sigma_{TPC}^{proton}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      //     histos.add("QAafter/TOF_Nsigma_pr_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{proton};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      //     histos.add("QAafter/TPC_Nsigma_pr_all", "TPC NSigma for proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{proton};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
    }

    // 3d histogram + Flags
    histos.add("h3Xi1530invmassDS", "Invariant mass of Xi(1530)0 differnt sign", kTHnSparseF, {centAxis, ptAxis, invMassAxis, FlagAxis});
    histos.add("h3Xi1530invmassLS", "Invariant mass of Xi(1530)0 same sign", kTHnSparseF, {centAxis, ptAxis, invMassAxis, FlagAxis});
    histos.add("h3Xi1530invmassME", "Invariant mass of Xi(1530)0 mixed event", kTHnSparseF, {centAxis, ptAxis, invMassAxis, FlagAxis});

    if (study_antiparticle) {
      histos.add("h3Xi1530invmassDSAnti", "Invariant mass of Anti-Xi(1530)0 differnt sign", kTHnSparseF, {centAxis, ptAxis, invMassAxis, FlagAxis});
      histos.add("h3Xi1530invmassLSAnti", "Invariant mass of Anti-Xi(1530)0 same sign", kTHnSparseF, {centAxis, ptAxis, invMassAxis, FlagAxis});
    }

    if (additionalMEPlots) {
      histos.add("h3Xi1530invmassME_DS", "Invariant mass of Xi(1530)0 mixed event DS", kTHnSparseF, {centAxis, ptAxis, invMassAxis, FlagAxis});
      histos.add("h3Xi1530invmassME_DSAnti", "Invariant mass of Xi(1530)0 mixed event DSAnti", kTHnSparseF, {centAxis, ptAxis, invMassAxis, FlagAxis});
    }

    if (doprocessMC) {
      // MC QA
      histos.add("QAMCTrue/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH1F, {dcaxyAxis});
      histos.add("QAMCTrue/trkDCAxy_Xi", "DCAxy distribution of Xi- track candidates", HistType::kTH1F, {dcaxyAxis});

      histos.add("QAMCTrue/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH1F, {dcazAxis});
      histos.add("QAMCTrue/trkDCAz_Xi", "DCAz distribution of Xi- track candidates", HistType::kTH1F, {dcazAxis});

      // Can't use PID for cascades !!
      if (PIDplots) {
        histos.add("QAMCTrue/TOF_Nsigma_pi_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
        histos.add("QAMCTrue/TPC_Nsigma_pi_all", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

        // histos.add("QAMCTrue/TOF_Nsigma_pr_all", "TOF NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Proton};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
        // histos.add("QAMCTrue/TPC_Nsigma_pr_all", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      }

      histos.add("h3RecXi1530invmass", "Invariant mass of Reconstructed MC Xi(1530)0", kTHnSparseF, {centAxis, ptAxis, invMassAxis, FlagAxis});
      histos.add("h3RecXi1530invmassAnti", "Invariant mass of Reconstructed MC Anti-Xi(1530)0", kTHnSparseF, {centAxis, ptAxis, invMassAxis, FlagAxis});

      histos.add("h3Xi1530Gen", "pT distribution of True MC Xi(1530)0", kTHnSparseF, {mcLabelAxis, ptAxis, centAxis, FlagAxis});
      histos.add("h3Xi1530GenAnti", "pT distribution of True MC Anti-Xi(1530)0", kTHnSparseF, {mcLabelAxis, ptAxis, centAxis, FlagAxis});

      histos.add("Xi1530Rec", "pT distribution of Reconstructed MC Xi(1530)0", kTH2F, {ptAxis, centAxis});
      histos.add("Xi1530RecAnti", "pT distribution of Reconstructed MC Anti-Xi(1530)0", kTH2F, {ptAxis, centAxis});
      histos.add("Xi1530Recinvmass", "Inv mass distribution of Reconstructed MC Xi(1530)0", kTH1F, {invMassAxis});
    }
  }

  double massPi = MassPionCharged;

  // Primary track selection for the first pion //
  template <typename TrackType>
  bool trackCut(const TrackType track)
  {
    if (std::abs(track.eta()) > cMaxetacut)
      return false;
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (std::abs(track.dcaXY()) > cMaxDCArToPVcut)
      return false;
    if (std::abs(track.dcaZ()) > cMaxDCAzToPVcut)
      return false;
    if (track.itsNCls() < cfgITScluster)
      return false;
    if (track.tpcNClsFound() < cfgTPCcluster)
      return false;
    if (track.tpcCrossedRowsOverFindableCls() < cfgRatioTPCRowsOverFindableCls)
      return false;
    if (track.itsChi2NCl() >= cfgITSChi2NCl)
      return false;
    if (track.tpcChi2NCl() >= cfgTPCChi2NCl)
      return false;
    if (cfgHasITS && !track.hasITS())
      return false;
    if (cfgHasTPC && !track.hasTPC())
      return false;
    if (cfgHasTOF && !track.hasTOF())
      return false;
    if (cfgUseITSRefit && !track.passedITSRefit())
      return false;
    if (cfgUseTPCRefit && !track.passedTPCRefit())
      return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    if (cfgGlobalWoDCATrack && !track.isGlobalTrackWoDCA())
      return false;
    if (cfgGlobalTrack && !track.isGlobalTrack())
      return false;

    return true;
  }

  // Primary track selection for cascades (Need to more informations for cascades!) //
  template <typename TracksTypeCasc>
  bool casctrackCut(const TracksTypeCasc track)
  {
    if (std::abs(track.eta()) > cMaxetacut)
      return false;
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (std::abs(track.dcaXYCascToPV()) > cMaxDCArToPVcut)
      return false;
    if (std::abs(track.dcaZCascToPV()) > cMaxDCAzToPVcut)
      return false;

    return true;
  }

  // Secondary track selection for cascades  //
  //(Not added yet!-> Need to more informations for cascades!)

  // Topological cuts for cascades
  template <typename TracksTypeCasc>
  bool casctopCut(const TracksTypeCasc track)
  {
    // Topological cuts for V0s
    if (std::abs(track.daughDCA()) > cDCALambdaDaugtherscut)
      return false;
    if (std::abs(track.dcav0topv()) < cDCALambdaToPVcut)
      return false;
    if (track.sign() < 0) {
      if (std::abs(track.dcanegtopv()) < cDCAPionToPVcut)
        return false;
      if (std::abs(track.dcapostopv()) < cDCAProtonToPVcut)
        return false;
    } else if (track.sign() > 0) {
      if (std::abs(track.dcanegtopv()) < cDCAProtonToPVcut)
        return false;
      if (std::abs(track.dcapostopv()) < cDCAPionToPVcut)
        return false;
    }
    if (track.v0CosPA() < cCosV0cut)
      return false;
    if (track.transRadius() > cMaxV0radiuscut || track.transRadius() < cMinV0radiuscut)
      return false;

    // Topological Cuts for Cascades
    if (track.dcabachtopv() < cDCABachlorToPVcut)
      return false;
    if (track.cascdaughDCA() > cDCAXiDaugtherscut)
      return false;
    if (track.cascCosPA() < cCosPACasc)
      return false;
    if (track.casctransRadius() > cMaxCascradiuscut || track.casctransRadius() < cMinCascradiuscut)
      return false;
    // if (std::abs(track.mXi() - pdgDB->Mass(3312)) > cMasswindowCasccut) // codes are not operated when using pdgDB->Mass() !!
    //     return false;
    if (std::abs(track.mXi() - cMassXiminus) > cMasswindowCasccut)
      return false;

    return true;
  }

  // PID selection for the First Pion //
  template <typename T>
  bool selectionPIDPionFirst(const T& candidate)
  {
    if (tof_at_high_pt) {
      if (candidate.hasTOF() && (std::abs(candidate.tofNSigmaPi()) < cMaxTOFnSigmaPionFirst)) {
        return true;
      }
      if (!candidate.hasTOF() && (std::abs(candidate.tpcNSigmaPi()) < cMaxTPCnSigmaPionFirst)) {
        return true;
      }
    } else {
      bool tpcPIDPassed{false}, tofPIDPassed{false};
      if (std::abs(candidate.tpcNSigmaPi()) < cMaxTPCnSigmaPionFirst) {
        tpcPIDPassed = true;
      }
      if (cByPassTOFPionFirst && tpcPIDPassed) {
        return true;
      }
      if (candidate.hasTOF()) {
        if (std::abs(candidate.tofNSigmaPi()) < cMaxTOFnSigmaPionFirst) {
          tofPIDPassed = true;
        }
        if ((nsigmaCutCombinedPionFirst > 0) && (candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi() + candidate.tofNSigmaPi() * candidate.tofNSigmaPi() < nsigmaCutCombinedPionFirst * nsigmaCutCombinedPionFirst)) {
          tofPIDPassed = true;
        }
      } else {
        tofPIDPassed = true;
      }
      if (tpcPIDPassed && tofPIDPassed) {
        return true;
      }
    }
    return true;
  }

  // PID selection for the Cascades // -> Does not use yet!
  // template <class TCascTracksTo, typename TCascade>
  // bool selectionPIDCascade(const TCascade& candidate)
  // {
  //     auto bachTrack = candidate.template bachelor_as<TCascTracksTo>();
  //     auto posTrack = candidate.template posTrack_as<TCascTracksTo>();
  //     auto negTrack = candidate.template negTrack_as<TCascTracksTo>();

  //     bool lConsistentWithLambdaPos = false;
  //     bool lConsistentWithLambdaNeg = false;
  //     bool lConsistentWithLambda = false;
  //     bool lConsistentWithXi = false;

  //     if (tof_at_high_pt) {
  //         if (bachTrack.hasTOF() && (std::abs(bachTrack.tofNSigmaPi()) < cMaxTOFnSigmaPionBachelor)) {
  //             lConsistentWithXi = true;
  //         }
  //         if (!bachTrack.hasTOF() && (std::abs(bachTrack.tpcNSigmaPi()) < cMaxTPCnSigmaPionBachelor)) {
  //             lConsistentWithXi = true;
  //         }
  //     } else {
  //         bool tpcPIDPassed{false}, tofPIDPassed{false};
  //         if (std::abs(bachTrack.tpcNSigmaPi()) < cMaxTPCnSigmaPionBachelor) {
  //             tpcPIDPassed = true;
  //         }
  //         if (cByPassTOF && tpcPIDPassed) {
  //             lConsistentWithXi = true;
  //         }
  //         if (bachTrack.hasTOF()) {
  //             if (std::abs(bachTrack.tofNSigmaPi()) < cMaxTPCnSigmaPionBachelor) {
  //                 tofPIDPassed = true;
  //             }
  //             if ((nsigmaCutCombinedPionBachelor > 0) && (bachTrack.tpcNSigmaPi() * bachTrack.tpcNSigmaPi() + bachTrack.tofNSigmaPi() * bachTrack.tofNSigmaPi() < nsigmaCutCombinedPionBachelor * nsigmaCutCombinedPionBachelor)) {
  //                 tofPIDPassed = true;
  //             }
  //         } else {
  //             tofPIDPassed = true;
  //         }
  //         if (tpcPIDPassed && tofPIDPassed) {
  //             lConsistentWithXi = true;
  //         }

  //     }

  //     if(candidate.sign() > 0) {
  //         bool lConsistentWithPion = false;
  //         bool lConsistentWithProton = false;

  //         if (tof_at_high_pt) {
  //             if (posTrack.hasTOF() && (std::abs(posTrack.tofNSigmaPi()) < cMaxTOFnSigmaPion)) {
  //                 lConsistentWithPion = true;
  //             }
  //             if (!posTrack.hasTOF() && (std::abs(posTrack.tpcNSigmaPi()) < cMaxTPCnSigmaPion)) {
  //                 lConsistentWithPion = true;
  //             }
  //         } else {
  //             bool tpcPIDPassed{false}, tofPIDPassed{false};
  //             if (std::abs(posTrack.tpcNSigmaPi()) < cMaxTPCnSigmaPion) {
  //                 tpcPIDPassed = true;
  //             }
  //             if (cByPassTOF && tpcPIDPassed) {
  //                 lConsistentWithPion = true;
  //             }
  //             if (posTrack.hasTOF()) {
  //                 if (std::abs(posTrack.tofNSigmaPi()) < cMaxTPCnSigmaPion) {
  //                     tofPIDPassed = true;
  //                 }
  //                 if ((nsigmaCutCombinedPion > 0) && (posTrack.tpcNSigmaPi() * posTrack.tpcNSigmaPi() + posTrack.tofNSigmaPi() * posTrack.tofNSigmaPi() < nsigmaCutCombinedPion * nsigmaCutCombinedPion)) {
  //                     tofPIDPassed = true;
  //                 }
  //             } else {
  //                 tofPIDPassed = true;
  //             }
  //             if (tpcPIDPassed && tofPIDPassed) {
  //                 lConsistentWithPion = true;
  //             }
  //         }

  //         if (tof_at_high_pt) {
  //             if (posTrack.hasTOF() && (std::abs(posTrack.tofNSigmaPr()) < cMaxTOFnSigmaProton)) {
  //                 lConsistentWithProton = true;
  //             }
  //             if (!posTrack.hasTOF() && (std::abs(posTrack.tpcNSigmaPr()) < cMaxTPCnSigmaProton)) {
  //                 lConsistentWithProton = true;
  //             }
  //         } else {
  //             bool tpcPIDPassed{false}, tofPIDPassed{false};
  //             if (std::abs(posTrack.tpcNSigmaPr()) < cMaxTPCnSigmaProton) {
  //                 tpcPIDPassed = true;
  //             }
  //             if (cByPassTOF && tpcPIDPassed) {
  //                 lConsistentWithProton = true;
  //             }
  //             if (posTrack.hasTOF()) {
  //                 if (std::abs(posTrack.tofNSigmaPr()) < cMaxTPCnSigmaProton) {
  //                     tofPIDPassed = true;
  //                 }
  //                 if ((nsigmaCutCombinedProton > 0) && (posTrack.tpcNSigmaPr() * posTrack.tpcNSigmaPr() + posTrack.tofNSigmaPr() * posTrack.tofNSigmaPr() < nsigmaCutCombinedProton * nsigmaCutCombinedProton)) {
  //                     tofPIDPassed = true;
  //                 }
  //             } else {
  //                 tofPIDPassed = true;
  //             }
  //             if (tpcPIDPassed && tofPIDPassed) {
  //                 lConsistentWithProton = true;
  //             }
  //         }
  //         lConsistentWithLambdaPos = lConsistentWithProton*lConsistentWithPion;
  //     }

  //     if(candidate.sign() < 0) {
  //         bool lConsistentWithPion = false;
  //         bool lConsistentWithProton = false;

  //         if (tof_at_high_pt) {
  //             if (negTrack.hasTOF() && (std::abs(negTrack.tofNSigmaPi()) < cMaxTOFnSigmaPion)) {
  //                 lConsistentWithPion = true;
  //             }
  //             if (!negTrack.hasTOF() && (std::abs(negTrack.tpcNSigmaPi()) < cMaxTPCnSigmaPion)) {
  //                 lConsistentWithPion = true;
  //             }
  //         } else {
  //             bool tpcPIDPassed{false}, tofPIDPassed{false};
  //             if (std::abs(negTrack.tpcNSigmaPi()) < cMaxTPCnSigmaPion) {
  //                 tpcPIDPassed = true;
  //             }
  //             if (cByPassTOF && tpcPIDPassed) {
  //                 lConsistentWithPion = true;
  //             }
  //             if (negTrack.hasTOF()) {
  //                 if (std::abs(negTrack.tofNSigmaPi()) < cMaxTPCnSigmaPion) {
  //                     tofPIDPassed = true;
  //                 }
  //                 if ((nsigmaCutCombinedPion > 0) && (negTrack.tpcNSigmaPi() * negTrack.tpcNSigmaPi() + negTrack.tofNSigmaPi() * negTrack.tofNSigmaPi() < nsigmaCutCombinedPion * nsigmaCutCombinedPion)) {
  //                     tofPIDPassed = true;
  //                 }
  //             } else {
  //                 tofPIDPassed = true;
  //             }
  //             if (tpcPIDPassed && tofPIDPassed) {
  //                 lConsistentWithPion = true;
  //             }
  //         }

  //         if (tof_at_high_pt) {
  //             if (negTrack.hasTOF() && (std::abs(negTrack.tofNSigmaPr()) < cMaxTOFnSigmaProton)) {
  //                 lConsistentWithProton = true;
  //             }
  //             if (!negTrack.hasTOF() && (std::abs(negTrack.tpcNSigmaPr()) < cMaxTPCnSigmaProton)) {
  //                 lConsistentWithProton = true;
  //             }
  //         } else {
  //             bool tpcPIDPassed{false}, tofPIDPassed{false};
  //             if (std::abs(negTrack.tpcNSigmaPr()) < cMaxTPCnSigmaProton) {
  //                 tpcPIDPassed = true;
  //             }
  //             if (cByPassTOF && tpcPIDPassed) {
  //                 lConsistentWithProton = true;
  //             }
  //             if (negTrack.hasTOF()) {
  //                 if (std::abs(negTrack.tofNSigmaPr()) < cMaxTPCnSigmaProton) {
  //                     tofPIDPassed = true;
  //                 }
  //                 if ((nsigmaCutCombinedProton > 0) && (negTrack.tpcNSigmaPr() * negTrack.tpcNSigmaPr() + negTrack.tofNSigmaPr() * negTrack.tofNSigmaPr() < nsigmaCutCombinedProton * nsigmaCutCombinedProton)) {
  //                     tofPIDPassed = true;
  //                 }
  //             } else {
  //                 tofPIDPassed = true;
  //             }
  //             if (tpcPIDPassed && tofPIDPassed) {
  //                 lConsistentWithProton = true;
  //             }
  //         }
  //         lConsistentWithLambdaNeg = lConsistentWithProton*lConsistentWithPion;
  //     }

  //     lConsistentWithLambda = lConsistentWithLambdaPos * lConsistentWithLambdaNeg;

  //     return lConsistentWithXi * lConsistentWithLambda;
  // }

  template <bool IsMC, bool IsMix, typename CollisionType, typename TracksType, typename TracksTypeCasc>
  void fillHistograms(const CollisionType& collision, const TracksType& dTracks1, const TracksTypeCasc& dTracks2) // Order: ResoColl, ResoTrack, ResoCascTrack
  {
    auto multiplicity = collision.cent();

    if (additionalQAplots) {
      histos.fill(HIST("MultCalib/centglopi_before"), multiplicity, dTracks1.size());            // centrality vs global tracks before the multiplicity calibration cuts
      histos.fill(HIST("MultCalib/GloPVpi_before"), dTracks1.size(), collision.multNTracksPV()); // global tracks vs PV tracks before the multiplicity calibration cuts
    }

    if (additionalQAeventPlots) {
      if constexpr (!IsMix) {
        histos.fill(HIST("QAevent/hVertexZSameE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentSameE"), collision.cent());
        histos.fill(HIST("TestME/hCollisionIndexSameE"), collision.globalIndex());
        histos.fill(HIST("TestME/hnTrksSameE"), dTracks1.size());
      } else {
        histos.fill(HIST("QAevent/hVertexZMixedE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentMixedE"), collision.cent());
        histos.fill(HIST("TestME/hCollisionIndexMixedE"), collision.globalIndex());
        histos.fill(HIST("TestME/hnTrksMixedE"), dTracks1.size());
      }
    }

    if (additionalQAplots) {
      histos.fill(HIST("MultCalib/centglopi_after"), multiplicity, dTracks1.size());            // centrality vs global tracks after the multiplicity calibration cuts
      histos.fill(HIST("MultCalib/GloPVpi_after"), dTracks1.size(), collision.multNTracksPV()); // global tracks vs PV tracks after the multiplicity calibration cuts
    }

    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;
    for (auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(dTracks1, dTracks2))) {

      if (additionalQAeventPlots) {
        if constexpr (!IsMix) {
          histos.fill(HIST("TestME/hPairsCounterSameE"), 1.0);
        } else {
          histos.fill(HIST("TestME/hPairsCounterMixedE"), 1.0);
        }
      }

      if (!trackCut(trk1) || !casctrackCut(trk2))
        continue;

      auto isTrk1hasTOF = trk1.hasTOF();
      // auto isTrk2hasTOF = trk2.hasTOF(); // Can't use TPC, TOF info. for cascades yet!!

      auto trk1ptPi = trk1.pt();
      auto trk1NSigmaPiTPC = trk1.tpcNSigmaPi();
      auto trk1NSigmaPiTOF = (isTrk1hasTOF) ? trk1.tofNSigmaPi() : -999.;

      auto trk2ptXi = trk2.pt();

      if constexpr (!IsMix) {
        //// QA plots before the selection
        //  --- PID QA Pion
        if (PIDplots) {
          histos.fill(HIST("QAbefore/TPC_Nsigma_pi_first_all"), multiplicity, trk1ptPi, trk1NSigmaPiTPC);
          if (isTrk1hasTOF) {
            histos.fill(HIST("QAbefore/TOF_Nsigma_pi_first_all"), multiplicity, trk1ptPi, trk1NSigmaPiTOF);
            histos.fill(HIST("QAbefore/TOF_TPC_Map_pi_first_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
          }
        }

        histos.fill(HIST("QAbefore/trkpT_pi"), trk1ptPi);
        histos.fill(HIST("QAbefore/trkpT_Xi"), trk2ptXi);
        histos.fill(HIST("QAbefore/trkDCAxy_pi"), trk1.dcaXY());
        histos.fill(HIST("QAbefore/trkDCAxy_Xi"), trk2.dcaXYCascToPV());
        histos.fill(HIST("QAbefore/trkDCAz_pi"), trk1.dcaZ());
        histos.fill(HIST("QAbefore/trkDCAz_Xi"), trk2.dcaZCascToPV());
      }

      if (cUseOnlyTOFTrackPionFirst && !isTrk1hasTOF)
        continue;

      // if (cUseOnlyTOFTrackPr && !isTrk2hasTOF) // Can't use TPC, TOF info. for cascades yet!!
      //     continue;
      if (!selectionPIDPionFirst(trk1))
        continue;

      // if (!selectionPIDCascades(trk2)) // Can't use TPC, TOF info. for cascades yet!!
      //     continue;
      if (!casctopCut(trk2))
        continue;

      if (additionalQAplots) {
        // TPCncluster distributions
        histos.fill(HIST("TPCncluster/TPCnclusterpifirst"), trk1.tpcNClsFound());
        histos.fill(HIST("TPCncluster/TPCnclusterPhipifirst"), trk1.tpcNClsFound(), trk1.phi());
      }

      if constexpr (!IsMix) {
        //// QA plots before the selection
        //  --- PID QA Pion
        if (PIDplots) {
          histos.fill(HIST("QAafter/TPC_Nsigma_pi_first_all"), multiplicity, trk1ptPi, trk1NSigmaPiTPC);
          if (isTrk1hasTOF) {
            histos.fill(HIST("QAafter/TOF_Nsigma_pi_first_all"), multiplicity, trk1ptPi, trk1NSigmaPiTOF);
            histos.fill(HIST("QAafter/TOF_TPC_Map_pi_first_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
          }
        }

        histos.fill(HIST("QAafter/trkpT_pi"), trk1ptPi);
        histos.fill(HIST("QAafter/trkpT_Xi"), trk2ptXi);
        histos.fill(HIST("QAafter/trkDCAxy_pi"), trk1.dcaXY());
        histos.fill(HIST("QAafter/trkDCAxy_Xi"), trk2.dcaXYCascToPV());
        histos.fill(HIST("QAafter/trkDCAz_pi"), trk1.dcaZ());
        histos.fill(HIST("QAafter/trkDCAz_Xi"), trk2.dcaZCascToPV());
      }

      lDecayDaughter1.SetPtEtaPhiM(trk1.pt(), trk1.eta(), trk1.phi(), massPi);
      lDecayDaughter2.SetPtEtaPhiM(trk2.pt(), trk2.eta(), trk2.phi(), trk2.mXi());
      lResonance = lDecayDaughter1 + lDecayDaughter2;

      if (abs(lResonance.Rapidity()) >= 0.5)
        continue;

      if (cfgCutsOnMother) {
        if (lResonance.Pt() >= cMaxPtMotherCut) // excluding candidates in overflow
          continue;
        if (lResonance.M() >= cMaxMinvMotherCut) // excluding candidates in overflow
          continue;
      }

      if (trk1.sign() * trk2.sign() < 0) {

        if constexpr (!IsMix) {
          if (study_antiparticle) {
            if (trk1.sign() > 0) {
              if (invmass1D)
                histos.fill(HIST("Xi1530invmassDS"), lResonance.M());
              histos.fill(HIST("h3Xi1530invmassDS"), multiplicity, lResonance.Pt(), lResonance.M(), kData);
            } else if (trk1.sign() < 0) {
              if (invmass1D)
                histos.fill(HIST("Xi1530invmassDSAnti"), lResonance.M());
              histos.fill(HIST("h3Xi1530invmassDSAnti"), multiplicity, lResonance.Pt(), lResonance.M(), kData);
            }
          } else {
            if (invmass1D)
              histos.fill(HIST("Xi1530invmassDS"), lResonance.M());
            histos.fill(HIST("h3Xi1530invmassDS"), multiplicity, lResonance.Pt(), lResonance.M(), kData);
          }
        } else {
          if (invmass1D)
            histos.fill(HIST("Xi1530invmassME"), lResonance.M());
          if (trk1.sign() > 0) {
            if (invmass1D)
              histos.fill(HIST("Xi1530invmassME_DS"), lResonance.M());
            histos.fill(HIST("h3Xi1530invmassME_DS"), multiplicity, lResonance.Pt(), lResonance.M(), kMixing);
          } else if (trk1.sign() < 0) {
            if (invmass1D)
              histos.fill(HIST("Xi1530invmassME_DSAnti"), lResonance.M());
            histos.fill(HIST("h3Xi1530invmassME_DSAnti"), multiplicity, lResonance.Pt(), lResonance.M(), kMixing);
          }
          histos.fill(HIST("h3Xi1530invmassME"), multiplicity, lResonance.Pt(), lResonance.M(), kMixing);
        }

        if constexpr (IsMC) {
          if (abs(trk1.pdgCode()) != 211 || abs(trk2.pdgCode()) != 3312)
            continue;
          if (trk1.motherId() != trk2.motherId())
            continue;
          if (abs(trk2.motherPDG()) != 3324)
            continue;

          histos.fill(HIST("QAMCTrue/trkDCAxy_pi"), trk1.dcaXY());
          histos.fill(HIST("QAMCTrue/trkDCAxy_xi"), trk2.dcaXYCascToPV());
          histos.fill(HIST("QAMCTrue/trkDCAz_pi"), trk1.dcaZ());
          histos.fill(HIST("QAMCTrue/trkDCAz_xi"), trk2.dcaZCascToPV());
          histos.fill(HIST("QAMCTrue/TPC_Nsigma_pi_all"), multiplicity, trk1ptPi, trk1NSigmaPiTPC);
          if (isTrk1hasTOF) {
            histos.fill(HIST("QAMCTrue/TOF_Nsigma_pi_all"), multiplicity, trk1ptPi, trk1NSigmaPiTOF);
          }

          // MC histograms
          if (trk2.motherPDG() > 0) {
            histos.fill(HIST("Xi1530Rec"), lResonance.Pt(), multiplicity);
            histos.fill(HIST("Xi1530Recinvmass"), lResonance.M());
            histos.fill(HIST("h3RecXi1530invmass"), multiplicity, lResonance.Pt(), lResonance.M(), kMCReco);
          } else {
            histos.fill(HIST("Xi1530RecAnti"), lResonance.Pt(), multiplicity);
            histos.fill(HIST("Xi1530Recinvmass"), lResonance.M());
            histos.fill(HIST("h3RecXi1530invmassAnti"), multiplicity, lResonance.Pt(), lResonance.M(), kMCReco);
          }
        }

      } else {
        if constexpr (!IsMix) {
          if (study_antiparticle) {
            if (trk1.sign() < 0) {
              if (invmass1D)
                histos.fill(HIST("Xi1530invmassLS"), lResonance.M());
              histos.fill(HIST("h3Xi1530invmassLS"), multiplicity, lResonance.Pt(), lResonance.M(), kLS);
            } else if (trk1.sign() > 0) {
              if (invmass1D)
                histos.fill(HIST("Xi1530invmassLSAnti"), lResonance.M());
              histos.fill(HIST("h3Xi1530invmassLSAnti"), multiplicity, lResonance.Pt(), lResonance.M(), kLS);
            }
          } else {
            if (invmass1D)
              histos.fill(HIST("Xi1530invmassLS"), lResonance.M());
            histos.fill(HIST("h3Xi1530invmassLS"), multiplicity, lResonance.Pt(), lResonance.M(), kLS);
          }
        }
      }
    }
  }

  void processData(aod::ResoCollision& resoCollision, aod::ResoTracks const& resoTracks, aod::ResoCascades const& cascTracks)
  {
    if (additionalQAeventPlots)
      histos.fill(HIST("QAevent/hEvtCounterSameE"), 1.0);
    fillHistograms<false, false>(resoCollision, resoTracks, cascTracks);
  }

  void processMC(ResoMCCols::iterator const& resoCollision,
                 soa::Join<aod::ResoCascades, aod::ResoMCCascades> const& cascTracks,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resoTracks)
  {
    if (!resoCollision.isInAfterAllCuts() || (abs(resoCollision.posZ()) > cZvertCutMC)) // MC event selection, all cuts missing vtx cut
      return;
    fillHistograms<true, false>(resoCollision, resoTracks, cascTracks);
  }

  void processMCTrue(ResoMCCols::iterator const& resoCollision, aod::ResoMCParents& resoParents)
  {
    auto multiplicity = resoCollision.cent();
    for (auto& part : resoParents) { // loop over all pre-filtered MC particles
      if (abs(part.pdgCode()) != 3324 || abs(part.y()) >= 0.5)
        continue;
      bool pass1 = abs(part.daughterPDG1()) == 211 || abs(part.daughterPDG2()) == 211;
      bool pass2 = abs(part.daughterPDG1()) == 3312 || abs(part.daughterPDG2()) == 3312;

      if (!pass1 || !pass2)
        continue;

      if (resoCollision.isVtxIn10()) // INEL10
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("h3Xi1530Gen"), 0, part.pt(), multiplicity, kINEL10);
        else
          histos.fill(HIST("h3Xi1530GenAnti"), 0, part.pt(), multiplicity, kINEL10);
      }
      if (resoCollision.isVtxIn10() && resoCollision.isInSel8()) // INEL>10, vtx10
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("h3Xi1530Gen"), 1, part.pt(), multiplicity, kINELg010);
        else
          histos.fill(HIST("h3Xi1530GenAnti"), 1, part.pt(), multiplicity, kINELg010);
      }
      if (resoCollision.isVtxIn10() && resoCollision.isTriggerTVX()) // vtx10, TriggerTVX
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("h3Xi1530Gen"), 2, part.pt(), multiplicity, kMCTruePS);
        else
          histos.fill(HIST("h3Xi1530GenAnti"), 2, part.pt(), multiplicity, kMCTruePS);
      }
      if (resoCollision.isInAfterAllCuts()) // after all event selection
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("h3Xi1530Gen"), 3, part.pt(), multiplicity, kAllType);
        else
          histos.fill(HIST("h3Xi1530GenAnti"), 3, part.pt(), multiplicity, kAllType);
      }
    }
  }

  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processME(aod::ResoCollisions const& resoCollisions, aod::ResoTracks const& resoTracks, aod::ResoCascades const& cascTracks)
  {
    auto tracksTuple = std::make_tuple(resoTracks, cascTracks);

    BinningTypeVtxZT0M colBinning{{CfgVtxBins, CfgMultBins}, true};
    Pair<aod::ResoCollisions, aod::ResoTracks, aod::ResoCascades, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, resoCollisions, tracksTuple, &cache};

    for (auto& [collision1, tracks1, collision2, tracks2] : pairs) {
      if (additionalQAeventPlots)
        histos.fill(HIST("QAevent/hEvtCounterMixedE"), 1.0);
      fillHistograms<false, true>(collision1, tracks1, tracks2);
    }
  }

  PROCESS_SWITCH(xi1530analysisqa, processData, "Process Event for Data", true);
  PROCESS_SWITCH(xi1530analysisqa, processMC, "Process Event for MC (Reconstructed)", true);
  PROCESS_SWITCH(xi1530analysisqa, processMCTrue, "Process Event for MC (Generated)", true);
  PROCESS_SWITCH(xi1530analysisqa, processME, "Process EventMixing light without partition", true);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<xi1530analysisqa>(cfgc)};
}
