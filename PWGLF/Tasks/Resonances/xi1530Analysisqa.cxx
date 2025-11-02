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

/// \file xi1530Analysisqa.cxx
/// \brief Reconstruction of Xi* resonance.
///
/// \author Min-jae Kim <minjae.kim@cern.ch>, Bong-Hwi Lim <bong-hwi.lim@cern.ch>
// #include <TLorentzVector.h>
#include "PWGLF/DataModel/LFResonanceTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/mcCentrality.h"
#include "PWGLF/Utils/inelGt.h"

#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"

#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/runDataProcessing.h"

#include "Math/Vector4D.h"
#include "TF1.h"
#include "TRandom3.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;
using LorentzVectorPtEtaPhiMass = ROOT::Math::PtEtaPhiMVector;
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
  kAllType,
  kXiStar = 3324
};

struct Xi1530Analysisqa {

  // Basic set-up //
  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  Preslice<aod::ResoMicroTracks> perResoCollision =
    aod::resodaughter::resoCollisionId;
  Preslice<aod::ResoCascades> perResoCollisionCasc =
    aod::resodaughter::resoCollisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  using ResoMCCols = soa::Join<aod::ResoCollisions, aod::ResoMCCollisions>;
  Configurable<float> cMassXiminus{"cMassXiminus", 1.32171, "Mass of Xi baryon"};

  // Associated with histograms
  ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
  ConfigurableAxis binsPtQA{"binsPtQA", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0}, "Binning of the pT axis"};
  ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0, 110.0}, "Binning of the centrality axis"};

  Configurable<float> cInvMassStart{"cInvMassStart", 1.4, "Invariant mass start"};
  Configurable<float> cInvMassEnd{"cInvMassEnd", 3.0, "Invariant mass end"};
  Configurable<int> cInvMassBins{"cInvMassBins", 800, "Invariant mass binning"};

  Configurable<bool> invMass1D{"invMass1D", true, "Invariant mass 1D"};
  Configurable<bool> studyAntiparticle{"studyAntiparticle", true, "Study anti-particles separately"};
  Configurable<bool> pidPlots{"pidPlots", true, "Make TPC and TOF PID plots"};
  Configurable<bool> additionalQAplots{"additionalQAplots", true, "Additional QA plots"};

  // Event Mixing
  Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of events to mix"};
  ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
  ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Mixing bins - z-vertex"};

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  // Track selections (Except DCA selelctions) //

  // Primary track selections
  Configurable<float> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};
  Configurable<float> cMaxetacut{"cMaxetacut", 0.8, "Track maximum eta cut"};

  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};
  Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};

  Configurable<bool> tofAtHighPt{"tofAtHighPt", false, "Use TOF at high pT"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 1, "Minimum Number of TPC cluster"}; // Minmimum

  Configurable<int> cfgTPCRows{"cfgTPCRows", 80, "Minimum Number of TPC Crossed Rows "};
  Configurable<float> cfgRatioTPCRowsOverFindableCls{"cfgRatioTPCRowsOverFindableCls", 0.8, "Minimum of TPC Crossed Rows to Findable Clusters"}; // Minmimum

  // Configurable<bool> cfgUseTPCRefit{"cfgUseTPCRefit", true, "Require TPC Refit"}; //refit is included in global track selection
  // Configurable<bool> cfgUseITSRefit{"cfgUseITSRefit", true, "Require ITS Refit"};

  Configurable<bool> cfgHasTOF{"cfgHasTOF", false, "Require TOF"};

  Configurable<float> cfgRapidityCut{"cfgRapidityCut", 0.5, "Rapidity cut for tracks"};

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  // Cascade and V0 selections //

  // Primary track DCAxy to PV
  Configurable<float> cDCAxytoPVByPtPiFirstP0{"cDCAxytoPVByPtPiFirstP0", 0.004, "Coeff. Track DCAxy cut to PV by pt for Pion First (p0)"};
  Configurable<float> cDCAxyToPVByPtPiFirstExp{"cDCAxyToPVByPtPiFirstExp", 0.013, "Coeff. Track DCAxy cut to PV by pt for Pion First (exp)"};

  Configurable<bool> cDCAxyToPVAsPtForCasc{"cDCAxyToPVAsPtForCasc", true, "Set DCAxy to PV selection as pt cut"};
  Configurable<float> cDCAxyToPVByPtCascP0{"cDCAxyToPVByPtCascP0", 999., "Coeff. for Track DCAxy cut to PV by pt for Cascade (p0)"};
  Configurable<float> cDCAxyToPVByPtCascExp{"cDCAxyToPVByPtCascExp", 1., "Coeff. Track DCAxy cut to PV by pt for Cascade (exp)"};

  // Primary track DCAz to PV
  Configurable<bool> cDCAzToPVAsPt{"cDCAzToPVAsPt", true, "DCAz to PV selection as pt"};
  Configurable<bool> cDCAzToPVAsPtForCasc{"cDCAzToPVAsPtForCasc", true, "Set DCA to PV selection as pt cut"};
  Configurable<float> cMaxDCAzToPVCut{"cMaxDCAzToPVCut", 0.5, "Track DCAz cut to PV Maximum"};
  Configurable<float> cMaxDCAzToPVCutCasc{"cMaxDCAzToPVCutCasc", 0.5, "Track DCAz cut to PV Maximum for casc"};

  // Topological selections for V0
  Configurable<float> cDCALambdaDaugtherscut{"cDCALambdaDaugtherscut", 0.7, "Lambda dauthers DCA cut"};
  Configurable<float> cDCALambdaToPVcut{"cDCALambdaToPVcut", 0.02, "Lambda DCA cut to PV"};
  Configurable<float> cDCAPionToPVcut{"cDCAPionToPVcut", 0.06, "pion DCA cut to PV"};
  Configurable<float> cDCAProtonToPVcut{"cDCAProtonToPVcut", 0.07, "proton DCA cut to PV"};

  Configurable<float> cV0CosPACutPtDepP0{"cV0CosPACutPtDepP0", 0.25, "Coeff. for Cosine Pointing angle for V0 as pt (p0)"};
  Configurable<float> cV0CosPACutPtDepP1{"cV0CosPACutPtDepP1", 0.022, "Coeff. for Cosine Pointing angle for V0 as pt (p1)"};

  Configurable<float> cMaxV0radiuscut{"cMaxV0radiuscut", 200., "V0 radius cut Maximum"};
  Configurable<float> cMinV0radiuscut{"cMinV0radiuscut", 2.5, "V0 radius cut Minimum"};
  Configurable<float> cMasswindowV0cut{"cMasswindowV0cut", 0.005, "V0 Mass window cut"};

  // Topological selections for Cascade

  Configurable<float> cDCABachlorToPVcut{"cDCABachlorToPVcut", 0.06, "Bachelor DCA cut to PV"};
  Configurable<float> cDCAXiDaugthersCutPtRangeLower{"cDCAXiDaugthersCutPtRangeLower", 1., "Xi- DCA cut to PV as pt range lower"};
  Configurable<float> cDCAXiDaugthersCutPtRangeUpper{"cDCAXiDaugthersCutPtRangeUpper", 4., "Xi- DCA cut to PV as pt range upper"};
  Configurable<float> cDCAXiDaugthersCutPtDepLower{"cDCAXiDaugthersCutPtDepLower", 0.8, "Xi- DCA cut to PV as pt Under 1 GeV/c"};
  Configurable<float> cDCAXiDaugthersCutPtDepMiddle{"cDCAXiDaugthersCutPtDepMiddle", 0.5, "Xi- DCA cut to PV as pt 1 - 4 GeV/c"};
  Configurable<float> cDCAXiDaugthersCutPtDepUpper{"cDCAXiDaugthersCutPtDepUpper", 0.2, "Xi- DCA cut to PV as pt Over 4 GeV/c"};

  Configurable<float> cCosPACascCutPtDepP0{"cCosPACascCutPtDepP0", 0.2, "Coeff. for Cosine Pointing angle for Cascade as pt (p0)"};
  Configurable<float> cCosPACascCutPtDepP1{"cCosPACascCutPtDepP1", 0.022, "Coeff. for Cosine Pointing angle for Cascade as pt (p1)"};

  Configurable<float> cMaxCascradiuscut{"cMaxCascradiuscut", 200., "Cascade radius cut Maximum"};
  Configurable<float> cMinCascradiuscut{"cMinCascradiuscut", 1.1, "Cascade radius cut Minimum"};
  Configurable<float> cMasswindowCasccut{"cMasswindowCasccut", 0.008, "Cascade Mass window cut"};

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  // PID Selections//

  Configurable<float> cPIDBound{"cPIDBound", 6.349, "configurable for replacing to .has"};

  Configurable<float> cMinTOFpt{"cMinTOFpt", 0.5, "Maximum TOF pt cut"};

  // PID Selections for Pion First
  Configurable<float> cMaxtpcnSigmaPionFirst{"cMaxtpcnSigmaPionFirst", 4.0, "TPC nSigma cut for Pion First"};
  Configurable<float> cMaxtofnSigmaPionFirst{"cMaxtofnSigmaPionFirst", 3.0, "TOF nSigma cut for Pion First"};

  Configurable<float> nsigmaCutCombinedPionFirst{"nsigmaCutCombinedPionFirst", -4.0, "Combined nSigma cut for Pion First"};

  Configurable<bool> cUseOnlyTOFTrackPionFirst{"cUseOnlyTOFTrackPionFirst", false, "Use only TOF track for PID selection Pion First"};
  Configurable<bool> cByPassTOFPionFirst{"cByPassTOFPionFirst", true, "By pass TOF Pion First PID selection"};

  // PID Selections for Pion Bachelor
  Configurable<float> cMaxtpcnSigmaPionBachelor{"cMaxtpcnSigmaPionBachelor", 4.0, "TPC nSigma cut for Pion Bachelor"};
  Configurable<float> cMaxtofnSigmaPionBachelor{"cMaxtofnSigmaPionBachelor", 3.0, "TOF nSigma cut for Pion Bachelor"};

  Configurable<float> nsigmaCutCombinedPionBachelor{"nsigmaCutCombinedPionBachelor", -4.0, "Combined nSigma cut for Pion Bachelor"};

  Configurable<bool> cUseOnlyTOFTrackPionBachelor{"cUseOnlyTOFTrackPionBachelor", false, "Use only TOF track for PID selection Pion Bachelor"};
  Configurable<bool> cByPassTOFPionBachelor{"cByPassTOFPionBachelor", true, "By pass TOF Pion Bachelor PID selection"};

  // PID Selections for Pion
  Configurable<float> cMaxtpcnSigmaPion{"cMaxtpcnSigmaPion", 4.0, "TPC nSigma cut for Pion"};
  Configurable<float> cMaxtofnSigmaPion{"cMaxtofnSigmaPion", 3.0, "TOF nSigma cut for Pion"};

  Configurable<float> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", -4.0, "Combined nSigma cut for Pion"};

  Configurable<bool> cUseOnlyTOFTrackPion{"cUseOnlyTOFTrackPion", false, "Use only TOF track for PID selection Pion"};
  Configurable<bool> cByPassTOFPion{"cByPassTOFPion", true, "By pass TOF Pion PID selection"};

  // PID Selections for Proton
  Configurable<float> cMaxtpcnSigmaProton{"cMaxtpcnSigmaProton", 4.0, "TPC nSigma cut for Proton"};
  Configurable<float> cMaxtofnSigmaProton{"cMaxtofnSigmaProton", 3.0, "TOF nSigma cut for Proton"};

  Configurable<float> nsigmaCutCombinedProton{"nsigmaCutCombinedProton", -4.0, "Combined nSigma cut for Proton"};

  Configurable<bool> cUseOnlyTOFTrackProton{"cUseOnlyTOFTrackProton", false, "Use only TOF track for PID selection Proton"};
  Configurable<bool> cByPassTOFProton{"cByPassTOFProton", true, "By pass TOF Proton PID selection"};

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  // MC Event selection //
  Configurable<float> cZvertCutMC{"cZvertCutMC", 10.0, "MC Z-vertex cut"};

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  // Cuts on mother particle and others
  Configurable<bool> cfgCutsOnMother{"cfgCutsOnMother", true, "Enamble additional cuts on mother"};
  Configurable<float> cMaxPtMotherCut{"cMaxPtMotherCut", 9.0, "Maximum pt of mother cut"};
  Configurable<float> cMaxMinvMotherCut{"cMaxMinvMotherCut", 3.0, "Maximum Minv of mother cut"};

  Configurable<bool> studyStableXi{"studyStableXi", false, "Study stable Xi"};

  Configurable<bool> cMCCent{"cMCCent", true, "Using calibrated MC centrality (for FT0M)"};
  Configurable<bool> cRecoINELgt0{"cRecoINELgt0", true, "check if INEL>0 for reco events"};
  TRandom* rn = new TRandom();

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  struct PidSelectionParam {
    double cMaxTPCnSigma;
    double cMaxTOFnSigma;
    bool cByPassTOF;
    double nsigmaCutCombined;
  };

  void init(o2::framework::InitContext&)
  {
    AxisSpec centAxis = {binsCent, "FT0M (%)"};
    AxisSpec dcaxyAxis = {1500, 0.0, 0.3, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {1500, 0.0, 0.3, "DCA_{#it{z}} (cm)"};
    AxisSpec dcaDaugAxis = {1000, 0.0, 1, "DCA_{#it{Daughter}} (cm)"};
    AxisSpec cosPAAxis = {3000, 0.0, 0.06, "1-cos(PA)"};
    AxisSpec mcLabelAxis = {6, -1.5, 4.5, "MC Label"};
    AxisSpec ptAxis = {binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {cInvMassBins, cInvMassStart, cInvMassEnd, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec invMassAxisCasc = {800, 1.25, 1.65, "Invariant Mass for Casc. (GeV/#it{c}^2)"};
    AxisSpec pidQAAxis = {65, -6.5, 6.5};
    AxisSpec flagAxis = {9, 0, 9, "Flags"};

    {
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
      histos.add("QAevent/hMultiplicityPercentSameE", "Multiplicity percentile of collision", HistType::kTH1F, {centAxis});

      histos.add("QAevent/hEvtCounterMixedE", "Number of analyzed Mixed Events", HistType::kTH1F, {{1, 0.5, 1.5}});
      histos.add("QAevent/hVertexZMixedE", "Collision Vertex Z position", HistType::kTH1F, {{100, -15., 15.}});
      histos.add("QAevent/hMultiplicityPercentMixedE", "Multiplicity percentile of collision", HistType::kTH1F, {centAxis});
    }

    if (invMass1D) {
      histos.add("Xi1530invmassDS", "Invariant mass of Xi(1530)0 differnt sign", kTH1F, {invMassAxis});
      histos.add("Xi1530invmassLS", "Invariant mass of Xi(1530)0 like sign", kTH1F, {invMassAxis});
      histos.add("Xi1530invmassME", "Invariant mass of Xi(1530)0 mixed event", kTH1F, {invMassAxis});

      if (studyAntiparticle) {
        histos.add("Xi1530invmassDSAnti", "Invariant mass of Anti-Xi(1530)0 differnt sign", kTH1F, {invMassAxis});
        histos.add("Xi1530invmassLSAnti", "Invariant mass of Anti-Xi(1530)0 like sign", kTH1F, {invMassAxis});
      }
    }

    if (doprocessMEDF || doprocessMEMicro) {
      histos.add("Xi1530invmassME_DS", "Invariant mass of Xi(1530)0 mixed event DS", kTH1F, {invMassAxis});
      histos.add("Xi1530invmassME_DSAnti", "Invariant mass of Xi(1530)0 mixed event DSAnti", kTH1F, {invMassAxis});
    }

    // DCA QA to candidates for first pion and Xi-
    histos.add("QAbefore/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH2F, {ptAxis, dcaxyAxis});
    histos.add("QAbefore/trkDCAxy_Xi", "DCAxy distribution of Xi- track candidates", HistType::kTH2F, {ptAxis, dcaxyAxis});

    histos.add("QAbefore/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH2F, {ptAxis, dcazAxis});
    histos.add("QAbefore/trkDCAz_Xi", "DCAz distribution of Xi- track candidates", HistType::kTH2F, {ptAxis, dcazAxis});

    histos.add("QAafter/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH2F, {ptAxis, dcaxyAxis});
    histos.add("QAafter/trkDCAxy_Xi", "DCAxy distribution of Xi- track candidates", HistType::kTH2F, {ptAxis, dcaxyAxis});

    histos.add("QAafter/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH2F, {ptAxis, dcazAxis});
    histos.add("QAafter/trkDCAz_Xi", "DCAz distribution of Xi- track candidates", HistType::kTH2F, {ptAxis, dcazAxis});

    if (pidPlots) {
      // Plots for QA before, Need to pt info. for the daugthers
      histos.add("QAbefore/TOF_TPC_Map_pi_first_all", "TOF + TPC Combined PID for Pion_{First};#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QAbefore/TOF_Nsigma_pi_first_all", "TOF NSigma for Pion_{First};#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAbefore/TPC_Nsigma_pi_first_all", "TPC NSigma for Pion_{First};#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

      histos.add("QAbefore/TOF_TPC_Map_pi_bachelor_all", "TOF + TPC Combined PID for Pion_{Bachelor};#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QAbefore/TPC_Nsigma_pi_bachelor_all", "TPC NSigma for Pion_{Bachelor};#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

      histos.add("QAbefore/TOF_TPC_Map_pr_all", "TOF + TPC Combined PID for Proton;#sigma_{TOF}^{Proton};#sigma_{TPC}^{Proton}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QAbefore/TOF_TPC_Map_antipr_all", "TOF + TPC Combined PID for Anti-Proton;#sigma_{TOF}^{Proton};#sigma_{TPC}^{Proton}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});

      histos.add("QAbefore/TPC_Nsigma_pr_all", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAbefore/TPC_Nsigma_antipr_all", "TPC NSigma for Anti-Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

      histos.add("QAbefore/TOF_TPC_Map_pi_all", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QAbefore/TOF_TPC_Map_piminus_all", "TOF + TPC Combined PID for Pion -;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});

      histos.add("QAbefore/TPC_Nsigma_pi_all", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAbefore/TPC_Nsigma_piminus_all", "TPC NSigma for Pion -;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

      // Plots for QA after
      histos.add("QAafter/TOF_TPC_Map_pi_first_all", "TOF + TPC Combined PID for Pion_{First};#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QAafter/TOF_Nsigma_pi_first_all", "TOF NSigma for Pion_{First};#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAafter/TPC_Nsigma_pi_first_all", "TPC NSigma for Pion_{First};#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

      histos.add("QAafter/TOF_TPC_Map_pi_bachelor_all", "TOF + TPC Combined PID for Pion_{Bachelor};#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QAafter/TPC_Nsigma_pi_bachelor_all", "TPC NSigma for Pion_{Bachelor};#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

      histos.add("QAafter/TOF_TPC_Map_pr_all", "TOF + TPC Combined PID for Proton;#sigma_{TOF}^{Proton};#sigma_{TPC}^{Proton}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QAafter/TOF_TPC_Map_antipr_all", "TOF + TPC Combined PID for Anti-Proton;#sigma_{TOF}^{Proton};#sigma_{TPC}^{Proton}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});

      histos.add("QAafter/TPC_Nsigma_pr_all", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAafter/TPC_Nsigma_antipr_all", "TPC NSigma for Anti-Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

      histos.add("QAafter/TOF_TPC_Map_pi_all", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
      histos.add("QAafter/TOF_TPC_Map_piminus_all", "TOF + TPC Combined PID for Pion -;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});

      histos.add("QAafter/TPC_Nsigma_pi_all", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      histos.add("QAafter/TPC_Nsigma_piminus_all", "TPC NSigma for Pion -;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
    }

    // 3d histogram + Flags
    histos.add("h3Xi1530invmassDS", "Invariant mass of Xi- differnt sign", kTHnSparseF, {centAxis, ptAxis, invMassAxis, flagAxis});
    histos.add("h3XiinvmassDS", "Invariant mass of Xi- differnt sign", kTHnSparseF, {centAxis, ptAxis, invMassAxisCasc, flagAxis});

    histos.add("h3Xi1530invmassLS", "Invariant mass of Xi(1530)0 same sign", kTHnSparseF, {centAxis, ptAxis, invMassAxis, flagAxis});
    histos.add("h3XiinvmassLS", "Invariant mass of Xi- same sign", kTHnSparseF, {centAxis, ptAxis, invMassAxisCasc, flagAxis});

    histos.add("h3Xi1530invmassME", "Invariant mass of Xi(1530)0 mixed event", kTHnSparseF, {centAxis, ptAxis, invMassAxis, flagAxis});
    histos.add("h3XiinvmassME", "Invariant mass of Xi- mixed event", kTHnSparseF, {centAxis, ptAxis, invMassAxisCasc, flagAxis});

    if (studyAntiparticle) {
      histos.add("h3Xi1530invmassDSAnti", "Invariant mass of Anti-Xi(1530)0 differnt sign", kTHnSparseF, {centAxis, ptAxis, invMassAxis, flagAxis});
      histos.add("h3XiinvmassDSAnti", "Invariant mass of Anti-Xi- differnt sign", kTHnSparseF, {centAxis, ptAxis, invMassAxisCasc, flagAxis});

      histos.add("h3Xi1530invmassLSAnti", "Invariant mass of Anti-Xi(1530)0 same sign", kTHnSparseF, {centAxis, ptAxis, invMassAxis, flagAxis});
      histos.add("h3XiinvmassLSAnti", "Invariant mass of Anti-Xi- same sign", kTHnSparseF, {centAxis, ptAxis, invMassAxisCasc, flagAxis});
    }

    if (doprocessMEDF || doprocessMEMicro) {
      histos.add("h3Xi1530invmassME_DS", "Invariant mass of Xi(1530)0 mixed event DS", kTHnSparseF, {centAxis, ptAxis, invMassAxis, flagAxis});
      histos.add("h3XiinvmassME_DS", "Invariant mass of Xi- mixed event DS", kTHnSparseF, {centAxis, ptAxis, invMassAxisCasc, flagAxis});

      histos.add("h3Xi1530invmassME_DSAnti", "Invariant mass of Xi(1530)0 mixed event DSAnti", kTHnSparseF, {centAxis, ptAxis, invMassAxis, flagAxis});
      histos.add("h3XiinvmassME_DSAnti", "Invariant mass of Xi- mixed event DSAnti", kTHnSparseF, {centAxis, ptAxis, invMassAxisCasc, flagAxis});
    }

    if (doprocessMC) {
      // MC QA
      histos.add("QAMCTrue/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH2F, {ptAxis, dcaxyAxis});
      histos.add("QAMCTrue/trkDCAxy_xi", "DCAxy distribution of Xi- track candidates", HistType::kTH2F, {ptAxis, dcaxyAxis});

      histos.add("QAMCTrue/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH2F, {ptAxis, dcazAxis});
      histos.add("QAMCTrue/trkDCAz_xi", "DCAz distribution of Xi- track candidates", HistType::kTH2F, {ptAxis, dcazAxis});

      if (pidPlots) {
        histos.add("QAMCTrue/TOF_TPC_Map_pi_first_all", "TOF + TPC Combined PID for Pion_{First};#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
        histos.add("QAMCTrue/TOF_Nsigma_pi_first_all", "TOF NSigma for Pion_{First};#it{p}_{T} (GeV/#it{c});#sigma_{TOF}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
        histos.add("QAMCTrue/TPC_Nsigma_pi_first_all", "TPC NSigma for Pion_{First};#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

        histos.add("QAMCTrue/TOF_TPC_Map_pi_bachelor_all", "TOF + TPC Combined PID for Pion_{Bachelor};#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
        histos.add("QAMCTrue/TPC_Nsigma_pi_bachelor_all", "TPC NSigma for Pion_{Bachelor};#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

        histos.add("QAMCTrue/TOF_TPC_Map_pr_all", "TOF + TPC Combined PID for Proton;#sigma_{TOF}^{Proton};#sigma_{TPC}^{Proton}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
        histos.add("QAMCTrue/TOF_TPC_Map_antipr_all", "TOF + TPC Combined PID for Anti-Proton;#sigma_{TOF}^{Proton};#sigma_{TPC}^{Proton}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});

        histos.add("QAMCTrue/TPC_Nsigma_pr_all", "TPC NSigma for Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
        histos.add("QAMCTrue/TPC_Nsigma_antipr_all", "TPC NSigma for Anti-Proton;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Proton};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});

        histos.add("QAMCTrue/TOF_TPC_Map_pi_all", "TOF + TPC Combined PID for Pion;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});
        histos.add("QAMCTrue/TOF_TPC_Map_piminus_all", "TOF + TPC Combined PID for Pion -;#sigma_{TOF}^{Pion};#sigma_{TPC}^{Pion}", {HistType::kTH2F, {pidQAAxis, pidQAAxis}});

        histos.add("QAMCTrue/TPC_Nsigma_pi_all", "TPC NSigma for Pion;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
        histos.add("QAMCTrue/TPC_Nsigma_piminus_all", "TPC NSigma for Pion -;#it{p}_{T} (GeV/#it{c});#sigma_{TPC}^{Pion};", {HistType::kTH3F, {centAxis, ptAxisQA, pidQAAxis}});
      }

      histos.add("h3RecXi1530invmass", "Invariant mass of Reconstructed MC Xi(1530)0", kTHnSparseF, {centAxis, ptAxis, invMassAxis, flagAxis});
      histos.add("h3RecXiinvmass", "Invariant mass of Reconstructed MC Xi-", kTHnSparseF, {centAxis, ptAxis, invMassAxisCasc, flagAxis});

      histos.add("h3RecXi1530invmassAnti", "Invariant mass of Reconstructed MC Anti-Xi(1530)0", kTHnSparseF, {centAxis, ptAxis, invMassAxis, flagAxis});
      histos.add("h3RecXiinvmassAnti", "Invariant mass of Reconstructed MC Anti-Xi-", kTHnSparseF, {centAxis, ptAxis, invMassAxisCasc, flagAxis});

      histos.add("h3Xi1530Gen", "pT distribution of True MC Xi(1530)0", kTHnSparseF, {mcLabelAxis, ptAxis, centAxis});
      histos.add("h3Xi1530GenAnti", "pT distribution of True MC Anti-Xi(1530)0", kTHnSparseF, {mcLabelAxis, ptAxis, centAxis});

      histos.add("Xi1530Rec", "pT distribution of Reconstructed MC Xi(1530)0", kTH2F, {ptAxis, centAxis});
      histos.add("Xi1530RecAnti", "pT distribution of Reconstructed MC Anti-Xi(1530)0", kTH2F, {ptAxis, centAxis});
      histos.add("Xi1530Recinvmass", "Inv mass distribution of Reconstructed MC Xi(1530)0", kTH1F, {invMassAxis});
    }

    if (additionalQAplots) {
      histos.add("QAbefore/V0sDCADoughter_aspt", "V0s DCA Doughter distribution as pt", HistType::kTH2F, {ptAxis, dcaDaugAxis});
      histos.add("QAbefore/CascDCADoughter_aspt", "Casc DCA Doughter distribution as pt", HistType::kTH2F, {ptAxis, dcaDaugAxis});
      histos.add("QAbefore/CascMass_aspt", "Casc DCA Bachlor distribution as pt", HistType::kTH2F, {ptAxis, invMassAxisCasc});
      histos.add("QAbefore/V0sCosPA_aspt", "V0s CosPA distribution as pt", HistType::kTH2F, {ptAxis, cosPAAxis});
      histos.add("QAbefore/CascCosPA_aspt", "Casc CosPA distribution as pt", HistType::kTH2F, {ptAxis, cosPAAxis});

      histos.add("QAafter/V0sDCADoughter_aspt", "V0s DCA Doughter distribution as pt", HistType::kTH2F, {ptAxis, dcaDaugAxis});
      histos.add("QAafter/CascDCADoughter_aspt", "Casc DCA Doughter distribution as pt", HistType::kTH2F, {ptAxis, dcaDaugAxis});
      histos.add("QAafter/CascMass_aspt", "Casc DCA Bachlor distribution as pt", HistType::kTH2F, {ptAxis, invMassAxisCasc});
      histos.add("QAafter/V0sCosPA_aspt", "V0s CosPA distribution as pt", HistType::kTH2F, {ptAxis, cosPAAxis});
      histos.add("QAafter/CascCosPA_aspt", "Casc CosPA distribution as pt", HistType::kTH2F, {ptAxis, cosPAAxis});

      histos.add("QAMCTrue/V0sDCADoughter_aspt", "V0s DCA Doughter distribution as pt", HistType::kTH2F, {ptAxis, dcaDaugAxis});
      histos.add("QAMCTrue/CascDCADoughter_aspt", "Casc DCA Doughter distribution as pt", HistType::kTH2F, {ptAxis, dcaDaugAxis});
      histos.add("QAMCTrue/CascMass_aspt", "Casc DCA Bachlor distribution as pt", HistType::kTH2F, {ptAxis, invMassAxisCasc});
      histos.add("QAMCTrue/V0sCosPA_aspt", "V0s CosPA distribution as pt", HistType::kTH2F, {ptAxis, cosPAAxis});
      histos.add("QAMCTrue/CascCosPA_aspt", "Casc CosPA distribution as pt", HistType::kTH2F, {ptAxis, cosPAAxis});
    }
  }

  double massPi = MassPionCharged;

  // Primary track selection for the first pion //
  template <bool IsResoMicrotrack, typename TrackType>
  bool primaryTrackCut(const TrackType track)
  {
    if (std::abs(track.eta()) > cMaxetacut)
      return false;
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if constexpr (IsResoMicrotrack) {
      if (std::abs(o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(track.trackSelectionFlags())) > (cDCAxytoPVByPtPiFirstP0 + cDCAxyToPVByPtPiFirstExp * std::pow(track.pt(), -1.1)))
        return false;
      if (cDCAzToPVAsPt) {
        if (std::abs(o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(track.trackSelectionFlags())) > (cDCAxytoPVByPtPiFirstP0 + cDCAxyToPVByPtPiFirstExp * std::pow(track.pt(), -1.1)))
          return false;
      } else {
        if (std::abs(o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(track.trackSelectionFlags())) > cMaxDCAzToPVCut)
          return false;
      }
    } else {
      if (std::abs(track.dcaXY()) > (cDCAxytoPVByPtPiFirstP0 + cDCAxyToPVByPtPiFirstExp * std::pow(track.pt(), -1.1)))
        return false;
      if (cDCAzToPVAsPt) {
        if (std::abs(track.dcaZ()) > (cDCAxytoPVByPtPiFirstP0 + cDCAxyToPVByPtPiFirstExp * std::pow(track.pt(), -1.1)))
          return false;
      } else {
        if (std::abs(track.dcaZ()) > cMaxDCAzToPVCut)
          return false;
      }
      if (track.tpcNClsFound() < cfgTPCcluster)
        return false;
      if (track.tpcNClsCrossedRows() < cfgTPCRows)
        return false;
    }
    if (cfgHasTOF && !track.hasTOF())
      return false;
    // if (cfgUseITSRefit && !track.passedITSRefit())
    //   return false;
    // if (cfgUseTPCRefit && !track.passedTPCRefit())
    //   return false;
    if (cfgPVContributor && !track.isPVContributor())
      return false;
    if (cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    return true;
  }

  bool hasSubsystemInfo(float Nsigma)
  {
    return std::abs(Nsigma) < cPIDBound;
  }

  // Primary track selection for cascades, Need to more informations for cascades //
  template <typename TracksTypeCasc>
  bool cascprimaryTrackCut(const TracksTypeCasc track)
  {
    if (std::abs(track.eta()) > cMaxetacut)
      return false;
    if (std::abs(track.pt()) < cMinPtcut)
      return false;
    if (cDCAxyToPVAsPtForCasc) {
      if (std::abs(track.dcaXYCascToPV()) > (cDCAxyToPVByPtCascP0 + cDCAxyToPVByPtCascExp * track.pt()))
        return false;
    }
    if (cDCAzToPVAsPtForCasc) {
      if (std::abs(track.dcaZCascToPV()) > (cDCAxyToPVByPtCascP0 + cDCAxyToPVByPtCascExp * std::pow(track.pt(), -1.1)))
        return false;
    }

    return true;
  }

  // Secondary track selection for cascades //

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
    } else {
      if (std::abs(track.dcanegtopv()) < cDCAProtonToPVcut)
        return false;
      if (std::abs(track.dcapostopv()) < cDCAPionToPVcut)
        return false;
    }
    if (track.v0CosPA() < std::cos(cV0CosPACutPtDepP0 - cV0CosPACutPtDepP1 * track.pt()))
      return false;
    if (track.transRadius() > cMaxV0radiuscut || track.transRadius() < cMinV0radiuscut)
      return false;
    if (std::abs(track.mLambda() - MassLambda) > cMasswindowV0cut)
      return false;

    // Topological Cuts for Cascades
    if (std::abs(track.dcabachtopv()) < cDCABachlorToPVcut)
      return false;
    if (track.pt() < cDCAXiDaugthersCutPtRangeLower) {
      if (track.cascDaughDCA() > cDCAXiDaugthersCutPtDepLower)
        return false;
    }
    if (track.pt() >= cDCAXiDaugthersCutPtRangeLower && track.pt() < cDCAXiDaugthersCutPtRangeUpper) {
      if (track.cascDaughDCA() > cDCAXiDaugthersCutPtDepMiddle)
        return false;
    }
    if (track.pt() >= cDCAXiDaugthersCutPtRangeUpper) {
      if (track.cascDaughDCA() > cDCAXiDaugthersCutPtDepUpper)
        return false;
    }
    if (track.cascCosPA() < std::cos(cCosPACascCutPtDepP0 - cCosPACascCutPtDepP1 * track.pt()))
      return false;

    if (track.cascTransRadius() > cMaxCascradiuscut || track.cascTransRadius() < cMinCascradiuscut)
      return false;
    if (std::abs(track.mXi() - cMassXiminus) > cMasswindowCasccut)
      return false;

    return true;
  }

  bool pidSelector(float TPCNsigma, float TOFNsigma, const PidSelectionParam& params, bool tofAtHighPt, float trackPt)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};

    if (tofAtHighPt && trackPt > cMinTOFpt) {
      if (std::abs(TPCNsigma) < params.cMaxTPCnSigma) {
        tpcPIDPassed = true;
      }

      if (params.cByPassTOF && tpcPIDPassed) {
        return true;
      }

      if (hasSubsystemInfo(TOFNsigma)) {
        if (std::abs(TOFNsigma) < params.cMaxTOFnSigma) {
          tofPIDPassed = true;
        }
        if ((params.nsigmaCutCombined > 0) &&
            (TPCNsigma * TPCNsigma + TOFNsigma * TOFNsigma < params.nsigmaCutCombined * params.nsigmaCutCombined)) {
          tofPIDPassed = true;
        }
      } else {
        tofPIDPassed = true;
      }
      return tpcPIDPassed && tofPIDPassed;
    } else {

      if (std::abs(TPCNsigma) < params.cMaxTPCnSigma) {
        tpcPIDPassed = true;
      }

      if (params.cByPassTOF && tpcPIDPassed) {
        return true;
      }

      if (hasSubsystemInfo(TOFNsigma)) {
        if (std::abs(TOFNsigma) < params.cMaxTOFnSigma) {
          tofPIDPassed = true;
        }
        if ((params.nsigmaCutCombined > 0) &&
            (TPCNsigma * TPCNsigma + TOFNsigma * TOFNsigma < params.nsigmaCutCombined * params.nsigmaCutCombined)) {
          tofPIDPassed = true;
        }
      } else {
        tofPIDPassed = true;
      }

      return tpcPIDPassed && tofPIDPassed;
    }
  }

  // PID selection for the First Pion //
  template <bool IsResoMicrotrack, typename T>
  bool selectionPIDPionFirst(const T& candidate)
  {

    float tpcNsigmaPionFirst, tofNsigmaPionFirst;
    float trackPt = candidate.pt();

    if constexpr (IsResoMicrotrack) {
      tpcNsigmaPionFirst = o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(candidate.pidNSigmaPiFlag());
      tofNsigmaPionFirst = o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(candidate.pidNSigmaPiFlag());
    } else {
      tpcNsigmaPionFirst = candidate.tpcNSigmaPi();
      tofNsigmaPionFirst = candidate.tofNSigmaPi();
    }

    PidSelectionParam pionFirstParams = {cMaxtpcnSigmaPionFirst, cMaxtofnSigmaPionFirst, cByPassTOFPionFirst, nsigmaCutCombinedPionFirst};

    return pidSelector(tpcNsigmaPionFirst, tofNsigmaPionFirst, pionFirstParams, tofAtHighPt, trackPt);
  }

  template <typename TCascade>
  bool selectionPIDCascades(const TCascade& candidate)
  {
    bool lConsistentWithXi{false}, lConsistentWithLambda{false}, lConsistentWithPion{false}, lConsistentWithProton{false};

    float tpcNsigmaBachelor, tofNsigmaBachelor;
    float tpcNsigmaPion, tofNsigmaPion;
    float tpcNsigmaProton, tofNsigmaProton;
    float trackPt = candidate.pt();

    if (candidate.sign() < 0) { // Xi- candidates
      tpcNsigmaBachelor = candidate.daughterTPCNSigmaBachPi();
      tofNsigmaBachelor = candidate.daughterTOFNSigmaBachPi();

      tpcNsigmaPion = candidate.daughterTPCNSigmaNegPi();
      tofNsigmaPion = candidate.daughterTOFNSigmaNegPi();

      tpcNsigmaProton = candidate.daughterTPCNSigmaPosPr();
      tofNsigmaProton = candidate.daughterTOFNSigmaPosPr();
    } else { // Anti-Xi- candidates

      tpcNsigmaBachelor = candidate.daughterTPCNSigmaBachPi();
      tofNsigmaBachelor = candidate.daughterTOFNSigmaBachPi();

      tpcNsigmaPion = candidate.daughterTPCNSigmaPosPi();
      tofNsigmaPion = candidate.daughterTOFNSigmaPosPi();

      tpcNsigmaProton = candidate.daughterTPCNSigmaNegPr();
      tofNsigmaProton = candidate.daughterTOFNSigmaNegPr();
    }

    PidSelectionParam bachelorParams = {cMaxtpcnSigmaPionBachelor, cMaxtofnSigmaPionBachelor, cByPassTOFPionBachelor, nsigmaCutCombinedPionBachelor};
    PidSelectionParam pionParams = {cMaxtpcnSigmaPion, cMaxtofnSigmaPion, cByPassTOFPion, nsigmaCutCombinedPion};
    PidSelectionParam protonParams = {cMaxtpcnSigmaProton, cMaxtofnSigmaProton, cByPassTOFProton, nsigmaCutCombinedProton};

    lConsistentWithXi = pidSelector(tpcNsigmaBachelor, tofNsigmaBachelor, bachelorParams, tofAtHighPt, trackPt);
    lConsistentWithPion = pidSelector(tpcNsigmaPion, tofNsigmaPion, pionParams, tofAtHighPt, trackPt);
    lConsistentWithProton = pidSelector(tpcNsigmaProton, tofNsigmaProton, protonParams, tofAtHighPt, trackPt);

    lConsistentWithLambda = lConsistentWithProton && lConsistentWithPion;

    return lConsistentWithXi && lConsistentWithLambda;
  }

  template <bool IsResoMicrotrack, bool IsMC, bool IsMix, typename CollisionType, typename CenMult, typename TracksType, typename TracksTypeCasc>
  void fillHistograms(const CollisionType& collision, const CenMult& multiplicity, const TracksType& dTracks1, const TracksTypeCasc& dTracks2) // Order: ResoColl, ResoTrack, ResoCascTrack
  {
    // auto multiplicity = collision.cent();

    {
      if constexpr (!IsMix) {
        histos.fill(HIST("QAevent/hVertexZSameE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentSameE"), multiplicity);
        histos.fill(HIST("TestME/hCollisionIndexSameE"), collision.globalIndex());
        histos.fill(HIST("TestME/hnTrksSameE"), dTracks1.size());
      } else {
        histos.fill(HIST("QAevent/hVertexZMixedE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentMixedE"), multiplicity);
        histos.fill(HIST("TestME/hCollisionIndexMixedE"), collision.globalIndex());
        histos.fill(HIST("TestME/hnTrksMixedE"), dTracks1.size());
      }
    }

    LorentzVectorPtEtaPhiMass lDecayDaughter1, lDecayDaughter2, lResonance; // It will be replaced to use RecoDecay (In fixing...)

    for (const auto& [trk1, trk2] : combinations(CombinationsFullIndexPolicy(dTracks1, dTracks2))) {

      {
        if constexpr (!IsMix) {
          histos.fill(HIST("TestME/hPairsCounterSameE"), 1.0);
        } else {
          histos.fill(HIST("TestME/hPairsCounterMixedE"), 1.0);
        }
      }

      auto trk1ptPi = trk1.pt();
      static float trk1DCAXY;
      static float trk1DCAZ;
      static float trk1NSigmaPiTPC;
      static float trk1NSigmaPiTOF;
      if constexpr (IsResoMicrotrack) {
        trk1DCAXY = o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(trk1.trackSelectionFlags());
        trk1DCAZ = o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(trk1.trackSelectionFlags());
        trk1NSigmaPiTPC = o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(trk1.pidNSigmaPiFlag());
        trk1NSigmaPiTOF = o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(trk1.pidNSigmaPiFlag());
      } else {
        trk1DCAXY = trk1.dcaXY();
        trk1DCAZ = trk1.dcaZ();
        trk1NSigmaPiTPC = trk1.tpcNSigmaPi();
        trk1NSigmaPiTOF = trk1.tofNSigmaPi();
      }

      auto trk2ptXi = trk2.pt();

      auto trk2DCAXY = trk2.dcaXYCascToPV();
      auto trk2DCAZ = trk2.dcaZCascToPV();

      auto trk2DCAV0sDougthers = trk2.daughDCA();
      auto trk2DCACascDougthers = trk2.cascDaughDCA();
      auto trk2Mass = trk2.mXi();
      auto trk2CascCosPA = trk2.cascCosPA();
      auto trk2V0sCosPA = trk2.v0CosPA();

      // Need to daughther's pt info. in the table
      // auto trk2ptPiBachelor = trk2.pt();
      float trk2NSigmaPiBachelorTPC = trk2.daughterTPCNSigmaBachPi();
      float trk2NSigmaPiBachelorTOF = trk2.daughterTOFNSigmaBachPi();

      // auto trk2ptPr = trk2.pt();
      float trk2NSigmaPrPosTPC = trk2.daughterTPCNSigmaPosPr();
      float trk2NSigmaPrNegTPC = trk2.daughterTPCNSigmaNegPr();

      float trk2NSigmaPrPosTOF = trk2.daughterTOFNSigmaPosPr();
      float trk2NSigmaPrNegTOF = trk2.daughterTOFNSigmaNegPr();

      // auto trk2ptPi = trk2.pt();
      float trk2NSigmaPiPosTPC = trk2.daughterTPCNSigmaPosPi();
      float trk2NSigmaPiNegTPC = trk2.daughterTPCNSigmaNegPi();

      float trk2NSigmaPiPosTOF = trk2.daughterTOFNSigmaPosPi();
      float trk2NSigmaPiNegTOF = trk2.daughterTOFNSigmaNegPi();

      if constexpr (!IsMix) {
        //// QA plots before the selection //
        //  --- PID QA
        if (pidPlots) {
          histos.fill(HIST("QAbefore/TPC_Nsigma_pi_first_all"), multiplicity, trk1ptPi, trk1NSigmaPiTPC);
          if (hasSubsystemInfo(trk1NSigmaPiTOF)) {
            histos.fill(HIST("QAbefore/TOF_Nsigma_pi_first_all"), multiplicity, trk1ptPi, trk1NSigmaPiTOF);
            histos.fill(HIST("QAbefore/TOF_TPC_Map_pi_first_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
          }

          histos.fill(HIST("QAbefore/TPC_Nsigma_pi_bachelor_all"), multiplicity, 0, trk2NSigmaPiBachelorTPC); // can't take pt information for the cascade secondary
          if (hasSubsystemInfo(trk2NSigmaPiBachelorTOF)) {
            histos.fill(HIST("QAbefore/TOF_TPC_Map_pi_bachelor_all"), trk2NSigmaPiBachelorTOF, trk2NSigmaPiBachelorTPC);
          }

          histos.fill(HIST("QAbefore/TPC_Nsigma_pr_all"), multiplicity, 0, trk2NSigmaPrPosTPC);
          if (hasSubsystemInfo(trk2NSigmaPrPosTOF)) {
            histos.fill(HIST("QAbefore/TOF_TPC_Map_pr_all"), trk2NSigmaPrPosTOF, trk2NSigmaPrPosTPC);
          }

          histos.fill(HIST("QAbefore/TPC_Nsigma_antipr_all"), multiplicity, 0, trk2NSigmaPrNegTPC);
          if (hasSubsystemInfo(trk2NSigmaPrNegTOF)) {
            histos.fill(HIST("QAbefore/TOF_TPC_Map_antipr_all"), trk2NSigmaPrNegTOF, trk2NSigmaPrNegTPC);
          }

          histos.fill(HIST("QAbefore/TPC_Nsigma_pi_all"), multiplicity, 0, trk2NSigmaPiPosTPC);
          if (hasSubsystemInfo(trk2NSigmaPiPosTOF)) {
            histos.fill(HIST("QAbefore/TOF_TPC_Map_pi_all"), trk2NSigmaPiPosTOF, trk2NSigmaPiPosTPC);
          }

          histos.fill(HIST("QAbefore/TPC_Nsigma_piminus_all"), multiplicity, 0, trk2NSigmaPiNegTPC);
          if (hasSubsystemInfo(trk2NSigmaPiNegTOF)) {
            histos.fill(HIST("QAbefore/TOF_TPC_Map_piminus_all"), trk2NSigmaPiNegTOF, trk2NSigmaPiNegTPC);
          }
        }

        histos.fill(HIST("QAbefore/trkDCAxy_pi"), trk1ptPi, trk1DCAXY);
        histos.fill(HIST("QAbefore/trkDCAxy_Xi"), trk2ptXi, trk2DCAXY);

        histos.fill(HIST("QAbefore/trkDCAz_pi"), trk1ptPi, trk1DCAZ);
        histos.fill(HIST("QAbefore/trkDCAz_Xi"), trk2ptXi, trk2DCAZ);

        if (additionalQAplots) {
          histos.fill(HIST("QAbefore/V0sDCADoughter_aspt"), trk2ptXi, trk2DCAV0sDougthers);
          histos.fill(HIST("QAbefore/CascDCADoughter_aspt"), trk2ptXi, trk2DCACascDougthers);
          histos.fill(HIST("QAbefore/CascMass_aspt"), trk2ptXi, trk2Mass);
          histos.fill(HIST("QAbefore/V0sCosPA_aspt"), trk2ptXi, 1. - trk2V0sCosPA);
          histos.fill(HIST("QAbefore/CascCosPA_aspt"), trk2ptXi, 1. - trk2CascCosPA);
        }
      }

      if (!primaryTrackCut<IsResoMicrotrack>(trk1) || !cascprimaryTrackCut(trk2)) // Primary track selections
        continue;

      // PID selection
      if (cUseOnlyTOFTrackPionFirst && !hasSubsystemInfo(trk1NSigmaPiTOF))
        continue;

      if (cUseOnlyTOFTrackPionBachelor && !hasSubsystemInfo(trk2NSigmaPiBachelorTOF))
        continue;

      if (cUseOnlyTOFTrackProton && !hasSubsystemInfo(trk2NSigmaPrPosTOF))
        continue;
      if (cUseOnlyTOFTrackProton && !hasSubsystemInfo(trk2NSigmaPrNegTOF))
        continue;

      if (cUseOnlyTOFTrackPion && !hasSubsystemInfo(trk2NSigmaPiPosTOF))
        continue;
      if (cUseOnlyTOFTrackPion && !hasSubsystemInfo(trk2NSigmaPiNegTOF))
        continue;

      if (!selectionPIDPionFirst<IsResoMicrotrack>(trk1) || !selectionPIDCascades(trk2))
        continue;

      if (!casctopCut(trk2))
        continue;

      if constexpr (!IsMix) {
        //// QA plots after the selection
        //  --- PID QA
        if (pidPlots) {
          histos.fill(HIST("QAafter/TPC_Nsigma_pi_first_all"), multiplicity, trk1ptPi, trk1NSigmaPiTPC);

          if (hasSubsystemInfo(trk1NSigmaPiTOF)) {
            histos.fill(HIST("QAafter/TOF_Nsigma_pi_first_all"), multiplicity, trk1ptPi, trk1NSigmaPiTOF);
            histos.fill(HIST("QAafter/TOF_TPC_Map_pi_first_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
          }

          if (trk2.sign() < 0) {
            histos.fill(HIST("QAafter/TPC_Nsigma_pi_bachelor_all"), multiplicity, 0, trk2NSigmaPiBachelorTPC); // not exist pt information in resocascade yet.
            if (hasSubsystemInfo(trk2NSigmaPiBachelorTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_pi_bachelor_all"), trk2NSigmaPiBachelorTOF, trk2NSigmaPiBachelorTPC);
            }

            histos.fill(HIST("QAafter/TPC_Nsigma_pr_all"), multiplicity, 0, trk2NSigmaPrPosTPC);
            if (hasSubsystemInfo(trk2NSigmaPrPosTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_pr_all"), trk2NSigmaPrPosTOF, trk2NSigmaPrPosTPC);
            }

            histos.fill(HIST("QAafter/TPC_Nsigma_piminus_all"), multiplicity, 0, trk2NSigmaPiNegTPC);
            if (hasSubsystemInfo(trk2NSigmaPiNegTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_piminus_all"), trk2NSigmaPiNegTOF, trk2NSigmaPiNegTPC);
            }

          } else {

            histos.fill(HIST("QAafter/TPC_Nsigma_pi_bachelor_all"), multiplicity, 0, trk2NSigmaPiBachelorTPC); // not exist pt information in resocascade yet.
            if (hasSubsystemInfo(trk2NSigmaPiBachelorTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_pi_bachelor_all"), trk2NSigmaPiBachelorTOF, trk2NSigmaPiBachelorTPC);
            }

            histos.fill(HIST("QAafter/TPC_Nsigma_antipr_all"), multiplicity, 0, trk2NSigmaPrNegTPC);
            if (hasSubsystemInfo(trk2NSigmaPrNegTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_antipr_all"), trk2NSigmaPrNegTOF, trk2NSigmaPrNegTPC);
            }

            histos.fill(HIST("QAafter/TPC_Nsigma_pi_all"), multiplicity, 0, trk2NSigmaPiPosTPC);
            if (hasSubsystemInfo(trk2NSigmaPiPosTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_pi_all"), trk2NSigmaPiPosTOF, trk2NSigmaPiPosTPC);
            }
          }
        }

        histos.fill(HIST("QAafter/trkDCAxy_pi"), trk1ptPi, trk1DCAXY);
        histos.fill(HIST("QAafter/trkDCAxy_Xi"), trk2ptXi, trk2DCAXY);

        histos.fill(HIST("QAafter/trkDCAz_pi"), trk1ptPi, trk1DCAZ);
        histos.fill(HIST("QAafter/trkDCAz_Xi"), trk2ptXi, trk2DCAZ);

        if (additionalQAplots) {
          histos.fill(HIST("QAafter/V0sDCADoughter_aspt"), trk2ptXi, trk2DCAV0sDougthers);
          histos.fill(HIST("QAafter/CascDCADoughter_aspt"), trk2ptXi, trk2DCACascDougthers);
          histos.fill(HIST("QAafter/CascMass_aspt"), trk2ptXi, trk2Mass);
          histos.fill(HIST("QAafter/V0sCosPA_aspt"), trk2ptXi, 1. - trk2V0sCosPA);
          histos.fill(HIST("QAafter/CascCosPA_aspt"), trk2ptXi, 1. - trk2CascCosPA);
        }
      }

      lDecayDaughter1 = LorentzVectorPtEtaPhiMass(trk1ptPi, trk1.eta(), trk1.phi(), massPi);
      lDecayDaughter2 = LorentzVectorPtEtaPhiMass(trk2ptXi, trk2.eta(), trk2.phi(), trk2.mXi());
      lResonance = lDecayDaughter1 + lDecayDaughter2;

      auto lResonancePt = lResonance.Pt();
      auto lResonanceInMass = lResonance.M();

      if (std::abs(lResonance.Rapidity()) >= cfgRapidityCut)
        continue;

      if (cfgCutsOnMother) {
        if (lResonancePt >= cMaxPtMotherCut) // excluding candidates in overflow
          continue;
        if (lResonanceInMass >= cMaxMinvMotherCut) // excluding candidates in overflow
          continue;
      }

      if (trk1.sign() * trk2.sign() < 0) {

        if constexpr (!IsMix) {
          if (studyAntiparticle) {
            if (trk1.sign() > 0) {
              if (invMass1D)
                histos.fill(HIST("Xi1530invmassDS"), lResonanceInMass);
              histos.fill(HIST("h3Xi1530invmassDS"), multiplicity, lResonancePt, lResonanceInMass, kData);
            } else if (trk1.sign() < 0) {
              if (invMass1D)
                histos.fill(HIST("Xi1530invmassDSAnti"), lResonanceInMass);
              histos.fill(HIST("h3Xi1530invmassDSAnti"), multiplicity, lResonancePt, lResonanceInMass, kData);
            }
          } else {
            if (invMass1D)
              histos.fill(HIST("Xi1530invmassDS"), lResonanceInMass);
            histos.fill(HIST("h3Xi1530invmassDS"), multiplicity, lResonancePt, lResonanceInMass, kData);
          }

          if (studyStableXi) {
            if (trk1.sign() > 0) {
              histos.fill(HIST("h3XiinvmassDS"), multiplicity, trk2ptXi, trk2Mass, kData);
            } else if (trk1.sign() < 0) {
              histos.fill(HIST("h3XiinvmassDSAnti"), multiplicity, trk2ptXi, trk2Mass, kData);
            }
          }
        } else {
          if (invMass1D)
            histos.fill(HIST("Xi1530invmassME"), lResonanceInMass);
          if (trk1.sign() > 0) {
            if (invMass1D)
              histos.fill(HIST("Xi1530invmassME_DS"), lResonanceInMass);
            histos.fill(HIST("h3Xi1530invmassME_DS"), multiplicity, lResonancePt, lResonanceInMass, kMixing);
          } else if (trk1.sign() < 0) {
            if (invMass1D)
              histos.fill(HIST("Xi1530invmassME_DSAnti"), lResonanceInMass);
            histos.fill(HIST("h3Xi1530invmassME_DSAnti"), multiplicity, lResonancePt, lResonanceInMass, kMixing);
          }
          histos.fill(HIST("h3Xi1530invmassME"), multiplicity, lResonancePt, lResonanceInMass, kMixing);

          if (studyStableXi) {
            if (trk1.sign() > 0) {
              histos.fill(HIST("h3XiinvmassME_DS"), multiplicity, trk2ptXi, trk2Mass, kMixing);
            } else if (trk1.sign() < 0) {
              histos.fill(HIST("h3XiinvmassME_DSAnti"), multiplicity, trk2ptXi, trk2Mass, kMixing);
            }
          }
        }
        if constexpr (IsMC) {
          if (std::abs(trk2.motherPDG()) != kXiStar)
            continue;
          if (std::abs(trk1.pdgCode()) != kPiPlus || std::abs(trk2.pdgCode()) != kXiMinus)
            continue;
          if (trk1.motherId() != trk2.motherId())
            continue;

          histos.fill(HIST("QAMCTrue/trkDCAxy_pi"), trk1ptPi, trk1DCAXY);
          histos.fill(HIST("QAMCTrue/trkDCAxy_xi"), trk2ptXi, trk2DCAXY);
          histos.fill(HIST("QAMCTrue/trkDCAz_pi"), trk1ptPi, trk1DCAZ);
          histos.fill(HIST("QAMCTrue/trkDCAz_xi"), trk2ptXi, trk2DCAZ);

          if (additionalQAplots) {
            histos.fill(HIST("QAMCTrue/V0sDCADoughter_aspt"), trk2ptXi, trk2DCAV0sDougthers);
            histos.fill(HIST("QAMCTrue/CascDCADoughter_aspt"), trk2ptXi, trk2DCACascDougthers);
            histos.fill(HIST("QAMCTrue/CascMass_aspt"), trk2ptXi, trk2Mass);
            histos.fill(HIST("QAMCTrue/V0sCosPA_aspt"), trk2ptXi, 1. - trk2V0sCosPA);
            histos.fill(HIST("QAMCTrue/CascCosPA_aspt"), trk2ptXi, 1. - trk2CascCosPA);
          }

          histos.fill(HIST("QAMCTrue/TPC_Nsigma_pi_first_all"), multiplicity, trk1ptPi, trk1NSigmaPiTPC);
          if (hasSubsystemInfo(trk1NSigmaPiTOF)) {
            histos.fill(HIST("QAMCTrue/TOF_Nsigma_pi_first_all"), multiplicity, trk1ptPi, trk1NSigmaPiTOF);
            histos.fill(HIST("QAMCTrue/TOF_TPC_Map_pi_first_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
          }

          if (trk2.sign() < 0) {
            histos.fill(HIST("QAMCTrue/TPC_Nsigma_pi_bachelor_all"), multiplicity, 0, trk2NSigmaPiBachelorTPC); // not exist pt information in resocascade yet.
            if (hasSubsystemInfo(trk2NSigmaPiBachelorTOF)) {
              histos.fill(HIST("QAMCTrue/TOF_TPC_Map_pi_bachelor_all"), trk2NSigmaPiBachelorTOF, trk2NSigmaPiBachelorTPC);
            }

            histos.fill(HIST("QAMCTrue/TPC_Nsigma_pr_all"), multiplicity, 0, trk2NSigmaPrPosTPC);
            if (hasSubsystemInfo(trk2NSigmaPrPosTOF)) {
              histos.fill(HIST("QAMCTrue/TOF_TPC_Map_pr_all"), trk2NSigmaPrPosTOF, trk2NSigmaPrPosTPC);
            }

            histos.fill(HIST("QAMCTrue/TPC_Nsigma_piminus_all"), multiplicity, 0, trk2NSigmaPiNegTPC);
            if (hasSubsystemInfo(trk2NSigmaPiNegTOF)) {
              histos.fill(HIST("QAMCTrue/TOF_TPC_Map_piminus_all"), trk2NSigmaPiNegTOF, trk2NSigmaPiNegTPC);
            }

          } else {

            histos.fill(HIST("QAMCTrue/TPC_Nsigma_pi_bachelor_all"), multiplicity, 0, trk2NSigmaPiBachelorTPC); // not exist pt information in resocascade yet.
            if (hasSubsystemInfo(trk2NSigmaPiBachelorTOF)) {
              histos.fill(HIST("QAMCTrue/TOF_TPC_Map_pi_bachelor_all"), trk2NSigmaPiBachelorTOF, trk2NSigmaPiBachelorTPC);
            }

            histos.fill(HIST("QAMCTrue/TPC_Nsigma_antipr_all"), multiplicity, 0, trk2NSigmaPrNegTPC);
            if (hasSubsystemInfo(trk2NSigmaPrNegTOF)) {
              histos.fill(HIST("QAMCTrue/TOF_TPC_Map_antipr_all"), trk2NSigmaPrNegTOF, trk2NSigmaPrNegTPC);
            }

            histos.fill(HIST("QAMCTrue/TPC_Nsigma_pi_all"), multiplicity, 0, trk2NSigmaPiPosTPC);
            if (hasSubsystemInfo(trk2NSigmaPiPosTOF)) {
              histos.fill(HIST("QAMCTrue/TOF_TPC_Map_pi_all"), trk2NSigmaPiPosTOF, trk2NSigmaPiPosTPC);
            }
          }
          // MC histograms
          if (trk2.motherPDG() > 0) {
            histos.fill(HIST("Xi1530Rec"), lResonancePt, multiplicity);
            histos.fill(HIST("Xi1530Recinvmass"), lResonanceInMass);
            histos.fill(HIST("h3RecXi1530invmass"), multiplicity, lResonancePt, lResonanceInMass, kMCReco);
            histos.fill(HIST("h3RecXiinvmass"), multiplicity, trk2ptXi, trk2Mass, kMCReco);
          } else {
            histos.fill(HIST("Xi1530RecAnti"), lResonancePt, multiplicity);
            histos.fill(HIST("Xi1530Recinvmass"), lResonanceInMass);
            histos.fill(HIST("h3RecXi1530invmassAnti"), multiplicity, lResonancePt, lResonanceInMass, kMCReco);
            histos.fill(HIST("h3RecXiinvmassAnti"), multiplicity, trk2ptXi, trk2Mass, kMCReco);
          }
        }

      } else {
        if constexpr (!IsMix) {
          if (studyAntiparticle) {
            if (trk1.sign() < 0) {
              if (invMass1D)
                histos.fill(HIST("Xi1530invmassLS"), lResonanceInMass);
              histos.fill(HIST("h3Xi1530invmassLS"), multiplicity, lResonancePt, lResonanceInMass, kLS);
              histos.fill(HIST("h3XiinvmassLS"), multiplicity, trk2ptXi, trk2Mass, kLS);
            } else if (trk1.sign() > 0) {
              if (invMass1D)
                histos.fill(HIST("Xi1530invmassLSAnti"), lResonanceInMass);
              histos.fill(HIST("h3Xi1530invmassLSAnti"), multiplicity, lResonancePt, lResonanceInMass, kLS);
              histos.fill(HIST("h3XiinvmassLSAnti"), multiplicity, trk2ptXi, trk2Mass, kLS);
            }
          } else {
            if (invMass1D)
              histos.fill(HIST("Xi1530invmassLS"), lResonanceInMass);
            histos.fill(HIST("h3Xi1530invmassLS"), multiplicity, lResonancePt, lResonanceInMass, kLS);
            histos.fill(HIST("h3XiinvmassLS"), multiplicity, trk2ptXi, trk2Mass, kLS);
          }
        }
      }
    }
  }

  void processData(aod::ResoCollision const& resoCollision,
                   aod::ResoCollisionColls const& collisionIndex,
                   soa::Join<aod::Collisions, aod::PVMults> const& collisions,
                   aod::ResoTracks const& resoTracks,
                   aod::ResoCascades const& cascTracks)
  {
    if (cRecoINELgt0) {
      auto linkRow = collisionIndex.iteratorAt(resoCollision.globalIndex());
      auto collId = linkRow.collisionId(); // Take original collision global index matched with resoCollision

      auto coll = collisions.iteratorAt(collId); // Take original collision matched with resoCollision

      if (!coll.isInelGt0()) // Check reco INELgt0 (at least one PV track in |eta| < 1) about the collision
        return;
    }

    histos.fill(HIST("QAevent/hEvtCounterSameE"), 1.0);
    auto multiplicity = resoCollision.cent();
    fillHistograms<false, false, false>(resoCollision, multiplicity, resoTracks, cascTracks);
  }

  // Reconstructed level MC for the track
  void processMC(ResoMCCols::iterator const& resoCollision,
                 aod::ResoCollisionColls const& collisionIndex,
                 soa::Join<aod::ResoCollisionCandidatesMC, aod::PVMults> const& collisionsMC,
                 soa::Join<aod::ResoCascades, aod::ResoMCCascades> const& cascTracks,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resoTracks,
                 soa::Join<aod::McCollisions, aod::McCentFT0Ms> const&)
  {
    float multiplicity;
    if (cMCCent && cRecoINELgt0) {
      auto linkRow = collisionIndex.iteratorAt(resoCollision.globalIndex());
      auto collId = linkRow.collisionId(); // Take original collision global index matched with resoCollision

      auto coll = collisionsMC.iteratorAt(collId); // Take original collision matched with resoCollision

      if (!coll.isInelGt0()) // Check reco INELgt0 (at least one PV track in |eta| < 1) about the collision
        return;

      auto mcColl = coll.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>();
      multiplicity = mcColl.centFT0M();
    } else if (!cMCCent && cRecoINELgt0) {
      auto linkRow = collisionIndex.iteratorAt(resoCollision.globalIndex());
      auto collId = linkRow.collisionId(); // Take original collision global index matched with resoCollision

      auto coll = collisionsMC.iteratorAt(collId); // Take original collision matched with resoCollision

      if (!coll.isInelGt0()) // Check reco INELgt0 (at least one PV track in |eta| < 1) about the collision
        return;

      multiplicity = resoCollision.cent();
    } else if (cMCCent && !cRecoINELgt0) {
      auto linkRow = collisionIndex.iteratorAt(resoCollision.globalIndex());
      auto collId = linkRow.collisionId(); // Take original collision global index matched with resoCollision

      auto coll = collisionsMC.iteratorAt(collId); // Take original collision matched with resoCollision

      auto mcColl = coll.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>();
      multiplicity = mcColl.centFT0M();
    } else {
      multiplicity = resoCollision.cent();
    }

    if (!resoCollision.isInAfterAllCuts() || (std::abs(resoCollision.posZ()) > cZvertCutMC)) // MC event selection, all cuts missing vtx cut
      return;

    fillHistograms<false, true, false>(resoCollision, multiplicity, resoTracks, cascTracks);
  }

  // Truth level MC for the track with reco event
  void processMCTrue(ResoMCCols::iterator const& resoCollision,
                     aod::ResoCollisionColls const& collisionIndex,
                     aod::ResoMCParents const& resoParents,
                     aod::ResoCollisionCandidatesMC const& collisionsMC,
                     soa::Join<aod::McCollisions, aod::McCentFT0Ms> const&)
  {
    float multiplicity;
    if (cMCCent && cRecoINELgt0) {
      auto linkRow = collisionIndex.iteratorAt(resoCollision.globalIndex());
      auto collId = linkRow.collisionId(); // Take original collision global index matched with resoCollision

      auto coll = collisionsMC.iteratorAt(collId); // Take original collision matched with resoCollision

      if (!coll.isInelGt0()) // Check reco INELgt0 (at least one PV track in |eta| < 1) about the collision
        return;

      auto mcColl = coll.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>();
      multiplicity = mcColl.centFT0M();
    } else if (!cMCCent && cRecoINELgt0) {
      auto linkRow = collisionIndex.iteratorAt(resoCollision.globalIndex());
      auto collId = linkRow.collisionId(); // Take original collision global index matched with resoCollision

      auto coll = collisionsMC.iteratorAt(collId); // Take original collision matched with resoCollision

      if (!coll.isInelGt0()) // Check reco INELgt0 (at least one PV track in |eta| < 1) about the collision
        return;

      multiplicity = resoCollision.cent();
    } else if (cMCCent && !cRecoINELgt0) {
      auto linkRow = collisionIndex.iteratorAt(resoCollision.globalIndex());
      auto collId = linkRow.collisionId(); // Take original collision global index matched with resoCollision

      auto coll = collisionsMC.iteratorAt(collId); // Take original collision matched with resoCollision

      auto mcColl = coll.mcCollision_as<soa::Join<aod::McCollisions, aod::McCentFT0Ms>>();
      multiplicity = mcColl.centFT0M();
    } else {
      multiplicity = resoCollision.cent();
    }

    for (const auto& part : resoParents) { // loop over all pre-filtered MC particles
      if (std::abs(part.pdgCode()) != kXiStar || std::abs(part.y()) >= cfgRapidityCut)
        continue;
      bool pass1 = std::abs(part.daughterPDG1()) == kPiPlus || std::abs(part.daughterPDG2()) == kPiPlus;
      bool pass2 = std::abs(part.daughterPDG1()) == kXiMinus || std::abs(part.daughterPDG2()) == kXiMinus;

      if (!pass1 || !pass2)
        continue;

      if (part.pdgCode() > 0) // INELt0 or INEL
        histos.fill(HIST("h3Xi1530Gen"), -1, part.pt(), multiplicity);
      else
        histos.fill(HIST("h3Xi1530GenAnti"), -1, part.pt(), multiplicity);

      if (resoCollision.isVtxIn10()) // vtx10
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("h3Xi1530Gen"), 0, part.pt(), multiplicity);
        else
          histos.fill(HIST("h3Xi1530GenAnti"), 0, part.pt(), multiplicity);
      }
      if (resoCollision.isVtxIn10() && resoCollision.isInSel8()) // vtx10 && sel8
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("h3Xi1530Gen"), 1, part.pt(), multiplicity);
        else
          histos.fill(HIST("h3Xi1530GenAnti"), 1, part.pt(), multiplicity);
      }
      if (resoCollision.isVtxIn10() && resoCollision.isTriggerTVX()) // vtx10 && TriggerTVX
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("h3Xi1530Gen"), 2, part.pt(), multiplicity);
        else
          histos.fill(HIST("h3Xi1530GenAnti"), 2, part.pt(), multiplicity);
      }
      if (resoCollision.isInAfterAllCuts()) // after all event selection
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("h3Xi1530Gen"), 3, part.pt(), multiplicity);
        else
          histos.fill(HIST("h3Xi1530GenAnti"), 3, part.pt(), multiplicity);
      }
    }
  }

  void processDataMicro(aod::ResoCollision const& resoCollision,
                        aod::ResoCollisionColls const& collisionIndex,
                        soa::Join<aod::Collisions, aod::PVMults> const& collisions,
                        aod::ResoMicroTracks const& resomicrotracks,
                        aod::ResoCascades const& cascTracks)
  {
    if (cRecoINELgt0) {
      auto linkRow = collisionIndex.iteratorAt(resoCollision.globalIndex());
      auto collId = linkRow.collisionId(); // Take original collision global index matched with resoCollision

      auto coll = collisions.iteratorAt(collId); // Take original collision matched with resoCollision

      if (!coll.isInelGt0()) // Check reco INELgt0 (at least one PV track in |eta| < 1) about the collision
        return;
    }

    histos.fill(HIST("QAevent/hEvtCounterSameE"), 1.0);
    auto multiplicity = resoCollision.cent();
    fillHistograms<true, false, false>(resoCollision, multiplicity, resomicrotracks, cascTracks);
  }

  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  Preslice<aod::ResoTrackDFs> perRColdf = aod::resodaughter::resoCollisionDFId;
  Preslice<aod::ResoCascadeDFs> perRColdfCasc = aod::resodaughter::resoCollisionDFId;

  void processMEDF(aod::ResoCollisionDFs const& resoCollisions, aod::ResoTrackDFs const& resotracks, aod::ResoCascadeDFs const& cascTracks)
  {

    auto tracksTuple = std::make_tuple(resotracks, cascTracks);

    BinningTypeVtxZT0M colBinning{{cfgVtxBins, cfgMultBins}, true};
    Pair<aod::ResoCollisionDFs, aod::ResoTrackDFs, aod::ResoCascadeDFs, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, resoCollisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {

      histos.fill(HIST("QAevent/hEvtCounterMixedE"), 1.0);
      auto multiplicity = collision1.cent();
      fillHistograms<false, false, true>(collision1, multiplicity, tracks1, tracks2);
    }
  }
  void processDataDF(aod::ResoCollisionDF const& resoCollision, aod::ResoTrackDFs const& resotracks, aod::ResoCascadeDFs const& cascTracks)
  {
    auto multiplicity = resoCollision.cent();
    fillHistograms<false, false, false>(resoCollision, multiplicity, resotracks, cascTracks);
  }

  void processMEMicro(aod::ResoCollisions const& resoCollisions,
                      aod::ResoCollisionColls const& collisionIndex,
                      soa::Join<aod::Collisions, aod::PVMults> const& collisions,
                      aod::ResoMicroTracks const& resomicrotracks,
                      aod::ResoCascades const& cascTracks)
  {
    auto tracksTuple = std::make_tuple(resomicrotracks, cascTracks);

    BinningTypeVtxZT0M colBinning{{cfgVtxBins, cfgMultBins}, true};
    Pair<aod::ResoCollisions, aod::ResoMicroTracks, aod::ResoCascades, BinningTypeVtxZT0M> pairs{colBinning, nEvtMixing, -1, resoCollisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {

      if (cRecoINELgt0) {
        const auto rcIdx = collision1.globalIndex();
        auto linkRow = collisionIndex.iteratorAt(rcIdx);
        auto collId = linkRow.collisionId(); // Take original collision global index matched with resoCollision

        auto coll = collisions.iteratorAt(collId); // Take original collision matched with resoCollision

        if (!coll.isInelGt0()) // Check reco INELgt0 (at least one PV track in |eta| < 1) about the collision
          continue;
      }
      histos.fill(HIST("QAevent/hEvtCounterMixedE"), 1.0);
      auto multiplicity = collision1.cent();
      fillHistograms<true, false, true>(collision1, multiplicity, tracks1, tracks2);
    }
  }

  PROCESS_SWITCH(Xi1530Analysisqa, processData, "Process Event for Data", false);
  PROCESS_SWITCH(Xi1530Analysisqa, processMC, "Process Event for MC (Reconstructed)", false);
  PROCESS_SWITCH(Xi1530Analysisqa, processMCTrue, "Process Event for MC (Generated)", false);
  PROCESS_SWITCH(Xi1530Analysisqa, processDataMicro, "Process Event for Data (MicroTrack)", false);
  PROCESS_SWITCH(Xi1530Analysisqa, processMEMicro, "Process EventMixing (MicroTrack) ", false);
  PROCESS_SWITCH(Xi1530Analysisqa, processMEDF, "Process EventMixing (DF) ", true);
  PROCESS_SWITCH(Xi1530Analysisqa, processDataDF, "Process Event for Data (DF) ", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Xi1530Analysisqa>(cfgc)};
}
