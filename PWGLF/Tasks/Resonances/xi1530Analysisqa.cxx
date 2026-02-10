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

#include <vector>

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
  // Associated with histograms
  struct : ConfigurableGroup {
    ConfigurableAxis binsPt{"binsPt", {VARIABLE_WIDTH, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9, 4.0, 4.1, 4.2, 4.3, 4.4, 4.5, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.6, 5.7, 5.8, 5.9, 6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 8.0, 8.1, 8.2, 8.3, 8.4, 8.5, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0, 10.1, 10.2, 10.3, 10.4, 10.5, 10.6, 10.7, 10.8, 10.9, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0, 12.1, 12.2, 12.3, 12.4, 12.5, 12.6, 12.7, 12.8, 12.9, 13.0, 13.1, 13.2, 13.3, 13.4, 13.5, 13.6, 13.7, 13.8, 13.9, 14.0, 14.1, 14.2, 14.3, 14.4, 14.5, 14.6, 14.7, 14.8, 14.9, 15.0}, "Binning of the pT axis"};
    ConfigurableAxis binsPtQA{"binsPtQA", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.2, 8.4, 8.6, 8.8, 9.0, 9.2, 9.4, 9.6, 9.8, 10.0}, "Binning of the pT axis"};
    ConfigurableAxis binsCent{"binsCent", {VARIABLE_WIDTH, 0.0, 1.0, 5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0, 101.0, 110.0}, "Binning of the centrality axis"};
    ConfigurableAxis multNTracksAxis{"multNTracksAxis", {500, 0, 500}, "N_{tracks}"};
    ConfigurableAxis impactParameterAxis{"impactParameterAxis", {500, 0, 50}, "IP (fm)"};

    Configurable<float> cInvMassStart{"cInvMassStart", 1.4, "Invariant mass start"};
    Configurable<float> cInvMassEnd{"cInvMassEnd", 3.0, "Invariant mass end"};
    Configurable<int> cInvMassBins{"cInvMassBins", 800, "Invariant mass binning"};

    Configurable<bool> invMass1D{"invMass1D", true, "Invariant mass 1D"};
    Configurable<bool> pidPlots{"pidPlots", true, "Make TPC and TOF PID plots"};
    Configurable<bool> additionalQAplots{"additionalQAplots", true, "Additional QA plots"};
    Configurable<bool> truthQA{"truthQA", true, "Truth particle's QA plots"};
    Configurable<bool> multQA{"multQA", true, "Multiplicity QA after all cuts"};
    Configurable<bool> eventQA{"eventQA", true, "Event QA after all cuts"};
  } histoConfig;

  // Event Mixing
  struct : ConfigurableGroup {
    Configurable<int> nEvtMixing{"nEvtMixing", 10, "Number of events to mix"};
    ConfigurableAxis cfgVtxBins{"cfgVtxBins", {VARIABLE_WIDTH, -10.0f, -8.f, -6.f, -4.f, -2.f, 0.f, 2.f, 4.f, 6.f, 8.f, 10.f}, "Mixing bins - z-vertex"};
    ConfigurableAxis cfgMultBins{"cfgMultBins", {VARIABLE_WIDTH, 0.0f, 1.0f, 5.0f, 10.0f, 20.0f, 30.0f, 40.0f, 50.0f, 60.0f, 70.0f, 80.0f, 90.0f, 100.0f, 110.0f}, "Mixing bins - z-vertex"};

  } mixingConfig;

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  // Track selections (Except DCA selelctions) //

  struct : ConfigurableGroup {
    // Primary track selections
    Configurable<float> cMinPtcut{"cMinPtcut", 0.15, "Track minium pt cut"};
    Configurable<float> cMaxetacut{"cMaxetacut", 0.8, "Track maximum eta cut"};

    Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", true, "Primary track selection"};
    Configurable<bool> cfgPVContributor{"cfgPVContributor", true, "PV contributor track selection"};

    Configurable<int> cfgTPCcluster{"cfgTPCcluster", 1, "Minimum Number of TPC cluster"}; // Minmimum

    Configurable<int> cfgTPCRows{"cfgTPCRows", 80, "Minimum Number of TPC Crossed Rows "};

    Configurable<bool> cfgHasTOF{"cfgHasTOF", false, "Require TOF"};

    Configurable<float> cfgRapidityCut{"cfgRapidityCut", 0.5, "Rapidity cut for tracks"};

    // Primary track DCAxy to PV
    ConfigurableAxis cDCAtoPVBins{"cDCAtoPVBins", {1500, 0, 0.3}, "Bins for track DCA to PV"};

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

  } primarytrackConfig;

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  // Cascade and V0 selections //

  struct : ConfigurableGroup {
    // Topological selections for V0s
    ConfigurableAxis cDCADaughtersBins{"cDCADaughtersBins", {1000, 0, 0.1}, "Bins for DCA to daughters"};
    Configurable<float> cDCALambdaDaugtherscut{"cDCALambdaDaugtherscut", 0.7, "Lambda dauthers DCA cut"};
    Configurable<float> cDCALambdaToPVcut{"cDCALambdaToPVcut", 0.02, "Lambda DCA cut to PV"};
    Configurable<float> cDCAPionToPVcut{"cDCAPionToPVcut", 0.06, "pion DCA cut to PV"};
    Configurable<float> cDCAProtonToPVcut{"cDCAProtonToPVcut", 0.07, "proton DCA cut to PV"};

    Configurable<float> cV0CosPACutPtDepP0{"cV0CosPACutPtDepP0", 0.25, "Coeff. for Cosine Pointing angle for V0 as pt (p0)"};
    Configurable<float> cV0CosPACutPtDepP1{"cV0CosPACutPtDepP1", 0.022, "Coeff. for Cosine Pointing angle for V0 as pt (p1)"};

    Configurable<float> cMaxV0radiuscut{"cMaxV0radiuscut", 200., "V0 radius cut Maximum"};
    Configurable<float> cMinV0radiuscut{"cMinV0radiuscut", 2.5, "V0 radius cut Minimum"};
    Configurable<float> cMasswindowV0cut{"cMasswindowV0cut", 0.005, "V0 Mass window cut"};

  } v0sConfig;

  struct : ConfigurableGroup {
    // Topological selections for Cascades
    Configurable<float> cDCABachlorToPVcut{"cDCABachlorToPVcut", 0.06, "Bachelor DCA cut to PV"};
    Configurable<float> cDCAXiDaugthersCutPtRangeLower{"cDCAXiDaugthersCutPtRangeLower", 1., "Xi- DCA cut to PV as pt range lower"};
    Configurable<float> cDCAXiDaugthersCutPtRangeUpper{"cDCAXiDaugthersCutPtRangeUpper", 4., "Xi- DCA cut to PV as pt range upper"};
    Configurable<float> cDCAXiDaugthersCutPtDepLower{"cDCAXiDaugthersCutPtDepLower", 0.8, "Xi- DCA cut to PV as pt Under 1 GeV/c"};
    Configurable<float> cDCAXiDaugthersCutPtDepMiddle{"cDCAXiDaugthersCutPtDepMiddle", 0.5, "Xi- DCA cut to PV as pt 1 - 4 GeV/c"};
    Configurable<float> cDCAXiDaugthersCutPtDepUpper{"cDCAXiDaugthersCutPtDepUpper", 0.2, "Xi- DCA cut to PV as pt Over 4 GeV/c"};

    ConfigurableAxis cCosPABins{"cCosPABins", {3000, 0, 0.06}, "Bins for Cosine Pointing Angle"};
    Configurable<float> cCosPACascCutPtDepP0{"cCosPACascCutPtDepP0", 0.2, "Coeff. for Cosine Pointing angle for Cascade as pt (p0)"};
    Configurable<float> cCosPACascCutPtDepP1{"cCosPACascCutPtDepP1", 0.022, "Coeff. for Cosine Pointing angle for Cascade as pt (p1)"};

    Configurable<float> cMaxCascradiuscut{"cMaxCascradiuscut", 200., "Cascade radius cut Maximum"};
    Configurable<float> cMinCascradiuscut{"cMinCascradiuscut", 1.1, "Cascade radius cut Minimum"};
    Configurable<float> cMasswindowCasccut{"cMasswindowCasccut", 0.008, "Cascade Mass window cut"};
    Configurable<float> cMassXiminus{"cMassXiminus", 1.32171, "Mass of Xi baryon"};

  } cascadeConfig;

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  // PID Selections//

  struct : ConfigurableGroup {

    // Genenral selections
    Configurable<float> cPIDBound{"cPIDBound", 6.349, "configurable for replacing to .has"};

    Configurable<bool> tofAtHighPt{"tofAtHighPt", false, "Use TOF at high pT"};
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

  } pidConfig;

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  struct : ConfigurableGroup {

    // MC Event selection //
    Configurable<float> cZvertCutMC{"cZvertCutMC", 10.0, "MC Z-vertex cut"};
    // Cuts on mother particle and others

    Configurable<bool> cfgCutsOnMother{"cfgCutsOnMother", true, "Enamble additional cuts on mother"};
    Configurable<float> cMaxPtMotherCut{"cMaxPtMotherCut", 9.0, "Maximum pt of mother cut"};
    Configurable<float> cMaxMinvMotherCut{"cMaxMinvMotherCut", 3.0, "Maximum Minv of mother cut"};

    Configurable<bool> studyStableXi{"studyStableXi", false, "Study stable Xi"};

    Configurable<bool> cMultNTracksPVFull{"cMultNTracksPVFull", false, "Use full PV track multiplicity"};
    Configurable<bool> cMultNTracksPVeta1{"cMultNTracksPVeta1", false, "Use PV track multiplicity within |eta|<1"};
    Configurable<bool> cMultNTracksPVetaHalf{"cMultNTracksPVetaHalf", true, "Use PV track multiplicity within |eta|<0.5"};

    Configurable<bool> cRecoINELgt0{"cRecoINELgt0", true, "check if INEL>0 for reco events"};

    Configurable<bool> cUseFixedMassXi{"cUseFixedMassXi", false, "Use fixed mass for Xi-"};
    Configurable<bool> cUseTruthRapidity{"cUseTruthRapidity", false, "Use truth rapidity for Xi*"};

    Configurable<bool> cConsiderPairOnly{"cConsiderPairOnly", true, "Consider only existing particle pairs in the event"};

  } additionalConfig;

  TRandom* rn = new TRandom();

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  struct PidSelectionParam {
    float cMaxTPCnSigma;
    float cMaxTOFnSigma;
    bool cByPassTOF;
    float nsigmaCutCombined;
  };

  void init(o2::framework::InitContext&)
  {
    AxisSpec centAxis = {histoConfig.binsCent, "FT0M (%)"};
    AxisSpec dcaxyAxis = {primarytrackConfig.cDCAtoPVBins, "DCA_{#it{xy}} (cm)"};
    AxisSpec dcazAxis = {primarytrackConfig.cDCAtoPVBins, "DCA_{#it{z}} (cm)"};
    AxisSpec dcaDaugAxis = {v0sConfig.cDCADaughtersBins, "DCA_{#it{Daughter}} (cm)"};
    AxisSpec cosPAAxis = {cascadeConfig.cCosPABins, "1-cos(PA)"};
    AxisSpec mcLabelAxis = {6, -1.5, 4.5, "MC Label"};
    AxisSpec ptAxis = {histoConfig.binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {histoConfig.binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {histoConfig.cInvMassBins, histoConfig.cInvMassStart, histoConfig.cInvMassEnd, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec invMassAxisCasc = {800, 1.25, 1.65, "Invariant Mass for Casc. (GeV/#it{c}^2)"};
    AxisSpec pidQAAxis = {65, -6.5, 6.5};
    AxisSpec flagAxis = {9, 0, 9, "Flags"};

    if (histoConfig.multQA) {
      // multiplicity histograms
      histos.add("multQA/h2MultCent", "Multiplicity vs Centrality", HistType::kTH2F, {centAxis, histoConfig.multNTracksAxis});
      histos.add("multQA/h2MultCentMC", "Multiplicity vs Centrality MC", HistType::kTH2F, {centAxis, histoConfig.multNTracksAxis});
    }

    // event histograms
    if (histoConfig.eventQA) {
      histos.add("QAevent/hCollisionIndexSameE", "coll index sameE", HistType::kTH1F, {{500, 0.0f, 500.0f}});
      histos.add("QAevent/hCollisionIndexMixedE", "coll index mixedE", HistType::kTH1F, {{500, 0.0f, 500.0f}});

      histos.add("QAevent/hnTrksSameE", "n tracks per event SameE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});
      histos.add("QAevent/hnTrksMixedE", "n tracks per event MixedE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});

      histos.add("QAevent/hnCascsSameE", "n cascs per event SameE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});
      histos.add("QAevent/hnCascsMixedE", "n cascs per event MixedE", HistType::kTH1F, {{1000, 0.0f, 1000.0f}});

      histos.add("QAevent/hVertexZSameE", "Collision Vertex Z position", HistType::kTH1F, {{100, -15., 15.}});
      histos.add("QAevent/hMultiplicityPercentSameE", "Multiplicity percentile of collision", HistType::kTH1F, {centAxis});

      histos.add("QAevent/hVertexZMixedE", "Collision Vertex Z position", HistType::kTH1F, {{100, -15., 15.}});
      histos.add("QAevent/hMultiplicityPercentMixedE", "Multiplicity percentile of collision", HistType::kTH1F, {centAxis});
    }

    if (histoConfig.invMass1D) {
      histos.add("Xi1530invmassDS", "Invariant mass of Xi(1530)0 differnt sign", kTH1F, {invMassAxis});
      histos.add("Xi1530invmassLS", "Invariant mass of Xi(1530)0 like sign", kTH1F, {invMassAxis});

      histos.add("Xi1530invmassDSAnti", "Invariant mass of Anti-Xi(1530)0 differnt sign", kTH1F, {invMassAxis});
      histos.add("Xi1530invmassLSAnti", "Invariant mass of Anti-Xi(1530)0 like sign", kTH1F, {invMassAxis});
    }

    if (histoConfig.additionalQAplots) {
      // DCA QA to candidates for first pion and Xi-
      histos.add("QAbefore/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH2F, {ptAxis, dcaxyAxis});
      histos.add("QAbefore/trkDCAxy_Xi", "DCAxy distribution of Xi- track candidates", HistType::kTH2F, {ptAxis, dcaxyAxis});

      histos.add("QAbefore/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH2F, {ptAxis, dcazAxis});
      histos.add("QAbefore/trkDCAz_Xi", "DCAz distribution of Xi- track candidates", HistType::kTH2F, {ptAxis, dcazAxis});

      histos.add("QAafter/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH2F, {ptAxis, dcaxyAxis});
      histos.add("QAafter/trkDCAxy_Xi", "DCAxy distribution of Xi- track candidates", HistType::kTH2F, {ptAxis, dcaxyAxis});

      histos.add("QAafter/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH2F, {ptAxis, dcazAxis});
      histos.add("QAafter/trkDCAz_Xi", "DCAz distribution of Xi- track candidates", HistType::kTH2F, {ptAxis, dcazAxis});
    }

    if (histoConfig.pidPlots) {
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

    histos.add("h3Xi1530invmassDSAnti", "Invariant mass of Anti-Xi(1530)0 differnt sign", kTHnSparseF, {centAxis, ptAxis, invMassAxis, flagAxis});
    histos.add("h3XiinvmassDSAnti", "Invariant mass of Anti-Xi- differnt sign", kTHnSparseF, {centAxis, ptAxis, invMassAxisCasc, flagAxis});

    histos.add("h3Xi1530invmassLSAnti", "Invariant mass of Anti-Xi(1530)0 same sign", kTHnSparseF, {centAxis, ptAxis, invMassAxis, flagAxis});

    if (doprocessMEDF || doprocessMEMicro) {
      histos.add("h3Xi1530invmassME_DS", "Invariant mass of Xi(1530)0 mixed event DS", kTHnSparseF, {centAxis, ptAxis, invMassAxis, flagAxis});
      histos.add("h3Xi1530invmassME_DSAnti", "Invariant mass of Xi(1530)0 mixed event DSAnti", kTHnSparseF, {centAxis, ptAxis, invMassAxis, flagAxis});
    }

    if (doprocessMC) {
      // MC QA
      histos.add("QAMCTrue/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH2F, {ptAxis, dcaxyAxis});
      histos.add("QAMCTrue/trkDCAxy_xi", "DCAxy distribution of Xi- track candidates", HistType::kTH2F, {ptAxis, dcaxyAxis});

      histos.add("QAMCTrue/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH2F, {ptAxis, dcazAxis});
      histos.add("QAMCTrue/trkDCAz_xi", "DCAz distribution of Xi- track candidates", HistType::kTH2F, {ptAxis, dcazAxis});

      if (histoConfig.pidPlots) {
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

      histos.add("h3RecXi1530invmass", "Invariant mass of Reconstructed MC Xi(1530)0", kTHnSparseF, {centAxis, ptAxis, invMassAxis, ptAxis});
      histos.add("h3RecXiinvmass", "Invariant mass of Reconstructed MC Xi-", kTHnSparseF, {centAxis, ptAxis, invMassAxisCasc, ptAxis});

      histos.add("h3RecXi1530invmassAnti", "Invariant mass of Reconstructed MC Anti-Xi(1530)0", kTHnSparseF, {centAxis, ptAxis, invMassAxis, ptAxis});
      histos.add("h3RecXiinvmassAnti", "Invariant mass of Reconstructed MC Anti-Xi-", kTHnSparseF, {centAxis, ptAxis, invMassAxisCasc, ptAxis});

      histos.add("h3Xi1530Gen", "pT distribution of True MC Xi(1530)0", kTHnSparseF, {ptAxis, centAxis, histoConfig.multNTracksAxis});
      histos.add("h3Xi1530GenAnti", "pT distribution of True MC Anti-Xi(1530)0", kTHnSparseF, {ptAxis, centAxis, histoConfig.multNTracksAxis});

      histos.add("Xi1530Rec", "pT distribution of Reconstructed MC Xi(1530)0", kTH2F, {ptAxis, centAxis});
      histos.add("Xi1530RecAnti", "pT distribution of Reconstructed MC Anti-Xi(1530)0", kTH2F, {ptAxis, centAxis});
      histos.add("Xi1530Recinvmass", "Inv mass distribution of Reconstructed MC Xi(1530)0", kTH1F, {invMassAxis});
    }

    if (histoConfig.additionalQAplots) {
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

  float massPi = MassPionCharged;

  // Primary track selection for the first pion //
  template <bool IsResoMicrotrack, typename TrackType>
  bool primaryTrackCut(const TrackType track)
  {
    if (std::abs(track.eta()) > primarytrackConfig.cMaxetacut)
      return false;
    if (std::abs(track.pt()) < primarytrackConfig.cMinPtcut)
      return false;
    if constexpr (IsResoMicrotrack) {
      if (std::abs(o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(track.trackSelectionFlags())) > (primarytrackConfig.cDCAxytoPVByPtPiFirstP0 + primarytrackConfig.cDCAxyToPVByPtPiFirstExp * std::pow(track.pt(), -1.)))
        return false;
      if (primarytrackConfig.cDCAzToPVAsPt) {
        if (std::abs(o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(track.trackSelectionFlags())) > (primarytrackConfig.cDCAxytoPVByPtPiFirstP0 + primarytrackConfig.cDCAxyToPVByPtPiFirstExp * std::pow(track.pt(), -1.)))
          return false;
      } else {
        if (std::abs(o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(track.trackSelectionFlags())) > primarytrackConfig.cMaxDCAzToPVCut)
          return false;
      }
    } else {
      if (std::abs(track.dcaXY()) > (primarytrackConfig.cDCAxytoPVByPtPiFirstP0 + primarytrackConfig.cDCAxyToPVByPtPiFirstExp * std::pow(track.pt(), -1.)))
        return false;
      if (primarytrackConfig.cDCAzToPVAsPt) {
        if (std::abs(track.dcaZ()) > (primarytrackConfig.cDCAxytoPVByPtPiFirstP0 + primarytrackConfig.cDCAxyToPVByPtPiFirstExp * std::pow(track.pt(), -1.)))
          return false;
      } else {
        if (std::abs(track.dcaZ()) > primarytrackConfig.cMaxDCAzToPVCut)
          return false;
      }
      if (track.tpcNClsFound() < primarytrackConfig.cfgTPCcluster)
        return false;
      if (track.tpcNClsCrossedRows() < primarytrackConfig.cfgTPCRows)
        return false;
    }
    if (primarytrackConfig.cfgHasTOF && !track.hasTOF())
      return false;
    if (primarytrackConfig.cfgPVContributor && !track.isPVContributor())
      return false;
    if (primarytrackConfig.cfgPrimaryTrack && !track.isPrimaryTrack())
      return false;
    return true;
  }

  bool hasSubsystemInfo(float Nsigma)
  {
    return std::abs(Nsigma) < pidConfig.cPIDBound;
  }

  // Primary track selection for cascades, Need to more informations for cascades //
  template <typename TracksTypeCasc>
  bool cascprimaryTrackCut(const TracksTypeCasc track)
  {
    if (std::abs(track.eta()) > primarytrackConfig.cMaxetacut)
      return false;
    if (std::abs(track.pt()) < primarytrackConfig.cMinPtcut)
      return false;
    if (primarytrackConfig.cDCAxyToPVAsPtForCasc) {
      if (std::abs(track.dcaXYCascToPV()) > (primarytrackConfig.cDCAxyToPVByPtCascP0 + primarytrackConfig.cDCAxyToPVByPtCascExp * track.pt()))
        return false;
    }
    if (primarytrackConfig.cDCAzToPVAsPtForCasc) {
      if (std::abs(track.dcaZCascToPV()) > (primarytrackConfig.cDCAxyToPVByPtCascP0 + primarytrackConfig.cDCAxyToPVByPtCascExp * std::pow(track.pt(), -1.)))
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
    if (std::abs(track.daughDCA()) > v0sConfig.cDCALambdaDaugtherscut)
      return false;
    if (std::abs(track.dcav0topv()) < v0sConfig.cDCALambdaToPVcut)
      return false;
    if (track.sign() < 0) {
      if (std::abs(track.dcanegtopv()) < v0sConfig.cDCAPionToPVcut)
        return false;
      if (std::abs(track.dcapostopv()) < v0sConfig.cDCAProtonToPVcut)
        return false;
    } else {
      if (std::abs(track.dcanegtopv()) < v0sConfig.cDCAProtonToPVcut)
        return false;
      if (std::abs(track.dcapostopv()) < v0sConfig.cDCAPionToPVcut)
        return false;
    }
    if (track.v0CosPA() < std::cos(v0sConfig.cV0CosPACutPtDepP0 - v0sConfig.cV0CosPACutPtDepP1 * track.pt()))
      return false;
    if (track.transRadius() > v0sConfig.cMaxV0radiuscut || track.transRadius() < v0sConfig.cMinV0radiuscut)
      return false;
    if (std::abs(track.mLambda() - MassLambda) > v0sConfig.cMasswindowV0cut)
      return false;

    // Topological Cuts for Cascades
    if (std::abs(track.dcabachtopv()) < cascadeConfig.cDCABachlorToPVcut)
      return false;
    if (track.pt() < cascadeConfig.cDCAXiDaugthersCutPtRangeLower) {
      if (track.cascDaughDCA() > cascadeConfig.cDCAXiDaugthersCutPtDepLower)
        return false;
    }
    if (track.pt() >= cascadeConfig.cDCAXiDaugthersCutPtRangeLower && track.pt() < cascadeConfig.cDCAXiDaugthersCutPtRangeUpper) {
      if (track.cascDaughDCA() > cascadeConfig.cDCAXiDaugthersCutPtDepMiddle)
        return false;
    }
    if (track.pt() >= cascadeConfig.cDCAXiDaugthersCutPtRangeUpper) {
      if (track.cascDaughDCA() > cascadeConfig.cDCAXiDaugthersCutPtDepUpper)
        return false;
    }
    if (track.cascCosPA() < std::cos(cascadeConfig.cCosPACascCutPtDepP0 - cascadeConfig.cCosPACascCutPtDepP1 * track.pt()))
      return false;

    if (track.cascTransRadius() > cascadeConfig.cMaxCascradiuscut || track.cascTransRadius() < cascadeConfig.cMinCascradiuscut)
      return false;
    if (std::abs(track.mXi() - cascadeConfig.cMassXiminus) > cascadeConfig.cMasswindowCasccut)
      return false;

    return true;
  }

  bool pidSelector(float TPCNsigma, float TOFNsigma, const PidSelectionParam& params, bool tofAtHighPt, float trackPt)
  {
    bool tpcPIDPassed{false}, tofPIDPassed{false};

    if (tofAtHighPt && trackPt > pidConfig.cMinTOFpt) {
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

    PidSelectionParam pionFirstParams = {pidConfig.cMaxtpcnSigmaPionFirst, pidConfig.cMaxtofnSigmaPionFirst, pidConfig.cByPassTOFPionFirst, pidConfig.nsigmaCutCombinedPionFirst};

    return pidSelector(tpcNsigmaPionFirst, tofNsigmaPionFirst, pionFirstParams, pidConfig.tofAtHighPt, trackPt);
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

    PidSelectionParam bachelorParams = {pidConfig.cMaxtpcnSigmaPionBachelor, pidConfig.cMaxtofnSigmaPionBachelor, pidConfig.cByPassTOFPionBachelor, pidConfig.nsigmaCutCombinedPionBachelor};
    PidSelectionParam pionParams = {pidConfig.cMaxtpcnSigmaPion, pidConfig.cMaxtofnSigmaPion, pidConfig.cByPassTOFPion, pidConfig.nsigmaCutCombinedPion};
    PidSelectionParam protonParams = {pidConfig.cMaxtpcnSigmaProton, pidConfig.cMaxtofnSigmaProton, pidConfig.cByPassTOFProton, pidConfig.nsigmaCutCombinedProton};

    lConsistentWithXi = pidSelector(tpcNsigmaBachelor, tofNsigmaBachelor, bachelorParams, pidConfig.tofAtHighPt, trackPt);
    lConsistentWithPion = pidSelector(tpcNsigmaPion, tofNsigmaPion, pionParams, pidConfig.tofAtHighPt, trackPt);
    lConsistentWithProton = pidSelector(tpcNsigmaProton, tofNsigmaProton, protonParams, pidConfig.tofAtHighPt, trackPt);

    lConsistentWithLambda = lConsistentWithProton && lConsistentWithPion;

    return lConsistentWithXi && lConsistentWithLambda;
  }

  template <bool IsResoMicrotrack, bool IsMC, bool IsMix, typename CollisionType, typename centType, typename TracksType, typename TracksTypeCasc>
  void fillHistograms(const CollisionType& collision, const centType& inCent, const TracksType& dTracks1, const TracksTypeCasc& dTracks2) // Order: ResoColl, ResoTrack, ResoCascTrack
  {
    auto Cent = inCent;

    if (histoConfig.eventQA) {
      if constexpr (!IsMix) {
        histos.fill(HIST("QAevent/hVertexZSameE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentSameE"), Cent);
        histos.fill(HIST("QAevent/hCollisionIndexSameE"), collision.globalIndex());
        histos.fill(HIST("QAevent/hnTrksSameE"), dTracks1.size());
        histos.fill(HIST("QAevent/hnCascsSameE"), dTracks2.size());
      } else {
        histos.fill(HIST("QAevent/hVertexZMixedE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentMixedE"), Cent);
        histos.fill(HIST("QAevent/hCollisionIndexMixedE"), collision.globalIndex());
        histos.fill(HIST("QAevent/hnTrksMixedE"), dTracks1.size());
        histos.fill(HIST("QAevent/hnCascsMixedE"), dTracks2.size());
      }
    }

    if (additionalConfig.cConsiderPairOnly && (dTracks2.size() < 1 || dTracks1.size() < 1))
      return;

    LorentzVectorPtEtaPhiMass lDecayDaughter1, lDecayDaughter2, lResonance;
    std::vector<int64_t> pionCandateIndicies = {};
    std::vector<int64_t> xiCandateIndicies = {};
    pionCandateIndicies.reserve(dTracks1.size());
    xiCandateIndicies.reserve(dTracks2.size());

    for (const auto& trk1 : dTracks1) {
      auto trk1ptPi = trk1.pt();
      float trk1DCAXY = -1.f;
      float trk1DCAZ = -1.f;
      float trk1NSigmaPiTPC = -999.f;
      float trk1NSigmaPiTOF = -999.f;
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

      // QA before
      if constexpr (!IsMix) {
        if (histoConfig.pidPlots) {
          histos.fill(HIST("QAbefore/TPC_Nsigma_pi_first_all"), Cent, trk1ptPi, trk1NSigmaPiTPC);
          if (hasSubsystemInfo(trk1NSigmaPiTOF)) {
            histos.fill(HIST("QAbefore/TOF_Nsigma_pi_first_all"), Cent, trk1ptPi, trk1NSigmaPiTOF);
            histos.fill(HIST("QAbefore/TOF_TPC_Map_pi_first_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
          }
        }
        if (histoConfig.additionalQAplots) {
          histos.fill(HIST("QAbefore/trkDCAxy_pi"), trk1ptPi, trk1DCAXY);
          histos.fill(HIST("QAbefore/trkDCAz_pi"), trk1ptPi, trk1DCAZ);
        }
      }
      if (pidConfig.cUseOnlyTOFTrackPionBachelor && hasSubsystemInfo(trk1NSigmaPiTOF))
        continue;
      if (!selectionPIDPionFirst<IsResoMicrotrack>(trk1)) // PID selection for the first pion
        continue;
      if (!primaryTrackCut<IsResoMicrotrack>(trk1)) // Primary track selections
        continue;

      if constexpr (!IsMix) {
        if (histoConfig.pidPlots) {
          histos.fill(HIST("QAafter/TPC_Nsigma_pi_first_all"), Cent, trk1ptPi, trk1NSigmaPiTPC);

          if (hasSubsystemInfo(trk1NSigmaPiTOF)) {
            histos.fill(HIST("QAafter/TOF_Nsigma_pi_first_all"), Cent, trk1ptPi, trk1NSigmaPiTOF);
            histos.fill(HIST("QAafter/TOF_TPC_Map_pi_first_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
          }
        }
        if (histoConfig.additionalQAplots) {
          histos.fill(HIST("QAafter/trkDCAxy_pi"), trk1ptPi, trk1DCAXY);
          histos.fill(HIST("QAafter/trkDCAz_pi"), trk1ptPi, trk1DCAZ);
        }
      }

      pionCandateIndicies.push_back(trk1.index());
    }

    for (const auto& trk2 : dTracks2) {
      auto trk2ptXi = trk2.pt();
      auto trk2InvMass = trk2.mXi();
      auto trk2DCAXY = trk2.dcaXYCascToPV();
      auto trk2DCAZ = trk2.dcaZCascToPV();

      auto trk2DCAV0sDougthers = trk2.daughDCA();
      auto trk2DCACascDougthers = trk2.cascDaughDCA();
      auto massXiCand = trk2.mXi();
      auto trk2CascCosPA = trk2.cascCosPA();
      auto trk2V0sCosPA = trk2.v0CosPA();
      // QA before selections
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
        if (histoConfig.pidPlots) {
          histos.fill(HIST("QAbefore/TPC_Nsigma_pi_bachelor_all"), Cent, 0, trk2NSigmaPiBachelorTPC); // can't take pt information for the cascade secondary
          if (hasSubsystemInfo(trk2NSigmaPiBachelorTOF)) {
            histos.fill(HIST("QAbefore/TOF_TPC_Map_pi_bachelor_all"), trk2NSigmaPiBachelorTOF, trk2NSigmaPiBachelorTPC);
          }

          histos.fill(HIST("QAbefore/TPC_Nsigma_pr_all"), Cent, 0, trk2NSigmaPrPosTPC);
          if (hasSubsystemInfo(trk2NSigmaPrPosTOF)) {
            histos.fill(HIST("QAbefore/TOF_TPC_Map_pr_all"), trk2NSigmaPrPosTOF, trk2NSigmaPrPosTPC);
          }

          histos.fill(HIST("QAbefore/TPC_Nsigma_antipr_all"), Cent, 0, trk2NSigmaPrNegTPC);
          if (hasSubsystemInfo(trk2NSigmaPrNegTOF)) {
            histos.fill(HIST("QAbefore/TOF_TPC_Map_antipr_all"), trk2NSigmaPrNegTOF, trk2NSigmaPrNegTPC);
          }

          histos.fill(HIST("QAbefore/TPC_Nsigma_pi_all"), Cent, 0, trk2NSigmaPiPosTPC);
          if (hasSubsystemInfo(trk2NSigmaPiPosTOF)) {
            histos.fill(HIST("QAbefore/TOF_TPC_Map_pi_all"), trk2NSigmaPiPosTOF, trk2NSigmaPiPosTPC);
          }

          histos.fill(HIST("QAbefore/TPC_Nsigma_piminus_all"), Cent, 0, trk2NSigmaPiNegTPC);
          if (hasSubsystemInfo(trk2NSigmaPiNegTOF)) {
            histos.fill(HIST("QAbefore/TOF_TPC_Map_piminus_all"), trk2NSigmaPiNegTOF, trk2NSigmaPiNegTPC);
          }
        }

        if (histoConfig.additionalQAplots) {
          histos.fill(HIST("QAbefore/trkDCAxy_Xi"), trk2ptXi, trk2DCAXY);
          histos.fill(HIST("QAbefore/trkDCAz_Xi"), trk2ptXi, trk2DCAZ);
          histos.fill(HIST("QAbefore/V0sDCADoughter_aspt"), trk2ptXi, trk2DCAV0sDougthers);
          histos.fill(HIST("QAbefore/CascDCADoughter_aspt"), trk2ptXi, trk2DCACascDougthers);
          histos.fill(HIST("QAbefore/CascMass_aspt"), trk2ptXi, massXiCand);
          histos.fill(HIST("QAbefore/V0sCosPA_aspt"), trk2ptXi, 1. - trk2V0sCosPA);
          histos.fill(HIST("QAbefore/CascCosPA_aspt"), trk2ptXi, 1. - trk2CascCosPA);
        }
      }

      if (pidConfig.cUseOnlyTOFTrackPionBachelor && !hasSubsystemInfo(trk2NSigmaPiBachelorTOF))
        continue;
      if (pidConfig.cUseOnlyTOFTrackProton && !hasSubsystemInfo(trk2NSigmaPrPosTOF))
        continue;
      if (pidConfig.cUseOnlyTOFTrackProton && !hasSubsystemInfo(trk2NSigmaPrNegTOF))
        continue;
      if (pidConfig.cUseOnlyTOFTrackPion && !hasSubsystemInfo(trk2NSigmaPiPosTOF))
        continue;
      if (pidConfig.cUseOnlyTOFTrackPion && !hasSubsystemInfo(trk2NSigmaPiNegTOF))
        continue;
      if (!selectionPIDCascades(trk2))
        continue;
      if (!cascprimaryTrackCut(trk2) || !casctopCut(trk2)) // Primary track selections
        continue;

      // QA after selections
      if constexpr (!IsMix) {
        //// QA plots after the selection
        //  --- PID QA
        if (histoConfig.pidPlots) {

          if (trk2.sign() < 0) {
            histos.fill(HIST("QAafter/TPC_Nsigma_pi_bachelor_all"), Cent, 0, trk2NSigmaPiBachelorTPC); // not exist pt information in resocascade yet.
            if (hasSubsystemInfo(trk2NSigmaPiBachelorTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_pi_bachelor_all"), trk2NSigmaPiBachelorTOF, trk2NSigmaPiBachelorTPC);
            }

            histos.fill(HIST("QAafter/TPC_Nsigma_pr_all"), Cent, 0, trk2NSigmaPrPosTPC);
            if (hasSubsystemInfo(trk2NSigmaPrPosTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_pr_all"), trk2NSigmaPrPosTOF, trk2NSigmaPrPosTPC);
            }

            histos.fill(HIST("QAafter/TPC_Nsigma_piminus_all"), Cent, 0, trk2NSigmaPiNegTPC);
            if (hasSubsystemInfo(trk2NSigmaPiNegTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_piminus_all"), trk2NSigmaPiNegTOF, trk2NSigmaPiNegTPC);
            }

          } else {

            histos.fill(HIST("QAafter/TPC_Nsigma_pi_bachelor_all"), Cent, 0, trk2NSigmaPiBachelorTPC); // not exist pt information in resocascade yet.
            if (hasSubsystemInfo(trk2NSigmaPiBachelorTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_pi_bachelor_all"), trk2NSigmaPiBachelorTOF, trk2NSigmaPiBachelorTPC);
            }

            histos.fill(HIST("QAafter/TPC_Nsigma_antipr_all"), Cent, 0, trk2NSigmaPrNegTPC);
            if (hasSubsystemInfo(trk2NSigmaPrNegTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_antipr_all"), trk2NSigmaPrNegTOF, trk2NSigmaPrNegTPC);
            }

            histos.fill(HIST("QAafter/TPC_Nsigma_pi_all"), Cent, 0, trk2NSigmaPiPosTPC);
            if (hasSubsystemInfo(trk2NSigmaPiPosTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_pi_all"), trk2NSigmaPiPosTOF, trk2NSigmaPiPosTPC);
            }
          }
        }
        if (histoConfig.additionalQAplots) {
          histos.fill(HIST("QAafter/trkDCAxy_Xi"), trk2ptXi, trk2DCAXY);
          histos.fill(HIST("QAafter/trkDCAz_Xi"), trk2ptXi, trk2DCAZ);
          histos.fill(HIST("QAafter/V0sDCADoughter_aspt"), trk2ptXi, trk2DCAV0sDougthers);
          histos.fill(HIST("QAafter/CascDCADoughter_aspt"), trk2ptXi, trk2DCACascDougthers);
          histos.fill(HIST("QAafter/CascMass_aspt"), trk2ptXi, massXiCand);
          histos.fill(HIST("QAafter/V0sCosPA_aspt"), trk2ptXi, 1. - trk2V0sCosPA);
          histos.fill(HIST("QAafter/CascCosPA_aspt"), trk2ptXi, 1. - trk2CascCosPA);
        }

        if (additionalConfig.studyStableXi) {
          if (trk2.sign() < 0) {
            histos.fill(HIST("h3XiinvmassDS"), Cent, trk2ptXi, trk2InvMass, kData);
          } else if (trk2.sign() > 0) {
            histos.fill(HIST("h3XiinvmassDSAnti"), Cent, trk2ptXi, trk2InvMass, kData);
          }
          if constexpr (IsMC) {
            if (trk2.motherPDG() > 0) {
              histos.fill(HIST("h3RecXiinvmass"), Cent, trk2ptXi, trk2InvMass, kMCReco);
            } else {
              histos.fill(HIST("h3RecXiinvmassAnti"), Cent, trk2ptXi, trk2InvMass, kMCReco);
            }
          }
        }
      }

      xiCandateIndicies.push_back(trk2.index());
    }

    for (const auto& trk1cand : pionCandateIndicies) {
      auto pionCand = dTracks1.iteratorAt(trk1cand);
      for (const auto& trk2cand : xiCandateIndicies) {
        auto xiCand = dTracks2.iteratorAt(trk2cand);
        auto pionCandPt = pionCand.pt();
        auto xiCandPt = xiCand.pt();
        float massXiCand = xiCand.mXi();
        if (additionalConfig.cUseFixedMassXi)
          massXiCand = cascadeConfig.cMassXiminus;

        lDecayDaughter1 = LorentzVectorPtEtaPhiMass(pionCandPt, pionCand.eta(), pionCand.phi(), massPi);
        lDecayDaughter2 = LorentzVectorPtEtaPhiMass(xiCandPt, xiCand.eta(), xiCand.phi(), massXiCand);

        lResonance = lDecayDaughter1 + lDecayDaughter2;
        auto lResonanceMass = lResonance.M();
        auto lResonancePt = lResonance.Pt();

        if (std::abs(lResonance.Rapidity()) >= primarytrackConfig.cfgRapidityCut)
          continue;

        if (additionalConfig.cfgCutsOnMother) {
          if (lResonancePt >= additionalConfig.cMaxPtMotherCut) // excluding candidates in overflow
            continue;
          if (lResonanceMass >= additionalConfig.cMaxMinvMotherCut) // excluding candidates in overflow
            continue;
        }

        if (pionCand.sign() * xiCand.sign() < 0) { // Signal candidates

          if constexpr (!IsMix) {
            if (pionCand.sign() > 0) {
              if (histoConfig.invMass1D)
                histos.fill(HIST("Xi1530invmassDS"), lResonanceMass);
              histos.fill(HIST("h3Xi1530invmassDS"), Cent, lResonancePt, lResonanceMass, kData);
            } else if (pionCand.sign() < 0) {
              if (histoConfig.invMass1D)
                histos.fill(HIST("Xi1530invmassDSAnti"), lResonanceMass);

              histos.fill(HIST("h3Xi1530invmassDSAnti"), Cent, lResonancePt, lResonanceMass, kData);
            }
          } else {
            if (pionCand.sign() > 0) {
              histos.fill(HIST("h3Xi1530invmassME_DS"), Cent, lResonancePt, lResonanceMass, kData);
            } else if (pionCand.sign() < 0) {
              histos.fill(HIST("h3Xi1530invmassME_DSAnti"), Cent, lResonancePt, lResonanceMass, kData);
            }
          }

          if constexpr (IsMC) {
            if (std::abs(xiCand.motherPDG()) != kXiStar)
              continue;
            if (std::abs(pionCand.pdgCode()) != kPiPlus || std::abs(xiCand.pdgCode()) != kXiMinus)
              continue;
            if (pionCand.motherId() != xiCand.motherId())
              continue;
            auto lResonancePtMC = xiCand.motherPt();
            if (additionalConfig.cUseTruthRapidity && std::abs(xiCand.motherRap()) >= primarytrackConfig.cfgRapidityCut)
              continue;
            if (histoConfig.truthQA) {
              float trk1DCAXY = -1.f;
              float trk1DCAZ = -1.f;
              float trk1NSigmaPiTPC = -999.f;
              float trk1NSigmaPiTOF = -999.f;
              if constexpr (IsResoMicrotrack) {
                trk1DCAXY = o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAxy(pionCand.trackSelectionFlags());
                trk1DCAZ = o2::aod::resomicrodaughter::ResoMicroTrackSelFlag::decodeDCAz(pionCand.trackSelectionFlags());
                trk1NSigmaPiTPC = o2::aod::resomicrodaughter::PidNSigma::getTPCnSigma(pionCand.pidNSigmaPiFlag());
                trk1NSigmaPiTOF = o2::aod::resomicrodaughter::PidNSigma::getTOFnSigma(pionCand.pidNSigmaPiFlag());
              } else {
                trk1DCAXY = pionCand.dcaXY();
                trk1DCAZ = pionCand.dcaZ();
                trk1NSigmaPiTPC = pionCand.tpcNSigmaPi();
                trk1NSigmaPiTOF = pionCand.tofNSigmaPi();
              }

              auto trk2DCAXY = xiCand.dcaXYCascToPV();
              auto trk2DCAZ = xiCand.dcaZCascToPV();

              auto trk2DCAV0sDougthers = xiCand.daughDCA();
              auto trk2DCACascDougthers = xiCand.cascDaughDCA();
              auto massXiCand = xiCand.mXi();
              auto trk2CascCosPA = xiCand.cascCosPA();
              auto trk2V0sCosPA = xiCand.v0CosPA();

              // auto trk2ptPiBachelor = xiCand.pt();
              float trk2NSigmaPiBachelorTPC = xiCand.daughterTPCNSigmaBachPi();
              float trk2NSigmaPiBachelorTOF = xiCand.daughterTOFNSigmaBachPi();
              // auto trk2ptPr = xiCand.pt();
              float trk2NSigmaPrPosTPC = xiCand.daughterTPCNSigmaPosPr();
              float trk2NSigmaPrNegTPC = xiCand.daughterTPCNSigmaNegPr();

              float trk2NSigmaPrPosTOF = xiCand.daughterTOFNSigmaPosPr();
              float trk2NSigmaPrNegTOF = xiCand.daughterTOFNSigmaNegPr();

              // auto trk2ptPi = xiCand.pt();
              float trk2NSigmaPiPosTPC = xiCand.daughterTPCNSigmaPosPi();
              float trk2NSigmaPiNegTPC = xiCand.daughterTPCNSigmaNegPi();

              float trk2NSigmaPiPosTOF = xiCand.daughterTOFNSigmaPosPi();
              float trk2NSigmaPiNegTOF = xiCand.daughterTOFNSigmaNegPi();

              if (histoConfig.additionalQAplots) {
                histos.fill(HIST("QAMCTrue/V0sDCADoughter_aspt"), xiCandPt, trk2DCAV0sDougthers);
                histos.fill(HIST("QAMCTrue/CascDCADoughter_aspt"), xiCandPt, trk2DCACascDougthers);
                histos.fill(HIST("QAMCTrue/CascMass_aspt"), xiCandPt, massXiCand);
                histos.fill(HIST("QAMCTrue/V0sCosPA_aspt"), xiCandPt, 1. - trk2V0sCosPA);
                histos.fill(HIST("QAMCTrue/CascCosPA_aspt"), xiCandPt, 1. - trk2CascCosPA);

                histos.fill(HIST("QAMCTrue/trkDCAxy_pi"), pionCandPt, trk1DCAXY);
                histos.fill(HIST("QAMCTrue/trkDCAxy_xi"), xiCandPt, trk2DCAXY);

                histos.fill(HIST("QAMCTrue/trkDCAz_pi"), pionCandPt, trk1DCAZ);
                histos.fill(HIST("QAMCTrue/trkDCAz_xi"), xiCandPt, trk2DCAZ);
              }

              histos.fill(HIST("QAMCTrue/TPC_Nsigma_pi_first_all"), Cent, pionCandPt, trk1NSigmaPiTPC);
              if (hasSubsystemInfo(trk1NSigmaPiTOF)) {
                histos.fill(HIST("QAMCTrue/TOF_Nsigma_pi_first_all"), Cent, pionCandPt, trk1NSigmaPiTOF);
                histos.fill(HIST("QAMCTrue/TOF_TPC_Map_pi_first_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
              }

              histos.fill(HIST("QAMCTrue/TPC_Nsigma_pi_bachelor_all"), Cent, 0, trk2NSigmaPiBachelorTPC); // not exist pt information in resocascade yet.
              if (hasSubsystemInfo(trk2NSigmaPiBachelorTOF)) {
                histos.fill(HIST("QAMCTrue/TOF_TPC_Map_pi_bachelor_all"), trk2NSigmaPiBachelorTOF, trk2NSigmaPiBachelorTPC);
              }

              histos.fill(HIST("QAMCTrue/TPC_Nsigma_pr_all"), Cent, 0, trk2NSigmaPrPosTPC);
              if (hasSubsystemInfo(trk2NSigmaPrPosTOF)) {
                histos.fill(HIST("QAMCTrue/TOF_TPC_Map_pr_all"), trk2NSigmaPrPosTOF, trk2NSigmaPrPosTPC);
              }
              histos.fill(HIST("QAMCTrue/TPC_Nsigma_antipr_all"), Cent, 0, trk2NSigmaPrNegTPC);
              if (hasSubsystemInfo(trk2NSigmaPrNegTOF)) {
                histos.fill(HIST("QAMCTrue/TOF_TPC_Map_antipr_all"), trk2NSigmaPrNegTOF, trk2NSigmaPrNegTPC);
              }

              histos.fill(HIST("QAMCTrue/TPC_Nsigma_piminus_all"), Cent, 0, trk2NSigmaPiNegTPC);
              if (hasSubsystemInfo(trk2NSigmaPiNegTOF)) {
                histos.fill(HIST("QAMCTrue/TOF_TPC_Map_piminus_all"), trk2NSigmaPiNegTOF, trk2NSigmaPiNegTPC);
              }
              histos.fill(HIST("QAMCTrue/TPC_Nsigma_pi_all"), Cent, 0, trk2NSigmaPiPosTPC);
              if (hasSubsystemInfo(trk2NSigmaPiPosTOF)) {
                histos.fill(HIST("QAMCTrue/TOF_TPC_Map_pi_all"), trk2NSigmaPiPosTOF, trk2NSigmaPiPosTPC);
              }

            } // truthQA

            // MC histograms
            if (xiCand.motherPDG() > 0) {
              histos.fill(HIST("Xi1530Rec"), lResonancePt, Cent);
              histos.fill(HIST("Xi1530Recinvmass"), lResonanceMass);
              histos.fill(HIST("h3RecXi1530invmass"), Cent, lResonancePt, lResonanceMass, lResonancePtMC);
            } else {
              histos.fill(HIST("Xi1530RecAnti"), lResonancePt, Cent);
              histos.fill(HIST("Xi1530Recinvmass"), lResonanceMass);
              histos.fill(HIST("h3RecXi1530invmassAnti"), Cent, lResonancePt, lResonanceMass, lResonancePtMC);
            }
          } // is MC
        } else { // Bkg candidates
          if constexpr (!IsMix) {
            if (pionCand.sign() < 0) {
              if (histoConfig.invMass1D)
                histos.fill(HIST("Xi1530invmassLS"), lResonanceMass);
              histos.fill(HIST("h3Xi1530invmassLS"), Cent, lResonancePt, lResonanceMass, kLS);
            } else if (pionCand.sign() > 0) {
              if (histoConfig.invMass1D)
                histos.fill(HIST("Xi1530invmassLSAnti"), lResonanceMass);
              histos.fill(HIST("h3Xi1530invmassLSAnti"), Cent, lResonancePt, lResonanceMass, kLS);
            }
          }
        } // -> End if signal or bkg
      } // -> End loop over xi candidates
    } // -> End loop over pion and xi candidates
  } // -> End fillHistograms

  void processData(aod::ResoCollision const& resoCollision,
                   aod::ResoTracks const& resoTracks,
                   aod::ResoCascades const& cascTracks)
  {
    auto inCent = resoCollision.cent();
    auto multiplicity = 0.f;

    if (additionalConfig.cRecoINELgt0 && !resoCollision.isRecINELgt0()) // Check reco INELgt0 (at least one PV track in |eta| < 1) about the collision
      return;
    if (additionalConfig.cMultNTracksPVFull)
      multiplicity = resoCollision.multNTracksPV();
    else if (additionalConfig.cMultNTracksPVeta1)
      multiplicity = resoCollision.multNTracksPVeta1();
    else if (additionalConfig.cMultNTracksPVetaHalf)
      multiplicity = resoCollision.multNTracksPVetaHalf();

    if (histoConfig.multQA) {
      histos.fill(HIST("multQA/h2MultCent"), inCent, multiplicity);
    }
    fillHistograms<false, false, false>(resoCollision, inCent, resoTracks, cascTracks);
  }

  // Calculate numerator for the Acceptance x Efficiency
  void processMC(ResoMCCols::iterator const& resoCollision,
                 soa::Join<aod::ResoCascades, aod::ResoMCCascades> const& cascTracks,
                 soa::Join<aod::ResoTracks, aod::ResoMCTracks> const& resoTracks)
  {
    auto inCent = resoCollision.cent();
    float multiplicity = 0.f;
    if (additionalConfig.cRecoINELgt0 && !resoCollision.isRecINELgt0()) // Check reco INELgt0 (at least one PV track in |eta| < 1) about the collision
      return;
    if (additionalConfig.cMultNTracksPVFull)
      multiplicity = resoCollision.multNTracksPV();
    else if (additionalConfig.cMultNTracksPVeta1)
      multiplicity = resoCollision.multNTracksPVeta1();
    else if (additionalConfig.cMultNTracksPVetaHalf)
      multiplicity = resoCollision.multNTracksPVetaHalf();

    if (!resoCollision.isInAfterAllCuts()) // MC event selection, all cuts missing vtx cut
      return;
    if (histoConfig.multQA) {
      histos.fill(HIST("multQA/h2MultCent"), inCent, multiplicity);
    }
    fillHistograms<false, true, false>(resoCollision, inCent, resoTracks, cascTracks);
  }

  // Calculate denominator for the Acceptance x Efficiency, actually it is not Trueth info...
  void processMCTrue(ResoMCCols::iterator const& resoCollision,
                     aod::ResoMCParents const& resoParents)
  {
    auto multiplicity = resoCollision.mcMultiplicity();
    auto inCent = resoCollision.cent();
    if (additionalConfig.cRecoINELgt0 && !resoCollision.isRecINELgt0()) // Check reco INELgt0 (at least one PV track in |eta| < 1) about the collision
      return;
    if (!resoCollision.isInAfterAllCuts())
      return;
    if (histoConfig.multQA) {
      histos.fill(HIST("multQA/h2MultCentMC"), inCent, multiplicity);
    }

    for (const auto& part : resoParents) { // loop over all pre-filtered MC particles
      if (std::abs(part.pdgCode()) != kXiStar || std::abs(part.y()) >= primarytrackConfig.cfgRapidityCut)
        continue;
      bool pass1 = std::abs(part.daughterPDG1()) == kPiPlus || std::abs(part.daughterPDG2()) == kPiPlus;
      bool pass2 = std::abs(part.daughterPDG1()) == kXiMinus || std::abs(part.daughterPDG2()) == kXiMinus;

      if (!pass1 || !pass2)
        continue;

      if (resoCollision.isInAfterAllCuts()) // after all event selection
      {
        if (part.pdgCode() > 0)
          histos.fill(HIST("h3Xi1530Gen"), part.pt(), inCent, multiplicity);
        else
          histos.fill(HIST("h3Xi1530GenAnti"), part.pt(), inCent, multiplicity);
      }
    }
  }

  void processDataMicro(aod::ResoCollision const& resoCollision,
                        aod::ResoMicroTracks const& resomicrotracks,
                        aod::ResoCascades const& cascTracks)
  {
    float multiplicity = 0.f;
    auto inCent = resoCollision.cent();
    if (additionalConfig.cRecoINELgt0 && !resoCollision.isRecINELgt0()) // Check reco INELgt0 (at least one PV track in |eta| < 1) about the collision
      return;
    if (additionalConfig.cMultNTracksPVFull)
      multiplicity = resoCollision.multNTracksPV();
    else if (additionalConfig.cMultNTracksPVeta1)
      multiplicity = resoCollision.multNTracksPVeta1();
    else if (additionalConfig.cMultNTracksPVetaHalf)
      multiplicity = resoCollision.multNTracksPVetaHalf();

    if (histoConfig.multQA) {
      histos.fill(HIST("multQA/h2MultCent"), inCent, multiplicity);
    }
    fillHistograms<true, false, false>(resoCollision, inCent, resomicrotracks, cascTracks);
  }

  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processMEMicro(aod::ResoCollisions const& resoCollisions,
                      aod::ResoMicroTracks const& resomicrotracks,
                      aod::ResoCascades const& cascTracks)
  {
    auto tracksTuple = std::make_tuple(resomicrotracks, cascTracks);

    BinningTypeVtxZT0M colBinning{{mixingConfig.cfgVtxBins, mixingConfig.cfgMultBins}, true};
    Pair<aod::ResoCollisions, aod::ResoMicroTracks, aod::ResoCascades, BinningTypeVtxZT0M> pairs{colBinning, mixingConfig.nEvtMixing, -1, resoCollisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {

      float multiplicity = 0.f;
      auto inCent = collision1.cent();
      if (additionalConfig.cRecoINELgt0 && !collision1.isRecINELgt0()) // Check reco INELgt0 (at least one PV track in |eta| < 1) about the collision
        continue;
      if (additionalConfig.cMultNTracksPVFull)
        multiplicity = collision1.multNTracksPV();
      else if (additionalConfig.cMultNTracksPVeta1)
        multiplicity = collision1.multNTracksPVeta1();
      else if (additionalConfig.cMultNTracksPVetaHalf)
        multiplicity = collision1.multNTracksPVetaHalf();

      if (histoConfig.multQA) {
        histos.fill(HIST("multQA/h2MultCent"), inCent, multiplicity);
      }
      fillHistograms<true, false, true>(collision1, inCent, tracks1, tracks2);
    }
  }
  void processMEDF(aod::ResoCollisionDFs const& resoCollisions, aod::ResoTrackDFs const& resotracks, aod::ResoCascadeDFs const& cascTracks)
  {

    auto tracksTuple = std::make_tuple(resotracks, cascTracks);

    BinningTypeVtxZT0M colBinning{{mixingConfig.cfgVtxBins, mixingConfig.cfgMultBins}, true};
    Pair<aod::ResoCollisionDFs, aod::ResoTrackDFs, aod::ResoCascadeDFs, BinningTypeVtxZT0M> pairs{colBinning, mixingConfig.nEvtMixing, -1, resoCollisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {

      float multiplicity = 0.f;
      auto inCent = collision1.cent();
      if (histoConfig.multQA) {
        histos.fill(HIST("multQA/h2MultCent"), inCent, multiplicity);
      }
      fillHistograms<false, false, true>(collision1, inCent, tracks1, tracks2);
    }
  }

  PROCESS_SWITCH(Xi1530Analysisqa, processData, "Process Event for Data", false);
  PROCESS_SWITCH(Xi1530Analysisqa, processMC, "Process Event for MC (Reconstructed)", false);
  PROCESS_SWITCH(Xi1530Analysisqa, processMCTrue, "Process Event for MC (Generated)", false);
  PROCESS_SWITCH(Xi1530Analysisqa, processDataMicro, "Process Event for Data (MicroTrack)", false);
  PROCESS_SWITCH(Xi1530Analysisqa, processMEMicro, "Process EventMixing (MicroTrack) ", false);
  PROCESS_SWITCH(Xi1530Analysisqa, processMEDF, "Process EventMixing (DataFrame) ", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Xi1530Analysisqa>(cfgc)};
}
