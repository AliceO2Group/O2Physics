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

#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/BinningPolicy.h>
#include <Framework/Configurable.h>
#include <Framework/GroupedCombinations.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/SliceCache.h>
#include <Framework/runDataProcessing.h>

#include <Math/GenVector/RotationZ.h>
#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TPDGCode.h>

#include <arrow/array/concatenate.h>
#include <arrow/chunked_array.h>
#include <arrow/result.h>
#include <arrow/table.h>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <map>
#include <memory>
#include <tuple>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using namespace o2::constants::physics;
using LorentzVectorPtEtaPhiMass = ROOT::Math::PtEtaPhiMVector;
using LorentzVectorSetXYZM = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<float>>;

enum {
  kData = 0,
  kLS = 1,
  kMixing = 2,
  kMCReco = 3,
  kMCTrue = 4,
  kMCTruePS = 5,
  kINEL10 = 6,
  kINELg010 = 7,
  kAllType = 8,
  kXiStar = 3324
};

struct Xi1530Analysisqa {

  // Module-initializer tables; full tracks and cascades retain their schema.
  using ResoCollisions = aod::ResoCollisions_001;
  // Read original IDs only; no track_as() or original AO2D is needed here.
  using ResoTracks = soa::Join<aod::ResoTracks, aod::ResoTrackTracks>;
  using ResoMicroTracks = aod::ResoMicroTracks_001;
  using ResoMCCols = soa::Join<ResoCollisions, aod::ResoMCCollisions_001>;

  // Own just the event rows: zero-copy slices would retain entire DF buffers.
  struct MixingCollision {
    float x, y, z, centrality, multiplicity;
    int64_t index;
    bool recINELgt0;
    [[nodiscard]] float posX() const { return x; }
    [[nodiscard]] float posY() const { return y; }
    [[nodiscard]] float posZ() const { return z; }
    [[nodiscard]] int64_t globalIndex() const { return index; }
  };
  struct MixingEvent {
    MixingCollision collision{};
    std::shared_ptr<arrow::Table> tracks;
  };
  std::map<int, std::deque<MixingEvent>> mixingPools;
  float mixingBField = 0.f;

  // Basic set-up //
  SliceCache cache;
  Preslice<aod::ResoTracks> perRCol = aod::resodaughter::resoCollisionId;
  Preslice<aod::Tracks> perCollision = aod::track::collisionId;
  Preslice<ResoMicroTracks> perResoCollision =
    aod::resodaughter::resoCollisionId;
  Preslice<aod::ResoCascades> perResoCollisionCasc =
    aod::resodaughter::resoCollisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

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

    Configurable<float> cfgRapidityMinCut{"cfgRapidityMinCut", -0.5, "Rapidity cut for tracks"};
    Configurable<float> cfgRapidityMaxCut{"cfgRapidityMaxCut", 0.5, "Rapidity cut for tracks"};

    // Primary track DCAxy to PV
    ConfigurableAxis cDCAtoPVBins{"cDCAtoPVBins", {1500, 0, 0.3}, "Bins for track DCA to PV"};

    Configurable<bool> cDCAxyToPVAsPt{"cDCAxyToPVAsPt", true, "Use the pT-dependent pion DCAxy cut; otherwise use cMaxDCAxyToPVCut"};
    Configurable<float> cMaxDCAxyToPVCut{"cMaxDCAxyToPVCut", 0.5f, "Constant maximum pion |DCAxy| (cm) when its pT-dependent cut is disabled"};
    Configurable<float> cDCAxytoPVByPtPiFirstP0{"cDCAxytoPVByPtPiFirstP0", 0.004, "Constant term of the pion pT-dependent DCAxy/z cut (cm)"};
    Configurable<float> cDCAxyToPVByPtPiFirstExp{"cDCAxyToPVByPtPiFirstExp", 0.013, "Coefficient in pion DCAxy/z cut = P0 + coefficient / pT^power (legacy key)"};
    Configurable<float> cDCAxyToPVByPtPiFirstPower{"cDCAxyToPVByPtPiFirstPower", 1.f, "Power in pion DCAxy/z cut = P0 + coefficient / pT^power"};

    Configurable<bool> cDCAxyToPVAsPtForCasc{"cDCAxyToPVAsPtForCasc", true, "Set DCAxy to PV selection as pt cut"};
    Configurable<float> cDCAxyToPVByPtCascP0{"cDCAxyToPVByPtCascP0", 999., "Coeff. for Track DCAxy cut to PV by pt for Cascade (p0)"};
    Configurable<float> cDCAxyToPVByPtCascExp{"cDCAxyToPVByPtCascExp", 1., "Coeff. Track DCAxy cut to PV by pt for Cascade (exp)"};

    // Primary track DCAz to PV
    Configurable<bool> cDCAzToPVAsPt{"cDCAzToPVAsPt", true, "Use the shared pion pT-dependent DCA cut for DCAz; otherwise use cMaxDCAzToPVCut"};
    Configurable<bool> cDCAzToPVAsPtForCasc{"cDCAzToPVAsPtForCasc", true, "Set DCA to PV selection as pt cut"};
    Configurable<float> cMaxDCAzToPVCut{"cMaxDCAzToPVCut", 0.5, "Track DCAz cut to PV Maximum"};
    Configurable<float> cMaxDCAzToPVCutCasc{"cMaxDCAzToPVCutCasc", 0.5, "Track DCAz cut to PV Maximum for casc"};

  } primarytrackConfig;

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  // Cascade and V0 selections //

  struct : ConfigurableGroup {
    // Topological selections for V0s
    ConfigurableAxis cDCADaughtersBins{"cDCADaughtersBins", {1000, 0, 0.1}, "Bins for DCA to daughters"};
    ConfigurableAxis cTransRadiusBins{"cTransRadiusBins", {500, 0, 100}, "Bins for transverse radius"};
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
    ConfigurableAxis cDCASecondaryBins{"cDCASecondaryBins", {1000, 0, 0.1}, "Bins for DCA to Cascade secondary"};
    ConfigurableAxis cProperLifetimeBins{"cProperLifetimeBins", {300, 0, 30}, "Bins for proper lifetime"};
    ConfigurableAxis cDCABachelorToPVcutBins{"cDCABachelorToPVcutBins", {1000, 0, 0.1}, "Bins for DCA bachelor to PV cut"};
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

    Configurable<float> cMaxProperLifetimeCut{"cMaxProperLifetimeCut", 4.7, "Maximum proper lifetime cut for Xi- candidates"};
    Configurable<float> cMinNCrossedRowsTPCPos{"cMinNCrossedRowsTPCPos", 50, "Minimum number of crossed rows in TPC for positive track in cascade"};
    Configurable<float> cMinNCrossedRowsTPCNeg{"cMinNCrossedRowsTPCNeg", 50, "Minimum number of crossed rows in TPC for negative track in cascade"};
    Configurable<float> cMinNCrossedRowsTPCBach{"cMinNCrossedRowsTPCBach", 50, "Minimum number of crossed rows in TPC for bachelor track in cascade"};

  } cascadeConfig;

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  // PID Selections//

  struct : ConfigurableGroup {

    // Genenral selections
    Configurable<float> cPIDBound{"cPIDBound", 6.349, "configurable for replacing to .has"};
    ConfigurableAxis cPIDnSigmaBins{"cPIDnSigmaBins", {131, -6.5, 6.5}, "Bins for nSigma PID"};

    Configurable<bool> tofAtHighPt{"tofAtHighPt", false, "Apply TOF PID only above cMinTOFpt (cascade pT for its daughters)"};
    Configurable<float> cMinTOFpt{"cMinTOFpt", 0.5, "Lower pT threshold for TOF PID when tofAtHighPt is enabled"};

    // PID Selections for Pion First
    Configurable<float> cMintpcnSigmaPionFirst{"cMintpcnSigmaPionFirst", -4.0, "Strict lower TPC (nSigma - mean) bound for Pion First"};
    Configurable<float> cMaxtpcnSigmaPionFirst{"cMaxtpcnSigmaPionFirst", 4.0, "Strict upper TPC (nSigma - mean) bound for Pion First"};
    Configurable<float> cMeantpcnSigmaPionFirst{"cMeantpcnSigmaPionFirst", 0.0, "TPC nSigma mean subtracted for Pion First PID"};
    Configurable<float> cMintofnSigmaPionFirst{"cMintofnSigmaPionFirst", -3.0, "Strict lower TOF (nSigma - mean) bound for Pion First"};
    Configurable<float> cMaxtofnSigmaPionFirst{"cMaxtofnSigmaPionFirst", 3.0, "Strict upper TOF (nSigma - mean) bound for Pion First"};
    Configurable<float> cMeantofnSigmaPionFirst{"cMeantofnSigmaPionFirst", 0.0, "TOF nSigma mean subtracted for Pion First PID"};

    Configurable<float> nsigmaCutCombinedPionFirst{"nsigmaCutCombinedPionFirst", -4.0, "Combined nSigma cut for Pion First"};

    Configurable<bool> cByPassTOFPionFirst{"cByPassTOFPionFirst", true, "By pass TOF Pion First PID selection"};

    // PID Selections for Pion Bachelor
    Configurable<float> cMintpcnSigmaPionBachelor{"cMintpcnSigmaPionBachelor", -4.0, "Strict lower TPC (nSigma - mean) bound for Pion Bachelor"};
    Configurable<float> cMaxtpcnSigmaPionBachelor{"cMaxtpcnSigmaPionBachelor", 4.0, "Strict upper TPC (nSigma - mean) bound for Pion Bachelor"};
    Configurable<float> cMeantpcnSigmaPionBachelor{"cMeantpcnSigmaPionBachelor", 0.0, "TPC nSigma mean subtracted for Pion Bachelor PID"};
    Configurable<float> cMintofnSigmaPionBachelor{"cMintofnSigmaPionBachelor", -3.0, "Strict lower TOF (nSigma - mean) bound for Pion Bachelor"};
    Configurable<float> cMaxtofnSigmaPionBachelor{"cMaxtofnSigmaPionBachelor", 3.0, "Strict upper TOF (nSigma - mean) bound for Pion Bachelor"};
    Configurable<float> cMeantofnSigmaPionBachelor{"cMeantofnSigmaPionBachelor", 0.0, "TOF nSigma mean subtracted for Pion Bachelor PID"};

    Configurable<float> nsigmaCutCombinedPionBachelor{"nsigmaCutCombinedPionBachelor", -4.0, "Combined nSigma cut for Pion Bachelor"};

    Configurable<bool> cByPassTOFPionBachelor{"cByPassTOFPionBachelor", true, "By pass TOF Pion Bachelor PID selection"};

    // PID Selections for Pion
    Configurable<float> cMintpcnSigmaPion{"cMintpcnSigmaPion", -4.0, "Strict lower TPC (nSigma - mean) bound for Pion"};
    Configurable<float> cMaxtpcnSigmaPion{"cMaxtpcnSigmaPion", 4.0, "Strict upper TPC (nSigma - mean) bound for Pion"};
    Configurable<float> cMeantpcnSigmaPion{"cMeantpcnSigmaPion", 0.0, "TPC nSigma mean subtracted for Pion PID"};
    Configurable<float> cMintofnSigmaPion{"cMintofnSigmaPion", -3.0, "Strict lower TOF (nSigma - mean) bound for Pion"};
    Configurable<float> cMaxtofnSigmaPion{"cMaxtofnSigmaPion", 3.0, "Strict upper TOF (nSigma - mean) bound for Pion"};
    Configurable<float> cMeantofnSigmaPion{"cMeantofnSigmaPion", 0.0, "TOF nSigma mean subtracted for Pion PID"};

    Configurable<float> nsigmaCutCombinedPion{"nsigmaCutCombinedPion", -4.0, "Combined nSigma cut for Pion"};

    Configurable<bool> cByPassTOFPion{"cByPassTOFPion", true, "By pass TOF Pion PID selection"};

    // PID Selections for Proton
    Configurable<float> cMintpcnSigmaProton{"cMintpcnSigmaProton", -4.0, "Strict lower TPC (nSigma - mean) bound for Proton"};
    Configurable<float> cMaxtpcnSigmaProton{"cMaxtpcnSigmaProton", 4.0, "Strict upper TPC (nSigma - mean) bound for Proton"};
    Configurable<float> cMeantpcnSigmaProton{"cMeantpcnSigmaProton", 0.0, "TPC nSigma mean subtracted for Proton PID"};
    Configurable<float> cMintofnSigmaProton{"cMintofnSigmaProton", -3.0, "Strict lower TOF (nSigma - mean) bound for Proton"};
    Configurable<float> cMaxtofnSigmaProton{"cMaxtofnSigmaProton", 3.0, "Strict upper TOF (nSigma - mean) bound for Proton"};
    Configurable<float> cMeantofnSigmaProton{"cMeantofnSigmaProton", 0.0, "TOF nSigma mean subtracted for Proton PID"};

    Configurable<float> nsigmaCutCombinedProton{"nsigmaCutCombinedProton", -4.0, "Combined nSigma cut for Proton"};

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

    Configurable<bool> cRecoINELgt0{"cRecoINELgt0", true, "check if INEL>0 for reco events"};
    Configurable<bool> cMCINELgt0{"cMCINELgt0", true, "Require generator INEL>0 in reconstructed MC and generated-parent processes"};
    Configurable<bool> cMCVtxIn10{"cMCVtxIn10", true, "Require generator |vertex z| < 10 cm using isVtxIn10 in reconstructed MC and generated-parent processes"};

    Configurable<bool> cUseFixedMassXi{"cUseFixedMassXi", false, "Use fixed mass for Xi-"};
    Configurable<bool> cUseTruthRapidity{"cUseTruthRapidity", true, "Use truth rapidity for Xi*"};

    Configurable<bool> cfgFillRotBkg{"cfgFillRotBkg", true, "Fill rotated background"};
    Configurable<float> cfgMinRot{"cfgMinRot", 5.0 * constants::math::PI / 6.0, "Minimum of rotation"};
    Configurable<float> cfgMaxRot{"cfgMaxRot", 7.0 * constants::math::PI / 6.0, "Maximum of rotation"};
    Configurable<bool> cfgRotPion{"cfgRotPion", true, "Rotate pion"};
    Configurable<int> cfgNrotBkg{"cfgNrotBkg", 4, "Number of rotated copies (background) per each original candidate"};
  } additionalConfig;

  //*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//*//

  struct PidSelectionParam {
    float cMinTPCnSigma;
    float cMaxTPCnSigma;
    float cMeanTPCnSigma;
    float cMinTOFnSigma;
    float cMaxTOFnSigma;
    float cMeanTOFnSigma;
    bool cByPassTOF;
    float nsigmaCutCombined;
  };

  void init(o2::framework::InitContext&)
  {
    if (doprocessMEDF && doprocessMEMicro) {
      LOGF(fatal, "Enable only one mixing process: processMEDF or processMEMicro");
    }
    AxisSpec centAxis = {histoConfig.binsCent, "FT0M (%)"};
    AxisSpec dcaxyAxis = {primarytrackConfig.cDCAtoPVBins, "|DCA_{#it{xy}}| (cm)"};
    AxisSpec dcazAxis = {primarytrackConfig.cDCAtoPVBins, "|DCA_{#it{z}}| (cm)"};
    AxisSpec dcaSecondaryAxis = {cascadeConfig.cDCASecondaryBins, "DCA_{#it{Secondary}} (cm)"};
    AxisSpec dcaBachAxis = {cascadeConfig.cDCABachelorToPVcutBins, "DCA_{#it{Bach}} (cm)"};
    AxisSpec dcaDaugAxis = {v0sConfig.cDCADaughtersBins, "DCA_{#it{Daughter}} (cm)"};
    AxisSpec cosPAAxis = {cascadeConfig.cCosPABins, "1-cos(PA)"};
    AxisSpec properLifetimeAxis = {cascadeConfig.cProperLifetimeBins, "c#tau (cm)"};
    AxisSpec mcLabelAxis = {6, -1.5, 4.5, "MC Label"};
    AxisSpec ptAxis = {histoConfig.binsPt, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec ptAxisQA = {histoConfig.binsPtQA, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invMassAxis = {histoConfig.cInvMassBins, histoConfig.cInvMassStart, histoConfig.cInvMassEnd, "Invariant Mass (GeV/#it{c}^2)"};
    AxisSpec invMassAxisCasc = {400, 1.121, 1.521, "Invariant Mass for Casc. (GeV/#it{c}^2)"};
    AxisSpec invMassAxisLambda = {400, 0.921, 1.321, "Invariant Mass for Lambda (GeV/#it{c}^2)"};
    AxisSpec transRadiusAxis = {v0sConfig.cTransRadiusBins, "Transverse radius (cm)"};
    AxisSpec pidQAAxis = {pidConfig.cPIDnSigmaBins, "nSigma PID"};
    AxisSpec flagAxis = {9, 0, 9, "Flags"};

    if (histoConfig.multQA) {
      // multiplicity histograms
      histos.add("multQA/h2MultCent", "Same-event multiplicity vs centrality (one entry per event)", HistType::kTH2F, {centAxis, histoConfig.multNTracksAxis});
      if (doprocessMEMicro || doprocessMEDF) {
        histos.add("multQA/h2MultCentME", "Mixed-event anchor multiplicity vs centrality (one entry per event pair)", HistType::kTH2F, {centAxis, histoConfig.multNTracksAxis});
      }
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

      histos.add("QAevent/hRotBkg", "Rotated angle of rotated background", HistType::kTH1F, {{360, 0.0, o2::constants::math::TwoPI}});
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
      histos.add("QAbefore/trkDCAxy_xi", "DCAxy distribution of Xi- track candidates", HistType::kTH2F, {ptAxis, dcaxyAxis});

      histos.add("QAbefore/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH2F, {ptAxis, dcazAxis});
      histos.add("QAbefore/trkDCAz_xi", "DCAz distribution of Xi- track candidates", HistType::kTH2F, {ptAxis, dcazAxis});

      histos.add("QAafter/trkDCAxy_pi", "DCAxy distribution of pion track candidates", HistType::kTH2F, {ptAxis, dcaxyAxis});
      histos.add("QAafter/trkDCAxy_xi", "DCAxy distribution of Xi- track candidates", HistType::kTH2F, {ptAxis, dcaxyAxis});

      histos.add("QAafter/trkDCAz_pi", "DCAz distribution of pion track candidates", HistType::kTH2F, {ptAxis, dcazAxis});
      histos.add("QAafter/trkDCAz_xi", "DCAz distribution of Xi- track candidates", HistType::kTH2F, {ptAxis, dcazAxis});
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
    histos.add("h3Xi1530invmassDS", "Invariant mass of Xi- differnt sign", kTHnSparseD, {centAxis, ptAxis, invMassAxis, flagAxis});
    histos.add("h3XiinvmassDS", "Invariant mass of Xi- differnt sign", kTHnSparseD, {centAxis, ptAxis, invMassAxisCasc, flagAxis});

    histos.add("h3Xi1530invmassLS", "Invariant mass of Xi(1530)0 same sign", kTHnSparseD, {centAxis, ptAxis, invMassAxis, flagAxis});
    histos.add("h3Xi1530invmassRotDS", "Invariant mass of Xi(1530)0 rotated DS", kTHnSparseD, {centAxis, ptAxis, invMassAxis, flagAxis});

    histos.add("h3Xi1530invmassDSAnti", "Invariant mass of Anti-Xi(1530)0 differnt sign", kTHnSparseD, {centAxis, ptAxis, invMassAxis, flagAxis});
    histos.add("h3XiinvmassDSAnti", "Invariant mass of Anti-Xi- differnt sign", kTHnSparseD, {centAxis, ptAxis, invMassAxisCasc, flagAxis});

    histos.add("h3Xi1530invmassLSAnti", "Invariant mass of Anti-Xi(1530)0 same sign", kTHnSparseD, {centAxis, ptAxis, invMassAxis, flagAxis});
    histos.add("h3Xi1530invmassRotDSAnti", "Invariant mass of Anti-Xi(1530)0 rotated DS", kTHnSparseD, {centAxis, ptAxis, invMassAxis, flagAxis});

    if (doprocessMEMicro || doprocessMEDF) {
      histos.add("h3Xi1530invmassME_DS", "Invariant mass of Xi(1530)0 mixed event DS", kTHnSparseD, {centAxis, ptAxis, invMassAxis, flagAxis});
      histos.add("h3Xi1530invmassME_DSAnti", "Invariant mass of Xi(1530)0 mixed event DSAnti", kTHnSparseD, {centAxis, ptAxis, invMassAxis, flagAxis});
    }

    if (doprocessMC || doprocessMCMicro) {
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

      histos.add("h3RecXi1530invmass", "Invariant mass of Reconstructed MC Xi(1530)0", kTHnSparseD, {centAxis, ptAxis, invMassAxis, ptAxis});
      histos.add("h3RecXiinvmass", "Invariant mass of Reconstructed MC Xi-", kTHnSparseD, {centAxis, ptAxis, invMassAxisCasc, flagAxis});

      histos.add("h3RecXi1530invmassAnti", "Invariant mass of Reconstructed MC Anti-Xi(1530)0", kTHnSparseD, {centAxis, ptAxis, invMassAxis, ptAxis});
      histos.add("h3RecXiinvmassAnti", "Invariant mass of Reconstructed MC Anti-Xi-", kTHnSparseD, {centAxis, ptAxis, invMassAxisCasc, flagAxis});

      histos.add("Xi1530Rec", "pT distribution of Reconstructed MC Xi(1530)0", kTH2D, {ptAxis, centAxis});
      histos.add("Xi1530RecAnti", "pT distribution of Reconstructed MC Anti-Xi(1530)0", kTH2D, {ptAxis, centAxis});
      histos.add("Xi1530Recinvmass", "Inv mass distribution of Reconstructed MC Xi(1530)0", kTH1F, {invMassAxis});
    }
    if (doprocessMC || doprocessMCMicro || doprocessMCTrue) {
      histos.add("h3Xi1530Gen", "pT distribution of True MC Xi(1530)0", kTHnSparseD, {ptAxis, centAxis, histoConfig.multNTracksAxis});
      histos.add("h3Xi1530GenAnti", "pT distribution of True MC Anti-Xi(1530)0", kTHnSparseD, {ptAxis, centAxis, histoConfig.multNTracksAxis});
    }
    // QA for topological, kinematical cut for cascades
    if (histoConfig.additionalQAplots) {
      histos.add("QAbefore/V0DCATopPV", "V0s DCA to PV distribution as pt", HistType::kTH2F, {ptAxis, dcaxyAxis});
      histos.add("QAbefore/V0DCADoughter", "V0s DCA Doughter distribution as pt", HistType::kTH2F, {ptAxis, dcaDaugAxis});
      histos.add("QAbefore/CascDCADoughter", "Casc DCA Doughter distribution as pt", HistType::kTH2F, {ptAxis, dcaDaugAxis});
      histos.add("QAbefore/CascDCABachPV", "Casc DCA Bachlor to PV distribution as pt", HistType::kTH2F, {ptAxis, dcaBachAxis});
      histos.add("QAbefore/CascDCAPosPV", "Casc DCA Positive to PV distribution as pt", HistType::kTH2F, {ptAxis, dcaSecondaryAxis});
      histos.add("QAbefore/CascDCANegPV", "Casc DCA Negative to PV distribution as pt", HistType::kTH2F, {ptAxis, dcaSecondaryAxis});
      histos.add("QAbefore/V0Mass", "V0 mass distribution as pt", HistType::kTH2F, {ptAxis, invMassAxisLambda});
      histos.add("QAbefore/CascMass", "Casc mass distribution as pt", HistType::kTH2F, {ptAxis, invMassAxisCasc});
      histos.add("QAbefore/V0CosPA", "V0s CosPA distribution as pt", HistType::kTH2F, {ptAxis, cosPAAxis});
      histos.add("QAbefore/CascCosPA", "Casc CosPA distribution as pt", HistType::kTH2F, {ptAxis, cosPAAxis});
      histos.add("QAbefore/V0Radius", "V0 Radius distribution as pt", HistType::kTH2F, {ptAxis, transRadiusAxis});
      histos.add("QAbefore/CascRadius", "Casc Radius distribution as pt", HistType::kTH2F, {ptAxis, transRadiusAxis});
      histos.add("QAbefore/ProperLifetime", "Proper Lifetime distribution as pt", HistType::kTH2F, {ptAxis, properLifetimeAxis});
      histos.add("QAbefore/NCrossedRowsPos", "Number of crossed rows in TPC for positive daughter in V0s", HistType::kTH2F, {ptAxis, {200, 0, 200, "N crossed rows"}});
      histos.add("QAbefore/NCrossedRowsNeg", "Number of crossed rows in TPC for negative daughter in V0s", HistType::kTH2F, {ptAxis, {200, 0, 200, "N crossed rows"}});
      histos.add("QAbefore/NCrossedRowsBach", "Number of crossed rows in TPC for bachelor in Cascades", HistType::kTH2F, {ptAxis, {200, 0, 200, "N crossed rows"}});

      histos.add("QAafter/V0DCATopPV", "V0s DCA to PV distribution as pt", HistType::kTH2F, {ptAxis, dcaxyAxis});
      histos.add("QAafter/V0DCADoughter", "V0s DCA Doughter distribution as pt", HistType::kTH2F, {ptAxis, dcaDaugAxis});
      histos.add("QAafter/CascDCADoughter", "Casc DCA Doughter distribution as pt", HistType::kTH2F, {ptAxis, dcaDaugAxis});
      histos.add("QAafter/CascDCABachPV", "Casc DCA Bachlor to PV distribution as pt", HistType::kTH2F, {ptAxis, dcaBachAxis});
      histos.add("QAafter/CascDCAPosPV", "Casc DCA Positive to PV distribution as pt", HistType::kTH2F, {ptAxis, dcaSecondaryAxis});
      histos.add("QAafter/CascDCANegPV", "Casc DCA Negative to PV distribution as pt", HistType::kTH2F, {ptAxis, dcaSecondaryAxis});
      histos.add("QAafter/V0Mass", "V0 mass distribution as pt", HistType::kTH2F, {ptAxis, invMassAxisLambda});
      histos.add("QAafter/CascMass", "Casc mass distribution as pt", HistType::kTH2F, {ptAxis, invMassAxisCasc});
      histos.add("QAafter/V0CosPA", "V0s CosPA distribution as pt", HistType::kTH2F, {ptAxis, cosPAAxis});
      histos.add("QAafter/CascCosPA", "Casc CosPA distribution as pt", HistType::kTH2F, {ptAxis, cosPAAxis});
      histos.add("QAafter/V0Radius", "V0 Radius distribution as pt", HistType::kTH2F, {ptAxis, transRadiusAxis});
      histos.add("QAafter/CascRadius", "Casc Radius distribution as pt", HistType::kTH2F, {ptAxis, transRadiusAxis});
      histos.add("QAafter/ProperLifetime", "Proper Lifetime distribution as pt", HistType::kTH2F, {ptAxis, properLifetimeAxis});
      histos.add("QAafter/NCrossedRowsPos", "Number of crossed rows in TPC for positive daughter in V0s", HistType::kTH2F, {ptAxis, {200, 0, 200, "N crossed rows"}});
      histos.add("QAafter/NCrossedRowsNeg", "Number of crossed rows in TPC for negative daughter in V0s", HistType::kTH2F, {ptAxis, {200, 0, 200, "N crossed rows"}});
      histos.add("QAafter/NCrossedRowsBach", "Number of crossed rows in TPC for bachelor in Cascades", HistType::kTH2F, {ptAxis, {200, 0, 200, "N crossed rows"}});

      histos.add("QAMCTrue/V0DCATopPV", "V0s DCA to PV distribution as pt", HistType::kTH2F, {ptAxis, dcaxyAxis});
      histos.add("QAMCTrue/V0DCADoughter", "V0s DCA Doughter distribution as pt", HistType::kTH2F, {ptAxis, dcaDaugAxis});
      histos.add("QAMCTrue/CascDCADoughter", "Casc DCA Doughter distribution as pt", HistType::kTH2F, {ptAxis, dcaDaugAxis});
      histos.add("QAMCTrue/CascDCABachPV", "Casc DCA Bachlor to PV distribution as pt", HistType::kTH2F, {ptAxis, dcaBachAxis});
      histos.add("QAMCTrue/CascDCAPosPV", "Casc DCA Positive to PV distribution as pt", HistType::kTH2F, {ptAxis, dcaSecondaryAxis});
      histos.add("QAMCTrue/CascDCANegPV", "Casc DCA Negative to PV distribution as pt", HistType::kTH2F, {ptAxis, dcaSecondaryAxis});
      histos.add("QAMCTrue/V0Mass", "V0 mass distribution as pt", HistType::kTH2F, {ptAxis, invMassAxisLambda});
      histos.add("QAMCTrue/CascMass", "Casc mass distribution as pt", HistType::kTH2F, {ptAxis, invMassAxisCasc});
      histos.add("QAMCTrue/V0CosPA", "V0s CosPA distribution as pt", HistType::kTH2F, {ptAxis, cosPAAxis});
      histos.add("QAMCTrue/CascCosPA", "Casc CosPA distribution as pt", HistType::kTH2F, {ptAxis, cosPAAxis});
      histos.add("QAMCTrue/V0Radius", "V0 Radius distribution as pt", HistType::kTH2F, {ptAxis, transRadiusAxis});
      histos.add("QAMCTrue/CascRadius", "Casc Radius distribution as pt", HistType::kTH2F, {ptAxis, transRadiusAxis});
      histos.add("QAMCTrue/ProperLifetime", "Proper Lifetime distribution as pt", HistType::kTH2F, {ptAxis, properLifetimeAxis});
      histos.add("QAMCTrue/NCrossedRowsPos", "Number of crossed rows in TPC for positive daughter in V0s", HistType::kTH2F, {ptAxis, {200, 0, 200, "N crossed rows"}});
      histos.add("QAMCTrue/NCrossedRowsNeg", "Number of crossed rows in TPC for negative daughter in V0s", HistType::kTH2F, {ptAxis, {200, 0, 200, "N crossed rows"}});
      histos.add("QAMCTrue/NCrossedRowsBach", "Number of crossed rows in TPC for bachelor in Cascades", HistType::kTH2F, {ptAxis, {200, 0, 200, "N crossed rows"}});
    }
  }

  float massPi = MassPionCharged;

  // Primary track selection for the first pion //
  template <bool IsResoMicrotrack, typename TrackType>
  bool primaryTrackCut(const TrackType& track)
  {
    if (!std::isfinite(track.pt()) || !std::isfinite(track.eta()) ||
        !std::isfinite(track.dcaXY()) || !std::isfinite(track.dcaZ())) {
      return false;
    }
    if (std::abs(track.eta()) >= primarytrackConfig.cMaxetacut) {
      return false;
    }
    if (std::abs(track.pt()) <= primarytrackConfig.cMinPtcut) {
      return false;
    }
    // Keep double precision and Power=1 equivalent to the previous pow(pt, -1.) cut.
    // Evaluate the common threshold only when at least one axis needs it.
    const double dcaPtCut = (primarytrackConfig.cDCAxyToPVAsPt || primarytrackConfig.cDCAzToPVAsPt)
                              ? primarytrackConfig.cDCAxytoPVByPtPiFirstP0.value + primarytrackConfig.cDCAxyToPVByPtPiFirstExp.value * std::pow(track.pt(), -static_cast<double>(primarytrackConfig.cDCAxyToPVByPtPiFirstPower.value))
                              : 0.;
    const double dcaXYCut = primarytrackConfig.cDCAxyToPVAsPt ? dcaPtCut : primarytrackConfig.cMaxDCAxyToPVCut.value;
    const double dcaZCut = primarytrackConfig.cDCAzToPVAsPt ? dcaPtCut : primarytrackConfig.cMaxDCAzToPVCut.value;
    if (std::abs(track.dcaXY()) >= dcaXYCut || std::abs(track.dcaZ()) >= dcaZCut) {
      return false;
    }
    if constexpr (!IsResoMicrotrack) {
      if (track.tpcNClsFound() <= primarytrackConfig.cfgTPCcluster) {
        return false;
      }
      if (track.tpcNClsCrossedRows() <= primarytrackConfig.cfgTPCRows) {
        return false;
      }
    }
    if (primarytrackConfig.cfgPVContributor && !track.isPVContributor()) {
      return false;
    }
    if (primarytrackConfig.cfgPrimaryTrack && !track.isPrimaryTrack()) {
      return false;
    }
    return true;
  }

  bool hasSubsystemInfo(float Nsigma)
  {
    return std::abs(Nsigma) < pidConfig.cPIDBound;
  }

  // Primary track selection for cascades, Need to more informations for cascades //
  template <typename TracksTypeCasc>
  bool cascprimaryTrackCut(const TracksTypeCasc& track)
  {
    if (std::abs(track.eta()) >= primarytrackConfig.cMaxetacut) {
      return false;
    }
    if (std::abs(track.pt()) <= primarytrackConfig.cMinPtcut) {
      return false;
    }
    if (track.nCrossedRowsPos() <= cascadeConfig.cMinNCrossedRowsTPCPos) {
      return false;
    }
    if (track.nCrossedRowsNeg() <= cascadeConfig.cMinNCrossedRowsTPCNeg) {
      return false;
    }
    if (track.nCrossedRowsBach() <= cascadeConfig.cMinNCrossedRowsTPCBach) {
      return false;
    }
    if (primarytrackConfig.cDCAxyToPVAsPtForCasc) {
      if (std::abs(track.dcaXYCascToPV()) >= (primarytrackConfig.cDCAxyToPVByPtCascP0 + primarytrackConfig.cDCAxyToPVByPtCascExp * track.pt())) {
        return false;
      }
    }
    if (primarytrackConfig.cDCAzToPVAsPtForCasc) {
      if (std::abs(track.dcaZCascToPV()) >= (primarytrackConfig.cDCAxyToPVByPtCascP0 + primarytrackConfig.cDCAxyToPVByPtCascExp * std::pow(track.pt(), -1.))) {
        return false;
      }
    } else if (std::abs(track.dcaZCascToPV()) >= primarytrackConfig.cMaxDCAzToPVCutCasc) {
      return false;
    }

    return true;
  }

  // Secondary track selection for cascades //

  // Topological cuts for cascades
  template <typename TracksTypeCasc>
  bool casctopCut(const TracksTypeCasc& track)
  {
    // Topological cuts for V0s
    if (std::abs(track.daughDCA()) >= v0sConfig.cDCALambdaDaugtherscut) {
      return false;
    }
    if (std::abs(track.dcav0topv()) <= v0sConfig.cDCALambdaToPVcut) {
      return false;
    }
    if (track.sign() < 0) {
      if (std::abs(track.dcanegtopv()) <= v0sConfig.cDCAPionToPVcut) {
        return false;
      }
      if (std::abs(track.dcapostopv()) <= v0sConfig.cDCAProtonToPVcut) {
        return false;
      }
    } else {
      if (std::abs(track.dcanegtopv()) <= v0sConfig.cDCAProtonToPVcut) {
        return false;
      }
      if (std::abs(track.dcapostopv()) <= v0sConfig.cDCAPionToPVcut) {
        return false;
      }
    }
    if (track.v0CosPA() <= std::cos(v0sConfig.cV0CosPACutPtDepP0 - v0sConfig.cV0CosPACutPtDepP1 * track.pt())) {
      return false;
    }
    if (track.transRadius() >= v0sConfig.cMaxV0radiuscut || track.transRadius() <= v0sConfig.cMinV0radiuscut) {
      return false;
    }
    if (std::abs(track.mLambda() - MassLambda) >= v0sConfig.cMasswindowV0cut) {
      return false;
    }

    // Topological Cuts for Cascades
    if (std::abs(track.dcabachtopv()) <= cascadeConfig.cDCABachlorToPVcut) {
      return false;
    }
    if (track.pt() <= cascadeConfig.cDCAXiDaugthersCutPtRangeLower) {
      if (track.cascDaughDCA() >= cascadeConfig.cDCAXiDaugthersCutPtDepLower) {
        return false;
      }
    }
    if (track.pt() >= cascadeConfig.cDCAXiDaugthersCutPtRangeLower && track.pt() <= cascadeConfig.cDCAXiDaugthersCutPtRangeUpper) {
      if (track.cascDaughDCA() >= cascadeConfig.cDCAXiDaugthersCutPtDepMiddle) {
        return false;
      }
    }
    if (track.pt() >= cascadeConfig.cDCAXiDaugthersCutPtRangeUpper) {
      if (track.cascDaughDCA() >= cascadeConfig.cDCAXiDaugthersCutPtDepUpper) {
        return false;
      }
    }
    if (track.cascCosPA() <= std::cos(cascadeConfig.cCosPACascCutPtDepP0 - cascadeConfig.cCosPACascCutPtDepP1 * track.pt())) {
      return false;
    }

    if (track.cascTransRadius() >= cascadeConfig.cMaxCascradiuscut || track.cascTransRadius() <= cascadeConfig.cMinCascradiuscut) {
      return false;
    }
    if (std::abs(track.mXi() - cascadeConfig.cMassXiminus) >= cascadeConfig.cMasswindowCasccut) {
      return false;
    }

    return true;
  }

  bool pidSelector(float TPCNsigma, float TOFNsigma, const PidSelectionParam& params, bool tofAtHighPt, float trackPt, bool hasTOF)
  {
    // Bounds are signed offsets from the detector-specific mean; exclude both edges.
    const float centeredTPC = TPCNsigma - params.cMeanTPCnSigma;
    const bool passesTPCWindow = params.cMinTPCnSigma < centeredTPC && centeredTPC < params.cMaxTPCnSigma;
    if (!passesTPCWindow) {
      return false;
    }
    // Missing or bypassed TOF never vetoes a track that passes the TPC window.
    if (params.cByPassTOF || (tofAtHighPt && trackPt <= pidConfig.cMinTOFpt) || !hasTOF) {
      return true;
    }
    // Retain TOF-window OR combined acceptance; center both detectors for the latter.
    // A present Micro001 TOF value at +/-infinity is overflow, not missing TOF.
    const float centeredTOF = TOFNsigma - params.cMeanTOFnSigma;
    return (params.cMinTOFnSigma < centeredTOF && centeredTOF < params.cMaxTOFnSigma) ||
           (params.nsigmaCutCombined > 0 &&
            centeredTPC * centeredTPC + centeredTOF * centeredTOF < params.nsigmaCutCombined * params.nsigmaCutCombined);
  }

  // PID selection for the First Pion //
  template <typename T>
  bool selectionPIDPionFirst(const T& candidate)
  {

    const float tpcNsigmaPionFirst = candidate.tpcNSigmaPi();
    const float tofNsigmaPionFirst = candidate.tofNSigmaPi();
    float trackPt = candidate.pt();

    PidSelectionParam pionFirstParams = {.cMinTPCnSigma = pidConfig.cMintpcnSigmaPionFirst, .cMaxTPCnSigma = pidConfig.cMaxtpcnSigmaPionFirst, .cMeanTPCnSigma = pidConfig.cMeantpcnSigmaPionFirst, .cMinTOFnSigma = pidConfig.cMintofnSigmaPionFirst, .cMaxTOFnSigma = pidConfig.cMaxtofnSigmaPionFirst, .cMeanTOFnSigma = pidConfig.cMeantofnSigmaPionFirst, .cByPassTOF = pidConfig.cByPassTOFPionFirst, .nsigmaCutCombined = pidConfig.nsigmaCutCombinedPionFirst};

    return pidSelector(tpcNsigmaPionFirst, tofNsigmaPionFirst, pionFirstParams, pidConfig.tofAtHighPt, trackPt, candidate.hasTOF());
  }

  template <typename TCascade>
  bool selectionPIDCascades(const TCascade& candidate)
  {
    bool lConsistentWithXi{false}, lConsistentWithLambda{false}, lConsistentWithPion{false}, lConsistentWithProton{false};

    float tpcNsigmaBachelor = 0.f, tofNsigmaBachelor = 0.f;
    float tpcNsigmaPion = 0.f, tofNsigmaPion = 0.f;
    float tpcNsigmaProton = 0.f, tofNsigmaProton = 0.f;
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

    PidSelectionParam bachelorParams = {.cMinTPCnSigma = pidConfig.cMintpcnSigmaPionBachelor, .cMaxTPCnSigma = pidConfig.cMaxtpcnSigmaPionBachelor, .cMeanTPCnSigma = pidConfig.cMeantpcnSigmaPionBachelor, .cMinTOFnSigma = pidConfig.cMintofnSigmaPionBachelor, .cMaxTOFnSigma = pidConfig.cMaxtofnSigmaPionBachelor, .cMeanTOFnSigma = pidConfig.cMeantofnSigmaPionBachelor, .cByPassTOF = pidConfig.cByPassTOFPionBachelor, .nsigmaCutCombined = pidConfig.nsigmaCutCombinedPionBachelor};
    PidSelectionParam pionParams = {.cMinTPCnSigma = pidConfig.cMintpcnSigmaPion, .cMaxTPCnSigma = pidConfig.cMaxtpcnSigmaPion, .cMeanTPCnSigma = pidConfig.cMeantpcnSigmaPion, .cMinTOFnSigma = pidConfig.cMintofnSigmaPion, .cMaxTOFnSigma = pidConfig.cMaxtofnSigmaPion, .cMeanTOFnSigma = pidConfig.cMeantofnSigmaPion, .cByPassTOF = pidConfig.cByPassTOFPion, .nsigmaCutCombined = pidConfig.nsigmaCutCombinedPion};
    PidSelectionParam protonParams = {.cMinTPCnSigma = pidConfig.cMintpcnSigmaProton, .cMaxTPCnSigma = pidConfig.cMaxtpcnSigmaProton, .cMeanTPCnSigma = pidConfig.cMeantpcnSigmaProton, .cMinTOFnSigma = pidConfig.cMintofnSigmaProton, .cMaxTOFnSigma = pidConfig.cMaxtofnSigmaProton, .cMeanTOFnSigma = pidConfig.cMeantofnSigmaProton, .cByPassTOF = pidConfig.cByPassTOFProton, .nsigmaCutCombined = pidConfig.nsigmaCutCombinedProton};

    lConsistentWithXi = pidSelector(tpcNsigmaBachelor, tofNsigmaBachelor, bachelorParams, pidConfig.tofAtHighPt, trackPt, hasSubsystemInfo(tofNsigmaBachelor));
    lConsistentWithPion = pidSelector(tpcNsigmaPion, tofNsigmaPion, pionParams, pidConfig.tofAtHighPt, trackPt, hasSubsystemInfo(tofNsigmaPion));
    lConsistentWithProton = pidSelector(tpcNsigmaProton, tofNsigmaProton, protonParams, pidConfig.tofAtHighPt, trackPt, hasSubsystemInfo(tofNsigmaProton));

    lConsistentWithLambda = lConsistentWithProton && lConsistentWithPion;

    return lConsistentWithXi && lConsistentWithLambda;
  }

  template <typename CollisionType, typename TCascade>
  float properLifetime(const CollisionType& collision, const TCascade& candidate)
  {
    float kSmallMomentumDenominator = 1e-6f; // To avoid division by zero, if momentum is extremely small, we consider it as 1e-6 GeV/c
    float dx = candidate.decayVtxX() - collision.posX();
    float dy = candidate.decayVtxY() - collision.posY();
    float dz = candidate.decayVtxZ() - collision.posZ();
    float l = std::sqrt(dx * dx + dy * dy + dz * dz);
    float p = std::sqrt(candidate.px() * candidate.px() + candidate.py() * candidate.py() + candidate.pz() * candidate.pz());
    auto properLifetime = (l / (p + kSmallMomentumDenominator)) * candidate.mXi();
    return properLifetime;
  }

  template <bool IsResoMicrotrack, bool IsMC, bool IsMix, typename CollisionType, typename centType, typename TracksType, typename TracksTypeCasc>
  void fillHistograms(const CollisionType& collision, const centType& inCent, const TracksType& dTracks1, const TracksTypeCasc& dTracks2,
                      const MixingCollision* cascadeCollision = nullptr) // Order: ResoColl, ResoTrack, ResoCascTrack
  {
    auto cent = inCent;

    if (histoConfig.eventQA) {
      if constexpr (!IsMix) {
        histos.fill(HIST("QAevent/hVertexZSameE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentSameE"), cent);
        histos.fill(HIST("QAevent/hCollisionIndexSameE"), collision.globalIndex());
        histos.fill(HIST("QAevent/hnTrksSameE"), dTracks1.size());
        histos.fill(HIST("QAevent/hnCascsSameE"), dTracks2.size());
      } else {
        histos.fill(HIST("QAevent/hVertexZMixedE"), collision.posZ());
        histos.fill(HIST("QAevent/hMultiplicityPercentMixedE"), cent);
        histos.fill(HIST("QAevent/hCollisionIndexMixedE"), collision.globalIndex());
        histos.fill(HIST("QAevent/hnTrksMixedE"), dTracks1.size());
        histos.fill(HIST("QAevent/hnCascsMixedE"), dTracks2.size());
      }
    }

    LorentzVectorPtEtaPhiMass lDecayDaughter1, lDecayDaughter2, lResonance, lDaughterRot, lResonanceRot;
    std::vector<int64_t> pionCandateIndicies = {};
    std::vector<int64_t> xiCandateIndicies = {};
    pionCandateIndicies.reserve(dTracks1.size());
    xiCandateIndicies.reserve(dTracks2.size());

    for (const auto& trk1 : dTracks1) {
      auto trk1ptPi = trk1.pt();
      const float trk1DCAXY = trk1.dcaXY();
      const float trk1DCAZ = trk1.dcaZ();
      const float trk1NSigmaPiTPC = trk1.tpcNSigmaPi();
      const float trk1NSigmaPiTOF = trk1.tofNSigmaPi();

      // QA before
      if constexpr (!IsMix) {
        if (histoConfig.pidPlots) {
          histos.fill(HIST("QAbefore/TPC_Nsigma_pi_first_all"), cent, trk1ptPi, trk1NSigmaPiTPC);
          if (trk1.hasTOF() && !std::isnan(trk1NSigmaPiTOF)) {
            histos.fill(HIST("QAbefore/TOF_Nsigma_pi_first_all"), cent, trk1ptPi, trk1NSigmaPiTOF);
            histos.fill(HIST("QAbefore/TOF_TPC_Map_pi_first_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
          }
        }
        if (histoConfig.additionalQAplots) {
          histos.fill(HIST("QAbefore/trkDCAxy_pi"), trk1ptPi, std::abs(trk1DCAXY));
          histos.fill(HIST("QAbefore/trkDCAz_pi"), trk1ptPi, std::abs(trk1DCAZ));
        }
      }
      if (!selectionPIDPionFirst(trk1)) { // PID selection for the first pion
        continue;
      }
      if (!primaryTrackCut<IsResoMicrotrack>(trk1)) { // Primary track selections
        continue;
      }

      if constexpr (!IsMix) {
        if (histoConfig.pidPlots) {
          histos.fill(HIST("QAafter/TPC_Nsigma_pi_first_all"), cent, trk1ptPi, trk1NSigmaPiTPC);

          if (trk1.hasTOF() && !std::isnan(trk1NSigmaPiTOF)) {
            histos.fill(HIST("QAafter/TOF_Nsigma_pi_first_all"), cent, trk1ptPi, trk1NSigmaPiTOF);
            histos.fill(HIST("QAafter/TOF_TPC_Map_pi_first_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
          }
        }
        if (histoConfig.additionalQAplots) {
          histos.fill(HIST("QAafter/trkDCAxy_pi"), trk1ptPi, std::abs(trk1DCAXY));
          histos.fill(HIST("QAafter/trkDCAz_pi"), trk1ptPi, std::abs(trk1DCAZ));
        }
      }

      pionCandateIndicies.push_back(trk1.index());
    }

    for (const auto& trk2 : dTracks2) {
      // Kinematic variables for cascades
      auto trk2ptXi = trk2.pt();
      auto massLambdaCand = trk2.mLambda();
      auto massXiCand = trk2.mXi();
      auto trk2ProperLifetime = cascadeCollision ? properLifetime(*cascadeCollision, trk2) : properLifetime(collision, trk2);

      // Topological variables for cascades
      auto trk2DCAV0TopPV = trk2.dcav0topv();
      auto trk2DCAXY = trk2.dcaXYCascToPV();
      auto trk2DCAZ = trk2.dcaZCascToPV();
      auto trk2DCAV0sDougthers = trk2.daughDCA(); // V0s daughter DCA
      auto trk2DCACascDougthers = trk2.cascDaughDCA();
      auto trk2DCABachPV = trk2.dcabachtopv();
      auto trk2DCAPosPV = trk2.dcapostopv();
      auto trk2DCANegPV = trk2.dcanegtopv();
      auto trk2V0CosPA = trk2.v0CosPA();
      auto trk2CascCosPA = trk2.cascCosPA();
      auto trk2V0Radius = trk2.transRadius();
      auto trk2CascRadius = trk2.cascTransRadius();
      auto trks2NCrossedRowsPos = trk2.nCrossedRowsPos();
      auto trks2NCrossedRowsNeg = trk2.nCrossedRowsNeg();
      auto trks2NCrossedRowsBach = trk2.nCrossedRowsBach();

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
          histos.fill(HIST("QAbefore/TPC_Nsigma_pi_bachelor_all"), cent, 0, trk2NSigmaPiBachelorTPC); // can't take pt information for the cascade secondary
          if (hasSubsystemInfo(trk2NSigmaPiBachelorTOF)) {
            histos.fill(HIST("QAbefore/TOF_TPC_Map_pi_bachelor_all"), trk2NSigmaPiBachelorTOF, trk2NSigmaPiBachelorTPC);
          }

          histos.fill(HIST("QAbefore/TPC_Nsigma_pr_all"), cent, 0, trk2NSigmaPrPosTPC);
          if (hasSubsystemInfo(trk2NSigmaPrPosTOF)) {
            histos.fill(HIST("QAbefore/TOF_TPC_Map_pr_all"), trk2NSigmaPrPosTOF, trk2NSigmaPrPosTPC);
          }

          histos.fill(HIST("QAbefore/TPC_Nsigma_antipr_all"), cent, 0, trk2NSigmaPrNegTPC);
          if (hasSubsystemInfo(trk2NSigmaPrNegTOF)) {
            histos.fill(HIST("QAbefore/TOF_TPC_Map_antipr_all"), trk2NSigmaPrNegTOF, trk2NSigmaPrNegTPC);
          }

          histos.fill(HIST("QAbefore/TPC_Nsigma_pi_all"), cent, 0, trk2NSigmaPiPosTPC);
          if (hasSubsystemInfo(trk2NSigmaPiPosTOF)) {
            histos.fill(HIST("QAbefore/TOF_TPC_Map_pi_all"), trk2NSigmaPiPosTOF, trk2NSigmaPiPosTPC);
          }

          histos.fill(HIST("QAbefore/TPC_Nsigma_piminus_all"), cent, 0, trk2NSigmaPiNegTPC);
          if (hasSubsystemInfo(trk2NSigmaPiNegTOF)) {
            histos.fill(HIST("QAbefore/TOF_TPC_Map_piminus_all"), trk2NSigmaPiNegTOF, trk2NSigmaPiNegTPC);
          }
        }

        if (histoConfig.additionalQAplots) {
          histos.fill(HIST("QAbefore/V0DCATopPV"), trk2ptXi, trk2DCAV0TopPV);
          histos.fill(HIST("QAbefore/trkDCAxy_xi"), trk2ptXi, std::abs(trk2DCAXY));
          histos.fill(HIST("QAbefore/trkDCAz_xi"), trk2ptXi, std::abs(trk2DCAZ));
          histos.fill(HIST("QAbefore/V0DCADoughter"), trk2ptXi, trk2DCAV0sDougthers);
          histos.fill(HIST("QAbefore/CascDCADoughter"), trk2ptXi, trk2DCACascDougthers);
          histos.fill(HIST("QAbefore/CascDCABachPV"), trk2ptXi, trk2DCABachPV);
          histos.fill(HIST("QAbefore/CascDCAPosPV"), trk2ptXi, trk2DCAPosPV);
          histos.fill(HIST("QAbefore/CascDCANegPV"), trk2ptXi, trk2DCANegPV);
          histos.fill(HIST("QAbefore/V0CosPA"), trk2ptXi, 1. - trk2V0CosPA);
          histos.fill(HIST("QAbefore/CascCosPA"), trk2ptXi, 1. - trk2CascCosPA);
          histos.fill(HIST("QAbefore/V0Radius"), trk2ptXi, trk2V0Radius);
          histos.fill(HIST("QAbefore/CascRadius"), trk2ptXi, trk2CascRadius);
          histos.fill(HIST("QAbefore/V0Mass"), trk2ptXi, massLambdaCand);
          histos.fill(HIST("QAbefore/CascMass"), trk2ptXi, massXiCand);
          histos.fill(HIST("QAbefore/ProperLifetime"), trk2ptXi, trk2ProperLifetime);
          histos.fill(HIST("QAbefore/NCrossedRowsPos"), trk2ptXi, trks2NCrossedRowsPos);
          histos.fill(HIST("QAbefore/NCrossedRowsNeg"), trk2ptXi, trks2NCrossedRowsNeg);
          histos.fill(HIST("QAbefore/NCrossedRowsBach"), trk2ptXi, trks2NCrossedRowsBach);
        }
      }

      if (!selectionPIDCascades(trk2)) {
        continue;
      }
      if (!cascprimaryTrackCut(trk2) || !casctopCut(trk2)) { // Primary track selections
        continue;
      }
      if (trk2ProperLifetime >= cascadeConfig.cMaxProperLifetimeCut) {
        continue;
      }

      // QA after selections
      if constexpr (!IsMix) {
        //// QA plots after the selection
        //  --- PID QA
        if (histoConfig.pidPlots) {

          if (trk2.sign() < 0) {
            histos.fill(HIST("QAafter/TPC_Nsigma_pi_bachelor_all"), cent, 0, trk2NSigmaPiBachelorTPC); // not exist pt information in resocascade yet.
            if (hasSubsystemInfo(trk2NSigmaPiBachelorTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_pi_bachelor_all"), trk2NSigmaPiBachelorTOF, trk2NSigmaPiBachelorTPC);
            }

            histos.fill(HIST("QAafter/TPC_Nsigma_pr_all"), cent, 0, trk2NSigmaPrPosTPC);
            if (hasSubsystemInfo(trk2NSigmaPrPosTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_pr_all"), trk2NSigmaPrPosTOF, trk2NSigmaPrPosTPC);
            }

            histos.fill(HIST("QAafter/TPC_Nsigma_piminus_all"), cent, 0, trk2NSigmaPiNegTPC);
            if (hasSubsystemInfo(trk2NSigmaPiNegTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_piminus_all"), trk2NSigmaPiNegTOF, trk2NSigmaPiNegTPC);
            }

          } else {

            histos.fill(HIST("QAafter/TPC_Nsigma_pi_bachelor_all"), cent, 0, trk2NSigmaPiBachelorTPC); // not exist pt information in resocascade yet.
            if (hasSubsystemInfo(trk2NSigmaPiBachelorTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_pi_bachelor_all"), trk2NSigmaPiBachelorTOF, trk2NSigmaPiBachelorTPC);
            }

            histos.fill(HIST("QAafter/TPC_Nsigma_antipr_all"), cent, 0, trk2NSigmaPrNegTPC);
            if (hasSubsystemInfo(trk2NSigmaPrNegTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_antipr_all"), trk2NSigmaPrNegTOF, trk2NSigmaPrNegTPC);
            }

            histos.fill(HIST("QAafter/TPC_Nsigma_pi_all"), cent, 0, trk2NSigmaPiPosTPC);
            if (hasSubsystemInfo(trk2NSigmaPiPosTOF)) {
              histos.fill(HIST("QAafter/TOF_TPC_Map_pi_all"), trk2NSigmaPiPosTOF, trk2NSigmaPiPosTPC);
            }
          }
        }
        if (histoConfig.additionalQAplots) {
          histos.fill(HIST("QAafter/V0DCATopPV"), trk2ptXi, trk2DCAV0TopPV);
          histos.fill(HIST("QAafter/trkDCAxy_xi"), trk2ptXi, std::abs(trk2DCAXY));
          histos.fill(HIST("QAafter/trkDCAz_xi"), trk2ptXi, std::abs(trk2DCAZ));
          histos.fill(HIST("QAafter/V0DCADoughter"), trk2ptXi, trk2DCAV0sDougthers);
          histos.fill(HIST("QAafter/CascDCADoughter"), trk2ptXi, trk2DCACascDougthers);
          histos.fill(HIST("QAafter/CascDCABachPV"), trk2ptXi, trk2DCABachPV);
          histos.fill(HIST("QAafter/CascDCAPosPV"), trk2ptXi, trk2DCAPosPV);
          histos.fill(HIST("QAafter/CascDCANegPV"), trk2ptXi, trk2DCANegPV);
          histos.fill(HIST("QAafter/V0CosPA"), trk2ptXi, 1. - trk2V0CosPA);
          histos.fill(HIST("QAafter/CascCosPA"), trk2ptXi, 1. - trk2CascCosPA);
          histos.fill(HIST("QAafter/V0Radius"), trk2ptXi, trk2V0Radius);
          histos.fill(HIST("QAafter/CascRadius"), trk2ptXi, trk2CascRadius);
          histos.fill(HIST("QAafter/V0Mass"), trk2ptXi, massLambdaCand);
          histos.fill(HIST("QAafter/CascMass"), trk2ptXi, massXiCand);
          histos.fill(HIST("QAafter/ProperLifetime"), trk2ptXi, trk2ProperLifetime);
          histos.fill(HIST("QAafter/NCrossedRowsPos"), trk2ptXi, trks2NCrossedRowsPos);
          histos.fill(HIST("QAafter/NCrossedRowsNeg"), trk2ptXi, trks2NCrossedRowsNeg);
          histos.fill(HIST("QAafter/NCrossedRowsBach"), trk2ptXi, trks2NCrossedRowsBach);
        }

        if (additionalConfig.studyStableXi) {
          if (trk2.sign() < 0) {
            histos.fill(HIST("h3XiinvmassDS"), cent, trk2ptXi, massXiCand, kData);
          } else if (trk2.sign() > 0) {
            histos.fill(HIST("h3XiinvmassDSAnti"), cent, trk2ptXi, massXiCand, kData);
          }
          if constexpr (IsMC) {
            if (trk2.pdgCode() == kXiMinus) {
              histos.fill(HIST("h3RecXiinvmass"), cent, trk2ptXi, massXiCand, kMCReco);
            } else if (trk2.pdgCode() == -kXiMinus) {
              histos.fill(HIST("h3RecXiinvmassAnti"), cent, trk2ptXi, massXiCand, kMCReco);
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
        if constexpr (!IsMix) {
          const auto trackId = pionCand.trackId();
          const auto& daughterIds = xiCand.cascadeIndices();
          if (trackId >= 0 && (trackId == daughterIds[0] || trackId == daughterIds[1] || trackId == daughterIds[2])) {
            continue;
          }
        }
        auto pionCandPt = pionCand.pt();
        auto xiCandPt = xiCand.pt();
        float massXiCand = xiCand.mXi();
        if (additionalConfig.cUseFixedMassXi) {
          massXiCand = cascadeConfig.cMassXiminus;
        }

        lDecayDaughter1 = LorentzVectorPtEtaPhiMass(pionCandPt, pionCand.eta(), pionCand.phi(), massPi);
        lDecayDaughter2 = LorentzVectorPtEtaPhiMass(xiCandPt, xiCand.eta(), xiCand.phi(), massXiCand);

        lResonance = lDecayDaughter1 + lDecayDaughter2;
        auto lResonanceMass = lResonance.M();
        auto lResonancePt = lResonance.Pt();

        if ((lResonance.Rapidity() <= primarytrackConfig.cfgRapidityMinCut) || (lResonance.Rapidity() >= primarytrackConfig.cfgRapidityMaxCut)) {
          continue;
        }

        if (additionalConfig.cfgCutsOnMother) {
          if (lResonancePt >= additionalConfig.cMaxPtMotherCut) { // excluding candidates in overflow
            continue;
          }
          if (lResonanceMass >= additionalConfig.cMaxMinvMotherCut) { // excluding candidates in overflow
            continue;
          }
        }

        if (pionCand.sign() * xiCand.sign() < 0) { // Signal candidates

          if constexpr (!IsMix) {
            if (pionCand.sign() > 0) {
              if (histoConfig.invMass1D) {
                histos.fill(HIST("Xi1530invmassDS"), lResonanceMass);
              }
              histos.fill(HIST("h3Xi1530invmassDS"), cent, lResonancePt, lResonanceMass, kData);
              if (additionalConfig.cfgFillRotBkg) {
                for (int i = 0; i < additionalConfig.cfgNrotBkg; i++) {
                  const auto lRotAngle = additionalConfig.cfgNrotBkg == 1
                                           ? 0.5f * (additionalConfig.cfgMinRot + additionalConfig.cfgMaxRot)
                                           : additionalConfig.cfgMinRot + i * ((additionalConfig.cfgMaxRot - additionalConfig.cfgMinRot) / (additionalConfig.cfgNrotBkg - 1));
                  if (histoConfig.eventQA) {
                    histos.fill(HIST("QAevent/hRotBkg"), lRotAngle);
                  }
                  if (additionalConfig.cfgRotPion) {
                    lDaughterRot = lDecayDaughter1;
                    ROOT::Math::RotationZ rot(lRotAngle);
                    auto p3 = rot * lDaughterRot.Vect();
                    lDaughterRot = LorentzVectorSetXYZM(p3.X(), p3.Y(), p3.Z(), lDaughterRot.M());
                    lResonanceRot = lDaughterRot + lDecayDaughter2;
                  } else {
                    lDaughterRot = lDecayDaughter2;
                    ROOT::Math::RotationZ rot(lRotAngle);
                    auto p3 = rot * lDaughterRot.Vect();
                    lDaughterRot = LorentzVectorSetXYZM(p3.X(), p3.Y(), p3.Z(), lDaughterRot.M());
                    lResonanceRot = lDecayDaughter1 + lDaughterRot;
                  }
                  if (additionalConfig.cfgCutsOnMother &&
                      (lResonanceRot.Pt() >= additionalConfig.cMaxPtMotherCut || lResonanceRot.M() >= additionalConfig.cMaxMinvMotherCut)) {
                    continue;
                  }
                  histos.fill(HIST("h3Xi1530invmassRotDS"), cent, lResonanceRot.Pt(), lResonanceRot.M(), kData);
                }
              }

            } else if (pionCand.sign() < 0) {
              if (histoConfig.invMass1D) {
                histos.fill(HIST("Xi1530invmassDSAnti"), lResonanceMass);
              }
              histos.fill(HIST("h3Xi1530invmassDSAnti"), cent, lResonancePt, lResonanceMass, kData);
              if (additionalConfig.cfgFillRotBkg) {
                for (int i = 0; i < additionalConfig.cfgNrotBkg; i++) {
                  const auto lRotAngle = additionalConfig.cfgNrotBkg == 1
                                           ? 0.5f * (additionalConfig.cfgMinRot + additionalConfig.cfgMaxRot)
                                           : additionalConfig.cfgMinRot + i * ((additionalConfig.cfgMaxRot - additionalConfig.cfgMinRot) / (additionalConfig.cfgNrotBkg - 1));
                  if (histoConfig.eventQA) {
                    histos.fill(HIST("QAevent/hRotBkg"), lRotAngle);
                  }
                  if (additionalConfig.cfgRotPion) {
                    lDaughterRot = lDecayDaughter1;
                    ROOT::Math::RotationZ rot(lRotAngle);
                    auto p3 = rot * lDaughterRot.Vect();
                    lDaughterRot = LorentzVectorSetXYZM(p3.X(), p3.Y(), p3.Z(), lDaughterRot.M());
                    lResonanceRot = lDaughterRot + lDecayDaughter2;
                  } else {
                    lDaughterRot = lDecayDaughter2;
                    ROOT::Math::RotationZ rot(lRotAngle);
                    auto p3 = rot * lDaughterRot.Vect();
                    lDaughterRot = LorentzVectorSetXYZM(p3.X(), p3.Y(), p3.Z(), lDaughterRot.M());
                    lResonanceRot = lDecayDaughter1 + lDaughterRot;
                  }
                  if (additionalConfig.cfgCutsOnMother &&
                      (lResonanceRot.Pt() >= additionalConfig.cMaxPtMotherCut || lResonanceRot.M() >= additionalConfig.cMaxMinvMotherCut)) {
                    continue;
                  }
                  histos.fill(HIST("h3Xi1530invmassRotDSAnti"), cent, lResonanceRot.Pt(), lResonanceRot.M(), kData);
                }
              }
            }
          } else {
            if (pionCand.sign() > 0) {
              histos.fill(HIST("h3Xi1530invmassME_DS"), cent, lResonancePt, lResonanceMass, kData);
            } else if (pionCand.sign() < 0) {
              histos.fill(HIST("h3Xi1530invmassME_DSAnti"), cent, lResonancePt, lResonanceMass, kData);
            }
          }

          if constexpr (IsMC) {
            if (std::abs(xiCand.motherPDG()) != kXiStar || pionCand.motherPDG() != xiCand.motherPDG()) {
              continue;
            }
            const int motherSign = xiCand.motherPDG() > 0 ? 1 : -1;
            if (pionCand.pdgCode() != motherSign * kPiPlus || xiCand.pdgCode() != motherSign * kXiMinus) {
              continue;
            }
            if (xiCand.motherId() < 0 || pionCand.motherId() != xiCand.motherId()) {
              continue;
            }
            auto lResonancePtMC = xiCand.motherPt();
            if (additionalConfig.cUseTruthRapidity) {
              const bool passesTruthRapidity = xiCand.motherRap() > primarytrackConfig.cfgRapidityMinCut && xiCand.motherRap() < primarytrackConfig.cfgRapidityMaxCut;
              if (!passesTruthRapidity) {
                continue;
              }
            }
            if (histoConfig.truthQA) {
              const float trk1DCAXY = pionCand.dcaXY();
              const float trk1DCAZ = pionCand.dcaZ();
              const float trk1NSigmaPiTPC = pionCand.tpcNSigmaPi();
              const float trk1NSigmaPiTOF = pionCand.tofNSigmaPi();

              auto trk2DCAXY = xiCand.dcaXYCascToPV();
              auto trk2DCAZ = xiCand.dcaZCascToPV();
              auto trk2ProperLifetime = properLifetime(collision, xiCand);
              auto massLambdaCand = xiCand.mLambda();

              auto trk2DCAV0TopPV = xiCand.dcav0topv();
              auto trk2DCAV0sDougthers = xiCand.daughDCA(); // V0s daughter DCA
              auto trk2DCACascDougthers = xiCand.cascDaughDCA();
              auto trk2DCABachPV = xiCand.dcabachtopv();
              auto trk2DCAPosPV = xiCand.dcapostopv();
              auto trk2DCANegPV = xiCand.dcanegtopv();
              auto trk2V0CosPA = xiCand.v0CosPA();
              auto trk2CascCosPA = xiCand.cascCosPA();
              auto trk2V0Radius = xiCand.transRadius();
              auto trk2CascRadius = xiCand.cascTransRadius();
              auto trks2NCrossedRowsPos = xiCand.nCrossedRowsPos();
              auto trks2NCrossedRowsNeg = xiCand.nCrossedRowsNeg();
              auto trks2NCrossedRowsBach = xiCand.nCrossedRowsBach();

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
                histos.fill(HIST("QAMCTrue/trkDCAxy_pi"), pionCandPt, std::abs(trk1DCAXY));
                histos.fill(HIST("QAMCTrue/trkDCAz_pi"), pionCandPt, std::abs(trk1DCAZ));
                histos.fill(HIST("QAMCTrue/V0DCATopPV"), xiCandPt, trk2DCAV0TopPV);
                histos.fill(HIST("QAMCTrue/trkDCAxy_xi"), xiCandPt, std::abs(trk2DCAXY));
                histos.fill(HIST("QAMCTrue/trkDCAz_xi"), xiCandPt, std::abs(trk2DCAZ));

                histos.fill(HIST("QAMCTrue/V0DCADoughter"), xiCandPt, trk2DCAV0sDougthers);
                histos.fill(HIST("QAMCTrue/CascDCADoughter"), xiCandPt, trk2DCACascDougthers);
                histos.fill(HIST("QAMCTrue/CascDCABachPV"), xiCandPt, trk2DCABachPV);
                histos.fill(HIST("QAMCTrue/CascDCAPosPV"), xiCandPt, trk2DCAPosPV);
                histos.fill(HIST("QAMCTrue/CascDCANegPV"), xiCandPt, trk2DCANegPV);
                histos.fill(HIST("QAMCTrue/V0CosPA"), xiCandPt, 1. - trk2V0CosPA);
                histos.fill(HIST("QAMCTrue/CascCosPA"), xiCandPt, 1. - trk2CascCosPA);
                histos.fill(HIST("QAMCTrue/V0Radius"), xiCandPt, trk2V0Radius);
                histos.fill(HIST("QAMCTrue/CascRadius"), xiCandPt, trk2CascRadius);
                histos.fill(HIST("QAMCTrue/V0Mass"), xiCandPt, massLambdaCand);
                histos.fill(HIST("QAMCTrue/CascMass"), xiCandPt, xiCand.mXi());
                histos.fill(HIST("QAMCTrue/ProperLifetime"), xiCandPt, trk2ProperLifetime);
                histos.fill(HIST("QAMCTrue/NCrossedRowsPos"), xiCandPt, trks2NCrossedRowsPos);
                histos.fill(HIST("QAMCTrue/NCrossedRowsNeg"), xiCandPt, trks2NCrossedRowsNeg);
                histos.fill(HIST("QAMCTrue/NCrossedRowsBach"), xiCandPt, trks2NCrossedRowsBach);
              }

              if (histoConfig.pidPlots) {
                histos.fill(HIST("QAMCTrue/TPC_Nsigma_pi_first_all"), cent, pionCandPt, trk1NSigmaPiTPC);
                if (pionCand.hasTOF() && !std::isnan(trk1NSigmaPiTOF)) {
                  histos.fill(HIST("QAMCTrue/TOF_Nsigma_pi_first_all"), cent, pionCandPt, trk1NSigmaPiTOF);
                  histos.fill(HIST("QAMCTrue/TOF_TPC_Map_pi_first_all"), trk1NSigmaPiTOF, trk1NSigmaPiTPC);
                }

                histos.fill(HIST("QAMCTrue/TPC_Nsigma_pi_bachelor_all"), cent, 0, trk2NSigmaPiBachelorTPC); // not exist pt information in resocascade yet.
                if (hasSubsystemInfo(trk2NSigmaPiBachelorTOF)) {
                  histos.fill(HIST("QAMCTrue/TOF_TPC_Map_pi_bachelor_all"), trk2NSigmaPiBachelorTOF, trk2NSigmaPiBachelorTPC);
                }

                if (xiCand.sign() < 0) {
                  histos.fill(HIST("QAMCTrue/TPC_Nsigma_pr_all"), cent, 0, trk2NSigmaPrPosTPC);
                  if (hasSubsystemInfo(trk2NSigmaPrPosTOF)) {
                    histos.fill(HIST("QAMCTrue/TOF_TPC_Map_pr_all"), trk2NSigmaPrPosTOF, trk2NSigmaPrPosTPC);
                  }
                  histos.fill(HIST("QAMCTrue/TPC_Nsigma_piminus_all"), cent, 0, trk2NSigmaPiNegTPC);
                  if (hasSubsystemInfo(trk2NSigmaPiNegTOF)) {
                    histos.fill(HIST("QAMCTrue/TOF_TPC_Map_piminus_all"), trk2NSigmaPiNegTOF, trk2NSigmaPiNegTPC);
                  }
                } else {
                  histos.fill(HIST("QAMCTrue/TPC_Nsigma_antipr_all"), cent, 0, trk2NSigmaPrNegTPC);
                  if (hasSubsystemInfo(trk2NSigmaPrNegTOF)) {
                    histos.fill(HIST("QAMCTrue/TOF_TPC_Map_antipr_all"), trk2NSigmaPrNegTOF, trk2NSigmaPrNegTPC);
                  }
                  histos.fill(HIST("QAMCTrue/TPC_Nsigma_pi_all"), cent, 0, trk2NSigmaPiPosTPC);
                  if (hasSubsystemInfo(trk2NSigmaPiPosTOF)) {
                    histos.fill(HIST("QAMCTrue/TOF_TPC_Map_pi_all"), trk2NSigmaPiPosTOF, trk2NSigmaPiPosTPC);
                  }
                }
              }

            } // truthQA

            // MC histograms
            if (xiCand.motherPDG() > 0) {
              histos.fill(HIST("Xi1530Rec"), lResonancePt, cent);
              histos.fill(HIST("Xi1530Recinvmass"), lResonanceMass);
              histos.fill(HIST("h3RecXi1530invmass"), cent, lResonancePt, lResonanceMass, lResonancePtMC);
            } else {
              histos.fill(HIST("Xi1530RecAnti"), lResonancePt, cent);
              histos.fill(HIST("Xi1530Recinvmass"), lResonanceMass);
              histos.fill(HIST("h3RecXi1530invmassAnti"), cent, lResonancePt, lResonanceMass, lResonancePtMC);
            }
          } // is MC
        } else { // Bkg candidates
          if constexpr (!IsMix) {
            if (pionCand.sign() < 0) {
              if (histoConfig.invMass1D) {
                histos.fill(HIST("Xi1530invmassLS"), lResonanceMass);
              }
              histos.fill(HIST("h3Xi1530invmassLS"), cent, lResonancePt, lResonanceMass, kLS);
            } else if (pionCand.sign() > 0) {
              if (histoConfig.invMass1D) {
                histos.fill(HIST("Xi1530invmassLSAnti"), lResonanceMass);
              }
              histos.fill(HIST("h3Xi1530invmassLSAnti"), cent, lResonancePt, lResonanceMass, kLS);
            }
          }
        } // -> End if signal or bkg
      } // -> End loop over xi candidates
    } // -> End loop over pion and xi candidates
  } // -> End fillHistograms

  void processData(ResoCollisions::iterator const& resoCollision,
                   ResoTracks const& resoTracks,
                   aod::ResoCascades const& cascTracks)
  {
    auto inCent = resoCollision.cent();
    // Version 001 stores the estimator chosen by cfgMultiplicityEstimator in the producer.
    const auto multiplicity = resoCollision.multiplicity();

    if (additionalConfig.cRecoINELgt0 && !resoCollision.isRecINELgt0()) { // Check reco INELgt0 (at least one PV track in |eta| < 1) about the collision
      return;
    }
    if (histoConfig.multQA) {
      histos.fill(HIST("multQA/h2MultCent"), inCent, multiplicity);
    }
    fillHistograms<false, false, false>(resoCollision, inCent, resoTracks, cascTracks);
  }

  // Calculate numerator for the Acceptance x Efficiency
  void processMC(ResoMCCols::iterator const& resoCollision,
                 soa::Join<aod::ResoCascades, aod::ResoMCCascades> const& cascTracks,
                 soa::Join<aod::ResoTracks, aod::ResoTrackTracks, aod::ResoMCTracks> const& resoTracks)
  {
    auto inCent = resoCollision.cent();
    const auto multiplicity = resoCollision.multiplicity();
    if ((additionalConfig.cRecoINELgt0 && !resoCollision.isRecINELgt0()) ||
        (additionalConfig.cMCINELgt0 && !resoCollision.isINELgt0()) ||
        (additionalConfig.cMCVtxIn10 && !resoCollision.isVtxIn10())) {
      return;
    }
    // The modular producer writes only selected reconstructed collisions.
    if (histoConfig.multQA) {
      histos.fill(HIST("multQA/h2MultCent"), inCent, multiplicity);
    }
    fillHistograms<false, true, false>(resoCollision, inCent, resoTracks, cascTracks);
  }

  void processMCMicro(ResoMCCols::iterator const& resoCollision,
                      soa::Join<aod::ResoCascades, aod::ResoMCCascades> const& cascTracks,
                      soa::Join<ResoMicroTracks, aod::ResoMCMicroTracks_001> const& resoTracks)
  {
    if ((additionalConfig.cRecoINELgt0 && !resoCollision.isRecINELgt0()) ||
        (additionalConfig.cMCINELgt0 && !resoCollision.isINELgt0()) ||
        (additionalConfig.cMCVtxIn10 && !resoCollision.isVtxIn10())) {
      return;
    }
    if (histoConfig.multQA) {
      histos.fill(HIST("multQA/h2MultCent"), resoCollision.cent(), resoCollision.multiplicity());
    }
    fillHistograms<true, true, false>(resoCollision, resoCollision.cent(), resoTracks, cascTracks);
  }

  // Calculate denominator for the Acceptance x Efficiency, actually it is not Trueth info...
  void processMCTrue(ResoMCCols::iterator const& resoCollision,
                     aod::ResoMCParents_001 const& resoParents)
  {
    auto multiplicity = resoCollision.mcMultiplicity();
    auto inCent = resoCollision.cent();
    if ((additionalConfig.cRecoINELgt0 && !resoCollision.isRecINELgt0()) ||
        (additionalConfig.cMCINELgt0 && !resoCollision.isINELgt0()) ||
        (additionalConfig.cMCVtxIn10 && !resoCollision.isVtxIn10())) {
      return;
    }
    if (histoConfig.multQA) {
      histos.fill(HIST("multQA/h2MultCentMC"), inCent, multiplicity);
    }

    for (const auto& part : resoParents) { // loop over all pre-filtered MC particles
      if (std::abs(part.pdgCode()) != kXiStar) {
        continue;
      }
      if ((part.y() <= primarytrackConfig.cfgRapidityMinCut) || (part.y() >= primarytrackConfig.cfgRapidityMaxCut)) {
        continue;
      }
      bool pass1 = std::abs(part.daughterPDG1()) == kPiPlus || std::abs(part.daughterPDG2()) == kPiPlus;
      bool pass2 = std::abs(part.daughterPDG1()) == kXiMinus || std::abs(part.daughterPDG2()) == kXiMinus;

      if (!pass1 || !pass2) {
        continue;
      }

      if (part.pdgCode() > 0) {
        histos.fill(HIST("h3Xi1530Gen"), part.pt(), inCent, multiplicity);
      } else {
        histos.fill(HIST("h3Xi1530GenAnti"), part.pt(), inCent, multiplicity);
      }
    }
  }

  void processDataMicro(ResoCollisions::iterator const& resoCollision,
                        ResoMicroTracks const& resomicrotracks,
                        aod::ResoCascades const& cascTracks)
  {
    const auto multiplicity = resoCollision.multiplicity();
    auto inCent = resoCollision.cent();
    if (additionalConfig.cRecoINELgt0 && !resoCollision.isRecINELgt0()) { // Check reco INELgt0 (at least one PV track in |eta| < 1) about the collision
      return;
    }
    if (histoConfig.multQA) {
      histos.fill(HIST("multQA/h2MultCent"), inCent, multiplicity);
    }
    fillHistograms<true, false, false>(resoCollision, inCent, resomicrotracks, cascTracks);
  }

  using BinningTypeVtxZT0M = ColumnBinningPolicy<aod::collision::PosZ, aod::resocollision::Cent>;
  void processMEMicro(ResoCollisions const& resoCollisions,
                      ResoMicroTracks const& resomicrotracks,
                      aod::ResoCascades const& cascTracks)
  {
    if (mixingConfig.nEvtMixing <= 0) {
      return;
    }
    auto tracksTuple = std::make_tuple(resomicrotracks, cascTracks);

    BinningTypeVtxZT0M colBinning{{mixingConfig.cfgVtxBins, mixingConfig.cfgMultBins}, true};
    Pair<ResoCollisions, ResoMicroTracks, aod::ResoCascades, BinningTypeVtxZT0M> pairs{colBinning, mixingConfig.nEvtMixing, -1, resoCollisions, tracksTuple, &cache};

    for (const auto& [collision1, tracks1, collision2, tracks2] : pairs) {

      const auto multiplicity = collision1.multiplicity();
      auto inCent = collision1.cent();
      if (additionalConfig.cRecoINELgt0 && (!collision1.isRecINELgt0() || !collision2.isRecINELgt0())) {
        continue;
      }
      if (!std::isfinite(collision1.bMagField()) || collision1.bMagField() != collision2.bMagField()) {
        continue;
      }
      if (histoConfig.multQA) {
        histos.fill(HIST("multQA/h2MultCentME"), inCent, multiplicity);
      }
      const MixingCollision cascadeCollision{.x = collision2.posX(), .y = collision2.posY(), .z = collision2.posZ(), .centrality = collision2.cent(), .multiplicity = collision2.multiplicity(), .index = collision2.globalIndex(), .recINELgt0 = collision2.isRecINELgt0()};
      fillHistograms<true, false, true>(collision1, inCent, tracks1, tracks2, &cascadeCollision);
    }
  }
  static std::shared_ptr<arrow::Table> copyMixingTracks(ResoMicroTracks const& tracks)
  {
    const auto source = tracks.asArrowTableConstrained();
    if (source->num_rows() == 0) {
      return arrow::Table::MakeEmpty(source->schema()).ValueOrDie();
    }
    std::vector<std::shared_ptr<arrow::ChunkedArray>> columns;
    columns.reserve(source->num_columns());
    for (const auto& column : source->columns()) {
      // Concatenate copies only the slice's rows, releasing the source DF buffers.
      columns.push_back(std::make_shared<arrow::ChunkedArray>(arrow::Concatenate(column->chunks()).ValueOrDie()));
    }
    return arrow::Table::Make(source->schema(), columns);
  }

  // Cross-DF mixing on module v001 input; no legacy Reso*DF merger is needed.
  // History is local to this task instance. Use inputs from one run/condition set.
  void processMEDF(ResoCollisions::iterator const& collision,
                   ResoMicroTracks const& tracks,
                   aod::ResoCascades const& cascades)
  {
    if (mixingConfig.nEvtMixing <= 0 || !std::isfinite(collision.bMagField())) {
      return;
    }
    if (mixingBField != collision.bMagField()) {
      mixingPools.clear();
      mixingBField = collision.bMagField();
    }
    BinningTypeVtxZT0M binning{{mixingConfig.cfgVtxBins, mixingConfig.cfgMultBins}, true};
    const int bin = binning.getBin(std::make_tuple(collision.posZ(), collision.cent()));
    if (bin < 0) {
      return;
    }
    const MixingCollision current{.x = collision.posX(), .y = collision.posY(), .z = collision.posZ(), .centrality = collision.cent(), .multiplicity = collision.multiplicity(), .index = collision.globalIndex(), .recINELgt0 = collision.isRecINELgt0()};
    auto& pool = mixingPools[bin];
    for (const auto& previous : pool) {
      const auto& anchor = previous.collision;
      if (additionalConfig.cRecoINELgt0 && (!anchor.recINELgt0 || !current.recINELgt0)) {
        continue;
      }
      if (histoConfig.multQA) {
        histos.fill(HIST("multQA/h2MultCentME"), anchor.centrality, anchor.multiplicity);
      }
      // Preserve the old pion -> new Cascade direction and the Cascade's own PV.
      fillHistograms<true, false, true>(anchor, anchor.centrality, ResoMicroTracks{previous.tracks}, cascades, &current);
    }
    // Insert after mixing: DF-local row IDs may repeat, but no event mixes with itself.
    pool.push_back({current, copyMixingTracks(tracks)});
    while (pool.size() > static_cast<std::size_t>(mixingConfig.nEvtMixing.value)) {
      pool.pop_front();
    }
  }

  PROCESS_SWITCH(Xi1530Analysisqa, processData, "Process Event for Data", false);
  PROCESS_SWITCH(Xi1530Analysisqa, processMC, "Process Event for MC (Reconstructed)", false);
  PROCESS_SWITCH(Xi1530Analysisqa, processMCMicro, "Process Event for MC (Reconstructed MicroTrack 001)", false);
  PROCESS_SWITCH(Xi1530Analysisqa, processMCTrue, "Process Event for MC (Generated)", false);
  PROCESS_SWITCH(Xi1530Analysisqa, processDataMicro, "Process Event for Data (MicroTrack)", false);
  PROCESS_SWITCH(Xi1530Analysisqa, processMEMicro, "Process EventMixing (MicroTrack) ", false);
  PROCESS_SWITCH(Xi1530Analysisqa, processMEDF, "Process cross-DF EventMixing (MicroTrack 001)", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& context)
{
  return WorkflowSpec{adaptAnalysisTask<Xi1530Analysisqa>(context)};
}
