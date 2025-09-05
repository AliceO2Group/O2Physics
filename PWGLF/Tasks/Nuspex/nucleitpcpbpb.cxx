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

/// \file nucleitpcpbpb.cxx
/// \brief nuclei analysis
/// \note under work
///
/// \author Jaideep Tanwar <jaideep.tanwar@cern.ch>, Panjab University

#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/PID/TPCPIDResponse.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsTPC/BetheBlochAleph.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include <Math/Vector4D.h>
#include <TRandom3.h>

#include <limits>
#include <string>
#include <vector>
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, o2::aod::TracksDCA, aod::TrackSelectionExtension, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCDe, aod::pidTPCTr, aod::pidTPCHe, aod::pidTPCAl, aod::TOFSignal, aod::TOFEvTime>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;
//---------------------------------------------------------------------------------------------------------------------------------
namespace
{
static const int nParticles = 6;
static const std::vector<std::string> particleNames{"pion", "proton", "deuteron", "triton", "helion", "alpha"};
static const std::vector<int> particlePdgCodes{211, 2212, o2::constants::physics::kDeuteron, o2::constants::physics::kTriton, o2::constants::physics::kHelium3, o2::constants::physics::kAlpha};
static const std::vector<double> particleMasses{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron, o2::constants::physics::MassTriton, o2::constants::physics::MassHelium3, o2::constants::physics::MassAlpha};
static const std::vector<int> particleCharge{1, 1, 1, 1, 2, 2};
const int nBetheParams = 6;
std::vector<int> hfMothCodes = {511, 521, 531, 541, 5122}; // b-mesons + Lambda_b
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
constexpr double kBetheBlochDefault[nParticles][nBetheParams]{
  {13.611469, 3.598765, -0.021138, 2.039562, 0.651040, 0.09},    // pion
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},      // proton
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},      // deuteron
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},      // triton
  {-126.557359, -0.858569, 1.111643, 1.210323, 2.656374, 0.09},  // helion
  {-126.557359, -0.858569, 1.111643, 1.210323, 2.656374, 0.09}}; // alpha
const int nTrkSettings = 13;
static const std::vector<std::string> trackPIDsettingsNames{"useBBparams", "minITSnCls", "minITSnClscos", "minTPCnCls", "maxTPCchi2", "minTPCchi2", "maxITSchi2", "maxTPCnSigma", "maxDcaXY", "maxDcaZ", "minITSclsSize", "minTPCnClsCrossedRows", "minReqClusterITSib"};
constexpr double kTrackPIDSettings[nParticles][nTrkSettings]{
  {0, 0, 4, 60, 4.0, 0.5, 100, 2.5, 2., 2., 0., 70, 1},
  {1, 0, 4, 70, 4.0, 0.5, 100, 3.0, 2., 2., 0., 70, 1},
  {1, 0, 4, 70, 4.0, 0.5, 100, 3.0, 2., 2., 0., 70, 1},
  {1, 0, 4, 70, 4.0, 0.5, 100, 3.0, 2., 2., 0., 70, 1},
  {1, 0, 4, 75, 4.0, 0.5, 100, 5.0, 2., 2., 0., 70, 1},
  {1, 0, 4, 70, 4.0, 0.5, 100, 5.0, 2., 2., 0., 70, 1}};

const int nTrkSettings2 = 6;
static const std::vector<std::string> trackPIDsettingsNames2{"useITSnsigma", "minITSnsigma", "maxITSnsigma", "fillsparsh", "useTPCnsigmaTOF", "maxTPCnsigmaTOF"};
constexpr double kTrackPIDSettings2[nParticles][nTrkSettings2]{
  {1, -5, 4, 0, 1, 2},
  {1, -5, 4, 0, 1, 2},
  {1, -5, 4, 0, 1, 2},
  {1, -5, 4, 1, 1, 2},
  {1, -5, 4, 1, 1, 2},
  {1, -5, 4, 1, 1, 2}};

struct PrimParticles {
  TString name;
  int pdgCode, charge;
  double mass, resolution;
  std::vector<double> betheParams;
  bool active;
  PrimParticles(std::string name_, int pdgCode_, double mass_, int charge_, LabeledArray<double> bethe) : name(name_), pdgCode(pdgCode_), charge(charge_), mass(mass_), active(false)
  {
    resolution = bethe.get(name, "resolution");
    betheParams.clear();
    constexpr unsigned int kNSpecies = 5;
    for (unsigned int i = 0; i < kNSpecies; i++)
      betheParams.push_back(bethe.get(name, i));
  }
}; // struct PrimParticles
//----------------------------------------------------------------------------------------------------------------
std::vector<std::shared_ptr<TH3>> hmass;
std::vector<std::shared_ptr<TH2>> hmassnsigma;
} // namespace
//----------------------------------------------------------------------------------------------------------------
struct NucleitpcPbPb {
  Preslice<aod::TrackAssoc> perCollision = aod::track_association::collisionId;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histomc{"histomc", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  Configurable<int> cfgDebug{"cfgDebug", 1, "debug level"};
  // event Selections cuts
  Configurable<bool> removeITSROFrameBorder{"removeITSROFrameBorder", false, "Remove TF border"};
  Configurable<bool> removeNoSameBunchPileup{"removeNoSameBunchPileup", false, "Remove TF border"};
  Configurable<bool> requireIsGoodZvtxFT0vsPV{"requireIsGoodZvtxFT0vsPV", false, "Remove TF border"};
  Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "Remove TF border"};
  Configurable<bool> removeNoTimeFrameBorder{"removeNoTimeFrameBorder", false, "Remove TF border"};
  Configurable<bool> cfgRigidityCorrection{"cfgRigidityCorrection", false, "apply rigidity correction"};
  // Track Selection Cuts
  Configurable<float> cfgCutEta{"cfgCutEta", 0.9f, "Eta range for tracks"};
  Configurable<bool> cfgetaRequire{"cfgetaRequire", true, "eta cut require"};
  Configurable<bool> cfgetaRequireMC{"cfgetaRequireMC", true, "eta cut require for generated particles"};
  Configurable<bool> cfgrapidityRequireMC{"cfgrapidityRequireMC", true, "rapidity cut require for generated particles"};
  Configurable<bool> cfgUsePVcontributors{"cfgUsePVcontributors", true, "use tracks that are PV contibutors"};
  Configurable<bool> cfgITSrequire{"cfgITSrequire", true, "Additional cut on ITS require"};
  Configurable<bool> cfgTPCrequire{"cfgTPCrequire", true, "Additional cut on TPC require"};
  Configurable<bool> cfgPassedITSRefit{"cfgPassedITSRefit", true, "Require ITS refit"};
  Configurable<bool> cfgPassedTPCRefit{"cfgPassedTPCRefit", true, "Require TPC refit"};
  Configurable<bool> cfgRapidityRequire{"cfgRapidityRequire", true, "Require Rapidity cut"};
  Configurable<bool> cfgTPCNClsfoundRequire{"cfgTPCNClsfoundRequire", true, "Require TPCNClsfound Cut"};
  Configurable<bool> cfgTPCNClsCrossedRowsRequire{"cfgTPCNClsCrossedRowsRequire", true, "Require TPCNClsCrossedRows Cut"};
  Configurable<bool> cfgmaxTPCchi2Require{"cfgmaxTPCchi2Require", true, "Require maxTPCchi2 Cut"};
  Configurable<bool> cfgminTPCchi2Require{"cfgminTPCchi2Require", true, "Require minTPCchi2 Cut"};
  Configurable<bool> cfgminITSnClsRequire{"cfgminITSnClsRequire", false, "Require minITSnCls Cut"};
  Configurable<bool> cfgmccorrectionhe4Require{"cfgmccorrectionhe4Require", true, "MC correction for pp he4 particle"};
  Configurable<bool> cfgminITSnClscosRequire{"cfgminITSnClscosRequire", true, "Require minITSnCls / cosh(eta) Cut"};
  Configurable<bool> cfgminReqClusterITSibRequire{"cfgminReqClusterITSibRequire", true, " Require min number of clusters required in ITS inner barrel"};
  Configurable<bool> cfgmaxITSchi2Require{"cfgmaxITSchi2Require", true, "Require maxITSchi2 Cut"};
  Configurable<bool> cfgmaxTPCnSigmaRequire{"cfgmaxTPCnSigmaRequire", true, "Require maxTPCnSigma Cut"};
  Configurable<bool> cfgminGetMeanItsClsSizeRequire{"cfgminGetMeanItsClsSizeRequire", true, "Require minGetMeanItsClsSize Cut"};
  Configurable<bool> cfgmaxGetMeanItsClsSizeRequire{"cfgmaxGetMeanItsClsSizeRequire", true, "Require maxGetMeanItsClsSize Cut"};
  Configurable<bool> cfgDCAwithptRequire{"cfgDCAwithptRequire", true, "Require DCA cuts with pt dependance"};
  Configurable<bool> cfgRequirebetaplot{"cfgRequirebetaplot", true, "Require beta plot"};
  Configurable<bool> cfgIncludeMaterialInEfficiency{"cfgIncludeMaterialInEfficiency", true, "Require from material in efficiency"};

  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {kBetheBlochDefault[0], nParticles, nBetheParams, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgTrackPIDsettings{"cfgTrackPIDsettings", {kTrackPIDSettings[0], nParticles, nTrkSettings, particleNames, trackPIDsettingsNames}, "track selection and PID criteria"};
  Configurable<LabeledArray<double>> cfgTrackPIDsettings2{"cfgTrackPIDsettings2", {kTrackPIDSettings2[0], nParticles, nTrkSettings2, particleNames, trackPIDsettingsNames2}, "track selection and PID criteria"};
  Configurable<bool> cfgFillhspectra{"cfgFillhspectra", true, "fill data sparsh"};
  Configurable<bool> cfgFillmass{"cfgFillmass", false, "Fill mass histograms"};
  Configurable<bool> cfgFillmassnsigma{"cfgFillmassnsigma", true, "Fill mass vs nsigma histograms"};
  Configurable<float> centcut{"centcut", 80.0f, "centrality cut"};
  Configurable<float> cfgCutRapidity{"cfgCutRapidity", 0.5f, "Rapidity range"};
  Configurable<float> cfgtpcNClsFindable{"cfgtpcNClsFindable", 0.8f, "tpcNClsFindable over crossedRows"}; /////////////
  Configurable<float> cfgZvertex{"cfgZvertex", 10, "Min Z Vertex"};
  Configurable<bool> cfgZvertexRequire{"cfgZvertexRequire", true, "Pos Z cut require"};
  Configurable<bool> cfgZvertexRequireMC{"cfgZvertexRequireMC", true, "Pos Z cut require for generated particles"};
  Configurable<bool> cfgsel8Require{"cfgsel8Require", true, "sel8 cut require"};
  o2::track::TrackParametrizationWithError<float> mTrackParCov;
  // Binning configuration
  ConfigurableAxis axisMagField{"axisMagField", {10, -10., 10.}, "magnetic field"};
  ConfigurableAxis axisNev{"axisNev", {5, 0., 5.}, "Number of events"};
  ConfigurableAxis axisRigidity{"axisRigidity", {4000, -10., 10.}, "#it{p}^{TPC}/#it{z}"};
  ConfigurableAxis axisdEdx{"axisdEdx", {4000, 0, 4000}, "d#it{E}/d#it{x}"};
  ConfigurableAxis axisCent{"axisCent", {100, 0, 100}, "centrality"};
  ConfigurableAxis axisVtxZ{"axisVtxZ", {120, -20, 20}, "z"};

  ConfigurableAxis ptAxis{"ptAxis", {200, 0, 10}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axiseta{"axiseta", {100, -1, 1}, "eta"};
  ConfigurableAxis axisrapidity{"axisrapidity", {100, -2, 2}, "rapidity"};
  ConfigurableAxis axismass{"axismass", {100, -10, 10}, "mass"};
  ConfigurableAxis axismassnsigma{"axismassnsigma", {100, 0, 20}, "mass"};
  ConfigurableAxis nsigmaAxis{"nsigmaAxis", {160, -10, 10}, "n#sigma_{#pi^{+}}"};
  ConfigurableAxis speciesBitAxis{"speciesBitAxis", {8, -0.5, 7.5}, "particle type 0: pion, 1: proton, 2: deuteron, 3: triton, 4:He3, 5:He4"};
  ConfigurableAxis speciesTrackingAxis{"speciesTrackingAxis", {11, -0.5, 10.5}, "particle type 0: pion, 1: proton, 2: deuteron, 3: triton, 4:He3, 5:He4"};
  ConfigurableAxis axisDCA{"axisDCA", {400, -10., 10.}, "DCA axis"};
  ConfigurableAxis particleAntiAxis{"particleAntiAxis", {2, 0, 2}, "Particle/Anti-particle"}; // 0 = particle, 1 = anti-particle
  ConfigurableAxis decayTypeAxis{"decayTypeAxis", {3, 0, 3}, "Decay type"};                   // 0 = primary, 1 = from decay

  // CCDB
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<double> bField{"bField", -999, "bz field, -999 is automatic"};
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> pidPath{"pidPath", "", "Path to the PID response object"};

  std::vector<PrimParticles> primaryParticles;
  std::vector<float> primVtx, cents;
  bool collHasCandidate, collPassedEvSel, collPassedEvSelMc;
  int mRunNumber, occupancy;
  float dBz, momn;
  TRandom3 rand;
  float he3 = 4;
  float he4 = 5;
  float processmaterial = 23;
  //----------------------------------------------------------------------------------------------------------------
  void init(InitContext const&)
  {
    mRunNumber = 0;
    dBz = 0;
    rand.SetSeed(0);
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    for (int i = 0; i < nParticles; i++) { // create primaryParticles
      primaryParticles.push_back(PrimParticles(particleNames.at(i), particlePdgCodes.at(i), particleMasses.at(i), particleCharge.at(i), cfgBetheBlochParams));
    }
    // create histograms
    if (doprocessData) {
      histos.add("histMagField", "histMagField", kTH1F, {axisMagField});
      histos.add("histNev", "histNev", kTH1F, {axisNev});
      histos.add("histVtxZ", "histVtxZ", kTH1F, {axisVtxZ});
      histos.add("histCentFT0C", "histCentFT0C", kTH1F, {axisCent});
      histos.add("histCentFT0M", "histCentFT0M", kTH1F, {axisCent});
      histos.add("histCentFTOC_cut", "histCentFTOC_cut", kTH1F, {axisCent});
      histos.add<THnSparse>("hSpectra", " ", HistType::kTHnSparseF, {speciesBitAxis, ptAxis, nsigmaAxis, {5, -2.5, 2.5}, axisCent, axisDCA, axisDCA});
    }
    histos.add("histeta", "histeta", kTH1F, {axiseta});
    histos.add("Tofsignal", "Tofsignal", kTH2F, {axisRigidity, {4000, 0.2, 1.2, "#beta"}});
    histos.add("Tpcsignal", "Tpcsignal", kTH2F, {axisRigidity, axisdEdx});

    hmass.resize(2 * nParticles + 2);
    hmassnsigma.resize(2 * nParticles + 2);

    for (int i = 0; i < nParticles; i++) {
      TString histName = primaryParticles[i].name;
      if (cfgFillmass) {
        hmass[2 * i] = histos.add<TH3>(Form("histmass_pt/histmass_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); mass^{2};  centrality(%)", HistType::kTH3F, {ptAxis, axismass, axisCent});
        hmass[2 * i + 1] = histos.add<TH3>(Form("histmass_ptanti/histmass_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); mass^{2}; centrality(%)", HistType::kTH3F, {ptAxis, axismass, axisCent});
      }
    }
    for (int i = 0; i < nParticles; i++) {
      TString histName = primaryParticles[i].name;
      if (cfgFillmassnsigma) {
        hmassnsigma[2 * i] = histos.add<TH2>(Form("histmass_nsigma/histmass_%s", histName.Data()), ";nsigma; mass^{2}", HistType::kTH2F, {nsigmaAxis, axismassnsigma});
        hmassnsigma[2 * i + 1] = histos.add<TH2>(Form("histmass_nsigmaanti/histmass_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); mass^{2}", HistType::kTH2F, {nsigmaAxis, axismassnsigma});
      }
    }

    if (doprocessMC) {
      histomc.add<THnSparse>("hSpectramc", " ", HistType::kTHnSparseF, {speciesBitAxis, {5, -2.5, 2.5}, axisCent, ptAxis, ptAxis});

      histomc.add("hDenomEffAcc", "Denominator for Efficiency x Acceptance",
                  {HistType::kTHnSparseF, {speciesBitAxis, ptAxis, axisrapidity, axisCent, particleAntiAxis, decayTypeAxis}});
      histomc.add("hNumerEffAcc", "Numerator for Efficiency x Acceptance",
                  {HistType::kTHnSparseF, {speciesBitAxis, ptAxis, axisrapidity, axisCent, particleAntiAxis, decayTypeAxis}});
      histomc.add("hDenomSignalLoss", "Denominator for Signal Loss",
                  {HistType::kTHnSparseF, {speciesBitAxis, ptAxis, axisrapidity, axisCent, particleAntiAxis, decayTypeAxis}});
      histomc.add("hNumerSignalLoss", "Numerator for Signal Loss",
                  {HistType::kTHnSparseF, {speciesBitAxis, ptAxis, axisrapidity, axisCent, particleAntiAxis, decayTypeAxis}});

      // Material secondary histograms
      histomc.add("histSecondaryMaterialHe3", "He3 from material", HistType::kTH2F, {ptAxis, axisCent});
      histomc.add("histSecondaryMaterialAntiHe3", "Anti-He3 from material", HistType::kTH2F, {ptAxis, axisCent});
      histomc.add("histSecondaryMaterialHe4", "He4 from material", HistType::kTH2F, {ptAxis, axisCent});
      histomc.add("histSecondaryMaterialAntiHe4", "Anti-He4 from material", HistType::kTH2F, {ptAxis, axisCent});

      histomc.add("histVtxZgen", "histVtxZgen", kTH1F, {axisVtxZ});
      histomc.add("histEtagen", "histEtagen", kTH1F, {axiseta});
      histomc.add("histPtgenHe3", "histPtgenHe3", kTH1F, {ptAxis});
      histomc.add("histPtgenAntiHe3", "histPtgenAntiHe3", kTH1F, {ptAxis});
      histomc.add("histPtgenHe4", "histPtgenHe4", kTH1F, {ptAxis});
      histomc.add("histPtgenAntiHe4", "histPtgenAntiHe4", kTH1F, {ptAxis});
      histomc.add("histNevReco", "histNevReco", kTH1F, {axisNev});
      histomc.add("histVtxZReco", "histVtxZReco", kTH1F, {axisVtxZ});
      histomc.add("histCentFT0CReco", "histCentFT0CReco", kTH1F, {axisCent});
      histomc.add("histCentFT0MReco", "histCentFT0MReco", kTH1F, {axisCent});
      histomc.add("histDeltaPtVsPtGen", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10}, {1000, -0.5, 0.5, "p_{T}(reco) - p_{T}(gen);p_{T}(reco)"}});
      histomc.add("histDeltaPtVsPtGenanti", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10}, {1000, -0.5, 0.5, "p_{T}(reco) - p_{T}(gen);p_{T}(reco)"}});
      histomc.add("histDeltaPtVsPtGenHe4", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10}, {1000, -0.5, 0.5, "p_{T}(reco) - p_{T}(gen);p_{T}(reco)"}});
      histomc.add("histDeltaPtVsPtGenHe4anti", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10}, {1000, -0.5, 0.5, "p_{T}(reco) - p_{T}(gen);p_{T}(reco)"}});
      histomc.add("histPIDtrack", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10, "p_{T}(reco)"}, {9, -0.5, 8.5, "p_{T}(reco) - p_{T}(gen)"}});
      histomc.add("histPIDtrackanti", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10, "p_{T}(reco)"}, {9, -0.5, 8.5, "p_{T}(reco) - p_{T}(gen)"}});
      histomc.add("histPIDtrackhe4", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10, "p_{T}(reco)"}, {9, -0.5, 8.5, "p_{T}(reco) - p_{T}(gen)"}});
      histomc.add("histPIDtrackantihe4", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10, "p_{T}(reco)"}, {9, -0.5, 8.5, "p_{T}(reco) - p_{T}(gen)"}});

      histomc.add("histWeakDecayPtHe3", "Pt distribution of He3 from weak decays", kTH2F, {ptAxis, axisCent});
      histomc.add("histWeakDecayPtAntiHe3", "Pt distribution of Anti-He3 from weak decays", kTH2F, {ptAxis, axisCent});
      histomc.add("histWeakDecayPtHe4", "Pt distribution of He4 from weak decays", kTH2F, {ptAxis, axisCent});
      histomc.add("histWeakDecayPtAntiHe4", "Pt distribution of Anti-He4 from weak decays", kTH2F, {ptAxis, axisCent});
      histomc.add("histRecoWeakDecayPtHe3", "Pt distribution of reconstructed He3 from weak decays", kTH2F, {ptAxis, axisCent});
      histomc.add("histRecoWeakDecayPtHe4", "Pt distribution of reconstructed He4 from weak decays", kTH2F, {ptAxis, axisCent});
      histomc.add("histRecoSecondaryMaterialHe3", "Reco He3 from material", HistType::kTH2F, {ptAxis, axisCent});
      histomc.add("histRecoSecondaryMaterialAntiHe3", "Reco Anti-He3 from material", HistType::kTH2F, {ptAxis, axisCent});
      histomc.add("histRecoSecondaryMaterialHe4", "Reco He4 from material", HistType::kTH2F, {ptAxis, axisCent});
      histomc.add("histRecoSecondaryMaterialAntiHe4", "Reco Anti-He4 from material", HistType::kTH2F, {ptAxis, axisCent});
      histomc.add("histProcessCodeHe", "Process codes for He isotopes", kTH1F, {{100, 0, 100, "process code"}});
      histomc.add("histProcess23Details", "Process 23 details", kTH2F, {{4, 0.5, 4.5, "particle type"}, {100, 0, 10, "p_{T}"}});
      histomc.add("histAllMaterialSecondariesGen", "All material secondaries (gen)", kTH3F, {{100, 0, 10, "p_{T}"}, {20, -1, 1, "y"}, {5, -0.5, 4.5, "type"}});
      histomc.add("histAllMaterialSecondariesReco", "All material secondaries (reco)", kTH3F, {{100, 0, 10, "p_{T}"}, {20, -1, 1, "y"}, {5, -0.5, 4.5, "type"}});

      // Add axis labels for type: 0=unknown, 1=He3, 2=anti-He3, 3=He4, 4=anti-He4
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------------------
  void processData(CollisionsFull const& collisions, TracksFull const& tracks, aod::BCsWithTimestamps const&, aod::TrackAssoc const& tracksColl)
  {
    for (const auto& collision : collisions) {
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      initCollision(collision);
      if (!collPassedEvSel)
        continue;
      if (collision.centFT0C() > centcut)
        continue;
      if (removeITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
        continue;
      if (removeNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup))
        continue;
      if (requireIsGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
        continue;
      if (requireIsVertexITSTPC && !collision.selection_bit(aod::evsel::kIsVertexITSTPC))
        continue;
      if (removeNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
        continue;
      histos.fill(HIST("histCentFTOC_cut"), collision.centFT0C());
      histos.fill(HIST("histNev"), 2.5);
      const uint64_t collIdx = collision.globalIndex();
      auto tracksByColl = tracksColl.sliceBy(perCollision, collIdx);
      /////////////////////////////////////////////////////////////////////////////////
      for (const auto& trackId : tracksByColl) {
        const auto& track = tracks.rawIteratorAt(trackId.trackId());
        if (!track.isPVContributor() && cfgUsePVcontributors)
          continue;
        if (!track.hasITS() && cfgITSrequire)
          continue;
        if (!track.hasTPC() && cfgTPCrequire)
          continue;
        if (!track.passedITSRefit() && cfgPassedITSRefit)
          continue;
        if (!track.passedTPCRefit() && cfgPassedTPCRefit)
          continue;
        if (std::abs(track.eta()) > cfgCutEta && cfgetaRequire)
          continue;
        for (size_t i = 0; i < primaryParticles.size(); i++) {

          float ptMomn;
          setTrackParCov(track, mTrackParCov);
          mTrackParCov.setPID(track.pidForTracking());
          ptMomn = (i == he3 || i == he4) ? 2 * mTrackParCov.getPt() : mTrackParCov.getPt();
          int sign = 0;
          if (track.sign() > 0) {
            sign = 1;
          }
          if (track.sign() < 0) {
            sign = -1;
          }
          if (std::abs(getRapidity(track, i)) > cfgCutRapidity && cfgRapidityRequire)
            continue;
          if (track.tpcNClsFound() < cfgTrackPIDsettings->get(i, "minTPCnCls") && cfgTPCNClsfoundRequire)
            continue;
          if (((track.tpcNClsCrossedRows() < cfgTrackPIDsettings->get(i, "minTPCnClsCrossedRows")) || track.tpcNClsCrossedRows() < cfgtpcNClsFindable * track.tpcNClsFindable()) && cfgTPCNClsCrossedRowsRequire)
            continue;
          if (track.tpcChi2NCl() > cfgTrackPIDsettings->get(i, "maxTPCchi2") && cfgmaxTPCchi2Require)
            continue;
          if (track.tpcChi2NCl() < cfgTrackPIDsettings->get(i, "minTPCchi2") && cfgminTPCchi2Require)
            continue;
          if (track.itsNCls() < cfgTrackPIDsettings->get(i, "minITSnCls") && cfgminITSnClsRequire)
            continue;
          double cosheta = std::cosh(track.eta());
          if ((track.itsNCls() / cosheta) < cfgTrackPIDsettings->get(i, "minITSnClscos") && cfgminITSnClscosRequire)
            continue;
          if ((track.itsNClsInnerBarrel() < cfgTrackPIDsettings->get(i, "minReqClusterITSib")) && cfgminReqClusterITSibRequire)
            continue;
          if (track.itsChi2NCl() > cfgTrackPIDsettings->get(i, "maxITSchi2") && cfgmaxITSchi2Require)
            continue;
          if (getMeanItsClsSize(track) < cfgTrackPIDsettings->get(i, "minITSclsSize") && cfgminGetMeanItsClsSizeRequire)
            continue;

          bool insideDCAxy = (std::abs(track.dcaXY()) <= (cfgTrackPIDsettings->get(i, "maxDcaXY") * (0.0105f + 0.0350f / std::pow(ptMomn, 1.1f)))); // o2-linter: disable=magic-number (To be checked)
          if ((!(insideDCAxy) || std::abs(track.dcaZ()) > dcazSigma(ptMomn, cfgTrackPIDsettings->get(i, "maxDcaZ"))) && cfgDCAwithptRequire)
            continue;

          float tpcNsigma = getTPCnSigma(track, primaryParticles.at(i));
          if ((std::abs(tpcNsigma) > cfgTrackPIDsettings->get(i, "maxTPCnSigma")) && cfgmaxTPCnSigmaRequire)
            continue;
          float itsSigma = getITSnSigma(track, primaryParticles.at(i));
          if (itsSigma < cfgTrackPIDsettings2->get(i, "minITSnsigma") && cfgTrackPIDsettings2->get(i, "useITSnsigma") < 1)
            continue;
          if (itsSigma > cfgTrackPIDsettings2->get(i, "maxITSnsigma") && cfgTrackPIDsettings2->get(i, "useITSnsigma") < 1)
            continue;
          histos.fill(HIST("Tpcsignal"), getRigidity(track) * track.sign(), track.tpcSignal());

          if (cfgFillhspectra && cfgTrackPIDsettings2->get(i, "fillsparsh") == 1) {
            histos.fill(HIST("hSpectra"), i, ptMomn, tpcNsigma, sign, collision.centFT0C(), track.dcaZ(), track.dcaXY());
          }
          fillhmassnsigma(track, i, tpcNsigma);
          if ((std::abs(tpcNsigma) > cfgTrackPIDsettings2->get(i, "maxTPCnsigmaTOF")) && cfgTrackPIDsettings2->get(i, "useTPCnsigmaTOF") < 1)
            continue;
          fillhmass(track, i, collision.centFT0C());

          if (cfgRequirebetaplot) {
            histos.fill(HIST("Tofsignal"), getRigidity(track) * track.sign(), o2::pid::tof::Beta::GetBeta(track));
          }
        }
        histos.fill(HIST("histeta"), track.eta());
      } // track loop
      ///////////////////////////////////////////////
    }
  }
  PROCESS_SWITCH(NucleitpcPbPb, processData, "data analysis", false);

  //----------------------------------------------------------------------------------------------------------------
  // MC particles - Efficiency x Acceptance and Signal Loss calculations
  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  struct McCollInfo {
    bool passedEvSel = false;
    float centrality = -1.0f;
  };
  std::vector<McCollInfo> mcCollInfos;

  void processMC(CollisionsFullMC const& collisions,
                 aod::McCollisions const& mcCollisions,
                 TracksFull const& tracks,
                 aod::BCsWithTimestamps const& bcs,
                 aod::McParticles const& particlesMC,
                 aod::McTrackLabels const& trackLabelsMC,
                 aod::McCollisionLabels const& collLabels,
                 aod::TrackAssoc const& tracksColl,
                 CollisionsFull const& colls)
  {

    (void)bcs;
    (void)collLabels;
    (void)colls;

    mcCollInfos.clear();
    mcCollInfos.resize(mcCollisions.size());
    for (auto const& collision : collisions) {
      int mcCollIdx = collision.mcCollisionId();
      if (mcCollIdx < 0 || mcCollIdx >= static_cast<int>(mcCollisions.size())) {
        continue;
      }
      if (std::abs(collision.posZ()) > cfgZvertex)
        continue;
      if (!collision.sel8())
        continue;
      if (collision.centFT0C() > centcut)
        continue;

      mcCollInfos[mcCollIdx].passedEvSel = true;
      mcCollInfos[mcCollIdx].centrality = collision.centFT0C();
    }
    for (auto const& mcCollision : mcCollisions) {
      size_t idx = mcCollision.globalIndex();

      if (idx >= mcCollInfos.size())
        continue;

      for (auto const& mcParticle : particlesMC) {

        if (mcParticle.mcCollisionId() != mcCollision.globalIndex())
          continue;

        int pdgCode = mcParticle.pdgCode();
        bool isHe3 = (std::abs(pdgCode) == particlePdgCodes.at(4));
        bool isHe4 = (std::abs(pdgCode) == particlePdgCodes.at(5));

        if (std::abs(mcParticle.eta()) > cfgCutEta)
          continue;
        if (std::abs(mcParticle.y()) > cfgCutRapidity)
          continue;
        bool isMaterialSecondary = false;
        if (!mcParticle.isPhysicalPrimary() && (isHe3 || isHe4)) {
          if (!mcParticle.has_mothers()) {
            isMaterialSecondary = true;
          } else {
            auto mother = mcParticle.mothers_as<aod::McParticles>().front();
            float notmother = 1000000000;
            if (std::abs(mother.pdgCode()) < notmother) {
              isMaterialSecondary = true;
            }
          }
        }
        if (isHe3 || isHe4) {
          histomc.fill(HIST("histProcessCodeHe"), mcParticle.getProcess());
          if (mcParticle.getProcess() == processmaterial) {
            histomc.fill(HIST("histProcess23Details"), std::abs(pdgCode), mcParticle.pt());
          }
        }

        if (isMaterialSecondary) {
          float centrality = mcCollInfos[idx].passedEvSel ? mcCollInfos[idx].centrality : -1.0f;

          if (pdgCode == particlePdgCodes.at(4)) {
            histomc.fill(HIST("histSecondaryMaterialHe3"), mcParticle.pt(), centrality);
          } else if (pdgCode == -particlePdgCodes.at(4)) {
            histomc.fill(HIST("histSecondaryMaterialAntiHe3"), mcParticle.pt(), centrality);
          } else if (pdgCode == particlePdgCodes.at(5)) {
            histomc.fill(HIST("histSecondaryMaterialHe4"), mcParticle.pt(), centrality);
          } else if (pdgCode == -particlePdgCodes.at(5)) {
            histomc.fill(HIST("histSecondaryMaterialAntiHe4"), mcParticle.pt(), centrality);
          }
          int type = 0;
          if (std::abs(pdgCode) == particlePdgCodes.at(4))
            type = (pdgCode > 0) ? 1 : 2;
          else if (std::abs(pdgCode) == particlePdgCodes.at(5))
            type = (pdgCode > 0) ? 3 : 4;
          histomc.fill(HIST("histAllMaterialSecondariesGen"), mcParticle.pt(), mcParticle.y(), type);
        }

        if (!isHe3 && !isHe4)
          continue;

        int decayType = 0;
        int particleAnti = (pdgCode > 0) ? 0 : 1;

        if (mcParticle.isPhysicalPrimary()) {
          decayType = 0;
          if (mcParticle.has_mothers()) {
            for (const auto& motherparticle : mcParticle.mothers_as<aod::McParticles>()) {
              if (std::find(hfMothCodes.begin(), hfMothCodes.end(),
                            std::abs(motherparticle.pdgCode())) != hfMothCodes.end()) {
                decayType = 1;
                break;
              }
            }
          }
        } else if (mcParticle.has_mothers()) {
          decayType = 1;
        } else {
          decayType = 2;
          continue;
        }
        bool isFromWeakDecay = (decayType == 1);

        if (!mcParticle.isPhysicalPrimary() && !isFromWeakDecay)
          continue;

        int particleType = -1;
        if (std::abs(pdgCode) == particlePdgCodes.at(4))
          particleType = he3;
        else if (std::abs(pdgCode) == particlePdgCodes.at(5))
          particleType = he4;

        if (particleType >= 0) {
          float centrality = mcCollInfos[idx].passedEvSel ? mcCollInfos[idx].centrality : -1.0f;
          histomc.fill(HIST("hDenomSignalLoss"), particleType, mcParticle.pt(), mcParticle.y(), centrality, particleAnti, decayType);

          if (mcCollInfos[idx].passedEvSel) {
            histomc.fill(HIST("hNumerSignalLoss"), particleType, mcParticle.pt(), mcParticle.y(), centrality, particleAnti, decayType);
          }
        }

        if (isFromWeakDecay) {
          float centrality = mcCollInfos[idx].passedEvSel ? mcCollInfos[idx].centrality : -1.0f;
          if (pdgCode == particlePdgCodes.at(4)) {
            histomc.fill(HIST("histWeakDecayPtHe3"), mcParticle.pt(), centrality);
          } else if (pdgCode == -particlePdgCodes.at(4)) {
            histomc.fill(HIST("histWeakDecayPtAntiHe3"), mcParticle.pt(), centrality);
          } else if (pdgCode == particlePdgCodes.at(5)) {
            histomc.fill(HIST("histWeakDecayPtHe4"), mcParticle.pt(), centrality);
          } else if (pdgCode == -particlePdgCodes.at(5)) {
            histomc.fill(HIST("histWeakDecayPtAntiHe4"), mcParticle.pt(), centrality);
          }
        }
        if (pdgCode == particlePdgCodes.at(4)) {
          histomc.fill(HIST("histPtgenHe3"), mcParticle.pt());
        } else if (pdgCode == -particlePdgCodes.at(4)) {
          histomc.fill(HIST("histPtgenAntiHe3"), mcParticle.pt());
        } else if (pdgCode == particlePdgCodes.at(5)) {
          histomc.fill(HIST("histPtgenHe4"), mcParticle.pt());
        } else if (pdgCode == -particlePdgCodes.at(5)) {
          histomc.fill(HIST("histPtgenAntiHe4"), mcParticle.pt());
        }
      }
    }
    for (auto const& mcCollision : mcCollisions) {
      size_t idx = mcCollision.globalIndex();

      if (idx >= mcCollInfos.size())
        continue;
      if (!mcCollInfos[idx].passedEvSel)
        continue;
      if (std::abs(mcCollision.posZ()) > cfgZvertex)
        continue;

      histomc.fill(HIST("histVtxZgen"), mcCollision.posZ());
      float centrality = mcCollInfos[idx].centrality;

      for (auto const& mcParticle : particlesMC) {
        if (mcParticle.mcCollisionId() != mcCollision.globalIndex())
          continue;

        int pdgCode = mcParticle.pdgCode();
        bool isHe3 = (std::abs(pdgCode) == particlePdgCodes.at(4));
        bool isHe4 = (std::abs(pdgCode) == particlePdgCodes.at(5));

        if (!isHe3 && !isHe4)
          continue;
        if (std::abs(mcParticle.eta()) > cfgCutEta)
          continue;
        if (std::abs(mcParticle.y()) > cfgCutRapidity)
          continue;

        int decayType = 0;
        int particleAnti = (pdgCode > 0) ? 0 : 1;

        if (mcParticle.isPhysicalPrimary()) {
          decayType = 0;
          if (mcParticle.has_mothers()) {
            for (const auto& motherparticle : mcParticle.mothers_as<aod::McParticles>()) {
              if (std::find(hfMothCodes.begin(), hfMothCodes.end(),
                            std::abs(motherparticle.pdgCode())) != hfMothCodes.end()) {
                decayType = 1;
                break;
              }
            }
          }
        } else if (mcParticle.has_mothers()) {
          decayType = 1;
        } else {
          decayType = 2;
          continue;
        }
        bool isFromWeakDecay = (decayType == 1);

        if (!mcParticle.isPhysicalPrimary() && !isFromWeakDecay)
          continue;

        int particleType = -1;
        if (std::abs(pdgCode) == particlePdgCodes.at(4))
          particleType = he3;
        else if (std::abs(pdgCode) == particlePdgCodes.at(5))
          particleType = he4;

        if (particleType >= 0) {
          histomc.fill(HIST("hDenomEffAcc"), particleType, mcParticle.pt(), mcParticle.y(), centrality, particleAnti, decayType);
        }
      }
    }
    for (auto const& collision : collisions) {
      auto mcCollIdx = collision.mcCollisionId();
      if (mcCollIdx < 0 || mcCollIdx >= static_cast<int>(mcCollisions.size()))
        continue;

      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      collHasCandidate = false;
      histomc.fill(HIST("histNevReco"), 0.5);

      if (std::abs(collision.posZ()) > cfgZvertex && cfgZvertexRequire)
        continue;
      collPassedEvSel = collision.sel8();
      if (!collPassedEvSel && cfgsel8Require)
        continue;
      histomc.fill(HIST("histNevReco"), 1.5);
      histomc.fill(HIST("histVtxZReco"), collision.posZ());
      histomc.fill(HIST("histCentFT0CReco"), collision.centFT0C());
      histomc.fill(HIST("histCentFT0MReco"), collision.centFT0M());
      if (collision.centFT0C() > centcut)
        continue;
      histomc.fill(HIST("histNevReco"), 2.5);
      if (removeITSROFrameBorder && !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))
        continue;
      if (removeNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup))
        continue;
      if (requireIsGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
        continue;
      if (requireIsVertexITSTPC && !collision.selection_bit(aod::evsel::kIsVertexITSTPC))
        continue;
      if (removeNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
        continue;
      for (const auto& assoc : tracksColl) {
        if (assoc.collisionId() != collision.globalIndex())
          continue;
        const auto& track = tracks.rawIteratorAt(assoc.trackId());
        auto labelRow = trackLabelsMC.iteratorAt(track.globalIndex());
        int label = labelRow.mcParticleId();
        if (label < 0 || label >= static_cast<int>(particlesMC.size()))
          continue;
        auto const& matchedMCParticle = particlesMC.iteratorAt(label);

        int pdg = matchedMCParticle.pdgCode();
        bool isHe3 = (std::abs(pdg) == particlePdgCodes.at(4));
        bool isHe4 = (std::abs(pdg) == particlePdgCodes.at(5));

        bool isMaterialSecondary = false;
        if ((isHe3 || isHe4) && !matchedMCParticle.isPhysicalPrimary()) {

          if (matchedMCParticle.getProcess() == processmaterial) { // kPMaterial
            isMaterialSecondary = true;
          } else if (!matchedMCParticle.has_mothers()) {
            isMaterialSecondary = true;
          }
        }

        if (isMaterialSecondary) {
          if (pdg == particlePdgCodes.at(4)) {
            histomc.fill(HIST("histRecoSecondaryMaterialHe3"), track.pt(), collision.centFT0C());
          } else if (pdg == -particlePdgCodes.at(4)) {
            histomc.fill(HIST("histRecoSecondaryMaterialAntiHe3"), track.pt(), collision.centFT0C());
          } else if (pdg == particlePdgCodes.at(5)) {
            histomc.fill(HIST("histRecoSecondaryMaterialHe4"), track.pt(), collision.centFT0C());
          } else if (pdg == -particlePdgCodes.at(5)) {
            histomc.fill(HIST("histRecoSecondaryMaterialAntiHe4"), track.pt(), collision.centFT0C());
          }

          int type = 0;
          if (std::abs(pdg) == particlePdgCodes.at(4))
            type = (pdg > 0) ? 1 : 2;
          else if (std::abs(pdg) == particlePdgCodes.at(5))
            type = (pdg > 0) ? 3 : 4;
          histomc.fill(HIST("histAllMaterialSecondariesReco"), track.pt(), getRapidity(track, he3), type);

          if (!cfgIncludeMaterialInEfficiency) {
            continue;
          }
        }
        int decayType = 0;
        bool isFromWeakDecay = false;

        if (matchedMCParticle.isPhysicalPrimary()) {
          decayType = 0;
          if (matchedMCParticle.has_mothers()) {
            for (const auto& motherparticle : matchedMCParticle.mothers_as<aod::McParticles>()) {
              if (std::find(hfMothCodes.begin(), hfMothCodes.end(),
                            std::abs(motherparticle.pdgCode())) != hfMothCodes.end()) {
                isFromWeakDecay = true;
                decayType = 1;
                break;
              }
            }
          }
        } else if (matchedMCParticle.has_mothers()) {
          isFromWeakDecay = true;
          decayType = 1;
        } else {
          decayType = 2;
        }

        if (!track.isPVContributor() && cfgUsePVcontributors)
          continue;
        if (!track.hasITS() && cfgITSrequire)
          continue;
        if (!track.hasTPC() && cfgTPCrequire)
          continue;
        if (!track.passedITSRefit() && cfgPassedITSRefit)
          continue;
        if (!track.passedTPCRefit() && cfgPassedTPCRefit)
          continue;
        if (std::abs(track.eta()) > cfgCutEta && cfgetaRequire)
          continue;
        if (!matchedMCParticle.isPhysicalPrimary() && !isFromWeakDecay && !isMaterialSecondary)
          continue;

        for (size_t i = 0; i < primaryParticles.size(); i++) {
          if (std::abs(pdg) != std::abs(particlePdgCodes.at(i)))
            continue;
          float ptMomn;
          setTrackParCov(track, mTrackParCov);
          mTrackParCov.setPID(track.pidForTracking());
          ptMomn = (i == he3 || i == he4) ? 2 * mTrackParCov.getPt() : mTrackParCov.getPt();
          int sign = 0;
          if (track.sign() > 0) {
            sign = 1;
          }
          if (track.sign() < 0) {
            sign = -1;
          }
          if (std::abs(getRapidity(track, i)) > cfgCutRapidity && cfgRapidityRequire)
            continue;
          if (track.tpcNClsFound() < cfgTrackPIDsettings->get(i, "minTPCnCls") && cfgTPCNClsfoundRequire)
            continue;
          if (((track.tpcNClsCrossedRows() < cfgTrackPIDsettings->get(i, "minTPCnClsCrossedRows")) || track.tpcNClsCrossedRows() < cfgtpcNClsFindable * track.tpcNClsFindable()) && cfgTPCNClsCrossedRowsRequire)
            continue;
          if (track.tpcChi2NCl() > cfgTrackPIDsettings->get(i, "maxTPCchi2") && cfgmaxTPCchi2Require)
            continue;
          if (track.tpcChi2NCl() < cfgTrackPIDsettings->get(i, "minTPCchi2") && cfgminTPCchi2Require)
            continue;
          if (track.itsNCls() < cfgTrackPIDsettings->get(i, "minITSnCls") && cfgminITSnClsRequire)
            continue;
          double cosheta = std::cosh(track.eta());
          if ((track.itsNCls() / cosheta) < cfgTrackPIDsettings->get(i, "minITSnClscos") && cfgminITSnClscosRequire)
            continue;
          if ((track.itsNClsInnerBarrel() < cfgTrackPIDsettings->get(i, "minReqClusterITSib")) && cfgminReqClusterITSibRequire)
            continue;
          if (track.itsChi2NCl() > cfgTrackPIDsettings->get(i, "maxITSchi2") && cfgmaxITSchi2Require)
            continue;
          if (getMeanItsClsSize(track) < cfgTrackPIDsettings->get(i, "minITSclsSize") && cfgminGetMeanItsClsSizeRequire)
            continue;

          bool insideDCAxy = (std::abs(track.dcaXY()) <= (cfgTrackPIDsettings->get(i, "maxDcaXY") * (0.0105f + 0.0350f / std::pow(ptMomn, 1.1f))));
          if ((!(insideDCAxy) || std::abs(track.dcaZ()) > dcazSigma(ptMomn, cfgTrackPIDsettings->get(i, "maxDcaZ"))) && cfgDCAwithptRequire)
            continue;

          float tpcNsigma = getTPCnSigma(track, primaryParticles.at(i));
          if ((std::abs(tpcNsigma) > cfgTrackPIDsettings->get(i, "maxTPCnSigma")) && cfgmaxTPCnSigmaRequire)
            continue;
          float itsSigma = getITSnSigma(track, primaryParticles.at(i));
          if (itsSigma < cfgTrackPIDsettings2->get(i, "minITSnsigma") && cfgTrackPIDsettings2->get(i, "useITSnsigma") < 1)
            continue;
          if (itsSigma > cfgTrackPIDsettings2->get(i, "maxITSnsigma") && cfgTrackPIDsettings2->get(i, "useITSnsigma") < 1)
            continue;
          float ptReco;
          mTrackParCov.setPID(track.pidForTracking());
          ptReco = (std::abs(pdg) == particlePdgCodes.at(4) || std::abs(pdg) == particlePdgCodes.at(5)) ? 2 * mTrackParCov.getPt() : mTrackParCov.getPt();

          int particleAnti = (pdg > 0) ? 0 : 1;

          if (i == he3 || i == he4) {
            histomc.fill(HIST("hNumerEffAcc"), i, matchedMCParticle.pt(), matchedMCParticle.y(), collision.centFT0C(), particleAnti, decayType);
          }
          if (isFromWeakDecay) {
            if (std::abs(pdg) == particlePdgCodes.at(4)) {
              histomc.fill(HIST("histRecoWeakDecayPtHe3"), ptReco, collision.centFT0C());
            } else if (std::abs(pdg) == particlePdgCodes.at(5)) {
              histomc.fill(HIST("histRecoWeakDecayPtHe4"), ptReco, collision.centFT0C());
            }
          }

          histos.fill(HIST("Tpcsignal"), getRigidity(track) * track.sign(), track.tpcSignal());

          if (pdg == -particlePdgCodes.at(5) && cfgmccorrectionhe4Require) {
            ptReco = ptReco + 0.00765 + 0.503791 * std::exp(-1.10517 * ptReco);
          }

          if (pdg == -particlePdgCodes.at(4) && cfgmccorrectionhe4Require) {
            int pidGuess = track.pidForTracking();
            int antitriton = 6;
            if (pidGuess == antitriton) {
              ptReco = ptReco - 0.464215 + 0.195771 * ptReco - 0.0183111 * ptReco * ptReco;
            }
          }
          float ptGen = matchedMCParticle.pt();

          float deltaPt = ptReco - ptGen;

          if (pdg == -particlePdgCodes.at(4)) {
            histomc.fill(HIST("histDeltaPtVsPtGenanti"), ptReco, deltaPt);
            histomc.fill(HIST("histPIDtrackanti"), ptReco, track.pidForTracking());
          }
          if (pdg == particlePdgCodes.at(4)) {
            histomc.fill(HIST("histDeltaPtVsPtGen"), ptReco, deltaPt);
            histomc.fill(HIST("histPIDtrack"), ptReco, track.pidForTracking());
          }
          if (pdg == -particlePdgCodes.at(5)) {
            histomc.fill(HIST("histDeltaPtVsPtGenHe4anti"), ptReco, deltaPt);
            histomc.fill(HIST("histPIDtrackantihe4"), ptReco, track.pidForTracking());
          }
          if (pdg == particlePdgCodes.at(5)) {
            histomc.fill(HIST("histDeltaPtVsPtGenHe4"), ptReco, deltaPt);
            histomc.fill(HIST("histPIDtrackhe4"), ptReco, track.pidForTracking());
          }

          // For TOF matching efficiency
          float ptTOF = -1.0; // Default: no TOF
          if (track.hasTOF()) {
            ptTOF = ptReco;
          }

          if (cfgTrackPIDsettings2->get(i, "fillsparsh") == 1) {
            histomc.fill(HIST("hSpectramc"), i, sign, collision.centFT0C(),
                         ptReco, ptTOF);
          }

          if (cfgRequirebetaplot) {
            histos.fill(HIST("Tofsignal"), getRigidity(track) * track.sign(), o2::pid::tof::Beta::GetBeta(track));
          }
        } // loop over primaryParticles types
        histos.fill(HIST("histeta"), track.eta());
      } // TrackAssoc loop
    } // Collision loop
  }

  PROCESS_SWITCH(NucleitpcPbPb, processMC, "MC reco+gen analysis with efficiency corrections", false);
  //=-=-=-==-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    constexpr float kInvalidBField = -990.f;
    auto run3grpTimestamp = bc.timestamp();
    dBz = 0;
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grpTimestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      if (bField < kInvalidBField) {
        // Fetch magnetic field from ccdb for current collision
        dBz = grpo->getNominalL3Field();
        LOG(info) << "Retrieved GRP for timestamp " << run3grpTimestamp << " with magnetic field of " << dBz << " kZG";
      } else {
        dBz = bField;
      }
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grpTimestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grpTimestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      if (bField < kInvalidBField) {
        // Fetch magnetic field from ccdb for current collision
        dBz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
        LOG(info) << "Retrieved GRP for timestamp " << run3grpTimestamp << " with magnetic field of " << dBz << " kZG";
      } else {
        dBz = bField;
      }
    }
    mRunNumber = bc.runNumber();
  }
  //----------------------------------------------------------------------------------------------------------------
  template <typename T>
  void initCollision(const T& collision)
  {
    collHasCandidate = false;
    histos.fill(HIST("histMagField"), dBz);
    histos.fill(HIST("histNev"), 0.5);
    collPassedEvSel = collision.sel8() && std::abs(collision.posZ()) < cfgZvertex;
    occupancy = collision.trackOccupancyInTimeRange();
    if (collPassedEvSel) {
      histos.fill(HIST("histNev"), 1.5);
      histos.fill(HIST("histVtxZ"), collision.posZ());
      histos.fill(HIST("histCentFT0C"), collision.centFT0C());
      histos.fill(HIST("histCentFT0M"), collision.centFT0M());
    }
    primVtx.assign({collision.posX(), collision.posY(), collision.posZ()});
    cents.assign({collision.centFT0A(), collision.centFT0C(), collision.centFT0M()});
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  void fillhmass(T const& track, int species, float cent)
  {
    if (!track.hasTOF() || !cfgFillmass)
      return;
    float beta{o2::pid::tof::Beta::GetBeta(track)};
    const float eps = 1e-6f;
    if (beta < eps || beta > 1.0f - eps)
      return;
    float charge = (species == he3 || species == he4) ? 2.f : 1.f;
    float p = getRigidity(track); // assuming this is the momentum from inner TPC
    float massTOF = p * charge * std::sqrt(1.f / (beta * beta) - 1.f);
    // get PDG mass
    float pdgMass = particleMasses[species];
    float massDiff = massTOF - pdgMass;
    float ptMomn;
    setTrackParCov(track, mTrackParCov);
    mTrackParCov.setPID(track.pidForTracking());
    ptMomn = (species == he3 || species == he4) ? 2 * mTrackParCov.getPt() : mTrackParCov.getPt();
    if (track.sign() > 0) {
      hmass[2 * species]->Fill(ptMomn, massDiff, cent);
    } else if (track.sign() < 0) {
      hmass[2 * species + 1]->Fill(ptMomn, massDiff, cent);
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  void fillhmassnsigma(T const& track, int species, float sigma)
  {
    if (!track.hasTOF() || !cfgFillmassnsigma)
      return;
    float beta{o2::pid::tof::Beta::GetBeta(track)};
    const float eps = 1e-6f;
    if (beta < eps || beta > 1.0f - eps)
      return;
    float charge = (species == he3 || species == he4) ? 2.f : 1.f;
    float p = getRigidity(track);
    float massTOF = p * charge * std::sqrt(1.f / (beta * beta) - 1.f);
    // get PDG mass
    float masssquare = massTOF * massTOF;
    if (track.sign() > 0) {
      hmassnsigma[2 * species]->Fill(sigma, masssquare);
    } else if (track.sign() < 0) {
      hmassnsigma[2 * species + 1]->Fill(sigma, masssquare);
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  float getTPCnSigma(T const& track, PrimParticles& particle)
  {
    const float rigidity = getRigidity(track);
    if (!track.hasTPC())
      return -999;
    if (particle.name == "pion" && cfgTrackPIDsettings->get("pion", "useBBparams") < 1)
      return cfgTrackPIDsettings->get("pion", "useBBparams") == 0 ? track.tpcNSigmaPi() : 0;
    if (particle.name == "proton" && cfgTrackPIDsettings->get("proton", "useBBparams") < 1)
      return cfgTrackPIDsettings->get("proton", "useBBparams") == 0 ? track.tpcNSigmaPr() : 0;
    if (particle.name == "deuteron" && cfgTrackPIDsettings->get("deuteron", "useBBparams") < 1)
      return cfgTrackPIDsettings->get("deuteron", "useBBparams") == 0 ? track.tpcNSigmaDe() : 0;
    if (particle.name == "triton" && cfgTrackPIDsettings->get("triton", "useBBparams") < 1)
      return cfgTrackPIDsettings->get("triton", "useBBparams") == 0 ? track.tpcNSigmaTr() : 0;
    if (particle.name == "helion" && cfgTrackPIDsettings->get("helion", "useBBparams") < 1)
      return cfgTrackPIDsettings->get("helion", "useBBparams") == 0 ? track.tpcNSigmaHe() : 0;
    if (particle.name == "alpha" && cfgTrackPIDsettings->get("alpha", "useBBparams") < 1)
      return cfgTrackPIDsettings->get("alpha", "useBBparams") == 0 ? track.tpcNSigmaAl() : 0;

    double expBethe{tpc::BetheBlochAleph(static_cast<double>(particle.charge * rigidity / particle.mass), particle.betheParams[0], particle.betheParams[1], particle.betheParams[2], particle.betheParams[3], particle.betheParams[4])};
    double expSigma{expBethe * particle.resolution};
    float sigmaTPC = static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
    return sigmaTPC;
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  float getITSnSigma(T const& track, PrimParticles& particle)
  {
    if (!track.hasITS())
      return -999;
    o2::aod::ITSResponse itsResponse;
    if (particle.name == "pion")
      return itsResponse.nSigmaITS<o2::track::PID::Pion>(track);
    if (particle.name == "proton")
      return itsResponse.nSigmaITS<o2::track::PID::Proton>(track);
    if (particle.name == "deuteron")
      return itsResponse.nSigmaITS<o2::track::PID::Deuteron>(track);
    if (particle.name == "triton")
      return itsResponse.nSigmaITS<o2::track::PID::Triton>(track);
    if (particle.name == "helion")
      return itsResponse.nSigmaITS<o2::track::PID::Helium3>(track);
    if (particle.name == "alpha")
      return itsResponse.nSigmaITS<o2::track::PID::Alpha>(track);
    return -999; // fallback if no match
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  float getMeanItsClsSize(T const& track)
  {
    int sum = 0, n = 0;
    constexpr int kNITSLayers = 8;
    for (int i = 0; i < kNITSLayers; i++) {
      sum += (track.itsClusterSizes() >> (4 * i) & 15);
      if (track.itsClusterSizes() >> (4 * i) & 15)
        n++;
    }
    return n > 0 ? static_cast<float>(sum) / n : 0.f;
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  float getRigidity(T const& track)
  {
    if (!cfgRigidityCorrection)
      return track.tpcInnerParam();
    bool hePID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
    return hePID ? track.tpcInnerParam() / 2 : track.tpcInnerParam();
  }
  //----------------------------------------------------------------------------------------------------------------
  float dcazSigma(double pt, float dcasigma)
  {
    float invPt = 1.f / pt;
    return (5.00000e-04 + 8.73690e-03 * invPt + 9.62329e-04 * invPt * invPt) * dcasigma; // o2-linter: disable=magic-number (To be checked)
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  float getRapidity(T const& track, int species)
  {
    using PtEtaPhiMVector = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;
    double momn;
    int speciesHe3 = 4;
    int speciesHe4 = 5;
    if (species == speciesHe3 || species == speciesHe4) {
      momn = 2 * track.pt();
    } else {
      momn = track.pt();
    }
    PtEtaPhiMVector lorentzVectorParticle(momn, track.eta(), track.phi(), particleMasses[species]);
    return lorentzVectorParticle.Rapidity();
  }
}; // end of the task here
//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NucleitpcPbPb>(cfgc)};
}
