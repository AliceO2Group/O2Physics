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
static const std::vector<std::string> correctedparticleNames{"helion", "antihelion", "alpha", "antialpha"};
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

static const int nfittingparticle = 4;
const int nfittingparameters = 4;
static const std::vector<std::string> trackcorrectionNames{"correctionneed", "a", "b", "c"};
constexpr double ktrackcorrection[nfittingparticle][nfittingparameters]{
  {1, 0.464215, 0.195771, 0.0183111}, // He3
  {1, 0.464215, 0.195771, 0.0183111}, // anti-He3
  {1, 0.00765, 0.503791, -1.10517},   // He4
  {1, 0.00765, 0.503791, -1.10517}};  // anti-He4

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

  Preslice<TracksFull> tracksPerCollision = aod::track::collisionId;

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
  Configurable<bool> cfgetaRequire{"cfgetaRequire", true, "eta cut require"};
  Configurable<bool> cfgetaRequireMC{"cfgetaRequireMC", true, "eta cut require for generated particles"};
  Configurable<bool> cfgRapidityRequireMC{"cfgRapidityRequireMC", true, "rapidity cut require for generated particles"};
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
  Configurable<bool> cfgRequirebetaplot{"cfgRequirebetaplot", true, "Require beta plot"};
  Configurable<bool> cfgdcaxynopt{"cfgdcaxynopt", true, "DCA xy cut without pT dependent"};
  Configurable<bool> cfgdcaznopt{"cfgdcaznopt", false, "DCA xy cut without pT dependent"};
  Configurable<bool> cfgmass2{"cfgmass2", true, "Fill mass square difference"};

  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {kBetheBlochDefault[0], nParticles, nBetheParams, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgTrackPIDsettings{"cfgTrackPIDsettings", {kTrackPIDSettings[0], nParticles, nTrkSettings, particleNames, trackPIDsettingsNames}, "track selection and PID criteria"};
  Configurable<LabeledArray<double>> cfgTrackPIDsettings2{"cfgTrackPIDsettings2", {kTrackPIDSettings2[0], nParticles, nTrkSettings2, particleNames, trackPIDsettingsNames2}, "track selection and PID criteria"};
  Configurable<LabeledArray<double>> cfgktrackcorrection{"cfgktrackcorrection", {ktrackcorrection[0], nfittingparticle, nfittingparameters, correctedparticleNames, trackcorrectionNames}, "fitting paramters"};
  Configurable<bool> cfgFillhspectra{"cfgFillhspectra", true, "fill data sparsh"};
  Configurable<bool> cfgFillmass{"cfgFillmass", false, "Fill mass histograms"};
  Configurable<bool> cfgFillmassnsigma{"cfgFillmassnsigma", true, "Fill mass vs nsigma histograms"};
  Configurable<float> centcut{"centcut", 80.0f, "centrality cut"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.9f, "Eta range for tracks"};
  Configurable<float> cfgCutRapidity{"cfgCutRapidity", 0.5f, "Rapidity range"};
  Configurable<float> cfgtpcNClsFindable{"cfgtpcNClsFindable", 0.8f, "tpcNClsFindable over crossedRows"};
  Configurable<float> cfgZvertex{"cfgZvertex", 10, "Min Z Vertex"};
  Configurable<bool> cfgZvertexRequireMC{"cfgZvertexRequireMC", true, "Pos Z cut in MC"};
  Configurable<bool> cfgsel8Require{"cfgsel8Require", true, "sel8 cut require"};
  Configurable<float> cfgminmassrejection{"cfgminmassrejection", 6.5, "Min side of He3 particle rejection"};
  Configurable<float> cfgmaxmassrejection{"cfgmaxmassrejection", 9.138, "Max side of He3 particle rejection"};
  Configurable<float> correctionsigma{"correctionsigma", 2, "Max sigma value outside which correction is require"};
  Configurable<bool> cfghe3massrejreq{"cfghe3massrejreq", true, "Require mass cut on He4 particles"};

  o2::track::TrackParametrizationWithError<float> mTrackParCov;
  // Binning configuration
  ConfigurableAxis axisMagField{"axisMagField", {10, -10., 10.}, "magnetic field"};
  ConfigurableAxis axisNev{"axisNev", {10, 0., 10.}, "Number of events"};
  ConfigurableAxis axisRigidity{"axisRigidity", {4000, -10., 10.}, "#it{p}^{TPC}/#it{z}"};
  ConfigurableAxis axisdEdx{"axisdEdx", {4000, 0, 4000}, "d#it{E}/d#it{x}"};
  ConfigurableAxis axisCent{"axisCent", {100, 0, 100}, "centrality"};
  ConfigurableAxis axisVtxZ{"axisVtxZ", {120, -20, 20}, "z"};

  ConfigurableAxis ptAxis{"ptAxis", {200, 0, 10}, "#it{p}_{T} (GeV/#it{c})"};

  ConfigurableAxis ptAxisa{"ptAxisa", {20, 0, 10}, "#it{p}_{T} (GeV/#it{c})"}; // just check

  ConfigurableAxis axiseta{"axiseta", {100, -1, 1}, "eta"};
  ConfigurableAxis axisrapidity{"axisrapidity", {100, -2, 2}, "rapidity"};
  ConfigurableAxis axismass{"axismass", {100, -10, 10}, "mass"};
  ConfigurableAxis axismassnsigma{"axismassnsigma", {100, 0, 20}, "nsigma mass"};
  ConfigurableAxis nsigmaAxis{"nsigmaAxis", {160, -10, 10}, "n#sigma_{#pi^{+}}"};
  ConfigurableAxis speciesBitAxis{"speciesBitAxis", {8, -0.5, 7.5}, "particle type 0: pion, 1: proton, 2: deuteron, 3: triton, 4:He3, 5:He4"};
  ConfigurableAxis speciesTrackingAxis{"speciesTrackingAxis", {11, -0.5, 10.5}, "particle type 0: pion, 1: proton, 2: deuteron, 3: triton, 4:He3, 5:He4"};
  ConfigurableAxis axisDCA{"axisDCA", {400, -10., 10.}, "DCA axis"};
  ConfigurableAxis particleAntiAxis{"particleAntiAxis", {2, -0.5, 1.5}, "Particle/Anti-particle"}; // 0 = particle, 1 = anti-particle
  ConfigurableAxis decayTypeAxis{"decayTypeAxis", {3, -0.5, 2.5}, "Decay type"};                   // 0 = primary, 1 = from decay, 2 = material

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
  bool collHasCandidate, collPassedEvSel;
  int mRunNumber, occupancy;
  float dBz;
  TRandom3 rand;
  float he3 = 4;
  float he4 = 5;
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
    histos.add("dcaZ", "dcaZ", kTH2F, {ptAxis, axisDCA});
    histos.add("dcaXY", "dcaXY", kTH2F, {ptAxis, axisDCA});
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

      // Efficiency x Acceptance
      histomc.add("hDenomEffAcc", "Denominator for Efficiency x Acceptance",
                  {HistType::kTHnSparseF, {speciesBitAxis, ptAxis, axisrapidity, axisCent, particleAntiAxis, decayTypeAxis}});
      histomc.add("hNumerEffAcc", "Numerator for Efficiency x Acceptance",
                  {HistType::kTHnSparseF, {speciesBitAxis, ptAxis, axisrapidity, axisCent, particleAntiAxis, decayTypeAxis}});

      // The Signal loss correction
      histomc.add("hHe3SignalLossDenom", "He3 Signal Loss Denominator", kTH1F, {axisCent});
      histomc.add("hHe3SignalLossNumer", "He3 Signal Loss Numerator", kTH1F, {axisCent});
      histomc.add("hHe4SignalLossDenom", "He4 Signal Loss Denominator", kTH1F, {axisCent});
      histomc.add("hHe4SignalLossNumer", "He4 Signal Loss Numerator", kTH1F, {axisCent});

      histomc.add("haHe3SignalLossDenom", "He3 Signal Loss Denominator", kTH1F, {axisCent});
      histomc.add("haHe3SignalLossNumer", "He3 Signal Loss Numerator", kTH1F, {axisCent});
      histomc.add("haHe4SignalLossDenom", "He4 Signal Loss Denominator", kTH1F, {axisCent});
      histomc.add("haHe4SignalLossNumer", "He4 Signal Loss Numerator", kTH1F, {axisCent});

      // The event loss correction
      histomc.add("hEventLossDenom", "Event loss denominator", kTH1F, {axisCent});
      histomc.add("hEventLossNumer", "Event loss numerator", kTH1F, {axisCent});

      histomc.add("histVtxZgen", "histVtxZgen", kTH1F, {axisVtxZ});
      histomc.add("histNevReco", "histNevReco", kTH1F, {axisNev});
      histomc.add("histVtxZReco", "histVtxZReco", kTH1F, {axisVtxZ});
      histomc.add("histCentFT0CReco", "histCentFT0CReco", kTH1F, {axisCent});
      histomc.add("histCentFT0MReco", "histCentFT0MReco", kTH1F, {axisCent});

      histomc.add("histdetapttriton", " delta pt vs pt rec for trition detected", HistType::kTH2F, {{1000, 0, 10}, {1000, -0.5, 0.5, "p_{T}(reco) - p_{T}(gen);p_{T}(reco)"}});
      histomc.add("histdetapttritonanti", " delta pt vs pt rec for trition detected", HistType::kTH2F, {{1000, 0, 10}, {1000, -0.5, 0.5, "p_{T}(reco) - p_{T}(gen);p_{T}(reco)"}});

      histomc.add("histDeltaPtVsPtGen", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10}, {1000, -0.5, 0.5, "p_{T}(reco) - p_{T}(gen);p_{T}(reco)"}});
      histomc.add("histDeltaPtVsPtGenanti", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10}, {1000, -0.5, 0.5, "p_{T}(reco) - p_{T}(gen);p_{T}(reco)"}});
      histomc.add("histDeltaPtVsPtGenHe4", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10}, {1000, -0.5, 0.5, "p_{T}(reco) - p_{T}(gen);p_{T}(reco)"}});
      histomc.add("histDeltaPtVsPtGenHe4anti", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10}, {1000, -0.5, 0.5, "p_{T}(reco) - p_{T}(gen);p_{T}(reco)"}});
      histomc.add("histPIDtrack", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10, "p_{T}(reco)"}, {9, -0.5, 8.5, "p_{T}(reco) - p_{T}(gen)"}});
      histomc.add("histPIDtrackanti", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10, "p_{T}(reco)"}, {9, -0.5, 8.5, "p_{T}(reco) - p_{T}(gen)"}});
      histomc.add("histPIDtrackhe4", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10, "p_{T}(reco)"}, {9, -0.5, 8.5, "p_{T}(reco) - p_{T}(gen)"}});
      histomc.add("histPIDtrackantihe4", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10, "p_{T}(reco)"}, {9, -0.5, 8.5, "p_{T}(reco) - p_{T}(gen)"}});
    }

    if (doprocessDCA) {

      histomc.add<THnSparse>("hSpectraDCA", " ", HistType::kTHnSparseF, {speciesBitAxis, {5, -2.5, 2.5}, axisCent, ptAxis, ptAxis, decayTypeAxis, axisDCA});
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------------------------------------
  void processData(CollisionsFull const& collisions,
                   TracksFull const& tracks,
                   aod::BCsWithTimestamps const&)
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
      histos.fill(HIST("histNev"), 2.5);
      if (removeNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup))
        continue;
      histos.fill(HIST("histNev"), 3.5);

      if (requireIsGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
        continue;
      histos.fill(HIST("histNev"), 4.5);
      if (requireIsVertexITSTPC && !collision.selection_bit(aod::evsel::kIsVertexITSTPC))
        continue;

      histos.fill(HIST("histNev"), 5.5);
      if (removeNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
        continue;
      histos.fill(HIST("histNev"), 6.5);
      histos.fill(HIST("histCentFTOC_cut"), collision.centFT0C());

      // new slicing
      auto tracksInColl = tracks.sliceBy(tracksPerCollision, collision.globalIndex());

      // loop over sliced tracks
      for (const auto& track : tracksInColl) {
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

          float tpcNsigma = getTPCnSigma(track, primaryParticles.at(i));
          if ((std::abs(tpcNsigma) > cfgTrackPIDsettings->get(i, "maxTPCnSigma")) && cfgmaxTPCnSigmaRequire)
            continue;
          if (std::abs(tpcNsigma) > correctionsigma) {
            double a = 0, b = 0, c = 0;

            int param = -1;
            if (i == he3) {
              param = (track.sign() > 0) ? 0 : 1;
            } else if (i == he4) {
              param = (track.sign() > 0) ? 2 : 3;
            }

            if (param >= 0) {
              a = cfgktrackcorrection->get(param, "a");
              b = cfgktrackcorrection->get(param, "b");
              c = cfgktrackcorrection->get(param, "c");
            }

            if (i == he4 && cfgmccorrectionhe4Require) {
              ptMomn = ptMomn + a + b * std::exp(c * ptMomn);
            }

            if (i == he3 && cfgmccorrectionhe4Require) {
              int pidGuess = track.pidForTracking();
              int antitriton = 6;
              if (pidGuess == antitriton) {
                ptMomn = ptMomn - a + b * ptMomn - c * ptMomn * ptMomn;
              }
            }
          }
          int sign = (track.sign() > 0) ? 1 : ((track.sign() < 0) ? -1 : 0);

          if (std::abs(getRapidity(track, i)) > cfgCutRapidity && cfgRapidityRequire)
            continue;
          if (track.tpcNClsFound() < cfgTrackPIDsettings->get(i, "minTPCnCls") && cfgTPCNClsfoundRequire)
            continue;
          if (((track.tpcNClsCrossedRows() < cfgTrackPIDsettings->get(i, "minTPCnClsCrossedRows")) ||
               track.tpcNClsCrossedRows() < cfgtpcNClsFindable * track.tpcNClsFindable()) &&
              cfgTPCNClsCrossedRowsRequire)
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

          // DCA XY cut
          bool insideDCAxy = cfgdcaxynopt ? (std::abs(track.dcaXY()) <= cfgTrackPIDsettings->get(i, "maxDcaXY")) : (std::abs(track.dcaXY()) <= (cfgTrackPIDsettings->get(i, "maxDcaXY") * (0.0105f + 0.0350f / std::pow(ptMomn, 1.1f))));

          // DCA Z cut
          bool insideDCAz = cfgdcaznopt ? (std::abs(track.dcaZ()) <= cfgTrackPIDsettings->get(i, "maxDcaZ")) : (std::abs(track.dcaZ()) <= dcazSigma(ptMomn, cfgTrackPIDsettings->get(i, "maxDcaZ")));

          if ((!insideDCAxy || !insideDCAz)) {
            continue;
          }

          float itsSigma = getITSnSigma(track, primaryParticles.at(i));
          if (itsSigma < cfgTrackPIDsettings2->get(i, "minITSnsigma") && cfgTrackPIDsettings2->get(i, "useITSnsigma") < 1)
            continue;
          if (itsSigma > cfgTrackPIDsettings2->get(i, "maxITSnsigma") && cfgTrackPIDsettings2->get(i, "useITSnsigma") < 1)
            continue;

          histos.fill(HIST("Tpcsignal"), getRigidity(track) * track.sign(), track.tpcSignal());

          histos.fill(HIST("dcaXY"), ptMomn, track.dcaXY());
          histos.fill(HIST("dcaZ"), ptMomn, track.dcaZ());

          if (cfgFillhspectra && cfgTrackPIDsettings2->get(i, "fillsparsh") == 1) {

            if (i != he4) {
              histos.fill(HIST("hSpectra"), i, ptMomn, tpcNsigma, sign, collision.centFT0C(), track.dcaZ(), track.dcaXY());
            } else {
              if (!track.hasTOF()) {
                // Fill without TOF
                histos.fill(HIST("hSpectra"), i, ptMomn, tpcNsigma, sign, collision.centFT0C(), track.dcaZ(), track.dcaXY());
              } else {
                // Has TOF - apply mass cut
                float beta = o2::pid::tof::Beta::GetBeta(track);
                const float eps = 1e-6f;
                if (beta < eps || beta > 1.0f - eps)
                  continue;

                float charge = 2.f; // he4 has charge 2
                float p = getRigidity(track);
                float massTOF = p * charge * std::sqrt(1.f / (beta * beta) - 1.f);

                // Apply mass cut for he4 (mass^2 around 3.73^2 = 13.9)
                if (cfghe3massrejreq && (massTOF * massTOF > cfgminmassrejection && massTOF * massTOF < cfgmaxmassrejection)) {
                  continue; // Skip if mass cut fails
                }

                histos.fill(HIST("hSpectra"), i, ptMomn, tpcNsigma, sign, collision.centFT0C(), track.dcaZ(), track.dcaXY());
              }
            }
          }
          fillhmassnsigma(track, i, tpcNsigma);

          if ((std::abs(tpcNsigma) > cfgTrackPIDsettings2->get(i, "maxTPCnsigmaTOF")) && cfgTrackPIDsettings2->get(i, "useTPCnsigmaTOF") < 1)
            continue;
          fillhmass(track, i, collision.centFT0C());

          if (cfgRequirebetaplot) {
            histos.fill(HIST("Tofsignal"), getRigidity(track) * track.sign(), o2::pid::tof::Beta::GetBeta(track));
          }
        } // loop primaryParticles

        histos.fill(HIST("histeta"), track.eta());
      } // loop sliced tracks
    } // collision loop
  }
  PROCESS_SWITCH(NucleitpcPbPb, processData, "data analysis", false);

  //----------------------------------------------------------------------------------------------------------------
  // MC particles - Efficiency x Acceptance and Signal Loss calculations
  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  struct McCollInfo {
    bool passedEvSel = false;
    float centrality = -1.0f;
    bool passedEvSelVtZ = false;
  };
  std::vector<McCollInfo> mcCollInfos;

  void processMC(CollisionsFullMC const& collisions,
                 aod::McCollisions const& mcCollisions,
                 soa::Join<TracksFull, aod::McTrackLabels> const& tracks,
                 aod::McParticles const& particlesMC,
                 aod::BCsWithTimestamps const&)

  {

    mcCollInfos.clear();
    mcCollInfos.resize(mcCollisions.size());

    // First pass: Store centrality and apply event selection
    for (auto const& collision : collisions) {
      int mcCollIdx = collision.mcCollisionId();
      if (mcCollIdx < 0 || mcCollIdx >= static_cast<int>(mcCollisions.size())) {
        continue;
      }

      // STORE CENTRALITY WITHOUt CUTS
      mcCollInfos[mcCollIdx].centrality = collision.centFT0C();

      if (!collision.sel8() && cfgsel8Require)
        continue;
      if (collision.centFT0C() > centcut)
        continue;

      // Additional cuts
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

      // Mark this MC collision as passing event selection
      mcCollInfos[mcCollIdx].passedEvSel = true;

      // Apply event selection cuts
      if (std::abs(collision.posZ()) > cfgZvertex && cfgZvertexRequireMC)
        continue;

      mcCollInfos[mcCollIdx].passedEvSelVtZ = true;
    }

    // FILL EVENT LOSS AND SIGNAL LOSS: Combined loop per MC collision
    for (size_t i = 0; i < mcCollInfos.size(); i++) {
      if (mcCollInfos[i].centrality >= 0) { // Only if we found a matching collision
        // Event loss denominator
        histomc.fill(HIST("hEventLossDenom"), mcCollInfos[i].centrality);

        // Event loss numerator (if passed selection)
        if (mcCollInfos[i].passedEvSel) {
          histomc.fill(HIST("hEventLossNumer"), mcCollInfos[i].centrality);
        }

        // Fill signal loss for all primary particles in this MC collision
        for (auto const& mcParticle : particlesMC) {
          if (mcParticle.mcCollisionId() != static_cast<int>(i)) {
            continue;
          }
          if (!mcParticle.isPhysicalPrimary()) {
            continue;
          }

          // Signal loss denominator
          if (mcParticle.pdgCode() == particlePdgCodes.at(4)) { // He3
            histomc.fill(HIST("hHe3SignalLossDenom"), mcCollInfos[i].centrality);

          } else if (mcParticle.pdgCode() == particlePdgCodes.at(5)) { // He4
            histomc.fill(HIST("hHe4SignalLossDenom"), mcCollInfos[i].centrality);
          } else if (mcParticle.pdgCode() == -particlePdgCodes.at(4)) { // anti-He3
            histomc.fill(HIST("haHe3SignalLossDenom"), mcCollInfos[i].centrality);

          } else if (mcParticle.pdgCode() == -particlePdgCodes.at(5)) { // He4
            histomc.fill(HIST("haHe4SignalLossDenom"), mcCollInfos[i].centrality);
          }

          // Signal loss numerator (if event passed selection)
          if (mcCollInfos[i].passedEvSel) {
            if (mcParticle.pdgCode() == particlePdgCodes.at(4)) { // He3
              histomc.fill(HIST("hHe3SignalLossNumer"), mcCollInfos[i].centrality);
            } else if (mcParticle.pdgCode() == particlePdgCodes.at(5)) { // He4
              histomc.fill(HIST("hHe4SignalLossNumer"), mcCollInfos[i].centrality);
            } else if (mcParticle.pdgCode() == -particlePdgCodes.at(4)) { // anti-He3
              histomc.fill(HIST("haHe3SignalLossNumer"), mcCollInfos[i].centrality);
            } else if (mcParticle.pdgCode() == -particlePdgCodes.at(5)) { // anti-He4
              histomc.fill(HIST("haHe4SignalLossNumer"), mcCollInfos[i].centrality);
            }
          }
        }
      }
    }

    // Process MC collisions for efficiency and reconstructed collisions
    for (auto const& mcCollision : mcCollisions) {
      size_t idx = mcCollision.globalIndex();
      if (idx >= mcCollInfos.size())
        continue;

      float centrality = mcCollInfos[idx].centrality;
      //      bool passedEvSel = mcCollInfos[idx].passedEvSel;
      bool passedEvSelVtZ = mcCollInfos[idx].passedEvSelVtZ;

      // Process generated particles for efficiency denominators
      histomc.fill(HIST("histVtxZgen"), mcCollision.posZ());

      for (auto const& mcParticle : particlesMC) {
        if (mcParticle.mcCollisionId() != mcCollision.globalIndex())
          continue;

        int pdgCode = mcParticle.pdgCode();
        bool isHe3 = (std::abs(pdgCode) == particlePdgCodes.at(4));
        bool isHe4 = (std::abs(pdgCode) == particlePdgCodes.at(5));

        if (!isHe3 && !isHe4)
          continue;

        if (std::abs(mcParticle.eta()) > cfgCutEta && cfgetaRequireMC)
          continue;
        if (std::abs(mcParticle.y()) > cfgCutRapidity && cfgRapidityRequireMC)
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
          //  if (!mcParticle.isPhysicalPrimary())
          continue;

        int particleType = -1;
        if (std::abs(pdgCode) == particlePdgCodes.at(4))
          particleType = he3;
        else if (std::abs(pdgCode) == particlePdgCodes.at(5))
          particleType = he4;

        if (particleType >= 0) {

          // Efficiency x Acceptance histograms
          if (passedEvSelVtZ) {
            histomc.fill(HIST("hDenomEffAcc"), particleType, mcParticle.pt(), mcParticle.y(), centrality, particleAnti, decayType);
          }
        }
      }

      // Process reconstructed collisions for this MC collision
      if (passedEvSelVtZ) {
        // Find the corresponding reconstructed collision
        for (auto const& collision : collisions) {
          if (collision.mcCollisionId() != static_cast<int>(idx))
            continue;

          auto bc = collision.bc_as<aod::BCsWithTimestamps>();
          initCCDB(bc);

          histomc.fill(HIST("histNevReco"), 0.5);
          histomc.fill(HIST("histVtxZReco"), collision.posZ());
          histomc.fill(HIST("histCentFT0CReco"), collision.centFT0C());
          histomc.fill(HIST("histCentFT0MReco"), collision.centFT0M());

          auto tracksInColl = tracks.sliceBy(tracksPerCollision, collision.globalIndex());

          for (auto const& track : tracksInColl) {
            if (!track.has_mcParticle())
              continue; // skip un-matched reco tracks

            auto const& matchedMCParticle = track.mcParticle_as<aod::McParticles>();

            // Only process particles from this MC collision
            if (matchedMCParticle.mcCollisionId() != mcCollision.globalIndex())
              continue;

            int pdg = matchedMCParticle.pdgCode();
            bool isHe3 = (std::abs(pdg) == particlePdgCodes.at(4));
            bool isHe4 = (std::abs(pdg) == particlePdgCodes.at(5));

            if (!isHe3 && !isHe4)
              continue;

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
            if (!matchedMCParticle.isPhysicalPrimary() && isFromWeakDecay)
              continue;

            for (size_t i = 0; i < primaryParticles.size(); i++) {
              if (std::abs(pdg) != std::abs(particlePdgCodes.at(i)))
                continue;

              float ptReco;
              setTrackParCov(track, mTrackParCov);
              mTrackParCov.setPID(track.pidForTracking());

              ptReco = (std::abs(pdg) == particlePdgCodes.at(4) || std::abs(pdg) == particlePdgCodes.at(5)) ? 2 * mTrackParCov.getPt() : mTrackParCov.getPt();

              int particleAnti = (pdg > 0) ? 0 : 1;

              double a = 0, b = 0, c = 0;

              int param = -1;
              if (i == he3) {
                param = (-particlePdgCodes.at(4) > 0) ? 0 : 1;
              } else if (i == he4) {
                param = (-particlePdgCodes.at(4) > 0) ? 2 : 3;
              }

              if (param >= 0) {
                a = cfgktrackcorrection->get(param, "a");
                b = cfgktrackcorrection->get(param, "b");
                c = cfgktrackcorrection->get(param, "c");
              }

              if (std::abs(pdg) == particlePdgCodes.at(5) && cfgmccorrectionhe4Require) {
                ptReco = ptReco + a + b * std::exp(c * ptReco);
              }

              if (std::abs(pdg) == particlePdgCodes.at(4) && cfgmccorrectionhe4Require) {
                int pidGuess = track.pidForTracking();
                int antitriton = 6;
                if (pidGuess == antitriton) {
                  ptReco = ptReco - a + b * ptReco - c * ptReco * ptReco;
                }
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

              // DCA XY cut
              bool insideDCAxy = cfgdcaxynopt ? (std::abs(track.dcaXY()) <= cfgTrackPIDsettings->get(i, "maxDcaXY")) : (std::abs(track.dcaXY()) <= (cfgTrackPIDsettings->get(i, "maxDcaXY") * (0.0105f + 0.0350f / std::pow(ptReco, 1.1f))));

              // DCA Z cut
              bool insideDCAz = cfgdcaznopt ? (std::abs(track.dcaZ()) <= cfgTrackPIDsettings->get(i, "maxDcaZ")) : (std::abs(track.dcaZ()) <= dcazSigma(ptReco, cfgTrackPIDsettings->get(i, "maxDcaZ")));

              if ((!insideDCAxy || !insideDCAz)) {
                continue;
              }

              float tpcNsigma = getTPCnSigma(track, primaryParticles.at(i));
              if ((std::abs(tpcNsigma) > cfgTrackPIDsettings->get(i, "maxTPCnSigma")) && cfgmaxTPCnSigmaRequire)
                continue;

              if (i == he3 || i == he4) {
                histomc.fill(HIST("hNumerEffAcc"), i, ptReco, getRapidity(track, i), collision.centFT0C(), particleAnti, decayType);
              }

              float ptTOF = -1.0; // Default: no TOF
              if (track.hasTOF()) {
                ptTOF = ptReco;
              }

              if (cfgTrackPIDsettings2->get(i, "fillsparsh") == 1) {
                histomc.fill(HIST("hSpectramc"), i, particleAnti, collision.centFT0C(),
                             ptReco, ptTOF);
              }

              fillhmassnsigma(track, i, tpcNsigma);
              histos.fill(HIST("dcaXY"), ptReco, track.dcaXY());
              histos.fill(HIST("dcaZ"), ptReco, track.dcaZ());

              histos.fill(HIST("Tpcsignal"), getRigidity(track) * track.sign(), track.tpcSignal());

              // Fill the requested histograms
              float ptGen = matchedMCParticle.pt();
              float deltaPt = ptReco - ptGen;

              if (pdg == -particlePdgCodes.at(4)) {
                histomc.fill(HIST("histDeltaPtVsPtGenanti"), ptReco, deltaPt);
                histomc.fill(HIST("histPIDtrackanti"), ptReco, track.pidForTracking());

                int pidGuess = track.pidForTracking();
                int antitriton = 6;
                if (pidGuess == antitriton) {
                  histomc.fill(HIST("histdetapttritonanti"), ptReco, deltaPt);
                }
              }
              if (pdg == particlePdgCodes.at(4)) {
                histomc.fill(HIST("histDeltaPtVsPtGen"), ptReco, deltaPt);
                histomc.fill(HIST("histPIDtrack"), ptReco, track.pidForTracking());

                int pidGuess = track.pidForTracking();
                int antitriton = 6;
                if (pidGuess == antitriton) {
                  histomc.fill(HIST("histdetapttriton"), ptReco, deltaPt);
                }
              }
              if (pdg == -particlePdgCodes.at(5)) {
                histomc.fill(HIST("histDeltaPtVsPtGenHe4anti"), ptReco, deltaPt);
                histomc.fill(HIST("histPIDtrackantihe4"), ptReco, track.pidForTracking());
              }
              if (pdg == particlePdgCodes.at(5)) {
                histomc.fill(HIST("histDeltaPtVsPtGenHe4"), ptReco, deltaPt);
                histomc.fill(HIST("histPIDtrackhe4"), ptReco, track.pidForTracking());
              }
            }
          }
          break; // Found the matching collision, break out of collision loop
        }
      }
    }
  }
  PROCESS_SWITCH(NucleitpcPbPb, processMC, "MC reco+gen analysis with efficiency corrections", false);
  //=-=-=-==-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  //----------------------------------------------------------------------------------------------------------------
  // MC particles - DCA secondary fraction
  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

  void processDCA(CollisionsFullMC const& collisions,
                  aod::McCollisions const& mcCollisions,
                  aod::McParticles const& particlesMC,
                  soa::Join<TracksFull, aod::McTrackLabels> const& tracks,
                  aod::BCsWithTimestamps const&)

  {
    (void)particlesMC;
    mcCollInfos.clear();
    mcCollInfos.resize(mcCollisions.size());

    // First pass: Store centrality and apply event selection
    for (auto const& collision : collisions) {
      int mcCollIdx = collision.mcCollisionId();
      if (mcCollIdx < 0 || mcCollIdx >= static_cast<int>(mcCollisions.size())) {
        continue;
      }

      // STORE CENTRALITY WITHOUt CUTS
      mcCollInfos[mcCollIdx].centrality = collision.centFT0C();

      if (!collision.sel8() && cfgsel8Require)
        continue;
      if (collision.centFT0C() > centcut)
        continue;

      // Additional cuts
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

      // Mark this MC collision as passing event selection
      mcCollInfos[mcCollIdx].passedEvSel = true;

      // Apply event selection cuts
      if (std::abs(collision.posZ()) > cfgZvertex && cfgZvertexRequireMC)
        continue;

      mcCollInfos[mcCollIdx].passedEvSelVtZ = true;
    }

    // Process MC collisions for efficiency and reconstructed collisions
    for (auto const& mcCollision : mcCollisions) {
      size_t idx = mcCollision.globalIndex();
      if (idx >= mcCollInfos.size())
        continue;

      //      bool passedEvSel = mcCollInfos[idx].passedEvSel;
      bool passedEvSelVtZ = mcCollInfos[idx].passedEvSelVtZ;

      // Process reconstructed collisions for this MC collision
      if (passedEvSelVtZ) {
        // Find the corresponding reconstructed collision
        for (auto const& collision : collisions) {
          if (collision.mcCollisionId() != static_cast<int>(idx))
            continue;

          auto bc = collision.bc_as<aod::BCsWithTimestamps>();
          initCCDB(bc);
          auto tracksInColl = tracks.sliceBy(tracksPerCollision, collision.globalIndex());

          for (auto const& track : tracksInColl) {
            if (!track.has_mcParticle())
              continue; // skip un-matched reco tracks

            auto const& matchedMCParticle = track.mcParticle_as<aod::McParticles>();

            // Only process particles from this MC collision
            if (matchedMCParticle.mcCollisionId() != mcCollision.globalIndex())
              continue;

            int pdg = matchedMCParticle.pdgCode();
            bool isHe3 = (std::abs(pdg) == particlePdgCodes.at(4));
            bool isHe4 = (std::abs(pdg) == particlePdgCodes.at(5));

            if (!isHe3 && !isHe4)
              continue;

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
              if (std::abs(pdg) != std::abs(particlePdgCodes.at(i)))
                continue;

              float ptReco;
              setTrackParCov(track, mTrackParCov);
              mTrackParCov.setPID(track.pidForTracking());

              ptReco = (std::abs(pdg) == particlePdgCodes.at(4) || std::abs(pdg) == particlePdgCodes.at(5)) ? 2 * mTrackParCov.getPt() : mTrackParCov.getPt();

              int particleAnti = (pdg > 0) ? 0 : 1;

              double a = 0, b = 0, c = 0;

              int param = -1;
              if (i == he3) {
                param = (-particlePdgCodes.at(4) > 0) ? 0 : 1;
              } else if (i == he4) {
                param = (-particlePdgCodes.at(4) > 0) ? 2 : 3;
              }

              if (param >= 0) {
                a = cfgktrackcorrection->get(param, "a");
                b = cfgktrackcorrection->get(param, "b");
                c = cfgktrackcorrection->get(param, "c");
              }

              if (std::abs(pdg) == particlePdgCodes.at(5) && cfgmccorrectionhe4Require) {
                ptReco = ptReco + a + b * std::exp(c * ptReco);
              }

              if (std::abs(pdg) == particlePdgCodes.at(4) && cfgmccorrectionhe4Require) {
                int pidGuess = track.pidForTracking();
                int antitriton = 6;
                if (pidGuess == antitriton) {
                  ptReco = ptReco - a + b * ptReco - c * ptReco * ptReco;
                }
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

              // DCA XY cut
              bool insideDCAxy = cfgdcaxynopt ? (std::abs(track.dcaXY()) <= cfgTrackPIDsettings->get(i, "maxDcaXY")) : (std::abs(track.dcaXY()) <= (cfgTrackPIDsettings->get(i, "maxDcaXY") * (0.0105f + 0.0350f / std::pow(ptReco, 1.1f))));

              // DCA Z cut
              bool insideDCAz = cfgdcaznopt ? (std::abs(track.dcaZ()) <= cfgTrackPIDsettings->get(i, "maxDcaZ")) : (std::abs(track.dcaZ()) <= dcazSigma(ptReco, cfgTrackPIDsettings->get(i, "maxDcaZ")));

              if ((!insideDCAxy || !insideDCAz)) {
                continue;
              }

              float ptTOF = -1.0; // Default: no TOF
              if (track.hasTOF()) {
                ptTOF = ptReco;
              }

              int decayType = 0; // 0 = primary, 1 = weak decay, 2 = material

              bool isProdByGen = false;
              isProdByGen = track.mcParticle().producedByGenerator();

              if (matchedMCParticle.isPhysicalPrimary()) {
                // ---- Primary particles ----
                decayType = 0;

              } else if (matchedMCParticle.getProcess() == TMCProcess::kPDecay && !isProdByGen) {
                // ---- Secondary from weak decay ----
                decayType = 1;

              } else if (matchedMCParticle.getProcess() == 23) {
                // ---- Secondary from material interaction ----
                decayType = 2;
              }

              if (cfgTrackPIDsettings2->get(i, "fillsparsh") == 1) {
                histomc.fill(HIST("hSpectraDCA"), i, particleAnti, collision.centFT0C(),
                             ptReco, ptTOF, decayType, track.dcaXY());
              }

              //
            }
          }
          break; // Found the matching collision, break out of collision loop
        }
      }
    }
  }
  PROCESS_SWITCH(NucleitpcPbPb, processDCA, "MC DCA analysis For secondary correction", false);
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
    float massDiff = 0.0;
    if (species != he4) {
      if (cfgmass2) {
        // Compare squared masses
        massDiff = massTOF * massTOF - pdgMass * pdgMass;
      } else {
        // Compare linear masses
        massDiff = massTOF - pdgMass;
      }
    }
    if (species == he4) {
      if (cfghe3massrejreq && (massTOF * massTOF > cfgminmassrejection && massTOF * massTOF < cfgmaxmassrejection))
        return;
      if (cfgmass2) {
        // Compare squared masses
        massDiff = massTOF * massTOF - pdgMass * pdgMass;
      } else {
        // Compare linear masses
        massDiff = massTOF - pdgMass;
      }
    }

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

    if (species != he4) {
      masssquare = massTOF * massTOF;
    }
    if (species == he4) {
      if (cfghe3massrejreq && (massTOF * massTOF > cfgminmassrejection && massTOF * massTOF < cfgmaxmassrejection))
        return;
      masssquare = massTOF * massTOF;
    }

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
