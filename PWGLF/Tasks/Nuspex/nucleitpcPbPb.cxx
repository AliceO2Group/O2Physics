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

/// \file nucleitpcPbPb.cxx
/// \brief Nuclei (pion, proton, deuteron, triton, He3, He4) TPC(+TOF) PID
///        analysis for Pb-Pb, including MC efficiency/acceptance, signal-loss,
///        event-loss and DCA-template corrections.
/// \note under work
///
/// \author Jaideep Tanwar <jaideep.tanwar@cern.ch>, Panjab University

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseITS.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DataFormatsParameters/GRPObject.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MathUtils/BetheBlochAleph.h>
#include <ReconstructionDataFormats/TrackParametrizationWithError.h>

#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <TH3.h>
#include <THnSparse.h>
#include <TMCProcess.h>
#include <TRandom3.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using CollisionsFull = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, o2::aod::TracksDCA, aod::TrackSelectionExtension, aod::pidTPCPi, aod::pidTPCPr, aod::pidTPCDe, aod::pidTPCTr, aod::pidTPCHe, aod::pidTPCAl, aod::TOFSignal, aod::TOFEvTime, aod::pidTOFFullPi, aod::pidTOFFullPr, aod::pidTOFFullDe, aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl>;
using CollisionsFullMC = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::EvSels, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As>;

//===================================================================================
// Static configuration tables (species definitions, PID/quality defaults, and
// track/DCA correction parameters). None of the values below are changed by
// this cleanup -- only the surrounding comments/formatting are.
//===================================================================================
namespace
{

// ---- Particle species definitions -------------------------------------------
constexpr int NParticles = 6; // pion, proton, deuteron, triton, helion(He3), alpha(He4)
const std::vector<std::string> particleNames{"pion", "proton", "deuteron", "triton", "helion", "alpha"};
const std::vector<std::string> correctedparticleNames{"helion", "antihelion", "alpha", "antialpha"};
const std::vector<int> particlePdgCodes{211, 2212, o2::constants::physics::kDeuteron, o2::constants::physics::kTriton, o2::constants::physics::kHelium3, o2::constants::physics::kAlpha};
const std::vector<double> particleMasses{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron, o2::constants::physics::MassTriton, o2::constants::physics::MassHelium3, o2::constants::physics::MassAlpha};
const std::vector<int> particleCharge{1, 1, 1, 1, 2, 2};

// ---- TPC Bethe-Bloch parameterisation defaults (per species) ---------------
constexpr int NBetheParams = 6;
const std::vector<int> hfMothCodes = {511, 521, 531, 541, 5122}; // b-mesons + Lambda_b
const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
const std::array<std::array<double, NBetheParams>, NParticles> kBetheBlochDefault{{
  {13.611469, 3.598765, -0.021138, 2.039562, 0.651040, 0.09},     // pion
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},       // proton
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},       // deuteron
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},       // triton
  {-126.557359, -0.858569, 1.111643, 1.210323, 2.656374, 0.09},   // helion
  {-126.557359, -0.858569, 1.111643, 1.210323, 2.656374, 0.09}}}; // alpha

// ---- Track PID/quality selection defaults (per species) --------------------
constexpr int NTrkSettings = 13;
const std::vector<std::string> trackPIDsettingsNames{"useBBparams", "minITSnCls", "minITSnClscos", "minTPCnCls", "maxTPCchi2", "minTPCchi2", "maxITSchi2", "maxTPCnSigma", "maxDcaXY", "maxDcaZ", "minITSclsSize", "minTPCnClsCrossedRows", "minReqClusterITSib"};
const std::array<std::array<double, NTrkSettings>, NParticles> kTrackPIDSettings{{
  {0, 0, 4, 60, 4.0, 0.5, 100, 2.5, 2., 2., 0., 70, 1},
  {1, 0, 4, 70, 4.0, 0.5, 100, 3.0, 2., 2., 0., 70, 1},
  {1, 0, 4, 70, 4.0, 0.5, 100, 3.0, 2., 2., 0., 70, 1},
  {1, 0, 4, 70, 4.0, 0.5, 100, 3.0, 2., 2., 0., 70, 1},
  {1, 0, 4, 75, 4.0, 0.5, 100, 5.0, 2., 2., 0., 70, 1},
  {1, 0, 4, 70, 4.0, 0.5, 100, 5.0, 2., 2., 0., 70, 1}}};

constexpr int NTrkSettings2 = 7;
const std::vector<std::string> trackPIDsettingsNames2{"useITSnsigma", "minITSnsigma", "maxITSnsigma", "fillsparsh", "useTPCnsigmaTOF", "maxTPCnsigmaTOF", "maxTOFnsigma"};
const std::array<std::array<double, NTrkSettings2>, NParticles> kTrackPIDSettings2{{
  {1, -5, 4, 0, 1, 2, 2},
  {1, -5, 4, 0, 1, 2, 2},
  {1, -5, 4, 0, 1, 2, 2},
  {1, -5, 4, 1, 1, 2, 2},
  {1, -5, 4, 1, 1, 2, 2},
  {1, -5, 4, 1, 1, 2, 2}}};

// ---- pT-dependent reconstruction correction parameters for He3/He4 ---------
constexpr int NFittingParticle = 4;
constexpr int NFittingParameters = 4;
const std::vector<std::string> trackcorrectionNames{"a", "b", "c", "d"};
const std::array<std::array<double, NFittingParameters>, NFittingParticle> ktrackcorrection{{
  {0.464215, 0.195771, 0.0183111, 0.0},  // He3
  {0.464215, 0.195771, 0.0183111, 0.0},  // anti-He3
  {0.00765, 0.503791, -1.10517, 0.0},    // He4
  {0.00765, 0.503791, -1.10517, 0.0}}};  // anti-He4

// ---- pT-dependent DCA cut parameters ----------------------------------------
//   DCAxy(pt) = p0 * exp( p1 * pt) + p2
//   DCAz(pt)  = p0 * exp(-p1 * pt) + p2
constexpr int NDCATypes = 2;
constexpr int NDCAParameters = 3;
const std::vector<std::string> dcaNames{"p0", "p1", "p2"};
const std::vector<std::string> dcaTypeNames{"DCAxy", "DCAz"};
const std::array<std::array<double, NDCAParameters>, NDCATypes> kDCAcorrection{{
  {0.0118, -0.6889, 0.0017},
  {0.1014, 1.7512, 0.0024}}};

// ---- MC decay-origin classification -----------------------------------------
constexpr int DecayTypePrimary = 0;
constexpr int DecayTypeWeak = 1;
constexpr int DecayTypeTransport = 2;

/// One row of the Bethe-Bloch/PID configuration table, resolved for a given species.
struct PrimParticles {
  TString name;
  int pdgCode;
  int charge;
  double mass;
  double resolution;
  std::vector<double> betheParams;
  bool active;

  PrimParticles(const std::string& name_, int pdgCode_, double mass_, int charge_, const LabeledArray<double>& bethe)
    : name(name_), pdgCode(pdgCode_), charge(charge_), mass(mass_), active(false)
  {
    resolution = bethe.get(name, "resolution");
    betheParams.clear();
    constexpr unsigned int NSpecies = 5;
    for (unsigned int i = 0; i < NSpecies; ++i) {
      betheParams.push_back(bethe.get(name, i));
    }
  }
}; // struct PrimParticles

std::vector<std::shared_ptr<TH3>> hmass;
std::vector<std::shared_ptr<TH3>> hmassnsigma;

} // namespace

//===================================================================================
struct NucleitpcPbPb {

  Preslice<TracksFull> tracksPerCollision = aod::track::collisionId;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histomc{"histomc", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // ---- Event selection configurables ------------------------------------
  Configurable<bool> removeNoSameBunchPileup{"removeNoSameBunchPileup", false, "Remove no same bunch pileup"};
  Configurable<bool> requireIsGoodZvtxFT0vsPV{"requireIsGoodZvtxFT0vsPV", false, "Require is good Zvtx FT0 vs PV"};
  Configurable<bool> requireIsVertexITSTPC{"requireIsVertexITSTPC", false, "Require is vertex ITS TPC"};
  Configurable<bool> removeNoTimeFrameBorder{"removeNoTimeFrameBorder", false, "Remove no time frame border"};
  Configurable<bool> cfgsel8Require{"cfgsel8Require", true, "sel8 cut require"};
  Configurable<float> cfgZvertex{"cfgZvertex", 10, "Min Z Vertex"};

  // ---- Track selection configurables ------------------------------------
  Configurable<bool> cfgUsePVcontributors{"cfgUsePVcontributors", true, "use tracks that are PV contibutors"};
  Configurable<bool> cfgPassedITSRefit{"cfgPassedITSRefit", true, "Require ITS refit"};
  Configurable<bool> cfgPassedTPCRefit{"cfgPassedTPCRefit", true, "Require TPC refit"};
  Configurable<bool> cfgetaRequire{"cfgetaRequire", true, "eta cut require"};
  Configurable<bool> cfgetaRequireMC{"cfgetaRequireMC", true, "eta cut require for generated particles"};
  Configurable<bool> cfgRapidityRequire{"cfgRapidityRequire", true, "Require Rapidity cut"};
  Configurable<bool> cfgRapidityRequireMC{"cfgRapidityRequireMC", true, "rapidity cut require for generated particles"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.9f, "Eta range for tracks"};
  Configurable<float> cfgCutRapiditymin{"cfgCutRapiditymin", -0.5f, "min Rapidity range"};
  Configurable<float> cfgCutRapiditymax{"cfgCutRapiditymax", 0.5f, "max Rapidity range"};
  Configurable<float> cfgtpcNClsFindable{"cfgtpcNClsFindable", 0.8f, "tpcNClsFindable over crossedRows"};
  Configurable<float> centcut{"centcut", 80.0f, "centrality cut"};
  Configurable<bool> cfgZvertexRequireMC{"cfgZvertexRequireMC", true, "Pos Z cut in MC"};

  // ---- Track quality-cut configurables -----------------------------------
  Configurable<bool> cfgTPCNClsCrossedRowsRequire{"cfgTPCNClsCrossedRowsRequire", true, "Require TPCNClsCrossedRows Cut"};
  Configurable<bool> cfgmaxTPCchi2Require{"cfgmaxTPCchi2Require", true, "Require maxTPCchi2 Cut"};
  Configurable<bool> cfgminTPCchi2Require{"cfgminTPCchi2Require", true, "Require minTPCchi2 Cut"};
  Configurable<bool> cfgminITSnClsRequire{"cfgminITSnClsRequire", false, "Require minITSnCls Cut"};
  Configurable<bool> cfgminITSnClscosRequire{"cfgminITSnClscosRequire", true, "Require minITSnCls / cosh(eta) Cut"};
  Configurable<bool> cfgminReqClusterITSibRequire{"cfgminReqClusterITSibRequire", true, " Require min number of clusters required in ITS inner barrel"};
  Configurable<bool> cfgmaxITSchi2Require{"cfgmaxITSchi2Require", true, "Require maxITSchi2 Cut"};
  Configurable<bool> cfgmaxTPCnSigmaRequire{"cfgmaxTPCnSigmaRequire", true, "Require maxTPCnSigma Cut"};
  Configurable<bool> cfgminGetMeanItsClsSizeRequire{"cfgminGetMeanItsClsSizeRequire", true, "Require minGetMeanItsClsSize Cut"};
  Configurable<bool> cfgmaxGetMeanItsClsSizeRequire{"cfgmaxGetMeanItsClsSizeRequire", true, "Require maxGetMeanItsClsSize Cut"};

  // ---- He4-specific configurables -----------------------------------------
  Configurable<bool> cfghe3massrejreq{"cfghe3massrejreq", true, "Require mass cut on He4 particles"};
  Configurable<float> cfgminmassrejection{"cfgminmassrejection", 6.5, "Min side of He3 particle rejection"};
  Configurable<float> cfgmaxmassrejection{"cfgmaxmassrejection", 9.138, "Max side of He3 particle rejection"};
  Configurable<float> deuteronsigmarejection{"deuteronsigmarejection", 2.0f, "Deuteron TPC nsigma rejection for He4"};

  // ---- MC-only configurables -----------------------------------------------
  Configurable<bool> cfgmccorrectionhe4Require{"cfgmccorrectionhe4Require", true, "MC correction for pp he4 particle"};

  // ---- DCA cut configurables -------------------------------------------------
  Configurable<bool> cfgUseDCAxyCorrection{"cfgUseDCAxyCorrection", true, "Use pT-dependent DCAxy cut"};
  Configurable<bool> cfgUseDCAzCorrection{"cfgUseDCAzCorrection", true, "Use pT-dependent DCAz cut"};
  Configurable<bool> cfgMultiplyTransportPt{"cfgMultiplyTransportPt", true, "In the DCA templates (processDCA), double the reconstructed pT for transport (material) secondaries; if false, use pT as reconstructed"};

  // ---- Misc configurables -----------------------------------------------------
  Configurable<int> cfgDebug{"cfgDebug", 1, "debug level"};
  Configurable<bool> cfgRigidityCorrection{"cfgRigidityCorrection", false, "apply rigidity correction"};
  Configurable<bool> cfgRequirebetaplot{"cfgRequirebetaplot", true, "Require beta plot"};
  Configurable<bool> cfgmass2{"cfgmass2", true, "Fill mass square difference"};
  Configurable<bool> cfgFillmass{"cfgFillmass", false, "Fill mass histograms"};
  Configurable<bool> cfgFillmassnsigma{"cfgFillmassnsigma", true, "Fill mass vs nsigma histograms"};

  // ---- Table-valued configurables (species x parameter) ------------------
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {kBetheBlochDefault[0].data(), NParticles, NBetheParams, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgTrackPIDsettings{"cfgTrackPIDsettings", {kTrackPIDSettings[0].data(), NParticles, NTrkSettings, particleNames, trackPIDsettingsNames}, "track selection and PID criteria"};
  Configurable<LabeledArray<double>> cfgTrackPIDsettings2{"cfgTrackPIDsettings2", {kTrackPIDSettings2[0].data(), NParticles, NTrkSettings2, particleNames, trackPIDsettingsNames2}, "track selection and PID criteria"};
  Configurable<LabeledArray<double>> cfgktrackcorrection{"cfgktrackcorrection", {ktrackcorrection[0].data(), NFittingParticle, NFittingParameters, correctedparticleNames, trackcorrectionNames}, "fitting parameters"};
  Configurable<LabeledArray<double>> cfgDCAcorrection{"cfgDCAcorrection", {kDCAcorrection[0].data(), NDCATypes, NDCAParameters, dcaTypeNames, dcaNames}, "DCA parameters: DCAxy: p0*exp(p1*pt)+p2, DCAz: p0*exp(-p1*pt)+p2"};

  o2::track::TrackParametrizationWithError<float> mTrackParCov;

  // ---- Binning configuration -----------------------------------------------
  ConfigurableAxis axisMagField{"axisMagField", {10, -10., 10.}, "magnetic field"};
  ConfigurableAxis axisNev{"axisNev", {10, 0., 10.}, "Number of events"};
  ConfigurableAxis axisRigidity{"axisRigidity", {4000, -10., 10.}, "#it{p}^{TPC}/#it{z}"};
  ConfigurableAxis axisdEdx{"axisdEdx", {4000, 0, 4000}, "d#it{E}/d#it{x}"};
  ConfigurableAxis axisCent{"axisCent", {100, 0, 100}, "centrality"};
  ConfigurableAxis axisOccupancy{"axisOccupancy", {5000, 0, 50000}, "axis for Occupancy of event"};
  ConfigurableAxis axisVtxZ{"axisVtxZ", {120, -20, 20}, "z"};
  ConfigurableAxis ptAxis{"ptAxis", {200, 0, 10}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axiseta{"axiseta", {100, -1, 1}, "eta"};
  ConfigurableAxis axisrapidity{"axisrapidity", {100, -2, 2}, "rapidity"};
  ConfigurableAxis axismass{"axismass", {1200, 0, 12}, "mass"};
  ConfigurableAxis axismassnsigma{"axismassnsigma", {1200, 0, 12}, "nsigma mass"};
  ConfigurableAxis nsigmaAxis{"nsigmaAxis", {160, -10, 10}, "n#sigma_{#pi^{+}}"};
  ConfigurableAxis axisDCA{"axisDCA", {400, -10., 10.}, "DCA axis"};
  ConfigurableAxis particleAntiAxis{"particleAntiAxis", {2, -0.5, 1.5}, "Particle/Anti-particle"}; // 0 = particle, 1 = anti-particle
  ConfigurableAxis decayTypeAxis{"decayTypeAxis", {3, -0.5, 2.5}, "Decay type"};                   // 0 = primary, 1 = from decay, 2 = material

  // ---- CCDB -------------------------------------------------------------------
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
  bool collHasCandidate = false;
  bool collPassedEvSel = false;
  int mRunNumber = 0;
  int occupancy = 0;
  float dBz = 0.0f;
  TRandom3 rand;

  // Species indices into primaryParticles / particlePdgCodes (order fixed by particleNames above).
  const float proton = 1;
  const float he3 = 4;
  const float he4 = 5;

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

    for (int i = 0; i < NParticles; i++) {
      primaryParticles.push_back(PrimParticles(particleNames.at(i), particlePdgCodes.at(i), particleMasses.at(i), particleCharge.at(i), cfgBetheBlochParams));
    }

    if (doprocessData) {
      histos.add("histNev", "histNev", kTH1F, {axisNev});
      histos.add("histVtxZ", "histVtxZ", kTH1F, {axisVtxZ});
      histos.add("histCentFT0C", "histCentFT0C", kTH1F, {axisCent});
      histos.add("histCentFT0M", "histCentFT0M", kTH1F, {axisCent});
      histos.add("histCentFTOC_cut", "histCentFTOC_cut", kTH1F, {axisCent});
      histos.add<THnSparse>("hSpectra", " ", HistType::kTHnSparseF, {ptAxis, nsigmaAxis, {5, -2.5, 2.5}, axisCent});
      histos.add<THnSparse>("hSpectratof", " ", HistType::kTHnSparseF, {ptAxis, nsigmaAxis, {5, -2.5, 2.5}, axisCent});
      histos.add("DCAxy_vs_pT_data", "DCA_{xy} vs p_{T} for He3 (Data);p_{T} (GeV/c);DCA_{xy} (cm)",
                 {HistType::kTH3F, {ptAxis, axisDCA, axisCent}});
      histos.add("DCAxy_vs_pT_anti_data", "DCA_{xy} vs p_{T} for anti-He3 (Data);p_{T} (GeV/c);DCA_{xy} (cm)",
                 {HistType::kTH3F, {ptAxis, axisDCA, axisCent}});

      hmass.resize(2 * NParticles + 2);
      hmassnsigma.resize(2 * NParticles + 2);

      for (int i = 0; i < NParticles; i++) {
        TString histName = primaryParticles[i].name;
        if (cfgFillmass) {
          hmass[2 * i] = histos.add<TH3>(Form("histmass_pt/histmass_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); mass^{2};  centrality(%)", HistType::kTH3F, {ptAxis, axismass, axisCent});
          hmass[2 * i + 1] = histos.add<TH3>(Form("histmass_ptanti/histmass_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); mass^{2}; centrality(%)", HistType::kTH3F, {ptAxis, axismass, axisCent});
        }
      }
      for (int i = 0; i < NParticles; i++) {
        TString histName = primaryParticles[i].name;
        if (cfgFillmassnsigma) {
          hmassnsigma[2 * i] = histos.add<TH3>(Form("histmass_nsigma/histmass_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); mass^{2}/z^{2}; Centrality(%)", HistType::kTH3F, {ptAxis, axismassnsigma, axisCent});
          hmassnsigma[2 * i + 1] = histos.add<TH3>(Form("histmass_nsigmaanti/histmass_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); mass^{2}/z^{2}; centrality(%)", HistType::kTH3F, {ptAxis, axismassnsigma, axisCent});
        }
      }

      histos.add("histeta", "histeta", kTH1F, {axiseta});
      histos.add("histEvents", "histEvents", kTH2F, {axisCent, axisOccupancy});
    }

    histos.add("dcaZ", "dcaZ", kTH2F, {ptAxis, axisDCA});
    histos.add("dcaXY", "dcaXY", kTH2F, {ptAxis, axisDCA});
    histos.add("Tofsignal", "Tofsignal", kTH2F, {axisRigidity, {4000, 0.2, 1.2, "#beta"}});
    histos.add("Tpcsignal", "Tpcsignal", kTH2F, {axisRigidity, axisdEdx});

    if (doprocessMC) {
      histomc.add<THnSparse>("hSpectramc", " ", HistType::kTHnSparseF, {particleAntiAxis, axisCent, ptAxis, ptAxis});

      // Efficiency x Acceptance
      histomc.add("hDenomEffAcc", "Denominator for Efficiency x Acceptance",
                  {HistType::kTHnSparseF, {ptAxis, axisCent, particleAntiAxis, decayTypeAxis}});
      histomc.add("hNumerEffAcc", "Numerator for Efficiency x Acceptance",
                  {HistType::kTHnSparseF, {ptAxis, axisCent, particleAntiAxis, decayTypeAxis}});

      // Signal loss correction
      histomc.add("hSignalLossDenom", "Signal Loss Denominator", kTH2F, {axisCent, ptAxis});
      histomc.add("hSignalLossNumer", "Signal Loss Numerator", kTH2F, {axisCent, ptAxis});
      histomc.add("haSignalLossDenom", "anti particle Signal Loss Denominator", kTH2F, {axisCent, ptAxis});
      histomc.add("haSignalLossNumer", "antiparticle Signal Loss Numerator", kTH2F, {axisCent, ptAxis});

      // Event loss correction
      histomc.add("hEventLossDenom", "Event loss denominator", kTH1F, {axisCent});
      histomc.add("hEventLossNumer", "Event loss numerator", kTH1F, {axisCent});

      histomc.add("histVtxZgen", "histVtxZgen", kTH1F, {axisVtxZ});
      histomc.add("histVtxZReco", "histVtxZReco", kTH1F, {axisVtxZ});

      histomc.add("histDeltaPtVsPtGenanti", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10}, {1000, -0.5, 0.5, "p_{T}(reco) - p_{T}(gen);p_{T}(reco)"}});
      histomc.add("histDeltaPtVsPtGen", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10}, {1000, -0.5, 0.5, "p_{T}(reco) - p_{T}(gen);p_{T}(reco)"}});
      histomc.add("histPIDtrack", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10, "p_{T}(reco)"}, {9, -0.5, 8.5, "p_{T}(reco) - p_{T}(gen)"}});
      histomc.add("histPIDtrackanti", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10, "p_{T}(reco)"}, {9, -0.5, 8.5, "p_{T}(reco) - p_{T}(gen)"}});

      // Mass^2/z^2 3D plots for MC (pT, mass^2/z^2, centrality) - separate for particles and anti-particles
      histomc.add("hMassVsPtMC", "mass^{2}/z^{2} vs p_{T} (MC Particles);p_{T} (GeV/c);mass^{2}/z^{2};Centrality",
                  {HistType::kTH3F, {ptAxis, axismassnsigma, axisCent}});
      histomc.add("hMassVsPtAntiMC", "mass^{2}/z^{2} vs p_{T} (MC Anti-particles);p_{T} (GeV/c);mass^{2}/z^{2};Centrality",
                  {HistType::kTH3F, {ptAxis, axismassnsigma, axisCent}});
    }

    if (doprocessDCA) {
      histomc.add("DCAxy_vs_pT_primary", "DCA_{xy} vs p_{T} (primary);p_{T} (GeV/c);DCA_{xy} (cm)",
                  {HistType::kTH3F, {ptAxis, axisDCA, axisCent}});
      histomc.add("DCAxy_vs_pT_weakdecay", "DCA_{xy} vs p_{T}  (Weak Decay);p_{T} (GeV/c);DCA_{xy} (cm)",
                  {HistType::kTH3F, {ptAxis, axisDCA, axisCent}});
      histomc.add("DCAxy_vs_pT_transport", "DCA_{xy} vs p_{T} (Transport);p_{T} (GeV/c);DCA_{xy} (cm)",
                  {HistType::kTH3F, {ptAxis, axisDCA, axisCent}});
      histomc.add("DCAxy_vs_pT_anti_primary", "DCA_{xy} vs p_{T} for anti (primary);p_{T} (GeV/c);DCA_{xy} (cm)",
                  {HistType::kTH3F, {ptAxis, axisDCA, axisCent}});
      histomc.add("DCAxy_vs_pT_anti_weakdecay", "DCA_{xy} vs p_{T} for anti (Weak Decay);p_{T} (GeV/c);DCA_{xy} (cm)",
                  {HistType::kTH3F, {ptAxis, axisDCA, axisCent}});
      histomc.add("DCAxy_vs_pT_anti_transport", "DCA_{xy} vs p_{T} for anti (Transport);p_{T} (GeV/c);DCA_{xy} (cm)",
                  {HistType::kTH3F, {ptAxis, axisDCA, axisCent}});
    }
  }

  //====================================================================================
  // processData: reconstructed data, fills TPC(+TOF) PID spectra and DCA/mass QA
  //====================================================================================
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
      if (removeNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup))
        continue;
      histos.fill(HIST("histNev"), 2.5);

      if (requireIsGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
        continue;
      histos.fill(HIST("histNev"), 3.5);
      if (requireIsVertexITSTPC && !collision.selection_bit(aod::evsel::kIsVertexITSTPC))
        continue;

      histos.fill(HIST("histNev"), 4.5);
      if (removeNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
        continue;
      histos.fill(HIST("histNev"), 5.5);
      histos.fill(HIST("histEvents"), collision.centFT0C(), occupancy);
      histos.fill(HIST("histVtxZ"), collision.posZ());
      histos.fill(HIST("histCentFT0C"), collision.centFT0C());
      histos.fill(HIST("histCentFT0M"), collision.centFT0M());
      if (collision.centFT0C() > centcut)
        continue;
      histos.fill(HIST("histNev"), 6.5);
      histos.fill(HIST("histCentFTOC_cut"), collision.centFT0C());

      auto tracksInColl = tracks.sliceBy(tracksPerCollision, collision.globalIndex());
      for (const auto& track : tracksInColl) {
        if (!passesBasicRecoTrackCuts(track))
          continue;

        for (size_t i = 0; i < primaryParticles.size(); i++) {
          if (cfgTrackPIDsettings2->get(i, "fillsparsh") != 1)
            continue;

          float ptMomn = computeSpeciesPt(track, i);
          float tpcNsigma = getTPCnSigma(track, primaryParticles.at(i));
          if ((std::abs(tpcNsigma) > cfgTrackPIDsettings->get(i, "maxTPCnSigma")) && cfgmaxTPCnSigmaRequire)
            continue;

          int sign = (track.sign() > 0) ? 1 : ((track.sign() < 0) ? -1 : 0);
          float rapidity = getRapidity(track, i);
          if ((rapidity < cfgCutRapiditymin || rapidity > cfgCutRapiditymax) && cfgRapidityRequire)
            continue;

          if (!passesTrackQualityCuts(track, i))
            continue;
          if (!passesDCACuts(track, i, ptMomn))
            continue;

          float itsSigma = getITSnSigma(track, primaryParticles.at(i));
          if (itsSigma < cfgTrackPIDsettings2->get(i, "minITSnsigma") && cfgTrackPIDsettings2->get(i, "useITSnsigma") < 1)
            continue;
          if (itsSigma > cfgTrackPIDsettings2->get(i, "maxITSnsigma") && cfgTrackPIDsettings2->get(i, "useITSnsigma") < 1)
            continue;

          histos.fill(HIST("Tpcsignal"), getRigidity(track) * track.sign(), track.tpcSignal());
          histos.fill(HIST("dcaXY"), ptMomn, track.dcaXY());
          histos.fill(HIST("dcaZ"), ptMomn, track.dcaZ());

          float tofnsigma = getTOFnSigma(track, primaryParticles.at(i));
          if (track.sign() > 0) {
            histos.fill(HIST("DCAxy_vs_pT_data"), ptMomn, track.dcaXY(), collision.centFT0C());
          } else if (track.sign() < 0) {
            histos.fill(HIST("DCAxy_vs_pT_anti_data"), ptMomn, track.dcaXY(), collision.centFT0C());
          }

          // He4 selection uses the deuteron TPC nsigma to reject deuterons.
          float tpcNsigmaDe = track.tpcNSigmaDe();

          if (i != he4) {
            histos.fill(HIST("hSpectra"), ptMomn, tpcNsigma, sign, collision.centFT0C());
            if (track.hasTOF()) {
              histos.fill(HIST("hSpectratof"), ptMomn, tofnsigma, sign, collision.centFT0C());
            }
          } else if (!track.hasTOF()) {
            if (std::abs(tpcNsigmaDe) > deuteronsigmarejection) {
              histos.fill(HIST("hSpectra"), ptMomn, tpcNsigma, sign, collision.centFT0C());
            }
          } else {
            // Has TOF - apply He4 mass rejection.
            float beta = o2::pid::tof::Beta::GetBeta(track);
            const float eps = 1e-6f;
            if (beta < eps || beta > 1.0f - eps)
              continue;

            float charge = 2.f; // he4 has charge 2
            float p = getRigidity(track);
            float massTOF = p * charge * std::sqrt(1.f / (beta * beta) - 1.f);
            if (cfghe3massrejreq && (massTOF * massTOF > cfgminmassrejection && massTOF * massTOF < cfgmaxmassrejection)) {
              continue; // reject He4 mass window
            }
            if (std::abs(tpcNsigmaDe) > deuteronsigmarejection) {
              histos.fill(HIST("hSpectra"), ptMomn, tpcNsigma, sign, collision.centFT0C());
            }
          }

          if (std::abs(tpcNsigma) > cfgTrackPIDsettings2->get(i, "maxTPCnsigmaTOF")) {
            if (i == he4) {
              if (std::abs(tpcNsigmaDe) > deuteronsigmarejection) {
                fillhmassnsigma(track, i, collision.centFT0C());
              }
            } else {
              fillhmassnsigma(track, i, collision.centFT0C());
            }
          }

          if ((std::abs(tpcNsigma) > cfgTrackPIDsettings2->get(i, "maxTPCnsigmaTOF")) && cfgTrackPIDsettings2->get(i, "useTPCnsigmaTOF") < 1)
            continue;

          if (i == he4) {
            if (std::abs(tpcNsigmaDe) > deuteronsigmarejection) {
              fillhmass(track, i, collision.centFT0C());
            }
          } else if (i == proton && (std::abs(tofnsigma) < cfgTrackPIDsettings2->get(i, "maxTOFnsigma"))) {
            fillhmassnsigma(track, i, collision.centFT0C());
          } else {
            fillhmass(track, i, collision.centFT0C());
          }

          if (cfgRequirebetaplot) {
            histos.fill(HIST("Tofsignal"), getRigidity(track) * track.sign(), o2::pid::tof::Beta::GetBeta(track));
          }
        } // species loop

        histos.fill(HIST("histeta"), track.eta());
      } // track loop
    } // collision loop
  }
  PROCESS_SWITCH(NucleitpcPbPb, processData, "data analysis", false);

  //====================================================================================
  // processMC: MC reco+gen, for efficiency x acceptance, signal-loss and event-loss
  //====================================================================================
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

    // ---- Pass 1: store centrality and event-selection flags per MC collision ----
    for (auto const& collision : collisions) {
      int mcCollIdx = collision.mcCollisionId();
      if (mcCollIdx < 0 || mcCollIdx >= static_cast<int>(mcCollisions.size())) {
        continue;
      }

      mcCollInfos[mcCollIdx].centrality = collision.centFT0C(); // stored without cuts

      if (!collision.sel8() && cfgsel8Require)
        continue;
      if (collision.centFT0C() > centcut)
        continue;
      if (removeNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup))
        continue;
      if (requireIsGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
        continue;
      if (requireIsVertexITSTPC && !collision.selection_bit(aod::evsel::kIsVertexITSTPC))
        continue;
      if (removeNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
        continue;

      mcCollInfos[mcCollIdx].passedEvSel = true;

      if (std::abs(collision.posZ()) > cfgZvertex && cfgZvertexRequireMC)
        continue;
      mcCollInfos[mcCollIdx].passedEvSelVtZ = true;
    }

    // ---- Event-loss and signal-loss denominators/numerators ----
    for (size_t i = 0; i < mcCollInfos.size(); i++) {
      if (mcCollInfos[i].centrality < 0)
        continue; // no matching reconstructed collision found

      histomc.fill(HIST("hEventLossDenom"), mcCollInfos[i].centrality);
      if (mcCollInfos[i].passedEvSel) {
        histomc.fill(HIST("hEventLossNumer"), mcCollInfos[i].centrality);
      }

      for (auto const& mcParticle : particlesMC) {
        if (mcParticle.mcCollisionId() != static_cast<int>(i))
          continue;
        if (!mcParticle.isPhysicalPrimary())
          continue;

        int particleType = findParticleTypeIndex(mcParticle.pdgCode());
        if (!speciesEnabled(particleType))
          continue;

        if (mcParticle.pdgCode() == particlePdgCodes.at(particleType)) {
          histomc.fill(HIST("hSignalLossDenom"), mcCollInfos[i].centrality, mcParticle.pt());
        } else if (mcParticle.pdgCode() == -particlePdgCodes.at(particleType)) {
          histomc.fill(HIST("haSignalLossDenom"), mcCollInfos[i].centrality, mcParticle.pt());
        }

        if (mcCollInfos[i].passedEvSel) {
          if (mcParticle.pdgCode() == particlePdgCodes.at(particleType)) {
            histomc.fill(HIST("hSignalLossNumer"), mcCollInfos[i].centrality, mcParticle.pt());
          } else if (mcParticle.pdgCode() == -particlePdgCodes.at(particleType)) {
            histomc.fill(HIST("haSignalLossNumer"), mcCollInfos[i].centrality, mcParticle.pt());
          }
        }
      }
    }

    // ---- Generated particles: efficiency x acceptance denominator, and reco pass ----
    for (auto const& mcCollision : mcCollisions) {
      size_t idx = mcCollision.globalIndex();
      if (idx >= mcCollInfos.size())
        continue;

      float centrality = mcCollInfos[idx].centrality;
      bool passedEvSelVtZ = mcCollInfos[idx].passedEvSelVtZ;

      histomc.fill(HIST("histVtxZgen"), mcCollision.posZ());

      for (auto const& mcParticle : particlesMC) {
        if (mcParticle.mcCollisionId() != mcCollision.globalIndex())
          continue;

        int pdgCode = mcParticle.pdgCode();
        int particleType = findParticleTypeIndex(pdgCode);
        if (!speciesEnabled(particleType))
          continue;

        if (std::abs(mcParticle.eta()) > cfgCutEta && cfgetaRequireMC)
          continue;
        float rapidity = mcParticle.y();
        if ((rapidity < cfgCutRapiditymin || rapidity > cfgCutRapiditymax) && cfgRapidityRequireMC)
          continue;

        int particleAnti = (pdgCode > 0) ? 0 : 1;
        int decayType = classifyDecayType(mcParticle);
        if (decayType == DecayTypeTransport)
          continue; // transported secondaries are not part of the denominator

        if (passedEvSelVtZ) {
          histomc.fill(HIST("hDenomEffAcc"), mcParticle.pt(), centrality, particleAnti, decayType);
        }
      }

      if (!passedEvSelVtZ)
        continue;

      // ---- Find the reconstructed collision matching this MC collision ----
      for (auto const& collision : collisions) {
        if (collision.mcCollisionId() != static_cast<int>(idx))
          continue;

        auto bc = collision.bc_as<aod::BCsWithTimestamps>();
        initCCDB(bc);
        histomc.fill(HIST("histVtxZReco"), collision.posZ());

        auto tracksInColl = tracks.sliceBy(tracksPerCollision, collision.globalIndex());
        for (auto const& track : tracksInColl) {
          if (!track.has_mcParticle())
            continue; // skip un-matched reco tracks

          auto const& matchedMCParticle = track.mcParticle_as<aod::McParticles>();
          if (matchedMCParticle.mcCollisionId() != mcCollision.globalIndex())
            continue;

          int pdg = matchedMCParticle.pdgCode();
          int decayType = classifyDecayType(matchedMCParticle);

          if (!passesBasicRecoTrackCuts(track))
            continue;

          for (size_t i = 0; i < primaryParticles.size(); i++) {
            if (std::abs(pdg) != std::abs(particlePdgCodes.at(i)))
              continue;
            if (cfgTrackPIDsettings2->get(i, "fillsparsh") != 1)
              continue;

            float ptReco = computeSpeciesPt(track, i);
            ptReco = applyHe34PtCorrection(track, pdg, i, ptReco);

            int particleAnti = (pdg > 0) ? 0 : 1;

            float rapidity = getRapidity(track, i);
            if ((rapidity < cfgCutRapiditymin || rapidity > cfgCutRapiditymax) && cfgRapidityRequire)
              continue;

            if (!passesTrackQualityCuts(track, i))
              continue;
            if (!passesDCACuts(track, i, ptReco))
              continue;

            float tpcNsigma = getTPCnSigma(track, primaryParticles.at(i));
            bool passesNsigma = !((std::abs(tpcNsigma) > cfgTrackPIDsettings->get(i, "maxTPCnSigma")) && cfgmaxTPCnSigmaRequire);

            // Single fill per track: axis 2 (pTReco) always gets the reco pT,
            // with all cuts applied except TPC nsigma and hasTOF -- this is
            // the denominator population for the TOF matching efficiency.
            // Axis 3 (pTOF) gets the same pT value only for tracks that both
            // pass the TPC nsigma cut and have TOF -- the numerator
            // population -- otherwise it gets a sentinel (-1, landing in the
            // pTOF axis underflow bin since ptAxis starts at 0) so it's
            // excluded from any real-pT projection of that axis. Filling
            // once (not twice) avoids double-counting TOF-matched tracks in
            // the pTReco axis relative to non-matched tracks.
            bool isTOFmatched = passesNsigma && track.hasTOF();
            float ptTOF = isTOFmatched ? ptReco : -1.0f;
            histomc.fill(HIST("hSpectramc"), particleAnti, collision.centFT0C(), ptReco, ptTOF);

            if (!passesNsigma)
              continue;

            // Efficiency x Acceptance numerator (only after the TPC nsigma cut).
            histomc.fill(HIST("hNumerEffAcc"), ptReco, collision.centFT0C(), particleAnti, decayType);

            if (decayType == DecayTypePrimary) {
              histos.fill(HIST("dcaXY"), ptReco, track.dcaXY());
              histos.fill(HIST("dcaZ"), ptReco, track.dcaZ());
              histos.fill(HIST("Tpcsignal"), getRigidity(track) * track.sign(), track.tpcSignal());
            }

            float ptGen = matchedMCParticle.pt();
            float deltaPt = ptReco - ptGen;
            if (pdg == -particlePdgCodes.at(i) && decayType == DecayTypePrimary) { // anti-particle
              histomc.fill(HIST("histDeltaPtVsPtGenanti"), ptReco, deltaPt);
              histomc.fill(HIST("histPIDtrackanti"), ptReco, track.pidForTracking());
            }
            if (pdg == particlePdgCodes.at(i) && decayType == DecayTypePrimary) { // particle
              histomc.fill(HIST("histDeltaPtVsPtGen"), ptReco, deltaPt);
              histomc.fill(HIST("histPIDtrack"), ptReco, track.pidForTracking());
            }

            // Mass^2/z^2 for MC, split by particle/anti-particle.
            if (track.hasTOF()) {
              float beta = o2::pid::tof::Beta::GetBeta(track);
              const float eps = 1e-6f;
              if (beta >= eps && beta <= 1.0f - eps) {
                float charge = (i == he3 || i == he4) ? 2.f : 1.f;
                float p = getRigidity(track);
                float massTOF = p * charge * std::sqrt(1.f / (beta * beta) - 1.f);
                float massSquareOverChargeSquare = (massTOF * massTOF) / (charge * charge);

                bool skipHe4 = (std::abs(pdg) == particlePdgCodes.at(5)) && cfghe3massrejreq &&
                               (massTOF * massTOF > cfgminmassrejection && massTOF * massTOF < cfgmaxmassrejection);

                if (!skipHe4 && decayType == DecayTypePrimary) {
                  if (pdg > 0) {
                    histomc.fill(HIST("hMassVsPtMC"), ptReco, massSquareOverChargeSquare, collision.centFT0C());
                  } else {
                    histomc.fill(HIST("hMassVsPtAntiMC"), ptReco, massSquareOverChargeSquare, collision.centFT0C());
                  }
                }
              }
            }
          } // species loop
        } // track loop
        break; // matching reconstructed collision found and processed
      } // collision loop
    } // mcCollision loop
  }
  PROCESS_SWITCH(NucleitpcPbPb, processMC, "MC reco+gen analysis with efficiency corrections", false);

  //====================================================================================
  // processDCA: MC reco, builds DCAxy templates (primary / weak-decay / transport)
  // used for the secondary-contamination correction. Note: no DCA cut is applied
  // here, since the DCA distribution itself is what is being measured.
  //====================================================================================
  void processDCA(CollisionsFullMC const& collisions,
                  aod::McCollisions const& mcCollisions,
                  aod::McParticles const& particlesMC,
                  soa::Join<TracksFull, aod::McTrackLabels> const& tracks,
                  aod::BCsWithTimestamps const&)
  {
    (void)particlesMC;
    mcCollInfos.clear();
    mcCollInfos.resize(mcCollisions.size());

    for (auto const& collision : collisions) {
      int mcCollIdx = collision.mcCollisionId();
      if (mcCollIdx < 0 || mcCollIdx >= static_cast<int>(mcCollisions.size())) {
        continue;
      }

      mcCollInfos[mcCollIdx].centrality = collision.centFT0C();

      if (!collision.sel8() && cfgsel8Require)
        continue;
      if (collision.centFT0C() > centcut)
        continue;
      if (removeNoSameBunchPileup && !collision.selection_bit(aod::evsel::kNoSameBunchPileup))
        continue;
      if (requireIsGoodZvtxFT0vsPV && !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))
        continue;
      if (requireIsVertexITSTPC && !collision.selection_bit(aod::evsel::kIsVertexITSTPC))
        continue;
      if (removeNoTimeFrameBorder && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder))
        continue;

      mcCollInfos[mcCollIdx].passedEvSel = true;

      if (std::abs(collision.posZ()) > cfgZvertex && cfgZvertexRequireMC)
        continue;
      mcCollInfos[mcCollIdx].passedEvSelVtZ = true;
    }

    for (auto const& mcCollision : mcCollisions) {
      size_t idx = mcCollision.globalIndex();
      if (idx >= mcCollInfos.size())
        continue;
      if (!mcCollInfos[idx].passedEvSelVtZ)
        continue;

      for (auto const& collision : collisions) {
        if (collision.mcCollisionId() != static_cast<int>(idx))
          continue;

        auto bc = collision.bc_as<aod::BCsWithTimestamps>();
        initCCDB(bc);
        auto tracksInColl = tracks.sliceBy(tracksPerCollision, collision.globalIndex());

        for (auto const& track : tracksInColl) {
          if (!track.has_mcParticle())
            continue;

          auto const& matchedMCParticle = track.mcParticle_as<aod::McParticles>();
          if (matchedMCParticle.mcCollisionId() != mcCollision.globalIndex())
            continue;

          int pdg = matchedMCParticle.pdgCode();
          if (!passesBasicRecoTrackCuts(track))
            continue;

          for (size_t i = 0; i < primaryParticles.size(); i++) {
            if (std::abs(pdg) != std::abs(particlePdgCodes.at(i)))
              continue;
            if (cfgTrackPIDsettings2->get(i, "fillsparsh") != 1)
              continue;

            float ptReco = computeSpeciesPt(track, i);
            ptReco = applyHe34PtCorrection(track, pdg, i, ptReco);

            float rapidity = getRapidity(track, i);
            if ((rapidity < cfgCutRapiditymin || rapidity > cfgCutRapiditymax) && cfgRapidityRequire)
              continue;

            if (!passesTrackQualityCuts(track, i))
              continue;

            float tpcNsigma = getTPCnSigma(track, primaryParticles.at(i));
            if ((std::abs(tpcNsigma) > cfgTrackPIDsettings->get(i, "maxTPCnSigma")) && cfgmaxTPCnSigmaRequire)
              continue;

            // Decay-origin classification specific to the DCA-template building:
            // physical primaries and weak decays (via getProcess()==kPDecay) are
            // distinguished from transport/material secondaries.
            int decayType = -999;
            if (matchedMCParticle.isPhysicalPrimary()) {
              decayType = DecayTypePrimary;
              if (matchedMCParticle.has_mothers()) {
                for (const auto& motherparticle : matchedMCParticle.mothers_as<aod::McParticles>()) {
                  if (std::find(hfMothCodes.begin(), hfMothCodes.end(), std::abs(motherparticle.pdgCode())) != hfMothCodes.end()) {
                    decayType = DecayTypeWeak;
                    break;
                  }
                }
              }
            } else if (matchedMCParticle.getProcess() == TMCProcess::kPDecay) {
              if (!matchedMCParticle.has_mothers()) {
                continue; // secondaries from weak decay without mothers are skipped
              }
              decayType = DecayTypeWeak;
            } else {
              decayType = DecayTypeTransport;
            }

            // Transport (material) secondaries optionally get their reconstructed
            // pT doubled for the DCA template, per cfgMultiplyTransportPt; primary
            // and weak-decay secondaries always use pT as reconstructed.
            float ptTransport = cfgMultiplyTransportPt ? 2.f * ptReco : ptReco;

            if (pdg == particlePdgCodes.at(i)) { // particle
              if (decayType == DecayTypePrimary) {
                histomc.fill(HIST("DCAxy_vs_pT_primary"), ptReco, track.dcaXY(), collision.centFT0C());
              } else if (decayType == DecayTypeWeak) {
                histomc.fill(HIST("DCAxy_vs_pT_weakdecay"), ptReco, track.dcaXY(), collision.centFT0C());
              } else if (decayType == DecayTypeTransport) {
                histomc.fill(HIST("DCAxy_vs_pT_transport"), ptTransport, track.dcaXY(), collision.centFT0C());
              }
            } else if (pdg == -particlePdgCodes.at(i)) { // anti-particle
              if (decayType == DecayTypePrimary) {
                histomc.fill(HIST("DCAxy_vs_pT_anti_primary"), ptReco, track.dcaXY(), collision.centFT0C());
              } else if (decayType == DecayTypeWeak) {
                histomc.fill(HIST("DCAxy_vs_pT_anti_weakdecay"), ptReco, track.dcaXY(), collision.centFT0C());
              } else if (decayType == DecayTypeTransport) {
                histomc.fill(HIST("DCAxy_vs_pT_anti_transport"), ptTransport, track.dcaXY(), collision.centFT0C());
              }
            }
          } // species loop
        } // track loop
        break; // matching reconstructed collision found and processed
      } // collision loop
    } // mcCollision loop
  }
  PROCESS_SWITCH(NucleitpcPbPb, processDCA, "MC DCA analysis For secondary correction", false);

  //====================================================================================
  // Shared per-track helper cuts (extracted from the near-identical blocks that used
  // to be duplicated across processData / processMC / processDCA).
  //====================================================================================

  /// Basic reconstructed-track quality gate: PV contributor, has TPC, ITS/TPC refit,
  /// eta acceptance. Identical requirement in all three process functions.
  template <typename T>
  bool passesBasicRecoTrackCuts(T const& track)
  {
    if (!track.isPVContributor() && cfgUsePVcontributors)
      return false;
    if (!track.hasTPC())
      return false;
    if (!track.passedITSRefit() && cfgPassedITSRefit)
      return false;
    if (!track.passedTPCRefit() && cfgPassedTPCRefit)
      return false;
    if (std::abs(track.eta()) > cfgCutEta && cfgetaRequire)
      return false;
    return true;
  }

  /// Species-dependent TPC/ITS cluster-quality cuts (no DCA cut here -- processDCA
  /// deliberately does not cut on DCA since it is measuring the DCA distribution).
  template <typename T>
  bool passesTrackQualityCuts(T const& track, size_t i)
  {
    if (track.tpcNClsFound() < cfgTrackPIDsettings->get(i, "minTPCnCls"))
      return false;
    if (((track.tpcNClsCrossedRows() < cfgTrackPIDsettings->get(i, "minTPCnClsCrossedRows")) ||
         track.tpcNClsCrossedRows() < cfgtpcNClsFindable * track.tpcNClsFindable()) &&
        cfgTPCNClsCrossedRowsRequire)
      return false;
    if (track.tpcChi2NCl() > cfgTrackPIDsettings->get(i, "maxTPCchi2") && cfgmaxTPCchi2Require)
      return false;
    if (track.tpcChi2NCl() < cfgTrackPIDsettings->get(i, "minTPCchi2") && cfgminTPCchi2Require)
      return false;
    if (track.itsNCls() < cfgTrackPIDsettings->get(i, "minITSnCls") && cfgminITSnClsRequire)
      return false;
    double cosheta = std::cosh(track.eta());
    if ((track.itsNCls() / cosheta) < cfgTrackPIDsettings->get(i, "minITSnClscos") && cfgminITSnClscosRequire)
      return false;
    if ((track.itsNClsInnerBarrel() < cfgTrackPIDsettings->get(i, "minReqClusterITSib")) && cfgminReqClusterITSibRequire)
      return false;
    if (track.itsChi2NCl() > cfgTrackPIDsettings->get(i, "maxITSchi2") && cfgmaxITSchi2Require)
      return false;
    if (getMeanItsClsSize(track) < cfgTrackPIDsettings->get(i, "minITSclsSize") && cfgminGetMeanItsClsSizeRequire)
      return false;
    return true;
  }

  /// pT-dependent DCAxy/DCAz cuts: sigma(pt) = p0*exp(+/-p1*pt) + p2, scaled by the
  /// per-species sigma factor. Used by processData and processMC (not processDCA).
  template <typename T>
  bool passesDCACuts(T const& track, size_t i, float pt)
  {
    double dcaXYp0 = cfgDCAcorrection->get(0U, "p0");
    double dcaXYp1 = cfgDCAcorrection->get(0U, "p1");
    double dcaXYp2 = cfgDCAcorrection->get(0U, "p2");
    double dcaZp0 = cfgDCAcorrection->get(1U, "p0");
    double dcaZp1 = cfgDCAcorrection->get(1U, "p1");
    double dcaZp2 = cfgDCAcorrection->get(1U, "p2");

    bool insideDCAxy;
    if (cfgUseDCAxyCorrection) {
      double sigma = cfgTrackPIDsettings->get(i, "maxDcaXY") * (dcaXYp0 * std::exp(dcaXYp1 * pt) + dcaXYp2);
      insideDCAxy = std::abs(track.dcaXY()) <= sigma;
    } else {
      insideDCAxy = std::abs(track.dcaXY()) <= cfgTrackPIDsettings->get(i, "maxDcaXY");
    }

    bool insideDCAz;
    if (cfgUseDCAzCorrection) {
      double sigma = cfgTrackPIDsettings->get(i, "maxDcaZ") * (dcaZp0 * std::exp(-dcaZp1 * pt) + dcaZp2);
      insideDCAz = std::abs(track.dcaZ()) <= sigma;
    } else {
      insideDCAz = std::abs(track.dcaZ()) <= cfgTrackPIDsettings->get(i, "maxDcaZ");
    }

    return insideDCAxy && insideDCAz;
  }

  /// Base candidate pT for species i: He3/He4 have charge 2, so the tracking pT
  /// (computed for charge 1) is doubled.
  template <typename T>
  float computeSpeciesPt(T const& track, size_t i)
  {
    setTrackParCov(track, mTrackParCov);
    mTrackParCov.setPID(track.pidForTracking());
    return (i == he3 || i == he4) ? 2.f * mTrackParCov.getPt() : mTrackParCov.getPt();
  }

  /// MC-only pT correction for He3/He4 reconstructed tracks (empirical a/b/c
  /// parameterisation per species/charge, plus an antitriton-mistag correction
  /// for He3). No-op when cfgmccorrectionhe4Require is false or species is
  /// neither He3 nor He4.
  template <typename T>
  float applyHe34PtCorrection(T const& track, int pdg, size_t i, float ptIn)
  {
    if (!cfgmccorrectionhe4Require)
      return ptIn;

    double a = 0, b = 0, c = 0;
    int param = -1;
    if (i == he3) {
      param = (pdg > 0) ? 0 : 1;
    } else if (i == he4) {
      param = (pdg > 0) ? 2 : 3;
    }
    if (param >= 0) {
      a = cfgktrackcorrection->get(param, "a");
      b = cfgktrackcorrection->get(param, "b");
      c = cfgktrackcorrection->get(param, "c");
    }

    float ptOut = ptIn;
    if (std::abs(pdg) == particlePdgCodes.at(5)) { // He4/alpha
      ptOut = ptOut + a + b * std::exp(c * ptOut);
    }
    if (std::abs(pdg) == particlePdgCodes.at(4)) { // He3/helion, antitriton-mistag correction
      constexpr int AntiTritonPid = 6;
      if (track.pidForTracking() == AntiTritonPid) {
        ptOut = ptOut - a + b * ptOut - c * ptOut * ptOut;
      }
    }
    return ptOut;
  }

  /// Look up the primaryParticles/particlePdgCodes index matching an (MC) pdg code,
  /// or -1 if none of the configured species match.
  int findParticleTypeIndex(int pdgCode)
  {
    for (size_t j = 0; j < primaryParticles.size(); j++) {
      if (std::abs(pdgCode) == std::abs(particlePdgCodes.at(j)))
        return static_cast<int>(j);
    }
    return -1;
  }

  /// True if particleType is a valid species index with fillsparsh enabled.
  bool speciesEnabled(int particleType)
  {
    return particleType >= 0 && cfgTrackPIDsettings2->get(particleType, "fillsparsh") == 1;
  }

  /// Decay-origin classification shared by processMC's generated-particle loop and
  /// reconstructed-track loop (identical logic in both places).
  template <typename McPart>
  int classifyDecayType(McPart const& mcPart)
  {
    int decayType;
    if (mcPart.isPhysicalPrimary()) {
      decayType = DecayTypePrimary;
      if (mcPart.has_mothers()) {
        for (const auto& motherparticle : mcPart.template mothers_as<aod::McParticles>()) {
          if (std::find(hfMothCodes.begin(), hfMothCodes.end(), std::abs(motherparticle.pdgCode())) != hfMothCodes.end()) {
            decayType = DecayTypeWeak;
            break;
          }
        }
      }
    } else if (mcPart.has_mothers()) {
      decayType = DecayTypeWeak;
    } else {
      decayType = DecayTypeTransport;
    }
    return decayType;
  }

  //====================================================================================
  // CCDB / event / mass-histogram helper functions
  //====================================================================================

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    constexpr float InvalidBField = -990.f;
    auto run3grpTimestamp = bc.timestamp();
    dBz = 0;
    auto* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grpTimestamp);
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      if (bField < InvalidBField) {
        dBz = grpo->getNominalL3Field();
        LOG(info) << "Retrieved GRP for timestamp " << run3grpTimestamp << " with magnetic field of " << dBz << " kZG";
      } else {
        dBz = bField;
      }
    } else {
      auto* grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grpTimestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grpTimestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      if (bField < InvalidBField) {
        dBz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
        LOG(info) << "Retrieved GRP for timestamp " << run3grpTimestamp << " with magnetic field of " << dBz << " kZG";
      } else {
        dBz = bField;
      }
    }
    mRunNumber = bc.runNumber();
  }

  template <typename T>
  void initCollision(const T& collision)
  {
    collHasCandidate = false;
    histos.fill(HIST("histNev"), 0.5);
    collPassedEvSel = collision.sel8() && std::abs(collision.posZ()) < cfgZvertex;
    occupancy = collision.trackOccupancyInTimeRange();
    if (collPassedEvSel) {
      histos.fill(HIST("histNev"), 1.5);
    }
    primVtx.assign({collision.posX(), collision.posY(), collision.posZ()});
    cents.assign({collision.centFT0A(), collision.centFT0C(), collision.centFT0M()});
  }

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
    float p = getRigidity(track);
    float massTOF = p * charge * std::sqrt(1.f / (beta * beta) - 1.f);

    if (species == he4 && cfghe3massrejreq && (massTOF * massTOF > cfgminmassrejection && massTOF * massTOF < cfgmaxmassrejection))
      return;

    float massDiff = cfgmass2 ? (massTOF * massTOF) / (charge * charge) : massTOF / charge;

    float ptMomn = computeSpeciesPt(track, species);
    if (track.sign() > 0) {
      hmass[2 * species]->Fill(ptMomn, massDiff, cent);
    } else if (track.sign() < 0) {
      hmass[2 * species + 1]->Fill(ptMomn, massDiff, cent);
    }
  }

  template <class T>
  void fillhmassnsigma(T const& track, int species, float cent)
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
    float massSquareOverChargeSquare = (massTOF * massTOF) / (charge * charge);

    if (species == he4 && cfghe3massrejreq && (massTOF * massTOF > cfgminmassrejection && massTOF * massTOF < cfgmaxmassrejection))
      return;

    float ptMomn = computeSpeciesPt(track, species);
    if (track.sign() > 0) {
      hmassnsigma[2 * species]->Fill(ptMomn, massSquareOverChargeSquare, cent);
    } else if (track.sign() < 0) {
      hmassnsigma[2 * species + 1]->Fill(ptMomn, massSquareOverChargeSquare, cent);
    }
  }

  //====================================================================================
  // PID / kinematics helper functions
  //====================================================================================

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

    double expBethe{common::BetheBlochAleph(static_cast<double>(particle.charge * rigidity / particle.mass), particle.betheParams[0], particle.betheParams[1], particle.betheParams[2], particle.betheParams[3], particle.betheParams[4])};
    double expSigma{expBethe * particle.resolution};
    return static_cast<float>((track.tpcSignal() - expBethe) / expSigma);
  }

  /// TOF n-sigma for the given species, mirroring getTPCnSigma's species
  /// dispatch. Returns -999 if the track has no TOF signal.
  template <class T>
  float getTOFnSigma(T const& track, PrimParticles& particle)
  {
    if (!track.hasTOF())
      return -999;
    if (particle.name == "pion")
      return track.tofNSigmaPi();
    if (particle.name == "proton")
      return track.tofNSigmaPr();
    if (particle.name == "deuteron")
      return track.tofNSigmaDe();
    if (particle.name == "triton")
      return track.tofNSigmaTr();
    if (particle.name == "helion")
      return track.tofNSigmaHe();
    if (particle.name == "alpha")
      return track.tofNSigmaAl();
    return -999; // fallback if no match
  }

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

  template <class T>
  float getMeanItsClsSize(T const& track)
  {
    int sum = 0, n = 0;
    constexpr int NITSLayers = 8;
    for (int i = 0; i < NITSLayers; i++) {
      sum += (track.itsClusterSizes() >> (4 * i) & 15);
      if (track.itsClusterSizes() >> (4 * i) & 15)
        n++;
    }
    return n > 0 ? static_cast<float>(sum) / n : 0.f;
  }

  template <class T>
  float getRigidity(T const& track)
  {
    if (!cfgRigidityCorrection)
      return track.tpcInnerParam();
    bool hePID = track.pidForTracking() == o2::track::PID::Helium3 || track.pidForTracking() == o2::track::PID::Alpha;
    return hePID ? track.tpcInnerParam() / 2 : track.tpcInnerParam();
  }

  template <class T>
  float getRapidity(T const& track, int species)
  {
    using PtEtaPhiMVector = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double>>;
    int speciesHe3 = 4;
    int speciesHe4 = 5;
    double momn = (species == speciesHe3 || species == speciesHe4) ? 2 * track.pt() : track.pt();
    PtEtaPhiMVector lorentzVectorParticle(momn, track.eta(), track.phi(), particleMasses[species]);
    return lorentzVectorParticle.Rapidity();
  }
}; // end of the task here

//----------------------------------------------------------------------------------------------------------------
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<NucleitpcPbPb>(cfgc)};
}
