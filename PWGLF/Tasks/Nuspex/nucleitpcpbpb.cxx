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
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
constexpr double kBetheBlochDefault[nParticles][nBetheParams]{
  {13.611469, 3.598765, -0.021138, 2.039562, 0.651040, 0.09},    // pion
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},      // proton
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},      // deuteron
  {5.393020, 7.859534, 0.004048, 2.323197, 1.609307, 0.09},      // triton
  {-126.557359, -0.858569, 1.111643, 1.210323, 2.656374, 0.09},  // helion
  {-126.557359, -0.858569, 1.111643, 1.210323, 2.656374, 0.09}}; // alpha
const int nTrkSettings = 19;
static const std::vector<std::string> trackPIDsettingsNames{"useBBparams", "minITSnCls", "minITSnClscos", "minTPCnCls", "maxTPCchi2", "minTPCchi2", "maxITSchi2", "minRigidity", "maxRigidity", "maxTPCnSigma", "TOFrequiredabove", "minTOFmass", "maxTOFmass", "maxDcaXY", "maxDcaZ", "minITSclsSize", "maxITSclsSize", "minTPCnClsCrossedRows", "minReqClusterITSib"};
constexpr double kTrackPIDSettings[nParticles][nTrkSettings]{
  {0, 0, 4, 60, 4.0, 0.5, 100, 0.15, 1.2, 2.5, -1, 0, 100, 2., 2., 0., 1000, 70, 1},
  {1, 0, 4, 70, 4.0, 0.5, 100, 0.20, 4.0, 3.0, -1, 0, 100, 2., 2., 0., 1000, 70, 1},
  {1, 0, 4, 70, 4.0, 0.5, 100, 0.50, 5.0, 3.0, -1, 0, 100, 2., 2., 0., 1000, 70, 1},
  {1, 0, 4, 70, 4.0, 0.5, 100, 0.50, 5.0, 3.0, -1, 0, 100, 2., 2., 0., 1000, 70, 1},
  {1, 0, 4, 75, 4.0, 0.5, 100, 0.50, 5.0, 5.0, -1, 0, 100, 2., 2., 0., 1000, 70, 1},
  {1, 0, 4, 70, 4.0, 0.5, 100, 0.50, 5.0, 5.0, -1, 0, 100, 2., 2., 0., 1000, 70, 1}};
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
std::vector<std::shared_ptr<TH2>> hDeDx;
std::vector<std::shared_ptr<TH2>> hNsigmaPt;
std::vector<std::shared_ptr<TH2>> hmass;
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
  Configurable<bool> cfgmaxITSnSigmaRequire{"cfgmaxITSnSigmaRequire", true, "Require maxITSnSigma Cut for helium"};
  Configurable<bool> cfgminGetMeanItsClsSizeRequire{"cfgminGetMeanItsClsSizeRequire", true, "Require minGetMeanItsClsSize Cut"};
  Configurable<bool> cfgmaxGetMeanItsClsSizeRequire{"cfgmaxGetMeanItsClsSizeRequire", true, "Require maxGetMeanItsClsSize Cut"};
  Configurable<bool> cfgRigidityCutRequire{"cfgRigidityCutRequire", true, "Require Rigidity Cut"};
  Configurable<bool> cfgmassRequire{"cfgmassRequire", true, "Require mass Cuts"};
  Configurable<bool> cfgDCAwithptRequire{"cfgDCAwithptRequire", true, "Require DCA cuts with pt dependance"};
  Configurable<bool> cfgDCAnopt{"cfgDCAnopt", true, "Require DCA cuts without pt dependance"};
  Configurable<bool> cfgTwicemass{"cfgTwicemass", true, "multiply mass by its charge"};
  Configurable<bool> cfgRequirebetaplot{"cfgRequirebetaplot", true, "Require beta plot"};
  Configurable<bool> cfgRequireMCposZ{"cfgRequireMCposZ", true, "Require beta plot"};

  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {kBetheBlochDefault[0], nParticles, nBetheParams, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for light nuclei"};
  Configurable<LabeledArray<double>> cfgTrackPIDsettings{"cfgTrackPIDsettings", {kTrackPIDSettings[0], nParticles, nTrkSettings, particleNames, trackPIDsettingsNames}, "track selection and PID criteria"};
  Configurable<bool> cfgFillDeDxWithCut{"cfgFillDeDxWithCut", true, "Fill with cut beth bloch"};
  Configurable<bool> cfgFillnsigma{"cfgFillnsigma", false, "Fill n-sigma histograms"};
  Configurable<bool> cfgFillmass{"cfgFillmass", false, "Fill mass histograms"};
  Configurable<float> centcut{"centcut", 80.0f, "centrality cut"};
  Configurable<float> cfgCutRapidity{"cfgCutRapidity", 0.5f, "Rapidity range"};
  Configurable<float> cfgZvertex{"cfgZvertex", 10, "Min Z Vertex"};
  Configurable<float> cfgITSnsigma{"cfgITSnsigma", 5, "Max ITS nsigma value"};
  Configurable<float> cfgtpcNClsFound{"cfgtpcNClsFound", 100.0f, "min. no. of tpcNClsFound"};
  Configurable<float> cfgitsNCls{"cfgitsNCls", 2.0f, "min. no. of itsNCls"};
  o2::track::TrackParametrizationWithError<float> mTrackParCov;
  // Binning configuration
  ConfigurableAxis axisMagField{"axisMagField", {10, -10., 10.}, "magnetic field"};
  ConfigurableAxis axisNev{"axisNev", {3, 0., 3.}, "Number of events"};
  ConfigurableAxis axisRigidity{"axisRigidity", {4000, -10., 10.}, "#it{p}^{TPC}/#it{z}"};
  ConfigurableAxis axisdEdx{"axisdEdx", {4000, 0, 4000}, "d#it{E}/d#it{x}"};
  ConfigurableAxis axisCent{"axisCent", {100, 0, 100}, "centrality"};
  ConfigurableAxis axisVtxZ{"axisVtxZ", {120, -20, 20}, "z"};
  ConfigurableAxis axisImpt{"axisImpt", {100, 0, 20}, "impact parameter"};
  ConfigurableAxis ptAxis{"ptAxis", {200, 0, 10}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis axiseta{"axiseta", {100, -1, 1}, "eta"};
  ConfigurableAxis axisrapidity{"axisrapidity", {100, -2, 2}, "rapidity"};
  ConfigurableAxis axismass{"axismass", {100, -10, 10}, "mass^{2}"};
  ConfigurableAxis nsigmaAxis{"nsigmaAxis", {160, -20, 20}, "n#sigma_{#pi^{+}}"};
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
    }
    histos.add("histeta", "histeta", kTH1F, {axiseta});
    histos.add("Tofsignal", "Tofsignal", kTH2F, {axisRigidity, {4000, 0.2, 1.2, "#beta"}});
    histos.add("histDcaZVsPtData_particle", "dcaZ vs Pt (particle)", HistType::kTH2F, {{1000, 0, 20}, {1000, -2.5, 2.5, "dca"}});
    histos.add("histDcaXYVsPtData_particle", "dcaXY vs Pt (particle)", HistType::kTH2F, {{1000, 0, 20}, {1000, -2.0, 2.0, "dca"}});
    histos.add("histDcaZVsPtData_antiparticle", "dcaZ vs Pt (antiparticle)", HistType::kTH2F, {{1000, 0, 20}, {1000, -2.5, 2.5, "dca"}});
    histos.add("histDcaXYVsPtData_antiparticle", "dcaXY vs Pt (antiparticle)", HistType::kTH2F, {{1000, 0, 20}, {1000, -2.0, 2.0, "dca"}});
    hDeDx.resize(2 * nParticles + 2);
    hNsigmaPt.resize(2 * nParticles + 2);
    hmass.resize(2 * nParticles + 2);
    for (int i = 0; i < nParticles + 1; i++) {
      TString histName = i < nParticles ? primaryParticles[i].name : "all";
      if (cfgFillDeDxWithCut) {
        hDeDx[2 * i] = histos.add<TH2>(Form("dedx/histdEdx_%s_Cuts", histName.Data()), ";p_{TPC}/z (GeV/#it{c}); d#it{E}/d#it{x}", HistType::kTH2F, {axisRigidity, axisdEdx});
      }
    }
    for (int i = 0; i < nParticles; i++) {
      TString histName = primaryParticles[i].name;
      if (cfgFillnsigma) {
        hNsigmaPt[2 * i] = histos.add<TH2>(Form("histnsigmaTPC_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); TPCnsigma", HistType::kTH2F, {ptAxis, nsigmaAxis});
        hNsigmaPt[2 * i + 1] = histos.add<TH2>(Form("histnsigmaTPC_anti_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); TPCnsigma", HistType::kTH2F, {ptAxis, nsigmaAxis});
      }
      if (cfgFillmass) {
        hmass[2 * i] = histos.add<TH2>(Form("histmass_pt/histmass_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); mass^{2}", HistType::kTH2F, {ptAxis, axismass});
        hmass[2 * i + 1] = histos.add<TH2>(Form("histmass_ptanti/histmass_%s", histName.Data()), ";p_T{TPC} (GeV/#it{c}); mass^{2}", HistType::kTH2F, {ptAxis, axismass});
      }
    }

    if (doprocessMC) {
      histomc.add("histVtxZgen", "histVtxZgen", kTH1F, {axisVtxZ});
      histomc.add("histEtagen", "histEtagen", kTH1F, {axiseta});
      histomc.add("histPtgenHe3", "histPtgenHe3", kTH1F, {ptAxis});
      histomc.add("histPtgenAntiHe3", "histPtgenAntiHe3", kTH1F, {ptAxis});
      histomc.add("histPtgenHe4", "histPtgenHe4", kTH1F, {ptAxis});
      histomc.add("histPtgenAntiHe4", "histPtgenAntiHe4", kTH1F, {ptAxis});
      // Reconstrcuted eta
      histomc.add("histNevReco", "histNevReco", kTH1F, {axisNev});
      histomc.add("histVtxZReco", "histVtxZReco", kTH1F, {axisVtxZ});
      histomc.add("histCentFT0CReco", "histCentFT0CReco", kTH1F, {axisCent});
      histomc.add("histCentFT0MReco", "histCentFT0MReco", kTH1F, {axisCent});
      histomc.add("histPtRecoHe3", "histPtgenHe3", kTH1F, {ptAxis});
      histomc.add("histPtRecoAntiHe3", "histPtgenAntiHe3", kTH1F, {ptAxis});
      histomc.add("histPtRecoHe4", "histPtgenHe4", kTH1F, {ptAxis});
      histomc.add("histPtRecoAntiHe4", "histPtgenAntiHe4", kTH1F, {ptAxis});
      histomc.add("histDeltaPtVsPtGen", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10}, {1000, -0.5, 0.5, "p_{T}(reco) - p_{T}(gen);p_{T}(reco)"}});
      histomc.add("histPIDtrack", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10, "p_{T}(reco)"}, {9, -0.5, 8.5, "p_{T}(reco) - p_{T}(gen)"}});
      histomc.add("histDeltaPtVsPtGenHe4", " delta pt vs pt rec", HistType::kTH2F, {{1000, 0, 10}, {1000, -0.5, 0.5, "p_{T}(reco) - p_{T}(gen);p_{T}(reco)"}});
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  void findprimaryParticles(aod::TrackAssoc const& tracksByColl, TracksFull const& tracks)
  {
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
      if (std::abs(track.eta()) > cfgCutEta)
        continue;
      for (size_t i = 0; i < primaryParticles.size(); i++) {
        if (std::abs(getRapidity(track, i)) > cfgCutRapidity && cfgRapidityRequire)
          continue;
        if (track.tpcNClsFound() < cfgTrackPIDsettings->get(i, "minTPCnCls") && cfgTPCNClsfoundRequire)
          continue;
        if (((track.tpcNClsCrossedRows() < cfgTrackPIDsettings->get(i, "minTPCnClsCrossedRows")) || track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable()) && cfgTPCNClsCrossedRowsRequire) // o2-linter: disable=magic-number (To be checked)
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
        if (getMeanItsClsSize(track) > cfgTrackPIDsettings->get(i, "maxITSclsSize") && cfgmaxGetMeanItsClsSizeRequire)
          continue;
        if ((getRigidity(track) < cfgTrackPIDsettings->get(i, "minRigidity") || getRigidity(track) > cfgTrackPIDsettings->get(i, "maxRigidity")) && cfgRigidityCutRequire)
          continue;
        float ptMomn;
        setTrackParCov(track, mTrackParCov);
        mTrackParCov.setPID(track.pidForTracking());
        ptMomn = (i == he3 || i == he4) ? 2 * mTrackParCov.getPt() : mTrackParCov.getPt();
        bool insideDCAxy = (std::abs(track.dcaXY()) <= (cfgTrackPIDsettings->get(i, "maxDcaXY") * (0.0105f + 0.0350f / std::pow(ptMomn, 1.1f)))); // o2-linter: disable=magic-number (To be checked)
        if ((!(insideDCAxy) || std::abs(track.dcaZ()) > dcazSigma(ptMomn, cfgTrackPIDsettings->get(i, "maxDcaZ"))) && cfgDCAwithptRequire)
          continue;
        if (track.sign() > 0) {
          histos.fill(HIST("histDcaZVsPtData_particle"), ptMomn, track.dcaZ());
          histos.fill(HIST("histDcaXYVsPtData_particle"), ptMomn, track.dcaXY());
        }
        if (track.sign() < 0) {
          histos.fill(HIST("histDcaZVsPtData_antiparticle"), ptMomn, track.dcaZ());
          histos.fill(HIST("histDcaXYVsPtData_antiparticle"), ptMomn, track.dcaXY());
        }
        float tpcNsigma = getTPCnSigma(track, primaryParticles.at(i));
        if ((std::abs(tpcNsigma) > cfgTrackPIDsettings->get(i, "maxTPCnSigma")) && cfgmaxTPCnSigmaRequire)
          continue;
        float itsSigma = getITSnSigma(track, primaryParticles.at(i));
        if ((std::abs(itsSigma) > cfgITSnsigma) && cfgmaxITSnSigmaRequire)
          continue;
        fillnsigma(track, i);
        filldedx(track, i);
        // TOF selection
        if (!track.hasTOF() && cfgTrackPIDsettings->get(i, "TOFrequiredabove") < 1)
          continue;
        float beta{o2::pid::tof::Beta::GetBeta(track)};
        float charge{1.f + static_cast<float>(i == he3 || i == he4)};
        float tofMasses = getRigidity(track) * charge * std::sqrt(1.f / (beta * beta) - 1.f);
        if ((getRigidity(track) > cfgTrackPIDsettings->get(i, "TOFrequiredabove") && (tofMasses < cfgTrackPIDsettings->get(i, "minTOFmass") || tofMasses > cfgTrackPIDsettings->get(i, "maxTOFmass"))) && cfgmassRequire)
          continue;
        fillhmass(track, i);
      }
      histos.fill(HIST("histeta"), track.eta());
      if (cfgRequirebetaplot) {
        histos.fill(HIST("Tofsignal"), getRigidity(track), o2::pid::tof::Beta::GetBeta(track));
      }
      filldedx(track, nParticles);
    } // track loop
  }
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
      const uint64_t collIdx = collision.globalIndex();
      auto tracksByColl = tracksColl.sliceBy(perCollision, collIdx);
      findprimaryParticles(tracksByColl, tracks);
      if (!collHasCandidate)
        continue;
    }
  }
  PROCESS_SWITCH(NucleitpcPbPb, processData, "data analysis", false);
  //----------------------------------------------------------------------------------------------------------------
  // MC particles
  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  struct McCollInfo {
    bool passedEvSel = false;
  };
  std::vector<McCollInfo> mcCollInfos;
  void processMC(CollisionsFullMC const& collisions, aod::McCollisions const& mcCollisions, TracksFull const& tracks, aod::BCsWithTimestamps const& bcs, aod::McParticles const& particlesMC, aod::McTrackLabels const& trackLabelsMC, aod::McCollisionLabels const& collLabels, aod::TrackAssoc const& tracksColl, CollisionsFull const& colls)
  {
    (void)colls;
    (void)collLabels;
    (void)bcs;
    mcCollInfos.clear();
    mcCollInfos.resize(mcCollisions.size());
    // ----------------------------- Generated particles loop -----------------------------
    for (auto const& mcColl : mcCollisions) {
      if (std::abs(mcColl.posZ()) > cfgZvertex)
        continue;
      histomc.fill(HIST("histVtxZgen"), mcColl.posZ());
      for (auto const& mcParticle : particlesMC) {
        if (mcParticle.mcCollisionId() != mcColl.globalIndex())
          continue;
        if (!mcParticle.isPhysicalPrimary())
          continue;
        int pdgCode = mcParticle.pdgCode();
        if (std::abs(mcParticle.y()) > cfgCutRapidity)
          continue;
        if (std::abs(mcParticle.eta()) > cfgCutEta)
          continue;
        histomc.fill(HIST("histEtagen"), mcParticle.eta());
        float ptScaled = mcParticle.pt();
        if (pdgCode == particlePdgCodes.at(4)) {
          histomc.fill(HIST("histPtgenHe3"), ptScaled);
        } else if (pdgCode == -particlePdgCodes.at(4)) {
          histomc.fill(HIST("histPtgenAntiHe3"), ptScaled);
        } else if (pdgCode == particlePdgCodes.at(5)) {
          histomc.fill(HIST("histPtgenHe4"), ptScaled);
        } else if (pdgCode == -particlePdgCodes.at(5)) {
          histomc.fill(HIST("histPtgenAntiHe4"), ptScaled);
        }
      } // mc track loop generated
    } // mc collision loop generated
    // ----------------------------- Reconstructed track loop -----------------------------
    for (auto const& collision : collisions) {
      auto mcCollIdx = collision.mcCollisionId();
      if (mcCollIdx < 0 || mcCollIdx >= mcCollisions.size())
        continue;
      auto bc = collision.bc_as<aod::BCsWithTimestamps>();
      initCCDB(bc);
      collHasCandidate = false;
      histomc.fill(HIST("histNevReco"), 0.5);
      collPassedEvSel = collision.sel8() && std::abs(collision.posZ()) < cfgZvertex;
      occupancy = collision.trackOccupancyInTimeRange();
      if (!collPassedEvSel)
        continue;
      histomc.fill(HIST("histNevReco"), 1.5);
      histomc.fill(HIST("histVtxZReco"), collision.posZ());
      histomc.fill(HIST("histCentFT0CReco"), collision.centFT0C());
      histomc.fill(HIST("histCentFT0MReco"), collision.centFT0M());
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
      // Loop through TrackAssoc and check if this track belongs to the current collision
      for (const auto& assoc : tracksColl) {
        if (assoc.collisionId() != collision.globalIndex())
          continue;
        const auto& track = tracks.rawIteratorAt(assoc.trackId());
        auto labelRow = trackLabelsMC.iteratorAt(track.globalIndex());
        int label = labelRow.mcParticleId();
        if (label < 0 || label >= particlesMC.size())
          continue;
        auto const& matchedMCParticle = particlesMC.iteratorAt(label);
        int pdg = matchedMCParticle.pdgCode();
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
        if (std::abs(track.eta()) > cfgCutEta)
          continue;
        if (!matchedMCParticle.isPhysicalPrimary())
          continue;
        for (size_t i = 0; i < primaryParticles.size(); i++) {
          if (std::abs(pdg) != std::abs(particlePdgCodes.at(i)))
            continue;
          if (std::abs(getRapidity(track, i)) > cfgCutRapidity && cfgRapidityRequire)
            continue;
          if (track.tpcNClsFound() < cfgTrackPIDsettings->get(i, "minTPCnCls") && cfgTPCNClsfoundRequire)
            continue;
          if (((track.tpcNClsCrossedRows() < cfgTrackPIDsettings->get(i, "minTPCnClsCrossedRows")) || track.tpcNClsCrossedRows() < 0.8 * track.tpcNClsFindable()) && cfgTPCNClsCrossedRowsRequire) // o2-linter: disable=magic-number (To be checked)
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
          if (getMeanItsClsSize(track) > cfgTrackPIDsettings->get(i, "maxITSclsSize") && cfgmaxGetMeanItsClsSizeRequire)
            continue;
          if ((getRigidity(track) < cfgTrackPIDsettings->get(i, "minRigidity") || getRigidity(track) > cfgTrackPIDsettings->get(i, "maxRigidity")) && cfgRigidityCutRequire)
            continue;
          float ptMomn;
          setTrackParCov(track, mTrackParCov);
          mTrackParCov.setPID(track.pidForTracking());
          ptMomn = (i == he3 || i == he4) ? 2 * mTrackParCov.getPt() : mTrackParCov.getPt();
          bool insideDCAxy = (std::abs(track.dcaXY()) <= (cfgTrackPIDsettings->get(i, "maxDcaXY") * (0.0105f + 0.0350f / std::pow(ptMomn, 1.1f)))); // o2-linter: disable=magic-number (To be checked)
          if ((!(insideDCAxy) || std::abs(track.dcaZ()) > dcazSigma(ptMomn, cfgTrackPIDsettings->get(i, "maxDcaZ"))) && cfgDCAwithptRequire)
            continue;
          if (track.sign() > 0) {
            histos.fill(HIST("histDcaZVsPtData_particle"), ptMomn, track.dcaZ());
            histos.fill(HIST("histDcaXYVsPtData_particle"), ptMomn, track.dcaXY());
          }
          if (track.sign() < 0) {
            histos.fill(HIST("histDcaZVsPtData_antiparticle"), ptMomn, track.dcaZ());
            histos.fill(HIST("histDcaXYVsPtData_antiparticle"), ptMomn, track.dcaXY());
          }
          float tpcNsigma = getTPCnSigma(track, primaryParticles.at(i));
          if ((std::abs(tpcNsigma) > cfgTrackPIDsettings->get(i, "maxTPCnSigma")) && cfgmaxTPCnSigmaRequire)
            continue;
          float itsSigma = getITSnSigma(track, primaryParticles.at(i));
          if ((std::abs(itsSigma) > cfgITSnsigma) && cfgmaxITSnSigmaRequire)
            continue;
          fillnsigma(track, i);
          filldedx(track, i);
          if (!track.hasTOF() && cfgTrackPIDsettings->get(i, "TOFrequiredabove") < 1)
            continue;
          float beta{o2::pid::tof::Beta::GetBeta(track)};
          float charge{1.f + static_cast<float>(i == he3 || i == he4)};
          float tofMasses = getRigidity(track) * charge * std::sqrt(1.f / (beta * beta) - 1.f);
          if ((getRigidity(track) > cfgTrackPIDsettings->get(i, "TOFrequiredabove") && (tofMasses < cfgTrackPIDsettings->get(i, "minTOFmass") || tofMasses > cfgTrackPIDsettings->get(i, "maxTOFmass"))) && cfgmassRequire)
            continue;
          fillhmass(track, i);
        }
        histos.fill(HIST("histeta"), track.eta());
        if (cfgRequirebetaplot) {
          histos.fill(HIST("Tofsignal"), getRigidity(track), o2::pid::tof::Beta::GetBeta(track));
        }
        filldedx(track, nParticles);
        /*----------------------------------------------------------------------------------------------------------------*/
        float ptReco;
        setTrackParCov(track, mTrackParCov);
        mTrackParCov.setPID(track.pidForTracking());
        ptReco = (std::abs(pdg) == particlePdgCodes.at(4) || std::abs(pdg) == particlePdgCodes.at(5)) ? 2 * mTrackParCov.getPt() : mTrackParCov.getPt();
        if (pdg == -particlePdgCodes.at(5) && cfgmccorrectionhe4Require) {
          ptReco = ptReco + 0.00765 + 0.503791 * std::exp(-1.10517 * ptReco);
        }

        if (pdg == -particlePdgCodes.at(4) && cfgmccorrectionhe4Require) {
          int pidGuess = track.pidForTracking();
          int antitriton = 6;
          if (pidGuess == antitriton) {
            ptReco = ptReco - 0.464215 + 0.195771 * ptReco - 0.0183111 * ptReco * ptReco;
            //  LOG(info) << "we have he3" << pidGuess;
          }
        }
        float ptGen = matchedMCParticle.pt();
        float deltaPt = ptReco - ptGen;

        if (pdg == -particlePdgCodes.at(4)) {
          histomc.fill(HIST("histDeltaPtVsPtGen"), ptReco, deltaPt);
          histomc.fill(HIST("histPIDtrack"), ptReco, track.pidForTracking());
        }
        if (pdg == -particlePdgCodes.at(5)) {
          histomc.fill(HIST("histDeltaPtVsPtGenHe4"), ptReco, deltaPt);
        }
        if (pdg == particlePdgCodes.at(4)) {
          histomc.fill(HIST("histPtRecoHe3"), ptReco);
        } else if (pdg == -particlePdgCodes.at(4)) {
          histomc.fill(HIST("histPtRecoAntiHe3"), ptReco);
        } else if (pdg == particlePdgCodes.at(5)) {
          histomc.fill(HIST("histPtRecoHe4"), ptReco);
        } else if (pdg == -particlePdgCodes.at(5)) {
          histomc.fill(HIST("histPtRecoAntiHe4"), ptReco);
        }
      } // Track loop
    } // Collision loop
  }
  PROCESS_SWITCH(NucleitpcPbPb, processMC, "MC reco+gen analysis", true);
  //-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
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
  void filldedx(T const& track, int species)
  {
    const float rigidity = getRigidity(track);
    if (track.tpcNClsFound() < cfgtpcNClsFound || track.itsNCls() < cfgitsNCls)
      return;
    if (cfgFillDeDxWithCut) {
      hDeDx[2 * species]->Fill(track.sign() * rigidity, track.tpcSignal());
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  void fillnsigma(T const& track, int species)
  {
    if (track.tpcNClsFound() < cfgtpcNClsFound || track.itsNCls() < cfgitsNCls)
      return;
    if (cfgFillnsigma) {
      int i = species;
      const float tpcNsigma = getTPCnSigma(track, primaryParticles.at(i));
      float ptMomn;
      setTrackParCov(track, mTrackParCov);
      mTrackParCov.setPID(track.pidForTracking());
      ptMomn = (i == he3 || i == he4) ? 2 * mTrackParCov.getPt() : mTrackParCov.getPt();
      if (track.sign() > 0) {
        hNsigmaPt[2 * species]->Fill(ptMomn, tpcNsigma);
      }
      if (track.sign() < 0) {
        hNsigmaPt[2 * species + 1]->Fill(ptMomn, tpcNsigma);
      }
    }
  }
  //----------------------------------------------------------------------------------------------------------------
  template <class T>
  void fillhmass(T const& track, int species)
  {
    if (track.tpcNClsFound() < cfgtpcNClsFound || track.itsNCls() < cfgitsNCls)
      return;
    if (!track.hasTOF() || !cfgFillmass)
      return;
    float beta{o2::pid::tof::Beta::GetBeta(track)};
    if (beta <= 0.f || beta >= 1.f)
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
      hmass[2 * species]->Fill(ptMomn, massDiff * massDiff);
    } else if (track.sign() < 0) {
      hmass[2 * species + 1]->Fill(ptMomn, massDiff * massDiff);
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
