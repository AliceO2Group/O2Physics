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

/// \file   kinkBuilder.cxx
/// \brief Builder task for kink decay topologies using ITS standalone tracks for the mother
/// \author Francesco Mazzaschi <francesco.mazzaschi@cern.ch>

#include "PWGLF/DataModel/LFKinkDecayTables.h"
#include "PWGLF/Utils/svPoolCreator.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DCAFitter/DCAFitterN.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/MatLayerCylSet.h>
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
#include <MathUtils/Primitive2D.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/Vertex.h>

#include <TH2.h>
#include <TPDGCode.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using std::array;
using VBracket = o2::math_utils::Bracket<int>;
using TracksFull = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU>;
using TracksFullMc = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::McTrackLabels>;
using CollisionsCent = soa::Join<aod::Collisions, aod::CentFT0Cs>;
using McRecoCollisions = soa::Join<aod::Collisions, aod::McCollisionLabels>;
using McRecoCollisionsCent = soa::Join<aod::Collisions, aod::McCollisionLabels, aod::CentFT0Cs>;

enum KinkDecayType { SigmaMinusToPiMinusNeutron = 0,
                     SigmaPlusToPiPlusNeutron,
                     SigmaPlusToProtonPi0,
                     NMatchedDecays };

namespace
{
constexpr std::array<float, 7> LayerRadii{2.33959f, 3.14076f, 3.91924f, 19.6213f, 24.5597f, 34.388f, 39.3329f};
constexpr double betheBlochDefault[1][6]{{-1.e32, -1.e32, -1.e32, -1.e32, -1.e32, -1.e32}};
static const std::vector<std::string> betheBlochParNames{"p0", "p1", "p2", "p3", "p4", "resolution"};
static const std::vector<std::string> particleNames{"Daughter"};
constexpr int DepthMcMatchMax = 1; // Max depth for MC matching

std::shared_ptr<TH2> h2ClsMapPtMoth;
std::shared_ptr<TH2> h2ClsMapPtDaug;
std::shared_ptr<TH2> h2DeDxDaugSel;
std::shared_ptr<TH2> h2KinkAnglePt;
std::shared_ptr<TH2> h2MothMassPt;
std::shared_ptr<TH2> hSelMotherQA;
std::shared_ptr<TH2> hSelDaugQA;
std::shared_ptr<TH2> hSelKinkedTrackQA;
std::shared_ptr<TH2> hMothDaughSignsInit;
std::shared_ptr<TH2> hMothDaughSignsFinal;
std::shared_ptr<TH2> h2ItsClsMothBeforeSel;
std::shared_ptr<TH2> h2ItsClsDaugBeforeSel;
std::shared_ptr<TH2> hZDiff;
std::shared_ptr<TH2> hPhiDiff;
std::shared_ptr<TH2> hDCAMothToPV;
std::shared_ptr<TH2> hDCADaugToPV;
std::shared_ptr<TH2> hMothDecRad2;
std::shared_ptr<TH2> hGenCandidates;
std::shared_ptr<TH2> hRecCandidates;
std::array<std::shared_ptr<TH2>, NMatchedDecays> hGenPtKinkAngle;
std::array<std::shared_ptr<TH2>, NMatchedDecays> hRecPtKinkAngle;
} // namespace

struct kinkCandidate {

  float recoPtMoth() const { return std::hypot(momMoth[0], momMoth[1]); }
  float recoPhiMoth() const { return std::atan2(momMoth[1], momMoth[0]); }
  float recoEtaMoth() const { return std::asinh(momMoth[2] / recoPtMoth()); }

  float recoPtDaug() const { return std::hypot(momDaug[0], momDaug[1]); }
  float recoPhiDaug() const { return std::atan2(momDaug[1], momDaug[0]); }
  float recoEtaDaug() const { return std::asinh(momDaug[2] / recoPtDaug()); }

  int mothTrackID;
  int daugTrackID;
  int collisionID;

  int mothSign;
  std::array<float, 3> momMoth = {-999, -999, -999};
  std::array<float, 3> momDaug = {-999, -999, -999};
  std::array<float, 3> primVtx = {-999, -999, -999};
  std::array<float, 3> decVtx = {-999, -999, -999};

  float dcaKinkTopo = -999;
  float nSigmaTPCDaug = -999;
  float nSigmaTOFDaug = -999;
  float dcaXYdaug = -999;
  float dcaXYmoth = -999;
  float kinkAngle = -999;
};

struct kinkBuilder {

  enum PartType { kSigmaMinus = 0,
                  kHypertriton,
                  kHyperhelium4sigma };

  Produces<aod::KinkCands> outputDataTable;
  Produces<aod::KinkCandsUnbound> outputDataTableUB;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  Configurable<int> hypoMoth{"hypoMoth", kSigmaMinus, "Mother particle hypothesis"};
  Configurable<bool> fillDebugTable{"fillDebugTable", false, "If true, fill the debug table with all candidates unbound"};
  // Selection criteria
  Configurable<float> maxDCAMothToPV{"maxDCAMothToPV", 0.1, "Max DCA of the mother to the PV"};
  Configurable<float> minDCADaugToPV{"minDCADaugToPV", 0., "Min DCA of the daughter to the PV"};
  Configurable<float> minPtMoth{"minPtMoth", 0.5, "Minimum pT of the hypercandidate"};
  Configurable<float> maxZDiff{"maxZDiff", 20., "Max z difference between the kink daughter and the mother"};
  Configurable<float> maxPhiDiff{"maxPhiDiff", 100, "Max phi difference between the kink daughter and the mother"};
  Configurable<float> timeMarginNS{"timeMarginNS", 600, "Additional time res tolerance in ns"};
  Configurable<float> etaMax{"etaMax", 1., "eta daughter"};
  Configurable<float> nTPCClusMinDaug{"nTPCClusMinDaug", 80, "daug NTPC clusters cut"};
  Configurable<bool> askTOFforDaug{"askTOFforDaug", false, "If true, ask for TOF signal"};
  Configurable<bool> doSVRadiusCut{"doSVRadiusCut", true, "If true, apply the cut on the radius of the secondary vertex and tracksIU"};
  Configurable<bool> updateMothTrackUsePV{"updateMothTrackUsePV", false, "If true, update the mother track parameters using the primary vertex"};

  o2::vertexing::DCAFitterN<2> fitter;
  o2::base::MatLayerCylSet* lut = nullptr;
  int nItsInnerBarrelLayers = 3;
  int nItsOuterBarrelLayers = 4;
  int nItsTotalLayers = 7;
  int nCoords = 7;
  int itsChi2NClMax = 36;

  // constants
  float radToDeg = o2::constants::math::Rad2Deg;
  svPoolCreator svCreator;

  // bethe bloch parameters
  Configurable<LabeledArray<double>> cfgBetheBlochParams{"cfgBetheBlochParams", {betheBlochDefault[0], 1, 6, particleNames, betheBlochParNames}, "TPC Bethe-Bloch parameterisation for charged daughter"};
  Configurable<int> cfgMaterialCorrection{"cfgMaterialCorrection", static_cast<int>(o2::base::Propagator::MatCorrType::USEMatCorrNONE), "Type of material correction"};
  Configurable<float> customVertexerTimeMargin{"customVertexerTimeMargin", 800, "Time margin for custom vertexer (ns)"};
  Configurable<bool> skipAmbiTracks{"skipAmbiTracks", false, "Skip ambiguous tracks"};
  Configurable<bool> unlikeSignBkg{"unlikeSignBkg", false, "Use unlike sign background"};
  Configurable<bool> skipBkgCands{"skipBkgCands", true, "Skip bkg candidates in MC process"};

  // Centrality selection
  Configurable<float> centralityMin{"centralityMin", 0., "Minimum centrality accepted in SP/EP computation (not applied in resolution process)"};
  Configurable<float> centralityMax{"centralityMax", 100., "Maximum centrality accepted in SP/EP computation (not applied in resolution process)"};

  // CCDB options
  Configurable<std::string> ccdbPath{"ccdbPath", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};

  // histogram axes
  ConfigurableAxis rigidityBins{"rigidityBins", {200, -10.f, 10.f}, "Binning for rigidity #it{p}^{TPC}/#it{z}"};
  ConfigurableAxis dedxBins{"dedxBins", {1000, 0.f, 1000.f}, "Binning for dE/dx"};
  ConfigurableAxis nSigmaBins{"nSigmaBins", {200, -5.f, 5.f}, "Binning for n sigma"};
  ConfigurableAxis zDiffBins{"zDiffBins", {300, 0.f, 30.f}, "Binning for z difference"};
  ConfigurableAxis phiDiffBins{"phiDiffBins", {360, 0.f, 360.f}, "Binning for phi difference"};
  ConfigurableAxis dcaMothToPVBins{"dcaMothToPVBins", {120, 0.f, 60.f}, "Binning for DCA of mother to PV"};
  ConfigurableAxis dcaDaugToPVBins{"dcaDaugToPVBins", {200, 0.f, 20.f}, "Binning for DCA of daughter to PV"};
  ConfigurableAxis mothDecRad2Bins{"mothDecRad2Bins", {100, 0, 100}, "Binning for mother decay radius squared"};

  // Filter collisions with centrality information (for the process with centrality selection)
  using CollisionsCentSel = soa::Filtered<CollisionsCent>;
  using McRecoCollisionsCentSel = soa::Filtered<McRecoCollisionsCent>;

  Filter filterCentrality = aod::cent::centFT0C >= centralityMin && aod::cent::centFT0C <= centralityMax;

  // std vector of candidates
  std::vector<kinkCandidate> kinkCandidates;

  HistogramRegistry qaRegistry{"QA", {}, OutputObjHandlingPolicy::AnalysisObject};

  int mRunNumber;
  float mBz;
  std::array<float, 6> mBBparamsDaug;

  // mother and daughter tracks' properties (absolute charge and mass)
  int charge = 1;
  std::vector<int> mothMatchPdgCodes{};
  float mothMass = o2::constants::physics::MassSigmaMinus;
  float chargedDauMass = o2::constants::physics::MassPiMinus;
  float neutDauMass = o2::constants::physics::MassNeutron;

  void init(InitContext const&)
  {
    if (hypoMoth == kSigmaMinus) {
      charge = 1;
      mothMatchPdgCodes = {PDG_t::kSigmaMinus, PDG_t::kSigmaPlus};
      mothMass = o2::constants::physics::MassSigmaMinus;
      chargedDauMass = o2::constants::physics::MassPiMinus;
      neutDauMass = o2::constants::physics::MassNeutron;
    } else if (hypoMoth == kHypertriton) {
      charge = 1;
      mothMass = o2::constants::physics::MassHyperTriton;
      chargedDauMass = o2::constants::physics::MassTriton;
      neutDauMass = o2::constants::physics::MassPi0;
    } else if (hypoMoth == kHyperhelium4sigma) {
      charge = 2;
      mothMass = o2::constants::physics::MassHyperHelium4;
      chargedDauMass = o2::constants::physics::MassAlpha;
      neutDauMass = o2::constants::physics::MassPi0;
    }

    // dummy values, 1 for mother, 0 for daughter
    svCreator.setPDGs(1, 0);

    mRunNumber = 0;
    mBz = 0;

    ccdb->setURL(ccdbPath);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    fitter.setPropagateToPCA(true);
    fitter.setMaxR(200.);
    fitter.setMinParamChange(1e-3);
    fitter.setMinRelChi2Change(0.9);
    fitter.setMaxDZIni(1e9);
    fitter.setMaxChi2(1e9);
    fitter.setUseAbsDCA(true);

    svCreator.setTimeMargin(customVertexerTimeMargin);
    if (skipAmbiTracks) {
      svCreator.setSkipAmbiTracks();
    }

    const AxisSpec itsClusterMapAxis(128, 0, 127, "ITS cluster map");
    const AxisSpec rigidityAxis{rigidityBins, "#it{p}^{TPC}/#it{z}"};
    const AxisSpec ptAxis{rigidityBins, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec absPtAxis{1000, 0, 10, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec kinkAngleAxis{100, 0, 180, "#theta_{kink} (deg)"};
    const AxisSpec dedxAxis{dedxBins, "d#it{E}/d#it{x}"};

    AxisSpec massAxis(100, 1.1, 1.4, "m (GeV/#it{c}^{2})");
    if (hypoMoth == kSigmaMinus) {
      massAxis = AxisSpec{100, 1.1, 1.4, "m (GeV/#it{c}^{2})"};
    } else if (hypoMoth == kHypertriton) {
      massAxis = AxisSpec{100, 2.94, 3.2, "m (GeV/#it{c}^{2})"};
    } else if (hypoMoth == kHyperhelium4sigma) {
      massAxis = AxisSpec{100, 3.85, 4.25, "m (GeV/#it{c}^{2})"};
    }

    h2DeDxDaugSel = qaRegistry.add<TH2>("h2DeDxDaugSel", "h2DeDxDaugSel; p_{TPC}/z (GeV/#it{c}); dE/dx", HistType::kTH2F, {rigidityAxis, dedxAxis});
    h2KinkAnglePt = qaRegistry.add<TH2>("h2KinkAnglePt", "h2KinkAnglePt; p_{T} (GeV/#it{c}); #theta_{kink} (deg)", HistType::kTH2F, {ptAxis, kinkAngleAxis});
    h2MothMassPt = qaRegistry.add<TH2>("h2MothMassPt", "h2MothMassPt; p_{T} (GeV/#it{c}); m (GeV/#it{c}^{2})", HistType::kTH2F, {ptAxis, massAxis});
    h2ClsMapPtMoth = qaRegistry.add<TH2>("h2ClsMapPtMoth", "h2ClsMapPtMoth; p_{T} (GeV/#it{c}); ITS cluster map", HistType::kTH2F, {ptAxis, itsClusterMapAxis});
    h2ClsMapPtDaug = qaRegistry.add<TH2>("h2ClsMapPtDaug", "h2ClsMapPtDaug; p_{T} (GeV/#it{c}); ITS cluster map", HistType::kTH2F, {ptAxis, itsClusterMapAxis});
    h2ItsClsMothBeforeSel = qaRegistry.add<TH2>("h2ItsClsMothBeforeSel", "h2ItsClsMothBeforeSel; Tot Cls; IB clusters", HistType::kTH2F, {{8, -0.5, 7.5}, {4, -0.5, 3.5}});
    h2ItsClsDaugBeforeSel = qaRegistry.add<TH2>("h2ItsClsDaugBeforeSel", "h2ItsClsDaugBeforeSel; Tot Cls; IB clusters", HistType::kTH2F, {{8, -0.5, 7.5}, {4, -0.5, 3.5}});

    // QA
    hSelMotherQA = qaRegistry.add<TH2>("hSelMotherQA", "hSelMotherQA;;Q_{Mother}", {HistType::kTH2F, {{9, -0.5, 8.5}, {2, -0.5, 1.5}}});
    hSelMotherQA->GetYaxis()->SetBinLabel(1, "-");
    hSelMotherQA->GetYaxis()->SetBinLabel(2, "+");
    hSelMotherQA->GetXaxis()->SetBinLabel(1, "All");
    hSelMotherQA->GetXaxis()->SetBinLabel(2, "Has Collision");
    hSelMotherQA->GetXaxis()->SetBinLabel(3, "Has ITS");
    hSelMotherQA->GetXaxis()->SetBinLabel(4, "Has NO TPC");
    hSelMotherQA->GetXaxis()->SetBinLabel(5, "Has NO TOF");
    hSelMotherQA->GetXaxis()->SetBinLabel(6, "ITS cls sel");
    hSelMotherQA->GetXaxis()->SetBinLabel(7, "ITS IB cls sel");
    hSelMotherQA->GetXaxis()->SetBinLabel(8, "ITS chi2 N cls sel");
    hSelMotherQA->GetXaxis()->SetBinLabel(9, "Pt sel");

    // QA
    hSelDaugQA = qaRegistry.add<TH2>("hSelDaugQA", "hSelDaugQA;;Q_{Daug}", {HistType::kTH2F, {{8, -0.5, 7.5}, {2, -0.5, 1.5}}});
    hSelDaugQA->GetYaxis()->SetBinLabel(1, "-");
    hSelDaugQA->GetYaxis()->SetBinLabel(2, "+");
    hSelDaugQA->GetXaxis()->SetBinLabel(1, "All");
    hSelDaugQA->GetXaxis()->SetBinLabel(2, "Has ITS");
    hSelDaugQA->GetXaxis()->SetBinLabel(3, "Has TPC");
    hSelDaugQA->GetXaxis()->SetBinLabel(4, "Has TOF");
    hSelDaugQA->GetXaxis()->SetBinLabel(5, "ITS tot cls sel");
    hSelDaugQA->GetXaxis()->SetBinLabel(6, "ITS IB cls sel");
    hSelDaugQA->GetXaxis()->SetBinLabel(7, "TPC Crossed rows vs findable sel");
    hSelDaugQA->GetXaxis()->SetBinLabel(8, "TPC cls sel");

    hSelKinkedTrackQA = qaRegistry.add<TH2>("hSelKinkedTrackQA", "hSelKinkedTrackQA;;(Q_{Mother}, Q_{Daughter})", {HistType::kTH2F, {{11, -0.5, 10.5}, {4, -0.5, 3.5}}});
    hSelKinkedTrackQA->GetYaxis()->SetBinLabel(1, "(+,+)");
    hSelKinkedTrackQA->GetYaxis()->SetBinLabel(2, "(-,-)");
    hSelKinkedTrackQA->GetYaxis()->SetBinLabel(3, "(+,-)");
    hSelKinkedTrackQA->GetYaxis()->SetBinLabel(4, "(-,+)");
    hSelKinkedTrackQA->GetXaxis()->SetBinLabel(1, "All");
    hSelKinkedTrackQA->GetXaxis()->SetBinLabel(2, "MaxDCAMothToPV");
    hSelKinkedTrackQA->GetXaxis()->SetBinLabel(3, "MaxZDiff");
    hSelKinkedTrackQA->GetXaxis()->SetBinLabel(4, "MaxPhiDiff");
    hSelKinkedTrackQA->GetXaxis()->SetBinLabel(5, "MinDCADaugToPV");
    hSelKinkedTrackQA->GetXaxis()->SetBinLabel(6, "UpdateMothTrackUsePV");
    hSelKinkedTrackQA->GetXaxis()->SetBinLabel(7, "DCAFitter");
    hSelKinkedTrackQA->GetXaxis()->SetBinLabel(8, "PropagateTracksToVertex");
    hSelKinkedTrackQA->GetXaxis()->SetBinLabel(9, "MothDecRad2InsideL4");
    hSelKinkedTrackQA->GetXaxis()->SetBinLabel(10, "MothDaugLastLayers");
    hSelKinkedTrackQA->GetXaxis()->SetBinLabel(11, "MothLastLayerRadiusCheck");
    hMothDaughSignsInit = qaRegistry.add<TH2>("hMothDaughSignsInit", "hMothDaughSignsInit; Sign Mother; Sign Daughter", {HistType::kTH2F, {{3, -1.5, 1.5}, {3, -1.5, 1.5}}});
    hMothDaughSignsFinal = qaRegistry.add<TH2>("hMothDaughSignsFinal", "hMothDaughSignsFinal; Sign Mother; Sign Daughter", {HistType::kTH2F, {{3, -1.5, 1.5}, {3, -1.5, 1.5}}});
    hZDiff = qaRegistry.add<TH2>("hZDiff", "hZDiff; #Delta z (#mu m);(Q_{Mother}, Q_{Daughter})", HistType::kTH2F, {zDiffBins, {4, -0.5, 3.5}});
    hZDiff->GetYaxis()->SetBinLabel(1, "(+,+)");
    hZDiff->GetYaxis()->SetBinLabel(2, "(-,-)");
    hZDiff->GetYaxis()->SetBinLabel(3, "(+,-)");
    hZDiff->GetYaxis()->SetBinLabel(4, "(-,+)");
    hPhiDiff = qaRegistry.add<TH2>("hPhiDiff", "hPhiDiff; #Delta #phi (rad);(Q_{Mother}, Q_{Daughter})", HistType::kTH2F, {phiDiffBins, {4, -0.5, 3.5}});
    hPhiDiff->GetYaxis()->SetBinLabel(1, "(+,+)");
    hPhiDiff->GetYaxis()->SetBinLabel(2, "(-,-)");
    hPhiDiff->GetYaxis()->SetBinLabel(3, "(+,-)");
    hPhiDiff->GetYaxis()->SetBinLabel(4, "(-,+)");
    hDCAMothToPV = qaRegistry.add<TH2>("hDCAMothToPV", "hDCAMothToPV; DCA moth to PV;(Q_{Mother}, Q_{Daughter})", HistType::kTH2F, {dcaMothToPVBins, {4, -0.5, 3.5}});
    hDCAMothToPV->GetYaxis()->SetBinLabel(1, "(+,+)");
    hDCAMothToPV->GetYaxis()->SetBinLabel(2, "(-,-)");
    hDCAMothToPV->GetYaxis()->SetBinLabel(3, "(+,-)");
    hDCAMothToPV->GetYaxis()->SetBinLabel(4, "(-,+)");
    hDCADaugToPV = qaRegistry.add<TH2>("hDCADaugToPV", "hDCADaugToPV; DCA daug to PV;(Q_{Mother}, Q_{Daughter})", HistType::kTH2F, {dcaDaugToPVBins, {4, -0.5, 3.5}});
    hDCADaugToPV->GetYaxis()->SetBinLabel(1, "(+,+)");
    hDCADaugToPV->GetYaxis()->SetBinLabel(2, "(-,-)");
    hDCADaugToPV->GetYaxis()->SetBinLabel(3, "(+,-)");
    hDCADaugToPV->GetYaxis()->SetBinLabel(4, "(-,+)");
    hMothDecRad2 = qaRegistry.add<TH2>("hMothDecRad2", "hMothDecRad2; Mother Dec. Radius Squared (cm^{2});(Q_{Mother}, Q_{Daughter})", HistType::kTH2F, {mothDecRad2Bins, {4, -0.5, 3.5}});
    hMothDecRad2->GetYaxis()->SetBinLabel(1, "(+,+)");
    hMothDecRad2->GetYaxis()->SetBinLabel(2, "(-,-)");
    hMothDecRad2->GetYaxis()->SetBinLabel(3, "(+,-)");
    hMothDecRad2->GetYaxis()->SetBinLabel(4, "(-,+)");

    mBBparamsDaug[0] = cfgBetheBlochParams->get("Daughter", Form("p%i", 0));
    mBBparamsDaug[1] = cfgBetheBlochParams->get("Daughter", Form("p%i", 1));
    mBBparamsDaug[2] = cfgBetheBlochParams->get("Daughter", Form("p%i", 2));
    mBBparamsDaug[3] = cfgBetheBlochParams->get("Daughter", Form("p%i", 3));
    mBBparamsDaug[4] = cfgBetheBlochParams->get("Daughter", Form("p%i", 4));
    mBBparamsDaug[5] = cfgBetheBlochParams->get("Daughter", "resolution");

    if (doprocessMc || doprocessMcWCent) {
      if (skipBkgCands) {
        hRecCandidates = qaRegistry.add<TH2>("hRecCandidates", "hRecCandidates;Counts;", {HistType::kTH2F, {{NMatchedDecays, -0.5, static_cast<float>(NMatchedDecays) - 0.5}, absPtAxis}});
        hRecCandidates->GetXaxis()->SetBinLabel(1, "#Sigma^{-} #rightarrow n#pi^{-}");
        hRecCandidates->GetXaxis()->SetBinLabel(2, "#Sigma^{+} #rightarrow n#pi^{+}");
        hRecCandidates->GetXaxis()->SetBinLabel(3, "#Sigma^{+} #rightarrow p#pi^{0}");
      }
      hGenCandidates = qaRegistry.add<TH2>("hGenCandidates", "hGenCandidates;Counts;", {HistType::kTH2F, {{NMatchedDecays, -0.5, static_cast<float>(NMatchedDecays) - 0.5}, absPtAxis}});
      hGenCandidates->GetXaxis()->SetBinLabel(1, "#Sigma^{-} #rightarrow n#pi^{-}");
      hGenCandidates->GetXaxis()->SetBinLabel(2, "#Sigma^{+} #rightarrow n#pi^{+}");
      hGenCandidates->GetXaxis()->SetBinLabel(3, "#Sigma^{+} #rightarrow p#pi^{0}");
      hRecPtKinkAngle[SigmaMinusToPiMinusNeutron] = qaRegistry.add<TH2>("hRecPtKinkAngleSigmaMinusToPiMinusNeutron", "Rec Pt vs KinkAngle #Sigma^{-} #rightarrow #pi^{-}n;Counts;", {HistType::kTH2F, {absPtAxis, kinkAngleAxis}});
      hRecPtKinkAngle[SigmaPlusToPiPlusNeutron] = qaRegistry.add<TH2>("hRecPtKinkAngleSigmaPlusToPiPlusNeutron", "Rec Pt vs KinkAngle #Sigma^{+} #rightarrow #pi^{+}n;Counts;", {HistType::kTH2F, {absPtAxis, kinkAngleAxis}});
      hRecPtKinkAngle[SigmaPlusToProtonPi0] = qaRegistry.add<TH2>("hRecPtKinkAngleSigmaPlusToProtonPi0", "Rec Pt vs KinkAngle #Sigma^{+} #rightarrow p#pi^{0};Counts;", {HistType::kTH2F, {absPtAxis, kinkAngleAxis}});
      hGenPtKinkAngle[SigmaMinusToPiMinusNeutron] = qaRegistry.add<TH2>("hGenPtKinkAngleSigmaMinusToPiMinusNeutron", "Gen Pt vs KinkAngle #Sigma^{-} #rightarrow #pi^{-}n;Counts;", {HistType::kTH2F, {absPtAxis, kinkAngleAxis}});
      hGenPtKinkAngle[SigmaPlusToPiPlusNeutron] = qaRegistry.add<TH2>("hGenPtKinkAngleSigmaPlusToPiPlusNeutron", "Gen Pt vs KinkAngle #Sigma^{+} #rightarrow #pi^{+}n;Counts;", {HistType::kTH2F, {absPtAxis, kinkAngleAxis}});
      hGenPtKinkAngle[SigmaPlusToProtonPi0] = qaRegistry.add<TH2>("hGenPtKinkAngleSigmaPlusToProtonPi0", "Gen Pt vs KinkAngle #Sigma^{+} #rightarrow p#pi^{0};Counts;", {HistType::kTH2F, {absPtAxis, kinkAngleAxis}});
    }
  }

  template <typename T>
  bool selectMothTrack(const T& candidate)
  {
    bool isPositive = candidate.sign() > 0;
    hSelMotherQA->Fill(0.f, isPositive);

    if (!candidate.has_collision())
      return false;
    hSelMotherQA->Fill(1.f, isPositive);

    if (!candidate.hasITS())
      return false;
    hSelMotherQA->Fill(2.f, isPositive);

    if (candidate.hasTPC())
      return false;
    hSelMotherQA->Fill(3.f, isPositive);

    if (candidate.hasTOF())
      return false;
    hSelMotherQA->Fill(4.f, isPositive);

    h2ItsClsMothBeforeSel->Fill(candidate.itsNCls(), candidate.itsNClsInnerBarrel());
    if (candidate.itsNCls() >= nItsTotalLayers - 1)
      return false;
    hSelMotherQA->Fill(5.f, isPositive);

    if (candidate.itsNClsInnerBarrel() != nItsInnerBarrelLayers)
      return false;
    hSelMotherQA->Fill(6.f, isPositive);

    if (candidate.itsChi2NCl() >= itsChi2NClMax)
      return false;
    hSelMotherQA->Fill(7.f, isPositive);

    if (candidate.pt() <= minPtMoth)
      return false;
    hSelMotherQA->Fill(8.f, isPositive);

    return true;
  }

  template <typename T>
  bool selectDaugTrack(const T& candidate)
  {
    bool isPositive = candidate.sign() > 0;
    hSelDaugQA->Fill(0.f, isPositive);

    if (!candidate.hasITS())
      return false;
    hSelDaugQA->Fill(1.f, isPositive);

    if (!candidate.hasTPC())
      return false;
    hSelDaugQA->Fill(2.f, isPositive);

    if (askTOFforDaug && !candidate.hasTOF())
      return false;
    hSelDaugQA->Fill(3.f, isPositive);

    h2ItsClsDaugBeforeSel->Fill(candidate.itsNCls(), candidate.itsNClsInnerBarrel());
    if (candidate.itsNClsInnerBarrel() != 0)
      return false;
    hSelDaugQA->Fill(4.f, isPositive);

    if (candidate.itsNCls() >= nItsOuterBarrelLayers)
      return false;
    hSelDaugQA->Fill(5.f, isPositive);

    if (candidate.tpcNClsCrossedRows() <= 0.8 * candidate.tpcNClsFindable())
      return false;
    hSelDaugQA->Fill(6.f, isPositive);

    if (candidate.tpcNClsFound() <= nTPCClusMinDaug)
      return false;
    hSelDaugQA->Fill(7.f, isPositive);

    return true;
  }

  template <class TSVCand, class TTracks, class TColls>
  void buildKinkCand(const TSVCand& svCand, const TTracks& tracks, bool& isAccepted, const TColls&)
  {
    isAccepted = false;
    kinkCandidate kinkCand;

    auto trackMoth = tracks.rawIteratorAt(svCand.tr0Idx);
    auto trackDaug = tracks.rawIteratorAt(svCand.tr1Idx);

    // Fill Selections QA histo
    int chargeCombSvCand = 2 * unlikeSignBkg + (trackMoth.sign() == -1 ? 1 : 0);
    hSelKinkedTrackQA->Fill(0.f, chargeCombSvCand); // all candidates bin
    hMothDaughSignsInit->Fill(trackMoth.sign(), trackDaug.sign());

    auto const& collision = trackMoth.template collision_as<TColls>();
    auto const& bc = collision.template bc_as<aod::BCs>();
    initCCDB(bc);

    o2::dataformats::VertexBase primaryVertex;
    primaryVertex.setPos({collision.posX(), collision.posY(), collision.posZ()});
    primaryVertex.setCov(collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ());
    kinkCand.primVtx = {primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()};

    o2::track::TrackParCov trackParCovMoth = getTrackParCov(trackMoth);
    o2::track::TrackParCov trackParCovMothPV{trackParCovMoth};
    o2::base::Propagator::Instance()->PropagateToXBxByBz(trackParCovMoth, LayerRadii[trackMoth.itsNCls() - 1]);

    std::array<float, 2> dcaInfoMoth;
    o2::base::Propagator::Instance()->propagateToDCABxByBz({primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, trackParCovMothPV, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfoMoth);

    hDCAMothToPV->Fill(std::abs(dcaInfoMoth[0]), chargeCombSvCand);
    if (std::abs(dcaInfoMoth[0]) > maxDCAMothToPV) {
      return;
    }
    hSelKinkedTrackQA->Fill(1.f, chargeCombSvCand); // MaxDCAMothToPV cut bin

    o2::track::TrackParCov trackParCovDaug = getTrackParCov(trackDaug);

    // check if the kink daughter is close to the mother
    hZDiff->Fill(std::abs(trackParCovMoth.getZ() - trackParCovDaug.getZ()), chargeCombSvCand);
    if (std::abs(trackParCovMoth.getZ() - trackParCovDaug.getZ()) > maxZDiff) {
      return;
    }
    hSelKinkedTrackQA->Fill(2.f, chargeCombSvCand); // MaxZDiff cut bin

    hPhiDiff->Fill(std::abs(trackParCovMoth.getPhi() - trackParCovDaug.getPhi()) * radToDeg, chargeCombSvCand);
    if ((std::abs(trackParCovMoth.getPhi() - trackParCovDaug.getPhi()) * radToDeg) > maxPhiDiff) {
      return;
    }
    hSelKinkedTrackQA->Fill(3.f, chargeCombSvCand); // MaxPhiDiff cut bin

    // propagate to PV
    std::array<float, 2> dcaInfoDaug;
    o2::base::Propagator::Instance()->propagateToDCABxByBz({primaryVertex.getX(), primaryVertex.getY(), primaryVertex.getZ()}, trackParCovDaug, 2.f, static_cast<o2::base::Propagator::MatCorrType>(cfgMaterialCorrection.value), &dcaInfoDaug);
    hDCADaugToPV->Fill(std::abs(dcaInfoDaug[0]), chargeCombSvCand);
    if (std::abs(dcaInfoDaug[0]) < minDCADaugToPV) {
      return;
    }
    hSelKinkedTrackQA->Fill(4.f, chargeCombSvCand); // MinDCADaugToPV cut bin

    if (updateMothTrackUsePV) {
      // update the mother track parameters using the primary vertex
      trackParCovMoth = trackParCovMothPV;
      if (!trackParCovMoth.update(primaryVertex)) {
        return;
      }
    }
    hSelKinkedTrackQA->Fill(5.f, chargeCombSvCand); // UpdateMothTrackUsePV cut bin

    int nCand = 0;
    try {
      nCand = fitter.process(trackParCovMoth, trackParCovDaug);
    } catch (...) {
      LOG(error) << "Exception caught in DCA fitter process call!";
      return;
    }
    if (nCand == 0) {
      return;
    }
    hSelKinkedTrackQA->Fill(6.f, chargeCombSvCand); // DCAFitter cut bin

    if (!fitter.propagateTracksToVertex()) {
      return;
    }
    hSelKinkedTrackQA->Fill(7.f, chargeCombSvCand); // PropagateTracksToVertex cut bin

    auto propMothTrack = fitter.getTrack(0);
    auto propDaugTrack = fitter.getTrack(1);
    kinkCand.decVtx = fitter.getPCACandidatePos();

    // cut on decay radius to 17 cm
    float decRad2 = kinkCand.decVtx[0] * kinkCand.decVtx[0] + kinkCand.decVtx[1] * kinkCand.decVtx[1];
    hMothDecRad2->Fill(decRad2, chargeCombSvCand);
    int idxFirstOuterBarrelLayer = nItsInnerBarrelLayers;
    if (doSVRadiusCut && decRad2 < LayerRadii[idxFirstOuterBarrelLayer] * LayerRadii[idxFirstOuterBarrelLayer]) {
      return;
    }
    hSelKinkedTrackQA->Fill(8.f, chargeCombSvCand); // MothDecRad2InsideL4 cut bin

    // get last layer hitted by the mother and the first layer hitted by the daughter
    int lastLayerMoth = 0, firstLayerDaug = 0;
    for (int i = 0; i < nItsTotalLayers; i++) {
      if (trackMoth.itsClusterMap() & (1 << i)) {
        lastLayerMoth = i;
      }
    }

    for (int i = 0; i < nItsTotalLayers; i++) {
      if (trackDaug.itsClusterMap() & (1 << i)) {
        firstLayerDaug = i;
        break;
      }
    }

    if (doSVRadiusCut && lastLayerMoth >= firstLayerDaug) {
      return;
    }
    hSelKinkedTrackQA->Fill(9.f, chargeCombSvCand); // MothDaugLastLayers cut bin

    if (doSVRadiusCut && decRad2 < LayerRadii[lastLayerMoth] * LayerRadii[lastLayerMoth]) {
      return;
    }
    hSelKinkedTrackQA->Fill(10.f, chargeCombSvCand); // MothLastLayerRadiusCheck cut bin

    for (int i = 0; i < nCoords; i++) {
      kinkCand.decVtx[i] -= kinkCand.primVtx[i];
    }

    propMothTrack.getPxPyPzGlo(kinkCand.momMoth);
    propDaugTrack.getPxPyPzGlo(kinkCand.momDaug);
    for (int i = 0; i < nCoords; i++) {
      kinkCand.momMoth[i] *= charge;
      kinkCand.momDaug[i] *= charge;
    }
    float pMoth = propMothTrack.getP() * charge;
    float pDaug = propDaugTrack.getP() * charge;
    float spKink = kinkCand.momMoth[0] * kinkCand.momDaug[0] + kinkCand.momMoth[1] * kinkCand.momDaug[1] + kinkCand.momMoth[2] * kinkCand.momDaug[2];
    kinkCand.kinkAngle = std::acos(spKink / (pMoth * pDaug));

    std::array<float, 3> neutDauMom{0.f, 0.f, 0.f};
    for (int i = 0; i < nCoords; i++) {
      neutDauMom[i] = kinkCand.momMoth[i] - kinkCand.momDaug[i];
    }

    float chargedDauE = std::sqrt(pDaug * pDaug + chargedDauMass * chargedDauMass);
    float neutE = std::sqrt(neutDauMom[0] * neutDauMom[0] + neutDauMom[1] * neutDauMom[1] + neutDauMom[2] * neutDauMom[2] + neutDauMass * neutDauMass);
    float invMass = std::sqrt((chargedDauE + neutE) * (chargedDauE + neutE) - (pMoth * pMoth));

    h2DeDxDaugSel->Fill(trackDaug.tpcInnerParam() * trackDaug.sign(), trackDaug.tpcSignal());
    h2KinkAnglePt->Fill(trackMoth.pt() * charge * trackMoth.sign(), kinkCand.kinkAngle * radToDeg);
    h2MothMassPt->Fill(trackMoth.pt() * charge * trackMoth.sign(), invMass);
    h2ClsMapPtMoth->Fill(trackMoth.pt() * charge * trackMoth.sign(), trackMoth.itsClusterMap());
    h2ClsMapPtDaug->Fill(trackDaug.pt() * charge * trackDaug.sign(), trackDaug.itsClusterMap());
    hMothDaughSignsFinal->Fill(trackMoth.sign(), trackDaug.sign());

    kinkCand.collisionID = collision.globalIndex();
    kinkCand.mothTrackID = trackMoth.globalIndex();
    kinkCand.daugTrackID = trackDaug.globalIndex();

    kinkCand.dcaXYmoth = dcaInfoMoth[0];
    kinkCand.mothSign = trackMoth.sign();
    kinkCand.dcaXYdaug = dcaInfoDaug[0];
    kinkCand.dcaKinkTopo = std::sqrt(fitter.getChi2AtPCACandidate());
    kinkCandidates.push_back(kinkCand);
    isAccepted = true;
    return;
  }

  void initCCDB(aod::BCs::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();
    LOG(info) << "Initializing CCDB for run " << mRunNumber;
    o2::parameters::GRPMagField* grpmag = ccdb->getForRun<o2::parameters::GRPMagField>(grpmagPath, mRunNumber);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    mBz = grpmag->getNominalL3Field();
    fitter.setBz(mBz);

    if (!lut) {
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
      int mat{static_cast<int>(cfgMaterialCorrection)};
      fitter.setMatCorrType(static_cast<o2::base::Propagator::MatCorrType>(mat));
    }
    o2::base::Propagator::Instance()->setMatLUT(lut);
    LOG(info) << "Task initialized for run " << mRunNumber << " with magnetic field " << mBz << " kZG";
  }

  template <typename TColls, typename TTracks, typename TAmbiTracks>
  void buildSvPool(const TColls& collisions, const TTracks& tracks, const TAmbiTracks& ambiguousTracks, const aod::BCs& bcs)
  {
    svCreator.clearPools();
    svCreator.fillBC2Coll(collisions, bcs);

    for (const auto& track : tracks) {
      if (std::abs(track.eta()) > etaMax) {
        continue;
      }

      bool isDaug = selectDaugTrack(track);
      bool isMoth = selectMothTrack(track);
      if (!isDaug && !isMoth) {
        continue;
      }
      int pdgHypo = isMoth ? 1 : 0;
      svCreator.appendTrackCand(track, collisions, pdgHypo, ambiguousTracks, bcs);
    }
  }

  template <typename TColls, typename TTracks>
  void fillOutputsData(const TColls& collisions, const TTracks& tracks, const aod::AmbiguousTracks& ambiTracks, const aod::BCs& bcs)
  {
    kinkCandidates.clear();

    buildSvPool(collisions, tracks, ambiTracks, bcs);
    auto& kinkPool = svCreator.getSVCandPool(collisions, !unlikeSignBkg);
    bool isAccepted = false;
    for (const auto& svCand : kinkPool) {
      buildKinkCand(svCand, tracks, isAccepted, collisions);
    }

    // sort kinkCandidates by collisionID to allow joining with collision table
    std::sort(kinkCandidates.begin(), kinkCandidates.end(), [](const kinkCandidate& a, const kinkCandidate& b) { return a.collisionID < b.collisionID; });

    for (const auto& kinkCand : kinkCandidates) {
      if (fillDebugTable) {
        outputDataTableUB(kinkCand.decVtx[0], kinkCand.decVtx[1], kinkCand.decVtx[2],
                          kinkCand.mothSign, kinkCand.momMoth[0], kinkCand.momMoth[1], kinkCand.momMoth[2],
                          kinkCand.momDaug[0], kinkCand.momDaug[1], kinkCand.momDaug[2],
                          kinkCand.dcaXYmoth, kinkCand.dcaXYdaug, kinkCand.dcaKinkTopo);
      } else {
        outputDataTable(kinkCand.collisionID, kinkCand.mothTrackID, kinkCand.daugTrackID,
                        kinkCand.decVtx[0], kinkCand.decVtx[1], kinkCand.decVtx[2],
                        kinkCand.mothSign, kinkCand.momMoth[0], kinkCand.momMoth[1], kinkCand.momMoth[2],
                        kinkCand.momDaug[0], kinkCand.momDaug[1], kinkCand.momDaug[2],
                        kinkCand.dcaXYmoth, kinkCand.dcaXYdaug, kinkCand.dcaKinkTopo);
      }
    }
  }

  void processData(aod::Collisions const& collisions, TracksFull const& tracks, aod::AmbiguousTracks const& ambiTracks, aod::BCs const& bcs)
  {
    fillOutputsData(collisions, tracks, ambiTracks, bcs);
  }
  PROCESS_SWITCH(kinkBuilder, processData, "Data processing", true);

  void processDataWCentSel(CollisionsCentSel const& collisions, TracksFull const& tracks, aod::AmbiguousTracks const& ambiTracks, aod::BCs const& bcs)
  {
    fillOutputsData(collisions, tracks, ambiTracks, bcs);
  }
  PROCESS_SWITCH(kinkBuilder, processDataWCentSel, "Data processing with centrality selection", false);

  template <bool checkKinkDaugPdg, typename TMother>
  int matchKinkDecay(const TMother& motherPart, const aod::McParticles& mcParticles)
  {
    int pdgMother = motherPart.pdgCode();
    int8_t sign = 0;
    int pdgCodeNeutralDaug{-1}, pdgCodeChargedDaug{-1};
    std::array<int, 2> finState = {-1, -1};
    switch (std::abs(pdgMother)) {
      case PDG_t::kSigmaMinus: {
        // Swap the sign of the neutral decay products in case of anti-particles
        pdgCodeNeutralDaug = (pdgMother > 0) ? +PDG_t::kNeutron : -PDG_t::kNeutron;
        pdgCodeChargedDaug = (pdgMother > 0) ? +PDG_t::kPiMinus : +PDG_t::kPiPlus;
        finState = {pdgCodeChargedDaug, pdgCodeNeutralDaug}; // Both decay channels have the same neutral daughter
        if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, motherPart, pdgMother, finState, true, &sign, DepthMcMatchMax)) {
          return SigmaMinusToPiMinusNeutron;
        }
        break;
      }
      case PDG_t::kSigmaPlus: {
        // Swap the sign of the neutral decay products in case of anti-particles
        pdgCodeNeutralDaug = (pdgMother > 0) ? +PDG_t::kNeutron : -PDG_t::kNeutron;
        pdgCodeChargedDaug = (pdgMother > 0) ? +PDG_t::kPiPlus : +PDG_t::kPiMinus;
        finState = {pdgCodeChargedDaug, pdgCodeNeutralDaug};
        if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, motherPart, pdgMother, finState, true, &sign, DepthMcMatchMax)) {
          return SigmaPlusToPiPlusNeutron;
        }

        // Swap the sign of the neutral decay products in case of anti-particles
        pdgCodeNeutralDaug = (pdgMother > 0) ? +PDG_t::kPi0 : -PDG_t::kPi0;
        pdgCodeChargedDaug = (pdgMother > 0) ? +PDG_t::kProton : -PDG_t::kProton;
        finState = {pdgCodeChargedDaug, pdgCodeNeutralDaug};
        if (RecoDecay::isMatchedMCGen<false, true>(mcParticles, motherPart, pdgMother, finState, true, &sign, DepthMcMatchMax)) {
          return SigmaPlusToProtonPi0;
        }
        break;
      }
      default:
        LOG(warning) << "No matching function implemented for the selected mother particle hypothesis. Returning -1.";
        return -1;
    }

    return -1;
  }

  template <typename TColls, typename TTracks>
  void fillOutputsMc(const TColls& mcRecoCollisions, const TTracks& tracksMc, const aod::AmbiguousTracks& ambiTracksMc, const aod::BCs& bcs, const aod::McParticles& mcParticles)
  {
    kinkCandidates.clear();

    buildSvPool(mcRecoCollisions, tracksMc, ambiTracksMc, bcs);
    auto& kinkPool = svCreator.getSVCandPool(mcRecoCollisions, !unlikeSignBkg);
    for (const auto& svCand : kinkPool) {
      // Perform matching of the kink candidate
      auto trackMoth = tracksMc.rawIteratorAt(svCand.tr0Idx);
      auto genMothPart = trackMoth.template mcParticle_as<aod::McParticles>();
      auto trackDaug = tracksMc.rawIteratorAt(svCand.tr1Idx);
      auto genDaugPart = trackDaug.template mcParticle_as<aod::McParticles>();
      int decayChannel{-1};
      if (skipBkgCands) {

        // Check mother PDG first
        if (std::find(mothMatchPdgCodes.begin(), mothMatchPdgCodes.end(), std::abs(genMothPart.pdgCode())) == mothMatchPdgCodes.end()) {
          continue;
        }

        // Check that the daughter was generated from the selected mother
        int genDaugMothIdx = RecoDecay::getMother(mcParticles, genDaugPart, genMothPart.pdgCode(), true);
        if (genDaugMothIdx != trackMoth.mcParticleId()) {
          continue; // Skip candidates where the daughter is not coming from the selected mother
        }
        decayChannel = matchKinkDecay<true>(genMothPart, mcParticles);
        if (decayChannel < 0) {
          continue; // Skip candidates that do not match the decay channels of interest
        }
        hRecCandidates->Fill(decayChannel, trackMoth.pt()); // Decay channel match bin
      }
      bool isAccepted = false;
      buildKinkCand(svCand, tracksMc, isAccepted, mcRecoCollisions);
      // If candidate passed all selections and was built, fill histogram
      if (isAccepted && skipBkgCands) {
        hRecPtKinkAngle[decayChannel]->Fill(trackMoth.pt(), std::abs(trackMoth.phi() - trackDaug.phi()) * radToDeg);
      }
    }

    // sort kinkCandidates by collisionID to allow joining with collision table
    std::sort(kinkCandidates.begin(), kinkCandidates.end(), [](const kinkCandidate& a, const kinkCandidate& b) { return a.collisionID < b.collisionID; });

    for (const auto& kinkCand : kinkCandidates) {
      if (fillDebugTable) {
        outputDataTableUB(kinkCand.decVtx[0], kinkCand.decVtx[1], kinkCand.decVtx[2],
                          kinkCand.mothSign, kinkCand.momMoth[0], kinkCand.momMoth[1], kinkCand.momMoth[2],
                          kinkCand.momDaug[0], kinkCand.momDaug[1], kinkCand.momDaug[2],
                          kinkCand.dcaXYmoth, kinkCand.dcaXYdaug, kinkCand.dcaKinkTopo);
      } else {
        outputDataTable(kinkCand.collisionID, kinkCand.mothTrackID, kinkCand.daugTrackID,
                        kinkCand.decVtx[0], kinkCand.decVtx[1], kinkCand.decVtx[2],
                        kinkCand.mothSign, kinkCand.momMoth[0], kinkCand.momMoth[1], kinkCand.momMoth[2],
                        kinkCand.momDaug[0], kinkCand.momDaug[1], kinkCand.momDaug[2],
                        kinkCand.dcaXYmoth, kinkCand.dcaXYdaug, kinkCand.dcaKinkTopo);
      }
    }

    // Generated kinks
    for (const auto& mcPart : mcParticles) {
      if (std::find(mothMatchPdgCodes.begin(), mothMatchPdgCodes.end(), std::abs(mcPart.pdgCode())) == mothMatchPdgCodes.end()) {
        continue; // Skip if mother particle does not match the PDG codes of interest
      }
      int decayChannel = matchKinkDecay<false>(mcPart, mcParticles); // Dummy daughter PDG code since we are not checking it in this case
      if (decayChannel < 0) {
        continue; // Skip if no matching decay channel is found
      }
      hGenCandidates->Fill(decayChannel, mcPart.pt());
      std::vector<int> arrDaughIdxs = {};
      RecoDecay::getDaughters<false, true>(mcPart, &arrDaughIdxs, std::array{0}, DepthMcMatchMax);
      for (auto iProng = 0u; iProng < arrDaughIdxs.size(); ++iProng) {
        auto daughI = mcParticles.rawIteratorAt(arrDaughIdxs[iProng]);
        if (std::abs(daughI.pdgCode()) == PDG_t::kPiPlus || std::abs(daughI.pdgCode()) == PDG_t::kProton) {
          hGenPtKinkAngle[decayChannel]->Fill(mcPart.pt(), std::abs(mcPart.phi() - daughI.phi()) * radToDeg);
        }
      }
    }
  }

  void processMc(McRecoCollisions const& mcRecoCollisions,
                 aod::McParticles const& mcParticles,
                 TracksFullMc const& tracksMc,
                 aod::AmbiguousTracks const& ambiTracksMc,
                 aod::BCs const& bcs)
  {
    fillOutputsMc(mcRecoCollisions, tracksMc, ambiTracksMc, bcs, mcParticles);
  }
  PROCESS_SWITCH(kinkBuilder, processMc, "MC processing", false);

  void processMcWCent(McRecoCollisionsCentSel const& mcRecoCollisions,
                      aod::McParticles const& mcParticles,
                      TracksFullMc const& tracksMc,
                      aod::AmbiguousTracks const& ambiTracksMc,
                      aod::BCs const& bcs)
  {
    fillOutputsMc(mcRecoCollisions, tracksMc, ambiTracksMc, bcs, mcParticles);
  }
  PROCESS_SWITCH(kinkBuilder, processMcWCent, "MC processing with centrality selection", false);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<kinkBuilder>(cfgc)};
}
