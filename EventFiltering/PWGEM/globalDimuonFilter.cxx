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

// \brief software trigger for global dimuons
// \author daiki.sekihata@cern.ch

#include "Common/CCDB/EventSelectionParams.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/fwdtrackUtilities.h"
#include "Common/DataModel/CollisionAssociationTables.h"
#include "Common/DataModel/EventSelection.h"
#include "EventFiltering/filterTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Field/MagFieldParam.h>
#include <Field/MagneticField.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <MCHTracking/TrackExtrap.h>
#include <ReconstructionDataFormats/GlobalFwdTrack.h>

#include <Math/Vector4D.h> // IWYU pragma: keep (do not replace with Math/Vector4Dfwd.h)
#include <Math/Vector4Dfwd.h>
#include <TGeoGlobalMagField.h>
#include <TH1.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

struct globalDimuonFilter {
  o2::framework::Produces<o2::aod::GlobalDimuonFilters> tags;

  // Configurables
  o2::framework::Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  o2::framework::Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  o2::framework::Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  o2::framework::Configurable<int> minNmuon{"minNmuon", 2, "min number of global muon candidates per collision"};

  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "eventCutGroup";
    o2::framework::Configurable<uint16_t> cfgNcontribMin{"cfgNcontribMin", 0, "min. Ncontributors to PV"};
    o2::framework::Configurable<float> cfgChi2PerNcontribMax{"cfgChi2PerNcontribMax", 1e+10, "max. chi2/Ncontrib of PV"};
    o2::framework::Configurable<float> minZvtx{"minZvtx", -11.f, "min. Zvtx of collision"};
    o2::framework::Configurable<float> maxZvtx{"maxZvtx", +11.f, "max. Zvtx of collision"};
    o2::framework::Configurable<bool> cfgRequireFT0AND{"cfgRequireFT0AND", true, "require FT0AND"};
    o2::framework::Configurable<bool> cfgRequireNoTFB{"cfgRequireNoTFB", true, "require No time frame border"};
    o2::framework::Configurable<bool> cfgRequireNoITSROFB{"cfgRequireNoITSROFB", false, "require no ITS readout frame border"};
    o2::framework::Configurable<bool> cfgRequireVertexITSTPC{"cfgRequireVertexITSTPC", false, "require Vertex ITSTPC in event cut"};             // ITS-TPC matched track contributes PV.
    o2::framework::Configurable<bool> cfgRequireVertexTOFmatched{"cfgRequireVertexTOFmatched", false, "require Vertex TOFmatched in event cut"}; // ITS-TPC-TOF matched track contributes PV.
    o2::framework::Configurable<bool> cfgRequireNoSameBunchPileup{"cfgRequireNoSameBunchPileup", false, "require no same bunch pileup in event cut"};
    o2::framework::Configurable<bool> cfgRequireGoodZvtxFT0vsPV{"cfgRequireGoodZvtxFT0vsPV", false, "require good Zvtx between FT0 vs. PV in event cut"};
  } eventCutGroup;

  struct : o2::framework::ConfigurableGroup {
    std::string prefix = "glMuonCutGroup";
    o2::framework::Configurable<float> minPt{"minPt", 0.55, "min pt for muon"};
    o2::framework::Configurable<float> maxPt{"maxPt", 1e+10, "max pt for muon"};
    o2::framework::Configurable<float> minEta{"minEta", -3.65, "min. eta acceptance for MFT-MCH-MID"};
    o2::framework::Configurable<float> maxEta{"maxEta", -2.45, "max. eta acceptance for MFT-MCH-MID"};
    o2::framework::Configurable<float> minRabs{"minRabs", 17.6, "min. R at absorber end for global muon"}; // std::tan(2.f * std::atan(std::exp(- -3.6)) ) * -505. = 27.6
    o2::framework::Configurable<float> maxRabs{"maxRabs", 89.5, "max. R at absorber end"};
    o2::framework::Configurable<float> maxDCAxy{"maxDCAxy", 0.5, "max. DCAxy for global muons"};
    o2::framework::Configurable<float> maxMatchingChi2MCHMFT{"maxMatchingChi2MCHMFT", 100.f, "max. chi2 for MCH-MFT matching"};
    o2::framework::Configurable<float> maxChi2{"maxChi2", 1e+10, "max. chi2/ndf for global muon"};
    o2::framework::Configurable<float> maxChi2MFT{"maxChi2MFT", 1e+10, "max. chi2/ndf for MFTsa"};
    o2::framework::Configurable<int> minNclsMFT{"minNclsMFT", 5, "min ncluster of MFT"};
    o2::framework::Configurable<int> minNclsMCH{"minNclsMCH", 5, "min ncluster of MCH"};
    o2::framework::Configurable<bool> refitGlobalMuon{"refitGlobalMuon", true, "flag to refit global muon"};
    o2::framework::Configurable<float> matchingZ{"matchingZ", -77.5, "z position where matching is performed"};
    o2::framework::Configurable<float> maxDEta{"maxDEta", 1e+10, "max. deta between MFT-MCH-MID and MCH-MID"};
    o2::framework::Configurable<float> maxDPhi{"maxDPhi", 1e+10, "max. dphi between MFT-MCH-MID and MCH-MID"};
  } glMuonCutGroup;

  struct : o2::framework::ConfigurableGroup { // tight cut
    std::string prefix = "tagMuonCutGroup";
    o2::framework::Configurable<float> minPt{"minPt", 0.8, "min pt for muon"};
    o2::framework::Configurable<float> maxPt{"maxPt", 1e+10, "max pt for muon"};
    o2::framework::Configurable<float> minEta{"minEta", -3.6, "min. eta acceptance for MFT-MCH-MID"};
    o2::framework::Configurable<float> maxEta{"maxEta", -2.5, "max. eta acceptance for MFT-MCH-MID"};
    o2::framework::Configurable<float> minRabs{"minRabs", 27.6, "min. R at absorber end for global muon (min. eta = -3.6)"}; // std::tan(2.f * std::atan(std::exp(- -3.6)) ) * -505. = 27.6
    o2::framework::Configurable<float> maxRabs{"maxRabs", 89.5, "max. R at absorber end"};
    o2::framework::Configurable<float> maxDCAxy{"maxDCAxy", 0.06, "max. DCAxy for global muons"};
    o2::framework::Configurable<float> maxMatchingChi2MCHMFT{"maxMatchingChi2MCHMFT", 40.f, "max. chi2 for MCH-MFT matching"};
    o2::framework::Configurable<float> maxChi2{"maxChi2", 4.f, "max. chi2/ndf for global muon"};
    o2::framework::Configurable<float> maxChi2MFT{"maxChi2MFT", 4.f, "max. chi2/ndf for MFTsa"};
    o2::framework::Configurable<int> minNclsMFT{"minNclsMFT", 5, "min ncluster of MFT"};
    o2::framework::Configurable<int> minNclsMCH{"minNclsMCH", 5, "min ncluster of MCH"};
    o2::framework::Configurable<bool> refitGlobalMuon{"refitGlobalMuon", true, "flag to refit global muon"};
    o2::framework::Configurable<float> matchingZ{"matchingZ", -77.5, "z position where matching is performed"};
    o2::framework::Configurable<float> maxDEta{"maxDEta", 0.08, "max. deta between MFT-MCH-MID and MCH-MID"};
    o2::framework::Configurable<float> maxDPhi{"maxDPhi", 0.08, "max. dphi between MFT-MCH-MID and MCH-MID"};
  } tagMuonCutGroup;

  struct : o2::framework::ConfigurableGroup { // loose cut
    std::string prefix = "probeMuonCutGroup";
    o2::framework::Configurable<float> minPt{"minPt", 0.01, "min pt for muon"};
    o2::framework::Configurable<float> maxPt{"maxPt", 1e+10, "max pt for muon"};
    o2::framework::Configurable<float> minEta{"minEta", -3.6, "min. eta acceptance for MFT-MCH-MID"};
    o2::framework::Configurable<float> maxEta{"maxEta", -2.5, "max. eta acceptance for MFT-MCH-MID"};
    o2::framework::Configurable<float> minRabs{"minRabs", 17.6, "min. R at absorber end for global muon (min. eta = -3.6)"}; // std::tan(2.f * std::atan(std::exp(- -3.6)) ) * -505. = 27.6
    o2::framework::Configurable<float> maxRabs{"maxRabs", 89.5, "max. R at absorber end"};
    o2::framework::Configurable<float> maxDCAxy{"maxDCAxy", 1e+10, "max. DCAxy for global muons"};
    o2::framework::Configurable<float> maxMatchingChi2MCHMFT{"maxMatchingChi2MCHMFT", 1e+10, "max. chi2 for MCH-MFT matching"};
    o2::framework::Configurable<float> maxChi2{"maxChi2", 1e+10, "max. chi2/ndf for global muon"};
    o2::framework::Configurable<float> maxChi2MFT{"maxChi2MFT", 1e+10, "max. chi2/ndf for MFTsa"};
    o2::framework::Configurable<int> minNclsMFT{"minNclsMFT", 5, "min ncluster of MFT"};
    o2::framework::Configurable<int> minNclsMCH{"minNclsMCH", 5, "min ncluster of MCH"};
    o2::framework::Configurable<bool> refitGlobalMuon{"refitGlobalMuon", true, "flag to refit global muon"};
    o2::framework::Configurable<float> matchingZ{"matchingZ", -77.5, "z position where matching is performed"};
    o2::framework::Configurable<float> maxDEta{"maxDEta", 1e+10, "max. deta between MFT-MCH-MID and MCH-MID"};
    o2::framework::Configurable<float> maxDPhi{"maxDPhi", 1e+10, "max. dphi between MFT-MCH-MID and MCH-MID"};
  } probeMuonCutGroup;

  // for z shift for propagation
  o2::framework::Configurable<bool> cfgApplyZShiftFromCCDB{"cfgApplyZShiftFromCCDB", false, "flag to apply z shift"};
  o2::framework::Configurable<std::string> cfgZShiftPath{"cfgZShiftPath", "Users/m/mcoquet/ZShift", "CCDB path for z shift to apply to forward tracks"};
  o2::framework::Configurable<float> cfgManualXShiftMFTtop{"cfgManualXShiftMFTtop", 0, "manual x shift for propagation of global muon to PV"};
  o2::framework::Configurable<float> cfgManualYShiftMFTtop{"cfgManualYShiftMFTtop", 0, "manual y shift for propagation of global muon to PV"};
  o2::framework::Configurable<float> cfgManualZShiftMFTtop{"cfgManualZShiftMFTtop", 0, "manual z shift for propagation of global muon to PV"};
  o2::framework::Configurable<float> cfgManualXShiftMFTbottom{"cfgManualXShiftMFTbottom", 0, "manual x shift for propagation of global muon to PV"};
  o2::framework::Configurable<float> cfgManualYShiftMFTbottom{"cfgManualYShiftMFTbottom", 0, "manual y shift for propagation of global muon to PV"};
  o2::framework::Configurable<float> cfgManualZShiftMFTbottom{"cfgManualZShiftMFTbottom", 0, "manual z shift for propagation of global muon to PV"};

  o2::framework::HistogramRegistry fRegistry{"output", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject, false, false};
  o2::ccdb::CcdbApi ccdbApi;
  o2::framework::Service<o2::ccdb::BasicCCDBManager> ccdb;
  int mRunNumber = 0;
  float mBz = 0;
  float mXShiftMFTtop = 0;
  float mYShiftMFTtop = 0;
  float mZShiftMFTtop = 0;
  float mXShiftMFTbottom = 0;
  float mYShiftMFTbottom = 0;
  float mZShiftMFTbottom = 0;

  void init(o2::framework::InitContext&)
  {
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);
    ccdbApi.init(ccdburl);
    mRunNumber = 0;
    mBz = 0;
    mXShiftMFTtop = 0;
    mYShiftMFTtop = 0;
    mZShiftMFTtop = 0;
    mXShiftMFTbottom = 0;
    mYShiftMFTbottom = 0;
    mZShiftMFTbottom = 0;

    addHistograms();
  }

  template <typename TBC>
  void initCCDB(TBC const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }
    mRunNumber = bc.runNumber();

    std::map<std::string, std::string> metadata;
    auto soreor = o2::ccdb::BasicCCDBManager::getRunDuration(ccdbApi, mRunNumber);
    auto ts = soreor.first;
    auto grpmag = ccdbApi.retrieveFromTFileAny<o2::parameters::GRPMagField>(grpmagPath, metadata, ts);
    o2::base::Propagator::initFieldFromGRP(grpmag);
    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      ccdb->get<TGeoManager>(geoPath);
    }
    o2::mch::TrackExtrap::setField();
    const double centerMFT[3] = {0, 0, -61.4};
    o2::field::MagneticField* field = dynamic_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());
    mBz = field->getBz(centerMFT); // Get field at centre of MFT
    LOGF(info, "Bz at center of MFT = %f kZG", mBz);

    if (cfgApplyZShiftFromCCDB) {
      auto* zShift = ccdb->getForTimeStamp<std::vector<float>>(cfgZShiftPath, bc.timestamp());
      if (zShift != nullptr && !zShift->empty()) {
        LOGF(info, "reading z shift %f from %s", (*zShift)[0], cfgZShiftPath.value);
        mZShiftMFTtop = (*zShift)[0];
        mZShiftMFTbottom = (*zShift)[0];
      } else {
        LOGF(info, "z shift is not found in ccdb path %s. set to 0 cm", cfgZShiftPath.value);
        mZShiftMFTtop = 0;
        mZShiftMFTbottom = 0;
      }
    } else {
      LOGF(info, "X shift for MFT top is manually set to %f cm", cfgManualXShiftMFTtop.value);
      LOGF(info, "X shift for MFT bottom is manually set to %f cm", cfgManualXShiftMFTbottom.value);
      LOGF(info, "Y shift for MFT top is manually set to %f cm", cfgManualYShiftMFTtop.value);
      LOGF(info, "Y shift for MFT bottom is manually set to %f cm", cfgManualYShiftMFTbottom.value);
      LOGF(info, "Z shift for MFT top is manually set to %f cm", cfgManualZShiftMFTtop.value);
      LOGF(info, "Z shift for MFT bottom is manually set to %f cm", cfgManualZShiftMFTbottom.value);
      mXShiftMFTtop = cfgManualXShiftMFTtop;
      mYShiftMFTtop = cfgManualYShiftMFTtop;
      mZShiftMFTtop = cfgManualZShiftMFTtop;
      mXShiftMFTbottom = cfgManualXShiftMFTbottom;
      mYShiftMFTbottom = cfgManualYShiftMFTbottom;
      mZShiftMFTbottom = cfgManualZShiftMFTbottom;
    }
  }

  void addHistograms()
  {
    auto hCollisionCounter{std::get<std::shared_ptr<TH1>>(fRegistry.add("hCollisionCounter", "Number of collisions", o2::framework::HistType::kTH1D, {{10, 0.5, 10.5}}))};
    hCollisionCounter->GetXaxis()->SetBinLabel(1, "all");
    hCollisionCounter->GetXaxis()->SetBinLabel(2, "TVX");
    hCollisionCounter->GetXaxis()->SetBinLabel(3, "No TFB");
    hCollisionCounter->GetXaxis()->SetBinLabel(4, "No ITS ROFB");
    hCollisionCounter->GetXaxis()->SetBinLabel(5, "|Z_{vtx}| < 10 cm");
    hCollisionCounter->GetXaxis()->SetBinLabel(6, "sel8");
    hCollisionCounter->GetXaxis()->SetBinLabel(7, "MB");
    hCollisionCounter->GetXaxis()->SetBinLabel(8, "MB && global dimuon");
    hCollisionCounter->GetXaxis()->SetBinLabel(9, "MB && global dimuon for QA");

    fRegistry.add("hNmu", "#mu multiplicity;N_{#mu}^{global} candidates per collision", o2::framework::HistType::kTH1D, {{11, -0.5, 10.5}}, false);

    fRegistry.add("Vertex/MB/hZvtx", "vertex z; Z_{vtx} (cm)", o2::framework::HistType::kTH1D, {{40, -20, +20}}, false);
    fRegistry.add("Vertex/MB/hChi2vsNContrib", "vertex #chi^{2}/N_{contrib} vs. N_{contrib};N_{contrib};#chi^{2}/N_{contrib}", o2::framework::HistType::kTH2D, {{200, 0.5, 200.5}, {100, 0, 10}}, false);
    fRegistry.add("Vertex/MB/hSigmaX", "vertex #sigma_{X} vs. N_{contrib};N_{contrib};#sigma_{X} (#mum)", o2::framework::HistType::kTH2D, {{200, 0.5, 200.5}, {1000, 0, 100}}, false);
    fRegistry.add("Vertex/MB/hSigmaY", "vertex #sigma_{Y} vs. N_{contrib};N_{contrib};#sigma_{Y} (#mum)", o2::framework::HistType::kTH2D, {{200, 0.5, 200.5}, {1000, 0, 100}}, false);
    fRegistry.add("Vertex/MB/hSigmaZ", "vertex #sigma_{Z} vs. N_{contrib};N_{contrib};#sigma_{Z} (#mum)", o2::framework::HistType::kTH2D, {{200, 0.5, 200.5}, {1000, 0, 100}}, false);
    fRegistry.addClone("Vertex/MB/", "Vertex/Triggered/");

    fRegistry.add("MFTMCHMID/positive/hPt", "pT;p_{T} (GeV/c)", o2::framework::HistType::kTH1D, {{200, 0.0f, 10}}, false);
    fRegistry.add("MFTMCHMID/positive/hPtResolution", "pT;p_{T} (GeV/c);#Deltap_{T}/p_{T}", o2::framework::HistType::kTH2D, {{200, 0.0f, 10}, {200, 0, 0.1}}, false);
    fRegistry.add("MFTMCHMID/positive/hEtaPhi", "#eta vs. #varphi;#varphi (rad.);#eta", o2::framework::HistType::kTH2D, {{180, 0, 2 * M_PI}, {200, -4.f, -2.f}}, false);
    fRegistry.add("MFTMCHMID/positive/hEtaPhi_MatchedMCHMID", "#eta vs. #varphi;#varphi (rad.);#eta", o2::framework::HistType::kTH2D, {{180, 0, 2 * M_PI}, {200, -4.f, -2.f}}, false);
    fRegistry.add("MFTMCHMID/positive/hDEtaDPhi", "#Delta#eta vs. #Delta#varphi;#Delta#varphi (rad.);#Delta#eta", o2::framework::HistType::kTH2D, {{90, -M_PI / 4, M_PI / 4}, {100, -0.5, +0.5}}, false);
    fRegistry.add("MFTMCHMID/positive/hSign", "sign;sign", o2::framework::HistType::kTH1D, {{3, -1.5, +1.5}}, false);
    fRegistry.add("MFTMCHMID/positive/hNclusters", "Nclusters;Nclusters", o2::framework::HistType::kTH1D, {{21, -0.5f, 20.5}}, false);
    fRegistry.add("MFTMCHMID/positive/hNclustersMFT", "NclustersMFT;Nclusters MFT", o2::framework::HistType::kTH1D, {{11, -0.5f, 10.5}}, false);
    fRegistry.add("MFTMCHMID/positive/hRatAbsorberEnd", "R at absorber end;R at absorber end (cm)", o2::framework::HistType::kTH1D, {{200, 0, 100}}, false);
    fRegistry.add("MFTMCHMID/positive/hPDCA_Rabs", "pDCA vs. Rabs;R at absorber end (cm);p #times DCA (GeV/c #upoint cm)", o2::framework::HistType::kTH2D, {{200, 0, 100}, {100, 0, 1000}}, false);
    fRegistry.add("MFTMCHMID/positive/hChi2_Pt", "chi2;p_{T,#mu} (GeV/c);#chi^{2}/ndf", o2::framework::HistType::kTH2D, {{200, 0, 10}, {100, 0.0f, 10}}, false);
    fRegistry.add("MFTMCHMID/positive/hChi2MFT_Pt", "chi2 MFT;p_{T,#mu} (GeV/c);#chi^{2} MFT/ndf", o2::framework::HistType::kTH2D, {{200, 0, 10}, {100, 0.0f, 10}}, false);
    fRegistry.add("MFTMCHMID/positive/hChi2MatchMCHMID_Pt", "chi2 match MCH-MID;p_{T,#mu} (GeV/c);matching #chi^{2}/ndf between MCH-MID", o2::framework::HistType::kTH2D, {{200, 0, 10}, {200, 0.0f, 20}}, false);
    fRegistry.add("MFTMCHMID/positive/hChi2MatchMCHMFT_Pt", "chi2 match MCH-MFT;p_{T,#mu} (GeV/c);matching #chi^{2}/ndf between MFT-MCH", o2::framework::HistType::kTH2D, {{200, 0, 10}, {100, 0.0f, 100}}, false);
    fRegistry.add("MFTMCHMID/positive/hDCAxy2D", "DCA x vs. y;DCA_{x} (cm);DCA_{y} (cm)", o2::framework::HistType::kTH2D, {{400, -1, 1}, {400, -1, +1}}, false);
    fRegistry.add("MFTMCHMID/positive/hDCAxy2DinSigma", "DCA x vs. y in sigma;DCA_{x} (#sigma);DCA_{y} (#sigma)", o2::framework::HistType::kTH2D, {{200, -10, 10}, {200, -10, +10}}, false);
    fRegistry.add("MFTMCHMID/positive/hDCAxy", "DCAxy;DCA_{xy} (cm);", o2::framework::HistType::kTH1D, {{1000, 0, 1}}, false);
    fRegistry.add("MFTMCHMID/positive/hDCAxyinSigma", "DCAxy in sigma;DCA_{xy} (#sigma);", o2::framework::HistType::kTH1D, {{100, 0, 10}}, false);
    fRegistry.add("MFTMCHMID/positive/hLog10Chi2IP", "chi2IP;log_{10}(#chi^{2}_{IP})", o2::framework::HistType::kTH1D, {{1000, -5, 5}}, false);
    fRegistry.add("MFTMCHMID/positive/hDCAxResolutionvsPt", "DCA_{x} resolution vs. p_{T};p_{T} (GeV/c);DCA_{x} resolution (#mum);", o2::framework::HistType::kTH2D, {{100, 0, 10.f}, {500, 0, 500}}, false);
    fRegistry.add("MFTMCHMID/positive/hDCAyResolutionvsPt", "DCA_{y} resolution vs. p_{T};p_{T} (GeV/c);DCA_{y} resolution (#mum);", o2::framework::HistType::kTH2D, {{100, 0, 10.f}, {500, 0, 500}}, false);
    fRegistry.add("MFTMCHMID/positive/hDCAxyResolutionvsPt", "DCA_{xy} resolution vs. p_{T};p_{T} (GeV/c);DCA_{xy} resolution (#mum);", o2::framework::HistType::kTH2D, {{100, 0, 10.f}, {500, 0, 500}}, false);
    fRegistry.add("MFTMCHMID/positive/hDCAx_Phi", "DCAx vs. #varphi;#varphi (rad.);DCA_{x} (cm)", o2::framework::HistType::kTH2D, {{180, -M_PI, M_PI}, {400, -0.2, +0.2}}, false);
    fRegistry.add("MFTMCHMID/positive/hDCAy_Phi", "DCAy vs. #varphi;#varphi (rad.);DCA_{y} (cm)", o2::framework::HistType::kTH2D, {{180, -M_PI, M_PI}, {400, -0.2, +0.2}}, false);
    fRegistry.add("MFTMCHMID/positive/hDCAx_PosZ_MFTtop", "DCAx vs. posZ for MFT top;Z_{vtx} (cm);DCA_{x} (cm)", o2::framework::HistType::kTH2D, {{200, -10, +10}, {400, -0.2, +0.2}}, false);
    fRegistry.add("MFTMCHMID/positive/hDCAy_PosZ_MFTtop", "DCAy vs. posZ for MFT top;Z_{vtx} (cm);DCA_{y} (cm)", o2::framework::HistType::kTH2D, {{200, -10, +10}, {400, -0.2, +0.2}}, false);
    fRegistry.add("MFTMCHMID/positive/hDCAx_PosZ_MFTbottom", "DCAx vs. posZ for MFT bottom;Z_{vtx} (cm);DCA_{x} (cm)", o2::framework::HistType::kTH2D, {{200, -10, +10}, {400, -0.2, +0.2}}, false);
    fRegistry.add("MFTMCHMID/positive/hDCAy_PosZ_MFTbottom", "DCAy vs. posZ for MFT bottom;Z_{vtx} (cm);DCA_{y} (cm)", o2::framework::HistType::kTH2D, {{200, -10, +10}, {400, -0.2, +0.2}}, false);
    fRegistry.addClone("MFTMCHMID/positive/", "MFTMCHMID/negative/");

    const o2::framework::AxisSpec axisMll{{0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.505, 0.51, 0.515, 0.52, 0.525, 0.53, 0.535, 0.54, 0.545, 0.55, 0.555, 0.56, 0.565, 0.57, 0.575, 0.58, 0.585, 0.59, 0.595, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.75, 2.80, 2.85, 2.90, 2.95, 3.00, 3.05, 3.10, 3.15, 3.20, 3.25, 3.30, 3.35, 3.40, 3.45, 3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85, 3.90, 3.95, 4.00, 4.50, 5.00, 5.50, 6.00, 6.50, 7.00, 7.50, 8.00, 8.1, 8.2, 8.3, 8.4, 8.50, 8.60, 8.70, 8.80, 8.90, 9.00, 9.10, 9.20, 9.30, 9.40, 9.50, 9.60, 9.70, 9.80, 9.90, 10.00, 10.10, 10.20, 10.30, 10.40, 10.50, 10.60, 10.70, 10.80, 10.90, 11.0, 11.1, 11.2, 11.3, 11.4, 11.5, 11.6, 11.7, 11.8, 11.9, 12.0}, "m_{#mu#mu} (GeV/c^{2})"};
    const o2::framework::AxisSpec axisPtll{{0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10}, "p_{T,#mu#mu} (GeV/c)"};
    const o2::framework::AxisSpec axisYll{40, -4.0, -2.0, "y_{#mu#mu}"};

    fRegistry.add("Pair/same/uls/hs", "dimuon", o2::framework::HistType::kTHnSparseD, {axisMll, axisPtll, axisYll}, false);
    fRegistry.add("Pair/same/lspp/hs", "dimuon", o2::framework::HistType::kTHnSparseD, {axisMll, axisPtll, axisYll}, false);
    fRegistry.add("Pair/same/lsmm/hs", "dimuon", o2::framework::HistType::kTHnSparseD, {axisMll, axisPtll, axisYll}, false);

    const o2::framework::AxisSpec axisMllLMR{{0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.01, 2.02, 2.03, 2.04, 2.05, 2.06, 2.07, 2.08, 2.09, 2.10, 2.11, 2.12, 2.13, 2.14, 2.15, 2.16, 2.17, 2.18, 2.19, 2.20, 2.21, 2.22, 2.23, 2.24, 2.25, 2.26, 2.27, 2.28, 2.29, 2.30, 2.31, 2.32, 2.33, 2.34, 2.35, 2.36, 2.37, 2.38, 2.39, 2.40, 2.41, 2.42, 2.43, 2.44, 2.45, 2.46, 2.47, 2.48, 2.49, 2.50, 2.51, 2.52, 2.53, 2.54, 2.55, 2.56, 2.57, 2.58, 2.59, 2.60, 2.61, 2.62, 2.63, 2.64, 2.65, 2.66, 2.67, 2.68, 2.69, 2.70, 2.71, 2.72, 2.73, 2.74, 2.75, 2.76, 2.77, 2.78, 2.79, 2.80, 2.81, 2.82, 2.83, 2.84, 2.85, 2.86, 2.87, 2.88, 2.89, 2.90, 2.91, 2.92, 2.93, 2.94, 2.95, 2.96, 2.97, 2.98, 2.99, 3.00, 3.01, 3.02, 3.03, 3.04, 3.05, 3.06, 3.07, 3.08, 3.09, 3.10, 3.11, 3.12, 3.13, 3.14, 3.15, 3.16, 3.17, 3.18, 3.19, 3.20, 3.21, 3.22, 3.23, 3.24, 3.25, 3.26, 3.27, 3.28, 3.29, 3.30, 3.31, 3.32, 3.33, 3.34, 3.35, 3.36, 3.37, 3.38, 3.39, 3.40, 3.41, 3.42, 3.43, 3.44, 3.45, 3.46, 3.47, 3.48, 3.49, 3.50, 3.51, 3.52, 3.53, 3.54, 3.55, 3.56, 3.57, 3.58, 3.59, 3.60, 3.61, 3.62, 3.63, 3.64, 3.65, 3.66, 3.67, 3.68, 3.69, 3.70, 3.71, 3.72, 3.73, 3.74, 3.75, 3.76, 3.77, 3.78, 3.79, 3.80, 3.81, 3.82, 3.83, 3.84, 3.85, 3.86, 3.87, 3.88, 3.89, 3.90, 3.91, 3.92, 3.93, 3.94, 3.95, 3.96, 3.97, 3.98, 3.99, 4.00}, "m_{#mu#mu} (GeV/c^{2})"};
    const o2::framework::AxisSpec axisPt{{0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.5, 3, 3.5, 4, 4.5, 5, 6, 7, 8, 9, 10}, "p_{T,#mu} (GeV/c)"};
    const o2::framework::AxisSpec axisEta{40, -4.0, -2.0, "#eta_{#mu}"};
    const o2::framework::AxisSpec axisPhi{36, 0, 2 * M_PI, "#varphi_{#mu} (rad.)"};
    const o2::framework::AxisSpec axisDEta{100, -0.5, 0.5, "#Delta#eta"};
    const o2::framework::AxisSpec axisDPhi{90, -M_PI / 4, M_PI / 4, "#Delta#varphi (rad.)"};
    const o2::framework::AxisSpec axisDCAxy{100, 0, 1, "DCA_{xy} (cm)"};
    const o2::framework::AxisSpec axisChi2{200, 0, 20, "global #chi^{2}/ndf"};
    const o2::framework::AxisSpec axisChi2MFT{100, 0, 10, "MFTsa #chi^{2}/ndf"};
    const o2::framework::AxisSpec axisChi2MFTMCH{100, 0, 100, "matching #chi^{2}/ndf between MFT-MCH"};
    const o2::framework::AxisSpec axisPDCA{100, 0, 1000, "p #times DCA_{xy} (GeV/c #upoint cm)"};
    const o2::framework::AxisSpec axisRabs{150, 15, 90, "R_{xy} at absorber end (cm)"};
    const o2::framework::AxisSpec axisNclsMCH{21, -0.5, 20.5, "N_{cls}^{MCH}"};
    const o2::framework::AxisSpec axisNclsMFT{11, -0.5, 10.5, "N_{cls}^{MFT}"};

    // for track cut efficiency with tag and probe
    fRegistry.add("TAP/same/uls/hs", "dimuon for T&P", o2::framework::HistType::kTHnSparseD, {axisMllLMR, axisPt, axisEta, axisPhi, axisDEta, axisDPhi, axisDCAxy, axisChi2, axisChi2MFT, axisChi2MFTMCH, axisPDCA, axisRabs, axisNclsMCH, axisNclsMFT}, false);
    fRegistry.add("TAP/same/lspp/hs", "dimuon for T&P", o2::framework::HistType::kTHnSparseD, {axisMllLMR, axisPt, axisEta, axisPhi, axisDEta, axisDPhi, axisDCAxy, axisChi2, axisChi2MFT, axisChi2MFTMCH, axisPDCA, axisRabs, axisNclsMCH, axisNclsMFT}, false);
    fRegistry.add("TAP/same/lsmm/hs", "dimuon for T&P", o2::framework::HistType::kTHnSparseD, {axisMllLMR, axisPt, axisEta, axisPhi, axisDEta, axisDPhi, axisDCAxy, axisChi2, axisChi2MFT, axisChi2MFTMCH, axisPDCA, axisRabs, axisNclsMCH, axisNclsMFT}, false);
  }

  std::vector<std::tuple<int, int, int>> vec_min_chi2MatchMCHMFT;

  template <typename TCollision, typename TFwdTrack, typename TFwdTracks, typename TMFTTracks>
  void findBestMatchPerMCHMID(TCollision const&, TFwdTrack const& fwdtrack, TFwdTracks const& fwdtracks, TMFTTracks const&)
  {
    if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
      return;
    }

    const auto& muons_per_MCHMID = fwdtracks.sliceBy(fwdtracksPerMCHTrack, fwdtrack.globalIndex());
    // LOGF(info, "stanadalone: muon.globalIndex() = %d, muon.chi2MatchMCHMFT() = %f", muon.globalIndex(), muon.chi2MatchMCHMFT());
    // LOGF(info, "muons_per_MCHMID.size() = %d", muons_per_MCHMID.size());

    float min_chi2MatchMCHMFT = FLT_MAX;
    std::tuple<int, int, int> tupleIds_at_min_chi2mftmch;
    for (const auto& muon_tmp : muons_per_MCHMID) {
      if (muon_tmp.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        if (0.f < muon_tmp.chi2MatchMCHMFT() && muon_tmp.chi2MatchMCHMFT() < min_chi2MatchMCHMFT) {
          min_chi2MatchMCHMFT = muon_tmp.chi2MatchMCHMFT();
          tupleIds_at_min_chi2mftmch = std::make_tuple(muon_tmp.globalIndex(), muon_tmp.matchMCHTrackId(), muon_tmp.matchMFTTrackId());
        }
      }
    }
    vec_min_chi2MatchMCHMFT.emplace_back(tupleIds_at_min_chi2mftmch);
  }

  template <bool fillHistograms, typename TCollision, typename TFwdTrack, typename TFwdTracks, typename TMFTTracks>
  int countMatchedCandidatesPerMCHMID(TCollision const& collision, TFwdTrack const& fwdtrack, TFwdTracks const& fwdtracks, TMFTTracks const& mfttracks)
  {
    if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::MuonStandaloneTrack) {
      return -1;
    }

    int nGoodGLPerSA = 0;
    auto muons_per_MCHMID = fwdtracks.sliceBy(fwdtracksPerMCHTrack, fwdtrack.globalIndex());
    for (const auto& muon_tmp : muons_per_MCHMID) {
      if (muon_tmp.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
        float pt = 999.f, eta = 999.f, phi = 999.f;
        if (isSelectedGlobalMuon<fillHistograms>(collision, muon_tmp, fwdtracks, mfttracks, pt, eta, phi)) {
          nGoodGLPerSA++;
        }
      }
    }
    return nGoodGLPerSA;
  }

  template <typename TCollision, typename TFwdTrack, typename TFwdTracks, typename TMFTTracks>
  bool isBestMatch(TCollision const& collision, TFwdTrack const& fwdtrack, TFwdTracks const& fwdtracks, TMFTTracks const& mfttracks)
  {
    if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
      std::map<int, float> map_chi2MCHMFT;
      map_chi2MCHMFT[fwdtrack.globalIndex()] = fwdtrack.chi2MatchMCHMFT(); // add myself
      // LOGF(info, "add myself: fwdtrack.globalIndex() = %d, fwdtrack.chi2MatchMCHMFT() = %f", fwdtrack.globalIndex(), fwdtrack.chi2MatchMCHMFT());

      auto candidates = fwdtracks.sliceBy(perMFTTrack, fwdtrack.matchMFTTrackId()); // global muon candidates including this fwdtrack.
      for (const auto& candidate : candidates) {
        if (candidate.globalIndex() == fwdtrack.globalIndex()) {
          continue; // don't add fwdtrack.globalIndex() again.
        }

        float pt = 999.f, eta = 999.f, phi = 999.f;
        if (isSelectedGlobalMuon<false>(collision, candidate, fwdtracks, mfttracks, pt, eta, phi)) {
          map_chi2MCHMFT[candidate.globalIndex()] = candidate.chi2MatchMCHMFT();
          // LOGF(info, "same MFT found: candidate.globalIndex() = %d, candidate.chi2MatchMCHMFT() = %f", candidate.globalIndex(), candidate.chi2MatchMCHMFT());
        }
      }

      auto it = std::min_element(map_chi2MCHMFT.begin(), map_chi2MCHMFT.end(), [](decltype(map_chi2MCHMFT)::value_type& l, decltype(map_chi2MCHMFT)::value_type& r) -> bool { return l.second < r.second; });

      // LOGF(info, "min: globalIndex = %d, chi2 = %f", it->first, it->second);
      // LOGF(info, "bool = %d", it->first == fwdtrack.globalIndex());

      if (it->first == fwdtrack.globalIndex()) { // search for minimum matching-chi2
        map_chi2MCHMFT.clear();
        return true;
      } else {
        map_chi2MCHMFT.clear();
        return false;
      }
    } else {
      return true;
    }
  }

  template <bool fillHistograms = false, typename TCollision, typename TFwdTrack, typename TFwdTracks, typename TMFTTracks>
  bool isSelectedGlobalMuon(TCollision const& collision, TFwdTrack const& fwdtrack, TFwdTracks const&, TMFTTracks const&, float& pt, float& eta, float& phi)
  {
    if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
      return false;
    }

    if (fwdtrack.chi2MatchMCHMFT() < 0.f || glMuonCutGroup.maxMatchingChi2MCHMFT < fwdtrack.chi2MatchMCHMFT()) {
      return false;
    }

    if (fwdtrack.chi2MatchMCHMID() < 0.f) { // this should never happen. only for protection.
      return false;
    }

    if (fwdtrack.rAtAbsorberEnd() < glMuonCutGroup.minRabs || glMuonCutGroup.maxRabs < fwdtrack.rAtAbsorberEnd()) {
      return false;
    }

    auto mchtrack = fwdtrack.template matchMCHTrack_as<TFwdTracks>(); // MCH-MID
    auto mfttrack = fwdtrack.template matchMFTTrack_as<TMFTTracks>(); // MFTsa
    float xMFT = mfttrack.x();
    float yMFT = mfttrack.y();
    float relPtResolution = 0.f;

    if (mfttrack.nClusters() < glMuonCutGroup.minNclsMFT) {
      return false;
    }
    if (mchtrack.nClusters() < glMuonCutGroup.minNclsMCH) {
      return false;
    }

    float rAtAbsorberEnd = fwdtrack.rAtAbsorberEnd(); // this works only for GlobalMuonTrack
    int nClustersMFT = mfttrack.nClusters();
    int ndf_mchmft = 2.f * (mchtrack.nClusters() + nClustersMFT) - 5.f;
    float chi2 = fwdtrack.chi2() / static_cast<float>(ndf_mchmft);
    if (chi2 < 0.f || glMuonCutGroup.maxChi2 < chi2) {
      return false;
    }

    int ndf_mft = 2.f * nClustersMFT - 5.f;
    float chi2mft = mfttrack.chi2() / static_cast<float>(ndf_mft);
    if (chi2mft < 0.f || glMuonCutGroup.maxChi2MFT < chi2mft) {
      return false;
    }

    o2::dataformats::GlobalFwdTrack propmuonAtPV_Matched = o2::aod::fwdtrackutils::propagateMuon(mchtrack, mchtrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToVertex, glMuonCutGroup.matchingZ, mBz, 0.f, 0.f, 0.f);
    float etaMatchedMCHMID = propmuonAtPV_Matched.getEta();
    float phiMatchedMCHMID = propmuonAtPV_Matched.getPhi();
    phiMatchedMCHMID = RecoDecay::constrainAngle(phiMatchedMCHMID, 0, 1U);

    o2::dataformats::GlobalFwdTrack propmuonAtPV = mfttrack.y() > 0.f ? o2::aod::fwdtrackutils::propagateMuon(fwdtrack, fwdtrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToVertex, glMuonCutGroup.matchingZ, mBz, mXShiftMFTtop, mYShiftMFTtop, mZShiftMFTtop) : o2::aod::fwdtrackutils::propagateMuon(fwdtrack, fwdtrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToVertex, glMuonCutGroup.matchingZ, mBz, mXShiftMFTbottom, mYShiftMFTbottom, mZShiftMFTbottom);
    pt = propmuonAtPV.getPt();
    eta = propmuonAtPV.getEta();
    phi = propmuonAtPV.getPhi();
    phi = RecoDecay::constrainAngle(phi, 0, 1U);
    if (glMuonCutGroup.refitGlobalMuon) {
      pt = propmuonAtPV_Matched.getP() * std::sin(2.f * std::atan(std::exp(-eta)));
    }
    relPtResolution = std::sqrt(propmuonAtPV.getSigma2InvQPt()) / std::fabs(propmuonAtPV.getInvQPt());

    float deta = etaMatchedMCHMID - eta;
    float dphi = phiMatchedMCHMID - phi;
    dphi = RecoDecay::constrainAngle(dphi, -M_PI, 1U);
    if (std::sqrt(std::pow(deta / glMuonCutGroup.maxDEta, 2) + std::pow(dphi / glMuonCutGroup.maxDPhi, 2)) > 1.f) {
      return false;
    }

    o2::dataformats::GlobalFwdTrack propmuonAtDCA = mfttrack.y() > 0.f ? o2::aod::fwdtrackutils::propagateMuon(fwdtrack, fwdtrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToDCA, glMuonCutGroup.matchingZ, mBz, mXShiftMFTtop, mYShiftMFTtop, mZShiftMFTtop) : o2::aod::fwdtrackutils::propagateMuon(fwdtrack, fwdtrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToDCA, glMuonCutGroup.matchingZ, mBz, mXShiftMFTbottom, mYShiftMFTbottom, mZShiftMFTbottom);
    float dcaX = propmuonAtDCA.getX() - collision.posX();
    float dcaY = propmuonAtDCA.getY() - collision.posY();
    float dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
    float cXXatDCA = propmuonAtDCA.getSigma2X();
    float cYYatDCA = propmuonAtDCA.getSigma2Y();
    float cXYatDCA = propmuonAtDCA.getSigmaXY();
    if (glMuonCutGroup.maxDCAxy < dcaXY) {
      return false;
    }

    // add convariance matrix of vertex
    float tx = std::cos(propmuonAtDCA.getPhi()) / propmuonAtDCA.getTgl(); // dx/dz
    float ty = std::sin(propmuonAtDCA.getPhi()) / propmuonAtDCA.getTgl(); // dy/dz

    // PV covariance projected into the (dx, dy) residual space.
    float resPVXX = collision.covXX() - 2.0 * tx * collision.covXZ() + tx * tx * collision.covZZ();
    float resPVXY = collision.covXY() - ty * collision.covXZ() - tx * collision.covYZ() + tx * ty * collision.covZZ();
    float resPVYY = collision.covYY() - 2.0 * ty * collision.covYZ() + ty * ty * collision.covZZ();

    // Full residual covariance.
    float cXX = cXXatDCA + resPVXX;
    float cXY = cXYatDCA + resPVXY;
    float cYY = cYYatDCA + resPVYY;

    float det = cXX * cYY - cXY * cXY; // determinanat
    float dcaXYinSigma = 999.f;
    if (det < 0) {
      dcaXYinSigma = 999.f;
    } else {
      dcaXYinSigma = std::sqrt(std::fabs((dcaX * dcaX * cYY + dcaY * dcaY * cXX - 2.f * dcaX * dcaY * cXY) / det / 2.f)); // dca xy in sigma
    }
    float sigma_dcaXY = dcaXY / dcaXYinSigma / std::sqrt(2.f);

    o2::dataformats::GlobalFwdTrack propmuonAtDCA_Matched = o2::aod::fwdtrackutils::propagateMuon(mchtrack, mchtrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToDCA, glMuonCutGroup.matchingZ, mBz, 0.f, 0.f, 0.f);
    float dcaX_Matched = propmuonAtDCA_Matched.getX() - collision.posX();
    float dcaY_Matched = propmuonAtDCA_Matched.getY() - collision.posY();
    float dcaXY_Matched = std::sqrt(dcaX_Matched * dcaX_Matched + dcaY_Matched * dcaY_Matched);
    float pDCA = mchtrack.p() * dcaXY_Matched;

    if (pt < glMuonCutGroup.minPt || glMuonCutGroup.maxPt < pt) {
      return false;
    }

    if (eta < glMuonCutGroup.minEta || glMuonCutGroup.maxEta < eta) {
      return false;
    }

    float chi2IP = mfttrack.y() > 0.f ? o2::aod::fwdtrackutils::getFwdChi2IP(fwdtrack, collision, mBz, mXShiftMFTtop, mYShiftMFTtop, mZShiftMFTtop) : o2::aod::fwdtrackutils::getFwdChi2IP(fwdtrack, collision, mBz, mXShiftMFTbottom, mYShiftMFTbottom, mZShiftMFTbottom);

    if constexpr (fillHistograms) {
      if (fwdtrack.sign() > 0) {
        fRegistry.fill(HIST("MFTMCHMID/positive/hPt"), pt);
        fRegistry.fill(HIST("MFTMCHMID/positive/hPtResolution"), pt, relPtResolution);
        fRegistry.fill(HIST("MFTMCHMID/positive/hEtaPhi"), phi, eta);
        fRegistry.fill(HIST("MFTMCHMID/positive/hEtaPhi_MatchedMCHMID"), phiMatchedMCHMID, etaMatchedMCHMID);
        fRegistry.fill(HIST("MFTMCHMID/positive/hDEtaDPhi"), dphi, deta);
        fRegistry.fill(HIST("MFTMCHMID/positive/hSign"), fwdtrack.sign());
        fRegistry.fill(HIST("MFTMCHMID/positive/hNclusters"), fwdtrack.nClusters());
        fRegistry.fill(HIST("MFTMCHMID/positive/hNclustersMFT"), nClustersMFT);
        fRegistry.fill(HIST("MFTMCHMID/positive/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
        fRegistry.fill(HIST("MFTMCHMID/positive/hRatAbsorberEnd"), rAtAbsorberEnd);
        fRegistry.fill(HIST("MFTMCHMID/positive/hChi2_Pt"), pt, chi2);
        fRegistry.fill(HIST("MFTMCHMID/positive/hChi2MFT_Pt"), pt, chi2mft);
        fRegistry.fill(HIST("MFTMCHMID/positive/hChi2MatchMCHMID_Pt"), pt, fwdtrack.chi2MatchMCHMID());
        fRegistry.fill(HIST("MFTMCHMID/positive/hChi2MatchMCHMFT_Pt"), pt, fwdtrack.chi2MatchMCHMFT());
        fRegistry.fill(HIST("MFTMCHMID/positive/hDCAxy2D"), dcaX, dcaY);
        fRegistry.fill(HIST("MFTMCHMID/positive/hDCAxy2DinSigma"), dcaX / std::sqrt(cXX), dcaY / std::sqrt(cYY));
        fRegistry.fill(HIST("MFTMCHMID/positive/hDCAxy"), dcaXY);
        fRegistry.fill(HIST("MFTMCHMID/positive/hDCAxyinSigma"), dcaXYinSigma);
        fRegistry.fill(HIST("MFTMCHMID/positive/hLog10Chi2IP"), std::log10(chi2IP));
        fRegistry.fill(HIST("MFTMCHMID/positive/hDCAxResolutionvsPt"), pt, std::sqrt(cXX) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/positive/hDCAyResolutionvsPt"), pt, std::sqrt(cYY) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/positive/hDCAxyResolutionvsPt"), pt, sigma_dcaXY * 1e+4);   // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/positive/hDCAx_Phi"), std::atan2(yMFT, xMFT), dcaX);
        fRegistry.fill(HIST("MFTMCHMID/positive/hDCAy_Phi"), std::atan2(yMFT, xMFT), dcaY);
        if (std::atan2(yMFT, xMFT) > 0.f) {
          fRegistry.fill(HIST("MFTMCHMID/positive/hDCAx_PosZ_MFTtop"), collision.posZ(), dcaX);
          fRegistry.fill(HIST("MFTMCHMID/positive/hDCAy_PosZ_MFTtop"), collision.posZ(), dcaY);
        } else {
          fRegistry.fill(HIST("MFTMCHMID/positive/hDCAx_PosZ_MFTbottom"), collision.posZ(), dcaX);
          fRegistry.fill(HIST("MFTMCHMID/positive/hDCAy_PosZ_MFTbottom"), collision.posZ(), dcaY);
        }
      } else {
        fRegistry.fill(HIST("MFTMCHMID/negative/hPt"), pt);
        fRegistry.fill(HIST("MFTMCHMID/negative/hPtResolution"), pt, relPtResolution);
        fRegistry.fill(HIST("MFTMCHMID/negative/hEtaPhi"), phi, eta);
        fRegistry.fill(HIST("MFTMCHMID/negative/hEtaPhi_MatchedMCHMID"), phiMatchedMCHMID, etaMatchedMCHMID);
        fRegistry.fill(HIST("MFTMCHMID/negative/hDEtaDPhi"), dphi, deta);
        fRegistry.fill(HIST("MFTMCHMID/negative/hSign"), fwdtrack.sign());
        fRegistry.fill(HIST("MFTMCHMID/negative/hNclusters"), fwdtrack.nClusters());
        fRegistry.fill(HIST("MFTMCHMID/negative/hNclustersMFT"), nClustersMFT);
        fRegistry.fill(HIST("MFTMCHMID/negative/hPDCA_Rabs"), rAtAbsorberEnd, pDCA);
        fRegistry.fill(HIST("MFTMCHMID/negative/hRatAbsorberEnd"), rAtAbsorberEnd);
        fRegistry.fill(HIST("MFTMCHMID/negative/hChi2_Pt"), pt, chi2);
        fRegistry.fill(HIST("MFTMCHMID/negative/hChi2MFT_Pt"), pt, chi2mft);
        fRegistry.fill(HIST("MFTMCHMID/negative/hChi2MatchMCHMID_Pt"), pt, fwdtrack.chi2MatchMCHMID());
        fRegistry.fill(HIST("MFTMCHMID/negative/hChi2MatchMCHMFT_Pt"), pt, fwdtrack.chi2MatchMCHMFT());
        fRegistry.fill(HIST("MFTMCHMID/negative/hDCAxy2D"), dcaX, dcaY);
        fRegistry.fill(HIST("MFTMCHMID/negative/hDCAxy2DinSigma"), dcaX / std::sqrt(cXX), dcaY / std::sqrt(cYY));
        fRegistry.fill(HIST("MFTMCHMID/negative/hDCAxy"), dcaXY);
        fRegistry.fill(HIST("MFTMCHMID/negative/hDCAxyinSigma"), dcaXYinSigma);
        fRegistry.fill(HIST("MFTMCHMID/negative/hLog10Chi2IP"), std::log10(chi2IP));
        fRegistry.fill(HIST("MFTMCHMID/negative/hDCAxResolutionvsPt"), pt, std::sqrt(cXX) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/negative/hDCAyResolutionvsPt"), pt, std::sqrt(cYY) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/negative/hDCAxyResolutionvsPt"), pt, sigma_dcaXY * 1e+4);   // convert cm to um
        fRegistry.fill(HIST("MFTMCHMID/negative/hDCAx_Phi"), std::atan2(yMFT, xMFT), dcaX);
        fRegistry.fill(HIST("MFTMCHMID/negative/hDCAy_Phi"), std::atan2(yMFT, xMFT), dcaY);
        if (std::atan2(yMFT, xMFT) > 0.f) {
          fRegistry.fill(HIST("MFTMCHMID/negative/hDCAx_PosZ_MFTtop"), collision.posZ(), dcaX);
          fRegistry.fill(HIST("MFTMCHMID/negative/hDCAy_PosZ_MFTtop"), collision.posZ(), dcaY);
        } else {
          fRegistry.fill(HIST("MFTMCHMID/negative/hDCAx_PosZ_MFTbottom"), collision.posZ(), dcaX);
          fRegistry.fill(HIST("MFTMCHMID/negative/hDCAy_PosZ_MFTbottom"), collision.posZ(), dcaY);
        }
      }
    }
    return true;
  }

  template <typename TCollision, typename TFwdTrack, typename TFwdTracks, typename TMFTTracks, typename TMFTCovs>
  bool isTaggedMuon(TCollision const& collision, TFwdTrack const& fwdtrack, TFwdTracks const&, TMFTTracks const&, TMFTCovs const&, float& pt, float& eta, float& phi)
  {
    if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
      return false;
    }

    if (std::find(vec_min_chi2MatchMCHMFT.begin(), vec_min_chi2MatchMCHMFT.end(), std::make_tuple(fwdtrack.globalIndex(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId())) == vec_min_chi2MatchMCHMFT.end()) {
      return false;
    }

    if (fwdtrack.chi2MatchMCHMFT() > tagMuonCutGroup.maxMatchingChi2MCHMFT) {
      return false;
    }

    if (fwdtrack.chi2MatchMCHMID() < 0.f) { // this should never happen. only for protection.
      return false;
    }

    if (fwdtrack.rAtAbsorberEnd() < tagMuonCutGroup.minRabs || tagMuonCutGroup.maxRabs < fwdtrack.rAtAbsorberEnd()) {
      return false;
    }

    auto mchtrack = fwdtrack.template matchMCHTrack_as<TFwdTracks>(); // MCH-MID
    auto mfttrack = fwdtrack.template matchMFTTrack_as<TMFTTracks>(); // MFTsa
    if (mfttrack.nClusters() < tagMuonCutGroup.minNclsMFT) {
      return false;
    }
    if (mchtrack.nClusters() < tagMuonCutGroup.minNclsMCH) {
      return false;
    }

    int nClustersMFT = mfttrack.nClusters();
    int ndf_mchmft = 2.f * (mchtrack.nClusters() + nClustersMFT) - 5.f;
    float chi2 = fwdtrack.chi2() / static_cast<float>(ndf_mchmft);
    if (tagMuonCutGroup.maxChi2 < chi2) {
      return false;
    }

    int ndf_mft = 2.f * nClustersMFT - 5.f;
    float chi2mft = mfttrack.chi2() / static_cast<float>(ndf_mft);
    if (tagMuonCutGroup.maxChi2MFT < chi2mft) {
      return false;
    }

    o2::dataformats::GlobalFwdTrack propmuonAtPV_Matched = o2::aod::fwdtrackutils::propagateMuon(mchtrack, mchtrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToVertex, tagMuonCutGroup.matchingZ, mBz, 0.f, 0.f, 0.f);
    float etaMatchedMCHMID = propmuonAtPV_Matched.getEta();
    float phiMatchedMCHMID = propmuonAtPV_Matched.getPhi();
    phiMatchedMCHMID = RecoDecay::constrainAngle(phiMatchedMCHMID, 0, 1U);

    o2::dataformats::GlobalFwdTrack propmuonAtPV = mfttrack.y() > 0.f ? o2::aod::fwdtrackutils::propagateMuon(fwdtrack, fwdtrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToVertex, tagMuonCutGroup.matchingZ, mBz, mXShiftMFTtop, mYShiftMFTtop, mZShiftMFTtop) : o2::aod::fwdtrackutils::propagateMuon(fwdtrack, fwdtrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToVertex, tagMuonCutGroup.matchingZ, mBz, mXShiftMFTbottom, mYShiftMFTbottom, mZShiftMFTbottom);
    pt = propmuonAtPV.getPt();
    eta = propmuonAtPV.getEta();
    phi = propmuonAtPV.getPhi();
    phi = RecoDecay::constrainAngle(phi, 0, 1U);
    if (tagMuonCutGroup.refitGlobalMuon) {
      pt = propmuonAtPV_Matched.getP() * std::sin(2.f * std::atan(std::exp(-eta)));
    }

    float deta = etaMatchedMCHMID - eta;
    float dphi = phiMatchedMCHMID - phi;
    dphi = RecoDecay::constrainAngle(dphi, -M_PI, 1U);

    if (std::sqrt(std::pow(deta / tagMuonCutGroup.maxDEta, 2) + std::pow(dphi / tagMuonCutGroup.maxDPhi, 2)) > 1.f) {
      return false;
    }

    float dcaX = propmuonAtPV.getX() - collision.posX();
    float dcaY = propmuonAtPV.getY() - collision.posY();
    float dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
    if (tagMuonCutGroup.maxDCAxy < dcaXY) {
      return false;
    }

    if (pt < tagMuonCutGroup.minPt || tagMuonCutGroup.maxPt < pt) {
      return false;
    }

    if (eta < tagMuonCutGroup.minEta || tagMuonCutGroup.maxEta < eta) {
      return false;
    }

    return true;
  }

  template <typename TCollision, typename TFwdTrack, typename TFwdTracks, typename TMFTTracks, typename TMFTCovs>
  bool isProbeMuon(TCollision const& collision, TFwdTrack const& fwdtrack, TFwdTracks const&, TMFTTracks const&, TMFTCovs const&, float& pt, float& eta, float& phi, float& deta, float& dphi, float& dcaXY, float& chi2, float& chi2mft, float& chi2MFTMCH, float& pDCA, float& rAtAbsorberEnd, int& nclsMCH, int& nclsMFT)
  {
    if (fwdtrack.trackType() != o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
      return false;
    }

    if (std::find(vec_min_chi2MatchMCHMFT.begin(), vec_min_chi2MatchMCHMFT.end(), std::make_tuple(fwdtrack.globalIndex(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId())) == vec_min_chi2MatchMCHMFT.end()) {
      return false;
    }

    if (fwdtrack.chi2MatchMCHMFT() > probeMuonCutGroup.maxMatchingChi2MCHMFT) {
      return false;
    }

    if (fwdtrack.chi2MatchMCHMID() < 0.f) { // this should never happen. only for protection.
      return false;
    }

    if (fwdtrack.rAtAbsorberEnd() < probeMuonCutGroup.minRabs || probeMuonCutGroup.maxRabs < fwdtrack.rAtAbsorberEnd()) {
      return false;
    }
    rAtAbsorberEnd = fwdtrack.rAtAbsorberEnd();
    chi2MFTMCH = fwdtrack.chi2MatchMCHMFT();

    auto mchtrack = fwdtrack.template matchMCHTrack_as<TFwdTracks>(); // MCH-MID
    auto mfttrack = fwdtrack.template matchMFTTrack_as<TMFTTracks>(); // MFTsa
    if (mfttrack.nClusters() < probeMuonCutGroup.minNclsMFT) {
      return false;
    }
    if (mchtrack.nClusters() < probeMuonCutGroup.minNclsMCH) {
      return false;
    }

    nclsMCH = mchtrack.nClusters();
    nclsMFT = mfttrack.nClusters();
    int ndf_mchmft = 2.f * (nclsMCH + nclsMFT) - 5.f;
    chi2 = fwdtrack.chi2() / static_cast<float>(ndf_mchmft);
    if (probeMuonCutGroup.maxChi2 < chi2) {
      return false;
    }

    int ndf_mft = 2.f * nclsMFT - 5.f;
    chi2mft = mfttrack.chi2() / static_cast<float>(ndf_mft);
    if (probeMuonCutGroup.maxChi2MFT < chi2mft) {
      return false;
    }

    o2::dataformats::GlobalFwdTrack propmuonAtPV_Matched = o2::aod::fwdtrackutils::propagateMuon(mchtrack, mchtrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToVertex, probeMuonCutGroup.matchingZ, mBz, 0.f, 0.f, 0.f);
    float etaMatchedMCHMID = propmuonAtPV_Matched.getEta();
    float phiMatchedMCHMID = propmuonAtPV_Matched.getPhi();
    phiMatchedMCHMID = RecoDecay::constrainAngle(phiMatchedMCHMID, 0, 1U);

    o2::dataformats::GlobalFwdTrack propmuonAtPV = mfttrack.y() > 0.f ? o2::aod::fwdtrackutils::propagateMuon(fwdtrack, fwdtrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToVertex, probeMuonCutGroup.matchingZ, mBz, mXShiftMFTtop, mYShiftMFTtop, mZShiftMFTtop) : o2::aod::fwdtrackutils::propagateMuon(fwdtrack, fwdtrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToVertex, probeMuonCutGroup.matchingZ, mBz, mXShiftMFTbottom, mYShiftMFTbottom, mZShiftMFTbottom);
    pt = propmuonAtPV.getPt();
    eta = propmuonAtPV.getEta();
    phi = propmuonAtPV.getPhi();
    phi = RecoDecay::constrainAngle(phi, 0, 1U);
    if (probeMuonCutGroup.refitGlobalMuon) {
      pt = propmuonAtPV_Matched.getP() * std::sin(2.f * std::atan(std::exp(-eta)));
    }

    deta = etaMatchedMCHMID - eta;
    dphi = phiMatchedMCHMID - phi;
    dphi = RecoDecay::constrainAngle(dphi, -M_PI, 1U);

    if (std::sqrt(std::pow(deta / probeMuonCutGroup.maxDEta, 2) + std::pow(dphi / probeMuonCutGroup.maxDPhi, 2)) > 1.f) {
      return false;
    }

    float dcaX = propmuonAtPV.getX() - collision.posX();
    float dcaY = propmuonAtPV.getY() - collision.posY();
    dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
    if (probeMuonCutGroup.maxDCAxy < dcaXY) {
      return false;
    }

    o2::dataformats::GlobalFwdTrack propmuonAtDCA_Matched = o2::aod::fwdtrackutils::propagateMuon(mchtrack, mchtrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToDCA, probeMuonCutGroup.matchingZ, mBz, 0.f, 0.f, 0.f);
    float dcaX_Matched = propmuonAtDCA_Matched.getX() - collision.posX();
    float dcaY_Matched = propmuonAtDCA_Matched.getY() - collision.posY();
    float dcaXY_Matched = std::sqrt(dcaX_Matched * dcaX_Matched + dcaY_Matched * dcaY_Matched);
    pDCA = mchtrack.p() * dcaXY_Matched;

    if (pt < probeMuonCutGroup.minPt || probeMuonCutGroup.maxPt < pt) {
      return false;
    }

    if (eta < probeMuonCutGroup.minEta || probeMuonCutGroup.maxEta < eta) {
      return false;
    }

    return true;
  }

  template <typename TCollision, typename TMuons, typename TFwdTracks, typename TMFTTracks, typename TMFTCovs>
  void runPairing(TCollision const& collision, TMuons const& posMuons, TMuons const& negMuons, TFwdTracks const& fwdtracks, TMFTTracks const& mfttracks, TMFTCovs const&)
  {
    for (const auto& pos : posMuons) {
      auto v1 = std::get<1>(pos);
      auto fwdtrack1 = fwdtracks.rawIteratorAt(std::get<0>(pos));
      if (!isBestMatch(collision, fwdtrack1, fwdtracks, mfttracks)) {
        continue;
      }

      for (const auto& neg : negMuons) {
        auto v2 = std::get<1>(neg);
        auto fwdtrack2 = fwdtracks.rawIteratorAt(std::get<0>(neg));
        if (!isBestMatch(collision, fwdtrack2, fwdtracks, mfttracks)) {
          continue;
        }

        auto v12 = v1 + v2;
        fRegistry.fill(HIST("Pair/same/uls/hs"), v12.M(), v12.Pt(), v12.Rapidity());

      } // end end of negative muon loop
    } // end end of positive muon loop

    for (int i1 = 0; i1 < static_cast<int>(posMuons.size()) - 1; i1++) {
      auto pos1 = posMuons[i1];
      auto v1 = std::get<1>(pos1);
      auto fwdtrack1 = fwdtracks.rawIteratorAt(std::get<0>(pos1));
      if (!isBestMatch(collision, fwdtrack1, fwdtracks, mfttracks)) {
        continue;
      }

      for (int i2 = i1 + 1; i2 < static_cast<int>(posMuons.size()); i2++) {
        auto pos2 = posMuons[i2];
        auto v2 = std::get<1>(pos2);
        auto fwdtrack2 = fwdtracks.rawIteratorAt(std::get<0>(pos2));
        // if (std::get<0>(pos1) == std::get<0>(pos2)) {
        //   continue;
        // }

        if (!isBestMatch(collision, fwdtrack2, fwdtracks, mfttracks)) {
          continue;
        }

        auto v12 = v1 + v2;
        fRegistry.fill(HIST("Pair/same/lspp/hs"), v12.M(), v12.Pt(), v12.Rapidity());
      } // end end of positive muon loop
    } // end end of positive muon loop

    for (int i1 = 0; i1 < static_cast<int>(negMuons.size()) - 1; i1++) {
      auto neg1 = negMuons[i1];
      auto v1 = std::get<1>(neg1);
      auto fwdtrack1 = fwdtracks.rawIteratorAt(std::get<0>(neg1));
      if (!isBestMatch(collision, fwdtrack1, fwdtracks, mfttracks)) {
        continue;
      }

      for (int i2 = i1 + 1; i2 < static_cast<int>(negMuons.size()); i2++) {
        auto neg2 = negMuons[i2];
        auto v2 = std::get<1>(neg2);
        auto fwdtrack2 = fwdtracks.rawIteratorAt(std::get<0>(neg2));
        // if (std::get<0>(neg1) == std::get<0>(neg2)) {
        //   continue;
        // }

        if (!isBestMatch(collision, fwdtrack2, fwdtracks, mfttracks)) {
          continue;
        }

        auto v12 = v1 + v2;
        fRegistry.fill(HIST("Pair/same/lsmm/hs"), v12.M(), v12.Pt(), v12.Rapidity());
      } // end end of negative muon loop
    } // end end of negative muon loop
  }

  template <typename TCollision, typename TMuons, typename TFwdTracks, typename TMFTTracks, typename TMFTCovs>
  void runTAP(TCollision const& collision, TMuons const& posMuons, TMuons const& negMuons, TFwdTracks const& fwdtracks, TMFTTracks const& mfttracks, TMFTCovs const& mftCovs)
  {
    // pos is tag, neg is probe
    for (const auto& pos : posMuons) {
      auto fwdtrack1 = fwdtracks.rawIteratorAt(pos);
      if (!isBestMatch(collision, fwdtrack1, fwdtracks, mfttracks)) {
        continue;
      }
      float pt1 = 999.f, eta1 = 999.f, phi1 = 999.f;
      if (!isTaggedMuon(collision, fwdtrack1, fwdtracks, mfttracks, mftCovs, pt1, eta1, phi1)) {
        continue;
      }
      ROOT::Math::PtEtaPhiMVector v1(pt1, eta1, phi1, o2::constants::physics::MassMuon);

      for (const auto& neg : negMuons) {
        auto fwdtrack2 = fwdtracks.rawIteratorAt(neg);
        if (!isBestMatch(collision, fwdtrack2, fwdtracks, mfttracks)) {
          continue;
        }

        float pt2 = 999.f, eta2 = 999.f, phi2 = 999.f;
        float deta = 999.f, dphi = 999.f, dcaXY = 999.f, chi2 = 999.f, chi2MFT = 999.f, chi2MFTMCH = 999.f, pDCA = 999.f, rAtAbsorberEnd = 999.f;
        int nclsMCH = 0, nclsMFT = 0;
        if (!isProbeMuon(collision, fwdtrack2, fwdtracks, mfttracks, mftCovs, pt2, eta2, phi2, deta, dphi, dcaXY, chi2, chi2MFT, chi2MFTMCH, pDCA, rAtAbsorberEnd, nclsMCH, nclsMFT)) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v2(pt2, eta2, phi2, o2::constants::physics::MassMuon);

        auto v12 = v1 + v2;
        fRegistry.fill(HIST("TAP/same/uls/hs"), v12.M(), pt2, eta2, phi2, deta, dphi, dcaXY, chi2, chi2MFT, chi2MFTMCH, pDCA, rAtAbsorberEnd, nclsMCH, nclsMFT);

      } // end end of negative muon loop
    } // end end of positive muon loop

    // neg is tag, pos is probe
    for (const auto& neg : negMuons) {
      auto fwdtrack1 = fwdtracks.rawIteratorAt(neg);
      if (!isBestMatch(collision, fwdtrack1, fwdtracks, mfttracks)) {
        continue;
      }
      float pt1 = 999, eta1 = 999, phi1 = 999;
      if (!isTaggedMuon(collision, fwdtrack1, fwdtracks, mfttracks, mftCovs, pt1, eta1, phi1)) {
        continue;
      }
      ROOT::Math::PtEtaPhiMVector v1(pt1, eta1, phi1, o2::constants::physics::MassMuon);

      for (const auto& pos : posMuons) {
        auto fwdtrack2 = fwdtracks.rawIteratorAt(pos);
        if (!isBestMatch(collision, fwdtrack2, fwdtracks, mfttracks)) {
          continue;
        }

        float pt2 = 999, eta2 = 999, phi2 = 999;
        float deta = 999.f, dphi = 999.f, dcaXY = 999.f, chi2 = 999.f, chi2MFT = 999.f, chi2MFTMCH = 999.f, pDCA = 999.f, rAtAbsorberEnd = 999.f;
        int nclsMCH = 0, nclsMFT = 0;
        if (!isProbeMuon(collision, fwdtrack2, fwdtracks, mfttracks, mftCovs, pt2, eta2, phi2, deta, dphi, dcaXY, chi2, chi2MFT, chi2MFTMCH, pDCA, rAtAbsorberEnd, nclsMCH, nclsMFT)) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v2(pt2, eta2, phi2, o2::constants::physics::MassMuon);

        auto v12 = v1 + v2;
        fRegistry.fill(HIST("TAP/same/uls/hs"), v12.M(), pt2, eta2, phi2, deta, dphi, dcaXY, chi2, chi2MFT, chi2MFTMCH, pDCA, rAtAbsorberEnd, nclsMCH, nclsMFT);

      } // end end of negative muon loop
    } // end end of positive muon loop

    for (const auto& pos1 : posMuons) {
      auto fwdtrack1 = fwdtracks.rawIteratorAt(pos1);
      if (!isBestMatch(collision, fwdtrack1, fwdtracks, mfttracks)) {
        continue;
      }
      float pt1 = 999, eta1 = 999, phi1 = 999;
      if (!isTaggedMuon(collision, fwdtrack1, fwdtracks, mfttracks, mftCovs, pt1, eta1, phi1)) {
        continue;
      }
      ROOT::Math::PtEtaPhiMVector v1(pt1, eta1, phi1, o2::constants::physics::MassMuon);

      for (const auto& pos2 : posMuons) {
        auto fwdtrack2 = fwdtracks.rawIteratorAt(pos2);
        if (pos1 == pos2) {
          continue;
        }

        if (!isBestMatch(collision, fwdtrack2, fwdtracks, mfttracks)) {
          continue;
        }
        float pt2 = 999, eta2 = 999, phi2 = 999;
        float deta = 999.f, dphi = 999.f, dcaXY = 999.f, chi2 = 999.f, chi2MFT = 999.f, chi2MFTMCH = 999.f, pDCA = 999.f, rAtAbsorberEnd = 999.f;
        int nclsMCH = 0, nclsMFT = 0;
        if (!isProbeMuon(collision, fwdtrack2, fwdtracks, mfttracks, mftCovs, pt2, eta2, phi2, deta, dphi, dcaXY, chi2, chi2MFT, chi2MFTMCH, pDCA, rAtAbsorberEnd, nclsMCH, nclsMFT)) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v2(pt2, eta2, phi2, o2::constants::physics::MassMuon);

        auto v12 = v1 + v2;
        fRegistry.fill(HIST("TAP/same/lspp/hs"), v12.M(), pt2, eta2, phi2, deta, dphi, dcaXY, chi2, chi2MFT, chi2MFTMCH, pDCA, rAtAbsorberEnd, nclsMCH, nclsMFT);
      } // end end of positive muon loop
    } // end end of positive muon loop

    for (const auto& neg1 : negMuons) {
      auto fwdtrack1 = fwdtracks.rawIteratorAt(neg1);
      if (!isBestMatch(collision, fwdtrack1, fwdtracks, mfttracks)) {
        continue;
      }
      float pt1 = 999, eta1 = 999, phi1 = 999;
      if (!isTaggedMuon(collision, fwdtrack1, fwdtracks, mfttracks, mftCovs, pt1, eta1, phi1)) {
        continue;
      }
      ROOT::Math::PtEtaPhiMVector v1(pt1, eta1, phi1, o2::constants::physics::MassMuon);

      for (const auto& neg2 : negMuons) {
        auto fwdtrack2 = fwdtracks.rawIteratorAt(neg2);
        if (neg1 == neg2) {
          continue;
        }

        if (!isBestMatch(collision, fwdtrack2, fwdtracks, mfttracks)) {
          continue;
        }
        float pt2 = 999, eta2 = 999, phi2 = 999;
        float deta = 999.f, dphi = 999.f, dcaXY = 999.f, chi2 = 999.f, chi2MFT = 999.f, chi2MFTMCH = 999.f, pDCA = 999.f, rAtAbsorberEnd = 999.f;
        int nclsMCH = 0, nclsMFT = 0;
        if (!isProbeMuon(collision, fwdtrack2, fwdtracks, mfttracks, mftCovs, pt2, eta2, phi2, deta, dphi, dcaXY, chi2, chi2MFT, chi2MFTMCH, pDCA, rAtAbsorberEnd, nclsMCH, nclsMFT)) {
          continue;
        }
        ROOT::Math::PtEtaPhiMVector v2(pt2, eta2, phi2, o2::constants::physics::MassMuon);

        auto v12 = v1 + v2;
        fRegistry.fill(HIST("TAP/same/lsmm/hs"), v12.M(), pt2, eta2, phi2, deta, dphi, dcaXY, chi2, chi2MFT, chi2MFTMCH, pDCA, rAtAbsorberEnd, nclsMCH, nclsMFT);
      } // end end of negative muon loop
    } // end end of negative muon loop
  }

  template <typename TCollision>
  bool isSelectedCollision(TCollision const& collision)
  {
    if (!(eventCutGroup.minZvtx < collision.posZ() && collision.posZ() < eventCutGroup.maxZvtx)) {
      return false;
    }
    if (eventCutGroup.cfgRequireFT0AND && !collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
      return false;
    }
    if (eventCutGroup.cfgRequireNoTFB && !collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
      return false;
    }
    if (eventCutGroup.cfgRequireNoITSROFB && !collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
      return false;
    }
    if (eventCutGroup.cfgRequireVertexITSTPC && !collision.selection_bit(o2::aod::evsel::kIsVertexITSTPC)) {
      return false;
    }
    if (eventCutGroup.cfgRequireVertexTOFmatched && !collision.selection_bit(o2::aod::evsel::kIsVertexTOFmatched)) {
      return false;
    }
    if (eventCutGroup.cfgRequireNoSameBunchPileup && !collision.selection_bit(o2::aod::evsel::kNoSameBunchPileup)) {
      return false;
    }
    if (eventCutGroup.cfgRequireGoodZvtxFT0vsPV && !collision.selection_bit(o2::aod::evsel::kIsGoodZvtxFT0vsPV)) {
      return false;
    }
    if (collision.chi2() / collision.numContrib() > eventCutGroup.cfgChi2PerNcontribMax) {
      return false;
    }
    if (collision.numContrib() < eventCutGroup.cfgNcontribMin) {
      return false;
    }
    return true;
  }

  using MyCollisions = o2::soa::Join<o2::aod::Collisions, o2::aod::EvSels>;
  using MyCollision = MyCollisions::iterator;

  using MyBCs = o2::soa::Join<o2::aod::BCsWithTimestamps, o2::aod::BcSels>;
  using MyBC = MyBCs::iterator;

  using MyFwdTracks = o2::soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov>; // muon tracks are repeated. i.e. not exclusive.
  using MyFwdTrack = MyFwdTracks::iterator;

  o2::framework::SliceCache cache;
  o2::framework::Preslice<o2::aod::FwdTracks> perCollision = o2::aod::fwdtrack::collisionId;
  o2::framework::Preslice<o2::aod::FwdTrackAssoc> fwdtrackIndicesPerCollision = o2::aod::track_association::collisionId;
  o2::framework::PresliceUnsorted<o2::aod::FwdTrackAssoc> fwdtrackIndicesPerFwdTrack = o2::aod::track_association::fwdtrackId;
  o2::framework::PresliceUnsorted<o2::aod::FwdTracks> fwdtracksPerMCHTrack = o2::aod::fwdtrack::matchMCHTrackId;
  o2::framework::PresliceUnsorted<o2::aod::FwdTracks> perMFTTrack = o2::aod::fwdtrack::matchMFTTrackId;

  // std::unordered_map<int, int> map_mfttrackcovs;
  std::unordered_map<int, bool> mapSelectedCollisions;

  void processSA(MyCollisions const& collisions, MyBCs const&, MyFwdTracks const& fwdtracks, o2::aod::MFTTracks const& mfttracks, o2::aod::MFTTracksCov const& mftCovs)
  {
    // for (const auto& mfttrackConv : mftCovs) {
    //   map_mfttrackcovs[mfttrackConv.matchMFTTrackId()] = mfttrackConv.globalIndex();
    // }

    vec_min_chi2MatchMCHMFT.reserve(fwdtracks.size());
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<MyBCs>();
      initCCDB(bc);
      const auto& fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      for (const auto& fwdtrack : fwdtracks_per_coll) {
        findBestMatchPerMCHMID(collision, fwdtrack, fwdtracks, mfttracks);
      } // end of fwdtrack loop
    } // end of collision loop

    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<MyBCs>();
      initCCDB(bc);
      fRegistry.fill(HIST("hCollisionCounter"), 1);

      if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        fRegistry.fill(HIST("hCollisionCounter"), 2);
      }
      if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        fRegistry.fill(HIST("hCollisionCounter"), 3);
      }
      if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        fRegistry.fill(HIST("hCollisionCounter"), 4);
      }
      if (eventCutGroup.minZvtx < collision.posZ() && collision.posZ() < eventCutGroup.maxZvtx) {
        fRegistry.fill(HIST("hCollisionCounter"), 5);
      }
      if (collision.sel8()) {
        fRegistry.fill(HIST("hCollisionCounter"), 6);
      }

      if (isSelectedCollision(collision)) {
        fRegistry.fill(HIST("hCollisionCounter"), 7);
        fRegistry.fill(HIST("Vertex/MB/hZvtx"), collision.posZ());
        fRegistry.fill(HIST("Vertex/MB/hChi2vsNContrib"), collision.numContrib(), collision.chi2() / collision.numContrib());
        fRegistry.fill(HIST("Vertex/MB/hSigmaX"), collision.numContrib(), std::sqrt(collision.covXX()) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("Vertex/MB/hSigmaY"), collision.numContrib(), std::sqrt(collision.covYY()) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("Vertex/MB/hSigmaZ"), collision.numContrib(), std::sqrt(collision.covZZ()) * 1e+4); // convert cm to um

        int nGlobalMuon = 0;
        // LOGF(info, "collision.globalIndex() = %d, fwdtrackIdsThisCollision = %d", collision.globalIndex(), fwdtrackIdsThisCollision.size());
        auto fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
        for (const auto& fwdtrack : fwdtracks_per_coll) {
          int nGoodGLPerSA = countMatchedCandidatesPerMCHMID<true>(collision, fwdtrack, fwdtracks, mfttracks);
          if (nGoodGLPerSA > 0) {
            nGlobalMuon++;
          }
        }

        // LOGF(info, "nGlobalMuon = %d", nGlobalMuon);
        fRegistry.fill(HIST("hNmu"), nGlobalMuon);

        if (nGlobalMuon >= minNmuon) {
          fRegistry.fill(HIST("hCollisionCounter"), 8);
          fRegistry.fill(HIST("Vertex/Triggered/hZvtx"), collision.posZ());
          fRegistry.fill(HIST("Vertex/Triggered/hChi2vsNContrib"), collision.numContrib(), collision.chi2() / collision.numContrib());
          fRegistry.fill(HIST("Vertex/Triggered/hSigmaX"), collision.numContrib(), std::sqrt(collision.covXX()) * 1e+4); // convert cm to um
          fRegistry.fill(HIST("Vertex/Triggered/hSigmaY"), collision.numContrib(), std::sqrt(collision.covYY()) * 1e+4); // convert cm to um
          fRegistry.fill(HIST("Vertex/Triggered/hSigmaZ"), collision.numContrib(), std::sqrt(collision.covZZ()) * 1e+4); // convert cm to um

          mapSelectedCollisions[collision.globalIndex()] = true;
          tags(true);
        } else {
          mapSelectedCollisions[collision.globalIndex()] = false;
          tags(false);
        }
      } else {
        mapSelectedCollisions[collision.globalIndex()] = false;
        tags(false);
      }
    } // end of collision loop

    // trigger QA
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<MyBCs>();
      initCCDB(bc);

      if (!mapSelectedCollisions[collision.globalIndex()]) {
        continue;
      }
      fRegistry.fill(HIST("hCollisionCounter"), 9);

      auto fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      std::vector<std::pair<int, ROOT::Math::PtEtaPhiMVector>> posMuons;
      std::vector<std::pair<int, ROOT::Math::PtEtaPhiMVector>> negMuons;
      posMuons.reserve(fwdtracks_per_coll.size());
      negMuons.reserve(fwdtracks_per_coll.size());

      for (const auto& fwdtrack : fwdtracks_per_coll) {
        float pt = 999.f, eta = 999.f, phi = 999.f;
        if (isSelectedGlobalMuon<false>(collision, fwdtrack, fwdtracks, mfttracks, pt, eta, phi)) {
          if (fwdtrack.sign() > 0) {
            posMuons.emplace_back(std::make_pair(fwdtrack.globalIndex(), ROOT::Math::PtEtaPhiMVector(pt, eta, phi, o2::constants::physics::MassMuon)));
          } else {
            negMuons.emplace_back(std::make_pair(fwdtrack.globalIndex(), ROOT::Math::PtEtaPhiMVector(pt, eta, phi, o2::constants::physics::MassMuon)));
          }
        }
      } // end of fwdtrack loop

      // make pairs
      runPairing(collision, posMuons, negMuons, fwdtracks, mfttracks, mftCovs);

      posMuons.clear();
      posMuons.shrink_to_fit();
      negMuons.clear();
      negMuons.shrink_to_fit();
    } // end of collision loop

    // for TAP
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<MyBCs>();
      initCCDB(bc);

      if (!isSelectedCollision(collision)) {
        continue;
      }

      auto fwdtracks_per_coll = fwdtracks.sliceBy(perCollision, collision.globalIndex());
      std::vector<int> posProbeMuons;
      std::vector<int> negProbeMuons;
      posProbeMuons.reserve(fwdtracks_per_coll.size());
      negProbeMuons.reserve(fwdtracks_per_coll.size());

      for (const auto& fwdtrack : fwdtracks_per_coll) {
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          if (fwdtrack.sign() > 0) {
            posProbeMuons.emplace_back(fwdtrack.globalIndex());
          } else {
            negProbeMuons.emplace_back(fwdtrack.globalIndex());
          }
        }
      } // end of fwdtrack loop

      // track cut efficiency
      runTAP(collision, posProbeMuons, negProbeMuons, fwdtracks, mfttracks, mftCovs); // currently, select global muons with min matching chi2 MCH-MFT

      posProbeMuons.clear();
      posProbeMuons.shrink_to_fit();
      negProbeMuons.clear();
      negProbeMuons.shrink_to_fit();
    } // end of collision loop

    // map_mfttrackcovs.clear();
    mapSelectedCollisions.clear();
    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
  }
  PROCESS_SWITCH(globalDimuonFilter, processSA, "process without FTTCA", false);

  void processFTTCA(MyCollisions const& collisions, MyBCs const&, MyFwdTracks const& fwdtracks, o2::aod::MFTTracks const& mfttracks, o2::aod::MFTTracksCov const& mftCovs, o2::aod::FwdTrackAssoc const& fwdtrackIndices)
  {
    // for (const auto& mfttrackConv : mftCovs) {
    //   map_mfttrackcovs[mfttrackConv.matchMFTTrackId()] = mfttrackConv.globalIndex();
    // }

    vec_min_chi2MatchMCHMFT.reserve(fwdtracks.size());
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<MyBCs>();
      initCCDB(bc);
      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        findBestMatchPerMCHMID(collision, fwdtrack, fwdtracks, mfttracks);
      } // end of fwdtrack loop
    } // end of collision loop

    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<MyBCs>();
      initCCDB(bc);

      fRegistry.fill(HIST("hCollisionCounter"), 1);

      if (collision.selection_bit(o2::aod::evsel::kIsTriggerTVX)) {
        fRegistry.fill(HIST("hCollisionCounter"), 2);
      }
      if (collision.selection_bit(o2::aod::evsel::kNoTimeFrameBorder)) {
        fRegistry.fill(HIST("hCollisionCounter"), 3);
      }
      if (collision.selection_bit(o2::aod::evsel::kNoITSROFrameBorder)) {
        fRegistry.fill(HIST("hCollisionCounter"), 4);
      }
      if (eventCutGroup.minZvtx < collision.posZ() && collision.posZ() < eventCutGroup.maxZvtx) {
        fRegistry.fill(HIST("hCollisionCounter"), 5);
      }
      if (collision.sel8()) {
        fRegistry.fill(HIST("hCollisionCounter"), 6);
      }

      if (isSelectedCollision(collision)) {
        fRegistry.fill(HIST("hCollisionCounter"), 7);
        fRegistry.fill(HIST("Vertex/MB/hZvtx"), collision.posZ());
        fRegistry.fill(HIST("Vertex/MB/hChi2vsNContrib"), collision.numContrib(), collision.chi2() / collision.numContrib());
        fRegistry.fill(HIST("Vertex/MB/hSigmaX"), collision.numContrib(), std::sqrt(collision.covXX()) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("Vertex/MB/hSigmaY"), collision.numContrib(), std::sqrt(collision.covYY()) * 1e+4); // convert cm to um
        fRegistry.fill(HIST("Vertex/MB/hSigmaZ"), collision.numContrib(), std::sqrt(collision.covZZ()) * 1e+4); // convert cm to um

        int nGlobalMuon = 0;
        auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
        // LOGF(info, "collision.globalIndex() = %d, fwdtrackIdsThisCollision = %d", collision.globalIndex(), fwdtrackIdsThisCollision.size());
        for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
          auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
          int nGoodGLPerSA = countMatchedCandidatesPerMCHMID<true>(collision, fwdtrack, fwdtracks, mfttracks);
          if (nGoodGLPerSA > 0) {
            nGlobalMuon++;
          }
        } // end of fwdtrack loop

        // LOGF(info, "nGlobalMuon = %d", nGlobalMuon);
        fRegistry.fill(HIST("hNmu"), nGlobalMuon);

        if (nGlobalMuon >= minNmuon) {
          fRegistry.fill(HIST("hCollisionCounter"), 8);
          fRegistry.fill(HIST("Vertex/Triggered/hZvtx"), collision.posZ());
          fRegistry.fill(HIST("Vertex/Triggered/hChi2vsNContrib"), collision.numContrib(), collision.chi2() / collision.numContrib());
          fRegistry.fill(HIST("Vertex/Triggered/hSigmaX"), collision.numContrib(), std::sqrt(collision.covXX()) * 1e+4); // convert cm to um
          fRegistry.fill(HIST("Vertex/Triggered/hSigmaY"), collision.numContrib(), std::sqrt(collision.covYY()) * 1e+4); // convert cm to um
          fRegistry.fill(HIST("Vertex/Triggered/hSigmaZ"), collision.numContrib(), std::sqrt(collision.covZZ()) * 1e+4); // convert cm to um

          mapSelectedCollisions[collision.globalIndex()] = true;
          tags(true);
        } else {
          mapSelectedCollisions[collision.globalIndex()] = false;
          tags(false);
        }
      } else {
        mapSelectedCollisions[collision.globalIndex()] = false;
        tags(false);
      }
    } // end of collision loop

    // trigger QA
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<MyBCs>();
      initCCDB(bc);

      if (!mapSelectedCollisions[collision.globalIndex()]) {
        continue;
      }
      fRegistry.fill(HIST("hCollisionCounter"), 9);

      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      std::vector<std::pair<int, ROOT::Math::PtEtaPhiMVector>> posMuons;
      std::vector<std::pair<int, ROOT::Math::PtEtaPhiMVector>> negMuons;
      posMuons.reserve(fwdtrackIdsThisCollision.size());
      negMuons.reserve(fwdtrackIdsThisCollision.size());

      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        if (std::find(vec_min_chi2MatchMCHMFT.begin(), vec_min_chi2MatchMCHMFT.end(), std::make_tuple(fwdtrack.globalIndex(), fwdtrack.matchMCHTrackId(), fwdtrack.matchMFTTrackId())) == vec_min_chi2MatchMCHMFT.end()) {
          continue;
        }

        float pt = 999.f, eta = 999.f, phi = 999.f;
        if (isSelectedGlobalMuon<false>(collision, fwdtrack, fwdtracks, mfttracks, pt, eta, phi)) {
          if (fwdtrack.sign() > 0) {
            posMuons.emplace_back(std::make_pair(fwdtrack.globalIndex(), ROOT::Math::PtEtaPhiMVector(pt, eta, phi, o2::constants::physics::MassMuon)));
          } else {
            negMuons.emplace_back(std::make_pair(fwdtrack.globalIndex(), ROOT::Math::PtEtaPhiMVector(pt, eta, phi, o2::constants::physics::MassMuon)));
          }
        }
      } // end of fwdtrack loop

      // make pairs
      runPairing(collision, posMuons, negMuons, fwdtracks, mfttracks, mftCovs);

      posMuons.clear();
      posMuons.shrink_to_fit();
      negMuons.clear();
      negMuons.shrink_to_fit();
    } // end of collision loop

    // for TAP
    for (const auto& collision : collisions) {
      auto bc = collision.template bc_as<MyBCs>();
      initCCDB(bc);

      if (!isSelectedCollision(collision)) {
        continue;
      }

      auto fwdtrackIdsThisCollision = fwdtrackIndices.sliceBy(fwdtrackIndicesPerCollision, collision.globalIndex());
      std::vector<int> posProbeMuons;
      std::vector<int> negProbeMuons;
      posProbeMuons.reserve(fwdtrackIdsThisCollision.size());
      negProbeMuons.reserve(fwdtrackIdsThisCollision.size());

      for (const auto& fwdtrackId : fwdtrackIdsThisCollision) {
        auto fwdtrack = fwdtrackId.template fwdtrack_as<MyFwdTracks>();
        if (fwdtrack.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack) {
          if (fwdtrack.sign() > 0) {
            posProbeMuons.emplace_back(fwdtrack.globalIndex());
          } else {
            negProbeMuons.emplace_back(fwdtrack.globalIndex());
          }
        }
      } // end of fwdtrack loop

      // track cut efficiency
      runTAP(collision, posProbeMuons, negProbeMuons, fwdtracks, mfttracks, mftCovs); // currently, select global muons with min matching chi2 MCH-MFT

      posProbeMuons.clear();
      posProbeMuons.shrink_to_fit();
      negProbeMuons.clear();
      negProbeMuons.shrink_to_fit();
    } // end of collision loop

    // map_mfttrackcovs.clear();
    mapSelectedCollisions.clear();
    vec_min_chi2MatchMCHMFT.clear();
    vec_min_chi2MatchMCHMFT.shrink_to_fit();
  }
  PROCESS_SWITCH(globalDimuonFilter, processFTTCA, "process with FTTCA", true);
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfg)
{
  return o2::framework::WorkflowSpec{
    adaptAnalysisTask<globalDimuonFilter>(cfg, o2::framework::TaskName{"em-global-dimuon-filter"})};
}
