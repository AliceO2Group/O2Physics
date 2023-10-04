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

// Builder task for hypertriton 3-body decay reconstruction
// author: yuanzhe.wang@cern.ch

#include <cmath>
#include <array>
#include <cstdlib>

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "DCAFitter/DCAFitterN.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/Vtx3BodyTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"

#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsTPC/BetheBlochAleph.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using FullTracksExtIU = soa::Join<aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPr, aod::pidTPCFullPi, aod::pidTPCFullDe>;
using FullTracksExtMCIU = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

using LabeledTracks = soa::Join<FullTracksExtIU, aod::McTrackLabels>;

struct hypertriton3bodyBuilder {

  Produces<aod::StoredVtx3BodyDatas> vtx3bodydata;
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // Configurables
  Configurable<bool> d_UseAbsDCA{"d_UseAbsDCA", true, "Use Abs DCAs"};

  HistogramRegistry registry{
    "registry",
    {
      {"hEventCounter", "hEventCounter", {HistType::kTH1F, {{1, 0.0f, 1.0f}}}},
      {"hVtx3BodyCounter", "hVtx3BodyCounter", {HistType::kTH1F, {{7, 0.0f, 7.0f}}}},
    },
  };

  // Selection criteria
  Configurable<double> d_bz_input{"d_bz", -999, "bz field, -999 is automatic"};
  Configurable<int> mincrossedrows{"mincrossedrows", 70, "min crossed rows"};
  Configurable<float> minCosPA3body{"minCosPA3body", 0.8, "minCosPA3body"};
  Configurable<float> dcavtxdau{"dcavtxdau", 1.0, "DCA Vtx Daughters"};

  Configurable<int> useMatCorrType{"useMatCorrType", 0, "0: none, 1: TGeo, 2: LUT"};
  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  int mRunNumber;
  float d_bz;
  float maxSnp;  // max sine phi for propagation
  float maxStep; // max step size (cm) for propagation
  o2::base::MatLayerCylSet* lut = nullptr;
  o2::vertexing::DCAFitterN<3> fitter3body;

  void init(InitContext& context)
  {
    mRunNumber = 0;
    d_bz = 0;
    maxSnp = 0.85f;  // could be changed later
    maxStep = 2.00f; // could be changed later
    fitter3body.setPropagateToPCA(true);
    fitter3body.setMaxR(200.); //->maxRIni3body
    fitter3body.setMinParamChange(1e-3);
    fitter3body.setMinRelChi2Change(0.9);
    fitter3body.setMaxDZIni(1e9);
    fitter3body.setMaxChi2(1e9);
    fitter3body.setUseAbsDCA(d_UseAbsDCA);

    // Material correction in the DCA fitter

    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setFatalWhenNull(false);

    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(1, "Total");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(2, "CollisionID");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(3, "TPCNcls");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(4, "HasSV");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(5, "HasSVAfterCorr");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(6, "DcaVtxDau");
    registry.get<TH1>(HIST("hVtx3BodyCounter"))->GetXaxis()->SetBinLabel(7, "CosPA");

    // Material correction in the DCA fitter
    o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
    if (useMatCorrType == 1) {
      LOGF(info, "TGeo correction requested, loading geometry");
      if (!o2::base::GeometryManager::isGeometryLoaded()) {
        ccdb->get<TGeoManager>(geoPath);
      }
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
    }
    if (useMatCorrType == 2) {
      LOGF(info, "LUT correction requested, loading LUT");
      lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(lutPath));
      matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
    }
    fitter3body.setMatCorrType(matCorr);
  }

  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (mRunNumber == bc.runNumber()) {
      return;
    }

    // In case override, don't proceed, please - no CCDB access required
    if (d_bz_input > -990) {
      d_bz = d_bz_input;
      fitter3body.setBz(d_bz);
      o2::parameters::GRPMagField grpmag;
      if (fabs(d_bz) > 1e-5) {
        grpmag.setL3Current(30000.f / (d_bz / 5.0f));
      }
      o2::base::Propagator::initFieldFromGRP(&grpmag);
      mRunNumber = bc.runNumber();
      return;
    }

    auto run3grp_timestamp = bc.timestamp();
    o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(grpPath, run3grp_timestamp);
    o2::parameters::GRPMagField* grpmag = 0x0;
    if (grpo) {
      o2::base::Propagator::initFieldFromGRP(grpo);
      // Fetch magnetic field from ccdb for current collision
      d_bz = grpo->getNominalL3Field();
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    } else {
      grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, run3grp_timestamp);
      if (!grpmag) {
        LOG(fatal) << "Got nullptr from CCDB for path " << grpmagPath << " of object GRPMagField and " << grpPath << " of object GRPObject for timestamp " << run3grp_timestamp;
      }
      o2::base::Propagator::initFieldFromGRP(grpmag);
      // Fetch magnetic field from ccdb for current collision
      d_bz = std::lround(5.f * grpmag->getL3Current() / 30000.f);
      LOG(info) << "Retrieved GRP for timestamp " << run3grp_timestamp << " with magnetic field of " << d_bz << " kZG";
    }
    mRunNumber = bc.runNumber();
    // Set magnetic field value once known
    fitter3body.setBz(d_bz);

    if (useMatCorrType == 2) {
      // setMatLUT only after magfield has been initalized
      // (setMatLUT has implicit and problematic init field call if not)
      o2::base::Propagator::Instance()->setMatLUT(lut);
    }
  }
  //------------------------------------------------------------------

  void process(aod::Collision const& collision, FullTracksExtIU const& tracks, aod::Decay3Bodys const& decay3bodys, aod::BCsWithTimestamps const&)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    initCCDB(bc);
    registry.fill(HIST("hEventCounter"), 0.5);

    for (auto& vtx3body : decay3bodys) {

      registry.fill(HIST("hVtx3BodyCounter"), 0.5);

      auto t0 = vtx3body.track0_as<FullTracksExtIU>();
      auto t1 = vtx3body.track1_as<FullTracksExtIU>();
      auto t2 = vtx3body.track2_as<FullTracksExtIU>();
      if (t0.collisionId() != t1.collisionId() || t0.collisionId() != t2.collisionId()) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), 1.5);

      if (t0.tpcNClsCrossedRows() < mincrossedrows && t1.tpcNClsCrossedRows() < mincrossedrows && t2.tpcNClsCrossedRows() < mincrossedrows) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), 2.5);

      auto Track0 = getTrackParCov(t0);
      auto Track1 = getTrackParCov(t1);
      auto Track2 = getTrackParCov(t2);
      int n3bodyVtx = fitter3body.process(Track0, Track1, Track2);
      if (n3bodyVtx == 0) { // discard this pair
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), 3.5);

      double finalXTrack0 = fitter3body.getTrack(0).getX();
      double finalXTrack1 = fitter3body.getTrack(1).getX();
      double finalXTrack2 = fitter3body.getTrack(2).getX();

      // Rotate to desired alpha
      Track0.rotateParam(fitter3body.getTrack(0).getAlpha());
      Track1.rotateParam(fitter3body.getTrack(1).getAlpha());
      Track2.rotateParam(fitter3body.getTrack(2).getAlpha());

      // Retry closer to minimum with material corrections
      o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;
      if (useMatCorrType == 1)
        matCorr = o2::base::Propagator::MatCorrType::USEMatCorrTGeo;
      if (useMatCorrType == 2)
        matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;

      o2::base::Propagator::Instance()->propagateToX(Track0, finalXTrack0, d_bz, maxSnp, maxStep, matCorr);
      o2::base::Propagator::Instance()->propagateToX(Track1, finalXTrack1, d_bz, maxSnp, maxStep, matCorr);
      o2::base::Propagator::Instance()->propagateToX(Track2, finalXTrack2, d_bz, maxSnp, maxStep, matCorr);

      n3bodyVtx = fitter3body.process(Track0, Track1, Track2);
      if (n3bodyVtx == 0) { // discard this pair
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), 4.5);

      std::array<float, 3> pos = {0.};
      const auto& vtxXYZ = fitter3body.getPCACandidate();
      for (int i = 0; i < 3; i++) {
        pos[i] = vtxXYZ[i];
      }

      std::array<float, 3> p0 = {0.}, p1 = {0.}, p2{0.};
      Track0.getPxPyPzGlo(p0);
      Track1.getPxPyPzGlo(p1);
      Track2.getPxPyPzGlo(p2);
      std::array<float, 3> p3B = {p0[0] + p1[0] + p2[0], p0[1] + p1[1] + p2[1], p0[2] + p1[2] + p2[2]};

      if (fitter3body.getChi2AtPCACandidate() > dcavtxdau) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), 5.5);

      float VtxcosPA = RecoDecay::cpa(array{collision.posX(), collision.posY(), collision.posZ()}, array{pos[0], pos[1], pos[2]}, array{p3B[0], p3B[1], p3B[2]});
      if (VtxcosPA < minCosPA3body) {
        continue;
      }
      registry.fill(HIST("hVtx3BodyCounter"), 6.5);

      vtx3bodydata(
        t0.globalIndex(), t1.globalIndex(), t2.globalIndex(), collision.globalIndex(), vtx3body.globalIndex(),
        pos[0], pos[1], pos[2],
        p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2],
        fitter3body.getChi2AtPCACandidate(),
        t0.dcaXY(), t1.dcaXY(), t2.dcaXY());
    }
  }
};

struct hypertriton3bodyDataLinkBuilder {
  Produces<aod::Decay3BodyDataLink> vtxdataLink;

  void init(InitContext const&) {}

  void process(aod::Decay3Bodys const& decay3bodytable, aod::Vtx3BodyDatas const& vtxdatatable)
  {
    std::vector<int> lIndices;
    lIndices.reserve(decay3bodytable.size());
    for (int ii = 0; ii < decay3bodytable.size(); ii++)
      lIndices[ii] = -1;
    for (auto& vtxdata : vtxdatatable) {
      lIndices[vtxdata.decay3bodyId()] = vtxdata.globalIndex();
    }
    for (int ii = 0; ii < decay3bodytable.size(); ii++) {
      vtxdataLink(lIndices[ii]);
    }
  }
};

struct hypertriton3bodyLabelBuilder {

  Produces<aod::McVtx3BodyLabels> vtxlabels;
  Produces<aod::McFullVtx3BodyLabels> vtxfulllabels;

  // for bookkeeping purposes: how many V0s come from same mother etc
  HistogramRegistry registry{
    "registry",
    {
      {"hLabelCounter", "hLabelCounter", {HistType::kTH1F, {{3, 0.0f, 3.0f}}}},
      {"hHypertritonCounter", "hHypertritonCounter", {HistType::kTH1F, {{4, 0.0f, 4.0f}}}},
      {"hPIDCounter", "hPIDCounter", {HistType::kTH1F, {{6, 0.0f, 6.0f}}}},
      {"hHypertritonMCPt", "hHypertritonMCPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
      {"hAntiHypertritonMCPt", "hAntiHypertritonMCPt", {HistType::kTH1F, {{100, 0.0f, 10.0f}}}},
      {"hHypertritonMCMass", "hHypertritonMCMass", {HistType::kTH1F, {{40, 2.95f, 3.05f}}}},
      {"hAntiHypertritonMCMass", "hAntiHypertritonMCMass", {HistType::kTH1F, {{40, 2.95f, 3.05f}}}},
      {"h3dTotalTrueHypertriton", "h3dTotalTrueHypertriton", {HistType::kTH3F, {{50, 0, 50, "ct(cm)"}, {200, 0.0f, 10.0f, "#it{p}_{T} (GeV/c)"}, {40, 2.95f, 3.05f, "Inv. Mass (GeV/c^{2})"}}}},
    },
  };

  void init(InitContext const&)
  {
    registry.get<TH1>(HIST("hLabelCounter"))->GetXaxis()->SetBinLabel(1, "Total");
    registry.get<TH1>(HIST("hLabelCounter"))->GetXaxis()->SetBinLabel(2, "Same MotherParticle");
    registry.get<TH1>(HIST("hLabelCounter"))->GetXaxis()->SetBinLabel(3, "True H3L");
    registry.get<TH1>(HIST("hHypertritonCounter"))->GetXaxis()->SetBinLabel(1, "H3L");
    registry.get<TH1>(HIST("hHypertritonCounter"))->GetXaxis()->SetBinLabel(2, "H3L daughters pass PID");
    registry.get<TH1>(HIST("hHypertritonCounter"))->GetXaxis()->SetBinLabel(3, "#bar{H3L}");
    registry.get<TH1>(HIST("hHypertritonCounter"))->GetXaxis()->SetBinLabel(4, "#bar{H3L} daughters pass PID");
    registry.get<TH1>(HIST("hPIDCounter"))->GetXaxis()->SetBinLabel(1, "H3L Proton PID > 5");
    registry.get<TH1>(HIST("hPIDCounter"))->GetXaxis()->SetBinLabel(2, "H3L Pion PID > 5");
    registry.get<TH1>(HIST("hPIDCounter"))->GetXaxis()->SetBinLabel(3, "H3L Deuteron PID > 5");
    registry.get<TH1>(HIST("hPIDCounter"))->GetXaxis()->SetBinLabel(4, "#bar{H3L} Proton PID > 5");
    registry.get<TH1>(HIST("hPIDCounter"))->GetXaxis()->SetBinLabel(5, "#bar{H3L} Pion PID > 5");
    registry.get<TH1>(HIST("hPIDCounter"))->GetXaxis()->SetBinLabel(6, "#bar{H3L} Deuteron PID > 5");
  }

  Configurable<float> TpcPidNsigmaCut{"TpcPidNsigmaCut", 5, "TpcPidNsigmaCut"};

  void processDoNotBuildLabels(aod::Collisions::iterator const& collision)
  {
    // dummy process function - should not be required in the future
  }
  PROCESS_SWITCH(hypertriton3bodyLabelBuilder, processDoNotBuildLabels, "Do not produce MC label tables", true);

  void processBuildLabels(aod::Collisions::iterator const& collision, aod::Decay3BodysLinked const& decay3bodys, aod::Vtx3BodyDatas const& vtx3bodydatas, LabeledTracks const&, aod::McParticles const& particlesMC)
  {
    std::vector<int> lIndices;
    lIndices.reserve(vtx3bodydatas.size());
    for (int ii = 0; ii < vtx3bodydatas.size(); ii++) {
      lIndices[ii] = -1;
    }

    for (auto& decay3body : decay3bodys) {

      int lLabel = -1;
      int lPDG = -1;
      float lPt = -1;
      double MClifetime = -1;
      bool is3bodyDecay = false;
      int lGlobalIndex = -1;

      auto lTrack0 = decay3body.track0_as<LabeledTracks>();
      auto lTrack1 = decay3body.track1_as<LabeledTracks>();
      auto lTrack2 = decay3body.track2_as<LabeledTracks>();
      registry.fill(HIST("hLabelCounter"), 0.5);

      // Association check
      // There might be smarter ways of doing this in the future
      if (!lTrack0.has_mcParticle() || !lTrack1.has_mcParticle() || !lTrack2.has_mcParticle()) {
        vtxfulllabels(-1);
        continue;
      }
      auto lMCTrack0 = lTrack0.mcParticle_as<aod::McParticles>();
      auto lMCTrack1 = lTrack1.mcParticle_as<aod::McParticles>();
      auto lMCTrack2 = lTrack2.mcParticle_as<aod::McParticles>();
      if (!lMCTrack0.has_mothers() || !lMCTrack1.has_mothers() || !lMCTrack2.has_mothers()) {
        vtxfulllabels(-1);
        continue;
      }

      for (auto& lMother0 : lMCTrack0.mothers_as<aod::McParticles>()) {
        for (auto& lMother1 : lMCTrack1.mothers_as<aod::McParticles>()) {
          for (auto& lMother2 : lMCTrack2.mothers_as<aod::McParticles>()) {
            if (lMother0.globalIndex() == lMother1.globalIndex() && lMother0.globalIndex() == lMother2.globalIndex()) {
              lGlobalIndex = lMother1.globalIndex();
              lPt = lMother1.pt();
              lPDG = lMother1.pdgCode();
              MClifetime = RecoDecay::sqrtSumOfSquares(lMCTrack2.vx() - lMother2.vx(), lMCTrack2.vy() - lMother2.vy(), lMCTrack2.vz() - lMother2.vz()) * o2::constants::physics::MassHyperTriton / lMother2.p();
              is3bodyDecay = true; // vtxs with the same mother
            }
          }
        }
      } // end association check
      if (!is3bodyDecay) {
        vtxfulllabels(-1);
        continue;
      }
      registry.fill(HIST("hLabelCounter"), 1.5);

      // Intended for cross-checks only
      // N.B. no rapidity cut!
      if (lPDG == 1010010030 && lMCTrack0.pdgCode() == 2212 && lMCTrack1.pdgCode() == -211 && lMCTrack2.pdgCode() == 1000010020) {
        lLabel = lGlobalIndex;
        double hypertritonMCMass = RecoDecay::m(array{array{lMCTrack0.px(), lMCTrack0.py(), lMCTrack0.pz()}, array{lMCTrack1.px(), lMCTrack1.py(), lMCTrack1.pz()}, array{lMCTrack2.px(), lMCTrack2.py(), lMCTrack2.pz()}}, array{o2::constants::physics::MassProton, o2::constants::physics::MassPionCharged, o2::constants::physics::MassDeuteron});
        registry.fill(HIST("hLabelCounter"), 2.5);
        registry.fill(HIST("hHypertritonCounter"), 0.5);
        registry.fill(HIST("hHypertritonMCPt"), lPt);
        registry.fill(HIST("hHypertritonMCMass"), hypertritonMCMass);
        registry.fill(HIST("h3dTotalTrueHypertriton"), MClifetime, lPt, hypertritonMCMass);
        if (TMath::Abs(lTrack0.tpcNSigmaPr()) > TpcPidNsigmaCut) {
          registry.fill(HIST("hPIDCounter"), 0.5);
        }
        if (TMath::Abs(lTrack1.tpcNSigmaPi()) > TpcPidNsigmaCut) {
          registry.fill(HIST("hPIDCounter"), 1.5);
        }
        if (TMath::Abs(lTrack2.tpcNSigmaDe()) > TpcPidNsigmaCut) {
          registry.fill(HIST("hPIDCounter"), 2.5);
        }
        if (TMath::Abs(lTrack0.tpcNSigmaPr()) < TpcPidNsigmaCut && TMath::Abs(lTrack1.tpcNSigmaPi()) < TpcPidNsigmaCut && TMath::Abs(lTrack2.tpcNSigmaDe()) < TpcPidNsigmaCut) {
          registry.fill(HIST("hHypertritonCounter"), 1.5);
        }
      }
      if (lPDG == -1010010030 && lMCTrack0.pdgCode() == 211 && lMCTrack1.pdgCode() == -2212 && lMCTrack2.pdgCode() == -1000010020) {
        lLabel = lGlobalIndex;
        double antiHypertritonMCMass = RecoDecay::m(array{array{lMCTrack0.px(), lMCTrack0.py(), lMCTrack0.pz()}, array{lMCTrack1.px(), lMCTrack1.py(), lMCTrack1.pz()}, array{lMCTrack2.px(), lMCTrack2.py(), lMCTrack2.pz()}}, array{o2::constants::physics::MassPionCharged, o2::constants::physics::MassProton, o2::constants::physics::MassDeuteron});
        registry.fill(HIST("hLabelCounter"), 2.5);
        registry.fill(HIST("hHypertritonCounter"), 2.5);
        registry.fill(HIST("hAntiHypertritonMCPt"), lPt);
        registry.fill(HIST("hAntiHypertritonMCMass"), antiHypertritonMCMass);
        registry.fill(HIST("h3dTotalTrueHypertriton"), MClifetime, lPt, antiHypertritonMCMass);
        if (TMath::Abs(lTrack0.tpcNSigmaPi()) > TpcPidNsigmaCut) {
          registry.fill(HIST("hPIDCounter"), 4.5);
        }
        if (TMath::Abs(lTrack1.tpcNSigmaPr()) > TpcPidNsigmaCut) {
          registry.fill(HIST("hPIDCounter"), 3.5);
        }
        if (TMath::Abs(lTrack2.tpcNSigmaDe()) > TpcPidNsigmaCut) {
          registry.fill(HIST("hPIDCounter"), 5.5);
        }
        if (TMath::Abs(lTrack0.tpcNSigmaPi()) < TpcPidNsigmaCut && TMath::Abs(lTrack1.tpcNSigmaPr()) < TpcPidNsigmaCut && TMath::Abs(lTrack2.tpcNSigmaDe()) < TpcPidNsigmaCut) {
          registry.fill(HIST("hHypertritonCounter"), 3.5);
        }
      }

      // Construct label table, only true hypertriton and true daughters with a specified order is labeled
      // for matter: track0->p, track1->pi, track2->d
      // for antimatter: track0->pi, track1->p, track2->d
      vtxfulllabels(lLabel);
      if (decay3body.vtx3BodyDataId() != -1) {
        lIndices[decay3body.vtx3BodyDataId()] = lLabel;
      }
    }
    for (int ii = 0; ii < vtx3bodydatas.size(); ii++) {
      vtxlabels(lIndices[ii]);
    }
  }
  PROCESS_SWITCH(hypertriton3bodyLabelBuilder, processBuildLabels, "Produce MC label tables", false);
};

struct hypertriton3bodyInitializer {
  Spawns<aod::Vtx3BodyDatas> vtx3bodydatas;
  void init(InitContext const&) {}
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<hypertriton3bodyBuilder>(cfgc),
    adaptAnalysisTask<hypertriton3bodyDataLinkBuilder>(cfgc),
    adaptAnalysisTask<hypertriton3bodyLabelBuilder>(cfgc),
    adaptAnalysisTask<hypertriton3bodyInitializer>(cfgc),
  };
}
