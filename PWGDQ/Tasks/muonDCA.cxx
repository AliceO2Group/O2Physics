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
//
/// \file muonDCA.cxx
/// \brief Task to compute and evaluate DCA quantities
/// \author Nicolas Bizé <nicolas.bize@cern.ch>, SUBATECH
//
#include "PWGDQ/Core/VarManager.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "GlobalTracking/MatchGlobalFwd.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::aod;

using MyMuons = soa::Join<aod::ReducedMuons, aod::ReducedMuonsExtra, aod::ReducedMuonsCov>;
using MyEvents = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended>;
using MyEventsVtxCov = soa::Join<aod::ReducedEvents, aod::ReducedEventsExtended, aod::ReducedEventsVtxCov>;

// constexpr static uint32_t gkMuonDCAFillMapWithCov = VarManager::ObjTypes::ReducedMuon | VarManager::ObjTypes::ReducedMuonExtra | VarManager::ObjTypes::ReducedMuonCov | VarManager::ObjTypes::MuonDCA;

constexpr static int toVertex = VarManager::kToVertex;
constexpr static int toDCA = VarManager::kToDCA;
constexpr static int toRabs = VarManager::kToRabs;

static o2::globaltracking::MatchGlobalFwd mExtrap;
template <typename T>
bool isSelected(const T& muon);

struct muonExtrap {
  Produces<ReducedMuonsDca> dcaTable;
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> fConfigCcdbUrl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::parameters::GRPMagField* grpmag = nullptr; // for run 3, we access GRPMagField from GLO/Config/GRPMagField
  int fCurrentRun;                               // needed to detect if the run changed and trigger update of magnetic field

  HistogramRegistry registry{
    "registry",
    {}};

  void init(o2::framework::InitContext&)
  {
    // Load geometry
    fCCDB->setURL(fConfigCcdbUrl);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      LOGF(info, "Load geometry from CCDB");
      fCCDB->get<TGeoManager>(geoPath);
    }

    AxisSpec pdcaAxis = {5000, 0.0, 5000.0, "p #times DCA"};
    AxisSpec dcaAxis = {200, 0.0, 200.0, "DCA"};
    AxisSpec dcaxAxis = {200, -100.0, 100.0, "DCA_x"};
    AxisSpec dcayAxis = {200, -100.0, 100.0, "DCA_y"};
    AxisSpec rabsAxis = {100, 0., 100.0, "R_{abs}"};
    AxisSpec xAxis = {200, -100., 100.0, "x (cm)"};
    AxisSpec yAxis = {200, -100., 100.0, "y (cm)"};
    AxisSpec zAxis = {200, -100., 100.0, "z (cm)"};

    HistogramConfigSpec pdcaSpec({HistType::kTH1F, {pdcaAxis}});
    HistogramConfigSpec dcaSpec({HistType::kTH1F, {dcaAxis}});
    HistogramConfigSpec dcaxSpec({HistType::kTH1F, {dcaxAxis}});
    HistogramConfigSpec dcaySpec({HistType::kTH1F, {dcayAxis}});
    HistogramConfigSpec rabsSpec({HistType::kTH1F, {rabsAxis}});
    HistogramConfigSpec xSpec({HistType::kTH1F, {xAxis}});
    HistogramConfigSpec ySpec({HistType::kTH1F, {yAxis}});
    HistogramConfigSpec zSpec({HistType::kTH1F, {zAxis}});

    registry.add("pdca", "pDCA", pdcaSpec);
    registry.add("dca", "DCA", dcaSpec);
    registry.add("dcax", "DCA_x", dcaxSpec);
    registry.add("dcay", "DCA_y", dcaySpec);
    registry.add("rabs", "R_{abs}", rabsSpec);
    registry.add("xAtVtx", "x at vertex", xSpec);
    registry.add("xAtDCA", "x at DCA", xSpec);
    registry.add("xAtRabs", "x at end abs", xSpec);
    registry.add("yAtVtx", "y at vertex", ySpec);
    registry.add("yAtDCA", "y at DCA", ySpec);
    registry.add("yAtRabs", "y at end abs", ySpec);
    registry.add("zAtVtx", "z at vertex", zSpec);
    registry.add("zAtDCA", "z at DCA", zSpec);
    registry.add("zAtRabs", "z at end abs", zSpec);
  }

  void processExtrapolation(MyEventsVtxCov::iterator const& collision, MyMuons const& muons)
  {
    if (fCurrentRun != collision.runNumber()) {
      grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, collision.timestamp());
      if (grpmag != nullptr) {
        LOGF(info, "Init field from GRP");
        o2::base::Propagator::initFieldFromGRP(grpmag);
      }
      LOGF(info, "Set field for muons");
      VarManager::SetupMuonMagField();
      fCurrentRun = collision.runNumber();
    }

    for (auto& muon : muons) {
      if (static_cast<int>(muon.trackType()) < 2) {
        continue; // Make sure to remove global muon tracks
      }
      // propagate muon track to vertex
      o2::dataformats::GlobalFwdTrack muonTrackAtVertex = VarManager::PropagateMuon(muon, collision, toVertex);

      // propagate muon track to DCA
      o2::dataformats::GlobalFwdTrack muonTrackAtDCA = VarManager::PropagateMuon(muon, collision, toDCA);

      // propagate to Rabs
      o2::dataformats::GlobalFwdTrack muonTrackAtRabs = VarManager::PropagateMuon(muon, collision, toRabs);

      // Calculate DCA quantities (preferable to do it with VarManager)
      double dcax = muonTrackAtDCA.getX() - collision.posX();
      double dcay = muonTrackAtDCA.getY() - collision.posY();
      double dca = std::sqrt(dcax * dcax + dcay * dcay);
      double pdca = muonTrackAtVertex.getP() * dca;
      double xAtVtx = muonTrackAtVertex.getX();
      double yAtVtx = muonTrackAtVertex.getY();
      double zAtVtx = muonTrackAtVertex.getZ();
      double xAtDCA = muonTrackAtDCA.getX();
      double yAtDCA = muonTrackAtDCA.getY();
      double zAtDCA = muonTrackAtDCA.getZ();
      double xAbs = muonTrackAtRabs.getX();
      double yAbs = muonTrackAtRabs.getY();
      double zAbs = muonTrackAtRabs.getZ();

      double rabs = std::sqrt(xAbs * xAbs + yAbs * yAbs);

      // QA histograms
      registry.get<TH1>(HIST("pdca"))->Fill(pdca);
      registry.get<TH1>(HIST("dca"))->Fill(dca);
      registry.get<TH1>(HIST("dcax"))->Fill(dcax);
      registry.get<TH1>(HIST("dcay"))->Fill(dcay);
      registry.get<TH1>(HIST("rabs"))->Fill(rabs);

      registry.get<TH1>(HIST("xAtDCA"))->Fill(xAtDCA);
      registry.get<TH1>(HIST("xAtRabs"))->Fill(xAbs);
      registry.get<TH1>(HIST("xAtVtx"))->Fill(xAtVtx);

      registry.get<TH1>(HIST("yAtDCA"))->Fill(yAtDCA);
      registry.get<TH1>(HIST("yAtRabs"))->Fill(yAbs);
      registry.get<TH1>(HIST("yAtVtx"))->Fill(yAtVtx);

      registry.get<TH1>(HIST("zAtDCA"))->Fill(zAtDCA);
      registry.get<TH1>(HIST("zAtRabs"))->Fill(zAbs);
      registry.get<TH1>(HIST("zAtVtx"))->Fill(zAtVtx);

      // Fill DCA table
      dcaTable(pdca,
               dca,
               dcax,
               dcay,
               rabs,
               muonTrackAtVertex.getPt(),
               muonTrackAtVertex.getEta(),
               muonTrackAtVertex.getPhi(),
               muon.sign(),
               muon.isAmbiguous(),
               muonTrackAtVertex.getPx(),
               muonTrackAtVertex.getPy(),
               muonTrackAtVertex.getPz());
    }
  }

  PROCESS_SWITCH(muonExtrap, processExtrapolation, "process extrapolation", false);

  void processDummy(MyEventsVtxCov&)
  {
    // do nothing
  }

  PROCESS_SWITCH(muonExtrap, processDummy, "do nothing", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<muonExtrap>(cfgc)};
};
