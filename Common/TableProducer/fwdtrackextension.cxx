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
// \file fwdtrackextension.cxx
// \brief Task performing forward track DCA computation.
// \author Maurice Coquet, maurice.louis.coquet@cern.ch
//

#include "Common/Core/fwdtrackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <GlobalTracking/MatchGlobalFwd.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/DataTypes.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/TrackFwd.h>

#include <Math/MatrixRepresentationsStatic.h>
#include <Math/SMatrix.h>

#include <vector>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

using MuonsWithCov = soa::Join<aod::FwdTracks, aod::FwdTracksCov>;

struct FwdTrackExtension {
  Produces<aod::FwdTracksDCA> fwdDCA;
  Configurable<std::string> fGeoPath{"fGeoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};
  Configurable<std::string> fGrpmagPath{"fGrpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> fConfigCcdbUrl{"fConfigCcdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<bool> fRefitGlobalMuon{"fRefitGlobalMuon", true, "Recompute parameters of global muons"};

  Service<o2::ccdb::BasicCCDBManager> fCCDB;
  o2::parameters::GRPMagField* grpmag = nullptr; // for run 3, we access GRPMagField from GLO/Config/GRPMagField
  int fCurrentRun;                               // needed to detect if the run changed and trigger update of magnetic field

  void init(o2::framework::InitContext&)
  {
    // Load geometry
    fCCDB->setURL(fConfigCcdbUrl);
    fCCDB->setCaching(true);
    fCCDB->setLocalObjectValidityChecking();

    if (!o2::base::GeometryManager::isGeometryLoaded()) {
      LOGF(info, "Load geometry from CCDB");
      fCCDB->get<TGeoManager>(fGeoPath);
    }
  }

  void process(aod::Collisions::iterator const& collision, o2::aod::BCsWithTimestamps const& /*...*/, MuonsWithCov const& tracks, aod::MFTTracks const& /*...*/)
  {
    auto bc = collision.template bc_as<o2::aod::BCsWithTimestamps>();
    if (fCurrentRun != bc.runNumber()) {
      grpmag = fCCDB->getForTimeStamp<o2::parameters::GRPMagField>(fGrpmagPath, bc.timestamp());
      if (grpmag != nullptr) {
        LOGF(info, "Init field from GRP");
        o2::base::Propagator::initFieldFromGRP(grpmag);
      }
      LOGF(info, "Set field for muons");
      o2::mch::TrackExtrap::setField();
      fCurrentRun = bc.runNumber();
    }
    const float zField = grpmag->getNominalL3Field();
    for (const auto& track : tracks) {
      const auto trackType = track.trackType();
      o2::dataformats::GlobalFwdTrack fwdtrack = o2::aod::fwdtrackutils::getTrackParCovFwd(track, track);
      if (fRefitGlobalMuon && (trackType == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack || trackType == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalForwardTrack)) {
        auto muontrack = track.template matchMCHTrack_as<MuonsWithCov>();
        auto mfttrack = track.template matchMFTTrack_as<aod::MFTTracks>();
        o2::dataformats::GlobalFwdTrack propmuon = o2::aod::fwdtrackutils::propagateMuon(muontrack, muontrack, collision, o2::aod::fwdtrackutils::propagationPoint::kToVertex, 0.f, zField);
        SMatrix5 tpars(mfttrack.x(), mfttrack.y(), mfttrack.phi(), mfttrack.tgl(), mfttrack.signed1Pt());
        SMatrix55 tcovs{};
        o2::track::TrackParCovFwd mft{mfttrack.z(), tpars, tcovs, mfttrack.chi2()};
        fwdtrack = o2::aod::fwdtrackutils::refitGlobalMuonCov(propmuon, mft);
      }
      const auto proptrack = o2::aod::fwdtrackutils::propagateTrackParCovFwd(fwdtrack, trackType, collision, o2::aod::fwdtrackutils::propagationPoint::kToDCA, 0.f, zField);
      const float dcaX = (proptrack.getX() - collision.posX());
      const float dcaY = (proptrack.getY() - collision.posY());
      fwdDCA(dcaX, dcaY);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FwdTrackExtension>(cfgc)};
  return workflow;
}
