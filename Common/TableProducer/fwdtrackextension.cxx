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
// Task performing forward track DCA computation
//

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"

#include "Math/SMatrix.h"
#include "ReconstructionDataFormats/TrackFwd.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<double, 5>;

struct FwdTrackExtension {

  Produces<aod::FwdTracksDCA> extendedTrackQuantities;

  float mMagField = 0.0;

  Configurable<bool> fUseRemoteField{"cfgUseRemoteField", true, "Chose whether to fetch the magnetic field from ccdb or set it manually"};
  Configurable<float> fConfigMagField{"cfgMagField", 5.0f, "Manually set magnetic field"};
  Configurable<std::string> ccdburl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::parameters::GRPMagField* grpmag = nullptr;

  int fCurrentRun;

  void init(o2::framework::InitContext& context)
  {
    fCurrentRun = 0;

    ccdb->setURL(ccdburl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
  }

  void process(aod::FwdTracks const& tracks, aod::BCsWithTimestamps const& bcs, aod::Collisions const&)
  {
    for (auto& track : tracks) {
      float dcaX = -999;
      float dcaY = -999;
      if (track.has_collision()) {
        if (track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalMuonTrack || track.trackType() == o2::aod::fwdtrack::ForwardTrackTypeEnum::GlobalForwardTrack) {

          auto const& collision = track.collision();
          auto bc = collision.bc_as<aod::BCsWithTimestamps>();

          if (fCurrentRun != bc.runNumber()) {
            if (fUseRemoteField.value) {
              grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
              if (grpmag != nullptr) {
                mMagField = grpmag->getNominalL3Field();
              } else {
                LOGF(fatal, "GRP object is not available in CCDB at timestamp=%llu", bc.timestamp());
              }
            } else {
              mMagField = fConfigMagField.value;
            }
            fCurrentRun = bc.runNumber();
          }

          double chi2 = track.chi2();
          SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
          std::vector<double> v1;
          SMatrix55 tcovs(v1.begin(), v1.end());
          o2::track::TrackParCovFwd pars1{track.z(), tpars, tcovs, chi2};

          pars1.propagateToZhelix(collision.posZ(), mMagField);

          dcaX = (pars1.getX() - collision.posX());
          dcaY = (pars1.getY() - collision.posY());
        }
      }
      extendedTrackQuantities(dcaX, dcaY);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{adaptAnalysisTask<FwdTrackExtension>(cfgc)};
  return workflow;
}
