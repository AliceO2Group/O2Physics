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

///
/// \file   pidITS.cxx
/// \since  2024-11-12
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \author Francesco Mazzaschi francesco.mazzaschi@cern.ch
/// \brief  Task to produce PID tables for ITS split for each particle.
///         Only the tables for the mass hypotheses requested are filled, the others are sent empty.
///

#include <utility>
#include <vector>
#include <string>

// O2 includes
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "ReconstructionDataFormats/Track.h"
#include "CCDB/BasicCCDBManager.h"
#include "TOFBase/EventTimeMaker.h"

// O2Physics includes
#include "Common/DataModel/PIDResponseITS.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

/// Task to produce the ITS PID information for each particle species
/// The parametrization is: [p0/(bg)**p1 + p2] * pow(q, p3), being bg = p/m and q the charge 
struct itsPid {

  Configurable<float> bb1{"bb1", 0.f, "Bethe Bloch parameter 1"};
  Configurable<float> bb2{"bb2", 0.f, "Bethe Bloch parameter 2"};
  Configurable<float> bb3{"bb3", 0.f, "Bethe Bloch parameter 3"};
  Configurable<float> chargeExponent{"chargeExponent", 0.f, "Charge exponent"};
  Configurable<float> resolution{"resolution", 0.f, "Charge exponent"};
  Configurable<bool> getFromCCDB{"getFromCCDB", false, "Get the parameters from CCDB"};

  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> paramfile{"param-file", "", "Path to the parametrization object, if empty the parametrization is not taken from file"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdbPath", "Analysis/PID/TPC/Response", "Path of the TPC parametrization on the CCDB"};
  Configurable<std::string> recoPass{"recoPass", "", "Reconstruction pass name for CCDB query (automatically takes latest object for timestamp if blank)"};
  Configurable<int64_t> ccdbTimestamp{"ccdb-timestamp", 0, "timestamp of the object used to query in CCDB the detector response. Exceptions: -1 gets the latest object, 0 gets the run dependent timestamp"};

  void init(o2::framework::InitContext&)
  {
    if (getFromCCDB) {
      ccdb->setURL(url.value);
      ccdb->setCaching(true);
      ccdb->setLocalObjectValidityChecking();
      ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    } else {
      o2::aod::mITSResponse.mITSRespParams[0] = bb1;
      o2::aod::mITSResponse.mITSRespParams[1] = bb2;
      o2::aod::mITSResponse.mITSRespParams[2] = bb3;
      o2::aod::mITSResponse.mChargeFactor = chargeExponent;
      o2::aod::mITSResponse.mResolution = resolution;
    }
  }

  /// Dummy process function for BCs, needed in case both Run2 and Run3 process functions are disabled
  void process(aod::BCs const&) {}

  void processTest(o2::soa::Join<aod::TracksIU, aod::TracksExtra> const& tracks)
  {
    auto tracksWithPid = soa::Attach<o2::soa::Join<aod::TracksIU, aod::TracksExtra>,
                                     aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi,
                                     aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe,
                                     aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe, aod::pidits::ITSNSigmaAl>(tracks);

    for (const auto& track : tracksWithPid) {
      LOG(info) << track.itsNSigmaPr();
    }
  }
  PROCESS_SWITCH(itsPid, processTest, "Produce a test", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<itsPid>(cfgc)};
  return workflow;
}
