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
/// \author Giorgio Alberto Lucia giorgio.alberto.lucia@cern.ch
/// \brief  Task to produce PID tables for ITS split for each particle.
///         Only the tables for the mass hypotheses requested are filled, the others are sent empty.
///

#include "Common/Core/MetadataHelper.h"
#include "Common/DataModel/PIDResponseITS.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/ASoA.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Array2D.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/runDataProcessing.h>

#include <chrono>
#include <cstdint>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

o2::common::core::MetadataHelper metadataInfo;

static constexpr int nCases = 2;
static constexpr int nParameters = 12;
static const std::vector<std::string> casesNames{"Data", "MC"};
static const std::vector<std::string> parameterNames{"RespITSPar1", "RespITSPar2", "RespITSPar3",
                                                     "RespITSPar1_Z2", "RespITSPar2_Z2", "RespITSPar3_Z2",
                                                     "ResolutionPar1", "ResolutionPar2", "ResolutionPar3",
                                                     "ResolutionPar1_Z2", "ResolutionPar2_Z2", "ResolutionPar3_Z2"};

static constexpr float defaultParameters[nCases][nParameters] = {
  {1.18941, 1.53792, 1.69961, 2.35117, 1.80347, 5.14355, 1.94669e-01, -2.08616e-01, 1.30753, 0.09, -999., -999.},
  {1.63806, 1.58847, 2.52275, 2.66505, 1.48405, 6.90453, 1.40487e-01, -4.31078e-01, 1.50052, 0.09, -999., -999.}};

/// Task to produce the ITS PID information for each particle species
/// The parametrization is: [p0/(bg)**p1 + p2] being bg = p/m. Different parametrizations are used for He3 and Alpha particles.
/// The resolution depends on the bg and is modelled with an erf function: p0*TMath::Erf((bg-p1)/p2). If p1/p2 is -999, the resolution is set to p0.
struct itsPid {

  Configurable<LabeledArray<float>> itsParams{"itsParams",
                                              {defaultParameters[0], nCases, nParameters, casesNames, parameterNames},
                                              "Response parameters"};
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
      LOG(fatal) << "Not implemented yet";
    } else {
      const char* dataType = metadataInfo.isMC() ? "MC" : "Data";
      o2::aod::ITSResponse::setParameters(itsParams->get(dataType, "RespITSPar1"),
                                          itsParams->get(dataType, "RespITSPar2"),
                                          itsParams->get(dataType, "RespITSPar3"),
                                          itsParams->get(dataType, "RespITSPar1_Z2"),
                                          itsParams->get(dataType, "RespITSPar2_Z2"),
                                          itsParams->get(dataType, "RespITSPar3_Z2"),
                                          itsParams->get(dataType, "ResolutionPar1"),
                                          itsParams->get(dataType, "ResolutionPar2"),
                                          itsParams->get(dataType, "ResolutionPar3"),
                                          itsParams->get(dataType, "ResolutionPar1_Z2"),
                                          itsParams->get(dataType, "ResolutionPar2_Z2"),
                                          itsParams->get(dataType, "ResolutionPar3_Z2"));
    }
  }

  /// Dummy process function for BCs, needed in case both Run2 and Run3 process functions are disabled
  void process(aod::Timestamps const&) {}

  void processTest(o2::soa::Join<aod::TracksIU, aod::TracksExtra> const& tracks)
  {
    auto tracksWithPid = soa::Attach<o2::soa::Join<aod::TracksIU, aod::TracksExtra>,
                                     aod::pidits::ITSNSigmaEl, aod::pidits::ITSNSigmaMu, aod::pidits::ITSNSigmaPi,
                                     aod::pidits::ITSNSigmaKa, aod::pidits::ITSNSigmaPr, aod::pidits::ITSNSigmaDe,
                                     aod::pidits::ITSNSigmaTr, aod::pidits::ITSNSigmaHe, aod::pidits::ITSNSigmaAl>(tracks);

    for (const auto& track : tracksWithPid) {
      LOG(info) << track.itsNSigmaEl();
      LOG(info) << track.itsNSigmaPi();
      LOG(info) << track.itsNSigmaPr();
    }
  }
  PROCESS_SWITCH(itsPid, processTest, "Produce a test", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  // Parse the metadata
  metadataInfo.initMetadata(cfgc);
  auto workflow = WorkflowSpec{adaptAnalysisTask<itsPid>(cfgc)};
  return workflow;
}
