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
#include "MetadataHelper.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

MetadataHelper metadataInfo;

static constexpr int nCases = 2;
static constexpr int nParameters = 12;
static const std::vector<std::string> casesNames{"Data", "MC"};
static const std::vector<std::string> parameterNames{"RespITSPar1", "RespITSPar2", "RespITSPar3",
                                                     "RespITSPar1_Z2", "RespITSPar2_Z2", "RespITSPar3_Z2",
                                                     "ResolutionPar1", "ResolutionPar2", "ResolutionPar3",
                                                     "ResolutionPar1_Z2", "ResolutionPar2_Z2", "ResolutionPar3_Z2"};

static constexpr float defaultParameters[nCases][nParameters] = {
  {1.18941, 1.53792, 1.69961, 2.35117, 1.80347, 5.14355, 1.94669e-01, -2.08616e-01, 1.30753, 8.74371e-02, -1.82804, 5.06449e-01},
  {1.18941, 1.53792, 1.69961, 2.35117, 1.80347, 5.14355, 1.94669e-01, -2.08616e-01, 1.30753, 8.74371e-02, -1.82804, 5.06449e-01}};

/// Task to produce the ITS PID information for each particle species
/// The parametrization is: [p0/(bg)**p1 + p2] being bg = p/m. Different parametrizations are used for He3 and Alpha particles.
/// The resolution depends on the bg and is modelled with an erf function: p0*TMath::Erf((bg-p1)/p2)
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
      const char* key = metadataInfo.isMC() ? "MC" : "Data";
      o2::aod::ITSResponse::setParameters(itsParams->get(key, "RespITSPar1"),
                                          itsParams->get(key, "RespITSPar2"),
                                          itsParams->get(key, "RespITSPar3"),
                                          itsParams->get(key, "RespITSPar1_Z2"),
                                          itsParams->get(key, "RespITSPar2_Z2"),
                                          itsParams->get(key, "RespITSPar3_Z2"),
                                          itsParams->get(key, "ResolutionPar1"),
                                          itsParams->get(key, "ResolutionPar2"),
                                          itsParams->get(key, "ResolutionPar3"),
                                          itsParams->get(key, "ResolutionPar1_Z2"),
                                          itsParams->get(key, "ResolutionPar2_Z2"),
                                          itsParams->get(key, "ResolutionPar3_Z2"));
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
