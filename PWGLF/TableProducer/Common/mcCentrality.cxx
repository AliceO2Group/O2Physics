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
/// \file   mcCentrality.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \author Francesca Ercolessi francesca.ercolessi@cern.ch
/// \since  2024-06-05
/// \brief  Task to produce the table for the equalized multiplicity into centrality bins
///

#include "PWGLF/DataModel/mcCentrality.h"

#include "PWGLF/Utils/inelGt.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TH1.h>
#include <TFile.h>

#include <chrono>
#include <cstdint>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

/// Task to produce the response table
struct McCentrality {

  // Tables to produce
  Produces<aod::McCentFT0Ms> centFT0M;
  Produces<aod::McCentFT0As> centFT0A;
  Produces<aod::McCentFT0Cs> centFT0C;
  
  // NOTE: Commented out unused produces to prevent garbage columns in AOD
  // Produces<aod::McCentFV0As> centFV0A;
  // Produces<aod::McCentFDDMs> centFDDM;
  // Produces<aod::McCentNTPVs> centNTPV;

  // Input parameters
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdbUrl{"ccdbUrl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> ccdbTimestamp{"ccdbTimestamp", -1, "timestamp of the object used to query in CCDB the detector response. If 0 the object corresponding to the run number is used, if < 0 the latest object is used"};
  Configurable<std::string> path{"path", "/tmp/InputCalibMC.root", "path to calib file or ccdb path if begins with ccdb://"};
  Configurable<bool> selectPrimaries{"selectPrimaries", true, "Select only primary particles"};
  Service<o2::framework::O2DatabasePDG> pdgDB;
  ConfigurableAxis binsPercentile{"binsPercentile", {VARIABLE_WIDTH, 0, 0.001, 0.01, 1.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0}, "Binning of the percentile axis"};
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {1000, 0, 5000}, "Binning of the multiplicity axis"};
  
  // Added fillFt0M to allow running strictly on A/C without crashing on missing M
  Configurable<bool> fillFt0M{"fillFt0M", true, "Fills the FT0M histogram"};
  Configurable<bool> fillFt0A{"fillFt0A", false, "Fills the FT0A histogram"};
  Configurable<bool> fillFt0C{"fillFt0C", false, "Fills the FT0C histogram"};
  Configurable<bool> doNotCrashOnNull{"doNotCrashOnNull", false, "If ccdb object does not exist, fill with dummy values"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  TH1F* h1dFT0M = nullptr;
  TH1F* h1dFT0A = nullptr;
  TH1F* h1dFT0C = nullptr;

  o2::pwglf::ParticleCounter<o2::framework::O2DatabasePDG> mCounter;

  void init(o2::framework::InitContext& /*initContext*/)
  {
    // Set up the CCDB
    ccdb->setURL(ccdbUrl.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdb->setFatalWhenNull(false);

    mCounter.mPdgDatabase = pdgDB.service;
    mCounter.mSelectPrimaries = selectPrimaries.value;

    if (fillFt0M) {
      histos.add("FT0M/percentile", "FT0M percentile.", HistType::kTH1D, {{binsPercentile, "FT0M percentile"}});
      histos.add("FT0M/percentilevsMult", "FT0M percentile.", HistType::kTH2D, {{binsPercentile, "FT0M percentile"}, {binsMultiplicity, "FT0M mult."}});
    }
    if (fillFt0A) {
      histos.add("FT0A/percentile", "FT0A percentile.", HistType::kTH1D, {{binsPercentile, "FT0A percentile"}});
      histos.add("FT0A/percentilevsMult", "FT0A percentile.", HistType::kTH2D, {{binsPercentile, "FT0A percentile"}, {binsMultiplicity, "FT0A mult."}});
    }
    if (fillFt0C) {
      histos.add("FT0C/percentile", "FT0C percentile.", HistType::kTH1D, {{binsPercentile, "FT0C percentile"}});
      histos.add("FT0C/percentilevsMult", "FT0C percentile.", HistType::kTH2D, {{binsPercentile, "FT0C percentile"}, {binsMultiplicity, "FT0C mult."}});
    }

    TList* lOfInput = nullptr;
    TFile* f = nullptr; // Track file to close later

    if (path.value.rfind("ccdb://", 0) == 0) { // Getting post calib. from CCDB
      path.value.replace(0, 7, "");
      lOfInput = ccdb->get<TList>(path);
      if (!lOfInput) {
        if (doNotCrashOnNull) {
          LOG(info) << "Could not find the calibration TList from CCDB in path " << path << ", will fill tables with dummy values";
        } else {
          LOG(fatal) << "Could not find the calibration TList from CCDB in path " << path;
          return;
        }
      }
    } else { // Getting post calib. from file
      f = TFile::Open(path.value.c_str(), "READ");
      if (!f || !f->IsOpen()) {
        LOG(fatal) << "The input file " << path << " is not valid or open";
      }
      lOfInput = static_cast<TList*>(f->Get("ccdb_object"));
      if (!lOfInput) {
        f->ls();
        LOG(fatal) << "The input file " << path.value << " does not contain the TList ccdb_object";
      }
    }

    auto getHist = [this, lOfInput](const char* name) -> TH1F* {
      if (!lOfInput) {
        return nullptr;
      }
      auto obj = lOfInput->FindObject(name);
      if (!obj) {
        if (this->doNotCrashOnNull) {
          LOG(info) << "Could not open histogram " << name << " from TList, will fill tables with dummy values";
        } else {
          LOG(fatal) << "Could not open histogram " << name << " from TList";
        }
        return nullptr;
      }
      // Clone to detach from TFile directory so we can safely close the file
      auto hist = static_cast<TH1F*>(obj->Clone(Form("%s_clone", name)));
      hist->SetDirectory(nullptr); 
      return hist;
    };

    if (fillFt0M) {
      h1dFT0M = getHist("h1dFT0M");
    }
    if (fillFt0A) {
      h1dFT0A = getHist("h1dFT0A");
    }
    if (fillFt0C) {
      h1dFT0C = getHist("h1dFT0C");
    }

    // Safely close the file to prevent memory leaks
    if (f) {
      f->Close();
      delete f;
    }
  }

  void process(aod::McCollision const& /*mcCollision*/,
               aod::McParticles const& mcParticles)
  {
    const float nFT0A = mCounter.countFT0A(mcParticles);
    const float nFT0C = mCounter.countFT0C(mcParticles);
    const float nFT0M = nFT0A + nFT0C;

    if (fillFt0M) {
      float valueCentFT0M = 105.0f;
      if (h1dFT0M) {
        valueCentFT0M = h1dFT0M->GetBinContent(h1dFT0M->FindBin(nFT0M));
      }
      centFT0M(valueCentFT0M);
      histos.fill(HIST("FT0M/percentile"), valueCentFT0M);
      histos.fill(HIST("FT0M/percentilevsMult"), valueCentFT0M, nFT0M);
    }

    if (fillFt0A) {
      float valueCentFT0A = 105.0f;
      if (h1dFT0A) {
        valueCentFT0A = h1dFT0A->GetBinContent(h1dFT0A->FindBin(nFT0A));
      }
      centFT0A(valueCentFT0A);
    }

    if (fillFt0C) {
      float valueCentFT0C = 105.0f;
      if (h1dFT0C) {
        valueCentFT0C = h1dFT0C->GetBinContent(h1dFT0C->FindBin(nFT0C));
      }
      centFT0C(valueCentFT0C);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<McCentrality>(cfgc)}; }
