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

/// \file FilterAndAnalysisFramework.cxx
/// \brief Implementation of the FilterAndAnalysisFramework class for
/// filtering collisions for derived data production
/// \author victor.gonzalez.sebastian@gmail.com

#include "FilterAndAnalysisFramework.h"

#include "PWGCF/TwoParticleCorrelations/Core/EventSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/PIDSelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/SelectionFilterAndAnalysis.h"
#include "PWGCF/TwoParticleCorrelations/Core/TrackSelectionFilterAndAnalysis.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <Framework/Logger.h>

#include <TList.h>
#include <TObjString.h>

#include <Rtypes.h>

#include <cstdint>
#include <ctime>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

using namespace o2;
using namespace o2::analysis::pwgcf;

ClassImp(FilterAndAnalysisFramework);

/// \brief Default constructor
FilterAndAnalysisFramework::FilterAndAnalysisFramework()

  = default;

/// \brief Named constructor for a concrete operating mode
void FilterAndAnalysisFramework::setConfiguration(EventSelectionConfigurable& evtf, TrackSelectionConfigurable& trkf, PIDSelectionConfigurable& pidf, SelectionFilterAndAnalysis::Selmodes mode)
{
  tmpTrackFilter = new pwgcf::TrackSelectionFilterAndAnalysis(trkf, mode);
  tmpEventFilter = new pwgcf::EventSelectionFilterAndAnalysis(evtf, mode);
  tmpPIDFilter = new pwgcf::PIDSelectionFilterAndAnalysis(pidf, mode);
}

void FilterAndAnalysisFramework::registerConfiguration()
{
  LOGF(info, "Registering configuration into CCDB for date %s", filterdate.data());

  /* sanity checks */
  if (tmpTrackFilter == nullptr || tmpEventFilter == nullptr || tmpPIDFilter == nullptr) {
    LOGF(fatal, "Configuration not stored yet, please use SetConfiguration(evtf, trkf, pidf, mode), or configuration already registered");
  }
  if (fTrackFilter != nullptr || fEventFilter != nullptr || fPIDFilter != nullptr) {
    LOGF(fatal, "Configuration already registered");
  }

  TList* signatures = nullptr;
  std::tm t{};
  /* we set noon for checking */
  std::istringstream ss(filterdate + "1200");
  /* let's extract the timestamp */
  ss >> std::get_time(&t, "%Y%m%d%H%M");
  if (ss.fail()) {
    LOGF(fatal, "Wrong skimming filter date format. Please use YYYYMMDD");
  }
  std::time_t timeStamp = mktime(&t);
  /* and get it in ms */
  uint64_t timestamp = 1000 * timeStamp;
  /* let's prepare the sot and eot */
  uint64_t sot = timestamp - (12 * 60 * 60 * 1000);        /* the beginning of the day */
  uint64_t eot = timestamp + (12 * 60 * 60 * 1000 - 1000); /* one second before the end of the day */

  /* get access to the CCDB */
  auto* ccdb = new o2::ccdb::CcdbApi();
  ccdb->init(ccdburl);
  o2::ccdb::BasicCCDBManager* ccdbm = &o2::ccdb::BasicCCDBManager::instance();
  ccdbm->setFatalWhenNull(false); /* we need to check the availability of objects */
  ccdbm->setURL(ccdburl);
  signatures = ccdbm->getForTimeStamp<TList>(ccdbpath, timestamp);

  auto storesignatures = [&]() {
    auto* siglst = new TList();
    siglst->SetOwner();
    siglst->SetName("Skimming signatures");
    siglst->Add(new TObjString(tmpEventFilter->getCutStringSignature().Data()));
    siglst->Add(new TObjString(tmpTrackFilter->getCutStringSignature().Data()));
    siglst->Add(new TObjString(tmpPIDFilter->getCutStringSignature().Data()));
    std::map<std::string, std::string> metadata = {{"FilterDate", filterdate}};
    ccdb->storeAsTFileAny<TList>(siglst, ccdbpath, metadata, sot, eot);
    fEventFilter = tmpEventFilter;
    fTrackFilter = tmpTrackFilter;
    fPIDFilter = tmpPIDFilter;
    tmpEventFilter = nullptr;
    tmpTrackFilter = nullptr;
    tmpPIDFilter = nullptr;
  };

  if (signatures != nullptr) {
    /* signatures already stored in the CCDB, check if they match with the current configuration */
    auto domatch = [&signatures](int index, auto filter) {
      return (strcmp(signatures->At(index)->GetName(), filter->getCutStringSignature().Data()) == 0);
    };
    if (domatch(0, tmpEventFilter) && domatch(1, tmpTrackFilter) && domatch(2, tmpPIDFilter)) {
      /* signatures match so configuration already registered and available for skimming and/or analysis */
      LOGF(info, "Filter signatures already registered in the CCDB. You can proceed with the skimming and/or analysis");
    } else {
      /* signatures don't match the version on the CCDB */
      if (forceupdate) {
        /* forced update required so store the configuration signatures in the CCDB */
        LOGF(info, "Filter configuration signatures registered into the CCDB");
        storesignatures();
      } else {
        /* no forced update required, inform the user */
        LOGF(error, "The filter configuration signatures for the requested date don't match the ones available at the CCDB");
        LOGF(error, "If the new filter configuration is registered the previous one will not be anylonger available");
        LOGF(error, "If you agree on that please invoke againg the registering with --forced=true");
        LOGF(error, "Nothing was done!!!");
      }
    }
  } else {
    /* signatures not stored in the CCDB, let's store them */
    LOGF(info, "Filter configuration signatures registered into the CCDB");
    storesignatures();
  }
}

void FilterAndAnalysisFramework::init()
{
  LOGF(info, "Initializing filter configuration from date %s", filterdate.data());

  /* sanity checks */
  if (tmpTrackFilter == nullptr || tmpEventFilter == nullptr || tmpPIDFilter == nullptr) {
    LOGF(fatal, "Configuration not stored yet, please use SetConfiguration(evtf, trkf, pidf, mode, force), or configuration already initalized");
  }
  if (fTrackFilter != nullptr || fEventFilter != nullptr || fPIDFilter != nullptr) {
    LOGF(fatal, "Configuration already initialized");
  }

  TList* signatures = nullptr;
  std::tm t{};
  /* we set noon for checking */
  std::istringstream ss(filterdate + "1200");
  /* let's extract the timestamp */
  ss >> std::get_time(&t, "%Y%m%d%H%M");
  if (ss.fail()) {
    LOGF(fatal, "Wrong skimming filter date format. Please use YYYYMMDD");
  }
  std::time_t timeStamp = mktime(&t);
  /* and get it in ms */
  uint64_t timestamp = 1000 * timeStamp;

  /* get access to the CCDB */
  o2::ccdb::BasicCCDBManager* ccdbm = &o2::ccdb::BasicCCDBManager::instance();
  ccdbm->setFatalWhenNull(false); /* we need to check the availability of objects */
  ccdbm->setURL(ccdburl);
  signatures = ccdbm->getForTimeStamp<TList>(ccdbpath, timestamp);

  if (signatures != nullptr) {
    /* signatures registered in the CCDB, check if they match with the current configuration */
    auto domatch = [&signatures](int index, auto filter) {
      LOGF(info,
           "Comparing signatures\n"
           "\t\t%s\n"
           "\tand"
           "\t\t%s",
           signatures->At(index)->GetName(), filter->getCutStringSignature().Data());
      return (strcmp(signatures->At(index)->GetName(), filter->getCutStringSignature().Data()) == 0);
    };
    if (domatch(0, tmpEventFilter) && domatch(1, tmpTrackFilter) && domatch(2, tmpPIDFilter)) {
      /* signatures match the configuration proceed with skimming and/or analysis */
      LOGF(info, "Filter signatures registered in the CCDB. Proceeding with the skimming and/or analysis");
      fEventFilter = tmpEventFilter;
      fTrackFilter = tmpTrackFilter;
      fPIDFilter = tmpPIDFilter;
      tmpEventFilter = nullptr;
      tmpTrackFilter = nullptr;
      tmpPIDFilter = nullptr;
    } else {
      /* signatures don't match the version on the CCDB */
      LOGF(fatal,
           "The stored configuration signatures don't match with the ones registered at the CCDB\n"
           "\t\tEither you are trying to use a new filter configuration which is not registerd, please register it,\n"
           "\t\tor you have changed the configuration for carrying your analysis, but it will not match the skimmed data,\n"
           "\t\tplease do not modify the configuration, just the selection \"yes\" \"not\" flags.");
    }
  } else {
    /* signatures not stored in the CCDB, abort! */
    LOGF(fatal, "Filter configuration signatures not stored in the CCDB!!! We cannot proceed!!!");
  }
}
