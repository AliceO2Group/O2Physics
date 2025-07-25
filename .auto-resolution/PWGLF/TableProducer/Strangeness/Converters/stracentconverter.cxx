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
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

using namespace o2;
using namespace o2::framework;

// Converts Stra Cents from 000 to 001
struct stracentconverter {
  Produces<aod::StraCents_001> straCents_001;

  void process(aod::StraCents_000 const& straCents_000)
  {
    for (auto& values : straCents_000) {
      straCents_001(values.centFT0M(),
                    values.centFT0A(),
                    values.centFT0C(),
                    values.centFV0A(),
                    -999., /*dummy FT0Cvariant1 value*/
                    -999., /*dummy MFT value*/
                    -999. /*dummy NGlobal value*/);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<stracentconverter>(cfgc)};
}
