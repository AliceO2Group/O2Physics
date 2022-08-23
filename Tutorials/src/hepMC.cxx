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
/// \brief Accessing HepMC information
/// \author jgrosseo
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;

// Simple access to collision
struct HepMC {
  void processXSection(aod::HepMCXSection const& hepMC)
  {
    LOGF(info, "X-section: Coll id: %d, generator: %d, accepted: %lld, attempted: %lld",
         hepMC.mcCollisionId(), hepMC.generatorsID(), hepMC.accepted(), hepMC.attempted());
  }
  PROCESS_SWITCH(HepMC, processXSection, "Process x-section", true);

  void processPdfInfo(aod::HepMCPdfInfo const& hepMC)
  {
    LOGF(info, "PDFInfo: Coll id: %d, generator: %d, id1: %lld, id2: %lld",
         hepMC.mcCollisionId(), hepMC.generatorsID(), hepMC.id1(), hepMC.id2());
  }
  PROCESS_SWITCH(HepMC, processPdfInfo, "Process PDF info", true);

  void processHeavyIon(aod::HepMCHeavyIon const& hepMC)
  {
    LOGF(info, "Heavy Ion: Coll id: %d, generator: %d, ncoll_hard: %lld, ncoll: %lld",
         hepMC.mcCollisionId(), hepMC.generatorsID(), hepMC.ncollHard(), hepMC.ncoll());
  }
  PROCESS_SWITCH(HepMC, processHeavyIon, "Process Heavy Ion info", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HepMC>(cfgc)};
}
