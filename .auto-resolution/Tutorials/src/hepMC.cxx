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

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;

// Simple access to collision
struct HepMc {
  OutputObj<TH2F> xsection{TH2F("xsection", "collision id vs. generator ID", 1000, -0.5, 999.5, 20, -0.5, 19.5)};
  OutputObj<TH2F> pdfinfo{TH2F("pdfinfo", "collision id vs. generator ID", 1000, -0.5, 999.5, 20, -0.5, 19.5)};
  OutputObj<TH2F> heavyion{TH2F("heavyion", "collision id vs. generator ID", 1000, -0.5, 999.5, 20, -0.5, 19.5)};

  void processXSection(aod::HepMCXSection const& hepMC)
  {
    LOGF(debug, "X-section: Coll id: %d, generator: %d, accepted: %lld, attempted: %lld",
         hepMC.mcCollisionId(), hepMC.generatorsID(), hepMC.accepted(), hepMC.attempted());
    xsection->Fill(hepMC.mcCollisionId(), hepMC.generatorsID());
  }
  PROCESS_SWITCH(HepMc, processXSection, "Process x-section", true);

  void processPdfInfo(aod::HepMCPdfInfo const& hepMC)
  {
    LOGF(debug, "PDFInfo: Coll id: %d, generator: %d, id1: %lld, id2: %lld",
         hepMC.mcCollisionId(), hepMC.generatorsID(), hepMC.id1(), hepMC.id2());
    pdfinfo->Fill(hepMC.mcCollisionId(), hepMC.generatorsID());
  }
  PROCESS_SWITCH(HepMc, processPdfInfo, "Process PDF info", true);

  void processHeavyIon(aod::HepMCHeavyIon const& hepMC)
  {
    LOGF(debug, "Heavy Ion: Coll id: %d, generator: %d, ncoll_hard: %lld, ncoll: %lld",
         hepMC.mcCollisionId(), hepMC.generatorsID(), hepMC.ncollHard(), hepMC.ncoll());
    heavyion->Fill(hepMC.mcCollisionId(), hepMC.generatorsID());
  }
  PROCESS_SWITCH(HepMc, processHeavyIon, "Process Heavy Ion info", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<HepMc>(cfgc)};
}
