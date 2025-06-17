#include <THashList.h>
#include <TList.h>
#include <TString.h>
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "PWGDQ/DataModel/ReducedInfoTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod;

struct eventExtendedConverter000_001 {
  Produces<aod::ReducedEventsExtended_001> eventExtended_001;

  void processExtendedConverting(aod::ReducedEventsExtended_000 const& events)
  {
    for (const auto& event : events) {
        eventExtended_001(event.globalBC(), event.alias_raw(), event.selection_raw(), event.timestamp(), event.centRun2V0M(),
            event.multTPC(), event.multFV0A(), event.multFV0C(), event.multFT0A(), event.multFT0C(),
            event.multFDDA(), event.multFDDC(), event.multZNA(), event.multZNC(), event.multTracklets(), event.multNTracksPV(),
            event.centFT0C(), -1.0f, -1.0f);
    }
  }

  void processDummy(aod::ReducedEvents& /*events*/)
  {
    // do nothing
  }

  PROCESS_SWITCH(eventExtendedConverter000_001, processExtendedConverting, "Convert Table EventsExtended_000 to Table EventsExtended_001", false);
  PROCESS_SWITCH(eventExtendedConverter000_001, processDummy, "do nothing", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<eventExtendedConverter000_001>(cfgc)};
}
