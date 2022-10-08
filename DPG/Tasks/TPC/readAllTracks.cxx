#include <CCDB/BasicCCDBManager.h>
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
/// O2Physics
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/EventSelection.h"

#include "tpcSkimsTableCreator.h"
/// ROOT
#include "TRandom3.h"
#include <cmath>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;
using namespace o2::dataformats;


using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::pidTPCFullEl, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullEl, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr, aod::TrackSelection>;
using Coll = soa::Join<aod::Collisions, aod::Mults, aod::EvSels>;

namespace o2::aod
{

namespace tableProducer
{
DECLARE_SOA_COLUMN(NormMultTPC, normMultTPC, float);
DECLARE_SOA_COLUMN(NormNClustersTPC, normNClustersTPC, float);
} // namespace tpcskims

DECLARE_SOA_TABLE(ALLTRACKSTABLE, "AOD", "ALLTRACKS",
                  o2::aod::track::TPCSignal,
                  o2::aod::track::TPCInnerParam,
                  o2::aod::track::Tgl,
                  o2::aod::track::Signed1Pt,
                  o2::aod::track::Eta,
                  o2::aod::tableProducer::NormMultTPC,
                  o2::aod::tableProducer::NormNClustersTPC);
};


struct tableWriter {

    Produces<o2::aod::ALLTRACKSTABLE> allTracksTable;

    Configurable<bool> isRun2{"isRun2", false, "Normalization for the number of clusters"};
    Configurable<bool> useEventSel{"useEventSel", true, "Whether or not to use the event selection"};
    Configurable<bool> useTrackSel{"useTrackSel", true, "Whether or not to use the track selection"};

    int nClNorm = 152;

    /// Function to fill trees
    template <typename T, typename C>
    void fillTracksTable(T const& track, C const& collision)
    {

        const double ncl = track.tpcNClsFound();
        const int multTPC = collision.multTPC();

        allTracksTable(track.tpcSignal(),
            track.tpcInnerParam(),
            track.tgl(),
            track.signed1Pt(),
            track.eta(),
            multTPC / 11000.,
            std::sqrt(nClNorm / ncl));

    };

    // Event selection
    template <typename CollisionType, typename TrackType>
    bool isEvSel(const CollisionType& collision, const TrackType& tracks)
    {
        if (isRun2) {
            if (!collision.sel7()) {
                return false;
            }
            else{
                return true;
            }
        } else {
            if (!collision.sel8()) {
                return false;
            }
            else{
                return true;
            }
        }
    };

    /// Track selection
    template <typename CollisionType, typename TrackType>
    bool isTrkSel(const CollisionType& collision, const TrackType& track)
    {
        if (!track.isGlobalTrack()) { // Skipping non global tracks
            return false;
        }
        if (!track.hasITS()) { // Skipping tracks without ITS
            return false;
        }
        if (!track.hasTPC()) { // Skipping tracks without TPC
            return false;
        }
        return true;
    };

    void init(o2::framework::InitContext& initContext)
    {
        if(isRun2){
            nClNorm = 159;
        }
    };
    void process(Coll::iterator const& collision, Trks const& tracks){

        for (auto const& track : tracks) {
            bool useTrack = true;

            if(useTrackSel){
                useTrack = useTrack && isTrkSel(collision, track);
            };
            if(useEventSel){
                useTrack = useTrack && isEvSel(collision, track);
            };

            if(useTrack){ fillTracksTable(track, collision); };
        };
    };

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<tableWriter>(cfgc)};
  return workflow;
}