#include "Framework/AnalysisTask.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

//_______________________________________________________________________
TrackSelection CreateTPCOnlyTrackCuts(float max_eta, int min_tpc_ncr, float max_chi2_tpc)
{
  TrackSelection selectedTracks;
  selectedTracks.SetPtRange(0.01f, 1e10f);
  selectedTracks.SetEtaRange(-max_eta, max_eta);
  selectedTracks.SetMinNCrossedRowsTPC(min_tpc_ncr);
  selectedTracks.SetMinNCrossedRowsOverFindableClustersTPC(0.6f);
  selectedTracks.SetMaxChi2PerClusterTPC(max_chi2_tpc);
  selectedTracks.SetMaxDcaXY(2.4f);
  // selectedTracks.SetMaxDcaZ(3.2f);
  return selectedTracks;
}
//_______________________________________________________________________
bool checkAP(float alpha, float qt)
{
  const float alpha_max = 0.95;
  const float qt_max = 0.05;
  float ellipse = pow(alpha / alpha_max, 2) + pow(qt / qt_max, 2);
  if (ellipse < 1.0) {
    return true;
  } else {
    return false;
  }
}
//_______________________________________________________________________
//_______________________________________________________________________
//_______________________________________________________________________
