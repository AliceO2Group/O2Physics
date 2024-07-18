#include "JEPFlowAnalysis.h"

using namespace o2;
using namespace o2::framework;
using namespace std;

void JEPFlowAnalysis::FillHistograms(const Int_t cBin, Float_t det, Float_t v2, Float_t v3, Float_t v4) {
  switch (cBin) {
    case 0: {
      mHistRegistry->fill(HIST(mCentClasses[0]) + HIST("fV2EP"), v2, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[0]) + HIST("fV3EP"), v3, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[0]) + HIST("fV4EP"), v4, det, 1.);
    } break;
    case 1: {
      mHistRegistry->fill(HIST(mCentClasses[1]) + HIST("fV2EP"), v2, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[1]) + HIST("fV3EP"), v3, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[1]) + HIST("fV4EP"), v4, det, 1.);
    } break;
    case 2: {
      mHistRegistry->fill(HIST(mCentClasses[2]) + HIST("fV2EP"), v2, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[2]) + HIST("fV3EP"), v3, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[2]) + HIST("fV4EP"), v4, det, 1.);
    } break;
    case 3: {
      mHistRegistry->fill(HIST(mCentClasses[3]) + HIST("fV2EP"), v2, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[3]) + HIST("fV3EP"), v3, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[3]) + HIST("fV4EP"), v4, det, 1.);
    } break;
    case 4: {
      mHistRegistry->fill(HIST(mCentClasses[4]) + HIST("fV2EP"), v2, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[4]) + HIST("fV3EP"), v3, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[4]) + HIST("fV4EP"), v4, det, 1.);
    } break;
    case 5: {
      mHistRegistry->fill(HIST(mCentClasses[5]) + HIST("fV2EP"), v2, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[5]) + HIST("fV3EP"), v3, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[5]) + HIST("fV4EP"), v4, det, 1.);
    } break;
    case 6: {
      mHistRegistry->fill(HIST(mCentClasses[6]) + HIST("fV2EP"), v2, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[6]) + HIST("fV3EP"), v3, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[6]) + HIST("fV4EP"), v4, det, 1.);
    } break;
    case 7: {
      mHistRegistry->fill(HIST(mCentClasses[7]) + HIST("fV2EP"), v2, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[7]) + HIST("fV3EP"), v3, det, 1.);
      mHistRegistry->fill(HIST(mCentClasses[7]) + HIST("fV4EP"), v4, det, 1.);
    } break;
    default:
      return;
  }
}