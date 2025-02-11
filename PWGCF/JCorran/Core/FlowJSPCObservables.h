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

// \brief   Helper class for the SPC-related analyses.
// \author  Maxim Virta (maxim.virta@cern.fi)

#ifndef PWGCF_JCORRAN_CORE_FLOWJSPCOBSERVABLES_H_
#define PWGCF_JCORRAN_CORE_FLOWJSPCOBSERVABLES_H_

// O2 headers. //
#include "Framework/HistogramRegistry.h"

const int maxNrComb = 12;
class FlowJSPCObservables
{
 public:
  FlowJSPCObservables() = default;

  int harmonicArray[maxNrComb][8] = {{0}};

  void setSPCObservables(int index)
  {
    // int *harmonicArray = (int*)malloc(sizeof(int)*maxNrComb*8);

    // Switch to set up correct symmetry plane combinations
    switch (index) {
      case 0: {
        LOGF(info, "Computing three harmonic SPC");
        int harmonicArray01[maxNrComb][8] = {
          {3, 6, -3, -3, 0, 0, 0, 0},
          {3, 4, -2, -2, 0, 0, 0, 0},
          {3, 8, -4, -4, 0, 0, 0, 0},
          {3, 2, 4, -6, 0, 0, 0, 0},
          {3, 2, 3, -5, 0, 0, 0, 0},
          {3, 3, 4, -7, 0, 0, 0, 0}, // These are three harmonic SPC!!
          {3, 2, 5, -7, 0, 0, 0, 0}, // These are three harmonic SPC!!
          {3, 3, 5, -8, 0, 0, 0, 0}, // These are three harmonic SPC!!
          {0, 6, -2, -2, -2, 0, 0, 0},
          {0, 2, -3, -4, 5, 0, 0, 0},
          {0, 2, -3, -3, 4, 0, 0, 0},
          {0, 3, 3, -2, -2, -2, 0, 0}};

        memcpy(harmonicArray, harmonicArray01, sizeof(int) * maxNrComb * 8);
      } break;
      case 1: {
        LOGF(info, "Computing four harmonic SPC");
        int harmonicArray02[maxNrComb][8] = {
          {4, 6, -2, -2, -2, 0, 0, 0},
          {4, 2, -3, -4, 5, 0, 0, 0},
          {4, 2, -3, -3, 4, 0, 0, 0},
          {4, 2, 2, 3, -7, 0, 0, 0}, // These are three harmonic SPC!!
          {4, 2, 2, 4, -8, 0, 0, 0}, // These are three harmonic SPC!!
          {4, 2, 7, -4, -5, 0, 0, 0},
          {4, 3, -4, -4, 5, 0, 0, 0},
          {0, 0, 0, 0, 0, 0, 0, 0},
          {0, 0, 0, 0, 0, 0, 0, 0},
          {0, 0, 0, 0, 0, 0, 0, 0},
          {0, 0, 0, 0, 0, 0, 0, 0},
          {0, 0, 0, 0, 0, 0, 0, 0}};
        memcpy(harmonicArray, harmonicArray02, sizeof(int) * maxNrComb * 8);
      } break;
      case 3: {
        LOGF(info, "Computing five and six harmonic SPC");
        int harmonicArray03[maxNrComb][8] = {
          {5, 3, 3, -2, -2, -2, 0, 0},
          {5, 2, 2, -3, 4, -5, 0, 0},
          {5, 2, 3, 3, -4, -4, 0, 0},
          {5, 3, 3, 3, -4, -5, 0, 0},
          {5, 2, 3, -4, 5, -6, 0, 0},
          {5, 8, -2, -2, -2, -2, 0, 0},
          {6, 2, 2, 2, 2, -4, -4, 0},
          {6, 2, 3, 4, 4, -6, -7, 0},
          {6, 2, 2, 2, 2, -3, -5, 0},
          {6, 2, 2, 2, 3, -4, -5, 0},
          {6, 2, 2, 3, 3, -4, -6, 0},
          {0, 0, 0, 0, 0, 0, 0, 0}};
        memcpy(harmonicArray, harmonicArray03, sizeof(int) * maxNrComb * 8);
      } break;
      default:
        LOGF(error, "ERROR: Invalid configuration index. Skipping this element.");
    }
  }

 private:
  ClassDefNV(FlowJSPCObservables, 1);
};
#endif // PWGCF_JCORRAN_CORE_FLOWJSPCOBSERVABLES_H_
