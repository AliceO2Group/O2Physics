// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoDreamEventHisto.h
/// \brief FemtoDreamEventHisto - Histogram class for tracks, V0s and cascades
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#ifndef PWGCF_FEMTODREAM_FEMTODREAMEVENTHISTO_H_
#define PWGCF_FEMTODREAM_FEMTODREAMEVENTHISTO_H_

#include "PWGCF/DataModel/FemtoDerived.h"
#include "Framework/HistogramRegistry.h"

using namespace o2::framework;
namespace o2::analysis::femtoDream
{
/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
class FemtoDreamEventHisto
{
 public:
  /// Destructor
  virtual ~FemtoDreamEventHisto() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  void init(HistogramRegistry* registry)
  {
    mHistogramRegistry = registry;
    mHistogramRegistry->add("Event/hZvtx", "; vtx_{z} (cm); Entries", kTH1F, {{300, -12.5, 12.5}});
    mHistogramRegistry->add("Event/hMultV0M", "; vMultV0M; Entries", kTH1F, {{16384, 0, 32768}});
    mHistogramRegistry->add("Event/hMultNTr", "; vMultNTr; Entries", kTH1F, {{200, 0, 200}});
    mHistogramRegistry->add("Event/hMultNTrVsZvtx", "; vMultNTr; vtx_{z} (cm)", kTH2F, {{200, 0, 200}, {300, -12.5, 12.5}});
    mHistogramRegistry->add("Event/hMultNTrVsMultV0M", "; vMultNTr; vMultV0M", kTH2F, {{200, 0, 200}, {16384, 0, 32768}});
  }

  /// Some basic QA of the event
  /// \tparam T type of the collision
  /// \param col Collision
  template <typename T>
  void fillQA(T const& col)
  {
    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST("Event/hZvtx"), col.posZ());
      mHistogramRegistry->fill(HIST("Event/hMultV0M"), col.multV0M());
      mHistogramRegistry->fill(HIST("Event/hMultNTr"), col.multNtr());
      mHistogramRegistry->fill(HIST("Event/hMultNTrVsZvtx"), col.multNtr(), col.posZ());
      mHistogramRegistry->fill(HIST("Event/hMultNTrVsMultV0M"), col.multNtr(), col.multV0M());
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry; ///< For QA output
};
} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_FEMTODREAMEVENTHISTO_H_
