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

/// \file AliFemtoEventHisto.h
/// \brief AliFemtoEventHisto - Histogram class for tracks, V0s and cascades
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#ifndef ANALYSIS_TASKS_PWGCF_O2ALIFEMTO_INCLUDE_O2ALIFEMTO_ALIFEMTOEVENTHISTO_H_
#define ANALYSIS_TASKS_PWGCF_O2ALIFEMTO_INCLUDE_O2ALIFEMTO_ALIFEMTOEVENTHISTO_H_

//#include "AliFemtoDerived.h"
#include "PWGCF/DataModel/FemtoDerived.h"

#include "Framework/HistogramRegistry.h"

using namespace o2::framework;
namespace o2::analysis::aliFemto
{
/// \class AliFemtoEventHisto
/// \brief Class for histogramming event properties
class AliFemtoEventHisto
{
 public:
  /// Destructor
  virtual ~AliFemtoEventHisto() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  void init(HistogramRegistry* registry)
  {
    mHistogramRegistry = registry;
    mHistogramRegistry->add("Event/zvtxhist", "; vtx_{z} (cm); Entries", kTH1F, {{300, -12.5, 12.5}});
    mHistogramRegistry->add("Event/MultV0M", "; vMultV0M; Entries", kTH1F, {{600, 0, 600}});
  }

  /// Some basic QA of the event
  /// \tparam T type of the collision
  /// \param col Collision
  template <typename T>
  void fillQA(T const& col)
  {
    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST("Event/zvtxhist"), col.posZ());
      mHistogramRegistry->fill(HIST("Event/MultV0M"), col.multV0M());
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry; ///< For QA output
};
} // namespace o2::analysis::aliFemto

#endif /* ANALYSIS_TASKS_PWGCF_O2ALIFEMTO_INCLUDE_O2ALIFEMTO_ALIFEMTOEVENTHISTO_H_ */
