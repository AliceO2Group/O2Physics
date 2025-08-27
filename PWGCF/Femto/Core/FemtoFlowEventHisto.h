// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoFlowEventHisto.h
/// \brief FemtoFlowEventHisto - Histogram class for tracks, V0s and cascades
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de
/// \author Zuzanna Chochulska, WUT Warsaw & CTU Prague, zchochul@cern.ch

#ifndef PWGCF_FEMTOUFLOW_CORE_FEMTOUFLOWEVENTHISTO_H_
#define PWGCF_FEMTOUFLOW_CORE_FEMTOUFLOWEVENTHISTO_H_

#include "PWGCF/Femto/DataModel/FemtoDerived.h"

#include "Framework/HistogramRegistry.h"

namespace o2::analysis::femto_flow
{

/// \class FemtoFlowEventHisto
/// \brief Class for histogramming event properties
class FemtoFlowEventHisto
{
 public:
  /// Destructor
  virtual ~FemtoFlowEventHisto() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  void init(o2::framework::HistogramRegistry* registry)
  {
    using namespace o2::framework;

    mHistogramRegistry = registry;
    mHistogramRegistry->add("Event/zvtxhist", "; vtx_{z} (cm); Entries", kTH1F, {{250, -12.5, 12.5}});
    mHistogramRegistry->add("Event/MultV0M", "; vMultV0M; Entries", kTH1F, {{2000, 0, 20000}});
    mHistogramRegistry->add("Event/MultNTr", "; vMultNTr; Entries", kTH1F, {{20, 0, 200}});
    mHistogramRegistry->add("Event/MultNTrVSMultV0M", "; vMultNTr; MultV0M", kTH2F, {{200, 0, 4000}, {2000, 0, 20000}});
    mHistogramRegistry->add("Event/zvtxhist_MultNTr", "; zvtxhist; MultNTr", kTH2F, {{250, -12.5, 12.5}, {20, 0, 200}});
  }

  /// Some basic QA of the event
  /// \tparam T type of the collision
  /// \param col Collision
  template <typename T>
  void fillQA(T const& col)
  {
    using namespace o2::framework;

    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST("Event/zvtxhist"), col.posZ());
      mHistogramRegistry->fill(HIST("Event/MultV0M"), col.multV0M());
      mHistogramRegistry->fill(HIST("Event/MultNTr"), col.multNtr());
      mHistogramRegistry->fill(HIST("Event/MultNTrVSMultV0M"), col.multNtr(), col.multV0M());
      mHistogramRegistry->fill(HIST("Event/zvtxhist_MultNTr"), col.posZ(), col.multNtr());
    }
  }

 private:
  o2::framework::HistogramRegistry* mHistogramRegistry; ///< For QA output
};
} // namespace o2::analysis::femto_flow

#endif // PWGCF_FEMTOUFLOW_CORE_FEMTOUFLOWEVENTHISTO_H_
