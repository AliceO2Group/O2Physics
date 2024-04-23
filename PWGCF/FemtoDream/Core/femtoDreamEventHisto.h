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

#ifndef PWGCF_FEMTODREAM_CORE_FEMTODREAMEVENTHISTO_H_
#define PWGCF_FEMTODREAM_CORE_FEMTODREAMEVENTHISTO_H_

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
  void init(HistogramRegistry* registry, bool isMC)
  {
    mHistogramRegistry = registry;
    mHistogramRegistry->add("Event/hZvtx", "; vtx_{z} (cm); Entries", kTH1F, {{300, -12.5, 12.5}});
    mHistogramRegistry->add("Event/hMultPercentile", "; Multiplicity Percentile (FT0M); Entries", kTH1F, {{110, 0, 110}});
    mHistogramRegistry->add("Event/hMultNTr", "; Multiplicity (MultNtr); Entries", kTH1F, {{200, 0, 200}});
    mHistogramRegistry->add("Event/hMultNTrVsZvtx", "; Multiplicity (MultNtr); vtx_{z} (cm)", kTH2F, {{200, 0, 200}, {300, -12.5, 12.5}});
    mHistogramRegistry->add("Event/hMultNTrVsMultPercentile", "; Multiplicity (MultNtr); Multiplicity Percentile (FT0M)", kTH2F, {{200, 0, 200}, {110, 0, 110}});
    mHistogramRegistry->add("Event/hMultPercentileVsZvtx", "; Multiplicity Percentile (FT0M); vtx_{z} (cm)", kTH2F, {{110, 0, 110}, {300, -12.5, 12.5}});

    if (isMC) {
      mHistogramRegistry->add("Event_MC/hGenMult08VsMultPercentile", "; generated MC multiplicity (#eta<0.8); Multiplicity Percentile (FT0M)", kTH2F, {{200, 0, 200}, {110, 0, 110}});
    }
  }

  /// Some basic QA of the event
  /// \tparam T type of the collision
  /// \param col Collision
  template <bool isMC, typename T>
  void fillQA(T const& col)
  {
    if (mHistogramRegistry) {
      mHistogramRegistry->fill(HIST("Event/hZvtx"), col.posZ());
      mHistogramRegistry->fill(HIST("Event/hMultPercentile"), col.multV0M());
      mHistogramRegistry->fill(HIST("Event/hMultNTr"), col.multNtr());
      mHistogramRegistry->fill(HIST("Event/hMultNTrVsZvtx"), col.multNtr(), col.posZ());
      mHistogramRegistry->fill(HIST("Event/hMultNTrVsMultPercentile"), col.multNtr(), col.multV0M());
      mHistogramRegistry->fill(HIST("Event/hMultPercentileVsZvtx"), col.multV0M(), col.posZ());

      if constexpr (isMC) {
        if (col.has_fdMCCollision()) {
          mHistogramRegistry->fill(HIST("Event_MC/hGenMult08VsMultPercentile"), col.fdMCCollision().multMCgenPartEta08(), col.multV0M());
        }
      }
    }
  }

 private:
  HistogramRegistry* mHistogramRegistry; ///< For QA output
};
} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_CORE_FEMTODREAMEVENTHISTO_H_
