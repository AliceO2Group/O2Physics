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
//
/// \author Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
/// \brief Implementation of AnalysisCompositeCut
/// \file AnalysisCompositeCut.h

#ifndef PWGDQ_CORE_ANALYSISCOMPOSITECUT_H_
#define PWGDQ_CORE_ANALYSISCOMPOSITECUT_H_

#include "PWGDQ/Core/AnalysisCut.h"

#include <RtypesCore.h>

#include <vector>

//_________________________________________________________________________
class AnalysisCompositeCut : public AnalysisCut
{
 public:
  explicit AnalysisCompositeCut(bool useAND = kTRUE);
  AnalysisCompositeCut(const char* name, const char* title, bool useAND = kTRUE);
  AnalysisCompositeCut(const AnalysisCompositeCut& c);
  AnalysisCompositeCut& operator=(const AnalysisCompositeCut& c);
  ~AnalysisCompositeCut() override;

  void addCut(AnalysisCut* cut)
  {
    if (auto* composite = dynamic_cast<AnalysisCompositeCut*>(cut)) {
      fCompositeCutList.push_back(*composite);
    } else {
      fCutList.push_back(*cut);
    }
  };

  bool getUseAND() const { return fOptionUseAND; }
  int getNCuts() const { return fCutList.size() + fCompositeCutList.size(); }

  bool isSelected(float* values) override;

 protected:
  bool fOptionUseAND;                                  // true (default): apply AND on all cuts; false: use OR
  std::vector<AnalysisCut> fCutList;                   // list of cuts
  std::vector<AnalysisCompositeCut> fCompositeCutList; // list of composite cuts
};

#endif // PWGDQ_CORE_ANALYSISCOMPOSITECUT_H_
