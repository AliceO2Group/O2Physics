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
#ifndef SELECTIONFILTERANDANALYSIS_H
#define SELECTIONFILTERANDANALYSIS_H

#include <Rtypes.h>
#include <TString.h>
#include <TObject.h>
#include <TNamed.h>
#include <TList.h>

#include "SkimmingConfigurableCuts.h"

namespace o2
{
namespace analysis
{
namespace PWGCF
{

/// \brief Base class for filter and selection once filetered
class SelectionFilterAndAnalysis : public TNamed
{
 public:
  /// \enum selmodes
  /// \brief the modes in which the selection can operate
  enum selmodes {
    kFilter = 0, ///< filter mode to produce skimmed data
    kAnalysis    ///< analysis mode over already skimmed data
  };

  SelectionFilterAndAnalysis();
  SelectionFilterAndAnalysis(const TString& name, selmodes mode);

 private:
  virtual int CalculateMaskLength() = 0;

 protected:
  selmodes mMode = kFilter;      /// the operating mode of the selection instance
  int mMaskLength = 0;           /// the length of the mask needed to filter the selection cuts
  ULong64_t mSelectedMask = 0UL; /// the selection mask for the current passed collision
  ULong64_t mArmedMask = 0UL;    /// the armed mask identifying the significative selection cuts

  ClassDefNV(SelectionFilterAndAnalysis, 1)
};

} // namespace PWGCF
} // namespace analysis
} // namespace o2

#endif // SELECTIONFILTERANDANALYSIS_H
