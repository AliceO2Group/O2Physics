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

  /// \brief get the valid (armed) mask associated to the configuration
  /// \return the armed mask
  uint64_t getMask() { return mArmedMask; }
  /// \brief get the valid (armed) optional part mask associated to the configuration
  /// \return the armed optional mask
  /// A clear example of the optional part mask is the mask of the multiplicity classes
  /// where only one of the available in the whole mask will be flagged
  std::vector<uint64_t>& getOptMask() { return mOptArmedMask; }
  TString printOptionalMasks() const;
  /// \brief get the valid (armed) mandatory part mask associated to the configuration
  /// \return the armed optional mask
  /// A clear example of the mandatory part mask is the mask of the zvertex and their
  /// alternatives where only one concrete one is required to be flagged
  uint64_t getForcedMask() { return mForcedArmedMask; }
  /// \brief get the cut string signature
  const TString& getCutStringSignature() { return mCutStringSignature; }

 private:
  virtual int CalculateMaskLength() = 0;
  virtual void StoreArmedMask() = 0;

 protected:
  selmodes mMode = kFilter;                 /// the operating mode of the selection instance
  TString mCutStringSignature;              /// the signature of the cut string, i.e. without the alternative selection tags
  int mMaskLength = 0;                      /// the length of the mask needed to filter the selection cuts
  uint64_t mSelectedMask = 0UL;             /// the selection mask for the current passed collision
  uint64_t mArmedMask = 0UL;                /// the complete armed mask identifying the applicable selection cuts
  std::vector<uint64_t> mOptArmedMask = {}; /// the list of armed masks for options of the applicable selection cuts
                                            /// all of them have to have at least one option active
  uint64_t mForcedArmedMask = 0UL;          /// the mandatory armed mask of the applicable selection cuts

  ClassDefNV(SelectionFilterAndAnalysis, 1)
};

} // namespace PWGCF
} // namespace analysis
} // namespace o2

#endif // SELECTIONFILTERANDANALYSIS_H
