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
#ifndef SKIMMING_CONFIGURABLE_CUTS_CLASSES_H
#define SKIMMING_CONFIGURABLE_CUTS_CLASSES_H

#include <Rtypes.h>
#include <TString.h>
#include <TObject.h>
#include <TNamed.h>
#include <TMath.h>
#include <TList.h>
#include <set>
#include <vector>
#include <regex>
#include <TObjArray.h>

#include "Framework/DataTypes.h"

namespace o2
{
namespace analysis
{
namespace PWGCF
{
/// \class CutBrick
/// \brief Virtual class which implements the base component of the selection cuts
///
template <typename TValueToFilter>
class CutBrick : public TNamed
{
 public:
  CutBrick();
  CutBrick(const char*, const char*);
  CutBrick(const CutBrick&) = delete;
  virtual ~CutBrick() override = default;
  CutBrick& operator=(const CutBrick&) = delete;

 public:
  /// Returns whether the cut brick is active alowing the selection
  /// \return true if the brick is active
  bool IsActive() { return mState == kACTIVE; }
  /// Pure virual function
  /// Returns whether the cut brick is incorporated in the selection chain
  virtual std::vector<bool> IsArmed() = 0;
  /// Pure virtual function. Filters the passed value
  /// The brick or brick components will change to active if the passed value
  /// fits within the brick or brick components scope
  /// \returns a vector of booleans with true on the component for which the value activated the component brick
  virtual std::vector<bool> Filter(const TValueToFilter&) = 0;
  /// Pure virtual function. Return the length needed to code the brick status
  /// The length is in brick units. The actual length is implementation dependent
  /// \returns Brick length in units of bricks
  virtual int Length() = 0;

  static CutBrick<TValueToFilter>* constructBrick(const char* name, const char* regex, const std::set<std::string>& allowed);
  static const char* mgImplementedbricks[];

  /// Set the status of the cut significative (or not) for the selection chain
  void Arm(bool doit = true) { doit ? mMode = kSELECTED : mMode = kUNSELECTED; }

 protected:
  /// \enum BrickStatus
  /// \brief The status of the brick
  enum BrickStatus {
    kPASSIVE, ///< if the passed value does not comply whith the brick condition
    kACTIVE   ///< if the passed value comply whith the brick condition
  };
  /// \enum BrickMode
  /// \brief The mode of operation of the brick
  enum BrickMode {
    kUNSELECTED, ///< if the status of the brick is not significative
    kSELECTED    ///< if the status of the brick is significative
  };

  BrickStatus mState = kPASSIVE;
  BrickMode mMode = kUNSELECTED;

  ClassDef(CutBrick, 1);
};

/// \class CutBrickLimit
/// \brief Class which implements a limiting cut brick.
/// The brick will be active if the filtered value is below the limit
template <typename TValueToFilter>
class CutBrickLimit : public CutBrick<TValueToFilter>
{
 public:
  CutBrickLimit();
  CutBrickLimit(const char*, const TValueToFilter&);
  CutBrickLimit(const TString&);
  virtual ~CutBrickLimit() override = default;
  CutBrickLimit(const CutBrickLimit&) = delete;
  CutBrickLimit& operator=(const CutBrickLimit&) = delete;

  virtual std::vector<bool> IsArmed();
  virtual std::vector<bool> Filter(const TValueToFilter&);
  virtual int Length() { return 1; }

 private:
  void ConstructCutFromString(const TString&);

  TValueToFilter mLimit; ///< the limiting upper value
  ClassDef(CutBrickLimit, 1);
};

/// \class CutBrickThreshold
/// \brief Class which implements a threshold cut brick.
/// The brick will be active if the filtered value is above the threshold
template <typename TValueToFilter>
class CutBrickThreshold : public CutBrick<TValueToFilter>
{
 public:
  CutBrickThreshold();
  CutBrickThreshold(const char*, const TValueToFilter&);
  CutBrickThreshold(const TString&);
  virtual ~CutBrickThreshold() override = default;
  CutBrickThreshold(const CutBrickThreshold&) = delete;
  CutBrickThreshold& operator=(const CutBrickThreshold&) = delete;

  virtual std::vector<bool> IsArmed();
  virtual std::vector<bool> Filter(const TValueToFilter&);
  virtual int Length() { return 1; }

 private:
  void ConstructCutFromString(const TString&);

  TValueToFilter mThreshold; ///< the threshold value
  ClassDef(CutBrickThreshold, 1);
};

/// \class CutBrickRange
/// \brief Class which implements a range cut brick.
/// The brick will be active if the filtered value is within the range
template <typename TValueToFilter>
class CutBrickRange : public CutBrick<TValueToFilter>
{
 public:
  CutBrickRange();
  CutBrickRange(const char*, const TValueToFilter&, const TValueToFilter&);
  CutBrickRange(const TString&);
  virtual ~CutBrickRange() override = default;
  CutBrickRange(const CutBrickRange&) = delete;
  CutBrickRange& operator=(const CutBrickRange&) = delete;

  virtual std::vector<bool> IsArmed();
  virtual std::vector<bool> Filter(const TValueToFilter&);
  virtual int Length() { return 1; }

 private:
  void ConstructCutFromString(const TString&);

  TValueToFilter mLow;  ///< the lower value of the range
  TValueToFilter mHigh; ///< the upper value of the range
  ClassDef(CutBrickRange, 1);
};

/// \class CutBrickExtToRange
/// \brief Class which implements an external to range cut brick.
/// The brick will be active if the filtered value is outside the range
template <typename TValueToFilter>
class CutBrickExtToRange : public CutBrick<TValueToFilter>
{
 public:
  CutBrickExtToRange();
  CutBrickExtToRange(const char*, const TValueToFilter&, const TValueToFilter&);
  CutBrickExtToRange(const TString&);
  virtual ~CutBrickExtToRange() override = default;
  CutBrickExtToRange(const CutBrickExtToRange&) = delete;
  CutBrickExtToRange& operator=(const CutBrickExtToRange&) = delete;

  virtual std::vector<bool> IsArmed();
  virtual std::vector<bool> Filter(const TValueToFilter&);
  virtual int Length() { return 1; }

 private:
  void ConstructCutFromString(const TString&);

  TValueToFilter mLow;  ///< the lower value of the range
  TValueToFilter mHigh; ///< the upper value of the range
  ClassDef(CutBrickExtToRange, 1);
};

/// \class CutBrickSelectorMultipleRanges
/// \brief Class which implements a string as selector an multiple ranges
/// The brick will be active if the filtered value is inside any of the multiple ranges
/// Otherwise it will be passive
/// Each range might be active in its own if the filtered value is within its scope
template <typename TValueToFilter>
class CutBrickSelectorMultipleRanges : public CutBrick<TValueToFilter>
{
 public:
  CutBrickSelectorMultipleRanges();
  CutBrickSelectorMultipleRanges(const char*, const std::vector<TValueToFilter>&);
  CutBrickSelectorMultipleRanges(const TString&);
  virtual ~CutBrickSelectorMultipleRanges() override = default;
  CutBrickSelectorMultipleRanges(const CutBrickSelectorMultipleRanges&) = delete;
  CutBrickSelectorMultipleRanges& operator=(const CutBrickSelectorMultipleRanges&) = delete;

  virtual std::vector<bool> IsArmed() override;
  virtual std::vector<bool> Filter(const TValueToFilter&) override;
  /// Return the length needed to code the brick status
  /// The length is in brick units. The actual length is implementation dependent
  /// \returns Brick length in units of bricks
  virtual int Length() override { return mActive.size(); }

 private:
  void ConstructCutFromString(const TString&);

  std::vector<TValueToFilter> mEdges; ///< the value of the ranges edges (len = nranges + 1)
  std::vector<bool> mActive;          ///< if the associated range is active with the passed value to filter (len = nranges)
  ClassDef(CutBrickSelectorMultipleRanges, 1);
};

/// \class CutWithVariations
/// \brief Class which implements a cut with its default configuration
/// and its potential variations to be used for systematic tests
template <typename TValueToFilter>
class CutWithVariations : public CutBrick<TValueToFilter>
{
 public:
  /// Default constructor
  CutWithVariations();
  CutWithVariations(const char*, const char*, bool);
  CutWithVariations(const TString&);
  virtual ~CutWithVariations() override = default;
  CutWithVariations(const CutWithVariations&) = delete;
  CutWithVariations& operator=(const CutWithVariations&) = delete;

  bool AddDefaultBrick(CutBrick<TValueToFilter>* brick);
  bool AddVariationBrick(CutBrick<TValueToFilter>* brick);
  TList& getDefaultBricks() { return mDefaultBricks; }
  TList& getVariantBricks() { return mVariationBricks; }
  virtual std::vector<bool> IsArmed();
  virtual std::vector<bool> Filter(const TValueToFilter&);
  virtual int Length();

 private:
  void ConstructCutFromString(const TString&);

  bool mAllowSeveralDefaults; ///< true if allows to store several cut default values
  TList mDefaultBricks;       ///< the list with the cut default values bricks
  TList mVariationBricks;     ///< the list with the cut variation values bricks
  ClassDef(CutWithVariations, 1);
};

/// \class SpecialCutBrick
/// \brief Virtual class which implements the base component of the special selection cuts
/// Special selection cuts are needed because the tables access seems cannot be
/// implemented in templated classes in the class dictionary. If this were found to be
/// feasible things will be easier deriving the special cut bricks directly from CutBrick
///
class SpecialCutBrick : public TNamed
{
 public:
  SpecialCutBrick();
  SpecialCutBrick(const char*, const char*);
  SpecialCutBrick(const SpecialCutBrick&) = delete;
  virtual ~SpecialCutBrick() override = default;
  SpecialCutBrick& operator=(const SpecialCutBrick&) = delete;

 public:
  /// Returns whether the cut brick is active alowing the selection
  /// \return true if the brick is active
  bool IsActive() { return mState == kACTIVE; }
  /// Pure virtual function
  /// Returns whether the cut brick is incorporated in the selection chain
  virtual std::vector<bool> IsArmed() = 0;
  /// Pure virtual function. Return the length needed to code the brick status
  /// The length is in brick units. The actual length is implementation dependent
  /// \returns Brick length in units of bricks
  virtual int Length() = 0;

 protected:
  /// Set the status of the cut significative (or not) for the selection chain
  void Arm(bool doit = true) { doit ? mMode = kSELECTED : mMode = kUNSELECTED; }

 protected:
  /// \enum BrickStatus
  /// \brief The status of the brick
  enum BrickStatus {
    kPASSIVE, ///< if the passed value does not comply whith the brick condition
    kACTIVE   ///< if the passed value comply whith the brick condition
  };
  /// \enum BrickMode
  /// \brief The mode of operation of the brick
  enum BrickMode {
    kUNSELECTED, ///< if the status of the brick is not significative
    kSELECTED    ///< if the status of the brick is significative
  };

  BrickStatus mState = kPASSIVE;
  BrickMode mMode = kUNSELECTED;

  ClassDef(SpecialCutBrick, 1);
};

class TrackSelectionBrick : public SpecialCutBrick
{
 public:
  TrackSelectionBrick() = default;
  TrackSelectionBrick(const TString&);
  virtual ~TrackSelectionBrick() override = default;

  enum class TrackCuts : int {
    kTrackType = 0,
    kTPCNCls,
    kTPCCrossedRows,
    kTPCCrossedRowsOverNCls,
    kTPCChi2NDF,
    kTPCRefit,
    kITSNCls,
    kITSChi2NDF,
    kITSRefit,
    kITSHits,
    kGoldenChi2,
    kDCAxy,
    kDCAz,
    kNCuts
  };

  static const std::string mCutNames[static_cast<int>(TrackCuts::kNCuts)];

  virtual std::vector<bool> IsArmed() override;
  template <typename TrackToFilter>
  bool Filter(TrackToFilter const& track)
  {
    if ((track.trackType() == mTrackType) &&
        (not mCheckNClustersTPC or (track.tpcNClsFound() >= mMinNClustersTPC)) &&
        (not mCheckNCrossedRowsTPC or (track.tpcNClsCrossedRows() >= mMinNCrossedRowsTPC)) &&
        (not mCheckMinNCrossedRowsOverFindableClustersTPC or (track.tpcCrossedRowsOverFindableCls() >= mMinNCrossedRowsOverFindableClustersTPC)) &&
        (not mCheckNClustersITS or (track.itsNCls() >= mMinNClustersITS)) &&
        (not mCheckMaxChi2PerClusterITS or (track.itsChi2NCl() <= mMaxChi2PerClusterITS)) &&
        (not mCheckMaxChi2PerClusterTPC or (track.tpcChi2NCl() <= mMaxChi2PerClusterTPC)) &&
        ((mRequireITSRefit) ? (track.flags() & o2::aod::track::ITSrefit) : true) &&
        ((mRequireTPCRefit) ? (track.flags() & o2::aod::track::TPCrefit) : true) &&
        ((mRequireGoldenChi2) ? (track.flags() & o2::aod::track::GoldenChi2) : true) &&
        FulfillsITSHitRequirements(track.itsClusterMap()) &&
        (mCheckMaxDcaXY && (abs(track.dcaXY()) <= ((mMaxDcaXYPtDep) ? mMaxDcaXYPtDep(track.pt()) : mMaxDcaXY))) &&
        (mCheckMaxDcaZ && (abs(track.dcaZ()) <= mMaxDcaZ))) {
      this->mState = this->kACTIVE;
      return true;
    } else {
      this->mState = this->kPASSIVE;
      return false;
    }
  }

  int Length() override { return 1; }

  void SetTrackType(o2::aod::track::TrackTypeEnum trackType) { mTrackType = trackType; }
  void SetRequireITSRefit(bool requireITSRefit = true) { mRequireITSRefit = requireITSRefit; }
  void SetRequireTPCRefit(bool requireTPCRefit = true) { mRequireTPCRefit = requireTPCRefit; }
  void SetRequireGoldenChi2(bool requireGoldenChi2 = true) { mRequireGoldenChi2 = requireGoldenChi2; }
  void SetMinNClustersTPC(int minNClustersTPC) { mMinNClustersTPC = minNClustersTPC; }
  void SetMinNCrossedRowsTPC(int minNCrossedRowsTPC) { mMinNCrossedRowsTPC = minNCrossedRowsTPC; }
  void SetMinNCrossedRowsOverFindableClustersTPC(float minNCrossedRowsOverFindableClustersTPC) { mMinNCrossedRowsOverFindableClustersTPC = minNCrossedRowsOverFindableClustersTPC; }
  void SetMinNClustersITS(int minNClustersITS) { mMinNClustersITS = minNClustersITS; }
  void SetMaxChi2PerClusterTPC(float maxChi2PerClusterTPC) { mMaxChi2PerClusterTPC = maxChi2PerClusterTPC; }
  void SetMaxChi2PerClusterITS(float maxChi2PerClusterITS) { mMaxChi2PerClusterITS = maxChi2PerClusterITS; }
  void SetMaxDcaXY(float maxDcaXY) { mMaxDcaXY = maxDcaXY; }
  void SetMaxDcaZ(float maxDcaZ) { mMaxDcaZ = maxDcaZ; }

  void SetMaxDcaXYPtDep(std::function<float(float)> ptDepCut) { mMaxDcaXYPtDep = ptDepCut; }
  void SetRequireHitsInITSLayers(int8_t minNRequiredHits, std::set<uint8_t> requiredLayers) { mRequiredITSHits.push_back(std::make_pair(minNRequiredHits, requiredLayers)); }
  void SetRequireNoHitsInITSLayers(std::set<uint8_t> excludedLayers) { mRequiredITSHits.push_back(std::make_pair(-1, excludedLayers)); }
  void ResetITSRequirements() { mRequiredITSHits.clear(); }

  void DisableNClustersTPCCheck(bool disable = true) { mCheckNClustersTPC = not disable; }
  void DisableNCrossedRowsTPCCheck(bool disable = true) { mCheckNCrossedRowsTPC = not disable; }
  void DisableNClustersITSCheck(bool disable = true) { mCheckNClustersITS = not disable; }
  void DisableMaxChi2PerClusterTPCCheck(bool disable = true) { mCheckMaxChi2PerClusterTPC = not disable; }
  void DisableMaxChi2PerClusterITSCheck(bool disable = true) { mCheckMaxChi2PerClusterITS = not disable; }
  void DisableMinNCrossedRowsOverFindableClustersTPCCheck(bool disable = true) { mCheckMinNCrossedRowsOverFindableClustersTPC = not disable; }
  void DisableMaxDcaXYCheck(bool disable = true) { mCheckMaxDcaXY = not disable; }
  void DisableMaxDcaZCheck(bool disable = true) { mCheckMaxDcaZ = not disable; }

 private:
  void constructFB1LHC2010();
  void constructFB1LHC2011();
  void constructFB32LHC2010();
  void constructFB32LHC2011();
  void constructFB64LHC2010();
  void constructFB64LHC2011();

  bool FulfillsITSHitRequirements(uint8_t itsClusterMap);

  o2::aod::track::TrackTypeEnum mTrackType{o2::aod::track::TrackTypeEnum::Track};

  // track quality cuts
  int mMinNClustersTPC{0};                            // min number of TPC clusters
  int mMinNCrossedRowsTPC{0};                         // min number of crossed rows in TPC
  int mMinNClustersITS{0};                            // min number of ITS clusters
  float mMaxChi2PerClusterTPC{1e10f};                 // max tpc fit chi2 per TPC cluster
  float mMaxChi2PerClusterITS{1e10f};                 // max its fit chi2 per ITS cluster
  float mMinNCrossedRowsOverFindableClustersTPC{0.f}; // min ratio crossed rows / findable clusters

  float mMaxDcaXY{1e10f};                       // max dca in xy plane
  float mMaxDcaZ{1e10f};                        // max dca in z direction
  std::function<float(float)> mMaxDcaXYPtDep{}; // max dca in xy plane as function of pT

  bool mRequireITSRefit{false};   // require refit in ITS
  bool mRequireTPCRefit{false};   // require refit in TPC
  bool mRequireGoldenChi2{false}; // require golden chi2 cut (Run 2 only)

  bool mCheckNClustersTPC{true};                           // check the number of TPC clusters
  bool mCheckNCrossedRowsTPC{true};                        // check the number of crossed rows int TPC
  bool mCheckNClustersITS{true};                           // check the number of ITS clusters
  bool mCheckMaxChi2PerClusterTPC{true};                   // check max tpc fit chi2 per TPC cluster
  bool mCheckMaxChi2PerClusterITS{true};                   // check max its fit chi2 per ITS cluster
  bool mCheckMinNCrossedRowsOverFindableClustersTPC{true}; // check min ratio crossed rows / findable clusters
  bool mCheckMaxDcaXY{true};                               // check max dca in xy plane
  bool mCheckMaxDcaZ{true};                                // check max dca in z direction

  // vector of ITS requirements (minNRequiredHits in specific requiredLayers)
  std::vector<std::pair<int8_t, std::set<uint8_t>>> mRequiredITSHits{};

  ClassDef(TrackSelectionBrick, 1);
};

} // namespace PWGCF
} // namespace analysis
} // namespace o2
#endif // SKIMMING_CONFIGURABLE_CUTS_CLASSES_H
