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

/// \file TListHandler.h
/// \brief TListHandler class to add, get and fill the histograms stored in a root TList
///        in similar way as HisgtogramRegistry works. most of code is written with help of
///        functions and methods defined in HistogramRegistry.h and HisgtogramRegistry.cxx
/// \author Rahul Verma (rahul.verma@cern.ch, rahul.verma@iitb.ac.in)

#ifndef COMMON_TOOLS_TLISTHANDLER_H_
#define COMMON_TOOLS_TLISTHANDLER_H_

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramSpec.h"

#include <TObjString.h>

#include <algorithm>
#include <functional>
#include <memory>
#include <regex>
#include <string>
#include <utility>
#include <vector>

namespace o2::framework
{

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
// TListHandler Class
// Class to handle the histograms stored in a TList in efficient way just like Histogram Registry.
// Many methods/definitions/parts of code are taken from HistogramRegisty.h and HistogramRegistry.cxx
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
// ToDo :: create templated class to keep the BitMask modifiable and declrable at compile time. // template <uint32_t BitMask>
class TListHandler
{

  struct HistName {
    // ctor for histogram names that are already hashed at compile time via HIST("myHistName")
    template <char... chars>
    constexpr HistName(const ConstStr<chars...>& hashedHistName); // NOLINT(runtime/explicit)
    char const* const str{};
    const uint32_t hash{};
    const uint32_t idx{};

   protected:
    friend class TListHandler;
    // ctor that does the hashing at runtime (for internal use only)
    constexpr HistName(char const* const name); // NOLINT(runtime/explicit)
  };

 public:
  TList* rootList;
  bool makeNestedList;
  TListHandler() : rootList(nullptr) {}                  // Default constructor
  explicit TListHandler(TList* list) : rootList(list) {} // Constructor accepting a TList*

  // operator() overload to assign a new list
  void operator()(TList* list, bool makeNestedListFlag)
  {
    rootList = list;
    makeNestedList = makeNestedListFlag;
    mName = rootList->GetName();
  }

  inline uint32_t probeDistance(uint32_t hash_idx, uint32_t slot_idx, uint32_t capacity)
  {
    return ((capacity + slot_idx - hash_idx) & (capacity - 1)); // simplification of logic  (slot_idx >= hash_idx) ? (slot_idx - hash_idx) : (capacity + slot_idx)
  }

  // functions to add histograms to the TList
  HistPtr add(const HistogramSpec& histSpec);
  HistPtr add(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2 = false);
  HistPtr add(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2 = false);
  HistPtr add(const std::string& name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2 = false);

  template <typename T>
  std::shared_ptr<T> add(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2 = false);
  template <typename T>
  std::shared_ptr<T> add(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2 = false);
  template <typename T>
  std::shared_ptr<T> add(const std::string& name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2 = false);

  void addClone(const std::string& source, const std::string& target);

  // get the underlying histogram pointer
  template <typename T>
  std::shared_ptr<T> get(const HistName& histName);

  // fill hist with values
  template <typename... Ts>
  void fill(const HistName& histName, Ts... positionAndWeight)
    requires(FillValue<Ts> && ...);

  // fill hist with content of (filtered) table columns
  template <typename... Cs, typename T>
  void fill(const HistName& histName, const T& table, const o2::framework::expressions::Filter& filter);

  // get rough estimate for size of histogram stored in TList
  double getSize(const HistName& histName, double fillFraction = 1.);

  // get rough estimate for size of all histograms stored in TList
  double getSize(double fillFraction = 1.);

  /// deletes all the histograms from the registry
  void clean();

  // print summary of the histograms stored in registry
  void print(bool showAxisDetails = false);

  void totalShifts(bool printDetails = false)
  {
    uint32_t totalLookups = 0;
    for (uint i = 0; i < kMaxTListSize; i++) {
      if (mTListKey[i] != 0) {
        if (printDetails) {
          LOG(info) << "DEBUG :: " << "i = " << i << " :: " << imask(mTListKey[i]) << " :: " << probeDistance(imask(mTListKey[i]), i, kMaxTListSize) << " ::  ";
        }
        totalLookups += probeDistance(imask(mTListKey[i]), i, kMaxTListSize);
      }
    }
    LOG(info) << "DEBUG :: " << "totalLookups = " << totalLookups << " ::  ";
  }

  static HistName makeHistName(const std::string& name);

  mutable uint32_t lookup = 0;

  std::vector<std::string> mRegisteredNames{};

 private:
  uint32_t performRobinhoodSwapAndGetIndex(const uint32_t histSpecHash, const HistPtr hist);
  void insertInNestedTList(const std::string& name, const HistPtr hist);

  // create histogram from specification and insert it into the TList
  HistPtr insert(const HistogramSpec& histSpec);

  // clone an existing histogram and insert it into the TList
  template <typename T>
  HistPtr insertClone(const HistName& histName, const std::shared_ptr<T> originalHist);

  void cloningInList(const TList* sourceList, const std::string& source, const std::string& target);

  // function to query if name is already in use
  bool contains(const HistName& histName);

  // helper function to find the histogram position in the TList
  template <typename T>
  uint32_t getHistIndex(const T& histName);

  // helper function that checks if histogram name can be used in TList
  void validateHistName(const std::string& name, const uint32_t hash);

  constexpr uint32_t imask(uint32_t i) const
  {
    return i & kTListBitmask;
  }

  // void insertInNestedTList(const std::string& name, const HistPtr hist){

  void registerName(const std::string& name);

  std::string mName{};
  bool mSortHistos{}; // Sorting to be implement to a TList in future.
  // std::vector<std::string> mRegisteredNames{};

  // To Do - make templated class with bitmask --> static constexpr uint32_t kTListBitmask = BitMask;.
  // {0x1FF = 511}; {0x3FF = 1023}; {0x7FF = 2047}; {0xFFF = 4095}
  static constexpr uint32_t kTListBitmask{0x7FF};             // o2-linter: disable=name/constexpr-constant,name/function-variable
  static constexpr uint32_t kMaxTListSize{kTListBitmask + 1}; // o2-linter: disable=name/constexpr-constant,name/function-variable
  static_assert((kMaxTListSize & (kMaxTListSize - 1)) == 0, "BitMask must be a power of 2^(n)-1");
  std::array<uint32_t, kMaxTListSize> mTListKey{};
  std::array<HistPtr, kMaxTListSize> mTListValue{};
  // HistPtr
};

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
// Implementation of TListHandler template functions.
// HistFiller template functions are already defined in the Histogramspec.h file
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

template <typename T>
concept IsTListHandler = std::same_as<T, TListHandler>;

template <typename T>
struct IsVariant : std::false_type {
};

template <typename... Args>
struct IsVariant<std::variant<Args...>> : std::true_type {
};

template <typename T>
constexpr bool IsVariantV = IsVariant<T>::value;

template <typename T>
constexpr bool AlwaysFalseV = false;

template <typename T>
const char* getTObjectName(const T& obj)
{
  if constexpr (std::is_base_of_v<TObject, std::remove_pointer_t<T>>) {
    return obj ? obj->GetName() : nullptr;
  } else if constexpr (IsVariantV<T>) {
    const char* name = nullptr;
    std::visit([&](auto&& val) {
      if (val)
        name = val->GetName();
    },
               obj);
    return name;
  } else {
    static_assert(AlwaysFalseV<T>, "Unsupported type in getTObjectName");
  }
}

TListHandler::HistName TListHandler::makeHistName(const std::string& name)
{
  return HistName{name.c_str()};
}

// function to query if name is already in use
bool TListHandler::contains(const HistName& histName)
{
  // check for all occurances of the hash
  auto iter = mTListKey.begin();
  while ((iter = std::find(iter, mTListKey.end(), histName.hash)) != mTListKey.end()) {
    LOG(fatal) << "DEBUG :: " << "Hash collision Occured " << " ::  ";
    auto& val = mTListValue[iter - mTListKey.begin()];
    const char* curName = getTObjectName(val);
    if (strcmp(curName, histName.str) == 0) {
      return true;
    }
    iter++;
  }
  return false;
}

template <char... chars>
constexpr TListHandler::HistName::HistName(const ConstStr<chars...>& hashedHistName)
  : str(hashedHistName.str),
    hash(hashedHistName.hash),
    idx(hash & kTListBitmask) // BITMASK changed
{
}

template <typename T>
std::shared_ptr<T> TListHandler::add(const std::string& name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2)
{
  return add<T>(name.c_str(), title, histType, axes, callSumw2);
}

template <typename T>
std::shared_ptr<T> TListHandler::get(const HistName& histName)
{
  if (auto histPtr = std::get_if<std::shared_ptr<T>>(&mTListValue[getHistIndex(histName)])) {
    return *histPtr;
  } else {
    throw runtime_error_f(R"(Histogram type specified in get<>(HIST("%s")) does not match the actual type of the histogram!)", histName.str);
  }
}

template <typename T>
uint32_t TListHandler::getHistIndex(const T& histName)
{
  if (O2_BUILTIN_LIKELY(histName.hash == mTListKey[histName.idx])) {
    return histName.idx;
  }
  for (auto i = 1u; i < kMaxTListSize; ++i) {
    if (histName.hash == mTListKey[imask(histName.idx + i)]) {
      return imask(histName.idx + i);
    }
  }
  throw runtime_error_f(R"(Could not find histogram "%s" in TListHandler "%s"!)", histName.str, mName.data());
}

template <typename... Ts>
void TListHandler::fill(const HistName& histName, Ts... positionAndWeight)
  requires(FillValue<Ts> && ...)
{
  std::visit([positionAndWeight...](auto&& hist) { HistFiller::fillHistAny(hist, positionAndWeight...); }, mTListValue[getHistIndex(histName)]);
}

extern template void TListHandler::fill(const HistName& histName, double);
extern template void TListHandler::fill(const HistName& histName, float);
extern template void TListHandler::fill(const HistName& histName, int);

extern template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<TH1>);
extern template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<TH2>);
extern template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<TH3>);
extern template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<TProfile>);
extern template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<TProfile2D>);
extern template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<TProfile3D>);
extern template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<THnSparse>);
extern template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<THn>);
extern template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<StepTHn>);

extern template std::shared_ptr<TH1> TListHandler::add<TH1>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
extern template std::shared_ptr<TH1> TListHandler::add<TH1>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);
extern template std::shared_ptr<TH2> TListHandler::add<TH2>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
extern template std::shared_ptr<TH2> TListHandler::add<TH2>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);
extern template std::shared_ptr<TH3> TListHandler::add<TH3>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
extern template std::shared_ptr<TH3> TListHandler::add<TH3>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);
extern template std::shared_ptr<TProfile> TListHandler::add<TProfile>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
extern template std::shared_ptr<TProfile> TListHandler::add<TProfile>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);
extern template std::shared_ptr<TProfile2D> TListHandler::add<TProfile2D>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
extern template std::shared_ptr<TProfile2D> TListHandler::add<TProfile2D>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);
extern template std::shared_ptr<TProfile3D> TListHandler::add<TProfile3D>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
extern template std::shared_ptr<TProfile3D> TListHandler::add<TProfile3D>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);
extern template std::shared_ptr<THn> TListHandler::add<THn>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
extern template std::shared_ptr<THn> TListHandler::add<THn>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);
extern template std::shared_ptr<THnSparse> TListHandler::add<THnSparse>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
extern template std::shared_ptr<THnSparse> TListHandler::add<THnSparse>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);
extern template std::shared_ptr<StepTHn> TListHandler::add<StepTHn>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
extern template std::shared_ptr<StepTHn> TListHandler::add<StepTHn>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);

template <typename... Cs, typename T>
void TListHandler::fill(const HistName& histName, const T& table, const o2::framework::expressions::Filter& filter)
{
  std::visit([&table, &filter](auto&& hist) { HistFiller::fillHistAny<Cs...>(hist, table, filter); }, mTListValue[getHistIndex(histName)]);
}

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
// Implementation Taken From HistogramRegistry.cxx File
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

template void TListHandler::fill(const HistName& histName, double);
template void TListHandler::fill(const HistName& histName, float);
template void TListHandler::fill(const HistName& histName, int);

// runtime  constructor
constexpr TListHandler::HistName::HistName(char const* const name)
  : str(name),
    hash(runtime_hash(name)),
    idx(hash & kTListBitmask)
{
}

// create histogram from specification and insert it into the TList
HistPtr TListHandler::insert(const HistogramSpec& histSpec)
{
  validateHistName(histSpec.name, histSpec.hash);
  HistPtr hist = HistFactory::createHistVariant(histSpec);
  performRobinhoodSwapAndGetIndex(histSpec.hash, hist);
  mRegisteredNames.push_back(histSpec.name);
  insertInNestedTList(histSpec.name, hist);
  return hist;
}

void TListHandler::insertInNestedTList(const std::string& name, const HistPtr hist)
{
  // created nested TList structue to easy reading of file.
  HistPtr tObj = hist;
  TList* parentList = rootList;
  if (makeNestedList) {
    std::string pathStr(name);
    std::string dirPath, histName;

    auto pos = pathStr.find_last_of('/');
    if (pos != std::string::npos) {
      dirPath = pathStr.substr(0, pos);
      histName = pathStr.substr(pos + 1);
    } else {
      histName = pathStr;
    }

    TList* subList = nullptr;

    if (!dirPath.empty()) {
      TString dirPathT(dirPath.c_str());
      TObjArray* folders = dirPathT.Tokenize('/');
      for (int i = 0; i < folders->GetEntries(); ++i) {
        TString subdir = (static_cast<TObjString*>(folders->At(i)))->GetString();
        TObject* existingObj = parentList->FindObject(subdir);
        if (existingObj && existingObj->InheritsFrom(TList::Class())) {
          subList = static_cast<TList*>(existingObj);
        } else {
          subList = new TList();
          subList->SetName(subdir);
          parentList->Add(subList);
        }
        parentList = subList;
      }
      delete folders;
    }

    TNamed* rawPtrToObj = nullptr;
    std::visit([&](const auto& sharedPtr) { rawPtrToObj = static_cast<TNamed*>(sharedPtr.get()); }, tObj);
    rawPtrToObj->SetName(histName.c_str());
  }

  TObject* rawPtrToObj = nullptr;
  std::visit([&](const auto& sharedPtr) { rawPtrToObj = sharedPtr.get(); }, tObj);
  parentList->Add(rawPtrToObj);
}

void sortTListRecursively(TList* list)
{
  if (!list)
    return;

  list->Sort(); // This uses TObject::Compare if implemented (e.g., TNamed compares by name)

  for (TObject* obj : *list) { // o2-linter: disable=const-ref-in-for-loop (not constant)
    if (obj->InheritsFrom(TList::Class())) {
      sortTListRecursively(static_cast<TList*>(obj));
    }
  }
}

uint32_t TListHandler::performRobinhoodSwapAndGetIndex(const uint32_t histSpecHash, const HistPtr hist)
{
  const uint32_t idx = imask(histSpecHash);
  uint32_t swappedHash = 0;
  HistPtr swappedPtr; // is place of std::remove_reference_t<decltype(mTListValue[0])> swappedPtr = nullptr;
  uint32_t pos = 0;
  uint32_t returnPos = kMaxTListSize;
  for (auto i = 0u; i < kMaxTListSize; ++i) {
    pos = imask(idx + i);
    TObject* rawPtr = nullptr;
    std::visit([&](const auto& sharedPtr) { rawPtr = sharedPtr.get(); }, mTListValue[pos]);
    if (!rawPtr) {
      if (swappedHash == 0) { // No Swapping Happened
        mTListKey[pos] = histSpecHash;
        mTListValue[pos] = hist;
        returnPos = pos;
      } else {
        std::swap(swappedPtr, mTListValue[pos]);
        std::swap(swappedHash, mTListKey[pos]);
      }

      lookup += i;
      break;
    } else {
      uint32_t activeHash = (swappedHash == 0) ? histSpecHash : swappedHash;
      uint32_t currntHashDist = probeDistance(imask(activeHash), pos, kMaxTListSize);
      uint32_t storedHashDist = probeDistance(imask(mTListKey[pos]), pos, kMaxTListSize);

      if (storedHashDist < currntHashDist) {
        if (swappedHash == 0) {
          swappedPtr = hist;
          swappedHash = histSpecHash;
          returnPos = pos;
        }
        std::swap(swappedPtr, mTListValue[pos]);
        std::swap(swappedHash, mTListKey[pos]);
        lookup++;
      }
    }
  }

  if (swappedHash != 0) {
    LOGF(fatal, R"(Internal array of TListHandler "%s" is full.)", mName);
    throw std::runtime_error("Insertion failed: Hashmap is full");
  }

  if (returnPos == kMaxTListSize) {
    LOGF(fatal, R"(TListHandler "%s" bad return value, histogram not assigned)", mName);
  }
  return returnPos;
}

// helper function that checks if histogram name can be used in TList
void TListHandler::validateHistName(const std::string& name, const uint32_t hash)
{
  // check that there are still slots left in the TList
  if (mRegisteredNames.size() == kMaxTListSize) {
    LOGF(fatal, R"(TListHandler "%s" is full! It can hold only %d histograms.)", mName, kMaxTListSize);
  }

  // validate that hash is unique
  auto it = std::find(mTListKey.begin(), mTListKey.end(), hash);
  if (it != mTListKey.end()) {
    auto idx = it - mTListKey.begin();
    std::string collidingName{};
    // std::visit([&](const auto& hist) { collidingName = hist->GetName(); }, mTListValue[idx]);
    auto& val = mTListValue[idx];
    collidingName = getTObjectName(val);
    LOGF(fatal, R"(Hash collision in TListHandler "%s"! Please rename histogram "%s" or "%s".)", mName, name, collidingName);
  }

  // validate that name contains only allowed characters
  if (!std::regex_match(name, std::regex("([a-zA-Z0-9])(([\\/_-])?[a-zA-Z0-9])*"))) {
    LOGF(fatal, R"(Histogram name "%s" contains invalid characters. Only letters, numbers, and (except for the beginning or end of the word) the special characters '/', '_', '-' are allowed.)", name);
  }
}

HistPtr TListHandler::add(const HistogramSpec& histSpec)
{
  return insert(histSpec);
}

HistPtr TListHandler::add(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2)
{
  return insert({name, title, histConfigSpec, callSumw2});
}

HistPtr TListHandler::add(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2)
{
  return insert({name, title, {histType, axes}, callSumw2});
}

HistPtr TListHandler::add(const std::string& name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2)
{
  return add(name.c_str(), title, histType, axes, callSumw2);
}

const TList* getRootSubList(TList* rootList, std::string dirPath)
{
  const TList* parentList = rootList;
  if (!dirPath.empty()) {
    TString dirPathT(dirPath.c_str());
    TObjArray* folders = dirPathT.Tokenize('/');
    for (int i = 0; i < folders->GetEntries(); ++i) {
      TString subdir = (static_cast<TObjString*>(folders->At(i)))->GetString();
      TObject* existingObj = parentList->FindObject(subdir);
      if (existingObj && existingObj->InheritsFrom(TList::Class())) {
        parentList = static_cast<const TList*>(existingObj);
      } else {
        delete folders;
        throw std::runtime_error(Form("In TList path %s, object %s is not a TList", dirPath.c_str(), subdir.Data()));
      }
    }
    delete folders;
  }
  return parentList;
}

void TListHandler::cloningInList(const TList* sourceList, const std::string& source, const std::string& target)
{
  TIter next(sourceList);
  TObject* obj = nullptr;

  sourceList->Print();
  while ((obj = next())) {
    std::string currentName = obj->GetName();
    std::string fullSourceName = source + currentName;
    std::string fullTargetName = target + currentName;

    if (obj->InheritsFrom(TList::Class())) {
      // If it's a nested directory, recursively process it
      const TList* subSourceList = dynamic_cast<const TList*>(obj);
      if (!subSourceList) {
        LOGF(fatal, "Failed to cast to TList");
        return;
      }

      // Append '/' to preserve nested directory structure
      std::string newSource = fullSourceName + "/";
      std::string newTarget = fullTargetName + "/";

      // Recursive call inside the function
      cloningInList(subSourceList, newSource, newTarget);
    } else {
      // It’s a histogram — clone it
      auto& sharedPtr = mTListValue[getHistIndex(makeHistName(fullSourceName))];
      HistName targetHistName = makeHistName(fullTargetName);
      std::visit([&](const auto& ptr) {
        insertClone(targetHistName, ptr);
      },
                 sharedPtr);
    }
  }
}

// store a copy of an existing histogram (or group of histograms) under a different name
void TListHandler::addClone(const std::string& source, const std::string& target)
{
  if (!makeNestedList) {
    auto doInsertClone = [&](const auto& sharedPtr) {
      if (!sharedPtr.get()) {
        return;
      }
      std::string sourceName{(static_cast<TNamed*>(sharedPtr.get()))->GetName()};
      // search for histograms starting with source_ substring
      if (sourceName.rfind(source, 0) == 0) {
        // when cloning groups of histograms source_ and target_ must end with "/"
        if (sourceName.size() != source.size() && (source.back() != '/' || target.back() != '/')) {
          return;
        }
        // when cloning a single histogram the specified target_ must not be a group name
        if (sourceName.size() == source.size() && target.back() == '/') {
          LOGF(fatal, "Cannot turn histogram into folder!");
        }
        std::string targetName{target};
        targetName += sourceName.substr(sourceName.find(source) + source.size());
        insertClone(targetName.data(), sharedPtr);
      }
    };

    for (auto& histVariant : mTListValue) { // o2-linter: disable=const-ref-in-for-loop (not constant)
      std::visit(doInsertClone, histVariant);
    }
  } else {
    const TList* sourceList = getRootSubList(rootList, source);
    cloningInList(sourceList, source, target);
  }
}

// get rough estimate for size of histogram stored in registry
double TListHandler::getSize(const HistName& histName, double fillFraction)
{
  double size{};
  std::visit([&fillFraction, &size](auto&& hist) { size = HistFiller::getSize(hist, fillFraction); }, mTListValue[getHistIndex(histName)]);
  return size;
}

// get rough estimate for size of all histograms stored in registry
double TListHandler::getSize(double fillFraction)
{
  double size{};
  for (auto j = 0u; j < kMaxTListSize; ++j) {
    std::visit([&fillFraction, &size](auto&& hist) {
      if (hist) {
        size += HistFiller::getSize(hist, fillFraction);
      }
    },
               mTListValue[j]);
  }
  return size;
}

void TListHandler::clean()
{
  for (auto& value : mTListValue) { // o2-linter: disable=const-ref-in-for-loop (not constant)
    std::visit([](auto&& hist) { hist.reset(); }, value);
  }
}

// print some useful meta-info about the stored histograms
void TListHandler::print(bool showAxisDetails)
{
  std::vector<double> fillFractions{0.1, 0.25, 0.5};
  std::vector<double> totalSizes(fillFractions.size());

  uint32_t nHistos{};
  bool containsSparseHist{};
  auto printHistInfo = [&](auto&& hist) {
    if (hist) {
      using T = std::decay_t<decltype(*hist)>;
      bool isSparse{};
      if (hist->InheritsFrom(THnSparse::Class())) {
        isSparse = true;
        containsSparseHist = true;
      }
      ++nHistos;
      std::vector<double> sizes;
      std::string sizeInfo{};
      if (isSparse) {
        std::transform(std::begin(fillFractions), std::end(fillFractions), std::back_inserter(sizes), [&hist](auto& fraction) { return HistFiller::getSize(hist, fraction); });
        for (uint i = 0; i < fillFractions.size(); ++i) {
          sizeInfo += fmt::format("{:.2f} kB ({:.0f} %)", sizes[i] * 1024, fillFractions[i] * 100);
          if (i != fillFractions.size() - 1) {
            sizeInfo += ", ";
          }
        }
      } else {
        double size = HistFiller::getSize(hist);
        sizes.resize(fillFractions.size(), size);
        sizeInfo = fmt::format("{:.2f} kB", sizes[0] * 1024);
      }
      std::transform(totalSizes.begin(), totalSizes.end(), sizes.begin(), totalSizes.begin(), std::plus<double>());
      LOGF(info, "Hist %03d: %-35s  %-19s [%s]", nHistos, hist->GetName(), hist->IsA()->GetName(), sizeInfo);

      if (showAxisDetails) {
        int nDim = 0;
        if constexpr (std::is_base_of_v<THnBase, T>) {
          nDim = hist->GetNdimensions();
        } else if constexpr (std::is_base_of_v<TH1, T>) {
          nDim = hist->GetDimension();
        }
        TAxis* axis{nullptr};
        for (int d = 0; d < nDim; ++d) {
          if constexpr (std::is_base_of_v<THnBase, T> || std::is_base_of_v<StepTHn, T>) {
            axis = hist->GetAxis(d);
          } else {
            if (d == 0) {
              axis = hist->GetXaxis();
            } else if (d == 1) {
              axis = hist->GetYaxis();
            } else if (d == 2) { // o2-linter: disable=magic-number (not constant in original histogram registry code)
              axis = hist->GetZaxis();
            }
          }
          LOGF(info, "- Axis %d: %-20s (%d bins)", d, axis->GetTitle(), axis->GetNbins());
        }
      }
    }
  };

  std::string titleString{"======================== TListHandler ========================"};
  LOGF(info, "");
  LOGF(info, "%s", titleString);
  LOGF(info, "%s\"%s\"", std::string(static_cast<int>(0.5 * titleString.size() - (1 + 0.5 * mName.size())), ' '), mName);
  for (auto& curHistName : mRegisteredNames) { // o2-linter: disable=const-ref-in-for-loop (not constant in original histogram regitry code)
    std::visit(printHistInfo, mTListValue[getHistIndex(HistName{curHistName.data()})]);
  }
  std::string totalSizeInfo{};
  if (containsSparseHist) {
    for (uint i = 0; i < totalSizes.size(); ++i) {
      totalSizeInfo += fmt::format("{:.2f} MB ({:.0f} %)", totalSizes[i], fillFractions[i] * 100);
      if (i != totalSizes.size() - 1) {
        totalSizeInfo += ", ";
      }
    }
  } else {
    totalSizeInfo = fmt::format("{:.2f} MB", totalSizes[0]);
  }
  LOGF(info, "%s", std::string(titleString.size(), '='), titleString);
  LOGF(info, "Total: %d histograms, ca. %s", nHistos, totalSizeInfo);
  if (lookup) {
    LOGF(info, "Due to index collisions, histograms were shifted by %d registry slots in total.", lookup);
  }
  LOGF(info, "%s", std::string(titleString.size(), '='), titleString);
  LOGF(info, "");
}

template <typename T>
HistPtr TListHandler::insertClone(const HistName& histName, const std::shared_ptr<T> originalHist)
{
  validateHistName(histName.str, histName.hash);
  HistPtr hist = std::shared_ptr<T>(static_cast<T*>(originalHist->Clone(histName.str)));
  performRobinhoodSwapAndGetIndex(histName.hash, hist);
  mRegisteredNames.push_back(histName.str);
  insertInNestedTList(histName.str, hist);
  return hist;
}

template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<TH1>);
template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<TH2>);
template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<TH3>);
template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<TProfile>);
template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<TProfile2D>);
template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<TProfile3D>);
template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<THnSparse>);
template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<THn>);
template HistPtr TListHandler::insertClone(const HistName&, const std::shared_ptr<StepTHn>);

template <typename T>
std::shared_ptr<T> TListHandler::add(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2)
{
  auto histVariant = add(name, title, histConfigSpec, callSumw2);
  if (auto histPtr = std::get_if<std::shared_ptr<T>>(&histVariant)) {
    return *histPtr;
  } else {
    throw runtime_error_f(R"(Histogram type specified in add<>("%s") does not match the actual type of the histogram!)", name);
  }
}

template <typename T>
std::shared_ptr<T> TListHandler::add(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2)
{
  auto histVariant = add(name, title, histType, axes, callSumw2);
  if (auto histPtr = std::get_if<std::shared_ptr<T>>(&histVariant)) {
    return *histPtr;
  } else {
    throw runtime_error_f(R"(Histogram type specified in add<>("%s") does not match the actual type of the histogram!)", name);
  }
}

template std::shared_ptr<TH1> TListHandler::add<TH1>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
template std::shared_ptr<TH1> TListHandler::add<TH1>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);
template std::shared_ptr<TH2> TListHandler::add<TH2>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
template std::shared_ptr<TH2> TListHandler::add<TH2>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);
template std::shared_ptr<TH3> TListHandler::add<TH3>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
template std::shared_ptr<TH3> TListHandler::add<TH3>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);
template std::shared_ptr<TProfile> TListHandler::add<TProfile>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
template std::shared_ptr<TProfile> TListHandler::add<TProfile>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);
template std::shared_ptr<TProfile2D> TListHandler::add<TProfile2D>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
template std::shared_ptr<TProfile2D> TListHandler::add<TProfile2D>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);
template std::shared_ptr<TProfile3D> TListHandler::add<TProfile3D>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
template std::shared_ptr<TProfile3D> TListHandler::add<TProfile3D>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);
template std::shared_ptr<THn> TListHandler::add<THn>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
template std::shared_ptr<THn> TListHandler::add<THn>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);
template std::shared_ptr<THnSparse> TListHandler::add<THnSparse>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
template std::shared_ptr<THnSparse> TListHandler::add<THnSparse>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);
template std::shared_ptr<StepTHn> TListHandler::add<StepTHn>(char const* const name, char const* const title, const HistogramConfigSpec& histConfigSpec, bool callSumw2);
template std::shared_ptr<StepTHn> TListHandler::add<StepTHn>(char const* const name, char const* const title, HistType histType, const std::vector<AxisSpec>& axes, bool callSumw2);

} // namespace o2::framework
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
// TListHandler Class Implementation Over
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

#endif // COMMON_TOOLS_TLISTHANDLER_H_
