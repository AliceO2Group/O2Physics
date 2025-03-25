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

/// \file RCTSelectionFlags.h
/// \brief RCT selection flags
///
/// \author Andrea Ferrero <andrea.ferrero@cern.ch> and Evgeny Kryshen <evgeny.kryshen@cern.ch>

#ifndef COMMON_CCDB_RCTSELECTIONFLAGS_H_
#define COMMON_CCDB_RCTSELECTIONFLAGS_H_

#include <CommonUtils/EnumFlags.h>
#include <Rtypes.h>
#include <TMath.h>

#include <stdexcept>
#include <algorithm>
#include <string>
#include <vector>

namespace o2::aod::rctsel
{
/*
 * Bit mapping used for populating the CCDB objects from the RCT flags
 * From https://github.com/JianLIUhep/RCTutils/blob/main/CCDB/process_and_upload.C
std::map<std::string, std::map<std::string, int>> detailedBitMapping = {
  {"CPV", { {"Bad", 0}, {"Invalid", 0} }},
  {"EMC", { {"Bad", 1}, {"NoDetectorData", 1}, {"BadEMCalorimetry", 1}, {"LimitedAcceptanceMCReproducible", 2} }},
  {"FDD", { {"Bad", 3}, {"Invalid", 3}, {"NoDetectorData", 3} }},
  {"FT0", { {"Bad", 4}, {"UnknownQuality", 4}, {"Unknown", 4} }},
  {"FV0", { {"Bad", 5} }},
  {"HMP", { {"Bad", 6}, {"NoDetectorData", 6} }},
  {"ITS", { {"Bad", 7}, {"UnknownQuality", 7}, {"BadTracking", 7}, {"LimitedAcceptanceMCReproducible", 8} }},
  {"MCH", { {"Bad", 9}, {"NoDetectorData", 9}, {"Unknown", 9}, {"LimitedAcceptanceMCReproducible", 10} }},
  {"MFT", { {"Bad", 11}, {"BadTracking", 11}, {"LimitedAcceptanceMCReproducible", 12} }},
  {"MID", { {"Bad", 13}, {"BadTracking", 13}, {"LimitedAcceptanceMCReproducible", 14} }},
  {"PHS", { {"Bad", 15}, {"Invalid", 15} }},
  {"TOF", { {"Bad", 16}, {"NoDetectorData", 16}, {"BadPID", 16}, {"LimitedAcceptanceMCReproducible", 17} }},
  {"TPC", { {"Bad", 18}, {"BadTracking", 18}, {"BadPID", 19}, {"LimitedAcceptanceMCNotReproducible", 18}, {"LimitedAcceptanceMCReproducible", 20} }},
  {"TRD", { {"Bad", 21}, {"BadTracking", 21} }},
  {"ZDC", { {"Bad", 22}, {"UnknownQuality", 22}, {"Unknown", 22}, {"NoDetectorData", 22} }}
};
*/

// RCT selection flags
enum RCTSelectionFlags {
  kCPVBad = 0,
  kEMCBad,
  kEMCLimAccMCRepr,
  kFDDBad,
  kFT0Bad,
  kFV0Bad,
  kHMPBad,
  kITSBad,
  kITSLimAccMCRepr,
  kMCHBad,
  kMCHLimAccMCRepr,
  kMFTBad,
  kMFTLimAccMCRepr,
  kMIDBad,
  kMIDLimAccMCRepr,
  kPHSBad,
  kTOFBad,
  kTOFLimAccMCRepr,
  kTPCBadTracking,
  kTPCBadPID,
  kTPCLimAccMCRepr,
  kTRDBad,
  kZDCBad,
  kNRCTSelectionFlags
};

template <typename T>
concept HasRCTFlags = requires(T a, int bit) {
  { a.rct_bit(bit) } -> std::convertible_to<bool>;
  { a.rct_raw() } -> std::convertible_to<uint64_t>;
};

class RCTFlagsChecker : public o2::utils::EnumFlags<RCTSelectionFlags>
{
 public:
  RCTFlagsChecker() = default;

  // Construct the object from an initializer list, like this:
  // RCTFlagsChecker qualityFlagsChecker{ kFT0Bad, kITSBad, kMFTBad, kMFTLimAccMCRepr };
  using o2::utils::EnumFlags<RCTSelectionFlags>::EnumFlags;

  // Construct the object from one of the pre-defined runlist selections.
  // The label parameter can take the following values:
  // - "CBT"
  // - "CBT_hadronPID"
  // - "CBT_electronPID"
  // - "CCBT_calo"
  // - "CBT_muon"
  // - "CBT_muon_glo"
  // The checkZDC boolean flag controls whether to iclude the ZDC quality in all the pre-defined selections (for Pb-Pb data)
  // The treatLimitedAcceptanceAsBad boolean flag controls whether "LimitedAcceptanceMCReproducible" flags should be
  // treated as Bad and the corresponding events excluded
  explicit RCTFlagsChecker(const std::string& label, bool checkZDC = false, bool treatLimitedAcceptanceAsBad = false)
  {
    init(label, checkZDC, treatLimitedAcceptanceAsBad);
  }

  // Initialize the object from an initializer list of RCTSelectionFlags values
  void init(std::initializer_list<RCTSelectionFlags> flags)
  {
    reset();
    *this = RCTFlagsChecker(flags);
  }

  // Initialize the object from one of the pre-defined runlist selections.
  // The label parameter can take the following values:
  // - "CBT"
  // - "CBT_hadronPID"
  // - "CBT_electronPID"
  // - "CCBT_calo"
  // - "CBT_muon"
  // - "CBT_muon_glo"
  // The checkZDC boolean flag controls whether to iclude the ZDC quality in all the pre-defined selections (for Pb-Pb data)
  // The treatLimitedAcceptanceAsBad boolean flag controls whether "LimitedAcceptanceMCReproducible" flags should be
  // treated as Bad and the corresponding events excluded
  void init(const std::string& label, bool checkZDC = false, bool treatLimitedAcceptanceAsBad = false)
  {
    auto setFlags = [this](std::initializer_list<RCTSelectionFlags> flags) {
      std::for_each(flags.begin(),
                    flags.end(),
                    [this](const RCTSelectionFlags f) noexcept { set(f); });
    };

    reset();

    if (label == "CBT") {
      setFlags({kFT0Bad, kITSBad, kTPCBadTracking, kTPCBadPID});
      if (treatLimitedAcceptanceAsBad) {
        setFlags({kITSLimAccMCRepr, kTPCLimAccMCRepr});
      }
    }

    if (label == "CBT_hadronPID") {
      setFlags({kFT0Bad, kITSBad, kTPCBadTracking, kTPCBadPID, kTOFBad});
      if (treatLimitedAcceptanceAsBad) {
        setFlags({kITSLimAccMCRepr, kTPCLimAccMCRepr, kTOFLimAccMCRepr});
      }
    }

    if (label == "CBT_electronPID") {
      setFlags({kFT0Bad, kITSBad, kTPCBadTracking, kTPCBadPID, kTRDBad});
      if (treatLimitedAcceptanceAsBad) {
        setFlags({kITSLimAccMCRepr, kTPCLimAccMCRepr});
      }
    }

    if (label == "CBT_calo") {
      setFlags({kFT0Bad, kITSBad, kTPCBadTracking, kTPCBadPID, kEMCBad});
      if (treatLimitedAcceptanceAsBad) {
        setFlags({kITSLimAccMCRepr, kTPCLimAccMCRepr, kEMCLimAccMCRepr});
      }
    }

    if (label == "CBT_muon") {
      setFlags({kFT0Bad, kITSBad, kTPCBadTracking, kMCHBad, kMIDBad});
      if (treatLimitedAcceptanceAsBad) {
        setFlags({kITSLimAccMCRepr, kTPCLimAccMCRepr, kMCHLimAccMCRepr, kMIDLimAccMCRepr});
      }
    }

    if (label == "CBT_muon_glo") {
      setFlags({kFT0Bad, kITSBad, kTPCBadTracking, kMCHBad, kMIDBad, kMFTBad});
      if (treatLimitedAcceptanceAsBad) {
        setFlags({kITSLimAccMCRepr, kTPCLimAccMCRepr, kMCHLimAccMCRepr, kMIDLimAccMCRepr, kMFTLimAccMCRepr});
      }
    }

    if (checkZDC) {
      set(kZDCBad);
    }
  }

  // Check the RCT column of a given event selection table.
  // The function returns true if none of the checked flags is set in the RCT column.
  bool checkTable(const HasRCTFlags auto& table)
  {
    if (!any()) {
      throw std::out_of_range("RCTFlagsCheckerAlt with empty RCTSelectionFlags bits mask");
    }

    // bitmask of the current table
    uint64_t tableBits = table.rct_raw();
    // bitmask of flags to be checked
    uint64_t flagsBits = value();

    // return true if none of the checked bits is set in the table bitmask
    return ((tableBits & flagsBits) == 0);
  }

  bool operator()(const HasRCTFlags auto& table)
  {
    return checkTable(table);
  }
};

} // namespace o2::aod::rctsel
#endif // COMMON_CCDB_RCTSELECTIONFLAGS_H_
