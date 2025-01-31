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

/// \file QualitySelectionParams.h
/// \brief Quality selection parameters
///
/// \author Andrea Ferrero <andrea.ferrero@cern.ch> and Evgeny Kryshen <evgeny.kryshen@cern.ch>

#ifndef COMMON_CCDB_QUALITYSELECTIONPARAMS_H_
#define COMMON_CCDB_QUALITYSELECTIONPARAMS_H_

#include <Rtypes.h>
#include <TMath.h>

#include <stdexcept>
#include <algorithm>
#include <string>
#include <vector>

namespace o2::aod::qualitysel
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

// Quality selection flags
enum QualitySelectionFlags {
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
  kNQualitySelectionFlags
};

template <typename T>
concept HasDetectorQuality = requires(T a, int bit)
{
  { a.qc_bit(bit) } -> std::convertible_to<bool>;
};

class QualityFlagsChecker
{
  // Helper class to allow appending new elements to a vector of QualitySelectionFlags
  class QualitySelection: public std::vector<QualitySelectionFlags>
  {
  public:
    QualitySelection& operator+=(const std::initializer_list<QualitySelectionFlags>& l)
    {
      insert(end(), l.begin(), l.end());
      return (*this);
    }

    QualitySelection& operator+=(const std::vector<QualitySelectionFlags>& v)
    {
      insert(end(), v.begin(), v.end());
      return (*this);
    }
  };

public:
  QualityFlagsChecker() = default;

  // construct the object from an initializer list, like this:
  // QualityFlagsChecker qualityFlagsChecker{ kFT0Bad, kITSBad, kMFTBad, kMFTLimAccMCRepr };
  QualityFlagsChecker(std::initializer_list<QualitySelectionFlags> bitsToCheck) : vBitsToCheck(bitsToCheck.size())
  {
    std::copy(bitsToCheck.begin(), bitsToCheck.end(), vBitsToCheck.begin());
  }

  // construct the object from a vector of QualitySelectionFlags
  explicit QualityFlagsChecker(const std::vector<QualitySelectionFlags> bitsToCheck) : vBitsToCheck(bitsToCheck) {}

  // Construct the object from one of the pre-defined official runlist selections.
  // The label parameter can take the following values:
  // - "CBT"
  // - "CBT_hadronPID"
  // - "CBT_electronPID"
  // - "CCBT_caloBT"
  // - "CBT_muon"
  // - "CBT_muon_glo"
  // The checkZDC boolean flag controls whether to iclude the ZDC quality in all the pre-defined selections (for Pb-Pb data)
  // The treatLimitedAcceptanceAsBad boolean flag controls whether "LimitedAcceptanceMCReproducible" flags should be
  // treated as Bad and the corresponding events excluded
  QualityFlagsChecker(const std::string& label, bool checkZDC = false, bool treatLimitedAcceptanceAsBad = false)  // o2-linter: disable=runtime/explicit
  {
    initialize(label, checkZDC, treatLimitedAcceptanceAsBad);
  }

  QualityFlagsChecker(const char* label, bool checkZDC = false, bool treatLimitedAcceptanceAsBad = false)  // o2-linter: disable=runtime/explicit
  {
    initialize(std::string(label), checkZDC, treatLimitedAcceptanceAsBad);
  }

  void initialize(const std::string& label, bool checkZDC = false, bool treatLimitedAcceptanceAsBad = false)
  {
    QualitySelection bitsToCheck;

    if (label == "CBT") {
      bitsToCheck += { kFT0Bad, kITSBad, kTPCBadTracking, kTPCBadPID };
      if (treatLimitedAcceptanceAsBad) {
        bitsToCheck += { kITSLimAccMCRepr, kTPCLimAccMCRepr };
      }
    }

    if (label == "CBT_hadronPID") {
      bitsToCheck += {kFT0Bad, kITSBad, kTPCBadTracking, kTPCBadPID, kTOFBad};
      if (treatLimitedAcceptanceAsBad) {
        bitsToCheck += {kITSLimAccMCRepr, kTPCLimAccMCRepr, kTOFLimAccMCRepr};
      }
    }

    if (label == "CBT_electronPID") {
      bitsToCheck += {kFT0Bad, kITSBad, kTPCBadTracking, kTPCBadPID, kTRDBad};
      if (treatLimitedAcceptanceAsBad) {
        bitsToCheck += {kITSLimAccMCRepr, kTPCLimAccMCRepr};
      }
    }

    if (label == "CBT_calo") {
      bitsToCheck += {kFT0Bad, kITSBad, kTPCBadTracking, kTPCBadPID, kEMCBad};
      if (treatLimitedAcceptanceAsBad) {
        bitsToCheck += {kITSLimAccMCRepr, kTPCLimAccMCRepr, kEMCLimAccMCRepr};
      }
    }

    if (label == "CBT_muon") {
      bitsToCheck += {kFT0Bad, kITSBad, kTPCBadTracking, kMCHBad, kMIDBad};
      if (treatLimitedAcceptanceAsBad) {
        bitsToCheck += {kITSLimAccMCRepr, kTPCLimAccMCRepr, kMCHLimAccMCRepr, kMIDLimAccMCRepr};
      }
    }

    if (label == "CBT_muon_glo") {
      bitsToCheck += {kFT0Bad, kITSBad, kTPCBadTracking, kMCHBad, kMIDBad, kMFTBad};
      if (treatLimitedAcceptanceAsBad) {
        bitsToCheck += {kITSLimAccMCRepr, kTPCLimAccMCRepr, kMCHLimAccMCRepr, kMIDLimAccMCRepr, kMFTLimAccMCRepr};
      }
    }

    if (checkZDC) {
      bitsToCheck += { kZDCBad };
    }

    vBitsToCheck = bitsToCheck;
  }

  bool checkTable(const HasDetectorQuality auto& table)
  {
    if (vBitsToCheck.empty()) {
      throw std::out_of_range ("QualityFlagsChecker with empty QualitySelectionFlags bits vector");
    }

    for (auto bit : vBitsToCheck) {
      if (table.qc_bit(bit)) {
        return false;
      }
    }
    return true;
  }

  bool operator ()(const HasDetectorQuality auto& table)
  {
    return checkTable(table);
  }

private:
  std::vector<QualitySelectionFlags> vBitsToCheck;
};

} // namespace o2::aod::qualitysel
#endif // COMMON_CCDB_QUALITYSELECTIONPARAMS_H_
