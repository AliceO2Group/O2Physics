// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICE3_CORE_FLATLUTROOTWRITER_H_
#define ALICE3_CORE_FLATLUTROOTWRITER_H_

#include "ALICE3/Core/FlatLutEntry.h"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace o2::fastsim
{

struct FlatLutRootRow {
  int entryIndex = -1;
  int nchBin = -1;
  int radBin = -1;
  int etaBin = -1;
  int ptBin = -1;

  float nch = 0.f;
  float eta = 0.f;
  float pt = 0.f;

  bool valid = false;
  float eff = 0.f;
  float eff2 = 0.f;
  float itof = 0.f;
  float otof = 0.f;

  std::array<float, kNumCovarianceTerms> covm{};
  std::array<float, kNumEigenModes> eigval{};
  EigenMatrix eigvec{};
  EigenMatrix eiginv{};
};

class FlatLutRootWriter
{
 public:
  FlatLutRootWriter(const lutHeader_t& header, bool covarianceStoresSigma);

  void reserve(std::size_t n);
  void addEntry(int entryIndex,
                int nchBin,
                int radBin,
                int etaBin,
                int ptBin,
                const lutEntry_t& entry);

  bool writeToFile(const std::string& rootFileName) const;

 private:
  lutHeader_t mHeader;
  bool mCovarianceStoresSigma = false;
  std::vector<FlatLutRootRow> mRows;
};

} // namespace o2::fastsim

#endif // ALICE3_CORE_FLATLUTROOTWRITER_H_
