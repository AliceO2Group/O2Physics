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

/// \file SecondaryLut.h
/// \brief ROOT-persistent LUT for secondary particles (weak decay products),
///        binned in nch, radius, eta, pt and the particle's number of hits.

#ifndef ALICE3_CORE_SECONDARYLUT_H_
#define ALICE3_CORE_SECONDARYLUT_H_

#include <TNamed.h>
#include <TObject.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <array>
#include <cstddef>
#include <vector>

class TH1D;
class TTree;

namespace o2::fastsim
{

/// Uniform (linear or log10) binning axis. Streamable, trivially copyable.
struct SecondaryLutAxis {
  int nbins = 1;
  float min = 0.f;
  float max = 1.f;
  bool log = false; // min/max are expressed in log10 of the variable

  /// Bin index of val, clamped to [0, nbins-1]
  [[nodiscard]] int find(float val) const;
  /// Bin centre in the variable's natural units
  [[nodiscard]] float center(int bin) const;
  /// Lower edge in the variable's natural units
  [[nodiscard]] float lowEdge(int bin) const;
  /// Upper edge in the variable's natural units
  [[nodiscard]] float upEdge(int bin) const;
  /// Number of bins in the axis
  [[nodiscard]] int getNbins() const { return nbins; }

  ClassDefNV(SecondaryLutAxis, 1);
};

/// Payload of one bin. Plain floats only for compact and fast I/O.
struct SecondaryLutEntry {
  bool valid = false;
  float eff = 0.f;
  float eff2 = 0.f;
  float itof = 0.f;
  float otof = 0.f;
  std::array<float, 15> covm{};
  std::array<float, 5> eigval{};
  std::array<float, 25> eigvec{}; // row-major 5x5
  std::array<float, 25> eiginv{}; // row-major 5x5

  ClassDefNV(SecondaryLutEntry, 1);
};

/// LUT for secondary particles. Entries are stored in one contiguous vector with
/// layout index = ((((nch * nRad + rad) * nEta + eta) * nHits + hits) * nPt + pt),
/// pt varying fastest. Can be written to / read from a ROOT file and inspected
/// with Print(), Browse(), toTree() and projectEfficiency().
class SecondaryLut : public TNamed
{
 public:
  enum Dim : int { kNch = 0,
                   kRad,
                   kEta,
                   kPt,
                   kNHits,
                   kNDims };

  SecondaryLut() = default;
  SecondaryLut(const char* name, int pdg, float mass, float field);

  /// Set the binning of one dimension. Must be called before init().
  void setAxis(Dim dim, int nbins, float min, float max, bool log = false);
  [[nodiscard]] const SecondaryLutAxis& getAxis(Dim dim) const { return mAxes[dim]; }

  /// Allocate (zero-initialised) storage for all bins
  void init();
  [[nodiscard]] bool isInitialized() const { return !mEntries.empty(); }

  [[nodiscard]] int getPdg() const { return mPdg; }
  [[nodiscard]] float getMass() const { return mMass; }
  [[nodiscard]] float getField() const { return mField; }
  [[nodiscard]] std::size_t size() const { return mEntries.size(); }

  [[nodiscard]] std::size_t index(int nch, int rad, int eta, int pt, int nHits) const
  {
    const auto& a = mAxes;
    return (((static_cast<std::size_t>(nch) * a[kRad].nbins + rad) * a[kEta].nbins + eta) * a[kNHits].nbins + nHits) * a[kPt].nbins + pt;
  }

  /// Direct bin access (no bounds clamping)
  SecondaryLutEntry& at(int nch, int rad, int eta, int pt, int nHits)
  {
    return mEntries[index(nch, rad, eta, pt, nHits)];
  }
  [[nodiscard]] const SecondaryLutEntry& at(int nch, int rad, int eta, int pt, int nHits) const
  {
    return mEntries[index(nch, rad, eta, pt, nHits)];
  }

  /// Lookup by physical values; out-of-range values are clamped to edge bins.
  /// Returns nullptr if the LUT is not initialised.
  [[nodiscard]] const SecondaryLutEntry* lookup(float nch, float rad, float eta, float pt, float nHits) const;

  [[nodiscard]] const std::vector<SecondaryLutEntry>& entries() const { return mEntries; }

  /// Write to a ROOT file (object stored under its own name)
  bool writeToFile(const char* fileName, const char* option = "RECREATE") const;
  /// Read from a ROOT file; returns nullptr on failure. Caller owns the result.
  static SecondaryLut* readFromFile(const char* fileName, const char* objName = "SecondaryLut");

  /// Flat TTree with one row per bin (bin indices, bin centres and payload)
  [[nodiscard]] TTree* toTree(const char* treeName = "secondaryLut") const;
  /// Efficiency (average over valid bins) projected on one dimension. Caller owns the result.
  [[nodiscard]] TH1D* projectEfficiency(Dim dim, const char* histName = nullptr) const;

  void Print(Option_t* option = "") const override;
  void Browse(TBrowser* b) override;
  [[nodiscard]] Bool_t IsFolder() const override { return kTRUE; }

 private:
  void unravel(std::size_t idx, std::array<int, kNDims>& bins) const;

  int mPdg = 0;
  float mMass = 0.f;
  float mField = 0.f;
  std::array<SecondaryLutAxis, kNDims> mAxes;
  std::vector<SecondaryLutEntry> mEntries;

  ClassDefOverride(SecondaryLut, 1);
};

} // namespace o2::fastsim

#endif // ALICE3_CORE_SECONDARYLUT_H_
