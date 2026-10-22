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

/// \file upcTestFitBitMapping.cxx
/// \brief Tutorial for accessing and mapping the FIT-bit table
/// \author Sandor Lokos, sandor.lokos@cern.ch
/// \since October 2026

#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

#include <CommonConstants/MathConstants.h>
#include <FT0Base/Geometry.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <Framework/runDataProcessing.h>

#include <array>

using namespace o2;
using namespace o2::framework;

struct UpcTestFitBitMapping {
  static constexpr int kThr1Selector = 1;
  static constexpr int kThr2Selector = 2;

  Configurable<int> whichThr{"whichThr", 1, "FIT amplitude threshold to use: 1 for Thr1 or 2 for Thr2"};

  // Minimal zero-offset container required by getPhiEtaFromFitBit().
  // Detector-alignment corrections are not applied in this tutorial.
  struct OffsetXYZ {
    double x{0.};
    double y{0.};
    double z{0.};

    double getX() const { return x; }
    double getY() const { return y; }
    double getZ() const { return z; }
  };

  std::array<OffsetXYZ, 1> offsetFT0{};
  static constexpr int kRunOffsetIndex = 0;

  o2::ft0::Geometry ft0Geometry{};

  HistogramRegistry registry{
    "registry",
    {
      {"association/hRowsPerCollision",
       "FIT-bit rows associated with each collision;rows;collisions",
       {HistType::kTH1F, {{4, -0.5, 3.5}}}},

      {"multiplicity/hFT0A",
       "FT0A fired-channel multiplicity;N_{fired}^{FT0A};collisions",
       {HistType::kTH1F, {{97, -0.5, 96.5}}}},
      {"multiplicity/hFT0C",
       "FT0C fired-channel multiplicity;N_{fired}^{FT0C};collisions",
       {HistType::kTH1F, {{113, -0.5, 112.5}}}},
      {"multiplicity/hFV0A",
       "FV0A fired-channel multiplicity;N_{fired}^{FV0A};collisions",
       {HistType::kTH1F, {{49, -0.5, 48.5}}}},
      {"multiplicity/hFT0AVsFT0C",
       "FT0A versus FT0C fired channels;N_{fired}^{FT0A};N_{fired}^{FT0C}",
       {HistType::kTH2F, {{97, -0.5, 96.5}, {113, -0.5, 112.5}}}},

      {"occupancy/hPackedBit",
       "Selected FIT-bit occupancy;packed bit;counts",
       {HistType::kTH1F, {{256, -0.5, 255.5}}}},
      {"occupancy/hFT0AChannel",
       "FT0A channel occupancy;FT0A channel;counts",
       {HistType::kTH1F, {{96, -0.5, 95.5}}}},
      {"occupancy/hFT0CChannel",
       "FT0C channel occupancy;FT0C channel;counts",
       {HistType::kTH1F, {{112, -0.5, 111.5}}}},
      {"occupancy/hFV0AChannel",
       "FV0A channel occupancy;FV0A channel;counts",
       {HistType::kTH1F, {{48, -0.5, 47.5}}}},

      {"mapping/hPhiA",
       "FT0A #varphi;#varphi;counts",
       {HistType::kTH1F, {{18, 0., o2::constants::math::TwoPI}}}},
      {"mapping/hEtaA",
       "FT0A #eta;#eta;counts",
       {HistType::kTH1F, {{8, 3.5, 5.0}}}},
      {"mapping/hEtaPhiA",
       "FT0A #eta versus #varphi;#eta;#varphi",
       {HistType::kTH2F, {{8, 3.5, 5.0}, {18, 0., o2::constants::math::TwoPI}}}},
      {"mapping/hXYA",
       "FT0A channel positions;x (cm);y (cm)",
       {HistType::kTH2F, {{24, -18., 18.}, {24, -18., 18.}}}},

      {"mapping/hPhiC",
       "FT0C #varphi;#varphi;counts",
       {HistType::kTH1F, {{18, 0., o2::constants::math::TwoPI}}}},
      {"mapping/hEtaC",
       "FT0C #eta;#eta;counts",
       {HistType::kTH1F, {{8, -3.5, -2.0}}}},
      {"mapping/hEtaPhiC",
       "FT0C #eta versus #varphi;#eta;#varphi",
       {HistType::kTH2F, {{8, -3.5, -2.0}, {18, 0., o2::constants::math::TwoPI}}}},
      {"mapping/hXYC",
       "FT0C channel positions;x (cm);y (cm)",
       {HistType::kTH2F, {{24, -18., 18.}, {24, -18., 18.}}}},
    }};

  void init(InitContext&)
  {
    if (whichThr.value != kThr1Selector && whichThr.value != kThr2Selector) {
      LOGF(fatal,
           "Invalid whichThr=%d: allowed values are 1 and 2",
           whichThr.value);
    }

    ft0Geometry.calculateChannelCenter();
  }

  template <typename T>
  udhelpers::Bits256 getSelectedBits(T const& row) const
  {
    if (whichThr.value == kThr2Selector) {
      return udhelpers::makeBits256(
        row.thr2W0(), row.thr2W1(), row.thr2W2(), row.thr2W3());
    }

    return udhelpers::makeBits256(
      row.thr1W0(), row.thr1W1(), row.thr1W2(), row.thr1W3());
  }

  void process(aod::UDCollisions::iterator const&,
               aod::UDCollisionFITBits const& fitBits)
  {
    // Because UDCollisionFITBits contains UDCollisionId, the framework
    // supplies only the FIT-bit rows associated with the current collision.
    registry.fill(HIST("association/hRowsPerCollision"), fitBits.size());

    // The producer is expected to write exactly one FIT-bit row per collision.
    if (fitBits.size() != 1) {
      return;
    }

    const auto row = fitBits.begin();
    const auto bits = getSelectedBits(row);

    int nFT0A = 0;
    int nFT0C = 0;
    int nFV0A = 0;

    // Packed layout:
    //   0--95:    FT0A
    //   96--207:  FT0C
    //   208--255: FV0A
    for (int bit = 0; bit < udhelpers::kTotalBits; ++bit) {
      if (!udhelpers::testBit(bits, bit)) {
        continue;
      }

      registry.fill(HIST("occupancy/hPackedBit"), bit);

      const auto decoded = udhelpers::decodeFitBit(bit);

      if (decoded.det == udhelpers::FitBitRef::Det::FV0) {
        ++nFV0A;
        registry.fill(HIST("occupancy/hFV0AChannel"), decoded.ch);
        continue;
      }

      if (decoded.det != udhelpers::FitBitRef::Det::FT0) {
        continue;
      }

      if (decoded.isC) {
        ++nFT0C;
        registry.fill(
          HIST("occupancy/hFT0CChannel"),
          decoded.ch - udhelpers::kFT0AChannels);
      } else {
        ++nFT0A;
        registry.fill(HIST("occupancy/hFT0AChannel"), decoded.ch);
      }

      double phi = 0.;
      double eta = 0.;
      if (!udhelpers::getPhiEtaFromFitBit(
            ft0Geometry, bit, offsetFT0, kRunOffsetIndex, phi, eta)) {
        continue;
      }

      const auto position = ft0Geometry.getChannelCenter(decoded.ch);
      const double x = position.X() + offsetFT0[kRunOffsetIndex].getX();
      const double y = position.Y() + offsetFT0[kRunOffsetIndex].getY();

      if (decoded.isC) {
        registry.fill(HIST("mapping/hPhiC"), phi);
        registry.fill(HIST("mapping/hEtaC"), eta);
        registry.fill(HIST("mapping/hEtaPhiC"), eta, phi);
        registry.fill(HIST("mapping/hXYC"), x, y);
      } else {
        registry.fill(HIST("mapping/hPhiA"), phi);
        registry.fill(HIST("mapping/hEtaA"), eta);
        registry.fill(HIST("mapping/hEtaPhiA"), eta, phi);
        registry.fill(HIST("mapping/hXYA"), x, y);
      }
    }

    registry.fill(HIST("multiplicity/hFT0A"), nFT0A);
    registry.fill(HIST("multiplicity/hFT0C"), nFT0C);
    registry.fill(HIST("multiplicity/hFV0A"), nFV0A);
    registry.fill(HIST("multiplicity/hFT0AVsFT0C"), nFT0A, nFT0C);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<UpcTestFitBitMapping>(cfgc)};
}
