// UpcTestFITBitMapping.cxx
// Minimal task to test FT0 bitset -> (eta,phi) mapping using udhelpers::getPhiEtaFromFitBit
//
// Run example (from O2Physics build dir):
//   ./stage/bin/o2-pwgud-ud-fitbit-mapping --aod-file AO2D_UD.root -b
//
// Note: Add this file to PWGUD/Tasks/CMakeLists.txt.

#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"

#include "PWGUD/DataModel/UDTables.h"     // aod::UDCollisionFT0Bits
#include "PWGUD/Core/UDHelpers.h"         // udhelpers::Bits256, makeBits256, testBit, getPhiEtaFromFitBit
#include "FT0Base/Geometry.h"             // o2::ft0::Geometry

#include <array>
#include <cmath>

using namespace o2;
using namespace o2::framework;

struct UpcTestFITBitMapping
{
  Configurable<int> whichThr{"whichThr", 1, "Use 1=Thr1 bits or 2=Thr2 bits"};
  Configurable<int> maxEvents{"maxEvents", -1, "Process at most this many rows (-1 = all)"};

  // Minimal offset container compatible with UDHelpers.h expectations: getX/getY/getZ
  struct OffsetXYZ {
    double x{0}, y{0}, z{0};
    double getX() const { return x; }
    double getY() const { return y; }
    double getZ() const { return z; }
  };

  std::array<OffsetXYZ, 1> offsetFT0{}; // iRunOffset = 0 for now
  int iRunOffset = 0;

  o2::ft0::Geometry ft0Det{};

  HistogramRegistry registry{
    "registry",
    {
      {"hPhiA",    "FT0A #varphi;#varphi;counts",                {HistType::kTH1F, {{18, 0.0,  2.*M_PI}}}},
      {"hEtaA",    "FT0A #eta;#eta;counts",                      {HistType::kTH1F, {{8, 3.5, 5.0}}}},
      {"hEtaPhiA", "FT0A #eta vs #varphi;#eta;#varphi",          {HistType::kTH2F, {{8, 3.5, 5.0}, {18, 0.0,  2.*M_PI}}}},

      {"hPhiC",    "FT0C #varphi;#varphi;counts",                {HistType::kTH1F, {{18, 0.0,  2.*M_PI}}}},
      {"hEtaC",    "FT0C #eta;#eta;counts",                      {HistType::kTH1F, {{8, -3.5, -2.0}}}},
      {"hEtaPhiC", "FT0C #eta vs #varphi;#eta;#varphi",          {HistType::kTH2F, {{8, -3.5, -2.0}, {18, 0.0,  2.*M_PI}}}},
    }};

  void init(InitContext&)
  {
    // UDHelpers calls calculateChannelCenter() inside, but doing it once here is fine.
    ft0Det.calculateChannelCenter();
  }

  void process(aod::UDCollisionFT0Bits const& bitsTable)
  {
    int64_t nProcessed = 0;

    for (auto const& row : bitsTable) {
      if (maxEvents >= 0 && nProcessed >= maxEvents) {
        break;
      }
      ++nProcessed;

      // Use udhelpers' canonical packed type + builder
      udhelpers::Bits256 w{};
      if (whichThr == 2) {
        w = udhelpers::makeBits256(row.thr2W0(), row.thr2W1(), row.thr2W2(), row.thr2W3());
      } else {
        w = udhelpers::makeBits256(row.thr1W0(), row.thr1W1(), row.thr1W2(), row.thr1W3());
      }

      // Loop FT0 bits only (0..207). FV0 starts at 208 but ignored here.
      for (int bit = 0; bit < udhelpers::kFT0Bits; ++bit) {
        if (!udhelpers::testBit(w, bit)) {
          continue;
        }

        double phi = 0., eta = 0.;
        const bool ok = udhelpers::getPhiEtaFromFitBit(ft0Det, bit, offsetFT0, iRunOffset, phi, eta);
        if (!ok) {
          continue;
        }

        if (bit < udhelpers::kFT0AChannels) {
          registry.fill(HIST("hPhiA"), phi);
          registry.fill(HIST("hEtaA"), eta);
          registry.fill(HIST("hEtaPhiA"), eta, phi);
        } else {
          registry.fill(HIST("hPhiC"), phi);
          registry.fill(HIST("hEtaC"), eta);
          registry.fill(HIST("hEtaPhiC"), eta, phi);
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcTestFITBitMapping>(cfgc, TaskName{"fitbit-mapping"})
  };
}
