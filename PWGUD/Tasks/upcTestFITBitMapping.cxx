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
// \brief FIT bits to phi, eta mapping + multiplicity from FIT bits
// \author Sandor Lokos, sandor.lokos@cern.ch
// \since March 2026

#include "PWGUD/Core/SGSelector.h"
#include "PWGUD/Core/UDHelpers.h"
#include "PWGUD/DataModel/UDTables.h"

#include "FT0Base/Geometry.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"
#include "Framework/runDataProcessing.h"

#include <array>
#include <cmath>
#include <cstdint>

using namespace o2;
using namespace o2::framework;

// using UDCollisionsFull = soa::Join<
  // aod::UDCollisions,
  // aod::UDCollisionsSels,
  // aod::UDCollisionSelExtras>;

using UDCollisionsFull = aod::UDCollisions;

struct UpcTestFITBitMapping {
	
	SGSelector sgSelector;

  Configurable<int>  whichThr{"whichThr", 1, "Use 1=Thr1 bits or 2=Thr2 bits"};
  Configurable<int>  maxEvents{"maxEvents", -1, "Process at most this many collisions (-1 = all)"};
  Configurable<int>  debugPrintEvery{"debugPrintEvery", 1000, "Print every N events (-1 disables)"};
  Configurable<int>  debugPrintFirst{"debugPrintFirst", 10, "Always print first N events"};
  Configurable<bool> useRctFlag{"useRctFlag", false, "use RCT flags for event selection"};
  Configurable<int>  cutRctFlag{"cutRctFlag", 0, "0 = off, 1 = CBT, 2 = CBT+ZDC, 3 = CBThadron, 4 = CBThadron+ZDC"};

  struct OffsetXYZ {
    double x{0.}, y{0.}, z{0.};
    double getX() const { return x; }
    double getY() const { return y; }
    double getZ() const { return z; }
  };

  std::array<OffsetXYZ, 1> offsetFT0{};
  int iRunOffset = 0;

  o2::ft0::Geometry ft0Det{};

  HistogramRegistry registry{
    "registry",
    {
      {"debug/hEventCounter", "Event counter;step;events", {HistType::kTH1F, {{6, -0.5, 5.5}}}},
      {"debug/hFitBitsSize", "fitBits.size() per collision;fitBits.size();events", {HistType::kTH1F, {{10, -0.5, 9.5}}}},
      {"debug/hCollisionIndexMod", "Collision index mod 100;collision.globalIndex() % 100;events", {HistType::kTH1F, {{100, -0.5, 99.5}}}},

      {"map/hPhiA", "FT0A #varphi;#varphi;counts", {HistType::kTH1F, {{18, 0.0, 2. * M_PI}}}},
      {"map/hEtaA", "FT0A #eta;#eta;counts", {HistType::kTH1F, {{8, 3.5, 5.0}}}},
      {"map/hEtaPhiA", "FT0A #eta vs #varphi;#eta;#varphi", {HistType::kTH2F, {{8, 3.5, 5.0}, {18, 0.0, 2. * M_PI}}}},

      {"map/hPhiC", "FT0C #varphi;#varphi;counts", {HistType::kTH1F, {{18, 0.0, 2. * M_PI}}}},
      {"map/hEtaC", "FT0C #eta;#eta;counts", {HistType::kTH1F, {{8, -3.5, -2.0}}}},
      {"map/hEtaPhiC", "FT0C #eta vs #varphi;#eta;#varphi", {HistType::kTH2F, {{8, -3.5, -2.0}, {18, 0.0, 2. * M_PI}}}},

      {"mult/hPnFT0A", "P(n): FT0A fired-channel multiplicity;N_{fired}^{FT0A};events", {HistType::kTH1F, {{97, -0.5, 96.5}}}},
      {"mult/hPnFT0C", "P(n): FT0C fired-channel multiplicity;N_{fired}^{FT0C};events", {HistType::kTH1F, {{113, -0.5, 112.5}}}},
      {"mult/hPnFT0", "P(n): FT0 fired-channel multiplicity;N_{fired}^{FT0};events", {HistType::kTH1F, {{209, -0.5, 208.5}}}},
      {"mult/hPnFV0A", "P(n): FV0A fired-channel multiplicity;N_{fired}^{FV0A};events", {HistType::kTH1F, {{49, -0.5, 48.5}}}},
      {"mult/hPnFIT", "P(n): FIT fired-channel multiplicity;N_{fired}^{FIT};events", {HistType::kTH1F, {{257, -0.5, 256.5}}}},

      {"mult/hNfiredA_vs_C", "FT0A vs FT0C fired channels;N_{fired}^{FT0A};N_{fired}^{FT0C}",
       {HistType::kTH2F, {{97, -0.5, 96.5}, {113, -0.5, 112.5}}}}
    }};

  int countBitsInRange(udhelpers::Bits256 const& bits, int first, int last) const
  {
    int n = 0;
    for (int i = first; i <= last; ++i) {
      if (udhelpers::testBit(bits, i)) {
        ++n;
				}
    }
    return n;
  }
	
	// template <typename C>
	// bool isGoodRctFlag(const C& collision)
	// {
		// switch (cutRctFlag) {
			// case 1: // CBT
				// return sgSelector.isCBTOk(collision);
			// case 2: // CBT + ZDC
				// return sgSelector.isCBTZdcOk(collision);
			// case 3: // CBT hadron
				// return sgSelector.isCBTHadronOk(collision);
			// case 4: // CBT hadron + ZDC
				// return sgSelector.isCBTHadronZdcOk(collision);
			// default: // no RCT cut applied
				// return true;
		// }
	// }
	
	
	// template <typename C>
	// bool collisionPassesCuts(const C& collision)
	// {
		// /* good vertex */
		// if (!collision.vtxITSTPC()) {
			// return false;
		// }
		// registry.fill(HIST("debug/hEventCounter"), 4.);

		// /* same bunch pile-up rejection */
		// if (!collision.sbp()) {
			// return false;
		// }
		// registry.fill(HIST("debug/hEventCounter"), 5.);

		// /* ITS ROF rejection */
		// if (!collision.itsROFb()) {
			// return false;
		// }
		// registry.fill(HIST("debug/hEventCounter"), 6.);

		// /* Timeframe border collision rejection */
		// if (!collision.tfb()) {
			// return false;
		// }
		// registry.fill(HIST("debug/hEventCounter"), 7.);

		// /* z-vertex cut */
		// if (std::abs(collision.posZ()) > 10.0) {
			// return false;
		// }
		// registry.fill(HIST("debug/hEventCounter"), 8.);
		
		// /* RCT flag check*/
		// if (useRctFlag) {
			// if (!isGoodRctFlag(collision)) {
				// return false;
			// }
			// registry.fill(HIST("debug/hEventCounter"), 9.); // passed RCT flag
		// }
		// return true;
	// }

  void init(InitContext&)
  {
    ft0Det.calculateChannelCenter();
    LOGF(info, "UpcTestFITBitMapping initialized. whichThr=%d", static_cast<int>(whichThr));
  }

  void process(UDCollisionsFull::iterator const& collision, aod::UDCollisionFITBits const& fitBits)
  {
    if (maxEvents >= 0 && collision.globalIndex() >= maxEvents) {
      return;
    }

    registry.fill(HIST("debug/hEventCounter"), 0.);
    registry.fill(HIST("debug/hFitBitsSize"), fitBits.size());
    registry.fill(HIST("debug/hCollisionIndexMod"), collision.globalIndex() % 100);

    if (fitBits.size() == 0) {
      if (collision.globalIndex() < debugPrintFirst ||
          (debugPrintEvery > 0 && collision.globalIndex() % debugPrintEvery == 0)) {
        LOGF(info, "collision %d: fitBits.size() = 0", collision.globalIndex());
      }
      registry.fill(HIST("debug/hEventCounter"), 1.);
      return;
			
    } else if (fitBits.size() == 1) {
      LOGF(debug, "collision %d: fitBits.size() = %d , as it should be", collision.globalIndex(), fitBits.size());
			registry.fill(HIST("debug/hEventCounter"), 2.);
			
    } else { // > 1 case
			if (collision.globalIndex() < debugPrintFirst ||
          (debugPrintEvery > 0 && collision.globalIndex() % debugPrintEvery == 0)) {
        LOGF(warn, "collision %d: fitBits.size() = %d (expected 0 or 1)", collision.globalIndex(), fitBits.size());
      }
			registry.fill(HIST("debug/hEventCounter"), 3.);
			return;
    }

		/* Checking collision level cuts */
		// if (!collisionPassesCuts(collision)) {
			// return;
		// }

		/* Only one row per collision is expected. */
    // auto const& row = *fitBits.begin();
    auto row = fitBits.begin();

    udhelpers::Bits256 w{};
    if (whichThr == 2) {
      w = udhelpers::makeBits256(row.thr2W0(), row.thr2W1(), row.thr2W2(), row.thr2W3());
    } else {
      w = udhelpers::makeBits256(row.thr1W0(), row.thr1W1(), row.thr1W2(), row.thr1W3());
    }

    const int nFT0A = countBitsInRange(w, 0, udhelpers::kFT0AChannels - 1);
    const int nFT0C = countBitsInRange(w, udhelpers::kFT0AChannels, udhelpers::kFT0Bits - 1);
    const int nFT0 = nFT0A + nFT0C;
    const int nFV0A = countBitsInRange(w, udhelpers::kFT0Bits, udhelpers::kTotalBits - 1);
    const int nFIT = nFT0 + nFV0A;

    registry.fill(HIST("mult/hPnFT0A"), nFT0A);
    registry.fill(HIST("mult/hPnFT0C"), nFT0C);
    registry.fill(HIST("mult/hPnFT0"), nFT0);
    registry.fill(HIST("mult/hPnFV0A"), nFV0A);
    registry.fill(HIST("mult/hPnFIT"), nFIT);
    registry.fill(HIST("mult/hNfiredA_vs_C"), nFT0A, nFT0C);

    if (collision.globalIndex() < debugPrintFirst ||
        (debugPrintEvery > 0 && collision.globalIndex() % debugPrintEvery == 0)) {
      LOGF(info,
           "collision %d: fitBits.size=%d  nFT0A=%d nFT0C=%d nFT0=%d nFV0A=%d nFIT=%d",
           collision.globalIndex(), fitBits.size(), nFT0A, nFT0C, nFT0, nFV0A, nFIT);
    }

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
        registry.fill(HIST("map/hPhiA"), phi);
        registry.fill(HIST("map/hEtaA"), eta);
        registry.fill(HIST("map/hEtaPhiA"), eta, phi);
      } else {
        registry.fill(HIST("map/hPhiC"), phi);
        registry.fill(HIST("map/hEtaC"), eta);
        registry.fill(HIST("map/hEtaPhiC"), eta, phi);
      }
    }

    registry.fill(HIST("debug/hEventCounter"), 4.);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcTestFITBitMapping>(cfgc, TaskName{"fitbit-mapping"})};
}