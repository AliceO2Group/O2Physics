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

using UDCollisionsFull = soa::Join<
  aod::UDCollisions,
  aod::UDCollisionsSels,
  aod::UDCollisionSelExtras,
	aod::UDZdcsReduced>;

// using UDCollisionsFull = aod::UDCollisions;

struct UpcTestFITBitMapping {
	
	SGSelector sgSelector;

  Configurable<int>  whichThr{"whichThr", 1, "Use 1=Thr1 bits or 2=Thr2 bits"};
  Configurable<int>  maxEvents{"maxEvents", -1, "Process at most this many collisions (-1 = all)"};
  Configurable<int>  debugPrintEvery{"debugPrintEvery", 1000, "Print every N events (-1 disables)"};
  Configurable<int>  debugPrintFirst{"debugPrintFirst", 10, "Always print first N events"};
  Configurable<bool> useRctFlag{"useRctFlag", true, "use RCT flags for event selection"};
  Configurable<int>  cutRctFlag{"cutRctFlag", 2, "0 = off, 1 = CBT, 2 = CBT+ZDC, 3 = CBThadron, 4 = CBThadron+ZDC"};
	
	Configurable<int> applyGapSideZDCVeto{"applyGapSideZDCVeto", true, "Apply ZDC neutron veto on the gap side"};
	Configurable<float> cutZNAnoNeutron{"cutZNAnoNeutron", 1.f, "Maximum allowed ZNA energy for zero-neutron condition"};
	Configurable<float> cutZNCnoNeutron{"cutZNCnoNeutron", 1.f, "Maximum allowed ZNC energy for zero-neutron condition"};
	
	Configurable<int> useFT0TimeCut{"useFT0TimeCut", true, "Apply direct FT0 timing cut"};
	Configurable<float> cutFT0TimeA{"cutFT0TimeA", 2., "FT0-A time cut [ns]"};
	Configurable<float> cutFT0TimeC{"cutFT0TimeC", 2., "FT0-C time cut [ns]"};
	
	Configurable<int> useFT0AVeto{"useFT0AVeto", true, "Apply FT0-A total amplitude veto"};
	Configurable<float> cutFT0AmpAVeto{"cutFT0AmpAVeto", 100., "Maximum allowed FT0-A total amplitude for veto"};
	
	Configurable<int> useFV0AVeto{"useFV0AVeto", true, "Apply FV0-A total amplitude veto"};
	Configurable<float> cutFV0AmpAVeto{"cutFV0AmpAVeto", 50., "Maximum allowed FV0-A total amplitude for veto"};
	
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
      {"debug/hEventCounter", "Event counter;step;events", {HistType::kTH1F, {{20, -0.5, 19.5}}}},
      {"debug/hFitBitsSize", "fitBits.size() per collision;fitBits.size();events", {HistType::kTH1F, {{10, -0.5, 9.5}}}},
      {"debug/hCollisionIndexMod", "Collision index mod 100;collision.globalIndex() % 100;events", {HistType::kTH1F, {{100, -0.5, 99.5}}}},
			{"debug/hFT0CChannelOccupancy", "FT0C fired channel occupancy;FT0C channel;counts",{HistType::kTH1F, {{112, 95.5, 207.5}}}},

      {"map/hPhiA", "FT0A #varphi;#varphi;counts", {HistType::kTH1F, {{9, 0, 2*M_PI}}}},
      {"map/hEtaA", "FT0A #eta;#eta;counts", {HistType::kTH1F, {{7, 3.4, 4.8}}}},
      {"map/hEtaPhiA", "FT0A #eta vs #varphi;#eta;#varphi", {HistType::kTH2F, {{7, 3.4, 4.8}, {9, 0, 2*M_PI}}}},

      {"map/hPhiC", "FT0C #varphi;#varphi;counts", {HistType::kTH1F, {{9, 0, 2*M_PI}}}},
      {"map/hEtaC", "FT0C #eta;#eta;counts", {HistType::kTH1F, {{5, -3.1, -2.1}}}},
      {"map/hEtaPhiC", "FT0C #eta vs #varphi;#eta;#varphi", {HistType::kTH2F, {{5, -3.1, -2.1}, {9, 0, 2*M_PI}}}},

      {"mult/hPnFT0A", "P(n): FT0A fired-channel multiplicity;N_{fired}^{FT0A};events", {HistType::kTH1F, {{97, -0.5, 96.5}}}},
      {"mult/hPnFT0C", "P(n): FT0C fired-channel multiplicity;N_{fired}^{FT0C};events", {HistType::kTH1F, {{50, -0.5, 49.5}}}},
      {"mult/hPnFT0", "P(n): FT0 fired-channel multiplicity;N_{fired}^{FT0};events", {HistType::kTH1F, {{209, -0.5, 208.5}}}},
      {"mult/hPnFV0A", "P(n): FV0A fired-channel multiplicity;N_{fired}^{FV0A};events", {HistType::kTH1F, {{49, -0.5, 48.5}}}},
      {"mult/hPnFIT", "P(n): FIT fired-channel multiplicity;N_{fired}^{FIT};events", {HistType::kTH1F, {{257, -0.5, 256.5}}}},
      // {"mult/hPnFT0A_sel", "Selected P(n): FT0A ;N_{fired}^{FT0A};events", {HistType::kTH1F, {{97, -0.5, 96.5}}}},
      // {"mult/hPnFT0C_sel", "Selected P(n): FT0C ;N_{fired}^{FT0C};events", {HistType::kTH1F, {{50, -0.5, 49.5}}}},
      // {"mult/hPnFT0_sel", "Selected P(n): FT0 ;N_{fired}^{FT0};events", {HistType::kTH1F, {{209, -0.5, 208.5}}}},
      // {"mult/hPnFV0A_sel", "Selected P(n): FV0A ;N_{fired}^{FV0A};events", {HistType::kTH1F, {{49, -0.5, 48.5}}}},
      // {"mult/hPnFIT_sel", "Selected P(n): FIT ;N_{fired}^{FIT};events", {HistType::kTH1F, {{257, -0.5, 256.5}}}},

      {"mult/hNfiredA_vs_C", "FT0A vs FT0C fired channels;N_{fired}^{FT0A};N_{fired}^{FT0C}", {HistType::kTH2F, {{97, -0.5, 96.5}, {113, -0.5, 112.5}}}},
		  // {"mult/hNfiredA_vs_C_sel", "Selected FT0A vs FT0C fired channels;N_{fired}^{FT0A};N_{fired}^{FT0C}",
       // {HistType::kTH2F, {{97, -0.5, 96.5}, {113, -0.5, 112.5}}}}
			 
			{"qaBeforeCuts/hFT0AAmplitudeVsTime","FT0A amplitude vs time before cuts;FT0A time (ns);FT0A total amplitude", {HistType::kTH2F, {{200, -40., 40.}, {300, 0., 3000.}}}},
			{"qaBeforeCuts/hFT0CAmplitudeVsTime","FT0C amplitude vs time before cuts;FT0C time (ns);FT0C total amplitude", {HistType::kTH2F, {{200, -40., 40.}, {300, 0., 3000.}}}},
			
			{"qaBeforeCuts/hZNEnergy","ZNA vs ZNC energy before cuts;ZNA energy;ZNC energy", {HistType::kTH2F, {{250, -5.0, 20.}, {250, -5.0, 20.}}}},
			
			{"map/hXYA", "FT0A fired channels in x-y;x [cm];y [cm]", {HistType::kTH2F, {{12, -18., 18.}, {12, -18., 18.}}}},
			{"map/hXYC", "FT0C fired channels in x-y;x [cm];y [cm]", {HistType::kTH2F, {{12, -18., 18.}, {12, -18., 18.}}}}
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
	
	template <typename C>
	bool isGoodRctFlag(const C& collision)
	{
		switch (cutRctFlag) {
			case 1: // CBT
				return sgSelector.isCBTOk(collision);
			case 2: // CBT + ZDC
				return sgSelector.isCBTZdcOk(collision);
			case 3: // CBT hadron
				return sgSelector.isCBTHadronOk(collision);
			case 4: // CBT hadron + ZDC
				return sgSelector.isCBTHadronZdcOk(collision);
			default: // no RCT cut applied
				return true;
		}
	}
	
	template <typename C>
	bool isZeroNeutronA(const C& collision)
	{
		const float eZNA = collision.energyCommonZNA(); // exact accessor to be checked
		return eZNA <= cutZNAnoNeutron;
	}

	template <typename C>
	bool isZeroNeutronC(const C& collision)
	{
		const float eZNC = collision.energyCommonZNC(); // exact accessor to be checked
		return eZNC <= cutZNCnoNeutron;
	}
	
	
	template <typename C>
	bool collisionPassesCuts(const C& collision)
	{
		/* good vertex */
		if (!collision.vtxITSTPC()) {
			return false;
		}
		registry.fill(HIST("debug/hEventCounter"), 4.);

		/* same bunch pile-up rejection */
		if (!collision.sbp()) {
			return false;
		}
		registry.fill(HIST("debug/hEventCounter"), 5.);

		/* ITS ROF rejection */
		if (!collision.itsROFb()) {
			return false;
		}
		registry.fill(HIST("debug/hEventCounter"), 6.);

		/* Timeframe border collision rejection */
		if (!collision.tfb()) {
			return false;
		}
		registry.fill(HIST("debug/hEventCounter"), 7.);

		/* z-vertex cut */
		if (std::abs(collision.posZ()) > 10.0) {
			return false;
		}
		registry.fill(HIST("debug/hEventCounter"), 8.);
		
		/* RCT flag check*/
		if (useRctFlag) {
			if (!isGoodRctFlag(collision)) {
				return false;
			}
			registry.fill(HIST("debug/hEventCounter"), 9.); // passed RCT flag
		}

		/* FT0-A veto on total amplitude */
		if (useFT0AVeto) {
			const float ampA = collision.totalFT0AmplitudeA();
			if (ampA > cutFT0AmpAVeto) {
				return false;
			}
			registry.fill(HIST("debug/hEventCounter"), 10.);
		}
		
		/* FV0-A veto on total amplitude */
		if (useFV0AVeto) {
			const float ampA = collision.totalFV0AmplitudeA();
			if (ampA > cutFV0AmpAVeto) {
				return false;
			}
			registry.fill(HIST("debug/hEventCounter"), 11.);
		}
		
		/* direct FT0 timing cuts */
		if (useFT0TimeCut) {
			const float tA = collision.timeFT0A();
			const float tC = collision.timeFT0C();

			if (abs(tA) > cutFT0TimeA) {
				return false;
			}
			if (abs(tC) > cutFT0TimeC) {
				return false;
			}
			registry.fill(HIST("debug/hEventCounter"), 12.);
		}
		
		
		int gapSide = 0 ; // Veto on A-side for now.
		if (applyGapSideZDCVeto) {
			if (gapSide == 0) { // gap-A
				if (!isZeroNeutronA(collision)) {
					return false;
				}
			} else if (gapSide == 1) { // gap-C
				if (!isZeroNeutronC(collision)) {
					return false;
				}
			}
			registry.fill(HIST("debug/hEventCounter"), 13.);
		}
		
		/*
		if (applyGapSideZDCVeto) {
			if (!isZeroNeutronA(collision) || !isZeroNeutronC(collision)) {
				return false;
			}
			registry.fill(HIST("debug/hEventCounter"), 13.);
		}
		*/
		return true;
	}

  void init(InitContext&)
	{
		ft0Det.calculateChannelCenter();

		LOGF(info, "UpcTestFITBitMapping initialized. whichThr=%d", static_cast<int>(whichThr));

		LOGF(info, "=== FT0A geometry dump ===");
		for (int ch = 0; ch < udhelpers::kFT0AChannels; ++ch) {
			auto pos = ft0Det.getChannelCenter(ch);

			const double x = pos.X();
			const double y = pos.Y();
			const double z = pos.Z();

			const double phi = RecoDecay::phi(x, y);
			const double r = std::sqrt(x * x + y * y);
			const double theta = std::atan2(r, z);
			const double eta = -std::log(std::tan(0.5 * theta));

			LOGF(info, "FT0A ch=%3d  x=%+8.4f y=%+8.4f z=%+8.4f  phi=%+8.4f eta=%+8.4f",
					 ch, x, y, z, phi, eta);
		}

		LOGF(info, "=== FT0C geometry dump ===");
		for (int ch = udhelpers::kFT0AChannels; ch < udhelpers::kFT0Bits; ++ch) {
			auto pos = ft0Det.getChannelCenter(ch);

			const double x = pos.X();
			const double y = pos.Y();
			const double z = pos.Z();

			const double phi = RecoDecay::phi(x, y);
			const double r = std::sqrt(x * x + y * y);
			const double theta = std::atan2(r, z);
			const double eta = -std::log(std::tan(0.5 * theta));

			LOGF(info, "FT0C ch=%3d  x=%+8.4f y=%+8.4f z=%+8.4f  phi=%+8.4f eta=%+8.4f",
					 ch, x, y, z, phi, eta);
		}

		LOGF(info, "=== end of FT0 geometry dump ===");
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
		
		const float ft0TimeA = collision.timeFT0A();
		const float ft0TimeC = collision.timeFT0C();
		const float ft0AmpA  = collision.totalFT0AmplitudeA();
		const float ft0AmpC  = collision.totalFT0AmplitudeC();

		registry.fill(HIST("qaBeforeCuts/hFT0AAmplitudeVsTime"), ft0TimeA, ft0AmpA);
		registry.fill(HIST("qaBeforeCuts/hFT0CAmplitudeVsTime"), ft0TimeC, ft0AmpC);
		registry.fill(HIST("qaBeforeCuts/hZNEnergy"), collision.energyCommonZNA(), collision.energyCommonZNC());

		/* Checking collision level cuts */
		if (!collisionPassesCuts(collision)) {
			return;
		}

		/* Only one row per collision is expected. */
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
		
		// /* A-side veto using detector-level multiplicity proxy */
		// if (applyFT0Aveto) {
			// if (nFT0A > maxFT0AforVeto) {
				// return;
			// }
			// registry.fill(HIST("debug/hEventCounter"), 13.);
		// }
		
		// /* Optional stronger veto with FV0A */
		// if (applyFV0Aveto) {
			// if (nFV0A > maxFV0AforVeto) {
				// return;
			// }
			// registry.fill(HIST("debug/hEventCounter"), 14.);
		// }
		
    // registry.fill(HIST("mult/hPnFT0A_sel"), nFT0A);
    // registry.fill(HIST("mult/hPnFT0C_sel"), nFT0C);
    // registry.fill(HIST("mult/hPnFT0_sel"), nFT0);
    // registry.fill(HIST("mult/hPnFV0A_sel"), nFV0A);
    // registry.fill(HIST("mult/hPnFIT_sel"), nFIT);
    // registry.fill(HIST("mult/hNfiredA_vs_C_sel"), nFT0A, nFT0C);

    if (collision.globalIndex() < debugPrintFirst ||
        (debugPrintEvery > 0 && collision.globalIndex() % debugPrintEvery == 0)) {
      LOGF(info,
           "collision %d: fitBits.size=%d  nFT0A=%d nFT0C=%d nFT0=%d nFV0A=%d nFIT=%d",
           collision.globalIndex(), fitBits.size(), nFT0A, nFT0C, nFT0, nFV0A, nFIT);
    }
		if (collision.globalIndex() < 5) {
			for (int bit = 96; bit < 208; ++bit) {
				if (udhelpers::testBit(w, bit)) {
					LOGF(info, "ANALYSIS sees fired bit %d", bit);
				}
			}
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

			auto pos = ft0Det.getChannelCenter(bit);

			const double x = pos.X() - offsetFT0[iRunOffset].getX();
			const double y = pos.Y() - offsetFT0[iRunOffset].getY();

      if (bit < udhelpers::kFT0AChannels) {
				// LOGF(info,"phi = %.4f, phiBin = %i",phi,phiBin);
				registry.fill(HIST("map/hPhiA"), phi);
        registry.fill(HIST("map/hEtaA"), eta);
        registry.fill(HIST("map/hEtaPhiA"), eta, phi);
				registry.fill(HIST("map/hXYA"), x, y);
      } else {
				registry.fill(HIST("debug/hFT0CChannelOccupancy"), bit);
        registry.fill(HIST("map/hPhiC"), phi);
        registry.fill(HIST("map/hEtaC"), eta);
        registry.fill(HIST("map/hEtaPhiC"), eta, phi);
			  registry.fill(HIST("map/hXYC"), x, y);
      }
    }

    registry.fill(HIST("debug/hEventCounter"), 4.);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<UpcTestFITBitMapping>(cfgc, TaskName{"fitbitmapping"})};
}