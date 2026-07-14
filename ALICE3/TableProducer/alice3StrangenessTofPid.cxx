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
///
/// \file alice3StrangenessTofPid.cxx
/// \brief This task produces tof strangeness pid tables
/// \author Jesper Karlsson Gumprecht
///

#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "ALICE3/DataModel/OTFStrangeness.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "ALICE3/Utils/a3StrangenessTofPidHelper.h"
#include "Common/Core/trackUtilities.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/TrackParametrization.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <string>

using namespace o2;
using namespace o2::framework;
using namespace o2::upgrade::stratofpid;

using CascTofResults = StrangenessTofResults<Topology::Cascade>;
using Alice3Tracks = soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels, aod::TracksDCA, aod::TracksExtraA3>;
using FullCascadeCandidates = soa::Join<aod::StoredCascCores, aod::CascIndices, aod::A3CascadeMcLabels>;
using Alice3Collision = soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator;

struct Alice3StrangenessTofPid {
  Produces<aod::A3XiInnerTofPid> tableA3XiInnerTofPid;                     // contains xi strangeness iTOF nsigma
  Produces<aod::A3XiOuterTofPid> tableA3XiOuterTofPid;                     // contains xi strangeness oTOF nsigma
  Produces<aod::A3OmegaInnerTofPid> tableA3OmegaInnerTofPid;               // contains omega strangeness iTOF nsigma
  Produces<aod::A3OmegaOuterTofPid> tableA3OmegaOuterTofPid;               // contains omega strangeness oTOF nsigma
  Produces<aod::A3XiExpectedInnerTimes> tableA3XiExpectedInnerTimes;       // contains xi strangeness iTOF info
  Produces<aod::A3XiExpectedOuterTimes> tableA3XiExpectedOuterTimes;       // contains xi strangeness oTOF info
  Produces<aod::A3OmegaExpectedInnerTimes> tableA3OmegaExpectedInnerTimes; // contains omega strangeness iTOF info
  Produces<aod::A3OmegaExpectedOuterTimes> tableA3OmegaExpectedOuterTimes; // contains omega strangeness oTOF info

  struct : ConfigurableGroup {
    ConfigurableAxis axisNSigma{"axisNSigma", {200, -10, 10}, "N sigma axis"};
    ConfigurableAxis axisCollisionTime{"axisCollisionTime", {200, -1e-3, 1e-3}, "Time delta axis"};
    ConfigurableAxis axisTimeDelta{"axisTimeDelta", {300, -1500, 1500}, "Time delta axis"};
    ConfigurableAxis axisTime{"axisTime", {1000, 0, 10000}, "Time axis"};
    ConfigurableAxis axisPt{"axisPt", {VARIABLE_WIDTH, 0.0f, 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f, 1.3f, 1.4f, 1.5f, 1.6f, 1.7f, 1.8f, 1.9f, 2.0f, 2.2f, 2.4f, 2.6f, 2.8f, 3.0f, 3.2f, 3.4f, 3.6f, 3.8f, 4.0f, 4.4f, 4.8f, 5.2f, 5.6f, 6.0f, 6.5f, 7.0f, 7.5f, 8.0f, 9.0f, 10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 17.0f, 19.0f, 21.0f, 23.0f, 25.0f, 30.0f, 35.0f, 40.0f, 50.0f}, "pt axis for QA histograms"};
  } axes;

  struct : ConfigurableGroup {
    std::string prefix = "cfgTof";

    // Time Of Flight
    Configurable<float> innerRadius{"innerRadius", 21., "Radius of iTOF (cm)"};
    Configurable<float> innerResolution{"innerResolution", 20., "Timing resolution of iTOF (ps)"};
    Configurable<float> outerRadius{"outerRadius", 92., "Radius of oTOF (cm)"};
    Configurable<float> outerResolution{"outerResolution", 20., "Timing resolution of oTOF (ps)"};
  } cfgTof;

  struct : ConfigurableGroup {
    std::string prefix = "produce";

    // Tables
    Configurable<bool> tableInnerNSigmaXi{"tableInnerNSigmaXi", true, "Produce NSigma tables for xi candidate"};
    Configurable<bool> tableOuterNSigmaXi{"tableOuterNSigmaXi", true, "Produce NSigma tables for xi candidate"};
    Configurable<bool> tableInnerExpectedTimesXi{"tableInnerExpectedTimesXi", false, "Produce expected and measured times tables for xi candidate"};
    Configurable<bool> tableOuterExpectedTimesXi{"tableOuterExpectedTimesXi", false, "Produce expected and measured times tables for xi candidate"};
    Configurable<bool> tableInnerNSigmaOmega{"tableInnerNSigmaOmega", true, "Produce NSigma tables for omega candidate"};
    Configurable<bool> tableOuterNSigmaOmega{"tableOuterNSigmaOmega", true, "Produce NSigma tables for omega candidate"};
    Configurable<bool> tableInnerExpectedTimesOmega{"tableInnerExpectedTimesOmega", false, "Produce expected and measured times tables for omega candidate"};
    Configurable<bool> tableOuterExpectedTimesOmega{"tableOuterExpectedTimesOmega", false, "Produce expected and measured times tables for omega candidate"};

    // Histograms
    Configurable<bool> histosXi{"histosXi", true, "Produce histograms for xi candidates"};
    Configurable<bool> histosAntiXi{"histosAntiXi", true, "Produce histograms for anti xi candidates"};
    Configurable<bool> histosOmega{"histosOmega", true, "Produce histograms for omega candidates"};
    Configurable<bool> histosAntiOmega{"histosAntiOmega", true, "Produce histograms for anti omega candidates"};

  } produce;

  Configurable<float> magneticField{"magneticField", 20.0f, "Magnetic field (in kilogauss)"};

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  StrangenessTofPid<Topology::Cascade> tofPidCasc;
  Service<o2::framework::O2DatabasePDG> pdgDB{};
  static constexpr float NanoToPico = 1e+3;
  static constexpr std::array<float, o2::track::kLabCovMatSize> ParentTrackCovMatrix{};

  void init(InitContext&)
  {
    // Todo: First try and get from tofpid task
    tofPidCasc.setResolution(cfgTof.innerResolution, cfgTof.outerResolution);
    tofPidCasc.setRadius(cfgTof.innerRadius, cfgTof.innerResolution);
    tofPidCasc.setMagneticField(magneticField);

    if (doprocessCascades) {
      if (produce.histosXi) {
        histos.add("Xi/hInnerArrivalTimeDeltaNeg", "hInnerArrivalTimeDeltaNeg", kTH1D, {{axes.axisTimeDelta}});
        histos.add("Xi/hInnerArrivalTimeDeltaPos", "hInnerArrivalTimeDeltaPos", kTH1D, {{axes.axisTimeDelta}});
        histos.add("Xi/hInnerArrivalTimeDeltaBach", "hInnerArrivalTimeDeltaBach", kTH1D, {{axes.axisTimeDelta}});
        histos.add("Xi/hOuterArrivalTimeDeltaNeg", "hOuterArrivalTimeDeltaNeg", kTH1D, {{axes.axisTimeDelta}});
        histos.add("Xi/hOuterArrivalTimeDeltaPos", "hOuterArrivalTimeDeltaPos", kTH1D, {{axes.axisTimeDelta}});
        histos.add("Xi/hOuterArrivalTimeDeltaBach", "hOuterArrivalTimeDeltaBach", kTH1D, {{axes.axisTimeDelta}});

        histos.add("Xi/hInnerNSigmaNeg", "hInnerNSigmaNeg", kTH1D, {{axes.axisNSigma}});
        histos.add("Xi/hInnerNSigmaPos", "hInnerNSigmaPos", kTH1D, {{axes.axisNSigma}});
        histos.add("Xi/hInnerNSigmaBach", "hInnerNSigmaBach", kTH1D, {{axes.axisNSigma}});
        histos.add("Xi/hOuterNSigmaNeg", "hOuterNSigmaNeg", kTH1D, {{axes.axisNSigma}});
        histos.add("Xi/hOuterNSigmaPos", "hOuterNSigmaPos", kTH1D, {{axes.axisNSigma}});
        histos.add("Xi/hOuterNSigmaBach", "hOuterNSigmaBach", kTH1D, {{axes.axisNSigma}});

        histos.add("Xi/hInnerExpectedArrivalTimeNeg", "hInnerExpectedArrivalTimeNeg", kTH1D, {{axes.axisTime}});
        histos.add("Xi/hInnerExpectedArrivalTimePos", "hInnerExpectedArrivalTimePos", kTH1D, {{axes.axisTime}});
        histos.add("Xi/hInnerExpectedArrivalTimeBach", "hInnerExpectedArrivalTimeBach", kTH1D, {{axes.axisTime}});
        histos.add("Xi/hOuterExpectedArrivalTimeNeg", "hOuterExpectedArrivalTimeNeg", kTH1D, {{axes.axisTime}});
        histos.add("Xi/hOuterExpectedArrivalTimePos", "hOuterExpectedArrivalTimePos", kTH1D, {{axes.axisTime}});
        histos.add("Xi/hOuterExpectedArrivalTimeBach", "hOuterExpectedArrivalTimeBach", kTH1D, {{axes.axisTime}});

        histos.add("Xi/hInnerMeasuredArrivalTimeNeg", "hInnerMeasuredArrivalTimeNeg", kTH1D, {{axes.axisTime}});
        histos.add("Xi/hInnerMeasuredArrivalTimePos", "hInnerMeasuredArrivalTimePos", kTH1D, {{axes.axisTime}});
        histos.add("Xi/hInnerMeasuredArrivalTimeBach", "hInnerMeasuredArrivalTimeBach", kTH1D, {{axes.axisTime}});
        histos.add("Xi/hOuterMeasuredArrivalTimeNeg", "hOuterMeasuredArrivalTimeNeg", kTH1D, {{axes.axisTime}});
        histos.add("Xi/hOuterMeasuredArrivalTimePos", "hOuterMeasuredArrivalTimePos", kTH1D, {{axes.axisTime}});
        histos.add("Xi/hOuterMeasuredArrivalTimeBach", "hOuterMeasuredArrivalTimeBach", kTH1D, {{axes.axisTime}});

        histos.add("Xi/hInnerMeasuredVsExpectedArrivalTimeNeg", "hInnerMeasuredVsExpectedArrivalTimeNeg;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("Xi/hInnerMeasuredVsExpectedArrivalTimePos", "hInnerMeasuredVsExpectedArrivalTimePos;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("Xi/hInnerMeasuredVsExpectedArrivalTimeBach", "hInnerMeasuredVsExpectedArrivalTimeBach;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("Xi/hOuterMeasuredVsExpectedArrivalTimeNeg", "hOuterMeasuredVsExpectedArrivalTimeNeg;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("Xi/hOuterMeasuredVsExpectedArrivalTimePos", "hOuterMeasuredVsExpectedArrivalTimePos;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("Xi/hOuterMeasuredVsExpectedArrivalTimeBach", "hOuterMeasuredVsExpectedArrivalTimeBach;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
      }

      if (produce.histosAntiXi) {
        histos.add("AntiXi/hInnerArrivalTimeDeltaNeg", "hInnerArrivalTimeDeltaNeg", kTH1D, {{axes.axisTimeDelta}});
        histos.add("AntiXi/hInnerArrivalTimeDeltaPos", "hInnerArrivalTimeDeltaPos", kTH1D, {{axes.axisTimeDelta}});
        histos.add("AntiXi/hInnerArrivalTimeDeltaBach", "hInnerArrivalTimeDeltaBach", kTH1D, {{axes.axisTimeDelta}});
        histos.add("AntiXi/hOuterArrivalTimeDeltaNeg", "hOuterArrivalTimeDeltaNeg", kTH1D, {{axes.axisTimeDelta}});
        histos.add("AntiXi/hOuterArrivalTimeDeltaPos", "hOuterArrivalTimeDeltaPos", kTH1D, {{axes.axisTimeDelta}});
        histos.add("AntiXi/hOuterArrivalTimeDeltaBach", "hOuterArrivalTimeDeltaBach", kTH1D, {{axes.axisTimeDelta}});

        histos.add("AntiXi/hInnerNSigmaNeg", "hInnerNSigmaNeg", kTH1D, {{axes.axisNSigma}});
        histos.add("AntiXi/hInnerNSigmaPos", "hInnerNSigmaPos", kTH1D, {{axes.axisNSigma}});
        histos.add("AntiXi/hInnerNSigmaBach", "hInnerNSigmaBach", kTH1D, {{axes.axisNSigma}});
        histos.add("AntiXi/hOuterNSigmaNeg", "hOuterNSigmaNeg", kTH1D, {{axes.axisNSigma}});
        histos.add("AntiXi/hOuterNSigmaPos", "hOuterNSigmaPos", kTH1D, {{axes.axisNSigma}});
        histos.add("AntiXi/hOuterNSigmaBach", "hOuterNSigmaBach", kTH1D, {{axes.axisNSigma}});

        histos.add("AntiXi/hInnerExpectedArrivalTimeNeg", "hInnerExpectedArrivalTimeNeg", kTH1D, {{axes.axisTime}});
        histos.add("AntiXi/hInnerExpectedArrivalTimePos", "hInnerExpectedArrivalTimePos", kTH1D, {{axes.axisTime}});
        histos.add("AntiXi/hInnerExpectedArrivalTimeBach", "hInnerExpectedArrivalTimeBach", kTH1D, {{axes.axisTime}});
        histos.add("AntiXi/hOuterExpectedArrivalTimeNeg", "hOuterExpectedArrivalTimeNeg", kTH1D, {{axes.axisTime}});
        histos.add("AntiXi/hOuterExpectedArrivalTimePos", "hOuterExpectedArrivalTimePos", kTH1D, {{axes.axisTime}});
        histos.add("AntiXi/hOuterExpectedArrivalTimeBach", "hOuterExpectedArrivalTimeBach", kTH1D, {{axes.axisTime}});

        histos.add("AntiXi/hInnerMeasuredArrivalTimeNeg", "hInnerMeasuredArrivalTimeNeg", kTH1D, {{axes.axisTime}});
        histos.add("AntiXi/hInnerMeasuredArrivalTimePos", "hInnerMeasuredArrivalTimePos", kTH1D, {{axes.axisTime}});
        histos.add("AntiXi/hInnerMeasuredArrivalTimeBach", "hInnerMeasuredArrivalTimeBach", kTH1D, {{axes.axisTime}});
        histos.add("AntiXi/hOuterMeasuredArrivalTimeNeg", "hOuterMeasuredArrivalTimeNeg", kTH1D, {{axes.axisTime}});
        histos.add("AntiXi/hOuterMeasuredArrivalTimePos", "hOuterMeasuredArrivalTimePos", kTH1D, {{axes.axisTime}});
        histos.add("AntiXi/hOuterMeasuredArrivalTimeBach", "hOuterMeasuredArrivalTimeBach", kTH1D, {{axes.axisTime}});

        histos.add("AntiXi/hInnerMeasuredVsExpectedArrivalTimeNeg", "hInnerMeasuredVsExpectedArrivalTimeNeg;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("AntiXi/hInnerMeasuredVsExpectedArrivalTimePos", "hInnerMeasuredVsExpectedArrivalTimePos;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("AntiXi/hInnerMeasuredVsExpectedArrivalTimeBach", "hInnerMeasuredVsExpectedArrivalTimeBach;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("AntiXi/hOuterMeasuredVsExpectedArrivalTimeNeg", "hOuterMeasuredVsExpectedArrivalTimeNeg;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("AntiXi/hOuterMeasuredVsExpectedArrivalTimePos", "hOuterMeasuredVsExpectedArrivalTimePos;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("AntiXi/hOuterMeasuredVsExpectedArrivalTimeBach", "hOuterMeasuredVsExpectedArrivalTimeBach;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
      }

      if (produce.histosOmega) {
        histos.add("Omega/hInnerArrivalTimeDeltaNeg", "hInnerArrivalTimeDeltaNeg", kTH1D, {{axes.axisTimeDelta}});
        histos.add("Omega/hInnerArrivalTimeDeltaPos", "hInnerArrivalTimeDeltaPos", kTH1D, {{axes.axisTimeDelta}});
        histos.add("Omega/hInnerArrivalTimeDeltaBach", "hInnerArrivalTimeDeltaBach", kTH1D, {{axes.axisTimeDelta}});
        histos.add("Omega/hOuterArrivalTimeDeltaNeg", "hOuterArrivalTimeDeltaNeg", kTH1D, {{axes.axisTimeDelta}});
        histos.add("Omega/hOuterArrivalTimeDeltaPos", "hOuterArrivalTimeDeltaPos", kTH1D, {{axes.axisTimeDelta}});
        histos.add("Omega/hOuterArrivalTimeDeltaBach", "hOuterArrivalTimeDeltaBach", kTH1D, {{axes.axisTimeDelta}});

        histos.add("Omega/hInnerNSigmaNeg", "hInnerNSigmaNeg", kTH1D, {{axes.axisNSigma}});
        histos.add("Omega/hInnerNSigmaPos", "hInnerNSigmaPos", kTH1D, {{axes.axisNSigma}});
        histos.add("Omega/hInnerNSigmaBach", "hInnerNSigmaBach", kTH1D, {{axes.axisNSigma}});
        histos.add("Omega/hOuterNSigmaNeg", "hOuterNSigmaNeg", kTH1D, {{axes.axisNSigma}});
        histos.add("Omega/hOuterNSigmaPos", "hOuterNSigmaPos", kTH1D, {{axes.axisNSigma}});
        histos.add("Omega/hOuterNSigmaBach", "hOuterNSigmaBach", kTH1D, {{axes.axisNSigma}});

        histos.add("Omega/hInnerExpectedArrivalTimeNeg", "hInnerExpectedArrivalTimeNeg", kTH1D, {{axes.axisTime}});
        histos.add("Omega/hInnerExpectedArrivalTimePos", "hInnerExpectedArrivalTimePos", kTH1D, {{axes.axisTime}});
        histos.add("Omega/hInnerExpectedArrivalTimeBach", "hInnerExpectedArrivalTimeBach", kTH1D, {{axes.axisTime}});
        histos.add("Omega/hOuterExpectedArrivalTimeNeg", "hOuterExpectedArrivalTimeNeg", kTH1D, {{axes.axisTime}});
        histos.add("Omega/hOuterExpectedArrivalTimePos", "hOuterExpectedArrivalTimePos", kTH1D, {{axes.axisTime}});
        histos.add("Omega/hOuterExpectedArrivalTimeBach", "hOuterExpectedArrivalTimeBach", kTH1D, {{axes.axisTime}});

        histos.add("Omega/hInnerMeasuredArrivalTimeNeg", "hInnerMeasuredArrivalTimeNeg", kTH1D, {{axes.axisTime}});
        histos.add("Omega/hInnerMeasuredArrivalTimePos", "hInnerMeasuredArrivalTimePos", kTH1D, {{axes.axisTime}});
        histos.add("Omega/hInnerMeasuredArrivalTimeBach", "hInnerMeasuredArrivalTimeBach", kTH1D, {{axes.axisTime}});
        histos.add("Omega/hOuterMeasuredArrivalTimeNeg", "hOuterMeasuredArrivalTimeNeg", kTH1D, {{axes.axisTime}});
        histos.add("Omega/hOuterMeasuredArrivalTimePos", "hOuterMeasuredArrivalTimePos", kTH1D, {{axes.axisTime}});
        histos.add("Omega/hOuterMeasuredArrivalTimeBach", "hOuterMeasuredArrivalTimeBach", kTH1D, {{axes.axisTime}});

        histos.add("Omega/hInnerMeasuredVsExpectedArrivalTimeNeg", "hInnerMeasuredVsExpectedArrivalTimeNeg;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("Omega/hInnerMeasuredVsExpectedArrivalTimePos", "hInnerMeasuredVsExpectedArrivalTimePos;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("Omega/hInnerMeasuredVsExpectedArrivalTimeBach", "hInnerMeasuredVsExpectedArrivalTimeBach;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("Omega/hOuterMeasuredVsExpectedArrivalTimeNeg", "hOuterMeasuredVsExpectedArrivalTimeNeg;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("Omega/hOuterMeasuredVsExpectedArrivalTimePos", "hOuterMeasuredVsExpectedArrivalTimePos;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("Omega/hOuterMeasuredVsExpectedArrivalTimeBach", "hOuterMeasuredVsExpectedArrivalTimeBach;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
      }

      if (produce.histosAntiOmega) {
        histos.add("AntiOmega/hInnerArrivalTimeDeltaNeg", "hInnerArrivalTimeDeltaNeg", kTH1D, {{axes.axisTimeDelta}});
        histos.add("AntiOmega/hInnerArrivalTimeDeltaPos", "hInnerArrivalTimeDeltaPos", kTH1D, {{axes.axisTimeDelta}});
        histos.add("AntiOmega/hInnerArrivalTimeDeltaBach", "hInnerArrivalTimeDeltaBach", kTH1D, {{axes.axisTimeDelta}});
        histos.add("AntiOmega/hOuterArrivalTimeDeltaNeg", "hOuterArrivalTimeDeltaNeg", kTH1D, {{axes.axisTimeDelta}});
        histos.add("AntiOmega/hOuterArrivalTimeDeltaPos", "hOuterArrivalTimeDeltaPos", kTH1D, {{axes.axisTimeDelta}});
        histos.add("AntiOmega/hOuterArrivalTimeDeltaBach", "hOuterArrivalTimeDeltaBach", kTH1D, {{axes.axisTimeDelta}});

        histos.add("AntiOmega/hInnerNSigmaNeg", "hInnerNSigmaNeg", kTH1D, {{axes.axisNSigma}});
        histos.add("AntiOmega/hInnerNSigmaPos", "hInnerNSigmaPos", kTH1D, {{axes.axisNSigma}});
        histos.add("AntiOmega/hInnerNSigmaBach", "hInnerNSigmaBach", kTH1D, {{axes.axisNSigma}});
        histos.add("AntiOmega/hOuterNSigmaNeg", "hOuterNSigmaNeg", kTH1D, {{axes.axisNSigma}});
        histos.add("AntiOmega/hOuterNSigmaPos", "hOuterNSigmaPos", kTH1D, {{axes.axisNSigma}});
        histos.add("AntiOmega/hOuterNSigmaBach", "hOuterNSigmaBach", kTH1D, {{axes.axisNSigma}});

        histos.add("AntiOmega/hInnerExpectedArrivalTimeNeg", "hInnerExpectedArrivalTimeNeg", kTH1D, {{axes.axisTime}});
        histos.add("AntiOmega/hInnerExpectedArrivalTimePos", "hInnerExpectedArrivalTimePos", kTH1D, {{axes.axisTime}});
        histos.add("AntiOmega/hInnerExpectedArrivalTimeBach", "hInnerExpectedArrivalTimeBach", kTH1D, {{axes.axisTime}});
        histos.add("AntiOmega/hOuterExpectedArrivalTimeNeg", "hOuterExpectedArrivalTimeNeg", kTH1D, {{axes.axisTime}});
        histos.add("AntiOmega/hOuterExpectedArrivalTimePos", "hOuterExpectedArrivalTimePos", kTH1D, {{axes.axisTime}});
        histos.add("AntiOmega/hOuterExpectedArrivalTimeBach", "hOuterExpectedArrivalTimeBach", kTH1D, {{axes.axisTime}});

        histos.add("AntiOmega/hInnerMeasuredArrivalTimeNeg", "hInnerMeasuredArrivalTimeNeg", kTH1D, {{axes.axisTime}});
        histos.add("AntiOmega/hInnerMeasuredArrivalTimePos", "hInnerMeasuredArrivalTimePos", kTH1D, {{axes.axisTime}});
        histos.add("AntiOmega/hInnerMeasuredArrivalTimeBach", "hInnerMeasuredArrivalTimeBach", kTH1D, {{axes.axisTime}});
        histos.add("AntiOmega/hOuterMeasuredArrivalTimeNeg", "hOuterMeasuredArrivalTimeNeg", kTH1D, {{axes.axisTime}});
        histos.add("AntiOmega/hOuterMeasuredArrivalTimePos", "hOuterMeasuredArrivalTimePos", kTH1D, {{axes.axisTime}});
        histos.add("AntiOmega/hOuterMeasuredArrivalTimeBach", "hOuterMeasuredArrivalTimeBach", kTH1D, {{axes.axisTime}});

        histos.add("AntiOmega/hInnerMeasuredVsExpectedArrivalTimeNeg", "hInnerMeasuredVsExpectedArrivalTimeNeg;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("AntiOmega/hInnerMeasuredVsExpectedArrivalTimePos", "hInnerMeasuredVsExpectedArrivalTimePos;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("AntiOmega/hInnerMeasuredVsExpectedArrivalTimeBach", "hInnerMeasuredVsExpectedArrivalTimeBach;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("AntiOmega/hOuterMeasuredVsExpectedArrivalTimeNeg", "hOuterMeasuredVsExpectedArrivalTimeNeg;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("AntiOmega/hOuterMeasuredVsExpectedArrivalTimePos", "hOuterMeasuredVsExpectedArrivalTimePos;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
        histos.add("AntiOmega/hOuterMeasuredVsExpectedArrivalTimeBach", "hOuterMeasuredVsExpectedArrivalTimeBach;Measured time (ps);Expected time (ps)", kTH2D, {{axes.axisTime}, {axes.axisTime}});
      }
    }
    histos.print();
  }

  void processCascades(Alice3Collision const& collision, FullCascadeCandidates const& cascadeCandidates, Alice3Tracks const&, aod::McParticles const&, aod::McCollisions const&)
  {
    const std::array<float, 3> vtx = {collision.posX(), collision.posY(), collision.posZ()};
    const auto mcCollision = collision.template mcCollision_as<aod::McCollisions>();
    for (const auto& casc : cascadeCandidates) {
      const std::array<float, 3> cascSV = {casc.x(), casc.y(), casc.z()};
      const std::array<float, 3> cascP = {casc.px(), casc.py(), casc.pz()};
      const int chargeCascade = (casc.sign() > 0) ? 1 : -1;
      o2::track::TrackParCov cascade(cascSV, cascP, ParentTrackCovMatrix, chargeCascade);

      const std::array<float, 3> v0SV = {casc.xlambda(), casc.ylambda(), casc.zlambda()};
      const std::array<float, 3> v0P = {casc.pxpos() + casc.pxneg(), casc.pypos() + casc.pyneg(), casc.pzpos() + casc.pzneg()};
      static constexpr int ChargeV0 = 0;
      o2::track::TrackParCov v0(v0SV, v0P, ParentTrackCovMatrix, ChargeV0);

      const auto bachTrack = casc.template bachelor_as<Alice3Tracks>();
      const auto bachParticle = bachTrack.template mcParticle_as<aod::McParticles>();
      o2::track::TrackParCov bachelor = getTrackParCov(bachTrack);

      const auto posTrack = casc.template posTrack_as<Alice3Tracks>();
      const auto posParticle = posTrack.template mcParticle_as<aod::McParticles>();
      o2::track::TrackParCov positive = getTrackParCov(posTrack);

      const auto negTrack = casc.template negTrack_as<Alice3Tracks>();
      const auto negParticle = negTrack.template mcParticle_as<aod::McParticles>();
      o2::track::TrackParCov negative = getTrackParCov(negTrack);

      tofPidCasc.reset();
      tofPidCasc.setTrack(TrackType::Positive, positive, mcCollision.t() + posParticle.vt(), pdgDB->Mass(posParticle.pdgCode()));
      tofPidCasc.setTrack(TrackType::Negative, negative, mcCollision.t() + negParticle.vt(), pdgDB->Mass(negParticle.pdgCode()));
      tofPidCasc.setTrack(TrackType::Bachelor, bachelor, mcCollision.t() + bachParticle.vt(), pdgDB->Mass(bachParticle.pdgCode()));
      tofPidCasc.setTrack(TrackType::V0, v0);           // no track time
      tofPidCasc.setTrack(TrackType::Cascade, cascade); // no track time

      CascTofResults straTofResultsXi{}, straTofResultsAntiXi{}, straTofResultsOmega{}, straTofResultsAntiOmega{};
      straTofResultsXi = tofPidCasc.findNSigmas<CandidateType::Xi>(vtx, cascSV, v0SV);
      straTofResultsAntiXi = tofPidCasc.findNSigmas<CandidateType::AntiXi>(vtx, cascSV, v0SV);
      straTofResultsOmega = tofPidCasc.findNSigmas<CandidateType::Omega>(vtx, cascSV, v0SV);
      straTofResultsAntiOmega = tofPidCasc.findNSigmas<CandidateType::AntiOmega>(vtx, cascSV, v0SV);

      if (produce.histosXi) {
        if (straTofResultsXi.pos.hasInnerTof) {
          histos.fill(HIST("Xi/hInnerArrivalTimeDeltaPos"), (straTofResultsXi.pos.measuredTimeInner - straTofResultsXi.pos.expectedTimeInner) * NanoToPico);
          histos.fill(HIST("Xi/hInnerNSigmaPos"), straTofResultsXi.pos.nSigmaInner);
          histos.fill(HIST("Xi/hInnerExpectedArrivalTimePos"), straTofResultsXi.pos.expectedTimeInner * NanoToPico);
          histos.fill(HIST("Xi/hInnerMeasuredArrivalTimePos"), straTofResultsXi.pos.measuredTimeInner * NanoToPico);
          histos.fill(HIST("Xi/hInnerMeasuredVsExpectedArrivalTimePos"), straTofResultsXi.pos.measuredTimeInner * NanoToPico, straTofResultsXi.pos.expectedTimeInner * NanoToPico);
        }
        if (straTofResultsXi.neg.hasInnerTof) {
          histos.fill(HIST("Xi/hInnerArrivalTimeDeltaNeg"), (straTofResultsXi.neg.measuredTimeInner - straTofResultsXi.neg.expectedTimeInner) * NanoToPico);
          histos.fill(HIST("Xi/hInnerNSigmaNeg"), straTofResultsXi.neg.nSigmaInner);
          histos.fill(HIST("Xi/hInnerExpectedArrivalTimeNeg"), straTofResultsXi.neg.expectedTimeInner * NanoToPico);
          histos.fill(HIST("Xi/hInnerMeasuredArrivalTimeNeg"), straTofResultsXi.neg.measuredTimeInner * NanoToPico);
          histos.fill(HIST("Xi/hInnerMeasuredVsExpectedArrivalTimeNeg"), straTofResultsXi.neg.measuredTimeInner * NanoToPico, straTofResultsXi.neg.expectedTimeInner * NanoToPico);
        }
        if (straTofResultsXi.bach.hasInnerTof) {
          histos.fill(HIST("Xi/hInnerArrivalTimeDeltaBach"), (straTofResultsXi.bach.measuredTimeInner - straTofResultsXi.bach.expectedTimeInner) * NanoToPico);
          histos.fill(HIST("Xi/hInnerNSigmaBach"), straTofResultsXi.bach.nSigmaInner);
          histos.fill(HIST("Xi/hInnerExpectedArrivalTimeBach"), straTofResultsXi.bach.expectedTimeInner * NanoToPico);
          histos.fill(HIST("Xi/hInnerMeasuredArrivalTimeBach"), straTofResultsXi.bach.measuredTimeInner * NanoToPico);
          histos.fill(HIST("Xi/hInnerMeasuredVsExpectedArrivalTimeBach"), straTofResultsXi.bach.measuredTimeInner * NanoToPico, straTofResultsXi.bach.expectedTimeInner * NanoToPico);
        }
        if (straTofResultsXi.pos.hasOuterTof) {
          histos.fill(HIST("Xi/hOuterArrivalTimeDeltaPos"), (straTofResultsXi.pos.measuredTimeOuter - straTofResultsXi.pos.expectedTimeOuter) * NanoToPico);
          histos.fill(HIST("Xi/hOuterNSigmaPos"), straTofResultsXi.pos.nSigmaOuter);
          histos.fill(HIST("Xi/hOuterExpectedArrivalTimePos"), straTofResultsXi.pos.expectedTimeOuter * NanoToPico);
          histos.fill(HIST("Xi/hOuterMeasuredArrivalTimePos"), straTofResultsXi.pos.measuredTimeOuter * NanoToPico);
          histos.fill(HIST("Xi/hOuterMeasuredVsExpectedArrivalTimePos"), straTofResultsXi.pos.measuredTimeOuter * NanoToPico, straTofResultsXi.pos.expectedTimeOuter * NanoToPico);
        }
        if (straTofResultsXi.neg.hasOuterTof) {
          histos.fill(HIST("Xi/hOuterArrivalTimeDeltaNeg"), (straTofResultsXi.neg.measuredTimeOuter - straTofResultsXi.neg.expectedTimeOuter) * NanoToPico);
          histos.fill(HIST("Xi/hOuterNSigmaNeg"), straTofResultsXi.neg.nSigmaOuter);
          histos.fill(HIST("Xi/hOuterExpectedArrivalTimeNeg"), straTofResultsXi.neg.expectedTimeOuter * NanoToPico);
          histos.fill(HIST("Xi/hOuterMeasuredArrivalTimeNeg"), straTofResultsXi.neg.measuredTimeOuter * NanoToPico);
          histos.fill(HIST("Xi/hOuterMeasuredVsExpectedArrivalTimeNeg"), straTofResultsXi.neg.measuredTimeOuter * NanoToPico, straTofResultsXi.neg.expectedTimeOuter * NanoToPico);
        }
        if (straTofResultsXi.bach.hasOuterTof) {
          histos.fill(HIST("Xi/hOuterArrivalTimeDeltaBach"), (straTofResultsXi.bach.measuredTimeOuter - straTofResultsXi.bach.expectedTimeOuter) * NanoToPico);
          histos.fill(HIST("Xi/hOuterNSigmaBach"), straTofResultsXi.bach.nSigmaOuter);
          histos.fill(HIST("Xi/hOuterExpectedArrivalTimeBach"), straTofResultsXi.bach.expectedTimeOuter * NanoToPico);
          histos.fill(HIST("Xi/hOuterMeasuredArrivalTimeBach"), straTofResultsXi.bach.measuredTimeOuter * NanoToPico);
          histos.fill(HIST("Xi/hOuterMeasuredVsExpectedArrivalTimeBach"), straTofResultsXi.bach.measuredTimeOuter * NanoToPico, straTofResultsXi.bach.expectedTimeOuter * NanoToPico);
        }
      }
      if (produce.histosAntiXi) {
        if (straTofResultsAntiXi.pos.hasInnerTof) {
          histos.fill(HIST("AntiXi/hInnerArrivalTimeDeltaPos"), (straTofResultsAntiXi.pos.measuredTimeInner - straTofResultsAntiXi.pos.expectedTimeInner) * NanoToPico);
          histos.fill(HIST("AntiXi/hInnerNSigmaPos"), straTofResultsAntiXi.pos.nSigmaInner);
          histos.fill(HIST("AntiXi/hInnerExpectedArrivalTimePos"), straTofResultsAntiXi.pos.expectedTimeInner * NanoToPico);
          histos.fill(HIST("AntiXi/hInnerMeasuredArrivalTimePos"), straTofResultsAntiXi.pos.measuredTimeInner * NanoToPico);
          histos.fill(HIST("AntiXi/hInnerMeasuredVsExpectedArrivalTimePos"), straTofResultsAntiXi.pos.measuredTimeInner * NanoToPico, straTofResultsAntiXi.pos.expectedTimeInner * NanoToPico);
        }
        if (straTofResultsAntiXi.neg.hasInnerTof) {
          histos.fill(HIST("AntiXi/hInnerArrivalTimeDeltaNeg"), (straTofResultsAntiXi.neg.measuredTimeInner - straTofResultsAntiXi.neg.expectedTimeInner) * NanoToPico);
          histos.fill(HIST("AntiXi/hInnerNSigmaNeg"), straTofResultsAntiXi.neg.nSigmaInner);
          histos.fill(HIST("AntiXi/hInnerExpectedArrivalTimeNeg"), straTofResultsAntiXi.neg.expectedTimeInner * NanoToPico);
          histos.fill(HIST("AntiXi/hInnerMeasuredArrivalTimeNeg"), straTofResultsAntiXi.neg.measuredTimeInner * NanoToPico);
          histos.fill(HIST("AntiXi/hInnerMeasuredVsExpectedArrivalTimeNeg"), straTofResultsAntiXi.neg.measuredTimeInner * NanoToPico, straTofResultsAntiXi.neg.expectedTimeInner * NanoToPico);
        }
        if (straTofResultsAntiXi.bach.hasInnerTof) {
          histos.fill(HIST("AntiXi/hInnerArrivalTimeDeltaBach"), (straTofResultsAntiXi.bach.measuredTimeInner - straTofResultsAntiXi.bach.expectedTimeInner) * NanoToPico);
          histos.fill(HIST("AntiXi/hInnerNSigmaBach"), straTofResultsAntiXi.bach.nSigmaInner);
          histos.fill(HIST("AntiXi/hInnerExpectedArrivalTimeBach"), straTofResultsAntiXi.bach.expectedTimeInner * NanoToPico);
          histos.fill(HIST("AntiXi/hInnerMeasuredArrivalTimeBach"), straTofResultsAntiXi.bach.measuredTimeInner * NanoToPico);
          histos.fill(HIST("AntiXi/hInnerMeasuredVsExpectedArrivalTimeBach"), straTofResultsAntiXi.bach.measuredTimeInner * NanoToPico, straTofResultsAntiXi.bach.expectedTimeInner * NanoToPico);
        }
        if (straTofResultsAntiXi.pos.hasOuterTof) {
          histos.fill(HIST("AntiXi/hOuterArrivalTimeDeltaPos"), (straTofResultsAntiXi.pos.measuredTimeOuter - straTofResultsAntiXi.pos.expectedTimeOuter) * NanoToPico);
          histos.fill(HIST("AntiXi/hOuterNSigmaPos"), straTofResultsAntiXi.pos.nSigmaOuter);
          histos.fill(HIST("AntiXi/hOuterExpectedArrivalTimePos"), straTofResultsAntiXi.pos.expectedTimeOuter * NanoToPico);
          histos.fill(HIST("AntiXi/hOuterMeasuredArrivalTimePos"), straTofResultsAntiXi.pos.measuredTimeOuter * NanoToPico);
          histos.fill(HIST("AntiXi/hOuterMeasuredVsExpectedArrivalTimePos"), straTofResultsAntiXi.pos.measuredTimeOuter * NanoToPico, straTofResultsAntiXi.pos.expectedTimeOuter * NanoToPico);
        }
        if (straTofResultsAntiXi.neg.hasOuterTof) {
          histos.fill(HIST("AntiXi/hOuterArrivalTimeDeltaNeg"), (straTofResultsAntiXi.neg.measuredTimeOuter - straTofResultsAntiXi.neg.expectedTimeOuter) * NanoToPico);
          histos.fill(HIST("AntiXi/hOuterNSigmaNeg"), straTofResultsAntiXi.neg.nSigmaOuter);
          histos.fill(HIST("AntiXi/hOuterExpectedArrivalTimeNeg"), straTofResultsAntiXi.neg.expectedTimeOuter * NanoToPico);
          histos.fill(HIST("AntiXi/hOuterMeasuredArrivalTimeNeg"), straTofResultsAntiXi.neg.measuredTimeOuter * NanoToPico);
          histos.fill(HIST("AntiXi/hOuterMeasuredVsExpectedArrivalTimeNeg"), straTofResultsAntiXi.neg.measuredTimeOuter * NanoToPico, straTofResultsAntiXi.neg.expectedTimeOuter * NanoToPico);
        }
        if (straTofResultsAntiXi.bach.hasOuterTof) {
          histos.fill(HIST("AntiXi/hOuterArrivalTimeDeltaBach"), (straTofResultsAntiXi.bach.measuredTimeOuter - straTofResultsAntiXi.bach.expectedTimeOuter) * NanoToPico);
          histos.fill(HIST("AntiXi/hOuterNSigmaBach"), straTofResultsAntiXi.bach.nSigmaOuter);
          histos.fill(HIST("AntiXi/hOuterExpectedArrivalTimeBach"), straTofResultsAntiXi.bach.expectedTimeOuter * NanoToPico);
          histos.fill(HIST("AntiXi/hOuterMeasuredArrivalTimeBach"), straTofResultsAntiXi.bach.measuredTimeOuter * NanoToPico);
          histos.fill(HIST("AntiXi/hOuterMeasuredVsExpectedArrivalTimeBach"), straTofResultsAntiXi.bach.measuredTimeOuter * NanoToPico, straTofResultsAntiXi.bach.expectedTimeOuter * NanoToPico);
        }
      }
      if (produce.histosOmega) {
        if (straTofResultsOmega.pos.hasInnerTof) {
          histos.fill(HIST("Omega/hInnerArrivalTimeDeltaPos"), (straTofResultsOmega.pos.measuredTimeInner - straTofResultsOmega.pos.expectedTimeInner) * NanoToPico);
          histos.fill(HIST("Omega/hInnerNSigmaPos"), straTofResultsOmega.pos.nSigmaInner);
          histos.fill(HIST("Omega/hInnerExpectedArrivalTimePos"), straTofResultsOmega.pos.expectedTimeInner * NanoToPico);
          histos.fill(HIST("Omega/hInnerMeasuredArrivalTimePos"), straTofResultsOmega.pos.measuredTimeInner * NanoToPico);
          histos.fill(HIST("Omega/hInnerMeasuredVsExpectedArrivalTimePos"), straTofResultsOmega.pos.measuredTimeInner * NanoToPico, straTofResultsOmega.pos.expectedTimeInner * NanoToPico);
        }
        if (straTofResultsOmega.neg.hasInnerTof) {
          histos.fill(HIST("Omega/hInnerArrivalTimeDeltaNeg"), (straTofResultsOmega.neg.measuredTimeInner - straTofResultsOmega.neg.expectedTimeInner) * NanoToPico);
          histos.fill(HIST("Omega/hInnerNSigmaNeg"), straTofResultsOmega.neg.nSigmaInner);
          histos.fill(HIST("Omega/hInnerExpectedArrivalTimeNeg"), straTofResultsOmega.neg.expectedTimeInner * NanoToPico);
          histos.fill(HIST("Omega/hInnerMeasuredArrivalTimeNeg"), straTofResultsOmega.neg.measuredTimeInner * NanoToPico);
          histos.fill(HIST("Omega/hInnerMeasuredVsExpectedArrivalTimeNeg"), straTofResultsOmega.neg.measuredTimeInner * NanoToPico, straTofResultsOmega.neg.expectedTimeInner * NanoToPico);
        }
        if (straTofResultsOmega.bach.hasInnerTof) {
          histos.fill(HIST("Omega/hInnerArrivalTimeDeltaBach"), (straTofResultsOmega.bach.measuredTimeInner - straTofResultsOmega.bach.expectedTimeInner) * NanoToPico);
          histos.fill(HIST("Omega/hInnerNSigmaBach"), straTofResultsOmega.bach.nSigmaInner);
          histos.fill(HIST("Omega/hInnerExpectedArrivalTimeBach"), straTofResultsOmega.bach.expectedTimeInner * NanoToPico);
          histos.fill(HIST("Omega/hInnerMeasuredArrivalTimeBach"), straTofResultsOmega.bach.measuredTimeInner * NanoToPico);
          histos.fill(HIST("Omega/hInnerMeasuredVsExpectedArrivalTimeBach"), straTofResultsOmega.bach.measuredTimeInner * NanoToPico, straTofResultsOmega.bach.expectedTimeInner * NanoToPico);
        }
        if (straTofResultsOmega.pos.hasOuterTof) {
          histos.fill(HIST("Omega/hOuterArrivalTimeDeltaPos"), (straTofResultsOmega.pos.measuredTimeOuter - straTofResultsOmega.pos.expectedTimeOuter) * NanoToPico);
          histos.fill(HIST("Omega/hOuterNSigmaPos"), straTofResultsOmega.pos.nSigmaOuter);
          histos.fill(HIST("Omega/hOuterExpectedArrivalTimePos"), straTofResultsOmega.pos.expectedTimeOuter * NanoToPico);
          histos.fill(HIST("Omega/hOuterMeasuredArrivalTimePos"), straTofResultsOmega.pos.measuredTimeOuter * NanoToPico);
          histos.fill(HIST("Omega/hOuterMeasuredVsExpectedArrivalTimePos"), straTofResultsOmega.pos.measuredTimeOuter * NanoToPico, straTofResultsOmega.pos.expectedTimeOuter * NanoToPico);
        }
        if (straTofResultsOmega.neg.hasOuterTof) {
          histos.fill(HIST("Omega/hOuterArrivalTimeDeltaNeg"), (straTofResultsOmega.neg.measuredTimeOuter - straTofResultsOmega.neg.expectedTimeOuter) * NanoToPico);
          histos.fill(HIST("Omega/hOuterNSigmaNeg"), straTofResultsOmega.neg.nSigmaOuter);
          histos.fill(HIST("Omega/hOuterExpectedArrivalTimeNeg"), straTofResultsOmega.neg.expectedTimeOuter * NanoToPico);
          histos.fill(HIST("Omega/hOuterMeasuredArrivalTimeNeg"), straTofResultsOmega.neg.measuredTimeOuter * NanoToPico);
          histos.fill(HIST("Omega/hOuterMeasuredVsExpectedArrivalTimeNeg"), straTofResultsOmega.neg.measuredTimeOuter * NanoToPico, straTofResultsOmega.neg.expectedTimeOuter * NanoToPico);
        }
        if (straTofResultsOmega.bach.hasOuterTof) {
          histos.fill(HIST("Omega/hOuterArrivalTimeDeltaBach"), (straTofResultsOmega.bach.measuredTimeOuter - straTofResultsOmega.bach.expectedTimeOuter) * NanoToPico);
          histos.fill(HIST("Omega/hOuterNSigmaBach"), straTofResultsOmega.bach.nSigmaOuter);
          histos.fill(HIST("Omega/hOuterExpectedArrivalTimeBach"), straTofResultsOmega.bach.expectedTimeOuter * NanoToPico);
          histos.fill(HIST("Omega/hOuterMeasuredArrivalTimeBach"), straTofResultsOmega.bach.measuredTimeOuter * NanoToPico);
          histos.fill(HIST("Omega/hOuterMeasuredVsExpectedArrivalTimeBach"), straTofResultsOmega.bach.measuredTimeOuter * NanoToPico, straTofResultsOmega.bach.expectedTimeOuter * NanoToPico);
        }
      }

      if (produce.histosAntiOmega) {
        if (straTofResultsAntiOmega.pos.hasInnerTof) {
          histos.fill(HIST("AntiOmega/hInnerArrivalTimeDeltaPos"), (straTofResultsAntiOmega.pos.measuredTimeInner - straTofResultsAntiOmega.pos.expectedTimeInner) * NanoToPico);
          histos.fill(HIST("AntiOmega/hInnerNSigmaPos"), straTofResultsAntiOmega.pos.nSigmaInner);
          histos.fill(HIST("AntiOmega/hInnerExpectedArrivalTimePos"), straTofResultsAntiOmega.pos.expectedTimeInner * NanoToPico);
          histos.fill(HIST("AntiOmega/hInnerMeasuredArrivalTimePos"), straTofResultsAntiOmega.pos.measuredTimeInner * NanoToPico);
          histos.fill(HIST("AntiOmega/hInnerMeasuredVsExpectedArrivalTimePos"), straTofResultsAntiOmega.pos.measuredTimeInner * NanoToPico, straTofResultsAntiOmega.pos.expectedTimeInner * NanoToPico);
        }
        if (straTofResultsAntiOmega.neg.hasInnerTof) {
          histos.fill(HIST("AntiOmega/hInnerArrivalTimeDeltaNeg"), (straTofResultsAntiOmega.neg.measuredTimeInner - straTofResultsAntiOmega.neg.expectedTimeInner) * NanoToPico);
          histos.fill(HIST("AntiOmega/hInnerNSigmaNeg"), straTofResultsAntiOmega.neg.nSigmaInner);
          histos.fill(HIST("AntiOmega/hInnerExpectedArrivalTimeNeg"), straTofResultsAntiOmega.neg.expectedTimeInner * NanoToPico);
          histos.fill(HIST("AntiOmega/hInnerMeasuredArrivalTimeNeg"), straTofResultsAntiOmega.neg.measuredTimeInner * NanoToPico);
          histos.fill(HIST("AntiOmega/hInnerMeasuredVsExpectedArrivalTimeNeg"), straTofResultsAntiOmega.neg.measuredTimeInner * NanoToPico, straTofResultsAntiOmega.neg.expectedTimeInner * NanoToPico);
        }
        if (straTofResultsAntiOmega.bach.hasInnerTof) {
          histos.fill(HIST("AntiOmega/hInnerArrivalTimeDeltaBach"), (straTofResultsAntiOmega.bach.measuredTimeInner - straTofResultsAntiOmega.bach.expectedTimeInner) * NanoToPico);
          histos.fill(HIST("AntiOmega/hInnerNSigmaBach"), straTofResultsAntiOmega.bach.nSigmaInner);
          histos.fill(HIST("AntiOmega/hInnerExpectedArrivalTimeBach"), straTofResultsAntiOmega.bach.expectedTimeInner * NanoToPico);
          histos.fill(HIST("AntiOmega/hInnerMeasuredArrivalTimeBach"), straTofResultsAntiOmega.bach.measuredTimeInner * NanoToPico);
          histos.fill(HIST("AntiOmega/hInnerMeasuredVsExpectedArrivalTimeBach"), straTofResultsAntiOmega.bach.measuredTimeInner * NanoToPico, straTofResultsAntiOmega.bach.expectedTimeInner * NanoToPico);
        }
        if (straTofResultsAntiOmega.pos.hasOuterTof) {
          histos.fill(HIST("AntiOmega/hOuterArrivalTimeDeltaPos"), (straTofResultsAntiOmega.pos.measuredTimeOuter - straTofResultsAntiOmega.pos.expectedTimeOuter) * NanoToPico);
          histos.fill(HIST("AntiOmega/hOuterNSigmaPos"), straTofResultsAntiOmega.pos.nSigmaOuter);
          histos.fill(HIST("AntiOmega/hOuterExpectedArrivalTimePos"), straTofResultsAntiOmega.pos.expectedTimeOuter * NanoToPico);
          histos.fill(HIST("AntiOmega/hOuterMeasuredArrivalTimePos"), straTofResultsAntiOmega.pos.measuredTimeOuter * NanoToPico);
          histos.fill(HIST("AntiOmega/hOuterMeasuredVsExpectedArrivalTimePos"), straTofResultsAntiOmega.pos.measuredTimeOuter * NanoToPico, straTofResultsAntiOmega.pos.expectedTimeOuter * NanoToPico);
        }
        if (straTofResultsAntiOmega.neg.hasOuterTof) {
          histos.fill(HIST("AntiOmega/hOuterArrivalTimeDeltaNeg"), (straTofResultsAntiOmega.neg.measuredTimeOuter - straTofResultsAntiOmega.neg.expectedTimeOuter) * NanoToPico);
          histos.fill(HIST("AntiOmega/hOuterNSigmaNeg"), straTofResultsAntiOmega.neg.nSigmaOuter);
          histos.fill(HIST("AntiOmega/hOuterExpectedArrivalTimeNeg"), straTofResultsAntiOmega.neg.expectedTimeOuter * NanoToPico);
          histos.fill(HIST("AntiOmega/hOuterMeasuredArrivalTimeNeg"), straTofResultsAntiOmega.neg.measuredTimeOuter * NanoToPico);
          histos.fill(HIST("AntiOmega/hOuterMeasuredVsExpectedArrivalTimeNeg"), straTofResultsAntiOmega.neg.measuredTimeOuter * NanoToPico, straTofResultsAntiOmega.neg.expectedTimeOuter * NanoToPico);
        }
        if (straTofResultsAntiOmega.bach.hasOuterTof) {
          histos.fill(HIST("AntiOmega/hOuterArrivalTimeDeltaBach"), (straTofResultsAntiOmega.bach.measuredTimeOuter - straTofResultsAntiOmega.bach.expectedTimeOuter) * NanoToPico);
          histos.fill(HIST("AntiOmega/hOuterNSigmaBach"), straTofResultsAntiOmega.bach.nSigmaOuter);
          histos.fill(HIST("AntiOmega/hOuterExpectedArrivalTimeBach"), straTofResultsAntiOmega.bach.expectedTimeOuter * NanoToPico);
          histos.fill(HIST("AntiOmega/hOuterMeasuredArrivalTimeBach"), straTofResultsAntiOmega.bach.measuredTimeOuter * NanoToPico);
          histos.fill(HIST("AntiOmega/hOuterMeasuredVsExpectedArrivalTimeBach"), straTofResultsAntiOmega.bach.measuredTimeOuter * NanoToPico, straTofResultsAntiOmega.bach.expectedTimeOuter * NanoToPico);
        }
      }

      if (produce.tableInnerNSigmaXi) {
        tableA3XiInnerTofPid(
          (straTofResultsXi.bach.hasInnerTof) ? straTofResultsXi.bach.nSigmaInner : o2::upgrade::pid::NoPidSignal,
          (straTofResultsXi.pos.hasInnerTof) ? straTofResultsXi.bach.nSigmaInner : o2::upgrade::pid::NoPidSignal,
          (straTofResultsAntiXi.pos.hasInnerTof) ? straTofResultsAntiXi.bach.nSigmaInner : o2::upgrade::pid::NoPidSignal,
          (straTofResultsAntiXi.neg.hasInnerTof) ? straTofResultsAntiXi.bach.nSigmaInner : o2::upgrade::pid::NoPidSignal,
          (straTofResultsXi.neg.hasInnerTof) ? straTofResultsXi.bach.nSigmaInner : o2::upgrade::pid::NoPidSignal);
      }
      if (produce.tableOuterNSigmaXi) {
        tableA3XiOuterTofPid(
          (straTofResultsXi.bach.hasOuterTof) ? straTofResultsXi.bach.nSigmaOuter : o2::upgrade::pid::NoPidSignal,
          (straTofResultsXi.pos.hasOuterTof) ? straTofResultsXi.bach.nSigmaOuter : o2::upgrade::pid::NoPidSignal,
          (straTofResultsAntiXi.pos.hasOuterTof) ? straTofResultsAntiXi.bach.nSigmaOuter : o2::upgrade::pid::NoPidSignal,
          (straTofResultsAntiXi.neg.hasOuterTof) ? straTofResultsAntiXi.bach.nSigmaOuter : o2::upgrade::pid::NoPidSignal,
          (straTofResultsXi.neg.hasOuterTof) ? straTofResultsXi.bach.nSigmaOuter : o2::upgrade::pid::NoPidSignal);
      }
      if (produce.tableInnerNSigmaOmega) {
        tableA3OmegaInnerTofPid(
          (straTofResultsOmega.bach.hasInnerTof) ? straTofResultsOmega.bach.nSigmaInner : o2::upgrade::pid::NoPidSignal,
          (straTofResultsOmega.pos.hasInnerTof) ? straTofResultsOmega.bach.nSigmaInner : o2::upgrade::pid::NoPidSignal,
          (straTofResultsAntiOmega.pos.hasInnerTof) ? straTofResultsAntiOmega.bach.nSigmaInner : o2::upgrade::pid::NoPidSignal,
          (straTofResultsAntiOmega.neg.hasInnerTof) ? straTofResultsAntiOmega.bach.nSigmaInner : o2::upgrade::pid::NoPidSignal,
          (straTofResultsOmega.neg.hasInnerTof) ? straTofResultsOmega.bach.nSigmaInner : o2::upgrade::pid::NoPidSignal);
      }
      if (produce.tableOuterNSigmaOmega) {
        tableA3OmegaOuterTofPid(
          (straTofResultsOmega.bach.hasOuterTof) ? straTofResultsOmega.bach.nSigmaOuter : o2::upgrade::pid::NoPidSignal,
          (straTofResultsOmega.pos.hasOuterTof) ? straTofResultsOmega.bach.nSigmaOuter : o2::upgrade::pid::NoPidSignal,
          (straTofResultsAntiOmega.pos.hasOuterTof) ? straTofResultsAntiOmega.bach.nSigmaOuter : o2::upgrade::pid::NoPidSignal,
          (straTofResultsAntiOmega.neg.hasOuterTof) ? straTofResultsAntiOmega.bach.nSigmaOuter : o2::upgrade::pid::NoPidSignal,
          (straTofResultsOmega.neg.hasOuterTof) ? straTofResultsOmega.bach.nSigmaOuter : o2::upgrade::pid::NoPidSignal);
      }

      if (produce.tableInnerExpectedTimesXi) {
        tableA3XiExpectedInnerTimes(
          straTofResultsXi.bach.expectedTimeInner,
          straTofResultsXi.pos.expectedTimeInner,
          straTofResultsAntiXi.pos.expectedTimeInner,
          straTofResultsAntiXi.neg.expectedTimeInner,
          straTofResultsXi.neg.expectedTimeInner,
          straTofResultsXi.bach.measuredTimeInner,
          straTofResultsXi.pos.measuredTimeInner,
          straTofResultsAntiXi.pos.measuredTimeInner,
          straTofResultsAntiXi.neg.measuredTimeInner,
          straTofResultsXi.neg.measuredTimeInner);
      }
      if (produce.tableOuterExpectedTimesXi) {
        tableA3XiExpectedOuterTimes(
          straTofResultsXi.bach.expectedTimeOuter,
          straTofResultsXi.pos.expectedTimeOuter,
          straTofResultsAntiXi.pos.expectedTimeOuter,
          straTofResultsAntiXi.neg.expectedTimeOuter,
          straTofResultsXi.neg.expectedTimeOuter,
          straTofResultsXi.bach.measuredTimeOuter,
          straTofResultsXi.pos.measuredTimeOuter,
          straTofResultsAntiXi.pos.measuredTimeOuter,
          straTofResultsAntiXi.neg.measuredTimeOuter,
          straTofResultsXi.neg.measuredTimeOuter);
      }
      if (produce.tableInnerExpectedTimesOmega) {
        tableA3OmegaExpectedInnerTimes(
          straTofResultsOmega.bach.expectedTimeInner,
          straTofResultsOmega.pos.expectedTimeInner,
          straTofResultsAntiOmega.pos.expectedTimeInner,
          straTofResultsAntiOmega.neg.expectedTimeInner,
          straTofResultsOmega.neg.expectedTimeInner,
          straTofResultsOmega.bach.measuredTimeInner,
          straTofResultsOmega.pos.measuredTimeInner,
          straTofResultsAntiOmega.pos.measuredTimeInner,
          straTofResultsAntiOmega.neg.measuredTimeInner,
          straTofResultsOmega.neg.measuredTimeInner);
      }
      if (produce.tableOuterExpectedTimesOmega) {
        tableA3OmegaExpectedOuterTimes(
          straTofResultsOmega.bach.expectedTimeOuter,
          straTofResultsOmega.pos.expectedTimeOuter,
          straTofResultsAntiOmega.pos.expectedTimeOuter,
          straTofResultsAntiOmega.neg.expectedTimeOuter,
          straTofResultsOmega.neg.expectedTimeOuter,
          straTofResultsOmega.bach.measuredTimeOuter,
          straTofResultsOmega.pos.measuredTimeOuter,
          straTofResultsAntiOmega.pos.measuredTimeOuter,
          straTofResultsAntiOmega.neg.measuredTimeOuter,
          straTofResultsOmega.neg.measuredTimeOuter);
      }
    }
  }

  PROCESS_SWITCH(Alice3StrangenessTofPid, processCascades, "", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<Alice3StrangenessTofPid>(cfgc)};
}
