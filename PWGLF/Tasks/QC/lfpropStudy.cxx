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
/// \file lfpropStudy.cxx
/// \since 27-11-2023
/// \author Carolina Reetz <c.reetz@cern.ch>
/// \brief QA task to study properties of propagated tracks

#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"

#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

using TracksType = soa::Join<aod::StoredTracks, aod::TracksDCA, aod::TracksExtra>;
using TracksCovMatType = soa::Join<TracksType, aod::StoredTracksCov, aod::TracksDCACov>;
using TracksLabeledType = soa::Join<TracksType, aod::McTrackLabels>;

struct lfpropStudy {

  ConfigurableAxis axisDCAxy{"axisDCAxy", {1000, -10.f, 10.f}, "DCA_{xy} (cm)"};
  ConfigurableAxis axisDCAz{"axisDCAz", {1000, -10.f, 10.f}, "DCA_{z} (cm)"};
  ConfigurableAxis axisMom{"axisMom", {1000, 0.0f, 10.f}, "p (GeV/c)"};

  Configurable<bool> d_noITS{"d_noITS", true, "flag to select tracks without ITS"};
  Configurable<float> d_pTMin{"d_pTMin", 0.3, "minimum track momentum"};
  Configurable<float> d_TPCrowsMin{"d_TPCrowsMin", 70, "minimum number of TPC crossed rows"};
  Configurable<float> d_TPCrowsOverFindMin{"d_TPCrowsOverFindMin", 0.8, "minimum for ratio of TPC crossed rows over findable"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  Filter collisionFilter = (aod::evsel::sel8 == true);

  void init(InitContext const&)
  {
    histos.add("hEventCounter", "hEventCounter", kTH1F, {{1, 0.0f, 1.0f}});

    const AxisSpec AxisPx{axisMom, "#it{p}_{x} (GeV/#it{c})"};
    const AxisSpec AxisPy{axisMom, "#it{p}_{y} (GeV/#it{c})"};
    const AxisSpec AxisPz{axisMom, "#it{p}_{z} (GeV/#it{c})"};
    const AxisSpec AxisPt{axisMom, "#it{p}_{T} (GeV/#it{c})"};

    // momentum
    histos.add("IU/hPx", "IU Px All", kTH1F, {AxisPx});
    histos.add("IU/hPy", "IU Py All", kTH1F, {AxisPy});
    histos.add("IU/hPz", "IU Pz All", kTH1F, {AxisPz});
    histos.add("IU/hPt", "IU Pt All", kTH1F, {AxisPt});
    histos.add("IU/El/hPx", "IU Px El", kTH1F, {AxisPx});
    histos.add("IU/El/hPy", "IU Py El", kTH1F, {AxisPy});
    histos.add("IU/El/hPz", "IU Pz El", kTH1F, {AxisPz});
    histos.add("IU/El/hPt", "IU Pt El", kTH1F, {AxisPt});
    histos.add("IU/Mu/hPx", "IU Px Mu", kTH1F, {AxisPx});
    histos.add("IU/Mu/hPy", "IU Py Mu", kTH1F, {AxisPy});
    histos.add("IU/Mu/hPz", "IU Pz Mu", kTH1F, {AxisPz});
    histos.add("IU/Mu/hPt", "IU Pt Mu", kTH1F, {AxisPt});
    histos.add("IU/Pi/hPx", "IU Px Pi", kTH1F, {AxisPx});
    histos.add("IU/Pi/hPy", "IU Py Pi", kTH1F, {AxisPy});
    histos.add("IU/Pi/hPz", "IU Pz Pi", kTH1F, {AxisPz});
    histos.add("IU/Pi/hPt", "IU Pt Pi", kTH1F, {AxisPt});
    histos.add("IU/Ka/hPx", "IU Px Ka", kTH1F, {AxisPx});
    histos.add("IU/Ka/hPy", "IU Py Ka", kTH1F, {AxisPy});
    histos.add("IU/Ka/hPz", "IU Pz Ka", kTH1F, {AxisPz});
    histos.add("IU/Ka/hPt", "IU Pt Ka", kTH1F, {AxisPt});
    histos.add("IU/Pr/hPx", "IU Px Pr", kTH1F, {AxisPx});
    histos.add("IU/Pr/hPy", "IU Py Pr", kTH1F, {AxisPy});
    histos.add("IU/Pr/hPz", "IU Pz Pr", kTH1F, {AxisPz});
    histos.add("IU/Pr/hPt", "IU Pt Pr", kTH1F, {AxisPt});
    histos.add("IU/De/hPx", "IU Px De", kTH1F, {AxisPx});
    histos.add("IU/De/hPy", "IU Py De", kTH1F, {AxisPy});
    histos.add("IU/De/hPz", "IU Pz De", kTH1F, {AxisPz});
    histos.add("IU/De/hPt", "IU Pt De", kTH1F, {AxisPt});
    histos.add("IU/Tr/hPx", "IU Px Tr", kTH1F, {AxisPx});
    histos.add("IU/Tr/hPy", "IU Py Tr", kTH1F, {AxisPy});
    histos.add("IU/Tr/hPz", "IU Pz Tr", kTH1F, {AxisPz});
    histos.add("IU/Tr/hPt", "IU Pt Tr", kTH1F, {AxisPt});
    histos.add("IU/He/hPx", "IU Px He", kTH1F, {AxisPx});
    histos.add("IU/He/hPy", "IU Py He", kTH1F, {AxisPy});
    histos.add("IU/He/hPz", "IU Pz He", kTH1F, {AxisPz});
    histos.add("IU/He/hPt", "IU Pt He", kTH1F, {AxisPt});
    histos.add("IU/Al/hPx", "IU Px Al", kTH1F, {AxisPx});
    histos.add("IU/Al/hPy", "IU Py Al", kTH1F, {AxisPy});
    histos.add("IU/Al/hPz", "IU Pz Al", kTH1F, {AxisPz});
    histos.add("IU/Al/hPt", "IU Pt Al", kTH1F, {AxisPt});

    histos.add("hPx", "Propagated Px All", kTH1F, {AxisPx});
    histos.add("hPy", "Propagated Py All", kTH1F, {AxisPy});
    histos.add("hPz", "Propagated Pz All", kTH1F, {AxisPz});
    histos.add("hPt", "Propagated Pt All", kTH1F, {AxisPt});
    histos.add("El/hPx", "Propagated Px El", kTH1F, {AxisPx});
    histos.add("El/hPy", "Propagated Py El", kTH1F, {AxisPy});
    histos.add("El/hPz", "Propagated Pz El", kTH1F, {AxisPz});
    histos.add("El/hPt", "Propagated Pt El", kTH1F, {AxisPt});
    histos.add("Mu/hPx", "Propagated Px Mu", kTH1F, {AxisPx});
    histos.add("Mu/hPy", "Propagated Py Mu", kTH1F, {AxisPy});
    histos.add("Mu/hPz", "Propagated Pz Mu", kTH1F, {AxisPz});
    histos.add("Mu/hPt", "Propagated Pt Mu", kTH1F, {AxisPt});
    histos.add("Pi/hPx", "Propagated Px Pi", kTH1F, {AxisPx});
    histos.add("Pi/hPy", "Propagated Py Pi", kTH1F, {AxisPy});
    histos.add("Pi/hPz", "Propagated Pz Pi", kTH1F, {AxisPz});
    histos.add("Pi/hPt", "Propagated Pt Pi", kTH1F, {AxisPt});
    histos.add("Ka/hPx", "Propagated Px Ka", kTH1F, {AxisPx});
    histos.add("Ka/hPy", "Propagated Py Ka", kTH1F, {AxisPy});
    histos.add("Ka/hPz", "Propagated Pz Ka", kTH1F, {AxisPz});
    histos.add("Ka/hPt", "Propagated Pt Ka", kTH1F, {AxisPt});
    histos.add("Pr/hPx", "Propagated Px Pr", kTH1F, {AxisPx});
    histos.add("Pr/hPy", "Propagated Py Pr", kTH1F, {AxisPy});
    histos.add("Pr/hPz", "Propagated Pz Pr", kTH1F, {AxisPz});
    histos.add("Pr/hPt", "Propagated Pt Pr", kTH1F, {AxisPt});
    histos.add("De/hPx", "Propagated Px De", kTH1F, {AxisPx});
    histos.add("De/hPy", "Propagated Py De", kTH1F, {AxisPy});
    histos.add("De/hPz", "Propagated Pz De", kTH1F, {AxisPz});
    histos.add("De/hPt", "Propagated Pt De", kTH1F, {AxisPt});
    histos.add("Tr/hPx", "Propagated Px Tr", kTH1F, {AxisPx});
    histos.add("Tr/hPy", "Propagated Py Tr", kTH1F, {AxisPy});
    histos.add("Tr/hPz", "Propagated Pz Tr", kTH1F, {AxisPz});
    histos.add("Tr/hPt", "Propagated Pt Tr", kTH1F, {AxisPt});
    histos.add("He/hPx", "Propagated Px He", kTH1F, {AxisPx});
    histos.add("He/hPy", "Propagated Py He", kTH1F, {AxisPy});
    histos.add("He/hPz", "Propagated Pz He", kTH1F, {AxisPz});
    histos.add("He/hPt", "Propagated Pt He", kTH1F, {AxisPt});
    histos.add("Al/hPx", "Propagated Px Al", kTH1F, {AxisPx});
    histos.add("Al/hPy", "Propagated Py Al", kTH1F, {AxisPy});
    histos.add("Al/hPz", "Propagated Pz Al", kTH1F, {AxisPz});
    histos.add("Al/hPt", "Propagated Pt Al", kTH1F, {AxisPt});

    const AxisSpec AxisDCAxy{axisDCAxy, "DCA_{xy} (cm)"};
    const AxisSpec AxisDCAz{axisDCAz, "DCA_{z} (cm)"};

    // DCA
    histos.add("IU/hDCAxy", "IU DCAxy All", kTH1F, {AxisDCAxy});
    histos.add("IU/hDCAz", "IU DCAz All", kTH1F, {AxisDCAz});
    histos.add("IU/El/hDCAxy", "IU DCAxy El", kTH1F, {AxisDCAxy});
    histos.add("IU/El/hDCAz", "IU DCAz El", kTH1F, {AxisDCAz});
    histos.add("IU/Mu/hDCAxy", "IU DCAxy Mu", kTH1F, {AxisDCAxy});
    histos.add("IU/Mu/hDCAz", "IU DCAz Mu", kTH1F, {AxisDCAz});
    histos.add("IU/Pi/hDCAxy", "IU DCAxy Pi", kTH1F, {AxisDCAxy});
    histos.add("IU/Pi/hDCAz", "IU DCAz Pi", kTH1F, {AxisDCAz});
    histos.add("IU/Ka/hDCAxy", "IU DCAxy Ka", kTH1F, {AxisDCAxy});
    histos.add("IU/Ka/hDCAz", "IU DCAz Ka", kTH1F, {AxisDCAz});
    histos.add("IU/Pr/hDCAxy", "IU DCAxy Pr", kTH1F, {AxisDCAxy});
    histos.add("IU/Pr/hDCAz", "IU DCAz Pr", kTH1F, {AxisDCAz});
    histos.add("IU/De/hDCAxy", "IU DCAxy De", kTH1F, {AxisDCAxy});
    histos.add("IU/De/hDCAz", "IU DCAz De", kTH1F, {AxisDCAz});
    histos.add("IU/Tr/hDCAxy", "IU DCAxy Tr", kTH1F, {AxisDCAxy});
    histos.add("IU/Tr/hDCAz", "IU DCAz Tr", kTH1F, {AxisDCAz});
    histos.add("IU/He/hDCAxy", "IU DCAxy He", kTH1F, {AxisDCAxy});
    histos.add("IU/He/hDCAz", "IU DCAz He", kTH1F, {AxisDCAz});
    histos.add("IU/Al/hDCAxy", "IU DCAxy Al", kTH1F, {AxisDCAxy});
    histos.add("IU/Al/hDCAz", "IU DCAz Al", kTH1F, {AxisDCAz});

    histos.add("hDCAxy", "Propagated DCAxy All", kTH1F, {AxisDCAxy});
    histos.add("hDCAz", "Propagated DCAz All", kTH1F, {AxisDCAz});
    histos.add("El/hDCAxy", "Propagated DCAxy El", kTH1F, {AxisDCAxy});
    histos.add("El/hDCAz", "Propagated DCAz El", kTH1F, {AxisDCAz});
    histos.add("Mu/hDCAxy", "Propagated DCAxy Mu", kTH1F, {AxisDCAxy});
    histos.add("Mu/hDCAz", "Propagated DCAz Mu", kTH1F, {AxisDCAz});
    histos.add("Pi/hDCAxy", "Propagated DCAxy Pi", kTH1F, {AxisDCAxy});
    histos.add("Pi/hDCAz", "Propagated DCAz Pi", kTH1F, {AxisDCAz});
    histos.add("Ka/hDCAxy", "Propagated DCAxy Ka", kTH1F, {AxisDCAxy});
    histos.add("Ka/hDCAz", "Propagated DCAz Ka", kTH1F, {AxisDCAz});
    histos.add("Pr/hDCAxy", "Propagated DCAxy Pr", kTH1F, {AxisDCAxy});
    histos.add("Pr/hDCAz", "Propagated DCAz Pr", kTH1F, {AxisDCAz});
    histos.add("De/hDCAxy", "Propagated DCAxy De", kTH1F, {AxisDCAxy});
    histos.add("De/hDCAz", "Propagated DCAz De", kTH1F, {AxisDCAz});
    histos.add("Tr/hDCAxy", "Propagated DCAxy Tr", kTH1F, {AxisDCAxy});
    histos.add("Tr/hDCAz", "Propagated DCAz Tr", kTH1F, {AxisDCAz});
    histos.add("He/hDCAxy", "Propagated DCAxy He", kTH1F, {AxisDCAxy});
    histos.add("He/hDCAz", "Propagated DCAz He", kTH1F, {AxisDCAz});
    histos.add("Al/hDCAxy", "Propagated DCAxy Al", kTH1F, {AxisDCAxy});
    histos.add("Al/hDCAz", "Propagated DCAz Al", kTH1F, {AxisDCAz});
  }

  template <typename T>
  void fillPerTrack(const T& track)
  {
    histos.fill(HIST("hPx"), track.px());
    histos.fill(HIST("hPy"), track.py());
    histos.fill(HIST("hPz"), track.pz());
    histos.fill(HIST("hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
    histos.fill(HIST("hDCAxy"), track.dcaXY());
    histos.fill(HIST("hDCAz"), track.dcaZ());

    switch (track.pidForTracking()) {
      case o2::track::PID::Electron:
        histos.fill(HIST("El/hPx"), track.px());
        histos.fill(HIST("El/hPy"), track.py());
        histos.fill(HIST("El/hPz"), track.pz());
        histos.fill(HIST("El/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("El/hDCAxy"), track.dcaXY());
        histos.fill(HIST("El/hDCAz"), track.dcaZ());
        break;
      case o2::track::PID::Muon:
        histos.fill(HIST("Mu/hPx"), track.px());
        histos.fill(HIST("Mu/hPy"), track.py());
        histos.fill(HIST("Mu/hPz"), track.pz());
        histos.fill(HIST("Mu/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("Mu/hDCAxy"), track.dcaXY());
        histos.fill(HIST("Mu/hDCAz"), track.dcaZ());
        break;
      case o2::track::PID::Pion:
        histos.fill(HIST("Pi/hPx"), track.px());
        histos.fill(HIST("Pi/hPy"), track.py());
        histos.fill(HIST("Pi/hPz"), track.pz());
        histos.fill(HIST("Pi/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("Pi/hDCAxy"), track.dcaXY());
        histos.fill(HIST("Pi/hDCAz"), track.dcaZ());
        break;
      case o2::track::PID::Kaon:
        histos.fill(HIST("Ka/hPx"), track.px());
        histos.fill(HIST("Ka/hPy"), track.py());
        histos.fill(HIST("Ka/hPz"), track.pz());
        histos.fill(HIST("Ka/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("Ka/hDCAxy"), track.dcaXY());
        histos.fill(HIST("Ka/hDCAz"), track.dcaZ());
        break;
      case o2::track::PID::Proton:
        histos.fill(HIST("Pr/hPx"), track.px());
        histos.fill(HIST("Pr/hPy"), track.py());
        histos.fill(HIST("Pr/hPz"), track.pz());
        histos.fill(HIST("Pr/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("Pr/hDCAxy"), track.dcaXY());
        histos.fill(HIST("Pr/hDCAz"), track.dcaZ());
        break;
      case o2::track::PID::Deuteron:
        histos.fill(HIST("De/hPx"), track.px());
        histos.fill(HIST("De/hPy"), track.py());
        histos.fill(HIST("De/hPz"), track.pz());
        histos.fill(HIST("De/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("De/hDCAxy"), track.dcaXY());
        histos.fill(HIST("De/hDCAz"), track.dcaZ());
        break;
      case o2::track::PID::Triton:
        histos.fill(HIST("Tr/hPx"), track.px());
        histos.fill(HIST("Tr/hPy"), track.py());
        histos.fill(HIST("Tr/hPz"), track.pz());
        histos.fill(HIST("Tr/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("Tr/hDCAxy"), track.dcaXY());
        histos.fill(HIST("Tr/hDCAz"), track.dcaZ());
        break;
      case o2::track::PID::Helium3:
        histos.fill(HIST("He/hPx"), track.px());
        histos.fill(HIST("He/hPy"), track.py());
        histos.fill(HIST("He/hPz"), track.pz());
        histos.fill(HIST("He/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("He/hDCAxy"), track.dcaXY());
        histos.fill(HIST("He/hDCAz"), track.dcaZ());
        break;
      case o2::track::PID::Alpha:
        histos.fill(HIST("Al/hPx"), track.px());
        histos.fill(HIST("Al/hPy"), track.py());
        histos.fill(HIST("Al/hPz"), track.pz());
        histos.fill(HIST("Al/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("Al/hDCAxy"), track.dcaXY());
        histos.fill(HIST("Al/hDCAz"), track.dcaZ());
        break;
      default:
        LOG(fatal) << "Unknown PID: " << track.pidForTracking();
    }
  }

  template <typename T>
  void fillPerTrackIU(const T& track)
  {
    histos.fill(HIST("IU/hPx"), track.px());
    histos.fill(HIST("IU/hPy"), track.py());
    histos.fill(HIST("IU/hPz"), track.pz());
    histos.fill(HIST("IU/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
    histos.fill(HIST("IU/hDCAxy"), track.dcaXY());
    histos.fill(HIST("IU/hDCAz"), track.dcaZ());

    switch (track.pidForTracking()) {
      case o2::track::PID::Electron:
        histos.fill(HIST("IU/El/hPx"), track.px());
        histos.fill(HIST("IU/El/hPy"), track.py());
        histos.fill(HIST("IU/El/hPz"), track.pz());
        histos.fill(HIST("IU/El/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("IU/El/hDCAxy"), track.dcaXY());
        histos.fill(HIST("IU/El/hDCAz"), track.dcaZ());
        break;
      case o2::track::PID::Muon:
        histos.fill(HIST("IU/Mu/hPx"), track.px());
        histos.fill(HIST("IU/Mu/hPy"), track.py());
        histos.fill(HIST("IU/Mu/hPz"), track.pz());
        histos.fill(HIST("IU/Mu/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("IU/Mu/hDCAxy"), track.dcaXY());
        histos.fill(HIST("IU/Mu/hDCAz"), track.dcaZ());
        break;
      case o2::track::PID::Pion:
        histos.fill(HIST("IU/Pi/hPx"), track.px());
        histos.fill(HIST("IU/Pi/hPy"), track.py());
        histos.fill(HIST("IU/Pi/hPz"), track.pz());
        histos.fill(HIST("IU/Pi/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("IU/Pi/hDCAxy"), track.dcaXY());
        histos.fill(HIST("IU/Pi/hDCAz"), track.dcaZ());
        break;
      case o2::track::PID::Kaon:
        histos.fill(HIST("IU/Ka/hPx"), track.px());
        histos.fill(HIST("IU/Ka/hPy"), track.py());
        histos.fill(HIST("IU/Ka/hPz"), track.pz());
        histos.fill(HIST("IU/Ka/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("IU/Ka/hDCAxy"), track.dcaXY());
        histos.fill(HIST("IU/Ka/hDCAz"), track.dcaZ());
        break;
      case o2::track::PID::Proton:
        histos.fill(HIST("IU/Pr/hPx"), track.px());
        histos.fill(HIST("IU/Pr/hPy"), track.py());
        histos.fill(HIST("IU/Pr/hPz"), track.pz());
        histos.fill(HIST("IU/Pr/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("IU/Pr/hDCAxy"), track.dcaXY());
        histos.fill(HIST("IU/Pr/hDCAz"), track.dcaZ());
        break;
      case o2::track::PID::Deuteron:
        histos.fill(HIST("IU/De/hPx"), track.px());
        histos.fill(HIST("IU/De/hPy"), track.py());
        histos.fill(HIST("IU/De/hPz"), track.pz());
        histos.fill(HIST("IU/De/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("IU/De/hDCAxy"), track.dcaXY());
        histos.fill(HIST("IU/De/hDCAz"), track.dcaZ());
        break;
      case o2::track::PID::Triton:
        histos.fill(HIST("IU/Tr/hPx"), track.px());
        histos.fill(HIST("IU/Tr/hPy"), track.py());
        histos.fill(HIST("IU/Tr/hPz"), track.pz());
        histos.fill(HIST("IU/Tr/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("IU/Tr/hDCAxy"), track.dcaXY());
        histos.fill(HIST("IU/Tr/hDCAz"), track.dcaZ());
        break;
      case o2::track::PID::Helium3:
        histos.fill(HIST("IU/He/hPx"), track.px());
        histos.fill(HIST("IU/He/hPy"), track.py());
        histos.fill(HIST("IU/He/hPz"), track.pz());
        histos.fill(HIST("IU/He/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("IU/He/hDCAxy"), track.dcaXY());
        histos.fill(HIST("IU/He/hDCAz"), track.dcaZ());
        break;
      case o2::track::PID::Alpha:
        histos.fill(HIST("IU/Al/hPx"), track.px());
        histos.fill(HIST("IU/Al/hPy"), track.py());
        histos.fill(HIST("IU/Al/hPz"), track.pz());
        histos.fill(HIST("IU/Al/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("IU/Al/hDCAxy"), track.dcaXY());
        histos.fill(HIST("IU/Al/hDCAz"), track.dcaZ());
        break;
      default:
        LOG(fatal) << "Unknown PID: " << track.pidForTracking();
    }
  }

  void processData(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& /*collision*/,
                   TracksType const& tracks)
  {
    histos.fill(HIST("hEventCounter"), 0.5);
    for (auto& track : tracks) {
      // track selection
      if (d_noITS && track.hasITS())
        continue; // optional: only look at tracks which start outside the ITS
      if (track.tpcNClsCrossedRows() < d_TPCrowsMin || track.tpcCrossedRowsOverFindableCls() < d_TPCrowsOverFindMin)
        continue;
      if (track.trackType() == aod::track::TrackIU) { // only look at tracks which were propagated
        fillPerTrackIU(track);
        continue;
      }
      fillPerTrack(track);
    }
  }
  PROCESS_SWITCH(lfpropStudy, processData, "process data", true);

  void processDataNoColl(TracksType const& tracks)
  {
    histos.fill(HIST("hEventCounter"), 0.5);
    for (auto& track : tracks) {
      // track selection
      if (d_noITS && track.hasITS())
        continue; // optional: only look at tracks which start outside the ITS
      if (track.tpcNClsCrossedRows() < d_TPCrowsMin || track.tpcCrossedRowsOverFindableCls() < d_TPCrowsOverFindMin)
        continue;
      if (track.trackType() == aod::track::TrackIU) { // only look at tracks which were propagated
        fillPerTrackIU(track);
        continue;
      }
      fillPerTrack(track);
    }
  }
  PROCESS_SWITCH(lfpropStudy, processDataNoColl, "process data without collisions association", false);

  void processDataCovMat(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& /*collision*/,
                         TracksCovMatType const& tracks)
  {
    histos.fill(HIST("hEventCounter"), 0.5);
    for (auto& track : tracks) {
      // track selection
      if (d_noITS && track.hasITS())
        continue; // optional: only look at tracks which start outside the ITS
      if (track.tpcNClsCrossedRows() < d_TPCrowsMin || track.tpcCrossedRowsOverFindableCls() < d_TPCrowsOverFindMin)
        continue;
      if (track.trackType() == aod::track::TrackIU) { // only look at tracks which were propagated
        fillPerTrackIU(track);
        continue;
      }
      fillPerTrack(track);
    }
  }
  PROCESS_SWITCH(lfpropStudy, processDataCovMat, "process data with Cov. Mat.", false);

  void processMC(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& /*collision*/,
                 TracksLabeledType const& tracks,
                 aod::McParticles const& particlesMC)
  {
    histos.fill(HIST("hEventCounter"), 0.5);
    for (auto& track : tracks) {
      if (!track.has_mcParticle() || track.mcParticleId() <= -1 || track.mcParticleId() > particlesMC.size())
        continue;
      if (d_noITS && track.hasITS())
        continue; // optional: only look at tracks which start outside the ITS
      if (track.tpcNClsCrossedRows() < d_TPCrowsMin || track.tpcCrossedRowsOverFindableCls() < d_TPCrowsOverFindMin)
        continue;
      if (track.trackType() == aod::track::TrackIU) { // only look at tracks which were propagated
        histos.fill(HIST("IU/hPx"), track.px());
        histos.fill(HIST("IU/hPy"), track.py());
        histos.fill(HIST("IU/hPz"), track.pz());
        histos.fill(HIST("IU/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
        histos.fill(HIST("IU/hDCAxy"), track.dcaXY());
        histos.fill(HIST("IU/hDCAz"), track.dcaZ());

        switch (std::abs(track.mcParticle().pdgCode())) {
          case kElectron:
            histos.fill(HIST("IU/El/hPx"), track.px());
            histos.fill(HIST("IU/El/hPy"), track.py());
            histos.fill(HIST("IU/El/hPz"), track.pz());
            histos.fill(HIST("IU/El/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
            histos.fill(HIST("IU/El/hDCAxy"), track.dcaXY());
            histos.fill(HIST("IU/El/hDCAz"), track.dcaZ());
            break;
          case kMuonPlus:
            histos.fill(HIST("IU/Mu/hPx"), track.px());
            histos.fill(HIST("IU/Mu/hPy"), track.py());
            histos.fill(HIST("IU/Mu/hPz"), track.pz());
            histos.fill(HIST("IU/Mu/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
            histos.fill(HIST("IU/Mu/hDCAxy"), track.dcaXY());
            histos.fill(HIST("IU/Mu/hDCAz"), track.dcaZ());
            break;
          case kPiPlus:
            histos.fill(HIST("IU/Pi/hPx"), track.px());
            histos.fill(HIST("IU/Pi/hPy"), track.py());
            histos.fill(HIST("IU/Pi/hPz"), track.pz());
            histos.fill(HIST("IU/Pi/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
            histos.fill(HIST("IU/Pi/hDCAxy"), track.dcaXY());
            histos.fill(HIST("IU/Pi/hDCAz"), track.dcaZ());
            break;
          case kKPlus:
            histos.fill(HIST("IU/Ka/hPx"), track.px());
            histos.fill(HIST("IU/Ka/hPy"), track.py());
            histos.fill(HIST("IU/Ka/hPz"), track.pz());
            histos.fill(HIST("IU/Ka/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
            histos.fill(HIST("IU/Ka/hDCAxy"), track.dcaXY());
            histos.fill(HIST("IU/Ka/hDCAz"), track.dcaZ());
            break;
          case kProton:
            histos.fill(HIST("IU/Pr/hPx"), track.px());
            histos.fill(HIST("IU/Pr/hPy"), track.py());
            histos.fill(HIST("IU/Pr/hPz"), track.pz());
            histos.fill(HIST("IU/Pr/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
            histos.fill(HIST("IU/Pr/hDCAxy"), track.dcaXY());
            histos.fill(HIST("IU/Pr/hDCAz"), track.dcaZ());
            break;
          case 1000010020:
            histos.fill(HIST("IU/De/hPx"), track.px());
            histos.fill(HIST("IU/De/hPy"), track.py());
            histos.fill(HIST("IU/De/hPz"), track.pz());
            histos.fill(HIST("IU/De/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
            histos.fill(HIST("IU/De/hDCAxy"), track.dcaXY());
            histos.fill(HIST("IU/De/hDCAz"), track.dcaZ());
            break;
          case 1000010030:
            histos.fill(HIST("IU/Tr/hPx"), track.px());
            histos.fill(HIST("IU/Tr/hPy"), track.py());
            histos.fill(HIST("IU/Tr/hPz"), track.pz());
            histos.fill(HIST("IU/Tr/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
            histos.fill(HIST("IU/Tr/hDCAxy"), track.dcaXY());
            histos.fill(HIST("IU/Tr/hDCAz"), track.dcaZ());
            break;
          case 1000020030:
            histos.fill(HIST("IU/He/hPx"), track.px());
            histos.fill(HIST("IU/He/hPy"), track.py());
            histos.fill(HIST("IU/He/hPz"), track.pz());
            histos.fill(HIST("IU/He/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
            histos.fill(HIST("IU/He/hDCAxy"), track.dcaXY());
            histos.fill(HIST("IU/He/hDCAz"), track.dcaZ());
            break;
          case 1000020040:
            histos.fill(HIST("IU/Al/hPx"), track.px());
            histos.fill(HIST("IU/Al/hPy"), track.py());
            histos.fill(HIST("IU/Al/hPz"), track.pz());
            histos.fill(HIST("IU/Al/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
            histos.fill(HIST("IU/Al/hDCAxy"), track.dcaXY());
            histos.fill(HIST("IU/Al/hDCAz"), track.dcaZ());
            break;
          default:
            continue;
        }
        continue;
      }

      histos.fill(HIST("hPx"), track.px());
      histos.fill(HIST("hPy"), track.py());
      histos.fill(HIST("hPz"), track.pz());
      histos.fill(HIST("hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
      histos.fill(HIST("hDCAxy"), track.dcaXY());
      histos.fill(HIST("hDCAz"), track.dcaZ());

      switch (std::abs(track.mcParticle().pdgCode())) {
        case kElectron:
          histos.fill(HIST("El/hPx"), track.px());
          histos.fill(HIST("El/hPy"), track.py());
          histos.fill(HIST("El/hPz"), track.pz());
          histos.fill(HIST("El/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("El/hDCAxy"), track.dcaXY());
          histos.fill(HIST("El/hDCAz"), track.dcaZ());
          break;
        case kMuonPlus:
          histos.fill(HIST("Mu/hPx"), track.px());
          histos.fill(HIST("Mu/hPy"), track.py());
          histos.fill(HIST("Mu/hPz"), track.pz());
          histos.fill(HIST("Mu/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("Mu/hDCAxy"), track.dcaXY());
          histos.fill(HIST("Mu/hDCAz"), track.dcaZ());
          break;
        case kPiPlus:
          histos.fill(HIST("Pi/hPx"), track.px());
          histos.fill(HIST("Pi/hPy"), track.py());
          histos.fill(HIST("Pi/hPz"), track.pz());
          histos.fill(HIST("Pi/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("Pi/hDCAxy"), track.dcaXY());
          histos.fill(HIST("Pi/hDCAz"), track.dcaZ());
          break;
        case kKPlus:
          histos.fill(HIST("Ka/hPx"), track.px());
          histos.fill(HIST("Ka/hPy"), track.py());
          histos.fill(HIST("Ka/hPz"), track.pz());
          histos.fill(HIST("Ka/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("Ka/hDCAxy"), track.dcaXY());
          histos.fill(HIST("Ka/hDCAz"), track.dcaZ());
          break;
        case kProton:
          histos.fill(HIST("Pr/hPx"), track.px());
          histos.fill(HIST("Pr/hPy"), track.py());
          histos.fill(HIST("Pr/hPz"), track.pz());
          histos.fill(HIST("Pr/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("Pr/hDCAxy"), track.dcaXY());
          histos.fill(HIST("Pr/hDCAz"), track.dcaZ());
          break;
        case 1000010020:
          histos.fill(HIST("De/hPx"), track.px());
          histos.fill(HIST("De/hPy"), track.py());
          histos.fill(HIST("De/hPz"), track.pz());
          histos.fill(HIST("De/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("De/hDCAxy"), track.dcaXY());
          histos.fill(HIST("De/hDCAz"), track.dcaZ());
          break;
        case 1000010030:
          histos.fill(HIST("Tr/hPx"), track.px());
          histos.fill(HIST("Tr/hPy"), track.py());
          histos.fill(HIST("Tr/hPz"), track.pz());
          histos.fill(HIST("Tr/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("Tr/hDCAxy"), track.dcaXY());
          histos.fill(HIST("Tr/hDCAz"), track.dcaZ());
          break;
        case 1000020030:
          histos.fill(HIST("He/hPx"), track.px());
          histos.fill(HIST("He/hPy"), track.py());
          histos.fill(HIST("He/hPz"), track.pz());
          histos.fill(HIST("He/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("He/hDCAxy"), track.dcaXY());
          histos.fill(HIST("He/hDCAz"), track.dcaZ());
          break;
        case 1000020040:
          histos.fill(HIST("Al/hPx"), track.px());
          histos.fill(HIST("Al/hPy"), track.py());
          histos.fill(HIST("Al/hPz"), track.pz());
          histos.fill(HIST("Al/hPt"), sqrt(track.px() * track.px() + track.py() * track.py()));
          histos.fill(HIST("Al/hDCAxy"), track.dcaXY());
          histos.fill(HIST("Al/hDCAz"), track.dcaZ());
          break;
        default:
          continue;
      }
    }
  }
  PROCESS_SWITCH(lfpropStudy, processMC, "process MC", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<lfpropStudy>(cfgc)}; }
