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
/// \author Alberto Caliva (alberto.caliva@cern.ch), Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \since May 31, 2024

#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TVector2.h>
#include <TVector3.h>
#include <cmath>
#include <vector>

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Framework/ASoA.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/Track.h"
#include "ReconstructionDataFormats/DCA.h"

using namespace std;
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using std::array;

using SelectedCollisions = soa::Join<aod::Collisions, aod::EvSels>;

using FullTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;

struct tracked_cascade_properties {

  // QC Histograms
  HistogramRegistry registryQC{
    "registryQC",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Analysis Histograms: Data
  HistogramRegistry registryData{
    "registryData",
    {},
    OutputObjHandlingPolicy::AnalysisObject,
    true,
    true};

  // Global Parameters
  Configurable<float> zVtx{"zVtx", 10.0f, "z vertex cut"};

  // Track Parameters
  Configurable<float> minITSnCls{"minITSnCls", 3.0f, "min number of ITS clusters"};
  Configurable<float> minTPCnClsFound{"minTPCnClsFound", 80.0f, "min number of found TPC clusters"};
  Configurable<float> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80.0f, "min number of TPC crossed rows"};
  Configurable<float> maxChi2TPC{"maxChi2TPC", 4.0f, "max chi2 per cluster TPC"};
  Configurable<float> etaMin{"etaMin", -0.8f, "eta min"};
  Configurable<float> etaMax{"etaMax", +0.8f, "eta max"};
  Configurable<float> ptMin_V0_proton{"ptMin_V0_proton", 0.3f, "pt min of proton from V0"};
  Configurable<float> ptMax_V0_proton{"ptMax_V0_proton", 10.0f, "pt max of proton from V0"};
  Configurable<float> ptMin_V0_pion{"ptMin_V0_pion", 0.1f, "pt min of pion from V0"};
  Configurable<float> ptMax_V0_pion{"ptMax_V0_pion", 1.5f, "pt max of pion from V0"};
  Configurable<float> nsigmaTPCmin{"nsigmaTPCmin", -3.0f, "Minimum nsigma TPC"};
  Configurable<float> nsigmaTPCmax{"nsigmaTPCmax", +3.0f, "Maximum nsigma TPC"};
  Configurable<float> nsigmaTOFmin{"nsigmaTOFmin", -3.0f, "Minimum nsigma TOF"};
  Configurable<float> nsigmaTOFmax{"nsigmaTOFmax", +3.0f, "Maximum nsigma TOF"};
  Configurable<float> dcanegtoPVmin{"dcanegtoPVmin", 0.0f, "Minimum DCA Neg To PV"};
  Configurable<float> dcapostoPVmin{"dcapostoPVmin", 0.0f, "Minimum DCA Pos To PV"};
  Configurable<bool> requireTOF{"requireTOF", false, "require TOF hit"};
  Configurable<bool> requireITS{"requireITS", false, "require ITS hit"};

  // V0 Parameters
  Configurable<float> yMin{"yMin", -0.5f, "minimum y"};
  Configurable<float> yMax{"yMax", +0.5f, "maximum y"};
  Configurable<float> minimumV0Radius{"minimumV0Radius", 0.0f, "Minimum V0 Radius"};
  Configurable<float> maximumV0Radius{"maximumV0Radius", 100.0f, "Maximum V0 Radius"};
  Configurable<float> v0cospaMin{"v0cospaMin", 0.99f, "Minimum V0 CosPA"};
  Configurable<float> dcaV0DaughtersMax{"dcaV0DaughtersMax", 0.5f, "Maximum DCA Daughters"};

  // Cascade Parameters
  Configurable<float> minimumCascRadius{"minimumCascRadius", 5.0f, "Minimum Cascade Radius"};
  Configurable<float> maximumCascRadius{"maximumCascRadius", 18.0f, "Maximum Cascade Radius"};
  Configurable<float> casccospaMin{"casccospaMin", 0.99f, "Minimum Cascade CosPA"};
  Configurable<float> dcabachtopvMin{"dcabachtopvMin", 0.0f, "Minimum DCA bachelor to PV"};
  Configurable<float> dcaV0topvMin{"dcaV0topvMin", 0.0f, "Minimum DCA V0 to PV"};
  Configurable<float> dcaCascDaughtersMax{"dcaCascDaughtersMax", 0.5f, "Maximum DCA Daughters"};

  // Mass Cuts
  Configurable<float> mMin_xi{"mMin_xi", 1.28f, "mMin Xi"};
  Configurable<float> mMax_xi{"mMax_xi", 1.36f, "mMax Xi"};
  Configurable<float> mMin_omega{"mMin_omega", 1.63f, "mMin Omega"};
  Configurable<float> mMax_omega{"mMax_omega", 1.71f, "mMax Omega"};

  void init(InitContext const&)
  {
    registryQC.add("number_of_events_data", "number of events in data", HistType::kTH1F, {{15, 0, 15, "Event Cuts"}});

    registryData.add("xi_pos", "xi_pos", HistType::kTH3F, {{100, 0.0, 20.0, "#LT ITS cluster size #GT"}, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}});
    registryData.add("xi_neg", "xi_neg", HistType::kTH3F, {{100, 0.0, 20.0, "#LT ITS cluster size #GT"}, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}});
    registryData.add("omega_pos", "omega_pos", HistType::kTH3F, {{100, 0.0, 20.0, "#LT ITS cluster size #GT"}, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}});
    registryData.add("omega_neg", "omega_neg", HistType::kTH3F, {{100, 0.0, 20.0, "#LT ITS cluster size #GT"}, {100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {16, -0.8, 0.8, "#eta"}});

    registryData.add("xi_mass_pos", "xi_mass_pos", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"}});
    registryData.add("xi_mass_neg", "xi_mass_neg", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"}});
    registryData.add("omega_mass_pos", "omega_mass_pos", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"}});
    registryData.add("omega_mass_neg", "omega_mass_neg", HistType::kTH2F, {{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"}, {200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"}});
  }

  // Hits on ITS Layer
  bool hasHitOnITSlayer(uint8_t itsClsmap, int layer)
  {
    unsigned char test_bit = 1 << layer;
    return (itsClsmap & test_bit);
  }

  // Xi Selections
  template <typename CascType, typename T1, typename T2, typename T3, typename C>
  bool passedXiSelection(const CascType& casc, const T1& ptrack, const T2& ntrack, const T3& btrack, const C& coll)
  {
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;
    if (!passedSingleTrackSelection(btrack))
      return false;

    // Xi+ Selection (Xi+ -> antiL + pi+)
    if (btrack.sign() > 0) {
      if (ntrack.pt() < ptMin_V0_proton)
        return false;
      if (ntrack.pt() > ptMax_V0_proton)
        return false;
      if (ptrack.pt() < ptMin_V0_pion)
        return false;
      if (ptrack.pt() > ptMax_V0_pion)
        return false;

      // PID Selections (TPC)
      if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID Selections (TOF)
      if (requireTOF) {
        if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }
    }

    // Xi- Selection (Xi- -> L + pi-)
    if (btrack.sign() < 0) {
      if (ptrack.pt() < ptMin_V0_proton)
        return false;
      if (ptrack.pt() > ptMax_V0_proton)
        return false;
      if (ntrack.pt() < ptMin_V0_pion)
        return false;
      if (ntrack.pt() > ptMax_V0_pion)
        return false;

      // PID Selections (TPC)
      if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID Selections (TOF)
      if (requireTOF) {
        if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }
    }

    // V0 Selections
    if (casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()) < v0cospaMin)
      return false;
    if (casc.v0radius() < minimumV0Radius || casc.v0radius() > maximumV0Radius)
      return false;
    if (casc.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    if (casc.dcapostopv() < dcapostoPVmin)
      return false;
    if (casc.dcanegtopv() < dcanegtoPVmin)
      return false;

    // Cascade Selections
    if (casc.cascradius() < minimumCascRadius || casc.cascradius() > maximumCascRadius)
      return false;
    if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) < casccospaMin)
      return false;
    if (casc.dcabachtopv() < dcabachtopvMin)
      return false;
    if (casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()) < dcaV0topvMin)
      return false;
    if (casc.dcacascdaughters() > dcaCascDaughtersMax)
      return false;

    // PID Selection on bachelor
    if (btrack.tpcNSigmaPi() < nsigmaTPCmin || btrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (requireTOF) {
      if (btrack.tofNSigmaPi() < nsigmaTOFmin || btrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }

    // Rapidity Selection
    /*
    TLorentzVector lorentzVect;
    lorentzVect.SetXYZM(casc.px(), casc.py(), casc.pz(), 1.32171);
    if (lorentzVect.Rapidity() < yMin || lorentzVect.Rapidity() > yMax)
      return false;
    */
    return true;
  }

  // Omega Selections
  template <typename CascType, typename T1, typename T2, typename T3, typename C>
  bool passedOmegaSelection(const CascType& casc, const T1& ptrack, const T2& ntrack, const T3& btrack, const C& coll)
  {
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;
    if (!passedSingleTrackSelection(btrack))
      return false;

    // Omega+ Selection (Omega+ -> antiL + K+)
    if (btrack.sign() > 0) {
      if (ntrack.pt() < ptMin_V0_proton)
        return false;
      if (ntrack.pt() > ptMax_V0_proton)
        return false;
      if (ptrack.pt() < ptMin_V0_pion)
        return false;
      if (ptrack.pt() > ptMax_V0_pion)
        return false;

      // PID Selections (TPC)
      if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID Selections (TOF)
      if (requireTOF) {
        if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }
    }

    // Omega- Selection (Omega- -> L + K-)
    if (btrack.sign() < 0) {
      if (ptrack.pt() < ptMin_V0_proton)
        return false;
      if (ptrack.pt() > ptMax_V0_proton)
        return false;
      if (ntrack.pt() < ptMin_V0_pion)
        return false;
      if (ntrack.pt() > ptMax_V0_pion)
        return false;

      // PID Selections (TPC)
      if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID Selections (TOF)
      if (requireTOF) {
        if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }
    }

    // V0 Selections
    if (casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()) < v0cospaMin)
      return false;
    if (casc.v0radius() < minimumV0Radius || casc.v0radius() > maximumV0Radius)
      return false;
    if (casc.dcaV0daughters() > dcaV0DaughtersMax)
      return false;
    if (casc.dcapostopv() < dcapostoPVmin)
      return false;
    if (casc.dcanegtopv() < dcanegtoPVmin)
      return false;

    // Cascade Selections
    if (casc.cascradius() < minimumCascRadius || casc.cascradius() > maximumCascRadius)
      return false;
    if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) < casccospaMin)
      return false;
    if (casc.dcabachtopv() < dcabachtopvMin)
      return false;
    if (casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ()) < dcaV0topvMin)
      return false;
    if (casc.dcacascdaughters() > dcaCascDaughtersMax)
      return false;

    // PID Selection on bachelor
    if (btrack.tpcNSigmaKa() < nsigmaTPCmin || btrack.tpcNSigmaKa() > nsigmaTPCmax)
      return false;

    // PID Selections (TOF)
    if (requireTOF) {
      if (btrack.tofNSigmaKa() < nsigmaTOFmin || btrack.tofNSigmaKa() > nsigmaTOFmax)
        return false;
    }

    // Rapidity Selection
    /*
    TLorentzVector lorentzVect;
    lorentzVect.SetXYZM(casc.px(), casc.py(), casc.pz(), 1.6724);
    if (lorentzVect.Rapidity() < yMin || lorentzVect.Rapidity() > yMax)
      return false;
    */
    return true;
  }

  // Single-Track Selection
  template <typename T1>
  bool passedSingleTrackSelection(const T1& track)
  {
    if (requireITS && (!track.hasITS()))
      return false;
    if (requireITS && track.itsNCls() < minITSnCls)
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsFound() < minTPCnClsFound)
      return false;
    if (track.tpcNClsCrossedRows() < minNCrossedRowsTPC)
      return false;
    if (track.tpcChi2NCl() > maxChi2TPC)
      return false;
    if (track.eta() < etaMin || track.eta() > etaMax)
      return false;
    if (requireTOF && (!track.hasTOF()))
      return false;
    return true;
  }

  void processData(SelectedCollisions::iterator const& collision, aod::AssignedTrackedCascades const& trackedCascades,
                   aod::Cascades const&, FullTracks const& tracks)
  {
    registryQC.fill(HIST("number_of_events_data"), 0.5);
    if (!collision.sel8())
      return;

    registryQC.fill(HIST("number_of_events_data"), 1.5);
    if (abs(collision.posZ()) > zVtx)
      return;

    registryQC.fill(HIST("number_of_events_data"), 2.5);

    for (const auto& trackedCascade : trackedCascades) {

      const auto track = trackedCascade.track_as<FullTracks>();
      const auto& casc = trackedCascade.cascade();
      const auto& btrack = casc.bachelor_as<FullTracks>();
      const auto& ptrack = casc.posTrack_as<FullTracks>();
      const auto& ntrack = casc.negTrack_as<FullTracks>();

      // Calculate Average Cluster Size
      int nITScls(0);
      int clusterSize[7];
      double averageClusterSize(0);
      for (int i = 0; i < 7; i++) {
        clusterSize[i] = (((1 << 4) - 1) & (track.itsClusterSizes() >> 4 * i));
        averageClusterSize += static_cast<double>(clusterSize[i]);

        if (hasHitOnITSlayer(track.itsClusterMap(), i))
          nITScls++;
      }
      averageClusterSize = averageClusterSize / static_cast<double>(nITScls);

      if (passedXiSelection(casc, ptrack, ntrack, btrack, collision) && casc.mXi() > mMin_xi && casc.mXi() < mMax_xi) {
        if (btrack.sign() > 0) {
          registryData.fill(HIST("xi_pos"), averageClusterSize, casc.pt(), casc.eta());
          registryData.fill(HIST("xi_mass_pos"), casc.pt(), casc.mXi());
        }
        if (btrack.sign() < 0) {
          registryData.fill(HIST("xi_neg"), averageClusterSize, casc.pt(), casc.eta());
          registryData.fill(HIST("xi_mass_neg"), casc.pt(), casc.mXi());
        }
      }

      if (passedOmegaSelection(casc, ptrack, ntrack, btrack, collision) && casc.mOmega() > mMin_omega && casc.mOmega() < mMax_omega) {
        if (btrack.sign() > 0) {
          registryData.fill(HIST("omega_pos"), averageClusterSize, casc.pt(), casc.eta());
          registryData.fill(HIST("omega_mass_pos"), casc.pt(), casc.mOmega());
        }
        if (btrack.sign() < 0) {
          registryData.fill(HIST("omega_neg"), averageClusterSize, casc.pt(), casc.eta());
          registryData.fill(HIST("omega_mass_neg"), casc.pt(), casc.mOmega());
        }
      }
    }
  }
  PROCESS_SWITCH(tracked_cascade_properties, processData, "Process data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tracked_cascade_properties>(cfgc)};
}
