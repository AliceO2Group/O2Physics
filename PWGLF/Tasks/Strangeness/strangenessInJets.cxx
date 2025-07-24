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
/// \file strangenessInJets.cxx
///
/// \brief task for analysis of strangeness in jets
/// \author Alberto Calivà (alberto.caliva@cern.ch)
/// \author Francesca Ercolessi (francesca.ercolessi@cern.ch)
/// \author Nicolò Jacazio (nicolo.jacazio@cern.ch)
/// \author Sara Pucillo (sara.pucillo@cern.ch)
///
/// \since May 22, 2024

#include "PWGJE/Core/JetBkgSubUtils.h"
#include "PWGJE/Core/JetDerivedDataUtilities.h"
#include "PWGJE/DataModel/Jet.h"
#include "PWGJE/DataModel/JetReducedData.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"

#include "Common/Core/RecoDecay.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "EventFiltering/Zorro.h"
#include "EventFiltering/ZorroSummary.h"

#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/Track.h"

#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include <TVector2.h>
#include <TVector3.h>

#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/Selector.hh>
#include <fastjet/tools/JetMedianBackgroundEstimator.hh>
#include <fastjet/tools/Subtractor.hh>

#include <cmath>
#include <string>
#include <vector>

using namespace std;
using namespace o2;
using namespace o2::soa;
using namespace o2::aod;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::constants::physics;
using namespace o2::constants::math;
using std::array;

// Define convenient aliases for joined AOD tables
using SelCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms>;
using SimCollisions = soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Ms, aod::McCollisionLabels>;
using DaughterTracks = soa::Join<aod::Tracks, aod::TracksIU, aod::TracksExtra, aod::TracksCovIU, aod::TracksDCA,
                                 aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr,
                                 aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
using DaughterTracksMC = soa::Join<DaughterTracks, aod::McTrackLabels>;

struct StrangenessInJets {

  // Instantiate the CCDB service and API interface
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::ccdb::CcdbApi ccdbApi;

  // Instantiate the Zorro processor for skimmed data and define an output object
  Zorro zorro;
  OutputObj<ZorroSummary> zorroSummary{"zorroSummary"};

  // Define histogram registries
  HistogramRegistry registryData{"registryData", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryMC{"registryMC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry registryQC{"registryQC", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Global analysis parameters
  Configurable<int> particleOfInterest{"particleOfInterest", 0, "0 = K0 and Lambda, 1 = Xi and Omega"};
  Configurable<double> minJetPt{"minJetPt", 10.0, "Minimum reconstructed pt of the jet (GeV/c)"};
  Configurable<double> rJet{"rJet", 0.3, "Jet resolution parameter (R)"};
  Configurable<double> zVtx{"zVtx", 10.0, "Maximum z-vertex position"};
  Configurable<double> deltaEtaEdge{"deltaEtaEdge", 0.05, "eta gap from detector edge"};
  Configurable<bool> cfgSkimmedProcessing{"cfgSkimmedProcessing", false, "Enable processing of skimmed data"};
  Configurable<std::string> triggerName{"triggerName", "fOmega", "Software trigger name"};

  // Track analysis parameters
  Configurable<int> minITSnCls{"minITSnCls", 4, "Minimum number of ITS clusters"};
  Configurable<int> minNCrossedRowsTPC{"minNCrossedRowsTPC", 80, "Minimum number of TPC crossed rows"};
  Configurable<double> maxChi2TPC{"maxChi2TPC", 4.0f, "Maximum chi2 per cluster TPC"};
  Configurable<double> etaMin{"etaMin", -0.8f, "Minimum eta"};
  Configurable<double> etaMax{"etaMax", +0.8f, "Maximum eta"};
  Configurable<double> ptMinV0Proton{"ptMinV0Proton", 0.3f, "Minimum pt of protons from V0"};
  Configurable<double> ptMaxV0Proton{"ptMaxV0Proton", 10.0f, "Maximum pt of protons from V0"};
  Configurable<double> ptMinV0Pion{"ptMinV0Pion", 0.1f, "Minimum pt of pions from V0"};
  Configurable<double> ptMaxV0Pion{"ptMaxV0Pion", 1.5f, "Maximum pt of pions from V0"};
  Configurable<double> ptMinK0Pion{"ptMinK0Pion", 0.3f, "Minimum pt of pions from K0"};
  Configurable<double> ptMaxK0Pion{"ptMaxK0Pion", 10.0f, "Maximum pt of pions from K0"};
  Configurable<double> nsigmaTPCmin{"nsigmaTPCmin", -3.0f, "Minimum nsigma TPC"};
  Configurable<double> nsigmaTPCmax{"nsigmaTPCmax", +3.0f, "Maximum nsigma TPC"};
  Configurable<double> nsigmaTOFmin{"nsigmaTOFmin", -3.0f, "Minimum nsigma TOF"};
  Configurable<double> nsigmaTOFmax{"nsigmaTOFmax", +3.0f, "Maximum nsigma TOF"};
  Configurable<bool> requireITS{"requireITS", false, "Require ITS hit"};
  Configurable<bool> requireTOF{"requireTOF", false, "Require TOF hit"};

  // V0 analysis parameters
  Configurable<double> minimumV0Radius{"minimumV0Radius", 0.5f, "Minimum V0 Radius"};
  Configurable<double> maximumV0Radius{"maximumV0Radius", 40.0f, "Maximum V0 Radius"};
  Configurable<double> dcanegtoPVmin{"dcanegtoPVmin", 0.1f, "Minimum DCA of negative track to primary vertex"};
  Configurable<double> dcapostoPVmin{"dcapostoPVmin", 0.1f, "Minimum DCA of positive track to primary vertex"};
  Configurable<double> v0cospaMin{"v0cospaMin", 0.99f, "Minimum V0 cosine of pointing angle"};
  Configurable<double> dcaV0DaughtersMax{"dcaV0DaughtersMax", 0.5f, "Maximum DCA between V0 daughters"};

  // Cascade analysis parameters
  Configurable<float> minimumCascRadius{"minimumCascRadius", 0.1f, "Minimum cascade radius"};
  Configurable<float> maximumCascRadius{"maximumCascRadius", 40.0f, "Maximum cascade radius"};
  Configurable<float> casccospaMin{"casccospaMin", 0.99f, "Minimum cascade cosine of pointing angle"};
  Configurable<float> dcabachtopvMin{"dcabachtopvMin", 0.1f, "Minimum DCA of bachelor to primary vertex"};
  Configurable<float> dcaV0topvMin{"dcaV0topvMin", 0.1f, "Minimum DCA of V0 to primary vertex"};
  Configurable<float> dcaCascDaughtersMax{"dcaCascDaughtersMax", 0.5f, "Maximum DCA between daughters"};
  Configurable<float> deltaMassXi{"deltaMassXi", 0.02f, "Mass window for Xi rejection"};
  Configurable<float> deltaMassOmega{"deltaMassOmega", 0.02f, "Mass window for Omega rejection"};
  Configurable<float> deltaMassLambda{"deltaMassLambda", 0.02f, "Mass window for Lambda inclusion"};

  // List of Particles
  enum Option { kV0Particles,
                kCascades };

  // Instantiate utility class for jet background subtraction
  JetBkgSubUtils backgroundSub;

  // Initialize CCDB access and histogram registry for Zorro processing
  void initCCDB(aod::BCsWithTimestamps::iterator const& bc)
  {
    if (cfgSkimmedProcessing) {
      zorro.initCCDB(ccdb.service, bc.runNumber(), bc.timestamp(), triggerName.value);
      zorro.populateHistRegistry(registryData, bc.runNumber());
    }
  }

  void init(InitContext const&)
  {
    if (cfgSkimmedProcessing) {
      zorroSummary.setObject(zorro.getZorroSummary());
    }

    // Define binning and axis specifications for multiplicity, eta, pT, PID, and invariant mass histograms
    std::vector<double> multBinning = {0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
    AxisSpec multAxis = {multBinning, "FT0C percentile"};
    const AxisSpec ptAxis{100, 0.0, 10.0, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec invMassK0sAxis{200, 0.44, 0.56, "m_{#pi#pi} (GeV/#it{c}^{2})"};
    const AxisSpec invMassLambdaAxis{200, 1.09, 1.14, "m_{p#pi} (GeV/#it{c}^{2})"};
    const AxisSpec invMassXiAxis{200, 1.28, 1.36, "m_{p#pi#pi} (GeV/#it{c}^{2})"};
    const AxisSpec invMassOmegaAxis{200, 1.63, 1.71, "m_{p#piK} (GeV/#it{c}^{2})"};

    // Histograms for real data
    if (doprocessData) {

      // Event counters
      registryData.add("number_of_events_data", "number of events in data", HistType::kTH1D, {{20, 0, 20, "Event Cuts"}});
      registryData.add("number_of_events_vsmultiplicity", "number of events in data vs multiplicity", HistType::kTH1D, {{101, 0, 101, "Multiplicity percentile"}});

      // Histograms for analysis of strange hadrons
      switch (particleOfInterest) {
        case kV0Particles:
          registryData.add("Lambda_in_jet", "Lambda_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassLambdaAxis});
          registryData.add("AntiLambda_in_jet", "AntiLambda_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassLambdaAxis});
          registryData.add("Lambda_in_ue", "Lambda_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassLambdaAxis});
          registryData.add("AntiLambda_in_ue", "AntiLambda_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassLambdaAxis});
          registryData.add("K0s_in_jet", "K0s_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassK0sAxis});
          registryData.add("K0s_in_ue", "K0s_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassK0sAxis});
          break;
        case kCascades:
          registryData.add("XiPos_in_jet", "XiPos_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassXiAxis});
          registryData.add("XiPos_in_ue", "XiPos_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassXiAxis});
          registryData.add("XiNeg_in_jet", "XiNeg_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassXiAxis});
          registryData.add("XiNeg_in_ue", "XiNeg_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassXiAxis});
          registryData.add("OmegaPos_in_jet", "OmegaPos_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassOmegaAxis});
          registryData.add("OmegaPos_in_ue", "OmegaPos_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassOmegaAxis});
          registryData.add("OmegaNeg_in_jet", "OmegaNeg_in_jet", HistType::kTH3F, {multBinning, ptAxis, invMassOmegaAxis});
          registryData.add("OmegaNeg_in_ue", "OmegaNeg_in_ue", HistType::kTH3F, {multBinning, ptAxis, invMassOmegaAxis});
          break;
        default:
          LOG(fatal) << "Cannot interpret particle " << particleOfInterest;
          break;
      }
    }

    // Histograms for mc generated
    if (doprocessMCgenerated) {

      // Event counter
      registryMC.add("number_of_events_mc_gen", "number of gen events in mc", HistType::kTH1D, {{10, 0, 10, "Event Cuts"}});

      // Histograms for analysis
      registryMC.add("K0s_generated_jet", "K0s_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("K0s_generated_ue", "K0s_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("Lambda_generated_jet", "Lambda_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("Lambda_generated_ue", "Lambda_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("AntiLambda_generated_jet", "AntiLambda_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("AntiLambda_generated_ue", "AntiLambda_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiPos_generated_jet", "XiPos_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiPos_generated_ue", "XiPos_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiNeg_generated_jet", "XiNeg_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiNeg_generated_ue", "XiNeg_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaPos_generated_jet", "OmegaPos_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaPos_generated_ue", "OmegaPos_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaNeg_generated_jet", "OmegaNeg_generated_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaNeg_generated_ue", "OmegaNeg_generated_ue", HistType::kTH2F, {multBinning, ptAxis});
    }

    // Histograms for mc reconstructed
    if (doprocessMCreconstructed) {

      // Event counter
      registryMC.add("number_of_events_mc_rec", "number of rec events in mc", HistType::kTH1D, {{10, 0, 10, "Event Cuts"}});

      // Histograms for analysis
      registryMC.add("K0s_reconstructed_jet", "K0s_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("K0s_reconstructed_ue", "K0s_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("Lambda_reconstructed_jet", "Lambda_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("Lambda_reconstructed_ue", "Lambda_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("AntiLambda_reconstructed_jet", "AntiLambda_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("AntiLambda_reconstructed_ue", "AntiLambda_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiPos_reconstructed_jet", "XiPos_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiPos_reconstructed_ue", "XiPos_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiNeg_reconstructed_jet", "XiNeg_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("XiNeg_reconstructed_ue", "XiNeg_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaPos_reconstructed_jet", "OmegaPos_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaPos_reconstructed_ue", "OmegaPos_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaNeg_reconstructed_jet", "OmegaNeg_reconstructed_jet", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("OmegaNeg_reconstructed_ue", "OmegaNeg_reconstructed_ue", HistType::kTH2F, {multBinning, ptAxis});

      // Histograms for secondary hadrons
      registryMC.add("K0s_reconstructed_jet_incl", "K0s_reconstructed_jet_incl", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("K0s_reconstructed_ue_incl", "K0s_reconstructed_ue_incl", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("Lambda_reconstructed_jet_incl", "Lambda_reconstructed_jet_incl", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("Lambda_reconstructed_ue_incl", "Lambda_reconstructed_ue_incl", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("AntiLambda_reconstructed_jet_incl", "AntiLambda_reconstructed_jet_incl", HistType::kTH2F, {multBinning, ptAxis});
      registryMC.add("AntiLambda_reconstructed_ue_incl", "AntiLambda_reconstructed_ue_incl", HistType::kTH2F, {multBinning, ptAxis});
    }
  }

  // Calculation of perpendicular axes
  void getPerpendicularAxis(TVector3 p, TVector3& u, double sign)
  {
    // initialization
    double ux(0), uy(0), uz(0);

    // components of vector p
    double px = p.X();
    double py = p.Y();
    double pz = p.Z();

    // protection 1
    if (px == 0 && py != 0) {
      uy = -(pz * pz) / py;
      ux = sign * std::sqrt(py * py - (pz * pz * pz * pz) / (py * py));
      uz = pz;
      u.SetXYZ(ux, uy, uz);
      return;
    }

    // protection 2
    if (py == 0 && px != 0) {
      ux = -(pz * pz) / px;
      uy = sign * std::sqrt(px * px - (pz * pz * pz * pz) / (px * px));
      uz = pz;
      u.SetXYZ(ux, uy, uz);
      return;
    }

    // equation parameters
    double a = px * px + py * py;
    double b = 2.0 * px * pz * pz;
    double c = pz * pz * pz * pz - py * py * py * py - px * px * py * py;
    double delta = b * b - 4.0 * a * c;

    // protection agains delta<0
    if (delta < 0) {
      return;
    }

    // solutions
    ux = (-b + sign * std::sqrt(delta)) / (2.0 * a);
    uy = (-pz * pz - px * ux) / py;
    uz = pz;
    u.SetXYZ(ux, uy, uz);
    return;
  }

  // Delta phi calculation
  double getDeltaPhi(double a1, double a2)
  {
    double deltaPhi(0);
    double phi1 = TVector2::Phi_0_2pi(a1);
    double phi2 = TVector2::Phi_0_2pi(a2);
    double diff = std::fabs(phi1 - phi2);

    if (diff <= PI)
      deltaPhi = diff;
    if (diff > PI)
      deltaPhi = TwoPI - diff;

    return deltaPhi;
  }

  // Find ITS hit
  template <typename TrackIts>
  bool hasITSHitOnLayer(const TrackIts& track, int layer)
  {
    int ibit = layer - 1;
    return (track.itsClusterMap() & (1 << ibit));
  }

  // Single-track selection for particles inside jets
  template <typename JetTrack>
  bool passedTrackSelectionForJetReconstruction(const JetTrack& track)
  {
    const int minTpcCr = 70;
    const double maxChi2Tpc = 4.0;
    const double maxChi2Its = 36.0;
    const double maxPseudorapidity = 0.8;
    const double minPtTrack = 0.1;
    const double dcaxyMaxTrackPar0 = 0.0105;
    const double dcaxyMaxTrackPar1 = 0.035;
    const double dcaxyMaxTrackPar2 = 1.1;
    const double dcazMaxTrack = 2.0;

    if (!track.hasITS())
      return false;
    if ((!hasITSHitOnLayer(track, 1)) && (!hasITSHitOnLayer(track, 2)) && (!hasITSHitOnLayer(track, 3)))
      return false;
    if (!track.hasTPC())
      return false;
    if (track.tpcNClsCrossedRows() < minTpcCr)
      return false;
    if (track.tpcChi2NCl() > maxChi2Tpc)
      return false;
    if (track.itsChi2NCl() > maxChi2Its)
      return false;
    if (std::fabs(track.eta()) > maxPseudorapidity)
      return false;
    if (track.pt() < minPtTrack)
      return false;
    if (std::fabs(track.dcaXY()) > (dcaxyMaxTrackPar0 + dcaxyMaxTrackPar1 / std::pow(track.pt(), dcaxyMaxTrackPar2)))
      return false;
    if (std::fabs(track.dcaZ()) > dcazMaxTrack)
      return false;
    return true;
  }

  // Lambda selections
  template <typename Lambda, typename TrackPos, typename TrackNeg>
  bool passedLambdaSelection(const Lambda& v0, const TrackPos& ptrack, const TrackNeg& ntrack)
  {
    // Single-track selections
    if (!passedSingleTrackSelection(ptrack) || !passedSingleTrackSelection(ntrack))
      return false;

    // Momentum of lambda daughters
    TVector3 proton(v0.pxpos(), v0.pypos(), v0.pzpos());
    TVector3 pion(v0.pxneg(), v0.pyneg(), v0.pzneg());

    // Selection on pt of Lambda daughters
    if (proton.Pt() < ptMinV0Proton || proton.Pt() > ptMaxV0Proton)
      return false;
    if (pion.Pt() < ptMinV0Pion || pion.Pt() > ptMaxV0Pion)
      return false;

    // V0 selections
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (std::fabs(v0.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;

    // PID selections (TPC): positive track = proton, negative track = pion
    if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
      return false;
    if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // PID selections (TOF): positive track = proton, negative track = pion
    if (requireTOF) {
      if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }
    return true;
  }

  // AntiLambda selections
  template <typename AntiLambda, typename TrackPos, typename TrackNeg>
  bool passedAntiLambdaSelection(const AntiLambda& v0, const TrackPos& ptrack, const TrackNeg& ntrack)
  {
    // Single-track selections
    if (!passedSingleTrackSelection(ptrack) || !passedSingleTrackSelection(ntrack))
      return false;

    // Momentum AntiLambda daughters
    TVector3 pion(v0.pxpos(), v0.pypos(), v0.pzpos());
    TVector3 proton(v0.pxneg(), v0.pyneg(), v0.pzneg());

    // Selections on pt of Antilambda daughters
    if (proton.Pt() < ptMinV0Proton || proton.Pt() > ptMaxV0Proton)
      return false;
    if (pion.Pt() < ptMinV0Pion || pion.Pt() > ptMaxV0Pion)
      return false;

    // V0 selections
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (std::fabs(v0.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;

    // PID selections (TPC): negative track = proton, positive track = pion
    if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;
    if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
      return false;

    // PID selections (TOF): negative track = proton, positive track = pion
    if (requireTOF) {
      if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
        return false;
    }
    return true;
  }

  // K0s selections
  template <typename K0short, typename TrackPos, typename TrackNeg>
  bool passedK0ShortSelection(const K0short& v0, const TrackPos& ptrack, const TrackNeg& ntrack)
  {
    // Single-Track Selections
    if (!passedSingleTrackSelection(ptrack) || !passedSingleTrackSelection(ntrack))
      return false;

    // Momentum of K0s daughters
    TVector3 pionPos(v0.pxpos(), v0.pypos(), v0.pzpos());
    TVector3 pionNeg(v0.pxneg(), v0.pyneg(), v0.pzneg());

    // Selections on pt of K0s daughters
    if (pionPos.Pt() < ptMinK0Pion || pionPos.Pt() > ptMaxK0Pion)
      return false;
    if (pionNeg.Pt() < ptMinK0Pion || pionNeg.Pt() > ptMaxK0Pion)
      return false;

    // V0 selections
    if (v0.v0cosPA() < v0cospaMin)
      return false;
    if (v0.v0radius() < minimumV0Radius || v0.v0radius() > maximumV0Radius)
      return false;
    if (std::fabs(v0.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    if (std::fabs(v0.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(v0.dcanegtopv()) < dcanegtoPVmin)
      return false;

    // PID selections (TPC)
    if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;
    if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // PID selections (TOF)
    if (requireTOF) {
      if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
      if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }
    return true;
  }

  // Xi Selections
  template <typename Xi, typename TrackPos, typename TrackNeg, typename TrackBac, typename Coll>
  bool passedXiSelection(const Xi& casc, const TrackPos& ptrack, const TrackNeg& ntrack, const TrackBac& btrack, const Coll& coll)
  {
    // Single-track selections on cascade daughters
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;
    if (!passedSingleTrackSelection(btrack))
      return false;

    // Xi+ selection (Xi+ -> antiL + pi+)
    if (btrack.sign() > 0) {
      if (ntrack.pt() < ptMinV0Proton || ntrack.pt() > ptMaxV0Proton)
        return false;
      if (ptrack.pt() < ptMinV0Pion || ptrack.pt() > ptMaxV0Pion)
        return false;

      // PID selections (TPC)
      if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID selections (TOF)
      if (requireTOF) {
        if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }

      // Require that V0 is compatible with Lambda
      ROOT::Math::PxPyPzMVector pProton;
      ROOT::Math::PxPyPzMVector pPion;
      pProton.SetCoordinates(ntrack.px(), ntrack.py(), ntrack.pz(), MassProton);
      pPion.SetCoordinates(ptrack.px(), ptrack.py(), ptrack.pz(), MassPionCharged);
      double mLambda = (pProton + pPion).M();
      if (std::fabs(mLambda - MassLambda0) > deltaMassLambda)
        return false;
    }

    // Xi- selection (Xi- -> L + pi-)
    if (btrack.sign() < 0) {
      if (ptrack.pt() < ptMinV0Proton || ptrack.pt() > ptMaxV0Proton)
        return false;
      if (ntrack.pt() < ptMinV0Pion || ntrack.pt() > ptMaxV0Pion)
        return false;

      // PID selections (TPC)
      if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID selections (TOF)
      if (requireTOF) {
        if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }

      // Require that V0 is compatible with Lambda
      ROOT::Math::PxPyPzMVector pProton;
      ROOT::Math::PxPyPzMVector pPion;
      pProton.SetCoordinates(ptrack.px(), ptrack.py(), ptrack.pz(), MassProton);
      pPion.SetCoordinates(ntrack.px(), ntrack.py(), ntrack.pz(), MassPionCharged);
      double mLambda = (pProton + pPion).M();
      if (std::fabs(mLambda - MassLambda0) > deltaMassLambda)
        return false;
    }

    // V0 selections
    if (casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()) < v0cospaMin)
      return false;
    if (casc.v0radius() < minimumV0Radius || casc.v0radius() > maximumV0Radius)
      return false;
    if (std::fabs(casc.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    if (std::fabs(casc.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(casc.dcanegtopv()) < dcanegtoPVmin)
      return false;

    // Cascade selections
    if (casc.cascradius() < minimumCascRadius || casc.cascradius() > maximumCascRadius)
      return false;
    if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) < casccospaMin)
      return false;
    if (std::fabs(casc.dcabachtopv()) < dcabachtopvMin)
      return false;
    if (std::fabs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ())) < dcaV0topvMin)
      return false;
    if (std::fabs(casc.dcacascdaughters()) > dcaCascDaughtersMax)
      return false;

    // PID selection on bachelor
    if (btrack.tpcNSigmaPi() < nsigmaTPCmin || btrack.tpcNSigmaPi() > nsigmaTPCmax)
      return false;

    // PID selections (TOF)
    if (requireTOF) {
      if (btrack.tofNSigmaPi() < nsigmaTOFmin || btrack.tofNSigmaPi() > nsigmaTOFmax)
        return false;
    }

    // Reject candidates compatible with Omega
    if (std::fabs(casc.mOmega() - MassOmegaMinus) < deltaMassOmega)
      return false;
    return true;
  }

  // Omega selections
  template <typename Omega, typename TrackPos, typename TrackNeg, typename TrackBac, typename Coll>
  bool passedOmegaSelection(const Omega& casc, const TrackPos& ptrack, const TrackNeg& ntrack, const TrackBac& btrack, const Coll& coll)
  {
    // Single-track selections on cascade daughters
    if (!passedSingleTrackSelection(ptrack))
      return false;
    if (!passedSingleTrackSelection(ntrack))
      return false;
    if (!passedSingleTrackSelection(btrack))
      return false;

    // Omega+ selection (Omega+ -> antiL + K+)
    if (btrack.sign() > 0) {
      if (ntrack.pt() < ptMinV0Proton || ntrack.pt() > ptMaxV0Proton)
        return false;
      if (ptrack.pt() < ptMinV0Pion || ptrack.pt() > ptMaxV0Pion)
        return false;

      // PID selections (TPC)
      if (ntrack.tpcNSigmaPr() < nsigmaTPCmin || ntrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ptrack.tpcNSigmaPi() < nsigmaTPCmin || ptrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID selections (TOF)
      if (requireTOF) {
        if (ntrack.tofNSigmaPr() < nsigmaTOFmin || ntrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ptrack.tofNSigmaPi() < nsigmaTOFmin || ptrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }

      // Require that V0 is compatible with Lambda
      ROOT::Math::PxPyPzMVector pProton;
      ROOT::Math::PxPyPzMVector pPion;
      pProton.SetCoordinates(ntrack.px(), ntrack.py(), ntrack.pz(), MassProton);
      pPion.SetCoordinates(ptrack.px(), ptrack.py(), ptrack.pz(), MassPionCharged);
      double mLambda = (pProton + pPion).M();
      if (std::fabs(mLambda - MassLambda0) > deltaMassLambda)
        return false;
    }

    // Omega- selection (Omega- -> L + K-)
    if (btrack.sign() < 0) {
      if (ptrack.pt() < ptMinV0Proton || ptrack.pt() > ptMaxV0Proton)
        return false;
      if (ntrack.pt() < ptMinV0Pion || ntrack.pt() > ptMaxV0Pion)
        return false;

      // PID selections (TPC)
      if (ptrack.tpcNSigmaPr() < nsigmaTPCmin || ptrack.tpcNSigmaPr() > nsigmaTPCmax)
        return false;
      if (ntrack.tpcNSigmaPi() < nsigmaTPCmin || ntrack.tpcNSigmaPi() > nsigmaTPCmax)
        return false;

      // PID selections (TOF)
      if (requireTOF) {
        if (ptrack.tofNSigmaPr() < nsigmaTOFmin || ptrack.tofNSigmaPr() > nsigmaTOFmax)
          return false;
        if (ntrack.tofNSigmaPi() < nsigmaTOFmin || ntrack.tofNSigmaPi() > nsigmaTOFmax)
          return false;
      }

      // Require that V0 is compatible with Lambda
      ROOT::Math::PxPyPzMVector pProton;
      ROOT::Math::PxPyPzMVector pPion;
      pProton.SetCoordinates(ptrack.px(), ptrack.py(), ptrack.pz(), MassProton);
      pPion.SetCoordinates(ntrack.px(), ntrack.py(), ntrack.pz(), MassPionCharged);
      double mLambda = (pProton + pPion).M();
      if (std::fabs(mLambda - MassLambda0) > deltaMassLambda)
        return false;
    }

    // V0 selections
    if (casc.v0cosPA(coll.posX(), coll.posY(), coll.posZ()) < v0cospaMin)
      return false;
    if (casc.v0radius() < minimumV0Radius || casc.v0radius() > maximumV0Radius)
      return false;
    if (std::fabs(casc.dcaV0daughters()) > dcaV0DaughtersMax)
      return false;
    if (std::fabs(casc.dcapostopv()) < dcapostoPVmin)
      return false;
    if (std::fabs(casc.dcanegtopv()) < dcanegtoPVmin)
      return false;

    // Cascade selections
    if (casc.cascradius() < minimumCascRadius || casc.cascradius() > maximumCascRadius)
      return false;
    if (casc.casccosPA(coll.posX(), coll.posY(), coll.posZ()) < casccospaMin)
      return false;
    if (std::fabs(casc.dcabachtopv()) < dcabachtopvMin)
      return false;
    if (std::fabs(casc.dcav0topv(coll.posX(), coll.posY(), coll.posZ())) < dcaV0topvMin)
      return false;
    if (std::fabs(casc.dcacascdaughters()) > dcaCascDaughtersMax)
      return false;

    // PID selection on bachelor
    if (btrack.tpcNSigmaKa() < nsigmaTPCmin || btrack.tpcNSigmaKa() > nsigmaTPCmax)
      return false;

    // PID selections (TOF)
    if (requireTOF) {
      if (btrack.tofNSigmaKa() < nsigmaTOFmin || btrack.tofNSigmaKa() > nsigmaTOFmax)
        return false;
    }

    // Reject candidates compatible with Xi
    if (std::fabs(casc.mXi() - MassXiMinus) < deltaMassXi)
      return false;
    return true;
  }

  // Single-track selection
  template <typename Track>
  bool passedSingleTrackSelection(const Track& track)
  {
    if (requireITS && (!track.hasITS()))
      return false;
    if (requireITS && track.itsNCls() < minITSnCls)
      return false;
    if (!track.hasTPC())
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

  // Process data
  void processData(SelCollisions::iterator const& collision, aod::V0Datas const& fullV0s,
                   aod::CascDataExt const& Cascades, DaughterTracks const& tracks,
                   aod::BCsWithTimestamps const&)
  {
    // Fill event counter before event selection
    registryData.fill(HIST("number_of_events_data"), 0.5);

    // Get the bunch crossing (BC) information associated with the collision
    auto bc = collision.template bc_as<aod::BCsWithTimestamps>();

    // Initialize CCDB objects using the BC info
    initCCDB(bc);

    // If skimmed processing is enabled, skip this event unless it passes Zorro selection
    if (cfgSkimmedProcessing && !zorro.isSelected(collision.template bc_as<aod::BCsWithTimestamps>().globalBC())) {
      return;
    }

    // Fill event counter after zorro selection
    registryData.fill(HIST("number_of_events_data"), 1.5);

    // Event selection
    if (!collision.sel8() || std::fabs(collision.posZ()) > zVtx)
      return;

    // Fill event counter after event selection
    registryData.fill(HIST("number_of_events_data"), 2.5);

    // Loop over reconstructed tracks
    std::vector<fastjet::PseudoJet> fjParticles;
    for (auto const& track : tracks) {

      // Require that tracks pass selection criteria
      if (!passedTrackSelectionForJetReconstruction(track))
        continue;

      // 4-momentum representation of a particle
      fastjet::PseudoJet fourMomentum(track.px(), track.py(), track.pz(), track.energy(MassPionCharged));
      fjParticles.emplace_back(fourMomentum);
    }

    // Reject empty events
    if (fjParticles.size() < 1)
      return;
    registryData.fill(HIST("number_of_events_data"), 3.5);

    // Cluster particles using the anti-kt algorithm
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
    fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
    fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
    auto [rhoPerp, rhoMPerp] = backgroundSub.estimateRhoPerpCone(fjParticles, jets);

    // Jet selection
    bool isAtLeastOneJetSelected = false;
    std::vector<TVector3> selectedJet;
    std::vector<TVector3> ue1;
    std::vector<TVector3> ue2;

    // Loop over reconstructed jets
    for (const auto& jet : jets) {

      // Jet must be fully contained in the acceptance
      if ((std::fabs(jet.eta()) + rJet) > (etaMax - deltaEtaEdge))
        continue;

      // Jet pt must be larger than threshold
      auto jetForSub = jet;
      fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
      if (jetMinusBkg.pt() < minJetPt)
        continue;
      isAtLeastOneJetSelected = true;

      // Calculation of perpendicular cones
      TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
      TVector3 ueAxis1(0, 0, 0);
      TVector3 ueAxis2(0, 0, 0);
      getPerpendicularAxis(jetAxis, ueAxis1, +1);
      getPerpendicularAxis(jetAxis, ueAxis2, -1);

      // Store jet and UE axes
      selectedJet.emplace_back(jetAxis);
      ue1.emplace_back(ueAxis1);
      ue2.emplace_back(ueAxis2);
    }
    if (!isAtLeastOneJetSelected)
      return;

    // Fill event counter with events with at least one jet
    registryData.fill(HIST("number_of_events_data"), 4.5);

    // Event multiplicity
    const float multiplicity = collision.centFT0M();

    // Fill event multiplicity
    registryData.fill(HIST("number_of_events_vsmultiplicity"), multiplicity);

    // Loop over selected jets
    for (int i = 0; i < static_cast<int>(selectedJet.size()); i++) {

      // kV0Particles
      if (particleOfInterest == Option::kV0Particles) {
        for (const auto& v0 : fullV0s) {

          // Get V0 daughters
          const auto& pos = v0.posTrack_as<DaughterTracks>();
          const auto& neg = v0.negTrack_as<DaughterTracks>();
          TVector3 v0dir(v0.px(), v0.py(), v0.pz());

          // Calculate distance from jet and UE axes
          float deltaEtaJet = v0dir.Eta() - selectedJet[i].Eta();
          float deltaPhiJet = getDeltaPhi(v0dir.Phi(), selectedJet[i].Phi());
          float deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
          float deltaEtaUe1 = v0dir.Eta() - ue1[i].Eta();
          float deltaPhiUe1 = getDeltaPhi(v0dir.Phi(), ue1[i].Phi());
          float deltaRue1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
          float deltaEtaUe2 = v0dir.Eta() - ue2[i].Eta();
          float deltaPhiUe2 = getDeltaPhi(v0dir.Phi(), ue2[i].Phi());
          float deltaRue2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

          // K0s
          if (passedK0ShortSelection(v0, pos, neg)) {
            if (deltaRjet < rJet) {
              registryData.fill(HIST("K0s_in_jet"), multiplicity, v0.pt(), v0.mK0Short());
            }
            if (deltaRue1 < rJet || deltaRue2 < rJet) {
              registryData.fill(HIST("K0s_in_ue"), multiplicity, v0.pt(), v0.mK0Short());
            }
          }
          // Lambda
          if (passedLambdaSelection(v0, pos, neg)) {
            if (deltaRjet < rJet) {
              registryData.fill(HIST("Lambda_in_jet"), multiplicity, v0.pt(), v0.mLambda());
            }
            if (deltaRue1 < rJet || deltaRue2 < rJet) {
              registryData.fill(HIST("Lambda_in_ue"), multiplicity, v0.pt(), v0.mLambda());
            }
          }
          // AntiLambda
          if (passedAntiLambdaSelection(v0, pos, neg)) {
            if (deltaRjet < rJet) {
              registryData.fill(HIST("AntiLambda_in_jet"), multiplicity, v0.pt(), v0.mAntiLambda());
            }
            if (deltaRue1 < rJet || deltaRue2 < rJet) {
              registryData.fill(HIST("AntiLambda_in_ue"), multiplicity, v0.pt(), v0.mAntiLambda());
            }
          }
        }
      }

      // Cascades
      if (particleOfInterest == Option::kCascades) {
        for (const auto& casc : Cascades) {

          // Get cascade daughters
          auto bach = casc.bachelor_as<DaughterTracks>();
          auto pos = casc.posTrack_as<DaughterTracks>();
          auto neg = casc.negTrack_as<DaughterTracks>();
          TVector3 cascadeDir(casc.px(), casc.py(), casc.pz());

          // Calculate distance from jet and UE axes
          double deltaEtaJet = cascadeDir.Eta() - selectedJet[i].Eta();
          double deltaPhiJet = getDeltaPhi(cascadeDir.Phi(), selectedJet[i].Phi());
          double deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
          double deltaEtaUe1 = cascadeDir.Eta() - ue1[i].Eta();
          double deltaPhiUe1 = getDeltaPhi(cascadeDir.Phi(), ue1[i].Phi());
          double deltaRue1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
          double deltaEtaUe2 = cascadeDir.Eta() - ue2[i].Eta();
          double deltaPhiUe2 = getDeltaPhi(cascadeDir.Phi(), ue2[i].Phi());
          double deltaRue2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

          // Xi+
          if (passedXiSelection(casc, pos, neg, bach, collision) && bach.sign() > 0) {
            if (deltaRjet < rJet) {
              registryData.fill(HIST("XiPos_in_jet"), multiplicity, casc.pt(), casc.mXi());
            }
            if (deltaRue1 < rJet || deltaRue2 < rJet) {
              registryData.fill(HIST("XiPos_in_ue"), multiplicity, casc.pt(), casc.mXi());
            }
          }
          // Xi-
          if (passedXiSelection(casc, pos, neg, bach, collision) && bach.sign() < 0) {
            if (deltaRjet < rJet) {
              registryData.fill(HIST("XiNeg_in_jet"), multiplicity, casc.pt(), casc.mXi());
            }
            if (deltaRue1 < rJet || deltaRue2 < rJet) {
              registryData.fill(HIST("XiNeg_in_ue"), multiplicity, casc.pt(), casc.mXi());
            }
          }
          // Omega+
          if (passedOmegaSelection(casc, pos, neg, bach, collision) && bach.sign() > 0) {
            if (deltaRjet < rJet) {
              registryData.fill(HIST("OmegaPos_in_jet"), multiplicity, casc.pt(), casc.mOmega());
            }
            if (deltaRue1 < rJet || deltaRue2 < rJet) {
              registryData.fill(HIST("OmegaPos_in_ue"), multiplicity, casc.pt(), casc.mOmega());
            }
          }
          // Omega-
          if (passedOmegaSelection(casc, pos, neg, bach, collision) && bach.sign() < 0) {
            if (deltaRjet < rJet) {
              registryData.fill(HIST("OmegaNeg_in_jet"), multiplicity, casc.pt(), casc.mOmega());
            }
            if (deltaRue1 < rJet || deltaRue2 < rJet) {
              registryData.fill(HIST("OmegaNeg_in_ue"), multiplicity, casc.pt(), casc.mOmega());
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(StrangenessInJets, processData, "Process data", true);

  Preslice<aod::V0Datas> perCollisionV0 = o2::aod::v0data::collisionId;
  Preslice<aod::CascDataExt> perCollisionCasc = o2::aod::cascade::collisionId;
  Preslice<aod::McParticles> perMCCollision = o2::aod::mcparticle::mcCollisionId;
  Preslice<DaughterTracksMC> perCollisionTrk = o2::aod::track::collisionId;

  // Generated MC events
  void processMCgenerated(aod::McCollisions const& collisions, aod::McParticles const& mcParticles)
  {
    // Loop over all simulated collision events
    for (const auto& collision : collisions) {

      // Fill event counter before any selection
      registryMC.fill(HIST("number_of_events_mc_gen"), 0.5);

      // Need to apply event selection to simulated events
      registryMC.fill(HIST("number_of_events_mc_gen"), 1.5);

      // Require vertex position within the allowed z range
      if (std::fabs(collision.posZ()) > zVtx)
        continue;

      // Fill event counter after selection on z-vertex
      registryMC.fill(HIST("number_of_events_mc_gen"), 2.5);

      // Multiplicity of generated event
      double genMultiplicity = 0.0;

      // MC particles per collision
      auto mcParticlesPerColl = mcParticles.sliceBy(perMCCollision, collision.globalIndex());

      // Loop over all MC particles and select physical primaries within acceptance
      std::vector<fastjet::PseudoJet> fjParticles;
      for (const auto& particle : mcParticlesPerColl) {
        if (!particle.isPhysicalPrimary())
          continue;
        double minPtParticle = 0.1;
        if (particle.eta() < etaMin || particle.eta() > etaMax || particle.pt() < minPtParticle)
          continue;

        // Build 4-momentum assuming charged pion mass
        double energy = std::sqrt(particle.p() * particle.p() + MassPionCharged * MassPionCharged);
        fastjet::PseudoJet fourMomentum(particle.px(), particle.py(), particle.pz(), energy);
        fourMomentum.set_user_index(particle.pdgCode());
        fjParticles.emplace_back(fourMomentum);
      }

      // Skip events with no particles
      if (fjParticles.size() < 1)
        continue;
      registryMC.fill(HIST("number_of_events_mc_gen"), 3.5);

      // Cluster MC particles into jets using anti-kt algorithm
      fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
      fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
      fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());

      // Estimate background energy density (rho) in perpendicular cone
      auto [rhoPerp, rhoMPerp] = backgroundSub.estimateRhoPerpCone(fjParticles, jets);

      // Loop over clustered jets
      for (const auto& jet : jets) {

        // Jet must be fully contained in acceptance
        if ((std::fabs(jet.eta()) + rJet) > (etaMax - deltaEtaEdge))
          continue;

        // Subtract background energy from jet
        auto jetForSub = jet;
        fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);

        // Apply jet pT threshold
        if (jetMinusBkg.pt() < minJetPt)
          continue;
        registryMC.fill(HIST("number_of_events_mc_gen"), 4.5);

        // Set up two perpendicular cone axes for underlying event estimation
        TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
        double coneRadius = std::sqrt(jet.area() / PI);
        TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
        getPerpendicularAxis(jetAxis, ueAxis1, +1);
        getPerpendicularAxis(jetAxis, ueAxis2, -1);

        // Loop over MC particles
        for (const auto& particle : mcParticlesPerColl) {
          if (!particle.isPhysicalPrimary())
            continue;
          double minPtParticle = 0.1;
          if (particle.eta() < etaMin || particle.eta() > etaMax || particle.pt() < minPtParticle)
            continue;

          // Compute distance of particles from jet and UE axes
          double deltaEtaJet = particle.eta() - jetAxis.Eta();
          double deltaPhiJet = getDeltaPhi(particle.phi(), jetAxis.Phi());
          double deltaRJet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
          double deltaEtaUe1 = particle.eta() - ueAxis1.Eta();
          double deltaPhiUe1 = getDeltaPhi(particle.phi(), ueAxis1.Phi());
          double deltaRUe1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
          double deltaEtaUe2 = particle.eta() - ueAxis2.Eta();
          double deltaPhiUe2 = getDeltaPhi(particle.phi(), ueAxis2.Phi());
          double deltaRUe2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

          // Select particles inside jet
          if (deltaRJet < coneRadius) {
            switch (particle.pdgCode()) {
              case kK0Short:
                registryMC.fill(HIST("K0s_generated_jet"), genMultiplicity, particle.pt());
                break;
              case kLambda0:
                registryMC.fill(HIST("Lambda_generated_jet"), genMultiplicity, particle.pt());
                break;
              case kLambda0Bar:
                registryMC.fill(HIST("AntiLambda_generated_jet"), genMultiplicity, particle.pt());
                break;
              case kXiMinus:
                registryMC.fill(HIST("XiNeg_generated_jet"), genMultiplicity, particle.pt());
                break;
              case kXiPlusBar:
                registryMC.fill(HIST("XiPos_generated_jet"), genMultiplicity, particle.pt());
                break;
              case kOmegaMinus:
                registryMC.fill(HIST("OmegaNeg_generated_jet"), genMultiplicity, particle.pt());
                break;
              case kOmegaPlusBar:
                registryMC.fill(HIST("OmegaPos_generated_jet"), genMultiplicity, particle.pt());
                break;
              default:
                break;
            }
          }

          // Select particles inside UE cones
          if (deltaRUe1 < coneRadius || deltaRUe2 < coneRadius) {
            switch (particle.pdgCode()) {
              case kK0Short:
                registryMC.fill(HIST("K0s_generated_ue"), genMultiplicity, particle.pt());
                break;
              case kLambda0:
                registryMC.fill(HIST("Lambda_generated_ue"), genMultiplicity, particle.pt());
                break;
              case kLambda0Bar:
                registryMC.fill(HIST("AntiLambda_generated_ue"), genMultiplicity, particle.pt());
                break;
              case kXiMinus:
                registryMC.fill(HIST("XiNeg_generated_ue"), genMultiplicity, particle.pt());
                break;
              case kXiPlusBar:
                registryMC.fill(HIST("XiPos_generated_ue"), genMultiplicity, particle.pt());
                break;
              case kOmegaMinus:
                registryMC.fill(HIST("OmegaNeg_generated_ue"), genMultiplicity, particle.pt());
                break;
              case kOmegaPlusBar:
                registryMC.fill(HIST("OmegaPos_generated_ue"), genMultiplicity, particle.pt());
                break;
              default:
                break;
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(StrangenessInJets, processMCgenerated, "process generated events", false);

  // Reconstructed MC events
  void processMCreconstructed(SimCollisions const& collisions, DaughterTracksMC const& mcTracks,
                              aod::V0Datas const& fullV0s, aod::CascDataExt const& Cascades,
                              const aod::McParticles&)
  {
    for (const auto& collision : collisions) {

      // Fill event counter before any selection
      registryMC.fill(HIST("number_of_events_mc_rec"), 0.5);
      if (!collision.sel8())
        continue;

      // Fill event counter after event selection
      registryMC.fill(HIST("number_of_events_mc_rec"), 1.5);
      if (std::fabs(collision.posZ()) > zVtx)
        continue;

      // Fill event counter after selection on z-vertex
      registryMC.fill(HIST("number_of_events_mc_rec"), 2.5);

      // Event multiplicity
      const float multiplicity = collision.centFT0M();

      // Number of V0 and cascades per collision
      auto v0sPerColl = fullV0s.sliceBy(perCollisionV0, collision.globalIndex());
      auto cascPerColl = Cascades.sliceBy(perCollisionCasc, collision.globalIndex());
      auto tracksPerColl = mcTracks.sliceBy(perCollisionTrk, collision.globalIndex());

      // Loop over reconstructed tracks
      std::vector<fastjet::PseudoJet> fjParticles;
      for (auto const& track : tracksPerColl) {
        if (!passedTrackSelectionForJetReconstruction(track))
          continue;

        // 4-momentum representation of a particle
        fastjet::PseudoJet fourMomentum(track.px(), track.py(), track.pz(), track.energy(MassPionCharged));
        fjParticles.emplace_back(fourMomentum);
      }

      // Reject empty events
      if (fjParticles.size() < 1)
        continue;
      registryMC.fill(HIST("number_of_events_mc_rec"), 3.5);

      // Cluster particles using the anti-kt algorithm
      fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rJet);
      fastjet::AreaDefinition areaDef(fastjet::active_area, fastjet::GhostedAreaSpec(1.0));
      fastjet::ClusterSequenceArea cs(fjParticles, jetDef, areaDef);
      std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets());
      auto [rhoPerp, rhoMPerp] = backgroundSub.estimateRhoPerpCone(fjParticles, jets);

      // Jet selection
      bool isAtLeastOneJetSelected = false;
      std::vector<TVector3> selectedJet;
      std::vector<TVector3> ue1;
      std::vector<TVector3> ue2;

      // Loop over clustered jets
      for (const auto& jet : jets) {

        // jet must be fully contained in the acceptance
        if ((std::fabs(jet.eta()) + rJet) > (etaMax - deltaEtaEdge))
          continue;

        // jet pt must be larger than threshold
        auto jetForSub = jet;
        fastjet::PseudoJet jetMinusBkg = backgroundSub.doRhoAreaSub(jetForSub, rhoPerp, rhoMPerp);
        if (jetMinusBkg.pt() < minJetPt)
          continue;
        isAtLeastOneJetSelected = true;

        // Perpendicular cones
        TVector3 jetAxis(jet.px(), jet.py(), jet.pz());
        TVector3 ueAxis1(0, 0, 0), ueAxis2(0, 0, 0);
        getPerpendicularAxis(jetAxis, ueAxis1, +1);
        getPerpendicularAxis(jetAxis, ueAxis2, -1);

        // Store selected jet and UE cone axes
        selectedJet.emplace_back(jetAxis);
        ue1.emplace_back(ueAxis1);
        ue2.emplace_back(ueAxis2);
      }
      if (!isAtLeastOneJetSelected)
        continue;

      // Fill event counter for events with at least one selected jet
      registryMC.fill(HIST("number_of_events_mc_rec"), 4.5);

      // Loop over selected jets
      for (int i = 0; i < static_cast<int>(selectedJet.size()); i++) {

        // V0 particles
        if (particleOfInterest == Option::kV0Particles) {
          for (const auto& v0 : v0sPerColl) {
            const auto& pos = v0.posTrack_as<DaughterTracksMC>();
            const auto& neg = v0.negTrack_as<DaughterTracksMC>();
            TVector3 v0dir(v0.px(), v0.py(), v0.pz());

            // Get MC particles
            if (!pos.has_mcParticle() || !neg.has_mcParticle())
              continue;
            auto posParticle = pos.mcParticle_as<aod::McParticles>();
            auto negParticle = neg.mcParticle_as<aod::McParticles>();
            if (!posParticle.has_mothers() || !negParticle.has_mothers())
              continue;

            // Select particles originating from the same parent
            int pdgParent(0);
            bool isPhysPrim = false;
            for (const auto& particleMotherOfNeg : negParticle.mothers_as<aod::McParticles>()) {
              for (const auto& particleMotherOfPos : posParticle.mothers_as<aod::McParticles>()) {
                if (particleMotherOfNeg == particleMotherOfPos) {
                  pdgParent = particleMotherOfNeg.pdgCode();
                  isPhysPrim = particleMotherOfNeg.isPhysicalPrimary();
                }
              }
            }
            if (pdgParent == 0)
              continue;

            // Compute distance from jet and UE axes
            double deltaEtaJet = v0dir.Eta() - selectedJet[i].Eta();
            double deltaPhiJet = getDeltaPhi(v0dir.Phi(), selectedJet[i].Phi());
            double deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
            double deltaEtaUe1 = v0dir.Eta() - ue1[i].Eta();
            double deltaPhiUe1 = getDeltaPhi(v0dir.Phi(), ue1[i].Phi());
            double deltaRue1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
            double deltaEtaUe2 = v0dir.Eta() - ue2[i].Eta();
            double deltaPhiUe2 = getDeltaPhi(v0dir.Phi(), ue2[i].Phi());
            double deltaRue2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

            // K0s
            if (passedK0ShortSelection(v0, pos, neg) && pdgParent == kK0Short && isPhysPrim) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("K0s_reconstructed_jet"), multiplicity, v0.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("K0s_reconstructed_ue"), multiplicity, v0.pt());
              }
            }
            // Lambda
            if (passedLambdaSelection(v0, pos, neg) && pdgParent == kLambda0 && isPhysPrim) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("Lambda_reconstructed_jet"), multiplicity, v0.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("Lambda_reconstructed_ue"), multiplicity, v0.pt());
              }
            }
            // AntiLambda
            if (passedAntiLambdaSelection(v0, pos, neg) && pdgParent == kLambda0Bar && isPhysPrim) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("AntiLambda_reconstructed_jet"), multiplicity, v0.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("AntiLambda_reconstructed_ue"), multiplicity, v0.pt());
              }
            }

            // Fill inclusive spectra
            // K0s
            if (passedK0ShortSelection(v0, pos, neg) && pdgParent == kK0Short) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("K0s_reconstructed_jet_incl"), multiplicity, v0.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("K0s_reconstructed_ue_incl"), multiplicity, v0.pt());
              }
            }
            // Lambda
            if (passedLambdaSelection(v0, pos, neg) && pdgParent == kLambda0) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("Lambda_reconstructed_jet_incl"), multiplicity, v0.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("Lambda_reconstructed_ue_incl"), multiplicity, v0.pt());
              }
            }
            // AntiLambda
            if (passedAntiLambdaSelection(v0, pos, neg) && pdgParent == kLambda0Bar) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("AntiLambda_reconstructed_jet_incl"), multiplicity, v0.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("AntiLambda_reconstructed_ue_incl"), multiplicity, v0.pt());
              }
            }
          }
        }

        // Cascades
        if (particleOfInterest == Option::kCascades) {
          for (const auto& casc : cascPerColl) {
            auto bach = casc.bachelor_as<DaughterTracksMC>();
            auto pos = casc.posTrack_as<DaughterTracksMC>();
            auto neg = casc.negTrack_as<DaughterTracksMC>();

            // Get MC particles
            if (!bach.has_mcParticle() || !pos.has_mcParticle() || !neg.has_mcParticle())
              continue;
            auto posParticle = pos.mcParticle_as<aod::McParticles>();
            auto negParticle = neg.mcParticle_as<aod::McParticles>();
            auto bachParticle = bach.mcParticle_as<aod::McParticles>();
            if (!posParticle.has_mothers() || !negParticle.has_mothers() || !bachParticle.has_mothers())
              continue;

            // Select particles originating from the same parent
            int pdgParent(0);
            bool isPhysPrim = false;
            for (const auto& particleMotherOfNeg : negParticle.mothers_as<aod::McParticles>()) {
              for (const auto& particleMotherOfPos : posParticle.mothers_as<aod::McParticles>()) {
                for (const auto& particleMotherOfBach : bachParticle.mothers_as<aod::McParticles>()) {
                  if (particleMotherOfNeg != particleMotherOfPos)
                    continue;
                  if (std::abs(particleMotherOfNeg.pdgCode()) != kLambda0)
                    continue;
                  isPhysPrim = particleMotherOfBach.isPhysicalPrimary();
                  pdgParent = particleMotherOfBach.pdgCode();
                }
              }
            }
            if (pdgParent == 0)
              continue;
            if (!isPhysPrim)
              continue;

            // Compute distances from jet and UE axes
            TVector3 cascadeDir(casc.px(), casc.py(), casc.pz());
            double deltaEtaJet = cascadeDir.Eta() - selectedJet[i].Eta();
            double deltaPhiJet = getDeltaPhi(cascadeDir.Phi(), selectedJet[i].Phi());
            double deltaRjet = std::sqrt(deltaEtaJet * deltaEtaJet + deltaPhiJet * deltaPhiJet);
            double deltaEtaUe1 = cascadeDir.Eta() - ue1[i].Eta();
            double deltaPhiUe1 = getDeltaPhi(cascadeDir.Phi(), ue1[i].Phi());
            double deltaRue1 = std::sqrt(deltaEtaUe1 * deltaEtaUe1 + deltaPhiUe1 * deltaPhiUe1);
            double deltaEtaUe2 = cascadeDir.Eta() - ue2[i].Eta();
            double deltaPhiUe2 = getDeltaPhi(cascadeDir.Phi(), ue2[i].Phi());
            double deltaRue2 = std::sqrt(deltaEtaUe2 * deltaEtaUe2 + deltaPhiUe2 * deltaPhiUe2);

            // Xi+
            if (passedXiSelection(casc, pos, neg, bach, collision) && bach.sign() > 0 && pdgParent == kXiPlusBar) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("XiPos_reconstructed_jet"), multiplicity, casc.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("XiPos_reconstructed_ue"), multiplicity, casc.pt());
              }
            }
            // Xi-
            if (passedXiSelection(casc, pos, neg, bach, collision) && bach.sign() < 0 && pdgParent == kXiMinus) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("XiNeg_reconstructed_jet"), multiplicity, casc.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("XiNeg_reconstructed_ue"), multiplicity, casc.pt());
              }
            }
            // Omega+
            if (passedOmegaSelection(casc, pos, neg, bach, collision) && bach.sign() > 0 && pdgParent == kOmegaPlusBar) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("OmegaPos_reconstructed_jet"), multiplicity, casc.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("OmegaPos_reconstructed_ue"), multiplicity, casc.pt());
              }
            }
            // Omega-
            if (passedOmegaSelection(casc, pos, neg, bach, collision) && bach.sign() < 0 && pdgParent == kOmegaMinus) {
              if (deltaRjet < rJet) {
                registryMC.fill(HIST("OmegaNeg_reconstructed_jet"), multiplicity, casc.pt());
              }
              if (deltaRue1 < rJet || deltaRue2 < rJet) {
                registryMC.fill(HIST("OmegaNeg_reconstructed_ue"), multiplicity, casc.pt());
              }
            }
          }
        }
      }
    }
  }
  PROCESS_SWITCH(StrangenessInJets, processMCreconstructed, "process reconstructed events", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<StrangenessInJets>(cfgc)};
}
