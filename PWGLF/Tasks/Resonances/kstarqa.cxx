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
/// \brief this is a code for the kstarqa resonance
/// \author prottay das, sawan
/// \since 13/03/2024

#include <TDatabasePDG.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPDGCode.h>
#include "TRandom3.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include "CCDB/BasicCCDBManager.h"
#include "CCDB/CcdbApi.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/runDataProcessing.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "ReconstructionDataFormats/Track.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::soa;
using std::array;

struct kstarqa {

  SliceCache cache;

  // Histograms are defined with HistogramRegistry
  HistogramRegistry rEventSelection{"eventSelection", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  // Confugrable for QA histograms
  Configurable<bool> CalcLikeSign{"CalcLikeSign", true, "Calculate Like Sign"};
  Configurable<bool> CalcRotational{"CalcRotational", true, "Calculate Rotational"};
  Configurable<bool> QA{"QA", false, "QA"};
  Configurable<bool> QAbefore{"QAbefore", true, "QAbefore"};
  Configurable<bool> QAafter{"QAafter", true, "QAafter"};
  Configurable<bool> QAevents{"QAevents", true, "Multiplicity dist, DCAxy, DCAz"};
  Configurable<bool> onlyTOF{"onlyTOF", false, "only TOF tracks"};
  Configurable<bool> onlyTOFHIT{"onlyTOFHIT", false, "accept only TOF hit tracks at high pt"};
  Configurable<bool> onlyTPC{"onlyTPC", true, "only TPC tracks"};
  Configurable<bool> cfgFT0M{"cfgFT0M", true, "1: pp, 0: PbPb"};

  // Configurables for track selections
  Configurable<int> rotational_cut{"rotational_cut", 10, "Cut value (Rotation angle pi - pi/cut and pi + pi/cut)"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2f, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPCPi{"nsigmacutTPCPi", 3.0, "Value of the TPC Nsigma cut for pions"};
  Configurable<float> nsigmaCutTPCKa{"nsigmacutTPCKa", 3.0, "Value of the TPC Nsigma cut for kaons"};
  Configurable<float> nsigmaCutTOFPi{"nsigmacutTOFPi", 3.0, "Value of the TOF Nsigma cut for pions"};
  Configurable<float> nsigmaCutTOFKa{"nsigmacutTOFKa", 3.0, "Value of the TOF Nsigma cut for kaons"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the Combined Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 5, "Number of mixed events per event"};
  Configurable<bool> ismanualDCAcut{"ismanualDCAcut", true, "ismanualDCAcut"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 120, "Number of TPC cluster"};
  Configurable<float> cfgRCRFC{"cfgRCRFC", 0.8f, "Crossed Rows to Findable Clusters"};
  Configurable<float> cfgITSChi2NCl{"cfgITSChi2NCl", 36.0, "ITS Chi2/NCl"};
  Configurable<float> cfgTPCChi2NCl{"cfgTPCChi2NCl", 4.0, "TPC Chi2/NCl"};
  Configurable<bool> cfgPVContributor{"cfgPVContributor", false, "PV contributor track selection"};           // PV Contriuibutor
  Configurable<bool> cfgPrimaryTrack{"cfgPrimaryTrack", false, "Primary track selection"};                    // kGoldenChi2 | kDCAxy | kDCAz
  Configurable<bool> cfgGlobalWoDCATrack{"cfgGlobalWoDCATrack", false, "Global track selection without DCA"}; // kQualityTracks (kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF | kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits) | kInAcceptanceTracks (kPtRange | kEtaRange)
  Configurable<bool> cfgGlobalTrack{"cfgGlobalTrack", false, "Global track selection"};                       // kGoldenChi2 | kDCAxy | kDCAz

  // Event selection configurables
  Configurable<bool> timFrameEvsel{"timFrameEvsel", false, "TPC Time frame boundary cut"};
  Configurable<bool> TVXEvsel{"TVXEvsel", false, "Triggger selection"};
  Configurable<float> cutzvertex{"cutzvertex", 10.0f, "Accepted z-vertex range (cm)"};
  Configurable<bool> MID{"MID", false, "Misidentification of tracks"};

  // Configurable for histograms
  Configurable<int> nBins{"nBins", 100, "N bins in all histos"};
  Configurable<int> nBinsinvMass{"nBinsinvMass", 180, "N bins in invMass histos"};
  Configurable<float> invMassbinlow{"invMassbinlow", 0.6, "invMass bin low"};
  Configurable<float> invMassbinhigh{"invMassbinhigh", 1.5, "invMass bin high"};
  Configurable<int> nBinspT{"nBinspT", 200, "N bins in pT histos"};
  Configurable<float> pTbinlow{"pTbinlow", 0.0, "pT bin low"};
  Configurable<float> pTbinhigh{"pTbinhigh", 20.0, "pT bin high"};
  ConfigurableAxis binsMultPlot{"binsCent", {201, -0.5f, 200.5f}, "Binning of the centrality axis for plots"};
  Configurable<bool> avoidsplitrackMC{"avoidsplitrackMC", true, "avoid split track in MC"};
  Configurable<bool> AllGenCollisions{"AllGenCollisions", true, "To fill all generated collisions for the signal loss calculations"};
  ConfigurableAxis axisdEdx{"axisdEdx", {20000, 0.0f, 200.0f}, "dE/dx (a.u.)"};
  ConfigurableAxis axisPtfordEbydx{"axisPtfordEbydx", {2000, 0, 20}, "pT (GeV/c)"};
  ConfigurableAxis axisMultdist{"axisMultdist", {3500, 0, 70000}, "Multiplicity distribution"};

  // Event plane configurables
  Configurable<bool> boost_daugter1{"boost_daugter1", false, "Boost daughter Kaon in the COM frame"};
  Configurable<bool> boost_daugter2{"boost_daugter2", true, "Boost daughter Pion in the COM frame"};
  Configurable<bool> activateTHnSparseCosThStarHelicity{"activateTHnSparseCosThStarHelicity", true, "Activate the THnSparse with cosThStar w.r.t. helicity axis"};
  Configurable<bool> activateTHnSparseCosThStarProduction{"activateTHnSparseCosThStarProduction", false, "Activate the THnSparse with cosThStar w.r.t. production axis"};
  Configurable<bool> activateTHnSparseCosThStarBeam{"activateTHnSparseCosThStarBeam", false, "Activate the THnSparse with cosThStar w.r.t. beam axis (Gottified jackson frame)"};
  Configurable<bool> activateTHnSparseCosThStarRandom{"activateTHnSparseCosThStarRandom", false, "Activate the THnSparse with cosThStar w.r.t. random axis"};
  Configurable<int> c_nof_rotations{"c_nof_rotations", 3, "Number of random rotations in the rotational background"};
  ConfigurableAxis configThnAxisPOL{"configThnAxisPOL", {20, -1.0, 1.0}, "Costheta axis"};

  void init(InitContext const&)
  {
    // Axes
    AxisSpec vertexZAxis = {nBins, -15., 15., "vrtx_{Z} [cm] for plots"};
    AxisSpec ptAxis = {nBinspT, pTbinlow, pTbinhigh, "#it{p}_{T} (GeV/#it{c})"};
    AxisSpec invmassAxis = {nBinsinvMass, invMassbinlow, invMassbinhigh, "Invariant mass (GeV/#it{c}^{2})"};
    AxisSpec thnAxisPOL{configThnAxisPOL, "cos(#theta)"};

    // Histograms
    // Event selection
    rEventSelection.add("hVertexZRec", "hVertexZRec", {HistType::kTH1F, {vertexZAxis}});
    rEventSelection.add("hmult", "Multiplicity percentile", kTH1F, {{binsMultPlot}});

    // for primary tracks
    if (QAbefore && QAafter) {
      histos.add("hNsigmaTPC_before", "NsigmaKaon TPC distribution before", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      histos.add("hNsigmaTOF_before", "NsigmaKaon TOF distribution before", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      histos.add("dE_by_dx_TPC", "dE/dx signal in the TPC as a function of pT", kTH2F, {axisPtfordEbydx, axisdEdx});
      histos.add("hphi", "Phi distribution", kTH1F, {{65, 0, 6.5}});

      histos.add("hEta_after", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
      histos.add("hCRFC_after", "CRFC after distribution", kTH1F, {{100, 0.0f, 10.0f}});
      histos.add("hCRFC_before", "CRFC before distribution", kTH1F, {{100, 0.0f, 10.0f}});
      // histos.add("hNsigmaPionTPC_after", "NsigmaPion TPC distribution", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      // histos.add("hNsigmaPionTOF_after", "NsigmaPion TOF distribution", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      // histos.add("hNsigmaKaonTPC_after", "NsigmaKaon TPC distribution", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
      // histos.add("hNsigmaKaonTOF_after", "NsigmaKaon TOF distribution", kTH2F, {{100, 0.0f, 10.0f}, {200, -10.0f, 10.0f}});
    }

    // KStar histograms
    histos.add("h3KstarInvMassUnlikeSign", "kstar Unlike Sign", kTHnSparseF, {binsMultPlot, ptAxis, invmassAxis, thnAxisPOL}, true);
    histos.add("h3KstarInvMassMixed", "kstar Mixed", kTHnSparseF, {binsMultPlot, ptAxis, invmassAxis, thnAxisPOL}, true);
    if (CalcLikeSign)
      histos.add("h3KstarInvMasslikeSign", "kstar like Sign", kTHnSparseF, {binsMultPlot, ptAxis, invmassAxis, thnAxisPOL}, true);
    if (CalcRotational)
      histos.add("h3KstarInvMassRotated", "kstar rotated", kTHnSparseF, {binsMultPlot, ptAxis, invmassAxis, thnAxisPOL}, true);

    // MC generated histograms
    histos.add("k892Gen", "pT distribution of True MC K(892)0", kTH1D, {ptAxis});
    // histos.add("k892GenAnti", "pT distribution of True MC Anti-K(892)0", kTH1D, {ptAxis});
    // Reconstructed MC histogram
    histos.add("h3KstarRec", "pT of reconstructed kstar", kTH1F, {ptAxis});
    histos.add("h1KstarRecMass", "Invariant mass of kstar meson", kTH1D, {invmassAxis});
    histos.add("h1KstarRecpt", "pT of kstar meson", kTH1F, {ptAxis});
    histos.add("h1genmass", "Invariant mass of generated kstar meson", kTH1D, {invmassAxis});
    histos.add("h1recpt", "pT of generated kstar meson", kTH1F, {ptAxis});
    histos.add("events_check_data", "No. of events in the data", kTH1I, {{20, 0, 20}});
    histos.add("events_check", "No. of events in the generated MC", kTH1I, {{20, 0, 20}});
    histos.add("events_checkrec", "No. of events in the reconstructed MC", kTH1I, {{20, 0, 20}});
    histos.add("h1KSRecsplit", "KS meson Rec split", kTH1F, {{100, 0.0f, 10.0f}});

    // Multplicity distribution
    if (QAevents) {
      histos.add("multdist_FT0M", "FT0M Multiplicity distribution", kTH1F, {axisMultdist});
      histos.add("multdist_FT0A", "FT0A Multiplicity distribution", kTH1F, {axisMultdist});
      histos.add("multdist_FT0C", "FT0C Multiplicity distribution", kTH1F, {axisMultdist});
      histos.add("hNcontributor", "Number of primary vertex contributor", kTH1F, {{2000, 0.0f, 10000.0f}});
      histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
      histos.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    }
  }

  double massPi = TDatabasePDG::Instance()->GetParticle(kPiPlus)->Mass(); // FIXME: Get from the common header
  double massKa = o2::constants::physics::MassKPlus;
  ROOT::Math::PtEtaPhiMVector CKSVector;
  ROOT::Math::PtEtaPhiMVector CKSVectorRot1;
  ROOT::Math::PtEtaPhiMVector CKSVectorRot2;

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (ismanualDCAcut && !(candidate.isGlobalTrackWoDCA() && candidate.isPVContributor() && std::abs(candidate.dcaXY()) < cfgCutDCAxy && std::abs(candidate.dcaZ()) < cfgCutDCAz && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster)) {
      return false;
    } else if (!ismanualDCAcut) {
      if (std::abs(candidate.pt()) < cfgCutPT)
        return false;
      if (std::abs(candidate.dcaXY()) > cfgCutDCAxy)
        return false;
      if (std::abs(candidate.dcaZ()) > cfgCutDCAz)
        return false;
      if (candidate.tpcCrossedRowsOverFindableCls() < cfgRCRFC)
        return false;
      if (candidate.itsNCls() < cfgITScluster)
        return false;
      if (candidate.tpcNClsFound() < cfgTPCcluster)
        return false;
      if (candidate.itsChi2NCl() >= cfgITSChi2NCl)
        return false;
      if (candidate.tpcChi2NCl() >= cfgTPCChi2NCl)
        return false;
      if (cfgPVContributor && !candidate.isPVContributor())
        return false;
      // if (cfgPrimaryTrack && !candidate.isPrimaryTrack())
      //   return false;
      if (cfgGlobalWoDCATrack && !candidate.isGlobalTrackWoDCA())
        return false;
      if (cfgGlobalTrack && !candidate.isGlobalTrack())
        return false;
    }

    return true;
  }

  template <typename T>
  bool selectionPID(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOFPi) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < nsigmaCutTOFPi) {
          return true;
        }
        if (!candidate.hasTOF() &&
            std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPCPi) {
          return true;
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPCPi) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (nsigmaCutCombined * nsigmaCutCombined)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < nsigmaCutTPCPi) {
          return true;
        }
      }
    } else if (PID == 1) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOFKa) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < nsigmaCutTOFKa) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPCKa) {
          return true;
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPCKa) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (nsigmaCutCombined * nsigmaCutCombined)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < nsigmaCutTPCKa) {
          return true;
        }
      }
    }
    return false;
  }

  template <typename T>
  bool MIDselectionPID(const T& candidate, int PID)
  {
    if (PID == 0) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < 3.0) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPi()) < 3.0) {
          return true;
        }
        if (!candidate.hasTOF() &&
            std::abs(candidate.tpcNSigmaPi()) < 3.0) {
          return true;
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaPi()) < 3.0) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaPi() * candidate.tofNSigmaPi() + candidate.tpcNSigmaPi() * candidate.tpcNSigmaPi()) < (3.0 * 3.0)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPi()) < 3.0) {
          return true;
        }
      }
    } else if (PID == 1) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < 3.0) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaKa()) < 3.0) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < 3.0) {
          return true;
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaKa()) < 3.0) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaKa() * candidate.tofNSigmaKa() + candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa()) < (3.0 * 3.0)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaKa()) < 3.0) {
          return true;
        }
      }
    } else if (PID == 2) {
      if (onlyTOF) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < 3.0) {
          return true;
        }
      } else if (onlyTOFHIT) {
        if (candidate.hasTOF() && std::abs(candidate.tofNSigmaPr()) < 3.0) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < 3.0) {
          return true;
        }
      } else if (onlyTPC) {
        if (std::abs(candidate.tpcNSigmaPr()) < 3.0) {
          return true;
        }
      } else {
        if (candidate.hasTOF() && (candidate.tofNSigmaPr() * candidate.tofNSigmaPr() + candidate.tpcNSigmaPr() * candidate.tpcNSigmaPr()) < (3.0 * 3.0)) {
          return true;
        }
        if (!candidate.hasTOF() && std::abs(candidate.tpcNSigmaPr()) < 3.0) {
          return true;
        }
      }
    }
    return false;
  }

  array<float, 3> pvec0;
  array<float, 3> pvec1;

  // Defining filters for events (event selection)
  // Processed events will be already fulfilling the event selection
  // requirements
  // Filter eventFilter = (o2::aod::evsel::sel8 == true);
  Filter posZFilter = (nabs(o2::aod::collision::posZ) < cutzvertex);

  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::MultZeqs, aod::CentFT0Ms, aod::CentFT0As, aod::CentFT0Cs>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPr, aod::pidTOFFullPr, aod::TrackSelectionExtension>>;
  using V0TrackCandidate = aod::V0Datas;
  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;

  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::pidTPCFullPi, aod::pidTOFFullPi, aod::McTrackLabels>>;

  template <typename T1, typename T2, typename T3, typename T4>
  void fillInvMass(const T1& track1, const T2& track2, const T3& lv2, const T4& lv3, float multiplicity, bool isMix)
  {
    ROOT::Math::PxPyPzMVector daughter1, daughter2, daughter_selected;
    daughter1 = ROOT::Math::PxPyPzMVector(track1.px(), track1.py(), track1.pz(), massKa); // Kaon
    daughter2 = ROOT::Math::PxPyPzMVector(track2.px(), track2.py(), track2.pz(), massPi); // Pion
    daughter_selected = (boost_daugter1) ? daughter1 : daughter2;
    auto selected_dau_mass = (boost_daugter1) ? massKa : massPi;

    TLorentzVector lv4, lv5;
    // polarization calculations

    auto phiRandom = gRandom->Uniform(0.f, constants::math::TwoPI);
    auto thetaRandom = gRandom->Uniform(0.f, constants::math::PI);
    ROOT::Math::PxPyPzMVector fourVecDau1 = ROOT::Math::PxPyPzMVector(daughter_selected.Px(), daughter_selected.Py(), daughter_selected.Pz(), selected_dau_mass); // Kaon or Pion

    ROOT::Math::PxPyPzMVector fourVecMother = ROOT::Math::PxPyPzMVector(lv3.Px(), lv3.Py(), lv3.Pz(), lv3.M()); // mass of KshortKshort pair
    ROOT::Math::Boost boost{fourVecMother.BoostToCM()};                                                         // boost mother to center of mass frame
    ROOT::Math::PxPyPzMVector fourVecDauCM = boost(fourVecDau1);                                                // boost the frame of daughter same as mother
    ROOT::Math::XYZVector threeVecDauCM = fourVecDauCM.Vect();                                                  // get the 3 vector of daughter in the frame of mother

    TRandom* rn = new TRandom();

    if (TMath::Abs(lv3.Rapidity() < 0.5)) {

      if (activateTHnSparseCosThStarHelicity) {
        ROOT::Math::XYZVector helicityVec = fourVecMother.Vect(); // 3 vector of mother in COM frame
        auto cosThetaStarHelicity = helicityVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(helicityVec.Mag2()));

        if (track1.sign() * track2.sign() < 0) {
          if (!isMix) {
            histos.fill(HIST("h3KstarInvMassUnlikeSign"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarHelicity);

            for (int i = 0; i < c_nof_rotations; i++) {
              float theta2 = rn->Uniform(TMath::Pi() - TMath::Pi() / rotational_cut, TMath::Pi() + TMath::Pi() / rotational_cut);
              lv4.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi() + theta2, massKa); // for rotated background
              lv5 = lv2 + lv4;
              if (CalcRotational)
                histos.fill(HIST("h3KstarInvMassRotated"), multiplicity, lv5.Pt(), lv5.M(), cosThetaStarHelicity);
            }
          } else {
            histos.fill(HIST("h3KstarInvMassMixed"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarHelicity);
          }
        } else {
          if (!isMix) {
            if (CalcLikeSign)
              histos.fill(HIST("h3KstarInvMasslikeSign"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarHelicity);
          }
        }

      } else if (activateTHnSparseCosThStarProduction) {
        ROOT::Math::XYZVector normalVec = ROOT::Math::XYZVector(lv3.Py(), -lv3.Px(), 0.f);
        auto cosThetaStarProduction = normalVec.Dot(threeVecDauCM) / (std::sqrt(threeVecDauCM.Mag2()) * std::sqrt(normalVec.Mag2()));

        if (track1.sign() * track2.sign() < 0) {
          if (!isMix) {
            histos.fill(HIST("h3KstarInvMassUnlikeSign"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarProduction);
            for (int i = 0; i < c_nof_rotations; i++) {
              float theta2 = rn->Uniform(0, TMath::Pi());
              lv4.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi() + theta2, massKa); // for rotated background
              lv5 = lv2 + lv4;
              if (CalcRotational)
                histos.fill(HIST("h3KstarInvMassRotated"), multiplicity, lv5.Pt(), lv5.M(), cosThetaStarProduction);
            }
          } else {
            histos.fill(HIST("h3KstarInvMassMixed"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarProduction);
          }
        } else {
          if (!isMix) {
            if (CalcLikeSign)
              histos.fill(HIST("h3KstarInvMasslikeSign"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarProduction);
          }
        }
      } else if (activateTHnSparseCosThStarBeam) {
        ROOT::Math::XYZVector beamVec = ROOT::Math::XYZVector(0.f, 0.f, 1.f);
        auto cosThetaStarBeam = beamVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());

        if (track1.sign() * track2.sign() < 0) {
          if (!isMix) {
            histos.fill(HIST("h3KstarInvMassUnlikeSign"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarBeam);
            for (int i = 0; i < c_nof_rotations; i++) {
              float theta2 = rn->Uniform(0, TMath::Pi());
              lv4.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi() + theta2, massKa); // for rotated background
              lv5 = lv2 + lv4;
              if (CalcRotational)
                histos.fill(HIST("h3KstarInvMassRotated"), multiplicity, lv5.Pt(), lv5.M(), cosThetaStarBeam);
            }
          } else {
            histos.fill(HIST("h3KstarInvMassMixed"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarBeam);
          }
        } else {
          if (CalcLikeSign)
            histos.fill(HIST("h3KstarInvMasslikeSign"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarBeam);
        }
      } else if (activateTHnSparseCosThStarRandom) {
        ROOT::Math::XYZVector randomVec = ROOT::Math::XYZVector(std::sin(thetaRandom) * std::cos(phiRandom), std::sin(thetaRandom) * std::sin(phiRandom), std::cos(thetaRandom));
        auto cosThetaStarRandom = randomVec.Dot(threeVecDauCM) / std::sqrt(threeVecDauCM.Mag2());

        if (track1.sign() * track2.sign() < 0) {
          if (!isMix) {
            histos.fill(HIST("h3KstarInvMassUnlikeSign"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarRandom);
            for (int i = 0; i < c_nof_rotations; i++) {
              float theta2 = rn->Uniform(0, TMath::Pi());
              lv4.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi() + theta2, massKa); // for rotated background
              lv5 = lv2 + lv4;
              if (CalcRotational)
                histos.fill(HIST("h3KstarInvMassRotated"), multiplicity, lv5.Pt(), lv5.M(), cosThetaStarRandom);
            }
          } else {
            histos.fill(HIST("h3KstarInvMassMixed"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarRandom);
          }
        } else {
          if (!isMix) {
            if (CalcLikeSign)
              histos.fill(HIST("h3KstarInvMasslikeSign"), multiplicity, lv3.Pt(), lv3.M(), cosThetaStarRandom);
          }
        }
      }
    }
  }

  void processSE(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)

  {
    histos.fill(HIST("events_check_data"), 0.5);

    if (TVXEvsel && (!collision.selection_bit(aod::evsel::kIsTriggerTVX))) {
      return;
    }
    histos.fill(HIST("events_check_data"), 1.5);

    if (timFrameEvsel && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }
    histos.fill(HIST("events_check_data"), 2.5);

    if (!collision.sel8()) {
      return;
    }
    histos.fill(HIST("events_check_data"), 3.5);

    float multiplicity = 0.0f;

    multiplicity = (cfgFT0M) ? collision.centFT0M() : collision.centFT0C();

    // Fill the event counter
    if (QAevents) {
      rEventSelection.fill(HIST("hVertexZRec"), collision.posZ());
      rEventSelection.fill(HIST("hmult"), multiplicity);
      histos.fill(HIST("multdist_FT0M"), collision.multFT0M());
      histos.fill(HIST("multdist_FT0A"), collision.multFT0A());
      histos.fill(HIST("multdist_FT0C"), collision.multFT0C());
      histos.fill(HIST("hNcontributor"), collision.numContrib());
    }

    for (auto& [track1, track2] : combinations(CombinationsFullIndexPolicy(tracks, tracks))) {
      if (QAbefore) {
        histos.fill(HIST("hNsigmaTPC_before"), track1.pt(), track1.tpcNSigmaKa());
        histos.fill(HIST("hNsigmaTOF_before"), track1.pt(), track1.tofNSigmaKa());
        histos.fill(HIST("hCRFC_before"), track1.tpcCrossedRowsOverFindableCls());
        histos.fill(HIST("dE_by_dx_TPC"), track1.p(), track1.tpcSignal());
        histos.fill(HIST("hphi"), track1.phi());
      }
      histos.fill(HIST("events_check_data"), 4.5);

      if (!selectionTrack(track1)) {
        continue;
      }
      if (!selectionTrack(track2)) {
        continue;
      }
      histos.fill(HIST("events_check_data"), 5.5);
      if (QAevents) {
        histos.fill(HIST("hDcaxy"), track1.dcaXY());
        histos.fill(HIST("hDcaz"), track1.dcaZ());
      }

      // since we are using combinations full index policy, so repeated pairs are allowed, so we can check one with Kaon and other with pion
      if (!selectionPID(track1, 1)) // Track 1 is checked with Kaon
        continue;
      if (!selectionPID(track2, 0)) // Track 2 is checked with Pion
        continue;

      histos.fill(HIST("events_check_data"), 6.5);

      if (MID) {
        if (MIDselectionPID(track1, 0)) // Kaon misidentified as pion
          continue;
        if (MIDselectionPID(track1, 2)) // Kaon misidentified as proton
          continue;
        if (MIDselectionPID(track2, 1)) // Pion misidentified as kaon
          continue;
      }

      histos.fill(HIST("events_check_data"), 7.5);

      if (QAafter) {
        histos.fill(HIST("hEta_after"), track1.eta());
        histos.fill(HIST("hCRFC_after"), track1.tpcCrossedRowsOverFindableCls());
        // histos.fill(HIST("hNsigmaKaonTPC_after"), track1.pt(), track1.tpcNSigmaKa());
        // histos.fill(HIST("hNsigmaKaonTOF_after"), track1.pt(), track1.tofNSigmaKa());
        // histos.fill(HIST("hNsigmaPionTPC_after"), track2.pt(), track2.tpcNSigmaPi());
        // histos.fill(HIST("hNsigmaPionTOF_after"), track2.pt(), track2.tofNSigmaPi());
      }

      if (track1.globalIndex() == track2.globalIndex())
        continue;

      histos.fill(HIST("events_check_data"), 8.5);

      TLorentzVector lv1, lv2, lv3;
      lv1.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi(), massKa);
      lv2.SetPtEtaPhiM(track2.pt(), track2.eta(), track2.phi(), massPi);
      lv3 = lv1 + lv2;
      bool isMix = false;
      fillInvMass(track1, track2, lv2, lv3, multiplicity, isMix);
    }
  }

  PROCESS_SWITCH(kstarqa, processSE, "Process Same event", true);

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for ME mixing"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {10, 0, 100}, "multiplicity percentile for ME mixing"};
  // ConfigurableAxis axisMultiplicity{"axisMultiplicity", {2000, 0, 10000}, "TPC multiplicity  for bin for ME mixing"};

  using BinningTypeTPCMultiplicity = ColumnBinningPolicy<aod::collision::PosZ, aod::mult::MultTPC>;
  using BinningTypeCentralityM = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0M>;

  BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicityClass}, true};
  BinningTypeCentralityM binningOnCentrality{{axisVertex, axisMultiplicityClass}, true};

  SameKindPair<EventCandidates, TrackCandidates, BinningTypeVertexContributor> pair1{binningOnPositions, cfgNoMixedEvents, -1, &cache};
  SameKindPair<EventCandidates, TrackCandidates, BinningTypeCentralityM> pair2{binningOnCentrality, cfgNoMixedEvents, -1, &cache};

  void processME(EventCandidates const&, TrackCandidates const&)
  {
    if (cfgFT0M) {
      for (auto& [c1, tracks1, c2, tracks2] : pair1) {

        if (!c1.sel8()) {
          continue;
        }
        if (!c2.sel8()) {
          continue;
        }

        if (timFrameEvsel && (!c1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c2.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !c2.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
          continue;
        }

        if (TVXEvsel && (!c1.selection_bit(aod::evsel::kIsTriggerTVX) || !c2.selection_bit(aod::evsel::kIsTriggerTVX))) {
          return;
        }

        auto multiplicity = c1.centFT0M();

        for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

          if (!selectionTrack(t1)) // Kaon
            continue;
          if (!selectionTrack(t2)) // Pion
            continue;
          if (!selectionPID(t1, 1)) // Kaon
            continue;
          if (!selectionPID(t2, 0)) // Pion
            continue;
          if (MID) {
            if (MIDselectionPID(t1, 0)) // misidentified as pion
              continue;
            if (MIDselectionPID(t1, 2)) // misidentified as proton
              continue;
            if (MIDselectionPID(t2, 1)) // misidentified as kaon
              continue;
          }

          TLorentzVector KAON;
          KAON.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), massKa);
          TLorentzVector PION;
          PION.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), massPi);

          TLorentzVector Kstar = KAON + PION;
          bool isMix = true;

          if (!QA) {
            if (TMath::Abs(Kstar.Rapidity()) < 0.5) {
              fillInvMass(t1, t2, PION, Kstar, multiplicity, isMix);
            }
          }
        }
      }
    } else {
      for (auto& [c1, tracks1, c2, tracks2] : pair2) {

        if (!c1.sel8()) {
          continue;
        }
        if (!c2.sel8()) {
          continue;
        }

        if (timFrameEvsel && (!c1.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c2.selection_bit(aod::evsel::kNoTimeFrameBorder) || !c1.selection_bit(aod::evsel::kNoITSROFrameBorder) || !c2.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
          continue;
        }

        if (TVXEvsel && (!c1.selection_bit(aod::evsel::kIsTriggerTVX) || !c2.selection_bit(aod::evsel::kIsTriggerTVX))) {
          return;
        }

        auto multiplicity = c1.centFT0M();

        for (auto& [t1, t2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(tracks1, tracks2))) {

          if (!selectionTrack(t1)) // Kaon
            continue;
          if (!selectionTrack(t2)) // Pion
            continue;
          if (!selectionPID(t1, 1)) // Kaon
            continue;
          if (!selectionPID(t2, 0)) // Pion
            continue;
          if (MID) {
            if (MIDselectionPID(t1, 0)) // misidentified as pion
              continue;
            if (MIDselectionPID(t1, 2)) // misidentified as proton
              continue;
            if (MIDselectionPID(t2, 1)) // misidentified as kaon
              continue;
          }

          TLorentzVector KAON;
          KAON.SetPtEtaPhiM(t1.pt(), t1.eta(), t1.phi(), massKa);
          TLorentzVector PION;
          PION.SetPtEtaPhiM(t2.pt(), t2.eta(), t2.phi(), massPi);

          TLorentzVector Kstar = KAON + PION;
          bool isMix = true;

          if (!QA) {
            if (TMath::Abs(Kstar.Rapidity()) < 0.5) {
              fillInvMass(t1, t2, PION, Kstar, multiplicity, isMix);
            }
          }
        }
      }
    }
  }

  PROCESS_SWITCH(kstarqa, processME, "Process Mixed event", true);

  void processGen(aod::McCollision const& mcCollision, aod::McParticles& mcParticles, const soa::SmallGroups<EventCandidatesMC>& collisions)
  {
    histos.fill(HIST("events_check"), 0.5);
    if (std::abs(mcCollision.posZ()) < cutzvertex) {
      histos.fill(HIST("events_check"), 1.5);
    }
    int Nchinel = 0;
    for (auto& mcParticle : mcParticles) {
      auto pdgcode = std::abs(mcParticle.pdgCode());
      if (mcParticle.isPhysicalPrimary() && (pdgcode == 211 || pdgcode == 321 || pdgcode == 2212 || pdgcode == 11 || pdgcode == 13)) {
        if (std::abs(mcParticle.eta()) < 1.0) {
          Nchinel = Nchinel + 1;
        }
      }
    }
    if (Nchinel > 0 && std::abs(mcCollision.posZ()) < cutzvertex)
      histos.fill(HIST("events_check"), 2.5);
    std::vector<int64_t> SelectedEvents(collisions.size());
    int nevts = 0;
    for (const auto& collision : collisions) {
      // if (!collision.sel8() || std::abs(collision.mcCollision().posZ()) > cutzvertex) {
      if (std::abs(collision.mcCollision().posZ()) > cutzvertex) {
        continue;
      }

      if (timFrameEvsel && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
        continue;
      }
      if (TVXEvsel && (!collision.selection_bit(aod::evsel::kIsTriggerTVX))) {
        continue;
      }

      SelectedEvents[nevts++] = collision.mcCollision_as<aod::McCollisions>().globalIndex();
    }
    SelectedEvents.resize(nevts);
    histos.fill(HIST("events_check"), 3.5);
    const auto evtReconstructedAndSelected = std::find(SelectedEvents.begin(), SelectedEvents.end(), mcCollision.globalIndex()) != SelectedEvents.end();

    if (!AllGenCollisions && !evtReconstructedAndSelected) { // Check that the event is reconstructed and that the reconstructed events pass the selection
      return;
    }
    histos.fill(HIST("events_check"), 4.5);
    for (auto& mcParticle : mcParticles) {
      if (std::abs(mcParticle.y()) >= 0.5) {
        continue;
      }
      histos.fill(HIST("events_check"), 5.5);
      if (abs(mcParticle.pdgCode()) != 313) {
        continue;
      }
      histos.fill(HIST("events_check"), 6.5);
      auto kDaughters = mcParticle.daughters_as<aod::McParticles>();
      if (kDaughters.size() != 2) {
        continue;
      }
      histos.fill(HIST("events_check"), 7.5);
      auto passkaon = false;
      auto passpion = false;
      for (auto kCurrentDaughter : kDaughters) {
        if (!kCurrentDaughter.isPhysicalPrimary()) {
          continue;
        }
        histos.fill(HIST("events_check"), 8.5);
        if (abs(kCurrentDaughter.pdgCode()) == 321) {
          // if (kCurrentDaughter.pdgCode() == +321) {
          passkaon = true;
          histos.fill(HIST("events_check"), 9.5);
        } else if (abs(kCurrentDaughter.pdgCode()) == 211) {
          //} else if (kCurrentDaughter.pdgCode() == -321) {
          passpion = true;
          // histos.fill(HIST("events_check"), 10.5);
        }
      }
      if (passkaon && passpion) {
        // if (mcParticle.pdgCode() > 0)
        histos.fill(HIST("k892Gen"), mcParticle.pt());
        // else
        //   histos.fill(HIST("k892GenAnti"), mcParticle.pt());
      }
    }
  }
  PROCESS_SWITCH(kstarqa, processGen, "Process Generated", false);

  void processRec(EventCandidatesMC::iterator const& collision, TrackCandidatesMC const& tracks, aod::McParticles const&, aod::McCollisions const& /*mcCollisions*/)
  {

    TLorentzVector lDecayDaughter1, lDecayDaughter2, lResonance;

    histos.fill(HIST("events_checkrec"), 0.5);
    if (!collision.has_mcCollision()) {
      return;
    }
    histos.fill(HIST("events_checkrec"), 1.5);
    // if (std::abs(collision.mcCollision().posZ()) > cutzvertex || !collision.sel8()) {
    if (std::abs(collision.mcCollision().posZ()) > cutzvertex) {
      return;
    }
    histos.fill(HIST("events_checkrec"), 2.5);

    if (timFrameEvsel && !collision.selection_bit(aod::evsel::kNoTimeFrameBorder)) {
      return;
    }
    histos.fill(HIST("events_checkrec"), 3.5);
    if (TVXEvsel && (!collision.selection_bit(aod::evsel::kIsTriggerTVX))) {
      return;
    }
    histos.fill(HIST("events_checkrec"), 4.5);

    auto oldindex = -999;
    for (auto track1 : tracks) {
      if (!selectionTrack(track1)) {
        continue;
      }
      histos.fill(HIST("events_checkrec"), 5.5);
      if (!track1.has_mcParticle()) {
        continue;
      }
      histos.fill(HIST("events_checkrec"), 6.5);
      auto track1ID = track1.index();
      for (auto track2 : tracks) {
        if (!track2.has_mcParticle()) {
          continue;
        }
        histos.fill(HIST("events_checkrec"), 7.5);
        if (!selectionTrack(track2)) {
          continue;
        }
        histos.fill(HIST("events_checkrec"), 8.5);
        auto track2ID = track2.index();
        if (track2ID <= track1ID) {
          continue;
        }
        // if (!selectionPair(track1, track2)) {
        //   continue;
        // }
        histos.fill(HIST("events_checkrec"), 9.5);
        if (track1.sign() * track2.sign() >= 0) {
          continue;
        }
        histos.fill(HIST("events_checkrec"), 10.5);
        const auto mctrack1 = track1.mcParticle();
        const auto mctrack2 = track2.mcParticle();
        int track1PDG = std::abs(mctrack1.pdgCode());
        int track2PDG = std::abs(mctrack2.pdgCode());
        if (!mctrack1.isPhysicalPrimary()) {
          continue;
        }
        histos.fill(HIST("events_checkrec"), 11.5);
        if (!mctrack2.isPhysicalPrimary()) {
          continue;
        }
        histos.fill(HIST("events_checkrec"), 12.5);
        if (!(selectionPID(track1, 1) && selectionPID(track2, 0))) { // pion and kaon
          continue;
        }
        histos.fill(HIST("events_checkrec"), 13.5);
        if (!(track1PDG == 321 && track2PDG == 211)) {
          continue;
        }
        // LOG(info) << "trackpdgs are:"<<track1PDG<<" "<<track2PDG;
        histos.fill(HIST("events_checkrec"), 14.5);
        for (auto& mothertrack1 : mctrack1.mothers_as<aod::McParticles>()) {
          for (auto& mothertrack2 : mctrack2.mothers_as<aod::McParticles>()) {
            if (mothertrack1.pdgCode() != mothertrack2.pdgCode()) {
              continue;
            }
            histos.fill(HIST("events_checkrec"), 15.5);
            if (mothertrack1.globalIndex() != mothertrack2.globalIndex()) {
              continue;
            }
            histos.fill(HIST("events_checkrec"), 16.5);
            if (!mothertrack1.producedByGenerator()) {
              continue;
            }
            histos.fill(HIST("events_checkrec"), 17.5);
            if (std::abs(mothertrack1.y()) >= 0.5) {
              continue;
            }
            histos.fill(HIST("events_checkrec"), 18.5);
            if (std::abs(mothertrack1.pdgCode()) != 313) {
              continue;
            }

            if (avoidsplitrackMC && oldindex == mothertrack1.globalIndex()) {
              histos.fill(HIST("h1KSRecsplit"), mothertrack1.pt());
              continue;
            }
            oldindex = mothertrack1.globalIndex();
            pvec0 = array{track1.px(), track1.py(), track1.pz()};
            pvec1 = array{track2.px(), track2.py(), track2.pz()};
            auto arrMomrec = array{pvec0, pvec1};
            auto motherP = mothertrack1.p();
            auto motherE = mothertrack1.e();
            auto genMass = std::sqrt(motherE * motherE - motherP * motherP);
            auto recMass = RecoDecay::m(arrMomrec, array{massKa, massPi});
            // auto recpt = TMath::Sqrt((track1.px() + track2.px()) * (track1.px() + track2.px()) + (track1.py() + track2.py()) * (track1.py() + track2.py()));
            //// Resonance reconstruction
            lDecayDaughter1.SetXYZM(track1.px(), track1.py(), track1.pz(), massKa);
            lDecayDaughter2.SetXYZM(track2.px(), track2.py(), track2.pz(), massPi);
            lResonance = lDecayDaughter1 + lDecayDaughter2;
            histos.fill(HIST("h3KstarRec"), motherP);
            histos.fill(HIST("h1KstarRecMass"), recMass);
            histos.fill(HIST("h1KstarRecpt"), mothertrack1.pt());
            histos.fill(HIST("h1genmass"), genMass);
            histos.fill(HIST("h1recpt"), lResonance.Pt());
          }
        }
      }
    }
  }
  PROCESS_SWITCH(kstarqa, processRec, "Process Reconstructed", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<kstarqa>(cfgc)};
}
