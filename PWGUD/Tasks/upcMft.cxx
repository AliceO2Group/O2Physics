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

// c++ headers
#include <iostream>

// root headers
#include <TH1F.h>
#include <TH2F.h>
#include "TMath.h"
#include "TLorentzVector.h"

// framework headers
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

// for track propagation
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/MatrixFunctions.h"
#include "Math/SMatrix.h"
#include "MathUtils/Utils.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DetectorsBase/Propagator.h"
#include "TGeoGlobalMagField.h"
#include "Field/MagneticField.h"
#include "DataFormatsMFT/TrackMFT.h"

// define namespaces
using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// namespace for track propagation

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

using ExtBCs = soa::Join<aod::BCs, aod::Timestamps>;

struct BcTracks {
  std::vector<aod::AmbiguousMFTTrack> tracks;
  float timeFT0C = -999; // Default is -999, because FT0 time can be negative, so the usual -1 is not ideal.
  float timeFT0A = -999;
  float ZposFT0 = -999; // Z position calculated from timeA and timeC in cm
  int nChannelsFT0C = 0;
  std::vector<o2::track::TrackParCovFwd> trackParams = {};
};

//-----------------------------------------------------------------------------------------
struct UpcMftRec {

  // defining histograms using histogram registry
  HistogramRegistry registry{"registry", {}, OutputObjHandlingPolicy::AnalysisObject};

  // tools for track propagation
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  int runNumber = -1;
  float Bz = 0;                                         // Magnetic field for MFT
  static constexpr double centerMFT[3] = {0, 0, -61.4}; // Field at center of MFT
  o2::parameters::GRPMagField* grpmag = nullptr;

  // initialisation
  // -- define histograms
  void init(o2::framework::InitContext&)
  {
    // to access CCDB
    ccdb->setURL(ccdburl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    // histograms for the reconstructed rho
    registry.add("hMassRhoRec", "Mass of reconstructed #rho^{0}; mass (GeV/c^{2}); Entries", kTH1F, {{150, 0., 1.5}});
    registry.add("hMassRhoRec4", "Mass of reconstructed #rho (4 tracks); mass (GeV/c^{2}); Entries", kTH1F, {{200, 0., 2.0}});
    registry.add("hPtRhoRec", "p_{T} of reconstructed #rho^{0};#it{p}_{T} (GeV/c); Entries", kTH1F, {{200, 0., 1.0}});
    registry.add("hPtRhoRec4", "p_{T} of reconstructed #rho (4 tracks);#it{p}_{T} (GeV/c); Entries", kTH1F, {{400, 0., 2.0}});
    registry.add("hMassPtRhoRec", "Mass vs. p_{T};Mass of reconstructed #rho^{0} (GeV/c^{2});p_{T} of reconstructed #rho^{0} (GeV/c)", {HistType::kTH2F, {{100, 0., 2.0}, {200, 0., 2.0}}});
    registry.add("hMassPtRhoRec4", "Mass vs. p_{T} (4 tracks);Mass of reconstructed #rho^{0} (GeV/c^{2});p_{T} of reconstructed #rho^{0} (GeV/c)", {HistType::kTH2F, {{100, 0., 2.0}, {200, 0., 2.0}}});
    registry.add("hMassPtRhoRecElectrons", "#gamma#gamma#rightarrow#it{e}^{+}#it{e}^{-} (2 tracks);mass (GeV/c^{2}) ;#it{p}_{T} (GeV/c)", {HistType::kTH2F, {{150, 0., 1.5}, {200, 0., 1.0}}});
    registry.add("hTimeFT0C", "FT0C time (ns); Entries", kTH1F, {{200, -1.0, 1.0}});
    registry.add("hTimeFT0C4", "FT0C time (ns) (4 tracks); Entries", kTH1F, {{200, -1.0, 1.0}});
    registry.add("hFT0CSignalsPerBC", "Number of FT0C channels with signal per BC; Entries", kTH1F, {{112, 0.5, 112.5}});
    registry.add("hTracksPerBC", "Tracks per BC; Tracks per BC; Entries", kTH1F, {{10000, 0, 10000}});
    registry.add("hChargeDistribution4Tracks", "Charge distribution of 4-track events; q_{1}+q_{2}+q_{3}+q_{4}; Entries", kTH1F, {{9, -4.5, 4.5}});
    // control plots for 2-track events
    registry.add("ControlHistos/hTrackChi2", "Track #chi^{2}; #chi^{2}; Entries", kTH1F, {{210, -0.5, 20.5}});
    registry.add("ControlHistos/hTrackNClusters", "Number of clusters per track; Number of clusters; Entries", kTH1F, {{10, 0.5, 10.5}});
    registry.add("ControlHistos/hTrackPt", "Track p_{T}; p_{T} (GeV/c); Entries", kTH1F, {{200, 0, 1.0}});
    registry.add("ControlHistos/hTrackEta", "Track #eta; #eta; Entries", kTH1F, {{50, -4, -2}});
    registry.add("ControlHistos/hTrackEtaPos", "Positive track #eta; #eta; Entries", kTH1F, {{50, -4, -2}});
    registry.add("ControlHistos/hTrackEtaNeg", "Negative track #eta; #eta; Entries", kTH1F, {{50, -4, -2}});
    registry.add("ControlHistos/hTrackPhi", "Track #phi; #phi; Entries", kTH1F, {{100, -3.2, 3.2}});
    registry.add("ControlHistos/hTrackPhiPos", "Positive track #phi; #phi; Entries", kTH1F, {{100, -3.2, 3.2}});
    registry.add("ControlHistos/hTrackPhiNeg", "Negative track #phi; #phi; Entries", kTH1F, {{100, -3.2, 3.2}});
    // propagated tracks control plots
    registry.add("PropControlHistos/hTrackPt", "Track p_{T}; p_{T} (GeV/c); Entries", kTH1F, {{200, 0, 1.0}});
    registry.add("PropControlHistos/hTrackEta", "Track #eta; #eta; Entries", kTH1F, {{50, -4, -2}});
    registry.add("PropControlHistos/hTrackEta2D", "Pseudorapidity of propagated track pairs; #eta^{+}; #eta^{-}", {HistType::kTH2F, {{50, -4, -2}, {50, -4, -2}}});
    registry.add("PropControlHistos/hTrackPhi", "Track #phi; #phi; Entries", kTH1F, {{100, -3.2, 3.2}});
    registry.add("PropControlHistos/hTrackPhi2D", "#phi of propagated track pairs; #phi^{+}; #phi^{-}", {HistType::kTH2F, {{100, -3.2, 3.2}, {100, -3.2, 3.2}}});
  }

  void initCCDB(ExtBCs::iterator const& bc)
  {
    if (runNumber == bc.runNumber()) {
      return;
    }
    grpmag = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(grpmagPath, bc.timestamp());
    LOG(info) << "Setting magnetic field to current " << grpmag->getL3Current()
              << " A for run " << bc.runNumber()
              << " from its GRPMagField CCDB object";
    o2::base::Propagator::initFieldFromGRP(grpmag); // for some reason this is necessary for the next next line
    runNumber = bc.runNumber();

    o2::field::MagneticField* field = static_cast<o2::field::MagneticField*>(TGeoGlobalMagField::Instance()->GetField());

    Bz = field->getBz(centerMFT); // gives error if the propagator is not initFielded
    LOG(info) << "The field at the center of the MFT is Bz = " << Bz;
  }

  void process(aod::BCs const& bcs,
               aod::Collisions const& /*collisions*/,
               aod::MFTTracks const& /*tracksMFT*/,
               aod::AmbiguousMFTTracks const& ambiguousTracksMFT,
               aod::FT0s const& ft0s,
               ExtBCs const& ebcs)
  {

    std::map<int, BcTracks> bcGroups; // Map of tracks that belong to a given bunch crossing. One track can be associated with multiple BCs.

    initCCDB(ebcs.begin()); // Bz is needed for propagation

    for (auto& ambiguousTrackMFT : ambiguousTracksMFT) {
      auto mftBCs = ambiguousTrackMFT.bc();
      for (auto& mftBC : mftBCs) {
        int bcId = mftBC.globalBC();        // Get global BC that will be used for grouping.
        auto element = bcGroups.find(bcId); // Find the iterator poiniting to value associated with this BCid.
        if (element == bcGroups.end()) {    // If this BCid was encountered for the first time, create a vector containing the current track
          BcTracks bct;
          bct.tracks.push_back(ambiguousTrackMFT); // was ambiguousTrackMFT

          auto track = ambiguousTrackMFT.mfttrack();
          std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
          SMatrix55 tcovs(v1.begin(), v1.end());
          SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
          o2::track::TrackParCovFwd trackPar{track.z(), tpars, tcovs, track.chi2()};
          bct.trackParams.push_back(trackPar);
          bcGroups.emplace(std::make_pair(bcId, bct));
        } else { // If this BcId was already encountered, add the track to the existing vector
          element->second.tracks.push_back(ambiguousTrackMFT);
          auto track = ambiguousTrackMFT.mfttrack();
          std::vector<double> v1; // Temporary null vector for the computation of the covariance matrix
          SMatrix55 tcovs(v1.begin(), v1.end());
          SMatrix5 tpars(track.x(), track.y(), track.phi(), track.tgl(), track.signed1Pt());
          o2::track::TrackParCovFwd trackPar{track.z(), tpars, tcovs, track.chi2()};
          element->second.trackParams.push_back(trackPar);
        }
      }

      for (auto& ft0 : ft0s) {
        int collisionBC = (bcs.iteratorAt(ft0.bcId())).globalBC(); // Get global BC that will be used for grouping.
        auto element = bcGroups.find(collisionBC);                 // Find the iterator poiniting to value associated with this BCid.
        if (element == bcGroups.end()) {
          // std::cout << " There is no global BC id for FT0. " << std::endl;
        } else {
          element->second.timeFT0C = ft0.timeC();
          element->second.timeFT0A = ft0.timeA();
          element->second.ZposFT0 = ft0.posZ();
          element->second.nChannelsFT0C = ft0.channelC().size();
        }
      }
    }
    BcTracks defaultBcTrack; // Struct with default values used later to filter out default FT0 time
    for (auto& bcGroup : bcGroups) {
      registry.fill(HIST("hTracksPerBC"), bcGroup.second.tracks.size());

      if (bcGroup.second.timeFT0C != defaultBcTrack.timeFT0C && std::abs(bcGroup.second.timeFT0C) < 1 && std::abs(bcGroup.second.timeFT0A) > 5 && bcGroup.second.tracks.size() == 2) {

        // Fill control histos for non-propagated tracks
        for (int i = 0; i < 2; i++) {
          auto mftTrack = bcGroup.second.tracks[i].mfttrack(); // Get MFT tracks for the two-track case
          registry.fill(HIST("ControlHistos/hTrackChi2"), mftTrack.chi2());
          registry.fill(HIST("ControlHistos/hTrackNClusters"), mftTrack.nClusters());
          registry.fill(HIST("ControlHistos/hTrackPt"), mftTrack.pt());
          registry.fill(HIST("ControlHistos/hTrackEta"), mftTrack.eta());
          registry.fill(HIST("ControlHistos/hTrackPhi"), mftTrack.phi());
          if (mftTrack.sign() > 0) {
            registry.fill(HIST("ControlHistos/hTrackEtaPos"), mftTrack.eta());
            registry.fill(HIST("ControlHistos/hTrackPhiPos"), mftTrack.phi());
          }
          if (mftTrack.sign() < 0) {
            registry.fill(HIST("ControlHistos/hTrackEtaNeg"), mftTrack.eta());
            registry.fill(HIST("ControlHistos/hTrackPhiNeg"), mftTrack.phi());
          }
        }

        // Get propagated MFT tracks (two-track case)
        auto mftPropTrack1 = bcGroup.second.trackParams[0];
        auto mftPropTrack2 = bcGroup.second.trackParams[1];

        mftPropTrack1.propagateToZhelix(bcGroup.second.ZposFT0, Bz); // track parameters propagation to the z position given by FT0
        mftPropTrack2.propagateToZhelix(bcGroup.second.ZposFT0, Bz); // track parameters propagation to the z position given by FT0

        // Fill control histos for the propagated tracks
        for (int i = 0; i < 2; i++) {
          auto mftPropTrack = bcGroup.second.trackParams[i];
          registry.fill(HIST("PropControlHistos/hTrackPt"), mftPropTrack.getPt());
          registry.fill(HIST("PropControlHistos/hTrackEta"), mftPropTrack.getEta());
          registry.fill(HIST("PropControlHistos/hTrackPhi"), mftPropTrack.getPhi()); // getCharge() is below
        }

        // Make a rho candidate ( using propagated tracks )
        if (mftPropTrack1.getCharge() + mftPropTrack2.getCharge() != 0)
          continue; // Skip if total charge is not 0 ; if using non propagated tracks use sign() instead of getCharge()
        if (mftPropTrack1.getCharge() > 0 && mftPropTrack2.getCharge() < 0) {
          registry.fill(HIST("PropControlHistos/hTrackEta2D"), mftPropTrack1.getEta(), mftPropTrack2.getEta());
          registry.fill(HIST("PropControlHistos/hTrackPhi2D"), mftPropTrack1.getPhi(), mftPropTrack2.getPhi());
        }
        if (mftPropTrack1.getCharge() < 0 && mftPropTrack2.getCharge() > 0) {
          registry.fill(HIST("PropControlHistos/hTrackEta2D"), mftPropTrack2.getEta(), mftPropTrack1.getEta());
          registry.fill(HIST("PropControlHistos/hTrackPhi2D"), mftPropTrack2.getPhi(), mftPropTrack1.getPhi());
        }

        TLorentzVector pi1LV;
        TLorentzVector pi2LV;
        pi1LV.SetPtEtaPhiM(mftPropTrack1.getPt(), mftPropTrack1.getEta(), mftPropTrack1.getPhi(), o2::constants::physics::MassPionCharged); // if using non propagated tracks use pt(),eta(),phi() instead
        pi2LV.SetPtEtaPhiM(mftPropTrack2.getPt(), mftPropTrack2.getEta(), mftPropTrack2.getPhi(), o2::constants::physics::MassPionCharged); // if using non propagated tracks use pt(),eta(),phi() instead
        TLorentzVector rhoLV = pi1LV + pi2LV;
        // Fill histos
        registry.fill(HIST("hMassRhoRec"), rhoLV.M());
        registry.fill(HIST("hPtRhoRec"), rhoLV.Pt());
        registry.fill(HIST("hMassPtRhoRec"), rhoLV.M(), rhoLV.Pt());
        registry.fill(HIST("hTimeFT0C"), bcGroup.second.timeFT0C);
        registry.fill(HIST("hFT0CSignalsPerBC"), bcGroup.second.nChannelsFT0C);
      }

      if (bcGroup.second.timeFT0C != defaultBcTrack.timeFT0C && std::abs(bcGroup.second.timeFT0C) < 1 && std::abs(bcGroup.second.timeFT0A) > 5 && bcGroup.second.tracks.size() == 4) { // 4 tracks events
        // Get tracks
        auto mftTrack1 = bcGroup.second.tracks[0].mfttrack();
        auto mftTrack2 = bcGroup.second.tracks[1].mfttrack();
        auto mftTrack3 = bcGroup.second.tracks[2].mfttrack();
        auto mftTrack4 = bcGroup.second.tracks[3].mfttrack();
        // Get propagated MFT tracks (four-track case)
        auto mftPropTrack1 = bcGroup.second.trackParams[0];
        auto mftPropTrack2 = bcGroup.second.trackParams[1];
        auto mftPropTrack3 = bcGroup.second.trackParams[2];
        auto mftPropTrack4 = bcGroup.second.trackParams[3];
        mftPropTrack1.propagateToZhelix(bcGroup.second.ZposFT0, Bz); // track parameters propagation to the z position given by FT0
        mftPropTrack2.propagateToZhelix(bcGroup.second.ZposFT0, Bz);
        mftPropTrack3.propagateToZhelix(bcGroup.second.ZposFT0, Bz); // track parameters propagation to the z position given by FT0
        mftPropTrack4.propagateToZhelix(bcGroup.second.ZposFT0, Bz); // track parameters propagation to the z position given by FT0
        // Make a rho candidate
        TLorentzVector pi1LV;
        TLorentzVector pi2LV;
        TLorentzVector pi3LV;
        TLorentzVector pi4LV;
        pi1LV.SetPtEtaPhiM(mftPropTrack1.getPt(), mftPropTrack1.getEta(), mftPropTrack1.getPhi(), o2::constants::physics::MassPionCharged);
        pi2LV.SetPtEtaPhiM(mftPropTrack2.getPt(), mftPropTrack2.getEta(), mftPropTrack2.getPhi(), o2::constants::physics::MassPionCharged);
        pi3LV.SetPtEtaPhiM(mftPropTrack3.getPt(), mftPropTrack3.getEta(), mftPropTrack3.getPhi(), o2::constants::physics::MassPionCharged);
        pi4LV.SetPtEtaPhiM(mftPropTrack4.getPt(), mftPropTrack4.getEta(), mftPropTrack4.getPhi(), o2::constants::physics::MassPionCharged);
        TLorentzVector rhoLV = pi1LV + pi2LV + pi3LV + pi4LV;
        // Fill histos
        registry.fill(HIST("hChargeDistribution4Tracks"), mftTrack1.sign() + mftTrack2.sign() + mftTrack3.sign() + mftTrack4.sign());
        if (mftTrack1.sign() + mftTrack2.sign() + mftTrack3.sign() + mftTrack4.sign() != 0)
          continue; // Skip if total charge is not 0
        registry.fill(HIST("hMassRhoRec4"), rhoLV.M());
        registry.fill(HIST("hPtRhoRec4"), rhoLV.Pt());
        registry.fill(HIST("hMassPtRhoRec4"), rhoLV.M(), rhoLV.Pt());
        registry.fill(HIST("hTimeFT0C4"), bcGroup.second.timeFT0C);
      }

      if (bcGroup.second.timeFT0C != defaultBcTrack.timeFT0C && std::abs(bcGroup.second.timeFT0C) < 1 && std::abs(bcGroup.second.timeFT0A) > 5 && bcGroup.second.tracks.size() == 2) { // electron mass hypothesis
        // Get tracks
        auto mftTrack1 = bcGroup.second.tracks[0].mfttrack();
        auto mftTrack2 = bcGroup.second.tracks[1].mfttrack();
        // Get propagated MFT tracks (two-track case)
        auto mftPropTrack1 = bcGroup.second.trackParams[0];
        auto mftPropTrack2 = bcGroup.second.trackParams[1];

        mftPropTrack1.propagateToZhelix(bcGroup.second.ZposFT0, Bz);
        mftPropTrack2.propagateToZhelix(bcGroup.second.ZposFT0, Bz);

        // Make a rho candidate
        if (mftTrack1.sign() + mftTrack2.sign() != 0)
          continue; // Skip if total charge is not 0
        TLorentzVector el1LV;
        TLorentzVector el2LV;
        el1LV.SetPtEtaPhiM(mftPropTrack1.getPt(), mftPropTrack1.getEta(), mftPropTrack1.getPhi(), o2::constants::physics::MassElectron);
        el2LV.SetPtEtaPhiM(mftPropTrack2.getPt(), mftPropTrack2.getEta(), mftPropTrack2.getPhi(), o2::constants::physics::MassElectron);
        TLorentzVector rhoLV = el1LV + el2LV;
        // Fill histos
        registry.fill(HIST("hMassPtRhoRecElectrons"), rhoLV.M(), rhoLV.Pt());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    o2::framework::adaptAnalysisTask<UpcMftRec>(cfgc)};
}
