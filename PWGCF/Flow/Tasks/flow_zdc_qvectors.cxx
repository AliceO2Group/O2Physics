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

// In this task the energy calibration and recentring of Q-vectors constructed in the ZDCs will be done 
// Start with step one and add other steps later.

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/CCDB/EventSelectionParams.h"
#include "Common/CCDB/TriggerAliases.h"
#include "CCDB/BasicCCDBManager.h"
#include <DataFormatsParameters/GRPObject.h>
#include <DataFormatsParameters/GRPMagField.h>
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/Multiplicity.h"
#include "ReconstructionDataFormats/GlobalTrackID.h"
#include "ReconstructionDataFormats/Track.h"

#include <cmath>
#include <algorithm>
#include <numeric>
#include <vector>
#include "Framework/ASoAHelpers.h"
#include "Framework/RunningWorkflowInfo.h"

#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/Centrality.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TObjArray.h"
#include <stdlib.h>
#include "TF1.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TSystem.h"
#include <TROOT.h>

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::aod::track;
using namespace o2::aod::evsel;

// define my.....
using myCollisions = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::CentFT0Cs>>;
using myTracks = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TrackSelection, aod::TracksDCA>>;
using BCsRun3 = soa::Join<aod::BCs, aod::Timestamps, aod::BcSels, aod::Run3MatchedToBCSparse>;

namespace o2::analysis::qvectortask
{

int counter = 0;

// step0 -> Energy calib
std::shared_ptr<TProfile2D> ZN_Energy[10] = {{nullptr}}; // fill for weights
std::shared_ptr<TH2> hQx_vs_Qy[6] = {{nullptr}};         // fill Qx vs Qy for each step in the recentring process.

// <XX> <XY> <YX> and <YY> for step 0 (and step 5)
std::shared_ptr<TProfile> COORD_correlations[2][4] = {{nullptr}};

std::vector<const char*> sides = {"A", "C"};
std::vector<const char*> coords = {"x", "y", "z"};
std::vector<const char*> COORDS = {"X", "Y"};


// Define histogrm names here to use same names for creating and later uploading and retrieving data from ccdb 
// Eneergy calibration: 
std::vector<TString> names_Ecal(10,"");


// https://alice-notes.web.cern.ch/system/files/notes/analysis/620/017-May-31-analysis_note-ALICE_analysis_note_v2.pdf
std::vector<double> ZDC_px = {-1.75, 1.75, -1.75, 1.75};
std::vector<double> ZDC_py = {-1.75, -1.75, 1.75, 1.75};
double alphaZDC = 0.395; // Dit is oud PAS OP!! 

// step 0 tm 5 A&C
std::vector<std::vector<double>> q(6,std::vector<double>(4,0));   // 6 steps, each with 4 values
std::vector<double> v(3);                                         // vx, vy, vz
std::vector<double> v_mean(3);                                    // vx_mean, vy_mean, vz_mean

// for energy calibration
std::vector<double> EZN(8);      // uncalibrated energy for the 2x4 towers (a1, a2, a3, a4, c1, c2, c3, c4)
std::vector<double> meanEZN(10); // mean energies from calibration histos (common A, t1-4 A,common C, t1-4C)
std::vector<double> e(8,0.);     // calibrated energies (a1, a2, a3, a4, c1, c2, c3, c4))

} // namespace o2::analysis::qvectortask

using namespace o2::analysis::qvectortask;

struct ZDCqvectors {
  ConfigurableAxis axisCent{"axisCent", {90, 0, 90}, "Centrality axis in 1% bins"};
  ConfigurableAxis axisCent10{"axisCent10", {9, 0, 90}, "Centrality axis in 10% bins"};
  ConfigurableAxis axisQ{"axisQ", {100, -2, 2}, "Q vector (xy) in ZDC"};
  ConfigurableAxis axisVx{"axisVx", {100, 0, 0.15}, "for Pos X of collision"};
  ConfigurableAxis axisVy{"axisVy", {100, 0.35, 0.41}, "for Pos Y of collision"};
  ConfigurableAxis axisVx_mean{"axisVx_mean", {100, -0.03, 0.03}, "for Pos X of collision"};
  ConfigurableAxis axisVy_mean{"axisVy_mean", {100, -0.03, 0.03}, "for Pos Y of collision"};
  ConfigurableAxis axisVz{"axisVz", {100, -12, 12}, "for vz of collision"}; // take 12 because we shift vi - <vi>
  ConfigurableAxis axisRun{"axisRun", {1e6, 0, 1e6}, "for runNumber in ThnSparse"};
  ConfigurableAxis axisPolarity{"axisPolarity", {2, -1, 1}, "Magnet Polarity"};
  ConfigurableAxis axisMult{"axisMult", {5000, 0, 5000}, "Track Multiplicity"};
  ConfigurableAxis axisMult_mean{"axisMult_mean", {5000, -2500, 2500}, "Track Multiplicity mean per run"};

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 10.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMin, float, 0.2f, "Minimal.q pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtPOIMax, float, 10.0f, "Maximal pT for poi tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMin, float, 0.2f, "Minimal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutPtMax, float, 3.0f, "Maximal pT for ref tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutChi2prTPCcls, float, 2.5, "Chi2 per TPC clusters")
  O2_DEFINE_CONFIGURABLE(cfgMagField, float, 99999, "Configurable magnetic field; default CCDB will be queried")
  O2_DEFINE_CONFIGURABLE(cfgEnergyCal, std::string, "", "ccdb path for energy calibration histos")

  //  Define variables needed to do the recentring steps.
  double centrality;
  double trackmultiplicity;
  double trackmultiplicity_mean;
  int runnumber;
  double polarity;

  //  Define output
  HistogramRegistry registry{"Registry"};

  //  Filters
  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter trackFilter = (nabs(aod::track::eta) < cfgCutEta) && (aod::track::pt > cfgCutPtMin) && (aod::track::pt < cfgCutPtMax) && ((requireGlobalTrackInFilter()) || (aod::track::isGlobalTrackSDD == (uint8_t) true)) && (aod::track::tpcChi2NCl < cfgCutChi2prTPCcls);

  Service<ccdb::BasicCCDBManager> ccdb;

  // keep calibration histos for each given step 
  struct Calib {
      std::vector<TList*> calibList = std::vector<TList*>(8, nullptr);
      std::vector<bool> calibfilesLoaded = std::vector<bool>(8, false);
  } cal;

  void init(InitContext const&)
  {
    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    int64_t now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);

    // Tower mean energies vs. centrality used for tower gain equalisation
    for (int tower = 0; tower < 10; tower++) {
      names_Ecal[tower] = TString::Format("hZN%s_mean_t%i_cent", sides[(tower < 5) ? 0 : 1], tower % 5);
      ZN_Energy[tower]  = registry.add<TProfile2D>(Form("Energy/%s", names_Ecal[tower].Data()), Form("%s", names_Ecal[tower].Data()), kTProfile2D, {{1, 0, 1}, axisCent});
    }
    
    // Qx_vs_Qy before recentering for ZNA and ZNC
    int step = 0;
    for (const char* side : sides) {
    hQx_vs_Qy[step] = registry.add<TH2>(Form("step%i/hZN%s_Qx_vs_Qy", step, side), Form("hZN%s_Qx_vs_Qy", side), kTH2F, {axisQ, axisQ});
    }
    for (const char* COORD1 : COORDS) {
        for (const char* COORD2 : COORDS) {
        // Now we get: <XX> <XY> & <YX> <YY> vs. Centrality : do only before (step0) and after (step5) of recentring.
        COORD_correlations[step][0] = registry.add<TProfile>(Form("step%i/hQ%sA_Q%sC_vs_cent", step, COORD1, COORD2), Form("hQ%sA_Q%sC_vs_cent", COORD1, COORD2), kTProfile, {axisCent10});
        }
    }

    // recentered q-vectors (to check what steps are finished in the end)
    registry.add("hStep", "hStep", {HistType::kTH1D, {{6, 0., 6.}}});
  }

  inline void fillRegistry(int step, double qxa, double qya, double qxc, double qyc)
  {
    if (std::isnan(qxa) || std::isnan(qya) || std::isnan(qxc) || std::isnan(qyc)) {
      return;
    }

    if (step == 0) {
      registry.get<TProfile>(HIST("step0/hQXA_QXC_vs_cent"))->Fill(centrality, qxa * qxc);
      registry.get<TProfile>(HIST("step0/hQYA_QYC_vs_cent"))->Fill(centrality, qya * qyc);
      registry.get<TProfile>(HIST("step0/hQYA_QXC_vs_cent"))->Fill(centrality, qxc * qya);
      registry.get<TProfile>(HIST("step0/hQXA_QYC_vs_cent"))->Fill(centrality, qyc * qxa);

      registry.fill(HIST("step0/hZNA_Qx_vs_Qy"), qxa, qya);
      registry.fill(HIST("step0/hZNC_Qx_vs_Qy"), qxc, qyc);

      registry.fill(HIST("hStep"), .5, 1);
    }
  }


  
  void loadCalibrations(int step, uint64_t timestamp)
  {
    std::string ccdb_dir = cfgEnergyCal;

    // for now return if step!=0 because only energy calibration is implemented
    if(step!=0) return; 

    if (cfgEnergyCal.value.empty() == false){ //only proceed if string is not empty
      cal.calibList[step] = ccdb->getForTimeStamp<TList>(ccdb_dir, timestamp);

      if (cal.calibList[step]){
        for (int i = 0; i < 10; i++) {
          TProfile2D* vec = nullptr;  
          vec = (TProfile2D*)cal.calibList[step]->FindObject(Form("%s", names_Ecal[i].Data()));
          if (vec == nullptr || vec->GetEntries() < 1) {
            LOGF(info, "%s not found or empty! Produce calibration file at given step", names_Ecal[i].Data());
            cal.calibfilesLoaded[step] =  false;
            break; 
          }
        }
          cal.calibfilesLoaded[step] = true; 
          } else {
          LOGF(warning, "Could not load TList with calibration histos from %s", ccdb_dir);
          cal.calibfilesLoaded[step] = false; }
    }
  }

  int getMagneticField(uint64_t timestamp)
  {
    static o2::parameters::GRPMagField* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>("GLO/Config/GRPMagField", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  
  double getCorrectionStep0(int runnumber, double centrality, int tower){ 

    TProfile2D* hist = nullptr;  

    hist = (TProfile2D*)cal.calibList[0]->FindObject(Form("%s",names_Ecal[tower].Data()));
    if(!hist) {
      LOGF(fatal, "Histo with calib. conatants for tower %i not available.. Abort..", tower);
    }
    int binrunnumber = hist->GetXaxis()->FindBin(TString::Format("%i",runnumber)); 
    if(binrunnumber!=1) LOGF(fatal, "bin runnumber is wrong.."); 
    double calibConstant = hist->GetBinContent(int(binrunnumber), int(centrality) + 1);

    return calibConstant; 
  }


  void process(myCollisions::iterator const& collision,
               BCsRun3 const& /*bcs*/,
               aod::Zdcs const& /*zdcs*/,
               myTracks const& tracks)
  {

    // for Q-vector calculation
    //  A[0] & C[1]
    std::vector<double> sumZN(2,0.);  
    std::vector<double> xEnZN(2,0.); 
    std::vector<double> yEnZN(2,0.); 

    if (!collision.sel8())
      return;
    auto cent = collision.centFT0C();
    if (cent < 0 || cent > 90)
      return;

    const auto& foundBC = collision.foundBC_as<BCsRun3>();
    auto field = (cfgMagField == 999999) ? getMagneticField(foundBC.timestamp()) : cfgMagField;
    if (foundBC.has_zdc()) {

      v[0] = collision.posX();
      v[1] = collision.posY();
      v[2] = collision.posZ();
      centrality = cent;
      trackmultiplicity = tracks.size();
      polarity = (field < 0) ? -.5 : .5;

      int runnumber = foundBC.runNumber();
      const auto& zdcCol = foundBC.zdc();

      // Get the raw energies EZN[8] (not the common A,C)
      for (int tower = 0; tower < 4; tower++) {
        EZN[tower] = zdcCol.energySectorZNA()[tower];
        EZN[tower + 4] = zdcCol.energySectorZNC()[tower];
      }

      // load the calibration histos for step 0 (Energy Calibration)
      loadCalibrations(0, foundBC.timestamp()); 

      if (!cal.calibfilesLoaded[0]) {
        if (counter < 1) {
          LOGF(info, "files to start recentering not found.. ");
          LOGF(info, "=====================> .....Start Energy Calibration..... <=====================");
        }
      }

      bool isZNAhit = true; 
      bool isZNChit = true; 

      for (int i = 0; i < 8; ++i) {
        if (i<4 && EZN[i] <= 0) isZNAhit = false;
        if (i>3 && EZN[i] <= 0) isZNChit = false;
      }
      
      if(zdcCol.energyCommonZNA() <= 0 ) isZNAhit = false; 
      if(zdcCol.energyCommonZNC() <= 0 ) isZNChit = false; 


      // Fill to get mean energy per tower in 1% centrality bins
      for (int tower = 0; tower < 5; tower++) {
        if (tower == 0) {
          if (isZNAhit)
            ZN_Energy[tower]->Fill(Form("%d", runnumber), cent, zdcCol.energyCommonZNA(), 1);
          if (isZNChit)
            ZN_Energy[tower + 5]->Fill(Form("%d", runnumber), cent, zdcCol.energyCommonZNC(), 1);
          LOGF(debug, "Common A tower filled with: %i, %.2f, %.2f", runnumber, cent,zdcCol.energyCommonZNA() );
        } else {
          if (isZNAhit)
            ZN_Energy[tower]->Fill(Form("%d", runnumber), cent, EZN[tower - 1], 1);
          if (isZNChit)
            ZN_Energy[tower + 5]->Fill(Form("%d", runnumber), cent, EZN[tower-1+4], 1);
          LOGF(debug, "Tower ZNC[%i] filled with: %i, %.2f, %.2f", tower, runnumber, cent,EZN[tower-1+4]);
        }
      }

      if(!isZNAhit || !isZNChit) {
        LOGF(info, "ZDC not hit (correctly) ... don't use event!"); 
        counter++;
        return; 
      }

      if (cal.calibfilesLoaded[0]) {
        if (counter < 1)
          LOGF(info, "files for step 0 (energy Calibraton) are open!");

        if (counter < 1) {
          LOGF(info, "=====================> .....Start Calculating Q-Vectors..... <=====================");
        }

        // Now start gain equalisation!
        //Fill the list with calibration constants. 
        for (int tower = 0; tower < 10; tower++) {
          meanEZN[tower] = getCorrectionStep0(runnumber, centrality, tower);
          }

       
        // Use the calibration constants but now only loop over towers 1-4
        int calibtower=0; 
        std::vector<int> towers_nocom = {1,2,3,4,6,7,8,9}; 

        for (int tower : towers_nocom) { 
          if (meanEZN[tower] > 0) {
            double ecommon = (tower>4) ? meanEZN[5] : meanEZN[0]; 
            e[calibtower] = EZN[calibtower] * (0.25 * ecommon) / meanEZN[tower];
          }
          calibtower++; 
        }

        // Now calculate Q-vector
        for (int tower = 0; tower < 8; tower++) {
          int side = (tower > 3) ? 1 : 0; 
          int sector = tower % 4;
          double energy = TMath::Power(e[tower], alphaZDC);
          sumZN[side] += energy; 
          xEnZN[side] += ZDC_px[sector] * energy; 
          yEnZN[side] += ZDC_py[sector] * energy;
        }

        // "QXA", "QYA", "QXC", "QYC"
        for (int i = 0; i < 2; ++i) {
          if (sumZN[i] > 0) {
            q[0][i * 2] = xEnZN[i] / sumZN[i];     // for QXA[0] and QXC[2]
            q[0][i * 2 + 1] = yEnZN[i] / sumZN[i]; // for QYA[1] and QYC[3]
          }
        }     
        

        if (!cal.calibfilesLoaded[1]) {
          if (counter < 1)
            LOGF(warning, "Calibation files for recentring does not exist (yet) | Output created with non-recentered q-vectors!!!!");
          // LOGF(info,"Start filling: "); 
          // for (double vecc : q[0]) LOGF(info,"%.2f", vecc); 
          fillRegistry(0, q[0][0], q[0][1], q[0][2], q[0][3]); 
          counter++;

          return;
        }
      } // end of cal.calibfilesLoaded[0]

    } // end collision found ZDC
    counter++;
  } // end of process
  
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ZDCqvectors>(cfgc)};
}
