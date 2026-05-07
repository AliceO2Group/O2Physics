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

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/TrackSelectionTables.h" // needed for aod::TracksDCA table

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/MathConstants.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>

#include <TCollection.h>
#include <TFile.h>
#include <TGrid.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

using namespace o2;
using namespace o2::framework;

// Definitions of join tables for Run 3 analysis:
using EventSelection = soa::Join<aod::EvSels, aod::Mults, aod::CentFT0As, aod::CentFT0Cs, aod::CentFT0Ms, aod::CentFV0As, aod::CentNTPVs>;
using CollisionRec = soa::Join<aod::Collisions, EventSelection>::iterator; // use in json "isMC": "true" for "event-selection-task"
using CollisionRecSim = soa::Join<aod::Collisions, aod::McCollisionLabels, EventSelection>::iterator;
using CollisionSim = aod::McCollision;
using TracksRec = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>;
using TrackRec = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>::iterator;
using TracksRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>; // + use in json "isMC" : "true"
using TrackRecSim = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::McTrackLabels>::iterator;
using TracksSim = aod::McParticles;
using TrackSim = aod::McParticles::iterator;

using namespace std;

// *) Define enums:
enum eRecSim { eRec = 0,
               eSim,
               eRecAndSim };

enum eProcess {
  eProcessRec = 0, // Run 3, only reconstructed
  eProcessRecSim,  // Run 3, both reconstructed and simulated
  eProcessSim,     // Run 3, only simulated
  eProcess_N
};

enum eEventHistograms {
  eVertexZ = 0,
  ePt,
  eEventHistograms_N
};

// *) Main task:
struct MultiharmonicCorrelations { // this name is used in lower-case format to name the TDirectoryFile in AnalysisResults.root
  Service<ccdb::BasicCCDBManager> ccdb;

  // *) Base TList to hold all output objects:
  TString sBaseListName = "Default list name";
  OutputObj<TList> fBaseList{sBaseListName.Data(),
                             OutputObjHandlingPolicy::AnalysisObject,
                             OutputObjSourceType::OutputObjSource};

  // *) Define configurables:
  Configurable<bool> cfDryRun{"cfDryRun", false, "book all histos and run without filling and calculating anything"}; // example for built-in type (float, string, etc.)
  Configurable<std::vector<float>> cf_pt_bins{"cf_pt_bins", {1000, 0., 100.}, "nPtBins, ptMin, ptMax"};               // example for an array
  Configurable<std::vector<float>> cf_phi_bins{"cf_phi_bins", {100, 0., o2::constants::math::TwoPI}, "nPhiBins, phiMin, phiMax"};
  Configurable<std::vector<float>> cf_centr_bins{"cf_centr_bins", {1000, 0., 100.}, "nCentrBins, centrMin, centrMax"};
  Configurable<std::vector<float>> cf_x_bins{"cf_x_bins", {1000, -100., 100.}, "nXBins, xMin, xMax"};
  Configurable<std::vector<float>> cf_y_bins{"cf_y_bins", {1000, -100., 100.}, "nYBins, yMin, yMax"};
  Configurable<std::vector<float>> cf_z_bins{"cf_z_bins", {1000, -100., 100.}, "nZBins, zMin, zMax"};
  Configurable<std::vector<float>> cf_mult_bins{"cf_mult_bins", {50, 0, 3e3}, "nMultBins, multMin, multMax"};

  Configurable<int> cfCent{"cfCent", 1, "centrality estimator"};
  Configurable<int> cfMult{"cfMult", 1, "multiplicity"};
  Configurable<bool> cfQA{"cfQA", true, "quality assurance"};

  Configurable<std::vector<float>> cfVertexZ{"cfVertexZ", {-10, 10.}, "vertex z position range: {min, max}[cm], with convention: min <= Vz < max"};
  Configurable<std::vector<float>> cfPt{"cfPt", {0.2, 5.0}, "transverse momentum range"};

  Configurable<std::string> cfFileWithWeights{"cfFileWithWeights", "~/O2/weights.root", "path to external ROOT file which holds all particle weights"};

  // *) Define and initialize all data members to be called in the main process* functions:
  // **) Task configuration:
  struct TaskConfiguration {
    bool fProcess[eProcess_N] = {false}; // Set what to process. See enum eProcess for full description. Set via implicit variables within a PROCESS_SWITCH clause.
    bool fDryRun = false;                // book all histos and run without filling and calculating anything
  } tc;                                  // you have to prepend "tc." for all objects name in this group later in the code

  // **) Particle histograms:
  struct ParticleHistograms {
    TList* fParticleHistogramsList = NULL; //!<! list to hold all control particle histograms
    TH1F* fHistPt[2] = {NULL};             // pt distribution of a particle [ 0 = rec, 1 = sim ]
    TH1F* fHistPhi[2] = {NULL};
    TH1F* histWeights = NULL;
  } pc; // you have to prepend "pc." for all objects name in this group later in the code

  struct EventHistograms {
    TList* fEventHistogramsList = NULL;
    TH1F* fHistCentr[2] = {NULL};
    TH1I* fHistMult[2] = {NULL};
    TH1F* fHistX[2] = {NULL};
    TH1F* fHistY[2] = {NULL};
    TH1F* fHistZ[2] = {NULL};
    TH2F* fQA = NULL;
    TH1F* fEventHistograms[eEventHistograms_N][2][2] = {{{NULL}}}; //! [ type - see enum eEventHistograms ][reco,sim][before, after event cuts]
  } event;

  TObject* GetObjectFromList(TList* list, const char* objectName)
  {
    // Get TObject pointer from TList, even if it's in some nested TList. Foreseen
    // to be used to fetch histograms or profiles from files directly.
    // Some ideas taken from TCollection::ls()
    // If you have added histograms directly to files (without TList's), then you can fetch them directly with
    // file->Get("hist-name").

    // Usage: TH1D *hist = (TH1D*) GetObjectFromList("some-valid-TList-pointer","some-object-name");

    // Insanity checks:
    if (!list) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    if (!objectName) {
      LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
    }
    if (0 == list->GetEntries()) {
      return NULL;
    }

    // The object is in the current base list:
    TObject* objectFinal = list->FindObject(objectName); // the final object I am after
    if (objectFinal) {
      return objectFinal;
    }

    // Otherwise, search for the object recursively in the nested lists:
    TObject* objectIter; // iterator object in the loop below
    TIter next(list);
    while ((objectIter = next())) // double round braces are to silence the warnings
    {
      if (TString(objectIter->ClassName()).EqualTo("TList")) {
        objectFinal = GetObjectFromList(reinterpret_cast<TList*>(objectIter), objectName);
        if (objectFinal)
          return objectFinal;
      }
    } // while(objectIter = next())

    return NULL;

  } // TObject* GetObjectFromList(TList *list, char *objectName)

  TH1F* GetHistogramWithWeights(const char* filePath, const char* runNumber)
  {
    // *) Return value:
    TH1F* hist = NULL;
    TList* baseList = NULL;     // base top-level list in the TFile, e.g. named "ccdb_object"
    TList* listWithRuns = NULL; // nested list with run-wise TList's holding run-specific weights

    // *) Determine from filePath if the file is on a local machine, or in home dir AliEn, or in CCDB:
    //    Algorithm: If filePath begins with "/alice/cern.ch/" then it's in the home dir AliEn;
    //                     If filePath begins with "/alice-ccdb.cern.ch/" then it's in CCDB. Therefore, files in AliEn and CCDB must be specified with abs path;
    //                     for local files both abs and relative paths are just fine.
    bool bFileIsInAliEn = false;
    bool bFileIsInCCDB = false;

    string pathstr = filePath;
    const string pathalien = "/alice/cern.ch/";
    const string pathccdb = "/alice-ccdb.cern.ch/";
    if (pathstr.find(pathalien) == 0) {
      bFileIsInAliEn = true;
    } else if (pathstr.find(pathccdb) == 0) {
      bFileIsInCCDB = true;
    }

    if (bFileIsInAliEn) {
      TGrid* alien = TGrid::Connect("alien", gSystem->Getenv("USER"), "", ""); // do not forget to add #include <TGrid.h> to the preamble of your analysis task
      if (!alien) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }
      TFile* weightsFile = TFile::Open(Form("alien://%s", filePath), "READ"); // yes, ROOT can open a file transparently, even if it's sitting in AliEn, with this specific syntax
      if (!weightsFile) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }
      weightsFile->GetObject("ccdb_object", baseList);
      if (!baseList) {
        // weightsFile->ls();
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      // Finally, from the top-level TList, get the desired nested TList => the technical problem here is that it can be nested at any level,
      // for thare there is a helper utility function GetObjectFromList(...) , see its implementation further below
      listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
        }
      }

    } else if (bFileIsInCCDB) {
      // File you want to access is in your home dir in CCDB:
      // Remember that here I do not access the file; instead, I directly access the object in that file.
      ccdb->setURL("http://alice-ccdb.cern.ch"); // to be able to use "ccdb" this object in your analysis task, see 4b/ below
      baseList = reinterpret_cast<TList*>(ccdb->get<TList>(TString(filePath).ReplaceAll("/alice-ccdb.cern.ch/", "").Data()));
      if (!baseList) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
        }
      }

      // OK, we got the desired TList with efficiency corrections, after that we can use the common code for all 3 cases (local, AliEn, CCDB, that common code is below)

    } else {
      // this is the local case:
      // Check if the external ROOT file exists at the specified path:
      if (gSystem->AccessPathName(filePath, kFileExists)) {
        LOGF(info, "\033[1;33m if(gSystem->AccessPathName(filePath,kFileExists)), filePath = %s \033[0m", filePath);
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      TFile* weightsFile = TFile::Open(filePath, "READ");
      if (!weightsFile) {
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      weightsFile->GetObject("ccdb_object", baseList);

      if (!baseList) {
        // weightsFile->ls();
        LOGF(fatal, "\033[1;31m%s at line %d\033[0m", __FUNCTION__, __LINE__);
      }

      listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumber));
      if (!listWithRuns) {
        TString runNumberWithLeadingZeroes = "000";
        runNumberWithLeadingZeroes += runNumber; // another try, with "000" prepended to run number
        listWithRuns = reinterpret_cast<TList*>(GetObjectFromList(baseList, runNumberWithLeadingZeroes.Data()));
        if (!listWithRuns) {
          // baseList->ls();
          LOGF(fatal, "\033[1;31m%s at line %d : this crash can happen if in the output file there is no list with weights for the current run number = %s\033[0m", __FUNCTION__, __LINE__ /*, tc.fRunNumber.Data()*/);
        }
      }

    } // end of  else

    // Here comes the common code for all three cases, where from "listWithRuns" you fetch the desired histogram with efficiency corrections:
    listWithRuns->ls();

    hist = reinterpret_cast<TH1F*>(listWithRuns->FindObject("histWithEfficiencyCorrections"));
    if (!hist) {
      LOGF(fatal, "no histweight");
    }

    // Once you have a valid pointer to "hist", add these technical lines:
    hist->SetDirectory(0);                                    // remove ownerhip from an external file
    TH1F* histClone = reinterpret_cast<TH1F*>(hist->Clone()); // yes, I have to clone here

    return histClone;

  } // end of TH1F* GetHistogramWithWeights(const char* filePath, const char* runNumber) {

  // ...

  // *) Define all member functions to be called in the main process* functions:
  template <eRecSim rs, typename T1>
  bool EventCuts(T1 const& collision)
  {
    vector<float> vertexZ = cfVertexZ.value;
    float vertexZmin = static_cast<float>(vertexZ[0]);
    float vertexZmax = static_cast<float>(vertexZ[1]);
    float posZ = collision.posZ();
    if (posZ < vertexZmin || posZ > vertexZmax)
      return false;
    return true;
  }

  template <eRecSim rs, typename T>
  bool ParticleCuts(T const& track)
  {
    vector<float> Pt = cfPt.value;
    float ptcutmin = static_cast<float>(Pt[0]);
    float ptcutmax = static_cast<float>(Pt[1]);
    float pt = track.pt();
    if (pt < ptcutmin || pt > ptcutmax)
      return false;
    return true;
  }

  template <eRecSim rs, typename T1, typename T2>
  void Steer(T1 const& collision, T2 const& tracks)
  {
    // Dry run:
    if (tc.fDryRun) {
      return;
    }
    // Print current run number:
    LOGF(info, "Run number: %d", collision.bc().runNumber());
    // Print centrality estimated with "FT0M" estimator:
    LOGF(info, "Centrality: %f", collision.centFT0M());
    // Print vertex X position:
    LOGF(info, "Vertex X position: %f", collision.posX());

    float zrec = 0., zsim = 0.;
    if constexpr (rs == eRec || rs == eRecAndSim) {
      event.fHistX[eRec]->Fill(collision.posX());
      event.fHistY[eRec]->Fill(collision.posY());
      event.fHistZ[eRec]->Fill(collision.posZ());
      event.fEventHistograms[eVertexZ][eRec][0]->Fill(collision.posZ());
      zrec = collision.posZ();
      float centr = 0;
      if (cfCent == 1)
        centr = collision.centFT0C();
      if (cfCent == 2)
        centr = collision.centFT0M();
      if (cfCent == 3)
        centr = collision.centFT0A();
      event.fHistCentr[eRec]->Fill(centr);
      if (cfMult == 1)
        event.fHistMult[eRec]->Fill(collision.multTPC());
      if (cfMult == 2)
        event.fHistMult[eRec]->Fill(collision.multFV0M());
      if (cfMult == 3)
        event.fHistMult[eRec]->Fill(collision.multFT0C());
      if (cfMult == 4)
        event.fHistMult[eRec]->Fill(collision.multFT0M());
      if (cfMult == 5)
        event.fHistMult[eRec]->Fill(collision.multNTracksPV());

      if constexpr (rs == eRecAndSim) {
        auto mccollision = collision.mcCollision();
        float b = mccollision.impactParameter();
        float centrsim = o2::constants::math::PI * b * b / 7.71;
        event.fHistX[eSim]->Fill(mccollision.posX());
        event.fHistY[eSim]->Fill(mccollision.posY());
        event.fHistZ[eSim]->Fill(mccollision.posZ());
        event.fEventHistograms[eVertexZ][eSim][0]->Fill(mccollision.posZ());
        zsim = mccollision.posZ();
        event.fHistCentr[eSim]->Fill(centrsim);
        event.fQA->Fill(centrsim, centr);
        // event.fHistMult[eSim]->Fill(tracks.size());
      }

      // *) Event cuts:
      if (!EventCuts<rs>(collision)) { // Main call for event cuts
        return;
      }
      event.fEventHistograms[eVertexZ][eRec][1]->Fill(zrec);
      if constexpr (rs == eRecAndSim)
        event.fEventHistograms[eVertexZ][eSim][1]->Fill(zsim);
    }

    // Main loop over particles:
    auto track = tracks.iteratorAt(0); // set the type and scope from one instance
    for (int64_t i = 0; i < tracks.size(); i++) {

      // Print track azimuthal angle:

      // Fill reconstructed ...:
      float ptrec = 0., ptsim = 0.;
      if constexpr (rs == eRec || rs == eRecAndSim) {
        // Fill track pt distribution:
        pc.fHistPt[eRec]->Fill(track.pt());
        event.fEventHistograms[ePt][eRec][0]->Fill(track.pt());
        ptrec = track.pt();
        pc.fHistPhi[eRec]->Fill(track.phi());

        // ... and corresponding MC truth simulated:
        // See https://github.com/AliceO2Group/O2Physics/blob/master/Tutorials/src/mcHistograms.cxx
        // See https://aliceo2group.github.io/analysis-framework/docs/datamodel/ao2dTables.html#montecarlo
        if constexpr (rs == eRecAndSim) {
          if (!track.has_mcParticle()) {
            LOGF(warning, "  No MC particle for this track, skip...");
            return;
          }
          auto mcparticle = track.mcParticle(); // corresponding MC truth simulated particle
          pc.fHistPt[eSim]->Fill(mcparticle.pt());
          event.fEventHistograms[ePt][eSim][0]->Fill(mcparticle.pt());
          ptsim = mcparticle.pt();
          pc.fHistPhi[eSim]->Fill(mcparticle.phi());
        } // end of if constexpr (rs == eRecAndSim) {

        // *) Particle cuts:
        if (!ParticleCuts<rs>(track)) { // Main call for particle cuts.
          continue;                     // not return!!
        }
        event.fEventHistograms[ePt][eRec][1]->Fill(ptrec);
        if constexpr (rs == eRecAndSim)
          event.fEventHistograms[ePt][eSim][0]->Fill(ptsim);

      } // if constexpr (rs == eRec || rs == eRecAndSim) {
    } // end of for (int64_t i = 0; i < tracks.size(); i++) {

  } // end of template <eRecSim rs, typename T1, typename T2> void Steer(T1 const& collision, T2 const& tracks) {

  // *) Initialize and book all objects:
  void init(InitContext&)
  {

    // ... code to book and initialize all analysis objects ...

    // *) Set automatically what to process, from an implicit variable "doprocessSomeProcessName" within a PROCESS_SWITCH clause:
    tc.fProcess[eProcessRec] = doprocessRec;
    tc.fProcess[eProcessRecSim] = doprocessRecSim;
    tc.fProcess[eProcessSim] = doprocessSim;

    // *) Configure your task using configurables in the json file:
    tc.fDryRun = cfDryRun;

    // *) Book base list:
    TList* temp = new TList();
    temp->SetOwner(true);
    fBaseList.setObject(temp);

    // *) Book and nest all other TLists:
    pc.fParticleHistogramsList = new TList();
    pc.fParticleHistogramsList->SetName("ParticleHistograms");
    pc.fParticleHistogramsList->SetOwner(true);
    fBaseList->Add(pc.fParticleHistogramsList); // any nested TList in the base TList appears as a subdir in the output ROOT file
    event.fEventHistogramsList = new TList();
    event.fEventHistogramsList->SetName("EventHistograms");
    event.fEventHistogramsList->SetOwner(true);
    fBaseList->Add(event.fEventHistogramsList);

    // *) Book pt distribution with binning defined through configurables in the json file:
    vector<float> l_pt_bins = cf_pt_bins.value; // define local array and initialize it from an array set in the configurables
    vector<float> l_phi_bins = cf_phi_bins.value;
    vector<float> l_centr_bins = cf_centr_bins.value;
    vector<float> l_x_bins = cf_x_bins.value;
    vector<float> l_y_bins = cf_y_bins.value;
    vector<float> l_z_bins = cf_z_bins.value;
    vector<float> l_mult_bins = cf_mult_bins.value;
    int nBins = static_cast<int>(l_pt_bins[0]);
    int nBinsphi = static_cast<int>(l_phi_bins[0]);
    int nBinscentr = static_cast<int>(l_centr_bins[0]);
    int nBinsx = static_cast<int>(l_x_bins[0]);
    int nBinsy = static_cast<int>(l_y_bins[0]);
    int nBinsz = static_cast<int>(l_z_bins[0]);
    int nBinsmult = static_cast<int>(l_mult_bins[0]);

    float min = l_pt_bins[1];
    float max = l_pt_bins[2];
    float minphi = l_phi_bins[1];
    float maxphi = l_phi_bins[2];
    float mincentr = l_centr_bins[1];
    float maxcentr = l_centr_bins[2];
    float minx = l_x_bins[1];
    float maxx = l_x_bins[2];
    float miny = l_y_bins[1];
    float maxy = l_y_bins[2];
    float minz = l_z_bins[1];
    float maxz = l_z_bins[2];
    float minmult = l_mult_bins[1];
    float maxmult = l_mult_bins[2];

    pc.fHistPt[eRec] = new TH1F("fHistPt[eRec]", "pt distribution for reconstructed particles", nBins, min, max);
    pc.fHistPt[eRec]->GetXaxis()->SetTitle("p_{T}");
    pc.fParticleHistogramsList->Add(pc.fHistPt[eRec]);
    pc.fHistPhi[eRec] = new TH1F("fHistPhi[eRec]", "phi distribution for reconstructed particles", nBinsphi, minphi, maxphi);
    pc.fHistPhi[eRec]->GetXaxis()->SetTitle("phi");
    pc.fParticleHistogramsList->Add(pc.fHistPhi[eRec]);

    pc.fHistPt[eSim] = new TH1F("fHistPt[eSim]", "pt distribution for simulated particles", nBins, min, max);
    pc.fHistPt[eSim]->GetXaxis()->SetTitle("p_{T}");
    pc.fParticleHistogramsList->Add(pc.fHistPt[eSim]);
    pc.fHistPhi[eSim] = new TH1F("fHistPhi[eSim]", "phi distribution for simulated particles", nBinsphi, minphi, maxphi);
    pc.fHistPhi[eSim]->GetXaxis()->SetTitle("phi");
    pc.fParticleHistogramsList->Add(pc.fHistPhi[eSim]);

    pc.histWeights = GetHistogramWithWeights(cfFileWithWeights.value.c_str(), "000123456");
    pc.fParticleHistogramsList->Add(pc.histWeights);

    event.fHistCentr[eRec] = new TH1F("fHistCentr[eRec]", "centrality distribution for reconstructed particles", nBinscentr, mincentr, maxcentr);
    event.fHistX[eRec] = new TH1F("fHistX[eRec]", "posX distribution for reconstructed particles", nBinsx, minx, maxx);
    event.fHistY[eRec] = new TH1F("fHistY[eRec]", "posY distribution for reconstructed particles", nBinsy, miny, maxy);
    event.fHistZ[eRec] = new TH1F("fHistZ[eRec]", "posZ distribution for reconstructed particles", nBinsz, minz, maxz);
    event.fHistMult[eRec] = new TH1I("fHistMult[eRec]", "mult distribution for reconstructed particles", nBinsmult, minmult, maxmult);
    event.fHistCentr[eRec]->GetXaxis()->SetTitle("centrality");
    event.fHistX[eRec]->GetXaxis()->SetTitle("x");
    event.fHistY[eRec]->GetXaxis()->SetTitle("y");
    event.fHistZ[eRec]->GetXaxis()->SetTitle("z");
    event.fHistMult[eRec]->GetXaxis()->SetTitle("multiplicity");
    event.fEventHistogramsList->Add(event.fHistCentr[eRec]);
    event.fEventHistogramsList->Add(event.fHistX[eRec]);
    event.fEventHistogramsList->Add(event.fHistY[eRec]);
    event.fEventHistogramsList->Add(event.fHistZ[eRec]);
    event.fEventHistogramsList->Add(event.fHistMult[eRec]);

    event.fHistCentr[eSim] = new TH1F("fHistCentr[eSim]", "centrality distribution for simulated particles", nBinscentr, mincentr, maxcentr);
    event.fHistX[eSim] = new TH1F("fHistX[eSim]", "posX distribution for simulated particles", nBinsx, minx, maxx);
    event.fHistY[eSim] = new TH1F("fHistY[eSim]", "posY distribution for simulated particles", nBinsy, miny, maxy);
    event.fHistZ[eSim] = new TH1F("fHistZ[eSim]", "posZ distribution for simulated particles", nBinsz, minz, maxz);
    event.fHistMult[eSim] = new TH1I("fHistMult[eSim]", "mult distribution for simulated particles", nBinsmult, minmult, maxmult);
    event.fHistCentr[eSim]->GetXaxis()->SetTitle("centrality");
    event.fHistX[eSim]->GetXaxis()->SetTitle("x");
    event.fHistY[eSim]->GetXaxis()->SetTitle("y");
    event.fHistZ[eSim]->GetXaxis()->SetTitle("z");
    event.fHistMult[eSim]->GetXaxis()->SetTitle("multiplicity");
    event.fEventHistogramsList->Add(event.fHistCentr[eSim]);
    event.fEventHistogramsList->Add(event.fHistX[eSim]);
    event.fEventHistogramsList->Add(event.fHistY[eSim]);
    event.fEventHistogramsList->Add(event.fHistZ[eSim]);
    event.fEventHistogramsList->Add(event.fHistMult[eSim]);

    const char* cevent[] = {"vertexZ", "Pt"};
    const char* cpro[] = {"rec", "sim"};
    const char* ccut[] = {"before", "after"};
    for (int i = 0; i < eEventHistograms_N; i++) {
      for (int j = 0; j < 2; j++) {
        for (int k = 0; k < 2; k++) {
          TString histname = Form("fEventHistograms[%s][%s][%s]", cevent[i], cpro[j], ccut[k]);
          TString histtitle = Form("%s distribution for %s, %s cut", cevent[i], cpro[j], ccut[k]);
          if (i == 0)
            event.fEventHistograms[i][j][k] = new TH1F(histname, histtitle, nBinsz, minz, maxz);
          if (i == 1)
            event.fEventHistograms[i][j][k] = new TH1F(histname, histtitle, nBins, min, max);
          event.fEventHistograms[i][j][k]->GetXaxis()->SetTitle(Form("%s", cevent[i]));
          event.fEventHistogramsList->Add(event.fEventHistograms[i][j][k]);
        }
      }
    }

    event.fQA = new TH2F("QA", "quality assurance", nBinscentr, mincentr, maxcentr, nBinscentr, mincentr, maxcentr);
    if (cfQA) {
      event.fEventHistogramsList->Add(event.fQA);
    }

  } // end of void init(InitContext&) {

  // A) Process only reconstructed data:
  void processRec(CollisionRec const& collision, aod::BCs const&, TracksRec const& tracks)
  {
    // ...

    // *) Steer all analysis steps:
    Steer<eRec>(collision, tracks);
  }
  PROCESS_SWITCH(MultiharmonicCorrelations, processRec, "process only reconstructed data", true); // yes, keep always one process switch "true", so that there is default running version

  // -------------------------------------------

  // B) Process both reconstructed and corresponding MC truth simulated data:
  void processRecSim(CollisionRecSim const& collision, aod::BCs const&, TracksRecSim const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {
    Steer<eRecAndSim>(collision, tracks);
  }
  PROCESS_SWITCH(MultiharmonicCorrelations, processRecSim, "process both reconstructed and corresponding MC truth simulated data", false);

  // -------------------------------------------

  // C) Process only simulated data:
  void processSim(CollisionSim const& /*collision*/, aod::BCs const&, TracksSim const& /*tracks*/)
  {
    // Steer<eSim>(collision, tracks); // TBI 20241105 not ready yet, but I do not really need this one urgently, since RecSim is working, and I need that one for efficiencies...
  }
  PROCESS_SWITCH(MultiharmonicCorrelations, processSim, "process only simulated data", false);

}; // struct MultiharmonicCorrelations {

// *) The final touch:
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<MultiharmonicCorrelations>(cfgc),
  };
} // WorkflowSpec...
