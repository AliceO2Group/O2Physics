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
#include <TComplex.h>
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

enum eCut {
  eBefore = 0,
  eAfter,
  eCut_N
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
  Configurable<std::vector<float>> cfPtBins{"cfPtBins", {1000, 0., 100.}, "nPtBins, ptMin, ptMax"};                   // example for an array
  Configurable<std::vector<float>> cfPhiBins{"cfPhiBins", {100, 0., o2::constants::math::TwoPI}, "nPhiBins, phiMin, phiMax"};
  Configurable<std::vector<float>> cfCentrBins{"cfCentrBins", {100, 0., 100.}, "nCentrBins, centrMin, centrMax"};
  Configurable<std::vector<float>> cfXBins{"cfXBins", {1000, -100., 100.}, "nXBins, xMin, xMax"};
  Configurable<std::vector<float>> cfYBins{"cfYBins", {1000, -100., 100.}, "nYBins, yMin, yMax"};
  Configurable<std::vector<float>> cfZBins{"cfZBins", {1000, -100., 100.}, "nZBins, zMin, zMax"};
  Configurable<std::vector<float>> cfMultBins{"cfMultBins", {50, 0, 3e3}, "nMultBins, multMin, multMax"};
  Configurable<std::vector<float>> cfTPCnclsBins{"cfTPCnclsBins", {100, 0., 1000.}, "ntpcnclsBins, tpnclsMin, tpcnclsMax"};
  Configurable<std::vector<float>> cfDCAxyBins{"cfDCAxyBins", {1000, -20., 20.}, "ndcaxyBins, dcaxyMin, dcaxyMax"};
  Configurable<std::vector<float>> cfDCAzBins{"cfDCAzBins", {1000, -10., 10.}, "ndcazBins, dcazMin, dcazMax"};
  Configurable<std::vector<float>> cfNcontrBins{"cfNcontrBins", {100, 0., 1000}, "nNContrBins, NContrMin, NContrMax"};

  Configurable<std::string> cfCent{"cfCent", "FT0C", "centrality estimator"};
  Configurable<std::string> cfMult{"cfMult", "TPC", "multiplicity"};
  Configurable<bool> cfQA{"cfQA", true, "quality assurance"};

  Configurable<std::vector<float>> cfVertexZ{"cfVertexZ", {-10, 10.}, "vertex z position range: {min, max}[cm], with convention: min <= Vz < max"};
  Configurable<std::vector<float>> cfPt{"cfPt", {0.2, 5.0}, "transverse momentum range"};
  Configurable<std::vector<float>> cfEta{"cfEta", {-0.8, 0.8}, "eta range"};

  Configurable<std::string> cfFileWithWeights{"cfFileWithWeights", "/alice-ccdb.cern.ch/Users/p/pengchon/test04", "path to external ROOT file which holds all particle weights"};

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
    TH1F* fHistCharge[2] = {NULL};  // charge distribution
    TH1F* fHistTPCncls[2] = {NULL}; // TPCNClsFindable
    TH1F* fHistTracksdcaXY[2] = {NULL};
    TH1F* fHistTracksdcaZ[2] = {NULL};
    TH1F* histWeights = NULL;
  } pc; // you have to prepend "pc." for all objects name in this group later in the code

  struct EventHistograms {
    TList* fEventHistogramsList = NULL;
    TH1F* fHistCentr[2] = {NULL};
    TH1I* fHistMult[2] = {NULL};
    TH1F* fHistX[2] = {NULL};
    TH1F* fHistY[2] = {NULL};
    TH1F* fHistZ[2] = {NULL};
    TH1I* fHistNContr = NULL;
    TH1F* fEventHistograms[eEventHistograms_N][2][2] = {{{NULL}}}; //! [ type - see enum eEventHistograms ][reco,sim][before, after event cuts]
  } event;

  struct QA {
    TList* fQAList = NULL;
    TH2F* fQA = NULL;
  } qa;

  static constexpr int maxHarmonic = 7;
  struct CorrelationVariables {
    TList* fCorrelationVariablesList = NULL;
    TProfile* pv22_centr = NULL;
    TProfile* pv32_centr = NULL;
    TProfile* pv42_centr = NULL;
    TProfile* pfour32_centr = NULL;
    TProfile* pfour42_centr = NULL;
    TComplex Qvector[maxHarmonic];
  } cor;

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
    LOGF(info, "bFileIsInCCDB= %d", bFileIsInCCDB);

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
    if (posZ < vertexZmin || posZ > vertexZmax || (!collision.sel8()))
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
    vector<float> Eta = cfEta.value;
    float etacutmin = static_cast<float>(Eta[0]);
    float etacutmax = static_cast<float>(Eta[1]);
    float eta = track.eta();
    if (pt < ptcutmin || pt > ptcutmax || eta < etacutmin || eta > etacutmax)
      return false;
    return true;
  }

  TComplex Q(Int_t n)
  {
    // Using the fact that Q{-n,p} = Q{n,p}^*.
    if (n >= 0) {
      return cor.Qvector[n];
    }
    return TComplex::Conjugate(cor.Qvector[-n]);
  }

  TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4)
  { // Generic four-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4)]>.
    TComplex four = Q(n1) * Q(n2) * Q(n3) * Q(n4) - Q(n1 + n2) * Q(n3) * Q(n4) - Q(n2) * Q(n1 + n3) * Q(n4) - Q(n1) * Q(n2 + n3) * Q(n4) + 2. * Q(n1 + n2 + n3) * Q(n4) - Q(n2) * Q(n3) * Q(n1 + n4) + Q(n2 + n3) * Q(n1 + n4) - Q(n1) * Q(n3) * Q(n2 + n4) + Q(n1 + n3) * Q(n2 + n4) + 2. * Q(n3) * Q(n1 + n2 + n4) - Q(n1) * Q(n2) * Q(n3 + n4) + Q(n1 + n2) * Q(n3 + n4) + 2. * Q(n2) * Q(n1 + n3 + n4) + 2. * Q(n1) * Q(n2 + n3 + n4) - 6. * Q(n1 + n2 + n3 + n4);

    return four;

  } // TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4)

  template <eRecSim rs, typename T1, typename T2>
  void Steer(T1 const& collision, T2 const& tracks)
  {
    // Dry run:
    if (tc.fDryRun) {
      return;
    }
    // Print current run number:
    // LOGF(info, "Run number: %d", collision.bc().runNumber());

    float zrec = 0., zsim = 0., centr = 0, M = 0.;

    if constexpr (rs == eRec || rs == eRecAndSim) {
      event.fHistX[eRec]->Fill(collision.posX());
      event.fHistY[eRec]->Fill(collision.posY());
      event.fHistZ[eRec]->Fill(collision.posZ());
      event.fEventHistograms[eVertexZ][eRec][0]->Fill(collision.posZ());
      zrec = collision.posZ();
      event.fHistNContr->Fill(collision.numContrib());
      if (cfCent.value == "FT0C")
        centr = collision.centFT0C();
      else if (cfCent.value == "FT0M")
        centr = collision.centFT0M();
      else if (cfCent.value == "FT0A")
        centr = collision.centFT0A();
      event.fHistCentr[eRec]->Fill(centr);

      std::string multType = "TPC";
      if (cfMult.value == "TPC")
        M = collision.multTPC();
      else if (cfMult.value == "FV0M")
        M = collision.multFV0M();
      else if (cfMult.value == "FT0C")
        M = collision.multFT0C();
      else if (cfMult.value == "FT0M")
        M = collision.multFT0M();
      else if (cfMult.value == "NTracksPV")
        M = collision.multNTracksPV();
      event.fHistMult[eRec]->Fill(M);

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
        qa.fQA->Fill(centrsim, centr);
        centr = centrsim;
      }

      // *) Event cuts:
      float centrcut = 80.;
      if (!EventCuts<rs>(collision) || centr > centrcut) { // Main call for event cuts
        return;
      }
      event.fEventHistograms[eVertexZ][eRec][1]->Fill(zrec);
      if constexpr (rs == eRecAndSim)
        event.fEventHistograms[eVertexZ][eSim][1]->Fill(zsim);
    }

    // before loop over particles
    float phi = 0;
    for (int ih = 0; ih < maxHarmonic; ih++) {
      cor.Qvector[ih] = TComplex(0., 0.);
    }

    // Main loop over particles:
    for (const auto& track : tracks) {

      // Fill reconstructed ...:
      float ptrec = 0., ptsim = 0.;
      if constexpr (rs == eRec || rs == eRecAndSim) {
        // Fill track pt distribution:
        pc.fHistPt[eRec]->Fill(track.pt());
        event.fEventHistograms[ePt][eRec][0]->Fill(track.pt());
        ptrec = track.pt();
        phi = track.phi();
        pc.fHistPhi[eRec]->Fill(track.phi());
        pc.fHistCharge[eRec]->Fill(track.sign());
        pc.fHistTPCncls[eRec]->Fill(track.tpcNClsFindable());
        pc.fHistTracksdcaXY[eRec]->Fill(track.dcaXY());
        pc.fHistTracksdcaZ[eRec]->Fill(track.dcaZ());

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
          phi = mcparticle.phi();
          pc.fHistPhi[eSim]->Fill(mcparticle.phi());
          // pc.fHistCharge[eSim]->Fill(mcparticle.sign());
          // pc.fHistTPCncls[eSim]->Fill(mcparticle.tpcNClsFindable());
          // pc.fHistTracksdcaXY[eSim]->Fill(mcparticle.dcaXY());
          // pc.fHistTracksdcaZ[eSim]->Fill(mcparticle.dcaZ());
        } // end of if constexpr (rs == eRecAndSim)

        // *) Particle cuts:
        if (!ParticleCuts<rs>(track)) { // Main call for particle cuts.
          continue;                     // not return!!
        }
        event.fEventHistograms[ePt][eRec][1]->Fill(ptrec);
        if constexpr (rs == eRecAndSim)
          event.fEventHistograms[ePt][eSim][0]->Fill(ptsim);

      } // if constexpr (rs == eRec || rs == eRecAndSim)

      // analysis in the loop over particle
      for (int ih = 0; ih < maxHarmonic; ih++) {
        cor.Qvector[ih] += TComplex(TMath::Cos(ih * phi), TMath::Sin(ih * phi));
      }
    } // end of for (auto track: tracks)
    // calculate correlations
    float four32 = Four(3, 2, -3, -2).Re() / Four(0, 0, 0, 0).Re();
    float four42 = Four(4, 2, -4, -2).Re() / Four(0, 0, 0, 0).Re();
    float v22 = (Q(2).Rho2() - M) / (M * (M - 1.));
    float v32 = (Q(3).Rho2() - M) / (M * (M - 1.));
    float v42 = (Q(4).Rho2() - M) / (M * (M - 1.));

    cor.pv22_centr->Fill(centr, v22, M * (M - 1));
    cor.pv32_centr->Fill(centr, v32, M * (M - 1));
    cor.pv42_centr->Fill(centr, v42, M * (M - 1));
    cor.pfour32_centr->Fill(centr, four32, M * (M - 1));
    cor.pfour42_centr->Fill(centr, four42, M * (M - 1));

  } // end of template <eRecSim rs, typename T1, typename T2> void Steer(T1 const& collision, T2 const& tracks)

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

    qa.fQAList = new TList();
    qa.fQAList->SetName("QAHistograms");
    qa.fQAList->SetOwner(true);
    fBaseList->Add(qa.fQAList);

    cor.fCorrelationVariablesList = new TList();
    cor.fCorrelationVariablesList->SetName("CorrelationVariables");
    cor.fCorrelationVariablesList->SetOwner(true);
    fBaseList->Add(cor.fCorrelationVariablesList);

    // *) Book pt distribution with binning defined through configurables in the json file:
    vector<float> l_pt_bins = cfPtBins.value; // define local array and initialize it from an array set in the configurables
    vector<float> l_phi_bins = cfPhiBins.value;
    vector<float> l_centr_bins = cfCentrBins.value;
    vector<float> l_x_bins = cfXBins.value;
    vector<float> l_y_bins = cfYBins.value;
    vector<float> l_z_bins = cfZBins.value;
    vector<float> l_mult_bins = cfMultBins.value;
    vector<float> l_tpcncls_bins = cfTPCnclsBins.value;
    vector<float> l_dcaxy_bins = cfDCAxyBins.value;
    vector<float> l_dcaz_bins = cfDCAzBins.value;
    vector<float> l_ncontr_bins = cfNcontrBins.value;
    int nBins = static_cast<int>(l_pt_bins[0]);
    int nBinsphi = static_cast<int>(l_phi_bins[0]);
    int nBinscentr = static_cast<int>(l_centr_bins[0]);
    int nBinsx = static_cast<int>(l_x_bins[0]);
    int nBinsy = static_cast<int>(l_y_bins[0]);
    int nBinsz = static_cast<int>(l_z_bins[0]);
    int nBinsmult = static_cast<int>(l_mult_bins[0]);
    int nBinscharge = 2;
    int nBinstpcncls = static_cast<int>(l_tpcncls_bins[0]);
    int nBinsdcaxy = static_cast<int>(l_dcaxy_bins[0]);
    int nBinsdcaz = static_cast<int>(l_dcaz_bins[0]);
    int nBinsncontr = static_cast<int>(l_ncontr_bins[0]);

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
    float mincharge = -2.;
    float maxcharge = 2.;
    float mintpcncls = l_tpcncls_bins[1];
    float maxtpcncls = l_tpcncls_bins[2];
    float mindcaxy = l_dcaxy_bins[1];
    float maxdcaxy = l_dcaxy_bins[2];
    float mindcaz = l_dcaz_bins[1];
    float maxdcaz = l_dcaz_bins[2];
    float maxncontr = l_ncontr_bins[1];
    float minncontr = l_ncontr_bins[2];

    pc.fHistPt[eRec] = new TH1F("fHistPt[eRec]", "pt distribution for reconstructed particles", nBins, min, max);
    pc.fHistPhi[eRec] = new TH1F("fHistPhi[eRec]", "phi distribution for reconstructed particles", nBinsphi, minphi, maxphi);
    pc.fHistCharge[eRec] = new TH1F("fHistCharge[eRec]", "charge distribution for reconstructed particles", nBinscharge, mincharge, maxcharge);
    pc.fHistTPCncls[eRec] = new TH1F("fHistTPCncls[eRec]", "tpcncls distribution for reconstructed particles", nBinstpcncls, mintpcncls, maxtpcncls);
    pc.fHistTracksdcaXY[eRec] = new TH1F("fHistTracksdcaXY[eRec]", "dcaxy distribution for reconstructed particles", nBinsdcaxy, mindcaxy, maxdcaxy);
    pc.fHistTracksdcaZ[eRec] = new TH1F("fHistTracksdcaZ[eRec]", "dcaz distribution for reconstructed particles", nBinsdcaz, mindcaz, maxdcaz);
    pc.fHistPt[eRec]->GetXaxis()->SetTitle("p_{T}");
    pc.fHistPhi[eRec]->GetXaxis()->SetTitle("phi");
    pc.fHistCharge[eRec]->GetXaxis()->SetTitle("charge");
    pc.fHistTPCncls[eRec]->GetXaxis()->SetTitle("TPCNClsFindable");
    pc.fHistTracksdcaXY[eRec]->GetXaxis()->SetTitle("DCA XY");
    pc.fHistTracksdcaZ[eRec]->GetXaxis()->SetTitle("DCA Z");
    pc.fParticleHistogramsList->Add(pc.fHistPt[eRec]);
    pc.fParticleHistogramsList->Add(pc.fHistPhi[eRec]);
    pc.fParticleHistogramsList->Add(pc.fHistCharge[eRec]);
    pc.fParticleHistogramsList->Add(pc.fHistTPCncls[eRec]);
    pc.fParticleHistogramsList->Add(pc.fHistTracksdcaXY[eRec]);
    pc.fParticleHistogramsList->Add(pc.fHistTracksdcaZ[eRec]);

    pc.fHistPt[eSim] = new TH1F("fHistPt[eSim]", "pt distribution for simulated particles", nBins, min, max);
    pc.fHistPhi[eSim] = new TH1F("fHistPhi[eSim]", "phi distribution for simulated particles", nBinsphi, minphi, maxphi);
    pc.fHistCharge[eSim] = new TH1F("fHistCharge[eSim]", "charge distribution for simulated particles", nBinscharge, mincharge, maxcharge);
    pc.fHistTPCncls[eSim] = new TH1F("fHistTPCncls[eSim]", "tpcncls distribution for simulated particles", nBinstpcncls, minphi, maxtpcncls);
    pc.fHistTracksdcaXY[eSim] = new TH1F("fHistTracksdcaXY[eSim]", "dcaxy distribution for simulated particles", nBinsdcaxy, mindcaxy, maxdcaxy);
    pc.fHistTracksdcaZ[eSim] = new TH1F("fHistTracksdcaZ[eSim]", "dcaz distribution for simulated particles", nBinsdcaz, mindcaz, maxdcaz);
    pc.fHistPt[eSim]->GetXaxis()->SetTitle("p_{T}");
    pc.fHistPhi[eSim]->GetXaxis()->SetTitle("phi");
    pc.fHistCharge[eSim]->GetXaxis()->SetTitle("charge");
    pc.fHistTPCncls[eSim]->GetXaxis()->SetTitle("TPCNClsFindable");
    pc.fHistTracksdcaXY[eSim]->GetXaxis()->SetTitle("DCA XY");
    pc.fHistTracksdcaZ[eSim]->GetXaxis()->SetTitle("DCA Z");
    pc.fParticleHistogramsList->Add(pc.fHistPt[eSim]);
    pc.fParticleHistogramsList->Add(pc.fHistPhi[eSim]);
    pc.fParticleHistogramsList->Add(pc.fHistCharge[eSim]);
    pc.fParticleHistogramsList->Add(pc.fHistTPCncls[eSim]);
    pc.fParticleHistogramsList->Add(pc.fHistTracksdcaXY[eSim]);
    pc.fParticleHistogramsList->Add(pc.fHistTracksdcaZ[eSim]);

    pc.histWeights = GetHistogramWithWeights(cfFileWithWeights.value.c_str(), "000123456");
    pc.fParticleHistogramsList->Add(pc.histWeights);

    event.fHistCentr[eRec] = new TH1F("fHistCentr[eRec]", "centrality distribution for reconstructed particles", nBinscentr, mincentr, maxcentr);
    event.fHistX[eRec] = new TH1F("fHistX[eRec]", "posX distribution for reconstructed particles", nBinsx, minx, maxx);
    event.fHistY[eRec] = new TH1F("fHistY[eRec]", "posY distribution for reconstructed particles", nBinsy, miny, maxy);
    event.fHistZ[eRec] = new TH1F("fHistZ[eRec]", "posZ distribution for reconstructed particles", nBinsz, minz, maxz);
    event.fHistMult[eRec] = new TH1I("fHistMult[eRec]", "mult distribution for reconstructed particles", nBinsmult, minmult, maxmult);
    event.fHistNContr = new TH1I("fHistNContr", "NContr distribution", nBinsncontr, minncontr, maxncontr);
    event.fHistCentr[eRec]->GetXaxis()->SetTitle("centrality");
    event.fHistX[eRec]->GetXaxis()->SetTitle("x");
    event.fHistY[eRec]->GetXaxis()->SetTitle("y");
    event.fHistZ[eRec]->GetXaxis()->SetTitle("z");
    event.fHistMult[eRec]->GetXaxis()->SetTitle("multiplicity");
    event.fHistNContr->GetXaxis()->SetTitle("numContrib");
    event.fEventHistogramsList->Add(event.fHistCentr[eRec]);
    event.fEventHistogramsList->Add(event.fHistX[eRec]);
    event.fEventHistogramsList->Add(event.fHistY[eRec]);
    event.fEventHistogramsList->Add(event.fHistZ[eRec]);
    event.fEventHistogramsList->Add(event.fHistMult[eRec]);
    event.fEventHistogramsList->Add(event.fHistNContr);

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
      for (int j = 0; j < eRecAndSim; j++) {
        for (int k = 0; k < eCut_N; k++) {
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

    qa.fQA = new TH2F("QA", "quality assurance", nBinscentr, mincentr, maxcentr, nBinscentr, mincentr, maxcentr);
    if (cfQA) {
      qa.fQAList->Add(qa.fQA);
    }

    // float quantiles[10] = {0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    float quantiles[10] = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80};
    cor.pv22_centr = new TProfile("pv22", "profile of v_{2}^{2}", 9, quantiles);
    cor.pv32_centr = new TProfile("pv32", "profile of v_{3}^{2}", 9, quantiles);
    cor.pv42_centr = new TProfile("pv42", "profile of v_{4}^{2}", 9, quantiles);
    cor.pfour32_centr = new TProfile("pfour32", "profile of v_{2}^{2}*v_{3}^{2}", 9, quantiles);
    cor.pfour42_centr = new TProfile("pfour42", "profile of v_{2}^{2}*v_{4}^{2}", 9, quantiles);
    cor.pv22_centr->GetYaxis()->SetTitle("v_{2}^{2}");
    cor.pv32_centr->GetYaxis()->SetTitle("v_{3}^{2}");
    cor.pv42_centr->GetYaxis()->SetTitle("v_{4}^{2}");
    cor.pv22_centr->GetXaxis()->SetTitle("centrality");
    cor.pv32_centr->GetXaxis()->SetTitle("centrality");
    cor.pv42_centr->GetXaxis()->SetTitle("centrality");
    cor.pfour32_centr->GetYaxis()->SetTitle("v_{2}^{2}v_{3}^{2}");
    cor.pfour42_centr->GetYaxis()->SetTitle("v_{2}^{2}v_{4}^{2}");
    cor.pfour32_centr->GetXaxis()->SetTitle("centrality");
    cor.pfour42_centr->GetXaxis()->SetTitle("centrality");
    cor.fCorrelationVariablesList->Add(cor.pv22_centr);
    cor.fCorrelationVariablesList->Add(cor.pv32_centr);
    cor.fCorrelationVariablesList->Add(cor.pv42_centr);
    cor.fCorrelationVariablesList->Add(cor.pfour32_centr);
    cor.fCorrelationVariablesList->Add(cor.pfour42_centr);

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
