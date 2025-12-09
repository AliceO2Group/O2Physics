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

/// \file alice3TrackingTranslator.cxx
///
/// \brief Translator task to convert tracking software to the AO2D format digestible with the O2Physics analysis framework
///
/// \author Nicolò Jacazio, Universita del Piemonte Orientale (IT)
///

#include "ALICE3/DataModel/collisionAlice3.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/TrackParametrizationWithError.h>

#include <TFile.h>
#include <TList.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTree.h>

#include <cmath>
#include <map>
#include <string>
#include <vector>

#ifdef __CLING__
#pragma link C++ class std::vector < std::vector < unsigned int>> + ;
#pragma link C++ class std::vector < std::vector < std::uint32_t>> + ;
#endif

TString inputPath;

struct Alice3TrackingTranslator {
  o2::framework::Produces<o2::aod::Collisions> tableCollisions;
  o2::framework::Produces<o2::aod::McCollisionLabels> tableMcCollisionLabels;
  o2::framework::Produces<o2::aod::StoredTracks> tableStoredTracks;
  o2::framework::Produces<o2::aod::TracksExtension> tableTracksExtension;
  o2::framework::Produces<o2::aod::StoredTracksCov> tableStoredTracksCov;
  o2::framework::Produces<o2::aod::TracksCovExtension> tableTracksCovExtension;
  o2::framework::Produces<o2::aod::McTrackLabels> tableMcTrackLabels;
  o2::framework::Produces<o2::aod::TracksDCA> tableTracksDCA;
  o2::framework::Produces<o2::aod::TracksDCACov> tableTracksDCACov;
  o2::framework::Produces<o2::aod::CollisionsAlice3> tableCollisionsAlice3;
  o2::framework::Produces<o2::aod::TracksAlice3> tableTracksAlice3;
  o2::framework::Produces<o2::aod::TracksExtraA3> tableTracksExtraA3;

  o2::framework::Produces<o2::aod::StoredTracksExtra_002> tableStoredTracksExtra;
  o2::framework::Produces<o2::aod::TrackSelection> tableTrackSelection;
  o2::framework::Produces<o2::aod::TrackSelectionExtension> tableTrackSelectionExtension;
  o2::framework::Produces<o2::aod::StoredMcParticles> tableStoredMcParticles;
  o2::framework::Produces<o2::aod::McCollisions> tableMcCollisions;

  void init(o2::framework::InitContext&)
  {
    // Initialization if needed
    LOG(info) << "Alice3TrackingTranslator init called";
    // Load dictionary for nested vector
    gInterpreter->GenerateDictionary("vector<vector<unsigned int>>", "vector");
  }

#define SETADDRESS(branchname, branchvar)                            \
  if (mTree->SetBranchAddress(branchname, &branchvar)) {             \
    LOG(fatal) << "Could not set branch address for " << branchname; \
  }
  struct FileStruct {
    FileStruct(std::string filename, std::string treename) : mFile(filename.c_str(), "READ")
    {
      if (mFile.IsZombie()) {
        LOG(fatal) << "Could not open file " << filename;
      }
      mFile.GetObject(treename.c_str(), mTree);
      if (mTree) {
        LOG(info) << "Found " << treename << " tree with " << mTree->GetEntries() << " entries.";
      } else {
        LOG(fatal) << treename << " tree not found in " << filename;
      }
    }
    void setEventEntry(Long64_t entry)
    {
      if (mTree->GetEntry(entry) < 0) {
        LOG(fatal) << "Could not read entry " << entry << " from tree.";
      }
    }
    Long64_t getEntries() const { return mTree->GetEntries(); }
    TFile mFile;
    TTree* mTree;
  };

  struct ParticleStruct : public FileStruct {
    ParticleStruct(std::string filename, std::string treename) : FileStruct(filename, treename)
    {
      // mTree->Print();
      SETADDRESS("particle_type", m_particle_type);
      SETADDRESS("vx", m_vx);
      SETADDRESS("vy", m_vy);
      SETADDRESS("vz", m_vz);
      SETADDRESS("vt", m_vt);
      SETADDRESS("px", m_px);
      SETADDRESS("py", m_py);
      SETADDRESS("pz", m_pz);
      SETADDRESS("m", m_m);
      SETADDRESS("p", m_p);
    }
    std::vector<int>* m_particle_type = nullptr;
    std::vector<float>* m_vx = nullptr;
    std::vector<float>* m_vy = nullptr;
    std::vector<float>* m_vz = nullptr;
    std::vector<float>* m_vt = nullptr;
    std::vector<float>* m_px = nullptr;
    std::vector<float>* m_py = nullptr;
    std::vector<float>* m_pz = nullptr;
    std::vector<float>* m_m = nullptr;
    std::vector<float>* m_p = nullptr;
  };

  struct TrackStruct : public FileStruct {
    TrackStruct(std::string filename, std::string treename) : FileStruct(filename, treename)
    {
      mTree->Print();
      // Set branch addresses for ACTS track parameters

      SETADDRESS("event_nr", m_event_nr);
      SETADDRESS("nMeasurements", m_nMeasurements);
      SETADDRESS("nStates", m_nStates);
      SETADDRESS("nHoles", m_nHoles);
      SETADDRESS("chi2Sum", m_chi2Sum);
      SETADDRESS("NDF", m_NDF);
      SETADDRESS("eLOC0_fit", m_eLOC0_fit);
      SETADDRESS("eLOC1_fit", m_eLOC1_fit);
      SETADDRESS("ePHI_fit", m_ePHI_fit);
      SETADDRESS("eTHETA_fit", m_eTHETA_fit);
      SETADDRESS("eQOP_fit", m_eQOP_fit);
      SETADDRESS("eT_fit", m_eT_fit);
      SETADDRESS("nMajorityHits", m_nMajorityHits);
      // SETADDRESS("majorityParticleId", m_majorityParticleId);
      mTree->SetBranchAddress("majorityParticleId", &m_majorityParticleId);
      SETADDRESS("t_charge", m_t_charge);
      SETADDRESS("t_vx", m_t_vx);
      SETADDRESS("t_vy", m_t_vy);
      SETADDRESS("t_vz", m_t_vz);
      SETADDRESS("t_time", m_t_time);
      SETADDRESS("t_px", m_t_px);
      SETADDRESS("t_py", m_t_py);
      SETADDRESS("t_pz", m_t_pz);
      SETADDRESS("t_theta", m_t_theta);
      SETADDRESS("t_phi", m_t_phi);
      SETADDRESS("t_pT", m_t_pT);
      SETADDRESS("t_eta", m_t_eta);
    }
    // Define track-related members here
    UInt_t* m_event_nr = nullptr;
    std::vector<uint32_t>* m_nMeasurements = nullptr;
    std::vector<uint32_t>* m_nStates = nullptr;
    std::vector<uint32_t>* m_nHoles = nullptr;
    std::vector<float>* m_chi2Sum = nullptr;
    std::vector<uint32_t>* m_NDF = nullptr;
    // Fitted track parameters
    std::vector<float>* m_eLOC0_fit = nullptr;  // local position 0 (typically y in local frame)
    std::vector<float>* m_eLOC1_fit = nullptr;  // local position 1 (typically z in local frame)
    std::vector<float>* m_ePHI_fit = nullptr;   // azimuthal angle
    std::vector<float>* m_eTHETA_fit = nullptr; // polar angle
    std::vector<float>* m_eQOP_fit = nullptr;   // q/m_p (charge over momentum)
    std::vector<float>* m_eT_fit = nullptr;     // time

    // The majority truth particle info
    std::vector<unsigned int>* m_nMajorityHits = nullptr;                    /// The number of hits from majority particle
    std::vector<std::vector<std::uint32_t>>* m_majorityParticleId = nullptr; /// The particle Id of the majority particle
    std::vector<int>* m_t_charge = nullptr;                                  /// Charge of majority particle
    std::vector<float>* m_t_time = nullptr;                                  /// Time of majority particle
    std::vector<float>* m_t_vx = nullptr;                                    /// Vertex x positions of majority particle
    std::vector<float>* m_t_vy = nullptr;                                    /// Vertex y positions of majority particle
    std::vector<float>* m_t_vz = nullptr;                                    /// Vertex z positions of majority particle
    std::vector<float>* m_t_px = nullptr;                                    /// Initial momenta m_px of majority particle
    std::vector<float>* m_t_py = nullptr;                                    /// Initial momenta m_py of majority particle
    std::vector<float>* m_t_pz = nullptr;                                    /// Initial momenta m_pz of majority particle
    std::vector<float>* m_t_theta = nullptr;                                 /// Initial momenta theta of majority particle
    std::vector<float>* m_t_phi = nullptr;                                   /// Initial momenta phi of majority particle
    std::vector<float>* m_t_pT = nullptr;                                    /// Initial momenta pT of majority particle
    std::vector<float>* m_t_eta = nullptr;                                   /// Initial momenta eta of majority particle
  };

  struct HitsStruct : public FileStruct {
    HitsStruct(std::string filename, std::string treename) : FileStruct(filename, treename)
    {
      mTree->Print();
      SETADDRESS("barcode", barcode);
    }
    std::vector<unsigned int>* barcode = nullptr;
  };

  void process(o2::aod::BCs const&)
  {
    LOG(info) << "Alice3TrackingTranslator process called";
    // Find all ROOT files in the folder
    std::vector<std::string> rootFiles;
    TSystemDirectory dir(inputPath.Data(), inputPath.Data());
    TList* filesList = dir.GetListOfFiles();
    if (filesList) {
      TIter next(filesList);
      TSystemFile* file;
      while ((file = static_cast<TSystemFile*>(next()))) {
        TString fname = file->GetName();
        if (!file->IsDirectory() && fname.EndsWith(".root")) {
          TString fullPath = TString::Format("%s/%s", inputPath.Data(), fname.Data());
          rootFiles.push_back(fullPath.Data());
        }
      }
      delete filesList;
    }
    // Open all found ROOT files
    std::map<std::string, std::string> files;
    for (const auto& filename : rootFiles) {
      LOG(info) << "Opened ROOT file: " << filename;
      // Extract just the filename without path
      TString tfilename(filename.c_str());
      TString justFilename = gSystem->BaseName(tfilename);
      LOG(info) << "Processing file: " << justFilename.Data();
      files[justFilename.Data()] = filename;
    }

    // Now open the files to translate and read the trees
    ParticleStruct fileParticles(files["particles_simulation.root"], "particles");
    // FileStruct fileVertices(files["performance_vertexing.root"], "vertexing");
    TrackStruct fileTracksummary(files["tracksummary_ckf.root"], "tracksummary");
    // HitsStruct fileHits(files["hits.root"], "hits");

    const Long64_t kEvents = fileParticles.getEntries();
    for (Long64_t iEvent = 0; iEvent < kEvents; ++iEvent) {
      fileParticles.setEventEntry(iEvent);
      // fileVertices.setEventEntry(iEvent);
      fileTracksummary.setEventEntry(iEvent);
      // fileHits.setEventEntry(iEvent);

      LOG(info) << "Processing event " << iEvent << "/" << kEvents;

      // Create collision entry for this event
      // TODO: Extract proper collision position from vertex file if available
      float collisionX = 0.0f;
      float collisionY = 0.0f;
      float collisionZ = 0.0f;

      tableCollisions(0,          // bcId
                      collisionX, // posX
                      collisionY, // posY
                      collisionZ, // posZ
                      0.0f,       // covXX
                      0.0f,       // covXY
                      0.0f,       // covXZ
                      0.0f,       // covYY
                      0.0f,       // covYZ
                      0.0f,       // covZZ
                      0,          // flags
                      0.0f,       // m_chi2Sum
                      0,          // numContrib
                      0.0f,       // collisionTime
                      0.0f);      // collisionTimeRes

      tableMcCollisionLabels(iEvent, // mcCollisionId
                             0);     // mcMask

      tableCollisionsAlice3(0.f); // multDensity

      // Fill MC particles
      int mothers[2] = {-1, -1};
      int daughters[2] = {-1, -1};
      const size_t nParticlesGen = fileParticles.m_vx->size();
      for (size_t iParticle = 0; iParticle < nParticlesGen; ++iParticle) {
        continue;
        if (iParticle == 0) {
          tableMcCollisions(0,                                 // mccollision::BCId,
                            0,                                 // mccollision::GeneratorsID,
                            fileParticles.m_vx->at(iParticle), // mccollision::PosX,
                            fileParticles.m_vy->at(iParticle), // mccollision::PosY,
                            fileParticles.m_vz->at(iParticle), // mccollision::PosZ
                            fileParticles.m_vt->at(iParticle), // mccollision::T
                            1.0f,                              // mccollision::Weight
                            0.0f,                              // mccollision::ImpactParameter,
                            0.f);                              // mccollision::EventPlaneAngle,
        }

        uint8_t flags = 0;
        flags |= o2::aod::mcparticle::enums::PhysicalPrimary;
        tableStoredMcParticles(tableMcCollisions.lastIndex(),                                                  // mcCollisionId
                               fileParticles.m_particle_type->at(iParticle),                                   // pdgCode
                               0,                                                                              // statusCode
                               flags,                                                                          // flags
                               mothers,                                                                        // mothersIds
                               daughters,                                                                      // daughtersIdSlice
                               1.0f,                                                                           // weight
                               fileParticles.m_px->at(iParticle),                                              // m_px
                               fileParticles.m_py->at(iParticle),                                              // m_py
                               fileParticles.m_pz->at(iParticle),                                              // m_pz
                               std::hypot(fileParticles.m_p->at(iParticle), fileParticles.m_m->at(iParticle)), // e
                               fileParticles.m_vx->at(iParticle),                                              // m_vx
                               fileParticles.m_vy->at(iParticle),                                              // m_vy
                               fileParticles.m_vz->at(iParticle),                                              // m_vz
                               fileParticles.m_vt->at(iParticle));                                             // m_vt
      }

      // Convert tracks from ACTS to ALICE format
      const size_t nParticles = fileTracksummary.m_t_vx->size();
      const size_t nTracks = fileTracksummary.m_eLOC0_fit->size();
      for (size_t iTrack = 0; iTrack < nTracks; ++iTrack) {
        LOG(info) << "Processing track " << iTrack << "/" << nTracks << " (nParticles=" << nParticles << ") nParticlesGen=" << nParticlesGen;
        const size_t iParticle = iTrack;
        if (iParticle == 0) {
          tableMcCollisions(0,                                        // mccollision::BCId,
                            0,                                        // mccollision::GeneratorsID,
                            fileTracksummary.m_t_vx->at(iParticle),   // mccollision::PosX,
                            fileTracksummary.m_t_vy->at(iParticle),   // mccollision::PosY,
                            fileTracksummary.m_t_vz->at(iParticle),   // mccollision::PosZ
                            fileTracksummary.m_t_time->at(iParticle), // mccollision::T
                            1.0f,                                     // mccollision::Weight
                            0.0f,                                     // mccollision::ImpactParameter,
                            0.f);                                     // mccollision::EventPlaneAngle,
        }
        uint8_t flags = 0;
        flags |= o2::aod::mcparticle::enums::PhysicalPrimary;

        // fileTracksummary.m_majorityParticleId->at(iParticle).at(2), // pdgCode
        const size_t iParticleGen = fileTracksummary.m_majorityParticleId->at(iParticle).empty() ? 0 : fileTracksummary.m_majorityParticleId->at(iParticle).at(0);
        tableStoredMcParticles(tableMcCollisions.lastIndex(),                   // mcCollisionId
                               fileParticles.m_particle_type->at(iParticleGen), // pdgCode
                               0,                                               // statusCode
                               flags,                                           // flags
                               mothers,                                         // mothersIds
                               daughters,                                       // daughtersIdSlice
                               1.0f,                                            // weight
                               fileTracksummary.m_t_px->at(iParticle),          // m_px
                               fileTracksummary.m_t_py->at(iParticle),          // m_py
                               fileTracksummary.m_t_pz->at(iParticle),          // m_pz
                               0,                                               // e
                               fileTracksummary.m_t_vx->at(iParticle),          // m_vx
                               fileTracksummary.m_t_vy->at(iParticle),          // m_vy
                               fileTracksummary.m_t_vz->at(iParticle),          // m_vz
                               fileTracksummary.m_t_time->at(iParticle));       // m_vt

        // Extract ACTS track parameters
        const float phi = fileTracksummary.m_ePHI_fit->at(iTrack);
        const float theta = fileTracksummary.m_eTHETA_fit->at(iTrack);
        const float qOverP = fileTracksummary.m_eQOP_fit->at(iTrack);
        const float loc0 = fileTracksummary.m_eLOC0_fit->at(iTrack);
        const float loc1 = fileTracksummary.m_eLOC1_fit->at(iTrack);

        // Convert to ALICE track parameters
        // ALICE uses: alpha, x, y, z, snp, tgl, signed1Pt
        float alpha = phi; // Track angle in global frame
        float x = loc0;    // Local x position
        float y = loc1;    // Local y position
        float z = 0.0f;    // Will be set from DCA or collision vertex

        // Calculate snp (sin of track momentum azimuthal angle)
        float snp = std::sin(phi);

        // Calculate tgl (tangent of track momentum dip angle)
        float tgl = 1.0f / std::tan(theta);

        // Calculate signed1Pt (charge/pt)
        const float m_p = (qOverP != 0) ? std::abs(1.0f / qOverP) : 0.0f;
        const float pt = m_p * std::sin(theta);
        int8_t charge = (qOverP > 0) ? 1 : -1;
        const float signed1Pt = (pt != 0) ? charge / pt : 0.0f;

        // Track quality
        float m_chi2Sum = fileTracksummary.m_chi2Sum->at(iTrack);
        uint32_t m_nMeasurements = fileTracksummary.m_nMeasurements->at(iTrack);
        uint32_t m_NDF = fileTracksummary.m_NDF->at(iTrack);

        // Fill covariance matrices (simplified - should be extracted from ACTS if available)
        float cYY = 0.1f;
        float cZY = 0.0f;
        float cZZ = 0.1f;
        float cSnpY = 0.0f;
        float cSnpZ = 0.0f;
        float cSnpSnp = 0.001f;
        float cTglY = 0.0f;
        float cTglZ = 0.0f;
        float cTglSnp = 0.0f;
        float cTglTgl = 0.001f;
        float c1PtY = 0.0f;
        float c1PtZ = 0.0f;
        float c1PtSnp = 0.0f;
        float c1PtTgl = 0.0f;
        float c1Pt21Pt2 = 0.001f * signed1Pt * signed1Pt;

        // Create TrackParCov object with dummy covariance matrix
        std::array<float, 5> trackParams = {y, z, snp, tgl, signed1Pt};
        std::array<float, 15> trackCov = {cYY, cZY, cZZ, cSnpY, cSnpZ, cSnpSnp,
                                          cTglY, cTglZ, cTglSnp, cTglTgl,
                                          c1PtY, c1PtZ, c1PtSnp, c1PtTgl, c1Pt21Pt2};
        o2::track::TrackParCov trackParCov(x, alpha, trackParams, trackCov, charge);

        // Fill StoredTracks table (basic track parameters)
        tableStoredTracks(tableCollisions.lastIndex(),          // collisionId
                          o2::aod::track::TrackTypeEnum::Track, // trackType
                          trackParCov.getX(),                   // x
                          trackParCov.getAlpha(),               // alpha
                          trackParCov.getY(),                   // y
                          trackParCov.getZ(),                   // z
                          trackParCov.getSnp(),                 // snp
                          trackParCov.getTgl(),                 // tgl
                          trackParCov.getQ2Pt());               // signed1Pt

        // Fill TracksExtension table
        tableTracksExtension(trackParCov.getPt(),
                             trackParCov.getP(),
                             trackParCov.getEta(),
                             trackParCov.getPhi());

        tableStoredTracksCov(std::sqrt(trackParCov.getSigmaY2()),   // SigmaY
                             std::sqrt(trackParCov.getSigmaZ2()),   // SigmaZ
                             std::sqrt(trackParCov.getSigmaSnp2()), // SigmaSnp
                             std::sqrt(trackParCov.getSigmaTgl2()), // SigmaTgl
                             std::sqrt(trackParCov.getSigma1Pt2()), // Sigma1Pt
                             0,                                     // RhoZY
                             0,                                     // RhoSnpY
                             0,                                     // RhoSnpZ
                             0,                                     // RhoTglY
                             0,                                     // RhoTglZ
                             0,                                     // RhoTglSnp
                             0,                                     // Rho1PtY
                             0,                                     // Rho1PtZ
                             0,                                     // Rho1PtSnp
                             0);                                    // Rho1PtTgl

        // covariance matrix at collision vertex
        tableTracksCovExtension(trackParCov.getSigmaY2(),     // sigmaY2
                                trackParCov.getSigmaZY(),     // sigmaZY
                                trackParCov.getSigmaZ2(),     // sigmaZ2
                                trackParCov.getSigmaSnpY(),   // sigmaSnpY
                                trackParCov.getSigmaSnpZ(),   // sigmaSnpZ
                                trackParCov.getSigmaSnp2(),   // sigmaSnp2
                                trackParCov.getSigmaTglY(),   // sigmaTglY
                                trackParCov.getSigmaTglZ(),   // sigmaTglZ
                                trackParCov.getSigmaTglSnp(), // sigmaTglSnp
                                trackParCov.getSigmaTgl2(),   // sigmaTgl2
                                trackParCov.getSigma1PtY(),   // sigma1PtY
                                trackParCov.getSigma1PtZ(),   // sigma1PtZ
                                trackParCov.getSigma1PtSnp(), // sigma1PtSnp
                                trackParCov.getSigma1PtTgl(), // sigma1PtTgl
                                trackParCov.getSigma1Pt2());  // sigma1Pt2

        // Fill MC track labels
        // Get particle linkage from hits using the majority hit index
        int32_t mcParticleId = -1;                         // Default to invalid particle ID
        mcParticleId = tableStoredMcParticles.lastIndex(); // Temporary: link all tracks to the last added MC particle
        // if (fileTracksummary.nMajorityHits && iTrack < fileTracksummary.nMajorityHits->size()) {
        //   unsigned int hitIndex = fileTracksummary.nMajorityHits->at(iTrack);
        //   if (fileHits.barcode && hitIndex < fileHits.barcode->size()) {
        //     mcParticleId = static_cast<int32_t>(fileHits.barcode->at(hitIndex));
        //     LOG(debug) << "Track " << iTrack << " linked to MC particle " << mcParticleId
        //                << " via hit index " << hitIndex;
        //   } else {
        //     LOG(warning) << "Hit index " << hitIndex << " out of range for track " << iTrack
        //                  << " (barcode vector size: " << (fileHits.barcode ? fileHits.barcode->size() : 0) << ")";
        //   }
        // } else {
        //   LOG(warning) << "No majority hit information available for track " << iTrack;
        // }
        // for ( const auto &vv : fileTracksummary.majorityParticleId->at(iTrack) ){
        //   LOG(info) << vv;
        // }
        tableMcTrackLabels(mcParticleId, // McParticleId
                           0);           // mcMask

        // Fill DCA info (simplified - should be calculated properly)
        tableTracksDCA(0.0f,  // dcaXY
                       0.0f); // dcaZ

        tableTracksDCACov(0.0f,  // sigmaDcaXY2
                          0.0f); // sigmaDcaZ2

        // Fill ALICE3 specific tables
        tableTracksAlice3(true); // isReconstructed

        tableTracksExtraA3(m_nMeasurements, // nSiliconHits (using m_nMeasurements as proxy)
                           0);              // nTPCHits

        // Fill extra track info
        tableStoredTracksExtra(0.f,                                 // TPCInnerParam
                               static_cast<uint32_t>(0),            // Flags
                               static_cast<uint8_t>(0),             // ITSClusterSizes
                               static_cast<uint8_t>(0),             // TPCNClsFindable
                               static_cast<uint8_t>(0),             // TPCNClsFindableMinusFound
                               static_cast<int8_t>(0),              // TPCNClsFindableMinusPID
                               static_cast<int8_t>(0),              // TPCNClsFindableMinusCrossedRows
                               static_cast<uint8_t>(0),             // TPCNClsShared
                               static_cast<uint8_t>(0),             // TRDPattern
                               m_chi2Sum / (m_NDF > 0 ? m_NDF : 1), // ITSChi2NCl
                               0.f,                                 // TPCChi2NCl
                               0.f,                                 // TRDChi2
                               0.f,                                 // TOFChi2
                               0.f,                                 // TPCSignal
                               0.f,                                 // TRDSignal
                               0.f,                                 // Length
                               0.f,                                 // TOFExpMom
                               0.f,                                 // TrackEtaEMCAL
                               0.f,                                 // TrackPhiEMCAL
                               0.f,                                 // TrackTime
                               0.f);                                // TrackTimeRes

        // Fill track selection
        tableTrackSelection(false,  // IsGlobalTrackSDD,
                            false,  // TrackCutFlag,
                            false,  // TrackCutFlagFb1,
                            false,  // TrackCutFlagFb2,
                            false,  // TrackCutFlagFb3,
                            false,  // TrackCutFlagFb4,
                            false); // TrackCutFlagFb5,

        tableTrackSelectionExtension(false,  // PassedTrackType,
                                     false,  // PassedPtRange,
                                     false,  // PassedEtaRange,
                                     false,  // PassedTPCNCls,
                                     false,  // PassedTPCCrossedRows,
                                     false,  // PassedTPCCrossedRowsOverNCls,
                                     false,  // PassedTPCChi2NDF,
                                     false,  // PassedTPCRefit,
                                     false,  // PassedITSNCls,
                                     false,  // PassedITSChi2NDF,
                                     false,  // PassedITSRefit,
                                     false,  // PassedITSHits,
                                     false,  // PassedGoldenChi2,
                                     false,  // PassedDCAxy,
                                     false,  // PassedDCAz,
                                     false,  // PassedITSHitsFB1,
                                     false); // PassedITSHitsFB2
      }

      LOG(info) << "Event " << iEvent << ": has " << nTracks << " tracks and " << nParticles << " particles.";
    }
  }
};

o2::framework::WorkflowSpec defineDataProcessing(o2::framework::ConfigContext const& cfgc)
{
  o2::framework::WorkflowSpec w;
  if (cfgc.options().hasOption("aod-file")) {
    std::string inputFile = cfgc.options().get<std::string>("aod-file");
    if (!inputFile.empty()) {
      LOG(info) << "  " << inputFile;
      TString tinputFile(inputFile.c_str());
      inputPath = gSystem->DirName(tinputFile);
    }
  }
  w.push_back(adaptAnalysisTask<Alice3TrackingTranslator>(cfgc));
  return w;
}
