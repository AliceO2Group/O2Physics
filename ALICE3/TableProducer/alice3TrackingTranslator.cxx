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

#include "ALICE3/DataModel/OTFCollision.h"
#include "ALICE3/DataModel/collisionAlice3.h"
#include "ALICE3/DataModel/tracksAlice3.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/DataTypes.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>

#include <TCollection.h>
#include <TFile.h>
#include <TH1.h>
#include <TList.h>
#include <TString.h>
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTree.h>

#include <Rtypes.h>
#include <RtypesCore.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <map>
#include <string>
#include <utility>
#include <vector>

TString inputPath;

using namespace o2::framework;

namespace
{
using Key = std::array<std::uint32_t, 5>; // barcode (vp, vs, particle, gen, sub)
Key keyOf(std::uint32_t vp, std::uint32_t vs, std::uint32_t pa,
          std::uint32_t ge, std::uint32_t sp) { return Key{vp, vs, pa, ge, sp}; }
// struct ParticleTruthInfo { int pdg; float px, py, pz, m; };
//   struct Seg { double px, py, pz, p; };
} // namespace

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
  o2::framework::Produces<o2::aod::TracksAlice3Pdg> tableTracksAlice3Pdg;
  o2::framework::Produces<o2::aod::TracksExtraA3> tableTracksExtraA3;
  o2::framework::Produces<o2::aod::MCParticlesExtraA3> tableMCParticlesExtraA3;

  o2::framework::Produces<o2::aod::StoredTracksExtra_002> tableStoredTracksExtra;
  o2::framework::Produces<o2::aod::TrackSelection> tableTrackSelection;
  o2::framework::Produces<o2::aod::TrackSelectionExtension> tableTrackSelectionExtension;
  o2::framework::Produces<o2::aod::StoredMcParticles> tableStoredMcParticles;
  o2::framework::Produces<o2::aod::McCollisions> tableMcCollisions;
  o2::framework::Produces<o2::aod::OTFLUTConfigId> tableOTFLUTConfigId;

  o2::framework::Configurable<int> maxCollisions{"maxCollisions", -1000, "Maximum number of collisions translated"};
  o2::framework::Configurable<bool> useTrueInfoForRecoTracks{"useTrueInfoForRecoTracks", false, "Use true information for reconstructed tracks"};
  // o2::framework::Configurable<bool> addDaughterInfo{"addDaughterInfo", false, "Add daughter particle information to the MC truth output tables"};

  o2::framework::HistogramRegistry histos{"histos", {}, o2::framework::OutputObjHandlingPolicy::AnalysisObject};
  void init(o2::framework::InitContext&)
  {
    // Initialization if needed
    LOG(info) << "Alice3TrackingTranslator init called";

    histos.add("hPx", "Px distribution;Px [GeV/c];Entries", kTH1F, {{100, -5.0, 5.0}});
    histos.add("hPy", "Py distribution;Py [GeV/c];Entries", kTH1F, {{100, -5.0, 5.0}});
    histos.add("hPz", "Pz distribution;Pz [GeV/c];Entries", kTH1F, {{100, -5.0, 5.0}});
  }

#define SETADDRESS(branchname, branchvar)                            \
  if (mTree->SetBranchAddress(branchname, &branchvar)) {             \
    LOG(fatal) << "Could not set branch address for " << branchname; \
  }

  struct FileStruct {
    FileStruct(const std::string& filename, const std::string& treename) : mFile(filename.c_str(), "READ")
    {
      if (mFile.IsZombie()) {
        LOG(fatal) << "Could not open file '" << filename << "'";
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
    ParticleStruct(const std::string& filename, const std::string& treename) : FileStruct(std::move(filename), std::move(treename))
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
      SETADDRESS("q", m_q);
      SETADDRESS("number_of_hits", m_number_of_hits);
      SETADDRESS("vertex_primary", m_vertex_primary);
      SETADDRESS("vertex_secondary", m_vertex_secondary);
      SETADDRESS("generation", m_generation);
      SETADDRESS("sub_particle", m_sub_particle);
      SETADDRESS("particle", m_particle);
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
    std::vector<int>* m_number_of_hits = nullptr;
    std::vector<ULong64_t>* m_particleId = nullptr;
    std::vector<float>* m_q = nullptr;
    std::vector<std::uint32_t>* m_vertex_primary = nullptr;
    std::vector<std::uint32_t>* m_vertex_secondary = nullptr;
    std::vector<std::uint32_t>* m_generation = nullptr;
    std::vector<std::uint32_t>* m_sub_particle = nullptr;
    std::vector<std::uint32_t>* m_particle = nullptr;
    // std::vector<ULong64_t>* m_motherId = nullptr;
  };

  struct VertexStruct : public FileStruct {
    VertexStruct(const std::string& filename, const std::string& treename) : FileStruct(filename, treename)
    {
      SETADDRESS("vx", m_x);
      SETADDRESS("vy", m_y);
      SETADDRESS("vz", m_z);
      SETADDRESS("vt", m_t);
      SETADDRESS("incoming_particles_vertex_primary", m_incoming_particles_vertex_primary);
      SETADDRESS("incoming_particles_vertex_secondary", m_incoming_particles_vertex_secondary);
      SETADDRESS("incoming_particles_generation", m_incoming_particles_generation);
      SETADDRESS("incoming_particles_sub_particle", m_incoming_particles_sub_particle);
      SETADDRESS("incoming_particles_particle", m_incoming_particles_particle);
      SETADDRESS("outgoing_particles_vertex_primary", m_outgoing_particles_vertex_primary);
      SETADDRESS("outgoing_particles_vertex_secondary", m_outgoing_particles_vertex_secondary);
      SETADDRESS("outgoing_particles_generation", m_outgoing_particles_generation);
      SETADDRESS("outgoing_particles_sub_particle", m_outgoing_particles_sub_particle);
      SETADDRESS("outgoing_particles_particle", m_outgoing_particles_particle);
    }
    std::vector<float>* m_x = nullptr;
    std::vector<float>* m_y = nullptr;
    std::vector<float>* m_z = nullptr;
    std::vector<float>* m_t = nullptr;
    std::vector<std::vector<std::uint32_t>>* m_incoming_particles_vertex_primary = nullptr;
    std::vector<std::vector<std::uint32_t>>* m_incoming_particles_vertex_secondary = nullptr;
    std::vector<std::vector<std::uint32_t>>* m_incoming_particles_generation = nullptr;
    std::vector<std::vector<std::uint32_t>>* m_incoming_particles_sub_particle = nullptr;
    std::vector<std::vector<std::uint32_t>>* m_incoming_particles_particle = nullptr;
    std::vector<std::vector<std::uint32_t>>* m_outgoing_particles_vertex_primary = nullptr;
    std::vector<std::vector<std::uint32_t>>* m_outgoing_particles_vertex_secondary = nullptr;
    std::vector<std::vector<std::uint32_t>>* m_outgoing_particles_generation = nullptr;
    std::vector<std::vector<std::uint32_t>>* m_outgoing_particles_sub_particle = nullptr;
    std::vector<std::vector<std::uint32_t>>* m_outgoing_particles_particle = nullptr;
  };

  struct TrackStruct : public FileStruct {
    TrackStruct(const std::string& filename, const std::string& treename) : FileStruct(std::move(filename), std::move(treename))
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
      SETADDRESS("majorityParticleId_vertex_primary", m_majorityParticleId_vertex_primary);
      SETADDRESS("majorityParticleId_vertex_secondary", m_majorityParticleId_vertex_secondary);
      SETADDRESS("majorityParticleId_generation", m_majorityParticleId_generation);
      SETADDRESS("majorityParticleId_sub_particle", m_majorityParticleId_sub_particle);
      SETADDRESS("majorityParticleId_particle", m_majorityParticleId_particle);
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
    std::vector<unsigned int>* m_nMajorityHits = nullptr; /// The number of hits from majority particle
    std::vector<int>* m_t_charge = nullptr;               /// Charge of majority particle
    std::vector<float>* m_t_time = nullptr;               /// Time of majority particle
    std::vector<float>* m_t_vx = nullptr;                 /// Vertex x positions of majority particle
    std::vector<float>* m_t_vy = nullptr;                 /// Vertex y positions of majority particle
    std::vector<float>* m_t_vz = nullptr;                 /// Vertex z positions of majority particle
    std::vector<float>* m_t_px = nullptr;                 /// Initial momenta m_px of majority particle
    std::vector<float>* m_t_py = nullptr;                 /// Initial momenta m_py of majority particle
    std::vector<float>* m_t_pz = nullptr;                 /// Initial momenta m_pz of majority particle
    std::vector<float>* m_t_theta = nullptr;              /// Initial momenta theta of majority particle
    std::vector<float>* m_t_phi = nullptr;                /// Initial momenta phi of majority particle
    std::vector<float>* m_t_pT = nullptr;                 /// Initial momenta pT of majority particle
    std::vector<float>* m_t_eta = nullptr;                /// Initial momenta eta of majority particle

    std::vector<std::uint32_t>* m_majorityParticleId_vertex_primary = nullptr;
    std::vector<std::uint32_t>* m_majorityParticleId_vertex_secondary = nullptr;
    std::vector<std::uint32_t>* m_majorityParticleId_generation = nullptr;
    std::vector<std::uint32_t>* m_majorityParticleId_sub_particle = nullptr;
    std::vector<std::uint32_t>* m_majorityParticleId_particle = nullptr;
  };

  struct HitsStruct : public FileStruct {
    HitsStruct(const std::string& filename, const std::string& treename) : FileStruct(std::move(filename), std::move(treename))
    {
      mTree->Print();
      SETADDRESS("barcode", barcode);
    }
    std::vector<unsigned int>* barcode = nullptr;
  };
  void addMCParticle(int collIndex, ParticleStruct& fileParticles, int iParticle, uint8_t flags, int firstMother, int firstDaughter, int secondDaughter, int numberOfHits)
  {
    int mothers[2] = {firstMother, firstMother};
    int daughters[2] = {firstDaughter, secondDaughter};
    tableStoredMcParticles(collIndex,                                                                      // mcCollisionId
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
    tableMCParticlesExtraA3(numberOfHits,                                                                  // number_of_hits
                            fileParticles.m_q->at(iParticle));                                             // charge
  };

  void process(o2::aod::BCs const&)
  {
    LOG(info) << "Alice3TrackingTranslator process called";
    // Find all ROOT files in the folder
    std::vector<std::string> rootFiles;
    LOG(info) << "Reading input files from path: " << inputPath.Data();
    TSystemDirectory dir(inputPath.Data(), inputPath.Data());
    TList* filesList = dir.GetListOfFiles();
    if (filesList) {
      TIter next(filesList);
      TSystemFile* file;
      while ((file = static_cast<TSystemFile*>(next()))) {
        TString fname = file->GetName();
        LOG(debug) << "Found file: " << fname.Data();
        if (!file->IsDirectory() && fname.EndsWith(".root")) {
          TString fullPath = TString::Format("%s/%s", inputPath.Data(), fname.Data());
          LOG(info) << "Adding file: " << fullPath.Data();
          rootFiles.push_back(fullPath.Data());
        }
      }
      delete filesList;
    } else {
      LOG(fatal) << "Could not open directory: " << inputPath.Data();
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
    LOG(info) << "All files loaded successfully";
    // Now open the files to translate and read the trees
    // ParticleStruct fileParticles(files["particles.root"], "particles");
    // LOG(info) << "Particles loaded successfully";
    ParticleStruct fileParticlesSim(files["particles_simulation.root"], "particles");
    LOG(info) << "Particles Sim loaded successfully";
    // std::string daughterFileName = addDaughterInfo ? "particles_decay.root" : "particles_simulation.root";
    // ParticleStruct fileDaughterParticles(files[daughterFileName], "particles");
    // LOG(info) << "Daughter particles loaded successfully from file " << daughterFileName;
    // FileStruct fileVertices(files["performance_vertexing.root"], "vertexing");
    TrackStruct fileTracksummary(files["tracksummary_ambi-tracks-merged.root"], "tracksummary");
    // HitsStruct fileHits(files["hits.root"], "hits");
    VertexStruct fileVertices(files["vertices_gen_and_geant.root"], "vertices");

    LOG(info) << "Tracks loaded successfully";
    const Long64_t kEvents = fileParticlesSim.getEntries();
    // int indexOfLastParticleAfterEvent = -1;
    for (Long64_t iEvent = 0; iEvent < kEvents; ++iEvent) {
      if (iEvent > 0 && maxCollisions.value > 0 && (iEvent % maxCollisions) == 0) {
        LOG(info) << "Stopping at event " << iEvent << "/" << kEvents;
        break;
      }
      // fileParticles.setEventEntry(iEvent);
      fileVertices.setEventEntry(iEvent);
      fileTracksummary.setEventEntry(iEvent);
      // fileHits.setEventEntry(iEvent);
      // if (addDaughterInfo)
      //   fileDaughterParticles.setEventEntry(iEvent);
      fileParticlesSim.setEventEntry(iEvent);

      LOG(info) << "Processing event " << iEvent << "/" << kEvents;

      // Create collision entry for this event
      // TODO: Extract proper collision position from vertex file if available
      float collisionX = 0.0f;
      float collisionY = 0.0f;
      float collisionZ = 0.0f;

      tableOTFLUTConfigId(0); // dummy for the moment

      // Determine the collision ID for the new entry.
      // If the table is empty, lastIndex() returns -1, so we start at 0.
      // If it has entries, lastIndex() returns the index of the last element, so we use lastIndex() + 1.
      int collisionId = tableCollisions.lastIndex() + 1;

      // Convert tracks from ACTS to ALICE format
      // const size_t nParticlesGen = fileParticles.m_vx->size();
      const size_t nParticlesSim = fileParticlesSim.m_vx->size();
      // const size_t nDaughterParticles = fileDaughterParticles.m_vx->size();
      const size_t nTracks = fileTracksummary.m_eLOC0_fit->size();
      std::vector<size_t> idMCparticles;

      // local index k within this event -> global AO2D index
      int firstIdxThisEvent = tableStoredMcParticles.lastIndex() + 1;

      // barcode -> global AO2D index (used later for track labels)
      std::map<Key, int> barcodeToGlobalIdx;
      std::map<Key, std::vector<float>> barcodeToVertexPosition;
      std::map<int, int> GlobalIdxToMotherIdx;
      std::map<int, int> GlobalIdxToPDGCode;
      std::map<int, int> GlobalToLocalIdx;
      std::map<int, std::vector<int>> GlobalIdxToDaughterIdxs;
      for (size_t iPart = 0; iPart < nParticlesSim; ++iPart) {
        Key key = keyOf(
          fileParticlesSim.m_vertex_primary->at(iPart),
          fileParticlesSim.m_vertex_secondary->at(iPart),
          fileParticlesSim.m_particle->at(iPart),
          fileParticlesSim.m_generation->at(iPart),
          fileParticlesSim.m_sub_particle->at(iPart));
        barcodeToGlobalIdx[key] = firstIdxThisEvent + (int)iPart;

        if (iPart == 0) {
          tableMcCollisions(0,                                // mccollision::BCId,
                            0,                                // mccollision::GeneratorsID,
                            fileParticlesSim.m_vx->at(iPart), // mccollision::PosX,
                            fileParticlesSim.m_vy->at(iPart), // mccollision::PosY,
                            fileParticlesSim.m_vz->at(iPart), // mccollision::PosZ
                            fileParticlesSim.m_vt->at(iPart), // mccollision::T
                            1.0f,                             // mccollision::Weight
                            0.0f,                             // mccollision::ImpactParameter,
                            0.f);                             // mccollision::EventPlaneAngle,
        }
      }

      for (size_t j = 0; j < fileVertices.m_incoming_particles_vertex_secondary->size(); ++j) {
        for (size_t d = 0; d < (*fileVertices.m_incoming_particles_vertex_secondary)[j].size(); ++d) {
          Key kd = keyOf(
            (*fileVertices.m_incoming_particles_vertex_primary)[j][d],
            (*fileVertices.m_incoming_particles_vertex_secondary)[j][d],
            (*fileVertices.m_incoming_particles_particle)[j][d],
            (*fileVertices.m_incoming_particles_generation)[j][d],
            (*fileVertices.m_incoming_particles_sub_particle)[j][d]);

          auto iteratorMother = barcodeToGlobalIdx.find(kd);
          if (iteratorMother != barcodeToGlobalIdx.end()) {
            int motherIdx = iteratorMother->second;
            for (size_t d2 = 0; d2 < (*fileVertices.m_outgoing_particles_vertex_secondary)[j].size(); ++d2) {
              Key kd2 = keyOf(
                (*fileVertices.m_outgoing_particles_vertex_primary)[j][d2],
                (*fileVertices.m_outgoing_particles_vertex_secondary)[j][d2],
                (*fileVertices.m_outgoing_particles_particle)[j][d2],
                (*fileVertices.m_outgoing_particles_generation)[j][d2],
                (*fileVertices.m_outgoing_particles_sub_particle)[j][d2]);
              barcodeToVertexPosition[kd2].push_back(fileVertices.m_x->at(j));
              barcodeToVertexPosition[kd2].push_back(fileVertices.m_y->at(j));
              barcodeToVertexPosition[kd2].push_back(fileVertices.m_z->at(j));
              auto iteratorDaughter = barcodeToGlobalIdx.find(kd2);
              if (iteratorDaughter != barcodeToGlobalIdx.end()) {
                int daughterIdx = iteratorDaughter->second;
                GlobalIdxToMotherIdx[daughterIdx] = motherIdx;
                GlobalIdxToDaughterIdxs[motherIdx].push_back(daughterIdx);
              }
            }
          }
        }
      }
      // int idPrimary = -1;
      for (size_t iPart = 0; iPart < nParticlesSim; ++iPart) {
        int globalIdx = firstIdxThisEvent + (int)iPart;
        int motherIdx = -1;
        if (GlobalIdxToMotherIdx.find(globalIdx) != GlobalIdxToMotherIdx.end()) {
          motherIdx = GlobalIdxToMotherIdx[globalIdx];
        }
        int firstDaughter = -1;
        int secondDaughter = -1;
        if (GlobalIdxToDaughterIdxs.find(globalIdx) != GlobalIdxToDaughterIdxs.end()) {
          const auto& daughters = GlobalIdxToDaughterIdxs[globalIdx];
          if (daughters.size() > 0) {
            firstDaughter = daughters[0];
          }
          if (daughters.size() > 1) {
            secondDaughter = daughters[1];
          }
        }
        // Determine flags for the MC particle
        // TODO: For now primary are only those without mother, but all based on PhysicalPrimary definition should be considered.
        uint8_t flags = 0;
        if (motherIdx == -1) {
          collisionX = fileVertices.m_x->at(iPart);
          collisionY = fileVertices.m_y->at(iPart);
          collisionZ = fileVertices.m_z->at(iPart);
          flags |= o2::aod::mcparticle::enums::PhysicalPrimary;
          // idPrimary += 1;
        }
        // if(motherIdx==-1 && std::abs(fileParticlesSim.m_particle_type->at(iPart)) == 211) {
        std::cout << "Adding MC particle " << iPart << ", with globalIdx " << globalIdx << ", with motherIdx " << motherIdx << ", firstDaughter " << firstDaughter << ", secondDaughter " << secondDaughter << ", pdg " << fileParticlesSim.m_particle_type->at(iPart) << std::endl;
        //  std::cout << "PV: " << fileVertices.m_x->at(iPart) << ", " << fileVertices.m_y->at(iPart) << ", " << fileVertices.m_z->at(iPart) << std::endl;
        // std::cout << "Daughter info: ";
        // if(firstDaughter != -1){
        //   std::cout << "firstDaughter pdg " << fileParticlesSim.m_particle_type->at(firstDaughter) << " ";
        // }
        // if(secondDaughter != -1){
        //   std::cout << "secondDaughter pdg " << fileParticlesSim.m_particle_type->at(secondDaughter) << " ";
        // }
        // std::cout << std::endl;
        // }
        GlobalIdxToPDGCode[globalIdx] = fileParticlesSim.m_particle_type->at(iPart);
        GlobalToLocalIdx[globalIdx] = iPart;
        addMCParticle(tableMcCollisions.lastIndex(), fileParticlesSim, iPart, flags, motherIdx, firstDaughter, secondDaughter, fileParticlesSim.m_number_of_hits->at(iPart));
      }

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

      // Reconstructed information
      for (size_t iTrack = 0; iTrack < nTracks; ++iTrack) {
        std::cout << "Processing track " << iTrack << "/" << nTracks << std::endl;

        Key key = keyOf(
          fileTracksummary.m_majorityParticleId_vertex_primary->at(iTrack),
          fileTracksummary.m_majorityParticleId_vertex_secondary->at(iTrack),
          fileTracksummary.m_majorityParticleId_particle->at(iTrack),
          fileTracksummary.m_majorityParticleId_generation->at(iTrack),
          fileTracksummary.m_majorityParticleId_sub_particle->at(iTrack));
        int mcParticleIdx = -1;
        auto iterator = barcodeToGlobalIdx.find(key);
        if (iterator != barcodeToGlobalIdx.end()) {
          mcParticleIdx = iterator->second;
        }
        std::cout << "Track " << iTrack << " is associated with MC particle index " << mcParticleIdx << std::endl;
        // Extract ACTS track parameters
        float phi = fileTracksummary.m_ePHI_fit->at(iTrack);
        float theta = fileTracksummary.m_eTHETA_fit->at(iTrack);
        float qOverP = fileTracksummary.m_eQOP_fit->at(iTrack);
        float loc0 = fileTracksummary.m_eLOC0_fit->at(iTrack);
        float loc1 = fileTracksummary.m_eLOC1_fit->at(iTrack);
        float vx = loc0;
        float vy = loc1;
        float vz = 0.0f;
        // auto iteratorVz = barcodeToVertexPosition.find(key);
        // if (iteratorVz != barcodeToVertexPosition.end()){
        //   vx = iteratorVz->second[0];
        //   vy = iteratorVz->second[1];
        //   vz = iteratorVz->second[2];
        // }

        if (useTrueInfoForRecoTracks) {
          if (mcParticleIdx != -1) {
            int localIdx = GlobalToLocalIdx[mcParticleIdx];
            float p = fileParticlesSim.m_p->at(localIdx);
            // float m = fileParticlesSim.m_m->at(localIdx);
            float pt = std::hypot(fileParticlesSim.m_px->at(localIdx), fileParticlesSim.m_py->at(localIdx));
            phi = std::acos(fileParticlesSim.m_px->at(localIdx) / pt);
            theta = std::acos(fileParticlesSim.m_pz->at(localIdx) / p);
            qOverP = (fileParticlesSim.m_q->at(localIdx) != 0) ? (fileParticlesSim.m_q->at(localIdx) / p) : 0.0f;
            loc0 = fileParticlesSim.m_vx->at(localIdx);
            loc1 = fileParticlesSim.m_vy->at(localIdx);
          }
        }

        std::cout << "Track parameters: phi=" << phi << ", theta=" << theta << ", qOverP=" << qOverP << ", loc0=" << loc0 << ", loc1=" << loc1 << std::endl;
        std::cout << "vx= " << vx << ", vy= " << vy << ", vz= " << vz << std::endl;
        // Convert to ALICE track parameters
        // ALICE uses: alpha, x, y, z, snp, tgl, signed1Pt
        float alpha = phi - M_PI / 2; // Track angle in global frame
        float x = vx;                 // Local x position
        float y = vy;                 // Local y position
        float z = vz;                 // Will be set from DCA or collision vertex

        // Calculate snp (sin of track momentum azimuthal angle)
        float snp = o2::constants::math::Almost1; // std::sin(phi);//

        // Calculate tgl (tangent of track momentum dip angle)
        float tgl = 1.0f / std::tan(theta);

        // Calculate signed1Pt (charge/pt)
        const float m_p = (qOverP != 0) ? std::abs(1.0f / qOverP) : 0.0f;
        const float px = m_p * std::cos(phi) * std::sin(theta);
        const float py = m_p * std::sin(phi) * std::sin(theta);
        const float pz = m_p * std::cos(theta);
        const float pt = m_p * std::sin(theta);
        int8_t charge = (qOverP > 0) ? 1 : -1;
        const float signed1Pt = (pt != 0) ? charge / pt : 0.0f;

        if (charge > 0) {
          histos.fill(HIST("hPx"), px);
          histos.fill(HIST("hPy"), py);
          histos.fill(HIST("hPz"), pz);
        }

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
        std::cout << "TrackParCov created with parameters:" << std::endl;
        // Fill StoredTracks table (basic track parameters)
        tableStoredTracks(collisionId,                          // collisionId
                          o2::aod::track::TrackTypeEnum::Track, // trackType
                          trackParCov.getX(),                   // x
                          trackParCov.getAlpha(),               // alpha
                          trackParCov.getY(),                   // y
                          trackParCov.getZ(),                   // z
                          trackParCov.getSnp(),                 // snp
                          trackParCov.getTgl(),                 // tgl
                          trackParCov.getQ2Pt());               // signed1Pt
        std::cout << "Filling StoredTracks table" << std::endl;
        // Fill TracksExtension table
        tableTracksExtension(trackParCov.getPt(),
                             trackParCov.getP(),
                             trackParCov.getEta(),
                             trackParCov.getPhi());
        std::cout << "Filling TracksExtension table" << std::endl;
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
        std::cout << "Filling StoredTracksCov table" << std::endl;
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
                                trackParCov.getSigma1Pt2());  // sigma1Pt
        std::cout << "Filling TracksCovExtension table" << std::endl;
        // Fill MC label
        tableMcTrackLabels(mcParticleIdx, // McParticleId
                           0);            // mcMask
        std::cout << "Filling McTrackLabels table" << std::endl;
        // Fill DCA info (simplified - should be calculated properly)
        tableTracksDCA(0.0f,  // dcaXY
                       0.0f); // dcaZ
        std::cout << "Filling TracksDCA table" << std::endl;
        tableTracksDCACov(0.0f,  // sigmaDcaXY2
                          0.0f); // sigmaDcaZ2
        std::cout << "Filling TracksDCACov table" << std::endl;
        // Fill ALICE3 specific tables
        tableTracksAlice3(true); // isReconstructed
        std::cout << GlobalIdxToPDGCode[mcParticleIdx] << std::endl;
        if (mcParticleIdx > 0)
          tableTracksAlice3Pdg(GlobalIdxToPDGCode[mcParticleIdx]); // PdgCode to the linked MC truth particle
        else
          tableTracksAlice3Pdg(0); // No linked MC truth particle
        std::cout << "Filling TracksAlice3 and TracksAlice3Pdg tables" << std::endl;
        tableTracksExtraA3(m_nMeasurements, // nSiliconHits (using m_nMeasurements as proxy)
                           0,               // nTPCHits
                           0,               // trackType
                           false);          // isPVContributor
        std::cout << "Filling TracksExtraA3 table" << std::endl;
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
        std::cout << "Filling StoredTracksExtra table" << std::endl;
        // Fill track selection
        tableTrackSelection(false,  // IsGlobalTrackSDD,
                            false,  // TrackCutFlag,
                            false,  // TrackCutFlagFb1,
                            false,  // TrackCutFlagFb2,
                            false,  // TrackCutFlagFb3,
                            false,  // TrackCutFlagFb4,
                            false); // TrackCutFlagFb5,
        std::cout << "Filling TrackSelection table" << std::endl;
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
        std::cout << "Filling TrackSelectionExtension table" << std::endl;
      }

      //   for (size_t iTrack = 0; iTrack < nTracks; ++iTrack) {
      //     LOG(info) << "Processing track " << iTrack << "/" << nTracks << " (nParticlesSim=" << nParticlesSim << ") nParticlesGen=" << nParticlesGen;
      //     const size_t iParticle = iTrack;

      //     if (iParticle == 0) {
      //       tableMcCollisions(0,                                 // mccollision::BCId,
      //                         0,                                 // mccollision::GeneratorsID,
      //                         fileParticles.m_vx->at(iParticle), // mccollision::PosX,
      //                         fileParticles.m_vy->at(iParticle), // mccollision::PosY,
      //                         fileParticles.m_vz->at(iParticle), // mccollision::PosZ
      //                         fileParticles.m_vt->at(iParticle), // mccollision::T
      //                         1.0f,                              // mccollision::Weight
      //                         0.0f,                              // mccollision::ImpactParameter,
      //                         0.f);                              // mccollision::EventPlaneAngle,
      //     }
      //     uint8_t flags = 0;

      //     Key majorityVertexInfo = keyOf(fileTracksummary.m_majorityParticleId_vertex_primary->at(iTrack), fileTracksummary.m_majorityParticleId_vertex_secondary->at(iTrack), (*ouPa)[j][d],
      //                    (*ouGe)[j][d], (*ouSp)[j][d]);

      //     // ULong64_t idMCTrueParticle = fileTracksummary.m_majorityParticleId->at(iParticle);
      //     // int32_t mcParticleId = -1;
      //     // int pdgCode = -1;

      //     // for (size_t iMC = 0; iMC < nParticlesGen; ++iMC) {
      //     //   if (fileParticles.m_particleId->at(iMC) == idMCTrueParticle) {
      //     //     if (count(idMCparticles.begin(), idMCparticles.end(), fileParticles.m_particleId->at(iMC)) > 0) {
      //     //       continue;
      //     //     }
      //     //     idMCparticles.push_back(fileParticles.m_particleId->at(iMC));
      //     //     flags |= o2::aod::mcparticle::enums::PhysicalPrimary;
      //     //     int nHits = 0;
      //     //     for (size_t iPartSim = 0; iPartSim < nParticlesSim; ++iPartSim) {
      //     //       if (fileParticlesSim.m_particleId->at(iPartSim) == fileParticles.m_particleId->at(iMC)) {
      //     //         nHits = fileParticlesSim.m_number_of_hits->at(iPartSim);
      //     //         break;
      //     //       }
      //     //     }
      //     //     addMCParticle(tableMcCollisions.lastIndex(), fileParticles, iMC, flags, -1, -1, nHits);
      //     //     mcParticleId = tableStoredMcParticles.lastIndex();
      //     //     pdgCode = fileParticles.m_particle_type->at(iMC);
      //     //     break;
      //     //   }
      //     // }
      //     // if (addDaughterInfo) {
      //       // for (size_t iMC = 0; iMC < nParticlesSim; ++iMC) {
      //     //     if (fileDaughterParticles.m_particleId->at(iMC) == idMCTrueParticle) {
      //     //       if (count(idMCparticles.begin(), idMCparticles.end(), fileDaughterParticles.m_particleId->at(iMC)) > 0) {
      //     //         break;
      //     //       }

      //     //       int nHits = 0;
      //     //       for (size_t iPartSim = 0; iPartSim < nParticlesSim; ++iPartSim) {
      //     //         if (fileParticlesSim.m_particleId->at(iPartSim) == fileDaughterParticles.m_particleId->at(iMC)) {
      //     //           nHits = fileParticlesSim.m_number_of_hits->at(iPartSim);
      //     //           break;
      //     //         }
      //     //       }
      //     //       for (size_t iMother = 0; iMother < nParticlesGen; ++iMother) {
      //     //         if (fileDaughterParticles.m_motherId->at(iMC) == fileParticles.m_particleId->at(iMother)) {
      //     //           if (count(idMCparticles.begin(), idMCparticles.end(), fileParticles.m_particleId->at(iMother)) > 0) {
      //     //             break;
      //     //           }
      //     //           idMCparticles.push_back(fileParticles.m_particleId->at(iMother));
      //     //           uint8_t flagsMother = o2::aod::mcparticle::enums::PhysicalPrimary;
      //     //           addMCParticle(tableMcCollisions.lastIndex(), fileParticles, iMother, flagsMother, -1, tableStoredMcParticles.lastIndex() + 2, 0);
      //     //           break;
      //     //         }
      //     //       }
      //     //       int motherId = -1;
      //     //       if (count(idMCparticles.begin(), idMCparticles.end(), fileDaughterParticles.m_motherId->at(iMC)) > 0) {
      //     //         auto it = find(idMCparticles.begin(), idMCparticles.end(), fileDaughterParticles.m_motherId->at(iMC));
      //     //         motherId = it - idMCparticles.begin() + indexOfLastParticleAfterEvent + 1;
      //     //       }
      //     //       idMCparticles.push_back(fileDaughterParticles.m_particleId->at(iMC));
      //     //       addMCParticle(tableMcCollisions.lastIndex(), fileDaughterParticles, iMC, flags, motherId, -1, nHits);
      //     //       mcParticleId = tableStoredMcParticles.lastIndex();
      //     //       pdgCode = fileDaughterParticles.m_particle_type->at(iMC);
      //     //       break;
      //     //     }
      //       // }
      //     // }
      //

      //

      //

      //     // Fill MC track labels
      //     // Get particle linkage from hits using the majority hit index
      //     // if (fileTracksummary.nMajorityHits && iTrack < fileTracksummary.nMajorityHits->size()) {
      //     //   unsigned int hitIndex = fileTracksummary.nMajorityHits->at(iTrack);
      //     //   if (fileHits.barcode && hitIndex < fileHits.barcode->size()) {
      //     //     mcParticleId = static_cast<int32_t>(fileHits.barcode->at(hitIndex));
      //     //     LOG(debug) << "Track " << iTrack << " linked to MC particle " << mcParticleId
      //     //                << " via hit index " << hitIndex;
      //     //   } else {
      //     //     LOG(warning) << "Hit index " << hitIndex << " out of range for track " << iTrack
      //     //                  << " (barcode vector size: " << (fileHits.barcode ? fileHits.barcode->size() : 0) << ")";
      //     //   }
      //     // } else {
      //     //   LOG(warning) << "No majority hit information available for track " << iTrack;
      //     // }
      //     // for ( const auto &vv : fileTracksummary.majorityParticleId->at(iTrack) ){
      //     //   LOG(info) << vv;
      //     // }
      //
      //   }

      //   for (size_t iParticle = 0; iParticle < nParticlesGen; ++iParticle) {
      //     if (iParticle == 0 && nTracks == 0) {
      //       tableMcCollisions(0,                                 // mccollision::BCId,
      //                         0,                                 // mccollision::GeneratorsID,
      //                         fileParticles.m_vx->at(iParticle), // mccollision::PosX,
      //                         fileParticles.m_vy->at(iParticle), // mccollision::PosY,
      //                         fileParticles.m_vz->at(iParticle), // mccollision::PosZ
      //                         fileParticles.m_vt->at(iParticle), // mccollision::T
      //                         1.0f,                              // mccollision::Weight
      //                         0.0f,                              // mccollision::ImpactParameter,
      //                         0.f);                              // mccollision::EventPlaneAngle,
      //     }
      //     if (idMCparticles.end() != std::find(idMCparticles.begin(), idMCparticles.end(), fileParticles.m_particleId->at(iParticle))) {
      //       // Already added via track
      //       continue;
      //     }
      //     uint8_t flags = 0;
      //     flags |= o2::aod::mcparticle::enums::PhysicalPrimary;

      //     int nHits = 0;
      //     for (size_t iPartSim = 0; iPartSim < nParticlesSim; ++iPartSim) {
      //       if (fileParticlesSim.m_particleId->at(iPartSim) == fileParticles.m_particleId->at(iParticle)) {
      //         nHits = fileParticlesSim.m_number_of_hits->at(iPartSim);
      //         break;
      //       }
      //     }
      //     addMCParticle(tableMcCollisions.lastIndex(), fileParticles, iParticle, flags, -1, -1, nHits);
      //     idMCparticles.push_back(fileParticles.m_particleId->at(iParticle));
      //   }
      //   // if (addDaughterInfo) {
      //   //   for (size_t iParticle = 0; iParticle < nDaughterParticles; ++iParticle) {
      //   //     if (idMCparticles.end() != std::find(idMCparticles.begin(), idMCparticles.end(), fileDaughterParticles.m_particleId->at(iParticle))) {
      //   //       // Already added via track
      //   //       continue;
      //   //     }
      //   //     uint8_t flags = 0;
      //   //     int nHits = 0;
      //   //     for (size_t iPartSim = 0; iPartSim < nParticlesSim; ++iPartSim) {
      //   //       if (fileParticlesSim.m_particleId->at(iPartSim) == fileDaughterParticles.m_particleId->at(iParticle)) {
      //   //         nHits = fileParticlesSim.m_number_of_hits->at(iPartSim);
      //   //         break;
      //   //       }
      //   //     }
      //   //     int motherId = -1;
      //   //     for (size_t iMother = 0; iMother < nParticlesGen; ++iMother) {
      //   //       if (fileDaughterParticles.m_motherId->at(iParticle) == fileParticles.m_particleId->at(iMother)) {
      //   //         if (count(idMCparticles.begin(), idMCparticles.end(), fileDaughterParticles.m_motherId->at(iParticle)) > 0) {
      //   //           auto it = find(idMCparticles.begin(), idMCparticles.end(), fileDaughterParticles.m_motherId->at(iParticle));
      //   //           motherId = it - idMCparticles.begin() + indexOfLastParticleAfterEvent + 1;
      //   //         }
      //   //       }
      //   //     }
      //   //     addMCParticle(tableMcCollisions.lastIndex(), fileDaughterParticles, iParticle, flags, motherId, -1, nHits);
      //   //   }
      //   // }

      //   LOG(info) << "Event " << iEvent << ": has " << nTracks << " tracks, " << nParticlesGen << " particles " << nDaughterParticles << " daughter particles, " << nParticlesSim << " propagated particles.";
      //   LOG(info) << "Total numbers of stored MC particles: " << tableStoredMcParticles.lastIndex() + 1;
      //   indexOfLastParticleAfterEvent = tableStoredMcParticles.lastIndex();
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
