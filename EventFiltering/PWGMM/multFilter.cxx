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

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "../filterTables.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include <cmath>
#include <string>
#include <TMath.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

static const std::vector<std::string> mmObjectsNames{"kHmTrk", "kHmFddFt0cMftFv0", "kHmFddMftFv0", "kHmFv0Mft", "kHmFv0", "kHmMft", "kHfFv0", "kHfMftTrk", "kHfMftFv0Trk", "kHfMftFv0", "kHmMftFt0a", "kHmFt0", "kHfFt0", "kHfMftFt0a", "kHmFt0cFv0", "kHfFt0cFv0", "kHtPt"};

struct multFilter {
  enum { kHighTrackMult = 0,
         kHighFddFt0cMftFv0Mult,
         kHighFddMftFv0Mult,
         kHighFv0MftMult,
         kHighFv0Mult,
         kHighMftMult,
         kHighFv0Flat,
         kHighMftTrkFlat,
         kHighMftFv0TrkFlat,
         kHighMftFv0Flat,
         kHighMftFt0aMult,
         kHighFt0Mult,
         kHighFt0Flat,
         kHighMftFt0aFlat,
         kHighFt0cFv0Mult,
         kHighFt0cFv0Flat,
         kLeadingPtTrack,
         kNtriggersMM };

  // event selection cuts
  Configurable<float> selHTrkMult{"selHTrkMult", 24., "global trk multiplicity threshold"};
  Configurable<float> selHMfddft0cfv0mft{"selHMfddft0cfv0mft", 138.6, "FDD+FV0+FT0C+MFT mult threshold"};
  Configurable<float> selHMfddfv0mft{"selHMfddfv0mft", 137.9, "FDD+FV0+MFT mult threshold"};
  Configurable<float> selHMfv0mft{"selHMfv0mft", 65.6, "FV0+MFT mult threshold"};
  Configurable<float> selHMFv0{"selHMFv0", 9991.85, "FV0-amplitude threshold"};
  Configurable<float> selHMMFTMult{"selHMMFTMult", 28.14, "MFT multilicity threshold"};
  Configurable<float> sel1Ffv0{"sel1Ffv0", 0.9, "1-flatenicity FV0  threshold"};
  Configurable<float> sel1Fmftglob{"sel1Fmftglob", 0.80, "1-flatenicity MFT+Globtracks threshold"};
  Configurable<float> sel1Fmftfv0glob{"sel1Fmftfv0glob", 0.81, "1-flatenicity FV0+MFT+Globtracks threshold"};
  Configurable<float> sel1Fmftfv0{"sel1Fmftfv0", 0.86, "1-flatenicity MFT+FV0 threshold"};
  Configurable<float> sel1Mmftft0a{"sel1Mmftft0a", 41.1, "1-flatenicity MFT+FT0A threshold"};
  Configurable<float> sel1Mft0{"sel1Mft0", 43.4, "FT0 mult threshold"};
  Configurable<float> sel1Fft0{"sel1Fft0", 0.78, "1-flatenicity FT0 threshold"};
  Configurable<float> sel1Fmftft0a{"sel1Fmftft0a", 0.8, "1-flatenicity MFT+FT0A threshold"};
  Configurable<float> sel1Mft0cFv0{"sel1Mft0cfv0", 70.0, "FT0C+FV0 mult threshold"};
  Configurable<float> sel1Fft0cFv0{"sel1Fft0cfv0", 0.84, "1-flatenicity FT0C+FV0 threshold"};
  Configurable<float> selPtTrig{"selPtTrig", 5., "track pT leading threshold"};

  Produces<aod::MultFilters> tags;

  // acceptance cuts
  Configurable<float> cfgTrkEtaCut{"cfgTrkEtaCut", 1.5f, "Eta range for tracks"};
  Configurable<float> cfgTrkLowPtCut{"cfgTrkLowPtCut", 0.15f, "Minimum constituent pT"};

  HistogramRegistry multiplicity{"multiplicity", {}, OutputObjHandlingPolicy::AnalysisObject, true, true};

  static constexpr std::string_view nhEst_before[17] = {
    "eGlobaltrack", "eFDDAFDDCFT0CFV0MFT", "eFDDAFDDCFV0MFT", "eFV0MFT", "eFV0", "eMFTmult",
    "e1flatencityFV0", "e1flatencitytrkMFT", "e1flatencitytrkMFTFV0", "e1flatencityMFTFV0",
    "eMFTmultFT0A", "eFT0", "e1flatencityFT0", "e1flatencityMFTFT0A", "eFT0CFV0", "e1flatencityFT0CFV0", "ePtl"};
  static constexpr std::string_view nhEst_after[17] = {
    "eGlobaltrack_selected", "eFDDAFDDCFT0CFV0MFT_selected", "eFDDAFDDCFV0MFT_selected", "eFV0MFT_selected", "eFV0_selected", "eMFTmult_selected", "e1flatencityFV0_selected", "e1flatencitytrkMFT_selected", "e1flatencitytrkMFTFV0_selected", "e1flatencityMFTFV0_selected", "eMFTmultFT0A_selected", "eFT0_selected", "e1flatencityFT0_selected", "e1flatencityMFTFT0A_selected", "eFT0CFV0_selected", "e1flatencityFT0CFV0_selected", "ePtl_selected"};
  static constexpr std::string_view npEst[17] = {
    "epGlobaltrack", "epFDDAFDDCFT0CFV0MFT", "epFDDAFDDCFV0MFT", "epFV0MFT", "epFV0", "epMFTmult",
    "ep1flatencityFV0", "ep1flatencitytrkMFT", "ep1flatencitytrkMFTFV0", "ep1flatencityMFTFV0",
    "epMFTmultFT0A", "epFT0", "ep1flatencityFT0", "e1pflatencityMFTFT0A", "epFT0CFV0", "ep1flatencityFT0CFV0", "epPtl"};

  static constexpr std::string_view tEst[17] = {
    "GlobalTrk", "FDDA_FDDC_FT0C_FV0_MFT", "FDDA_FDDC_FV0_MFT", "FV0_MFT", "FV0",
    "MFTTrk", "1-flatencity_FV0", "1-flatencity_trk_MFT", "1-flatencity_trk_MFT_FV0",
    "1-flatencity_MFT_FV0", "MFT_FT0A", "FT0", "1-flatencityFT0", "1-flatencity_MFT_FT0A", "FT0C_FV0", "1-flatencity_FT0C_FV0", "pT^{trig} (GeV/#it{c})"};

  void init(o2::framework::InitContext&)
  {
    int nBinsEst[17] = {100, 400, 400, 400, 500, 200, 102, 102, 102, 102, 400, 100, 102, 102, 200, 102, 150};
    float lowEdgeEst[17] = {-0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
                            -0.01, -0.01, -0.01, -0.01, -0.5, -0.5, -0.01, -0.01, -0.5, -0.01, .0};
    float upEdgeEst[17] = {99.5, 399.5, 399.5, 399.5, 39999.5, 199.5,
                           1.01, 1.01, 1.01, 1.01, 399.5, 99.5, 1.01, 1.01, 199.5, 1.01, 150.0};

    // QA event level
    multiplicity.add("fCollZpos", "Vtx_z", HistType::kTH1F, {{200, -20., +20., "#it{z}_{vtx} position (cm)"}});

    // QA FT0
    multiplicity.add("hAmpT0AVsCh", "", HistType::kTH2F, {{24, -0.5, 23.5, "ch"}, {600, -0.5, +5999.5, "FT0A amplitude"}});
    multiplicity.add("hAmpT0CVsCh", "", HistType::kTH2F, {{28, -0.5, 27.5, "ch"}, {600, -0.5, +5999.5, "FT0C amplitude"}});
    multiplicity.add("hFT0C", "FT0C", HistType::kTH1F, {{600, -0.5, 599.5, "FT0C amplitudes"}});
    multiplicity.add("hFT0A", "FT0A", HistType::kTH1F, {{600, -0.5, 599.5, "FT0A amplitudes"}});
    multiplicity.add("hAmpT0AvsVtx", "", HistType::kTH2F, {{30, -15.0, +15.0, "Vtx_z"}, {600, -0.5, +5999.5, "FT0A amplitude"}});
    multiplicity.add("hAmpT0CvsVtx", "", HistType::kTH2F, {{30, -15.0, +15.0, "Vtx_z"}, {600, -0.5, +5999.5, "FT0C amplitude"}});

    // QA FDD
    multiplicity.add("hAmpFDAVsCh", "", HistType::kTH2F, {{8, -0.5, 7.5, "channel"}, {600, -0.5, +599.5, "FDA amplitude"}});
    multiplicity.add("hAmpFDCVsCh", "", HistType::kTH2F, {{8, -0.5, 7.5, "channel"}, {600, -0.5, +599.5, "FDC amplitude"}});
    multiplicity.add("hAmpFDAvsVtx", "", HistType::kTH2F, {{30, -15.0, +15.0, "Vtx_z"}, {6000, -0.5, 7999.5, "Ampl. FDA"}});
    multiplicity.add("hAmpFDCvsVtx", "", HistType::kTH2F, {{30, -15.0, +15.0, "Vtx_z"}, {6000, -0.5, 7999.5, "Ampl. FDC"}});

    // QA global tracks
    multiplicity.add("hdNdetaGlobal", "dNdeta", HistType::kTH1F, {{50, -5.0, 5.0, " "}});
    multiplicity.add("hPhiGlobal", "Phi", HistType::kTH1F, {{64, 0., 2.0 * M_PI, " "}});
    multiplicity.add("hMFTvsVtx", "", HistType::kTH2F, {{30, -15.0, +15.0, "Vtx_z"}, {200, -0.5, +199.5, "MFT mult (-3.6<#eta<-2.5)"}});

    // QA MFT tracks
    multiplicity.add("hdNdetaMFT", "dNdeta (MFTtracks)", HistType::kTH1F, {{50, -5.0, 5.0, " "}});
    multiplicity.add("hPhiMFT", "Phi (MFTtracks)", HistType::kTH1F, {{64, 0., 2.0 * M_PI, " "}});

    // estimators
    for (int i_e = 0; i_e < 17; ++i_e) {
      multiplicity.add(
        npEst[i_e].data(), "", HistType::kTProfile, {{nBinsEst[i_e], lowEdgeEst[i_e], upEdgeEst[i_e], tEst[i_e].data()}});
    }
    for (int i_e = 0; i_e < 17; ++i_e) {
      multiplicity.add(
        nhEst_before[i_e].data(), "", HistType::kTH1F, {{nBinsEst[i_e], lowEdgeEst[i_e], upEdgeEst[i_e], tEst[i_e].data()}});
      multiplicity.add(
        nhEst_after[i_e].data(), "", HistType::kTH1F, {{nBinsEst[i_e], lowEdgeEst[i_e], upEdgeEst[i_e], tEst[i_e].data()}});
    }

    std::array<std::string, 2> eventTitles = {"all", "rejected"};

    auto scalers{std::get<std::shared_ptr<TH1>>(multiplicity.add("fProcessedEvents", "Multiplicity - event filtered;;events", HistType::kTH1F, {{kNtriggersMM + 2, -0.5, kNtriggersMM + 2 - 0.5}}))};
    for (size_t iBin = 0; iBin < eventTitles.size() + mmObjectsNames.size(); iBin++) {
      if (iBin < 2)
        scalers->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
      else
        scalers->GetXaxis()->SetBinLabel(iBin + 1, mmObjectsNames[iBin - 2].data());
    }
    // overlap with HtrackMult
    auto scalerso1{std::get<std::shared_ptr<TH1>>(multiplicity.add("fProcessedEvents_overlap1", "Multiplicity - event filtered;;events", HistType::kTH1F, {{kNtriggersMM + 2, -0.5, kNtriggersMM + 2 - 0.5}}))};
    for (size_t iBin = 0; iBin < eventTitles.size() + mmObjectsNames.size(); iBin++) {
      if (iBin < 2)
        scalerso1->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
      else
        scalerso1->GetXaxis()->SetBinLabel(iBin + 1, mmObjectsNames[iBin - 2].data());
    }
    auto scalerso2{std::get<std::shared_ptr<TH1>>(multiplicity.add("fProcessedEvents_overlap2", "Multiplicity - event filtered;;events", HistType::kTH1F, {{kNtriggersMM + 2, -0.5, kNtriggersMM + 2 - 0.5}}))};
    for (size_t iBin = 0; iBin < eventTitles.size() + mmObjectsNames.size(); iBin++) {
      if (iBin < 2)
        scalerso2->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
      else
        scalerso2->GetXaxis()->SetBinLabel(iBin + 1, mmObjectsNames[iBin - 2].data());
    }
    auto scalerso3{std::get<std::shared_ptr<TH1>>(multiplicity.add("fProcessedEvents_overlap3", "Multiplicity - event filtered;;events", HistType::kTH1F, {{kNtriggersMM + 2, -0.5, kNtriggersMM + 2 - 0.5}}))};
    for (size_t iBin = 0; iBin < eventTitles.size() + mmObjectsNames.size(); iBin++) {
      if (iBin < 2)
        scalerso3->GetXaxis()->SetBinLabel(iBin + 1, eventTitles[iBin].data());
      else
        scalerso3->GetXaxis()->SetBinLabel(iBin + 1, mmObjectsNames[iBin - 2].data());
    }
  }

  Filter trackFilter = (nabs(aod::track::eta) < cfgTrkEtaCut) && (aod::track::pt > cfgTrkLowPtCut);
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection>>;

  int getT0ASector(int i_ch)
  {
    int i_sec_t0a = -1;
    for (int i_sec = 0; i_sec < 24; ++i_sec) {
      if (i_ch >= 4 * i_sec && i_ch <= 3 + 4 * i_sec) {
        i_sec_t0a = i_sec;
        break;
      }
    }
    return i_sec_t0a;
  }
  int getT0CSector(int i_ch)
  {
    int i_sec_t0c = -1;
    for (int i_sec = 0; i_sec < 28; ++i_sec) {
      if (i_ch >= 4 * i_sec && i_ch <= 3 + 4 * i_sec) {
        i_sec_t0c = i_sec;
        break;
      }
    }
    return i_sec_t0c;
  }
  float GetFlatenicity(float signals[], int entries)
  {
    float flat = 9999;
    float mRho = 0;
    for (int iCell = 0; iCell < entries; ++iCell) {
      mRho += 1.0 * signals[iCell];
    }
    // average activity per cell
    mRho /= (1.0 * entries);
    if (mRho <= 0) {
      return 9999;
    }
    // get sigma
    float sRho_tmp = 0;
    for (int iCell = 0; iCell < entries; ++iCell) {
      sRho_tmp += TMath::Power(1.0 * signals[iCell] - mRho, 2);
    }
    sRho_tmp /= (1.0 * entries * entries);
    float sRho = TMath::Sqrt(sRho_tmp);
    if (mRho > 0) {
      flat = sRho / mRho;
    }
    return flat;
  }
  void process(soa::Join<aod::Collisions, aod::EvSels>::iterator const& collision, TrackCandidates const& tracks, aod::MFTTracks const& mfttracks, aod::FT0s const& ft0s, aod::FV0As const& fv0s, aod::FDDs_001 const& fdds_001)
  {

    bool keepEvent[kNtriggersMM]{false};
    auto vtxZ = collision.posZ();
    multiplicity.fill(HIST("fProcessedEvents"), 0);
    multiplicity.fill(HIST("fCollZpos"), collision.posZ());
    // global observables
    int multTrack = 0;
    float flPt = 0; // leading pT
    double flatenicity_fv0 = 9999;

    float sumAmpFV0 = 0;
    float sumAmpFV01to4Ch = 0;
    int innerFV0 = 32;
    const int nCells = 48; // 48 sectors in FV0

    if (collision.has_foundFV0()) {
      float RhoLattice[nCells];
      for (Int_t iCh = 0; iCh < nCells; iCh++) {
        RhoLattice[iCh] = 0;
      }

      auto fv0 = collision.foundFV0();
      // LOGP(info, "amplitude.size()={}", fv0.amplitude().size());
      for (std::size_t ich = 0; ich < fv0.amplitude().size(); ich++) {

        int channelv0 = fv0.channel()[ich];
        float ampl_ch = fv0.amplitude()[ich];
        sumAmpFV0 += ampl_ch;
        if (channelv0 >= 8) { // exclude the 1st ch, eta 2.2,4.52
          sumAmpFV01to4Ch += ampl_ch;
        }
        if (channelv0 < innerFV0) {
          RhoLattice[channelv0] = ampl_ch;
        } else {
          RhoLattice[channelv0] = ampl_ch / 2.0; // two channels per bin
        }
      }

      flatenicity_fv0 = GetFlatenicity(RhoLattice, nCells);
    }

    // FT0
    float sumAmpFT0A = 0.f;
    float sumAmpFT0C = 0.f;
    const int nCellsT0A = 24;
    float RhoLatticeT0A[nCellsT0A];
    for (int iCh = 0; iCh < nCellsT0A; iCh++) {
      RhoLatticeT0A[iCh] = 0.0;
    }
    const int nCellsT0C = 28;
    float RhoLatticeT0C[nCellsT0C];
    for (int iCh = 0; iCh < nCellsT0C; iCh++) {
      RhoLatticeT0C[iCh] = 0.0;
    }

    if (collision.has_foundFT0()) {
      auto ft0 = collision.foundFT0();
      for (std::size_t i_a = 0; i_a < ft0.amplitudeA().size(); i_a++) {
        float amplitude = ft0.amplitudeA()[i_a];
        uint8_t channel = ft0.channelA()[i_a];
        int sector = getT0ASector(channel);
        if (sector >= 0 && sector < 24) {
          RhoLatticeT0A[sector] += amplitude;
          multiplicity.fill(HIST("hAmpT0AVsCh"), sector, amplitude);
        }
        sumAmpFT0A += amplitude;
        multiplicity.fill(HIST("hFT0A"), amplitude);
      }
      for (std::size_t i_c = 0; i_c < ft0.amplitudeC().size(); i_c++) {
        float amplitude = ft0.amplitudeC()[i_c];
        sumAmpFT0C += amplitude;
        uint8_t channel = ft0.channelC()[i_c];
        int sector = getT0CSector(channel);
        if (sector >= 0 && sector < 28) {
          RhoLatticeT0C[sector] += amplitude;
          multiplicity.fill(HIST("hAmpT0CVsCh"), sector, amplitude);
        }
        multiplicity.fill(HIST("hFT0C"), amplitude);
      }
      multiplicity.fill(HIST("hAmpT0AvsVtx"), vtxZ, sumAmpFT0A);
      multiplicity.fill(HIST("hAmpT0CvsVtx"), vtxZ, sumAmpFT0C);
    }
    float flatenicity_t0a = GetFlatenicity(RhoLatticeT0A, nCellsT0A);
    float flatenicity_t0c = GetFlatenicity(RhoLatticeT0C, nCellsT0C);

    // FDD
    float sumAmpFDDA = 0;
    float sumAmpFDDC = 0;
    if (collision.has_foundFDD()) {
      auto fdd = collision.foundFDD();
      for (std::size_t ich = 0; ich < 8; ich++) {
        float amplitude = fdd.chargeA()[ich];
        sumAmpFDDA += amplitude;
        multiplicity.fill(HIST("hAmpFDAVsCh"), ich, amplitude);
      }
      for (std::size_t ich = 0; ich < 8; ich++) {
        float amplitude = fdd.chargeC()[ich];
        sumAmpFDDC += amplitude;
        multiplicity.fill(HIST("hAmpFDCVsCh"), ich, amplitude);
      }
      multiplicity.fill(HIST("hAmpFDAvsVtx"), vtxZ, sumAmpFDDA);
      multiplicity.fill(HIST("hAmpFDCvsVtx"), vtxZ, sumAmpFDDC);
    }

    // MFTtracks
    float multMFTTrackParc = 0;
    float multMFTTrack = 0;

    const int nRings1 = 2;
    const int nSectors1 = 8;
    const int nCells1 = nRings1 * nSectors1;
    float maxEta1[nRings1] = {-3.05, -2.50};
    float minEta1[nRings1] = {-3.60, -3.05};
    float maxPhi1[nSectors1] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    float minPhi1[nSectors1] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    float RhoLattice1[nCells1];
    for (int iCh = 0; iCh < nCells1; iCh++) {
      RhoLattice1[iCh] = 0.0;
    }

    for (auto& track : mfttracks) {
      float eta_a = track.eta();
      float phi_a = track.phi();
      o2::math_utils::bringTo02Pi(phi_a);
      if (eta_a > -2.5 || eta_a < -3.6) { // the MFT eta coverage
        continue;
      }
      multiplicity.fill(HIST("hdNdetaMFT"), eta_a);
      multiplicity.fill(HIST("hPhiMFT"), phi_a);

      int i_ch = 0;
      for (int ir = 0; ir < nRings1; ir++) {
        for (int is = 0; is < nSectors1; is++) {
          if (eta_a >= minEta1[ir] && eta_a < maxEta1[ir] &&
              phi_a >= minPhi1[is] * 2.0 * M_PI / (1.0 * nSectors1) &&
              phi_a < maxPhi1[is] * 2.0 * M_PI / (1.0 * nSectors1)) {
            RhoLattice1[i_ch]++;
          }
          i_ch++;
        }
      }

      multMFTTrack++;

      if (eta_a > -3.4 || eta_a < -3.6) {
        continue;
      }
      multMFTTrackParc++;
    }
    multiplicity.fill(HIST("hMFTvsVtx"), vtxZ, multMFTTrack);
    float flatenicity_mft = GetFlatenicity(RhoLattice1, nCells1);
    // Globaltracks
    const int nRings2 = 4;
    const int nSectors2 = 8;
    const int nCells2 = nRings2 * nSectors2;
    float maxEta2[nRings2] = {-0.4, 0.0, +0.4, +0.8};
    float minEta2[nRings2] = {-0.8, -0.4, +0.0, +0.4};
    float maxPhi2[nSectors2] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
    float minPhi2[nSectors2] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    float RhoLattice2[nCells2];
    for (int iCh = 0; iCh < nCells2; iCh++) {
      RhoLattice2[iCh] = 0.0;
    }

    for (auto& track : tracks) {
      if (!track.isGlobalTrack()) {
        continue;
      }
      float eta_a = track.eta();
      float phi_a = track.phi();
      float pt_a = track.pt();
      multTrack++;
      multiplicity.fill(HIST("hdNdetaGlobal"), eta_a);
      multiplicity.fill(HIST("hPhiGlobal"), phi_a);
      if (flPt < pt_a) {
        flPt = pt_a;
      }
      int i_ch = 0;
      for (int ir = 0; ir < nRings2; ir++) {
        for (int is = 0; is < nSectors2; is++) {
          if (eta_a >= minEta2[ir] && eta_a < maxEta2[ir] &&
              phi_a >= minPhi2[is] * 2.0 * M_PI / (1.0 * nSectors2) &&
              phi_a < maxPhi2[is] * 2.0 * M_PI / (1.0 * nSectors2)) {
            RhoLattice2[i_ch]++;
          }
          i_ch++;
        }
      }
    }
    float flatenicity_glob = GetFlatenicity(RhoLattice2, nCells2);

    float combined_estimator1 = 0;
    float combined_estimator2 = 0;
    float combined_estimator3 = 0;
    float combined_estimator4 = 0;
    float combined_estimator5 = 0;
    float combined_estimator6 = 0;
    float estimator[17];
    float cut[17];
    for (int i_e = 0; i_e < 17; ++i_e) {
      estimator[i_e] = 0;
      cut[i_e] = 0;
    }
    cut[0] = selHTrkMult;
    cut[1] = selHMfddft0cfv0mft;
    cut[2] = selHMfddfv0mft;
    cut[3] = selHMfv0mft;
    cut[4] = selHMFv0;
    cut[5] = selHMMFTMult;
    cut[6] = sel1Ffv0;
    cut[7] = sel1Fmftglob;
    cut[8] = sel1Fmftfv0glob;
    cut[9] = sel1Fmftfv0;
    cut[10] = sel1Mmftft0a;
    cut[11] = sel1Mft0;
    cut[12] = sel1Fft0;
    cut[13] = sel1Fmftft0a;
    cut[14] = sel1Mft0cFv0;
    cut[15] = sel1Fft0cFv0;
    cut[16] = selPtTrig;
    // option 1
    const int nEta1 = 5; // FDDC + MFTparc + FT0C + FV0 (rings 1-4) + FDDA
    float weigthsEta1[nEta1] = {0.0117997, 1.66515, 0.0569502, 0.00548221, 0.0037175};
    float ampl1[nEta1] = {0, 0, 0, 0, 0};
    ampl1[0] = sumAmpFDDC;
    ampl1[1] = multMFTTrackParc;
    ampl1[2] = sumAmpFT0C;
    ampl1[3] = sumAmpFV01to4Ch;
    ampl1[4] = sumAmpFDDA;

    for (int i_1 = 0; i_1 < nEta1; ++i_1) {
      combined_estimator1 += ampl1[i_1] * weigthsEta1[i_1];
    }
    // option 2
    const int nEta2 = 4; // FDDC + MFT + FV0 (rings 1-4) + FDDA
    float weigthsEta2[nEta2] = {0.0117997, 1.05258, 0.00548221, 0.0037175};
    float ampl2[nEta2] = {0, 0, 0, 0};
    ampl2[0] = sumAmpFDDC;
    ampl2[1] = multMFTTrack;
    ampl2[2] = sumAmpFV01to4Ch;
    ampl2[3] = sumAmpFDDA;
    for (int i_2 = 0; i_2 < nEta2; ++i_2) {
      combined_estimator2 += ampl2[i_2] * weigthsEta2[i_2];
    }

    // option 3
    const int nEta3 = 2; // MFT + FV0
    float weigthsEta3[nEta3] = {1.05258, 0.00535717};
    float ampl3[nEta3] = {0, 0};
    ampl3[0] = multMFTTrack;
    ampl3[1] = sumAmpFV0;
    for (int i_3 = 0; i_3 < nEta3; ++i_3) {
      combined_estimator3 += ampl3[i_3] * weigthsEta3[i_3];
    }

    // option 4
    const int nEta4 = 2; // MFT + FT0A
    float weigthsEta4[nEta4] = {1.05258, 0.014552069};
    float ampl4[nEta4] = {0, 0};
    ampl4[0] = multMFTTrack;
    ampl4[1] = sumAmpFT0A;
    for (int i_4 = 0; i_4 < nEta4; ++i_4) {
      combined_estimator4 += ampl4[i_4] * weigthsEta4[i_4];
    }

    // option 5
    const int nEta5 = 2; // FT0C + FT0A
    float weigthsEta5[nEta5] = {0.0569502, 0.014552069};
    float ampl5[nEta5] = {0, 0};
    ampl5[0] = sumAmpFT0C;
    ampl5[1] = sumAmpFT0A;
    for (int i_5 = 0; i_5 < nEta5; ++i_5) {
      combined_estimator5 += ampl5[i_5] * weigthsEta5[i_5];
    }

    // option 6
    const int nEta6 = 2; //  FT0C + FV0
    float weigthsEta6[nEta6] = {0.0569502, 0.00535717};
    float ampl6[nEta6] = {0, 0};
    ampl6[0] = sumAmpFT0C;
    ampl6[1] = sumAmpFV0;
    for (int i_6 = 0; i_6 < nEta6; ++i_6) {
      combined_estimator6 += ampl6[i_6] * weigthsEta6[i_6];
    }

    float flatenicity_mft_glob = (flatenicity_mft + flatenicity_glob) / 2.0;
    float flatenicity_mft_fv0 = (flatenicity_mft + flatenicity_fv0) / 2.0;
    float flatenicity_mft_glob_fv0 =
      (flatenicity_mft + flatenicity_glob + flatenicity_fv0) / 3.0;

    estimator[0] = multTrack;
    estimator[1] = combined_estimator1;
    estimator[2] = combined_estimator2;
    estimator[3] = combined_estimator3;
    estimator[4] = sumAmpFV0;
    estimator[5] = multMFTTrack;
    estimator[6] = 1.0 - flatenicity_fv0;
    estimator[7] = 1.0 - flatenicity_mft_glob;
    estimator[8] = 1.0 - flatenicity_mft_glob_fv0;
    estimator[9] = 1.0 - flatenicity_mft_fv0;
    estimator[10] = combined_estimator4;
    estimator[11] = combined_estimator5;
    float flatenicity_ft0 = (flatenicity_t0a + flatenicity_t0c) / 2.0;
    estimator[12] = 1.0 - flatenicity_ft0;
    float flatenicity_mft_ft0a = (flatenicity_mft + flatenicity_t0a) / 2.0;
    estimator[13] = 1.0 - flatenicity_mft_ft0a;
    estimator[14] = combined_estimator6;
    estimator[15] = 1.0 - (flatenicity_fv0 + flatenicity_t0c) / 2.0;
    estimator[16] = flPt;

    static_for<0, 16>([&](auto i) {
      constexpr int index = i.value;
      multiplicity.fill(HIST(npEst[index]), estimator[index], estimator[0]);
      multiplicity.fill(HIST(nhEst_before[index]), estimator[index]);
    });

    static_for<0, 16>([&](auto i) {
      constexpr int index = i.value;
      if (estimator[index] > cut[index]) {
        multiplicity.fill(HIST(nhEst_after[index]), estimator[index]);
        keepEvent[index] = true;
      }
    });

    tags(keepEvent[kHighTrackMult], keepEvent[kHighFddFt0cMftFv0Mult], keepEvent[kHighFddMftFv0Mult], keepEvent[kHighFv0MftMult], keepEvent[kHighFv0Mult], keepEvent[kHighMftMult], keepEvent[kHighMftMult], keepEvent[kHighMftTrkFlat], keepEvent[kHighMftFv0TrkFlat], keepEvent[kHighMftFv0Flat], keepEvent[kHighMftFt0aMult], keepEvent[kHighFt0Mult], keepEvent[kHighFt0Flat], keepEvent[kHighMftFt0aFlat], keepEvent[kHighFt0cFv0Mult], keepEvent[kHighFt0cFv0Flat], keepEvent[kLeadingPtTrack]);

    if (!keepEvent[kHighTrackMult] && !keepEvent[kHighFddFt0cMftFv0Mult] && !keepEvent[kHighFddMftFv0Mult] && !keepEvent[kHighFv0MftMult] && !keepEvent[kHighFv0Mult] && !keepEvent[kHighMftMult] && !keepEvent[kHighMftMult] && !keepEvent[kHighMftTrkFlat] && !keepEvent[kHighMftFv0TrkFlat] && !keepEvent[kHighMftFv0Flat] && !keepEvent[kHighMftFt0aMult] && !keepEvent[kHighFt0Mult] && !keepEvent[kHighFt0Flat] && !keepEvent[kHighMftFt0aFlat] && !keepEvent[kHighFt0cFv0Mult] && !keepEvent[kHighFt0cFv0Flat] && !keepEvent[kLeadingPtTrack]) {
      multiplicity.fill(HIST("fProcessedEvents"), 1);
    } else {
      for (int iTrigger{0}; iTrigger < kNtriggersMM; iTrigger++) {
        if (keepEvent[iTrigger]) {
          multiplicity.fill(HIST("fProcessedEvents"), iTrigger + 2);
        }
      }
      if (keepEvent[kHighTrackMult]) {
        for (int iTrigger{0}; iTrigger < kNtriggersMM; iTrigger++) {
          if (keepEvent[iTrigger]) {
            multiplicity.fill(HIST("fProcessedEvents_overlap1"), iTrigger + 2);
          }
        }
      }
      if (keepEvent[kHighFv0MftMult]) {
        for (int iTrigger{0}; iTrigger < kNtriggersMM; iTrigger++) {
          if (keepEvent[iTrigger]) {
            multiplicity.fill(HIST("fProcessedEvents_overlap2"), iTrigger + 2);
          }
        }
      }
      if (keepEvent[kHighFt0Mult]) {
        for (int iTrigger{0}; iTrigger < kNtriggersMM; iTrigger++) {
          if (keepEvent[iTrigger]) {
            multiplicity.fill(HIST("fProcessedEvents_overlap3"), iTrigger + 2);
          }
        }
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfg)
{
  return WorkflowSpec{adaptAnalysisTask<multFilter>(cfg)};
}
