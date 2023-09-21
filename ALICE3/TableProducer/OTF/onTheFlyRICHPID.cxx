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

//
// Task to add a table of track parameters propagated to the primary vertex
//

#include <utility>
#include <cmath>
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Framework/ASoAHelpers.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "ReconstructionDataFormats/DCA.h"
#include "DetectorsBase/Propagator.h"
#include "DetectorsBase/GeometryManager.h"
#include "CommonUtils/NameConf.h"
#include "CCDB/CcdbApi.h"
#include "CCDB/BasicCCDBManager.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsCalibration/MeanVertexObject.h"
#include "CommonConstants/GeomConstants.h"
#include "CommonConstants/PhysicsConstants.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TString.h"
#include "ALICE3/DataModel/OTFRICH.h"
#include "DetectorsVertexing/HelixHelper.h"

#include "TableHelper.h"
#include "ALICE3/Core/DelphesO2TrackSmearer.h"

/// \file onTheFlyRichPid.cxx
///
/// \brief This task goes straight from a combination of track table and mcParticles
/// and a projective bRICH configuration to a table of TOF NSigmas for the particles
/// being analysed. It currently contemplates 5 particle types:
/// electrons, pions, kaons, protons and muons.
///
/// More particles could be added but would have to be added to the LUT
/// being used in the onTheFly tracker task.
///
/// \warning Geometry parameters are configurable, but resolution values should be adapted.
/// Since angular resolution depends on the specific geometric details, it is better to
/// calculate it from full simulation and add new input. Alternatively, an analytical
/// expression can be provided as a function of the main parameters.
///
/// \author David Dobrigkeit Chinellato, UNICAMP, Nicola Nicassio, University and INFN Bari

using namespace o2;
using namespace o2::framework;

struct OnTheFlyRichPid {
  Produces<aod::UpgradeRich> upgradeRich;

  // necessary for particle charges
  Service<o2::framework::O2DatabasePDG> pdg;

  // master setting: magnetic field
  Configurable<float> dBz{"dBz", 20, "magnetic field (kilogauss)"};

  // add rich-specific configurables here
  Configurable<int> bRichNumberOfSectors{"bRichNumberOfSectors", 21, "barrel RICH number of sectors"};
  Configurable<float> bRichPhotodetectorCentralModuleHalfLength{"bRichPhotodetectorCentralModuleHalfLength", 18.4 / 2.0, "barrel RICH photodetector central module half length (cm)"};
  Configurable<float> bRichPhotodetectorOtherModuleLength{"bRichPhotodetectorOtherModuleLength", 18.4, "barrel RICH photodetector other module length (cm)"};
  Configurable<float> bRichRadiatorInnerRadius{"bRichRadiatorInnerRadius", 85., "barrel RICH radiator inner radius (cm)"};
  Configurable<float> bRichPhotodetectorOuterRadius{"bRichPhotodetectorOuterRadius", 112., "barrel RICH photodetector outer radius (cm)"};
  Configurable<float> bRichRadiatorThickness{"bRichRadiatorThickness", 2., "barrel RICH radiator thickness (cm)"};
  Configurable<float> bRichRefractiveIndex{"bRichRefractiveIndex", 1.03, "barrel RICH refractive index"};
  Configurable<bool> bRichFlagAbsorbingWalls{"bRichFlagAbsorbingWalls", false, "barrel RICH flag absorbing walls between sectors"};
  Configurable<int> nStepsLIntegrator{"nStepsLIntegrator", 200, "number of steps in length integrator"};
  Configurable<bool> doQAplots{"doQAplots", true, "do basic velocity plot qa"};
  Configurable<int> nBinsThetaRing{"nBinsThetaRing", 3000, "number of bins in theta ring"};
  Configurable<int> nBinsP{"nBinsP", 400, "number of bins in momentum"};
  Configurable<int> nBinsNsigmaCorrectSpecies{"nBinsNsigmaCorrectSpecies", 200, "number of bins in Nsigma plot (correct speies)"};
  Configurable<int> nBinsNsigmaWrongSpecies{"nBinsNsigmaWrongSpecies", 200, "number of bins in Nsigma plot (wrong species)"};
  Configurable<int> nBinsAngularRes{"nBinsAngularRes", 400, "number of bins plots angular resolution"};
  Configurable<int> nBinsRelativeEtaPt{"nBinsRelativeEtaPt", 400, "number of bins plots pt and eta relative errors"};
  Configurable<int> nBinsEta{"nBinsEta", 400, "number of bins plot relative eta error"};
  Configurable<bool> flagIncludeTrackAngularRes{"flagIncludeTrackAngularRes", true, "flag to include or exclude track time resolution"};
  Configurable<float> multiplicityEtaRange{"multiplicityEtaRange", 0.800000012, "eta range to compute the multiplicity"};
  Configurable<bool> flagRICHLoadDelphesLUTs{"flagRICHLoadDelphesLUTs", false, "flag to load Delphes LUTs for tracking correction (use recoTrack parameters if false)"};

  Configurable<std::string> lutEl{"lutEl", "lutCovm.el.dat", "LUT for electrons"};
  Configurable<std::string> lutMu{"lutMu", "lutCovm.mu.dat", "LUT for muons"};
  Configurable<std::string> lutPi{"lutPi", "lutCovm.pi.dat", "LUT for pions"};
  Configurable<std::string> lutKa{"lutKa", "lutCovm.ka.dat", "LUT for kaons"};
  Configurable<std::string> lutPr{"lutPr", "lutCovm.pr.dat", "LUT for protons"};

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // Track smearer (here used to get relative pt and eta uncertainties)
  o2::delphes::DelphesO2TrackSmearer mSmearer;

  // needed: random number generator for smearing
  TRandom3 pRandomNumberGenerator;

  // for handling basic QA histograms if requested
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  /// Flag unphysical and unavailable values
  float error_value = -1000;

  // Variables projective/hybrid layout
  int mNumberSectors = bRichNumberOfSectors;                            // 21;
  float mTileLength = bRichPhotodetectorOtherModuleLength;              // 18.4; // [cm]
  float mTileLengthCentral = bRichPhotodetectorCentralModuleHalfLength; // 18.4 / 2.0; // [cm]
  float mProjectiveLengthInner = mTileLengthCentral;
  float mRadiusProjIn = bRichRadiatorInnerRadius;       // 85.;   // [cm]
  float mRadiusProjOut = bRichPhotodetectorOuterRadius; // 112.; // [cm]
  std::vector<TVector3> det_centers;
  std::vector<TVector3> rad_centers;
  std::vector<float> angle_centers;
  std::vector<float> theta_min;
  std::vector<float> theta_max;

  // Update projective geometry
  void updateProjectiveParameters()
  {
    const int number_of_sectors_in_z = mNumberSectors;
    det_centers.resize(number_of_sectors_in_z);
    rad_centers.resize(number_of_sectors_in_z);
    angle_centers.resize(number_of_sectors_in_z);
    theta_min.resize(number_of_sectors_in_z);
    theta_max.resize(number_of_sectors_in_z);
    float square_size_barrel_cylinder = 2.0 * mTileLengthCentral;
    float square_size_z = mTileLength;
    float R_min = mRadiusProjIn;
    float R_max = mRadiusProjOut;
    std::vector<float> theta_bi;
    std::vector<float> R0_tilt;
    std::vector<float> z0_tilt;
    std::vector<float> T_r_plus_g;
    std::vector<float> l_aerogel_z;
    std::vector<float> l_detector_z;
    theta_bi.resize(number_of_sectors_in_z);
    R0_tilt.resize(number_of_sectors_in_z);
    z0_tilt.resize(number_of_sectors_in_z);
    T_r_plus_g.resize(number_of_sectors_in_z);
    l_aerogel_z.resize(number_of_sectors_in_z);
    l_detector_z.resize(number_of_sectors_in_z);

    // Central sector
    int i_central_mirror = static_cast<int>((number_of_sectors_in_z) / 2.0);
    float m_val = std::tan(0.0);
    theta_bi[i_central_mirror] = std::atan(m_val);
    R0_tilt[i_central_mirror] = R_max;
    z0_tilt[i_central_mirror] = R0_tilt[i_central_mirror] * std::tan(theta_bi[i_central_mirror]);
    l_detector_z[i_central_mirror] = square_size_barrel_cylinder;
    l_aerogel_z[i_central_mirror] = std::sqrt(1.0 + m_val * m_val) * R_min * square_size_barrel_cylinder / (std::sqrt(1.0 + m_val * m_val) * R_max - m_val * square_size_barrel_cylinder);
    T_r_plus_g[i_central_mirror] = R_max - R_min;
    float t = std::tan(atan(m_val) + std::atan(square_size_barrel_cylinder / (2.0 * R_max * std::sqrt(1.0 + m_val * m_val) - square_size_barrel_cylinder * m_val)));
    theta_max[i_central_mirror] = M_PI / 2.0 - std::atan(t);
    theta_min[i_central_mirror] = M_PI / 2.0 + std::atan(t);
    mProjectiveLengthInner = R_min * t;
    for (int i = static_cast<int>((number_of_sectors_in_z) / 2.0) + 1; i < number_of_sectors_in_z; i++) {
      float par_a = t;
      float par_b = 2.0 * R_max / square_size_z;
      m_val = (std::sqrt(par_a * par_a * par_b * par_b + par_b * par_b - 1.0) + par_a * par_b * par_b) / (par_b * par_b - 1.0);
      theta_min[i] = M_PI / 2.0 - std::atan(t);
      theta_max[2 * i_central_mirror - i] = M_PI / 2.0 + std::atan(t);
      t = std::tan(std::atan(m_val) + std::atan(square_size_z / (2.0 * R_max * std::sqrt(1.0 + m_val * m_val) - square_size_z * m_val)));
      theta_max[i] = M_PI / 2.0 - std::atan(t);
      theta_min[2 * i_central_mirror - i] = M_PI / 2.0 + std::atan(t);
      // Forward sectors
      theta_bi[i] = std::atan(m_val);
      R0_tilt[i] = R_max - square_size_z / 2.0 * std::sin(std::atan(m_val));
      z0_tilt[i] = R0_tilt[i] * std::tan(theta_bi[i]);
      l_detector_z[i] = square_size_z;
      l_aerogel_z[i] = std::sqrt(1.0 + m_val * m_val) * R_min * square_size_z / (std::sqrt(1.0 + m_val * m_val) * R_max - m_val * square_size_z);
      T_r_plus_g[i] = std::sqrt(1.0 + m_val * m_val) * (R_max - R_min) - m_val / 2.0 * (square_size_z + l_aerogel_z[i]);
      // Backword sectors
      theta_bi[2 * i_central_mirror - i] = -std::atan(m_val);
      R0_tilt[2 * i_central_mirror - i] = R_max - square_size_z / 2.0 * std::sin(std::atan(m_val));
      z0_tilt[2 * i_central_mirror - i] = -R0_tilt[i] * std::tan(theta_bi[i]);
      l_detector_z[2 * i_central_mirror - i] = square_size_z;
      l_aerogel_z[2 * i_central_mirror - i] = std::sqrt(1.0 + m_val * m_val) * R_min * square_size_z / (std::sqrt(1.0 + m_val * m_val) * R_max - m_val * square_size_z);
      T_r_plus_g[2 * i_central_mirror - i] = std::sqrt(1.0 + m_val * m_val) * (R_max - R_min) - m_val / 2.0 * (square_size_z + l_aerogel_z[i]);
      mProjectiveLengthInner = R_min * t; // <-- At the end of the loop this will be the maximum Z
    }
    // Coordinate radiali layer considerati
    std::vector<float> R0_detector;
    std::vector<float> R0_aerogel;
    R0_detector.resize(number_of_sectors_in_z);
    R0_aerogel.resize(number_of_sectors_in_z);
    for (int i = 0; i < number_of_sectors_in_z; i++) {
      R0_detector[i] = R0_tilt[i];
      det_centers[i].SetXYZ(R0_detector[i], 0, R0_detector[i] * std::tan(theta_bi[i]));
      R0_aerogel[i] = R0_tilt[i] - (T_r_plus_g[i]) * std::cos(theta_bi[i]);
      rad_centers[i].SetXYZ(R0_aerogel[i], 0, R0_aerogel[i] * std::tan(theta_bi[i]));
      angle_centers[i] = theta_bi[i];
    }
  }

  void init(o2::framework::InitContext& initContext)
  {
    pRandomNumberGenerator.SetSeed(0); // fully randomize

    // Load LUT for pt and eta smearing
    if (flagIncludeTrackAngularRes && flagRICHLoadDelphesLUTs) {
      std::map<int, const char*> mapPdgLut;
      const char* lutElChar = lutEl->c_str();
      const char* lutMuChar = lutMu->c_str();
      const char* lutPiChar = lutPi->c_str();
      const char* lutKaChar = lutKa->c_str();
      const char* lutPrChar = lutPr->c_str();

      LOGF(info, "Will load electron lut file ..: %s for RICH PID", lutElChar);
      LOGF(info, "Will load muon lut file ......: %s for RICH PID", lutMuChar);
      LOGF(info, "Will load pion lut file ......: %s for RICH PID", lutPiChar);
      LOGF(info, "Will load kaon lut file ......: %s for RICH PID", lutKaChar);
      LOGF(info, "Will load proton lut file ....: %s for RICH PID", lutPrChar);

      mapPdgLut.insert(std::make_pair(11, lutElChar));
      mapPdgLut.insert(std::make_pair(13, lutMuChar));
      mapPdgLut.insert(std::make_pair(211, lutPiChar));
      mapPdgLut.insert(std::make_pair(321, lutKaChar));
      mapPdgLut.insert(std::make_pair(2212, lutPrChar));

      for (auto e : mapPdgLut) {
        if (!mSmearer.loadTable(e.first, e.second)) {
          LOG(fatal) << "Having issue with loading the LUT " << e.first << " " << e.second;
        }
      }
    }

    if (doQAplots) {
      const AxisSpec axisMomentum{static_cast<int>(nBinsP), 0.0f, +20.0f, "#it{p} (GeV/#it{c})"};
      const AxisSpec axisAngle{static_cast<int>(nBinsThetaRing), 0.0f, +0.30f, "Measured Cherenkov angle (rad)"};
      histos.add("h2dAngleVsMomentumBarrelRICH", "h2dAngleVsMomentumBarrelRICH", kTH2F, {axisMomentum, axisAngle});

      std::string particle_names1[5] = {"#it{e}", "#it{#mu}", "#it{#pi}", "#it{K}", "#it{p}"};
      std::string particle_names2[5] = {"Elec", "Muon", "Pion", "Kaon", "Prot"};
      for (int i_true = 0; i_true < 5; i_true++) {
        std::string name_title_barrel_track_res = "h2dBarrelAngularResTrack" + particle_names2[i_true] + "VsP";
        std::string name_title_barrel_total_res = "h2dBarrelAngularResTotal" + particle_names2[i_true] + "VsP";
        const AxisSpec axisTrackAngularRes{static_cast<int>(nBinsAngularRes), 0.0f, +5.0f, "Track angular resolution - " + particle_names1[i_true] + " (mrad)"};
        const AxisSpec axisTotalAngularRes{static_cast<int>(nBinsAngularRes), 0.0f, +5.0f, "Total angular resolution - " + particle_names1[i_true] + " (mrad)"};
        histos.add(name_title_barrel_track_res.c_str(), name_title_barrel_track_res.c_str(), kTH2F, {axisMomentum, axisTrackAngularRes});
        histos.add(name_title_barrel_total_res.c_str(), name_title_barrel_total_res.c_str(), kTH2F, {axisMomentum, axisTotalAngularRes});
      }
      for (int i_true = 0; i_true < 5; i_true++) {
        for (int i_hyp = 0; i_hyp < 5; i_hyp++) {
          std::string name_title = "h2dBarrelNsigmaTrue" + particle_names2[i_true] + "Vs" + particle_names2[i_hyp] + "Hypothesis";
          if (i_true == i_hyp) {
            const AxisSpec axisNsigmaCorrect{static_cast<int>(nBinsNsigmaCorrectSpecies), -10.0f, +10.0f, "N#sigma - True " + particle_names1[i_true] + " vs " + particle_names1[i_hyp] + " hypothesis"};
            histos.add(name_title.c_str(), name_title.c_str(), kTH2F, {axisMomentum, axisNsigmaCorrect});
          } else {
            const AxisSpec axisNsigmaWrong{static_cast<int>(nBinsNsigmaWrongSpecies), -10.0f, +10.0f, "N#sigma -  True " + particle_names1[i_true] + " vs " + particle_names1[i_hyp] + " hypothesis"};
            histos.add(name_title.c_str(), name_title.c_str(), kTH2F, {axisMomentum, axisNsigmaWrong});
          }
        }
      }
    }

    // Update projective parameters
    updateProjectiveParameters();
  }

  /// Function to convert a McParticle into a perfect Track
  /// \param particle the particle to convert (mcParticle)
  /// \param o2track the address of the resulting TrackParCov
  template <typename McParticleType>
  o2::track::TrackParCov convertMCParticleToO2Track(McParticleType& particle)
  {
    // FIXME: this is a fundamentally important piece of code.
    // It could be placed in a utility file instead of here.
    auto pdgInfo = pdg->GetParticle(particle.pdgCode());
    int charge = 0;
    if (pdgInfo != nullptr) {
      charge = pdgInfo->Charge() / 3;
    }
    std::array<float, 5> params;
    std::array<float, 15> covm = {0.};
    float s, c, x;
    o2::math_utils::sincos(particle.phi(), s, c);
    o2::math_utils::rotateZInv(particle.vx(), particle.vy(), x, params[0], s, c);
    params[1] = particle.vz();
    params[2] = 0.; // since alpha = phi
    auto theta = 2. * std::atan(std::exp(-particle.eta()));
    params[3] = 1. / std::tan(theta);
    params[4] = charge / particle.pt();

    // Return TrackParCov
    return o2::track::TrackParCov(x, particle.phi(), params, covm);
  }

  /// check if particle reaches radiator
  /// \param track the input track
  /// \param radius the radius of the layer you're calculating the length to
  /// \param magneticField the magnetic field to use when propagating
  bool checkMagfieldLimit(o2::track::TrackParCov track, float radius, float magneticField)
  {
    o2::math_utils::CircleXYf_t trcCircle;
    float sna, csa;
    track.getCircleParams(magneticField, trcCircle, sna, csa);

    // distance between circle centers (one circle is at origin -> easy)
    float centerDistance = std::hypot(trcCircle.xC, trcCircle.yC);

    // condition of circles touching - if not satisfied returned value if false
    if (centerDistance < trcCircle.rC + radius && centerDistance > std::fabs(trcCircle.rC - radius)) {
      return true;
    } else {
      return false;
    }
  }

  /// returns radiator radius in cm at considered eta (accounting for sector inclination)
  /// \param eta the pseudorapidity of the tarck (assuming primary vertex at origin)
  float radiusRipple(float eta)
  {
    float polar = 2.0 * std::atan(std::exp(-eta));
    float i_sector = 0;
    bool flag_sector = false;
    for (int j_sec = 0; j_sec < mNumberSectors; j_sec++) {
      if (polar > theta_max[j_sec] && polar < theta_min[j_sec]) {
        flag_sector = true;
        i_sector = j_sec;
      }
    }
    if (flag_sector) {
      float R_sec_rich = rad_centers[i_sector].X();
      float z_sec_rich = rad_centers[i_sector].Z();
      return (std::pow(R_sec_rich, 2) + std::pow(z_sec_rich, 2)) / (R_sec_rich + z_sec_rich / std::tan(polar));
    } else {
      return error_value;
    }
  }

  /// returns Cherenkov angle in rad (above threshold) or bad flag (below threshold)
  /// \param momentum the momentum of the tarck
  /// \param mass the mass of the particle
  float CherenkovAngle(float momentum, float mass)
  {
    // Check if particle is above the threshold
    if (momentum > mass / std::sqrt(bRichRefractiveIndex * bRichRefractiveIndex - 1.0)) {
      // Calculate angle
      float angle = std::acos(std::sqrt(momentum * momentum + mass * mass) / (momentum * bRichRefractiveIndex));

      // Mean number of detected photons (SiPM PDE is accountd in multiplicative factor!!!)
      float meanNumberofDetectedPhotons = 230. * std::sin(angle) * std::sin(angle) * bRichRadiatorThickness;

      // Require at least 3 photons on average for real angle reconstruction
      if (meanNumberofDetectedPhotons > 3.) {
        return angle;
      } else {
        return error_value;
      }
    } else {
      return error_value;
    }
  }

  /// returns linear interpolation
  /// \param x the eta we want the resolution for
  /// \param x0 the closest smaller available eta
  /// \param x1 the closest larger available eta
  /// \param y0 the resolution corresponding to x0
  /// \param y1 the resolution corresponding to x1
  float interpolate(double x, double x0, double x1, double y0, double y1)
  {
    float y = y0 + ((y1 - y0) / (x1 - x0)) * (x - x0);
    return y;
  }

  /// returns angular resolution for considered track eta
  /// \param eta the pseudorapidity of the tarck (assuming primary vertex at origin)
  float AngularResolution(float eta)
  {
    // Vectors for sampling (USE ANALYTICAL EXTRAPOLATION FOR BETTER RESULTS)
    float eta_sampling[] = {-2.000000, -1.909740, -1.731184, -1.552999, -1.375325, -1.198342, -1.022276, -0.847390, -0.673976, -0.502324, -0.332683, -0.165221, 0.000000, 0.165221, 0.332683, 0.502324, 0.673976, 0.847390, 1.022276, 1.198342, 1.375325, 1.552999, 1.731184, 1.909740, 2.000000};
    float res_ring_sampling_with_abs_walls[] = {0.0009165, 0.000977, 0.001098, 0.001198, 0.001301, 0.001370, 0.001465, 0.001492, 0.001498, 0.001480, 0.001406, 0.001315, 0.001241, 0.001325, 0.001424, 0.001474, 0.001480, 0.001487, 0.001484, 0.001404, 0.001273, 0.001197, 0.001062, 0.000965, 0.0009165};
    float res_ring_sampling_without_abs_walls[] = {0.0009165, 0.000977, 0.001095, 0.001198, 0.001300, 0.001369, 0.001468, 0.001523, 0.001501, 0.001426, 0.001299, 0.001167, 0.001092, 0.001179, 0.001308, 0.001407, 0.001491, 0.001508, 0.001488, 0.001404, 0.001273, 0.001196, 0.001061, 0.000965, 0.0009165};
    int size_res_vector = sizeof(eta_sampling) / sizeof(eta_sampling[0]);
    // Use binary search to find the lower and upper indices
    int lowerIndex = std::lower_bound(eta_sampling, eta_sampling + size_res_vector, eta) - eta_sampling - 1;
    int upperIndex = lowerIndex + 1;
    if (lowerIndex >= 0 && upperIndex < size_res_vector) {
      // Resolution interpolation
      if (bRichFlagAbsorbingWalls) {
        float interpolatedResRing = interpolate(eta, eta_sampling[lowerIndex], eta_sampling[upperIndex], res_ring_sampling_with_abs_walls[lowerIndex], res_ring_sampling_with_abs_walls[upperIndex]);
        // std::cout << "Interpolated y value: " << interpolatedY << std::endl;
        return interpolatedResRing;
      } else {
        float interpolatedResRing = interpolate(eta, eta_sampling[lowerIndex], eta_sampling[upperIndex], res_ring_sampling_without_abs_walls[lowerIndex], res_ring_sampling_without_abs_walls[upperIndex]);
        // std::cout << "Interpolated y value: " << interpolatedY << std::endl;
        return interpolatedResRing;
      }
    } else {
      // std::cout << "Unable to interpolate. Target x value is outside the range of available data." << std::endl;
      return error_value;
    }
  }

  /// returns track angular resolution
  /// \param pt the transverse momentum of the tarck
  /// \param eta the pseudorapidity of the tarck
  /// \param track_pt_resolution the absolute resolution on pt
  /// \param track_pt_resolution the absolute resolution on eta
  /// \param mass the mass of the particle
  /// \param refractive_index the refractive index of the radiator
  double calculate_track_time_resolution_advanced(float pt, float eta, float track_pt_resolution, float track_eta_resolution, float mass, float refractive_index)
  {
    // Compute tracking contribution to timing using the error propagation formula
    // Uses light speed in m/ps, magnetic field in T (*0.1 for conversion kGauss -> T)
    double a0 = mass * mass;
    double a1 = refractive_index;
    double dtheta_on_dpt = a0 / (pt * std::sqrt(a0 + std::pow(pt * std::cosh(eta), 2)) * std::sqrt(std::pow(pt * std::cosh(eta), 2) * (std::pow(a1, 2) - 1.0) - a0));
    double dtheta_on_deta = (a0 * std::tanh(eta)) / (std::sqrt(a0 + std::pow(pt * std::cosh(eta), 2)) * std::sqrt(std::pow(pt * std::cosh(eta), 2) * (std::pow(a1, 2) - 1.0) - a0));
    double track_angular_resolution = std::hypot(std::fabs(dtheta_on_dpt) * track_pt_resolution, std::fabs(dtheta_on_deta) * track_eta_resolution);
    return track_angular_resolution;
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels>::iterator const& collision, soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels> const& tracks, aod::McParticles const&, aod::McCollisions const&)
  {

    o2::dataformats::VertexBase pvVtx({collision.posX(), collision.posY(), collision.posZ()},
                                      {collision.covXX(), collision.covXY(), collision.covYY(), collision.covXZ(), collision.covYZ(), collision.covZZ()});

    std::array<float, 6> mcPvCov = {0.};
    o2::dataformats::VertexBase mcPvVtx({0.0f, 0.0f, 0.0f}, mcPvCov);
    if (collision.has_mcCollision()) {
      auto mcCollision = collision.mcCollision();
      mcPvVtx.setX(mcCollision.posX());
      mcPvVtx.setY(mcCollision.posY());
      mcPvVtx.setZ(mcCollision.posZ());
    } // else remains untreated for now

    // First we compute the number of charged particles in the event
    float dNdEta = 0.f;
    if (flagRICHLoadDelphesLUTs) {
      for (const auto& track : tracks) {
        if (!track.has_mcParticle())
          continue;
        auto mcParticle = track.mcParticle();
        if (std::abs(mcParticle.eta()) > multiplicityEtaRange) {
          continue;
        }
        if (mcParticle.has_daughters()) {
          continue;
        }
        const auto& pdgInfo = pdg->GetParticle(mcParticle.pdgCode());
        if (!pdgInfo) {
          // LOG(warning) << "PDG code " << mcParticle.pdgCode() << " not found in the database";
          continue;
        }
        if (pdgInfo->Charge() == 0) {
          continue;
        }
        dNdEta += 1.f;
      }
    }

    for (const auto& track : tracks) {
      // first step: find precise arrival time (if any)
      // --- convert track into perfect track
      if (!track.has_mcParticle()) // should always be OK but check please
        continue;

      auto mcParticle = track.mcParticle();
      o2::track::TrackParCov o2track = convertMCParticleToO2Track(mcParticle);

      // float xPv = error_value;
      if (o2track.propagateToDCA(mcPvVtx, dBz)) {
        // xPv = o2track.getX();
      }

      // get particle to calculate Cherenkov angle and resolution
      auto pdgInfo = pdg->GetParticle(mcParticle.pdgCode());
      if (pdgInfo == nullptr) {
        continue;
      }
      float expectedAngleBarrelRich = CherenkovAngle(o2track.getP(), pdgInfo->Mass());
      float barrelRICHAngularResolution = AngularResolution(o2track.getEta());
      float projectiveRadiatorRadius = radiusRipple(o2track.getEta());
      bool flagReachesRadiator = false;
      if (projectiveRadiatorRadius > error_value + 1.) {
        flagReachesRadiator = checkMagfieldLimit(o2track, projectiveRadiatorRadius, dBz);
      }
      /// DISCLAIMER: Exact extrapolation of angular resolution would require track propagation
      ///             to the RICH radiator (accounting sector inclination) in terms of (R,z).
      ///             The extrapolation with Eta is correct only if the primary vertex is at origin.
      ///             Discrepancies may be negligible, but would be more rigorous if propagation tool is available

      // Smear with expected resolutions
      float measuredAngleBarrelRich = pRandomNumberGenerator.Gaus(expectedAngleBarrelRich, barrelRICHAngularResolution);

      // Now we calculate the expected arrival time following certain mass hypotheses
      // and the (imperfect!) reconstructed track parametrizations
      auto recoTrack = getTrackParCov(track);
      if (recoTrack.propagateToDCA(pvVtx, dBz)) {
        // xPv = recoTrack.getX();
      }

      // Straight to Nsigma
      float deltaThetaBarrelRich[5], nSigmaBarrelRich[5];
      int lpdg_array[5] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton};
      float masses[5];

      for (int ii = 0; ii < 5; ii++) {
        nSigmaBarrelRich[ii] = error_value;

        auto pdgInfoThis = pdg->GetParticle(lpdg_array[ii]);
        masses[ii] = pdgInfoThis->Mass();
        float hypothesisAngleBarrelRich = CherenkovAngle(recoTrack.getP(), masses[ii]);

        // Evaluate total sigma (layer + tracking resolution)
        float barrelTotalAngularReso = barrelRICHAngularResolution;
        if (flagIncludeTrackAngularRes) {
          double pt_resolution = std::pow(recoTrack.getP() / std::cosh(recoTrack.getEta()), 2) * std::sqrt(recoTrack.getSigma1Pt2());
          double eta_resolution = std::fabs(std::sin(2.0 * std::atan(std::exp(-recoTrack.getEta())))) * std::sqrt(recoTrack.getSigmaTgl2());
          if (flagRICHLoadDelphesLUTs) {
            pt_resolution = mSmearer.getAbsPtRes(pdgInfoThis->PdgCode(), dNdEta, recoTrack.getEta(), recoTrack.getP() / std::cosh(recoTrack.getEta()));
            eta_resolution = mSmearer.getAbsEtaRes(pdgInfoThis->PdgCode(), dNdEta, recoTrack.getEta(), recoTrack.getP() / std::cosh(recoTrack.getEta()));
          }
          // cout << endl <<  "Pt resolution: " << pt_resolution << ", Eta resolution: " << eta_resolution << endl << endl;
          float barrelTrackAngularReso = calculate_track_time_resolution_advanced(recoTrack.getP() / std::cosh(recoTrack.getEta()), recoTrack.getEta(), pt_resolution, eta_resolution, masses[ii], bRichRefractiveIndex);
          barrelTotalAngularReso = std::hypot(barrelRICHAngularResolution, barrelTrackAngularReso);
          if (doQAplots && hypothesisAngleBarrelRich > error_value + 1. && measuredAngleBarrelRich > error_value + 1. && barrelRICHAngularResolution > error_value + 1. && flagReachesRadiator) {
            float momentum = recoTrack.getP();
            // float pseudorapidity = recoTrack.getEta();
            // float transverse_momentum = momentum / std::cosh(pseudorapidity);
            if (ii == 0 && std::fabs(mcParticle.pdgCode()) == pdg->GetParticle(lpdg_array[0])->PdgCode()) {
              histos.fill(HIST("h2dBarrelAngularResTrackElecVsP"), momentum, 1000.0 * barrelTrackAngularReso);
              histos.fill(HIST("h2dBarrelAngularResTotalElecVsP"), momentum, 1000.0 * barrelTotalAngularReso);
            }
            if (ii == 1 && std::fabs(mcParticle.pdgCode()) == pdg->GetParticle(lpdg_array[1])->PdgCode()) {
              histos.fill(HIST("h2dBarrelAngularResTrackMuonVsP"), momentum, 1000.0 * barrelTrackAngularReso);
              histos.fill(HIST("h2dBarrelAngularResTotalMuonVsP"), momentum, 1000.0 * barrelTotalAngularReso);
            }
            if (ii == 2 && std::fabs(mcParticle.pdgCode()) == pdg->GetParticle(lpdg_array[2])->PdgCode()) {
              histos.fill(HIST("h2dBarrelAngularResTrackPionVsP"), momentum, 1000.0 * barrelTrackAngularReso);
              histos.fill(HIST("h2dBarrelAngularResTotalPionVsP"), momentum, 1000.0 * barrelTotalAngularReso);
            }
            if (ii == 3 && std::fabs(mcParticle.pdgCode()) == pdg->GetParticle(lpdg_array[3])->PdgCode()) {
              histos.fill(HIST("h2dBarrelAngularResTrackKaonVsP"), momentum, 1000.0 * barrelTrackAngularReso);
              histos.fill(HIST("h2dBarrelAngularResTotalKaonVsP"), momentum, 1000.0 * barrelTotalAngularReso);
            }
            if (ii == 4 && std::fabs(mcParticle.pdgCode()) == pdg->GetParticle(lpdg_array[4])->PdgCode()) {
              histos.fill(HIST("h2dBarrelAngularResTrackProtVsP"), momentum, 1000.0 * barrelTrackAngularReso);
              histos.fill(HIST("h2dBarrelAngularResTotalProtVsP"), momentum, 1000.0 * barrelTotalAngularReso);
            }
          }
        }

        /// DISCLAIMER: here tracking is accounted only for momentum value, but not for track parameters at impact point on the
        ///             RICH radiator, since exact resolution would require photon generation and transport to photodetector.
        ///             Effects are expected to be negligible (a few tenths of a milliradian) but further studies are required !
        if (hypothesisAngleBarrelRich > error_value + 1. && measuredAngleBarrelRich > error_value + 1. && barrelRICHAngularResolution > error_value + 1. && flagReachesRadiator) {
          deltaThetaBarrelRich[ii] = hypothesisAngleBarrelRich - measuredAngleBarrelRich;
          nSigmaBarrelRich[ii] = deltaThetaBarrelRich[ii] / barrelTotalAngularReso;
        }
      }

      // Fill histograms
      if (doQAplots) {
        float momentum = recoTrack.getP();
        float barrelRichTheta = measuredAngleBarrelRich;

        if (barrelRichTheta > error_value + 1. && barrelRICHAngularResolution > error_value + 1. && flagReachesRadiator) {
          histos.fill(HIST("h2dAngleVsMomentumBarrelRICH"), momentum, barrelRichTheta);

          if (std::fabs(mcParticle.pdgCode()) == pdg->GetParticle(lpdg_array[0])->PdgCode()) {
            histos.fill(HIST("h2dBarrelNsigmaTrueElecVsElecHypothesis"), momentum, nSigmaBarrelRich[0]);
            histos.fill(HIST("h2dBarrelNsigmaTrueElecVsMuonHypothesis"), momentum, nSigmaBarrelRich[1]);
            histos.fill(HIST("h2dBarrelNsigmaTrueElecVsPionHypothesis"), momentum, nSigmaBarrelRich[2]);
            histos.fill(HIST("h2dBarrelNsigmaTrueElecVsKaonHypothesis"), momentum, nSigmaBarrelRich[3]);
            histos.fill(HIST("h2dBarrelNsigmaTrueElecVsProtHypothesis"), momentum, nSigmaBarrelRich[4]);
          }
          if (std::fabs(mcParticle.pdgCode()) == pdg->GetParticle(lpdg_array[1])->PdgCode()) {
            histos.fill(HIST("h2dBarrelNsigmaTrueMuonVsElecHypothesis"), momentum, nSigmaBarrelRich[0]);
            histos.fill(HIST("h2dBarrelNsigmaTrueMuonVsMuonHypothesis"), momentum, nSigmaBarrelRich[1]);
            histos.fill(HIST("h2dBarrelNsigmaTrueMuonVsPionHypothesis"), momentum, nSigmaBarrelRich[2]);
            histos.fill(HIST("h2dBarrelNsigmaTrueMuonVsKaonHypothesis"), momentum, nSigmaBarrelRich[3]);
            histos.fill(HIST("h2dBarrelNsigmaTrueMuonVsProtHypothesis"), momentum, nSigmaBarrelRich[4]);
          }
          if (std::fabs(mcParticle.pdgCode()) == pdg->GetParticle(lpdg_array[2])->PdgCode()) {
            histos.fill(HIST("h2dBarrelNsigmaTruePionVsElecHypothesis"), momentum, nSigmaBarrelRich[0]);
            histos.fill(HIST("h2dBarrelNsigmaTruePionVsMuonHypothesis"), momentum, nSigmaBarrelRich[1]);
            histos.fill(HIST("h2dBarrelNsigmaTruePionVsPionHypothesis"), momentum, nSigmaBarrelRich[2]);
            histos.fill(HIST("h2dBarrelNsigmaTruePionVsKaonHypothesis"), momentum, nSigmaBarrelRich[3]);
            histos.fill(HIST("h2dBarrelNsigmaTruePionVsProtHypothesis"), momentum, nSigmaBarrelRich[4]);
          }
          if (std::fabs(mcParticle.pdgCode()) == pdg->GetParticle(lpdg_array[3])->PdgCode()) {
            histos.fill(HIST("h2dBarrelNsigmaTrueKaonVsElecHypothesis"), momentum, nSigmaBarrelRich[0]);
            histos.fill(HIST("h2dBarrelNsigmaTrueKaonVsMuonHypothesis"), momentum, nSigmaBarrelRich[1]);
            histos.fill(HIST("h2dBarrelNsigmaTrueKaonVsPionHypothesis"), momentum, nSigmaBarrelRich[2]);
            histos.fill(HIST("h2dBarrelNsigmaTrueKaonVsKaonHypothesis"), momentum, nSigmaBarrelRich[3]);
            histos.fill(HIST("h2dBarrelNsigmaTrueKaonVsProtHypothesis"), momentum, nSigmaBarrelRich[4]);
          }
          if (std::fabs(mcParticle.pdgCode()) == pdg->GetParticle(lpdg_array[4])->PdgCode()) {
            histos.fill(HIST("h2dBarrelNsigmaTrueProtVsElecHypothesis"), momentum, nSigmaBarrelRich[0]);
            histos.fill(HIST("h2dBarrelNsigmaTrueProtVsMuonHypothesis"), momentum, nSigmaBarrelRich[1]);
            histos.fill(HIST("h2dBarrelNsigmaTrueProtVsPionHypothesis"), momentum, nSigmaBarrelRich[2]);
            histos.fill(HIST("h2dBarrelNsigmaTrueProtVsKaonHypothesis"), momentum, nSigmaBarrelRich[3]);
            histos.fill(HIST("h2dBarrelNsigmaTrueProtVsProtHypothesis"), momentum, nSigmaBarrelRich[4]);
          }
        }
      }

      // Sigmas have been fully calculated. Please populate the NSigma helper table (once per track)
      upgradeRich(nSigmaBarrelRich[0], nSigmaBarrelRich[1], nSigmaBarrelRich[2], nSigmaBarrelRich[3], nSigmaBarrelRich[4]);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyRichPid>(cfgc)};
}
