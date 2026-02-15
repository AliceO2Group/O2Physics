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
/// \file onTheFlyRichPid.cxx
///
/// \brief This task goes straight from a combination of track table and mcParticles
///        and a projective bRICH configuration to a table of TOF NSigmas for the particles
///        being analysed. It currently contemplates 5 particle types:
///        electrons, pions, kaons, protons and muons
///
///        More particles could be added but would have to be added to the LUT
///        being used in the onTheFly tracker task.
///
/// \warning Geometry parameters are configurable, but resolution values should be adapted.
///          Since angular resolution depends on the specific geometric details, it is better to
///          calculate it from full simulation and add new input. Alternatively, an analytical
///          expression can be provided as a function of the main parameters.
///          Latest version: analytical parametrization of angular resolution !!!
///
/// \author David Dobrigkeit Chinellato, UNICAMP
/// \author Nicola Nicassio, University and INFN Bari
/// \since  May 22, 2024
///

#include "ALICE3/Core/DelphesO2TrackSmearer.h"
#include "ALICE3/Core/FastTracker.h"
#include "ALICE3/Core/TrackUtilities.h"
#include "ALICE3/DataModel/OTFCollision.h"
#include "ALICE3/DataModel/OTFRICH.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <CCDB/BasicCCDBManager.h>
#include <CCDB/CcdbApi.h>
#include <CommonConstants/GeomConstants.h>
#include <CommonConstants/MathConstants.h>
#include <CommonConstants/PhysicsConstants.h>
#include <CommonUtils/NameConf.h>
#include <DataFormatsCalibration/MeanVertexObject.h>
#include <DataFormatsParameters/GRPMagField.h>
#include <DetectorsBase/GeometryManager.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/ASoAHelpers.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/RunningWorkflowInfo.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/DCA.h>
#include <ReconstructionDataFormats/HelixHelper.h>
#include <ReconstructionDataFormats/PID.h>

#include <TPDGCode.h>
#include <TRandom3.h>
#include <TString.h>
#include <TVector3.h>

#include <cmath>
#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::constants::math;

struct OnTheFlyRichPid {
  Produces<aod::UpgradeRich> upgradeRich;
  Produces<aod::UpgradeRichSignal> upgradeRichSignal;

  // necessary for particle charges
  Service<o2::framework::O2DatabasePDG> pdg;
  // Necessary for LUTs
  Service<o2::ccdb::BasicCCDBManager> ccdb;

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
  Configurable<float> gasRadiatorRindex{"gasRadiatorRindex", 1.0006f, "gas radiator refractive index"};
  Configurable<float> gasRichRadiatorThickness{"gasRichRadiatorThickness", 25.f, "gas radiator thickness (cm)"};
  Configurable<float> bRichRefractiveIndexSector0{"bRichRefractiveIndexSector0", 1.03, "barrel RICH refractive index central(s)"};                        // central(s)
  Configurable<float> bRichRefractiveIndexSector1{"bRichRefractiveIndexSector1", 1.03, "barrel RICH refractive index central(s)-1 and central(s)+1"};     // central(s)-1 and central(s)+1
  Configurable<float> bRichRefractiveIndexSector2{"bRichRefractiveIndexSector2", 1.03, "barrel RICH refractive index central(s)-2 and central(s)+2"};     // central(s)-2 and central(s)+2
  Configurable<float> bRichRefractiveIndexSector3{"bRichRefractiveIndexSector3", 1.03, "barrel RICH refractive index central(s)-3 and central(s)+3"};     // central(s)-3 and central(s)+3
  Configurable<float> bRichRefractiveIndexSector4{"bRichRefractiveIndexSector4", 1.03, "barrel RICH refractive index central(s)-4 and central(s)+4"};     // central(s)-4 and central(s)+4
  Configurable<float> bRichRefractiveIndexSector5{"bRichRefractiveIndexSector5", 1.03, "barrel RICH refractive index central(s)-5 and central(s)+5"};     // central(s)-5 and central(s)+5
  Configurable<float> bRichRefractiveIndexSector6{"bRichRefractiveIndexSector6", 1.03, "barrel RICH refractive index central(s)-6 and central(s)+6"};     // central(s)-6 and central(s)+6
  Configurable<float> bRichRefractiveIndexSector7{"bRichRefractiveIndexSector7", 1.03, "barrel RICH refractive index central(s)-7 and central(s)+7"};     // central(s)-7 and central(s)+7
  Configurable<float> bRichRefractiveIndexSector8{"bRichRefractiveIndexSector8", 1.03, "barrel RICH refractive index central(s)-8 and central(s)+8"};     // central(s)-8 and central(s)+8
  Configurable<float> bRichRefractiveIndexSector9{"bRichRefractiveIndexSector9", 1.03, "barrel RICH refractive index central(s)-9 and central(s)+9"};     // central(s)-9 and central(s)+9
  Configurable<float> bRichRefractiveIndexSector10{"bRichRefractiveIndexSector10", 1.03, "barrel RICH refractive index central(s)-10 and central(s)+10"}; // central(s)-10 and central(s)+10
  Configurable<float> bRichRefractiveIndexSector11{"bRichRefractiveIndexSector11", 1.03, "barrel RICH refractive index central(s)-11 and central(s)+11"}; // central(s)-11 and central(s)+11
  Configurable<float> bRichRefractiveIndexSector12{"bRichRefractiveIndexSector12", 1.03, "barrel RICH refractive index central(s)-12 and central(s)+12"}; // central(s)-12 and central(s)+12
  Configurable<float> bRichRefractiveIndexSector13{"bRichRefractiveIndexSector13", 1.03, "barrel RICH refractive index central(s)-13 and central(s)+13"}; // central(s)-13 and central(s)+13
  Configurable<float> bRichRefractiveIndexSector14{"bRichRefractiveIndexSector14", 1.03, "barrel RICH refractive index central(s)-14 and central(s)+14"}; // central(s)-14 and central(s)+14
  Configurable<float> bRichRefractiveIndexSector15{"bRichRefractiveIndexSector15", 1.03, "barrel RICH refractive index central(s)-15 and central(s)+15"}; // central(s)-15 and central(s)+15
  Configurable<float> bRichRefractiveIndexSector16{"bRichRefractiveIndexSector16", 1.03, "barrel RICH refractive index central(s)-16 and central(s)+16"}; // central(s)-16 and central(s)+16
  Configurable<float> bRichRefractiveIndexSector17{"bRichRefractiveIndexSector17", 1.03, "barrel RICH refractive index central(s)-17 and central(s)+17"}; // central(s)-17 and central(s)+17
  Configurable<float> bRichRefractiveIndexSector18{"bRichRefractiveIndexSector18", 1.03, "barrel RICH refractive index central(s)-18 and central(s)+18"}; // central(s)-18 and central(s)+18
  Configurable<float> bRichRefractiveIndexSector19{"bRichRefractiveIndexSector19", 1.03, "barrel RICH refractive index central(s)-19 and central(s)+19"}; // central(s)-19 and central(s)+19
  Configurable<float> bRichRefractiveIndexSector20{"bRichRefractiveIndexSector20", 1.03, "barrel RICH refractive index central(s)-20 and central(s)+20"}; // central(s)-20 and central(s)+20
  Configurable<float> bRICHPixelSize{"bRICHPixelSize", 0.1, "barrel RICH pixel size (cm)"};
  Configurable<float> bRichGapRefractiveIndex{"bRichGapRefractiveIndex", 1.000283, "barrel RICH gap refractive index"};

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // Track smearer array, one per geometry
  std::vector<std::unique_ptr<o2::delphes::DelphesO2TrackSmearer>> mSmearer;

  // needed: random number generator for smearing
  TRandom3 pRandomNumberGenerator;

  // for handling basic QA histograms if requested
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  /// Flag unphysical and unavailable values (must be negative)
  static constexpr float kErrorValue = -1000;

  // Variables projective/hybrid layout
  std::vector<TVector3> detCenters;
  std::vector<TVector3> radCenters;
  std::vector<float> angleCenters;
  std::vector<float> thetaMin;
  std::vector<float> thetaMax;
  std::vector<float> bRichRefractiveIndexSector;
  std::vector<float> aerogelRindex;
  std::vector<float> photodetrctorLength;
  std::vector<float> gapThickness;

  // Update projective geometry
  void updateProjectiveParameters()
  {
    const int numberOfSectorsInZ = bRichNumberOfSectors;
    detCenters.resize(numberOfSectorsInZ);
    radCenters.resize(numberOfSectorsInZ);
    angleCenters.resize(numberOfSectorsInZ);
    thetaMin.resize(numberOfSectorsInZ);
    thetaMax.resize(numberOfSectorsInZ);
    aerogelRindex.resize(numberOfSectorsInZ);
    photodetrctorLength.resize(numberOfSectorsInZ);
    gapThickness.resize(numberOfSectorsInZ);
    float squareSizeBarrelCylinder = 2.0 * bRichPhotodetectorCentralModuleHalfLength;
    float squareSizeZ = bRichPhotodetectorOtherModuleLength;
    float rMin = bRichRadiatorInnerRadius;
    float rMax = bRichPhotodetectorOuterRadius;
    std::vector<float> thetaBi;
    std::vector<float> r0Tilt;
    std::vector<float> z0Tilt;
    std::vector<float> tRPlusG;
    std::vector<float> lAerogelZ;
    std::vector<float> lDetectorZ;
    thetaBi.resize(numberOfSectorsInZ);
    r0Tilt.resize(numberOfSectorsInZ);
    z0Tilt.resize(numberOfSectorsInZ);
    tRPlusG.resize(numberOfSectorsInZ);
    lAerogelZ.resize(numberOfSectorsInZ);
    lDetectorZ.resize(numberOfSectorsInZ);

    // Odd number of sectors
    static constexpr int kTwo = 2;
    if (numberOfSectorsInZ % kTwo != 0) {
      int iCentralMirror = static_cast<int>((numberOfSectorsInZ) / 2.0);
      float mVal = std::tan(0.0);
      thetaBi[iCentralMirror] = std::atan(mVal);
      r0Tilt[iCentralMirror] = rMax;
      z0Tilt[iCentralMirror] = r0Tilt[iCentralMirror] * std::tan(thetaBi[iCentralMirror]);
      lDetectorZ[iCentralMirror] = squareSizeBarrelCylinder;
      lAerogelZ[iCentralMirror] = std::sqrt(1.0 + mVal * mVal) * rMin * squareSizeBarrelCylinder / (std::sqrt(1.0 + mVal * mVal) * rMax - mVal * squareSizeBarrelCylinder);
      tRPlusG[iCentralMirror] = rMax - rMin;
      float t = std::tan(std::atan(mVal) + std::atan(squareSizeBarrelCylinder / (2.0 * rMax * std::sqrt(1.0 + mVal * mVal) - squareSizeBarrelCylinder * mVal)));
      thetaMax[iCentralMirror] = o2::constants::math::PIHalf - std::atan(t);
      thetaMin[iCentralMirror] = o2::constants::math::PIHalf + std::atan(t);
      bRichPhotodetectorCentralModuleHalfLength.value = rMin * t;
      aerogelRindex[iCentralMirror] = bRichRefractiveIndexSector[0];
      for (int i = iCentralMirror + 1; i < numberOfSectorsInZ; i++) {
        float parA = t;
        float parB = 2.0 * rMax / squareSizeZ;
        mVal = (std::sqrt(parA * parA * parB * parB + parB * parB - 1.0) + parA * parB * parB) / (parB * parB - 1.0);
        thetaMin[i] = o2::constants::math::PIHalf - std::atan(t);
        thetaMax[2 * iCentralMirror - i] = o2::constants::math::PIHalf + std::atan(t);
        t = std::tan(std::atan(mVal) + std::atan(squareSizeZ / (2.0 * rMax * std::sqrt(1.0 + mVal * mVal) - squareSizeZ * mVal)));
        thetaMax[i] = o2::constants::math::PIHalf - std::atan(t);
        thetaMin[2 * iCentralMirror - i] = o2::constants::math::PIHalf + std::atan(t);
        // Forward sectors
        thetaBi[i] = std::atan(mVal);
        r0Tilt[i] = rMax - squareSizeZ / 2.0 * std::sin(std::atan(mVal));
        z0Tilt[i] = r0Tilt[i] * std::tan(thetaBi[i]);
        lDetectorZ[i] = squareSizeZ;
        lAerogelZ[i] = std::sqrt(1.0 + mVal * mVal) * rMin * squareSizeZ / (std::sqrt(1.0 + mVal * mVal) * rMax - mVal * squareSizeZ);
        tRPlusG[i] = std::sqrt(1.0 + mVal * mVal) * (rMax - rMin) - mVal / 2.0 * (squareSizeZ + lAerogelZ[i]);
        aerogelRindex[i] = bRichRefractiveIndexSector[i - iCentralMirror];
        // Backword sectors
        thetaBi[2 * iCentralMirror - i] = -std::atan(mVal);
        r0Tilt[2 * iCentralMirror - i] = rMax - squareSizeZ / 2.0 * std::sin(std::atan(mVal));
        z0Tilt[2 * iCentralMirror - i] = -r0Tilt[i] * std::tan(thetaBi[i]);
        lDetectorZ[2 * iCentralMirror - i] = squareSizeZ;
        lAerogelZ[2 * iCentralMirror - i] = std::sqrt(1.0 + mVal * mVal) * rMin * squareSizeZ / (std::sqrt(1.0 + mVal * mVal) * rMax - mVal * squareSizeZ);
        tRPlusG[2 * iCentralMirror - i] = std::sqrt(1.0 + mVal * mVal) * (rMax - rMin) - mVal / 2.0 * (squareSizeZ + lAerogelZ[i]);
        aerogelRindex[2 * iCentralMirror - i] = bRichRefractiveIndexSector[i - iCentralMirror];
        bRichPhotodetectorCentralModuleHalfLength.value = rMin * t; // <-- At the end of the loop this will be the maximum z
      }
    } else { // Even number of sectors
      float twoHalfGap = 1.0;
      int iCentralMirror = static_cast<int>((numberOfSectorsInZ) / 2.0);
      float mVal = std::tan(0.0);
      float t = std::tan(std::atan(mVal) + std::atan(twoHalfGap / (2.0 * rMax * std::sqrt(1.0 + mVal * mVal) - twoHalfGap * mVal)));
      bRichPhotodetectorCentralModuleHalfLength.value = rMin * t;
      for (int i = iCentralMirror; i < numberOfSectorsInZ; i++) {
        float parA = t;
        float parB = 2.0 * rMax / squareSizeZ;
        mVal = (std::sqrt(parA * parA * parB * parB + parB * parB - 1.0) + parA * parB * parB) / (parB * parB - 1.0);
        thetaMin[i] = o2::constants::math::PIHalf - std::atan(t);
        thetaMax[2 * iCentralMirror - i - 1] = o2::constants::math::PIHalf + std::atan(t);
        t = std::tan(std::atan(mVal) + std::atan(squareSizeZ / (2.0 * rMax * std::sqrt(1.0 + mVal * mVal) - squareSizeZ * mVal)));
        thetaMax[i] = o2::constants::math::PIHalf - std::atan(t);
        thetaMin[2 * iCentralMirror - i - 1] = o2::constants::math::PIHalf + std::atan(t);
        // Forward sectors
        thetaBi[i] = std::atan(mVal);
        r0Tilt[i] = rMax - squareSizeZ / 2.0 * std::sin(std::atan(mVal));
        z0Tilt[i] = r0Tilt[i] * std::tan(thetaBi[i]);
        lDetectorZ[i] = squareSizeZ;
        lAerogelZ[i] = std::sqrt(1.0 + mVal * mVal) * rMin * squareSizeZ / (std::sqrt(1.0 + mVal * mVal) * rMax - mVal * squareSizeZ);
        tRPlusG[i] = std::sqrt(1.0 + mVal * mVal) * (rMax - rMin) - mVal / 2.0 * (squareSizeZ + lAerogelZ[i]);
        aerogelRindex[i] = bRichRefractiveIndexSector[i - iCentralMirror];
        // Backword sectors
        thetaBi[2 * iCentralMirror - i - 1] = -std::atan(mVal);
        r0Tilt[2 * iCentralMirror - i - 1] = rMax - squareSizeZ / 2.0 * std::sin(std::atan(mVal));
        z0Tilt[2 * iCentralMirror - i - 1] = -r0Tilt[i] * std::tan(thetaBi[i]);
        lDetectorZ[2 * iCentralMirror - i - 1] = squareSizeZ;
        lAerogelZ[2 * iCentralMirror - i - 1] = std::sqrt(1.0 + mVal * mVal) * rMin * squareSizeZ / (std::sqrt(1.0 + mVal * mVal) * rMax - mVal * squareSizeZ);
        tRPlusG[2 * iCentralMirror - i - 1] = std::sqrt(1.0 + mVal * mVal) * (rMax - rMin) - mVal / 2.0 * (squareSizeZ + lAerogelZ[i]);
        aerogelRindex[2 * iCentralMirror - i - 1] = bRichRefractiveIndexSector[i - iCentralMirror];
        bRichPhotodetectorCentralModuleHalfLength.value = rMin * t; // <-- At the end of the loop this will be the maximum z
      }
    }
    // Coordinate radiali layer considerati
    std::vector<float> r0Detector;
    std::vector<float> r0Aerogel;
    r0Detector.resize(numberOfSectorsInZ);
    r0Aerogel.resize(numberOfSectorsInZ);
    for (int i = 0; i < numberOfSectorsInZ; i++) {
      r0Detector[i] = r0Tilt[i]; // + (detector_thickness / 2.0) * std::cos(thetaBi[i]) NEGLIGIBLE
      detCenters[i].SetXYZ(r0Detector[i], 0, r0Detector[i] * std::tan(thetaBi[i]));
      r0Aerogel[i] = r0Tilt[i] - (tRPlusG[i] - bRichRadiatorThickness / 2.0) * std::cos(thetaBi[i]);
      radCenters[i].SetXYZ(r0Aerogel[i], 0, r0Aerogel[i] * std::tan(thetaBi[i]));
      angleCenters[i] = thetaBi[i];
      photodetrctorLength[i] = lDetectorZ[i];
      gapThickness[i] = (detCenters[i] - radCenters[i]).Mag() - bRichRadiatorThickness / 2.0;
    }
    // DEBUG
    // std::cout << std::endl << std::endl;
    // for (int i = 0; i < numberOfSectorsInZ; i++) {
    //  std::cout << "Sector " << i << ": Gap = " << gapThickness[i] << " cm, Index aerogel = " << aerogelRindex[i] << ", Gap thickness = " << gapThickness[i] << " cm" << std::endl;
    // }
    // std::cout << std::endl << std::endl;
  }

  // Configuration defined at init time
  o2::fastsim::GeometryContainer mGeoContainer;
  float mMagneticField = 0.0f;
  void init(o2::framework::InitContext& initContext)
  {
    mGeoContainer.init(initContext);

    const int nGeometries = mGeoContainer.getNumberOfConfigurations();
    pRandomNumberGenerator.SetSeed(0); // fully randomize
    mMagneticField = mGeoContainer.getFloatValue(0, "global", "magneticfield");

    if (flagRICHLoadDelphesLUTs) {
      for (int icfg = 0; icfg < nGeometries; ++icfg) {
        mSmearer.emplace_back(std::make_unique<o2::delphes::DelphesO2TrackSmearer>());
        mSmearer[icfg]->setCcdbManager(ccdb.operator->());
        std::map<std::string, std::string> globalConfiguration = mGeoContainer.getConfiguration(icfg, "global");
        for (const auto& entry : globalConfiguration) {
          int pdg = 0;
          if (entry.first.find("lut") != 0) {
            continue;
          }
          if (entry.first.find("lutEl") != std::string::npos) {
            pdg = kElectron;
          } else if (entry.first.find("lutMu") != std::string::npos) {
            pdg = kMuonMinus;
          } else if (entry.first.find("lutPi") != std::string::npos) {
            pdg = kPiPlus;
          } else if (entry.first.find("lutKa") != std::string::npos) {
            pdg = kKPlus;
          } else if (entry.first.find("lutPr") != std::string::npos) {
            pdg = kProton;
          } else if (entry.first.find("lutDe") != std::string::npos) {
            pdg = o2::constants::physics::kDeuteron;
          } else if (entry.first.find("lutTr") != std::string::npos) {
            pdg = o2::constants::physics::kTriton;
          } else if (entry.first.find("lutHe3") != std::string::npos) {
            pdg = o2::constants::physics::kHelium3;
          } else if (entry.first.find("lutAl") != std::string::npos) {
            pdg = o2::constants::physics::kAlpha;
          }

          std::string filename = entry.second;
          if (pdg == 0) {
            LOG(fatal) << "Unknown LUT entry " << entry.first << " for global configuration";
          }
          LOG(info) << "Loading LUT for pdg " << pdg << " for config " << icfg << " from provided file '" << filename << "'";
          if (filename.empty()) {
            LOG(warning) << "No LUT file passed for pdg " << pdg << ", skipping.";
          }
          // strip from leading/trailing spaces
          filename.erase(0, filename.find_first_not_of(" "));
          filename.erase(filename.find_last_not_of(" ") + 1);
          if (filename.empty()) {
            LOG(warning) << "No LUT file passed for pdg " << pdg << ", skipping.";
          }
          bool success = mSmearer[icfg]->loadTable(pdg, filename.c_str());
          if (!success) {
            LOG(fatal) << "Having issue with loading the LUT " << pdg << " " << filename;
          }
        }
      }
    }

    if (doQAplots) {
      const AxisSpec axisMomentum{static_cast<int>(nBinsP), 0.0f, +20.0f, "#it{p} (GeV/#it{c})"};
      const AxisSpec axisAngle{static_cast<int>(nBinsThetaRing), 0.0f, +0.30f, "Measured Cherenkov angle (rad)"};
      const AxisSpec axisEta{static_cast<int>(nBinsEta), -2.0f, +2.0f, "#eta"};
      const AxisSpec axisSector{static_cast<int>(bRichNumberOfSectors), -0.5f, static_cast<double>(bRichNumberOfSectors) - 0.5f, "Sector ID number"};
      const AxisSpec axisRingAngularResolution{static_cast<int>(nBinsAngularRes), 0.0f, +5.0f, "Ring angular resolution - rec. only (mrad)"};
      histos.add("h2dAngleVsMomentumBarrelRICH", "h2dAngleVsMomentumBarrelRICH", kTH2F, {axisMomentum, axisAngle});
      histos.add("h2dAngleVsEtaBarrelRICH", "h2dAngleVsEtaBarrelRICH", kTH2F, {axisEta, axisAngle});
      histos.add("h2dAngularResolutionVsMomentumBarrelRICH", "h2dAngularResolutionVsMomentumBarrelRICH", kTH2F, {axisMomentum, axisRingAngularResolution});
      histos.add("h2dAngularResolutionVsEtaBarrelRICH", "h2dAngularResolutionVsEtaBarrelRICH", kTH2F, {axisEta, axisRingAngularResolution});
      histos.add("hSectorID", "hSectorID", kTH1F, {axisSector});

      const int kNspec = 9; // electron, muon, pion, kaon, proton, deuteron, triton, helium3, alpha
      std::string particleNames1[kNspec] = {"#it{e}", "#it{#mu}", "#it{#pi}", "#it{K}", "#it{p}", "#it{d}", "#it{t}", "^{3}He", "#it{#alpha}"};
      std::string particleNames2[kNspec] = {"Elec", "Muon", "Pion", "Kaon", "Prot", "Deut", "Trit", "He3", "Al"};
      for (int iTrue = 0; iTrue < kNspec; iTrue++) {
        std::string nameTitleBarrelTrackRes = "h2dBarrelAngularResTrack" + particleNames2[iTrue] + "VsP";
        std::string nameTitleBarrelTotalRes = "h2dBarrelAngularResTotal" + particleNames2[iTrue] + "VsP";
        const AxisSpec axisTrackAngularRes{static_cast<int>(nBinsAngularRes), 0.0f, +5.0f, "Track angular resolution - " + particleNames1[iTrue] + " (mrad)"};
        const AxisSpec axisTotalAngularRes{static_cast<int>(nBinsAngularRes), 0.0f, +5.0f, "Total angular resolution - " + particleNames1[iTrue] + " (mrad)"};
        histos.add(nameTitleBarrelTrackRes.c_str(), nameTitleBarrelTrackRes.c_str(), kTH2F, {axisMomentum, axisTrackAngularRes});
        histos.add(nameTitleBarrelTotalRes.c_str(), nameTitleBarrelTotalRes.c_str(), kTH2F, {axisMomentum, axisTotalAngularRes});
      }
      for (int iTrue = 0; iTrue < kNspec; iTrue++) {
        for (int iHyp = 0; iHyp < kNspec; iHyp++) {
          std::string nameTitle = "h2dBarrelNsigmaTrue" + particleNames2[iTrue] + "Vs" + particleNames2[iHyp] + "Hypothesis";
          if (iTrue == iHyp) {
            const AxisSpec axisNsigmaCorrect{static_cast<int>(nBinsNsigmaCorrectSpecies), -10.0f, +10.0f, "N#sigma - True " + particleNames1[iTrue] + " vs " + particleNames1[iHyp] + " hypothesis"};
            histos.add(nameTitle.c_str(), nameTitle.c_str(), kTH2F, {axisMomentum, axisNsigmaCorrect});
          } else {
            const AxisSpec axisNsigmaWrong{static_cast<int>(nBinsNsigmaWrongSpecies), -10.0f, +10.0f, "N#sigma -  True " + particleNames1[iTrue] + " vs " + particleNames1[iHyp] + " hypothesis"};
            histos.add(nameTitle.c_str(), nameTitle.c_str(), kTH2F, {axisMomentum, axisNsigmaWrong});
          }
        }
      }
    }
    // Load all sector refractive indices
    /*bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector0AndNMinus1);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector1AndNMinus2);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector2AndNMinus3);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector3AndNMinus4);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector4AndNMinus5);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector5AndNMinus6);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector6AndNMinus7);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector7AndNMinus8);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector8AndNMinus9);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector9AndNMinus10);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector10AndNMinus11);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector11AndNMinus12);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector12AndNMinus13);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector13AndNMinus14);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector14AndNMinus15);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector15AndNMinus16);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector16AndNMinus17);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector17AndNMinus18);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector18AndNMinus19);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector19AndNMinus20);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector20AndNMinus21);*/
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector0);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector1);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector2);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector3);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector4);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector5);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector6);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector7);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector8);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector9);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector10);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector11);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector12);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector13);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector14);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector15);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector16);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector17);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector18);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector19);
    bRichRefractiveIndexSector.push_back(bRichRefractiveIndexSector20);

    // Update projective parameters
    updateProjectiveParameters();
  }

  /// check if particle reaches radiator
  /// \param track the input track
  /// \param radius the radius of the layer you're calculating the length to
  /// \param magneticField the magnetic field to use when propagating
  bool checkMagfieldLimit(o2::track::TrackParCov track, const float radius, const float magneticField)
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

  /// returns sector hit by the track (if any), -1 otherwise
  /// \param eta the pseudorapidity of the tarck (assuming primary vertex at origin)
  int findSector(const float eta)
  {
    float polar = 2.0 * std::atan(std::exp(-eta));
    for (int jSector = 0; jSector < bRichNumberOfSectors; jSector++) {
      if (polar > thetaMax[jSector] && polar < thetaMin[jSector]) {
        return jSector;
      }
    }
    return -1;
  }

  /// returns radiator radius in cm at considered eta (accounting for sector inclination)
  /// \param eta the pseudorapidity of the tarck (assuming primary vertex at origin)
  /// \param iSecor the index of the track RICH sector
  float radiusRipple(const float eta, const int iSecor)
  {
    if (iSecor > 0) {
      const float polar = 2.0 * std::atan(std::exp(-eta));
      const float rSecRich = radCenters[iSecor].X();
      const float zSecRich = radCenters[iSecor].Z();
      const float rSecRichSquared = rSecRich * rSecRich;
      const float zSecRichSquared = zSecRich * zSecRich;
      return (rSecRichSquared + zSecRichSquared) / (rSecRich + zSecRich / std::tan(polar));
    } else {
      return kErrorValue;
    }
  }

  /// returns Cherenkov angle in rad (above threshold) or bad flag (below threshold)
  /// \param momentum the momentum of the tarck
  /// \param mass the mass of the particle
  /// \param n the refractive index of the considered sector
  /// \param angle the output angle (passed by reference)
  /// \return true if particle is above the threshold and enough photons, false otherwise
  bool cherenkovAngle(const float momentum, const float mass, const float n, float& angle)
  {
    if (momentum > mass / std::sqrt(n * n - 1.0)) { // Check if particle is above the threshold
      // Calculate angle
      angle = std::acos(std::sqrt(momentum * momentum + mass * mass) / (momentum * n));

      // Mean number of detected photons (SiPM PDE is accountd in multiplicative factor!!!)
      const float meanNumberofDetectedPhotons = 230. * std::sin(angle) * std::sin(angle) * bRichRadiatorThickness;

      // Require at least 3 photons on average for real angle reconstruction
      static constexpr float kMinPhotons = 3.f;
      if (meanNumberofDetectedPhotons > kMinPhotons) {
        return true; // Particle is above the threshold and enough photons
      }
      return false; // Particle is above the threshold, but not enough photons
    }
    return false; // Particle is below the threshold
  }

  bool isOverTrhesholdInGasRadiator(const float momentum, const float mass)
  {
    if (momentum < mass / std::sqrt(gasRadiatorRindex * gasRadiatorRindex - 1.0)) { // Check if particle is above the threshold
      return false;
    }
    const float angle = std::acos(std::sqrt(momentum * momentum + mass * mass) / (momentum * gasRadiatorRindex));
    const float meanNumberofDetectedPhotons = 230. * std::sin(angle) * std::sin(angle) * gasRichRadiatorThickness;

    // Require at least 3 photons on average for real angle reconstruction
    static constexpr float kMinPhotons = 3.f;
    if (meanNumberofDetectedPhotons <= kMinPhotons) {
      return false;
    }
    return true;
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
  float angularResolution(float eta)
  {
    // Vectors for sampling (USE ANALYTICAL EXTRAPOLATION FOR BETTER RESULTS)
    static constexpr float kEtaSampling[] = {-2.000000, -1.909740, -1.731184, -1.552999, -1.375325, -1.198342, -1.022276, -0.847390, -0.673976, -0.502324, -0.332683, -0.165221, 0.000000, 0.165221, 0.332683, 0.502324, 0.673976, 0.847390, 1.022276, 1.198342, 1.375325, 1.552999, 1.731184, 1.909740, 2.000000};
    static constexpr float kResRingSamplingWithAbsWalls[] = {0.0009165, 0.000977, 0.001098, 0.001198, 0.001301, 0.001370, 0.001465, 0.001492, 0.001498, 0.001480, 0.001406, 0.001315, 0.001241, 0.001325, 0.001424, 0.001474, 0.001480, 0.001487, 0.001484, 0.001404, 0.001273, 0.001197, 0.001062, 0.000965, 0.0009165};
    static constexpr float kResRingSamplingWithoutAbsWalls[] = {0.0009165, 0.000977, 0.001095, 0.001198, 0.001300, 0.001369, 0.001468, 0.001523, 0.001501, 0.001426, 0.001299, 0.001167, 0.001092, 0.001179, 0.001308, 0.001407, 0.001491, 0.001508, 0.001488, 0.001404, 0.001273, 0.001196, 0.001061, 0.000965, 0.0009165};
    static constexpr int kSizeResVector = sizeof(kEtaSampling) / sizeof(kEtaSampling[0]);
    // Use binary search to find the lower and upper indices
    const int lowerIndex = std::lower_bound(kEtaSampling, kEtaSampling + kSizeResVector, eta) - kEtaSampling - 1;
    const int upperIndex = lowerIndex + 1;
    if (lowerIndex >= 0 && upperIndex < kSizeResVector) {
      // Resolution interpolation
      if (bRichFlagAbsorbingWalls) {
        float interpolatedResRing = interpolate(eta, kEtaSampling[lowerIndex], kEtaSampling[upperIndex], kResRingSamplingWithAbsWalls[lowerIndex], kResRingSamplingWithAbsWalls[upperIndex]);
        // std::cout << "Interpolated y value: " << interpolatedY << std::endl;
        return interpolatedResRing;
      } else {
        float interpolatedResRing = interpolate(eta, kEtaSampling[lowerIndex], kEtaSampling[upperIndex], kResRingSamplingWithoutAbsWalls[lowerIndex], kResRingSamplingWithoutAbsWalls[upperIndex]);
        // std::cout << "Interpolated y value: " << interpolatedY << std::endl;
        return interpolatedResRing;
      }
    } else {
      // std::cout << "Unable to interpolate. Target x value is outside the range of available data." << std::endl;
      return kErrorValue;
    }
  }

  /// To account border effects in bRICH
  /// Approximation for analytic calculation: track along projective sector normal (reasonable)
  /// \param eta the pseudorapidity of the tarck (assuming primary vertex at origin)
  /// \param tileZlength the Z-length of the photodetector plane of the considered sector
  /// \param radius the nominal radius of the Cherenkov ring
  float fractionPhotonsProjectiveRICH(float eta, float tileZlength, float radius)
  {
    const float polar = 2.0 * std::atan(std::exp(-eta));
    int iSecor = 0;
    bool flagSector = false;
    for (int jSector = 0; jSector < bRichNumberOfSectors; jSector++) {
      if (polar > thetaMax[jSector] && polar < thetaMin[jSector]) {
        flagSector = true;
        iSecor = jSector;
      }
    }
    if (flagSector == false) {
      return kErrorValue; // <-- Returning negative value
    }
    float rSecRich = radCenters[iSecor].X();
    float zSecRich = radCenters[iSecor].Z();
    // float rSecTof = detCenters[iSecor].X();
    // float zSecTof = detCenters[iSecor].Z();
    const float rSecRichSquared = rSecRich * rSecRich;
    const float zSecRichSquared = zSecRich * zSecRich;
    const float radiusRipple = (rSecRichSquared + zSecRichSquared) / (rSecRich + zSecRich / std::tan(polar));
    const float zRipple = radiusRipple / std::tan(polar);
    const float absZ = std::hypot(radiusRipple - rSecRich, zRipple - zSecRich);
    float fraction = 1.;
    if (tileZlength / 2. - absZ < radius) {
      fraction = fraction - (1. / o2::constants::math::PI) * std::acos((tileZlength / 2. - absZ) / radius);
      if (tileZlength / 2. + absZ < radius) {
        fraction = fraction - (1. / o2::constants::math::PI) * std::acos((tileZlength / 2. + absZ) / radius);
      }
    }
    return fraction;
  }

  /// returns refined ring angular resolution
  /// \param eta the pseudorapidity of the tarck (assuming primary vertex at origin)
  /// \param n the refractive index of the considered sector
  /// \param nGas the refractive index of the gas filling the RICH proximity gap
  /// \param thicknessRad the radiator thickness in cm
  /// \param thicknessGas the proximity gap thickness of the considered sector in cm
  /// \param pixelSize the SiPM pixel size in cm
  /// \param thetaCherenkov the Cherenkov angle for the considered track
  /// \param tileZlength the Z-length of the photodetector plane of the considered sector
  float extractRingAngularResolution(const float eta, const float n, const float nGas, const float thicknessRad, const float thicknessGas, const float pixelSize, const float thetaCherenkov, const float tileZlength)
  {
    // Check if input angle is error value
    if (thetaCherenkov <= kErrorValue + 1)
      return kErrorValue;
    // Parametrization variables (notation from https://doi.org/10.1016/0168-9002(94)90532-0)
    const float phiC = 0.;
    const float thetaP = 0.;
    const float sinThetaCherenkov = std::sin(thetaCherenkov);
    const float cosThetaCherenkov = std::cos(thetaCherenkov);
    const float xP = std::sin(thetaP);
    const float zP = std::cos(thetaP);
    const float sinPhi = std::sin(phiC);
    const float cosPhi = std::cos(phiC);
    // const float ze = thicknessRad / (2. * zP);
    const float aZ = zP * cosThetaCherenkov - xP * sinThetaCherenkov * cosPhi;
    const float e3z = std::sqrt(aZ * aZ + (nGas / n) * (nGas / n) - 1.);
    const float z = thicknessGas;
    const float alpha = e3z / aZ;
    const float etac = e3z * n * n;
    const float k = thicknessRad / (2. * z);
    const float m = 1. / (n * n);
    const float lambda = (1. + k * alpha * alpha * alpha) / (1. + k * alpha);
    // Derivative d(thetaCherenkov)/dx
    const float temp1 = etac / z;
    const float temp2 = alpha * e3z * cosPhi;
    const float temp3 = xP * sinThetaCherenkov * sinPhi * sinPhi;
    const float dThetaX = temp1 * (temp2 - temp3);
    // Derivative d(thetaCherenkov)/dy
    const float temp4 = etac * sinPhi / z;
    const float temp5 = cosThetaCherenkov - zP * (1 - m) / aZ;
    const float dThetaY = temp4 * temp5;
    // Derivative d(thetaCherenkov)/dze
    const float temp8 = etac * sinThetaCherenkov;
    const float temp9 = z * (1.0 + k * alpha * alpha * alpha * n * n);
    const float temp10 = alpha * alpha * (1.0 - xP * xP * sinPhi * sinPhi) + lambda * xP * xP * sinPhi * sinPhi;
    const float dThetaZe = temp8 * temp10 / temp9;
    // Derivative d(thetaCherenkov)/dn
    const float dThetaN = zP / (n * aZ * sinThetaCherenkov);
    // RMS wavelength (using Sellmeier formula with measured aerogel parameters)
    const float a0 = 0.0616;
    const float w0 = 56.5;
    const float n0 = 1.0307;
    const float wMean = 436.0; // 450.0;//484.9;//426.0; //450;//485;//450; // nm  //450   426 493.5 // 426.0
    const float scaling = n / n0;
    const float temp11 = a0 * w0 * w0 * wMean;
    const float temp12 = std::pow(wMean * wMean - w0 * w0, 1.5);
    const float temp13 = std::sqrt(wMean * wMean * (1.0 + a0) - w0 * w0);
    const float dNw = temp11 / (temp12 * temp13) * scaling;
    const float sigmaW = 115.0;
    // RMS x y
    const float sigmaThetaX = pixelSize / std::sqrt(12.0);
    const float sigmaThetaY = pixelSize / std::sqrt(12.0);
    // RMS emission point
    const float sigmaThetaZe = thicknessRad / (zP * std::sqrt(12.0));
    // Single photon angular resolution
    const float resX = dThetaX * sigmaThetaX;
    const float resY = dThetaY * sigmaThetaY;
    const float resZe = dThetaZe * sigmaThetaZe;
    const float resN = dThetaN * dNw * sigmaW;
    const float resXsquared = resX * resX;
    const float resYsquared = resY * resY;
    const float resZesquared = resZe * resZe;
    const float resNsquared = resN * resN;
    const float singlePhotonAngularResolution = std::sqrt(resXsquared + resYsquared + resZesquared + resNsquared);
    // Fraction of photons in rings (to account loss close to sector boundary in projective geometry)
    const float sinThetaCherenkovSquared = sinThetaCherenkov * sinThetaCherenkov;
    const float radius = (thicknessRad / 2.0) * std::tan(thetaCherenkov) + thicknessGas * std::tan(std::asin((n / nGas) * sinThetaCherenkov));
    const float n0Photons = 24. * thicknessRad / 2.; // photons for N = 1.03 at saturation ( 24/2 factor per radiator cm )
    const float sinAngle = std::sin(std::acos(1. / 1.03));
    const float sinAngleSquared = sinAngle * sinAngle;
    const float multiplicitySpectrumFactor = sinThetaCherenkovSquared / sinAngleSquared; // scale multiplicity w.r.t. N = 1.03 at saturation
    // Considering average resolution (integrated over the sector)
    // float nPhotons = (tileZlength / 2.0 > radius) ? n0Photons * multiplicitySpectrumFactor * (1.-(2.0*radius)/(o2::constants::math::PI*tileZlength)) : n0Photons * multiplicitySpectrumFactor * (1.-(2.0*radius)/(o2::constants::math::PI*tileZlength) - (2.0/(tileZlength*o2::constants::math::PI))*(-(tileZlength/(2.0))*std::acos(tileZlength/(2.0*radius)) + radius*std::sqrt(1.-std::pow(tileZlength/(2.0*radius),2.0))));
    // Considering "exact" resolution (eta by eta)
    const float nPhotons = n0Photons * multiplicitySpectrumFactor * fractionPhotonsProjectiveRICH(eta, tileZlength, radius);
    if (nPhotons <= kErrorValue + 1)
      return kErrorValue;
    // Ring angular resolution
    const float ringAngularResolution = singlePhotonAngularResolution / std::sqrt(nPhotons);
    return ringAngularResolution;
  }

  /// returns track angular resolution
  /// \param pt the transverse momentum of the tarck
  /// \param eta the pseudorapidity of the tarck
  /// \param trackPtResolution the absolute resolution on pt
  /// \param trackPtResolution the absolute resolution on eta
  /// \param mass the mass of the particle
  /// \param refractiveIndex the refractive index of the radiator
  double calculateTrackAngularResolutionAdvanced(const float pt, const float eta, const float trackPtResolution, const float trackEtaResolution, const float mass, const float refractiveIndex)
  {
    // Compute tracking contribution to timing using the error propagation formula
    // Uses light speed in m/ps, magnetic field in T (*0.1 for conversion kGauss -> T)
    const double a0 = mass * mass;
    const double a1 = refractiveIndex;
    const double a1Squared = a1 * a1;
    const float ptCoshEta = pt * std::cosh(eta);
    const float ptCoshEtaSquared = ptCoshEta * ptCoshEta;
    const double dThetaOndPt = a0 / (pt * std::sqrt(a0 + ptCoshEtaSquared) * std::sqrt(ptCoshEtaSquared * (a1Squared - 1.0) - a0));
    const double dThetaOndEta = (a0 * std::tanh(eta)) / (std::sqrt(a0 + ptCoshEtaSquared) * std::sqrt(ptCoshEtaSquared * (a1Squared - 1.0) - a0));
    const double trackAngularResolution = std::hypot(std::fabs(dThetaOndPt) * trackPtResolution, std::fabs(dThetaOndEta) * trackEtaResolution);
    return trackAngularResolution;
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::OTFLUTConfigId>::iterator const& collision,
               soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels> const& tracks,
               aod::McParticles const&,
               aod::McCollisions const&)
  {
    const o2::dataformats::VertexBase pvVtx({collision.posX(), collision.posY(), collision.posZ()},
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
      auto fillDummyValues = [&](bool gasRich = false) {
        upgradeRich(kErrorValue, kErrorValue, kErrorValue, kErrorValue, kErrorValue, kErrorValue, kErrorValue, kErrorValue, kErrorValue);
        upgradeRichSignal(false, false, false, false, false, false, false, false, false, false, gasRich);
      };

      // first step: find precise arrival time (if any)
      // --- convert track into perfect track
      if (!track.has_mcParticle()) { // should always be OK but check please
        fillDummyValues();
        continue;
      }

      const auto& mcParticle = track.mcParticle();
      o2::track::TrackParCov o2track = o2::upgrade::convertMCParticleToO2Track(mcParticle, pdg);

      // float xPv = kErrorValue;
      if (o2track.propagateToDCA(mcPvVtx, mMagneticField)) {
        // xPv = o2track.getX();
      }

      // get particle to calculate Cherenkov angle and resolution
      auto pdgInfo = pdg->GetParticle(mcParticle.pdgCode());
      if (pdgInfo == nullptr) {
        fillDummyValues();
        continue;
      }

      // find track bRICH sector
      const int iSecor = findSector(o2track.getEta());
      if (iSecor < 0) {
        fillDummyValues();
        continue;
      }

      const bool expectedAngleBarrelGasRichOk = isOverTrhesholdInGasRadiator(o2track.getP(), pdgInfo->Mass());
      float expectedAngleBarrelRich = kErrorValue;
      const bool expectedAngleBarrelRichOk = cherenkovAngle(o2track.getP(), pdgInfo->Mass(), aerogelRindex[iSecor], expectedAngleBarrelRich);
      if (!expectedAngleBarrelRichOk) {
        fillDummyValues(expectedAngleBarrelGasRichOk);
        continue; // Particle is below the threshold or not enough photons
      }
      // float barrelRICHAngularResolution = angularResolution(o2track.getEta());
      const float barrelRICHAngularResolution = extractRingAngularResolution(o2track.getEta(), aerogelRindex[iSecor], bRichGapRefractiveIndex, bRichRadiatorThickness, gapThickness[iSecor], bRICHPixelSize, expectedAngleBarrelRich, photodetrctorLength[iSecor]);
      const float projectiveRadiatorRadius = radiusRipple(o2track.getEta(), iSecor);
      bool flagReachesRadiator = false;
      if (projectiveRadiatorRadius > kErrorValue + 1.) {
        flagReachesRadiator = checkMagfieldLimit(o2track, projectiveRadiatorRadius, mMagneticField);
      }
      /// DISCLAIMER: Exact extrapolation of angular resolution would require track propagation
      ///             to the RICH radiator (accounting sector inclination) in terms of (R,z).
      ///             The extrapolation with Eta is correct only if the primary vertex is at origin.
      ///             Discrepancies may be negligible, but would be more rigorous if propagation tool is available

      // Smear with expected resolutions
      const float measuredAngleBarrelRich = pRandomNumberGenerator.Gaus(expectedAngleBarrelRich, barrelRICHAngularResolution);

      // Now we calculate the expected Cherenkov angle following certain mass hypotheses
      // and the (imperfect!) reconstructed track parametrizations
      auto recoTrack = getTrackParCov(track);
      if (recoTrack.propagateToDCA(pvVtx, mMagneticField)) {
        // xPv = recoTrack.getX();
      }

      // Straight to Nsigma
      static constexpr int kNspecies = 9;
      static constexpr int kEl = 0;
      static constexpr int kMu = 1;
      static constexpr int kPi = 2;
      static constexpr int kKa = 3;
      static constexpr int kPr = 4;
      static constexpr int kDe = 5;
      static constexpr int kTr = 6;
      static constexpr int kHe3 = 7;
      static constexpr int kAl = 8;
      float nSigmaBarrelRich[kNspecies] = {kErrorValue, kErrorValue, kErrorValue, kErrorValue, kErrorValue, kErrorValue, kErrorValue, kErrorValue, kErrorValue};
      bool signalBarrelRich[kNspecies] = {false, false, false, false, false, false, false, false, false};
      float deltaThetaBarrelRich[kNspecies]; //, nSigmaBarrelRich[kNspecies];
      static constexpr int kParticlePdgs[kNspecies] = {kElectron,
                                                       kMuonMinus,
                                                       kPiPlus,
                                                       kKPlus,
                                                       kProton,
                                                       o2::constants::physics::kDeuteron,
                                                       o2::constants::physics::kTriton,
                                                       o2::constants::physics::kHelium3,
                                                       o2::constants::physics::kAlpha};
      static constexpr float kParticleMasses[kNspecies] = {o2::track::pid_constants::sMasses[o2::track::PID::Electron],
                                                           o2::track::pid_constants::sMasses[o2::track::PID::Muon],
                                                           o2::track::pid_constants::sMasses[o2::track::PID::Pion],
                                                           o2::track::pid_constants::sMasses[o2::track::PID::Kaon],
                                                           o2::track::pid_constants::sMasses[o2::track::PID::Proton],
                                                           o2::track::pid_constants::sMasses[o2::track::PID::Deuteron],
                                                           o2::track::pid_constants::sMasses[o2::track::PID::Triton],
                                                           o2::track::pid_constants::sMasses[o2::track::PID::Helium3],
                                                           o2::track::pid_constants::sMasses[o2::track::PID::Alpha]};

      for (int ii = 0; ii < kNspecies; ii++) { // Loop on the particle hypotheses

        float hypothesisAngleBarrelRich = kErrorValue;
        const bool hypothesisAngleBarrelRichOk = cherenkovAngle(recoTrack.getP(), kParticleMasses[ii], aerogelRindex[iSecor], hypothesisAngleBarrelRich);
        signalBarrelRich[ii] = hypothesisAngleBarrelRichOk; // Particle is above the threshold and enough photons

        // Evaluate total sigma (layer + tracking resolution)
        float barrelTotalAngularReso = barrelRICHAngularResolution;
        if (flagIncludeTrackAngularRes) {
          const float transverseMomentum = recoTrack.getP() / std::cosh(recoTrack.getEta());
          double ptResolution = transverseMomentum * transverseMomentum * std::sqrt(recoTrack.getSigma1Pt2());
          double etaResolution = std::fabs(std::sin(2.0 * std::atan(std::exp(-recoTrack.getEta())))) * std::sqrt(recoTrack.getSigmaTgl2());
          if (flagRICHLoadDelphesLUTs) {
            if (mSmearer[collision.lutConfigId()]->hasTable(kParticlePdgs[ii])) {
              ptResolution = mSmearer[collision.lutConfigId()]->getAbsPtRes(kParticlePdgs[ii], dNdEta, recoTrack.getEta(), transverseMomentum);
              etaResolution = mSmearer[collision.lutConfigId()]->getAbsEtaRes(kParticlePdgs[ii], dNdEta, recoTrack.getEta(), transverseMomentum);
            }
          }
          // cout << endl <<  "Pt resolution: " << ptResolution << ", Eta resolution: " << etaResolution << endl << endl;
          const float barrelTrackAngularReso = calculateTrackAngularResolutionAdvanced(recoTrack.getP() / std::cosh(recoTrack.getEta()), recoTrack.getEta(), ptResolution, etaResolution, kParticleMasses[ii], aerogelRindex[iSecor]);
          barrelTotalAngularReso = std::hypot(barrelRICHAngularResolution, barrelTrackAngularReso);
          if (doQAplots &&
              hypothesisAngleBarrelRich > kErrorValue + 1. &&
              measuredAngleBarrelRich > kErrorValue + 1. &&
              barrelRICHAngularResolution > kErrorValue + 1. &&
              flagReachesRadiator) {
            switch (mcParticle.pdgCode()) {
              case kParticlePdgs[kEl]:  // Electron
              case -kParticlePdgs[kEl]: // Positron
                if (ii == kEl) {
                  histos.fill(HIST("h2dBarrelAngularResTrackElecVsP"), recoTrack.getP(), 1000.0 * barrelTrackAngularReso);
                  histos.fill(HIST("h2dBarrelAngularResTotalElecVsP"), recoTrack.getP(), 1000.0 * barrelTotalAngularReso);
                }
                break;
              case kParticlePdgs[kMu]:  // Muon
              case -kParticlePdgs[kMu]: // AntiMuon
                if (ii == kMu) {
                  histos.fill(HIST("h2dBarrelAngularResTrackMuonVsP"), recoTrack.getP(), 1000.0 * barrelTrackAngularReso);
                  histos.fill(HIST("h2dBarrelAngularResTotalMuonVsP"), recoTrack.getP(), 1000.0 * barrelTotalAngularReso);
                }
                break;
              case kParticlePdgs[kPi]:  // Pion
              case -kParticlePdgs[kPi]: // AntiPion
                if (ii == kPi) {
                  histos.fill(HIST("h2dBarrelAngularResTrackPionVsP"), recoTrack.getP(), 1000.0 * barrelTrackAngularReso);
                  histos.fill(HIST("h2dBarrelAngularResTotalPionVsP"), recoTrack.getP(), 1000.0 * barrelTotalAngularReso);
                }
                break;
              case kParticlePdgs[kKa]:  // Kaon
              case -kParticlePdgs[kKa]: // AntiKaon
                if (ii == kKa) {
                  histos.fill(HIST("h2dBarrelAngularResTrackKaonVsP"), recoTrack.getP(), 1000.0 * barrelTrackAngularReso);
                  histos.fill(HIST("h2dBarrelAngularResTotalKaonVsP"), recoTrack.getP(), 1000.0 * barrelTotalAngularReso);
                }
                break;
              case kParticlePdgs[kPr]:  // Proton
              case -kParticlePdgs[kPr]: // AntiProton
                if (ii == kPr) {
                  histos.fill(HIST("h2dBarrelAngularResTrackProtVsP"), recoTrack.getP(), 1000.0 * barrelTrackAngularReso);
                  histos.fill(HIST("h2dBarrelAngularResTotalProtVsP"), recoTrack.getP(), 1000.0 * barrelTotalAngularReso);
                }
                break;
              case kParticlePdgs[kDe]:  // Deuteron
              case -kParticlePdgs[kDe]: // AntiDeuteron
                if (ii == kDe) {
                  histos.fill(HIST("h2dBarrelAngularResTrackDeutVsP"), recoTrack.getP(), 1000.0 * barrelTrackAngularReso);
                  histos.fill(HIST("h2dBarrelAngularResTotalDeutVsP"), recoTrack.getP(), 1000.0 * barrelTotalAngularReso);
                }
                break;
              case kParticlePdgs[kTr]:  // Triton
              case -kParticlePdgs[kTr]: // AntiTriton
                if (ii == kTr) {
                  histos.fill(HIST("h2dBarrelAngularResTrackTritVsP"), recoTrack.getP(), 1000.0 * barrelTrackAngularReso);
                  histos.fill(HIST("h2dBarrelAngularResTotalTritVsP"), recoTrack.getP(), 1000.0 * barrelTotalAngularReso);
                }
                break;
              case kParticlePdgs[kHe3]:  // Helium3
              case -kParticlePdgs[kHe3]: // AntiHelium3
                if (ii == kHe3) {
                  histos.fill(HIST("h2dBarrelAngularResTrackHe3VsP"), recoTrack.getP(), 1000.0 * barrelTrackAngularReso);
                  histos.fill(HIST("h2dBarrelAngularResTotalHe3VsP"), recoTrack.getP(), 1000.0 * barrelTotalAngularReso);
                }
                break;
              case kParticlePdgs[kAl]:  // Alpha
              case -kParticlePdgs[kAl]: // AntiAlpha
                if (ii == kAl) {
                  histos.fill(HIST("h2dBarrelAngularResTrackAlVsP"), recoTrack.getP(), 1000.0 * barrelTrackAngularReso);
                  histos.fill(HIST("h2dBarrelAngularResTotalAlVsP"), recoTrack.getP(), 1000.0 * barrelTotalAngularReso);
                }
                break;
              default:
                break;
            }
          }
        }

        /// DISCLAIMER: here tracking is accounted only for momentum value, but not for track parameters at impact point on the
        ///             RICH radiator, since exact resolution would require photon generation and transport to photodetector.
        ///             Effects are expected to be negligible (a few tenths of a milliradian) but further studies are required !
        if (hypothesisAngleBarrelRich > kErrorValue + 1. &&
            measuredAngleBarrelRich > kErrorValue + 1. &&
            barrelRICHAngularResolution > kErrorValue + 1. &&
            flagReachesRadiator) {
          deltaThetaBarrelRich[ii] = hypothesisAngleBarrelRich - measuredAngleBarrelRich;
          nSigmaBarrelRich[ii] = deltaThetaBarrelRich[ii] / barrelTotalAngularReso;
        }
      }

      // Fill histograms
      if (doQAplots) {
        if (measuredAngleBarrelRich > kErrorValue + 1. &&
            barrelRICHAngularResolution > kErrorValue + 1. &&
            flagReachesRadiator) {
          histos.fill(HIST("h2dAngleVsMomentumBarrelRICH"), recoTrack.getP(), measuredAngleBarrelRich);
          histos.fill(HIST("h2dAngleVsEtaBarrelRICH"), recoTrack.getEta(), measuredAngleBarrelRich);
          histos.fill(HIST("h2dAngularResolutionVsMomentumBarrelRICH"), recoTrack.getP(), 1000 * barrelRICHAngularResolution);
          histos.fill(HIST("h2dAngularResolutionVsEtaBarrelRICH"), recoTrack.getEta(), 1000 * barrelRICHAngularResolution);
          histos.fill(HIST("hSectorID"), iSecor);

          switch (mcParticle.pdgCode()) {
            case kParticlePdgs[kEl]:  // Electron
            case -kParticlePdgs[kEl]: // Positron
              histos.fill(HIST("h2dBarrelNsigmaTrueElecVsElecHypothesis"), recoTrack.getP(), nSigmaBarrelRich[0]);
              histos.fill(HIST("h2dBarrelNsigmaTrueElecVsMuonHypothesis"), recoTrack.getP(), nSigmaBarrelRich[1]);
              histos.fill(HIST("h2dBarrelNsigmaTrueElecVsPionHypothesis"), recoTrack.getP(), nSigmaBarrelRich[2]);
              histos.fill(HIST("h2dBarrelNsigmaTrueElecVsKaonHypothesis"), recoTrack.getP(), nSigmaBarrelRich[3]);
              histos.fill(HIST("h2dBarrelNsigmaTrueElecVsProtHypothesis"), recoTrack.getP(), nSigmaBarrelRich[4]);
              break;
            case kParticlePdgs[kMu]:  // Muon
            case -kParticlePdgs[kMu]: // AntiMuon
              histos.fill(HIST("h2dBarrelNsigmaTrueMuonVsElecHypothesis"), recoTrack.getP(), nSigmaBarrelRich[0]);
              histos.fill(HIST("h2dBarrelNsigmaTrueMuonVsMuonHypothesis"), recoTrack.getP(), nSigmaBarrelRich[1]);
              histos.fill(HIST("h2dBarrelNsigmaTrueMuonVsPionHypothesis"), recoTrack.getP(), nSigmaBarrelRich[2]);
              histos.fill(HIST("h2dBarrelNsigmaTrueMuonVsKaonHypothesis"), recoTrack.getP(), nSigmaBarrelRich[3]);
              histos.fill(HIST("h2dBarrelNsigmaTrueMuonVsProtHypothesis"), recoTrack.getP(), nSigmaBarrelRich[4]);
              break;
            case kParticlePdgs[kPi]:  // Pion
            case -kParticlePdgs[kPi]: // AntiPion
              histos.fill(HIST("h2dBarrelNsigmaTruePionVsElecHypothesis"), recoTrack.getP(), nSigmaBarrelRich[0]);
              histos.fill(HIST("h2dBarrelNsigmaTruePionVsMuonHypothesis"), recoTrack.getP(), nSigmaBarrelRich[1]);
              histos.fill(HIST("h2dBarrelNsigmaTruePionVsPionHypothesis"), recoTrack.getP(), nSigmaBarrelRich[2]);
              histos.fill(HIST("h2dBarrelNsigmaTruePionVsKaonHypothesis"), recoTrack.getP(), nSigmaBarrelRich[3]);
              histos.fill(HIST("h2dBarrelNsigmaTruePionVsProtHypothesis"), recoTrack.getP(), nSigmaBarrelRich[4]);
              break;
            case kParticlePdgs[kKa]:  // Kaon
            case -kParticlePdgs[kKa]: // AntiKaon
              histos.fill(HIST("h2dBarrelNsigmaTrueKaonVsElecHypothesis"), recoTrack.getP(), nSigmaBarrelRich[0]);
              histos.fill(HIST("h2dBarrelNsigmaTrueKaonVsMuonHypothesis"), recoTrack.getP(), nSigmaBarrelRich[1]);
              histos.fill(HIST("h2dBarrelNsigmaTrueKaonVsPionHypothesis"), recoTrack.getP(), nSigmaBarrelRich[2]);
              histos.fill(HIST("h2dBarrelNsigmaTrueKaonVsKaonHypothesis"), recoTrack.getP(), nSigmaBarrelRich[3]);
              histos.fill(HIST("h2dBarrelNsigmaTrueKaonVsProtHypothesis"), recoTrack.getP(), nSigmaBarrelRich[4]);
              break;
            case kParticlePdgs[kPr]:  // Proton
            case -kParticlePdgs[kPr]: // AntiProton
              histos.fill(HIST("h2dBarrelNsigmaTrueProtVsElecHypothesis"), recoTrack.getP(), nSigmaBarrelRich[0]);
              histos.fill(HIST("h2dBarrelNsigmaTrueProtVsMuonHypothesis"), recoTrack.getP(), nSigmaBarrelRich[1]);
              histos.fill(HIST("h2dBarrelNsigmaTrueProtVsPionHypothesis"), recoTrack.getP(), nSigmaBarrelRich[2]);
              histos.fill(HIST("h2dBarrelNsigmaTrueProtVsKaonHypothesis"), recoTrack.getP(), nSigmaBarrelRich[3]);
              histos.fill(HIST("h2dBarrelNsigmaTrueProtVsProtHypothesis"), recoTrack.getP(), nSigmaBarrelRich[4]);
              histos.fill(HIST("h2dBarrelNsigmaTrueProtVsDeutHypothesis"), recoTrack.getP(), nSigmaBarrelRich[5]);
              break;
            case kParticlePdgs[kDe]:  // Deuteron
            case -kParticlePdgs[kDe]: // AntiDeuteron
              histos.fill(HIST("h2dBarrelNsigmaTrueDeutVsProtHypothesis"), recoTrack.getP(), nSigmaBarrelRich[4]);
              histos.fill(HIST("h2dBarrelNsigmaTrueDeutVsDeutHypothesis"), recoTrack.getP(), nSigmaBarrelRich[5]);
              histos.fill(HIST("h2dBarrelNsigmaTrueDeutVsTritHypothesis"), recoTrack.getP(), nSigmaBarrelRich[6]);
              histos.fill(HIST("h2dBarrelNsigmaTrueDeutVsHe3Hypothesis"), recoTrack.getP(), nSigmaBarrelRich[7]);
              histos.fill(HIST("h2dBarrelNsigmaTrueDeutVsAlHypothesis"), recoTrack.getP(), nSigmaBarrelRich[8]);
              break;
            case kParticlePdgs[kTr]:  // Triton
            case -kParticlePdgs[kTr]: // AntiTriton
              histos.fill(HIST("h2dBarrelNsigmaTrueTritVsProtHypothesis"), recoTrack.getP(), nSigmaBarrelRich[4]);
              histos.fill(HIST("h2dBarrelNsigmaTrueTritVsDeutHypothesis"), recoTrack.getP(), nSigmaBarrelRich[5]);
              histos.fill(HIST("h2dBarrelNsigmaTrueTritVsTritHypothesis"), recoTrack.getP(), nSigmaBarrelRich[6]);
              histos.fill(HIST("h2dBarrelNsigmaTrueTritVsHe3Hypothesis"), recoTrack.getP(), nSigmaBarrelRich[7]);
              histos.fill(HIST("h2dBarrelNsigmaTrueTritVsAlHypothesis"), recoTrack.getP(), nSigmaBarrelRich[8]);
              break;
            case kParticlePdgs[kHe3]:  // Helium3
            case -kParticlePdgs[kHe3]: // AntiHelium3
              histos.fill(HIST("h2dBarrelNsigmaTrueHe3VsDeutHypothesis"), recoTrack.getP(), nSigmaBarrelRich[5]);
              histos.fill(HIST("h2dBarrelNsigmaTrueHe3VsTritHypothesis"), recoTrack.getP(), nSigmaBarrelRich[6]);
              histos.fill(HIST("h2dBarrelNsigmaTrueHe3VsHe3Hypothesis"), recoTrack.getP(), nSigmaBarrelRich[7]);
              histos.fill(HIST("h2dBarrelNsigmaTrueHe3VsAlHypothesis"), recoTrack.getP(), nSigmaBarrelRich[8]);
              break;
            case kParticlePdgs[kAl]:  // Alpha
            case -kParticlePdgs[kAl]: // AntiAlpha
              histos.fill(HIST("h2dBarrelNsigmaTrueAlVsDeutHypothesis"), recoTrack.getP(), nSigmaBarrelRich[5]);
              histos.fill(HIST("h2dBarrelNsigmaTrueAlVsTritHypothesis"), recoTrack.getP(), nSigmaBarrelRich[6]);
              histos.fill(HIST("h2dBarrelNsigmaTrueAlVsHe3Hypothesis"), recoTrack.getP(), nSigmaBarrelRich[7]);
              histos.fill(HIST("h2dBarrelNsigmaTrueAlVsAlHypothesis"), recoTrack.getP(), nSigmaBarrelRich[8]);
              break;
            default:
              break;
          }
        }
      }

      // Sigmas have been fully calculated. Please populate the NSigma helper table (once per track)
      upgradeRich(nSigmaBarrelRich[0], nSigmaBarrelRich[1], nSigmaBarrelRich[2], nSigmaBarrelRich[3], nSigmaBarrelRich[4], nSigmaBarrelRich[5], nSigmaBarrelRich[6], nSigmaBarrelRich[7], nSigmaBarrelRich[8]);
      upgradeRichSignal(expectedAngleBarrelRichOk, signalBarrelRich[0], signalBarrelRich[1], signalBarrelRich[2], signalBarrelRich[3], signalBarrelRich[4], signalBarrelRich[5], signalBarrelRich[6], signalBarrelRich[7], signalBarrelRich[8], expectedAngleBarrelGasRichOk);
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyRichPid>(cfgc)};
}
