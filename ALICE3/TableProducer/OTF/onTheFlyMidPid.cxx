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
/// \file onTheFlyMidPid.cxx
///
/// \brief This task goes straight from a combination of track table and mcParticles
///        and a custom MID configuration to a table of muon id for the particles
///        being analysed.
///
/// \author David Dobrigkeit Chinellato, UNICAMP
/// \author Antonio Ortiz, UNAM
/// \since  May 22, 2024
///

#include "GeometryContainer.h"

#include "ALICE3/Core/FlatTrackSmearer.h"
#include "ALICE3/Core/TrackUtilities.h"
#include "ALICE3/DataModel/OTFCollision.h"
#include "ALICE3/DataModel/OTFMID.h"
#include "Common/Core/trackUtilities.h"

#include <CCDB/BasicCCDBManager.h>
#include <CommonConstants/PhysicsConstants.h>
#include <DetectorsBase/Propagator.h>
#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisHelpers.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/HistogramRegistry.h>
#include <Framework/HistogramSpec.h>
#include <Framework/InitContext.h>
#include <Framework/O2DatabasePDGPlugin.h>
#include <Framework/OutputObjHeader.h>
#include <Framework/runDataProcessing.h>
#include <ReconstructionDataFormats/Track.h>
#include <ReconstructionDataFormats/TrackParametrization.h>

#include <Math/GenVector/PositionVector3D.h>
#include <TAxis.h>
#include <TEfficiency.h>
#include <TH2.h>
#include <THashList.h>
#include <TPDGCode.h>
#include <TRandom3.h>
#include <TString.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <math.h>

using namespace o2;
using namespace o2::framework;

struct OnTheFlyMidPid {
  Produces<aod::UpgradeMidMC> upgradeMidMC;
  Produces<aod::UpgradeMid> upgradeMid;
  Produces<aod::UpgradeMidExpected> upgradeMidExpected;

  // necessary for particle charges
  Service<o2::framework::O2DatabasePDG> pdgDatabase;
  // Necessary for LUTs
  Service<o2::ccdb::BasicCCDBManager> ccdb;

  // these are the settings governing the MID layer to be used
  // in the evolution of this effort
  struct : ConfigurableGroup {
    Configurable<bool> considerEventTime{"considerEventTime", true, "flag to consider event time"};
    // Parametrization of muon/pion absortion in iron
    Configurable<float> magRadiusC{"magRadiusC", 160., "Magnet coil radius (cm)"};
    Configurable<float> absRadius1{"absRadius1", 220., "Absorber inner radius (cm)"};
    Configurable<float> absRadius2{"absRadius2", 290., "Absorber outer radius (cm)"};
    Configurable<float> absDensity{"absDensity", 7.874, "iron density (g/cm3)"};
    Configurable<float> absLambdaPi{"absLambdaPi", 20.41, "illustrative effective interaction length (cm)"};
    Configurable<float> absDeDx{"absDeDxMIP", 1.7, "dE/(rho dx) MIP MeV/g/cm2"};
    Configurable<float> absFakeAt3{"absFakeAt3", 0.032, "fake-muon efficiency at 3 GeV/c from GEANT4+BDT abs(eta)<1.24, 70 cm"};
    Configurable<float> absFakeAt10{"absFakeAt10", 0.08, "fake-muon efficiency at 10 GeV/c from GEANT4+BDT abs(eta)<1.24, 70 cm"};
    Configurable<float> absSigmaEta{"absSigmaEta", 10.65, "user-configurable phenomenological eta width"};
    Configurable<float> absPionFractionAt3{"absPionFractionAt3", 0.5, "parameter from GEANT4+BDT"}; // assumed, NOT measured
    Configurable<float> absEtaMax{"absEtaMax", 1.24, "max eta covered by detector"};
    Configurable<float> midRadius{"midRadius", 302.5, "MID radius (cm)"};
    Configurable<float> midLength{"midLength", 1000.0, "MID length (cm)"};
    Configurable<std::array<float, 2>> midCellDimension{"midCellDimension", {5.1, 5.1}, "MID cell dimension in Z and RPhi (cm)"};
    Configurable<float> midFractionOfInactiveArea{"midFractionOfInactiveArea", (0.1 / 5.1), "MID fraction of inactive area"};
    Configurable<float> midTimeReso{"midTimeReso", 1500, "MID time error (ps)"};
    Configurable<int> nStepsLIntegrator{"nStepsLIntegrator", 200, "number of steps in length integrator"};
    Configurable<float> multiplicityEtaRange{"multiplicityEtaRange", 0.800000012, "eta range to compute the multiplicity"};
    Configurable<bool> flagIncludeTrackTimeRes{"flagIncludeTrackTimeRes", true, "flag to include or exclude track time resolution"};
    Configurable<bool> flagMIDLoadDelphesLUTs{"flagMIDLoadDelphesLUTs", false, "flag to load Delphes LUTs for tracking correction (use recoTrack parameters if false)"};
  } simConfig;

  struct : ConfigurableGroup {
    Configurable<bool> doQAplots{"doQAplots", true, "do basic velocity plot qa"};
    Configurable<bool> doSeparationVsPt{"doSeparationVsPt", true, "Produce plots vs pt or p"};
    Configurable<int> nBinsBeta{"nBinsBeta", 2200, "number of bins in beta"};
    Configurable<int> nBinsP{"nBinsP", 80, "number of bins in momentum/pT depending on doSeparationVsPt"};
    Configurable<int> nBinsTrackLength{"nBinsTrackLength", 500, "number of bins in track length"};
    Configurable<int> nBinsNsigmaWrongSpecies{"nBinsNsigmaWrongSpecies", 200, "number of bins in Nsigma plot (wrong species)"};
    Configurable<float> minNsigmaRange{"minNsigmaRange", -10, "lower limit for the nsigma axis"};
    Configurable<float> maxNsigmaRange{"maxNsigmaRange", +10, "upper limit for the nsigma axis"};
    Configurable<int> nBinsDeltaCorrectSpecies{"nBinsDeltaCorrectSpecies", 200, "number of bins in Delta plot (correct speies)"};
    Configurable<int> nBinsDeltaWrongSpecies{"nBinsDeltaWrongSpecies", 200, "number of bins in Delta plot (wrong species)"};
    Configurable<int> nBinsTrackDeltaLength{"nBinsTrackDeltaLength", 100, "number of bins in delta track length"};
    Configurable<int> nBinsNsigmaCorrectSpecies{"nBinsNsigmaCorrectSpecies", 200, "number of bins in Nsigma plot (correct speies)"};
    Configurable<float> minDeltaRange{"minDeltaRange", -100, "lower limit for the nsigma axis"};
    Configurable<float> maxDeltaRange{"maxDeltaRange", +100, "upper limit for the nsigma axis"};
    Configurable<int> nBinsTimeRes{"nBinsTimeRes", 400, "number of bins plots time resolution"};
    Configurable<int> nBinsRelativeEtaPt{"nBinsRelativeEtaPt", 400, "number of bins plots pt and eta relative errors"};
    Configurable<int> nBinsEta{"nBinsEta", 400, "number of bins plot relative eta error"};
    Configurable<int> nBinsMult{"nBinsMult", 200, "number of bins in multiplicity"};
    Configurable<float> maxMultRange{"maxMultRange", 1000.f, "upper limit in multiplicity plots"};
  } plotsConfig;

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // Track smearer array, one per geometry
  // std::vector<std::unique_ptr<o2::delphes::TrackSmearer>> mSmearer;
  std::vector<std::unique_ptr<o2::fastsim::TrackSmearer>> mSmearer;
  // needed: random number generator for smearing
  TRandom3 pRandomNumberGenerator;

  // for handling basic QA histograms if requested
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<THashList> listEfficiency{"efficiency"};

  // Configuration defined at init time
  o2::fastsim::GeometryContainer mGeoContainer;
  float mMagneticField = 0.0f;
  void init(o2::framework::InitContext& initContext)
  {
    mGeoContainer.setCcdbManager(ccdb.operator->());
    mGeoContainer.init(initContext);

    const int nGeometries = mGeoContainer.getNumberOfConfigurations();
    mMagneticField = mGeoContainer.getFloatValue(0, "global", "magneticfield");

    pRandomNumberGenerator.SetSeed(0); // fully randomize
    if (simConfig.flagMIDLoadDelphesLUTs) {
      for (int icfg = 0; icfg < nGeometries; ++icfg) {
        mSmearer.emplace_back(std::make_unique<o2::fastsim::TrackSmearer>());
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

    if (plotsConfig.doQAplots) {
      const AxisSpec axisdNdeta{plotsConfig.nBinsMult, 0.0f, plotsConfig.maxMultRange, Form("dN/d#eta in |#eta| < %f", simConfig.multiplicityEtaRange.value)};

      histos.add("h1dNdeta", "h2dNdeta", kTH1F, {axisdNdeta});
      histos.add("h2dEventTime", "h2dEventTime", kTH2F, {{2000, -10000, 10000, "computed"}, {2000, -10000, 10000, "generated"}});
      histos.add("h1dEventTimegen", "h1dEventTimegen", kTH1F, {{200, -1000, 1000, "generated"}});
      histos.add("h1dEventTimerec", "h1dEventTimerec", kTH1F, {{2000, -10000, 10000, "computed"}});
      histos.add("h1dEventTimedelta", "h1dEventTimedelta", kTH1F, {{200, -1000, 1000, "generated - computed"}});
      histos.add("h2dEventTimeres", "h2dEventTimeres", kTH2F, {axisdNdeta, {300, 0, 300, "resolution"}});
      listEfficiency.setObject(new THashList);
      listEfficiency->Add(new TEfficiency("effEventTime", "effEventTime", plotsConfig.nBinsMult, 0.0f, plotsConfig.maxMultRange));

      const AxisSpec axisMomentum{static_cast<int>(plotsConfig.nBinsP), 0.0f, +10.0f, "#it{p} (GeV/#it{c})"};
      const AxisSpec axisRigidity{static_cast<int>(plotsConfig.nBinsP), 0.0f, +10.0f, "#it{p} / |#it{z}| (GeV/#it{c})"};
      const AxisSpec axisVelocity{static_cast<int>(plotsConfig.nBinsBeta), 0.0f, +1.1f, "Measured #beta"};
      const AxisSpec axisMomentumSmall{static_cast<int>(plotsConfig.nBinsP), 0.0f, +1.0f, "#it{p} (GeV/#it{c})"};
      const AxisSpec axisTrackLength{static_cast<int>(plotsConfig.nBinsTrackLength), 0.0f, 1000.0f, "Track length (cm)"};
      const AxisSpec axisTrackDeltaLength{static_cast<int>(plotsConfig.nBinsTrackDeltaLength), 0.0f, 30.0f, "Delta Track length (cm)"};

      histos.add("h2dVelocityVsMomentumInner", "h2dVelocityVsMomentumInner", kTH2F, {axisMomentum, axisVelocity});
      histos.add("h2dVelocityVsRigidityInner", "h2dVelocityVsRigidityInner", kTH2F, {axisRigidity, axisVelocity});
      histos.add("h2dTrackLengthInnerVsPt", "h2dTrackLengthInnerVsPt", kTH2F, {axisMomentum, axisTrackLength});
      histos.add("h2dTrackLengthInnerRecoVsPt", "h2dTrackLengthInnerRecoVsPt", kTH2F, {axisMomentum, axisTrackLength});

      histos.add("h2dTrackLengthB2vsB0", "h2dDeltaTrackLengthInnerVsPt; for B=0; for B=2 T", kTH2F, {axisTrackLength, axisTrackLength});
      histos.add("h2HitMap", "h2HitMap", kTH2F, {{200, -simConfig.midLength / 2, simConfig.midLength / 2}, {1000, 0, simConfig.midRadius * 2 * M_PI}});
      histos.add("h2HitMapRec", "h2HitMapRec", kTH2F, {{200, -simConfig.midLength / 2, simConfig.midLength / 2}, {1000, 0, simConfig.midRadius * 2 * M_PI}});
      histos.add("h2HitDeltaXY", "h2HitDeltaXY", kTH2F, {{200, -100, 100, "#Delta#it{x} (cm)"}, {200, -100, 100, "#Delta#it{y} (cm)"}});

      const AxisSpec axisPt{static_cast<int>(plotsConfig.nBinsP), 0.0f, +10.0f, "#it{p}_{T} (GeV/#it{c})"};
      const AxisSpec axisPtRes{1000, -1.0f, +1.0f, "(#it{p}_{T}^{gen}-#it{p}_{T}^{rec})/#it{p}_{T}^{gen}"};
      const AxisSpec axisEta{static_cast<int>(plotsConfig.nBinsEta), -2.0f, +2.0f, "#eta"};
      const AxisSpec axisRelativePt{static_cast<int>(plotsConfig.nBinsRelativeEtaPt), 0.0f, +10.0f, "#it{#sigma_{p_{T}}} / #it{p_{T}} (%)"};
      const AxisSpec axisRelativeEta{static_cast<int>(plotsConfig.nBinsRelativeEtaPt), 0.0f, +10.0f, "#it{#sigma_{#eta}} / #it{#eta} (%)"};
      histos.add("h2dEtaVsPt", "h2dTrackLengthInnerRecoVsPt", kTH2F, {axisMomentum, axisEta});
      histos.add("h2dEtaVsP", "h2dTrackLengthInnerRecoVsP", kTH2F, {axisMomentum, axisEta});

      const AxisSpec axisDelta{plotsConfig.nBinsDeltaCorrectSpecies, plotsConfig.minDeltaRange, plotsConfig.maxDeltaRange, "#Delta#{t} (ps)"};
      const AxisSpec axisNsigmaCorrect{plotsConfig.nBinsNsigmaCorrectSpecies, plotsConfig.minNsigmaRange, plotsConfig.maxNsigmaRange, "N#sigma - True "};
      const AxisSpec axisDeltaWrong{plotsConfig.nBinsDeltaWrongSpecies, plotsConfig.minDeltaRange, plotsConfig.maxDeltaRange, "#Delta - True "};
      const AxisSpec axisX{plotsConfig.doSeparationVsPt.value ? axisPt : axisMomentum};
      // const AxisSpec axisNSigma{axisNsigmaCorrect};

      const AxisSpec axisTrackTimeRes{plotsConfig.nBinsTimeRes, 0.0f, +3000.0f, "Track time resolution (ps)"};
      const AxisSpec axisTotalTimeRes{plotsConfig.nBinsTimeRes, 0.0f, +3000.0f, "Total time resolution (ps)"};

      histos.add("h1PtRes", "", kTH2F, {axisMomentum, axisPtRes});
      histos.add("h2EkinVsEloss", "", kTH2F, {{200, 0.0f, 2000.f, "#it{E}_{loss} (MeV)"}, {200, 0.0f, 20000.f, "KE (MeV)"}});
    }
  }

  struct MIDEfficiency {
   private:
    const float layerRadius;
    const float layerLength;
    const float cellDimensionZ;
    const float cellDimensionRPhi;
    const float fractionInactive;
    const float magField;

    TAxis* axisZ = nullptr;
    TAxis* axisRPhi = nullptr;
    TAxis* axisInCellZ = nullptr;
    TAxis* axisInCellRPhi = nullptr;

    TH2F* hHitMapInCell = nullptr;
    TH2F* hHitMapInCellBefore = nullptr;
    TH2F* hHitMap = nullptr;

   public:
    ~MIDEfficiency()
    {
      if (0) {
        hHitMap->SaveAs(Form("/tmp/%s.png", hHitMap->GetName()));
        hHitMapInCell->SaveAs(Form("/tmp/%s.png", hHitMapInCell->GetName()));
        hHitMapInCellBefore->SaveAs(Form("/tmp/%s.png", hHitMapInCellBefore->GetName()));
      }

      delete axisZ;
      delete axisRPhi;
      delete axisInCellZ;
      delete axisInCellRPhi;
      delete hHitMap;
      delete hHitMapInCell;
      delete hHitMapInCellBefore;
    }

    MIDEfficiency(const MIDEfficiency&) = delete;
    MIDEfficiency& operator=(const MIDEfficiency&) = delete;

    MIDEfficiency(float r, float l, std::array<float, 2> pDimensions, float fIA, float m)
      : layerRadius(r), layerLength(l), cellDimensionZ(pDimensions[0]), cellDimensionRPhi(pDimensions[1]), fractionInactive(fIA), magField(m)
    {
      // Assuming square cells for simplicity
      const float circumference = 2.0f * M_PI * layerRadius;
      axisZ = new TAxis(static_cast<int>(layerLength / cellDimensionZ), -layerLength / 2, layerLength);
      axisRPhi = new TAxis(static_cast<int>(circumference / cellDimensionRPhi), 0.f, circumference);

      const float inactiveBorderRPhi = cellDimensionRPhi * std::sqrt(fractionInactive) / 2;
      const float inactiveBorderZ = cellDimensionZ * std::sqrt(fractionInactive) / 2;
      const double arrayRPhi[4] = {-cellDimensionRPhi / 2, -cellDimensionRPhi / 2 + inactiveBorderRPhi, cellDimensionRPhi / 2 - inactiveBorderRPhi, cellDimensionRPhi / 2};
      for (int i = 0; i < 4; i++) {
        LOG(info) << "arrayRPhi[" << i << "] = " << arrayRPhi[i];
      }
      axisInCellRPhi = new TAxis(3, arrayRPhi);
      const double arrayZ[4] = {-cellDimensionZ / 2, -cellDimensionZ / 2 + inactiveBorderZ, cellDimensionZ / 2 - inactiveBorderZ, cellDimensionZ / 2};
      for (int i = 0; i < 4; i++) {
        LOG(info) << "arrayZ[" << i << "] = " << arrayZ[i];
      }
      axisInCellZ = new TAxis(3, arrayZ);

      hHitMap = new TH2F(Form("hHitMap_R%.0f", layerRadius), "HitMap;z (cm); r#phi (cm)", 1000, -1000, 1000, 1000, -1000, 1000);
      hHitMapInCell = new TH2F(Form("hHitMapInCell_R%.0f", layerRadius), "HitMapInCell;z (cm); r#phi (cm)", 1000, -10, 10, 1000, -10, 10);
      hHitMapInCellBefore = new TH2F(Form("hHitMapInCellBefore_R%.0f", layerRadius), "HitMapInCell;z (cm); r#phi (cm)", 1000, -10, 10, 1000, -10, 10);
    }

    bool isInMIDActiveArea(std::array<float, 3> hitPosition)
    {
      if (fractionInactive <= 0.0f) {
        return true;
      }
      if (fractionInactive >= 1.0f) {
        return false;
      }

      // Convert 3D position to cylindrical coordinates for area calculation
      const float phi = std::atan2(hitPosition[1], hitPosition[0]);
      const float rphi = phi * layerRadius;
      const float z = hitPosition[2];
      const float r = std::sqrt(hitPosition[0] * hitPosition[0] + hitPosition[1] * hitPosition[1]);

      // Check if hit is within layer geometric acceptance
      if (std::abs(layerRadius - r) > 10.f) {
        LOG(debug) << "Hit out of MID layer acceptance: r=" << r << " cm with respect to the layer radius " << layerRadius;
        return false;
      }
      if (std::abs(z) > layerLength / 2.0f) {
        LOG(debug) << "Hit out of MID layer acceptance: z=" << z << " cm with respect to the layer length " << layerLength;
        return false;
      }

      const int cellIndexPhi = axisRPhi->FindBin(rphi);
      const int cellIndexZ = axisZ->FindBin(z);

      if (cellIndexPhi <= 0 || cellIndexPhi > axisRPhi->GetNbins() || cellIndexZ <= 0 || cellIndexZ > axisZ->GetNbins()) {
        return false;
      }
      // Calculate local position within the pixel
      const float localRPhi = (rphi - axisRPhi->GetBinCenter(cellIndexPhi));
      const float localZ = (z - axisZ->GetBinCenter(cellIndexZ));

      // The difference between the hit and the pixel position cannot be greater than the size of the pixel
      if (std::abs(localRPhi - axisRPhi->GetBinCenter(cellIndexPhi)) > axisRPhi->GetBinWidth(cellIndexPhi)) {
        // LOG(warning) << "Local hit difference in phi is bigger than the pixel size";
      }
      if (std::abs(localZ - axisZ->GetBinCenter(cellIndexZ)) > axisZ->GetBinWidth(cellIndexZ)) {
        // LOG(warning) << "Local hit difference in z is bigger than the pixel size";
      }
      hHitMapInCellBefore->Fill(localZ, localRPhi);
      switch (axisInCellRPhi->FindBin(localRPhi)) {
        case 0:
        case 1:
        case 3:
        case 4:
          return false;
      }
      switch (axisInCellZ->FindBin(localZ)) {
        case 0:
        case 1:
        case 3:
        case 4:
          return false;
      }
      hHitMapInCell->Fill(localZ, localRPhi);
      hHitMap->Fill(z, rphi);
      return true;
    }

    /// Check if a track hits the active area (convenience function)
    /// \param track the track parameters for automatic hit position calculation
    /// \return true if the hit is in the active area
    bool isTrackInActiveArea(const o2::track::TrackParCov& track)
    {
      if (fractionInactive <= 0.f) {
        return true;
      }
      float x, y, z;
      if (!track.getXatLabR(layerRadius, x, magField)) {
        LOG(debug) << "Could not propagate track to TOF layer at radius " << layerRadius << " cm";
        return false;
      }
      bool b;
      ROOT::Math::PositionVector3D hit = track.getXYZGloAt(x, magField, b);
      if (!b) {
        LOG(debug) << "Could not get hit position at radius " << layerRadius << " cm";
        return false;
      }
      hit.GetCoordinates(x, y, z);
      return isInMIDActiveArea(std::array<float, 3>{x, y, z});
    }
  };

  bool isInMIDActiveArea(const o2::track::TrackParCov& track)
  {
    static MIDEfficiency midLayer(simConfig.midRadius, simConfig.midLength, simConfig.midCellDimension, simConfig.midFractionOfInactiveArea, mMagneticField);
    return midLayer.isTrackInActiveArea(track);
  }

  struct TracksWithTime {
    TracksWithTime(int pdgCode,
                   std::pair<float, float> midTime,
                   std::pair<float, float> trackLengthMID,
                   std::pair<float, float> trackLengthInIron,
                   std::pair<std::array<float, 3>, std::array<float, 3>> posInMid,
                   std::pair<float, float> momentum,
                   std::pair<float, float> pseudorapidity,
                   float noSmearingPt)
      : mPdgCode(pdgCode),
        mMIDTime(midTime),
        mTrackLengthMID(trackLengthMID),
        mTrackLengthInIron(trackLengthInIron),
        mPosInMid(posInMid),
        mMomentum(momentum),
        mPseudorapidity{pseudorapidity},
        mNoSmearingPt{noSmearingPt} {}

    int mPdgCode;
    std::pair<float, float> mMIDTime;           // Measured time and expected resolution at the MID [ps]
    std::pair<float, float> mTrackLengthMID;    // Measured and generated length at the MID [cm]
    std::pair<float, float> mTrackLengthInIron; // Measured and generated length within absorber [cm]
    std::pair<std::array<float, 3>, std::array<float, 3>> mPosInMid;
    std::pair<float, float> mMomentum;       // Measured momentum and uncertainty on the momentum [GeV/c]
    std::pair<float, float> mPseudorapidity; // Measured pseudorapidity and uncertainty on the pseudorapidity
    float mNoSmearingPt;                     // No smearing pt
  };

  std::vector<TracksWithTime> tracksWithTime;
  bool eventTime(const std::vector<TracksWithTime>& tracks,
                 std::array<float, 2>& tzero)
  {

    float sum = 0.;
    float sumw = 0.;

    // Todo: check the different mass hypothesis iteratively
    for (const auto& track : tracks) {
      auto pdgInfo = pdgDatabase->GetParticle(track.mPdgCode);
      if (pdgInfo == nullptr) {
        continue;
      }
      const float mass = pdgInfo->Mass();
      const float mass2 = mass * mass;
      const float tof = track.mMIDTime.first; // [ps]
      if (tof <= 0) {
        continue;
      }
      const float etof = track.mMIDTime.second;         // [ps]
      const float length = track.mTrackLengthMID.first; // [cm]
      float p = track.mMomentum.first;                  // [GeV/c]
      p *= std::abs(pdgInfo->Charge()) / 3.;            // Total momentum
      const float ep = track.mMomentum.second;          // [GeV/c]
      const float p2 = p * p;
      const float lengthOverC = length * o2::constants::physics::invLightSpeedCm2PS;
      const float texp = lengthOverC / p * std::sqrt(mass2 + p2);
      // LOG(info) << "TOF: " << tof << " " << etof << " vs " << texp;
      const float etexp = lengthOverC * mass2 / p2 / std::sqrt(mass2 + p2) * ep;
      const float sigma = std::sqrt(etexp * etexp + etof * etof);
      const float deltat = tof - texp;

      const float w = 1. / (sigma * sigma);

      sum += w * deltat;
      sumw += w;
    }

    static constexpr float kMaxEventTimeResolution = 200.f;
    if (sumw <= 0. || tracks.size() <= 1 || std::sqrt(1. / sumw) > kMaxEventTimeResolution) {
      tzero[0] = 0.;                      // [ps]
      tzero[1] = kMaxEventTimeResolution; // [ps]
      return false;
    }

    tzero[0] = sum / sumw;
    tzero[1] = std::sqrt(1. / sumw);
    return true;
  }

  /// returns track time resolution
  /// \param pt the transverse momentum of the tarck
  /// \param eta the pseudorapidity of the tarck
  /// \param trackPtResolution the absolute resolution on pt
  /// \param trackEtaResolution the absolute resolution on eta
  /// \param mass the mass of the particle
  /// \param detRadius the radius of the cylindrical layer
  /// \param magneticField the magnetic field (along Z)
  double calculateTrackTimeResolutionAdvanced(const float pt,
                                              const float eta,
                                              const float trackPtResolution,
                                              const float trackEtaResolution,
                                              const float mass,
                                              const float detRadius,
                                              const float magneticField)
  {
    // Compute tracking contribution to timing using the error propagation formula
    // Uses light speed in m/ps, magnetic field in T (*0.1 for conversion kGauss -> T)
    float a0 = mass * mass;
    float a1 = 0.299792458 * (0.1 * magneticField) * (0.01 * o2::constants::physics::LightSpeedCm2NS / 1e+3);
    float a2 = (detRadius * 0.01) * (detRadius * 0.01) * (0.299792458) * (0.299792458) * (0.1 * magneticField) * (0.1 * magneticField) / 2.0;
    float dtofOndPt = (std::pow(pt, 4) * std::pow(std::cosh(eta), 2) * std::acos(1.0 - a2 / std::pow(pt, 2)) - 2.0 * a2 * std::pow(pt, 2) * (a0 + std::pow(pt * std::cosh(eta), 2)) / std::sqrt(a2 * (2.0 * std::pow(pt, 2) - a2))) / (a1 * std::pow(pt, 3) * std::sqrt(a0 + std::pow(pt * std::cosh(eta), 2)));
    float dtofOndEta = std::pow(pt, 2) * std::sinh(eta) * std::cosh(eta) * std::acos(1.0 - a2 / std::pow(pt, 2)) / (a1 * std::sqrt(a0 + std::pow(pt * std::cosh(eta), 2)));
    float trackTimeResolution = std::hypot(std::fabs(dtofOndPt) * trackPtResolution, std::fabs(dtofOndEta) * trackEtaResolution);
    return trackTimeResolution;
  }

  /*
  // ============================================================
  // Track length in a B = 0 region
  //
  // Calculates the distance from the CURRENT track position
  // to the cylinder of radius "radius".
  // ============================================================
   */
  float computeTrackLengthB0(const o2::track::TrackParCov& track,
                             const float radius)
  {
    // Current global position
    std::array<float, 3> pos{};
    track.getXYZGlo(pos);

    // Current global momentum
    std::array<float, 3> mom{};
    if (!track.getPxPyPzGlo(mom)) {
      return -100.f;
    }

    const float px = mom[0];
    const float py = mom[1];
    // const float pz = mom[2];

    const float pt = std::hypot(px, py);

    if (pt <= 0.f) {
      return -100.f;
    }

    // Unit vector in transverse plane
    const float ux = px / pt;
    const float uy = py / pt;

    const float b = pos[0] * ux + pos[1] * uy;
    const float c = pos[0] * pos[0] +
                    pos[1] * pos[1] -
                    radius * radius;

    const float discriminant = b * b - c;

    if (discriminant < 0.f) {
      return -100.f;
    }

    const float sqrtDisc = std::sqrt(discriminant);

    const float t1 = -b + sqrtDisc;
    const float t2 = -b - sqrtDisc;

    // We want the outward/forward intersection
    float t = -1.f;

    if (t1 > 0.f && t2 > 0.f) {
      t = std::min(t1, t2);
    } else if (t1 > 0.f) {
      t = t1;
    } else if (t2 > 0.f) {
      t = t2;
    } else {
      return -100.f;
    }

    return t * std::sqrt(1.f + track.getTgl() * track.getTgl());
  }

  /*
  // ============================================================
  // MID piecewise track length
  //
  // B = +20.0 kG   R < R0
  // B =   0.0 kG   R0 < R < R1
  // B = -14.3 kG   R1 < R < R2
  //
  // Returns total track length from current track position to R3.
  // ============================================================
   */
  float computeTrackLengthMID(const o2::track::TrackParCov& inputTrack,
                              float R0,
                              float R1,
                              float R2,
                              float R3,
                              std::array<float, 3>& posMID,
                              std::array<float, 3>& deltaMID,
                              float& deltaposInIron)
  {
    float B0_m = +20.0f; // kG
    float B1_m = 0.0f;   // kG
    float B2_m = 0.0f;
    B2_m = -B0_m * pow(R0, 2) / (pow(R2, 2) - pow(R1, 2)); // kG
    float B3_m = 0.0f;                                     // kG

    // Work with a copy
    auto track = inputTrack;

    float totalLength = 0.f;
    float xAtR = 0.f;

    // ==========================================================
    // REGION 1
    // current position -> R0
    // B = +20 kG
    // ==========================================================

    const float L0 =
      o2::upgrade::computeTrackLength(track, R0, B0_m);

    if (L0 < 0.f) {
      return -100.f;
    }

    totalLength += L0;

    // Find local x corresponding to global radius R0
    if (!track.getXatLabR(R0,
                          xAtR,
                          B0_m,
                          o2::track::DirOutward)) {
      return -100.f;
    }

    // Propagate track to R0
    if (!track.propagateTo(xAtR, B0_m)) {
      return -100.f;
    }

    // ==========================================================
    // REGION 2
    // R0 -> R1
    // B = 0
    // ==========================================================

    const float L01 = computeTrackLengthB0(track, R1);

    if (L01 < 0.f) {
      return -100.f;
    }

    totalLength += L01;

    // Find local x at R1 assuming B = 0
    if (!track.getXatLabR(R1,
                          xAtR,
                          B1_m,
                          o2::track::DirOutward)) {
      return -100.f;
    }

    // Straight-line propagation to R1
    if (!track.propagateTo(xAtR, B1_m)) {
      return -100.f;
    }

    // ==========================================================
    // REGION 3 (absorber)
    // R1 -> R2
    // B = -14.3 kG
    // ==========================================================

    std::array<float, 3> pos1{};
    track.getXYZGlo(pos1);

    float L12 =
      o2::upgrade::computeTrackLength(track, R2, B2_m);

    if (L12 < 0.f) {
      return -100.f;
    }
    totalLength += L12;
    deltaposInIron = L12;
    // Find local x corresponding to global radius R2
    if (!track.getXatLabR(R2,
                          xAtR,
                          B2_m,
                          o2::track::DirOutward)) {
      return -100.f;
    }

    // Propagate track to R2
    if (!track.propagateTo(xAtR, B2_m)) {
      return -100.f;
    }

    // ==========================================================
    // REGION 2
    // R2 -> R3
    // B = 0
    // ==========================================================

    const float L23 = computeTrackLengthB0(track, R3);

    if (L23 < 0.f) {
      return -100.f;
    }

    totalLength += L23;

    // Find local x corresponding to global radius R3
    if (!track.getXatLabR(R3,
                          xAtR,
                          B3_m,
                          o2::track::DirOutward)) {
      return -100.f;
    }

    // Propagate track to R3
    if (!track.propagateTo(xAtR, B3_m)) {
      return -100.f;
    }

    // ==========================================================
    // GLOBAL POSITION AT R3
    // ==========================================================

    std::array<float, 3> pos3{};
    track.getXYZGlo(pos3);

    posMID[0] = pos3[0];
    posMID[1] = pos3[1];
    posMID[2] = pos3[2];

    // Displacement projections
    deltaMID[0] = posMID[0] - pos1[0];
    deltaMID[1] = posMID[1] - pos1[1];
    deltaMID[2] = posMID[2] - pos1[2];

    return totalLength;
  }
  // Multiple Coulomb scattering angle in iron
  //
  // momentum  : particle momentum [GeV/c]
  // thickness : iron thickness [cm]
  //
  // returns theta0 [rad]

  float moliereThetaIron(float momentum, float mass, float thickness)
  {
    float X0 = 1.757f; // radiation length of iron [cm]
                       //  constexpr float mMuon = 0.105658f; // muon mass [GeV/c^2]

    if (momentum <= 0.f || thickness <= 0.f) {
      return 0.f;
    }

    // beta = p/E
    const float energy = std::sqrt(momentum * momentum +
                                   mass * mass);

    const float beta = momentum / energy;

    // Thickness in radiation lengths
    const float xOverX0 = thickness / X0;

    // Highland approximation
    const float theta0 = (0.0136f / (beta * momentum)) * sqrt(xOverX0) * (1.f + 0.038f * log(xOverX0));

    return theta0;
  }
  float fPi(double p)
  {
    return p > 1. ? 1. - std::exp(-(p - 1.) / 1.5) : 0.;
  }
  float fMu(double p)
  {
    return p > 0.7 ? (1. - std::exp(-(p - 0.7) / 0.8)) * std::exp(-p / 3.) : 0.;
  }
  float refPi(double p)
  {
    if (p <= 1.)
      return 0.;
    float f3 = fPi(3.);
    float mu10 = simConfig.absFakeAt3 * (1. - simConfig.absPionFractionAt3) * fMu(10.) / fMu(3.);
    float slope = ((simConfig.absFakeAt10 - mu10) / (simConfig.absFakeAt3 * simConfig.absPionFractionAt3) / (fPi(10.) / f3) - 1.) / 7.;
    return std::max<double>(0., simConfig.absFakeAt3 * simConfig.absPionFractionAt3) * fPi(p) / f3 * (1. + slope * (p - 3.));
  }
  float refMu(double p)
  {
    return simConfig.absFakeAt3 * (1. - simConfig.absPionFractionAt3) * fMu(p) / fMu(3.);
  }
  float decayAir(double p, double airCm)
  {
    if (p <= 0.)
      return 0.;
    return -std::expm1(-airCm * 0.13957039 / (p * 780.45)); // assume pion mass, ctau=780.45 cm
  }
  float survivalAir(double p, double airCm)
  {
    return 1. - decayAir(p, airCm);
  }
  // Range proxy: kinetic energy minus mean ionization loss in iron.
  // Smoothed only to represent missing fluctuations; NOT a Geant4 range calculation.
  float rangeTurnOn(double p, double ironCm, double mass)
  {
    if (p <= 0.)
      return 0.;
    float kinetic = std::sqrt(p * p + mass * mass) - mass;
    float deficit = kinetic - (simConfig.absDensity * simConfig.absDeDx / 1000.) * ironCm;
    float width = 0.10; // GeV, phenomenological
    return 1. / (1. + std::exp(-std::max<double>(-60., std::min<double>(60., deficit / width))));
  }
  float etaShape(double eta)
  {
    return std::exp(-0.5 * eta * eta / (simConfig.absSigmaEta * simConfig.absSigmaEta));
  }
  float pion(float p, float eta, float airCm, float ironCm)
  {
    if (p <= 0.)
      return 0.;
    const float c = std::cosh(eta);
    const float survivalRatio = survivalAir(p, airCm * c) / survivalAir(p, simConfig.absRadius1 * c);
    const float hadronRatio = std::exp(-(ironCm - (simConfig.absRadius2 - simConfig.absRadius1)) * c / simConfig.absLambdaPi);
    // Ratio to reference range at the SAME eta; reference shape is phenomenological.
    const float rangeRatio = rangeTurnOn(p, ironCm * c, 0.13957039) / rangeTurnOn(p, (simConfig.absRadius2 - simConfig.absRadius1) * c, 0.13957039); // assume pion mass
    return std::max<double>(0., refPi(p) * etaShape(eta) * survivalRatio * hadronRatio * rangeRatio);
  }
  float muon(float p, float eta, float airCm, float ironCm)
  {
    if (p <= 0.)
      return 0.;
    const float c = std::cosh(eta);
    const float d0 = decayAir(p, simConfig.absRadius1 * c);
    const float decayRatio = d0 > 0. ? decayAir(p, airCm * c) / d0 : 0.;
    // Effective secondary-muon momentum proxy: 0.8 * pion momentum.
    // Secondary muons born in iron are NOT modeled; this is a phenomenological air-decay term.
    const float pMu = 0.8 * p;
    const float rangeRatio = rangeTurnOn(pMu, ironCm * c, 0.1056583755) /
                             rangeTurnOn(pMu, (simConfig.absRadius2 - simConfig.absRadius1) * c, 0.1056583755);
    return std::max<double>(0., refMu(p) * etaShape(eta) * decayRatio * rangeRatio);
  }
  float total(float p, float eta, float airCm, float ironCm)
  {
    return std::min<double>(1., pion(p, eta, airCm, ironCm) + muon(p, eta, airCm, ironCm));
  }
  float meanEta(float p, float airCm, float ironCm)
  {
    const int n = 1000;
    float s = 0.;
    for (int i = 0; i < n; ++i) {
      float eta = -simConfig.absEtaMax + (i + .5) * 2. * simConfig.absEtaMax / n;
      s += total(p, eta, airCm, ironCm);
    }
    return s / n;
  }

  void process(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::OTFLUTConfigId>::iterator const& collision,
               soa::Join<aod::Tracks, aod::TracksCov, aod::McTrackLabels> const& tracks,
               aod::McParticles const&,
               aod::McCollisions const&)
  {
    o2::dataformats::VertexBase pvVtx({collision.posX(), collision.posY(), collision.posZ()},
                                      {collision.covXX(), collision.covXY(), collision.covYY(),
                                       collision.covXZ(), collision.covYZ(), collision.covZZ()});

    std::array<float, 6> mcPvCov = {0.};
    o2::dataformats::VertexBase mcPvVtx({0.0f, 0.0f, 0.0f}, mcPvCov);
    const float eventCollisionTimePS = (simConfig.considerEventTime.value ? collision.collisionTime() * 1e3 : 0.f); // convert ns to ps
    if (collision.has_mcCollision()) {
      auto mcCollision = collision.mcCollision();
      mcPvVtx.setX(mcCollision.posX());
      mcPvVtx.setY(mcCollision.posY());
      mcPvVtx.setZ(mcCollision.posZ());
    } // else remains untreated for now

    // First we compute the number of charged particles in the event if LUTs are loaded
    float dNdEta = 0.f;
    for (const auto& track : tracks) {
      if (!track.has_mcParticle())
        continue;
      auto mcParticle = track.mcParticle();
      if (std::abs(mcParticle.eta()) > simConfig.multiplicityEtaRange) {
        continue;
      }
      if (mcParticle.has_daughters()) {
        continue;
      }
      const auto& pdgInfo = pdgDatabase->GetParticle(mcParticle.pdgCode());
      if (!pdgInfo) {
        // LOG(warning) << "PDG code " << mcParticle.pdgCode() << " not found in the database";
        continue;
      }
      if (pdgInfo->Charge() == 0) {
        continue;
      }
      dNdEta += 1.f;
    }
    if (plotsConfig.doQAplots) {
      histos.fill(HIST("h1dNdeta"), dNdEta);
    }

    tracksWithTime.clear(); // clear the vector of tracks with time to prepare the cache for the next event
    tracksWithTime.reserve(tracks.size());
    // First loop to generate the arrival time of the particles to the MID layers
    for (const auto& track : tracks) {
      // first step: find precise arrival time (if any)
      // --- convert track into perfect track
      if (!track.has_mcParticle()) { // should always be OK but check please
        upgradeMidMC(-999.f, -999.f);
        continue;
      } else {
        LOG(debug) << "Track without mcParticle found!";
      }
      const auto& mcParticle = track.mcParticle();
      float true_pt = mcParticle.pt();
      o2::track::TrackParCov o2track = o2::upgrade::convertMCParticleToO2Track(mcParticle, pdgDatabase);
      // get mass to calculate velocity
      auto pdgInfo = pdgDatabase->GetParticle(mcParticle.pdgCode());
      if (pdgInfo == nullptr) {
        LOG(error) << "PDG code " << mcParticle.pdgCode() << " not found in the database";
        upgradeMidMC(-999.f, -999.f);
        continue;
      }

      float xPv = -100.f;
      static constexpr float kTrkXThreshold = -99.f; // Threshold to consider a good propagation of the track
      if (o2track.propagateToDCA(mcPvVtx, mMagneticField)) {
        xPv = o2track.getX();
      }
      float trackLengthMID = -1;
      float trackLengthMIDB2 = -1;
      std::array<float, 3> posMid{-999., -999., -999.};
      std::array<float, 3> delaPosMid{-999., -999., -999.};
      float trackLengthInIron = -999.;
      if (xPv > kTrkXThreshold) {
        trackLengthMIDB2 = o2::upgrade::computeTrackLength(o2track, simConfig.magRadiusC, mMagneticField);
        trackLengthMID = computeTrackLengthMID(o2track, simConfig.magRadiusC, simConfig.absRadius1, simConfig.absRadius2, simConfig.midRadius, posMid, delaPosMid, trackLengthInIron);
        // thetamolInIron = moliereThetaIron(o2track.getP(), pdgInfo->Mass(), trackLengthInIron);
      }
      if (plotsConfig.doQAplots) {
        histos.fill(HIST("h2dTrackLengthB2vsB0"), trackLengthMIDB2, trackLengthMID);
      }

      // Check if the track hit a sensitive area of the MID
      const bool activeMID = isInMIDActiveArea(o2track);
      if (!activeMID) {
        trackLengthMID = -999.f;
      } else {
        if (plotsConfig.doQAplots) {
          float phi_tmp = std::atan2(posMid[1], posMid[0]);
          if (phi_tmp < 0.f) {
            phi_tmp += 2.f * M_PI;
          }
          histos.fill(HIST("h2HitMap"), posMid[2], simConfig.midRadius * phi_tmp);
        }
      }

      const float v = o2::upgrade::computeParticleVelocity(o2track.getP(), pdgInfo->Mass());
      const float expectedTimeMID = trackLengthMID > 0 ? trackLengthMID / v + eventCollisionTimePS : -999.f; // arrival time to the MID in ps

      upgradeMidMC(expectedTimeMID, trackLengthMID);

      // Smear with expected resolutions
      const float measuredTimeMID = pRandomNumberGenerator.Gaus(expectedTimeMID, simConfig.midTimeReso);

      // Now we calculate the expected arrival time following certain mass hypotheses
      // and the (imperfect!) reconstructed track parametrizations
      float trackLengthRecoMID = -1;
      auto recoTrack = getTrackParCov(track);
      xPv = -100.f;
      std::array<float, 3> posRecMid{-999., -999., -999.};
      std::array<float, 3> deltaPosRecMid{-999., -999., -999.};
      float trackLengthRecoInIron = -1.;
      trackLengthRecoInIron = -1.;
      trackLengthRecoMID = -1.;
      if (recoTrack.propagateToDCA(pvVtx, mMagneticField)) {
        xPv = recoTrack.getX();
      }
      if (xPv > kTrkXThreshold) {
        trackLengthRecoMID = computeTrackLengthMID(recoTrack, simConfig.magRadiusC, simConfig.absRadius1, simConfig.absRadius2, simConfig.midRadius, posRecMid, deltaPosRecMid, trackLengthRecoInIron);
        //  mimic scattering in iron
        std::array<float, 3> mom{-999., -999., -999.};
        recoTrack.getPxPyPzGlo(mom);
        for (int jj = 0; jj < 3; ++jj) {
          float posRectmp = 0.0f;
          float molieretmp = 0.0f;
          molieretmp = moliereThetaIron(mom[jj], pdgInfo->Mass(), deltaPosRecMid[jj]);
          posRectmp = pRandomNumberGenerator.Gaus(posRecMid[jj], molieretmp);
          posRecMid[jj] = posRectmp;
        }
        if (plotsConfig.doQAplots) {
          float phi_tmp = std::atan2(posRecMid[1], posRecMid[0]);
          if (phi_tmp < 0.f) {
            phi_tmp += 2.f * M_PI;
          }
          histos.fill(HIST("h2HitMapRec"), posRecMid[2], simConfig.midRadius * phi_tmp);
          histos.fill(HIST("h2HitDeltaXY"), posRecMid[0] - posMid[0], posRecMid[1] - posMid[1]);
          histos.fill(HIST("h1PtRes"), true_pt, (true_pt - recoTrack.getPt()) / true_pt);
        }
      }
      //  cache the track info needed for the event time calculation
      //  Reuse or emplace a new object in the vector
      tracksWithTime.emplace_back(TracksWithTime{mcParticle.pdgCode(),
                                                 {measuredTimeMID, expectedTimeMID},
                                                 {trackLengthRecoMID, trackLengthMID},
                                                 {trackLengthRecoInIron, trackLengthInIron},
                                                 {posRecMid, posMid},
                                                 {recoTrack.getP(), recoTrack.getSigma1Pt2()},
                                                 {recoTrack.getEta(), recoTrack.getSigmaTgl2()},
                                                 o2track.getPt()});
    }
    //  Now we compute the event time for the tracks

    std::array<float, 2> tzero = {0.f, 0.f};
    bool etStatus = false;
    if (simConfig.considerEventTime.value) {
      etStatus = eventTime(tracksWithTime, tzero);
      if (!etStatus) {
        LOG(debug) << "Event time calculation failed with " << tracksWithTime.size() << " tracks with time and " << dNdEta << " charged particles";
      }
    }

    if (plotsConfig.doQAplots) {
      histos.fill(HIST("h2dEventTime"), tzero[0], eventCollisionTimePS);
      histos.fill(HIST("h1dEventTimegen"), eventCollisionTimePS);
      histos.fill(HIST("h1dEventTimerec"), tzero[0]);
      histos.fill(HIST("h2dEventTimeres"), dNdEta, tzero[1]);
      if (etStatus) {
        histos.fill(HIST("h1dEventTimedelta"), eventCollisionTimePS - tzero[0]);
      }
      static_cast<TEfficiency*>(listEfficiency->At(0))->Fill(etStatus, dNdEta);
    }

    // Then we do a second loop to compute the measured quantities with the measured event time
    int trackWithTimeIndex = 0;
    for (const auto& track : tracks) {
      if (!track.has_mcParticle()) { // should always be OK but check please
        continue;
      }
      const auto& mcParticle = track.mcParticle();

      const auto& trkWithTime = tracksWithTime[trackWithTimeIndex++];
      const float trackLengthRecoMID = trkWithTime.mTrackLengthMID.first;
      const float trackLengthMID = trkWithTime.mTrackLengthMID.second;
      const float measuredTimeMID = trkWithTime.mMIDTime.first;
      const float expecTimeMID = trkWithTime.mMIDTime.second;
      const float momentum = trkWithTime.mMomentum.first;
      const float pseudorapidity = trkWithTime.mPseudorapidity.first;
      const float noSmearingPt = trkWithTime.mNoSmearingPt;
      std::array<float, 3> expectedPosMid = trkWithTime.mPosInMid.second;
      std::array<float, 3> measuredPosMid = trkWithTime.mPosInMid.first;
      float trackLengthInIronRec = trkWithTime.mTrackLengthInIron.first;
      bool isMuCandidate = false;
      if (abs(pseudorapidity) < simConfig.absEtaMax) {
        if (trackLengthInIronRec > 0) {
          float kineticEnergy = std::sqrt(pow(momentum, 2) + pow(0.1057, 2)) - 0.1057;
          float energyloss = trackLengthInIronRec * simConfig.absDeDx * simConfig.absDensity / 1000.;
          if (abs(trkWithTime.mPdgCode) == 13) {
            if ((kineticEnergy - energyloss) > 0.1) {
              if (plotsConfig.doQAplots) {
                histos.fill(HIST("h2EkinVsEloss"), energyloss * 1000., kineticEnergy * 1000.);
              }
              isMuCandidate = true;
            }
          } else if (abs(trkWithTime.mPdgCode) == 11) {
            isMuCandidate = false;
          } else { // mostly hadrons, parametrization from simulations using pions, GEANT4+BDT
            float eff_rndm = pRandomNumberGenerator.Uniform(0., 1.);
            float eff_had = total(momentum, pseudorapidity, trackLengthRecoMID - trackLengthInIronRec, trackLengthInIronRec);
            if (eff_rndm < eff_had && (kineticEnergy - energyloss) > 0.1) {
              if (plotsConfig.doQAplots) {
                histos.fill(HIST("h2EkinVsEloss"), energyloss * 1000., kineticEnergy * 1000.);
              }
              isMuCandidate = true;
            }
          }
        }
      }
      auto truePdgInfo = pdgDatabase->GetParticle(mcParticle.pdgCode());
      float rigidity = momentum; // fallback to momentum if charge unknown

      // Use MC truth charge for rigidity calculation
      if (truePdgInfo) {
        const float trueCharge = std::abs(truePdgInfo->Charge()) / 3.0f;
        if (trueCharge > 0) {
          rigidity = momentum / trueCharge;
        }
      }

      if (plotsConfig.doQAplots) {
        // unit conversion: length in cm, time in ps
        const float midBeta = (trackLengthRecoMID / measuredTimeMID) / o2::constants::physics::LightSpeedCm2PS;
        if (trackLengthRecoMID > 0) {
          histos.fill(HIST("h2dVelocityVsMomentumInner"), momentum, midBeta);
          histos.fill(HIST("h2dVelocityVsRigidityInner"), rigidity, midBeta);
          histos.fill(HIST("h2dTrackLengthInnerVsPt"), noSmearingPt, trackLengthMID);
          histos.fill(HIST("h2dTrackLengthInnerRecoVsPt"), noSmearingPt, trackLengthRecoMID);
          histos.fill(HIST("h2dEtaVsPt"), noSmearingPt, pseudorapidity);
          histos.fill(HIST("h2dEtaVsP"), momentum, pseudorapidity);
        }
      }
      std::vector<float> measuredPosMidVec(measuredPosMid.begin(), measuredPosMid.end());
      upgradeMid(tzero[0], tzero[1], measuredPosMidVec, measuredTimeMID, trackLengthRecoMID, isMuCandidate);

      std::vector<float> expectedPosMidVec(expectedPosMid.begin(), expectedPosMid.end());
      upgradeMidExpected(expectedPosMidVec, expecTimeMID, trackLengthMID);
    }

    if (trackWithTimeIndex != tracks.size()) {
      LOG(fatal) << "Track with time index " << trackWithTimeIndex << " does not match the number of tracks " << tracks.size();
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<OnTheFlyMidPid>(cfgc)}; }
