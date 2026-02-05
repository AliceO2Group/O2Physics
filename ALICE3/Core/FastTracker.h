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

#ifndef ALICE3_CORE_FASTTRACKER_H_
#define ALICE3_CORE_FASTTRACKER_H_

#include "DetLayer.h"

#include <CCDB/BasicCCDBManager.h>
#include <Framework/InitContext.h>
#include <Framework/Logger.h>
#include <ReconstructionDataFormats/Track.h>

#include <map>
#include <string>
#include <vector>

namespace o2
{
namespace fastsim
{

class GeometryContainer
{
 public:
  GeometryContainer() = default;
  virtual ~GeometryContainer() = default;

  void init(o2::framework::InitContext& initContext);

  /**
   * @brief Parses a TEnv configuration file and returns the key-value pairs split per entry
   * @param filename Path to the TEnv configuration file
   * @param layers Vector to store the order of the layers as they appear in the file
   * @return A map where each key is a layer name and the value is another map of key-value pairs for that layer
   */
  static std::map<std::string, std::map<std::string, std::string>> parseTEnvConfiguration(std::string& filename, std::vector<std::string>& layers);

  // A container for the geometry info
  struct GeometryEntry {
    // Default constructor
    GeometryEntry() = default;
    explicit GeometryEntry(std::string filename)
    {
      mFileName = filename;
      mConfigurations = GeometryContainer::parseTEnvConfiguration(mFileName, mLayerNames);
      LOG(info) << "Loaded geometry configuration from file: " << filename << " with " << mLayerNames.size() << " layers.";
    }
    std::map<std::string, std::map<std::string, std::string>> getConfigurations() const { return mConfigurations; }
    std::map<std::string, std::string> getConfiguration(const std::string& layerName) const;
    std::vector<std::string> getLayerNames() const { return mLayerNames; }
    bool hasValue(const std::string& layerName, const std::string& key) const;
    std::string getValue(const std::string& layerName, const std::string& key, bool require = true) const;
    void setValue(const std::string& layerName, const std::string& key, const std::string& value) { mConfigurations[layerName][key] = value; }
    void replaceValue(const std::string& layerName, const std::string& key, const std::string& value);
    float getFloatValue(const std::string& layerName, const std::string& key) const { return std::stof(getValue(layerName, key)); }
    int getIntValue(const std::string& layerName, const std::string& key) const { return std::stoi(getValue(layerName, key)); }

   private:
    std::string mFileName;                                                     // Filename of the geometry
    std::map<std::string, std::map<std::string, std::string>> mConfigurations; // Layer configurations
    std::vector<std::string> mLayerNames;                                      // Ordered names of the layers
  };

  // Add a geometry entry from a configuration file
  void addEntry(const std::string& filename) { entries.emplace_back(filename); }

  // Getters
  int getNumberOfConfigurations() const { return entries.size(); }
  const std::vector<GeometryEntry>& getEntries() const { return entries; }
  const GeometryEntry& getEntry(const int id) const { return entries.at(id); }
  GeometryEntry getGeometryEntry(const int id) const { return entries.at(id); }

  // Get configuration maps
  std::map<std::string, std::map<std::string, std::string>> getConfigurations(const int id) const { return entries.at(id).getConfigurations(); }
  std::map<std::string, std::string> getConfiguration(const int id, const std::string& layerName) const { return entries.at(id).getConfiguration(layerName); }

  // Get specific values
  std::string getValue(const int id, const std::string& layerName, const std::string& key, bool require = true) const { return entries.at(id).getValue(layerName, key, require); };
  float getFloatValue(const int id, const std::string& layerName, const std::string& key) const { return entries.at(id).getFloatValue(layerName, key); }

 private:
  std::vector<GeometryEntry> entries;
};

// +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+

// this class implements a synthetic smearer that allows
// for on-demand smearing of TrackParCovs in a certain flexible t
// detector layout.
class FastTracker
{
 public:
  // Constructor/destructor
  FastTracker() = default;
  // Destructor
  virtual ~FastTracker() {}

  // Layer and layer configuration
  DetLayer* AddLayer(TString name, float r, float z, float x0, float xrho, float resRPhi = 0.0f, float resZ = 0.0f, float eff = 0.0f, int type = 0);

  /// Add a dead region in phi for a specific layer
  /// \param layerName Name of the layer to modify
  /// \param phiStart Start angle of the dead region (in radians)
  /// \param phiEnd End angle of the dead region (in radians)
  void addDeadPhiRegionInLayer(const std::string& layerName, float phiStart, float phiEnd);
  DetLayer GetLayer(const int layer) const { return layers[layer]; }
  std::vector<DetLayer> GetLayers() const { return layers; }
  int GetLayerIndex(const std::string& name) const;
  size_t GetNLayers() const { return layers.size(); }
  bool IsLayerInert(const int layer) const { return layers[layer].isInert(); }
  void ClearLayers() { layers.clear(); }
  void SetRadiationLength(const std::string layerName, float x0) { layers[GetLayerIndex(layerName)].setRadiationLength(x0); }
  void SetRadius(const std::string layerName, float r) { layers[GetLayerIndex(layerName)].setRadius(r); }
  void SetResolutionRPhi(const std::string layerName, float resRPhi) { layers[GetLayerIndex(layerName)].setResolutionRPhi(resRPhi); }
  void SetResolutionZ(const std::string layerName, float resZ) { layers[GetLayerIndex(layerName)].setResolutionZ(resZ); }
  void SetResolution(const std::string layerName, float resRPhi, float resZ)
  {
    SetResolutionRPhi(layerName, resRPhi);
    SetResolutionZ(layerName, resZ);
  }

  void AddSiliconALICE3v4(std::vector<float> pixelResolution);
  void AddSiliconALICE3v2(std::vector<float> pixelResolution);
  void AddSiliconALICE3(float scaleX0VD, std::vector<float> pixelResolution);
  void AddTPC(float phiResMean, float zResMean);

  /**
   * @brief Adds a generic detector configuration from the specified file.
   *
   * This function loads and integrates a detector configuration into the tracker
   * using the provided filename. The file should contain the necessary parameters
   * and settings for the detector to be added.
   *
   * @param configMap Configuration map describing the detector.
   */
  void AddGenericDetector(GeometryContainer::GeometryEntry configMap, o2::ccdb::BasicCCDBManager* ccdbManager = nullptr);

  void Print();

  /**
   * @brief Performs fast tracking on the input track parameters.
   *
   * Propagates the given input track through the detector layers, applying
   * relevant corrections and updates, and stores the result in outputTrack.
   *
   * @param inputTrack The input track parameters and covariance (const, by value).
   * @param outputTrack Reference to the output track parameters and covariance, to be filled.
   * @param nch Charged particle multiplicity (used for hit density calculations).
   * @return int i.e. number of intercepts (implementation-defined).
   */
  int FastTrack(o2::track::TrackParCov inputTrack, o2::track::TrackParCov& outputTrack, const float nch);

  // For efficiency calculation
  float Dist(float z, float radius);
  float OneEventHitDensity(float multiplicity, float radius);
  float IntegratedHitDensity(float multiplicity, float radius);
  float UpcHitDensity(float radius);
  float HitDensity(float radius);
  float ProbGoodChiSqHit(float radius, float searchRadiusRPhi, float searchRadiusZ);

  // Setters and getters for configuration
  void SetIntegrationTime(float t) { integrationTime = t; }
  void SetMaxRadiusOfSlowDetectors(float r) { maxRadiusSlowDet = r; }
  void SetAvgRapidity(float y) { avgRapidity = y; }
  void SetdNdEtaCent(int d) { dNdEtaCent = d; }
  void SetLhcUPCscale(float s) { lhcUPCScale = s; }
  void SetBField(float b) { magneticField = b; }
  void SetMinRadTrack(float r) { fMinRadTrack = r; }
  void SetMagneticField(float b) { magneticField = b; }
  void SetApplyZacceptance(bool b) { mApplyZacceptance = b; }
  void SetApplyMSCorrection(bool b) { mApplyMSCorrection = b; }
  void SetApplyElossCorrection(bool b) { mApplyElossCorrection = b; }
  void SetApplyEffCorrection(bool b) { mApplyEffCorrection = b; }

  // Getters for the last track
  int GetNIntercepts() const { return nIntercepts; }
  int GetNSiliconPoints() const { return nSiliconPoints; }
  int GetNGasPoints() const { return nGasPoints; }
  float GetGoodHitProb(int layer) const
  {
    return (layer >= 0 && static_cast<size_t>(layer) < goodHitProbability.size()) ? goodHitProbability[layer] : 0.0f;
  }
  std::size_t GetNHits() const { return hits.size(); }
  float GetHitX(const int i) const { return hits[i][0]; }
  float GetHitY(const int i) const { return hits[i][1]; }
  float GetHitZ(const int i) const { return hits[i][2]; }
  uint64_t GetCovMatOK() const { return covMatOK; }
  uint64_t GetCovMatNotOK() const { return covMatNotOK; }

 private:
  // Definition of detector layers
  std::vector<DetLayer> layers;
  std::vector<std::vector<float>> hits; // bookkeep last added hits

  /// configuration parameters
  bool mApplyZacceptance = false;       /// check z acceptance or not
  bool mApplyMSCorrection = true;       /// Apply correction for multiple scattering
  bool mApplyElossCorrection = true;    /// Apply correction for eloss (requires MS correction)
  bool mApplyEffCorrection = true;      /// Apply correction for hit efficiency
  int mVerboseLevel = 0;                /// 0: not verbose, >0 more verbose
  const float mCrossSectionMinB = 8;    /// Minimum bias Cross section for event under study (PbPb MinBias ~ 8 Barns)
  int dNdEtaCent = 2200;                /// dN/deta e.g. at centrality 0-5% (for 5 TeV PbPb)
  int dNdEtaMinB = 1;                   /// dN/deta for minimum bias events
  float integrationTime = 0.02f;        /// Integration time in ms
  float magneticField = 20.f;           /// Magnetic field in kiloGauss (5 = 0.5T, 20 = 2T, etc)
  float covMatFactor = 0.99f;           /// covmat off-diagonal factor to use for covmat fix (negative: no factor)
  float sigmaD = 6.0f;                  /// sigma for the detector resolution in cm
  float luminosity = 1.e27f;            /// luminosity in cm^-2 s^-1 (e.g. 1.e27 for PbPb at 5 TeV)
  float otherBackground = 0.0f;         /// background from other sources, e.g. pileup, in [0, 1]
  float maxRadiusSlowDet = 10.f;        /// maximum radius of slow detectors in cm
  float avgRapidity = 0.45f;            /// average rapidity for hit density calculation
  float lhcUPCScale = 1.0f;             /// scale factor for LHC UPC events
  float upcBackgroundMultiplier = 1.0f; /// multiplier for UPC background
  float fMinRadTrack = 132.f;           /// minimum radius for track propagation in cm

  /// counters for covariance matrix statuses
  uint64_t covMatOK = 0;    /// cov mat has positive eigenvals
  uint64_t covMatNotOK = 0; /// cov mat has negative eigenvals

  /// last track information
  int nIntercepts = 0;    /// found in first outward propagation
  int nSiliconPoints = 0; /// silicon-based space points added to track
  int nGasPoints = 0;     /// tpc-based space points added to track
  std::vector<float> goodHitProbability;

  ClassDef(FastTracker, 1);
};

// +-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+-~-<*>-~-+

} // namespace fastsim
} // namespace o2

#endif // ALICE3_CORE_FASTTRACKER_H_
