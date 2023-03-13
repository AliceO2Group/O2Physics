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

/// \file FemtoDreamParticleHisto.h
/// \brief FemtoDreamParticleHisto - Histogram class for tracks, V0s and cascades
/// \author Andi Mathis, TU München, andreas.mathis@ph.tum.de

#ifndef PWGCF_FEMTODREAM_FEMTODREAMPARTICLEHISTO_H_
#define PWGCF_FEMTODREAM_FEMTODREAMPARTICLEHISTO_H_

#include <string>
#include "PWGCF/DataModel/FemtoDerived.h"
#include "Framework/HistogramRegistry.h"

using namespace o2::framework;

namespace o2::analysis::femtoDream
{

/// \class FemtoDreamParticleHisto
/// \brief Class for histogramming particle properties
/// \tparam particleType Type of the particle (Track/V0/Cascade/...)
/// \tparam suffixType (optional) Takes care of the suffix for the folder name in case of analyses of pairs of the same kind (T-T, V-V, C-C)
template <o2::aod::femtodreamparticle::ParticleType particleType, int suffixType = 0>
class FemtoDreamParticleHisto
{
 public:
  /// Destructor
  virtual ~FemtoDreamParticleHisto() = default;

  template <o2::aod::femtodreamMCparticle::MCType mc, typename T>
  void init_base(std::string folderName, std::string tempFitVarAxisTitle, T& tempFitVarAxis, T& tempFitVarpTAxis, bool isMC)
  {
      std::string folderSuffix = static_cast<std::string>(o2::aod::femtodreamMCparticle::MCTypeName[mc]).c_str(); 
      /// Histograms of the kinematic properties
      mHistogramRegistry->add((folderName + folderSuffix + "/hPt").c_str(), "; #it{p}_{T} (GeV/#it{c}); Entries", kTH1F, {{240, 0, 6}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hEta").c_str(), "; #eta; Entries", kTH1F, {{200, -1.5, 1.5}});
      mHistogramRegistry->add((folderName + folderSuffix + "/hPhi").c_str(), "; #phi; Entries", kTH1F, {{200, 0, 2. * M_PI}});
      
      /// particle specific histogramms for the TempFitVar column in FemtoDreamParticles
      if constexpr (o2::aod::femtodreamMCparticle::MCType::kRecon == mc){
        mHistogramRegistry->add((folderName + folderSuffix + static_cast<std::string>(o2::aod::femtodreamparticle::TempFitVarName[mParticleType])).c_str(), ("; #it{p}_{T} (GeV/#it{c}); " + tempFitVarAxisTitle).c_str(), kTH2F, {{tempFitVarpTAxis}, {tempFitVarAxis}});
      }

  }
  
  template <typename T>
  void init_MC(std::string folderName, std::string tempFitVarAxisTitle, T& tempFitVarAxis, T& tempFitVarpTAxis)
  {

    /// Particle-type specific histograms
    if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kTrack) {
      /// Track histograms
      mHistogramRegistry->add((folderName + static_cast<std::string>(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]).c_str() + "/hDCAxy_Primary").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
      mHistogramRegistry->add((folderName + static_cast<std::string>(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]).c_str() + "/hDCAxy_Daughter").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
      mHistogramRegistry->add((folderName + static_cast<std::string>(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]).c_str() + "/hDCAxy_Material").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
      mHistogramRegistry->add((folderName + static_cast<std::string>(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]).c_str() + "/hDCAxy_NotPrimary").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
      mHistogramRegistry->add((folderName + static_cast<std::string>(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]).c_str() + "/hDCAxy_Fake").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
      mHistogramRegistry->add((folderName + static_cast<std::string>(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]).c_str() + "/hDCAxy_DaughterLambda").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
      mHistogramRegistry->add((folderName + static_cast<std::string>(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth]).c_str() + "/hDCAxy_DaughterSigmaplus").c_str(), "; #it{p}_{T} (GeV/#it{c}); DCA_{xy} (cm)", kTH2F, {tempFitVarpTAxis, tempFitVarAxis});
    } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0) {
      /// V0 histograms
    } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
      /// Cascade histograms
    } else {
      LOG(fatal) << "FemtoDreamParticleHisto: Histogramming for requested object not defined - quitting!";
    }

    /*
    for (int theorigin = o2::aod::femtodreamMCparticle::ParticleOriginMCTruth::kPrimary; theorigin < o2::aod::femtodreamMCparticle::ParticleOriginMCTruth::kNOriginMCTruthTypes; theorigin++){
      mHistogramRegistry->add(folderName + (static_cast<std::string>(o2::aod::femtodreamparticle::TempFitVarName[mParticleType])
                                          + static_cast<std::string>(o2::aod::femtodreamMCparticle::ParticleOriginMCTruthName[theorigin])).c_str(),
                              ("; #it{p}_{T} (GeV/#it{c}); " + tempFitVarAxisTitle).c_str(), kTH2F, {{tempFitVarpTAxis}, {tempFitVarAxis}});
    }
    */


  }

  /// Initialization of the QA histograms
  /// \param registry HistogramRegistry
  template <typename T>
  void init(HistogramRegistry* registry, T& tempFitVarpTBins, T& tempFitVarBins, bool isMC)
  {
    if (registry) {
      mHistogramRegistry = registry;
      /// The folder names are defined by the type of the object and the suffix (if applicable)
      std::string tempFitVarAxisTitle;
      if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kTrack) {
        /// Track histograms
        tempFitVarAxisTitle = "DCA_{xy} (cm)";
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0) {
        /// V0 histograms
        tempFitVarAxisTitle = "cos#alpha";
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
        /// Cascade histograms
        tempFitVarAxisTitle = "cos#alpha";
      } else {
        LOG(fatal) << "FemtoDreamParticleHisto: Histogramming for requested object not defined - quitting!";
      }

      framework::AxisSpec tempFitVarpTAxis = {tempFitVarpTBins, "#it{p}_{T} (GeV/#it{c})"}; // the pT binning may vary
      framework::AxisSpec tempFitVarAxis = {tempFitVarBins, tempFitVarAxisTitle};
      
      std::string folderName = (static_cast<std::string>(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]).c_str() + static_cast<std::string>(mFolderSuffix[mFolderSuffixType])).c_str();

      //Fill here the actual histogramms by calling init_base and init_MC
      init_base<o2::aod::femtodreamMCparticle::MCType::kRecon>(folderName, tempFitVarAxisTitle, tempFitVarAxis, tempFitVarpTAxis, isMC);
      if(isMC){
        init_base<o2::aod::femtodreamMCparticle::MCType::kTruth>(folderName, tempFitVarAxisTitle, tempFitVarAxis, tempFitVarpTAxis, isMC);
        init_MC(folderName, tempFitVarAxisTitle, tempFitVarAxis, tempFitVarpTAxis);
      }

    }
  }

  /// Filling of the histograms
  /// \tparam T Data type of the particle
  /// \param part Particle
  template <o2::aod::femtodreamMCparticle::MCType mc, typename T>
  void fillQA_base(T const& part)
  {
    /// Histograms of the kinematic properties
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hPt"), part.pt());
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hEta"), part.eta());
    mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) + HIST("/hPhi"), part.phi());
        
    /// particle specific histogramms for the TempFitVar column in FemtoDreamParticles
    if constexpr (mc == o2::aod::femtodreamMCparticle::MCType::kRecon){
      mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST(o2::aod::femtodreamMCparticle::MCTypeName[mc]) 
                             + HIST(o2::aod::femtodreamparticle::TempFitVarName[mParticleType]), part.pt(), part.tempFitVar());
    }
  } 


  /// Filling of the histograms for MC
  /// \tparam T Data type of the particle
  /// \param part Particle
  /// \param mctruthorigin Origin of the associated mc Truth particle
  template <typename T>
  void fillQA_MC(T const& part, int mctruthorigin)
  {
    if (mHistogramRegistry) {
      

      if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kTrack) {
        /// Track histograms
        switch (mctruthorigin) {
          case (o2::aod::femtodreamMCparticle::kPrimary):
            mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Primary"),
                                     part.pt(), part.tempFitVar());
            break;
          case (o2::aod::femtodreamMCparticle::kDaughter):
            mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Daughter"),
                                     part.pt(), part.tempFitVar());
            break;
          case (o2::aod::femtodreamMCparticle::kMaterial):
            mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Material"),
                                     part.pt(), part.tempFitVar());
            break;
          case (o2::aod::femtodreamMCparticle::kNotPrimary):
            mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_NotPrimary"),
                                     part.pt(), part.tempFitVar());
            break;
          case (o2::aod::femtodreamMCparticle::kFake):
            mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_Fake"),
                                     part.pt(), part.tempFitVar());
            break;
          case (o2::aod::femtodreamMCparticle::kDaughterLambda):
            mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_DaughterLambda"),
                                     part.pt(), part.tempFitVar());
            break;
          case (o2::aod::femtodreamMCparticle::kDaughterSigmaplus):
            mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType]) + HIST(mFolderSuffix[mFolderSuffixType]) + HIST("_MC/hDCAxy_DaughterSigmaplus"),
                                     part.pt(), part.tempFitVar());
            break;
          default:
            LOG(fatal) << "femtodreamparticleMC: not known value for ParticleOriginMCTruth - please check. Quitting!";
        }
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kV0) {
        /// V0 histograms
      } else if constexpr (mParticleType == o2::aod::femtodreamparticle::ParticleType::kCascade) {
        /// Cascade histograms
      } else {
        LOG(fatal) << "FemtoDreamParticleHisto: Histogramming for requested object not defined - quitting!";
      }


      /*
      for (int theorigin = o2::aod::femtodreamMCparticle::ParticleOriginMCTruth::kPrimary; theorigin < o2::aod::femtodreamMCparticle::ParticleOriginMCTruth::kNOriginMCTruthTypes; theorigin++){
        if (theorigin == mctruthorigin){
          mHistogramRegistry->fill(HIST(o2::aod::femtodreamparticle::ParticleTypeName[mParticleType])
                                 + HIST(mFolderSuffix[mFolderSuffixType]) 
                                 + HIST(o2::aod::femtodreamMCparticle::MCTypeName[o2::aod::femtodreamMCparticle::MCType::kTruth])
                                 + HIST(o2::aod::femtodreamparticle::TempFitVarName[mParticleType])
                                 + HIST(o2::aod::femtodreamMCparticle::ParticleOriginMCTruthName[mctruthorigin]),
                                 part.pt(), part.tempFitVar());
        }
      }
      */

    }
  }

  
  template <bool isMC, typename T>
  void fillQA(T const& part)
  {
    std::string tempFitVarName;
    if (mHistogramRegistry) {
      
      fillQA_base<o2::aod::femtodreamMCparticle::MCType::kRecon>(part);
      if constexpr (isMC){
        if( part.has_femtoDreamMCParticle()){
          fillQA_base<o2::aod::femtodreamMCparticle::MCType::kTruth>(part.femtoDreamMCParticle());
          fillQA_MC(part, (part.femtoDreamMCParticle()).partOriginMCTruth());
        } 
      } 

    }
  }

 private:
  HistogramRegistry* mHistogramRegistry;                                                   ///< For QA output
  static constexpr o2::aod::femtodreamparticle::ParticleType mParticleType = particleType; ///< Type of the particle under analysis
  static constexpr int mFolderSuffixType = suffixType;                                     ///< Counter for the folder suffix specified below
  static constexpr std::string_view mFolderSuffix[3] = {"", "_one", "_two"};               ///< Suffix for the folder name in case of analyses of pairs of the same kind (T-T, V-V, C-C)
};
} // namespace o2::analysis::femtoDream

#endif // PWGCF_FEMTODREAM_FEMTODREAMPARTICLEHISTO_H_
