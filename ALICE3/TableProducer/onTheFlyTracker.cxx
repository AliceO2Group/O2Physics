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

/// \file onTheFlyTracker.cxx
///
/// \brief LUT-based on-the-fly analysis task-level tracking
///
/// This task allows for the calculation of aod::collisions and aod::Tracks in a synthetic manner,
/// smearing MC particles with very configurable settings. This will allow for the usage of
/// custom LUTs (obtained through separate studies) and the subsequent estimate of the performance
/// of a future detector even in very statistics-hungry analyses.
///
/// \author David Dobrigkeit Chinellato, UNICAMP
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>, UniBo
/// \author Roberto Preghenella preghenella@bo.infn.it
///

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

#include "TRandom.h"

using namespace o2;
using namespace o2::framework;

// Imported from DelphesO2
/// @author: Roberto Preghenella
/// @email: preghenella@bo.infn.it

#define LUTCOVM_VERSION 20210801

struct map_t {
  int nbins = 1;
  float min = 0.;
  float max = 1.e6;
  bool log = false;
  float eval(int bin)
  {
    float width = (max - min) / nbins;
    float val = min + (bin + 0.5) * width;
    if (log)
      return pow(10., val);
    return val;
  };
  int find(float val)
  {
    float width = (max - min) / nbins;
    int bin;
    if (log)
      bin = (int)((log10(val) - min) / width);
    else
      bin = (int)((val - min) / width);
    if (bin < 0)
      return 0;
    if (bin > nbins - 1)
      return nbins - 1;
    return bin;
  };
  void print() { printf("nbins = %d, min = %f, max = %f, log = %s \n", nbins, min, max, log ? "on" : "off"); };
};

struct lutHeader_t {
  int version = LUTCOVM_VERSION;
  int pdg = 0;
  float mass = 0.;
  float field = 0.;
  map_t nchmap;
  map_t radmap;
  map_t etamap;
  map_t ptmap;
  bool check_version()
  {
    return (version == LUTCOVM_VERSION);
  };
  void print()
  {
    printf(" version: %d \n", version);
    printf("     pdg: %d \n", pdg);
    printf("   field: %f \n", field);
    printf("  nchmap: ");
    nchmap.print();
    printf("  radmap: ");
    radmap.print();
    printf("  etamap: ");
    etamap.print();
    printf("   ptmap: ");
    ptmap.print();
  };
};

struct lutEntry_t {
  float nch = 0.;
  float eta = 0.;
  float pt = 0.;
  bool valid = false;
  float eff = 0.;
  float eff2 = 0.;
  float itof = 0.;
  float otof = 0.;
  float covm[15] = {0.};
  float eigval[5] = {0.};
  float eigvec[5][5] = {0.};
  float eiginv[5][5] = {0.};
  void print()
  {
    printf(" --- lutEntry: pt = %f, eta = %f (%s)\n", pt, eta, valid ? "valid" : "not valid");
    printf("     efficiency: %f\n", eff);
    printf("     covMatix: ");
    int k = 0;
    for (int i = 0; i < 5; ++i) {
      for (int j = 0; j < i + 1; ++j)
        printf("% e ", covm[k++]);
      printf("\n               ");
    }
    printf("\n");
  }
};

class TrackSmearer
{

 public:
  TrackSmearer() = default;
  ~TrackSmearer() = default;

  /** LUT methods **/
  bool loadTable(int pdg, const char* filename, bool forceReload = false)
  {
    auto ipdg = getIndexPDG(pdg);
    if (mLUTHeader[ipdg] && !forceReload) {
      std::cout << " --- LUT table for PDG " << pdg << " has been already loaded with index " << ipdg << std::endl;
      return false;
    }
    mLUTHeader[ipdg] = new lutHeader_t;

    std::ifstream lutFile(filename, std::ifstream::binary);
    if (!lutFile.is_open()) {
      std::cout << " --- cannot open covariance matrix file for PDG " << pdg << ": " << filename << std::endl;
      delete mLUTHeader[ipdg];
      mLUTHeader[ipdg] = nullptr;
      return false;
    }
    lutFile.read(reinterpret_cast<char*>(mLUTHeader[ipdg]), sizeof(lutHeader_t));
    if (lutFile.gcount() != sizeof(lutHeader_t)) {
      std::cout << " --- troubles reading covariance matrix header for PDG " << pdg << ": " << filename << std::endl;
      delete mLUTHeader[ipdg];
      mLUTHeader[ipdg] = nullptr;
      return false;
    }
    if (mLUTHeader[ipdg]->version != LUTCOVM_VERSION) {
      std::cout << " --- LUT header version mismatch: expected/detected = " << LUTCOVM_VERSION << "/" << mLUTHeader[ipdg]->version << std::endl;
      delete mLUTHeader[ipdg];
      mLUTHeader[ipdg] = nullptr;
      return false;
    }
    if (mLUTHeader[ipdg]->pdg != pdg) {
      std::cout << " --- LUT header PDG mismatch: expected/detected = " << pdg << "/" << mLUTHeader[ipdg]->pdg << std::endl;
      delete mLUTHeader[ipdg];
      mLUTHeader[ipdg] = nullptr;
      return false;
    }
    const int nnch = mLUTHeader[ipdg]->nchmap.nbins;
    const int nrad = mLUTHeader[ipdg]->radmap.nbins;
    const int neta = mLUTHeader[ipdg]->etamap.nbins;
    const int npt = mLUTHeader[ipdg]->ptmap.nbins;
    mLUTEntry[ipdg] = new lutEntry_t****[nnch];
    for (int inch = 0; inch < nnch; ++inch) {
      mLUTEntry[ipdg][inch] = new lutEntry_t***[nrad];
      for (int irad = 0; irad < nrad; ++irad) {
        mLUTEntry[ipdg][inch][irad] = new lutEntry_t**[neta];
        for (int ieta = 0; ieta < neta; ++ieta) {
          mLUTEntry[ipdg][inch][irad][ieta] = new lutEntry_t*[npt];
          for (int ipt = 0; ipt < npt; ++ipt) {
            mLUTEntry[ipdg][inch][irad][ieta][ipt] = new lutEntry_t;
            lutFile.read(reinterpret_cast<char*>(mLUTEntry[ipdg][inch][irad][ieta][ipt]), sizeof(lutEntry_t));
            if (lutFile.gcount() != sizeof(lutEntry_t)) {
              std::cout << " --- troubles reading covariance matrix entry for PDG " << pdg << ": " << filename << std::endl;
              return false;
            }
          }
        }
      }
    }
    std::cout << " --- read covariance matrix table for PDG " << pdg << ": " << filename << std::endl;
    mLUTHeader[ipdg]->print();

    lutFile.close();
    return true;
  }
  void useEfficiency(bool val) { mUseEfficiency = val; };
  void setWhatEfficiency(int val) { mWhatEfficiency = val; };
  lutHeader_t* getLUTHeader(int pdg) { return mLUTHeader[getIndexPDG(pdg)]; };
  lutEntry_t* getLUTEntry(int pdg, float nch, float radius, float eta, float pt)
  {
    auto ipdg = getIndexPDG(pdg);
    if (!mLUTHeader[ipdg])
      return nullptr;
    auto inch = mLUTHeader[ipdg]->nchmap.find(nch);
    auto irad = mLUTHeader[ipdg]->radmap.find(radius);
    auto ieta = mLUTHeader[ipdg]->etamap.find(eta);
    auto ipt = mLUTHeader[ipdg]->ptmap.find(pt);
    return mLUTEntry[ipdg][inch][irad][ieta][ipt];
  }

  using O2Track = o2::track::TrackParCov;
  bool smearTrack(O2Track& o2track, lutEntry_t* lutEntry)
  {
    // generate efficiency
    if (mUseEfficiency) {
      auto eff = 0.;
      if (mWhatEfficiency == 1)
        eff = lutEntry->eff;
      if (mWhatEfficiency == 2)
        eff = lutEntry->eff2;
      if (gRandom->Uniform() > eff)
        return false;
    }
    // transform params vector and smear
    double params_[5];
    for (int i = 0; i < 5; ++i) {
      double val = 0.;
      for (int j = 0; j < 5; ++j)
        val += lutEntry->eigvec[j][i] * o2track.getParam(j);
      params_[i] = gRandom->Gaus(val, sqrt(lutEntry->eigval[i]));
    }
    // transform back params vector
    for (int i = 0; i < 5; ++i) {
      double val = 0.;
      for (int j = 0; j < 5; ++j)
        val += lutEntry->eiginv[j][i] * params_[j];
      o2track.setParam(val, i);
    }
    // should make a sanity check that par[2] sin(phi) is in [-1, 1]
    if (fabs(o2track.getParam(2)) > 1.) {
      std::cout << " --- smearTrack failed sin(phi) sanity check: " << o2track.getParam(2) << std::endl;
    }
    // set covariance matrix
    for (int i = 0; i < 15; ++i)
      o2track.setCov(lutEntry->covm[i], i);
    return true;
  }
  bool smearTrack(O2Track& o2track, int pid, float nch)
  {

    auto pt = o2track.getPt();
    if (abs(pid) == 1000020030) {
      pt *= 2.f;
    }
    auto eta = o2track.getEta();
    auto lutEntry = getLUTEntry(pid, nch, 0., eta, pt);
    if (!lutEntry || !lutEntry->valid)
      return false;
    return smearTrack(o2track, lutEntry);
  }

  int getIndexPDG(int pdg)
  {
    switch (abs(pdg)) {
      case 11:
        return 0; // Electron
      case 13:
        return 1; // Muon
      case 211:
        return 2; // Pion
      case 321:
        return 3; // Kaon
      case 2212:
        return 4; // Proton
      case 1000010020:
        return 5; // Deuteron
      case 1000010030:
        return 6; // Triton
      case 1000020030:
        return 7; // Helium3
      default:
        return 2; // Default: pion
    };
  };

  void setdNdEta(float val) { mdNdEta = val; };

 protected:
  static constexpr unsigned int nLUTs = 8; // Number of LUT available
  lutHeader_t* mLUTHeader[nLUTs] = {nullptr};
  lutEntry_t***** mLUTEntry[nLUTs] = {nullptr};
  bool mUseEfficiency = true;
  int mWhatEfficiency = 1;
  float mdNdEta = 1600.;
};

struct OnTheFlyTracker {
  Produces<aod::Collisions> collisions;
  Produces<aod::McCollisionLabels> collLabels;
  Produces<aod::StoredTracks> tracksPar;
  Produces<aod::TracksExtension> tracksParExtension;
  Produces<aod::StoredTracksCov> tracksParCov;
  Produces<aod::TracksCovExtension> tracksParCovExtension;
  Produces<aod::McTrackLabels> tracksLabels;
  Produces<aod::TracksDCA> tracksDCA;

  Configurable<float> maxEta{"maxEta", 1.5, "maximum eta to consider viable"};
  Configurable<float> multEtaRange{"multEtaRange", 0.8, "eta range to compute the multiplicity"};
  Configurable<float> minPt{"minPt", 0.1, "minimum pt to consider viable"};
  Configurable<bool> enableLUT{"enableLUT", false, "Enable track smearing"};
  Configurable<bool> enableNucleiSmearing{"enableNucleiSmearing", false, "Enable smearing of nuclei"};

  bool fillTracksDCA = false;

  // necessary for particle charges
  Service<O2DatabasePDG> pdgDB;

  // for handling basic QA histograms if requested
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrNONE;

  // Track smearer
  TrackSmearer mSmearer;

  void init(o2::framework::InitContext& initContext)
  {
    // Checking if the tables are requested in the workflow and enabling them
    auto& workflows = initContext.services().get<RunningWorkflowInfo const>();
    for (DeviceSpec const& device : workflows.devices) {
      for (auto const& input : device.inputs) {
        if (input.matcher.binding == "TracksDCA") {
          fillTracksDCA = true;
        }
      }
    }

    if (enableLUT) {
      std::map<int, const char*> mapPdgLut;
      mapPdgLut.insert(std::make_pair(11, "lutCovm.el.dat"));
      mapPdgLut.insert(std::make_pair(13, "lutCovm.mu.dat"));
      mapPdgLut.insert(std::make_pair(211, "lutCovm.pi.dat"));
      mapPdgLut.insert(std::make_pair(321, "lutCovm.ka.dat"));
      mapPdgLut.insert(std::make_pair(2212, "lutCovm.pr.dat"));
      if (enableNucleiSmearing) {
        mapPdgLut.insert(std::make_pair(1000010020, "lutCovm.de.dat"));
        mapPdgLut.insert(std::make_pair(1000010030, "lutCovm.tr.dat"));
        mapPdgLut.insert(std::make_pair(1000020030, "lutCovm.he3.dat"));
      }
      for (auto e : mapPdgLut) {
        if (!mSmearer.loadTable(e.first, e.second)) {
          LOG(fatal) << "Having issue with loading the LUT " << e.first << " " << e.second;
        }
      }
    }

    // Basic QA
    const AxisSpec axisMomentum{static_cast<int>(100), 0.0f, +10.0f, "#it{p} (GeV/#it{c})"};
    histos.add("hPt", "hPt", kTH1F, {axisMomentum});
  }

  /// Function to convert a McParticle into a perfect Track
  /// \param particle the particle to convert (mcParticle)
  /// \param o2track the address of the resulting TrackParCov
  template <typename McParticleType>
  void convertMCParticleToO2Track(McParticleType& particle, o2::track::TrackParCov& o2track)
  {
    auto pdgInfo = pdgDB->GetParticle(particle.pdgCode());
    int charge = 0;
    if (pdgInfo != nullptr) {
      charge = pdgInfo->Charge();
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

    // Initialize TrackParCov in-place
    new (&o2track)(o2::track::TrackParCov)(x, particle.phi(), params, covm);
  }

  /// Function to fill track parameter table
  /// \param coll collision (for index)
  /// \param trackType type of created track
  /// \param trackPar track for parameters
  template <typename CollType, typename TTrackPar>
  void fillTracksPar(CollType& coll, aod::track::TrackTypeEnum trackType, TTrackPar& trackPar)
  {
    tracksPar(coll.globalIndex(), trackType, trackPar.getX(), trackPar.getAlpha(), trackPar.getY(), trackPar.getZ(), trackPar.getSnp(), trackPar.getTgl(), trackPar.getQ2Pt());
    tracksParExtension(trackPar.getPt(), trackPar.getP(), trackPar.getEta(), trackPar.getPhi());
  }

  float dNdEta = 0.f; // Charged particle multiplicity to use in the efficiency evaluation
  void process(aod::McCollision const& mcCollision, aod::McParticles const& mcParticles)
  {
    o2::dataformats::DCA dcaInfoCov;
    o2::dataformats::VertexBase vtx;
    // First we compute the number of charged particles in the event
    dNdEta = 0.f;
    for (const auto& mcParticle : mcParticles) {
      if (TMath::Abs(mcParticle.eta()) > multEtaRange) {
        continue;
      }
      if (mcParticle.has_daughters()) {
        continue;
      }
      const auto& pdgInfo = pdgDB->GetParticle(mcParticle.pdgCode());
      if (!pdgInfo) {
        LOG(warning) << "PDG code " << mcParticle.pdgCode() << " not found in the database";
        continue;
      }
      if (pdgInfo->Charge() == 0) {
        continue;
      }
      dNdEta += 1.f;
    }

    for (const auto& mcParticle : mcParticles) {
      const auto pdg = std::abs(mcParticle.pdgCode());
      if (pdg != kElectron && pdg != kMuonMinus && pdg != kPiPlus && pdg != kKPlus && pdg != kProton) {
        continue;
      }
      if (std::fabs(mcParticle.eta()) > maxEta) {
        continue;
      }
      if (mcParticle.pt() < minPt) {
        continue;
      }
      o2::track::TrackParCov trackParCov;
      convertMCParticleToO2Track(mcParticle, trackParCov);

      if (!mSmearer.smearTrack(trackParCov, mcParticle.pdgCode(), dNdEta)) {
        continue;
      }

      // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*
      // Calculate primary vertex
      // To be added once smeared tracks are in place
      // *+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*+~+*

      // Base QA
      histos.fill(HIST("hPt"), trackParCov.getPt());

      // Fixme: collision index could be changeable
      aod::track::TrackTypeEnum trackType = aod::track::Track;
      fillTracksPar(mcCollision, trackType, trackParCov);
      if (fillTracksDCA) {
        tracksDCA(1e-3, 1e-3);
      }
      // TODO do we keep the rho as 0? Also the sigma's are duplicated information
      tracksParCov(std::sqrt(trackParCov.getSigmaY2()), std::sqrt(trackParCov.getSigmaZ2()), std::sqrt(trackParCov.getSigmaSnp2()),
                   std::sqrt(trackParCov.getSigmaTgl2()), std::sqrt(trackParCov.getSigma1Pt2()), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
      tracksParCovExtension(trackParCov.getSigmaY2(), trackParCov.getSigmaZY(), trackParCov.getSigmaZ2(), trackParCov.getSigmaSnpY(),
                            trackParCov.getSigmaSnpZ(), trackParCov.getSigmaSnp2(), trackParCov.getSigmaTglY(), trackParCov.getSigmaTglZ(), trackParCov.getSigmaTglSnp(),
                            trackParCov.getSigmaTgl2(), trackParCov.getSigma1PtY(), trackParCov.getSigma1PtZ(), trackParCov.getSigma1PtSnp(), trackParCov.getSigma1PtTgl(),
                            trackParCov.getSigma1Pt2());
      tracksLabels(mcParticle.globalIndex(), 0);
    }
    collisions(-1, // BC is irrelevant in synthetic MC tests for now, could be adjusted in future
               mcCollision.posX(), mcCollision.posY(), mcCollision.posZ(),
               1e-3, 0.0, 1e-3, 0.0, 0.0, 1e-3,
               0, 1e-3, mcParticles.size(),
               0, 0);
    collLabels(mcCollision.globalIndex(), 0);
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<OnTheFlyTracker>(cfgc)};
}
