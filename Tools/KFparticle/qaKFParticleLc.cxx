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
/// \file   qaKFParticleLc.cxx
/// \author Annalena Kalteyer <annalena.sophie.kalteyer@cern.ch>, GSI Darmstadt
/// \brief  Task to test the performance of the KFParticle package on the Lc to pKpi decay
///

#include "Tools/KFparticle/qaKFParticleLc.h"

#include "TableHelper.h"

#include <CCDB/BasicCCDBManager.h>

#include <TDatabasePDG.h>
#include <TPDGCode.h>

#include <string>

/// includes O2
#include "DataFormatsParameters/GRPMagField.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/Propagator.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/runDataProcessing.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"

/// includes O2Physics
#include "Common/Core/RecoDecay.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/Core/trackUtilities.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Tools/KFparticle/KFUtilities.h"

/// includes KFParticle
#include "KFPTrack.h"
#include "KFPVertex.h"
#include "KFParticle.h"
#include "KFParticleBase.h"
#include "KFVertex.h"

#ifndef HomogeneousField

#define HomogeneousField

#endif

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::dataformats;

struct qaKFParticleLc {

  /// general steering settings
  Configurable<bool> isRun3{"isRun3", true, "Is Run3 dataset"};
  Configurable<std::string> ccdbUrl{"ccdburl", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPathLut{"ccdbPathLut", "GLO/Param/MatLUT", "Path for LUT parametrization"};
  Configurable<std::string> ccdbPathGrp{"ccdbPathGrp", "GLO/GRP/GRP", "Path of the grp file (Run 2)"};
  Configurable<std::string> ccdbPathGrpMag{"ccdbPathGrpMag", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object (Run 3)"};
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  o2::base::MatLayerCylSet* lut;
  o2::base::Propagator::MatCorrType matCorr = o2::base::Propagator::MatCorrType::USEMatCorrLUT;
  int runNumber;
  double magneticField = 0.;
  KFParticle KFPion, KFKaon, KFProton, KFLc, KFLc_PV;

  /// option to select good events
  Configurable<bool> eventSelection{"eventSelection", true, "select good events"}; // currently only sel8 is defined for run3
  /// options to select only specific tracks
  Configurable<int> trackSelection{"trackSelection", 1, "Track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks"};

  /// Particle Identification
  // TPC PID
  Configurable<double> nSigmaTpcMaxPi{"nSigmaTpcMaxPi", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcMaxKa{"nSigmaTpcMaxKa", 3., "Nsigma cut on TPC only"};
  Configurable<double> nSigmaTpcMaxPr{"nSigmaTpcMaxPr", 3., "Nsigma cut on TPC only"};
  // TOF PID
  Configurable<double> ptPidTofMinPi{"ptPidTofMinPi", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMaxPi{"nSigmaTofMaxPi", 3., "Nsigma cut on TOF only"};
  Configurable<double> ptPidTofMinKa{"ptPidTofMinKa", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMaxKa{"nSigmaTofMaxKa", 3., "Nsigma cut on TOF only"};
  Configurable<double> ptPidTofMinPr{"ptPidTofMinPr", 0.15, "Lower bound of track pT for TOF PID"};
  Configurable<double> nSigmaTofMaxPr{"nSigmaTofMaxPr", 3., "Nsigma cut on TOF only"};
  // TPC & TOF Combined
  Configurable<double> nSigmaCombMaxPi{"nSigmaCombMaxPi", 3., "Nsigma cut on TPC & TOF"};
  Configurable<double> nSigmaCombMaxKa{"nSigmaCombMaxKa", 3., "Nsigma cut on TPC & TOF"};
  Configurable<double> nSigmaCombMaxPr{"nSigmaCombMaxPr", 3., "Nsigma cut on TPC & TOF"};

  /// singe track selections
  Configurable<float> d_pTMin{"d_pTMin", 0.3, "minimum momentum for tracks"};
  Configurable<float> d_etaRange{"d_etaRange", 0.8, "eta Range for tracks"};
  Configurable<float> d_dcaXYTrackPV{"d_dcaXYTrackPV", 2., "DCA XY of the daughter tracks to the PV"};
  Configurable<float> d_dcaZTrackPV{"d_dcaZTrackPV", 10., "DCA Z of the daughter tracks to the PV"};
  /// Lc Daughter selections
  Configurable<float> d_PtMinPi{"d_PtMinPi", 0., "minimum momentum for Pi from Lc candidates"};
  Configurable<float> d_PtMinKa{"d_PtMinKa", 0., "minimum momentum for Ka from Lc candidates"};
  Configurable<float> d_PtMinPr{"d_PtMinPr", 0., "minimum momentum for Pr from Lc candidates"};
  Configurable<float> d_dist3DSVDau{"d_dist3DSVDau", 1000., "maximum geometrical distance 3D daughter tracks at the SV"};
  /// Lc selection after geometrical fitting
  Configurable<float> d_pTMinLc{"d_pTMinLc", 0., "minimum momentum for Lc candidates"};
  Configurable<float> d_pTMaxLc{"d_pTMaxLc", 36., "maximum momentum for Lc candidates"};
  Configurable<float> d_massMin{"d_massMin", 1.65, "minimum mass"};
  Configurable<float> d_massMax{"d_massMax", 2.08, "minimum mass"};
  Configurable<float> d_cosPA{"d_cosPA", -1., "minimum cosine Pointing angle"};
  Configurable<float> d_cosPAXY{"d_cosPAXY", -1., "minimum cosine Pointing angle"};
  Configurable<float> d_distPVSV{"d_distPVSV", -1., "minimum distance between PV and SV"};
  Configurable<float> d_chi2geo{"d_chi2geo", 1000., "maximum chi2 geometrical"};
  /// Lc selection after topological constrain to the PV
  Configurable<bool> applySelectionWithTopoConst{"applySelectionWithTopoConst", true, "Apply selections constraining the mother to the PV"};
  Configurable<float> d_decayLength{"d_decayLength", 0., "minimum decay length"};
  Configurable<float> d_normdecayLength{"d_normdecayLength", 100., "minimum normalised decay length"};
  Configurable<float> d_chi2topo{"d_chi2topo", 1000., "maximum chi2 topological"};
  /// Option to write D0 variables in a tree
  Configurable<double> d_DwnSmplFact{"d_DwnSmplFact", 1., "Downsampling factor for tree"};

  // Define which track selection should be used:
  // 0 -> No track selection is applied
  // 1 kGlobalTrack = kQualityTracks | kPrimaryTracks | kInAcceptanceTracks
  //        kQualityTracks = kTrackType | kTPCNCls | kTPCCrossedRows | kTPCCrossedRowsOverNCls | kTPCChi2NDF |
  //                         kTPCRefit | kITSNCls | kITSChi2NDF | kITSRefit | kITSHits
  //        kPrimaryTracks = kGoldenChi2 | kDCAxy | kDCAz
  //        kInAcceptanceTracks = kPtRange | kEtaRange
  // 2 kGlobalTrackWoPtEta = kQualityTracks | kPrimaryTracks
  // 3 kGlobalTrackWoDCA = kQualityTracks | kInAcceptanceTracks
  // 4 kQualityTracks
  // 5 kInAcceptanceTracks
  Filter trackFilter = (trackSelection.node() == 0) ||
                       ((trackSelection.node() == 1) && requireGlobalTrackInFilter()) ||
                       ((trackSelection.node() == 2) && requireGlobalTrackWoPtEtaInFilter()) ||
                       ((trackSelection.node() == 3) && requireGlobalTrackWoDCAInFilter()) ||
                       ((trackSelection.node() == 4) && requireQualityTracksInFilter()) ||
                       ((trackSelection.node() == 5) && requireTrackCutInFilter(TrackSelectionFlags::kInAcceptanceTracks));

  Filter trackFilterEta = (nabs(aod::track::eta) < 0.8f);
  Filter trackFilterPTMin = (aod::track::pt > d_pTMin);

  Filter eventFilter = (o2::aod::evsel::sel8 == true);

  using CollisionTableData = soa::Join<aod::Collisions, aod::EvSels>;
  using BigTracks = soa::Join<aod::Tracks, aod::TracksCov, aod::TracksExtra>;
  using BigTracksExtended = soa::Join<BigTracks, aod::TracksDCA>;
  using BigTracksPID = soa::Join<BigTracksExtended, aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi, aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi, aod::pidTOFFullKa, aod::pidTOFFullPr>;
  using TrackTableData = soa::Join<BigTracksPID, aod::TrackSelection>;

  /// Table to be produced
  Produces<o2::aod::TreeKFLc> rowKFLc;

  void initMagneticFieldCCDB(o2::aod::BCsWithTimestamps::iterator const& bc, int& mRunNumber,
                             o2::framework::Service<o2::ccdb::BasicCCDBManager> const& ccdb, std::string ccdbPathGrp, o2::base::MatLayerCylSet* lut,
                             bool isRun3)
  {

    if (mRunNumber != bc.runNumber()) {

      LOGF(info, "====== initCCDB function called (isRun3==%d)", isRun3);
      if (!isRun3) { // Run 2 GRP object
        o2::parameters::GRPObject* grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>(ccdbPathGrp, bc.timestamp());
        if (grpo == nullptr) {
          LOGF(fatal, "Run 2 GRP object (type o2::parameters::GRPObject) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        }
        o2::base::Propagator::initFieldFromGRP(grpo);
        o2::base::Propagator::Instance()->setMatLUT(lut);
        LOGF(info, "Setting magnetic field to %d kG for run %d from its GRP CCDB object (type o2::parameters::GRPObject)", grpo->getNominalL3Field(), bc.runNumber());
      } else { // Run 3 GRP object
        o2::parameters::GRPMagField* grpo = ccdb->getForTimeStamp<o2::parameters::GRPMagField>(ccdbPathGrp, bc.timestamp());
        if (grpo == nullptr) {
          LOGF(fatal, "Run 3 GRP object (type o2::parameters::GRPMagField) is not available in CCDB for run=%d at timestamp=%llu", bc.runNumber(), bc.timestamp());
        }
        o2::base::Propagator::initFieldFromGRP(grpo);
        o2::base::Propagator::Instance()->setMatLUT(lut);
        LOGF(info, "Setting magnetic field to current %f A for run %d from its GRP CCDB object (type o2::parameters::GRPMagField)", grpo->getL3Current(), bc.runNumber());
      }
      mRunNumber = bc.runNumber();
    }
  } /// end initMagneticFieldCCDB

  void init(InitContext const&)
  {
    ccdb->setURL(ccdbUrl);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    lut = o2::base::MatLayerCylSet::rectifyPtrFromFile(ccdb->get<o2::base::MatLayerCylSet>(ccdbPathLut));
    runNumber = 0;
  } /// End init

  /// Function for single track selection
  template <typename T>
  bool isSelectedTracks(const T& track1, const T& track2, const T& track3)
  {
    /// DCA XY of the daughter tracks to the primaty vertex
    if ((fabs(track1.dcaXY()) > d_dcaXYTrackPV) || (fabs(track2.dcaXY()) > d_dcaXYTrackPV) || (fabs(track3.dcaXY()) > d_dcaXYTrackPV)) {
      return false;
    }
    /// DCA Z of the daughter tracks to the primaty vertex
    if ((fabs(track1.dcaZ()) > d_dcaZTrackPV) || (fabs(track2.dcaZ()) > d_dcaZTrackPV) || (fabs(track3.dcaZ()) > d_dcaZTrackPV)) {
      return false;
    }
    /// reject if the tracks with sum of charge != 1
    if (abs(track1.sign() + track2.sign() + track3.sign()) != 1) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isSelectedLc(const T& trackKaon, const T& trackPion, const T& trackProton)
  {
    if (!(trackKaon.sign() == -1 && trackPion.sign() == 1 && trackProton.sign() == 1)) {
      return false;
    }
    bool pidKaon = SelectPIDCombined(trackKaon, kKPlus);
    bool pidPion = SelectPIDCombined(trackPion, kPiPlus);
    bool pidProton = SelectPIDCombined(trackProton, kProton);
    if (!(pidKaon && pidPion && pidProton)) {
      return false;
    }
    return true;
  }
  template <typename T>
  bool isSelectedLcBar(const T& trackKaon, const T& trackPion, const T& trackProton)
  {
    if (!(trackKaon.sign() == 1 && trackPion.sign() == -1 && trackProton.sign() == -1)) {
      return false;
    }
    bool pidKaon = SelectPIDCombined(trackKaon, kKPlus);
    bool pidPion = SelectPIDCombined(trackPion, kPiPlus);
    bool pidProton = SelectPIDCombined(trackProton, kProton);
    if (!(pidKaon && pidPion && pidProton)) {
      return false;
    }
    return true;
  }

  template <typename T, typename T2>
  bool ReconstructLc(const T& trackKaon, const T& trackPion, const T& trackProton, const T2& KFPV)
  {
    KFPTrack kfpTrackKa = createKFPTrackFromTrack(trackKaon);
    KFPTrack kfpTrackPi = createKFPTrackFromTrack(trackPion);
    KFPTrack kfpTrackPr = createKFPTrackFromTrack(trackProton);
    KFParticle KFKa(kfpTrackKa, 321);
    KFKaon = KFKa;
    KFParticle KFPi(kfpTrackPi, 211);
    KFPion = KFPi;
    KFParticle KFPr(kfpTrackPr, 2212);
    KFProton = KFPr;
    const KFParticle* LcDaughters[3] = {&KFKaon, &KFPion, &KFProton};
    KFLc.SetConstructMethod(2);
    KFLc.Construct(LcDaughters, 3);
    if (!isSelectedDaughters(KFPion, KFKaon, KFProton)) {
      return false;
    }
    if (!isSelectedGeo(KFLc, KFPV)) {
      return false;
    }
    KFLc_PV = KFLc;
    KFLc_PV.SetProductionVertex(KFPV);
    if (applySelectionWithTopoConst) {
      if (!isSelectedLcTopo(KFLc_PV, KFPV)) {
        return false;
      }
    }
    return true;
  }

  template <typename T>
  bool isSelectedDaughters(const T& KFPion, const T& KFKaon, const T& KFProton)
  {
    /// Pt of daughters
    if ((KFPion.GetPt() < d_PtMinPi) || (KFKaon.GetPt() < d_PtMinKa) || (KFProton.GetPt() < d_PtMinPr)) {
      return false;
    }
    /// distance 3D daughter tracks at the secondary vertex
    if ((KFPion.GetDistanceFromParticle(KFKaon) > d_dist3DSVDau) || (KFPion.GetDistanceFromParticle(KFProton) > d_dist3DSVDau) || (KFKaon.GetDistanceFromParticle(KFProton) > d_dist3DSVDau)) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isSelectedGeo(const T& KFLc, const T& KFPV)
  {
    /// Pt selection
    if (KFLc.GetPt() < d_pTMinLc || KFLc.GetPt() > d_pTMaxLc) {
      return false;
    }
    /// Mass window selection
    if (KFLc.GetMass() < d_massMin || KFLc.GetMass() > d_massMax) {
      return false;
    }
    /// cosine pointing angle selection
    if (cpaFromKF(KFLc, KFPV) < d_cosPA) {
      return false;
    }
    /// cosine pointing XY angle selection
    if (cpaXYFromKF(KFLc, KFPV) < d_cosPAXY) {
      return false;
    }
    /// Minimum distance between PV and SV
    if (KFLc.GetDistanceFromVertex(KFPV) < d_distPVSV) {
      return false;
    }
    /// chi2 geometrical
    float chi2geo = KFLc.GetChi2() / KFLc.GetNDF();
    if (chi2geo > d_chi2geo) {
      return false;
    }
    return true;
  }

  template <typename T>
  bool isSelectedLcTopo(const T& KFLc_PV, const T& /*KFPV*/)
  {
    /// Pt selection
    if (KFLc_PV.GetPt() < d_pTMinLc || KFLc_PV.GetPt() > d_pTMaxLc) {
      return false;
    }
    /// Mass window selection
    if (KFLc_PV.GetMass() < d_massMin || KFLc_PV.GetMass() > d_massMax) {
      return false;
    }
    /// decay length selection
    if (KFLc_PV.GetDecayLength() < d_decayLength) {
      return false;
    }
    /// decay length error selection
    float normdecayLength = KFLc_PV.GetDecayLength() / KFLc_PV.GetErrDecayLength();
    if (normdecayLength < d_normdecayLength) {
      return false;
    }
    /// chi2 topological
    float chi2topo = KFLc_PV.GetChi2() / KFLc_PV.GetNDF();
    if (chi2topo > d_chi2topo) {
      return false;
    }
    return true;
  }

  template <typename T1>
  bool SelectPIDCombined(const T1& track, int particle)
  {
    switch (particle) {
      case kPiPlus: {
        if ((track.pt() <= ptPidTofMinPi) && track.hasTPC() && (abs(track.tpcNSigmaPi()) < nSigmaTpcMaxPi)) {
          return true;
        } else if ((track.pt() > ptPidTofMinPi) && track.hasTPC() && !track.hasTOF() && (abs(track.tpcNSigmaPi()) < nSigmaTpcMaxPi)) {
          return true;
        } else if ((track.pt() > ptPidTofMinPi) && !track.hasTPC() && track.hasTOF() && (abs(track.tofNSigmaPi()) < nSigmaTofMaxPi)) {
          return true;
        } else if ((track.pt() > ptPidTofMinPi) && track.hasTPC() && track.hasTOF()) {
          float CombinednSigma = 1. / sqrt(2) * sqrt((track.tpcNSigmaPi() * track.tpcNSigmaPi()) + (track.tofNSigmaPi() * track.tofNSigmaPi()));
          if (abs(CombinednSigma) < nSigmaCombMaxPi) {
            return true;
          } else {
            return false;
          }
        } else {
          return false;
        }
        break;
      }
      case kKPlus: {
        if ((track.pt() <= ptPidTofMinKa) && track.hasTPC() && (abs(track.tpcNSigmaKa()) < nSigmaTpcMaxKa)) {
          return true;
        } else if ((track.pt() > ptPidTofMinKa) && track.hasTPC() && !track.hasTOF() && (abs(track.tpcNSigmaKa()) < nSigmaTpcMaxKa)) {
          return true;
        } else if ((track.pt() > ptPidTofMinKa) && !track.hasTPC() && track.hasTOF() && (abs(track.tofNSigmaKa()) < nSigmaTofMaxKa)) {
          return true;
        } else if ((track.pt() > ptPidTofMinKa) && track.hasTPC() && track.hasTOF()) {
          float CombinednSigma = 1. / sqrt(2) * sqrt((track.tpcNSigmaKa() * track.tpcNSigmaKa()) + (track.tofNSigmaKa() * track.tofNSigmaKa()));
          if (abs(CombinednSigma) < nSigmaCombMaxKa) {
            return true;
          } else {
            return false;
          }
        } else {
          return false;
        }
        break;
      }
      case kProton: {
        if ((track.pt() <= ptPidTofMinPr) && track.hasTPC() && (abs(track.tpcNSigmaPr()) < nSigmaTpcMaxPr)) {
          return true;
        } else if ((track.pt() > ptPidTofMinPr) && track.hasTPC() && !track.hasTOF() && (abs(track.tpcNSigmaPr()) < nSigmaTpcMaxPr)) {
          return true;
        } else if ((track.pt() > ptPidTofMinPr) && !track.hasTPC() && track.hasTOF() && (abs(track.tofNSigmaPr()) < nSigmaTofMaxPr)) {
          return true;
        } else if ((track.pt() > ptPidTofMinPr) && track.hasTPC() && track.hasTOF()) {
          float CombinednSigma = 1. / sqrt(2) * sqrt((track.tpcNSigmaPr() * track.tpcNSigmaPr()) + (track.tofNSigmaPr() * track.tofNSigmaPr()));
          if (abs(CombinednSigma) < nSigmaCombMaxPr) {
            return true;
          } else {
            return false;
          }
        } else {
          return false;
        }
        break;
      }
      default: {
        LOGF(error, "ERROR: Species is not implemented");
        return false;
      }
    }
  }

  template <typename T1, typename T2>
  void writeVarTree(const T2& KFPion, const T2& KFKaon, const T2& KFProton, const T2& KFLc, const T2& KFLc_PV, const T2& KFPV, const T1& trackKa, const T1& trackPi, const T1& trackPr, const int source)
  {

    float chi2geo = KFLc.GetChi2() / KFLc.GetNDF();
    float normdecayLength = KFLc_PV.GetDecayLength() / KFLc_PV.GetErrDecayLength();
    float chi2topo = KFLc_PV.GetChi2() / KFLc_PV.GetNDF();
    const double pseudoRndm = trackKa.pt() * 1000. - (int64_t)(trackKa.pt() * 1000);
    if (pseudoRndm < d_DwnSmplFact) {
      /// Filling the D0 tree
      rowKFLc(runNumber,
              KFPion.GetPt(),
              KFKaon.GetPt(),
              KFProton.GetPt(),
              KFPion.GetDistanceFromVertexXY(KFPV),
              KFKaon.GetDistanceFromVertexXY(KFPV),
              KFProton.GetDistanceFromVertexXY(KFPV),
              KFPion.GetDistanceFromVertex(KFPV),
              KFKaon.GetDistanceFromVertex(KFPV),
              KFProton.GetDistanceFromVertex(KFPV),
              trackPi.tpcNSigmaPi(),
              trackKa.tpcNSigmaKa(),
              trackPr.tpcNSigmaPr(),
              trackPi.tofNSigmaPi(),
              trackKa.tofNSigmaKa(),
              trackPr.tofNSigmaPr(),
              KFPion.GetDistanceFromVertexXY(KFLc),
              KFKaon.GetDistanceFromVertexXY(KFLc),
              KFProton.GetDistanceFromVertexXY(KFLc),
              KFLc.GetPt(),
              KFLc.GetMass(),
              cpaFromKF(KFLc, KFPV),
              cpaXYFromKF(KFLc, KFPV),
              KFLc.GetDistanceFromVertex(KFPV),
              KFLc.GetDistanceFromVertexXY(KFPV),
              chi2geo,
              KFLc_PV.GetPt(),
              KFLc_PV.GetMass(),
              KFLc_PV.GetDecayLength(),
              KFLc_PV.GetDecayLengthXY(),
              cpaFromKF(KFLc_PV, KFPV),
              KFLc_PV.GetLifeTime(),
              normdecayLength,
              chi2topo,
              source);
    }
  }

  /// Process function for data
  void processData(soa::Filtered<CollisionTableData>::iterator const& collision, soa::Filtered<TrackTableData> const& tracks, aod::BCsWithTimestamps const&)
  {

    auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    if (runNumber != bc.runNumber()) {
      initMagneticFieldCCDB(bc, runNumber, ccdb, isRun3 ? ccdbPathGrpMag : ccdbPathGrp, lut, isRun3);
      magneticField = o2::base::Propagator::Instance()->getNominalBz();
/// Set magnetic field for KF vertexing
#ifdef HomogeneousField
      KFParticle::SetField(magneticField);
#endif
    }
    /// set KF primary vertex
    KFPVertex kfpVertex = createKFPVertexFromCollision(collision);
    KFParticle KFPV(kfpVertex);

    for (auto& [track1, track2, track3] : combinations(soa::CombinationsStrictlyUpperIndexPolicy(tracks, tracks, tracks))) {

      /// Apply single track selection
      if (!isSelectedTracks(track1, track2, track3)) {
        continue;
      }

      bool CandLc = false;
      bool CandLcbar = false;
      int source = 0;

      if (isSelectedLc(track1, track2, track3)) {
        CandLc = true;
        source = 1;
        bool LcReconstructed = ReconstructLc(track1, track2, track3, KFPV);
        if (LcReconstructed) {
          writeVarTree(KFPion, KFKaon, KFProton, KFLc, KFLc_PV, KFPV, track1, track2, track3, source);
        }
      }
      if (isSelectedLc(track1, track3, track2)) {
        CandLc = true;
        source = 1;
        bool LcReconstructed = ReconstructLc(track1, track3, track2, KFPV);
        if (LcReconstructed) {
          writeVarTree(KFPion, KFKaon, KFProton, KFLc, KFLc_PV, KFPV, track1, track3, track2, source);
        }
      }
      if (isSelectedLc(track2, track1, track3)) {
        CandLc = true;
        source = 1;
        bool LcReconstructed = ReconstructLc(track2, track1, track3, KFPV);
        if (LcReconstructed) {
          writeVarTree(KFPion, KFKaon, KFProton, KFLc, KFLc_PV, KFPV, track2, track1, track3, source);
        }
      }
      if (isSelectedLc(track2, track3, track1)) {
        CandLc = true;
        source = 1;
        bool LcReconstructed = ReconstructLc(track2, track3, track1, KFPV);
        if (LcReconstructed) {
          writeVarTree(KFPion, KFKaon, KFProton, KFLc, KFLc_PV, KFPV, track2, track3, track1, source);
        }
      }
      if (isSelectedLc(track3, track1, track2)) {
        CandLc = true;
        source = 1;
        bool LcReconstructed = ReconstructLc(track3, track1, track2, KFPV);
        if (LcReconstructed) {
          writeVarTree(KFPion, KFKaon, KFProton, KFLc, KFLc_PV, KFPV, track3, track1, track2, source);
        }
      }
      if (isSelectedLc(track3, track2, track1)) {
        CandLc = true;
        source = 1;
        bool LcReconstructed = ReconstructLc(track3, track2, track1, KFPV);
        if (LcReconstructed) {
          writeVarTree(KFPion, KFKaon, KFProton, KFLc, KFLc_PV, KFPV, track3, track2, track1, source);
        }
      }
      if (isSelectedLcBar(track1, track2, track3)) {
        CandLcbar = true;
        source = 2;
        bool LcReconstructed = ReconstructLc(track1, track2, track3, KFPV);
        if (LcReconstructed) {
          writeVarTree(KFPion, KFKaon, KFProton, KFLc, KFLc_PV, KFPV, track1, track2, track3, source);
        }
      }
      if (isSelectedLcBar(track1, track3, track2)) {
        CandLcbar = true;
        source = 2;
        bool LcReconstructed = ReconstructLc(track1, track3, track2, KFPV);
        if (LcReconstructed) {
          writeVarTree(KFPion, KFKaon, KFProton, KFLc, KFLc_PV, KFPV, track1, track3, track2, source);
        }
      }
      if (isSelectedLcBar(track2, track1, track3)) {
        CandLcbar = true;
        source = 2;
        bool LcReconstructed = ReconstructLc(track2, track1, track3, KFPV);
        if (LcReconstructed) {
          writeVarTree(KFPion, KFKaon, KFProton, KFLc, KFLc_PV, KFPV, track2, track1, track3, source);
        }
      }
      if (isSelectedLcBar(track2, track3, track1)) {
        CandLcbar = true;
        source = 2;
        bool LcReconstructed = ReconstructLc(track2, track3, track1, KFPV);
        if (LcReconstructed) {
          writeVarTree(KFPion, KFKaon, KFProton, KFLc, KFLc_PV, KFPV, track2, track3, track1, source);
        }
      }
      if (isSelectedLcBar(track3, track1, track2)) {
        CandLcbar = true;
        source = 2;
        bool LcReconstructed = ReconstructLc(track3, track1, track2, KFPV);
        if (LcReconstructed) {
          writeVarTree(KFPion, KFKaon, KFProton, KFLc, KFLc_PV, KFPV, track3, track1, track2, source);
        }
      }
      if (isSelectedLcBar(track3, track2, track1)) {
        CandLcbar = true;
        source = 2;
        bool LcReconstructed = ReconstructLc(track3, track2, track1, KFPV);
        if (LcReconstructed) {
          writeVarTree(KFPion, KFKaon, KFProton, KFLc, KFLc_PV, KFPV, track3, track2, track1, source);
        }
      }
      if (!CandLc && !CandLcbar) {
        continue;
      }
    }
  }
  PROCESS_SWITCH(qaKFParticleLc, processData, "process data", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<qaKFParticleLc>(cfgc)};
}
