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
// Lambda polarisation task
// prottay.das@cern.ch

#include <TH1F.h>
#include <TDirectory.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <cmath>
#include <array>
#include <cstdlib>
#include <algorithm>
#include "Common/Core/RecoDecay.h"
#include <string>

#include "TRandom3.h"
#include "Math/Vector3D.h"
#include "Math/Vector4D.h"
#include "Math/GenVector/Boost.h"
#include "TF1.h"

#include "PWGLF/DataModel/SPCalibrationTables.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"
#include "PWGLF/DataModel/LFStrangenessTables.h"
#include "PWGLF/DataModel/LFStrangenessPIDTables.h"
#include "Common/DataModel/FT0Corrected.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;

using dauTracks = soa::Join<aod::DauTrackExtras, aod::DauTrackTPCPIDs>;
using v0Candidates = soa::Join<aod::V0CollRefs, aod::V0Cores, aod::V0Extras>;

struct lambdapolsp {

  int mRunNumber;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  // fill output
  Configurable<bool> additionalEvSel{"additionalEvSel", false, "additionalEvSel"};
  Configurable<bool> additionalEvSel2{"additionalEvSel2", false, "additionalEvSel2"};
  Configurable<bool> additionalEvSel3{"additionalEvSel3", false, "additionalEvSel3"};
  Configurable<bool> additionalEvSel4{"additionalEvSel4", false, "additionalEvSel4"};
  Configurable<bool> globalpt{"globalpt", true, "select tracks based on pt global vs tpc"};
  Configurable<bool> cqvas{"cqvas", false, "change q vectors after shift correction"};
  Configurable<int> useprofile{"useprofile", 3, "flag to select profile vs Sparse"};
  Configurable<int> QxyNbins{"QxyNbins", 100, "Number of bins in QxQy histograms"};
  Configurable<float> lbinQxy{"lbinQxy", -5.0, "lower bin value in QxQy histograms"};
  Configurable<float> hbinQxy{"hbinQxy", 5.0, "higher bin value in QxQy histograms"};
  Configurable<int> cfgMaxOccupancy{"cfgMaxOccupancy", 1000, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> cfgMinOccupancy{"cfgMinOccupancy", 0, "maximum occupancy of tracks in neighbouring collisions in a given time range"};
  Configurable<int> sys{"sys", 1, "flag to select systematic source"};
  Configurable<bool> dosystematic{"dosystematic", false, "flag to perform systematic study"};
  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutCentralityMax{"cfgCutCentralityMax", 50.0f, "Accepted maximum Centrality"};
  Configurable<float> cfgCutCentralityMin{"cfgCutCentralityMin", 30.0f, "Accepted minimum Centrality"};
  // proton track cut
  Configurable<float> cfgCutPT{"cfgCutPT", 0.15, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 0.1f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 0.1f, "DCAz range for tracks"};
  Configurable<int> cfgITScluster{"cfgITScluster", 5, "Number of ITS cluster"};
  Configurable<int> cfgTPCcluster{"cfgTPCcluster", 70, "Number of TPC cluster"};
  Configurable<bool> isPVContributor{"isPVContributor", true, "is PV contributor"};
  Configurable<bool> checkwithpub{"checkwithpub", true, "checking results with published"};

  // Configs for V0
  Configurable<float> ConfV0PtMin{"ConfV0PtMin", 0.f, "Minimum transverse momentum of V0"};
  Configurable<float> ConfV0Rap{"ConfV0Rap", 0.8f, "Rapidity range of V0"};
  Configurable<double> ConfV0DCADaughMax{"ConfV0DCADaughMax", 0.2f, "Maximum DCA between the V0 daughters"};
  Configurable<double> ConfV0CPAMin{"ConfV0CPAMin", 0.9998f, "Minimum CPA of V0"};
  Configurable<float> ConfV0TranRadV0Min{"ConfV0TranRadV0Min", 1.5f, "Minimum transverse radius"};
  Configurable<float> ConfV0TranRadV0Max{"ConfV0TranRadV0Max", 100.f, "Maximum transverse radius"};
  Configurable<double> cMaxV0DCA{"cMaxV0DCA", 1.2, "Maximum V0 DCA to PV"};
  Configurable<double> cMinV0DCAPr{"cMinV0DCAPr", 0.05, "Minimum V0 daughters DCA to PV for Pr"};
  Configurable<double> cMinV0DCAPi{"cMinV0DCAPi", 0.05, "Minimum V0 daughters DCA to PV for Pi"};
  Configurable<float> cMaxV0LifeTime{"cMaxV0LifeTime", 20, "Maximum V0 life time"};

  // config for V0 daughters
  Configurable<float> ConfDaughEta{"ConfDaughEta", 0.8f, "V0 Daugh sel: max eta"};
  Configurable<float> cfgDaughPrPt{"cfgDaughPrPt", 0.4, "minimum daughter proton pt"};
  Configurable<float> cfgDaughPiPt{"cfgDaughPiPt", 0.2, "minimum daughter pion pt"};
  Configurable<float> ConfDaughTPCnclsMin{"ConfDaughTPCnclsMin", 50.f, "V0 Daugh sel: Min. nCls TPC"};
  Configurable<float> ConfDaughPIDCuts{"ConfDaughPIDCuts", 3, "PID selections for Lambda daughters"};
  Configurable<bool> usesubdet{"usesubdet", false, "use subdet"};
  Configurable<bool> useAccCorr{"useAccCorr", false, "use acceptance correction"};
  Configurable<std::string> ConfAccPathL{"ConfAccPathL", "Users/p/prottay/My/Object/From379780/Fulldata/NewPbPbpass4_28032025/acccorrL", "Path to acceptance correction for Lambda"};
  Configurable<std::string> ConfAccPathAL{"ConfAccPathAL", "Users/p/prottay/My/Object/From379780/Fulldata/NewPbPbpass4_28032025/acccorrAL", "Path to acceptance correction for AntiLambda"};

  Configurable<int> CentNbins{"CentNbins", 16, "Number of bins in cent histograms"};
  Configurable<float> lbinCent{"lbinCent", 0.0, "lower bin value in cent histograms"};
  Configurable<float> hbinCent{"hbinCent", 80.0, "higher bin value in cent histograms"};
  Configurable<int> SANbins{"SANbins", 20, "Number of bins in costhetastar"};
  Configurable<float> lbinSA{"lbinSA", -1.0, "lower bin value in costhetastar histograms"};
  Configurable<float> hbinSA{"hbinSA", 1.0, "higher bin value in costhetastar histograms"};
  Configurable<int> PolNbins{"PolNbins", 20, "Number of bins in polarisation"};
  Configurable<float> lbinPol{"lbinPol", -1.0, "lower bin value in #phi-#psi histograms"};
  Configurable<float> hbinPol{"hbinPol", 1.0, "higher bin value in #phi-#psi histograms"};
  Configurable<int> IMNbins{"IMNbins", 100, "Number of bins in invariant mass"};
  Configurable<float> lbinIM{"lbinIM", 1.0, "lower bin value in IM histograms"};
  Configurable<float> hbinIM{"hbinIM", 1.2, "higher bin value in IM histograms"};
  Configurable<int> resNbins{"resNbins", 50, "Number of bins in reso"};
  Configurable<float> lbinres{"lbinres", 0.0, "lower bin value in reso histograms"};
  Configurable<float> hbinres{"hbinres", 10.0, "higher bin value in reso histograms"};
  Configurable<int> spNbins{"spNbins", 2000, "Number of bins in sp"};
  Configurable<float> lbinsp{"lbinsp", -1.0, "lower bin value in sp histograms"};
  Configurable<float> hbinsp{"hbinsp", 1.0, "higher bin value in sp histograms"};

  ConfigurableAxis configcentAxis{"configcentAxis", {VARIABLE_WIDTH, 0.0, 10.0, 40.0, 80.0}, "Cent V0M"};
  ConfigurableAxis configthnAxispT{"configthnAxisPt", {VARIABLE_WIDTH, 0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.5, 8.0, 10.0, 100.0}, "#it{p}_{T} (GeV/#it{c})"};
  ConfigurableAxis configetaAxis{"configetaAxis", {VARIABLE_WIDTH, -0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.8}, "Eta"};
  ConfigurableAxis configthnAxisPol{"configthnAxisPol", {VARIABLE_WIDTH, -1.0, -0.6, -0.2, 0, 0.2, 0.4, 0.8}, "Pol"};
  ConfigurableAxis configbinAxis{"configbinAxis", {VARIABLE_WIDTH, -0.8, -0.4, -0.2, 0, 0.2, 0.4, 0.8}, "BA"};
  // ConfigurableAxis configphiAxis{"configphiAxis", {VARIABLE_WIDTH, 0.0, 0.2, 0.4, 0.8, 1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 5.5, 6.28}, "PhiAxis"};

  SliceCache cache;
  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext&)
  {
    AxisSpec thnAxisres{resNbins, lbinres, hbinres, "Reso"};
    AxisSpec thnAxisInvMass{IMNbins, lbinIM, hbinIM, "#it{M} (GeV/#it{c}^{2})"};
    AxisSpec thnAxisCosThetaStar{SANbins, lbinSA, hbinSA, "SA"};
    AxisSpec centAxis = {CentNbins, lbinCent, hbinCent, "V0M (%)"};
    AxisSpec spAxis = {spNbins, lbinsp, hbinsp, "Sp"};
    AxisSpec qxZDCAxis = {QxyNbins, lbinQxy, hbinQxy, "Qx"};

    if (checkwithpub) {
      if (useprofile == 2) {
        histos.add("hpuxQxpvscentpteta", "hpuxQxpvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuyQypvscentpteta", "hpuyQypvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuxQxtvscentpteta", "hpuxQxtvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuyQytvscentpteta", "hpuyQytvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuxyQxytvscentpteta", "hpuxyQxytvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuxyQxypvscentpteta", "hpuxyQxypvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpoddv1vscentpteta", "hpoddv1vscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpevenv1vscentpteta", "hpevenv1vscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpv21", "hpv21", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpv22", "hpv22", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpv23", "hpv23", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpx2Tx1Ax1Cvscentpteta", "hpx2Tx1Ax1Cvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpx2Ty1Ay1Cvscentpteta", "hpx2Ty1Ay1Cvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpy2Tx1Ay1Cvscentpteta", "hpy2Tx1Ay1Cvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpy2Ty1Ax1Cvscentpteta", "hpy2Ty1Ax1Cvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpx1Ax1Cvscentpteta", "hpx1Ax1Cvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpy1Ay1Cvscentpteta", "hpy1Ay1Cvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpx1Avscentpteta", "hpx1Avscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpx1Cvscentpteta", "hpx1Cvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpy1Avscentpteta", "hpy1Avscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpy1Cvscentpteta", "hpy1Cvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);

        histos.add("hpx2Tx1Avscentpteta", "hpx2Tx1Avscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpx2Tx1Cvscentpteta", "hpx2Tx1Cvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpx2Ty1Avscentpteta", "hpx2Ty1Avscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpx2Ty1Cvscentpteta", "hpx2Ty1Cvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpy2Tx1Avscentpteta", "hpy2Tx1Avscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpy2Ty1Cvscentpteta", "hpy2Ty1Cvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpy2Ty1Avscentpteta", "hpy2Ty1Avscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpy2Tx1Cvscentpteta", "hpy2Tx1Cvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpx1Ay1Cvscentpteta", "hpx1Ay1Cvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpy1Ax1Cvscentpteta", "hpy1Ax1Cvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpx2Tvscentpteta", "hpx2Tvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpy2Tvscentpteta", "hpy2Tvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);

        histos.add("hpuxvscentpteta", "hpuxvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuyvscentpteta", "hpuyvscentpteta", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        /*
              histos.add("hpuxvscentptetaneg", "hpuxvscentptetaneg", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
              histos.add("hpuyvscentptetaneg", "hpuyvscentptetaneg", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);

              histos.add("hpuxQxpvscentptetaneg", "hpuxQxpvscentptetaneg", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
              histos.add("hpuyQypvscentptetaneg", "hpuyQypvscentptetaneg", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
              histos.add("hpuxQxtvscentptetaneg", "hpuxQxtvscentptetaneg", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
              histos.add("hpuyQytvscentptetaneg", "hpuyQytvscentptetaneg", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
              histos.add("hpuxyQxytvscentptetaneg", "hpuxyQxytvscentptetaneg", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
              histos.add("hpuxyQxypvscentptetaneg", "hpuxyQxypvscentptetaneg", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
              histos.add("hpoddv1vscentptetaneg", "hpoddv1vscentptetaneg", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
              histos.add("hpevenv1vscentptetaneg", "hpevenv1vscentptetaneg", HistType::kTHnSparseF, {centAxis, configthnAxispT, configetaAxis, spAxis}, true);
        */

        histos.add("hpQxtQxpvscent", "hpQxtQxpvscent", HistType::kTHnSparseF, {centAxis, spAxis}, true);
        histos.add("hpQytQypvscent", "hpQytQypvscent", HistType::kTHnSparseF, {centAxis, spAxis}, true);
        histos.add("hpQxytpvscent", "hpQxytpvscent", HistType::kTHnSparseF, {centAxis, spAxis}, true);
        histos.add("hpQxtQypvscent", "hpQxtQypvscent", HistType::kTHnSparseF, {centAxis, spAxis}, true);
        histos.add("hpQxpQytvscent", "hpQxpQytvscent", HistType::kTHnSparseF, {centAxis, spAxis}, true);

        histos.add("hpQxpvscent", "hpQxpvscent", HistType::kTHnSparseF, {centAxis, spAxis}, true);
        histos.add("hpQxtvscent", "hpQxtvscent", HistType::kTHnSparseF, {centAxis, spAxis}, true);
        histos.add("hpQypvscent", "hpQypvscent", HistType::kTHnSparseF, {centAxis, spAxis}, true);
        histos.add("hpQytvscent", "hpQytvscent", HistType::kTHnSparseF, {centAxis, spAxis}, true);
      } else {
        histos.add("hpuxQxpvscentpteta", "hpuxQxpvscentpteta", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuyQypvscentpteta", "hpuyQypvscentpteta", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuxQxtvscentpteta", "hpuxQxtvscentpteta", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuyQytvscentpteta", "hpuyQytvscentpteta", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuxyQxytvscentpteta", "hpuxyQxytvscentpteta", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuxyQxypvscentpteta", "hpuxyQxypvscentpteta", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpoddv1vscentpteta", "hpoddv1vscentpteta", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpevenv1vscentpteta", "hpevenv1vscentpteta", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);

        histos.add("hpuxvscentpteta", "hpuxvscentpteta", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuyvscentpteta", "hpuyvscentpteta", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        /*histos.add("hpuxvscentptetaneg", "hpuxvscentptetaneg", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuyvscentptetaneg", "hpuyvscentptetaneg", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);

        histos.add("hpuxQxpvscentptetaneg", "hpuxQxpvscentptetaneg", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuyQypvscentptetaneg", "hpuyQypvscentptetaneg", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuxQxtvscentptetaneg", "hpuxQxtvscentptetaneg", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuyQytvscentptetaneg", "hpuyQytvscentptetaneg", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuxyQxytvscentptetaneg", "hpuxyQxytvscentptetaneg", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpuxyQxypvscentptetaneg", "hpuxyQxypvscentptetaneg", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpoddv1vscentptetaneg", "hpoddv1vscentptetaneg", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);
        histos.add("hpevenv1vscentptetaneg", "hpevenv1vscentptetaneg", HistType::kTHnSparseF, {configcentAxis, configthnAxispT, configetaAxis, spAxis}, true);*/

        histos.add("hpQxtQxpvscent", "hpQxtQxpvscent", HistType::kTHnSparseF, {configcentAxis, spAxis}, true);
        histos.add("hpQytQypvscent", "hpQytQypvscent", HistType::kTHnSparseF, {configcentAxis, spAxis}, true);
        histos.add("hpQxytpvscent", "hpQxytpvscent", HistType::kTHnSparseF, {configcentAxis, spAxis}, true);
        histos.add("hpQxtQypvscent", "hpQxtQypvscent", HistType::kTHnSparseF, {configcentAxis, spAxis}, true);
        histos.add("hpQxpQytvscent", "hpQxpQytvscent", HistType::kTHnSparseF, {configcentAxis, spAxis}, true);

        histos.add("hpQxpvscent", "hpQxpvscent", HistType::kTHnSparseF, {configcentAxis, spAxis}, true);
        histos.add("hpQxtvscent", "hpQxtvscent", HistType::kTHnSparseF, {configcentAxis, spAxis}, true);
        histos.add("hpQypvscent", "hpQypvscent", HistType::kTHnSparseF, {configcentAxis, spAxis}, true);
        histos.add("hpQytvscent", "hpQytvscent", HistType::kTHnSparseF, {configcentAxis, spAxis}, true);
      }
    }

    histos.add("hCentrality", "Centrality distribution", kTH1F, {{centAxis}});
    // histos.add("hpsiApsiC", "hpsiApsiC", kTHnSparseF, {psiACAxis, psiACAxis});
    //  histos.add("hpsiApsiC", "hpsiApsiC", kTH2F, {psiACAxis, psiACAxis});
    // histos.add("hphiminuspsiA", "hphiminuspisA", kTH1F, {{50, 0, 6.28}}, true);
    // histos.add("hphiminuspsiC", "hphiminuspisC", kTH1F, {{50, 0, 6.28}}, true);
    //  histos.add("hCentrality0", "Centrality distribution0", kTH1F, {{centAxis}});
    //  histos.add("hCentrality1", "Centrality distribution1", kTH1F, {{centAxis}});
    //  histos.add("hCentrality2", "Centrality distribution2", kTH1F, {{centAxis}});
    //  histos.add("hCentrality3", "Centrality distribution3", kTH1F, {{centAxis}});

    if (!checkwithpub) {
      // histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{20, -10.0, 10.0}});
      histos.add("hpRes", "hpRes", HistType::kTHnSparseF, {configcentAxis, thnAxisres});
      histos.add("hpResSin", "hpResSin", HistType::kTHnSparseF, {configcentAxis, thnAxisres});
      /*histos.add("hpCosPsiA", "hpCosPsiA", HistType::kTHnSparseF, {configcentAxis, thnAxisres});
      histos.add("hpCosPsiC", "hpCosPsiC", HistType::kTHnSparseF, {configcentAxis, thnAxisres});
      histos.add("hpSinPsiA", "hpSinPsiA", HistType::kTHnSparseF, {configcentAxis, thnAxisres});
      histos.add("hpSinPsiC", "hpSinPsiC", HistType::kTHnSparseF, {configcentAxis, thnAxisres});*/
      /*histos.add("hcentQxZDCA", "hcentQxZDCA", kTH2F, {{centAxis}, {qxZDCAxis}});
      histos.add("hcentQyZDCA", "hcentQyZDCA", kTH2F, {{centAxis}, {qxZDCAxis}});
      histos.add("hcentQxZDCC", "hcentQxZDCC", kTH2F, {{centAxis}, {qxZDCAxis}});
      histos.add("hcentQyZDCC", "hcentQyZDCC", kTH2F, {{centAxis}, {qxZDCAxis}});*/

      if (usesubdet) {
        histos.add("hSparseLambdaCosPsiA", "hSparseLambdaCosPsiA", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
        histos.add("hSparseLambdaSinPsiA", "hSparseLambdaSinPsiA", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
        histos.add("hSparseLambdaCosPsiC", "hSparseLambdaCosPsiC", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
        histos.add("hSparseLambdaSinPsiC", "hSparseLambdaSinPsiC", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
      }
      histos.add("hSparseLambdaCosPsi", "hSparseLambdaCosPsi", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
      histos.add("hSparseLambdaSinPsi", "hSparseLambdaSinPsi", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
      if (usesubdet) {
        histos.add("hSparseAntiLambdaCosPsiA", "hSparseAntiLambdaCosPsiA", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
        histos.add("hSparseAntiLambdaSinPsiA", "hSparseAntiLambdaSinPsiA", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
        histos.add("hSparseAntiLambdaCosPsiC", "hSparseAntiLambdaCosPsiC", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
        histos.add("hSparseAntiLambdaSinPsiC", "hSparseAntiLambdaSinPsiC", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
      }
      histos.add("hSparseAntiLambdaCosPsi", "hSparseAntiLambdaCosPsi", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
      histos.add("hSparseAntiLambdaSinPsi", "hSparseAntiLambdaSinPsi", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);

      histos.add("hSparseLambdaPol", "hSparseLambdaPol", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
      if (usesubdet) {
        histos.add("hSparseLambdaPolA", "hSparseLambdaPolA", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
        histos.add("hSparseLambdaPolC", "hSparseLambdaPolC", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
      }
      histos.add("hSparseAntiLambdaPol", "hSparseAntiLambdaPol", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
      if (usesubdet) {
        histos.add("hSparseAntiLambdaPolA", "hSparseAntiLambdaPolA", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
        histos.add("hSparseAntiLambdaPolC", "hSparseAntiLambdaPolC", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
      }
      histos.add("hSparseLambda_corr1a", "hSparseLambda_corr1a", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
      histos.add("hSparseLambda_corr1b", "hSparseLambda_corr1b", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
      // histos.add("hSparseLambda_corr1c", "hSparseLambda_corr1c", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configphiAxis, configcentAxis, configbinAxis}, true);
      histos.add("hSparseAntiLambda_corr1a", "hSparseAntiLambda_corr1a", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
      histos.add("hSparseAntiLambda_corr1b", "hSparseAntiLambda_corr1b", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
      // histos.add("hSparseAntiLambda_corr1c", "hSparseAntiLambda_corr1c", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configphiAxis, configcentAxis, configbinAxis}, true);

      histos.add("hSparseLambda_corr2a", "hSparseLambda_corr2a", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
      // histos.add("hSparseLambda_corr2b", "hSparseLambda_corr2b", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
      histos.add("hSparseAntiLambda_corr2a", "hSparseAntiLambda_corr2a", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
      // histos.add("hSparseAntiLambda_corr2b", "hSparseAntiLambda_corr2b", HistType::kTHnSparseF, {thnAxisInvMass, configthnAxispT, configthnAxisPol, configcentAxis, configbinAxis}, true);
    }
  }

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (!(candidate.isGlobalTrack() && candidate.isPVContributor() && candidate.itsNCls() > cfgITScluster && candidate.tpcNClsFound() > cfgTPCcluster && candidate.itsNClsInnerBarrel() >= 1)) {
      return false;
    }
    return true;
  }

  template <typename Collision, typename V0>
  bool SelectionV0(Collision const& collision, V0 const& candidate)
  {
    if (TMath::Abs(candidate.dcav0topv()) > cMaxV0DCA) {
      return false;
    }
    const float pT = candidate.pt();
    const float tranRad = candidate.v0radius();
    const float dcaDaughv0 = TMath::Abs(candidate.dcaV0daughters());
    const float cpav0 = candidate.v0cosPA();

    float CtauLambda = candidate.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massLambda;
    // float lowmasscutlambda = cMinLambdaMass;
    // float highmasscutlambda = cMaxLambdaMass;

    if (pT < ConfV0PtMin) {
      return false;
    }
    if (dcaDaughv0 > ConfV0DCADaughMax) {
      return false;
    }
    if (cpav0 < ConfV0CPAMin) {
      return false;
    }
    if (tranRad < ConfV0TranRadV0Min) {
      return false;
    }
    if (tranRad > ConfV0TranRadV0Max) {
      return false;
    }
    if (TMath::Abs(CtauLambda) > cMaxV0LifeTime) {
      return false;
    }
    if (TMath::Abs(candidate.yLambda()) > ConfV0Rap) {
      return false;
    }
    return true;
  }
  template <typename V0, typename T>
  bool isSelectedV0Daughter(V0 const& candidate, T const& track, int pid)
  {
    // const auto eta = track.eta();
    // const auto pt = track.pt();
    const auto tpcNClsF = track.tpcNClsFound();
    if (track.tpcNClsCrossedRows() < 70) {
      return false;
    }
    /*if (TMath::Abs(eta) > ConfDaughEta) {
      return false;
      }*/
    if (tpcNClsF < ConfDaughTPCnclsMin) {
      return false;
    }
    if (track.tpcCrossedRowsOverFindableCls() < 0.8) {
      return false;
    }

    if (pid == 0 && TMath::Abs(track.tpcNSigmaPr()) > ConfDaughPIDCuts) {
      return false;
    }
    if (pid == 1 && TMath::Abs(track.tpcNSigmaPi()) > ConfDaughPIDCuts) {
      return false;
    }
    if (pid == 0 && (candidate.positivept() < cfgDaughPrPt || candidate.negativept() < cfgDaughPiPt)) {
      return false; // doesn´t pass lambda pT sels
    }
    if (pid == 1 && (candidate.positivept() < cfgDaughPiPt || candidate.negativept() < cfgDaughPrPt)) {
      return false; // doesn´t pass antilambda pT sels
    }
    if (std::abs(candidate.positiveeta()) > ConfDaughEta || std::abs(candidate.negativeeta()) > ConfDaughEta) {
      return false;
    }

    if (pid == 0 && (TMath::Abs(candidate.dcapostopv()) < cMinV0DCAPr || TMath::Abs(candidate.dcanegtopv()) < cMinV0DCAPi)) {
      return false;
    }
    if (pid == 1 && (TMath::Abs(candidate.dcapostopv()) < cMinV0DCAPi || TMath::Abs(candidate.dcanegtopv()) < cMinV0DCAPr)) {
      return false;
    }

    return true;
  }

  template <typename TV0>
  bool isCompatible(TV0 const& v0, int pid /*0: lambda, 1: antilambda*/)
  {
    // checks if this V0 is compatible with the requested hypothesis

    // de-ref track extras
    auto posTrackExtra = v0.template posTrackExtra_as<dauTracks>();
    auto negTrackExtra = v0.template negTrackExtra_as<dauTracks>();

    // check for desired kinematics
    if (pid == 0 && (v0.positivept() < cfgDaughPrPt || v0.negativept() < cfgDaughPiPt)) {
      return false; // doesn´t pass lambda pT sels
    }
    if (pid == 1 && (v0.positivept() < cfgDaughPiPt || v0.negativept() < cfgDaughPrPt)) {
      return false; // doesn´t pass antilambda pT sels
    }
    if (std::abs(v0.positiveeta()) > ConfDaughEta || std::abs(v0.negativeeta()) > ConfDaughEta) {
      return false;
    }

    // check TPC tracking properties
    if (posTrackExtra.tpcNClsCrossedRows() < 70 || negTrackExtra.tpcNClsCrossedRows() < 70) {
      return false;
    }
    if (posTrackExtra.tpcNClsFound() < ConfDaughTPCnclsMin || negTrackExtra.tpcNClsFound() < ConfDaughTPCnclsMin) {
      return false;
    }
    if (posTrackExtra.tpcCrossedRowsOverFindableCls() < 0.8 || negTrackExtra.tpcCrossedRowsOverFindableCls() < 0.8) {
      return false;
    }

    // check TPC PID
    if (pid == 0 && ((std::abs(posTrackExtra.tpcNSigmaPr()) > ConfDaughPIDCuts) || (std::abs(negTrackExtra.tpcNSigmaPi()) > ConfDaughPIDCuts))) {
      return false;
    }
    if (pid == 1 && ((std::abs(posTrackExtra.tpcNSigmaPi()) > ConfDaughPIDCuts) || (std::abs(negTrackExtra.tpcNSigmaPr()) > ConfDaughPIDCuts))) {
      return false;
    }

    if (pid == 0 && (TMath::Abs(v0.dcapostopv()) < cMinV0DCAPr || TMath::Abs(v0.dcanegtopv()) < cMinV0DCAPi)) {
      return false;
    }
    if (pid == 1 && (TMath::Abs(v0.dcapostopv()) < cMinV0DCAPi || TMath::Abs(v0.dcanegtopv()) < cMinV0DCAPr)) {
      return false;
    }

    // if we made it this far, it's good
    return true;
  }

  double GetPhiInRange(double phi)
  {
    double result = RecoDecay::constrainAngle(phi);

    /*
      double result = phi;
      while (result < 0) {
      // result = result + 2. * TMath::Pi();
      result = result + 2. * o2::constants::math::PI;
      }
      while (result > 2. * TMath::Pi()) {
      // result = result - 2. * TMath::Pi();
      result = result - 2. * o2::constants::math::PI;
      }*/
    return result;
  }

  bool shouldReject(bool LambdaTag, bool aLambdaTag,
                    const ROOT::Math::PxPyPzMVector& Lambdadummy,
                    const ROOT::Math::PxPyPzMVector& AntiLambdadummy)
  {
    const double minMass = 1.105;
    const double maxMass = 1.125;
    return (LambdaTag && aLambdaTag &&
            (Lambdadummy.M() > minMass && Lambdadummy.M() < maxMass) &&
            (AntiLambdadummy.M() > minMass && AntiLambdadummy.M() < maxMass));
  }

  void fillHistograms(bool tag1, bool tag2, const ROOT::Math::PxPyPzMVector& particle,
                      const ROOT::Math::PxPyPzMVector& daughter,
                      double psiZDCC, double psiZDCA, double psiZDC, double centrality,
                      double candmass, double candpt, float desbinvalue, double acvalue)
  {

    ROOT::Math::Boost boost{particle.BoostToCM()};
    auto fourVecDauCM = boost(daughter);
    auto phiangle = TMath::ATan2(fourVecDauCM.Py(), fourVecDauCM.Px());
    auto phiminuspsiC = GetPhiInRange(phiangle - psiZDCC);
    auto phiminuspsiA = GetPhiInRange(phiangle - psiZDCA);
    auto phiminuspsi = GetPhiInRange(phiangle - psiZDC);
    auto cosThetaStar = fourVecDauCM.Pz() / fourVecDauCM.P();
    auto sinThetaStar = TMath::Sqrt(1 - (cosThetaStar * cosThetaStar));
    auto PolC = TMath::Sin(phiminuspsiC);
    auto PolA = TMath::Sin(phiminuspsiA);
    auto Pol = TMath::Sin(phiminuspsi);

    auto sinPhiStar = TMath::Sin(GetPhiInRange(phiangle));
    auto cosPhiStar = TMath::Cos(GetPhiInRange(phiangle));
    // auto sinThetaStarcosphiphiStar = sinThetaStar * TMath::Cos(2 * GetPhiInRange(particle.Phi() - phiangle));
    // auto phiphiStar = GetPhiInRange(particle.Phi() - phiangle);

    if (useAccCorr) {
      acvalue = (4 / 3.14) * acvalue;
      PolC = PolC / acvalue;
      PolA = PolA / acvalue;
      Pol = Pol / acvalue;
    }

    // Fill histograms using constructed names
    if (tag2) {
      if (usesubdet) {
        histos.fill(HIST("hSparseAntiLambdaCosPsiA"), candmass, candpt, (TMath::Cos(GetPhiInRange(psiZDCA))), centrality, desbinvalue);
        histos.fill(HIST("hSparseAntiLambdaCosPsiC"), candmass, candpt, (TMath::Cos(GetPhiInRange(psiZDCC))), centrality, desbinvalue);
        histos.fill(HIST("hSparseAntiLambdaSinPsiA"), candmass, candpt, (TMath::Sin(GetPhiInRange(psiZDCA))), centrality, desbinvalue);
        histos.fill(HIST("hSparseAntiLambdaSinPsiC"), candmass, candpt, (TMath::Sin(GetPhiInRange(psiZDCC))), centrality, desbinvalue);
      }
      histos.fill(HIST("hSparseAntiLambdaCosPsi"), candmass, candpt, (TMath::Cos(GetPhiInRange(psiZDC))), centrality, desbinvalue);
      histos.fill(HIST("hSparseAntiLambdaSinPsi"), candmass, candpt, (TMath::Sin(GetPhiInRange(psiZDC))), centrality, desbinvalue);

      if (usesubdet) {
        histos.fill(HIST("hSparseAntiLambdaPolA"), candmass, candpt, PolA, centrality, desbinvalue);
        histos.fill(HIST("hSparseAntiLambdaPolC"), candmass, candpt, PolC, centrality, desbinvalue);
      }
      histos.fill(HIST("hSparseAntiLambdaPol"), candmass, candpt, Pol, centrality, desbinvalue);
      histos.fill(HIST("hSparseAntiLambda_corr1a"), candmass, candpt, sinPhiStar, centrality, desbinvalue);
      histos.fill(HIST("hSparseAntiLambda_corr1b"), candmass, candpt, cosPhiStar, centrality, desbinvalue);
      // histos.fill(HIST("hSparseAntiLambda_corr1c"), candmass, candpt, phiphiStar, centrality, desbinvalue);
      histos.fill(HIST("hSparseAntiLambda_corr2a"), candmass, candpt, sinThetaStar, centrality, desbinvalue);
      // histos.fill(HIST("hSparseAntiLambda_corr2b"), candmass, candpt, sinThetaStarcosphiphiStar, centrality, desbinvalue);
    }
    if (tag1) {
      if (usesubdet) {
        histos.fill(HIST("hSparseLambdaCosPsiA"), candmass, candpt, (TMath::Cos(GetPhiInRange(psiZDCA))), centrality, desbinvalue);
        histos.fill(HIST("hSparseLambdaCosPsiC"), candmass, candpt, (TMath::Cos(GetPhiInRange(psiZDCC))), centrality, desbinvalue);
        histos.fill(HIST("hSparseLambdaSinPsiA"), candmass, candpt, (TMath::Sin(GetPhiInRange(psiZDCA))), centrality, desbinvalue);
        histos.fill(HIST("hSparseLambdaSinPsiC"), candmass, candpt, (TMath::Sin(GetPhiInRange(psiZDCC))), centrality, desbinvalue);
      }
      histos.fill(HIST("hSparseLambdaCosPsi"), candmass, candpt, (TMath::Cos(GetPhiInRange(psiZDC))), centrality, desbinvalue);
      histos.fill(HIST("hSparseLambdaSinPsi"), candmass, candpt, (TMath::Sin(GetPhiInRange(psiZDC))), centrality, desbinvalue);

      if (usesubdet) {
        histos.fill(HIST("hSparseLambdaPolA"), candmass, candpt, PolA, centrality, desbinvalue);
        histos.fill(HIST("hSparseLambdaPolC"), candmass, candpt, PolC, centrality, desbinvalue);
      }
      histos.fill(HIST("hSparseLambdaPol"), candmass, candpt, Pol, centrality, desbinvalue);
      histos.fill(HIST("hSparseLambda_corr1a"), candmass, candpt, sinPhiStar, centrality, desbinvalue);
      histos.fill(HIST("hSparseLambda_corr1b"), candmass, candpt, cosPhiStar, centrality, desbinvalue);
      // histos.fill(HIST("hSparseLambda_corr1c"), candmass, candpt, phiphiStar, centrality, desbinvalue);
      histos.fill(HIST("hSparseLambda_corr2a"), candmass, candpt, sinThetaStar, centrality, desbinvalue);
      // histos.fill(HIST("hSparseLambda_corr2b"), candmass, candpt, sinThetaStarcosphiphiStar, centrality, desbinvalue);
    }
  }

  ROOT::Math::PxPyPzMVector Lambda, AntiLambda, Lambdadummy, AntiLambdadummy, Proton, Pion, AntiProton, AntiPion, fourVecDauCM;
  ROOT::Math::XYZVector threeVecDauCM, threeVecDauCMXY;
  double phiangle = 0.0;
  // double angleLambda=0.0;
  // double angleAntiLambda=0.0;
  double massLambda = o2::constants::physics::MassLambda;
  double massPr = o2::constants::physics::MassProton;
  double massPi = o2::constants::physics::MassPionCharged;

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter centralityFilter = (nabs(aod::cent::centFT0C) < cfgCutCentralityMax && nabs(aod::cent::centFT0C) > cfgCutCentralityMin);
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter dcaCutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::SPCalibrationTables, aod::Mults>>;
  using AllTrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullPi, aod::pidTPCFullPr, aod::pidTPCFullKa>>;
  using ResoV0s = aod::V0Datas;

  TProfile2D* accprofileL;
  TProfile2D* accprofileAL;
  int currentRunNumber = -999;
  int lastRunNumber = -999;

  using BCsRun3 = soa::Join<aod::BCsWithTimestamps, aod::Run3MatchedToBCSparse>;

  void processData(EventCandidates::iterator const& collision, AllTrackCandidates const& tracks, ResoV0s const& V0s, BCsRun3 const&)
  {

    if (!collision.sel8()) {
      return;
    }
    auto centrality = collision.centFT0C();
    // histos.fill(HIST("hCentrality0"), centrality);
    if (!collision.triggereventsp()) {
      return;
    }
    // histos.fill(HIST("hCentrality1"), centrality);

    if (additionalEvSel && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }
    // histos.fill(HIST("hCentrality2"), centrality);
    // if (additionalEvSel2 && (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
    if (additionalEvSel2 && (collision.trackOccupancyInTimeRange() > cfgMaxOccupancy || collision.trackOccupancyInTimeRange() < cfgMinOccupancy)) {
      return;
    }
    // histos.fill(HIST("hCentrality3"), centrality);
    if (additionalEvSel3 && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }

    if (additionalEvSel4 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return;
    }

    currentRunNumber = collision.foundBC_as<BCsRun3>().runNumber();
    auto bc = collision.foundBC_as<BCsRun3>();

    auto qxZDCA = collision.qxZDCA();
    auto qxZDCC = collision.qxZDCC();
    auto qyZDCA = collision.qyZDCA();
    auto qyZDCC = collision.qyZDCC();
    auto psiZDCC = collision.psiZDCC();
    auto psiZDCA = collision.psiZDCA();

    double modqxZDCA;
    double modqyZDCA;
    double modqxZDCC;
    double modqyZDCC;

    if (cqvas) {
      modqxZDCA = TMath::Sqrt((qxZDCA * qxZDCA) + (qyZDCA * qyZDCA)) * TMath::Cos(psiZDCA);
      modqyZDCA = TMath::Sqrt((qxZDCA * qxZDCA) + (qyZDCA * qyZDCA)) * TMath::Sin(psiZDCA);
      modqxZDCC = TMath::Sqrt((qxZDCC * qxZDCC) + (qyZDCC * qyZDCC)) * TMath::Cos(psiZDCC);
      modqyZDCC = TMath::Sqrt((qxZDCC * qxZDCC) + (qyZDCC * qyZDCC)) * TMath::Sin(psiZDCC);
    } else {
      modqxZDCA = qxZDCA;
      modqyZDCA = qyZDCA;
      modqxZDCC = qxZDCC;
      modqyZDCC = qyZDCC;
    }

    auto psiZDC = TMath::ATan2((modqyZDCC - modqyZDCA), (modqxZDCC - modqxZDCA)); // full event plane
    /*if (useonlypsis) {
      psiZDC = psiZDCC - psiZDCA;
      }*/

    histos.fill(HIST("hCentrality"), centrality);
    if (!checkwithpub) {
      // histos.fill(HIST("hVtxZ"), collision.posZ());
      histos.fill(HIST("hpRes"), centrality, (TMath::Cos(GetPhiInRange(psiZDCA - psiZDCC))));
      histos.fill(HIST("hpResSin"), centrality, (TMath::Sin(GetPhiInRange(psiZDCA - psiZDCC))));
      /*histos.fill(HIST("hpCosPsiA"), centrality, (TMath::Cos(GetPhiInRange(psiZDCA))));
      histos.fill(HIST("hpCosPsiC"), centrality, (TMath::Cos(GetPhiInRange(psiZDCC))));
      histos.fill(HIST("hpSinPsiA"), centrality, (TMath::Sin(GetPhiInRange(psiZDCA))));
      histos.fill(HIST("hpSinPsiC"), centrality, (TMath::Sin(GetPhiInRange(psiZDCC))));*/
      /*histos.fill(HIST("hcentQxZDCA"), centrality, qxZDCA);
        histos.fill(HIST("hcentQyZDCA"), centrality, qyZDCA);
        histos.fill(HIST("hcentQxZDCC"), centrality, qxZDCC);
        histos.fill(HIST("hcentQyZDCC"), centrality, qyZDCC);*/
    }

    ///////////checking v1////////////////////////////////
    if (checkwithpub) {

      auto QxtQxp = modqxZDCA * modqxZDCC;
      auto QytQyp = modqyZDCA * modqyZDCC;
      auto Qxytp = QxtQxp + QytQyp;
      auto QxpQyt = modqxZDCA * modqyZDCC;
      auto QxtQyp = modqxZDCC * modqyZDCA;

      histos.fill(HIST("hpQxtQxpvscent"), centrality, QxtQxp);
      histos.fill(HIST("hpQytQypvscent"), centrality, QytQyp);
      histos.fill(HIST("hpQxytpvscent"), centrality, Qxytp);
      histos.fill(HIST("hpQxpQytvscent"), centrality, QxpQyt);
      histos.fill(HIST("hpQxtQypvscent"), centrality, QxtQyp);

      histos.fill(HIST("hpQxpvscent"), centrality, modqxZDCA);
      histos.fill(HIST("hpQxtvscent"), centrality, modqxZDCC);
      histos.fill(HIST("hpQypvscent"), centrality, modqyZDCA);
      histos.fill(HIST("hpQytvscent"), centrality, modqyZDCC);

      for (const auto& track : tracks) {
        if (!selectionTrack(track)) {
          continue;
        }

        float sign = track.sign();
        if (sign == 0.0) // removing neutral particles
          continue;

        auto ux = TMath::Cos(GetPhiInRange(track.phi()));
        auto uy = TMath::Sin(GetPhiInRange(track.phi()));
        // auto py=track.py();

        auto uxQxp = ux * modqxZDCA;
        auto uyQyp = uy * modqyZDCA;
        auto uxyQxyp = uxQxp + uyQyp;
        auto uxQxt = ux * modqxZDCC;
        auto uyQyt = uy * modqyZDCC;
        auto uxyQxyt = uxQxt + uyQyt;
        auto oddv1 = ux * (modqxZDCA - modqxZDCC) + uy * (modqyZDCA - modqyZDCC);
        auto evenv1 = ux * (modqxZDCA + modqxZDCC) + uy * (modqyZDCA + modqyZDCC);
        auto v21 = TMath::Cos(2 * (GetPhiInRange(track.phi()) - psiZDCA - psiZDCC));
        auto v22 = TMath::Cos(2 * (GetPhiInRange(track.phi()) + psiZDCA - psiZDCC));
        auto v23 = TMath::Cos(2 * (GetPhiInRange(track.phi()) - psiZDC));

        auto x2Tx1Ax1C = TMath::Cos(2 * GetPhiInRange(track.phi())) * modqxZDCA * modqxZDCC;
        auto x2Ty1Ay1C = TMath::Cos(2 * GetPhiInRange(track.phi())) * modqyZDCA * modqyZDCC;
        auto y2Tx1Ay1C = TMath::Sin(2 * GetPhiInRange(track.phi())) * modqxZDCA * modqyZDCC;
        auto y2Ty1Ax1C = TMath::Sin(2 * GetPhiInRange(track.phi())) * modqyZDCA * modqxZDCC;
        auto x1Ax1C = modqxZDCA * modqxZDCC;
        auto y1Ay1C = modqyZDCA * modqyZDCC;
        auto x1Ay1C = modqxZDCA * modqyZDCC;
        auto x1Cy1A = modqxZDCC * modqyZDCA;

        // detector acceptance corrections to match v2{ZDC}
        auto x1A = modqxZDCA;
        auto x1C = modqxZDCC;
        auto y1A = modqyZDCA;
        auto y1C = modqyZDCC;
        auto x2T = TMath::Cos(2 * GetPhiInRange(track.phi()));
        auto y2T = TMath::Sin(2 * GetPhiInRange(track.phi()));
        auto x2Tx1A = TMath::Cos(2 * GetPhiInRange(track.phi())) * modqxZDCA;
        auto x2Tx1C = TMath::Cos(2 * GetPhiInRange(track.phi())) * modqxZDCC;
        auto x2Ty1A = TMath::Cos(2 * GetPhiInRange(track.phi())) * modqyZDCA;
        auto x2Ty1C = TMath::Cos(2 * GetPhiInRange(track.phi())) * modqyZDCC;
        auto y2Tx1A = TMath::Sin(2 * GetPhiInRange(track.phi())) * modqxZDCA;
        auto y2Tx1C = TMath::Sin(2 * GetPhiInRange(track.phi())) * modqxZDCC;
        auto y2Ty1A = TMath::Sin(2 * GetPhiInRange(track.phi())) * modqyZDCA;
        auto y2Ty1C = TMath::Sin(2 * GetPhiInRange(track.phi())) * modqyZDCC;

        if (globalpt) {
          // if (sign > 0) {
          histos.fill(HIST("hpuxQxpvscentpteta"), centrality, track.pt(), track.eta(), uxQxp);
          histos.fill(HIST("hpuyQypvscentpteta"), centrality, track.pt(), track.eta(), uyQyp);
          histos.fill(HIST("hpuxQxtvscentpteta"), centrality, track.pt(), track.eta(), uxQxt);
          histos.fill(HIST("hpuyQytvscentpteta"), centrality, track.pt(), track.eta(), uyQyt);

          histos.fill(HIST("hpuxvscentpteta"), centrality, track.pt(), track.eta(), ux);
          histos.fill(HIST("hpuyvscentpteta"), centrality, track.pt(), track.eta(), uy);

          histos.fill(HIST("hpuxyQxytvscentpteta"), centrality, track.pt(), track.eta(), uxyQxyt);
          histos.fill(HIST("hpuxyQxypvscentpteta"), centrality, track.pt(), track.eta(), uxyQxyp);
          histos.fill(HIST("hpoddv1vscentpteta"), centrality, track.pt(), track.eta(), oddv1);
          histos.fill(HIST("hpevenv1vscentpteta"), centrality, track.pt(), track.eta(), evenv1);

          histos.fill(HIST("hpv21"), centrality, track.pt(), track.eta(), v21);
          histos.fill(HIST("hpv22"), centrality, track.pt(), track.eta(), v22);
          histos.fill(HIST("hpv23"), centrality, track.pt(), track.eta(), v23);

          histos.fill(HIST("hpx2Tx1Ax1Cvscentpteta"), centrality, track.pt(), track.eta(), x2Tx1Ax1C);
          histos.fill(HIST("hpx2Ty1Ay1Cvscentpteta"), centrality, track.pt(), track.eta(), x2Ty1Ay1C);
          histos.fill(HIST("hpy2Tx1Ay1Cvscentpteta"), centrality, track.pt(), track.eta(), y2Tx1Ay1C);
          histos.fill(HIST("hpy2Ty1Ax1Cvscentpteta"), centrality, track.pt(), track.eta(), y2Ty1Ax1C);
          histos.fill(HIST("hpx2Tvscentpteta"), centrality, track.pt(), track.eta(), x2T);
          histos.fill(HIST("hpy2Tvscentpteta"), centrality, track.pt(), track.eta(), y2T);
          histos.fill(HIST("hpx2Tx1Avscentpteta"), centrality, track.pt(), track.eta(), x2Tx1A);
          histos.fill(HIST("hpx2Tx1Cvscentpteta"), centrality, track.pt(), track.eta(), x2Tx1C);
          histos.fill(HIST("hpx2Ty1Avscentpteta"), centrality, track.pt(), track.eta(), x2Ty1A);
          histos.fill(HIST("hpx2Ty1Cvscentpteta"), centrality, track.pt(), track.eta(), x2Ty1C);
          histos.fill(HIST("hpy2Tx1Avscentpteta"), centrality, track.pt(), track.eta(), y2Tx1A);
          histos.fill(HIST("hpy2Ty1Cvscentpteta"), centrality, track.pt(), track.eta(), y2Ty1C);
          histos.fill(HIST("hpy2Ty1Avscentpteta"), centrality, track.pt(), track.eta(), y2Ty1A);
          histos.fill(HIST("hpy2Tx1Cvscentpteta"), centrality, track.pt(), track.eta(), y2Tx1C);
          histos.fill(HIST("hpx1Ax1Cvscentpteta"), centrality, track.pt(), track.eta(), x1Ax1C);
          histos.fill(HIST("hpy1Ay1Cvscentpteta"), centrality, track.pt(), track.eta(), y1Ay1C);
          histos.fill(HIST("hpx1Ay1Cvscentpteta"), centrality, track.pt(), track.eta(), x1Ay1C);
          histos.fill(HIST("hpy1Ax1Cvscentpteta"), centrality, track.pt(), track.eta(), x1Cy1A);
          histos.fill(HIST("hpx1Avscentpteta"), centrality, track.pt(), track.eta(), x1A);
          histos.fill(HIST("hpx1Cvscentpteta"), centrality, track.pt(), track.eta(), x1C);
          histos.fill(HIST("hpy1Avscentpteta"), centrality, track.pt(), track.eta(), y1A);
          histos.fill(HIST("hpy1Cvscentpteta"), centrality, track.pt(), track.eta(), y1C);

          /*} else {
            histos.fill(HIST("hpuxQxpvscentptetaneg"), centrality, track.pt(), track.eta(), uxQxp);
            histos.fill(HIST("hpuyQypvscentptetaneg"), centrality, track.pt(), track.eta(), uyQyp);
            histos.fill(HIST("hpuxQxtvscentptetaneg"), centrality, track.pt(), track.eta(), uxQxt);
            histos.fill(HIST("hpuyQytvscentptetaneg"), centrality, track.pt(), track.eta(), uyQyt);

            histos.fill(HIST("hpuxvscentptetaneg"), centrality, track.pt(), track.eta(), ux);
            histos.fill(HIST("hpuyvscentptetaneg"), centrality, track.pt(), track.eta(), uy);

            histos.fill(HIST("hpuxyQxytvscentptetaneg"), centrality, track.pt(), track.eta(), uxyQxyt);
            histos.fill(HIST("hpuxyQxypvscentptetaneg"), centrality, track.pt(), track.eta(), uxyQxyp);
            histos.fill(HIST("hpoddv1vscentptetaneg"), centrality, track.pt(), track.eta(), oddv1);
            histos.fill(HIST("hpevenv1vscentptetaneg"), centrality, track.pt(), track.eta(), evenv1);
            }*/
        } else {
          histos.fill(HIST("hpuxQxpvscentpteta"), centrality, track.tpcInnerParam(), track.eta(), uxQxp);
          histos.fill(HIST("hpuyQypvscentpteta"), centrality, track.tpcInnerParam(), track.eta(), uyQyp);
          histos.fill(HIST("hpuxQxtvscentpteta"), centrality, track.tpcInnerParam(), track.eta(), uxQxt);
          histos.fill(HIST("hpuyQytvscentpteta"), centrality, track.tpcInnerParam(), track.eta(), uyQyt);

          histos.fill(HIST("hpuxvscentpteta"), centrality, track.pt(), track.eta(), ux);
          histos.fill(HIST("hpuyvscentpteta"), centrality, track.pt(), track.eta(), uy);

          histos.fill(HIST("hpuxyQxytvscentpteta"), centrality, track.tpcInnerParam(), track.eta(), uxyQxyt);
          histos.fill(HIST("hpuxyQxypvscentpteta"), centrality, track.tpcInnerParam(), track.eta(), uxyQxyp);
          histos.fill(HIST("hpoddv1vscentpteta"), centrality, track.pt(), track.eta(), oddv1);
          histos.fill(HIST("hpevenv1vscentpteta"), centrality, track.pt(), track.eta(), evenv1);
        }
      }
    } else {
      for (const auto& v0 : V0s) {

        auto postrack = v0.template posTrack_as<AllTrackCandidates>();
        auto negtrack = v0.template negTrack_as<AllTrackCandidates>();

        int LambdaTag = 0;
        int aLambdaTag = 0;

        const auto signpos = postrack.sign();
        const auto signneg = negtrack.sign();

        if (signpos < 0 || signneg > 0) {
          continue;
        }

        if (isSelectedV0Daughter(v0, postrack, 0) && isSelectedV0Daughter(v0, negtrack, 1)) {
          LambdaTag = 1;
        }
        if (isSelectedV0Daughter(v0, negtrack, 0) && isSelectedV0Daughter(v0, postrack, 1)) {
          aLambdaTag = 1;
        }

        if (!LambdaTag && !aLambdaTag)
          continue;

        if (!SelectionV0(collision, v0)) {
          continue;
        }

        if (LambdaTag) {
          Proton = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPr);
          AntiPion = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPi);
          Lambdadummy = Proton + AntiPion;
          // angleLambda = calculateAngleBetweenLorentzVectors(Proton, AntiPion);
        }
        if (aLambdaTag) {
          AntiProton = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPr);
          Pion = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPi);
          AntiLambdadummy = AntiProton + Pion;
          // angleAntiLambda = calculateAngleBetweenLorentzVectors(AntiProton, Pion);
        }

        if (shouldReject(LambdaTag, aLambdaTag, Lambdadummy, AntiLambdadummy)) {
          continue;
        }

        if (TMath::Abs(v0.eta()) > 0.8)
          continue;

        int taga = LambdaTag;
        int tagb = aLambdaTag;

        if (useAccCorr && (currentRunNumber != lastRunNumber)) {
          accprofileL = ccdb->getForTimeStamp<TProfile2D>(ConfAccPathL.value, bc.timestamp());
          accprofileAL = ccdb->getForTimeStamp<TProfile2D>(ConfAccPathAL.value, bc.timestamp());
        }

        float desbinvalue = 0.0;
        if (dosystematic) {
          ////////////////////////////////////////////////////
          float LTsys = TMath::Abs(v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massLambda);
          float CPAsys = v0.v0cosPA();
          float DCADaughsys = TMath::Abs(v0.dcaV0daughters());
          float DCApossys = TMath::Abs(v0.dcapostopv());
          float DCAnegsys = TMath::Abs(v0.dcanegtopv());
          float sysvar = -999.9;
          double syst[10];
          if (sys == 1) {
            double temp[10] = {26, 27, 28, 29, 30, 31, 32, 33, 34, 35};
            std::copy(std::begin(temp), std::end(temp), std::begin(syst));
            sysvar = LTsys;
          }
          if (sys == 2) {
            double temp[10] = {0.992, 0.993, 0.9935, 0.994, 0.9945, 0.995, 0.9955, 0.996, 0.9965, 0.997};
            std::copy(std::begin(temp), std::end(temp), std::begin(syst));
            sysvar = CPAsys;
          }
          if (sys == 3) {
            double temp[10] = {0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25};
            std::copy(std::begin(temp), std::end(temp), std::begin(syst));
            sysvar = DCADaughsys;
          }
          if (sys == 4) {
            double temp[10] = {0.05, 0.07, 0.1, 0.15, 0.18, 0.2, 0.22, 0.25, 0.28, 0.3};
            std::copy(std::begin(temp), std::end(temp), std::begin(syst));
            sysvar = DCApossys;
          }
          if (sys == 5) {
            double temp[10] = {0.05, 0.07, 0.1, 0.15, 0.18, 0.2, 0.22, 0.25, 0.28, 0.3};
            std::copy(std::begin(temp), std::end(temp), std::begin(syst));
            sysvar = DCAnegsys;
          }

          for (int i = 0; i < 10; i++) {
            if (sys == 1 || sys == 3) {
              if (sysvar < syst[i])
                desbinvalue = i + 0.5;
              else
                continue;
            }
            if (sys == 2 || sys == 4 || sys == 5) {
              if (sysvar > syst[i])
                desbinvalue = i + 0.5;
              else
                continue;
            }

            ///////////////////////////////////////////////////
            if (LambdaTag) {
              Lambda = Proton + AntiPion;
              tagb = 0;
              double acvalue = accprofileL->GetBinContent(accprofileL->FindBin(v0.eta(), v0.pt()));
              fillHistograms(taga, tagb, Lambda, Proton, psiZDCC, psiZDCA, psiZDC, centrality, v0.mLambda(), v0.pt(), desbinvalue, acvalue);
            }

            tagb = aLambdaTag;
            if (aLambdaTag) {
              AntiLambda = AntiProton + Pion;
              taga = 0;
              double acvalue = accprofileAL->GetBinContent(accprofileAL->FindBin(v0.eta(), v0.pt()));
              fillHistograms(taga, tagb, AntiLambda, AntiProton, psiZDCC, psiZDCA, psiZDC, centrality, v0.mAntiLambda(), v0.pt(), desbinvalue, acvalue);
            }
          }
        } else {
          if (LambdaTag) {
            Lambda = Proton + AntiPion;
            tagb = 0;
            double acvalue = accprofileL->GetBinContent(accprofileL->FindBin(v0.eta(), v0.pt()));
            fillHistograms(taga, tagb, Lambda, Proton, psiZDCC, psiZDCA, psiZDC, centrality, v0.mLambda(), v0.pt(), v0.eta(), acvalue);
          }

          tagb = aLambdaTag;
          if (aLambdaTag) {
            AntiLambda = AntiProton + Pion;
            taga = 0;
            double acvalue = accprofileAL->GetBinContent(accprofileAL->FindBin(v0.eta(), v0.pt()));
            fillHistograms(taga, tagb, AntiLambda, AntiProton, psiZDCC, psiZDCA, psiZDC, centrality, v0.mAntiLambda(), v0.pt(), v0.eta(), acvalue);
          }
        }
      }
    }
    lastRunNumber = currentRunNumber;
  }
  PROCESS_SWITCH(lambdapolsp, processData, "Process data", true);

  // process function for derived data - mimics the functionality of the original data
  void processDerivedData(soa::Join<aod::StraCollisions, aod::StraCents, aod::StraEvSels, aod::StraStamps, aod::StraZDCSP>::iterator const& collision, v0Candidates const& V0s, dauTracks const&, BCsRun3 const&)
  {
    //___________________________________________________________________________________________________
    // event selection
    if (!collision.sel8()) {
      return;
    }
    auto centrality = collision.centFT0C();
    if (!collision.triggereventsp()) { // provided by StraZDCSP
      return;
    }

    if (additionalEvSel && (!collision.selection_bit(aod::evsel::kNoSameBunchPileup) || !collision.selection_bit(aod::evsel::kIsGoodZvtxFT0vsPV))) {
      return;
    }
    // histos.fill(HIST("hCentrality2"), centrality);
    //  if (additionalEvSel2 && (!collision.selection_bit(o2::aod::evsel::kNoCollInTimeRangeStandard))) {
    if (additionalEvSel2 && (collision.trackOccupancyInTimeRange() > cfgMaxOccupancy || collision.trackOccupancyInTimeRange() < cfgMinOccupancy)) {
      return;
    }
    // histos.fill(HIST("hCentrality3"), centrality);
    if (additionalEvSel3 && (!collision.selection_bit(aod::evsel::kNoTimeFrameBorder) || !collision.selection_bit(aod::evsel::kNoITSROFrameBorder))) {
      return;
    }

    if (additionalEvSel4 && !collision.selection_bit(o2::aod::evsel::kIsGoodITSLayersAll)) {
      return;
    }

    /*currentRunNumber = collision.foundBC_as<BCsRun3>().runNumber();
    auto bc = collision.foundBC_as<BCsRun3>();

    if (useAccCorr && (currentRunNumber != lastRunNumber)) {
      accprofileL = ccdb->getForTimeStamp<TProfile2D>(ConfAccPathL.value, bc.timestamp());
      accprofileAL = ccdb->getForTimeStamp<TProfile2D>(ConfAccPathAL.value, bc.timestamp());
    }
    */
    //___________________________________________________________________________________________________
    // retrieve further info provided by StraZDCSP
    auto qxZDCA = collision.qxZDCA();
    auto qxZDCC = collision.qxZDCC();
    auto qyZDCA = collision.qyZDCA();
    auto qyZDCC = collision.qyZDCC();
    auto psiZDCC = collision.psiZDCC();
    auto psiZDCA = collision.psiZDCA();
    double modqxZDCA;
    double modqyZDCA;
    double modqxZDCC;
    double modqyZDCC;

    if (cqvas) {
      modqxZDCA = TMath::Sqrt((qxZDCA * qxZDCA) + (qyZDCA * qyZDCA)) * TMath::Cos(psiZDCA);
      modqyZDCA = TMath::Sqrt((qxZDCA * qxZDCA) + (qyZDCA * qyZDCA)) * TMath::Sin(psiZDCA);
      modqxZDCC = TMath::Sqrt((qxZDCC * qxZDCC) + (qyZDCC * qyZDCC)) * TMath::Cos(psiZDCC);
      modqyZDCC = TMath::Sqrt((qxZDCC * qxZDCC) + (qyZDCC * qyZDCC)) * TMath::Sin(psiZDCC);
    } else {
      modqxZDCA = qxZDCA;
      modqyZDCA = qyZDCA;
      modqxZDCC = qxZDCC;
      modqyZDCC = qyZDCC;
    }

    auto psiZDC = TMath::ATan2((modqyZDCC - modqyZDCA), (modqxZDCC - modqxZDCA)); // full event plane

    // fill histograms
    histos.fill(HIST("hCentrality"), centrality);
    if (!checkwithpub) {
      // histos.fill(HIST("hVtxZ"), collision.posZ());
      histos.fill(HIST("hpRes"), centrality, (TMath::Cos(GetPhiInRange(psiZDCA - psiZDCC))));
      histos.fill(HIST("hpResSin"), centrality, (TMath::Sin(GetPhiInRange(psiZDCA - psiZDCC))));
      /*histos.fill(HIST("hpCosPsiA"), centrality, (TMath::Cos(GetPhiInRange(psiZDCA))));
      histos.fill(HIST("hpCosPsiC"), centrality, (TMath::Cos(GetPhiInRange(psiZDCC))));
      histos.fill(HIST("hpSinPsiA"), centrality, (TMath::Sin(GetPhiInRange(psiZDCA))));
      histos.fill(HIST("hpSinPsiC"), centrality, (TMath::Sin(GetPhiInRange(psiZDCC))));*/
    }

    //___________________________________________________________________________________________________
    // loop over V0s as necessary
    for (const auto& v0 : V0s) {
      bool LambdaTag = isCompatible(v0, 0);
      bool aLambdaTag = isCompatible(v0, 1);

      if (!LambdaTag && !aLambdaTag)
        continue;

      if (!SelectionV0(collision, v0)) {
        continue;
      }

      if (LambdaTag) {
        Proton = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPr);
        AntiPion = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPi);
        Lambdadummy = Proton + AntiPion;
        // angleLambda = calculateAngleBetweenLorentzVectors(Proton, AntiPion);
      }
      if (aLambdaTag) {
        AntiProton = ROOT::Math::PxPyPzMVector(v0.pxneg(), v0.pyneg(), v0.pzneg(), massPr);
        Pion = ROOT::Math::PxPyPzMVector(v0.pxpos(), v0.pypos(), v0.pzpos(), massPi);
        AntiLambdadummy = AntiProton + Pion;
        // angleAntiLambda = calculateAngleBetweenLorentzVectors(AntiProton, Pion);
      }

      if (shouldReject(LambdaTag, aLambdaTag, Lambdadummy, AntiLambdadummy)) {
        continue;
      }

      if (TMath::Abs(v0.eta()) > 0.8)
        continue;

      int taga = LambdaTag;
      int tagb = aLambdaTag;

      float desbinvalue = 0.0;
      if (dosystematic) {
        ////////////////////////////////////////////////////
        float LTsys = TMath::Abs(v0.distovertotmom(collision.posX(), collision.posY(), collision.posZ()) * massLambda);
        float CPAsys = v0.v0cosPA();
        float DCADaughsys = TMath::Abs(v0.dcaV0daughters());
        float DCApossys = TMath::Abs(v0.dcapostopv());
        float DCAnegsys = TMath::Abs(v0.dcanegtopv());
        float sysvar = -999.9;
        double syst[10];
        if (sys == 1) {
          double temp[10] = {26, 27, 28, 29, 30, 31, 32, 33, 34, 35};
          std::copy(std::begin(temp), std::end(temp), std::begin(syst));
          sysvar = LTsys;
        }
        if (sys == 2) {
          double temp[10] = {0.992, 0.993, 0.9935, 0.994, 0.9945, 0.995, 0.9955, 0.996, 0.9965, 0.997};
          std::copy(std::begin(temp), std::end(temp), std::begin(syst));
          sysvar = CPAsys;
        }
        if (sys == 3) {
          double temp[10] = {0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25};
          std::copy(std::begin(temp), std::end(temp), std::begin(syst));
          sysvar = DCADaughsys;
        }
        if (sys == 4) {
          double temp[10] = {0.05, 0.07, 0.1, 0.15, 0.18, 0.2, 0.22, 0.25, 0.28, 0.3};
          std::copy(std::begin(temp), std::end(temp), std::begin(syst));
          sysvar = DCApossys;
        }
        if (sys == 5) {
          double temp[10] = {0.05, 0.07, 0.1, 0.15, 0.18, 0.2, 0.22, 0.25, 0.28, 0.3};
          std::copy(std::begin(temp), std::end(temp), std::begin(syst));
          sysvar = DCAnegsys;
        }

        for (int i = 0; i < 10; i++) {
          if (sys == 1 || sys == 3) {
            if (sysvar < syst[i])
              desbinvalue = i + 0.5;
            else
              continue;
          }
          if (sys == 2 || sys == 4 || sys == 5) {
            if (sysvar > syst[i])
              desbinvalue = i + 0.5;
            else
              continue;
          }

          ///////////////////////////////////////////////////
          if (LambdaTag) {
            Lambda = Proton + AntiPion;
            tagb = 0;
            double acvalue = 1.0;
            fillHistograms(taga, tagb, Lambda, Proton, psiZDCC, psiZDCA, psiZDC, centrality, v0.mLambda(), v0.pt(), desbinvalue, acvalue);
          }

          tagb = aLambdaTag;
          if (aLambdaTag) {
            AntiLambda = AntiProton + Pion;
            taga = 0;
            double acvalue = 1.0;
            fillHistograms(taga, tagb, AntiLambda, AntiProton, psiZDCC, psiZDCA, psiZDC, centrality, v0.mAntiLambda(), v0.pt(), desbinvalue, acvalue);
          }
        }
      } else {
        if (LambdaTag) {
          Lambda = Proton + AntiPion;
          tagb = 0;
          double acvalue = 1.0;
          fillHistograms(taga, tagb, Lambda, Proton, psiZDCC, psiZDCA, psiZDC, centrality, v0.mLambda(), v0.pt(), v0.eta(), acvalue);
        }

        tagb = aLambdaTag;
        if (aLambdaTag) {
          AntiLambda = AntiProton + Pion;
          taga = 0;
          double acvalue = 1.0;
          fillHistograms(taga, tagb, AntiLambda, AntiProton, psiZDCC, psiZDCA, psiZDC, centrality, v0.mAntiLambda(), v0.pt(), v0.eta(), acvalue);
        }
      }
    }
    // lastRunNumber = currentRunNumber;
  }
  PROCESS_SWITCH(lambdapolsp, processDerivedData, "Process derived data", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<lambdapolsp>(cfgc, TaskName{"lambdapolsp"})};
}
