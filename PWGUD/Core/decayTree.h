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

#ifndef PWGUD_CORE_DECAYTREE_H_
#define PWGUD_CORE_DECAYTREE_H_

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/Logger.h"
#include "Framework/O2DatabasePDGPlugin.h"

#include "TDatabasePDG.h"
#include "TLorentzVector.h"

#include <map>
#include <string>
#include <utility>
#include <vector>

// -----------------------------------------------------------------------------
class pidSelector
{
 public:
  // constructor/destructor
  pidSelector() {}
  explicit pidSelector(std::vector<std::vector<double>>& pidcuts);
  ~pidSelector() {}

  // setter
  void clear();

  // getters
  void Print();

  // templated functions
  template <typename TTs>
  double getTPCnSigma(TTs track, int pid)
  {
    auto hypo = pid2ind(pid);
    switch (hypo) {
      case 0:
        return track.tpcNSigmaEl();
      case 1:
        return track.tpcNSigmaPi();
      case 2:
        return track.tpcNSigmaMu();
      case 3:
        return track.tpcNSigmaKa();
      case 4:
        return track.tpcNSigmaPr();
      default:
        return 0.;
    }
  };

  template <typename TTs>
  double getTOFnSigma(TTs track, int pid)
  {
    auto hypo = pid2ind(pid);
    switch (hypo) {
      case 0:
        return track.tofNSigmaEl();
      case 1:
        return track.tofNSigmaPi();
      case 2:
        return track.tofNSigmaMu();
      case 3:
        return track.tofNSigmaKa();
      case 4:
        return track.tofNSigmaPr();
      default:
        return 0.;
    }
  };

  template <typename TTs>
  bool goodTrack(TTs track)
  {
    // loop over pidcuts
    for (const auto& pidcut : fpidCuts) {

      float mom = 0.;
      float detValue = 0.;
      if (pidcut[1] == 1) {
        // TPC
        if (track.hasTPC()) {
          mom = track.tpcInnerParam();
          if (mom < pidcut[4] || mom > pidcut[5]) {
            // not in relevant momentum range
            continue;
          }
          if (std::abs(pidcut[2]) == 1) {
            // nSigma
            detValue = getTPCnSigma(track, pidcut[0]);
          } else {
            // signal
            detValue = track.tpcSignal();
          }
        } else {
          if (pidcut[3] == 2) {
            // TPC is required
            return false;
          } else {
            continue;
          }
        }
      } else {
        // TOF
        if (track.hasTOF()) {
          mom = track.tofExpMom();
          if (mom < pidcut[4] || mom > pidcut[5]) {
            // not in relevant momentum range
            continue;
          }
          if (std::abs(pidcut[2]) == 1) {
            // nSigma
            detValue = getTOFnSigma(track, pidcut[0]);
          } else {
            // signal
            detValue = track.tofSignal();
          }
        } else {
          if (pidcut[3] == 2) {
            // TOF is required
            return false;
          } else {
            continue;
          }
        }
      }

      // inclusive / exclusive
      if (pidcut[2] > 0 && (detValue < pidcut[6] || detValue > pidcut[7])) {
        return false;
      } else if (pidcut[2] < 0 && (detValue > pidcut[6] && detValue < pidcut[7])) {
        return false;
      }
    }

    // the track is good if we arrive here
    return true;
  }

 private:
  std::vector<std::vector<double>> fpidCuts;
  int pid2ind(int pid);

  // ClassDefNV(pidSelector, 1);
};

// -----------------------------------------------------------------------------
class angleCut
{
 public:
  // constructor/destructor
  angleCut() {}
  explicit angleCut(std::pair<std::string, std::string> rNames, double angleMin, double angleMax);
  ~angleCut() {}

  std::pair<std::string, std::string> rNames() { return fRnames; }
  std::pair<double, double> angleRange() { return std::pair<double, double>{fAngleMin, fAngleMax}; }

  void Print();

 private:
  std::pair<std::string, std::string> fRnames;
  double fAngleMin;
  double fAngleMax;

  // ClassDefNV(angleCut, 1);
};

// -----------------------------------------------------------------------------
class reconstructedParticle
{
 public:
  // constructor/destructor
  reconstructedParticle() {}
  explicit reconstructedParticle(std::string name, TLorentzVector ivm, std::vector<int>& comb)
  {
    fName = name;
    fIVM = ivm;
    fComb = comb;
  };
  ~reconstructedParticle() {}

  std::string name() { return fName; }
  TLorentzVector lv() { return fIVM; }
  std::vector<int> comb() { return fComb; }

 private:
  std::string fName;
  TLorentzVector fIVM;
  std::vector<int> fComb;
};

using recResType = std::vector<std::map<std::string, reconstructedParticle>>;

class reconstructedEvent
{
 public:
  // constructor/destructor
  reconstructedEvent() {}
  explicit reconstructedEvent(recResType recs, int chargeState, std::vector<int>& comb)
  {
    fRecs = recs;
    fComb = comb;
    fChargeState = chargeState;
  };
  ~reconstructedEvent() {}

  recResType recResonances() { return fRecs; }
  std::vector<int> comb() { return fComb; }

 private:
  recResType fRecs;
  std::vector<int> fComb;
  int fChargeState;
};

using decayTreeResType = std::map<std::string, recResType>;

// -----------------------------------------------------------------------------
class resonance
{
 public:
  // constructor/destructor
  resonance();
  ~resonance() {}

  // setters
  void init();
  void reset();

  void setisFinal() { fisFinal = true; }
  void setCounter(int counter) { fCounter = counter; }
  void setName(std::string name) { fName = name; }
  void setStatus(int status) { fStatus = status; }
  void setPID(int pid) { fPID = pid; }
  void setPIDFun(int pidfun) { fPIDfun = pidfun; }
  void setDetectorHits(int its, int tpc, int trd, int tof)
  {
    fdetectorHits = std::vector<int>{its, tpc, trd, tof};
  }
  void clearParents() { fParents.clear(); }
  void addParent(std::string parent) { fParents.push_back(parent); }
  void setDaughters(std::vector<std::string>& daughters) { fDaughters = daughters; }
  void setIVM(TLorentzVector ivm)
  {
    fIVM = ivm;
    fStatus = 1;
  }
  void setCharge(int charge) { fCharge = charge; }

  // selections
  void setMassRange(double mmin, double mmax)
  {
    fmassMin = mmin;
    fmassMax = mmax;
  }
  void setPtRange(double ptmin, double ptmax)
  {
    fptMin = ptmin;
    fptMax = ptmax;
  }
  void setEtaRange(double etamin, double etamax)
  {
    fetaMin = etamin;
    fetaMax = etamax;
  }
  void setNcltpcRange(int ncltpcmin, int ncltpcmax)
  {
    fncltpcMin = ncltpcmin;
    fncltpcMax = ncltpcmax;
  }
  void setChi2ncltpcRange(double chi2ncltpcmin, double chi2ncltpcmax)
  {
    fchi2ncltpcMin = chi2ncltpcmin;
    fchi2ncltpcMax = chi2ncltpcmax;
  }
  void setDCAxyzMax(double dcaxymax, double dcazmax)
  {
    fdcaxyMax = dcaxymax;
    fdcazMax = dcazmax;
  }
  void setPIDSelector(pidSelector pidcuts) { fpidSelector = pidcuts; }
  void setAngleCuts(std::vector<angleCut*> anglecuts) { fangleCuts = anglecuts; }

  // histograms
  void setMassHistAxis(int nbins, double binmin, double binmax)
  {
    fnmassBins = nbins;
    fmassHistMin = binmin;
    fmassHistMax = binmax;
  }
  void setMomHistAxis(int nbins, double binmin, double binmax)
  {
    fnmomBins = nbins;
    fmomHistMin = binmin;
    fmomHistMax = binmax;
  }
  void updateStatus();

  // getters
  bool isFinal() { return fisFinal; }
  int counter() { return fCounter; }
  std::string name() { return fName; }
  int status() { return fStatus; }
  int pid() { return fPID; }
  int pidFun() { return fPIDfun; }
  std::vector<int> detectorHits() { return fdetectorHits; }
  std::vector<std::string> getParents() { return fParents; }
  std::vector<std::string> getDaughters() { return fDaughters; }
  double massMin() { return fmassMin; }
  double massMax() { return fmassMax; }
  double ptMin() { return fptMin; }
  double ptMax() { return fptMax; }
  double etaMin() { return fetaMin; }
  double etaMax() { return fetaMax; }
  TLorentzVector IVM() { return fIVM; }
  int charge() { return fCharge; }
  pidSelector getPIDSelector() { return fpidSelector; }
  std::vector<angleCut*> getAngleCuts() { return fangleCuts; }

  // histograms
  int nmassBins() { return fnmassBins; }
  std::vector<double> massHistRange() { return std::vector<double>({fmassHistMin, fmassHistMax}); }
  int nmomBins() { return fnmomBins; }
  std::vector<double> momHistRange() { return std::vector<double>({fmomHistMin, fmomHistMax}); }

  void Print();

  // templated functions
  // resonance status
  //  0: unset
  //  1: IVM calculated
  //  2: not accepted
  //  3: accepted
  template <typename TTs>
  void updateStatus(TTs const& track)
  {
    // IVM has to be computed
    if (fStatus == 0) {
      return;
    }

    // check mass, pt, eta range and charge
    updateStatus();

    // check detector hits, track cuts
    if (fStatus >= 3) {
      // detector hits
      if (fdetectorHits[0] >= 0) {
        if ((fdetectorHits[0] == 0 && track.hasITS()) || (fdetectorHits[0] > 0 && !track.hasITS())) {
          fStatus = 2;
        }
      }
      if (fdetectorHits[1] >= 0) {
        if ((fdetectorHits[1] == 0 && track.hasTPC()) || (fdetectorHits[1] > 0 && !track.hasTPC())) {
          fStatus = 2;
        }
      }
      if (fdetectorHits[2] >= 0) {
        if ((fdetectorHits[2] == 0 && track.hasTRD()) || (fdetectorHits[2] > 0 && !track.hasTRD())) {
          fStatus = 2;
        }
      }
      if (fdetectorHits[3] >= 0) {
        if ((fdetectorHits[3] == 0 && track.hasTOF()) || (fdetectorHits[3] > 0 && !track.hasTOF())) {
          fStatus = 2;
        }
      }
      // PID cuts
      if (!fpidSelector.goodTrack(track)) {
        fStatus = 2;
      }
      // nclTPC
      auto nclTPC = track.tpcNClsFindable() - track.tpcNClsFindableMinusFound();
      if (nclTPC < fncltpcMin || nclTPC > fncltpcMax) {
        fStatus = 2;
      }
      // chi2nclTPC
      if (track.tpcChi2NCl() < fchi2ncltpcMin || track.tpcChi2NCl() > fchi2ncltpcMax) {
        fStatus = 2;
      }
      // dcaxyz
      auto lim = fdcaxyMax + std::pow(0.0350 / track.pt(), 1.1);
      if (std::abs(track.dcaXY()) > lim || std::abs(track.dcaZ()) > fdcazMax) {
        fStatus = 2;
      }
    }
  }

 private:
  bool fisFinal;
  int fCounter;

  // resonance name
  std::string fName;
  int fStatus;

  // nominal pid
  int fPID;
  int fPIDfun;
  std::vector<int> fdetectorHits;

  // name of parents and daughters
  std::vector<std::string> fParents;
  std::vector<std::string> fDaughters;
  void updateParents();

  // mass, pT, , eta range
  double fmassMin;
  double fmassMax;
  double fptMin;
  double fptMax;
  double fetaMin;
  double fetaMax;
  int fncltpcMin;
  int fncltpcMax;
  double fchi2ncltpcMin;
  double fchi2ncltpcMax;
  double fdcaxyMax;
  double fdcazMax;

  // histogram axes
  int fnmassBins;
  double fmassHistMax;
  double fmassHistMin;
  int fnmomBins;
  double fmomHistMax;
  double fmomHistMin;

  // invariant mass
  TLorentzVector fIVM;
  int fCharge;

  // pidcuts, anglecuts
  pidSelector fpidSelector;
  std::vector<angleCut*> fangleCuts;

  // ClassDefNV(resonance, 1);
};

// -----------------------------------------------------------------------------
class decayTree
{
 public:
  // constructor/destructor
  decayTree();
  ~decayTree() {}

  // setters
  // read decay tree from json file
  bool init(std::string const& filename, o2::framework::HistogramRegistry& registry);

  // reset status of all resonances to 0
  void reset();
  void updateStatus();

  // getters
  int nFinals() { return fnFinals; }
  std::vector<resonance*> getResonances() { return fResonances; }
  resonance* getResonance(std::string name);
  resonance* getFinal(int counter);
  std::vector<resonance*> getFinals(resonance* res);
  std::vector<int> ntrackRange() { return std::vector<int>{fnTracksMin, fnTracksMax}; }
  double rgtrTOFMin() { return frgtwtofMin; }
  std::vector<int> dBCRange() { return std::vector<int>{fdBCMin, fdBCMax}; }
  std::vector<int> FITvetos() { return fFITvetos; }

  void Print();

  template <typename TTs>
  decayTreeResType processTree(TTs const& tracks, bool withFill = true)
  {
    recResType ULSresults;
    recResType LSresults;

    // return if nFinals > tracks.size()
    if (fnFinals > tracks.size()) {
      LOGF(info, "Number of tracks (%d) is smaller than the number of finals (%d)", tracks.size(), fnFinals);
      return decayTreeResType{{"ULS", ULSresults}, {"LS", LSresults}};
    }

    // create all possible track combinations including permutations
    auto combs = combinations(tracks.size());

    // a vector to keep track of successful combinations
    std::vector<std::size_t> goodCombs;

    // loop over possible combinations
    LOGF(debug, "New event");
    for (auto& comb : combs) {
      std::string scomb("");
      for (const auto& i : comb) {
        scomb.append(" ").append(std::to_string(i));
      }
      LOGF(debug, "  combination:%s", scomb);

      // has an equivalent combination been accepted already?
      auto newHash = combHash(comb);
      if (std::find(goodCombs.begin(), goodCombs.end(), newHash) != goodCombs.end()) {
        LOGF(debug, "    Equivalent combination is already accepted!");
        continue;
      }

      // loop over resonances and compute
      reset();
      for (auto res : fResonances) {
        computeResonance(res, tracks, comb);
      }

      // check angles between daughters of all resonances
      checkAngles();

      // check status of all resonances
      updateStatus();
      if (fStatus >= 2) {
        goodCombs.push_back(newHash);
        std::map<std::string, reconstructedParticle> recResonances;
        for (const auto& res : fResonances) {
          recResonances.insert({res->name(), reconstructedParticle(res->name(), res->IVM(), comb)});
        }

        if (fStatus == 2) {
          ULSresults.push_back(recResonances);
        } else {
          LSresults.push_back(recResonances);
        }
      }
    }
    auto results = decayTreeResType{{"ULS", ULSresults}, {"LS", LSresults}};
    if (withFill) {
      fillHistograms(results, tracks);
    }
    return results;
  }

#define getHist(type, name) std::get<std::shared_ptr<type>>(fhistPointers[name])
  template <typename TTs>
  void fillHistograms(decayTreeResType results, TTs const& tracks)
  {
    // fill the histograms
    std::string base;
    std::string hname;

    // results["ULS"] contains the ULS results
    // results["LS"] contains the LS results
    for (const auto& cc : fccs) {
      // result is a std::vector<std::map<std::string, reconstructedParticle>>
      for (auto result : results[cc]) {

        // loop over the reconstructed particles
        //  rec.first:  name of the reconstructed particle
        //  rec.second: reconstructed particle
        for (auto rec : result) {
          auto lv = rec.second.lv();
          base = cc;
          base.append("/").append(rec.first).append("/");

          hname = base + "mpt";
          getHist(TH2, hname)->Fill(lv.M(), lv.Perp(), 1.);
          hname = base + "meta";
          getHist(TH2, hname)->Fill(lv.M(), lv.Eta(), 1.);
          hname = base + "pteta";
          getHist(TH2, hname)->Fill(lv.Perp(), lv.Eta(), 1.);

          // M vs daughters
          auto res = getResonance(rec.first);
          auto daughs = res->getDaughters();
          auto ndaughs = daughs.size();
          for (auto i = 0; i < static_cast<int>(ndaughs); i++) {
            auto d1 = getResonance(daughs[i]);

            // M vs pT daughter
            hname = base;
            hname.append("MvspT_").append(rec.first).append(d1->name());
            getHist(TH2, hname)->Fill(lv.M(), result[d1->name()].lv().Perp(), 1.);

            // M vs eta daughter
            hname = base;
            hname.append("Mvseta_").append(rec.first).append(d1->name());
            getHist(TH2, hname)->Fill(lv.M(), result[d1->name()].lv().Eta(), 1.);

            if (d1->isFinal()) {
              auto tr = tracks.begin() + result[d1->name()].comb()[d1->counter()];

              // M vs dca
              hname = base;
              hname.append("MvsdcaXY_").append(rec.first).append(d1->name());
              getHist(TH2, hname)->Fill(lv.M(), tr.dcaXY(), 1.);
              hname = base;
              hname.append("MvsdcaZ_").append(rec.first).append(d1->name());
              getHist(TH2, hname)->Fill(lv.M(), tr.dcaZ(), 1.);

              // M vs chi2 track
              hname = base;
              hname.append("Mvschi2_").append(rec.first).append(d1->name());
              getHist(TH2, hname)->Fill(lv.M(), tr.tpcChi2NCl(), 1.);

              // M vs nCl track
              hname = base;
              hname.append("MvsnCl_").append(rec.first).append(d1->name());
              getHist(TH2, hname)->Fill(lv.M(), tr.tpcNClsFindable() - tr.tpcNClsFindableMinusFound(), 1.);

              // M versus detector hits
              hname = base;
              hname.append("MvsdetHits_").append(rec.first).append(d1->name());
              auto ind = tr.hasITS() + tr.hasTPC() * 2 + tr.hasTRD() * 4 + tr.hasTOF() * 8;
              getHist(TH2, hname)->Fill(lv.M(), ind, 1.);
            } else {
              // M vs Mi
              hname = base;
              hname.append("MvsM_").append(rec.first).append(d1->name());
              getHist(TH2, hname)->Fill(lv.M(), result[d1->name()].lv().M(), 1.);
            }
          }

          // daughters vs daughters
          for (auto i = 0; i < static_cast<int>(ndaughs - 1); i++) {
            auto d1 = getResonance(daughs[i]);
            auto ivm1 = result[daughs[i]].lv();
            for (auto j = i + 1; j < static_cast<int>(ndaughs); j++) {
              auto d2 = getResonance(daughs[j]);
              auto ivm2 = result[daughs[j]].lv();

              // M1 vs M2
              hname = base;
              hname.append("MvsM_").append(d1->name()).append(d2->name());
              getHist(TH2, hname)->Fill(ivm1.M(), ivm2.M(), 1.);

              // angle(d1, d2)
              auto ang = ivm1.Angle(ivm2.Vect());
              hname = base;
              hname.append("angle_").append(d1->name()).append(d2->name());
              getHist(TH1, hname)->Fill(ang, 1.);

              // M vs angle(d1, d2)
              hname = base;
              hname.append("Mvsangle_").append(d1->name()).append(d2->name());
              getHist(TH2, hname)->Fill(lv.M(), ang, 1.);

              // both daughters are finals
              if (d1->isFinal() && d2->isFinal()) {
                auto tr1 = tracks.begin() + result[d1->name()].comb()[d1->counter()];
                auto tr2 = tracks.begin() + result[d2->name()].comb()[d2->counter()];

                // TPC signal vs TPC signal
                hname = base;
                hname.append("TPCsignal_").append(d1->name()).append(d2->name());
                getHist(TH2, hname)->Fill(tr1.tpcSignal(), tr2.tpcSignal(), 1.);
              }
            }
          }

          // finals specific histograms
          if (res->isFinal()) {
            auto tr = tracks.begin() + rec.second.comb()[res->counter()];

            // dca XYZ
            hname = base;
            hname.append("dcaXY");
            getHist(TH1, hname)->Fill(tr.dcaXY(), 1.);
            hname = base;
            hname.append("dcaZ");
            getHist(TH1, hname)->Fill(tr.dcaZ(), 1.);

            // TPC
            hname = base;
            hname.append("nS").append(fparts[0]).append(fdets[0]);
            getHist(TH2, hname)->Fill(tr.tpcInnerParam(), tr.tpcNSigmaEl(), 1.);
            hname = base;
            hname.append("nS").append(fparts[1]).append(fdets[0]);
            getHist(TH2, hname)->Fill(tr.tpcInnerParam(), tr.tpcNSigmaPi(), 1.);
            hname = base;
            hname.append("nS").append(fparts[2]).append(fdets[0]);
            getHist(TH2, hname)->Fill(tr.tpcInnerParam(), tr.tpcNSigmaMu(), 1.);
            hname = base;
            hname.append("nS").append(fparts[3]).append(fdets[0]);
            getHist(TH2, hname)->Fill(tr.tpcInnerParam(), tr.tpcNSigmaKa(), 1.);
            hname = base;
            hname.append("nS").append(fparts[4]).append(fdets[0]);
            getHist(TH2, hname)->Fill(tr.tpcInnerParam(), tr.tpcNSigmaPr(), 1.);

            // TOF
            if (tr.hasTOF()) {
              hname = base;
              hname.append("nS").append(fparts[0]).append(fdets[1]);
              getHist(TH2, hname)->Fill(tr.tofExpMom(), tr.tofNSigmaEl(), 1.);
              hname = base;
              hname.append("nS").append(fparts[1]).append(fdets[1]);
              getHist(TH2, hname)->Fill(tr.tofExpMom(), tr.tofNSigmaPi(), 1.);
              hname = base;
              hname.append("nS").append(fparts[2]).append(fdets[1]);
              getHist(TH2, hname)->Fill(tr.tofExpMom(), tr.tofNSigmaMu(), 1.);
              hname = base;
              hname.append("nS").append(fparts[3]).append(fdets[1]);
              getHist(TH2, hname)->Fill(tr.tofExpMom(), tr.tofNSigmaKa(), 1.);
              hname = base;
              hname.append("nS").append(fparts[4]).append(fdets[1]);
              getHist(TH2, hname)->Fill(tr.tofExpMom(), tr.tofNSigmaPr(), 1.);
            }

            // detector hits
            hname = base;
            hname.append("detectorHits");
            if (tr.hasITS()) {
              getHist(TH1, hname)->Fill(1, 1.);
            }
            if (tr.hasTPC()) {
              getHist(TH1, hname)->Fill(2, 1.);
            }
            if (tr.hasTRD()) {
              getHist(TH1, hname)->Fill(3, 1.);
            }
            if (tr.hasTOF()) {
              getHist(TH1, hname)->Fill(4, 1.);
            }
          }
        }
      }
    }
  }

 private:
  // decayTree status
  //  0: unset
  //  1: not accepted
  //  2: ULS accepted
  //  3: LS accepted
  int fStatus;
  TDatabasePDG* fPDG;

  // event requierements
  int fnTracksMin;
  int fnTracksMax;
  double frgtwtofMin;
  int fdBCMin;
  int fdBCMax;
  std::vector<int> fFITvetos;
  std::vector<int> fULSstates;
  std::vector<int> fLSstates;

  // vectors of Resonances
  std::vector<resonance*> fResonances;
  int fChargeState;

  // number of finals
  int fnFinals;
  std::vector<std::vector<int>> fPermutations;

  // histogram registry
  std::vector<std::string> fccs;
  std::vector<std::string> fdets;
  std::vector<std::string> fparts;
  std::map<std::string, o2::framework::HistPtr> fhistPointers;

  // generate parent information for all resonances
  void updateParents();

  // helper functions to compute combinations and permutations
  //  combination:  selection of n out of N
  //  permutation:  order of n selected items
  // create all permutations of all combinations
  std::size_t combHash(std::vector<int>& comb);
  void permutations(std::vector<int>& ref, int n0, int np, std::vector<std::vector<int>>& perms);
  int permutations(int n0, std::vector<std::vector<int>>& perms);
  void combinations(int n0, std::vector<int>& pool, int np, std::vector<int>& inds, int n,
                    std::vector<std::vector<int>>& combs);
  int combinations(int n0, int np, std::vector<std::vector<int>>& combs);
  std::vector<std::vector<int>> combinations(int nPool);

  // check all angle requirements
  void checkAngles();

  // compute the charge state
  int chargeState(std::vector<int> chs);
  void updateChargeState();

  // templated functions
  template <typename TTs>
  void computeResonance(resonance* res, TTs const& tracks, std::vector<int>& comb)
  {
    // if status > 0 then return
    if (res->status() > 0) {
      return;
    }

    // initialisations
    TLorentzVector ivm{0., 0., 0., 0.};
    int charge = 0;

    // is this a final state or a resonance
    if (res->isFinal()) {
      // is a final
      auto pdgparticle = fPDG->GetParticle(res->pid());
      auto track = (tracks.begin() + comb[res->counter()]);
      ivm.SetXYZM(track.px(), track.py(), track.pz(), pdgparticle->Mass());
      res->setIVM(ivm);
      res->setCharge(track.sign());
      res->setStatus(1);

      // apply cuts
      res->updateStatus(track);

    } else {
      // is a resonance
      // loop over daughters
      for (const auto& daughName : res->getDaughters()) {
        auto daugh = getResonance(daughName);
        computeResonance(daugh, tracks, comb);
        ivm += daugh->IVM();
        charge += daugh->charge();
      }
      res->setIVM(ivm);
      res->setCharge(charge);
      res->setStatus(1);

      // apply cuts
      res->updateStatus();
    }
  }

  // create histograms
  void createHistograms(o2::framework::HistogramRegistry& registry)
  {
    // definitions
    auto etax = o2::framework::AxisSpec(100, -1.5, 1.5);
    auto nSax = o2::framework::AxisSpec(300, -15.0, 15.0);
    auto chi2ax = o2::framework::AxisSpec(100, 0.0, 5.0);
    auto nClax = o2::framework::AxisSpec(170, 0.0, 170.0);
    auto angax = o2::framework::AxisSpec(315, 0.0, 3.15);
    auto dcaxyax = o2::framework::AxisSpec(400, -0.2, 0.2);
    auto dcazax = o2::framework::AxisSpec(600, -0.3, 0.3);
    auto sTPCax = o2::framework::AxisSpec(1000, 0., 1000.);

    std::string base;
    std::string hname;
    std::string annot;
    fhistPointers.clear();
    for (const auto& res : getResonances()) {
      auto max = o2::framework::AxisSpec(res->nmassBins(), res->massHistRange()[0], res->massHistRange()[1]);
      auto momax = o2::framework::AxisSpec(res->nmomBins(), res->momHistRange()[0], res->momHistRange()[1]);

      // M-pT, M-eta, pT-eta
      for (const auto& cc : fccs) {
        base = cc;
        base.append("/").append(res->name()).append("/");
        hname = base + "mpt";
        annot = "M versus pT; M (" + res->name() + ") GeV/c^{2}; pT (" + res->name() + ") GeV/c";
        fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH2F, {max, momax}})});
        hname = base + "meta";
        annot = "M versus eta; M (" + res->name() + ") GeV/c^{2}; eta (" + res->name() + ")";
        fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH2F, {max, etax}})});
        hname = base + "pteta";
        annot = "pT versus eta; pT (" + res->name() + ") GeV/c; eta (" + res->name() + ")";
        fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH2F, {momax, etax}})});

        // M versus daughters
        auto daughs = res->getDaughters();
        auto ndaughs = daughs.size();
        for (auto i = 0; i < static_cast<int>(ndaughs); i++) {
          auto d1 = getResonance(daughs[i]);

          // M vs pT daughter
          hname = base;
          hname.append("MvspT_").append(res->name()).append(d1->name());
          annot = "M versus pT; M (" + res->name() + ") GeV/c^{2}; pT (" + d1->name() + ") GeV/c";
          auto momax1 = o2::framework::AxisSpec(d1->nmomBins(), d1->momHistRange()[0], d1->momHistRange()[1]);
          fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH2F, {max, momax1}})});

          // M vs eta daughter
          hname = base;
          hname.append("Mvseta_").append(res->name()).append(d1->name());
          annot = "M versus eta; M (" + res->name() + ") GeV/c^{2}; eta (" + d1->name() + ")";
          fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH2F, {max, etax}})});

          if (d1->isFinal()) {
            // M vs dcaXYZ
            hname = base;
            hname.append("MvsdcaXY_").append(res->name()).append(d1->name());
            annot = "M versus dcaXY; M (" + res->name() + ") GeV/c^{2}; dca_{XY} (" + d1->name() + ") #mu m";
            fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH2F, {max, dcaxyax}})});
            hname = base;
            hname.append("MvsdcaZ_").append(res->name()).append(d1->name());
            annot = "M versus dcaZ; M (" + res->name() + ") GeV/c^{2}; dca_{Z} (" + d1->name() + ") #mu m";
            fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH2F, {max, dcazax}})});

            // M vs chi2 track
            hname = base;
            hname.append("Mvschi2_").append(res->name()).append(d1->name());
            annot = "M versus chi2; M (" + res->name() + ") GeV/c^{2}; chi2 (" + d1->name() + ")";
            fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH2F, {max, chi2ax}})});

            // M vs nCl track
            hname = base;
            hname.append("MvsnCl_").append(res->name()).append(d1->name());
            annot = "M versus nCl; M (" + res->name() + ") GeV/c^{2}; nCl (" + d1->name() + ")";
            fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH2F, {max, nClax}})});

            // M versus detector hits
            hname = base;
            hname.append("MvsdetHits_").append(res->name()).append(d1->name());
            annot = "M versus detector hits; M (" + res->name() + ") GeV/c^{2}; ITS + 2*TPC + 4*TRD + 8*TOF (" + d1->name() + ")";
            fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH2F, {max, {16, -0.5, 15.5}}})});
          } else {
            // M vs Mi
            hname = base;
            hname.append("MvsM_").append(res->name()).append(d1->name());
            annot = "M versus M; M (" + res->name() + ") GeV/c^{2}; M (" + d1->name() + ") GeV/c^{2}";
            auto max1 = o2::framework::AxisSpec(res->nmassBins(), d1->massHistRange()[0], d1->massHistRange()[1]);
            fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH2F, {max, max1}})});
          }
        }

        // daughters vs daughters
        for (auto i = 0; i < static_cast<int>(ndaughs - 1); i++) {
          auto d1 = getResonance(daughs[i]);
          auto max1 = o2::framework::AxisSpec(d1->nmassBins(), d1->massHistRange()[0], d1->massHistRange()[1]);
          for (auto j = i + 1; j < static_cast<int>(ndaughs); j++) {
            auto d2 = getResonance(daughs[j]);
            auto max2 = o2::framework::AxisSpec(d2->nmassBins(), d2->massHistRange()[0], d2->massHistRange()[1]);

            // M1 vs M2
            hname = base;
            hname.append("MvsM_").append(d1->name()).append(d2->name());
            annot = std::string("M versus M; M (").append(d1->name()).append(") GeV/c^{2}; M (").append(d2->name()).append(") GeV/c^{2}");
            fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH2F, {max1, max2}})});

            // angle(d1, d2)
            hname = base;
            hname.append("angle_").append(d1->name()).append(d2->name());
            annot = std::string("angle; Angle (").append(d1->name()).append(", ").append(d2->name()).append(")");
            fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH1F, {angax}})});

            // M vs angle(d1, d2)
            hname = base;
            hname.append("Mvsangle_").append(d1->name()).append(d2->name());
            annot = std::string("M versus angle; M (").append(res->name()).append(") GeV/c^{2}; Angle (").append(d1->name()).append(", ").append(d2->name()).append(")");
            fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH2F, {max, angax}})});

            // both daughters are finals
            if (d1->isFinal() && d2->isFinal()) {
              hname = base;
              hname.append("TPCsignal_").append(d1->name()).append(d2->name());
              annot = std::string("TPC signal of both tracks; TPCsignal (").append(d1->name()).append("); TPCsignal (").append(d2->name()).append(")");
              fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH2F, {sTPCax, sTPCax}})});
            }
          }
        }

        // for finals only
        if (res->isFinal()) {
          // dca
          hname = base;
          hname.append("dcaXY");
          annot = std::string("dcaXY; dca_{XY}(").append(res->name()).append(")");
          fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH1F, {dcaxyax}})});
          hname = base;
          hname.append("dcaZ");
          annot = std::string("dcaZ; dca_{Z}(").append(res->name()).append(")");
          fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH1F, {dcazax}})});

          // nSIgma[TPC, TOF] vs pT
          for (const auto& det : fdets) {
            for (const auto& part : fparts) {
              hname = base;
              hname.append("nS").append(part).append(det);
              annot = std::string("nSigma_").append(det).append(" versus p; p (").append(res->name()).append(") GeV/c; nSigma_{").append(det).append(", ").append(part).append("} (").append(res->name()).append(")");
              fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH2F, {momax, nSax}})});
            }
          }

          // detector hits
          hname = base;
          hname.append("detectorHits");
          annot = std::string("detectorHits; Detector(").append(res->name()).append(")");
          fhistPointers.insert({hname, registry.add(hname.c_str(), annot.c_str(), {o2::framework::HistType::kTH1F, {{4, 0.5, 4.5}}})});
        }
      }
    }
  }

  // ClassDefNV(decayTree, 1);
};

#endif // PWGUD_CORE_DECAYTREE_H_
