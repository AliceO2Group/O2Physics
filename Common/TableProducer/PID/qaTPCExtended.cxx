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
/// \file qaTPCExtended.cxx
/// \author Jeremy Wilkinson <jeremy.wilkinson@cern.ch>, GSI Darmstadt
/// \brief Task to give QA output for TPC PID response based on external classes (TOF cut, V0s, etc)


#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Common/Core/trackUtilities.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{         //runtime customisation goes here
  };
  std::swap(workflowOptions, options);
     
}

#include "Framework/runDataProcessing.h"

using namespace o2::dataformats;


struct QaTpcPid {
  //Configurables
  Configurable<float> cutTOF{"cutTOF", 3.f, "TOF nsigma cut for TPC-TOF PID"};
  Configurable<int> pBins{"pBins", 400, "Number of momentum bins"};
  Configurable<float> pMin{"pMin", 0.f, "Lower limit in momentum"};
  Configurable<float> pMax{"pMax", 20.f, "Upper limit in momentum"};
  Configurable<int> nBinsNSigma{"nBinsNSigma",200, "Number of bins for TPC nSigma"};
  Configurable<float> minNSigma{"minNSigma", -10.f, "Lower limit for TPC nSigma"};
  Configurable<float> maxNSigma{"maxNSigma",  10.f, "Upper limit for TPC nSigma"};
  
  HistogramRegistry hists{"HistogramsTPCPIDQA"};
  
  static constexpr int Np = 9;
  static constexpr std::string_view hnsigmaTPC[Np] = {"nsigmaTPC/El", "nsigmaTPC/Mu", "nsigmaTPC/Pi",
                                                   "nsigmaTPC/Ka", "nsigmaTPC/Pr", "nsigmaTPC/De",
                                                   "nsigmaTPC/Tr", "nsigmaTPC/He", "nsigmaTPC/Al"};
  static constexpr std::string_view hnsigmaTPCTOF[Np] = {"nsigmaTPCTOF/El", "nsigmaTPCTOF/Mu", "nsigmaTPCTOF/Pi",
                                                   "nsigmaTPCTOF/Ka", "nsigmaTPCTOF/Pr", "nsigmaTPCTOF/De",
                                                   "nsigmaTPCTOF/Tr", "nsigmaTPCTOF/He", "nsigmaTPCTOF/Al"};
                                                   
  static constexpr std::string_view hnsigmaTOF[Np] = {"nsigmaTOF/El", "nsigmaTOF/Mu", "nsigmaTOF/Pi",
                                                   "nsigmaTOF/Ka", "nsigmaTOF/Pr", "nsigmaTOF/De",
                                                   "nsigmaTOF/Tr", "nsigmaTOF/He", "nsigmaTOF/Al"};
  static constexpr std::string_view hnsigmaTOFAfter[Np] = {"nsigmaTOFAfter/El", "nsigmaTOFAfter/Mu", "nsigmaTOFAfter/Pi",
                                                   "nsigmaTOFAfter/Ka", "nsigmaTOFAfter/Pr", "nsigmaTOFAfter/De",
                                                   "nsigmaTOFAfter/Tr", "nsigmaTOFAfter/He", "nsigmaTOFAfter/Al"};
  static constexpr const char* partName[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};

  template <uint8_t i>
  void addTPCQAParticleHistos()
  {
    AxisSpec pAxis{pBins, pMin, pMax, "#it{p} [GeV/#it{c}]"};
    AxisSpec nSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, "TPC n_{#sigma}"};
    AxisSpec tofnSigmaAxis{nBinsNSigma, minNSigma, maxNSigma, "TOF n_{#sigma}"};
    
    hists.add(hnsigmaTPC[i].data(), Form("TPC signal (%s) without TOF cut",partName[i]), kTH2F, {pAxis, nSigmaAxis}); 
    hists.add(hnsigmaTPCTOF[i].data(), Form("TPC signal (%s) after TOF cut",partName[i],cutTOF), kTH2F, {pAxis, nSigmaAxis}); 
    hists.add(hnsigmaTOF[i].data(), "TOF signal", kTH2F, {pAxis, tofnSigmaAxis}); 
    hists.add(hnsigmaTOFAfter[i].data(), "TOF signal after TOF cut", kTH2F, {pAxis, tofnSigmaAxis}); 
    
  }//addParticleHistos
  

  void init(InitContext&)
  {
  
    addTPCQAParticleHistos<0>();
    addTPCQAParticleHistos<1>();
    addTPCQAParticleHistos<2>();
    addTPCQAParticleHistos<3>();
    addTPCQAParticleHistos<4>();
    addTPCQAParticleHistos<5>();
    addTPCQAParticleHistos<6>();
    addTPCQAParticleHistos<7>();
    addTPCQAParticleHistos<8>();
  
  }//init
  
  template <uint8_t i, typename T>
  void fillTPCQAParticleHistos(const T& t, const float mom, const float tofNSigma, const float tpcNSigma)
  {
    hists.fill(HIST(hnsigmaTPC[i]),mom, tpcNSigma);
    if (abs(tofNSigma) < cutTOF) hists.fill(HIST(hnsigmaTPCTOF[i]),mom, tpcNSigma);
    
    hists.fill(HIST(hnsigmaTOF[i]),mom, tofNSigma);
    if (abs(tofNSigma) < cutTOF) hists.fill(HIST(hnsigmaTOFAfter[i]),mom, tofNSigma);
  }//fillParticleHistos
  
  void process(aod::Collision const& collision, soa::Join<aod::Tracks, aod::TracksExtra,
                                                          aod::pidTPCFullEl, aod::pidTPCFullMu, aod::pidTPCFullPi,
                                                          aod::pidTPCFullKa, aod::pidTPCFullPr, aod::pidTPCFullDe,
                                                          aod::pidTPCFullTr, aod::pidTPCFullHe, aod::pidTPCFullAl,
                                                          aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                                          aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe,
                                                          aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl,
                                                          aod::TrackSelection> const& tracks)
  {
    for (auto t : tracks) {
      const float mom = t.tpcInnerParam();
      fillTPCQAParticleHistos<0>(t, mom, t.tofNSigmaEl(), t.tpcNSigmaEl());
      fillTPCQAParticleHistos<1>(t, mom, t.tofNSigmaMu(), t.tpcNSigmaMu());
      fillTPCQAParticleHistos<2>(t, mom, t.tofNSigmaPi(), t.tpcNSigmaPi());
      fillTPCQAParticleHistos<3>(t, mom, t.tofNSigmaKa(), t.tpcNSigmaKa());
      fillTPCQAParticleHistos<4>(t, mom, t.tofNSigmaPr(), t.tpcNSigmaPr());
      fillTPCQAParticleHistos<5>(t, mom, t.tofNSigmaDe(), t.tpcNSigmaDe());
      fillTPCQAParticleHistos<6>(t, mom, t.tofNSigmaTr(), t.tpcNSigmaTr());
      fillTPCQAParticleHistos<7>(t, mom, t.tofNSigmaHe(), t.tpcNSigmaHe());
      fillTPCQAParticleHistos<8>(t, mom, t.tofNSigmaAl(), t.tpcNSigmaAl());
      
    }//for
  }//process
};//struct QaTPCPID


WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec w;
  w.push_back(adaptAnalysisTask<QaTpcPid>(cfgc));
  return w;
}