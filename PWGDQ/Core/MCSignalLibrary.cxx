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
// Contact: iarsene@cern.ch, i.c.arsene@fys.uio.no
//
#include "PWGDQ/Core/MCSignalLibrary.h"

#include "MCProng.h"
#include "MCSignal.h"

#include <CommonConstants/PhysicsConstants.h>
#include <Framework/Logger.h>

#include <TPDGCode.h>
#include <TString.h>

#include <rapidjson/document.h>
#include <rapidjson/error/error.h>

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

using namespace o2::constants::physics;
// using std::cout;
// using std::endl;

MCSignal* o2::aod::dqmcsignals::GetMCSignal(const char* name)
{
  std::string nameStr = name;
  MCSignal* signal = nullptr;
  // 1-prong signals
  if (nameStr == "alicePrimary") {
    MCProng prong(1);                                  // 1-generation prong
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);  // set source to be ALICE primary particles
    signal = new MCSignal(1, name, "ALICE primaries"); // define a signal with one prong
    signal->AddProng(prong);                           // add the previously defined prong to the signal
    return signal;
  }
  if (nameStr == "electron") {
    MCProng prong(1, {11}, {true}, {false}, {0}, {0}, {false});        // define 1-generation prong using the full constructor
    signal = new MCSignal(name, "Inclusive electrons", {prong}, {-1}); // define the signal using the full constructor
    return signal;
  }
  if (nameStr == "electronPrimary") {
    MCProng prong(1, {11}, {true}, {false}, {0}, {0}, {false});      // define 1-generation prong using the full constructor
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);                // set source to be ALICE primary particles
    signal = new MCSignal(name, "Primary electrons", {prong}, {-1}); // define the signal using the full constructor
    return signal;
  }
  if (nameStr == "pionPrimary") {
    MCProng prong(1, {211}, {true}, {false}, {0}, {0}, {false}); // define 1-generation prong using the full constructor
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);            // set source to be ALICE primary particles
    signal = new MCSignal(name, "Primary pions", {prong}, {-1}); // define the signal using the full constructor
    return signal;
  }
  if (nameStr == "pionPrimaryFromHc") {
    MCProng prong(2, {211, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);
    signal = new MCSignal(name, "Primary pions from open charmed hadron decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "pionPrimaryFromHb") {
    MCProng prong(2, {211, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);
    signal = new MCSignal(name, "Primary pions from open beauty hadron decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "photon") {
    MCProng prong(1, {22}, {true}, {false}, {0}, {0}, {false}); // define 1-generation prong using the full constructor
    signal = new MCSignal(name, "Photon", {prong}, {-1});       // define the signal using the full constructor
    return signal;
  }
  if (nameStr == "muonPrimary") {
    MCProng prong(1, {13}, {true}, {false}, {0}, {0}, {false});  // define 1-generation prong using the full constructor
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);            // set source to be ALICE primary particles
    signal = new MCSignal(name, "Primary Muons", {prong}, {-1}); // define the signal using the full constructor
    return signal;
  }
  if (nameStr == "kaonFromPhi") {
    MCProng prong(2, {321, 333}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); // define 2-generation prong using the full constructor
    signal = new MCSignal(name, "Kaons from phi-mesons", {prong}, {-1});                        // define the signal using the full constructor
    return signal;
  }
  if (nameStr == "kaonPrimary") {
    MCProng prong(1, {321}, {true}, {false}, {0}, {0}, {false}); // define 1-generation prong using the full constructor
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);            // set source to be ALICE primary particles
    signal = new MCSignal(name, "Primary Kaons", {prong}, {-1}); // define the signal using the full constructor
    return signal;
  }
  if (nameStr == "Lambda0Baryon") {
    MCProng prong(1, {3122}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Lambda0 Baryon", {prong}, {-1});
    return signal;
  }
  if (nameStr == "SigmaPlusBaryon") {
    MCProng prong(1, {3222}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "SigmaPlus Baryon", {prong}, {-1});
    return signal;
  }
  if (nameStr == "proton") {
    MCProng prong(1, {2212}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "proton", {prong}, {-1});
    return signal;
  }
  if (nameStr == "protonPrimary") {
    MCProng prong(1, {2212}, {true}, {false}, {0}, {0}, {false}); // define 1-generation prong using the full constructor
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);             // set source to be ALICE primary particles
    signal = new MCSignal(name, "Primary Proton", {prong}, {-1}); // define the signal using the full constructor
    return signal;
  }
  if (nameStr == "protonFromTransport") {
    MCProng prong(1, {2212}, {true}, {false}, {0}, {0}, {false});
    prong.SetSourceBit(0, MCProng::kProducedInTransport);
    signal = new MCSignal(name, "ProtonFromTransport", {prong}, {-1});
    return signal;
  }
  if (nameStr == "protonFromLambda0") {
    MCProng prong(2, {2212, 3122}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Proton from Lambda0 decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "protonFromSigmaPlus") {
    MCProng prong(2, {2212, 3222}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Proton from Sigma+ decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "phiMeson") {
    MCProng prong(1, {333}, {true}, {false}, {0}, {0}, {false}); // define 1-generation prong using the full constructor
    signal = new MCSignal(name, "phi meson", {prong}, {-1});     // define the signal using the full constructor
    return signal;
  }
  if (nameStr == "muon") {
    MCProng prong(1, {13}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive muons", {prong}, {-1});
    return signal;
  }
  if (nameStr == "electronNOTfromTransport") {
    MCProng prong(1);
    prong.SetPDGcode(0, 11, true);
    prong.SetSourceBit(0, MCProng::kProducedInTransport, true); // exclude particles produces in transport
    signal = new MCSignal(name, "Electrons which are not produced during transport in detector", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromNonpromptJpsi") {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {false});
    signal = new MCSignal(name, "Electrons from non-prompt jpsi decays with beauty in decay chain", {prong}, {-1});
    return signal;
  }
  if (nameStr == "ePrimaryFromPromptJpsi") {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from prompt jpsi decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "ePrimaryFromNonpromptJpsi") {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from non-prompt jpsi decays with beauty in decay chain", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Jpsi") {
    MCProng prong(1, {443}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive jpsi", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Helium3") {
    MCProng prong(1, {1000020030}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Helium3", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Helium3Primary") {
    MCProng prong(1, {1000020030}, {true}, {false}, {0}, {0}, {false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Helium3Primary", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Helium3FromTransport") {
    MCProng prong(1, {1000020030}, {true}, {false}, {0}, {0}, {false});
    prong.SetSourceBit(0, MCProng::kProducedInTransport);
    signal = new MCSignal(name, "Helium3FromTransport", {prong}, {-1});
    return signal;
  }
  if (nameStr == "promptJpsi") {
    MCProng prong(1, {443}, {true}, {false}, {0}, {0}, {false}, false, {503}, {true});
    signal = new MCSignal(name, "Prompt jpsi (not from beauty)", {prong}, {-1});
    return signal;
  }
  if (nameStr == "nonPromptJpsi") {
    MCProng prong(1, {443}, {true}, {false}, {0}, {0}, {false}, false, {503}, {false});
    signal = new MCSignal(name, "Non-prompt jpsi (from beauty)", {prong}, {-1});
    return signal;
  }
  if (nameStr == "nonPromptJpsiFromBeauty") {
    MCProng prong(2, {503, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true);
    signal = new MCSignal(name, "Non-prompt jpsi directly from beauty", {prong}, {-1});
    return signal;
  }
  if (nameStr == "nonPromptJpsiNotDirectlyFromBeauty") {
    MCProng prong(2, {443, 503}, {true, true}, {false, true}, {0, 0}, {0, 0}, {false, false}, false, {503}, {false});
    signal = new MCSignal(name, "Non-prompt jpsi from other but with beauty in decay chain", {prong}, {-1});
    return signal;
  }
  if (nameStr == "AnythingDecayToJpsi") {
    MCProng prong(2, {MCProng::kPDGCodeNotAssigned, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true);
    signal = new MCSignal(name, "Decay of anything into J/psi", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eeFromNonpromptPsi2S") {
    MCProng prong(2, {11, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {false});
    signal = new MCSignal(name, "ee pairs from non-prompt psi2s decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eeFromPromptPsi2S") {
    MCProng prong(2, {11, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {true});
    signal = new MCSignal(name, "ee pairs from prompt psi2s decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eFromNonpromptPsi2S") {
    MCProng prong(2, {11, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {false});
    signal = new MCSignal(name, "Electrons from beauty psi2s decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromPromptPsi2S") {
    MCProng prong(2, {11, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {true});
    signal = new MCSignal(name, "Electrons from prompt psi2s decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Psi2S") {
    MCProng prong(1, {100443}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive psi2s", {prong}, {-1});
    return signal;
  }
  if (nameStr == "nonPromptPsi2S") {
    MCProng prong(1, {100443}, {true}, {false}, {0}, {0}, {false}, false, {503}, {false});
    signal = new MCSignal(name, "Non-prompt psi2s", {prong}, {-1});
    return signal;
  }
  if (nameStr == "promptPsi2S") {
    MCProng prong(1, {100443}, {true}, {false}, {0}, {0}, {false}, false, {503}, {true});
    signal = new MCSignal(name, "Prompt psi2s (not from beauty)", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Chic0") {
    MCProng prong(1, {10441}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive Chic0", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Chic1") {
    MCProng prong(1, {20443}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive Chic1", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Chic2") {
    MCProng prong(1, {445}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive Chic2", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Chic012") {
    MCProng prong(1, {904}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive Chic0, Chic1 and Chic2", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Upsilon1S") {
    MCProng prong(1, {553}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive Upsilon1S", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Upsilon2S") {
    MCProng prong(1, {100553}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive Upsilon2S", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Upsilon3S") {
    MCProng prong(1, {200553}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive Upsilon3S", {prong}, {-1});
    return signal;
  }
  if (nameStr == "allBeautyHadrons") {
    MCProng prong(1, {503}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "All beauty hadrons", {prong}, {-1});
    return signal;
  }
  if (nameStr == "allBeautyHadronsFS") {
    MCProng prong(1, {503}, {true}, {false}, {0}, {0}, {false});
    prong.SetSourceBit(0, MCProng::kHEPMCFinalState);
    signal = new MCSignal(name, "All beauty hadrons", {prong}, {-1});
    return signal;
  }
  if (nameStr == "allOpenBeautyHadrons") {
    MCProng prong(1, {502}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "All open beauty hadrons", {prong}, {-1});
    return signal;
  }
  if (nameStr == "allOpenBeautyHadronsFS") {
    MCProng prong(1, {502}, {true}, {false}, {0}, {0}, {false});
    prong.SetSourceBit(0, MCProng::kHEPMCFinalState);
    signal = new MCSignal(name, "All open beauty hadrons", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Bc") {
    MCProng prong(1, {541}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Bc", {prong}, {-1});
    return signal;
  }
  if (nameStr == "mumuFromJpsiFromBc") {
    MCProng prong(3, {13, 443, 541}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Muon pair from jpsi from Bc decays", {prong, prong}, {1, 1});
    return signal;
  }
  if (nameStr == "muFromBc") {
    MCProng prong(2, {13, 541}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Muon from Bc decays", {prong}, {1});
    return signal;
  }
  if (nameStr == "mumumuFromBc") {
    MCProng prongMuFromJpsi(3, {13, 443, 541}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongMuFromBc(2, {13, 541}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Trimuon from Bc decays", {prongMuFromJpsi, prongMuFromJpsi, prongMuFromBc}, {2, 2, 1});
    return signal;
  }
  if (nameStr == "everythingFromBeauty") {
    MCProng prong(2, {0, 503}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Everything from beauty", {prong}, {-1});
    return signal;
  }
  if (nameStr == "everythingFromBeautyFS") {
    MCProng prong(2, {0, 503}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(1, MCProng::kHEPMCFinalState);
    signal = new MCSignal(name, "Everything from beauty", {prong}, {-1});
    return signal;
  }
  if (nameStr == "everythingFromEverythingFromBeauty") {
    MCProng prong(3, {0, 0, 503}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Everything from everything from beauty", {prong}, {-1});
    return signal;
  }
  if (nameStr == "everythingFromEverythingFromBeautyFS") {
    MCProng prong(3, {0, 0, 503}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSourceBit(2, MCProng::kHEPMCFinalState);
    signal = new MCSignal(name, "Everything from everything from beauty", {prong}, {-1});
    return signal;
  }
  if (nameStr == "allCharmHadrons") {
    MCProng prong(1, {403}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "All charm hadrons", {prong}, {-1});
    return signal;
  }
  if (nameStr == "allOpenCharmHadrons") {
    MCProng prong(1, {402}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "All open charm hadrons", {prong}, {-1});
    return signal;
  }
  if (nameStr == "allCharmFromBeauty") {
    MCProng prong(2, {403, 503}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "All charm hadrons from beauty", {prong}, {-1});
    return signal;
  }
  if (nameStr == "allPromptCharm") {
    MCProng prong(2, {403, 503}, {true, true}, {false, true}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "All prompt charm hadrons (not from beauty)", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Pi0DecayToe") {
    MCProng prong(2, {111, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Pi0 decays into an electron", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Pi0DecayTog") {
    MCProng prong(2, {111, 22}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Pi0 decays into an gamma", {prong}, {1});
    return signal;
  }
  if (nameStr == "Pi0DecayTogg") {
    MCProng prong(2, {111, 22}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Pi0 decays into an gamma, gamma", {prong, prong}, {1, 1});
    return signal;
  }
  if (nameStr == "PromptPi0DecayToe") {
    MCProng prong(2, {111, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, true, {403, 503}, {true, true});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Pi0 decays into an electron", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Pi0") {
    MCProng prong(1, {111}, {true}, {false}, {0}, {0}, {false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
    signal = new MCSignal(name, "Pi0", {prong}, {-1});
    return signal;
  }
  if (nameStr == "LMeeLFQ") {
    MCProng prong(1, {900}, {true}, {false}, {0}, {0}, {false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
    signal = new MCSignal(name, "light flavor mesons + quarkonia", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi,jpsi,psi2s
    return signal;
  }
  if (nameStr == "LMeeLF") {
    MCProng prong(1, {901}, {true}, {false}, {0}, {0}, {false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
    signal = new MCSignal(name, "ligh flavor mesons", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi
    return signal;
  }
  if (nameStr == "PromptJpsiDecayToe") {
    MCProng prong(2, {443, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, true, {503}, {true});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Prompt jpsi (not from beauty) decay to electron", {prong}, {-1});
    return signal;
  }
  if (nameStr == "electronFromDs") {
    MCProng prong(2, {11, 431}, {true, true}, {false, false}, {0, 0}, {0, 0}, {true, true});
    signal = new MCSignal(name, "Electrons from Ds decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "dsMeson") {
    MCProng prong(1, {431}, {true}, {false}, {0}, {0}, {true});
    signal = new MCSignal(name, "Ds mesons", {prong}, {-1});
    return signal;
  }
  if (nameStr == "electronFromPC") {
    MCProng prong(2, {11, 22}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "electron from a photon conversion", {prong}, {-1});
    return signal;
  }
  if (nameStr == "PowhegDYMuon1") {
    MCProng prong(1, {13}, {true}, {false}, {0}, {0}, {false});        // selecting muons
    prong.SetSourceBit(0, MCProng::kIsPowhegDYMuon);                   // set source to be Muon from POWHEG
    signal = new MCSignal(name, "POWHEG Muon singles", {prong}, {-1}); // define a signal with 1-prong
    return signal;
  }

  // 2-prong signals
  if (nameStr == "dielectron") {
    MCProng prong(1, {11}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Electron pair", {prong, prong}, {-1, -1});
    return signal;
  }
  if (nameStr == "dimuon") {
    MCProng prong(1, {13}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Muon pair", {prong, prong}, {-1, -1});
    return signal;
  }
  if (nameStr == "electronMuonPair") {
    MCProng electron(1, {11}, {true}, {false}, {0}, {0}, {false});
    MCProng muon(1, {13}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Electron-muon pair", {electron, muon}, {-1, -1});
    return signal;
  }
  if (nameStr == "emuFromOpenHFhadron") {
    MCProng electron(2, {11, 902}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    electron.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng muon(2, {13, 902}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    muon.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "e and mu each from an open charm or beauty hadron decay", {electron, muon}, {-1, -1});
    return signal;
  }
  if (nameStr == "emuFromOpenCharmHadron") {
    MCProng electron(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    electron.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng muon(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    muon.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "e and mu each from an open charm hadron decay", {electron, muon}, {-1, -1});
    return signal;
  }
  if (nameStr == "emuFromOpenBeautyHadron") {
    MCProng electron(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    electron.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng muon(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    muon.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "e and mu each from an open beauty hadron decay", {electron, muon}, {-1, -1});
    return signal;
  }
  if (nameStr == "dielectronFromPC") {
    MCProng prong(2, {11, 22}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "dielectron from a photon conversion", {prong, prong}, {1, 1});
    return signal;
  }
  if (nameStr == "dielectronFromAllPC") {
    MCProng prong(2, {11, 22}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "dielectron from a photon conversion", {prong, prong}, {-1, -1});
    return signal;
  }
  if (nameStr == "dielectronPCPi0") {
    MCProng prong(3, {11, 22, 111}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "dielectron from a photon conversion from a pi0", {prong, prong}, {1, 1});
    return signal;
  }
  if (nameStr == "PowhegDYMuon2") {
    MCProng prong(1, {13}, {true}, {false}, {0}, {0}, {false});                // selecting muons
    prong.SetSourceBit(0, MCProng::kIsPowhegDYMuon);                           // set source to be Muon from POWHEG
    signal = new MCSignal(name, "POWHEG Muon pair", {prong, prong}, {-1, -1}); // define a signal with 2-prong
    return signal;
  }

  // LMEE single signals
  // electron signals with mother X: e from mother X
  if (nameStr == "eFromAnything") {
    MCProng prong(2, {11, MCProng::kPDGCodeNotAssigned}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from any mother", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromPhoton") {
    MCProng prong(2, {11, 22}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from photon conversion", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromPi0") {
    MCProng prong(2, {11, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from pi0 decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "ePrimaryFromPromptPi0") {
    MCProng prong(2, {11, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502, 402}, {true, true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from prompt pi0 decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromEta") {
    MCProng prong(2, {11, 221}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from eta decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromEtaPrime") {
    MCProng prong(2, {11, 331}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from eta' decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromRho") {
    MCProng prong(2, {11, 113}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from rho decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromOmega") {
    MCProng prong(2, {11, 223}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from omega decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromPhi") {
    MCProng prong(2, {11, 333}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from phi decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromJpsi") {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from jpsi decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "anythingFromJpsi") {
    MCProng prong(2, {MCProng::kPDGCodeNotAssigned, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Anything from jpsi decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromPromptJpsi") {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {true});
    signal = new MCSignal(name, "Electrons from jpsi decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromPsi2S") {
    MCProng prong(2, {11, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from psi2s decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromLMeeLF") {
    MCProng prong(2, {11, 901}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from LF meson decays", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi mesons
    return signal;
  }
  if (nameStr == "ePrimaryFromLMeeLF") {
    MCProng prong(2, {11, 901}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);                      // set source to be ALICE primary particles
    signal = new MCSignal(name, "Electrons from LF meson decays", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi mesons
    return signal;
  }
  if (nameStr == "eFromLMeeLFQ") {
    MCProng prong(2, {11, 900}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(false);                                                             // set direction to check generation in time (true) or back in time (false)
    signal = new MCSignal(name, "Electrons from LF meson + quarkonia decays", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi,jpsi,psi2s mesons
    return signal;
  }
  if (nameStr == "ePrimaryFromLMeeLFQ") {
    MCProng prong(2, {11, 900}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);                                  // set source to be ALICE primary particles
    prong.SetSignalInTime(false);                                                             // set direction to check generation in time (true) or back in time (false)
    signal = new MCSignal(name, "Electrons from LF meson + quarkonia decays", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi,jpsi,psi2s mesons
    return signal;
  }
  if (nameStr == "eFromHc") {
    MCProng prong(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open charmed hadron decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromAnyHc") {
    MCProng prong(1, {11}, {true}, {false}, {0}, {0}, {false}, false, {402}, {false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);
    signal = new MCSignal(name, "Electrons from any open charm hadron decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromD0") {
    MCProng prong(2, {11, Pdg::kD0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from D0 decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromChargedD") {
    MCProng prong(2, {11, Pdg::kDPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from D+/- decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromDs") {
    MCProng prong(2, {11, Pdg::kDS}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from Ds +/- decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromLambdaC") {
    MCProng prong(2, {11, Pdg::kLambdaCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from Lambda_c decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromXiC0") {
    MCProng prong(2, {11, Pdg::kXiC0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from Xi_c_0 decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromXiCPlus") {
    MCProng prong(2, {11, Pdg::kXiCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from Xi_c_+ decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromXiCPlusPlus") {
    MCProng prong(2, {11, Pdg::kXiCCPlusPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from Xi_c_++ decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromHb") {
    MCProng prong(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open beauty hadron decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromAnyHb") {
    MCProng prong(1, {11}, {true}, {false}, {0}, {0}, {false}, false, {502}, {false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);
    signal = new MCSignal(name, "Electrons from any open beauty hadron decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromHbc") {
    MCProng prong(2, {11, 902}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open charm or beauty hadron decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromMc") {
    MCProng prong(2, {11, 401}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open charmed meson decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromMb") {
    MCProng prong(2, {11, 501}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open beauty meson decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromBc") {
    MCProng prong(2, {11, 4001}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open charmed baryon decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromBb") {
    MCProng prong(2, {11, 5001}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open beauty baryon decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromPromptHc") {
    MCProng prong(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open charmed hadron decays without beauty in decay history", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromHbtoHc") {
    MCProng prong(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open charmed hadron decays with b hadron in decay history", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromPromptLM") {
    MCProng prong(2, {11, 101}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502, 402}, {true, true});
    signal = new MCSignal(name, "Electrons from light mesons without B/D in decay history", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromHbtoLM") {
    MCProng prong(2, {11, 101}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {false});
    signal = new MCSignal(name, "Electrons from light mesons with B hadron in decay history", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromHctoLM") {
    MCProng prong(2, {11, 101, 402}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false, {502}, {true});
    signal = new MCSignal(name, "Electrons from light mesons from D hadron decays and no B in decay history", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromUpsilon1S") {
    MCProng prong(2, {11, 553}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from Upsilon1S decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromUpsilon2S") {
    MCProng prong(2, {11, 100553}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from Upsilon2S decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "eFromUpsilon3S") {
    MCProng prong(2, {11, 200553}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from Upsilon3S decays", {prong}, {-1});
    return signal;
  }

  // muon signals with mother X: mu from mother X
  if (nameStr == "muFromJpsi") {
    MCProng prong(2, {13, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "muons from jpsi decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "muFromPsi2S") {
    MCProng prong(2, {13, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "muons from psi2s decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "muFromHb") {
    MCProng prong(2, {13, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "muons from b->mu", {prong}, {-1});
    return signal;
  }
  if (nameStr == "muFromPromptHc") {
    MCProng prong(2, {13, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    signal = new MCSignal(name, "muons from c->mu, without beauty in decay history", {prong}, {-1});
    return signal;
  }
  if (nameStr == "muFromHbtoHc") {
    MCProng prong(3, {13, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "muons from b->c->mu", {prong}, {-1});
    return signal;
  }
  if (nameStr == "secondaryMuon") {
    MCProng prong(1, {13}, {true}, {false}, {0}, {0}, {false});
    prong.SetSourceBit(0, MCProng::kProducedInTransport);
    signal = new MCSignal(name, "muons produced during transport in detector", {prong}, {-1});
    return signal;
  }
  if (nameStr == "muFromPromptLM") {
    MCProng prong(2, {13, 101}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502, 402}, {true, true});
    signal = new MCSignal(name, "muons from light mesons without B/D in decay history", {prong}, {-1});
    return signal;
  }
  if (nameStr == "muFromHbtoLM") {
    MCProng prong(2, {13, 101}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {false});
    signal = new MCSignal(name, "muons from light mesons with B hadron in decay history", {prong}, {-1});
    return signal;
  }
  if (nameStr == "muFromHctoLM") {
    MCProng prong(2, {13, 101, 402}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false, {502}, {true});
    signal = new MCSignal(name, "muons from light mesons from D hadron decays and no B in decay history", {prong}, {-1});
    return signal;
  }
  if (nameStr == "muFromUpsilon1S") {
    MCProng prong(2, {13, 553}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "muons from Upsilon1S decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "muFromUpsilon2S") {
    MCProng prong(2, {13, 100553}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "muons from Upsilon2S decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "muFromUpsilon3S") {
    MCProng prong(2, {13, 200553}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "muons from Upsilon3S decays", {prong}, {-1});
    return signal;
  }

  // Decay signal: Mother to electron: X -> e
  if (nameStr == "AnythingToE") {
    MCProng prong(2, {MCProng::kPDGCodeNotAssigned, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
    prong.SetSignalInTime(true);                                            // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Decay of anything into e", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi,jpsi,psi2s mesons
    return signal;
  }
  if (nameStr == "LFQdecayToE") {
    MCProng prong(2, {900, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
    prong.SetSignalInTime(true);                                                       // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "LF meson  + quarkonia decays into e", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi,jpsi,psi2s mesons
    return signal;
  }
  if (nameStr == "HcToE") {
    MCProng prong(2, {402, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open charmed hadron decays into e", {prong}, {-1});
    return signal;
  }
  if (nameStr == "HbToE") {
    MCProng prong(2, {502, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open beauty hadron decays into e", {prong}, {-1});
    return signal;
  }
  if (nameStr == "McToE") {
    MCProng prong(2, {401, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open charmed meson decays into e", {prong}, {-1});
    return signal;
  }
  if (nameStr == "MbToE") {
    MCProng prong(2, {501, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open beauty meson decays into e", {prong}, {-1});
    return signal;
  }
  if (nameStr == "BcToE") {
    MCProng prong(2, {4001, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open charmed baryon decays into e", {prong}, {-1});
    return signal;
  }
  if (nameStr == "BbToE") {
    MCProng prong(2, {5001, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open beauty baryon decay into e", {prong}, {-1});
    return signal;
  }
  if (nameStr == "PromptHcToE") {
    MCProng prong(3, {502, 402, 11}, {true, true, true}, {true, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open charmed hadron decays into e", {prong}, {-1});
    return signal;
  }
  if (nameStr == "NonPromptHcToE") {
    MCProng prong(3, {502, 402, 11}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "b hadron decays to open charmed hadron decays to e", {prong}, {-1});
    return signal;
  }
  if (nameStr == "HFdecayToE") {
    MCProng prong(2, {902, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open charm and beauty to electrons", {prong}, {-1});
    return signal;
  }
  if (nameStr == "AnyHFdecayToE") {
    MCProng prong(1, {902}, {true}, {false}, {0}, {0}, {false}, true, {11}, {false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open charm and beauty to electrons", {prong}, {-1});
    return signal;
  }
  //   if (nameStr == "LFQtoPC") {
  //   MCProng prong(3, {900, 22, 11}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
  //   prong.SetSignalInTime(true);  // set direction to check for daughters (true, in time) or for mothers (false, back in time)
  //   signal = new MCSignal(name, "LF meson + quarkonia decays into photon conversion electron", {prong}, {-1});
  //   return signal;
  // }

  //_________________________________________________________________________________________________________________________
  // LMEE pair signals for LF, same mother
  if (nameStr == "eeFromAnything") {
    MCProng prong(2, {11, MCProng::kPDGCodeNotAssigned}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from any decay", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eeFromPi0") {
    MCProng prong(2, {11, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from pi0 decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eePrimaryFromPromptPi0") {
    MCProng prong(2, {11, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502, 402}, {true, true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from prompt pi0 decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eeFromEta") {
    MCProng prong(2, {11, 221}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from eta decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eeFromEtaprime") {
    MCProng prong(2, {11, 331}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from eta' decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eeFromRho") {
    MCProng prong(2, {11, 113}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from rho decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eeFromOmega") {
    MCProng prong(2, {11, 223}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from omega decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eeFromPhi") {
    MCProng prong(2, {11, 333}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from phi decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eeFromJpsi") {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from j/psi decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eeFromJpsiExclusive") {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from j/psi decays", {prong, prong}, {1, 1}); // signal at pair level
    signal->SetDecayChannelIsExclusive(2, true);
    return signal;
  }
  if (nameStr == "eeFromJpsiNotExclusive") {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from j/psi decays", {prong, prong}, {1, 1}); // signal at pair level
    signal->SetDecayChannelIsNotExclusive(2, true);
    return signal;
  }
  if (nameStr == "eePrimaryFromPromptJPsi") {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from prompt j/psi decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eePrimaryFromNonPromptJPsi") {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from non-prompt j/psi decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "mumuFromPhi") {
    MCProng prong(2, {13, 333}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "mumu pairs from phi decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "mumuFromJpsi") {
    MCProng prong(2, {13, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "mumu pairs from j/psi decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "mumuFromPromptJpsi") {
    MCProng prong(2, {13, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "mumu pairs from prompt j/psi decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "mumuFromNonPromptJpsi") {
    MCProng prong(2, {13, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "mumu pairs from non-prompt j/psi decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eeFromPsi2S") {
    MCProng prong(2, {11, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from psi2s decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "mumuFromPsi2S") {
    MCProng prong(2, {13, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "mumu pairs from psi2s decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "mumuFromPromptPsi2S") {
    MCProng prong(2, {13, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "mumu pairs from prompt psi2s decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "mumuFromNonPromptPsi2S") {
    MCProng prong(2, {13, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "mumu pairs from non-prompt psi2s decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "mumuFromUpsilon1S") {
    MCProng prong(2, {13, 553}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "mumu pairs from upsilon1s decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "mumuFromUpsilon2S") {
    MCProng prong(2, {13, 100553}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "mumu pairs from upsilon2s decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "mumuFromUpsilon3S") {
    MCProng prong(2, {13, 200553}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "mumu pairs from upsilon3s decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eeFromLMeeLFQ") {
    MCProng prong(2, {11, 900}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from light flavor meson + quarkonia decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eeFromLMeeLF") {
    MCProng prong(2, {11, 901}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from light flavor meson decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eeFromLMeeNoHFLFQ") {
    MCProng prong(2, {11, 900}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502, 402}, {true, true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from light flavor meson + quarkonia decays not from open-HF decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (nameStr == "eeFromLMeeNoHFLF") {
    MCProng prong(2, {11, 901}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502, 402}, {true, true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from light flavor meson decays not from open-HF decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }

  // LMEE pair signals for HF
  // D0->e and D0->e
  if (nameStr == "eeFromD0") {
    MCProng prong(2, {kElectron, Pdg::kD0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from D0 decays, no beauty in history", {prong, prong}, {-1, -1});
    return signal;
  }

  // D0->e and D0->e
  if (nameStr == "eeFromPi0FromD0") {
    MCProng prong(2, {kElectron, kPi0, Pdg::kD0}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from D0 to Pi0 decays, no beauty in history", {prong, prong}, {1, 1});
    return signal;
  }

  // D+/- -> e and D+/- -> e
  if (nameStr == "eeFromChargedD") {
    MCProng prong(2, {kElectron, Pdg::kDPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from D+/- decays, no beauty in history", {prong, prong}, {-1, -1});
    return signal;
  }

  // D0 -> e and D+/- -> e
  if (nameStr == "eeFromD0andChargedD") {
    MCProng prongD0(2, {kElectron, Pdg::kD0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    MCProng prongDch(2, {kElectron, Pdg::kDPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongD0.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongDch.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "one e from D0 and one e from D+/- decays, no beauty in history", {prongD0, prongDch}, {-1, -1});
    return signal;
  }

  // D+/- -> e and D0 -> e
  if (nameStr == "eeFromD0andChargedDBis") {
    MCProng prongD0(2, {kElectron, Pdg::kD0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    MCProng prongDch(2, {kElectron, Pdg::kDPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongD0.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongDch.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "one e from D+/- and one e from D0 decays, no beauty in history (inverse signal)", {prongDch, prongD0}, {-1, -1});
    return signal;
  }

  // D_s->e and D_s->e
  if (nameStr == "eeFromDs") {
    MCProng prong(2, {kElectron, Pdg::kDS}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from Ds +/- decays, no beauty in history", {prong, prong}, {-1, -1});
    return signal;
  }

  // Lambda_c->e and Lambda_c->e
  if (nameStr == "eeFromLambdaC") {
    MCProng prong(2, {kElectron, Pdg::kLambdaCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from Lambda_c, no beauty in history", {prong, prong}, {-1, -1});
    return signal;
  }

  // Lambda_c->e and D0->e
  if (nameStr == "eeFromLambdaCandD0") {
    MCProng prongLc(2, {kElectron, Pdg::kLambdaCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    MCProng prongD0(2, {kElectron, Pdg::kD0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongLc.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongD0.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "one e from Lambda_c and one e from D0 decays, no beauty in history", {prongLc, prongD0}, {-1, -1});
    return signal;
  }

  // D0->e and Lambda_c->e
  if (nameStr == "eeFromLambdaCandD0Bis") {
    MCProng prongLc(2, {kElectron, Pdg::kLambdaCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    MCProng prongD0(2, {kElectron, Pdg::kD0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongLc.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongD0.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "one e from Lambda_c and one e from D0 decays, no beauty in history (inverse signal)", {prongD0, prongLc}, {-1, -1});
    return signal;
  }

  // Lambda_c->e and D+/- -> e
  if (nameStr == "eeFromLambdaCandChargedD") {
    MCProng prongLc(2, {kElectron, Pdg::kLambdaCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    MCProng prongDch(2, {kElectron, Pdg::kDPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongLc.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongDch.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "one e from Lambda_c and one e from D+/- decays, no beauty in history", {prongLc, prongDch}, {-1, -1});
    return signal;
  }

  // D+/- -> e and Lambda_c->e
  if (nameStr == "eeFromLambdaCandChargedDBis") {
    MCProng prongLc(2, {kElectron, Pdg::kLambdaCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    MCProng prongDch(2, {kElectron, Pdg::kDPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongLc.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongDch.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "one e from Lambda_c and one e from D+/- decays, no beauty in history (inverse signal)", {prongDch, prongLc}, {-1, -1});
    return signal;
  }

  // Xic0 ->e and Xic0 ->e
  if (nameStr == "eeFromXiC0") {
    MCProng prong(2, {kElectron, Pdg::kXiC0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from Xi_c0, no beauty in history", {prong, prong}, {-1, -1});
    return signal;
  }

  // Xi_c+ ->e and Xi_c+ ->e
  if (nameStr == "eeFromXiCPlus") {
    MCProng prong(2, {kElectron, Pdg::kXiCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from Xi_c+, no beauty in history", {prong, prong}, {-1, -1});
    return signal;
  }

  // Xi_c0 ->e and Xi_c+ ->e
  if (nameStr == "eeFromXiC0andXiCPlus") {
    MCProng prongXiCPlus(2, {kElectron, Pdg::kXiCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    MCProng prongXiC0(2, {kElectron, Pdg::kXiC0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongXiCPlus.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongXiC0.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "one e from Xi_c+ and one e from Xi_c0 decays, no beauty in history", {prongXiCPlus, prongXiC0}, {-1, -1});
    return signal;
  }

  // Xi_c+ ->e and Xi_c0 ->e
  if (nameStr == "eeFromXiC0andXiCPlusBis") {
    MCProng prongXiCPlus(2, {kElectron, Pdg::kXiCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    MCProng prongXiC0(2, {kElectron, Pdg::kXiC0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongXiCPlus.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongXiC0.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "one e from Xi_c+ and one e from Xi_c0 decays, no beauty in history (inverse signal)", {prongXiC0, prongXiCPlus}, {-1, -1});
    return signal;
  }

  // Xi_cc++ ->e and Xi_cc++ ->e
  if (nameStr == "eeFromXiCPlusPlus") {
    MCProng prong(2, {kElectron, Pdg::kXiCCPlusPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from Xi_cc++, no beauty in history", {prong, prong}, {-1, -1});
    return signal;
  }

  // c->e and c->e (no check)
  if (nameStr == "eeFromCCNoCheck") {
    MCProng prong(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from c->e and c->e without check", {prong, prong}, {-1, -1}); // signal at pair level
    return signal;
  }

  // ee from HF in general
  if (nameStr == "eeFromHF") {
    MCProng prong(2, {11, 902}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b,c->e and b,c->e without check", {prong, prong}, {-1, -1}); // signal at pair level
    return signal;
  }

  // Any c in history but no b -> c -> e
  if (nameStr == "eeFromPromptCandPromptC") {
    MCProng prong(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true}); // check if mother pdg code is in history
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs with any charm but no beauty in decay chain", {prong, prong}, {-1, -1}); // signal at pair level
    return signal;
  }

  // Any b to any c in history b -> c -> e
  if (nameStr == "eeFromAnyBtoCandAnyBtoC") {
    MCProng prong(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {false}); // check if mother pdg code is in history
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs with any beauty to charm in decay chain", {prong, prong}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->c->e, b->c->e
  if (nameStr == "eeFromBtoCandBtoC") {
    MCProng prong(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false); // check if mother pdg code is in history
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs with any beauty to charm in decay chain", {prong, prong}, {-1, -1}); // signal at pair level
    return signal;
  }

  // Any b->e and Any b->X->c->e
  // Looking at such decays: B -> (e) D -> (e)e and bar{B} -> e
  // Signal allows combinations of ee from the same B meson
  // + the combination of e fom B and e from bar{B}
  if (nameStr == "eeFromBandAnyBtoC") {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {false}); // check if mother pdg code is in history
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->e and b->X->c->e", {prongB, prongBtoC}, {-1, -1}); // signal at pair level
    return signal;
  }

  // Any b->e and Any b->X->c->e
  if (nameStr == "eeFromBandAnyBtoCBis") {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {false}); // check if mother pdg code is in history
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->X->c->e and b->e", {prongBtoC, prongB}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->e and b->c->e
  if (nameStr == "eeFromBandBtoC") {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false); // check if mother pdg code is in history
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "direkt ee pairs from b->e and b->c->e", {prongB, prongBtoC}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->e and b->c->e
  if (nameStr == "eeFromBandBtoCBis") {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false); // check if mother pdg code is in history
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->c->e and b->e", {prongBtoC, prongB}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->e and b->c->e (same mother/grandmother)
  // require that the mother is the grandmother of the other electron
  if (nameStr == "eeFromBandBtoCsameGM") {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false); // check if mother pdg code is in history
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->e and b->c->e, mother = grandmother", {prongB, prongBtoC}, {1, 2}, false); // signal at pair level, accept commonAncestor Pairs
    return signal;
  }

  // b->e and b->c->e (same mother/grandmother)
  // require that the mother is the grandmother of the other electron
  if (nameStr == "eeFromBandBtoCsameGMBis") {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false); // check if mother pdg code is in history
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->c->e and b->e, mother = grandmother", {prongBtoC, prongB}, {2, 1}, false); // signal at pair level, accept commonAncestor Pairs
    return signal;
  }

  // b->e and b->c->e (different mother/grandmother)
  // require that the mother is not the grandmother of the other electron
  if (nameStr == "eeFromBandBtoCdiffGM") {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false); // check if mother pdg code is in history
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->e and b->c->e, mother != grandmother", {prongB, prongBtoC}, {1, 2}, true); // signal at pair level, exclude commonAncestor Pairs
    return signal;
  }

  // b->e and b->c->e (different mother/grandmother)
  // require that the mother is not the grandmother of the other electron
  if (nameStr == "eeFromBandBtoCdiffGMBis") {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false); // check if mother pdg code is in history
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->c->e and b->e, mother != grandmother", {prongBtoC, prongB}, {2, 1}, true); // signal at pair level, exclude commonAncestor Pairs
    return signal;
  }

  // b->e and b->e
  if (nameStr == "eeFromBB") {
    MCProng prong(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->e and b->e", {prong, prong}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->e and b->e (commonAncestors)
  if (nameStr == "eeFromSameB") {
    MCProng prong(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->e and b->e", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }

  // b->e and c->e no check
  if (nameStr == "eeFromBandFromC") {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongC(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->e and c->e", {prongB, prongC}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->e and c->e no check
  if (nameStr == "eeFromBandFromCBis") {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongC(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->e and c->e", {prongC, prongB}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->e and b->c->e (single b)
  if (nameStr == "eeFromSingleBandBtoC") {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->e and b->c->e (single b)", {prongB, prongBtoC}, {1, 2}); // signal at pair level
    return signal;
  }

  //_________________________________________________________________________________________________________________________

  if (nameStr == "kaonFromBplus") {
    MCProng prong(2, {321, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaons from B+ decays", {prong}, {1});
    return signal;
  }

  if (nameStr == "kaonFromBplusHistory") {
    MCProng prong(1, {321}, {true}, {false}, {0}, {0}, {false}, false, {521}, {false});
    signal = new MCSignal(name, "Kaons from B+ decays", {prong}, {-1});
    return signal;
  }

  if (nameStr == "kaonPrimaryFromBplusHistory") {
    MCProng prong(1, {321}, {true}, {false}, {0}, {0}, {false}, false, {521}, {false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Kaons from B+ decays", {prong}, {-1});
    return signal;
  }

  if (nameStr == "kaonPrimaryFromBplusFS") {
    MCProng prong(2, {321, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prong.SetSourceBit(1, MCProng::kHEPMCFinalState);
    signal = new MCSignal(name, "Kaons from B+ decays", {prong}, {-1});
    return signal;
  }

  if (nameStr == "kaonFromAnyBHistory") {
    MCProng prong(1, {321}, {true}, {false}, {0}, {0}, {false}, false, {503}, {false});
    signal = new MCSignal(name, "Kaons from B+ decays", {prong}, {-1});
    return signal;
  }

  if (nameStr == "JpsiFromBplus") {
    MCProng prong(2, {443, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Jpsi from B+ decays", {prong}, {1});
    return signal;
  }

  if (nameStr == "eFromJpsiFromBplus") {
    MCProng prong(3, {11, 443, 521}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Electrons from Jpsi from B+ decays", {prong}, {1});
    return signal;
  }

  if (nameStr == "electronFromJpsiFromBplus") {
    MCProng prong(3, {11, 443, 521}, {false, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Electrons from Jpsi from B+ decays", {prong}, {1});
    return signal;
  }

  if (nameStr == "positronFromJpsiFromBplus") {
    MCProng prong(3, {-11, 443, 521}, {false, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Positrons from Jpsi from B+ decays", {prong}, {1});
    return signal;
  }

  if (nameStr == "eeFromJpsiFromBplus") {
    MCProng prong(3, {11, 443, 521}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Electron pair from Jpsi from B+ decays", {prong, prong}, {1, 1});
    return signal;
  }

  if (nameStr == "eeKaonFromBplus") {
    MCProng pronge(3, {11, 443, 521}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongKaon(2, {321, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon and electron pair from B+", {pronge, pronge, prongKaon}, {2, 2, 1});
    return signal;
  }

  if (nameStr == "eeFromJpsiKaonAny") {
    MCProng pronge(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongKaon(1, {321}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Kaon and electron pair", {pronge, pronge, prongKaon}, {-1, -1, -1});
    return signal;
  }

  if (nameStr == "eeKaonFromBplusExclusive") {
    MCProng pronge(3, {11, 443, 521}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongKaon(2, {321, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon and electron pair from B+", {pronge, pronge, prongKaon}, {2, 2, 1});
    signal->SetDecayChannelIsExclusive(2, true);
    return signal;
  }

  if (nameStr == "eeKaonFromBplusNotExclusive") {
    MCProng pronge(3, {11, 443, 521}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongKaon(2, {321, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon and electron pair from B+", {pronge, pronge, prongKaon}, {2, 2, 1});
    signal->SetDecayChannelIsNotExclusive(2, true);
    return signal;
  }

  // correlated background
  if (nameStr == "eePionFromBplus") {
    MCProng pronge(3, {11, 443, 521}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPion(2, {211, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Pion and electron pair from B+", {pronge, pronge, prongPion}, {2, 2, 1});
    return signal;
  }

  if (nameStr == "eeKaonFromBplusViaEverything") {
    MCProng pronge(3, {11, 443, 521}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongKaon(3, {321, 0, 521}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Kaon and electron pair from B+ via everything", {pronge, pronge, prongKaon}, {2, 2, 2});
    return signal;
  }

  if (nameStr == "eeKaonFromB0") {
    MCProng pronge(3, {11, 443, 511}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongKaon(2, {321, 511}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon and electron pair from B0", {pronge, pronge, prongKaon}, {2, 2, 1});
    return signal;
  }

  if (nameStr == "eePionFromB0ViaEverything") { // catching feed-down for B0
    MCProng pronge(3, {11, 443, 511}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPion(3, {211, 0, 511}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Pion and electron pair from B0", {pronge, pronge, prongPion}, {2, 2, 1});
    return signal;
  }

  if (nameStr == "eeKaonFromOpenBeautyMesons") {
    MCProng pronge(3, {11, 443, 501}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongKaon(2, {321, 501}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Excited kaon and electron pair from B0", {pronge, pronge, prongKaon}, {2, 2, 2});
    return signal;
  }

  if (nameStr == "eeKaonFromOpenBeautyHadrons") {
    MCProng pronge(3, {11, 443, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongKaon(2, {321, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon and electron pair from open beauty hadrons", {pronge, pronge, prongKaon}, {2, 2, 1});
    return signal;
  }

  if (nameStr == "eeKaonFromLambdaB") {
    MCProng pronge(3, {11, 443, 5122}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongKaon(2, {321, 5122}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon and electron pair from lambda B", {pronge, pronge, prongKaon}, {2, 2, 1});
    return signal;
  }

  if (nameStr == "eeKaonPion0FromBplus") {
    MCProng pronge(3, {11, 443, 521}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongKaon(2, {321, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongPion(2, {111, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon, pi0 and electron pair from B+", {pronge, pronge, prongKaon, prongPion}, {2, 2, 1, 1});
    return signal;
  }

  if (nameStr == "eeKaonEtaFromBplus") {
    MCProng pronge(3, {11, 443, 521}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongKaon(2, {321, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongEta(2, {221, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon, eta and electron pair from B+", {pronge, pronge, prongKaon, prongEta}, {2, 2, 1, 1});
    return signal;
  }

  if (nameStr == "eeKaonOmegaFromBplus") {
    MCProng pronge(3, {11, 443, 521}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongKaon(2, {321, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongOmega(2, {223, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon, omega and electron pair from B+", {pronge, pronge, prongKaon, prongOmega}, {2, 2, 1, 1});
    return signal;
  }

  if (nameStr == "eeKaonPionFromBplus") {
    MCProng pronge(3, {11, 443, 521}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongKaon(2, {321, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongPion(2, {211, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon, pion and electron pair from B+", {pronge, pronge, prongKaon, prongPion}, {2, 2, 1, 1});
    return signal;
  }

  if (nameStr == "Bplus") {
    MCProng prong(1, {521}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "B+", {prong}, {-1});
    return signal;
  }

  if (nameStr == "BplusFS") {
    MCProng prong(1, {521}, {true}, {false}, {0}, {0}, {false});
    prong.SetSourceBit(0, MCProng::kHEPMCFinalState);
    signal = new MCSignal(name, "B+", {prong}, {-1});
    return signal;
  }

  if (nameStr == "beautyPairs") {
    MCProng prong(1, {503}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal("beautyPairs", "Beauty hadron pair", {prong, prong}, {-1, -1});
    return signal;
  }
  if (nameStr == "everythingFromBeautyPairs") {
    MCProng prong(2, {0, 503}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal("everythingFromBeautyPairs", "Everything from beauty hadrons pair", {prong, prong}, {-1, -1});
    return signal;
  }
  if (nameStr == "everythingFromEverythingFromBeautyPairsCM") {
    MCProng prong(3, {0, 0, 503}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal("everythingFromEverythingFromBeautyPairs", "Everything from everything from beauty hadrons pair with common grand-mother", {prong, prong}, {2, 2});
    return signal;
  }
  if (nameStr == "everythingFromBeautyANDeverythingFromEverythingFromBeautyPairs") {
    MCProng prong1(3, {0, 0, 503}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prong2(2, {0, 503}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal("everythingFromBeautyANDeverythingFromEverythingFromBeautyPairs", "Everything beauty and everything from everything from beauty hadrons pair", {prong1, prong2}, {2, 1});
    return signal;
  }

  //------------------------------------------------------------------------------------

  if (nameStr == "D0") {
    MCProng prong(1, {Pdg::kD0}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "D0", {prong}, {-1});
    return signal;
  }
  if (nameStr == "nonPromptD0") {
    MCProng prong(2, {Pdg::kD0, 503}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Non-prompt D0", {prong}, {-1});
    return signal;
  }
  if (nameStr == "D0FS") {
    MCProng prong(1, {Pdg::kD0}, {true}, {false}, {0}, {0}, {false});
    prong.SetSourceBit(0, MCProng::kHEPMCFinalState);
    signal = new MCSignal(name, "D0", {prong}, {-1});
    return signal;
  }
  if (nameStr == "KPiFromD0") {
    MCProng prongKaon(2, {321, Pdg::kD0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongPion(2, {211, Pdg::kD0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon and pion pair from D0", {prongKaon, prongPion}, {1, 1});
    return signal;
  }
  if (nameStr == "KPiFromD0Reflected") {
    MCProng prongFalseKaon(2, {211, Pdg::kD0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongFalsePion(2, {321, Pdg::kD0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon and pion pair from D0 with reflected mass assumption", {prongFalseKaon, prongFalsePion}, {1, 1});
    return signal;
  }
  if (nameStr == "Dcharged") {
    MCProng prong(1, {Pdg::kDPlus}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "D+/-", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Dplus") {
    MCProng prong(1, {Pdg::kDPlus}, {false}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "D+", {prong}, {-1});
    return signal;
  }
  if (nameStr == "Dminus") {
    MCProng prong(1, {-Pdg::kDPlus}, {false}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "D+", {prong}, {-1});
    return signal;
  }
  if (nameStr == "KPiPiFromDcharged") {
    MCProng prongKaon(2, {321, Pdg::kDPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongPion(2, {211, Pdg::kDPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon pion pion triplet from D+/-", {prongKaon, prongPion, prongPion}, {1, 1, 1});
    return signal;
  }
  if (nameStr == "KPiPiFromDplus") {
    MCProng prongKaon(2, {-321, Pdg::kDPlus}, {false, false}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongPion(2, {211, Pdg::kDPlus}, {false, false}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon pion pion triplet from D+", {prongKaon, prongPion, prongPion}, {1, 1, 1});
    return signal;
  }
  if (nameStr == "KPiPiFromDminus") {
    MCProng prongKaon(2, {321, -Pdg::kDPlus}, {false, false}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongPion(2, {-211, -Pdg::kDPlus}, {false, false}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon pion pion triplet from D-", {prongKaon, prongPion, prongPion}, {1, 1, 1});
    return signal;
  }
  if (nameStr == "Dstar") {
    MCProng prong(1, {Pdg::kDStar}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "D*", {prong}, {-1});
    return signal;
  }
  if (nameStr == "DstarPlus") {
    MCProng prong(1, {Pdg::kDStar}, {false}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "D*+", {prong}, {-1});
    return signal;
  }
  if (nameStr == "DstarMinus") {
    MCProng prong(1, {-Pdg::kDStar}, {false}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "D*-", {prong}, {-1});
    return signal;
  }
  if (nameStr == "pionFromDstar") {
    MCProng prong(2, {211, Pdg::kDStar}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Pions from D* decays", {prong}, {1});
    return signal;
  }
  if (nameStr == "D0FromDstar") {
    MCProng prong(2, {Pdg::kD0, Pdg::kDStar}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "D0 from D* decays", {prong}, {1});
    return signal;
  }
  if (nameStr == "KFromD0FromDstar") {
    MCProng prong(3, {321, Pdg::kD0, Pdg::kDStar}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Kaons from D0 from D* decays", {prong}, {1});
    return signal;
  }
  if (nameStr == "PiFromD0FromDstar") {
    MCProng prong(3, {211, Pdg::kD0, Pdg::kDStar}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Pions from D0 from D* decays", {prong}, {1});
    return signal;
  }
  if (nameStr == "KPiFromD0FromDstar") {
    MCProng prongKaon(3, {321, Pdg::kD0, Pdg::kDStar}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPion(3, {321, Pdg::kD0, Pdg::kDStar}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Kaon and pion pair from D0 from D* decay", {prongKaon, prongPion}, {1, 1});
    return signal;
  }
  if (nameStr == "KPiPiFromD0FromDstar") {
    MCProng prongKaon(3, {321, Pdg::kD0, Pdg::kDStar}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPionSecondary(3, {211, Pdg::kD0, Pdg::kDStar}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPion(2, {211, Pdg::kDStar}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon pion pion triplet from D*", {prongKaon, prongPionSecondary, prongPion}, {2, 2, 1});
    return signal;
  }
  if (nameStr == "KPiPiFromD0FromDstarPlus") {
    MCProng prongKaon(3, {-321, Pdg::kD0, Pdg::kDStar}, {false, false, false}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPionSecondary(3, {211, Pdg::kD0, Pdg::kDStar}, {false, false, false}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPion(2, {211, Pdg::kDStar}, {false, false}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon pion pion triplet from D*+", {prongKaon, prongPionSecondary, prongPion}, {2, 2, 1});
    return signal;
  }
  if (nameStr == "KPiPiFromD0FromDstarMinus") {
    MCProng prongKaon(3, {321, Pdg::kD0, Pdg::kDStar}, {false, false, false}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPionSecondary(3, {-211, Pdg::kD0, Pdg::kDStar}, {false, false, false}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPion(2, {-211, Pdg::kDStar}, {false, false}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon pion pion triplet from D*-", {prongKaon, prongPionSecondary, prongPion}, {2, 2, 1});
    return signal;
  }
  if (nameStr == "KFromDplus") {
    MCProng prong(2, {321, Pdg::kDPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Kaons from D+/- decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "LambdaC") {
    MCProng prong(1, {Pdg::kLambdaCPlus}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Lambda_c", {prong}, {-1});
    return signal;
  }

  //--------------------------------------------------------------------------------

  if (nameStr == "JpsiFromChic0") {
    MCProng prong(2, {443, 10441}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Jpsi from Chic0 decays", {prong}, {1});
    return signal;
  }
  if (nameStr == "eFromJpsiFromChic0") {
    MCProng prong(3, {11, 443, 10441}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from Jpsi from Chic0 decays", {prong}, {1});
    return signal;
  }
  if (nameStr == "eeFromJpsiFromChic0") {
    MCProng prong(3, {11, 443, 10441}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electron pair from Jpsi from Chic0 decays", {prong, prong}, {2, 2});
    return signal;
  }
  if (nameStr == "JpsiFromChic1") {
    MCProng prong(2, {443, 20443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Jpsi from Chic1 decays", {prong}, {1});
    return signal;
  }
  if (nameStr == "eFromJpsiFromChic1") {
    MCProng prong(3, {11, 443, 20443}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from Jpsi from Chic1 decays", {prong}, {2});
    return signal;
  }
  if (nameStr == "eeFromJpsiFromChic1") {
    MCProng prong(3, {11, 443, 20443}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electron pair from Jpsi from Chic1 decays", {prong, prong}, {2, 2});
    return signal;
  }
  if (nameStr == "JpsiFromChic2") {
    MCProng prong(2, {443, 445}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Jpsi from Chic2 decays", {prong}, {1});
    return signal;
  }
  if (nameStr == "JpsiFromChic2") {
    MCProng prong(2, {443, 904}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Jpsi from Chic0, Chic1 or Chic2 decays", {prong}, {1});
    return signal;
  }
  if (nameStr == "eFromJpsiFromChic2") {
    MCProng prong(3, {11, 443, 445}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from Jpsi from Chic2 decays", {prong}, {2});
    return signal;
  }
  if (nameStr == "eeFromJpsiFromChic2") {
    MCProng prong(3, {11, 443, 445}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electron pair from Jpsi from Chic2 decays", {prong, prong}, {2, 2});
    return signal;
  }
  if (nameStr == "eeFromJpsiFromChic012") {
    MCProng prong(3, {11, 443, 904}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electron pair from Jpsi from Chic0, Chic1 or Chic2 decays", {prong, prong}, {2, 2});
    return signal;
  }
  if (nameStr == "PhotonFromChic0") {
    MCProng prong(2, {22, 10441}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon from Chic0 decays", {prong}, {1});
    return signal;
  }
  if (nameStr == "PhotonFromChic1") {
    MCProng prong(2, {22, 20443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon from Chic1 decays", {prong}, {1});
    return signal;
  }
  if (nameStr == "PhotonFromChic2") {
    MCProng prong(2, {22, 445}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon from Chic2 decays", {prong}, {1});
    return signal;
  }
  if (nameStr == "PhotonFromChic012") {
    MCProng prong(2, {22, 904}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon from Chic0, Chic1, and Chic2 decays", {prong}, {-1});
    return signal;
  }
  if (nameStr == "PhotonFromPi0") {
    MCProng prong(2, {22, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon from Pi0 decays", {prong}, {1});
    return signal;
  }
  if (nameStr == "PhotonPhotonFromPi0") {
    MCProng prong(2, {22, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon Photon from Pi0 decays", {prong, prong}, {1, 1});
    return signal;
  }

  if (nameStr == "eePhotonFromChic1") {
    MCProng pronge(3, {11, 443, 20443}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPhoton(2, {22, 20443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    pronge.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongPhoton.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon and electron pair from Chic1", {pronge, pronge, prongPhoton}, {2, 2, 1});
    return signal;
  }

  if (nameStr == "eePhotonFromChic2") {
    MCProng pronge(3, {11, 443, 445}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPhoton(2, {22, 445}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    pronge.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongPhoton.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon and electron pair from Chic2", {pronge, pronge, prongPhoton}, {2, 2, 1});
    return signal;
  }

  if (nameStr == "eePhotonFromChic12") {
    MCProng pronge(3, {11, 443, MCProng::kPDGCodeNotAssigned}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPhoton(2, {22, MCProng::kPDGCodeNotAssigned}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    pronge.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongPhoton.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon and electron pair from Chic0, Chic1 and Chic2", {pronge, pronge, prongPhoton}, {2, 2, 1});
    return signal;
  }

  if (nameStr == "eePhotonFromPi0") {
    MCProng pronge(2, {11, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongPhoton(2, {22, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    pronge.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongPhoton.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon and electron pair from Pi0", {pronge, pronge, prongPhoton}, {1, 1, 1});
    return signal;
  }

  //--------------------------------------------------------------------------------

  if (nameStr == "X3872") {
    MCProng prong(1, {9920443}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive X(3872)", {prong}, {-1});
    return signal;
  }

  if (nameStr == "JpsiFromX3872") {
    MCProng prong(1, {443}, {true}, {false}, {0}, {0}, {false}, false, {9920443}, {false});
    signal = new MCSignal(name, "Jpsi from X3872", {prong}, {-1});
    return signal;
  }

  if (nameStr == "eFromX3872") {
    MCProng prong(3, {11, 443, 9920443}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Electron from Jpsi from X3872", {prong}, {1});
    return signal;
  }

  if (nameStr == "PionFromX3872") {
    MCProng prong(1, {211}, {true}, {false}, {0}, {0}, {false}, false, {9920443}, {false});
    signal = new MCSignal(name, "Pion from Jpsi from X3872", {prong}, {-1});
    return signal;
  }

  if (nameStr == "JpsiFromPsi2S") {
    MCProng prong(1, {443}, {true}, {false}, {0}, {0}, {false}, false, {100443}, {false});
    signal = new MCSignal(name, "Jpsi from Psi2S", {prong}, {-1});
    return signal;
  }

  if (nameStr == "JpsiFromPromptPsi2S") {
    MCProng prong(2, {443, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {true});
    signal = new MCSignal(name, "Jpsi from prompt Psi2S", {prong}, {-1});
    return signal;
  }

  if (nameStr == "JpsiFromNonpromptPsi2S") {
    MCProng prong(2, {443, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {false});
    signal = new MCSignal(name, "Jpsi from non-prompt Psi2S", {prong}, {-1});
    return signal;
  }

  if (nameStr == "eFromJpsiFromPsi2S") {
    MCProng prong(3, {11, 443, 100443}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Electron from Jpsi from Psi2S", {prong}, {1});
    return signal;
  }

  if (nameStr == "eFromJpsiFromPromptPsi2S") {
    MCProng prong(3, {11, 443, 100443}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false, {503}, {true});
    signal = new MCSignal(name, "Electron from Jpsi from prompt Psi2S", {prong}, {1});
    return signal;
  }

  if (nameStr == "eFromJpsiFromNonpromptPsi2S") {
    MCProng prong(3, {11, 443, 100443}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false, {503}, {false});
    signal = new MCSignal(name, "Electron from Jpsi from non-prompt Psi2S", {prong}, {1});
    return signal;
  }

  if (nameStr == "PionFromPsi2S") {
    MCProng prong(1, {211}, {true}, {false}, {0}, {0}, {false}, false, {100443}, {false});
    signal = new MCSignal(name, "Pion from Jpsi from Psi2S", {prong}, {-1});
    return signal;
  }

  if (nameStr == "PionFromPromptPsi2S") {
    MCProng prong(2, {211, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {true});
    signal = new MCSignal(name, "Pion from prompt Psi2S", {prong}, {-1});
    return signal;
  }

  if (nameStr == "PionFromNonpromptPsi2S") {
    MCProng prong(2, {211, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {false});
    signal = new MCSignal(name, "Pion from non-prompt Psi2S", {prong}, {-1});
    return signal;
  }

  if (nameStr == "eeFromJpsiFromX3872") {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {9920443}, {false});
    signal = new MCSignal(name, "Electron pair from Jpsi from X3872", {prong, prong}, {1, 1});
    return signal;
  }

  if (nameStr == "JpsiPiPiFromX3872") {
    MCProng prongJpsi(2, {443, 9920443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongPi(2, {211, 9920443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Jpsi and pion pair from X3872", {prongJpsi, prongPi, prongPi}, {1, 1, 1});
    return signal;
  }

  if (nameStr == "eePiPiFromX3872") {
    MCProng pronge(3, {11, 443, 9920443}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPi(2, {211, 9920443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electron pair and pion pair from X3872", {pronge, pronge, prongPi, prongPi}, {2, 2, 1, 1});
    return signal;
  }

  if (nameStr == "eeFromJpsiFromPsi2S") {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {100443}, {false});
    signal = new MCSignal(name, "Electron pair from Jpsi from Psi2S", {prong, prong}, {1, 1});
    return signal;
  }

  if (nameStr == "eeFromJpsiFromPromptPsi2S") {
    MCProng prong(3, {11, 443, 100443}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false, {503}, {true});
    signal = new MCSignal(name, "Electron pair from Jpsi from prompt Psi2S", {prong, prong}, {1, 1});
    return signal;
  }

  if (nameStr == "eeFromJpsiFromNonpromptPsi2S") {
    MCProng prong(3, {11, 443, 100443}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false, {503}, {false});
    signal = new MCSignal(name, "Electron pair from Jpsi from non-prompt Psi2S", {prong, prong}, {1, 1});
    return signal;
  }

  if (nameStr == "JpsiPiPiFromPsi2S") {
    MCProng prongJpsi(2, {443, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongPi(2, {211, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Jpsi and pion pair from Psi2S", {prongJpsi, prongPi, prongPi}, {1, 1, 1});
    return signal;
  }

  if (nameStr == "eePiPiFromPsi2S") {
    MCProng pronge(3, {11, 443, 100443}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPi(2, {211, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electron pair and pion pair from Psi2S", {pronge, pronge, prongPi, prongPi}, {2, 2, 1, 1});
    return signal;
  }

  if (nameStr == "eePiPiFromPromptPsi2S") {
    MCProng pronge(3, {11, 443, 100443}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false, {503}, {true});
    MCProng prongPi(2, {211, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {true});
    signal = new MCSignal(name, "Electron pair and pion pair from prompt Psi2S", {pronge, pronge, prongPi, prongPi}, {2, 2, 1, 1});
    return signal;
  }

  if (nameStr == "eePiPiFromNonpromptPsi2S") {
    MCProng pronge(3, {11, 443, 100443}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false, {503}, {false});
    MCProng prongPi(2, {211, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {false});
    signal = new MCSignal(name, "Electron pair and pion pair from non-prompt Psi2S", {pronge, pronge, prongPi, prongPi}, {2, 2, 1, 1});
    return signal;
  }

  if (nameStr == "eeFromPromptJpsiAnyPrimary") {
    MCProng pronge(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {true});
    pronge.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongPrimary(1);
    prongPrimary.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "anyprimary and electron pair from prompt jpsi", {pronge, pronge, prongPrimary}, {1, 1, -1});
    return signal;
  }

  if (nameStr == "eeFromJpsiAnyPrimary") {
    MCProng pronge(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    pronge.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongPrimary(1);
    prongPrimary.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "anyprimary and electron pair from prompt jpsi", {pronge, pronge, prongPrimary}, {1, 1, -1});
    return signal;
  }

  if (nameStr == "eeFromNonPromptJpsiAnyPrimary") {
    MCProng pronge(3, {11, 443, 503}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    pronge.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongPrimary(1);
    prongPrimary.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "anyprimary and electron pair from non-prompt jpsi", {pronge, pronge, prongPrimary}, {1, 1, -1});
    return signal;
  }
  return nullptr;
}

//_______________________________________________________________________________________________
std::vector<MCSignal*> o2::aod::dqmcsignals::GetMCSignalsFromJSON(const char* json)
{
  //
  // configure MC signals using a json file
  //
  std::vector<MCSignal*> signals;
  LOG(info) << "========================================== interpreting JSON for MC signals";
  LOG(info) << "JSON string: " << json;
  //
  // Create a vector of MCSignal from a JSON formatted string
  //   The JSON is expected to contain a list of objects, with each object containing the fields needed
  //    to define an MCSignal
  rapidjson::Document document;

  // Check that the json is parsed correctly
  rapidjson::ParseResult ok = document.Parse(json);
  if (!ok) {
    TString str = "";
    for (int i = ok.Offset() - 30; i < static_cast<int>(ok.Offset()) + 50; i++) {
      if ((i >= 0) && (i < static_cast<int>(strlen(json)))) {
        str += json[i];
      }
    }
    LOG(fatal) << "JSON parse error: " << rapidjson::GetParseErrorFunc(ok.Code()) << " (" << ok.Offset() << ")" << " **** Parsing error is somewhere here: " << str.Data();
    return signals;
  }

  // loop over the top level objects in the json
  for (rapidjson::Value::ConstMemberIterator it = document.MemberBegin(); it != document.MemberEnd(); it++) {

    const char* sigName = it->name.GetString();
    LOG(info) << "=================================================== Configuring MC signal " << sigName;
    const auto& signal = it->value;

    // Validate the entry for this MCSignal
    if (!ValidateJSONMCSignal(&signal, sigName)) {
      LOG(fatal) << "MCSignal JSON not properly defined for " << sigName << ". Skipping";
      continue;
    } else {
      LOG(debug) << "MCSignal validated";
    }

    // Get the signal title
    const char* title = (signal.HasMember("title") ? signal.FindMember("title")->value.GetString() : "");
    LOG(info) << "Title is: " << title;

    // Get the exclude common ancestor
    bool excludeCommonAncestor = false;
    if (signal.HasMember("excludeCommonAncestor")) {
      excludeCommonAncestor = signal.FindMember("excludeCommonAncestor")->value.GetBool();
    }
    LOG(debug) << "exclude common ancestor " << excludeCommonAncestor;

    // Check for MCProng objects in the json
    std::vector<MCProng> prongs;
    for (rapidjson::Value::ConstMemberIterator prongIt = signal.MemberBegin(); prongIt != signal.MemberEnd(); prongIt++) {

      // If the name is not MCProng, continue
      TString prongName = prongIt->name.GetString();
      if (!prongName.Contains("MCProng")) {
        continue;
      }

      // Call the function to parse the MCProng object and get the pointer to the created MCProng
      MCProng* prong = ParseJSONMCProng(&prongIt->value, prongName.Data());
      if (prong == nullptr) {
        LOG(fatal) << "MCProng not built! MCSignal not configured";
        return signals;
      }
      LOG(debug) << "MCProng defined";
      // Print the contents of the configured prong
      prong->Print();
      // push the prong to the vector
      prongs.push_back(*prong);
    }

    // Get the common ancestors array
    std::vector<int8_t> commonAncestors;
    if (signal.HasMember("commonAncestors")) {
      for (auto& v : signal.FindMember("commonAncestors")->value.GetArray()) {
        commonAncestors.push_back(v.GetInt());
        LOG(debug) << "common ancestor " << v.GetInt();
      }
    } else {
      for (uint32_t i = 0; i < prongs.size(); i++) {
        commonAncestors.push_back(-1);
      }
    }

    if (prongs.size() == 0) {
      LOG(fatal) << "No prongs were defined for this MCSignal!";
      return signals;
    }

    // Check that we have as many prongs defined as the size of the common ancestors array
    if (prongs.size() != commonAncestors.size()) {
      LOG(fatal) << "Number of defined prongs and size of commonAncestors array must coincide in MCSignal definition";
      return signals;
    }

    // Create the signal and add it to the output vector
    auto* mcSignal = new MCSignal(sigName, title, prongs, commonAncestors, excludeCommonAncestor);
    LOG(debug) << "MCSignal defined, adding to the output vector";
    mcSignal->PrintConfig();
    signals.push_back(mcSignal);
  }

  return signals;
}

//_______________________________________________________________________________________________
template <typename T>
bool o2::aod::dqmcsignals::ValidateJSONMCProng(T prongJSON, const char* prongName)
{

  // Check that the json entry for this prong is correctly given
  LOG(debug) << "Validating the prong " << prongName;

  // The fields for the number of generations, pdg codes and checkBothCharges are required
  if (!prongJSON->HasMember("n") || !prongJSON->HasMember("pdgs") || !prongJSON->HasMember("checkBothCharges")) {
    LOG(fatal) << "Missing either n, pdgs or checkBothCharges fields in MCProng JSON definition";
    return false;
  }
  // the size of the pdgs array must be equal to n
  int n = prongJSON->FindMember("n")->value.GetInt();
  uint32_t nSigned = static_cast<uint64_t>(n);
  if (prongJSON->FindMember("pdgs")->value.GetArray().Size() != nSigned) {
    LOG(fatal) << "Size of the pdgs array must be equal to n in MCProng JSON definition";
    return false;
  }
  // the size of the checkBothCharges array must be equal to n
  if (prongJSON->FindMember("checkBothCharges")->value.GetArray().Size() != nSigned) {
    LOG(fatal) << "Size of the checkBothCharges array must be equal to n in MCProng JSON definition";
    return false;
  }
  // the size of the exclude pdg array must be equal to n
  if (prongJSON->HasMember("excludePDG")) {
    if (prongJSON->FindMember("excludePDG")->value.GetArray().Size() != nSigned) {
      LOG(fatal) << "Size of the excludePDG array must be equal to n in MCProng JSON definition";
      return false;
    }
  }

  // Check the corectness of the source bits fields, if these are specified
  // The sourceBits field should be an array of size n, with each element being another array containing an array of sources specified as strings
  //    The source strings have to be the ones specified in MCProng::Source
  //  If the excludeSource is specified in addition, then this field should be an array of size n, with each element being another array of booleans
  //     corresponding to each source specified in sourceBits
  if (prongJSON->HasMember("sourceBits")) {
    if (prongJSON->FindMember("sourceBits")->value.GetArray().Size() != nSigned) {
      LOG(fatal) << "Size of the sourceBits array must be equal to n in MCProng JSON definition";
      return false;
    }
    std::vector<uint32_t> nSourceBits;
    for (auto& ii : prongJSON->FindMember("sourceBits")->value.GetArray()) {
      if (!ii.IsArray()) {
        LOG(fatal) << "The sourceBits field should be an array of arrays of MCProng::Source";
        return false;
      }
      nSourceBits.push_back(ii.GetArray().Size());
      for (auto& iii : ii.GetArray()) {
        if (MCProng::fgSourceNames.find(iii.GetString()) == MCProng::fgSourceNames.end()) {
          LOG(fatal) << "Source " << iii.GetString() << " not implemented in MCProng";
          return false;
        }
      }
    }
    if (prongJSON->HasMember("excludeSource")) {
      if (prongJSON->FindMember("excludeSource")->value.GetArray().Size() != nSigned) {
        LOG(fatal) << "Size of the excludeSource array must be equal to n in MCProng JSON definition";
        return false;
      }
      int iElem = 0;
      for (auto& ii : prongJSON->FindMember("excludeSource")->value.GetArray()) {
        if (!ii.IsArray()) {
          LOG(fatal) << "The excludeSource field should be an array of arrays of bool";
          return false;
        }
        if (ii.GetArray().Size() != nSourceBits[iElem]) {
          LOG(fatal) << "The size of excludeSource arrays does not match the size of the arrays in sourceBits";
          return false;
        }
        iElem++;
      }
    }
    // Check the useAND on source bit map
    if (prongJSON->HasMember("useANDonSourceBitMap")) {
      if (prongJSON->FindMember("useANDonSourceBitMap")->value.GetArray().Size() != nSigned) {
        LOG(fatal) << "Size of the useANDonSourceBitMap array must be equal to n in MCProng JSON definition";
        return false;
      }
    }
  }

  // sourceBits is needed in case the other source related fields are specified
  if ((prongJSON->HasMember("excludeSource") || prongJSON->HasMember("useANDonSourceBitMap")) && !prongJSON->HasMember("sourceBits")) {
    LOG(fatal) << "Field sourceBits is needed when specifying excludeSource or useANDonSourceBitMap";
    return false;
  }

  // check checkGenerationsInTime
  if (prongJSON->HasMember("checkGenerationsInTime")) {
    if (!prongJSON->FindMember("checkGenerationsInTime")->value.IsBool()) {
      LOG(fatal) << "Field checkGeneretionsInTime must be boolean";
      return false;
    }
  }

  if (prongJSON->HasMember("checkIfPDGInHistory")) {
    if (!prongJSON->FindMember("checkIfPDGInHistory")->value.IsArray()) {
      LOG(fatal) << "Field checkGeneretionsInTime must be an array of integers";
      return false;
    }
    uint32_t vecSize = prongJSON->FindMember("checkIfPDGInHistory")->value.GetArray().Size();
    if (prongJSON->HasMember("excludePDGInHistory")) {
      if (!prongJSON->FindMember("excludePDGInHistory")->value.IsArray()) {
        LOG(fatal) << "Field excludePDGInHistory must be an array of booleans";
        return false;
      }
      if (prongJSON->FindMember("excludePDGInHistory")->value.GetArray().Size() != vecSize) {
        LOG(fatal) << "Field excludePDGInHistory must be an array of equal size with the array specified by checkIfPDGInHistory";
        return false;
      }
    }
  } else {
    if (prongJSON->HasMember("excludePDGInHistory")) {
      LOG(fatal) << "Field checkIfPDGInHistory is required when excludePDGInHistory is specified";
      return false;
    }
  }

  return true;
}

//_______________________________________________________________________________________________
template <typename T>
MCProng* o2::aod::dqmcsignals::ParseJSONMCProng(T prongJSON, const char* prongName)
{

  // Check that the entry for this prong is validated
  LOG(debug) << "Parsing the prong " << prongName;
  if (!ValidateJSONMCProng(prongJSON, prongName)) {
    LOG(fatal) << "MCProng not properly defined in the JSON file.";
    return nullptr;
  }

  // Get the number of generations
  int n = prongJSON->FindMember("n")->value.GetInt();
  LOG(debug) << "n: " << n;
  // Get the array of PDG codes
  std::vector<int> pdgs;
  for (auto& pdg : prongJSON->FindMember("pdgs")->value.GetArray()) {
    pdgs.push_back(pdg.GetInt());
    LOG(debug) << "pdgs: " << pdg.GetInt();
  }
  // get the array of booleans for check both charges option
  std::vector<bool> checkBothCharges;
  for (auto& ii : prongJSON->FindMember("checkBothCharges")->value.GetArray()) {
    checkBothCharges.push_back(ii.GetBool());
    LOG(debug) << "check both charges " << ii.GetBool();
  }

  // get the array of booleans for the excludePDG option, defaults to false
  std::vector<bool> excludePDG;
  if (prongJSON->HasMember("excludePDG")) {
    for (auto& ii : prongJSON->FindMember("excludePDG")->value.GetArray()) {
      excludePDG.push_back(ii.GetBool());
      LOG(debug) << "exclude pdg " << ii.GetBool();
    }
  } else {
    for (int i = 0; i < n; i++) {
      excludePDG.push_back(false);
    }
  }

  // get the source bits, and transform from string to int
  std::vector<std::vector<int>> sourceBitsVec;
  if (prongJSON->HasMember("sourceBits")) {
    for (auto& ii : prongJSON->FindMember("sourceBits")->value.GetArray()) {
      std::vector<int> sourceBits;
      for (auto& iii : ii.GetArray()) {
        sourceBits.push_back(MCProng::fgSourceNames[iii.GetString()]);
        LOG(debug) << "source bit " << iii.GetString();
      }
      sourceBitsVec.push_back(sourceBits);
    }
  }
  // prepare the exclusion source options if specified
  std::vector<std::vector<bool>> excludeSourceVec;
  if (prongJSON->HasMember("excludeSource")) {
    for (auto& ii : prongJSON->FindMember("excludeSource")->value.GetArray()) {
      std::vector<bool> excludeSource;
      for (auto& iii : ii.GetArray()) {
        excludeSource.push_back(iii.GetBool());
        LOG(debug) << "exclude source bit " << iii.GetBool();
      }
      excludeSourceVec.push_back(excludeSource);
    }
  }

  // prepare the useANDonSourceBitMap vector, defaults to true for each generation
  std::vector<bool> useANDonSourceBitMap;
  if (prongJSON->HasMember("useANDonSourceBitMap")) {
    for (auto& ii : prongJSON->FindMember("useANDonSourceBitMap")->value.GetArray()) {
      useANDonSourceBitMap.push_back(ii.GetBool());
      LOG(debug) << "use AND on source map " << ii.GetBool();
    }
  } else {
    for (int i = 0; i < n; i++) {
      useANDonSourceBitMap.push_back(true);
    }
  }

  // prepare the bit maps suitable for the MCProng constructor
  bool hasExclude = prongJSON->HasMember("excludeSource");
  int igen = 0;
  std::vector<uint64_t> sBitsVec;
  std::vector<uint64_t> sBitsExcludeVec;
  for (auto& itgen : sourceBitsVec) {
    int is = 0;
    uint64_t sBits = 0;
    uint64_t sBitsExclude = 0;
    auto excludeVec = (hasExclude ? excludeSourceVec[igen] : std::vector<bool>{});
    for (auto& s : itgen) {
      bool exclude = (hasExclude ? excludeVec[is] : false);
      if (s != MCProng::kNothing) {
        sBits |= (static_cast<uint64_t>(1) << s);
        if (exclude) {
          sBitsExclude |= (static_cast<uint64_t>(1) << s);
        }
      }
      is++;
    }
    sBitsVec.push_back(sBits);
    sBitsExcludeVec.push_back(sBitsExclude);
    LOG(debug) << "igen " << igen;
    LOG(debug) << "igen sBits " << sBits;
    LOG(debug) << "igen exclude " << sBitsExclude;
    igen++;
  }

  // check that the sourceBits has the size of n generations
  if (prongJSON->HasMember("sourceBits")) {
    if (sBitsVec.size() != static_cast<uint32_t>(n)) {
      LOG(fatal) << "sourceBits array should have a size equal to n";
    }
  } else {
    sBitsVec.clear();
    for (int i = 0; i < n; i++) {
      sBitsVec.push_back(0);
    }
  }
  // check that the sourceBits exclude has the size of n generations
  if (prongJSON->HasMember("excludeSource")) {
    if (sBitsExcludeVec.size() != static_cast<uint32_t>(n)) {
      LOG(fatal) << "sourceBits exclude array should have a size equal to n";
    }
  } else {
    sBitsExcludeVec.clear();
    for (int i = 0; i < n; i++) {
      sBitsExcludeVec.push_back(0);
    }
  }

  bool checkGenerationsInTime = false;
  if (prongJSON->HasMember("checkGenerationsInTime")) {
    checkGenerationsInTime = prongJSON->FindMember("checkGenerationsInTime")->value.GetBool();
  }
  LOG(debug) << "checkGenerationsInTime: " << checkGenerationsInTime;

  std::vector<int> checkIfPDGInHistory = {};
  if (prongJSON->HasMember("checkIfPDGInHistory")) {
    for (auto& ii : prongJSON->FindMember("checkIfPDGInHistory")->value.GetArray()) {
      checkIfPDGInHistory.push_back(ii.GetInt());
      LOG(debug) << "checkIfPDGInHistory: " << ii.GetInt();
    }
  }

  std::vector<bool> excludePDGInHistory = {};
  if (prongJSON->HasMember("excludePDGInHistory")) {
    for (auto& ii : prongJSON->FindMember("excludePDGInHistory")->value.GetArray()) {
      excludePDGInHistory.push_back(ii.GetBool());
      LOG(debug) << "excludePDGInHistory: " << ii.GetBool();
    }
  }

  // Calling the MCProng constructor
  auto* prong = new MCProng(n, pdgs, checkBothCharges, excludePDG, sBitsVec, sBitsExcludeVec, useANDonSourceBitMap,
                            checkGenerationsInTime, checkIfPDGInHistory, excludePDGInHistory);
  // Print the configuration
  prong->Print();
  return prong;
}

//_______________________________________________________________________________________________
template <typename T>
bool o2::aod::dqmcsignals::ValidateJSONMCSignal(T sigJSON, const char* sigName)
{

  LOG(info) << "Validating MC signal " << sigName;
  if (sigJSON->HasMember("commonAncestors")) {
    if (!sigJSON->FindMember("commonAncestors")->value.IsArray()) {
      LOG(fatal) << "In MCSignal definition, commonAncestors must be an array";
      return false;
    }
  }
  if (sigJSON->HasMember("excludeCommonAncestor") && !sigJSON->HasMember("commonAncestors")) {
    LOG(fatal) << "In MCSignal definition, commonAncestors field is needed if excludeCommonAncestor is specified";
    return false;
  }

  return true;
}
