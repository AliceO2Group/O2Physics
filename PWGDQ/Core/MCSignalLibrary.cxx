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
#include <TPDGCode.h>
#include "CommonConstants/PhysicsConstants.h"
#include "PWGDQ/Core/MCSignalLibrary.h"

using namespace o2::constants::physics;

MCSignal* o2::aod::dqmcsignals::GetMCSignal(const char* name)
{
  std::string nameStr = name;
  MCSignal* signal;
  // 1-prong signals
  if (!nameStr.compare("alicePrimary")) {
    MCProng prong(1);                                  // 1-generation prong
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);  // set source to be ALICE primary particles
    signal = new MCSignal(1, name, "ALICE primaries"); // define a signal with one prong
    signal->AddProng(prong);                           // add the previously defined prong to the signal
    return signal;
  }
  if (!nameStr.compare("electron")) {
    MCProng prong(1, {11}, {true}, {false}, {0}, {0}, {false});        // define 1-generation prong using the full constructor
    signal = new MCSignal(name, "Inclusive electrons", {prong}, {-1}); // define the signal using the full constructor
    return signal;
  }
  if (!nameStr.compare("electronPrimary")) {
    MCProng prong(1, {11}, {true}, {false}, {0}, {0}, {false});      // define 1-generation prong using the full constructor
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);                // set source to be ALICE primary particles
    signal = new MCSignal(name, "Primary electrons", {prong}, {-1}); // define the signal using the full constructor
    return signal;
  }
  if (!nameStr.compare("pionPrimary")) {
    MCProng prong(1, {211}, {true}, {false}, {0}, {0}, {false}); // define 1-generation prong using the full constructor
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);            // set source to be ALICE primary particles
    signal = new MCSignal(name, "Primary pions", {prong}, {-1}); // define the signal using the full constructor
    return signal;
  }
  if (!nameStr.compare("pionPrimaryFromHc")) {
    MCProng prong(2, {211, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);
    signal = new MCSignal(name, "Primary pions from open charmed hadron decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("pionPrimaryFromHb")) {
    MCProng prong(2, {211, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);
    signal = new MCSignal(name, "Primary pions from open beauty hadron decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("photon")) {
    MCProng prong(1, {22}, {true}, {false}, {0}, {0}, {false}); // define 1-generation prong using the full constructor
    signal = new MCSignal(name, "Photon", {prong}, {-1});       // define the signal using the full constructor
    return signal;
  }
  if (!nameStr.compare("muonPrimary")) {
    MCProng prong(1, {13}, {true}, {false}, {0}, {0}, {false});  // define 1-generation prong using the full constructor
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);            // set source to be ALICE primary particles
    signal = new MCSignal(name, "Primary Muons", {prong}, {-1}); // define the signal using the full constructor
    return signal;
  }
  if (!nameStr.compare("kaonFromPhi")) {
    MCProng prong(2, {321, 333}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}); // define 2-generation prong using the full constructor
    signal = new MCSignal(name, "Kaons from phi-mesons", {prong}, {-1});                        // define the signal using the full constructor
    return signal;
  }
  if (!nameStr.compare("kaonPrimary")) {
    MCProng prong(1, {321}, {true}, {false}, {0}, {0}, {false}); // define 1-generation prong using the full constructor
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);            // set source to be ALICE primary particles
    signal = new MCSignal(name, "Primary Kaons", {prong}, {-1}); // define the signal using the full constructor
    return signal;
  }
  if (!nameStr.compare("protonPrimary")) {
    MCProng prong(1, {2212}, {true}, {false}, {0}, {0}, {false}); // define 1-generation prong using the full constructor
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);             // set source to be ALICE primary particles
    signal = new MCSignal(name, "Primary Proton", {prong}, {-1}); // define the signal using the full constructor
    return signal;
  }
  if (!nameStr.compare("phiMeson")) {
    MCProng prong(1, {333}, {true}, {false}, {0}, {0}, {false}); // define 1-generation prong using the full constructor
    signal = new MCSignal(name, "phi meson", {prong}, {-1});     // define the signal using the full constructor
    return signal;
  }
  if (!nameStr.compare("muon")) {
    MCProng prong(1, {13}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive muons", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("electronNOTfromTransport")) {
    MCProng prong(1);
    prong.SetPDGcode(0, 11, true);
    prong.SetSourceBit(0, MCProng::kProducedInTransport, true); // exclude particles produces in transport
    signal = new MCSignal(name, "Electrons which are not produced during transport in detector", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromNonpromptJpsi")) {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {false});
    signal = new MCSignal(name, "Electrons from non-prompt jpsi decays with beauty in decay chain", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("ePrimaryFromPromptJpsi")) {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from prompt jpsi decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("Jpsi")) {
    MCProng prong(1, {443}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive jpsi", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("Helium3")) {
    MCProng prong(1, {1000020030}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Helium3", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("nonPromptJpsi")) {
    MCProng prong(2, {443, 503}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Non-prompt jpsi", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("promptJpsi")) {
    MCProng prong(1, {443}, {true}, {false}, {0}, {0}, {false}, false, {503}, {true});
    signal = new MCSignal(name, "Prompt jpsi (not from beauty)", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromNonpromptPsi2S")) {
    MCProng prong(3, {11, 100443, 503}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Electrons from beauty psi2s decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromPromptPsi2S")) {
    MCProng prong(3, {11, 100443, 503}, {true, true, true}, {false, false, true}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Electrons from prompt psi2s decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("Psi2S")) {
    MCProng prong(1, {100443}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive psi2s", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("nonPromptPsi2S")) {
    MCProng prong(2, {100443, 503}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Non-prompt psi2s", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("promptPsi2S")) {
    MCProng prong(2, {100443, 503}, {true, true}, {false, true}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Prompt psi2s (not from beauty)", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("Chic0")) {
    MCProng prong(1, {10441}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive Chic0", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("Chic1")) {
    MCProng prong(1, {20443}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive Chic1", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("Chic2")) {
    MCProng prong(1, {445}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive Chic2", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("Chic012")) {
    MCProng prong(1, {904}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive Chic0, Chic1 and Chic2", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("Upsilon1S")) {
    MCProng prong(1, {553}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive Upsilon1S", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("Upsilon2S")) {
    MCProng prong(1, {100553}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive Upsilon2S", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("Upsilon3S")) {
    MCProng prong(1, {200553}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive Upsilon3S", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("allBeautyHadrons")) {
    MCProng prong(1, {503}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "All beauty hadrons", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("allOpenBeautyHadrons")) {
    MCProng prong(1, {502}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "All open beauty hadrons", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("Bc")) {
    MCProng prong(1, {541}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Bc", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("mumuFromJpsiFromBc")) {
    MCProng prong(3, {13, 443, 541}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Muon pair from jpsi from Bc decays", {prong, prong}, {1, 1});
    return signal;
  }
  if (!nameStr.compare("muFromBc")) {
    MCProng prong(2, {13, 541}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Muon from Bc decays", {prong}, {1});
    return signal;
  }
  if (!nameStr.compare("mumumuFromBc")) {
    MCProng prongMuFromJpsi(3, {13, 443, 541}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongMuFromBc(2, {13, 541}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Trimuon from Bc decays", {prongMuFromJpsi, prongMuFromJpsi, prongMuFromBc}, {2, 2, 1});
    return signal;
  }
  if (!nameStr.compare("everythingFromBeauty")) {
    MCProng prong(2, {0, 503}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Everything from beauty", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("everythingFromEverythingFromBeauty")) {
    MCProng prong(3, {0, 0, 503}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Everything from everything from beauty", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("allCharmHadrons")) {
    MCProng prong(1, {403}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "All charm hadrons", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("allOpenCharmHadrons")) {
    MCProng prong(1, {402}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "All open charm hadrons", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("allCharmFromBeauty")) {
    MCProng prong(2, {403, 503}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "All charm hadrons from beauty", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("allPromptCharm")) {
    MCProng prong(2, {403, 503}, {true, true}, {false, true}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "All prompt charm hadrons (not from beauty)", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("Pi0DecayToe")) {
    MCProng prong(2, {111, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Pi0 decays into an electron", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("Pi0DecayTog")) {
    MCProng prong(2, {111, 22}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Pi0 decays into an gamma", {prong}, {1});
    return signal;
  }
  if (!nameStr.compare("Pi0DecayTogg")) {
    MCProng prong(2, {111, 22}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Pi0 decays into an gamma, gamma", {prong, prong}, {1, 1});
    return signal;
  }
  if (!nameStr.compare("PromptPi0DecayToe")) {
    MCProng prong(2, {111, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, true, {403, 503}, {true, true});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Pi0 decays into an electron", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("Pi0")) {
    MCProng prong(1, {111}, {true}, {false}, {0}, {0}, {false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
    signal = new MCSignal(name, "Pi0", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("LMeeLFQ")) {
    MCProng prong(1, {900}, {true}, {false}, {0}, {0}, {false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
    signal = new MCSignal(name, "light flavor mesons + quarkonia", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi,jpsi,psi2s
    return signal;
  }
  if (!nameStr.compare("LMeeLF")) {
    MCProng prong(1, {901}, {true}, {false}, {0}, {0}, {false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
    signal = new MCSignal(name, "ligh flavor mesons", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi
    return signal;
  }
  if (!nameStr.compare("PromptJpsiDecayToe")) {
    MCProng prong(2, {443, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, true, {503}, {true});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Prompt jpsi (not from beauty) decay to electron", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("electronFromDs")) {
    MCProng prong(2, {11, 431}, {true, true}, {false, false}, {0, 0}, {0, 0}, {true, true});
    signal = new MCSignal(name, "Electrons from Ds decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("dsMeson")) {
    MCProng prong(1, {431}, {true}, {false}, {0}, {0}, {true});
    signal = new MCSignal(name, "Ds mesons", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("electronFromPC")) {
    MCProng prong(2, {11, 22}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "electron from a photon conversion", {prong}, {-1});
    return signal;
  }

  // 2-prong signals
  if (!nameStr.compare("dielectron")) {
    MCProng prong(1, {11}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Electron pair", {prong, prong}, {-1, -1});
    return signal;
  }
  if (!nameStr.compare("dimuon")) {
    MCProng prong(1, {13}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Muon pair", {prong, prong}, {-1, -1});
    return signal;
  }
  if (!nameStr.compare("electronMuonPair")) {
    MCProng electron(1, {11}, {true}, {false}, {0}, {0}, {false});
    MCProng muon(1, {13}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Electron-muon pair", {electron, muon}, {-1, -1});
    return signal;
  }
  if (!nameStr.compare("dielectronFromPC")) {
    MCProng prong(2, {11, 22}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "dielectron from a photon conversion", {prong, prong}, {1, 1});
    return signal;
  }
  if (!nameStr.compare("dielectronFromAllPC")) {
    MCProng prong(2, {11, 22}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "dielectron from a photon conversion", {prong, prong}, {-1, -1});
    return signal;
  }
  if (!nameStr.compare("dielectronPCPi0")) {
    MCProng prong(3, {11, 22, 111}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "dielectron from a photon conversion from a pi0", {prong, prong}, {1, 1});
    return signal;
  }

  // LMEE single signals
  // electron signals with mother X: e from mother X
  if (!nameStr.compare("eFromAnything")) {
    MCProng prong(2, {11, MCProng::kPDGCodeNotAssigned}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from any mother", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromPhoton")) {
    MCProng prong(2, {11, 22}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from photon conversion", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromPi0")) {
    MCProng prong(2, {11, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from pi0 decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("ePrimaryFromPromptPi0")) {
    MCProng prong(2, {11, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502, 402}, {true, true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from prompt pi0 decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromEta")) {
    MCProng prong(2, {11, 221}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from eta decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromEtaPrime")) {
    MCProng prong(2, {11, 331}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from eta' decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromRho")) {
    MCProng prong(2, {11, 113}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from rho decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromOmega")) {
    MCProng prong(2, {11, 223}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from omega decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromPhi")) {
    MCProng prong(2, {11, 333}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from phi decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromJpsi")) {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from jpsi decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromPromptJpsi")) {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from jpsi decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromPsi2S")) {
    MCProng prong(2, {11, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from psi2s decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromLMeeLF")) {
    MCProng prong(2, {11, 901}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from LF meson decays", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi mesons
    return signal;
  }
  if (!nameStr.compare("ePrimaryFromLMeeLF")) {
    MCProng prong(2, {11, 901}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);                      // set source to be ALICE primary particles
    signal = new MCSignal(name, "Electrons from LF meson decays", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi mesons
    return signal;
  }
  if (!nameStr.compare("eFromLMeeLFQ")) {
    MCProng prong(2, {11, 900}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(false);                                                             // set direction to check generation in time (true) or back in time (false)
    signal = new MCSignal(name, "Electrons from LF meson + quarkonia decays", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi,jpsi,psi2s mesons
    return signal;
  }
  if (!nameStr.compare("ePrimaryFromLMeeLFQ")) {
    MCProng prong(2, {11, 900}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);                                  // set source to be ALICE primary particles
    prong.SetSignalInTime(false);                                                             // set direction to check generation in time (true) or back in time (false)
    signal = new MCSignal(name, "Electrons from LF meson + quarkonia decays", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi,jpsi,psi2s mesons
    return signal;
  }
  if (!nameStr.compare("eFromHc")) {
    MCProng prong(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open charmed hadron decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromAnyHc")) {
    MCProng prong(1, {11}, {true}, {false}, {0}, {0}, {false}, false, {402}, {false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);
    signal = new MCSignal(name, "Electrons from any open charm hadron decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromD0")) {
    MCProng prong(2, {11, Pdg::kD0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from D0 decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromChargedD")) {
    MCProng prong(2, {11, Pdg::kDPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from D+/- decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromDs")) {
    MCProng prong(2, {11, Pdg::kDS}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from Ds +/- decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromLambdaC")) {
    MCProng prong(2, {11, Pdg::kLambdaCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from Lambda_c decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromXiC0")) {
    MCProng prong(2, {11, Pdg::kXiC0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from Xi_c_0 decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromXiCPlus")) {
    MCProng prong(2, {11, Pdg::kXiCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from Xi_c_+ decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromXiCPlusPlus")) {
    MCProng prong(2, {11, Pdg::kXiCCPlusPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from Xi_c_++ decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromHb")) {
    MCProng prong(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open beauty hadron decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromAnyHb")) {
    MCProng prong(1, {11}, {true}, {false}, {0}, {0}, {false}, false, {502}, {false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);
    signal = new MCSignal(name, "Electrons from any open beauty hadron decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromHbc")) {
    MCProng prong(2, {11, 902}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open charm or beauty hadron decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromMc")) {
    MCProng prong(2, {11, 401}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open charmed meson decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromMb")) {
    MCProng prong(2, {11, 501}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open beauty meson decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromBc")) {
    MCProng prong(2, {11, 4001}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open charmed baryon decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromBb")) {
    MCProng prong(2, {11, 5001}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open beauty baryon decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromPromptHc")) {
    MCProng prong(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open charmed hadron decays without beauty in decay history", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromHbtoHc")) {
    MCProng prong(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from open charmed hadron decays with b hadron in decay history", {prong}, {-1});
    return signal;
  }

  // muon signals with mother X: mu from mother X
  if (!nameStr.compare("muFromJpsi")) {
    MCProng prong(2, {13, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "muons from jpsi decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("muFromPsi2S")) {
    MCProng prong(2, {13, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "muons from psi2s decays", {prong}, {-1});
    return signal;
  }

  // Decay signal: Mother to electron: X -> e
  if (!nameStr.compare("AnythingToE")) {
    MCProng prong(2, {MCProng::kPDGCodeNotAssigned, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
    prong.SetSignalInTime(true);                                            // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Decay of anything into e", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi,jpsi,psi2s mesons
    return signal;
  }
  if (!nameStr.compare("LFQdecayToE")) {
    MCProng prong(2, {900, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
    prong.SetSignalInTime(true);                                                       // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "LF meson  + quarkonia decays into e", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi,jpsi,psi2s mesons
    return signal;
  }
  if (!nameStr.compare("HcToE")) {
    MCProng prong(2, {402, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open charmed hadron decays into e", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("HbToE")) {
    MCProng prong(2, {502, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open beauty hadron decays into e", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("McToE")) {
    MCProng prong(2, {401, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open charmed meson decays into e", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("MbToE")) {
    MCProng prong(2, {501, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open beauty meson decays into e", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("BcToE")) {
    MCProng prong(2, {4001, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open charmed baryon decays into e", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("BbToE")) {
    MCProng prong(2, {5001, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open beauty baryon decay into e", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("PromptHcToE")) {
    MCProng prong(3, {502, 402, 11}, {true, true, true}, {true, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open charmed hadron decays into e", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("NonPromptHcToE")) {
    MCProng prong(3, {502, 402, 11}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "b hadron decays to open charmed hadron decays to e", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("HFdecayToE")) {
    MCProng prong(2, {902, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open charm and beauty to electrons", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("AnyHFdecayToE")) {
    MCProng prong(1, {902}, {true}, {false}, {0}, {0}, {false}, true, {11}, {false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
    prong.SetSignalInTime(true); // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "Open charm and beauty to electrons", {prong}, {-1});
    return signal;
  }
  //   if (!nameStr.compare("LFQtoPC")) {
  //   MCProng prong(3, {900, 22, 11}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
  //   prong.SetSignalInTime(true);  // set direction to check for daughters (true, in time) or for mothers (false, back in time)
  //   signal = new MCSignal(name, "LF meson + quarkonia decays into photon conversion electron", {prong}, {-1});
  //   return signal;
  // }

  //_________________________________________________________________________________________________________________________
  // LMEE pair signals for LF, same mother
  if (!nameStr.compare("eeFromAnything")) {
    MCProng prong(2, {11, MCProng::kPDGCodeNotAssigned}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from any decay", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromPi0")) {
    MCProng prong(2, {11, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from pi0 decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eePrimaryFromPromptPi0")) {
    MCProng prong(2, {11, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502, 402}, {true, true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from prompt pi0 decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromEta")) {
    MCProng prong(2, {11, 221}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from eta decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromEtaprime")) {
    MCProng prong(2, {11, 331}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from eta' decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromRho")) {
    MCProng prong(2, {11, 113}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from rho decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromOmega")) {
    MCProng prong(2, {11, 223}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from omega decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromPhi")) {
    MCProng prong(2, {11, 333}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from phi decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromJpsi")) {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from j/psi decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eePrimaryFromPromptJPsi")) {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from prompt j/psi decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eePrimaryFromNonPromptJPsi")) {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {503}, {false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from non-prompt j/psi decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("mumuFromJpsi")) {
    MCProng prong(2, {13, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "mumu pairs from j/psi decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromPsi2S")) {
    MCProng prong(2, {11, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from psi2s decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("mumuFromPsi2S")) {
    MCProng prong(2, {13, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "mumu pairs from psi2s decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("mumuFromUpsilon1S")) {
    MCProng prong(2, {13, 553}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "mumu pairs from upsilon1s decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("mumuFromUpsilon2S")) {
    MCProng prong(2, {13, 100553}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "mumu pairs from upsilon2s decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("mumuFromUpsilon3S")) {
    MCProng prong(2, {13, 200553}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "mumu pairs from upsilon3s decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromLMeeLFQ")) {
    MCProng prong(2, {11, 900}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from light flavor meson + quarkonia decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromLMeeLF")) {
    MCProng prong(2, {11, 901}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from light flavor meson decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromLMeeNoHFLFQ")) {
    MCProng prong(2, {11, 900}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502, 402}, {true, true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from light flavor meson + quarkonia decays not from open-HF decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromLMeeNoHFLF")) {
    MCProng prong(2, {11, 901}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502, 402}, {true, true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from light flavor meson decays not from open-HF decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }

  // LMEE pair signals for HF
  // D0->e and D0->e
  if (!nameStr.compare("eeFromD0")) {
    MCProng prong(2, {kElectron, Pdg::kD0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from D0 decays, no beauty in history", {prong, prong}, {-1, -1});
    return signal;
  }

  // D0->e and D0->e
  if (!nameStr.compare("eeFromPi0FromD0")) {
    MCProng prong(2, {kElectron, kPi0, Pdg::kD0}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from D0 to Pi0 decays, no beauty in history", {prong, prong}, {1, 1});
    return signal;
  }

  // D+/- -> e and D+/- -> e
  if (!nameStr.compare("eeFromChargedD")) {
    MCProng prong(2, {kElectron, Pdg::kDPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from D+/- decays, no beauty in history", {prong, prong}, {-1, -1});
    return signal;
  }

  // D0 -> e and D+/- -> e
  if (!nameStr.compare("eeFromD0andChargedD")) {
    MCProng prongD0(2, {kElectron, Pdg::kD0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    MCProng prongDch(2, {kElectron, Pdg::kDPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongD0.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongDch.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "one e from D0 and one e from D+/- decays, no beauty in history", {prongD0, prongDch}, {-1, -1});
    return signal;
  }

  // D+/- -> e and D0 -> e
  if (!nameStr.compare("eeFromD0andChargedDBis")) {
    MCProng prongD0(2, {kElectron, Pdg::kD0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    MCProng prongDch(2, {kElectron, Pdg::kDPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongD0.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongDch.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "one e from D+/- and one e from D0 decays, no beauty in history (inverse signal)", {prongDch, prongD0}, {-1, -1});
    return signal;
  }

  // D_s->e and D_s->e
  if (!nameStr.compare("eeFromDs")) {
    MCProng prong(2, {kElectron, Pdg::kDS}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from Ds +/- decays, no beauty in history", {prong, prong}, {-1, -1});
    return signal;
  }

  // Lambda_c->e and Lambda_c->e
  if (!nameStr.compare("eeFromLambdaC")) {
    MCProng prong(2, {kElectron, Pdg::kLambdaCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from Lambda_c, no beauty in history", {prong, prong}, {-1, -1});
    return signal;
  }

  // Lambda_c->e and D0->e
  if (!nameStr.compare("eeFromLambdaCandD0")) {
    MCProng prongLc(2, {kElectron, Pdg::kLambdaCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    MCProng prongD0(2, {kElectron, Pdg::kD0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongLc.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongD0.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "one e from Lambda_c and one e from D0 decays, no beauty in history", {prongLc, prongD0}, {-1, -1});
    return signal;
  }

  // D0->e and Lambda_c->e
  if (!nameStr.compare("eeFromLambdaCandD0Bis")) {
    MCProng prongLc(2, {kElectron, Pdg::kLambdaCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    MCProng prongD0(2, {kElectron, Pdg::kD0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongLc.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongD0.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "one e from Lambda_c and one e from D0 decays, no beauty in history (inverse signal)", {prongD0, prongLc}, {-1, -1});
    return signal;
  }

  // Lambda_c->e and D+/- -> e
  if (!nameStr.compare("eeFromLambdaCandChargedD")) {
    MCProng prongLc(2, {kElectron, Pdg::kLambdaCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    MCProng prongDch(2, {kElectron, Pdg::kDPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongLc.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongDch.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "one e from Lambda_c and one e from D+/- decays, no beauty in history", {prongLc, prongDch}, {-1, -1});
    return signal;
  }

  // D+/- -> e and Lambda_c->e
  if (!nameStr.compare("eeFromLambdaCandChargedDBis")) {
    MCProng prongLc(2, {kElectron, Pdg::kLambdaCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    MCProng prongDch(2, {kElectron, Pdg::kDPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongLc.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongDch.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "one e from Lambda_c and one e from D+/- decays, no beauty in history (inverse signal)", {prongDch, prongLc}, {-1, -1});
    return signal;
  }

  // Xic0 ->e and Xic0 ->e
  if (!nameStr.compare("eeFromXiC0")) {
    MCProng prong(2, {kElectron, Pdg::kXiC0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from Xi_c0, no beauty in history", {prong, prong}, {-1, -1});
    return signal;
  }

  // Xi_c+ ->e and Xi_c+ ->e
  if (!nameStr.compare("eeFromXiCPlus")) {
    MCProng prong(2, {kElectron, Pdg::kXiCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from Xi_c+, no beauty in history", {prong, prong}, {-1, -1});
    return signal;
  }

  // Xi_c0 ->e and Xi_c+ ->e
  if (!nameStr.compare("eeFromXiC0andXiCPlus")) {
    MCProng prongXiCPlus(2, {kElectron, Pdg::kXiCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    MCProng prongXiC0(2, {kElectron, Pdg::kXiC0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongXiCPlus.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongXiC0.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "one e from Xi_c+ and one e from Xi_c0 decays, no beauty in history", {prongXiCPlus, prongXiC0}, {-1, -1});
    return signal;
  }

  // Xi_c+ ->e and Xi_c0 ->e
  if (!nameStr.compare("eeFromXiC0andXiCPlusBis")) {
    MCProng prongXiCPlus(2, {kElectron, Pdg::kXiCPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    MCProng prongXiC0(2, {kElectron, Pdg::kXiC0}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prongXiCPlus.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongXiC0.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "one e from Xi_c+ and one e from Xi_c0 decays, no beauty in history (inverse signal)", {prongXiC0, prongXiCPlus}, {-1, -1});
    return signal;
  }

  // Xi_cc++ ->e and Xi_cc++ ->e
  if (!nameStr.compare("eeFromXiCPlusPlus")) {
    MCProng prong(2, {kElectron, Pdg::kXiCCPlusPlus}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from Xi_cc++, no beauty in history", {prong, prong}, {-1, -1});
    return signal;
  }

  // c->e and c->e (no check)
  if (!nameStr.compare("eeFromCCNoCheck")) {
    MCProng prong(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from c->e and c->e without check", {prong, prong}, {-1, -1}); // signal at pair level
    return signal;
  }

  // ee from HF in general
  if (!nameStr.compare("eeFromHF")) {
    MCProng prong(2, {11, 902}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b,c->e and b,c->e without check", {prong, prong}, {-1, -1}); // signal at pair level
    return signal;
  }

  // Any c in history but no b -> c -> e
  if (!nameStr.compare("eeFromPromptCandPromptC")) {
    MCProng prong(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {true}); // check if mother pdg code is in history
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs with any charm but no beauty in decay chain", {prong, prong}, {-1, -1}); // signal at pair level
    return signal;
  }

  // Any b to any c in history b -> c -> e
  if (!nameStr.compare("eeFromAnyBtoCandAnyBtoC")) {
    MCProng prong(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {false}); // check if mother pdg code is in history
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs with any beauty to charm in decay chain", {prong, prong}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->c->e, b->c->e
  if (!nameStr.compare("eeFromBtoCandBtoC")) {
    MCProng prong(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false); // check if mother pdg code is in history
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs with any beauty to charm in decay chain", {prong, prong}, {-1, -1}); // signal at pair level
    return signal;
  }

  // Any b->e and Any b->X->c->e
  // Looking at such decays: B -> (e) D -> (e)e and bar{B} -> e
  // Signal allows combinations of ee from the same B meson
  // + the combination of e fom B and e from bar{B}
  if (!nameStr.compare("eeFromBandAnyBtoC")) {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {false}); // check if mother pdg code is in history
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->e and b->X->c->e", {prongB, prongBtoC}, {-1, -1}); // signal at pair level
    return signal;
  }

  // Any b->e and Any b->X->c->e
  if (!nameStr.compare("eeFromBandAnyBtoCBis")) {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false}, false, {502}, {false}); // check if mother pdg code is in history
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->X->c->e and b->e", {prongBtoC, prongB}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->e and b->c->e
  if (!nameStr.compare("eeFromBandBtoC")) {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false); // check if mother pdg code is in history
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "direkt ee pairs from b->e and b->c->e", {prongB, prongBtoC}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->e and b->c->e
  if (!nameStr.compare("eeFromBandBtoCBis")) {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false); // check if mother pdg code is in history
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->c->e and b->e", {prongBtoC, prongB}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->e and b->c->e (same mother/grandmother)
  // require that the mother is the grandmother of the other electron
  if (!nameStr.compare("eeFromBandBtoCsameGM")) {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false); // check if mother pdg code is in history
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->e and b->c->e, mother = grandmother", {prongB, prongBtoC}, {1, 2}, false); // signal at pair level, accept commonAncestor Pairs
    return signal;
  }

  // b->e and b->c->e (same mother/grandmother)
  // require that the mother is the grandmother of the other electron
  if (!nameStr.compare("eeFromBandBtoCsameGMBis")) {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false); // check if mother pdg code is in history
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->c->e and b->e, mother = grandmother", {prongBtoC, prongB}, {2, 1}, false); // signal at pair level, accept commonAncestor Pairs
    return signal;
  }

  // b->e and b->c->e (different mother/grandmother)
  // require that the mother is not the grandmother of the other electron
  if (!nameStr.compare("eeFromBandBtoCdiffGM")) {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false); // check if mother pdg code is in history
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->e and b->c->e, mother != grandmother", {prongB, prongBtoC}, {1, 2}, true); // signal at pair level, exclude commonAncestor Pairs
    return signal;
  }

  // b->e and b->c->e (different mother/grandmother)
  // require that the mother is not the grandmother of the other electron
  if (!nameStr.compare("eeFromBandBtoCdiffGMBis")) {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false}, false); // check if mother pdg code is in history
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->c->e and b->e, mother != grandmother", {prongBtoC, prongB}, {2, 1}, true); // signal at pair level, exclude commonAncestor Pairs
    return signal;
  }

  // b->e and b->e
  if (!nameStr.compare("eeFromBB")) {
    MCProng prong(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->e and b->e", {prong, prong}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->e and b->e (commonAncestors)
  if (!nameStr.compare("eeFromSameB")) {
    MCProng prong(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->e and b->e", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }

  // b->e and c->e no check
  if (!nameStr.compare("eeFromBandFromC")) {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongC(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->e and c->e", {prongB, prongC}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->e and c->e no check
  if (!nameStr.compare("eeFromBandFromCBis")) {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongC(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->e and c->e", {prongC, prongB}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->e and b->c->e (single b)
  if (!nameStr.compare("eeFromSingleBandBtoC")) {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prongB.SetSourceBit(0, MCProng::kPhysicalPrimary);
    MCProng prongBtoC(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prongBtoC.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "ee pairs from b->e and b->c->e (single b)", {prongB, prongBtoC}, {1, 2}); // signal at pair level
    return signal;
  }

  //_________________________________________________________________________________________________________________________

  if (!nameStr.compare("kaonFromBplus")) {
    MCProng prong(2, {321, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaons from B+ decays", {prong}, {1});
    return signal;
  }

  if (!nameStr.compare("JpsiFromBplus")) {
    MCProng prong(2, {443, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Jpsi from B+ decays", {prong}, {1});
    return signal;
  }

  if (!nameStr.compare("eFromJpsiFromBplus")) {
    MCProng prong(3, {11, 443, 521}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Electrons from Jpsi from B+ decays", {prong}, {1});
    return signal;
  }

  if (!nameStr.compare("eeFromJpsiFromBplus")) {
    MCProng prong(3, {11, 443, 521}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Electron pair from Jpsi from B+ decays", {prong, prong}, {1, 1});
    return signal;
  }

  if (!nameStr.compare("eeKaonFromBplus")) {
    MCProng pronge(3, {11, 443, 521}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongKaon(2, {321, 521}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Kaon and electron pair from B+", {pronge, pronge, prongKaon}, {2, 2, 1});
    return signal;
  }

  if (!nameStr.compare("Bplus")) {
    MCProng prong(1, {521}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "B+", {prong}, {-1});
    return signal;
  }

  if (!nameStr.compare("beautyPairs")) {
    MCProng prong(1, {503}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal("beautyPairs", "Beauty hadron pair", {prong, prong}, {-1, -1});
    return signal;
  }
  if (!nameStr.compare("everythingFromBeautyPairs")) {
    MCProng prong(2, {0, 503}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal("everythingFromBeautyPairs", "Everything from beauty hadrons pair", {prong, prong}, {-1, -1});
    return signal;
  }
  if (!nameStr.compare("everythingFromEverythingFromBeautyPairsCM")) {
    MCProng prong(3, {0, 0, 503}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal("everythingFromEverythingFromBeautyPairs", "Everything from everything from beauty hadrons pair with common grand-mother", {prong, prong}, {2, 2});
    return signal;
  }
  if (!nameStr.compare("everythingFromBeautyANDeverythingFromEverythingFromBeautyPairs")) {
    MCProng prong1(3, {0, 0, 503}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prong2(2, {0, 503}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal("everythingFromBeautyANDeverythingFromEverythingFromBeautyPairs", "Everything beauty and everything from everything from beauty hadrons pair", {prong1, prong2}, {2, 1});
    return signal;
  }

  //--------------------------------------------------------------------------------

  if (!nameStr.compare("JpsiFromChic0")) {
    MCProng prong(2, {443, 10441}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Jpsi from Chic0 decays", {prong}, {1});
    return signal;
  }
  if (!nameStr.compare("eFromJpsiFromChic0")) {
    MCProng prong(3, {11, 443, 10441}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from Jpsi from Chic0 decays", {prong}, {1});
    return signal;
  }
  if (!nameStr.compare("eeFromJpsiFromChic0")) {
    MCProng prong(3, {11, 443, 10441}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electron pair from Jpsi from Chic0 decays", {prong, prong}, {2, 2});
    return signal;
  }
  if (!nameStr.compare("JpsiFromChic1")) {
    MCProng prong(2, {443, 20443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Jpsi from Chic1 decays", {prong}, {1});
    return signal;
  }
  if (!nameStr.compare("eFromJpsiFromChic1")) {
    MCProng prong(3, {11, 443, 20443}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from Jpsi from Chic1 decays", {prong}, {2});
    return signal;
  }
  if (!nameStr.compare("eeFromJpsiFromChic1")) {
    MCProng prong(3, {11, 443, 20443}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electron pair from Jpsi from Chic1 decays", {prong, prong}, {2, 2});
    return signal;
  }
  if (!nameStr.compare("JpsiFromChic2")) {
    MCProng prong(2, {443, 445}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Jpsi from Chic2 decays", {prong}, {1});
    return signal;
  }
  if (!nameStr.compare("JpsiFromChic2")) {
    MCProng prong(2, {443, 904}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Jpsi from Chic0, Chic1 or Chic2 decays", {prong}, {1});
    return signal;
  }
  if (!nameStr.compare("eFromJpsiFromChic2")) {
    MCProng prong(3, {11, 443, 445}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electrons from Jpsi from Chic2 decays", {prong}, {2});
    return signal;
  }
  if (!nameStr.compare("eeFromJpsiFromChic2")) {
    MCProng prong(3, {11, 443, 445}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electron pair from Jpsi from Chic2 decays", {prong, prong}, {2, 2});
    return signal;
  }
  if (!nameStr.compare("eeFromJpsiFromChic012")) {
    MCProng prong(3, {11, 443, 904}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Electron pair from Jpsi from Chic0, Chic1 or Chic2 decays", {prong, prong}, {2, 2});
    return signal;
  }
  if (!nameStr.compare("PhotonFromChic0")) {
    MCProng prong(2, {22, 10441}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon from Chic0 decays", {prong}, {1});
    return signal;
  }
  if (!nameStr.compare("PhotonFromChic1")) {
    MCProng prong(2, {22, 20443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon from Chic1 decays", {prong}, {1});
    return signal;
  }
  if (!nameStr.compare("PhotonFromChic2")) {
    MCProng prong(2, {22, 445}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon from Chic2 decays", {prong}, {1});
    return signal;
  }
  if (!nameStr.compare("PhotonFromChic012")) {
    MCProng prong(2, {22, 904}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon from Chic0, Chic1, and Chic2 decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("PhotonFromPi0")) {
    MCProng prong(2, {22, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon from Pi0 decays", {prong}, {1});
    return signal;
  }
  if (!nameStr.compare("PhotonPhotonFromPi0")) {
    MCProng prong(2, {22, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon Photon from Pi0 decays", {prong, prong}, {1, 1});
    return signal;
  }

  if (!nameStr.compare("eePhotonFromChic1")) {
    MCProng pronge(3, {11, 443, 20443}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPhoton(2, {22, 20443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    pronge.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongPhoton.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon and electron pair from Chic1", {pronge, pronge, prongPhoton}, {2, 2, 1});
    return signal;
  }

  if (!nameStr.compare("eePhotonFromChic2")) {
    MCProng pronge(3, {11, 443, 445}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPhoton(2, {22, 445}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    pronge.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongPhoton.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon and electron pair from Chic2", {pronge, pronge, prongPhoton}, {2, 2, 1});
    return signal;
  }

  if (!nameStr.compare("eePhotonFromChic12")) {
    MCProng pronge(3, {11, 443, MCProng::kPDGCodeNotAssigned}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    MCProng prongPhoton(2, {22, MCProng::kPDGCodeNotAssigned}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    pronge.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongPhoton.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon and electron pair from Chic0, Chic1 and Chic2", {pronge, pronge, prongPhoton}, {2, 2, 1});
    return signal;
  }

  if (!nameStr.compare("eePhotonFromPi0")) {
    MCProng pronge(2, {11, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongPhoton(2, {22, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    pronge.SetSourceBit(0, MCProng::kPhysicalPrimary);
    prongPhoton.SetSourceBit(0, MCProng::kPhysicalPrimary);
    signal = new MCSignal(name, "Photon and electron pair from Pi0", {pronge, pronge, prongPhoton}, {1, 1, 1});
    return signal;
  }
  return nullptr;
}
