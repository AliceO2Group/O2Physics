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

#ifndef PWGDQ_CORE_MCSIGNALLIBRARY_H_
#define PWGDQ_CORE_MCSIGNALLIBRARY_H_

#include <string>
#include "PWGDQ/Core/MCProng.h"
#include "PWGDQ/Core/MCSignal.h"

namespace o2::aod
{
namespace dqmcsignals
{
MCSignal* GetMCSignal(const char* signalName);
}
} // namespace o2::aod

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
    MCProng prong(3, {11, 443, 503}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Electrons from beauty jpsi decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromPromptJpsi")) {
    MCProng prong(3, {11, 443, 503}, {true, true, true}, {false, false, true}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Electrons from prompt jpsi decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("Jpsi")) {
    MCProng prong(1, {443}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "Inclusive jpsi", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("nonPromptJpsi")) {
    MCProng prong(2, {443, 503}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Non-prompt jpsi", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("promptJpsi")) {
    MCProng prong(2, {443, 503}, {true, true}, {false, true}, {0, 0}, {0, 0}, {false, false});
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
  if (!nameStr.compare("anyBeautyHadron")) {
    MCProng prong(1, {503}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal(name, "All beauty hadrons", {prong}, {-1});
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
  if (!nameStr.compare("electronFromPi0")) {
    MCProng prong(2, {11, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from pi0 decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("Pi0decayTOe")) {
    MCProng prong(2, {111, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
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
    signal = new MCSignal("dielectron", "Electron pair", {prong, prong}, {-1, -1});
    return signal;
  }
  if (!nameStr.compare("dimuon")) {
    MCProng prong(1, {13}, {true}, {false}, {0}, {0}, {false});
    signal = new MCSignal("dielectron", "Electron pair", {prong, prong}, {-1, -1});
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
  if (!nameStr.compare("dielectronPCPi0")) {
    MCProng prong(3, {11, 22, 111}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "dielectron from a photon conversion from a pi0", {prong, prong}, {1, 1});
    return signal;
  }

  // LMEE single signals
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
  if (!nameStr.compare("muFromJpsi")) {
    MCProng prong(2, {13, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "muons from jpsi decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromPsi2S")) {
    MCProng prong(2, {11, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from psi2s decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("muFromPsi2S")) {
    MCProng prong(2, {13, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "muons from psi2s decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromLMeeLFQ")) {
    MCProng prong(2, {11, 900}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
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
  if (!nameStr.compare("LFQdecayToE")) {
    MCProng prong(2, {900, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
    prong.SetSignalInTime(true);                                                       // set direction to check for daughters (true, in time) or for mothers (false, back in time)
    signal = new MCSignal(name, "LF meson  + quarkonia decays into e", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi,jpsi,psi2s mesons
    return signal;
  }
  if (!nameStr.compare("eFromLMeeLF")) {
    MCProng prong(2, {11, 901}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    // prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);  // set source to be ALICE primary particles
    signal = new MCSignal(name, "Electrons from LF meson decays", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi mesons
    return signal;
  }
  if (!nameStr.compare("ePrimaryFromLMeeLF")) {
    MCProng prong(2, {11, 901}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary, false);                      // set source to be ALICE primary particles
    signal = new MCSignal(name, "Electrons from LF meson decays", {prong}, {-1}); // pi0,eta,eta',rho,omega,phi mesons
    return signal;
  }
  if (!nameStr.compare("eFromHc")) {
    MCProng prong(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from open charmed hadron decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromHb")) {
    MCProng prong(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from open beauty hadron decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromMc")) {
    MCProng prong(2, {11, 401}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from open charmed meson decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromMb")) {
    MCProng prong(2, {11, 501}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from open beauty meson decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromBc")) {
    MCProng prong(2, {11, 4001}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from open charmed baryon decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromBb")) {
    MCProng prong(2, {11, 5001}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "Electrons from open beauty baryon decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromPromptHc")) {
    MCProng prong(3, {11, 402, 502}, {true, true, true}, {false, false, true}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Electrons from open charmed hadron decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("eFromNonPromptHc")) {
    MCProng prong(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "Electrons from open charmed hadron decays from b hadron decays", {prong}, {-1});
    return signal;
  }
  if (!nameStr.compare("HFdecayToE")) {
    MCProng prong(2, {902, 11}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
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
  if (!nameStr.compare("eeFromPi0")) {
    MCProng prong(2, {11, 111}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "ee pairs from pi0 decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromEta")) {
    MCProng prong(2, {11, 221}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "ee pairs from eta decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromEtaprime")) {
    MCProng prong(2, {11, 331}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "ee pairs from eta' decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromRho")) {
    MCProng prong(2, {11, 113}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "ee pairs from rho decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromOmega")) {
    MCProng prong(2, {11, 223}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "ee pairs from omega decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromPhi")) {
    MCProng prong(2, {11, 333}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "ee pairs from phi decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromJpsi")) {
    MCProng prong(2, {11, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "ee pairs from j/psi decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("mumuFromJpsi")) {
    MCProng prong(2, {13, 443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "mumu pairs from j/psi decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromPsi2S")) {
    MCProng prong(2, {11, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "ee pairs from psi2s decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("mumuFromPsi2S")) {
    MCProng prong(2, {13, 100443}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "mumu pairs from psi2s decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromLMeeLFQ")) {
    MCProng prong(2, {11, 900}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "ee pairs from light flavor meson + quarkonia decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eeFromLMeeLF")) {
    MCProng prong(2, {11, 901}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "ee pairs from light flavor meson decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eePrimaryFromLMeeLFQ")) {
    MCProng prong(2, {11, 900}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);                                                           // set source to be ALICE primary particles
    signal = new MCSignal(name, "ee pairs from light flavor meson + quarkonia decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }
  if (!nameStr.compare("eePrimaryFromLMeeLF")) {
    MCProng prong(2, {11, 901}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    prong.SetSourceBit(0, MCProng::kPhysicalPrimary);                                               // set source to be ALICE primary particles
    signal = new MCSignal(name, "ee pairs from light flavor meson decays", {prong, prong}, {1, 1}); // signal at pair level
    return signal;
  }

  // LMEE pair signals for HF
  // c->e and c->e (no check)
  if (!nameStr.compare("eeFromCCNoCheck")) {
    MCProng prong(2, {11, 402}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "ee pairs from c->e and c->e without check", {prong, prong}, {-1, -1}); // signal at pair level
    return signal;
  }

  // c->e and c->e (prompt)
  if (!nameStr.compare("eeFromCC")) {
    MCProng prong(3, {11, 402, 502}, {true, true, true}, {false, false, true}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "ee pairs from c->e and c->e", {prong, prong}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->e and b->e
  if (!nameStr.compare("eeFromBB")) {
    MCProng prong(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    signal = new MCSignal(name, "ee pairs from b->e and b->e", {prong, prong}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->c->e and b->c->e
  if (!nameStr.compare("eeFromBtoC")) {
    MCProng prong(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "ee pairs from b->c->e and b->c->e", {prong, prong}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->e and b->c->e
  if (!nameStr.compare("eeFromBandBtoC")) {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongBtoC(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "ee pairs from b->e and b->c->e", {prongB, prongBtoC}, {-1, -1}); // signal at pair level
    return signal;
  }

  // b->e and b->c->e (single b)
  if (!nameStr.compare("eeFromSingleBandBtoC")) {
    MCProng prongB(2, {11, 502}, {true, true}, {false, false}, {0, 0}, {0, 0}, {false, false});
    MCProng prongBtoC(3, {11, 402, 502}, {true, true, true}, {false, false, false}, {0, 0, 0}, {0, 0, 0}, {false, false, false});
    signal = new MCSignal(name, "ee pairs from b->e and b->c->e (single b)", {prongB, prongBtoC}, {1, 2}); // signal at pair level
    return signal;
  }

  //_________________________________________________________________________________________________________________________

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
  return nullptr;
}
#endif // PWGDQ_CORE_MCSIGNALLIBRARY_H_
