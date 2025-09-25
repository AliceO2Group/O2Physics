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
/// \file   checkPidPacking.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  exec to check the packing and depacking of PID signals (nsigmas)
/// \since  03/05/2024
///

#include "Common/DataModel/PIDResponseTOF.h"
#include "Common/DataModel/PIDResponseTPC.h"

#include <Framework/Logger.h>

#include <TCanvas.h>
#include <TH1.h>
#include <TRandom.h>
#include <TString.h>

#include <cstdint>
#include <string>

using namespace o2;

template <typename T>
bool process(std::string outputName, int nevents = 100000)
{
  class NsigmaContainer
  {
   public:
    NsigmaContainer() {}
    void operator()(const int8_t& packed) { mPacked = packed; }
    int8_t mPacked = 0;
    float unpack() { return T::unPackInTable(mPacked); }
  } container;

  TH1F* hgaus = new TH1F("hgaus", "", 20 / T::bin_width,
                         -10 + T::bin_width * 0.5,
                         10 + T::bin_width * 0.5);
  LOG(info) << "Bin width = " << T::bin_width << " vs histo " << hgaus->GetXaxis()->GetBinWidth(1);
  hgaus->SetLineColor(2);
  hgaus->SetLineStyle(1);
  TH1F* hgausPacked = static_cast<TH1F*>(hgaus->Clone("hgausPacked"));
  hgausPacked->SetLineColor(4);
  hgausPacked->SetLineStyle(2);

  TH1F* huniform = static_cast<TH1F*>(hgaus->Clone("huniform"));
  huniform->SetLineColor(2);
  huniform->SetLineStyle(1);
  TH1F* huniformPacked = static_cast<TH1F*>(hgaus->Clone("huniformPacked"));
  huniformPacked->SetLineColor(4);
  huniformPacked->SetLineStyle(2);

  for (int i = 0; i < nevents; i++) {
    float nsigma = gRandom->Gaus(0, 1);
    hgaus->Fill(nsigma);
    T::packInTable(nsigma, container);
    hgausPacked->Fill(container.unpack());

    nsigma = gRandom->Uniform(-10, 10);
    huniform->Fill(nsigma);
    T::packInTable(nsigma, container);
    huniformPacked->Fill(container.unpack());
  }

  TCanvas* can = new TCanvas("can");
  hgaus->Draw();
  hgausPacked->Draw("same");
  outputName = "/tmp/" + outputName + ".pdf";
  can->SaveAs(Form("%s[", outputName.c_str()));
  can->SaveAs(outputName.c_str());

  huniform->Draw();
  huniformPacked->Draw("same");
  can->SaveAs(outputName.c_str());
  can->SaveAs(Form("%s]", outputName.c_str()));
  const bool gausOk = (hgaus->GetBinContent(hgaus->FindBin(0)) == hgausPacked->GetBinContent(hgausPacked->FindBin(0)));
  const bool uniformOk = (huniform->GetBinContent(huniform->FindBin(0)) == huniformPacked->GetBinContent(huniformPacked->FindBin(0)));
  return gausOk && uniformOk;
}

int main(int /*argc*/, char* /*argv*/[])
{

  LOG(info) << "Checking the packing and unpacking of PID signals (nsigmas) in the TPC PID response.";
  if (process<aod::pidtpc_tiny::binning>("TPC", 100000)) {
    LOG(info) << "Packing and unpacking of PID signals (nsigmas) in the TPC PID response is correct.";
  } else {
    LOG(fatal) << "Packing and unpacking of PID signals (nsigmas) in the TPC PID response is incorrect.";
  }

  LOG(info) << "Checking the packing and unpacking of PID signals (nsigmas) in the TOF PID response.";
  if (process<aod::pidtof_tiny::binning>("TOF", 100000)) {
    LOG(info) << "Packing and unpacking of PID signals (nsigmas) in the TOF PID response is correct.";
  } else {
    LOG(fatal) << "Packing and unpacking of PID signals (nsigmas) in the TOF PID response is incorrect.";
  }

} // main
