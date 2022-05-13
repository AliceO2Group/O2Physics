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
/// \file LFNucleiBATask.cxx
///
/// \author Rutuparna Rath <rutuparna.rath@cern.ch> and Giovanni Malfattore <giovanni.malfattore@cern.ch>
///

#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/HistogramRegistry.h"

#include "Common/Core/PID/PIDResponse.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/trackUtilities.h"
#include "PWGLF/DataModel/LFNucleiTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
//class o2::aod::LfTreeCreatorNuclei;
struct LFNucleiDeuteronTask {
  OutputObj<TH1F> h1VtxZ{TH1F("h1VtxZ", "V_{z};V_{z} (in cm); counts", 3000, -15, 15)};
  OutputObj<TH1F> h1CentV0M{TH1F("h1CentV0M", "V0M; Multiplicity; counts", 110, 0, 110)};
  OutputObj<TH2F> h2PionVspNSigmaTPC{TH2F(" h2PionVspNSigmaTPC", "NSigmaTPC(pi) vs p; Momentum(p) in GeV; NSigmaTPC", 200, 0, 20., 200, -10, 10.)};
  OutputObj<TH2F> h2KaonVspNSigmaTPC{TH2F(" h2KaonVspNSigmaTPC", "NSigmaTPC(Ka) vs p; Momentum(p) in GeV; NSigmaTPC", 200, 0, 20., 200, -10, 10.)};
  OutputObj<TH2F> h2ProtonVspNSigmaTPC{TH2F(" h2ProtonVspNSigmaTPC", "NSigmaTPC(proton) vs p; Momentum(p) in GeV; NSigmaTPC", 200, 0, 20., 200, -10, 10.)};
  OutputObj<TH2F> h2DeuteronVspNSigmaTPC{TH2F(" h2DeuteronVspNSigmaTPC", "NSigmaTPC(D) vs p; Momentum(p) in GeV; NSigmaTPC", 200, 0, 20., 200, -10, 10.)};
  OutputObj<TH2F> h2HeliumVspNSigmaTPC{TH2F(" h2HeliumVspNSigmaTPC", "NSigmaTPC(He) vs p; Momentum(p) in GeV; NSigmaTPC", 200, 0, 20., 200, -10, 10.)};
  OutputObj<TH2F> h2PionVspNSigmaTOF{TH2F(" h2PionVspNSigmaTOF", "NSigmaTOF(pi) vs p; Momentum(p) in GeV; NSigmaTOF", 200, 0, 20., 200, -10, 10.)};
  OutputObj<TH2F> h2KaonVspNSigmaTOF{TH2F(" h2KaonVspNSigmaTOF", "NSigmaTOF(Ka) vs p; Momentum(p) in GeV; NSigmaTOF", 200, 0, 20., 200, -10, 10.)};
  OutputObj<TH2F> h2ProtonVspNSigmaTOF{TH2F(" h2ProtonVspNSigmaTOF", "NSigmaTOF(proton) vs p; Momentum(p) in GeV; NSigmaTOF", 200, 0, 20., 200, -10, 10.)};
  OutputObj<TH2F> h2DeuteronVspNSigmaTOF{TH2F(" h2DeuteronVspNSigmaTOF", "NSigmaTOF(D) vs p; Momentum(p) in GeV; NSigmaTOF", 200, 0, 20., 200, -10, 10.)};
  OutputObj<TH2F> h2HeliumVspNSigmaTOF{TH2F(" h2HeliumVspNSigmaTOF", "NSigmaTOF(He) vs p; Momentum(p) in GeV; NSigmaTOF", 200, 0, 20., 200, -10, 10.)};
  OutputObj<TH1F> h1DeuteronSpectra{TH1F("h1DeuteronSpectra", "pT; counts", 200, 0, 20)};

  void init(o2::framework::InitContext&)
  {
  }
  //void process(soa::Join<o2::aod::Collisions, o2::aod::LfCandNucleusFullEvents>::iterator const& events)
  void process(o2::aod::LfCandNucleusFullEvents::iterator const& event, o2::aod::LfCandNucleusFull const& tracks)
  {

    //for(auto& coll:events)
    h1VtxZ->Fill(event.vz());
    h1CentV0M->Fill(event.v0m());
    for (auto& track : tracks) {
      h2PionVspNSigmaTPC->Fill(track.p(), track.nsigTPCPi());
      h2KaonVspNSigmaTPC->Fill(track.p(), track.nsigTPCKa());
      h2ProtonVspNSigmaTPC->Fill(track.p(), track.nsigTPCPr());
      h2DeuteronVspNSigmaTPC->Fill(track.p(), track.nsigTPCD());
      h2HeliumVspNSigmaTPC->Fill(track.p(), track.nsigTPC3He());
      h2PionVspNSigmaTOF->Fill(track.p(), track.nsigTOFPi());
      h2KaonVspNSigmaTOF->Fill(track.p(), track.nsigTOFKa());
      h2ProtonVspNSigmaTOF->Fill(track.p(), track.nsigTOFPr());
      h2DeuteronVspNSigmaTOF->Fill(track.p(), track.nsigTOFD());
      h2HeliumVspNSigmaTOF->Fill(track.p(), track.nsigTOF3He());
      //if(std::abs(track.nsigTPCD()) < 5. && std::abs(track.nsigTOFD()) < 3.)
      if (std::abs(track.nsigTPCD()) < 5.)
        h1DeuteronSpectra->Fill(track.pt());
    }
    //LOG(info)<<"Vertex Z ==="<<coll.posZ();
  }
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<LFNucleiDeuteronTask>(cfgc)};
}
