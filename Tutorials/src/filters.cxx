// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"

namespace o2::aod
{
namespace etaphi
{
DECLARE_SOA_INDEX_COLUMN(Collision, collision);
DECLARE_SOA_COLUMN(Eta, eta, float, "fEta");
DECLARE_SOA_COLUMN(Phi, phi, float, "fPhi");
} // namespace etaphi
DECLARE_SOA_TABLE(EtaPhi, "AOD", "ETAPHI",
                  etaphi::Eta, etaphi::Phi);

DECLARE_SOA_TABLE(EtaPhiCol, "AOD", "ETAPHICOL",
                  etaphi::CollisionId, etaphi::Eta, etaphi::Phi);
} // namespace o2::aod

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// This is a very simple example showing how to iterate over tracks
// and create a new collection for them.
// FIXME: this should really inherit from AnalysisTask but
//        we need GCC 7.4+ for that
struct ATask {
  Produces<aod::EtaPhi> etaphi;

  void process(aod::Tracks const& tracks)
  {
    for (auto& track : tracks) {
      etaphi(track.eta(), track.phi());
    }
  }
};

struct BTask {
  Filter ptFilter = (aod::etaphi::phi > 1.f) && (aod::etaphi::phi < 2.f);

  void process(soa::Filtered<aod::EtaPhi> const& etaPhis)
  {
    for (auto& etaPhi : etaPhis) {
      LOGF(INFO, "(%f, 1 < %f < 2)", etaPhi.eta(), etaPhi.phi());
    }
  }
};

struct CTask {
  Produces<aod::EtaPhiCol> epc;

  void process(aod::Collision const& collision, aod::Tracks const& tracks)
  {
    for (auto& track : tracks) {
      epc(collision, track.eta(), track.phi());
    }
  }
};

struct DTask {
  Filter flt = (aod::etaphi::eta < 1.0f) && (aod::etaphi::eta > -1.0f);

  void process(aod::Collision const& collision, soa::Filtered<aod::EtaPhiCol> const& epc)
  {
    uint32_t count = 0;
    LOG(INFO) << "===========================================================";
    for (auto& entry : epc) {
      LOG(INFO) << "[" << count++ << "] " << collision.globalIndex() << " = " << entry.collisionId() << "\n"
                << "(-1 < " << entry.eta() << " < 1)";
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const&)
{
  return WorkflowSpec{
    adaptAnalysisTask<ATask>("produce-etaphi"),
    adaptAnalysisTask<BTask>("consume-etaphi"),
    adaptAnalysisTask<CTask>("produce-etaphi-col"),
    adaptAnalysisTask<DTask>("consume-etaphi-col")};
}
