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

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;

struct MFTtoFwd {

  HistogramRegistry registry{
    "registry",
    {
      {"N_FwdForOneMFT", "; N_{trk}^{fwd}; count", {HistType::kTH1F, {{10, -0.5, 10.5}}}},
      {"chi2_for_mult_Fwd", "; chi^2; count", {HistType::kTH1F, {{101, -0.5, 100.5}}}},

    }};

  void init(InitContext&)
  {
      LOGP(info, " --> Init ");

  }

  void process(aod::MFTTrack const& mfttrack, soa::SmallGroups<aod::FwdTracks> const& fwdtracks)
  {

    registry.fill(HIST("N_FwdForOneMFT"), fwdtracks.size());

    if (fwdtracks.size() > 1) {
      for (auto& fwdtrack : fwdtracks) {
        registry.fill(HIST("chi2_for_mult_Fwd"), fwdtrack.chi2MatchMCHMFT());
      }
    }
  }

  // PROCESS_SWITCH(MFTtoFwd, processIndexingFwd, "Create reverse index from particles to tracks", false);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<MFTtoFwd>(cfgc)};
}
