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
#include "ReconstructionDataFormats/TrackFwd.h"

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::track;

struct dummymc {

  void process(soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov, aod::McFwdTrackLabels> const& tracks, aod::McParticles const& mcParticles) // AmbiguousMFTTracks and fwd doesn't work yet
  {
    for (auto& track : tracks) {
      LOGF(info, "CI %d,  tracks.size() %d \n", track.globalIndex(), tracks.size());
      auto extAmbiTrack = tracks.iteratorAt(track.globalIndex());
      //auto extAmbiTrack = tracks.fwdtrack();
      if (!extAmbiTrack.phi()) {
        continue;
      }
      LOGF(info, "extAmbiTrack %d \n", extAmbiTrack.phi());
      //      if (!extAmbiTrack.has_mcParticle()) {
      //        LOGF(warning, "No MC particle for ambiguous track, skip...");
      //        continue;
      //      }
      auto particle = extAmbiTrack.mcParticle();
      int mcCollAmbiID = particle.mcCollisionId();

    }   // ambitracks
  }

};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<dummymc>(cfgc)};
}
