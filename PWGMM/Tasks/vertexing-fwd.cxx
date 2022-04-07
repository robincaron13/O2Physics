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
#include <cmath>
#include "ReconstructionDataFormats/TrackFwd.h"
#include "Math/MatrixFunctions.h"
#include "Math/SMatrix.h"
#include "MathUtils/Utils.h"
#include "Framework/Configurable.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "CommonConstants/MathConstants.h"
#include "CommonConstants/LHCConstants.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::track;

using SMatrix55 = ROOT::Math::SMatrix<double, 5, 5, ROOT::Math::MatRepSym<double, 5>>;
using SMatrix5 = ROOT::Math::SVector<Double_t, 5>;

struct vertexingfwd {

  // Configurable<int> rangeBC{"rangeBC", 10, "Range for collision BCId and BCglobalIndex correspondance"};
  std::vector<int> vecCollForAmb;       // vector for collisions associated to an ambiguous track
  std::vector<double> vecDCACollForAmb; // vector for dca collision associated to an ambiguous track
  std::vector<double> vecAmbTrack;      // vector for dca collision associated to an ambiguous track

  HistogramRegistry registry{
    "registry",
    {{"EventsNtrkZvtx", "; N_{trk}; Z_{vtx}; events", {HistType::kTH2F, {{301, -0.5, 300.5}, {201, -20.1, 20.1}}}}, //
     {"EventSelection", "; status; events", {HistType::kTH1F, {{4, 0.5, 4.5}}}},                                    //
     {"TrackDCAxy", "; DCA_{xy}; counts", {HistType::kTH1F, {{100, -10, 10}}}},                                     //
     {"TrackDCAx", "; DCA_{x}; counts", {HistType::kTH1F, {{100, -10, 10}}}},                                       //
     {"TrackDCAy", "; DCA_{y}; counts", {HistType::kTH1F, {{100, -10, 10}}}},                                       //
     {"CollisionTime", "; Collisions time (ns); counts", {HistType::kTH1F, {{100, 0, 0.1}}}},
     {"CollisionTimeRes", "; Resolution of collision time (ns); counts", {HistType::kTH1F, {{100, 0, 0.1}}}},
     {"Chi2", "; #chi^{2}; counts", {HistType::kTH1F, {{100, 0, 10}}}},
     {"vecCollSize", "; coll size", {HistType::kTH1F, {{100, 0, 100}}}},
     {"IndicesCollMC", "; i CollMC; i CollAmbi", {HistType::kTH2F, {{401, -0.5, 1100.5}, {401, -0.5, 1100.5}}}},
     {"NumContrib", "; N_{tr} used for the vertex; counts", {HistType::kTH1F, {{100, 0, 100}}}}} //
  };

  void process(aod::AmbiguousMFTTracks const& ambitracks, aod::BCs const& bcs, soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels> const& tracks, soa::Join<aod::Collisions, aod::McCollisionLabels> const& collisions, aod::McParticles const& mcParticles, aod::McCollisions const& mcCollisions) // AmbiguousMFTTracks and fwd doesn't work yet
  {
    for (auto& ambitrack : ambitracks) {
      vecCollForAmb.clear();
      vecDCACollForAmb.clear();
      vecAmbTrack.clear();
      if (!ambitrack.bc().size()) {
        continue;
      }
      if (!ambitrack.globalIndex()) {
        continue;
      }
      if (!tracks.size()) {
        continue;
      }
      LOGF(info, "------------------------------------ We look at ambitrack %d which has %d possible BCs", ambitrack.globalIndex(), ambitrack.bc().size());
      if (tracks.size() <= 1)
        continue;
      if (ambitrack.globalIndex() >= tracks.size()) {
        continue;
      }
      if (ambitrack.globalIndex() < 0) {
        continue;
      }
      LOGF(info, "CI %d,  tracks.size() %d \n", ambitrack.globalIndex(), tracks.size());
      int mcCollAmbiID = -1;
      // soa::Join<o2::aod::FwdTracks, o2::aod::FwdTracksCov, aod::McFwdTrackLabels>::iterator extAmbiTrack;
      soa::Join<o2::aod::MFTTracks, aod::McMFTTrackLabels>::iterator extAmbiTrack;

      // auto extAmbiTrack = tracks.iteratorAt(ambitrack.globalIndex());
      auto extAmbiTrackid = ambitrack.mfttrackId();
      bool itIsDef = false;
      for (auto& track : tracks) {
        if (extAmbiTrackid == track.globalIndex()) {
          if (!track.has_mcParticle()) {
            LOGF(warning, "No MC particle for ambiguous track, skip...");
            continue;
          }
          auto particle = track.mcParticle();
          mcCollAmbiID = particle.mcCollisionId();
          vecAmbTrack.push_back(track.x());
          vecAmbTrack.push_back(track.y());
          vecAmbTrack.push_back(track.phi());
          vecAmbTrack.push_back(track.tgl());
          vecAmbTrack.push_back(track.signed1Pt());
          vecAmbTrack.push_back(track.z());
          vecAmbTrack.push_back(track.chi2());

          extAmbiTrack = track; // can not access to the extAmbiTrack .x, .y(), .phi() ...
          itIsDef = true;
          break;
        }
      }
      // if( track.globalIndex() == (tracks.size()-1) ) continue;
      if (!extAmbiTrack.x()) {
        continue;
      }
      // printf("extAmbiTrack %d \n", ambitrack.globalIndex());
      // LOGF(info, "extAmbiTrack %d \n", extAmbiTrack.phi());
      //      if (!extAmbiTrack.has_mcParticle()) {
      //        LOGF(warning, "No MC particle for ambiguous track, skip...");
      //        continue;
      //      }
      //      auto particle = extAmbiTrack.mcParticle();
      //      int mcCollAmbiID = particle.mcCollisionId();

      // printf("phi = %f, mcCollAmbiID  %d \n", extAmbiTrack.phi(), mcCollAmbiID);

      for (auto& bc : ambitrack.bc()) {
        // LOGF(info, "       5) BC %d with global BC %lld", bc.globalIndex(), bc.globalBC());

        for (auto& collision : collisions) {
          int rangeBC = 100; // retrive the global bc for the indices given in the bcslice vector
          registry.fill(HIST("EventSelection"), 1.);

          if ((collision.bc().globalBC() > (bc.globalBC() - rangeBC)) && (collision.bc().globalBC() < (bc.globalBC() + rangeBC))) {
            // registry.fill(HIST("EventsNtrkZvtx"), tracks.size(), collision.posZ());
            // registry.fill(HIST("EventSelection"), 2.);
            //  printf("collision BC ID = %d, bc.globalIndex() %lld\n", collision.bcId(), bc.globalIndex());
            // LOGF(info,"collision pos Z %f\n", collision.posZ());

            SMatrix5 tpars(vecAmbTrack[0], vecAmbTrack[1], vecAmbTrack[2], vecAmbTrack[3], vecAmbTrack[4]);
            // printf("SMatrix5  extAmbiTrack.signed1Pt() %f\n", extAmbiTrack.signed1Pt());
            //            std::vector<double> v1{extAmbiTrack.cXX(), extAmbiTrack.cXY(), extAmbiTrack.cYY(), extAmbiTrack.cPhiX(), extAmbiTrack.cPhiY(),
            //                                   extAmbiTrack.cPhiPhi(), extAmbiTrack.cTglX(), extAmbiTrack.cTglY(), extAmbiTrack.cTglPhi(), extAmbiTrack.cTglTgl(),
            //                                   extAmbiTrack.c1PtX(), extAmbiTrack.c1PtY(), extAmbiTrack.c1PtPhi(), extAmbiTrack.c1PtTgl(), extAmbiTrack.c1Pt21Pt2()};
            std::vector<double> v1;
            // printf("vector<double>  extAmbiTrack.cXX() %f\n", extAmbiTrack.cXX());

            SMatrix55 tcovs(v1.begin(), v1.end());

            o2::track::TrackParCovFwd pars1{vecAmbTrack[5], tpars, tcovs, vecAmbTrack[6]};
            // double chi2 = 1.0;
            // o2::track::TrackParCovFwd pars1{extAmbiTrack.z(), tpars, tcovs, chi2};
            // printf("SMatrix55  extAmbiTrack.chi2() %f\n", extAmbiTrack.chi2());
            pars1.propagateToZlinear(collision.posZ());
            // printf("propagateToZlinear  collision.posZ() %f\n", collision.posZ());
            const auto dcaX(pars1.getX() - collision.posX());
            const auto dcaY(pars1.getY() - collision.posY());
            auto dcaXY = std::sqrt(dcaX * dcaX + dcaY * dcaY);
            // printf("dcaXY  %f\n", dcaXY);
            registry.fill(HIST("TrackDCAxy"), dcaXY);
            registry.fill(HIST("TrackDCAx"), dcaX);
            registry.fill(HIST("TrackDCAy"), dcaY);
            registry.fill(HIST("CollisionTime"), collision.collisionTime());
            registry.fill(HIST("CollisionTimeRes"), collision.collisionTimeRes());
            registry.fill(HIST("Chi2"), collision.chi2());
            registry.fill(HIST("NumContrib"), collision.numContrib());

            int mcCollindex = collision.mcCollision().globalIndex();
            vecCollForAmb.push_back(mcCollindex);
            vecDCACollForAmb.push_back(dcaXY);
          } // condition BC
        }   // collisions

        int indexMinDCA = std::distance(vecDCACollForAmb.begin(), std::min_element(vecDCACollForAmb.begin(), vecDCACollForAmb.end()));
        int indexMCcoll = vecCollForAmb[indexMinDCA];
        registry.fill(HIST("IndicesCollMC"), mcCollAmbiID, indexMCcoll);
        registry.fill(HIST("vecCollSize"), vecCollForAmb.size());
        if (mcCollAmbiID == indexMCcoll) {
          LOGF(info, " ------> ambitrack correctly associated to collision \n");
        }
      } // bc of ambitracks
    }   // ambitracks
  }

  PROCESS_SWITCH(vertexingfwd, process, "Process ambiguous track DCA", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<vertexingfwd>(cfgc)};
}
