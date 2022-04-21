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
/// \brief By default the histogram is saved to file
///        AnalysisResults.root.
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
//#include <vector>
#include <utility>
#include "TComplex.h"
#include "Framework/Configurable.h"
#include "PWGCF/GenericFramework/GFW.h"
#include "PWGCF/GenericFramework/GFWCumulant.h"
#include "PWGCF/GenericFramework/FlowContainer.h"
#include "PWGCF/GenericFramework/GFWWeights.h"
#include <TProfile.h>
#include <TRandom3.h>

//#include "Common/DataModel/EventSelection.h"
//#include "Common/DataModel/TrackSelectionTables.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using std::vector;
using MyMCTracks = soa::Join<aod::Tracks, aod::McTrackLabels>;
using MyTracks = aod::Tracks;
using MyTracksSelected = aod::FwdTracks;
using MyCollisions = aod::Collisions::iterator;

const Int_t maxHarmonic = 4;
const Int_t maxPower = 1;
TComplex Qvector[maxHarmonic][maxPower];    // Q-vector components
TComplex QvectorPos[maxHarmonic][maxPower]; // Q-vector components with positive eta range
TComplex QvectorNeg[maxHarmonic][maxPower]; // Q-vector components with negative eta range
// const float etaLimit = -3.05;

struct RootHistograms {

  // normal creation of a histogram
  TH1F* phiHA = new TH1F("phiA", "phiA", 100, 0., 2. * M_PI);
  TH1F* etaHA = new TH1F("etaA", "etaA", 102, -4.01, 4.01);

  void process(aod::Tracks const& tracks)
  {
    for (auto& track : tracks) {
      phiHA->Fill(track.phi());
      etaHA->Fill(track.eta());
    }
  }
};

struct OutputObjects {

  // histogram created with OutputObj<TH1F>
  OutputObj<TH1F> phiB{TH1F("phiB", "phiB", 100, 0., 2. * M_PI), OutputObjHandlingPolicy::QAObject};
  OutputObj<TH1F> Nch{TH1F("Nch", "Nch", 100, 0., 100.), OutputObjHandlingPolicy::QAObject};

  OutputObj<TH2F> etaptB{TH2F("etaptB", "etaptB", 102, -2.01, 2.01, 100, 0.0, 5.0), OutputObjHandlingPolicy::AnalysisObject};
  OutputObj<TH2F> trXY{TH2F("trXY", "trXY", 100, -40, 40, 100, -40, 40), OutputObjHandlingPolicy::QAObject};
  OutputObj<TH2F> trXY0{TH2F("trXY0", "trXY0", 200, -100, 100, 100, -20, 20), OutputObjHandlingPolicy::QAObject};
  OutputObj<TH2F> trXY4{TH2F("trXY4", "trXY4", 110, -55, 55, 100, -20, 20), OutputObjHandlingPolicy::QAObject};

  void process(aod::Tracks const& tracks)
  {
    for (auto& track : tracks) {
      phiB->Fill(track.phi());
      etaptB->Fill(track.eta(), track.pt());
      trXY->Fill(track.x(), track.y());
      trXY0->Fill(track.z(), track.x());
      trXY4->Fill(track.z(), track.y());
    }
    Nch->Fill(tracks.size());
  }
};

// Simple access to collision
struct VertexDistribution {
  OutputObj<TH1F> vertex{TH1F("vertex", "vertex", 100, -10, 10)};

  // loop over MC truth McCollisions
  void process(aod::McCollision const& mcCollision)
  {
    LOGF(info, "MC. vtx-z = %f", mcCollision.posZ());
    vertex->Fill(mcCollision.posZ());
  }
};

/*
struct MultiplicityEventTrackSelection {

  OutputObj<TH1F> multiplicity{TH1F("multiplicity", "multiplicity", 5000, -0.5, 4999.5), OutputObjHandlingPolicy::AnalysisObject};

  Filter collisionZFilter = nabs(aod::collision::posZ) < 10.0f;
  Filter trackFilter = (nabs(aod::track::eta) < 0.8f) && (aod::track::pt > 0.15f) && (aod::track::isGlobalTrack == (uint8_t) true);


  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
               soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>> const& tracks)
  {
    if (!collision.alias()[kINT7]) {
      return;
    }
    if (!collision.sel7()) {
      return;
    }

    LOGP(INFO, "Collision with {} tracks", tracks.size());
    multiplicity->Fill(tracks.size());
  }
};

struct MultiplicityEventTrackSelectionMFT {

  OutputObj<TH1F> multiplicityMFT{TH1F("multiplicityMFT", "multiplicityMFT", 5000, -0.5, 4999.5), OutputObjHandlingPolicy::AnalysisObject};

  Filter collisionZFilter = nabs(aod::collision::posZ) < 30.0f;
  Filter trackFilter = (aod::track::eta > -3.6f) && (aod::track::eta < -2.5f) && (aod::track::pt > 0.15f) && (aod::track::isGlobalTrack == (uint8_t) true);

  void process(soa::Filtered<soa::Join<aod::Collisions, aod::EvSels>>::iterator const& collision,
               soa::Filtered<soa::Join<aod::Tracks, aod::TrackSelection>> const& tracks)
  {
    if (!collision.alias()[kINT7]) {
      return;
    }
    if (!collision.sel7()) {
      return;
    }

    LOGP(INFO, "Collision with {} MFT tracks", tracks.size());
    multiplicityMFT->Fill(tracks.size());
  }
};
*/

struct OutputObjSet {
  // incomplete definition of an OutputObj
  OutputObj<TH1F> trZ{"trZ", OutputObjHandlingPolicy::QAObject};
  OutputObj<TH2F> phieta{"phieta", OutputObjHandlingPolicy::QAObject};
  OutputObj<TH2F> ZvtxEta{"ZvtxEta", OutputObjHandlingPolicy::QAObject};
  OutputObj<TH2F> trackZeta{"trackZeta", OutputObjHandlingPolicy::QAObject};

  // Filter ptfilter = aod::track::pt > 0.000001f;

  // Filter trackFilter = (aod::track::eta > -3.6f) && (aod::track::eta < -2.5f) ;

  void init(InitContext const&)
  {
    // complete the definition of the OutputObj
    trZ.setObject(new TH1F("Z", "Z", 100, -100., 100.));
    phieta.setObject(new TH2F("phieta", "MFT tracks; #varphi; #eta;", 102, -2 * M_PI, 2 * M_PI, 102, -4.51, 4.51));
    ZvtxEta.setObject(new TH2F("ZvtxEta", "MFT tracks; z_{vtx}; #eta ", 102, -50., 50., 102, -4.51, 4.51));
    trackZeta.setObject(new TH2F("trackZeta", "MFT tracks; track.z; #eta ", 100, -100., 100., 102, -4.51, 4.51));

    // other options:
    // TH1F* t = new TH1F(); trZ.setObject(t); <- resets content!
    // TH1F t(); trZ.setObject(t) <- makes a copy
    // trZ.setObject({"Z","Z",100,-10.,10.}); <- creates new
  }

  void process(o2::aod::Collision const& collision, aod::Tracks const& tracks)
  {
    for (auto& track : tracks) {
      ZvtxEta->Fill(collision.posZ(), track.eta());
      trackZeta->Fill(track.z(), track.eta());
      trZ->Fill(track.z());
      phieta->Fill(track.phi(), track.eta());
    }
  }
};

// Iterate on muon using the collision iterator in the dq-analysis style
struct IterateMuons {

  // histogram defined with HistogramRegistry
  HistogramRegistry registryMuons{
    "registryMuons",
    {{"phiMuons", "phiMuons", {HistType::kTH1F, {{100, 0., 2. * M_PI}}}},
     {"etaptMuons", "etaptMuons", {HistType::kTH2F, {{102, -4.51, 4.51}, {100, 0.0, 5.0}}}},
     {"ptMuons", "ptMuons", {HistType::kTH1F, {{102, 0., 40.}}}}}};

  void process(aod::Collisions::iterator const& collision, aod::FwdTracks const& muons)
  {
    // LOGF(info, "Vertex = %f has %d muons", collision.posZ(), muons.size());
    for (auto& muon : muons) {
      // LOGF(info, "  pT = %.2f", muon.pt());
      registryMuons.get<TH1>(HIST("phiMuons"))->Fill(muon.phi());
      registryMuons.get<TH2>(HIST("etaptMuons"))->Fill(muon.pt(), muon.phi());
      registryMuons.get<TH1>(HIST("ptMuons"))->Fill(muon.pt());
    }
  }
};

struct HistRegistry {

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {{"phiC", "; #varphi", {HistType::kTH1F, {{100, -M_PI, M_PI}}}},
     {"etaptC", "etaptC", {HistType::kTH2F, {{102, -4.51, 4.51}, {100, 0.0, 5.0}}}},
     {"ptC", ";p_{T}", {HistType::kTH1F, {{102, -1., 20.}}}},
     {"phietaC", "MFT tracks; #varphi; #eta", {HistType::kTH2F, {{200, -M_PI, M_PI}, {200, -4.51, 4.51}}}},
     {"Nch", "; N_{ch}", {HistType::kTH1F, {{100, 0., 100.}}}},
     {"EventsNtrkZvtx", "MFT tracks; N_{trk}; z_{vtx}", {HistType::kTH2F, {{50, 0, 50}, {100, -20.0, 20.0}}}}}};

  void process(aod::Collisions::iterator const& collision, MyTracks const& tracks)
  {
    for (auto& track : tracks) {
      registry.get<TH1>(HIST("phiC"))->Fill(track.phi());
      registry.get<TH2>(HIST("etaptC"))->Fill(track.eta(), track.pt());
      registry.get<TH1>(HIST("ptC"))->Fill(track.pt());
      registry.fill(HIST("phietaC"), track.phi(), track.eta());
      // registry.fill(HIST("trZeta"), track.z(), track.eta());
      // registry.fill(HIST("trZphi"), track.z(), track.phi());
    }
    registry.get<TH2>(HIST("EventsNtrkZvtx"))->Fill(tracks.size(), collision.posZ());
    registry.get<TH1>(HIST("Nch"))->Fill(tracks.size());
  }
};

struct QvectorAnalysis {

  TRandom3* fRndm = new TRandom3(0);
  Configurable<std::string> fConfigEventCuts{"cfgEventCuts", "eventStandard", "Event selection"};
  Configurable<std::string> fConfigTrackCuts{"cfgTrackCuts", "", "Comma separated list of barrel track cuts"};
  Configurable<std::string> fConfigMuonCuts{"cfgMuonCuts", "", "Comma separated list of muon cuts"};
  Configurable<bool> fConfigQA{"cfgQA", false, "If true, fill QA histograms"};
  Configurable<int> nHarm{"nHarm", 2, "Number of harmonics"};
  // Configurable<int> fmaxHarmonic{"maxHarmonic", 3, "Maximum harmonic to be computed"};
  // Configurable<int> fmaxPower{"maxPower", 1, "Maximum power to be computed"};
  Configurable<bool> bUseWeights{"UseWeights", false, "If true, fill Q vectors with weights for phi and p_T"};
  Configurable<bool> bsubEvents{"subEvents", true, "If true, fill use sub-events methods with different detector gaps"};
  Configurable<float> fetaLimit{"etaLimit", 0.0, "Eta gap separation (e.g ITS=0.0, MFT=-3.05), only if subEvents=true"};

  Configurable<std::string> url{"ccdb-url", "http://ccdb-test.cern.ch:8080", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdb-path", "Users/lm", "base path to the ccdb object"};
  Configurable<long> nolaterthan{"ccdb-no-later-than", std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count(), "latest acceptable timestamp of creation for the object"};

  HistogramRegistry registryQ{
    "registryQ",
    {{"hmult", "; N_{ch}", {HistType::kTH1F, {{100, 0, 200}}}},
     {"hpT", "; p_{T}", {HistType::kTH1F, {{100, 0, 20}}}},
     {"hpT_0", "; p_{T}", {HistType::kTH1F, {{100, 0, 20}}}},
     {"hpT_4", "; p_{T}", {HistType::kTH1F, {{100, 0, 20}}}},
     {"htracketa", "; #eta", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"htracketa_0", "; #eta", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"htracketa_4", "; #eta", {HistType::kTH1F, {{100, -5., 5.}}}},
     {"h2QnX", "; Q_{n,x}", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"h2QnY", "; Q_{n,y}", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"h2Psin", "; raw #Psi_{n}", {HistType::kTH1F, {{100, -2., 2.}}}},
     {"Zvtx_h2QnX", ";  z_{vtx}; raw Q_{n,x}", {HistType::kTProfile, {{40, -20.0f, 20.0f}}}},
     {"Zvtx_h2QnY", ";  z_{vtx}; raw Q_{n,y}", {HistType::kTProfile, {{40, -20.0f, 20.0f}}}},
     {"Mult_h2QnX", "; N_{ch}; raw Q_{n,x}", {HistType::kTProfile, {{20, -1, 99}}}},
     {"Mult_h2QnY", "; N_{ch}; raw Q_{n,y}", {HistType::kTProfile, {{20, -1, 99}}}},
     {"Mult_resGap", "; N_{ch}; resGap", {HistType::kTProfile, {{20, -0.5f, 99.5f}}}},
     {"h2Vn", "; raw v_{n}", {HistType::kTH1F, {{50, -3., 3.}}}},
     {"Mult_Qn0Qn1", "; N_{ch}; raw Q_{n,A}Q_{n,B}^{*}", {HistType::kTProfile, {{100, -1, 99}}}},
     {"Mult_Qn0Qn2", "; N_{ch}; raw Q_{n,A}Q_{n,C}^{*}", {HistType::kTProfile, {{100, -1, 99}}}},
     {"Mult_Qn1Qn2", "; N_{ch}; raw Q_{n,B}Q_{n,C}^{*}", {HistType::kTProfile, {{100, -1, 99}}}},
     {"pT_h2VnSP", "; p_{T}; raw v_{n} {SP}", {HistType::kTProfile, {{20, -0.1, 9.9}}}},
     {"pT_h2VnEP", "; p_{T}; raw v_{n} {EP}", {HistType::kTProfile, {{20, -0.1, 9.9}}}}

    }};
  int nMult = 0;    // event multiplicity
  int nMultPos = 0; // event multiplicity positive eta gap
  int nMultNeg = 0; // event multiplicity negative eta gap

  double dPhi = 0., wPhi = 1., wPhiToPowerP = 1.;    // azimuthal angle and corresponding weight
  double vnrawSP = 0.0, vnrawEP = 0.0, resGap = 0.0; // flow coefficients and resolution

  void process(aod::Collisions::iterator const& collision, MyTracks const& tracks, MyTracksSelected const& muons)
  {

    nMult = tracks.size();
    // auto bc = collision.bc_as<aod::BCsWithTimestamps>();
    nMultPos = 0;
    nMultNeg = 0;

    registryQ.get<TH1>(HIST("hmult"))->Fill(nMult);

    for (auto& track : tracks) {
      dPhi = track.phi();
      // determine corresponding weight:
      if (bUseWeights) {
        wPhi = track.pt();
      }
      // Calculate Q-vector components:
      for (Int_t h = 0; h < maxHarmonic; h++) {
        for (Int_t p = 0; p < maxPower; p++) {
          // if(bUseWeights){wPhiToPowerP = pow(wPhi,p);}
          Qvector[h][p] += TComplex(wPhiToPowerP * TMath::Cos(h * dPhi), wPhiToPowerP * TMath::Sin(h * dPhi));
          if (bsubEvents) {
            if (track.eta() > fetaLimit) {
              QvectorPos[h][p] += TComplex(wPhiToPowerP * TMath::Cos(h * dPhi), wPhiToPowerP * TMath::Sin(h * dPhi));
              nMultPos++;
            }
            if (track.eta() < fetaLimit) {
              QvectorNeg[h][p] += TComplex(wPhiToPowerP * TMath::Cos(h * dPhi), wPhiToPowerP * TMath::Sin(h * dPhi));
              nMultNeg++;
            }
          }
          registryQ.get<TH1>(HIST("htracketa"))->Fill(track.eta());
          registryQ.get<TH1>(HIST("hpT"))->Fill(track.pt());
        }
      }
    } // loop over tracks

    if (nMult > 1.0) {
      // double normQ = TMath::Sqrt(Qvector[nHarm][0].Re() * Qvector[nHarm][0].Re() + Qvector[nHarm][0].Im() * Qvector[nHarm][0].Im());
      double normFactor = nMult; // or normQ;
      TComplex QvectorNormalized = TComplex(Qvector[nHarm][0].Re() / normFactor, Qvector[nHarm][0].Im() / normFactor);

      if (bsubEvents && (nMultPos>0 && nMultNeg>0) ) {
        // double normQPos = TMath::Sqrt(QvectorPos[nHarm][0].Re() * QvectorPos[nHarm][0].Re() + QvectorPos[nHarm][0].Im() * QvectorPos[nHarm][0].Im());
        // double normQNeg = TMath::Sqrt(QvectorNeg[nHarm][0].Re() * QvectorNeg[nHarm][0].Re() + QvectorNeg[nHarm][0].Im() * QvectorNeg[nHarm][0].Im());

        // double normFactorPos = nMultPos; // or normQ;
        // double normFactorNeg = nMultNeg; // or normQ;

        TComplex QvectorNormalizedPos = TComplex(QvectorPos[nHarm][0].Re() / normFactor, QvectorPos[nHarm][0].Im() / nMultPos);
        TComplex QvectorNormalizedNeg = TComplex(QvectorNeg[nHarm][0].Re() / normFactor, QvectorNeg[nHarm][0].Im() / nMultNeg);

        resGap = (QvectorPos[nHarm][0].Re() * QvectorNeg[nHarm][0].Re() + QvectorPos[nHarm][0].Im() * QvectorNeg[nHarm][0].Im()) / (nMultPos * nMultNeg);
      }

      double rawPsin = (1.0 / nHarm) * TMath::ATan2(Qvector[nHarm][0].Re() / normFactor, Qvector[nHarm][0].Im() / normFactor);

      // double refFlowRn = TMath::Sqrt( TMath::Abs((QvectorNormalized.Re()*QvectorNormalizedNeg.Im() + QvectorNormalized.Im()*QvectorNormalizedNeg.Re() )*( QvectorNormalized.Re()*QvectorNormalizedPos.Im() + QvectorNormalized.Im()*QvectorNormalizedPos.Re() )/(QvectorNormalizedPos.Re()*QvectorNormalizedNeg.Im() + QvectorNormalizedPos.Im()*QvectorNormalizedNeg.Re() ) ) );
      // double normQ0 =  QvectorNormalized.Re()*QvectorNormalized.Re() + QvectorNormalized.Im()*QvectorNormalized.Im()  ;
      // double normQ1 =  QvectorNormalizedPos.Re()*QvectorNormalizedPos.Re() + QvectorNormalizedPos.Im()*QvectorNormalizedPos.Im()  ;
      // double normQ2 =  QvectorNormalizedNeg.Re()*QvectorNormalizedNeg.Re() + QvectorNormalizedNeg.Im()*QvectorNormalizedNeg.Im()  ;

      // registryQ.get<TH1>(HIST("hEvperRun"))->Fill(bc.runNumber());

      registryQ.get<TH1>(HIST("h2QnX"))->Fill(QvectorNormalized.Re());
      registryQ.get<TH1>(HIST("h2QnY"))->Fill(QvectorNormalized.Im());
      registryQ.get<TH1>(HIST("h2Psin"))->Fill(rawPsin);

      registryQ.get<TProfile>(HIST("Zvtx_h2QnX"))->Fill(collision.posZ(), QvectorNormalized.Re());
      registryQ.get<TProfile>(HIST("Zvtx_h2QnY"))->Fill(collision.posZ(), QvectorNormalized.Im());
      registryQ.get<TProfile>(HIST("Mult_h2QnX"))->Fill(nMult, QvectorNormalized.Re());
      registryQ.get<TProfile>(HIST("Mult_h2QnY"))->Fill(nMult, QvectorNormalized.Im());

      registryQ.get<TProfile>(HIST("Mult_resGap"))->Fill(nMult, resGap);

      //            for (auto& track : tracks) {
      //                dPhi = track.phi();
      //
      //                //LOGF(info, "     nMult = %.2f", nMult );
      //                // Calculate vn = uQ products:
      //                rawSP =  (TMath::Cos(nHarm*dPhi)*QvectorNormalized.Re() + TMath::Sin(nHarm*dPhi)*QvectorNormalized.Im()) ;
      //                rawEP =  TMath::Cos(nHarm*(dPhi-rawPsin) ) ;
      //
      //                if(rawSP) registryQ.get<TH1>(HIST("h2Vn"))->Fill(rawSP);
      //                if(normQ0 && normQ1) registryQ.get<TProfile>(HIST("Mult_Qn0Qn1"))->Fill(nMult, normQ0*normQ1 );
      //                if(normQ0 && normQ2) registryQ.get<TProfile>(HIST("Mult_Qn0Qn2"))->Fill(nMult, normQ0*normQ2 );
      //                if(normQ1 && normQ2) registryQ.get<TProfile>(HIST("Mult_Qn1Qn2"))->Fill(nMult, normQ1*normQ2 );
      //                if( rawSP ) registryQ.get<TProfile>(HIST("Mult_h2VnSP"))->Fill(nMult, rawSP );
      //                if( rawEP ) registryQ.get<TProfile>(HIST("Mult_h2VnEP"))->Fill(nMult, rawEP );
      //
      //            } // loop over tracks

      for (auto& muon : muons) {
        dPhi = muon.phi();
        vnrawSP = (TMath::Cos(nHarm * dPhi) * QvectorNormalized.Re() + TMath::Sin(nHarm * dPhi) * QvectorNormalized.Im()) / nMult;
        vnrawEP = TMath::Cos(nHarm * (dPhi - rawPsin));

        if (muon.trackType() == 0) {
          registryQ.get<TH1>(HIST("hpT_0"))->Fill(muon.pt());
          registryQ.get<TH1>(HIST("htracketa_0"))->Fill(muon.pt());
        }
        if (muon.trackType() == 4) {
          registryQ.get<TH1>(HIST("hpT_4"))->Fill(muon.pt());
          registryQ.get<TH1>(HIST("htracketa_4"))->Fill(muon.pt());
        }

        if (vnrawSP)
          registryQ.get<TProfile>(HIST("pT_h2VnSP"))->Fill(muon.pt(), vnrawSP);
        if (vnrawEP)
          registryQ.get<TProfile>(HIST("pT_h2VnEP"))->Fill(muon.pt(), vnrawEP);

      } // loop over tracks
    }

    ResetQvector();
  }

  void ResetQvector()
  {
    // Reset all Q-vector components to zero before starting a new event.
    for (Int_t h = 0; h < maxHarmonic; h++) {
      for (Int_t p = 0; p < maxPower; p++) {
        Qvector[h][p] = TComplex(0., 0.);
        QvectorPos[h][p] = TComplex(0., 0.);
        QvectorNeg[h][p] = TComplex(0., 0.);
      }
    }
  }
};

struct LoopOverMcMatched {
  OutputObj<TH1F> etaDiff{TH1F("etaDiff", ";eta_{MC} - eta_{Rec}", 100, -2, 2)};
  void process(aod::McCollision const& mcCollision, soa::SmallGroups<soa::Join<aod::McCollisionLabels, aod::Collisions>> const& collisions,
               soa::Join<aod::Tracks, aod::McTrackLabels> const& tracks, aod::McParticles_001 const& mcParticles)
  {
    // access MC truth information with mcCollision() and mcParticle() methods
    // LOGF(info, "MC collision at vtx-z = %f with %d mc particles and %d reconstructed collisions", mcCollision.posZ(), mcParticles.size(), collisions.size());
    for (auto& collision : collisions) {
      // LOGF(info, "  Reconstructed collision at vtx-z = %f", collision.posZ());

      // NOTE this will be replaced by a improved grouping in the future
      auto groupedTracks = tracks.sliceBy(aod::track::collisionId, collision.globalIndex());
      // LOGF(info, "  which has %d tracks", groupedTracks.size());
      for (auto& track : groupedTracks) {
        if (!track.has_mcParticle()) {
          // LOGF(warning, "No MC particle for track, skip...");
          continue;
        }
        etaDiff->Fill(track.mcParticle().eta() - track.eta());
      }
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    // adaptAnalysisTask<RootHistograms>(cfgc),
    // adaptAnalysisTask<OutputObjects>(cfgc),
    // adaptAnalysisTask<OutputObjSet>(cfgc),
    // adaptAnalysisTask<HistRegistry>(cfgc),
    // adaptAnalysisTask<LoopOverMcMatched>(cfgc),
    // adaptAnalysisTask<MultiplicityEventTrackSelectionMFT>(cfgc),
    adaptAnalysisTask<IterateMuons>(cfgc),
    adaptAnalysisTask<QvectorAnalysis>(cfgc)};
}
