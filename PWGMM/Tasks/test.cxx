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
/// \brief Both tasks, ATask and BTask create two histograms. But whereas in
///        the first case (ATask) the histograms are not saved to file, this
///        happens automatically if OutputObj<TH1F> is used to create a
///        histogram. By default the histogram is saved to file
///        AnalysisResults.root. HistogramRegistry is yet an other possibility
///        to deal with histograms. See tutorial example histogramRegistery.cxx
///        for details.
/// \author
/// \since

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

struct RootHistograms {

  // normal creation of a histogram
  TH1F* phiHA = new TH1F("phiA", "phiA", 100, 0., 2. * M_PI);
  TH1F* etaHA = new TH1F("etaA", "etaA", 102, -4.01, 4.01);
    TH2F* trXYHA = new TH2F("trXYHA", "trXYHA", 100, -20, 20, 100, -20, 20);
    TH2F* trXY0HA = new TH2F("trXY0HA", "trXY0HA", 100, -20, 20, 100, -20, 20);
    TH2F* trXY4HA = new TH2F("trXY4HA", "trXY4HA", 100, -20, 20, 100, -20, 20);

  void process(aod::Tracks const& tracks)
  {
    for (auto& track : tracks) {
      phiHA->Fill(track.phi());
      etaHA->Fill(track.eta());
        trXYHA->Fill(track.x(),track.y());
        if(track.z() > 46.) {  trXY0HA->Fill(track.x(),track.y()); }
        if(track.z() < 49.) { trXY4HA->Fill(track.x(),track.y()); }
    }
  }
};

struct OutputObjects {

  // histogram created with OutputObj<TH1F>
  OutputObj<TH1F> phiB{TH1F("phiB", "phiB", 100, 0., 2. * M_PI), OutputObjHandlingPolicy::QAObject};
  OutputObj<TH2F> etaptB{TH2F("etaptB", "etaptB", 102, -2.01, 2.01, 100, 0.0, 5.0), OutputObjHandlingPolicy::AnalysisObject};
    OutputObj<TH2F> trXY{TH2F("trXY", "trXY", 100, -20, 20, 100, -20, 20), OutputObjHandlingPolicy::QAObject};
    OutputObj<TH2F> trXY0{TH2F("trXY0", "trXY0",  200, -100, 100, 100, -20, 20), OutputObjHandlingPolicy::QAObject};
    OutputObj<TH2F> trXY4{TH2F("trXY4", "trXY4",  110, -55, 55, 100, -20, 20), OutputObjHandlingPolicy::QAObject};

  void process(aod::Tracks const& tracks)
  {
    for (auto& track : tracks) {
      phiB->Fill(track.phi());
      etaptB->Fill(track.eta(), track.pt());
        trXY->Fill(track.x(),track.y());
        trXY0->Fill(track.z(),track.x());
        trXY4->Fill(track.z(),track.y());
    }
  }
};

struct OutputObjSet {
  // incomplete definition of an OutputObj
  OutputObj<TH1F> trZ{"trZ", OutputObjHandlingPolicy::QAObject};

  Filter ptfilter = aod::track::pt > 0.5f;

  void init(InitContext const&)
  {
    // complete the definition of the OutputObj
    trZ.setObject(new TH1F("Z", "Z", 100, -20., 20.));

    // other options:
    // TH1F* t = new TH1F(); trZ.setObject(t); <- resets content!
    // TH1F t(); trZ.setObject(t) <- makes a copy
    // trZ.setObject({"Z","Z",100,-10.,10.}); <- creates new
  }

  void process(soa::Filtered<aod::Tracks> const& tracks)
  {
    for (auto& track : tracks) {
      trZ->Fill(track.z());
      
    }
  }
};

struct HistRegistry {

  // histogram defined with HistogramRegistry
  HistogramRegistry registry{
    "registry",
    {{"phiC", "phiC", {HistType::kTH1F, {{100, 0., 2. * M_PI}}}},
     {"etaptC", "etaptC", {HistType::kTH2F, {{102, -4.51, 4.51}, {100, 0.0, 5.0}}}},
    {"ptC", "ptC", {HistType::kTH1F, {{102, 0., 20.}}}}  }};

  void process(aod::Tracks const& tracks)
  {
    for (auto& track : tracks) {
      registry.get<TH1>(HIST("phiC"))->Fill(track.phi());
      registry.get<TH2>(HIST("etaptC"))->Fill(track.eta(), track.pt());
        registry.get<TH1>(HIST("ptC"))->Fill(track.pt());

    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<RootHistograms>(cfgc),
    adaptAnalysisTask<OutputObjects>(cfgc),
    adaptAnalysisTask<OutputObjSet>(cfgc),
    adaptAnalysisTask<HistRegistry>(cfgc),
  };
}
