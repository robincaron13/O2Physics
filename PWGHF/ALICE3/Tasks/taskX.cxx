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

/// \file taskX.cxx
/// \brief X(3872) analysis task
///
/// \author Gian Michele Innocenti <gian.michele.innocenti@cern.ch>, CERN
/// \author Rik Spijkers <r.spijkers@students.uu.nl>, Utrecht University
/// \author Luca Micheletti <luca.micheletti@to.infn.it>, INFN

#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"

#include "PWGHF/Core/SelectorCuts.h"
#include "PWGHF/DataModel/CandidateReconstructionTables.h"
#include "PWGHF/DataModel/CandidateSelectionTables.h"

using namespace o2;
using namespace o2::aod;
using namespace o2::analysis;
using namespace o2::analysis::hf_cuts_x_to_jpsi_pi_pi;
using namespace o2::framework;
using namespace o2::aod::hf_cand_x;
using namespace o2::framework::expressions;
using namespace o2::aod::hf_cand_2prong;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  ConfigParamSpec optionDoMC{"doMC", VariantType::Bool, true, {"Fill MC histograms."}};
  workflowOptions.push_back(optionDoMC);
}

#include "Framework/runDataProcessing.h"

/// X(3872) analysis task
struct HfTaskX {
  Configurable<int> selectionFlagX{"selectionFlagX", 1, "Selection Flag for X"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<bool> modeXToJpsiToMuMuPiPi{"modeXToJpsiToMuMuPiPi", false, "Perform Jpsi to mu+mu- analysis"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_x_to_jpsi_pi_pi::vecBinsPt}, "pT bin limits"};

  Filter filterSelectCandidates = (aod::hf_sel_candidate_x::isSelXToJpsiToEEPiPi >= selectionFlagX || aod::hf_sel_candidate_x::isSelXToJpsiToMuMuPiPi >= selectionFlagX);

  HistogramRegistry registry{
    "registry",
    {{"hPtProng0", "3-prong candidates;prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtProng1", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtProng2", "3-prong candidates;prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{100, 0., 10.}}}},
     {"hPtCand", "3-prong candidates;candidate #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{150, 0., 15.}}}}}};

  void init(o2::framework::InitContext&)
  {
    registry.add("hMass", "3-prong candidates;inv. mass (J/#psi #pi+ #pi-) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 3., 4.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLength", "3-prong candidates;decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLengthxy", "3-prong candidates;decay length xy (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0", "3-prong candidates;prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong2", "3-prong candidates;prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPA", "3-prong candidates;cosine of pointing angle;entries", {HistType::kTH2F, {{500, 0.9, 1.0}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEta", "3-prong candidates;candidate #it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hImpParErr", "3-prong candidates;impact parameter error (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenErr", "3-prong candidates;decay length error (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDecLenXYErr", "3-prong candidates;decay length xy error (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  void process(soa::Filtered<soa::Join<aod::HfCandX, aod::HfSelXToJpsiPiPi>> const& candidates)
  {
    int decayMode = modeXToJpsiToMuMuPiPi ? hf_cand_x::DecayType::XToJpsiToMuMuPiPi : hf_cand_x::DecayType::XToJpsiToEEPiPi;
    for (const auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << decayMode)) {
        continue;
      }
      if (yCandMax >= 0. && std::abs(yX(candidate)) > yCandMax) {
        continue;
      }

      registry.fill(HIST("hMass"), invMassXToJpsiPiPi(candidate), candidate.pt());
      registry.fill(HIST("hPtCand"), candidate.pt());
      registry.fill(HIST("hPtProng0"), candidate.ptProng0());
      registry.fill(HIST("hPtProng1"), candidate.ptProng1());
      registry.fill(HIST("hPtProng2"), candidate.ptProng2());
      registry.fill(HIST("hDecLength"), candidate.decayLength(), candidate.pt());
      registry.fill(HIST("hDecLengthxy"), candidate.decayLengthXY(), candidate.pt());
      registry.fill(HIST("hd0Prong0"), candidate.impactParameter0(), candidate.pt());
      registry.fill(HIST("hd0Prong1"), candidate.impactParameter1(), candidate.pt());
      registry.fill(HIST("hd0Prong2"), candidate.impactParameter2(), candidate.pt());
      registry.fill(HIST("hCPA"), candidate.cpa(), candidate.pt());
      registry.fill(HIST("hEta"), candidate.eta(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter0(), candidate.pt());
      registry.fill(HIST("hImpParErr"), candidate.errorImpactParameter1(), candidate.pt());
      registry.fill(HIST("hDecLenErr"), candidate.errorDecayLength(), candidate.pt());
      registry.fill(HIST("hDecLenXYErr"), candidate.errorDecayLengthXY(), candidate.pt());
    } // candidate loop
  }   // process
};    // struct

struct HfTaskXMc {
  Configurable<int> selectionFlagX{"selectionFlagX", 1, "Selection Flag for X"};
  Configurable<double> yCandMax{"yCandMax", -1., "max. cand. rapidity"};
  Configurable<bool> modeXToJpsiToMuMuPiPi{"modeXToJpsiToMuMuPiPi", false, "Perform Jpsi to mu+mu- analysis"};
  Configurable<std::vector<double>> binsPt{"binsPt", std::vector<double>{hf_cuts_x_to_jpsi_pi_pi::vecBinsPt}, "pT bin limits"};

  Filter filterSelectCandidates = (aod::hf_sel_candidate_x::isSelXToJpsiToEEPiPi >= selectionFlagX || aod::hf_sel_candidate_x::isSelXToJpsiToMuMuPiPi >= selectionFlagX);

  HistogramRegistry registry{
    "registry",
    {{"hPtRecSig", "3-prong candidates (rec. matched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{150, 0., 15.}}}},
     {"hPtRecBg", "3-prong candidates (rec. unmatched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{150, 0., 15.}}}},
     {"hPtGen", "3-prong candidates (gen. matched);#it{p}_{T} (GeV/#it{c});entries", {HistType::kTH1F, {{150, 0., 15.}}}},
     {"hPtGenSig", "3-prong candidates (rec. matched);#it{p}_{T}^{gen.} (GeV/#it{c});entries", {HistType::kTH1F, {{150, 0., 15.}}}}}};

  void init(o2::framework::InitContext&)
  {
    registry.add("hCPARecSig", "3-prong candidates (rec. matched);cosine of pointing angle;entries", {HistType::kTH2F, {{500, 0.9, 1.0}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCPARecBg", "3-prong candidates (rec. unmatched);cosine of pointing angle;entries", {HistType::kTH2F, {{500, 0.9, 1.0}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecSig", "3-prong candidates (rec. matched);#it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaRecBg", "3-prong candidates (rec. unmatched);#it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hEtaGen", "3-prong candidates (gen. matched);#it{#eta};entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hPtProng0RecSig", "3-prong candidates (rec. matched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecSig", "3-prong candidates (rec. matched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng2RecSig", "3-prong candidates (rec. matched);prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng0RecBg", "3-prong candidates (rec. unmatched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng1RecBg", "3-prong candidates (rec. unmatched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtProng2RecBg", "3-prong candidates (rec. unmatched);prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenProng0", "3-prong candidates (gen. matched);prong 0 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenProng1", "3-prong candidates (gen. matched);prong 1 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hPtGenProng2", "3-prong candidates (gen. matched);prong 2 #it{p}_{T} (GeV/#it{c});entries", {HistType::kTH2F, {{100, 0., 10.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hMassRecSig", "3-prong candidates (rec. matched);inv. mass (J/#psi #pi+ #pi-) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 3., 4.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hMassRecBg", "3-prong candidates (rec. unmatched);inv. mass (J/#psi #pi+ #pi-) (GeV/#it{c}^{2});entries", {HistType::kTH2F, {{200, 3., 4.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecSig", "3-prong candidates (rec. matched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecSig", "3-prong candidates (rec. matched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong2RecSig", "3-prong candidates (rec. matched);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong0RecBg", "3-prong candidates (rec. unmatched);prong 0 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong1RecBg", "3-prong candidates (rec. unmatched);prong 1 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hd0Prong2RecBg", "3-prong candidates (rec. unmatched);prong 2 DCAxy to prim. vertex (cm);entries", {HistType::kTH2F, {{400, -0.002, 0.002}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDeclengthRecSig", "3-prong candidates (rec. matched);decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hDeclengthRecBg", "3-prong candidates (rec. unmatched);decay length (cm);entries", {HistType::kTH2F, {{100, 0., 0.01}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});

    registry.add("hChi2PCASig", "3-prong candidates (rec. matched);chi2 PCA (cm);entries", {HistType::kTH2F, {{1000, 0., 0.0001}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hChi2PCABg", "3-prong candidates (rec. unmatched);chi2 PCA (cm);entries", {HistType::kTH2F, {{1000, 0., 0.0001}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtSig", "3-prong candidates (rec. matched);proper lifetime X(3872) * #it{c} (cm);entries", {HistType::kTH2F, {{400, 0., 0.001}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hCtBg", "3-prong candidates (rec. unmatched);proper lifetime X(3872) * #it{c} (cm);entries", {HistType::kTH2F, {{400, 0., 0.001}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYSig", "3-prong candidates (rec. matched);candidate rapidity;entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
    registry.add("hYBg", "3-prong candidates (rec. unmatched);candidate rapidity;entries", {HistType::kTH2F, {{100, -2., 2.}, {(std::vector<double>)binsPt, "#it{p}_{T} (GeV/#it{c})"}}});
  }

  void process(soa::Filtered<soa::Join<aod::HfCandX, aod::HfSelXToJpsiPiPi, aod::HfCandXMcRec>> const& candidates,
               soa::Join<aod::McParticles, aod::HfCandXMcGen> const& particlesMC,
               aod::TracksWMc const& tracks)
  {
    // MC rec.
    // Printf("MC Candidates: %d", candidates.size());
    int decayMode = modeXToJpsiToMuMuPiPi ? hf_cand_x::DecayType::XToJpsiToMuMuPiPi : hf_cand_x::DecayType::XToJpsiToEEPiPi;
    for (const auto& candidate : candidates) {
      if (!(candidate.hfflag() & 1 << decayMode)) {
        continue;
      }
      if (yCandMax >= 0. && std::abs(yX(candidate)) > yCandMax) {
        continue;
      }
      if (candidate.flagMcMatchRec() == 1 << decayMode) {
        auto indexMother = RecoDecay::getMother(particlesMC, candidate.prong1_as<aod::TracksWMc>().mcParticle_as<soa::Join<aod::McParticles, aod::HfCandXMcGen>>(), 9920443, true);
        auto particleMother = particlesMC.rawIteratorAt(indexMother);
        registry.fill(HIST("hPtGenSig"), particleMother.pt());
        registry.fill(HIST("hPtRecSig"), candidate.pt());
        registry.fill(HIST("hCPARecSig"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecSig"), candidate.eta(), candidate.pt());

        registry.fill(HIST("hDeclengthRecSig"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hMassRecSig"), invMassXToJpsiPiPi(candidate), candidate.pt());
        registry.fill(HIST("hd0Prong0RecSig"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecSig"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hd0Prong2RecSig"), candidate.impactParameter2(), candidate.pt());
        registry.fill(HIST("hPtProng0RecSig"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecSig"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("hPtProng2RecSig"), candidate.ptProng2(), candidate.pt());
        registry.fill(HIST("hChi2PCASig"), candidate.chi2PCA(), candidate.pt());
        registry.fill(HIST("hCtSig"), ctX(candidate), candidate.pt());
        registry.fill(HIST("hYSig"), yX(candidate), candidate.pt());
      } else {
        registry.fill(HIST("hPtRecBg"), candidate.pt());
        registry.fill(HIST("hCPARecBg"), candidate.cpa(), candidate.pt());
        registry.fill(HIST("hEtaRecBg"), candidate.eta(), candidate.pt());

        registry.fill(HIST("hDeclengthRecBg"), candidate.decayLength(), candidate.pt());
        registry.fill(HIST("hMassRecBg"), invMassXToJpsiPiPi(candidate), candidate.pt());
        registry.fill(HIST("hd0Prong0RecBg"), candidate.impactParameter0(), candidate.pt());
        registry.fill(HIST("hd0Prong1RecBg"), candidate.impactParameter1(), candidate.pt());
        registry.fill(HIST("hd0Prong2RecBg"), candidate.impactParameter2(), candidate.pt());
        registry.fill(HIST("hPtProng0RecBg"), candidate.ptProng0(), candidate.pt());
        registry.fill(HIST("hPtProng1RecBg"), candidate.ptProng1(), candidate.pt());
        registry.fill(HIST("hPtProng2RecBg"), candidate.ptProng2(), candidate.pt());
        registry.fill(HIST("hChi2PCABg"), candidate.chi2PCA(), candidate.pt());
        registry.fill(HIST("hCtBg"), ctX(candidate), candidate.pt());
        registry.fill(HIST("hYBg"), yX(candidate), candidate.pt());
      }
    } // rec
    // MC gen.
    // Printf("MC Particles: %d", particlesMC.size());
    for (const auto& particle : particlesMC) {
      if (particle.flagMcMatchGen() == 1 << decayMode) {
        // TODO: add X(3872) mass such that we can use the getMassPDG function instead of hardcoded mass
        if (yCandMax >= 0. && std::abs(RecoDecay::y(std::array{particle.px(), particle.py(), particle.pz()}, 3.87168)) > yCandMax) {
          // Printf("MC Gen.: Y rejection: %g", RecoDecay::Y(std::array{particle.px(), particle.py(), particle.pz()}, 3.87168));
          continue;
        }
        registry.fill(HIST("hPtGen"), particle.pt());
        registry.fill(HIST("hEtaGen"), particle.eta(), particle.pt());

        // properties of gen matched X(3872), to get a first look at some cuts
        float ptProngs[3];
        int counter = 0;
        for (const auto& daugh : particle.daughters_as<aod::McParticles>()) {
          ptProngs[counter] = daugh.pt();
          counter++;
        }
        registry.fill(HIST("hPtGenProng0"), ptProngs[0], particle.pt());
        registry.fill(HIST("hPtGenProng1"), ptProngs[1], particle.pt());
        registry.fill(HIST("hPtGenProng2"), ptProngs[2], particle.pt());
      }
    } // gen
  }   // process
};    // struct

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  WorkflowSpec workflow{
    adaptAnalysisTask<HfTaskX>(cfgc)};
  const bool doMC = cfgc.options().get<bool>("doMC");
  if (doMC) {
    workflow.push_back(adaptAnalysisTask<HfTaskXMc>(cfgc));
  }
  return workflow;
}
