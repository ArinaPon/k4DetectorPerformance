/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include "TrackingValidationHelpers.h"
#include "TrackingValidationPlots.h"

// k4FWCore
#include "k4FWCore/Consumer.h"

// Gaudi
#include "Gaudi/Property.h"
#include "GaudiKernel/MsgStream.h"

// EDM4hep
#include "edm4hep/MCParticleCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TCanvas.h"


// STL
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>


/** @struct TrackingValidation
 *
 *  Gaudi Consumer that validates the performance of track finding and track fitting
 *  by comparing reconstructed tracks with Monte Carlo truth information and,
 *  optionally, with perfectly associated fitted tracks.
 *
 *  The consumer writes several ROOT TTrees containing finder-level associations,
 *  fitter residuals with respect to MC truth, and fitter residuals with respect
 *  to perfectly fitted reference tracks. In addition, summary performance plots
 *  are produced in finalize() and written to the same ROOT file.
 *
 *  The supported validation modes are:
 *    - full pipeline validation,
 *    - finder-only validation,
 *    - fitter-only validation.
 *
 *  input:
 *    - MC particle collection : edm4hep::MCParticleCollection
 *    - planar digi-to-sim link collections : std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>
 *    - drift-chamber digi-to-sim link collections : std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>
 *    - finder track collection : edm4hep::TrackCollection
 *    - fitted track collection : edm4hep::TrackCollection
 *    - optional perfect fitted-track collections : std::vector<const edm4hep::TrackCollection*>
 *
 *  output:
 *    - ROOT file containing validation TTrees
 *    - summary performance plots written to the same ROOT file
 *
 *  @author Arina Ponomareva
 *  @date   2026-03
 *
 */


// ---------- CONSUMER ----------
struct TrackingValidation final
    : k4FWCore::Consumer<void(
          const edm4hep::MCParticleCollection&,
          const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&,  // planar links
          const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&,  // optional DCH links
          const edm4hep::TrackCollection&,                                            // finder tracks
          const edm4hep::TrackCollection&,                                            // fitted tracks (reco)
          const std::vector<const edm4hep::TrackCollection*>&                         // optional perfect fitted tracks
          )> {

  TrackingValidation(const std::string& name, ISvcLocator* svcLoc)
      : Consumer(
            name, svcLoc,
            {
                KeyValues("MCParticles", {"MCParticles"}),

                KeyValues("PlanarLinks",
                          {"SiWrBSimDigiLinks", "SiWrDSimDigiLinks", "VTXBSimDigiLinks", "VTXDSimDigiLinks"}),

                
                KeyValues("DCHLinks", {"DCH_DigiSimAssociationCollection"}),

                KeyValues("FinderTracks", {"GGTFTracks"}),
                KeyValues("FittedTracks", {"FittedTracks"}),

                
                KeyValues("PerfectFittedTracks", {"PerfectFittedTracks"}),
            }) {}

  StatusCode initialize() override {
    info() << "Initializing TrackingValidationConsumer" << endmsg;

    m_outFile = std::make_unique<TFile>(m_outputFile.value().c_str(), "RECREATE");
    if (!m_outFile || m_outFile->IsZombie()) {
      error() << "Cannot open output file: " << m_outputFile.value() << endmsg;
      return StatusCode::FAILURE;
    }

    bookAssocTree(m_finder_p2t, "finder_particle_to_tracks");
    bookAssocTree(m_finder_t2p, "finder_track_to_particles");
    bookAssocTree(m_perf_p2t, "perfect_particle_to_tracks");
    bookAssocTree(m_perf_t2p, "perfect_track_to_particles");

    bookFitterTree(m_fit_vs_mc, "fitter_vs_mc");
    bookFitterTree(m_fit_vs_perfect, "fitter_vs_perfect");

    return StatusCode::SUCCESS;
  }

  void operator()(const edm4hep::MCParticleCollection& mcParts,
                  const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& planarLinksVec,
                  const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& dchLinksVec,
                  const edm4hep::TrackCollection& finderTracks,
                  const edm4hep::TrackCollection& fittedTracks,
                  const std::vector<const edm4hep::TrackCollection*>& perfectFittedTracksVec) const override {

    const int event = m_evt++;
    const int mode = m_mode.value();  // 0 full, 1 finder-only, 2 fitter-only

    // ---------- Build truth maps: hit -> particle, particle -> hits ----------
    std::unordered_map<int, std::vector<uint64_t>> hitsPerParticle;
    hitsPerParticle.reserve(mcParts.size());

    std::unordered_map<uint64_t, int> hitToParticle;
    hitToParticle.reserve(200000);

    // planar links
    for (const auto* links : planarLinksVec) {
      if (!links) continue;
      for (const auto& link : *links) {
        const auto digi = link.getFrom();
        const auto sim  = link.getTo();
        const auto mc   = sim.getParticle();
        if (!digi.isAvailable() || !mc.isAvailable()) continue;

        const int pid = mc.getObjectID().index;
        const uint64_t key = TrackingValidationHelpers::oidKey(digi.getObjectID());
        hitsPerParticle[pid].push_back(key);
        hitToParticle[key] = pid;
      }
    }

    // DCH links
    for (const auto* links : dchLinksVec) {
      if (!links) continue;
      for (const auto& link : *links) {
        const auto digi = link.getFrom();
        const auto sim  = link.getTo();
        const auto mc   = sim.getParticle();
        if (!digi.isAvailable() || !mc.isAvailable()) continue;

        const int pid = mc.getObjectID().index;
        const uint64_t key = TrackingValidationHelpers::oidKey(digi.getObjectID());
        hitsPerParticle[pid].push_back(key);
        hitToParticle[key] = pid;
      }
    }

    // ---------- Finder & Perfect association trees ----------
    if (mode == 0 || mode == 1) {
      fillPerfectAssoc(event, mcParts, hitsPerParticle);
      fillFinderAssoc(event, mcParts, finderTracks, hitToParticle, hitsPerParticle);
    }

    // ---------- Build pid -> best perfect-fitted AtIP state ----------
    struct StateWithNHits {
      edm4hep::TrackState st;
      int nHits = 0;
    };

    std::unordered_map<int, StateWithNHits> perfectAtIPByPid;

    const bool wantPerfect = m_doPerfectFit.value();
    const bool havePerfectCollections = !perfectFittedTracksVec.empty();
    const bool doPerfect = wantPerfect && havePerfectCollections;

    if (wantPerfect && !havePerfectCollections && !m_warnedMissingPerfectInput) {
      warning() << "DoPerfectFit=true but no PerfectFittedTracks collection was provided. "
                << "fitter_vs_perfect will be filled with NaNs/empty content." << endmsg;
      m_warnedMissingPerfectInput = true;
    }

    if (doPerfect) {
      size_t nPerfectTracks = 0;
      for (const auto* coll : perfectFittedTracksVec) {
        if (!coll) continue;
        nPerfectTracks += coll->size();
      }
      perfectAtIPByPid.reserve(nPerfectTracks);

      for (const auto* coll : perfectFittedTracksVec) {
        if (!coll) continue;

        for (const auto& trk : *coll) {
          edm4hep::TrackState st;
          if (!TrackingValidationHelpers::getAtIPState(trk, st)) continue;

          const int pid = majorityParticleForTrack(trk, hitToParticle);
          if (pid < 0 || pid >= (int)mcParts.size()) continue;

          const int nHits = (int)trk.getTrackerHits().size();
          auto it = perfectAtIPByPid.find(pid);
          if (it == perfectAtIPByPid.end() || nHits > it->second.nHits) {
            perfectAtIPByPid[pid] = StateWithNHits{st, nHits};
          }
        }
      }
    }

    // ---------- Fitter trees ----------
    if (mode == 0 || mode == 2) {
      fillFitterTrees(event, mcParts, fittedTracks, hitToParticle, perfectAtIPByPid, doPerfect);
    }
  }

  StatusCode finalize() override {
    info() << "Finalizing TrackingValidation, wrote " << m_evt << " events" << endmsg;

    if (m_outFile) {
      m_outFile->cd();

      //write trees
      if (m_finder_p2t.tree) m_finder_p2t.tree->Write();
      if (m_finder_t2p.tree) m_finder_t2p.tree->Write();
      if (m_perf_p2t.tree) m_perf_p2t.tree->Write();
      if (m_perf_t2p.tree) m_perf_t2p.tree->Write();

      if (m_fit_vs_mc.tree) m_fit_vs_mc.tree->Write();
      if (m_fit_vs_perfect.tree) m_fit_vs_perfect.tree->Write();

      //fitter summary plots
      // d0 resolution vs momentum from fitter_vs_mc
      TGraphErrors* g_d0_vs_p = TrackingValidationPlots::makeD0ResolutionVsMomentum(m_fit_vs_mc.tree,
                                                         "g_d0_resolution_vs_p",
                                                         0.1, 100.0, 0.15);
      if (g_d0_vs_p) {
        TCanvas* c_d0_vs_p = TrackingValidationPlots::drawD0ResolutionCanvas(g_d0_vs_p,
                                                  "c_d0_resolution_vs_p",
                                                  0.1, 100.0);
        g_d0_vs_p->Write();
        if (c_d0_vs_p) c_d0_vs_p->Write();
      }
      // p resolution vs momentum
      TGraphErrors* g_p_vs_p = TrackingValidationPlots::makeMomentumResolutionVsMomentum(m_fit_vs_mc.tree,
                                                          "g_p_resolution_vs_p",
                                                          0.1, 100.0, 0.15);
      if (g_p_vs_p) {
        TCanvas* c_p_vs_p = TrackingValidationPlots::drawResolutionCanvas(g_p_vs_p,
                                           "c_p_resolution_vs_p",
                                           "momentum resolution vs momentum;p_{ref} [GeV];#sigma((p_{reco}-p_{ref})/p_{ref})",
                                           0.1, 100.0);
        g_p_vs_p->Write();
        if (c_p_vs_p) c_p_vs_p->Write();
      }

      // pT resolution vs momentum
      TGraphErrors* g_pt_vs_p = TrackingValidationPlots::makePtResolutionVsMomentum(m_fit_vs_mc.tree,
                                                     "g_pt_resolution_vs_p",
                                                     0.1, 100.0, 0.15);
      if (g_pt_vs_p) {
        TCanvas* c_pt_vs_p = TrackingValidationPlots::drawResolutionCanvas(g_pt_vs_p,
                                            "c_pt_resolution_vs_p",
                                            "pT resolution vs momentum;p_{ref} [GeV];#sigma((pT_{reco}-pT_{ref})/pT_{ref})",
                                            0.1, 100.0);
        g_pt_vs_p->Write();
        if (c_pt_vs_p) c_pt_vs_p->Write();
      }

      // finder summary plot
      TGraphErrors* g_eff_vs_p = TrackingValidationPlots::makeEfficiencyVsMomentum(
          m_finder_p2t.tree,
          "g_efficiency_vs_p",
          m_finderEfficiencyDefinition.value(),
          m_finderPurityThreshold.value(),
          0.1, 100.0, 0.15);
      if (g_eff_vs_p) {
        TCanvas* c_eff_vs_p = TrackingValidationPlots::drawEfficiencyCanvas(
            g_eff_vs_p,
            "c_efficiency_vs_p",
            "tracking efficiency vs momentum;p [GeV];Efficiency",
            0.1, 100.0);
        g_eff_vs_p->Write();
        if (c_eff_vs_p) c_eff_vs_p->Write();
        
      }
    
      m_outFile->Close();
    }
    return StatusCode::SUCCESS;
  }

private:
  // ---------- properties ----------
  Gaudi::Property<std::string> m_outputFile{this, "OutputFile", "validation.root", "Output ROOT file (TTrees)"};

  // 0 full, 1 finder-only, 2 fitter-only
  Gaudi::Property<int> m_mode{this, "Mode", 0, "Validation mode: 0=Full, 1=FinderOnly, 2=FitterOnly"};

  Gaudi::Property<float> m_Bz{this, "Bz", 2.f, "Magnetic field Bz [T] used in omega convention (GenFit-style)"};

  Gaudi::Property<float> m_refX{this, "RefPointX", 0.f, "Reference point X [mm] (must match fitter m_VP_referencePoint)"};
  Gaudi::Property<float> m_refY{this, "RefPointY", 0.f, "Reference point Y [mm] (must match fitter m_VP_referencePoint)"};
  Gaudi::Property<float> m_refZ{this, "RefPointZ", 0.f, "Reference point Z [mm] (must match fitter m_VP_referencePoint)"};
  
// Definition used for the summary tracking-efficiency plot.
// Default = 1 keeps the current behaviour unchanged.
  Gaudi::Property<int> m_finderEfficiencyDefinition{
  this, "FinderEfficiencyDefinition", 1,
  "Definition used for the tracking-efficiency summary plot: "
  "1 = require purity >= FinderPurityThreshold; "
  "2 = require purity >= 0.5 and efficiency >= 0.5"
      };

  Gaudi::Property<float> m_finderPurityThreshold{
  this, "FinderPurityThreshold", 0.75f,
  "Minimum purity for a particle-track match to count in tracking efficiency "
  "when FinderEfficiencyDefinition = 1"
      };

  Gaudi::Property<bool> m_doPerfectFit{
      this, "DoPerfectFit", false,
      "If true: fill fitter_vs_perfect using PerfectFitted_tracks if available. "
      "If false: tree exists but is empty per event."};

  // ---------- output structs ----------
  struct AssocTree {
    TTree* tree = nullptr;
    int event = 0;
    std::vector<int> index;
    std::vector<float> p;
    std::vector<float> pT;
    std::vector<int> nTrueHits; 
    std::vector<std::vector<int>> assoc;

    // only really used for finder_particle_to_tracks
    std::vector<std::vector<int>> sharedHits;
    std::vector<std::vector<float>> matchEfficiency;
    std::vector<std::vector<float>> matchPurity;

    void clear() { 
        index.clear(); 
        p.clear();
        pT.clear();
        nTrueHits.clear();
        assoc.clear();
        sharedHits.clear();
        matchEfficiency.clear();
        matchPurity.clear();
     }
  };

  struct FitterTree {
    TTree* tree = nullptr;
    int event = 0;

    std::vector<int> track_index;
    std::vector<int> track_location;

    std::vector<float> resD0, resZ0, resPhi, resOmega, resTanL;
    std::vector<float> p_reco, p_ref;
    std::vector<float> pT_reco, pT_ref;

    void clear() {
      track_index.clear();
      track_location.clear();
      resD0.clear();
      resZ0.clear();
      resPhi.clear();
      resOmega.clear();
      resTanL.clear();
      p_reco.clear();
      p_ref.clear();
      pT_reco.clear();
      pT_ref.clear();
    }
  };

  static void bookAssocTree(AssocTree& t, const char* name) {
    t.tree = new TTree(name, name);
    t.tree->Branch("event", &t.event);
    t.tree->Branch("index", &t.index);
    t.tree->Branch("p", &t.p);
    t.tree->Branch("pT", &t.pT);
    t.tree->Branch("nTrueHits", &t.nTrueHits);
    t.tree->Branch("assoc", &t.assoc);
    t.tree->Branch("sharedHits", &t.sharedHits);
    t.tree->Branch("matchEfficiency", &t.matchEfficiency);
    t.tree->Branch("matchPurity", &t.matchPurity);
  }

  static void bookFitterTree(FitterTree& t, const char* name) {
    t.tree = new TTree(name, name);
    t.tree->Branch("event", &t.event);
    t.tree->Branch("track_index", &t.track_index);
    t.tree->Branch("track_location", &t.track_location);
    t.tree->Branch("resD0", &t.resD0);
    t.tree->Branch("resZ0", &t.resZ0);
    t.tree->Branch("resPhi", &t.resPhi);
    t.tree->Branch("resOmega", &t.resOmega);
    t.tree->Branch("resTanLambda", &t.resTanL);
    t.tree->Branch("p_reco", &t.p_reco);
    t.tree->Branch("p_ref", &t.p_ref);
    t.tree->Branch("pT_reco", &t.pT_reco);
    t.tree->Branch("pT_ref", &t.pT_ref);
  }

  // ---------- association trees ----------
  void fillPerfectAssoc(int event, const edm4hep::MCParticleCollection& mcParts,
                        const std::unordered_map<int, std::vector<uint64_t>>& hitsPerParticle) const {

    m_perf_p2t.clear();
    m_perf_t2p.clear();
    m_perf_p2t.event = event;
    m_perf_t2p.event = event;

    for (int i = 0; i < (int)mcParts.size(); ++i) {
      const auto& mc = mcParts[i];
      if (mc.getGeneratorStatus() != 1) continue;

      auto it = hitsPerParticle.find(i);
      if (it == hitsPerParticle.end() || it->second.empty()) continue;
      const auto& mom = mc.getMomentum();
      const float px = float(mom.x);
      const float py = float(mom.y);
      const float pz = float(mom.z);
      const float p  = std::sqrt(px*px + py*py + pz*pz);
      const float pT = std::sqrt(px*px + py*py);
      const int nHits = (int)it->second.size();

      m_perf_p2t.index.push_back(i);
      m_perf_p2t.p.push_back(p);
      m_perf_p2t.pT.push_back(pT);
      m_perf_p2t.nTrueHits.push_back(nHits);
      m_perf_p2t.assoc.push_back({i});
      m_perf_p2t.sharedHits.push_back({});
      m_perf_p2t.matchEfficiency.push_back({});
      m_perf_p2t.matchPurity.push_back({});

      m_perf_t2p.index.push_back(i);
      m_perf_t2p.p.push_back(p);
      m_perf_t2p.pT.push_back(pT);
      m_perf_t2p.nTrueHits.push_back(nHits);
      m_perf_t2p.assoc.push_back({i});
      m_perf_t2p.sharedHits.push_back({});
      m_perf_t2p.matchEfficiency.push_back({});
      m_perf_t2p.matchPurity.push_back({});
    }

    if (m_perf_p2t.tree) m_perf_p2t.tree->Fill();
    if (m_perf_t2p.tree) m_perf_t2p.tree->Fill();
  }

  void fillFinderAssoc(int event, const edm4hep::MCParticleCollection& mcParts,
                     const edm4hep::TrackCollection& finderTracks,
                     const std::unordered_map<uint64_t, int>& hitToParticle,
                     const std::unordered_map<int, std::vector<uint64_t>>& hitsPerParticle) const {

  m_finder_p2t.clear();
  m_finder_t2p.clear();
  m_finder_p2t.event = event;
  m_finder_t2p.event = event;

  std::vector<std::unordered_map<int, int>> trackParticleCounts;
  std::vector<int> trackNHits;
  trackParticleCounts.resize(finderTracks.size());
  trackNHits.resize(finderTracks.size(), 0);

  int tIdx = 0;
  for (const auto& trk : finderTracks) {
    trackNHits[tIdx] = (int)trk.getTrackerHits().size();
    for (const auto& h : trk.getTrackerHits()) {
      const uint64_t hk = TrackingValidationHelpers::oidKey(h.getObjectID());
      auto it = hitToParticle.find(hk);
      if (it == hitToParticle.end()) continue;
      trackParticleCounts[tIdx][it->second] += 1;
    }
    ++tIdx;
  }

  // track -> particles
  for (int t = 0; t < (int)finderTracks.size(); ++t) {
    m_finder_t2p.index.push_back(t);
    m_finder_t2p.p.push_back(-1.f);
    m_finder_t2p.pT.push_back(-1.f);
    m_finder_t2p.nTrueHits.push_back(trackNHits[t]);

    std::vector<int> parts;
    std::vector<int> sh;
    std::vector<float> effs;
    std::vector<float> purs;

    for (const auto& kv : trackParticleCounts[t]) {
      const int pid = kv.first;
      if (pid < 0) continue;

      const int shared = kv.second;
      const int nTrackHits = trackNHits[t];
      const int nParticleHits =
          hitsPerParticle.count(pid) ? (int)hitsPerParticle.at(pid).size() : 0;

      const float eff = (nParticleHits > 0) ? float(shared) / float(nParticleHits) : 0.f;
      const float pur = (nTrackHits > 0) ? float(shared) / float(nTrackHits) : 0.f;

      parts.push_back(pid);
      sh.push_back(shared);
      effs.push_back(eff);
      purs.push_back(pur);
    }

    m_finder_t2p.assoc.push_back(parts);
    m_finder_t2p.sharedHits.push_back(sh);
    m_finder_t2p.matchEfficiency.push_back(effs);
    m_finder_t2p.matchPurity.push_back(purs);
  }

  // particle -> tracks
  for (int p = 0; p < (int)mcParts.size(); ++p) {
    const auto& mc = mcParts[p];
    if (mc.getGeneratorStatus() != 1) continue;

    auto itHits = hitsPerParticle.find(p);
    if (itHits == hitsPerParticle.end() || itHits->second.empty()) continue;

    const auto& mom = mc.getMomentum();
    const float px = float(mom.x);
    const float py = float(mom.y);
    const float pz = float(mom.z);
    const float pAbs = std::sqrt(px*px + py*py + pz*pz);
    const float pT = std::sqrt(px*px + py*py);
    const int nParticleHits = (int)itHits->second.size();

    std::vector<int> tracks;
    std::vector<int> sh;
    std::vector<float> effs;
    std::vector<float> purs;

    for (int t = 0; t < (int)finderTracks.size(); ++t) {
      auto it = trackParticleCounts[t].find(p);
      if (it == trackParticleCounts[t].end()) continue;

      const int shared = it->second;
      const int nTrackHits = trackNHits[t];

      const float eff = (nParticleHits > 0) ? float(shared) / float(nParticleHits) : 0.f;
      const float pur = (nTrackHits > 0) ? float(shared) / float(nTrackHits) : 0.f;

      tracks.push_back(t);
      sh.push_back(shared);
      effs.push_back(eff);
      purs.push_back(pur);
    }

    m_finder_p2t.index.push_back(p);
    m_finder_p2t.p.push_back(pAbs);
    m_finder_p2t.pT.push_back(pT);
    m_finder_p2t.nTrueHits.push_back(nParticleHits);
    m_finder_p2t.assoc.push_back(tracks);
    m_finder_p2t.sharedHits.push_back(sh);
    m_finder_p2t.matchEfficiency.push_back(effs);
    m_finder_p2t.matchPurity.push_back(purs);
  }

  if (m_finder_p2t.tree) m_finder_p2t.tree->Fill();
  if (m_finder_t2p.tree) m_finder_t2p.tree->Fill();
}

  // ---------- matching helper ----------
  int majorityParticleForTrack(const edm4hep::Track& trk,
                               const std::unordered_map<uint64_t, int>& hitToParticle) const {
    std::unordered_map<int, int> counts;
    for (const auto& h : trk.getTrackerHits()) {
      const uint64_t hk = TrackingValidationHelpers::oidKey(h.getObjectID());
      auto it = hitToParticle.find(hk);
      if (it == hitToParticle.end()) continue;
      counts[it->second] += 1;
    }
    if (counts.empty()) return -1;

    int bestP = -1;
    int bestN = -1;
    for (const auto& kv : counts) {
      if (kv.second > bestN) {
        bestN = kv.second;
        bestP = kv.first;
      }
    }
    return bestP;
  }

  // ---------- fitter trees ----------
  template <typename PerfectMapT>
  void fillFitterTrees(int event,
                       const edm4hep::MCParticleCollection& mcParts,
                       const edm4hep::TrackCollection& fittedTracks,
                       const std::unordered_map<uint64_t, int>& hitToParticle,
                       const PerfectMapT& perfectAtIPByPid,
                       bool doPerfect) const {

    m_fit_vs_mc.clear();
    m_fit_vs_perfect.clear();
    m_fit_vs_mc.event = event;
    m_fit_vs_perfect.event = event;

    const float NaN = std::numeric_limits<float>::quiet_NaN();

    int tIdx = 0;
    for (const auto& trk : fittedTracks) {
      edm4hep::TrackState stReco;
      if (!TrackingValidationHelpers::getAtIPState(trk, stReco)) {
        ++tIdx;
        continue;
      }

      const int pid = majorityParticleForTrack(trk, hitToParticle);
      if (pid < 0 || pid >= (int)mcParts.size()) {
        ++tIdx;
        continue;
      }

      const auto& mc = mcParts[pid];

      // reco params (already in fitter convention)
      TrackingValidationHelpers::HelixParams reco;
      reco.D0 = float(stReco.D0);
      reco.Z0 = float(stReco.Z0);
      reco.phi = float(stReco.phi);
      reco.omega = float(stReco.omega);
      reco.tanLambda = float(stReco.tanLambda);
      reco.pT = TrackingValidationHelpers::ptFromState(stReco, m_Bz.value());
      reco.p = TrackingValidationHelpers::momentumFromState(stReco, m_Bz.value());

      // ref from MC using the SAME convention as fitter (PCA + phi0 + ZPCA + omega=a*B/pT)
      const TrackingValidationHelpers::HelixParams refMC = TrackingValidationHelpers::truthFromMC_GenfitConvention(mc, m_Bz.value(), m_refX.value(), m_refY.value(), m_refZ.value());
  
      
      // --- vs MC  ---
      m_fit_vs_mc.track_index.push_back(tIdx);
      m_fit_vs_mc.track_location.push_back(int(stReco.location));
      m_fit_vs_mc.resD0.push_back(reco.D0 - refMC.D0);
      m_fit_vs_mc.resZ0.push_back(reco.Z0 - refMC.Z0);
      m_fit_vs_mc.resPhi.push_back(TrackingValidationHelpers::wrapDeltaPhi(reco.phi, refMC.phi));
      m_fit_vs_mc.resOmega.push_back(reco.omega - refMC.omega);
      m_fit_vs_mc.resTanL.push_back(reco.tanLambda - refMC.tanLambda);
      m_fit_vs_mc.p_reco.push_back(reco.p);
      m_fit_vs_mc.p_ref.push_back(refMC.p);
      m_fit_vs_mc.pT_reco.push_back(reco.pT);
      m_fit_vs_mc.pT_ref.push_back(refMC.pT);

      // --- vs perfect-fitted  ---
      if (doPerfect) {
        auto it = perfectAtIPByPid.find(pid);
        if (it != perfectAtIPByPid.end()) {
          const auto& stPerf = it->second.st;

          TrackingValidationHelpers::HelixParams refP;
          refP.D0 = float(stPerf.D0);
          refP.Z0 = float(stPerf.Z0);
          refP.phi = float(stPerf.phi);
          refP.omega = float(stPerf.omega);
          refP.tanLambda = float(stPerf.tanLambda);
          refP.pT = TrackingValidationHelpers::ptFromState(stPerf, m_Bz.value());
          refP.p = TrackingValidationHelpers::momentumFromState(stPerf, m_Bz.value());
          
          m_fit_vs_perfect.track_index.push_back(tIdx);
          m_fit_vs_perfect.track_location.push_back(int(stReco.location));
          m_fit_vs_perfect.resD0.push_back(reco.D0 - refP.D0);
          m_fit_vs_perfect.resZ0.push_back(reco.Z0 - refP.Z0);
          m_fit_vs_perfect.resPhi.push_back(TrackingValidationHelpers::wrapDeltaPhi(reco.phi, refP.phi));
          m_fit_vs_perfect.resOmega.push_back(reco.omega - refP.omega);
          m_fit_vs_perfect.resTanL.push_back(reco.tanLambda - refP.tanLambda);
          m_fit_vs_perfect.p_reco.push_back(reco.p);
          m_fit_vs_perfect.p_ref.push_back(refP.p);
          m_fit_vs_perfect.pT_reco.push_back(reco.pT);
          m_fit_vs_perfect.pT_ref.push_back(refP.pT);
        } else {
          m_fit_vs_perfect.track_index.push_back(tIdx);
          m_fit_vs_perfect.track_location.push_back(int(stReco.location));
          m_fit_vs_perfect.resD0.push_back(NaN);
          m_fit_vs_perfect.resZ0.push_back(NaN);
          m_fit_vs_perfect.resPhi.push_back(NaN);
          m_fit_vs_perfect.resOmega.push_back(NaN);
          m_fit_vs_perfect.resTanL.push_back(NaN);
          m_fit_vs_perfect.p_reco.push_back(reco.p);
          m_fit_vs_perfect.p_ref.push_back(NaN);
          m_fit_vs_perfect.pT_reco.push_back(reco.pT);
          m_fit_vs_perfect.pT_ref.push_back(NaN);
        }
      }

      ++tIdx;
    }

    if (m_fit_vs_mc.tree) m_fit_vs_mc.tree->Fill();
    if (m_fit_vs_perfect.tree) m_fit_vs_perfect.tree->Fill();
  }

private:
  mutable int m_evt = 0;
  mutable bool m_warnedMissingPerfectInput = false;

  std::unique_ptr<TFile> m_outFile;

  mutable AssocTree m_finder_p2t;
  mutable AssocTree m_finder_t2p;
  mutable AssocTree m_perf_p2t;
  mutable AssocTree m_perf_t2p;

  mutable FitterTree m_fit_vs_mc;
  mutable FitterTree m_fit_vs_perfect;
};

DECLARE_COMPONENT(TrackingValidation)