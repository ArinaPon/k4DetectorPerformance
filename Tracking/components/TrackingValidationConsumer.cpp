// TrackingValidationConsumer
//
// Validation consumer that writes the following TTrees:
//   1) finder_particle_to_tracks
//   2) finder_track_to_particles
//   3) perfect_particle_to_tracks
//   4) perfect_track_to_particles
//   5) fitter_vs_mc
//   6) fitter_vs_perfect
//
// In finalize(), the consumer also produces summary plots written to the same ROOT file:
//   - tracking efficiency vs momentum
//   - d0 resolution vs momentum
//   - momentum resolution vs momentum
//   - transverse-momentum resolution vs momentum
//
// The fitter-vs-perfect tree is filled only when perfect-fitted tracks are provided
// and DoPerfectFit is enabled

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

// podio
#include "podio/ObjectID.h"

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"

// STL
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

// ---------- helpers ----------
static inline uint64_t oidKey(const podio::ObjectID& id) {
  return (uint64_t(id.collectionID) << 32) | uint64_t(uint32_t(id.index));
}

static inline float safeAtan2(float y, float x) { return std::atan2(y, x); }

static inline float wrapDeltaPhi(float a, float b) {
  float d = a - b;
  while (d >  M_PI) d -= 2.f * M_PI;
  while (d < -M_PI) d += 2.f * M_PI;
  return d;
}

struct HelixParams {
  float D0 = 0.f;        // mm
  float Z0 = 0.f;        // mm
  float phi = 0.f;       // rad
  float omega = 0.f;     // 1/mm 
  float tanLambda = 0.f; // unitless
  float p = 0.f;         // GeV
  float pT = 0.f;        // GeV
};

// Constants matching the fitter code 
static constexpr float c_mm_s = 2.998e11f;
static constexpr float a_genfit = 1e-15f * c_mm_s; // ~2.998e-4

// --- GenFit-like PCAInfo in mm, ported from GenfitTrack::PCAInfo ---
// Returns PCA point (x,y,z) and Phi0 (tangent angle at PCA).
struct PCAInfoHelper {
  float pcaX = 0.f;
  float pcaY = 0.f;
  float pcaZ = 0.f;
  float phi0 = 0.f;
  bool ok = false;
};

// position (x,y,z) in mm, momentum (px,py,pz) in GeV, refPoint in mm
static PCAInfoHelper PCAInfo_mm(float x, float y, float z,
                                float px, float py, float pz,
                                int chargeSign,
                                float refX, float refY,
                                float Bz) {
  PCAInfoHelper out;

  const float pt = std::sqrt(px*px + py*py);
  if (pt == 0.f) return out;
  if (chargeSign == 0) chargeSign = 1;
  if (Bz == 0.f) return out;

  // Radius in mm:
  // GenfitTrack::PCAInfo uses R = pt/(0.3*|q|*Bz)*100  [cm]
  // -> multiply by 10 to get mm: *1000
  const float R = pt / (0.3f * std::abs(chargeSign) * Bz) * 1000.f;

  const float tx = px / pt;
  const float ty = py / pt;

  const float nx = float(chargeSign) * (ty);
  const float ny = float(chargeSign) * (-tx);

  const float xc = x + R * nx;
  const float yc = y + R * ny;

  const float vx = refX - xc;
  const float vy = refY - yc;
  const float vxy = std::sqrt(vx*vx + vy*vy);
  if (vxy == 0.f) return out;

  const float ux = vx / vxy;
  const float uy = vy / vxy;

  const float pcaX = xc + R * ux;
  const float pcaY = yc + R * uy;

  // tangent direction at PCA (same as GenfitTrack::PCAInfo)
  const float rx = pcaX - xc;
  const float ry = pcaY - yc;

  const int sign = (chargeSign > 0) ? 1 : -1;
  float tanX = -sign * ry;
  float tanY =  sign * rx;

  const float tnorm = std::sqrt(tanX*tanX + tanY*tanY);
  if (tnorm == 0.f) return out;

  tanX /= tnorm;
  tanY /= tnorm;

  const float phi0 = std::atan2(tanY, tanX);

  // ZPCA approximation from GenfitTrack::PCAInfo (ported)
  // Uses a straight-line minimization in (R,z) with pR=pt, pZ=pz.
  const float pR = pt;
  const float pZ = pz;
  const float R0 = std::sqrt(x*x + y*y);
  const float Z0 = z;

  const float denom = (pR*pR + pZ*pZ);
  if (denom == 0.f) return out;

  const float tPCA = -(R0*pR + Z0*pZ) / denom;
  const float ZPCA = Z0 + pZ * tPCA;

  out.pcaX = pcaX;
  out.pcaY = pcaY;
  out.pcaZ = ZPCA;
  out.phi0 = phi0;
  out.ok = true;
  return out;
}

// Build MC truth helix parameters using the fitter convention
static HelixParams truthFromMC_GenfitConvention(const edm4hep::MCParticle& mc,
                                                float Bz,
                                                float refX, float refY, float refZ) {
  HelixParams hp;

  const auto& mom = mc.getMomentum();
  const float px = float(mom.x);
  const float py = float(mom.y);
  const float pz = float(mom.z);

  const float pT = std::sqrt(px*px + py*py);
  const float p  = std::sqrt(px*px + py*py + pz*pz);

  hp.pT = pT;
  hp.p  = p;

  // Charge sign consistent with fitter usage (sign matters for omega)
  int qSign = 1;
  if (mc.getCharge() < 0.f) qSign = -1;

  const auto& v = mc.getVertex();
  const float x = float(v.x); // mm
  const float y = float(v.y); // mm
  const float z = float(v.z); // mm

  const auto info = PCAInfo_mm(x, y, z, px, py, pz, qSign, refX, refY, Bz);
  if (!info.ok) {
    const float NaN = std::numeric_limits<float>::quiet_NaN();
    hp.D0 = NaN;
    hp.Z0 = NaN;
    hp.phi = NaN;
    hp.omega = NaN;
    hp.tanLambda = NaN;
    return hp;
  }

  // D0/Z0 in the same convention as the fitter (mm)
  hp.D0 = ( (-(refX - info.pcaX)) * std::sin(info.phi0) + (refY - info.pcaY) * std::cos(info.phi0) ); // mm
  hp.Z0 = (info.pcaZ - refZ); // mm

  // phi in fitter is taken from momentum.Phi() at the evaluated state.
  // For MC we use the momentum direction.
  hp.phi = safeAtan2(py, px);

  hp.tanLambda = (pT > 0.f) ? (pz / pT) : 0.f;

  // omega convention matches fitter: omega = +/- |a * Bz / pT|
  hp.omega = (pT > 0.f) ? (std::abs(a_genfit * Bz / pT) * float(qSign)) : 0.f;

  return hp;
}

static bool getAtIPState(const edm4hep::Track& trk, edm4hep::TrackState& out) {
  for (const auto& st : trk.getTrackStates()) {
    if (st.location == edm4hep::TrackState::AtIP) {
      out = st;
      return true;
    }
  }
  return false;
}




static float ptFromState(const edm4hep::TrackState& st, float Bz) {
  const float omega = std::abs(float(st.omega));
  if (omega == 0.f) return 0.f;
  return a_genfit * std::abs(Bz) / omega;
}

static float momentumFromState(const edm4hep::TrackState& st, float Bz) {
  const float pT = ptFromState(st, Bz);
  const float tl = float(st.tanLambda);
  return pT * std::sqrt(1.f + tl * tl);
}
// Helper functions for plotting

static std::vector<double> makeLogBins(double min, double max, double step) {
  std::vector<double> bins;
  for (double x = std::log10(min); x <= std::log10(max); x += step) {
    bins.push_back(std::pow(10., x));
  }
  if (bins.empty() || bins.back() < max) bins.push_back(max);
  return bins;
}

static TF1* fitGaussianCore(TH1F* h, const std::string& name) {
  if (!h || h->GetEntries() < 10) return nullptr;

  const double mean = h->GetMean();
  const double rms  = h->GetRMS();
  if (rms <= 0.) return nullptr;

  TF1* g1 = new TF1((name + "_g1").c_str(), "gaus", mean - 2.0 * rms, mean + 2.0 * rms);
  h->Fit(g1, "RQ0");

  double m1 = g1->GetParameter(1);
  double s1 = std::abs(g1->GetParameter(2));
  if (s1 <= 0.) s1 = rms;

  TF1* g2 = new TF1(name.c_str(), "gaus", m1 - 1.5 * s1, m1 + 1.5 * s1);
  h->Fit(g2, "RQ0");

  return g2;
}

static TGraphErrors* makeD0ResolutionVsMomentum(TTree* tree,
                                                const char* graphName = "g_d0_resolution_vs_p",
                                                double pMin = 0.1,
                                                double pMax = 100.0,
                                                double logStep = 0.15) {
  if (!tree) return nullptr;

  std::vector<double> bins = makeLogBins(pMin, pMax, logStep);
  const int nBins = bins.size() - 1;

  std::vector<std::unique_ptr<TH1F>> hists;
  hists.reserve(nBins);

  for (int i = 0; i < nBins; ++i) {
    hists.emplace_back(std::make_unique<TH1F>(
        Form("h_d0_bin_%d", i),
        Form("d0 residual bin %d;resD0 [#mum];Entries", i),
        120, -20.0, 20.0));
  }

  
  std::vector<float>* resD0 = nullptr;
  std::vector<float>* p_ref_vec = nullptr;

  // Tree stores vectors per event
  tree->SetBranchAddress("p_ref", &p_ref_vec);
  tree->SetBranchAddress("resD0", &resD0);

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t ievt = 0; ievt < nEntries; ++ievt) {
    tree->GetEntry(ievt);

    if (!p_ref_vec || !resD0) continue;
    if (p_ref_vec->size() != resD0->size()) continue;

    for (size_t i = 0; i < p_ref_vec->size(); ++i) {
      const double p = (*p_ref_vec)[i];
      const double d0_um = (*resD0)[i] * 1000.0; // mm -> um

      if (!std::isfinite(p) || !std::isfinite(d0_um)) continue;
      if (p < pMin || p >= pMax) continue;

      int bin = -1;
      for (int b = 0; b < nBins; ++b) {
        if (p >= bins[b] && p < bins[b + 1]) {
          bin = b;
          break;
        }
      }
      if (bin < 0) continue;

      hists[bin]->Fill(d0_um);
    }
  }

  TGraphErrors* g = new TGraphErrors();
  g->SetName(graphName);
  g->SetTitle(";p_{ref} [GeV];#sigma(d_{0}) [#mum]");

  int ip = 0;
  for (int b = 0; b < nBins; ++b) {
    if (hists[b]->GetEntries() < 20) continue;

    TF1* fit = fitGaussianCore(hists[b].get(), Form("fit_d0_bin_%d", b));
    if (!fit) continue;

    const double sigma = std::abs(fit->GetParameter(2));
    const double sigmaErr = fit->GetParError(2);
    const double pCenter = std::sqrt(bins[b] * bins[b + 1]);

    g->SetPoint(ip, pCenter, sigma);
    g->SetPointError(ip, 0.0, sigmaErr);
    ++ip;
  }

  return g;
}

static TCanvas* drawD0ResolutionCanvas(TGraphErrors* g,
                                       const char* canvasName = "c_d0_resolution_vs_p",
                                       double xMin = 0.1,
                                       double xMax = 100.0) {
  if (!g) return nullptr;

  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas(canvasName, "d0 resolution vs momentum", 800, 600);
  c->SetLogx();

  g->SetMarkerStyle(20);
  g->SetLineWidth(2);
  g->GetXaxis()->SetLimits(xMin, xMax);
  g->Draw("AP");

  return c;
}

static TGraphErrors* makeMomentumResolutionVsMomentum(TTree* tree,
                                                      const char* graphName = "g_p_resolution_vs_p",
                                                      double pMin = 0.1,
                                                      double pMax = 100.0,
                                                      double logStep = 0.15) {
  if (!tree) return nullptr;

  std::vector<double> bins = makeLogBins(pMin, pMax, logStep);
  const int nBins = bins.size() - 1;

  std::vector<std::unique_ptr<TH1F>> hists;
  hists.reserve(nBins);

  for (int i = 0; i < nBins; ++i) {
    hists.emplace_back(std::make_unique<TH1F>(
        Form("h_pres_bin_%d", i),
        Form("p resolution bin %d;(p_{reco}-p_{ref})/p_{ref};Entries", i),
        120, -0.2, 0.2));
  }

  std::vector<float>* p_ref_vec = nullptr;
  std::vector<float>* p_reco_vec = nullptr;

  tree->SetBranchAddress("p_ref", &p_ref_vec);
  tree->SetBranchAddress("p_reco", &p_reco_vec);

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t ievt = 0; ievt < nEntries; ++ievt) {
    tree->GetEntry(ievt);

    if (!p_ref_vec || !p_reco_vec) continue;
    if (p_ref_vec->size() != p_reco_vec->size()) continue;

    for (size_t i = 0; i < p_ref_vec->size(); ++i) {
      const double pRef = (*p_ref_vec)[i];
      const double pReco = (*p_reco_vec)[i];

      if (!std::isfinite(pRef) || !std::isfinite(pReco)) continue;
      if (pRef <= 0.) continue;
      if (pRef < pMin || pRef >= pMax) continue;

      const double res = (pReco - pRef) / pRef;

      int bin = -1;
      for (int b = 0; b < nBins; ++b) {
        if (pRef >= bins[b] && pRef < bins[b + 1]) {
          bin = b;
          break;
        }
      }
      if (bin < 0) continue;

      hists[bin]->Fill(res);
    }
  }

  TGraphErrors* g = new TGraphErrors();
  g->SetName(graphName);
  g->SetTitle(";p_{ref} [GeV];#sigma((p_{reco}-p_{ref})/p_{ref})");

  int ip = 0;
  for (int b = 0; b < nBins; ++b) {
    if (hists[b]->GetEntries() < 20) continue;

    TF1* fit = fitGaussianCore(hists[b].get(), Form("fit_pres_bin_%d", b));
    if (!fit) continue;

    const double sigma = std::abs(fit->GetParameter(2));
    const double sigmaErr = fit->GetParError(2);
    const double pCenter = std::sqrt(bins[b] * bins[b + 1]);

    g->SetPoint(ip, pCenter, sigma);
    g->SetPointError(ip, 0.0, sigmaErr);
    ++ip;
  }

  return g;
}

static TGraphErrors* makePtResolutionVsMomentum(TTree* tree,
                                                const char* graphName = "g_pt_resolution_vs_p",
                                                double pMin = 0.1,
                                                double pMax = 100.0,
                                                double logStep = 0.15) {
  if (!tree) return nullptr;

  std::vector<double> bins = makeLogBins(pMin, pMax, logStep);
  const int nBins = bins.size() - 1;

  std::vector<std::unique_ptr<TH1F>> hists;
  hists.reserve(nBins);

  for (int i = 0; i < nBins; ++i) {
    hists.emplace_back(std::make_unique<TH1F>(
        Form("h_ptres_bin_%d", i),
        Form("pT resolution bin %d;(pT_{reco}-pT_{ref})/pT_{ref};Entries", i),
        120, -0.2, 0.2));
  }

  std::vector<float>* p_ref_vec = nullptr;
  std::vector<float>* pt_ref_vec = nullptr;
  std::vector<float>* pt_reco_vec = nullptr;

  tree->SetBranchAddress("p_ref", &p_ref_vec);
  tree->SetBranchAddress("pT_ref", &pt_ref_vec);
  tree->SetBranchAddress("pT_reco", &pt_reco_vec);

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t ievt = 0; ievt < nEntries; ++ievt) {
    tree->GetEntry(ievt);

    if (!p_ref_vec || !pt_ref_vec || !pt_reco_vec) continue;
    if (p_ref_vec->size() != pt_ref_vec->size()) continue;
    if (pt_ref_vec->size() != pt_reco_vec->size()) continue;

    for (size_t i = 0; i < p_ref_vec->size(); ++i) {
      const double pRef = (*p_ref_vec)[i];
      const double ptRef = (*pt_ref_vec)[i];
      const double ptReco = (*pt_reco_vec)[i];

      if (!std::isfinite(pRef) || !std::isfinite(ptRef) || !std::isfinite(ptReco)) continue;
      if (ptRef <= 0.) continue;
      if (pRef < pMin || pRef >= pMax) continue;

      const double res = (ptReco - ptRef) / ptRef;

      int bin = -1;
      for (int b = 0; b < nBins; ++b) {
        if (pRef >= bins[b] && pRef < bins[b + 1]) {
          bin = b;
          break;
        }
      }
      if (bin < 0) continue;

      hists[bin]->Fill(res);
    }
  }

  TGraphErrors* g = new TGraphErrors();
  g->SetName(graphName);
  g->SetTitle(";p_{ref} [GeV];#sigma((pT_{reco}-pT_{ref})/pT_{ref})");

  int ip = 0;
  for (int b = 0; b < nBins; ++b) {
    if (hists[b]->GetEntries() < 20) continue;

    TF1* fit = fitGaussianCore(hists[b].get(), Form("fit_ptres_bin_%d", b));
    if (!fit) continue;

    const double sigma = std::abs(fit->GetParameter(2));
    const double sigmaErr = fit->GetParError(2);
    const double pCenter = std::sqrt(bins[b] * bins[b + 1]);

    g->SetPoint(ip, pCenter, sigma);
    g->SetPointError(ip, 0.0, sigmaErr);
    ++ip;
  }

  return g;
}

static TCanvas* drawResolutionCanvas(TGraphErrors* g,
                                     const char* canvasName,
                                     const char* title,
                                     double xMin = 0.1,
                                     double xMax = 100.0) {
  if (!g) return nullptr;

  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas(canvasName, title, 800, 600);
  c->SetLogx();

  g->SetMarkerStyle(20);
  g->SetLineWidth(2);
  g->SetTitle(title);
  g->GetXaxis()->SetLimits(xMin, xMax);
  g->Draw("AP");

  return c;
}
// A truth particle is counted as reconstructed if at least one associated finder track
// has purity above the configured threshold
static TGraphErrors* makeEfficiencyVsMomentum(TTree* finderTree,
                                              const char* graphName,
                                              double purityThreshold,
                                              double pMin = 0.1,
                                              double pMax = 100.0,
                                              double logStep = 0.15) {
  if (!finderTree) return nullptr;

  std::vector<double> bins = makeLogBins(pMin, pMax, logStep);
  const int nBins = bins.size() - 1;

  std::vector<int> nDen(nBins, 0);
  std::vector<int> nNum(nBins, 0);

  std::vector<float>* pVec = nullptr;
  std::vector<std::vector<float>>* purVec = nullptr;

  finderTree->SetBranchAddress("p", &pVec);
  finderTree->SetBranchAddress("matchPurity", &purVec);

  const Long64_t nEntries = finderTree->GetEntries();
  for (Long64_t ievt = 0; ievt < nEntries; ++ievt) {
    finderTree->GetEntry(ievt);

    if (!pVec || !purVec) continue;
    if (pVec->size() != purVec->size()) continue;

    for (size_t i = 0; i < pVec->size(); ++i) {
      const double p = (*pVec)[i];
      if (!std::isfinite(p) || p < pMin || p >= pMax) continue;

      int bin = -1;
      for (int b = 0; b < nBins; ++b) {
        if (p >= bins[b] && p < bins[b + 1]) {
          bin = b;
          break;
        }
      }
      if (bin < 0) continue;

      // denominator: all truth particles present in finder_particle_to_tracks tree
      // (genStatus1 and with at least one true hit, as enforced in fillFinderAssoc)
      nDen[bin]++;

      bool isMatched = false;
      for (size_t j = 0; j < (*purVec)[i].size(); ++j) {
        if ((*purVec)[i][j] >= purityThreshold) {
          isMatched = true;
          break;
      }
  }
      if (isMatched) nNum[bin]++;
    }
  }

  TGraphErrors* g = new TGraphErrors();
  g->SetName(graphName);
  g->SetTitle(";p [GeV];Tracking efficiency");

  int ip = 0;
  for (int b = 0; b < nBins; ++b) {
    if (nDen[b] == 0) continue;

    const double eff = double(nNum[b]) / double(nDen[b]);
    const double err = std::sqrt(eff * (1.0 - eff) / double(nDen[b]));
    const double pCenter = std::sqrt(bins[b] * bins[b + 1]);

    g->SetPoint(ip, pCenter, eff);
    g->SetPointError(ip, 0.0, err);
    ++ip;
  }

  return g;
}

static TCanvas* drawEfficiencyCanvas(TGraphErrors* g,
                                     const char* canvasName,
                                     const char* title,
                                     double xMin = 0.1,
                                     double xMax = 100.0) {
  if (!g) return nullptr;

  gStyle->SetOptStat(0);

  TCanvas* c = new TCanvas(canvasName, title, 800, 600);
  c->SetLogx();

  g->SetMarkerStyle(20);
  g->SetLineWidth(2);
  g->SetTitle(title);
  g->GetYaxis()->SetRangeUser(0.0, 1.05);
  g->GetXaxis()->SetLimits(xMin, xMax);
  g->Draw("AP");

  return c;
}

// ---------- CONSUMER ----------
struct TrackingValidationConsumer final
    : k4FWCore::Consumer<void(
          const edm4hep::MCParticleCollection&,
          const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&,  // planar links
          const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&,  // optional DCH links
          const edm4hep::TrackCollection&,                                            // finder tracks
          const edm4hep::TrackCollection&,                                            // fitted tracks (reco)
          const std::vector<const edm4hep::TrackCollection*>&                         // optional perfect fitted tracks
          )> {

  TrackingValidationConsumer(const std::string& name, ISvcLocator* svcLoc)
      : Consumer(
            name, svcLoc,
            {
                KeyValues("MCParticles", {"MCParticles"}),

                KeyValues("PlanarLinks",
                          {"SiWrBSimDigiLinks", "SiWrDSimDigiLinks", "VTXBSimDigiLinks", "VTXDSimDigiLinks"}),

                
                KeyValues("DCHLinks", {"DCH_DigiSimAssociationCollection"}),

                KeyValues("FinderTracks", {"GGTFTracks"}),
                KeyValues("FittedTracks", {"FittedTracks"}),

                
                KeyValues("PerfectFittedTracks", {"PerfectFitted_tracks"}),
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
        const uint64_t key = oidKey(digi.getObjectID());
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
        const uint64_t key = oidKey(digi.getObjectID());
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
          if (!getAtIPState(trk, st)) continue;

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
    info() << "Finalizing TrackingValidationConsumer, wrote " << m_evt << " events" << endmsg;

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
      TGraphErrors* g_d0_vs_p = makeD0ResolutionVsMomentum(m_fit_vs_mc.tree,
                                                         "g_d0_resolution_vs_p",
                                                         0.1, 100.0, 0.15);
      if (g_d0_vs_p) {
        TCanvas* c_d0_vs_p = drawD0ResolutionCanvas(g_d0_vs_p,
                                                  "c_d0_resolution_vs_p",
                                                  0.1, 100.0);
        g_d0_vs_p->Write();
        if (c_d0_vs_p) c_d0_vs_p->Write();
      }
      // p resolution vs momentum
      TGraphErrors* g_p_vs_p = makeMomentumResolutionVsMomentum(m_fit_vs_mc.tree,
                                                          "g_p_resolution_vs_p",
                                                          0.1, 100.0, 0.15);
      if (g_p_vs_p) {
        TCanvas* c_p_vs_p = drawResolutionCanvas(g_p_vs_p,
                                           "c_p_resolution_vs_p",
                                           "momentum resolution vs momentum;p_{ref} [GeV];#sigma((p_{reco}-p_{ref})/p_{ref})",
                                           0.1, 100.0);
        g_p_vs_p->Write();
        if (c_p_vs_p) c_p_vs_p->Write();
      }

      // pT resolution vs momentum
      TGraphErrors* g_pt_vs_p = makePtResolutionVsMomentum(m_fit_vs_mc.tree,
                                                     "g_pt_resolution_vs_p",
                                                     0.1, 100.0, 0.15);
      if (g_pt_vs_p) {
        TCanvas* c_pt_vs_p = drawResolutionCanvas(g_pt_vs_p,
                                            "c_pt_resolution_vs_p",
                                            "pT resolution vs momentum;p_{ref} [GeV];#sigma((pT_{reco}-pT_{ref})/pT_{ref})",
                                            0.1, 100.0);
        g_pt_vs_p->Write();
        if (c_pt_vs_p) c_pt_vs_p->Write();
      }

      // finder summary plot
      TGraphErrors* g_eff_vs_p = makeEfficiencyVsMomentum(
          m_finder_p2t.tree,
          "g_efficiency_vs_p",
          m_finderPurityThreshold.value(),
          0.1, 100.0, 0.15);
      if (g_eff_vs_p) {
        TCanvas* c_eff_vs_p = drawEfficiencyCanvas(
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
  
  Gaudi::Property<float> m_finderPurityThreshold{
    this, "FinderPurityThreshold", 0.75f,
    "Minimum purity for a particle-track match to count in tracking efficiency"};

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
      const uint64_t hk = oidKey(h.getObjectID());
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
      const uint64_t hk = oidKey(h.getObjectID());
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
      if (!getAtIPState(trk, stReco)) {
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
      HelixParams reco;
      reco.D0 = float(stReco.D0);
      reco.Z0 = float(stReco.Z0);
      reco.phi = float(stReco.phi);
      reco.omega = float(stReco.omega);
      reco.tanLambda = float(stReco.tanLambda);
      reco.pT = ptFromState(stReco, m_Bz.value());
      reco.p = momentumFromState(stReco, m_Bz.value());

      // ref from MC using the SAME convention as fitter (PCA + phi0 + ZPCA + omega=a*B/pT)
      const HelixParams refMC = truthFromMC_GenfitConvention(mc, m_Bz.value(), m_refX.value(), m_refY.value(), m_refZ.value());
  
      
      // --- vs MC  ---
      m_fit_vs_mc.track_index.push_back(tIdx);
      m_fit_vs_mc.track_location.push_back(int(stReco.location));
      m_fit_vs_mc.resD0.push_back(reco.D0 - refMC.D0);
      m_fit_vs_mc.resZ0.push_back(reco.Z0 - refMC.Z0);
      m_fit_vs_mc.resPhi.push_back(wrapDeltaPhi(reco.phi, refMC.phi));
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

          HelixParams refP;
          refP.D0 = float(stPerf.D0);
          refP.Z0 = float(stPerf.Z0);
          refP.phi = float(stPerf.phi);
          refP.omega = float(stPerf.omega);
          refP.tanLambda = float(stPerf.tanLambda);
          refP.pT = ptFromState(stPerf, m_Bz.value());
          refP.p = momentumFromState(stPerf, m_Bz.value());
          
          m_fit_vs_perfect.track_index.push_back(tIdx);
          m_fit_vs_perfect.track_location.push_back(int(stReco.location));
          m_fit_vs_perfect.resD0.push_back(reco.D0 - refP.D0);
          m_fit_vs_perfect.resZ0.push_back(reco.Z0 - refP.Z0);
          m_fit_vs_perfect.resPhi.push_back(wrapDeltaPhi(reco.phi, refP.phi));
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

DECLARE_COMPONENT(TrackingValidationConsumer)