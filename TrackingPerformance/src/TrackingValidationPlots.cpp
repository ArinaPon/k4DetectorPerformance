#include "TrackingValidationPlots.h"
#include "TStyle.h"
#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

// makeLogBins
//fitGaussianCore
//makeD0ResolutionVsMomentum
//drawD0ResolutionCanvas
//makeMomentumResolutionVsMomentum
//makePtResolutionVsMomentum
//drawResolutionCanvas
//makeEfficiencyVsMomentum
//drawEfficiencyCanvas

namespace TrackingValidationPlots {
std::vector<double> makeLogBins(double min, double max, double step) {
  std::vector<double> bins;
  for (double x = std::log10(min); x <= std::log10(max); x += step) {
    bins.push_back(std::pow(10., x));
  }
  if (bins.empty() || bins.back() < max) bins.push_back(max);
  return bins;
}

TF1* fitGaussianCore(TH1F* h, const std::string& name) {
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

TGraphErrors* makeD0ResolutionVsMomentum(TTree* tree,
                                         const char* graphName,
                                         double pMin,
                                         double pMax,
                                         double logStep) {
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

  tree->SetBranchAddress("p_ref", &p_ref_vec);
  tree->SetBranchAddress("resD0", &resD0);

  const Long64_t nEntries = tree->GetEntries();
  for (Long64_t ievt = 0; ievt < nEntries; ++ievt) {
    tree->GetEntry(ievt);

    if (!p_ref_vec || !resD0) continue;
    if (p_ref_vec->size() != resD0->size()) continue;

    for (size_t i = 0; i < p_ref_vec->size(); ++i) {
      const double p = (*p_ref_vec)[i];
      const double d0_um = (*resD0)[i] * 1000.0;

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

TCanvas* drawD0ResolutionCanvas(TGraphErrors* g,
                                const char* canvasName,
                                double xMin,
                                double xMax) {
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

TGraphErrors* makeMomentumResolutionVsMomentum(TTree* tree,
                                               const char* graphName,
                                               double pMin,
                                               double pMax,
                                               double logStep) {
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

TGraphErrors* makePtResolutionVsMomentum(TTree* tree,
                                         const char* graphName,
                                         double pMin,
                                         double pMax,
                                         double logStep) {
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

TCanvas* drawResolutionCanvas(TGraphErrors* g,
                              const char* canvasName,
                              const char* title,
                              double xMin,
                              double xMax) {
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

TGraphErrors* makeEfficiencyVsMomentum(TTree* finderTree,
                                       const char* graphName,
                                       int efficiencyDefinition,
                                       double purityThreshold,
                                       double pMin,
                                       double pMax,
                                       double logStep) {
  if (!finderTree) return nullptr;

  std::vector<double> bins = makeLogBins(pMin, pMax, logStep);
  const int nBins = bins.size() - 1;

  std::vector<int> nDen(nBins, 0);
  std::vector<int> nNum(nBins, 0);

  std::vector<float>* pVec = nullptr;
  std::vector<std::vector<float>>* purVec = nullptr;
  std::vector<std::vector<float>>* effVec = nullptr;

  finderTree->SetBranchAddress("p", &pVec);
  finderTree->SetBranchAddress("matchPurity", &purVec);
  finderTree->SetBranchAddress("matchEfficiency", &effVec);

  const Long64_t nEntries = finderTree->GetEntries();
  for (Long64_t ievt = 0; ievt < nEntries; ++ievt) {
    finderTree->GetEntry(ievt);

    if (!pVec || !purVec || !effVec) continue;
    if (pVec->size() != purVec->size()) continue;
    if (pVec->size() != effVec->size()) continue;

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

      nDen[bin]++;

      bool isMatched = false;

      const auto& purities = (*purVec)[i];
      const auto& efficiencies = (*effVec)[i];
      const size_t nMatches = std::min(purities.size(), efficiencies.size());

      for (size_t j = 0; j < nMatches; ++j) {
        const float purity = purities[j];
        const float efficiency = efficiencies[j];

        if (efficiencyDefinition == 2) {
          if (purity >= 0.5f && efficiency >= 0.5f) {
            isMatched = true;
            break;
          }
        } else {
          if (purity >= purityThreshold) {
            isMatched = true;
            break;
          }
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

TCanvas* drawEfficiencyCanvas(TGraphErrors* g,
                              const char* canvasName,
                              const char* title,
                              double xMin,
                              double xMax) {
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
}
