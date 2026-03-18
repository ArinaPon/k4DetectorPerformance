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

#ifndef TRACKINGVALIDATIONPLOTS_H
#define TRACKINGVALIDATIONPLOTS_H

#include "TCanvas.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TTree.h"

#include <memory>
#include <string>
#include <vector>

namespace TrackingValidationPlots {

std::vector<double> makeLogBins(double min, double max, double step);
TF1* fitGaussianCore(TH1F* h, const std::string& name);

TGraphErrors* makeD0ResolutionVsMomentum(TTree* tree,
                                         const char* graphName = "g_d0_resolution_vs_p",
                                         double pMin = 0.1,
                                         double pMax = 100.0,
                                         double logStep = 0.15);

TCanvas* drawD0ResolutionCanvas(TGraphErrors* g,
                                const char* canvasName = "c_d0_resolution_vs_p",
                                double xMin = 0.1,
                                double xMax = 100.0);

TGraphErrors* makeMomentumResolutionVsMomentum(TTree* tree,
                                               const char* graphName = "g_p_resolution_vs_p",
                                               double pMin = 0.1,
                                               double pMax = 100.0,
                                               double logStep = 0.15);

TGraphErrors* makePtResolutionVsMomentum(TTree* tree,
                                         const char* graphName = "g_pt_resolution_vs_p",
                                         double pMin = 0.1,
                                         double pMax = 100.0,
                                         double logStep = 0.15);

TCanvas* drawResolutionCanvas(TGraphErrors* g,
                              const char* canvasName,
                              const char* title,
                              double xMin = 0.1,
                              double xMax = 100.0);

TGraphErrors* makeEfficiencyVsMomentum(TTree* finderTree,
                                       const char* graphName,
                                       int efficiencyDefinition,
                                       double purityThreshold,
                                       double pMin = 0.1,
                                       double pMax = 100.0,
                                       double logStep = 0.15);

TCanvas* drawEfficiencyCanvas(TGraphErrors* g,
                              const char* canvasName,
                              const char* title,
                              double xMin = 0.1,
                              double xMax = 100.0);

} // namespace TrackingValidationPlots

#endif
