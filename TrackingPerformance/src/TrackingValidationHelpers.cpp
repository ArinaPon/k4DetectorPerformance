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

#include <cmath>
#include <limits>
#include <algorithm>

namespace TrackingValidationHelpers {

// Conversion constant for pT in GeV, B in T, and curvature in mm^-1.
static constexpr float c_mm_s = 2.998e11f;
static constexpr float kBFieldToCurvature = 1e-15f * c_mm_s;

float wrapDeltaPhi(float a, float b) {
  constexpr float pi = 3.14159265358979323846f;
  constexpr float twoPi = 2.f * pi;

  float d = a - b;

  while (d > pi)
    d -= twoPi;

  while (d < -pi)
    d += twoPi;

  return d;
}

PCAInfoHelper PCAInfo_mm(float x, float y, float z, float px, float py, float pz, int chargeSign, float refX,
                         float refY, float Bz) {
  PCAInfoHelper out;

  const float pt = std::sqrt(px * px + py * py);
  if (pt == 0.f)
    return out;

  if (chargeSign == 0)
    return out;

  // This helper currently assumes that Bz is supplied as a positive magnitude.
  if (Bz <= 0.f)
    return out;

  const float R = pt / (kBFieldToCurvature * Bz);

  const float tx = px / pt;
  const float ty = py / pt;

  const float nx = float(chargeSign) * ty;
  const float ny = float(chargeSign) * (-tx);

  const float xc = x + R * nx;
  const float yc = y + R * ny;

  const float vx = refX - xc;
  const float vy = refY - yc;
  const float vxy = std::sqrt(vx * vx + vy * vy);
  if (vxy == 0.f)
    return out;

  const float ux = vx / vxy;
  const float uy = vy / vxy;

  const float pcaX = xc + R * ux;
  const float pcaY = yc + R * uy;

  const float rx = pcaX - xc;
  const float ry = pcaY - yc;

  const int sign = (chargeSign > 0) ? 1 : -1;

  // Tangent at the PCA, oriented along the particle momentum.
  float tanX = sign * ry;
  float tanY = -sign * rx;

  const float tnorm = std::sqrt(tanX * tanX + tanY * tanY);
  if (tnorm == 0.f)
    return out;

  tanX /= tnorm;
  tanY /= tnorm;

  const float phi0 = std::atan2(tanY, tanX);

  // Signed transverse arc length from the production point to the PCA.
  const float startRx = x - xc;
  const float startRy = y - yc;

  const float cross = startRx * ry - startRy * rx;
  const float dot = startRx * rx + startRy * ry;
  const float deltaAlpha = std::atan2(cross, dot);

  const float signedArcLength = -sign * R * deltaAlpha;

  // Propagate z consistently along the helix.
  const float ZPCA = z + signedArcLength * pz / pt;

  out.pcaX = pcaX;
  out.pcaY = pcaY;
  out.pcaZ = ZPCA;
  out.phi0 = phi0;
  out.ok = true;
  return out;
}

HelixParams truthHelixParamsFromMC(const edm4hep::MCParticle& mc, float Bz, float refX, float refY, float refZ) {
  HelixParams hp;

  const auto& mom = mc.getMomentum();
  const float px = float(mom.x);
  const float py = float(mom.y);
  const float pz = float(mom.z);

  const float pT = std::sqrt(px * px + py * py);
  const float p = std::sqrt(px * px + py * py + pz * pz);

  hp.pT = pT;
  hp.p = p;

  int qSign = 0;
  if (mc.getCharge() > 0.f)
    qSign = 1;
  else if (mc.getCharge() < 0.f)
    qSign = -1;

  const auto& v = mc.getVertex();
  const float x = float(v.x);
  const float y = float(v.y);
  const float z = float(v.z);

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

  const float dx = info.pcaX - refX;
  const float dy = info.pcaY - refY;

  hp.D0 = dx * std::sin(info.phi0) - dy * std::cos(info.phi0);
  hp.Z0 = info.pcaZ - refZ;
  hp.phi = info.phi0;
  hp.tanLambda = (pT > 0.f) ? (pz / pT) : 0.f;
  hp.omega = (pT > 0.f) ? (kBFieldToCurvature * Bz / pT * float(qSign)) : 0.f;

  return hp;
}

std::vector<float> computeDeltaMC(
    const edm4hep::MCParticleCollection& mcParts,
    const std::vector<int>& particleIndices) {

  const float NaN = std::numeric_limits<float>::quiet_NaN();
  const float infinity = std::numeric_limits<float>::infinity();

  // These vectors use the same indexing as mcParts.
  std::vector<float> eta(mcParts.size(), NaN);
  std::vector<float> phi(mcParts.size(), NaN);
  std::vector<float> deltaMC(mcParts.size(), NaN);

  std::vector<int> validIndices;
  validIndices.reserve(particleIndices.size());

  // Compute eta and phi for the requested MC particles.
  for (const int index : particleIndices) {
    if (index < 0 ||
        index >= static_cast<int>(mcParts.size())) {
      continue;
    }

    const auto& momentum = mcParts[index].getMomentum();

    const float px = static_cast<float>(momentum.x);
    const float py = static_cast<float>(momentum.y);
    const float pz = static_cast<float>(momentum.z);

    const float pT = std::hypot(px, py);

    if (!std::isfinite(pT) || !(pT > 0.f)) {
      continue;
    }

    eta[index] = std::asinh(pz / pT);
    phi[index] = std::atan2(py, px);

    // A valid particle starts with no finite neighbour.
    deltaMC[index] = infinity;
    validIndices.push_back(index);
  }

  // Calculate each pair only once and update both particles.
  for (std::size_t i = 0; i < validIndices.size(); ++i) {
    const int first = validIndices[i];

    for (std::size_t j = i + 1;
         j < validIndices.size();
         ++j) {
      const int second = validIndices[j];

      // Protect against repeated indices in particleIndices.
      if (first == second) {
        continue;
      }

      const float deltaEta = eta[first] - eta[second];
      const float deltaPhi =
          wrapDeltaPhi(phi[first], phi[second]);

      const float deltaR =
          std::hypot(deltaEta, deltaPhi);

      if (!std::isfinite(deltaR)) {
        continue;
      }

      deltaMC[first] =
          std::min(deltaMC[first], deltaR);

      deltaMC[second] =
          std::min(deltaMC[second], deltaR);
    }
  }

  return deltaMC;
}


std::optional<edm4hep::TrackState> getAtIPState(const edm4hep::Track& trk) {
  for (const auto& st : trk.getTrackStates()) {
    if (st.location == edm4hep::TrackState::AtIP) {
      return st;
    }
  }
  return std::nullopt;
}

float ptFromState(const edm4hep::TrackState& st, float Bz) {
  const float omega = std::abs(float(st.omega));
  if (omega == 0.f)
    return 0.f;
  return kBFieldToCurvature * std::abs(Bz) / omega;
}

float momentumFromState(const edm4hep::TrackState& st, float Bz) {
  const float pT = ptFromState(st, Bz);
  const float tl = float(st.tanLambda);
  return pT * std::sqrt(1.f + tl * tl);
}

} // namespace TrackingValidationHelpers