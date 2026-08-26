<!--
Copyright (c) 2020-2024 Key4hep-Project.

This file is part of Key4hep.
See https://key4hep.github.io/key4hep-doc/ for further info.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
-->
# TrackingValidation

## Overview

`TrackingValidation` is a validation algorithm for studying the performance (efficiency, purity, residuals, resolutions and pulls) of track finding and track fitting in the tracking reconstruction.
It is designed to compare reconstructed and fitted tracks with Monte Carlo truth information and, when enabled, with tracks obtained from perfect tracking, i.e. tracks fitted using the correct simHits from the particle truth information. The algorithm writes validation information to a ROOT output file containing TTrees and summary plots that can be used later for performance studies and plotting.

Typical use cases include:
- validation of track-finder performance,
- validation of fitted-track parameters against MC truth,
- comparison between standard reconstructed tracks and perfectly associated reference tracks.

---

## Inputs

`TrackingValidation` expects EDM4hep event content in which the relevant collections have already been produced by the preceding steps of the reconstruction chain.

### Input collection types

`TrackingValidation` consumes the following input collections:

- **MC particle collection**
  Type: `edm4hep::MCParticleCollection`

  Used as the truth reference for particle-level validation.

- **Hit-to-sim link collections**
  Type: `std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>`

  Collections used to associate reconstructed tracker hits with the originating simulated particles. The steering combines the available detector-specific link collections (e.g. silicon and drift-chamber links) into a single input vector passed to the validation algorithm.

- **Finder track collection**
  Type: `edm4hep::TrackCollection`

  Collection of tracks produced by the track-finding stage.

- **Fitted track collection**
  Type: `edm4hep::TrackCollection`

  Collection of tracks produced by the standard fitting stage.

- **Perfect fitted-track collections (optional)**
  Type: `std::vector<const edm4hep::TrackCollection*>`

  Optional reference collections produced from perfect truth-based associations and fitting, used when perfect-fit validation is enabled.

---

## Outputs

The algorithm writes a ROOT file specified by `OutputFile`.

The file contains validation TTrees for finder-level and fitter-level studies, together with summary performance plots produced in `finalize()`. The fitter validation trees store the reconstructed and reference track parameters, their residuals, the corresponding parameter uncertainties extracted from the fitted covariance matrix, and the resulting pull values for the five helix parameters.

### Output content by mode

The exact content filled in the output depends on the validation mode selected through `Mode`:

- **`Mode = 0` (full-pipeline mode)**
  Both finder-level and fitter-level validation are performed.
  The output includes the association trees and the fitter residual trees.

- **`Mode = 1` (finder-only mode)**
  Only the finder-level validation is performed.
  The finder and perfect-association trees are filled, while the fitter trees are booked in the file but are not filled.

- **`Mode = 2` (fitter-only mode)**
  Only the fitter-level validation is performed.
  The fitter trees are filled, while the finder and perfect-association trees are booked in the file but are not filled.

### Effect of `DoPerfectFit`

The flag `DoPerfectFit` controls the handling of the `fitter_vs_perfect` output:

- if **`DoPerfectFit = true`** and perfect fitted-track collections are provided, the fitter-to-perfect comparison is filled;
- if **`DoPerfectFit = false`**, the `fitter_vs_perfect` tree is still created but its per-event content remains empty;
- if **`DoPerfectFit = true`** but no perfect fitted-track collection is provided, the tree is still written and a warning is issued.

### Summary plots

In `finalize()`, the algorithm also writes summary plots to the same ROOT file, including:

- tracking efficiency vs momentum;
- helix-parameter resolutions (`d0`, `z0`, `phi`, `omega`, `tanLambda`) as a function of momentum;
- total-momentum resolution vs momentum;
- transverse-momentum resolution vs momentum;
- pull distributions for the five helix parameters (`d0`, `z0`, `phi`, `omega`, `tanLambda`).

The pull distributions are computed as
pull = (reconstructed − reference) / σ
where σ is taken from the corresponding diagonal element of the fitted track covariance matrix. Gaussian fits are performed for sufficiently populated pull distributions to facilitate validation of the covariance estimates.
Additional plots may be added in future developments.
---

## Finder validation: efficiency and purity

To evaluate finder performance, each reconstructed track is matched to the truth particle with which it shares the largest number of hits.

For each particle -track pair, the algorithm stores two standard hit-based quantities:

- **track hit purity**: the fraction of hits on the reconstructed track that originate from the matched truth particle;
- **track hit efficiency**: the fraction of the truth-particle hits that are recovered in the reconstructed track.

The summary **tracking efficiency** can then be defined in more than one way.

- **`FinderEfficiencyDefinition = 1`**
  A truth particle is counted as reconstructed if it is associated to at least one finder track with
  `purity >= FinderPurityThreshold`.
  In the default configuration, `FinderPurityThreshold = 0.75`, following the CMS association convention in which a reconstructed track is associated to a simulated particle if more than 75% of its hits originate from that particle. The tracking efficiency is then defined as the fraction of simulated tracks associated to at least one reconstructed track. :contentReference[oaicite:0]{index=0}

- **`FinderEfficiencyDefinition = 2`**
  A truth particle is counted as reconstructed if it is associated to at least one finder track with
  `purity >= 0.5` **and** `efficiency >= 0.5`.
  This corresponds to the stricter two-ratio variant, where both the purity of the reconstructed track and the fraction of recovered truth hits must exceed 50%.

In the current implementation, the denominator of the efficiency plot includes generator-level particles with status 1 and at least one truth-linked hit.
For more details on the CMS association convention and the related definitions of tracking efficiency, fake rate, and duplicate rate, see the CMS performance note [*Performance of the track selection DNN in Run 3*](https://cds.cern.ch/record/2854696/files/DP2023_009.pdf).




---

## How to run

`TrackingValidation` is tested through a lightweight `ctest` workflow. The CI test does **not** run the full reconstruction chain. Instead, it starts from a small pre-produced CLD reconstruction file and runs only the validation algorithm.

In the current setup:

- the input CLD reconstruction file is retrieved through CMake `ExternalData`;
- the validation steering is controlled by `runTrackingValidation.py`;
- the shell test `testTrackingValidation.sh` runs only `TrackingValidation`;

### Running the CI test

From the build directory run the validation test:

```bash
ctest -V -R testTrackingValidation
```

The test writes the validation output file to the test working directory, which is the build directory configured by CMake. The output file is named:

```text
validation.root
```

### CI test workflow

The registered CI test runs the following reduced workflow:

```text
pre-produced CLD reco file → TrackingValidation → validation.root
```

It intentionally does **not** run:

- DDSim;
- digitization;
- track finding;
- track fitting;
- perfect tracking.

The test uses the CLD geometry file from `k4geo`:

```text
${K4GEO}/FCCee/CLD/compact/CLD_o3_v01/CLD_o3_v01.xml
```

and the pre-produced input file:

```text
MuGuns_CLD_o3_v01_2026_07_01.root
```

which is provided through CMake `ExternalData`.

### CI test configuration

The validation-only CI test uses the following steering configuration:

- `runDigi = false`
  Do not run digitization.

- `runFinder = false`
  Do not run the track finder. Finder tracks are read from the input file.

- `runFitter = false`
  Do not run the track fitter. Fitted tracks are read from the input file.

- `runPerfectTracking = false`
  Do not run perfect tracking or perfect fitting.

- `runValidation = true`
  Run the validation algorithm.

- `useDCH = false`
  Do not use drift-chamber collections. The CI input is a CLD reconstruction file.

- `mode = 0`
  Run full validation, using both finder-level and fitter-level inputs from the pre-produced file.

- `doPerfectFit = false`
  Do not fill the fitter-versus-perfect-track comparison.

- `finderEfficiencyDefinition = 1`
  Use the purity-based tracking-efficiency definition.

- `finderPurityThreshold = 0.75`
  Use a purity threshold of 0.75 when `FinderEfficiencyDefinition = 1`.

The CLD collection names are passed explicitly to the steering file:

```text
mcParticles  = MCPhysicsParticles
hitSimLinks  = VXDTrackerHitRelations
finderTracks = SiTracks
fittedTracks = FittedTracks
```

The boolean steering options accept both `true/false` and `1/0` inputs.

---

## Plot recipes

The summary plots produced by `TrackingValidation` are written directly to the validation ROOT file during `finalize()`.

The recipes below describe how to reproduce the different plot families and which ROOT objects are produced.

### Common validation configuration

For an already reconstructed EDM4hep sample, the validation can be run without repeating digitization, track finding, or track fitting:

```bash
k4run runTrackingValidation.py \
  --inputFile <reconstructed.root> \
  --validationFile validation.root \
  --compactFile <detector-geometry.xml> \
  --runDigi false \
  --runFinder false \
  --runFitter false \
  --runPerfectTracking false \
  --runValidation true \
  --mode 0 \
  --doPerfectFit false \
  --mcParticles <MC-particle-collection> \
  --hitSimLinks <hit-to-sim-link-collections> \
  --finderTracks <finder-track-collection> \
  --fittedTracks <fitted-track-collection>
```

For the CLD sample used in the CI test, the corresponding collection names are:

```text
mcParticles  = MCPhysicsParticles
hitSimLinks  = VXDTrackerHitRelations
finderTracks = SiTracks
fittedTracks = FittedTracks
```

The validation mode determines which plot families are produced:

```text
Mode = 0   finder and fitter validation
Mode = 1   finder validation only
Mode = 2   fitter validation only
```

---

### Tracking efficiency vs momentum

The tracking-efficiency-versus-momentum plot is produced from the `finder_particle_to_tracks` tree.

It requires finder validation, i.e. `Mode = 0` or `Mode = 1`.

The matching definition is selected with:

```bash
--finderEfficiencyDefinition 1 \
--finderPurityThreshold 0.75
```

For `FinderEfficiencyDefinition = 1`, a truth particle is counted as reconstructed if at least one associated finder track satisfies

```text
purity >= FinderPurityThreshold
```

Alternatively,

```bash
--finderEfficiencyDefinition 2
```

requires at least one associated track satisfying

```text
purity >= 0.5
efficiency >= 0.5
```

Example:

```bash
k4run runTrackingValidation.py \
  --inputFile <reconstructed.root> \
  --validationFile validation.root \
  --compactFile <detector-geometry.xml> \
  --runDigi false \
  --runFinder false \
  --runFitter false \
  --runValidation true \
  --mode 1 \
  --finderEfficiencyDefinition 1 \
  --finderPurityThreshold 0.75 \
  --mcParticles <MC-particle-collection> \
  --hitSimLinks <hit-to-sim-link-collections> \
  --finderTracks <finder-track-collection>
```

The following objects are written to the validation ROOT file:

```text
g_efficiency_vs_p
c_efficiency_vs_p
```

The underlying information is stored in:

```text
finder_particle_to_tracks
```

---

### Tracking efficiency vs production radius

The tracking-efficiency-versus-production-radius plot is intended for studies of displaced-track reconstruction.

Enable the plot with:

```bash
--makeEfficiencyVsVertexR true
```

It requires finder validation (`Mode = 0` or `Mode = 1`).

A meaningful efficiency-versus-radius curve requires an input sample containing particles produced at non-zero transverse displacement. A prompt sample produced at the interaction point will populate only the first production-radius bin.

The current plotting range is

```text
0 <= vertexR < 1000 mm
```

with a bin width of `25 mm`.

#### Optional transverse-momentum selection

Enable the selection with:

```bash
--vertexREfficiencyApplyPtCut true \
--vertexREfficiencyMinPt 1.0
```

which applies

```text
pT > 1 GeV
```

#### Optional polar-angle selection

Enable the selection with:

```bash
--vertexREfficiencyApplyThetaCut true \
--vertexREfficiencyMinThetaDeg 10.0 \
--vertexREfficiencyMaxThetaDeg 170.0
```

which applies

```text
10 deg < theta < 170 deg
```

#### Optional MC-particle separation selection

Enable the selection with:

```bash
--vertexREfficiencyApplyDeltaMCCut true \
--vertexREfficiencyMinDeltaMC 0.02
```

which applies

```text
deltaMC > 0.02
```

Here, `deltaMC` is the minimum angular distance in eta-phi space between the selected MC particle and another selected MC particle.

#### Optional production-vertex-z selection

Enable the selection with:

```bash
--vertexREfficiencyApplyVertexZCut true \
--vertexREfficiencyMaxAbsVertexZ 30.0
```

which applies

```text
|vertexZ| <= 30 mm
```

A complete example with all selections enabled is:

```bash
k4run runTrackingValidation.py \
  --inputFile <reconstructed-displaced-sample.root> \
  --validationFile validation_displaced.root \
  --compactFile <detector-geometry.xml> \
  --runDigi false \
  --runFinder false \
  --runFitter false \
  --runValidation true \
  --mode 1 \
  --finderEfficiencyDefinition 1 \
  --finderPurityThreshold 0.75 \
  --makeEfficiencyVsVertexR true \
  --vertexREfficiencyApplyPtCut true \
  --vertexREfficiencyMinPt 1.0 \
  --vertexREfficiencyApplyThetaCut true \
  --vertexREfficiencyMinThetaDeg 10.0 \
  --vertexREfficiencyMaxThetaDeg 170.0 \
  --vertexREfficiencyApplyDeltaMCCut true \
  --vertexREfficiencyMinDeltaMC 0.02 \
  --vertexREfficiencyApplyVertexZCut true \
  --vertexREfficiencyMaxAbsVertexZ 30.0 \
  --mcParticles <MC-particle-collection> \
  --hitSimLinks <hit-to-sim-link-collections> \
  --finderTracks <finder-track-collection>
```

The following objects are written:

```text
g_efficiency_vs_vertexR
c_efficiency_vs_vertexR
```

The truth-level quantities used for this plot are stored in:

```text
finder_particle_to_tracks/pT
finder_particle_to_tracks/theta
finder_particle_to_tracks/vertexR
finder_particle_to_tracks/vertexZ
finder_particle_to_tracks/deltaMC
```

---

### Track-parameter resolutions vs momentum

Fitter validation automatically produces resolution plots for the five helix parameters:

```text
d0
z0
phi
omega
tanLambda
```

These plots require fitter validation (`Mode = 0` or `Mode = 2`). No additional plot-specific steering option is required.

A fitter-only validation run can be performed with:

```bash
k4run runTrackingValidation.py \
  --inputFile <reconstructed.root> \
  --validationFile validation.root \
  --compactFile <detector-geometry.xml> \
  --runDigi false \
  --runFinder false \
  --runFitter false \
  --runValidation true \
  --mode 2 \
  --doPerfectFit false \
  --mcParticles <MC-particle-collection> \
  --hitSimLinks <hit-to-sim-link-collections> \
  --fittedTracks <fitted-track-collection>
```

The following graphs and canvases are written:

```text
g_d0_resolution_vs_p
c_d0_resolution_vs_p

g_z0_resolution_vs_p
c_z0_resolution_vs_p

g_phi_resolution_vs_p
c_phi_resolution_vs_p

g_omega_resolution_vs_p
c_omega_resolution_vs_p

g_tanlambda_resolution_vs_p
c_tanlambda_resolution_vs_p
```

The corresponding residual branches are:

```text
fitter_vs_mc/resD0
fitter_vs_mc/resZ0
fitter_vs_mc/resPhi
fitter_vs_mc/resOmega
fitter_vs_mc/resTanLambda
```

The resolution in each momentum bin is defined using the effective sigma, calculated as half the width of the narrowest interval containing `68.27%` of the residual distribution.

---

### Momentum and transverse-momentum resolutions

The fitter-validation run also produces the relative total-momentum and transverse-momentum resolutions.

The corresponding residual quantities are

```text
(p_reco - p_ref) / p_ref
```

and

```text
(pT_reco - pT_ref) / pT_ref
```

respectively.

No additional steering option is required.

The following objects are written:

```text
g_p_resolution_vs_p
c_p_resolution_vs_p

g_pt_resolution_vs_p
c_pt_resolution_vs_p
```

The underlying quantities are stored in:

```text
fitter_vs_mc/p_reco
fitter_vs_mc/p_ref
fitter_vs_mc/pT_reco
fitter_vs_mc/pT_ref
```

---

### Pull distributions

Fitter validation automatically produces pull distributions for the five helix parameters:

```text
d0
z0
phi
omega
tanLambda
```

For a parameter `x`, the pull is defined as

```text
(x_reco - x_ref) / sigma_x
```

where `sigma_x` is obtained from the corresponding diagonal element of the fitted track covariance matrix.

No additional steering option is required.

The following histograms are written:

```text
h_pull_d0
h_pull_z0
h_pull_phi
h_pull_omega
h_pull_tanlambda
```

with the corresponding canvases:

```text
c_h_pull_d0_vs_mc
c_h_pull_z0_vs_mc
c_h_pull_phi_vs_mc
c_h_pull_omega_vs_mc
c_h_pull_tanlambda_vs_mc
```

The underlying branches are:

```text
fitter_vs_mc/pullD0
fitter_vs_mc/pullZ0
fitter_vs_mc/pullPhi
fitter_vs_mc/pullOmega
fitter_vs_mc/pullTanLambda
```

For sufficiently populated pull distributions, a Gaussian fit is performed in the central region.

---

### Track-fit chi2/ndf distribution

Fitter validation also produces the reduced-chi2 distribution for MC-associated fitted tracks with a valid track state at the interaction point.

For each track,

```text
chi2Ndf = chi2 / ndf
```

is calculated when `chi2` is finite and non-negative and `ndf > 0`.

No additional steering option is required. The plot is produced when fitter validation is enabled (`Mode = 0` or `Mode = 2`).

The following objects are written:

```text
h_chi2_ndf
c_chi2_ndf
```

The underlying fit-quality quantities are stored in:

```text
fitter_vs_mc/chi2
fitter_vs_mc/ndf
fitter_vs_mc/chi2Ndf
```

The histogram currently contains `100` bins in the range

```text
0 <= chi2/ndf < 10
```

---

### Perfect-track comparison

Residuals with respect to perfectly associated fitted tracks can additionally be stored by enabling perfect tracking and perfect-fit validation:

```bash
--runPerfectTracking true \
--doPerfectFit true
```

The resulting comparison is stored in:

```text
fitter_vs_perfect
```

When `runPerfectTracking = false`, an already reconstructed input sample may instead provide the perfect fitted-track collection directly through the corresponding steering option.

---

### Inspecting the produced plots

The contents of the validation ROOT file can be listed with:

```bash
rootls -r validation.root
```

To list the main summary plots:

```bash
rootls -r validation.root | grep -Ei "efficiency|resolution|pull|chi2"
```


