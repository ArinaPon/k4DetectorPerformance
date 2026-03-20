# TrackingValidation

## Overview

`TrackingValidation` is a validation algorithm for studying the performance (efficiency, purity, residuals, resolutions) of track finding and track fitting in the tracking reconstruction.
It is desighned to compare reconstructed and fitted tracks with Monte Carlo truth information and, when enabled, with tracks obtained from perfect tracking i.e. tracks fitted using the correct simhits from the particle truth information. The algorithm writes validation information to a ROOT output file containing TTrees and summary plots that can be used later for performance studies and plotting.

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

- **Planar digi-to-sim link collections**  
  Type: `std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>`  
  Used to connect reconstructed planar hits to the originating simulated particles.

- **Drift-chamber digi-to-sim link collections**  
  Type: `std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>`  
  Used to connect reconstructed drift-chamber hits to the originating simulated particles.

- **Finder track collection**  
  Type: `edm4hep::TrackCollection`  
  Collection of tracks produced by the track-finding stage.

- **Fitted track collection**  
  Type: `edm4hep::TrackCollection`  
  Collection of tracks produced by the standard fitting stage.

- **Perfect fitted-track collections (optional)**  
  Type: `std::vector<const edm4hep::TrackCollection*>`  
  Optional reference collections produced from perfect truth-based associations, used when perfect-fit validation is enabled.

---

## Outputs

The algorithm writes a ROOT file specified by `OutputFile`.

The file contains validation TTrees for finder-level and fitter-level studies, together with summary performance plots produced in `finalize()`. The fitter validation trees store residuals of the reconstructed track parameters with respect to the chosen reference.

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

- tracking efficiency vs momentum,
- `d0` resolution vs momentum,
- momentum resolution vs momentum,
- transverse-momentum resolution vs momentum.  

---

## Finder validation: efficiency and purity

To evaluate finder performance, each reconstructed track is matched to the truth particle with which it shares the largest number of hits.

For each particle-track pair, the algorithm stores two standard hit-based quantities:

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

`TrackingValidation` is tested through a small end-to-end workflow driven by `ctest`. The test starts from a simulated EDM4hep file, runs the reconstruction and validation steering, and writes the final validation ROOT output.

In the current setup:

- the simulation step is performed with `ddsim` in the shell test,
- the reconstruction and validation steps are controlled by `runTrackingValidation.py`,
- the full test is launched through `ctest`.


### Build and run the test

Set up the Key4hep environment:

```bash
cd k4RecTracker
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
k4_local_repo
```

Then build and install `k4RecTracker`:

```bash

mkdir build install
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
```
Next, build and install `k4DetectorPerformance`, and expose its local install:

```bash
cd k4DetectorPerformance
k4_local_repo
mkdir build install
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
```

Finally, run the test from the `k4DetectorPerformance`` build directory:

```bash
cd k4DetectorPerformance/build
ctest -V -R testTrackingValidation
```

The validation output is written to:

```text
k4DetectorPerformance/TrackingPerformance/test/validation_output.root
```
### Test configuration and steering options

The current test runs the full reconstruction and validation chain with the following settings:

- `RUN_SIM = 1`  
  simulation step in the shell test (`0` = skip simulation and use an existing EDM4hep input file via `INPUT_FILE_OVERRIDE`, `1` = run simulation with `ddsim`);
  
- `runDigi = 1`  
  digitization step (`0` = skip digitization, `1` = run digitization);

- `runFinder = 1`  
  track-finder step (`0` = skip track finder, `1` = run track finder);

- `runFitter = 1`  
  reconstructed-track fitter (`0` = skip reco fitter, `1` = run reco fitter);

- `runPerfectTracking = 1`  
  perfect-tracking and perfect-fitter chain (`0` = skip perfect tracking/perfect fitter, `1` = run them);

- `runValidation = 1`  
  validation step (`0` = skip validation, `1` = run validation);

- `useDCH = 1`  
  drift-chamber collections (`0` = disable DCH, `1` = use DCH);

- `mode = 0`  
  validation mode (`0` = full validation, `1` = finder-only validation, `2` = fitter-only validation);

- `doPerfectFit = 1`  
  fitter-versus-perfect comparison (`0` = disable fitter-vs-perfect filling, `1` = enable it);

- `finderEfficiencyDefinition = 1`  
  tracking-efficiency definition (`1` = purity-based definition, `2` = purity >= 0.5 and efficiency >= 0.5);

- `finderPurityThreshold = 0.75`  
  purity threshold used when `FinderEfficiencyDefinition = 1`.

These command-line flags are defined in `runTrackingValidation.py`, which allows the same steering file to be used either for the full chain or for reduced workflows in which some reconstruction steps are skipped and only the validation is run.
