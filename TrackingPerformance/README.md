# TrackingValidation

## Overview

`TrackingValidation` is a validation algorithm for studying the performance of track finding and track fitting in the tracking reconstruction.
It is desighned to compare reconstructed and fitted tracks with Monte Carlo truth information and, when enabled, with tracks obtained from perfect tracking. The algorithm writes validation information to a ROOT output file containing TTrees that can be used later for performance studies and plotting.

Typical use cases include:
- validation of track-finder performance,
- validation of fitted-track parameters against MC truth,
- comparison between standard reconstructed tracks and perfectly associated reference tracks.

---

## Inputs

`TrackingValidation` expects an EDM4hep event content in which the relevant collections have already been produced by the preceding steps of the reconstruction chain.

### Input collection types

`TrackingValidation` consumes the following types of event collections:

- **MC particle collection**  
  Used as the truth reference for particle-level validation.

- **Planar digi-to-sim association collections**  
  Used to connect reconstructed planar hits to the originating simulated particles.

- **Wire-hit digi-to-sim association collection**  
  Used to connect reconstructed drift-chamber hits to the originating simulated particles.

- **Finder track collection**  
  Collection of tracks produced by the track-finding stage.

- **Fitted track collection**  
  Collection of tracks produced by the standard fitting stage.

- **Perfect fitted-track collection (optional)**  
  Collection of fitted tracks produced from perfect truth-based associations, used as an additional reference when perfect-fit validation is enabled.

  ---

  ## Outputs

  The algorithm writes a ROOT file specified by the `OutputFile`.

The exact content of the output depends mainly on the validation mode selected through `Mode`:

- **`Mode = 0` (full-pipeline mode)**  
  Produces both finder-level and fitter-level validation trees.

- **`Mode = 1` (finder-only mode)**  
  Produces the trees related to track-finder validation.

- **`Mode = 2` (fitter-only mode)**  
  Produces the trees related to fitted-track validation.

In addition, `DoPerfectFit` controls whether the comparison to perfectly associated fitted tracks is filled. When enabled, the output also includes the fitter-versus-perfect validation information.

---

## Finder validation: efficiency and purity

To evaluate finder performance, each reconstructed track is matched to the truth particle with which it shares the largest number of hits.

For each particle–track pair, the algorithm stores two standard hit-based quantities:

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

For more details on the CMS association convention and the related definitions of tracking efficiency, fake rate, and duplicate rate, see the CMS performance note *Performance of the track selection DNN in Run 3*. :contentReference[oaicite:1]{index=1}

---

## How to run

To run the `TrackingValidation` test locally, first set up the Key4hep environment:

```bash
cd k4RecTracker
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
```

Then build and install `k4RecTracker`:

```bash
k4_local_repo
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
```

Next, build and install `k4DetectorPerformance`, and expose its local install:

```bash
cd k4DetectorPerformance
k4_local_repo
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j 8
```

Finally, run the test from the `k4DetectorPerformance` build directory:

```bash
cd k4DetectorPerformance/build
ctest -V -R testTrackingValidation
```

The validation output is written to:

```text
k4DetectorPerformance/TrackingPerformance/test/validation_output_test.root
```


