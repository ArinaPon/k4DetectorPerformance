#
# Copyright (c) 2020-2024 Key4hep-Project.
#
# This file is part of Key4hep.
# See https://key4hep.github.io/key4hep-doc/ for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
import os
import math

from Gaudi.Configuration import INFO
from Configurables import EventDataSvc, GeoSvc, UniqueIDGenSvc, RndmGenSvc
from k4FWCore import IOSvc, ApplicationMgr
from k4FWCore.parseArgs import parser

# --------------------
# Arguments
# --------------------
parser.add_argument("--inputFile", required=True,
                    help="Input EDM4hep ROOT file")
parser.add_argument("--modelPath", default="",
                    help="Path to the GGTF ONNX model, required only if --runFinder 1")
parser.add_argument("--outputFile", default="out_reco.root",
                    help="Output EDM4hep ROOT file with reconstructed collections")
parser.add_argument("--validationFile", default="validation.root",
                    help="Output ROOT file written by TrackingValidation")

parser.add_argument("--compactFile",
                    default=os.path.join(os.environ["K4GEO"], "FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml"),
                    help="Detector geometry XML file")

# Pipeline control
parser.add_argument("--runDigi", type=int, default=1, choices=[0, 1],
                    help="0=skip digitization, 1=run digitization")
parser.add_argument("--runFinder", type=int, default=1, choices=[0, 1],
                    help="0=skip track finder, 1=run track finder")
parser.add_argument("--runFitter", type=int, default=1, choices=[0, 1],
                    help="0=skip reco fitter, 1=run reco fitter")
parser.add_argument("--runPerfectTracking", type=int, default=1, choices=[0, 1],
                    help="0=skip perfect tracking/perfect fitter, 1=run them")
parser.add_argument("--runValidation", type=int, default=1, choices=[0, 1],
                    help="0=skip validation, 1=run TrackingValidation")
parser.add_argument("--useDCH", type=int, default=1, choices=[0, 1],
                    help="0=disable DCH collections, 1=use DCH collections")

# Validation control
parser.add_argument("--mode", type=int, default=0, choices=[0, 1, 2],
                    help="Validation mode: 0=Full, 1=FinderOnly, 2=FitterOnly")
parser.add_argument("--doPerfectFit", type=int, default=1, choices=[0, 1],
                    help="0=do not fill fitter_vs_perfect, 1=fill fitter_vs_perfect")
parser.add_argument("--finderEfficiencyDefinition", type=int, default=1, choices=[1, 2],
                    help="1=purity-based definition, 2=purity+efficiency >= 0.5 definition")
parser.add_argument("--finderPurityThreshold", type=float, default=0.75,
                    help="Purity threshold used when FinderEfficiencyDefinition = 1")

args = parser.parse_args()

if args.runFinder == 1 and not args.modelPath:
    parser.error("--modelPath is required when --runFinder 1")

if args.runValidation == 0 and args.doPerfectFit == 1:
    print("WARNING: --doPerfectFit is ignored when --runValidation 0")

if args.runValidation == 1 and args.mode == 1 and args.doPerfectFit == 1:
    print("WARNING: --doPerfectFit is ignored in finder-only validation mode")

if args.runDigi == 0 and args.runFinder == 1:
    print("WARNING: --runFinder 1 with --runDigi 0 assumes digi collections are already present in the input file")

if args.runFinder == 0 and args.runFitter == 1:
    print("WARNING: --runFitter 1 with --runFinder 0 assumes finder-track collections are already present in the input file")

if args.runPerfectTracking == 0 and args.doPerfectFit == 1 and args.runValidation == 1 and args.mode in [0, 2]:
    print("WARNING: --doPerfectFit 1 with --runPerfectTracking 0 assumes PerfectFittedTracks is already present in the input file")

if args.runValidation == 1 and args.mode in [0, 1] and args.runFinder == 0:
    print("WARNING: Validation mode requires FinderTracks, but --runFinder 0. "
          "Assuming finder tracks are already present in the input file.")

if args.runValidation == 1 and args.mode in [0, 2] and args.runFitter == 0:
    print("WARNING: Validation mode requires FittedTracks, but --runFitter 0. "
          "Assuming fitted tracks are already present in the input file.")

if all(flag == 0 for flag in [args.runDigi, args.runFinder, args.runFitter,
                              args.runPerfectTracking, args.runValidation]):
    parser.error("Nothing to do: all run flags are set to 0")

# --------------------
# Fixed collection names
# --------------------
MC_COLLECTION = "MCParticles"

PLANAR_LINK_COLLECTIONS = [
    "SiWrBSimDigiLinks",
    "SiWrDSimDigiLinks",
    "VTXBSimDigiLinks",
    "VTXDSimDigiLinks",
]

DCH_LINK_COLLECTIONS = ["DCHDigiSimAssociationCollection"] if args.useDCH == 1 else []

HIT_SIM_LINK_COLLECTIONS = PLANAR_LINK_COLLECTIONS + DCH_LINK_COLLECTIONS

PLANAR_DIGI_COLLECTIONS = [
    "VTXBDigis",
    "VTXDDigis",
    "SiWrBDigis",
    "SiWrDDigis",
]

DCH_DIGI_COLLECTIONS = ["DCHDigis"] if args.useDCH == 1 else []

FINDER_TRACK_COLLECTION = "GGTFTracks"
FITTED_TRACK_COLLECTION = "Fitted_tracks_with_filtered_hits"
PERFECT_TRACK_COLLECTION = "PerfectTracks"
PERFECT_FITTED_TRACK_COLLECTION = "PerfectFittedTracks"

# --------------------
# IO
# --------------------
io = IOSvc("IOSvc")
io.Input = args.inputFile
io.Output = args.outputFile

# --------------------
# Geometry
# --------------------
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [args.compactFile]
geoservice.EnableGeant4Geo = False
geoservice.OutputLevel = INFO

# --------------------
# Algorithm sequence
# --------------------
TopAlg = []

# --------------------
# Digitizers
# --------------------
if args.runDigi == 1:
    from Configurables import DDPlanarDigi, DCHdigi_v02

    innerVertexResolution_x = 0.003
    innerVertexResolution_y = 0.003
    innerVertexResolution_t = 1000

    outerVertexResolution_x = 0.050 / math.sqrt(12)
    outerVertexResolution_y = 0.150 / math.sqrt(12)
    outerVertexResolution_t = 1000

    siWrapperResolution_x = 0.050 / math.sqrt(12)
    siWrapperResolution_y = 1.0 / math.sqrt(12)
    siWrapperResolution_t = 0.040

    vtxb_digitizer = DDPlanarDigi("VTXBdigitizer")
    vtxb_digitizer.SubDetectorName = "Vertex"
    vtxb_digitizer.IsStrip = False
    vtxb_digitizer.ResolutionU = [
        innerVertexResolution_x,
        innerVertexResolution_x,
        innerVertexResolution_x,
        outerVertexResolution_x,
        outerVertexResolution_x,
    ]
    vtxb_digitizer.ResolutionV = [
        innerVertexResolution_y,
        innerVertexResolution_y,
        innerVertexResolution_y,
        outerVertexResolution_y,
        outerVertexResolution_y,
    ]
    vtxb_digitizer.ResolutionT = [
        innerVertexResolution_t,
        innerVertexResolution_t,
        innerVertexResolution_t,
        outerVertexResolution_t,
        outerVertexResolution_t,
    ]
    vtxb_digitizer.SimTrackHitCollectionName = ["VertexBarrelCollection"]
    vtxb_digitizer.SimTrkHitRelCollection = ["VTXBSimDigiLinks"]
    vtxb_digitizer.TrackerHitCollectionName = ["VTXBDigis"]
    vtxb_digitizer.ForceHitsOntoSurface = True

    vtxd_digitizer = DDPlanarDigi("VTXDdigitizer")
    vtxd_digitizer.SubDetectorName = "Vertex"
    vtxd_digitizer.IsStrip = False
    vtxd_digitizer.ResolutionU = [
        outerVertexResolution_x,
        outerVertexResolution_x,
        outerVertexResolution_x,
    ]
    vtxd_digitizer.ResolutionV = [
        outerVertexResolution_y,
        outerVertexResolution_y,
        outerVertexResolution_y,
    ]
    vtxd_digitizer.ResolutionT = [
        outerVertexResolution_t,
        outerVertexResolution_t,
        outerVertexResolution_t,
    ]
    vtxd_digitizer.SimTrackHitCollectionName = ["VertexEndcapCollection"]
    vtxd_digitizer.SimTrkHitRelCollection = ["VTXDSimDigiLinks"]
    vtxd_digitizer.TrackerHitCollectionName = ["VTXDDigis"]
    vtxd_digitizer.ForceHitsOntoSurface = True

    siwrb_digitizer = DDPlanarDigi("SiWrBdigitizer")
    siwrb_digitizer.SubDetectorName = "SiWrB"
    siwrb_digitizer.IsStrip = False
    siwrb_digitizer.ResolutionU = [siWrapperResolution_x, siWrapperResolution_x]
    siwrb_digitizer.ResolutionV = [siWrapperResolution_y, siWrapperResolution_y]
    siwrb_digitizer.ResolutionT = [siWrapperResolution_t, siWrapperResolution_t]
    siwrb_digitizer.SimTrackHitCollectionName = ["SiWrBCollection"]
    siwrb_digitizer.SimTrkHitRelCollection = ["SiWrBSimDigiLinks"]
    siwrb_digitizer.TrackerHitCollectionName = ["SiWrBDigis"]
    siwrb_digitizer.ForceHitsOntoSurface = True

    siwrd_digitizer = DDPlanarDigi("SiWrDdigitizer")
    siwrd_digitizer.SubDetectorName = "SiWrD"
    siwrd_digitizer.IsStrip = False
    siwrd_digitizer.ResolutionU = [siWrapperResolution_x, siWrapperResolution_x]
    siwrd_digitizer.ResolutionV = [siWrapperResolution_y, siWrapperResolution_y]
    siwrd_digitizer.ResolutionT = [siWrapperResolution_t, siWrapperResolution_t]
    siwrd_digitizer.SimTrackHitCollectionName = ["SiWrDCollection"]
    siwrd_digitizer.SimTrkHitRelCollection = ["SiWrDSimDigiLinks"]
    siwrd_digitizer.TrackerHitCollectionName = ["SiWrDDigis"]
    siwrd_digitizer.ForceHitsOntoSurface = True

    TopAlg += [vtxb_digitizer, vtxd_digitizer, siwrb_digitizer, siwrd_digitizer]

    if args.useDCH == 1:
        dch_digitizer = DCHdigi_v02(
            "DCHdigi2",
            InputSimHitCollection=["DCHCollection"],
            OutputDigihitCollection=["DCH_DigiCollection"],
            OutputLinkCollection=["DCH_DigiSimAssociationCollection"],
            DCH_name="DCH_v2",
            zResolution_mm=30.0,
            xyResolution_mm=0.1,
            Deadtime_ns=400.0,
            GasType=0,
            ReadoutWindowStartTime_ns=1.0,
            ReadoutWindowDuration_ns=450.0,
            DriftVelocity_um_per_ns=-1.0,
            SignalVelocity_mm_per_ns=200.0,
            OutputLevel=INFO,
        )
        TopAlg += [dch_digitizer]

# --------------------
# Track finder
# --------------------
if args.runFinder == 1:
    from Configurables import GGTFTrackFinder

    ggtf = GGTFTrackFinder(
        "GGTFTrackFinder",
        InputPlanarHitCollections=PLANAR_DIGI_COLLECTIONS,
        InputWireHitCollections=DCH_DIGI_COLLECTIONS,
        OutputTracksGGTF=[FINDER_TRACK_COLLECTION],
        ModelPath=args.modelPath,
        Tbeta=0.6,
        Td=0.3,
        OutputLevel=INFO,
    )
    TopAlg += [ggtf]

# --------------------
# Reco fitter
# --------------------
if args.runFitter == 1:
    from Configurables import GenfitTrackFitter

    reco_fitter = GenfitTrackFitter("RecoTrackFitter")
    reco_fitter.InputTracks = [FINDER_TRACK_COLLECTION]

    #reco_fitter.OutputFittedTracks = [FITTED_TRACK_COLLECTION]
    reco_fitter.OutputFittedTracksWithFilteredHits = ["Fitted_tracks_with_filtered_hits"]
    reco_fitter.OutputFittedHits = ["FittedHits"]

    reco_fitter.RunSingleEvaluation = True
    reco_fitter.UseBrems = False
    reco_fitter.BetaInit = 100.0
    reco_fitter.BetaFinal = 0.1
    reco_fitter.BetaSteps = 15
    reco_fitter.InitializationType = 1

    reco_fitter.SkipTrackOrdering = False
    reco_fitter.ListOfTypesToSkip = [0]
    reco_fitter.FilterTrackHits = True
    reco_fitter.RunCalorimeterExtrapolation = False

    reco_fitter.OutputLevel = INFO

    TopAlg += [reco_fitter]

# --------------------
# Perfect tracking + perfect fitter
# --------------------
if args.runPerfectTracking == 1:
    from Configurables import PerfectTrackFinder, GenfitTrackFitter

    perfect = PerfectTrackFinder("PerfectTrackFinder")
    perfect.InputMCParticles = [MC_COLLECTION]
    perfect.InputPlanarHitCollections = PLANAR_LINK_COLLECTIONS
    perfect.InputWireHitCollections = DCH_LINK_COLLECTIONS
    perfect.OutputPerfectTracks = [PERFECT_TRACK_COLLECTION]
    perfect.OutputLevel = INFO

    perfect_fitter = GenfitTrackFitter("PerfectTrackFitter")
    perfect_fitter.InputTracks = [PERFECT_TRACK_COLLECTION]

    perfect_fitter.OutputFittedTracks = [PERFECT_FITTED_TRACK_COLLECTION]
    perfect_fitter.OutputFittedTracksWithFilteredHits = [PERFECT_FITTED_TRACK_COLLECTION + "_FilteredHits"]
    perfect_fitter.OutputFittedHits = ["PerfectFittedHits"]

    perfect_fitter.RunSingleEvaluation = True
    perfect_fitter.UseBrems = False
    perfect_fitter.BetaInit = 100.0
    perfect_fitter.BetaFinal = 0.1
    perfect_fitter.BetaSteps = 15
    perfect_fitter.InitializationType = 0

    perfect_fitter.SkipTrackOrdering = False
    perfect_fitter.ListOfTypesToSkip = [1]
    perfect_fitter.FilterTrackHits = True

    perfect_fitter.OutputLevel = INFO

    TopAlg += [perfect, perfect_fitter]

# --------------------
# Validation consumer
# --------------------
if args.runValidation == 1:
    from Configurables import TrackingValidation

    val = TrackingValidation("TrackingValidation")
    val.OutputFile = args.validationFile
    val.Mode = args.mode

    val.Bz = 2.0
    val.RefPointX = 0.0
    val.RefPointY = 0.0
    val.RefPointZ = 0.0

    val.DoPerfectFit = bool(args.doPerfectFit and args.mode in [0, 2])
    val.FinderEfficiencyDefinition = args.finderEfficiencyDefinition
    val.FinderPurityThreshold = args.finderPurityThreshold

    val.MCParticles = [MC_COLLECTION]
    val.HitSimLinks = HIT_SIM_LINK_COLLECTIONS

    # Mode-dependent validation inputs.
    # mode 0: full validation
    # mode 1: finder-only validation
    # mode 2: fitter-only validation
    val.FinderTracks = [FINDER_TRACK_COLLECTION] if args.mode in [0, 1] else []
    val.FittedTracks = [FITTED_TRACK_COLLECTION] if args.mode in [0, 2] else []

    val.PerfectFittedTracks = (
        [PERFECT_FITTED_TRACK_COLLECTION]
        if args.doPerfectFit == 1 and args.mode in [0, 2]
        else []
    )

    val.OutputLevel = INFO

    TopAlg += [val]

# --------------------
# AppMgr
# --------------------
ApplicationMgr(
    TopAlg=TopAlg,
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[
        geoservice,
        EventDataSvc("EventDataSvc"),
        UniqueIDGenSvc("uidSvc"),
        RndmGenSvc(),
    ],
    OutputLevel=INFO,
)