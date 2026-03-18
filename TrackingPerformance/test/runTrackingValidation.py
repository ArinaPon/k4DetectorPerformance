import os
import math

from Gaudi.Configuration import INFO
from Configurables import EventDataSvc, GeoSvc, UniqueIDGenSvc, RndmGenSvc
from k4FWCore import IOSvc, ApplicationMgr
from k4FWCore.parseArgs import parser

# --------------------
# Arguments
# --------------------
parser.add_argument("--inputFile", required=True, help="Input simulated EDM4hep ROOT file")
parser.add_argument("--modelPath", required=True, help="Path to the GGTF ONNX model")
parser.add_argument("--outputFile", default="out_reco.root",
                    help="Output EDM4hep ROOT file with reconstructed collections")
parser.add_argument("--validationFile", default="validation.root",
                    help="Output ROOT file written by TrackingValidationConsumer")
args = parser.parse_args()

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
geoservice.detectors = [
    os.path.join(os.environ["K4GEO"], "FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml")
]
geoservice.EnableGeant4Geo = False
geoservice.OutputLevel = INFO

# --------------------
# Digitizers
# --------------------
from Configurables import DDPlanarDigi, DCHdigi_v02

innerVertexResolution_x = 0.003
innerVertexResolution_y = 0.003
innerVertexResolution_t = 1000

outerVertexResolution_x = 0.050 / math.sqrt(12)
outerVertexResolution_y = 0.150 / math.sqrt(12)
outerVertexResolution_t = 1000

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

siWrapperResolution_x = 0.050 / math.sqrt(12)
siWrapperResolution_y = 1.0 / math.sqrt(12)
siWrapperResolution_t = 0.040

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

# --------------------
# Track finder
# --------------------
from Configurables import GGTFTrackFinder


ggtf = GGTFTrackFinder(
    "GGTFTrackFinder",
    InputPlanarHitCollections=["VTXBDigis", "VTXDDigis", "SiWrBDigis", "SiWrDDigis"],
    InputWireHitCollections=["DCH_DigiCollection"],
    OutputTracksGGTF=["GGTFTracks"],
    ModelPath=args.modelPath,
    Tbeta=0.6,
    Td=0.3,
    OutputLevel=INFO,
)

# --------------------
# Reco fitter
# --------------------
from Configurables import GenfitTrackFitter

reco_fitter = GenfitTrackFitter("RecoTrackFitter")
reco_fitter.InputTracks = ["GGTFTracks"]
reco_fitter.OutputFittedTracks = ["FittedTracks"]
reco_fitter.RunSingleEvaluation = True
reco_fitter.UseBrems = False
reco_fitter.BetaInit = 100.0
reco_fitter.BetaFinal = 0.1
reco_fitter.BetaSteps = 10
reco_fitter.InitializationType = 1
reco_fitter.SkipTrackOrdering = False
reco_fitter.SkipUnmatchedTracks = False
reco_fitter.OutputLevel = INFO

# --------------------
# Perfect track finder
# --------------------
from Configurables import PerfectTrackFinder

perfect = PerfectTrackFinder("PerfectTrackFinder")
perfect.InputMCParticles = ["MCParticles"]
perfect.InputPlanarHitCollections = [
    "SiWrBSimDigiLinks",
    "SiWrDSimDigiLinks",
    "VTXBSimDigiLinks",
    "VTXDSimDigiLinks",
]
perfect.InputWireHitCollections = ["DCH_DigiSimAssociationCollection"]
perfect.OutputPerfectTracks = ["PerfectTracks"]
perfect.OutputLevel = INFO

# --------------------
# Perfect fitter
# --------------------
perfect_fitter = GenfitTrackFitter("PerfectTrackFitter")
perfect_fitter.InputTracks = ["PerfectTracks"]
perfect_fitter.OutputFittedTracks = ["PerfectFittedTracks"]
perfect_fitter.RunSingleEvaluation = True
perfect_fitter.UseBrems = False
perfect_fitter.BetaInit = 100.0
perfect_fitter.BetaFinal = 0.1
perfect_fitter.BetaSteps = 10
perfect_fitter.InitializationType = 1
perfect_fitter.SkipTrackOrdering = False
perfect_fitter.SkipUnmatchedTracks = False
perfect_fitter.OutputLevel = INFO

# --------------------
# Validation consumer
# --------------------
from Configurables import TrackingValidation

val = TrackingValidation("TrackingValidation")
val.OutputFile = args.validationFile
val.Mode = 0
val.Bz = 2.0
val.RefPointX = 0.0
val.RefPointY = 0.0
val.RefPointZ = 0.0
val.DoPerfectFit = True
val.MCParticles = ["MCParticles"]
val.PlanarLinks = [
    "SiWrBSimDigiLinks",
    "SiWrDSimDigiLinks",
    "VTXBSimDigiLinks",
    "VTXDSimDigiLinks",
]
val.DCHLinks = ["DCH_DigiSimAssociationCollection"]
val.FinderTracks = ["GGTFTracks"]
val.FittedTracks = ["FittedTracks"]
val.PerfectFittedTracks = ["PerfectFittedTracks"]
val.OutputLevel = INFO

# --------------------
# AppMgr
# --------------------
ApplicationMgr(
    TopAlg=[
        dch_digitizer,
        vtxb_digitizer,
        vtxd_digitizer,
        siwrb_digitizer,
        siwrd_digitizer,
        ggtf,
        reco_fitter,
        perfect,
        perfect_fitter,
        val,
    ],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[geoservice, EventDataSvc("EventDataSvc"), UniqueIDGenSvc("uidSvc"), RndmGenSvc()],
    OutputLevel=INFO,
)
