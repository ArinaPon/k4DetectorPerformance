#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
MODEL_FILE="${1:-}"

if [ -z "${MODEL_FILE}" ]; then
  echo "ERROR: missing ONNX model path argument"
  echo "Usage: $0 /full/path/to/SimpleGatrIDEAv3o1.onnx"
  exit 1
fi

if [ ! -f "${MODEL_FILE}" ]; then
  echo "ERROR: ONNX model file not found: ${MODEL_FILE}"
  exit 1
fi

if [ -z "${K4GEO:-}" ]; then
  echo "ERROR: K4GEO is not set"
  exit 1
fi

if [ ! -d "${K4GEO}" ]; then
  echo "ERROR: K4GEO does not point to a valid directory: ${K4GEO}"
  exit 1
fi

if ! command -v ddsim >/dev/null 2>&1; then
  echo "ERROR: ddsim not found in PATH"
  exit 1
fi

if ! command -v k4run >/dev/null 2>&1; then
  echo "ERROR: k4run not found in PATH"
  exit 1
fi

if ! command -v curl >/dev/null 2>&1; then
  echo "ERROR: curl not found in PATH"
  exit 1
fi

XML_FILE="${K4GEO}/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml"
RUN_FILE="${SCRIPT_DIR}/runTrackingValidation.py"
VAL_FILE="${SCRIPT_DIR}/validation_output_test.root"

if [ ! -f "${XML_FILE}" ]; then
  echo "ERROR: geometry XML file not found: ${XML_FILE}"
  exit 1
fi

if [ ! -f "${RUN_FILE}" ]; then
  echo "ERROR: tracking validation run file not found: ${RUN_FILE}"
  exit 1
fi

TMPDIR="$(mktemp -d "${TMPDIR:-/tmp}/tracking_validation.XXXXXX")"
trap 'rm -rf "${TMPDIR}"' EXIT

STEERING_FILE="${TMPDIR}/SteeringFile_IDEA_o1_v03.py"
SIM_FILE="${TMPDIR}/out_sim_edm4hep.root"
RECO_FILE="${TMPDIR}/out_reco.root"

rm -f "${VAL_FILE}"

N_EVENTS=5
SEED=42

echo "=== Downloading DDSim steering file ==="
curl -L \
  -o "${STEERING_FILE}" \
  https://raw.githubusercontent.com/key4hep/k4geo/master/example/SteeringFile_IDEA_o1_v03.py

if [ ! -f "${STEERING_FILE}" ]; then
  echo "ERROR: failed to download DDSim steering file"
  exit 1
fi

echo "=== Test configuration ==="
echo "Script dir:                 ${SCRIPT_DIR}"
echo "Temporary dir:              ${TMPDIR}"
echo "Geometry XML:               ${XML_FILE}"
echo "DDSim steering:             ${STEERING_FILE}"
echo "Run script:                 ${RUN_FILE}"
echo "ONNX model:                 ${MODEL_FILE}"
echo "Simulation file:            ${SIM_FILE}"
echo "Reco file:                  ${RECO_FILE}"
echo "Validation file:            ${VAL_FILE}"
echo "Events:                     ${N_EVENTS}"
echo "Seed:                       ${SEED}"
echo "runDigi:                    1"
echo "runFinder:                  1"
echo "runFitter:                  1"
echo "runPerfectTracking:         1"
echo "runValidation:              1"
echo "useDCH:                     1"
echo "mode:                       0"
echo "doPerfectFit:               1"
echo "finderEfficiencyDefinition: 1"
echo "finderPurityThreshold:      0.75"

echo "=== Step 1: DDSim ==="
ddsim \
  --steeringFile "${STEERING_FILE}" \
  --compactFile "${XML_FILE}" \
  -G \
  --gun.particle mu- \
  --gun.energy "5*GeV" \
  --gun.distribution uniform \
  --gun.thetaMin "89*deg" \
  --gun.thetaMax "89*deg" \
  --gun.phiMin "0*deg" \
  --gun.phiMax "0*deg" \
  --random.seed "${SEED}" \
  --numberOfEvents "${N_EVENTS}" \
  --outputFile "${SIM_FILE}"

if [ ! -f "${SIM_FILE}" ]; then
  echo "ERROR: simulation output was not created: ${SIM_FILE}"
  exit 1
fi

echo "=== Step 2: digitization + tracking + validation ==="
k4run "${RUN_FILE}" \
  --inputFile "${SIM_FILE}" \
  --modelPath "${MODEL_FILE}" \
  --outputFile "${RECO_FILE}" \
  --validationFile "${VAL_FILE}" \
  --geom "${XML_FILE}" \
  --runDigi 1 \
  --runFinder 1 \
  --runFitter 1 \
  --runPerfectTracking 1 \
  --runValidation 1 \
  --useDCH 1 \
  --mode 0 \
  --doPerfectFit 1 \
  --finderEfficiencyDefinition 1 \
  --finderPurityThreshold 0.75

if [ ! -f "${RECO_FILE}" ]; then
  echo "ERROR: reconstruction output was not created: ${RECO_FILE}"
  exit 1
fi

if [ ! -f "${VAL_FILE}" ]; then
  echo "ERROR: validation output was not created: ${VAL_FILE}"
  exit 1
fi

echo "=== Step 3: check outputs ==="
test -f "${SIM_FILE}"
test -f "${RECO_FILE}"
test -f "${VAL_FILE}"

echo "Test completed successfully."
echo "Validation file:  ${VAL_FILE}"