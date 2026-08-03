#!/usr/bin/env bash
set -euo pipefail

if [[ "${WM_PROJECT_VERSION:-}" != "14" ]]; then
    echo "ERROR: source OpenFOAM Foundation v14 before running this test." >&2
    exit 1
fi

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
tutorial="$WM_PROJECT_DIR/tutorials/multicomponentFluid/verticalChannel"
test_root="$(mktemp -d /tmp/cigma-single-phase-smoke.XXXXXX)"

prepare_case()
{
    local case_dir="$1"
    local solver_name="$2"

    mkdir -p "$case_dir"
    cp -a "$tutorial/." "$case_dir/"
    cp "$case_dir/0/U.orig" "$case_dir/0/U"

    foamDictionary "$case_dir/system/controlDict" \
        -entry solver -set "$solver_name" >/dev/null
    foamDictionary "$case_dir/system/controlDict" \
        -entry startFrom -set startTime >/dev/null
    foamDictionary "$case_dir/system/controlDict" \
        -entry endTime -set 1e-5 >/dev/null
    foamDictionary "$case_dir/system/controlDict" \
        -entry writeInterval -set 1e-5 >/dev/null

    blockMesh -case "$case_dir" >"$case_dir/log.blockMesh" 2>&1
}

run_case()
{
    local case_dir="$1"
    foamRun -case "$case_dir" >"$case_dir/log.foamRun" 2>&1
    grep -q '^End$' "$case_dir/log.foamRun"
}

reference_case="$test_root/reference"
cigma_case="$test_root/cigma"
condensation_case="$test_root/condensation"

prepare_case "$reference_case" multicomponentFluid
prepare_case "$cigma_case" cigmaSinglePhaseFluid
run_case "$reference_case"
run_case "$cigma_case"

for field in H2O T U air alphat k nut omega p
do
    cmp "$reference_case/1e-05/$field" "$cigma_case/1e-05/$field"
done

prepare_case "$condensation_case" cigmaSinglePhaseFluid

foamDictionary "$condensation_case/0/H2O" \
    -entry boundaryField/walls/type -set saturatedSteam >/dev/null
foamDictionary "$condensation_case/0/H2O" \
    -entry boundaryField/walls/value -set 'uniform 0.01' >/dev/null
foamDictionary "$condensation_case/0/U" \
    -entry boundaryField/walls/type -set condensingWallVelocity >/dev/null
foamDictionary "$condensation_case/0/U" \
    -entry boundaryField/walls/value -set 'uniform (0 0 0)' >/dev/null
foamDictionary "$condensation_case/0/T" \
    -entry boundaryField/walls -set \
    '{ type externalCondensationTemperature; Ta constant 273.15; h uniform 4000; relaxation 1; qcond qcond; qcondRelaxation 1; value uniform 273.15; }' \
    >/dev/null

cp "$project_dir/tests/singlePhase/condensationProperties" \
    "$condensation_case/constant/condensationProperties"

run_case "$condensation_case"

grep -q 'Single-phase wall condensation is ON' \
    "$condensation_case/log.foamRun"
grep -q 'Updating single-phase diffusion-layer wall condensation' \
    "$condensation_case/log.foamRun"
test -f "$condensation_case/1e-05/massTransferRate"
test -f "$condensation_case/1e-05/qcond"
grep -q '^-' "$condensation_case/1e-05/massTransferRate"

echo "PASS: no-condensation fields match the native multicomponentFluid module."
echo "PASS: diffusionLayer and saturatedSteam completed an active one-step test."
echo "PASS: qcond was coupled through externalCondensationTemperature."
echo "Test cases retained at: $test_root"
