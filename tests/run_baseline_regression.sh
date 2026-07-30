#!/usr/bin/env bash
set -euo pipefail

if [[ "${WM_PROJECT_VERSION:-}" != "14" ]]; then
    echo "ERROR: source OpenFOAM Foundation v14 before running this test." >&2
    exit 1
fi

for command_name in foamDictionary foamListTimes foamMultiRun sha256sum; do
    if ! command -v "$command_name" >/dev/null 2>&1; then
        echo "ERROR: required command '$command_name' was not found." >&2
        exit 1
    fi
done

if [[ ! -x /usr/bin/time ]]; then
    echo "ERROR: GNU time is required at /usr/bin/time." >&2
    exit 1
fi

project_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
case_root="$(cd "$project_dir/../CIGMA/2Dmodel" && pwd)"

standard_source="${1:-$case_root/CCSJ01_multiphaseEuler_wallCondensation_CHT}"
modular_source="${2:-$case_root/CCSJ01_cigmaMPE_wallCondensation_CHT}"

for source_case in "$standard_source" "$modular_source"; do
    for required_path in 0 constant system system/controlDict; do
        if [[ ! -e "$source_case/$required_path" ]]; then
            echo "ERROR: missing '$source_case/$required_path'." >&2
            exit 1
        fi
    done
done

test_root="$(mktemp -d "${TMPDIR:-/tmp}/cigmaMPEFoam-baseline.XXXXXX")"

cleanup()
{
    if [[ "${KEEP_TEST_CASES:-0}" == "1" ]]; then
        echo "Temporary cases retained at: $test_root"
    else
        rm -rf -- "$test_root"
    fi
}
trap cleanup EXIT

standard_case="$test_root/standard"
modular_case="$test_root/modular"
mkdir "$standard_case" "$modular_case"

copy_case_inputs()
{
    local source_case="$1"
    local target_case="$2"

    cp -a \
        "$source_case/0" \
        "$source_case/constant" \
        "$source_case/system" \
        "$target_case/"
}

configure_one_step()
{
    local test_case="$1"

    foamDictionary -case "$test_case" system/controlDict \
        -entry startFrom -set startTime >/dev/null
    foamDictionary -case "$test_case" system/controlDict \
        -entry startTime -set 0 >/dev/null
    foamDictionary -case "$test_case" system/controlDict \
        -entry stopAt -set endTime >/dev/null
    foamDictionary -case "$test_case" system/controlDict \
        -entry endTime -set 1e-6 >/dev/null
    foamDictionary -case "$test_case" system/controlDict \
        -entry deltaT -set 1e-6 >/dev/null
    foamDictionary -case "$test_case" system/controlDict \
        -entry adjustTimeStep -set no >/dev/null
    foamDictionary -case "$test_case" system/controlDict \
        -entry writeControl -set timeStep >/dev/null
    foamDictionary -case "$test_case" system/controlDict \
        -entry writeInterval -set 1 >/dev/null
}

copy_case_inputs "$standard_source" "$standard_case"
copy_case_inputs "$modular_source" "$modular_case"
configure_one_step "$standard_case"
configure_one_step "$modular_case"

run_case()
{
    local name="$1"
    local test_case="$2"

    /usr/bin/time -v \
        foamMultiRun -case "$test_case" -noFunctionObjects \
        >"$test_root/$name.log" \
        2>"$test_root/$name.time"
}

run_case standard "$standard_case"
run_case modular "$modular_case"

standard_time="$(foamListTimes -case "$standard_case" -latestTime)"
modular_time="$(foamListTimes -case "$modular_case" -latestTime)"

if [[ "$standard_time" != "1e-06" || "$modular_time" != "1e-06" ]]; then
    echo "ERROR: one or both cases did not write time 1e-06." >&2
    exit 1
fi

standard_output="$standard_case/$standard_time"
modular_output="$modular_case/$modular_time"

find "$standard_output" -type f -printf '%P\n' | sort \
    >"$test_root/standard.files"
find "$modular_output" -type f -printf '%P\n' | sort \
    >"$test_root/modular.files"

if ! diff -u "$test_root/standard.files" "$test_root/modular.files"; then
    echo "ERROR: the written file lists differ." >&2
    exit 1
fi

normalise_field()
{
    sed \
        -e 's/libmultiphaseEulerFvModels.so/libCONDENSATION.so/g' \
        -e 's/libmultiphaseEulerFoamFvModels.so/libCONDENSATION.so/g' \
        -e 's/libcigmaMPECondensationModels.so/libCONDENSATION.so/g' \
        "$1"
}

raw_difference_count=0
normalised_difference_count=0

while IFS= read -r relative_file; do
    if ! cmp -s \
        "$standard_output/$relative_file" \
        "$modular_output/$relative_file"; then
        case "$relative_file" in
            fluid/alphat.gas|fluid/alphat.liquid)
                raw_difference_count=$((raw_difference_count + 1))
                ;;
            *)
                echo "ERROR: unexpected raw difference in '$relative_file'." >&2
                exit 1
                ;;
        esac
    fi

    if ! diff -q \
        <(normalise_field "$standard_output/$relative_file") \
        <(normalise_field "$modular_output/$relative_file") \
        >/dev/null; then
        echo "ERROR: numerical difference in '$relative_file'." >&2
        normalised_difference_count=$((normalised_difference_count + 1))
    fi
done <"$test_root/standard.files"

if [[ "$raw_difference_count" -ne 2 ]]; then
    echo "ERROR: expected two library-metadata differences, found $raw_difference_count." >&2
    exit 1
fi

if [[ "$normalised_difference_count" -ne 0 ]]; then
    exit 1
fi

normalise_log()
{
    sed -E \
        -e '/^(Exec|Time|PID|Case)[[:space:]]*:/d' \
        -e 's/Selecting solver (multiphaseEuler|cigmaMPEFluid)/Selecting solver FLUID_SOLVER/' \
        -e 's/Selecting solver (solid|cigmaMPESolid)/Selecting solver SOLID_SOLVER/' \
        -e 's/(ExecutionTime = )[0-9.]+ s  ClockTime = [0-9]+ s/\1<TIME> s  ClockTime = <TIME> s/' \
        "$1"
}

normalise_log "$test_root/standard.log" >"$test_root/standard.normalised.log"
normalise_log "$test_root/modular.log" >"$test_root/modular.normalised.log"

if ! diff -u \
    "$test_root/standard.normalised.log" \
    "$test_root/modular.normalised.log"; then
    echo "ERROR: solver traces differ after approved metadata normalisation." >&2
    exit 1
fi

if grep -q 'FOAM FATAL' "$test_root/standard.log" "$test_root/modular.log"; then
    echo "ERROR: a fatal OpenFOAM message was detected." >&2
    exit 1
fi

standard_warnings="$(grep -c 'FOAM Warning' "$test_root/standard.log" || true)"
modular_warnings="$(grep -c 'FOAM Warning' "$test_root/modular.log" || true)"

if [[ "$standard_warnings" != "$modular_warnings" ]]; then
    echo "ERROR: warning counts differ." >&2
    exit 1
fi

aggregate_hash()
{
    local output_dir="$1"

    while IFS= read -r relative_file; do
        normalise_field "$output_dir/$relative_file" \
            | sha256sum \
            | awk -v file="$relative_file" '{print $1, file}'
    done <"$test_root/standard.files" | sha256sum | awk '{print $1}'
}

standard_hash="$(aggregate_hash "$standard_output")"
modular_hash="$(aggregate_hash "$modular_output")"

standard_wall_time="$(awk -F': ' '/Elapsed \(wall clock\) time/{print $2}' "$test_root/standard.time")"
modular_wall_time="$(awk -F': ' '/Elapsed \(wall clock\) time/{print $2}' "$test_root/modular.time")"
standard_rss="$(awk -F': ' '/Maximum resident set size/{print $2}' "$test_root/standard.time")"
modular_rss="$(awk -F': ' '/Maximum resident set size/{print $2}' "$test_root/modular.time")"
written_files="$(wc -l <"$test_root/standard.files")"
openfoam_build="$(awk '$1 == "Build" && $2 == ":" {print $3; exit}' "$test_root/standard.log")"

echo "PASS: cigmaMPEFoam baseline regression"
echo "OpenFOAM build: $openfoam_build"
echo "Written files per case: $written_files"
echo "Expected raw metadata differences: $raw_difference_count"
echo "Normalised field differences: $normalised_difference_count"
echo "Normalised aggregate hash: $standard_hash"
echo "Standard aggregate hash: $standard_hash"
echo "Modular aggregate hash: $modular_hash"
echo "Standard warnings: $standard_warnings"
echo "Modular warnings: $modular_warnings"
echo "Standard wall time: $standard_wall_time"
echo "Modular wall time: $modular_wall_time"
echo "Standard peak RSS: $standard_rss kB"
echo "Modular peak RSS: $modular_rss kB"
