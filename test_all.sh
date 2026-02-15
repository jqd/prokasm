#!/bin/bash
#
# Test script for ProkASM pipeline
# Tests all four data types: Illumina, Nanopore, PacBio, and Hybrid
#

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}╔════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║          ProkASM Pipeline - Complete Test Suite           ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Track test results
TESTS_PASSED=0
TESTS_FAILED=0
FAILED_TESTS=()

# Function to run a test
run_test() {
    local test_name=$1
    shift
    local cmd="$@"
    
    echo -e "\n${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${YELLOW}Testing: ${test_name}${NC}"
    echo -e "${YELLOW}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "Command: ${cmd}"
    echo ""
    
    if eval $cmd; then
        echo -e "\n${GREEN}✓ ${test_name} - PASSED${NC}"
        ((TESTS_PASSED++))
    else
        echo -e "\n${RED}✗ ${test_name} - FAILED${NC}"
        ((TESTS_FAILED++))
        FAILED_TESTS+=("$test_name")
    fi
}

# Test 1: Illumina workflow
run_test "Illumina Workflow" \
    "nextflow run main.nf \
        --data_type illumina \
        --samplesheet illumina.csv \
        --outdir test_results/illumina \
        -profile cluster \
        -resume"

# Test 2: Nanopore workflow
run_test "Nanopore Workflow" \
    "nextflow run main.nf \
        --data_type nanopore \
        --samplesheet nanopore.csv \
        --longread_format fastq \
        --flye_mode nano-hq \
        --outdir test_results/nanopore \
        -profile cluster \
        -resume"

# Test 3: PacBio workflow
run_test "PacBio Workflow" \
    "nextflow run main.nf \
        --data_type pacbio \
        --samplesheet pacbio.csv \
        --longread_format fastq \
        --pacbio_assembler hifiasm \
        --outdir test_results/pacbio \
        -profile cluster \
        -resume"

# Test 4: Hybrid workflow
run_test "Hybrid Workflow" \
    "nextflow run main.nf \
        --data_type hybrid \
        --samplesheet hybrid.csv \
        --longread_type nanopore \
        --longread_format fastq \
        --outdir test_results/hybrid \
        -profile cluster \
        -resume"

# Summary
echo -e "\n${BLUE}╔════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║                      Test Summary                          ║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "Total Tests: $((TESTS_PASSED + TESTS_FAILED))"
echo -e "${GREEN}Passed: ${TESTS_PASSED}${NC}"
echo -e "${RED}Failed: ${TESTS_FAILED}${NC}"
echo ""

if [ $TESTS_FAILED -gt 0 ]; then
    echo -e "${RED}Failed tests:${NC}"
    for test in "${FAILED_TESTS[@]}"; do
        echo -e "  ${RED}✗${NC} $test"
    done
    echo ""
    exit 1
else
    echo -e "${GREEN}✓ All tests passed successfully!${NC}"
    echo ""
    exit 0
fi
