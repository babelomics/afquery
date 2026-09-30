#!/usr/bin/env bash

set -euo pipefail

# -----------------------------
# Usage:
# normalize_vcf.sh <vcf_id> <vcf> <ref_genome> <gender(M|F|male|female)> [threads]
#
# Example:
# normalize_vcf.sh sample1 input.vcf.gz ref.fa male 8
# -----------------------------

# -----------------------------
# Check input parameters
# -----------------------------
if [[ $# -lt 4 ]]
then
    echo "USAGE: $(basename "$0") <vcf_id> <vcf> <ref_genome> <gender(M|F|male|female)> [threads(default=4)]"
    exit 1
fi

VCF_ID="$1"
VCF="$2"
REF="$3"
GENDER="$4"
THREADS="${5:-4}"

# -----------------------------
# bcftools version (< 1.20 refuses to combine -m and -d in bcftools norm)
# -----------------------------
MIN_BCFTOOLS="1.20"
BCFTOOLS_VERSION=$(bcftools --version | head -n1 | awk '{print $2}')
if [[ "$(printf '%s\n' "${MIN_BCFTOOLS}" "${BCFTOOLS_VERSION}" | sort -V | head -n1)" != "${MIN_BCFTOOLS}" ]]
then
    echo "ERROR: bcftools >= ${MIN_BCFTOOLS} is required (found ${BCFTOOLS_VERSION})"
    exit 1
fi

OUT="${VCF_ID}_norm.vcf.gz"
CHR_MAP="chr_list.txt"
GENDER_FILE="${VCF_ID}_gender.txt"

# -----------------------------
# Chromosome rename table
# -----------------------------
if [[ ! -e "${CHR_MAP}" ]]
then
    {
        echo "chrM MT"
        for i in {1..22}; do
            echo "chr${i} ${i}"
        done
        echo "chrX X"
        echo "chrY Y"
    } > "${CHR_MAP}"
fi

# -----------------------------
# Detect sample ID
# -----------------------------
NSAMPLES=$(bcftools query -l "${VCF}" | wc -l)
ID=$(bcftools query -l "${VCF}" | head -n1)

if [[ "${NSAMPLES}" -gt 1 ]]
then
    echo "WARNING: VCF contains more than one sample. Using first sample: ${ID}"
fi

# -----------------------------
# Normalize gender input
# -----------------------------
case "${GENDER,,}" in
    m|male)
        GENDER="M"
        ;;
    f|female)
        GENDER="F"
        ;;
    *)
        echo "ERROR: gender must be M/F or male/female"
        exit 1
        ;;
esac

echo "${ID} ${GENDER}" > "${GENDER_FILE}"

# -----------------------------
# Chromosome targets
# -----------------------------
TARGETS=$(printf "%s," {1..22} X Y MT)
TARGETS="${TARGETS%,}"

# -----------------------------
# Normalize VCF
# -----------------------------
bcftools annotate \
    -x INFO,^FORMAT/GT,^FORMAT/DP,^FORMAT/GQ \
    --force \
    --rename-chrs "${CHR_MAP}" \
    "${VCF}" | \
bcftools norm \
    -m -both \
    -d exact \
    -f "${REF}" \
    --check-ref ws \
    --targets "${TARGETS}" | \
bcftools +setGT -- \
    -n . \
    -i 'FILTER!="PASS"' \
    -t q | \
bcftools +fixploidy -- \
    -s "${GENDER_FILE}" | \
bcftools view \
    -e 'GT="0" || GT="0/0"' \
    --threads "${THREADS}" \
    -Oz \
    -o "${OUT}"

# -----------------------------
# Index
# -----------------------------
bcftools index --threads "${THREADS}" "${OUT}"

echo "Output:"
echo "  VCF: ${OUT}"
echo "  Index: ${OUT}.csi"
