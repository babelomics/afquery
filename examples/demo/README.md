# AFQuery Demo Dataset

This directory contains a script to generate synthetic demo data for the AFQuery tutorial.

## Usage

```bash
python create_demo_data.py
```

This creates a `demo_output/` directory with:

- **10 single-sample VCFs** (`demo_output/vcfs/`) — synthetic genotype data across autosomes, chrX, chrY, and chrM. Genotypes are haploid (`0`/`1`) where the sample carries one copy: males on chrX (outside the PARs) and chrY, and every sample on chrM
- **2 BED files** (`demo_output/beds/`) — capture regions for `wes_v1` and `wes_v2` technologies
- **1 manifest** (`demo_output/manifest.tsv`) — sample metadata with sex, technology, and phenotype codes

## Sample Composition

| Samples | Sex | Technology | Phenotype codes |
|---------|-----|------------|-----------------|
| DEMO_001–DEMO_004 | 2F, 2M | WGS | E11.9, I10, control |
| DEMO_005–DEMO_007 | 2F, 1M | wes_v1 | E11.9, I10, control |
| DEMO_008–DEMO_010 | 1F, 2M | wes_v2 | E11.9, I10, control |

## Next Steps

After generating the data, follow the [Tutorial](../../docs/getting-started/tutorial.md) for a complete walkthrough.
