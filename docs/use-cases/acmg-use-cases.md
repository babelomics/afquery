# ACMG Criteria: BA1, PM2, and PS4

AFQuery enables evidence-based ACMG/AMP variant classification using allele frequencies computed on your own cohort. This page shows how to apply the three AF-dependent ACMG criteria — BA1, PM2, and PS4 — with AFQuery commands and Python examples.

!!! note "AFQuery supplements population databases"
    Local cohort AF provides complementary evidence to reference population databases like gnomAD. Use both: gnomAD for global population context, AFQuery for local cohort context.

---

## BA1 — Stand-Alone Benign

**Definition**: Allele frequency > 5% in any population → variant is benign (stand-alone evidence).

### When Local AF Matters for BA1

- Your cohort contains a population underrepresented in gnomAD
- You need to confirm gnomAD AF in your local sequencing pipeline
- A variant shows discrepant AF between gnomAD and your cohort

AFQuery Command:

```bash
afquery query \
  --db ./db/ \
  --locus chr7:117559590 \
  --ref C --alt T \
  --phenotype controls
```

Example output:

```
chr7:117559590 C>T  AC=160  AN=2000  AF=0.0800  n_eligible=1000  N_HET=140  N_HOM_ALT=10  N_HOM_REF=850  N_FAIL=0  N_NO_COVERAGE=0
```

**Interpretation**: AF=0.08 (8%) with AN=2000 → BA1 criterion met. This variant is benign in the controls subgroup of your cohort.

**Minimum AN for BA1**: An AF estimate is only reliable with sufficient AN. With AN=20, two carriers would give AF=0.10 — this is noise, not evidence.


---

## PM2 — Supporting Pathogenic

**Definition**: Absent or extremely low frequency in population databases → supporting evidence for pathogenicity.

!!! important "PM2 strength update (ClinGen SVI, 2020)"
    The ClinGen Sequence Variant Interpretation Working Group [recommends](https://clinicalgenome.org/site/assets/files/5182/pm2_-_svi_recommendation_-_approved_sept2020.pdf) applying PM2 at **supporting strength** (PM2_Supporting) rather than the original moderate level. Most ClinGen Variant Curation Expert Panels have adopted this downgrade. When using AFQuery results to apply PM2, treat the evidence as supporting unless your laboratory's classification framework specifies otherwise.

### Distinguishing "Absent" from "Unknown"

This is the most critical distinction when applying PM2:

| Result | Meaning | PM2 Applicable? |
|--------|---------|-----------------|
| AC=0, AN=2000 | Variant absent in well-covered cohort | Yes |
| AC=0, AN=50 | Low coverage — absence is not meaningful | No |
| AC=0, AN=0 | Position not covered by any eligible sample | No |

### Python Example: PM2 with AN Validation

```python
from afquery import Database

db = Database("./db/")
results = db.query("chr12", pos=49416565, alt="A")

if not results:
    print("Variant not in database — cannot apply PM2 (no coverage data)")
else:
    r = results[0]
    if r.AN >= 2000 and r.AC == 0:
        print(f"PM2 applicable: absent in cohort (AN={r.AN})")
    elif r.AN >= 2000 and r.AF < 0.0001:
        print(f"PM2 applicable: extremely rare (AF={r.AF:.6f}, AN={r.AN})")
    elif r.AN < 2000:
        print(f"PM2 not reliable: insufficient AN ({r.AN})")
```

---

## PS4 — Strong Pathogenic (Case Enrichment)

**Definition**: Significantly increased prevalence in affected individuals compared to controls → strong evidence for pathogenicity.

### Pseudo-Control Comparison

AFQuery's phenotype filtering enables direct comparison between cases and controls without separate databases:

```bash
# AF in affected samples (tagged with the disease phenotype)
afquery query \
  --db ./db/ \
  --locus chr2:166845670 \
  --ref G --alt A \
  --phenotype rare_disease

# AF in unaffected samples (everything except rare_disease)
afquery query \
  --db ./db/ \
  --locus chr2:166845670 \
  --ref G --alt A \
  --phenotype ^rare_disease
```

### Enrichment Ratio Calculation

```python
from afquery import Database

db = Database("./db/")

# Query affected group
cases = db.query("chr2", pos=166845670, alt="A", phenotype=["rare_disease"])
# Query unaffected group (pseudo-controls)
controls = db.query("chr2", pos=166845670, alt="A", phenotype=["^rare_disease"])

if cases and controls:
    c = cases[0]
    ctrl = controls[0]

    if ctrl.AN >= 500 and c.AN >= 100:
        case_af = c.AF if c.AF else 0.0
        ctrl_af = ctrl.AF if ctrl.AF else 0.0

        if ctrl_af > 0:
            enrichment = case_af / ctrl_af
            print(f"Case AF:    {case_af:.4f} (AN={c.AN})")
            print(f"Control AF: {ctrl_af:.4f} (AN={ctrl.AN})")
            print(f"Enrichment: {enrichment:.1f}x")
        else:
            print(f"Variant absent in controls (AN={ctrl.AN}), present in cases (AF={case_af:.4f})")
            print("Strong enrichment — PS4 applicable if case count is sufficient")
```

!!! tip "Verify carriers"
    Use `afquery variant-info` to list the specific samples carrying the variant and confirm their phenotype assignments match expectations. See [Variant Info](../guides/variant-info.md).

!!! warning "Statistical considerations"
    The original ACMG/AMP framework (Richards et al., 2015) does not specify a minimum odds ratio for PS4. The OR > 5.0 threshold used here is a widely adopted convention, but ClinGen expert panels have proposed different thresholds depending on the gene/disease context (e.g., OR ≥ 6 for strong, OR ≥ 3 for moderate in hearing loss). The enrichment ratio above is a screening tool. For formal PS4 application, compute a Fisher's exact test on the 2×2 table of (carrier/non-carrier) × (case/control) and apply the thresholds relevant to your clinical context.

For a full pseudo-control workflow, see [Pseudo-controls](../use-cases/pseudo-controls.md).

---

## AN Threshold Guidance

| Criterion | Minimum AN | Recommended AN | Rationale |
|-----------|-----------|---------------|-----------|
| BA1 | ≥500 | ≥1000 | 5% AF requires reliable denominator; AN=500 gives 95% CI of ±1.9% |
| PM2 | ≥1000 | ≥2000 | Absence is only meaningful with sufficient power to detect rare variants |
| PS4 | ≥100 per group | ≥500 per group | Enrichment tests need adequate sample size in both arms |

---

## Worked Example: Evaluating One Variant

Consider variant **chr15:48762884 C>T** in a cohort of 2500 samples:

### Step 1: Check BA1

```bash
afquery query --db ./db/ --locus chr15:48762884 --ref C --alt T
# Result: AC=3, AN=4800, AF=0.000625
```

AF=0.000625 < 0.05 → **BA1 not met** (variant is not common enough to be stand-alone benign).

### Step 2: Check PM2

AC=3, AF=0.000625 → the variant is present, so strict PM2 (absent) does not apply. However, AF < 0.001 with AN=4800 means the variant is extremely rare — PM2_Supporting may apply depending on disease prevalence and inheritance model (see PM2 strength note above).

### Step 3: Check PS4

```bash
# Cases: 300 samples tagged with the disease
afquery query --db ./db/ --locus chr15:48762884 --ref C --alt T --phenotype DISEASE_X
# chr15:48762884 C>T  AC=3  AN=580  AF=0.005172  n_eligible=290  N_HET=3  N_HOM_ALT=0  N_HOM_REF=287  N_FAIL=0  N_NO_COVERAGE=0

# Controls: remaining 2200 samples
afquery query --db ./db/ --locus chr15:48762884 --ref C --alt T --phenotype ^DISEASE_X
# chr15:48762884 C>T  AC=0  AN=4220  AF=0.000000  n_eligible=2110  N_HET=0  N_HOM_ALT=0  N_HOM_REF=2110  N_FAIL=0  N_NO_COVERAGE=0
```

All 3 carriers are in the disease group. Enrichment is significant (AF=0.005 vs AF=0.000, Fisher's exact p < 0.001). **PS4 applicable** — the variant is enriched in affected individuals.

### Summary

| Criterion | Result | Evidence |
|-----------|--------|----------|
| BA1 | Not met | AF=0.000625, well below 5% |
| PM2 | Borderline | Present but extremely rare (AF < 0.001) |
| PS4 | Met | All carriers are cases; Fisher's exact p < 0.001 |

---

## Common Pitfalls


### Low AN Masking as Absence

AC=0 with AN=50 does **not** mean the variant is absent — it means you have no statistical power. Always check AN before interpreting AC=0 as evidence for PM2.

### chrX Ploidy Effects

On chrX non-PAR regions, males contribute AN=1 and females AN=2. A cohort with 80% males has lower AN than expected from sample count alone, and hemizygous males who carry the variant contribute AC=1 as homozygotes. See [Ploidy & Sex Chromosomes](../advanced/ploidy-and-sex-chroms.md).

### High N_FAIL at the Variant Site

`N_FAIL` counts eligible samples whose call at the site failed quality filters (`FILTER≠PASS`). These samples are left out of both AC and AN, so AF is computed over the samples with a usable call. A high N_FAIL means that AF rests on fewer samples and that the site may have systematic QC problems, which warrants caution before applying ACMG criteria.

| N_FAIL / n_eligible | Interpretation |
|---|---|
| < 5% | Isolated low-quality calls — typically not concerning |
| 5–15% | Elevated failure rate — review site-level QC metrics |
| > 15% | Systematic artifact likely — treat AF with caution; do not apply BA1/PM2 without manual site inspection |

Check `N_FAIL` in the query output or use `AFQUERY_N_FAIL` in annotated VCFs:

```bash
# Flag sites with high failure rate after annotation
bcftools filter -i 'AFQUERY_N_FAIL > 0' annotated.vcf | head
```

---

## Next Steps

- [Variant Info](../guides/variant-info.md) — list carriers with metadata to verify phenotype assignments
- [Clinical Prioritization](../use-cases/clinical-prioritization.md) — full annotation and filtering workflow
- [Pseudo-controls](../use-cases/pseudo-controls.md) — case vs. control AF comparison
- [Population-Specific AF](../use-cases/population-specific-af.md) — ancestry-stratified queries
- [Understanding Output](../getting-started/understanding-output.md) — what each output field means
