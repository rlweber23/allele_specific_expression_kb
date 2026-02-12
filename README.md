# Allele-Specific Expression KB

Nextflow pipelines for allele-specific expression (ASE) mapping of mouse snRNA-seq data from the IGVF consortium. Reads are mapped using kb-python against strain-specific concatenated genomes, followed by ambient RNA removal with CellBender.

---

## Repository Structure

```
.
├── pipeline.nf             # Main ASE mapping pipeline
├── nextflow.config         # Config for main pipeline
├── ASE_samplesheets/       # Sample sheet inputs
├── references/             # Barcode whitelists
├── scripts/                # Downstream processing scripts
└── make_references/        # Reference genome generation pipeline
    ├── pipeline.nf
    ├── nextflow.config
    └── references/         # Generated reference outputs
```

---

## Pipelines

### 1. `make_references/` — Reference Generation *(run first)*

Downloads genome, annotation, and SNP files and builds strain-specific concatenated references for allele-specific mapping.

**Downloads:**
- Mouse SNP VCF: `https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-2112-v8-SNPs_Indels/mgp_REL2021_snps.vcf.gz`
- Reference FASTA (IGVF): `IGVFFI9282QLXO.fasta.gz`
- Reference GTF (IGVF): `IGVFFI4777RDZK.gtf.gz`

**Outputs:** Strain-specific FASTA/GTF files and kb-python indices (e.g., under `references/A_J/`)

```bash
cd make_references
nextflow run pipeline.nf -c nextflow.config
```

### 2. `pipeline.nf` — ASE Mapping

Maps FASTQ files against the strain-specific references generated above using kb-python, then runs CellBender for ambient RNA removal.

```bash
nextflow run pipeline.nf -c nextflow.config
```

---

## Dependencies

| Package | Version |
|---|---|
| `kb-python` | 0.28.2 |
| `cellbender` | 0.3.2 |
| `g2gtools` | 0.2.7 |

---

## Input

Sample sheets are provided in CSV format. See `ASE_samplesheets/test_sample_sheet.csv` for formatting.
