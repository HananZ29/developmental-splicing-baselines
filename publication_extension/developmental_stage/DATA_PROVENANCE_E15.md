# E15 data provenance

All files below are public and were downloaded on 2026-11-05 with the user's approval. Nothing here
is an outcome: no sQTL ground truth has been downloaded or read for E15.

| Local path (`data/external/…`) | Source URL | Size (B) | MD5 / check | Status |
|---|---|---|---|---|
| `brainspan_pec/meta/mRNA-seq_Sample_metadata.xlsx` | http://development.psychencode.org/files/raw_data/mRNA-seq_Sample%20metadata.xlsx | 31,834 | 54bcf31901aab52651f4a78fb8c44014 | done |
| `brainspan_pec/meta/mRNA-seq_QC.xlsx` | http://development.psychencode.org/files/processed_data/RNA-seq/mRNA-seq_QC.xlsx | 213,488 | 0240a62bf4ddcf8327903af906f37860 | done |
| `brainspan_pec/psi/brainSpan.hg38.gencode21.STAR.DEXSeq.exonPSI.txt` | http://development.psychencode.org/files/processed_data/RNA-seq/brainSpan.hg38.gencode21.STAR.DEXSeq.exonPSI.txt | 1,206,684,334 expected | size checked on completion | downloading |
| `gtex_junctions/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz` | https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz | 4,245,392,642 expected | size checked on completion | downloading |
| `brainspan_mrf/neocortex/*.zipentry` + `manifest.tsv` | http://download.alleninstitute.org/brainspan/MRF_BigWig_Gencode_v10/mrf/<DONOR>.zip, neocortex entries only, fetched by byte range with `scripts/E15_fetch_brainspan_mrf.py` | 20 GB total; 362 entries, 42 donors (completed 2026-11-06) | Byte length checked per entry: all 362 match the manifest. CRC-32 from the zip central directory is checked on decompression (first entry, HSB_103_A1C: CRC OK, uncompressed size OK) | done |

## Sample metadata (from `mRNA-seq_Sample_metadata.xlsx`)

**41 donors with PSI-matrix samples:**
- **18 prenatal** (8–37 pcw):
  - HSB112, 148, 153, 150, 113, 103, 149, 114, 178, 154, 96, 97, 98, 107, 92, 159;
  - two late-fetal donors: HSB155 (35 pcw) and HSB194 (37 pcw).
- **23 postnatal** (4 months–40 years). Of these, **8 are ≥ 18 years**: HSB105, 127, 130, 136, 126, 145, 123, 135.
- Covariates per donor: sex, ethnicity, PMI, pH, RIN, hemisphere and **sequencing site (USC or Yale)**.
  - Sequencing site is a potential batch factor. It must be checked for balance between the fetal
    and adult baseline groups before E15 is pre-registered.

**Excluded from any fetal baseline evaluated against Wen sQTLs** (they are also in the Wen cohort):
HSB107, HSB148, HSB150, HSB154, HSB194.

**MRF archives vs metadata:**
- There are 42 MRF donor archives.
- HSB168 and HSB169 have MRF archives but are not in the PSI/metadata table, so their ages are
  unknown here. They are not used unless their ages are documented.
- HSB131 is in the metadata but has no MRF archive.

**Region codes in the MRF entry names** (checked): adult A1C, DFC, IPC, ITC, M1C, MFC, OFC, S1C,
STC, V1C, VFC; early prenatal M1C-S1C, Ocx, PCx.

## Not obtained

- **Wen et al. 2024 sQTL files** (Synapse syn50897018): need a Synapse login by the user.

  **Synapse `file_description` annotations, read via the API on 2026-11-05:**

  | ID | File | Annotation |
  |---|---|---|
  | syn51183885 | `mixed_cissqtl_40hcp_perm_purity_filtered.txt.gz` | "SuSiE fine-mapping results of cis-sQTL" (full multi-ancestry cohort) |
  | syn51183855 | `conditional_top_variants.txt` | "Sum stats of conditional cis-sQTL mapping in the full multi-ancestry dataset" |
  | syn51183854 | `leafcutter_clusters_to_genes.txt` | "Intron annotation from GTEx script" |
  | syn51183859 | `all.introns.tested.tsv` | "Intron annotation from customized leafviz" |
  | syn51184552 | `mixed_nominal_sqtl_40HCP_all_assoc.txt.gz` | Nominal, full cohort. Optional; only for the cross-stage lookup |
  | syn51177112, syn51177113 | Tri1/Tri2 permutation files | "…permutation pass for the **EUR** Tri1 / Tri2 sQTL". The trimester mapping used EUR-ancestry samples |
  | syn51183886, syn51183884 | Tri1/Tri2 `perm_purity_filtered` | SuSiE fine-mapping in Tri1 / Tri2 |

  - **Correction to the earlier feasibility text:** the `*_perm_purity_filtered` files are
    **SuSiE fine-mapping** results, not permutation leads.
  - **Correction to `UPDATED_DATA_INVENTORY.md`:** it says fine-mapped sQTL sets are "not public".
    They are public on Synapse.
  - **Trimester definitions** (local ST1-1 metadata, 654 samples):
    - Tri1: 4–13 pcw, n = 216.
    - Tri2: 13.4–26 pcw, n = 433.
    - Tri3: 28–39 pcw, n = 5.

    Both Tri1 and Tri2 are **prenatal**, so neither is an adult stage. The per-trimester sample
    counts after the EUR restriction have not been verified.
  - Trimester permutation leads are already local (ST3-5 Tri1: 14,193 introns / 5,312 sGenes; ST3-6
    Tri2: 5,348 / 2,318). Full-cohort group-permutation leads are local (ST1-4).
- **dbGaP phs000755** (raw reads plus genotypes, for E16): controlled access.

## Junction counts from MRF (`scripts/E15_mrf_junction_counts.py`)

**Rules:**
- A junction is a gap of ≥ 20 nt between consecutive alignment blocks of one mate.
- Coordinates are 1-based intron first and last bases.
- Counts are unstranded and taken once per read line.
- The CRC-32 must match the zip directory.

**Validation on HSB_103_A1C** (13 pcw), checked before any downstream use:

| Measure | Value |
|---|---|
| Reads or pairs | 29,098,744 |
| Spliced reads | 4,798,096 |
| Junctions | 151,190 |
| Runtime | ~75 s per sample |

Agreement with GENCODE v24lift37 annotated introns (hg19):

| Junction filter | Junctions | Annotated | Read-weighted annotated |
|---|---|---|---|
| ≥ 1 read | 151,190 | 90.4% | 99.1% |
| ≥ 5 reads | 92,202 | 98.3% | 99.5% |
| ≥ 20 reads | 46,614 | 99.6% | 99.7% |

The coordinate convention is therefore correct.

## Annotations downloaded 2026-11-05 (approved)

| File (`data/external/gencode/`) | Source | Size (B) | MD5 (matches GENCODE `MD5SUMS`) |
|---|---|---|---|
| `gencode.v19.annotation.gtf.gz` | ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/ | 37,991,892 | bd83e28270e595d3bde6bfcb21c9748f ✓ |
| `gencode.v10.annotation.gtf.gz` | …/release_10/ | 25,072,189 | 3f8404dbaf488370a8d5e06e9583cad5 ✓ |

- Official checksum lists are saved as `MD5SUMS.v19` and `MD5SUMS.v10`.
- Purpose:
  - v19: PsychENCODE DER-10 transcript IDs;
  - v10: the BrainSpan MRF annotation.
- **Check:** all 381 sampled PsychENCODE DER-10a transcript IDs, with version, exist in GENCODE v19
  (100%) [V]. v19 is therefore PEC's annotation.
- DER-10a and the SNP table are downloading to `data/external/psychencode_qtl/`. Sizes are verified
  on completion. **No slope will be read before the E15 pre-registration is frozen.**

## Completed downloads (2026-11-05)

No publisher checksum is available for any of these four files. Integrity is checked by exact byte
size, gzip test and per-row field counts. The MD5s below are ours, for future re-verification.

| File (`data/external/…`) | Size (B) | Checks | MD5 (ours) |
|---|---|---|---|
| `brainspan_pec/psi/brainSpan.hg38.gencode21.STAR.DEXSeq.exonPSI.txt` | 1,206,684,334 ✓ | 575,341 lines, all with 608 fields (1 ID + 607 samples) | 45bfa99db0416fc5ece5f1749b31c9e3 |
| `gtex_junctions/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz` | 4,245,392,642 ✓ | `gzip -t` OK | 69b59d8ec1aa4adeec354ada094cd1ea |
| `psychencode_qtl/DER-10a_hg19_isoQTL.significant.txt` | 390,983,712 ✓ | 2,628,260 lines, all with 15 fields | e6145312c9759c60ee91b980c7d86fba |
| `psychencode_qtl/SNP_Information_Table_with_Alleles.txt` | 217,980,289 ✓ | 5,312,509 lines, all with 6 fields | 0cb74c8e8a07eda64f9bfa5114b60017 |

The PSI and GTEx downloads were interrupted by connection resets and completed with byte-range resume.

## Wen Synapse files: access check (2026-11-05, Synapse REST `entity/{id}/bundle2`, anonymous)

The four E15 inputs are fixed by the plan (AUDIT_GAP_ANALYSIS_AND_PLAN.md G2; E15_PREFREEZE_VERIFICATION.md gap table).
All four are in folder syn51177091.

| ID (version) | File | Synapse description | Restriction | Anonymous download |
|---|---|---|---|---|
| syn51183885 (v3) | `mixed_cissqtl_40hcp_perm_purity_filtered.txt.gz` | SuSiE fine-mapping results of cis-sQTL | OPEN, no unmet requirement | not allowed (login required) |
| syn51183855 (v4) | `conditional_top_variants.txt` | Conditional cis-sQTL, full multi-ancestry dataset | OPEN | login required |
| syn51183854 (v4) | `leafcutter_clusters_to_genes.txt` | Intron annotation from GTEx script | OPEN | login required |
| syn51183859 (v4) | `all.introns.tested.tsv` | Intron annotation from customized leafviz | OPEN | login required |

Not E15 inputs (checked because they were proposed):
- syn51184373 (v2) `tri1_nominal_sqtl_15HCP_all_assoc.txt.gz`
- syn51184386 (v2) `tri2_nominal_sqtl_10HCP_all_assoc.txt.gz`

Both are OPEN. Their Synapse description is "Results from sQTL nominal mapping for Tri1/Tri2".

File sizes and publisher MD5s are not returned to anonymous requests. They will be recorded from
Synapse after an authenticated download.

### File sizes and type (2026-11-06, anonymous)

**Method:** Synapse Dataset syn50897018 (version 6, 68 items) was queried with
`SELECT id FROM syn50897018 WHERE id='<syn>'` and QueryBundleRequest partMask = 64 (sumFileSizes).

**Type:** all four inputs are ordinary `FileEntity` objects. They are not Synapse Tables, so their rows
cannot be queried or subset through the API.

The Dataset columns are file metadata only:
name, parentId, file_description, resourceType, id, description, modifiedOn, modifiedBy, type,
currentVersion, benefactorId, projectId.

The parent folder syn51177091 holds 18 files totalling 162,856,945,323 B.

| ID | File | Size (B) | Needed for E15 |
|---|---|---|---|
| syn51183885 | `mixed_cissqtl_40hcp_perm_purity_filtered.txt.gz` | 19,113,819 | yes |
| syn51183855 | `conditional_top_variants.txt` | 5,888,581 | yes |
| syn51183854 | `leafcutter_clusters_to_genes.txt` | 3,476,917 | yes |
| syn51183859 | `all.introns.tested.tsv` | 13,471,910 | yes |
| syn51184373 | `tri1_nominal_sqtl_15HCP_all_assoc.txt.gz` | 21,716,056,762 | no |
| syn51184386 | `tri2_nominal_sqtl_10HCP_all_assoc.txt.gz` | 21,804,119,064 | no |
| syn51184552 | `mixed_nominal_sqtl_40HCP_all_assoc.txt.gz` | 35,249,355,917 | no (optional; not requested) |

The four E15 inputs total 41,951,227 B (about 42 MB). Downloading them whole is cheaper than any subsetting.

**Still needed, and only available with an authenticated request:**
- the download URL;
- the file-handle MD5.

### Downloaded (2026-11-06, by the user from the Synapse website) and verified

Verification script: `scripts/E15_verify_wen_files.py`. It checks exact size against Synapse,
full-read integrity and line count, and records our MD5. Publisher MD5s were not supplied, so the MD5s
below are ours.

| ID | File (`data/external/wen_synapse/`) | Size (B) | Lines | MD5 (ours) | Columns |
|---|---|---|---|---|---|
| syn51183885 | `mixed_cissqtl_40hcp_perm_purity_filtered.txt.gz` | 19,113,819 ✓ | 560,689 (incl. header) | 9b5072076ae30fbc58638e8a30c34e82 | phenotype_id, variant_id, chr, pos, ref, alt, cs_id, cs_index, finemapped_region, pip, z, cs_min_r2, cs_avg_r2, cs_size, posterior_mean, posterior_sd, cs_log10bf |
| syn51183855 | `conditional_top_variants.txt` | 5,888,581 ✓ | 38,777 (no header) | 14ed16b51adaa4ffd979d075549c6d5e | 20 fields, unlabelled (QTLtools conditional-pass layout; to be confirmed from QTLtools documentation before use) |
| syn51183854 | `leafcutter_clusters_to_genes.txt` | 3,476,917 ✓ | 71,430 (incl. header) | 08527a522e62f6a39f943a3e53e18279 | clu, genes, names |
| syn51183859 | `all.introns.tested.tsv` | 13,471,910 ✓ | 273,167 (no header) | d13437680f3420ad20b57cf15a6f9100 | 2 fields: intron ID (`chr:start:end:clu_N_strand`), Ensembl gene ID |

**Pre-freeze disclosure.** The first version of the verification script printed the first line of each
file as "columns".
- `all.introns.tested.tsv` and `conditional_top_variants.txt` have no header row, so one data row of
  each was displayed.
- For `conditional_top_variants.txt` this was one intron–variant record, including its nominal
  p-value and slope.
- This single row was not used for any decision, and the script was fixed to print field counts only
  for files without a header.
- No other values have been read.
