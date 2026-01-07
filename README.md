# Variant Calling Pipeline using GATK

A comprehensive variant calling pipeline for genomic data analysis using GATK (Genome Analysis Toolkit) best practices. This pipeline processes raw FASTQ files through alignment, duplicate marking, variant calling, and quality filtering to produce high-quality variant calls suitable for downstream GWAS analysis.

## Pipeline Overview

This pipeline implements the GATK best practices workflow for germline variant discovery:

```
Raw FASTQ → Alignment → Mark Duplicates → Variant Calling → Joint Genotyping → Hard Filtering
```

## Pipeline Steps

### 1. Data Integrity Check (`1_interity`)
Validates the integrity of raw FASTQ files before processing.

**Features:**
- Checks gzip file integrity
- Identifies corrupted files early in the pipeline

**Usage:**
```bash
bash 1_interity
```

### 2. BWA-MEM Alignment (`2_alignment_script`)
Aligns paired-end reads to a reference genome using BWA-MEM.

**Features:**
- Automatic read group extraction from FASTQ headers
- Parallel processing with configurable threads (default: 32)
- Automatic SAM to BAM conversion and sorting
- BAM indexing
- Comprehensive alignment statistics

**Configuration:**
- Reference: `references/genome.fa`
- Input: `raw_data/*_R1.fastq.gz` and `*_R2.fastq.gz`
- Output: `alignment/*.bam`
- Threads: 32

**Usage:**
```bash
bash 2_alignment_script
```

**Output:**
- Sorted and indexed BAM files
- Alignment summary with mapping rates
- Log files for each sample

### 3. Mark Duplicates (`3_mark_dup`)
Identifies and marks PCR/optical duplicates using GATK MarkDuplicates.

**Features:**
- Automatic BAM validation
- Index creation
- Duplicate metrics generation
- Resume capability (skips already processed files)

**Configuration:**
- Input: `alignment/*.bam`
- Output: `alignment/marked/*_marked.bam`
- Threads: 32

**Usage:**
```bash
bash 3_mark_dup
```

**Output:**
- Marked BAM files with duplicates flagged
- Duplicate metrics for each sample
- BAM indices

### 4. HaplotypeCaller (`4_haplotype_caller`)
Calls variants per sample in GVCF mode for joint genotyping.

**Features:**
- GVCF output for efficient joint genotyping
- Enhanced GQ bands (1-60 individually, then 70, 80, 90, 99)
- Configurable ploidy (default: 2)
- Quality filtering (min confidence: 30.0, min mapping quality: 20)

**Configuration:**
- Reference: `references/genome.fa`
- Input: `alignment/marked/*.bam`
- Output: `gatk/haplotypecaller/*.g.vcf.gz`
- Ploidy: 2
- Threads: 32

**Usage:**
```bash
bash 4_haplotype_caller
```

**Output:**
- GVCF files for each sample
- Variant statistics
- Processing logs

### 5. GenomicsDB Import (`5_GenomicDB`)
Consolidates individual GVCFs into a GenomicsDB workspace for efficient joint genotyping.

**Features:**
- Automatic sample map creation
- Interval-based processing
- Parallel interval import (configurable)
- Database validation

**Configuration:**
- Input: `gatk/haplotypecaller/*.g.vcf.gz`
- Output: `gatk/genomicsdb/genomicsdb_database`
- Parallel intervals: 4

**Usage:**
```bash
bash 5_GenomicDB
```

**Output:**
- GenomicsDB workspace
- Sample name map
- Genome intervals list
- Import summary

### 6. Joint Genotyping (`6_joint_genotype`)
Performs joint genotyping across all samples to generate a multi-sample VCF.

**Features:**
- Joint genotyping from GenomicsDB
- Configurable heterozygosity rates
- Ti/Tv ratio calculation for quality assessment
- Comprehensive variant statistics

**Configuration:**
- Input: `gatk/genomicsdb/genomicsdb_database`
- Output: `gatk/joint_genotyping/cohort_raw.vcf.gz`
- Ploidy: 2
- Min confidence: 30.0
- Heterozygosity: 0.001
- Indel heterozygosity: 0.000125

**Usage:**
```bash
bash 6_joint_genotype
```

**Output:**
- Raw multi-sample VCF
- Variant statistics (SNPs, indels, Ti/Tv ratio)
- Quality assessment report

### 7. Hard Filtering (`7_hard_filter`)
Applies GATK hard filters to remove low-quality variants.

**Features:**
- Multiple quality filters based on GATK best practices
- Separate PASS-only VCF generation
- Filter breakdown statistics
- Ti/Tv ratio on filtered variants

**Filter Criteria:**
- `QD < 2.0`: Quality by Depth
- `FS > 60.0`: Fisher Strand bias
- `MQ < 40.0`: Mapping Quality
- `MQRankSum < -12.5`: Mapping Quality Rank Sum
- `ReadPosRankSum < -8.0`: Read Position Rank Sum
- `SOR > 3.0`: Strand Odds Ratio
- Clustering: 3 variants in 0 bp window

**Usage:**
```bash
bash 7_hard_filter
```

**Output:**
- Filtered VCF with filter tags: `gatk/filtered/cohort_hard_filtered.vcf.gz`
- PASS-only VCF: `gatk/filtered/cohort_PASS_only.vcf.gz`
- Filtration summary report
- Quality metrics

## Prerequisites

### Software Requirements
- **BWA** (>= 0.7.17): For read alignment
- **SAMtools** (>= 1.10): For BAM manipulation
- **GATK** (>= 4.0): For variant calling and filtering
- **BCFtools** (>= 1.10): For VCF statistics and filtering

### Reference Genome
- Reference FASTA: `references/genome.fa`
- BWA index: `references/genome.fa.bwt` (and other index files)
- GATK dictionary: `references/genome.dict`
- FAI index: `references/genome.fa.fai`

**Create required indices:**
```bash
# BWA index
bwa index references/genome.fa

# SAMtools index
samtools faidx references/genome.fa

# GATK dictionary
gatk CreateSequenceDictionary -R references/genome.fa
```

## Directory Structure

```
variant_calling/
├── raw_data/                    # Input FASTQ files
│   ├── sample1_R1.fastq.gz
│   └── sample1_R2.fastq.gz
├── references/                  # Reference genome and indices
│   ├── genome.fa
│   ├── genome.fa.fai
│   ├── genome.dict
│   └── genome.fa.bwt (and other BWA indices)
├── alignment/                   # Aligned BAM files
│   └── marked/                  # Duplicate-marked BAM files
├── gatk/                        # GATK outputs
│   ├── haplotypecaller/        # Per-sample GVCFs
│   ├── genomicsdb/             # GenomicsDB workspace
│   ├── joint_genotyping/       # Multi-sample VCF
│   ├── filtered/               # Filtered VCFs
│   ├── logs/                   # Processing logs
│   └── tmp/                    # Temporary files
└── logs/                        # Pipeline logs
```

## Running the Complete Pipeline

Execute scripts in order:

```bash
# 1. Check data integrity
bash 1_interity

# 2. Align reads
bash 2_alignment_script

# 3. Mark duplicates
bash 3_mark_dup

# 4. Call variants (per sample)
bash 4_haplotype_caller

# 5. Import to GenomicsDB
bash 5_GenomicDB

# 6. Joint genotyping
bash 6_joint_genotype

# 7. Apply hard filters
bash 7_hard_filter
```

## Quality Metrics

### Ti/Tv Ratio (Transition/Transversion)
A key quality indicator for SNP calls:
- **≥ 2.0**: Excellent quality
- **1.8-2.0**: Good quality
- **1.5-1.8**: Fair quality (consider additional filtering)
- **< 1.5**: Poor quality (strong filtering recommended)

### Mapping Rate
Expected mapping rates for good quality data:
- **> 95%**: Excellent
- **90-95%**: Good
- **< 90%**: May indicate contamination or reference mismatch

## Downstream Analysis

### Additional Filtering for GWAS

Apply MAF and missingness filters:
```bash
bcftools view -i 'INFO/AF>=0.05 & F_MISSING<=0.1' \
    gatk/filtered/cohort_PASS_only.vcf.gz \
    -Oz -o cohort_gwas_ready.vcf.gz
```

### Convert to PLINK Format

```bash
plink --vcf gatk/filtered/cohort_PASS_only.vcf.gz \
      --make-bed --out cohort_plink \
      --allow-extra-chr
```

## Performance Considerations

- **Threads**: Adjust `THREADS` variable in each script based on available CPU cores
- **Memory**: GATK steps may require 16-32 GB RAM depending on genome size
- **Storage**: Ensure sufficient disk space (BAM files can be large)
- **Runtime**: 
  - Alignment: 15-30 minutes per sample
  - HaplotypeCaller: 30-60 minutes per sample
  - GenomicsDB Import: 1-3 hours (depends on sample count)
  - Joint Genotyping: 2-4 hours (for ~300 samples)

## Troubleshooting

### Common Issues

1. **Missing BAM index**: Scripts automatically create indices, but if errors occur:
   ```bash
   samtools index alignment/sample.bam
   ```

2. **GATK memory errors**: Increase Java heap size:
   ```bash
   export JAVA_OPTS="-Xmx32g"
   ```

3. **Reference dictionary not found**: Create using:
   ```bash
   gatk CreateSequenceDictionary -R references/genome.fa
   ```

## Citation

If you use this pipeline, please cite:

- **GATK**: Van der Auwera GA & O'Connor BD. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition). O'Reilly Media.
- **BWA**: Li H. and Durbin R. (2009) Fast and accurate short read alignment with Burrows-Wheeler Transform. Bioinformatics, 25:1754-60.
- **SAMtools**: Li H., Handsaker B., Wysoker A., et al. (2009) The Sequence Alignment/Map format and SAMtools. Bioinformatics, 25:2078-9.

## License

This pipeline is provided as-is for academic and research purposes.

## Author

Sudhakar Basu

## Repository

https://github.com/SudhakarBasu/Variant_calling_GATK
