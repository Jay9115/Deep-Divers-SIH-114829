# Deep-Sea eDNA Analysis Pipeline

**AI-Driven Environmental DNA Analysis for Deep-Sea Biodiversity Assessment**

[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://python.org)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![SIH 2025](https://img.shields.io/badge/SIH-2025-orange.svg)](https://www.sih.gov.in)

---

## 🌊 Project Overview

This project addresses the critical challenge of analyzing environmental DNA (eDNA) from deep-sea ecosystems, where traditional bioinformatics approaches fail due to poor database representation of novel taxa. Our AI-driven pipeline provides a comprehensive solution for processing raw Illumina sequencing data from deep-sea samples, enabling rapid taxonomic classification and biodiversity assessment.

**Problem Statement ID:** SIH25042  
**Team ID:** 002200

### 🎯 Key Innovation

Unlike traditional pipelines that rely heavily on reference databases (QIIME2, mothur), our approach combines:
- **3-Tier Classification System**: Exact matching → ML embeddings → Novel characterization
- **Marker-Aware Processing**: Optimized workflows for 18S rRNA and COI genes
- **Deep Learning Integration**: Fine-tuned DNABERT-2/FCCGR-CNN models for taxonomic inference

---

## 📊 System Architecture

```
INPUT LAYER → QC & ASV → FEATURE EXTRACTION → CLASSIFICATION → ABUNDANCE → ECOLOGY → OUTPUT
     ↓             ↓              ↓               ↓            ↓          ↓         ↓
  FASTQ Files   Module 1       Module 2       Module 3     Module 4   Module 5  Reports
  Metadata                                                                        & Viz
  Controls                                    (GPU)
```

### Processing Modules

1. **Module 1: Marker-Aware QC & ASV Generation**
   - Automatic primer detection and marker identification
   - Parallel 18S (protist-optimized) and COI (metazoan-optimized) pipelines
   - Quality control, trimming, merging, and ASV generation using DADA2
   - 3-stage chimera removal and contamination filtering

2. **Module 2: Dual Feature Extraction**
   - K-mer spectrum analysis with MinHash LSH (fast path)
   - FCGR (Frequency Chaos Game Representation) generation for AI models
   - Parallel processing for computational efficiency

3. **Module 3: Hybrid 3-Tier Classification**
   - **Tier 1**: Exact matching for 80-85% of ASVs using MinHash LSH
   - **Tier 2**: ML embeddings via fine-tuned DNABERT-2 for 10-12% of ASVs
   - **Tier 3**: Novel characterization with phylogenetic placement for 5-8% of ASVs

4. **Module 4: Abundance Estimation & Correction**
   - rRNA copy number correction and PCR bias modeling
   - Optional spike-in calibration for absolute abundance
   - Rare taxa validation and quality control

5. **Module 5: Ecological Analysis**
   - Alpha and beta diversity calculations
   - Taxonomic composition and co-occurrence networks
   - Deep-sea specific validation and biogeographic checks

---

## 🔬 Technical Specifications

### AI Models & Algorithms
- **DNABERT-2**: Fine-tuned on 85K 18S + 120K COI sequences
- **FAISS HNSW**: High-dimensional similarity search
- **MinHash LSH**: Fast approximate matching
- **EPA-ng**: Phylogenetic placement for novel taxa

### Reference Databases (3.5+ GB)
- **PR2 v5.0.0**: 235K 18S protist sequences
- **SILVA 138.1**: 510K eukaryote SSU sequences
- **MIDORI2 v253**: 1.8M COI metazoan sequences
- **DeepSeaDB v1.0**: 12K curated deep-sea sequences

### System Requirements

**Minimum:**
- 16-core CPU, 32GB RAM
- 100GB SSD storage
- Python 3.10+, R 4.3+

**Recommended:**
- 32-core CPU, 64GB RAM
- RTX 4090 GPU (for 3x speed improvement)
- 200GB SSD storage

---

## 🚀 Installation & Setup

### 1. Clone Repository
```bash
git clone https://github.com/Jay9115/Deep-Divers-SIH-114829.git
cd Deep-Divers-SIH-114829
```

### 2. Environment Setup
```bash
# Using Conda (recommended)
conda create -n deepsea-edna python=3.10
conda activate deepsea-edna
pip install -r requirements.txt

# Or using Docker
docker build -t deepsea-edna .
docker run -it deepsea-edna
```

### 3. Install Dependencies
```bash
# R packages for DADA2
Rscript -e "install.packages(c('dada2', 'phyloseq', 'BiocManager'))"

# Download reference databases
python setup.py download_databases
```

---

## 💻 Usage

### Command Line Interface

```bash
# Quick start with demo data
deepsea-edna --input demo_fastq/ --output results/ --marker auto

# Full pipeline with custom parameters
deepsea-edna \
  --input /path/to/fastq/ \
  --metadata samples.tsv \
  --output results/ \
  --marker 18S \
  --threads 16 \
  --gpu \
  --confidence 0.75
```

### Python API

```python
from deepsea_edna_pipeline import DeepSeaAnalyzer

# Initialize pipeline
analyzer = DeepSeaAnalyzer(
    input_dir="fastq_files/",
    output_dir="results/",
    use_gpu=True
)

# Run analysis
results = analyzer.run_full_pipeline()



### Web Interface

```bash
# Start web server
deepsea-edna serve --port 8080

# Access at http://localhost:8080
```

---

## 📁 Project Structure

```
deepsea_edna_pipeline/
├── cli/                    # Command-line interface
│   ├── asv_cli.py         # ASV generation commands
│   └── features_cli.py    # Feature extraction commands
├── module1_qc_asv/        # Quality control & ASV generation
│   ├── pipeline_asv.py    # Main ASV pipeline
│   ├── primers/           # Primer sequences
│   └── scripts/           # Processing scripts
├── module2_features/      # Feature extraction
│   ├── fcgr/             # FCGR generation
│   └── kmer/             # K-mer analysis
├── utils/                 # Utility functions
├── data-DB/              # Demo data
├── tests/                # Unit tests
└── UI/                   # Web interface
```

---


### Primary Outputs
- `taxonomy.tsv` - Complete taxonomic assignments with confidence scores
- `asv_table.biom` - ASV abundance table in BIOM format
- `sequences.fasta` - Representative ASV sequences
- `novel_taxa.md` - Novel taxa characterization report

### Visualization
- `interactive.html` - Interactive dashboard with UMAP plots
- `figures/` - Publication-ready plots (300 DPI PNG/PDF)
- `diversity_plots.pdf` - Alpha/beta diversity visualizations

### Integration
- `phyloseq.rds` - R phyloseq object for further analysis
- `processing.log` - Complete processing log with parameters

---

## 🎥 Demo Video

[![Deep-Sea eDNA Pipeline Demo](https://img.youtube.com/vi/YOUR_VIDEO_ID/0.jpg)](https://www.youtube.com/watch?v=YOUR_VIDEO_ID)

*Click above to watch our pipeline demonstration*

---

## 👥 Team Information

**Team Name:** Deep Divers  
**Team ID:** 114829
**Problem Statement:** SIH25042

### Team Members

| Role | Name | College ID |
|------|------|-----------|
| **Team Leader** | Mahi Patel | 23DCS076 |
| Member | Jay Patel | 23DCS076 |
| Member | Samarth Patel | 23DCS089 |
| Member | Shubh Patel | 23DCS092 |
| Member | Jaimin Raval | 23DCS111 | 
| Member | Ashok Suthar | 23DCS130 | 

---

## 🔬 Scientific Background

### The Deep-Sea Challenge

The deep ocean harbors immense biodiversity, much of it unknown to science. Traditional biodiversity assessment methods are limited by:

1. **Accessibility**: Deep-sea environments are remote and difficult to sample
2. **Database Gaps**: Reference databases lack deep-sea organism sequences
3. **Computational Limitations**: Existing tools are slow and database-dependent
4. **Novel Taxa**: High proportion of unknown species in deep-sea samples

### Our Solution

Environmental DNA (eDNA) analysis provides a non-invasive method to assess biodiversity by detecting genetic material in water and sediment samples. Our AI-driven approach overcomes traditional limitations by:

- Minimizing reliance on incomplete reference databases
- Using machine learning to identify novel taxa
- Providing rapid, accurate taxonomic assignments
- Enabling discovery of previously unknown species

---

## 📚 Dependencies

### Core Libraries
```
Python >= 3.10
torch >= 2.0.0
transformers >= 4.35.0
faiss-cpu >= 1.7.4
numpy >= 1.24.0
pandas >= 2.0.0
biopython >= 1.81
scikit-learn >= 1.3.0
```

### Bioinformatics Tools
```
cutadapt >= 4.6
fastp >= 0.23.0
vsearch >= 2.22.0
dada2 >= 1.28.0 (R)
EPA-ng >= 0.3.8
```

---

## 🤝 Contributing

We welcome contributions from the scientific community! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details on:

- Code standards and testing
- Documentation requirements
- Submission process
- Bug reporting

---

## 🙏 Acknowledgments

- Centre for Marine Living Resources and Ecology (CMLRE)
- Smart India Hackathon 2025
- Hugging Face for DNABERT-2 pre-trained models
- PR2, SILVA, and MIDORI database maintainers

---

*Advancing deep-sea biodiversity research through AI-driven eDNA analysis*