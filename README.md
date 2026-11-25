# GeneWeaver – GPU-Accelerated CRISPR Alignment Engine

## Project Overview
GeneWeaver is a Python-based bioinformatics project designed to identify potential CRISPR off-target DNA sequence matches using sequence alignment techniques.

The project aims to explore high-performance computing through CPU-based alignment, GPU acceleration, and distributed computing.

## Current Development Status

**Stage 1 – CPU-Based DNA Alignment (Completed)**

Implemented features:
- FASTA file parsing using BioPython.
- DNA sequence validation.
- Sliding-window sequence alignment.
- Mismatch detection.
- Identification of exact and approximate DNA matches.
- Execution-time measurement.
- Initial chunk-based alignment implementation.

## Technologies Used
- Python
- BioPython
- VS Code
- Git and GitHub

## Project Structure
```text
GeneWeaver/
├── data/
│   └── sample.fasta
├── src/
│   └── main.py
├── results/
├── requirements.txt
├── .gitignore
└── README.md
```

## Installation

Clone the repository:

```bash
git clone YOUR_GITHUB_REPOSITORY_URL
cd GeneWeaver
```

Create a virtual environment:

```bash
python -m venv venv
```

Activate it on Windows:

```bash
venv\Scripts\activate
```

Install dependencies:

```bash
pip install -r requirements.txt
```

## Run the Project

```bash
python src/main.py
```

Enter a DNA target sequence containing A, T, C and G when prompted.

## Current Testing

Test dataset: Artificial FASTA sample.

Target DNA: `ATCGGCTA`

Maximum allowed mismatches: 2

Initial test results:
- Total potential matches: 6
- Exact matches: 3
- One-mismatch matches: 3

## Future Development

- Stage 2: Numba CUDA GPU acceleration.
- Stage 3: Dask distributed computing.
- Stage 4: Biological scoring and PAM recognition.
- Stage 5: Performance optimization and Textual terminal dashboard.

## Disclaimer
This is an educational and experimental bioinformatics project. Its current sequence-matching results are not validated predictions of CRISPR editing outcomes.

**Project Status:** Stage 1 – CPU Alignment Completed.# -GeneWeaver-GPU-Accelerated-CRISPR-Alignment-Engine
