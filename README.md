# Dock Analysis

> **Universal docking result analyzer** — HDOCK, HADDOCK 2.4, Vina, ClusPro output to reports, visualizations, and binding analysis.

[![Version](https://img.shields.io/badge/version-2.0.0-blue.svg)](https://github.com/CliffVale/dock_analysis)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

---

## What's New in v2.0

- **HADDOCK 2.4 support** — parse HADDOCK output (score lists, cluster analysis)
- **Multi-server comparison** — analyze HDOCK and HADDOCK results side-by-side
- **Enhanced cluster analysis** — automatic clustering by score proximity
- **H-bond detection** — identify hydrogen bonds at protein-DNA interface

---

## Quick Start

```bash
# Analyze HDOCK results
python scripts/hdock_analyzer.py --input hdock_results.tar.gz --output analysis/

# Analyze HADDOCK results
python scripts/haddock_analyzer.py --input run1/ --output analysis/

# Compare multiple HADDOCK runs
python scripts/haddock_analyzer.py --input run1/ run2/ run3/ --output comparison/
```

---

## Analyzer Scripts

| Script | Input | Server |
|--------|-------|--------|
| `hdock_analyzer.py` | HDOCK tar.gz or extracted folder | HDOCK |
| `haddock_analyzer.py` | HADDOCK output directory | HADDOCK 2.4 |
| `compare_docks.py` | Multiple docking results | Any |

---

## What Each Script Does

### `hdock_analyzer.py` — HDOCK Analysis

- Parses HDOCK `.out` score files
- Extracts interface contact files
- Generates binding residue frequency analysis
- Produces PyMOL visualization scripts
- Outputs Markdown + HTML reports

### `haddock_analyzer.py` — HADDOCK Analysis

- Parses `the_score.list` files
- Analyzes individual model PDBs for contacts
- Automatic cluster detection (by score proximity)
- Decomposes HADDOCK score into components (vdW, elec, desolv, BSA)
- Generates PyMOL scripts with cluster visualization
- Outputs Markdown + HTML reports + JSON

### `compare_docks.py` — Multi-Result Comparison

- Compare scores across different docking runs
- Identify common binding residues
- Generate comparative reports

---

## Output Files

### HDOCK Analysis Output
```
analysis/
├── report.md                    # Markdown report
├── report.html                  # HTML report (open in browser)
├── hdock_analysis.json          # Machine-readable results
├── pymol_visualize.py           # PyMOL visualization script
└── binding_residues.csv         # Contact residue details
```

### HADDOCK Analysis Output
```
analysis/
├── report.md                    # Markdown report
├── report.html                  # HTML report
├── haddock_analysis.json        # Machine-readable results
├── pymol_visualize.py           # PyMOL visualization script
└── binding_residues.csv         # Contact residue details
```

---

## HADDOCK Score Components

HADDOCK decomposes the total score into physically meaningful terms:

| Component | Description | Weight |
|-----------|-------------|--------|
| vdw | van der Waals interactions | 0.1 |
| elec | Electrostatic interactions | 1.0 |
| desolv | Desolvation energy | 1.0 |
| AIR | Ambiguous Interaction Restraints | 0.01 |
| BSA | Buried Surface Area | -0.005 |

**Total HADDOCK score** = 0.1×vdW + 1.0×elec + 1.0×desolv + 0.01×AIR − 0.005×BSA

---

## Usage Examples

### Analyze HDOCK Tarball

```bash
python scripts/hdock_analyzer.py \
    --input results/hdock_job12345.tar.gz \
    --output analysis/hdock_job12345
```

### Analyze HADDOCK Run

```bash
python scripts/haddock_analyzer.py \
    --input /path/to/haddock_run1/ \
    --output analysis/haddock_run1
```

### Compare Multiple HADDOCK Runs

```bash
python scripts/haddock_analyzer.py \
    --input run1/ run2/ run3/ \
    --output comparison/
```

### Generate Only HTML Report

```bash
python scripts/haddock_analyzer.py \
    --input run1/ \
    --output analysis/ \
    --format html
```

---

## PyMOL Visualization

After analysis, run the generated PyMOL script:

```bash
# View in PyMOL
pymol analysis/pymol_visualize.py

# Or from PyMOL GUI
File → Run Script → analysis/pymol_visualize.py
```

The script shows:
- Protein receptor (gray cartoon + surface)
- DNA/RNA ligand (orange sticks + phosphorus spheres)
- Binding site residues (yellow sticks)
- Hydrogen bonds (dashed lines)
- Top cluster models (colored by cluster)

---

## Installation

### Prerequisites

- Python 3.9+
- Conda/Mamba

### Setup

```bash
# Clone repository
git clone https://github.com/CliffVale/dock_analysis.git
cd dock_analysis

# Create environment
conda env create -f environment.yml
conda activate dock_analysis

# Verify
python scripts/hdock_analyzer.py --help
python scripts/haddock_analyzer.py --help
```

---

## Project Structure

```
dock_analysis/
├── README.md                    # This file
├── environment.yml              # Conda environment
├── scripts/
│   ├── hdock_analyzer.py        # HDOCK output parser
│   ├── haddock_analyzer.py      # HADDOCK output parser
│   └── compare_docks.py         # Multi-result comparison
├── data/                        # Sample data
├── results/                     # Example outputs
├── examples/                    # Usage examples
└── docs/
    └── TUTORIAL.md              # Step-by-step guide
```

---

## Citation

If you use this tool, please cite:

1. **HDOCK**: Yan et al. (2017) "HDOCK: a web server for protein-protein and protein-DNA/RNA docking" *Nucleic Acids Research*
2. **HADDOCK**: van Zundert et al. (2016) "The HADDOCK2.4 web server" *Nucleic Acids Research*

---

## License

MIT License

---

## Contributing

Contributions welcome! Please open an issue or PR.
