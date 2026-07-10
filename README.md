# Docking Result Analyzer

> **Universal tool for analyzing molecular docking outputs** — HDOCK, AutoDock Vina, ClusPro, and more.

> 👤 **Author**: CliffVale (Bhrigu)
> 
> 📅 **Version**: 1.0.0 (July 2026)

---

## What Does This Do?

```
┌─────────────────────────────────────────────────────────────────────┐
│  📦 INPUT: Docking results (tar.gz or folder)                       │
│                                                                     │
│  🔧 PROCESSING:                                                     │
│    • Parse scores and rankings                                      │
│    • Analyze binding interfaces                                     │
│    • Identify key contact residues                                  │
│    • Generate visualization scripts                                 │
│    • Create reports (HTML/Markdown)                                 │
│                                                                     │
│  📊 OUTPUT:                                                         │
│    • Ranked models with scores                                      │
│    • Binding residue analysis                                       │
│    • PyMOL visualization scripts                                    │
│    • Interactive HTML report                                        │
│    • Machine-readable JSON data                                     │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Quick Start (For Beginners)

### Step 1: Install

```bash
# Clone the repository
git clone https://github.com/CliffVale/dock_analysis.git
cd dock_analysis

# Create environment (optional, works with Python 3.9+)
conda env create -f environment.yml
conda activate dock_analysis

# Or just install dependencies
pip install numpy matplotlib
```

### Step 2: Get Your Docking Results

Download results from your docking server:

| Server | Download Format |
|--------|-----------------|
| **HDOCK** | Click "All results" → Downloads `all_results_*.tar.gz` |
| **ClusPro** | Download job zip file |
| **AutoDock Vina** | Save `.pdbqt` output file |

### Step 3: Analyze

```bash
# Analyze HDOCK results (tar.gz)
python scripts/hdock_analyzer.py --input all_results.tar.gz --output my_analysis/

# Analyze extracted folder
python scripts/hdock_analyzer.py --input hdock_output_folder/ --output my_analysis/
```

### Step 4: View Results

```bash
# Open HTML report in browser
xdg-open my_analysis/report.html  # Linux
open my_analysis/report.html      # macOS

# Open PyMOL visualization
pymol my_analysis/pymol_visualize.py
```

---

## What You Get

### 📄 Reports

| File | Description | How to Use |
|------|-------------|------------|
| `report.html` | **Interactive HTML report** | Open in any web browser |
| `report.md` | Markdown report | View on GitHub or any markdown viewer |
| `EXECUTIVE_SUMMARY.md` | Quick summary | For presentations/papers |

### 📊 Data Files

| File | Description | How to Use |
|------|-------------|------------|
| `hdock_analysis.json` | Machine-readable results | Import in Python, R, or any tool |
| `binding_residues.csv` | Contact residue details | Open in Excel/Google Sheets |
| `comparison.json` | Multi-dock comparison | For comparing aptamers |

### 🎨 Visualization

| File | Description | How to Use |
|------|-------------|------------|
| `pymol_visualize.py` | PyMOL script | `pymol pymol_visualize.py` |
| `score_chart.md` | ASCII bar chart | View in any text editor |
| `*.pse` | PyMOL session | Double-click to open |

---

## Understanding the Output

### HDOCK Scores

```
┌─────────────────────────────────────────────────────────────┐
│  HDOCK SCORE INTERPRETATION                                 │
│                                                             │
│  Score Range        │  Binding Strength   │  Confidence     │
│  ─────────────────  │  ─────────────────  │  ────────────   │
│  -300 to -200       │  Excellent          │  High           │
│  -200 to -150       │  Good               │  Moderate       │
│  -150 to -100       │  Moderate           │  Low            │
│  > -100             │  Weak/None          │  Very Low       │
│                                                             │
│  ⚠️  Scores are relative, not absolute binding affinities   │
│  ⚠️  Always validate experimentally (SPR, ITC, BLI)         │
└─────────────────────────────────────────────────────────────┘
```

### Binding Interface

```
┌─────────────────────────────────────────────────────────────┐
│  BINDING RESIDUE ANALYSIS                                   │
│                                                             │
│  What are "contact residues"?                               │
│  → Amino acids on the protein surface that touch the DNA    │
│                                                             │
│  How are they identified?                                   │
│  → Count how often each residue appears in top 10 models    │
│  → High frequency = likely real binding site                │
│                                                             │
│  Example output:                                            │
│    TRP67: 10/10 models  ← Key binding residue               │
│    PHE39: 8/10 models   ← Important for stacking            │
│    TYR73: 7/10 models   ← Contributing residue              │
│                                                             │
│  ⚠️  Residues in >5/10 models are reliable binding sites    │
└─────────────────────────────────────────────────────────────┘
```

### Visualization Colors

```
┌─────────────────────────────────────────────────────────────┐
│  PYMOL COLOR SCHEME                                         │
│                                                             │
│  🔵 Gray/White  →  Protein (receptor)                       │
│  🟠 Orange      →  DNA (ligand)                             │
│  🟡 Yellow      →  Binding site residues                    │
│  🟢 Green       →  Hydrogen bonds                           │
│  🔴 Red spheres →  Phosphate backbone (DNA)                 │
│                                                             │
│  To customize: Edit pymol_visualize.py                      │
└─────────────────────────────────────────────────────────────┘
```

---

## Advanced Usage

### Compare Multiple Docking Runs

```bash
# Compare different aptamers
python scripts/compare_docks.py \
    --results results_20mer/ results_50mer/ results_trail_seq10/ \
    --output comparison/
```

### Batch Processing

```bash
# Analyze multiple tar.gz files
for file in *.tar.gz; do
    name="${file%.tar.gz}"
    python scripts/hdock_analyzer.py --input "$file" --output "analysis/$name"
done

# Compare all
python scripts/compare_docks.py --results analysis/*/ --output comparison/
```

### Custom Analysis

```python
from scripts.hdock_analyzer import HDOCKAnalyzer

# Load and analyze
analyzer = HDOCKAnalyzer("results.tar.gz", "output/")
results = analyzer.parse_hdock_output()

# Access data programmatically
for model in results.models[:5]:
    print(f"Model {model.model_id}: {model.score}")
```

---

## Troubleshooting

### "No .out file found"
- **Cause**: Wrong input format
- **Fix**: Make sure you're pointing to the HDOCK output folder or tar.gz

### "No hdock_analysis.json"
- **Cause**: Haven't run analyzer yet
- **Fix**: Run `hdock_analyzer.py` first

### PyMOL script fails
- **Cause**: PyMOL not installed or not in PATH
- **Fix**: Install PyMOL: `conda install -c conda-forge pymol`

### Bad scores (all similar)
- **Cause**: Docking failed or poor quality input
- **Fix**: Check input PDB quality, re-run docking

---

## Example Walkthrough

### Input: HDOCK Results for CRP Aptamer

```
all_results_mCRP-vs-trail_seq10.tar.gz
├── hdock_6a50e897df44c.out      ← Scores
├── interface_1.txt              ← Contacts (model 1)
├── interface_2.txt              ← Contacts (model 2)
├── ...
├── model_1.pdb                  ← Best docking pose
├── model_2.pdb
├── ...
├── rec_*.pdb                    ← Receptor
└── lig_*.pdb                    ← Ligand
```

### Command

```bash
python scripts/hdock_analyzer.py \
    --input all_results_mCRP-vs-trail_seq10.tar.gz \
    --output analysis_trail_seq10/
```

### Output

```
analysis_trail_seq10/
├── report.html                  ← Open in browser
├── report.md                    ← Markdown version
├── hdock_analysis.json          ← For further analysis
├── pymol_visualize.py           ← Run in PyMOL
├── binding_residues.csv         ← Open in Excel
└── score_chart.md               ← ASCII visualization
```

---

## Technical Details

### Parsing Algorithm

1. **Extract** tar.gz if compressed
2. **Parse** `.out` file for scores and metadata
3. **Parse** `interface_*.txt` for contact residues
4. **Aggregate** contacts across top 10 models
5. **Rank** residues by frequency
6. **Generate** scripts and reports

### Scoring Metrics

| Metric | Source | Meaning |
|--------|--------|---------|
| `score` | Column 7 in .out | HDOCK docking score (kcal/mol) |
| `confidence` | Column 2 | Confidence score (0-1) |
| `rmsd` | Column 3 | Ligand RMSD from reference |

### Contact Detection

- Distance cutoff: 4.0 Å (van der Waals contact)
- H-bond cutoff: 3.5 Å (hydrogen bond)
- Stacking: Parallel/displaced < 4.5 Å

---

## FAQ

**Q: Are HDOCK scores binding affinities?**
A: No. HDOCK scores are relative rankings. Actual Kd requires experimental measurement.

**Q: How many models should I analyze?**
A: Top 10-20. Beyond rank 20, scores plateau and poses become unreliable.

**Q: Can I use this for Vina results?**
A: Currently optimized for HDOCK. Vina support coming in v2.0.

**Q: Why do some residues appear in all models?**
A: They're likely genuine binding site residues. >5/10 models = high confidence.

---

## Citation

If you use this tool, please cite:

```
CliffVale (2026). Docking Result Analyzer: Universal tool for molecular 
docking analysis. GitHub repository: https://github.com/CliffVale/dock_analysis
```

---

## License

MIT License

---

## Credits

**Author**: CliffVale (Bhrigu)

**Built with**: Python, NumPy, Matplotlib

**Tested on**: HDOCK, AutoDock Vina outputs

**Part of**: CRP Aptamer Wearable Biosensor Project, IIT Delhi

---

## Contributing

Contributions welcome! Please open an issue or PR.

### Ideas for Improvement
- [ ] Support for ClusPro results
- [ ] Support for HADDOCK results
- [ ] Binding free energy calculations
- [ ] MD trajectory analysis
- [ ] Interactive Jupyter notebook
- [ ] Web interface
