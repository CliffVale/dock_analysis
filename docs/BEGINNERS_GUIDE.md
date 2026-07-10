# Beginner's Guide to Docking Analysis

> **No experience needed!** This guide walks you through analyzing your first docking result.

---

## What is Molecular Docking?

Molecular docking predicts how two molecules (like a protein and DNA) bind together. It's like a molecular puzzle:

```
┌─────────────────────────────────────────────────────────────┐
│  🧩 THE DOCKING PUZZLE                                      │
│                                                             │
│  Protein (Receptor)  +  DNA (Ligand)  →  Bound Complex     │
│                                                             │
│  ┌─────────┐          ┌─────┐         ┌─────────┐          │
│  │         │    +     │ ═══ │    =    │    ═══  │          │
│  │  (●)    │          └─────┘         │   (●)═══│          │
│  │         │                          │         │          │
│  └─────────┘                          └─────────┘          │
│                                                             │
│  🔵 Protein    🟠 DNA    🟡 Binding site                    │
└─────────────────────────────────────────────────────────────┘
```

---

## Step-by-Step Tutorial

### Step 1: Get Your Docking Results

After docking (e.g., on HDOCK server), download the results:

```
📦 What you download:
├── all_results.tar.gz     ← This is your results file
└── OR
└── results_folder/        ← Already extracted folder
```

**How to download from HDOCK:**
1. Go to your HDOCK results page
2. Click "All results in a package"
3. Save the `.tar.gz` file

---

### Step 2: Install This Tool

Open terminal and run:

```bash
# Clone this tool
git clone https://github.com/CliffVale/dock_analysis.git
cd dock_analysis

# Install Python (if you don't have it)
# On Linux:
sudo pacman -S python    # Arch/CachyOS
sudo apt install python3 # Ubuntu/Debian

# Install required packages
pip install numpy matplotlib
```

---

### Step 3: Run the Analyzer

```bash
# Basic analysis (replace with your file path)
python scripts/hdock_analyzer.py --input /path/to/all_results.tar.gz --output my_results/
```

**What happens:**
```
┌─────────────────────────────────────────────────────────────┐
│  ⚙️  WHAT THE TOOL DOES                                     │
│                                                             │
│  1. 📂 Extracts your tar.gz file                            │
│  2. 📊 Reads all 100 docking models                         │
│  3. 🏆 Ranks them by score (best first)                     │
│  4. 🔬 Identifies which protein residues contact the DNA    │
│  5. 📝 Generates reports and visualization scripts          │
│                                                             │
│  ⏱️  Takes about 5-10 seconds                               │
└─────────────────────────────────────────────────────────────┘
```

---

### Step 4: View Your Results

After analysis, you'll have:

```
my_results/
├── 📄 report.html              ← OPEN THIS FIRST
├── 📄 report.md                ← For GitHub/documentation
├── 📊 hdock_analysis.json      ← For further analysis
├── 🎨 pymol_visualize.py       ← For 3D visualization
└── 📈 score_chart.md           ← Quick score overview
```

#### View HTML Report (Recommended)

```bash
# Linux
xdg-open my_results/report.html

# macOS
open my_results/report.html

# Windows
start my_results/report.html
```

**What you'll see:**
- Job information (receptor, ligand, etc.)
- Score ranking (best model highlighted)
- Binding residue analysis
- Visual charts

#### View in PyMOL (3D Visualization)

```bash
# Install PyMOL (if you don't have it)
conda install -c conda-forge pymol

# Open visualization
pymol my_results/pymol_visualize.py
```

**What you'll see:**
- Protein in gray/white
- DNA in orange
- Binding site highlighted in yellow
- Hydrogen bonds shown as dashes

---

## Understanding Your Results

### The Most Important Number: HDOCK Score

```
┌─────────────────────────────────────────────────────────────┐
│  📊 HDOCK SCORE GUIDE                                       │
│                                                             │
│  Score         │  What it means                             │
│  ──────────────│─────────────────────────────────────────── │
│  -300 to -200  │  🟢 EXCELLENT - Strong binding predicted   │
│  -200 to -150  │  🟡 GOOD - Likely real interaction         │
│  -150 to -100  │  🟠 MODERATE - May need validation         │
│  > -100        │  🔴 WEAK - Probably not binding             │
│                                                             │
│  ⚠️  IMPORTANT: These are PREDICTIONS, not measurements!    │
│     Always validate with experiments (SPR, ITC, etc.)       │
└─────────────────────────────────────────────────────────────┘
```

### What Are "Contact Residues"?

```
┌─────────────────────────────────────────────────────────────┐
│  🔬 BINDING RESIDUES EXPLAINED                              │
│                                                             │
│  When DNA binds to protein, some amino acids touch it:      │
│                                                             │
│      Protein surface        DNA                             │
│      ┌─────────────┐        ┌──────┐                        │
│      │  TRP67 ●────│────────│────  │  ← Contact residue    │
│      │             │        │      │                        │
│      │  PHE39 ●────│────────│────  │  ← Contact residue    │
│      │             │        │      │                        │
│      │  ALA52      │        └──────┘  ← Not a contact      │
│      └─────────────┘                                        │
│                                                             │
│  In top 10 models:                                          │
│    TRP67 appears 10/10 times → DEFINITELY binds             │
│    PHE39 appears 8/10 times  → LIKELY binds                 │
│    ALA52 appears 0/10 times  → Does NOT bind                │
└─────────────────────────────────────────────────────────────┘
```

---

## Common Tasks

### Compare Multiple Aptamers

```bash
# Analyze each aptamer separately
python scripts/hdock_analyzer.py --input 20mer_results.tar.gz --output analysis_20mer/
python scripts/hdock_analyzer.py --input 50mer_results.tar.gz --output analysis_50mer/
python scripts/hdock_analyzer.py --input trail_seq10.tar.gz --output analysis_seq10/

# Compare all
python scripts/compare_docks.py --results analysis_20mer/ analysis_50mer/ analysis_seq10/ --output comparison/
```

### Extract Specific Information

```python
# Python script to extract scores
import json

with open('my_results/hdock_analysis.json') as f:
    data = json.load(f)

print(f"Best score: {data['best_score']}")
print(f"Number of models: {data['total_models']}")
```

### Create Publication Figures

```bash
# Generate high-resolution PyMOL image
pymol -d "run my_results/pymol_visualize.py; save figure.png, dpi=300"
```

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| "command not found: python" | Install Python: `sudo pacman -S python` |
| "No module named numpy" | Run: `pip install numpy matplotlib` |
| "No .out file found" | Check file path, make sure you're pointing to the right folder |
| PyMOL won't open | Install: `conda install -c conda-forge pymol` |
| HTML report blank | Try different browser, or check file isn't corrupted |

---

## Getting Help

- 📧 Email: bhriguz6@gmail.com
- 🐛 Issues: https://github.com/CliffVale/dock_analysis/issues
- 📖 Docs: See full README.md

---

## Next Steps After Analysis

1. **Identify top candidate** → Look at best score
2. **Check binding residues** → Are they biologically relevant?
3. **Visualize in PyMOL** → Does the pose make sense?
4. **Validate experimentally** → Test binding affinity (SPR, ITC, fluorescence)
5. **Optimize** → Modify aptamer based on binding site insights

---

*Created by CliffVale | Part of CRP Aptamer Biosensor Project*
