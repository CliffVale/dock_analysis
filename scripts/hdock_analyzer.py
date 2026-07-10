#!/usr/bin/env python3
"""
hdock_analyzer.py — HDOCK Result Analyzer

Parses HDOCK output (tar.gz or extracted folder) and generates:
- Score rankings
- Binding interface analysis
- Contact residue identification
- Visualization scripts
- HTML/Markdown reports

Usage:
    python hdock_analyzer.py --input results.tar.gz --output analysis/
    python hdock_analyzer.py --input extracted_folder/ --output analysis/
"""

import argparse
import tarfile
import os
import sys
import json
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Tuple
import re


@dataclass
class DockingModel:
    """Single docking model data."""
    model_id: int
    score: float
    confidence: float = 0.0
    rmsd: float = 0.0
    ligand_rmsd: float = 0.0
    interface_residues_protein: List[str] = field(default_factory=list)
    interface_residues_dna: List[str] = field(default_factory=list)
    contacts: List[Dict] = field(default_factory=list)
    hbonds: List[Dict] = field(default_factory=list)
    stacking: List[Dict] = field(default_factory=list)
    

@dataclass
class DockingResults:
    """Complete docking results."""
    job_id: str
    receptor_file: str
    ligand_file: str
    receptor_size: Tuple[int, int, int] = (0, 0, 0)
    ligand_size: Tuple[int, int, int] = (0, 0, 0)
    models: List[DockingModel] = field(default_factory=list)
    best_model: DockingModel = None
    score_range: Tuple[float, float] = (0, 0)
    

class HDOCKAnalyzer:
    """Analyze HDOCK docking results."""
    
    def __init__(self, input_path: str, output_dir: str = "analysis"):
        self.input_path = Path(input_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results = None
        self.temp_dir = None
        
    def extract_tarball(self) -> Path:
        """Extract HDOCK tar.gz to temporary directory."""
        import tempfile
        
        if self.input_path.suffix == '.gz':
            self.temp_dir = Path(tempfile.mkdtemp())
            with tarfile.open(self.input_path, 'r:gz') as tar:
                tar.extractall(self.temp_dir)
            # Find the extracted folder
            extracted = list(self.temp_dir.glob('*'))[0]
            return extracted
        else:
            return self.input_path
    
    def parse_hdock_output(self, folder: Path) -> DockingResults:
        """Parse HDOCK output files."""
        print(f"\n[PARSING] {folder}")
        
        # Find .out file
        out_files = list(folder.glob('*.out'))
        if not out_files:
            raise FileNotFoundError(f"No .out file found in {folder}")
        
        out_file = out_files[0]
        job_id = out_file.stem.replace('hdock_', '')
        
        results = DockingResults(
            job_id=job_id,
            receptor_file="",
            ligand_file=""
        )
        
        # Parse .out file
        with open(out_file, 'r') as f:
            lines = f.readlines()
        
        # Parse header (grid info, receptor/ligand sizes)
        for i, line in enumerate(lines):
            if 'rec_' in line:
                parts = line.split()
                results.receptor_file = parts[0]
                results.receptor_size = (float(parts[1]), float(parts[2]), float(parts[3]))
            elif 'lig_' in line:
                parts = line.split()
                results.ligand_file = parts[0]
                results.ligand_size = (float(parts[1]), float(parts[2]), float(parts[3]))
            elif re.match(r'^\s*[\d\-]', line) and len(line.split()) >= 7:
                # This is a model line
                parts = line.split()
                if len(parts) >= 7:
                    try:
                        model = DockingModel(
                            model_id=len(results.models) + 1,
                            confidence=float(parts[1]) if len(parts) > 1 else 0,
                            rmsd=float(parts[2]) if len(parts) > 2 else 0,
                            score=float(parts[6]) if len(parts) > 6 else 0
                        )
                        results.models.append(model)
                    except (ValueError, IndexError):
                        continue
        
        if results.models:
            results.best_model = results.models[0]
            scores = [m.score for m in results.models]
            results.score_range = (min(scores), max(scores))
        
        print(f"  Job ID: {job_id}")
        print(f"  Models: {len(results.models)}")
        if results.best_model:
            print(f"  Best score: {results.best_model.score}")
            print(f"  Score range: {results.score_range}")
        
        self.results = results
        return results
    
    def parse_interfaces(self, folder: Path):
        """Parse interface_*.txt files for contact details."""
        print(f"\n[INTERFACES] Parsing contact files...")
        
        interface_files = sorted(folder.glob('interface_*.txt'))
        
        for i, iface_file in enumerate(interface_files[:10]):  # Top 10 models
            model_id = i + 1
            if model_id > len(self.results.models):
                break
            
            model = self.results.models[model_id - 1]
            
            with open(iface_file, 'r') as f:
                content = f.read()
            
            # Parse contacts
            for line in content.split('\n'):
                if line.strip() and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 4:
                        contact = {
                            'receptor_atom': parts[0],
                            'receptor_residue': parts[1] if len(parts) > 1 else '',
                            'dna_atom': parts[2] if len(parts) > 2 else '',
                            'dna_residue': parts[3] if len(parts) > 3 else '',
                            'distance': float(parts[4]) if len(parts) > 4 else 0
                        }
                        model.contacts.append(contact)
            
            print(f"  Model {model_id}: {len(model.contacts)} contacts")
    
    def analyze_binding_interface(self):
        """Identify key binding residues and interactions."""
        print(f"\n[ANALYSIS] Binding interface analysis...")
        
        if not self.results or not self.results.models:
            print("  No models to analyze")
            return
        
        # Aggregate contacts across top models
        residue_frequency = {}
        atom_contacts = {}
        
        for model in self.results.models[:10]:  # Top 10
            for contact in model.contacts:
                res = contact.get('receptor_residue', '')
                if res:
                    residue_frequency[res] = residue_frequency.get(res, 0) + 1
        
        # Sort by frequency
        sorted_residues = sorted(residue_frequency.items(), key=lambda x: x[1], reverse=True)
        
        print(f"\n  Top binding residues (frequency in top 10 models):")
        for res, freq in sorted_residues[:10]:
            print(f"    {res}: {freq}/10 models")
        
        return sorted_residues
    
    def generate_pymol_script(self):
        """Generate PyMOL visualization script."""
        print(f"\n[VIZ] Generating PyMOL script...")
        
        script_content = f"""
# Auto-generated PyMOL script for HDOCK results
# Job ID: {self.results.job_id}

# Load structures
load hdock_receptor.pdb, receptor
load hdock_ligand.pdb, ligand
load model_1.pdb, best_model

# Style receptor
hide all
show cartoon, receptor
color gray80, receptor
show surface, receptor
set transparency, 0.7, receptor

# Style ligand (DNA)
show sticks, ligand
color orange, ligand
show spheres, ligand and name P

# Style best model
show cartoon, best_model
color spectrum, best_model

# Highlight binding site
select binding_site, best_model and resi {"+".join(self._get_binding_residues())}
show sticks, binding_site
color yellow, binding_site

# Distance contacts
distance hbonds, best_model, receptor, mode=2
hide labels, hbonds

# View
center best_model
zoom best_model
orient

# Save
save {self.results.job_id}_visualization.pse
save {self.results.job_id}_best_model.png, dpi=300
"""
        
        script_path = self.output_dir / "pymol_visualize.py"
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        print(f"  Saved: {script_path}")
    
    def _get_binding_residues(self) -> List[str]:
        """Get binding residue numbers from best model."""
        residues = set()
        if self.results.best_model:
            for contact in self.results.best_model.contacts:
                res = contact.get('receptor_residue', '')
                match = re.search(r'(\d+)', res)
                if match:
                    residues.add(match.group(1))
        return list(residues)[:20]  # Limit for PyMOL
    
    def generate_report(self, format='markdown'):
        """Generate analysis report."""
        print(f"\n[REPORT] Generating {format} report...")
        
        if format == 'markdown':
            self._generate_markdown_report()
        elif format == 'html':
            self._generate_html_report()
    
    def _generate_markdown_report(self):
        """Generate markdown report."""
        report = f"""# HDOCK Docking Analysis Report

## Job Information
| Property | Value |
|----------|-------|
| Job ID | `{self.results.job_id}` |
| Receptor | {self.results.receptor_file} |
| Ligand | {self.results.ligand_file} |
| Total Models | {len(self.results.models)} |

## Scoring Summary
| Metric | Value |
|--------|-------|
| Best Score | **{self.results.best_model.score if self.results.best_model else 'N/A'}** |
| Score Range | {self.results.score_range[0]} to {self.results.score_range[1]} |
| Mean Score | {sum(m.score for m in self.results.models) / len(self.results.models) if self.results.models else 0:.2f} |

## Top 10 Models
| Rank | Model | Score | Confidence |
|------|-------|-------|------------|
"""
        
        for model in self.results.models[:10]:
            report += f"| {model.model_id} | model_{model.model_id}.pdb | {model.score} | {model.confidence} |\n"
        
        report += f"""
## Binding Interface Analysis

### Key Contact Residues
Top residues identified across top 10 models:
"""
        
        sorted_residues = self.analyze_binding_interface()
        if sorted_residues:
            for res, freq in sorted_residues[:10]:
                report += f"- **{res}**: Present in {freq}/10 models\n"
        
        report += f"""
## Visualization

PyMOL script saved: `pymol_visualize.py`

To view:
```bash
pymol pymol_visualize.py
```

## Files Generated
- `hdock_analysis.json` — Machine-readable results
- `pymol_visualize.py` — PyMOL visualization script
- `binding_residues.csv` — Binding residue details
- `report.md` — This report

---
*Generated by dock_analysis*
*CliffVale (bhriguz6@gmail.com)*
"""
        
        report_path = self.output_dir / "report.md"
        with open(report_path, 'w') as f:
            f.write(report)
        
        print(f"  Saved: {report_path}")
    
    def _generate_html_report(self):
        """Generate HTML report with embedded styling."""
        html = """<!DOCTYPE html>
<html>
<head>
    <title>HDOCK Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; background: #f5f5f5; }}
        .container {{ max-width: 1000px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
        th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background: #3498db; color: white; }}
        tr:hover {{ background: #f5f5f5; }}
        .best {{ background: #d4edda; font-weight: bold; }}
        .score {{ font-family: monospace; font-size: 1.1em; }}
        code {{ background: #f4f4f4; padding: 2px 6px; border-radius: 3px; }}
        .info-box {{ background: #e8f4fc; border-left: 4px solid #3498db; padding: 15px; margin: 15px 0; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>HDOCK Docking Analysis Report</h1>
        
        <div class="info-box">
            <strong>What is this?</strong> This report analyzes results from HDOCK protein-DNA docking.
            HDOCK generates 100 docking models ranked by score. Lower (more negative) scores indicate better binding.
        </div>
        
        <h2>Job Information</h2>
        <table>
            <tr><th>Property</th><th>Value</th></tr>
            <tr><td>Job ID</td><td><code>JOB_ID</code></td></tr>
            <tr><td>Receptor</td><td>RECEPTOR</td></tr>
            <tr><td>Ligand</td><td>LIGAND</td></tr>
            <tr><td>Total Models</td><td>MODELS_COUNT</td></tr>
        </table>
        
        <h2>Scoring Summary</h2>
        <div class="info-box">
            <strong>Understanding Scores:</strong> HDOCK scores are binding energies (kcal/mol). 
            More negative = stronger predicted binding. Typical range: -100 to -300.
        </div>
        <table>
            <tr><th>Metric</th><th>Value</th></tr>
            <tr class="best"><td>Best Score</td><td class="score">BEST_SCORE</td></tr>
            <tr><td>Score Range</td><td class="score">RANGE</td></tr>
            <tr><td>Mean Score</td><td class="score">MEAN</td></tr>
        </table>
        
        <h2>Top 10 Models</h2>
        <table>
            <tr><th>Rank</th><th>Model</th><th>Score</th></tr>
            TOP_MODELS
        </table>
        
        <h2>Binding Interface</h2>
        <div class="info-box">
            <strong>What are contact residues?</strong> These are amino acids on the protein surface 
            that directly interact with the DNA. They're identified by counting how often each residue 
            appears in the top 10 docking models.
        </div>
        BINDING_RESIDUES
        
        <hr>
        <p><em>Generated by dock_analysis | CliffVale (bhriguz6@gmail.com)</em></p>
    </div>
</body>
</html>"""
        
        # Replace placeholders
        html = html.replace('JOB_ID', self.results.job_id)
        html = html.replace('RECEPTOR', self.results.receptor_file)
        html = html.replace('LIGAND', self.results.ligand_file)
        html = html.replace('MODELS_COUNT', str(len(self.results.models)))
        
        if self.results.best_model:
            html = html.replace('BEST_SCORE', str(self.results.best_model.score))
        html = html.replace('RANGE', f"{self.results.score_range[0]} to {self.results.score_range[1]}")
        
        mean_score = sum(m.score for m in self.results.models) / len(self.results.models) if self.results.models else 0
        html = html.replace('MEAN', f"{mean_score:.2f}")
        
        # Top models table
        top_models_html = ""
        for model in self.results.models[:10]:
            top_models_html += f"""
            <tr>
                <td>{model.model_id}</td>
                <td>model_{model.model_id}.pdb</td>
                <td class="score">{model.score}</td>
            </tr>"""
        html = html.replace('TOP_MODELS', top_models_html)
        
        # Binding residues
        sorted_residues = self.analyze_binding_interface()
        residues_html = "<ul>"
        for res, freq in sorted_residues[:10]:
            residues_html += f"<li><strong>{res}</strong> — present in {freq}/10 models</li>"
        residues_html += "</ul>"
        html = html.replace('BINDING_RESIDUES', residues_html)
        
        html_path = self.output_dir / "report.html"
        with open(html_path, 'w') as f:
            f.write(html)
        
        print(f"  Saved: {html_path}")
    
    def save_json(self):
        """Save results as JSON for programmatic access."""
        data = {
            'job_id': self.results.job_id,
            'receptor': self.results.receptor_file,
            'ligand': self.results.ligand_file,
            'total_models': len(self.results.models),
            'best_score': self.results.best_model.score if self.results.best_model else None,
            'score_range': list(self.results.score_range),
            'models': [
                {
                    'id': m.model_id,
                    'score': m.score,
                    'confidence': m.confidence,
                    'rmsd': m.rmsd,
                    'num_contacts': len(m.contacts)
                }
                for m in self.results.models
            ]
        }
        
        json_path = self.output_dir / "hdock_analysis.json"
        with open(json_path, 'w') as f:
            json.dump(data, f, indent=2)
        
        print(f"  Saved: {json_path}")
    
    def run(self):
        """Run complete analysis pipeline."""
        print("\n" + "="*60)
        print("HDOCK RESULT ANALYZER")
        print("="*60)
        
        # Extract if tarball
        folder = self.extract_tarball()
        
        # Parse results
        self.parse_hdock_output(folder)
        self.parse_interfaces(folder)
        
        # Analyze
        self.analyze_binding_interface()
        
        # Generate outputs
        self.save_json()
        self.generate_pymol_script()
        self.generate_report('markdown')
        self.generate_report('html')
        
        print("\n" + "="*60)
        print("ANALYSIS COMPLETE!")
        print("="*60)
        print(f"\nResults in: {self.output_dir}/")
        print("\nNext steps:")
        print("  1. Open report.html in browser")
        print("  2. Run pymol_visualize.py in PyMOL")
        print("  3. Use hdock_analysis.json for further analysis")


def main():
    parser = argparse.ArgumentParser(
        description="HDOCK Result Analyzer",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Analyze tar.gz from HDOCK
    python hdock_analyzer.py --input results.tar.gz --output analysis/
    
    # Analyze extracted folder
    python hdock_analyzer.py --input hdock_output/ --output analysis/
    
    # Generate only specific outputs
    python hdock_analyzer.py --input results.tar.gz --output analysis/ --format html
        """
    )
    
    parser.add_argument("--input", required=True,
                        help="HDOCK output (tar.gz or folder)")
    parser.add_argument("--output", default="analysis",
                        help="Output directory (default: analysis)")
    parser.add_argument("--format", choices=['markdown', 'html', 'both'],
                        default='both', help="Report format")
    
    args = parser.parse_args()
    
    try:
        analyzer = HDOCKAnalyzer(args.input, args.output)
        analyzer.run()
    except Exception as e:
        print(f"\n[ERROR] {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
