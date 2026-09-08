#!/usr/bin/env python3
"""
haddock_analyzer.py — HADDOCK Result Analyzer

Parses HADDOCK 2.4 output and generates:
- Score rankings and cluster analysis
- Binding interface analysis
- Contact residue identification
- Visualization scripts
- HTML/Markdown reports

Usage:
    python haddock_analyzer.py --input haddock_output/ --output analysis/
    python haddock_analyzer.py --input run1/ run2/ --output comparison/
"""

import argparse
import os
import sys
import json
import re
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Tuple, Optional
import shutil


@dataclass
class HADDOCKModel:
    """Single HADDOCK docking model."""
    model_id: int
    cluster_id: int = 0
    score: float = 0.0
    vdw: float = 0.0
    elec: float = 0.0
    desolv: float = 0.0
    bsa: float = 0.0
    rmsd: float = 0.0
    air: float = 0.0
    interface_residues_protein: List[str] = field(default_factory=list)
    interface_residues_dna: List[str] = field(default_factory=list)
    contacts: List[Dict] = field(default_factory=list)
    hbonds: List[Dict] = field(default_factory=list)


@dataclass
class HADDOCKCluster:
    """Cluster of HADDOCK models."""
    cluster_id: int
    models: List[HADDOCKModel] = field(default_factory=list)
    best_score: float = 0.0
    mean_score: float = 0.0
    size: int = 0


@dataclass
class HADDOCKResults:
    """Complete HADDOCK results."""
    run_id: str
    receptor_file: str = ""
    ligand_file: str = ""
    models: List[HADDOCKModel] = field(default_factory=list)
    clusters: List[HADDOCKCluster] = field(default_factory=list)
    best_model: Optional[HADDOCKModel] = None
    score_range: Tuple[float, float] = (0, 0)


class HADDOCKAnalyzer:
    """Analyze HADDOCK 2.4 docking results."""

    def __init__(self, input_path: str, output_dir: str = "analysis"):
        self.input_path = Path(input_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results = None

    def parse_haddock_output(self, folder: Path) -> HADDOCKResults:
        """
        Parse HADDOCK output directory.
        Looks for: structures/struc_*.pdb, the_score.list, etc.
        """
        print(f"\n[PARSING] {folder}")

        run_id = folder.name

        results = HADDOCKResults(run_id=run_id)

        # Look for the_score.list (main score file)
        score_file = folder / "the_score.list"
        if score_file.exists():
            self._parse_score_list(score_file, results)
        else:
            # Try to find score files in subdirectories
            for sub in folder.rglob("the_score.list"):
                self._parse_score_list(sub, results)
                break

        # Look for individual model PDBs
        pdb_files = sorted(folder.rglob("struc_*.pdb"))
        if not pdb_files:
            pdb_files = sorted(folder.rglob("model_*.pdb"))

        # Parse PDB files for contact information
        self._parse_model_contacts(folder, results, pdb_files)

        # Cluster models
        self._cluster_models(results)

        if results.models:
            results.best_model = results.models[0]
            scores = [m.score for m in results.models]
            results.score_range = (min(scores), max(scores))

        print(f"  Run ID: {run_id}")
        print(f"  Models: {len(results.models)}")
        print(f"  Clusters: {len(results.clusters)}")
        if results.best_model:
            print(f"  Best score: {results.best_model.score}")

        self.results = results
        return results

    def _parse_score_list(self, score_file: Path, results: HADDOCKResults):
        """Parse the_score.list file."""
        print(f"  Parsing {score_file}")

        with open(score_file) as f:
            lines = f.readlines()

        # HADDOCK score file format:
        # Rank  Model  Score  vdw  elec  desolv  air  bsa  ...
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith('-----'):
                continue

            parts = line.split()
            if len(parts) >= 6:
                try:
                    # Try to identify fields by position
                    # Common format: rank, model_name, score, vdw, elec, desolv, air, bsa
                    rank = int(parts[0]) if parts[0].isdigit() else len(results.models) + 1

                    # Find numeric values
                    numbers = []
                    for p in parts:
                        try:
                            numbers.append(float(p))
                        except ValueError:
                            continue

                    if len(numbers) >= 3:
                        model = HADDOCKModel(
                            model_id=rank,
                            score=numbers[2] if len(numbers) > 2 else numbers[-1],
                            vdw=numbers[3] if len(numbers) > 3 else 0,
                            elec=numbers[4] if len(numbers) > 4 else 0,
                            desolv=numbers[5] if len(numbers) > 5 else 0,
                            air=numbers[6] if len(numbers) > 6 else 0,
                            bsa=numbers[7] if len(numbers) > 7 else 0,
                        )
                        results.models.append(model)
                except (ValueError, IndexError):
                    continue

    def _parse_model_contacts(self, folder: Path, results: HADDOCKResults,
                               pdb_files: List[Path]):
        """Parse PDB files for contact information."""
        print(f"  Parsing {len(pdb_files)} model PDBs for contacts...")

        cutoff = 4.5  # Angstroms for contact detection

        for pdb_file in pdb_files[:20]:  # Top 20 models
            try:
                # Extract model number from filename
                match = re.search(r'struc_(\d+)|model_(\d+)', pdb_file.name)
                if not match:
                    continue
                model_id = int(match.group(1) or match.group(2))

                # Find matching model in results
                model = None
                for m in results.models:
                    if m.model_id == model_id:
                        model = m
                        break

                if not model:
                    # Create model entry if not from score list
                    model = HADDOCKModel(model_id=model_id, score=0)
                    results.models.append(model)

                # Parse contacts from PDB
                protein_atoms = []
                dna_atoms = []

                with open(pdb_file) as f:
                    for line in f:
                        if line.startswith(('ATOM', 'HETATM')):
                            x = float(line[30:38])
                            y = float(line[38:46])
                            z = float(line[46:54])
                            atom_name = line[12:16].strip()
                            resname = line[17:20].strip()
                            resid = line[22:26].strip()
                            chain = line[21]

                            atom_info = {
                                'name': atom_name,
                                'resname': resname,
                                'resid': resid,
                                'chain': chain,
                                'coords': [x, y, z]
                            }

                            # Separate protein vs DNA
                            if resname in ('DA', 'DT', 'DC', 'DG', 'A', 'T', 'C', 'G',
                                          'RA', 'RU', 'RC', 'RG', 'U'):
                                dna_atoms.append(atom_info)
                            else:
                                protein_atoms.append(atom_info)

                # Find contacts between protein and DNA
                import numpy as np
                contacts = []
                hbonds = []

                for p_atom in protein_atoms:
                    for d_atom in dna_atoms:
                        dist = np.sqrt(sum(
                            (a - b) ** 2
                            for a, b in zip(p_atom['coords'], d_atom['coords'])
                        ))

                        if dist < cutoff:
                            contact = {
                                'protein_atom': p_atom['name'],
                                'protein_residue': f"{p_atom['resname']}{p_atom['resid']}",
                                'protein_chain': p_atom['chain'],
                                'dna_atom': d_atom['name'],
                                'dna_residue': f"{d_atom['resname']}{d_atom['resid']}",
                                'dna_chain': d_atom['chain'],
                                'distance': round(dist, 2)
                            }
                            contacts.append(contact)

                            # H-bond detection (donor-acceptor distance < 3.5A)
                            if dist < 3.5 and (
                                p_atom['name'] in ('N', 'O', 'NZ', 'ND1', 'NE2', 'OG', 'OD1', 'OD2') or
                                d_atom['name'] in ('O4', 'O6', 'N3', 'N2', 'N4')
                            ):
                                hbonds.append(contact)

                model.contacts = contacts
                model.hbonds = hbonds
                model.interface_residues_protein = list(set(
                    c['protein_residue'] for c in contacts
                ))
                model.interface_residues_dna = list(set(
                    c['dna_residue'] for c in contacts
                ))

            except Exception as e:
                print(f"    [!] Error parsing {pdb_file.name}: {e}")
                continue

    def _cluster_models(self, results: HADDOCKResults):
        """Cluster models by score proximity."""
        if not results.models:
            return

        # Simple clustering: group by score proximity
        sorted_models = sorted(results.models, key=lambda m: m.score)
        clusters = []
        current_cluster = []
        threshold = 5.0  # Score difference threshold

        for model in sorted_models:
            if not current_cluster:
                current_cluster.append(model)
            elif abs(model.score - current_cluster[0].score) < threshold:
                current_cluster.append(model)
            else:
                clusters.append(current_cluster)
                current_cluster = [model]

        if current_cluster:
            clusters.append(current_cluster)

        # Create cluster objects
        for i, cluster_models in enumerate(clusters, 1):
            scores = [m.score for m in cluster_models]
            cluster = HADDOCKCluster(
                cluster_id=i,
                models=cluster_models,
                best_score=min(scores),
                mean_score=sum(scores) / len(scores),
                size=len(cluster_models)
            )
            for m in cluster_models:
                m.cluster_id = i
            results.clusters.append(cluster)

    def analyze_binding_interface(self):
        """Identify key binding residues across all clusters."""
        print(f"\n[ANALYSIS] Binding interface analysis...")

        if not self.results or not self.results.models:
            print("  No models to analyze")
            return

        # Aggregate contacts across top models (top 10 per cluster)
        residue_frequency = {}
        hbond_frequency = {}

        for model in self.results.models[:10]:
            for contact in model.contacts:
                res = contact.get('protein_residue', '')
                if res:
                    residue_frequency[res] = residue_frequency.get(res, 0) + 1

            for hb in model.hbonds:
                key = f"{hb['protein_residue']}-{hb['dna_residue']}"
                hbond_frequency[key] = hbond_frequency.get(key, 0) + 1

        # Sort by frequency
        sorted_residues = sorted(
            residue_frequency.items(), key=lambda x: x[1], reverse=True
        )
        sorted_hbonds = sorted(
            hbond_frequency.items(), key=lambda x: x[1], reverse=True
        )

        print(f"\n  Top binding residues (frequency in top models):")
        for res, freq in sorted_residues[:10]:
            print(f"    {res}: {freq}/{min(10, len(self.results.models))} models")

        print(f"\n  Top H-bonds:")
        for hb, freq in sorted_hbonds[:5]:
            print(f"    {hb}: {freq}/{min(10, len(self.results.models))} models")

        return sorted_residues, sorted_hbonds

    def generate_pymol_script(self):
        """Generate PyMOL visualization script."""
        print(f"\n[VIZ] Generating PyMOL script...")

        binding_residues = self._get_binding_residues()

        script_content = f"""#!/usr/bin/env python3
# Auto-generated PyMOL script for HADDOCK results
# Run ID: {self.results.run_id}

# Load structures
load haddock_receptor.pdb, receptor
load haddock_ligand.pdb, ligand

# Load top cluster models
"""
        # Load first model from each cluster
        for cluster in self.results.clusters[:3]:
            if cluster.models:
                best = cluster.models[0]
                script_content += f"load struc_{best.model_id}.pdb, cluster{cluster.cluster_id}\n"

        script_content += """
# Style receptor
hide all
show cartoon, receptor
color gray80, receptor
show surface, receptor
set transparency, 0.7, receptor

# Style ligand (DNA/RNA)
show sticks, ligand
color orange, ligand
show spheres, ligand and name P

# Style cluster models
"""
        colors = ["cyan", "magenta", "yellow"]
        for i, cluster in enumerate(self.results.clusters[:3]):
            color = colors[i % len(colors)]
            script_content += f"""
show cartoon, cluster{cluster.cluster_id}
color {color}, cluster{cluster.cluster_id}
"""

        if binding_residues:
            script_content += f"""
# Highlight binding site
select binding_site, receptor and resi {"+".join(binding_residues[:20])}
show sticks, binding_site
color yellow, binding_site
"""

        script_content += """
# Distance contacts
distance hbonds, receptor, ligand, mode=2
hide labels, hbonds

# View
center receptor
zoom receptor
orient

# Save session
save haddock_analysis.pse
save haddock_best_model.png, dpi=300
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
                res = contact.get('protein_residue', '')
                match = re.search(r'(\d+)', res)
                if match:
                    residues.add(match.group(1))
        return list(residues)[:20]

    def generate_report(self, format: str = 'markdown'):
        """Generate analysis report."""
        print(f"\n[REPORT] Generating {format} report...")

        if format == 'markdown':
            self._generate_markdown_report()
        elif format == 'html':
            self._generate_html_report()

    def _generate_markdown_report(self):
        """Generate markdown report."""
        report = f"""# HADDOCK Docking Analysis Report

## Run Information
| Property | Value |
|----------|-------|
| Run ID | `{self.results.run_id}` |
| Receptor | {self.results.receptor_file} |
| Ligand | {self.results.ligand_file} |
| Total Models | {len(self.results.models)} |
| Clusters | {len(self.results.clusters)} |

## Scoring Summary
| Metric | Value |
|--------|-------|
| Best Score | **{self.results.best_model.score if self.results.best_model else 'N/A'}** |
| Score Range | {self.results.score_range[0]:.2f} to {self.results.score_range[1]:.2f} |
| Mean Score | {sum(m.score for m in self.results.models) / len(self.results.models) if self.results.models else 0:.2f} |

## HADDOCK Score Components
| Component | Best Model |
|-----------|-----------|
| Total Score | {self.results.best_model.score if self.results.best_model else 'N/A'} |
| van der Waals | {self.results.best_model.vdw if self.results.best_model else 'N/A'} |
| Electrostatics | {self.results.best_model.elec if self.results.best_model else 'N/A'} |
| Desolvation | {self.results.best_model.desolv if self.results.best_model else 'N/A'} |
| AIR (restraints) | {self.results.best_model.air if self.results.best_model else 'N/A'} |
| BSA | {self.results.best_model.bsa if self.results.best_model else 'N/A'} |

## Cluster Analysis
| Cluster | Size | Best Score | Mean Score |
|---------|------|------------|------------|
"""
        for cluster in self.results.clusters:
            report += f"| {cluster.cluster_id} | {cluster.size} | {cluster.best_score:.2f} | {cluster.mean_score:.2f} |\n"

        report += f"""
## Top 10 Models
| Rank | Cluster | Model | Score | vdw | elec | desolv |
|------|---------|-------|-------|-----|------|--------|
"""
        for model in self.results.models[:10]:
            report += (
                f"| {model.model_id} | {model.cluster_id} | "
                f"struc_{model.model_id}.pdb | {model.score:.2f} | "
                f"{model.vdw:.2f} | {model.elec:.2f} | {model.desolv:.2f} |\n"
            )

        report += f"""
## Binding Interface Analysis

### Key Contact Residues
Top residues identified across top models:
"""
        result = self.analyze_binding_interface()
        if result:
            sorted_residues = result[0]
            for res, freq in sorted_residues[:10]:
                report += f"- **{res}**: Present in {freq}/{min(10, len(self.results.models))} models\n"

        report += f"""
## Visualization

PyMOL script saved: `pymol_visualize.py`

To view:
```bash
pymol pymol_visualize.py
```

## Files Generated
- `haddock_analysis.json` — Machine-readable results
- `pymol_visualize.py` — PyMOL visualization script
- `binding_residues.csv` — Binding residue details
- `report.md` — This report

---
*Generated by dock_analysis (HADDOCK analyzer)*
*CliffVale (bhriguz6@gmail.com)*
"""

        report_path = self.output_dir / "report.md"
        with open(report_path, 'w') as f:
            f.write(report)

        print(f"  Saved: {report_path}")

    def _generate_html_report(self):
        """Generate HTML report."""
        html = f"""<!DOCTYPE html>
<html>
<head>
    <title>HADDOCK Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; background: #f5f5f5; }}
        .container {{ max-width: 1000px; margin: 0 auto; background: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; border-bottom: 2px solid #e74c3c; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        table {{ width: 100%; border-collapse: collapse; margin: 20px 0; }}
        th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background: #e74c3c; color: white; }}
        tr:hover {{ background: #f5f5f5; }}
        .best {{ background: #d4edda; font-weight: bold; }}
        .score {{ font-family: monospace; font-size: 1.1em; }}
        code {{ background: #f4f4f4; padding: 2px 6px; border-radius: 3px; }}
        .info-box {{ background: #fce4ec; border-left: 4px solid #e74c3c; padding: 15px; margin: 15px 0; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>HADDOCK Docking Analysis Report</h1>

        <div class="info-box">
            <strong>HADDOCK vs HDOCK:</strong> This report analyzes HADDOCK 2.4 results.
            HADDOCK uses data-driven docking with restraints, producing HADDOCK scores
            (incorporating vdw, elec, desolvation, BSA). Lower = better binding.
        </div>

        <h2>Run Information</h2>
        <table>
            <tr><th>Property</th><th>Value</th></tr>
            <tr><td>Run ID</td><td><code>{self.results.run_id}</code></td></tr>
            <tr><td>Total Models</td><td>{len(self.results.models)}</td></tr>
            <tr><td>Clusters</td><td>{len(self.results.clusters)}</td></tr>
        </table>

        <h2>Scoring Summary</h2>
        <table>
            <tr><th>Metric</th><th>Value</th></tr>
            <tr class="best"><td>Best Score</td><td class="score">{self.results.best_model.score if self.results.best_model else 'N/A'}</td></tr>
            <tr><td>Score Range</td><td class="score">{self.results.score_range[0]:.2f} to {self.results.score_range[1]:.2f}</td></tr>
        </table>

        <h2>Cluster Analysis</h2>
        <table>
            <tr><th>Cluster</th><th>Size</th><th>Best Score</th><th>Mean Score</th></tr>
"""
        for cluster in self.results.clusters:
            html += f"            <tr><td>{cluster.cluster_id}</td><td>{cluster.size}</td><td>{cluster.best_score:.2f}</td><td>{cluster.mean_score:.2f}</td></tr>\n"

        html += """        </table>

        <h2>Top 10 Models</h2>
        <table>
            <tr><th>Rank</th><th>Cluster</th><th>Score</th><th>vdw</th><th>elec</th></tr>
"""
        for model in self.results.models[:10]:
            html += f"            <tr><td>{model.model_id}</td><td>{model.cluster_id}</td><td class='score'>{model.score:.2f}</td><td>{model.vdw:.2f}</td><td>{model.elec:.2f}</td></tr>\n"

        html += """        </table>

        <hr>
        <p><em>Generated by dock_analysis (HADDOCK analyzer) | CliffVale (bhriguz6@gmail.com)</em></p>
    </div>
</body>
</html>"""

        html_path = self.output_dir / "report.html"
        with open(html_path, 'w') as f:
            f.write(html)

        print(f"  Saved: {html_path}")

    def save_json(self):
        """Save results as JSON."""
        data = {
            'run_id': self.results.run_id,
            'receptor': self.results.receptor_file,
            'ligand': self.results.ligand_file,
            'total_models': len(self.results.models),
            'total_clusters': len(self.results.clusters),
            'best_score': self.results.best_model.score if self.results.best_model else None,
            'score_range': list(self.results.score_range),
            'clusters': [
                {
                    'id': c.cluster_id,
                    'size': c.size,
                    'best_score': c.best_score,
                    'mean_score': c.mean_score,
                }
                for c in self.results.clusters
            ],
            'models': [
                {
                    'id': m.model_id,
                    'cluster': m.cluster_id,
                    'score': m.score,
                    'vdw': m.vdw,
                    'elec': m.elec,
                    'desolv': m.desolv,
                    'bsa': m.bsa,
                    'num_contacts': len(m.contacts),
                    'num_hbonds': len(m.hbonds),
                }
                for m in self.results.models
            ]
        }

        json_path = self.output_dir / "haddock_analysis.json"
        with open(json_path, 'w') as f:
            json.dump(data, f, indent=2)

        print(f"  Saved: {json_path}")

    def save_contact_csv(self):
        """Save contact residues as CSV."""
        csv_path = self.output_dir / "binding_residues.csv"
        with open(csv_path, 'w') as f:
            f.write("Residue,Frequency,Models\n")

            residue_freq = {}
            for model in self.results.models[:10]:
                for contact in model.contacts:
                    res = contact.get('protein_residue', '')
                    if res:
                        residue_freq[res] = residue_freq.get(res, 0) + 1

            for res, freq in sorted(residue_freq.items(), key=lambda x: x[1], reverse=True):
                f.write(f"{res},{freq},{freq}/{min(10, len(self.results.models))}\n")

        print(f"  Saved: {csv_path}")

    def run(self):
        """Run complete analysis pipeline."""
        print("\n" + "="*60)
        print("HADDOCK RESULT ANALYZER")
        print("="*60)

        # Handle multiple input directories
        if self.input_path.is_dir():
            folders = [self.input_path]
        else:
            folders = [self.input_path]

        # Parse results from each folder
        for folder in folders:
            self.parse_haddock_output(folder)

        # Analyze
        self.analyze_binding_interface()

        # Generate outputs
        self.save_json()
        self.save_contact_csv()
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
        print("  3. Use haddock_analysis.json for further analysis")


def main():
    parser = argparse.ArgumentParser(
        description="HADDOCK 2.4 Result Analyzer",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Analyze single HADDOCK run
    python haddock_analyzer.py --input run1/ --output analysis/

    # Compare multiple runs
    python haddock_analyzer.py --input run1/ run2/ run3/ --output comparison/

    # Generate only HTML report
    python haddock_analyzer.py --input run1/ --output analysis/ --format html
        """
    )

    parser.add_argument("--input", "-i", required=True, nargs='+',
                        help="HADDOCK output directory(ies)")
    parser.add_argument("--output", "-o", default="analysis",
                        help="Output directory (default: analysis)")
    parser.add_argument("--format", choices=['markdown', 'html', 'both'],
                        default='both', help="Report format")

    args = parser.parse_args()

    try:
        # Analyze each input directory
        for input_dir in args.input:
            analyzer = HADDOCKAnalyzer(input_dir, args.output)
            analyzer.run()

    except Exception as e:
        print(f"\n[ERROR] {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
