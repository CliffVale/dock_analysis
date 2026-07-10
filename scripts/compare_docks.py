#!/usr/bin/env python3
"""
compare_docks.py — Compare Multiple Docking Results

Compare HDOCK scores across multiple aptamers/targets.
Generates comparison tables, bar charts, and summary reports.

Usage:
    python compare_docks.py --results results_20mer/ results_50mer/ --output comparison/
"""

import argparse
import os
import sys
import json
from pathlib import Path
from typing import List, Dict
from dataclasses import dataclass


@dataclass
class DockingSummary:
    """Summary of a single docking run."""
    name: str
    job_id: str
    best_score: float
    mean_score: float
    score_range: tuple
    num_models: int
    top_residues: List[str]
    

class DockingComparator:
    """Compare multiple docking results."""
    
    def __init__(self, result_dirs: List[str], output_dir: str = "comparison"):
        self.result_dirs = [Path(d) for d in result_dirs]
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.summaries = []
        
    def load_results(self):
        """Load JSON results from each docking run."""
        print("\n[LOADING] Loading docking results...")
        
        for result_dir in self.result_dirs:
            json_file = result_dir / "hdock_analysis.json"
            if json_file.exists():
                with open(json_file, 'r') as f:
                    data = json.load(f)
                
                summary = DockingSummary(
                    name=result_dir.name,
                    job_id=data.get('job_id', 'unknown'),
                    best_score=data.get('best_score', 0),
                    mean_score=sum(m['score'] for m in data.get('models', [])) / max(len(data.get('models', [])), 1),
                    score_range=tuple(data.get('score_range', [0, 0])),
                    num_models=data.get('total_models', 0),
                    top_residues=[]
                )
                self.summaries.append(summary)
                print(f"  ✓ {summary.name}: {summary.best_score}")
            else:
                print(f"  ✗ {result_dir.name}: No hdock_analysis.json found")
    
    def compare_scores(self):
        """Generate score comparison table."""
        print("\n[COMPARE] Score comparison...")
        
        # Sort by best score (most negative = best)
        sorted_results = sorted(self.summaries, key=lambda x: x.best_score)
        
        table = """
# Docking Score Comparison

## Ranking (Best to Worst)

| Rank | Aptamer/Target | Best Score | Mean Score | Range |
|------|----------------|------------|------------|-------|
"""
        
        for i, s in enumerate(sorted_results, 1):
            table += f"| {i} | {s.name} | **{s.best_score}** | {s.mean_score:.2f} | {s.score_range[0]} to {s.score_range[1]} |\n"
        
        # Add per-nucleotide efficiency if applicable
        table += """
## Per-Nucleotide Efficiency

| Aptamer | Length | Score | Score/nt |
|---------|--------|-------|----------|
"""
        
        for s in sorted_results:
            # Try to extract length from name
            length = self._extract_length(s.name)
            if length:
                score_per_nt = s.best_score / length
                table += f"| {s.name} | {length} nt | {s.best_score} | {score_per_nt:.2f} |\n"
        
        # Save
        report_path = self.output_dir / "comparison_report.md"
        with open(report_path, 'w') as f:
            f.write(table)
        
        print(f"  Saved: {report_path}")
        
        # Also save as JSON
        comparison_data = {
            'results': [
                {
                    'name': s.name,
                    'best_score': s.best_score,
                    'mean_score': s.mean_score,
                    'score_range': list(s.score_range),
                    'length': self._extract_length(s.name),
                    'score_per_nt': s.best_score / self._extract_length(s.name) if self._extract_length(s.name) else None
                }
                for s in sorted_results
            ]
        }
        
        json_path = self.output_dir / "comparison.json"
        with open(json_path, 'w') as f:
            json.dump(comparison_data, f, indent=2)
        
        print(f"  Saved: {json_path}")
    
    def _extract_length(self, name: str) -> int:
        """Extract aptamer length from name (e.g., '20mer' → 20)."""
        import re
        match = re.search(r'(\d+)\s*mer', name, re.IGNORECASE)
        if match:
            return int(match.group(1))
        return None
    
    def generate_bar_chart(self):
        """Generate ASCII bar chart for quick visualization."""
        print("\n[CHART] Generating bar chart...")
        
        sorted_results = sorted(self.summaries, key=lambda x: x.best_score)
        
        # Normalize scores for display
        max_abs_score = max(abs(s.best_score) for s in sorted_results) if sorted_results else 1
        
        chart = "\n# Score Comparison Chart\n\n"
        chart += "```\n"
        chart += f"{'Aptamer':<20} {'Score':>10} {'Bar':<40}\n"
        chart += "-" * 70 + "\n"
        
        for s in sorted_results:
            bar_length = int((abs(s.best_score) / max_abs_score) * 35)
            bar = "█" * bar_length
            chart += f"{s.name:<20} {s.best_score:>10.2f} {bar}\n"
        
        chart += "```\n"
        chart += "\n*More negative (longer bar) = better binding*\n"
        
        chart_path = self.output_dir / "score_chart.md"
        with open(chart_path, 'w') as f:
            f.write(chart)
        
        print(f"  Saved: {chart_path}")
    
    def generate_summary(self):
        """Generate executive summary."""
        print("\n[SUMMARY] Generating summary...")
        
        if not self.summaries:
            print("  No results to summarize")
            return
        
        sorted_results = sorted(self.summaries, key=lambda x: x.best_score)
        best = sorted_results[0]
        
        summary = f"""
# Executive Summary

## Best Performer: {best.name}

| Metric | Value |
|--------|-------|
| Best Score | **{best.best_score}** |
| Mean Score | {best.mean_score:.2f} |
| Score Range | {best.score_range[0]} to {best.score_range[1]} |

## All Results

| Aptamer | Best Score | Status |
|---------|------------|--------|
"""
        
        for s in sorted_results:
            status = "🏆 BEST" if s == best else "✓"
            summary += f"| {s.name} | {s.best_score} | {status} |\n"
        
        summary += f"""
## Recommendation

Based on docking scores, **{best.name}** shows the strongest predicted binding 
with a score of **{best.best_score}**.

### Next Steps
1. Validate top candidate experimentally
2. Test binding affinity (SPR, ITC, or fluorescence)
3. Characterize binding kinetics (kon, koff)
4. Test in relevant biological context

---
*Generated by dock_analysis | CliffVale (bhriguz6@gmail.com)*
"""
        
        summary_path = self.output_dir / "EXECUTIVE_SUMMARY.md"
        with open(summary_path, 'w') as f:
            f.write(summary)
        
        print(f"  Saved: {summary_path}")
    
    def run(self):
        """Run complete comparison pipeline."""
        print("\n" + "="*60)
        print("DOCKING RESULT COMPARATOR")
        print("="*60)
        
        self.load_results()
        
        if len(self.summaries) < 2:
            print("\n[WARNING] Need at least 2 results for comparison")
            print("Running single-result summary instead...")
        
        self.compare_scores()
        self.generate_bar_chart()
        self.generate_summary()
        
        print("\n" + "="*60)
        print("COMPARISON COMPLETE!")
        print("="*60)
        print(f"\nResults in: {self.output_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description="Compare Multiple Docking Results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Compare two docking runs
    python compare_docks.py --results results_20mer/ results_50mer/ --output comparison/
    
    # Compare multiple runs
    python compare_docks.py --results dock1/ dock2/ dock3/ --output comparison/
        """
    )
    
    parser.add_argument("--results", nargs='+', required=True,
                        help="Directories containing hdock_analysis.json")
    parser.add_argument("--output", default="comparison",
                        help="Output directory (default: comparison)")
    
    args = parser.parse_args()
    
    try:
        comparator = DockingComparator(args.results, args.output)
        comparator.run()
    except Exception as e:
        print(f"\n[ERROR] {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
