#!/usr/bin/env python3
"""
Generate comprehensive workflow report
Combines MultiQC and DESeq2 results into a final report
"""
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import base64
import io

def main():
    # Get inputs from Snakemake
    multiqc_file = snakemake.input["multiqc"]
    deseq2_file = snakemake.input["deseq2"]
    output_file = snakemake.output[0]
    print("Generating comprehensive workflow report...")

    # Load DESeq2 results
    de_results = pd.read_csv(deseq2_file)

    # Generate summary statistics
    total_genes = len(de_results)
    significant_genes = len(de_results[
        (de_results['padj'] < 0.05) &
        (abs(de_results['log2FoldChange']) > 1)
    ])

    upregulated = len(de_results[
        (de_results['padj'] < 0.05) &
        (de_results['log2FoldChange'] > 1)
    ])

    downregulated = len(de_results[
        (de_results['padj'] < 0.05) &
        (de_results['log2FoldChange'] < -1)
    ])

    # Create summary plot
    summary_plot = create_summary_plot(upregulated, downregulated, total_genes)

    # Generate HTML report
    html_content = f"""


    <div>
        <h1>RNA-seq Analysis Report</h1>
        <p>Comprehensive analysis results from raw reads to differential expression</
    </div>

    <div>
        <h2>Analysis Summary</h2>
        <div>
            <div>
                <div>{total_genes:,}</div>
                <div>Total Genes</div>
            </div>

            <div>
                <div>{significant_genes:,}</div>
                <div>Significant DEGs</div>
            </div>

            <div>
                <div>{upregulated:,}</div>
                <div>Upregulated</div>
            </div>

            <div>
                <div>{downregulated:,}</div>
                <div>Downregulated</div>
            </div>
        </div>
    </div>

    <div>
        <h2>Differential Expression Overview</h2>
        <img>
    </div>

    <div>
        <h2>Quality Control</h2>
        <p>Detailed quality control metrics are available in the <a href="../qc/multi
    </div>

    <div>
        <h2>Analysis Details</h2>
        <ul>
            <li>Alignment performed with STAR</li>
            <li>Gene quantification with featureCounts</li>
            <li>Differential expression analysis with DESeq2</li>
            <li>Significance criteria: padj &lt; 0.05 and |log2FC| &gt; 1</li>
        </ul>
    </div>

    &lt;/body&gt;
    &lt;/html&gt;
    """

    # Save report
    with open(output_file, 'w') as f:
        f.write(html_content)
        print(f"Report generated successfully: {output_file}")


def create_summary_plot(upregulated, downregulated, total_genes):
    """Create a summary plot of differential expression results"""
    fig, ax = plt.subplots(figsize=(8, 6))
    categories = ['Upregulated', 'Downregulated', 'Not Significant']
    values = [upregulated, downregulated, total_genes - upregulated - downregulated]
    colors = ['#ff9999', '#66b3ff', '#99ff99']
    wedges, texts, autotexts = ax.pie(values, labels=categories, colors=colors,
    autopct='%1.1f%%', startangle=90)
    ax.set_title('Differential Expression Results', fontsize=14, fontweight='bold')

    # Convert plot to base64 string
    buffer = io.BytesIO()
    plt.savefig(buffer, format='png', bbox_inches='tight', dpi=300)
    buffer.seek(0)
    plot_data = buffer.getvalue()
    buffer.close()
    plt.close()
    return base64.b64encode(plot_data).decode()

if __name__ == "__main__":
    main()
