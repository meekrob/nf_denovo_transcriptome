import os 
import pandas as pd 
import numpy as np 
from datetime import datetime
import matplotlib.pyplot as plt
from tabulate import tabulate

# Configuration
path = "/nfs/home/jlamb/Projects/"
dir_list = [
    "01_CulexTarsalisLifeStages",
    "02_CxTarMidgutSingleCellReads",
    "03_CxTarOvarySingleCellReads",
    "04_emily_ebel_bulk_midgut_data_eg_RNAseq",
]

# Initialize DataFrames to store results
file_stats_df = pd.DataFrame(columns=['Dataset', 'Stage', 'Files', 'Size (GB)'])
busco_summary_df = pd.DataFrame(columns=['Dataset', 'Complete%', 'Single Copy%', 'Duplicated%', 'Fragmented%', 'Missing%', 'Total BUSCOs'])
busco_details_df = pd.DataFrame(columns=['Dataset', 'Complete', 'Single Copy', 'Duplicated', 'Fragmented', 'Missing', 'Total'])

# Function to convert bytes to GB
def bytes_to_gb(size_bytes):
    return round(size_bytes / 1_000_000_000, 2)

# Function to count paired-end files
def count_paired_files(directory, extension=".fastq.gz"):
    return len([f for f in os.listdir(directory) if f.endswith(extension)]) / 2

# Function to get total size of files
def get_total_size(directory, extension=".fastq.gz"):
    return sum([os.path.getsize(os.path.join(directory, f)) for f in os.listdir(directory) if f.endswith(extension)])

# Process each dataset
for directory in dir_list:
    dataset_path = os.path.join(path, directory)
    results_path = os.path.join(dataset_path, "results")
    trimming_path = os.path.join(results_path, "trimming")
    seqkit_qc_path = os.path.join(results_path, "seqkit_qc")
    merging_path = os.path.join(results_path, "merging")
    normalization_path = os.path.join(results_path, "normalization")
    assembly_path = os.path.join(results_path, "assembly", "rnaspades", "rnaspades")
    busco_path = os.path.join(results_path, "busco", "busco_output")

    # Extract directory shortname for display
    dir_short = directory.split("_")[0]
    
    # Collect file statistics for raw files
    try:
        raw_files = count_paired_files(results_path)
        raw_size = bytes_to_gb(get_total_size(results_path))
        file_stats_df = file_stats_df.append({
            'Dataset': dir_short,
            'Stage': 'Raw',
            'Files': int(raw_files),
            'Size (GB)': raw_size
        }, ignore_index=True)
    except:
        pass

    # Collect file statistics for QC files
    try:
        qc_files = count_paired_files(seqkit_qc_path)
        qc_size = bytes_to_gb(get_total_size(seqkit_qc_path))
        file_stats_df = file_stats_df.append({
            'Dataset': dir_short,
            'Stage': 'QC',
            'Files': int(qc_files),
            'Size (GB)': qc_size
        }, ignore_index=True)
    except:
        pass

    # Collect file statistics for merged files
    try:
        merged_files = []
        for f in os.listdir(merging_path):
            if f.endswith(".fastq.gz"):
                file_size = bytes_to_gb(os.path.getsize(os.path.join(merging_path, f)))
                merged_files.append({
                    'Dataset': dir_short,
                    'Stage': 'Merged: ' + f,
                    'Files': 1,
                    'Size (GB)': file_size
                })
        for entry in merged_files:
            file_stats_df = file_stats_df.append(entry, ignore_index=True)
    except:
        pass

    # Collect file statistics for normalized files
    try:
        norm_files = []
        for f in os.listdir(normalization_path):
            if f.endswith(".fastq.gz"):
                file_size = bytes_to_gb(os.path.getsize(os.path.join(normalization_path, f)))
                norm_files.append({
                    'Dataset': dir_short,
                    'Stage': 'Normalized: ' + f,
                    'Files': 1,
                    'Size (GB)': file_size
                })
        for entry in norm_files:
            file_stats_df = file_stats_df.append(entry, ignore_index=True)
    except:
        pass

    # Collect assembly statistics
    try:
        assembly_size = bytes_to_gb(os.path.getsize(os.path.join(assembly_path, "transcripts.fasta")))
        file_stats_df = file_stats_df.append({
            'Dataset': dir_short,
            'Stage': 'Assembly',
            'Files': 1,
            'Size (GB)': assembly_size
        }, ignore_index=True)
    except:
        pass

    # Extract BUSCO statistics
    try:
        busco_stats = {}
        with open(os.path.join(busco_path, "short_summary.specific.diptera_odb10.busco_output.txt"), "r") as f:
            for line in f:
                line = line.strip()
                # Extract lineage dataset info
                if "The lineage dataset is:" in line:
                    busco_stats["lineage"] = line.split("The lineage dataset is:")[1].strip()
                
                # Extract overall metrics
                elif line.startswith("C:") and "]" in line:
                    parts = line.split(",")
                    complete_part = parts[0]  # C:41.4%[S:36.5%,D:4.9%]
                    fragmented_part = parts[1]  # F:15.6%
                    missing_part = parts[2].split(",n:")[0]  # M:43.0%
                    total_buscos = parts[2].split("n:")[1].strip()  # 3285
                    
                    busco_stats["complete_pct"] = float(complete_part.split(":")[1].split("[")[0])
                    busco_stats["single_copy_pct"] = float(complete_part.split("S:")[1].split("%")[0])
                    busco_stats["duplicated_pct"] = float(complete_part.split("D:")[1].split("%")[0])
                    busco_stats["fragmented_pct"] = float(fragmented_part.split(":")[1].split("%")[0])
                    busco_stats["missing_pct"] = float(missing_part.split(":")[1].split("%")[0])
                    busco_stats["total_buscos"] = int(total_buscos)
                
                # Extract detailed counts
                elif "Complete BUSCOs" in line and "(" in line:
                    busco_stats["complete"] = int(line.split()[0])
                elif "Complete and single-copy BUSCOs" in line:
                    busco_stats["single_copy"] = int(line.split()[0])
                elif "Complete and duplicated BUSCOs" in line:
                    busco_stats["duplicated"] = int(line.split()[0])
                elif "Fragmented BUSCOs" in line:
                    busco_stats["fragmented"] = int(line.split()[0])
                elif "Missing BUSCOs" in line:
                    busco_stats["missing"] = int(line.split()[0])
                elif "Total BUSCO groups searched" in line:
                    busco_stats["total"] = int(line.split()[0])

        # Add BUSCO summary to DataFrames
        if busco_stats:
            busco_summary_df = busco_summary_df.append({
                'Dataset': dir_short,
                'Complete%': busco_stats.get("complete_pct", "N/A"),
                'Single Copy%': busco_stats.get("single_copy_pct", "N/A"),
                'Duplicated%': busco_stats.get("duplicated_pct", "N/A"),
                'Fragmented%': busco_stats.get("fragmented_pct", "N/A"),
                'Missing%': busco_stats.get("missing_pct", "N/A"),
                'Total BUSCOs': busco_stats.get("total_buscos", "N/A")
            }, ignore_index=True)
            
            busco_details_df = busco_details_df.append({
                'Dataset': dir_short,
                'Complete': busco_stats.get("complete", "N/A"),
                'Single Copy': busco_stats.get("single_copy", "N/A"),
                'Duplicated': busco_stats.get("duplicated", "N/A"),
                'Fragmented': busco_stats.get("fragmented", "N/A"),
                'Missing': busco_stats.get("missing", "N/A"),
                'Total': busco_stats.get("total", "N/A")
            }, ignore_index=True)
    except:
        pass

# Generate the report
def generate_report():
    # Title and metadata
    now = datetime.now()
    report = []
    report.append("=" * 80)
    report.append(f"TRANSCRIPTOME ASSEMBLY PIPELINE SUMMARY REPORT")
    report.append(f"Generated on: {now.strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("=" * 80)
    report.append("")
    
    # File Statistics Table
    report.append("-" * 80)
    report.append("FILE STATISTICS")
    report.append("-" * 80)
    if not file_stats_df.empty:
        report.append(tabulate(file_stats_df, headers='keys', tablefmt='grid', showindex=False))
    else:
        report.append("No file statistics available.")
    report.append("")
    
    # BUSCO Summary Table
    report.append("-" * 80)
    report.append("BUSCO ASSESSMENT SUMMARY (PERCENTAGES)")
    report.append("-" * 80)
    if not busco_summary_df.empty:
        # Format percentages
        summary_formatted = busco_summary_df.copy()
        for col in ['Complete%', 'Single Copy%', 'Duplicated%', 'Fragmented%', 'Missing%']:
            summary_formatted[col] = summary_formatted[col].apply(lambda x: f"{x:.1f}%" if isinstance(x, (int, float)) else x)
        report.append(tabulate(summary_formatted, headers='keys', tablefmt='grid', showindex=False))
    else:
        report.append("No BUSCO summary statistics available.")
    report.append("")
    
    # BUSCO Details Table
    report.append("-" * 80)
    report.append("BUSCO ASSESSMENT DETAILS (COUNTS)")
    report.append("-" * 80)
    if not busco_details_df.empty:
        report.append(tabulate(busco_details_df, headers='keys', tablefmt='grid', showindex=False))
    else:
        report.append("No BUSCO detailed statistics available.")
    
    return "\n".join(report)

# Print the report
print(generate_report())

# Write report to file
with open('transcriptome_assembly_report.txt', 'w') as f:
    f.write(generate_report())

# Generate visualizations
try:
    if not busco_summary_df.empty:
        # BUSCO completeness bar chart
        plt.figure(figsize=(10, 8))
        busco_plot_data = busco_summary_df.set_index('Dataset')
        busco_plot_data[['Complete%', 'Fragmented%', 'Missing%']].plot(kind='bar', stacked=True, 
                                                                         color=['green', 'orange', 'red'])
        plt.title('BUSCO Assessment by Dataset')
        plt.ylabel('Percentage')
        plt.tight_layout()
        plt.savefig('busco_assessment.png')
        
    if not file_stats_df.empty:
        # File size comparison
        plt.figure(figsize=(12, 8))
        size_data = file_stats_df.pivot(index='Dataset', columns='Stage', values='Size (GB)')
        size_data = size_data.loc[:, ~size_data.columns.str.contains('Merged|Normalized')]
        size_data.plot(kind='bar')
        plt.title('File Size Comparison by Processing Stage')
        plt.ylabel('Size (GB)')
        plt.tight_layout()
        plt.savefig('file_size_comparison.png')
except Exception as e:
    print(f"Error generating plots: {e}")







