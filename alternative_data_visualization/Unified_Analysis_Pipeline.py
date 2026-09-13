#!/usr/bin/env python3
"""
Unified Analysis Pipeline
A comprehensive analysis tool for biological image data processing.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch
from scipy import stats
from scipy.stats import ttest_ind
import numpy as np
import os
import re
import json
import glob
import argparse
from datetime import datetime
from openpyxl import Workbook
from openpyxl.utils.dataframe import dataframe_to_rows

# =========================
# === HELPER FUNCTIONS ====
# =========================

def normalize_path(path_str):
    """Normalize file paths to work cross-platform (Windows/Linux/macOS)
    Accepts paths with forward slashes, backslashes, or mixed, and converts to OS-appropriate format
    """
    if path_str is None:
        return None
    # os.path.normpath automatically handles path separators for the current OS
    # It converts forward slashes and backslashes to the appropriate separator
    return os.path.normpath(path_str)

def is_valid_input_path(input_path):
    """Check if an input path is valid (not None, empty, or the string "None")
    
    Args:
        input_path: The path value to check (can be string, None, or empty)
    
    Returns:
        True if the path is valid, False otherwise
    """
    if input_path is None:
        return False
    if isinstance(input_path, str):
        # Check if it's empty or the string "None" (case-insensitive)
        if input_path.strip() == "" or input_path.strip().lower() == "none":
            return False
    return True


def parse_label(filename, file_col_type="standard"):
    """Parse filename to extract CellLine, Treatment, Replicate, Photo"""
    s = os.path.basename(str(filename)).strip()
    
    if file_col_type == "alignment":
        if s.startswith("MAX_"):
            s = s[4:]
    
    s = re.sub(r'^MAX_', '', s, flags=re.I)
    s = re.sub(r'(\.(nd2|csv|tif|tiff|txt|xlsx)$|:.+)$', '', s, flags=re.I)
    s = os.path.splitext(s)[0]
    parts = s.split("_")
    
    cell_line = treatment = replicate = photo = None
    if len(parts) >= 2:
        cell_line = parts[0]
        treatment = parts[1]
        if file_col_type == "alignment":
            if len(parts) >= 4:
                try:
                    replicate = int(parts[2])
                    photo = int(parts[3])
                except ValueError:
                    replicate = photo = None
        else:
            nums = [int(p) for p in parts if re.match(r'^\d+$', p)]
            if len(nums) >= 2:
                replicate, photo = nums[-2], nums[-1]
            elif len(nums) == 1:
                replicate = nums[0]
    
    return pd.Series([cell_line, treatment, replicate, photo])

def load_data(input_file):
    """Load CSV or Excel file - handles both file types and cross-platform paths
    
    Args:
        input_file: Path to CSV or Excel file (supports .csv, .xls, .xlsx)
                   Can use forward slashes, backslashes, or mixed separators
    """
    # Normalize the path first to handle any path separator issues
    input_file = normalize_path(input_file)
    
    if input_file.lower().endswith(('.xls', '.xlsx')):
        df = pd.read_excel(input_file)
    else:
        try:
            df = pd.read_csv(input_file)
        except UnicodeDecodeError:
            df = pd.read_csv(input_file, encoding="latin1")
    return df

def get_file_column(df):
    """Determine which column contains filenames"""
    if "Label" in df.columns:
        return "Label"
    elif "File_Name" in df.columns:
        return "File_Name"
    else:
        raise ValueError("No 'Label' or 'File_Name' column found in CSV.")

def perform_statistical_test(g1_vals, g2_vals):
    """Perform independent t-test"""
    t_stat, p_val = stats.ttest_ind(g1_vals, g2_vals, equal_var=False)
    return t_stat, p_val

def get_significance_stars(p_val, alpha):
    """Get significance stars based on p-value"""
    if p_val < alpha / 100: return "****"
    elif p_val < alpha / 10: return "***"
    elif p_val < alpha / 5: return "**"
    elif p_val < alpha: return "*"
    else: return "ns"

def generate_group_colors(groups, custom_colors=None):
    """Generate distinct colors for groups automatically
    
    Args:
        groups: List of unique group names
        custom_colors: Optional dict of {group_name: color} for manual override
    
    Returns:
        Dict mapping group names to colors (hex format)
    """
    num_groups = len(groups)
    
    # Use custom colors if provided
    if custom_colors is None:
        custom_colors = {}
    
    # Select colormap based on number of groups
    # Try to use matplotlib's tab10 colormap (10 distinct colors)
    # If more than 10 groups, use Set3 (12 colors) or generate more
    if num_groups <= 10:
        cmap_name = 'tab10'
        max_colors = 10
    elif num_groups <= 12:
        cmap_name = 'Set3'
        max_colors = 12
    elif num_groups <= 20:
        cmap_name = 'tab20'
        max_colors = 20
    else:
        # For many groups, generate colors in HSV space for maximum distinctness
        cmap_name = 'hsv'
        max_colors = num_groups
    
    try:
        # Try to get colormap (works for both old and new matplotlib)
        if hasattr(plt, 'colormaps'):
            cmap = plt.colormaps.get_cmap(cmap_name)
        else:
            cmap = plt.cm.get_cmap(cmap_name)
    except:
        # Fallback to tab10 if colormap not found
        cmap = plt.cm.get_cmap('tab10')
        max_colors = 10
    
    # Generate colors for each group
    group_colors = {}
    color_idx = 0
    
    # Sort groups for consistent color assignment
    sorted_groups = sorted(groups)
    
    for group in sorted_groups:
        if group in custom_colors:
            # Use custom color if provided
            group_colors[group] = custom_colors[group]
        else:
            # Generate color from colormap
            if cmap_name == 'hsv':
                # For HSV, spread colors evenly around the color wheel (0 to 1)
                color_value = color_idx / num_groups
            else:
                # For discrete colormaps, use indices 0 to max_colors-1
                # Normalize to 0-1 range for colormap
                color_value = (color_idx % max_colors) / max(max_colors - 1, 1)
            
            color_rgba = cmap(color_value)
            
            # Convert RGBA to hex (handle both tuple and array formats)
            if isinstance(color_rgba, np.ndarray):
                r, g, b = color_rgba[0], color_rgba[1], color_rgba[2]
            else:
                r, g, b = color_rgba[0], color_rgba[1], color_rgba[2]
            
            color_hex = '#{:02x}{:02x}{:02x}'.format(
                int(r * 255),
                int(g * 255),
                int(b * 255)
            )
            group_colors[group] = color_hex
            color_idx += 1
    
    return group_colors

# =========================
# === INTENSITY COMBINING ===
# =========================
def run_intensity_combining(cfg_intensity, report_subdir=None):
    """Combine intensity ROI measurement files
    
    Args:
        cfg_intensity: Configuration dictionary for intensity combining
        report_subdir: Report subdirectory (unused, kept for compatibility)
    
    Returns:
        Path to the combined CSV file, or None if skipped
    """
    if not is_valid_input_path(cfg_intensity.get("input_folder")):
        print(f"⚠️ Intensity combining skipped: 'input_folder' not specified or set to None/blank")
        return None
    
    input_folder = normalize_path(cfg_intensity.get("input_folder"))
    output_filename = cfg_intensity.get("output_filename", "Combined_Intensity_Results.csv")
    
    if not os.path.exists(input_folder):
        print(f"⚠️ Intensity combining skipped: folder not found: {input_folder}")
        return None
    
    csv_files = glob.glob(os.path.join(input_folder, "*_ROI_measurements.csv"))
    
    if not csv_files:
        print(f"⚠️ No ROI measurement files found in {input_folder}")
        return None
    
    combined_df = pd.DataFrame()
    for f in csv_files:
        df = pd.read_csv(f)
        if df.columns[0] != "Nucleus":
            df.rename(columns={df.columns[0]: "Nucleus"}, inplace=True)
        combined_df = pd.concat([combined_df, df], ignore_index=True)
    
    output_path = os.path.join(input_folder, output_filename)
    combined_df.to_csv(output_path, index=False)
    print(f"✅ Combined intensity CSV saved to: {output_path}")
    return output_path

# =========================
# === ALIGNMENT ANALYSIS ===
# =========================
def run_alignment_analysis(cfg_analysis, report_subdir, excel_writer, analysis_name="Alignment", global_group_colors=None, global_group_assignments=None, global_stat_comparisons=None):
    """Run alignment analysis (boxplot)"""
    # Check if required fields exist and are valid
    if "input_csv" not in cfg_analysis or not is_valid_input_path(cfg_analysis.get("input_csv")):
        print(f"⚠️  Skipping {analysis_name} analysis: 'input_csv' not specified or set to None/blank")
        return None
    if "y_col" not in cfg_analysis:
        print(f"⚠️  Skipping {analysis_name} analysis: 'y_col' not specified")
        return None
    
    input_csv = normalize_path(cfg_analysis["input_csv"])
    if not os.path.exists(input_csv):
        print(f"⚠️  Skipping {analysis_name} analysis: input file not found: {input_csv}")
        return None
    
    y_col = cfg_analysis["y_col"]
    # Use section-specific group_assignments if provided, otherwise use global
    group_assignments = cfg_analysis.get("group_assignments", global_group_assignments or {})
    # Use section-specific stat_comparisons if provided, otherwise use global
    stat_comparisons = cfg_analysis.get("stat_comparisons", global_stat_comparisons or [])
    
    print(f"   📊 Input file: {input_csv}")
    print(f"   📊 Y column: {y_col}")
    if cfg_analysis.get("group_assignments"):
        print(f"   📊 Using section-specific group_assignments (overriding global)")
    else:
        print(f"   📊 Using global group_assignments")
    if cfg_analysis.get("stat_comparisons"):
        print(f"   📊 Using section-specific stat_comparisons (overriding global)")
    else:
        print(f"   📊 Using global stat_comparisons")
    
    # Load data
    df = load_data(input_csv)
    file_col = get_file_column(df)
    
    # Parse labels
    df[["CellLine","Treatment","Replicate","Photo"]] = df[file_col].apply(
        lambda x: parse_label(x, "alignment")
    )
    df = df.dropna(subset=["Treatment","Replicate"])
    df[y_col] = pd.to_numeric(df[y_col], errors="coerce")
    df = df.dropna(subset=[y_col])
    df["Group"] = df["Treatment"].map(group_assignments)
    
    # Max per image
    sum_data = df.groupby(
        ["CellLine","Treatment","Replicate","Photo"], as_index=False
    )[y_col].max()
    sum_data["Group"] = sum_data["Treatment"].map(group_assignments)
    
    # Plotting
    treatment_order = list(group_assignments.keys())
    
    # Use global color mapping if provided, otherwise generate locally
    if global_group_colors is not None:
        group_colors = global_group_colors
    else:
        treatment_colors_cfg = cfg_analysis.get("treatment_colors", {})
        unique_groups = list(set(group_assignments.values()))
        group_colors = generate_group_colors(unique_groups, treatment_colors_cfg)
    
    # Map treatments to their group colors
    treatment_colors = {
        t: group_colors.get(group_assignments[t], "gray")
        for t in group_assignments
    }
    
    plot_data = [sum_data[sum_data["Treatment"] == t][y_col].values for t in treatment_order]
    plot_colors = [treatment_colors.get(t, "gray") for t in treatment_order]
    
    plt.figure(figsize=(12, 6))
    box = plt.boxplot(plot_data, patch_artist=True, showfliers=True, tick_labels=treatment_order)
    for patch, color in zip(box["boxes"], plot_colors):
        patch.set_facecolor(color)
    
    for i, yvals in enumerate(plot_data):
        if len(yvals) > 0:
            x = np.random.normal(i + 1, 0.08, size=len(yvals))
            plt.scatter(x, yvals, color="black", s=25, zorder=2)
    
    # Statistics
    stats_cfg = cfg_analysis.get("statistics", {})
    alpha = stats_cfg.get("alpha", 0.05)
    stat_line_offset = cfg_analysis.get("stat_line_offset", 3)
    star_offset = cfg_analysis.get("star_offset", 1)
    font_size = cfg_analysis.get("font_size", 10)
    y_axis_max = cfg_analysis.get("y_axis_max", None)
    plot_title = cfg_analysis.get("plot_title", f"{analysis_name} Analysis")
    
    print(f"   ⚙️  Statistical test: Independent t-test (alpha={alpha})")
    if y_axis_max:
        print(f"   ⚙️  Y-axis max: {y_axis_max} (from config)")
    else:
        print(f"   ⚙️  Y-axis max: Auto-calculated")
    
    stats_results = []
    y_max = max([max(d) if len(d) > 0 else 0 for d in plot_data])
    
    # Check if data exceeds configured y_axis_max
    if y_axis_max and y_max > y_axis_max:
        print(f"   ⚠️  WARNING: Data maximum ({y_max:.2f}) exceeds configured y_axis_max ({y_axis_max})")
        print(f"   ⚠️  Recommendation: Increase y_axis_max to at least {y_max * 1.2:.2f} to show all data")
    
    line_offsets = np.linspace(0, stat_line_offset*y_max, num=len(stat_comparisons)+1)[1:]
    
    for i, (g1, g2) in enumerate(stat_comparisons):
        g1_vals = sum_data[sum_data["Group"] == g1][y_col].values
        g2_vals = sum_data[sum_data["Group"] == g2][y_col].values
        
        if len(g1_vals) == 0 or len(g2_vals) == 0:
            continue
        
        t_stat, p_val = perform_statistical_test(g1_vals, g2_vals)
        stars = get_significance_stars(p_val, alpha)
        
        x1 = treatment_order.index(
            sum_data[sum_data["Group"] == g1]["Treatment"].iloc[0]
        ) + 1
        x2 = treatment_order.index(
            sum_data[sum_data["Group"] == g2]["Treatment"].iloc[0]
        ) + 1
        
        y = y_max + line_offsets[i]
        plt.plot([x1,x1,x2,x2], [y, y+star_offset, y+star_offset, y], c="black")
        plt.text((x1+x2)/2, y + 0.5*star_offset, stars,
                 ha="center", va="bottom", fontsize=font_size)
        
        stats_results.append({
            "Group1": g1,
            "Group2": g2,
            "t_stat": t_stat,
            "p_value": p_val,
            "significance": stars
        })
    
    # Formatting
    plt.xticks(rotation=45, ha="right", fontsize=font_size)
    plt.ylabel(y_col)
    plt.title(plot_title)
    plt.ylim(0, y_axis_max if y_axis_max else y_max * 1.2)
    
    legend_handles = [
        Patch(color=treatment_colors[t], label=f"{t} ({group_assignments[t]})")
        for t in treatment_order
    ]
    plt.legend(handles=legend_handles, bbox_to_anchor=(1.05,1), loc="upper left")
    plt.tight_layout()
    
    # Save plot
    plot_path = os.path.join(report_subdir, f"{analysis_name}_Plot.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ {analysis_name} plot saved to: {plot_path}")
    
    # Save stats to Excel
    if stats_results:
        stats_df = pd.DataFrame(stats_results)
        stats_df.to_excel(excel_writer, sheet_name=analysis_name, index=False)
    
    return stats_results

# =========================
# === THICKNESS ANALYSIS ===
# =========================
def run_thickness_analysis(cfg_analysis, report_subdir, excel_writer, global_group_colors=None, global_group_assignments=None, global_stat_comparisons=None):
    """Run thickness analysis (boxplot)"""
    # Check if required fields exist and are valid
    if "input_csv" not in cfg_analysis or not is_valid_input_path(cfg_analysis.get("input_csv")):
        print(f"⚠️  Skipping Thickness analysis: 'input_csv' not specified or set to None/blank")
        return None
    if "y_col" not in cfg_analysis:
        print(f"⚠️  Skipping Thickness analysis: 'y_col' not specified")
        return None
    
    input_csv = normalize_path(cfg_analysis["input_csv"])
    if not os.path.exists(input_csv):
        print(f"⚠️  Skipping Thickness analysis: input file not found: {input_csv}")
        return None
    
    y_col = cfg_analysis["y_col"]
    # Use section-specific group_assignments if provided, otherwise use global
    group_assignments = cfg_analysis.get("group_assignments", global_group_assignments or {})
    # Use section-specific stat_comparisons if provided, otherwise use global
    stat_comparisons = cfg_analysis.get("stat_comparisons", global_stat_comparisons or [])
    
    # Check if group_assignments is available
    if not group_assignments:
        print(f"⚠️  Skipping Thickness analysis: 'group_assignments' not specified (needed globally or in section)")
        return None
    
    print(f"   📊 Input file: {input_csv}")
    print(f"   📊 Y column: {y_col}")
    if cfg_analysis.get("group_assignments"):
        print(f"   📊 Using section-specific group_assignments (overriding global)")
    else:
        print(f"   📊 Using global group_assignments")
    if cfg_analysis.get("stat_comparisons"):
        print(f"   📊 Using section-specific stat_comparisons (overriding global)")
    else:
        print(f"   📊 Using global stat_comparisons")
    
    df = load_data(input_csv)
    file_col = get_file_column(df)
    
    df[["CellLine","Treatment","Replicate","Photo"]] = df[file_col].apply(parse_label)
    df = df.dropna(subset=["Treatment","Replicate"])
    df[y_col] = pd.to_numeric(df[y_col], errors="coerce")
    df = df.dropna(subset=[y_col])
    df["Group"] = df["Treatment"].map(group_assignments)
    
    max_data = df.groupby(["Treatment","Group","Replicate","Photo"], as_index=False)[y_col].max()
    
    treatment_order = list(group_assignments.keys())
    
    # Use global color mapping if provided, otherwise generate locally
    if global_group_colors is not None:
        group_colors = global_group_colors
    else:
        treatment_colors_cfg = cfg_analysis.get("treatment_colors", {})
        unique_groups = list(set(group_assignments.values()))
        group_colors = generate_group_colors(unique_groups, treatment_colors_cfg)
    
    plot_data = [max_data[max_data["Treatment"]==t][y_col].values for t in treatment_order]
    plot_colors = [group_colors.get(group_assignments[t], "gray") for t in treatment_order]
    
    plt.figure(figsize=(12,6))
    box = plt.boxplot(plot_data, tick_labels=treatment_order, patch_artist=True, showfliers=True)
    for patch, color in zip(box['boxes'], plot_colors):
        patch.set_facecolor(color)
    
    for i, yvals in enumerate(plot_data):
        if len(yvals)==0: continue
        x = np.random.normal(i+1, 0.05, size=len(yvals))
        plt.scatter(x, yvals, color='black', s=25, zorder=2)
    
    # Statistics
    statistics_cfg = cfg_analysis.get("statistics", {})
    alpha = statistics_cfg.get("alpha", 0.05)
    stat_line_offset = cfg_analysis.get("stat_line_offset", 3)
    star_offset = cfg_analysis.get("star_offset", 1)
    font_size = cfg_analysis.get("font_size", 10)
    y_axis_max = cfg_analysis.get("y_axis_max", None)
    
    print(f"   ⚙️  Statistical test: Independent t-test (alpha={alpha})")
    if y_axis_max:
        print(f"   ⚙️  Y-axis max: {y_axis_max} (from config)")
    else:
        print(f"   ⚙️  Y-axis max: Auto-calculated")
    
    stats_results = []
    y_max = max([max(d) if len(d)>0 else 0 for d in plot_data])
    
    # Check if data exceeds configured y_axis_max
    if y_axis_max and y_max > y_axis_max:
        print(f"   ⚠️  WARNING: Data maximum ({y_max:.2f}) exceeds configured y_axis_max ({y_axis_max})")
        print(f"   ⚠️  Recommendation: Increase y_axis_max to at least {y_max * 1.2:.2f} to show all data")
    
    # Use a more reasonable stat_line_offset (percentage-based instead of multiplier)
    # Default to 0.1 (10% of y_max) if stat_line_offset is the old large value
    original_offset = stat_line_offset
    if stat_line_offset >= 1:
        stat_line_offset = 0.1  # Use 10% of y_max as default
        print(f"   ⚙️  Auto-adjusted stat_line_offset from {original_offset} to {stat_line_offset} (10% of y_max) for better spacing")
    
    line_offsets = np.linspace(0, stat_line_offset*y_max, num=len(stat_comparisons)+1)[1:]
    
    group_order = [group_assignments[t] for t in treatment_order]
    group_to_pos = {grp: i+1 for i, grp in enumerate(group_order)}
    
    # Calculate maximum y position needed for statistical elements
    max_stat_y = y_max
    if len(stat_comparisons) > 0:
        max_stat_y = y_max + line_offsets[-1] + star_offset
    
    for i, (g1, g2) in enumerate(stat_comparisons):
        g1_vals = max_data[max_data["Group"]==g1][y_col].values
        g2_vals = max_data[max_data["Group"]==g2][y_col].values
        
        if len(g1_vals)==0 or len(g2_vals)==0:
            continue
        
        t_stat, p_val = perform_statistical_test(g1_vals, g2_vals)
        stars = get_significance_stars(p_val, alpha)
        
        x1, x2 = group_to_pos[g1], group_to_pos[g2]
        y = y_max + line_offsets[i]
        
        plt.plot([x1,x1,x2,x2], [y, y+star_offset, y+star_offset, y], lw=1.5, c='black')
        plt.text((x1+x2)/2, y+0.5*star_offset, stars, ha='center', va='bottom', fontsize=font_size)
        
        stats_results.append({"Control": g1, "Treatment": g2, "t_stat": t_stat, "p_value": p_val, "Significance": stars})
    
    plt.xticks(rotation=45, ha='right', fontsize=font_size)
    plt.ylabel(f"{y_col} per Image")
    plt.title("Fibronectin Thickness")
    
    # Set ylim to account for statistical elements, but use y_axis_max if provided
    if y_axis_max:
        plt.ylim(0, max(y_axis_max, max_stat_y * 1.1))
    else:
        plt.ylim(0, max_stat_y * 1.1)
    
    legend_handles = [Patch(color=group_colors.get(group_assignments[t],"gray"),
                            label=f"{t} ({group_assignments[t]})") for t in treatment_order]
    plt.legend(handles=legend_handles, title="Treatment → Group", bbox_to_anchor=(1.05,1), loc="upper left")
    plt.tight_layout()
    
    plot_path = os.path.join(report_subdir, "Thickness_Plot.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Thickness plot saved to: {plot_path}")
    
    if stats_results:
        stats_df = pd.DataFrame(stats_results)
        stats_df.to_excel(excel_writer, sheet_name="Thickness", index=False)
    
    return stats_results

# =========================
# === INTENSITY ANALYSIS ===
# =========================
def run_intensity_analysis(cfg_analysis, report_subdir, excel_writer, global_group_colors=None, global_group_assignments=None, global_stat_comparisons=None):
    """Run intensity analysis (violin plot)"""
    # Check if required fields exist and are valid
    if "input_csv" not in cfg_analysis or not is_valid_input_path(cfg_analysis.get("input_csv")):
        print(f"⚠️  Skipping Intensity analysis: 'input_csv' not specified or set to None/blank")
        return None
    if "y_col" not in cfg_analysis:
        print(f"⚠️  Skipping Intensity analysis: 'y_col' not specified")
        return None
    
    input_csv = normalize_path(cfg_analysis["input_csv"])
    if not os.path.exists(input_csv):
        print(f"⚠️  Skipping Intensity analysis: input file not found: {input_csv}")
        return None
    
    y_col = cfg_analysis["y_col"]
    # Use section-specific group_assignments if provided, otherwise use global
    group_assignments = cfg_analysis.get("group_assignments", global_group_assignments or {})
    # Use section-specific stat_comparisons if provided, otherwise use global
    stat_comparisons = cfg_analysis.get("stat_comparisons", global_stat_comparisons or [])
    
    # Check if group_assignments is available
    if not group_assignments:
        print(f"⚠️  Skipping Intensity analysis: 'group_assignments' not specified (needed globally or in section)")
        return None
    
    print(f"   📊 Input file: {input_csv}")
    print(f"   📊 Y column: {y_col}")
    if cfg_analysis.get("group_assignments"):
        print(f"   📊 Using section-specific group_assignments (overriding global)")
    else:
        print(f"   📊 Using global group_assignments")
    if cfg_analysis.get("stat_comparisons"):
        print(f"   📊 Using section-specific stat_comparisons (overriding global)")
    else:
        print(f"   📊 Using global stat_comparisons")
    
    df = load_data(input_csv)
    file_col = get_file_column(df)
    
    df[["CellLine","Treatment","Replicate","Photo"]] = df[file_col].apply(parse_label)
    df = df.dropna(subset=["Treatment","Replicate"])
    df[y_col] = pd.to_numeric(df[y_col], errors="coerce")
    df = df.dropna(subset=[y_col])
    df["Group"] = df["Treatment"].map(group_assignments)
    
    treatment_order = list(group_assignments.keys())
    
    # Use global color mapping if provided, otherwise generate locally
    if global_group_colors is not None:
        group_colors = global_group_colors
    else:
        treatment_colors_cfg = cfg_analysis.get("treatment_colors", {})
        unique_groups = list(set(group_assignments.values()))
        group_colors = generate_group_colors(unique_groups, treatment_colors_cfg)
    
    # Map treatments to their group colors
    treatment_colors = {t: group_colors.get(group_assignments[t], "gray") for t in group_assignments}
    
    plot_data = [df[df["Treatment"]==t][y_col].values for t in treatment_order]
    plot_colors = [treatment_colors.get(t, "gray") for t in treatment_order]
    
    sns.set(style="whitegrid", context="talk")
    plt.figure(figsize=(12,8))
    
    for i, d in enumerate(plot_data):
        if len(d) == 0:
            continue
        violin_width = cfg_analysis.get("violin_width", 0.7)
        scatter_jitter = cfg_analysis.get("scatter_jitter", 0.045)
        parts = plt.violinplot(d, positions=[i+1], widths=violin_width,
                               showmeans=False, showmedians=False, showextrema=False)
        for patch in parts['bodies']:
            patch.set_facecolor(plot_colors[i])
            patch.set_edgecolor('black')
            patch.set_alpha(0.7)
        x = np.random.normal(i+1, scatter_jitter, size=len(d))
        plt.scatter(x, d, color='black', s=0.75, zorder=2)
    
    # Statistics
    statistics_cfg = cfg_analysis.get("statistics", {})
    alpha = statistics_cfg.get("alpha", 0.05)
    stat_line_offset = cfg_analysis.get("stat_line_offset", 0.03)
    font_size = cfg_analysis.get("font_size", 10)
    y_axis_max = cfg_analysis.get("y_axis_max", None)
    y_axis_increment = cfg_analysis.get("y_axis_increment", None)
    plot_title = cfg_analysis.get("plot_title", f"{y_col} per Image")
    
    print(f"   ⚙️  Statistical test: Independent t-test (alpha={alpha})")
    if y_axis_max:
        print(f"   ⚙️  Y-axis max: {y_axis_max} (from config)")
    else:
        print(f"   ⚙️  Y-axis max: Auto-calculated")
    
    stats_results = []
    y_max = max([max(d) if len(d) > 0 else 0 for d in plot_data])
    
    # Check if data exceeds configured y_axis_max
    if y_axis_max and y_max > y_axis_max:
        print(f"   ⚠️  WARNING: Data maximum ({y_max:.2f}) exceeds configured y_axis_max ({y_axis_max})")
        print(f"   ⚠️  Recommendation: Increase y_axis_max to at least {y_max * 1.2:.2f} to show all data")
    
    line_offsets = np.linspace(0, stat_line_offset*y_max, num=len(stat_comparisons)+1)[1:]
    
    group_order = [group_assignments[t] for t in treatment_order]
    group_to_pos = {grp: i+1 for i, grp in enumerate(group_order)}
    
    # Use y_max for vert_bar_height calculation, not y_axis_max (which can be much larger)
    vert_bar_height = 0.05 * y_max
    
    # Calculate maximum y position needed for statistical elements
    max_stat_y = y_max
    if len(stat_comparisons) > 0:
        max_stat_y = y_max + line_offsets[-1] + vert_bar_height + 0.01*y_max
    
    for i, (treat1, treat2) in enumerate(stat_comparisons):
        g1_vals = df[df["Group"]==treat1][y_col].values
        g2_vals = df[df["Group"]==treat2][y_col].values
        
        if len(g1_vals) == 0 or len(g2_vals) == 0:
            continue
        
        t_stat, p_val = perform_statistical_test(g1_vals, g2_vals)
        stars = get_significance_stars(p_val, alpha)
        
        x1 = group_to_pos[treat1]
        x2 = group_to_pos[treat2]
        y = y_max + line_offsets[i]
        
        plt.plot([x1, x1, x2, x2], [y, y+vert_bar_height, y+vert_bar_height, y], lw=1.5, c='black')
        plt.text((x1+x2)/2, y + vert_bar_height + 0.01*y_max, stars,
                 ha='center', va='bottom', fontsize=font_size)
        
        stats_results.append({"Control": treat1, "Treatment": treat2, "t_stat": t_stat, "p_value": p_val, "Significance": stars})
    
    plt.xticks(range(1, len(treatment_order)+1), treatment_order, rotation=45, ha='right', fontsize=font_size)
    plt.ylabel(y_col)
    plt.title(plot_title)
    
    legend_handles = [Patch(color=treatment_colors[t], label=f"{t} ({group_assignments[t]})") for t in treatment_order]
    plt.legend(handles=legend_handles, title="Treatment → Group", bbox_to_anchor=(1.05,1), loc="upper left")
    
    # Set ylim to account for statistical elements, but use y_axis_max if provided
    if y_axis_max:
        plt.ylim(0, max(y_axis_max, max_stat_y * 1.1))
    else:
        plt.ylim(0, max_stat_y * 1.1)
    if y_axis_increment:
        plt.yticks(np.arange(0, (y_axis_max or max_stat_y*1.1)+1, y_axis_increment))
    plt.tight_layout()
    
    plot_path = os.path.join(report_subdir, "Intensity_Plot.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Intensity plot saved to: {plot_path}")
    
    if stats_results:
        stats_df = pd.DataFrame(stats_results)
        stats_df.to_excel(excel_writer, sheet_name="Intensity", index=False)
    
    return stats_results

# =========================
# === NUCLEI COUNT ANALYSIS ===
# =========================
def run_nuclei_analysis(cfg_analysis, report_subdir, excel_writer, global_group_colors=None, global_group_assignments=None, global_stat_comparisons=None):
    """Run nuclei count analysis (boxplot)"""
    # Check if required fields exist and are valid
    if "input_csv" not in cfg_analysis or not is_valid_input_path(cfg_analysis.get("input_csv")):
        print(f"⚠️  Skipping Nuclei Counts analysis: 'input_csv' not specified or set to None/blank")
        return None
    if "y_col" not in cfg_analysis:
        print(f"⚠️  Skipping Nuclei Counts analysis: 'y_col' not specified")
        return None
    
    input_csv = normalize_path(cfg_analysis["input_csv"])
    if not os.path.exists(input_csv):
        print(f"⚠️  Skipping Nuclei Counts analysis: input file not found: {input_csv}")
        return None
    
    y_col = cfg_analysis["y_col"]
    # Use section-specific group_assignments if provided, otherwise use global
    group_assignments = cfg_analysis.get("group_assignments", global_group_assignments or {})
    # Use section-specific stat_comparisons if provided, otherwise use global
    stat_comparisons = cfg_analysis.get("stat_comparisons", global_stat_comparisons or [])
    
    # Check if group_assignments is available
    if not group_assignments:
        print(f"⚠️  Skipping Nuclei Counts analysis: 'group_assignments' not specified (needed globally or in section)")
        return None
    
    print(f"   📊 Input file: {input_csv}")
    print(f"   📊 Y column: {y_col}")
    if cfg_analysis.get("group_assignments"):
        print(f"   📊 Using section-specific group_assignments (overriding global)")
    else:
        print(f"   📊 Using global group_assignments")
    if cfg_analysis.get("stat_comparisons"):
        print(f"   📊 Using section-specific stat_comparisons (overriding global)")
    else:
        print(f"   📊 Using global stat_comparisons")
    
    df = load_data(input_csv)
    file_col = get_file_column(df)
    
    df[["CellLine","Treatment","Replicate","Photo"]] = df[file_col].apply(parse_label)
    df = df.dropna(subset=["Treatment","Replicate"])
    df[y_col] = pd.to_numeric(df[y_col], errors="coerce")
    df = df.dropna(subset=[y_col])
    df["Group"] = df["Treatment"].map(group_assignments)
    
    max_nuclei = df.groupby(["CellLine","Treatment","Replicate","Photo"], as_index=False)[y_col].max()
    max_nuclei["Group"] = max_nuclei["Treatment"].map(group_assignments)
    
    treatment_order = list(group_assignments.keys())
    
    # Use global color mapping if provided, otherwise generate locally
    if global_group_colors is not None:
        group_colors = global_group_colors
    else:
        treatment_colors_cfg = cfg_analysis.get("treatment_colors", {})
        unique_groups = list(set(group_assignments.values()))
        group_colors = generate_group_colors(unique_groups, treatment_colors_cfg)
    
    # Map treatments to their group colors
    treatment_colors = {t: group_colors.get(group_assignments[t], "gray") for t in group_assignments}
    
    plot_data = [max_nuclei[max_nuclei["Treatment"]==t][y_col].values for t in treatment_order]
    plot_colors = [treatment_colors.get(t, "gray") for t in treatment_order]
    
    plt.figure(figsize=(12,6))
    box = plt.boxplot(plot_data, patch_artist=True, showfliers=True, tick_labels=treatment_order)
    for patch, color in zip(box["boxes"], plot_colors):
        patch.set_facecolor(color)
    
    for i, yvals in enumerate(plot_data):
        if len(yvals) == 0: continue
        x = np.random.normal(i+1, 0.08, size=len(yvals))
        plt.scatter(x, yvals, color="black", s=25, zorder=2)
    
    # Statistics
    statistics_cfg = cfg_analysis.get("statistics", {})
    alpha = statistics_cfg.get("alpha", 0.05)
    stat_line_offset = cfg_analysis.get("stat_line_offset", 3)
    font_size = cfg_analysis.get("font_size", 10)
    y_axis_max = cfg_analysis.get("y_axis_max", None)
    plot_title = cfg_analysis.get("plot_title", "Nuclei per Image")
    
    print(f"   ⚙️  Statistical test: Independent t-test (alpha={alpha})")
    if y_axis_max:
        print(f"   ⚙️  Y-axis max: {y_axis_max} (from config)")
    else:
        print(f"   ⚙️  Y-axis max: Auto-calculated")
    
    stats_results = []
    y_max = max([max(d) if len(d)>0 else 0 for d in plot_data])
    
    # Check if data exceeds configured y_axis_max
    if y_axis_max and y_max > y_axis_max:
        print(f"   ⚠️  WARNING: Data maximum ({y_max:.2f}) exceeds configured y_axis_max ({y_axis_max})")
        print(f"   ⚠️  Recommendation: Increase y_axis_max to at least {y_max * 1.2:.2f} to show all data")
    
    # Use a more reasonable stat_line_offset for nuclei counts (percentage-based instead of multiplier)
    # Default to 0.1 (10% of y_max) if stat_line_offset is the old large value
    original_offset = stat_line_offset
    if stat_line_offset >= 1:
        stat_line_offset = 0.1  # Use 10% of y_max as default
        print(f"   ⚙️  Auto-adjusted stat_line_offset from {original_offset} to {stat_line_offset} (10% of y_max) for better spacing")
    
    line_offsets = np.linspace(0, stat_line_offset*y_max, num=len(stat_comparisons)+1)[1:]
    group_order = [group_assignments[t] for t in treatment_order]
    group_to_pos = {grp: i+1 for i, grp in enumerate(group_order)}
    vert_bar_height = 0.05 * (y_axis_max or y_max)
    
    # Calculate maximum y position needed for statistical elements
    max_stat_y = y_max
    if len(stat_comparisons) > 0:
        max_stat_y = y_max + line_offsets[-1] + vert_bar_height + 0.01*(y_axis_max or y_max)
    
    for i, (treat1, treat2) in enumerate(stat_comparisons):
        g1_vals = max_nuclei[max_nuclei["Group"]==treat1][y_col].values
        g2_vals = max_nuclei[max_nuclei["Group"]==treat2][y_col].values
        if len(g1_vals) == 0 or len(g2_vals) == 0:
            continue
        
        t_stat, p_val = perform_statistical_test(g1_vals, g2_vals)
        stars = get_significance_stars(p_val, alpha)
        
        x1 = group_to_pos[treat1]
        x2 = group_to_pos[treat2]
        y = y_max + line_offsets[i]
        
        plt.plot([x1, x1, x2, x2], [y, y+vert_bar_height, y+vert_bar_height, y], lw=1.5, c='black')
        plt.text((x1+x2)/2, y+vert_bar_height+0.01*(y_axis_max or y_max), stars, ha='center', va='bottom', fontsize=font_size)
        stats_results.append({"Control": treat1, "Treatment": treat2, "t_stat": t_stat, "p_value": p_val, "Significance": stars})
    
    plt.xticks(rotation=45, ha="right", fontsize=font_size)
    plt.ylabel("Nuclei per Image")
    plt.title(plot_title)
    
    # Set ylim to account for statistical elements, but use y_axis_max if provided
    if y_axis_max:
        plt.ylim(0, max(y_axis_max, max_stat_y * 1.1))
    else:
        plt.ylim(0, max_stat_y * 1.1)
    
    legend_handles = [Patch(color=treatment_colors[t], label=f"{t} ({group_assignments[t]})") for t in treatment_order]
    plt.legend(handles=legend_handles, title="Treatment → Group", bbox_to_anchor=(1.05,1), loc="upper left")
    plt.tight_layout()
    
    plot_path = os.path.join(report_subdir, "NucleiCounts_Plot.png")
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Nuclei counts plot saved to: {plot_path}")
    
    if stats_results:
        stats_df = pd.DataFrame(stats_results)
        stats_df.to_excel(excel_writer, sheet_name="NucleiCounts", index=False)
    
    return stats_results

# =========================
# === MAIN FUNCTION =======
# =========================
def main():
    """Main execution function for the Unified Analysis Pipeline"""
    parser = argparse.ArgumentParser(
        description="Unified Analysis Pipeline - Comprehensive analysis tool for biological image data processing",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python Unified_Analysis_Pipeline.py
  python Unified_Analysis_Pipeline.py --config my_config.json
  python Unified_Analysis_Pipeline.py -c custom_config.json
        """
    )
    parser.add_argument(
        '-c', '--config',
        type=str,
        default='UnifiedAnalysisInput.json',
        help='Path to configuration JSON file (default: UnifiedAnalysisInput.json)'
    )
    
    args = parser.parse_args()
    config_path = args.config
    
    # Validate config file exists
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    print("\n" + "=" * 60)
    print("UNIFIED ANALYSIS PIPELINE")
    print("=" * 60)
    
    print(f"\n📄 Loading configuration from: {config_path}")
    with open(config_path, "r") as f:
        cfg = json.load(f)
    print(f"✅ Configuration loaded successfully")
    
    # Create report directory
    report_dir = normalize_path(cfg.get("report_directory", "Analysis_Report"))
    os.makedirs(report_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_subdir = os.path.join(report_dir, f"Report_{timestamp}")
    os.makedirs(report_subdir, exist_ok=True)
    
    # Excel writer for all stats
    excel_path = os.path.join(report_subdir, "All_Statistics.xlsx")
    excel_writer = pd.ExcelWriter(excel_path, engine='openpyxl')
    
    # Run intensity combining if configured
    if "intensity_combining" in cfg:
        print("\n[1/5] Running Intensity Combining...")
        combined_intensity_path = run_intensity_combining(cfg["intensity_combining"], report_subdir)
        if combined_intensity_path and "intensity" in cfg:
            # Update intensity analysis input_csv if combining was successful
            cfg["intensity"]["input_csv"] = combined_intensity_path
    
    # =========================
    # === GET GLOBAL SETTINGS ===
    # =========================
    # Get global group_assignments if provided, otherwise collect from all sections
    global_group_assignments = cfg.get("group_assignments", None)
    
    # If global group_assignments not provided, try to get from first analysis section
    if global_group_assignments is None:
        analysis_sections = ["alignment", "thickness", "intensity", "nuclei_counts"]
        for section in analysis_sections:
            if section in cfg and "group_assignments" in cfg[section]:
                global_group_assignments = cfg[section]["group_assignments"]
                print(f"⚠️  No global 'group_assignments' found. Using assignments from '{section}' section.")
                break
    else:
        print(f"✅ Using global group_assignments: {list(global_group_assignments.keys())} → {list(global_group_assignments.values())}")
    
    # Get global stat_comparisons if provided, otherwise collect from first section
    global_stat_comparisons = cfg.get("stat_comparisons", None)
    
    # If global stat_comparisons not provided, try to get from first analysis section
    if global_stat_comparisons is None:
        analysis_sections = ["alignment", "thickness", "intensity", "nuclei_counts"]
        for section in analysis_sections:
            if section in cfg and "stat_comparisons" in cfg[section]:
                global_stat_comparisons = cfg[section]["stat_comparisons"]
                print(f"⚠️  No global 'stat_comparisons' found. Using comparisons from '{section}' section.")
                break
    else:
        print(f"✅ Using global stat_comparisons: {global_stat_comparisons}")
    
    # =========================
    # === GENERATE GLOBAL COLORS ===
    # =========================
    # Collect all unique groups across all analyses for consistent color mapping
    all_groups = set()
    all_custom_colors = {}
    
    # Use global group_assignments to get all groups
    if global_group_assignments:
        all_groups.update(global_group_assignments.values())
    
    # Also check individual sections for any additional groups or overrides
    analysis_sections = ["alignment", "thickness", "intensity", "nuclei_counts"]
    for section in analysis_sections:
        if section in cfg:
            # Check for section-specific group_assignments (may override global)
            section_assignments = cfg[section].get("group_assignments", global_group_assignments or {})
            all_groups.update(section_assignments.values())
            # Collect custom colors from all analyses (later sections override earlier ones)
            section_colors = cfg[section].get("treatment_colors", {})
            if section_colors:
                all_custom_colors.update(section_colors)
    
    # Generate global color mapping once for all groups
    if all_groups:
        global_group_colors = generate_group_colors(list(all_groups), all_custom_colors if all_custom_colors else None)
        print(f"\n🎨 Generated colors for {len(all_groups)} unique groups:")
        for group in sorted(all_groups):
            print(f"   {group}: {global_group_colors[group]}")
    else:
        global_group_colors = None
    
    # Run analyses
    analyses_run = []
    analysis_count = 0
    
    # Count total analyses to run for step numbering
    total_analyses = sum(1 for section in ["alignment", "thickness", "intensity", "nuclei_counts"] if section in cfg)
    
    if "alignment" in cfg:
        analysis_count += 1
        step_num = f"[{analysis_count}/{total_analyses}]" if total_analyses > 0 else "[1]"
        print(f"\n{step_num} Running Alignment Analysis...")
        result = run_alignment_analysis(cfg["alignment"], report_subdir, excel_writer, "Alignment", global_group_colors, global_group_assignments, global_stat_comparisons)
        if result is not None:
            analyses_run.append("Alignment")
    
    if "thickness" in cfg:
        analysis_count += 1
        step_num = f"[{analysis_count}/{total_analyses}]" if total_analyses > 0 else "[1]"
        print(f"\n{step_num} Running Thickness Analysis...")
        result = run_thickness_analysis(cfg["thickness"], report_subdir, excel_writer, global_group_colors, global_group_assignments, global_stat_comparisons)
        if result is not None:
            analyses_run.append("Thickness")
    
    if "intensity" in cfg:
        analysis_count += 1
        step_num = f"[{analysis_count}/{total_analyses}]" if total_analyses > 0 else "[1]"
        print(f"\n{step_num} Running Intensity Analysis...")
        result = run_intensity_analysis(cfg["intensity"], report_subdir, excel_writer, global_group_colors, global_group_assignments, global_stat_comparisons)
        if result is not None:
            analyses_run.append("Intensity")
    
    if "nuclei_counts" in cfg:
        analysis_count += 1
        step_num = f"[{analysis_count}/{total_analyses}]" if total_analyses > 0 else "[1]"
        print(f"\n{step_num} Running Nuclei Counts Analysis...")
        result = run_nuclei_analysis(cfg["nuclei_counts"], report_subdir, excel_writer, global_group_colors, global_group_assignments, global_stat_comparisons)
        if result is not None:
            analyses_run.append("NucleiCounts")
    
    # Close Excel writer
    excel_writer.close()
    print(f"\n✅ All statistics exported to: {excel_path}")
    
    print("\n" + "=" * 60)
    print(f"REPORT GENERATED SUCCESSFULLY")
    print("=" * 60)
    print(f"📁 Report directory: {report_subdir}")
    print(f"📊 Analyses completed: {', '.join(analyses_run) if analyses_run else 'None'}")
    print(f"📈 Statistics file: All_Statistics.xlsx")
    print(f"🎨 Color scheme: Auto-generated for {len(all_groups)} unique groups")
    print("=" * 60)


if __name__ == "__main__":
    main()
