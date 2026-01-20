import os
import re
import shutil
import json

# ============================
# Load Configuration
# ============================

config_path = "RawFileRenamingInput.json"
if not os.path.exists(config_path):
    raise FileNotFoundError(f"Config file not found: {config_path}")

print(f"📄 Loading configuration from: {config_path}")
with open(config_path, "r") as f:
    cfg = json.load(f)
print(f"✅ Configuration loaded successfully\n")

# Normalize paths for cross-platform compatibility
def normalize_path(path_str):
    """Convert path string to use os.path.normpath for cross-platform compatibility"""
    if path_str is None:
        return None
    return os.path.normpath(path_str)

# Extract settings from config
input_root = normalize_path(cfg.get("input_root"))
output_root = normalize_path(cfg.get("output_root"))
well_map_raw = cfg.get("well_map", {})

# Convert well_map from JSON format (lists) to tuples for compatibility
well_map = {}
for well_id, mapping in well_map_raw.items():
    if isinstance(mapping, list) and len(mapping) >= 3:
        well_map[well_id] = (mapping[0], mapping[1], mapping[2])
    elif isinstance(mapping, dict):
        well_map[well_id] = (
            mapping.get("cell_line", ""),
            mapping.get("treatment", ""),
            mapping.get("replicate", "")
        )
    else:
        print(f"⚠️  Warning: Invalid mapping format for {well_id}, skipping")

if not input_root:
    raise ValueError("'input_root' not specified in config file")
if not output_root:
    raise ValueError("'output_root' not specified in config file")
if not os.path.exists(input_root):
    raise FileNotFoundError(f"Input directory not found: {input_root}")

os.makedirs(output_root, exist_ok=True)

print(f"📁 Input directory: {input_root}")
print(f"📁 Output directory: {output_root}")
print(f"📋 Well mappings: {len(well_map)} wells\n")

# ============================
# File Renaming Logic
# ============================

renamed_count = 0
skipped_no_well = 0
skipped_unmapped_well = 0

for root, dirs, files in os.walk(input_root):

    for file in files:
        old_path = os.path.join(root, file)
        ext = os.path.splitext(file)[1]

        # --- Extract well ID (WellA1, WellC04, etc.) ---
        well_match = re.search(r"(Well[A-H]\d{1,2})", file)
        if not well_match:
            skipped_no_well += 1
            print(f"⚠️  Skipped (no well ID): {file}")
            continue

        well_id = well_match.group(1)

        # Normalize WellA1 ➜ WellA01
        if re.match(r"Well[A-H]\d$", well_id):
            letter = well_id[4]
            number = well_id[5]
            well_id = f"Well{letter}0{number}"

        if well_id not in well_map:
            skipped_unmapped_well += 1
            print(f"⚠️  Skipped (well not in map: {well_id}): {file}")
            continue

        cell_line, treatment, replicate = well_map[well_id]

        # --- Extract image number ---
        img_match = re.search(r"_(\d{4})_", file)
        image_num = img_match.group(1) if img_match else "0000"

        # --- Build new filename ---
        new_name = f"{cell_line}_{treatment}_{replicate}_{image_num}{ext}"

        # --- Rebuild directory inside output_root (Option E) ---
        relative_path = os.path.relpath(root, input_root)
        new_dir = os.path.join(output_root, relative_path)
        os.makedirs(new_dir, exist_ok=True)

        new_path = os.path.join(new_dir, new_name)

        # --- Copy (safe) rather than overwrite ---
        shutil.copy2(old_path, new_path)
        renamed_count += 1

        print(f"✅ {file} → {new_name}")

# ============================
# Summary
# ============================

print("\n================ Summary ================")
print(f"✔ Renamed files: {renamed_count}")
print(f"⚠️ Skipped (no well ID): {skipped_no_well}")
print(f"⚠️ Skipped (well not in map): {skipped_unmapped_well}")
print("Done!")
