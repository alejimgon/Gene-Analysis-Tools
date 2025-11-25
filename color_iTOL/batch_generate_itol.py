#!/usr/bin/env python3
"""
Batch script to generate iTOL files for all tree files in a directory.

Usage:
    1. Ensure you have the required script `generate_itol_taxonomy_colors.py` in the same directory.
    2. Update the paths in the script as necessary.
    3. Run this script to process all `.treefile` files in the specified directory.
"""

import os
import glob
import subprocess
import sys

def main():
    # Configuration
    trees_dir = "BASE_DIR"  # /path/to/your/trees_directory
    taxonomy_file = "BASE_DIR" # /path/to/your/taxonomy_file.txt
    output_dir = "BASE_DIR" # /path/to/your/output_directory
    
    # Script path - should be in the same directory as this batch script
    script_path = os.path.join(output_dir, "generate_itol_taxonomy_colors.py")
    
    # Check if required files exist
    if not os.path.exists(script_path):
        print(f"Error: Script not found at {script_path}")
        sys.exit(1)
    
    if not os.path.exists(taxonomy_file):
        print(f"Error: Taxonomy file not found at {taxonomy_file}")
        sys.exit(1)
    
    if not os.path.exists(trees_dir):
        print(f"Error: Trees directory not found at {trees_dir}")
        sys.exit(1)
    
    # Find all .treefile files
    tree_files = glob.glob(os.path.join(trees_dir, "*.treefile"))
    
    if not tree_files:
        print(f"No .treefile files found in {trees_dir}")
        sys.exit(1)
    
    print(f"Found {len(tree_files)} tree files to process")
    print(f"Trees directory: {trees_dir}")
    print(f"Taxonomy file: {taxonomy_file}")
    print(f"Output directory: {output_dir}")
    print(f"Script: {script_path}")
    print()
    
    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    success_count = 0
    
    for tree_file in sorted(tree_files):  # Sort for consistent ordering
        try:
            # Extract gene name from filename
            basename = os.path.basename(tree_file)
            # Handle different naming patterns
            if '_verification_' in basename:
                gene_name = basename.split('_verification_')[0]
            elif '_filtered_' in basename:
                gene_name = basename.split('_filtered_')[0] 
            else:
                gene_name = basename.split('_')[0]  # fallback
            
            output_file = os.path.join(output_dir, f"{gene_name}_itol_colors.txt")
            labels_file = os.path.join(output_dir, f"{gene_name}_itol_labels.txt")
            
            print(f"Processing {gene_name}...")
            
            # Run the script with all the options from your successful run
            cmd = [
                sys.executable, script_path,
                tree_file,
                taxonomy_file,
                "-o", output_file,
                "-l", "phylum",
                "--labels",  # Generate labels file
                "--save-rooted",  # Save rooted trees
                "-v"  # Verbose output
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, cwd=output_dir)
            
            if result.returncode == 0:
                print(f"  ✓ Generated {os.path.basename(output_file)}")
                print(f"  ✓ Generated {os.path.basename(labels_file)}")
                if "--save-rooted" in cmd:
                    rooted_file = os.path.join(output_dir, f"{gene_name}_midpoint.treefile")
                    if os.path.exists(rooted_file):
                        print(f"  ✓ Generated {os.path.basename(rooted_file)}")
                success_count += 1
            else:
                print(f"  ✗ Failed processing {gene_name}:")
                if result.stderr:
                    print(f"    Error: {result.stderr.strip()}")
                if result.stdout:
                    print(f"    Output: {result.stdout.strip()}")
                
        except Exception as e:
            print(f"  ✗ Error processing {os.path.basename(tree_file)}: {e}")
    
    print(f"\nCompleted: {success_count}/{len(tree_files)} files processed successfully")
    
    if success_count > 0:
        print(f"\nGenerated files are in: {output_dir}")
        print("Files generated for each gene:")
        print("  - {gene}_itol_colors.txt (main coloring file)")
        print("  - {gene}_itol_labels.txt (external labels)")
        print("  - {gene}_midpoint.treefile (rooted tree)")

if __name__ == "__main__":
    main()