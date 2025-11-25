#!/usr/bin/env python3
"""
Script to automatically generate iTOL DATASET_RANGE files based on taxonomy and phylogenetic trees.
This script identifies taxonomic groups and creates consistent range-based coloring with automatic legends.

Usage:
    python generate_itol_taxonomy_colors.py <tree_file> <taxonomy_file> -o <output_file> -l <taxonomic_level> [-v] [--labels]

    - iTOL troubleshooting problem with the legend: even when specifying legend in vertical, it sometimes ignores it.
        - When the dataset is loaded, manually set legend to vertical in the iTOL interface.
            - Dataset -> Editing functions -> Legend -> Advanced options -> Direction -> Vertical -> Update dataset legend

Arguments:
    tree_file: Path to the Newick tree file.
    taxonomy_file: Path to the taxonomy TSV file.
    -o, --output: Output iTOL file path (default: auto-generated).
    -l, --level: Taxonomic level for grouping (default: phylum). Options: superkingdom, phylum, class, order, family.
    -v, --verbose: Verbose output.
    --labels: Generate external text labels file.
    --root-method: Tree rooting method. Options: midpoint, outgroup, none (default: midpoint).
    --outgroup: Outgroup sequence name for outgroup rooting (required if --root-method outgroup).
    --save-rooted: Save the rooted tree to a new file for reference.

Requires:
    - pandas
    - biopython
"""

import pandas as pd
import argparse
import os
import sys
import re
from collections import defaultdict
from Bio import Phylo

def apply_tree_rooting(tree, root_method='midpoint', outgroup=None, verbose=False):
    """
    Apply rooting to the tree to ensure consistent topology for monophyletic analysis.
    """
    if not tree:
        return tree
    
    if root_method == 'none':
        if verbose:
            print("  No rooting applied - using tree as-is")
        return tree
    
    elif root_method == 'midpoint':
        if verbose:
            print("  Applying midpoint rooting...")
        try:
            tree.root_at_midpoint()
            if verbose:
                print("  Midpoint rooting successful")
        except Exception as e:
            if verbose:
                print(f"  Warning: Midpoint rooting failed: {e}")
        return tree
    
    elif root_method == 'outgroup':
        if not outgroup:
            raise ValueError("Outgroup sequence name required for outgroup rooting")
        
        if verbose:
            print(f"  Applying outgroup rooting with: {outgroup}")
        
        # Find the outgroup terminal
        outgroup_terminal = None
        for terminal in tree.get_terminals():
            if terminal.name == outgroup:
                outgroup_terminal = terminal
                break
        
        if not outgroup_terminal:
            if verbose:
                print(f"  Warning: Outgroup '{outgroup}' not found in tree, using midpoint rooting instead")
            tree.root_at_midpoint()
        else:
            try:
                tree.root_with_outgroup(outgroup_terminal)
                if verbose:
                    print("  Outgroup rooting successful")
            except Exception as e:
                if verbose:
                    print(f"  Warning: Outgroup rooting failed: {e}, using midpoint rooting instead")
                tree.root_at_midpoint()
        
        return tree
    
    else:
        raise ValueError(f"Unknown rooting method: {root_method}")

def save_rooted_tree(tree, original_file, root_method, output_dir=".", verbose=False):
    """Save the rooted tree to a new file for reference."""
    if not tree or root_method == 'none':
        return None
    
    # Extract gene name from original filename for shorter name
    original_basename = os.path.basename(original_file)
    
    # Try to extract gene name
    if '_verification_' in original_basename:
        gene_name = original_basename.split('_verification_')[0]
    elif '_' in original_basename:
        gene_name = original_basename.split('_')[0]
    else:
        gene_name = os.path.splitext(original_basename)[0]
    
    # Create short filename: GeneName_rootmethod.treefile
    rooted_filename = f"{gene_name}_{root_method}.treefile"
    rooted_file = os.path.join(output_dir, rooted_filename)
    
    try:
        os.makedirs(output_dir, exist_ok=True)
        Phylo.write(tree, rooted_file, "newick")
        if verbose:
            print(f"  Rooted tree saved to: {rooted_file}")
        return rooted_file
    except Exception as e:
        if verbose:
            print(f"  Warning: Could not save rooted tree: {e}")
        return None

def parse_tree_file(tree_file):
    """Parse a Newick tree file and extract all leaf names."""
    try:
        tree = Phylo.read(tree_file, "newick")
        leaves = [terminal.name for terminal in tree.get_terminals()]
        return set(leaves), tree
    except Exception as e:
        print(f"Warning: Could not parse tree with BioPython: {e}")
        return set(), None

def load_taxonomy(taxonomy_file):
    """Load taxonomy file and create mapping from sequence ID to taxonomic info."""
    df = pd.read_csv(taxonomy_file, sep='\t')
    
    taxonomy_map = {}
    for _, row in df.iterrows():
        taxonomy_map[row['full_id']] = {
            'superkingdom': row['superkingdom'],
            'phylum': row['phylum'],
            'class': row['class'] if 'class' in row else None,
            'order': row['order'] if 'order' in row else None,
            'family': row['family'] if 'family' in row else None,
            'genus': row['genus'] if 'genus' in row else None,
            'species': row['species'] if 'species' in row else None
        }
    
    return taxonomy_map

def match_tree_leaves_to_taxonomy(tree_leaves, taxonomy_map):
    """Match tree leaf names to taxonomy entries."""
    matched = {}
    unmatched = []
    
    for leaf in tree_leaves:
        if not leaf or leaf.isdigit():
            continue
            
        # Try exact match first
        if leaf in taxonomy_map:
            matched[leaf] = taxonomy_map[leaf]
        else:
            # Try to find partial matches
            found = False
            for tax_id in taxonomy_map:
                if leaf in tax_id or tax_id in leaf:
                    matched[leaf] = taxonomy_map[tax_id]
                    found = True
                    break
                
                # Try matching without GCF prefix and version numbers
                clean_tax_id = re.sub(r'^GCF_\d+\.\d+-', '', tax_id)
                clean_leaf = re.sub(r'^GCF_\d+\.\d+-', '', leaf)
                if clean_tax_id == clean_leaf:
                    matched[leaf] = taxonomy_map[tax_id]
                    found = True
                    break
                    
                # Try matching just the WP identifier
                tax_wp = re.search(r'WP_\d+\.\d+', tax_id)
                leaf_wp = re.search(r'WP_\d+\.\d+', leaf)
                if tax_wp and leaf_wp and tax_wp.group() == leaf_wp.group():
                    matched[leaf] = taxonomy_map[tax_id]
                    found = True
                    break
            
            if not found:
                unmatched.append(leaf)
    
    return matched, unmatched

def group_by_taxonomy(matched_taxonomy, level='phylum'):
    """Group sequences by taxonomic level."""
    groups = defaultdict(list)
    
    for seq_id, tax_info in matched_taxonomy.items():
        if level == 'phylum':
            superkingdom = tax_info.get('superkingdom', 'Unknown')
            phylum = tax_info.get('phylum', 'Unknown')
            
            if superkingdom == 'Archaea':
                tax_group = 'Archaea'
            elif superkingdom == 'Bacteria':
                tax_group = phylum if phylum and phylum != 'Unknown' else 'Unknown'
            else:
                tax_group = 'Unknown'
        else:
            tax_group = tax_info.get(level, 'Unknown')
        
        if tax_group and tax_group != 'Unknown':
            groups[tax_group].append(seq_id)
    
    return groups

def define_taxonomic_colors():
    """Define colors for different taxonomic groups - ordered by abundance in database."""
    colors = {
        # Tier 1: Super abundant (>10,000) - Most distinct colors
        'Pseudomonadota': '#4682b4',         # 28,618 - Steel blue
        'Bacillota': '#7df9ff',              # 13,144 - Electric blue
        'Actinomycetota': '#32cd32',         # 8,770 - Lime green
        
        # Tier 2: Very abundant (1,000-10,000) - Clear, distinguishable colors
        'Bacteroidota': '#ff1493',           # 4,688 - Deep pink
        'Cyanobacteriota': '#00bfff',        # 2,101 - Deep sky blue
        'Euryarchaeota': '#00a86b',          # 1,720 - Sea green
        'Planctomycetota': '#ffb90f',        # 1,545 - Dark goldenrod
        'Thermodesulfobacteriota': '#87ceeb', # 1,411 - Sky blue
        
        # Tier 3: Abundant (500-1,000) - Good visibility colors
        'Verrucomicrobiota': '#9370db',      # 602 - Medium orchid
        'Campylobacterota': '#6495ed',       # 530 - Cornflower blue
        
        # Tier 4: Moderate (100-500) - Distinct but can be similar families
        'Myxococcota': '#ff8c00',            # 382 - Dark orange
        'Deinococcota': '#dc143c',           # 308 - Crimson
        'Acidobacteriota': '#ff6347',        # 306 - Tomato
        'Chlamydiota': '#8a2be2',            # 293 - Blue violet
        'Chloroflexota': '#00ced1',          # 178 - Dark turquoise
        'Spirochaetota': '#9932cc',          # 164 - Dark orchid
        'Thermoproteota': '#20b2aa',         # 160 - Light sea green (new, for archaea)
        'Aquificota': '#40e0d0',             # 152 - Turquoise
        'Nitrospirota': '#98fb98',           # 127 - Pale green
        'Fusobacteriota': '#ffd700',         # 123 - Gold
        'Nitrososphaerota': '#ff69b4',       # 114 - Hot pink
        'Chlorobiota': '#dda0dd',            # 101 - Plum
                
        # Special Archaea grouping (separate from abundance-based bacterial ordering)
        'Archaea': '#00a86b',                # Fallback for other archaea
        
        # Fallback categories for rare phyla (< 26 occurrences)
        'Other Phylum': '#b0e0e6',          # Light steel blue - for all phyla not listed above
        'Unknown': '#d3d3d3'                # Light gray
    }
    return colors

def find_monophyletic_clusters(tree, group_sequences, verbose=False):
    """
    Find monophyletic clusters within a taxonomic group using proper BioPython methods.
    Returns clusters of sequences that are monophyletic, plus isolated sequences.
    """
    if not group_sequences or not tree:
        return [group_sequences] if group_sequences else []
    
    # Get terminals that belong to this group
    group_terminals = []
    group_terminal_names = set(group_sequences)
    
    for terminal in tree.get_terminals():
        if terminal.name in group_terminal_names:
            group_terminals.append(terminal)
    
    if len(group_terminals) <= 1:
        return [[terminal.name for terminal in group_terminals]]
    
    if verbose:
        print(f"    Found {len(group_terminals)} terminals for this group")
    
    # First check if entire group is monophyletic
    try:
        if hasattr(tree, 'is_monophyletic') and tree.is_monophyletic(group_terminals):
            if verbose:
                print(f"    Entire group is monophyletic!")
            return [[terminal.name for terminal in group_terminals]]
    except:
        pass
    
    # Find clusters by examining the tree structure
    clusters = []
    processed_terminals = set()
    
    # Use a recursive approach to find monophyletic clusters
    def find_clusters_recursive(clade, current_cluster):
        """Recursively traverse tree to find monophyletic clusters."""
        if clade.is_terminal():
            if clade.name in group_terminal_names:
                current_cluster.append(clade)
            return
        
        # Check children
        group_children = []
        non_group_children = []
        
        for child in clade:
            child_terminals = child.get_terminals() if not child.is_terminal() else [child]
            child_group_terminals = [t for t in child_terminals if t.name in group_terminal_names]
            
            if child_group_terminals:
                if len(child_group_terminals) == len(child_terminals):
                    # This child contains only our group sequences (monophyletic)
                    group_children.append(child)
                else:
                    # This child contains mixed sequences, need to go deeper
                    find_clusters_recursive(child, current_cluster)
            else:
                non_group_children.append(child)
        
        # If we have group-only children, they form monophyletic clusters
        for group_child in group_children:
            child_terminals = group_child.get_terminals()
            cluster_names = [t.name for t in child_terminals if t.name in group_terminal_names]
            if cluster_names:
                clusters.append(cluster_names)
                processed_terminals.update(child_terminals)
    
    try:
        find_clusters_recursive(tree.root, [])
        
        # Add any remaining individual terminals
        for terminal in group_terminals:
            if terminal not in processed_terminals:
                clusters.append([terminal.name])
        
        if verbose:
            print(f"    Found {len(clusters)} clusters with sizes: {[len(cluster) for cluster in clusters]}")
        
    except Exception as e:
        if verbose:
            print(f"    Error in cluster detection: {e}, using fallback")
        # Fallback: treat each sequence individually
        clusters = [[terminal.name] for terminal in group_terminals]
    
    return clusters if clusters else [[terminal.name for terminal in group_terminals]]

def process_taxonomic_groups(groups, tree, verbose=False):
    """Process each taxonomic group to find monophyletic clusters."""
    all_clusters = {}
    
    for group_name, sequences in groups.items():
        if verbose:
            print(f"  Processing {group_name} with {len(sequences)} sequences")
        
        # Find monophyletic clusters within this group
        clusters = find_monophyletic_clusters(tree, sequences, verbose)
        
        # Store clusters with appropriate names
        if len(clusters) == 1:
            # Single monophyletic group
            all_clusters[group_name] = clusters[0]
        else:
            # Multiple clusters (polyphyletic group)
            for i, cluster in enumerate(clusters, 1):
                cluster_name = f"{group_name}_clade{i}" if len(cluster) > 1 else f"{group_name}_isolated{i}"
                all_clusters[cluster_name] = cluster
    
    return all_clusters

def generate_itol_file(clusters, colors, output_file, tree_name, rooting_info=None):
    """Generate iTOL DATASET_RANGE file with individual sequence coloring and dynamic legend."""
    
    # Color mapping with HEX colors
    color_mapping = {}
    
    # Get unique taxonomic groups actually present in the tree
    unique_groups = set()
    for cluster_name in clusters.keys():
        base_name = cluster_name.split('_clade')[0].split('_isolated')[0]
        unique_groups.add(base_name)
    
    # Create color mapping for each cluster and track which colors are used
    used_colors = {}  # Maps display name to color
    
    for cluster_name in clusters.keys():
        base_name = cluster_name.split('_clade')[0].split('_isolated')[0]
        
        if base_name in colors:
            # Use the defined color for known phyla
            color_mapping[cluster_name] = colors[base_name]
            used_colors[base_name] = colors[base_name]
        else:
            # Use "Other Phylum" color for unknown phyla
            color_mapping[cluster_name] = colors['Other Phylum']
            used_colors['Other Phylum'] = colors['Other Phylum']
            print(f"Info: {base_name} grouped as 'Other Phylum' in legend")
    
    # Sort the used colors by the order they appear in the original colors dictionary
    # This preserves the abundance-based ordering for phyla that are actually present
    legend_items = []
    
    # First, add all the defined phyla that are actually present (in abundance order)
    for phylum_name in colors.keys():
        if phylum_name in used_colors and phylum_name != 'Other Phylum' and phylum_name != 'Unknown':
            legend_items.append((phylum_name, colors[phylum_name]))
    
    # Then add "Other Phylum" if it was used
    if 'Other Phylum' in used_colors:
        legend_items.append(('Other Phylum', colors['Other Phylum']))
    
    # Finally add "Unknown" if it was used
    if 'Unknown' in used_colors:
        legend_items.append(('Unknown', colors['Unknown']))
    
    # Extract colors and labels for legend
    legend_colors = [item[1] for item in legend_items]
    legend_labels = [item[0] for item in legend_items]
    
    with open(output_file, 'w') as f:
        f.write("DATASET_RANGE\n")
        f.write("#Automatically generated iTOL file for taxonomic range coloring\n")
        f.write(f"#Tree: {tree_name}\n")
        f.write("#Generated by generate_itol_taxonomy_colors.py (DATASET_RANGE format)\n")
        f.write("#Individual sequence coloring (each sequence colored separately)\n")
        f.write("#Legend shows only taxonomic groups present in this tree\n")
        f.write("#Rare phyla are grouped as 'Other Phylum'\n")
        
        if rooting_info:
            f.write(f"#Rooting: {rooting_info}\n")
        
        f.write("\n")
        
        # Mandatory settings
        f.write("SEPARATOR COMMA\n")
        f.write("\n")
        f.write("DATASET_LABEL,Taxonomic_Ranges\n")
        f.write("COLOR,#ffff00\n")
        f.write("\n")
        
        # Optional settings
        f.write("RANGE_TYPE,box\n")
        f.write("RANGE_COVER,label\n")  # Individual sequence coloring
        f.write("UNROOTED_SMOOTH,simplify\n")
        f.write("COVER_LABELS,1\n")     
        f.write("COVER_DATASETS,0\n")
        f.write("FIT_LABELS,0\n")
        f.write("\n")
        
        # Legend settings - only show groups actually present in tree
        if legend_colors:
            f.write("LEGEND_TITLE,Taxonomy\n")
            f.write("LEGEND_POSITION_X,100\n")
            f.write("LEGEND_POSITION_Y,100\n")
            f.write("LEGEND_HORIZONTAL,0\n")  # Vertical legend
            
            # All shapes are squares (1)
            legend_shapes = ','.join(['1'] * len(legend_colors))
            f.write(f"LEGEND_SHAPES,{legend_shapes}\n")
            
            # Colors
            legend_color_str = ','.join(legend_colors)
            f.write(f"LEGEND_COLORS,{legend_color_str}\n")
            
            # Labels
            legend_label_str = ','.join(legend_labels)
            f.write(f"LEGEND_LABELS,{legend_label_str}\n")
            
            # All shapes have scale 1
            legend_scales = ','.join(['1'] * len(legend_colors))
            f.write(f"LEGEND_SHAPE_SCALES,{legend_scales}\n")
            
        f.write("\n")
        
        # Data section
        f.write("DATA\n")
        
        # Process each cluster - COLOR INDIVIDUAL SEQUENCES, NOT CLADES
        total_sequences = 0
        for cluster_name, sequence_list in clusters.items():
            color = color_mapping.get(cluster_name, colors['Unknown'])
            
            # Color each sequence individually (seq,seq,color format)
            for seq_id in sequence_list:
                f.write(f"{seq_id},{seq_id},{color}\n")
                total_sequences += 1
        
        f.write("\n# End of range data\n")
    
    # Enhanced summary
    print(f"Generated DATASET_RANGE file with {total_sequences} individually colored sequences")
    print(f"Legend contains {len(legend_items)} taxonomic groups actually present in tree:")
    for label, color in legend_items:
        if label == 'Other Phylum':
            # Count how many different phyla are grouped under "Other Phylum"
            other_phyla = []
            for cluster_name in clusters.keys():
                base_name = cluster_name.split('_clade')[0].split('_isolated')[0]
                if base_name not in colors or base_name in ['Other Phylum', 'Unknown']:
                    if base_name not in ['Other Phylum', 'Unknown'] and base_name not in other_phyla:
                        other_phyla.append(base_name)
            if other_phyla:
                print(f"  - {label}: {color} (includes: {', '.join(other_phyla)})")
            else:
                print(f"  - {label}: {color}")
        else:
            count = sum(1 for cluster_name in clusters.keys() 
                       if cluster_name.split('_clade')[0].split('_isolated')[0] == label)
            print(f"  - {label}: {color} ({count} group{'s' if count != 1 else ''})")

def generate_text_labels_file(clusters, colors, output_file, tree_name, rooting_info=None):
    """Generate iTOL DATASET_TEXT file for external taxonomic labels with real phylum names."""
    
    # Color mapping - use actual colors but show real taxonomic names
    color_mapping = {}
    for cluster_name in clusters.keys():
        base_name = cluster_name.split('_clade')[0].split('_isolated')[0]
        
        if base_name in colors:
            # Use the defined color for known phyla
            color_mapping[cluster_name] = colors[base_name]
        else:
            # Use "Other Phylum" color but keep the real taxonomic name in the label
            color_mapping[cluster_name] = colors['Other Phylum']
    
    with open(output_file, 'w') as f:
        f.write("DATASET_TEXT\n")
        f.write("#External text labels for taxonomic groups\n")
        f.write("#Shows actual phylum names only (no cluster/isolated info)\n")
        f.write(f"#Tree: {tree_name}\n")
        f.write("#Generated by generate_itol_taxonomy_colors.py\n")
        if rooting_info:
            f.write(f"#Rooting: {rooting_info}\n")
        f.write("\n")
        
        f.write("#=================================================================#\n")
        f.write("#                    MANDATORY SETTINGS                           #\n")
        f.write("#=================================================================#\n")
        f.write("SEPARATOR TAB\n\n")
        f.write("DATASET_LABEL\tTaxonomic Groups\n")
        f.write("COLOR\t#000000\n\n")
        
        f.write("#=================================================================#\n")
        f.write("#                    OPTIONAL SETTINGS                           #\n")
        f.write("#=================================================================#\n")
        f.write("MARGIN\t20\n")
        f.write("SHOW_INTERNAL\t0\n")
        f.write("LABEL_ROTATION\t0\n")
        f.write("ALIGN_TO_TREE\t0\n")
        f.write("SIZE_FACTOR\t1.3\n\n")
        
        f.write("#=================================================================#\n")
        f.write("#       Actual data follows after the \"DATA\" keyword              #\n")
        f.write("#=================================================================#\n")
        f.write("DATA\n")
        f.write("#ID\tlabel\tposition\tcolor\tstyle\tsize_factor\trotation\n\n")
        
        # Add labels for each cluster with SIMPLE taxonomic names only
        for cluster_name, sequence_list in clusters.items():
            color = color_mapping.get(cluster_name, colors['Unknown'])
            representative = sequence_list[0]  # Use first sequence as representative
            
            # Extract the real taxonomic name (base_name) from cluster_name
            base_name = cluster_name.split('_clade')[0].split('_isolated')[0]
            
            # Simple display label - just the phylum name
            display_label = base_name
            
            # Use a different text style for "Other Phylum" colored groups to distinguish them
            if base_name in colors:
                # Major phyla get bold style
                style = "bold"
            else:
                # Minor phyla (colored as "Other Phylum") get italic style to show they're grouped
                style = "italic"
            
            f.write(f"{representative}\t{display_label}\t-1\t{color}\t{style}\t1.3\t0\n")

def main():
    parser = argparse.ArgumentParser(description='Generate iTOL taxonomy coloring files with improved monophyly detection')
    parser.add_argument('tree_file', help='Path to Newick tree file')
    parser.add_argument('taxonomy_file', help='Path to taxonomy TSV file')
    parser.add_argument('-o', '--output', help='Output iTOL file path (default: auto-generated)')
    parser.add_argument('-l', '--level', default='phylum', 
                       choices=['superkingdom', 'phylum', 'class', 'order', 'family'],
                       help='Taxonomic level for grouping (default: phylum)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Verbose output')
    parser.add_argument('--labels', action='store_true', 
                       help='Generate external text labels file')
    parser.add_argument('--root-method', default='midpoint',
                       choices=['midpoint', 'outgroup', 'none'],
                       help='Tree rooting method (default: midpoint)')
    parser.add_argument('--outgroup', 
                       help='Outgroup sequence name (required for --root-method outgroup)')
    parser.add_argument('--save-rooted', action='store_true',
                       help='Save the rooted tree to a new file for reference')
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.root_method == 'outgroup' and not args.outgroup:
        parser.error("--outgroup is required when using --root-method outgroup")
    
    # Check files exist
    if not os.path.exists(args.tree_file):
        print(f"Error: Tree file {args.tree_file} not found")
        sys.exit(1)
    
    if not os.path.exists(args.taxonomy_file):
        print(f"Error: Taxonomy file {args.taxonomy_file} not found")
        sys.exit(1)
    
    # Generate output filename if not provided
    if not args.output:
        tree_basename = os.path.splitext(os.path.basename(args.tree_file))[0]
        args.output = f"{tree_basename}_itol_colors.txt"
    
    labels_output = args.output.replace('_colors.txt', '_labels.txt')
    
    if args.verbose:
        print(f"Processing tree: {args.tree_file}")
        print(f"Using taxonomy: {args.taxonomy_file}")
        print(f"Grouping by: {args.level}")
        print(f"Rooting method: {args.root_method}")
        if args.outgroup:
            print(f"Outgroup: {args.outgroup}")
        print(f"Output file: {args.output}")
        if args.labels:
            print(f"Labels file: {labels_output}")
        if args.save_rooted:
            print(f"Save rooted tree: Enabled")
    
    # Load data
    tree_leaves, tree = parse_tree_file(args.tree_file)
    taxonomy_map = load_taxonomy(args.taxonomy_file)
    
    if args.verbose:
        print(f"Found {len(tree_leaves)} leaves in tree")
        print(f"Loaded {len(taxonomy_map)} taxonomy entries")
    
    # Apply rooting
    rooting_info = None
    if tree:
        print(f"\nApplying tree rooting...")
        tree = apply_tree_rooting(tree, args.root_method, args.outgroup, args.verbose)
        
        if args.root_method == 'midpoint':
            rooting_info = "Midpoint rooting applied"
        elif args.root_method == 'outgroup':
            rooting_info = f"Outgroup rooting with {args.outgroup}"
        else:
            rooting_info = "No rooting applied"
        
        # Save rooted tree if requested
        if args.save_rooted:
            save_rooted_tree(tree, args.tree_file, args.root_method, ".", args.verbose)
    
    # Match sequences to taxonomy
    matched_taxonomy, unmatched = match_tree_leaves_to_taxonomy(tree_leaves, taxonomy_map)
    
    if args.verbose:
        print(f"Matched {len(matched_taxonomy)} sequences to taxonomy")
        if unmatched:
            print(f"Warning: {len(unmatched)} sequences could not be matched")
    
    # Group by taxonomy
    groups = group_by_taxonomy(matched_taxonomy, args.level)
    
    if args.verbose:
        print(f"\nTaxonomic groups found ({args.level}):")
        for group_name, sequences in groups.items():
            print(f"  {group_name}: {len(sequences)} sequences")
    
    # Find monophyletic clusters
    print(f"\nIdentifying monophyletic clusters...")
    clusters = process_taxonomic_groups(groups, tree, args.verbose)
    
    # Define colors
    colors = define_taxonomic_colors()
    
    # Generate iTOL file
    tree_name = os.path.basename(args.tree_file)
    generate_itol_file(clusters, colors, args.output, tree_name, rooting_info)
    
    print(f"iTOL file generated: {args.output}")
    
    # Generate text labels file if requested
    if args.labels:
        generate_text_labels_file(clusters, colors, labels_output, tree_name, rooting_info)
        print(f"iTOL labels file generated: {labels_output}")
    
    if args.verbose:
        print("\nSummary:")
        print(f"  Total sequences in tree: {len(tree_leaves)}")
        print(f"  Sequences matched to taxonomy: {len(matched_taxonomy)}")
        print(f"  Taxonomic groups: {len(groups)}")
        print(f"  Monophyletic clusters found: {len(clusters)}")
        print(f"  Tree rooting: {rooting_info}")
        
        # Report cluster information
        monophyletic_clusters = [name for name, seqs in clusters.items() if len(seqs) > 1 and '_isolated' not in name]
        isolated_sequences = [name for name in clusters.keys() if '_isolated' in name]
        
        if monophyletic_clusters:
            print(f"  Monophyletic clusters: {len(monophyletic_clusters)}")
        if isolated_sequences:
            print(f"  Isolated sequences: {len(isolated_sequences)}")

if __name__ == "__main__":
    main()