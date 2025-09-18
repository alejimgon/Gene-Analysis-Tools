#!/usr/bin/env python3
"""
Extracts a region from GenBank files based on a table of gene pairs and accession numbers.
For each line in the input TSV:
- Finds the two reference genes in the GenBank file.
- Extracts all features between them (inclusive), plus two genes upstream and two downstream (if available).
- Writes the region to a new GenBank file.
"""

import sys
import csv
from pathlib import Path
from Bio import SeqIO

# Paths (edit if needed)
BASE_DIR = Path(__file__).resolve().parent.parent
TSV_PATH = BASE_DIR / "path/to/gene_pairs.tsv"
GENBANK_DIR = BASE_DIR / "path/to/genbank_files"
OUTPUT_DIR = BASE_DIR / "path/to/short_genbanks"
OUTPUT_DIR.mkdir(exist_ok=True)

def get_gene_feature_indices(features, gene_ids):
    idxs = []
    for i, feat in enumerate(features):
        if feat.type == "CDS":
            # Try both 'locus_tag' and 'protein_id' and 'gene'
            for key in ("locus_tag", "protein_id", "gene"):
                if key in feat.qualifiers and feat.qualifiers[key][0] in gene_ids:
                    idxs.append(i)
    return idxs

def extract_region(features, idx1, idx2, flank=2):
    start = max(0, min(idx1, idx2) - flank)
    end = min(len(features), max(idx1, idx2) + flank + 1)
    return features[start:end]

def main():
    from Bio.SeqRecord import SeqRecord
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    with open(TSV_PATH) as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            gene1, gene2, acc = row
            gbk_path = GENBANK_DIR / f"{acc}.gbff"
            if not gbk_path.exists():
                print(f"GenBank file not found: {gbk_path}", file=sys.stderr)
                continue
            for record in SeqIO.parse(str(gbk_path), "genbank"):
                cds_features = [f for f in record.features if f.type == "CDS"]
                idxs = get_gene_feature_indices(cds_features, {gene1, gene2})
                if len(idxs) < 2:
                    print(f"Could not find both genes {gene1}, {gene2} in {acc}", file=sys.stderr)
                    continue
                region_features = extract_region(cds_features, idxs[0], idxs[1], flank=2)
                # Find min and max coordinates for the region
                region_starts = [int(f.location.start) for f in region_features]
                region_ends = [int(f.location.end) for f in region_features]
                region_start = min(region_starts)
                region_end = max(region_ends)
                # Extract subsequence
                sub_seq = record.seq[region_start:region_end]
                # Adjust features to new coordinates
                new_features = []
                for f in region_features:
                    # Shift feature location
                    loc = f.location
                    new_loc = FeatureLocation(
                        int(loc.start) - region_start,
                        int(loc.end) - region_start,
                        strand=loc.strand
                    )
                    new_feat = SeqFeature(
                        location=new_loc,
                        type=f.type,
                        qualifiers=f.qualifiers
                    )
                    new_features.append(new_feat)
                # Copy record and update
                new_record = SeqRecord(
                    sub_seq,
                    id=record.id,
                    name=record.name,
                    description=record.description,
                    dbxrefs=record.dbxrefs,
                    annotations=record.annotations
                )
                new_record.features = new_features
                # Write to file
                out_name = f"{acc}_short_{gene1}_{gene2}.gbk"
                out_path = OUTPUT_DIR / out_name
                SeqIO.write(new_record, str(out_path), "genbank")
                print(f"Wrote {out_path}")

if __name__ == "__main__":
    main()
