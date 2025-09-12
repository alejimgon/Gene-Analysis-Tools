# Taxonomy Retrieval Scripts

Scripts for retrieving taxonomy information using NCBI Taxonomy.

## Usage

```bash
python taxonomy_ncbitaxa.py path/to/taxids.tsv path/to/blastp.txt path/to/output_file.tsv
```

- **Taxonomy File:** TSV with `assembly_accession` and `taxid` columns.
- **BLASTP Results:** Text file with lines like `assembly_accession-unique_id`.

## Output

- Taxonomy results TSV
- Error log TSV
- Log file