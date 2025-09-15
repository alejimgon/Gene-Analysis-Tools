# Sequence Extraction Scripts

Perl scripts for extracting sequences from a FASTA file based on a list of sequence IDs.

## Scripts

- **extract_seq.pl**  
  Prints sequences present in a reference list of IDs.

- **extract_seq_not_equal.pl**  
  Writes sequences **not** present in a reference list of IDs to an output file.

## Usage

### Extract sequences present in a list

```bash
perl extract_seq.pl path/to/input.fasta path/to/id_list.txt > output.fasta
```

- **input.fasta:** Reference FASTA file.
- **id_list.txt:** Text file with one sequence ID per line (if using a FASTA, extract IDs with `grep '>' fasta_file | sed 's/>//'`).
- **output.fasta:** Output FASTA file (redirected from STDOUT).

### Extract sequences **not** present in a list

```bash
perl extract_seq_not_equal.pl path/to/input.fasta path/to/id_list.txt path/to/output.fasta
```

- **input.fasta:** Reference FASTA file.
- **id_list.txt:** Text file with one sequence ID per line.
- **output.fasta:** Output FASTA file.

## Output

- FASTA file containing the extracted sequences.

## Requirements

- Perl (tested with Perl 5)