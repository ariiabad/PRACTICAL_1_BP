# PRACTICAL_1_BP
Practical 1 — Introduction to protein structure manipulation in Python

##  Create the environment and install Biopython

In the terminal execute:

```bash
conda create -n biophysics python=3.10
conda activate biophysics
conda install -c conda-forge biopython
```

## How to run the scripts


Every script:

- accepts a local file (`--pdb_file`) **or** a PDB code to download (`--pdb_id`);
- prints its results to the console **and** saves them automatically to a
  `.txt` file (`--output` to choose the path);
- handles edge cases cleanly (missing chains, ambiguous residues, structures
  with no CYS, failed downloads) instead of crashing with a traceback.


## Usage

Use `--help` on any script to see all its parameters. A few examples using
the two structures:
```bash
# Exercise 1: pairs of residues with CA closer than 8 A
python ex1_ca_distances.py --pdb_file examples/real_pdbs/1UBQ.pdb --distance 8.0

# Exercise 2: atoms of a residue (NUM or CHAIN:NUM)
python ex2_residue_atoms.py --pdb_file examples/real_pdbs/4HHB.pdb --residue B:92

# Exercise 3: possible hydrogen bonds (3.5 A cutoff by default)
python ex3_hbonds.py --pdb_file examples/real_pdbs/1UBQ.pdb --cutoff 3.2

# Exercise 4: CA of all residues of a given type (1 or 3 letters)
python ex4_ca_by_restype.py --pdb_file examples/real_pdbs/4HHB.pdb --restype ARG

# Exercise 5: backbone connectivity (peptide bonds)
python ex5_backbone_connectivity.py --pdb_file examples/real_pdbs/1UBQ.pdb

# Exercise 6: disulphide bonds
python ex6_disulphide_bonds.py --pdb_file examples/real_pdbs/4HHB.pdb

# Exercise 7: distances between all atoms of two residues
python ex7_two_residue_distances.py --pdb_file examples/real_pdbs/1UBQ.pdb --residue1 A:1 --residue2 A:63

# Exercise 8: use a PDB id instead of a local file, in any script
python ex1_ca_distances.py --pdb_id 1AKI --distance 8.0
```

Results print to the console and are also written to a `.txt` file named
`<exercise>_<structure_id>.txt` in the current directory (e.g.
`ex1_ca_distances_1UBQ.txt`), unless `--output` is given.

#

## Example structures

The pdbs we are working with are two real structures from the PDB, used to
generate the sample output in `examples/outputs/`:

- **4HHB** — deoxygenated human hemoglobin, 4 chains (A/C: alpha, 141
  residues; B/D: beta, 146 residues).
- **1UBQ** — ubiquitin, 1 chain, 76 residues, with no CYS residues.

## Notes

- `common.py` centralizes structure loading (local file or `--pdb_id`),
  readable residue/atom formatting, residue-argument parsing
  (`CHAIN:NUMBER` or `NUMBER`), the 1-/3-letter amino acid code tables, and
  `emit_report`, which prints and saves results in one step.

