#!/usr/bin/env python3
"""
Exercise 7
----------
Prints the distances between all pairs of atoms of two given residues.

Parameters: PDB file (or --pdb_id), residue 1, residue 2 (format NUMBER or
CHAIN:NUMBER, same as in exercise 2).

Usage:
    python ex7_two_residue_distances.py --pdb_file structure.pdb \
        --residue1 A:2 --residue2 A:7

The results are printed to the console and also saved automatically to a
text file (see --output).
"""
import sys

from common import (
    add_structure_source_args,
    atom_id,
    base_parser,
    emit_report,
    find_residue,
    load_structure,
    residue_id,
)


def all_atom_distances(residue1, residue2):
    results = []
    for atom1 in residue1.get_atoms():
        for atom2 in residue2.get_atoms():
            d = atom1 - atom2
            results.append((atom1, atom2, d))
    results.sort(key=lambda item: (item[0].get_serial_number(), item[1].get_serial_number()))
    return results


def main():
    parser = base_parser(
        "Exercise 7: distances between all pairs of atoms of two residues"
    )
    add_structure_source_args(parser)
    parser.add_argument("--residue1", required=True, help="First residue (NUM or CHAIN:NUM)")
    parser.add_argument("--residue2", required=True, help="Second residue (NUM or CHAIN:NUM)")
    args = parser.parse_args()

    structure = load_structure(args)
    residue1 = find_residue(structure, args.residue1)
    residue2 = find_residue(structure, args.residue2)

    distances = all_atom_distances(residue1, residue2)

    lines = [
        f"# Distances between {residue_id(residue1)} and {residue_id(residue2)}",
        f"# {'Atom 1':<18}{'Atom 2':<18}{'Distance (A)'}",
    ]
    for atom1, atom2, d in distances:
        lines.append(f"{atom_id(atom1):<18}{atom_id(atom2):<18}{d:.2f}")
    lines.append(f"# Total: {len(distances)} atom pairs")

    emit_report(lines, args, "ex7_two_residue_distances")


if __name__ == "__main__":
    sys.exit(main())
