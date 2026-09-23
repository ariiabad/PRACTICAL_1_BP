#!/usr/bin/env python3
"""
Exercise 2
----------
Generates the list of all atoms of a given residue (with their
coordinates).

Parameters: PDB file (or --pdb_id), residue number 

Usage:
    python ex2_residue_atoms.py --pdb_file structure.pdb --residue A:315
    

The results are printed to the console and also saved automatically to a
text file (see --output).
"""
import sys

from common import (
    add_structure_source_args,
    base_parser,
    emit_report,
    find_residue,
    load_structure,
    residue_id,
)


def list_residue_atoms(residue):
    return sorted(residue.get_atoms(), key=lambda a: a.get_serial_number())


def main():
    parser = base_parser(
        "Exercise 2: list of atoms (and coordinates) of a given residue"
    )
    add_structure_source_args(parser)
    parser.add_argument(
        "--residue", required=True,
        help="Residue number, optionally with chain: '315' or 'A:315'",
    )
    args = parser.parse_args()

    structure = load_structure(args)
    residue = find_residue(structure, args.residue)
    atoms = list_residue_atoms(residue)

    lines = [
        f"# Residue {residue_id(residue)}",
        f"# {'Atom':<10}{'X':>10}{'Y':>10}{'Z':>10}",
    ]
    for atom in atoms:
        x, y, z = atom.get_coord()
        lines.append(f"{atom.get_name():<10}{x:>10.3f}{y:>10.3f}{z:>10.3f}")
    lines.append(f"# Total: {len(atoms)} atoms")

    emit_report(lines, args, "ex2_residue_atoms")


if __name__ == "__main__":
    sys.exit(main())
