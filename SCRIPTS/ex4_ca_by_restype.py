#!/usr/bin/env python3
"""
Exercise 4
----------
Generates the list of all CA atoms of a given residue type, with their
coordinates. Accepts the residue code in 1- or 3-letter format
automatically 

Parameters: PDB file (or --pdb_id), residue type.

Usage:
    python ex4_ca_by_restype.py --pdb_file structure.pdb --restype ARG
    python ex4_ca_by_restype.py --pdb_file structure.pdb --restype R

"""
import sys

from common import (
    add_structure_source_args,
    base_parser,
    emit_report,
    load_structure,
    normalize_resname,
    residue_id,
)


def find_ca_by_restype(structure, resname, model_id=0):
    residues = [
        residue
        for residue in structure[model_id].get_residues()
        if residue.resname == resname and "CA" in residue
    ]
    residues.sort(key=lambda r: (r.get_parent().id, r.id[1]))
    return residues


def main():
    parser = base_parser(
        "Exercise 4: CA atoms (with coordinates) of all residues of a "
        "given type"
    )
    add_structure_source_args(parser)
    parser.add_argument(
        "--restype", required=True,
        help="Residue type, either 3-letter (ARG) or 1-letter (R) format",
    )
    args = parser.parse_args()

    try:
        resname = normalize_resname(args.restype)
    except ValueError as e:
        sys.exit(f"Error: {e}")

    structure = load_structure(args)
    residues = find_ca_by_restype(structure, resname)

    lines = [
        f"# CA atoms of {resname} ({args.restype}) residues",
        f"# {'Residue':<15}{'X':>10}{'Y':>10}{'Z':>10}",
    ]
    for residue in residues:
        x, y, z = residue["CA"].get_coord()
        lines.append(f"{residue_id(residue):<15}{x:>10.3f}{y:>10.3f}{z:>10.3f}")
    lines.append(f"# Total: {len(residues)} residues")

    emit_report(lines, args, "ex4_ca_by_restype")


if __name__ == "__main__":
    sys.exit(main())
