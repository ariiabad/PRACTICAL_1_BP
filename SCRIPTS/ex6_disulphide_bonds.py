#!/usr/bin/env python3
"""
Exercise 6
----------
Same idea as exercise 5, but for disulphide bonds: they form between the
SG atoms of two CYS residues when they are at the appropriate distance.

Parameters: PDB file (or --pdb_id). Optional: cutoff distance (default
2.5 A).

Usage:
    python ex6_disulphide_bonds.py --pdb_file structure.pdb
    python ex6_disulphide_bonds.py --pdb_file structure.pdb --cutoff 2.2

"""
import sys

from Bio.PDB import NeighborSearch

from common import add_structure_source_args, base_parser, emit_report, load_structure, residue_id


def find_disulphide_bonds(structure, cutoff=2.5, model_id=0):
    sg_atoms = [
        residue["SG"]
        for residue in structure[model_id].get_residues()
        if residue.resname == "CYS" and "SG" in residue
    ]

    if len(sg_atoms) < 2:
        # no CYS, or a single one: there cannot be a disulphide bond
        return []

    ns = NeighborSearch(sg_atoms)
    pairs = ns.search_all(cutoff, level="A")

    results = []
    for atom1, atom2 in pairs:
        res1, res2 = atom1.get_parent(), atom2.get_parent()
        if res1 is res2:
            continue
        d = atom1 - atom2
        results.append((res1, res2, d))

    def sort_key(item):
        r1, r2, _ = item
        return ((r1.get_parent().id, r1.id[1]), (r2.get_parent().id, r2.id[1]))

    results.sort(key=sort_key)
    return results


def main():
    parser = base_parser(
        "Exercise 6: disulphide bonds (SG-SG contacts between CYS residues)"
    )
    add_structure_source_args(parser)
    parser.add_argument(
        "--cutoff", type=float, default=2.5,
        help="Cutoff distance in Angstrom (default 2.5; the typical S-S "
        "bond is ~2.05 A)",
    )
    args = parser.parse_args()

    structure = load_structure(args)
    bonds = find_disulphide_bonds(structure, args.cutoff)

    lines = [
        f"# Disulphide bonds SG-SG (< {args.cutoff} A)",
        f"# {'Residue 1':<15}{'Residue 2':<15}{'Distance (A)'}",
    ]
    for res1, res2, d in bonds:
        lines.append(f"{residue_id(res1):<15}{residue_id(res2):<15}{d:.2f}")
    lines.append(f"# Total: {len(bonds)} disulphide bonds")

    emit_report(lines, args, "ex6_disulphide_bonds")


if __name__ == "__main__":
    sys.exit(main())
