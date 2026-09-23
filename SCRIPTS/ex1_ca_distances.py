#!/usr/bin/env python3
"""
Exercise 1
----------
Determines the list of pairs of residues whose CA atoms are closer than a
given distance.

Parameters: PDB file (or --pdb_id), distance.

Usage:
    python ex1_ca_distances.py --pdb_file structure.pdb --distance 8.0
    python ex1_ca_distances.py --pdb_id 1AKI --distance 8.0
"""
import sys

from Bio.PDB import NeighborSearch

from common import (
    add_structure_source_args,
    base_parser,
    emit_report,
    load_structure,
    residue_id,
)


def get_ca_atoms(structure, model_id=0):
    return [
        residue["CA"]
        for residue in structure[model_id].get_residues()
        if "CA" in residue
    ]


def find_close_residue_pairs(structure, distance, model_id=0):
    ca_atoms = get_ca_atoms(structure, model_id)
    if len(ca_atoms) < 2:
        return []
    ns = NeighborSearch(ca_atoms)
    pairs = ns.search_all(distance, level="A")

    results = []
    for atom1, atom2 in pairs:
        res1 = atom1.get_parent()
        res2 = atom2.get_parent()
        if res1 is res2:
            continue
        d = atom1 - atom2
        results.append((res1, res2, d))

    # stable order by (chain, number) of the first residue, then the second
    def sort_key(item):
        r1, r2, _ = item
        k1 = (r1.get_parent().id, r1.id[1])
        k2 = (r2.get_parent().id, r2.id[1])
        return (k1, k2)

    results.sort(key=sort_key)
    return results


def main():
    parser = base_parser(
        "Exercise 1: pairs of residues whose CA atoms are closer than a "
        "given distance"
    )
    add_structure_source_args(parser)
    parser.add_argument(
        "--distance", type=float, required=True,
        help="Distance threshold in Angstrom",
    )
    args = parser.parse_args()

    structure = load_structure(args)
    pairs = find_close_residue_pairs(structure, args.distance)

    lines = [
        f"# Residue pairs with CA-CA < {args.distance} A",
        f"# {'Residue 1':<15}{'Residue 2':<15}{'Distance (A)'}",
    ]
    for res1, res2, d in pairs:
        lines.append(f"{residue_id(res1):<15}{residue_id(res2):<15}{d:.2f}")
    lines.append(f"# Total: {len(pairs)} pairs")

    emit_report(lines, args, "ex1_ca_distances")


if __name__ == "__main__":
    sys.exit(main())
