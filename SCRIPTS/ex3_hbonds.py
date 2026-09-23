#!/usr/bin/env python3
"""
Exercise 3
----------
Determines all possible hydrogen bonds: pairs of polar atoms (O, N, S)
closer than a cutoff distance (default 3.5 A).

Parameters: PDB file (or --pdb_id). Optional: cutoff distance (default
3.5).

Usage:
    python ex3_hbonds.py --pdb_file structure.pdb
    python ex3_hbonds.py --pdb_file structure.pdb --cutoff 3.2

The results are printed to the console and also saved automatically to a
text file (see --output).
"""
import sys

from Bio.PDB import NeighborSearch

from common import (
    POLAR_ELEMENTS,
    add_structure_source_args,
    atom_id,
    base_parser,
    emit_report,
    load_structure,
)


def get_polar_atoms(structure, model_id=0):
    return [
        atom
        for atom in structure[model_id].get_atoms()
        if atom.element in POLAR_ELEMENTS
    ]


def find_hbonds(structure, cutoff=3.5, model_id=0):
    polar_atoms = get_polar_atoms(structure, model_id)
    if len(polar_atoms) < 2:
        return []
    ns = NeighborSearch(polar_atoms)
    pairs = ns.search_all(cutoff, level="A")

    results = []
    for atom1, atom2 in pairs:
        res1 = atom1.get_parent()
        res2 = atom2.get_parent()
        if res1 is res2:
            # skip pairs within the same residue (not hydrogen bonds)
            continue
        d = atom1 - atom2
        results.append((atom1, atom2, d))

    def sort_key(item):
        a1, a2, _ = item
        r1, r2 = a1.get_parent(), a2.get_parent()
        return ((r1.get_parent().id, r1.id[1], a1.get_name()),
                (r2.get_parent().id, r2.id[1], a2.get_name()))

    results.sort(key=sort_key)
    return results


def main():
    parser = base_parser(
        "Exercise 3: possible hydrogen bonds (polar atoms O/N/S closer "
        "than a cutoff distance)"
    )
    add_structure_source_args(parser)
    parser.add_argument(
        "--cutoff", type=float, default=3.5,
        help="Cutoff distance in Angstrom (default 3.5)",
    )
    args = parser.parse_args()

    structure = load_structure(args)
    hbonds = find_hbonds(structure, args.cutoff)

    lines = [
        f"# Possible hydrogen bonds (< {args.cutoff} A)",
        f"# {'Atom 1':<20}{'Atom 2':<20}{'Distance (A)'}",
    ]
    for atom1, atom2, d in hbonds:
        lines.append(f"{atom_id(atom1):<20}{atom_id(atom2):<20}{d:.2f}")
    lines.append(f"# Total: {len(hbonds)} pairs")

    emit_report(lines, args, "ex3_hbonds")


if __name__ == "__main__":
    sys.exit(main())
