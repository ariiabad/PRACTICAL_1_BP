#!/usr/bin/env python3
"""
Exercise 5
----------
Generates the backbone connectivity list, i.e. which residues are linked
by an ordinary peptide bond (atom C of one residue with atom N of the
following residue).

Parameters: PDB file (or --pdb_id). Optional: cutoff distance for the
peptide bond (default 2.5 A).

Usage:
    python ex5_backbone_connectivity.py --pdb_file structure.pdb
    python ex5_backbone_connectivity.py --pdb_file structure.pdb --cutoff 2.0

The results are printed to the console and also saved automatically to a
text file (see --output).
"""
import sys

from Bio.PDB import NeighborSearch

from common import add_structure_source_args, base_parser, emit_report, load_structure, residue_id


def find_peptide_bonds(structure, cutoff=2.5, model_id=0):
    c_atoms = [
        residue["C"] for residue in structure[model_id].get_residues() if "C" in residue
    ]
    n_atoms = [
        residue["N"] for residue in structure[model_id].get_residues() if "N" in residue
    ]

    if not n_atoms or not c_atoms:
        return []

    ns = NeighborSearch(n_atoms)
    results = []
    for c_atom in c_atoms:
        res_c = c_atom.get_parent()
        neighbours = ns.search(c_atom.get_coord(), cutoff, level="A")
        for n_atom in neighbours:
            res_n = n_atom.get_parent()
            if res_c is res_n:
                # exclude the N of the same residue (not a peptide bond)
                continue
            d = c_atom - n_atom
            results.append((res_c, res_n, d))

    def sort_key(item):
        r1, r2, _ = item
        return ((r1.get_parent().id, r1.id[1]), (r2.get_parent().id, r2.id[1]))

    results.sort(key=sort_key)
    return results


def main():
    parser = base_parser(
        "Exercise 5: backbone connectivity (peptide bonds C-N between "
        "residues)"
    )
    add_structure_source_args(parser)
    parser.add_argument(
        "--cutoff", type=float, default=2.5,
        help="Cutoff distance in Angstrom for the peptide bond (default 2.5)",
    )
    args = parser.parse_args()

    structure = load_structure(args)
    bonds = find_peptide_bonds(structure, args.cutoff)

    lines = [
        f"# Peptide bonds C-N (< {args.cutoff} A)",
        f"# {'Residue (C)':<15}{'Residue (N)':<15}{'Distance (A)'}",
    ]
    for res_c, res_n, d in bonds:
        lines.append(f"{residue_id(res_c):<15}{residue_id(res_n):<15}{d:.2f}")
    lines.append(f"# Total: {len(bonds)} bonds")

    emit_report(lines, args, "ex5_backbone_connectivity")


if __name__ == "__main__":
    sys.exit(main())
