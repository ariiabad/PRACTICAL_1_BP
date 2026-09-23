"""
common.py
---------
Shared functions used by all exercise scripts (ex1 ... ex7).
Not meant to be run directly.
"""
import argparse
import os
import re
import sys

from Bio.PDB import PDBParser

# 3-letter <-> 1-letter amino acid codes (defined by hand, the 20 standard
# residues), used to accept both formats in exercise 4.
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
ONE_TO_THREE = {one: three for three, one in THREE_TO_ONE.items()}

# Atoms considered "polar" for exercise 3 (possible hydrogen bonds)
POLAR_ELEMENTS = {"O", "N", "S"}


def normalize_resname(user_input):
    """Accepts 'ARG' or 'R' (any case) and returns the 3-letter code
    ('ARG'). Raises ValueError if not a standard amino acid."""
    code = user_input.strip().upper()
    if len(code) == 3 and code in THREE_TO_ONE:
        return code
    if len(code) == 1 and code in ONE_TO_THREE:
        return ONE_TO_THREE[code]
    raise ValueError(
        f"Residue '{user_input}' not recognized. Use the 3-letter code "
        f"(e.g. ARG) or the 1-letter code (e.g. R)."
    )


def add_structure_source_args(parser):
    """Adds the structure source arguments (a local PDB file, or a PDB id
    to download automatically — exercise 8) plus the output file
    argument, to an ArgumentParser."""
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--pdb_file", help="Path to a local PDB file")
    group.add_argument(
        "--pdb_id",
        help="PDB code (e.g. 1AKI) to download automatically instead of "
        "using a local file (exercise 8, requires internet)",
    )
    parser.add_argument(
        "--download_dir",
        default="pdb_cache",
        help="Folder where files downloaded with --pdb_id are saved "
        "(default: ./pdb_cache)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Path of the output text file. If omitted, a name is "
        "generated automatically (e.g. ex1_ca_distances_1UBQ.txt)",
    )
    return parser


def download_structure_file(pdb_id, download_dir="pdb_cache"):
    """Downloads a PDB by its id and returns the path of the downloaded
    file. Exercise 8."""
    from Bio.PDB import PDBList

    os.makedirs(download_dir, exist_ok=True)
    pdbl = PDBList(verbose=False)
    fname = pdbl.retrieve_pdb_file(pdb_id, pdir=download_dir, file_format="pdb")
    if not fname or not os.path.isfile(fname):
        raise FileNotFoundError(
            f"Could not download structure '{pdb_id}'. Check the PDB code "
            f"and your internet connection."
        )
    return fname


def structure_id_from_args(args):
    """Short id for the structure (the PDB id, or the local file's
    basename without extension), used to build the default output file
    name."""
    if args.pdb_id:
        return args.pdb_id.upper()
    return os.path.splitext(os.path.basename(args.pdb_file))[0]


def load_structure(args):
    """Loads and returns the Biopython Structure object, either from
    args.pdb_file or by downloading args.pdb_id (exercise 8)."""
    if args.pdb_id:
        try:
            pdb_path = download_structure_file(args.pdb_id, args.download_dir)
        except FileNotFoundError as e:
            sys.exit(f"Error: {e}")
    else:
        pdb_path = args.pdb_file
        if not os.path.isfile(pdb_path):
            sys.exit(f"Error: PDB file '{pdb_path}' not found")

    parser = PDBParser(QUIET=True)
    return parser.get_structure(structure_id_from_args(args), pdb_path)


def residue_id(residue):
    """Readable string such as 'A:ARG315' (chain:name+number)."""
    chain = residue.get_parent()
    chain_id = chain.id if chain is not None else "?"
    _, resnum, icode = residue.id
    return f"{chain_id}:{residue.resname}{resnum}{icode.strip()}"


def atom_id(atom):
    """Readable string such as 'A:ARG315.CA'."""
    return f"{residue_id(atom.get_parent())}.{atom.get_name()}"


# Parses a residue argument such as "A:315", "A:315A" or "315" (chain only
# needed if the residue number is ambiguous across chains).
_RESIDUE_RE = re.compile(r"^(?:([A-Za-z0-9]):)?(-?\d+)([A-Za-z]?)$")


def parse_residue_spec(spec):
    """'A:315', '315A' or '315' -> (chain_id_or_None, resnum, icode)."""
    m = _RESIDUE_RE.match(spec.strip())
    if not m:
        raise ValueError(
            f"Invalid residue format: '{spec}'. Use NUMBER or "
            f"CHAIN:NUMBER (e.g. 315 or A:315)."
        )
    chain_id, resnum, icode = m.groups()
    return chain_id, int(resnum), (icode if icode else " ")


def find_residue(structure, spec, model_id=0):
    """Looks up a residue from its text specification. If no chain is
    given and the number is ambiguous across chains, asks for one."""
    chain_id, resnum, icode = parse_residue_spec(spec)
    model = structure[model_id]
    res_key = (" ", resnum, icode)

    if chain_id is not None:
        if chain_id not in model:
            sys.exit(f"Error: chain '{chain_id}' does not exist in the structure")
        chain = model[chain_id]
        if res_key not in chain:
            sys.exit(f"Error: residue {resnum}{icode.strip()} does not exist in chain {chain_id}")
        return chain[res_key]

    matches = [(chain, chain[res_key]) for chain in model if res_key in chain]
    if not matches:
        sys.exit(f"Error: no residue with number {resnum}{icode.strip()} exists")
    if len(matches) > 1:
        chains_found = ", ".join(c.id for c, _ in matches)
        sys.exit(
            f"Error: residue {resnum}{icode.strip()} exists in several "
            f"chains ({chains_found}). Specify the chain, e.g. "
            f"{matches[0][0].id}:{resnum}{icode.strip()}"
        )
    return matches[0][1]


def emit_report(lines, args, default_stem):
    """Prints `lines` to the console and saves them to a text file (either
    args.output, or an auto-generated '<default_stem>_<structure_id>.txt').
    The file is written in full before printing, so it's always complete
    even if console output gets interrupted."""
    output_path = args.output or f"{default_stem}_{structure_id_from_args(args)}.txt"

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n".join(str(line) for line in lines) + "\n")

    for line in lines:
        print(line)
    print(f"\n[Results written to {output_path}]")
    return output_path


def base_parser(description):
    return argparse.ArgumentParser(description=description)
