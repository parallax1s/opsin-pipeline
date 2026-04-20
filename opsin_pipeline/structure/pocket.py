"""Distance-based retinal-pocket computation per spec §4 and §7.

Given parsed PDB atoms and identified ligand matches, compute each residue's
minimum distance to any ligand heavy atom and assign a band. Output: a
``PocketMap`` that downstream code (position_map annotation, scaffold ingest)
can consume.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path

from .ligands import (
    LYR_LYS_ATOMS,
    LigandMatch,
    extract_ligand_heavy_atoms,
)
from .pdb import PDBAtom


# Spec §7.1 defaults; exposed as module constants so they're easy to tune.
DEFAULT_CUTOFF_A = 6.0
DEFAULT_STRONG_MAX_A = 4.0
DEFAULT_MEDIUM_MAX_A = 5.5

# ASCII three-letter -> one-letter. Unknown residues get "X".
_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # Common non-standard / PDB-specific residues that still sit on the protein chain
    "MSE": "M",  # selenomethionine
    "SEC": "U",
    "PYL": "O",
    "LYR": "K",  # Lys side of the retinal Schiff base
}


@dataclass(frozen=True)
class PocketResidue:
    pdb_chain: str
    pdb_resnum: int
    pdb_ins_code: str
    res_name_three: str
    res_name_one: str
    min_distance_A: float
    closest_ligand_atom: str
    closest_ligand_id: str
    closest_ligand_resnum: int
    band: str  # "strong" | "medium" | "none"
    role: str | None = None  # "schiff_base_linkage" for LYR Lys center, else None


@dataclass(frozen=True)
class PocketMap:
    scaffold_name: str
    pdb_path: str
    pdb_sha256: str
    pdb_chain: str
    ligand_matches: list[LigandMatch]
    pocket_residues: list[PocketResidue]
    cutoff_A: float
    thresholds_A: tuple[float, float]   # (strong_max, medium_max)


class EmptyPocketError(ValueError):
    """Raised when computation completes but no residues fall within cutoff."""


def compute_pocket(
    atoms: list[PDBAtom],
    ligand_matches: list[LigandMatch],
    *,
    scaffold_name: str,
    pdb_path: str | Path,
    pdb_chain: str | None = None,
    cutoff_A: float = DEFAULT_CUTOFF_A,
    strong_max_A: float = DEFAULT_STRONG_MAX_A,
    medium_max_A: float = DEFAULT_MEDIUM_MAX_A,
) -> PocketMap:
    """Build a PocketMap from parsed atoms + identified ligand matches.

    If ``pdb_chain`` is given, only protein residues on that chain are considered
    pocket candidates. The ligand itself may live on any chain already filtered
    upstream by ``identify_ligands(preferred_chain=...)``.
    """
    if not ligand_matches:
        raise ValueError("compute_pocket requires at least one ligand match")
    if strong_max_A > medium_max_A:
        raise ValueError("strong_max_A must be <= medium_max_A")

    ligand_heavy_atoms = _collect_ligand_heavy_atoms(atoms, ligand_matches)
    protein_residues = _group_protein_residues(atoms, pdb_chain=pdb_chain)

    pocket_residues: list[PocketResidue] = []
    for key, residue_atoms in sorted(protein_residues.items()):
        chain, resnum, ins_code = key
        min_distance, closest_lig_atom = _min_distance_to_ligand(residue_atoms, ligand_heavy_atoms)
        if min_distance > cutoff_A:
            continue
        band = _assign_band(min_distance, strong_max_A=strong_max_A, medium_max_A=medium_max_A)
        three = residue_atoms[0].res_name
        pocket_residues.append(
            PocketResidue(
                pdb_chain=chain,
                pdb_resnum=resnum,
                pdb_ins_code=ins_code,
                res_name_three=three,
                res_name_one=_THREE_TO_ONE.get(three, "X"),
                min_distance_A=round(min_distance, 3),
                closest_ligand_atom=closest_lig_atom.atom_name,
                closest_ligand_id=closest_lig_atom.res_name,
                closest_ligand_resnum=closest_lig_atom.res_num,
                band=band,
                role=None,
            )
        )

    # §5.3 / §12.4 — LYR Lys-side contributes a Schiff-base pocket residue regardless of cutoff
    for lyr_residue in _lyr_schiff_base_residues(atoms, ligand_matches):
        pocket_residues = [r for r in pocket_residues if r.pdb_resnum != lyr_residue.pdb_resnum or r.pdb_chain != lyr_residue.pdb_chain]
        pocket_residues.append(lyr_residue)

    if not pocket_residues:
        raise EmptyPocketError(
            f"No residues within {cutoff_A} A of retinal in {pdb_path}. "
            "Check chain selection or cutoff."
        )

    pocket_residues.sort(key=lambda r: (r.pdb_chain, r.pdb_resnum, r.pdb_ins_code))

    return PocketMap(
        scaffold_name=scaffold_name,
        pdb_path=str(pdb_path),
        pdb_sha256=_sha256_file(pdb_path),
        pdb_chain=pdb_chain or "",
        ligand_matches=list(ligand_matches),
        pocket_residues=pocket_residues,
        cutoff_A=cutoff_A,
        thresholds_A=(strong_max_A, medium_max_A),
    )


def _collect_ligand_heavy_atoms(
    atoms: list[PDBAtom], ligand_matches: list[LigandMatch]
) -> list[PDBAtom]:
    wanted_keys = {m.residue_key + (m.res_name,) for m in ligand_matches}
    collected: list[PDBAtom] = []
    for atom in atoms:
        # LYR may sit on ATOM records, not HETATM; match on (chain, resnum, ins, res_name)
        atom_key = (atom.chain, atom.res_num, atom.ins_code, atom.res_name)
        if atom_key not in wanted_keys:
            continue
        if atom.res_name == "LYR" and atom.atom_name in LYR_LYS_ATOMS:
            continue  # exclude Lys-side atoms from ligand heavy-atom set
        if atom.element.upper() == "H":
            continue
        collected.append(atom)
    return collected


def _group_protein_residues(
    atoms: list[PDBAtom], *, pdb_chain: str | None
) -> dict[tuple[str, int, str], list[PDBAtom]]:
    grouped: dict[tuple[str, int, str], list[PDBAtom]] = {}
    for atom in atoms:
        if atom.record != "ATOM":
            continue
        if atom.res_name == "LYR":
            continue  # LYR is handled separately as a Schiff-base contribution
        if atom.res_name not in _THREE_TO_ONE:
            continue
        if pdb_chain is not None and atom.chain != pdb_chain:
            continue
        grouped.setdefault((atom.chain, atom.res_num, atom.ins_code), []).append(atom)
    return grouped


def _min_distance_to_ligand(
    residue_atoms: list[PDBAtom], ligand_atoms: list[PDBAtom]
) -> tuple[float, PDBAtom]:
    best = float("inf")
    best_atom: PDBAtom | None = None
    for res_atom in residue_atoms:
        for lig_atom in ligand_atoms:
            d = _distance(res_atom, lig_atom)
            if d < best:
                best = d
                best_atom = lig_atom
    assert best_atom is not None  # ligand_atoms non-empty guaranteed by caller
    return best, best_atom


def _distance(a: PDBAtom, b: PDBAtom) -> float:
    dx = a.x - b.x
    dy = a.y - b.y
    dz = a.z - b.z
    return (dx * dx + dy * dy + dz * dz) ** 0.5


def _assign_band(distance: float, *, strong_max_A: float, medium_max_A: float) -> str:
    if distance <= strong_max_A:
        return "strong"
    if distance <= medium_max_A:
        return "medium"
    return "none"


def _lyr_schiff_base_residues(
    atoms: list[PDBAtom], ligand_matches: list[LigandMatch]
) -> list[PocketResidue]:
    """Emit a PocketResidue for the Lys side of every matched LYR.

    Per spec §12.4: band="strong", distance=0.0, role="schiff_base_linkage".
    """
    results: list[PocketResidue] = []
    for match in ligand_matches:
        if match.res_name != "LYR":
            continue
        results.append(
            PocketResidue(
                pdb_chain=match.chain,
                pdb_resnum=match.res_num,
                pdb_ins_code=match.ins_code,
                res_name_three="LYS",
                res_name_one="K",
                min_distance_A=0.0,
                closest_ligand_atom="N1",  # the Schiff-base nitrogen, nominal
                closest_ligand_id="LYR",
                closest_ligand_resnum=match.res_num,
                band="strong",
                role="schiff_base_linkage",
            )
        )
    return results


def _sha256_file(path: str | Path) -> str:
    p = Path(path)
    if not p.exists():
        return ""
    return hashlib.sha256(p.read_bytes()).hexdigest()
