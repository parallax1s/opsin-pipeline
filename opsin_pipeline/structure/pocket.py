"""Distance-based retinal-pocket computation per spec §4 and §7.

Given parsed PDB atoms and identified ligand matches, compute each residue's
minimum distance to any ligand heavy atom and assign a band. Output: a
``PocketMap`` that downstream code (position_map annotation, scaffold ingest)
can consume.
"""
from __future__ import annotations

import hashlib
import json
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
    role: str | None = None                 # "schiff_base_linkage" for LYR Lys center
    seq_index: int | None = None            # scaffold-numbered; only set when mapping applied
    mapping_note: str | None = None         # e.g. "PDB K85 vs scaffold position 87 (K)"


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


# ---- JSON I/O ----

def pocket_map_to_dict(pocket_map: PocketMap) -> dict:
    return {
        "scaffold_name": pocket_map.scaffold_name,
        "pdb_path": pocket_map.pdb_path,
        "pdb_sha256": pocket_map.pdb_sha256,
        "pdb_chain": pocket_map.pdb_chain,
        "cutoff_A": pocket_map.cutoff_A,
        "thresholds_A": list(pocket_map.thresholds_A),
        "ligand_matches": [
            {
                "res_name": m.res_name,
                "chain": m.chain,
                "res_num": m.res_num,
                "ins_code": m.ins_code,
                "heavy_atom_count": m.heavy_atom_count,
                "matched_by": m.matched_by,
            }
            for m in pocket_map.ligand_matches
        ],
        "pocket_residues": [
            {
                "pdb_chain": r.pdb_chain,
                "pdb_resnum": r.pdb_resnum,
                "pdb_ins_code": r.pdb_ins_code,
                "res_name_three": r.res_name_three,
                "res_name_one": r.res_name_one,
                "min_distance_A": r.min_distance_A,
                "closest_ligand_atom": r.closest_ligand_atom,
                "closest_ligand_id": r.closest_ligand_id,
                "closest_ligand_resnum": r.closest_ligand_resnum,
                "band": r.band,
                "role": r.role,
                "seq_index": r.seq_index,
                "mapping_note": r.mapping_note,
            }
            for r in pocket_map.pocket_residues
        ],
    }


def pocket_map_from_dict(data: dict) -> PocketMap:
    ligand_matches = [
        LigandMatch(
            res_name=str(m["res_name"]),
            chain=str(m["chain"]),
            res_num=int(m["res_num"]),
            ins_code=str(m.get("ins_code", "")),
            heavy_atom_count=int(m["heavy_atom_count"]),
            matched_by=str(m["matched_by"]),
        )
        for m in data.get("ligand_matches", [])
    ]
    pocket_residues = [
        PocketResidue(
            pdb_chain=str(r["pdb_chain"]),
            pdb_resnum=int(r["pdb_resnum"]),
            pdb_ins_code=str(r.get("pdb_ins_code", "")),
            res_name_three=str(r["res_name_three"]),
            res_name_one=str(r["res_name_one"]),
            min_distance_A=float(r["min_distance_A"]),
            closest_ligand_atom=str(r["closest_ligand_atom"]),
            closest_ligand_id=str(r["closest_ligand_id"]),
            closest_ligand_resnum=int(r["closest_ligand_resnum"]),
            band=str(r["band"]),
            role=r.get("role"),
            seq_index=r.get("seq_index"),
            mapping_note=r.get("mapping_note"),
        )
        for r in data.get("pocket_residues", [])
    ]
    thresholds = data.get("thresholds_A", [DEFAULT_STRONG_MAX_A, DEFAULT_MEDIUM_MAX_A])
    return PocketMap(
        scaffold_name=str(data["scaffold_name"]),
        pdb_path=str(data.get("pdb_path", "")),
        pdb_sha256=str(data.get("pdb_sha256", "")),
        pdb_chain=str(data.get("pdb_chain", "")),
        ligand_matches=ligand_matches,
        pocket_residues=pocket_residues,
        cutoff_A=float(data.get("cutoff_A", DEFAULT_CUTOFF_A)),
        thresholds_A=(float(thresholds[0]), float(thresholds[1])),
    )


def write_pocket_map(pocket_map: PocketMap, path: str | Path) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(pocket_map_to_dict(pocket_map), indent=2) + "\n", encoding="utf-8")
    return p


def read_pocket_map(path: str | Path) -> PocketMap:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return pocket_map_from_dict(data)


# ---- Mapping application (PDB numbering -> scaffold seq_index) ----

def apply_offset_mapping(
    pocket_map: PocketMap,
    *,
    offset: int,
    scaffold_sequence: str | None = None,
    strict: bool = False,
) -> PocketMap:
    """Return a copy with ``seq_index = pdb_resnum + offset`` on every residue.

    If ``scaffold_sequence`` is provided:

    - the resulting ``seq_index`` must land inside ``[1, len(scaffold_sequence)]``.
      Out-of-bounds indices are flagged via ``mapping_note`` + ``review_needed=true``
      (or raised under ``strict=True``) — never silently written.
    - the residue's one-letter code is checked against the scaffold at that index;
      mismatches are handled the same way.

    Insertion codes are rejected outright: the scalar offset mode cannot represent
    non-scalar numbering, so any pocket residue carrying an insertion code raises
    ``ValueError``. A per-residue mapping table is deferred — see spec §12b.
    """
    new_residues: list[PocketResidue] = []
    for r in pocket_map.pocket_residues:
        if r.pdb_ins_code:
            raise ValueError(
                f"apply_offset_mapping cannot handle insertion codes "
                f"(chain {r.pdb_chain} resnum {r.pdb_resnum}{r.pdb_ins_code}). "
                "Non-scalar numbering requires a per-residue mapping table, "
                "which is deferred — see spec §12b."
            )
        seq_index = r.pdb_resnum + offset
        note: str | None = None
        if scaffold_sequence is not None:
            if seq_index < 1 or seq_index > len(scaffold_sequence):
                msg = (
                    f"seq_index {seq_index} (from PDB chain {r.pdb_chain} resnum "
                    f"{r.pdb_resnum} + offset {offset}) is outside scaffold length "
                    f"{len(scaffold_sequence)}"
                )
                if strict:
                    raise ValueError(msg)
                note = msg
            else:
                scaffold_aa = scaffold_sequence[seq_index - 1].upper()
                if r.res_name_one.upper() != scaffold_aa and r.res_name_one != "X":
                    msg = (
                        f"PDB chain {r.pdb_chain} resnum {r.pdb_resnum} is "
                        f"{r.res_name_one} ({r.res_name_three}) but scaffold position "
                        f"{seq_index} is {scaffold_aa}"
                    )
                    if strict:
                        raise ValueError(msg)
                    note = msg
        new_residues.append(
            PocketResidue(
                pdb_chain=r.pdb_chain,
                pdb_resnum=r.pdb_resnum,
                pdb_ins_code=r.pdb_ins_code,
                res_name_three=r.res_name_three,
                res_name_one=r.res_name_one,
                min_distance_A=r.min_distance_A,
                closest_ligand_atom=r.closest_ligand_atom,
                closest_ligand_id=r.closest_ligand_id,
                closest_ligand_resnum=r.closest_ligand_resnum,
                band=r.band,
                role=r.role,
                seq_index=seq_index,
                mapping_note=note,
            )
        )
    return PocketMap(
        scaffold_name=pocket_map.scaffold_name,
        pdb_path=pocket_map.pdb_path,
        pdb_sha256=pocket_map.pdb_sha256,
        pdb_chain=pocket_map.pdb_chain,
        ligand_matches=list(pocket_map.ligand_matches),
        pocket_residues=new_residues,
        cutoff_A=pocket_map.cutoff_A,
        thresholds_A=pocket_map.thresholds_A,
    )
