"""Retinal-like ligand identification per spec §5.

Allowlist-driven, with an explicit deny-list to document common false positives.
LYR (covalent Lys+retinal Schiff base) is split into Lys-side and retinal-side
atoms so downstream distance math sees only the retinal heavy atoms.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from .pdb import PDBAtom


# Spec §5.1 — allow-list. Unknown HETATM codes are silently ignored.
RETINAL_WHITELIST: frozenset[str] = frozenset(
    {
        "RET",  # all-trans retinal (most common)
        "RTN",  # retinal, alt code
        "LYR",  # retinal Schiff base to Lys (covalent; split in §5.3)
        "A1H",  # all-trans retinal variant (e.g. ChRmine)
        "A1C",  # 13-cis retinal
        "13M",  # 13-cis retinal (alt)
        "9CR",  # 9-cis retinal
    }
)

# Catches vendor/refinement variants like RET_A, RT0, RET2.
# Exclusion guard: BCR is explicitly deny-listed below and must never match here.
_RETINAL_VARIANT_PATTERN = re.compile(r"^RE[T0-9A-Z]$")

# Spec §5.2 — deny-list. Never treated as the pocket ligand.
RETINAL_DENY_LIST: frozenset[str] = frozenset(
    {
        "BCR",   # β-carotene — similar atom count to retinal, common false positive
        "PLM", "CLR", "PC1", "POPC", "OLC",  # lipids
        "OGA", "LMT", "LDA", "LUT",          # detergents
        "ACT", "GOL", "PEG", "TRS",          # buffers
    }
)

# Spec §5.3 — atoms that belong to the Lys side of an LYR covalent residue.
LYR_LYS_ATOMS: frozenset[str] = frozenset({"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"})

# Retinal is C₂₀H₂₈O → 21 heavy atoms. Accept 15-25 as a sanity band.
MIN_HEAVY_ATOMS = 15
MAX_HEAVY_ATOMS = 25


@dataclass(frozen=True)
class LigandMatch:
    res_name: str
    chain: str
    res_num: int
    ins_code: str
    heavy_atom_count: int          # for LYR, this is the retinal-side count after split
    matched_by: str                # "whitelist" | "lyr_schiff_base" | "variant_regex"

    @property
    def residue_key(self) -> tuple[str, int, str]:
        return (self.chain, self.res_num, self.ins_code)


class NoRetinalLigandError(ValueError):
    """Raised when no retinal-like ligand is found. Carries observed HETATM codes."""

    def __init__(self, observed_codes: list[str], path_hint: str = ""):
        self.observed_codes = observed_codes
        hint = f" in {path_hint}" if path_hint else ""
        super().__init__(
            "No retinal-like ligand found"
            f"{hint}. Observed HETATM residue names: "
            + (", ".join(sorted(set(observed_codes))) if observed_codes else "<none>")
            + ". Extend RETINAL_WHITELIST if this is a legitimate retinal variant."
        )


class LigandShapeError(ValueError):
    """Raised when a whitelisted ligand has a heavy-atom count outside the sanity band."""


class MissingChainError(ValueError):
    """Raised when --pdb-chain was requested but that chain holds no ligand matches."""


def identify_ligands(
    atoms: list[PDBAtom],
    *,
    preferred_chain: str | None = None,
    allow_weird_ligand: bool = False,
    path_hint: str = "",
) -> list[LigandMatch]:
    """Return all retinal-like ligand groups in ``atoms``.

    Applies the whitelist/deny-list and heavy-atom sanity band (§5.1, §5.2, §5.4).
    LYR bypasses the sanity band until after splitting (§5.4 bullet 2).

    Raises NoRetinalLigandError if the final list is empty.
    Raises LigandShapeError on heavy-atom-count failures unless allow_weird_ligand.
    Raises MissingChainError if preferred_chain is given but matches none.
    """
    groups = _group_hetatm_residues(atoms)
    observed_codes = [res_name for (res_name, _chain, _num, _ins) in groups.keys()]

    matches: list[LigandMatch] = []
    for (res_name, chain, res_num, ins), group_atoms in sorted(groups.items()):
        if res_name in RETINAL_DENY_LIST:
            continue
        matched_by = _whitelist_match(res_name)
        if matched_by is None:
            continue

        if res_name == "LYR":
            retinal_side = extract_ligand_heavy_atoms(group_atoms, ligand_match_name="LYR")
            heavy_count = len(retinal_side)
        else:
            heavy_count = len(group_atoms)

        if not allow_weird_ligand and not (MIN_HEAVY_ATOMS <= heavy_count <= MAX_HEAVY_ATOMS):
            raise LigandShapeError(
                f"{res_name} at chain {chain!r} resnum {res_num} has {heavy_count} heavy "
                f"atoms after split; expected {MIN_HEAVY_ATOMS}-{MAX_HEAVY_ATOMS}. "
                "Pass allow_weird_ligand=True to skip this check."
            )

        matches.append(
            LigandMatch(
                res_name=res_name,
                chain=chain,
                res_num=res_num,
                ins_code=ins,
                heavy_atom_count=heavy_count,
                matched_by=matched_by,
            )
        )

    if preferred_chain is not None:
        matches_on_chain = [m for m in matches if m.chain == preferred_chain]
        if not matches_on_chain:
            raise MissingChainError(
                f"No retinal-like ligand on chain {preferred_chain!r}. "
                f"Matches were on chains: {sorted({m.chain for m in matches})}"
            )
        matches = matches_on_chain

    if not matches:
        raise NoRetinalLigandError(observed_codes, path_hint=path_hint)

    return matches


def extract_ligand_heavy_atoms(
    ligand_atoms: list[PDBAtom],
    *,
    ligand_match_name: str | None = None,
) -> list[PDBAtom]:
    """Heavy (non-H) ligand atoms, with LYR-side filtering applied if applicable."""
    heavy = [a for a in ligand_atoms if a.element.upper() != "H"]
    res_name = ligand_match_name or (ligand_atoms[0].res_name if ligand_atoms else "")
    if res_name == "LYR":
        return [a for a in heavy if a.atom_name not in LYR_LYS_ATOMS]
    return heavy


def extract_lyr_lys_atoms(ligand_atoms: list[PDBAtom]) -> list[PDBAtom]:
    """Lys-side atoms of an LYR residue (used to locate the Schiff-base Lys Cα)."""
    return [a for a in ligand_atoms if a.atom_name in LYR_LYS_ATOMS and a.element.upper() != "H"]


def _group_hetatm_residues(
    atoms: list[PDBAtom],
) -> dict[tuple[str, str, int, str], list[PDBAtom]]:
    # LYR can appear as ATOM or HETATM across different entries; match both.
    groups: dict[tuple[str, str, int, str], list[PDBAtom]] = {}
    for atom in atoms:
        is_candidate = atom.record == "HETATM" or atom.res_name == "LYR"
        if not is_candidate:
            continue
        key = (atom.res_name, atom.chain, atom.res_num, atom.ins_code)
        groups.setdefault(key, []).append(atom)
    return groups


def _whitelist_match(res_name: str) -> str | None:
    if res_name in RETINAL_DENY_LIST:
        return None
    if res_name == "LYR":
        return "lyr_schiff_base"
    if res_name in RETINAL_WHITELIST:
        return "whitelist"
    if _RETINAL_VARIANT_PATTERN.match(res_name):
        return "variant_regex"
    return None
