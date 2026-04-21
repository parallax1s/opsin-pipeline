"""Minimal PDB ATOM/HETATM parser.

Pure-Python, stdlib-only. Parses fixed-column PDB records per the legacy format
specification. Scope limits (§2 of the spec):
- no mmCIF
- take model 1 on NMR multi-model files
- skip hydrogens
- keep highest-occupancy alt-loc, ties break on alt_loc == "A"
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PDBAtom:
    record: str            # "ATOM" | "HETATM"
    serial: int
    atom_name: str         # stripped
    alt_loc: str           # usually ""
    res_name: str          # 3-letter, stripped
    chain: str             # single char, may be ""
    res_num: int
    ins_code: str          # usually ""
    x: float
    y: float
    z: float
    occupancy: float
    element: str           # right-justified 2-char; stripped

    @property
    def residue_key(self) -> tuple[str, int, str]:
        return (self.chain, self.res_num, self.ins_code)


class PDBParseError(ValueError):
    """Raised on malformed PDB lines; message carries the 1-indexed line number."""


def parse_pdb(source: str | Path) -> list[PDBAtom]:
    """Parse a PDB file path or raw PDB text into a list of non-H atoms.

    Only model 1 is returned (if MODEL/ENDMDL records are present).
    Alt-loc conflicts are resolved by keeping the highest-occupancy record
    per (chain, res_num, ins_code, atom_name); ties go to alt_loc == "A",
    then the lexicographically smallest alt_loc, then file order.
    """
    text = _read_text(source)
    atoms = _parse_lines(text)
    return _resolve_alt_locs(atoms)


def _read_text(source: str | Path) -> str:
    if isinstance(source, Path):
        return source.read_text(encoding="utf-8")
    # heuristic: path-like single-line strings vs actual PDB content
    if "\n" not in source and len(source) < 4096:
        candidate = Path(source)
        if candidate.exists():
            return candidate.read_text(encoding="utf-8")
    return source


def _parse_lines(text: str) -> list[PDBAtom]:
    atoms: list[PDBAtom] = []
    in_model = True  # treat files without MODEL records as "model 1"
    saw_model = False
    for lineno, raw in enumerate(text.splitlines(), start=1):
        if not raw:
            continue
        record = raw[0:6].rstrip()
        if record == "MODEL":
            saw_model = True
            # model number is in cols 11-14; keep only model 1
            try:
                model_num = int(raw[10:14].strip())
            except ValueError as exc:
                raise PDBParseError(f"line {lineno}: malformed MODEL record") from exc
            in_model = model_num == 1
            continue
        if record == "ENDMDL":
            in_model = False
            continue
        if record not in ("ATOM", "HETATM"):
            continue
        if saw_model and not in_model:
            continue
        atoms.append(_parse_atom_line(raw, lineno))
    return atoms


def _parse_atom_line(line: str, lineno: int) -> PDBAtom:
    # Pad short lines so slicing doesn't index-error on files with trimmed columns.
    if len(line) < 80:
        line = line.ljust(80)
    try:
        record = line[0:6].rstrip()
        serial = _int(line[6:11], lineno, "serial")
        atom_name = line[12:16].strip()
        alt_loc = line[16].strip()
        res_name = line[17:20].strip()
        chain = line[21].strip()
        res_num = _int(line[22:26], lineno, "res_num")
        ins_code = line[26].strip()
        x = _float(line[30:38], lineno, "x")
        y = _float(line[38:46], lineno, "y")
        z = _float(line[46:54], lineno, "z")
        occupancy = _float_default(line[54:60], 1.0)
        element = line[76:78].strip() or _guess_element(atom_name)
    except PDBParseError:
        raise
    except Exception as exc:  # defensive: slicing errors, etc.
        raise PDBParseError(f"line {lineno}: malformed {record or 'ATOM/HETATM'} record") from exc

    return PDBAtom(
        record=record,
        serial=serial,
        atom_name=atom_name,
        alt_loc=alt_loc,
        res_name=res_name,
        chain=chain,
        res_num=res_num,
        ins_code=ins_code,
        x=x,
        y=y,
        z=z,
        occupancy=occupancy,
        element=element,
    )


def _int(field: str, lineno: int, label: str) -> int:
    s = field.strip()
    if not s:
        raise PDBParseError(f"line {lineno}: missing integer field {label!r}")
    try:
        return int(s)
    except ValueError as exc:
        raise PDBParseError(f"line {lineno}: {label}={field!r} is not an integer") from exc


def _float(field: str, lineno: int, label: str) -> float:
    s = field.strip()
    if not s:
        raise PDBParseError(f"line {lineno}: missing float field {label!r}")
    try:
        return float(s)
    except ValueError as exc:
        raise PDBParseError(f"line {lineno}: {label}={field!r} is not a float") from exc


def _float_default(field: str, default: float) -> float:
    s = field.strip()
    if not s:
        return default
    try:
        return float(s)
    except ValueError:
        return default


def _guess_element(atom_name: str) -> str:
    """Fallback when the element column is blank (common in older PDBs)."""
    name = atom_name.strip()
    if not name:
        return ""
    # Handle two-letter elements that commonly show up right-justified
    if len(name) >= 2 and name[:2].upper() in {"FE", "MG", "ZN", "CA", "CL", "NA", "BR", "MN", "CU"}:
        return name[:2].upper()
    return name[0].upper()


def _resolve_alt_locs(atoms: list[PDBAtom]) -> list[PDBAtom]:
    # skip hydrogens up front
    heavy = [atom for atom in atoms if atom.element.upper() != "H"]
    # group by (chain, res_num, ins_code, atom_name); pick best alt-loc per group
    groups: dict[tuple[str, int, str, str], list[tuple[int, PDBAtom]]] = {}
    for index, atom in enumerate(heavy):
        key = (atom.chain, atom.res_num, atom.ins_code, atom.atom_name)
        groups.setdefault(key, []).append((index, atom))

    kept: list[PDBAtom] = []
    kept_indices: set[int] = set()
    for _key, entries in groups.items():
        if len(entries) == 1:
            index, atom = entries[0]
            kept_indices.add(index)
            continue
        # sort: highest occupancy, then alt_loc="A" preferred, then smallest alt_loc string,
        # then file order
        def sort_key(item: tuple[int, PDBAtom]) -> tuple[float, int, str, int]:
            idx, a = item
            return (-a.occupancy, 0 if a.alt_loc == "A" else 1, a.alt_loc, idx)

        chosen = sorted(entries, key=sort_key)[0]
        kept_indices.add(chosen[0])

    # emit in original file order
    for index, atom in enumerate(heavy):
        if index in kept_indices:
            kept.append(atom)
    return kept
