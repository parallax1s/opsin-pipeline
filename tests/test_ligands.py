import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from opsin_pipeline.structure.ligands import (
    LigandShapeError,
    LYR_LYS_ATOMS,
    MissingChainError,
    NoRetinalLigandError,
    extract_ligand_heavy_atoms,
    extract_lyr_lys_atoms,
    identify_ligands,
)
from opsin_pipeline.structure.pdb import parse_pdb


FIXTURE_DIR = Path(__file__).parent / "fixtures" / "structures"
SYNTHETIC = FIXTURE_DIR / "synthetic_mini.pdb"
SYNTHETIC_LYR = FIXTURE_DIR / "synthetic_lyr.pdb"


def _pdb_atom(
    record: str,
    serial: int,
    atom_name: str,
    res_name: str,
    chain: str,
    res_num: int,
    x: float,
    y: float = 0.0,
    z: float = 0.0,
    element: str = "C",
) -> str:
    """Build a single legal PDB ATOM/HETATM line at fixed columns."""
    assert len(record) <= 6
    name_field = atom_name.ljust(4) if len(atom_name) <= 3 else atom_name[:4]
    return (
        f"{record:<6}"
        f"{serial:>5}"
        f" "
        f"{name_field:<4}"
        f" "
        f"{res_name:<3}"
        f" "
        f"{chain:<1}"
        f"{res_num:>4}"
        f" "
        f"   "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
        f"  1.00 20.00          "
        f"{element:>2}"
    )


def _ligand_block(res_name: str, chain: str, res_num: int, n_atoms: int) -> str:
    return (
        "\n".join(
            _pdb_atom("HETATM", i + 1, f"C{i+1}", res_name, chain, res_num, float(i))
            for i in range(n_atoms)
        )
        + "\n"
    )


class WhitelistDenylistTests(unittest.TestCase):
    def test_ret_is_matched_by_whitelist(self):
        atoms = parse_pdb(SYNTHETIC)

        matches = identify_ligands(atoms)

        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0].res_name, "RET")
        self.assertEqual(matches[0].matched_by, "whitelist")
        self.assertEqual(matches[0].heavy_atom_count, 20)

    def test_glycerol_is_not_matched(self):
        atoms = parse_pdb(SYNTHETIC)

        matches = identify_ligands(atoms)

        self.assertFalse(any(m.res_name == "GOL" for m in matches))

    def test_bcr_is_hard_denied_even_if_heavy_atom_count_looks_ok(self):
        # 21 heavy atoms (inside the sanity band) but BCR is explicitly deny-listed
        text = _ligand_block("BCR", "A", 100, n_atoms=21)

        atoms = parse_pdb(text)
        with self.assertRaises(NoRetinalLigandError) as ctx:
            identify_ligands(atoms)
        self.assertIn("BCR", str(ctx.exception))

    def test_no_matching_ligand_lists_observed_codes(self):
        atoms = parse_pdb(
            "HETATM    1  O1  GOL A 201      50.000  50.000  50.000  1.00 40.00           O\n"
        )

        with self.assertRaises(NoRetinalLigandError) as ctx:
            identify_ligands(atoms)
        self.assertIn("GOL", str(ctx.exception))
        self.assertIn("RETINAL_WHITELIST", str(ctx.exception))

    def test_variant_regex_catches_unknown_retinal_suffix(self):
        # "REA" isn't in the explicit whitelist but matches ^RE[T0-9A-Z]$
        text = _ligand_block("REA", "A", 100, n_atoms=20)

        matches = identify_ligands(parse_pdb(text))

        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0].matched_by, "variant_regex")


class HeavyAtomSanityTests(unittest.TestCase):
    def test_too_few_atoms_raises_shape_error(self):
        with self.assertRaises(LigandShapeError):
            identify_ligands(parse_pdb(_ligand_block("RET", "A", 100, n_atoms=5)))

    def test_too_many_atoms_raises_shape_error(self):
        with self.assertRaises(LigandShapeError):
            identify_ligands(parse_pdb(_ligand_block("RET", "A", 100, n_atoms=30)))

    def test_allow_weird_ligand_bypasses_sanity_band(self):
        text = _ligand_block("RET", "A", 100, n_atoms=5)
        matches = identify_ligands(parse_pdb(text), allow_weird_ligand=True)
        self.assertEqual(len(matches), 1)


class LYRHandlingTests(unittest.TestCase):
    def test_lyr_matches_as_schiff_base_and_split_count_is_retinal_side(self):
        atoms = parse_pdb(SYNTHETIC_LYR)

        matches = identify_ligands(atoms)

        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0].res_name, "LYR")
        self.assertEqual(matches[0].matched_by, "lyr_schiff_base")
        # LYR fixture has 9 Lys-side atoms + 20 retinal-side atoms; after split, the
        # reported heavy count is the retinal side only
        self.assertEqual(matches[0].heavy_atom_count, 20)

    def test_lyr_split_retinal_heavy_atoms_excludes_lys_side(self):
        atoms = parse_pdb(SYNTHETIC_LYR)
        lyr_atoms = [a for a in atoms if a.res_name == "LYR"]

        retinal_side = extract_ligand_heavy_atoms(lyr_atoms, ligand_match_name="LYR")
        lys_side = extract_lyr_lys_atoms(lyr_atoms)

        retinal_names = {a.atom_name for a in retinal_side}
        lys_names = {a.atom_name for a in lys_side}
        self.assertTrue(retinal_names.isdisjoint(LYR_LYS_ATOMS))
        self.assertTrue(lys_names <= LYR_LYS_ATOMS)
        self.assertEqual(len(retinal_side), 20)
        self.assertEqual(len(lys_side), 9)

    def test_lyr_with_unsplit_count_above_25_still_passes_because_of_bypass(self):
        # combined LYR with 9 Lys + 20 retinal = 29 heavy atoms; the bypass means the
        # sanity band only evaluates the 20 retinal atoms, which is in range
        atoms = parse_pdb(SYNTHETIC_LYR)

        # should not raise
        matches = identify_ligands(atoms)

        self.assertEqual(matches[0].heavy_atom_count, 20)


class ChainSelectionTests(unittest.TestCase):
    def test_preferred_chain_filters_to_that_chain(self):
        # two RETs on different chains
        lines = []
        for chain, resnum, x_offset in (("A", 101, 0.0), ("B", 201, 100.0)):
            for i in range(20):
                lines.append(
                    _pdb_atom(
                        "HETATM",
                        i + 1,
                        f"C{i+1}",
                        "RET",
                        chain,
                        resnum,
                        x_offset + i,
                    )
                )
        atoms = parse_pdb("\n".join(lines) + "\n")

        on_b = identify_ligands(atoms, preferred_chain="B")

        self.assertEqual(len(on_b), 1)
        self.assertEqual(on_b[0].chain, "B")

    def test_missing_chain_raises(self):
        atoms = parse_pdb(SYNTHETIC)
        with self.assertRaises(MissingChainError):
            identify_ligands(atoms, preferred_chain="Z")


if __name__ == "__main__":
    unittest.main()
