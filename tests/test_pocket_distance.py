import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from opsin_pipeline.structure.ligands import identify_ligands
from opsin_pipeline.structure.pdb import parse_pdb
from opsin_pipeline.structure.pocket import (
    DEFAULT_MEDIUM_MAX_A,
    DEFAULT_STRONG_MAX_A,
    EmptyPocketError,
    compute_pocket,
)


FIXTURE_DIR = Path(__file__).parent / "fixtures" / "structures"
SYNTHETIC = FIXTURE_DIR / "synthetic_mini.pdb"
SYNTHETIC_LYR = FIXTURE_DIR / "synthetic_lyr.pdb"


class SyntheticDistanceTests(unittest.TestCase):
    def _pocket(self, **overrides):
        atoms = parse_pdb(SYNTHETIC)
        matches = identify_ligands(atoms)
        return compute_pocket(
            atoms,
            matches,
            scaffold_name="synthetic",
            pdb_path=str(SYNTHETIC),
            **overrides,
        )

    def test_band_counts_match_fixture_design(self):
        pocket = self._pocket()

        bands = {}
        for residue in pocket.pocket_residues:
            bands.setdefault(residue.band, []).append(residue)
        # fixture has 3 strong (<=4 A), 2 medium (4-5.5 A), 1 none (6.0 A exactly,
        # within the default 6.0 A cutoff but beyond medium_max 5.5), and 2 residues
        # outside the cutoff (7.0 and 10.0 A) which are excluded entirely.
        self.assertEqual(len(bands.get("strong", [])), 3)
        self.assertEqual(len(bands.get("medium", [])), 2)
        self.assertEqual(len(bands.get("none", [])), 1)

    def test_minimum_distances_match_hand_computation(self):
        pocket = self._pocket()
        by_resnum = {r.pdb_resnum: r for r in pocket.pocket_residues}

        self.assertAlmostEqual(by_resnum[1].min_distance_A, 2.0, places=3)
        self.assertAlmostEqual(by_resnum[2].min_distance_A, 2.5, places=3)
        self.assertAlmostEqual(by_resnum[3].min_distance_A, 3.5, places=3)
        self.assertAlmostEqual(by_resnum[4].min_distance_A, 4.5, places=3)
        self.assertAlmostEqual(by_resnum[5].min_distance_A, 5.0, places=3)

    def test_residues_outside_cutoff_are_dropped(self):
        pocket = self._pocket()
        # residues 7 (7.0 A) and 8 (10.0 A) exceed the default 6.0 A cutoff and are
        # dropped entirely; residue 6 (6.0 A exactly) is kept with band="none".
        resnums = {r.pdb_resnum for r in pocket.pocket_residues}
        self.assertNotIn(7, resnums)
        self.assertNotIn(8, resnums)
        self.assertIn(6, resnums)

    def test_cutoff_can_be_raised_to_include_outer_residues(self):
        pocket = self._pocket(cutoff_A=8.0)

        resnums = {r.pdb_resnum for r in pocket.pocket_residues}
        self.assertIn(6, resnums)
        self.assertIn(7, resnums)
        self.assertNotIn(8, resnums)  # residue 8 is 10 A away, still outside

    def test_band_thresholds_are_tunable(self):
        pocket = self._pocket(strong_max_A=3.0, medium_max_A=4.6)

        by_resnum = {r.pdb_resnum: r for r in pocket.pocket_residues}
        # now resnum 3 (3.5 A) should be medium, not strong
        self.assertEqual(by_resnum[3].band, "medium")
        # resnum 1 (2.0 A) still strong
        self.assertEqual(by_resnum[1].band, "strong")

    def test_defaults_sanity(self):
        # guard against accidental constant changes
        self.assertEqual(DEFAULT_STRONG_MAX_A, 4.0)
        self.assertEqual(DEFAULT_MEDIUM_MAX_A, 5.5)

    def test_pocket_map_carries_sha256_and_ligand_match(self):
        pocket = self._pocket()

        self.assertTrue(pocket.pdb_sha256)
        self.assertEqual(len(pocket.ligand_matches), 1)
        self.assertEqual(pocket.ligand_matches[0].res_name, "RET")


class EmptyPocketTests(unittest.TestCase):
    def test_zero_residues_within_cutoff_raises(self):
        atoms = parse_pdb(SYNTHETIC)
        matches = identify_ligands(atoms)
        with self.assertRaises(EmptyPocketError):
            compute_pocket(
                atoms,
                matches,
                scaffold_name="synthetic",
                pdb_path=str(SYNTHETIC),
                cutoff_A=0.5,
            )


class LYRSchiffBaseTests(unittest.TestCase):
    def test_lyr_residue_always_emitted_as_schiff_base_strong(self):
        atoms = parse_pdb(SYNTHETIC_LYR)
        matches = identify_ligands(atoms)

        pocket = compute_pocket(
            atoms,
            matches,
            scaffold_name="bacteriorhodopsin_test",
            pdb_path=str(SYNTHETIC_LYR),
        )

        lyr_residues = [r for r in pocket.pocket_residues if r.role == "schiff_base_linkage"]
        self.assertEqual(len(lyr_residues), 1)
        self.assertEqual(lyr_residues[0].pdb_resnum, 216)
        self.assertEqual(lyr_residues[0].band, "strong")
        self.assertEqual(lyr_residues[0].min_distance_A, 0.0)
        self.assertEqual(lyr_residues[0].res_name_one, "K")

    def test_lyr_lys_side_atoms_are_not_counted_as_ligand(self):
        # ALA 85 in the LYR fixture is at (3, 0, 0); the Lys-side CA of LYR 216 is at
        # (1.5, 0, 0). If Lys-side atoms were wrongly counted as the ligand, ALA's
        # nearest ligand distance would be 1.5 A. In reality the retinal C1 starts at
        # (10, 0, 0), so ALA-to-retinal distance must be at least 7 A (well outside
        # the 6.0 A cutoff) and ALA must not appear in the pocket.
        atoms = parse_pdb(SYNTHETIC_LYR)
        matches = identify_ligands(atoms)

        pocket = compute_pocket(
            atoms,
            matches,
            scaffold_name="br_test",
            pdb_path=str(SYNTHETIC_LYR),
        )

        self.assertFalse(any(r.res_name_three == "ALA" for r in pocket.pocket_residues))


if __name__ == "__main__":
    unittest.main()
