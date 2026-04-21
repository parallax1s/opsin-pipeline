"""End-to-end integration test against a stripped real bacteriorhodopsin PDB.

The fixture ``br_1c3w_strip.pdb`` is PDB 1C3W reduced to residues within 8 A of
the retinal (HETATM code RET). Expected pocket residues come from the published
BR retinal-binding site: K216 (Schiff-base Lys), D85 (primary counterion), D212
(proton shuttle), W86 / W182 / W189 (cage tryptophans), Y185 (counterion-
adjacent tyrosine), T89 (spectral tuning threonine). Tests assert the parser
plus pocket computation recover these residues with sensible band assignments.
"""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from opsin_pipeline.structure.ligands import identify_ligands
from opsin_pipeline.structure.pdb import parse_pdb
from opsin_pipeline.structure.pocket import compute_pocket


FIXTURE = Path(__file__).parent / "fixtures" / "structures" / "br_1c3w_strip.pdb"


class BRPocketIntegrationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        atoms = parse_pdb(FIXTURE)
        matches = identify_ligands(atoms)
        cls.pocket = compute_pocket(
            atoms,
            matches,
            scaffold_name="bacteriorhodopsin",
            pdb_path=str(FIXTURE),
            pdb_chain="A",
        )
        cls.by_resnum = {r.pdb_resnum: r for r in cls.pocket.pocket_residues}

    def test_retinal_is_identified(self):
        matches = self.pocket.ligand_matches
        self.assertEqual(len(matches), 1)
        self.assertEqual(matches[0].res_name, "RET")
        self.assertEqual(matches[0].heavy_atom_count, 20)  # retinal C1..C20

    def test_schiff_base_lysine_is_in_strong_band(self):
        self.assertIn(216, self.by_resnum)
        k216 = self.by_resnum[216]
        self.assertEqual(k216.res_name_three, "LYS")
        self.assertEqual(k216.band, "strong")
        self.assertLess(k216.min_distance_A, 2.0)  # Schiff base C-N distance

    def test_primary_counterion_d85_is_in_pocket(self):
        self.assertIn(85, self.by_resnum)
        d85 = self.by_resnum[85]
        self.assertEqual(d85.res_name_three, "ASP")
        self.assertIn(d85.band, {"strong", "medium"})

    def test_proton_shuttle_d212_is_strong(self):
        self.assertIn(212, self.by_resnum)
        d212 = self.by_resnum[212]
        self.assertEqual(d212.res_name_three, "ASP")
        self.assertEqual(d212.band, "strong")

    def test_cage_tryptophans_are_all_pocket(self):
        for resnum in (86, 182, 189):
            self.assertIn(resnum, self.by_resnum, f"W{resnum} missing from pocket")
            self.assertEqual(self.by_resnum[resnum].res_name_three, "TRP")
            self.assertIn(self.by_resnum[resnum].band, {"strong", "medium"})

    def test_y185_is_pocket(self):
        self.assertIn(185, self.by_resnum)
        y185 = self.by_resnum[185]
        self.assertEqual(y185.res_name_three, "TYR")
        self.assertEqual(y185.band, "strong")

    def test_t89_spectral_tuning_is_pocket(self):
        self.assertIn(89, self.by_resnum)
        t89 = self.by_resnum[89]
        self.assertEqual(t89.res_name_three, "THR")
        self.assertEqual(t89.band, "strong")

    def test_no_lipid_or_buffer_residues_appear_as_ligand(self):
        # 1C3W has SQU (squalene) and LI1 lipids nearby; they must not be picked
        # up as retinal ligands under any circumstance.
        for match in self.pocket.ligand_matches:
            self.assertNotIn(match.res_name, {"SQU", "LI1", "HOH"})


if __name__ == "__main__":
    unittest.main()
