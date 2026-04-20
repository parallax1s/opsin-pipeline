import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from opsin_pipeline.structure.pdb import PDBParseError, parse_pdb


FIXTURE_DIR = Path(__file__).parent / "fixtures" / "structures"
SYNTHETIC = FIXTURE_DIR / "synthetic_mini.pdb"


class PDBParserTests(unittest.TestCase):
    def test_parses_synthetic_fixture(self):
        atoms = parse_pdb(SYNTHETIC)

        # 8 protein CA + 20 retinal C + 1 glycerol O = 29 heavy atoms (hydrogen dropped)
        self.assertEqual(len(atoms), 29)
        protein = [a for a in atoms if a.record == "ATOM"]
        hetatms = [a for a in atoms if a.record == "HETATM"]
        self.assertEqual(len(protein), 8)
        self.assertEqual(len(hetatms), 21)

    def test_skips_hydrogens(self):
        atoms = parse_pdb(SYNTHETIC)
        self.assertFalse(any(a.element.upper() == "H" for a in atoms))
        # the explicit H on SER 7 in the fixture must not survive
        self.assertFalse(any(a.atom_name == "H" for a in atoms))

    def test_coordinates_and_residues_parse_correctly(self):
        atoms = parse_pdb(SYNTHETIC)
        by_key = {(a.chain, a.res_num, a.atom_name): a for a in atoms}

        ala = by_key[("A", 1, "CA")]
        self.assertEqual(ala.res_name, "ALA")
        self.assertAlmostEqual(ala.x, 2.0)
        self.assertAlmostEqual(ala.y, 2.0)

        ret_c13 = by_key[("A", 101, "C13")]
        self.assertEqual(ret_c13.record, "HETATM")
        self.assertEqual(ret_c13.res_name, "RET")
        self.assertAlmostEqual(ret_c13.x, 12.0)

    def test_accepts_text_input(self):
        text = SYNTHETIC.read_text(encoding="utf-8")
        from_text = parse_pdb(text)
        from_path = parse_pdb(SYNTHETIC)
        self.assertEqual(
            [(a.record, a.serial) for a in from_text],
            [(a.record, a.serial) for a in from_path],
        )

    def test_malformed_line_raises_with_line_number(self):
        text = (
            "ATOM      1  CA  ALA A   1       2.000   2.000   0.000  1.00 20.00           C\n"
            "ATOM   BAD                                                                    C\n"
        )
        with self.assertRaises(PDBParseError) as ctx:
            parse_pdb(text)
        self.assertIn("line 2", str(ctx.exception))

    def test_nmr_multi_model_keeps_model_one(self):
        text = (
            "MODEL        1\n"
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C\n"
            "ENDMDL\n"
            "MODEL        2\n"
            "ATOM      1  CA  ALA A   1      99.000  99.000  99.000  1.00 20.00           C\n"
            "ENDMDL\n"
        )
        atoms = parse_pdb(text)
        self.assertEqual(len(atoms), 1)
        self.assertAlmostEqual(atoms[0].x, 0.0)

    def test_alt_loc_keeps_highest_occupancy(self):
        text = (
            "ATOM      1  CA AALA A   1       0.000   0.000   0.000  0.30 20.00           C\n"
            "ATOM      2  CA BALA A   1       9.000   9.000   9.000  0.70 20.00           C\n"
        )
        atoms = parse_pdb(text)
        self.assertEqual(len(atoms), 1)
        self.assertAlmostEqual(atoms[0].x, 9.0)
        self.assertEqual(atoms[0].alt_loc, "B")

    def test_alt_loc_tiebreak_prefers_alt_a(self):
        text = (
            "ATOM      1  CA AALA A   1       1.000   1.000   1.000  0.50 20.00           C\n"
            "ATOM      2  CA BALA A   1       2.000   2.000   2.000  0.50 20.00           C\n"
        )
        atoms = parse_pdb(text)
        self.assertEqual(len(atoms), 1)
        self.assertEqual(atoms[0].alt_loc, "A")


if __name__ == "__main__":
    unittest.main()
