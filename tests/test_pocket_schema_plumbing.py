"""Ingest + generate plumbing: pocket_map_path on a scaffold flows distance
data through MutablePosition to Mutation with no changes to scoring.
"""
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from opsin_pipeline.generate import generate_candidates
from opsin_pipeline.ingest import load_scaffolds
from opsin_pipeline.structure.ligands import identify_ligands
from opsin_pipeline.structure.pdb import parse_pdb
from opsin_pipeline.structure.pocket import (
    apply_offset_mapping,
    compute_pocket,
    pocket_map_from_dict,
    pocket_map_to_dict,
    read_pocket_map,
    write_pocket_map,
)


FIXTURE_DIR = Path(__file__).parent / "fixtures" / "structures"
SYNTHETIC = FIXTURE_DIR / "synthetic_mini.pdb"


class PocketMapJSONTests(unittest.TestCase):
    def test_roundtrip_via_json(self):
        atoms = parse_pdb(SYNTHETIC)
        pocket = compute_pocket(
            atoms,
            identify_ligands(atoms),
            scaffold_name="synthetic",
            pdb_path=str(SYNTHETIC),
        )

        serialized = pocket_map_to_dict(pocket)
        restored = pocket_map_from_dict(serialized)

        self.assertEqual(restored.scaffold_name, pocket.scaffold_name)
        self.assertEqual(
            [r.pdb_resnum for r in restored.pocket_residues],
            [r.pdb_resnum for r in pocket.pocket_residues],
        )

    def test_write_and_read_roundtrip(self):
        atoms = parse_pdb(SYNTHETIC)
        pocket = compute_pocket(
            atoms,
            identify_ligands(atoms),
            scaffold_name="synthetic",
            pdb_path=str(SYNTHETIC),
        )

        with tempfile.TemporaryDirectory() as tmp:
            path = write_pocket_map(pocket, Path(tmp) / "pocket.json")
            restored = read_pocket_map(path)

        self.assertEqual(restored.cutoff_A, pocket.cutoff_A)
        self.assertEqual(restored.thresholds_A, pocket.thresholds_A)


class OffsetMappingTests(unittest.TestCase):
    def _pocket(self):
        atoms = parse_pdb(SYNTHETIC)
        return compute_pocket(
            atoms,
            identify_ligands(atoms),
            scaffold_name="synthetic",
            pdb_path=str(SYNTHETIC),
        )

    def test_offset_sets_seq_index(self):
        mapped = apply_offset_mapping(self._pocket(), offset=10)

        for r in mapped.pocket_residues:
            self.assertEqual(r.seq_index, r.pdb_resnum + 10)

    def test_aa_match_verification_flags_mismatch(self):
        pocket = self._pocket()
        # make a scaffold sequence where position 11 is M but residue 1 in fixture is ALA (A)
        scaffold_seq = "X" * 20
        scaffold_seq = scaffold_seq[:10] + "M" + scaffold_seq[11:]

        mapped = apply_offset_mapping(
            pocket, offset=10, scaffold_sequence=scaffold_seq, strict=False
        )
        ala_row = next(r for r in mapped.pocket_residues if r.pdb_resnum == 1)

        self.assertIsNotNone(ala_row.mapping_note)
        self.assertIn("A", ala_row.mapping_note)
        self.assertIn("M", ala_row.mapping_note)

    def test_strict_mismatch_raises(self):
        pocket = self._pocket()
        scaffold_seq = "M" * 20

        with self.assertRaises(ValueError):
            apply_offset_mapping(
                pocket, offset=10, scaffold_sequence=scaffold_seq, strict=True
            )

    def test_out_of_bounds_seq_index_is_flagged_not_silent(self):
        """Regression: apply_offset_mapping used to silently accept seq_index outside
        the scaffold sequence. Now it records a mapping_note so downstream reviewers
        notice (and with strict=True it becomes a hard error)."""
        pocket = self._pocket()
        # fixture has pdb_resnums 1..6; offset 100 pushes them to 101..106.
        # Against a 10-residue scaffold, every residue is out of bounds.
        scaffold_seq = "MKTAYIAKQR"

        mapped = apply_offset_mapping(
            pocket, offset=100, scaffold_sequence=scaffold_seq, strict=False
        )

        self.assertTrue(
            all(r.mapping_note is not None for r in mapped.pocket_residues),
            "every out-of-bounds residue must carry a mapping_note",
        )
        for r in mapped.pocket_residues:
            self.assertIn("outside scaffold length", r.mapping_note)

    def test_out_of_bounds_strict_raises(self):
        pocket = self._pocket()
        scaffold_seq = "MKTAYIAKQR"

        with self.assertRaises(ValueError) as ctx:
            apply_offset_mapping(
                pocket, offset=100, scaffold_sequence=scaffold_seq, strict=True
            )
        self.assertIn("outside scaffold length", str(ctx.exception))

    def test_without_scaffold_sequence_no_bounds_check(self):
        """Back-compat: callers who don't pass scaffold_sequence still get raw
        seq_index = pdb_resnum + offset with no bounds checking or AA verification."""
        pocket = self._pocket()

        mapped = apply_offset_mapping(pocket, offset=1000)

        for r in mapped.pocket_residues:
            self.assertIsNone(r.mapping_note)
            self.assertEqual(r.seq_index, r.pdb_resnum + 1000)


class ScaffoldIngestMergeTests(unittest.TestCase):
    def _write_pocket_map(self, directory: Path) -> Path:
        # build the pocket map for our synthetic, apply offset=0 so seq_index == pdb_resnum
        atoms = parse_pdb(SYNTHETIC)
        pocket = compute_pocket(
            atoms,
            identify_ligands(atoms),
            scaffold_name="demo",
            pdb_path=str(SYNTHETIC),
        )
        mapped = apply_offset_mapping(pocket, offset=0)
        return write_pocket_map(mapped, directory / "pocket.json")

    def test_load_scaffolds_merges_pocket_distances(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            pocket_path = self._write_pocket_map(tmp_path)
            # synthetic scaffold whose seq positions align 1:1 with the fixture PDB's resnums
            scaffold_payload = {
                "scaffolds": [
                    {
                        "name": "demo",
                        "family": "Synthetic",
                        "sequence": "AFYKWVSG",  # positions 1-8 match fixture residues
                        "pocket_map_path": pocket_path.name,
                        "mutable_positions": [
                            {"position": 1, "allowed": ["F"], "reason": "pocket test"},
                            {"position": 4, "allowed": ["F"], "reason": "medium test"},
                            {"position": 8, "allowed": ["F"], "reason": "far away"},
                        ],
                    }
                ]
            }
            scaffolds_path = tmp_path / "scaffolds.json"
            scaffolds_path.write_text(json.dumps(scaffold_payload), encoding="utf-8")

            scaffolds = load_scaffolds(scaffolds_path)

        positions_by_index = {p.position: p for p in scaffolds[0].mutable_positions}
        # position 1 -> strong band, distance ~2.0
        self.assertAlmostEqual(positions_by_index[1].distance_to_retinal, 2.0, places=3)
        # position 4 -> medium band, distance ~4.5
        self.assertAlmostEqual(positions_by_index[4].distance_to_retinal, 4.5, places=3)
        # position 8 -> outside cutoff (10 A > 6 A) -> no pocket entry -> distance stays None
        self.assertIsNone(positions_by_index[8].distance_to_retinal)

    def test_generated_mutations_carry_distance(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            pocket_path = self._write_pocket_map(tmp_path)
            payload = {
                "scaffolds": [
                    {
                        "name": "demo",
                        "family": "Synthetic",
                        "sequence": "AFYKWVSG",
                        "pocket_map_path": pocket_path.name,
                        "mutable_positions": [
                            {"position": 1, "allowed": ["F"], "reason": "pocket"},
                            {"position": 8, "allowed": ["F"], "reason": "far"},
                        ],
                    }
                ]
            }
            scaffolds_path = tmp_path / "scaffolds.json"
            scaffolds_path.write_text(json.dumps(payload), encoding="utf-8")
            scaffolds = load_scaffolds(scaffolds_path)

        candidates, _stats = generate_candidates(scaffolds, max_mutations=1)

        by_position = {m.position: m for c in candidates for m in c.mutations}
        self.assertAlmostEqual(by_position[1].distance_to_retinal, 2.0, places=3)
        self.assertIsNone(by_position[8].distance_to_retinal)


class CandidateCSVDistanceColumnTests(unittest.TestCase):
    def _write_pocket_map(self, directory: Path) -> Path:
        atoms = parse_pdb(SYNTHETIC)
        pocket = compute_pocket(
            atoms,
            identify_ligands(atoms),
            scaffold_name="demo",
            pdb_path=str(SYNTHETIC),
        )
        mapped = apply_offset_mapping(pocket, offset=0)
        return write_pocket_map(mapped, directory / "pocket.json")

    def test_csv_has_min_distance_column_populated_when_pocket_data_present(self):
        import csv as _csv

        from opsin_pipeline.report import write_candidate_csv
        from opsin_pipeline.score import rank_candidates

        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            pocket_path = self._write_pocket_map(tmp_path)
            payload = {
                "scaffolds": [
                    {
                        "name": "demo",
                        "family": "Synthetic",
                        "sequence": "AFYKWVSG",
                        "pocket_map_path": pocket_path.name,
                        "mutable_positions": [
                            {"position": 1, "allowed": ["F"], "reason": "pocket"},
                            {"position": 8, "allowed": ["F"], "reason": "far"},
                        ],
                    }
                ]
            }
            scaffolds_path = tmp_path / "scaffolds.json"
            scaffolds_path.write_text(json.dumps(payload), encoding="utf-8")
            scaffolds = load_scaffolds(scaffolds_path)
            candidates, _ = generate_candidates(scaffolds, max_mutations=1)
            ranked = rank_candidates(candidates, scaffolds)

            csv_path = write_candidate_csv(ranked, scaffolds, tmp_path / "ranked.csv")
            with csv_path.open(encoding="utf-8") as handle:
                rows = list(_csv.DictReader(handle))

        self.assertIn("min_distance_to_retinal_A", rows[0].keys())
        by_id = {row["candidate_id"]: row for row in rows}
        # position 1 has pocket distance 2.0 A
        self.assertEqual(by_id["demo_p1AtoF"]["min_distance_to_retinal_A"], "2.00")
        # position 8 is outside the cutoff -> no distance -> blank
        self.assertEqual(by_id["demo_p8GtoF"]["min_distance_to_retinal_A"], "")


if __name__ == "__main__":
    unittest.main()
