from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))
SCRIPT_PATH = SCRIPTS_DIR / "25_stage4_rosetta_interface_scoring.py"
SPEC = importlib.util.spec_from_file_location("stage4_scoring", SCRIPT_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Could not import {SCRIPT_PATH}")
stage4 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage4)


STAGE3D1_TEST_FIELDS = [
    "stage3d1_design_id",
    "repack_status",
    "backbone_geometry_status",
    "post_target_site_recovery_status",
    "post_hotspot_recovery_status",
    "post_macrocycle_geometry_status",
    "post_clash_status",
    "peptide_sequence",
    "repacked_peptide_sequence",
    "target_sequence",
    "repacked_target_sequence",
    "stage3d1_failure_reasons",
    "pass_stage3d1_qc",
]


def _stage3d1_pass_row(design_id: str) -> dict[str, str]:
    return {
        "stage3d1_design_id": design_id,
        "repack_status": "success",
        "backbone_geometry_status": "pass_fixed_backbone",
        "post_target_site_recovery_status": "site_contact_pass",
        "post_hotspot_recovery_status": "hotspot_contact_pass",
        "post_macrocycle_geometry_status": "pass_head_to_tail_macrocycle",
        "post_clash_status": "pass_no_severe_clash",
        "peptide_sequence": "AGK",
        "repacked_peptide_sequence": "AGK",
        "target_sequence": "AAAA",
        "repacked_target_sequence": "AAAA",
        "stage3d1_failure_reasons": "",
        "pass_stage3d1_qc": "true",
    }


def _priority_row(
    design_id: str,
    backbone_id: str,
    ddg: float,
) -> dict[str, object]:
    return {
        "stage4_design_id": design_id,
        "global_backbone_id": backbone_id,
        "macrocycle_geometry_status": "pass_head_to_tail_macrocycle",
        "target_site_recovery_status": "site_contact_pass",
        "hotspot_recovery_status": "hotspot_contact_pass",
        "clash_status": "pass_no_severe_clash",
        "detached_or_collapsed_flag": "pass_basic_pose_geometry",
        "ddg_proxy_no_repack": ddg,
        "ddg_proxy_status": "favorable_or_neutral_interface_proxy",
        "complex_rosetta_total_score": ddg,
        "num_hotspot_contacts": 2,
        "num_target_site_contacts": 3,
        "sequence_liability_count": 0,
        "sequence_liability_notes": "none",
    }


class Stage4ProvenanceTests(unittest.TestCase):
    def test_stage4_output_field_lists_have_no_duplicates(self) -> None:
        self.assertEqual(
            len(stage4.STAGE4_FIELDS),
            len(set(stage4.STAGE4_FIELDS)),
        )
        self.assertEqual(
            len(stage4.TOP_VALIDATION_FIELDS),
            len(set(stage4.TOP_VALIDATION_FIELDS)),
        )

    def test_stage3d1_pass_table_must_be_exact_subset(self) -> None:
        first = _stage3d1_pass_row("design_a")
        second = _stage3d1_pass_row("design_b")
        validated = stage4._validate_stage3d1_pass_subset(
            stage3d1_rows=[first, second],
            stage3d1_pass_rows=[first, second],
            stage3d1_fields=STAGE3D1_TEST_FIELDS,
        )
        self.assertEqual(
            [row["stage3d1_design_id"] for row in validated],
            ["design_a", "design_b"],
        )

        with self.assertRaisesRegex(RuntimeError, "not the exact pass subset"):
            stage4._validate_stage3d1_pass_subset(
                stage3d1_rows=[first, second],
                stage3d1_pass_rows=[first],
                stage3d1_fields=STAGE3D1_TEST_FIELDS,
            )

    def test_stage3d1_pass_flag_must_match_hard_qc(self) -> None:
        bad = _stage3d1_pass_row("design_bad")
        bad["post_clash_status"] = "fail_severe_clash"
        with self.assertRaisesRegex(RuntimeError, "pass flag disagrees"):
            stage4._validate_stage3d1_pass_subset(
                stage3d1_rows=[bad],
                stage3d1_pass_rows=[bad],
                stage3d1_fields=STAGE3D1_TEST_FIELDS,
            )

    def test_top_validation_selection_is_backbone_diverse(self) -> None:
        rows = [
            _priority_row("a1", "backbone_a", -20.0),
            _priority_row("a2", "backbone_a", -19.0),
            _priority_row("b1", "backbone_b", -18.0),
            _priority_row("c1", "backbone_c", -17.0),
        ]
        stage4._assign_priority_ranks(rows)
        selected = stage4._select_validation_source_rows(
            rows,
            top_validation_count=3,
            max_per_global_backbone=1,
        )
        self.assertEqual(
            [row["stage4_design_id"] for row in selected],
            ["a1", "b1", "c1"],
        )
        self.assertEqual(
            [row["stage4_validation_selection_rank"] for row in selected],
            [1, 2, 3],
        )

    def test_ddg_normalization_handles_zero_denominators(self) -> None:
        metrics = stage4._normalized_ddg_metrics(
            ddg_proxy=-12.0,
            peptide_length=12,
            target_contacts=6,
            site_contacts=0,
        )
        self.assertEqual(metrics["ddg_proxy_per_peptide_residue"], -1.0)
        self.assertEqual(metrics["ddg_proxy_per_target_contact"], -2.0)
        self.assertEqual(metrics["ddg_proxy_per_site_contact"], "")


if __name__ == "__main__":
    unittest.main()
