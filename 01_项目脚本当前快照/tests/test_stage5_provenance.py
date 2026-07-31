from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import stage5_contract as contract  # noqa: E402


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


stage5a_prep = _load("stage5a_prep", SCRIPTS / "26_prepare_afcycdesign_jobs.py")
stage5a_collect = _load("stage5a_collect", SCRIPTS / "27_collect_afcycdesign_validation.py")
stage5b_prep = _load("stage5b_prep", SCRIPTS / "30_prepare_stage5b_target_conditioned_jobs.py")
stage5b_collect = _load("stage5b_collect", SCRIPTS / "31_collect_stage5b_validation.py")
stage5a_runner = _load(
    "stage5a_runner",
    SCRIPTS / "external" / "run_afcycdesign_independent_recovery.py",
)
stage5b_runner = _load(
    "stage5b_runner",
    SCRIPTS / "external" / "run_afcycdesign_target_conditioned_recovery.py",
)


def _identity_row(candidate_id: str = "candidate_1") -> dict[str, str]:
    row = {
        "stage5_candidate_id": candidate_id,
        "stage5B_candidate_id": candidate_id,
        "peptide_sequence_hash": "seqhash1",
        "protocol_hash": "protocol1",
    }
    for field in contract.STAGE4_IDENTITY_FIELDS:
        row[field] = f"value_{field}"
    return row


class Stage5ProvenanceTests(unittest.TestCase):
    def test_output_field_lists_have_no_duplicates(self) -> None:
        for module, field_name in (
            (stage5a_prep, "MANIFEST_FIELDS"),
            (stage5a_prep, "JOB_FIELDS"),
            (stage5a_collect, "MODEL_FIELDS"),
            (stage5a_collect, "CANDIDATE_FIELDS"),
            (stage5b_prep, "MANIFEST_FIELDS"),
            (stage5b_prep, "JOB_FIELDS"),
            (stage5b_collect, "MODEL_FIELDS"),
            (stage5b_collect, "CANDIDATE_FIELDS"),
        ):
            fields = getattr(module, field_name)
            self.assertEqual(len(fields), len(set(fields)), f"{module.__name__}.{field_name}")

    def test_stage5_job_identity_must_match_candidate(self) -> None:
        candidate = _identity_row()
        job = dict(candidate)
        contract.validate_stage5_identity_link(
            candidate=candidate,
            job=job,
            candidate_id_field="stage5_candidate_id",
            job_candidate_id_field="stage5_candidate_id",
            label="test job",
        )
        job["global_backbone_id"] = "wrong_global_id"
        with self.assertRaisesRegex(RuntimeError, "global_backbone_id"):
            contract.validate_stage5_identity_link(
                candidate=candidate,
                job=job,
                candidate_id_field="stage5_candidate_id",
                job_candidate_id_field="stage5_candidate_id",
                label="test job",
            )

    def test_stage4_hard_gate_rejects_clash(self) -> None:
        row = {
            "pass_stage4_qc": "true",
            "pyrosetta_score_status": "success",
            "target_site_recovery_status": "site_contact_pass",
            "hotspot_recovery_status": "hotspot_contact_pass",
            "macrocycle_geometry_status": "pass_head_to_tail_macrocycle",
            "clash_status": "pass_no_severe_clash",
            "detached_or_collapsed_flag": "pass_basic_pose_geometry",
            "stage4_failure_reasons": "",
            "target_chain": "A",
            "peptide_chain": "B",
            "num_target_site_contacts": "2",
            "num_hotspot_contacts": "1",
        }
        contract._validate_stage4_hard_gate(row, "test row")
        row["clash_status"] = "fail_severe_clash"
        with self.assertRaisesRegex(RuntimeError, "clash_status"):
            contract._validate_stage4_hard_gate(row, "test row")

    def test_runners_cache_on_stage4_identity(self) -> None:
        required = {
            "global_backbone_id",
            "backbone_family_id",
            "source_stage4_design_id",
            "source_stage4_scored_pdb_sha256",
            "stage4_run_id",
            "stage4_protocol_identity_sha256",
            "stage4_scores_csv_sha256",
            "stage4_top_candidates_csv_sha256",
        }
        self.assertTrue(required.issubset(set(stage5a_runner.CACHE_IDENTITY_FIELDS)))
        self.assertTrue(required.issubset(set(stage5b_runner.CACHE_IDENTITY_FIELDS)))

    def test_preparation_protocol_versions_are_active_route_versions(self) -> None:
        self.assertEqual(
            stage5a_prep.PROTOCOL_VERSION,
            "stage5A_v3_active_route_single_sequence_mlm015",
        )
        self.assertEqual(
            stage5b_prep.PROTOCOL_VERSION,
            "stage5B_v2_active_route_target_only_template_single_sequence_no_mlm_dropout",
        )

    def test_collectors_require_complete_authoritative_outputs(self) -> None:
        stage5a_text = (SCRIPTS / "27_collect_afcycdesign_validation.py").read_text(encoding="utf-8")
        stage5b_text = (SCRIPTS / "31_collect_stage5b_validation.py").read_text(encoding="utf-8")
        self.assertIn("Missing authoritative Stage 5A output directory", stage5a_text)
        self.assertIn("has {len(metrics)} model records", stage5a_text)
        self.assertIn("Missing authoritative Stage 5B output", stage5b_text)
        self.assertIn("has {len(metrics)} model records", stage5b_text)


if __name__ == "__main__":
    unittest.main()
