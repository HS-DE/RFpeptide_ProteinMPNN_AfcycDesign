from __future__ import annotations

import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))


def _load_script(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS_DIR / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {filename}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


stage24 = _load_script(
    "stage24_under_test",
    "24_stage3d1_sidechain_repack.py",
)


def _atom(
    serial: int,
    atom: str,
    chain: str,
    resi: int,
    x: float,
    resname: str = "GLY",
) -> str:
    return (
        f"ATOM  {serial:5d} {atom:>4s} {resname:>3s} {chain}{resi:4d}    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C"
    )


def _write_backbone(
    path: Path,
    *,
    ca_shift: float = 0.0,
    resname: str = "GLY",
) -> None:
    lines = []
    serial = 1
    for resi, base in [(1, 0.0), (2, 4.0)]:
        for atom, offset in [("N", 0.0), ("CA", 1.0), ("C", 2.0), ("O", 3.0)]:
            shift = ca_shift if atom == "CA" and resi == 2 else 0.0
            lines.append(
                _atom(
                    serial,
                    atom,
                    "B",
                    resi,
                    base + offset + shift,
                    resname=resname,
                )
            )
            serial += 1
    path.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")


def _stage3c_row(
    backbone_id: str,
    sequence_id: str,
    sequence: str,
) -> dict[str, str]:
    return {
        "global_backbone_id": backbone_id,
        "backbone_id": backbone_id,
        "sequence_design_id": sequence_id,
        "peptide_sequence": sequence,
        "file_status": "pass",
        "parse_status": "pass",
        "sequence_status": "pass_sequence",
        "target_site_recovery_status": "site_contact_pass",
        "hotspot_recovery_status": "hotspot_contact_pass",
        "macrocycle_geometry_status": "pass_head_to_tail_macrocycle",
    }


class Stage3D1SelectionTests(unittest.TestCase):
    def test_sequence_deduplication_is_scoped_to_global_backbone(self) -> None:
        rows = [
            _stage3c_row("run_a__bb_1", "seq_1", "AAAA"),
            _stage3c_row("run_a__bb_1", "seq_2", "AAAA"),
            _stage3c_row("run_b__bb_1", "seq_3", "AAAA"),
        ]
        selected, skipped = stage24._selected_stage3c_rows(
            rows,
            set(),
            include_duplicate_sequences=False,
        )

        self.assertEqual(
            [row["sequence_design_id"] for row in selected],
            ["seq_1", "seq_3"],
        )
        self.assertEqual(len(skipped), 1)
        self.assertEqual(
            skipped[0]["skip_reason"],
            "duplicate_sequence_within_global_backbone",
        )

    def test_clash_failure_is_eligible_for_repack_when_other_gates_pass(self) -> None:
        row = _stage3c_row("run_a__bb_1", "seq_1", "AAAA")
        row["clash_status"] = "fail_severe_clash"
        row["pass_stage3c_qc"] = "false"
        selected, skipped = stage24._selected_stage3c_rows(
            [row],
            set(),
            include_duplicate_sequences=False,
        )

        self.assertEqual(len(selected), 1)
        self.assertFalse(skipped)


class Stage3D1BackboneTests(unittest.TestCase):
    def test_fixed_backbone_metrics_pass_for_identical_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            before = root / "before.pdb"
            after = root / "after.pdb"
            _write_backbone(before)
            _write_backbone(after)
            metrics = stage24._backbone_displacement_metrics(
                input_pdb=before,
                output_pdb=after,
                chain_id="B",
                chain_label="peptide",
                max_displacement_tolerance=0.002,
            )

        self.assertEqual(metrics["peptide_backbone_rmsd_A"], 0.0)
        self.assertEqual(
            metrics["peptide_backbone_max_displacement_A"],
            0.0,
        )
        self.assertTrue(metrics["peptide_backbone_fixed"])

    def test_fixed_backbone_metrics_detect_coordinate_movement(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            before = root / "before.pdb"
            after = root / "after.pdb"
            _write_backbone(before)
            _write_backbone(after, ca_shift=0.01)
            metrics = stage24._backbone_displacement_metrics(
                input_pdb=before,
                output_pdb=after,
                chain_id="B",
                chain_label="peptide",
                max_displacement_tolerance=0.002,
            )

        self.assertAlmostEqual(
            metrics["peptide_backbone_max_displacement_A"],
            0.01,
            places=6,
        )
        self.assertFalse(metrics["peptide_backbone_fixed"])

    def test_fixed_backbone_metrics_reject_sequence_change(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            before = root / "before.pdb"
            after = root / "after.pdb"
            _write_backbone(before)
            _write_backbone(after, resname="ALA")
            with self.assertRaisesRegex(RuntimeError, "sequence changed"):
                stage24._backbone_displacement_metrics(
                    input_pdb=before,
                    output_pdb=after,
                    chain_id="B",
                    chain_label="peptide",
                    max_displacement_tolerance=0.002,
                )


class Stage3D1TableContractTests(unittest.TestCase):
    def _fixture(self, root: Path):
        stage2_selection_csv = root / "selection.csv"
        stage3_jobs_csv = root / "jobs.csv"
        input_pdb = root / "input.pdb"
        source_pdb = root / "source.pdb"
        output_pdb = root / "global_a_dldesign_0.pdb"
        for path in [
            stage2_selection_csv,
            stage3_jobs_csv,
            input_pdb,
            source_pdb,
            output_pdb,
        ]:
            path.write_text(f"{path.name}\n", encoding="utf-8")

        backbone_id = "rfp_run_a__RFpep_Site_2_0001"
        sequence_id = output_pdb.stem
        selection_row = {
            "global_backbone_id": backbone_id,
            "source_local_design_id": "RFpep_Site_2_0001",
            "source_batch_label": "batch_a",
            "backbone_family_id": "family_1",
            "site_label": "RFpep_Site_2",
            "site_id": "site_2",
            "target_chain": "A",
            "peptide_chain": "B",
        }
        job_row = {
            **selection_row,
            "stage3_job_id": "job_a",
            "source_backbone_pdb": str(source_pdb),
            "source_backbone_pdb_sha256": stage24._sha256_file(source_pdb),
            "source_run_id": "rfp_run_a",
            "source_batch_id": "batch_a",
            "source_route_manifest": str(root / "route_manifest.json"),
            "source_route_manifest_sha256": "a" * 64,
        }
        row = {
            **selection_row,
            "backbone_id": backbone_id,
            "sequence_design_id": sequence_id,
            "relaxed_pdb": str(output_pdb),
            "relaxed_pdb_sha256": stage24._sha256_file(output_pdb),
            "stage3_job_id": "job_a",
            "run_group_id": "stage3_1bp_" + "b" * 12,
            "protocol_identity_sha256": "b" * 64,
            "stage3_mode": "proteinmpnn_only",
            "stage3_jobs_csv": str(stage3_jobs_csv),
            "stage3_jobs_csv_sha256": stage24._sha256_file(stage3_jobs_csv),
            "stage2_selection_csv": str(stage2_selection_csv),
            "stage2_selection_csv_sha256": stage24._sha256_file(
                stage2_selection_csv
            ),
            "stage3_input_pdb": str(input_pdb),
            "stage3_input_pdb_sha256": stage24._sha256_file(input_pdb),
            "source_backbone_pdb": str(source_pdb),
            "source_backbone_pdb_sha256": stage24._sha256_file(source_pdb),
            "expected_stage3b_outputs_for_backbone": "1",
            "observed_stage3b_outputs_for_backbone": "1",
            "source_run_id": "rfp_run_a",
            "source_batch_id": "batch_a",
            "source_route_manifest": str(root / "route_manifest.json"),
            "source_route_manifest_sha256": "a" * 64,
        }
        expected = {
            backbone_id: {
                output_pdb.resolve(): (input_pdb.resolve(), "input")
            }
        }
        return (
            row,
            selection_row,
            job_row,
            expected,
            stage2_selection_csv,
            stage3_jobs_csv,
        )

    def test_stage3c_table_rejects_local_backbone_id_substitution(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (
                row,
                selection_row,
                job_row,
                expected,
                selection_csv,
                jobs_csv,
            ) = self._fixture(root)
            qchelper = mock.Mock()
            stage24._validate_stage3c_table(
                qchelper=qchelper,
                stage3c_rows=[row],
                stage2_selection_lookup={
                    row["global_backbone_id"]: selection_row
                },
                stage3_job_lookup={row["global_backbone_id"]: job_row},
                expected_by_backbone=expected,
                stage3_route_provenance={},
                aggregate_sources={},
                stage2_selection_csv=selection_csv,
                stage2_selection_csv_sha256=stage24._sha256_file(
                    selection_csv
                ),
                stage3_jobs_csv=jobs_csv,
                stage3_jobs_csv_sha256=stage24._sha256_file(jobs_csv),
                run_group_id=row["run_group_id"],
                protocol_identity_sha256=row[
                    "protocol_identity_sha256"
                ],
                stage3_mode="proteinmpnn_only",
            )
            bad_row = dict(row)
            bad_row["backbone_id"] = "RFpep_Site_2_0001"
            with self.assertRaisesRegex(RuntimeError, "backbone_id"):
                stage24._validate_stage3c_table(
                    qchelper=qchelper,
                    stage3c_rows=[bad_row],
                    stage2_selection_lookup={
                        row["global_backbone_id"]: selection_row
                    },
                    stage3_job_lookup={
                        row["global_backbone_id"]: job_row
                    },
                    expected_by_backbone=expected,
                    stage3_route_provenance={},
                    aggregate_sources={},
                    stage2_selection_csv=selection_csv,
                    stage2_selection_csv_sha256=stage24._sha256_file(
                        selection_csv
                    ),
                    stage3_jobs_csv=jobs_csv,
                    stage3_jobs_csv_sha256=stage24._sha256_file(jobs_csv),
                    run_group_id=row["run_group_id"],
                    protocol_identity_sha256=row[
                        "protocol_identity_sha256"
                    ],
                    stage3_mode="proteinmpnn_only",
                )


if __name__ == "__main__":
    unittest.main()
