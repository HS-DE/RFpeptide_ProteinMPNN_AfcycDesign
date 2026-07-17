from __future__ import annotations

import builtins
import copy
import importlib.util
import json
import logging
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


SNAPSHOT_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = SNAPSHOT_ROOT.parent
SCRIPTS = SNAPSHOT_ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))

import common  # noqa: E402


def _load_stage26():
    path = SCRIPTS / "26_prepare_afcycdesign_jobs.py"
    spec = importlib.util.spec_from_file_location("stage26_route_tests", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load Stage 26")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_stage20():
    path = SCRIPTS / "20_make_rfpeptides_article_jobs.py"
    spec = importlib.util.spec_from_file_location("stage20_route_tests", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load Stage 20")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ACTIVE_CONFIG = REPO_ROOT / "05_配置快照" / "config" / "rfpeptides_head_to_tail.yaml"
LEGACY_CONFIG = REPO_ROOT / "05_配置快照" / "config" / "project_legacy_disulfide.yaml"
PROJECT_CONFIG = REPO_ROOT / "05_配置快照" / "config" / "project.yaml"


class RouteManifestTests(unittest.TestCase):
    def _payload(self, root: Path) -> dict:
        target = root / "target.pdb"
        mapping = root / "mapping.csv"
        target.write_text("ATOM\nEND\n", encoding="utf-8")
        mapping.write_text("rfpeptides_chain,rfpeptides_residue_number\nA,1\n", encoding="utf-8")
        _, config_sha, effective_sha = common.load_active_route_config(ACTIVE_CONFIG)
        hotspots = ["A82", "A84", "A85", "A86"]
        return {
            "batch_id": "N1_v3_fixture",
            "site_labels": ["RFpep_Site_2"],
            "protocol_peptide_length_min": 12,
            "protocol_peptide_length_max": 24,
            "run_peptide_length_min": 17,
            "run_peptide_length_max": 17,
            "num_designs_requested": 1,
            "project_config": str(ACTIVE_CONFIG),
            "project_config_sha256": config_sha,
            "effective_project_config_sha256": effective_sha,
            "stage0_sites": [
                {
                    "site_label": "RFpep_Site_2",
                    "target_pdb": str(target),
                    "target_pdb_sha256": common.sha256_file(target),
                    "mapping_csv": str(mapping),
                    "mapping_csv_sha256": common.sha256_file(mapping),
                    "normalized_hotspots": hotspots,
                    "normalized_hotspots_sha256": common.canonical_json_sha256(hotspots),
                }
            ],
        }

    def test_n1_manifest_round_trip_and_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path, manifest, digest = common.write_route_manifest(root, self._payload(root))
            loaded_path, loaded, loaded_digest = common.load_route_manifest(root)
            self.assertEqual(path, loaded_path)
            self.assertEqual(manifest["run_id"], loaded["run_id"])
            self.assertEqual(digest, loaded_digest)
            provenance = common.route_provenance_fields(path, manifest, digest)
            common.validate_row_route_provenance(provenance, provenance, "fixture")

            mirror_path, mirror_manifest, mirror_digest = common.write_route_manifest(root / "mirror", manifest)
            mirror_provenance = common.route_provenance_fields(mirror_path, mirror_manifest, mirror_digest)
            common.validate_row_route_provenance(provenance, mirror_provenance, "mirrored fixture")

    def test_run_identity_uses_content_hashes_not_absolute_input_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            first = self._payload(root)
            second = copy.deepcopy(first)
            second["project_config"] = str(root / "alternate" / "config.yaml")
            second["stage0_sites"][0]["target_pdb"] = str(root / "alternate" / "target.pdb")
            second["stage0_sites"][0]["mapping_csv"] = str(root / "alternate" / "mapping.csv")
            self.assertEqual(
                common.finalize_route_manifest(first)["run_id"],
                common.finalize_route_manifest(second)["run_id"],
            )

    def test_manifest_rejects_tampered_file_and_bad_length(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = self._payload(root)
            common.write_route_manifest(root, payload)
            Path(payload["stage0_sites"][0]["target_pdb"]).write_text("changed\n", encoding="utf-8")
            with self.assertRaises((RuntimeError, FileNotFoundError)):
                common.load_route_manifest(root)

        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            payload = self._payload(root)
            payload["run_peptide_length_max"] = 25
            common.write_route_manifest(root, payload)
            with self.assertRaises(RuntimeError):
                common.load_route_manifest(root)

    def test_manifest_rejects_route_identity_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            _, manifest, _ = common.write_route_manifest(root, self._payload(root))
            manifest["cyclization"] = "Cys-Cys disulfide"
            with self.assertRaises(RuntimeError):
                common.validate_route_manifest(root / "route_manifest.json", manifest)

    def test_manifest_missing_or_incomplete_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            with self.assertRaises((RuntimeError, FileNotFoundError)):
                common.load_route_manifest(root)

            path, _, _ = common.write_route_manifest(root, self._payload(root))
            payload = json.loads(path.read_text(encoding="utf-8"))
            payload.pop("hotspot_mapping_version")
            path.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaises(RuntimeError):
                common.load_route_manifest(root)

    def test_stage20_dry_run_writes_locked_manifest_and_job_provenance(self) -> None:
        stage20 = _load_stage20()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stage0_root = root / "stage0"
            inputs = stage0_root / "00_target_inputs"
            inputs.mkdir(parents=True)
            target = inputs / "RFpep_Site_2_target.pdb"
            mapping = inputs / "RFpep_Site_2_crop_renumbering_mapping.csv"
            target.write_text("ATOM\nEND\n", encoding="utf-8")
            mapping.write_text(
                "rfpeptides_chain,rfpeptides_residue_number,is_selected_hotspot\n"
                "A,1,true\nA,2,true\nA,3,true\nA,4,true\n",
                encoding="utf-8",
            )
            summary = inputs / "FGA_rfpeptides_stage0_target_inputs_summary.csv"
            summary.write_text(
                "site_label,site_id,site_quality_tier,target_pdb,crop_renumbering_mapping_csv,"
                "rfpeptides_hotspots,rfpeptides_residue_range,hotspots_txt\n"
                f"RFpep_Site_2,site2,high,{target},{mapping},\"A1,A2,A3,A4\",A1-A4,hotspots.txt\n",
                encoding="utf-8",
            )
            runtime = root / "rfd_macro"
            for rel in [
                "rfdiffusion/inference/utils.py",
                "rfdiffusion/inference/model_runners.py",
                "scripts/run_inference.py",
                "rfdiffusion/util.py",
            ]:
                path = runtime / rel
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("# fixture\n", encoding="utf-8")
            output = root / "stage1"
            argv = [
                "20_make_rfpeptides_article_jobs.py",
                "--input-root", str(stage0_root),
                "--output-root", str(output),
                "--selected-sites", "RFpep_Site_2",
                "--batch-id", "fixture_batch",
                "--project-config", str(ACTIVE_CONFIG),
                "--rfpeptides-root", str(runtime),
                "--num-designs", "1",
                "--length-min", "17",
                "--length-max", "17",
            ]
            with mock.patch("sys.argv", argv), mock.patch.object(
                stage20, "_git_head", return_value="f" * 40
            ), mock.patch.object(stage20, "setup_logger") as logger_factory:
                logger_factory.return_value = mock.Mock()
                self.assertEqual(stage20.main(), 0)
            _, manifest, manifest_sha = common.load_route_manifest(output)
            self.assertEqual(manifest["run_peptide_length_min"], 17)
            rows = common.read_csv(output / "01_rfpeptides_jobs" / "FGA_rfpeptides_stage1_jobs.csv")
            self.assertEqual(len(rows), 1)
            expected = common.route_provenance_fields(output / "route_manifest.json", manifest, manifest_sha)
            common.validate_row_route_provenance(rows[0], expected, "Stage 20 dry-run job")


class ActiveConfigTests(unittest.TestCase):
    def test_active_config_is_strict_and_hashes_are_stable(self) -> None:
        first = common.load_active_route_config(ACTIVE_CONFIG)
        second = common.load_active_route_config(ACTIVE_CONFIG)
        self.assertEqual(first[1:], second[1:])
        self.assertEqual(first[0]["peptide_design"]["cyclization"], "head_to_tail_amide")

    def test_legacy_and_project_configs_are_rejected(self) -> None:
        with self.assertRaises(RuntimeError):
            common.load_active_route_config(LEGACY_CONFIG)
        with self.assertRaises(RuntimeError):
            common.load_active_route_config(PROJECT_CONFIG)

    def test_missing_pyyaml_is_a_hard_failure(self) -> None:
        real_import = builtins.__import__

        def blocked_import(name, *args, **kwargs):
            if name == "yaml":
                raise ModuleNotFoundError("blocked for test")
            return real_import(name, *args, **kwargs)

        with mock.patch("builtins.__import__", side_effect=blocked_import):
            with self.assertRaises(RuntimeError):
                common.load_active_route_config(ACTIVE_CONFIG)

    def test_missing_shared_config_field_is_rejected(self) -> None:
        import yaml

        config = yaml.safe_load(ACTIVE_CONFIG.read_text(encoding="utf-8"))
        del config["structures"]["cleaned_pdb_file"]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "incomplete.yaml"
            path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
            with self.assertRaises(RuntimeError):
                common.load_active_route_config(path)


class MultiSourceRouteTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.stage26 = _load_stage26()

    def _manifest(self, run_id: str) -> dict:
        return {
            "run_id": run_id,
            "route_name": common.ROUTE_NAME,
            "route_protocol_version": common.ROUTE_PROTOCOL_VERSION,
            "hotspot_mapping_version": common.HOTSPOT_MAPPING_VERSION,
            "cyclization": common.ROUTE_CYCLIZATION,
            "project_config_sha256": "a" * 64,
            "effective_project_config_sha256": "b" * 64,
            "site_labels": ["RFpep_Site_2"],
            "stage0_sites": [{"target_pdb_sha256": "c" * 64}],
        }

    def test_duplicate_run_and_incompatible_sources_fail(self) -> None:
        one = self._manifest("run_one")
        duplicate = copy.deepcopy(one)
        with self.assertRaises(RuntimeError):
            self.stage26._validate_source_route_set(
                [
                    {"manifest": one, "manifest_sha256": "1" * 64},
                    {"manifest": duplicate, "manifest_sha256": "2" * 64},
                ]
            )

        two = self._manifest("run_two")
        two["stage0_sites"] = [{"target_pdb_sha256": "d" * 64}]
        with self.assertRaises(RuntimeError):
            self.stage26._validate_source_route_set(
                [
                    {"manifest": one, "manifest_sha256": "1" * 64},
                    {"manifest": two, "manifest_sha256": "2" * 64},
                ]
            )

    @staticmethod
    def _atom(serial: int, resname: str, chain: str, resi: int, x: float) -> str:
        return (
            f"ATOM  {serial:5d}  CA  {resname:>3s} {chain}{resi:4d}    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C"
        )

    def _write_stage5_fixture_pdb(self, path: Path, peptide_sequence: str = "") -> None:
        aa1_to_aa3 = {
            "A": "ALA", "V": "VAL", "I": "ILE", "L": "LEU", "M": "MET", "F": "PHE",
            "W": "TRP", "Y": "TYR", "K": "LYS", "R": "ARG", "D": "ASP", "E": "GLU",
        }
        lines = [self._atom(idx, "ALA", "A", idx, float(idx)) for idx in range(1, 87)]
        serial = 87
        for idx, aa in enumerate(peptide_sequence, start=1):
            lines.append(self._atom(serial, aa1_to_aa3[aa], "B", idx, float(idx)))
            serial += 1
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines + ["END"]) + "\n", encoding="utf-8")

    def test_stage26_merges_two_source_runs_with_per_candidate_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stage0_root = root / "stage0"
            stage0_inputs = stage0_root / "00_target_inputs"
            target_pdb = stage0_inputs / "RFpep_Site_2_target.pdb"
            mapping_csv = stage0_inputs / "RFpep_Site_2_crop_renumbering_mapping.csv"
            self._write_stage5_fixture_pdb(target_pdb)
            mapping_csv.parent.mkdir(parents=True, exist_ok=True)
            mapping_csv.write_text(
                "rfpeptides_residue_number,rfpeptides_residue_name,is_target_site_residue,is_selected_hotspot\n"
                + "".join(
                    f"{idx},ALA,{'true' if idx >= 80 else 'false'},{'true' if idx in {82, 84, 85, 86} else 'false'}\n"
                    for idx in range(1, 87)
                ),
                encoding="utf-8",
            )
            _, config_sha, effective_sha = common.load_active_route_config(ACTIVE_CONFIG)
            hotspots = ["A82", "A84", "A85", "A86"]
            source_roots = []
            sequences = ["AVILMFWYKRDE", "DEKRWYFMILVA"]
            for idx, sequence in enumerate(sequences, start=1):
                source_root = root / f"source_{idx}"
                source_roots.append(source_root)
                manifest_path, manifest, manifest_sha = common.write_route_manifest(
                    source_root,
                    {
                        "batch_id": f"fixture_batch_{idx}",
                        "site_labels": ["RFpep_Site_2"],
                        "protocol_peptide_length_min": 12,
                        "protocol_peptide_length_max": 24,
                        "run_peptide_length_min": 12,
                        "run_peptide_length_max": 12,
                        "num_designs_requested": 1,
                        "project_config": str(ACTIVE_CONFIG),
                        "project_config_sha256": config_sha,
                        "effective_project_config_sha256": effective_sha,
                        "stage0_sites": [
                            {
                                "site_label": "RFpep_Site_2",
                                "target_pdb": str(target_pdb),
                                "target_pdb_sha256": common.sha256_file(target_pdb),
                                "mapping_csv": str(mapping_csv),
                                "mapping_csv_sha256": common.sha256_file(mapping_csv),
                                "normalized_hotspots": hotspots,
                                "normalized_hotspots_sha256": common.canonical_json_sha256(hotspots),
                            }
                        ],
                    },
                )
                route = common.route_provenance_fields(manifest_path, manifest, manifest_sha)
                design_pdb = root / f"design_{idx}.pdb"
                self._write_stage5_fixture_pdb(design_pdb, sequence)
                row = {
                    "stage4_design_id": f"stage4_{idx}",
                    "backbone_id": f"RFpep_Site_2_fixture_{idx}",
                    "site_label": "RFpep_Site_2",
                    "site_id": "site2",
                    "scored_pdb": str(design_pdb),
                    "peptide_sequence": sequence,
                    "peptide_length": len(sequence),
                    "pass_stage4_qc": "true",
                    "macrocycle_geometry_status": "pass_head_to_tail_macrocycle",
                    "target_site_recovery_status": "site_contact_pass",
                    "hotspot_recovery_status": "hotspot_contact_pass",
                    "clash_status": "pass_no_severe_clash",
                    "detached_or_collapsed_flag": "pass_basic_pose_geometry",
                    "ddg_proxy_no_repack": -float(idx),
                    "num_target_contacts": 10,
                    "num_target_site_contacts": 4,
                    "num_hotspot_contacts": 2,
                    "stage4_priority_rank": 1,
                    **route,
                }
                common.write_csv(
                    source_root / "06_rosetta_scoring" / "FGA_rfpeptides_stage4_rosetta_interface_scores_pass.csv",
                    [row],
                    list(row),
                )

            output_root = root / "stage5"
            logger = logging.getLogger(f"stage26_fixture_{id(root)}")
            logger.handlers.clear()
            logger.addHandler(logging.NullHandler())
            argv = [
                "26_prepare_afcycdesign_jobs.py",
                "--source-run-root", str(source_roots[0]),
                "--source-run-root", str(source_roots[1]),
                "--stage0-root", str(stage0_root),
                "--output-root", str(output_root),
                "--project-config", str(ACTIVE_CONFIG),
                "--candidate-count", "2",
                "--seeds-per-candidate", "1",
            ]
            with mock.patch("sys.argv", argv), mock.patch.object(
                self.stage26, "setup_logger", return_value=logger
            ), mock.patch.object(self.stage26, "append_run_header", return_value=None):
                self.assertEqual(self.stage26.main(), 0)

            _, merged_manifest, _ = common.load_route_manifest(output_root)
            self.assertEqual(len(merged_manifest["source_route_manifests"]), 2)
            rows = common.read_csv(
                output_root / "07_structure_validation" / "FGA_rfpeptides_stage5_candidate_manifest.csv"
            )
            self.assertEqual(len(rows), 2)
            self.assertEqual({row["source_batch_id"] for row in rows}, {"fixture_batch_1", "fixture_batch_2"})
            self.assertTrue(all(Path(row["source_route_manifest"]).is_file() for row in rows))
            self.assertEqual(len({row["source_route_manifest_sha256"] for row in rows}), 2)


if __name__ == "__main__":
    unittest.main()
