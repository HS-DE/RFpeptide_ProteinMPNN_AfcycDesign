from __future__ import annotations

import csv
import importlib.util
import logging
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


stage22 = _load_script("stage22_under_test", "22_prepare_proteinmpnn_jobs.py")
stage23 = _load_script("stage23_under_test", "23_collect_proteinmpnn_sequences.py")


def _atom(serial: int, chain: str, resi: int, x: float) -> str:
    return (
        f"ATOM  {serial:5d}  CA  GLY {chain}{resi:4d}    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C"
    )


def _write_complex(path: Path, target_chain: str = "A", peptide_chain: str = "B") -> None:
    path.write_text(
        "\n".join(
            [
                "HEADER    SYNTHETIC TEST",
                "TER",
                _atom(1, target_chain, 1, 0.0),
                "TER",
                _atom(2, peptide_chain, 1, 2.0),
                "TER",
                "END",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


class Stage3ProvenanceTests(unittest.TestCase):
    def _prepare_fixture(self, root: Path) -> tuple[Path, Path, list[str]]:
        stage2_root = root / "stage2"
        stage2_root.mkdir(parents=True, exist_ok=True)
        ids = ["RFpep_Site_2_0002", "RFpep_Site_2_0017", "RFpep_Site_2_0018"]
        rows = []
        for idx, design_id in enumerate(ids, start=1):
            pdb_path = root / f"{design_id}.pdb"
            _write_complex(pdb_path)
            rows.append(
                {
                    "design_id": design_id,
                    "site_label": "RFpep_Site_2",
                    "site_id": "site2",
                    "rf_pdb": str(pdb_path),
                    "peptide_chain": "B",
                    "target_chain": "A",
                    "pass_backbone_qc": "true",
                    "peptide_length": str(12 + idx),
                }
            )
        pass_csv = root / "stage2_pass.csv"
        with pass_csv.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

        tool_root = root / "dl_binder_design"
        runner = tool_root / "mpnn_fr" / "dl_interface_design.py"
        checkpoint = tool_root / "mpnn_fr" / "ProteinMPNN" / "vanilla_model_weights" / "v_48_020.pt"
        runner.parent.mkdir(parents=True, exist_ok=True)
        checkpoint.parent.mkdir(parents=True, exist_ok=True)
        runner.write_text("print('fixture')\n", encoding="utf-8")
        checkpoint.write_bytes(b"fixture checkpoint")
        return stage2_root, pass_csv, ids

    def _run_stage22(self, root: Path, output_name: str, temperature: float = 0.1) -> list[dict[str, str]]:
        stage2_root, pass_csv, ids = self._prepare_fixture(root)
        output_root = root / output_name
        logger = logging.getLogger(f"stage22_fixture_{id(root)}")
        logger.handlers.clear()
        logger.addHandler(logging.NullHandler())
        argv = [
            "22_prepare_proteinmpnn_jobs.py",
            "--stage2-root",
            str(stage2_root),
            "--output-root",
            str(output_root),
            "--stage2-pass-csv",
            str(pass_csv),
            "--selected-backbones",
            ",".join(ids),
            "--dl-binder-design-root",
            str(root / "dl_binder_design"),
            "--temperature",
            str(temperature),
        ]
        with mock.patch("sys.argv", argv), mock.patch.object(stage22, "setup_logger", return_value=logger), mock.patch.object(
            stage22, "append_run_header", return_value=None
        ):
            self.assertEqual(stage22.main(), 0)
        jobs_files = list((output_root / "04_proteinmpnn_inputs").glob("*_jobs.csv"))
        self.assertEqual(len(jobs_files), 1)
        with jobs_files[0].open("r", encoding="utf-8-sig", newline="") as handle:
            return list(csv.DictReader(handle))

    def test_multi_backbone_jobs_are_one_to_one_and_protocol_aware(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            rows = self._run_stage22(root, "output_a", temperature=0.1)
            same_rows = self._run_stage22(root, "output_b", temperature=0.1)
            changed_rows = self._run_stage22(root, "output_c", temperature=0.2)

        self.assertEqual(len(rows), 3)
        self.assertEqual({row["design_id"] for row in rows}, {"RFpep_Site_2_0002", "RFpep_Site_2_0017", "RFpep_Site_2_0018"})
        self.assertTrue(all("," not in row["design_id"] for row in rows))
        self.assertTrue(all(row["design_id"] == row["backbone_id"] for row in rows))
        self.assertTrue(all("0007" not in row["stage3_job_id"] for row in rows))
        self.assertEqual({row["run_group_id"] for row in rows}, {same_rows[0]["run_group_id"]})
        self.assertNotEqual(rows[0]["run_group_id"], changed_rows[0]["run_group_id"])
        self.assertTrue(all(Path(row["source_backbone_pdb"]).stem == row["design_id"] for row in rows))

    def test_stage23_rejects_duplicate_and_aggregate_rows(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "Duplicate Stage 3 job"):
            stage23._strict_lookup_rows(
                [{"design_id": "A"}, {"design_id": "A"}], "design_id", "Stage 3 job"
            )
        with self.assertRaisesRegex(RuntimeError, "aggregate"):
            stage23._strict_lookup_rows([{"design_id": "A,B"}], "design_id", "Stage 3 job")


class PdbNormalizationTests(unittest.TestCase):
    def test_target_first_is_normalized_to_peptide_first(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.pdb"
            output = Path(tmp) / "output.pdb"
            _write_complex(source)
            note = stage22._copy_or_reorder_pdb(source, output, peptide_chain="B", target_chain="A")
            lines = output.read_text(encoding="utf-8").splitlines()

        self.assertEqual(note, "normalized_peptide_chain_first")
        atom_lines = [line for line in lines if line.startswith("ATOM  ")]
        self.assertEqual([line[21] for line in atom_lines], ["B", "A"])
        self.assertEqual(sum(line == "TER" for line in lines), 2)
        self.assertEqual(sum(line == "END" for line in lines), 1)
        self.assertEqual(lines[-1], "END")
        self.assertNotIn("TER", lines[: lines.index(atom_lines[0])])

    def test_multimodel_input_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "source.pdb"
            output = Path(tmp) / "output.pdb"
            source.write_text(
                "MODEL        1\n" + _atom(1, "B", 1, 0.0) + "\n" + _atom(2, "A", 1, 1.0) + "\nENDMDL\nEND\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(RuntimeError, "Multi-model"):
                stage22._copy_or_reorder_pdb(source, output, peptide_chain="B", target_chain="A")


if __name__ == "__main__":
    unittest.main()
