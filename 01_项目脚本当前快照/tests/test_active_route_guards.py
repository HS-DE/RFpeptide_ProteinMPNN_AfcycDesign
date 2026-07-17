from __future__ import annotations

import ast
from pathlib import Path
import sys
import tempfile
import unittest


SCRIPTS_DIR = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

from common import assert_active_route_path  # noqa: E402


class ActiveRoutePathTests(unittest.TestCase):
    def test_archive_markers_are_blocked_case_insensitively(self) -> None:
        markers = [
            "_ARCHIVED_INVALID_",
            "03_旧错误路线生成脚本",
            "04_旧N20假preflight样例_禁止运行",
            "旧错误路线",
            "禁止运行",
        ]
        for marker in markers:
            with self.subTest(marker=marker), self.assertRaisesRegex(RuntimeError, "blocked"):
                assert_active_route_path(Path("C:/safe") / marker / "input.csv", "fixture", must_exist=False)

    def test_allowed_existing_path_passes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "active" / "input.csv"
            path.parent.mkdir()
            path.write_text("ok\n", encoding="utf-8")
            self.assertEqual(assert_active_route_path(path, "fixture"), path.resolve())

    def test_symlink_into_forbidden_tree_is_blocked(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            forbidden = root / "_archived_invalid_" / "input.csv"
            forbidden.parent.mkdir()
            forbidden.write_text("bad\n", encoding="utf-8")
            link = root / "active_link.csv"
            try:
                link.symlink_to(forbidden)
            except OSError as exc:
                self.skipTest(f"Symlink creation is not available: {exc}")
            with self.assertRaisesRegex(RuntimeError, "blocked"):
                assert_active_route_path(link, "symlink fixture")


class RequiredCliInputTests(unittest.TestCase):
    REQUIRED_OPTIONS = {
        "20_make_rfpeptides_article_jobs.py": ["--input-root", "--output-root", "--selected-sites", "--rfpeptides-root", "--batch-id", "--project-config"],
        "21_collect_rfpeptides_backbones.py": ["--stage0-root", "--stage1-root", "--selected-sites", "--project-config"],
        "22_prepare_proteinmpnn_jobs.py": ["--stage2-root", "--selected-backbones", "--dl-binder-design-root", "--project-config"],
        "23_collect_proteinmpnn_sequences.py": ["--stage0-root", "--stage3-root", "--selected-backbones", "--stage3-jobs-csv", "--project-config"],
        "24_stage3d1_sidechain_repack.py": ["--stage0-root", "--stage3-root", "--selected-backbones", "--project-config"],
        "25_stage4_rosetta_interface_scoring.py": ["--stage0-root", "--stage3-root", "--selected-backbones", "--project-config"],
        "26_prepare_afcycdesign_jobs.py": ["--source-run-root", "--stage0-root", "--output-root", "--project-config"],
        "27_collect_afcycdesign_validation.py": ["--stage5-root", "--stage0-root", "--project-config"],
        "28_prepare_stage5_target_controls.py": ["--source-run-root", "--stage0-root", "--output-root", "--project-config"],
        "29_collect_stage5_target_controls.py": ["--control-root", "--stage0-root", "--project-config"],
        "30_prepare_stage5b_target_conditioned_jobs.py": ["--stage5a-root", "--stage0-root", "--output-root", "--project-config"],
        "31_collect_stage5b_validation.py": ["--stage5b-root", "--stage0-root", "--project-config"],
    }

    def test_all_production_entrypoints_require_upstream_roots(self) -> None:
        for filename, required_options in self.REQUIRED_OPTIONS.items():
            with self.subTest(script=filename):
                tree = ast.parse((SCRIPTS_DIR / filename).read_text(encoding="utf-8"), filename=filename)
                required_found: set[str] = set()
                for node in ast.walk(tree):
                    if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
                        continue
                    if node.func.attr != "add_argument" or not node.args:
                        continue
                    option_node = node.args[0]
                    if not isinstance(option_node, ast.Constant) or not isinstance(option_node.value, str):
                        continue
                    required_keyword = next((item for item in node.keywords if item.arg == "required"), None)
                    if (
                        required_keyword is not None
                        and isinstance(required_keyword.value, ast.Constant)
                        and required_keyword.value.value is True
                    ):
                        required_found.add(option_node.value)
                for option in required_options:
                    self.assertIn(option, required_found)

    def test_all_production_entrypoints_use_shared_route_guards(self) -> None:
        for filename in self.REQUIRED_OPTIONS:
            with self.subTest(script=filename):
                tree = ast.parse((SCRIPTS_DIR / filename).read_text(encoding="utf-8"), filename=filename)
                called_names = {
                    node.func.id
                    for node in ast.walk(tree)
                    if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
                }
                self.assertIn("assert_active_route_path", called_names)
                self.assertTrue(
                    {"load_active_route_config", "validate_route_project_config"} & called_names,
                    f"{filename} does not use the strict active-route config loader",
                )
                self.assertTrue(
                    {"load_route_manifest", "write_route_manifest"} & called_names,
                    f"{filename} does not read or create a route manifest",
                )

    def test_no_historical_run_identity_remains_in_active_scripts(self) -> None:
        forbidden = ["n1000", "rfpep_site_2_0007", "batch01", "batch02", "202606", "5529"]
        for filename in self.REQUIRED_OPTIONS:
            text = (SCRIPTS_DIR / filename).read_text(encoding="utf-8").casefold()
            for token in forbidden:
                with self.subTest(script=filename, token=token):
                    self.assertNotIn(token, text)


if __name__ == "__main__":
    unittest.main()
