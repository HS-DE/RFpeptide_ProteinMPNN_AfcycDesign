# FGA 环肽路线脚本人工检查包

快照日期：2026-07-16

这个目录用于人工复核，不是新的运行目录。请先阅读本文件，再查看
`脚本逐项清单.csv` 和 `SHA256SUMS.csv`。

## 状态标签

| 标签 | 含义 |
| --- | --- |
| `CURRENT_REVIEW_SNAPSHOT` | 当前修复后源码快照，仍需人工检查 |
| `CURRENT_CORRECTED_SMOKE_REFERENCE` | runtime-locked N1-v2 runner 样例 |
| `CURRENT_END_TO_END_PROVENANCE_EVIDENCE` | N1-v3 Stage 0/ContigMap/tensor/TRB/PDB 端到端闭环证据 |
| `ACTUAL_RUNTIME_SNAPSHOT` | 实际被推理进程加载的 home RFdiffusion 源码 |
| `NON_RUNTIME_REFERENCE` | 可见但未被旧 N20 inference 导入的本地参考副本 |
| `ARCHIVED_INVALID_DO_NOT_RUN` | 从错误路线归档复制，只能审计，禁止执行 |
| `REFERENCE_ONLY` | 配置、任务书或归档证据，仅供核对 |

## 目录

```text
01_项目脚本当前快照/       全部 46 个项目脚本
02_RFdiffusion关键运行文件/ 实际 home runtime 6 个文件，并保留本地非运行副本用于差异检查
03_旧错误路线生成脚本/     归档中实际存在的 100 个 job/runner，禁止执行
04_旧N20假preflight样例_禁止运行/ 旧 runner；preflight 与 inference 使用不同源码
04_修复后运行脚本样例/     runtime-locked N1-v2 Stage 1 runner
05_配置快照/              project.yaml
06_任务书快照/            当前任务书副本
07_归档审计参考/          根因说明、禁止使用标记和移动清单
08_运行时分叉与修复证据/   false-preflight 根因、新 runtime audit 和严格 Stage 2 证据
09_Stage1到Stage2端到端闭环证据_N1_v3/ 完整 N1-v3 产物、mapping export、QC 和日志
脚本逐项清单.csv           每一个复制文件的来源、用途、状态和 SHA-256
SHA256SUMS.csv             副本完整性校验表
```

## RFpeptides 主路线脚本

| 脚本 | Stage | 用途 |
| --- | --- | --- |
| `18_discover_rfpeptides_target_sites.py` | Stage -1 | native context 中重新发现位点，计算 RSA、几何、遮挡、chemical anchor 和 optional fpocket 旁证 |
| `19_prepare_rfpeptides_article_inputs.py` | Stage 0 | 生成 target crop、hotspot 定义和 crop renumbering mapping |
| `20_make_rfpeptides_article_jobs.py` | Stage 1 准备 | 生成 RFpeptides shell runner，并在启动前验证 hotspot index 已加 binder length |
| `21_collect_rfpeptides_backbones.py` | Stage 2 | 解析 PDB/TRB、按 target-chain order 映射 Site_2/hotspot、检查 direct contact、宏环几何和 clash |
| `22_prepare_proteinmpnn_jobs.py` | Stage 3A | 从 Stage 2 pass backbone 准备 ProteinMPNN-only 或 ProteinMPNN-FastRelax jobs |
| `23_collect_proteinmpnn_sequences.py` | Stage 3C | 收集 ProteinMPNN 输出并重新检查序列、Site_2/hotspot、宏环和 clash |
| `24_stage3d1_sidechain_repack.py` | Stage 3D-1 | PyRosetta side-chain repack-only，不移动 backbone |
| `25_stage4_rosetta_interface_scoring.py` | Stage 4A-v2 | no-repack/no-minimization Rosetta score proxy、序列性质和 validation priority |
| `26_prepare_afcycdesign_jobs.py` | Stage 5A 准备 | 合并 Stage 4 候选并准备 sequence-based independent-recovery jobs |
| `27_collect_afcycdesign_validation.py` | Stage 5A 收集 | 解析 independent-recovery 模型并计算位点、姿态、拓扑和置信度 |
| `28_prepare_stage5_target_controls.py` | Stage 5 control | 准备 target-only single-sequence/MLM/MSA 控制 |
| `29_collect_stage5_target_controls.py` | Stage 5 control | 收集 target-only recovery control |
| `30_prepare_stage5b_target_conditioned_jobs.py` | Stage 5B 准备 | 准备 target-only-template conditioned recovery；不提供 peptide template/initial guess |
| `31_collect_stage5b_validation.py` | Stage 5B 收集 | 解析 Stage 5B、统一坐标系并核验 Stage 4 reference 是否真实命中修正后的 Site_2 |

## Stage 5 外部 runner

| 脚本 | 用途 |
| --- | --- |
| `external/run_afcycdesign_independent_recovery.py` | 执行单个 Stage 5A candidate/seed job |
| `external/run_afcycdesign_target_recovery_control.py` | 执行单个 target-only control job |
| `external/run_afcycdesign_target_conditioned_recovery.py` | 执行单个 Stage 5B target-conditioned candidate/seed job |

## 公共依赖

| 脚本 | 用途 |
| --- | --- |
| `common.py` | 项目路径、CSV、Markdown、YAML 配置和日志函数 |
| `pdb_utils.py` | PDB residue/atom 解析和坐标距离函数 |
| `region_utils.py` | FGA 区域定义和序列切片 |
| `sequence_filters.py` | 旧路线序列硬过滤函数 |
| `ranking.py` | 旧路线候选排序函数 |
| `__init__.py` | Python package marker |

## 旧 ColabDesign/Boltz/ColabFold 路线脚本

这些也按“全部项目脚本”要求复制，但不属于当前 RFpeptides 主路线。

| 脚本 | 用途 |
| --- | --- |
| `00_check_environment.py` | 检查项目环境和依赖 |
| `01_extract_fga_sequence.py` | 提取 FGA/P02671 序列 |
| `02_prepare_fga_regions.py` | 生成不同 FGA 区域 FASTA |
| `03_prepare_structures.py` | 下载并清理 native fibrinogen PDB |
| `04_map_fga_structure.py` | 建立 PDB-UniProt residue mapping |
| `05_select_surface_patches.py` | 旧路线选择表面 patches |
| `06_make_design_jobs.py` | 旧 ColabDesign 路线准备设计 jobs |
| `07_collect_raw_designs.py` | 收集旧路线 raw candidates |
| `08_filter_sequences.py` | 旧 Cys-Cys 环肽序列过滤 |
| `09_prepare_complex_prediction_jobs.py` | 准备旧路线复合物预测 |
| `10_score_complex_predictions.py` | 解析和评分旧路线复合物预测 |
| `11_negative_screen.py` | 汇总旧路线负筛选 |
| `12_rank_candidates.py` | 旧路线综合排名 |
| `13_export_final_report.py` | 旧路线报告导出 |
| `14_prepare_boltz_prediction_jobs.py` | 准备 Boltz jobs |
| `15_parse_boltz_predictions.py` | 解析 Boltz 输出 |
| `16_prepare_colabfold_prediction_jobs.py` | 准备 ColabFold/AF-Multimer jobs |
| `17_parse_colabfold_predictions.py` | 解析 ColabFold 输出 |
| `external/run_colabdesign_cyclic_binder.py` | 旧 ColabDesign cyclic binder runner |
| `external/run_colabdesign_chunk_batch.sh` | 旧 ColabDesign 分块 runner |
| `external/run_colabdesign_safe_batch.sh` | 旧 ColabDesign 安全批处理封装 |
| `external/run_boltz_batch.sh` | Boltz 批处理 runner |
| `external/run_colabfold_batch.sh` | ColabFold 批处理 runner |

## RFdiffusion/RFpeptides 关键 runtime

以下应优先检查 `实际运行时_home_RFdiffusion/`。`本地参考副本_rfd_macro_非N20实际推理包/`
只能用于解释旧运行时分叉，不能代表实际 inference。

| 文件 | 用途和检查重点 |
| --- | --- |
| `scripts/run_inference.py` | Hydra inference 入口；记录 runtime audit，并把 `idx_pdb`、chain 和 hotspot provenance 写入输出 |
| `rfdiffusion/contigs.py` | 构建 receptor-local 与 complex-global mapping |
| `rfdiffusion/inference/utils.py` | `get_idx0_hotspots()` 必须执行 `receptor_local_index + binderlen` |
| `rfdiffusion/inference/model_runners.py` | hotspot feature 必须复用修正后的 complex-global `self.hotspot_0idx` |
| `rfdiffusion/util.py` | 实际 PDB writer；确认使用传入的 `chain_idx` 和 `idx_pdb`，此前审计包遗漏了此文件 |
| `config/inference/base.yaml` | 基础 inference 配置上下文，不是独立运行脚本 |

## 旧生成脚本

`03_旧错误路线生成脚本` 中的 100 个文件来自已隔离的 N1000、batch01-05、
Stage 5A 和 Stage 5B 输出目录。它们保留原始绝对路径和旧参数，只用于回答：

1. 当时实际运行了什么；
2. 哪些脚本继承了错误 hotspot mapping；
3. 下游使用了哪些候选、seed、template 和 protocol 参数。

每一个旧脚本都在 `脚本逐项清单.csv` 中单独列出并带有用途，状态统一为
`ARCHIVED_INVALID_DO_NOT_RUN`。

## 建议人工检查顺序

1. `19_prepare_rfpeptides_article_inputs.py`：确认 Stage 0 hotspot 是 crop sequence position。
2. `rfdiffusion/contigs.py`：确认 receptor-local mapping 的定义。
3. `实际运行时_home_RFdiffusion/rfdiffusion/inference/utils.py`：确认 hotspot index 加 `binderlen`。
4. `实际运行时_home_RFdiffusion/rfdiffusion/inference/model_runners.py`：确认模型 feature 使用同一修正 index。
5. `实际运行时_home_RFdiffusion/rfdiffusion/util.py` 与 `scripts/run_inference.py`：确认实际 chain/numbering 写出路径。
6. `20_make_rfpeptides_article_jobs.py`：确认 preflight 与 inference 锁定同一 runtime root、commit 和哈希。
7. `04_修复后运行脚本样例`：核对实际 contig、target PDB、hotspot 参数和 runtime identity。
8. `21_collect_rfpeptides_backbones.py`：确认 runtime audit、target sequence identity、cyclic mask 和 direct contact 全部过门槛。
9. `23_collect_proteinmpnn_sequences.py`：确认 Stage 3 继续使用每个结构自己的 mapped residue number。
10. `25_stage4_rosetta_interface_scoring.py`：确认不移动 backbone，分数不被描述为实验结合能。
11. `26-31`：确认 Stage 5 identity/cache、template coverage、坐标对齐和 reference-site premise。

## 当前 smoke 状态

旧 N20 已证实为 false-preflight 并移入错误归档，0 条可用。当前只有一个
runtime-locked N1-v2 工程 smoke 完成闭环：

```text
1 parsed
1 provenance pass
1 strict contact/macrocycle/clash pass: RFpep_Site_2_0000
hotspot_0idx observed/expected: 98,100,101,102
target chain A = 86 aa exact Stage 0 sequence
peptide chain B = 17 aa and cyclic mask covers peptide only
```

对应的修复后复合物结构已收录为
`08_运行时分叉与修复证据/RFpep_Site_2_L17_17_0.pdb`，用于与 runtime audit
和严格 Stage 2 QC 表逐项交叉核验。

这条只证明修复后的工程链路可以闭环，不能据此直接扩大生产。

## N1-v3 端到端 provenance closure

`09_Stage1到Stage2端到端闭环证据_N1_v3/` 保存了新的独立 smoke。该版本在
N1-v2 基础上补齐：

- Stage 0 target PDB、mapping CSV、规范化 hotspot list 的 SHA256 锁定；
- 真实 Stage 0 PDB/contig/hotspots 构造的 ContigMap preflight；
- helper、独立 ContigMap 推导和实际 model hotspot tensor 的运行时硬断言；
- model tensor 建立后才最终写出的 JSON/TRB identical runtime audit；
- Stage 2 对 JSON、TRB mapping、writer chain/numbering 和真实 PDB 的逐热点闭环。

N1-v3 的 provenance 全部通过，但随机 backbone 没有接触 target/Site_2/hotspot，
所以 `pass_backbone_qc=false`，不得进入 Stage 3。这一区分是预期行为：
provenance pass 只证明生成条件与输出身份可追溯，不代表结构设计质量通过。
