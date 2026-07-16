# Stage 1 到 Stage 2 端到端 hotspot provenance closure

日期：2026-07-16

## 结论

`N1-v3` 完成了全新的单条 RFpeptides smoke run，没有复用 N1-v2 输出。

本次结果分成两个互相独立的结论：

1. hotspot provenance closure：通过。
2. 该随机 backbone 的结构质量 QC：未通过。

因此，这条结构证明新的生成与审计链能够闭环，但它本身不是可进入 Stage 3 的
backbone。Stage 2 明确记录的结构失败原因为：

```text
detached_from_target_crop
site_missed
hotspot_missed
```

不能因为 provenance 通过而把这条结构当成设计候选，也不能为了得到一个结构
pass 而删除或重抽本次 smoke 结果。

## 锁定的 Stage 0 输入

```text
target PDB SHA256:
4a2253962e5b380215180c2f463af9d13ffaf2bd9470f092a19e36b6fce87ccd

mapping CSV SHA256:
b407e30325721318f18d1bd9dbbbc4ebe65f225916ab7c7eb436c71aa9bc5cd3

normalized hotspots:
A82,A84,A85,A86

normalized hotspot-list SHA256:
85ad077878621b8709a23cb8d0c86fe201d8889ae3645ccb9909c5912f05f3dc
```

## 端到端核验结果

真实 Stage 0 PDB 和真实 contig `17-17 A1-86/0` 构造出的 ContigMap 为：

```text
binderlen: 17
receptor-local target indices: 0..85
complex-global target indices: 17..102
expected hotspot complex-global indices: 98,100,101,102
```

运行时三组 hotspot 索引完全一致：

```text
get_idx0_hotspots:              98,100,101,102
ContigMap independently derived: 98,100,101,102
actual model hotspot tensor:     98,100,101,102
```

每个 hotspot 的完整交叉映射为：

```text
Stage0 A82 -> TRB global 98  -> tensor 98  -> writer A82 -> PDB A82
Stage0 A84 -> TRB global 100 -> tensor 100 -> writer A84 -> PDB A84
Stage0 A85 -> TRB global 101 -> tensor 101 -> writer A85 -> PDB A85
Stage0 A86 -> TRB global 102 -> tensor 102 -> writer A86 -> PDB A86
```

Stage 2 还确认：

```text
runtime JSON == TRB runtime_audit: exact match
Stage 0 locked hashes: pass
TRB mapping == Stage 0 mapping: pass
runtime chain_idx == PDB chain order: pass
runtime idx_pdb == PDB residue-number order: pass
cyclic indices cover peptide only: pass
hotspot provenance closure: pass
```

## 文件说明

| 文件 | 用途 |
| --- | --- |
| `RFpep_Site_2_target.pdb` | 被 SHA256 锁定的 Stage 0 target 输入 |
| `RFpep_Site_2_crop_renumbering_mapping.csv` | 被锁定的 Stage 0 residue mapping |
| `RFpep_Site_2_hotspots.txt` | Stage 0 hotspot 说明 |
| `FGA_rfpeptides_stage1_jobs.csv` | 含 Stage 0 与 runtime 哈希的 Stage 1 任务表 |
| `run_rfpeptides_stage1_site2_runtimefix_smoke_N1_v3.sh` | N1-v3 实际运行脚本 |
| `RFpep_Site_2_L17_17_0.pdb` | 本次真实输出 PDB |
| `RFpep_Site_2_L17_17_0.trb` | 本次真实原始 TRB |
| `RFpep_Site_2_L17_17_0.runtime_audit.json` | 模型 tensor 建立后最终重写的 runtime audit |
| `RFpep_Site_2_L17_17_0.trb_mapping.json` | Stage 2 导出的可读 TRB mapping 与 hotspot crosswalk |
| `FGA_rfpeptides_backbones_qc.csv` | Stage 2 全量 QC；provenance pass、结构 QC fail |
| `FGA_rfpeptides_backbones_qc_pass.csv` | 空表；该 smoke 结构不得进入 Stage 3 |
| `run_rfpeptides_stage1_site2_runtimefix_smoke_N1_v3.complete.log` | 包含真实 ContigMap preflight 和完整推理输出的日志 |
| `stage2_hotspot_provenance_qc.complete.log` | Stage 2 完整运行日志 |

## 当前边界

本目录证明的是工程 provenance closure 已能够自动执行和硬失败，不证明每次
RFpeptides 采样都会命中 Site_2。恢复小规模生产前仍应先运行少量新批次，并同时
报告 provenance pass 数和结构 QC pass 数。
