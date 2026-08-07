# FGA RFpeptides Stage 5 流程记录

更新时间：2026-08-06

## 1. Stage 5 的目的

Stage 5 不再生成 peptide 序列，而是检查 Stage 4 留下的序列能否在
AfCycDesign/AlphaFold 模型中形成合理的头尾环化 peptide，并回到 FGA 的
RFpep_Site_2 与 hotspot 附近。

本阶段属于计算结构恢复测试，不是实验结合验证。任何 Stage 5 支持都不能单独
把 peptide 定义为 final candidate；后续仍需要正交结构预测、负筛选与实验验证。

当前有效上游输入来自修复后的 head-to-tail 路线：

```text
Stage 4 hard-QC pass: 2372
Stage 4 official top validation candidates: 5
top-5 backbone families: 5
target crop length: 86 aa
hotspots: A82, A84, A85, A86
```

## 2. Stage 5A：sequence-only independent recovery

Stage 5A 的定义是只给 target sequence 与 peptide sequence，不给 target template、
peptide template 或 initial guess。peptide 仅通过 cyclic positional offset 表示为
头尾环化链。

初期 single-sequence + MLM 0.15 的 125-model 运行属于旧路线历史诊断，不能作为
修复后 head-to-tail 路线的候选证据。随后完成了缓存身份、sequence hash、协议
hash、pLDDT 尺度、PAE 字段和输入序列交叉核验等工程修复。

修复后 active-route top5 的 Stage 5A 任务已经准备在：

```text
results/rfpeptides_head_to_tail_v1_20260731_stage5A_top5/
```

该目录目前只有 manifest、protocol audit、job specs 和运行脚本；没有执行 active-route
Stage 5A 的 125 个模型。由于 86 aa target crop 的 sequence-only 恢复能力很弱，
该路线暂不作为当前优先判断依据。

## 3. target-only sequence controls

为区分“peptide 不对”和“target 自身无法由当前协议恢复”，准备了四组 target-only
sequence recovery controls：

```text
A: single sequence + MLM 0.15
B: single sequence + no MLM
C: real unpaired homolog MSA + MLM 0.15
D: real unpaired homolog MSA + no MLM
```

A/B 共完成 30 个模型，均未达到 target/Site_2/hotspot 的内部恢复阈值：

| group | models | median target RMSD (A) | median Site_2 RMSD (A) | median hotspot RMSD (A) |
| --- | ---: | ---: | ---: | ---: |
| A, single sequence + MLM | 15 | 36.525 | 65.135 | 81.704 |
| B, single sequence + no MLM | 15 | 10.984 | 18.817 | 26.015 |

C/D 因没有真实 FGA homolog MSA 而保持阻断，没有用伪造 paired MSA 替代。

结论：single-sequence target-only 基线不足以稳定恢复当前 86 aa crop，因此
sequence-only complex failure 不能直接解释为 peptide failure。

## 4. 旧 Stage 5B：masked target-template recovery

旧 Stage 5B 给模型加载了 86 aa target PDB，但实际模板选项为：

```text
template_sequence_masked: true
template_sidechains_masked: true
template_interchain_features_masked: true
```

因此它不是“给出完整 target 结构”的测试。模型可见 target backbone template 的一部分
几何信息，但 target residue identity、侧链和模板链间特征被屏蔽。peptide 仍没有
template 或 initial guess。

两轮运行结果：

```text
top5: 25 seed jobs, 125 models, 0/5 candidate support
all Stage 4 pass: 11860 seed jobs, 59300 models, 0/2372 candidate support
```

全量运行中：

```text
target recovery count: 0/59300
Site_2 contact count: 15865/59300
hotspot/same-site contact count: 2769/59300
strong or moderate peptide pose recovery: 0/59300
macrocycle geometry pass: 59300/59300
```

这些 contact 数不能单独证明结合位点恢复，因为对应模型的 target 并没有先通过结构
恢复门槛。该 59,300-model 结果应保留为 masked-template 协议诊断，不应据此宣布
2,372 条 peptide 全部失败。

## 5. C0-C3 target-context controls

随后用 target-only、无 peptide 的小规模控制检查 target template 的实际作用。每组
使用 1 个 seed 和全部 5 个 multimer-v3 model parameter sets，共 20 个模型。

| control | target context | template information | median FGA RMSD (A) | median Site_2 RMSD (A) | median hotspot RMSD (A) |
| --- | --- | --- | ---: | ---: | ---: |
| C0 | G115-G200, 86 aa | masked | 7.950 | 11.303 | 16.594 |
| C1 | G115-G200, 86 aa | full | 1.383 | 2.393 | 3.500 |
| C2 | G27-G200, 174 aa | full | 53.051 | 27.317 | 28.979 |
| C3 | G27-G200 + H118-H184 + I68-I127, 301 aa | full | 3.920 | 5.390 | 7.573 |

没有一组通过原先的全部严格阈值，但不能把四组简单归为同一种失败：

- C0 显示 masked template 下 86 aa crop 明显漂移。
- C1 保留完整 template 后，全局 target RMSD 明显改善到约 1.3-2.5 A；Site_2 与
  hotspot 仍略高于原先 2 A 阈值，局部置信度也偏低。
- C2 的孤立 174 aa FGA 片段出现整体构象崩塌，说明扩展单链不等于提供稳定的原生
  结构环境。
- C3 中 3/5 个模型保留了近似原生复合物架构：FGA RMSD 3.227-3.920 A，FGA 对齐后
  partner RMSD 2.670-3.582 A；另外 2 个模型明显崩塌。

因此，下一步应比较完整 86 aa crop template 的 C1 与包含邻近 FGB/FGG 片段的 C3，
而不是继续把旧 masked-template Stage 5B 当作完整 target-conditioned 验证。

## 6. 当前 Stage 5B-v2：top5 × C1/C3

新增脚本：

```text
scripts/34_prepare_stage5b_v2_context_jobs.py
scripts/external/run_afcycdesign_stage5b_v2_context_recovery.py
scripts/35_collect_stage5b_v2_context_validation.py
```

本轮固定使用相同的 Stage 4 official top5，并分别在 C1 与 C3 中预测：

```text
candidates: 5
contexts per candidate: 2
seed jobs: 10
model parameter sets per job: 5
planned model predictions: 50
requested recycles: 6
forward passes: 7
dropout: false
MLM: false
target MSA mode: single_sequence
peptide MSA mode: single_sequence
```

硬性协议：

- C1 target template 覆盖 86/86，peptide 是链 B，cyclic chain index=1。
- C3 target template 覆盖 301/301，target 是 A/B/C，peptide 是链 D，cyclic chain index=3。
- target template 的 sequence、backbone、sidechain 与 target-target interchain features
  全部保留。
- peptide template coverage=0；不加载 peptide 原设计坐标。
- use_initial_guess=false。
- Stage 4 PDB 只允许在预测完成后用于 target 对齐、peptide RMSD 与 contact recovery，
  runner 不读取其坐标。
- candidate/job/cache identity 同时绑定 peptide sequence hash、context hash、Stage 4
  identity、seed、recycles、commit 与 template protocol。

准备结果与真实输入预检：

```text
candidate-context rows validated: 10/10
seed jobs validated: 10/10
C1 target template coverage: 86/86
C3 target template coverage: 301/301
peptide template coverage: 0/10 jobs
ColabDesign commit: 5ab4efaba2321a6c3c314b82d2fff8e0241f5c2d
prediction run status: not started
```

输出目录：

```text
results/rfpeptides_head_to_tail_v1_20260806_stage5B_v2_top5_C1_C3/
  07_structure_validation_target_context_conditioned/
```

## 7. 正式运行命令

先运行 10 个任务的静态与模板张量预检：

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design

bash results/rfpeptides_head_to_tail_v1_20260806_stage5B_v2_top5_C1_C3/07_structure_validation_target_context_conditioned/jobs/check_stage5B_v2_C1_C3_protocol.sh
```

确认后，在一个 GPU 上依次运行 C1 与 C3：

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design

CUDA_VISIBLE_DEVICES=0 RUN_STAGE5B_V2_PREDICTIONS=YES bash \
  results/rfpeptides_head_to_tail_v1_20260806_stage5B_v2_top5_C1_C3/07_structure_validation_target_context_conditioned/jobs/run_stage5B_v2_C1_C3_all.sh
```

也可以在两个 GPU 上分别运行 C1 和 C3；不要再同时运行总脚本：

```bash
# GPU 0: C1
CUDA_VISIBLE_DEVICES=0 RUN_STAGE5B_V2_PREDICTIONS=YES bash \
  results/rfpeptides_head_to_tail_v1_20260806_stage5B_v2_top5_C1_C3/07_structure_validation_target_context_conditioned/jobs/run_stage5B_v2_C1_crop86_full_template_all.sh

# GPU 1: C3
CUDA_VISIBLE_DEVICES=1 RUN_STAGE5B_V2_PREDICTIONS=YES bash \
  results/rfpeptides_head_to_tail_v1_20260806_stage5B_v2_top5_C1_C3/07_structure_validation_target_context_conditioned/jobs/run_stage5B_v2_C3_native_GHI301_full_template_all.sh
```

全部预测结束后收集：

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate fga_stage1_fpocket

python scripts/35_collect_stage5b_v2_context_validation.py \
  --stage5b-v2-root results/rfpeptides_head_to_tail_v1_20260806_stage5B_v2_top5_C1_C3/07_structure_validation_target_context_conditioned \
  --project-config config/rfpeptides_head_to_tail.yaml
```

## 8. 本轮结果如何解释

1. C1 成功、C3 成功：说明 peptide pose 在受限 crop 与原生邻近链背景下都获得支持。
2. C1 成功、C3 失败：只能称为 crop-only support，可能受空间限制影响。
3. C1/C3 target 均恢复，但 peptide 不回到 Site_2：才可把证据更明确地指向 peptide
   pose 未恢复。
4. C3 target 自身未恢复：C3 对该 peptide 的否定结果仍不可判定，不能据此淘汰。
5. 即使获得 strong support，也只是进入下一轮正交预测/筛选，不是 final peptide。

## 9. 全量 C1/C3 campaign 准备，2026-08-06

在 top5 C1/C3 协议通过输入检查后，同一协议被扩展到全部 2,372 条
Stage 4 hard-QC pass。全量模式不增加 seed 数，仍然是每个
candidate-context 1 个 seed、5 个 AlphaFold multimer-v3 parameter sets：

```text
selection_mode: all_stage4_pass
candidates: 2372
contexts: C1_crop86_full_template, C3_native_GHI301_full_template
candidate-context pairs: 4744
seed jobs: 4744
model predictions: 23720
requested recycles: 6
forward passes per model: 7
job shards: 6
representative preflight specs: 26
```

输出目录：

```text
results/rfpeptides_head_to_tail_v1_20260806_stage5B_v2_all2372_C1_C3/
  07_structure_validation_target_context_conditioned/
```

准备阶段已经完成以下检查：

1. 2,372 条 Stage 4 candidate 的序列、PDB、链角色与 SHA-256 全部通过。
2. 4,744 个 candidate-context ID 和 job ID 全部唯一。
3. 六个 shard 均同时包含 C1 和 C3：shard 1/2 各含 396+396 个 jobs，
   shard 3-6 各含 395+395 个 jobs。
4. 26 个代表性 preflight 覆盖 C1/C3 和实际出现的 12-24 aa 肽长；C1
   template coverage 为 86/86，C3 为 301/301，peptide coverage 始终为 0。
5. Stage 35 对完整 4,744 行 manifest/job/spec 输入执行只读校验并通过。

全量预测尚未启动。推荐先按每张 GPU 一个进程顺序跑三个 shard，稳定性优先：

```bash
# GPU 0 窗口
cd /mnt/c/SH/fga_cyclic_peptide_design
ROOT=results/rfpeptides_head_to_tail_v1_20260806_stage5B_v2_all2372_C1_C3
JOBDIR="$ROOT/07_structure_validation_target_context_conditioned/jobs"

for s in 01 03 05; do
  CUDA_VISIBLE_DEVICES=0 RUN_STAGE5B_V2_PREDICTIONS=YES bash \
    "$JOBDIR/run_stage5B_v2_C1_C3_shard_${s}_of_06.sh" || break
done
```

```bash
# GPU 1 窗口
cd /mnt/c/SH/fga_cyclic_peptide_design
ROOT=results/rfpeptides_head_to_tail_v1_20260806_stage5B_v2_all2372_C1_C3
JOBDIR="$ROOT/07_structure_validation_target_context_conditioned/jobs"

for s in 02 04 06; do
  CUDA_VISIBLE_DEVICES=1 RUN_STAGE5B_V2_PREDICTIONS=YES bash \
    "$JOBDIR/run_stage5B_v2_C1_C3_shard_${s}_of_06.sh" || break
done
```

这些脚本可以重复执行；runner 只有在 metadata、协议身份、模型记录及 PDB/NPZ
全部匹配时才跳过已完成 job。不要同时运行总脚本和 shard 脚本，以免两个进程
竞争同一输出目录。

查看进度：

```bash
ROOT=results/rfpeptides_head_to_tail_v1_20260806_stage5B_v2_all2372_C1_C3/07_structure_validation_target_context_conditioned
echo "completed jobs: $(find "$ROOT/predictions" -type f -name run_metadata.json | wc -l) / 4744"
echo "written model PDBs: $(find "$ROOT/predictions" -type f -name '*.pdb' | wc -l) / 23720"
```

全部完成后收集：

```bash
cd /mnt/c/SH/fga_cyclic_peptide_design
source ~/fga_model_envs/miniforge3/etc/profile.d/conda.sh
conda activate fga_stage1_fpocket

python scripts/35_collect_stage5b_v2_context_validation.py \
  --stage5b-v2-root results/rfpeptides_head_to_tail_v1_20260806_stage5B_v2_all2372_C1_C3/07_structure_validation_target_context_conditioned \
  --project-config config/rfpeptides_head_to_tail.yaml
```

基于 C0-C3 control 的单 job 时间，本轮在每张 GPU 单进程时只能粗略估计为约
两天量级；实际时间取决于 GPU、并行度和 C3 301-aa context 的开销。正式输出
预计需要数 GB 磁盘空间，运行前应预留至少 10 GB。该 campaign 是大范围恢复
信号搜索，不是 final peptide selection。

## 10. 旧 top5 任务兼容修复，2026-08-07

top5 C1/C3 运行在完成全部 5 个 C1 和前 3 个 C3 后中断。报错不是模型预测失败，
而是旧 top5 job spec 生成时尚无 `stage5_selection_mode` 字段，之后为全量 campaign
增强的共享 runner 将该字段纳入必需协议身份，因而拒绝了剩余旧 spec。

兼容策略严格限定为旧协议版本
`stage5B_v2_C1_C3_full_target_template_top5_v1`：缺失字段只读解释为
`top_validation`。全量协议
`stage5B_v2_C1_C3_full_target_template_allpass_v1` 仍必须显式记录
`stage5_selection_mode=all_stage4_pass`，不允许使用兼容默认值。

runner 与 collector 同时兼容旧 top5 spec、metadata 和 model metrics。已完成的
8/8 个 top5 job 均重新通过缓存身份检查，因此重新运行原 master script 时会跳过
这 8 个任务，只补跑最后 2 个 C3 job。top5 10/10 preflight 重新通过；全量任务的
4,744 份 spec 均显式包含正确字段，C1/C3 各 2,372 份，代表性 26/26 preflight
重新通过。全量任务不需要重新生成。
