# Stage 1 运行时分叉根因

日期：2026-07-16

## 结论

旧 `20260715_hotspotfix_smoke` 不是修复后结果，不能继续使用。

旧 runner 在 `/mnt/c/SH/peptide_str/rfd_macro` 中运行。其静态 preflight 使用
`python -`，因此导入该当前目录下已经修改的 `rfdiffusion`。但真正执行
`./scripts/run_inference.py` 时，Python 从下面的 editable/install 路径导入包：

```text
/home/luomi/fga_model_envs/rfpeptides/RFdiffusion
```

home runtime 的 `get_idx0_hotspots()` 当时仍返回 receptor-local index，没有
加 `binderlen`。因此出现“preflight PASS，但实际 inference 仍使用错误 hotspot”
的 false positive。

## 链角色

ChatGPT 根据本地参考副本推断 `A=peptide, B=target`，这个推断不适用于实际
home runtime。实际 home `model_runners.py` 保留固定 target 的原 chain A，给
新设计 peptide 分配 chain B。因此实际输出是：

```text
chain A = 86-aa FGA target
chain B = 12-24-aa cyclic peptide
```

旧 Stage 2 把 A 识别为 target，在“链角色”这一点上恰好符合实际 PDB；但这
不能挽救错误的 hotspot conditioning。旧 N20 的 3 条 strict-pass 只是事后
碰巧接近真实 hotspot，不能证明生成时使用了正确热点。

## 修复

1. 修正实际 home runtime 的 `get_idx0_hotspots()`，使用
   `receptor_con_hal_idx0 + binderlen`。
2. `model_runners.py` 对 hotspot 数量和 complex-global index 做硬校验。
3. `run_inference.py` 为每个设计写入 `runtime_audit.json`，并在 TRB 中保存
   同一审计信息。
4. Stage 1 job maker 锁定单一 runtime root、Git commit 和四个源码 SHA-256。
5. Stage 2 要求 runtime audit 通过，并核验 target 序列、chain order、
   peptide-only cyclic mask 和实际 hotspot indices。

## 首个闭环 smoke

`20260716_runtimefix_smoke_N1_v2` 的检查结果：

```text
runtime_audit: pass
target: chain A, 86 aa, exact Stage 0 sequence
peptide: chain B, 17 aa
hotspot_0idx observed/expected: 98,100,101,102
cyclic mask: peptide only
direct Site_2 contacts: 5
direct hotspot contacts: 4
terminal C-N distance: 1.146 A
severe clash: no
strict Stage 2: pass
```

本目录中的 `RFpep_Site_2_L17_17_0.pdb` 是该 N1-v2 smoke 的实际输出结构，
用于核验 target/peptide 链身份、Site_2/hotspot 接触、cyclic geometry 和
runtime audit 是否彼此一致。

这只是工程 N1 smoke，不是生产候选，也不足以直接恢复 ProteinMPNN 或 Stage 5。

## 后续 N1-v3 闭环

N1-v2 之后发现静态 preflight 和 Stage 2 尚未自动读取真实 TRB mapping。该缺口
已在 `../09_Stage1到Stage2端到端闭环证据_N1_v3/` 中继续修复：模型输入 tensor
位置被直接读取并硬断言，最终 JSON/TRB audit 必须完全一致，Stage 2 逐热点验证
`Stage 0 -> TRB global index -> tensor -> writer -> PDB`。

N1-v3 的 provenance closure 通过，但随机结构的 contact QC 失败，因此不能进入
Stage 3。请勿用 N1-v2 的结构 pass 代替 N1-v3 的 provenance 证据，也不要用
N1-v3 的 provenance pass 宣称该 backbone 设计成功。
