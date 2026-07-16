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

这只是工程 N1 smoke，不是生产候选，也不足以直接恢复 ProteinMPNN 或 Stage 5。

