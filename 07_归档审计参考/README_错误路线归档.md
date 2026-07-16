# RFpep_Site_2 错误路线隔离归档

归档日期：2026-07-16

状态：**INVALID FOR INTENDED SITE_2 DESIGN / 禁止重新接入当前生产路线**

## 为什么归档

修复前的 RFpeptides Stage 1 hotspot conditioning 将 86-aa target crop 中的
hotspot 序列位置直接作为复合物全局索引使用，没有加入前置 peptide 的长度。

对于长度为 `L` 的 peptide 和 crop position `p`：

```text
正确 complex 0-based index = L + (p - 1)
正确输出 PDB residue number = L + p
```

例如 `L=17` 时，Stage 0 hotspot `82,84,85,86` 应映射为输出 PDB residue
`99,101,102,103`，旧流程却使用了 `82,84,85,86`。因此 Stage 1 从生成时就
条件化到了错误 target 区域。

旧 Stage 2 及下游 QC 又直接用 crop position 匹配输出 PDB residue number，
重复了同一种偏移错误，造成错误区域上的结构被标记为 Site_2/hotspot pass。

## 影响边界

- Stage -1 site rediscovery：保留，不是本次错误起点。
- Stage 0 target crop、序列和 hotspot 定义：保留，可作为修复后 Stage 1 输入。
- 修复前 Stage 1 RFpeptides backbones：对 intended Site_2 conditioning 无效。
- 旧 Stage 2 contact/hotspot pass：假阳性风险，不能继续使用。
- Stage 3 ProteinMPNN/side-chain repack：继承 off-target backbone，不能作为
  Site_2 序列设计证据。
- Stage 4 Rosetta score：只描述旧错误 pose 的界面，不能作为 Site_2 排名证据。
- Stage 5A/5B：候选输入前提无效；预测仅保留为 protocol engineering record。

Stage 5B 的 25/25 jobs 和 125/125 models 技术上完成，但 5/5 Stage 4 reference
designs 在修正映射后均为 Site_2 contacts=0、hotspot contacts=0，因此不能把
Stage 5B 结果解释为对 5 个正确 Site_2 候选的验证或否定。

## 本次归档规模

归档前盘点：

```text
affected_entries: 14 directories/subdirectories
files: 99,199
size: 3,813,640,162 bytes (about 3.55 GiB)
```

`archive_manifest.csv` 记录每个条目的原始路径、归档路径、文件数和字节数。
文件采用同一磁盘内移动，没有删除，也没有创建一份大体积重复副本。

## 被隔离的内容

```text
rfpeptides_article_route_clean_20260612/01_rfpeptides_jobs
rfpeptides_article_route_clean_20260615_fpocket/01_rfpeptides_jobs
rfpeptides_article_route_clean_20260615_fpocket/02_rfpeptides_backbones
rfpeptides_article_route_clean_20260615_fpocket/logs
rfpeptides_article_route_clean_20260615_fpocket_stage1_N1000_no_traj
rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch01
rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch02
rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch03
rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch04
rfpeptides_article_route_clean_20260623_stage1_N10000_L12_24_batch05
rfpeptides_article_route_clean_20260623_stage5_batch01_batch02
rfpeptides_article_route_clean_20260623_stage5A_v2_batch01_batch02
rfpeptides_article_route_clean_20260623_stage5B_batch01_batch02
rfpeptides_article_route_clean_20260715_hotspotfix_smoke_false_preflight
```

## 明确保留在活动结果区的内容

```text
results/rfpeptides_article_route_clean_20260612/00_site_discovery
results/rfpeptides_article_route_clean_20260612/00_target_inputs
results/rfpeptides_article_route_clean_20260615_fpocket/00_site_discovery
results/rfpeptides_article_route_clean_20260615_fpocket/00_target_inputs
results/rfpeptides_article_route_clean_20260623_stage5_target_controls_v1
results/rfpeptides_article_route_clean_20260716_runtimefix_smoke_N1_v2
```

其中 target-only controls 不依赖旧 peptide reference pose，可作为协议诊断记录。

2026-07-16 复核发现，旧 `20260715_hotspotfix_smoke` 的静态 preflight 从
`C:/SH/peptide_str/rfd_macro` 导入了已修改 helper，但实际 inference 从
`/home/luomi/fga_model_envs/rfpeptides/RFdiffusion` 导入了另一份未修复 package。
因此它是 false-preflight 结果，20 条 backbone 和 3 条 strict-pass 记录全部
转入本归档，不得作为修复后 smoke。新的 `20260716_runtimefix_smoke_N1_v2`
强制锁定单一 runtime root，并在 TRB/JSON 中记录源码路径、哈希、hotspot index、
chain identity 和 cyclic mask。

## 重新开始边界

路线从修复后的 Stage 1 重新开始。扩大生成前必须依次通过：

1. 静态 hotspot offset preflight；
2. TRB input-to-complex mapping 核对；
3. 输出 PDB 中实际 Site_2/hotspot residue mapping 核对；
4. 少量结构的人工 PyMOL 检查；
5. 修正版 Stage 2 必须证明 peptide 接触正确 Site_2/hotspot。

上述检查未全部通过前，不得进入 ProteinMPNN、Rosetta Stage 4 或 Stage 5。
