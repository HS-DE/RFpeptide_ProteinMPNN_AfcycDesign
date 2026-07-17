# 路线配置说明

- `project.yaml`：Stage 00-19 的原有兼容配置，不能用于 Stage 20-31 生产运行。
- `project_legacy_disulfide.yaml`：旧 Cys-Cys 二硫键路线的完整归档副本，`production_use=false`。
- `rfpeptides_head_to_tail.yaml`：Stage 20-31 唯一允许的活跃生产配置。

Stage 20-31 的每个命令都必须显式传入：

```text
--project-config config/rfpeptides_head_to_tail.yaml
```

活跃配置由独立的严格 loader 解析，不会与旧 `DEFAULT_CONFIG` 或
`project.yaml` 深度合并。每个 Stage 20 run 还会把配置原始 SHA-256 和解析后
canonical SHA-256 写入 `route_manifest.json`，下游不一致时立即停止。

活跃路线中的 `omit_aas: "CX"` 表示当前项目主动排除 Cys 和未知残基 X。这是本路线的序列空间选择，不是 head-to-tail 酰胺环化在物理上的必需条件。
