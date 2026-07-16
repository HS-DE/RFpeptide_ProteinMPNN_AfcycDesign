from __future__ import annotations

from typing import Any, Dict, List


def build_region_rows(sequence: str, config: Dict[str, Any]) -> List[Dict[str, Any]]:
    regions = config["target_regions"]
    ordered = [
        ("full_length_1_866", regions["full_length"]),
        ("extracellular_20_866", regions["extracellular"]),
        ("chain_36_866", regions["main_chain"]),
    ]
    rows: List[Dict[str, Any]] = []
    for name, spec in ordered:
        start = int(spec["start"])
        end = int(spec["end"])
        region_seq = sequence[start - 1 : end]
        rows.append(
            {
                "region_name": name,
                "start": start,
                "end": end,
                "length": len(region_seq),
                "use_for_design": bool(spec.get("use_for_design", False)),
                "priority": spec.get("priority", "preferred" if spec.get("preferred") else "secondary"),
                "note": spec.get("note", ""),
                "sequence": region_seq,
            }
        )
    return rows


def fasta_header(region_name: str, config: Dict[str, Any]) -> str:
    return f"{config['project']['target_gene']}_{config['project']['target_uniprot']}_{region_name}"
