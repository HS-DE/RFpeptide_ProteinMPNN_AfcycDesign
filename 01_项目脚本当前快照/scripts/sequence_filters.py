from __future__ import annotations

from typing import Any, Dict, Optional


HYDROPHOBIC = set("AVILMFWY")
LOW_COMPLEXITY_MOTIFS = ("AAAA", "KKKK", "RRRR", "DDDD", "EEEE", "GGGG", "SSSS", "GSGSGS", "PAPAPA")


def net_charge(seq: str) -> int:
    return sum(1 for aa in seq if aa in "KR") - sum(1 for aa in seq if aa in "DE")


def max_hydrophobic_run(seq: str) -> int:
    best = 0
    current = 0
    for aa in seq:
        if aa in HYDROPHOBIC:
            current += 1
            best = max(best, current)
        else:
            current = 0
    return best


def has_low_complexity(seq: str) -> bool:
    if any(motif in seq for motif in LOW_COMPLEXITY_MOTIFS):
        return True
    for size in (2, 3):
        for i in range(0, len(seq) - size * 3 + 1):
            motif = seq[i : i + size]
            if motif * 3 in seq:
                return True
    return False


def filter_candidate_sequence(seq: str, config: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    cfg = (config or {}).get("sequence_filters", {}) if config else {}
    peptide_cfg = (config or {}).get("peptide_design", {}) if config else {}
    seq = "".join(ch for ch in str(seq).strip().upper() if ch.isalpha())
    length_min = int(peptide_cfg.get("core_length_min", 10))
    length_max = int(peptide_cfg.get("core_length_max", 18))
    charge_min = int(cfg.get("net_charge_min", -3))
    charge_max = int(cfg.get("net_charge_max", 3))
    max_hydro = int(cfg.get("max_hydrophobic_run", 4))
    max_w = int(cfg.get("max_w_count", 1))
    max_m = int(cfg.get("max_m_count", 1))

    starts_with_cys = seq.startswith("C")
    ends_with_cys = seq.endswith("C")
    internal_cys_count = seq[1:-1].count("C") if len(seq) >= 2 else seq.count("C")
    core_length_pass = length_min <= len(seq) <= length_max
    charge = net_charge(seq)
    charge_pass = charge_min <= charge <= charge_max
    hydrophobic_run = max_hydrophobic_run(seq)
    hydrophobicity_pass = hydrophobic_run <= max_hydro
    w_count = seq.count("W")
    m_count = seq.count("M")
    oxidation_pass = w_count <= max_w and m_count <= max_m
    low_complexity_flag = has_low_complexity(seq)
    poly_ionic_flag = any(motif in seq for motif in ("KKK", "RRR", "DDD", "EEE"))

    notes = []
    if not starts_with_cys:
        notes.append("首位不是 Cys")
    if not ends_with_cys:
        notes.append("末位不是 Cys")
    if internal_cys_count:
        notes.append("内部含额外 Cys")
    if not core_length_pass:
        notes.append(f"长度 {len(seq)} 不在 {length_min}-{length_max}")
    if not charge_pass:
        notes.append(f"净电荷 {charge} 超出范围")
    if not hydrophobicity_pass:
        notes.append(f"连续强疏水残基 {hydrophobic_run} > {max_hydro}")
    if not oxidation_pass:
        notes.append("W 或 M 数量超限")
    if low_complexity_flag:
        notes.append("低复杂度序列")
    if poly_ionic_flag:
        notes.append("poly-basic/poly-acidic 风险")

    sequence_filter_pass = all(
        [
            starts_with_cys,
            ends_with_cys,
            internal_cys_count == 0,
            core_length_pass,
            charge_pass,
            hydrophobicity_pass,
            oxidation_pass,
            not low_complexity_flag,
            not poly_ionic_flag,
        ]
    )

    return {
        "core_sequence": seq,
        "starts_with_cys": starts_with_cys,
        "ends_with_cys": ends_with_cys,
        "internal_cys_count": internal_cys_count,
        "core_length_pass": core_length_pass,
        "net_charge": charge,
        "charge_pass": charge_pass,
        "hydrophobic_run_max": hydrophobic_run,
        "hydrophobicity_pass": hydrophobicity_pass,
        "w_count": w_count,
        "m_count": m_count,
        "low_complexity_flag": low_complexity_flag,
        "sequence_filter_pass": sequence_filter_pass,
        "filter_notes": "; ".join(notes) if notes else "pass",
    }
