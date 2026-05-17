#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import math
import statistics
from collections import defaultdict
from pathlib import Path


PSEUDO = 0.01
COUNT_PSEUDO = 1.0
SIM_DETECTION_THRESHOLD = 1.0
SIM_SUPPLEMENT_DETECTION_THRESHOLDS = (2.0, 5.0)
SEQ_CONDS = ("A", "B", "C", "D")
SEQ_REPS = (1, 2, 3, 4)


def pearson(xs: list[float], ys: list[float]) -> float:
    if not xs:
        return math.nan
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    dx = sum((x - mx) ** 2 for x in xs)
    dy = sum((y - my) ** 2 for y in ys)
    if dx == 0.0 or dy == 0.0:
        return math.nan
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys, strict=True)) / math.sqrt(dx * dy)


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    result = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            result[order[k]] = rank
        i = j
    return result


def spearman(xs: list[float], ys: list[float]) -> float:
    return pearson(ranks(xs), ranks(ys))


def slope(xs: list[float], ys: list[float]) -> float:
    if not xs:
        return math.nan
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    den = sum((x - mx) ** 2 for x in xs)
    if den == 0.0:
        return math.nan
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys, strict=True)) / den


def sample_cv(values: list[float]) -> float:
    mean = sum(values) / len(values)
    return statistics.stdev(values) / (mean + 1e-6)


def short_id(target: str) -> str:
    return target.split()[0]


def parse_gene(target: str | None) -> str | None:
    if target is None:
        return None
    fields = target.split("|")
    return fields[5] if len(fields) >= 6 else None


def parse_fasta_classes(fasta: Path) -> tuple[dict[str, str], dict[str, str]]:
    target_to_class: dict[str, str] = {}
    class_to_representative: dict[str, str] = {}
    members_by_class: dict[str, list[str]] = defaultdict(list)
    current_name: str | None = None
    current_hash = hashlib.sha256()
    current_len = 0

    def flush() -> None:
        nonlocal current_name, current_hash, current_len
        if current_name is None:
            return
        members_by_class[f"{current_len}:{current_hash.hexdigest()}"].append(current_name)

    with fasta.open() as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                flush()
                current_name = short_id(line[1:])
                current_hash = hashlib.sha256()
                current_len = 0
            else:
                seq = line.upper().encode("ascii")
                current_hash.update(seq)
                current_len += len(seq)
        flush()

    for class_id, members in members_by_class.items():
        class_to_representative[class_id] = members[0]
        for member in members:
            target_to_class[member] = class_id
    return target_to_class, class_to_representative


def class_of(target_to_class: dict[str, str], target: str) -> str:
    return target_to_class.get(short_id(target), short_id(target))


def read_piscem(path: Path, target_to_class: dict[str, str] | None = None) -> dict[str, tuple[float, float]]:
    if not path.exists():
        return {}
    result: dict[str, tuple[float, float]] = {}
    with path.open() as handle:
        handle.readline()
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            target = class_of(target_to_class, parts[0]) if target_to_class else parts[0]
            old_tpm, old_count = result.get(target, (0.0, 0.0))
            result[target] = (old_tpm + float(parts[3]), old_count + float(parts[4]))
    return result


def read_salmon(path: Path, target_to_class: dict[str, str] | None = None) -> dict[str, tuple[float, float]]:
    if not path.exists():
        return {}
    result: dict[str, tuple[float, float]] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            target = class_of(target_to_class, row["Name"]) if target_to_class else row["Name"]
            old_tpm, old_count = result.get(target, (0.0, 0.0))
            result[target] = (old_tpm + float(row["TPM"]), old_count + float(row["NumReads"]))
    return result


def read_kallisto(path: Path, target_to_class: dict[str, str] | None = None) -> dict[str, tuple[float, float]]:
    if not path.exists():
        return {}
    result: dict[str, tuple[float, float]] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            target = class_of(target_to_class, row["target_id"]) if target_to_class else row["target_id"]
            old_tpm, old_count = result.get(target, (0.0, 0.0))
            result[target] = (old_tpm + float(row["tpm"]), old_count + float(row["est_counts"]))
    return result


class Method:
    def __init__(self, name: str, reader, path_fn):
        self.name = name
        self.reader = reader
        self.path_fn = path_fn

    def read(self, sample: str, target_to_class: dict[str, str] | None = None) -> dict[str, tuple[float, float]]:
        return self.reader(self.path_fn(sample), target_to_class)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def format_float(value: object, digits: int = 2) -> str:
    return f"{float(value):.{digits}f}"


def format_int(value: object) -> str:
    return str(round(float(value)))


def format_approx_k(value: object) -> str:
    return f"≈{round(float(value) / 1000):.0f}K"


def write_display_csv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def write_display_tables(outdir: Path, computed: dict[str, list[dict[str, object]]]) -> None:
    display = outdir / "display"

    write_display_csv(
        display / "gencode_simulation_count.csv",
        ["Method", "TP", "FP", "Prec", "Rec", "F1", "Pearson", "Spear", "FC slope"],
        [
            [
                str(row["Method"]),
                format_int(row["TP"]),
                format_int(row["FP"]),
                format_float(row["Precision"]),
                format_float(row["Recall"]),
                format_float(row["F1"]),
                format_float(row["Pearson"]),
                format_float(row["Spearman"]),
                format_float(row["FC slope"]),
            ]
            for row in computed["gencode_simulation_count"]
        ],
    )
    for threshold in SIM_SUPPLEMENT_DETECTION_THRESHOLDS:
        write_display_csv(
            display / f"gencode_simulation_count_threshold_{threshold:g}.csv",
            ["Method", "TP", "FP", "Prec", "Rec", "F1", "Pearson", "Spear", "FC slope"],
            [
                [
                    str(row["Method"]),
                    format_int(row["TP"]),
                    format_int(row["FP"]),
                    format_float(row["Precision"]),
                    format_float(row["Recall"]),
                    format_float(row["F1"]),
                    format_float(row["Pearson"]),
                    format_float(row["Spearman"]),
                    format_float(row["FC slope"]),
                ]
                for row in computed[f"gencode_simulation_count_threshold_{threshold:g}"]
            ],
        )
    write_display_csv(
        display / "airway_replicate_cv.csv",
        ["Method", "Med CV (UNION)", "Med CV (INTER)"],
        [
            [str(row["Method"]), format_float(row["Med CV (UNION)"]), format_float(row["Med CV (INTER)"])]
            for row in computed["airway_replicate_cv"]
        ],
    )
    write_display_csv(
        display / "seqc_titration_full.csv",
        ["Method", "Pe(C)", "Pe(D)", "Sp(C)", "Sp(D)", "FC slope", "Det (A)"],
        [
            [
                str(row["Method"]),
                format_float(row["Pe(C)"]),
                format_float(row["Pe(D)"]),
                format_float(row["Sp(C)"]),
                format_float(row["Sp(D)"]),
                format_float(row["FC slope"]),
                format_approx_k(row["Det(A)"]),
            ]
            for row in computed["seqc_titration_full"]
        ],
    )
    write_display_csv(
        display / "seqc_taqman.csv",
        ["Method", "Pearson(A)", "Spearman(A)", "FC Pearson", "Det A&B"],
        [
            [
                str(row["Method"]),
                format_float(row["Pearson(A)"]),
                format_float(row["Spearman(A)"]),
                format_float(row["FC Pearson"]),
                format_int(row["Det A&B"]),
            ]
            for row in computed["seqc_taqman"]
        ],
    )

    full = {str(row["Method"]): row for row in computed["seqc_titration_full"]}
    low = {str(row["Method"]): row for row in computed["seqc_titration_2m"]}
    write_display_csv(
        display / "seqc_titration_depth.csv",
        ["Method", "Pe(C) full", "Pe(C) 2M", "FC full", "FC 2M"],
        [
            [
                method,
                format_float(full[method]["Pe(C)"]),
                format_float(low[method]["Pe(C)"]),
                format_float(full[method]["FC slope"]),
                format_float(low[method]["FC slope"]),
            ]
            for method in full
        ],
    )

    full_cv = {str(row["Method"]): row for row in computed["seqc_cv_full"]}
    low_cv = {str(row["Method"]): row for row in computed["seqc_cv_2m"]}
    write_display_csv(
        display / "seqc_cv_depth.csv",
        ["Method", "Full UNION", "Full INTER", "2M UNION", "2M INTER"],
        [
            [
                method,
                format_float(full_cv[method]["UNION"]),
                format_float(full_cv[method]["INTER"]),
                format_float(low_cv[method]["UNION"]),
                format_float(low_cv[method]["INTER"]),
            ]
            for method in full_cv
        ],
    )


def mean_for_samples(method: Method, samples: list[str], target_to_class: dict[str, str]) -> dict[str, float]:
    sums: dict[str, float] = defaultdict(float)
    counts: dict[str, int] = defaultdict(int)
    for sample in samples:
        quant = method.read(sample, target_to_class)
        for target, (tpm, _count) in quant.items():
            sums[target] += tpm
            counts[target] += 1
    return {target: sums[target] / counts[target] for target in sums}


def quant_maps(method: Method, samples: list[str], target_to_class: dict[str, str]) -> dict[str, dict[str, tuple[float, float]]]:
    return {sample: method.read(sample, target_to_class) for sample in samples}


def eval_sim(
    simdir: Path,
    target_to_class: dict[str, str],
    detection_threshold: float = SIM_DETECTION_THRESHOLD,
) -> list[dict[str, object]]:
    methods = [
        Method("piscem-infer (single)", read_piscem, lambda s: simdir / "quant_em" / s / f"{s}.quant"),
        Method("salmon (VBEM)", read_salmon, lambda s: simdir / "quant_salmon" / s / "quant.sf"),
        Method("salmon (EM)", read_salmon, lambda s: simdir / "quant_salmon_em" / s / "quant.sf"),
        Method("kallisto", read_kallisto, lambda s: simdir / "quant_kallisto" / s / "abundance.tsv"),
        Method("pims (strict)", read_piscem, lambda s: simdir / "quant_consensus_sel_support" / s / f"{s}.quant"),
        Method(
            "pims (+ rescue)",
            read_piscem,
            lambda s: simdir / "rerun_pims_rescue_enrichment_credible_stability_075_cv1_mean5_20260516" / "rescue" / s / f"{s}.quant",
        ),
    ]
    gt: dict[str, dict[str, float | bool]] = {}
    with (simdir / "ground_truth.csv").open() as handle:
        for row in csv.DictReader(handle):
            class_id = class_of(target_to_class, row["transcript_id"])
            entry = gt.setdefault(
                class_id,
                {"control": 0.0, "treatment": 0.0, "is_de": False},
            )
            entry["control"] = float(entry["control"]) + float(row["expected_reads_control"])
            entry["treatment"] = float(entry["treatment"]) + float(row["expected_reads_treatment"])
            entry["is_de"] = bool(entry["is_de"]) or row.get("is_de", "").lower() in {"true", "1"}

    with (simdir / "sample_info.csv").open() as handle:
        samples = list(csv.DictReader(handle))

    rows = []
    for method in methods:
        all_p = []
        all_s = []
        all_rmse = []
        tp = fp = fn = 0
        sample_counts: dict[str, dict[str, float]] = {}
        for sample in samples:
            sn = sample["sample_name"]
            cond = sample["condition"]
            quant = method.read(sn, target_to_class)
            sample_counts[sn] = {target: count for target, (_tpm, count) in quant.items()}
            est = []
            true = []
            for target in set(gt).union(quant):
                est_count = quant.get(target, (0.0, 0.0))[1]
                true_count = float(gt.get(target, {"control": 0.0, "treatment": 0.0})[cond])
                est.append(est_count)
                true.append(true_count)
                pred_detected = est_count > detection_threshold
                true_detected = true_count > detection_threshold
                tp += pred_detected and true_detected
                fp += pred_detected and not true_detected
                fn += not pred_detected and true_detected
            log_est = [math.log2(v + COUNT_PSEUDO) for v in est]
            log_true = [math.log2(v + COUNT_PSEUDO) for v in true]
            all_p.append(pearson(log_est, log_true))
            all_s.append(spearman(est, true))
            all_rmse.append(math.sqrt(sum((e - t) ** 2 for e, t in zip(log_est, log_true, strict=True)) / len(log_est)))

        ctrl = [sample["sample_name"] for sample in samples if sample["condition"] == "control"]
        trt = [sample["sample_name"] for sample in samples if sample["condition"] == "treatment"]
        est_lfc = []
        true_lfc = []
        for target, info in gt.items():
            truth_lfc = math.log2(float(info["treatment"]) + COUNT_PSEUDO) - math.log2(float(info["control"]) + COUNT_PSEUDO)
            if not bool(info["is_de"]) or abs(truth_lfc) <= 1.0:
                continue
            ctrl_count = sum(sample_counts[s].get(target, 0.0) for s in ctrl) / len(ctrl)
            trt_count = sum(sample_counts[s].get(target, 0.0) for s in trt) / len(trt)
            est_lfc.append(math.log2(trt_count + COUNT_PSEUDO) - math.log2(ctrl_count + COUNT_PSEUDO))
            true_lfc.append(truth_lfc)

        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        f1 = 2 * precision * recall / (precision + recall)
        rows.append(
            {
                "Method": method.name,
                "Detection threshold": detection_threshold,
                "TP": tp / len(samples),
                "FP": fp / len(samples),
                "FN": fn / len(samples),
                "Precision": precision,
                "Recall": recall,
                "F1": f1,
                "Pearson": sum(all_p) / len(all_p),
                "Spearman": sum(all_s) / len(all_s),
                "RMSE": sum(all_rmse) / len(all_rmse),
                "FC slope": slope(true_lfc, est_lfc),
            }
        )
    return rows


def eval_airway(benchdir: Path, target_to_class: dict[str, str]) -> list[dict[str, object]]:
    samples = ["SRR1039508", "SRR1039512", "SRR1039516", "SRR1039520"]
    methods = [
        Method("piscem-infer (single)", read_piscem, lambda s: benchdir / "quant/em" / f"{s}.quant"),
        Method("salmon (VBEM)", read_salmon, lambda s: benchdir / "quant/salmon" / s / "quant.sf"),
        Method("salmon (EM)", read_salmon, lambda s: benchdir / "quant/salmon_em" / s / "quant.sf"),
        Method("kallisto", read_kallisto, lambda s: benchdir / "quant/kallisto" / s / "abundance.tsv"),
        Method("pims (strict)", read_piscem, lambda s: benchdir / "quant/sel_support_all8_adaptive" / s / f"{s}.quant"),
        Method(
            "pims (+ rescue)",
            read_piscem,
            lambda s: benchdir / "quant/putative_locked_rescue_enrichment_credible_stability_075_cv1_mean5_20260516" / s / f"{s}.quant",
        ),
    ]
    means = {method.name: mean_for_samples(method, samples, target_to_class) for method in methods}
    all_targets = sorted({target for method_means in means.values() for target in method_means})
    union = [target for target in all_targets if any(method_means.get(target, 0.0) >= 1.0 for method_means in means.values())]
    inter = [target for target in all_targets if all(method_means.get(target, 0.0) >= 1.0 for method_means in means.values())]
    rows = []
    for method in methods:
        maps = quant_maps(method, samples, target_to_class)
        rows.append(
            {
                "Method": method.name,
                "UnionN": len(union),
                "InterN": len(inter),
                "Med CV (UNION)": statistics.median(sample_cv([maps[s].get(target, (0.0, 0.0))[0] for s in samples]) for target in union),
                "Med CV (INTER)": statistics.median(sample_cv([maps[s].get(target, (0.0, 0.0))[0] for s in samples]) for target in inter),
            }
        )
    return rows


def seqc_samples() -> list[str]:
    return [f"{cond}_{rep}" for cond in SEQ_CONDS for rep in SEQ_REPS]


def seqc_condition_samples(cond: str) -> list[str]:
    return [f"{cond}_{rep}" for rep in SEQ_REPS]


def seqc_methods(benchdir: Path, suffix: str = "") -> list[Method]:
    if suffix == "_2M":
        return [
            Method("piscem-infer (single)", read_piscem, lambda s: benchdir / "quant_piscem_single_2M" / s / f"{s}.quant"),
            Method("salmon (VBEM)", read_salmon, lambda s: benchdir / "quant_salmon_2M" / s / "quant.sf"),
            Method("salmon (EM)", read_salmon, lambda s: benchdir / "quant_salmon_em_2M" / s / "quant.sf"),
            Method("kallisto", read_kallisto, lambda s: benchdir / "quant_kallisto_2M" / s / "abundance.tsv"),
            Method("pims (strict)", read_piscem, lambda s: benchdir / "quant_2M_strict" / s / f"{s}.quant"),
            Method(
                "pims (+ rescue)",
                read_piscem,
                lambda s: benchdir / "quant_2M_putative_locked_rescue_enrichment_credible_stability_20260516" / s / f"{s}.quant",
            ),
        ]
    return [
        Method("piscem-infer (single)", read_piscem, lambda s: benchdir / "quant_piscem_single" / s / f"{s}.quant"),
        Method("salmon (VBEM)", read_salmon, lambda s: benchdir / "quant_salmon" / s / "quant.sf"),
        Method("salmon (EM)", read_salmon, lambda s: benchdir / "quant_salmon_em" / s / "quant.sf"),
        Method("kallisto", read_kallisto, lambda s: benchdir / "quant_kallisto" / s / "abundance.tsv"),
        Method("pims (strict)", read_piscem, lambda s: benchdir / "quant_sel_adapt_ecgraph" / s / f"{s}.quant"),
        Method(
            "pims (+ rescue)",
            read_piscem,
            lambda s: benchdir / "quant_putative_locked_rescue_enrichment_credible_stability_075_cv1_mean5_20260516" / s / f"{s}.quant",
        ),
    ]


def seqc_means(method: Method, target_to_class: dict[str, str]) -> dict[str, dict[str, float]]:
    return {cond: mean_for_samples(method, seqc_condition_samples(cond), target_to_class) for cond in SEQ_CONDS}


def seqc_ground_sets(methods: list[Method], target_to_class: dict[str, str]) -> tuple[list[str], list[str]]:
    means = {method.name: seqc_means(method, target_to_class) for method in methods}
    all_targets = sorted({t for method_means in means.values() for cond in ("A", "B") for t in method_means[cond]})
    union = [
        target
        for target in all_targets
        if any(means[method.name]["A"].get(target, 0.0) >= 1.0 or means[method.name]["B"].get(target, 0.0) >= 1.0 for method in methods)
    ]
    inter = [
        target
        for target in all_targets
        if all(means[method.name]["A"].get(target, 0.0) >= 1.0 and means[method.name]["B"].get(target, 0.0) >= 1.0 for method in methods)
    ]
    return union, inter


def eval_seqc_titration(benchdir: Path, target_to_class: dict[str, str], suffix: str = "") -> list[dict[str, object]]:
    ground_methods = seqc_methods(benchdir)
    union, _inter = seqc_ground_sets(ground_methods, target_to_class)
    rows = []
    for method in seqc_methods(benchdir, suffix):
        means = seqc_means(method, target_to_class)
        a = [means["A"].get(t, 0.0) for t in union]
        b = [means["B"].get(t, 0.0) for t in union]
        c = [means["C"].get(t, 0.0) for t in union]
        d = [means["D"].get(t, 0.0) for t in union]
        c_exp = [0.75 * av + 0.25 * bv for av, bv in zip(a, b, strict=True)]
        d_exp = [0.25 * av + 0.75 * bv for av, bv in zip(a, b, strict=True)]
        obs_lfc = [math.log2(cv + PSEUDO) - math.log2(dv + PSEUDO) for cv, dv in zip(c, d, strict=True)]
        exp_lfc = [math.log2(cv + PSEUDO) - math.log2(dv + PSEUDO) for cv, dv in zip(c_exp, d_exp, strict=True)]
        big = [idx for idx, value in enumerate(exp_lfc) if abs(value) > 0.5]
        rows.append(
            {
                "Method": method.name,
                "GroundN": len(union),
                "Pe(C)": pearson([math.log2(v + PSEUDO) for v in c], [math.log2(v + PSEUDO) for v in c_exp]),
                "Pe(D)": pearson([math.log2(v + PSEUDO) for v in d], [math.log2(v + PSEUDO) for v in d_exp]),
                "Sp(C)": spearman(c, c_exp),
                "Sp(D)": spearman(d, d_exp),
                "FC slope": slope([exp_lfc[i] for i in big], [obs_lfc[i] for i in big]),
                "Det(A)": sum(value > 0 for value in means["A"].values()),
            }
        )
    return rows


def eval_seqc_cv(benchdir: Path, target_to_class: dict[str, str], suffix: str = "") -> list[dict[str, object]]:
    ground_methods = seqc_methods(benchdir)
    union, inter = seqc_ground_sets(ground_methods, target_to_class)
    rows = []
    for method in seqc_methods(benchdir, suffix):
        maps = quant_maps(method, seqc_samples(), target_to_class)
        union_cvs = []
        inter_cvs = []
        for cond in SEQ_CONDS:
            samples = seqc_condition_samples(cond)
            union_cvs.append(statistics.median(sample_cv([maps[s].get(target, (0.0, 0.0))[0] for s in samples]) for target in union))
            inter_cvs.append(statistics.median(sample_cv([maps[s].get(target, (0.0, 0.0))[0] for s in samples]) for target in inter))
        rows.append(
            {
                "Method": method.name,
                "UnionN": len(union),
                "InterN": len(inter),
                "UNION": sum(union_cvs) / len(union_cvs),
                "INTER": sum(inter_cvs) / len(inter_cvs),
            }
        )
    return rows


def load_taqman(benchdir: Path) -> list[tuple[str, float, float]]:
    id_to_gene: dict[str, str] = {}
    with (benchdir / "GPL4097.annot").open(errors="ignore") as handle:
        header_seen = False
        for line in handle:
            if line.startswith(("#", "!", "^")) or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if not header_seen:
                header_seen = True
                continue
            if len(parts) >= 3 and parts[2] and parts[2] != "---":
                id_to_gene[parts[0].strip('"')] = parts[2]

    lines = (benchdir / "taqman_raw.txt").read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith('"ID_REF"'))
    end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("!")), len(lines))
    reader = csv.reader(lines[start:end], delimiter="\t")
    header = next(reader)
    labels = [f"{cond}_{rep}" for cond in SEQ_CONDS for rep in SEQ_REPS]
    rows = []
    seen = set()
    for row in reader:
        taq_id = row[0].strip('"')
        gene = id_to_gene.get(taq_id)
        if not gene or gene in seen:
            continue
        values = [float(value) if value else 0.0 for value in row[1 : 1 + len(labels)]]
        by_sample = dict(zip(labels, values, strict=True))
        rows.append(
            (
                gene,
                sum(by_sample[f"A_{rep}"] for rep in SEQ_REPS) / len(SEQ_REPS),
                sum(by_sample[f"B_{rep}"] for rep in SEQ_REPS) / len(SEQ_REPS),
            )
        )
        seen.add(gene)
    return rows


def gene_means(
    method: Method,
    cond: str,
    target_to_class: dict[str, str],
    class_to_representative: dict[str, str],
) -> dict[str, float]:
    sums: dict[str, float] = defaultdict(float)
    counts: dict[str, int] = defaultdict(int)
    for sample in seqc_condition_samples(cond):
        gene_tpm: dict[str, float] = defaultdict(float)
        quant = method.read(sample, target_to_class)
        for target, (tpm, _count) in quant.items():
            gene = parse_gene(class_to_representative.get(target, target))
            if gene is not None:
                gene_tpm[gene] += tpm
        for gene, tpm in gene_tpm.items():
            sums[gene] += tpm
            counts[gene] += 1
    return {gene: sums[gene] / counts[gene] for gene in sums}


def eval_taqman(
    benchdir: Path,
    target_to_class: dict[str, str],
    class_to_representative: dict[str, str],
) -> list[dict[str, object]]:
    taq = load_taqman(benchdir)
    ground_a = [gene for gene, a, _b in taq if a > 0.001]
    ground_b = [gene for gene, _a, b in taq if b > 0.001]
    ground_ab = [gene for gene in ground_a if gene in set(ground_b)]
    taq_a = {gene: a for gene, a, _b in taq}
    taq_b = {gene: b for gene, _a, b in taq}
    rows = []
    for method in seqc_methods(benchdir):
        a_means = gene_means(method, "A", target_to_class, class_to_representative)
        b_means = gene_means(method, "B", target_to_class, class_to_representative)
        a_values = [a_means.get(gene, 0.0) for gene in ground_a]
        b_values = [b_means.get(gene, 0.0) for gene in ground_b]
        ab_a = [a_means.get(gene, 0.0) for gene in ground_ab]
        ab_b = [b_means.get(gene, 0.0) for gene in ground_ab]
        rna_lfc = [math.log2(a + PSEUDO) - math.log2(b + PSEUDO) for a, b in zip(ab_a, ab_b, strict=True)]
        taq_lfc = [math.log2(taq_a[gene] + PSEUDO) - math.log2(taq_b[gene] + PSEUDO) for gene in ground_ab]
        rows.append(
            {
                "Method": method.name,
                "Pearson(A)": pearson(
                    [math.log2(value + PSEUDO) for value in a_values],
                    [math.log2(taq_a[gene] + PSEUDO) for gene in ground_a],
                ),
                "Spearman(A)": spearman(a_values, [taq_a[gene] for gene in ground_a]),
                "Pearson(B)": pearson(
                    [math.log2(value + PSEUDO) for value in b_values],
                    [math.log2(taq_b[gene] + PSEUDO) for gene in ground_b],
                ),
                "FC Pearson": pearson(rna_lfc, taq_lfc),
                "Det A&B": sum(a > 0 and b > 0 for a, b in zip(ab_a, ab_b, strict=True)),
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--simdir", type=Path, default=Path("sim_data_gencode"))
    parser.add_argument("--airwaydir", type=Path, default=Path("airway_benchmark"))
    parser.add_argument("--seqcdir", type=Path, default=Path("seqc_benchmark"))
    parser.add_argument("--fasta", type=Path, default=Path("sim_data_gencode/gencode.v49.pc_transcripts.fa"))
    parser.add_argument("--outdir", type=Path, default=Path("paper_results_export"))
    args = parser.parse_args()

    target_to_class, class_to_representative = parse_fasta_classes(args.fasta)
    computed = {
        "gencode_simulation_count": eval_sim(args.simdir, target_to_class),
        **{
            f"gencode_simulation_count_threshold_{threshold:g}": eval_sim(
                args.simdir,
                target_to_class,
                detection_threshold=threshold,
            )
            for threshold in SIM_SUPPLEMENT_DETECTION_THRESHOLDS
        },
        "airway_replicate_cv": eval_airway(args.airwaydir, target_to_class),
        "seqc_titration_full": eval_seqc_titration(args.seqcdir, target_to_class),
        "seqc_titration_2m": eval_seqc_titration(args.seqcdir, target_to_class, "_2M"),
        "seqc_cv_full": eval_seqc_cv(args.seqcdir, target_to_class),
        "seqc_cv_2m": eval_seqc_cv(args.seqcdir, target_to_class, "_2M"),
        "seqc_taqman": eval_taqman(args.seqcdir, target_to_class, class_to_representative),
    }
    for name, rows in computed.items():
        write_csv(args.outdir / f"{name}.csv", rows)
    write_display_tables(args.outdir, computed)


if __name__ == "__main__":
    main()
