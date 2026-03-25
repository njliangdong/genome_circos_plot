#!/usr/bin/env python3
import argparse
import csv
import json
import math
from collections import defaultdict, Counter

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Circle

try:
    import pandas as pd
except Exception:  # pragma: no cover
    pd = None


def read_fasta(path):
    contigs = []
    name = None
    seq_parts = []
    chrom_idx = 0
    scaffold_idx = 0
    contig_idx = 0
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seq = "".join(seq_parts)
                    display = name
                    lname = name.lower()
                    if "chromosome" in lname:
                        chrom_idx += 1
                        display = f"C{chrom_idx}"
                    elif "scaffold" in lname:
                        scaffold_idx += 1
                        display = f"S{scaffold_idx}"
                    elif "contig" in lname:
                        contig_idx += 1
                        display = f"c{contig_idx}"
                    contigs.append({"name": name, "display_name": display, "seq": seq, "length": len(seq)})
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
        if name is not None:
            seq = "".join(seq_parts)
            display = name
            lname = name.lower()
            if "chromosome" in lname:
                chrom_idx += 1
                display = f"C{chrom_idx}"
            elif "scaffold" in lname:
                scaffold_idx += 1
                display = f"S{scaffold_idx}"
            elif "contig" in lname:
                contig_idx += 1
                display = f"c{contig_idx}"
            contigs.append({"name": name, "display_name": display, "seq": seq, "length": len(seq)})
    return contigs


def parse_gff3(path):
    features = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            try:
                start = int(start)
                end = int(end)
            except ValueError:
                continue
            attr_map = {}
            for item in attrs.split(";"):
                if not item:
                    continue
                if "=" in item:
                    k, v = item.split("=", 1)
                    attr_map[k] = v
            fid = attr_map.get("ID", "")
            if "," in fid:
                fid = fid.split(",", 1)[0]
            features.append({
                "seqid": seqid,
                "type": ftype,
                "start": start,
                "end": end,
                "strand": strand,
                "attrs": attr_map,
                "id": fid,
            })
    return features


def build_contig_offsets(contigs, gap):
    offset = 0
    for c in contigs:
        c["offset"] = offset
        offset += c["length"] + gap
    total = offset - gap if contigs else 0
    return total


def windows_for_contig(length, window):
    for i in range(0, length, window):
        yield i, min(length, i + window)


def window_midpoints(length, window):
    for s, e in windows_for_contig(length, window):
        yield (s + e) / 2


def wave_scale(values):
    if not values:
        return []
    vmin = min(values)
    vmax = max(values)
    if vmax == vmin:
        return [0.0 for _ in values]
    mid = (vmax + vmin) / 2.0
    half = (vmax - vmin) / 2.0
    return [(v - mid) / half for v in values]


def plot_wave_track(ax, contigs, total_len, values_by_contig, window, base, amp, color, linewidth=1.1):
    for c in contigs:
        vals = values_by_contig.get(c["name"], [])
        if not vals:
            continue
        scaled = wave_scale(vals)
        thetas = []
        radii = []
        for mid, sc in zip(window_midpoints(c["length"], window), scaled):
            pos = c["offset"] + mid
            thetas.append(to_angle(pos, total_len))
            radii.append(base + amp * sc)
        ax.plot(thetas, radii, color=color, linewidth=linewidth, solid_joinstyle="round")


def gc_content(seq):
    if not seq:
        return 0.0
    s = seq.upper()
    gc = s.count("G") + s.count("C")
    return (gc / len(s)) * 100.0


def to_angle(pos, total):
    return (pos / total) * 2 * math.pi


def parse_table(path, sep):
    with open(path) as f:
        reader = csv.DictReader(f, delimiter=sep)
        return list(reader)


def auto_window(total_len):
    target_bins = 800
    base = max(1000.0, total_len / target_bins)
    if base < 2000:
        step = 500
    elif base < 10000:
        step = 1000
    elif base < 50000:
        step = 5000
    else:
        step = 10000
    return int(math.ceil(base / step) * step)


def auto_gap(total_len, n_contigs):
    if n_contigs <= 1:
        return 0
    avg = total_len / max(1, n_contigs)
    return int(max(200, min(5000, avg * 0.01)))


def label_contigs(
    ax,
    contigs,
    total_len,
    radius,
    label_max,
    label_min_len,
    fontsize=13,
    color="#8b5a2b",
    style="tangent",
    angle_offset_deg=0.0,
    keep_upright=True,
):
    if not contigs:
        return
    candidates = [c for c in contigs if c["length"] >= label_min_len]
    if len(candidates) > label_max:
        candidates = sorted(candidates, key=lambda x: x["length"], reverse=True)[:label_max]
    for c in candidates:
        mid = c["offset"] + (c["length"] / 2)
        theta = to_angle(mid, total_len) + math.radians(angle_offset_deg)
        angle_deg = math.degrees(theta)
        if style == "radial":
            rotation = angle_deg
        else:
            # tangent to the circle (parallel to ring); with theta_direction=-1 this is -angle_deg
            rotation = -angle_deg
        if keep_upright:
            rotation = (rotation + 360) % 360
            if rotation > 180:
                rotation -= 360
            if rotation > 90:
                rotation -= 180
            elif rotation < -90:
                rotation += 180
        ha = "center"
        label_text = c.get("display_name", c["name"])
        ax.text(
            theta,
            radius,
            label_text,
            fontsize=fontsize,
            fontweight="bold",
            rotation=rotation,
            rotation_mode="anchor",
            ha=ha,
            va="center",
            color=color,
        )


def add_track_label(ax, radius, text, color="#8b5a2b", angle=0.0, fontsize=10):
    # place label tangent to circle
    angle_deg = math.degrees(angle)
    rotation = angle_deg - 90
    if rotation > 90:
        rotation -= 180
    elif rotation < -90:
        rotation += 180
    ax.text(
        angle,
        radius,
        text,
        fontsize=fontsize,
        fontweight="bold",
        rotation=rotation,
        rotation_mode="anchor",
        ha="center",
        va="center",
        color=color,
    )


def add_track_label_radial(ax, radius, text, color="#8b5a2b", angle=0.0, offset=0.08, fontsize=10):
    # place label radially (like reference image)
    ax.plot([angle, angle], [radius, radius + offset], color=color, lw=0.9)
    ax.text(
        angle,
        radius + offset * 1.2,
        text,
        fontsize=fontsize,
        fontweight="bold",
        rotation=0,
        rotation_mode="anchor",
        ha="left",
        va="center",
        color=color,
    )


def add_bgc_arrow(ax, theta_start, theta_end, r0, r1, strand, color):
    if theta_end <= theta_start:
        return
    width = theta_end - theta_start
    tip = min(width * 0.35, 0.08)
    if strand == "+":
        body_end = max(theta_start, theta_end - tip)
        points = [
            (theta_start, r0),
            (theta_start, r1),
            (body_end, r1),
            (theta_end, (r0 + r1) / 2),
            (body_end, r0),
        ]
    else:
        body_start = min(theta_end, theta_start + tip)
        points = [
            (body_start, r1),
            (theta_end, r1),
            (theta_end, r0),
            (body_start, r0),
            (theta_start, (r0 + r1) / 2),
        ]
    poly = plt.Polygon(points, closed=True, facecolor=color, edgecolor="none", transform=ax.transData)
    ax.add_patch(poly)


def plot_heatmap_ring(ax, contigs, total_len, window, values_by_contig, base, height, cmap, vmin, vmax, empty_color="#e5e7eb"):
    def norm(v):
        if vmax == vmin:
            return 0.5
        return (v - vmin) / (vmax - vmin)

    for c in contigs:
        vals = values_by_contig.get(c["name"], [])
        if not vals:
            continue
        for idx, (s, e) in enumerate(windows_for_contig(c["length"], window)):
            val = vals[idx] if idx < len(vals) else None
            gstart = c["offset"] + s
            gend = c["offset"] + e
            theta1 = to_angle(gstart, total_len)
            theta2 = to_angle(gend, total_len)
            if val is None:
                color = empty_color
            else:
                color = cmap(norm(val))
            ax.bar(
                x=(theta1 + theta2) / 2,
                height=height,
                width=(theta2 - theta1),
                bottom=base,
                color=color,
                edgecolor=None,
                linewidth=0,
            )


def build_category_colors(categories):
    palette = []
    for cmap_name in ("tab20", "tab20b", "tab20c"):
        cmap = plt.get_cmap(cmap_name)
        palette.extend([cmap(i) for i in range(cmap.N)])
    colors = {}
    for idx, cat in enumerate(categories):
        colors[cat] = palette[idx % len(palette)]
    return colors


def load_layout(path):
    with open(path) as f:
        data = json.load(f)
    if not isinstance(data, dict) or "tracks" not in data:
        raise ValueError("layout json must be a dict with a 'tracks' list")
    if not isinstance(data["tracks"], list):
        raise ValueError("'tracks' must be a list")
    return data


def main():
    p = argparse.ArgumentParser(
        description="Genome circos plot for bacteria/fungi with optional annotations",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "Required: --fasta, --gff3. "
            "Optional: --eggnog, --expr, --bgc, --splice (see flags for column mappings)."
        ),
    )
    p.add_argument("--fasta", required=True, help="Genome assembly FASTA")
    p.add_argument("--gff3", required=True, help="Genome annotation GFF3")
    p.add_argument("--window", type=int, default=None, help="Window size for tracks (auto if omitted)")
    p.add_argument("--gap", type=int, default=None, help="Gap size between contigs (auto if omitted)")
    p.add_argument("--order", choices=["input", "length"], default="input", help="Contig order")
    p.add_argument("--gene-feature", default="gene", help="GFF3 feature type for gene list")
    p.add_argument("--label-contigs", action="store_true", default=True, help="Add contig labels (default on)")
    p.add_argument("--no-label-contigs", action="store_false", dest="label_contigs", help="Disable contig labels")
    p.add_argument("--label-max", type=int, default=50, help="Max contig labels to draw")
    p.add_argument("--label-min-len", type=int, default=0, help="Only label contigs >= this length")
    p.add_argument("--layout", default=None, help="Layout JSON to control track order/positions")
    p.add_argument("--conf", default=None, help="Alias of --layout (track/label config JSON)")

    p.add_argument("--repeat-top", type=int, default=3, help="Repeat types to plot (top N by count)")
    p.add_argument("--repeat-types", default=None, help="Comma-separated repeat types to force plot")

    p.add_argument("--eggnog", default=None, help="EggNOG annotation file (optional)")
    p.add_argument("--eggnog-sep", default="\t", help="EggNOG separator")
    p.add_argument("--eggnog-id-col", default="query", help="EggNOG gene id column")
    p.add_argument("--eggnog-cog-col", default="COG_category", help="EggNOG COG column")

    p.add_argument("--expr", default=None, help="Expression table (long format)")
    p.add_argument("--expr-sep", default="\t", help="Expression separator")
    p.add_argument("--expr-id-col", default="gene_id", help="Expression gene id column")
    p.add_argument("--expr-cond-col", default="condition", help="Expression condition column")
    p.add_argument("--expr-log2fc-col", default="log2fc", help="Expression log2fc column")
    p.add_argument("--expr-padj-col", default="padj", help="Expression padj column")
    p.add_argument("--expr-log2fc-th", type=float, default=1.0, help="log2FC threshold")
    p.add_argument("--expr-padj-th", type=float, default=0.05, help="padj threshold")

    p.add_argument("--bgc", default=None, help="BGC regions file (GECCO/antiSMASH/bed-like)")
    p.add_argument("--bgc-sep", default="\t", help="BGC separator")
    p.add_argument("--bgc-contig-col", default="sequence_id", help="BGC contig column")
    p.add_argument("--bgc-start-col", default="start", help="BGC start column")
    p.add_argument("--bgc-end-col", default="end", help="BGC end column")
    p.add_argument("--bgc-type-col", default="type", help="BGC type column")
    p.add_argument("--bgc-strand-col", default="strand", help="BGC strand column (optional)")

    p.add_argument("--splice", default=None, help="Alternative splicing file (bed-like)")
    p.add_argument("--splice-sep", default="\t", help="Splice separator")
    p.add_argument("--splice-contig-col", default="seqid", help="Splice contig column")
    p.add_argument("--splice-start-col", default="start", help="Splice start column")
    p.add_argument("--splice-end-col", default="end", help="Splice end column")
    p.add_argument("--splice-type-col", default="event", help="Splice event column")

    p.add_argument("--output", default="genome_circos.png", help="Output image")
    p.add_argument("--excel", default="genome_circos.xlsx", help="Output Excel")

    args = p.parse_args()

    contigs = read_fasta(args.fasta)
    if args.order == "length":
        contigs.sort(key=lambda x: x["length"], reverse=True)

    raw_total_len = sum(c["length"] for c in contigs)
    if raw_total_len <= 0:
        raise SystemExit("No contigs found in FASTA")

    window = args.window if args.window is not None else auto_window(raw_total_len)
    gap = args.gap if args.gap is not None else auto_gap(raw_total_len, len(contigs))

    total_len = build_contig_offsets(contigs, gap)

    features = parse_gff3(args.gff3)
    feat_by_type = defaultdict(list)
    for f in features:
        feat_by_type[f["type"]].append(f)

    gene_feature = args.gene_feature if args.gene_feature in feat_by_type else "CDS"
    genes = feat_by_type.get(gene_feature, [])

    # ---- GC content & gene density & repeat density ----
    gc_rows = []
    gene_density_rows = []
    repeat_density_rows = []

    gene_mids_by_contig = defaultdict(list)
    for g in genes:
        mid = (g["start"] + g["end"]) // 2
        gene_mids_by_contig[g["seqid"]].append(mid)

    repeat_types_all = [t for t in feat_by_type if "repeat" in t.lower() or "transpos" in t.lower()]
    repeat_mids_by_contig = defaultdict(list)
    for t in repeat_types_all:
        for f in feat_by_type[t]:
            mid = (f["start"] + f["end"]) // 2
            repeat_mids_by_contig[f["seqid"]].append(mid)

    gc_by_contig = {}
    gene_density_by_contig = {}
    repeat_density_by_contig = {}

    for c in contigs:
        seq = c["seq"]
        nwin = (c["length"] + window - 1) // window
        gene_counts = [0] * nwin
        repeat_counts = [0] * nwin
        for mid in gene_mids_by_contig.get(c["name"], []):
            idx = min(mid // window, nwin - 1)
            gene_counts[idx] += 1
        for mid in repeat_mids_by_contig.get(c["name"], []):
            idx = min(mid // window, nwin - 1)
            repeat_counts[idx] += 1

        gc_vals = []
        gene_vals = []
        repeat_vals = []
        for idx, (s, e) in enumerate(windows_for_contig(c["length"], window)):
            win_seq = seq[s:e]
            gc = gc_content(win_seq)
            gcount = gene_counts[idx] if idx < len(gene_counts) else 0
            rcount = repeat_counts[idx] if idx < len(repeat_counts) else 0
            gc_vals.append(gc)
            gene_vals.append(gcount)
            repeat_vals.append(rcount)
            gc_rows.append({
                "contig": c["name"],
                "start": s + 1,
                "end": e,
                "gc": gc,
            })
            gene_density_rows.append({
                "contig": c["name"],
                "start": s + 1,
                "end": e,
                "gene_count": gcount,
            })
            repeat_density_rows.append({
                "contig": c["name"],
                "start": s + 1,
                "end": e,
                "repeat_count": rcount,
            })
        gc_by_contig[c["name"]] = gc_vals
        gene_density_by_contig[c["name"]] = gene_vals
        repeat_density_by_contig[c["name"]] = repeat_vals

    # ---- Repeat distribution ----
    repeat_types = repeat_types_all
    if args.repeat_types:
        top_repeat_types = [t.strip() for t in args.repeat_types.split(",") if t.strip()]
    else:
        repeat_counts = Counter({t: len(feat_by_type[t]) for t in repeat_types})
        top_repeat_types = [t for t, _ in repeat_counts.most_common(args.repeat_top)]

    repeat_rows = []
    for t in top_repeat_types:
        for f in feat_by_type[t]:
            repeat_rows.append({
                "type": t,
                "contig": f["seqid"],
                "start": f["start"],
                "end": f["end"],
            })

    # ---- EggNOG annotation ----
    eggnog_map = {}
    eggnog_counts = Counter()
    if args.eggnog:
        for row in parse_table(args.eggnog, args.eggnog_sep):
            gid = row.get(args.eggnog_id_col, "").strip()
            cog = row.get(args.eggnog_cog_col, "").strip()
            if gid and cog:
                eggnog_map[gid] = cog
                eggnog_counts[cog] += 1

    # ---- Expression ----
    expr_map = defaultdict(list)
    if args.expr:
        for row in parse_table(args.expr, args.expr_sep):
            gid = row.get(args.expr_id_col, "").strip()
            cond = row.get(args.expr_cond_col, "").strip()
            if not gid or not cond:
                continue
            try:
                log2fc = float(row.get(args.expr_log2fc_col, ""))
            except ValueError:
                log2fc = None
            try:
                padj = float(row.get(args.expr_padj_col, ""))
            except ValueError:
                padj = None
            status = "ns"
            if log2fc is not None and padj is not None:
                if padj <= args.expr_padj_th and log2fc >= args.expr_log2fc_th:
                    status = "up"
                elif padj <= args.expr_padj_th and log2fc <= -args.expr_log2fc_th:
                    status = "down"
            expr_map[gid].append({"condition": cond, "log2fc": log2fc, "padj": padj, "status": status})

    expr_conditions = sorted({e["condition"] for items in expr_map.values() for e in items})
    expr_window_values = {cond: {} for cond in expr_conditions}
    expr_max_abs = 0.0
    if expr_conditions:
        # build per-window mean log2fc per condition
        gene_info = []
        for g in genes:
            gid = g["id"] or g["attrs"].get("ID", "")
            if not gid:
                continue
            mid = (g["start"] + g["end"]) // 2
            gene_info.append((g["seqid"], mid, gid))

        gene_by_contig = defaultdict(list)
        for seqid, mid, gid in gene_info:
            gene_by_contig[seqid].append((mid, gid))

        for c in contigs:
            nwin = (c["length"] + window - 1) // window
            for cond in expr_conditions:
                buckets = [[] for _ in range(nwin)]
                for mid, gid in gene_by_contig.get(c["name"], []):
                    idx = min(mid // window, nwin - 1)
                    for e in expr_map.get(gid, []):
                        if e["condition"] == cond and e["log2fc"] is not None:
                            buckets[idx].append(e["log2fc"])
                vals = []
                for b in buckets:
                    if b:
                        m = sum(b) / len(b)
                        vals.append(m)
                        expr_max_abs = max(expr_max_abs, abs(m))
                    else:
                        vals.append(None)
                expr_window_values[cond][c["name"]] = vals

    # ---- BGC regions ----
    bgc_rows = []
    if args.bgc:
        for row in parse_table(args.bgc, args.bgc_sep):
            contig = row.get(args.bgc_contig_col, "")
            start = row.get(args.bgc_start_col, "")
            end = row.get(args.bgc_end_col, "")
            btype = row.get(args.bgc_type_col, "")
            strand = row.get(args.bgc_strand_col, "+") or "+"
            if strand not in ("+", "-"):
                strand = "+"
            if contig and start and end:
                bgc_rows.append({
                    "contig": contig,
                    "start": int(float(start)),
                    "end": int(float(end)),
                    "type": btype or "BGC",
                    "strand": strand,
                })

    # ---- Splicing ----
    splice_rows = []
    if args.splice:
        for row in parse_table(args.splice, args.splice_sep):
            contig = row.get(args.splice_contig_col, "")
            start = row.get(args.splice_start_col, "")
            end = row.get(args.splice_end_col, "")
            stype = row.get(args.splice_type_col, "")
            if contig and start and end:
                splice_rows.append({
                    "contig": contig,
                    "start": int(float(start)),
                    "end": int(float(end)),
                    "type": stype or "AS",
                })

    splice_types = []
    splice_colors = {}
    if splice_rows:
        splice_types = [t for t, _ in Counter([s["type"] for s in splice_rows]).most_common()]
        splice_colors = build_category_colors(splice_types)

    eggnog_categories = []
    eggnog_colors = {}
    if eggnog_map:
        eggnog_categories = [c for c, _ in eggnog_counts.most_common()]
        eggnog_colors = build_category_colors(eggnog_categories)

    # ---- Plot ----
    image_cfg = {}
    legend_cfg = {}
    expr_cbar_cfg = {}
    if args.layout:
        image_cfg = layout.get("image", {})
        legend_cfg = layout.get("legend", {})
        expr_cbar_cfg = layout.get("expression_colorbar", {})
        # auto layout for expression/splicing: place outside other tracks without compression
        auto_expr_splice = bool(layout.get("auto_expr_splice", True))
        expr_track_width = float(layout.get("expression_track_width", 0.04))
        splicing_width = float(layout.get("splicing_width", 0.06))
        expr_gap = float(layout.get("expr_gap", 0.02))
        base_outer = outer_radius

        auto_bounds = {}
        if auto_expr_splice:
            inner_max = 0.0
            for t in layout.get("tracks", []):
                t_type = t.get("type")
                if t_type in ("expression", "splicing"):
                    continue
                r0 = t.get("r0")
                r1 = t.get("r1")
                if r0 is None or r1 is None:
                    continue
                r1 = float(r1)
                if r1 <= 1.5:
                    r1 = r1 * base_outer
                inner_max = max(inner_max, r1)
            if expr_conditions:
                expr_total = expr_track_width * max(1, len(expr_conditions))
                ex_r0 = inner_max + expr_gap
                ex_r1 = ex_r0 + expr_total
                auto_bounds["expression"] = (ex_r0, ex_r1)
                if splice_rows:
                    sp_r0 = ex_r1
                    sp_r1 = sp_r0 + splicing_width
                    auto_bounds["splicing"] = (sp_r0, sp_r1)
                    outer_radius = max(outer_radius, sp_r1)
                else:
                    outer_radius = max(outer_radius, ex_r1)
            elif splice_rows:
                sp_r0 = inner_max + expr_gap
                sp_r1 = sp_r0 + splicing_width
                auto_bounds["splicing"] = (sp_r0, sp_r1)
                outer_radius = max(outer_radius, sp_r1)
        else:
            base_outer = outer_radius
            auto_bounds = {}

        layout["_auto_bounds"] = auto_bounds
        layout["_base_outer"] = base_outer

        auto_size = bool(image_cfg.get("auto_size", True))
        base_fig = image_cfg.get("base_figsize", image_cfg.get("figsize", [15, 15]))
        expr_baseline = int(image_cfg.get("expr_baseline", 3))
        expr_step = float(image_cfg.get("expr_size_step", 0.6))
        extra = 0.0
        if auto_size and expr_conditions:
            extra = max(0, len(expr_conditions) - expr_baseline) * expr_step
        size_scale = outer_radius / base_outer if base_outer > 0 else 1.0
        figsize = [base_fig[0] * size_scale + extra, base_fig[1] * size_scale + extra]
    else:
        figsize = image_cfg.get("figsize", [15, 15])

    fig = plt.figure(figsize=tuple(figsize))
    ax = plt.subplot(111, polar=True)
    ax.set_theta_direction(-1)
    ax.set_theta_offset(math.pi / 2.0)
    ax.set_axis_off()

    if args.conf and not args.layout:
        args.layout = args.conf

    drawn_types = set()

    if args.layout:
        layout = load_layout(args.layout)
        outer_radius = float(layout.get("outer_radius", 1.7))
        contig_label_radius = float(layout.get("contig_label_radius", outer_radius * 1.08))
        contig_label_size = float(layout.get("contig_label_size", 14))
        contig_label_color = layout.get("contig_label_color", "#8b5a2b")
        contig_label_style = layout.get("contig_label_style", "tangent")
        contig_label_angle_deg = float(layout.get("contig_label_angle_deg", 0.0))
        contig_label_keep_upright = bool(layout.get("contig_label_keep_upright", True))
        contig_label_max = int(layout.get("contig_label_max", args.label_max))
        contig_label_min_len = int(layout.get("contig_label_min_len", args.label_min_len))
        contig_label_enable = bool(layout.get("contig_label_enable", True))
        label_mode = layout.get("label_mode", "text")
        legend_cfg = layout.get("legend", {})
        expr_cbar_cfg = layout.get("expression_colorbar", {})
        auto_expr_splice = bool(layout.get("auto_expr_splice", True))
        number_label_angle_deg = float(layout.get("number_label_angle_deg", 0.0))
        number_label_angle_step = float(layout.get("number_label_angle_step", 6.0))
        number_line_length = float(layout.get("number_line_length", 0.18))
        number_line_width = float(layout.get("number_line_width", 2.0))
        number_circle_size = float(layout.get("number_circle_size", 220))
        number_modes = {"number", "number_circle"}
        auto_bounds = layout.get("_auto_bounds", {})
        base_outer = layout.get("_base_outer", outer_radius)

        if args.label_contigs and contig_label_enable:
            label_contigs(
                ax,
                contigs,
                total_len,
                contig_label_radius,
                contig_label_max,
                contig_label_min_len,
                fontsize=contig_label_size,
                color=contig_label_color,
                style=contig_label_style,
                angle_offset_deg=contig_label_angle_deg,
                keep_upright=contig_label_keep_upright,
            )

        track_number = 0
        track_map = []
        for t in layout.get("tracks", []):
            t_type = t.get("type")
            r0 = t.get("r0")
            r1 = t.get("r1")
            if t_type in auto_bounds:
                r0, r1 = auto_bounds[t_type]
            if r0 is None or r1 is None:
                continue
            r0 = float(r0)
            r1 = float(r1)
            if r1 <= 1.5 and t_type not in auto_bounds:
                r0 = r0 * base_outer
                r1 = r1 * base_outer
            height = r1 - r0
            label = t.get("label", "")
            label_style = t.get("label_style", "radial")
            label_angle = math.radians(float(t.get("label_angle_deg", 0.0)))
            label_offset = float(t.get("label_offset", height * 0.6))
            label_radius = float(t.get("label_radius", r1))
            if label_radius <= 1.5 and t_type not in auto_bounds:
                label_radius = label_radius * base_outer
            label_size = float(t.get("label_size", 10))
            label_color = t.get("label_color", t.get("color", None))

            def default_color_for_type(tp):
                return {
                    "scaffolds": "#1e3a8a",
                    "gc_wave": "#16a34a",
                    "gene_wave": "#2563eb",
                    "repeat_wave": "#f97316",
                    "repeat_tiles": "#0ea5e9",
                    "bgc": "#fb7185",
                    "eggnog": "#06b6d4",
                    "expression": "#ef4444",
                    "splicing": "#f59e0b",
                }.get(tp, "#8b5a2b")

            def put_label(track_desc):
                nonlocal track_number
                track_number += 1
                color = label_color or default_color_for_type(t_type)
                track_map.append((track_number, track_desc, color, r1))
                if label_mode in number_modes:
                    label_text = str(track_number)
                else:
                    label_text = label or track_desc
                if not label_text:
                    return
                if label_mode != "none":
                    if label_style == "tangent":
                        add_track_label(ax, label_radius, label_text, color="#8b5a2b", angle=label_angle, fontsize=label_size)
                    else:
                        add_track_label_radial(ax, label_radius, label_text, color="#8b5a2b", angle=label_angle, offset=label_offset, fontsize=label_size)

            if t_type == "scaffolds":
                fill_color = t.get("fill_color", t.get("color", "#1e3a8a"))
                if isinstance(fill_color, str) and fill_color.lower() in ("none", "transparent"):
                    fill_color = "none"
                edge_color = t.get("edge_color", "white")
                edge_width = float(t.get("edge_width", 0.3))
                for c in contigs:
                    start = c["offset"]
                    end = c["offset"] + c["length"]
                    theta1 = to_angle(start, total_len)
                    theta2 = to_angle(end, total_len)
                    ax.bar(
                        x=(theta1 + theta2) / 2,
                        height=height,
                        width=(theta2 - theta1),
                        bottom=r0,
                        color=fill_color,
                        edgecolor=edge_color,
                        linewidth=edge_width,
                        align="center",
                    )
                put_label("Scaffolds")
                drawn_types.add("scaffolds")
            elif t_type == "gc_wave":
                plot_wave_track(ax, contigs, total_len, gc_by_contig, window, r0, amp=height * 0.5, color=t.get("color", "#16a34a"))
                put_label("GC content")
                drawn_types.add("gc")
            elif t_type == "gene_wave":
                plot_wave_track(ax, contigs, total_len, gene_density_by_contig, window, r0, amp=height * 0.5, color=t.get("color", "#2563eb"))
                put_label("Gene density")
                drawn_types.add("gene")
            elif t_type == "repeat_wave":
                plot_wave_track(ax, contigs, total_len, repeat_density_by_contig, window, r0, amp=height * 0.5, color=t.get("color", "#f97316"))
                put_label("Repeat density")
                drawn_types.add("repeat")
            elif t_type == "repeat_tiles":
                if repeat_rows:
                    repeat_colors = ["#06b6d4", "#38bdf8", "#7dd3fc", "#a5f3fc"]
                    for idx, r in enumerate(repeat_rows):
                        contig = next((c for c in contigs if c["name"] == r["contig"]), None)
                        if not contig:
                            continue
                        gstart = contig["offset"] + (r["start"] - 1)
                        gend = contig["offset"] + r["end"]
                        theta1 = to_angle(gstart, total_len)
                        theta2 = to_angle(gend, total_len)
                        ax.bar(
                            x=(theta1 + theta2) / 2,
                            height=height * 0.6,
                            width=(theta2 - theta1),
                            bottom=r0,
                            color=repeat_colors[idx % len(repeat_colors)],
                            edgecolor=None,
                            linewidth=0,
                        )
                    put_label("Repeat tiles")
                    drawn_types.add("repeat_tiles")
            elif t_type == "bgc":
                if bgc_rows:
                    for b in bgc_rows:
                        contig = next((c for c in contigs if c["name"] == b["contig"]), None)
                        if not contig:
                            continue
                        gstart = contig["offset"] + (b["start"] - 1)
                        gend = contig["offset"] + b["end"]
                        theta1 = to_angle(gstart, total_len)
                        theta2 = to_angle(gend, total_len)
                        add_bgc_arrow(ax, theta1, theta2, r0, r1, b.get("strand", "+"), t.get("color", "#fb7185"))
                    put_label("BGC")
                    drawn_types.add("bgc")
            elif t_type == "eggnog":
                if eggnog_map:
                    empty_color = t.get("empty_color", "#f1f5f9")
                    for c in contigs:
                        start = c["offset"]
                        end = c["offset"] + c["length"]
                        theta1 = to_angle(start, total_len)
                        theta2 = to_angle(end, total_len)
                        ax.bar(
                            x=(theta1 + theta2) / 2,
                            height=height,
                            width=(theta2 - theta1),
                            bottom=r0,
                            color=empty_color,
                            edgecolor=None,
                            linewidth=0,
                        )
                    for g in genes:
                        gid = g["id"] or g["attrs"].get("ID", "")
                        cog = eggnog_map.get(gid, "")
                        if not cog:
                            continue
                        contig = next((c for c in contigs if c["name"] == g["seqid"]), None)
                        if not contig:
                            continue
                        gstart = contig["offset"] + (g["start"] - 1)
                        gend = contig["offset"] + g["end"]
                        theta1 = to_angle(gstart, total_len)
                        theta2 = to_angle(gend, total_len)
                        ax.bar(
                            x=(theta1 + theta2) / 2,
                            height=height * 0.9,
                            width=(theta2 - theta1),
                            bottom=r0,
                            color=eggnog_colors.get(cog, "#06b6d4"),
                            edgecolor=None,
                            linewidth=0,
                        )
                    put_label("EggNOG")
                    drawn_types.add("eggnog")
            elif t_type == "expression":
                if expr_conditions:
                    cmap = plt.cm.RdYlGn_r
                    vmax = expr_max_abs if expr_max_abs > 0 else 1.0
                    cond = t.get("condition")
                    if cond and cond != "all":
                        conds = [cond] if cond in expr_conditions else []
                    else:
                        conds = expr_conditions
                    if conds:
                        sub_height = height / max(1, len(conds))
                        for idx, cnd in enumerate(conds):
                            base = r0 + idx * sub_height
                            plot_heatmap_ring(
                                ax, contigs, total_len, window, expr_window_values.get(cnd, {}),
                                base=base, height=sub_height * 0.9, cmap=cmap, vmin=-vmax, vmax=vmax, empty_color="#f1f5f9"
                            )
                            if label_mode == "text" and label:
                                add_track_label_radial(ax, base + sub_height * 0.6, f"Expr {cnd}", color="#8b5a2b", angle=label_angle, offset=label_offset)
                    put_label("Expression")
                    drawn_types.add("expression")
            elif t_type == "splicing":
                if splice_rows:
                    empty_color = t.get("empty_color", "#f1f5f9")
                    for c in contigs:
                        start = c["offset"]
                        end = c["offset"] + c["length"]
                        theta1 = to_angle(start, total_len)
                        theta2 = to_angle(end, total_len)
                        ax.bar(
                            x=(theta1 + theta2) / 2,
                            height=height,
                            width=(theta2 - theta1),
                            bottom=r0,
                            color=empty_color,
                            edgecolor=None,
                            linewidth=0,
                        )
                    for s in splice_rows:
                        contig = next((c for c in contigs if c["name"] == s["contig"]), None)
                        if not contig:
                            continue
                        gstart = contig["offset"] + (s["start"] - 1)
                        gend = contig["offset"] + s["end"]
                        theta1 = to_angle(gstart, total_len)
                        theta2 = to_angle(gend, total_len)
                        ax.bar(
                            x=(theta1 + theta2) / 2,
                            height=height * 0.9,
                            width=(theta2 - theta1),
                            bottom=r0,
                            color=splice_colors.get(s["type"], "#e5e7eb"),
                            edgecolor=None,
                            linewidth=0,
                        )
                    put_label("Alternative splicing")
                    drawn_types.add("splicing")
        # track legend for number mode
        if label_mode in number_modes and track_map:
            for idx, (num, name, color, r1) in enumerate(track_map):
                angle = math.radians(number_label_angle_deg + idx * number_label_angle_step)
                r_start = r1
                r_end = r1 + number_line_length
                ax.plot([angle, angle], [r_start, r_end], color=color, lw=number_line_width)
                ax.scatter([angle], [r_end], s=number_circle_size, facecolor="white", edgecolor=color, linewidth=number_line_width)
                ax.text(
                    angle,
                    r_end,
                    str(num),
                    fontsize=10,
                    fontweight="bold",
                    ha="center",
                    va="center",
                    color=color,
                )
            track_handles = [Line2D([0], [0], color="none", label=f"{num}. {name}") for num, name, _, _ in track_map]
            legend_tracks = fig.legend(handles=track_handles, loc="upper left", bbox_to_anchor=(-0.18, 1.02), frameon=False, fontsize=8)
            for text in legend_tracks.get_texts():
                text.set_color("#8b5a2b")
                text.set_fontweight("bold")
        # end layout mode
    else:
        track_total = 1 + 3  # contigs + GC/gene/repeat waves
        if repeat_rows:
            track_total += 1
        if bgc_rows:
            track_total += 1
        if eggnog_map:
            track_total += 1
        if expr_conditions:
            track_total += len(expr_conditions)
        if splice_rows:
            track_total += 1

        outer_radius = 1.7
        inner_min = 0.2
        gap_ratio = 0.4
        ring_step = (outer_radius - inner_min) / (track_total + (track_total - 1) * gap_ratio)
        ring_gap = ring_step * gap_ratio
        ring = outer_radius

    if not args.layout:
        # Splicing track (categorical heatmap) - outermost
        if splice_rows:
            empty_color = "#f1f5f9"
            for c in contigs:
                start = c["offset"]
                end = c["offset"] + c["length"]
                theta1 = to_angle(start, total_len)
                theta2 = to_angle(end, total_len)
                ax.bar(
                    x=(theta1 + theta2) / 2,
                    height=ring_step * 0.8,
                    width=(theta2 - theta1),
                    bottom=ring,
                    color=empty_color,
                    edgecolor=None,
                    linewidth=0,
                )
            for s in splice_rows:
                contig = next((c for c in contigs if c["name"] == s["contig"]), None)
                if not contig:
                    continue
                gstart = contig["offset"] + (s["start"] - 1)
                gend = contig["offset"] + s["end"]
                theta1 = to_angle(gstart, total_len)
                theta2 = to_angle(gend, total_len)
                ax.bar(
                    x=(theta1 + theta2) / 2,
                    height=ring_step * 0.8,
                    width=(theta2 - theta1),
                    bottom=ring,
                    color=splice_colors.get(s["type"], "#e5e7eb"),
                    edgecolor=None,
                    linewidth=0,
                )
            add_track_label_radial(ax, ring + ring_step * 0.35, "Alternative splicing", color="#8b5a2b", offset=ring_step * 0.7)
            ring -= (ring_step + ring_gap)
    
        # Expression heatmap tracks (one ring per condition)
        if expr_conditions:
            cmap = plt.cm.RdYlGn_r  # low=green, high=red
            vmax = expr_max_abs if expr_max_abs > 0 else 1.0
            for cond in expr_conditions:
                plot_heatmap_ring(
                    ax,
                    contigs,
                    total_len,
                    window,
                    expr_window_values.get(cond, {}),
                    base=ring,
                    height=ring_step * 0.8,
                    cmap=cmap,
                    vmin=-vmax,
                    vmax=vmax,
                    empty_color="#f1f5f9",
                )
                add_track_label_radial(ax, ring + ring_step * 0.35, f"Expr {cond}", color="#8b5a2b", offset=ring_step * 0.7)
                ring -= (ring_step + ring_gap)
    
        # contig blocks
        for c in contigs:
            start = c["offset"]
            end = c["offset"] + c["length"]
            theta1 = to_angle(start, total_len)
            theta2 = to_angle(end, total_len)
            ax.bar(
                x=(theta1 + theta2) / 2,
                height=ring_step * 0.6,
                width=(theta2 - theta1),
                bottom=ring,
                color="#1e3a8a",
                edgecolor="white",
                linewidth=0.3,
                align="center",
            )
    
        if args.label_contigs:
            label_contigs(ax, contigs, total_len, outer_radius + ring_step * 6.0, args.label_max, args.label_min_len, fontsize=14)
        add_track_label_radial(ax, ring + ring_step * 0.35, "Scaffolds", color="#8b5a2b", offset=ring_step * 0.7, fontsize=10)
    
        ring -= (ring_step + ring_gap)
    
        # GC content wave
        plot_wave_track(
            ax,
            contigs,
            total_len,
            gc_by_contig,
            window,
            base=ring,
            amp=ring_step * 0.45,
            color="#16a34a",
            linewidth=1.1,
        )
        add_track_label_radial(ax, ring, "GC content", color="#8b5a2b", offset=ring_step * 0.7, fontsize=10)
        ring -= (ring_step + ring_gap)
    
        # gene density wave
        plot_wave_track(
            ax,
            contigs,
            total_len,
            gene_density_by_contig,
            window,
            base=ring,
            amp=ring_step * 0.45,
            color="#2563eb",
            linewidth=1.1,
        )
        add_track_label_radial(ax, ring, "Gene density", color="#8b5a2b", offset=ring_step * 0.7, fontsize=10)
        ring -= (ring_step + ring_gap)
    
        # repeat density wave
        plot_wave_track(
            ax,
            contigs,
            total_len,
            repeat_density_by_contig,
            window,
            base=ring,
            amp=ring_step * 0.45,
            color="#f97316",
            linewidth=1.1,
        )
        add_track_label_radial(ax, ring, "Repeat density", color="#8b5a2b", offset=ring_step * 0.7, fontsize=10)
        ring -= (ring_step + ring_gap)
    
        # repeat tracks
        repeat_colors = ["#f97316", "#8b5cf6", "#06b6d4", "#84cc16", "#ec4899"]
        for idx, t in enumerate(top_repeat_types):
            for r in feat_by_type[t]:
                contig = next(c for c in contigs if c["name"] == r["seqid"])
                gstart = contig["offset"] + (r["start"] - 1)
                gend = contig["offset"] + r["end"]
                theta1 = to_angle(gstart, total_len)
                theta2 = to_angle(gend, total_len)
                ax.bar(
                    x=(theta1 + theta2) / 2,
                    height=ring_step * 0.5,
                    width=(theta2 - theta1),
                    bottom=ring,
                    color=repeat_colors[idx % len(repeat_colors)],
                    edgecolor=None,
                    linewidth=0,
                )
            ring -= (ring_step + ring_gap)
    
        # BGC track
        if bgc_rows:
            for b in bgc_rows:
                contig = next((c for c in contigs if c["name"] == b["contig"]), None)
                if not contig:
                    continue
                gstart = contig["offset"] + (b["start"] - 1)
                gend = contig["offset"] + b["end"]
                theta1 = to_angle(gstart, total_len)
                theta2 = to_angle(gend, total_len)
                add_bgc_arrow(
                    ax,
                    theta1,
                    theta2,
                    r0=ring,
                    r1=ring + ring_step * 0.6,
                    strand=b.get("strand", "+"),
                    color="#fb7185",
                )
            add_track_label_radial(ax, ring, "BGC (+/-)", color="#8b5a2b", offset=ring_step * 0.7)
            ring -= (ring_step + ring_gap)
    
        # EggNOG track (categorical heatmap aligned to genes)
        if eggnog_map:
            for g in genes:
                gid = g["id"] or g["attrs"].get("ID", "")
                cog = eggnog_map.get(gid, "")
                if not cog:
                    continue
                contig = next((c for c in contigs if c["name"] == g["seqid"]), None)
                if not contig:
                    continue
                gstart = contig["offset"] + (g["start"] - 1)
                gend = contig["offset"] + g["end"]
                theta1 = to_angle(gstart, total_len)
                theta2 = to_angle(gend, total_len)
                ax.bar(
                    x=(theta1 + theta2) / 2,
                    height=ring_step * 0.8,
                    width=(theta2 - theta1),
                    bottom=ring,
                    color=eggnog_colors.get(cog, "#06b6d4"),
                    edgecolor=None,
                    linewidth=0,
                )
            add_track_label_radial(ax, ring, "EggNOG", color="#8b5a2b", offset=ring_step * 0.7)
            ring -= (ring_step + ring_gap)
    
        # (splicing moved to outermost)
    
    if not args.layout:
        if True:
            drawn_types.update({"gc", "gene", "repeat", "scaffolds"})
        if repeat_rows:
            drawn_types.add("repeat_tiles")
        if bgc_rows:
            drawn_types.add("bgc")
        if eggnog_map:
            drawn_types.add("eggnog")
        if expr_conditions:
            drawn_types.add("expression")
        if splice_rows:
            drawn_types.add("splicing")

    labels_cfg = legend_cfg.get("labels", {}) if args.layout else {}
    label_gc = labels_cfg.get("gc", "GC content")
    label_gene = labels_cfg.get("gene", "Gene density")
    label_repeat = labels_cfg.get("repeat", "Repeat density")
    label_bgc = labels_cfg.get("bgc", "BGC (+/-)")
    label_expr = labels_cfg.get("expression", "")

    legend_handles = []
    if "gc" in drawn_types:
        legend_handles.append(Line2D([0], [0], color="#16a34a", lw=2, label=label_gc))
    if "gene" in drawn_types:
        legend_handles.append(Line2D([0], [0], color="#2563eb", lw=2, label=label_gene))
    if "repeat" in drawn_types:
        legend_handles.append(Line2D([0], [0], color="#f97316", lw=2, label=label_repeat))
    if "bgc" in drawn_types:
        legend_handles.append(Patch(facecolor="#fb7185", edgecolor="none", label=label_bgc))
    if "expression" in drawn_types and label_expr:
        legend_handles.append(Patch(facecolor="#6b7280", edgecolor="none", label=label_expr))

    if args.layout:
        main_cfg = legend_cfg.get("main", {})
        legend_main = fig.legend(
            handles=legend_handles,
            loc=main_cfg.get("loc", "upper right"),
            bbox_to_anchor=tuple(main_cfg.get("bbox", [1.28, 0.95])),
            frameon=False,
            fontsize=main_cfg.get("fontsize", 11),
            handlelength=main_cfg.get("handlelength", 2.6),
            handleheight=main_cfg.get("handleheight", 1.9),
            borderpad=main_cfg.get("borderpad", 0.8),
            labelspacing=main_cfg.get("labelspacing", 0.7),
            markerscale=main_cfg.get("markerscale", 1.6),
        )
    else:
        legend_main = fig.legend(
            handles=legend_handles,
            loc="upper right",
            bbox_to_anchor=(1.28, 0.95),
            frameon=False,
            fontsize=11,
            handlelength=2.6,
            handleheight=1.9,
            borderpad=0.8,
            labelspacing=0.7,
            markerscale=1.6,
        )
    legend_text_color = legend_cfg.get("text_color", "#8b5a2b") if args.layout else "#8b5a2b"
    legend_text_weight = legend_cfg.get("fontweight", "bold") if args.layout else "bold"
    for text in legend_main.get_texts():
        text.set_color(legend_text_color)
        text.set_fontweight(legend_text_weight)

    if eggnog_categories and "eggnog" in drawn_types:
        egg_prefix = legend_cfg.get("eggnog_prefix", "COG:") if args.layout else "COG:"
        egg_handles = [Patch(facecolor=eggnog_colors[c], edgecolor="none", label=f"{egg_prefix}{c}") for c in eggnog_categories]
        egg_cols = 2 if len(egg_handles) > 8 else 1
        if args.layout:
            egg_cfg = legend_cfg.get("eggnog", {})
            legend_egg = fig.legend(
                handles=egg_handles,
                loc=egg_cfg.get("loc", "upper right"),
                bbox_to_anchor=tuple(egg_cfg.get("bbox", [1.28, 0.70])),
                frameon=False,
                fontsize=egg_cfg.get("fontsize", 10),
                ncol=egg_cfg.get("ncol", egg_cols),
                handlelength=egg_cfg.get("handlelength", 2.6),
                handleheight=egg_cfg.get("handleheight", 1.9),
                borderpad=egg_cfg.get("borderpad", 0.8),
                labelspacing=egg_cfg.get("labelspacing", 0.7),
            )
        else:
            legend_egg = fig.legend(
                handles=egg_handles,
                loc="upper right",
                bbox_to_anchor=(1.28, 0.70),
                frameon=False,
                fontsize=10,
                ncol=egg_cols,
                handlelength=2.6,
                handleheight=1.9,
                borderpad=0.8,
                labelspacing=0.7,
            )
        for text in legend_egg.get_texts():
            text.set_color(legend_text_color)
            text.set_fontweight(legend_text_weight)

    if splice_types and "splicing" in drawn_types:
        sp_prefix = legend_cfg.get("splicing_prefix", "Splicing ") if args.layout else "Splicing "
        splice_handles = [Patch(facecolor=splice_colors[t], edgecolor="none", label=f"{sp_prefix}{t}") for t in splice_types]
        splice_cols = 2 if len(splice_handles) > 8 else 1
        if args.layout:
            sp_cfg = legend_cfg.get("splicing", {})
            legend_splice = fig.legend(
                handles=splice_handles,
                loc=sp_cfg.get("loc", "upper right"),
                bbox_to_anchor=tuple(sp_cfg.get("bbox", [1.28, 0.38])),
                frameon=False,
                fontsize=sp_cfg.get("fontsize", 10),
                ncol=sp_cfg.get("ncol", splice_cols),
                handlelength=sp_cfg.get("handlelength", 2.6),
                handleheight=sp_cfg.get("handleheight", 1.9),
                borderpad=sp_cfg.get("borderpad", 0.8),
                labelspacing=sp_cfg.get("labelspacing", 0.7),
            )
        else:
            legend_splice = fig.legend(
                handles=splice_handles,
                loc="upper right",
                bbox_to_anchor=(1.28, 0.38),
                frameon=False,
                fontsize=10,
                ncol=splice_cols,
                handlelength=2.6,
                handleheight=1.9,
                borderpad=0.8,
                labelspacing=0.7,
            )
        for text in legend_splice.get_texts():
            text.set_color(legend_text_color)
            text.set_fontweight(legend_text_weight)

    if expr_conditions and "expression" in drawn_types and expr_cbar_cfg.get("show", True):
        from matplotlib.colors import Normalize
        sm = plt.cm.ScalarMappable(cmap=plt.cm.RdYlGn_r, norm=Normalize(vmin=-expr_max_abs if expr_max_abs else -1, vmax=expr_max_abs if expr_max_abs else 1))
        sm.set_array([])
        cpos = expr_cbar_cfg.get("pos", [0.84, 0.78, 0.14, 0.035]) if args.layout else [0.84, 0.78, 0.14, 0.035]
        cax = fig.add_axes(cpos)
        cbar = fig.colorbar(sm, cax=cax, orientation="horizontal")
        cbar_label = expr_cbar_cfg.get("label", label_expr) if args.layout else "Expression (heatmap)"
        cbar_font = expr_cbar_cfg.get("fontsize", 9) if args.layout else 9
        cbar_tick = expr_cbar_cfg.get("ticksize", 8) if args.layout else 8
        cbar.set_label(cbar_label, fontsize=cbar_font)
        cbar.ax.xaxis.label.set_color(legend_text_color)
        cbar.ax.xaxis.label.set_fontweight(legend_text_weight)
        cbar.ax.tick_params(colors=legend_text_color, labelsize=cbar_tick)
    # Splicing is categorical; legend handles colors

    dpi = image_cfg.get("dpi", 300) if args.layout else 300
    bbox = image_cfg.get("bbox_inches", "tight") if args.layout else "tight"
    plt.savefig(args.output, dpi=dpi, bbox_inches=bbox)

    # ---- Excel output ----
    if pd is None:
        print("pandas not available; skipping Excel output")
        return

    with pd.ExcelWriter(args.excel) as writer:
        pd.DataFrame(gc_rows).to_excel(writer, index=False, sheet_name="gc_content")
        pd.DataFrame(gene_density_rows).to_excel(writer, index=False, sheet_name="gene_density")
        pd.DataFrame(repeat_density_rows).to_excel(writer, index=False, sheet_name="repeat_density")
        pd.DataFrame(repeat_rows).to_excel(writer, index=False, sheet_name="repeats")
        pd.DataFrame(bgc_rows).to_excel(writer, index=False, sheet_name="bgc")
        pd.DataFrame(splice_rows).to_excel(writer, index=False, sheet_name="splicing")
        if eggnog_map:
            egg_rows = []
            for g in genes:
                gid = g["id"] or g["attrs"].get("ID", "")
                cog = eggnog_map.get(gid, "")
                if cog:
                    egg_rows.append({
                        "gene_id": gid,
                        "contig": g["seqid"],
                        "start": g["start"],
                        "end": g["end"],
                        "cog_category": cog,
                    })
            pd.DataFrame(egg_rows).to_excel(writer, index=False, sheet_name="eggnog")
        if expr_map:
            expr_rows = []
            for gid, items in expr_map.items():
                for e in items:
                    expr_rows.append({
                        "gene_id": gid,
                        "condition": e["condition"],
                        "log2fc": e["log2fc"],
                        "padj": e["padj"],
                        "status": e["status"],
                    })
            pd.DataFrame(expr_rows).to_excel(writer, index=False, sheet_name="expression")

    print(f"Saved {args.output} and {args.excel}")


if __name__ == "__main__":
    main()
