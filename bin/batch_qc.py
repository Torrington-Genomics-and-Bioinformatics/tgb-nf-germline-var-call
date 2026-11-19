#!/usr/bin/env python3
"""Batch QC table and report (pipeline v3.1-qc).

  batch_qc.py <qc_dir> <out.tsv> <out.html> <title> <pipeline_version>

Reads, per sample, the files the pipeline writes to <outdir>/qc and <outdir>/bqsr:
  <s>.fastp.json              fastp (raw reads, Q30, GC)
  <s>.metrics.txt             MarkDuplicates (duplication, estimated library size)
  <s>.hs_metrics.txt          CollectHsMetrics on the run's target intervals (coverage, on-target, fold-80, dropout)
  <s>.per_target_coverage.txt.gz  CollectHsMetrics PER_TARGET_COVERAGE (used for the sex check)
  <s>.samtools.stats          samtools stats (mapped, properly paired, error rate, insert size)
  <s>.selfSM                  VerifyBamID2 (FREEMIX contamination)
  <s>.crosscheck_metrics      CrosscheckFingerprints by read group (lane agreement)
  <s>.recal.table, <s>.after.recal.table   BaseRecalibrator before / after ApplyBQSR
Writes one row per sample with pass / review / FAIL and the reasons, and a self-contained HTML report.
Python standard library only.
"""
import csv, glob, gzip, html, json, math, os, sys, datetime

qd, out_tsv, out_html, title, version = sys.argv[1:6]

# Thresholds (cohort QC, Nov 2026). Values beyond FAIL make a sample unusable without review; REVIEW = look at it.
T = dict(contam_fail=0.05, contam_review=0.02, cov_fail=20, cov_review=30, t20_fail=0.70, t20_review=0.80,
         dup_review=0.30, paired_review=0.97, fold80_review=2.0)


def p(s, ext):
    return os.path.join(qd, f"{s}.{ext}")


def picard(path, cls):
    if not os.path.exists(path):
        return None
    lines = open(path).read().split("\n")
    for i, l in enumerate(lines):
        if l.startswith("## METRICS CLASS") and cls in l:
            return dict(zip(lines[i + 1].split("\t"), lines[i + 2].split("\t")))
    return None


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def recal1(path):
    """RecalTable1 mismatch rows: (reported Q, empirical Q, observations)."""
    rows, on = [], False
    for l in open(path):
        if l.startswith("#:GATKTable:RecalTable1"):
            on = True
            continue
        if on:
            if not l.strip():
                break
            f = l.split()
            if f[0] == "ReadGroup" or f[2] != "M":
                continue
            rows.append((int(f[1]), float(f[3]), int(f[4])))
    return rows


def rmse(rows):
    n = sum(o for _, _, o in rows)
    return math.sqrt(sum(o * (e - q) ** 2 for q, e, o in rows) / n) if n else None


# PAR regions (GRCh38) are excluded from the sex-chromosome coverage
PAR = {"chrX": [(10001, 2781479), (155701383, 156030895)], "chrY": [(10001, 2781479), (56887903, 57217415)]}
AUTO = {f"chr{i}" for i in range(1, 23)}


def sex_check(path):
    """Length-weighted mean target coverage of chrX and chrY (non-PAR), relative to the autosomes."""
    if not os.path.exists(path):
        return None
    cov = {"auto": [0.0, 0], "chrX": [0.0, 0], "chrY": [0.0, 0]}
    with gzip.open(path, "rt") as f:
        h = f.readline().rstrip("\n").split("\t")
        ci, si, ei, mi = h.index("chrom"), h.index("start"), h.index("end"), h.index("mean_coverage")
        for l in f:
            v = l.rstrip("\n").split("\t")
            c, a, b, m = v[ci], int(v[si]), int(v[ei]), float(v[mi])
            if c in AUTO:
                k = "auto"
            elif c in PAR:
                if any(a <= pe and b >= ps for ps, pe in PAR[c]):
                    continue
                k = c
            else:
                continue
            n = b - a + 1
            cov[k][0] += m * n
            cov[k][1] += n
    au = cov["auto"][0] / cov["auto"][1] if cov["auto"][1] else 0
    if not au:
        return None
    x = cov["chrX"][0] / cov["chrX"][1] / au if cov["chrX"][1] else float("nan")
    y = cov["chrY"][0] / cov["chrY"][1] / au if cov["chrY"][1] else float("nan")
    # chrY decides (female exomes ~0-0.05, male ~0.3-0.5); chrX in male exomes runs ~0.5-0.75 (capture bias)
    if x > 0.8 and y < 0.1:
        call = "XX"
    elif x < 0.85 and y > 0.2:
        call = "XY"
    else:
        call = "check"
    return x, y, call


samples = sorted({os.path.basename(f)[:-len(".fastp.json")] for f in glob.glob(os.path.join(qd, "*.fastp.json"))}
                 | {os.path.basename(f)[:-len(".hs_metrics.txt")] for f in glob.glob(os.path.join(qd, "*.hs_metrics.txt"))})
cols = ["sample", "raw_reads_M", "q30_pct", "gc_pct", "dup_pct", "est_library_size_M", "mapped_pct",
        "properly_paired_pct", "error_rate_pct", "insert_mean", "mean_target_cov", "median_target_cov",
        "pct_targets_20x", "pct_targets_30x", "pct_on_target", "fold80", "at_dropout", "gc_dropout", "freemix",
        "read_groups", "lane_mismatch", "bqsr_rmse_before", "bqsr_rmse_after", "chrX_ratio", "chrY_ratio",
        "sex_genetic", "status", "flags"]
rows = []
for s in samples:
    r = dict(sample=s)
    fail, rev, miss = [], [], []
    if os.path.exists(p(s, "fastp.json")):
        b = json.load(open(p(s, "fastp.json")))["summary"]["before_filtering"]
        r.update(raw_reads_M=f"{b['total_reads'] / 1e6:.1f}", q30_pct=f"{100 * b['q30_rate']:.1f}",
                 gc_pct=f"{100 * b['gc_content']:.1f}")
    else:
        miss.append("fastp")
    md = picard(p(s, "metrics.txt"), "DuplicationMetrics")
    if md:
        d = num(md.get("PERCENT_DUPLICATION")) or 0
        lib = num(md.get("ESTIMATED_LIBRARY_SIZE"))
        r.update(dup_pct=f"{100 * d:.1f}", est_library_size_M=f"{lib / 1e6:.0f}" if lib else "")
        if d > T["dup_review"]:
            rev.append("dup>30%")
    else:
        miss.append("markdup")
    if os.path.exists(p(s, "samtools.stats")):
        v = {}
        for l in open(p(s, "samtools.stats")):
            if l.startswith("SN\t"):
                f = l.rstrip("\n").split("\t")
                v[f[1].rstrip(":")] = float(f[2])
        n = v.get("raw total sequences") or 0
        if n:
            pp = v["reads properly paired"] / n
            r.update(mapped_pct=f"{100 * v['reads mapped'] / n:.2f}", properly_paired_pct=f"{100 * pp:.2f}",
                     error_rate_pct=f"{100 * v['error rate']:.3f}", insert_mean=f"{v['insert size average']:.0f}")
            if pp < T["paired_review"]:
                rev.append("paired<97%")
    else:
        miss.append("samtools_stats")
    hs = picard(p(s, "hs_metrics.txt"), "HsMetrics")
    if hs:
        mc, t20, f80 = num(hs["MEAN_TARGET_COVERAGE"]), num(hs["PCT_TARGET_BASES_20X"]), num(hs["FOLD_80_BASE_PENALTY"])
        r.update(mean_target_cov=f"{mc:.1f}", median_target_cov=hs["MEDIAN_TARGET_COVERAGE"],
                 pct_targets_20x=f"{100 * t20:.1f}", pct_targets_30x=f"{100 * num(hs['PCT_TARGET_BASES_30X']):.1f}",
                 pct_on_target=f"{100 * num(hs['PCT_SELECTED_BASES']):.1f}",
                 fold80=f"{f80:.2f}" if f80 is not None else "", at_dropout=f"{num(hs['AT_DROPOUT']):.2f}",
                 gc_dropout=f"{num(hs['GC_DROPOUT']):.2f}")
        if mc < T["cov_fail"]:
            fail.append("cov<20x")
        elif mc < T["cov_review"]:
            rev.append("cov<30x")
        if t20 < T["t20_fail"]:
            fail.append("20x<70%")
        elif t20 < T["t20_review"]:
            rev.append("20x<80%")
        if f80 is not None and f80 > T["fold80_review"]:
            rev.append("fold80>2")
    else:
        miss.append("hsmetrics")
    if os.path.exists(p(s, "selfSM")):
        v = open(p(s, "selfSM")).read().split("\n")[1].split("\t")
        fm = num(v[6]) or 0
        r["freemix"] = f"{fm:.4f}"
        if fm > T["contam_fail"]:
            fail.append("contam>5%")
        elif fm > T["contam_review"]:
            rev.append("contam>2%")
    else:
        miss.append("verifybamid")
    if os.path.exists(p(s, "crosscheck_metrics")):
        body = [x.rstrip("\n").split("\t") for x in open(p(s, "crosscheck_metrics")) if x.strip() and not x.startswith("#")]
        h, body = body[0], body[1:]
        ri, li = h.index("RESULT"), h.index("LEFT_GROUP_VALUE")
        bad = [b for b in body if b[ri] == "UNEXPECTED_MISMATCH"]
        inc = sum(b[ri] == "INCONCLUSIVE" for b in body)
        r["read_groups"] = str(len({b[li] for b in body}))
        r["lane_mismatch"] = "yes" if bad else ("no" if not inc else f"no ({inc} inconclusive)")
        if bad:
            fail.append("lane_mismatch")
    else:
        miss.append("crosscheck")
    if os.path.exists(p(s, "recal.table")) and os.path.exists(p(s, "after.recal.table")):
        a, b = rmse(recal1(p(s, "recal.table"))), rmse(recal1(p(s, "after.recal.table")))
        r.update(bqsr_rmse_before=f"{a:.2f}" if a is not None else "", bqsr_rmse_after=f"{b:.2f}" if b is not None else "")
    else:
        miss.append("bqsr_tables")
    sx = sex_check(p(s, "per_target_coverage.txt.gz"))
    if sx:
        r.update(chrX_ratio=f"{sx[0]:.2f}", chrY_ratio=f"{sx[1]:.3f}", sex_genetic=sx[2])
        if sx[2] == "check":
            rev.append("sex_unclear")
    else:
        miss.append("per_target_coverage")
    if miss:
        rev.append("missing:" + "+".join(miss))
    r["flags"] = ",".join(fail + rev) or "-"
    r["status"] = "FAIL" if fail else ("review" if rev else "pass")
    rows.append(r)

with open(out_tsv, "w") as o:
    o.write("\t".join(cols) + "\n")
    for r in rows:
        o.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")

# ---- HTML report -------------------------------------------------------------------------------------------
n = {k: sum(r["status"] == k for r in rows) for k in ("pass", "review", "FAIL")}


def med(c):
    x = sorted(float(r[c]) for r in rows if r.get(c) not in (None, ""))
    if not x:
        return "-"
    m = len(x) // 2
    return f"{(x[m] if len(x) % 2 else (x[m - 1] + x[m]) / 2):g}"


show = [("sample", "Sample"), ("status", "Status"), ("raw_reads_M", "Reads (M)"), ("mean_target_cov", "Mean target cov"),
        ("pct_targets_20x", "% targets ≥20x"), ("pct_on_target", "% on target"), ("fold80", "Fold-80"),
        ("dup_pct", "Dup %"), ("freemix", "Contamination (FREEMIX)"), ("properly_paired_pct", "Properly paired %"),
        ("error_rate_pct", "Error rate %"), ("insert_mean", "Insert (bp)"), ("q30_pct", "Q30 %"),
        ("bqsr_rmse_before", "BQSR RMSE before"), ("bqsr_rmse_after", "after"), ("sex_genetic", "Genetic sex"),
        ("read_groups", "Lanes"), ("lane_mismatch", "Lane mismatch"), ("flags", "Flags")]
summ = [("mean_target_cov", "Mean target coverage"), ("pct_targets_20x", "% target bases ≥20x"),
        ("pct_on_target", "% on target"), ("fold80", "Fold-80 penalty"), ("dup_pct", "Duplication %"),
        ("freemix", "Contamination (FREEMIX)"), ("properly_paired_pct", "Properly paired %"), ("q30_pct", "Q30 %")]
e = html.escape
rowhtml = "\n".join(
    "<tr class='%s'>" % r["status"].lower() + "".join(f"<td>{e(str(r.get(c, '')))}</td>" for c, _ in show) + "</tr>"
    for r in rows)
page = f"""<!doctype html><html><head><meta charset="utf-8"><title>{e(title)} QC</title>
<style>
body{{font-family:Arial,Helvetica,sans-serif;margin:24px;color:#222}} h1{{font-size:20px}} h2{{font-size:16px;margin-top:28px}}
table{{border-collapse:collapse;font-size:12px}} th,td{{border:1px solid #ccc;padding:4px 6px;text-align:right}}
th{{background:#f2f2f2}} td:first-child,td:last-child{{text-align:left}}
tr.pass td:nth-child(2){{background:#d9f2d9}} tr.review td:nth-child(2){{background:#fff2cc}}
tr.fail td:nth-child(2){{background:#f8d0d0}} .note{{font-size:12px;color:#555;max-width:900px}}
.cards span{{display:inline-block;margin-right:18px;font-size:14px}}
</style></head><body>
<h1>{e(title)}: sample QC</h1>
<p class="note">Internal report. Pipeline {e(version)}; generated {datetime.datetime.now(datetime.timezone.utc):%Y-%m-%d %H:%M} UTC.
{len(rows)} samples.</p>
<p class="cards"><span>pass <b>{n['pass']}</b></span><span>review <b>{n['review']}</b></span><span>FAIL <b>{n['FAIL']}</b></span></p>
<h2>Batch medians</h2><table><tr>{''.join(f'<th>{e(l)}</th>' for _, l in summ)}</tr>
<tr>{''.join(f'<td>{med(c)}</td>' for c, _ in summ)}</tr></table>
<h2>Per sample</h2><table><tr>{''.join(f'<th>{e(l)}</th>' for _, l in show)}</tr>
{rowhtml}</table>
<h2>Thresholds</h2><p class="note">
FAIL: contamination FREEMIX &gt; {T['contam_fail']}; mean target coverage &lt; {T['cov_fail']}x; target bases ≥20x &lt; {int(100*T['t20_fail'])}%;
any read-group pair of the sample reported UNEXPECTED_MISMATCH by CrosscheckFingerprints.<br>
Review: FREEMIX &gt; {T['contam_review']}; mean target coverage &lt; {T['cov_review']}x; target bases ≥20x &lt; {int(100*T['t20_review'])}%;
duplication &gt; {int(100*T['dup_review'])}%; properly paired &lt; {int(100*T['paired_review'])}%; fold-80 &gt; {T['fold80_review']};
genetic sex unclear; any QC file missing.</p>
<h2>Definitions</h2><p class="note">
Coverage, on-target, fold-80 and AT/GC dropout: GATK CollectHsMetrics on the capture targets used for calling (bait = target).
Fold-80 penalty: how much more sequencing would be needed for 80% of target bases to reach the current mean coverage (1 = perfectly even).
Contamination: VerifyBamID2 FREEMIX (fraction of reads from another person). Mapped, properly paired, error rate (mismatches per
mapped base) and insert size: samtools stats on the final BAM. BQSR RMSE: observation-weighted root-mean-square difference between
reported and empirical base quality, before and after recalibration. Genetic sex: mean target coverage of chrX and chrY (PAR excluded)
relative to the autosomes (XX: X &gt; 0.8 and Y &lt; 0.1; XY: X &lt; 0.85 and Y &gt; 0.2; anything else: check). Compare it with the recorded sex.
Lanes: read groups (flowcell-lanes) in the BAM; lane mismatch: CrosscheckFingerprints between lanes of the same sample.
Detailed plots: multiqc_report.html in the same folder.</p>
</body></html>"""
open(out_html, "w").write(page)
print(f"{len(rows)} samples: {n}")
