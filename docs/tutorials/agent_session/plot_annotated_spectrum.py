"""Publication figure: annotated MS2 spectrum of FFLVLR (scan 21989, PXD044934).

Peaks are matched to theoretical b, y, a2 and phenylalanine-immonium ions at 10 ppm; the
first 13C isotope of a matched ion takes its colour but no label. A fragment ladder over the
sequence shows which b/y ions were observed. The y3 ion expected for the competing
Proteome Discoverer peptide FFVLLR (401.287) is marked where it is absent.

Usage: python plot_annotated_spectrum.py OBSERVED_JSON OUTPUT_SVG
"""

import json
import sys

import matplotlib

matplotlib.use("svg")
import matplotlib.pyplot as plt  # noqa: E402

PEPTIDE = "FFLVLR"
ALT_Y3 = ("FFVLLR", 401.28708)  # y3 (LLR+) of the Proteome Discoverer peptide
TOL_PPM = 10.0
ISOTOPE_SPACING = 1.003355

MONO = {"F": 147.06841, "L": 113.08406, "V": 99.06841, "R": 156.10111}
PROTON, WATER, CO = 1.007276, 18.010565, 27.994915

# Chart chrome and categorical slots from the dataviz reference palette (validated, light mode).
INK, INK_2, MUTED, GRID = "#0b0b0b", "#52514e", "#b5b4ae", "#e4e3df"
COLOURS = {"b": "#2a78d6", "y": "#eb6834", "other": "#1baf7a"}


def theoretical_ions(peptide: str) -> list[dict]:
    ions, n = [], len(peptide)
    for i in range(1, n):
        b = sum(MONO[a] for a in peptide[:i]) + PROTON
        y = sum(MONO[a] for a in peptide[n - i:]) + WATER + PROTON
        ions.append({"kind": "b", "name": f"b{i}", "seq": peptide[:i], "mz": b, "ordinal": i})
        ions.append({"kind": "y", "name": f"y{i}", "seq": peptide[n - i:], "mz": y, "ordinal": i})
    b2 = sum(MONO[a] for a in peptide[:2]) + PROTON
    ions.append({"kind": "other", "name": "a$_2$", "seq": peptide[:2], "mz": b2 - CO})
    ions.append({"kind": "other", "name": "Im(F)", "seq": "", "mz": MONO["F"] - CO + PROTON})
    return ions


def match(mz: list[float], ions: list[dict]) -> dict[int, tuple[dict, int]]:
    """Observed peak index -> (ion, isotope) for peaks within TOL_PPM of an ion or its M+1."""
    hits = {}
    for ion in ions:
        for iso in (0, 1):
            target = ion["mz"] + iso * ISOTOPE_SPACING
            best = min(range(len(mz)), key=lambda k: abs(mz[k] - target))
            if abs(mz[best] - target) / target * 1e6 <= TOL_PPM:
                hits[best] = (ion, iso)
    return hits


def main(observed_json: str, output_svg: str) -> None:
    obs = json.load(open(observed_json))
    mz, inten = obs["mz"], obs["intensity"]
    top = max(inten)
    rel = [100 * v / top for v in inten]
    hits = match(mz, theoretical_ions(PEPTIDE))

    plt.rcParams.update({
        "svg.fonttype": "none", "font.family": "DejaVu Sans", "font.size": 8,
        "axes.edgecolor": INK_2, "axes.labelcolor": INK, "xtick.color": INK_2,
        "ytick.color": INK_2, "axes.linewidth": 0.8,
    })
    fig, ax = plt.subplots(figsize=(7.0, 4.3))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.12, top=0.71)

    for k, (x, h) in enumerate(zip(mz, rel)):
        colour = COLOURS[hits[k][0]["kind"]] if k in hits else MUTED
        width = 1.4 if k in hits else 0.9
        ax.vlines(x, 0, h, color=colour, lw=width, zorder=3 if k in hits else 2)

    # Labels for monoisotopic matches; nudges keep the a2/y2/b2 cluster legible.
    nudge = {"a$_2$": (-9, 0), "y2": (-5, 3), "b2": (6, 3)}
    for k, (ion, iso) in sorted(hits.items(), key=lambda kv: mz[kv[0]]):
        if iso:
            continue
        dx, dy = nudge.get(ion["name"], (0, 0))
        name = ion["name"]
        if ion["kind"] in "by":
            name = f"{ion['kind']}$_{{{ion['ordinal']}}}$"
        label = name if not ion["seq"] else f"{name}\n{ion['seq']}"
        ax.annotate(label, (mz[k], rel[k]), xytext=(mz[k] + dx, rel[k] + 3 + dy),
                    ha="center", va="bottom", fontsize=7.5, color=INK, linespacing=1.05,
                    textcoords="data", zorder=4)

    # Where the competing peptide's diagnostic y3 should be: nothing observed.
    alt_seq, alt_mz = ALT_Y3
    ax.vlines(alt_mz, 0, 12, color=INK_2, lw=0.8, ls=(0, (2, 2)), zorder=1)
    ax.annotate(f"no peak at {alt_mz:.2f}\n(y$_3$ of {alt_seq})", (alt_mz, 6),
                xytext=(alt_mz + 5, 6), ha="left", va="center", fontsize=7, color=INK_2)

    ax.set_xlim(100, 700)
    ax.set_ylim(0, 112)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Relative intensity (%)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.yaxis.grid(True, color=GRID, lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(length=3, width=0.8)

    draw_ladder(fig, hits)
    fig.text(0.085, 0.965, f"{PEPTIDE}, 2+", fontsize=10, fontweight="bold", color=INK, va="top")
    fig.text(0.085, 0.925,
             f"PXD044934 · scan {obs['scan_number']} · RT {obs['rt_min']:.2f} min · "
             f"precursor m/z {obs['precursor_mz']:.4f} · HCD NCE 30 · fragments matched at {TOL_PPM:g} ppm",
             fontsize=7, color=INK_2, va="top")
    legend = [plt.Line2D([], [], color=COLOURS[k], lw=2, label=lbl)
              for k, lbl in (("b", "b ions"), ("y", "y ions"), ("other", "a$_2$ / immonium"))]
    legend.append(plt.Line2D([], [], color=MUTED, lw=2, label="unassigned"))
    ax.legend(handles=legend, loc="upper right", frameon=False, fontsize=7, ncol=4,
              bbox_to_anchor=(1.0, 1.02), handlelength=1.2, columnspacing=1.2)

    fig.savefig(output_svg)
    fig.savefig(output_svg.rsplit(".", 1)[0] + ".pdf")


def draw_ladder(fig, hits: dict) -> None:
    """Sequence with b (top-left hooks) and y (bottom-right hooks) cleavages observed."""
    seen = {(ion["kind"], ion.get("ordinal")) for ion, iso in hits.values() if iso == 0}
    n = len(PEPTIDE)
    lad = fig.add_axes([0.32, 0.745, 0.42, 0.13])
    lad.set_xlim(-0.5, n - 0.5)
    lad.set_ylim(-1, 1)
    lad.axis("off")
    for i, aa in enumerate(PEPTIDE):
        lad.text(i, 0, aa, ha="center", va="center", fontsize=12, color=INK, family="DejaVu Sans Mono")
    for i in range(1, n):
        x = i - 0.5
        if ("b", i) in seen:
            lad.plot([x, x, x - 0.3], [0, 0.6, 0.6], color=COLOURS["b"], lw=1.4,
                     solid_capstyle="round")
        if ("y", n - i) in seen:
            lad.plot([x, x, x + 0.3], [0, -0.6, -0.6], color=COLOURS["y"], lw=1.4,
                     solid_capstyle="round")


if __name__ == "__main__":
    main(*sys.argv[1:3])
