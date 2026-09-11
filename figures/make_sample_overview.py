"""
Target-sample overview figure (clean rebuild).

  x axis  = metallicity [Fe/H]
  y axis  = concentration c = log10(r_t / r_c)
  bubble  = half-light radius r_h (arcmin)
  colour  = none; every cluster drawn in the same light blue

Writes a dark version (matches the app) and a light version (for a
white-background slide deck) into this folder.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).parent
XLSX = HERE.parent / "membership_filter" / "gc_target_list.xlsx"

THEMES = {
    "dark":  dict(bg="#0a0e18", panel="#101827", ink="#e9edf6", dim="#93a0b8",
                  faint="#6b7893", grid="#1c2540", dot="#8fc2ff", edge="#cfe2ff"),
    "light": dict(bg="#ffffff", panel="#ffffff", ink="#1a2233", dim="#44506b",
                  faint="#6b7893", grid="#d7dEEA", dot="#7fb3f0", edge="#2f6db5"),
}

# bubble area scales with the half-light radius
def bubble(r_h):
    return 30 + 115 * np.asarray(r_h, float)


def declutter(ax, xs, ys, labels, colour, fs=8.5, rounds=120, avoid=None):
    """Place each label beside its point, then nudge overlapping labels
    apart vertically until they stop colliding (cheap text repulsion).
    `avoid` is an optional extra box (the legend) labels must stay out of."""
    fig = ax.figure
    fig.canvas.draw()
    # start each label to the right of its point, in offset (point) space
    texts, offs = [], []
    for x, y, l in zip(xs, ys, labels):
        texts.append(ax.annotate(l, (x, y), xytext=(8, 4), textcoords="offset points",
                                 fontsize=fs, color=colour, zorder=5))
        offs.append([8.0, 4.0])

    for _ in range(rounds):
        fig.canvas.draw()
        boxes = [t.get_window_extent() for t in texts]
        moved = False
        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                a, b = boxes[i], boxes[j]
                if not a.overlaps(b):
                    continue
                moved = True
                step = 1.6
                if a.y0 <= b.y0:
                    offs[i][1] -= step; offs[j][1] += step
                else:
                    offs[i][1] += step; offs[j][1] -= step
                texts[i].set_position(tuple(offs[i]))
                texts[j].set_position(tuple(offs[j]))
                boxes[i] = texts[i].get_window_extent()
                boxes[j] = texts[j].get_window_extent()
        # push any label off the legend by flipping it to the other side
        if avoid is not None:
            lb = avoid.get_window_extent()
            for i, t in enumerate(texts):
                if t.get_window_extent().overlaps(lb):
                    moved = True
                    offs[i][0] = -offs[i][0] - t.get_window_extent().width * 0.0
                    texts[i].set_ha("right" if offs[i][0] < 0 else "left")
                    texts[i].set_position(tuple(offs[i]))
        if not moved:
            break
    return texts


def build(df, theme_name):
    t = THEMES[theme_name]
    fig, ax = plt.subplots(figsize=(11, 7.6), dpi=160)
    fig.patch.set_facecolor(t["bg"])
    ax.set_facecolor(t["panel"])
    for s in ax.spines.values():
        s.set_color(t["grid"])
    ax.tick_params(colors=t["dim"], labelsize=10)
    ax.grid(color=t["grid"], lw=0.6, alpha=0.7, linestyle="--")
    ax.set_axisbelow(True)

    ax.scatter(df.feh, df.conc, s=bubble(df.rh),
               c=t["dot"], alpha=0.72, edgecolors=t["edge"], linewidths=0.9, zorder=3)

    ax.set_xlabel("Metallicity  [Fe/H]   (dex)", color=t["ink"], fontsize=12,
                  labelpad=9, fontweight="medium")
    ax.set_ylabel("Concentration  c = log$_{10}$($r_t$ / $r_c$)", color=t["ink"],
                  fontsize=12, labelpad=9, fontweight="medium")
    ax.set_title("The 24-cluster sample: metallicity, concentration and size",
                 color=t["ink"], fontsize=14, pad=14, fontweight="semibold")

    # the only legend: what the bubble size means
    ref = [1, 2, 3, 5]
    handles = [plt.scatter([], [], s=bubble(r), c=t["dot"], alpha=0.72,
                           edgecolors=t["edge"], linewidths=0.9,
                           label=f"{r} arcmin") for r in ref]
    # lower-right is the one genuinely empty corner of this data, so the
    # legend never sits on top of a cluster
    leg = ax.legend(handles=handles, title="Half-light radius $r_h$  (bubble size)",
                    loc="lower right", labelspacing=1.5, borderpad=1.1,
                    handletextpad=1.5, frameon=True, fontsize=9.5)
    leg.get_frame().set_facecolor(t["panel"])
    leg.get_frame().set_edgecolor(t["grid"])
    leg.get_title().set_color(t["ink"])
    leg.get_title().set_fontsize(10)
    for txt in leg.get_texts():
        txt.set_color(t["dim"])

    # extra room on the right so the outermost label (Terzan 5) isn't clipped
    ax.margins(x=0.09, y=0.13)
    x0, x1 = ax.get_xlim()
    ax.set_xlim(x0, x1 + 0.16 * (x1 - x0))
    fig.tight_layout()
    # labels last, once the axes limits and the legend box are final, so the
    # collision test sees the geometry the reader will actually see
    declutter(ax, df.feh.to_numpy(), df.conc.to_numpy(), df.name.tolist(),
              t["faint"], avoid=leg.get_frame())
    out = HERE / f"sample_overview_{theme_name}.png"
    fig.savefig(out, facecolor=t["bg"], bbox_inches="tight")
    plt.close(fig)
    print("wrote", out.name)


def main():
    HERE.mkdir(exist_ok=True)
    raw = pd.read_excel(XLSX, sheet_name="GC Target List", header=3)
    raw = raw.dropna(subset=["Cluster ID", "[Fe/H] (dex)"])
    df = pd.DataFrame({
        "name": raw["Cluster ID"].astype(str).str.strip(),
        "feh": pd.to_numeric(raw["[Fe/H] (dex)"], errors="coerce"),
        "conc": pd.to_numeric(raw["Concentration c = log10(r_t/r_c)"], errors="coerce"),
        "rh": pd.to_numeric(raw["Half-light radius r_h (arcmin)"], errors="coerce"),
    }).dropna()
    print(f"{len(df)} clusters:  [Fe/H] {df.feh.min():.2f}..{df.feh.max():.2f}   "
          f"c {df.conc.min():.2f}..{df.conc.max():.2f}   r_h {df.rh.min():.2f}..{df.rh.max():.2f}'")
    for th in THEMES:
        build(df, th)


if __name__ == "__main__":
    main()
