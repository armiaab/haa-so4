"""Original design geometry for theory.tex (all coordinates in millimetres).

The command-line entry point now delegates to gambar_photocraft.py, which
authors Figure 4 as gate circuits and exports the editable PhotoCraft sources.
The original drawing functions remain here as a record of the initial geometry.
"""
from argparse import ArgumentParser
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import to_rgba
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, Polygon, Rectangle
from matplotlib.path import Path as MplPath
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "figures"
INK, GRAY, PALE = "#1E293B", "#64748B", "#F1F5F9"
BLUE, OLIVE = "#0288D1", "#689F38"
EMERALD, GREEN = "#059669", "#2E7D32"
WINE, CRIMSON = "#881337", "#DC2626"
AMBER, CORAL = "#F57C00", "#FCA5A5"
DASH = (0, (3, 2))


def configure_fonts():
    # TeX Gyre Heros is the Helvetica-compatible family distributed with TeX.
    available = {f.name for f in font_manager.fontManager.ttflist}
    if not {"Inter", "Helvetica", "TeX Gyre Heros"} & available:
        for filename in ("texgyreheros-regular.otf", "texgyreheros-bold.otf"):
            try:
                located = subprocess.run(
                    ["kpsewhich", filename], capture_output=True, text=True, check=True
                ).stdout.strip()
                if located:
                    font_manager.fontManager.addfont(located)
            except (OSError, subprocess.CalledProcessError):
                pass
    available = {f.name for f in font_manager.fontManager.ttflist}
    family = next((n for n in ("Inter", "Helvetica", "TeX Gyre Heros") if n in available), None)
    if family is None:
        raise RuntimeError("Install Inter, Helvetica, or TeX Gyre Heros before rendering.")
    plt.rcParams.update({
        "font.family": family, "font.size": 8, "mathtext.fontset": "cm",
        "text.color": INK, "axes.labelcolor": INK, "axes.edgecolor": INK,
        "xtick.color": INK, "ytick.color": INK, "axes.linewidth": 1.2,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "path",
        "svg.hashsalt": "haa-so4-theory-2026", "savefig.facecolor": "white",
        "path.simplify": False,
    })
    print(f"Labels: {family}; mathematics: Computer Modern")


def canvas(width, height):
    fig = plt.figure(figsize=(width / 25.4, height / 25.4))
    ax = fig.add_axes([0, 0, 1, 1], xlim=(0, width), ylim=(0, height))
    ax.set_aspect("equal")
    ax.set_axis_off()
    return fig, ax


def text(ax, x, y, label, size=8, color=INK, ha="center", **kwargs):
    return ax.text(x, y, label, fontsize=size, color=color, ha=ha, va="center",
                   zorder=8, **kwargs)


def line(ax, points, color=INK, lw=1.2, ls="-", zorder=2):
    points = np.asarray(points)
    return ax.plot(points[:, 0], points[:, 1], color=color, lw=lw, ls=ls,
                   solid_capstyle="round", zorder=zorder)[0]


def arrow(ax, start, end, color=INK, lw=1.2, style="-|>", scale=8, **kwargs):
    patch = FancyArrowPatch(start, end, arrowstyle=style, mutation_scale=scale,
                            color=color, lw=lw, shrinkA=0, shrinkB=0,
                            zorder=5, **kwargs)
    ax.add_patch(patch)
    return patch


def rect(ax, x, y, w, h, fill="white", edge=INK, lw=1.2, ls="-", zorder=4):
    patch = Rectangle((x, y), w, h, facecolor=fill, edgecolor=edge,
                      lw=lw, ls=ls, zorder=zorder)
    ax.add_patch(patch)
    return patch


def node(ax, x, y, r, fill, edge, lw=1.2, zorder=5):
    ax.add_patch(Circle((x, y), r, fc=fill, ec=edge, lw=lw, zorder=zorder))


def heading(ax, x, y, panel, title, size=8.8):
    text(ax, x, y, f"({panel})", size, weight="bold", ha="left")
    text(ax, x + 6.5, y, title, size, weight="bold", ha="left")


def save(fig, name, preview):
    # Do not use bbox_inches='tight': the requested physical canvas is fixed.
    fig.savefig(OUT / f"{name}.pdf", metadata={
        "Creator": "src/gambar_teori.py", "Title": name,
        "Subject": "Analytic quantum-chemistry methods schematic",
        "CreationDate": None, "ModDate": None,
    })
    fig.savefig(OUT / f"{name}.svg", metadata={
        "Creator": "src/gambar_teori.py", "Title": name, "Date": None,
        "Description": "Vector schematic; mathematical labels outlined as paths.",
    })
    if preview:
        preview.mkdir(parents=True, exist_ok=True)
        fig.savefig(preview / f"{name}.png", dpi=220)
    plt.close(fig)
    print(name)


def workflow():
    fig, ax = canvas(170, 50)
    heading(ax, 2, 46, "a", "Molecular model", 8.3)
    heading(ax, 48, 46, "b", r"Qubit mapping & $\mathbb{Z}_2$ tapering", 8.1)
    heading(ax, 104, 46, "c", "Variational optimization", 8.1)
    xs = [6, 15, 24, 33]
    arrow(ax, (1.5, 27), (39, 27), BLUE, scale=7)
    text(ax, 40, 29, "$z$", 8, BLUE)
    for x in xs:
        node(ax, x, 27, 2.65, "#E1F5FE", BLUE)
        text(ax, x, 27, "H", 8.5, BLUE)
    for a, b in zip(xs[:-1], xs[1:]):
        line(ax, [(a, 31), (a, 36)], BLUE, 0.8)
        arrow(ax, (a + .4, 34), (b - .4, 34), BLUE, .8, "<->", 6)
        text(ax, (a + b) / 2, 37, "$R$", 8)
    line(ax, [(33, 31), (33, 36)], BLUE, .8)
    text(ax, 19.5, 18, r"$\hat{H}$ real", 9)
    text(ax, 19.5, 12, r"Linear H$_4$ / STO-3G", 7.5, GRAY)
    arrow(ax, (40.5, 24.5), (48.5, 24.5), GRAY, .8, scale=7)
    text(ax, 43.5, 17, "Jordan-\nWigner", 7.1, GRAY, linespacing=1.15)

    top = np.linspace(52, 96, 8)
    bottom = np.linspace(59, 91, 5)
    active = (2, 3, 5, 6, 7)
    for q, bx in zip(active, bottom):
        line(ax, [(top[q], 30.3), (bx, 24)], "#9E9E9E", 1, DASH)
    for q, x in enumerate(top):
        removed = q in (0, 1, 4)
        node(ax, x, 32.5, 2.15, "#FFEBEE" if removed else "#F0F4C3",
             "#D32F2F" if removed else OLIVE)
        text(ax, x, 38, str(q), 7.3)
        if removed:
            for sign in (-1, 1):
                line(ax, [(x - 1.2, 32.5 - sign * 1.2),
                          (x + 1.2, 32.5 + sign * 1.2)], "#D32F2F", 1.5, zorder=6)
    for q, x in enumerate(bottom):
        node(ax, x, 22, 2.15, "#F0F4C3", OLIVE)
        text(ax, x, 22, str(q), 6.5)
    text(ax, 75, 16.4, "5 qubits", 7.5)
    text(ax, 75, 11, r"Tapering $\mathbb{Z}_2$", 8)
    text(ax, 75, 5.5, r"$H=\sum_k h_kP_k$ (remains real)", 7.6)

    text(ax, 118, 39, r"$|\Phi_{\mathrm{HF}}\rangle=|00100\rangle$", 8.2)
    rect(ax, 108, 21, 16, 13, "#E8F5E9", GREEN, zorder=2)
    for y in np.linspace(23, 32, 5):
        line(ax, [(104, y), (126, y)], GREEN, 1.2, zorder=3)
    text(ax, 116, 27.5, r"$U(\theta)$", 9, GREEN,
         bbox={"facecolor": "#E8F5E9", "edgecolor": "none", "pad": 1.2})
    arrow(ax, (101.5, 27.5), (106.8, 27.5), GREEN, scale=6)
    arrow(ax, (127, 27.5), (138, 27.5), GRAY, .8, scale=7)
    text(ax, 132.5, 31, r"$\langle H\rangle$", 8)
    text(ax, 153, 40, r"$\min_{\theta}\ E(\theta)$", 8.6)
    text(ax, 153, 35, r"$E(\theta)\geq E_0$", 8.2)
    arrow(ax, (141, 18), (141, 32.5), scale=6)
    arrow(ax, (141, 18), (165, 18), scale=6)
    line(ax, [(141, 20), (165, 20)], GRAY, .8, DASH)
    text(ax, 167, 20.5, "$E_0$", 7)
    t = np.linspace(0, 1, 120)
    line(ax, np.c_[142 + 21 * t, 20.35 + 10.5 * np.exp(-5 * t)], AMBER, 1.5)
    text(ax, 153, 14.5, "Iteration", 7, GRAY)
    feedback = MplPath([(162, 16), (158, 5), (115, 4), (115, 20)],
                       [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
    ax.add_patch(FancyArrowPatch(path=feedback, arrowstyle="-|>", mutation_scale=8,
                                 color=GREEN, lw=1.2, zorder=5))
    text(ax, 137, 4.7, r"L-BFGS-B: $\theta$", 7.6, GREEN)
    return fig


def gate(ax, x, y, label, width=5.6, height=5.5, complex_gate=False, edge=None, size=8):
    color = edge or (GRAY if complex_gate else INK)
    rect(ax, x - width / 2, y - height / 2, width, height,
         PALE if complex_gate else "white", color,
         .8 if complex_gate else 1.2, DASH if complex_gate else "-")
    text(ax, x, y, label, size)


def layer(ax, left, right, top, bottom):
    rect(ax, left, bottom, right - left, top - bottom, "none", GRAY, .8, DASH, zorder=1)
    text(ax, (left + right) / 2, top + 3.4, r"$\times L$", 9)


def register(ax, x0, x1, ys, prep_x, ancilla=False):
    for q, y in enumerate(ys):
        label = r"$a:\,|0\rangle$" if q == 5 else rf"$q_{q}$"
        text(ax, x0 - 2, y, label, 8, ha="right")
        line(ax, [(x0, y), (x1, y)])
    gate(ax, prep_x, ys[2], "$X$", width=4.8, height=5.5)
    if ancilla:
        gate(ax, x1 - 3.3, ys[-1], r"$\mathrm{Tr}_a$", width=6.6, height=6, size=7.3)


def cnot(ax, x, yc, yt):
    line(ax, [(x, yc), (x, yt)], zorder=3)
    node(ax, x, yc, .7, INK, INK, .8)
    node(ax, x, yt, 1.65, "white", INK)
    line(ax, [(x - 1.6, yt), (x + 1.6, yt)], zorder=6)
    line(ax, [(x, yt - 1.6), (x, yt + 1.6)], zorder=6)


def circuits():
    fig, ax = canvas(170, 150)
    heading(ax, 2, 145, "a", "HEA")
    heading(ax, 88, 145, "b", "HAA")
    text(ax, 42, 138.5, r"Secondary $R_YR_Z$ variant; primary baseline: $U_3$", 7.5, GRAY)
    text(ax, 128, 138.5, "Common structure for CAN and SO(4)", 7.5, GRAY)
    ys = np.array([123, 116, 109, 102, 95, 88])
    register(ax, 9, 81, ys[:5], 12)
    layer(ax, 17, 62, 130, 90)
    for y in ys[:5]:
        gate(ax, 22, y, "$R_Y$")
        gate(ax, 30, y, "$R_Z$", complex_gate=True)
        gate(ax, 68, y, "$R_Y$")
        gate(ax, 76, y, "$R_Z$", complex_gate=True)
    for i, x in enumerate((39, 45.5, 52, 58.5)):
        cnot(ax, x, ys[i], ys[i + 1])
    text(ax, 42, 83, r"$X$ on $q_2$ prepares $|\Phi_{\mathrm{HF}}\rangle=|00100\rangle$.", 7.5, GRAY)
    text(ax, 42, 77.5, "The final rotation pair is outside the repeated layer.", 7.2, GRAY)

    register(ax, 97, 168, ys, 100, True)
    layer(ax, 116, 151, 130, 83)
    for y in ys[:5]:
        gate(ax, 109, y, "$V$", width=6)
        gate(ax, 157, y, "$V$", width=6)
    for y in ys:
        gate(ax, 123, y, "$W$", width=6)
    gate(ax, 140, (ys[0] + ys[-1]) / 2, r"$\mathcal{B}_{\mathcal{P}}$",
         width=14, height=ys[0] - ys[-1] + 5.5, size=13)
    text(ax, 128, 77.5, r"CAN: $V=U_3,\ W=I$; SO(4): $V=I,\ W=R_Y$", 8)
    text(ax, 128, 71.5, r"$\mathcal{B}_{\mathcal{P}}$: blocks on the pairs selected below", 7.5, GRAY)
    text(ax, 128, 66, "End rotations have independent parameters.", 7, GRAY)

    heading(ax, 2, 58, "c", "Brick topology")
    heading(ax, 88, 58, "d", "Star topology")
    text(ax, 42, 51.5, r"$n_a=0$ / omit the ancilla from HAA", 7.7, GRAY)
    text(ax, 128, 51.5, r"$n_a=1$ / retain the ancilla in HAA", 7.7, GRAY)
    brick_x = np.array([12, 27, 42, 57, 72])
    for i in range(4):
        first = i % 2 == 0
        line(ax, [(brick_x[i], 36), (brick_x[i + 1], 36)],
             WINE if first else GRAY, 1.2, "-" if first else DASH)
        text(ax, (brick_x[i] + brick_x[i + 1]) / 2, 42 if first else 30,
             "1" if first else "2", 8, WINE if first else GRAY)
    for i, x in enumerate(brick_x):
        node(ax, x, 36, 3.2, "#E1F5FE", BLUE)
        text(ax, x, 36, rf"$q_{i}$", 9)
    text(ax, 42, 22, r"1: $(q_0,q_1),\ (q_2,q_3)$ in parallel", 8.1, WINE)
    text(ax, 42, 16.5, r"2: $(q_1,q_2),\ (q_3,q_4)$ in parallel", 8.1, GRAY)
    text(ax, 42, 10.5, "Two half-layers of neighboring pairs", 7.7)

    system = np.array([[101, 39], [113, 43], [128, 44], [143, 43], [155, 39]])
    ancilla = np.array([128, 24])
    for i, xy in enumerate(system):
        line(ax, [xy, ancilla], INK, 1.2)
        node(ax, *xy, 3.2, "#E1F5FE", BLUE)
        text(ax, *xy, rf"$q_{i}$", 9)
    node(ax, *ancilla, 3.5, "#F0F4C3", OLIVE)
    text(ax, *ancilla, "$a$", 10)
    text(ax, 128, 15, r"$(q_0,a)\ \to\ (q_1,a)\ \to\ \cdots\ \to\ (q_4,a)$", 8.1)
    text(ax, 128, 9.5, "Five sequential blocks on the shared ancilla", 7.7)
    return fig


def weyl():
    fig, ax = canvas(120, 75)
    text(ax, 60, 69.5, r"Folded Weyl chamber: $\pi/4\geq\alpha\geq\beta\geq|\gamma|$", 9.2)
    # Linear projection of the exact vertices in units of pi/4.
    origin = np.array([9., 35.])
    projection = np.array([[30., -22.], [58., 24.], [0., 24.]])
    coords = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 1], [1, 1, -1]])
    vertices = origin + coords @ projection
    identity, cx, swap_plus, swap_minus = vertices
    iswap = origin + np.array([1, 1, 0]) @ projection
    ax.add_patch(Polygon([identity, swap_plus, swap_minus], fc=to_rgba(PALE, .6), ec="none", zorder=0))
    for i, j in ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)):
        if (i, j) == (0, 3):
            # Leave space around the floor annotation on this hidden edge.
            delta = vertices[j] - vertices[i]
            for a, b in ((0, .16), (.58, 1)):
                line(ax, [vertices[i] + a * delta, vertices[i] + b * delta], INK, .8, DASH)
        else:
            line(ax, vertices[[i, j]], INK, 1.2)
    ax.add_patch(Polygon([identity, cx, iswap], fc=to_rgba(CORAL, .3),
                         ec=CRIMSON, lw=1.2, zorder=3))
    for a, b in (([84, 36.7], [91, 47.5]), ([76, 28.3], [88, 22])):
        line(ax, [a, b], GRAY, .8, DASH, zorder=4)
    text(ax, 58, 44, r"Generic $\gamma\ne 0$", 8.3, GRAY)
    text(ax, 58, 39.5, "3 CNOTs [CAN]", 8, GRAY)
    text(ax, 34, 31.4, r"$\gamma=0$", 9, WINE)
    text(ax, 34, 26.8, "at most 2 CNOTs", 8, WINE)
    text(ax, 34, 22.5, "[SO(4)]", 8, WINE)
    for v in vertices:
        node(ax, *v, .8, INK, INK, .8)
    node(ax, *iswap, 1, "white", CRIMSON)
    text(ax, 9, 40, r"$I\otimes I$", 9)
    text(ax, 39, 8.5, "CNOT", 8.5)
    text(ax, 100, 61, "SWAP", 8.5, ha="left")
    text(ax, 100, 13, "SWAP", 8.5, ha="left")
    text(ax, 100, 37, r"$i$SWAP", 8.5, ha="left")
    text(ax, 60, 2.8, "Analytic schematic; no sampled gates", 7.5, GRAY)
    return fig


def tangents():
    fig, ax = canvas(120, 65)
    ax.add_patch(Ellipse((61, 28), 81, 26, angle=7, fc="#E0F2FE", ec=BLUE, lw=1.2, zorder=1))
    start = (39, 30)
    arrow(ax, start, (39, 48), "#334155", scale=10)
    arrow(ax, start, (87, 34), CRIMSON, scale=10)
    node(ax, *start, 1.3, "white", INK, zorder=7)
    text(ax, 39, 59, "Even-Y tangent: imaginary", 8.7, "#334155")
    text(ax, 39, 53, r"$\partial E/\partial\theta_e=0$", 10, "#334155")
    text(ax, 86, 49.5, "Odd-Y tangent: real", 8.6, CRIMSON)
    text(ax, 86, 44.3, "Gradient may be nonzero", 7.8, CRIMSON)
    text(ax, 61, 20, r"$\mathcal{M}_0=\{\boldsymbol{\theta}_e=0\}$", 10, BLUE)
    line(ax, [(23, 31), (37.4, 30)], GRAY, .8)
    text(ax, 13.5, 33.5, "Real circuit", 7.7)
    text(ax, 13.5, 28, r"at $\mathcal{M}_0$", 8)
    line(ax, [(4, 10), (116, 10)], "#CBD5E1", .8)
    text(ax, 31, 5.5, "CAN: 3 CNOTs at generic angles", 7.7, GRAY)
    text(ax, 90, 5.5, "SO(4): at most 2 CNOTs", 7.7, WINE)
    return fig


def symmetry():
    fig, ax = canvas(170, 60)
    heading(ax, 3, 55.5, "a", "All other even-Y angles zero", 8.6)
    heading(ax, 89, 55.5, "b", "Invariant set in exact arithmetic", 8.5)
    line(ax, [(14, 48.5), (20, 48.5)], GRAY, 1.2, DASH)
    text(ax, 22, 48.5, r"Even Y ($\partial E/\partial\theta_e=0$)", 7.8, GRAY, ha="left")
    line(ax, [(14, 43.5), (20, 43.5)], WINE, 1.2)
    text(ax, 22, 43.5, "Odd Y (example)", 7.8, WINE, ha="left")
    chart = fig.add_axes([14 / 170, 13 / 60, 63 / 170, 27 / 60])
    x = np.linspace(-np.pi, np.pi, 400)
    chart.plot(x, np.cos(x), color=GRAY, lw=1.2, ls=DASH)
    chart.plot(x, np.cos(x - .7), color=WINE, lw=1.5)
    chart.axvline(0, color=GRAY, lw=.8, ls=":")
    chart.plot([-.55, .55], [1, 1], color=GRAY, lw=.8)
    chart.plot(0, 1, "o", ms=3.2, mfc="white", mec=GRAY, mew=.8)
    chart.set(xlim=(-np.pi, np.pi), ylim=(-1.15, 1.2),
              xticks=[-np.pi, 0, np.pi], xticklabels=[r"$-\pi$", "$0$", r"$\pi$"], yticks=[])
    chart.spines[["top", "right"]].set_visible(False)
    chart.tick_params(axis="x", width=.8, length=3, pad=2, labelsize=8)
    chart.set_xlabel(r"Angle $\theta_k$", fontsize=8, labelpad=1)
    chart.set_ylabel("Schematic energy", fontsize=8, labelpad=4)

    plane = np.array([[98, 18], [144, 18], [163, 36], [117, 36]])
    ax.add_patch(Polygon(plane, fc="#E0F2FE", ec=BLUE, lw=1.2, zorder=1))
    arrow(ax, (105, 22), (157, 22), GRAY, .8, scale=7)
    arrow(ax, (105, 22), (119, 36), GRAY, .8, scale=7)
    arrow(ax, (105, 22), (105, 48), GRAY, 1.2, scale=8)
    text(ax, 101, 48, r"$\theta_e$", 10)
    text(ax, 160, 21, r"$\boldsymbol{\theta}_o$", 9)
    text(ax, 138, 33, r"$\mathcal{M}_0:\ \boldsymbol{\theta}_e=0$", 9, BLUE)
    points = np.array([[111, 23], [118, 26], [125, 25.2], [132, 27], [139, 28], [145, 28.6]])
    line(ax, points, EMERALD, 1.5, zorder=6)
    for p in points:
        node(ax, *p, .85, "white", EMERALD, 1.2, zorder=7)
    arrow(ax, points[-2], points[-1], EMERALD, 1.2, scale=7)
    text(ax, 132, 10, "HAA-CAN: real product states", 8.3)
    return fig


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--preview", type=Path, help="Optional directory for inspection PNGs")
    args = parser.parse_args()
    configure_fonts()
    OUT.mkdir(exist_ok=True)
    for draw, name in (
        (workflow, "fig1_vqe_pipeline"),
        (circuits, "fig4_ansatze_circuits"),
        (weyl, "fig5_weyl_chamber"),
        (tangents, "fig6_tangent_directions"),
        (symmetry, "fig7_symmetry_manifold"),
    ):
        save(draw(), name, args.preview)


if __name__ == "__main__":
    from gambar_photocraft import main as photocraft_main
    photocraft_main()
