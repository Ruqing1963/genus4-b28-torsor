"""
make_figures_b28.py -- vector figures for paper/b28_torsor.tex.

Figure 1 (fig_torsor.pdf): the 64-point coset C of lambda = (1,1,1) over F_2 (C = trace-one set of F_128) in
    symplectic coordinates adapted to a plus-type point P: the point P + sum(a_i e_i + b_i f_i) is drawn at row
    a = (a_1 a_2 a_3), column b = (b_1 b_2 b_3), so that Q_P = sum a_i b_i; the 28 roots of B_28 (eps = +1, minus type)
    are the filled dark discs, the 36 roots of x^8 B_28 + 1 (eps = -1, plus type) the light discs.  Left panel: the
    quadratic form and a symplectic transvection t_v (v a duad = root difference) drawn as arrows k -> k + v on the
    points with omega(k, v) = 1.  Right panel: the Frobenius x -> x^2 (type 7 1) on the same 64 points: one fixed point
    and nine 7-cycles, four inside the roots of B_28 and five inside the plus-type points.
Figure 2 (fig_phases.pdf): the eight normalised Frobenius eigenvalues omega_i / sqrt(q) of the rank-8 Witt sheaf on the
    unit circle for six Frobenius types at q = 4 (data: witt_lpoly_g4.json), lambda_1 = 0 (filled) and lambda_1 = 1 (open).
Usage: python make_figures_b28.py      (after witt_lpoly_g4.py)
"""
import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

from coset_signs import analyse

plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm", "pdf.fonttype": 42,
                     "ps.fonttype": 42})
DARK, LIGHT, GRID = "#1f3a5f", "#d9b99b", "#bbbbbb"
ORB = ["#c0392b", "#2980b9", "#27ae60", "#8e44ad", "#d35400", "#16a085", "#7f8c8d", "#f39c12", "#2c3e50"]


def symplectic_basis(K, omega, Q):
    """basis e1,f1,e2,f2,e3,f3 with omega(e_i,f_j)=delta, other pairings 0, and Q(e_i)=Q(f_i)=0 (Arf 0)"""
    remaining = [k for k in K if k]
    basis = []

    def perp(vs):
        return [k for k in K if all(omega[(k, v)] == 0 for v in vs)]

    space = list(K)
    for _ in range(3):
        cands = [k for k in space if k and Q[k] == 0]
        e = cands[0]
        f = next(k for k in space if omega[(e, k)] == 1 and Q[k] == 0)
        basis += [e, f]
        space = [k for k in space if omega[(k, e)] == 0 and omega[(k, f)] == 0]
    return basis


def coords(k, basis, omega):
    """coordinates (a1,a2,a3,b1,b2,b3) of k: a_i = omega(k, f_i), b_i = omega(k, e_i)"""
    a = [omega[(k, basis[2 * i + 1])] for i in range(3)]
    b = [omega[(k, basis[2 * i])] for i in range(3)]
    return a, b


def figure1():
    res = analyse(4, 1, 1, 1, verbose=False)
    E, C, K, eps, omega = res["E"], res["C"], res["K"], res["eps"], res["omega"]
    P = next(w for w in C if eps[w] == -1)  # a plus-type point
    Q = {k: int(eps[P ^ k] != eps[P]) for k in K}
    basis = symplectic_basis(K, omega, Q)
    # check Q = sum a_i b_i in these coordinates
    pos = {}
    for k in K:
        a, b = coords(k, basis, omega)
        assert Q[k] == (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) % 2
        pos[P ^ k] = (int("".join(map(str, b)), 2), int("".join(map(str, a)), 2))
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.7))
    for ax in axes:
        ax.set_aspect("equal")
        ax.set_xlim(-0.8, 7.8)
        ax.set_ylim(-0.8, 7.8)
        ax.set_xticks(range(8))
        ax.set_yticks(range(8))
        ax.set_xticklabels([format(i, "03b") for i in range(8)], fontsize=7)
        ax.set_yticklabels([format(i, "03b") for i in range(8)], fontsize=7)
        ax.set_xlabel(r"$b=(b_1b_2b_3)$", fontsize=8)
        ax.tick_params(length=0)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.grid(True, color=GRID, lw=0.4, alpha=0.6)
    axes[0].set_ylabel(r"$a=(a_1a_2a_3)$", fontsize=8)

    def draw_points(ax):
        for w, (x, y) in pos.items():
            if eps[w] == 1:
                ax.scatter([x], [y], s=70, c=DARK, edgecolors="black", linewidths=0.4, zorder=3)
            else:
                ax.scatter([x], [y], s=70, c=LIGHT, edgecolors="black", linewidths=0.4, zorder=3)
        x0, y0 = pos[P]
        ax.scatter([x0], [y0], s=160, facecolors="none", edgecolors="black", linewidths=1.0, zorder=4)

    # left: Q = sum a_i b_i and a transvection at a duad v (Q(v)=1)
    ax = axes[0]
    draw_points(ax)
    v = next(k for k in K if Q[k] == 1 and sum(coords(k, basis, omega)[0]) + sum(coords(k, basis, omega)[1]) == 2)
    moved = [k for k in K if omega[(k, v)] == 1]
    drawn = set()
    for k in moved:
        if k in drawn:
            continue
        k2 = k ^ v
        drawn |= {k, k2}
        (x1, y1), (x2, y2) = pos[P ^ k], pos[P ^ k2]
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="<->", mutation_scale=7, lw=0.6, color="#555555",
                                     connectionstyle="arc3,rad=0.25", zorder=2, alpha=0.9))
    ax.set_title(r"(a) $Q_P(k)=\sum a_ib_i$: $\bullet$ roots of $B_{28}$ ($\varepsilon=+1$), $\circ$ roots of $x^8B_{28}+1$;"
                 "\n" + r"arrows: the transvection $t_v$, $v$ a duad", fontsize=7.5)
    # right: Frobenius orbits
    ax = axes[1]
    draw_points(ax)
    seen, orbits = set(), []
    for w in C:
        if w in seen:
            continue
        orb, y = [], w
        while y not in orb:
            orb.append(y)
            y = E.mul(y, y)
        seen |= set(orb)
        orbits.append(orb)
    ci = 0
    for orb in orbits:
        if len(orb) == 1:
            x, y = pos[orb[0]]
            ax.annotate("fixed", (x, y), textcoords="offset points", xytext=(6, 6), fontsize=7)
            continue
        col = ORB[ci % len(ORB)]
        ci += 1
        pts = [pos[w] for w in orb] + [pos[orb[0]]]
        for (x1, y1), (x2, y2) in zip(pts, pts[1:]):
            ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=6, lw=0.7, color=col,
                                         connectionstyle="arc3,rad=0.15", zorder=2, alpha=0.85))
    nroot = sum(1 for orb in orbits if len(orb) == 7 and eps[orb[0]] == 1)
    ax.set_title(f"(b) Frobenius $x\\mapsto x^2$ (type $7\\,1$): one fixed point and nine 7-cycles,\n"
                 f"{nroot} among the roots of $B_{{28}}$, {9 - nroot} among the plus-type points", fontsize=7.5)
    fig.tight_layout()
    fig.savefig("paper/fig_torsor.pdf")
    plt.close(fig)
    print("fig_torsor.pdf written")


def figure2():
    data = json.load(open("witt_lpoly_g4.json", encoding="utf-8"))
    want = ["7 1", "8", "5 2 1", "non-octad(9)", "5 3", "6 1 1"]
    titles = {"non-octad(9)": r"$\Phi_9$ (order 9, no octad)"}
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 5.0))
    for ax, cls in zip(axes.ravel(), want):
        rows = [o for o in data if o["q"] == 4 and o["type"] == cls]
        if not rows:
            ax.set_visible(False)
            continue
        l5, l3 = rows[0]["l5"], rows[0]["l3"]
        rows = [o for o in rows if (o["l5"], o["l3"]) == (l5, l3)]
        ax.set_aspect("equal")
        th = [2 * math.pi * t / 400 for t in range(401)]
        ax.plot([math.cos(t) for t in th], [math.sin(t) for t in th], color=GRID, lw=0.8)
        for d in (8, 12, 16, 24):
            pass
        for o in rows:
            if o["l1"] not in (0, 1):
                continue
            filled = o["l1"] == 0
            for ph in o["phases"]:
                x, y = math.cos(2 * math.pi * ph), math.sin(2 * math.pi * ph)
                ax.scatter([x], [y], s=46 if filled else 60, c=DARK if filled else "none", edgecolors=DARK,
                           linewidths=1.0, zorder=3)
        N = sorted({o["unity_order"] for o in rows if o["l1"] in (0, 1)})
        ax.set_xlim(-1.25, 1.25)
        ax.set_ylim(-1.25, 1.25)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        name = titles.get(cls, f"type ${cls.replace(' ', chr(92) + ',')}$")
        ax.set_title(f"{name}, $\\lambda=(1,{l5},{l3})$ over $\\mathbb{{F}}_4$\n$(\\omega_i/2)^N=1$, $N={','.join(map(str, N))}$",
                     fontsize=7.5)
    fig.suptitle(r"Normalised Frobenius eigenvalues $\omega_i/\sqrt{q}$ of the rank-8 Witt sheaf, $q=4$: "
                 r"$\lambda_1=0$ filled, $\lambda_1=1$ open", fontsize=8)
    fig.tight_layout()
    fig.savefig("paper/fig_phases.pdf")
    plt.close(fig)
    print("fig_phases.pdf written")


if __name__ == "__main__":
    figure1()
    if os.path.exists("witt_lpoly_g4.json"):
        figure2()
    else:
        print("run witt_lpoly_g4.py first for figure 2")
