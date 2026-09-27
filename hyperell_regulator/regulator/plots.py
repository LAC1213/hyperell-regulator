r"""
Drawing the cut system the Stokes integrator works on.

One picture of the `w`-plane: the critical values of `\varphi` in the order the
integrator puts them, the path `\Gamma` joining `\varphi(P) = 0` through all of
them to `\varphi(Q) = \infty`, the gluing permutation of the `N` sheets across
each of its segments, and the racetrack that reads those permutations off.

`\Gamma` is also the branch cut of `\log\varphi`: it joins the zero of
`\varphi` to its pole, so the slit plane is simply connected and misses both,
and the logarithm is single valued on each sheet with a jump of `2\pi i`
across every edge.  The picture therefore shows the cut system and the log
branch at once, which is why they are not drawn separately.
"""

import numpy as np

from .cover import Cover
from .stokes import cut_data, gluing_permutations, racetrack

__all__ = ["plot_cut_system"]

INK = "#1c1c1c"
MUTED = "#8a8a8a"
CUT = "#1f6fb4"        # Gamma, and so the branch cut of log(phi)
TRACK = "#c4622d"      # the racetrack
CRIT = "#1c1c1c"


def _cycles(perm):
    seen, out = set(), []
    for s in range(len(perm)):
        if s in seen:
            continue
        c, t = [s], perm[s]
        seen.add(s)
        while t != s:
            c.append(t)
            seen.add(t)
            t = perm[t]
        if len(c) > 1:
            out.append("(" + " ".join(str(v + 1) for v in c) + ")")
    return "".join(out) or "id"


def plot_cut_system(F, phi, path=None, reach=1.35, show_racetrack=True,
                    figsize=(8.4, 7.2), dpi=160, eps=0.02):
    r"""
    Draw the critical values of `\varphi`, the cut system, and the log branch.

    ``F`` is a polynomial over `\QQ` and ``phi`` the pair `(A, B)` of rational
    functions with `\varphi = A(x) + B(x)y`.  With ``path``
    the figure is written there and the path returned; otherwise the figure is
    returned.

    EXAMPLES::

        sage: from hyperell_regulator.regulator import plot_cut_system
        sage: from hyperell_regulator.regulator.polynomials import ring
        sage: import os, tempfile
        sage: x = ring().gen()
        sage: out = os.path.join(tempfile.mkdtemp(), "cut.png")
        sage: bool(os.path.getsize(plot_cut_system(x^5 + x^2 + 2*x + 1,
        ....:                                      (x + 1, 1), path=out)) > 5000)
        True
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cov = Cover(F, *phi)
    nodes, kinds, d = cov.cut_path()
    upper, lower, _, dlam = cut_data(cov, nodes, d, eps=eps)
    tau = gluing_permutations(upper, lower)
    track, _, _ = racetrack(nodes, d, eps=eps)
    crit = [i for i, k in enumerate(kinds) if k != "none"]
    z = nodes[crit]

    span = max(np.abs(z - z.mean()).max(), 1e-3)
    lo = z.mean() - reach * span * (1 + 1j)
    hi = z.mean() + reach * span * (1 + 1j)
    far = z[-1] + d * 4 * span

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.set_facecolor("white")
    ax.axhline(0, color=MUTED, lw=0.6, alpha=0.35, zorder=0)
    ax.axvline(0, color=MUTED, lw=0.6, alpha=0.35, zorder=0)

    if show_racetrack:
        inside = ((track.real > lo.real) & (track.real < hi.real)
                  & (track.imag > lo.imag) & (track.imag < hi.imag))
        t = np.where(inside, track, np.nan)
        ax.plot(t.real, t.imag, color=TRACK, lw=1.0, ls=(0, (4, 3)),
                label="racetrack (reads off the sheet labels)", zorder=2)

    seg = np.append(nodes, far)
    ax.plot(seg.real, seg.imag, color=CUT, lw=2.0, solid_capstyle="round",
            label=r"$\Gamma$ = cut system = branch cut of $\log\varphi$",
            zorder=3)
    ax.annotate("", xy=(far.real, far.imag),
                xytext=(z[-1].real, z[-1].imag),
                arrowprops=dict(arrowstyle="-|>", color=CUT, lw=2.0), zorder=3)

    ax.scatter(z.real, z.imag, s=52, facecolor="white", edgecolor=CRIT,
               linewidth=1.6, zorder=5)
    ax.scatter([z[0].real], [z[0].imag], s=52, facecolor=CRIT,
               edgecolor=CRIT, linewidth=1.6, zorder=6)

    off = 0.055 * span
    for i, v in enumerate(z):
        lab = r"$z_1=\varphi(P)=0$" if i == 0 else r"$z_%d$" % (i + 1)
        ax.annotate(lab, (v.real, v.imag), xytext=(v.real + off, v.imag + off),
                    fontsize=10, color=INK, zorder=7,
                    path_effects=_halo(plt))
    for i in range(len(z)):
        a = z[i]
        b = far if i == len(z) - 1 else z[i + 1]
        # kept off the middle and pushed well clear, so it does not land on a
        # vertex label when consecutive critical values are close
        m = a + (0.42 if i < len(z) - 1 else 0.3) * (b - a)
        # opposite side from the vertex labels, which sit up and to the right
        n = -1j * (b - a) / abs(b - a) * 0.11 * span
        ax.annotate(r"$\tau_%d=$%s" % (i + 1, _cycles(tau[i])),
                    (m.real + n.real, m.imag + n.imag), fontsize=9,
                    color=CUT, ha="center", zorder=7,
                    path_effects=_halo(plt))

    ax.set_xlim(lo.real, hi.real)
    ax.set_ylim(lo.imag, hi.imag)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$\mathrm{Re}\,w$", color=INK)
    ax.set_ylabel(r"$\mathrm{Im}\,w$", color=INK)
    ax.set_title(r"Cut system of $\varphi$: degree %d, %d finite critical "
                 r"values%s, jump of $\log\varphi$ = %+.4f$\,i$"
                 % (cov.N, len(z),
                    "" if len(nodes) == len(z) else " (closing step routed)",
                    dlam.imag), fontsize=11, color=INK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.legend(loc="lower left", frameon=False, fontsize=9)
    fig.tight_layout()
    if path is None:
        return fig
    fig.savefig(path, facecolor="white")
    plt.close(fig)
    return path


def _halo(plt):
    import matplotlib.patheffects as pe
    return [pe.withStroke(linewidth=2.6, foreground="white")]
