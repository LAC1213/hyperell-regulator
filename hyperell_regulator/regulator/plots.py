"""Plot the shared x-plane cut of a hyperelliptic curve."""
from .stokes import Cut


def plot_cut_system(F, centers=(0,), path=None):
    """Draw one branch diagram, marking optional logarithmic zeros."""
    import matplotlib.pyplot as plt
    cut = Cut(F, centers)
    fig, ax = plt.subplots(figsize=(8, 6))
    for i, (a, b) in enumerate(cut.edges):
        start = complex(cut.nodes[a])
        end = complex(cut.nodes[b]) if b is not None else start+1.3*float(cut.span)*complex(cut.direction)
        ax.plot([start.real, end.real], [start.imag, end.imag], '-', color='#267f91',
                label='Integration cut' if i == 0 else None)
    roots = [complex(r) for r in cut.roots]
    ax.scatter([z.real for z in roots], [z.imag for z in roots], color='#bd5034', label='Branch values', zorder=3)
    for a in cut.centers:
        ax.plot(float(a), 0, 'o', mfc='white', mec='#267f91')
        ax.annotate(str(a), (float(a), 0), xytext=(5, 5), textcoords='offset points')
    ax.set(xlabel='Re(x)', ylabel='Im(x)', title='Hyperelliptic map x: two sheets')
    ax.set_aspect('equal', adjustable='datalim')
    ax.legend()
    if path is not None:
        fig.savefig(path, bbox_inches='tight')
    return fig, ax
