"""Plot the shared x-plane cut of a hyperelliptic curve."""
from .stokes import Cut


def plot_cut_system(F, centers=(0,), path=None):
    """Draw one branch diagram, marking optional logarithmic zeros."""
    import matplotlib.pyplot as plt
    cut = Cut(F, centers)
    nodes = [complex(z) for z in cut.nodes]
    nodes.append(nodes[-1]+1.3*float(cut.span)*complex(cut.direction))
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot([z.real for z in nodes], [z.imag for z in nodes], '-', color='#267f91', label='Integration cut')
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
