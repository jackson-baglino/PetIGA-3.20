from pathlib import Path
import matplotlib.pyplot as plt

# Path to your style file, relative to this module
_STYLE_PATH = Path(__file__).with_name("plot_style.mplstyle")

def apply_style():
    """Load the project-wide Matplotlib style."""
    plt.style.use(_STYLE_PATH)

def new_fig_ax(aspect=None, nrows=1, ncols=1, sharex=False, sharey=False):
    """Create a new figure + axes with consistent defaults."""
    fig, ax = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        sharex=sharex,
        sharey=sharey,
        constrained_layout=True,
    )
    if aspect is not None:
        # Works whether ax is a single Axes or array of Axes
        try:
            for a in ax.ravel():
                a.set_aspect(aspect)
        except AttributeError:
            ax.set_aspect(aspect)
    return fig, ax

def label_axes(ax, xlabel=None, ylabel=None, title=None):
    """Apply consistent labels and title to an Axes."""
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

def save_figure(fig, path, formats=("png", "pdf")):
    """Save a figure to one or more formats with consistent settings."""
    path = Path(path)
    for ext in formats:
        fig.savefig(path.with_suffix("." + ext))