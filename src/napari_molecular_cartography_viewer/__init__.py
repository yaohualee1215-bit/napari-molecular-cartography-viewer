from ._version import __version__

# Export the napari widget factory referenced in napari.yaml
from ._widget import make_molecular_cartography_viewer_widget

__all__ = [
    "__version__",
    "make_molecular_cartography_viewer_widget",
]
