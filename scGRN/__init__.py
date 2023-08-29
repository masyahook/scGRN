"""scGRN package."""
import sys

from . import network_analysis as ana
from . import utils as util
from .config import _GGPLOT_COLORS as ggcolors

# adding the ability to import using, shortened alias: e.g. import sc_grn.ana
sys.modules.update({f"{__name__}.{m}": globals()[m] for m in ["ana", "util"]})
