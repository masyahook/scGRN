import sys

from . import network_analysis as ana
from . import utils as util
from .utils import _GGPLOT_COLORS as ggcolors

sys.modules.update({f'{__name__}.{m}': globals()[m] for m in ['ana', 'util']})  # adding the ability to import using
# shortened alias: e.g. import sc_grn.ana
