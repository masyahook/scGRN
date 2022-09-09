import sys

from . import network_analysis as ana

sys.modules.update({f'{__name__}.{m}': globals()[m] for m in ['ana']})  # adding the ability to import using
# shortened alias: e.g. import sc_grn.ana
