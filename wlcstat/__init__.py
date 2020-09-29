"""Initialization of the wlcstat module

.. moduleauthor:: Andrew Spakowitz <ajspakow@stanford.edu>

"""

import wlcstat.wlcave
import wlcstat.wlcgreen
import wlcstat.wlcstruc
import wlcstat.poly_dyn

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
