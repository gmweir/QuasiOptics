# This section is to improve python compataibilty between py27 and py3
from __future__ import absolute_import

import sys
from . import ColorPy
sys.modules[__package__] = ColorPy