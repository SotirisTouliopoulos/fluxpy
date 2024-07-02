__author__ = "Haris Zafeiropoulos"
__version__ = "0.0.1"

from fluxpy import utils
from fluxpy import stats
from fluxpy import illustrations

from .auto_decorator import auto_decorate

auto_decorate()

__all__ = [
    'illustrations',
    'mappings',
    'stats',
    'utils',
]