import warnings

try:
    from . import betafab2, dihedraleditor
except ImportError as ie:
    warnings.warn('GUI is not supported due to an import error: {}'.format(ie.args))
