import logging

logging.basicConfig()
import warnings

try:
    from pymol import cmd
    from pymol.plugins import addmenuitemqt
except ImportError:
    warnings.warn(
        'Cannot import PyMOL: functionality will suffer (you can ignore this if you are just building the documentation).')
    cmd = None

from . import setbetahelix, utils, savegro, hbond, betafab2, betafabgui2, gmxselections


def __init_plugin__(self):
    if cmd is not None:
        addmenuitemqt('BetaFab2', command=betafabgui2.betafab2.run)
        addmenuitemqt('BetaFab2 dihedral editor', command=betafabgui2.dihedraleditor.run)


if cmd is not None:
    cmd.extend('fold_bp', setbetahelix.fold_bp)
    cmd.extend('select_bbb', utils.select_bbb)
    cmd.extend('save_gro', savegro.save_gro)
    cmd.extend('save_g96', savegro.save_g96)
    cmd.extend('save_crd', savegro.save_crd)
    cmd.extend('restrain_hbonds_gmx', hbond.restrain_hbonds_gmx)
    cmd.extend('gmx_beta_backbone_dihedrals_selection', gmxselections.gmx_beta_backbone_dihedrals_selection)
    cmd.extend('label_chains', utils.label_chains)
    cmd.extend('restrain_beta_backbone_dihedrals', utils.restrain_beta_backbone_dihedrals)
