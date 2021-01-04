"""Utilities for creating selection files for GROMACS"""


def gmx_beta_backbone_dihedrals_selection(outputfile: str, startresi: int, endresi: int, nterm_ace: bool = True,
                                          cterm_nme: bool = True):
    """Write a selection file for GROMACS containing the backbone torsion angles of a beta-peptide

    :param outputfile: output file name
    :type outputfile: str
    :param startresi: first residue number
    :type startresi: int
    :param endresi: last residue number
    :type endresi: int
    """
    startresi = int(startresi)
    endresi = int(endresi)
    with open(outputfile, "wt") as f:
        for i in range(startresi, endresi + 1):
            if i > startresi:
                f.write('"phi({0})" '
                        '(resnr {1} and name C) plus '
                        '(resnr {0} and name N) plus '
                        '(resnr {0} and (name CB or name CB1)) plus '
                        '(resnr {0} and name CA)\n'.format(i, i - 1))
            f.write('"theta({0})" '
                    '(resnr {0} and name N) plus '
                    '(resnr {0} and (name CB or name CB1)) plus '
                    '(resnr {0} and name CA) plus '
                    '(resnr {0} and name C)\n'.format(i, i - 1))
            if i < endresi:
                f.write('"psi({0})" '
                        '(resnr {0} and (name CB or name CB1)) plus '
                        '(resnr {0} and name CA) plus '
                        '(resnr {0} and name C) plus '
                        '(resnr {1} and name N)\n'.format(i, i + 1))
                f.write('"omega({0})" '
                        '(resnr {0} and name CA) plus '
                        '(resnr {0} and name C) plus '
                        '(resnr {1} and name N) plus '
                        '(resnr {1} and (name CB or name CB1))\n'.format(i, i + 1))
