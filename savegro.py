# coding: utf-8

import math
import warnings

try:
    from pymol import cmd
    import chempy
except ImportError:
    warnings.warn(
        'Cannot import PyMOL: functionality will suffer (you can ignore this if you are just building the documentation).')


def get_box_params(selection):
    try:
        a, b, c, alpha, beta, gamma, spacegroup = cmd.get_symmetry(selection)
    except TypeError:
        # this happens if get_symmetry() is None
        return 0, 0, 0, 0, 0, 0, 0, 0, 0
    alpha *= math.pi / 180.
    beta *= math.pi / 180.
    gamma *= math.pi / 180.
    a /= 10
    b /= 10
    c /= 10
    v1x = a
    v2x = b * math.cos(gamma)
    v2y = b * math.sin(gamma)
    v3x = c * math.cos(beta)
    v3y = c * (math.cos(alpha) - math.cos(gamma) * math.cos(beta)) / (math.sin(gamma))
    v3z = (c ** 2 - v3x ** 2 - v3y ** 2) ** 0.5
    return (v1x, v2y, v3z, 0, 0, v2x, 0, v3x, v3y)


def save_gro(filename, selection='(all)', bx=None, by=None, bz=None, sort=False):
    """
    DESCRIPTION
    
        "save_gro" saves content to a .gro file in a similar fashion as "save"

    USAGE

        save_gro filename [, selection [, bx [, by [, bz ]]]]

    ARGUMENTS
        
        filename = string: file path to be written

        selection = string: atoms to save {default: (all)}
        
        bx = float: rectangular box size X dimension in Angströms
        
        by = float: rectangular box size Y dimension in Angströms

        bz = float: rectangular box size Z dimension in Angströms
        
    NOTES

        In contrast to "save", only the current state can be saved
        as the GRO format is a single-structure format.
        
        If any of bx, by and bz is not supplied, the box will be 
        determined from the get_symmetry() command.
        
        Atoms are written with increasing ID.

    SEE ALSO

        load, save

    """
    if isinstance(sort, str):
        sort = eval(sort)
    with open(filename, 'wt') as f:
        f.write('{}\n'.format(selection))
        model = cmd.get_model(selection)
        f.write('{:>5d}\n'.format(len(model.atom)))
        if sort:
            sortfcn = lambda a: (a.chain, a.resi_number, model.atom.index(a))
        else:
            sortfcn = lambda a: model.atom.index(a)
        for i, atom in enumerate(sorted(model.atom, key=sortfcn)):
            assert isinstance(atom, chempy.Atom)
            f.write('{:>5d}{:<5s}{:>5s}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}\n'.format(atom.resi_number, atom.resn, atom.name,
                                                                                i + 1,
                                                                                float(atom.coord[0]) * 0.1,
                                                                                float(atom.coord[1]) * 0.1,
                                                                                float(atom.coord[2]) * 0.1))
        if bx is None and by is None and bz is None:
            f.write(' {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f}\n'.format(
                *get_box_params(selection)
            ))
        else:
            f.write(' {} {} {}\n'.format(bx, by, bz))
    print('Saved {} atoms.'.format(len(model.atom)))


def save_g96(filename, selection='(all)', bx=None, by=None, bz=None, sort=False):
    """
    DESCRIPTION

        "save_g96" saves content to a .g96 file in a similar fashion as "save"

    USAGE

        save_g96 filename [, selection [, bx [, by [, bz ]]]]

    ARGUMENTS

        filename = string: file path to be written

        selection = string: atoms to save {default: (all)}

        bx = float: rectangular box size X dimension in Angströms

        by = float: rectangular box size Y dimension in Angströms

        bz = float: rectangular box size Z dimension in Angströms

    NOTES

        In contrast to "save", only the current state can be saved
        as the GROMOS-96 format is a single-structure format.

        If any of bx, by and bz is not supplied, the box will be
        determined from the get_symmetry() command.

        Atoms are written with increasing ID.

    SEE ALSO

        load, save

    """
    if isinstance(sort, str):
        sort = eval(sort)
    with open(filename, 'wt') as f:
        f.write('TITLE\n')
        f.write(selection + '\n')
        f.write('END\n')
        f.write('POSITION\n')
        model = cmd.get_model(selection)
        if sort:
            sortfcn = lambda a: (a.chain, a.resi_number, model.atom.index(a))
        else:
            sortfcn = lambda a: model.atom.index(a)
        for i, atom in enumerate(sorted(model.atom, key=sortfcn)):
            f.write('{:>05d} {:<5s} {:<5s} {:>06d} {:>14.9f} {:>14.9f} {:>14.9f}\n'.format(
                atom.resi_number, atom.resn, atom.name, i + 1,
                                                        atom.coord[0] * 0.1, atom.coord[1] * 0.1, atom.coord[2] * 0.1))
        f.write('END\n')
        f.write('BOX\n')
        if bx is None and by is None and bz is None:
            f.write(
                ' {:>14.9f} {:>14.9f} {:>14.9f} {:>14.9f} {:>14.9f} {:>14.9f} {:>14.9f} {:>14.9f} {:>14.9f}\n'.format(
                    *get_box_params(selection)
                ))
        else:
            f.write(' {:>14.9f} {:>14.9f} {:>14.9f}\n'.format(bx, by, bz))
        f.write('END\n')


def save_crd(filename, selection='(all)', sort=True):
    """
    DESCRIPTION

        "save_crd" saves content to an expanded CHARMM CRD file in a similar fashion as "save"

    USAGE

        save_crd filename [, selection ]

    ARGUMENTS

        filename = string: file path to be written

        selection = string: atoms to save {default: (all)}

    NOTES

        In contrast to "save", only the current state can be saved
        as the CRD format is a single-structure format.

        Atoms are written with increasing ID.

    SEE ALSO

        load, save

    """
    if isinstance(sort, str):
        sort = eval(sort)
    with open(filename, 'wt') as f:
        f.write('* COORDINATES\n')
        f.write('* {}\n'.format(selection))
        model = cmd.get_model(selection)
        f.write('{:>10d}\n'.format(len(model.atom)))
        if sort:
            sortfcn = lambda a: (a.chain, a.resi_number, model.atom.index(a))
        else:
            sortfcn = lambda a: model.atom.index(a)
        for i, atom in enumerate(sorted(model.atom, key=sortfcn)):
            assert isinstance(atom, chempy.Atom)
            f.write('{:>10d}{:>10d}  {:<8s}  {:<8s}{:>20.10f}{:>20.10f}{:>20.10f}  {:<8s}  {:<8s}{:>20.10f}\n'.format(
                i, atom.resi_number, atom.resn, atom.name,
                float(atom.coord[0]),
                float(atom.coord[1]),
                float(atom.coord[2]),
                atom.chain, str(atom.resi_number), 0.0
            ))
    print('Saved {} atoms.'.format(len(model.atom)))
