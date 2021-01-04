"""Implement a database of secondary structures using the plugin preference store of PyMOL"""

import typing
import warnings

try:
    from pymol import cmd
    import pymol.plugins
except ImportError:
    warnings.warn(
        'Cannot import PyMOL: functionality will suffer (you can ignore this if you are just building the documentation).')
    cmd = None


class SecondaryStructureDB:
    """A database of the known secondary structures and the corresponding dihedral angles"""
    _instance = None

    DEFAULT_HELIXTYPES = {
        'Z6M': (126.7, 62.6, 152.7),
        'Z6P': (-126.7, -62.6, -152.7),
        'Z8M': (47.5, 53.5, -104.3),
        'Z8P': (-47.5, -53.5, 104.3),
        'H8M': (76.8, -120.6, 52.7),
        'H8P': (-76.8, 120.6, -52.7),
        'H10M': (-77.5, -51.8, -75.1),
        'H10P': (77.5, 51.8, 75.1),
        'H12M': (92.3, -90.0, 104.6),
        'H12P': (-92.3, 90.0, -104.6),
        'H14M': (-140.3, 66.5, -136.8),
        'H14P': (140.3, -66.5, 136.8),
        'SM': (70.5, 176.2, 168.9),
        'SP': (-70.5, -176.2, -168.9),
        'Straight': (180, 180, 180),
        'Straight alpha': (180, None, 180),
        'Alpha-helix': (-57, None, -47),
        '3_10-helix': (-49, None, -26),
        'P-beta-sheet': (-119, None, 113),
        'AP-beta-sheet': (-139, None, 135),
    }

    def __new__(cls):
        if cls._instance is None:
            obj = super().__new__()
            cls._instance = obj
        return cls._instance

    @classmethod
    def getAll(cls, alpha: bool = True, beta: bool = True) -> typing.Dict[str, typing.Tuple[float, float, float]]:
        """Get the list of secondary structures from the PyMOL plugin preferences

        :param alpha: if secondary structures for alpha-peptides are requested
        :type alpha: bool
        :param beta: if secondary structures for beta-peptides are requested
        :type beta: bool
        :return: the known secondary structures and their backbone torsion angles
        :rtype: a dict mapping names to triplets of floats
        """
        ss = pymol.plugins.pref_get('BETAFAB_HELIXTYPES', cls.DEFAULT_HELIXTYPES)
        ssbeta = {k: v for k, v in ss.items() if v[1] is not None}
        ssalpha = {k: v for k, v in ss.items() if v[1] is None}
        dic = {}
        dic.update(ssalpha if alpha else {})
        dic.update(ssbeta if beta else {})
        return dic

    @classmethod
    def add(cls, name: str, phi: float, theta: float, psi: float):
        """Add or edit an entry in the secondary structure torsion angle database

        :param name: the name of the secondary structure
        :type name: str
        :param phi: the first torsion angle
        :type phi: float
        :param theta: the second torsion angle (None for alpha-amino acid secondary structures)
        :type theta: float or None
        :param psi: the third torsion angle
        :type psi: float
        """
        dic = cls.getAll()
        dic[name] = (phi, theta, psi)
        pymol.plugins.pref_set('BETAFAB_HELIXTYPES', dic)
        pymol.plugins.pref_save(quiet=True)

    @classmethod
    def remove(cls, name: str):
        """Remove an entry in the secondary structure torsion angle database

        :param name: the name of the secondary structure
        :type name: str
        """
        dic = cls.getAll()
        try:
            del dic[name]
        except KeyError:
            pass
        pymol.plugins.pref_set('BETAFAB_HELIXTYPES', dic)
        pymol.plugins.pref_save(quiet=True)

    @classmethod
    def find(cls, phi: float, theta: typing.Optional[float], psi: float, tolerance: float = 0.5) -> typing.Optional[
        str]:
        """Try to find a secondary structure in the database which matches the given torsion angles

        :param phi: the first torsion angle
        :type phi: float
        :param theta: the second torsion angle (None for alpha-amino acids)
        :type theta: float or None
        :param psi: the third torsion angle
        :type psi: float
        :param tolerance: absolute tolerance in each angle
        :type tolerance: float
        :return: the name of the closest matching secondary structure or None if no match
        :rtype: str or None
        """
        if theta is None:
            # find only alpha-amino acid sec.structures
            diff = {k: abs(angles[0] - phi) + abs(angles[2] - psi)
                    for k, angles in cls.getAll().items()
                    if angles[1] is None and abs(angles[0] - phi) < tolerance
                    and abs(angles[2] - psi) < tolerance}
        else:
            # find only beta-amino acid sec.structures
            diff = {k: abs(angles[0] - phi) + abs(angles[1] - theta) + abs(angles[2] - psi)
                    for k, angles in cls.getAll().items()
                    if angles[1] is not None and abs(angles[0] - phi) < tolerance
                    and abs(angles[1] - theta) < tolerance and abs(angles[2] - psi) < tolerance}
        if not diff:
            return None
        mindiff = min(diff.values())
        return [k for k in diff if diff[k] == mindiff][0]

    @classmethod
    def dihedrals(cls, name: str) -> typing.Tuple[float, typing.Optional[float], float]:
        """Get the dihedral angles corresponding to a secondary structure.

        :param name: the name of the secondary structure
        :type name: str
        :return: the dihedral angles corresponding to that secondary structure, in degrees
        :rtype: a tuple of three floats, the central one can be None
        :raises KeyError: if the named entry does not exist
        """
        return cls.getAll(alpha=True, beta=True)[name]

    @classmethod
    def addDefaults(cls):
        """Add the default secondary structure types to the list"""
        dic = cls.getAll()
        dic.update(cls.DEFAULT_HELIXTYPES)
        pymol.plugins.pref_set('BETAFAB_HELIXTYPES', dic)
        pymol.plugins.pref_save(quiet=True)


def ssdb_add(entryname: str, phi: typing.Union[str, float], theta: typing.Union[str, float],
             psi: typing.Union[str, float, None] = None, _self=None):
    """
    DESCRIPTION

        Add/modify an entry in the secondary structure database

    USAGE

        ssdb_add entryname, phi, theta, psi   # for beta-amino acids

        or

        ssdb_add entryname, phi, psi  # for alpha-amino acids

    ARGUMENTS

        entryname = string: the name of the entry

        phi = float: the first (N-terminal) backbone dihedral angle: (C-, N, CA, C) for alpha-amino acids and
            (C-, N, CB, CA) for beta-amino acids

        theta = float: the middle backbone dihedral angle for beta-amino acids: (N, CB, CA, C)

        psi = float: the last (C-terminal) backbone dihedral angle: (N, CA, C, N+) for alpha-amino acids and
            (CB, CA, C, N+) for beta-amino acids

    NOTES

        changes are automatically saved for further PyMOL sessions.

    SEE ALSO

        ssdb_del, ssdb_list, ssdb_resetdefaults, ssdb_dihedrals
    """
    if psi is None:
        phi = float(phi)
        psi = float(theta)
        theta = None
    else:
        phi = float(phi)
        theta = float(theta)
        psi = float(psi)
    SecondaryStructureDB.add(entryname, phi, theta, psi)


def ssdb_del(entryname: str, _self=None):
    """
    DESCRIPTION

        Remove an entry from the secondary structure database

    USAGE

        ssdb_del entryname


    ARGUMENTS

        entryname = string: the name of the entry

    NOTES

        changes are automatically saved for further PyMOL sessions.

    SEE ALSO

        ssdb_add, ssdb_list, ssdb_resetdefaults, ssdb_dihedrals
    """
    SecondaryStructureDB.remove(entryname)


def ssdb_list(_self=None):
    """
    DESCRIPTION

        Lists the entries in the secondary structure database

    USAGE

        ssdb_list


    SEE ALSO

        ssdb_add, ssdb_del, ssdb_resetdefaults, ssdb_dihedrals
    """
    dic = SecondaryStructureDB.getAll()
    if not dic:
        print('No entries in the secondary structure database.')
        return
    namelen = max(max([len(e) for e in dic]), len('Secondary structure') + 2)
    numberformat = ' {:7.1f}° '
    numberlen = len(numberformat.format(-360.0))
    sepline = '+' + '-' * (namelen) + '+' + (numberlen * '-' + '+') * 3
    headerline = '|{{:^{}}}|'.format(namelen).format('Secondary structure') + '{{:^{}}}|'.format(numberlen).format(
        'φ') + '{{:^{}}}|'.format(numberlen).format('ϑ') + '{{:^{}}}|'.format(numberlen).format('ψ')
    print(sepline)
    print(headerline)
    print(sepline.replace('-', '='))
    for entry in sorted(dic):
        if dic[entry][1] is not None:
            print('|{{:^{}}}|'.format(namelen).format(entry) + numberformat.format(
                dic[entry][0]) + '|' + numberformat.format(dic[entry][1]) + '|' + numberformat.format(
                dic[entry][2]) + '|')
        else:
            print('|{{:^{}}}|'.format(namelen).format(entry) + numberformat.format(
                dic[entry][0]) + '|' + '{{:^{}}}'.format(numberlen).format('--') + '|' + numberformat.format(
                dic[entry][2]) + '|')
        print(sepline)


def ssdb_resetdefaults(_self=None):
    """
    DESCRIPTION

        Add back the default entries to the secondary structure database

    USAGE

        ssdb_resetdefaults


    NOTES

        changes are automatically saved for further PyMOL sessions.

    SEE ALSO

        ssdb_add, ssdb_del, ssdb_list, ssdb_dihedrals
    """
    SecondaryStructureDB.addDefaults()


def ssdb_dihedrals(entryname: str, _self=None):
    """
    DESCRIPTION

        Print the dihedral angles corresponding to a given entry
            secondary structure database

    USAGE

        ssdb_dihedrals entryname

    ARGUMENTS

        entryname = string: the name of the entry


    SEE ALSO

        ssdb_add, ssdb_del, ssdb_list, ssdb_resetdefaults
    """
    try:
        dihedrals = SecondaryStructureDB.dihedrals(entryname)
    except KeyError:
        print('Entry {} does not exist in the secondary structure database!'.format(entryname))
        return
    print('Dihedral angles for {}:'.format(entryname))
    for anglename, value in zip('φϑψ', dihedrals):
        if value is None:
            continue
        print('    {}: {:8.2f}°'.format(anglename, value))


if cmd is not None:
    cmd.extend('ssdb_add', ssdb_add)
    cmd.extend('ssdb_del', ssdb_del)
    cmd.extend('ssdb_list', ssdb_list)
    cmd.extend('ssdb_dihedrals', ssdb_dihedrals)
    cmd.extend('ssdb_resetdefaults', ssdb_resetdefaults)
