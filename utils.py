import warnings

try:
    from pymol import cmd
    import chempy.models
except ImportError:
    warnings.warn(
        'Cannot import PyMOL: functionality will suffer (you can ignore this if you are just building the documentation).')
    chempy = None

import random
import logging
import time

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


class tempObjectName:
    """Create a temporary object in PyMOL and delete it after use.

    Use this as a context manager, i.e.:

    >>> with tempObjectName('_temporary') as name:
    ...     # do something with `name`
    ...     pass
    >>>

    Name clashes with already existing objects are avoided.
    """
    tempobjnames = []

    def __init__(self, prefix: str = '_temporary', nocreate: bool = True):
        self._prefix = prefix
        self._objname = None
        self._nocreate = nocreate

    def __enter__(self) -> str:
        i = 0
        while True:
            objname = self._prefix + str(i) + str(time.monotonic())
            if objname not in cmd.get_object_list() and objname not in type(self).tempobjnames:
                break
            i += 1
        self._objname = objname
        type(self).tempobjnames.append(self._objname)

        if not self._nocreate:
            logger.debug('Creating model: {}'.format(self._objname))
            cmd.create(self._objname, 'none')
            logger.debug('Object list: {}'.format(cmd.get_object_list()))
            logger.debug('Done.')
        return self._objname

    def __exit__(self, exc_type, exc_val, exc_tb):
        cmd.delete('model {}'.format(self._objname))
        try:
            type(self).tempobjnames.remove(self._objname)
        except ValueError:
            pass
        self._objname = None


class fetchModel:
    """Fetch a model from PyMOL for low-level manipulation

    This is a context manager, i.e.:

    >>> with fetchModel('modelname') as model:
    ...     # manipulate model
    ...     pass
    >>>

    Your updates to the model are imported back into PyMOL. The mechanism uses :func:`pymol.cmd.get_model()` and
    :func:`pymol.cmd.load_model()` for fetching and updating the model, respectively.
    """

    def __init__(self, modelname: str):
        self._modelname = modelname
        self._model = None

    def __enter__(self) -> "chempy.models.Indexed":
        self._model = cmd.get_model('model {}'.format(self._modelname))
        logger.debug('Fetched model {}, containing {} atoms:'.format(self._modelname, len(self._model.atom)))
        for at in self._model.atom:
            logger.debug('   Name: {}\tResi: {}\tResn:{}'.format(at.name, at.resi_number, at.resn))
        return self._model

    def __exit__(self, exc_type, exc_val, exc_tb):
        cmd.delete('model {}'.format(self._modelname))
        cmd.load_model(self._model, self._modelname)
        # read back model for control
        self._model = cmd.get_model('model {}'.format(self._modelname))
        logger.debug('Putting back model {}, containing {} atoms:'.format(self._modelname, len(self._model.atom)))
        for at in self._model.atom:
            logger.debug('   Name: {}\tResi: {}\tResn:{}'.format(at.name, at.resi_number, at.resn))
        self._model = None


def getAtom(model, **kwargs):
    atoms = [a for a in model.atom if all([getattr(a, k) == v for k, v in kwargs.items()])]
    if len(atoms) > 1:
        raise ValueError('More than one atoms match the criteria')
    elif not atoms:
        raise ValueError('No atom matches the criteria')
    else:
        return atoms[0]


def get_atom(model, atom_or_index_or_name):
    #    print('get_atom: atom_or_index_or_name is ',atom_or_index_or_name)
    if isinstance(atom_or_index_or_name, chempy.Atom):
        return atom_or_index_or_name
    elif isinstance(atom_or_index_or_name, str):
        atoms = [a for a in model.atom if a.name == atom_or_index_or_name]
    else:
        atoms = [a for a in model.atom if a.index == atom_or_index_or_name]
    assert len(atoms) <= 1
    if not atoms:
        raise ValueError('No atom with index or name {} in the current model'.format(atom_or_index_or_name))
    return atoms[0]


def getBond(model, atom1, atom2):
    idx1 = model.atom.index(atom1)
    idx2 = model.atom.index(atom2)
    b = [b for b in model.bond if b.index == [idx1, idx2] or b.index == [idx2, idx1]][0]
    return b


def getBondedNeighbours(model, atom=None, **kwargs):
    if atom is None:
        atom = getAtom(model, **kwargs)
    idx = model.atom.index(atom)
    neighborindices = [[at for at in b.index if at != idx][0] for b in model.bond if idx in b.index]
    return [model.atom[ni] for ni in neighborindices]


def newAtom(**kwargs):
    atom = chempy.Atom()
    return updateAtom(atom, **kwargs)


def updateAtom(atom, **kwargs):
    for k, v in kwargs.items():
        setattr(atom, k, v)
    return atom


def addHydrogens(model):
    objname = '__tmp{}'.format(random.randint(0, 100000))
    cmd.delete(objname)
    cmd.load_model(model, objname)
    cmd.h_add('model {}'.format(objname))
    model = cmd.get_model(objname)
    cmd.delete(objname)
    return model


def set_dihedral(modelname, nameandresi1, nameandresi2, nameandresi3, nameandresi4, value):
    """A safer version of cmd.set_dihedral(): doesn't choke on non-singleton selections"""
    selections = [
        '({}) and name {} and resi {}'.format(modelname, *nr)
        for nr in [nameandresi1, nameandresi2, nameandresi3, nameandresi4]
    ]
    if not all([cmd.count_atoms(s) == 1 for s in selections]):
        return False
    cmd.set_dihedral(*(selections + [value]))
    cmd.unpick()
    return True


def get_dihedral(modelname, nameandresi1, nameandresi2, nameandresi3, nameandresi4):
    """A safer version of cmd.get_dihedral(): doesn't choke on non-singleton selections"""
    selections = [
        '({}) and name {} and resi {}'.format(modelname, *nr)
        for nr in [nameandresi1, nameandresi2, nameandresi3, nameandresi4]
    ]
    if not all([cmd.count_atoms(s) == 1 for s in selections]):
        return None
    return cmd.get_dihedral(*selections)


def planarize_peptide_bond(model, resi_name_center, resi_name_fix1, resi_name_fix2, resi_name_movable):
    try:
        ca = getAtom(model, resi_number=resi_name_fix1[0], name=resi_name_fix1[1])
        c = getAtom(model, resi_number=resi_name_center[0], name=resi_name_center[1])
        o = getAtom(model, resi_number=resi_name_movable[0], name=resi_name_movable[1])
        n = getAtom(model, resi_number=resi_name_fix2[0], name=resi_name_fix2[1])
    except ValueError:
        return
    v_c_ca = [x - y for x, y in zip(ca.coord, c.coord)]
    l_c_ca = sum([v ** 2 for v in v_c_ca]) ** 0.5
    v_c_ca_0 = [v / l_c_ca for v in v_c_ca]
    v_c_n = [x - y for x, y in zip(n.coord, c.coord)]
    l_c_n = sum([v ** 2 for v in v_c_n]) ** 0.5
    v_c_n_0 = [v / l_c_n for v in v_c_n]
    l_c_o_old = sum([(x - y) ** 2 for x, y in zip(o.coord, c.coord)]) ** 0.5
    v_c_o_new = [-(x + y) for x, y in zip(v_c_n_0, v_c_ca_0)]
    l_c_o_new = sum([v ** 2 for v in v_c_o_new]) ** 0.5
    v_c_o_new = [v / l_c_o_new * l_c_o_old for v in v_c_o_new]
    o.coord = tuple([x + y for x, y in zip(v_c_o_new, c.coord)])


def extend_peptide(model):
    objname = '__tmp{}'.format(random.randint(0, 100000))
    minresi = min([a.resi_number for a in model.atom])
    maxresi = max([a.resi_number for a in model.atom])

    def is_beta_residue(model, r):
        resn = getAtom(model, resi_number=r, name='N').resn
        return resn.startswith('B2') or resn.startswith('B3')

    # first planarize the atoms around C and N in the peptide bonds.
    for r in range(minresi, maxresi + 1):
        # now align the oxygens: CA, C, O and N must be planar
        planarize_peptide_bond(model, (r, 'C'), (r, 'CA'), (r + 1, 'N'), (r, 'O'))
        if is_beta_residue(model, r):
            # align the hydrogens: N, CB, C and H must be planar
            planarize_peptide_bond(model, (r, 'N'), (r, 'CB'), (r - 1, 'C'), (r, 'H'))
            # try it with CB1
            planarize_peptide_bond(model, (r, 'N'), (r, 'CB1'), (r - 1, 'C'), (r, 'H'))
        else:
            planarize_peptide_bond(model, (r, 'N'), (r, 'CA'), (r - 1, 'C'), (r, 'H'))
    cmd.delete(objname)
    cmd.load_model(model, objname)
    # set all torsions to straight
    for r in range(minresi, maxresi + 1):
        # planarize the peptide bond
        # print('Setting dihedrals of residue #{}'.format(r))
        set_dihedral(objname, ('O', r - 1), ('C', r - 1), ('N', r), ('H', r), 180)
        if is_beta_residue(model, r):
            # the "phi" torsion angle
            # print('This is a beta residue')
            set_dihedral(objname, ('H', r), ('N', r), ('CB+CB1', r), ('CA', r), 0)
            # the "theta" torsion angle
            set_dihedral(objname, ('N', r), ('CB+CB1', r), ('CA', r), ('C', r), 180)
            # the "psi" torsion angle
            set_dihedral(objname, ('CB+CB1', r), ('CA', r), ('C', r), ('O', r), 0)
        else:
            #            print('This is an alpha residue')
            # the "phi" dihedral
            set_dihedral(objname, ('H', r), ('N', r), ('CA', r), ('C', r), 0)
            # the "psi" dihedral
            set_dihedral(objname, ('N', r), ('CA', r), ('C', r), ('O', r), 0)
        # fix the hydrogens
        for atom in ['CA', 'CB', 'CB1']:
            if cmd.count_atoms('model {} and resi {} and name {}'.format(objname, r, atom)) == 1:
                cmd.h_fix('model {} and resi {} and name {}'.format(objname, r, atom))
        cmd.unpick()
    model = cmd.get_model(objname)
    cmd.delete(objname)
    return model


def select_bbb(selectionname='bbone', originalselection='all'):
    """
    DESCRIPTION

    Select the backbone of beta-peptides

    USAGE

    select_beta_backbone name [, originalselection]

    ARGUMENTS

    name = name of the new selection

    originalselection = superset in which the backbone will be searched
    """
    cmd.select(selectionname, '({}) and name CA+CB+CB1+CC+C+O+N+HN'.format(originalselection))


def label_chains(selection):
    """
    DESCRIPTION

    Detect and label chains

    USAGE

    label_chains selection

    ARGUMENTS

    selection = name of the selection on which to operate. Atoms outside the selection
                won't be changed
    """
    chainlabels = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    cmd.alter(selection, 'chain=""')
    for label in chainlabels:
        space = {'stored': []}
        cmd.iterate('({}) and (chain "")'.format(selection), 'stored.append(index)', space=space)
        unnamedindices = space['stored']
        if not unnamedindices:
            break
        cmd.alter('({}) and (bymol (({}) and index {}))'.format(selection, selection, unnamedindices[0]),
                  'chain="{}"'.format(label))
        cmd.sort()


def restrain_beta_backbone_dihedrals(restraintfile, selection='all', fc=10000):
    """
    DESCRIPTION

    Write a dihedral restraint .itp file for GROMACS to restraint the backbone dihedrals of a beta-peptide

    USAGE

    restrain_beta_backbone_dihedrals restraintfile [, selection]

    ARGUMENTS

    restraintfile = the name of the restraint file (e.g. 'dihre.itp')

    selection = the selection of the beta-peptide

    fc = force constant (defaults to 10000)
    """
    fc = float(fc)

    # first get the residue numbers
    space = {'residues': set()}
    cmd.iterate(selection, "residues.add(resv)", space=space)

    def presentandunique(selection, resi, atomname):
        count = cmd.count_atoms('({}) and (resi {}) and (name {})'.format(selection, resi, atomname))
        if count < 1:
            raise ValueError('Atom with name {} does not exist in residue {} of selection "{}"'.format(
                atomname, resi, selection))
        elif count > 1:
            raise ValueError('Atom name {} is not unique (found {}) in residue {} of selection "{}"'.format(
                atomname, count, resi, selection))
        return True

    def getdihedral(selection, resids, atomnames):
        # first check if all atoms are present and unique:
        for an, resi in zip(atomnames, resids):
            presentandunique(selection, resi, an)
        return cmd.get_dihedral(*['({0}) and (resi {1}) and (name {2})'.format(
            selection, resi, an) for resi, an in zip(resids, atomnames)],
                                quiet=True)

    def getindex(selection, resi, atomname):
        # check if all atoms are present and unique
        presentandunique(selection, resi, atomname)
        space = {'indices': []}
        cmd.iterate('({}) and (resi {}) and (name {})'.format(selection, resi, atomname),
                    'indices.append(index)', space=space)
        assert len(space['indices']) == 1
        return space['indices'][0]

    def dihedralline(selection, resids, atomnames, label=''):
        try:
            dihedral = getdihedral(selection, resids, atomnames)
            indices = [getindex(selection, resi, an) for resi, an in zip(resids, atomnames)]
        except ValueError as ve:
            return '; no dihedral restraint {}: {}\n'.format(label, ve.args[0])
        return '{:>8d} {:>8d} {:>8d} {:>8d}     1     {:8.2f}     0 {:8.2f}; {} ({}\t{}\t{}\t{})\n'.format(
            indices[0],
            indices[1],
            indices[2],
            indices[3],
            dihedral, fc, label,
            *['{}({})'.format(an, resi) for an, resi in zip(atomnames, resids)])

    with open(restraintfile, 'wt') as f:
        f.write('[ dihedral_restraints ]\n')
        for resi in sorted(space['residues']):
            # check each residue if they are alpha- or beta-ones
            thisresidue = '({}) and (resi {})'.format(selection, resi)
            space['atomnames'] = []
            cmd.iterate(thisresidue, "atomnames.append(name)", space=space)
            # in order for a residue to be recognized as an amino acid, it has to have the following atoms:
            # 'N', 'CA', 'C', ('O' or 'O1' and 'O2')
            ans = space['atomnames']  # abbreviation
            if not (('C' in ans)
                    and ('CA' in ans)
                    and ('N' in ans)
                    and (('O' in ans)
                         or (('O1' in ans)
                             and ('O2' in ans)))):
                # this is not an amino acid
                continue
            # see if 'N' is directly bonded to 'CA'
            if cmd.count_atoms(
                    '(neigh (({0}) and (resi {1}) and (name N))) and (({0}) and (resi {1}) and (name CA))'.format(
                        selection, resi)) == 1:
                # directly bonded, this is an alpha-amino acid
                f.write(dihedralline(selection, [resi - 1, resi, resi, resi], ['C', 'N', 'CA', 'C'],
                                     'phi{}'.format(resi)))  # phi
                f.write(dihedralline(selection, [resi, resi, resi, resi + 1], ['N', 'CA', 'C', 'N'],
                                     'psi{}'.format(resi)))  # psi
                f.write(dihedralline(selection, [resi - 1, resi - 1, resi, resi], ['O', 'C', 'N', 'HN+H'],
                                     'omega{}'.format(resi)))  # omega
            elif cmd.count_atoms(
                    '(neigh (({0}) and (resi {1}) and (name N))) and (neigh (({0}) and (resi {1}) and (name CA))) and (({0}) and (resi {1}) and name CB+CB1)'.format(
                        selection, resi)) == 1:
                # bonded through an atom which is called either 'CB' or 'CB1': this is a beta-amino acid
                f.write(dihedralline(selection, [resi - 1, resi, resi, resi], ['C', 'N', 'CB+CB1', 'CA'],
                                     'phi{}'.format(resi)))  # phi
                f.write(dihedralline(selection, [resi, resi, resi, resi], ['N', 'CB+CB1', 'CA', 'C'],
                                     'theta{}'.format(resi)))  # theta
                f.write(dihedralline(selection, [resi, resi, resi, resi + 1], ['CB+CB1', 'CA', 'C', 'N'],
                                     'psi{}'.format(resi)))  # psi
                f.write(dihedralline(selection, [resi - 1, resi - 1, resi, resi], ['O', 'C', 'N', 'HN+H'],
                                     'omega{}'.format(resi)))  # omega
            else:
                continue
