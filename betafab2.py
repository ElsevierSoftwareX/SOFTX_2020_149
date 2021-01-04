"""Fabricate a peptide from alpha- or beta-amino acid residues."""
import warnings

try:
    from pymol import cmd
    import chempy.models
    import pymol.plugins
except ImportError:
    warnings.warn(
        'Cannot import PyMOL: functionality will suffer (you can ignore this if you are just building the documentation).')
    cmd = None
import re
import os
from .utils import tempObjectName, fetchModel, set_dihedral
from typing import Optional, List, Dict, Union
import logging
from .secstructdb import SecondaryStructureDB
import itertools
import random

logger = logging.getLogger()
logger.setLevel(logging.INFO)


class BetaPeptideException(Exception):
    pass


class BetaPeptideConsistencyError(BetaPeptideException):
    pass


class BetaPeptide:
    """A class for constructing alpha/beta peptides.
    """
    GREEKLETTERS = 'ABGDEZHTIKLMNP'
    SIDECHAINS = {
        'A': 'ala',
        'C': 'cys',
        'CM': 'cys',  # this will be treated differently in the code, since PyMOL does not have a ready-made fragment
        'D': 'asp',
        'DH': 'asph',
        'E': 'glu',
        'EH': 'gluh',
        'F': 'phe',
        'G': 'gly',
        'HD': 'hid',
        'HE': 'hie',
        'HH': 'hip',
        'H': 'his',
        'I': 'ile',
        'K': 'lys',
        'O': 'orp',  # this will be treated differently in the code, since PyMOL does not have a ready-made fragment
        'KN': 'lysn',
        'ON': 'orn',  # this will be treated differently in the code, since PyMOL does not have a ready-made fragment
        'L': 'leu',
        'M': 'met',
        'N': 'asn',
        'P': 'pro',  # proline, only available for alpha-peptides!
        'Q': 'gln',
        'R': 'arg',
        'S': 'ser',
        'T': 'thr',
        'V': 'val',
        'W': 'trp',
        'Y': 'tyr',
    }

    def __init__(self, modelname: Optional[str], isbeta: bool = True, alphasidechain: str = 'G', alphastereo: str = 'S',
                 betasidechain: str = 'G', betastereo: str = 'S'):
        """Create a single amino-acid. Peptides can be built by *adding* them

        :param modelname: the name under which this peptide will be known by PyMOL
        :type modelname: str
        :param isbeta: True if this is to be a beta-amino acid, False if an alpha
        :type isbeta: bool
        :param alphasidechain: the sidechain on the C-alpha
        :type alphasidechain: str
        :param alphastereo: chirality of the C-alpha
        :type alphastereo: str, "R" or "S" (or "D" or "L" if `isbeta` is False, i.e. alpha-amino acids)
        :param betasidechain: the sidechain on the C-beta
        :type betasidechain: str
        :param betastereo: chirality of the C-beta
        :type betastereo: str, "R" or "S"
        """
        self._modelname = modelname
        if modelname is None:
            return
        if isbeta:
            self.initializeBetaResidue(alphasidechain, alphastereo, betasidechain, betastereo)
        else:
            self.initializeAlphaResidue(alphasidechain, alphastereo)

    def copy(self, newname: str) -> "BetaPeptide":
        """Make a deep copy of ourselves
        
        :param newname: the new object name
        :type newname: str
        :return: the new BetaPeptide object
        :rtype: BetaPeptide
        """
        model = cmd.get_model('model {}'.format(self._modelname))
        cmd.delete('model {}'.format(newname))
        cmd.load_model(model, newname)
        obj = BetaPeptide(None)
        obj._modelname = newname
        return obj

    def initializeAlphaResidue(self, sidechain: str = 'G', stereo: str = 'S'):
        """Initialize this residue to be an alpha-amino acid

        :param sidechain: designation of the sidechain
        :type sidechain: str
        :param stereo: absolute chirality of the sidechain ('R' or 'S')
        :type stereo: str
        :return: None
        :rtype: NoneType
        """
        cmd.delete('model {}'.format(self._modelname))
        try:
            cmd.fragment(self.SIDECHAINS[sidechain], self._modelname)
        except pymol.CmdException as exc:
            # could not load fragment: this might be a special one. Try to load it from our resource folder
            resourcedir = os.path.join(os.path.split(pymol.plugins.plugins['pmlbeta'].filename)[0], 'resource')
            cmd.load(os.path.join(resourcedir, '{}.pkl'.format(self.SIDECHAINS[sidechain])), self._modelname)
        if sidechain == 'CM':
            cmd.remove('model {} and name HG'.format(self._modelname))
            cmd.alter('model {} and name SG'.format(self._modelname), 'formal_charge=-1')
            cmd.alter('model {} and name SG'.format(self._modelname), 'partial_charge=-1')
            cmd.alter('model {}'.format(self._modelname), 'resn="CYM"')
        cmd.alter('model {}'.format(self._modelname), 'resv=1')
        cmd.alter('model {}'.format(self._modelname), 'chain="A"')

        # alpha-residues are by default all have L chirality. Now let's fix the chirality. The case of cisteine is
        # complicated...

        if sidechain != 'G' and (  # skip the achiral glycine
                (stereo == 'D') or  # all amino-acids if D is requested
                (stereo == 'R' and sidechain not in ['C',
                                                     'CM']) or  # for most amino-acids, L is S, thus for R, we have to flip
                (stereo == 'S' and sidechain in ['C', 'CM'])  # for Cys, L is R, thus for S, we have to flip
        ):
            cmd.edit('model {} and resi {} and name CA'.format(self._modelname, 1),
                     'model {} and resi {} and name HA'.format(self._modelname, 1),
                     'model {} and resi {} and name C'.format(self._modelname, 1))
            cmd.invert()
            cmd.unpick()

        # now rename the atoms
        with fetchModel(self._modelname) as model:
            for a in model.atom:
                if a.name in ['N', 'H', 'C', 'O']:
                    a.name = a.name + '_sc0'
                else:
                    a.name = a.name + '_sc2'
        self.renameAtomsInResidue(1)

    def checkPeptide(self):
        """Check this peptide for validity

        The following tests are done:
            1. each non-N-terminal peptide nitrogen must have one neighbour named "H" or "HN", except alpha-proline
            2. each peptide nitrogen must have only one neighbour named "CB" or "CB1"
            3. each peptide carbon must have only one neighbour named "CA"
            4. each non-C-terminal peptide carbon must have exactly one neighbour named "N" with resi=resi+1
            5. each non-C-terminal peptide carbon must have exactly one neighbour named "O"
            6. each non-N-terminal peptide nitrogen must have exactly one neighbour named "C" with resi=resi-1
            7. residue numbering must be continuous and start from one.

        Capping groups are exempt from the above rules:
            . ACE from rule 3.
            . NME from rule 2.

        :raises BetaPeptideConsistencyError: when one of the above tests fails
        """
        with fetchModel(self._modelname) as model:
            minresidue = min([a.resi_number for a in model.atom])
            maxresidue = max([a.resi_number for a in model.atom])

        if minresidue != 1:
            raise BetaPeptideConsistencyError(
                'Lowest residue number is {} instead of the expected 1.'.format(minresidue))

        for i in range(1, maxresidue + 1):
            resname = cmd.get_model('model {} and resi {}'.format(self._modelname, i)).atom[0].resn
            logger.debug('Residue {}: {}'.format(i, resname))
            if (i > 1) and resname not in ['PRO', 'DPRO']:  # not the N-terminal and not proline
                # count amide hydrogens on the nitrogen: must be one on a non-N-terminal residue
                hcount = cmd.count_atoms('model {} and resi {} and (neighbor name N ) and (name H+HN)'.format(
                    self._modelname, i))
                if hcount != 1:
                    raise BetaPeptideConsistencyError('Amide N in residue {} has {} hydrogens.'.format(i, hcount))
                # count amide carbons bonded to the nitrogen: must be one on a non-N-terminal residue. This also ensures
                # continuity of residue numbering (partly)
                ccount = cmd.count_atoms('(neighbor (model {} and resi {} and name N)) and (resi {} and name C)'.format(
                    self._modelname, i, i - 1))
                if ccount != 1:
                    raise BetaPeptideConsistencyError(
                        'Amide N in residue {} has {} amide carbons from the previous residue.'.format(i, ccount))
            if i < maxresidue:
                # not the C-terminal 
                # count amide oxygens on the carbon: must be one on a non-C-terminal residue
                ocount = cmd.count_atoms('model {} and resi {} and (neighbor name C) and name O'.format(
                    self._modelname, i))
                if ocount != 1:
                    raise BetaPeptideConsistencyError('Amide C in residue {} has {} oxygens.'.format(i, ocount))
                # count amide nitrogens bonded to the carbon: must be one on a non-C-terminal residue. Together with 
                # the previous similar check, this ensures the continuity of residue numbering.
                ncount = cmd.count_atoms(
                    '(neighbor (model {} and resi {} and name C)) and (resi {} and name N)'.format(
                        self._modelname, i, i + 1))
                if ncount != 1:
                    raise BetaPeptideConsistencyError(
                        'Amide C in residue {} has {} amide nitrogens from the next residue'.format(i, ncount)
                    )
            # all residues, including terminals
            cbcount = cmd.count_atoms('model {} and resi {} and (neighbor name N) and (name CB+CB1+CA)'.format(
                self._modelname, i))
            if (cbcount != 1) and (resname != 'NME') and (resname != 'ACE') and (resname != 'BUT'):
                raise BetaPeptideConsistencyError('Amide N in residue {} has {} CA or CB neighbours.'.format(
                    i, cbcount))
            cacount = cmd.count_atoms('model {} and resi {} and (neighbor name C) and (name CA)'.format(
                self._modelname, i))
            if (cacount != 1) and (resname != 'ACE') and (resname != 'NME'):
                raise BetaPeptideConsistencyError('Amide C in residue {} has {} CA neighbours.'.format(i, cacount))

    @property
    def maxResidue(self) -> int:
        """Get the highest residue number
        
        :return: the highest residue number
        :rtype: int
        """
        with fetchModel(self._modelname) as model:
            return max([a.resi_number for a in model.atom])

    def isNProtected(self) -> bool:
        """Check if the N-terminal is protected
        
        :return: True if the N-terminal is protected 
        :rtype: bool
        """
        return cmd.count_atoms(
            '(neighbor (model {} and name N and resi 1)) and not name CB+CB1+CA+H+HN+CH3'.format(self._modelname)) > 0

    def isCProtected(self) -> bool:
        """Check if the C-terminal is protected
        
        :return: True if the C-terminal is protected
        :rtype: bool
        """
        maxresidue = self.maxResidue
        return cmd.count_atoms(
            '(neighbor (model {} and name C and resi {})) and (not name O+CA+CH3)'.format(self._modelname,
                                                                                          maxresidue)) > 0

    def startsWithProline(self) -> bool:
        """Check if the N-terminus is a proline residue

        :return: True if the N-terminus is an alpha-proline
        :rtype: bool
        """
        with fetchModel(self._modelname) as m:
            minresi = min({a.resi_number for a in m.atom})
            firstresidue = {a.resn for a in m.atom if a.resi_number == minresi}
            if len(firstresidue) > 1:
                raise BetaPeptideConsistencyError('Atoms in the same residue must have the same residue names.')
            else:
                assert len(firstresidue) == 1  # len()==0 cannot happen.
                return firstresidue.pop() in ['PRO', 'DPRO']

    def __iadd__(self, right: "BetaPeptide") -> "BetaPeptide":
        """Append another peptide to our C-terminal

        :param right: the C-terminal part
        :type right: BetaPeptide
        :return: the concatenated beta-peptide (the present one is overwritten).
        :rtype: BetaPeptide
        """
        left = self  # the code below will be more logical this way.
        if not isinstance(left, BetaPeptide) or not isinstance(right, BetaPeptide):
            raise TypeError('Only peptides can be concatenated!')
        if left.isCProtected():
            raise ValueError('The N-terminal part is C-protected!')
        if right.isNProtected() and (not right.startsWithProline()):
            raise ValueError('The C-terminal part is N-protected!')
        # do some validation
        left.checkPeptide()
        right.checkPeptide()
        with tempObjectName() as mouldname, tempObjectName() as lefttemp, tempObjectName() as righttemp:
            self.loadPeptideBond(mouldname)
            mould = BetaPeptide(None)
            mould._modelname = mouldname
            left = left.copy(lefttemp)
            right = right.copy(righttemp)
            # increment the residue numbers in right
            maxleft = left.maxResidue
            with fetchModel(righttemp) as model:
                for a in model.atom:
                    a.resi_number += maxleft
            # mould is a peptide bond, with the following layout:
            #
            #   O        Cnext
            #   \\      /
            #    C --- N
            #   /       \
            #  Cprev     H
            #
            # align the following:
            #   C of the last residue of "left" to C of mould
            #   O of the last residue of "left" to O of mould
            #   CA or CH3 of the last residue of "left" to Cprev of mould

            # in case of alpha-proline:
            #
            #   O           CA
            #   \\        /    \
            #    C --- N        CB
            #   /       \      /
            #  Cprev     CD - CG
            logger.debug('Fitting left side')
            logger.debug('Models: \n' + '\n'.join(cmd.get_object_list()))
            logger.debug('Mould name: {}'.format(mouldname))
            logger.debug('Lefttemp: {}'.format(lefttemp))
            logger.debug('Righttemp: {}'.format(righttemp))
            self.safe_pairfit('model {}'.format(lefttemp),
                              'model {}'.format(mouldname),
                              maxleft,
                              1,
                              ['CA+CH3', 'C', 'O'],
                              ['Cprev', 'C', 'O'],
                              0.1)
            # find the neighbour of N
            cbname = cmd.get_model(
                '(neighbor (model {} and resi {} and name N)) and (symbol C) and (not name C)'.format(righttemp,
                                                                                                      1 + maxleft)).atom[
                0].name

            #   N of the first residue of "right" to N of mould
            #   H of the first residue of "right" to H of mould
            #   CB, CB1 or CA or CH3 of the first residue of "right" to Cnext of mould
            logger.debug('Fitting right side')
            self.safe_pairfit('model {}'.format(righttemp),
                              'model {}'.format(mouldname),
                              1 + maxleft,
                              1,
                              ['N', 'H' if not right.startsWithProline() else 'CD', cbname],
                              ['N', 'H', 'Cnext'],
                              0.1 if not right.startsWithProline() else 0.3)
            cmd.fuse("model {} and resi {} and name C".format(lefttemp, maxleft),
                     "model {} and resi {} and name N".format(righttemp, 1 + maxleft),
                     mode=3)
            cmd.bond("model {} and resi {} and name C".format(righttemp, maxleft),
                     "model {} and resi {} and name N".format(righttemp, 1 + maxleft), quiet=True)
            obj = right.copy(self._modelname)
        return obj

    @staticmethod
    def safe_pairfit(model1: str, model2: str, resi1: int, resi2: int, atoms1: List[str], atoms2: List[str],
                     rmstolerance: float = 0.1, ntries: int = 10):
        """Align model1 and model2 by the given atoms. Try to find an optimal tolerance.

        The problem with the original pair_fit of PyMOL was that sometimes it found a local minimum, not
        the desired one. If it converged into a local minimum and moved out by a random rotation, it
        can find the global minimum.

        :param model1: the name of the model to be moved
        :type model1: str
        :param model2: the name of the stationary model
        :type model2: str
        :param resi1: residue in model 1
        :type resi1: int
        :param resi2: residue in model 2
        :type resi2: int
        :param atoms1: name of atoms in the movable model
        :type atoms1: list of str
        :param atoms2: name of atoms in the stationary model
        :type atoms2: list of str
        :param rmstolerance: try to fit the models until the RMS is less than this
        :type rmstolerance: float
        :param ntries: the number of trials
        :type ntries: int
        :return: the final rms
        :rtype: float
        """
        selections1 = ['({}) and resi {} and name {}'.format(model1, resi1, name) for name in atoms1]
        selections2 = ['({}) and resi {} and name {}'.format(model2, resi2, name) for name in atoms2]
        selections = list(itertools.chain.from_iterable(zip(selections1, selections2)))
        for i in range(ntries):
            rms = cmd.pair_fit(*selections, quiet=True)
            if rms < rmstolerance:
                return rms
            else:
                # it may not have converged: do a random rotation and try to fit
                cmd.rotate('xyz'[random.randint(0, 2)], random.randrange(0, 360), model1)
        else:
            raise RuntimeError('Could not fit with RMS less than {}. RMS is: {}'.format(rmstolerance, rms))

    @staticmethod
    def betaResidueName(sidechain2: Optional[str] = None, sidechain3: Optional[str] = None) -> str:
        """Get a residue name usable in GROMACS/CHARMM input files for a beta-amino acid

        :param sidechain2: one (or two-)-letter abbreviation of the alpha-sidechain
        :type sidechain2: str
        :param sidechain3: one (or two-)-letter abbreviation of the beta-sidechain
        :type sidechain3: str
        :return: the residue name
        :rtype: str
        """
        if sidechain2 in [None, 'G'] and sidechain3 in [None, 'G']:
            return 'BALA'
        elif sidechain2 in [None, 'G']:
            return 'B3' + sidechain3
        elif sidechain3 in [None, 'G']:
            return 'B2' + sidechain2
        else:
            return 'B0' + sidechain2 + sidechain3

    def initializeBetaResidue(self, sidechain2: str = 'G', stereo2: str = 'S', sidechain3: str = 'G',
                              stereo3: str = 'S'):
        """Initialize this to be a beta-amino acid residue

        :param sidechain2: designation of the side chain of the alpha carbon
        :type sidechain2: str
        :param stereo2: absolute conformation of the alpha side chain
        :type stereo2: str, either 'R' or 'S'
        :param sidechain3: designation of the side chain of the beta carbon
        :type sidechain3: str
        :param stereo3: absolute conformation of the beta side chain
        :type stereo3: str, either 'R' or 'S'
        :return: None
        :rtype: NoneType
        :raises ValueError: on invalid input
        """
        residuename = self.betaResidueName(sidechain2, sidechain3)
        if sidechain2 == 'P' or sidechain3 == 'P':
            raise NotImplementedError('Beta-Proline residues not yet supported')
        # delete the model
        cmd.delete('model {}'.format(self._modelname))

        # prime us with a homo-glycine
        self.loadBetaBackbone(self._modelname)
        # set the residue numbers to 1
        with fetchModel(self._modelname) as model:
            for a in model.atom:
                a.resi_number = 1
                a.resn = residuename
                a.chain = 'A'
                a.name = a.name + '_sc0'  # add a tag identifying this atom to be part of the backbone

        # attach alpha sidechain if needed
        if sidechain2 != 'G':
            self.attachBetaSideChain(1, 2, stereo2, sidechain2)
            # rename the atoms
        # attach beta sidechain if needed
        if sidechain3 != 'G':
            self.attachBetaSideChain(1, 3, stereo3, sidechain3)
            # rename the atoms
        self.renameAtomsInResidue(1)

        cmd.alter('model {} and resi {}'.format(self._modelname, 1), 'chain="A"')
        # Adjustment of the chirality is needed in the following special cases:
        if sidechain2 in ['C', 'CM', 'S', 'T']:
            # switch chirality of the alpha carbon
            cmd.edit('model {} and resi {} and name CA'.format(self._modelname, 1),
                     'model {} and resi {} and name CB+CB1'.format(self._modelname, 1),
                     'model {} and resi {} and name C'.format(self._modelname, 1))
            cmd.invert()
            cmd.unpick()
        if sidechain3 in ['C', 'CM', 'S', 'T']:
            # switch chirality of the beta carbon
            cmd.edit('model {} and resi {} and name CB+CB1'.format(self._modelname, 1),
                     'model {} and resi {} and name CA'.format(self._modelname, 1),
                     'model {} and resi {} and name N'.format(self._modelname, 1))
            cmd.invert()
            cmd.unpick()
        if sidechain3 in ['M', 'D', 'DH'] and sidechain2 in ['G']:
            # switch chirality of the beta carbon
            cmd.edit('model {} and resi {} and name CB+CB1'.format(self._modelname, 1),
                     'model {} and resi {} and name CA'.format(self._modelname, 1),
                     'model {} and resi {} and name N'.format(self._modelname, 1))
            cmd.invert()
            cmd.unpick()

    def attachBetaSideChain(self, residue: int, site: int, stereo: str, sidechain: str):
        """Attach a side-chain to the beta-backbone. The atoms in the side-chain (including the attach site and its other
        hydrogen) will be correctly named.

        :param residue: number of the residue to be modified
        :type residue: int
        :param site: the site to attach the sidechain: either 2 or 3
        :type site: int
        :param stereo: absolute conformation
        :type stereo: str, either 'R' or 'S'
        :param sidechain: sidechain designation
        :type sidechain: str
        :return: None
        :rtype: NoneType
        :raises ValueError: if `site` or `stereo` is incorrect
        :raises KeyError: if `sidechain` is unknown
        """
        with tempObjectName() as fragmentname:
            # load the appropriate amino-acid fragment
            try:
                cmd.fragment(self.SIDECHAINS[sidechain], fragmentname)
            except pymol.CmdException as exc:
                # could not load fragment: this might be a special one. Try to load it from our resource folder
                resourcedir = os.path.join(os.path.split(pymol.plugins.plugins['pmlbeta'].filename)[0], 'resource')
                cmd.load(os.path.join(resourcedir, '{}.pkl'.format(self.SIDECHAINS[sidechain])), fragmentname)
            # keep only the side-chain and Calpha
            cmd.remove('model {} and (name C+O+N+H+HA)'.format(fragmentname))
            # hack for residue 'CM'
            if sidechain == 'CM':
                cmd.remove('model {} and name HG'.format(fragmentname))
                cmd.alter('model {} and name SG'.format(fragmentname), 'formal_charge=-1')

            has_beta2sidechain = cmd.count_atoms(
                'model {} and resi {} and name CB1'.format(self._modelname, residue)) == 1
            # find the name of the backbone atoms where the sidechain will be attached
            if site == 2 and stereo == 'S':
                attachatoms = ('CA_sc0', 'HA1_sc0')
            elif site == 2 and stereo == 'R':
                attachatoms = ('CA_sc0', 'HA2_sc0')
            elif site == 3 and stereo == 'S':
                attachatoms = (
                    'CB1_sc0' if has_beta2sidechain else 'CB_sc0', 'HB11_sc0' if has_beta2sidechain else 'HB1_sc0')
            elif site == 3 and stereo == 'R':
                attachatoms = (
                    'CB1_sc0' if has_beta2sidechain else 'CB_sc0', 'HB12_sc0' if has_beta2sidechain else 'HB2_sc0')
            else:
                raise ValueError('Invalid site/stereo combination: {}/{}'.format(site, stereo))

            # relabel the backbone atoms to avoid name clashes: simply add a "_sc0" tag at the end.
            with fetchModel(self._modelname) as model:
                for a in model.atom:
                    assert isinstance(a, chempy.Atom)
                    if a.resi_number != residue:
                        continue
                    if '_' not in a.name:  # don't add two tags if one is already present.
                        a.name = a.name + '_sc0'
            # relabel the sidechain atoms to avoid name clashes: add a "_sc<i>" tag at the end, where <i> is 2 or 3
            with fetchModel(self._modelname) as model:
                residuename = [a for a in model.atom if a.resi_number == residue][0].resn

            with fetchModel(fragmentname) as model:
                for a in model.atom:
                    assert isinstance(a, chempy.Atom)
                    a.name = a.name + '_sc{}'.format(site)
                    a.resi_number = residue
                    a.resn = residuename
            # select the first heavy atom (i.e. not hydrogen) attached to the CA in the backbone-stripped sidechain
            # model (typically CB)
            with tempObjectName(nocreate=True) as tempsel:
                cmd.select(tempsel,
                           'neighbor (model {0} and name CA_sc{1}) and not symbol H'.format(fragmentname, site))
                caneighbour = cmd.get_model('%{}'.format(tempsel)).atom[0].name
            # align the two models: fit CA and CB (or the first non-hydrogen attached to CA) in the backbone-stripped
            # sidechain to the carbon and hydrogen on the beta-backbone chosen above for attaching.
            cmd.pair_fit('(model {0} and name CA_sc{1}+{2})'.format(fragmentname, site, caneighbour),
                         'model {} and resi {} and name {}+{}'.format(self._modelname, residue, *attachatoms),
                         quiet=True)
            cmd.remove('model {} and name CA_sc{}'.format(fragmentname, site))
            cmd.remove('model {} and resi {} and name {}'.format(self._modelname, residue, attachatoms[1]))
            cmd.fuse('model {} and name {}'.format(fragmentname, caneighbour),
                     'model {} and resi {} and name {}'.format(self._modelname, residue, attachatoms[0]))
        # we are done with fusing, only renaming is needed.

    def renameAtomsInResidue(self, residue: int) -> None:
        """Rename atoms in a residue

        :param int residue: the residue number.
        """
        with fetchModel(self._modelname) as model:
            # first deconstruct the names of the atoms in the sidechains and the backbone into the form:
            #   <symbol>:<greek>:<index>:<sidechain>
            #
            # for example:
            #    CD1 in the alpha-sidechain -> C:D:1:2
            #    NE2 in the beta-sidechain -> N:E:2:3
            # hydrogens are skipped for the time being.
            for a in model.atom:
                assert isinstance(a, chempy.Atom)
                if a.resi_number != residue:
                    # skip atoms which do not belong to this residue
                    continue
                if a.symbol == 'H':
                    # this is a hydrogen, skip it for now.
                    continue
                m = re.match(
                    '^(?P<symbol>[BCNOFPS])(?P<greek>[{}])(?P<index>[123456789])?_sc(?P<sidechain>[023])$'.format(
                        self.GREEKLETTERS), a.name)
                if m is None:
                    if a.name in ['C_sc0', 'O_sc0', 'N_sc0', 'H_sc0']:
                        # these are legal
                        continue
                    raise ValueError('Cannot match atom name: {}'.format(a.name))
                m = m.groupdict()
                if m['index'] is None:
                    m['index'] = 0
                if m['sidechain'] == '3':
                    # increment the greek letters by one
                    m['greek'] = self.GREEKLETTERS[self.GREEKLETTERS.index(m['greek']) + 1]
                a.name = '{0[symbol]}_{0[greek]}_{0[index]}_{0[sidechain]}'.format(m)
            # now fix the numbering of heavy atoms. Sorting categories are (from strongest to weakest):
            #    - backbone, sidechain3, sidechain2
            #    - original index

            # the names of the heavy atoms at this point are: <symbol>_<greek>_<index>_<sidechain>, where:
            #    - symbol is the element (e.g. C, N, O, S...)
            #    - greek is the greek letter (e.g A, B, G, D, E, Z, H, T ...)
            #    - index is the number of this atom in the original residue (e.g. the "2" from CD2)
            #    - sidechain is the sidechain designation, '0' being the backbone, '2' the alpha- and '3' the beta-sidechain
            scsort = {'0': 0, '3': 1, '2': 2}
            for g in self.GREEKLETTERS:
                atoms = [a for a in model.atom if
                         a.symbol != 'H' and '_' in a.name and a.resi_number == residue and a.name.split('_')[1] == g]
                if not atoms:
                    continue
                elif len(atoms) == 1:
                    # this is the simplest case, no indexing required
                    a = atoms[0]
                    symbol, greek, index, sidechain = a.name.split('_')
                    a.name = symbol + greek
                else:
                    logger.debug('Sorting heavy atoms for greek letter {}:'.format(g))
                    for i, a in enumerate(
                            sorted(atoms, key=lambda a: (scsort[a.name.split('_')[3]], a.name.split('_')[2]))):
                        logger.debug('{}: {}'.format(i + 1, a.name))
                        symbol, greek, index, sidechain = a.name.split('_')
                        a.name = symbol + greek + str(i + 1)
        # re-fetch it from PyMOL, because we will use cmd.select for finding hydrogen neighbours of heavy atoms
        with fetchModel(self._modelname) as model:
            # now fix the hydrogens
            for a in model.atom:
                assert isinstance(a, chempy.Atom)
                if a.symbol == 'H':
                    continue
                if a.resi_number != residue:
                    continue
                # find the atom to which this is bonded
                atomindex = model.atom.index(a)
                pairindices = [b.index for b in model.bond if atomindex in b.index]
                neighbourindices = [x[0] if x[1] == atomindex else x[1] for x in pairindices]
                hydrogens = [model.atom[idx] for idx in neighbourindices if model.atom[idx].symbol == 'H']
                if not hydrogens:
                    continue
                elif len(hydrogens) == 1:
                    hydrogens[0].name = 'H' + a.name[1:]
                else:
                    # hydrogens are numbered as 1HG1, 2HG1, 3HG1 originally, i.e. the first character is the index.
                    # Another nomenclature is HG11, HG12, HG13. Either way, an alphabetical sort should do the trick of
                    # preserving the original order.
                    for i, h in enumerate(sorted(hydrogens, key=lambda h: h.name)):
                        h.name = 'H' + a.name[1:] + str(i + 1)

        with fetchModel(self._modelname) as model:
            # see if some atoms are still having some :tags at the end of their names. Must be the N:sc0, H:sc0, O:sc0, C:sc0
            for a in model.atom:
                if a.resi_number != residue:
                    continue
                if a.name.endswith('_sc0'):
                    assert a.name.split('_')[0] in ['N', 'H', 'O', 'C']
                    a.name = a.name.split('_')[0]
                assert '_' not in a.name

    def isBetaResidue(self, residue: int) -> bool:
        """Check if the residue is a beta-amino acid

        :param residue: residue index
        :type residue: int
        :return: True if it is a beta-amino acid
        :rtype: bool
        """
        if cmd.count_atoms('model {0} and resi {1} and name N+C+H+CA+O'.format(self._modelname, residue)) != 5:
            # this might not even be an amino-acid!
            return False
        if cmd.count_atoms(
                '(neighbor (model {0} and resi {1} and name N)) and resi {1} and name H'.format(self._modelname,
                                                                                                residue)) != 1:
            # N-H bond missing
            return False
        if cmd.count_atoms(
                '(neighbor (model {0} and resi {1} and name C)) and resi {1} and name CA'.format(self._modelname,
                                                                                                 residue)) != 1:
            # N-CA bond missing
            return False
        if cmd.count_atoms(
                '(neighbor (model {0} and resi {1} and name C)) and (neighbor (model {0} and resi {1} and name N)) and (resi {1} and name CA)'.format(
                    self._modelname, residue)) == 1:
            # N and C have a common neighbour, named CA. This is an alpha-amino acid
            return False
        if cmd.count_atoms(
                '(neighbor (model {0} and resi {1} and name CA)) and (neighbor (model {0} and resi {1} and name N)) and (resi {1} and name CB+CB1)'.format(
                    self._modelname, residue)) == 1:
            # N and CA have a common neighbour, named either CB or CB1. This is a beta-amino acid.
            return True
        return False

    def isAlphaResidue(self, residue: int) -> bool:
        """Check if the residue is an alpha-amino acid

        :param residue: residue index
        :type residue: int
        :return: True if it is an alpha-amino acid
        :rtype: bool
        """
        if cmd.count_atoms('model {0} and resi {1} and name N+C+H+CA+O'.format(self._modelname, residue)) != 5:
            # this might not even be an amino-acid!
            return False
        if cmd.count_atoms(
                '(neighbor (model {0} and resi {1} and name N)) and resi {1} and name H'.format(self._modelname,
                                                                                                residue)) != 1:
            # N-H bond missing
            return False
        if cmd.count_atoms(
                '(neighbor (model {0} and resi {1} and name C)) and resi {1} and name CA'.format(self._modelname,
                                                                                                 residue)) != 1:
            # N-CA bond missing
            return False
        if cmd.count_atoms(
                '(neighbor (model {0} and resi {1} and name C)) and (neighbor (model {0} and resi {1} and name N)) and (resi {1} and name CA)'.format(
                    self._modelname, residue)) == 1:
            # N and C have a common neighbour, named CA. This is an alpha-amino acid
            return True
        if cmd.count_atoms(
                '(neighbor (model {0} and resi {1} and name CA)) and (neighbor (model {0} and resi {1} and name N)) and (resi {1} and name CB+CB1)'.format(
                    self._modelname, residue)) == 1:
            # N and CA have a common neighbour, named either CB or CB1. This is a beta-amino acid.
            return False
        return False

    def fold(self, residue: int, dihedrals: Union[List[float], str]):
        """Fold a residue into the desired secondary structure

        The secondary structure can be given in `dihedrals` in two ways:
            1. a list of the torsion angles in degrees (2 or 3, depending on the type
                of the amino-acid
            2. the name of a secondary structure in the BETAFAB_HELIXTYPES PyMOL Plugin
                variable

        :param residue: the number of the residue
        :type residue: int
        :param dihedrals: secondary structure designation
        :type dihedrals: list of floats or a str
        :raises ValueError: if the residue is neither an alpha-, nor a beta-amino acid
        """
        if isinstance(dihedrals, str):
            dihedrals = SecondaryStructureDB.dihedrals(dihedrals)
        logger.debug('Folding alpha-residue {} to {}'.format(residue, dihedrals))
        if self.isAlphaResidue(residue):
            logger.debug('Folding alpha-residue {} to {}'.format(residue, dihedrals))
            logger.debug(
                set_dihedral(self._modelname, ('C', residue - 1), ('N', residue), ('CA', residue), ('C', residue),
                             dihedrals[0]))
            logger.debug(
                set_dihedral(self._modelname, ('N', residue), ('CA', residue), ('C', residue), ('N', residue + 1),
                             dihedrals[1]))
        elif self.isBetaResidue(residue):
            logger.debug('Folding beta-residue {} to {}'.format(residue, dihedrals))
            logger.debug(
                set_dihedral(self._modelname, ('C', residue - 1), ('N', residue), ('CB+CB1', residue), ('CA', residue),
                             dihedrals[0]))
            logger.debug(
                set_dihedral(self._modelname, ('N', residue), ('CB+CB1', residue), ('CA', residue), ('C', residue),
                             dihedrals[1]))
            logger.debug(
                set_dihedral(self._modelname, ('CB+CB1', residue), ('CA', residue), ('C', residue), ('N', residue + 1),
                             dihedrals[2]))
        else:
            raise ValueError('Residue {} is neither an alpha-, nor a beta-amino acid'.format(residue))

    @staticmethod
    def loadBetaBackbone(objectname: str) -> str:
        """Create a new bare beta-backbone as a new object

        :param objectname: desired name of the new object
        :type objectname: str
        :return: the name of the new object
        :rtype: str
        """
        resourcedir = os.path.join(os.path.split(pymol.plugins.plugins['pmlbeta'].filename)[0], 'resource')
        cmd.load(os.path.join(resourcedir, 'betabackbone.pkl'), objectname)
        return objectname

    @staticmethod
    def loadPeptideBond(objectname: str) -> str:
        """Create a new object: a bare peptide bond
        :param objectname: desired name of the new object
        :type objectname: str
        :return: the name of the new object
        :rtype: str
        """
        resourcedir = os.path.join(os.path.split(pymol.plugins.plugins['pmlbeta'].filename)[0], 'resource')
        cmd.load(os.path.join(resourcedir, 'peptidebond.pkl'), objectname)
        return objectname

    def _debugModel(self, model=None):
        if model is None:
            logger.debug('Fetching model')
            with fetchModel(self._modelname) as model:
                logger.debug('Fetched model.')
                self._debugModel(model)
                logger.debug('End of recursive call.')
        else:
            logger.debug('Atoms in the model:')
            for a in model.atom:
                logger.debug('    - name: {}, resi: {}, resn: {}'.format(a.name, a.resi_number, a.resn))
            logger.debug('Bonds in the model:')
            for b in model.bond:
                logger.debug('    - {} - {} (order {})'.format(model.atom[b.index[0]].name, model.atom[b.index[1]].name,
                                                               b.order))

    @classmethod
    def But(cls, modelname=None):
        obj = cls(None)
        resourcedir = os.path.join(os.path.split(pymol.plugins.plugins['pmlbeta'].filename)[0], 'resource')
        cmd.load(os.path.join(resourcedir, 'butyrate.pkl'), modelname)
        cmd.alter('model {}'.format(modelname), 'resv=1')
        cmd.alter('model {}'.format(modelname), 'chain="A"')
        obj._modelname = modelname
        return obj

    @classmethod
    def Ace(cls, modelname=None):
        obj = cls(None)
        cmd.fragment('ace', modelname)
        for i in range(1, 4):
            cmd.alter('model {} and name {}HH3'.format(modelname, i), 'name="HH3{}"'.format(i))
        cmd.alter('model {}'.format(modelname), 'resv=1')
        cmd.alter('model {}'.format(modelname), 'chain="A"')
        obj._modelname = modelname
        return obj

    @classmethod
    def NMe(cls, modelname=None):
        obj = cls(None)
        cmd.fragment('nme', modelname)
        for i in range(1, 4):
            cmd.alter('model {} and name {}HH3'.format(modelname, i), 'name="HH3{}"'.format(i))
        cmd.alter('model {}'.format(modelname), 'resv=1')
        cmd.alter('model {}'.format(modelname), 'chain="A"')
        obj._modelname = modelname
        return obj

    @classmethod
    def ACHC(cls, stereo2: str, stereo3: str, modelname: str = None):
        obj = cls(None)
        resourcedir = os.path.join(os.path.split(pymol.plugins.plugins['pmlbeta'].filename)[0], 'resource')
        try:
            cmd.load(os.path.join(resourcedir, 'ACHC_2{}3{}.pkl'.format(stereo2.upper(), stereo3.upper())), modelname)
        except pymol.CmdException:
            raise ValueError('Invalid stereochemistry: {} and {}'.format(stereo2, stereo3))
        cmd.alter('model {}'.format(modelname), 'resv=1')
        cmd.alter('model {}'.format(modelname), 'chain="A"')
        cmd.alter('model {}'.format(modelname), 'resn="ACHC"')
        obj._modelname = modelname
        return obj

    @classmethod
    def ACPC(cls, stereo2: str, stereo3: str, modelname: str = None):
        obj = cls(None)
        resourcedir = os.path.join(os.path.split(pymol.plugins.plugins['pmlbeta'].filename)[0], 'resource')
        try:
            cmd.load(os.path.join(resourcedir, 'ACPC_2{}3{}.pkl'.format(stereo2.upper(), stereo3.upper())), modelname)
        except pymol.CmdException:
            raise ValueError('Invalid stereochemistry: {} and {}'.format(stereo2, stereo3))
        cmd.alter('model {}'.format(modelname), 'resv=1')
        cmd.alter('model {}'.format(modelname), 'chain="A"')
        cmd.alter('model {}'.format(modelname), 'resn="ACPC"')
        obj._modelname = modelname
        return obj

    @staticmethod
    def parseBetaPeptideSequence(sequencestr: str) -> List[Dict[str, str]]:
        """Parse a beta-peptide sequence

        The markup is the following (case sensitive):

            1. beta-2 amino-acid (kind: B2):
                ({S|R})B2h<sc>   e.g. (S)B2hV, (R)B2hA

            2. beta-3 amino-acid (kind: B3):
                ({S|R})B3h<sc>   e.g. (S)B3hL, (R)B3hR

            3. beta-2,3 amino-acid (kind: B23):
                (2{S|R}3{S,R})B23h(2<sc1>3<sc2>)    e.g. (2S3R)B23h(2A3L)

            4. bare beta-backbone (kind: BA):
                BA or (2S3R)B23h(2G3G) (with any combination of R and S
                    in the first pair of parentheses)

            5. alpha-amino acid (kind: A):
                ({S|R|L|D})A<sc>     e.g. (S)AL, (R)AV, (L)AA, (D)AC

            6. capping groups:
                ACE: acetyl on the N-terminus (kind: ACE)
                NME: N-methylamide on the C-terminus (kind: NME)
                BUT: butyl on the N-terminus (kind: BUT)

        Subsequent residues must be separated with comma. Whitespace around commas
        and at the ends of the sequence string will be ignored.

        After each residue a set of backbone torsion angles can be optionally
        supplied in the form:
        [ <angle1> <angle2> <angle3>]
        i.e. a list enclosed in square brackets, separated by whitespace only. The angles
        must be expressed in degrees, in normal notation (no exponential!) For example:

        [-140.0 66 -136]

        For alpha-amino acids two angles must be given, for betas three.

        Employing the secondary structure database, the name of the desired fold can also be given
        in curly braces, e.g.:
        {H14M} or {Straight}

        :param sequence: a string representation of the peptide sequence, following the above markup
        :type sequence: str
        :return: parsed alpha- or beta-amino acids
        :rtype: a list of dicts having the following keys: 'kind', 'stereo2', 'stereo3', 'sidechain2', 'sidechain3'
        :raises ValueError: if an element of a sequence cannot be interpreted
        """
        float_regex = r"""([+-]\s*)?  # optional sign and whitespace
                        (  # start of the mantissa
                        (\d+(\.\d*)?)   # one option: some digits, then optionally some decimals
                        | or
                        (\.\d+)  # second option: only the decimals, led in by a decimal point
                        ) # end of the mantissa
                        ([eE][+-]?\d+)? # optionally, an exponent
                        """

        regex = re.compile(r"""^( # begin
                              (?P<ace>ACE) # acetyl group on the N-terminus
                              |
                              (?P<but>BUT) # butyl group on the N-terminus
                              |
                              (?P<nme>NME) # N-methylamide group on the C-terminus
                              |
                              \(2(?P<achcstereo2>[RS])3(?P<achcstereo3>[RS])\)ACHC  # 2-aminocyclohexanecarboxylic acid
                              | 
                              \(2(?P<acpcstereo2>[RS])3(?P<acpcstereo3>[RS])\)ACPC  # 2-aminocyclopentanecarboxylic acid
                              | 
                              (  # begin beta23-amino acid regex
                              \(2(?P<stereo2>[RS])3(?P<stereo3>[RS])\)  # chirality for beta2,3
                              B23h
                              \(2(?P<sc2>[A-Z]{{1,2}})3(?P<sc3>[A-Z]{{1,2}})\)  # sidechain designation
                              ) # end beta23-amino acid regex
                              | # or
                              ( # begin monosubstituted  beta-amino acid regex
                              \((?P<monostereo>[RS])\) # chirality for mono-substituted beta-amino acids
                              B(?P<monokind>[23])h
                              (?P<monosc>[A-Z]{{1,2}}) # sidechain designation
                              ) # end monostubstituted beta-amino acid regex
                              | # or
                              (?P<bare>BA) # bare beta-amino acid (homo-glycine) 
                              | # or
                              ( # begin alpha-amino acid regex
                              \((?P<alphastereo>[RSLD])\) # chirality for alpha-amino acids
                              A
                              (?P<alphasc>[A-Z]{{1,2}}) # sidechain designation
                              ))
                              \s* # allow for whitespace between the amino-acid designation and the secondary structure info
                              ( # begin optional torsion angle declaration
                              (\[\s*(?P<phi>{0})\s+(?P<theta>{0})(\s+(?P<psi>{0}))?\s*\]) # two or three torsion angles in square brackets 
                              | # or
                              (\{{\s*(?P<ssname>[-a-zA-Z0-9_+ ]+)\s*\}}) # name of a secondary structure database entry in curly braces 
                              )? # end optional torsion angle declaration 
                              $ # end alpha-amino acid regex
                              """.format(float_regex), re.X)

        sequence = []
        for aminoacid in [a.strip() for a in sequencestr.split(',')]:
            m = regex.match(aminoacid)
            if not m:
                raise ValueError('Invalid amino-acid designation: {}'.format(aminoacid))
            if m['phi'] is not None and m['theta'] is not None and m['psi'] is not None:
                dihedrals = [float(m['phi']), float(m['theta']), float(m['psi'])]
            elif m['phi'] is not None and m['theta'] is not None:  # implies that 'psi' not in m
                dihedrals = [float(m['phi']), float(m['theta'])]
            elif m['ssname'] is not None:
                dihedrals = SecondaryStructureDB.dihedrals(m['ssname'])
            else:
                dihedrals = None
            if m['bare'] is not None:
                if dihedrals is not None and len(dihedrals) != 3:
                    raise ValueError('Not enough dihedral angles specified in residue {}'.format(aminoacid))
                sequence.append(
                    {'kind': 'BA', 'sidechain2': 'G', 'stereo2': '', 'sidechain3': 'G', 'stereo3': '',
                     'dihedrals': dihedrals})
            elif m['alphastereo'] is not None:
                if dihedrals is not None and len(dihedrals) != 2:
                    raise ValueError('Too many dihedral angles specified in residue {}'.format(aminoacid))
                sequence.append(
                    {'kind': 'A', 'sidechain2': m['alphasc'], 'stereo2': m['alphastereo'],
                     'sidechain3': '', 'stereo3': '', 'dihedrals': dihedrals})
            elif m['monostereo'] is not None:
                if dihedrals is not None and len(dihedrals) != 3:
                    raise ValueError('Not enough dihedral angles specified in residue {}'.format(aminoacid))
                if m['monokind'] == '2':
                    sequence.append(
                        {'kind': 'B2', 'sidechain2': m['monosc'], 'stereo2': m['monostereo'],
                         'sidechain3': 'G', 'stereo3': '', 'dihedrals': dihedrals})
                else:
                    assert m['monokind'] == '3'
                    sequence.append(
                        {'kind': 'B3', 'sidechain2': 'G', 'stereo2': '',
                         'sidechain3': m['monosc'], 'stereo3': m['monostereo'], 'dihedrals': dihedrals})
            elif m['ace'] is not None:
                if dihedrals is not None:
                    raise ValueError('Residue ACE does not supports custom dihedrals.')
                sequence.append({'kind': 'ACE', 'sidechain2': '', 'stereo2': '', 'sidechain3': '', 'stereo3': '',
                                 'dihedrals': dihedrals})
            elif m['but'] is not None:
                if dihedrals is not None:
                    raise ValueError('Residue BUT does not supports custom dihedrals.')
                sequence.append({'kind': 'BUT', 'sidechain2': '', 'stereo2': '', 'sidechain3': '', 'stereo3': '',
                                 'dihedrals': dihedrals})
            elif m['nme'] is not None:
                if dihedrals is not None:
                    raise ValueError('Residue NME does not supports custom dihedrals.')
                sequence.append({'kind': 'NME', 'sidechain2': '', 'stereo2': '', 'sidechain3': '', 'stereo3': '',
                                 'dihedrals': dihedrals})
            elif m['acpcstereo2'] is not None and m['acpcstereo3'] is not None:
                sequence.append({'kind': 'ACPC', 'sidechain2': '', 'stereo2': m['acpcstereo2'], 'sidechain3': '',
                                 'stereo3': m['acpcstereo3'], 'dihedrals': dihedrals})
            elif m['achcstereo2'] is not None and m['achcstereo3'] is not None:
                sequence.append({'kind': 'ACHC', 'sidechain2': '', 'stereo2': m['achcstereo2'], 'sidechain3': '',
                                 'stereo3': m['achcstereo3'], 'dihedrals': dihedrals})
            else:
                if dihedrals is not None and len(dihedrals) != 3:
                    raise ValueError('Not enough dihedral angles specified in residue {}'.format(aminoacid))
                assert m['stereo2'] is not None
                sequence.append(
                    {'kind': 'B23', 'sidechain2': m['sc2'], 'stereo2': m['stereo2'],
                     'sidechain3': m['sc3'], 'stereo3': m['stereo3'], 'dihedrals': dihedrals})
        return sequence


def betafab2(objname, *args):
    """ Construct an alpha/beta peptide

        Amino acids can be given in the following way:

            1) beta-2 amino-acid:
                ({S|R})B2h<sc>   e.g. (S)B2hV, (R)B2hA

            2) beta-3 amino-acid:
                ({S|R})B3h<sc>   e.g. (S)B3hL, (R)B3hR

            3) beta-2,3 amino-acid:
                (2{S|R}3{S,R})B23h(2<sc1>3<sc2>)    e.g. (2S3R)B23h(2A3L)

            4) bare beta-backbone:
                BA or (2S3R)B23h(2G3G) (with any combination of R and S
                    in the first pair of parentheses)

            5) alpha-amino acid:
                ({S|R|L|D})A<sc>     e.g. (S)AL, (R)AV, (L)AA, (D)AC

            6) capping groups:
                ACE: acetyl on the N-terminus (must be first)
                BUT: butyl on the N-terminus (must be first)
                NME: N-methylamide on the C-terminus (must be last)

            7) cyclic beta-residues:
                (2{S|R}3{S|R})ACHC (2-aminocyclohexanecarboxylic acid) e.g. (2S3R)ACHC
                (2{S|R}3{S|R})ACPC (2-aminocyclopentanecarboxylic acid) e.g. (2S3R)ACPC

        where <sc> is the single-letter abbreviation of a natural
        peptide sidechain. The following are supported:

            A, C, D, E, F, G, I, K, L, M, N, O, Q, R, S, T, V, W, Y.

        In addition, the following variants are also understood:
            CM :
                deprotonated cysteine
            DH :
                aspartic acid
            EH :
                glutamic acid
            HD :
                histidine protonated on the delta-nitrogen
            HE :
                histidine protonated on the epsilon-nitrogen
            HH :
                histidine protonated on both nitrogens
            KN :
                neutral lysine
            ON :
                neutral ornithine

        The desired secondary structure can be given as well after
        each amino-acid in two possible ways:

            1. three (or two) floating point numbers in square brackets, e.g.
               (S)AQ[-57 -47]  or  (S)B3hV[-140.3 66.5 -136.8]
            2. referencing an entry of the secondary structure database in
               curly braces, e.g. (S)AQ{Alpha-helix}   or  (S)B3hV{H14M}

    :param objname: the name PyMOL will know about the resulting peptide
    :type objname: str
    :param args: amino-acid abbreviations
    :type args: str
    :return: the constructed peptide
    :rtype: BetaPeptide
    """
    sequence = [BetaPeptide.parseBetaPeptideSequence(a)[0] for a in args]

    betapeptide = None
    for ires, residue in enumerate(sequence):
        with tempObjectName() as nextname:
            if residue['kind'] == 'ACE':
                nextresidue = BetaPeptide.Ace(nextname)
            elif residue['kind'] == 'BUT':
                nextresidue = BetaPeptide.But(nextname)
            elif residue['kind'] == 'NME':
                nextresidue = BetaPeptide.NMe(nextname)
            elif residue['kind'] == 'ACHC':
                nextresidue = BetaPeptide.ACHC(residue['stereo2'], residue['stereo3'], nextname)
            elif residue['kind'] == 'ACPC':
                nextresidue = BetaPeptide.ACPC(residue['stereo2'], residue['stereo3'], nextname)
            else:
                nextresidue = BetaPeptide(
                    nextname, residue['kind'] != 'A', residue['sidechain2'], residue['stereo2'],
                    residue['sidechain3'], residue['stereo3'])
            if betapeptide is None:
                betapeptide = nextresidue
            else:
                betapeptide += nextresidue
            betapeptide = betapeptide.copy(objname)
    # we are ready. Now we need to fold it. Folding is skipped to this point, because folding a residue needs some
    # atoms from the previous and the next residues, respectively.
    for ires, residue in enumerate(sequence, start=1):
        if residue['dihedrals'] is not None:
            betapeptide.fold(ires, residue['dihedrals'])

    cmd.show_as('sticks', 'model {}'.format(objname))
    return


def betafab2cmd(objname, *args, **kwargs):  # kwargs is needed to swallow the _self argument given by PyMol
    """
    DESCRIPTION

        Build a beta peptide

    USAGE

        betafab2 objname [, aa1 [, aa2 [, aa3 [,...]]]]

    ARGUMENTS

        objname = str: target object name. Will be overwritten!

        aa1, aa2 etc. = str: descriptors of the alpha/beta amino
            acid residues. The following cases are understood:

            1) beta-2 amino-acid:
                ({S|R})B2h<sc>   e.g. (S)B2hV, (R)B2hA

            2) beta-3 amino-acid:
                ({S|R})B3h<sc>   e.g. (S)B3hL, (R)B3hR

            3) beta-2,3 amino-acid:
                (2{S|R}3{S,R})B23h(2<sc1>3<sc2>)    e.g. (2S3R)B23h(2A3L)

            4) bare beta-backbone:
                BA or (2S3R)B23h(2G3G) (with any combination of R and S
                    in the first pair of parentheses)

            5) alpha-amino acid:
                ({S|R|L|D})A<sc>     e.g. (S)AL, (R)AV, (L)AA, (D)AC

            6) capping groups:
                ACE: acetyl on the N-terminus (must be first)
                BUT: butyl on the N-terminus (must be first)
                NME: N-methylamide on the C-terminus (must be last)

            7) cyclic beta-residues:
                (2{S|R}3{S|R})ACHC (2-aminocyclohexanecarboxylic acid) e.g. (2S3R)ACHC
                (2{S|R}3{S|R})ACPC (2-aminocyclopentanecarboxylic acid) e.g. (2S3R)ACPC

        where <sc> is the single-letter abbreviation of a natural
        peptide sidechain. The following are supported:

            A, C, D, E, F, G, I, K, L, M, N, O, Q, R, S, T, V, W, Y.

        In addition, the following variants are also understood:
            CM : deprotonated cysteine
            DH : aspartic acid
            EH : glutamic acid
            HD : histidine protonated on the delta-nitrogen
            HE : histidine protonated on the epsilon-nitrogen
            HH : histidine protonated on both nitrogens
            KN : neutral lysine
            ON : neutral ornithine

        The desired secondary structure can be given as well after
        each amino-acid in two possible ways:

            1. three (or two) floating point numbers in square brackets, e.g.
               (S)AQ[-57 -47]  or  (S)B3hV[-140.3 66.5 -136.8]
            2. referencing an entry of the secondary structure database in
               curly braces, e.g.
               (S)AQ{Alpha-helix}   or  (S)B3hV{H14M}

    EXAMPLE

         betafab2 test, (S)B3hV, (S)B3hA, (R)AK, (S)B3hL, (2S3S)B23h(2A3A)

         betafab2 helix, (S)B3hV{H14M}, (S)B3hA{H14M}, (S)B3hL{H14M}, (2S3S)B23h(2A3A){H14M}, (S)B3hV{H14M}, (S)B3hA{H14M}, (S)B3hL{H14M}

    """
    return betafab2(objname, *args)


if cmd is not None:
    cmd.extend("betafab2", betafab2cmd)
