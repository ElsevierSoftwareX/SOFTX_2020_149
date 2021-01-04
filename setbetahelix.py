import warnings

try:
    from pymol import cmd
except ImportError:
    warnings.warn(
        'Cannot import PyMOL: functionality will suffer (you can ignore this if you are just building the documentation).')
import itertools
from .utils import set_dihedral
from typing import Union, List, Tuple, Optional, Iterable
from .secstructdb import SecondaryStructureDB
import re


def is_amino_acid(selection: str) -> bool:
    """Decide if the selection is a single amino-acid

    :param selection: the selection
    :type selection: str
    :return: True or False
    :rtype: bool
    """
    if not cmd.count_atoms('({}) and name N and symbol N'.format(selection)) == 1:
        # either does not contain a nitrogen called "N" (no amino-acid), or
        # contains more of them (not a single residue)
        return False
    elif not cmd.count_atoms('({}) and name C and symbol C'.format(selection)) == 1:
        # either does not contain a carbon called "C" (no amino-acid),
        # or contains more of them (not a single residue)
        return False
    elif not cmd.count_atoms(
            '({}) and (neighbor (name C and symbol C)) and (name CA and symbol C)'.format(selection)) == 1:
        # the "C" atom must have a "CA" neighbour
        return False
    else:
        return True


def is_beta(selection: str) -> bool:
    """Decide if a residue is a beta-amino acid

    :param selection: the selection to analyze
    :type selection: str
    :return: True or False
    :rtype: bool
    """
    if not is_amino_acid(selection):
        return False
    # see if "CA" has a neighbour which is also the neighbour of "N" and has the name "CB" or "CB1"
    if cmd.count_atoms(
            "(neighbor (({0}) and name N)) and (neighbor (({0}) and name CA)) and ({0}) and name CB+CB1".format(
                selection)) == 1:
        return True
    return False


def is_alpha(selection: str) -> bool:
    """Decide if a residue is an alpha-amino acid

    :param selection: the selection to analyze
    :type selection: str
    :return: True or False
    :rtype: bool
    """
    if not is_amino_acid(selection):
        return False
    # see if "CA" and "N" are neighbours
    if cmd.count_atoms("(neighbor (({0}) and name N)) and ({0}) and name CA".format(selection)) == 1:
        return True
    return False


def fold_bp(sstype: str, selection: str = '(all)'):
    """
    DESCRIPTION

    Adjust the torsion angles of an alpha- or beta-peptide to
    fold it into the desired secondary structure

    USAGE

    fold_bp sstype [, selection]

    ARGUMENTS

    sstype = string: the desired secondary structure. The following possibilities exist:
        1. a name from the secondary structure database (see command ssdb_list)
        2. two or three (for alpha- and beta-amino acids, respectively) space-separated
           floating point numbers in parentheses, corresponding to the backbone dihedral angles

        The such given secondary structure will be applied to all residues in the selection.

        Additionally, a space-separated list of the same length as the residues in the
        selection can be given in square brackets, containing entry names or angle tuples
        (or a mix of them), corresponding to the residues.

    selection = the selection to operate on. Must be a single peptide chain with
        unique and consecutive residue IDs (default: all)

    EXAMPLES

        fold_bp H14M, model valxval

        fold_bp (-140.3 66.5 -136.8), model valxval

        fold_bp [(-140.3 66.5 -136.8) (180 180 180) H14M], model tripeptide

    SEE ALSO

        ssdb_add, ssdb_del, ssdb_dihedrals, ssdb_resetdefaults, ssdb_list
    """
    if isinstance(sstype, str):
        # parse the string into a list of tuples of floats or Nones
        sstype = parse_sstype(sstype)
    elif isinstance(sstype, tuple):
        # we should have a tuple of floats or Nones
        if not all([isinstance(x, float) or x is None for x in sstype]):
            raise ValueError('Error in secondary structure type: {}'.format(sstype))
        sstype = itertools.cycle([sstype])  # use this for all residues
    # otherwise try to use sstype as an iterable, producing tuples of floats or Nones

    residues = list(sorted({a.resi_number for a in cmd.get_model(selection).atom}))
    r = None
    for r, angles in zip(residues, sstype):
        phi, theta, psi = angles
        if is_beta('({}) and resi {}'.format(selection, r)):
            set_dihedral(selection, ('C', r - 1), ('N', r), ('CB+CB1', r), ('CA', r), phi)
            set_dihedral(selection, ('N', r), ('CB+CB1', r), ('CA', r), ('C', r), theta)
            set_dihedral(selection, ('CB+CB1', r), ('CA', r), ('C', r), ('N', r + 1), psi)
        elif is_alpha('({}) and resi {}'.format(selection, r)):
            set_dihedral(selection, ('C', r - 1), ('N', r), ('CA', r), ('C', r), phi)
            set_dihedral(selection, ('N', r), ('CA', r), ('C', r), ('N', r + 1), psi)
        else:
            # not an amino acid, do nothing with this.
            continue
    if not r == max(residues):
        # after the for loop, r must be the largest residue number. If this is not the case,
        # too few angle triplets were given. Do nothing at present, just warn the user
        print('Warning: not all residues have been processed (too few secondary structures given)')

    cmd.unpick()
    cmd.orient(selection)


def parse_sstype(sstype: Union[
    str, Tuple[float, Optional[float], float], List[Union[str, Tuple[float, Optional[float], float]]]]) -> Iterable[
    Tuple[float, Optional[float], float]]:
    """Parse the secondary structure information for beta- and alpha-peptides

    The desired secondary structure type can be given in two ways:
        1. three floating point numbers in parentheses, separated by space:
           the backbone torsion angles. The middle one can be None for alpha-amino acids.

        2. a name of an entry in the secondary structure database (simple string)

    Additionally, a list can be given from any of the above, inside square brackets and
    separated by whitespace. E.g.:

        [ (-120 80 -136) H14M (130 50 45) ... ]

    The thing is more complicated because PyMOL gives all parameters in string format.
    """

    # now define some regular expressions which we will use

    # a regular expression for floating point numbers, including the exponential form
    float_regex = r"""([+-]\s*)?  # optional sign and whitespace
                    (  # start of the mantissa
                    (\d+(\.\d*)?)   # one option: some digits, then optionally some decimals
                    | or
                    (\.\d+)  # second option: only the decimals, led in by a decimal point
                    ) # end of the mantissa
                    ([eE][+-]?\d+)? # optionally, an exponent
                    """

    # a regular expression of a parenthesized part (without nested parentheses)
    paren_regex = re.compile(r"\([^()]*\)")

    sstuple_regex = re.compile(r""" #regular expression for dihedral angle tuples
    (
    \(\s*(?P<angle1>{0})\s*(?P<angle2>{0})\s*((?P<angle3>{0})\s*)?\)   # two or three floats in parentheses, whitespace separated 
    )
    | # or 
    (
    \(\s*(?P<phi>{0})\s*(?P<theta>{0}|None)\s*(?P<psi>{0})\s*\) # three floats, the middle can be None, in parentheses, whitespace separated
    )
    """.format(float_regex),
                               re.VERBOSE)

    def parse_parentuple(parenpart: str) -> Tuple[float, Optional[float], float]:
        """Parse a parenthesized part of a sec.structure definition string"""
        m = sstuple_regex.match(parenpart)
        if m is None:
            raise ValueError('Invalid dihedral angle tuple: {}'.format(parenpart))
        if m['angle1'] is not None and m['angle2'] is not None and m['angle3'] is not None:
            # three angles
            return (float(m['angle1']), float(m['angle2']), float(m['angle3']))
        elif m['angle1'] is not None and m['angle2'] is not None and m['angle3'] is None:
            # two angles
            return (float(m['angle1']), None, float(m['angle2']))
        elif m['angle1'] is not None or m['angle2'] is not None or m['angle3'] is not None:
            raise ValueError('Invalid dihedral angle tuple: {}'.format(parenpart))
        elif m['phi'] is not None and m['theta'] == 'None' and m['psi'] is not None:
            return (float(m['phi']), None, float(m['psi']))
        elif m['phi'] is not None and m['theta'] is not None and m['psi'] is not None:
            return (float(m['phi']), float(m['theta']), float(m['psi']))
        else:
            raise ValueError('Invalid dihedral angle tuple: {}'.format(parenpart))

    # now start working in earnest.

    if isinstance(sstype, tuple):
        # this must be a tuple of three floats, the middle one can be None
        if len(sstype) != 3:
            raise ValueError('Invalid tuple length: {}'.format(len(sstype)))
        elif not (isinstance(sstype[0], float) and
                  (isinstance(sstype[1], float) or sstype[1] is None) and
                  (isinstance(sstype[2], float))):
            raise ValueError(
                'Incorrect torsion angle type(s): {}, {}, {}'.format(type(sstype[0]), type(sstype[1]), type(sstype[2])))
        else:
            return itertools.cycle([sstype])  # one-applies-for-all case
    elif isinstance(sstype, list):
        # see if every element of the list is a correct tuple or a string. If a string, find the corresponding entry
        # in the secondary structure database.
        sstype_parsed = []
        for item in sstype:
            if isinstance(item, str):
                # this can fail with KeyError
                item = SecondaryStructureDB.dihedrals(item)
            elif isinstance(item, tuple):
                # validate it by simply calling ourselves again
                item = parse_sstype(item)[0]
            else:
                raise TypeError('Invalid type: {}, ({})'.format(item, type(item)))
            sstype_parsed.append(item)
        return sstype_parsed  # list of tuples
    elif isinstance(sstype, str):
        # this is the most complicated part.
        sstype = sstype.strip()  # strip leading whitespaces

        # some validation (not complete, i.e. mismatched/mixed up parentheses or brackets are not detected
        if sstype.count('[') != sstype.count(']'):  # 0==0 can also happen
            raise ValueError('Unmatched brackets in secondary structure type {}'.format(sstype))
        if sstype.count('(') != sstype.count(')'):  # 0==0 can also happen
            raise ValueError('Unmatched parentheses in secondary structure type {}'.format(sstype))
        if sstype.startswith('[') and sstype.endswith(']'):  # parse it as a list
            sstype = sstype[1:-1].strip()  # cut the brackets
            # first find the parts in parentheses
            parentheses = [(x.start(), x.end()) for x in paren_regex.finditer(sstype)]
            items = []

            def parse_outsideofparentheses_part(line: str):
                # If "x" does not contain parentheses:
                #   1. strip it
                #   2. split it up at whitespaces
                #   3. find the dihedrals in the SSDB
                # print('Parsing outside parentheses: *{}*'.format(line))
                return [SecondaryStructureDB.dihedrals(x) for x in line.strip().split()]

            if not parentheses:
                # parse the whole string as a list of entry names
                # print('No parentheses found.')
                items = parse_outsideofparentheses_part(sstype)
            else:
                # parse each pair of parentheses and the parts inbetween them.
                for ip in range(len(parentheses)):
                    parenpart = sstype[parentheses[ip][0]:parentheses[ip][1]]
                    # print('This parenthesis pair: {} to {} = *{}*'.format(*parentheses[ip], parenpart))
                    if ip == 0:
                        # parse the part before the first pair of parentheses
                        items.extend(parse_outsideofparentheses_part(sstype[0:parentheses[0][0]]))
                    else:
                        # parse the part before this pair of parentheses and after the previous pair of parentheses
                        items.extend(parse_outsideofparentheses_part(sstype[parentheses[ip - 1][1]:parentheses[ip][0]]))
                    # now parse the parenthesized part
                    items.append(parse_parentuple(parenpart))
                # parse the end of the string after the last pair of parentheses
                items.extend(parse_outsideofparentheses_part(sstype[parentheses[-1][1]:]))
            return items  # list of tuples
        elif '[' in sstype:
            raise ValueError(
                'Opening bracket must be the first character in the secondary structure type {}'.format(sstype))
        elif sstype.startswith('(') and sstype.endswith(')'):
            return itertools.cycle([parse_parentuple(sstype)])  # one-applies-for-all case
        elif '(' in sstype:
            raise ValueError(
                'Opening parenthesis must be the first character in the secondary structure type {}'.format(sstype))
        else:
            # no [, no ], no ( and no ) in the string, it must be an entry in the secondary structure database
            return itertools.cycle([SecondaryStructureDB.dihedrals(sstype)])  # one-applies-for-all case
