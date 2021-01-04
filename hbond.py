import warnings

try:
    from pymol import cmd
    import chempy
except ImportError:
    warnings.warn(
        'Cannot import PyMOL: functionality will suffer (you can ignore this if you are just building the documentation).')


def restrain_hbonds_gmx(selection, filename, strength=1000, distance=None,
                        mindist_detect=1.1, maxdist_detect=2.5, anglemin_detect=130, verbose=False):
    """
    DESCRIPTION

        Generate harmonic distance restraints for hydrogen bonds

    USAGE

        restrain_hbonds_gmx selection, filename [, strength [, distance [, mindist_detect [, maxdist_detect [, anglemin_detect ]]]]]]

    ARGUMENTS

        selection: the selection to operate on, containing the donors, the acceptors and the hydrogens

        filename: the file name to write the constraints to (a GROMACS .itp file)

        strength: bond strength (kJ mol-1 nm-2) (default: 1000)

        distance: the hydrogen-acceptor distance to restrain to (default: the current value)

        mindist_detect: minimum hydrogen-acceptor distance to consider (default: 1.1 A)

        maxdist_detect: maximum hydrogen-acceptor distance to consider (default: 2.5 A)

        anglemin_detect: minimum donor-hydrogen-acceptor angle to consider (default: 130°)

    NOTES

        A GROMACS .itp file will be written. The restraint potential has the form:

        V(r) = strength * ( r - distance )^2

    """
    mindist = float(mindist_detect)
    maxdist = float(maxdist_detect)
    anglemin = float(anglemin_detect)
    strength = float(strength)
    if distance is not None:
        distance = float(distance)

    with open(filename, 'wt') as f:
        f.write('; hydrogen bond distance restraints\n')
        f.write('[ bonds ]\n')
        f.write(';       ai        aj  funct        b0     kb       b0_B   kb_B\n')
        for donor, hydrogen, acceptor, dist, angle in find_beta_hbonds(selection, mindist, maxdist, anglemin, verbose):
            f.write(
                '{:10d}{:10d}   6   {:10.4f} {:.6f} {:10.4f} {:.6f}; original distance: {:.6f} nm, original angle: {:.2f} degrees\n'.format(
                    hydrogen, acceptor, distance / 10. if distance is not None else dist / 10., strength,
                    distance / 10. if distance is not None else dist / 10., 0.0, dist / 10., angle))


def find_beta_hbonds(selection, mindist, maxdist, anglemin, verbose: bool = False, nitrogen_atomnames='N+NT',
                     oxygen_atomnames='O+OT1+OT2', hydrogen_atomnames='H+HN+HT1+HT2+HT+H1+H2+H3'):
    """
    Detect hydrogen bonds in beta peptides (and probably other peptides, too)

    Inputs:
        selection: a selection containing everything,

    """
    donors = cmd.get_model('({}) and name {}'.format(selection, nitrogen_atomnames))
    acceptors = cmd.get_model('({}) and name {}'.format(selection, oxygen_atomnames))
    if verbose:
        print('Finding hydrogen bonds in selection {}'.format(selection))
        print('Donors: {}'.format(len(donors.atom)))
        print('Acceptors: {}'.format(len(acceptors.atom)))
    for d in donors.atom:
        assert isinstance(d, chempy.Atom)
        hydrogens = cmd.get_model(
            '(neighbor (({}) and index {})) and name {}'.format(selection, d.index, hydrogen_atomnames))
        if verbose:
            print('Looking at donor {}, having {} hydrogens'.format(d.index, len(hydrogens.atom)))
        for h in hydrogens.atom:
            assert isinstance(h, chempy.Atom)
            for a in acceptors.atom:
                assert isinstance(a, chempy.Atom)
                # check hydrogen-acceptor distance
                dist = cmd.get_distance('({}) and name {} and index {}'.format(selection, oxygen_atomnames, a.index),
                                        '({}) and name {} and index {}'.format(selection, hydrogen_atomnames, h.index))
                cmd.unpick()
                angle = cmd.get_angle(
                    '({}) and name {} and index {}'.format(selection, oxygen_atomnames, a.index),
                    '({}) and name {} and index {}'.format(selection, hydrogen_atomnames, h.index),
                    '({}) and name {} and index {}'.format(selection, nitrogen_atomnames, d.index))
                if verbose:
                    print('Distance between hydrogen {} and acceptor {}: {}'.format(h.index, a.index, dist))
                    print('Angle between donor {}, hydrogen {} and acceptor {}: {}'.format(d.index, h.index, a.index,
                                                                                           angle))
                if (dist <= float(maxdist)) and (dist >= float(mindist)) and (angle >= float(anglemin)):
                    yield d.index, h.index, a.index, dist, angle
