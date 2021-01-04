import typing

from PyQt5 import QtCore, QtGui

from ..betafab2 import BetaPeptide, BetaPeptideConsistencyError
from ..secstructdb import SecondaryStructureDB


class AminoAcidResidue:
    def __init__(self, kind: str, stereo2: str, sidechain2: str, stereo3: str, sidechain3: str,
                 dihedrals: typing.Union[typing.List[float], str, None]):
        if isinstance(dihedrals, str):
            dihedrals = [d for d in SecondaryStructureDB.dihedrals(dihedrals) if d is not None]
        if kind not in ['BA', 'B2', 'B3', 'B23', 'A', 'ACE', 'NME', 'BUT', 'ACPC', 'ACHC']:
            raise ValueError('Unknown residue kind: {}'.format(kind))
        self._kind = kind
        if self._kind in ['ACE', 'NME', 'BUT']:
            self._stereo2 = ''
            self._stereo3 = ''
            self._sidechain3 = ''
            self._sidechain2 = ''
            self._dihedrals = None
        else:
            if stereo2 not in ['R', 'S', 'D', 'L'] and kind in ['A']:
                raise ValueError('Unknown absolute conformation designation: {}'.format(stereo2))
            elif stereo2 not in ['R', 'S'] and kind in ['B2', 'B23', 'ACPC', 'ACHC']:
                raise ValueError('Unknown absolute conformation designation: {}'.format(stereo2))
            self._stereo2 = stereo2
            if stereo3 not in ['R', 'S'] and kind in ['B3', 'B23', 'ACPC', 'ACHC']:
                raise ValueError('Unknown absolute conformation designation: {}'.format(stereo3))
            self._stereo3 = stereo3
            if sidechain2 not in BetaPeptide.SIDECHAINS and kind in ['B2', 'B23', 'A']:
                raise ValueError('Unknown sidechain: {}'.format(sidechain2))
            self._sidechain2 = sidechain2
            if sidechain3 not in BetaPeptide.SIDECHAINS and kind in ['B3', 'B23']:
                raise ValueError('Unknown sidechain: {}'.format(sidechain3))
            self._sidechain3 = sidechain3
            if kind in ['B2', 'B23', 'B3', 'BA', 'ACPC', 'ACHC'] and (
                    (dihedrals is not None) and (len(dihedrals) != 3)):
                raise ValueError('Residue of kind {} needs 3 dihedrals!'.format(self._kind))
            if kind in ['A'] and ((dihedrals is not None) and (len(dihedrals) != 2)):
                raise ValueError('Residue of kind {} needs 2 dihedrals!'.format(self._kind))
            self._dihedrals = dihedrals

    @property
    def stereo2(self) -> str:
        """Chirality of the alpha-sidechain"""
        return self._stereo2 if self._kind in ['B2', 'B23', 'A', 'ACPC', 'ACHC'] else ''

    @stereo2.setter
    def stereo2(self, value: str):
        if value not in (['R', 'S', 'L', 'D', ''] if self._kind is 'A' else ['R', 'S', '']):
            raise ValueError('Invalid chirality: {}'.format(value))
        self._stereo2 = value

    @property
    def stereo3(self) -> str:
        """Chirality of the beta-sidechain"""
        return self._stereo3 if self._kind in ['B3', 'B23', 'ACPC', 'ACHC'] else ''

    @stereo3.setter
    def stereo3(self, value: str):
        if value not in ['R', 'S', '']:
            raise ValueError('Invalid chirality: {}'.format(value))
        self._stereo3 = value

    @property
    def sidechain2(self) -> str:
        """Type of the alpha-sidechain"""
        return self._sidechain2 if self._kind in ['B2', 'B23', 'A'] else ''

    @sidechain2.setter
    def sidechain2(self, value: str):
        if (value not in BetaPeptide.SIDECHAINS) and (value):
            raise ValueError('Unknown sidechain: {}'.format(value))
        self._sidechain2 = value

    @property
    def sidechain3(self) -> str:
        """Type of the beta-sidechain"""
        return self._sidechain3 if self._kind in ['B3', 'B23'] else ''

    @sidechain3.setter
    def sidechain3(self, value: str):
        if (value not in BetaPeptide.SIDECHAINS) and (value):
            raise ValueError('Unknown sidechain: {}'.format(value))
        self._sidechain3 = value

    @property
    def kind(self) -> str:
        return self._kind

    @kind.setter
    def kind(self, value: str):
        if value not in ['BA', 'B2', 'B3', 'B23', 'A', 'ACE', 'NME', 'BUT', 'ACPC', 'ACHC']:
            raise ValueError('Unknown residue type: {}'.format(value))
        self._kind = value
        if self._kind in ['B2', 'A', 'B23', 'ACPC', 'ACHC']:
            if not self._stereo2:
                self._stereo2 = 'S'
            if not self._sidechain2:
                self._sidechain2 = 'A'
        if self._kind in ['B3', 'B23', 'ACPC', 'ACHC']:
            if not self._stereo3:
                self._stereo3 = 'S'
            if not self._sidechain3:
                self._sidechain3 = 'A'
        if self._kind in ['B3', 'B23', 'ACPC', 'ACHC', 'B2']:
            # if this residue was an alpha-amino acid before and the chirality was specified with D/L, change it to R/S
            if self._stereo2 == 'L':
                self._stereo2 = 'S'
            elif self._stereo2 == 'D':
                self._stereo2 = 'R'
        if self._kind in ['A']:
            if not self._dihedrals or len(self._dihedrals) != 2:
                self._dihedrals = [180, 180]
        elif self._kind in ['BA', 'B2', 'B3', 'B23', 'ACPC', 'ACHC']:
            if not self._dihedrals or len(self._dihedrals) != 3:
                self._dihedrals = [180, 180, 180]
        else:
            if self._dihedrals:
                self._dihedrals = None

    @property
    def dihedrals(self):
        return self._dihedrals

    @dihedrals.setter
    def dihedrals(self, value: typing.Optional[typing.List[float]]):
        if self._kind in ['A'] and (value is None or len(value) != 2):
            raise ValueError('Alpha-amino acids need two dihedral angles')
        if self._kind in ['BA', 'B2', 'B3', 'B23', 'ACHC', 'ACPC'] and (value is None or len(value) != 3):
            raise ValueError('Beta-amino acids need three dihedral angles')
        if self._kind in ['ACE', 'NME', 'BUT'] and (value is not None):
            raise ValueError('ACE, BUT and NME residues do not need dihedral angles.')
        self._dihedrals = value

    def __str__(self, with_dih: bool = True) -> str:
        if self._dihedrals is None:
            dih = ''
        else:
            dih = '[' + ' '.join(['{:f}'.format(d) for d in self._dihedrals if d is not None]) + ']'
        if not with_dih:
            dih = ''
        if self.kind == 'BA':
            return 'BA' + dih
        elif self.kind == 'B2':
            return '({0.stereo2})B2h{0.sidechain2}'.format(self) + dih
        elif self.kind == 'B3':
            return '({0.stereo3})B3h{0.sidechain3}'.format(self) + dih
        elif self.kind == 'B23':
            return '(2{0.stereo2}3{0.stereo3})B23h(2{0.sidechain2}3{0.sidechain3})'.format(self) + dih
        elif self.kind == 'ACPC':
            return '(2{0.stereo2}3{0.stereo3})ACPC'.format(self) + dih
        elif self.kind == 'ACHC':
            return '(2{0.stereo2}3{0.stereo3})ACHC'.format(self) + dih
        elif self.kind == 'A':
            return '({0.stereo2})A{0.sidechain2}'.format(self) + dih
        elif self.kind == 'ACE':
            return 'ACE'
        elif self.kind == 'NME':
            return 'NME'
        elif self.kind == 'BUT':
            return 'BUT'
        else:
            raise ValueError('Invalid kind: {}'.format(self.kind))


class SequenceModel(QtCore.QAbstractItemModel):
    """A model for storing and manipulating alpha/beta-peptide sequence

    The columns are the following:

    0. index (integer), not user-editable
    1. abbreviation (string), not user-editable
    2. type (string), user-editable (combo box)
    3. stereo2 (string), user-editable (combo box)
    4. sidechain2 (string), user-editable (combo box)
    5. stereo3 (string), user-editable (combo box)
    6. sidechain3 (string), user-editable (combo box)
    7. secondary structure, user-editable
    """

    validationChanged = QtCore.pyqtSignal(bool)

    def __init__(self):
        super().__init__(parent=None)
        self._residues = []
        self._failed = []

    def rowCount(self, parent: QtCore.QModelIndex = None):
        """Return the number of rows"""
        if parent is None:
            parent = QtCore.QModelIndex()
        return len(self._residues) if not parent.isValid() else 0

    def columnCount(self, parent: QtCore.QModelIndex = None):
        """Return the number of columns"""
        if parent is None:
            parent = QtCore.QModelIndex()
        return 8 if not parent.isValid() else 0

    def flags(self, index: QtCore.QModelIndex):
        """Flags influencing the display and behaviour of the field"""
        basicflags = QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemNeverHasChildren | QtCore.Qt.ItemIsEnabled
        if index.column() in [2, 3, 4, 5, 6, 7] and self._residues[index.row()].kind in ['A', 'BA', 'B2', 'B3', 'B23']:
            return basicflags | QtCore.Qt.ItemIsEditable
        if index.column() in [2, 3, 5, 7] and self._residues[index.row()].kind in ['ACPC', 'ACHC']:
            return basicflags | QtCore.Qt.ItemIsEditable
        else:
            return basicflags

    def data(self, index: QtCore.QModelIndex, role: int = ...):
        """Return data for displaying items"""
        if (role == QtCore.Qt.DisplayRole) or (
                role == QtCore.Qt.EditRole):  # control the display behaviour of the items.
            # incidentally, the same information has to be given in the edit role.
            if index.column() == 0:
                # residue number
                return str(index.row() + 1)
            elif index.column() == 1:
                # residue name
                return self._residues[index.row()].__str__(with_dih=False)
            elif index.column() == 2:
                # residue kind
                return self._residues[index.row()].kind
            elif index.column() == 3:
                # stereo2
                return self._residues[index.row()].stereo2
            elif index.column() == 4:
                # sidechain2
                return self._residues[index.row()].sidechain2
            elif index.column() == 5:
                # stereo3
                return self._residues[index.row()].stereo3
            elif index.column() == 6:
                # sidechain3
                return self._residues[index.row()].sidechain3
            elif index.column() == 7:
                # sec.structure
                dihs = self._residues[index.row()].dihedrals
                if dihs is None:
                    return ''
                elif len(dihs) == 2:
                    return SecondaryStructureDB.find(dihs[0], None, dihs[1])
                else:
                    assert len(dihs) == 3
                    return SecondaryStructureDB.find(dihs[0], dihs[1], dihs[2])
            else:
                raise ValueError('Invalid column: {}'.format(index.column()))
        elif role == QtCore.Qt.TextColorRole:
            if index.row() in [f[0] for f in self._failed]:
                return QtGui.QColor('red')
            else:
                return None
        elif role == QtCore.Qt.ToolTipRole:
            if index.row() in [f[0] for f in self._failed]:
                return "This row is incorrect for the following reason(s):\n{}".format(
                    '\n'.join(['   - ' + f[1] for f in self._failed if f[0] == index.row()]))
            else:
                return "This row is validated."
        else:
            # do nothing for the other roles
            return None

    def setData(self, index: QtCore.QModelIndex, value: typing.Any, role: int = ...) -> bool:
        if role == QtCore.Qt.EditRole:
            if index.column() == 2:
                self._residues[index.row()].kind = value
            elif index.column() == 3:
                self._residues[index.row()].stereo2 = value
            elif index.column() == 4:
                self._residues[index.row()].sidechain2 = value
            elif index.column() == 5:
                self._residues[index.row()].stereo3 = value
            elif index.column() == 6:
                self._residues[index.row()].sidechain3 = value
            elif index.column() == 7:
                self._residues[index.row()].dihedrals = [d for d in SecondaryStructureDB.dihedrals(value) if
                                                         d is not None]
            self.dataChanged.emit(self.index(index.row(), index.column()), self.index(index.row(), index.column()))
            self.validate(noraise=True)
            return True
        return False

    def headerData(self, section: int, orientation: QtCore.Qt.Orientation, role: int = ...):
        if (role == QtCore.Qt.DisplayRole) and orientation == QtCore.Qt.Horizontal:
            # only the horizontal header is of matter to us, and only the display role.
            return \
                ['#', 'Residue', 'Type', 'Cα chirality', 'α side-chain', 'Cβ chirality', 'β side-chain', 'Sec.struct.'][
                    section]
        else:
            return None

    def parent(self, child: QtCore.QModelIndex) -> QtCore.QModelIndex:
        """Return the parent item for an index

        Because this is a flat model, this function always returns an invalid index.

        :param child:
        :type child:
        :return: an invalid QModelIndex index
        :rtype: QtCore.QModelIndex
        """
        return QtCore.QModelIndex()

    def index(self, row: int, column: int, parent: QtCore.QModelIndex = None) -> QtCore.QModelIndex:
        """Create an index for the item in (row, column)

        :param row: the row number of the item
        :type row: int
        :param column: the column number of the item
        :type column: int
        :param parent: index of the parent item, not used in this model
        :type parent: QtCore.QModelIndex or None
        :return: the index of the item
        :rtype: QtCore.QModelIndex
        """
        if parent is None:
            parent = QtCore.QModelIndex()
        return self.createIndex(row, column, None) if not parent.isValid() else QtCore.QModelIndex()

    def append(self, kind: str, stereo2: str, sidechain2: str, stereo3: str, sidechain3: str,
               dihedrals: typing.Union[typing.List[float], str, None]):
        """Append a new residue to the model.

        :param kind: amino-acid kind
        :type kind: str
        :param stereo2: chirality of the alpha sidechain
        :type stereo2: str
        :param sidechain2: alpha sidechain type
        :type sidechain2: str
        :param stereo3: chirality of the beta sidechain
        :type stereo3: str
        :param sidechain3: beta sidechain type
        :type sidechain3: str
        :param dihedrals: dihedral angles
        :type dihedrals: list of floats or a string or None
        """
        self.beginInsertRows(QtCore.QModelIndex(), len(self._residues), len(self._residues) + 1)
        self._residues.append(AminoAcidResidue(kind, stereo2, sidechain2, stereo3, sidechain3, dihedrals))
        self.endInsertRows()
        self.validate(noraise=True)

    def insert(self, row: int, kind: str, stereo2: str, sidechain2: str, stereo3: str, sidechain3: str,
               dihedrals: typing.Union[typing.List[float], str, None]):
        """Insert a new residue to the model before a given row

        :param row: index of the row. The new residue will be on this row at the end.
        :type row: int
        :param kind: amino-acid kind
        :type kind: str
        :param stereo2: chirality of the alpha sidechain
        :type stereo2: str
        :param sidechain2: alpha sidechain type
        :type sidechain2: str
        :param stereo3: chirality of the beta sidechain
        :type stereo3: str
        :param sidechain3: beta sidechain type
        :type sidechain3: str
        :param dihedrals: dihedral angles
        :type dihedrals: list of floats or a string or None
        """
        self.beginInsertRows(QtCore.QModelIndex(), row, row + 1)
        self._residues.insert(row, AminoAcidResidue(kind, stereo2, sidechain2, stereo3, sidechain3, dihedrals))
        self.endInsertRows()
        self.validate(noraise=True)

    def removeRow(self, row: int, parent: QtCore.QModelIndex = None):
        """Remove a single row from the model

        :param row: the row to remove
        :type row: int
        :param parent: index of the parent item; not used in this model
        :type parent: QtCore.QModelIndex or None
        """
        if parent is None:
            parent = QtCore.QModelIndex()
        self.beginRemoveRows(QtCore.QModelIndex(), row, row)
        del self._residues[row]
        self.endRemoveRows()
        self.validate(noraise=True)

    def removeRows(self, row: int, count: int, parent: QtCore.QModelIndex = None):
        """Remove a continuous range of rows from the model

        :param row: the first row to remove
        :type row: int
        :param count: the number of rows to remove
        :type count: int
        :param parent: index of the parent item; not used in this model
        :type parent: QtCore.QModelIndex or None
        """
        if parent is None:
            parent = QtCore.QModelIndex()
        self.beginRemoveRows(QtCore.QModelIndex(), row, row + count - 1)
        del self._residues[row:row + count]
        self.endRemoveRows()
        self.validate(noraise=True)

    def clear(self):
        """Clear the model: remove all residues"""
        self.beginResetModel()
        self._residues = []
        self.endResetModel()
        self.validate(noraise=True)

    def moveDown(self, row: int) -> int:
        """Move row up by one, if possible.

        :param row: the index of the row to move
        :type row: int
        :return: the index of the moved residue after moving.
        :rtype: int
        """
        if row >= len(self._residues) - 1:
            return row
        residue = self._residues[row]
        self.beginResetModel()
        del self._residues[row]
        self._residues.insert(row + 1, residue)
        self.endResetModel()
        self.validate(noraise=True)
        return row + 1

    def moveUp(self, row: int) -> int:
        """Move row down by one, if possible.

        :param row: the index of the row to move
        :type row: int
        :return: the index of the moved residue after moving.
        :rtype: int
        """
        if row <= 0:
            return row
        residue = self._residues[row]
        self.beginResetModel()
        del self._residues[row]
        self._residues.insert(row - 1, residue)
        self.endResetModel()
        self.validate(noraise=True)
        return row - 1

    def moveToTop(self, row: int) -> int:
        """Move row to the top, if possible.

        :param row: the index of the row to move
        :type row: int
        :return: the index of the moved residue after moving.
        :rtype: int
        """
        residue = self._residues[row]
        self.beginResetModel()
        del self._residues[row]
        self._residues.insert(0, residue)
        self.endResetModel()
        self.validate(noraise=True)
        return 0

    def moveToBottom(self, row: int) -> int:
        """Move row to the bottom, if possible.

        :param row: the index of the row to move
        :type row: int
        :return: the index of the moved residue after moving.
        :rtype: int
        """
        residue = self._residues[row]
        self.beginResetModel()
        del self._residues[row]
        self._residues.append(residue)
        self.endResetModel()
        self.validate(noraise=True)
        return len(self._residues) - 1

    def residues(self) -> typing.List[str]:
        """Get the list of residues in this model

        :return: the list of residues in a format acceptable by :func:`pmlbeta.betafab2.betafab2`
        :rtype: list of strings
        """
        return [str(res) for res in self._residues]

    def validate(self, noraise: bool = False):
        """Validate the sequence. Currently the following are tested:

        1. ACE and BUT can only be the first residue
        2. NME can only be the last residue

        :param noraise: if raising an exception is to be avoided
        :type noraise: bool
        :return: a list of failing rows and the reasons. One row may append multiple times with different reasons
        :rtype: list of tuples of (row, reason)
        :raise BetaPeptideConsistencyError: if the validation fails
        """
        failed = []
        for idx in [i for i in range(1, len(self._residues)) if self._residues[i].kind == 'ACE']:
            failed.append((idx, 'ACE not at the N terminus'))
        for idx in [i for i in range(1, len(self._residues)) if self._residues[i].kind == 'BUT']:
            failed.append((idx, 'BUT not at the N terminus'))
        for idx in [i for i in range(0, len(self._residues) - 1) if self._residues[i].kind == 'NME']:
            failed.append((idx, 'NME not at the C terminus'))
        self._failed = failed
        for idx in [f[0] for f in self._failed]:
            self.dataChanged.emit(self.index(idx, 0), self.index(idx, self.columnCount()),
                                  (QtCore.Qt.TextColorRole, QtCore.Qt.ToolTipRole))
        self.validationChanged.emit(self.isValid())
        if failed and not noraise:
            raise BetaPeptideConsistencyError(failed[0][1])
        return self._failed

    def isValid(self) -> bool:
        """If the sequence is valid

        :return: True if the sequence is buildable
        :rtype: bool
        """
        return not self._failed

    def residue(self, index: int) -> AminoAcidResidue:
        """Get the `index`-th residue

        :param index: the 0-based index of the residue
        :type index: int
        :return: the residue
        :rtype: AminoAcidResidue
        """
        return self._residues[index]
