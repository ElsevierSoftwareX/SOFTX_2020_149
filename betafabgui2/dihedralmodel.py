"""Model for storing secondary structures and the corresponding dihedral angles"""

import typing

from PyQt5 import QtCore

from ..secstructdb import SecondaryStructureDB


class DihedralModel(QtCore.QAbstractItemModel):
    """A model for storing secondary structures with their dihedral angles

    Columns:
        1. name
        2. type (alpha/beta)
        3. phi
        4. theta
        5. psi
    """

    def __init__(self):
        super().__init__(None)
        self._data = {}

    def columnCount(self, parent: typing.Optional[QtCore.QModelIndex] = None) -> int:
        """Return the number of columns

        :param parent: the parent index (not used, this is a flat model)
        :type parent: QtCore.QModelIndex or None
        :return: the number of columns
        :rtype: int
        """
        if parent is None:
            parent = QtCore.QModelIndex()
        return 5 if not parent.isValid() else 0

    def rowCount(self, parent: QtCore.QModelIndex = None) -> int:
        """Return the number of rows

        :param parent: the parent index (not used, this is a flat model)
        :type parent: QtCore.QModelIndex or None
        :return: the number of rows
        :rtype: int
        """
        if parent is None:
            parent = QtCore.QModelIndex()
        return len(self._data) if not parent.isValid() else 0

    def parent(self, child: QtCore.QModelIndex) -> QtCore.QModelIndex:
        """Get the parent of a child index

        :param child: index
        :type child: QtCore.QModelIndex
        :return: always an invalid QModelIndex, because this is a flat model
        :rtype: QtCore.QModelIndex
        """
        return QtCore.QModelIndex()

    def flags(self, index: QtCore.QModelIndex) -> int:
        """Get flags which influence the display of cells

        :param index: index referencing the cell.
        :type index: QModelIndex
        :return: flags
        :rtype: int
        """
        flags = QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemNeverHasChildren | QtCore.Qt.ItemIsSelectable
        key = self.nameForRow(index.row())
        if self._data[key][1] is None and index.column() == 3:
            # the theta dihedral is not editable for alpha-amino acids
            return flags
        else:
            return QtCore.Qt.ItemIsEditable | flags

    def data(self, index: QtCore.QModelIndex, role: int = None) -> typing.Any:
        """Return data needed for displaying/interacting with the cells

        :param index: index referencing the cell
        :type index: QModelIndex
        :param role: kind of data to be returned
        :type role: int
        :return: display/edit data according to `role`
        :rtype: various
        """
        key = self.nameForRow(index.row())
        if role == QtCore.Qt.DisplayRole:
            if index.column() == 0:
                return key
            elif index.column() == 1:
                return 'alpha' if self._data[key][1] is None else 'beta'
            elif index.column() == 2:
                return '{:.2f} °'.format(self._data[key][0])
            elif index.column() == 3:
                return '{:.2f} °'.format(self._data[key][1]) if self._data[key][1] is not None else '--'
            elif index.column() == 4:
                return '{:.2f} °'.format(self._data[key][2])
            else:
                raise ValueError('Invalid column: {}'.format(index.column()))
        elif role == QtCore.Qt.EditRole:
            if index.column() == 0:
                return key
            elif index.column() == 1:
                return 'alpha' if self._data[key][1] is None else 'beta'
            elif index.column() == 2:
                return self._data[key][0]
            elif index.column() == 3:
                return self._data[key][1]
            elif index.column() == 4:
                return self._data[key][2]
            else:
                raise ValueError('Invalid column: {}'.format(index.column()))
        return None

    def setData(self, index: QtCore.QModelIndex, value: typing.Any, role: int = None) -> bool:
        """Update data in a cell, possibly because of user interaction

        :param index: the index to update
        :type index: QModelIndex
        :param value: the new value
        :type value: various
        :param role: the role to update (typically QtCore.Qt.EditRole)
        :type role: int
        :return: True if updating succeeded
        :rtype: bool
        """
        key = self.nameForRow(index.row())
        if role == QtCore.Qt.EditRole:
            if index.column() == 0:
                # edit the name
                self.beginResetModel()
                ss = self._data.pop(key)
                self._data[value] = ss
                self.endResetModel()
                return True
            elif index.column() == 1:
                # change alpha/beta type
                if value == 'alpha':
                    self._data[key] = (self._data[key][0], None, self._data[key][2])
                elif value == 'beta':
                    self._data[key] = (self._data[key][0], 180., self._data[key][2])
                else:
                    raise ValueError(value)
                self.dataChanged.emit(self.index(index.row(), index.column()),
                                      self.index(index.row(), index.column()))
                return True
            elif index.column() == 2:
                self._data[key] = (value, self._data[key][1], self._data[key][2])
                self.dataChanged.emit(self.index(index.row(), index.column()),
                                      self.index(index.row(), index.column()))
                return True
            elif index.column() == 3:
                self._data[key] = (self._data[key][0], value, self._data[key][2])
                self.dataChanged.emit(self.index(index.row(), index.column()),
                                      self.index(index.row(), index.column()))
                return True
            elif index.column() == 4:
                self._data[key] = (self._data[key][0], self._data[key][1], value)
                self.dataChanged.emit(self.index(index.row(), index.column()),
                                      self.index(index.row(), index.column()))
                return True

    def headerData(self, section: int, orientation: QtCore.Qt.Orientation, role: int = None) -> typing.Any:
        """Get data for displaying the header (only column headers are implemented)

        :param section: column number
        :type section: int
        :param orientation: orientation. Qt.Horizontal is implemented
        :type orientation: int
        :param role: data role. Qt.DisplayRole is implemented only
        :type role: int
        :return: the caption of the column
        :rtype: str
        """
        if role == QtCore.Qt.DisplayRole and orientation == QtCore.Qt.Horizontal:
            return ['Name', 'α/β', 'φ (°)', 'ϑ (°)', 'ψ (°)'][section]
        return None

    def index(self, row: int, column: int, parent: QtCore.QModelIndex = None) -> QtCore.QModelIndex:
        """Create a model index for a cell

        :param row: row index
        :type row: int
        :param column: colum index
        :type column: int
        :param parent: parent (not used, this is a flat model)
        :type parent: QModelIndex or None
        :return: a model index
        :rtype: QModelIndex
        """
        if parent is None:
            parent = QtCore.QModelIndex()
        return self.createIndex(row, column, None) if not parent.isValid() else QtCore.QModelIndex()

    def fetchFromPyMOL(self):
        """Fetch secondary structure database information from PyMOL and update this model accordingly."""
        self.beginResetModel()
        self._data = SecondaryStructureDB.getAll()
        self.endResetModel()

    def applyToPyMOL(self):
        """Push the changes in this model to PyMOL"""
        allnames = list(SecondaryStructureDB.getAll().keys())
        for name in allnames:
            SecondaryStructureDB.remove(name)
        for name in self._data:
            SecondaryStructureDB.add(name, self._data[name][0], self._data[name][1], self._data[name][2])

    def add(self, nameprefix='Untitled') -> QtCore.QModelIndex:
        """Add a new secondary structure

        :param nameprefix: the prefix of the new name (e.g 'Untitled')
        :type nameprefix: str
        :return: the model index of the freshly added entry
        :rtype: QModelIndex
        """
        i = 0
        while nameprefix + str(i) in self._data:
            i += 1
        self.beginResetModel()
        self._data[nameprefix + str(i)] = (0, 0, 0)
        self.endResetModel()
        return self.index(self.rowForName(nameprefix + str(i)), 0, None)

    def remove(self, name: str):
        """Remove the named entry

        :param name: the name of the entry to be removed
        :type name: str
        """
        self.beginResetModel()
        del self._data[name]
        self.endResetModel()

    def addDefaults(self):
        """Update entries to their default values
        """
        self.beginResetModel()
        self._data.update(SecondaryStructureDB.DEFAULT_HELIXTYPES)
        self.endResetModel()

    def duplicate(self, name):
        """Create a duplicate of the named entry

        The new entry will be named "`name`_copy"

        :param name: the name of the entry to be duplicated
        :type name: str
        """
        self.beginResetModel()
        self._data[name + '_copy'] = self._data[name]
        self.endResetModel()
        return self.index(self.rowForName(name + '_copy'), 0, None)

    def nameForRow(self, row: int) -> str:
        """Get the name of the entry on row `row`

        :param row: row index
        :type row: int
        :return: the name of the entry
        :rtype: str
        """
        return list(sorted(self._data.keys()))[row]

    def rowForName(self, name: str) -> int:
        """Get the row index for entry named `name`

        :param name: the name of the entry
        :type name: str
        :return: the row index
        :rtype: int
        """
        return list(sorted(self._data.keys())).index(name)
