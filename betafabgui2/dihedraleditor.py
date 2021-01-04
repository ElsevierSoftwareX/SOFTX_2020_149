from PyQt5 import QtCore, QtWidgets

from .comboboxdelegate import ComboBoxDelegate
from .dihedraleditor_ui import Ui_Dialog
from .dihedralmodel import DihedralModel
from .spinboxdelegate import DoubleSpinBoxDelegate


class DihedralEditor(QtWidgets.QDialog, Ui_Dialog):
    """Editor widget for the secondary structure database"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)

    def setupUi(self, Form):
        """Initialize the user interface"""
        super().setupUi(Form)

        # add callbacks
        self.addDefaultsToolButton.clicked.connect(self.addDefaults)
        self.removeToolButton.clicked.connect(self.removeSecondaryStructure)
        self.addToolButton.clicked.connect(self.addSecondaryStructure)
        self.duplicateToolButton.clicked.connect(self.duplicateSecondaryStructure)

        self.buttonBox.button(self.buttonBox.Apply).clicked.connect(self.pushToPyMOL)
        self.buttonBox.button(self.buttonBox.Ok).clicked.connect(self.pushToPyMOL)

        self._delegates = {'type': ComboBoxDelegate(self.treeView, 'alpha/beta'),
                           'dihedral': DoubleSpinBoxDelegate(self.treeView, -180, 180, 1, ' °', True)}
        self.model = DihedralModel()
        self.treeView.setModel(self.model)
        self.treeView.setItemDelegateForColumn(1, self._delegates['type'])
        self.treeView.setItemDelegateForColumn(2, self._delegates['dihedral'])
        self.treeView.setItemDelegateForColumn(3, self._delegates['dihedral'])
        self.treeView.setItemDelegateForColumn(4, self._delegates['dihedral'])
        self.setWindowTitle('Dihedral editor for secondary structures')
        self.model.fetchFromPyMOL()
        self.model.dataChanged.connect(self.resizeColumns)
        self.model.modelReset.connect(self.resizeColumns)
        self.resizeColumns()

    def resizeColumns(self):
        """Resize the columns of the treeview. This is called whenever the underlying model changes."""
        for c in range(self.model.columnCount()):
            self.treeView.resizeColumnToContents(c)

    def addDefaults(self):
        """Add back defaults if they were deleted somehow
        """
        self.model.addDefaults()

    def removeSecondaryStructure(self):
        """Remove the selected secondary structure and keep the selection on the same line."""
        try:
            row = self.treeView.selectionModel().selectedRows(0)[0].row()
            self.model.remove(self.model.nameForRow(row))
            if row > self.model.rowCount() - 1:
                row = self.model.rowCount() - 1
            if row >= 0:
                self.treeView.selectionModel().select(
                    self.model.index(row, 0),
                    QtCore.QItemSelectionModel.Rows | QtCore.QItemSelectionModel.ClearAndSelect)
        except IndexError:
            # no selection
            return

    def addSecondaryStructure(self):
        """Add a new secondary structure and select it"""
        idx = self.model.add()
        self.treeView.selectionModel().select(
            idx, QtCore.QItemSelectionModel.Rows | QtCore.QItemSelectionModel.ClearAndSelect)

    def duplicateSecondaryStructure(self):
        """Duplicate the selected secondary structure and select the duplicate"""
        try:
            row = self.treeView.selectionModel().selectedRows(0)[0].row()
            idx = self.model.duplicate(self.model.nameForRow(row))
            self.treeView.selectionModel().select(
                idx,
                QtCore.QItemSelectionModel.Rows | QtCore.QItemSelectionModel.ClearAndSelect)
        except IndexError:
            # no selection
            return

    def pushToPyMOL(self):
        """Apply changes in PyMOL"""
        self.model.applyToPyMOL()


diheditwidget = None


def run():
    global diheditwidget
    if not diheditwidget:
        diheditwidget = DihedralEditor()
    diheditwidget.show()
    diheditwidget.raise_()
