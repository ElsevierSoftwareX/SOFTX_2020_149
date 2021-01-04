import warnings
from typing import List

from PyQt5 import QtCore, QtGui, QtWidgets

from .betafab2_ui import Ui_Form
from .sequencemodel import SequenceModel
from ..betafab2 import BetaPeptide, betafab2
from ..utils import select_bbb

try:
    from pymol import cmd
except ImportError:
    warnings.warn(
        'Cannot import PyMOL: functionality will suffer (you can ignore this if you are just building the documentation).')
from .comboboxdelegate import ComboBoxDelegate
from ..secstructdb import SecondaryStructureDB
from ..savegro import save_g96, save_crd


class Betafab2GUI(QtWidgets.QWidget, Ui_Form):
    """Graphical user interface for betafab2"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setupUi(self)

    def setupUi(self, Form):
        """Initialize the user interface"""
        super().setupUi(Form)

        # initialize callbacks
        self.addFromCommandLinePushButton.clicked.connect(self.addFromCommandLine)
        self.addFromComboBoxesPushButton.clicked.connect(self.addFromComboBoxes)
        self.AcePushButton.clicked.connect(self.addAce)
        self.ButPushButton.clicked.connect(self.addBut)
        self.NMePushButton.clicked.connect(self.addNMe)
        self.removeToolButton.clicked.connect(self.removeSelectedResidues)
        self.clearToolButton.clicked.connect(self.removeAllResidues)
        self.moveToTopToolButton.clicked.connect(self.moveResidueToTop)
        self.moveUpToolButton.clicked.connect(self.moveResidueUp)
        self.moveDownToolButton.clicked.connect(self.moveResidueDown)
        self.moveToBottomToolButton.clicked.connect(self.moveResidueToBottom)
        self.generatePushButton.clicked.connect(self.build)
        self.saveCRDPushButton.clicked.connect(self.saveCRD)
        self.saveg96PushButton.clicked.connect(self.saveG96)
        self.aminoAcidKindComboBox.currentIndexChanged.connect(self.aminoAcidKindChanged)
        self.commandLineLineEdit.returnPressed.connect(self.addFromCommandLine)
        self.targetObjectNameLineEdit.textChanged.connect(self.targetObjectNameChanged)
        self.targetObjectNameChanged()  # adjust the enabled/disabled state of the "Build" button
        self.copyToolButton.clicked.connect(self.copyResidues)
        self.pasteToolButton.clicked.connect(self.pasteResidues)
        self.cutToolButton.clicked.connect(self.cutResidues)

        # initialize the entry comboboxes
        self.aminoAcidKindChanged()
        self.sidechain2ComboBox.addItems(sorted(list(BetaPeptide.SIDECHAINS.keys())))
        self.sidechain3ComboBox.addItems(sorted(list(BetaPeptide.SIDECHAINS.keys())))
        self.sidechain2ComboBox.setCurrentIndex(0)
        self.sidechain3ComboBox.setCurrentIndex(0)

        # initialize the treeview
        self.model = SequenceModel()
        self.sequenceTreeView.setModel(self.model)
        self.resizeTreeViewColumns()
        self.sequenceTreeView.selectionModel().selectionChanged.connect(self.residueSelectionChanged)
        self.model.dataChanged.connect(self.modelChanged)
        self.model.modelReset.connect(self.modelChanged)
        self.model.rowsInserted.connect(self.modelChanged)
        self.model.rowsRemoved.connect(self.modelChanged)
        self.model.validationChanged.connect(self.modelValidationChanged)
        self._delegates = {
            'chirality': ComboBoxDelegate(self.sequenceTreeView, 'chirality'),
            'sidechain': ComboBoxDelegate(self.sequenceTreeView, 'sidechain'),
            'kind': ComboBoxDelegate(self.sequenceTreeView, 'kind'),
            'ss': ComboBoxDelegate(self.sequenceTreeView, 'ss')
        }
        self.sequenceTreeView.setItemDelegateForColumn(2, self._delegates['kind'])
        self.sequenceTreeView.setItemDelegateForColumn(3, self._delegates['chirality'])
        self.sequenceTreeView.setItemDelegateForColumn(4, self._delegates['sidechain'])
        self.sequenceTreeView.setItemDelegateForColumn(5, self._delegates['chirality'])
        self.sequenceTreeView.setItemDelegateForColumn(6, self._delegates['sidechain'])
        self.sequenceTreeView.setItemDelegateForColumn(7, self._delegates['ss'])

        # testing
        self.VALXVALPushButton.clicked.connect(self.addVALXVAL)

    def modelChanged(self):
        """Called whenever the sequence model changed"""
        self.resizeTreeViewColumns()

    def modelValidationChanged(self, valid: bool):
        """Called whenever the sequence model has been validated: either successfully or not

        :param valid: if the validation succeeded
        :type valid: bool
        """
        self.generatePushButton.setEnabled(valid and bool(self.targetObjectNameLineEdit.text().strip()))

    def resizeTreeViewColumns(self):
        """Resize the columns of the sequence treeview to the content"""
        for c in range(self.model.columnCount()):
            self.sequenceTreeView.resizeColumnToContents(c)

    def addVALXVAL(self):
        """Create a VALXVAL peptide, mostly for testing purposes"""
        self.model.append('ACE', '', '', '', '', None)
        self.model.append('B3', '', '', 'S', 'V', 'H14M')
        self.model.append('B3', '', '', 'S', 'A', 'H14M')
        self.model.append('B3', '', '', 'S', 'L', 'H14M')
        self.model.append('B23', 'S', 'A', 'S', 'A', 'H14M')
        self.model.append('B3', '', '', 'S', 'V', 'H14M')
        self.model.append('B3', '', '', 'S', 'A', 'H14M')
        self.model.append('B3', '', '', 'S', 'L', 'H14M')
        self.model.append('NME', '', '', '', '', None)

    def addFromCommandLine(self):
        """Add a new residue from the "command-line"."""
        try:
            sequence_parsed = BetaPeptide.parseBetaPeptideSequence(self.commandLineLineEdit.text())
            for residue in sequence_parsed:
                self.model.append(
                    residue['kind'], residue['stereo2'], residue['sidechain2'],
                    residue['stereo3'], residue['sidechain3'], residue['dihedrals'])
            self.commandLineLineEdit.clear()
        except ValueError as ve:
            QtWidgets.QMessageBox.critical(self, 'Error while parsing sequence', ve.args[0], QtWidgets.QMessageBox.Ok)

    def addFromComboBoxes(self):
        """Add a new residue from the combo boxes."""
        kinds = ['A', 'BA', 'B2', 'B3', 'B23', 'ACPC', 'ACHC']
        try:
            dihedrals = SecondaryStructureDB.dihedrals(self.secondaryStructureComboBox.currentText())
        except KeyError:
            raise ValueError('Unknown secondary structure: {}'.format(self.secondaryStructureComboBox.currentText()))
        dihedrals = [d for d in dihedrals if d is not None]
        self.model.append(
            kinds[self.aminoAcidKindComboBox.currentIndex()],
            self.stereo2ComboBox.currentText(), self.sidechain2ComboBox.currentText(),
            self.stereo3ComboBox.currentText(), self.sidechain3ComboBox.currentText(), dihedrals)

    def addAce(self):
        """Add an acetyl group at the N terminus"""
        self.model.append('ACE', '', '', '', '', None)

    def addBut(self):
        """Add a butyl group at the N terminus"""
        self.model.append('BUT', '', '', '', '', None)

    def addNMe(self):
        """Add an N-methylamide group at the C terminus"""
        self.model.append('NME', '', '', '', '', None)

    def updateSecondaryStructureEntryComboBox(self, alpha: bool, beta: bool):
        """Update the list of secondary structures in the corresponding entry combo box

        :param alpha: if alpha-peptidic secondary structures are requested
        :type alpha: bool
        :param beta: if beta-peptidic secondary structures are requested
        :type beta: bool
        """
        ssenabled = SecondaryStructureDB.getAll(alpha, beta)
        current = self.secondaryStructureComboBox.currentText()
        self.secondaryStructureComboBox.blockSignals(True)
        self.secondaryStructureComboBox.clear()
        self.secondaryStructureComboBox.addItems(ssenabled)
        self.secondaryStructureComboBox.blockSignals(False)
        index = self.secondaryStructureComboBox.findText(current)
        if index < 0: index = 0
        self.secondaryStructureComboBox.setCurrentIndex(index)

    def aminoAcidKindChanged(self):
        """The amino acid king is changed"""
        # the entry comboboxes must be enabled/disabled accordingly.
        # also the secondary structure combo box must be changed

        self.stereo2ComboBox.clear()
        self.stereo2ComboBox.addItems(
            ['D', 'L', 'R', 'S'] if self.aminoAcidKindComboBox.currentIndex() == 0 else ['R', 'S'])

        if self.aminoAcidKindComboBox.currentIndex() == 0:
            # alpha-amino acid
            self.stereo2ComboBox.setEnabled(True)
            if self.stereo2ComboBox.currentIndex() < 0:
                self.stereo2ComboBox.setCurrentIndex(0)
            self.stereo3ComboBox.setCurrentIndex(-1)
            self.stereo3ComboBox.setEnabled(False)
            self.sidechain2ComboBox.setEnabled(True)
            if self.sidechain2ComboBox.currentIndex() < 0:
                self.sidechain2ComboBox.setCurrentIndex(0)
            self.sidechain3ComboBox.setEnabled(False)
            self.sidechain3ComboBox.setCurrentIndex(-1)
            self.updateSecondaryStructureEntryComboBox(alpha=True, beta=False)
        elif self.aminoAcidKindComboBox.currentIndex() == 1:
            # beta-homo-glycine
            self.stereo2ComboBox.setEnabled(False)
            self.stereo2ComboBox.setCurrentIndex(-1)
            self.stereo3ComboBox.setEnabled(False)
            self.stereo3ComboBox.setCurrentIndex(-1)
            self.sidechain2ComboBox.setEnabled(False)
            self.sidechain2ComboBox.setCurrentIndex(-1)
            self.sidechain3ComboBox.setEnabled(False)
            self.sidechain3ComboBox.setCurrentIndex(-1)
            self.updateSecondaryStructureEntryComboBox(alpha=False, beta=True)
        elif self.aminoAcidKindComboBox.currentIndex() == 2:
            # beta2-amino acid
            self.stereo2ComboBox.setEnabled(True)
            if self.stereo2ComboBox.currentIndex() < 0:
                self.stereo2ComboBox.setCurrentIndex(0)
            self.stereo3ComboBox.setEnabled(False)
            self.stereo3ComboBox.setCurrentIndex(-1)
            self.sidechain2ComboBox.setEnabled(True)
            if self.sidechain2ComboBox.currentIndex() < 0:
                self.sidechain2ComboBox.setCurrentIndex(0)
            self.sidechain3ComboBox.setEnabled(False)
            self.sidechain3ComboBox.setCurrentIndex(-1)
            self.updateSecondaryStructureEntryComboBox(alpha=False, beta=True)
        elif self.aminoAcidKindComboBox.currentIndex() == 3:
            # beta3-amino acid
            self.stereo2ComboBox.setEnabled(False)
            self.stereo2ComboBox.setCurrentIndex(-1)
            self.stereo3ComboBox.setEnabled(True)
            if self.stereo3ComboBox.currentIndex() < 0:
                self.stereo3ComboBox.setCurrentIndex(0)
            self.sidechain2ComboBox.setEnabled(False)
            self.sidechain2ComboBox.setCurrentIndex(-1)
            self.sidechain3ComboBox.setEnabled(True)
            if self.sidechain3ComboBox.currentIndex() < 0:
                self.sidechain3ComboBox.setCurrentIndex(0)
            self.updateSecondaryStructureEntryComboBox(alpha=False, beta=True)
        elif self.aminoAcidKindComboBox.currentIndex() == 4:
            # beta23-amino acid
            self.stereo2ComboBox.setEnabled(True)
            if self.stereo2ComboBox.currentIndex() < 0:
                self.stereo2ComboBox.setCurrentIndex(0)
            self.stereo3ComboBox.setEnabled(True)
            if self.stereo3ComboBox.currentIndex() < 0:
                self.stereo3ComboBox.setCurrentIndex(0)
            self.sidechain2ComboBox.setEnabled(True)
            if self.sidechain2ComboBox.currentIndex() < 0:
                self.sidechain2ComboBox.setCurrentIndex(0)
            self.sidechain3ComboBox.setEnabled(True)
            if self.sidechain3ComboBox.currentIndex() < 0:
                self.sidechain3ComboBox.setCurrentIndex(0)
            self.updateSecondaryStructureEntryComboBox(alpha=False, beta=True)
        elif self.aminoAcidKindComboBox.currentIndex() in [5, 6]:
            # ACPC and ACHC
            self.stereo2ComboBox.setEnabled(True)
            if self.stereo2ComboBox.currentIndex() < 0:
                self.stereo2ComboBox.setCurrentIndex(0)
            self.stereo3ComboBox.setEnabled(True)
            if self.stereo3ComboBox.currentIndex() < 0:
                self.stereo3ComboBox.setCurrentIndex(0)
            self.sidechain2ComboBox.setEnabled(False)
            self.sidechain2ComboBox.setCurrentIndex(-1)
            self.sidechain3ComboBox.setEnabled(False)
            self.sidechain3ComboBox.setCurrentIndex(-1)
            self.updateSecondaryStructureEntryComboBox(alpha=False, beta=True)

    def removeSelectedResidues(self):
        """Remove the currently selected residues"""
        selectedindexes = self.sequenceTreeView.selectionModel().selectedIndexes()
        while selectedindexes:
            self.model.removeRow(selectedindexes[0].row())
            selectedindexes = self.sequenceTreeView.selectionModel().selectedIndexes()

    def removeAllResidues(self):
        """Remove all residues"""
        self.model.clear()

    def moveResidueToTop(self):
        """Move the selected residue to the top of the list"""
        try:
            row = self.sequenceTreeView.selectionModel().selectedIndexes()[0].row()
        except IndexError:
            return
        rowafter = self.model.moveToTop(row)
        self.sequenceTreeView.selectionModel().select(
            self.model.index(rowafter, 0),
            QtCore.QItemSelectionModel.ClearAndSelect | QtCore.QItemSelectionModel.Rows)

    def moveResidueUp(self):
        """Move the selected residue up in the list"""
        try:
            row = self.sequenceTreeView.selectionModel().selectedIndexes()[0].row()
        except IndexError:
            return
        rowafter = self.model.moveUp(row)
        self.sequenceTreeView.selectionModel().select(
            self.model.index(rowafter, 0),
            QtCore.QItemSelectionModel.ClearAndSelect | QtCore.QItemSelectionModel.Rows)

    def moveResidueDown(self):
        """Move the selected residue down in the list"""
        try:
            row = self.sequenceTreeView.selectionModel().selectedIndexes()[0].row()
        except IndexError:
            return
        rowafter = self.model.moveDown(row)
        self.sequenceTreeView.selectionModel().select(
            self.model.index(rowafter, 0),
            QtCore.QItemSelectionModel.ClearAndSelect | QtCore.QItemSelectionModel.Rows)

    def moveResidueToBottom(self):
        """Move the selected residue to the bottom"""
        try:
            row = self.sequenceTreeView.selectionModel().selectedIndexes()[0].row()
        except IndexError:
            return
        rowafter = self.model.moveToBottom(row)
        self.sequenceTreeView.selectionModel().select(
            self.model.index(rowafter, 0),
            QtCore.QItemSelectionModel.ClearAndSelect | QtCore.QItemSelectionModel.Rows)

    def build(self):
        """Build the peptide"""
        name = self.targetObjectNameLineEdit.text()
        if name in cmd.get_object_list():
            if QtWidgets.QMessageBox.question(
                    self, 'Overwrite existing model?',
                    'Model {} already exists. Do you want to overwrite it?'.format(name),
                    QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                    QtWidgets.QMessageBox.No) == QtWidgets.QMessageBox.NoRole:
                # "No" was selected by the user, return without building
                return
            cmd.delete('model {}'.format(name))
        residues = self.model.residues()
        try:
            betafab2(name, *residues)
        except Exception as exc:
            QtWidgets.QMessageBox.critical(
                self, 'Error while building peptide',
                'Error while building peptide:\n{}'.format(exc.args[0]))
        if self.cleanupCheckBox.isChecked():
            # try to energy minimize the structure
            select_bbb('_bbone_of_{}'.format(name), 'model {}'.format(name))
            cmd.flag('fix', '_bbone_of_{}'.format(name), 'set')
            try:
                import freemol.mengine
                cmd.clean('model {}'.format(name))
            except:
                warnings.warn('Cannot run clean(): This PyMOL build appears not to include full modeling capabilities '
                              '(could not import the freemol.mengine module). Using the inferior sculpting facility.')
                cmd.sculpt_activate(name)
                cmd.sculpt_iterate(name, cycles=1000)
                cmd.sculpt_deactivate(name)
            cmd.flag('fix', '_bbone_of_{}'.format(name), 'clear')
            cmd.delete('_bbone_of_{}'.format(name))

    def saveCRD(self):
        """Save the peptide to a CHARMM CRD file"""
        filename, fltr = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save to CHARMM CRD...', '', 'CHARMM CRD files (*.crd);;All files (*)', 'CHARMM CRD files (*.crd)')
        if not filename:
            return
        if not filename.endswith('.crd'):
            filename = filename + '.crd'
        save_crd(filename, 'model {}'.format(self.targetObjectNameLineEdit.text()))

    def saveG96(self):
        """Save the peptide to a GROMOS96 file"""
        filename, fltr = QtWidgets.QFileDialog.getSaveFileName(
            self, 'Save to GROMOS96...', '', 'GROMOS96 coordinate files (*.g96);;All files (*)',
            'GROMOS96 coordinate files (*.g96)')
        if not filename:
            return
        if not filename.endswith('.g96'):
            filename = filename + '.g96'
        save_g96(filename, 'model {}'.format(self.targetObjectNameLineEdit.text()))

    def residueSelectionChanged(self, selected: List[QtCore.QModelIndex], deselected: List[QtCore.QModelIndex]):
        """The selection in the treeview is changed. Enable/disable tool buttons accordingly"""
        selectedindexes = self.sequenceTreeView.selectionModel().selectedIndexes()
        selectedrows = {i.row() for i in selectedindexes}
        for movebutton in [self.moveToTopToolButton, self.moveUpToolButton, self.moveDownToolButton,
                           self.moveToBottomToolButton]:
            movebutton.setEnabled(len(selectedrows) == 1)
        self.removeToolButton.setEnabled(len(selectedrows) > 0)

    def targetObjectNameChanged(self):
        """The target object name field is changed. If it is empty, disable the "Build" button."""
        self.generatePushButton.setEnabled(bool(self.targetObjectNameLineEdit.text().strip()) and self.model.isValid())

    def cutResidues(self):
        """Cut the selected residues to the clipboard"""
        self.copyResidues()
        self.removeSelectedResidues()

    def copyResidues(self):
        """Copy the selected residues to the clipboard"""
        QtGui.QGuiApplication.instance().clipboard().setText(
            ', '.join(
                [str(self.model.residue(sel.row())) for sel in self.sequenceTreeView.selectionModel().selectedRows(0)]),
        )

    def pasteResidues(self):
        """Paste residues from the clipboard"""
        try:
            selrow = self.sequenceTreeView.selectionModel().selectedRows(0)[0].row()
        except IndexError:
            # no row is selected
            return
        try:
            sequence_parsed = BetaPeptide.parseBetaPeptideSequence(QtGui.QGuiApplication.instance().clipboard().text())
        except ValueError as ve:
            QtWidgets.QMessageBox.critical(self, 'Error while parsing sequence', ve.args[0], QtWidgets.QMessageBox.Ok)
            return
        for destrow, residue in enumerate(sequence_parsed, start=selrow):
            self.model.insert(destrow,
                              residue['kind'], residue['stereo2'], residue['sidechain2'],
                              residue['stereo3'], residue['sidechain3'], residue['dihedrals'])


betafabwidget = None


def run():
    global betafabwidget
    if not betafabwidget:
        betafabwidget = Betafab2GUI()
    betafabwidget.show()
    betafabwidget.raise_()
