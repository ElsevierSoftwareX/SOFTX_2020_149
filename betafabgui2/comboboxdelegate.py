import re

from PyQt5 import QtWidgets, QtCore

from ..betafab2 import BetaPeptide
from ..secstructdb import SecondaryStructureDB


def isAlphaResidue(resname: str) -> bool:
    return re.match(r'^\([RSLD]\)A[A-Z]{1,2}', resname) is not None


def isBetaResidue(resname: str) -> bool:
    return re.match(
        r'^(\(2[RS]3[RS]\)(ACPC|ACHC))|(\([RS]\)B[23]h[A-Z]{1,2})|BA|\(2[RS]3[RS]\)B23h\(2[A-Z]{1,2}3[A-Z]{1,2}\)',
        resname) is not None


class ComboBoxDelegate(QtWidgets.QStyledItemDelegate):
    def createEditor(self, parent: QtWidgets.QWidget, option: QtWidgets.QStyleOptionViewItem,
                     index: QtCore.QModelIndex) -> QtWidgets.QComboBox:
        """Create an editor widget

        :param parent: the parent widget
        :type parent: QtWidgets.QWidget
        :param option: style options
        :type option: QtWidgets.QStyleOptionViewItem
        :param index: the model index
        :type index: QtCore.QModelIndex
        :return: the combo box
        :rtype: QtWidgets.QComboBox
        """
        combobox = QtWidgets.QComboBox(parent)
        if self._kind == 'chirality':
            model = index.model()
            row = index.row()
            resname = str(model.residue(row))
            if isBetaResidue(resname):
                combobox.addItems(['R', 'S'])
            if isAlphaResidue(resname):
                combobox.addItems(['D', 'L', 'R', 'S'])
        elif self._kind == 'sidechain':
            combobox.addItems(sorted(BetaPeptide.SIDECHAINS))
        elif self._kind == 'kind':
            combobox.addItems(['A', 'BA', 'B2', 'B3', 'B23', 'ACPC', 'ACHC'])
        elif self._kind == 'ss':
            model = index.model()
            row = index.row()
            resname = str(model.residue(row))
            if isBetaResidue(resname):
                combobox.addItems(sorted(SecondaryStructureDB.getAll(alpha=False, beta=True)))
            if isAlphaResidue(resname):
                combobox.addItems(sorted(SecondaryStructureDB.getAll(alpha=True, beta=False)))
        elif self._kind == 'alpha/beta':
            combobox.addItems(['alpha', 'beta'])
        return combobox

    def setEditorData(self, editor: QtWidgets.QComboBox, index: QtCore.QModelIndex):
        """Update the editor state from the model.

        :param editor: the editor widget
        :type editor: QtWidgets.QComboBox
        :param index: the model index
        :type index: QtCore.QModelIndex
        """
        data = index.data(QtCore.Qt.EditRole)
        idx = editor.findText(data)
        if idx < 0:
            idx = 0
        editor.setCurrentIndex(idx)

    def setModelData(self, editor: QtWidgets.QComboBox, model: QtCore.QAbstractItemModel, index: QtCore.QModelIndex):
        """Update the model from the editor.

        :param editor: the editor widget
        :type editor: QtWidgets.QComboBox
        :param model: the model
        :type model: SequenceModel
        :param index: the model index
        :type index: QtCore.QModelIndex
        """
        model.setData(index, editor.currentText(), QtCore.Qt.EditRole)

    def updateEditorGeometry(self, editor: QtWidgets.QComboBox, option: QtWidgets.QStyleOptionViewItem,
                             index: QtCore.QModelIndex):
        """Update the editor geometry

        :param editor: the editor widget
        :type editor: QtWidgets.QComboBox
        :param option: options
        :type option: QtWidgets.QStyleOptionViewItem
        :param index: the model index
        :type index: QtCore.QModelIndex
        """
        editor.setGeometry(option.rect)

    def __init__(self, parent: QtWidgets.QWidget, kind: str):
        """Create a new combo box delegate

        :param parent: the parent widget
        :type parent: QtWidgets.QWidget
        :param kind: the content type of this delegate. Can be 'chirality' or 'sidechain'
        :type kind: str
        """
        self._kind = kind
        super().__init__(parent)
