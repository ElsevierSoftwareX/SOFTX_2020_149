from PyQt5 import QtWidgets, QtCore


class DoubleSpinBoxDelegate(QtWidgets.QStyledItemDelegate):
    def createEditor(self, parent: QtWidgets.QWidget, option: QtWidgets.QStyleOptionViewItem,
                     index: QtCore.QModelIndex) -> QtWidgets.QDoubleSpinBox:
        """Create an editor widget

        :param parent: the parent widget
        :type parent: QtWidgets.QWidget
        :param option: style options
        :type option: QtWidgets.QStyleOptionViewItem
        :param index: the model index
        :type index: QtCore.QModelIndex
        :return: the spin box
        :rtype: QtWidgets.QDoubleSpinBox
        """
        spinbox = QtWidgets.QDoubleSpinBox(parent)
        spinbox.setRange(self._minimum, self._maximum)
        spinbox.setSingleStep(self._step)
        if self._wrap is not None:
            spinbox.setWrapping(self._wrap)
        if self._suffix is not None:
            spinbox.setSuffix(self._suffix)
        return spinbox

    def setEditorData(self, editor: QtWidgets.QDoubleSpinBox, index: QtCore.QModelIndex):
        """Update the editor state from the model.

        :param editor: the editor widget
        :type editor: QtWidgets.QDoubleSpinBox
        :param index: the model index
        :type index: QtCore.QModelIndex
        """
        data = index.data(QtCore.Qt.EditRole)
        editor.setValue(data)

    def setModelData(self, editor: QtWidgets.QDoubleSpinBox, model: QtCore.QAbstractItemModel,
                     index: QtCore.QModelIndex):
        """Update the model from the editor.

        :param editor: the editor widget
        :type editor: QtWidgets.QDoubleSpinBox
        :param model: the model
        :type model: SequenceModel
        :param index: the model index
        :type index: QtCore.QModelIndex
        """
        model.setData(index, editor.value(), QtCore.Qt.EditRole)

    def updateEditorGeometry(self, editor: QtWidgets.QDoubleSpinBox, option: QtWidgets.QStyleOptionViewItem,
                             index: QtCore.QModelIndex):
        """Update the editor geometry

        :param editor: the editor widget
        :type editor: QtWidgets.QDoubleSpinBox
        :param option: options
        :type option: QtWidgets.QStyleOptionViewItem
        :param index: the model index
        :type index: QtCore.QModelIndex
        """
        editor.setGeometry(option.rect)

    def __init__(self, parent: QtWidgets.QWidget, minimum: float, maximum: float, step: float, suffix: str = None,
                 wrap: bool = None):
        """Create a new spin box delegate

        :param parent: the parent widget
        :type parent: QtWidgets.QWidget
        :param minimum: the lowest value in the spin box
        :type minimum: float
        :param maximum: the highets value in the spin box
        :type maximum: float
        :param step: step value of the spin box
        :type step: float
        :param suffix: what to display after the number in the spin box
        :type suffix: str
        :param wrap: if the spin box should wrap around
        :type wrap: bool
        """
        self._minimum = minimum
        self._maximum = maximum
        self._step = step
        self._suffix = suffix
        self._wrap = wrap
        super().__init__(parent)
