# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'src/pmlbeta/betafabgui2/dihedraleditor.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(471, 322)
        self.verticalLayout = QtWidgets.QVBoxLayout(Dialog)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.addToolButton = QtWidgets.QToolButton(Dialog)
        icon = QtGui.QIcon.fromTheme("list-add")
        self.addToolButton.setIcon(icon)
        self.addToolButton.setObjectName("addToolButton")
        self.horizontalLayout.addWidget(self.addToolButton)
        self.removeToolButton = QtWidgets.QToolButton(Dialog)
        icon = QtGui.QIcon.fromTheme("list-remove")
        self.removeToolButton.setIcon(icon)
        self.removeToolButton.setObjectName("removeToolButton")
        self.horizontalLayout.addWidget(self.removeToolButton)
        self.duplicateToolButton = QtWidgets.QToolButton(Dialog)
        icon = QtGui.QIcon.fromTheme("edit-copy")
        self.duplicateToolButton.setIcon(icon)
        self.duplicateToolButton.setObjectName("duplicateToolButton")
        self.horizontalLayout.addWidget(self.duplicateToolButton)
        self.addDefaultsToolButton = QtWidgets.QToolButton(Dialog)
        icon = QtGui.QIcon.fromTheme("view-refresh")
        self.addDefaultsToolButton.setIcon(icon)
        self.addDefaultsToolButton.setObjectName("addDefaultsToolButton")
        self.horizontalLayout.addWidget(self.addDefaultsToolButton)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.treeView = QtWidgets.QTreeView(Dialog)
        self.treeView.setObjectName("treeView")
        self.verticalLayout.addWidget(self.treeView)
        self.buttonBox = QtWidgets.QDialogButtonBox(Dialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Apply|QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(Dialog)
        self.buttonBox.accepted.connect(Dialog.accept)
        self.buttonBox.rejected.connect(Dialog.reject)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        self.addToolButton.setToolTip(_translate("Dialog", "Add a new entry"))
        self.addToolButton.setText(_translate("Dialog", "Add"))
        self.removeToolButton.setToolTip(_translate("Dialog", "Remove the selected entry"))
        self.removeToolButton.setText(_translate("Dialog", "Remove"))
        self.duplicateToolButton.setToolTip(_translate("Dialog", "Duplicate the selected entry"))
        self.duplicateToolButton.setText(_translate("Dialog", "Duplicate"))
        self.addDefaultsToolButton.setToolTip(_translate("Dialog", "Add defaults"))
        self.addDefaultsToolButton.setText(_translate("Dialog", "Add/update defaults"))

