# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'src/pmlbeta/recognition/recognition.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName("Form")
        Form.resize(498, 504)
        self.verticalLayout = QtWidgets.QVBoxLayout(Form)
        self.verticalLayout.setObjectName("verticalLayout")
        self.groupBox = QtWidgets.QGroupBox(Form)
        self.groupBox.setObjectName("groupBox")
        self.gridLayout = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout.setObjectName("gridLayout")
        self.browseRTPToolButton = QtWidgets.QToolButton(self.groupBox)
        icon = QtGui.QIcon.fromTheme("document-open")
        self.browseRTPToolButton.setIcon(icon)
        self.browseRTPToolButton.setObjectName("browseRTPToolButton")
        self.gridLayout.addWidget(self.browseRTPToolButton, 2, 2, 1, 1)
        self.label = QtWidgets.QLabel(self.groupBox)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.label_2 = QtWidgets.QLabel(self.groupBox)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 2, 0, 1, 1)
        self.modelComboBox = QtWidgets.QComboBox(self.groupBox)
        self.modelComboBox.setEditable(True)
        self.modelComboBox.setObjectName("modelComboBox")
        self.gridLayout.addWidget(self.modelComboBox, 0, 1, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.groupBox)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 1, 0, 1, 1)
        self.atomCountLabel = QtWidgets.QLabel(self.groupBox)
        self.atomCountLabel.setObjectName("atomCountLabel")
        self.gridLayout.addWidget(self.atomCountLabel, 1, 1, 1, 1)
        self.rtpLineEdit = QtWidgets.QLineEdit(self.groupBox)
        self.rtpLineEdit.setObjectName("rtpLineEdit")
        self.gridLayout.addWidget(self.rtpLineEdit, 2, 1, 1, 1)
        self.reloadModelsToolButton = QtWidgets.QToolButton(self.groupBox)
        icon = QtGui.QIcon.fromTheme("view-refresh")
        self.reloadModelsToolButton.setIcon(icon)
        self.reloadModelsToolButton.setObjectName("reloadModelsToolButton")
        self.gridLayout.addWidget(self.reloadModelsToolButton, 0, 2, 1, 1)
        self.label_4 = QtWidgets.QLabel(self.groupBox)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 3, 0, 1, 1)
        self.rtpCountLabel = QtWidgets.QLabel(self.groupBox)
        self.rtpCountLabel.setObjectName("rtpCountLabel")
        self.gridLayout.addWidget(self.rtpCountLabel, 3, 1, 1, 1)
        self.verticalLayout.addWidget(self.groupBox)
        self.groupBox_2 = QtWidgets.QGroupBox(Form)
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox_2)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.choicesTreeView = QtWidgets.QTreeView(self.groupBox_2)
        self.choicesTreeView.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.choicesTreeView.setSizeAdjustPolicy(QtWidgets.QAbstractScrollArea.AdjustToContents)
        self.choicesTreeView.setAlternatingRowColors(True)
        self.choicesTreeView.setRootIsDecorated(False)
        self.choicesTreeView.setUniformRowHeights(True)
        self.choicesTreeView.setObjectName("choicesTreeView")
        self.verticalLayout_2.addWidget(self.choicesTreeView)
        self.verticalLayout.addWidget(self.groupBox_2)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.findMatchingPushButton = QtWidgets.QPushButton(Form)
        self.findMatchingPushButton.setObjectName("findMatchingPushButton")
        self.horizontalLayout.addWidget(self.findMatchingPushButton)
        self.AssignPushButton = QtWidgets.QPushButton(Form)
        self.AssignPushButton.setObjectName("AssignPushButton")
        self.horizontalLayout.addWidget(self.AssignPushButton)
        self.verticalLayout.addLayout(self.horizontalLayout)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        _translate = QtCore.QCoreApplication.translate
        Form.setWindowTitle(_translate("Form", "Molecule recognition"))
        self.groupBox.setTitle(_translate("Form", "Input"))
        self.browseRTPToolButton.setText(_translate("Form", "Browse"))
        self.label.setText(_translate("Form", "Model / selection:"))
        self.label_2.setText(_translate("Form", "RTP file:"))
        self.label_3.setText(_translate("Form", "Atom count:"))
        self.atomCountLabel.setText(_translate("Form", "--"))
        self.reloadModelsToolButton.setText(_translate("Form", "Reload"))
        self.label_4.setText(_translate("Form", "Residue topologies:"))
        self.rtpCountLabel.setText(_translate("Form", "--"))
        self.groupBox_2.setTitle(_translate("Form", "Choices"))
        self.findMatchingPushButton.setText(_translate("Form", "Find matching"))
        self.AssignPushButton.setText(_translate("Form", "Assign"))