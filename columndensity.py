# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'columndensity.ui'
#
# Created by: PyQt5 UI code generator 5.12.3
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_ColumnDensityDialog(object):
    def setupUi(self, ColumnDensityDialog):
        ColumnDensityDialog.setObjectName("ColumnDensityDialog")
        ColumnDensityDialog.resize(288, 221)
        self.verticalLayoutWidget = QtWidgets.QWidget(ColumnDensityDialog)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(9, 9, 270, 201))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label.setWordWrap(True)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.columnDensity_edit = QtWidgets.QLineEdit(self.verticalLayoutWidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.columnDensity_edit.sizePolicy().hasHeightForWidth())
        self.columnDensity_edit.setSizePolicy(sizePolicy)
        self.columnDensity_edit.setObjectName("columnDensity_edit")
        self.horizontalLayout.addWidget(self.columnDensity_edit)
        self.label_2 = QtWidgets.QLabel(self.verticalLayoutWidget)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout.addWidget(self.label_2)
        self.exponent_rollbox = QtWidgets.QSpinBox(self.verticalLayoutWidget)
        self.exponent_rollbox.setMaximum(20)
        self.exponent_rollbox.setObjectName("exponent_rollbox")
        self.horizontalLayout.addWidget(self.exponent_rollbox)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.columnDensityAccept_button = QtWidgets.QPushButton(self.verticalLayoutWidget)
        self.columnDensityAccept_button.setObjectName("columnDensityAccept_button")
        self.verticalLayout.addWidget(self.columnDensityAccept_button)

        self.retranslateUi(ColumnDensityDialog)
        QtCore.QMetaObject.connectSlotsByName(ColumnDensityDialog)

    def retranslateUi(self, ColumnDensityDialog):
        _translate = QtCore.QCoreApplication.translate
        ColumnDensityDialog.setWindowTitle(_translate("ColumnDensityDialog", "Initial Column Density"))
        self.label.setText(_translate("ColumnDensityDialog", "Enter an estimate of the column density."))
        self.label_2.setText(_translate("ColumnDensityDialog", "x10^"))
        self.columnDensityAccept_button.setText(_translate("ColumnDensityDialog", "Accept"))
