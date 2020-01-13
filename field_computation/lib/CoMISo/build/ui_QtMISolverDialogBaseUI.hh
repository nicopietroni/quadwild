/********************************************************************************
** Form generated from reading UI file 'QtMISolverDialogBaseUI.ui'
**
** Created: Wed Aug 20 13:00:01 2014
**      by: Qt User Interface Compiler version 4.8.4
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_QTMISOLVERDIALOGBASEUI_H
#define UI_QTMISOLVERDIALOGBASEUI_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QDialog>
#include <QtGui/QDoubleSpinBox>
#include <QtGui/QGridLayout>
#include <QtGui/QHBoxLayout>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QVBoxLayout>

QT_BEGIN_NAMESPACE

class Ui_QtMISolverDialogBaseUI
{
public:
    QVBoxLayout *verticalLayout;
    QGridLayout *gridLayout;
    QCheckBox *cplexRoundingCB;
    QLabel *label_4;
    QLabel *label_6;
    QLabel *label_5;
    QDoubleSpinBox *localErrorDSB;
    QSpinBox *localItersSB;
    QLabel *label_3;
    QDoubleSpinBox *cgErrorDSB;
    QSpinBox *cgItersSB;
    QLabel *label_7;
    QHBoxLayout *horizontalLayout;
    QCheckBox *initialFullCB;
    QCheckBox *iterFullCB;
    QCheckBox *finalFullCB;
    QLabel *label_2;
    QHBoxLayout *horizontalLayout_3;
    QCheckBox *directRoundingCB;
    QCheckBox *noRoundingCB;
    QCheckBox *multipleRoundingCB;
    QDoubleSpinBox *multipleRoundingDSB;
    QLabel *label;
    QHBoxLayout *horizontalLayout_4;
    QSpinBox *infoSB;
    QCheckBox *solverStatsCheckBox;
    QLabel *label_8;
    QCheckBox *use_reordering_cb;
    QLabel *label_9;
    QCheckBox *gurobiRoundingCB;
    QDoubleSpinBox *gurobiMaxTimeDSB;
    QLabel *label_10;
    QHBoxLayout *horizontalLayout_2;
    QPushButton *okPB;
    QPushButton *cancelPB;
    QSpacerItem *verticalSpacer;

    void setupUi(QDialog *QtMISolverDialogBaseUI)
    {
        if (QtMISolverDialogBaseUI->objectName().isEmpty())
            QtMISolverDialogBaseUI->setObjectName(QString::fromUtf8("QtMISolverDialogBaseUI"));
        QtMISolverDialogBaseUI->setWindowModality(Qt::ApplicationModal);
        QtMISolverDialogBaseUI->resize(385, 273);
        verticalLayout = new QVBoxLayout(QtMISolverDialogBaseUI);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        gridLayout = new QGridLayout();
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        cplexRoundingCB = new QCheckBox(QtMISolverDialogBaseUI);
        cplexRoundingCB->setObjectName(QString::fromUtf8("cplexRoundingCB"));

        gridLayout->addWidget(cplexRoundingCB, 8, 1, 1, 1);

        label_4 = new QLabel(QtMISolverDialogBaseUI);
        label_4->setObjectName(QString::fromUtf8("label_4"));

        gridLayout->addWidget(label_4, 0, 1, 1, 1);

        label_6 = new QLabel(QtMISolverDialogBaseUI);
        label_6->setObjectName(QString::fromUtf8("label_6"));

        gridLayout->addWidget(label_6, 0, 2, 1, 1);

        label_5 = new QLabel(QtMISolverDialogBaseUI);
        label_5->setObjectName(QString::fromUtf8("label_5"));
        QSizePolicy sizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(label_5->sizePolicy().hasHeightForWidth());
        label_5->setSizePolicy(sizePolicy);

        gridLayout->addWidget(label_5, 1, 0, 1, 1);

        localErrorDSB = new QDoubleSpinBox(QtMISolverDialogBaseUI);
        localErrorDSB->setObjectName(QString::fromUtf8("localErrorDSB"));
        localErrorDSB->setDecimals(2);
        localErrorDSB->setMinimum(-1e+09);
        localErrorDSB->setMaximum(1e+09);
        localErrorDSB->setValue(6);

        gridLayout->addWidget(localErrorDSB, 1, 1, 1, 1);

        localItersSB = new QSpinBox(QtMISolverDialogBaseUI);
        localItersSB->setObjectName(QString::fromUtf8("localItersSB"));
        localItersSB->setMaximum(1000000000);
        localItersSB->setValue(10000);

        gridLayout->addWidget(localItersSB, 1, 2, 1, 1);

        label_3 = new QLabel(QtMISolverDialogBaseUI);
        label_3->setObjectName(QString::fromUtf8("label_3"));
        sizePolicy.setHeightForWidth(label_3->sizePolicy().hasHeightForWidth());
        label_3->setSizePolicy(sizePolicy);

        gridLayout->addWidget(label_3, 2, 0, 1, 1);

        cgErrorDSB = new QDoubleSpinBox(QtMISolverDialogBaseUI);
        cgErrorDSB->setObjectName(QString::fromUtf8("cgErrorDSB"));
        cgErrorDSB->setDecimals(2);
        cgErrorDSB->setMinimum(-1e+09);
        cgErrorDSB->setMaximum(1e+09);
        cgErrorDSB->setValue(6);

        gridLayout->addWidget(cgErrorDSB, 2, 1, 1, 1);

        cgItersSB = new QSpinBox(QtMISolverDialogBaseUI);
        cgItersSB->setObjectName(QString::fromUtf8("cgItersSB"));
        cgItersSB->setMaximum(1000000000);
        cgItersSB->setValue(20);

        gridLayout->addWidget(cgItersSB, 2, 2, 1, 1);

        label_7 = new QLabel(QtMISolverDialogBaseUI);
        label_7->setObjectName(QString::fromUtf8("label_7"));
        sizePolicy.setHeightForWidth(label_7->sizePolicy().hasHeightForWidth());
        label_7->setSizePolicy(sizePolicy);

        gridLayout->addWidget(label_7, 3, 0, 1, 1);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        initialFullCB = new QCheckBox(QtMISolverDialogBaseUI);
        initialFullCB->setObjectName(QString::fromUtf8("initialFullCB"));
        initialFullCB->setChecked(true);

        horizontalLayout->addWidget(initialFullCB);

        iterFullCB = new QCheckBox(QtMISolverDialogBaseUI);
        iterFullCB->setObjectName(QString::fromUtf8("iterFullCB"));
        iterFullCB->setChecked(true);

        horizontalLayout->addWidget(iterFullCB);

        finalFullCB = new QCheckBox(QtMISolverDialogBaseUI);
        finalFullCB->setObjectName(QString::fromUtf8("finalFullCB"));
        finalFullCB->setChecked(true);
        finalFullCB->setTristate(false);

        horizontalLayout->addWidget(finalFullCB);


        gridLayout->addLayout(horizontalLayout, 3, 1, 1, 2);

        label_2 = new QLabel(QtMISolverDialogBaseUI);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout->addWidget(label_2, 4, 0, 1, 1);

        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        directRoundingCB = new QCheckBox(QtMISolverDialogBaseUI);
        directRoundingCB->setObjectName(QString::fromUtf8("directRoundingCB"));

        horizontalLayout_3->addWidget(directRoundingCB);

        noRoundingCB = new QCheckBox(QtMISolverDialogBaseUI);
        noRoundingCB->setObjectName(QString::fromUtf8("noRoundingCB"));

        horizontalLayout_3->addWidget(noRoundingCB);

        multipleRoundingCB = new QCheckBox(QtMISolverDialogBaseUI);
        multipleRoundingCB->setObjectName(QString::fromUtf8("multipleRoundingCB"));

        horizontalLayout_3->addWidget(multipleRoundingCB);

        multipleRoundingDSB = new QDoubleSpinBox(QtMISolverDialogBaseUI);
        multipleRoundingDSB->setObjectName(QString::fromUtf8("multipleRoundingDSB"));
        multipleRoundingDSB->setDecimals(2);
        multipleRoundingDSB->setMinimum(0);
        multipleRoundingDSB->setMaximum(0.5);
        multipleRoundingDSB->setValue(0.5);

        horizontalLayout_3->addWidget(multipleRoundingDSB);


        gridLayout->addLayout(horizontalLayout_3, 4, 1, 1, 2);

        label = new QLabel(QtMISolverDialogBaseUI);
        label->setObjectName(QString::fromUtf8("label"));
        sizePolicy.setHeightForWidth(label->sizePolicy().hasHeightForWidth());
        label->setSizePolicy(sizePolicy);

        gridLayout->addWidget(label, 5, 0, 1, 1);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        infoSB = new QSpinBox(QtMISolverDialogBaseUI);
        infoSB->setObjectName(QString::fromUtf8("infoSB"));

        horizontalLayout_4->addWidget(infoSB);

        solverStatsCheckBox = new QCheckBox(QtMISolverDialogBaseUI);
        solverStatsCheckBox->setObjectName(QString::fromUtf8("solverStatsCheckBox"));

        horizontalLayout_4->addWidget(solverStatsCheckBox);


        gridLayout->addLayout(horizontalLayout_4, 5, 1, 1, 2);

        label_8 = new QLabel(QtMISolverDialogBaseUI);
        label_8->setObjectName(QString::fromUtf8("label_8"));

        gridLayout->addWidget(label_8, 6, 0, 1, 1);

        use_reordering_cb = new QCheckBox(QtMISolverDialogBaseUI);
        use_reordering_cb->setObjectName(QString::fromUtf8("use_reordering_cb"));
        use_reordering_cb->setChecked(true);

        gridLayout->addWidget(use_reordering_cb, 6, 1, 1, 1);

        label_9 = new QLabel(QtMISolverDialogBaseUI);
        label_9->setObjectName(QString::fromUtf8("label_9"));

        gridLayout->addWidget(label_9, 7, 0, 1, 1);

        gurobiRoundingCB = new QCheckBox(QtMISolverDialogBaseUI);
        gurobiRoundingCB->setObjectName(QString::fromUtf8("gurobiRoundingCB"));

        gridLayout->addWidget(gurobiRoundingCB, 7, 1, 1, 1);

        gurobiMaxTimeDSB = new QDoubleSpinBox(QtMISolverDialogBaseUI);
        gurobiMaxTimeDSB->setObjectName(QString::fromUtf8("gurobiMaxTimeDSB"));
        gurobiMaxTimeDSB->setMaximum(1e+09);
        gurobiMaxTimeDSB->setValue(60);

        gridLayout->addWidget(gurobiMaxTimeDSB, 7, 2, 1, 1);

        label_10 = new QLabel(QtMISolverDialogBaseUI);
        label_10->setObjectName(QString::fromUtf8("label_10"));

        gridLayout->addWidget(label_10, 7, 3, 1, 1);


        verticalLayout->addLayout(gridLayout);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        okPB = new QPushButton(QtMISolverDialogBaseUI);
        okPB->setObjectName(QString::fromUtf8("okPB"));

        horizontalLayout_2->addWidget(okPB);

        cancelPB = new QPushButton(QtMISolverDialogBaseUI);
        cancelPB->setObjectName(QString::fromUtf8("cancelPB"));

        horizontalLayout_2->addWidget(cancelPB);


        verticalLayout->addLayout(horizontalLayout_2);

        verticalSpacer = new QSpacerItem(20, 14, QSizePolicy::Minimum, QSizePolicy::Expanding);

        verticalLayout->addItem(verticalSpacer);


        retranslateUi(QtMISolverDialogBaseUI);

        QMetaObject::connectSlotsByName(QtMISolverDialogBaseUI);
    } // setupUi

    void retranslateUi(QDialog *QtMISolverDialogBaseUI)
    {
        QtMISolverDialogBaseUI->setWindowTitle(QApplication::translate("QtMISolverDialogBaseUI", "Dialog", 0, QApplication::UnicodeUTF8));
        cplexRoundingCB->setText(QApplication::translate("QtMISolverDialogBaseUI", "use cplex", 0, QApplication::UnicodeUTF8));
        label_4->setText(QApplication::translate("QtMISolverDialogBaseUI", "Max Error (1e-x)", 0, QApplication::UnicodeUTF8));
        label_6->setText(QApplication::translate("QtMISolverDialogBaseUI", "Iterations", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("QtMISolverDialogBaseUI", "Local", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("QtMISolverDialogBaseUI", "CG", 0, QApplication::UnicodeUTF8));
        label_7->setText(QApplication::translate("QtMISolverDialogBaseUI", "Full", 0, QApplication::UnicodeUTF8));
        initialFullCB->setText(QApplication::translate("QtMISolverDialogBaseUI", "initial", 0, QApplication::UnicodeUTF8));
        iterFullCB->setText(QApplication::translate("QtMISolverDialogBaseUI", "iter not converged", 0, QApplication::UnicodeUTF8));
        finalFullCB->setText(QApplication::translate("QtMISolverDialogBaseUI", "final", 0, QApplication::UnicodeUTF8));
        label_2->setText(QApplication::translate("QtMISolverDialogBaseUI", "Rounding", 0, QApplication::UnicodeUTF8));
        directRoundingCB->setText(QApplication::translate("QtMISolverDialogBaseUI", "direct", 0, QApplication::UnicodeUTF8));
        noRoundingCB->setText(QApplication::translate("QtMISolverDialogBaseUI", "no", 0, QApplication::UnicodeUTF8));
        multipleRoundingCB->setText(QApplication::translate("QtMISolverDialogBaseUI", "multiple", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("QtMISolverDialogBaseUI", "Info Level", 0, QApplication::UnicodeUTF8));
        solverStatsCheckBox->setText(QApplication::translate("QtMISolverDialogBaseUI", "Output solver statistics", 0, QApplication::UnicodeUTF8));
        label_8->setText(QApplication::translate("QtMISolverDialogBaseUI", "Constraints", 0, QApplication::UnicodeUTF8));
        use_reordering_cb->setText(QApplication::translate("QtMISolverDialogBaseUI", "reordering", 0, QApplication::UnicodeUTF8));
        label_9->setText(QApplication::translate("QtMISolverDialogBaseUI", "GUROBI", 0, QApplication::UnicodeUTF8));
        gurobiRoundingCB->setText(QApplication::translate("QtMISolverDialogBaseUI", "use gurobi", 0, QApplication::UnicodeUTF8));
        label_10->setText(QApplication::translate("QtMISolverDialogBaseUI", "s", 0, QApplication::UnicodeUTF8));
        okPB->setText(QApplication::translate("QtMISolverDialogBaseUI", "Ok", 0, QApplication::UnicodeUTF8));
        cancelPB->setText(QApplication::translate("QtMISolverDialogBaseUI", "Cancel", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class QtMISolverDialogBaseUI: public Ui_QtMISolverDialogBaseUI {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_QTMISOLVERDIALOGBASEUI_H
