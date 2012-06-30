#!/usr/bin/env python

# Copyright (C) 2012 Victor Semionov
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#  * Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#  * Neither the name of the copyright holder nor the names of the contributors
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


import sys
import os
import re

import multiprocessing

from PySide import QtCore, QtGui

import meta
import pamp
import npamp
import unitconv

import mainwin
import outwin


defaults = npamp.core.copy_conf(npamp.params.__dict__)


def repr_classless(v):
    return re.sub("<class '([^']+?)'>", "\\1", repr(v))

def file2conf(path):
    conf = defaults.copy()
    execfile(path, conf)
    conf = npamp.core.copy_conf(conf)
    return conf

def conf2file(conf, path):
    with open(path, "w") as fp:
        fp.write("import pamp\n")
        for param, value in sorted(conf.items()):
            fp.write("%s = %s\n" % (param, repr_classless(value)))


class AppWindow(QtGui.QMainWindow, mainwin.Ui_MainWindow):
    
    untitled_name = "untitled.%s" % meta.file_extension
    home_dir = QtCore.QDir.toNativeSeparators(QtCore.QDir.homePath())
    file_filter = "%s Files (*.%s)" % (meta.app_name, meta.file_extension)
    
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.old_excepthook = None
        self.widget_module_map = dict()
        self.working_conf = npamp.core.copy_conf(defaults)
        self.working_path = None
        self.output_path = None
        self.setupUi(self)
        self.initWidgets()
        self.connectSignals()
        self.updateUI()
        self.conf2gui(defaults)
        self.resize(self.sizeHint())
    
    def add_excepthook(self):
        self.old_excepthook = sys.excepthook
        sys.excepthook = self.gui_excepthook
    
    def gui_excepthook(self, exc_type, value, traceback):
        if self.old_excepthook:
            self.old_excepthook(exc_type, value, traceback)
        QtGui.QMessageBox.critical(self, "Error", str(value))
    
    def initWidgets(self):
        for name in dir(pamp.beam):
            obj = getattr(pamp.beam, name)
            if type(obj) is type and issubclass(obj, pamp.beam.BeamProfile) and obj is not pamp.beam.BeamProfile:
                self.comboBox_beam_class.addItem(name)
        self.widget_module_map[self.comboBox_beam_class.objectName()] = pamp.beam
        for name in dir(pamp.pulse):
            obj = getattr(pamp.pulse, name)
            if type(obj) is type and issubclass(obj, pamp.pulse.SinglePulse) and obj is not pamp.pulse.SinglePulse:
                self.comboBox_pulse_class.addItem(name)
        self.widget_module_map[self.comboBox_pulse_class.objectName()] = pamp.pulse
        for name in dir(pamp.loss):
            obj = getattr(pamp.loss, name)
            if type(obj) is type and issubclass(obj, pamp.loss.LossModel) and hasattr(obj, "concrete") and obj.concrete is True:
                self.comboBox_loss_model_class.addItem(name)
                self.comboBox_loss_dependence_alt_model.addItem(name)
                self.listWidget_compared_loss_models.addItem(name)
        self.widget_module_map[self.comboBox_loss_model_class.objectName()] = pamp.loss
        self.widget_module_map[self.comboBox_loss_dependence_alt_model.objectName()] = pamp.loss
        self.widget_module_map[self.listWidget_compared_loss_models.objectName()] = pamp.loss
        for name in dir(pamp.inverter):
            obj = getattr(pamp.inverter, name)
            if type(obj) is type and issubclass(obj, pamp.inverter.PopulationInverter) and obj is not pamp.inverter.PopulationInverter:
                self.comboBox_inverter_class.addItem(name)
        self.widget_module_map[self.comboBox_inverter_class.objectName()] = pamp.inverter
        for name in dir(pamp.integral):
            obj = getattr(pamp.integral, name)
            if type(obj) is type and issubclass(obj, pamp.integral.SampleIntegrator) and obj is not pamp.integral.SampleIntegrator:
                self.listWidget_integrator_classes.addItem(name)
        self.widget_module_map[self.listWidget_integrator_classes.objectName()] = pamp.integral
        self.listWidget_integrator_classes.resize(self.listWidget_integrator_classes.sizeHint())
        for name in dir(pamp.amplifier):
            obj = getattr(pamp.amplifier, name)
            if type(obj) is type and issubclass(obj, pamp.amplifier.NumericalAmplifier) and obj is not pamp.amplifier.NumericalAmplifier:
                self.listWidget_amplifier_classes.addItem(name)
        self.widget_module_map[self.listWidget_amplifier_classes.objectName()] = pamp.amplifier
        self.listWidget_amplifier_classes.resize(self.listWidget_amplifier_classes.sizeHint())
    
    def shortFilename(self):
        return os.path.splitext(os.path.basename(self.working_path or self.untitled_name))[0]
    
    def updateUI(self):
        filename = self.shortFilename()
        self.setWindowTitle("%s - %s" % (meta.gui_app_name, filename))
        self.actionClose.setEnabled(self.working_path is not None)
    
    def conf2gui(self, conf):
        def set_widget_value(label, widget, value):
            if type(widget) is QtGui.QLineEdit:
                if type(value) is float:
                    value = unitconv.convert_to_input(label.text(), value)
                widget.setText(repr_classless(value))
            elif type(widget) is QtGui.QSpinBox:
                widget.setValue(value)
            elif type(widget) is QtGui.QCheckBox:
                widget.setChecked(value)
            elif type(widget) is QtGui.QComboBox:
                widget.setCurrentIndex(widget.findText(value.__name__))
            elif type(widget) is QtGui.QListWidget:
                widget.setCurrentItem(None, QtGui.QItemSelectionModel.Clear)
                for v in value:
                    item = widget.findItems(v.__name__, QtCore.Qt.MatchExactly|QtCore.Qt.MatchCaseSensitive)[0]
                    widget.setCurrentItem(item, QtGui.QItemSelectionModel.Select)
            elif type(widget) is QtGui.QWidget:
                children = [w for w in widget.findChildren(QtGui.QWidget) if w.parent() is widget]
                assert len(children) == len(value), "wrong data length"
                children.sort(key = lambda item: item.objectName())
                for w, v in zip(children, value):
                    set_widget_value(label, w, v)
            else:
                assert False, "unhandled widget type"
        for parameter, value in conf.items():
            children = self.centralWidget().findChildren(QtGui.QWidget, QtCore.QRegExp("^[^_]+_%s$" % parameter))
            widgets = [w for w in children if type(w) is not QtGui.QLabel]
            labels = [w for w in children if type(w) is QtGui.QLabel]
            assert len(widgets) == 1, "none or more than one widget matches parameter name \"%s\"" % parameter
            assert len(labels) == 1, "none or more than one label matches parameter name \"%s\"" % parameter
            set_widget_value(labels[0], widgets[0], value)
        
        if self.gui2conf() != conf:
            raise npamp.core.ConfigurationError("invalid parameter value(s)")
    
    def gui2conf(self):
        def get_widget_value(label, widget, defval):
            def item2class(w, t):
                module = self.widget_module_map[w.objectName()]
                cls = getattr(module, t)
                return cls
            if type(widget) is QtGui.QLineEdit:
                if type(defval) is float:
                    value = float(widget.text())
                    value = unitconv.convert_from_input(label.text(), value)
                    return value
                else:
                    return eval(widget.text())
            elif type(widget) is QtGui.QSpinBox:
                return widget.value()
            elif type(widget) is QtGui.QCheckBox:
                return widget.isChecked()
            elif type(widget) is QtGui.QComboBox:
                return item2class(widget, widget.currentText())
            elif type(widget) is QtGui.QListWidget:
                clsnames = map(lambda item: item.text(), widget.selectedItems())
                classes = map(lambda clsname: item2class(widget, clsname), clsnames)
                return list(classes)
            elif type(widget) is QtGui.QWidget:
                children = [w for w in widget.findChildren(QtGui.QWidget) if w.parent() is widget]
                assert len(children) == len(defval), "wrong data length"
                children.sort(key = lambda item: item.objectName())
                return type(defval)(map(lambda wv: get_widget_value(label, *wv), zip(children, defval)))
            else:
                assert False, "unhandled widget type"
        conf = dict()
        for parameter, default in defaults.items():
            children = self.centralWidget().findChildren(QtGui.QWidget, QtCore.QRegExp("^[^_]+_%s$" % parameter))
            widgets = [w for w in children if type(w) is not QtGui.QLabel]
            labels = [w for w in children if type(w) is QtGui.QLabel]
            assert len(widgets) == 1, "none or more than one widget matches parameter name \"%s\"" % parameter
            assert len(labels) == 1, "none or more than one label matches parameter name \"%s\"" % parameter
            value = get_widget_value(labels[0], widgets[0], default)
            conf[parameter] = value
        return conf
    
    def connectSignals(self):
        self.actionNew.triggered.connect(self.onNew)
        self.actionOpen.triggered.connect(self.onOpen)
        self.actionSave.triggered.connect(self.onSave)
        self.actionSaveAs.triggered.connect(self.onSaveAs)
        self.actionClose.triggered.connect(self.onCloseFile)
        self.actionExecute.triggered.connect(self.onExecute)
        self.actionQuit.triggered.connect(self.onQuit)
        self.actionAbout.triggered.connect(self.onAbout)
    
    def checkSave(self):
        if self.gui2conf() != self.working_conf:
            choice = QtGui.QMessageBox.question(self, "Changes Made", "Do you want to save changes to %s?" % self.shortFilename(), QtGui.QMessageBox.Save|QtGui.QMessageBox.Discard|QtGui.QMessageBox.Cancel)
            if choice is QtGui.QMessageBox.Save:
                return True
            elif choice is QtGui.QMessageBox.Discard:
                return False
            else:
                return None
        else:
            return False
    
    def confirmProceed(self):
        save = self.checkSave()
        if save is True:
            self.save()
            return True
        elif save is False:
            return True
        else:
            return False
    
    def save(self):
        if self.working_path is not None:
            conf = self.gui2conf()
            conf2file(conf, self.working_path)
            self.working_conf = conf
            return True
        else:
            return self.saveAs()
    
    def saveAs(self):
        default = self.working_path or os.path.join(self.home_dir, self.untitled_name)
        path, _ = QtGui.QFileDialog.getSaveFileName(self, "Save As", default, self.file_filter)
        if not path:
            return False
        conf = self.gui2conf()
        conf2file(conf, path)
        self.working_conf = conf
        self.working_path = path
        self.output_path = None
        self.updateUI()
        return True
    
    def onNew(self):
        if self.confirmProceed():
            self.conf2gui(defaults)
            self.working_conf = npamp.core.copy_conf(defaults)
            self.working_path = None
            self.output_path = None
            self.updateUI()
    
    def openFile(self, path):
        conf = file2conf(path)
        old_conf = self.gui2conf()
        try:
            self.conf2gui(conf)
        except Exception:
            self.conf2gui(old_conf)
            raise
        self.working_conf = conf
        self.working_path = path
        self.output_path = None
        self.updateUI()
    
    def onOpen(self):
        if self.confirmProceed():
            default = os.path.dirname(self.working_path) if self.working_path else self.home_dir
            path, _ = QtGui.QFileDialog.getOpenFileName(self, "Open", default, self.file_filter)
            if not path:
                return
            self.openFile(path)
    
    def onSave(self):
        self.save()
    
    def onSaveAs(self):
        self.saveAs()
    
    def onCloseFile(self):
        if self.working_path is not None:
            if self.confirmProceed():
                self.working_path = None
                self.output_path = None
        self.updateUI()
    
    def onExecute(self):
        conf = self.gui2conf()
        
        if conf["graphs"]:
            if self.output_path is not None:
                output_path = self.output_path
            else:
                default = os.path.dirname(self.working_path) if self.working_path else self.home_dir
                output_path = QtGui.QFileDialog.getExistingDirectory(self, "Output Directory", default)
                if not output_path:
                    return
                if os.listdir(output_path):
                    choice = QtGui.QMessageBox.warning(self, "Directory Not Empty", "The selected directory is not empty. Do you want to continue?", QtGui.QMessageBox.Yes|QtGui.QMessageBox.No)
                    if choice is not QtGui.QMessageBox.Yes:
                        return
                self.output_path = output_path
        else:
            output_path = None
        
        in_conn, out_conn = multiprocessing.Pipe(False)
        proc = multiprocessing.Process(name=meta.app_name, target=worker, args=(out_conn, conf, output_path))
        
        thr = InputThread(in_conn)
        out = OutputWindow(self, proc, thr)
        thr.received.connect(out.onOutputRead)
        thr.finished.connect(out.onWorkerFinished)
        proc.start()
        out_conn.close()
        thr.start()
        out.exec_()
    
    def onQuit(self):
        self.close()
    
    def onAbout(self):
        text = "%s %s\n%s\n\n%s\n\n%s\n%s\n\n%s" % (meta.app_name, meta.app_version, meta.app_copyright, meta.app_description, meta.app_author_msg, meta.app_coauthors_msg, meta.app_website_msg)
        QtGui.QMessageBox.about(self, "About %s" % meta.app_name, text)
    
    def closeEvent(self, event):
        if self.confirmProceed():
            event.accept()
        else:
            event.ignore()

class OutputWindow(QtGui.QDialog, outwin.Ui_Dialog):
    
    def __init__(self, parent, proc, thr):
        QtGui.QDialog.__init__(self, parent)
        self.proc = proc
        self.thr = thr
        self.setupUi(self)
        self.connectSignals()
        self.setWindowTitle(meta.app_name + " output")
        width = max(self.sizeHint().width(), self.width())
        height = max(self.sizeHint().height(), self.height())
        self.resize(width, height)
        self.setCursor(QtCore.Qt.BusyCursor)
    
    def connectSignals(self):
        self.buttonBox.clicked.connect(self.onButtonClicked)
    
    def onButtonClicked(self, button):
        buttons = self.buttonBox.buttons()
        if button is buttons[0]:
            self.stopWorker()
        elif button is buttons[1]:
            self.close()
        else:
            assert False, "unhandled button"
    
    def reject(self):
        self.close()
    
    def stopWorker(self):
        self.proc.terminate()
        self.proc.join()
        self.thr.wait()
        self.buttonBox.buttons()[0].setEnabled(False)
    
    def onOutputRead(self, s):
        self.plainTextEdit.appendPlainText(s)
        self.plainTextEdit.ensureCursorVisible()
    
    def onWorkerFinished(self):
        self.stopWorker()
        self.unsetCursor()
    
    def closeEvent(self, event):
        self.stopWorker()
        event.accept()

class InputThread(QtCore.QThread):
    
    received = QtCore.Signal(str)
    finished = QtCore.Signal()
    
    def __init__(self, in_conn):
        super(InputThread, self).__init__()
        self.in_conn = in_conn
    
    def run(self):
        while True:
            try:
                s = self.in_conn.recv()
            except EOFError:
                break
            self.received.emit(s)
        self.in_conn.close()
        self.finished.emit()

def worker(out_conn, conf, output_path):
    mpout = npamp.mp.MPOutput(out_conn.send)
    sys.stdout = mpout
    sys.stderr = mpout
    npamp.params.__dict__.update(conf)
    npamp.run(None, output_path)


def main():
    multiprocessing.freeze_support()
    app = QtGui.QApplication(sys.argv)
    win = AppWindow()
    win.show()
    win.add_excepthook()
    
    args = sys.argv[1:]
    if len(args) > 1:
        raise npamp.InvocationError("too many arguments")
    if args:
        win.openFile(args[0])
    
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
