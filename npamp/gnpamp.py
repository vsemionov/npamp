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

import signal
import threading
import multiprocessing

import getopt
import re

from PySide import QtCore, QtGui

import model

import meta
import npamp
import unitconv

import mainwin
import outwin


old_excepthook = None

debug_mode = False


predef_attr_reprs = [("inf", "float(\"inf\")"), ("nan", "float(\"nan\")")]
predef_attrs = [(name, eval(value)) for (name, value) in predef_attr_reprs]

defaults = npamp.cfg.load_conf(npamp.params.__dict__, None)


def boot_excepthook(exc_type, value, traceback):
    if old_excepthook:
        old_excepthook(exc_type, value, traceback)
    QtGui.QMessageBox.critical(None, "%s Error" % meta.app_name, str(value))


def init_app_dir():
    try:
        os.mkdir(npamp.app_dir)
    except OSError:
        pass


def repr_cfg(v):
    return re.sub("<class '([^']+?)'>", "\\1", repr(v))

def file2conf(path):
    conf = npamp.cfg.load_conf(defaults, path)
    return conf

def conf2file(conf, path):
    with open(path, "w") as fp:
        fp.write("import model\n")
        fp.write("\n")
        fp.write("%s = %s\n" % ("version", npamp.params.version))
        fp.write("\n")
        for name, value in predef_attr_reprs:
            fp.write("%s = %s\n" % (name, value))
        fp.write("\n")
        for param, value in sorted(conf.items()):
            fp.write("%s = %s\n" % (param, repr_cfg(value)))
        fp.write("\n")
        for name, _ in predef_attr_reprs:
            fp.write("del %s\n" % name)


class AppWindow(QtGui.QMainWindow, mainwin.Ui_MainWindow):
    
    untitled_name = "untitled.%s" % meta.file_extension
    file_filter = "%s Files (*.%s)" % (meta.app_name, meta.file_extension)
    
    def __init__(self, extensions):
        QtGui.QMainWindow.__init__(self)
        self.extensions = extensions
        self.monitor_pipe = multiprocessing.Pipe(False)
        self.old_excepthook = None
        self.widget_module_map = dict()
        self.working_conf = npamp.cfg.copy_conf(defaults)
        self.working_path = None
        self.output_path = None
        self.default_directory = self.get_default_directory()
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
    
    def get_default_directory(self):
        if os.name == "nt":
            import winshell
            return winshell.desktop()
        else:
            home_dir = QtCore.QDir.toNativeSeparators(QtCore.QDir.homePath())
            return home_dir
    
    def initWidgets(self):
        for name in dir(model.beam):
            obj = getattr(model.beam, name)
            if type(obj) is type and issubclass(obj, model.beam.BeamProfile) and obj is not model.beam.BeamProfile:
                self.comboBox_beam_class.addItem(name)
        self.widget_module_map[self.comboBox_beam_class.objectName()] = model.beam
        
        for name in dir(model.pulse):
            obj = getattr(model.pulse, name)
            if type(obj) is type and issubclass(obj, model.pulse.SinglePulse) and obj is not model.pulse.SinglePulse:
                self.comboBox_pulse_class.addItem(name)
        self.widget_module_map[self.comboBox_pulse_class.objectName()] = model.pulse
        
        for name in dir(model.depop):
            obj = getattr(model.depop, name)
            if type(obj) is type and issubclass(obj, model.depop.DepopulationModel) and hasattr(obj, "concrete") and obj.concrete is True:
                self.comboBox_depop_model_class.addItem(name)
                self.comboBox_ext_alt_depop_model.addItem(name)
                self.listWidget_ext_depop_models.addItem(name)
        self.widget_module_map[self.comboBox_depop_model_class.objectName()] = model.depop
        self.widget_module_map[self.comboBox_ext_alt_depop_model.objectName()] = model.depop
        self.widget_module_map[self.listWidget_ext_depop_models.objectName()] = model.depop
        
        for name in dir(model.inverter):
            obj = getattr(model.inverter, name)
            if type(obj) is type and issubclass(obj, model.inverter.PopulationInverter) and obj is not model.inverter.PopulationInverter:
                self.comboBox_inverter_class.addItem(name)
        self.widget_module_map[self.comboBox_inverter_class.objectName()] = model.inverter
        
        for name in dir(model.integrator):
            obj = getattr(model.integrator, name)
            if type(obj) is type and issubclass(obj, model.integrator.NumericalIntegrator) and obj is not model.integrator.NumericalIntegrator:
                self.listWidget_integrator_classes.addItem(name)
        self.widget_module_map[self.listWidget_integrator_classes.objectName()] = model.integrator
        self.listWidget_integrator_classes.resize(self.listWidget_integrator_classes.sizeHint())
        
        for name in dir(model.amplifier):
            obj = getattr(model.amplifier, name)
            if type(obj) is type and issubclass(obj, model.amplifier.NumericalAmplifier) and obj is not model.amplifier.NumericalAmplifier:
                self.listWidget_amplifier_classes.addItem(name)
        self.widget_module_map[self.listWidget_amplifier_classes.objectName()] = model.amplifier
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
                widget.setText(repr_cfg(value))
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
            raise npamp.cfg.ConfigurationError("invalid parameter value(s)")
    
    def gui2conf(self):
        def get_widget_value(label, widget, defval):
            def item2class(w, t):
                module = self.widget_module_map[w.objectName()]
                cls = getattr(module, t)
                return cls
            
            if type(widget) is QtGui.QLineEdit:
                if type(defval) is float:
                    value = widget.text()
                    value = unitconv.convert_from_input(label.text(), value)
                    return value
                else:
                    glb = dict(predef_attrs)
                    glb["model"] = npamp.model
                    return eval(widget.text(), glb)
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
    
    def openFile(self, path):
        conf = file2conf(path)
        old_conf = self.gui2conf()
        try:
            self.conf2gui(conf)
        except:
            self.conf2gui(old_conf)
            raise
        self.working_conf = conf
        self.working_path = path
        self.output_path = None
        self.updateUI()
    
    def saveAs(self):
        default = self.working_path or os.path.join(self.default_directory, self.untitled_name)
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
    
    def save(self):
        if self.working_path is not None:
            conf = self.gui2conf()
            conf2file(conf, self.working_path)
            self.working_conf = conf
            return True
        else:
            return self.saveAs()
    
    def onNew(self):
        if self.confirmProceed():
            self.conf2gui(defaults)
            self.working_conf = npamp.cfg.copy_conf(defaults)
            self.working_path = None
            self.output_path = None
            self.updateUI()
    
    def onOpen(self):
        if self.confirmProceed():
            default = os.path.dirname(self.working_path) if self.working_path else self.default_directory
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
                self.working_conf = self.gui2conf()
                self.working_path = None
                self.output_path = None
        self.updateUI()
    
    def onExecute(self):
        conf = self.gui2conf()
        
        if conf["graphs"]:
            if self.output_path is not None:
                output_path = self.output_path
            else:
                default = os.path.dirname(self.working_path) if self.working_path else self.default_directory
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
        proc = multiprocessing.Process(name=meta.app_name, target=worker, args=(self.monitor_pipe, out_conn, conf, output_path, debug_mode))
        
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
        text = "%s %s\n\n%s\n%s\n\n%s\n\n%s\n%s\n\n%s" % (meta.app_name, meta.app_version, meta.app_copyright, meta.app_rights, meta.app_description, meta.app_author_msg, meta.app_coauthors_msg, meta.app_website_msg)
        
        text += "\n\n"
        if not self.extensions:
            text += "No extensions installed."
        else:
            text += "Installed extensions (name: description):"
            for extension in self.extensions:
                name, doc = extension.__name__, extension.__doc__
                text += "\n%s: %s" % (name, doc)
        
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
        self.buttonBox.buttons()[0].setDefault(True)
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

def worker((monitor_in, monitor_out), out_conn, conf, output_path, debug_mode):
    mpout = npamp.mp.MPOutput(out_conn.send)
    sys.stdout = sys.stderr = mpout
    
    monitor_out.close()
    thr = threading.Thread(target=npamp.mp.monitor_thread, args=(monitor_in,))
    thr.daemon = True
    thr.start()
    
    npamp.debug_mode = debug_mode
    npamp.params.__dict__.update(conf)
    npamp.run(None, output_path, None)


def main():
    multiprocessing.freeze_support()
    
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    
    app = QtGui.QApplication(sys.argv)
    
    global old_excepthook
    old_excepthook = sys.excepthook
    sys.excepthook = boot_excepthook
    
    init_app_dir()
    
    opts, args = getopt.getopt(sys.argv[1:], "g")
    for opt, _ in opts:
        if opt == "-g":
            global debug_mode
            debug_mode = True
            npamp.debug_mode = True
        else:
            assert False, "unhandled option"
    
    if len(args) > 1:
        raise npamp.InvocationError("too many arguments")
    
    extensions = npamp.load_extensions()
    
    win = AppWindow(extensions)
    
    if args:
        win.openFile(args[0])
    
    win.show()
    
    sys.excepthook = old_excepthook
    win.add_excepthook()
    
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
