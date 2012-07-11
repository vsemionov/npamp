
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


.PHONY: app exe installer clean

app:
	python setup.py build_ext
	find build -type f -and \( -name native.so -or -name native.pyd \) -exec cp \{\} npamp/model/ \;
	pyside-rcc -o npamp/main_rc.py res/main.qrc
	pyside-uic -o npamp/mainwin.py npamp/mainwin.ui
	pyside-uic -o npamp/outwin.py npamp/outwin.ui

exe: app
	python setup.py py2exe

installer: exe
	ISCC installer/npamp.iss

clean:
	rm -rf output
	rm -f installer/*.exe
	rm -rf dist
	rm -rf build
	find . -type f -name \*.pyc -exec rm -f \{\} \;
	rm -f npamp/outwin.py npamp/mainwin.py npamp/main_rc.py
	rm -f npamp/model/native.so npamp/model/native.pyd npamp/model/native.c
