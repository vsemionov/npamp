
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

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy
import matplotlib

sys.path.append(os.path.join(os.getcwd(), "npamp"))
import meta


ext_modules = [Extension("pamp.native", ["npamp/pamp/native.pyx", "npamp/pamp/randomkit.c"], include_dirs = ["npamp/pamp", numpy.get_include()], extra_compile_args=["-DRK_NO_WINCRYPT"])]

setup_args = dict(
    name = meta.app_name,
    version = meta.app_version,
    author = meta.app_author_name,
    author_email = meta.app_author_email,
    url = meta.app_url,
    description = meta.app_description,
    cmdclass = {"build_ext": build_ext},
    ext_modules = ext_modules,
)

if os.name == "nt":
    import py2exe
    py2exe.__name__ #suppress warning for unused import
    
    data_files_mpl = matplotlib.get_py2exe_datafiles()
    data_files_npamp = [("examples", [os.path.join("examples", name) for name in os.listdir("examples")]), ]
    data_files = data_files_mpl + data_files_npamp
    
    options = {
        "py2exe": {
            "includes": [
                "pamp.native", # not required, but used because it generates an error if the module is not found
            ],
            "excludes": ["Tkinter", ],
            "bundle_files" : 2, # unreliable
        },
    }
    
    setup_args.update(dict(
        console = [{"script": "npamp/%s.py" % meta.app_name, "copyright": meta.app_copyright, }, ],
        windows = [{"script": "npamp/%s.py" % meta.gui_app_name, "icon_resources": [(1, "res/app.ico"), (2, "res/doc.ico"), ], "copyright": meta.app_copyright, }, ],
        data_files = data_files,
        options = options,
    ))


setup(**setup_args)
