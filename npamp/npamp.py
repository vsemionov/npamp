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

import time
import glob
import getopt

import signal
import multiprocessing

import model

import meta
import params
import output
import core
import ext
import cfg
import mp
import svc


usage_help = \
"""Usage: {app_name} {{-h | -v | -e | [-D definition [...]] [-o output_dir] conf_file}}

Options:
 -h             print this help message and exit
 -v             print version information and exit
 -e             list installed extensions
 -D             define or override parameters
 -o output_dir  set output directory path (required for graphs)""".format(app_name=meta.app_name)

help_hint = "Try \"{app_name} -h\" for more information.".format(app_name=meta.app_name)


class InvocationError(Exception):
    pass

def print_help():
    print usage_help

def print_info():
    print meta.app_name, meta.app_version
    print meta.app_copyright
    print
    print meta.app_description
    print
    print meta.app_author_msg
    print meta.app_coauthors_msg
    print
    print meta.app_website_msg

def print_extensions(extensions):
    if not extensions:
        print "no extensions installed"
    else:
        print "installed extensions (name: description):"
        for extension in extensions:
            name, doc = extension.__name__, extension.__doc__
            print "%s: %s" % (name, doc)

def load_extensions():
    extension_path = os.path.normpath(os.path.expanduser("~/.%s/extensions" % meta.app_name.lower()))
    sys.path.append(extension_path)
    extension_pathnames = glob.glob(os.path.join(extension_path, "*.py"))
    extensions = []
    for pathname in extension_pathnames:
        _, name = os.path.split(pathname)
        name, _ = os.path.splitext(name)
        print "loading extension \"%s\"" % name
        extension = __import__(name)
        extensions.append(extension)
    return extensions

def execute(task_pool):
    dirname = "."
    
    if not params.initial_inversion:
        ref_inversion, inversion_rel_error = core.compute_inversion(dirname)
    else:
        ref_inversion, inversion_rel_error = params.initial_inversion, 0.0
    
    int_types, amp_types = params.integrator_classes, params.amplifier_classes
    
    numerics = None
    if params.amplification:
        numerics, rel_errors = core.select_methods((int_types, amp_types), ref_inversion, ret_rel_errors=True)
        num_types, counts = numerics
        core.amplify_ref_pulse(dirname, num_types, counts, ref_inversion)
        max_output_fluence, output_photon_counts, output_energy, rel_gain_decrease = core.amplify_train(dirname, num_types, counts, ref_inversion)
        core.report_results(ref_inversion, max_output_fluence, output_photon_counts, output_energy, rel_gain_decrease, inversion_rel_error, rel_errors)
    
    if params.extended_mode:
        ext.extended_mode(task_pool, dirname, ref_inversion, (int_types, amp_types), numerics)

def run(conf_path, output_path, definitions):
    try:
        print "configuring"
        if conf_path is not None:
            conf = cfg.load_conf(params.__dict__, conf_path)
            params.__dict__.update(conf)
        if definitions is not None:
            diff = dict()
            for definition in definitions:
                exec definition in diff
            defaults = cfg.copy_conf(params.__dict__)
            diff = cfg.copy_conf(diff)
            cfg.check_conf(defaults, diff)
            params.__dict__.update(diff)
        
        if params.verbose:
            print "verbose output enabled"
        
        if params.lower_process_priority:
            if params.verbose:
                print "reducing process priority"
            warning = svc.reduce_process_priority()
            if warning:
                output.warn(warning)
        
        if params.graphs:
            if not output_path:
                raise InvocationError("no output directory")
            output_path = os.path.normpath(output_path)
            if params.verbose:
                print "graphs will be written to:", output_path
        output.output_dir = output_path
        
        print "validating"
        cfg.validate()
        
        with mp.TaskPool(params.num_tasks, params.extended_mode, params.__dict__) as task_pool:
            print "executing"
            
            if params.verbose:
                print "number of parallel tasks:", task_pool.num_tasks
            
            start_time = time.time()
            execute(task_pool)
            end_time = time.time()
            elapsed_time = end_time - start_time
            
            print output.div_line
            print "done"
            
            if params.verbose:
                print "finished in %.2f s" % elapsed_time
    except (cfg.ConfigurationError, core.ComputationError, model.exc.ModelError):
        output.print_error()

def process(extensions):
    try:
        conf_path = None
        output_path = None
        definitions = []
        
        help_flag = False
        version_flag = False
        extensions_flag = False
        
        opts, args = getopt.getopt(sys.argv[1:], "hveD:o:")
        for opt, arg in opts:
            if opt == "-h":
                help_flag = True
            elif opt == "-v":
                version_flag = True
            elif opt == "-e":
                extensions_flag = True
            elif opt == "-D":
                definitions.append(arg)
            elif opt == "-o":
                output_path = arg
            else:
                assert False, "unhandled option"
        if args:
            conf_path = args[0]
        
        if help_flag:
            print_help()
            sys.exit()
        if version_flag:
            print_info()
            sys.exit()
        if extensions_flag:
            print_extensions(extensions)
            sys.exit()
        
        if len(args) < 1:
            raise InvocationError("no input file")
        elif len(args) > 1:
            raise InvocationError("too many arguments")
        
        run(conf_path, output_path, definitions)
    except InvocationError:
        output.print_error(help_hint)

def main():
    multiprocessing.freeze_support()
    
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    
    extensions = load_extensions()
    process(extensions)


if __name__ == "__main__":
    main()
