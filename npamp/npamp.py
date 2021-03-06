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
import multiprocessing

import time
import glob
import getopt

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
"""Usage: {app_name} {{-h | -v | -e | [-g] [-D definition [...]] [-o output_dir] conf_file}}

Options:
 -h             print this help message and exit
 -v             print version information and exit
 -e             list installed extensions and exit
 -g             debug mode
 -D             define or override parameters
 -o output_dir  set output directory path (required for graphs)

Arguments:
 conf_file      path to configuration file ("-" for standard input)""".format(app_name=meta.app_name)

help_hint = "Try \"{app_name} -h\" for more information.".format(app_name=meta.app_name)


debug_mode = False


app_dir = os.path.join(os.path.expanduser("~"), ".%s" % meta.app_name.lower())


class InvocationError(Exception):
    pass


def print_help():
    print usage_help

def print_info():
    print meta.app_name, meta.app_version
    print
    print meta.app_copyright
    print meta.app_rights
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
    extension_path = os.path.join(app_dir, "extensions")
    sys.path.append(extension_path)
    extension_pathnames = glob.glob(os.path.join(extension_path, "*.py"))
    extensions = []
    for pathname in extension_pathnames:
        _, name = os.path.split(pathname)
        name, _ = os.path.splitext(name)
        if debug_mode:
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
        numerics, rel_errors = core.select_methods((int_types, amp_types), ref_inversion)
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
            cfg.update_conf(params.__dict__, params.__dict__, definitions)
        
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
    except (cfg.ConfigurationError, core.ComputationError, model.exc.ModelError) as exc:
        output.print_error(str(exc))
        sys.exit(1)
    except MemoryError:
        output.print_error("out of memory")
        sys.exit(1)
    except:
        if not debug_mode:
            output.print_exception()
            sys.exit(1)
        else:
            raise

def process():
    try:
        conf_path = None
        output_path = None
        definitions = []
        
        help_flag = False
        version_flag = False
        extensions_flag = False
        
        opts, args = getopt.getopt(sys.argv[1:], "hvegD:o:")
        for opt, arg in opts:
            if opt == "-h":
                help_flag = True
            elif opt == "-v":
                version_flag = True
            elif opt == "-e":
                extensions_flag = True
            elif opt == "-g":
                global debug_mode
                debug_mode = True
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
        
        extensions = load_extensions()
        
        if extensions_flag:
            print_extensions(extensions)
            sys.exit()
        
        if conf_path is None:
            raise InvocationError("no input file")
        if len(args) > 1:
            raise InvocationError("too many arguments")
        
        run(conf_path, output_path, definitions)
    except (getopt.GetoptError, InvocationError) as ie:
        output.print_error(str(ie), help_hint)
        sys.exit(2)

def main():
    multiprocessing.freeze_support()
    
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    
    process()


if __name__ == "__main__":
    main()
