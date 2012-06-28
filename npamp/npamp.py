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
import getopt
import time

import multiprocessing

import meta
import params
import core
import mp
import output
import ext
import svc



usage_help = \
"""Usage: {app_name} {{-h | -v | [-o output_dir] conf_file}}

Options:
 -h             print this help message and exit
 -v             print version information and exit
 -o output_dir  output directory path (required when graphs are enabled)""".format(app_name=meta.app_name)

help_hint = "Try \"{app_name} -h\" for more information.".format(app_name=meta.app_name)



class InvocationError(Exception):
    pass


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

def run(conf_path, output_path):
    start_time = time.clock()
    
    print "configuring"
    if conf_path is not None:
        if params.verbose:
            print "reading configuration from:", conf_path
        execfile(conf_path, params.__dict__)
    output.output_dir = output_path
    if params.graphs:
        if not output.output_dir:
            raise InvocationError("no output directory")
        if params.verbose:
            print "result graphs will be saved in:", output_path
    else:
        if params.verbose:
            print "result graphs will not be saved"
    
    dirname = "."
    
    int_types = params.integrator_classes
    amp_types = params.amplifier_classes
    
    if params.lower_process_priority:
        svc.lower_process_priority()
    
    with mp.TaskPool(params.num_tasks, params.extended_mode, params.__dict__) as task_pool:
        if params.verbose:
            print "number of parallel tasks:", task_pool.num_tasks
        
        print "validating"
        core.validate()
        
        print "running"
        
        if params.extended_mode:
            ext.extended_mode(task_pool, dirname, (int_types, amp_types))
        
        if not params.initial_inversion:
            ref_inversion, ref_inversion_rel_error = core.compute_inversion(dirname)
        else:
            ref_inversion, ref_inversion_rel_error = params.initial_inversion, 0.0
        
        if params.amplification:
            num_types, counts = core.setup_methods(dirname, (int_types, amp_types), ref_inversion)
            energy_rel_error = core.compute_energy_rel_error(ref_inversion, ref_inversion_rel_error)
            max_output_fluence, output_photon_counts, output_energy, rel_gain_reduction = core.amplify_train(dirname, num_types, counts, ref_inversion)
            core.report_output_characteristics(ref_inversion, max_output_fluence, output_photon_counts, output_energy, rel_gain_reduction, ref_inversion_rel_error, energy_rel_error)
    
    print output.div_line
    print "done"
    
    end_time = time.clock()
    elapsed_time = end_time - start_time
    if params.verbose:
        print "finished in %.3f seconds" % elapsed_time

def print_help():
    print usage_help

def main():
    multiprocessing.freeze_support()
    
    try:
        conf_path = None
        output_path = None
        help_flag = False
        version_flag = False
        opts, args = getopt.getopt(sys.argv[1:], "hvo:")
        for opt, arg in opts:
            if opt == "-o":
                output_path = arg
            elif opt == "-h":
                help_flag = True
            elif opt == "-v":
                version_flag = True
            else:
                assert False, "unhandled option"
        if help_flag:
            print_help()
            sys.exit()
        if version_flag:
            print_info()
            sys.exit()
        if len(args) < 1:
            raise InvocationError("no input file")
        elif len(args) > 1:
            raise InvocationError("too many arguments")
        if args:
            conf_path = args[0]
        run(conf_path, output_path)
    except InvocationError as ie:
        print >>sys.stderr, "%s: %s" % (meta.app_name, ie.message)
        print >>sys.stderr, help_hint


if __name__ == "__main__":
    main()
