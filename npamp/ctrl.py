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


import params
import core
import ext


def execute(task_pool):
    dirname = "."
    
    if not params.initial_inversion:
        ref_inversion, inversion_rel_error = core.compute_inversion(dirname)
    else:
        ref_inversion, inversion_rel_error = params.initial_inversion, 0.0
    
    numerics = None
    if params.amplification:
        int_types, amp_types = params.integrator_classes, params.amplifier_classes
        numerics, rel_errors = core.select_methods((int_types, amp_types), ref_inversion, ret_rel_errors=True)
        num_types, counts = numerics
        core.amplify_ref_pulse(dirname, num_types, counts, ref_inversion)
        max_output_fluence, output_photon_counts, output_energy, rel_gain_decrease = core.amplify_train(dirname, num_types, counts, ref_inversion)
        core.report_results(ref_inversion, max_output_fluence, output_photon_counts, output_energy, rel_gain_decrease, inversion_rel_error, rel_errors)
    
    if params.extended_mode:
        ext.extended_mode(task_pool, dirname, ref_inversion, (int_types, amp_types), numerics)
