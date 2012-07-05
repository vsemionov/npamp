
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


import math
import warnings


def pulse_scale(pulse, rel_fluence_loss):
    total_fluence = pulse.density_integral(float("inf"))
    scale = 1.0
    offset = pulse.offset
    half_duration = pulse.duration/2.0
    while True:
        integral_t0 = pulse.density_integral(offset - half_duration * scale)
        integral_t1 = pulse.density_integral(offset + half_duration * scale)
        fluence = integral_t1 - integral_t0
        trunc_fluence = total_fluence - fluence
        if trunc_fluence <= total_fluence * rel_fluence_loss:
            break
        scale *= 2.0
    return scale

def steps(divs):
    steps = 2**divs + 1
    return steps

def divs(min_steps):
    divs = int(math.ceil(math.log(min_steps - 1, 2.0)))
    return divs

def min_steps(min_counts, varspace, rtol, compute_result, compute_rdiff, retextra=False):
    min_count_x, min_count_y = min_counts
    var_x, var_y = varspace
    
    nvars = sum(varspace)
    max_divs_sums = [0, 15, 24 - 1]
    max_divs_sum = max_divs_sums[nvars]
    
    if not var_x and min_count_x == 1:
        divs_x = 0
        count_x = 1
    else:
        min_count_x = max(min_count_x, 3)
        divs_x = divs(min_count_x)
        count_x = steps(divs_x)
    if not var_y and min_count_y == 1:
        divs_y = 0
        count_y = 1
    else:
        min_count_y = max(min_count_y, 3)
        divs_y = divs(min_count_y)
        count_y = steps(divs_y)
    
    counts = count_x, count_y
    
    if divs_x + divs_y > max_divs_sum:
        warnings.warn("min. step counts (%d, %d) and corresponding min. divs (%d, %d) too large; max. divs sum: %d" % (count_x, count_y, divs_x, divs_y, max_divs_sum), stacklevel=2)
        return None
    
    result, rel_error = compute_result(*counts)
    
    rdiff = float("inf")
    
    while True:
        last_counts, last_result, last_rel_error = counts, result, rel_error
        last_divs_x, last_divs_y = divs_x, divs_x
        count_x, count_y = counts
        
        if var_x:
            counts_x = steps(divs_x+1), count_y
            result_x, rel_error_x = compute_result(*counts_x)
            rdiff_x = compute_rdiff(result, result_x)
        if var_y:
            counts_y = count_x, steps(divs_y+1)
            result_y, rel_error_y = compute_result(*counts_y)
            rdiff_y = compute_rdiff(result, result_y)
        
        if var_x and (not var_y or rdiff_x >= rdiff_y):
            divs_x += 1
            counts = counts_x
            result, rel_error = result_x, rel_error_x
            rdiff = rdiff_x
        elif var_y:
            divs_y += 1
            counts = counts_y
            result, rel_error = result_y, rel_error_y
            rdiff = rdiff_y
        
        rdiff *= 2.0
        
        if max(rdiff, last_rel_error) < rtol:
            break
        elif divs_x + divs_y > max_divs_sum or math.isinf(rdiff):
            warnings.warn("unable to reach rtol (%f); latest counts: (%d, %d); latest divs: (%d, %d); max. divs sum: %d; latest rel. error: %f, latest rel. difference: %f" % (rtol, count_x, count_y, last_divs_x, last_divs_y, max_divs_sum, rel_error, rdiff), stacklevel=2)
            break
    
    if retextra:
        return last_counts, rdiff, last_result
    return last_counts
