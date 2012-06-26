
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


def pulse_scale(pulse, rel_energy_loss):
    half_energy = pulse.density_integral(pulse.offset)
    scale = 1.0
    while pulse.density_integral(pulse.offset - pulse.duration/2.0 * scale) > half_energy * rel_energy_loss:
        scale *= 2.0
    return scale

def steps(divs):
    steps = 2**divs + 1
    return steps

def divs(min_steps):
    divs = int(math.ceil(math.log(min_steps - 1, 2.0)))
    return divs

def min_steps(min_count_x, min_count_y, var_x, var_y, rtol, compute_result, compute_rdiff, retextra=False):
    max_divs = 12
    
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
    
    if divs_x > max_divs or divs_y > max_divs:
        warnings.warn("minimum step count(s) (%d, %d) too large" % counts, stacklevel=2)
        return None
    
    result, rel_error = compute_result(*counts)
    
    rdiff = float("inf")
    
    while True:
        last_counts, last_result, last_rel_error = counts, result, rel_error
        
        count_x, count_y = counts
        if var_x and divs_x < max_divs:
            counts_x = steps(divs_x+1), count_y
            result_x, rel_error_x = compute_result(*counts_x)
        else:
            counts_x = counts
            result_x, rel_error_x = result, rel_error
        if var_y and divs_y < max_divs:
            counts_y = count_x, steps(divs_y+1)
            result_y, rel_error_y = compute_result(*counts_y)
        else:
            counts_y = counts
            result_y, rel_error_y = result, rel_error
        
        rdiff_x = compute_rdiff(result, result_x)
        rdiff_y = compute_rdiff(result, result_y)
        
        if result_x is not result and rdiff_x >= rdiff_y:
            divs_x += 1
            counts = counts_x
            result, rel_error = result_x, rel_error_x
            rdiff = rdiff_x
        elif result_y is not result:
            divs_y += 1
            counts = counts_y
            result, rel_error = result_y, rel_error_y
            rdiff = rdiff_y
        else:
            warnings.warn("unable to reach rtol (%f); latest counts: (%d, %d); max count: %d; latest rel. error: %f, latest rdiff: %f" % (rtol, count_x, count_y, steps(max_divs), rel_error, rdiff), stacklevel=2)
            break
        
        if max(rdiff, last_rel_error) < rtol:
            break
    
    if retextra:
        return last_counts, rdiff, last_result
    return last_counts
