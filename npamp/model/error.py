
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

import copy

import numpy as np

import amplifier
import discrete
import util
import exc


def pulse_scale(pulse, trunc_rtol):
    if pulse.ref_density == 0.0:
        return 1.0, 0.0
    total_fluence = pulse.density_integral(float("inf"))
    scale = 1.0
    offset = pulse.offset
    half_duration = pulse.duration/2.0
    while True:
        integral_t0 = pulse.density_integral(offset - half_duration * scale)
        integral_t1 = pulse.density_integral(offset + half_duration * scale)
        fluence = integral_t1 - integral_t0
        trunc_rel_error = (total_fluence - fluence) / total_fluence
        if trunc_rel_error <= trunc_rtol:
            break
        scale *= 2.0
    return scale, trunc_rel_error

def min_steps(min_counts, varspace, rtol, compute_result, compute_rdiff, opname, varnames):
    min_count_x, min_count_y = min_counts
    var_x, var_y = varspace
    
    nvars = sum(varspace)
    max_divs_sums = [0, 15, 24 - 1]
    max_divs_sum = max_divs_sums[nvars]
    
    if not var_x and min_count_x <= 1:
        divs_x = 0
        count_x = 1
    else:
        min_count_x = max(min_count_x, 3)
        divs_x = discrete.divs(min_count_x)
        count_x = discrete.steps(divs_x)
    if not var_y and min_count_y <= 1:
        divs_y = 0
        count_y = 1
    else:
        min_count_y = max(min_count_y, 3)
        divs_y = discrete.divs(min_count_y)
        count_y = discrete.steps(divs_y)
    
    counts = count_x, count_y
    
    if divs_x + divs_y > max_divs_sum:
        raise exc.NumericalError("min. %s %s step counts (%s, %s) and corresponding min. divs (%s, %s) too large; max. divs sum: %s" % (opname, varnames, count_x, count_y, divs_x, divs_y, max_divs_sum))
    
    result, rel_error = compute_result(*counts)
    
    rdiff = float("inf")
    
    while True:
        last_counts, last_rel_error = counts, rel_error
        last_divs_x, last_divs_y = divs_x, divs_y
        count_x, count_y = counts
        
        if var_x:
            counts_x = discrete.steps(divs_x+1), count_y
            result_x, rel_error_x = compute_result(*counts_x)
            rdiff_x = compute_rdiff(result, result_x)
        if var_y:
            counts_y = count_x, discrete.steps(divs_y+1)
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
        max_rel_error = max(rdiff, last_rel_error)
        
        if max(rdiff, last_rel_error) < rtol:
            break
        elif divs_x + divs_y > max_divs_sum or math.isinf(rdiff):
            util.warn("max. %s %s divs sum (%s) exceeded; rtol: %s; latest counts: (%s, %s); latest divs: (%s, %s); latest rel. error: %s, latest rel. difference: %s" % (opname, varnames, max_divs_sum, rtol, count_x, count_y, last_divs_x, last_divs_y, rel_error, rdiff), stacklevel=2)
            break
    
    return last_counts, max_rel_error

def min_integration_steps(integrator, input_beam, int_rtol, (min_count_rho, min_count_phi)):
    def converge_steps(func, a, b, rtol, min_steps, varname):
        max_divs = 15
        
        min_steps = max(min_steps, 3)
        divs = discrete.divs(min_steps)
        if divs > max_divs:
            raise exc.NumericalError("min. integration %s step count (%s) and corresponding min. divs (%s) too large; max. divs: %s" % (varname, min_steps, divs, max_divs))
        
        divs = max(divs-1, 1)
        last_res = None
        while True:
            steps = discrete.steps(divs)
            X = np.linspace(a, b, steps)
            Y = np.vectorize(func)(X)
            res = integrator.integrate(X, Y)
            
            diff = None
            if last_res is not None:
                diff = abs(res - last_res)
                if diff <= rtol * abs(res):
                    break
            
            divs += 1
            if divs > max_divs:
                if diff is not None:
                    util.warn("max. integration %s divs (%s) exceeded; rtol: %s; latest step count: %s; latest difference: %s (current: %s; last: %s)" % (varname, max_divs, rtol, steps, diff, res, last_res), stacklevel=3)
                else:
                    util.warn("max. integration %s divs (%s) exceeded; rtol: %s; latest step count: %s" % (varname, max_divs, rtol, steps), stacklevel=3)
                break
            last_res = res
        
        return steps
    
    def fluence_integrals(steps_rho, steps_phi):
        Rho = np.linspace(0.0, active_medium.radius, steps_rho)
        Phi = np.linspace(0.0, 2.0*math.pi, steps_phi)
        
        fluence_min = np.empty((steps_rho, steps_phi))
        fluence_max = np.empty((steps_rho, steps_phi))
        fluence_max_3 = np.empty((steps_rho, steps_phi))
        fluence_max_4 = np.empty((steps_rho, steps_phi))
        
        for m, rho in enumerate(Rho):
            for n, phi in enumerate(Phi):
                beam_fluence = input_beam.fluence(rho, phi)
                inversion_fluence = active_medium.initial_inversion.inversion_integral(rho, phi, active_medium.length)
                fluence_min[m, n] = beam_fluence
                fluence_max[m, n] = beam_fluence * gain
                fluence_max_3[m, n] = beam_fluence + inversion_fluence / 2.0
                fluence_max_4[m, n] = beam_fluence + inversion_fluence
        
        photon_count_in = input_beam.fluence_integral(active_medium.radius)
        inversion_fluence_integral = active_medium.initial_inversion.inversion_fluence_integral(active_medium.radius, active_medium.length)
        
        num0 = integrator.integrate_base(Rho, Phi, fluence_min)
        exact0 = photon_count_in
        rel_error_0 = abs((num0 - exact0) / exact0)
        
        num1 = integrator.integrate_base(Rho, Phi, fluence_max)
        exact1 = photon_count_in * gain
        rel_error_1 = abs((num1 - exact1) / exact1)
        
        num3 = integrator.integrate_base(Rho, Phi, fluence_max_3)
        exact3 = photon_count_in + inversion_fluence_integral / 2.0
        rel_error_3 = abs((num3 - exact3) / exact3)
        
        num4 = integrator.integrate_base(Rho, Phi, fluence_max_4)
        exact4 = photon_count_in + inversion_fluence_integral
        rel_error_4 = abs((num4 - exact4) / exact4)
        
        result = (num0, num1, num3, num4)
        rel_error = max(rel_error_0, rel_error_1, rel_error_3, rel_error_4)
        return result, rel_error
    
    active_medium = integrator.active_medium
    gain = math.exp(active_medium.doping_agent.xsection * active_medium.initial_inversion.ref_inversion * active_medium.length)
    
    beam_rho, beam_phi = input_beam.xcoords
    
    if not beam_rho and not beam_phi:
        steps_rho, steps_phi = 1, 1
        rel_error = 0.0
    else:
        compute_rdiff = lambda last_res, res: max([abs((new - old) / new) for (old, new) in zip(last_res, res)])
        (steps_rho, steps_phi), rel_error = min_steps((min_count_rho, min_count_phi), (beam_rho, beam_phi), int_rtol, fluence_integrals, compute_rdiff, "integration", "(rho, phi)")
    
    return (steps_rho, steps_phi), rel_error

def min_amplification_steps(amp_type, active_medium, xverse, pulse_train, (min_count_z, min_count_t), integrator, amp_rtol):
    compute_rel_error = lambda num_fluence, exact_fluence: abs((exact_fluence - num_fluence) / exact_fluence)
    
    def amplify_signal(count_z, count_t):
        ref_pulse = pulse_train.pulse
        pulse_count = pulse_train.count
        lower_decay = amplifier.lower_state_decay(active_medium, pulse_train)
        
        rel_error = None
        if pulse_count == 1:
            for lower_lifetime in amplifier.ExactAmplifier.analytical_lower_lifetimes:
                test_active_medium = copy.deepcopy(active_medium)
                test_active_medium.doping_agent.lower_lifetime = lower_lifetime
                
                amp = amp_type(test_active_medium, count_z)
                num_density_out, num_population_final = amp.amplify(rho, phi, ref_pulse, count_t)
                num_fluence = integrator.integrate(amp.T, num_density_out)
                del amp, num_density_out, num_population_final
                
                exact_amp = amplifier.ExactOutputAmplifier(test_active_medium, count_z)
                exact_density_out, exact_population_final = exact_amp.amplify(rho, phi, ref_pulse, count_t)
                exact_fluence = integrator.integrate(exact_amp.T, exact_density_out)
                del exact_amp, exact_density_out, exact_population_final
                
                test_rel_error = compute_rel_error(num_fluence, exact_fluence)
                rel_error = max(test_rel_error, rel_error)
        
        amp = amp_type(active_medium, count_z)
        
        upper = np.vectorize(initial_inversion.inversion)(rho, phi, amp.Z)
        lower = np.zeros(count_z)
        population = (upper, lower)
        
        pulse_fluences = np.empty(pulse_count)
        
        for pnum in range(pulse_count):
            density_out, population_final = amp.amplify(rho, phi, ref_pulse, count_t, initial_population=population)
            
            upper = np.copy(population_final[0])
            lower = population_final[1] * lower_decay
            population = (upper, lower)
            
            fluence_out = integrator.integrate(amp.T, density_out)
            pulse_fluences[pnum] = fluence_out
        
        fluence_out = pulse_fluences[::-1].sum()
        
        return fluence_out, rel_error
    
    initial_inversion = active_medium.initial_inversion
    
    if xverse:
        rho, phi = xverse
    else:
        rho, phi = initial_inversion.rho_ref, initial_inversion.phi_ref
    
    data = min_steps((min_count_z, min_count_t), (True, True), amp_rtol, amplify_signal, compute_rel_error, "amplification", "(z, t)")
    
    return data

def perturbed_inversion_rel_error(ref_inversion, perturb_ref_inversion, inversion_rtol):
    abs_error = abs(perturb_ref_inversion - ref_inversion) + (ref_inversion + perturb_ref_inversion) * inversion_rtol
    rel_error = abs_error / ref_inversion
    return rel_error

def energy_rel_error(active_medium, ref_inversion_rel_error, (time_trunc_rel_error, amp_rtol, int_rtol)):
    rel_error_inversion = math.exp(active_medium.doping_agent.xsection * ref_inversion_rel_error * active_medium.initial_inversion.ref_inversion * active_medium.length) - 1.0
    rel_error_energy = time_trunc_rel_error + amp_rtol + int_rtol
    energy_rel_error = rel_error_inversion + rel_error_energy
    return energy_rel_error
