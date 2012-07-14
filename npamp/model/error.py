
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




def min_steps(min_counts, varspace, rtol, limit, compute_result, compute_rdiff, integrator, opname, varnames):
    min_count_x, min_count_y = min_counts
    var_x, var_y = varspace
    
    max_divs_sums = [0, 16 - 1, 24 - 1]
    nvars = 0
    
    min_count_int = integrator.num_integrator.min_count
    
    if not var_x and min_count_x <= 1:
        divs_x = 0
        count_x = 1
    else:
        nvars += 1
        min_count_x = max(min_count_x, min_count_int)
        divs_x = discrete.divs(min_count_x)
        count_x = discrete.steps(divs_x)
    
    if not var_y and min_count_y <= 1:
        divs_y = 0
        count_y = 1
    else:
        nvars += 1
        min_count_y = max(min_count_y, min_count_int)
        divs_y = discrete.divs(min_count_y)
        count_y = discrete.steps(divs_y)
    
    max_divs_sum = max_divs_sums[nvars]
    
    if limit and (count_x * count_y) > limit:
        raise exc.SoftLimitError()
    
    if (divs_x + divs_y) > max_divs_sum:
        raise exc.NumericalError("min. %s %s step counts (%s, %s) and corresponding min. divs (%s, %s) too large; max. divs sum: %s" % (opname, varnames, count_x, count_y, divs_x, divs_y, max_divs_sum))
    
    if not var_x and not var_y:
        return (count_x, count_y), 0.0
    
    result, rel_error = compute_result(count_x, count_y)
    rdiff = float("inf")
    
    while True:
        last_divs_x, last_divs_y = divs_x, divs_y
        last_count_x, last_count_y = count_x, count_y
        last_rel_error = rel_error
        
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
            count_x, count_y = counts_x
            result, rel_error = result_x, rel_error_x
            rdiff = rdiff_x
        elif var_y:
            divs_y += 1
            count_x, count_y = counts_y
            result, rel_error = result_y, rel_error_y
            rdiff = rdiff_y
        else:
            assert False, "unhandled case"
        
        rdiff *= 2.0
        max_rel_error = max(last_rel_error, rdiff)
        
        if max_rel_error <= rtol:
            break
        
        if limit and (count_x * count_y) > limit:
            raise exc.SoftLimitError()
        
        if (divs_x + divs_y) > max_divs_sum:
            util.warn("max. %s %s divs sum (%s) exceeded; rtol: %s; latest counts: (%s, %s); latest divs: (%s, %s); latest rel. error: %s, latest rel. difference: %s" % (opname, varnames, max_divs_sum, rtol, last_count_x, last_count_y, last_divs_x, last_divs_y, last_rel_error, rdiff), stacklevel=2)
            break
    
    return (last_count_x, last_count_y), max_rel_error


_compute_seq_rdiff = lambda approx_res, exact_res: max([abs((exact - approx) / exact) for (approx, exact) in zip(approx_res, exact_res)])


def min_integration_steps(integrator, (min_count_rho, min_count_phi), int_rtol, limit, active_medium, input_beam):
    def fluence_integrals(count_rho, count_phi):
        photon_count_in = input_beam.fluence_integral(active_medium.radius)
        inversion_fluence_integral = active_medium.initial_inversion.inversion_fluence_integral(active_medium.radius, active_medium.length)
        exact_results = (
            photon_count_in * gain,
            photon_count_in + inversion_fluence_integral / 2.0,
            photon_count_in + inversion_fluence_integral,
        )
        
        beam_fluence = lambda rho, phi: input_beam.fluence(rho, phi)
        inversion_fluence = lambda rho, phi: active_medium.initial_inversion.inversion_integral(rho, phi, active_medium.length)
        output_fluences = (
            lambda rho, phi: beam_fluence(rho, phi) * gain,
            lambda rho, phi: beam_fluence(rho, phi) + inversion_fluence(rho, phi) / 2.0,
            lambda rho, phi: beam_fluence(rho, phi) + inversion_fluence(rho, phi),
        )
        
        Rho = np.linspace(0.0, active_medium.radius, count_rho)
        Phi = np.linspace(0.0, 2.0*math.pi, count_phi)
        discrete_fluence = np.empty((count_rho, count_phi))
        
        num_results = []
        for output_fluence in output_fluences:
            for m, rho in enumerate(Rho):
                for n, phi in enumerate(Phi):
                    discrete_fluence[m, n] = output_fluence(rho, phi)
            num_result = integrator.integrate_base(Rho, Phi, discrete_fluence)
            num_results.append(num_result)
        
        rel_error = _compute_seq_rdiff(num_results, exact_results)
        return num_results, rel_error
    
    gain = math.exp(active_medium.doping_agent.xsection * active_medium.initial_inversion.ref_inversion * active_medium.length)
    
    (count_rho, count_phi), rel_error = min_steps((min_count_rho, min_count_phi), input_beam.xcoords, int_rtol, limit, fluence_integrals, _compute_seq_rdiff, integrator, "integration", "(rho, phi)")
    
    return (count_rho, count_phi), rel_error


def min_amplification_steps(amp_type, (min_count_z, min_count_t), amp_rtol, limit, active_medium, pulse_train, xverse, integrator):
    def amplify_signal(count_z, count_t):
        pulse_count = pulse_train.count
        lower_decay = amplifier.lower_state_decay(active_medium, pulse_train)
        
        num_fluences = []
        exact_fluences = []
        for lower_lifetime in amplifier.ExactAmplifier.analytical_lower_lifetimes:
            test_active_medium = copy.deepcopy(active_medium)
            test_active_medium.doping_agent.lower_lifetime = lower_lifetime
            
            amp = amp_type(test_active_medium, count_z)
            num_density_out, num_population_final = amp.amplify(rho, phi, ref_pulse, count_t)
            num_fluence = integrator.integrate(amp.T, num_density_out)
            num_fluences.append(num_fluence)
            del amp, num_density_out, num_population_final
            
            exact_amp = amplifier.ExactOutputAmplifier(test_active_medium, count_z)
            exact_density_out, exact_population_final = exact_amp.amplify(rho, phi, ref_pulse, count_t)
            exact_fluence = integrator.integrate(exact_amp.T, exact_density_out)
            exact_fluences.append(exact_fluence)
            del exact_amp, exact_density_out, exact_population_final
        
        rel_error = _compute_seq_rdiff(num_fluences, exact_fluences)
        
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
        
        fluence_out_first = pulse_fluences[0]
        fluence_out_total = pulse_fluences[::-1].sum()
        
        return (fluence_out_first, fluence_out_total), rel_error
    
    initial_inversion = active_medium.initial_inversion
    ref_pulse = pulse_train.pulse
    
    if xverse:
        rho, phi = xverse
    else:
        rho, phi = initial_inversion.rho_ref, initial_inversion.phi_ref
    
    min_count_z = max(min_count_z, amp_type.min_steps_z(active_medium))
    min_count_t = max(min_count_t, amp_type(active_medium, min_count_z).min_steps_t(ref_pulse))
    
    (count_z, count_t), rel_error = min_steps((min_count_z, min_count_t), (True, True), amp_rtol, limit, amplify_signal, _compute_seq_rdiff, integrator, "amplification", "(z, t)")
    
    return (count_z, count_t), rel_error




def perturbed_inversion_rel_error(ref_inversion, perturb_ref_inversion, inversion_rtol):
    abs_error = abs(perturb_ref_inversion - ref_inversion) + (ref_inversion + perturb_ref_inversion) * inversion_rtol
    rel_error = abs_error / ref_inversion
    return rel_error

def energy_rel_error(active_medium, ref_inversion_rel_error, (time_trunc_rel_error, amp_rtol, int_rtol)):
    rel_error_inversion = math.exp(active_medium.doping_agent.xsection * ref_inversion_rel_error * active_medium.initial_inversion.ref_inversion * active_medium.length) - 1.0
    rel_error_energy = time_trunc_rel_error + amp_rtol + int_rtol
    energy_rel_error = rel_error_inversion + rel_error_energy
    return energy_rel_error
