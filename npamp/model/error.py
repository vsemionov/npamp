
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
import warnings

import numpy as np

import amplifier
import util


def pulse_scale(pulse, trunc_rtol):
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

def min_steps(min_counts, varspace, rtol, compute_result, compute_rdiff, ret_extra=False):
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
        divs_x = util.divs(min_count_x)
        count_x = util.steps(divs_x)
    if not var_y and min_count_y == 1:
        divs_y = 0
        count_y = 1
    else:
        min_count_y = max(min_count_y, 3)
        divs_y = util.divs(min_count_y)
        count_y = util.steps(divs_y)
    
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
            counts_x = util.steps(divs_x+1), count_y
            result_x, rel_error_x = compute_result(*counts_x)
            rdiff_x = compute_rdiff(result, result_x)
        if var_y:
            counts_y = count_x, util.steps(divs_y+1)
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
    
    if ret_extra:
        return last_counts, rdiff, last_result
    return last_counts

def min_integration_steps(integrator, beam_profile, single_pulses, energy_rtol, fluence_rtol):
    def converge_steps(func, a, b):
        max_divs = 12
        divs = 1
        last_res = None
        while True:
            steps = util.steps(divs)
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
                    warnings.warn("max divs (%d) exceeded; latest difference: %f (current: %f; last: %f)" % (max_divs, diff, res, last_res), stacklevel=3)
                    break
                else:
                    raise ValueError("max divs (%d) exceeded" % max_divs)
            last_res = res
        return steps
    def fluence_integrals(steps_rho, steps_phi):
        Rho = np.linspace(0.0, active_medium.radius, steps_rho)
        Phi = np.linspace(0.0, 2.0*math.pi, steps_phi)
        fluence_min = np.empty((steps_rho, steps_phi))
        fluence_max = np.empty((steps_rho, steps_phi))
        fluence_max_3 = np.empty((steps_rho, steps_phi))
        fluence_max_4 = np.empty((steps_rho, steps_phi))
        inversion_fluence = active_medium.initial_inversion.inversion_integral(0.0, 0.0, L)
        for m, rho in enumerate(Rho):
            for n, phi in enumerate(Phi):
                beam_fluence = beam_profile.fluence(rho, phi)
                fluence_min[m, n] = beam_fluence
                fluence_max[m, n] = beam_fluence * gain
                fluence_max_3[m, n] = beam_fluence + inversion_fluence / 2.0
                fluence_max_4[m, n] = beam_fluence + inversion_fluence
        num0 = integrator.integrate_base(Rho, Phi, fluence_min)
        exact0 = beam_profile.fluence_integral(active_medium.radius)
        rel_error_0 = abs((num0 - exact0) / exact0)
        num1 = integrator.integrate_base(Rho, Phi, fluence_max)
        exact1 = beam_profile.fluence_integral(active_medium.radius) * gain
        rel_error_1 = abs((num1 - exact1) / exact1)
        num3 = integrator.integrate_base(Rho, Phi, fluence_max_3)
        exact3 = exact0 + inversion_fluence / 2.0 * active_medium.aperture
        rel_error_3 = abs((num3 - exact3) / exact3)
        num4 = integrator.integrate_base(Rho, Phi, fluence_max_4)
        exact4 = exact0 + inversion_fluence * active_medium.aperture
        rel_error_4 = abs((num4 - exact4) / exact4)
        result = (num0, num1, num3, num4)
        rel_error = max(rel_error_0, rel_error_1, rel_error_3, rel_error_4)
        return result, rel_error
    
    active_medium = integrator.active_medium
    L = active_medium.length
    gain = math.exp(active_medium.doping_agent.xsection * active_medium.initial_inversion.ref_inversion * L)
    
    rtol = energy_rtol
    
    beam_rho, beam_phi = beam_profile.xcoords
    
    if not beam_rho and not beam_phi:
        steps_rho = 1
        steps_phi = 1
    else:
        compute_rdiff = lambda last_res, res: max([abs((res[i] - last_res[i]) / res[i]) for i in range(len(res))])
        steps_rho, steps_phi = min_steps((1, 1), (beam_rho, beam_phi), rtol, fluence_integrals, compute_rdiff)
    
    active_medium_3 = copy.deepcopy(active_medium)
    active_medium_3.doping_agent.lower_lifetime = float("inf")
    exact_amp_3 = amplifier.ExactAmplifier(active_medium_3, 2)
    exact_amp_3.rho = 0.0
    exact_amp_3.phi = 0.0
    
    active_medium_4 = copy.deepcopy(active_medium)
    active_medium_4.doping_agent.lower_lifetime = 0.0
    exact_amp_4 = amplifier.ExactAmplifier(active_medium_4, 2)
    exact_amp_4.rho = 0.0
    exact_amp_4.phi = 0.0
    
    rtol = fluence_rtol
    steps_t = 0
    for single_pulse in single_pulses:
        t0 = single_pulse.t0
        T = single_pulse.duration
        exact_amp_3.input_pulse = single_pulse
        exact_amp_4.input_pulse = single_pulse
        density_in = lambda t: single_pulse.density(t)
        density_out = lambda t: single_pulse.density(t) * gain
        density_out_3 = lambda t: exact_amp_3.exact_density(L, t)
        density_out_4 = lambda t: exact_amp_4.exact_density(L, t)
        steps_t_0 = converge_steps(density_in, t0, t0 + T)
        steps_t_1 = converge_steps(density_out, t0, t0 + T)
        steps_t_3 = converge_steps(density_out_3, t0, t0 + T)
        steps_t_4 = converge_steps(density_out_4, t0, t0 + T)
        steps_t = max(steps_t, steps_t_0, steps_t_1, steps_t_3, steps_t_4)
    
    rtol = fluence_rtol
    population_inversion = lambda population: population[0] - population[1]
    steps_z = 0
    for single_pulse in single_pulses:
        T = single_pulse.duration
        exact_amp_3.input_pulse = single_pulse
        exact_amp_4.input_pulse = single_pulse
        inversion_out_3 = lambda z: population_inversion(exact_amp_3.exact_population(z, T))
        inversion_out_4 = lambda z: population_inversion(exact_amp_4.exact_population(z, T))
        steps_z_3 = converge_steps(inversion_out_3, 0.0, L)
        steps_z_4 = converge_steps(inversion_out_4, 0.0, L)
        steps_z = max(steps_z, steps_z_3, steps_z_4)
    
    return steps_rho, steps_phi, steps_z, steps_t

def min_amplification_steps(amp_type, active_medium, pulse, decay, (min_count_z, min_count_t), integrator, fluence_rtol, amp_rtol, ret_extra=False):
    def get_rel_error((num_density_out, num_population_out, num_T, num_Z), (exact_density_out, exact_population_out, exact_T, exact_Z), active_medium):
        num_inversion = num_population_out[0] - num_population_out[1] * decay
        exact_inversion = exact_population_out[0] - exact_population_out[1] * decay
        num_density_integral = integrator.integrate(num_T, num_density_out)
        exact_density_integral = integrator.integrate(exact_T, exact_density_out)
        rel_error_density = abs((num_density_integral - exact_density_integral) / exact_density_integral)
        rel_error_density += fluence_rtol * (1.0 + num_density_integral / exact_density_integral)
        num_inversion_integral = integrator.integrate(num_Z, num_inversion)
        exact_inversion_integral = integrator.integrate(exact_Z, exact_inversion)
        inversion_abs_error = abs(num_inversion_integral - exact_inversion_integral) + fluence_rtol * (num_inversion_integral + exact_inversion_integral)
        rel_error_inversion = math.exp(active_medium.doping_agent.xsection * inversion_abs_error) - 1.0
        rel_error = rel_error_density + rel_error_inversion
        return rel_error
    def get_rel_error_amp(num_amp, exact_amp):
        assert num_amp.active_medium is exact_amp.active_medium
        num_density_out = num_amp.density[-1]
        num_population_out = (num_amp.population[0].T[-1], num_amp.population[1].T[-1])
        num_T = num_amp.T
        num_Z = num_amp.Z
        exact_density_out = exact_amp.density[-1]
        exact_population_out = (exact_amp.population[0].T[-1], exact_amp.population[1].T[-1])
        exact_T = exact_amp.T
        exact_Z = exact_amp.Z
        active_medium = num_amp.active_medium
        rel_error = get_rel_error((num_density_out, num_population_out, num_T, num_Z), (exact_density_out, exact_population_out, exact_T, exact_Z), active_medium)
        return rel_error
    def amplify_pulse(count_z, count_t):
        rel_error = 0.0
        analytical_lower_lifetimes = [float("inf"), 0.0]
        for lower_lifetime in analytical_lower_lifetimes:
            test_active_medium = copy.deepcopy(active_medium)
            test_active_medium.doping_agent.lower_lifetime = lower_lifetime
            amp = amp_type(test_active_medium, count_z)
            num_density_out, num_population_out = amp.amplify(0.0, 0.0, pulse, count_t)
            exact = amplifier.ExactOutputAmplifier(test_active_medium, count_z)
            exact_density_out, exact_population_out = exact.amplify(0.0, 0.0, pulse, count_t)
            test_rel_error = get_rel_error((num_density_out, num_population_out, amp.T, amp.Z), (exact_density_out, exact_population_out, exact.T, exact.Z), test_active_medium)
            rel_error = max(test_rel_error, rel_error)
        amp = amp_type(active_medium, count_z)
        amp.amplify(0.0, 0.0, pulse, count_t)
        if active_medium.doping_agent.lower_lifetime in analytical_lower_lifetimes:
            exact = amplifier.ExactOutputAmplifier(active_medium, count_z)
            exact_density_out, exact_population_out = exact.amplify(0.0, 0.0, pulse, count_t)
        else:
            exact_density_out, exact_population_out = None, None
        results = (amp, exact_density_out, exact_population_out)
        return results, rel_error
    
    min_count_z = max(min_count_z, 3)
    min_count_t = max(min_count_t, 3)
    compute_rdiff = lambda last_res, res: get_rel_error_amp(last_res[0], res[0])
    data = min_steps((min_count_z, min_count_t), (True, True), amp_rtol, amplify_pulse, compute_rdiff, ret_extra)
    
    return data