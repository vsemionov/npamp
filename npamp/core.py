
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


import os

import math
import copy
import warnings

import types

import numpy as np

import params
import pamp
import output
import unitconv


class ConfigurationError(Exception):
    pass

class ComputationError(Exception):
    pass


copy_conf = lambda conf: {k: v for k, v in conf.items() if not k.startswith('_') and type(v) is not types.ModuleType}

def create_pulse(active_medium, beam, rho, phi):
    fluence = beam.fluence(rho, phi)
    ref_density = params.pulse_class.ref_density(active_medium.light_speed, params.pulse_duration, fluence)
    pulse = params.pulse_class(-params.pulse_duration/2.0, params.pulse_duration, ref_density)
    scale = pamp.util.pulse_scale(pulse, params.pulse_rel_energy_trunc)
    pulse = pamp.pulse.ExtendedPulse(pulse, scale)
    pulse = pamp.pulse.TruncatedPulse(pulse)
    return pulse

def mangle_count_z(count_z):
    count_z = max(count_z, 3)
    count_z = pamp.util.steps(pamp.util.divs(count_z))
    return count_z

def mangle_count_t(count_t):
    count_t = max(count_t, 3)
    count_t = pamp.util.steps(pamp.util.divs(count_t))
    return count_t

def mangle_count_cyl(count_cyl):
    count_cyl = max(count_cyl, 1)
    if count_cyl > 1:
        count_cyl = max(count_cyl, 3)
        count_cyl = pamp.util.steps(pamp.util.divs(count_cyl))
    return count_cyl

def compute_inversion(dirname):
    print output.div_line
    print "computing initial inversion"
    
    is_numerical = issubclass(params.loss_model_class, pamp.loss.NumericalLossModel)
    doping_agent = pamp.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = pamp.medium.ActiveMedium(None, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    pump_system = pamp.pump.PumpSystem(params.pump_wavelen, params.pump_duration, params.pump_power, params.pump_efficiency)
    loss_model_kwargs = dict(rtol=params.loss_rtol) if is_numerical else {}
    if is_numerical:
        loss_model_kwargs = dict(rtol=params.loss_rtol, min_count=params.loss_min_count)
    else:
        params.loss_rtol = 0.0
        loss_model_kwargs = {}
    loss_model_kwargs.update(params.loss_model_extra_args)
    loss_model = params.loss_model_class(active_medium, params.lasing_wavelen, **loss_model_kwargs)
    inv = params.inverter_class(active_medium, pump_system, loss_model)
    ref_inversion = inv.invert(params.initial_inversion_rtol, params.initial_inversion_min_count_t)
    rate_evals = (len(inv.inversion) - 1) * inv.evals_per_step
    pump_energy = params.pump_duration * params.pump_power
    stored_energy = pamp.energy.energy(params.lasing_wavelen, ref_inversion * active_medium.volume)
    
    if params.initial_inversion_validate:
        print "validating uniform ASE loss approximation"
        ross_num_model = loss_model if isinstance(loss_model, pamp.loss.RossNumericalASEModel) else pamp.loss.RossNumericalASEModel(active_medium, params.lasing_wavelen, params.loss_rtol, params.loss_min_count)
        rate_rel_stddev = ross_num_model.rate_rel_stddev(ref_inversion)
        unitconv.print_result("depopulation rate rel. std. deviation [{}]: {}", ("%",), (rate_rel_stddev,))
        if rate_rel_stddev > 10.0e-2:
            warnings.warn("uniform ASE loss approximation is invalid", stacklevel=2)
    
    if is_numerical:
        print "perturbing initial inversion"
        perturb_loss_model = pamp.loss.PerturbedLossModel(loss_model)
        perturb_inv = params.inverter_class(active_medium, pump_system, perturb_loss_model)
        perturb_ref_inversion = perturb_inv.invert(params.initial_inversion_rtol, params.initial_inversion_min_count_t)
        abs_error = abs(perturb_ref_inversion - ref_inversion) + (ref_inversion + perturb_ref_inversion) * params.initial_inversion_rtol
        rel_error = abs_error / ref_inversion
    else:
        rel_error = params.initial_inversion_rtol
    
    gain_coef = ref_inversion * active_medium.doping_agent.xsection
    gain = math.exp(gain_coef * active_medium.length)
    
    ref_inversion_atol = ref_inversion * rel_error
    gain_atol = gain * (math.exp(ref_inversion_atol * active_medium.doping_agent.xsection * active_medium.length) - 1.0)
    stored_energy_atol = stored_energy * rel_error
    
    if params.verbose:
        print "last number of depopulation rate evaluations:", rate_evals
    unitconv.print_result("pump energy [{}]: {}", ("mJ",), (pump_energy,))
    unitconv.print_result("initial inversion [{}]: {} ~ {}", ("cm^-3",), (ref_inversion, ref_inversion_atol))
    unitconv.print_result("small signal gain: {} ~ {}", (), (gain, gain_atol))
    unitconv.print_result("stored energy [{}]: {} ~ {}", ("mJ",), (stored_energy, stored_energy_atol))
    if params.graphs:
        print "generating output"
        dirname = output.init_dir(dirname)
        output.plot_inversion(dirname, inv)
    
    return ref_inversion, rel_error

def get_decay(active_medium, pulse):
    separation = params.train_pulse_period - pulse.duration
    lower_lifetime = active_medium.doping_agent.lower_lifetime
    if lower_lifetime == 0.0:
        decay = 0.0
    elif math.isinf(lower_lifetime):
        decay = 1.0
    elif params.train_pulse_count == 1:
        decay = 0.0
    else:
        decay = math.exp(- separation / lower_lifetime)
    return decay

def test_pulse(amp_type, active_medium, pulse, (min_count_z, min_count_t), integrator):
    def get_rel_error((num_density_out, num_population_out, num_T, num_Z), (exact_density_out, exact_population_out, exact_T, exact_Z), active_medium):
        decay = get_decay(active_medium, pulse)
        num_inversion = num_population_out[0] - num_population_out[1] * decay
        exact_inversion = exact_population_out[0] - exact_population_out[1] * decay
        num_density_integral = integrator.integral(num_T, num_density_out)
        exact_density_integral = integrator.integral(exact_T, exact_density_out)
        rel_error_density = abs((num_density_integral - exact_density_integral) / exact_density_integral)
        rel_error_density += params.fluence_rtol * (1.0 + num_density_integral / exact_density_integral)
        num_inversion_integral = integrator.integral(num_Z, num_inversion)
        exact_inversion_integral = integrator.integral(exact_Z, exact_inversion)
        inversion_abs_error = abs(num_inversion_integral - exact_inversion_integral) + params.fluence_rtol * (num_inversion_integral + exact_inversion_integral)
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
            exact = pamp.amplifier.ExactOutputAmplifier(test_active_medium, count_z)
            exact_density_out, exact_population_out = exact.amplify(0.0, 0.0, pulse, count_t)
            test_rel_error = get_rel_error((num_density_out, num_population_out, amp.T, amp.Z), (exact_density_out, exact_population_out, exact.T, exact.Z), test_active_medium)
            rel_error = max(test_rel_error, rel_error)
        amp = amp_type(active_medium, count_z)
        amp.amplify(0.0, 0.0, pulse, count_t)
        if active_medium.doping_agent.lower_lifetime in analytical_lower_lifetimes:
            exact = pamp.amplifier.ExactOutputAmplifier(active_medium, count_z)
            exact_density_out, exact_population_out = exact.amplify(0.0, 0.0, pulse, count_t)
        else:
            exact_density_out, exact_population_out = None, None
        results = (amp, exact_density_out, exact_population_out)
        return results, rel_error
    
    min_count_z = max(min_count_z, 2)
    min_count_t = max(min_count_t, 3)
    compute_rdiff = lambda last_res, res: get_rel_error_amp(last_res[0], res[0])
    data = pamp.util.min_steps(min_count_z, min_count_t, True, True, params.amp_rtol, amplify_pulse, compute_rdiff, retextra=True)
    
    return data

def most_efficient_method(dirname, active_medium, beam_profile, ref_pulse, int_types, amp_types, quiet=False):
    if not quiet:
        print output.div_line
    if not quiet or params.verbose:
        print "determining most efficient method combination"
    
    method_results = []
    for int_type in int_types:
        int_name = int_type.__name__
        integrator = pamp.energy.PhotonCountIntegrator(int_type, active_medium, beam_profile)
        count_rho, count_phi, min_count_z, min_count_t = integrator.min_steps((ref_pulse, ), params.energy_rtol, params.fluence_rtol)
        count_rho = max(count_rho, mangle_count_cyl(params.min_count_rho))
        count_phi = max(count_phi, mangle_count_cyl(params.min_count_phi))
        min_count_z = max(min_count_z, mangle_count_z(params.min_count_z))
        min_count_t = max(min_count_t, mangle_count_t(params.min_count_t))
        for amp_type in amp_types:
            amp_name = amp_type.__name__
            if params.verbose:
                print "%s, %s" % (int_name, amp_name)
            test_min_count_z = max(min_count_z, amp_type.min_steps_z(active_medium))
            test_min_count_t = max(min_count_t, amp_type(active_medium, test_min_count_z).min_steps_t(ref_pulse))
            data = test_pulse(amp_type, active_medium, ref_pulse, (test_min_count_z, test_min_count_t), integrator)
            if data is None:
                continue
            (count_z, count_t), rel_error, results = data
            count_z, count_t = mangle_count_z(count_z), mangle_count_t(count_t)
            count = count_rho * count_phi * count_z * count_t
            method_results.append(((int_type, amp_type), (count_rho, count_phi, count_z, count_t), results, count, rel_error))
    if not method_results:
        raise ComputationError("no suitable method combination found")
    min_count = min([r[3] for r in method_results])
    min_rel_error = min([r[4] for r in method_results if r[3] == min_count])
    ref_result = [r[:3] for r in method_results if r[3] == min_count and r[4] == min_rel_error][0]
    (int_type, amp_type), (count_rho, count_phi, count_z, count_t), results = ref_result
    amp, exact_density_out, exact_population_out = results
    
    if not quiet:
        if params.graphs:
            integrator = pamp.energy.PhotonCountIntegrator(int_type, active_medium, beam_profile)
            fluences = np.empty(amp.count_z)
            for l in range(amp.count_z):
                fluences[l] = integrator.integral(amp.T, amp.density[l]) * active_medium.light_speed
            print "generating output"
            dirname = output.init_dir(dirname)
            output.plot_output(dirname, beam_profile, ref_pulse, params.pulse_duration, amp, fluences, exact_density_out, exact_population_out)
    
    return (int_type, amp_type), (count_rho, count_phi, count_z, count_t)

def setup_methods(dirname, (int_types, amp_types), ref_inversion, quiet=False):
    ref_pulse_dir = os.path.join(dirname, output.ref_pulse_rel_path)
    
    initial_inversion = pamp.inversion.UniformInversion(ref_inversion)
    
    doping_agent = pamp.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = pamp.medium.ActiveMedium(initial_inversion, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    
    pulse_photon_count = pamp.energy.photon_count(params.lasing_wavelen, params.pulse_energy)
    ref_fluence = params.beam_class.ref_fluence(params.beam_radius, pulse_photon_count)
    beam_profile = params.beam_class(params.beam_radius, ref_fluence)
    
    ref_pulse = create_pulse(active_medium, beam_profile, beam_profile.rho_ref, beam_profile.phi_ref)
    
    (int_type, amp_type), (count_rho, count_phi, count_z, count_t) = most_efficient_method(ref_pulse_dir, active_medium, beam_profile, ref_pulse, int_types, amp_types, quiet)
    
    if params.verbose:
        print "int_type: %s; amp_type: %s" % (int_type.__name__, amp_type.__name__, )
        print "count_rho: %d; count_phi: %d" % (count_rho, count_phi)
        print "count_z: %d; count_t: %d" % (count_z, count_t)
    
    return (int_type, amp_type), (count_rho, count_phi, count_z, count_t)

def amplify_train(dirname, num_types, counts, ref_inversion, quiet=False):
    if not quiet:
        print output.div_line
    if not quiet or params.verbose:
        print "amplifying pulse train"
    
    initial_inversion = pamp.inversion.UniformInversion(ref_inversion)
    
    doping_agent = pamp.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = pamp.medium.ActiveMedium(initial_inversion, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    
    pulse_photon_count = pamp.energy.photon_count(params.lasing_wavelen, params.pulse_energy)
    ref_fluence = params.beam_class.ref_fluence(params.beam_radius, pulse_photon_count)
    beam_profile = params.beam_class(params.beam_radius, ref_fluence)
    
    ref_pulse = create_pulse(active_medium, beam_profile, beam_profile.rho_ref, beam_profile.phi_ref)
    
    int_type, amp_type = num_types
    count_rho, count_phi, count_z, count_t = counts
    integrator = pamp.energy.PhotonCountIntegrator(int_type, active_medium, beam_profile)
    amp = amp_type(active_medium, count_z)
    
    Rho = np.linspace(0.0, active_medium.radius, count_rho)
    Phi = np.linspace(0.0, 2.0*math.pi, count_phi)
    
    output_fluence = np.empty((count_rho, count_phi))
    output_photon_counts = np.empty(params.train_pulse_count)
    
    empty_mdlist = [[None for n in range(count_phi)] for m in range(count_rho)]
    population = copy.deepcopy(empty_mdlist)
    pulses = copy.deepcopy(empty_mdlist)
    for m, rho in enumerate(Rho):
        for n, phi in enumerate(Phi):
            upper = np.vectorize(active_medium.initial_inversion.inversion)(rho, phi, amp.Z)
            lower = np.zeros(count_z)
            population[m][n] = (upper, lower)
            pulse = create_pulse(active_medium, beam_profile, rho, phi)
            pulses[m][n] = pulse
    
    decay = get_decay(active_medium, ref_pulse)
    
    pulse_num_stride = params.pulse_num_stride
    for pnum in range(params.train_pulse_count):
        if not quiet or params.verbose:
            output.show_status((pnum, None), (pulse_num_stride, None), False)
        
        for m, rho in enumerate(Rho):
            for n, phi in enumerate(Phi):
                pulse = pulses[m][n]
                density_out, population_out = amp.amplify(rho, phi, pulse, count_t, initial_population=population[m][n])
                
                upper = np.copy(population_out[0])
                lower = population_out[1] * decay
                population[m][n] = (upper, lower)
                
                fluence_out = integrator.integral(amp.T, density_out) * active_medium.light_speed
                output_fluence[m, n] = fluence_out
        
        if pnum == 0:
            ref_output_fluence = np.copy(output_fluence)
        output_photon_counts[pnum] = integrator.integrate_base(Rho, Phi, output_fluence)
    
    if not quiet or params.verbose:
        output.show_status((pnum+1, None), (pulse_num_stride, None), True)
    
    if not quiet or params.verbose:
        print "processing results"
    max_output_fluence = pamp.energy.energy(params.lasing_wavelen, np.amax(ref_output_fluence))
    train_output_photon_count = sum(reversed(output_photon_counts))
    train_output_energy = pamp.energy.energy(params.lasing_wavelen, train_output_photon_count)
    if not quiet:
        if params.graphs:
            print "generating output"
            ref_pulse_dir = os.path.join(dirname, output.ref_pulse_rel_path)
            dirname = output.init_dir(dirname)
            ref_pulse_dir = output.init_dir(ref_pulse_dir)
            output.plot_beam(ref_pulse_dir, beam_profile, Rho, Phi, ref_output_fluence)
            output.plot_train(dirname, beam_profile, active_medium, output_photon_counts)
    
    return max_output_fluence, train_output_energy

def validate():
    doping_agent = pamp.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = pamp.medium.ActiveMedium(None, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    
    pulse_photon_count = pamp.energy.photon_count(params.lasing_wavelen, params.pulse_energy)
    ref_fluence = params.beam_class.ref_fluence(params.beam_radius, pulse_photon_count)
    beam_profile = params.beam_class(params.beam_radius, ref_fluence)
    
    ref_pulse = create_pulse(active_medium, beam_profile, beam_profile.rho_ref, beam_profile.phi_ref)
    
    if not (params.train_pulse_period >= ref_pulse.duration or params.train_pulse_count == 1):
        raise ConfigurationError("invalid parameters: pulse period is less than (extended) pulse duration")
    
    train_duration = params.train_pulse_period * (params.train_pulse_count - 1) + params.pulse_duration
    if not train_duration <= min(params.pump_duration, params.dopant_upper_lifetime) / 10.0:
        warnings.warn("approximation validity condition violated: pulse train duration is not much shorter than than upper level lifetime and pump duration", stacklevel=2)

def compute_rel_error(ref_inversion, ref_inversion_rel_error):
    rel_error_inversion = math.exp(params.dopant_xsection * ref_inversion_rel_error * ref_inversion * params.medium_length) - 1.0
    rel_error_energy = params.pulse_rel_energy_trunc + params.amp_rtol + params.energy_rtol
    return rel_error_inversion + rel_error_energy

def report_output_characteristics(ref_inversion, max_output_fluence, output_energy, inversion_rel_error, energy_rel_error):
    print output.div_line
    print "results:"
    
    doping_agent = pamp.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = pamp.medium.ActiveMedium(None, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    
    pump_energy = params.pump_duration * params.pump_power
    stored_energy = pamp.energy.energy(params.lasing_wavelen, ref_inversion * active_medium.volume)
    
    pulse_photon_count = pamp.energy.photon_count(params.lasing_wavelen, params.pulse_energy)
    ref_fluence = params.beam_class.ref_fluence(params.beam_radius, pulse_photon_count)
    beam_profile = params.beam_class(params.beam_radius, ref_fluence)
    
    input_photon_count = beam_profile.fluence_integral(active_medium.radius)
    input_energy = pamp.energy.energy(params.lasing_wavelen, input_photon_count)
    input_energy *= params.train_pulse_count
    
    added_energy = output_energy - input_energy
    extraction_eff = added_energy / stored_energy
    total_eff = added_energy / pump_energy
    
    output_energy_abs_error = energy_rel_error * output_energy
    stored_energy_abs_error = inversion_rel_error * stored_energy
    
    extraction_eff_abs_error = (added_energy + output_energy_abs_error) / (stored_energy - stored_energy_abs_error) - extraction_eff
    total_eff_abs_error = output_energy_abs_error / pump_energy
    
    max_output_fluence_abs_error = max_output_fluence * energy_rel_error # rtol is specified for the energy, not density/fluence, but use it nevertheless
    
    unitconv.print_result("maximum output fluence [{}]: {} ~ {}", ("J/cm^2",), (max_output_fluence, max_output_fluence_abs_error))
    unitconv.print_result("total output energy [{}]: {} ~ {}", ("mJ",), (output_energy, output_energy_abs_error))
    unitconv.print_result("extraction efficiency [{}]: {} ~ {}", ("%",), (extraction_eff, extraction_eff_abs_error))
    unitconv.print_result("opt-opt efficiency [{}]: {} ~ {}", ("%",), (total_eff, total_eff_abs_error))