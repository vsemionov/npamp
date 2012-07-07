
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

import math
import copy
import warnings
import traceback

import types

import numpy as np

import params
import model
import output
import unitconv


class ConfigurationError(Exception):
    pass

class ComputationError(Exception):
    pass


copy_conf = lambda conf: {k: v for k, v in conf.items() if not k.startswith('_') and type(v) is not types.ModuleType}

def create_pulse(active_medium, beam, rho, phi, ret_time_trunc_rel_error=False):
    fluence = beam.fluence(rho, phi)
    ref_density = params.pulse_class.ref_density(active_medium.light_speed, params.pulse_duration, fluence)
    pulse = params.pulse_class(-params.pulse_duration/2.0, params.pulse_duration, ref_density)
    scale, time_trunc_rel_error = model.error.pulse_scale(pulse, params.time_trunc_rtol)
    pulse = model.pulse.ExtendedPulse(pulse, scale)
    pulse = model.pulse.TruncatedPulse(pulse)
    return (pulse, time_trunc_rel_error) if ret_time_trunc_rel_error else pulse

def mangle_count_z(count_z):
    count_z = max(count_z, 3)
    count_z = model.util.steps(model.util.divs(count_z))
    return count_z

def mangle_count_t(count_t):
    count_t = max(count_t, 3)
    count_t = model.util.steps(model.util.divs(count_t))
    return count_t

def mangle_count_cyl(count_cyl):
    count_cyl = max(count_cyl, 1)
    if count_cyl > 1:
        count_cyl = max(count_cyl, 3)
        count_cyl = model.util.steps(model.util.divs(count_cyl))
    return count_cyl

def compute_inversion(dirname):
    print output.div_line
    print "computing population inversion"
    
    is_numerical = issubclass(params.depop_model_class, model.depop.NumericalDepopulationModel)
    doping_agent = model.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = model.medium.ActiveMedium(None, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    pump_system = model.pump.PumpSystem(params.pump_wavelen, params.pump_duration, params.pump_power, params.pump_efficiency)
    depop_model_kwargs = dict(rtol=params.depop_rate_rtol) if is_numerical else {}
    if is_numerical:
        depop_model_kwargs = dict(rtol=params.depop_rate_rtol, min_count=params.depop_rate_min_count)
    else:
        depop_model_kwargs = {}
    depop_model_kwargs.update(params.depop_model_extra_args)
    depop_model = params.depop_model_class(active_medium, params.lasing_wavelen, **depop_model_kwargs)
    inv = params.inverter_class(active_medium, pump_system, depop_model)
    ref_inversion = inv.invert(params.inversion_rtol, params.inversion_min_count_t)
    rate_evals = (len(inv.inversion) - 1) * inv.evals_per_step
    pump_energy = params.pump_duration * params.pump_power
    stored_energy = model.energy.energy(params.lasing_wavelen, ref_inversion * active_medium.volume)
    
    if params.verbose:
        print "count_t:", len(inv.T)
        print "number of depopulation rate evaluations:", rate_evals
    
    if params.inversion_validate:
        print "validating uniform ASE-induced depopulation rate approximation"
        ross_num_model = depop_model if isinstance(depop_model, model.depop.RossNumericalASEModel) else model.depop.RossNumericalASEModel(active_medium, params.lasing_wavelen, params.depop_rate_rtol, params.depop_rate_min_count)
        rate_rel_stddev = ross_num_model.rate_rel_stddev(ref_inversion)
        unitconv.print_result("depopulation rate rel. std. deviation [{}]: {}", ("%",), (rate_rel_stddev,))
        if rate_rel_stddev > 10.0e-2:
            warnings.warn("uniform ASE-induced depopulation rate approximation is invalid", stacklevel=2)
    
    if is_numerical:
        print "perturbing population inversion"
        perturb_depop_model = model.depop.PerturbedDepopulationModel(depop_model)
        perturb_inv = params.inverter_class(active_medium, pump_system, perturb_depop_model)
        perturb_ref_inversion = perturb_inv.invert(params.inversion_rtol, params.inversion_min_count_t)
        abs_error = abs(perturb_ref_inversion - ref_inversion) + (ref_inversion + perturb_ref_inversion) * params.inversion_rtol
        rel_error = abs_error / ref_inversion
    else:
        rel_error = params.inversion_rtol
    
    gain_coef = ref_inversion * active_medium.doping_agent.xsection
    gain = math.exp(gain_coef * active_medium.length)
    
    ref_inversion_atol = ref_inversion * rel_error
    gain_atol = gain * (math.exp(ref_inversion_atol * active_medium.doping_agent.xsection * active_medium.length) - 1.0)
    stored_energy_atol = stored_energy * rel_error
    
    unitconv.print_result("pump energy [{}]: {}", ("mJ",), (pump_energy,))
    unitconv.print_result("population inversion [{}]: {} ~ {}", ("cm^-3",), (ref_inversion, ref_inversion_atol))
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

def most_efficient_method(active_medium, beam_profile, ref_pulse, (int_types, amp_types), quiet=False):
    if not quiet:
        print output.div_line
    if not quiet or params.verbose:
        print "determining most efficient method combination"
    
    decay = get_decay(active_medium, ref_pulse)
    
    size_exc_types = (ValueError, MemoryError)
    best_method = None
    for int_type in int_types:
        int_name = int_type.__name__
        if params.verbose:
            print int_name
        integrator = model.integrator.DomainIntegrator(int_type, active_medium)
        try:
            count_rho, count_phi, min_count_z, min_count_t = model.error.min_integration_steps(integrator, beam_profile, (ref_pulse,), params.energy_rtol, params.fluence_rtol)
        except size_exc_types:
            traceback.print_exc()
            print "attempting to recover"
            sys.exc_clear()
            continue
        count_rho = max(count_rho, mangle_count_cyl(params.min_count_rho))
        count_phi = max(count_phi, mangle_count_cyl(params.min_count_phi))
        min_count_z = max(min_count_z, mangle_count_z(params.min_count_z))
        min_count_t = max(min_count_t, mangle_count_t(params.min_count_t))
        for amp_type in amp_types:
            amp_name = amp_type.__name__
            if params.verbose:
                print amp_name
            test_min_count_z = max(min_count_z, amp_type.min_steps_z(active_medium))
            test_min_count_t = max(min_count_t, amp_type(active_medium, test_min_count_z).min_steps_t(ref_pulse))
            try:
                data = model.error.min_amplification_steps(amp_type, active_medium, ref_pulse, decay, (test_min_count_z, test_min_count_t), integrator, params.fluence_rtol, params.amp_rtol, ret_extra=True)
            except size_exc_types:
                traceback.print_exc()
                print "attempting to recover"
                sys.exc_clear()
                continue
            if data is None:
                continue
            (count_z, count_t), rel_error, results = data
            count_z, count_t = mangle_count_z(count_z), mangle_count_t(count_t)
            count = count_rho * count_phi * count_z * count_t
            is_best = False
            if best_method is None:
                is_best = True
            else:
                _, _, _, best_count, best_rel_error = best_method
                if count < best_count:
                    is_best = True
                elif count == best_count:
                    if rel_error < best_rel_error:
                        is_best = True
            if is_best:
                best_method = ((int_type, amp_type), (count_rho, count_phi, count_z, count_t), results, count, rel_error)
    
    if best_method is None:
        raise ComputationError("no suitable method combination found")
    (int_type, amp_type), (count_rho, count_phi, count_z, count_t), results, _, _ = best_method
    
    return (int_type, amp_type), (count_rho, count_phi, count_z, count_t)

def setup_methods((int_types, amp_types), ref_inversion, ret_rel_errors=False, quiet=False):
    initial_inversion = model.inversion.UniformInversion(ref_inversion)
    
    doping_agent = model.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = model.medium.ActiveMedium(initial_inversion, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    
    pulse_photon_count = model.energy.photon_count(params.lasing_wavelen, params.pulse_energy)
    ref_fluence = params.beam_class.ref_fluence(params.beam_radius, pulse_photon_count)
    beam_profile = params.beam_class(params.beam_radius, ref_fluence)
    
    ref_pulse, time_trunc_rel_error = create_pulse(active_medium, beam_profile, beam_profile.rho_ref, beam_profile.phi_ref, ret_time_trunc_rel_error=True)
    
    (int_type, amp_type), (count_rho, count_phi, count_z, count_t) = most_efficient_method(active_medium, beam_profile, ref_pulse, (int_types, amp_types), quiet)
    
    if params.verbose:
        print "int_type: %s; amp_type: %s" % (int_type.__name__, amp_type.__name__, )
        print "count_rho: %d; count_phi: %d" % (count_rho, count_phi)
        print "count_z: %d; count_t: %d" % (count_z, count_t)
    
    numerics = (int_type, amp_type), (count_rho, count_phi, count_z, count_t)
    rel_errors = (time_trunc_rel_error,)
    return (numerics, rel_errors) if ret_rel_errors else numerics

def amplify_ref_pulse(dirname, num_types, counts, ref_inversion):
    print output.div_line
    print "amplifying ref. pulse"
    
    dirname = os.path.join(dirname, output.ref_pulse_rel_path)
    
    (int_type, amp_type), (_, _, count_z, count_t) = num_types, counts
    
    initial_inversion = model.inversion.UniformInversion(ref_inversion)
    
    doping_agent = model.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = model.medium.ActiveMedium(initial_inversion, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    
    pulse_photon_count = model.energy.photon_count(params.lasing_wavelen, params.pulse_energy)
    ref_fluence = params.beam_class.ref_fluence(params.beam_radius, pulse_photon_count)
    beam_profile = params.beam_class(params.beam_radius, ref_fluence)
    
    ref_pulse = create_pulse(active_medium, beam_profile, beam_profile.rho_ref, beam_profile.phi_ref)
    
    integrator = model.integrator.DomainIntegrator(int_type, active_medium)
    amp = amp_type(active_medium, count_z)
    
    num_density_out, _ = amp.amplify(0.0, 0.0, ref_pulse, count_t)
    
    if active_medium.doping_agent.lower_lifetime in model.amplifier.ExactAmplifier.analytical_lower_lifetimes:
        exact = model.amplifier.ExactOutputAmplifier(active_medium, count_z)
        exact_density_out, exact_population_out = exact.amplify(0.0, 0.0, ref_pulse, count_t)
    else:
        exact_density_out, exact_population_out = None, None
    
    fluence_out = integrator.integrate(amp.T, num_density_out) * active_medium.light_speed
    fluence_gain = fluence_out / ref_fluence
    unitconv.print_result("fluence gain: {}", (), (fluence_gain,))
    
    if params.graphs:
        integrator = model.integrator.DomainIntegrator(int_type, active_medium)
        fluences = np.empty(amp.count_z)
        for l in range(amp.count_z):
            fluences[l] = integrator.integrate(amp.T, amp.density[l]) * active_medium.light_speed
        print "generating output"
        dirname = output.init_dir(dirname)
        output.plot_output(dirname, beam_profile, ref_pulse, params.pulse_duration, amp, fluences, exact_density_out, exact_population_out)

def amplify_train(dirname, num_types, counts, ref_inversion, quiet=False):
    if not quiet:
        print output.div_line
    if not quiet or params.verbose:
        print "amplifying pulse train"
    
    initial_inversion = model.inversion.UniformInversion(ref_inversion)
    
    doping_agent = model.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = model.medium.ActiveMedium(initial_inversion, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    
    pulse_photon_count = model.energy.photon_count(params.lasing_wavelen, params.pulse_energy)
    ref_fluence = params.beam_class.ref_fluence(params.beam_radius, pulse_photon_count)
    beam_profile = params.beam_class(params.beam_radius, ref_fluence)
    
    ref_pulse = create_pulse(active_medium, beam_profile, beam_profile.rho_ref, beam_profile.phi_ref)
    
    int_type, amp_type = num_types
    count_rho, count_phi, count_z, count_t = counts
    integrator = model.integrator.DomainIntegrator(int_type, active_medium)
    amp = amp_type(active_medium, count_z)
    
    Rho = np.linspace(0.0, active_medium.radius, count_rho)
    Phi = np.linspace(0.0, 2.0*math.pi, count_phi)
    
    output_fluence = np.empty((count_rho, count_phi))
    max_fluences = np.empty(params.train_pulse_count)
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
                density_out, population_final = amp.amplify(rho, phi, pulse, count_t, initial_population=population[m][n])
                
                upper = np.copy(population_final[0])
                lower = population_final[1] * decay
                population[m][n] = (upper, lower)
                
                fluence_out = integrator.integrate(amp.T, density_out) * active_medium.light_speed
                output_fluence[m, n] = fluence_out
        
        if pnum == 0:
            ref_idx = np.unravel_index(output_fluence.argmax(), output_fluence.shape)
            ref_output_fluence = np.copy(output_fluence)
        max_fluences[pnum] = output_fluence[ref_idx]
        output_photon_counts[pnum] = integrator.integrate_base(Rho, Phi, output_fluence)
    
    if not quiet or params.verbose:
        output.show_status((pnum+1, None), (pulse_num_stride, None), True)
    
    if not quiet or params.verbose:
        print "processing results"
    max_output_fluence = max_fluences[::-1].sum()
    max_output_fluence = model.energy.energy(params.lasing_wavelen, max_output_fluence)
    train_output_photon_count = output_photon_counts[::-1].sum()
    train_output_energy = model.energy.energy(params.lasing_wavelen, train_output_photon_count)
    rel_gain_reduction = 1.0 - output_photon_counts[-1] / output_photon_counts[0]
    if not quiet:
        if params.graphs:
            print "generating output"
            ref_pulse_dir = os.path.join(dirname, output.ref_pulse_rel_path)
            dirname = output.init_dir(dirname)
            ref_pulse_dir = output.init_dir(ref_pulse_dir)
            output.plot_beam(ref_pulse_dir, beam_profile, Rho, Phi, ref_output_fluence)
            output.plot_train(dirname, beam_profile, active_medium, output_photon_counts)
    
    return max_output_fluence, output_photon_counts, train_output_energy, rel_gain_reduction

def validate():
    def validate_param_bounds(name, value):
        vtype = type(value)
        if vtype in (int, float):
            if vtype is float and name not in param_zero_allowed:
                if value <= 0.0:
                    raise ConfigurationError("parameter \"%s\" has (or contains) a zero value" % name)
            elif value < 0.0:
                raise ConfigurationError("parameter \"%s\" has (or contains) a negative value" % name)
            if vtype is float and name not in param_infty_allowed:
                if math.isinf(value):
                    raise ConfigurationError("parameter \"%s\" has (or contains) an infinite value" % name)
            if name in param_bounds:
                min_value, max_value = param_bounds[name]
                if min_value is not None:
                    if value < min_value:
                        raise ConfigurationError("parameter \"%s\" has (or contains) a value less than %s" % (name, min_value))
                if max_value is not None:
                    if value > max_value:
                        raise ConfigurationError("parameter \"%s\" has (or contains) a value greater than %s" % (name, max_value))
        elif vtype in (tuple, list):
            for element in value:
                validate_param_bounds(name, element)
    
    def validate_interval(name, values):
        if values[0] > values[1]:
            raise ConfigurationError("parameter \"%s\" has a lower bound greater than the upper bound" % name)
    
    def validate_list_uniques(name, values):
        if len(values) != len(set(values)):
            raise ConfigurationError("parameter \"%s\" contains duplicate values" % name)
    
    param_zero_allowed = set([
        "dopant_branching_ratio",
        "dopant_lower_lifetime",
        "initial_inversion",
    ])
    param_infty_allowed = set([
        "dopant_upper_lifetime",
        "dopant_lower_lifetime",
        "ext_opt_inversion_rdiff_max",
        "ext_opt_fluence_max",
    ])
    param_bounds = {
        "train_pulse_count": (1, None),
        "depop_rate_min_count": (16, None),
        "out_markers_step_divisor": (1, None),
        "out_rho_steps_divisor": (1, None),
        "out_phi_steps_divisor": (1, None),
        "out_z_steps_divisor": (1, None),
        "out_t_steps_divisor": (1, None),
        "ext_opt_pump_resolution": (2, None),
        "ext_opt_geom_resolution": (2, None),
        "dopant_branching_ratio": (None, 1.0),
        "pump_efficiency": (None, 1.0),
    }
    param_intervals = set([
        "ext_opt_pump_duration",
        "ext_opt_pump_power",
        "ext_opt_geom_mediumradius",
        "ext_opt_geom_beamradius",
    ])
    
    conf = copy_conf(params.__dict__)
    
    for parameter, value in conf.items():
        validate_param_bounds(parameter, value)
    
    for parameter in param_intervals:
        values = conf[parameter]
        validate_interval(parameter, values)
    
    doping_agent = model.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = model.medium.ActiveMedium(None, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    
    pulse_photon_count = model.energy.photon_count(params.lasing_wavelen, params.pulse_energy)
    ref_fluence = params.beam_class.ref_fluence(params.beam_radius, pulse_photon_count)
    beam_profile = params.beam_class(params.beam_radius, ref_fluence)
    
    ref_pulse = create_pulse(active_medium, beam_profile, beam_profile.rho_ref, beam_profile.phi_ref)
    
    if not (params.train_pulse_period >= ref_pulse.duration or params.train_pulse_count == 1):
        raise ConfigurationError("invalid parameters: pulse repetition period is less than (extended) pulse duration")
    
    if not params.pump_wavelen < params.lasing_wavelen:
        warnings.warn("pump wavelength is not less than lasing wavelength", stacklevel=2)
    
    if not params.dopant_lower_lifetime <= params.dopant_upper_lifetime / 10.0:
        warnings.warn("approximation validity condition violated: lower state lifetime is not much shorter than upper state (fluorescence) lifetime", stacklevel=2)
    
    train_duration = params.train_pulse_period * (params.train_pulse_count - 1) + params.pulse_duration
    if not train_duration <= params.dopant_upper_lifetime / 10.0:
        warnings.warn("approximation validity condition violated: pulse train duration is not much shorter than upper state (fluorescence) lifetime", stacklevel=2)

def compute_energy_rel_error(ref_inversion, ref_inversion_rel_error, rel_errors):
    (time_trunc_rel_error,) = rel_errors
    rel_error_inversion = math.exp(params.dopant_xsection * ref_inversion_rel_error * ref_inversion * params.medium_length) - 1.0
    rel_error_energy = time_trunc_rel_error + params.amp_rtol + params.energy_rtol
    energy_rel_error = rel_error_inversion + rel_error_energy
    return energy_rel_error

def report_output_characteristics(ref_inversion, max_output_fluence, output_photon_counts, output_energy, rel_gain_reduction, inversion_rel_error, energy_rel_error):
    print output.div_line
    print "results:"
    
    doping_agent = model.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = model.medium.ActiveMedium(None, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    
    pump_energy = params.pump_duration * params.pump_power
    stored_energy = model.energy.energy(params.lasing_wavelen, ref_inversion * active_medium.volume)
    
    pulse_photon_count = model.energy.photon_count(params.lasing_wavelen, params.pulse_energy)
    ref_fluence = params.beam_class.ref_fluence(params.beam_radius, pulse_photon_count)
    beam_profile = params.beam_class(params.beam_radius, ref_fluence)
    
    input_photon_count = beam_profile.fluence_integral(active_medium.radius)
    input_energy = model.energy.energy(params.lasing_wavelen, input_photon_count)
    input_energy *= params.train_pulse_count
    
    energy_gain = output_energy / input_energy
    added_energy = output_energy - input_energy
    extraction_eff = added_energy / stored_energy
    total_eff = added_energy / pump_energy
    
    stored_energy_abs_error = inversion_rel_error * stored_energy
    output_energy_abs_error = energy_rel_error * output_energy
    energy_gain_abs_error = energy_rel_error * energy_gain
    
    extraction_eff_abs_error = (added_energy + output_energy_abs_error) / max(stored_energy - stored_energy_abs_error, 0.0) - extraction_eff
    total_eff_abs_error = output_energy_abs_error / pump_energy
    
    max_output_fluence_abs_error = max_output_fluence * energy_rel_error # rtol is specified for the energy, not density/fluence, but use it nevertheless
    
    photon_count_first, photon_count_last = output_photon_counts[0], output_photon_counts[-1]
    photon_count_first_abs_error, photon_count_last_abs_error = photon_count_first * energy_rel_error, photon_count_last * energy_rel_error
    rel_gain_reduction_abs_error = 0.0
    if params.train_pulse_count > 1:
        rel_gain_reduction_abs_error = (photon_count_last + photon_count_last_abs_error) / max(photon_count_first - photon_count_first_abs_error, 0.0) - photon_count_last / photon_count_first
    
    unitconv.print_result("input energy [{}]: {}", ("mJ",), (input_energy,))
    unitconv.print_result("output energy [{}]: {} ~ {}", ("mJ",), (output_energy, output_energy_abs_error))
    unitconv.print_result("energy gain: {} ~ {}", (), (energy_gain, energy_gain_abs_error))
    unitconv.print_result("extraction efficiency [{}]: {} ~ {}", ("%",), (extraction_eff, extraction_eff_abs_error))
    unitconv.print_result("opt-opt efficiency [{}]: {} ~ {}", ("%",), (total_eff, total_eff_abs_error))
    unitconv.print_result("max output fluence [{}]: {} ~ {}", ("J/cm^2",), (max_output_fluence, max_output_fluence_abs_error))
    unitconv.print_result("rel gain reduction [{}]: {} ~ {}", ("%",), (rel_gain_reduction, rel_gain_reduction_abs_error))
