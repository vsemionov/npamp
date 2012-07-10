
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
import traceback

import numpy as np

import model

import params
import output
import unitconv


class ComputationError(Exception):
    pass


class PulseTrainMock(object):
    
    def __init__(self, pulse, count, period):
        self.pulse = pulse
        self.count = count
        self.period = period


def create_medium(ref_inversion):
    initial_inversion = None
    if ref_inversion is not None:
        initial_inversion = initial_inversion = model.inversion.UniformInversion(ref_inversion)
    doping_agent = model.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    medium = model.medium.ActiveMedium(initial_inversion, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    return medium

def create_beam():
    pulse_photon_count = model.energy.photon_count(params.lasing_wavelen, params.pulse_energy)
    ref_fluence = params.beam_class.ref_fluence(params.beam_radius, pulse_photon_count)
    beam = params.beam_class(params.beam_radius, ref_fluence)
    return beam

def create_pulse(active_medium, beam, rho, phi, ret_time_trunc_rel_error=False):
    fluence = beam.fluence(rho, phi)
    ref_density = params.pulse_class.ref_density(active_medium.light_speed, params.pulse_duration, fluence)
    pulse = params.pulse_class(-params.pulse_duration/2.0, params.pulse_duration, ref_density)
    scale, time_trunc_rel_error = model.error.pulse_scale(pulse, params.time_trunc_rtol)
    pulse = model.pulse.ExtendedPulse(pulse, scale)
    pulse = model.pulse.TruncatedPulse(pulse)
    return (pulse, time_trunc_rel_error) if ret_time_trunc_rel_error else pulse

def create_train(pulse):
    train = PulseTrainMock(pulse, params.train_pulse_count, params.train_pulse_period)
    return train

def create_depop_model(active_medium, depop_model_class):
    is_numerical = issubclass(depop_model_class, model.depop.NumericalDepopulationModel)
    if is_numerical:
        depop_model_kwargs = dict(rtol=params.depop_rate_rtol, min_samples=params.depop_rate_min_samples)
    else:
        depop_model_kwargs = {}
    if depop_model_class is params.depop_model_class:
        depop_model_kwargs.update(params.depop_model_extra_args)
    depop_model = depop_model_class(active_medium, params.lasing_wavelen, **depop_model_kwargs)
    return depop_model

def compute_inversion(dirname):
    print output.div_line
    print "computing population inversion"
    
    active_medium = create_medium(None)
    pump_system = model.pump.PumpSystem(params.pump_wavelen, params.pump_duration, params.pump_power, params.pump_efficiency)
    depop_model = create_depop_model(active_medium, params.depop_model_class)
    inv = params.inverter_class(active_medium, pump_system, depop_model)
    
    ref_inversion = inv.invert(params.inversion_rtol, params.inversion_min_count_t)
    rate_evals = (len(inv.inversion) - 1) * inv.evals_per_step
    pump_energy = params.pump_duration * params.pump_power
    stored_energy = model.energy.energy(params.lasing_wavelen, ref_inversion * active_medium.volume)
    
    if params.verbose:
        print "count_t:", len(inv.T)
        print "population rate evaluation count:", rate_evals
    
    if params.inversion_validate:
        print "validating uniform ASE-induced depopulation rate approximation"
        ross_num_model = depop_model if isinstance(depop_model, model.depop.RossNumericalASEModel) else model.depop.RossNumericalASEModel(active_medium, params.lasing_wavelen, params.depop_rate_rtol, params.depop_rate_min_samples)
        rate_rel_stddev = ross_num_model.rate_rel_stddev(ref_inversion)
        unitconv.print_result("depopulation rate rel. std. deviation [{}]: {}", ("%",), (rate_rel_stddev,))
        if rate_rel_stddev > 10.0e-2:
            output.warn("uniform ASE-induced depopulation rate approximation is invalid")
    
    if isinstance(depop_model, model.depop.NumericalDepopulationModel):
        print "perturbing population inversion"
        perturb_depop_model = model.depop.PerturbedDepopulationModel(depop_model)
        perturb_inv = params.inverter_class(active_medium, pump_system, perturb_depop_model)
        perturb_ref_inversion = perturb_inv.invert(params.inversion_rtol, params.inversion_min_count_t)
        rel_error = model.error.perturbed_inversion_rel_error(ref_inversion, perturb_ref_inversion, params.inversion_rtol)
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
        print output.status_writing
        dirname = output.init_dir(dirname)
        output.plot_inversion(dirname, inv)
    
    return ref_inversion, rel_error

def most_efficient_method((int_types, amp_types), active_medium, input_beam, ref_pulse, quiet=False):
    if not quiet:
        print output.div_line
    if not quiet or params.verbose:
        print "determining most efficient method combination"
    
    pulse_train = create_train(ref_pulse)
    
    size_exc_types = (ValueError, MemoryError)
    best_method = None
    for int_type in int_types:
        int_name = int_type.__name__
        if params.verbose:
            print int_name
        integrator = model.integrator.DomainIntegrator(int_type, active_medium)
        try:
            count_rho, count_phi, min_count_z, min_count_t = model.error.min_integration_steps(integrator, input_beam, (ref_pulse,), params.energy_rtol, params.fluence_rtol)
        except size_exc_types:
            output.print_exception()
            print "attempting to recover"
            sys.exc_clear()
            continue
        count_rho = max(count_rho, model.discrete.mangle_count_tv(params.min_count_rho))
        count_phi = max(count_phi, model.discrete.mangle_count_tv(params.min_count_phi))
        min_count_z = max(min_count_z, params.min_count_z)
        min_count_t = max(min_count_t, params.min_count_t)
        for amp_type in amp_types:
            amp_name = amp_type.__name__
            if params.verbose:
                print amp_name
            test_min_count_z = max(min_count_z, amp_type.min_steps_z(active_medium))
            test_min_count_t = max(min_count_t, amp_type(active_medium, test_min_count_z).min_steps_t(ref_pulse))
            try:
                data = model.error.min_amplification_steps(amp_type, active_medium, pulse_train, (test_min_count_z, test_min_count_t), integrator, params.fluence_rtol, params.amp_rtol, ret_extra=True)
            except size_exc_types:
                output.print_exception()
                print "attempting to recover"
                sys.exc_clear()
                continue
            if data is None:
                continue
            (count_z, count_t), rel_error, results = data
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
        raise ComputationError("no suitable numerical method combination found")
    (int_type, amp_type), (count_rho, count_phi, count_z, count_t), results, _, _ = best_method
    
    return (int_type, amp_type), (count_rho, count_phi, count_z, count_t)

def select_methods((int_types, amp_types), ref_inversion, ret_rel_errors=False, quiet=False):
    active_medium = create_medium(ref_inversion)
    input_beam = create_beam()
    ref_pulse, time_trunc_rel_error = create_pulse(active_medium, input_beam, input_beam.rho_ref, input_beam.phi_ref, ret_time_trunc_rel_error=True)
    
    (int_type, amp_type), (count_rho, count_phi, count_z, count_t) = most_efficient_method((int_types, amp_types), active_medium, input_beam, ref_pulse, quiet)
    
    if params.verbose:
        print "int_type: %s; amp_type: %s" % (int_type.__name__, amp_type.__name__, )
        print "count_rho: %d; count_phi: %d" % (count_rho, count_phi)
        print "count_z: %d; count_t: %d" % (count_z, count_t)
    
    numerics = (int_type, amp_type), (count_rho, count_phi, count_z, count_t)
    rel_errors = (time_trunc_rel_error, params.amp_rtol, params.energy_rtol)
    return (numerics, rel_errors) if ret_rel_errors else numerics

def amplify_ref_pulse(dirname, num_types, counts, ref_inversion):
    print output.div_line
    print "amplifying ref. pulse"
    
    dirname = os.path.join(dirname, output.ref_pulse_rel_path)
    
    (int_type, amp_type), (_, _, count_z, count_t) = num_types, counts
    
    active_medium = create_medium(ref_inversion)
    input_beam = create_beam()
    ref_pulse = create_pulse(active_medium, input_beam, input_beam.rho_ref, input_beam.phi_ref)
    
    integrator = model.integrator.DomainIntegrator(int_type, active_medium)
    amp = amp_type(active_medium, count_z)
    
    num_density_out, _ = amp.amplify(0.0, 0.0, ref_pulse, count_t)
    
    if active_medium.doping_agent.lower_lifetime in model.amplifier.ExactAmplifier.analytical_lower_lifetimes:
        exact = model.amplifier.ExactOutputAmplifier(active_medium, count_z)
        exact_density_out, exact_population_out = exact.amplify(0.0, 0.0, ref_pulse, count_t)
    else:
        exact_density_out, exact_population_out = None, None
    
    fluence_out = integrator.integrate(amp.T, num_density_out) * active_medium.light_speed
    fluence_gain = fluence_out / input_beam.ref_fluence
    unitconv.print_result("fluence gain: {}", (), (fluence_gain,))
    
    if params.graphs:
        integrator = model.integrator.DomainIntegrator(int_type, active_medium)
        fluences = np.empty(amp.count_z)
        for l in range(amp.count_z):
            fluences[l] = integrator.integrate(amp.T, amp.density[l]) * active_medium.light_speed
        print output.status_writing
        dirname = output.init_dir(dirname)
        output.plot_output(dirname, input_beam, ref_pulse, params.pulse_duration, amp, fluences, exact_density_out, exact_population_out)

def amplify_train(dirname, num_types, counts, ref_inversion, quiet=False):
    if not quiet:
        print output.div_line
    if not quiet or params.verbose:
        print "amplifying pulse train"
    
    active_medium = create_medium(ref_inversion)
    input_beam = create_beam()
    ref_pulse = create_pulse(active_medium, input_beam, input_beam.rho_ref, input_beam.phi_ref)
    pulse_train = create_train(ref_pulse)
    
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
    populations = copy.deepcopy(empty_mdlist)
    input_pulses = copy.deepcopy(empty_mdlist)
    for m, rho in enumerate(Rho):
        for n, phi in enumerate(Phi):
            upper = np.vectorize(active_medium.initial_inversion.inversion)(rho, phi, amp.Z)
            lower = np.zeros(count_z)
            populations[m][n] = (upper, lower)
            input_pulse = create_pulse(active_medium, input_beam, rho, phi)
            input_pulses[m][n] = input_pulse
    
    lower_decay = model.amplifier.lower_state_decay(active_medium, pulse_train)
    
    pulse_num_stride = params.pulse_num_stride
    for pnum in range(params.train_pulse_count):
        if not quiet or params.verbose:
            output.show_status((pnum, None), (pulse_num_stride, None), False)
        
        for m, rho in enumerate(Rho):
            for n, phi in enumerate(Phi):
                input_pulse = input_pulses[m][n]
                density_out, population_final = amp.amplify(rho, phi, input_pulse, count_t, initial_population=populations[m][n])
                
                upper = np.copy(population_final[0])
                lower = population_final[1] * lower_decay
                populations[m][n] = (upper, lower)
                
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
    rel_gain_decrease = 1.0 - output_photon_counts[-1] / output_photon_counts[0]
    
    if not quiet:
        if params.graphs:
            print output.status_writing
            ref_pulse_dir = os.path.join(dirname, output.ref_pulse_rel_path)
            dirname = output.init_dir(dirname)
            ref_pulse_dir = output.init_dir(ref_pulse_dir)
            output.plot_beam(ref_pulse_dir, input_beam, Rho, Phi, ref_output_fluence)
            output.plot_train(dirname, input_beam, active_medium, output_photon_counts)
    
    return max_output_fluence, output_photon_counts, train_output_energy, rel_gain_decrease

def report_results(ref_inversion, max_output_fluence, output_photon_counts, output_energy, rel_gain_decrease, inversion_rel_error, rel_errors):
    print output.div_line
    print "results:"
    
    active_medium = create_medium(ref_inversion)
    
    energy_rel_error = model.error.energy_rel_error(active_medium, inversion_rel_error, rel_errors)
    
    pump_energy = params.pump_duration * params.pump_power
    stored_energy = model.energy.energy(params.lasing_wavelen, ref_inversion * active_medium.volume)
    
    input_beam = create_beam()
    
    input_photon_count = input_beam.fluence_integral(active_medium.radius)
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
    rel_gain_decrease_abs_error = 0.0
    if params.train_pulse_count > 1:
        rel_gain_decrease_abs_error = (photon_count_last + photon_count_last_abs_error) / max(photon_count_first - photon_count_first_abs_error, 0.0) - photon_count_last / photon_count_first
    
    unitconv.print_result("input energy [{}]: {}", ("mJ",), (input_energy,))
    unitconv.print_result("output energy [{}]: {} ~ {}", ("mJ",), (output_energy, output_energy_abs_error))
    unitconv.print_result("energy gain: {} ~ {}", (), (energy_gain, energy_gain_abs_error))
    unitconv.print_result("extraction efficiency [{}]: {} ~ {}", ("%",), (extraction_eff, extraction_eff_abs_error))
    unitconv.print_result("opt.-opt. efficiency [{}]: {} ~ {}", ("%",), (total_eff, total_eff_abs_error))
    unitconv.print_result("max. output fluence [{}]: {} ~ {}", ("J/cm^2",), (max_output_fluence, max_output_fluence_abs_error))
    unitconv.print_result("rel. gain decrease [{}]: {} ~ {}", ("%",), (rel_gain_decrease, rel_gain_decrease_abs_error))
