
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

import numpy as np

import params
import model
import plot
import core
import output
import unitconv


def compare_depop_models(dirname):
    filename = lambda name: os.path.join(dirname, name)
    
    print output.div_line
    print "comparing depopulation models"
    
    if not params.ext_depop_models:
        return
    
    doping_agent = model.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = model.medium.ActiveMedium(None, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    pump_system = model.pump.PumpSystem(params.pump_wavelen, params.pump_duration, params.pump_power, params.pump_efficiency)
    data = []
    for depop_model_class in params.ext_depop_models:
        depop_model_label = depop_model_class.descr
        print depop_model_label
        depop_model_kwargs = dict(rtol=params.depop_rate_rtol, min_count=params.depop_rate_min_count) if issubclass(depop_model_class, model.depop.NumericalDepopulationModel) else {}
        if depop_model_class is params.depop_model_class:
            depop_model_kwargs.update(params.depop_model_extra_args)
        depop_model = depop_model_class(active_medium, params.lasing_wavelen, **depop_model_kwargs)
        inv = params.inverter_class(active_medium, pump_system, depop_model)
        inv.invert(params.inversion_rtol, params.inversion_min_count_t)
        depop_rate = np.vectorize(depop_model.rate)(inv.inversion) / active_medium.volume
        data.append((inv.T, inv.inversion, depop_rate, depop_model_label, depop_model_class))
        ref_inversion = inv.inversion[-1]
        unitconv.print_result("population inversion [{}]: {}", ("cm^-3",), (ref_inversion,))
        unitconv.print_result("depopulation rate [{}]: {}", ("cm^-3 s^-1",), (depop_rate[-1],))
    
    if params.graphs:
        print "generating output"
        dirname = os.path.join(dirname, output.models_rel_path)
        dirname = output.init_dir(dirname)
        data.sort(key = lambda x: x[1][-1], reverse=True)
        Ts, inversions, depop_rates, labels, depop_model_classes = zip(*data)
        plot.plot_data(filename("inversions_evo"), "Population Inversion Evolution", (Ts, None, None, output.t_pump_label), (inversions, None, None, output.inversion_abs_label), labels)
        pump_rate = pump_system.effective_pump_rate / active_medium.volume
        abs_rate_ylim = None #(0.0, pump_rate * 1.25)
        non_zero_Ts = [T[1:] for T in Ts]
        non_zero_inversions = [inversion[1:] for inversion in inversions]
        non_zero_rates = [depop_rate[1:] for depop_rate in depop_rates]
        rel_depop_rates = [depop_rate / inversion for depop_rate, inversion in zip(non_zero_rates, non_zero_inversions)]
        plot.plot_data(filename("depop_rates"), "Depopulation Rate", (inversions, None, None, output.inversion_abs_label), (depop_rates, None, abs_rate_ylim, output.rate_label), labels, yvals=[(pump_rate, "pump rate")])
        plot.plot_data(filename("depop_rates_alt"), "Depopulation Rate to Inversion Ratio", (non_zero_inversions, None, None, output.inversion_abs_label), (rel_depop_rates, None, None, output.rate_rel_label), labels)
        plot.plot_data(filename("depop_rates_evo"), "Depopulation Rate Evolution", (Ts, None, None, output.t_pump_label), (depop_rates, None, abs_rate_ylim, output.rate_label), labels, yvals=[(pump_rate, "pump rate")])
        plot.plot_data(filename("depop_rates_alt_evo"), "Depopulation Rate to Inversion Ratio Evolution", (non_zero_Ts, None, None, output.t_pump_label), (rel_depop_rates, None, None, output.rate_rel_label), labels)
        
        if params.ext_alt_depop_model not in depop_model_classes:
            return
        alt_model_idx = depop_model_classes.index(params.ext_alt_depop_model)
        alt_T = Ts[alt_model_idx]
        alt_inversion = inversions[alt_model_idx]
        altinvs = []
        aTs = []
        altinv_inversion_rdiffs = []
        for cls, T, inversion in zip(depop_model_classes, Ts, inversions):
            if cls is params.ext_alt_depop_model:
                continue
            uT = set(list(T) + list(alt_T))
            aT = np.array(sorted(list(uT)))
            altinv = np.interp(aT, alt_T, alt_inversion)[1:]
            inv = np.interp(aT, T, inversion)[1:]
            rdiff = np.fabs(inv - altinv) / np.fmin(inv, altinv)
            aTs.append(aT[1:])
            altinvs.append(altinv)
            altinv_inversion_rdiffs.append(rdiff)
        non_alt_labels = [label for i, label in enumerate(labels) if i != alt_model_idx]
        plot.plot_data(filename("inversions_rdiff_inv"), "Inversion Relative Difference", (altinvs, None, None, output.inversion_abs_label), (altinv_inversion_rdiffs, None, None, output.inversion_rdiff_label), non_alt_labels)
        plot.plot_data(filename("inversions_rdiff_evo"), "Inversion Relative Difference Evolution", (aTs, None, None, output.t_pump_label), (altinv_inversion_rdiffs, None, None, output.inversion_rdiff_label), non_alt_labels)

def compute_inversion_pump_dependence_task((i, j), (tau, pwr), active_medium, (depop_model1, depop_model2)):
    output.show_status((i, j), params.extended_status_strides, False)
    pump_system = model.pump.PumpSystem(params.pump_wavelen, tau, pwr, params.pump_efficiency)
    inv1 = params.inverter_class(active_medium, pump_system, depop_model1)
    inv2 = params.inverter_class(active_medium, pump_system, depop_model2)
    inversion1 = inv1.invert(params.inversion_rtol, params.inversion_min_count_t)
    inversion2 = inv2.invert(params.inversion_rtol, params.inversion_min_count_t)
    gain_coef1 = inversion1 * active_medium.doping_agent.xsection
    gain_coef2 = inversion2 * active_medium.doping_agent.xsection
    gain1 = math.exp(gain_coef1 * active_medium.length)
    gain2 = math.exp(gain_coef2 * active_medium.length)
    inversion_rdiff = abs((inversion1 - inversion2) / min(inversion1, inversion2))
    gain_rdiff = abs((gain1 - gain2) / min(gain1, gain2))
    return inversion1, inversion2, gain1, gain2, inversion_rdiff, gain_rdiff

def compute_inversion_pump_dependence(task_pool, dirname):
    filename = lambda name: os.path.join(dirname, name)
    
    print output.div_line
    print "computing inversion dependence on pumping parameters"
    
    doping_agent = model.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = model.medium.ActiveMedium(None, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    
    count_tau = params.ext_opt_pump_resolution[0]
    count_pwr = params.ext_opt_pump_resolution[1]
    
    Tau = np.linspace(params.ext_opt_pump_duration[0], params.ext_opt_pump_duration[1], count_tau)
    Pwr = np.linspace(params.ext_opt_pump_power[0], params.ext_opt_pump_power[1], count_pwr)
    
    num_model_kwargs = dict(rtol=params.depop_rate_rtol, min_count=params.depop_rate_min_count)
    depop_model_class1 = params.depop_model_class
    depop_model_class2 = params.ext_alt_depop_model
    depop_model_kwargs1 = num_model_kwargs if issubclass(depop_model_class1, model.depop.NumericalDepopulationModel) else {}
    depop_model_kwargs2 = num_model_kwargs if issubclass(depop_model_class2, model.depop.NumericalDepopulationModel) else {}
    if depop_model_class1 is params.depop_model_class:
        depop_model_kwargs1.update(params.depop_model_extra_args)
    if depop_model_class2 is params.depop_model_class:
        depop_model_kwargs2.update(params.depop_model_extra_args)
    depop_model1 = depop_model_class1(active_medium, params.lasing_wavelen, **depop_model_kwargs1)
    depop_model2 = depop_model_class2(active_medium, params.lasing_wavelen, **depop_model_kwargs2)
    
    inversions1, inversions2, gains1, gains2, inversion_rdiffs, gain_rdiffs = task_pool.parallel_task(compute_inversion_pump_dependence_task, (Tau, Pwr), (), (active_medium, (depop_model1, depop_model2)))
    
    output.show_status((count_tau, count_pwr), params.extended_status_strides, True)
    
    pump_energies = np.prod(np.array(np.meshgrid(Tau, Pwr)), axis=0).T
    stored_energies1 = model.energy.energy(params.lasing_wavelen, inversions1 * active_medium.volume)
    stored_energies2 = model.energy.energy(params.lasing_wavelen, inversions2 * active_medium.volume)
    
    if params.graphs:
        print "generating output"
        dirname = os.path.join(dirname, output.opt_pump_rel_path)
        inversion_rdiff_max = params.ext_opt_inversion_rdiff_max
        zlim = None #(0.0, inversion_rdiff_max)
        depop_contours = [inversion_rdiff_max]
        ref_pump_energy = params.pump_duration * params.pump_power
        ref_pump_contours = [(pump_energies.T, ref_pump_energy, "const. pump energy")]
        graph_types = [
            (dirname, Pwr, output.pump_power_label),
            (os.path.join(dirname, output.alt_plot_rel_path), Pwr * params.pump_efficiency / active_medium.volume, output.eff_power_density_label),
        ]
        for dirname, Y, ylabel in graph_types:
            dirname = output.init_dir(dirname)
            plot.plot_color(filename("energy_pump"), "Pump Energy", (Tau, None, None, output. pump_duration_label), (Y, None, None, ylabel), (pump_energies.T, None, output.energy_abs_pump_label), params.num_auto_contours)
            plot.plot_color(filename("inversion"), "Inversion (%s)" % depop_model_class1.descr, (Tau, None, None, output.pump_duration_label), (Y, None, None, ylabel), (inversions1.T, None, output.inversion_abs_label), params.num_auto_contours, extra_contours=ref_pump_contours)
            plot.plot_color(filename("inversion_alt"), "Inversion (%s)" % depop_model_class2.descr, (Tau, None, None, output.pump_duration_label), (Y, None, None, ylabel), (inversions2.T, None, output.inversion_abs_label), params.num_auto_contours, extra_contours=ref_pump_contours)
            plot.plot_color(filename("ss_gain"), "Small Signal Gain (%s)" % depop_model_class1.descr, (Tau, None, None, output.pump_duration_label), (Y, None, None, ylabel), (gains1.T, None, output.gain_label), params.num_auto_contours, extra_contours=ref_pump_contours)
            plot.plot_color(filename("ss_gain_alt"), "Small Signal Gain (%s)" % depop_model_class2.descr, (Tau, None, None, output.pump_duration_label), (Y, None, None, ylabel), (gains2.T, None, output.gain_label), params.num_auto_contours, extra_contours=ref_pump_contours)
            plot.plot_color(filename("energy_stored"), "Stored Energy (%s)" % depop_model_class1.descr, (Tau, None, None, output.pump_duration_label), (Y, None, None, ylabel), (stored_energies1.T, None, output.energy_abs_stored_label), params.num_auto_contours, extra_contours=ref_pump_contours)
            plot.plot_color(filename("energy_stored_alt"), "Stored Energy (%s)" % depop_model_class2.descr, (Tau, None, None, output.pump_duration_label), (Y, None, None, ylabel), (stored_energies2.T, None, output.energy_abs_stored_label), params.num_auto_contours, extra_contours=ref_pump_contours)
            plot.plot_color(filename("inversion_rdiff"), "Inversion Relative Difference", (Tau, None, None, output.pump_duration_label), (Y, None, None, ylabel), (inversion_rdiffs.T, zlim, output.inversion_rdiff_label), params.num_auto_contours, depop_contours)
            plot.plot_color(filename("ss_gain_rdiff"), "Small Signal Gain Relative Difference", (Tau, None, None, output.pump_duration_label), (Y, None, None, ylabel), (gain_rdiffs.T, zlim, output.gain_rdiff_label), params.num_auto_contours, depop_contours)
    
    return inversions1, inversion_rdiffs

def compute_inversion_geom_dependence_task(i, rm, doping_agent, (depop_model_class1, depop_model_class2), (depop_model_kwargs1, depop_model_kwargs2)):
    output.show_status((i, None), params.extended_status_strides, False)
    medium_radius_orig = params.medium_radius
    params.medium_radius = rm
    active_medium = model.medium.ActiveMedium(None, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    depop_model1 = depop_model_class1(active_medium, params.lasing_wavelen, **depop_model_kwargs1)
    depop_model2 = depop_model_class2(active_medium, params.lasing_wavelen, **depop_model_kwargs2)
    pump_system = model.pump.PumpSystem(params.pump_wavelen, params.pump_duration, params.pump_power, params.pump_efficiency)
    inv1 = params.inverter_class(active_medium, pump_system, depop_model1)
    inv2 = params.inverter_class(active_medium, pump_system, depop_model2)
    inversion1 = inv1.invert(params.inversion_rtol, params.inversion_min_count_t)
    inversion2 = inv2.invert(params.inversion_rtol, params.inversion_min_count_t)
    gain_coef1 = inversion1 * active_medium.doping_agent.xsection
    gain_coef2 = inversion2 * active_medium.doping_agent.xsection
    gain1 = math.exp(gain_coef1 * active_medium.length)
    gain2 = math.exp(gain_coef2 * active_medium.length)
    stored_energy1 = model.energy.energy(params.lasing_wavelen, inversion1 * active_medium.volume)
    stored_energy2 = model.energy.energy(params.lasing_wavelen, inversion2 * active_medium.volume)
    inversion_rdiff = abs((inversion1 - inversion2) / min(inversion1, inversion2))
    gain_rdiff = abs((gain1 - gain2) / min(gain1, gain2))
    params.medium_radius = medium_radius_orig
    return inversion1, inversion2, gain1, gain2, stored_energy1, stored_energy2, inversion_rdiff, gain_rdiff

def compute_inversion_geom_dependence(task_pool, dirname):
    filename = lambda name: os.path.join(dirname, name)
    
    print output.div_line
    print "computing inversion dependence on geometry parameters"
    
    doping_agent = model.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    
    min_medium_radius = params.ext_opt_geom_mediumradius[0]
    
    count_rm = params.ext_opt_geom_resolution[0]
    
    Rm = np.linspace(min_medium_radius, params.ext_opt_geom_mediumradius[1], count_rm)
    
    num_model_kwargs = dict(rtol=params.depop_rate_rtol, min_count=params.depop_rate_min_count)
    depop_model_class1 = params.depop_model_class
    depop_model_class2 = params.ext_alt_depop_model
    depop_model_kwargs1 = num_model_kwargs if issubclass(depop_model_class1, model.depop.NumericalDepopulationModel) else {}
    depop_model_kwargs2 = num_model_kwargs if issubclass(depop_model_class2, model.depop.NumericalDepopulationModel) else {}
    if depop_model_class1 is params.depop_model_class:
        depop_model_kwargs1.update(params.depop_model_extra_args)
    if depop_model_class2 is params.depop_model_class:
        depop_model_kwargs2.update(params.depop_model_extra_args)
    
    inversions1, inversions2, gains1, gains2, stored_energies1, stored_energies2, inversion_rdiffs, gain_rdiffs = task_pool.parallel_task(compute_inversion_geom_dependence_task, (Rm,), (), (doping_agent, (depop_model_class1, depop_model_class2), (depop_model_kwargs1, depop_model_kwargs2)))
    
    output.show_status((count_rm, None), params.extended_status_strides, True)
    
    if params.graphs:
        print "generating output"
        dirname = os.path.join(dirname, output.opt_geom_rel_path)
        dirname = output.init_dir(dirname)
        inversion_rdiff_max = params.ext_opt_inversion_rdiff_max
        pump_energy = params.pump_duration * params.pump_power
        energy_ylim = (0.0, pump_energy * 1.25)
        labels = [cls.descr for cls in [depop_model_class1, depop_model_class2]]
        plot.plot_data(filename("inversion"), "Inversion (%s)" % depop_model_class1.descr, (Rm, None, None, output.medium_radius_label), (inversions1, None, None, output.inversion_abs_label))
        plot.plot_data(filename("inversion_alt"), "Inversion (%s)" % depop_model_class2.descr, (Rm, None, None, output.medium_radius_label), (inversions2, None, None, output.inversion_abs_label))
        plot.plot_data(filename("inversions"), "Population Inversion", ([Rm]*2, None, None, output.medium_radius_label), ([inversions1, inversions2], None, None, output.inversion_abs_label), legend=labels)
        plot.plot_data(filename("ss_gain"), "Small Signal Gain (%s)" % depop_model_class1.descr, (Rm, None, None, output.medium_radius_label), (gains1, None, None, output.gain_label))
        plot.plot_data(filename("ss_gain_alt"), "Small Signal Gain (%s)" % depop_model_class2.descr, (Rm, None, None, output.medium_radius_label), (gains2, None, None, output.gain_label))
        plot.plot_data(filename("ss_gains"), "Small Signal Gain", ([Rm]*2, None, None, output.medium_radius_label), ([gains1, gains2], None, None, output.gain_label), legend=labels)
        plot.plot_data(filename("energy_stored"), "Stored Energy (%s)" % depop_model_class1.descr, (Rm, None, None, output.medium_radius_label), (stored_energies1, None, energy_ylim, output.energy_abs_stored_label), yvals=[(pump_energy, "pump energy")])
        plot.plot_data(filename("energy_stored_alt"), "Stored Energy (%s)" % depop_model_class2.descr, (Rm, None, None, output.medium_radius_label), (stored_energies2, None, energy_ylim, output.energy_abs_stored_label), yvals=[(pump_energy, "pump energy")])
        plot.plot_data(filename("energies_stored"), "Stored Energy", ([Rm]*2, None, None, output.medium_radius_label), ([stored_energies1, stored_energies2], None, energy_ylim, output.energy_abs_stored_label), legend=labels, yvals=[(pump_energy, "pump energy")])
        plot.plot_data(filename("inversion_rdiff"), "Inversion Relative Difference", (Rm, None, None, output.medium_radius_label), (inversion_rdiffs, None, None, output.inversion_rdiff_label), yvals=[(inversion_rdiff_max, None)])
        plot.plot_data(filename("ss_gain_rdiff"), "Small Signal Gain Relative Difference", (Rm, None, None, output.medium_radius_label), (gain_rdiffs, None, None, output.gain_rdiff_label), yvals=[(inversion_rdiff_max, None)])
    
    return inversions1, inversion_rdiffs

def compute_fluence_pump_dependence_task((i, j), (tau, pwr), inversion, active_medium, (rho, phi), integrator, amp, (count_z, count_t), ref_pulse, decay):
    output.show_status((i, j), params.extended_status_strides, False)
    
    initial_inversion = model.inversion.UniformInversion(inversion)
    active_medium.initial_inversion = initial_inversion
    
    upper = np.vectorize(initial_inversion.inversion)(rho, phi, amp.Z)
    lower = np.zeros(count_z)
    population = (upper, lower)
    
    pulse_fluences = np.empty(params.train_pulse_count)
    
    for pnum in range(params.train_pulse_count):
        density_out, population_out = amp.amplify(rho, phi, ref_pulse, count_t, initial_population=population)
        
        upper = np.copy(population_out[0])
        lower = population_out[1] * decay
        population = (upper, lower)
        
        fluence_out = integrator.integral(amp.T, density_out) * active_medium.light_speed
        pulse_fluences[pnum] = fluence_out
    
    fluence_out = pulse_fluences[::-1].sum()
    fluence_out = model.energy.energy(params.lasing_wavelen, fluence_out)
    return fluence_out

def compute_fluence_pump_dependence(task_pool, dirname, (int_types, amp_types), inversions, (int_type, amp_type), (count_z, count_t)):
    filename = lambda name: os.path.join(dirname, name)
    
    print output.div_line
    print "computing fluence dependence on pumping parameters"
    
    doping_agent = model.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = model.medium.ActiveMedium(None, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    
    pulse_photon_count = model.energy.photon_count(params.lasing_wavelen, params.pulse_energy)
    ref_fluence = params.beam_class.ref_fluence(params.beam_radius, pulse_photon_count)
    beam_profile = params.beam_class(params.beam_radius, ref_fluence)
    
    ref_pulse = core.create_pulse(active_medium, beam_profile, beam_profile.rho_ref, beam_profile.phi_ref)
    
    decay = core.get_decay(active_medium, ref_pulse)
    
    integrator = model.energy.PhotonCountIntegrator(int_type, active_medium, beam_profile)
    amp = amp_type(active_medium, count_z)
    
    count_tau = params.ext_opt_pump_resolution[0]
    count_pwr = params.ext_opt_pump_resolution[1]
    
    Tau = np.linspace(params.ext_opt_pump_duration[0], params.ext_opt_pump_duration[1], count_tau)
    Pwr = np.linspace(params.ext_opt_pump_power[0], params.ext_opt_pump_power[1], count_pwr)
    
    rho, phi = beam_profile.rho_ref, beam_profile.phi_ref
    fluences = task_pool.parallel_task(compute_fluence_pump_dependence_task, (Tau, Pwr), (inversions,), (active_medium, (rho, phi), integrator, amp, (count_z, count_t), ref_pulse, decay))
    
    output.show_status((count_tau, count_pwr), params.extended_status_strides, True)
    
    if params.graphs:
        print "generating output"
        dirname = os.path.join(dirname, output.opt_pump_rel_path)
        fluence_max = params.ext_opt_fluence_max
        zlim = None #(0.0, fluence_max)
        contours = [fluence_max]
        graph_types = [
            (dirname, Pwr, output.pump_power_label),
            (os.path.join(dirname, output.alt_plot_rel_path), Pwr * params.pump_efficiency / active_medium.volume, output.eff_power_density_label),
        ]
        for dirname, Y, ylabel in graph_types:
            dirname = output.init_dir(dirname)
            plot.plot_color(filename("fluence_out_max"), "Maximum Output Fluence", (Tau, None, None, output.pump_duration_label), (Y, None, None, ylabel), (fluences.T, zlim, output.fluence_abs_label_energy), params.num_auto_contours, contours)
    
    return fluences

def compute_fluence_geom_dependence_task((i, j), (rm, rb), inversion, doping_agent, pulse_photon_count, (int_type, amp_type), (count_z, count_t), decay):
    output.show_status((i, j), params.extended_status_strides, False)
    
    medium_radius_orig = params.medium_radius
    beam_radius_orig = params.beam_radius
    params.medium_radius = rm
    params.beam_radius = rb
    
    ref_inversion = inversion
    initial_inversion = model.inversion.UniformInversion(ref_inversion)
    active_medium = model.medium.ActiveMedium(initial_inversion, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    
    ref_fluence = params.beam_class.ref_fluence(params.beam_radius, pulse_photon_count)
    beam_profile = params.beam_class(params.beam_radius, ref_fluence)
    rho, phi = beam_profile.rho_ref, beam_profile.phi_ref
    ref_pulse = core.create_pulse(active_medium, beam_profile, rho, phi)
    
    integrator = model.energy.PhotonCountIntegrator(int_type, active_medium, beam_profile)
    amp = amp_type(active_medium, count_z)
    
    upper = np.vectorize(initial_inversion.inversion)(rho, phi, amp.Z)
    lower = np.zeros(count_z)
    population = (upper, lower)
    
    pulse_fluences = np.empty(params.train_pulse_count)
    
    for pnum in range(params.train_pulse_count):
        density_out, population_out = amp.amplify(rho, phi, ref_pulse, count_t, initial_population=population)
        
        upper = np.copy(population_out[0])
        lower = population_out[1] * decay
        population = (upper, lower)
        
        fluence_out = integrator.integral(amp.T, density_out) * active_medium.light_speed
        pulse_fluences[pnum] = fluence_out
    
    fluence_out = pulse_fluences[::-1].sum()
    fluence_out = model.energy.energy(params.lasing_wavelen, fluence_out)
    
    params.medium_radius = medium_radius_orig
    params.beam_radius = beam_radius_orig
    
    return fluence_out

def compute_fluence_geom_dependence(task_pool, dirname, (int_types, amp_types), inversions, (int_type, amp_type), (count_z, count_t)):
    filename = lambda name: os.path.join(dirname, name)
    
    print output.div_line
    print "computing fluence dependence on geometry parameters"
    
    doping_agent = model.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    pulse_photon_count = model.energy.photon_count(params.lasing_wavelen, params.pulse_energy)
    
    active_medium = model.medium.ActiveMedium(None, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    ref_fluence = params.beam_class.ref_fluence(params.beam_radius, pulse_photon_count)
    beam_profile = params.beam_class(params.beam_radius, ref_fluence)
    ref_pulse = core.create_pulse(active_medium, beam_profile, beam_profile.rho_ref, beam_profile.phi_ref)
    decay = core.get_decay(active_medium, ref_pulse)
    
    min_medium_radius = params.ext_opt_geom_mediumradius[0]
    min_beam_radius = params.ext_opt_geom_beamradius[0]
    
    count_rm = params.ext_opt_geom_resolution[0]
    count_rb = params.ext_opt_geom_resolution[1]
    
    Rm = np.linspace(min_medium_radius, params.ext_opt_geom_mediumradius[1], count_rm)
    Rb = np.linspace(min_beam_radius, params.ext_opt_geom_beamradius[1], count_rb)
    
    inversions2d = np.meshgrid(inversions, Rb)[0].T
    fluences = task_pool.parallel_task(compute_fluence_geom_dependence_task, (Rm, Rb), (inversions2d,), (doping_agent, pulse_photon_count, (int_type, amp_type), (count_z, count_t), decay))
    
    output.show_status((count_rm, count_rb), params.extended_status_strides, True)
    
    if params.graphs:
        print "generating output"
        dirname = os.path.join(dirname, output.opt_geom_rel_path)
        dirname = output.init_dir(dirname)
        fluence_max = params.ext_opt_fluence_max
        zlim = None #(0.0, fluence_max)
        contours = [fluence_max]
        plot.plot_color(filename("fluence_out_max"), "Maximum Output Fluence", (Rm, None, None, output.medium_radius_label), (Rb, None, None, output.beam_radius_label), (fluences.T, zlim, output.fluence_abs_label_energy), params.num_auto_contours, contours)
    
    return fluences

def output_pump_constraints(dirname, inversion_rdiffs, max_fluences):
    filename = lambda name: os.path.join(dirname, name)
    
    print output.div_line
    print "computing pumping parameters domain from constraints"
    
    doping_agent = model.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = model.medium.ActiveMedium(None, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    
    count_tau = params.ext_opt_pump_resolution[0]
    count_pwr = params.ext_opt_pump_resolution[1]
    
    Tau = np.linspace(params.ext_opt_pump_duration[0], params.ext_opt_pump_duration[1], count_tau)
    Pwr = np.linspace(params.ext_opt_pump_power[0], params.ext_opt_pump_power[1], count_pwr)
    
    inversion_rdiff_max = params.ext_opt_inversion_rdiff_max
    fluence_max = params.ext_opt_fluence_max
    contours = [
        (inversion_rdiffs.T, inversion_rdiff_max, "depopulation"),
        (max_fluences.T, fluence_max, "damage"),
    ]
    contour_comps = [1.0, 1.0]
    
    if params.graphs:
        print "generating output"
        dirname = os.path.join(dirname, output.opt_pump_rel_path)
        graph_types = [
            (dirname, Pwr,output. pump_power_label),
            (os.path.join(dirname, output.alt_plot_rel_path), Pwr * params.pump_efficiency / active_medium.volume, output.eff_power_density_label),
        ]
        for dirname, Y, ylabel in graph_types:
            dirname = output.init_dir(dirname)
            plot.plot_contour(filename("constraints"), "Pumping Parameters Domain Constraints", (Tau, None, None, output.pump_duration_label), (Y, None, None, ylabel), contours)
    
    return (contours, None, None), (contour_comps, None, None)

def output_geom_constraints(dirname, inversion_rdiffs, max_fluences):
    filename = lambda name: os.path.join(dirname, name)
    
    print output.div_line
    print "computing geometry parameters domain from constraints"
    
    count_rm = params.ext_opt_geom_resolution[0]
    count_rb = params.ext_opt_geom_resolution[1]
    
    Rm = np.linspace(params.ext_opt_geom_mediumradius[0], params.ext_opt_geom_mediumradius[1], count_rm)
    Rb = np.linspace(params.ext_opt_geom_beamradius[0], params.ext_opt_geom_beamradius[1], count_rb)
    
    inversion_rdiff_max = params.ext_opt_inversion_rdiff_max
    fluence_max = params.ext_opt_fluence_max
    contours = [
        (max_fluences.T, fluence_max, "damage"),
    ]
    contour_comps = [1.0]
    
    rm_depop = np.interp(inversion_rdiff_max, inversion_rdiffs[::-1], Rm[::-1])
    xvals = [(rm_depop, "depopulation")]
    xval_comps = [-1.0]
    
    if params.graphs:
        print "generating output"
        dirname = os.path.join(dirname, output.opt_geom_rel_path)
        dirname = output.init_dir(dirname)
        plot.plot_contour(filename("constraints"), "Geometry Parameter Domain Constraints", (Rm, None, None, output.medium_radius_label), (Rb, None, None, output.beam_radius_label), contours, xvals=xvals)
    
    return (contours, xvals, None), (contour_comps, xval_comps, None)

def optimize_output(domain, outputs, limits, comparisons, price):
    params1, params2 = domain
    contours, xvals, yvals = limits
    contour_comps, xval_comps, yval_comps = comparisons
    optimum = None
    for i, p1 in enumerate(params1):
        for j, p2 in enumerate(params2):
            safe = True
            if contours:
                for limit, comp in zip(contours, contour_comps):
                    Z, zmax, _ = limit
                    Z = Z.T
                    if Z[i, j] * comp > zmax * comp:
                        safe = False
                        break
            if xvals:
                for limit, comp in zip(xvals, xval_comps):
                    val, _ = limit
                    if p1 * comp > val * comp:
                        safe = False
                        break
            if yvals:
                for limit, comp in zip(yvals, yval_comps):
                    val, _ = limit
                    if p2 *comp > val * comp:
                        safe = False
                        break
            if not safe:
                continue
            if not optimum:
                optimum = (i, j)
            else:
                output = outputs[i, j]
                optimum_output = outputs[optimum]
                if output > optimum_output:
                    optimum = (i, j)
                elif output == optimum_output:
                    current_price = price(p1, p2)
                    optimum_price = price(params1[optimum[0]], params2[optimum[1]])
                    if current_price < optimum_price:
                        optimum = (i, j)
    return optimum

def compute_energy_pump_dependence_task((i, j), (tau, pwr), inversion, num_types, counts):
    output.show_status((i, j), params.extended_status_strides, False)
    _, _, output_energy, rel_gain_reduction = core.amplify_train(None, num_types, counts, inversion, quiet=True)
    return output_energy, rel_gain_reduction

def compute_energy_pump_dependence(task_pool, dirname, (int_types, amp_types), inversions, constraints, num_types, counts):
    filename = lambda name: os.path.join(dirname, name)
    
    print output.div_line
    print "computing energy dependence on pumping parameters"
    
    doping_agent = model.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = model.medium.ActiveMedium(None, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    
    pulse_photon_count = model.energy.photon_count(params.lasing_wavelen, params.pulse_energy)
    ref_fluence = params.beam_class.ref_fluence(params.beam_radius, pulse_photon_count)
    beam_profile = params.beam_class(params.beam_radius, ref_fluence)
    
    input_photon_count = beam_profile.fluence_integral(active_medium.radius)
    input_energy = model.energy.energy(params.lasing_wavelen, input_photon_count)
    input_energy *= params.train_pulse_count
    
    count_tau = params.ext_opt_pump_resolution[0]
    count_pwr = params.ext_opt_pump_resolution[1]
    
    Tau = np.linspace(params.ext_opt_pump_duration[0], params.ext_opt_pump_duration[1], count_tau)
    Pwr = np.linspace(params.ext_opt_pump_power[0], params.ext_opt_pump_power[1], count_pwr)
    
    output_energies, rel_gain_reductions = task_pool.parallel_task(compute_energy_pump_dependence_task, (Tau, Pwr), (inversions,), (num_types, counts))
    
    output.show_status((count_tau, count_pwr), params.extended_status_strides, True)
    
    pump_energies = np.prod(np.array(np.meshgrid(Tau, Pwr)), axis=0).T
    stored_energies = model.energy.energy(params.lasing_wavelen, inversions * active_medium.volume)
    energy_gains = output_energies / input_energy
    added_energies = output_energies - input_energy
    extraction_effs = added_energies / stored_energies
    total_effs = added_energies / pump_energies
    
    limits, comparisons = constraints
    price = lambda tau, pwr: tau * pwr
    
    unitconv.print_result("input energy [{}]: {}", ("mJ",), (input_energy,))
    
    optimum = optimize_output((Tau, Pwr), output_energies, limits, comparisons, price)
    output_energy_optimum_params = (Tau[optimum[0]], Pwr[optimum[1]]) if optimum else (None, None)
    unitconv.print_result("max output energy [{}]: {}", ("mJ",), (output_energies[optimum] if optimum else None,))
    unitconv.print_result("optimum pumping parameters (duration [{}], power [{}]): ({}, {})", ("us", "W"), output_energy_optimum_params)
    
    optimum = optimize_output((Tau, Pwr), energy_gains, limits, comparisons, price)
    energy_gain_optimum_params = (Tau[optimum[0]], Pwr[optimum[1]]) if optimum else (None, None)
    unitconv.print_result("max energy gain: {}", (), (energy_gains[optimum] if optimum else None,))
    unitconv.print_result("optimum pumping parameters (duration [{}], power [{}]): ({}, {})", ("us", "W"), energy_gain_optimum_params)
    
    optimum = optimize_output((Tau, Pwr), extraction_effs, limits, comparisons, price)
    extraction_eff_optimum_params = (Tau[optimum[0]], Pwr[optimum[1]]) if optimum else (None, None)
    unitconv.print_result("max extraction efficiency [{}]: {}", ("%",), (extraction_effs[optimum] if optimum else None,))
    unitconv.print_result("optimum pumping parameters (duration [{}], power [{}]): ({}, {})", ("us", "W"), extraction_eff_optimum_params)
    
    optimum = optimize_output((Tau, Pwr), total_effs, limits, comparisons, price)
    total_eff_optimum_params = (Tau[optimum[0]], Pwr[optimum[1]]) if optimum else (None, None)
    unitconv.print_result("max opt-opt efficiency [{}]: {}", ("%",), (total_effs[optimum] if optimum else None,))
    unitconv.print_result("optimum pumping parameters (duration [{}], power [{}]): ({}, {})", ("us", "W"), total_eff_optimum_params)
    
    optimum = optimize_output((Tau, Pwr), -rel_gain_reductions, limits, comparisons, price)
    rel_gain_reduction_optimum_params = (Tau[optimum[0]], Pwr[optimum[1]]) if optimum else (None, None)
    unitconv.print_result("min rel gain reduction [{}]: {}", ("%",), (rel_gain_reductions[optimum] if optimum else None,))
    unitconv.print_result("optimum pumping parameters (duration [{}], power [{}]): ({}, {})", ("us", "W"), rel_gain_reduction_optimum_params)
    
    if params.graphs:
        print "generating output"
        extra_contours, xvals, yvals = limits
        dirname = os.path.join(dirname, output.opt_pump_rel_path)
        graph_types = [
            (dirname, Pwr, output.pump_power_label),
            (os.path.join(dirname, output.alt_plot_rel_path), Pwr * params.pump_efficiency / active_medium.volume, output.eff_power_density_label),
        ]
        for dirname, Y, ylabel in graph_types:
            dirname = output.init_dir(dirname)
            plot.plot_color(filename("energy_out"), "Output Energy", (Tau, None, None, output.pump_duration_label), (Y, None, None, ylabel), (output_energies.T, None, output.energy_abs_pulse_label), params.num_auto_contours, extra_contours=extra_contours, xvals=xvals, yvals=yvals)
            plot.plot_color(filename("energy_gain"), "Energy Gain", (Tau, None, None, output.pump_duration_label), (Y, None, None, ylabel), (energy_gains.T, None, output.energy_rel_label), params.num_auto_contours, extra_contours=extra_contours, xvals=xvals, yvals=yvals)
            plot.plot_color(filename("efficiency_extr"), "Extraction Efficiency", (Tau, None, None, output.pump_duration_label), (Y, None, None, ylabel), (extraction_effs.T, None, output.extraction_eff_label), params.num_auto_contours, extra_contours=extra_contours, xvals=xvals, yvals=yvals)
            plot.plot_color(filename("efficiency_opt2"), "Optical to Optical Efficiency", (Tau, None, None, output.pump_duration_label), (Y, None, None, ylabel), (total_effs.T, None, output.total_eff_label), params.num_auto_contours, extra_contours=extra_contours, xvals=xvals, yvals=yvals)
            if params.train_pulse_count > 1:
                plot.plot_color(filename("gain_reduction"), "Gain Reduction", (Tau, None, None, output.pump_duration_label), (Y, None, None, ylabel), (rel_gain_reductions.T, None, output.rel_gain_reduction_label), params.num_auto_contours, extra_contours=extra_contours, xvals=xvals, yvals=yvals)

def compute_energy_geom_dependence_task((i, j), (rm, rb), inversion, doping_agent, pulse_photon_count, num_types, counts):
    output.show_status((i, j), params.extended_status_strides, False)
    medium_radius_orig = params.medium_radius
    beam_radius_orig = params.beam_radius
    params.medium_radius = rm
    params.beam_radius = rb
    ref_inversion = inversion
    active_medium = model.medium.ActiveMedium(None, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    stored_energy = model.energy.energy(params.lasing_wavelen, ref_inversion * active_medium.volume)
    ref_fluence = params.beam_class.ref_fluence(params.beam_radius, pulse_photon_count)
    beam_profile = params.beam_class(params.beam_radius, ref_fluence)
    input_photon_count = beam_profile.fluence_integral(active_medium.radius)
    input_energy = model.energy.energy(params.lasing_wavelen, input_photon_count)
    input_energy *= params.train_pulse_count
    _, _, output_energy, rel_gain_reduction = core.amplify_train(None, num_types, counts, ref_inversion, quiet=True)
    params.medium_radius = medium_radius_orig
    params.beam_radius = beam_radius_orig
    return stored_energy, input_energy, output_energy, rel_gain_reduction

def compute_energy_geom_dependence(task_pool, dirname, (int_types, amp_types), inversions, constraints, num_types, counts):
    filename = lambda name: os.path.join(dirname, name)
    
    print output.div_line
    print "computing energy dependence on geometry parameters"
    
    doping_agent = model.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, params.dopant_lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    pulse_photon_count = model.energy.photon_count(params.lasing_wavelen, params.pulse_energy)
    
    min_medium_radius = params.ext_opt_geom_mediumradius[0]
    min_beam_radius = params.ext_opt_geom_beamradius[0]
    
    count_rm = params.ext_opt_geom_resolution[0]
    count_rb = params.ext_opt_geom_resolution[1]
    
    Rm = np.linspace(min_medium_radius, params.ext_opt_geom_mediumradius[1], count_rm)
    Rb = np.linspace(min_beam_radius, params.ext_opt_geom_beamradius[1], count_rb)
    
    inversions2d = np.meshgrid(inversions, Rb)[0].T
    stored_energies, input_energies, output_energies, rel_gain_reductions = task_pool.parallel_task(compute_energy_geom_dependence_task, (Rm, Rb), (inversions2d,), (doping_agent, pulse_photon_count, num_types, counts))
    
    output.show_status((count_rm, count_rb), params.extended_status_strides, True)
    
    pump_energy = params.pump_duration * params.pump_power
    energy_gains = output_energies / input_energies
    added_energies = output_energies - input_energies
    extraction_effs = added_energies / stored_energies
    total_effs = added_energies / pump_energy
    
    limits, comparisons = constraints
    price = lambda rm, rb: 1.0 / (rm * rb)
    
    optimum = optimize_output((Rm, Rb), output_energies, limits, comparisons, price)
    output_energy_optimum_params = (Rm[optimum[0]], Rb[optimum[1]]) if optimum else (None, None)
    unitconv.print_result("max output energy [{}]: {}", ("mJ",), (output_energies[optimum] if optimum else None,))
    unitconv.print_result("optimum geometry parameters (medium diameter [{}], beam diameter [{}]): ({}, {})", ("mm", "mm"), output_energy_optimum_params)
    
    optimum = optimize_output((Rm, Rb), energy_gains, limits, comparisons, price)
    energy_gain_optimum_params = (Rm[optimum[0]], Rb[optimum[1]]) if optimum else (None, None)
    unitconv.print_result("max energy gain: {}", (), (energy_gains[optimum] if optimum else None,))
    unitconv.print_result("optimum geometry parameters (medium diameter [{}], beam diameter [{}]): ({}, {})", ("mm", "mm"), energy_gain_optimum_params)
    
    optimum = optimize_output((Rm, Rb), extraction_effs, limits, comparisons, price)
    extraction_eff_optimum_params = (Rm[optimum[0]], Rb[optimum[1]]) if optimum else (None, None)
    unitconv.print_result("max extraction efficiency [{}]: {}", ("%",), (extraction_effs[optimum] if optimum else None,))
    unitconv.print_result("optimum geometry parameters (medium diameter [{}], beam diameter [{}]): ({}, {})", ("mm", "mm"), extraction_eff_optimum_params)
    
    optimum = optimize_output((Rm, Rb), total_effs, limits, comparisons, price)
    total_eff_optimum_params = (Rm[optimum[0]], Rb[optimum[1]]) if optimum else (None, None)
    unitconv.print_result("max opt-opt efficiency [{}]: {}", ("%",), (total_effs[optimum] if optimum else None,))
    unitconv.print_result("optimum geometry parameters (medium diameter [{}], beam diameter [{}]): ({}, {})", ("mm", "mm"), total_eff_optimum_params)
    
    optimum = optimize_output((Rm, Rb), -rel_gain_reductions, limits, comparisons, price)
    rel_gain_reduction_optimum_params = (Rm[optimum[0]], Rb[optimum[1]]) if optimum else (None, None)
    unitconv.print_result("min rel gain reduction [{}]: {}", ("%",), (rel_gain_reductions[optimum] if optimum else None,))
    unitconv.print_result("optimum geometry parameters (medium diameter [{}], beam diameter [{}]): ({}, {})", ("mm", "mm"), rel_gain_reduction_optimum_params)
    
    if params.graphs:
        print "generating output"
        extra_contours, xvals, yvals = limits
        dirname = os.path.join(dirname, output.opt_geom_rel_path)
        dirname = output.init_dir(dirname)
        plot.plot_color(filename("energy_in"), "Input Energy", (Rm, None, None, output.medium_radius_label), (Rb, None, None, output.beam_radius_label), (input_energies.T, None, output.energy_abs_pulse_label), params.num_auto_contours, extra_contours=extra_contours, xvals=xvals, yvals=yvals)
        plot.plot_color(filename("energy_out"), "Output Energy", (Rm, None, None, output.medium_radius_label), (Rb, None, None, output.beam_radius_label), (output_energies.T, None, output.energy_abs_pulse_label), params.num_auto_contours, extra_contours=extra_contours, xvals=xvals, yvals=yvals)
        plot.plot_color(filename("energy_gain"), "Energy Gain", (Rm, None, None, output.medium_radius_label), (Rb, None, None, output.beam_radius_label), (energy_gains.T, None, output.energy_rel_label), params.num_auto_contours, extra_contours=extra_contours, xvals=xvals, yvals=yvals)
        plot.plot_color(filename("efficiency_extr"), "Extraction Efficiency", (Rm, None, None, output.medium_radius_label), (Rb, None, None, output.beam_radius_label), (extraction_effs.T, None, output.extraction_eff_label), params.num_auto_contours, extra_contours=extra_contours, xvals=xvals, yvals=yvals)
        plot.plot_color(filename("efficiency_opt2"), "Optical to Optical Efficiency", (Rm, None, None, output.medium_radius_label), (Rb, None, None, output.beam_radius_label), (total_effs.T, None, output.total_eff_label), params.num_auto_contours, extra_contours=extra_contours, xvals=xvals, yvals=yvals)
        if params.train_pulse_count > 1:
            plot.plot_color(filename("gain_reduction"), "Gain Reduction", (Rm, None, None, output.medium_radius_label), (Rb, None, None, output.beam_radius_label), (rel_gain_reductions.T, None, output.rel_gain_reduction_label), params.num_auto_contours, extra_contours=extra_contours, xvals=xvals, yvals=yvals)

def compare_lower_lifetimes(dirname, ref_inversion, (int_types, amp_types), numerics):
    filename = lambda name: os.path.join(dirname, name)
    
    print output.div_line
    print "comparing lower state lifetimes"
    
    if numerics is None:
        numerics = core.setup_methods(dirname, (int_types, amp_types), ref_inversion, quiet=True)
    num_types, counts = numerics
    
    lower_lifetime = params.dopant_lower_lifetime
    
    initial_inversion = model.inversion.UniformInversion(ref_inversion)
    
    doping_agent = model.dopant.DopingAgent(params.dopant_xsection, params.dopant_upper_lifetime, lower_lifetime, params.dopant_branching_ratio, params.dopant_concentration)
    active_medium = model.medium.ActiveMedium(initial_inversion, doping_agent, params.medium_radius, params.medium_length, params.medium_refr_idx)
    active_medium_3 = copy.deepcopy(active_medium)
    active_medium_3.doping_agent.lower_lifetime = float("inf")
    active_medium_4 = copy.deepcopy(active_medium)
    active_medium_4.doping_agent.lower_lifetime = 0.0
    
    pulse_photon_count = model.energy.photon_count(params.lasing_wavelen, params.pulse_energy)
    ref_fluence = params.beam_class.ref_fluence(params.beam_radius, pulse_photon_count)
    beam_profile = params.beam_class(params.beam_radius, ref_fluence)
    
    ref_pulse = core.create_pulse(active_medium, beam_profile, beam_profile.rho_ref, beam_profile.phi_ref)
    
    (int_type, amp_type), (_, _, count_z, count_t) = num_types, counts
    
    integrator = model.energy.PhotonCountIntegrator(int_type, active_medium, beam_profile)
    
    amp = amp_type(active_medium, count_z)
    amp_3 = model.amplifier.ExactOutputAmplifier(active_medium_3, count_z)
    amp_4 = model.amplifier.ExactOutputAmplifier(active_medium_4, count_z)
    amplify_args = (beam_profile.rho_ref, beam_profile.phi_ref, ref_pulse, count_t)
    
    print "zero"
    density_out_4, _ = amp_4.amplify(*amplify_args)
    fluence_out_4 = integrator.integral(amp_4.T, density_out_4) * active_medium_4.light_speed
    fluence_gain_4 = fluence_out_4 / ref_fluence
    unitconv.print_result("fluence gain: {}", (), (fluence_gain_4,))
    
    lsl_output_label = "finite"
    if lower_lifetime == 0.0:
        lsl_output_label = "zero"
    elif math.isinf(lower_lifetime):
        lsl_output_label = "infinite"
    print lsl_output_label
    density_out, _ = amp.amplify(*amplify_args)
    fluence_out = integrator.integral(amp.T, density_out) * active_medium.light_speed
    fluence_gain = fluence_out / ref_fluence
    unitconv.print_result("fluence gain: {}", (), (fluence_gain,))
    
    print "infinite"
    density_out_3, _ = amp_3.amplify(*amplify_args)
    fluence_out_3 = integrator.integral(amp_3.T, density_out_3) * active_medium_3.light_speed
    fluence_gain_3 = fluence_out_3 / ref_fluence
    unitconv.print_result("fluence gain: {}", (), (fluence_gain_3,))
    
    if params.graphs:
        print "generating output"
        dirname = os.path.join(dirname, output.ref_pulse_rel_path)
        dirname = output.init_dir(dirname)
        T = amp.T
        if params.output_rel_time:
            T = T / params.pulse_duration
        out_t_label = output.norm_t_label if params.output_rel_time else output.t_amp_label
        tlim = (T[0], T[-1])
        Ts = (T,) * 3
        ref_density = ref_pulse.ref_density
        densities = (density_out_4 / ref_density, density_out / ref_density, density_out_3 / ref_density)
        lifetime_scale = unitconv.units[output.lower_lifetime_unit]
        lsl_graph_label_fmt = (r"%g \, %s" % (lower_lifetime/lifetime_scale, output.lower_lifetime_unit))
        if lower_lifetime == 0.0:
            lsl_graph_label_fmt = "0"
        elif math.isinf(lower_lifetime):
            lsl_graph_label_fmt = "\\infty"
        labels = (output.lower_lifetime_legend % "0", output.lower_lifetime_legend % lsl_graph_label_fmt, output.lower_lifetime_legend % "\\infty")
        plot.plot_data(filename("lsl_effects"), "Effects of Lower State Lifetime", (Ts, None, tlim, out_t_label), (densities, None, None, output.density_rel_label), labels)

def extended_mode(task_pool, dirname, ref_inversion, (int_types, amp_types), numerics):
    if params.amplification:
        compare_lower_lifetimes(dirname, ref_inversion, (int_types, amp_types), numerics)
    
    if not params.initial_inversion:
        compare_depop_models(dirname)
        inversions_pump, inversion_rdiffs_pump = compute_inversion_pump_dependence(task_pool, dirname)
        inversions_geom, inversion_rdiffs_geom = compute_inversion_geom_dependence(task_pool, dirname)
        if params.amplification:
            print output.div_line
            print "determining extended mode method combinations"
            
            print "pumping"
            max_inversion_pump = inversions_pump[-1, -1]
            num_types_pump, counts_pump = core.setup_methods(dirname, (int_types, amp_types), max_inversion_pump, quiet=True)
            _, _, count_z_pump, count_t_pump = counts_pump
            
            print "geometry"
            min_medium_radius = params.ext_opt_geom_mediumradius[0]
            min_beam_radius = params.ext_opt_geom_beamradius[0]
            medium_radius_orig = params.medium_radius
            beam_radius_orig = params.beam_radius
            params.medium_radius = min_medium_radius
            params.beam_radius = min_beam_radius
            max_inversion_geom = inversions_geom[0]
            num_types_geom, counts_geom = core.setup_methods(dirname, (int_types, amp_types), max_inversion_geom, quiet=True)
            _, _, count_z_geom, count_t_geom = counts_geom
            params.medium_radius = medium_radius_orig
            params.beam_radius = beam_radius_orig
            
            max_fluences_pump = compute_fluence_pump_dependence(task_pool, dirname, (int_types, amp_types), inversions_pump, num_types_pump, (count_z_pump, count_t_pump))
            max_fluences_geom = compute_fluence_geom_dependence(task_pool, dirname, (int_types, amp_types), inversions_geom, num_types_geom, (count_z_geom, count_t_geom))
            pump_constraints = output_pump_constraints(dirname, inversion_rdiffs_pump, max_fluences_pump)
            geom_constraints = output_geom_constraints(dirname, inversion_rdiffs_geom, max_fluences_geom)
            compute_energy_pump_dependence(task_pool, dirname, (int_types, amp_types), inversions_pump, pump_constraints, num_types_pump, counts_pump)
            compute_energy_geom_dependence(task_pool, dirname, (int_types, amp_types), inversions_geom, geom_constraints, num_types_geom, counts_geom)
