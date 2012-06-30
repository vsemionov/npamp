
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

import numpy as np

import params
import plot


div_line = "=" * 32

output_dir = None
models_rel_path = "pumping"
ref_pulse_rel_path = "ref_pulse"
optimization_rel_path = "optimization"
pumpdep_rel_path = os.path.join(optimization_rel_path, "pump")
geomdep_rel_path = os.path.join(optimization_rel_path, "geom")
alt_plot_rel_path = "alt"

x_label = "x [mm]"
y_label = "y [mm]"
z_label = "z [mm]"
rho_label = "r [mm]"
t_amp_label = "t [ns]"
i_label = "pulse num."
norm_t_label = "t/T"
density_rel_label = "rel. photon density"
density_norm_rel_label = "norm. photon density"
upper_rel_label = "rel. upper state population"
lower_rel_label = "rel. lower state population"
inversion_rel_label = "rel. population inversion"
inversion_abs_label = "population inversion [cm^-3]"
t_pump_label = "t [us]"
pump_duration_label = "pump duration [us]"
pump_power_label = "pump power [W]"
eff_power_density_label = "absorbed power density [W/cm^3]"
rate_label = "depopulation rate [cm^-3 s^-1]"
rate_rel_label = "Rase/n [s^-1]"
gain_label = "small-signal gain"
error_label = "rel. error"
inversion_rdiff_label = "inversion rel. difference [%]"
gain_rdiff_label = "gain rel. difference [%]"
energy_rel_label = "energy gain"
energy_abs_pump_label = "optical pump energy [J]"
energy_abs_stored_label = "stored energy [J]"
energy_abs_pulse_label = "output energy [mJ]"
rel_gain_reduction_label = "rel. gain reduction [%]"
fluence_rel_label = "rel. fluence"
fluence_norm_rel_label = "norm. fluence"
fluence_abs_label_energy = "max. output fluence [J/cm^2]"
medium_radius_label = "medium diameter [mm]"
beam_radius_label = "beam diameter [mm]"
extraction_eff_label = "extraction efficiency [%]"
total_eff_label = "optical to optical efficiency [%]"

lower_lifetime_legend = r"$\tau_1 \, = \, %s$"
lower_lifetime_unit = "ps"


def init_dir(name):
    dirname = os.path.join(output_dir, name)
    if not os.path.isdir(dirname):
        os.makedirs(dirname)
    return dirname


def plot_inversion(dirname, inv):
    filename = lambda name: os.path.join(dirname, name)
    
    T = inv.T
    inversion = inv.inversion
    tlim = (T[0], T[-1])
    plot.plot_data(filename("inversion_evo"), "Population Inversion Evolution", (T, None, tlim, t_pump_label), (inversion, None, None, inversion_abs_label))

def plot_output(dirname, beam_profile, input_pulse, fwhm, amp, fluences, exact_density_out=None, exact_population_out=None):
    filename = lambda name: os.path.join(dirname, name)
    
    density = amp.density
    population = amp.population
    upper = population[0]
    lower = population[1]
    inversion = upper - lower
    
    Z = amp.Z
    T = amp.T
    
    if params.output_rel_time:
        T = T / fwhm
    
    TZ, ZT = np.meshgrid(T, Z)
    
    zlim = (Z[0], Z[-1])
    tlim = (T[0], T[-1])
    
    ref_density = input_pulse.ref_density
    ref_inversion = amp.active_medium.initial_inversion.ref_inversion
    
    out_t_label = norm_t_label if params.output_rel_time else t_amp_label
    
    stride_z = max(amp.count_z // params.out_z_steps_divisor, 1)
    stride_t = max(amp.count_t // params.out_t_steps_divisor, 1)
    
    plot.plot_data(filename("density_in"), "Input Photon Density", (T, None, tlim, out_t_label), (density[0]/ref_density, None, None, density_rel_label))
    plot.plot_data(filename("density_out"), "Output Photon Density", (T, None, tlim, out_t_label), (density[-1]/ref_density, None, None, density_rel_label))
    plot.plot_data(filename("densities"), "Input and Output Photon Densities", ((T, ) * 2, None, tlim, out_t_label), ((density[0]/ref_density, density[-1]/ref_density), None, None, density_rel_label), ("input pulse", "output pulse"))
    plot.plot_data(filename("densities_norm"), "Normalized Input and Output Photon Densities", ((T, ) * 2, None, tlim, out_t_label), ((density[0]/ref_density, density[-1]/np.amax(density[-1])), None, None, density_norm_rel_label), ("input pulse", "output pulse"))
    
    plot.plot_data(filename("upper_init"), "Initial Upper State Population", (Z, None, zlim, z_label), (upper.T[0]/ref_inversion, None, None, upper_rel_label))
    plot.plot_data(filename("upper_final"), "Final Upper State Population", (Z, None, zlim, z_label), (upper.T[-1]/ref_inversion, None, None, upper_rel_label))
    plot.plot_data(filename("lower_init"), "Initial Lower State Population", (Z, None, zlim, z_label), (lower.T[0]/ref_inversion, None, None, lower_rel_label))
    plot.plot_data(filename("lower_final"), "Final Lower State Population", (Z, None, zlim, z_label), (lower.T[-1]/ref_inversion, None, None, lower_rel_label))
    plot.plot_data(filename("inversion_init"), "Initial Population Inversion", (Z, None, zlim, z_label), (inversion.T[0]/ref_inversion, None, None, inversion_rel_label))
    plot.plot_data(filename("inversion_final"), "Final Population Inversion", (Z, None, zlim, z_label), (inversion.T[-1]/ref_inversion, None, None, inversion_rel_label))
    
    plot.plot_projection(filename("density_evo"), "Photon Density Evolution", (ZT, None, z_label), (TZ, None, out_t_label), (density/ref_density, None, density_rel_label), (30, -30), (stride_z, stride_t))
    plot.plot_projection(filename("upper_evo"), "Upper State Population Evolution", (ZT, None, z_label), (TZ, None, out_t_label), (upper/ref_inversion, None, upper_rel_label), (30, 30), (stride_z, stride_t))
    plot.plot_projection(filename("lower_evo"), "Lower State Population Evolution", (ZT, None, z_label), (TZ, None, out_t_label), (lower/ref_inversion, None, lower_rel_label), (30, 30), (stride_z, stride_t))
    plot.plot_projection(filename("inversion_evo"), "Population Inversion Evolution", (ZT, None, z_label), (TZ, None, out_t_label), (inversion/ref_inversion, None, inversion_rel_label), (30, 30), (stride_z, stride_t))
    
    if exact_density_out is not None:
        plot.plot_error(filename("density_err"), "Photon Density Relative Error", (T, None, tlim, out_t_label), ((exact_density_out, density[-1]), None, None, error_label))
    if exact_population_out is not None:
        plot.plot_error(filename("upper_err"), "Upper State Population Relative Error", (Z, None, zlim, z_label), ((exact_population_out[0], upper.T[-1]), None, None, error_label))
        plot.plot_error(filename("lower_err"), "Lower State Population Relative Error", (Z, None, zlim, z_label), ((exact_population_out[1], lower.T[-1]), None, None, error_label))
        plot.plot_error(filename("inversion_err"), "Population Inversion Relative Error", (Z, None, zlim, z_label), ((exact_population_out[0] - exact_population_out[1], inversion.T[-1]), None, None, error_label))
    
    norm_fluences = fluences / beam_profile.ref_fluence
    plot.plot_data(filename("fluence"), "Pulse Fluence", (Z, None, zlim, z_label), (norm_fluences, None, None, fluence_rel_label))

def plot_train(dirname, beam_profile, active_medium, output_photon_counts):
    filename = lambda name: os.path.join(dirname, name)
    
    pulse_count = len(output_photon_counts)
    pulse_nums = np.arange(1, pulse_count + 1)
    nlim = (pulse_nums[0] - 1, pulse_nums[-1] + 1)
    extra_args = dict(style="o", vlines=True, grid="y") if pulse_count <= 32 else {}
    input_photon_count = beam_profile.fluence_integral(active_medium.radius)
    plot.plot_data(filename("train_energy_gain"), "Pulse Train Energy Gain", (pulse_nums, None, nlim, i_label), (output_photon_counts/input_photon_count, None, None, energy_rel_label), **extra_args)

def plot_beam(dirname, beam_profile, Rho, Phi, ref_output_fluence):
    filename = lambda name: os.path.join(dirname, name)
    
    ref_input_fluence = np.empty(ref_output_fluence.shape)
    for m, rho in enumerate(Rho):
        for n, phi in enumerate(Phi):
            ref_input_fluence[m, n] = beam_profile.fluence(rho, phi)
    FR, RF = np.meshgrid(Phi, Rho)
    XY, YX = RF * np.cos(FR), RF * np.sin(FR)
    norm_input_fluence = ref_input_fluence / beam_profile.ref_fluence
    norm_output_fluence = ref_output_fluence / beam_profile.ref_fluence
    scale = np.amax(norm_output_fluence)
    stride_rho = max(len(Rho) // params.out_rho_steps_divisor, 1)
    stride_phi = max(len(Phi) // params.out_phi_steps_divisor, 1)
    plot.plot_projection(filename("fluence_in"), "Input Pulse Fluence", (XY, None, x_label), (YX, None, y_label), (norm_input_fluence, None, fluence_rel_label), (30, -60), (stride_rho, stride_phi))
    plot.plot_projection(filename("fluence_out"), "Output Pulse Fluence", (XY, None, x_label), (YX, None, y_label), (norm_output_fluence, None, fluence_rel_label), (30, -60), (stride_rho, stride_phi))
    
    n_ref = -1
    for n, phi in enumerate(Phi):
        if n_ref < 0 or abs(phi - beam_profile.phi_ref) < abs(Phi[n_ref] - beam_profile.phi_ref):
            n_ref = n
    rholim = (Rho[0], Rho[-1])
    plot.plot_data(filename("fluences"), "Input and Output Pulse Fluences", ((Rho,)*2, None, rholim, rho_label), ((norm_input_fluence[:, n_ref], norm_output_fluence[:, n_ref]), None, None, fluence_rel_label), ("input beam", "output beam"))
    plot.plot_data(filename("fluences_norm"), "Normalized Input and Output Pulse Fluences", ((Rho,)*2, None, rholim, rho_label), ((norm_input_fluence[:, n_ref], norm_output_fluence[:, n_ref] / scale), None, None, fluence_norm_rel_label), ("input beam", "output beam"))


def show_status((i, j), (si, sj), done):
    def print_status():
        if j is not None:
            print "%d, %d" % (i, j)
        else:
            print i
    if si != 0:
        if done:
            print_status()
        else:
            if i % si == 0:
                if j is None:
                    print_status()
                else:
                    if sj == 0:
                        if j == 0:
                            print_status()
                    else:
                        if j % sj == 0:
                            print_status()
