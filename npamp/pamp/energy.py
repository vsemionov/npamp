
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

import const
import amplifier
import util


def photon_energy(wavelen_vacuum):
    freq = const.light_speed_vacuum / wavelen_vacuum
    energy = const.planck_constant * freq
    return energy

def photon_count(wavelen_vacuum, energy):
    single_photon_energy = photon_energy(wavelen_vacuum)
    photon_count = energy / single_photon_energy
    return photon_count

def energy(wavelen_vacuum, photon_count):
    single_photon_energy = photon_energy(wavelen_vacuum)
    energy = photon_count * single_photon_energy
    return energy


class PhotonCountIntegrator(object):
    
    def __init__(self, integrator_type, active_medium, beam_profile):
        self.integrator_type = integrator_type
        self.active_medium = active_medium
        self.beam_profile = beam_profile
        
        self.integrator = integrator_type()
    
    def min_steps(self, single_pulses, energy_rtol, fluence_rtol):
        def converge_steps(func, a, b):
            max_divs = 12
            divs = 1
            last_res = None
            while True:
                steps = util.steps(divs)
                X = np.linspace(a, b, steps)
                Y = np.vectorize(func)(X)
                res = self.integral(X, Y)
                diff = None
                if last_res is not None:
                    diff = abs(res - last_res)
                    if diff <= rtol * abs(res):
                        break
                last_res = res
                divs += 1
                if divs > max_divs:
                    if diff is not None:
                        warnings.warn("max divs (%d) exceeded; latest difference: %f" % (max_divs, diff), stacklevel=3)
                        break
                    else:
                        raise ValueError("max divs (%d) exceeded" % max_divs)
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
            num0 = self.integrate_base(Rho, Phi, fluence_min)
            exact0 = beam_profile.fluence_integral(active_medium.radius)
            rel_error_0 = abs((num0 - exact0) / exact0)
            num1 = self.integrate_base(Rho, Phi, fluence_max)
            exact1 = beam_profile.fluence_integral(active_medium.radius) * gain
            rel_error_1 = abs((num1 - exact1) / exact1)
            num3 = self.integrate_base(Rho, Phi, fluence_max_3)
            exact3 = exact0 + inversion_fluence / 2.0 * active_medium.aperture
            rel_error_3 = abs((num3 - exact3) / exact3)
            num4 = self.integrate_base(Rho, Phi, fluence_max_4)
            exact4 = exact0 + inversion_fluence * active_medium.aperture
            rel_error_4 = abs((num4 - exact4) / exact4)
            result = (num0, num1, num3, num4)
            rel_error = max(rel_error_0, rel_error_1, rel_error_3, rel_error_4)
            return result, rel_error
        
        beam_profile = self.beam_profile
        active_medium = self.active_medium
        L = active_medium.length
        gain = math.exp(active_medium.doping_agent.xsection * active_medium.initial_inversion.ref_inversion * L)
        
        rtol = energy_rtol
        
        beam_rho, beam_phi = self.beam_profile.xcoords
        
        if not beam_rho and not beam_phi:
            steps_rho = 1
            steps_phi = 1
        else:
            compute_rdiff = lambda last_res, res: max([abs((res[i] - last_res[i]) / res[i]) for i in range(len(res))])
            steps_rho, steps_phi = util.min_steps((1, 1), (beam_rho, beam_phi), rtol, fluence_integrals, compute_rdiff)
        
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
    
    def integral(self, X, Y):
        assert X.ndim == Y.ndim == 1
        assert X.shape == Y.shape
        assert len(X) > 2 # scipy.integrate.romb() seems to require 3 or more samples
        
        divs_f = math.log(len(X) - 1, 2.0)
        divs = int(divs_f)
        assert divs == divs_f
        
        dx = (X[-1] - X[0]) / (len(X) - 1)
        I = self.integrator(Y, dx)
        return I
    
    def integrate_base(self, Rho, Phi, fluence):
        assert Rho.ndim == Phi.ndim == 1
        assert fluence.shape == (Rho.shape + Phi.shape)
        
        integrator = lambda Y, X, xmax: self.integral(X, Y) if len(X) > 1 else xmax*Y[0]
        
        phi_integrals = np.apply_along_axis(integrator, 1, fluence, Phi, 2.0*math.pi)
        radius = self.active_medium.radius
        phi_integrals *= Rho if len(Rho) > 1 else radius/2.0
        rho_intergral = integrator(phi_integrals, Rho, radius)
        
        return rho_intergral
