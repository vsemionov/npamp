
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

import numpy as np

import native
import exc


def lower_state_decay(active_medium, pulse_train):
    separation = pulse_train.period - pulse_train.pulse.duration
    lower_lifetime = active_medium.doping_agent.lower_lifetime
    if lower_lifetime == 0.0:
        decay = 0.0
    elif math.isinf(lower_lifetime):
        decay = 1.0
    elif pulse_train.count == 1:
        decay = 0.0
    else:
        decay = math.exp(- separation / lower_lifetime)
    return decay


class PulseAmplifier(object):
    
    def __init__(self, active_medium, count_z):
        self.active_medium = active_medium
        
        self.count_z = count_z
        self.Z, self.dz = np.linspace(0.0, active_medium.length, self.count_z, retstep=True)
        
        self.rho = None
        self.phi = None
        
        self.input_pulse = None
        
        self.sizet = None
        self.count_t = None
        self.T = None
        self.dt = None
        
        self.density = None
        self.population = None
    
    @staticmethod
    def min_steps_z(active_medium):
        raise NotImplementedError()
    
    def min_steps_t(self, input_pulse):
        raise NotImplementedError()
    
    def _calc_population(self, l, k):
        raise NotImplementedError()
    
    def _calc_density(self, l, k):
        raise NotImplementedError()
    
    def _init_arrays(self):
        new_array = lambda old, shape: np.empty(shape) if old is None or old.shape != shape else old
        shape = (self.Z.shape + self.T.shape)
        self.density = new_array(self.density, shape)
        self.population = self.population or (None, None)
        self.population = tuple(new_array(state, shape) for state in self.population)
    
    def _solve(self):
        for l in range(self.count_z):
            for k in range(self.count_t):
                idxs = (l, k)
                self.population[0][idxs], self.population[1][idxs] = self._calc_population(*idxs)
                self.density[idxs] = self._calc_density(*idxs)
    
    def _compute(self):
        self._init_arrays()
        self._solve()
        
        density_out = self.density[-1]
        population_final = tuple(state.T[-1] for state in self.population)
        return density_out, population_final
    
    def amplify(self, rho, phi, input_pulse, count_t):
        self.rho = rho
        self.phi = phi
        
        self.input_pulse = input_pulse
        
        if self.sizet != input_pulse.duration or self.count_t != count_t or self.T is None or self.T[-1] != input_pulse.t0 + self.sizet:
            self.sizet = input_pulse.duration
            self.count_t = count_t
            self.T, self.dt = np.linspace(input_pulse.t0, input_pulse.t0 + self.sizet, count_t, retstep=True)
        
        data_out = self._compute()
        return data_out

class ExactAmplifier(PulseAmplifier):
    
    analytical_lower_lifetimes = [0.0, float("inf")]
    
    def __init__(self, *args, **kwargs):
        super(ExactAmplifier, self).__init__(*args, **kwargs)
        
        lower_lifetime = self.active_medium.doping_agent.lower_lifetime
        if lower_lifetime == 0.0:
            four_level = True
        elif math.isinf(lower_lifetime):
            four_level = False
        else:
            raise exc.ModelError("%s only accepts zero or infinite lower state lifetimes" % self.__class__.__name__)
        self.four_level = four_level
    
    @staticmethod
    def min_steps_z(*args, **kwargs):
        return 2
    
    def min_steps_t(self, *args, **kwargs):
        return 2
    
    def _calc_population(self, l, k):
        z = self.Z[l]
        t = self.T[k]
        population = self._exact_population(z, t)
        return population
    
    def _calc_density(self, l, k):
        z = self.Z[l]
        t = self.T[k]
        density = self._exact_density(z, t)
        return density
    
    def _exact_population(self, z, t):
        active_medium = self.active_medium
        light_speed = active_medium.light_speed
        xsection = active_medium.doping_agent.xsection
        four_level = self.four_level
        inversion_integral = active_medium.initial_inversion.inversion_integral(self.rho, self.phi, z)
        density_integral = self.input_pulse.density_integral(t)
        inversion_exp = math.exp(- xsection * inversion_integral)
        inversion_coef = 2.0 if not four_level else 1.0
        density_exp_reciproc = math.exp(- inversion_coef * xsection * light_speed * density_integral)
        initial_inversion = active_medium.initial_inversion.inversion(self.rho, self.phi, z)
        inversion = initial_inversion * inversion_exp * density_exp_reciproc / (1.0 + density_exp_reciproc * (inversion_exp - 1.0))
        upper = (initial_inversion + inversion) / 2.0 if not four_level else inversion
        lower = (initial_inversion - inversion) / 2.0 if not four_level else 0.0
        return (upper, lower)
    
    def _exact_density(self, z, t):
        active_medium = self.active_medium
        input_pulse = self.input_pulse
        light_speed = active_medium.light_speed
        xsection = active_medium.doping_agent.xsection
        inversion_integral = active_medium.initial_inversion.inversion_integral(self.rho, self.phi, z)
        density_integral = input_pulse.density_integral(t)
        input_density = input_pulse.density(t)
        inversion_coef = 2.0 if not self.four_level else 1.0
        density_exp = math.exp(- inversion_coef * xsection * light_speed * density_integral)
        density = input_density / (1.0 - (1.0 - math.exp(- xsection * inversion_integral)) * density_exp)
        return density

class ExactOutputAmplifier(ExactAmplifier):
    
    def _compute(self):
        density_out = np.vectorize(self._exact_density)(self.Z[-1], self.T)
        population_final = np.vectorize(self._exact_population)(self.Z, self.T[-1])
        
        return density_out, population_final

class NumericalAmplifier(PulseAmplifier):
    
    def __init__(self, active_medium, count_z):
        super(NumericalAmplifier, self).__init__(active_medium, count_z)
        
        self.initial_population = None
        self.input_density = None
    
    @staticmethod
    def _min_steps(total_size, max_step_size):
        return int(math.ceil(total_size / max_step_size)) + 1
    
    @staticmethod
    def min_steps_z(active_medium, initial_population=None):
        raise NotImplementedError()
    
    def min_steps_t(self, input_pulse, initial_population=None):
        raise NotImplementedError()
    
    def _solve(self):
        raise NotImplementedError()
    
    def amplify(self, rho, phi, input_pulse, count_t, initial_population=None, input_density=None):
        assert initial_population is None or (len(initial_population) == 2 and [pop.shape == (self.count_z,) for pop in initial_population] == [True] * 2)
        assert input_density is None or input_density.shape == (count_t,)
        
        self.initial_population = initial_population
        self.input_density = input_density
        
        return super(NumericalAmplifier, self).amplify(rho, phi, input_pulse, count_t)

class HybridAmplifier(NumericalAmplifier):
    
    @staticmethod
    def min_steps_z(active_medium, initial_population=None):
        # trapezoid method
        # stability condition: unconditionally unstable for the solved equation
        # positiveness condition: 1.0 - xsection * max_inversion * (step/2.0) > 0
        # monotonically increasing condition: always true when the positiveness condition is true
        if initial_population is None:
            max_inversion = active_medium.initial_inversion.ref_inversion
        else:
            max_inversion = np.amax(initial_population[0])
        xsection = active_medium.doping_agent.xsection
        length = active_medium.length
        dz_sing = 2.0 / (xsection * max_inversion) if max_inversion != 0.0 else 2.0 * length
        dz_max = dz_sing / 2.0
        min_count_z = NumericalAmplifier._min_steps(length, dz_max)
        return min_count_z
    
    def min_steps_t(self, input_pulse, initial_population=None):
        # euler method
        # stability condition: step = -2.0/lambda, where lambda is the minimum (negative) eigenvalue of the Jacobian for the maximum photon density; this condition is weaker than the following
        # positiveness condition: (1.0 - xsection * light_speed * max_density * step >= 0) and (1.0 - (1.0/lower_lifetime) * step >= 0)
        # lower <= upper condition: 1.0 - 2.0 * xsection * light_speed * max_density * step >= 0 (stronger than one of the positiveness conditions, but taken into account under an extra condition)
        # upper monotonically decreasing condition: lower <= upper
        
        active_medium = self.active_medium
        if initial_population is None:
            max_inversion = active_medium.initial_inversion.ref_inversion
        else:
            max_inversion = np.amax(initial_population[0])
        light_speed = active_medium.light_speed
        xsection = active_medium.doping_agent.xsection
        duration = input_pulse.duration
        max_density = input_pulse.ref_density * ((1.0 + xsection*max_inversion*self.dz/2.0) / (1.0 - xsection*max_inversion*self.dz/2.0))**(self.count_z-1)
        
        # n2 positiveness
        dt_neg_n2 = 1.0 / (xsection * light_speed * max_density) if max_density != 0.0 else 2.0 * duration
        
        # n1 positiveness
        lower_lifetime = active_medium.doping_agent.lower_lifetime
        dt_neg_n1 = lower_lifetime if lower_lifetime != 0.0 else 2.0 * duration
        
        # n1 <= n2
        dt_ratio = 1.0 / (2.0 * xsection * light_speed * max_density) if lower_lifetime != 0.0 and max_density != 0.0 else 2.0 * duration
        
        dt_max = min(dt_neg_n2, dt_neg_n1, dt_ratio) / 2.0
        min_count_t = NumericalAmplifier._min_steps(duration, dt_max)
        return min_count_t
    
    def _solve(self):
        native._hybrid_amplifier_solve(self)

class NSFDAmplifier(NumericalAmplifier):
    
    @staticmethod
    def min_steps_z(*args, **kwargs):
        return 2
    
    def min_steps_t(self, *args, **kwargs):
        return 2
    
    def _solve(self):
        native._nsfd_amplifier_solve(self)
