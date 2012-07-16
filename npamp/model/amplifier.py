
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
        
        self.Z = np.linspace(0.0, active_medium.length, count_z)
        self.T = None
        
        self.density = None
        self.population = None
        
        self._input_pulse = None
        self._duration = None
    
    @staticmethod
    def min_steps_z(active_medium):
        raise NotImplementedError()
    
    def min_steps_t(self, input_pulse):
        raise NotImplementedError()
    
    def _calc_population(self, rho, phi, l, k):
        raise NotImplementedError()
    
    def _calc_density(self, rho, phi, l, k):
        raise NotImplementedError()
    
    def _solve(self, rho, phi):
        for l in range(len(self.Z)):
            for k in range(len(self.T)):
                self.population[0][l, k], self.population[1][l, k] = self._calc_population(rho, phi, l, k)
                self.density[l, k] = self._calc_density(rho, phi, l, k)
    
    def _init_arrays(self):
        new_array = lambda old, shape: np.empty(shape) if old is None or old.shape != shape else old
        shape = (self.Z.shape + self.T.shape)
        self.density = new_array(self.density, shape)
        self.population = self.population or (None, None)
        self.population = tuple(new_array(state, shape) for state in self.population)
    
    def _compute(self, rho, phi):
        self._init_arrays()
        self._solve(rho, phi)
        
        density_out = self.density[-1]
        population_final = tuple(state.T[-1] for state in self.population)
        return density_out, population_final
    
    def amplify(self, rho, phi, input_pulse, count_t):
        self._input_pulse = input_pulse
        
        T = self.T
        if T is None or len(T) != count_t or T[0] != input_pulse.t0 or self._duration != input_pulse.duration:
            self.T = np.linspace(input_pulse.t0, input_pulse.t0 + input_pulse.duration, count_t)
            self._duration = input_pulse.duration
        
        data_out = self._compute(rho, phi)
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
    
    def _calc_population(self, rho, phi, l, k):
        z, t = self.Z[l], self.T[k]
        population = self._exact_population(rho, phi, z, t)
        return population
    
    def _calc_density(self, rho, phi, l, k):
        z, t = self.Z[l], self.T[k]
        density = self._exact_density(rho, phi, z, t)
        return density
    
    def _exact_population(self, rho, phi, z, t):
        active_medium = self.active_medium
        light_speed = active_medium.light_speed
        xsection = active_medium.doping_agent.xsection
        four_level = self.four_level
        inversion_integral = active_medium.initial_inversion.inversion_integral(rho, phi, z)
        density_integral = self._input_pulse.density_integral(t)
        inversion_exp = math.exp(- xsection * inversion_integral)
        inversion_coef = 2.0 if not four_level else 1.0
        density_exp_reciproc = math.exp(- inversion_coef * xsection * light_speed * density_integral)
        initial_inversion = active_medium.initial_inversion.inversion(rho, phi, z)
        inversion = initial_inversion * inversion_exp * density_exp_reciproc / (1.0 + density_exp_reciproc * (inversion_exp - 1.0))
        upper = (initial_inversion + inversion) / 2.0 if not four_level else inversion
        lower = (initial_inversion - inversion) / 2.0 if not four_level else 0.0
        return (upper, lower)
    
    def _exact_density(self, rho, phi, z, t):
        active_medium = self.active_medium
        input_pulse = self._input_pulse
        light_speed = active_medium.light_speed
        xsection = active_medium.doping_agent.xsection
        inversion_integral = active_medium.initial_inversion.inversion_integral(rho, phi, z)
        density_integral = input_pulse.density_integral(t)
        input_density = input_pulse.density(t)
        inversion_coef = 2.0 if not self.four_level else 1.0
        density_exp = math.exp(- inversion_coef * xsection * light_speed * density_integral)
        density = input_density / (1.0 - (1.0 - math.exp(- xsection * inversion_integral)) * density_exp)
        return density

class ExactOutputAmplifier(ExactAmplifier):
    
    def _compute(self, rho, phi):
        Z, T = self.Z, self.T
        density_out = np.vectorize(self._exact_density)(rho, phi, Z[-1], T)
        population_final = np.vectorize(self._exact_population)(rho, phi, Z, T[-1])
        
        return density_out, population_final

class NumericalAmplifier(PulseAmplifier):
    
    def __init__(self, active_medium, count_z):
        super(NumericalAmplifier, self).__init__(active_medium, count_z)
        
        self._initial_population = None
        self._input_density = None
    
    @staticmethod
    def _min_steps(total_size, max_step_size):
        return int(max(math.ceil(total_size / max_step_size), 1.0)) + 1
    
    @staticmethod
    def min_steps_z(active_medium, initial_population=None):
        raise NotImplementedError()
    
    def min_steps_t(self, input_pulse, initial_population=None):
        raise NotImplementedError()
    
    def _solve(self, rho, phi):
        raise NotImplementedError()
    
    def amplify(self, rho, phi, input_pulse, count_t, initial_population=None, input_density=None):
        assert initial_population is None or (len(initial_population) == 2 and [state.shape == self.Z.shape for state in initial_population] == [True] * 2)
        assert input_density is None or input_density.shape == (count_t,)
        
        self._initial_population = initial_population
        self._input_density = input_density
        
        return super(NumericalAmplifier, self).amplify(rho, phi, input_pulse, count_t)

class HybridAmplifier(NumericalAmplifier):
    
    @staticmethod
    def min_steps_z(active_medium, initial_population=None):
        # trapezoid method
        # stability condition: unconditionally unstable for the solved equation
        # positiveness condition: 1.0 - xsection * max_inversion * (step/2.0) > 0
        # monotonically increasing condition: always true when the positiveness condition is true
        
        xsection = active_medium.doping_agent.xsection
        length = active_medium.length
        
        if initial_population is None:
            max_inversion = active_medium.initial_inversion.ref_inversion
        else:
            max_inversion = np.amax(initial_population[0])
        
        dz_sing = 2.0 / (xsection * max_inversion) if max_inversion != 0.0 else float("inf")
        dz_max = dz_sing / 2.0
        min_count_z = NumericalAmplifier._min_steps(length, dz_max)
        
        return min_count_z
    
    def min_steps_t(self, input_pulse, initial_population=None):
        # euler method
        # stability condition: step = -2.0/lambda, where lambda is the minimum (negative) eigenvalue of the Jacobian for the maximum photon density; this condition is weaker than the following set
        # positiveness conditions: (1.0 - xsection * light_speed * max_density * step >= 0) and (1.0 - (1.0/lower_lifetime) * step >= 0)
        # upper >= lower condition: 1.0 - 2.0 * xsection * light_speed * max_density * step >= 0 (stronger than one of the positiveness conditions, but taken into account under an extra condition)
        # upper monotonically decreasing condition: upper >= lower
        
        active_medium = self.active_medium
        doping_agent = active_medium.doping_agent
        light_speed = active_medium.light_speed
        xsection = doping_agent.xsection
        duration = input_pulse.duration
        
        Z = self.Z
        count_z = len(Z)
        dz = (Z[-1] - Z[0]) / (count_z - 1)
        
        if initial_population is None:
            max_inversion = active_medium.initial_inversion.ref_inversion
        else:
            max_inversion = np.amax(initial_population[0])
        
        max_density = input_pulse.ref_density * ((1.0 + xsection*max_inversion*dz/2.0) / (1.0 - xsection*max_inversion*dz/2.0))**(count_z-1)
        
        # n2 positiveness
        dt_neg_n2 = 1.0 / (xsection * light_speed * max_density) if max_density != 0.0 else float("inf")
        
        # n1 positiveness
        lower_lifetime = doping_agent.lower_lifetime
        dt_neg_n1 = lower_lifetime if lower_lifetime != 0.0 else float("inf")
        
        # n positiveness
        dt_neg_n = 1.0 / (2.0 * xsection * light_speed * max_density) if lower_lifetime != 0.0 and max_density != 0.0 else float("inf")
        
        dt_max = min(dt_neg_n2, dt_neg_n1, dt_neg_n) / 2.0
        min_count_t = NumericalAmplifier._min_steps(duration, dt_max)
        
        return min_count_t
    
    def _solve(self, rho, phi):
        native._hybrid_amplifier_solve(self, rho, phi)

class NSFDAmplifier(NumericalAmplifier):
    
    @staticmethod
    def min_steps_z(*args, **kwargs):
        return 2
    
    def min_steps_t(self, *args, **kwargs):
        return 2
    
    def _solve(self, rho, phi):
        native._nsfd_amplifier_solve(self, rho, phi)
