
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


class InputPulse(object):
    
    def __init__(self, t0, duration, ref_density):
        self.t0 = t0
        self.duration = duration
        self.ref_density = ref_density
        
        self.offset = t0 + duration/2.0
    
    def density(self, t):
        raise NotImplementedError()
    
    def density_integral(self, t):
        raise NotImplementedError()

class SinglePulse(InputPulse):
    
    @classmethod
    def ref_density(cls, light_speed, duration, fluence):
        ref_density = 1.0
        test_pulse = cls(-duration/2.0, duration, ref_density)
        test_fluence = light_speed * test_pulse.density_integral(float("inf"))
        ref_density *= fluence / test_fluence
        return ref_density

class SquarePulse(SinglePulse):
    
    def density(self, t):
        t -= self.t0
        return self.ref_density if 0.0 <= t <= self.duration else 0.0
    
    def density_integral(self, t):
        t -= self.t0
        t = max(t, 0.0)
        t = min(t, self.duration)
        return self.ref_density * t

class LorentzianPulse(SinglePulse):
    
    def __init__(self, *args, **kwargs):
        super(LorentzianPulse, self).__init__(*args, **kwargs)
        self.gamma = self.duration/2.0
        self.coef = math.pi * self.gamma
    
    def density(self, t):
        return self.ref_density * self.coef * (1.0/math.pi) * (self.gamma / ((t - self.offset)**2.0 + self.gamma**2.0))
    
    def density_integral(self, t):
        return self.ref_density * self.coef * ((1.0/math.pi) * math.atan((t - self.offset)/self.gamma) + 1.0/2.0)

class GaussianPulse(SinglePulse):
    
    def __init__(self, *args, **kwargs):
        super(GaussianPulse, self).__init__(*args, **kwargs)
        self.sigma = self.duration/(2.0*math.sqrt(2.0*math.log(2.0)))
    
    def density(self, t):
        return self.ref_density * math.exp(- (t - self.offset)**2.0 / (2.0 * self.sigma**2.0))
    
    def density_integral(self, t):
        return self.ref_density * self.sigma * math.sqrt(math.pi/2.0) * (math.erf((t - self.offset) / (math.sqrt(2.0)*self.sigma)) + 1.0)

class TransformedPulse(InputPulse):
    
    def __init__(self, pulse):
        super(TransformedPulse, self).__init__(pulse.t0, pulse.duration, pulse.ref_density)
        self.pulse = pulse
    
    def density(self, t):
        return self.pulse.density(t)
    
    def density_integral(self, t):
        return self.pulse.density_integral(t)

class ScaledPulse(TransformedPulse):
    
    def __init__(self, pulse, ref_density):
        super(ScaledPulse, self).__init__(pulse)
        self.ref_density = ref_density
        self.density_scale = self.ref_density / pulse.ref_density
    
    def density(self, t):
        return self.density_scale * self.pulse.density(t)
    
    def density_integral(self, t):
        return self.density_scale * self.pulse.density_integral(t)

class ExtendedPulse(TransformedPulse):
    
    def __init__(self, pulse, scale):
        super(ExtendedPulse, self).__init__(pulse)
        self.duration_scale = scale
        self.duration *= scale
        self.t0 = self.offset - self.duration/2.0

class TruncatedPulse(TransformedPulse):
    
    def __init__(self, pulse):
        super(TruncatedPulse, self).__init__(pulse)
    
    def density(self, t):
        return self.pulse.density(t) if 0.0 <= (t - self.t0) <= self.duration else 0.0
    
    def density_integral(self, t):
        t = min(max(t, self.t0), self.t0 + self.duration)
        return self.pulse.density_integral(t) - self.pulse.density_integral(self.t0)

class ShiftedPulse(TransformedPulse):
    
    def __init__(self, pulse, t0):
        super(ShiftedPulse, self).__init__(pulse)
        self.shift = t0 - pulse.t0
    
    def density(self, t):
        return self.pulse.density(t - self.shift)
    
    def density_integral(self, t):
        return self.pulse.density_integral(t - self.shift)

class PartialPulse(TransformedPulse):
    
    def __init__(self, pulse, t0_integral):
        super(PartialPulse, self).__init__(pulse)
        self.t0_integral = t0_integral
    
    def density_integral(self, t):
        return self.t0_integral + self.pulse.density_integral(t) - self.pulse.density_integral(self.t0)

class CompositePulse(InputPulse):
    
    def __init__(self, pulses):
        t0 = pulses[0].t0
        duration = sum(map(lambda p: p.duration, pulses))
        ref_density = max(map(lambda p: p.ref_density, pulses))
        super(CompositePulse, self).__init__(t0, duration, ref_density)
        
        self._pmap_t = np.empty(len(pulses))
        self._pmap_n = np.array(range(len(pulses)))
        
        self.pulses = [pulses[0]]
        self.shifted_pulses = [ShiftedPulse(pulses[0], self.t0)]
        self._pmap_t[0] = pulses[0].t0
        for i in range(1, len(pulses)):
            p = pulses[i]
            ppold = self.pulses[i-1]
            t0_integral = ppold.density_integral(ppold.t0 + ppold.duration)
            pulse = PartialPulse(p, t0_integral)
            self.pulses.append(pulse)
            t0 = self._pmap_t[i-1] + ppold.duration
            shifted_pulse = ShiftedPulse(pulse, t0)
            self.shifted_pulses.append(shifted_pulse)
            self._pmap_t[i] = t0
    
    def _pmap(self, t):
        n = int(math.floor(np.interp(t, self._pmap_t, self._pmap_n)))
        p = self.shifted_pulses[n]
        return p
    
    def density(self, t):
        p = self._pmap(t)
        return p.density(t)
    
    def density_integral(self, t):
        p = self._pmap(t)
        return p.density_integral(t)

class PulseTrain(CompositePulse):
    
    def __init__(self, pulse, count, period):
        self.pulse = pulse
        self.count = count
        self.period = period
        
        self.separator = SquarePulse(0.0, period - pulse.duration, 0.0)
        
        pulses = []
        for _ in range(count-1):
            pulses.append(pulse)
            pulses.append(self.separator)
        pulses.append(pulse)
        
        super(PulseTrain, self).__init__(pulses)
