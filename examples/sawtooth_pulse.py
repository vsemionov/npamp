
"""Provides a sawtooth temporal pulse shape"""


import model


class SawtoothPulse(model.pulse.SinglePulse):
    
    def density(self, t):
        t -= self.t0
        return (self.ref_density / self.duration) * t if 0.0 <= t <= self.duration else 0.0
    
    def density_integral(self, t):
        t -= self.t0
        t = max(t, 0.0)
        t = min(t, self.duration)
        return (self.ref_density / self.duration) * t**2.0 / 2.0


model.pulse.SawtoothPulse = SawtoothPulse
