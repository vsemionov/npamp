
"""Provides a sawtooth temporal pulse shape"""


import pamp


class SawtoothPulse(pamp.pulse.SinglePulse):
    
    def __init__(self, *args, **kwargs):
        super(SawtoothPulse, self).__init__(*args, **kwargs)
        self.slope = self.ref_density / self.duration
    
    def density(self, t):
        t -= self.t0
        return self.slope * t if 0.0 <= t <= self.duration else 0.0
    
    def density_integral(self, t):
        t -= self.t0
        t = max(t, 0.0)
        t = min(t, self.duration)
        return self.slope * t**2.0 / 2.0


setattr(pamp.pulse, SawtoothPulse.__name__, SawtoothPulse)
