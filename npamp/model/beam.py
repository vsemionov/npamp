
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


class BeamProfile(object):
    
    xcoords = False, False
    rho_ref, phi_ref = 0.0, 0.0
    rho_trunc = float("inf")
    
    def __init__(self, radius, ref_fluence):
        self.radius = radius
        self.ref_fluence = ref_fluence
    
    @classmethod
    def ref_fluence(cls, radius, photon_count):
        ref_fluence = 1.0
        test_beam = cls(radius, ref_fluence)
        test_photon_count = test_beam.fluence_integral(float("inf"))
        ref_fluence *= photon_count / test_photon_count
        return ref_fluence
    
    def fluence(self, rho, phi):
        raise NotImplementedError()
    
    def fluence_integral(self, rho):
        raise NotImplementedError()

class TophatBeam(BeamProfile):
    
    xcoords = False, False
    
    def __init__(self, *args, **kwargs):
        super(TophatBeam, self).__init__(*args, **kwargs)
        self.rho_trunc = self.radius
    
    def fluence(self, rho, phi):
        return self.ref_fluence if rho <= self.radius else 0.0
    
    def fluence_integral(self, rho):
        rho = min(rho, self.radius)
        return self.ref_fluence * math.pi * rho**2.0

class GaussianBeam(BeamProfile):
    
    xcoords = True, False
    
    def __init__(self, *args, **kwargs):
        super(GaussianBeam, self).__init__(*args, **kwargs)
        self.sigma = self.radius/2.0
    
    def fluence(self, rho, phi):
        return self.ref_fluence * math.exp(- rho**2.0 / (2.0 * self.sigma**2.0))
    
    def fluence_integral(self, rho):
        return self.ref_fluence * self.sigma**2.0 * 2.0*math.pi * (1.0 - math.exp(- rho**2.0 / (2.0 * self.sigma**2.0)))
