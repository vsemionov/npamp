
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
import scipy.integrate


class NumericalIntegrator(object):
    min_count = 2
    method = None
    
    def __call__(self, Y, dx):
        return self.method(Y, dx=dx)

class TrapezoidIntegrator(NumericalIntegrator):
    method = staticmethod(scipy.integrate.trapz)

class SimpsonIntegrator(NumericalIntegrator):
    method = staticmethod(scipy.integrate.simps)

class RombergIntegrator(NumericalIntegrator):
    min_count = 3 # scipy.integrate.romb() seems to require 3 or more samples
    method = staticmethod(scipy.integrate.romb)


class DomainIntegrator(object):
    
    def __init__(self, int_type, active_medium):
        self.active_medium = active_medium
        self.num_integrator = int_type()
    
    def integrate(self, X, Y):
        assert X.ndim == Y.ndim == 1
        assert X.shape == Y.shape
        assert len(X) >= self.num_integrator.min_count
        
        divs_f = math.log(len(X) - 1, 2.0)
        divs = int(divs_f)
        assert divs == divs_f
        
        dx = (X[-1] - X[0]) / (len(X) - 1)
        I = self.num_integrator(Y, dx)
        return I
    
    def integrate_base(self, Rho, Phi, fluence):
        assert Rho.ndim == Phi.ndim == 1
        assert fluence.shape == (Rho.shape + Phi.shape)
        
        integrate = lambda Y, X, xmax: self.integrate(X, Y) if len(X) > 1 else xmax*Y[0]
        
        phi_integrals = np.apply_along_axis(integrate, 1, fluence, Phi, 2.0*math.pi)
        radius = self.active_medium.radius
        phi_integrals *= Rho if len(Rho) > 1 else radius/2.0
        rho_integral = integrate(phi_integrals, Rho, radius)
        
        return rho_integral
