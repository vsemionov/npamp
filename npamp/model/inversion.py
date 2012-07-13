
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


class InitialInversion(object):
    
    xcoords = True, True, True
    rho_ref, phi_ref = 0.0, 0.0
    
    def __init__(self, ref_inversion):
        self.ref_inversion = ref_inversion
    
    def inversion(self, rho, phi, z):
        raise NotImplementedError()
    
    def inversion_integral(self, rho, phi, z):
        raise NotImplementedError()
    
    def inversion_fluence_integral(self, rho, z):
        if self.xcoords[0:2] == (False, False):
            return self.inversion_integral(rho, 0.0, z) * math.pi * rho**2.0
        else:
            raise NotImplementedError()

class UniformInversion(InitialInversion):
    
    xcoords = False, False, False
    
    def inversion(self, rho, phi, z):
        return self.ref_inversion
    
    def inversion_integral(self, rho, phi, z):
        return self.ref_inversion * z

class SquareInversion(InitialInversion):
    
    xcoords = False, False, True
    
    def __init__(self, ref_inversion, z0, length):
        super(SquareInversion, self).__init__(ref_inversion)
        assert z0 >= 0.0
        self.z0 = z0
        self.length = length
    
    def inversion(self, rho, phi, z):
        z -= self.z0
        return self.ref_inversion if 0.0 <= z <= self.length else 0.0
    
    def inversion_integral(self, rho, phi, z):
        z -= self.z0
        z = max(z, 0.0)
        z = min(z, self.length)
        return self.ref_inversion * z
