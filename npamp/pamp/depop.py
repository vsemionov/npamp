
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
import warnings

import energy
import native
import error


# All the ASE models assume a four-level system by considering the inversion equivalent to the upper state population, which is valid only when the lower state population is zero.
# A negative inversion is thus equivalent to a negative upper state population, which is not physically sane.
# A zero inversion is equivalent to a zero upper state population, which naturally yields zero fluorescence and ASE depopulation rate.
# Having these in mind, the depopulation rate functions, as per these ASE models, are only defined in the domain of non-negative inversions.


def decay_rate(inversion, volume, probability):
    rate = inversion * volume * probability
    return rate


def concrete(cls):
    setattr(cls, "concrete", True)
    return cls

class DepopulationModel(object):
    descr = "depopulation model"
    
    def __init__(self, active_medium, wavelen):
        self.active_medium = active_medium
        self.wavelen = wavelen
    
    def rate(self, inversion):
        if inversion < 0.0:
            raise ValueError("negative population inversion")
        return self.calc_rate(inversion)
    
    def calc_rate(self, inversion):
        raise NotImplementedError()

class ExactDepopulationModel(DepopulationModel):
    descr = "exact depopulation model"
    
    def __init__(self, active_medium, wavelen):
        super(ExactDepopulationModel, self).__init__(active_medium, wavelen)
        
        self.photon_energy = energy.photon_energy(self.wavelen)
        self.active_solid_angle = self.active_medium.aperture / self.active_medium.length**2.0
        self.upper_probability = 1.0 / active_medium.doping_agent.upper_lifetime
        self.radiative_probability = self.upper_probability * active_medium.doping_agent.branching_ratio
        self.nonradiative_probability = self.upper_probability * (1.0 - active_medium.doping_agent.branching_ratio)
        self.saturation_intensity = self.photon_energy * self.radiative_probability / self.active_medium.doping_agent.xsection

class NumericalDepopulationModel(DepopulationModel):
    descr = "numerical depopulation model"
    
    def __init__(self, active_medium, wavelen, rtol, min_count):
        super(NumericalDepopulationModel, self).__init__(active_medium, wavelen)
        
        self.rtol = rtol
        self.min_count = min_count

class PerturbedDepopulationModel(DepopulationModel):
    descr = "perturbed depopulation model"
    
    def __init__(self, numerical_depop_model):
        super(PerturbedDepopulationModel, self).__init__(numerical_depop_model.active_medium, numerical_depop_model.wavelen)
        
        self.numerical_depop_model = numerical_depop_model
        self.perturb = numerical_depop_model.rtol
    
    def calc_rate(self, inversion):
        rate = self.numerical_depop_model.calc_rate(inversion)
        rate *= (1.0 - self.perturb)
        return rate

@concrete
class NullDepopulation(ExactDepopulationModel):
    descr = "null depopulation"
    
    def calc_rate(self, inversion):
        return 0.0

@concrete
class FluorescenceModel(ExactDepopulationModel):
    descr = "fluorescence"
    
    def calc_rate(self, inversion):
        rate = decay_rate(inversion, self.active_medium.volume, self.upper_probability)
        return rate

@concrete
class LinfordASEModel(ExactDepopulationModel):
    descr = "Linford ASE model"
    
    def calc_rate(self, inversion):
        # at zero inversion, depopulation rate is undefined but approaches zero
        if inversion == 0.0:
            return 0.0
        active_medium = self.active_medium
        ln_gain = active_medium.doping_agent.xsection * inversion * active_medium.length
        gain = math.exp(ln_gain)
        intensity = self.saturation_intensity * (self.active_solid_angle/4.0) * (gain-1.0)**1.5/(gain*ln_gain)**0.5
        power = 2.0 * active_medium.aperture * intensity
        rate = energy.photon_count(self.wavelen, power)
        # fluorescence in directions outside of the two active solid angles:
        rate += (1.0 - 2.0 * self.active_solid_angle / (4.0 * math.pi)) * decay_rate(inversion, active_medium.volume, self.radiative_probability)
        # fluorescence to other states:
        rate += decay_rate(inversion, active_medium.volume, self.nonradiative_probability)
        return rate

@concrete
class SchulzASEModel(ExactDepopulationModel):
    descr = "Schulz ASE model"
    
    def calc_rate(self, inversion):
        active_medium = self.active_medium
        xsection = active_medium.doping_agent.xsection
        ln_gain = xsection * inversion * active_medium.length
        gain = math.exp(ln_gain)
        gamma_prime = ln_gain * self.upper_probability + (self.active_solid_angle * (gain - ln_gain - 1.0) / (2.0 * math.pi)) * self.radiative_probability
        rate = gamma_prime * active_medium.aperture / xsection
        return rate

@concrete
class RossApproximateASEModel(ExactDepopulationModel):
    descr = "Ross ASE model (approximate)"
    
    _Nd_YLF_branching_ratio = 0.56 * 0.525 / 0.500
    
    def _B(self, inversion):
        active_medium = self.active_medium
        length = active_medium.length
        ln_gain = active_medium.doping_agent.xsection * inversion * length
        gain = math.exp(ln_gain)
        B = 0.05 * (2.0*active_medium.radius/length)**0.3 * gain
        B *= self.active_medium.doping_agent.branching_ratio / self._Nd_YLF_branching_ratio
        return B
    
    def calc_rate(self, inversion):
        B = self._B(inversion)
        rate_coef = 1.0 + B
        rate = inversion * rate_coef * self.upper_probability
        rate *= self.active_medium.volume
        return rate

@concrete
class RossNumericalASEModel(NumericalDepopulationModel):
    descr = "Ross ASE model (numerical)"
    
    _integrate_B = lambda self, *args: native._ross_ase_model_integrate_B(self, *args)
    
    def __init__(self, active_medium, wavelen, rtol, min_count, sample_count_multiplier=16):
        super(RossNumericalASEModel, self).__init__(active_medium, wavelen, rtol, min_count)
        
        self.sample_count_multiplier = sample_count_multiplier
        
        self.max_nsamples = 16 * 1024**2
    
    def _rate_coef(self, inversion):
        if self.min_count > self.max_nsamples:
            raise error.ModelError("min. sample count (%d) is greater than max. number of samples (%d)" % (self.min_count, self.max_nsamples))
        nsamples = self.min_count
        while True:
            B, B_rel_error = self._integrate_B(inversion, nsamples)
            B_abs_error = B_rel_error * B
            rate_coef = 1.0 + B
            rate_rel_error = B_abs_error / rate_coef
            if rate_rel_error < self.rtol:
                break
            nsamples *= self.sample_count_multiplier
            if nsamples > self.max_nsamples:
                warnings.warn("max_nsamples (%d) exceeded; latest rel. error: %f" % (self.max_nsamples, rate_rel_error), stacklevel=2)
                break
        return rate_coef
    
    def calc_rate(self, inversion):
        rate_coef = self._rate_coef(inversion)
        rate = inversion * rate_coef * self.upper_probability
        rate *= self.active_medium.volume
        return rate
    
    def rate_rel_stddev(self, inversion):
        nsamples = nsubsamples = int(math.ceil(math.sqrt(self.max_nsamples)))
        B, _ = self._integrate_B(inversion, self.max_nsamples)
        B_stddev = native._ross_ase_model_B_stddev(self, inversion, nsamples, nsubsamples)
        rate_coef = 1.0 + B
        rate_rel_stddev = B_stddev / rate_coef
        return rate_rel_stddev

@concrete
class RossHybridASEModel(RossNumericalASEModel):
    descr = "Ross ASE model (hybrid)"
    
    _integrate_B = lambda self, *args: native._ross_ase_model_integrate_BP(self, 0.0, 0.0, self.z_rel * self.active_medium.length, *args)
    
    def __init__(self, active_medium, wavelen, rtol, min_count, sample_count_multiplier=16, z_rel=0.0):
        assert 0.0 <= z_rel <= 1.0, "z_rel is not in the interval [0.0, 1.0]"
        
        super(RossHybridASEModel, self).__init__(active_medium, wavelen, rtol, min_count, sample_count_multiplier)
        
        self.z_rel = z_rel
