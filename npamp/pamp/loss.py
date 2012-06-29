
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


# All the ASE models assume a four-level system by considering the inversion equivalent to the upper state population, which is valid only when the lower level population is zero.
# A negative inversion is thus equivalent to a negative upper state population, which is not physically sane.
# A zero inversion is equivalent to a zero upper state population, which naturally yields zero fluorescence and ASE loss rate.
# Having these in mind, the loss rate functions, as per these ASE models, are only defined in the domain of non-negative inversions.


def decay_rate(inversion, volume, lifetime):
    rate = inversion * volume / lifetime
    return rate


def concrete(cls):
    setattr(cls, "concrete", True)
    return cls

class LossModel(object):
    descr = "loss model"
    
    def __init__(self, active_medium, wavelen):
        self.active_medium = active_medium
        self.wavelen = wavelen
    
    def rate(self, inversion):
        if inversion < 0.0:
            raise ValueError("negative population inversion")
        return self.calc_rate(inversion)
    
    def calc_rate(self, inversion):
        raise NotImplementedError()

class ExactLossModel(LossModel):
    descr = "exact loss model"
    
    def __init__(self, active_medium, wavelen):
        super(ExactLossModel, self).__init__(active_medium, wavelen)
        
        self.photon_energy = energy.photon_energy(self.wavelen)
        self.active_solid_angle = self.active_medium.aperture / self.active_medium.length**2.0
        self.radiative_lifetime = active_medium.doping_agent.upper_lifetime / active_medium.doping_agent.branching_ratio
        self.saturation_intensity = self.photon_energy / (self.active_medium.doping_agent.xsection * self.radiative_lifetime)

class NumericalLossModel(LossModel):
    descr = "numerical loss model"
    
    def __init__(self, active_medium, wavelen, rtol, min_count):
        super(NumericalLossModel, self).__init__(active_medium, wavelen)
        
        self.rtol = rtol
        self.min_count = min_count

class PerturbedLossModel(LossModel):
    descr = "perturbed loss model"
    
    def __init__(self, numerical_loss_model):
        super(PerturbedLossModel, self).__init__(numerical_loss_model.active_medium, numerical_loss_model.wavelen)
        
        self.numerical_loss_model = numerical_loss_model
        self.perturb = numerical_loss_model.rtol
    
    def calc_rate(self, inversion):
        rate = self.numerical_loss_model.calc_rate(inversion)
        rate *= (1.0 - self.perturb)
        return rate

@concrete
class NullLossModel(ExactLossModel):
    descr = "null loss model"
    
    def calc_rate(self, inversion):
        return 0.0

@concrete
class FluorescenceLossModel(ExactLossModel):
    descr = "fluorescence"
    
    def calc_rate(self, inversion):
        active_medium = self.active_medium
        rate = decay_rate(inversion, active_medium.volume, active_medium.doping_agent.upper_lifetime)
        return rate

@concrete
class LinfordASEModel(ExactLossModel):
    descr = "Linford ASE model"
    
    def calc_rate(self, inversion):
        # at zero inversion, loss rate is undefined but approaches zero
        if inversion == 0.0:
            return 0.0
        active_medium = self.active_medium
        ln_gain = active_medium.doping_agent.xsection * inversion * active_medium.length
        gain = math.exp(ln_gain)
        intensity = self.saturation_intensity * (self.active_solid_angle/4.0) * (gain-1.0)**1.5/(gain*ln_gain)**0.5
        power = 2.0 * active_medium.aperture * intensity
        rate = energy.photon_count(self.wavelen, power)
        # fluorescence in directions outside of the two active solid angles:
        rate += (1.0 - 2.0 * self.active_solid_angle / (4.0 * math.pi)) * decay_rate(inversion, active_medium.volume, self.radiative_lifetime)
        # fluorescence to other levels:
        rate += decay_rate(inversion, active_medium.volume, active_medium.doping_agent.upper_lifetime / (1.0 - active_medium.doping_agent.branching_ratio))
        return rate

@concrete
class SchulzASEModel(ExactLossModel):
    descr = "Schulz ASE model"
    
    def calc_rate(self, inversion):
        active_medium = self.active_medium
        xsection = active_medium.doping_agent.xsection
        ln_gain = xsection * inversion * active_medium.length
        gain = math.exp(ln_gain)
        gamma_prime = ln_gain / active_medium.doping_agent.upper_lifetime + (self.active_solid_angle * (gain - ln_gain - 1.0) / (2.0 * math.pi)) / self.radiative_lifetime
        rate = gamma_prime * active_medium.aperture / xsection
        return rate

@concrete
class RossApproximateASEModel(ExactLossModel):
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
        active_medium = self.active_medium
        B = self._B(inversion)
        loss_coef = 1.0 + B
        rate = inversion * loss_coef / active_medium.doping_agent.upper_lifetime
        rate *= active_medium.volume
        return rate

@concrete
class RossNumericalASEModel(NumericalLossModel):
    descr = "Ross ASE model (numerical)"
    
    _integrate_B = lambda self, *args: native._ross_ase_model_integrate_B(self, *args)
    
    def __init__(self, active_medium, wavelen, rtol, min_count, sample_count_multiplier=16):
        super(RossNumericalASEModel, self).__init__(active_medium, wavelen, rtol, min_count)
        
        self.sample_count_multiplier = sample_count_multiplier
        
        self.max_nsamples = 16 * 1024**2
    
    def _loss_coef(self, inversion):
        nsamples = self.min_count
        while True:
            B, B_rel_error = self._integrate_B(inversion, nsamples)
            B_abs_error = B_rel_error * B
            loss_coef = 1.0 + B
            loss_rel_error = B_abs_error / loss_coef
            if loss_rel_error < self.rtol:
                break
            nsamples *= self.sample_count_multiplier
            if nsamples > self.max_nsamples:
                warnings.warn("max_nsamples (%d) exceeded; latest rel. error: %f" % (self.max_nsamples, loss_rel_error), stacklevel=2)
                break
        return loss_coef
    
    def calc_rate(self, inversion):
        active_medium = self.active_medium
        loss_coef = self._loss_coef(inversion)
        rate = inversion * loss_coef / active_medium.doping_agent.upper_lifetime
        rate *= active_medium.volume
        return rate
    
    def rate_rel_stddev(self, inversion):
        nsamples = nsubsamples = int(math.ceil(math.sqrt(self.max_nsamples)))
        B, _ = self._integrate_B(inversion, self.max_nsamples)
        B_stddev = native._ross_ase_model_B_stddev(self, inversion, nsamples, nsubsamples)
        loss_coef = 1.0 + B
        rate_rel_stddev = B_stddev / loss_coef
        return rate_rel_stddev

@concrete
class RossHybridASEModel(RossNumericalASEModel):
    descr = "Ross ASE model (hybrid)"
    
    _integrate_B = lambda self, *args: native._ross_ase_model_integrate_BP(self, 0.0, 0.0, 0.0, *args)
