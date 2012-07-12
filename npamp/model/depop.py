
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

import native
import energy
import util
import exc


# All the ASE models assume a four-level system by considering the inversion equivalent to the upper state population, which is valid only when the lower state population is zero.
# A negative inversion is thus equivalent to a negative upper state population, which is not physically sane.
# A zero inversion is equivalent to a zero upper state population, which naturally yields zero fluorescence and ASE depopulation rate.
# Having these in mind, the depopulation rate functions, as per these ASE models, are only defined in the domain of non-negative inversions.


def concrete(cls):
    setattr(cls, "concrete", True)
    return cls

class DepopulationModel(object):
    descr = "depopulation model"
    
    def __init__(self, active_medium):
        self.active_medium = active_medium
    
    def _decay_rate(self, inversion, probability):
        rate = inversion * self.active_medium.volume * probability
        return rate
    
    def _rate(self, inversion):
        raise NotImplementedError()
    
    def rate(self, inversion):
        if inversion < 0.0:
            raise exc.ModelError("negative population inversion")
        return self._rate(inversion)

class AnalyticalDepopulationModel(DepopulationModel):
    descr = "analytical depopulation model"

class NumericalDepopulationModel(DepopulationModel):
    descr = "numerical depopulation model"
    
    def __init__(self, active_medium, rtol, min_samples, prng_seed=-1):
        super(NumericalDepopulationModel, self).__init__(active_medium)
        
        self.max_samples = 16 * 1024**2
        
        self.rtol = rtol
        self.min_samples = min_samples
        self.prng_seed = prng_seed
        
        if self.min_samples > self.max_samples:
            raise exc.NumericalError("min. depop. rate sample size (%s) is greater than max. sample size (%s)" % (self.min_samples, self.max_samples))
        
        native._seed_prng(prng_seed)

class PerturbedDepopulationModel(DepopulationModel):
    descr = "perturbed depopulation model"
    
    def __init__(self, numerical_depop_model):
        super(PerturbedDepopulationModel, self).__init__(numerical_depop_model.active_medium)
        
        self.numerical_depop_model = numerical_depop_model
        self.rel_perturbation = numerical_depop_model.rtol
    
    def _rate(self, inversion):
        rate = self.numerical_depop_model._rate(inversion)
        rate *= (1.0 - self.rel_perturbation)
        return rate

@concrete
class NullDepopulation(AnalyticalDepopulationModel):
    descr = "null depopulation"
    
    def _rate(self, inversion):
        return 0.0

@concrete
class FluorescenceModel(AnalyticalDepopulationModel):
    descr = "fluorescence"
    
    def _rate(self, inversion):
        rate = self._decay_rate(inversion, self.active_medium.doping_agent.transition_probability)
        return rate

@concrete
class LinfordASEModel(AnalyticalDepopulationModel):
    descr = "Linford ASE model"
    
    def __init__(self, *args, **kwargs):
        super(LinfordASEModel, self).__init__(*args, **kwargs)
        doping_agent = self.active_medium.doping_agent
        photon_energy = energy.photon_energy(doping_agent.lasing_wavelen)
        self._saturation_intensity = photon_energy * doping_agent.laser_transition_probability / doping_agent.xsection
    
    def _rate(self, inversion):
        # at zero inversion, depopulation rate is undefined but approaches zero
        if inversion == 0.0:
            return 0.0
        active_medium = self.active_medium
        doping_agent = active_medium.doping_agent
        active_solid_angle = active_medium.active_solid_angle
        ln_gain = doping_agent.xsection * inversion * active_medium.length
        gain = math.exp(ln_gain)
        intensity = self._saturation_intensity * (active_solid_angle/4.0) * (gain-1.0)**1.5/(gain*ln_gain)**0.5
        power = 2.0 * active_medium.aperture * intensity
        rate = energy.photon_count(doping_agent.lasing_wavelen, power)
        # fluorescence in directions outside of the two active solid angles:
        rate += (1.0 - 2.0 * active_solid_angle / (4.0 * math.pi)) * self._decay_rate(inversion, doping_agent.laser_transition_probability)
        # fluorescence to other states:
        rate += self._decay_rate(inversion, doping_agent.nonlaser_transition_probability)
        return rate

@concrete
class SchulzASEModel(AnalyticalDepopulationModel):
    descr = "Schulz ASE model"
    
    def _rate(self, inversion):
        active_medium = self.active_medium
        doping_agent = active_medium.doping_agent
        xsection = doping_agent.xsection
        ln_gain = xsection * inversion * active_medium.length
        gain = math.exp(ln_gain)
        gamma_prime = ln_gain * doping_agent.transition_probability + (active_medium.active_solid_angle * (gain - ln_gain - 1.0) / (2.0 * math.pi)) * doping_agent.laser_transition_probability
        rate = gamma_prime * active_medium.aperture / xsection
        return rate

@concrete
class RossApproximateASEModel(AnalyticalDepopulationModel):
    descr = "Ross ASE model (approximate)"
    
    _Nd_YLF_branching_ratio = 0.56 * 0.525 / 0.500
    
    def _B(self, inversion):
        active_medium = self.active_medium
        length = active_medium.length
        ln_gain = active_medium.doping_agent.xsection * inversion * length
        gain = math.exp(ln_gain)
        B = 0.05 * (2.0*active_medium.radius/length)**0.3 * gain
        B *= active_medium.doping_agent.branching_ratio / self._Nd_YLF_branching_ratio
        return B
    
    def _rate(self, inversion):
        active_medium = self.active_medium
        B = self._B(inversion)
        rate_coef = 1.0 + B
        rate = inversion * rate_coef * active_medium.doping_agent.transition_probability
        rate *= active_medium.volume
        return rate

@concrete
class RossNumericalASEModel(NumericalDepopulationModel):
    descr = "Ross ASE model (numerical)"
    
    _integrate_B = lambda self, *args: native._ross_ase_model_integrate_B(self, *args)
    
    def __init__(self, active_medium, rtol, min_samples, prng_seed=-1, sample_count_multiplier=16):
        super(RossNumericalASEModel, self).__init__(active_medium, rtol, min_samples, prng_seed)
        
        self.sample_count_multiplier = sample_count_multiplier
    
    def _rate_coef(self, inversion):
        nsamples = self.min_samples
        while True:
            B, B_rel_error = self._integrate_B(inversion, nsamples)
            B_abs_error = B_rel_error * B
            rate_coef = 1.0 + B
            
            rate_rel_error = B_abs_error / rate_coef
            if rate_rel_error < self.rtol:
                break
            
            nsamples *= self.sample_count_multiplier
            if nsamples > self.max_samples:
                util.warn("max. depop. rate sample size (%s) exceeded; latest rel. error: %s" % (self.max_samples, rate_rel_error), stacklevel=2)
                break
        
        return rate_coef
    
    def _rate(self, inversion):
        active_medium = self.active_medium
        rate_coef = self._rate_coef(inversion)
        rate = inversion * rate_coef * active_medium.doping_agent.transition_probability
        rate *= active_medium.volume
        return rate
    
    def rate_rel_stddev(self, inversion):
        nsamples = nsubsamples = int(math.ceil(math.sqrt(self.max_samples)))
        B, _ = self._integrate_B(inversion, self.max_samples)
        B_stddev = native._ross_ase_model_B_stddev(self, inversion, nsamples, nsubsamples)
        rate_coef = 1.0 + B
        rate_rel_stddev = B_stddev / rate_coef
        return rate_rel_stddev

@concrete
class RossHybridASEModel(RossNumericalASEModel):
    descr = "Ross ASE model (hybrid)"
    
    _integrate_B = lambda self, *args: native._ross_ase_model_integrate_BP(self, 0.0, 0.0, self.z_rel * self.active_medium.length, *args)
    
    def __init__(self, active_medium, rtol, min_samples, prng_seed=-1, sample_count_multiplier=16, z_rel=0.0):
        assert 0.0 <= z_rel <= 1.0, "z_rel is not in the interval [0.0, 1.0]"
        
        super(RossHybridASEModel, self).__init__(active_medium, rtol, min_samples, prng_seed, sample_count_multiplier)
        
        self.z_rel = z_rel
