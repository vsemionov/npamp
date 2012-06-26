
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

import numpy as np

cimport cython
cimport numpy as np
cimport libc.math

cimport randomkit


DTYPE = np.double
ctypedef np.double_t DTYPE_t


@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def _hybrid_amplifier_calc_output(self):
    active_medium = self.active_medium
    input_pulse = self.input_pulse
    
    cdef DTYPE_t rho = self.rho
    cdef DTYPE_t phi = self.phi
    
    cdef np.ndarray[DTYPE_t, ndim=1] Z = self.Z
    cdef np.ndarray[DTYPE_t, ndim=1] T = self.T
    
    cdef np.ndarray[DTYPE_t, ndim=2] density = self.density
    cdef np.ndarray[DTYPE_t, ndim=2] upper = self.population[0]
    cdef np.ndarray[DTYPE_t, ndim=2] lower = self.population[1]
    
    cdef int have_initial_population = self.initial_population is not None
    cdef np.ndarray[DTYPE_t, ndim=1] initial_upper = self.initial_population[0] if have_initial_population else None
    cdef np.ndarray[DTYPE_t, ndim=1] initial_lower = self.initial_population[1] if have_initial_population else None
    cdef np.ndarray[DTYPE_t, ndim=1] initial_density = self.initial_density
    
    cdef unsigned int count_z = self.count_z
    cdef unsigned int count_t = self.count_t
    
    cdef DTYPE_t dz = self.dz
    cdef DTYPE_t dt = self.dt
    
    cdef DTYPE_t light_speed = active_medium.light_speed
    cdef DTYPE_t xsection = active_medium.doping_agent.xsection
    cdef DTYPE_t lower_lifetime = active_medium.doping_agent.lower_lifetime
    
    cdef DTYPE_t _sigma_c_dt = xsection * light_speed * dt
    cdef DTYPE_t _dt_over_tau = dt / lower_lifetime
    
    cdef DTYPE_t _n1, _n2, _n, _nn, _phi, _sigma_c_n_phi_dt
    
    cdef DTYPE_t density_coef = xsection * dz / 2.0
    
    cdef unsigned int l, k
    with nogil:
        for l in range(count_z):
            for k in range(count_t):
                # population
                if k > 0:
                    _phi = density[l, k-1]
                    _n1 = lower[l, k-1]
                    _n2 = upper[l, k-1]
                    _n = _n2 - _n1
                    _sigma_c_n_phi_dt = _sigma_c_dt * _n * _phi
                    
                    upper[l, k] = _n2 - _sigma_c_n_phi_dt
                    
                    if lower_lifetime != 0.0:
                        lower[l, k] = _n1 + _sigma_c_n_phi_dt - _n1 * _dt_over_tau
                    else:
                        lower[l, k] = 0.0
                else:
                    if have_initial_population:
                        upper[l, 0] = initial_upper[l]
                        lower[l, 0] = initial_lower[l]
                    else:
                        with gil:
                            upper[l, 0] = active_medium.initial_inversion.inversion(rho, phi, Z[l])
                        lower[l, 0] = 0.0
                
                # density
                if l > 0:
                    _n = upper[l-1, k] - lower[l-1, k]
                    _nn = upper[l, k] - lower[l, k]
                    _phi = density[l-1, k]
                    density[l, k] = _phi * (1.0 + density_coef * _n) / (1.0 - density_coef * _nn)
                else:
                    if initial_density is not None:
                        density[0, k] = initial_density[k]
                    else:
                        with gil:
                            density[0, k] = input_pulse.density(T[k])

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def _nsfd_amplifier_calc_output(self):
    # for the analytical solution of the rate equations, see http://eqworld.ipmnet.ru/en/solutions/sysode/sode0101.pdf
    active_medium = self.active_medium
    input_pulse = self.input_pulse
    
    cdef DTYPE_t rho = self.rho
    cdef DTYPE_t phi = self.phi
    
    cdef np.ndarray[DTYPE_t, ndim=1] Z = self.Z
    cdef np.ndarray[DTYPE_t, ndim=1] T = self.T
    
    cdef np.ndarray[DTYPE_t, ndim=2] density = self.density
    cdef np.ndarray[DTYPE_t, ndim=2] upper = self.population[0]
    cdef np.ndarray[DTYPE_t, ndim=2] lower = self.population[1]
    
    cdef int have_initial_population = self.initial_population is not None
    cdef np.ndarray[DTYPE_t, ndim=1] initial_upper = self.initial_population[0] if have_initial_population else None
    cdef np.ndarray[DTYPE_t, ndim=1] initial_lower = self.initial_population[1] if have_initial_population else None
    cdef np.ndarray[DTYPE_t, ndim=1] initial_density = self.initial_density
    
    cdef unsigned int count_z = self.count_z
    cdef unsigned int count_t = self.count_t
    
    cdef DTYPE_t dz = self.dz
    cdef DTYPE_t dt = self.dt
    
    cdef DTYPE_t light_speed = active_medium.light_speed
    cdef DTYPE_t xsection = active_medium.doping_agent.xsection
    cdef DTYPE_t lower_lifetime = active_medium.doping_agent.lower_lifetime
    
    cdef DTYPE_t _sigma_c = xsection * light_speed
    cdef DTYPE_t _tau_inv = 1.0 / lower_lifetime
    cdef DTYPE_t _neg_tau_inv_dt = _tau_inv * dt
    
    cdef DTYPE_t _n1, _n2, _n, _phi, _sigma_c_phi, _2_sigma_c_phi, _n1_over_sigma_c_phi, _sum, _sqrt_D, _lambda1, _lambda2, _lambda1_dt, _lambda2_dt, _c0, _c1, _c2, _C1, _C2
    
    cdef DTYPE_t density_coef = xsection * dz
    
    cdef unsigned int l, k
    with nogil:
        for l in range(count_z):
            for k in range(count_t):
                # population
                if k > 0:
                    _n2 = upper[l, k-1]
                    _phi = density[l, k-1]
                    _sigma_c_phi = _sigma_c * _phi
                    if lower_lifetime != 0.0:
                        _n1 = lower[l, k-1]
                        if _phi != 0.0:
                            _2_sigma_c_phi = 2.0 * _sigma_c_phi
                            _n1_over_sigma_c_phi = _n1 / _sigma_c_phi
                            _sum = - (_tau_inv + _2_sigma_c_phi)
                            _sqrt_D = libc.math.sqrt(_tau_inv*_tau_inv + _2_sigma_c_phi*_2_sigma_c_phi)
                            _lambda1 = (_sum + _sqrt_D) / 2.0
                            _lambda2 = (_sum - _sqrt_D) / 2.0
                            _lambda1_dt = _lambda1 * dt
                            _lambda2_dt = _lambda2 * dt
                            _c0 = _sigma_c_phi
                            _c1 = (_tau_inv + _sqrt_D) / 2.0
                            _c2 = (_tau_inv - _sqrt_D) / 2.0
                            _C1 = (_n2 - _n1_over_sigma_c_phi * _c2) / _sqrt_D
                            _C2 = _n1_over_sigma_c_phi - _C1
                            lower[l, k] = _C1 * _c0 * libc.math.exp(_lambda1_dt) + _C2 * _c0 * libc.math.exp(_lambda2_dt)
                            upper[l, k] = _C1 * _c1 * libc.math.exp(_lambda1_dt) + _C2 * _c2 * libc.math.exp(_lambda2_dt)
                        else:
                            lower[l, k] = _n1 * libc.math.exp(_neg_tau_inv_dt)
                            upper[l, k] = _n2
                    else:
                        lower[l, k] = 0.0
                        upper[l, k] = _n2 * libc.math.exp(- _sigma_c_phi * dt)
                else:
                    if have_initial_population:
                        lower[l, 0] = initial_lower[l]
                        upper[l, 0] = initial_upper[l]
                    else:
                        lower[l, 0] = 0.0
                        with gil:
                            upper[l, 0] = active_medium.initial_inversion.inversion(rho, phi, Z[l])
                
                # density
                if l > 0:
                    _n = upper[l-1, k] - lower[l-1, k]
                    _phi = density[l-1, k]
                    density[l, k] = _phi * libc.math.exp(density_coef * _n)
                else:
                    if initial_density is not None:
                        density[0, k] = initial_density[k]
                    else:
                        with gil:
                            density[0, k] = input_pulse.density(T[k])


cdef randomkit.rk_state _rk_state
randomkit.rk_randomseed(&_rk_state)

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def _ross_ase_model_integrate_B(self, inversion, nsamples):
    active_medium = self.active_medium
    
    branching_ratio = active_medium.doping_agent.branching_ratio
    
    volume = active_medium.volume
    integration_volume = volume**2
    
    cdef DTYPE_t radius = active_medium.radius
    cdef DTYPE_t length = active_medium.length
    cdef DTYPE_t spacing = active_medium.doping_agent.spacing
    cdef DTYPE_t gain_coef = inversion * active_medium.doping_agent.xsection
    
    F = np.empty(nsamples)
    cdef np.ndarray[DTYPE_t, ndim=1] Fbuff = F
    
    cdef DTYPE_t diameter = 2.0 * radius
    cdef DTYPE_t radius2 = radius * radius
    cdef DTYPE_t x1, x2, y1, y2, z1, z2
    cdef DTYPE_t dx, dy, dz, r2, r
    
    cdef unsigned int N = nsamples
    cdef unsigned int n = 0
    with nogil:
        while n < N:
            x1 = diameter * randomkit.rk_double(&_rk_state) - radius
            y1 = diameter * randomkit.rk_double(&_rk_state) - radius
            if x1*x1 + y1*y1 > radius2:
                continue
            x2 = diameter * randomkit.rk_double(&_rk_state) - radius
            y2 = diameter * randomkit.rk_double(&_rk_state) - radius
            if x2*x2 + y2*y2 > radius2:
                continue
            z1 = randomkit.rk_double(&_rk_state) * length
            z2 = randomkit.rk_double(&_rk_state) * length
            dx = x2 - x1
            dy = y2 - y1
            dz = z2 - z1
            r2 = dx*dx + dy*dy + dz*dz
            r = libc.math.sqrt(r2)
            if r < spacing:
                continue
            Fbuff[n] = libc.math.exp(gain_coef * r) / (4.0*libc.math.M_PI * r2)
            n += 1
    
    F.sort()
    total = np.sum(F)
    Fmean = total / nsamples
    integral = integration_volume * Fmean
    
    F2mean = np.sum(F**2) / nsamples
    Fmean2 = Fmean**2
    stddev2 = abs(F2mean - Fmean2)
    stddev = math.sqrt(stddev2)
    abs_error = integration_volume * stddev / math.sqrt(nsamples)
    rel_error = abs_error / integral
    
    B = branching_ratio * gain_coef * integral / volume
    
    return B, rel_error

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def _ross_ase_model_integrate_BP(self, x1, y1, z1, inversion, nsamples):
    active_medium = self.active_medium
    
    branching_ratio = active_medium.doping_agent.branching_ratio
    
    volume = active_medium.volume
    
    cdef DTYPE_t radius = active_medium.radius
    cdef DTYPE_t length = active_medium.length
    cdef DTYPE_t spacing = active_medium.doping_agent.spacing
    cdef DTYPE_t gain_coef = inversion * active_medium.doping_agent.xsection
    
    F = np.empty(nsamples)
    cdef np.ndarray[DTYPE_t, ndim=1] Fbuff = F
    
    cdef DTYPE_t _x1 = x1
    cdef DTYPE_t _y1 = y1
    cdef DTYPE_t _z1 = z1
    
    cdef DTYPE_t diameter = 2.0 * radius
    cdef DTYPE_t radius2 = radius * radius
    cdef DTYPE_t x2, y2, z2
    cdef DTYPE_t dx, dy, dz, r2, r
    
    cdef unsigned int N = nsamples
    cdef unsigned int n = 0
    with nogil:
        while n < N:
            x2 = diameter * randomkit.rk_double(&_rk_state) - radius
            y2 = diameter * randomkit.rk_double(&_rk_state) - radius
            if x2*x2 + y2*y2 > radius2:
                continue
            z2 = randomkit.rk_double(&_rk_state) * length
            dx = x2 - _x1
            dy = y2 - _y1
            dz = z2 - _z1
            r2 = dx*dx + dy*dy + dz*dz
            r = libc.math.sqrt(r2)
            if r < spacing:
                continue
            Fbuff[n] = libc.math.exp(gain_coef * r) / (4.0*libc.math.M_PI * r2)
            n += 1
    
    F.sort()
    total = np.sum(F)
    Fmean = total / nsamples
    integral = volume * Fmean
    
    F2mean = np.sum(F**2) / nsamples
    Fmean2 = Fmean**2
    stddev2 = abs(F2mean - Fmean2)
    stddev = math.sqrt(stddev2)
    abs_error = volume * stddev / math.sqrt(nsamples)
    rel_error = abs_error / integral
    
    B = branching_ratio * gain_coef * integral
    
    return B, rel_error

@cython.cdivision(True)
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def _ross_ase_model_B_stddev(self, inversion, nsamples, nsubsamples):
    active_medium = self.active_medium
    
    cdef DTYPE_t branching_ratio = active_medium.doping_agent.branching_ratio
    
    cdef DTYPE_t volume = active_medium.volume
    
    cdef DTYPE_t radius = active_medium.radius
    cdef DTYPE_t length = active_medium.length
    cdef DTYPE_t spacing = active_medium.doping_agent.spacing
    cdef DTYPE_t gain_coef = inversion * active_medium.doping_agent.xsection
    
    cdef DTYPE_t integral_coef = branching_ratio * gain_coef * volume / nsubsamples
    
    F1 = np.empty(nsamples)
    F2 = np.empty(nsubsamples)
    cdef np.ndarray[DTYPE_t, ndim=1] F1buff = F1
    cdef np.ndarray[DTYPE_t, ndim=1] F2buff = F2
    
    cdef DTYPE_t diameter = 2.0 * radius
    cdef DTYPE_t radius2 = radius * radius
    cdef DTYPE_t x1, x2, y1, y2, z1, z2
    cdef DTYPE_t dx, dy, dz, r2, r
    
    cdef unsigned int N1 = nsamples
    cdef unsigned int N2 = nsubsamples
    cdef unsigned int n1
    cdef unsigned int n2
    n1 = 0
    while n1 < N1:
        x1 = diameter * randomkit.rk_double(&_rk_state) - radius
        y1 = diameter * randomkit.rk_double(&_rk_state) - radius
        if x1*x1 + y1*y1 > radius2:
            continue
        z1 = randomkit.rk_double(&_rk_state) * length
        n2 = 0
        while n2 < N2:
            x2 = diameter * randomkit.rk_double(&_rk_state) - radius
            y2 = diameter * randomkit.rk_double(&_rk_state) - radius
            if x2*x2 + y2*y2 > radius2:
                continue
            z2 = randomkit.rk_double(&_rk_state) * length
            dx = x2 - x1
            dy = y2 - y1
            dz = z2 - z1
            r2 = dx*dx + dy*dy + dz*dz
            r = libc.math.sqrt(r2)
            if r < spacing:
                continue
            F2buff[n2] = libc.math.exp(gain_coef * r) / (4.0*libc.math.M_PI * r2)
            n2 += 1
        F2.sort()
        F1[n1] = np.sum(F2) * integral_coef
        n1 += 1
    
    F1.sort()
    F1mean = np.sum(F1) / nsamples
    squares = (F1 - F1mean)**2
    squares.sort()
    stddev2 = np.sum(squares) / (nsamples - 1)
    stddev = math.sqrt(stddev2)
    
    return stddev
