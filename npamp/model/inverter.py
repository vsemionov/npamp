
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


import numpy as np

import discrete
import util
import exc


class PopulationInverter(object):
    
    evals_per_step = None
    
    def __init__(self, active_medium, pump_system, depop_model):
        self.active_medium = active_medium
        self.pump_system = pump_system
        self.depop_model = depop_model
        
        self.T = None
        self.inversion = None
    
    def _deriv(self, inversion):
        active_medium = self.active_medium
        pump_system = self.pump_system
        volume = active_medium.volume
        pump_rate_eff = pump_system.effective_pump_rate
        depop_rate_total = self.depop_model.rate(inversion)
        deriv = (pump_rate_eff - depop_rate_total) / volume
        return deriv
    
    @staticmethod
    def _step(f, dt, n):
        raise NotImplementedError()
    
    def _integrate(self, rtol, count_t):
        # Adaptive step size is not necessary.
        # It would be useful only in cases when the majority of the time is spent in the steady state.
        # However, that is not a practical case -- any pumping, applied after reaching the steady state, is wasted.
        def is_stable(n, dn): # oscillation detection
            mid = n + dn/2.0
            return dn >= -2.0*rtol * mid
        
        f = self._deriv
        
        T, dt = np.linspace(0.0, self.pump_system.duration, count_t, retstep=True)
        inversion = np.empty(count_t)
        inversion[0] = 0.0
        
        for k in range(count_t-1):
            n = inversion[k]
            try:
                dn = self._step(f, dt, n)
            except (ValueError, ArithmeticError):
                return None
            if not is_stable(n, dn):
                return None
            inversion[k+1] = n + dn
        
        self.T = T
        self.inversion = inversion
        return inversion[-1]
    
    def invert(self, rtol, min_count_t):
        max_divs = 15
        
        min_count_t = max(min_count_t, 2)
        divs = discrete.divs(min_count_t)
        if divs > max_divs:
            raise exc.NumericalError("min. inversion time step count (%s) and corresponding min. divs (%s) too large; max. divs: %s" % (min_count_t, divs, max_divs))
        
        divs = max(divs-1, 0)
        last_res = None
        while True:
            count_t = discrete.steps(divs)
            res = self._integrate(rtol, count_t)
            
            diff = None
            if None not in (last_res, res):
                diff = abs(res - last_res)
                if diff <= rtol * abs(res):
                    break
            
            divs += 1
            if divs > max_divs:
                if diff is not None:
                    util.warn("max. inversion time divs (%s) exceeded; rtol: %s; latest step count: %s; latest difference: %s (current: %s; last: %s)" % (max_divs, rtol, count_t, diff, res, last_res), stacklevel=2)
                    break
                else:
                    msg = "max. inversion time divs (%s) exceeded; rtol: %s; latest step count: %s" % (max_divs, rtol, count_t)
                    if res is not None:
                        util.warn(msg, stacklevel=2)
                        break
                    else:
                        raise exc.NumericalError(msg)
            last_res = res
        
        return res

class EulerInverter(PopulationInverter):
    
    evals_per_step = 1
    
    @staticmethod
    def _step(f, dt, n):
        dn = dt * f(n)
        return dn

class RungeKuttaInverter(PopulationInverter):
    
    evals_per_step = 4
    
    @staticmethod
    def _step(f, dt, n):
        k1 = dt * f(n)
        k2 = dt * f(n + k1/2.0)
        k3 = dt * f(n + k2/2.0)
        k4 = dt * f(n + k3)
        dn = (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0
        return dn
