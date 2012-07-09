
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

import types

import params
import output
import core


class ConfigurationError(Exception):
    pass


copy_conf = lambda conf: {k: v for k, v in conf.items() if not k.startswith('_') and type(v) is not types.ModuleType}


def load_conf(defaults, path):
    conf = defaults.copy()
    diff = dict()
    if path:
        execfile(path, diff)
    conf.update(diff)
    if conf.pop("version", params.version) != params.version:
        raise ConfigurationError("unsupported configuration file format version")
    conf = copy_conf(conf)
    return conf


def validate():
    def validate_param_bounds(name, value):
        vtype = type(value)
        if vtype in (int, float):
            if vtype is float and name not in param_zero_allowed:
                if value <= 0.0:
                    raise ConfigurationError("parameter \"%s\" has (or contains) a zero value" % name)
            elif value < 0.0:
                raise ConfigurationError("parameter \"%s\" has (or contains) a negative value" % name)
            if vtype is float and name not in param_infty_allowed:
                if math.isinf(value):
                    raise ConfigurationError("parameter \"%s\" has (or contains) an infinite value" % name)
            if name in param_bounds:
                min_value, max_value = param_bounds[name]
                if min_value is not None:
                    if value < min_value:
                        raise ConfigurationError("parameter \"%s\" has (or contains) a value less than %s" % (name, min_value))
                if max_value is not None:
                    if value > max_value:
                        raise ConfigurationError("parameter \"%s\" has (or contains) a value greater than %s" % (name, max_value))
        elif vtype in (tuple, list):
            for element in value:
                validate_param_bounds(name, element)
    
    def validate_interval(name, values):
        if values[0] >= values[1]:
            raise ConfigurationError("parameter \"%s\" has an upper bound not greater than its lower bound" % name)
    
    def validate_list_uniques(name, values):
        if len(values) != len(set(values)):
            raise ConfigurationError("parameter \"%s\" contains duplicate values" % name)
    
    param_zero_allowed = set([
        "dopant_branching_ratio",
        "dopant_lower_lifetime",
        "initial_inversion",
    ])
    param_infty_allowed = set([
        "dopant_upper_lifetime",
        "dopant_lower_lifetime",
        "ext_opt_inversion_rdiff_max",
        "ext_opt_fluence_max",
    ])
    param_bounds = {
        "train_pulse_count": (1, None),
        "depop_rate_min_samples": (16, None),
        "out_markers_step_divisor": (1, None),
        "out_rho_steps_divisor": (1, None),
        "out_phi_steps_divisor": (1, None),
        "out_z_steps_divisor": (1, None),
        "out_t_steps_divisor": (1, None),
        "ext_opt_pump_resolution": (2, None),
        "ext_opt_geom_resolution": (2, None),
        "dopant_branching_ratio": (None, 1.0),
        "pump_efficiency": (None, 1.0),
    }
    param_intervals = set([
        "ext_opt_pump_duration",
        "ext_opt_pump_power",
        "ext_opt_geom_mediumradius",
        "ext_opt_geom_beamradius",
    ])
    
    conf = copy_conf(params.__dict__)
    
    for parameter, value in conf.items():
        validate_param_bounds(parameter, value)
    
    for parameter in param_intervals:
        values = conf[parameter]
        validate_interval(parameter, values)
    
    active_medium = core.create_medium(None)
    input_beam = core.create_beam()
    ref_pulse = core.create_pulse(active_medium, input_beam, input_beam.rho_ref, input_beam.phi_ref)
    
    if not (params.pump_wavelen < params.lasing_wavelen or params.initial_inversion):
        output.warn("pump wavelength is not less than lasing wavelength")
    
    if not (params.dopant_lower_lifetime <= params.dopant_upper_lifetime / 10.0 or params.initial_inversion):
        output.warn("approximation validity condition violated: lower state lifetime is not much shorter than upper state (fluorescence) lifetime")
    
    train_duration = params.pulse_duration + (params.train_pulse_count - 1) * params.train_pulse_period
    if not (train_duration <= params.dopant_upper_lifetime / 10.0 or not params.amplification):
        output.warn("approximation validity condition violated: signal duration is not much shorter than upper state (fluorescence) lifetime")
    
    if not (params.train_pulse_period >= ref_pulse.duration or params.train_pulse_count == 1 or not params.amplification):
        raise ConfigurationError("invalid parameters: pulse repetition period is less than (extended) pulse duration")
