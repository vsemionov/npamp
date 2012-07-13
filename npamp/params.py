
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


import model


version = 1.0


# physical parameters:

lasing_wavelen = 1064.0e-9

dopant_xsection = 2.8e-23
dopant_upper_lifetime = 230.0e-6
dopant_lower_lifetime = (115.0e-12 + 225.0e-12) / 2.0
dopant_branching_ratio = 0.56
dopant_concentration = 1.38e26

medium_radius = 2.0e-3
medium_length = 100.0e-3
medium_refr_idx = 1.82

initial_inversion = 0.0

pump_wavelen = 808.0e-9
pump_duration = 250.0e-6
pump_power = 4000.0
pump_efficiency = 0.5

depop_model_class = model.depop.RossApproximateASEModel

beam_class = model.beam.GaussianBeam
beam_radius = 2.0e-3

pulse_class = model.pulse.GaussianPulse
pulse_energy = 1.0e-6
pulse_duration = 6.0e-12

train_pulse_count = 1
train_pulse_period = 1.0e-9


# numerical parameters:

depop_rate_rtol = 5.0e-2
depop_rate_min_samples = 4096
depop_model_extra_args = {}

inverter_class = model.inverter.RungeKuttaInverter # RungeKuttaInverter usually makes fewer derivative evaluations for RossNumericalASEModel; for other depopulation models EulerInverter seems faster
inversion_rtol = 1.0e-2
inversion_min_count_t = 0
inversion_validate = False

ext_depop_models = [model.depop.FluorescenceModel, model.depop.RossApproximateASEModel, model.depop.RossNumericalASEModel]
ext_alt_depop_model = model.depop.FluorescenceModel
ext_opt_inversion_rdiff_max = 0.1
ext_opt_fluence_max = 1.0e4
ext_opt_pump_duration = (dopant_upper_lifetime / 2.0, dopant_upper_lifetime * 2.0)
ext_opt_pump_power = (2.0e3, 8.0e3)
ext_opt_pump_resolution = (24, 24)
ext_opt_geom_mediumradius = (medium_radius / 2.0, medium_radius * 2.0)
ext_opt_geom_beamradius = (beam_radius / 2.0, beam_radius * 2.0)
ext_opt_geom_resolution = (24, 24)

time_trunc_rtol = 1.0e-2
int_rtol = 1.0e-2
amp_rtol = 2.0e-2

min_count_rho = 0
min_count_phi = 0
min_count_z = 0
min_count_t = 0

integrator_classes = [model.integrator.SimpsonIntegrator]
amplifier_classes = [model.amplifier.HybridAmplifier, model.amplifier.NSFDAmplifier] # HybridAmplifier is faster and may produce less error in some cases; NSFDAmplifier is unconditionally stable and produces qualitatively correct solutions (with all expected properties)


# output parameters:

amplification = True

extended_mode = False
graphs = False
output_data = False

verbose = False

num_tasks = 0

lower_process_priority = True

pulse_num_stride = 100
extended_status_strides = (4, 0)

out_num_markers = 5
out_count_rho = 33
out_count_phi = 33
out_count_z = 25
out_count_t = 33

out_num_auto_contours = 5

output_rel_time = False
output_color = True
output_styled_lines = False
output_plot_titles = True

output_file_format = "png"
