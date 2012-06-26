
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


import pamp


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

loss_model_class = pamp.loss.RossApproximateASEModel

beam_class = pamp.beam.GaussianBeam
beam_radius = 2.0e-3

pulse_class = pamp.pulse.GaussianPulse
pulse_energy = 1.0e-6
pulse_duration = 6.0e-12

train_pulse_count = 1
train_pulse_period = 1.0e-9


# numerical parameters:

loss_rtol = 5.0e-2
loss_min_count = 4096
loss_model_extra_args = dict()

inverter_class = pamp.inverter.RungeKuttaInverter # RungeKuttaInverter usually makes fewer derivative evaluations for RossNumericalASEModel; for other loss models EulerInverter seems faster
initial_inversion_rtol = 1.0e-2
initial_inversion_validate = False
initial_inversion_min_count_t = 0

loss_dependence_alt_model = pamp.loss.FluorescenceLossModel
loss_dependence_inversion_rdiff_max = 0.1
fluence_dependence_fluence_max = 1.0e4
pumpdep_duration_interval = (dopant_upper_lifetime / 2.0, dopant_upper_lifetime * 2.0)
pumpdep_power_interval = (2.0e3, 8.0e3)
pumpdep_step_counts = (24, 24)
geomdep_mediumradius_interval = (medium_radius / 2.0, medium_radius * 2.0)
geomdep_beamradius_interval = (beam_radius / 2.0, beam_radius * 2.0)
geomdep_step_counts = (24, 24)

pulse_rel_energy_trunc = 1.0e-3

energy_rtol = 1.0e-6
fluence_rtol = 1.0e-6
min_count_rho = 0
min_count_phi = 0

amp_rtol = 1.0e-2
min_count_z = 0
min_count_t = 0

integrator_classes = (pamp.integral.TrapezoidIntegrator, pamp.integral.SimpsonIntegrator, pamp.integral.RombergIntegrator, )
amplifier_classes = (pamp.amplifier.HybridAmplifier, pamp.amplifier.NSFDAmplifier, ) # HybridAmplifier is faster and may produce less error in some cases; NSFDAmplifier is unconditionally stable and produces qualitatively correct solutions (with all expected properties)
compared_loss_models = (pamp.loss.FluorescenceLossModel, pamp.loss.RossApproximateASEModel, pamp.loss.RossNumericalASEModel, )


# output parameters:

amplification = True

extended_mode = False
graphs = False
output_data = False

verbose = False

num_tasks = 0

lower_process_priority = True

pulse_num_stride = 100
extended_status_strides = (0, 0)

out_markers_step_divisor = 5
out_rho_steps_divisor = 33
out_phi_steps_divisor = 65
out_z_steps_divisor = 33
out_t_steps_divisor = 33

num_contours = 5

output_rel_time = False
output_color = True
output_styled_lines = False
output_graph_titles = True

output_file_format = "png"
