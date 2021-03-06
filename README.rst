=====
NPAmp
=====
------------------------------------------------------------------------------------------------------------------------------------
Numerical model of a laser amplifier for pulse trains, accounting for amplified spontaneous emission and finite lower state lifetime
------------------------------------------------------------------------------------------------------------------------------------


NPAmp is an implementation of a numerical model of a laser amplifier. It is designed for practical and educational purposes.


.. image:: https://travis-ci.org/vsemionov/npamp.svg?branch=master
    :target: https://travis-ci.org/vsemionov/npamp


Screenshot
==========
.. image:: https://bitbucket.org/vsemionov/npamp/wiki/images/screenshots/basic.png


Examples of Results
===================
The following are examples of results from research, conducted with NPAmp.

Laser Amplifier Optimization
----------------------------
.. image:: https://bitbucket.org/vsemionov/npamp/wiki/images/examples/optimization/optimization/geometry/energy_out.png

Effects of Lower State Lifetime
-------------------------------
.. image:: https://bitbucket.org/vsemionov/npamp/wiki/images/examples/lower_state_lifetime/ref_pulse/lsl_effects.png

Effects of ASE
--------------
.. image:: https://bitbucket.org/vsemionov/npamp/wiki/images/examples/ase/pumping/inversions_evo.png

Effects of Gain Saturation
--------------------------
.. image:: https://bitbucket.org/vsemionov/npamp/wiki/images/examples/gain_saturation/ref_pulse/densities.png

Pulse Train Amplification
-------------------------
.. image:: https://bitbucket.org/vsemionov/npamp/wiki/images/examples/pulse_train/pulse_energy_gain.png


Overview
========

Mathematical Model Features
---------------------------
* Describes the propagation of single pulses or pulse trains in a laser amplifier.
* Accounts for the amplified spontaneous emission during pumping (by means of various ASE models).
* Accounts for the finite lower state lifetime during amplification.

Numerical Model Features
------------------------
* Produces physically consistent and, where possible, numerically stable results.
* Automatically determines discretization parameters from specified numerical error tolerances.
* Contains alternative numerical methods for the various intermediate computations.
* Automatically selects the most efficient combination of numerical methods.

Implementation Features
-----------------------
* Computes a wide variety of results.
* Generates publication-quality figures.
* Performs computations in parallel, on multiple CPUs, where possible.
* Operates on Windows and various Unix-like systems (Linux, BSD, etc.).
* Flexible and configurable operation and output.
* Simple graphical user interface.
* Command-line interface, facilitating research automation.
* Extensible by means of Python modules.
* Freely distributable under the terms of the BSD license.

Characteristics of the Modeled Amplifier and Input Signal
---------------------------------------------------------
* The amplifier is single-pass.
* The active medium is solid-state, four-level, with cylindrical geometry.
* The pumping is pulsed, homogeneous, and its optical power is a rectangular function of time.
* The active medium has a diffuse side surface and is water-cooled.
* There are no significant temperature gradients and mechanical stresses inside the active medium, therefore it is homogeneous.
* The active medium's temperature, and therefore its refractive index, is constant in time.
* The laser beam, propagating in the amplifier, is collimated.
* The input signal enters the amplifier immediately after the pumping; its power is null beforehand.
* The input signal consists of either a single pulse or a finite-length train of equally-separated in time pulses with arbitrary, but uniform, spatial and temporal shapes.
* The signal's frequency spectrum width is much smaller than its carrier frequency.
* The signal's total duration is much shorter than the upper state lifetime.
* The signal's optical intensity is proportional to its power and fluence.
* The time between amplification and the subsequent pumping is much longer than the upper state lifetime.
* During pumping, population inversion does not change significantly in the transit time of the amplifier.
* The pump band has a lifetime, much shorter than the upper state lifetime, and therefore its population is negligible; the ground state population is much larger than the sum of the laser state populations, and can therefore be considered constant in time; nonradiative transitions, originating from the upper state, can be neglected; the frequency dependence of the stimulated emission cross section of the laser transition can be neglected.
* During amplification, spontaneous emission is negligible; the transit time of the amplifier is much shorter than the upper state lifetime; absorption and scattering losses of signal radiation are negligible.


Installation
============

Windows
-------
The recommended way of installing NPAmp on Windows is to download and run the provided installer.

Unix
----
On Unix, NPAmp needs to be built from source. First, the required build and runtime dependencies need to be present and accessible on the system. The required packages, and their tested versions/variants, are:

* python (2.7.3)
* cython (0.15.1)
* numpy (1.6.1)
* scipy (0.9.0)
* matplotlib (1.1.1)
* pyside (1.1.0)
* make (3.81)
* C compiler (gcc 4.6.3)

Next, download the source archive, extract it into the target directory, and execute ``make`` inside that directory. The executable files are <target>/npamp/npamp.py and <target>/npamp/gnpamp.py, where <target> denotes the path to the target directory.


Operation
=========
NPAmp has a graphical interface (GNPAmp) and a command-line interface (NPAmp).

The model accepts a set of physical, numerical, operational, and output parameters. All of these parameters can be adjusted from the graphical interface. The parameters can optionally be written to and loaded from configuration files. 

There are two different modes of operation -- basic mode and extended mode. In basic mode, the results are computed only for fixed values of the input parameters. In extended mode, the results are computed for varying values of some parameters and/or for different models of some processes. This facilitates the studying of the effects of the various physical parameters on the results, as well as the effects of the different factors, taken into account by the different models.

Execution
---------
The most important parameters, affecting the execution of the model, are:

* extended mode: enables extended mode
* initial inversion: specifies the initial population inversion for the amplification stage; if non-zero, pumping-stage computations will not be performed
* perform amplification: specifies whether amplification-stage computations will be performed
* generate graphs: enables the generation of result graphs

Basic Mode
----------
The most important results, currently computed in basic mode, are:

* pumping stage:

 + numerical results:

  - pump energy
  - final population inversion
  - final small signal gain 
  - final stored energy

 + graphical results:

  - temporal evolution of population inversion

* amplification stage:

 + numerical results:

  - input energy
  - output energy
  - energy gain
  - extraction efficiency
  - optical to optical efficiency
  - maximum output fluence
  - relative gain decrease (the difference between the gain of the first and last pulses, divided by the gain of the first pulse)

 + graphical results:

  - energy gain of individual pulses in the train
  - longitudinal and temporal evolution of upper state population, lower state population, population inversion, and photon density
  - temporal shape of the input and output photon density
  - longitudinal evolution of the signal's fluence
  - radial shape of the input and output fluence
  - transverse shape of the input and output fluence

Extended Mode
-------------
The following computations are currently performed in extended mode:

* effects of lower state lifetime on the temporal shape and fluence of the output pulse

 + computed results:

  - output fluence for zero, finite, and infinite lifetimes
  - temporal shape of the output photon density for the above lifetimes

* comparison of various ASE models:

 + most important computed results:

  - population inversion and depopulation rate immediately after pumping, for each model
  - temporal evolution of population inversion and depopulation rate during pumping, for all models
  - dependence of depopulation rate on population inversion, for all models
  - temporal evolution of the relative difference between the population inversion, yielded by the compared models and by a reference model (example use: evaluation of the effect of ASE on population inversion by comparing results from ASE and fluorescence models)

* amplifier optimization:

 + optimized sets of parameters:

  - pumping system -- pump power and pump duration
  - geometry -- active medium diameter and input beam diameter

 + maximized output characteristics:

  - output energy
  - energy gain
  - extraction efficiency
  - optical to optical efficiency

 + minimized output characteristics:

  - relative gain decrease

 + imposed restrictions:

  - limited effect of ASE on population inversion during pumping
  - limited maximum output fluence (under the damage threshold)

 + most important computed results:

  - values of the optimized parameters, for which the output characteristics are extremal, while all restrictions are met
  - constraints on the domain of the optimized parameters
  - dependence of the following quantities on the optimized parameters: pump energy, population inversion, small signal gain, stored energy, effect of ASE on population inversion, input energy, output energy, energy gain, extraction efficiency, optical to optical efficiency, maximum output fluence, relative gain decrease


Practical Use
=============

ASE Models
----------
Of all implemented depopulation models, the Ross ASE model accounts for the highest number of contributing factors, and is expected to yield the most accurate results. In the graphical interface, this model is labeled as "RossApproximateASEModel", and is selected by default.

Numerical Methods
-----------------
The method for computing the depopulation rate depends on the selected ASE model. The numerical method for computing pumping-stage population inversion is selected manually.

Amplification-stage computations utilize a combination of methods (an integration method and an "amplification" method). More than one method can be selected for each of the two types. In that case, NPAmp will probe all possible combinations of the selected methods of the two types and will use the most efficient one (requiring the least total number of discretization steps, depending on physical and numerical parameters).

Generally, it is advisable to select both implemented amplification methods, since each one has advantages and disadvantages over the other one, and is more efficient in different cases. HybridAmplifier is first-order accurate in time and second-order accurate in space, but the physical consistency of its results is conditional and depends on the physical parameters and initial conditions. NSFDAmplifier is first-order accurate, but its results are unconditionally physically consistent. Both methods have a computational complexity of O(MxN), but HybridAmplifier generally executes faster for the same discretization parameters.

Numerical Parameters
--------------------
In most cases, the default values of the numerical parameters (in the Numerics tab of the graphical interface) are an acceptable compromise between accuracy and efficiency. In the absence of specific requirements, these parameters do not need to be modified.

In cases when examination of the temporal or spatial evolution of the computed quantities is required, the minimum numbers of discretization steps could need adjustment in order to increase the smoothness of the generated graphs.


Extensions
==========
NPAmp is extensible by means of Python modules. To install an extension:

1. Create the <home>/.npamp/extensions directory (if it doesn't exist), where <home> denotes the path to the user's home directory.
2. Copy the extension module to the above directory.

The extension will be loaded the next time NPAmp is run. To view the list of currently installed extensions from the graphical interface, use "About NPAmp" from the Help menu.

For an example of an extension and a configuration which uses it, see examples/sawtooth_pulse.py and examples/extensions.npc, respectively, in the installation directory.


Common Problems
===============
* Warning messages in the output. NPAmp could issue warnings about exceeding discretization limits while trying to ensure the required accuracy. This could occur during pumping-stage computations (depopulation rate and population inversion) or during the process of probing different combinations of amplification-stage methods. In the former case, select a different inversion method or reduce the corresponding error tolerance (depopulation rate or inversion). In the latter case, the warnings could be unrelated to the method combination that will be used, and therefore computations could still be performed with the required accuracy. To verify this, enable verbose output to find out which method combination causes the warnings and which one is actually used, or select different (sets of) integration and amplification methods. In cases of extreme values of the physical or numerical parameters, all method combinations could require unreasonably large numbers of discretization steps. In that case, the only solution is to reduce the corresponding error tolerance (amplification or integration).
* Unordered task progress indication numbers in the output. When computations are performed in parallel, the numbers, indicating the progress of the current task, could be unordered. This is normal and does not affect the results.


Contacts and Support
====================
Developed and maintained by Victor Semionov <vsemionov@gmail.com>.

Studied under the guidance of Assoc. Prof. Ivan Buchvarov and Assistant Prof. Alexander Gaydardzhiev at the Faculty of Physics, Sofia University.

For the latest information and downloads, visit the main website:
http://www.vsemionov.org/npamp/

For the source code and other technical resources, visit the development website:
https://bitbucket.org/vsemionov/npamp

Website of the Non-linear Optics and Solid State Lasers Laboratory:
http://www.phys.uni-sofia.bg/~ibuch/

For feedback, problem reports, or any other queries, email:
Victor Semionov <vsemionov@gmail.com>
