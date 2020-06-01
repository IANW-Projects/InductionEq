#!/usr/bin/env python3
# -*- coding: <utf8> -*-

#This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

#Use this script to specify the general conditions and the testcase you
#want to investigate. It calls functions to initializes all variables and fields relevant for
#the testcase and advance the numerical solution in time. Depending on the
#properties that shall be investigated additional functionalities can be
#enabled, e.g. divergence cleaning.
import InductionEq as induction_eq
from math import pi
import sys
#prepare_vars() initializes containers (maps) in which all variables (label, value) are
#stored. Variables are grouped together by purpose:
#I_Mesh: Contains all variables concerning the discrete mesh, e.g.
#stepsize, box size
#I_TI: Contains all variables concerning time integration, e.g. time
#integrator, timestep
#I_IEq: Contains all variables relevant for the induction equation, e.g.
#discretization form, switch for additional artificial dissipation
#I_DC: Contains all variables concerning divergence cleaning, e.g.
#discretization, error threshold
#I_Tech: Contains all variables of technical nature, e.g. the device,
#computation precision
#I_RunOps: Specifies the parameters of a computation, e.g. which variables
#will be saved, the testcase

#Variables in capital letters are program specific defines which are set by
#openCL compile settings. If the value is a string in captial letters
#beginning with USE it acts as a switch, e.g. to enable artificial dissipation
#or to switch between different discretizations. In case the key is written
#in capital letters the program define will be set to the corresponding
#value of the key.
if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")


induction_eq.preparse_vars()


#global I_Mesh I_TI I_IEq I_DC I_Tech I_RunOps I_Results

I_Mesh = induction_eq.I_Mesh
I_TI = induction_eq.I_TI
I_IEq = induction_eq.I_IEq
I_DC = induction_eq.I_DC
I_Tech = induction_eq.I_Tech
I_RunOps = induction_eq.I_RunOps
I_Results = induction_eq.I_Results

# Mesh related parameters. The minimum number of nodes per direction is
# constrained by the number of boundary nodes, e.g. the 4th order method
# has 2*4 boundary nodes which means the minimum number of nodes amounts to 8.

N = 40
I_Mesh['NODES_X'],I_Mesh['NODES_Y'],I_Mesh['NODES_Z'] = (N, N, 12)
I_Mesh['XMIN'], I_Mesh['XMAX'] = (-1.0, 1.0)
I_Mesh['YMIN'], I_Mesh['YMAX'] = (-1.0, 1.0)
I_Mesh['ZMIN'], I_Mesh['ZMAX'] = (-1.0, 1.0)

# Time integration related variables
I_TI['cfl'] = 0.95  #Define the Courant–Friedrichs–Lewy condition
I_TI['final_time'] = 2*pi
#Chose the time integrator. Below is a list of up to date available
#options:
# SSPRK33, SSPRK104,
# KennedyCarpenterLewis2R54C, CalvoFrancoRandez2R64,
# CarpenterKennedy2N54, ToulorgeDesmet2N84F
I_TI['time_integrator'] = 'SSPRK33' #'KennedyCarpenterLewis2R54C';


# Induction equation related variables
#Specify how the three part of the linear induction equation shall be
#computed.
I_IEq['form_uibj'] = 'USE_UIBJ_PRODUCT' # PRODUCT, SPLIT, CENTRAL
I_IEq['form_source'] = 'USE_SOURCE_CENTRAL' # CENTRAL, SPLIT, ZERO
I_IEq['form_ujbi'] = 'USE_UJBI_CENTRAL' # SPLIT, CENTRAL, PRODUCT
#Enable or disable Hall-Term
I_IEq['hall_term'] = 'NONE' # NONE, USE_HALL
#Enable or disable artificial dissipation
I_IEq['dissipation'] = 'NONE' # NONE, USE_ARTIFICIAL_DISSIPATION
#Specify what kind of artifical dissipation shall be used
    #USE_ADAPTIVE_DISSIPATION, USE_FIRST_ORDER_DISSIPATION, USE_HIGH_ORDER_DISSIPATION
I_IEq['dissipation_form'] = 'USE_HIGH_ORDER_DISSIPATION'
I_IEq['HO_DISSIPATION_FACTOR'] = 1 # a constant (non-negative) factor adapting the influence of the high order dissipation
#Additional parameters needed for adaptive dissipation. For typical values
#see Svärd et al. (2009).
I_IEq['MP2MIN'] = 1
I_IEq['MP2MAX'] = 10
I_IEq['CMIN'] = 1
I_IEq['CMAX'] = 10


# Divergence cleaning related variables
#Enable divergence cleaning (True) or disable divergence cleaning (False)
I_DC['divergence_cleaning'] = False
#Choose how the laplace operator will be discretized
I_DC['divergence_cleaning_form'] = 'USE_LAPLACE_WIDE_STENCIL_LNS'
#USE_LAPLACE_WIDE_STENCIL_LNS USE_LAPLACE_WIDE_STENCIL_DIRICHLET USE_LAPLACE_NARROW_STENCIL_DIRICHLET
#Specify a threshold for the divergence norm
I_DC['absolute_error_threshold'] = 1e-3
#The divergence cleaner will exit after max_iterations even if the error
#threshold is not reached
I_DC['max_iterations'] = 50


# Technical details
# Specify the OpenCL device that will be used for computation
I_Tech['device'] = 1

# Switch between single floating point and double floating point precision
I_Tech['REAL'] = 'double' # float
I_Tech['REAL4'] = I_Tech['REAL'] + '4' #Vector datatype

#Compiler based optimizations
if I_Tech['REAL'] == 'float':
    I_Tech['optimizations'] = ' -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only -cl-single-precision-constant'
else:
    I_Tech['optimizations'] = ' -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only'


# Options relevant for a run
# Defines the order
I_RunOps['order'] = 4 # 2, 4, 6
I_RunOps['operator_form'] = 'classical' # 'classical' or 'extended' operators
# Specify the testcase. The name of the testcase has to be equal to the
# name of a header file which contains a function describing the initial state.
# Example testcases are:
# rotation_2D, rotation_3D, alfven_periodic_2D, hall_travelling_wave, hall_periodic
induction_eq.I_RunOps['testcase'] = 'rotation_2D';
I_RunOps['variable_u'] = False; # must be set to True if a variable velocity is used
I_RunOps['periodic'] = 'NONE'; # 'NONE', 'USE_PERIODIC'; must be set to 'USE_PERIODIC'
                                       # if periodic boundary conditions should be used
# Optional plotting parameters (2D plots).
# Choose the cross section with
# 'x', 'y', 'z'
# If you want to plot multiple cross sections use
# 'xy', 'xz', 'yz', 'xyz'
I_RunOps['plot_numerical_solution'] = ''
I_RunOps['plot_analytical_solution'] = ''
I_RunOps['plot_difference'] = ''
I_RunOps['plot_divergence'] = ''
#If set to 1 the magnetic field will be saved to I_Results('field_b')
I_RunOps['save_fields'] = False
#If set to true the magnetic energy (L^2 norm) will be saved to I_Results('energy_over_time'), the
#L2 errors (if available) to I_Results('L2error_B_over_time') & I_Results('L2error_divB_over_time')
#and the time to I_Results('time')
I_RunOps['save_integrals_over_time'] = False

#Initialize the magnetic field, the velocity field and the density field
#according to the specified testcase. Also calculates and sets additional
#variables, e.g. stepsize, timestep, local work group size, etc...
(field_b_init, field_u_init, field_rho_init) = induction_eq.initialize()

print('Testcase: ' + I_RunOps['testcase'])
print('Order:' + str(I_RunOps['order']))
print('Time integrator:' + I_TI['time_integrator'])
print('DT: '+ str(I_TI['DT']) + '   N_STEPS: ' + str(I_TI['num_steps']) + '   FINAL_TIME: '+str(I_TI['final_time']))
print('DX: '+ str(I_Mesh['DX']) + '   NODES_X: '+str(I_Mesh['NODES_X']))
print('DY: '+ str(I_Mesh['DY']) + '   NODES_Y: '+str(I_Mesh['NODES_Y']))
print('DZ: '+ str(I_Mesh['DZ']) + '   NODES_Z: '+str(I_Mesh['NODES_Z']))
print('Dissipation: '+I_IEq['dissipation'])
print('Divergence cleaning: '+str(I_DC['divergence_cleaning']))
print('REAL:'+ I_Tech['REAL'])

# Takes the initialized fields and advances the solution in time
induction_eq.compute_numerical_solution(field_b_init, field_u_init, field_rho_init);

print('Divergence Norm: ' + str(I_Results['divergence_norm']) + ' Total Energy: '+str(I_Results['energy']))
print('Relative Error:' + str(100*I_Results['rel_err']) + '%\n\n')
