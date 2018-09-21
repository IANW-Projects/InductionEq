%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

%clc;
clear;
close all;

addpath('../matlab')

induction_eq.prepare_vars();

global I_Mesh I_TI I_IEq I_DC I_Tech I_RunOps I_Results

% Mesh related parameters. The minimum number of nodes per direction is
% constrained by the number of boundary nodes, e.g. the 4th order method
% has 2*4 boundary nodes which means the minimum number of nodes amounts to 8.
N = uint32(40);
I_Mesh('NODES_X') = N; I_Mesh('NODES_Y') = N; I_Mesh('NODES_Z') = N;
I_Mesh('XMIN') = -1.0; I_Mesh('XMAX') = 1.0;
I_Mesh('YMIN') = -1.0; I_Mesh('YMAX') = 1.0;
I_Mesh('ZMIN') = -1.0; I_Mesh('ZMAX') = 1.0;

% Time integration related variables
I_TI('cfl') = 0.95; %Define the Courant–Friedrichs–Lewy condition
I_TI('final_time') = 2*pi;
%Chose the time integrator. Below is a list of up to date available
%options:
% SSPRK33, SSPRK104,
% KennedyCarpenterLewis2R54C, CalvoFrancoRandez2R64,
% CarpenterKennedy2N54, ToulorgeDesmet2N84F
I_TI('time_integrator') = 'CarpenterKennedy2N54';

% Induction equation related variables
%Specify how the three part of the linear induction equation shall be
%computed.
I_IEq('form_uibj') = 'USE_UIBJ_PRODUCT'; % PRODUCT, SPLIT, CENTRAL
I_IEq('form_source') = 'USE_SOURCE_CENTRAL'; % CENTRAL, SPLIT, ZERO
I_IEq('form_ujbi') = 'USE_UJBI_CENTRAL'; % SPLIT, CENTRAL, PRODUCT
%Enable or disable Hall-Term
I_IEq('hall_term') = 'NONE'; % NONE, USE_HALL
%Enable or disable artificial dissipation
I_IEq('dissipation') = 'NONE'; %NONE, USE_ARTIFICIAL_DISSIPATION
%Specify what kind of artifical dissipation shall be used
%USE_ADAPTIVE_DISSIPATION, USE_FIRST_ORDER_DISSIPATION, USE_HIGH_ORDER_DISSIPATION
I_IEq('dissipation_form') = 'USE_ADAPTIVE_DISSIPATION';
%Additional parameters needed for adaptive dissipation. For typical values
%see Svärd et al. (2009).
I_IEq('MP2MIN') = 1;
I_IEq('MP2MAX') = 10;
I_IEq('CMIN') = 1;
I_IEq('CMAX') = 10;

% Divergence cleaning related variables
%Enable divergence cleaning (true) or disable divergence cleaning (false)
I_DC('divergence_cleaning') = false;
%Choose how the laplace operator will be discretized
I_DC('divergence_cleaning_form') = 'USE_LAPLACE_WIDE_STENCIL_LNS';
%USE_LAPLACE_WIDE_STENCIL_LNS USE_LAPLACE_WIDE_STENCIL_DIRICHLET USE_LAPLACE_NARROW_STENCIL_DIRICHLET
%Specify a threshold for the divergence norm
I_DC('absolute_error_threshold') = 1e-3;
%The divergence cleaner will exit after max_iterations even if the error
%threshold is not reached
I_DC('max_iterations') = 50;

% Technical details
% Specify the OpenCL device that will be used for computation
I_Tech('device') = 1;

%Switch between single floating point and double floating point precision
I_Tech('REAL') = 'double'; % float
I_Tech('REAL4') = sprintf('%s4',I_Tech('REAL')); %Vector datatype

%Compiler based optimizations
if strcmp(I_Tech('REAL'),'float')
    I_Tech('optimizations') = ' -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only -cl-single-precision-constant';
else
    I_Tech('optimizations') = ' -cl-mad-enable -cl-no-signed-zeros -cl-finite-math-only';
end

% Options relevant for a run
% Defines the order
I_RunOps('order') = 4; % 2, 4, 6
I_RunOps('operator_form') = 'classical'; % 'classical' or 'extended' operators
% Specify the testcase. The name of the testcase has to be equal to the
% name of a header file which contains a function describing the initial state.
% Example testcases are:
% rotation_2D, rotation_3D, alfven_periodic_2D, hall_travelling_wave, hall_periodic
I_RunOps('testcase') = 'rotation_2D';
I_RunOps('variable_u') = false; % must be set to true if a variable velocity is used
I_RunOps('periodic') = 'NONE'; % 'NONE', 'USE_PERIODIC'; must be set to 'USE_PERIODIC'
                                       % if periodic boundary conditions should be used
% Optional plotting parameters (2D plots).
% Choose the cross section with
% 'x', 'y', 'z'
% If you want to plot multiple cross sections use
% 'xy', 'xz', 'yz', 'xyz'
I_RunOps('plot_numerical_solution') = '';
I_RunOps('plot_analytical_solution') = '';
I_RunOps('plot_difference') = '';
I_RunOps('plot_divergence') = '';
%If set to 1 the magnetic field will be saved to I_Results('field_b')
I_RunOps('save_fields') = false;


Ns = [uint32(40), uint32(80), uint32(160), uint32(320), ];
orders = [2, 4, 6];
forms_uiBj = {'USE_UIBJ_PRODUCT', 'USE_UIBJ_SPLIT', 'USE_UIBJ_CENTRAL'};
forms_source = {'USE_SOURCE_ZERO', 'USE_SOURCE_SPLIT', 'USE_SOURCE_CENTRAL'};
forms_ujBi = {'USE_UJBI_PRODUCT', 'USE_UJBI_SPLIT', 'USE_UJBI_CENTRAL'};

forms = { ...
    {'USE_UIBJ_CENTRAL', 'USE_SOURCE_ZERO', 'USE_UJBI_CENTRAL'}, ...
    {'USE_UIBJ_CENTRAL', 'USE_SOURCE_CENTRAL', 'USE_UJBI_CENTRAL'}, ...
    {'USE_UIBJ_SPLIT', 'USE_SOURCE_CENTRAL', 'USE_UJBI_SPLIT'}, ...
    {'USE_UIBJ_PRODUCT', 'USE_SOURCE_CENTRAL', 'USE_UJBI_PRODUCT'}, ...
    {'USE_UIBJ_PRODUCT', 'USE_SOURCE_CENTRAL', 'USE_UJBI_SPLIT'}, ...
    {'USE_UIBJ_PRODUCT', 'USE_SOURCE_CENTRAL', 'USE_UJBI_CENTRAL'}, ...
};

io = fopen('rotation_3D_convergence.txt', 'w');
fprintf(io, '# N, order, form_uiBj, form_source, form_ujBi, runtime, energy, error in B, error in div B \n');

for N = Ns
  for order = orders
    for form = forms
      form_uiBj = form{1}{1};
      form_source = form{1}{2};
      form_ujBi = form{1}{3};
      
      fprintf('N = %4d, order = %d, form_uiBj = %16s, form_source = %18s, form_ujBi = %16s\n', ...
              N, order, char(form_uiBj), char(form_source), char(form_ujBi));

      I_Mesh('NODES_X') = N; I_Mesh('NODES_Y') = N; I_Mesh('NODES_Z') = N;
      I_IEq('form_uibj') = char(form_uiBj);
      I_IEq('form_source') = char(form_source);
      I_IEq('form_ujbi') = char(form_ujBi);
      I_RunOps('order') = order;

      [field_b_init, field_u_init, field_rho_init] = induction_eq.initialize();
      induction_eq.compute_numerical_solution(field_b_init, field_u_init, field_rho_init);
          
      fprintf(io, '%d, %d, %s, %s, %s, %.15e, %.15e, %.15e, %.15e \n', ...
              N, order, char(form_uiBj), char(form_source), char(form_ujBi), ...
              I_Results('runtime'), I_Results('energy'), I_Results('abs_err'), I_Results('divergence_norm'));

      fprintf('  Runtime: %8.2f s, Error in B: %.3e, Error in div B: %.3e \n\n', ...
              I_Results('runtime'), I_Results('abs_err'), I_Results('divergence_norm'));
    end
  end
end

fclose(io);

