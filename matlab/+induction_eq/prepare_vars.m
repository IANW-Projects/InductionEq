%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

function [] = prepare_vars()

    %prepare_vars() initializes containers (maps) in which all variables (label, value) are
    %stored. Variables are grouped together by purpose:
    %I_Mesh: Contains all variables concerning the discrete mesh, e.g.
    %stepsize, box size
    %I_TI: Contains all variables concerning time integration, e.g. time
    %integrator, timestep
    %I_IEq: Contains all variables relevant for the induction equation, e.g.
    %discretization form, switch for additional artificial dissipation
    %I_DC: Contains all variables concerning divergence cleaning, e.g.
    %discretization, error threshold
    %I_Tech: Contains all variables of technical nature, e.g. the device,
    %computation precision
    %I_RunOps: Specifies the parameters of a computation, e.g. which variables
    %will be saved, the testcase

    %Variables in capital letters are program specific defines which are set by
    %openCL compile settings. If the value is a string in captial letters
    %beginning with USE it acts as a switch, e.g. to enable artificial dissipation
    %or to switch between different discretizations. In case the key is written
    %in capital letters the program define will be set to the corresponding
    %value of the key.

    global I_Mesh I_TI I_IEq I_DC I_Tech I_RunOps I_Results

    keySet = {'NODES_X', 'NODES_Y', 'NODES_Z', 'DX', 'DY', 'DZ', 'XMIN','YMIN', 'ZMIN', 'XMAX', 'YMAX','ZMAX'};
    valueSet = {uint16(0) uint16(0) uint16(0) 0 0 0 0 0 0 0 0 0};
    I_Mesh = containers.Map(keySet, valueSet, 'UniformValues',false);

    keySet = {'cfl', 'final_time', 'time_integrator', 'DT', 'num_steps'};
    valueSet = {0 0 '' 0 0};
    I_TI = containers.Map(keySet, valueSet,'UniformValues',false);

    keySet = {'form_uibj', 'form_source', 'form_ujbi', 'dissipation', 'dissipation_form', 'HO_DISSIPATION_FACTOR', ...
              'hall_term', 'g_range', 'l_range', 'MP2MIN', 'MP2MAX', 'CMIN', 'CMAX'};
    valueSet = {'' '' '' 'NONE' '' 1 'NONE' 0 0 0 0 0 0};
    I_IEq = containers.Map(keySet, valueSet,'UniformValues',false);

    keySet = {'use', 'form', 'BNODES', 'error_threshold', 'g_range', 'l_range','max_iterations'};
    valueSet = {0 '' 0 0 0 0 0};
    I_DC = containers.Map(keySet, valueSet,'UniformValues',false);

    keySet = {'device', 'REAL','num_nodes_pad', 'num_groups','W_SIZE', 'g_range', 'l_range', 'optimizations'};
    valueSet = {0 '' 0 0 0 0 0 ''};
    I_Tech = containers.Map(keySet, valueSet,'UniformValues',false);

    keySet = {'order', 'operator_form', 'testcase', 'variable_u', 'periodic', 'plot_energy', 'plot_numerical_solution', 'plot_analytical_solution', 'plot_difference', 'plot_divergence', 'save_fields'};
    valueSet = {0 'classical' '' false 'NONE' 0 0 0 0 0 false};
    I_RunOps = containers.Map(keySet, valueSet,'UniformValues',false);

    keySet = {'abs_err', 'rel_err', 'field_b', 'field_b_ana', 'energy', 'divergence_norm', 'runtime', 'kernel_runtime'};
    valueSet = {0 0 0 0 0 0 0 0};
    I_Results = containers.Map(keySet, valueSet,'UniformValues',false);

end
