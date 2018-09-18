%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

function [field_b_init, field_u_init, field_rho_init] = initialize()

    %Initialize the magnetic field, the velocity field and the density field
    %according to the specified testcase. Also calculates and sets additional
    %variables, e.g. stepsize, timestep, local work group size, etc...

    global I_Mesh I_TI I_IEq I_DC I_Tech I_RunOps

    % Calculate the stepsize for each dimension
    if strcmp(I_RunOps('periodic'), 'USE_PERIODIC')
        % the nodes at the right boundary are not included
        DX = double(I_Mesh('XMAX') - I_Mesh('XMIN')) / double(I_Mesh('NODES_X'));
        DY = double(I_Mesh('YMAX') - I_Mesh('YMIN')) / double(I_Mesh('NODES_Y'));
        DZ = double(I_Mesh('ZMAX') - I_Mesh('ZMIN')) / double(I_Mesh('NODES_Z'));
    else
        DX = double(I_Mesh('XMAX') - I_Mesh('XMIN')) / (double(I_Mesh('NODES_X'))-1);
        DY = double(I_Mesh('YMAX') - I_Mesh('YMIN')) / (double(I_Mesh('NODES_Y'))-1);
        DZ = double(I_Mesh('ZMAX') - I_Mesh('ZMIN')) / (double(I_Mesh('NODES_Z'))-1);
    end

    I_Mesh('DX') = DX; I_Mesh('DY') = DY; I_Mesh('DZ') = DZ;

    num_nodes = I_Mesh('NODES_X')*I_Mesh('NODES_Y')*I_Mesh('NODES_Z');

    % Depending on device type chose the local work group size. If the
    % device is a gpu the work group size is taken from the device
    % specifications. In case of a cpu the optimal work group size was determined
    % through testing. Depending on the architecture this value may change.
    [~, dev_type, ~, ~, lw_size, ~] = cl_get_devices;
    type = dev_type(I_Tech('device'));
    if (strcmp(type{1},'CPU'))
        group_size = 16;
    else
        if (num_nodes > lw_size(I_Tech('device')))
            group_size = lw_size(I_Tech('device'));
        else
            group_size = 2^floor(log(num_nodes) / log(2));
        end
    end
    group_size = double(group_size);
    num_nodes_pad = ceil(double(num_nodes)/group_size)*group_size;
    num_groups = ceil(num_nodes_pad/group_size);

    I_Tech('num_nodes_pad') = num_nodes_pad;
    I_Tech('num_groups') = num_groups;
    I_Tech('W_SIZE') = uint16(group_size);

    % Initialise fields: magnetic, veclocity, and charge density
    if strcmp(I_Tech('REAL'),'float')
        field_u_init = single(zeros(4, num_nodes_pad));
        field_b_init = single(zeros(4, num_nodes_pad));
        field_rho_init = single(ones(1, num_nodes_pad));
    else
        field_u_init = double(zeros(4, num_nodes_pad));
        field_b_init = double(zeros(4, num_nodes_pad));
        field_rho_init = double(ones(1, num_nodes_pad));
    end

    % Generate settings needed for computation of the initial state for the
    % specified testcase.
    settings_tech = generate_settings(I_Tech, {'REAL'; 'REAL4'; 'optimizations'});
    settings_mesh = generate_settings(I_Mesh, {'DX'; 'DY'; 'DZ';...
                                               'NODES_X'; 'NODES_Y'; 'NODES_Z'; ...
                                               'XMIN'; 'XMAX'; 'YMIN'; 'YMAX';'ZMIN'; 'ZMAX'});

    settings = strcat(settings_tech, settings_mesh);

    % Include all OpenCL and header files needed for compilation.
    kernel_path_list =  {};
    kernel_path_list = [kernel_path_list, {sprintf('../include/%s.h', I_RunOps('testcase'))}];
    kernel_path_list = [kernel_path_list, {'../include/utils.h'}];
    kernel_path_list = [kernel_path_list, {'../kernel/kernel_init.cl'}];

    % Comile kernel with setting on device
    cl_run_kernel(I_Tech('device'), kernel_path_list, settings);

    % Set global and local range for calculation of the inital state
    % The global and local range can either be a vector of type uint32 or
    % double. Double is supported due to the fact that by default Matlab
    % stores all variables as double.
    I_IEq('g_range') = uint32([I_Mesh('NODES_X'), I_Mesh('NODES_Y'), I_Mesh('NODES_Z')]);
    I_IEq('l_range') = uint32([0]);

    % Compute initial state for magnetic field, velocity field and density
    % field
    cl_run_kernel(I_Tech('device'), 'init_b', I_IEq('g_range'), I_IEq('l_range'), field_b_init, 0);
    cl_run_kernel(I_Tech('device'), 'init_u', I_IEq('g_range'), I_IEq('l_range'), field_u_init, 0);
    cl_run_kernel(I_Tech('device'), 'init_rho', I_IEq('g_range'), I_IEq('l_range'), field_rho_init, 0);

    % Compute time step from CFL number, min step size and maximal velocity
    u_max = max(sqrt(field_u_init(1,:).^2 + field_u_init(2,:).^2 + field_u_init(3,:).^2));
    dt = I_TI('cfl') * min([DX DY DZ]) / double(u_max);
    num_steps = ceil(I_TI('final_time')/dt);
    dt = I_TI('final_time') / num_steps;

    I_TI('DT') = dt;
    I_TI('num_steps') = num_steps;

    % Global and local range for dot product and norm
    I_Tech('g_range') = uint32([I_Tech('num_nodes_pad'), 1, 1]);
    I_Tech('l_range') = uint32([I_Tech('W_SIZE'), 1, 1]);

    % Depending on the discretization of the laplace operator for
    % divergence cleaning set the number of boundary nodes
    if strcmp(I_DC('divergence_cleaning_form'), 'USE_LAPLACE_WIDE_STENCIL_LNS')
        I_DC('BNODES') = 0;
    elseif strcmp(I_DC('divergence_cleaning_form'), 'USE_LAPLACE_WIDE_STENCIL_DIRICHLET') || ...
           strcmp(I_DC('divergence_cleaning_form'), 'USE_LAPLACE_NARROW_STENCIL_DIRICHLET')
        I_DC('BNODES') = 1;
    else
        error('Option I_DC(''divergence_cleaning_form'') = %s unknown.', I_DC('divergence_cleaning_form'))
    end
    % Set global and local range for divergence cleaning
    I_DC('g_range') = uint32([I_Mesh('NODES_X')-2*I_DC('BNODES'), ...
                              I_Mesh('NODES_Y')-2*I_DC('BNODES'), ...
                              I_Mesh('NODES_Z')-2*I_DC('BNODES')]);
    I_DC('l_range') = uint32([0]);
end
