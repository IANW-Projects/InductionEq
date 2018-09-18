%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

function [field_b_init, DC_fields] = initialize()

    %Initialize the magnetic field and additional fields needed for divergence
    %cleaning.Also calculates and sets additional variables, e.g. stepsize,
    %local work group size, etc...

    global I_Mesh I_DC I_Tech I_RunOps

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
    I_Tech('W_SIZE') = uint32(group_size);

    % initialise fields: magnetic, divergence cleaning fields
    if strcmp(I_Tech('REAL'),'float')
        field_b_init = single(zeros(4, I_Tech('num_nodes_pad')));
        field_divB = single(zeros(1, I_Tech('num_nodes_pad')));
        field_phi  = single(zeros(1, I_Tech('num_nodes_pad')));

        %Fields needed for CG
        field_laplace = single(zeros(1, I_Tech('num_nodes_pad')));
        field_r = single(zeros(1, I_Tech('num_nodes_pad')));
        field_u = single(zeros(1, I_Tech('num_nodes_pad')));
        output  = single(zeros(1, I_Tech('num_groups')));
    else
        field_b_init = double(zeros(4, I_Tech('num_nodes_pad')));
        field_divB = double(zeros(1, I_Tech('num_nodes_pad')));
        field_phi  = double(zeros(1, I_Tech('num_nodes_pad')));

        %Fields needed for CG
        field_laplace = double(zeros(1, I_Tech('num_nodes_pad')));
        field_r = double(zeros(1, I_Tech('num_nodes_pad')));
        field_u = double(zeros(1, I_Tech('num_nodes_pad')));
        output  = double(zeros(1, I_Tech('num_groups')));
    end

    %TODO
    rng(12345);
    fac = 1.e-3;
    for iz=1:I_Mesh('NODES_Z')
        for iy=1:I_Mesh('NODES_Y')
            for ix=1:I_Mesh('NODES_X')
                idx=(ix + (iy -1)*I_Mesh('NODES_X') + (iz -1)*I_Mesh('NODES_Y')*I_Mesh('NODES_X'));

                x = I_Mesh('XMIN') + I_Mesh('DX')*(ix-1);
                y = I_Mesh('YMIN') + I_Mesh('DY')*(iy-1);
                z = I_Mesh('ZMIN') + I_Mesh('DZ')*(iz-1);

                field_b_init(:,idx) = [-y+fac*rand(), ...
                                        x+fac*rand(), ...
                                        0+fac*rand(), ...
                                        0];
            end
        end
    end


    % Set number of boundary nodes. Depends on the discretization of the
    % laplace operator.
    if strcmp(I_DC('divergence_cleaning_form'), 'USE_LAPLACE_WIDE_STENCIL_LNS')
        I_DC('BNODES') = 0;
    elseif strcmp(I_DC('divergence_cleaning_form'), 'USE_LAPLACE_WIDE_STENCIL_DIRICHLET') || ...
           strcmp(I_DC('divergence_cleaning_form'), 'USE_LAPLACE_NARROW_STENCIL_DIRICHLET')
        I_DC('BNODES') = 1;
    else
        error('Option I_DC(''divergence_cleaning_form'') = %s unknown.', I_DC('divergence_cleaning_form'))
    end

    %Generate settings
    settings_tech = generate_settings(I_Tech, {'REAL'; 'REAL4'; 'W_SIZE'; 'optimizations'});
    settings_mesh = generate_settings(I_Mesh, {'DX'; 'DY'; 'DZ';...
                                               'NODES_X'; 'NODES_Y'; 'NODES_Z'});
    settings_dc = generate_settings(I_DC, {'BNODES', 'divergence_cleaning_form'});

    settings = strcat(settings_tech, settings_mesh, settings_dc);

    % Fields for divergence cleaning
    DC_fields = struct('divB',field_divB,...
                       'phi',field_phi,...
                       'laplace', field_laplace,...
                       'r', field_r, ...
                       'u', field_u',...
                       'output', output);

    kernel_path_list =  {};
    if I_RunOps('order') == 2
        kernel_path_list = [kernel_path_list, {'../include/2ndOrder.h'}];
    elseif I_RunOps('order') == 4
        kernel_path_list = [kernel_path_list, {'../include/4thOrder.h'}];
    elseif I_RunOps('order') == 6
        kernel_path_list = [kernel_path_list, {'../include/6thOrder.h'}];
    else
        fprintf('Specify order \n')
    end
    kernel_path_list = [kernel_path_list, {'../include/utils.h'}];
    kernel_path_list = [kernel_path_list, {'../kernel/SBP_operator.cl'}];
    kernel_path_list = [kernel_path_list, {'../kernel/kernel_operator.cl'}];
    kernel_path_list = [kernel_path_list, {'../kernel/kernel_div_cleaning.cl'}];
    kernel_path_list = [kernel_path_list, {'../kernel/kernel_CG.cl'}];
    kernel_path_list = [kernel_path_list, {'../kernel/kernel_dot_product.cl'}];

    % Compile OpenCL kernels
    cl_run_kernel(I_Tech('device'), kernel_path_list, settings);

    % Define global and local range for divergence cleaning and dot product
    I_DC('g_range') = uint32([I_Mesh('NODES_X')-2*I_DC('BNODES'), ...
                              I_Mesh('NODES_Y')-2*I_DC('BNODES'), ...
                              I_Mesh('NODES_Z')-2*I_DC('BNODES')]);
    I_DC('l_range') = uint32([0]);
    I_Tech('g_range') = uint32([I_Tech('num_nodes_pad'), 1, 1]);
    I_Tech('l_range') = uint32([I_Tech('W_SIZE'), 1, 1]);
end
