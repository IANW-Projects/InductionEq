%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

function [kernel_path_list, settings, IEq_fields, time_integrator_num_fields, DC_fields, RK_Step] = compute_numerical_solution_setup(field_b, field_u, field_rho)

global I_Mesh I_TI I_IEq I_DC I_Tech I_RunOps

% Memory allocation for fields
if strcmp(I_Tech('REAL'),'float')
    % Initialize auxillary fields needed for time stepping
    field_b2 = single(zeros(4, I_Tech('num_nodes_pad')));
    field_b3 = single(zeros(4, I_Tech('num_nodes_pad')));
    field_divB = single(zeros(1, I_Tech('num_nodes_pad')));
    field_curlB_rho = single(zeros(4, I_Tech('num_nodes_pad')));
    current_time = single(zeros(2));
    norm2_output = single(zeros(1, I_Tech('num_groups')));

    % Allocate memory for field used for divergence cleaning
    if I_DC('divergence_cleaning')
        %Field for divergence of b-field and phi
        field_divB = single(zeros(1, I_Tech('num_nodes_pad')));
        field_phi = single(zeros(1, I_Tech('num_nodes_pad')));

        %Fields needed for CG
        field_laplace = single(zeros(1, I_Tech('num_nodes_pad')));
        field_r = single(zeros(1, I_Tech('num_nodes_pad')));
        field_u_div_cleaning = single(zeros(1, I_Tech('num_nodes_pad')));
    end
else
    % Initialize auxillary fields needed for time stepping
    field_b2 = double(zeros(4, I_Tech('num_nodes_pad')));
    field_b3 = double(zeros(4, I_Tech('num_nodes_pad')));
    field_divB = double(zeros(1, I_Tech('num_nodes_pad')));
    field_curlB_rho = double(zeros(4, I_Tech('num_nodes_pad')));
    current_time = double(zeros(2));
    norm2_output = double(zeros(1, I_Tech('num_groups')));

    % Allocate memory for field used for divergence cleaning
    if I_DC('divergence_cleaning')
        %Field for divergence of b-field and phi
        field_divB = double(zeros(1, I_Tech('num_nodes_pad')));
        field_phi = double(zeros(1, I_Tech('num_nodes_pad')));

        %Fields needed for CG
        field_laplace = double(zeros(1, I_Tech('num_nodes_pad')));
        field_r = double(zeros(1, I_Tech('num_nodes_pad')));
        field_u_div_cleaning = double(zeros(1, I_Tech('num_nodes_pad')));
    end
end

kernel_path_list = {};

% Include header file containing the coefficients of the respective order
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
kernel_path_list = [kernel_path_list, {sprintf('../include/%s.h', I_RunOps('testcase'))}];
kernel_path_list = [kernel_path_list, {'../kernel/kernel_init.cl'}];

kernel_path_list = [kernel_path_list, {'../kernel/SBP_operator.cl'}];
kernel_path_list = [kernel_path_list, {'../kernel/kernel_operator.cl'}];

% Include files for divergence cleaning if desired
if I_DC('divergence_cleaning')
    %Prepare temporary fields for divergence cleaning
    DC_fields = struct('divB',field_divB, ...
                       'phi',field_phi, ...
                       'laplace', field_laplace, ...
                       'r', field_r, ...
                       'u', field_u_div_cleaning', ...
                       'output', norm2_output);
    kernel_path_list = [kernel_path_list, {'../kernel/kernel_div_cleaning.cl'}];
    kernel_path_list = [kernel_path_list, {'../kernel/kernel_CG.cl'}];
    settings_dc = generate_settings(I_DC, {'BNODES', 'divergence_cleaning_form'});
else
    DC_fields = 0;
    settings_dc ='';
end

% Include files for artificial dissipation if desired
if strcmp(I_IEq('dissipation'), 'USE_ARTIFICIAL_DISSIPATION')
    kernel_path_list = [kernel_path_list, {'../include/artificial_dissipation.h'}];
    kernel_path_list = [kernel_path_list, {'../kernel/artificial_dissipation.cl'}];
    settings_dissipation = generate_settings(I_IEq, {'dissipation', 'dissipation_form'});

    if strcmp(I_IEq('dissipation_form'),'USE_ADAPTIVE_DISSIPATION')
        settings_dissipation = strcat(settings_dissipation, generate_settings(I_IEq, {'MP2MIN'; 'MP2MAX'; 'CMIN'; 'CMAX'}));
    end
else
    settings_dissipation = '';
end

kernel_path_list = [kernel_path_list, {'../kernel/induction_Eq_volume.cl'}];
kernel_path_list = [kernel_path_list, {'../kernel/induction_Eq_surface.cl'}];
kernel_path_list = [kernel_path_list, {'../kernel/kernel_dot_product.cl'}];
kernel_path_list = [kernel_path_list, {'../kernel/kernel_time_integrator.cl'}];

% Specify compile settings
settings_tech = generate_settings(I_Tech, {'REAL'; 'REAL4'; 'W_SIZE'; 'optimizations'});

settings_mesh = generate_settings(I_Mesh, {'DX'; 'DY'; 'DZ';...
                                           'NODES_X'; 'NODES_Y'; 'NODES_Z'; ...
                                           'XMIN'; 'XMAX'; 'YMIN'; 'YMAX';'ZMIN'; 'ZMAX'});

settings_induction_eq = generate_settings(I_IEq, {'form_uibj'; 'form_source'; 'form_ujbi'; 'hall_term'});

settings_time_integration = generate_settings(I_TI, {'DT'});

settings = strcat(settings_induction_eq, settings_time_integration, settings_tech, settings_mesh, settings_dissipation, settings_dc);


switch I_TI('time_integrator')

    case 'SSPRK33'
        RK_Step = {'SSPRK33_1', 'SSPRK33_2', 'SSPRK33_3'};
        if strcmp(I_IEq('hall_term'), 'USE_HALL')
            calc_curlB_rho = {'calc_curlB_rho_1_3args', 'calc_curlB_rho_2_3args', 'calc_curlB_rho_3_3args'};
        else
            calc_curlB_rho = [];
        end
        if I_RunOps('variable_u') == 1
            calc_u = repmat({'calc_u_3_args'}, size(RK_Step));
        else
            calc_u = [];
        end
        tmp = [calc_curlB_rho; calc_u; RK_Step];
        RK_Step = tmp(:)';
        RK_Step{end+1} = 'calc_time_3_args';
        time_integrator_num_fields = 3;

    case 'SSPRK104'

        RK_Step_1 = {'SSPRK104_01', 'SSPRK104_02', 'SSPRK104_03', 'SSPRK104_04', 'SSPRK104_05'};
        RK_Step_2 = {'SSPRK104_07', 'SSPRK104_08', 'SSPRK104_09', 'SSPRK104_10', 'SSPRK104_11'};
        if strcmp(I_IEq('hall_term'), 'USE_HALL')
            calc_curlB_rho_1 = {'calc_curlB_rho_1_3args', 'calc_curlB_rho_2_3args', 'calc_curlB_rho_3_3args', ...
                                'calc_curlB_rho_2_3args', 'calc_curlB_rho_3_3args'};
            calc_curlB_rho_2 = {'calc_curlB_rho_2_3args', 'calc_curlB_rho_1_3args', 'calc_curlB_rho_2_3args', ...
                                'calc_curlB_rho_1_3args', 'calc_curlB_rho_2_3args'};
        else
            calc_curlB_rho_1 = [];
            calc_curlB_rho_2 = [];
        end
        if I_RunOps('variable_u')
            calc_u = repmat({'calc_u_3_args'}, size(RK_Step_1));
        else
            calc_u = [];
        end
        tmp = [calc_curlB_rho_1; calc_u; RK_Step_1];
        RK_Step = tmp(:)';
        RK_Step{end+1} = 'SSPRK104_06';
        tmp = [calc_curlB_rho_2; calc_u; RK_Step_2];
        RK_Step = [RK_Step, tmp(:)'];
        RK_Step{end+1} = 'calc_time_3_args';
        time_integrator_num_fields = 3;

    case 'KennedyCarpenterLewis2R54C'

        RK_Step_a = {'KennedyCarpenterLewis2R54C_1a', 'KennedyCarpenterLewis2R54C_2a', ...
                     'KennedyCarpenterLewis2R54C_3a', 'KennedyCarpenterLewis2R54C_4a', ...
                     'KennedyCarpenterLewis2R54C_5'};
        RK_Step_b = {'KennedyCarpenterLewis2R54C_1b', 'KennedyCarpenterLewis2R54C_2b', ...
                     'KennedyCarpenterLewis2R54C_3b', 'KennedyCarpenterLewis2R54C_4b', ...
                     'calc_time_3_args'};
        if strcmp(I_IEq('hall_term'), 'USE_HALL')
            calc_curlB_rho = {'calc_curlB_rho_1_3args', 'calc_curlB_rho_1_3args', 'calc_curlB_rho_2_3args', ...
                              'calc_curlB_rho_1_3args', 'calc_curlB_rho_2_3args'};
        else
            calc_curlB_rho = [];
        end
        if I_RunOps('variable_u')
            calc_u = repmat({'calc_u_3_args'}, size(RK_Step_a));
        else
            calc_u = [];
        end
        tmp = [calc_curlB_rho; calc_u; RK_Step_a; RK_Step_b];
        RK_Step = tmp(:)';
        time_integrator_num_fields = 3;

    case 'CalvoFrancoRandez2R64'

        RK_Step_a = {'CalvoFrancoRandez2R64_1a', 'CalvoFrancoRandez2R64_2a', ...
                     'CalvoFrancoRandez2R64_3a', 'CalvoFrancoRandez2R64_4a', ...
                     'CalvoFrancoRandez2R64_5a', 'CalvoFrancoRandez2R64_6'};
        RK_Step_b = {'CalvoFrancoRandez2R64_1b', 'CalvoFrancoRandez2R64_2b', ...
                     'CalvoFrancoRandez2R64_3b', 'CalvoFrancoRandez2R64_4b', ...
                     'CalvoFrancoRandez2R64_5b', 'calc_time_3_args'};
        if strcmp(I_IEq('hall_term'), 'USE_HALL')
            calc_curlB_rho = {'calc_curlB_rho_1_3args', 'calc_curlB_rho_2_3args', 'calc_curlB_rho_2_3args', ...
                              'calc_curlB_rho_2_3args', 'calc_curlB_rho_2_3args', 'calc_curlB_rho_2_3args'};
        else
            calc_curlB_rho = [];
        end
        if I_RunOps('variable_u')
            calc_u = repmat({'calc_u_3_args'}, size(RK_Step_a));
        else
            calc_u = [];
        end
        tmp = [calc_curlB_rho; calc_u; RK_Step_a; RK_Step_b];
        RK_Step = tmp(:)';
        time_integrator_num_fields = 3;

    case 'CarpenterKennedy2N54'

        RK_Step_a = {'CarpenterKennedy2N54_1a', 'CarpenterKennedy2N54_2a', ...
                     'CarpenterKennedy2N54_3a', 'CarpenterKennedy2N54_4a', ...
                     'CarpenterKennedy2N54_5a'};
        RK_Step_b = {'CarpenterKennedy2N54_1b', 'CarpenterKennedy2N54_2b', ...
                     'CarpenterKennedy2N54_3b', 'CarpenterKennedy2N54_4b', ...
                     'CarpenterKennedy2N54_5b'};
        if strcmp(I_IEq('hall_term'), 'USE_HALL')
            calc_curlB_rho = repmat({'calc_curlB_rho_1_2args'}, size(RK_Step_a));
        else
            calc_curlB_rho = [];
        end
        if I_RunOps('variable_u')
            calc_u = repmat({'calc_u_2_args'}, size(RK_Step_a));
        else
            calc_u = [];
        end
        tmp = [calc_curlB_rho; calc_u; RK_Step_a; RK_Step_b];
        RK_Step = tmp(:)';
        RK_Step{end+1} = 'calc_time_2_args';
        time_integrator_num_fields = 2;

    case 'ToulorgeDesmet2N84F'

        RK_Step_a = {'ToulorgeDesmet2N84F_1a', 'ToulorgeDesmet2N84F_2a', ...
                     'ToulorgeDesmet2N84F_3a', 'ToulorgeDesmet2N84F_4a', ...
                     'ToulorgeDesmet2N84F_5a', 'ToulorgeDesmet2N84F_6a', ...
                     'ToulorgeDesmet2N84F_7a', 'ToulorgeDesmet2N84F_8a'};
        RK_Step_b = {'ToulorgeDesmet2N84F_1b', 'ToulorgeDesmet2N84F_2b', ...
                     'ToulorgeDesmet2N84F_3b', 'ToulorgeDesmet2N84F_4b', ...
                     'ToulorgeDesmet2N84F_5b', 'ToulorgeDesmet2N84F_6b', ...
                     'ToulorgeDesmet2N84F_7b', 'ToulorgeDesmet2N84F_8b'};
        if strcmp(I_IEq('hall_term'), 'USE_HALL')
            calc_curlB_rho = repmat({'calc_curlB_rho_1_2args'}, size(RK_Step_a));
        else
            calc_curlB_rho = [];
        end
        if I_RunOps('variable_u')
            calc_u = repmat({'calc_u_2_args'}, size(RK_Step_a));
        else
            calc_u = [];
        end
        tmp = [calc_curlB_rho; calc_u; RK_Step_a; RK_Step_b];
        RK_Step = tmp(:)';
        RK_Step{end+1} = 'calc_time_2_args';
        time_integrator_num_fields = 2;

      otherwise

        error('Unkown time integrator %s!', I_TI('time_integrator'))

end


switch time_integrator_num_fields
    case 2
        IEq_fields = struct('field_b', field_b, ...
                            'field_b2', field_b2, ...
                            'field_curlB_rho', field_curlB_rho, ...
                            'field_u', field_u, ...
                            'field_rho', field_rho, ...
                            'current_time', current_time, ...
                            'field_divB', field_divB, ...
                            'norm2_output', norm2_output);

    case 3
        IEq_fields = struct('field_b', field_b, ...
                            'field_b2', field_b2, ...
                            'field_b3', field_b3, ...
                            'field_curlB_rho', field_curlB_rho, ...
                            'field_u', field_u, ...
                            'field_rho', field_rho, ...
                            'current_time', current_time, ...
                            'field_divB', field_divB, ...
                            'norm2_output', norm2_output);

    otherwise
        error('Wrong value of time_integrator_num_fields = %d', time_integrator_num_fields)
end

end
