%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 4.0 license.

function [settings] = generate_settings(map, keys)

    % Atomatically generates a formatted settings string that contains all
    % compiler optimizations and OpenCL defines of map with the given keys where 
    % keys is a cell array, e.g. settings_mesh = generate_settings(I_Mesh, {'DX'; 'DY'; 'DZ'}

    settings = '';
    for i=1:length(keys)       
        if isKey(map, keys{i})
            if ischar(map(keys{i})) 
                if strcmp(keys{i}, 'optimizations')
                    settings = strcat(settings, map(keys{i}));
                elseif contains(map(keys{i}),'USE')
                    settings = strcat(settings, sprintf(' -D%s=1', map(keys{i})));
                else
                    settings = strcat(settings, sprintf(' -D%s=%s', keys{i}, map(keys{i})));
                end
            else
                if isinteger(map(keys{i}))
                    settings = strcat(settings, sprintf(' -D%s=%d', keys{i}, map(keys{i})));
                else
                    settings = strcat(settings, sprintf(' -D%s=%.16e%', keys{i}, map(keys{i})));
                end
            end
        else
            fprintf('Unknown identifier %s\n', keys{i})
        end    
    end
    
    
end
    
