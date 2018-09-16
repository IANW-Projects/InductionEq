%This project is licensed under the terms of the Creative Commons CC BY-NC-ND 3.0 license.

function [field_cs] = plot_2D(field, CS, NODES_X, NODES_Y, NODES_Z, plot_title)

% Reshapes a one dimensional field array into a 3 dimensional matrix with 
% length (NODES_X, NODES_Y, NODES_Z and generates a 2D color plot for a 
% chosen cross-section (CS). 

map = 'colormaps.plasma';

numNodes = NODES_X*NODES_Y*NODES_Z;

%Calc magnitude of specified field
for i=1:numNodes
    mag(i) = norm(field(:,i));
end
%Reshape 1D array to 3D array for plotting
field_cs = permute(reshape(mag, NODES_X, NODES_Y, NODES_Z),[2,1,3]);

FigX = figure;

%Choose cross section for plotting
if CS == 'X'
    plot_title = strcat(plot_title, ' (X-Cross-Section)');
    v(:,:) = field_cs(floor(NODES_X)/2,:,:);
    x_range = 1:NODES_Z;
    y_range = 1:NODES_Y;
elseif CS == 'Y'
    plot_title = strcat(plot_title, ' (Y-Cross-Section)');
    v(:,:) = field_cs(:,floor(NODES_Y)/2,:);
    x_range = 1:NODES_Z;
    y_range = 1:NODES_X;
else
    plot_title = strcat(plot_title, ' (Z-Cross-Section)');
    v(:,:) = field_cs(:,:,floor(NODES_Z)/2);
    x_range = 1:NODES_Y;
    y_range = 1:NODES_X;
end

hX = pcolor(x_range,y_range,v);

set(hX, 'EdgeColor', 'none');
c = colorbar;
colormap(map)
xlabel('Nodes')
ylabel('Nodes')
title(plot_title)
pbaspect([1 1 1])

shading interp;

end
