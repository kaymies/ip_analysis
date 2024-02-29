function fignum = drawStiffnessEllipseAsym(K, fignum, option, matrixID, lineSpecs, bool_arrow, num_color)

% Given the stiffness and / or damping matrix, draws an ellipse or
% epitrochoid (peanut) to qualitatively represent the stiffness directions
% and magnitude
%
% Inputs:
% K - matrix
%   - K, B, or [K, B]
% fignum - (optional) figure number on which to plot ellipse
% option - (optional) plot options:
%          1) plot ellipse using eigenvalues for magnitudes [default]
%          2) plot ellipse normalized to the principal eigenvalue
%          3) plot peanut plots
%          4) plot peanut and ellipse
% matrixID - (optional) 'K' or 'B', to identify whether the matrix is a 
%             stiffness or damping matrix
% lineSpecs - (optional) struct specifying the following line parameters:
%              - color: color of plot lines
%              - width: width of plot lines
% bool_arrow - true: draw arrow to indicate curl direction of asymmetry
%            - false: don't draw arrow

%NOTE: for the optional inputs, pass an empty matrix to leave them
%unspecified (see example below)
%
%---- Example usage ------------------------------------------------------- 
%     K = [2 -1; 1 5];
%     drawStiffnessEllipse(K);
%     option = 4;
%     lineSpecs.width = 2;
%     lineSpecs.color = [0.8,0,0];
%     drawStiffnessEllipse(K,2,option,[],lineSpecs);
% -------------------------------------------------------------------------
%
% Rika Sugimoto Dimitrova
% 2021-03-14
%
% Modified
% 2021-06-22 - included plotting of ellise
%              added axis labels
%              added option to pass [K, B] together
% 2021-06-29 - corrected ellipse minor axis calculations
%              added check for complex eigenvalues
%              added lineSpecs as input
% 2021-08-23 - corrected so that if D(2,2) > D(1,1) that will be the major
%              axis
% 2021-08-27 - plotting symmetric component as ellipse,
%              antisymmetric as circle
% 2021-11-27 - adding peanut plots
% 2022-02-28 - text interpreter to latex
% 2022-04-12 - include asymmetry arrow


if nargin > 1 && ~isempty(fignum)
    figH = figure(fignum);hold on;
else
    figH = figure;hold on;
    fignum = figH.Number;
end

if nargin < 3 || isempty(option)
    option = 1;
end

color_alpha = 0.8;
linecolor_alpha = color_alpha*ones(1,3);
% if nargin < 5 || isempty(lineSpecs)
%     linecolor = zeros(1,3);
%     linewidth = 1.2;
% else
%     linecolor = lineSpecs.color;
%     linewidth = lineSpecs.width;
% end

if nargin < 6 || isempty(bool_arrow)
    bool_arrow = false;
end


matID{1} = [];
if nargin > 3 && ~isempty(matrixID)
    matID{1} = matrixID;
end
    
K_x{1} = K(1:2,1:2);
if length(K) ~= min(size(K))
    K_x{2} = K(1:2,3:4);
    matID{1} = 'K'; matID{2} = 'B';    
end
N_plots = length(K_x);

if N_plots == 2
    set(figH,'Units','inches');    
    set(figH,'Position',[0,4.7,8.5,3.3]);
    %set(figH,'Units','normalized');    
    %set(figH,'Position',[0 0 0.62 0.42]);
end

for i = 1:N_plots
    subplot(1,N_plots,i); hold on;
    if num_color ~isempty(num_color)
        mycolors = parula(num_color);
        set(gca, 'ColorOrder', mycolors);
    end
    K_symmetric = (K_x{i}+K_x{i}')/2;
    K_asymmetric = (K_x{i}-K_x{i}')/2;
    [V_sym,D_sym] = eig(K_symmetric);
    D_asym = K_asymmetric(1,2);

    disp(D_sym);
    if ~isreal(D_sym)
        warning('Complex eigenvalues');
        V_sym = real(V_sym);
        D_sym = real(D_sym);
    end

    [eigval_max,max_index] = max([D_sym(1,1),D_sym(2,2)]);
    [~,min_index] = min([D_sym(1,1),D_sym(2,2)]);
    axis1_x = [(-1)*V_sym(1,max_index)*D_sym(max_index,max_index), V_sym(1,max_index)*D_sym(max_index,max_index)];
    axis1_y = [(-1)*V_sym(2,max_index)*D_sym(max_index,max_index), V_sym(2,max_index)*D_sym(max_index,max_index)];            
    axis2_x = [(-1)*V_sym(1,min_index)*D_sym(min_index,min_index), V_sym(1,min_index)*D_sym(min_index,min_index)];
    axis2_y = [(-1)*V_sym(2,min_index)*D_sym(min_index,min_index), V_sym(2,min_index)*D_sym(min_index,min_index)];
    switch option
        case 1
            if ~isempty(matID{i})
               switch matID{i}
                   case 'K'
                       title_str = 'Stiffness (Nm/rad)';
                       xlabel_str = '$$q_{1}$$';
                       ylabel_str = '$$q_{2}$$';
                   case 'B'
                       title_str = 'Damping (Nms/rad)';
                       xlabel_str = '$$\dot{q}_{1}$$';
                       ylabel_str = '$$\dot{q}_{2}$$';
               end 
            end

        case 2
            axis1_x = axis1_x/eigval_max;
            axis1_y = axis1_y/eigval_max;
            axis2_x = axis2_x/eigval_max;
            axis2_y = axis2_y/eigval_max;
            
            D_asym = D_asym/eigval_max;

            if ~isempty(matID{i})
               switch matID{i}
                   case 'K'
                       title_str = 'Stiffness';
                       xlabel_str = '$$q_{1}$$';
                       ylabel_str = '$$q_{2}$$';
                   case 'B'
                       title_str = 'Damping';
                       xlabel_str = '$$\dot{q}_{1}$$';
                       ylabel_str = '$$\dot{q}_{2}$$';
               end 
            end      
        case 3 
            if ~isempty(matID{i})
               switch matID{i}
                   case 'K'
                       title_str = 'Stiffness (Nm/rad)';
                       xlabel_str = '$$q_{1}$$';
                       ylabel_str = '$$q_{2}$$';
                   case 'B'
                       title_str = 'Damping (Nms/rad)';
                       xlabel_str = '$$\dot{q}_{1}$$';
                       ylabel_str = '$$\dot{q}_{2}$$';
               end 
            end
        case 4 
            if ~isempty(matID{i})
               switch matID{i}
                   case 'K'
                       title_str = 'Stiffness';
                       xlabel_str = '$$q_{1}$$';
                       ylabel_str = '$$q_{2}$$';
                   case 'B'
                       title_str = 'Damping';
                       xlabel_str = '$$\dot{q}_{1}$$';
                       ylabel_str = '$$\dot{q}_{2}$$';
               end 
            end
    end

    %%% Uncomment to plot major and minor axes %%%
    %plot(axis1_x,axis1_y,'Color',linecolor);
    %plot(axis2_x,axis2_y,'Color',linecolor);

    % Plot ellipse
    a = norm([axis1_x(2),axis1_y(2)]);
    b = norm([axis2_x(2),axis2_y(2)]);
    
    th = 0:0.01:(2*pi);
    
    %%% For "evenly spaced" (by angle) points long the ellipse %%%
    %     r_sym = sqrt(1 ./ ((cos(th)/a).^2 + (sin(th)/b).^2));
    %     y_sym = r_sym.*cos(th);
    %     y_sym = r_sym.*sin(th);

    x_sym = a.*cos(th);
    y_sym = b.*sin(th);
    
    r_asym = abs(D_asym(1,1));
    x_asym = r_asym.*cos(th);
    y_asym = r_asym.*sin(th);
    arrow_x = r_asym.*cos(pi/4);
    arrow_y = r_asym.*sin(pi/4);
    
    x_peanut = (a.*cos(th).^2 + b.*sin(th).^2).*cos(th);
    y_peanut = (a.*cos(th).^2 + b.*sin(th).^2).*sin(th);
    
    th_rot = atan2(axis1_y(2), axis1_x(2));
    R = [cos(th_rot) -sin(th_rot); sin(th_rot) cos(th_rot)];
    
    ellipse_points = R*[x_sym;y_sym];
    peanut_points = R*[x_peanut;y_peanut];
    switch option
        case 3
            plot(peanut_points(1,:),peanut_points(2,:),...
                'Color',linecolor,'Linewidth',linewidth);
        case 4
            plot(ellipse_points(1,:),ellipse_points(2,:),...
                'Color',linecolor_alpha,'Linewidth',linewidth);
            plot(peanut_points(1,:),peanut_points(2,:),...
                'Color',linecolor,'Linewidth',linewidth);
        otherwise
%             plot(ellipse_points(1,:),ellipse_points(2,:),...
%                 'Color',linecolor,'Linewidth',linewidth);
    end
            
    plot(x_asym,y_asym);%,...
        %'Color',linecolor,'Linewidth',linewidth);

    if option == 2
        xlim([-1,1]);
        ylim([-1,1]);
    end
    axis equal
    
    xticks(yticks);

    if ~isempty(matID{i})
        xlabel(xlabel_str,'Interpreter','latex');
        ylabel(ylabel_str,'Interpreter','latex');  
        title(title_str,'Interpreter','latex');
    end
    
    if bool_arrow
    % Get normalized units of arrow position (https://www.mathworks.com/matlabcentral/answers/166823-converting-plot-coordinates-to-normalized-coordinates)
    subaxis = subplot(1,2,i);
    AX=axis(gca); 
    Xrange=AX(2)-AX(1);
    Yrange=AX(4)-AX(3);       
    X = (arrow_x-AX(1))/Xrange*(subaxis.Position(3)) + subaxis.Position(1);
    Y = (arrow_y-AX(3))/Yrange*(subaxis.Position(4)) + subaxis.Position(2);
    annotation('arrow', ...
        [(X+0.005*subaxis.Position(3)) (X)], [(Y+0.01*subaxis.Position(4)) (Y+0.02*subaxis.Position(4))]);
    end
    
end

end