function [p_mean, p_SD] = ...
    plotStdArea(x,y_mean,y_std,ifig,option,lineSpecs,min_lb,bool_pow2dB)
% Plots a shaded area showing 1 STD from the mean.
% 
% Inputs:
% x      - independent variable
% y_mean - dependent variable, mean value
% y_std  - dependent variable, standard deviation
% ifig   - figure number to plot on; if empty, create a new figure
% option - plotting options:
%          0: plot mean+/-std as a shaded area only (default)
%          1: plot mean as a line only
%          2: plot mean+/-std area and mean line
% lineSpecs - line specifications
%           color: 1x3 vector for color
%           width: line width of mean
%           alpha: transparency of SD area
% min_lb - minimum lower bound;
%          if we are plotting log-scale plots, we would want min_lb > 0
% bool_pow2dB - true: convert y_mean and lb,ub to dB via 10*log10()
%             - false: leave y_mean and y_std as is.
%
% 2021-03-19
% Rika Sugimoto Dimitrova
% 
% Modified:
% 2021-04-12 - added plotting options
% 2022-03-02 - added lineSpecs

p_mean = [];
p_SD = [];

if nargin > 3 && ~isempty(ifig)
    figure(ifig);hold on;
else
    figure;hold on;
end

if nargin < 5 || isempty(option)
    if nargin > 2 && ~isempty(y_std)
        option = 0;
    else
        option = 1;
    end
end

if nargin < 6 || isempty(lineSpecs)
    lineSpecs.color = [0 0 0];
    lineSpecs.width = 1.2;
    lineSpecs.alpha = 0.5;
    lineSpecs.style = '-';
end

colors = get(gca,'colororder');

if option == 0 || option == 2
    
    if size(x,1) > size(x,2)
        x = x';
        y_mean = y_mean';
        y_std = y_std';
    end
    
    if size(y_std,1) == 2
        if y_std(1,1) < y_std(2,1)
            lb = y_std(1,:);
            ub = y_std(2,:);
        else
            lb = y_std(2,:);
            ub = y_std(1,:);   
        end
    else   
        lb = y_mean-y_std;
        ub = y_mean+y_std;
    end
    
    if nargin > 7 && min_lb > -inf && ~isempty(bool_pow2dB)
        lb = max( y_mean-y_std , abs(y_mean)*min_lb ); % confine to positive values if plotting CPSD
        if bool_pow2dB
            y_mean = 10*log10(y_mean);
            lb = 10*log10(lb);
            ub = 10*log10(ub);
        end
    end
    
    for i = 1:size(lb,1)
        if nargin > 6 && min_lb > 0
            p_SD = plot(x,lb,'-','Linewidth',0.5*lineSpecs.width,...
                'Color',lineSpecs.color);
            p_SD = plot(x,ub,'-','Linewidth',0.5*lineSpecs.width,...
                'Color',lineSpecs.color);
            %p_SD = fill([x,flip(x)],[ub(i,:),flip(lb(i,:))],lineSpecs.color);
        else
            p_SD = fill([x,flip(x)],[ub(i,:),flip(lb(i,:))],lineSpecs.color,...
                'FaceAlpha',lineSpecs.alpha,'EdgeAlpha',0);
        end
        %fill([x,flip(x)],[ub(i,:),flip(lb(i,:))],colors(i,:),...
        %    'FaceAlpha',lineSpecs.alpha,'EdgeAlpha',0);
    end
    
    if option == 2
        ax = gca;
        colors = get(gca,'colororder');
        %ax.ColorOrder(1:2*size(y_mean,1),:) = ...
        %    [colors(1:size(y_mean,1),:); colors(1:size(y_mean,1),:)];
        p_mean = plot(x,y_mean,lineSpecs.style,'Linewidth',lineSpecs.width,...
            'Color',lineSpecs.color);
    end
    
else
    p_mean = plot(x,y_mean,'Linewidth',lineSpecs.width);
end

end