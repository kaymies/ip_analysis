function postproc_animation(t_, x_, model_param, model_type, saveanim, anim_title)
% http://stackoverflow.com/questions/18776172/in-matlab-how-do-i-change-the-arrow-head-style-in-quiver-plot
% http://stackoverflow.com/questions/25729784/how-to-draw-an-arrow-in-matlab

switch model_type
    case 'SIP'
        n_q = 1;
    case 'DIP'
        n_q = 2;
end
%% 
t_window = 2;
isTerminate = 0;
di_anim = 4;
di_anim_default  = di_anim;
%% Initialize the animation figure and axes
set(0,'Units','pixels')
screenSize = get(0, 'ScreenSize');
figWidth = 600; figHeight = 600; 
X0 = (screenSize(3)-figWidth)/2;
Y0 = (screenSize(4)-figHeight)/2;

% Set figure
fig = figure(100);
set(gcf, 'Position', [X0,Y0, figWidth, figHeight]);
set(gcf,'KeyPressFcn',@printfig);
set(gcf, 'color','w'); % none: terrible result
% set(gcf, 'menubar','none')

% Set axis
animWidth  = 1; animHeight  = 1;
animX = 0;  animY = 0;
axes_anim = subplot('position', [animX animY animWidth animHeight]);
grid on;
set(axes_anim, 'xtick', []); set(axes_anim,'xticklabel',[]); 

%% set up movie
if saveanim
    writerObj = VideoWriter([anim_title,'MPEG-4']);
    open(writerObj)
end
%% Animation
i_anim= 1;
n_data = length(t_);
q_ = x_(1:n_q,:);
while (1)
    %% Draw
    draw_anim(t_(i_anim),q_(:,i_anim), i_anim, axes_anim, model_param );

    % video instruction
    getxlim = get(axes_anim, 'xlim');
    getylim = get(axes_anim, 'ylim');
    xpos = getxlim(1);    ypos = getylim(2)*0.95;
    txt21 = 'Press F: forward, R: reward, S: stop, E: terminate, 1-5: Control Speed';
    hText = text(xpos,ypos,txt21); 
    drawnow;    
    
    % settings    
    if saveanim
        frame = getframe(fig);
        writeVideo(writerObj,frame);
    end
    
    if abs(di_anim) < 1e-4
        % stop
    else
        i_anim = i_anim + di_anim;
    end
    
    if i_anim >= n_data
        i_anim = n_data;
        isTerminate = 1;
    elseif i_anim < 1
        i_anim = 1;
    end
    if isTerminate
        if saveanim
            close(writerObj);
        end
        break;
    end
end

function draw_anim(t, q, i_anim, axes, param)
    % cleaning before drawing
    tem = get(axes, 'Children');
    delete(tem);
    hold (axes, 'on');
    
    % draw ground
    hGround = line([-2 2], [0 0], 'linewidth',2,'color','k');
    set(hGround, 'parent', axes);
    
    % draw link    
    [p_O_CF,p_O_CH,p_O_CK,p_O_NE,p_O_HE,p_O_RS,p_O_RE,p_O_RW,p_O_RH,p_O_LS,p_O_LE,p_O_LW,p_O_LH] = auto_JointPosition(q,param);
    [p_O_CoM,~,~] = auto_COMinformation(q,zeros(size(q)),param);
    
    
    draw_segment(axes, p_O_CF, p_O_CK, 'b', 10);
    draw_segment(axes, p_O_CK, p_O_CH, 'b', 10);
    draw_segment(axes, p_O_CH, p_O_HE, 'r', 10);
    
    draw_segment(axes, p_O_RS, p_O_RH, 'g', 5);
    draw_segment(axes, p_O_LS, p_O_LH, 'g', 5);
    
    
    % draw joints
    draw_joint(axes, p_O_CF, 'k', 0.05/2);
    draw_joint(axes, p_O_CH, 'k', 0.05/2);
    draw_joint(axes, p_O_CoM, 'k', 0.05);
    
    % draw texts
    txt1 = 'Center of mass';
    y_txt1 = 1; 
    x_txt1 = 0.2;
    hText1 = text(x_txt1,y_txt1,txt1); 
    hText1.FontSize = 15; 
    hText1.VerticalAlignment = 'top';
    hText1.HorizontalAlignment = 'left';
    set(hText1, 'parent', axes);
    qx = 0.2; qy = 0.97;
%     hQuiver = quiver(axes, qx, qy, (-qx+centerX)*.7, (-qy+centerY)*.7,0,'k');
    
    % indicate time
    txt2 = ['t= ', num2str(t),' [s]'];
    y_txt2 = 1.5; 
    x_txt2 = -0.5;
    hText2 = text(x_txt2,y_txt2,txt2); 
    set(hText2, 'parent', axes);
    
    % set axis
    ylim(axes,[-0.1 2-.1]);
    xlim(axes, [-1 1]) ;   
    xlabel(axes,'distance [m]','fontsize',12); 
    ylabel(axes,'height [m]','fontsize',12);
end
    function handle = draw_segment(axes, p1, p2, color, width)
       handle = plot(axes, [p1(1) p2(1)], [p1(2) p2(2)],color,'linewidth',width);
    end
    function handle = draw_joint(axes, p, color, radius) 
        J1X = p(1);
        J1Y = p(2);        
        handle  = rectangle('Position',[J1X - radius, J1Y - radius, radius*2, radius*2],...
        'Curvature',[1,1],'FaceColor',color);
        set(handle,'parent',axes);
    end
% function draw_phase(t,q, Dq, i_anim, axes)
%     tem = get(axes, 'Children');
%     delete(tem);
%             
%     plot(axes,rad2deg(q), rad2deg(Dq),'r'); hold on
%     i_del = find (t>t_window , 1);
%     if i_anim > i_del
%         plot(axes,rad2deg(q(i_anim-i_del:i_anim,:)), rad2deg(Dq(i_anim-i_del:i_anim,:)),'k:','linewidth',1.5); hold (axes,'on');
%     else
%         plot(axes,rad2deg(q(1:i_anim,:)), rad2deg(Dq(1:i_anim,:)),'k:','linewidth',1.5); hold (axes,'on');
%     end
%     xlabel(axes, '$\theta$ [deg]','fontsize',12); 
%     ylabel(axes,'$\dot{\theta}$ [deg/s]','fontsize',12)
%     grid (axes,'on');
% end

function printfig(src,evnt)
    if evnt.Character == '1'
    di_anim = di_anim_default/4;
    end
    if evnt.Character == '2'
    di_anim = di_anim_default/2;
    end
    if evnt.Character == '3'
    di_anim = di_anim_default;
    end
    if evnt.Character == '4'
    di_anim = di_anim_default*2;
    end
    if evnt.Character == '5'
    di_anim = di_anim_default*4;
    end
    if evnt.Character == 'r'
    di_anim = -1*di_anim;
    di_anim_default = -di_anim_default;
    end
    if evnt.Character == 's'
    di_anim = eps;
    end
    if evnt.Character == 'e'
    isTerminate = 1;
    end 
end 

end