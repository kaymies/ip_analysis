% All parameters are based on 
% De Leva, Paolo. "Adjustments to Zatsiorsky-Seluyanov's segment inertia parameters." Journal of biomechanics 29.9 (1996): 1223-1230.
% Shoulder width for frontal plane model parameters were obtained from 
% NASA man-systems integration standards, https://msis.jsc.nasa.gov/sections/section03.htm

% m_bar: mass ratio
% r_bar: radii of gyration rel. segment length
% l_bar: segment length ratio

% input:
% - 'gender': 'F', 'M'
% - 'plane': 'sgt', 'frt'

% note: in de Leva, 'sagittal' of radii of gyration means radii of gyration
% about sagittal axis. Therefore the moment of inertia calculated by that 
% value should be used for frontal plane model. 
% likewise, 'transverse' should be used for sagittal plane model.
%% if symbolic parameter is needed
% syms g positive
% syms m_l L_1 c_1 j_1  positive  % lower body: *two* legs (thigh+shank)
% syms m_2 L_2 c_2 j_2  positive  % upper body: trunk + head
% syms m_a L_a c_a j_a  positive  % arm: *one* arm, upper arm + forearm + hand



%% Define class
ClassUpperarm = ClassLink(275.1, 281.7, 2.55, 2.71, 57.54, 57.72, 27.8, 28.5, 26.0, 26.9);
ClassForearm = ClassLink(264.3, 268.9, 1.38, 1.62, 45.59, 45.74, 26.1, 27.6, 25.7, 26.5);
ClassHand = ClassLink(78.0, 86.2, 0.56, 0.61, 74.74, 79.00, 53.1, 62.8, 45.4, 51.3);
% d= ClassLink(78.0, 86.2, 0.56, 0.61, 74.74, 79.00, 53.1, 62.8, 45.4, 51.3);

ClassThigh = ClassLink(368.5, 422.2, 14.78, 14.16, 100-36.12, 100-40.95, 36.9, 32.9, 36.4, 32.9);

ClassHead = ClassLink(243.7, 242.9, 6.68, 6.94, 100-48.41, 100-50.02, 27.1, 30.3, 29.5, 31.5);
ClassTrunk = ClassLink(614.8, 603.3, 42.57, 43.46, 100-49.64, 100-51.38, 30.7, 32.8, 29.2, 30.6);

ClassShank = ClassLink(438.6, 440.3, 4.81, 4.33, 100-43.52, 100-43.95, 26.7, 25.1, 26.3, 24.6);
%% Gravity
g = 9.81; 


%% Lump Lower Body Segments
Thigh = ClassThigh.assignLink(totalMass,totalHeight,gender,plane);
Shank = ClassShank.assignLink(totalMass,totalHeight,gender,plane);

%% Lump HAT
Head = ClassHead.assignLink(totalMass,totalHeight,gender,plane);
Trunk = ClassTrunk.assignLink(totalMass,totalHeight,gender,plane);

%% Lump Arm
Upperarm = ClassUpperarm.assignLink(totalMass,totalHeight,gender,plane);
Forearm = ClassForearm.assignLink(totalMass,totalHeight,gender,plane);
Hand = ClassHand.assignLink(totalMass,totalHeight,gender,plane);

%% Lump Upper Body (adducted body)
l_sjc = L_SJC(totalHeight, gender, plane);
w_sjc = w_SJC(totalHeight, gender, plane);


%%
function l_sjc = L_SJC(totalHeight, gender, plane)
    switch gender
        case 'F'
            l_sjc = 497.9/1735 * totalHeight;
        case 'M'
            l_sjc = 515.5/1741 * totalHeight;
    end
end
function w_sjc = w_SJC(totalHeight, gender, plane)
% obtained from NASA standards, averaged 5th, 50th, 95th percentile data
    switch plane
        case 'sgt'
            w_sjc = 0;
        case 'frt'
            switch gender
                case 'F'
                    w_sjc = 0.23/2 * totalHeight;
                    % Japanese, 40 yo, 2000
                case 'M'
                    w_sjc = 0.22/2 * totalHeight;
                    % American, 40 yo, 2000
            end
    end
    
end



