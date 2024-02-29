function [f_IP, IP_ratio] = getIPvsFreq(data,params)
% Obtain IP at different frequencies
%
% Inputs:
% data = struct containing the following info
%      - COM = 2xN COM position in the 2D plane of interest
%      - COP = 2xN COP position in the 2D plane of interest
%      - ForcePlateDynamics = 2xN ground reaction forces
%      - FreqSamp = sampling frequency of data, in Hz
% params = struct containing parameters for IP analysis
%      - method = 'bpf','cpsd' % band-pass-filter method vs. co-spectral
%      density method
%
% Rika Sugimoto Dimitrova
% 2020-12-13
% Modified:
% 2021-05-27 - cleaning up IP calc code (deleting unused plotting code)
% 2023-04-04 - updated to use cpsd
% 2023-05-31 - adding params.method as input
% 2023-06-27 - updated to take 2D data only (select plane outside function)

window_size = 2^9;
nfft = window_size*2;
noverlap = 0.5*window_size;

iX = 1; iZ = 2;

Fs_Hz = data.FreqSamp;
dt = 1/Fs_Hz;

COP = data.COP(iX,:); % 1D CoP
Fxy = data.ForcePlateDynamics(iX,:); % horizontal ground reaction force
COM_z = data.COM(iZ,:); % vertical CoM
Fz = data.ForcePlateDynamics(iZ,:);       % vertical ground reaction force

theta = -Fxy./Fz;

switch params.method
    case 'cpsd'
        %y = [detrend(COP,0)' detrend(theta,0)'];
        y = [detrend(COP,1)' detrend(theta,1)']; % 2023-04-06 linear detrend
                                                  % better for some cases 
                                                  % (e.g. simulation of integrated white noise)
        [Gyy, f_IP] = cpsd(y,y,hamming(window_size),noverlap,nfft,'mimo',Fs_Hz);
        
        IP_ratio = zeros(window_size+1,1);
        for i = 1:(window_size+1)
            [V,D] = eig(squeeze(real(Gyy(i,:,:))));
            [~,i_max] = max(diag(abs(D)));
            IP_ratio(i) = V(1,i_max)/V(2,i_max) / mean(COM_z);
        end

    case 'bpf'

        %% Kaymie's code
        % Obtain IP at different frequencies (inspired by Marta's code)
        f_i = 0.4; f_int = 0.2; f_range = 7.4; % starting freq, interval between tested frequencies, frequency range
        f = f_i:f_int:f_i+f_range; % frequencies
        IP = zeros(1,length(f)); % IP
        theta = -Fxy./Fz;
        %theta = atan2(Fz,Fxy); % Basically gives the same result as above after subtracting the mean (Rika: 2021-05-27)
        theta = theta - mean(theta);
        cop_ham = detrend(COP,0)'.*hamming(length(COP)); % apply Hamming window to CoP data
        theta_ham = theta'.*hamming(length(theta)); % apply Hamming window to the force angle
        for n = 1:length(f) % for all frequencies designated above:
            [B,A] = butter(2,[f(n) f(n)+f_int]/(1/dt/2)); % define the Butterworth filter
            COP_f = filtfilt(B,A,cop_ham); % apply filter to CoP data
            theta_f = filtfilt(B,A,theta_ham); % apply filter to force angle data
            [coeff,~,~,~,~,~] = pca([COP_f theta_f]); % PCA
            IP(n) = coeff(1,1)/coeff(2,1); % compute IP and store
        end
        IP_ratio = IP/mean(COM_z);
        f_IP = f+0.1;
end

% figure;plot(f+0.1,IP_ratio1);
% hold on;plot(f_IP,IP_ratio)
% xlim([0,8]);
% ylim([-0.5,2.5]);
% legend('Band-pass','CPSD');

end