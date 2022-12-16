function [] = generate_plot(human,sim)
%% Load Data
% Human Data
human_struct = human;
% Best-Fit Simulation Data
sim_struct = sim;
%% Plot human vs. simulation
figure();
plot(human_struct.Frequency,human_struct.IPDataAverage)
hold on
plot(sim_struct.Frequency,sim_struct.DataWithoutOutliers)
legend('human', 'simdata')
xlabel('Frequency [Hz]')
ylabel('IP (Fraction of CoM)')
xlim([0.5 5.3])
ylim([0 2.5])
end