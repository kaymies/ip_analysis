young = load("duarte_young.mat");
youngfit = load('duarte_old_fitted.mat');
numSubplots = 38;
gridRows = 5;
gridCols = 8;
% Create a figure
figure;
% Iterate over each subplot
for i = 1:numSubplots
% Select the current subplot
subplot(gridRows, gridCols, i);
% Plot the data in the current subplot
plot(youngfit.Frequency, squeeze(youngfit.IP(i,:,:)),'r-');
hold on
plot(young.Frequency,squeeze(young.IP(i,:,:)),'b-')
% Add labels and title to the subplot
xlabel('Frequency [Hz]')
ylabel('zIP/zCoM')
title(['Old Sub', num2str(i)]);
end
legend('Fitted (Gruben)','','','Raw (Duarte)')

% IP_reshape = zeros(38,3,29);
% for col = 1:38
%     % Extract every three rows from the original matrix
%     extractedRows = IP((col-1)*3+1:col*3, :);
%     
%     % Assign the extracted rows as a column in the reshaped matrix
%     IP_reshape(col, :, :) = extractedRows;
% end