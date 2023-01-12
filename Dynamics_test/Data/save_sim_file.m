function [sim_data] = save_sim_file(input,filename,data,meansim,numtrial)
folder_name = sprintf('%s_%s',input.Controller.type,input.CoordinateFrame);
folder = fullfile('Data','Simulation',folder_name);
addpath(folder);
file = sprintf('%s.mat',filename);
folder = fullfile(folder,file);
%Save data as struct
field1 = 'Alpha';  value1 = input.Controller.alpha;
field2 = 'Beta';  value2 = input.Controller.beta;
field3 = 'NoiseRatio';  value3 = input.NoiseRatio;
field4 = 'CoordinateFrame';  value4 = input.CoordinateFrame;
field5 = 'Plane'; value5 = input.plane;
field6 = 'Controller'; value6 = input.Controller.type;
field7 = 'NumTrials'; value7 = numtrial;
field8 = 'SamplingFrequency'; value8 = input.FreqSampKin;
field9 = 'TrialDuration'; value9 = input.trialDuration;
field10 = 'Mass'; value10 = input.TotalMass;
field11 = 'Height'; value11 = input.TotalHeight;
field12 = 'Frequency'; value12 = input.Frequency;
field13 = 'FrequencyWindow'; value13 = input.FrequencyWindow;
field14 = 'Data'; value14 = data;
field15 = 'DataWithoutOutliers'; value15 = meansim;
sim_data = struct(field1,value1,field2,value2,field3,value3,field4,value4,...
    field5,value5,field6,value6,field7,value7,field8,value8,field9,value9,...
    field10,value10,field11,value11,field12,value12,field13,value13,...
    field14,value14,field15,value15);
save(folder, '-struct', 'sim_data');
end

