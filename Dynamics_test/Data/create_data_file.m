%% Saving Human Data as Struct
setpath
filename = 'duarte_young_paper';
folder = fullfile('Data','Human');
file = sprintf('%s.mat',filename);
folder = fullfile(folder,file);

% sub_data = old_ave;

field1 = 'SubjectType';  value1 = 'Duarte Young Not Fitted from Paper';
field2 = 'NumSubjects';  value2 = 70.3;
% field8 = 'MeanAge'; value8 = mean(sub_data(:,2));
field9 = 'MeanHeight_m'; value9 = 1.71;
field10 = 'MeanMass_kg'; value10 = 68.5;
field3 = 'Frequency';  value3 = freq_Hz_bpf;
field4 = 'IPDataAverage';  value4 = mean(zIP_young_bpf);
field5 = 'Plane'; value5 = 'sgt';
% field6 = 'WhichLeg'; value6 = 'Paretic';
% field7 = 'StandardDeviation'; value7 = std_frt_nobeam;
% field14 = 'SubNumber'; value14 = num_sub_frt_nobeam;
% field15 = 'Age'; value15 = sub_data(:,2);
% field16 = 'Parameters'; value16 = sub_data(:,3:8);
field17 = 'IP'; value17 = zIP_young_bpf;
% field18 = 'IPAveSubject'; 
% B = reshape(value17,3,[]);
% C = mean(B,1,'omitnan');
% value18 = reshape(C,[], size(value17,2));

% field11 = 'ForceAP'; value11 = F_ap_N_nonparetic;
% field12 = 'ForceVert'; value12 = F_vert_N_nonparetic;
% field13 = 'CoPAP'; value13 = CoP_ap_m_nonparetic;
field19 = 'Pose'; value19 = 'pose_I';
duarte_young = struct(field1,value1,field5,value5,field2,value2,field9,value9,...
    field10,value10,field3,value3,field4,value4,...
    field17,value17,field19,value19);
save(folder, '-struct', filename);
%% Updating Data File
filename = 'duarte_young';
duarte_young = load(sprintf('%s.mat',filename));

duarte_young.Frequency_bpf = Frequency;
duarte_young.IP_bpf = IP;
% duarte_old.Frequency_cpsd = freq_Hz_cpsd(6:82);
% % duarte_old_fitted.IP = IP_reshape;
duarte_young.IPAveSubject_bpf = IPAveSubject;
duarte_young.IPDataAverage_bpf = IPDataAverage;
% duarte_old_fitted.Frequency = freq_Hz_bpf;
% duarte_young.MeanHeight_m = MeanHeight_m;
% duarte_young.Plane = 'sgt';
% duarte_young.Pose = 'pose_I';
% duarte_young.LowFreqRange = 5:11;
% duarte_young.HighFreqRange = 14:29;

folder = fullfile('Data','Human');
file = sprintf('%s.mat',filename);
folder = fullfile(folder,file);
save(folder, '-struct', filename);
%% Extracting data from Duarte dataset
old_data = [];
young_data = [];
for i = 1:height(data)
    if data(i,2) > 60
        old_new = data(i,:);
        old_data = [old_data; old_new];
    else
        young = data(i,:);
        young_new = data(i,:);
        young_data = [young_data; young_new];
    end
end