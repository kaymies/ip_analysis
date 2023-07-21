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

duarte_young.DataInfo = young_data_n;
% duarte_old.MeanHeight_m = MeanHeight_m;
% duarte_old.Plane = 'sgt';
% duarte_old.Pose = 'pose_I';
% duarte_old.LowFreqRange = 5:11;
% duarte_old.HighFreqRange = 14:29;

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
%% Extracting data info by age
count_old = 1;
count_young = 1;
for i = 1:height(datainfo_bySubj_GrubenSet)
    if datainfo_bySubj_GrubenSet.Age(i) > 60
        datainfo_bySubj_GrubenSet.Age(i)
        old_data(count_old,:) = datainfo_bySubj_GrubenSet(i,:);
        count_old = count_old + 1;
    else
        young_data(count_young,:) = datainfo_bySubj_GrubenSet(i,:);
        count_young = count_young + 1;
    end
end
%% Take every 3 rows
i = 1;
for col = 1:3:195
    young_data_n(i,:) = young_data(col, :);
    i = i+1;
end