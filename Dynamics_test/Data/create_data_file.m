%% Saving Human Data as Struct
setpath
filename = 'duarte_young';
folder = fullfile('Data','Human');
file = sprintf('%s.mat',filename);
folder = fullfile(folder,file);

sub_data = young_data;

field1 = 'SubjectType';  value1 = 'Old Duarte';
field2 = 'NumSubjects';  value2 = height(sub_data)/3;
field8 = 'MeanAge'; value8 = mean(sub_data(:,2));
% field9 = 'MeanHeight_m'; value9 = 1.67;
% field10 = 'MeanMass_kg'; value10 = 87.5;
field3 = 'Frequency';  value3 = Frequency;
field4 = 'IPDataAverage';  value4 = mean(sub_data(:,9:end),'omitnan');
field5 = 'Plane'; value5 = 'Sagittal';
% field6 = 'WhichLeg'; value6 = 'Paretic';
field7 = 'StandardDeviation'; value7 = std(sub_data(:,9:end),'omitnan');
field14 = 'TrialNumber'; value14 = sub_data(:,1);
field15 = 'Age'; value15 = sub_data(:,2);
field16 = 'Parameters'; value16 = sub_data(:,3:8);
field17 = 'IP'; value17 = sub_data(:,9:end);
field18 = 'IPAveSubject'; 
B = reshape(value17,3,[]);
C = mean(B,1,'omitnan');
value18 = reshape(C,[], size(value17,2));

% field11 = 'ForceAP'; value11 = F_ap_N_nonparetic;
% field12 = 'ForceVert'; value12 = F_vert_N_nonparetic;
% field13 = 'CoPAP'; value13 = CoP_ap_m_nonparetic;
duarte_young = struct(field1,value1,field5,value5,field2,value2,field14,value14,...
    field8,value8,field15,value15,...
    field3,value3,field17,value17,field18,value18,field4,value4,...
    field7,value7,field16,value16);
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