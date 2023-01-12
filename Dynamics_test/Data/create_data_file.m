%% Saving Human Data as Struct
setpath
filename = 'nonpar';
folder = fullfile('Data','Human');
file = sprintf('%s.mat',filename);
folder = fullfile(folder,file);

field1 = 'SubjectType';  value1 = 'Paretic';
field2 = 'NumSubjects';  value2 = 9;
field8 = 'MeanAge'; value8 = 59.8;
field9 = 'MeanHeight_m'; value9 = 1.67;
field10 = 'MeanMass_kg'; value10 = 87.5;
field3 = 'Frequency';  value3 = Frequency;
field4 = 'IPDataAverage';  value4 = IPDataAverage;
field5 = 'Plane'; value5 = 'Sagittal';
field6 = 'WhichLeg'; value6 = 'Paretic';
field7 = 'StandardDeviation'; value7 = StandardDeviation;
field11 = 'ForceAP'; value11 = F_ap_N_nonparetic;
field12 = 'ForceVert'; value12 = F_vert_N_nonparetic;
field13 = 'CoPAP'; value13 = CoP_ap_m_nonparetic;
nonpar = struct(field1,value1,field2,value2,field8,value8,field9,value9,...
    field10,value10,field3,value3,field4,value4,field5,value5,field6,value6,...
    field7,value7,field11,value11,field12,value12,field13,value13);
save(folder, '-struct', filename);