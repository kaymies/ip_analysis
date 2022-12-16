%% Saving Human Data as Struct
setpath
filename = 'young_dom';
folder = fullfile('Data','Human');
file = sprintf('%s.mat',filename);
folder = fullfile(folder,file);

field1 = 'SubjectType';  value1 = 'Old';
field2 = 'NumSubjects';  value2 = [];
field3 = 'Frequency';  value3 = Frequency;
field4 = 'IPDataAverage';  value4 = mean';
field5 = 'Plane'; value5 = 'Sagittal';
field6 = 'WhichLeg'; value6 = 'Dominant';
field7 = 'StandardDeviation'; value7 = sd';
old_dom = struct(field1,value1,field2,value2,field3,value3,field4,value4,...
    field5,value5,field6,value6,field7,value7);
save(folder, '-struct', filename);