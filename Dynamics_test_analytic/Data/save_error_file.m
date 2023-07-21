function [] = save_error_file(input,subject_type,params,filename,error,human_struct)
% folder_name = sprintf('%s_%s',input.Controller.type,input.CoordinateFrame);
folder = fullfile('Data','Error');
addpath(folder);
file = sprintf('%s.mat',filename);
folder = fullfile(folder,file);
%Save data as struct
field1 = 'Input';  value1 = input;
field2 = 'SubjectType'; value2 = subject_type;
field3 = 'NumSubjects'; value3 = human_struct.NumSubjects;
field4 = 'Parameters'; value4 = params;
field5 = 'Error'; value5 = error;
error_data = struct(field1,value1,field2,value2,field3,value3, ...
    field4,value4,...
    field5,value5);
save(folder, '-struct', 'error_data');
end

