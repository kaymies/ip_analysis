clear

%% Subject parameters
datainfo_byTrial = readtable('../Gruben_Data/BDSinfo.xlsx');
% datainfo_Gruben = readtable('./ProcessedData/Santos2016_zIPfits_Gruben.xlsx');

%% Get indices of subjects included in Prof. Gruben's set (assumed to be healthy adults)

iS = 0;
count_subject = 0;
for idatafile = 1:1930
% 
    if datainfo_byTrial.Subject(idatafile) ~= iS 
% 
%         iS = datainfo_byTrial.Subject(idatafile);
%         datainfo_bySubj(iS,:) = datainfo_byTrial(idatafile,:);

        % Check whether subject is included in Prof. Gruben's set
        if ismember(str2double(datainfo_byTrial.Trial{idatafile}(4:end)),...
                GTrial)
            count_subject = count_subject + 1;
            subjects_GrubenSet(count_subject) = iS;
            datainfo_bySubj_GrubenSet(count_subject,:) = ... 
                datainfo_byTrial(idatafile,:); 
        end

    end
end

% All subject data
filename = './ProcessedData/Santos2016_BDSinfo_bySubj.xlsx';
writetable(datainfo_bySubj,filename,'Sheet',1);

% Only subjects that are also included in Prof. Gruben's set
filename = './ProcessedData/Santos2016_BDSinfo_bySubj_GrubenSet.xlsx';
writetable(datainfo_bySubj_GrubenSet,filename,'Sheet',1);