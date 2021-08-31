function [] = trust_group()
%% Grouping data from all subjects that have been pre-processed and
% saved by trust behavior as a .mat file

%column 1 = subject;
%column 2 = trial number;
%column 3 = condition;
%column 4 = exchange number;
%column 5 = reward schedule
%column 6 = exchange within a reward schedule within a condition (1-16);
%column 7 = subject decides;
%column 8 = trustee decides;
%column 9 = decision onset time;
%column 10 = decision response time;
%column 11 = feedback onset time;
%column 12 = feedback offset time;

clear all;

%data_location = '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior';
data_location = 'C:\Users\timot\Box\skinner\personal_folders\Tim\trust_bsocial\data\processed';
%data_location = 'C:\Users\timot\Desktop\s1\';
cd(data_location)

files_mat = dir('trust*.mat');
num_of_subjects = length(files_mat);

%output_filename = sprintf(data_dir_str,'/forspss');


c.Names = {'subject','trialnum','trustee','exchange','reward_schedule',...
    'exch2','s_decision','t_decision','decision_Onset','decision_RT',...
    'feedback_Onset','feedback_Offset'};
c.Data = zeros([0,12]);

%% subject loop
for index=1:num_of_subjects
    file_name = files_mat(index).name;
    load(file_name);
    subject_ID = str2double(file_name(isstrprop(file_name,'digit')));
    fprintf('File processing: %s\n', file_name);
    start = 0;
    trials = length(b.TrialNumber);
    conditions= length(unique(b.ConditionOrder));
       
    %% identifying the beginning of nonpractice trials
    if conditions > 3 % for subs with 4 conditions instead of 3
        if length(b.TrialNumber) > 192
        start = length(b.TrialNumber) - 192;
        end
    else % for subs with 3 conditions
        if length(b.TrialNumber) > 144
        start = length(b.TrialNumber) - 144;
        end
    end
      
    c.s(index).trialnum = transpose([1:trials]);
    c.s(index).exchange = transpose(repmat([1:48], 1, conditions));
    
    %% correcting for incomplete experimental sessions
    if conditions > 3 % for subs with 4 conditions instead of 3
        if length(b.Reversal)<192
            c.s(index).reward_schedule(1:length(b.Reversal)) = b.Reversal(start+1:length(b.Reversal));
            c.s(index).reward_schedule(length(b.Reversal)+1:trials)=99;
            c.s(index).reward_schedule = c.s(index).reward_schedule';
        else
            c.s(index).reward_schedule = b.Reversal(start+1:trials);
        end
    else %for subs with 3 conditions
        if length(b.Reversal)<144
            c.s(index).reward_schedule(1:length(b.Reversal)) = b.Reversal(start+1:length(b.Reversal));
            c.s(index).reward_schedule(length(b.Reversal)+1:trials)=99;
            c.s(index).reward_schedule = c.s(index).reward_schedule';
        else
            c.s(index).reward_schedule = b.Reversal(start+1:trials);
        end
    end
    
    %% Condition = good = 1; bad = 2; neutral = 3; computer = 4;
    good =~cellfun(@isempty,strfind(b.identity(start+1:trials),'good'));
    bad =~cellfun(@isempty,strfind(b.identity(start+1:trials),'bad'));
    neutral = ~cellfun(@isempty,strfind(b.identity(start+1:trials),'neutral'));
    comp =~cellfun(@isempty,strfind(b.identity(start+1:trials),'computer'));
    
    c.s(index).condition = zeros(trials-start, 1);
    c.s(index).condition(good) = 1;
    c.s(index).condition(bad) = 2;
    c.s(index).condition(neutral) = 3;
    c.s(index).condition(comp) = 4;
    
    %% Reward schedule = 50% = 0; 25% = -1; 88% = 1;
    fifty = find(b.Reversal(start+1:length(b.Reversal)),50);
    eightyeight = find(b.Reversal(start+1:length(b.Reversal)),88);
    twentyfive = find(b.Reversal(start+1:length(b.Reversal)),25);
    
    c.s(index).reversal=zeros(trials-start,1);
    c.s(index).reversal(fifty) = 0;
    c.s(index).reversal(eightyeight)=1;
    c.s(index).reversal(twentyfive)=-1;
    
    c.s(index).exch2 = transpose(repmat([1:16], 1, conditions*3));
   
    
    %% Decisions share = 1; keep = -1; no reponse = 0;
    share =strcmp(b.PartDecides(start+1:trials),'share');
    keep =strcmp(b.PartDecides(start+1:trials),'keep');
    noresponse = strcmp(b.PartDecides(start+1:trials),'noresponse');;
    
    c.s(index).s_decides = zeros(trials-start, 1);
    c.s(index).s_decides(share) = 1;
    c.s(index).s_decides(keep) = -1;
    c.s(index).s_decides(noresponse) = 0; 
    
    %% Decisions share = 1; keep = -1; no reponse = 0;
    share =~cellfun(@isempty,strfind(b.TrusteeDecides(start+1:trials),'share'));
    keep =~cellfun(@isempty,strfind(b.TrusteeDecides(start+1:trials),'keep'));
    
    c.s(index).t_decides = zeros(trials-start, 1);
    c.s(index).t_decides(share) = 1;
    c.s(index).t_decides(keep) = -1;
    
    %% Timing of events
    c.s(index).decision_Onset = b.partnerchoice_OnsetTime(start+1:trials);
    c.s(index).decision_RT = b.partnerchoice_RT(start+1:trials);
    c.s(index).feedback_Onset = b.outcome_OnsetTime(start+1:trials);
    c.s(index).feedback_Offset = b.outcome_OffsetTime(start+1:trials);
    
    file=file_name(length(file_name)-9:length(file_name)-4);
    if str2num(file)==[]
        file=file_name(length(file_name)-8:length(file_name)-4);
    end
%     c.s(index).subject = repmat(str2num(file),trials-start,1);
      c.s(index).subject = repmat(subject_ID,trials-start,1);
  
    c.Data = vertcat(c.Data,(horzcat(c.s(index).subject,c.s(index).trialnum,...
        c.s(index).condition,c.s(index).exchange, c.s(index).reward_schedule,...
        c.s(index).exch2,c.s(index).s_decides,...
        c.s(index).t_decides, c.s(index).decision_Onset,c.s(index).decision_RT,...
        c.s(index).feedback_Onset, c.s(index).feedback_Offset)));
    
end
%save4spss_mac(c.Names, c.Data, output_filename);
mkdir group_data
cd group_data
%save scan_behavioral_group c; %everything
save beha_behavioral_group c; %everything
x = c.Data;
% data_dir_str = '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior/group_data';
% save(sprintf(strcat(data_dir_str,'/scan_behavior_data.mat')),'x'); %data only

%data_dir_str = '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/beha_behavior/group_data';
data_dir_str = 'C:\Users\timot\Box\skinner\personal_folders\Tim\trust_bsocial\data\processed\group_data';
%data_dir_str = 'C:\Users\timot\Desktop\s1\gd\';
save(sprintf(strcat(data_dir_str,'/beha_behavior_data.mat')),'x'); %data only



end

