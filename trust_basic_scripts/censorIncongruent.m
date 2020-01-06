function b = censorIncongruent()
%% Grabbing all subjects that have already been processed and turned 
%into .mat file. Then making regressors for each one.


%% Grabbing files
data_dir_str= '/Users/polinavanyukov/Box Sync/Project Trust Game/data/processed/scan_behavior';
data_dump_str = '/Users/polinavanyukov/Box Sync/Project Trust Game/regs/';

if ~exist(data_dump_str,'file')
    mkdir(data_dump_str)
    fprintf('Creating id specific reg folder in: %s\n\n',data_dump_str);
end

cd(data_dir_str)
files = dir('trust*.mat');
num_of_subjects = length(files);

%hard coded times
dur_choice_display = 300;
dur_feedback = 1200; %1100msec in the file?

%scanning parameters
scan_tr = 1.67;
block_length = 155;
block_end = block_length*scan_tr*1000; %in msec
hemoir = spm_hrf(scan_tr, [6,16,1,1,6,0,32]); % better than resampling and smoothing 
frequency_scale_hz = 10;
% this scale is in msec, but it is separated into bins of X
% Hz (defined by 'frequency_scale' above. the resulting
% output will be in the scale of X Hz.
bin_size = 1/frequency_scale_hz*1000; % convert Hz to ms

%% Subject loop
%for index = 8
for index=37
%for index=1:num_of_subjects
    filename = files(index).name;
    fprintf('File processing: %s\n', filename);
    id = filename(isstrprop(filename,'digit'));
    if not(str2double(id)==219956||str2double(id)==220017||str2double(id)==208572)
        load(filename);
        data_dump_str = strcat('/Users/polinavanyukov/Box Sync/Project Trust Game/regs/', num2str(id));
        b.regs = [];    
        
    %Subjects may be using different trial versions of the task, 
    %which only affects how the timings for the scanner should be calculated. 
    %For some participants first fixation onset is the ITI fixation(1). 
    %For others, first fixation onset is labeled as such, 2nd fixation onset 
    %= ITI fixation (1), etc.
        if not(exist('decisions'))
            share =~cellfun(@isempty,strfind(b.PartDecides,'share'));
            keep =~cellfun(@isempty,strfind(b.PartDecides,'keep'));
            noresponse = ~cellfun(@isempty,strfind(b.PartDecides,'noresponse'));
            b.decisions = zeros(192, 1);
            b.decisions(share) = 1;
            b.decisions(keep) = -1;
            b.decisions(noresponse) = 0;
        end
        
        if length(b.ITIfixation_OnsetTime) < 192
            b.ITIfixation_OnsetTime(length(b.ITIfixation_OnsetTime)+1:192)=-999;
            b.ITIfixation_OffsetTime(145:192)=num2cell(-999);
        end
        
        if b.ITIfixation_OnsetTime(192) == -999 || not(iscell(b.firstFixation_OnsetTime))
            firstfix_Onset = b.firstFixation_OnsetTime(1);
            trial2_ITI = 1;
            trial48_ITI = 47;
        else firstfix_Onset = b.ITIfixation_OnsetTime(1);
            trial2_ITI = 2;
            trial48_ITI = 48;
        end

        %for those participants for whom the partnerchoice offset time was not
        %recorded by e-prime
        if iscell(b.partnerchoice_OffsetTime)
            partnerchoice_OffsetTime=b.displaychoice_OnsetTime-1;
        else
            partnerchoice_OffsetTime=b.partnerchoice_OffsetTime;
        end

        fixations = [];
        fixations.event_beg=zeros(48,4);
        fixations.event_end=zeros(48,4); 
        taskness.event_beg =zeros(48,4);
        taskness.event_end=zeros(48,4);
        decision.event_beg=zeros(48,4);
        decision.event_end=zeros(48,4);
        feedback.event_beg=zeros(48,4);
        feedback.event_end=zeros(48,4);
        if iscell(b.partnerchoice_RESP)
            b.missed_trials = b.decisions == 0;
        else       
            b.missed_trials = (b.partnerchoice_RESP==-999);
        end
        b.notmissed_trials = ~(b.missed_trials); 
        
        % trustee decisions
        share =~cellfun(@isempty,strfind(b.TrusteeDecides(1:192),'share'));
        keep =~cellfun(@isempty,strfind(b.TrusteeDecides(1:192),'keep'));
        b.t_shareVSkeep = zeros(192,1);
        b.t_shareVSkeep(share) = 1;
        b.t_shareVSkeep(keep) = -1;
                
        % subject's decisions, share/keep;
        b.shareVSkeep = zeros(size(b.partnerchoice_RESP));
        b.shareVSkeep(b.decisions==1 & b.partnerchoice_RESP ~= -999) = 1;
        b.shareVSkeep(b.decisions==-1 & b.partnerchoice_RESP ~= -999) = -1;
        b.congruent1inNEG1 = b.shareVSkeep.*b.t_shareVSkeep;
        b.congruent1in0 = b.congruent1inNEG1;
        b.congruent1in0(b.congruent1inNEG1==-1)=0;
        b.missNincon = b.missed_trials.*~b.congruent1in0;
        b.notmissNcon = b.notmissed_trials.*b.congruent1in0;
        
%         plot(b.shareVSkeep(1:48));
%         hold on;
%         plot(b.t_shareVSkeep(1:48),'r');
%         figure(2);
%         plot(b.congruent1in0(1:48));

        %plot(b.notmissed_trials);
        trial1_index = 1;
        trial48_index = 48;
        
        for block= 1:4
            %for fixation screens
            if block == 1 && (b.ITIfixation_OnsetTime(192) == -999 || not(isempty([b.firstFixation_OnsetTime])))
                fixations.event_beg(:,block) = [0; b.ITIfixation_OnsetTime(trial2_ITI:trial48_ITI)-firstfix_Onset];
            else
                %fixations.event_beg(:,block) = b.ITIfixation_OnsetTime(trial2_ITI-1:trial48_ITI)-firstfix_Onset+block_end*(block-1);
                fixations.event_beg(:,block) = [155*1670*(block-1); b.ITIfixation_OnsetTime(trial2_ITI:trial48_ITI)-firstfix_Onset+155*1670*(block-1)];
            end        
            %fixations.event_end(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial48_index) - firstfix_Onset+block_end*(block-1);
            fixations.event_end(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial48_index) - firstfix_Onset+155*1670*(block-1);

            %for trial onset to offset; Taskness
            taskness.event_beg(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial48_index)-firstfix_Onset+155*1670*(block-1);
            taskness.event_end(:,block) = b.outcome_OffsetTime(trial1_index:trial48_index)-firstfix_Onset+155*1670*(block-1);

            %for decision onset to response (motor response)
            decision.event_beg(:,block) = b.partnerchoice_OnsetTime(trial1_index:trial48_index)-firstfix_Onset+155*1670*(block-1);
            decision.event_end(:,block) = partnerchoice_OffsetTime(trial1_index:trial48_index)-firstfix_Onset+155*1670*(block-1); 

            %for feedback onset to offset
            feedback.event_beg(:,block) = b.outcome_OnsetTime(trial1_index:trial48_index)-firstfix_Onset+155*1670*(block-1);
            feedback.event_end(:,block) = b.outcome_OffsetTime(trial1_index:trial48_index)-firstfix_Onset+155*1670*(block-1); 
           
            %% New for aligning value: RT aligned to feedback
%             decision.event_beg(:,block) = b.partnerchoice_RTTime(trial1_index:trial48_index)-firstfix_Onset;
%             decision.event_end(:,block) = partnerchoice_OffsetTime(trial1_index:trial48_index)-firstfix_Onset; 

            %epoch window + missed trials + to censor regressor 
            epoch_window = 0:bin_size:taskness.event_end(48, block);
            %% There is a strange glitch, likely owing to the createSimpleRegressor code;
            %% modified by AD to include ITI
            if(sum(b.missed_trials(trial1_index:trial48_index))) > 0
                % write 1s for missed trials to censor
                tmp_reg.(['regressors' num2str(block)]).to_censor = createSimpleRegressor(taskness.event_beg,taskness.event_end, epoch_window, b.missNincon(trial1_index:trial48_index));
            else
                % write a vector of 0s the size of regressors
                tmp_reg.(['regressors' num2str(block)]).to_censor = zeros(size(createSimpleRegressor(taskness.event_beg,taskness.event_end, epoch_window, b.notmissNcon(trial1_index:trial48_index))));
            end
            % flip to_censor
            tmp_reg.(['regressors' num2str(block)]).to_censor = 1-tmp_reg.(['regressors' num2str(block)]).to_censor;

            %% GSR resample
            tmp = gsresample( ...
             [zeros(50,1)' tmp_reg.(['regressors' num2str(block)]).to_censor(1:end-51)], ...
             10,1./scan_tr);
            tmp = ceil(tmp);
            tmp = [tmp ones(1, (block_length-1)-length(tmp))];
            tmp = [tmp zeros(1,155-length(tmp))];
            tmp_reg.(['regressors' num2str(block)]).to_censor = tmp;
            plot(tmp_reg.(['regressors' num2str(block)]).to_censor);
            
            %% Censoring blocks w/ movement
            tmp_reg=censorMovement(id, tmp_reg, block); 
            tmp_reg=censorComputer(id, tmp_reg, block);
            plot(tmp_reg.(['regressors' num2str(block)]).to_censor);
            
            %for next loop iteration, reinitilize variables
            if block < 4
                firstfix_Onset = b.ITIfixation_OnsetTime(trial2_ITI-1+48);
            end
            trial2_ITI=trial2_ITI+48;
            trial48_ITI=trial48_ITI+48;
            trial1_index = trial1_index+48;
            trial48_index = trial48_index+48;       
        end
        
        %plotting to censor regressors
%         plot(tmp_reg.regressors3.to_censor);
%         figure(2);
%         plot(b.notmissed_trials(192-48:192));
        %concatenating
        b.to_censor = [tmp_reg.regressors1.to_censor tmp_reg.regressors2.to_censor tmp_reg.regressors3.to_censor tmp_reg.regressors4.to_censor]; 
 %       plot(b.to_censor);
        b.to_censor = transpose(b.to_censor);
        gdlmwrite(strcat(data_dump_str, 'to_censored_mov_comp_incongr'),[b.to_censor],'\t');


    end
end



return


