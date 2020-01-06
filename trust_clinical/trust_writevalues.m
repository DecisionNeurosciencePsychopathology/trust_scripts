clear variables;
close all;

%needs glob.m (in Project Trust Game/scripts/temporal_instrumental_agent folder

%Quick username check, and path setting, this may have to change depending
%on the machine you are currently working on!
os = computer;
if strcmp(os(1:end-2),'PCWIN')
    values_location = glob('?');
else
    [~, me] = system('whoami');
    me = strtrim(me);
    if strcmp(me,'polinavanyukov')==1
        models_location = glob('/Users/polinavanyukov/Scripts/Trust/trust_rl_VBA/clinical_scan_models/05-Apr-2018f_trust_Qlearn_policy');
        pst_regs_location= glob('/Users/polinavanyukov/Regs/clinical_basic_regs/');
        write_location = strcat('/Users/polinavanyukov/Regs/clinical_policy_values/',date);

    else
        pe_location = glob('?');
        values_location = glob('?');
        pst_regs_location= glob('?');
        write_location = glob('?');
    end
end

cd(models_location{1})
files = dir(strcat('*.mat'));
num_of_subjects = length(files);

%files = dir(strcat('*',sprintf('cntr%d_mltrun%d_kappa%d_censor%d', counter, mltrun, kappa, censor),'.mat'));
if not(exist(write_location, 'dir'))
    mkdir(write_location);
end

hallquist = 0;
counter = 2;
mltrun = 1;
kappa = 2;
censor = 0;


%for index = 9
for index = 1:num_of_subjects    
    filename=files(index).name;
    fprintf('File processing: %s\n', filename);
    %subject_id = filename(isstrprop(filename,'digit'));
    load(filename);
    if ischar(id)
        subject_id = str2double(id);
    else
        subject_id = id;
    end
    blocks = 4;
    if id == 219471 %has only 2 blocks
        blocks = 2;
    elseif hallquist == 1
            blocks = 3;
    end
    if not(exist('b.decisions'))
        share =~cellfun(@isempty,strfind(b.PartDecides,'share'));
        keep =~cellfun(@isempty,strfind(b.PartDecides,'keep'));
        noresponse = ~cellfun(@isempty,strfind(b.PartDecides,'noresponse'));
        b.decisions = zeros(192, 1);
        b.decisions(share) = 1;
        b.decisions(keep) = -1;
        b.decisions(noresponse) = 0;
    end

    share =~cellfun(@isempty,strfind(b.TrusteeDecides,'share'));
    keep =~cellfun(@isempty,strfind(b.TrusteeDecides,'keep'));

    trustee.decisions = zeros(192, 1);
    trustee.decisions(share) = 1;
    trustee.decisions(keep) = -1;

    kept = b.decisions==-1;
    for block =1:blocks
        block_beg = 1+48*(block-1);
        block_end = 48*(block);
        values(block_beg:block_end) = posterior.muX(block*2-1,block_beg:block_end);      
    end
        
         figure(2); clf;
         subplot(2,1,1);
%         hold on;
         plot(values(1:192));
%         subplot(4,1,2);
%         plot(posterior.muX(1,(1:48)));
%         subplot(4,1,3);
%         plot(b.decisions(1:48));
         subplot(2,1,2);
         plot(trustee.decisions(1:192),'red');
        
        values(kept) = values(kept)*-1; %flipping values     
        values = values(2:blocks*48);
        %values = circpshift(values, total_trials-1);
        %values = abs(values);
        values(192) = 0;
        stdev = std(values(1:end-1));
        values(1:end-1) = values(1:end-1) - mean(values(1:end-1)); %z-scoring
        values(1:end-1) = values(1:end-1)./stdev;
        
        %cd(pst_regs_location{1});
        reg_DT = load(strcat(pst_regs_location{:}, num2str(id),'decision_Times.dat'));
        reg_FT = load(strcat(pst_regs_location{:}, num2str(id),'feedback_Times.dat'));

        values_DT = [reg_DT(:,1:2), values'];
        values_FT = [reg_FT(:,1:2), values'];
        dlmwrite([write_location '/trust' num2str(id) 'values_policy_flipped_DT.dat'],values_DT,'delimiter','\t','precision','%.6f');
        dlmwrite([write_location '/trust' num2str(id) 'values_policy_flipped_FT.dat'],values_FT,'delimiter','\t','precision','%.6f');
        close all; 
        
end
