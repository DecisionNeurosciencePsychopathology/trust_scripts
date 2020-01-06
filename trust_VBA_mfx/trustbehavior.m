function [b, filepath] = trustbehavior(data_dir,id,force)
% Polina Vanyukov, Jon Wilson & Alex Dombrovski
% 2014-05: revision of the program for the trust game paradigm
% Program that reads in a single data file
% 2019-09-09: Jiazhou modified it to be more universal...
filepath = fullfile(data_dir,id,sprintf('trust%s.mat', id));
if exist(filepath) && ~force
    load(filepath,'b'); 
    b.id=id;
    return
end
file = dir(fullfile(data_dir,id,'trust*_scan*.txt'));
if isempty(file) || length(file)>1
    b=[];
    filepath=[];
    return
end
try
fname = fullfile(file.folder,file.name);
%fname = sprintf('trust_05222015_scan_rs3-%s-1.txt', id5);
%fname = 'trust_05192015_scan_rs4-0064-1.txt';
%Added more fields for timing checking
debug_fields = {'displaychoice.OnsetTime', 'displaychoice.OffsetTime', ...
    'firstFixation.OnsetTime'};

%JW: Might be better to put these in order from txt file
fields = {'ConditionOrder','identity','Running','partnerchoice.RESP',...
     'TrusteeDecides','TrialNumber','exchangeNum','Reversal', ...
     'partnerchoice.OnsetTime','partnerchoice.RTTime','partnerchoice.RT'...
     ,'partnerchoice.Duration','partnerchoice.RESP','PartDecides','displaychoice.Duration',...
     'outcome.OnsetTime','outcome.OffsetTime','Order_RS', ...
     'ITIfixation.OnsetTime','ITIfixation.Duration',...
     'ITIfixation.OffsetTime', 'ITIjitter','partnerchoice.OffsetTime', debug_fields{:}};
startoffset = 13; % trust: from TrialProc occurrence to the beginning of the block containing relevant info
endoffset = 35;   % trust: from TrialProc occurrence to the end of the block containing relevant info
b = read_in_trust(fname,'TrialProc', fields, startoffset, endoffset);

start = 0;
trials = max(b.TrialNumber);
%Decisions share = 1; keep = -1; no response = 0;
share =contains(b.PartDecides(start+1:trials),'share');
keep =contains(b.PartDecides(start+1:trials),'keep');
noresponse = contains(b.PartDecides(start+1:trials),'noresponse');
b.decisions = zeros(trials-start, 1);
b.decisions(share) = 1;
b.decisions(keep) = -1;
b.decisions(noresponse) = 0;

trustee1 = 'none';
trustee2 = 'none';
trustee3 = 'none';
trustee4 = 'none';


%get identity for each phase of the task
for ind=start+1:trials
    if strcmp(trustee1, 'none')|| strcmp(trustee1, b.identity(ind))
        trustee1 = b.identity(ind);
    elseif strcmp(trustee2, 'none')||strcmp(trustee2, b.identity(ind))
        trustee2 = b.identity(ind);
    elseif strcmp(trustee3, 'none')||strcmp(trustee3, b.identity(ind))
        trustee3 = b.identity(ind);
    else
        trustee4 = b.identity(ind);
    end
end

share1_50 = 0;
share2_50 = 0;
share3_50 = 0;
share4_50 = 0;

share1_88 = 0;
share2_88 = 0;
share3_88 = 0;
share4_88 = 0;

share1_25 = 0;
share2_25 = 0;
share3_25 = 0;
share4_25 = 0;

share1Total = 0;
share2Total = 0;
share3Total = 0;
share4Total = 0;

b.trustee1 = {trustee1, share1_50, share1_88, share1_25, share1Total};
b.trustee2 = {trustee2, share2_50, share2_88, share2_25, share2Total};
b.trustee3 = {trustee3, share3_50, share3_88, share3_25, share3Total};
b.trustee4 = {trustee4, share4_50, share4_88, share4_25, share4Total};

% %share responses for different phases of the task
for ind=start+1:trials
    if strcmp(trustee1, b.identity(ind)) && strcmp(b.PartDecides(ind), 'share')
        share1Total = share1Total + 1;
        if b.Reversal(ind) == 50
            share1_50 = share1_50 + 1;
        end
        if b.Reversal(ind) == 25
            share1_25 = share1_25+1;
        end
        if b.Reversal(ind) == 88
            share1_88 = share1_88+1;
        end
        b.trustee1 = {trustee1, share1_50, share1_88, share1_25, share1Total};
    elseif strcmp(trustee2, b.identity(ind)) && strcmp(b.PartDecides(ind), 'share')
        share2Total = share2Total + 1;
        if b.Reversal(ind) == 50
            share2_50 = share2_50 + 1;
        end
        if b.Reversal(ind) == 25
            share2_25 = share2_25+1;
        end
        if b.Reversal(ind) == 88
            share2_88 = share2_88+1;
        end
        b.trustee2 = {trustee2, share2_50, share2_88, share2_25, share2Total};
    elseif strcmp(trustee3, b.identity(ind)) && strcmp(b.PartDecides(ind), 'share')
        share3Total = share3Total + 1;
        if b.Reversal(ind) == 50
            share3_50 = share3_50 + 1;
        end
        if b.Reversal(ind) == 25
            share3_25 = share3_25+1;
        end
        if b.Reversal(ind) == 88
            share3_88 = share3_88+1;
        end
        b.trustee3 = {trustee3, share3_50, share3_88, share3_25, share3Total};
     elseif strcmp(trustee4, b.identity(ind)) && strcmp(b.PartDecides(ind), 'share')
         share4Total = share4Total + 1;
         if b.Reversal(ind) == 50
             share4_50 = share4_50 + 1;
         end
         if b.Reversal(ind) == 25
             share4_25 = share4_25+1;
         end
         if b.Reversal(ind) == 88
             share4_88 = share4_88+1;
         end
         b.trustee4 = {trustee4, share4_50, share4_88, share4_25, share4Total};
    end
end
b.id=id;
save(filepath, 'b');

%%
return
catch
    b=[];
    filepath=[];
end
end


%JW: Functions
%If erroneous subjs creep up
function throw_trust_error(id)
global skip_subjs
    skip_subjs = [skip_subjs  id];
    fprintf(['\n\nTo many files for this subj review for later!',...
        'OR to something is wrong with the trials!\n'])
    save skip_subjs skip_subjs
end


