function [y, u, noresponse, hasComputer, trustees] = trust_get_data(datafile, vo)

if ~exist(datafile, 'file'), error('Cannot file data file: %s', datafile); end

df = load(datafile);

b=df.b; %output structure from task

ntrials = length(b.TrialNumber);

if ~isfield(df, 'decisions')
    share =~cellfun(@isempty,strfind(b.PartDecides,'share'));
    keep =~cellfun(@isempty,strfind(b.PartDecides,'keep'));
    missed = ~cellfun(@isempty,strfind(b.PartDecides,'noresponse'));
    b.decisions = zeros(ntrials, 1);
    b.decisions(share) = 1;
    b.decisions(keep) = -1;
    b.decisions(missed) = 0;
end

if isfield(df, 'share')
    actions = df.share(1:ntrials)'; %subject's actions
else
    share = zeros(length(b.decisions),1);
    share(b.decisions==1) = 1;
    actions = share(1:ntrials);
end

if ~isfield(df, 'noresponse')
    noresponse = zeros(length(b.decisions),1);
    noresponse(b.decisions==0|b.decisions==-999)=1;
end

actions = double(actions);
u(1,:) = actions;

%rewards = double(strcmp(b.TrusteeDecides(b.Order_RS>-999),'share')); %rewards including counterfactual ones (trustee's actions)
rewards = double(strcmp(b.TrusteeDecides,'share')); %rewards including counterfactual ones (trustee's actions)
rewards(rewards==0) = -1;
%rewards(noresponse'==1) = -999;
u(2,:) = rewards(1:ntrials)';

trustees = unique(b.identity);
if ismember('neutral', trustees)
    reference = 'neutral'; %kappaS becomes intercept in dummy coding scheme
else
    reference = trustees{1}; %first level
end

%move reference to front of vector
trustees = [trustees(strcmpi(reference, trustees)); trustees(~strcmpi(reference, trustees))];
[~, u(9,:)] = ismember(b.identity, trustees); %position in trustees vector of trustee on each trial (for indexing kappa)

% Experimental design
%reputation vector
if strcmpi(vo.dataset,'learn') || (exist('vo.fourTrustee') && vo.fourTrustee)
    %This is lina's learn dataset's coding...
    trustee_ID = zeros(length(b.identity),1);
    trustee_ID(strcmp(b.identity,'computer')) = 1; 
    trustee_ID(strcmp(b.identity,'good')) = 2;  
    trustee_ID(strcmp(b.identity,'bad')) = 3;
    trustee_ID(strcmp(b.identity,'neutral')) = 0; 
    u(3,:) =  trustee_ID(1:ntrials);
else
    trustee_ID = zeros(length(b.identity),1);
    trustee_ID(strcmp(b.identity,'good')) = 1;  % Lina's initial coding = 2
    trustee_ID(strcmp(b.identity,'bad')) = -1;
    trustee_ID(strcmp(b.identity,'neutral')) = 0; % Lina: -1
    u(3,:) =  trustee_ID(1:ntrials);
end

%reward contingency vector
reversal = zeros(length(b.identity),1);
fifty = find(b.Reversal(0+1:length(b.Reversal)),50);
eightyeight = find(b.Reversal(0+1:length(b.Reversal)),88);
twentyfive = find(b.Reversal(0+1:length(b.Reversal)),25);
reversal(fifty) = 0;
reversal(eightyeight)=1;
reversal(twentyfive)=-1;
u(8,:) = reversal(1:ntrials);

% humanity vector
human = ones(length(b.identity),1);
human(strcmp(b.identity,'computer')) = 0;
u(5,:) = human(1:ntrials);

% positive/negative valence vector
valence_pos = zeros(length(b.identity),1);
valence_neg = zeros(length(b.identity),1);
valence_pos(strcmp(b.identity,'good')) = 1;
valence_neg(strcmp(b.identity,'bad')) = 1;
u(6,:) = valence_pos(1:ntrials);
u(7,:) = valence_neg(1:ntrials);

%initial state value only to be sensitive
index_vector = zeros(length(b.identity),1);
blocklength = 48;
index_vector([1,1+blocklength,1+2*blocklength, 1+3*blocklength],1) = 1;
u(4,:) = index_vector(1:ntrials);

%if ntrials > 144 %include practice trials
%    u = u(:,ntrials - 143:end);
%    noresponse = noresponse(ntrials - 143:end);
%    b.identity = b.identity(ntrials - 143:end);
%end

y = u(1,:); %the subject's actions

hasComputer = any(strcmp(b.identity,'computer'));
% Optional censoring a block of trials
if vo.censor == 1
    u = u(:,~strcmp(b.identity,'computer'));
    y = y(:,~strcmp(b.identity,'computer'));
    noresponse = noresponse(~strcmp(b.identity,'computer'));
end

% shifting U
u = [zeros(9,1) u(:,1:end-1)];

end
