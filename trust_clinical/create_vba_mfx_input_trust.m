function vba_df=create_vba_mfx_input_trust()
%Create the input for VBA_MFX analysis on the trust game clinical data

%END HLEP

models = {'f_trust_Qlearn_policy_clinical'};

for model = models
    
    %Initialize the tbale
    vba_df = table();
    
    %Use glob to pull all the vba files
        data_path = 'E:/Box Sync/skinner/projects_analyses/Project Trust/data/model-derived/clinical_scan_models/';
        vba_files = glob([data_path model{:} '/*.mat']);

    
    %File loop
    for vba_file = vba_files'
        load(vba_file{:}, 'out') %Load in the file should contain b,out,posterior
        close all; %Remove figures
        id = regexp(vba_file{:}, '[0-9]{4,6}', 'match');
        id = str2double(id{:});
        
        %Reorder the multisession y and u inputs for MFX
        %Note I am going by out.design and assuming that was the order the
        %subject recieved the contingency!
        out=remap_trust_inputs(out);
        
        %Initialize temporay dataframe
        tmp_table = table();
        
        %Grab id
        tmp_table.ID = id;
        
        %Grab y & u
        tmp_table.y = {out.y};
        tmp_table.u = {out.u};
        
        %Grab the options used or perhaps create a sub function to create this
        tmp_table.options = {out.options};
        
        vba_df = [vba_df; tmp_table];
    end
    
    %Save the data
    if ~exist([pwd filesep 'vba_mfx_input/'], 'dir')
        mkdir([pwd filesep 'vba_mfx_input/'])
    end
    
    save([pwd filesep 'vba_mfx_input/' sprintf('vba_mfx_input_%s',model{:})], 'vba_df');
end



function out = remap_trust_inputs(in)
%Reorganize the y and u to have the same trustee order for every
%subject the order will be alphabetical b, c, g, n.
subj_order = in.design;
alphabet_order = sort(in.design);
out=in;
num_trials = 48;
seq_idx = 1:48:192;

%Do we have any bad subjects?
if length(in.y)~=192
    error('Figure out what to do in these situations')
else
    %Remap y and u 
    out.y =nan(size(in.y));
    out.u =nan(size(in.u));
    for i = 1:length(alphabet_order)
       %Find where the contingecy happened
       idx=find(ismember(subj_order,alphabet_order{i}));
       start_idx = num_trials*(idx-1)+1;
       end_idx = num_trials*idx;
       out.y(:,seq_idx(i):i*num_trials)=in.y(:,start_idx:end_idx);
       out.u(:,seq_idx(i):i*num_trials)=in.u(:,start_idx:end_idx);
    end
end
       
        
       
       




        


