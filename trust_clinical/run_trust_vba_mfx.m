function [p_sub,o_sub,p_group,o_group] = run_trust_vba_mfx(vba_df)

%Ensure the vba evo and observation functions are on path
addpath('vba_funcitons\')

%Parse out the vba dataframe
y = vba_df.y;
u = vba_df.u;

models = {'f_trust_Qlearn_policy'};

for model = models
    
    close all; %Get rid of extensive gpu memory figs
    
    %Load in the options per model -- this should be changed to the file
    %name?
    load(sprintf('vba_mfx_input/vba_mfx_input_%s_reordered.mat',model{:}))
    options = vba_df.options;
        
    %Initialize the output matrix
    vba_mfx_df = table();
    
    %Load in the options 
    dim = vba_df.options{1}.dim;
    
    %Designate the f and g function handles
    f_name = vba_df.options{1}.f_fname; %Evolution function
    g_name = vba_df.options{1}.g_fname; %Observation function

    %clear vars for new output
    clearvars p_sub o_sub p_group o_group
    
    %Run the main MFX function
    [p_sub,o_sub,p_group,o_group] = VBA_MFX(y,u,f_name,g_name,dim,options);
    
    %Save uncompressed albeit they will be large
    save(sprintf('E:/data/trust/vba_mfx_input/vba_mfx_output_p_sub_reodered_%s',model{:}),'p_sub', '-v7.3')
    save(sprintf('E:/data/trust/vba_mfx_input/vba_mfx_output_o_sub_reodered_%s',model{:}),'o_sub', '-v7.3')
    save(sprintf('E:/data/trust/vba_mfx_input/vba_mfx_output_p_group_reodered_%s',model{:}),'p_group', '-v7.3')
    save(sprintf('E:/data/trust/vba_mfx_input/vba_mfx_output_o_group_reodered_%s',model{:}),'o_group', '-v7.3')
    
end