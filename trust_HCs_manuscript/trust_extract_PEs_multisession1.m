function trust_extract_PEs_multisession1()
close all;

hallquist = 0;
counter = 2;
mltrun = 1;
kappa = 2;
censor = 0;
files = dir(strcat('*',sprintf('cntr%d_mltrun%d_kappa%d_censor%d', counter, mltrun, kappa, censor),'.mat'));
num_of_subjects = length(files);
PEs=zeros(num_of_subjects,193); %MODEL PEs
value=zeros(num_of_subjects,193); %value
for ct = 1:num_of_subjects
    filename=files(ct).name;
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
    for block=1:blocks
        jump = (block-1)*48;
        if block == 1
            PEs(ct,1:49) = [subject_id, out.suffStat.muX(2,1:48)]; %even lines are PEs
            value(ct,1:49) = [subject_id, out.suffStat.muX(1,1:48)];%odd lines are values
        else
            PEs(ct,2+jump:49+jump) = out.suffStat.muX(2+(block-1)*2,(1+jump):(48+jump));
            value(ct,2+jump:49+jump) = out.suffStat.muX(1+(block-1)*2,(1+jump):(48+jump));
        end
    end
    close all;
end
    
    if ~exist('PEsandValues', 'dir')
        mkdir('PEsandValues')
    end
    
% M_name = sprintf('modelPEs_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d_regret%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);
M_name = ['PEsandValues' filesep sprintf('modelPEs_cntr%d_mltrun%d_kappa%d_censor%d',counter, mltrun, kappa, censor)];
filename = char(M_name);
save(filename,'PEs');

%N_name = sprintf('values_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices);
N_name = ['PEsandValues' filesep sprintf('values_cntr%d_mltrun%d%_kappa%d_censor%d',counter, mltrun, kappa, censor)];
save(char(N_name), 'value');
close all;
