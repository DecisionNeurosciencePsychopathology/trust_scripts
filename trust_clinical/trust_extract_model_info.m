function trust_extract_model_info()
close all;

hallquist = 0;
counter = 2;
mltrun = 1;
kappa = 2;
censor = 0;
files = dir(strcat('*',sprintf('cntr%d_mltrun%d_kappa%d_censor%d', counter, mltrun, kappa, censor),'.mat'));
num_of_subjects = length(files);
%PEs=zeros(num_of_subjects,193); %MODEL PEs
%value=zeros(num_of_subjects,193); %value

%for ct = 1
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
    condition_order = strings(size(ConditionOrder));
    [condition_order{:}] = ConditionOrder{:};
    if blocks == 4
        k_good = condition_order == "good";
        k_bad = condition_order == "bad";
        k_neutral = condition_order == "neutral";
        k_computer = condition_order == "computer";
    else %hallquist dataset
        k_good = condition_order == "good";
        k_bad = condition_order == "bad";
        k_neutral = condition_order == "neutral";
    end
    %saving ID, learning_rate, temperature, subject_specific bias, 4
    %trustee_specific biases, model fits (AIC, BIC, R2)
    model_diag(ct,1:12) = [subject_id, posterior.muTheta, posterior.muPhi(1:2)',posterior.muPhi(find(k_good, 1)), posterior.muPhi(find(k_bad,1)), posterior.muPhi(find(k_neutral, 1)), posterior.muPhi(find(k_computer, 1)), out.fit.AIC, out.fit.BIC, out.fit.R2, out.fit.LL];
    close all;
end    
   
% M_name = sprintf('modelPEs_counter%d_multisession%d_fixed%d_SigmaKappa%d_reputation%d_humanity%d_valence_p%d_valence_n%d_assymetry_choice%d_regret%d',counter, multisession, fixed_params_across_runs, sigma_kappa, reputation_sensitive, humanity, valence_p, valence_n, assymetry_choices, regret);
M_name = ['model_diagnostics'];
filename = char(M_name);
save(filename,'model_diag');

close all;
