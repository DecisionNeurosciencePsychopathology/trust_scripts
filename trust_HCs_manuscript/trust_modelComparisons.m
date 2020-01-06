L_p0n1 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n1.mat');
L_p1n0 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p1_valence_n0.mat');
L_p1n1 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p1_valence_n1.mat');

L_m0r0h0p0n0a0 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0.mat');
L_m0r0h0p0n1a0 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n1_assymetry_choice0.mat');
L_m0r0h0p1n0a0 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p1_valence_n0_assymetry_choice0.mat');
L_m0r1h0p0n0a0 = load('L_multisession0_fixed1_SigmaKappa1_reputation1_humanity0_valence_p0_valence_n0_assymetry_choice0.mat');
L_m0r0h1p0n0a0 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity1_valence_p0_valence_n0_assymetry_choice0.mat');
L_m0r0h0p1n1a0 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p1_valence_n1_assymetry_choice0.mat');
L_m0r0h0p0n0a1 = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice1.mat');
L_k0m0r0h0p0n0a1 = load('L_multisession0_fixed1_SigmaKappa0_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice1.mat');



Lgood = [L_m0r0h0p0n0a0.L;L_m0r1h0p0n0a0.L;L_m0r0h1p0n0a0.L;L_m0r0h0p0n1a0.L;L_m0r0h0p1n0a0.L;L_m0r0h0p1n1a0.L; L_m0r0h0p0n0a1.L; L_k0m0r0h0p0n0a1.L];


%%
clear all; close all;
new = load('Ushifted2_L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_beta0');
Lnew = new.L;
old = load('L_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_beta0');
Lold = old.L;

[posterior,out] = VBA_groupBMC([Lnew; Lold]);

clear all; close all;
old = load('L_counter0_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_beta0');
Lold = old.L;

new = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_beta0');
Lnew = new.L;
Lgood = [Lnew; Lold];
Lgood = Lgood(:,[1:10,12:28,30:47,49]);

[posterior,out] = VBA_groupBMC([Lnew; Lold]);

L_c1m0r0h0p0n0a0 = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_beta0');
L_c1m1r0h0p0n0a0 = load('L_counter1_multisession1_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_beta0');
L_c1m0r1h0p0n0a0 = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation1_humanity0_valence_p0_valence_n0_assymetry_choice0_beta0');
L_c1m0r0h1p0n0a0 = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation0_humanity1_valence_p0_valence_n0_assymetry_choice0_beta0');
L_c1m0r0h0p1n0a0 = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p1_valence_n0_assymetry_choice0_beta0');
L_c1m0r0h0p0n1a0 = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n1_assymetry_choice0_beta0');
L_c1m0r0h0p0n0a1 = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice1_beta0');

Lgood = [L_c1m0r0h0p0n0a0.L;L_c1m0r1h0p0n0a0.L;L_c1m0r0h1p0n0a0.L;L_c1m0r0h0p1n0a0.L;L_c1m0r0h0p0n1a0.L;L_c1m0r0h0p0n0a1.L];

[posterior,out] = VBA_groupBMC(Lgood);

Lgood = [L_c1m0r0h0p0n0a0.L;L_c1m0r1h0p0n0a0.L;L_c1m0r0h1p0n0a0.L];

[posterior,out] = VBA_groupBMC(Lgood);

Lgood = [L_c1m0r1h0p0n0a0.L;L_c1m0r0h1p0n0a0.L];

[posterior,out] = VBA_groupBMC(Lgood([1,2],:));

Lgood = [L_c1m0r1h0p0n0a0.L;L_c1m0r0h1p0n0a0.L;L_c1m0r0h0p1n0a0.L;L_c1m0r0h0p0n1a0.L;L_c1m0r0h0p0n0a1.L];

L_c1m0r0h0p0n0a0u1 = load('L_counter1_multisession0_fixed1_SigmaKappa1_reputation0_humanity0_valence_p0_valence_n0_assymetry_choice0_regret1');
Lgood = [L_c1m0r0h0p0n0a0.L;L_c1m0r0h0p0n0a0u1.L];
[posterior,out] = VBA_groupBMC(Lgood(:,[1:10,12:28,30:47,49])); %L_c1m0r0h0p0n0a0.L > L_c1m0r0h0p0n0a0u1.L

%01/30/2017: comparing mixed counterfactuals+null w/rest
%c0t0r0 = load('L_cntr0_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0'); %null
c0t0r0_pmv = load('L_cntr0_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0_pmv'); %null corrected
c1t0r0 = load('L_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0'); %subject-counterfactual or policy model
%c0t0r1 = load('L_cntr0_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg1');%old regret
newnull_original_choice_rule = load('newnull_original_choice_rule');
c0t0r1_pmv = load('L_cntr0_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg1_pmv'); %regret
c1t1r0 = load('L_cntr1_mltrun0_fixed1_kappa1_rep1_hum0_val_p0_val_n0_as_choices0_reg0'); %trustee-counterfactual
c2t0r0 = load('L_cntr2_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0'); %mixed null+subject-counterfactual
SVM = load('L_SV');

Lgood = [c0t0r0.L;c1t0r0.L;c2t0r0.L]; 
[posterior,out] = VBA_groupBMC(Lgood);
Lgood = [c0t0r0.L;c2t0r0.L]; 
Lgood = [c0t0r0.L;c1t0r0.L;c1t1r0.L;c0t0r1.L]; 
Lgood = [c0t0r0.L;c1t0r0.L;c0t0r1_pmv.L;c1t1r0.L]; 
Lgood = [c2t0r0.L;c0t0r1_pmv.L;c1t1r0.L]; 
Lgood = [c0t0r0_pmv.L;c0t0r1_pmv.L;c1t1r0.L;c1t0r0.L]; %manuscript main comparison
Lgood = [c0t0r0_pmv.L;c0t0r1_pmv.L;c1t1r0.L;c1t0r0.L; SVM.L]; %main comparison


good_index_SVM = [1:6 8:14]; %behavioral
good_index_SVM = [1:16 18:22]; %scan
Lgood = [SVM.L(good_index_SVM); c0t0r0_pmv.L(good_index_SVM);c0t0r1_pmv.L(good_index_SVM);c1t1r0.L(good_index_SVM);c1t0r0.L(good_index_SVM)]; % main comparison w/ SVM
Lgood = [SVM.L(good_index_SVM); c0t0r0_pmv.L(good_index_SVM)]; % svm VS null comparison w/ SVM
Lgood = [SVM.L(good_index_SVM); c1t0r0.L(good_index_SVM)]; % svm VS policy comparison w/ SVM

%08/07/2017: comparing policy model w/ multisession kappa parameter
m0c1t0r0 = load('L_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat'); %policy
m1c1t0r0k4 = load('L_cntr2_mltrun1_fixed1_kappa4_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat'); %multisession kappa parameter

Lgood=[m0c1t0r0.L; m1c1t0r0k4.L];
[posterior,out] = VBA_groupBMC(Lgood);
%08/07/2017: comparing policy model w/ multisession kappa parameter +
%subject-specific kappa parameter
m0c1t0r0 = load('L_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat'); %policy
m1c1t0r0k4 = load('L_cntr2_mltrun1_fixed1_kappa2_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat'); %multisession kappa parameter

Lgood=[m0c1t0r0.L; m1c1t0r0k4.L];
[posterior,out] = VBA_groupBMC(Lgood);

kappa_varying = load('L_cntr2_mltrun1_fixed1_kappa1_free.mat');
kappa_fixed = load('L_cntr2_mltrun1_fixed1_kappa1_multisession_fixed_param.mat');
Lgood=[kappa_fixed.L; kappa_varying.L];

theta_free = load('L_cntr2_mltrun1_fixed1_kappa1_theta_free');

two_kappas_free = load('L_cntr2_mltrun1_fixed1_kappa2_x0_fixed');

two_kappas_free_SigmaKappa1 = load('L_cntr2_mltrun1_fixed1_kappa2_censor0');

%% a more systematic comparison of models, where sigma of kappa has been set to either 1/3 or 1
policy_multi0_kappa_sigmathird = load('L_cntr1_mltrun0_fixed1_kappa1_rep0_hum0_val_p0_val_n0_as_choices0_reg0.mat');
policy_multi0_kappa_sigma1 = load('L_cntr2_mltrun0_fixed1_kappa1_sigma1_censor0');
policy_multi1_kappaS_sigma1 = load('L_cntr2_mltrun1_fixed1_kappa1_free1_x0_fixed_censor0');
policy_multi1_kappaS_kappaT_sigma1 = load('L_cntr2_mltrun1_fixed1_kappa2_free1_x0_fixed_censor0');

Lgood = [policy_multi0_kappa_sigmathird.L,policy_multi0_kappa_sigma1.L];

policy_multi1_kappaS_sigma_third;
policy_multi1_kappaS_kappaT_sigma_third;

%% comparison of policy w/ old choice rule (subject-level kappa only) w/ new choice rule (equivocal results)
old_policy = load('L_cntr2_mltrun0_kappa1_censor0.mat');
actual = load('L_cntr0_mltrun1_kappa2_censor0');
regret = load('L_cntr1_mltrun1_kappa2_censor0');
policy = load('L_cntr2_mltrun1_kappa2_censor0');
trustee = load('L_cntr3_mltrun1_kappa2_censor0');
Lgood = [old_policy.L; actual.L; regret.L; policy.L; trustee.L]; %beha: old_policy > rest, ns; scan: old_policy > rest, ns; hall: old_policy > rest (s)


%% (not censored) a comparison of actual, regret, policy, trustee-counterfactual, all w/ new choice rule (k_s + k_t)

SV = load('L_SVM');
actual = load('L_cntr0_mltrun1_kappa2_censor0');
regret = load('L_cntr1_mltrun1_kappa2_censor0');
policy = load('L_cntr2_mltrun1_kappa2_censor0');
trustee = load('L_cntr3_mltrun1_kappa2_censor0');
regret2 = load('L_cntr4_mltrun1_kappa2_censor0'); %scan: p < .001; beha: p = .888; hallquist: p = .062 (but regret2 won)
Lgood = [actual.L; regret.L; policy.L; trustee.L];
Lgood = [actual.L; regret.L; policy.L; trustee.L; regret2.L];
Lgood = [actual.L; regret.L; policy.L; trustee.L; SV.L];
[posterior,out] = VBA_groupBMC(Lgood);

save(char('vba_4mc_outstats'),'out');

%% comparisons of alternative multisession models

ss = load('L_cntr2_mltrun0_fixed0_kappa1_censor0'); %single session
X0_v_KS_f_KT_0 = load('L_cntr2_mltrun1_fixed1_kappa1_censor0');
X0_f_KS_v_KT_0 = load('L_cntr2_mltrun1_fixed2_kappa1_censor0');
X0_f_KS_f_KT_f = load('L_cntr2_mltrun1_fixed3_kappa2_censor0');
Lgood = [ss.L;X0_v_KS_f_KT_0.L; X0_f_KS_v_KT_0.L; X0_f_KS_f_KT_f.L];
[posterior,out] = VBA_groupBMC(Lgood); %beha: ss > rest (BOR = 0.858, ep = .9998); scan: ss > rest (BOR = 0.982; ep = 0.792); hallquist: x varying > rest (BOR = .435; ep = 1)

Lgood = [X0_v_KS_f_KT_0.L; X0_f_KS_v_KT_0.L; X0_f_KS_f_KT_f.L];
[posterior,out] = VBA_groupBMC(Lgood); %beha: x0 varying > rest (BOR = .729; ep = .979); scan: x0 varying > rest (BOR = 0.916; ep = 0.928); hallquist:x varying > rest (BOR = .008; ep = 1)

Lgood = [X0_f_KS_v_KT_0.L; X0_f_KS_f_KT_f.L];
[posterior,out] = VBA_groupBMC(Lgood); %beha: ks varying < ks + kt (BOR = .524; ep = 1.0); scan: ks varying < ks + kt (BOR = .640; ep = 1.0); hallquist: ks varying < ks + kt (BOR = .718; ep = 1)

%% comparison of all three datasets pooled (computer block censored)

scan_ks1kt0 = load('L_cntr2_mltrun1_fixed2_kappa1_censor1');
scan_ks1kt1 = load('L_cntr2_mltrun1_fixed3_kappa2_censor1');

beha_ks1kt0 = load('L_cntr2_mltrun1_fixed2_kappa1_censor1');
beha_ks1kt1 = load('L_cntr2_mltrun1_fixed3_kappa2_censor1');

hallquist_ks1kt0 = load('L_cntr2_mltrun1_fixed2_kappa1_censor0');
hallquist_ks1kt1 = load('L_cntr2_mltrun1_fixed3_kappa2_censor0');

ks1kt0 = [scan_ks1kt0.L beha_ks1kt0.L hallquist_ks1kt0.L];
ks1kt1 = [scan_ks1kt1.L beha_ks1kt1.L hallquist_ks1kt1.L];
Lgood = [ks1kt0; ks1kt1];
[posterior,out] = VBA_groupBMC(Lgood);

%% comparison of policy w/ one kappa and free x0 and policy w/ 2 kappas (also w/ free learning rate)
scan_ks1kt0x1 = load('L_cntr2_mltrun1_fixed1_kappa1_censor0');
scan_ks1kt1x0 = load('L_cntr2_mltrun1_kappa2_censor0');
scan_ks1kt0x0l1 = load('L_cntr2_mltrun1_fixed4_kappa1_censor0');
scan_ks1kt0x0b1 = load('L_cntr2_mltrun1_fixed5_kappa1_censor0');

beha_ks1kt0x1 = load('L_cntr2_mltrun1_fixed1_kappa1_censor0');
beha_ks1kt1x0 = load('L_cntr2_mltrun1_kappa2_censor0');
beha_ks1kt0x0l1 = load('L_cntr2_mltrun1_fixed4_kappa1_censor0');
beha_ks1kt0x0b1 = load('L_cntr2_mltrun1_fixed5_kappa1_censor0');

hallquist_ks1kt0x1 = load('L_cntr2_mltrun1_fixed1_kappa1_censor0');
hallquist_ks1kt1x0 = load('L_cntr2_mltrun1_kappa2_censor0');
hallquist_ks1kt0x0l1 = load('L_cntr2_mltrun1_fixed4_kappa1_censor0');
hallquist_ks1kt0x0b1 = load('L_cntr2_mltrun1_fixed5_kappa1_censor0');

ks1kt0x1 = [scan_ks1kt0x1.L beha_ks1kt0x1.L hallquist_ks1kt0x1.L];
ks1kt1x0 = [scan_ks1kt1x0.L beha_ks1kt1x0.L hallquist_ks1kt1x0.L];
ks1kt0L1 = [scan_ks1kt0x0l1.L beha_ks1kt0x0l1.L hallquist_ks1kt0x0l1.L];
ks1kt1B1 = [scan_ks1kt0x0b1.L beha_ks1kt0x0b1.L hallquist_ks1kt0x0b1.L];

Lgood = [ks1kt0x0; ks1kt1x1];
Lgood = [ks1kt0x1; ks1kt1x0; ks1kt0L1; ks1kt1B1];
[posterior,out] = VBA_groupBMC(Lgood);

% comparison of new policy w/ SVM
scan_policy = load('L_cntr2_mltrun1_kappa2_censor0');
scan_SVM = load('L_cntr6_mltrun0_fixed1_kappa0_censor0');
beha_policy = load('L_cntr2_mltrun1_kappa2_censor0');
beha_SVM = load('L_cntr6_mltrun0_fixed1_kappa0_censor0');
hallquist_policy = load('L_cntr2_mltrun1_kappa2_censor0');
hallquist_SVM = load('L_cntr6_mltrun0_fixed1_kappa0_censor0');
index = [1:6, 8:15]; %215129, i = 7 is missing from the SV modeling results because does not have ratings.
beha_policy.L = beha_policy.L(index); 

Lgood = [scan_policy.L; scan_SVM.L];
Lgood = [beha_policy.L; beha_SVM.L];
Lgood = [hallquist_policy.L; hallquist_SVM.L];
Lgood = [scan_policy.L beha_policy.L hallquist_policy.L; scan_SVM.L beha_SVM.L hallquist_SVM.L];
[posterior,out] = VBA_groupBMC(Lgood);

% comparison of clinical manuscript dataset

actual = load('L_cntr0_mltrun1_fixed3_kappa2_censor0');
regret = load('L_cntr1_mltrun1_fixed3_kappa2_censor0');
policy = load('L_cntr2_mltrun1_fixed3_kappa2_censor0');
trustee =load('L_cntr3_mltrun1_fixed3_kappa2_censor0');

Lgood = [actual.L; regret.L; policy.L; trustee.L];
[posterior,out] = VBA_groupBMC(Lgood); %bor < .001 eps[0, 0, 1, 0]
