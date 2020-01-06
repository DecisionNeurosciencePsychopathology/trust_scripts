mfx_policy_mult_l = load('/Users/polinavanyukov/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/clinical_mfx/policy/HC/learn_policy_policy_no_multisession_FourTrusteekappaST_vba_mfx_L.mat');
mfx_policy_mixed_l = load('/Users/polinavanyukov/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/clinical_mfx/mixed/HC/learn_mixed_mixed_no_multisession_FourTrusteekappaST_vba_mfx_L.mat');
mfx_rcounter_l=load('/Users/polinavanyukov/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/clinical_mfx/rcounter/HC/learn_rcounter_rcounter_no_multisession_FourTrusteekappaST_vba_mfx_L.mat');
mfx_rcounter2_l=load('/Users/polinavanyukov/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/clinical_mfx/rcounte2/HC/learn_rcounter2_rcounter2_no_multisession_FourTrusteekappaST_vba_mfx_L.mat');
mfx_rcounter3_l=load('/Users/polinavanyukov/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/clinical_mfx/rcounte3/HC/learn_rcounter3_rcounter3_no_multisession_FourTrusteekappaST_vba_mfx_L.mat');
mfx_mixed2_l=load('/Users/polinavanyukov/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/clinical_mfx/mixe2/HC/learn_mixed2_mixed2_no_multisession_FourTrusteekappaST_vba_mfx_L.mat');


ll_hc = vertcat(mfx_policy_mult_l.F,mfx_policy_mixed_l.F);
[posterior, out] =VBA_groupBMC(ll_hc); %policy > mixed

ll_hc = vertcat(mfx_rcounter_l.F,mfx_rcounter2_l.F,mfx_rcounter3_l.F);
[posterior, out] =VBA_groupBMC(ll_hc); %rcounter > rcounter2, rcounter3

ll_hc = vertcat(mfx_rcounter_l.F,mfx_policy_mult_l.F);
[posterior, out] =VBA_groupBMC(ll_hc); %ns, but policy is marginally better

ll_hc = vertcat(mfx_policy_mult_l.F,mfx_mixed2_l.F);
[posterior, out] =VBA_groupBMC(ll_hc); %ns, but policy is somewhat better


mfx_policy_mult_l = load('/Users/polinavanyukov/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/clinical_mfx/policy/all/learn_policy_policy_no_multisession_FourTrusteekappaST_vba_mfx_L.mat');
mfx_policy_mixed_l = load('/Users/polinavanyukov/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/clinical_mfx/mixed/all/learn_mixed_mixed_no_multisession_FourTrusteekappaST_vba_mfx_L.mat');
mfx_rcounter_l=load('/Users/polinavanyukov/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/clinical_mfx/rcounter/all/learn_rcounter_rcounter_no_multisession_FourTrusteekappaST_vba_mfx_L.mat');
mfx_rcounter2_l=load('/Users/polinavanyukov/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/clinical_mfx/rcounte2/all/learn_rcounter2_rcounter2_no_multisession_FourTrusteekappaST_vba_mfx_L.mat');
mfx_rcounter3_l=load('/Users/polinavanyukov/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/clinical_mfx/rcounte3/all/learn_rcounter3_rcounter3_no_multisession_FourTrusteekappaST_vba_mfx_L.mat');
mfx_mixed2_l=load('/Users/polinavanyukov/OneDrive/Documents/Project Trust Game/analyses/manuscript clinical/clinical_mfx/mixe2/all/learn_mixed2_mixed2_no_multisession_FourTrusteekappaST_vba_mfx_L.mat');


ll_all = vertcat(mfx_policy_mult_l.F,mfx_policy_mixed_l.F);
[posterior, out] = VBA_groupBMC(ll_all); %policy > mixed

ll_all = vertcat(mfx_rcounter_l.F,mfx_rcounter2_l.F,mfx_rcounter3_l.F);
[posterior, out] =VBA_groupBMC(ll_all); %rcounter > rcounter2, rcounter3

ll_all = vertcat(mfx_rcounter_l.F,mfx_policy_mult_l.F);
[posterior, out] =VBA_groupBMC(ll_all); %qualtiatively, policy > rcounter, but ns

ll_all = vertcat(mfx_policy_mult_l.F,mfx_mixed2_l.F);
[posterior, out] =VBA_groupBMC(ll_all); %ns, but policy is somewhat better

ll_ffx = ll_all([1 2 3], :); 
VBA_groupBMC(ll_ffx); %3 preferred
VBA_groupBMC(ll_all([1 2],:));
ll_mfx = ll_all([4 5],:);
VBA_groupBMC(ll_mfx); % not conclusive
ll_ffx_v_mfx_st = ll_all([2 5],:);
VBA_groupBMC(ll_ffx_v_mfx_st); %mfx preferred
ll_ffx_v_mfx_s = ll_all([3 4],:);
VBA_groupBMC(ll_ffx_v_mfx_s); %mfx preferred



% now testing kappaST + nomult + different evo fx
null_nomult_kappaST_ffx_l = load('/Users/alisonmarie526/Box/DEPENd/Projects/Social_Trust/output/vba_mfx/matfiles_null_null_no_multisession_kappaST_vba_ffx_L.mat');
null_nomult_kappaST_mfx_l = load('/Users/alisonmarie526/Box/DEPENd/Projects/Social_Trust/output/vba_mfx/matfiles_null_null_no_multisession_kappaST_vba_mfx_L.mat');
regret_nomult_kappaST_ffx_l = load('/Users/alisonmarie526/Box/DEPENd/Projects/Social_Trust/output/vba_mfx/matfiles_regret_regret_no_multisession_kappaST_vba_ffx_L.mat');
regret_nomult_kappaST_mfx_l = load('/Users/alisonmarie526/Box/DEPENd/Projects/Social_Trust/output/vba_mfx/matfiles_regret_regret_no_multisession_kappaST_vba_mfx_L.mat');
countertrustee_nomult_kappaST_ffx_l = load('/Users/alisonmarie526/Box/DEPENd/Projects/Social_Trust/output/vba_mfx/matfiles_countertrustee_countertrustee_no_multisession_kappaST_vba_ffx_L.mat');
countertrustee_nomult_kappaST_mfx_l = load('/Users/alisonmarie526/Box/DEPENd/Projects/Social_Trust/output/vba_mfx/matfiles_countertrustee_countertrustee_no_multisession_kappaST_vba_mfx_L.mat');
policy_nomult_kappaST_ffx_l = load('/Users/alisonmarie526/Box/DEPENd/Projects/Social_Trust/output/vba_mfx/matfiles_policy_policy_no_multisession_kappaST_nomult_vba_ffx_L.mat');
policy_nomult_kappaST_mfx_l = load('/Users/alisonmarie526/Box/DEPENd/Projects/Social_Trust/output/vba_mfx/matfiles_policy_policy_no_multisession_kappaST_nomult_vba_mfx_L.mat');
ptrustee_nomult_kappaST_mfx_l = load('/Users/alisonmarie526/Box/DEPENd/Projects/Social_Trust/output/vba_mfx/matfiles_ptrustee_ptrustee_no_multisession_kappaST_vba_mfx_L.mat');
prs_nomult_kappaST_mfx_l = load('/Users/alisonmarie526/Box/DEPENd/Projects/Social_Trust/output/vba_mfx/matfiles_prs_prs_no_multisession_kappaST_vba_mfx_L.mat');

ll_mfx = vertcat(null_nomult_kappaST_mfx_l.F,regret_nomult_kappaST_mfx_l.F, countertrustee_nomult_kappaST_mfx_l.F,policy_nomult_kappaST_mfx_l.F, ptrustee_nomult_kappaST_mfx_l.F, prs_nomult_kappaST_mfx_l.F);
[a.mfx, r.mfx] = VBA_groupBMC(ll_mfx)
ll_ffx = vertcat(null_nomult_kappaST_ffx_l.negF', regret_nomult_kappaST_ffx_l.negF',countertrustee_nomult_kappaST_ffx_l.negF',  policy_nomult_kappaST_ffx_l.negF');
[a, rr] = VBA_groupBMC(ll_ffx)
