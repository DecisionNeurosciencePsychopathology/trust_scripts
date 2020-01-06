# trust_HCs_manuscript
VBA-style Q-learning and other models of modified Trust Game.
A procedure to run programs in this collection can go as follows. (1) run trust_fit_masterlist_vba. (2) compare models with the trust_modelComparison script.

1. trust_fit_masterlist_vba.m: gets data from all subjects in masterlist (.mat), runs a Qlearning (trust_Qlearning_ushifted_v2) algorithm w/ preset parameters, saves F(it)-values in table.

2a. trust_Qlearning_ushifted_v2.m: VBA fitting of Qlearning to a single subject trust data. Note that the u vector (which includes the subject's behavior +outcomes+experimental design) is shifted by one trial ahead, because enivronmental feedback on trial t-1 is what determines the value of the hidden states on trial t. Calls to f() and g() (obervation/evolution functions, respectively). Saves the figure and the results of the VBA inversion routine. Also is used to run with Social Value Model aka Fareri et al., 2015.

3. Evolution functions:
  3a. f_trust_Qlearn_policy.m: subject-counterfactual *counterfactuals are only calculated when participant's action = KEEP*.
  3b. f_trust_Qlearn_null_pmv.m: simpler evolution function where reinforcement for keep/share is considered separately (actual rewards).
  3c. f_trust_Qlearn_counter_trustee.m: trustee-counterfactual = actual rewards - would-be outcome from trustee's alternative action. 
  3d. f_trust_Qlearn_counter_hybrid_regret.m: regret = actual rewards - maximum outcome (1.5).
  3e. f_trust_SVM1: Social Value Model aka Fareri et al., 2015.

4. Observation functions:
  4a. g_trust_softmax_ED1: observation function that includes the kappa parameter (bias for keep/share actions) *for actual rewards*.
  4b. g_trust_softmax_ED2: observation function that includes the kappa parameter (bias for keep/share actions) *for all other evolution fns*.
  4b. g_trust_softmax_1: simpler observation function.
  4c. g_trust_SVM1.m: Social Value Model aka Fareri et al., 2015 *also softmax*.

5. trust_modelComparisons: helper code to compares F(its) between models. Uses saved tables from trust_fit_group_vba.m.
