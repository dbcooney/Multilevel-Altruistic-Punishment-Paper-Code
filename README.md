# Multilevel-Altruistic-Punishment-Paper-Code

This repository accompanies the preprint "Exploring the Evolution of Altruistic Punishment with a PDE Model of Cultural Multilevel Selection", by Daniel B. Cooney. It will include all of the scripts used to generate the figures in the paper.

The repository is organized into three folders: Scripts, Figures, and Simulation Outputs. All of the scripts can be run using Python 3.10.

For reference, below is a list of figures and the scripts that were used to generate each figure.

- Figure 2.1: Run within_trimorphic_altruistic.py
- Figure 4.1: Run group_local_density.py
- Figure 4.2: Run Gdiff_threshold.py
- Figure 4.3: Run Gdiff_threshold.py
- Figure 4.4: Run DPvsDCedgepayoff.py
- Figure 5.1: Run altruistic_fvpairwisegroup.py with "switch_prob = Fermi"
- Figure 5.2: Run altruistic_steady_fv_pairwise.py with "switch_prob = Fermi"
- Figure 5.3: Run loop_altruistic_fvpairwise.py with "switch_prob = Fermi" to generate data, then run q_model_plots_nonlinear.py to produce the figure
- Figure 5.4: Run loop_altruistic_fvpairwise.py with "switch_prob = Fermi" to generate data, then run k_model_plots_nonlinear.py to produce the figure
- Figure 5.5: Run loop_altruistic_fvpairwise.py with "switch_prob = local" to generate data, then run q_model_local_nonlinear_plots.py to produce the figure
- Figure 5.6: Run loop_altruistic_fvpairwise.py with "switch_prob = Tullock" to generate data, then run q_model_plots_Tullock.py to produce the figure
- Figure 6.1: Run trimorphic_altruistic.py with "quantity = trajectory"
- Figure 6.2: Run loop_lambda_trimorphic_altruistic.py with "group_rate_type = fraction cooperating" to generate data, then run cooperation_steady_trimorphic_plots.py to produce the figure
- Figure 6.3: Run trimorphic_altruistic.py with "quantity = steady"
- Figure 6.4: Run loop_lambda_trimorphic_altruistic.py with "group_rate_type = average payoff" to generate data, then run average_payoff_steady_state_trimorphic.py to produce the figure
