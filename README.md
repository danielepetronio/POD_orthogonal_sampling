# POD_orthogonal_sampling
The code 'predictive_analysis.m' performs the POD-orthogonal sampling, starting from a corner of the test matrix (default).
It loads 'all_data.mat' and it uses the functions: 'POD.m' 'predictor_matrix.m' 'figure_cpt.m'.
It is possible to select GP or Lasso as regression methods.
The code creates a .jpg of the cpt distribution modelled with a reduced number of tests (Fig. 12 of the paper) and a .gif of the test matrix progressive filling.
