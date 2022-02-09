####################
# B4
run_number_counts(flist=np.logspace(-1.1, -0.4, 3),
  detections_file='B4/B4_SMGs_run3.0.3arcsec.flux.fovscale1.5.txt',
  effective_area_file='B4/B4_effarea_1.5FWHM.txt',
  band='B4',
  default_simulation='./simulations/default.json',
  simulation_folder='simulations/B4',
  flux_mode='gaussian',
  n_points=3,
  savefile='B4/B4_number_counts_fovscale1.5.txt')

run_number_counts(flist=np.logspace(-1.1, -0.4, 3),
  detections_file='B4/B4_SMGs_run3.0.3arcsec.flux.fovscale1.5.txt',
  effective_area_file='B4/B4_effarea_1.5FWHM.txt',
  band='B4',
  default_simulation='./simulations/default.json',
  simulation_folder='simulations/B4',
  flux_mode='gaussian',
  perturbation=0.1,
  n_points=3,
  savefile='B4/B4_number_counts_fovscale1.5_withperturbation.txt')

run_number_counts(flist=np.logspace(-1.1, -0.4, 3),
  detections_file='B4/B4_SMGs_run3.0.3arcsec.flux.fovscale1.5.txt',
  effective_area_file='B4/B4_effarea_1.5FWHM.txt',
  band='B4',
  default_simulation='./simulations/default.json',
  simulation_folder='simulations/B4',
  flux_mode='gaussian',
  n_points=3,
  mode='differential',
  savefile='B4/B4_number_counts_fovscale1.5_differential.txt')

run_number_counts(flist=np.logspace(-1.1, -0.4, 3),
  detections_file='B4/B4_SMGs_run3.0.3arcsec.flux.fovscale1.5.txt',
  effective_area_file='B4/B4_effarea_1.5FWHM.txt',
  band='B4',
  default_simulation='./simulations/default.json',
  simulation_folder='simulations/B4',
  flux_mode='gaussian',
  perturbation=0.1,
  n_points=3,
  mode='differential',
  savefile='B4/B4_number_counts_fovscale1.5_differential_withperturbation.txt')



######################
# B5
run_number_counts(flist=[0.3,0.6],\
  detections_file='B5/B5_SMGs_run3.0.3arcsec.flux.fovscale1.5.txt',
  effective_area_file='B5/B5_effarea_1.5FWHM.txt',
  band='B5',
  default_simulation='./simulations/default.json',
  simulation_folder='simulations/B5',
  flux_mode='gaussian',
  mode='cumulative',
  savefile='B5/B5_number_counts_fovscale1.5.txt')

run_number_counts(flist=[0.3,0.6],\
  detections_file='B5/B5_SMGs_run3.0.3arcsec.flux.fovscale1.5.txt',
  effective_area_file='B5/B5_effarea_1.5FWHM.txt',
  band='B5',
  default_simulation='./simulations/default.json',
  simulation_folder='simulations/B5',
  flux_mode='gaussian',
  mode='cumulative',
  perturbation=0.1,
  savefile='B5/B5_number_counts_fovscale1.5_withperturbation.txt')

run_number_counts(flist=np.array([0.3,0.6]),\
  detections_file='B5/B5_SMGs_run3.0.3arcsec.flux.fovscale1.5.txt',
  effective_area_file='B5/B5_effarea_1.5FWHM.txt',
  band='B5',
  default_simulation='./simulations/default.json',
  simulation_folder='simulations/B5',
  flux_mode='gaussian',
  mode='differential',
  savefile='B5/B5_number_counts_fovscale1.5_differential.txt')

run_number_counts(flist=np.array([0.3,0.6]),\
  detections_file='B5/B5_SMGs_run3.0.3arcsec.flux.fovscale1.5.txt',
  effective_area_file='B5/B5_effarea_1.5FWHM.txt',
  band='B5',
  default_simulation='./simulations/default.json',
  simulation_folder='simulations/B5',
  flux_mode='gaussian',
  mode='differential',
  perturbation=0.1,
  savefile='B5/B5_number_counts_fovscale1.5_differential_withperturbation.txt')



######################
# B6
run_number_counts(#flist=np.logspace(-1.1, 0.3, 8),\
  detections_file='B6/B6_SMGs_run3.0.3arcsec.flux.txt',
  effective_area_file='B6/B6_effarea_1.5FWHM.txt',
  band='B6',
  default_simulation='./simulations/default.json',
  simulation_folder='simulations/B6',
  flux_mode='gaussian',
  mode='cumulative',
  n_points=11,
  savefile='B6/B6_number_counts_fovscale1.5.txt')

run_number_counts(#flist=np.logspace(-1.1, 0.3, 8),\
  detections_file='B6/B6_SMGs_run3.0.3arcsec.flux.txt',
  effective_area_file='B6/B6_effarea_1.5FWHM.txt',
  band='B6',
  default_simulation='./simulations/default.json',
  simulation_folder='simulations/B6',
  flux_mode='gaussian',
  mode='cumulative',
  n_points=11,
  perturbation=0.1,
  savefile='B6/B6_number_counts_fovscale1.5_withperturbation.txt')

run_number_counts(#flist=np.logspace(-1.1, 0.3, 8),\
  detections_file='B6/B6_SMGs_run3.0.3arcsec.flux.txt',
  effective_area_file='B6/B6_effarea_1.5FWHM.txt',
  band='B6',
  default_simulation='./simulations/default.json',
  simulation_folder='simulations/B6',
  flux_mode='gaussian',
  mode='differential',
  n_points=11,
  savefile='B6/B6_number_counts_fovscale1.5_differential.txt')

run_number_counts(#flist=np.logspace(-1.1, 0.3, 8),\
  detections_file='B6/B6_SMGs_run3.0.3arcsec.flux.txt',
  effective_area_file='B6/B6_effarea_1.5FWHM.txt',
  band='B6',
  default_simulation='./simulations/default.json',
  simulation_folder='simulations/B6',
  flux_mode='gaussian',
  mode='differential',
  n_points=11,
  perturbation=0.1,
  savefile='B6/B6_number_counts_fovscale1.5_differential_withperturbation.txt')



######################
# B7
run_number_counts(#flist=np.logspace(-0.7, 0.7, 10),
  detections_file='B7/B7_SMGs_run3.0.3arcsec.flux.fovscale1.5.txt',
  effective_area_file='B7/B7_effarea_1.5FWHM.txt',
  band='B7',
  default_simulation='./simulations/default.json',
  simulation_folder='simulations/B7',
  flux_mode='gaussian',
  n_points=9,
  mode='cumulative',
  savefile='B7/B7_number_counts_fovscale1.5.txt')

run_number_counts(#flist=np.logspace(-0.7, 0.7, 10),
  detections_file='B7/B7_SMGs_run3.0.3arcsec.flux.fovscale1.5.txt',
  effective_area_file='B7/B7_effarea_1.5FWHM.txt',
  band='B7',
  default_simulation='./simulations/default.json',
  simulation_folder='simulations/B7',
  flux_mode='gaussian',
  n_points=9,
  perturbation=0.1,
  mode='cumulative',
  savefile='B7/B7_number_counts_fovscale1.5_perturbation.txt')

run_number_counts(#flist=np.logspace(-0.7, 0.7, 10),
  detections_file='B7/B7_SMGs_run3.0.3arcsec.flux.fovscale1.5.txt',
  effective_area_file='B7/B7_effarea_1.5FWHM.txt',
  band='B7',
  default_simulation='./simulations/default.json',
  simulation_folder='simulations/B7',
  flux_mode='gaussian',
  n_points=9,
  mode='differential',
  savefile='B7/B7_number_counts_fovscale1.5_differential.txt')

run_number_counts(#flist=np.logspace(-0.7, 0.7, 10),
  detections_file='B7/B7_SMGs_run3.0.3arcsec.flux.fovscale1.5.txt',
  effective_area_file='B7/B7_effarea_1.5FWHM.txt',
  band='B7',
  default_simulation='./simulations/default.json',
  simulation_folder='simulations/B7',
  flux_mode='gaussian',
  n_points=9,
  perturbation=0.1,
  mode='differential',
  savefile='B7/B7_number_counts_fovscale1.5_differential_withperturbation.txt')

