default_simulation = './simulations_ext/ext_0.0arcsec/ext_0.0arcsec'
simulation_folder='./simulations_main_run3/'
version = 'v5'
outdir = 'number_counts/'

if False:
    # generate the effective area for all the bands
    almacal_catalogue='./check_SMGs_sextractor2.txt'
    dets_all = Table.read(almacal_catalogue, format='ascii')
    for band in ['B3','B4','B5','B6','B7','B8']:
        band_good = dets_all['goodfield_{}'.format(band)] != 0
        run_calculate_effarea_py3(
            imagedir='make_all_good_images/{}'.format(band),
            flux=np.logspace(-2, 1.2, 100),
            objs=dets_all['obj'][band_good],
            band=band,
            fov_scale=1.8,
            central_mask=2.0, #arcsec
            almacal_catalogue=almacal_catalogue,
            savefile='{0}/{0}_effarea_py3_fov1.8_v4.txt'.format(band))

if True:
    for bus in [False, True]:
        # B3 
        run_number_counts_py3(flist=[0.07,0.2],#, 0.6],
          detections_file='B3/B3_catalogue_0.6arcsec_fov1.8_py3_{v}.txt'.format(v=version),
          effective_area_file='B3/B3_effarea_py3_fov1.8_v4.txt',
          band='B3',
          default_simulation=default_simulation,
          simulation_folder=simulation_folder+'B3',
          flux_mode='aperture',
          plot=True,
          rbins=np.arange(0,60,10),
          n_bootstrap=1000,
          boostrap_uncertain_sources=bus,
          savefile=os.path.join(outdir, 'B3_bus-{}'.format(bus)),
         )

        # B4
        run_number_counts_py3(np.logspace(-1.1, -0.3, 4),
          detections_file='B4/B4_catalogue_0.6arcsec_fov1.8_py3_{v}.txt'.format(v=version),
          effective_area_file='B4/B4_effarea_py3_fov1.8_v4.txt',
          band='B4',
          rbins=np.arange(0,41,5),
          default_simulation=default_simulation,
          simulation_folder=simulation_folder+'B4',
          flux_mode='aperture',
          boostrap_uncertain_sources=bus,
          n_bootstrap=1000,
          savefile=os.path.join(outdir, 'B4_bus-{}'.format(bus)),
         )

        # B5 cumu
        run_number_counts_py3(np.logspace(-0.5, -0.1, 2),
          detections_file='B5/B5_catalogue_0.6arcsec_fov1.8_py3_{v}.txt'.format(v=version),
          effective_area_file='B5/B5_effarea_py3_fov1.8_v4.txt',
          band='B5',
          rbins=np.arange(0,31,10),
          default_simulation=default_simulation,
          simulation_folder=simulation_folder+'B5',
          flux_mode='aperture',
          boostrap_uncertain_sources=bus,
          n_bootstrap=1000,
          savefile=os.path.join(outdir, 'B5_bus-{}'.format(bus)),
         )


        # B6 cumulative number counts
        run_number_counts_py3(flist=np.logspace(-1.1, 0.35, 8),
          detections_file='B6/B6_catalogue_0.6arcsec_fov1.8_py3_{v}.txt'.format(v=version),
          effective_area_file='B6/B6_effarea_py3_fov1.8_v4.txt',
          band='B6',
          rbins=np.arange(0,26,4),
          default_simulation=default_simulation,
          simulation_folder=simulation_folder+'B6',
          flux_mode='aperture',
          boostrap_uncertain_sources=bus,
          n_bootstrap=1000, 
          savefile=os.path.join(outdir, 'B6_bus-{}'.format(bus)),
         )

        # B7 cumulative
        run_number_counts_py3(flist=np.logspace(-0.75, 0.78, 8),
          detections_file='B7/B7_catalogue_0.6arcsec_fov1.8_py3_{v}.txt'.format(v=version),
          effective_area_file='B7/B7_effarea_py3_fov1.8_v4.txt',
          band='B7',
          rbins=np.arange(0,21,4),
          default_simulation=default_simulation,
          simulation_folder=simulation_folder+'B7',
          flux_mode='aperture',
          boostrap_uncertain_sources=bus,
          n_bootstrap=1000,
          savefile=os.path.join(outdir, 'B7_bus-{}'.format(bus)),
         )


if False:
    # update the catalogue with accurate flux measurements
    for band in ['B3',]:#['B3','B4','B5','B6','B7','B8']:
        run_update_catalague('{0}/{0}_catalogue_0.6arcsec_fov1.8_py3_v4.txt'.format(band), basedir='make_all_good_images/', band=band,
                              flux_deboosting=True, default_deboosting_simu='simulations/quick_test2_boosting.dat',
                              savefile='{0}/{0}_catalogue_0.6arcsec_fov1.8_py3_v5.txt'.format(band))
    # check the radio images
    for obj in B3_dets2[B3_dets2['type']==1]['obj']:
        print(obj)
        B3_image = 'make_all_good_images/B3/{}/{}_B3_combine.ms.auto.cont.0.6arcsec.image.fits'.format(obj,obj)
        obj_outdir = os.path.join('radio_images/', obj)
        check_radio_images(B3_image, plot=True, savefig=True,
        outdir=obj_outdir, basename=obj, check_NVAS=True, check_FIRST=True, check_NVSS=True)
    

if False:
    """ Notes: the second version (with a suffix of 2) of this code is used to calculate the r_sep for 5 in band 7
    and r_sep=5 in band 6 to have almost equal number of detections in the inner and outer region
    The following code compute the r_sep=5" 
    """
    almacal_catalogue='./check_SMGs_sextractor2.txt'
    # define the seperation radius
    r_sep = 8 # in arcsec
    # seperate the catalogue into inner and outer regions
    B7_cat = Table.read('B7/B7_catalogue_0.6arcsec_fov1.8_py3_v5.txt', format='ascii')
    B7_effarea = Table.read('B7/B7_effarea_py3_fov1.8_v4.txt', format='ascii')
    outer_select = B7_cat['radial_distance'] > r_sep
    inner_select = B7_cat['radial_distance'] < r_sep
    B7_cat_inner = B7_cat[inner_select]
    B7_cat_outer = B7_cat[outer_select]
    B7_cat_inner.write('test_results/B7_catalogue_inner_{}.txt'.format(r_sep), format='ascii', overwrite=True)
    B7_cat_outer.write('test_results/B7_catalogue_outer_{}.txt'.format(r_sep), format='ascii', overwrite=True)
    # calculate the effective area of the inner and outer regions
    dets_all = Table.read(almacal_catalogue, format='ascii')
    for band in ['B7',]:
        band_good = dets_all['goodfield_{}'.format(band)] != 0
        run_calculate_effarea_py3(
            imagedir='make_all_good_images/{}'.format(band),
            flux=np.logspace(-2, 1.2, 100),
            objs=dets_all['obj'][band_good],
            band=band,
            fov_scale=1.8,
            central_mask=r_sep, #arcsec
            almacal_catalogue=almacal_catalogue,
            savefile='{}/{}_effarea_py3_fov1.8_outer_{}.txt'.format('test_results', band, r_sep))
    B7_effarea_of = Table.read('test_results/B7_effarea_py3_fov1.8_outer_{}.txt'.format(r_sep), format='ascii')
    B7_effarea_if = B7_effarea['B7'] - B7_effarea_of['B7']
    B7_effarea_iftab = Table([B7_effarea_of['flux'], B7_effarea_if], names=['flux','B7'])
    B7_effarea_iftab.write('test_results/B7_effarea_py3_fov1.8_inner_{}.txt'.format(r_sep), format='ascii')
    
    # check the number counts in the outer and inner region of each fields,
    # testing the effects of strong gravitional lensing and contamination from the calibartor
    run_number_counts_py3(flist=np.logspace(-0.75, 0.78, 6),
      detections_file='test_results/B7_catalogue_inner_{}.txt'.format(r_sep),
      effective_area_file='test_results/B7_effarea_py3_fov1.8_inner_{}.txt'.format(r_sep),
      band='B7',
      default_simulation='./simulations_ext/ext_0.0arcsec/ext_0.0arcsec',
      simulation_folder='./simulations_main_run3/B7/',
      flux_mode='aperture',
      plot=True,
      perturbation=0.0,
      #mode='cumulative',
      n_points=8,
      n_bootstrap=1000,
      #boostrap_uncertain_sources=True,
      #rbins=np.arange(0,20,4),
      savefile='tests/B7_number_counts_inner_{}.txt'.format(r_sep)
     )
    run_number_counts_py3(flist=np.logspace(-0.75, 0.78, 6),
      detections_file='test_results/B7_catalogue_outer_{}.txt'.format(r_sep),
      effective_area_file='test_results/B7_effarea_py3_fov1.8_outer_{}.txt'.format(r_sep),
      band='B7',
      default_simulation='./simulations_ext/ext_0.0arcsec/ext_0.0arcsec',
      simulation_folder='./simulations_main_run3/B7/',
      flux_mode='aperture',
      plot=True,
      perturbation=0.0,
      #mode='cumulative',
      n_points=8,
      n_bootstrap=1000,
      #boostrap_uncertain_sources=True,
      #rbins=np.arange(0,20,4),
      savefile='tests/B7_number_counts_outer_{}.txt'.format(r_sep)
     )
    
    #
    # boostraping the outer region to get the variaion of number counts under the same area as the inner region
    n_boostraping = 100
    n_fields = 200
    B7_cat = Table.read('B7/B7_catalogue_0.6arcsec_fov1.8_py3_v5.txt', format='ascii')
    dets_all = Table.read(almacal_catalogue, format='ascii')
    idx_select = list(range(len(dets_all)))
    for i in range(n_boostraping):
        dets_select = dets_all[np.random.choice(idx_select, n_fields, replace=False)]
        for band in ['B7',]:
            band_good = dets_select['goodfield_{}'.format(band)] != 0
            run_calculate_effarea_py3(
                imagedir='make_all_good_images/{}'.format(band),
                flux=np.logspace(-2, 1.2, 100),
                objs=dets_select['obj'][band_good],
                band=band,
                fov_scale=1.8,
                central_mask=r_sep, #arcsec
                almacal_catalogue=almacal_catalogue,
                savefile='{}/boostraping/{}_effarea_py3_fov1.8_outer_{}_run{}.txt'.format('test_results', band, r_sep, i))
        cat_selection = []
        for j,obj in enumerate(B7_cat['obj']):
            if obj in dets_select['obj']:
                cat_selection.append(j)
        cat_selected = B7_cat[cat_selection]
        cat_selected.write('test_results/boostraping/B7_catalogue_outer_{}_run{}.txt'.format(r_sep, i), format='ascii', overwrite=True)
        run_number_counts_py3(flist=np.logspace(-0.75, 0.78, 6),
          detections_file='test_results/boostraping/B7_catalogue_outer_{}_run{}.txt'.format(r_sep, i),
          effective_area_file='test_results/boostraping//B7_effarea_py3_fov1.8_outer_{}_run{}.txt'.format(r_sep, i),
          band='B7',
          default_simulation='./simulations_ext/ext_0.0arcsec/ext_0.0arcsec',
          simulation_folder='./simulations_main_run3/B7/',
          flux_mode='aperture',
          plot=False,
          perturbation=0.0,
          n_points=8,
          n_bootstrap=1000,
          #boostrap_uncertain_sources=True,
          #rbins=np.arange(0,20,4),
          savefile='tests/boostraping/B7_number_counts_outer_{}_run{}'.format(r_sep,i)
         )    
