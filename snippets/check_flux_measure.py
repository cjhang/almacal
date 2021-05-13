os.system('pwd')
mode = 'uv'
vis='./imaging/J1303-5540_combine_uvscale.ms' 
imagefile='./imaging/J1303-5540_combine_uvscale.ms.auto.cont.image.fits'
n=5
repeat=4
snr=(20,20)
fov_scale=1.5
outdir='./test1'
basename = os.path.basename(vis)
uvtaper_scale=None
budget=None
debug=False

# make_cont_img(vis, outdir=outdir, clean=True, niter=1000, suffix='',
              # only_fits=True, uvtaper_scale=uvtaper_scale, pblimit=-0.01,
              # fov_scale=fov_scale, basename=basename)
if True:
    # adding source in uv has not been test for new scheme
    md = msmdtool()
    if not md.open(vis):
        raise ValueError("Failed to open {}".format(vis))
    phasecenter = md.phasecenter()
    mydirection = phasecenter['refer'] +' '+ SkyCoord(phasecenter['m0']['value'], 
                    phasecenter['m1']['value'], unit="rad").to_string('hmsdms')
    freq_mean = np.mean(read_spw(vis))
    myfreq = "{:.2f}GHz".format(freq_mean)
    tb.open(vis + '/ANTENNA')
    antenna_diameter_list = tb.getcol('DISH_DIAMETER')
    tb.close()
    antenna_diameter = np.max(antenna_diameter_list) * u.m
    wavelength = const.c / (freq_mean * u.GHz) # in um
    fov = (fov_scale * 1.22 * wavelength / antenna_diameter * 206265).decompose().value

    im_info = imstat(imagefile)
    sensitivity = im_info['rms'] * 1000 # in mJy/pixel

    fluxrange = np.array(snr) * sensitivity
    known_sources = source_finder(imagefile, fov_scale=fov_scale)
    print(known_sources)
    clname_fullpath = os.path.join(outdir, basename+'.cl')
    for r in np.arange(0.3, 1, 0.1):
        print('radius', r)
        for i in range(repeat):
            print('run {}'.format(i))
            basename_repeat = basename + 'radius{}.run{}'.format(r, i)
            complist_file = os.path.join(outdir, basename_repeat+'.txt')
            if debug:
                print('basename_repeat', basename_repeat)
                print('complist_file', complist_file)
            mycomplist = make_random_source(mydirection, freq=myfreq, 
                    radius=[r*0.45*fov_scale*fov, r*0.45*fov*fov_scale], 
                    # radius=0.9*0.5*fov_scale*fov, # 0.9 is to compensate the optimal imsize 
                    debug=debug, fluxrange=fluxrange, savefile=complist_file, n=n, 
                    sampler=np.random.uniform, sampler_params={}, clname=clname_fullpath,
                    known_sources=known_sources, budget=budget) 
            add_random_sources(vis=vis, mycomplist=mycomplist,
                    outdir=outdir, outname=basename_repeat, debug=debug)
            rmtables(clname_fullpath)

