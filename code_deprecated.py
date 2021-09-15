def gen_image():
    # the code for detecting the central bright sources and do the point source subtraction
    if check_calibrator:
        imstat_info = imstat(myimagename+'.image')
        if imstat_info['flux'][0] > 0.05: #flux/flux density large than 0.05 Jy
            print("Maybe point source subtraction is failed")
            outfile = os.path.join(outdir, os.path.basename(vis)+".point.cl")
            uvmodelfit(vis=vis, niter=5, comptype="P", sourcepar=[1.0, 0.0, 0.0], varypar=[True,False,False], outfile=outfile)
            # uvmodelfit(vis=vis, niter=5, comptype="G", outfile=outfile,
                       # sourcepar=[1.0, 0.0, 0.0, 1.0, 0.9, 0], 
                       # varypar=[True, False, False, True, True, True])
            ft(vis=vis, complist=outfile)
            uvsub(vis=vis,reverse=False)
            outputvis = os.path.join(outdir, 'updated_data', os.path.basename(vis))
            os.system('mkdir -p {}'.format(os.path.dirname(outputvis)))
            split(vis=vis, outputvis=outputvis)
            if update_raw:
                print('removing', vis)
                rmtables(vis)
                print('copying', outputvis, vis)
                os.system('mv {} {}'.format(outputvis, vis))
                outputvis = vis
            rmtables(myimagename+'.*')
            make_cont_img(vis=outputvis, clean=True, myimagename=myimagename, outdir=outdir, niter=0, only_fits=True, **kwargs)
            return 0

