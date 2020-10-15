def gen_image(obj, band=None):
    for obs in os.listdir(obj):
        if band is not None:
            band_match = re.compile('_(?P<band>B\d{1,2})')
            if band_match.search(obs):
                obs_band = band_match.search(obs).groupdict()['band']
                if obs_band != band:
                    continue
        try:
            make_cont_img(vis=obj+'/'+obs, dirty_image=True)
        except:
            print("Error in imaging {}".format(obj))
        basename = os.path.basename(obs)
        myimagename = basename + '.cont.auto'
        exportfits(imagename=myimagename+'.image', fitsimage=myimagename+'.fits')
        rmtables(tablenames=myimagename+'.*')


if __name__ == '__main__':
    obj = '../1924-2914'
    gen_image(obj)
