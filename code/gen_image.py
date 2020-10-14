def gen_image(obj):
    for obs in os.listdir(obj):
        make_cont_img(vis=obs, dirty_image=True)
        basename = os.path.basename(obs)
        myimagename = basename + '.cont.auto'
        exportfits(imagename=myimagename+'.image', fitsimage=myimagename+'.fits')
        rmtables(tablenames=myimagename+'.*')


if __name__ == '__main__':
    obj = ''
    gen_image(obj)
