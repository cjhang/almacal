# this snippets plot the SMG surveys conducted to submillimeter number counts
# Data format: ['wavelength'[um], 'resolution'[arcsec], 'sensitivity[mJy]', 'area'[deg2],]

# 5sigma sensitivity and deepest sensitivity
#

# SCUBA
# 1sigma sensitivity, Geach+2017
# scuba = {'name':'', 'data':[850, , 400,]}
scuba2 = {'name':'S2CLS', 'data':[850., 15, 6,4.73], 'instrument':'SCUBA2', 
          'reference':'Geach et al. (2017)'}

# SPT
SPT = {'name':'SPT-SZ', 'data':[1400,63,11.0,771.2], 'instrument':'SPT', 
        'reference':'Mocanu et al. (2013)'}

# Herschel
H_ATLAS = {'name':'H-ATLAS', 'data':[500.,36.,53,510], 'instrument':'Herschel'}
HERMES = {'name':'HERMES', 'data':[14.0,36,20,12.8], 'instrument':'Herschel'}

# ALMA projects
LESS = {'name':'ALMA-LESS', 'data':[870, 1.5, 2.0, 0.35], 'reference':'Karin et al. (2013)'}
Ono2014 = {'name':'Ono+2014', 'data':[1200, 1.0, 0.1, 0.0008], 'reference':'Ono et al. (2014)'}
Fujimoto2016 = {'name':'Fujimoto+2016', 'data':[1200,0.5,0.055,0.002]}
Oteo2016 = {'name':'Oteo+2016', 'data':[1200,0.3,0.125,0.004]}
Simpson2015 = {'name':'Simpson2015', 'data':[870,0.5,1.5,0.8]}
Munoz2018 = {'name':'Munoz Arancibia2018','data':[1100,0.23,0.3,0.00058]}
HUDF_ALMA = {'name':'ASPECS', 'data':[1200,0.2,0.05,0.00116], 'reference':'Gonzalez-Lopez+2020'}

LABOCA_ECDFS = {'name':'LABOCA-ECDFS', 'data':[870,19,6.,0.25]}

ALMACAL = {'name':'ALMACAL', 'data':[870,0.2,0.03, 0.017]}
ALMARED = {'name':'ALMARED', 'data':[1200,0.7,1.0,510]}

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

plt.style.use('default')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xscale('log')

# circ = patches.Circle((10, 10), 10, alpha=0.8, fc='yellow')
# ax.add_patch(circ)
def scale(a):
    # make the size scal
    return np.log(3600*a)**2/5+2

single_dish = 'C0'
alma = 'C1'
satellite = 'C4'

ax.plot(scuba2['data'][2], scuba2['data'][1], 'o', color=single_dish, markersize=scale(scuba2['data'][-1]), label=scuba2['name'])
ax.text(scuba2['data'][2]/2, scuba2['data'][1], scuba2['name'], horizontalalignment='center', verticalalignment='center')

ax.plot(SPT['data'][2], SPT['data'][1], 'o', color=single_dish, markersize=scale(SPT['data'][-1]), label=SPT['name'])
ax.text(SPT['data'][2], SPT['data'][1], SPT['name'], horizontalalignment='center',verticalalignment='center')

ax.plot(LABOCA_ECDFS['data'][2], LABOCA_ECDFS['data'][1], 'o', color=single_dish, markersize=scale(LABOCA_ECDFS['data'][-1]), label=LABOCA_ECDFS['name'])
ax.text(LABOCA_ECDFS['data'][2], LABOCA_ECDFS['data'][1], LABOCA_ECDFS['name'], horizontalalignment='center',verticalalignment='center')


ax.plot(H_ATLAS['data'][2], H_ATLAS['data'][1], 'o', color=satellite, markersize=scale(H_ATLAS['data'][-1]), label=H_ATLAS['name'])
ax.text(H_ATLAS['data'][2], H_ATLAS['data'][1], H_ATLAS['name'], horizontalalignment='center',verticalalignment='center')

ax.plot(HERMES['data'][2], HERMES['data'][1], 'o', color=satellite, markersize=scale(HERMES['data'][-1]), label=HERMES['name'])
ax.text(HERMES['data'][2]/1.5, HERMES['data'][1], HERMES['name'], horizontalalignment='center',verticalalignment='center')

ax.plot(LESS['data'][2], LESS['data'][1], 'o', color=alma, markersize=scale(LESS['data'][-1]), label=Oteo2016['name'])
ax.text(LESS['data'][2], LESS['data'][1]+2, LESS['name'], horizontalalignment='center')

ax.plot(Ono2014['data'][2], Ono2014['data'][1], 'o', color=alma, markersize=scale(Ono2014['data'][-1]), label=Oteo2016['name'])
ax.text(Ono2014['data'][2], Ono2014['data'][1], Ono2014['name'], horizontalalignment='center')

ax.plot(Fujimoto2016['data'][2], Fujimoto2016['data'][1], 'o', color=alma, markersize=scale(Fujimoto2016['data'][-1]), label=Fujimoto2016['name'])
ax.text(Fujimoto2016['data'][2], Fujimoto2016['data'][1], Fujimoto2016['name'], horizontalalignment='center')

ax.plot(Oteo2016['data'][2], Oteo2016['data'][1], 'o', color=alma, markersize=scale(Oteo2016['data'][-1]), label=Oteo2016['name'])
ax.text(Oteo2016['data'][2], Oteo2016['data'][1], Oteo2016['name'], horizontalalignment='center')

ax.plot(Simpson2015['data'][2], Simpson2015['data'][1], 'o', color=alma, markersize=scale(Simpson2015['data'][-1]), label=Simpson2015['name'])
ax.text(Simpson2015['data'][2], Simpson2015['data'][1], Simpson2015['name'], horizontalalignment='center')

ax.plot(HUDF_ALMA['data'][2], HUDF_ALMA['data'][1], 'o', color=alma, markersize=scale(HUDF_ALMA['data'][-1]), label=HUDF_ALMA['name'])
ax.text(HUDF_ALMA['data'][2], HUDF_ALMA['data'][1], HUDF_ALMA['name'], horizontalalignment='center')

ax.plot(Munoz2018['data'][2], Munoz2018['data'][1], 'o', color=alma, markersize=scale(Munoz2018['data'][-1]), label=Munoz2018['name'])
ax.text(Munoz2018['data'][2], Munoz2018['data'][1], Munoz2018['name'], horizontalalignment='center')

ax.plot(ALMACAL['data'][2], ALMACAL['data'][1], 'o', color=alma, markersize=scale(ALMACAL['data'][-1]), label=ALMACAL['name'])
ax.text(ALMACAL['data'][2], ALMACAL['data'][1], ALMACAL['name'], horizontalalignment='center')

# ax.plot(ALMARED['data'][2], ALMARED['data'][1], 'o', color=alma, markersize=scale(ALMARED['data'][-1]), label=ALMARED['name'],zorder=4)
# ax.text(ALMARED['data'][2], ALMARED['data'][1]+1, ALMARED['name'], horizontalalignment='center', verticalalignment='center',zorder=5)
# ax.text(ALMARED['data'][2], ALMARED['data'][1]-1, "???", horizontalalignment='center', verticalalignment='center',zorder=5)


ax.set_xlabel('Sensitivity [mJy]')
ax.set_ylabel('Resolution [arcsec]')
ax.set_xlim(100, 0.005)
ax.set_ylim(70, -8)

# ax.legend()
plt.show()
fig.savefig('./surveys.png')
