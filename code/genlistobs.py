import os
import re

base_dir = '/science-ALMACAL/data'
p_obj = re.compile('J\d+[+-]\d')
p_obs = re.compile('uid___')

out_dir = '/scratch-ssd/jchen/listobs'

bad_obs = ['uid___A002_Xb21481_Xda94.ms.split.cal.J1743-0350_B6','uid___A002_X87f18c_X8ff.ms.split.cal.J1733-1304_B6','uid___A002_Xa2ea64_X18a0.ms.split.cal.J1733-1304_B7','uid___A002_X7d76cc_X1698.ms.split.cal.J1733-1304_B6', 'uid___A002_X7d2ec1_Xa80.ms.split.cal.J1733-1304_B7','uid___A002_X9aca45_X1d8b.ms.split.cal.J1733-1304_B3','uid___A002_X88063e_X190.ms.split.cal.J1733-1304_B6','uid___A002_Xb24884_X27e7.ms.split.cal.J1733-1304_B3','uid___A002_Xb73164_X5eb9.ms.split.cal.J1733-1304_B7','uid___A002_X8a951a_X355.ms.split.cal.J1733-1304_B6','uid___A002_X95e355_X220a.ms.split.cal.J1733-1304_B3','uid___A002_Xa309c2_X254f.ms.split.cal.J1733-1304_B6,uid___A002_Xbab09c_X4cef.ms.split.cal.J1733-1304_B3','uid___A002_Xa309c2_X254f.ms.split.cal.J1733-1304_B6','uid___A002_Xbab09c_X4cef.ms.split.cal.J1733-1304_B3','uid___A002_Xb24884_X1f46.ms.split.cal.J1733-1304_B3','uid___A002_X87f18c_Xed4.ms.split.cal.J1733-1304_B6','uid___A002_Xa2d681_X1ad9.ms.split.cal.J1733-1304_B6','uid___A002_X837c61_Xffd.ms.split.cal.J1733-1304_B6','uid___A002_X71cec2_X40.ms.split.cal.J1733-1304_B6','uid___A002_X8970e9_X2181.ms.split.cal.J1733-1304_B6','uid___A002_X8a951a_X69b.ms.split.cal.J1733-130_B6','uid___A002_Xb03fa0_X3715.ms.split.cal.J1733-1304_B6','uid___A002_Xaf985b_Xf3f.ms.split.cal.J1651+0129_B6','uid___A002_X9db5ce_X1dd1.ms.split.cal.J1550+0527_B7','uid___A002_Xa4b3e9_X828.ms.split.cal.J1550+054_B7','uid___A002_Xa4b3e9_X828.ms.split.cal.J1550+0527_B7','uid___A002_Xbb9cd7_Xae3b.ms.split.cal.J1058+0133_B6','uid___A002_X71a45c_X1a41.ms.split.cal.J1058+0133_B6','uid___A002_Xbbadbe_X88ec.ms.split.cal.J1058+0133_B6','uid___A002_Xbb85b6_X8556.ms.split.cal.J1058+0133_B6','uid___A002_Xbbadbe_X8f99.ms.split.cal.J1058+0133_B6','uid___A002_Xbabd27_X4714.ms.split.cal.J1058+0133_B6','uid___A002_Xb6a8c1_X47f7.ms.split.cal.J0431+2037_B7','uid___A002_Xa2300a_X1bbf.ms.split.cal.J1751+0939_B6','uid___A002_Xb148a2_X4fc6.ms.split.cal.J1751+0939_B6','uid___A002_Xb1586e_X37e3.ms.split.cal.J1751+0939_B6','uid___A002_Xb1cc39_X2f8d.ms.split.cal.J1751+0939_B6','uid___A002_Xb0be8b_X5366.ms.split.cal.J1751+0939_B6','uid___A002_Xa24618_X1c44.ms.split.cal.J1751+0939_B6','uid___A002_X822d50_X157b.ms.split.cal.J0006-0623_B7']



for obj in os.listdir(base_dir):
    if p_obj.match(obj):
        os.system('mkdir -p {}/{}'.format(out_dir, obj))
#if True:
#    if True:
#        obj = 'J1751+0939'
        for obs in os.listdir(base_dir +'/'+ obj):
            print(obs)
            if obs in bad_obs:
                continue
            if p_obs.match(obs):
                obs_filename = base_dir +'/'+ obj+'/'+obs
                listobs_filename = out_dir +'/'+  obj +'/'+ obs
                try:
                    listobs(vis=obs_filename, listfile=listobs_filename+'.listobs.txt')
                except:
                    print("Error: in", obs_filename)
                #print(obs_filename)
