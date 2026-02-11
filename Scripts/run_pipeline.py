import os
import astropy.io.ascii
import numpy as np
import glob
import pipeline_functions

#"The code to run it"
data = astropy.io.ascii.read('/home/user/spxmv2/LAGAF/clump_data_infor.csv')
clump_id = np.array(data['source_id'])
loc_7m = np.array(data['7m_loc'])
loc_12m = np.array(data['12m_loc'])
loc_tp = np.array(data['tp_loc'])
loc_folder = np.array(data['loc'])
subloc_tp = np.array(data['tp_sub_loc'])


clump_name = 'C1'
clump_loc = np.where(clump_id == clump_name)[0]
clump_loc_7m = loc_7m[clump_loc]
clump_loc_12m = loc_12m[clump_loc]
clump_loc_tp = loc_tp[clump_loc]
clump_loc_folder = loc_folder[clump_loc]
clump_subloc_tp = subloc_tp[clump_loc]

basic_folder = '/export/home/spxnp1/LAGAF/2022.1.01204.S/'
storage_folder = '/export/home/spxnp1/LAGAF/data_reduction/'
####go to the folder and make a clump folder
os.chdir(storage_folder)
os.makedirs(clump_name+'_combine', exist_ok=True)
new_folder_loc = storage_folder + clump_name + '_combine/'

######build clump ms in 7m
os.chdir(basic_folder + clump_loc_folder[0] + clump_loc_7m[0] + 'calibrated/working/')
file_name = 'calibrated_final.ms'
if os.path.exists(file_name):
    print(f"The file '{file_name}' exists. Skipping this part of the code.")
else:
    print(f"The file '{file_name}' does not exist. Running the code...")
    spw_target()
    print(f"Created '{file_name}'.")

split(vis = file_name, field = clump_name, outputvis = new_folder_loc+f'{clump_name}_7m.ms', datacolumn = 'data')
print('7m clump ms. build and put in the new folder')

######build clump ms in 12m
os.chdir(basic_folder + clump_loc_folder[0]+ clump_loc_12m[0] + 'calibrated/working/')
file_name = 'calibrated_final.ms'
if os.path.exists(file_name):
    print(f"The file '{file_name}' exists. Skipping this part of the code.")
else:
    print(f"The file '{file_name}' does not exist. Running the code...")
    spw_target()
    print(f"Created '{file_name}'.")

split(vis = file_name, field = clump_name, outputvis = new_folder_loc+f'{clump_name}_12m.ms', datacolumn = 'data')
print('12m clump ms. build and put in the new folder')

#####cd to the new clump folder
os.chdir(new_folder_loc)

#####(1) imaging continnum iamge for 12m and 7m
image_continuum(f'{clump_name}_12m.ms/', clump_name, 12)
image_continuum(f'{clump_name}_7m.ms/', clump_name, 7)
image_7m_12m_continuum(f'{clump_name}_7m.ms/', f'{clump_name}_12m.ms/', clump_name)
#####(2) imaging N2H+ spectral image for 12m and 7m
v_range = image_N2H(f'{clump_name}_12m.ms/', clump_name, 12, 0, 0, 0)
image_N2H(f'{clump_name}_7m.ms/', clump_name, 7, v_range[0], v_range[1], v_range[2])
print('7m and 12m N2H+ is done')

#####(3)combine 7m+12m, contsub, imaging spectral line
combine_and_image(f'{clump_name}_7m.ms/', f'{clump_name}_12m.ms/', 0, clump_name, v_range[0], v_range[1], v_range[2])
print('7m+12m N2H+ is done!')


#####(4) imcontsub
run_imcontsub(f'{clump_name}_12m.ms_N2H_aspclean.image')
print('imcontsub for 12m done')
run_imcontsub(f'{clump_name}_7m.ms_N2H_aspclean.image')
print('imcontsub for 7m done')
run_imcontsub('7m+12m.ms_N2H_aspclean.image')
print('imcontsub for 7m+12m done')



#####(5) feathering
tp_loc = basic_folder + clump_loc_folder[0] + clump_loc_tp[0] + f'product/{clump_subloc_tp[0][:-1]}.{clump_name}_sci.spw19.cube.I.sd.fits'
run_feather(tp_loc, '7m+12m.ms_N2H_aspclean_imcontsub.image/', clump_name, 'asp', '7m12m')
run_feather(tp_loc, f'{clump_name}_7m.ms_N2H_aspclean_imcontsub.image/', clump_name, 'asp', '7m')
print('two sets of feather are finished!')

####(6) export data
produce_moment0_fits()
print(f'finish running pipeline for {clump_name}')

####(7) remove all the dirty files except dirty.image
os.system('rm -rf *dirty*gridwt')
os.system('rm -rf *dirty*model')
os.system('rm -rf *dirty*pb')
os.system('rm -rf *dirty*psf')
os.system('rm -rf *dirty*residual')
os.system('rm -rf *dirty*sumwt')
os.system('rm -rf *dirty*weight')
print('remove part of the dirty image files')
