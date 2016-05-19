from necos_extract_vals import extract_values_mask
import sys

'''
Meta script to run necos local metrics
-------------------------------------------------
Can run in two modes:
python necos_meta.py s {subject_id}
python necos_meta.py f {text file containing list of subjects}
'''

##########################
### excluded Subjects ###
##########################
excluded_subjects = ['NECOS008', 'NECOS018','NECOS059']

################################
### Scans to iterate over ###
################################
scan_id_list = ['rest1', 'rest2', 'rest3', 'rest4','rest5', 'rest6']


################################
### Scans groups ###
################################
scan_group_list = ['rest_group_1', 'rest_group_2', 'rest_group_3']

#scan_group_id_list = [1, 2, 3]
#scan_group_id_dict = {1: ['rest1', 'rest2'], 2: ['rest3', 'rest4'], 3: ['rest5', 'rest6']}
#scan_group_name_dict = {1: 'rest_group_1', 2: 'rest_group_2', 3: 'rest_group_3'}

#########################
### Smoothing Kernel ###
########################
fwhm_int = int(4)

################################
### Reho Cluster Size ###
################################
reho_cluster_size = int(27)

##########################
### mask list ###
##########################
mask_list = ['mask_intact_bl2vs1_str_vs_con_le_amy.nii.gz', 
             'mask_intact_bl2vs1_str_vs_con_ri_amy.nii.gz']




# Subjects to iterate over
mode=sys.argv[1]            #takes the arguments in command line, pass them to python sys.argv[0] is the name of the script

if mode == 's':             #if the script is executed in the command line like "python dummy.py s NECOS001", sys.argv[0] --> 'dummy.py', sys.argv[1]--> s
    subject_list=[sys.argv[2]]
elif mode == 'f':
    with open(sys.argv[2], 'r') as f:
        full_list = [line.strip() for line in f]
        # Remove subjects from subjects (list)        
        subject_list = [part for part in full_list if part not in excluded_subjects]

for subject in subject_list:
    
    print 'Running subject '+subject

    working_dir = '/scr/kennedy2/reinelt/NECOS/wd_extract_vals_wf/'+subject+'/' 
    data_dir = '/scr/nil3/reinelt/NECOS/'
    out_dir = '/scr/nil2/reinelt/NECOS/sandkasten/extract_vals'
    
    preproc_dir = '/scr/kennedy2/reinelt/NECOS/MRI/'+subject+'/preprocessed/functional/'
            
    extract_values_mask(working_dir = working_dir, data_dir = data_dir, 
                        out_dir = out_dir, preproc_dir = preproc_dir,
                        subject_id = subject, mask_list = mask_list,
                        scan_group_list = scan_group_list)
                        
