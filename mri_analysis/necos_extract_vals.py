'''
Lifted from calc_metrics; https://github.com/fliem/LeiCA_LIFE/blob/master/metrics/calc_metrics.py#L70
'''
import os
from nipype.pipeline.engine import Node, Workflow
import nipype.interfaces.utility as util
import nipype.interfaces.fsl as fsl
import nipype.interfaces.afni as afni
import nipype.interfaces.io as nio

 
def extract_values_mask(working_dir, data_dir, out_dir, preproc_dir,
                            subject_id, mask_list, scan_group_list):
                        
                              
    afni.base.AFNICommand.set_default_output_type('NIFTI_GZ')
    fsl.FSLCommand.set_default_output_type('NIFTI_GZ')

    ########################
    ### SPECIFY WORKFLOW ###
    ########################

    # Create a local_metrics workflow
    extract_vals_wf = Workflow(name='extract_vals_wf')
    extract_vals_wf.base_dir =  working_dir
    extract_vals_wf.config['execution']['crashdump_dir'] = os.path.join(working_dir, 'crash')
    
    
    #####################
    ### SPECIFY NODES ###
    #####################
    
    ### INPUT ###
    
    # Infosource - a function free node to iterate over the list scans 
    scan_group_infosource = Node(util.IdentityInterface(fields=['scan_group']), 
                               name='scan_group_infosource')
    scan_group_infosource.iterables= [('scan_group', scan_group_list)]             
                                    
                          
    # select Z maps, 2 versions of func (bpf: band pass filtered; no_bpf:no temporal filtering)
    templates={'z_map_le_amy': 'MRI/' + subject_id + '/analysis/sca/nilearn/{scan_group}/left_amygdala/fwhm_6/corr_map_le_amy_z.nii.gz',
               'z_map_ri_amy': 'MRI/' + subject_id + '/analysis/sca/nilearn/{scan_group}/right_amygdala/fwhm_6/corr_map_ri_amy_z.nii.gz',
               }
    select_zmaps = Node(nio.SelectFiles(templates,
                                       base_directory = data_dir),
                       name="select_zmaps")
    extract_vals_wf.connect(scan_group_infosource, 'scan_group', select_zmaps, 'scan_group')
        
    
#    roi_infosource = Node(util.IdentityInterface(fields=['epi_mask']), 
#                               name='roi_infosource')
#    roi_infosource.iterables= [('epi_mask',roi_list)]                                  
    
    
#    # select files, 2 versions of func (bpf: band pass filtered; no_bpf:no temporal filtering)
#    templates={'epi_mask' : 'templates/harvard_oxford/{epi_mask}'
#               }
#    selectmasks = Node(nio.SelectFiles(templates,
#                                       base_directory = project_dir),
#                       name="selectmasks")
                     

    
    ### ANALYSIS NODES ###
    
    def get_vals(le_amy_nii):
                   
        from nipype.interfaces.afni import ROIStats
        import numpy as np
        import os        
           #from nipype.interfaces.fsl import ImageMeants        
            #from nilearn.input_data import NiftiMasker, NiftiMapsMasker
            #import numpy as np
            #import os
            
        ROIStats(mask = '/scr/nil3/reinelt/NECOS/templates/interaction_masks/mask_intact_bl2vs1_str_vs_con_le_amy.nii.gz',
                 quiet = True)
        mean_le_amy = ROIStats.run(le_amy_nii)
                                  
                                  
            # initialize  an empty file & "fill" it with the calculated value, necessary becaus nipype need file types or so... aehm hmm 
        mean_le_amy_file = os.path.abspath('mean_le_amy.txt')
        np.savetxt(mean_le_amy_file, np.atleast_1d(mean_le_amy))         
                      
        return mean_le_amy_file
                              
 

    vals_inside_mask = Node(util.Function(input_names = ['le_amy_nii'],
                                           output_names = ['mean_le_amy_file'],      
                                           function = get_vals),
                    name='sca_prob_seeds')                      
    extract_vals_wf.connect(select_zmaps, 'z_map_le_amy', vals_inside_mask, 'le_amy_nii')     
                        
    
      
#    vals_inside_mask = Node(afni.ROIStats(),
#                            name = 'vals_inside_mask')
#    vals_inside_mask.inputs.mask = '/scr/nil3/reinelt/NECOS/templates/interaction_masks/mask_intact_bl2vs1_str_vs_con_le_amy.nii.gz'
#                                 
#    extract_vals_wf.connect(select_zmaps, 'z_map_le_amy', vals_inside_mask, 'in_file')     
  
    ### OUTPUT ###
        
    # Node to sink data (function DataSink from nipype.interfaces.io)
    datasink = Node(nio.DataSink(base_directory = out_dir,
                                 substitutions=[('_scan_id_', ''),
                                                ('_epi_mask_',''),
                                                ('rest_preprocessed2mni_2mm_maths.nii.gz','temporal_std_2mni_2mm.nii.gz')    
                                                 ]),
                      name="datasink")    
    extract_vals_wf.connect(vals_inside_mask, 'mean_le_amy_file', datasink, 'mean_le_amy')

    
  
                   
    ####################
    ### RUN WORKFLOW ###
    ####################    
    
    #extract_vals_wf.write_graph(dotfilename='sca', graph2use='flat', format='pdf')
    #extract_vals_wf.run()rest_preprocessed2mn
    #extract_vals_wf.run(plugin='CondorDAGMan')    
    extract_vals_wf.run(plugin='MultiProc', plugin_args={'n_procs' : 30})

