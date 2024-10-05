

import dill
import scanpy as sc
import os
import warnings
warnings.filterwarnings("ignore")
import pandas
import pyranges
# Set stderr to null to avoid strange messages from ray
import sys
_stderr = sys.stderr
null = open(os.devnull,'wb')
import pandas as pd

import sys
args = sys.argv
work_dir=str(args[1])
tmp_dir = '/user1/tanpaku/okada/k-mrkm'

biomart_host = "http://sep2019.archive.ensembl.org/"



from scenicplus.wrappers.run_scenicplus import run_scenicplus

scplus_obj = dill.load(open(os.path.join(work_dir, 'scplus_obj.pkl'),'rb'))

try:
    run_scenicplus(
        scplus_obj = scplus_obj,
        variable = ['GEX_celltype'],
        species = 'hsapiens',
        assembly = 'hg38',
        tf_file = '../../TF_names_v_1.01.txt',
        save_path = os.path.join(work_dir, ''),
        biomart_host = biomart_host,
        upstream = [1000, 150000],
        downstream = [1000, 150000],
        calculate_TF_eGRN_correlation = True,
        calculate_DEGs_DARs = True,
        export_to_loom_file = True,
        export_to_UCSC_file = True,
        path_bedToBigBed = '../../',
        n_cpu = 12,
        _temp_dir = os.path.join(tmp_dir, 'ray_spill'))
except Exception as e:
    #in case of failure, still save the object
    dill.dump(scplus_obj, open(os.path.join(work_dir, 'scplus_obj.pkl'), 'wb'), protocol=-1)
    raise(e)

dill.dump(scplus_obj, open(os.path.join(work_dir, 'scplus_obj.pkl'), 'wb'), protocol=-1)
    

scplus_obj.uns['eRegulon_metadata'].head()


