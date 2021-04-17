# sensi_test.py
import os
import shutil
import multiprocessing as mp
from pathlib import Path
import joblib
from post_processing import write_h5
from h52nc import h52nc

run_path = Path("RUN_Pheno_Module_NEWMASK_.pkz")
pls_path = Path("pls_attrs.csv")

# Experiment - No eCO2 - HISTORICAL and  red precipitation

# new outputs folder
dump_folder = Path("RedPREC_HIST")

with open(run_path, 'rb') as fh:
    init_conditions = joblib.load(fh)

for gridcell in init_conditions:
    gridcell.clean_run(dump_folder, "init_cond")
    gridcell.pr = gridcell.pr - (gridcell.pr * 0.3)
    #gridcell.pr -= (gridcell.pr * 0.9)


def fun(gridcell):
    gridcell.run_caete("19790101", "20161231")

#def fun(gridcell):
#    gridcell.run_caete("19790101", "20161231"), fix_co2="1983")


n_proc = mp.cpu_count() // 2

with mp.Pool(processes=n_proc) as p:
    result = p.map(fun, init_conditions)

to_write = Path(os.path.join(Path("../outputs"), dump_folder)).resolve()
attrs = Path(os.path.join(to_write, Path("pls_attrs.csv"))).resolve()
h5path = Path(os.path.join(to_write, Path('CAETE.h5'))).resolve()
nc_outputs = Path(os.path.join(to_write, Path('nc_outputs')))

shutil.copy(pls_path, attrs)
write_h5(to_write)
h52nc(h5path, nc_outputs)