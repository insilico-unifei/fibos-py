import ctypes as ct
import numpy as np
import os
import main75
import ds75
import surfcal76

def call_main(iresf, iresl, maxres, maxat, meth):
    resnum = (ct.c_int*maxres)()
    x = (ct.c_double*maxat)()
    y = (ct.c_double*maxat)()
    z = (ct.c_double*maxat)()
    natm = np.array([1], dtype=np.int_)
    main75.main(resnum,natm,x,y,z,iresf,iresl)
    for ires in range(1, iresl+1):
        main75.main_intermediate(x,y,z,ires,resnum,natm[0])
        main75.main_intermediate01(x,y,z,ires,resnum,natm[0])
        ds75.runsims(meth)
        surfcal76.surfcal()
    os.rename("file.srf","prot.srf")
