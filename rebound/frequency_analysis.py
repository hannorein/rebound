import ctypes

FREQUENCY_ANALYSIS_TYPES = {"mft": 0, "fmft": 1, "fmft2": 2}
FREQUENCY_ANALYSIS_ERRORS = [  
    (-1, "Frequency analysis error: minfreq must be smaller than maxfreq."),
    (-2, "Frequency analysis error: nfreq must be larger than 0."),
    (-3, "Frequency analysis error: ndata must be power of 2."),
    (-4, "Frequency analysis error: input array is NULL."),
    (-5, "Frequency analysis error: pointer to output array is NULL."),
]

#DLLEXPORT int reb_frequency_analysis(double **output, int nfreq, double minfreq, double maxfreq, enum REB_FREQUENCY_ANALYSIS_TYPE type, double *input, unsigned long ndata);
def frequency_analysis(inp, type=0, nfreq=10, minfreq=-1e-3, maxfreq=1e-3):
    import numpy as np
    ndata = len(inp)//2
    inp_cont = np.ascontiguousarray(inp)
    inp_ptr = inp_cont.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    out = np.zeros(nfreq*3,dtype=np.double)
    out_ptr = out.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    clibrebound.reb_frequency_analysis.restype = ctypes.c_int
    ret = clibrebound.reb_frequency_analysis(out_ptr, ctypes.c_int(nfreq), ctypes.c_double(minfreq), ctypes.c_double(maxfreq), ctypes.c_int(type), inp_ptr, ctypes.c_uint(ndata))
    for value, message in FREQUENCY_ANALYSIS_ERRORS:
        if ret & value:
            raise RuntimeError(message)
    return out

from . import clibrebound
