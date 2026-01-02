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
def frequency_analysis(type, input, nfreq=10, minfreq=-1e-3, maxfreq=1e-3):

    ndata = len(input)//2

    out_ptr = ctypes.POINTER(ctypes.c_double)()
    clibrebound.reb_frequency_analysis.restype = ctypes.c_int
    ret = clibrebound.reb_frequency_analysis(ctypes.byref(out_ptr), ctypes.c_int(nfreq), ctypes.c_double(minfreq), ctypes.c_double(maxfreq), type, in_ptr, ctypes.c_uint(ndata))
    for value, message in FREQUENCY_ANALYSIS_ERRORS:
        if ret.value & value:
            raise RuntimeError(message)
