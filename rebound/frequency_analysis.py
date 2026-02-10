import ctypes

FREQUENCY_ANALYSIS_TYPES = {"mft": 0, "fmft": 1, "fmft2": 2}
FREQUENCY_ANALYSIS_ERRORS = [  
    (-1, "Frequency analysis error: minfreq must be smaller than maxfreq."),
    (-2, "Frequency analysis error: nfreq must be larger than 0."),
    (-3, "Frequency analysis error: ndata must be power of 2."),
    (-4, "Frequency analysis error: input array is NULL."),
    (-5, "Frequency analysis error: pointer to output array is NULL."),
]

def frequency_analysis(inp, type=0, nfreq=10, minfreq=-1e-3, maxfreq=1e-3):
    """Performs a frequency analysis on the timeseries data inp. Returns
    dominant modes (frequency, amplitude, phase). 
    
    Arguments
    ---------
    inp: numpy.array
        Input data in the order x[0], y[0], x[1], y[1], ... where
        x and y are the real and imaginary components of the signal.
    type: string
        Determines the type of the frequency analysis:
        "mft" = Modified Fourier Transform. See Laskar (1988).
                https://ui.adsabs.harvard.edu/abs/1988A%26A...198..341L/abstract
        "fmft" = Frequency Modified Transform. See Sidlichovsky and Nesvorny (1996).
                https://ui.adsabs.harvard.edu/abs/1996CeMDA..65..137S/abstract
        "fmft2" = Frequency Modified Transform with additional corrections. Most
                accurate but also slowest. See Sidlichovsky and Nesvorny (1996).
    nfreq: Int
        The number of frequencies to find.
    minfreq: Float
        The minimum frequency to consider. Units are [radians/datasep] where
        datasep is the timeinterval between sampling points.
    maxfreq: Float
        The maximim frequency to consider. Units are [radians/datasep] where
        datasep is the timeinterval between sampling points.

    Returns
    -------
    A list of lists, containing the frequencies, amplitudes, and phases of the 
    nfreq most dominant modes.
    The output units for the frequencies are [radians/datasep] where datasep is
    the timeinterval between sampling points. The output units for the phases
    are radians. 
    """
    import numpy as np
    if not isinstance(inp, np.ndarray):
        raise ValueError("Input array must be a numpy array")
    if inp.dtype != np.float64:
        raise ValueError("Input array must be have datatype np.float64")

    if isinstance(type, int):
        type = ctypes.c_int(type)
    elif isinstance(type, basestring):
        type = type.lower()
        if type in FREQUENCY_ANALYSIS_TYPES:
            type = FREQUENCY_ANALYSIS_TYPES[type]
        else:
            raise ValueError("Frequency Analysis Type not found.")
    else:
        raise ValueError("Frequency Analysis Type not found.")

    ndata = len(inp)//2
    inp_cont = np.ascontiguousarray(inp)
    inp_ptr = inp_cont.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    out = np.zeros(nfreq*3,dtype=np.double)
    out_ptr = out.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    clibrebound.reb_frequency_analysis.restype = ctypes.c_int
    ret = clibrebound.reb_frequency_analysis(out_ptr, ctypes.c_int(nfreq), ctypes.c_double(minfreq), ctypes.c_double(maxfreq), type, inp_ptr, ctypes.c_uint(ndata))
    for value, message in FREQUENCY_ANALYSIS_ERRORS:
        if ret & value:
            raise RuntimeError(message)
    return np.split(out,3)

from . import clibrebound
