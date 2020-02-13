cimport cython
cimport numpy as np
from libc.stdint cimport int64_t



import os
import time
import ctypes
from libcpp cimport bool
import numpy as np



from topological_thinning.utilities import dataIO



cdef extern from 'cpp-seg2seg.h':
    void CppDownsampleMapping(const char *prefix, int64_t *segmentation, float input_resolution[3], int64_t output_resolution[3], int64_t input_grid_size[3])




def DownsampleMapping(prefix, segmentation, output_resolution=(80, 80, 80)):
    # everything needs to be long ints to work with c++
    assert (segmentation.dtype == np.int64)

    if not os.path.isdir('skeletons'): os.mkdir('skeletons')
    if not os.path.isdir('skeletons/{}'.format(prefix)): os.mkdir('skeletons/{}'.format(prefix))

    start_time = time.time()

    # convert numpy arrays to c++ format
    cdef np.ndarray[int64_t, ndim=3, mode='c'] cpp_segmentation = np.ascontiguousarray(segmentation, dtype=ctypes.c_int64)
    cdef np.ndarray[float, ndim=1, mode='c'] cpp_input_resolution = np.ascontiguousarray(dataIO.Resolution(prefix), dtype=ctypes.c_float)
    cdef np.ndarray[int64_t, ndim=1, mode='c'] cpp_output_resolution = np.ascontiguousarray(output_resolution, dtype=ctypes.c_int64)
    cdef np.ndarray[int64_t, ndim=1, mode='c'] cpp_input_grid_size = np.ascontiguousarray(segmentation.shape, dtype=ctypes.c_int64)

    # call c++ function
    CppDownsampleMapping(prefix.encode('utf-8'), &(cpp_segmentation[0,0,0]), &(cpp_input_resolution[0]), &(cpp_output_resolution[0]), &(cpp_input_grid_size[0]))

    # free memory
    del cpp_segmentation
    del cpp_input_resolution
    del cpp_output_resolution
    del cpp_input_grid_size

    print ('Downsampled {} to resolution {} in {:0.2f} seconds.'.format(prefix, output_resolution, time.time() - start_time))
