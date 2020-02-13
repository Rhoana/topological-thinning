import os
import time
import struct


cimport cython
cimport numpy as np
from libcpp cimport bool
import ctypes
import numpy as np
from libc.stdint cimport int64_t



from topological_thinning.utilities import dataIO
from topological_thinning.utilities.constants import *



cdef extern from 'cpp-generate_skeletons.h':
    void CppTopologicalThinning(const char *prefix, int64_t skeleton_resolution[3], const char *lookup_table_directory)
    void CppFindEndpointVectors(const char *prefix, int64_t skeleton_resolution[3], float output_resolution[3])
    void CppApplyUpsampleOperation(const char *prefix, int64_t *input_segmentation, int64_t skeleton_resolution[3], float output_resolution[3])



# generate skeletons for this volume
def TopologicalThinning(prefix, input_segmentation, skeleton_resolution=(80, 80, 80)):
    # everything needs to be long ints to work with c++
    assert (input_segmentation.dtype == np.int64)

    start_time = time.time()

    # convert the numpy arrays to c++
    cdef np.ndarray[int64_t, ndim=1, mode='c'] cpp_skeleton_resolution = np.ascontiguousarray(skeleton_resolution, dtype=ctypes.c_int64)
    lut_directory = os.path.dirname(__file__)

    # call the topological skeleton algorithm
    CppTopologicalThinning(prefix.encode('utf-8'), &(cpp_skeleton_resolution[0]), lut_directory.encode('utf-8'))

    # call the upsampling operation
    cdef np.ndarray[int64_t, ndim=3, mode='c'] cpp_input_segmentation = np.ascontiguousarray(input_segmentation, dtype=ctypes.c_int64)
    cdef np.ndarray[float, ndim=1, mode='c'] cpp_output_resolution = np.ascontiguousarray(dataIO.Resolution(prefix), dtype=ctypes.c_float)

    CppApplyUpsampleOperation(prefix.encode('utf-8'), &(cpp_input_segmentation[0,0,0]), &(cpp_skeleton_resolution[0]), &(cpp_output_resolution[0]))

    print ('Generated skeletons for {} in {:0.2f} seconds.'.format(prefix, time.time() - start_time))



# find endpoint vectors for this skeleton
def FindEndpointVectors(prefix, skeleton_resolution=(80, 80, 80)):
    start_time = time.time()

    # convert to numpy array for c++ call
    cdef np.ndarray[int64_t, ndim=1, mode='c'] cpp_skeleton_resolution = np.ascontiguousarray(skeleton_resolution, dtype=ctypes.c_int64)
    cdef np.ndarray[float, ndim=1, mode='c'] cpp_output_resolution = np.ascontiguousarray(dataIO.Resolution(prefix), dtype=ctypes.c_float)

    CppFindEndpointVectors(prefix.encode('utf-8'), &(cpp_skeleton_resolution[0]), &(cpp_output_resolution[0]))

    print ('Found endpoint vectors for {} in {:0.2f} seconds.'.format(prefix, time.time() - start_time))
