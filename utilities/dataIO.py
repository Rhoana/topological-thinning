import h5py
import struct



import numpy as np



from topological_thinning.data_structures import skeleton_points, meta_data
from topological_thinning.utilities.constants import *



def ReadMetaData(prefix):
    # return the meta data for this prefix
    return meta_data.MetaData(prefix)



def ReadH5File(filename, dataset=None):
    # read the h5py file
    with h5py.File(filename, 'r') as hf:
        # read the first dataset if none given
        if dataset == None: data = np.array(hf[hf.keys()[0]])
        else: data = np.array(hf[dataset])

        return data.astype(np.int64)



def Resolution(prefix):
    # return the resolution for this prefix
    return meta_data.MetaData(prefix).Resolution()



def ReadSegmentationData(prefix):
    filename, dataset = meta_data.MetaData(prefix).SegmentationFilename()

    return ReadH5File(filename, dataset)



def ReadSkeletons(prefix, skeleton_algorithm='thinning', downsample_resolution=(80, 80, 80)):
    # read in all of the skeleton points
    skeleton_filename = 'skeletons/{}/{}-{:03d}x{:03d}x{:03d}-upsample-skeleton.pts'.format(prefix, skeleton_algorithm, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z])
    endpoint_filename = 'skeletons/{}/{}-{:03d}x{:03d}x{:03d}-endpoint-vectors.vec'.format(prefix, skeleton_algorithm, downsample_resolution[IB_X], downsample_resolution[IB_Y], downsample_resolution[IB_Z])

    # read the joints file and the vector file
    with open(skeleton_filename, 'rb') as sfd, open(endpoint_filename, 'rb') as efd:
        skel_zres, skel_yres, skel_xres, skel_max_label, = struct.unpack('qqqq', sfd.read(32))
        end_zres, end_yres, end_xres, end_max_label, = struct.unpack('qqqq', efd.read(32))
        assert (skel_zres == end_zres and skel_yres == end_yres and skel_xres == end_xres and skel_max_label == end_max_label)

        # create an array of skeletons
        skeletons = []
        resolution = Resolution(prefix)
        grid_size = GridSize(prefix)

        for label in range(skel_max_label):
            joints = []
            endpoints = []
            vectors = {}

            nelements, = struct.unpack('q', sfd.read(8))
            for _ in range(nelements):
                index, = struct.unpack('q', sfd.read(8))
                if (index < 0): endpoints.append(-1 * index)
                else: joints.append(index)

            nendpoints, = struct.unpack('q', efd.read(8))
            assert (len(endpoints) == nendpoints)
            for _ in range(nendpoints):
                endpoint, vz, vy, vx, = struct.unpack('qddd', efd.read(32))

                vectors[endpoint] = (vz, vy, vx)

            skeletons.append(skeleton_points.Skeleton(label, joints, endpoints, vectors, resolution, grid_size))

    return skeletons
