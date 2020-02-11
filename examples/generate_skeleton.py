from topological_thinning.utilities.dataIO import ReadSegmentationData
from topological_thinning.transforms.seg2seg import DownsampleMapping
from topological_thinning.skeletonization.generate_skeletons import TopologicalThinning, FindEndpointVectors



# each dataset is referenced by a unique identifier
# the unique identifer corresponds to a file in meta/{PREFIX}.meta
# that shows locations of various files and gives data attributes
prefix = 'SNEMI3D'



# read the segmentation
seg = ReadSegmentationData(prefix)

# downsample the data for smoother skeletons
DownsampleMapping(prefix, seg)

# call topological thinning function
TopologicalThinning(prefix, seg)
FindEndpointVectors(prefix)