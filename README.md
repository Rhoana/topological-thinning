## Installation

```
git clone https://github.com/Rhoana/topological_thinning.git
cd topological_thinning
conda create -n skeleton_env --file requirements.txt
source activate skeleton_env
cd skeletonization
python setup.py build_ext --inplace
cd ../transforms
python setup.py build_ext --inplace
```

## Meta Files

All datasets are referenced using a meta file. The meta file should have the format meta/{PREFIX}.meta where {PREFIX} is a unique identifier per dataset. The meta file needs to have the following format:

```
# resolution in nm
6x6x30
# segmentation filename
segmentations/SNEMI3D.h5 main
```
where the resolution represents the imaging resolution of the dataset and the segmentation filename is a path to the dataset with the accompanying h5 dataset name.


