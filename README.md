# slice2vol
Slice-to-Volume Registration for UK Biobank data

_Needs SPM12 in path_

Example: 

```
slcfnm = './example_images/shmolli_t1.nii.gz';
volfnm = './example_images/volume.nii.gz';

out = slice2vol(slcfnm, volfnm)
```
