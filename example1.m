clear; close all; clc
slcfnm = './example_images/shmolli_t1.nii.gz';
volfnm = './example_images/volume.nii.gz';
out = slice2vol(slcfnm, volfnm)