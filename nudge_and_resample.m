function slicedvol = nudge_and_resample(mapvol, volvol, translation_shift, interp_method)
% Resample volume (VIBE, WB, seg, etc) after a nudge by translation_shift
% amount using mapvol as reference coordinates
%   translation_shift is a 3x1 vec of shifts in mm for all directions
% 
%   Inspired by GR hepat1ca_PDFF_Sg code
% 
%   Assumes spm12 is in path
% 
%   Alex Bagur, 2021
arguments
	mapvol
	volvol
	translation_shift = zeros(3,1)
	interp_method = 0 % 0: NN, 1: (tri)linear
end

% addpath(genpath('./spm12'))
if numel(translation_shift)~=3
    error('translation_shift needs to be 3x1')
end

if ischar(mapvol)
	mapvol = spm_vol( mapvol );
end
if ischar(volvol)
	volvol = spm_vol( volvol );
end

% interpolate over larger bounds than what spm12 does, considering slice
% thickness information and slice profile
slice_increment = mapvol.mat * [1 1 2 1]' - mapvol.mat * [1 1 1 1]';
slice_thickness = slice_increment(3);
if interp_method == 0
zz=0;
elseif interp_method == 1
zz = linspace(-slice_thickness/2, +slice_thickness/2, 10);
end
w = gausswin(numel(zz));
% figure, plot(w,'-*'), xlabel('z'), title('Slice profile')

[~, XYZ] = spm_read_vols(mapvol);
slicedvol = zeros([mapvol.dim numel(zz)]);
for ii=1:numel(zz)
tmpXYZ = XYZ; % real-world coordinates of mapvol
% later could replace following with more flexible modification of
% world coordinates before subsequent conversion to voxel coordinates
tmpXYZ(1, :) = tmpXYZ(1, :) + translation_shift(1);
tmpXYZ(2, :) = tmpXYZ(2, :) + translation_shift(2);
tmpXYZ(3, :) = tmpXYZ(3, :) + translation_shift(3) + zz(ii);
% Find coordinates in voxel space of vibvol that correspond to
% real-world coordinates tmpXYZ (of parametric map)
tmpxyz = volvol.mat \ [tmpXYZ; ones(1, size(tmpXYZ, 2))];
% Sample vibe at those voxel coordinates
% Use (tri)linear interpolations
slicedvol_linearised = spm_sample_vol(volvol, tmpxyz(1, :), tmpxyz(2, :), tmpxyz(3, :), interp_method);
% Reshape sampled points (linear indexing) to an image with map dims
slicedvol(:,:,:,ii) = reshape(slicedvol_linearised, mapvol.dim(1), mapvol.dim(2)) .* w(ii);
end
if numel(zz)>1
% if interp1 == 0
% slicedvol = median(slicedvols,4);
% elseif interp1 == 1
slicedvol = mean(slicedvol,4);
% end
end

% figure, montage(slicedvols,'DisplayRange',[],'Size',[1 3])
% figure, imshow(slicedvol,[])

%%%%%%
% [~, vXYZ] = spm_read_vols(volvol);
% figure, 
% hold on
% plot3(XYZ(1,:),XYZ(2,:),XYZ(3,:),'.')
% plot3(vXYZ(1,:),vXYZ(2,:),vXYZ(3,:),'.')
% xlabel('X')
% ylabel('Y')
% zlabel('Z')
% legend('slice','vol')
% title('Scans in real-world coordinates')

end
