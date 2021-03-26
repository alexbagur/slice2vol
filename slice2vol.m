function out_reg = slice2vol(slcfnm,volfnm,similarity_metric,method,doFIGURES)
%SLICE2VOL Slice-to-volume registration.
% Affine in xy plane, Exhaustive search in z.
% Images are nifti filenames.
% 
% Fixed image is slice (since it can be quantitative)
% 
% Assumes spm12 is in path.
% 
%   out.tform can be used to resample a segmentation aligned with volume
%       posteriorly to obtain ROI on the slice.
%   Provide example
% 
%   method
%       * ZExhaustive exhaustive search over Z, no XY registration
%       * XYAffineZExhaustive affine in XY at z=0 then applies to all z,
%           choose z with highest similarity
%       * (No other methods supported at the moment)
% 
%   In order to resample the volume corresponding segmentation, use the
%   following:
% 
%       mask2d = nudge_and_resample(mapfnm, sgifnm, [0,0,out_reg.z_shift_opt]);
%       mask2d = imwarp(mask2d,affine2d(out_reg.xy_tform.T),'OutputView',imref2d(out_reg.slcvol_dim),'interp','nearest');
%       mask2d = round(rot90(mask2d));
% 
% Alex Bagur, 2021
arguments
    slcfnm
    volfnm
    similarity_metric = 'ssc'
    method = 'XYAffineZExhaustive'
    doFIGURES = true
end

if isempty(method)
    method = 'XYAffineZExhaustive';
end
    
switch method
    case 'ZExhaustive'
        doXYREG = 0;
    case 'XYAffineZExhaustive'
        doXYREG = 1;
    otherwise
        error('Registration method not supported yet')
end

disp(['Running slice-to-volume registration with similarity ' upper(similarity_metric) '...'])
if ischar(slcfnm)
    disp(['Slice: ' slcfnm])
    disp(['Volume: ' volfnm])
end
out_reg.slcfnm = slcfnm;
out_reg.volfnm = volfnm;
out_reg.method = method;

%%%%%%%% Load in nifti data
addpath(genpath('./spm12'))
% Read in parametric map
slcvol = spm_vol( slcfnm );
slcvol.dat(isnan(slcvol.dat))=0;
slcvoldat = slcvol.dat;

% Read in 3d volume
volvol = spm_vol( volfnm );

% Create Volume slice candidates
% Compute exhaustive search range based on volume slab height
% Maximum of 50 mm either side (-50mm to +50mm)
vol_increment = volvol.mat * [1 1 2 1]' - volvol.mat * [1 1 1 1]';
vol_thickness_mm = min((vol_increment(3) * volvol.dim(3))/4, 50);
z_shifts_mm = -vol_thickness_mm:1:+vol_thickness_mm;
disp(['z_shifts_mm = ' num2str(z_shifts_mm(1)) 'mm to ' num2str(z_shifts_mm(end)) 'mm relative to DICOM header reference'])
disp('Computing Volume candidate slices ...')
t=tic;
volarr = cell(1,numel(z_shifts_mm));
for n=1:numel(volarr)
    volarr{n} = nudge_and_resample(slcvol, volvol, [0,0,z_shifts_mm(n)], 1);
end
toc(t);

if doXYREG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Affine transformation in XY
%   Do it once at z=0 then apply to all subsequent z
%   reg params that do not change 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Calculating initial affine XY transformation for volume slice at z=0 ...')
disp('This transformation will be used for subsequent z values')
rng('default')
% tformarr = cell(1, N);
volarr_0 = nudge_and_resample(slcvol, volvol, zeros(3,1), 1);
[optimizer, metric] = imregconfig('multimodal');
optimizer.InitialRadius = 0.0002;
optimizer.Epsilon = 1.5e-4;
optimizer.GrowthFactor = 1.01;
optimizer.MaximumIterations = 50;
tform_0 = imregtform(volarr_0, slcvoldat, 'affine', optimizer, metric);

if doFIGURES
% Plot XY alignment after affine
volarr_0_warped = imwarp(volarr_0,tform_0,'OutputView',imref2d(slcvol.dim));
figure('Position',[32 283 1416 589]);
nexttile; imshowpair_checkerboard(rot90(slcvoldat),rot90(volarr_0)), 
title('Before XY registration')
nexttile; imshowpair_checkerboard(rot90(slcvoldat),rot90(volarr_0_warped))
title('After XY registration')
drawnow
end

% optimizer.GrowthFactor = 1.005; % lowered as less work needed with good starting point
% Initial tform_0 works well enough in most cases with no refinement
for n=1:numel(volarr)
%     tformarr{n} = imregtform(vibarr{n}, slcvoldat, 'affine', optimizer, metric, 'InitialTransformation', tform_0);
%     tformarr{n} = tform_0;
%     volarr{n} = imwarp(volarr{n},tformarr{n},'OutputView',imref2d(slcvol.dim));
    volarr{n} = imwarp(volarr{n},tform_0,'OutputView',imref2d(slcvol.dim));
end
out_reg.xy_tform = tform_0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Rigid transformations (shifts) along Z
%   Exhaustive search over predefined range of z
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=tic;
cost_mask = calculate_cost_evaluation_mask( volarr );
[z_cost_function, min_I] = optimisation_loop( slcvoldat, volarr, similarity_metric, cost_mask );
% opt_translation_shifts = fminsearch(@(x) slice2vol_cost_fun(x,volvol,slcvol,similarity_metric,cost_mask),[0 0 0])
toc(t)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save output
out_reg.similarity_metric = similarity_metric;
out_reg.slcvol_dim = slcvol.dim;
out_reg.z_shifts_mm = [num2str(z_shifts_mm(1)) ':' num2str(z_shifts_mm(2)-z_shifts_mm(1)) ':' num2str(z_shifts_mm(end))];
out_reg.min_I = min_I;
out_reg.z_shift_opt = z_shifts_mm(min_I);
out_reg.z_cost_function = z_cost_function;
out_reg.slcvoldat = slcvoldat;
out_reg.volarr = volarr;
out_reg.volarr_0 = volarr_0;
out_reg.vol_aligned_slice = volarr{min_I};
out_reg.cost_mask = cost_mask;
out_reg.cost_mask_indices = find(rot90(cost_mask));

if doFIGURES
figure
montage(volarr,'DisplayRange',[])
title('Volume slice candidates along z')

figure
imshow(cost_mask)

figure
plot(eval(out_reg.z_shifts_mm),out_reg.z_cost_function,'LineWidth',2)
hold on
xline(out_reg.z_shift_opt,'--');
legend('cost','optimum')
xlabel('z translations (mm)')
ylabel('Cost function to minimise (a.u.)')
title('Registration cost for z shifts from DICOM header reference')

figure('Position',[529 93 1096 467])
subplot(121)
imagesc(out_reg.slcvoldat)
colormap(gray)
title('Slice')
subplot(122)
imagesc(out_reg.vol_aligned_slice)
colormap(gray)
title(['Vol chosen slice at z=' num2str(out_reg.z_shift_opt) 'mm'])
drawnow
end

end

function [cost, min_I] = optimisation_loop( slcvoldat, volarr, similarity_metric, cost_mask )
%%%%%%%% OPTIMISATION: Maximise image similarity [Minimise 1/(image similarity)]
disp('Running optimisation ...')
disp('Calculating image similarity for all volume slice candidates ...')
N = numel(volarr);
cost = nan(1,N);
for n=1:N
    cost(n) = 1/similarity(slcvoldat, volarr{n}, similarity_metric, cost_mask);
end
[~,min_I]=min(cost,[],'omitnan');
end

function cost_evaluation_mask = calculate_cost_evaluation_mask( volarr )
% Cost mask is where image similarity will be evaluated
% Cost mask is the same for all candidates
cost_evaluation_mask = logical(ones(size(volarr{1})));
for n=1:numel(volarr)
    % Compute body mask from neighbourhood
%     outside_patch = volarr{n}(1:40,1:40);
%     noise_thr = mean(outside_patch(:),'omitnan')+2*std(outside_patch(:),'omitnan');
%     mask_n = imfill(volarr{n}>noise_thr,'holes');
    mask_n = volarr{n}>10;
    if ~any(mask_n(:))
        continue
    end
    cost_evaluation_mask = cost_evaluation_mask & mask_n;
end
% cost_evaluation_mask = bwareafilt(cost_mask, 1);
% figure; imshow(cost_evaluation_mask); drawnow
end
