function sim = similarity( IMG1, IMG2, similarity_metric, mask )
%SIMILARITY Calculates image similarity
%   sim = SIMILARITY( IMG1, IMG2 ) measures image similarity sim between
%   images IMG1 and IMG2.
% 
%   sim = SIMILARITY( IMG1, IMG2, matching_criteria ) uses similarity
%   metric specified. Available metrics are:
%       * 'monogenic'
%       * 'ssc'
%       * 'mind'
%       * 'nmi'
%       * 'mi'
%       * 'corr' [Default]
% 
%   sim = SIMILARITY( IMG1, IMG2, matching_criteria, mask ) evaluates
%   similarity inside the given binary mask ([all ones]).
% 
%   Alex Bagur, 2020
arguments
    IMG1
    IMG2
    similarity_metric = 'corr'
    mask = ones(size(IMG1))
end
if ~isequal(size(IMG1),size(IMG2))
    error('Images must have equal sizes')
end
mdir = fileparts(mfilename('fullpath'));
switch similarity_metric
    case 'monogenic'
        addpath(genpath('../lsn-monogenic/'))
        [EvenStore,OddStore,LocalPhaseStore1,LocalEnergyStore,OrientationStore_AzimuthAngle,OrientationStore_InclinationAngle, Tstore, OddStoreX, OddStoreY, OddStoreZ] = ...
            MonogenicSignal_MellorBrady(IMG1);
        [EvenStore,OddStore,LocalPhaseStore2,LocalEnergyStore,OrientationStore_AzimuthAngle,OrientationStore_InclinationAngle, Tstore, OddStoreX, OddStoreY, OddStoreZ] = ...
            MonogenicSignal_MellorBrady(IMG2);
        sqd = (LocalPhaseStore1 - LocalPhaseStore2).^2;
        masks = repmat(mask,1,1,1,5);
        assert(isequal(size(sqd),size(masks)),'in similarity:ssc unequal sizes (sqd, masks)')
        sim = 1./mean(sqd(masks));
        
%         set(0,'DefaultFigureWindowStyle','docked')
% %         figure; montage(EvenStore,'DisplayRange',[])
% %         figure; montage(OddStore,'DisplayRange',[])
%         figure; montage(LocalPhaseStore1,'DisplayRange',[])
%         figure; montage(LocalPhaseStore2,'DisplayRange',[])
% %         figure; montage(LocalEnergyStore,'DisplayRange',[])
% %         figure; montage(OrientationStore_AzimuthAngle,'DisplayRange',[])
% %         figure; montage(OrientationStore_InclinationAngle,'DisplayRange',[])
% %         figure; montage(OddStoreX,'DisplayRange',[])
% %         figure; montage(OddStoreY,'DisplayRange',[])
%         figure; montage(sqd,'DisplayRange',[]); title('SQD')
%         figure; montage(masks,'DisplayRange',[]); title('Masks')
    
    case 'ssc'
        addpath([mdir '/similarity_metrics/ssc/'])
        ssc1 = SSC_descriptor(repmat(IMG1,1,1,2));
        ssc2 = SSC_descriptor(repmat(IMG2,1,1,2));
        sqd = (ssc1 - ssc2).^2;
        masks = repmat(mask,1,1,2,12);
        assert(isequal(size(sqd),size(masks)),'in similarity:ssc unequal sizes (sqd, masks)')
        sim = 1./mean(sqd(masks));
%         set(0,'DefaultFigureWindowStyle','docked')
%         figure; montage(reshape(ssc1,size(ssc1,1),size(ssc1,2),size(ssc1,3)*size(ssc1,4)),'DisplayRange',[]); title('SSC IMG1')
%         figure; montage(reshape(ssc2,size(ssc2,1),size(ssc2,2),size(ssc2,3)*size(ssc2,4)),'DisplayRange',[]); title('SSC IMG2')
%         figure; montage(reshape(sqd,size(sqd,1),size(sqd,2),size(sqd,3)*size(sqd,4)),'DisplayRange',[]); title('SQD')
%         figure; montage(reshape(masks,size(masks,1),size(masks,2),size(masks,3)*size(masks,4)),'DisplayRange',[]); title('Masks')
        
    case 'mind'
        addpath([mdir '/similarity_metrics/mind/'])
        mind1 = MIND_descriptor2D(IMG1);
        mind2 = MIND_descriptor2D(IMG2);
        sqd = (mind1 - mind2).^2;
        masks = repmat(mask,1,1,4);
        assert(isequal(size(sqd),size(masks)),'in similarity:mind unequal sizes (sqd, masks)')
        sim = 1./mean(sqd(masks));
%         set(0,'DefaultFigureWindowStyle','docked')
%         figure; montage(mind1,'DisplayRange',[]); title('MIND IMG1')
%         figure; montage(mind2,'DisplayRange',[]); title('MIND IMG2')
%         figure; montage(sqd,'DisplayRange',[]); title('SQD')
%         figure; montage(masks,'DisplayRange',[]); title('Masks')
        
    case 'mind3'
        addpath([mdir '/similarity_metrics/mind/'])
        mind1 = MIND_descriptor(IMG1);
        mind2 = MIND_descriptor(IMG2);
        sqd = (mind1 - mind2).^2;
        masks = repmat(mask,1,1,2,6);
        assert(isequal(size(sqd),size(masks)),'in similarity:mind unequal sizes (sqd, masks)')
        sim = 1./mean(sqd(masks));
%         set(0,'DefaultFigureWindowStyle','docked')
%         figure; montage(reshape(mind1,size(mind1,1),size(mind1,2),size(mind1,3)*size(mind1,4)),'DisplayRange',[]); title('mind1 IMG1')
%         figure; montage(reshape(mind2,size(mind2,1),size(mind2,2),size(mind2,3)*size(mind2,4)),'DisplayRange',[]); title('mind1 IMG2')
%         figure; montage(reshape(sqd,size(sqd,1),size(sqd,2),size(sqd,3)*size(sqd,4)),'DisplayRange',[]); title('SQD')
%         figure; montage(reshape(masks,size(masks,1),size(masks,2),size(masks,3)*size(masks,4)),'DisplayRange',[]); title('Masks')
        
    case 'nmi'
        addpath([mdir '/similarity_metrics/nmi/'])
        sim = nmi(round(IMG1(mask)),round(IMG2(mask)));
        
    case 'mi'
        addpath([mdir '/similarity_metrics/minf/'])
        sim = mi(IMG1(mask),IMG2(mask));
        
    case 'corr'
        sim = corr(IMG1(mask),IMG2(mask)).^2;
        
    otherwise
        error('matching_criteria given not supported. See function help.')
end    
end

