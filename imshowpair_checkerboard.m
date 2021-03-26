function imshowpair_checkerboard(IMG1,IMG2)
%IMSHOWPAIR_CHECKERBOARD Display registration output using checkerboard
%plot
% Normalise images
IMG1 = (IMG1-mean(IMG1(:)))/std(IMG1(:))+1;
IMG2 = (IMG2-mean(IMG2(:)))/std(IMG2(:))+1;

% Checkerboard
gca;
imagesc(IMG1,[-5 5]);
% colormap(parula)
hold on 
% Save the handle for later use 
h = imagesc(-IMG2,[-5 5]);
hold off
colormap([flipud(gray); parula])
colorbar
I = IMG2;
[M,N] = size(I); 
block_size = 40;
P = ceil(M / block_size); 
Q = ceil(N / block_size); 
alpha = checkerboard(block_size, P, Q) > 0; 
alpha = alpha(1:M, 1:N); 
set(h, 'AlphaData', alpha);
end

