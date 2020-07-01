
function [outC center_TIFF ] = getFoveaCenter(Floresence_fileName, X_DMD, Y_DMD, DMD)

I = imread(Floresence_fileName);
I = squeeze(I(:,:,1)); % ROmain's setup
I = 2^16 - I; % ROmain's setup

I = I';

imshow(I); % permet � l'utilisateur de visualiser la photo TIF avec les �lectrodes


title('\bf Click the pixel of center of fovea.');
center_TIFF = ginput(1);

hold on; scatter(center_TIFF(1), center_TIFF(2), 30, 'o', 'Filled', 'MarkerFaceColor', [0 1 0]);

outC = coordTIFF2DMD_v2(X_DMD, Y_DMD, center_TIFF', DMD);

% outC = outC';
close all;

end