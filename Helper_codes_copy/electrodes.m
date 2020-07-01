%% Retourne une structure "electrode" contenant toutes les �lectrodes avec :
% - leur identifiant sur MEA (nombre de 1 � 252)
% - leurs coordonn�es dans la matrice 16x16 (ligne-colonne)
% - leurs coordon�es dans image TIFF
% - leurs coordon�es dans image DMD
function [electrode, distanceBetweenElectrodes, I] = electrodes(ligne_R_TIFF,colonne_R_TIFF,pictureNameCalibrationElectrods)
%% Initialisation :
electrode=struct('ID',{},'coord_mat',{},'coord_TIFF',{},'coord_DMD',{});
I=imread(pictureNameCalibrationElectrods);

I = 0.5*(squeeze(I(:,:,1)) + squeeze(I(:,:,2))); % ROmain's setup
I = 2^16 - I; % ROmain's setup

I = I';

imshow(I); % permet � l'utilisateur de visualiser la photo TIF avec les �lectrodes
title('\bf Cliquez sur les �lectrodes : haut-gauche, bas-gauche, haut-droite');
elecHG_TIFF=ginput(1); %elecHG = electrode Haut Gauche ; [colonne ligne]
elecBG_TIFF=ginput(1); %elecBG = electrode Haut Droite ; [colonne ligne]
elecHD_TIFF=ginput(1); %elecHD = electrode Haut Droite ; [colonne ligne]
close;



%% Construction des coordon�es TIFF pour toutes les �lectrodes
cran_horyzontal=(elecHD_TIFF-elecHG_TIFF)/13;
cran_vertical=(elecBG_TIFF-elecHG_TIFF)/13;
electrodeOrigine_TIFF=elecHG_TIFF-cran_horyzontal-cran_vertical; %On calcul les coordon�es de la vraie �lectrode HG (fictive)
num_elec=1;

load('Electrode_Retina_Grid.mat'); %ElecGrid = matrice contenant les �lectrodes
for i=1:16
    for j=1:16
        electrode(num_elec).coord_mat=[i j];
        electrode(num_elec).coord_TIFF=electrodeOrigine_TIFF+(i-1)*cran_vertical + (j-1)*cran_horyzontal;
        electrode(num_elec).coord_DMD=coordTIFF2DMD_v2(ligne_R_TIFF,colonne_R_TIFF,electrode(num_elec).coord_TIFF', 1);
        electrode(num_elec).ID=ElecGrid(i,j);
        num_elec=num_elec+1;
    end
end

x_h = electrode(1).coord_DMD(1)-electrode(2).coord_DMD(1);
y_h = electrode(1).coord_DMD(2)-electrode(2).coord_DMD(2);


distanceBetweenElectrodes = hypot(x_h, y_h);

% imshow(I);
% for k=1:256
%     hold on;
%     scatter(electrode(k).coord_TIFF(1),electrode(k).coord_TIFF(2));
% end
end

