%% Retourne une structure "electrode" contenant toutes les �lectrodes avec :
% - leur identifiant sur MEA (nombre de 1 � 252)
% - leurs coordonn�es dans la matrice 16x16 (ligne-colonne)
% - leurs coordon�es dans image TIFF
% - leurs coordon�es dans image DMD
function [electrode] = electrodes_packed(elecDMD, ligne_R_TIFF,colonne_R_TIFF, distanceBetweenElectrodes)
%% Initialisation :
electrode=struct('ID',{},'coord_mat',{},'coord_TIFF',{},'coord_DMD',{});

elecHG_DMD = elecDMD(1,:);
elecBG_DMD = elecDMD(2,:);
elecHD_DMD = elecDMD(3,:);

%% Construction des coordon�es TIFF pour toutes les �lectrodes

numHorizontal = round((elecHD_DMD-elecHG_DMD)/distanceBetweenElectrodes)+1;
numHorizontal = numHorizontal(1);

numVertical= round((elecBG_DMD-elecHG_DMD)/distanceBetweenElectrodes)+1;
numVertical = numVertical(2);


electrodeOrigine_DMD = elecHG_DMD - [distanceBetweenElectrodes, distanceBetweenElectrodes]; %On calcul les coordon�es de la vraie �lectrode HG (fictive)
num_elec=1;

load('Electrode_Retina_Grid.mat'); %ElecGrid = matrice contenant les �lectrodes

for i=1:numHorizontal
    for j=1:numVertical
        electrode(num_elec).coord_mat=[i j];
        
        electrode(num_elec).coord_DMD = electrodeOrigine_DMD + [(i-1)*distanceBetweenElectrodes, (j-1)*distanceBetweenElectrodes];
        
        electrode(num_elec).coord_TIFF = coordDMD2TIFF(ligne_R_TIFF,colonne_R_TIFF,electrode(num_elec).coord_DMD', 2);

        %         electrode(num_elec).coord_DMD=coordTIFF2DMD_v2(ligne_R_TIFF,colonne_R_TIFF,electrode(num_elec).coord_TIFF', 2);
%         electrode(num_elec).ID=ElecGrid(i,j);
        num_elec=num_elec+1;
    end
end

x_h = electrode(1).coord_DMD(1)-electrode(2).coord_DMD(1);
y_h = electrode(1).coord_DMD(2)-electrode(2).coord_DMD(2);


% distanceBetweenElectrodes = hypot(x_h, y_h)

% imshow(I);
% for k=1:256
%     hold on;
%     scatter(electrode(k).coord_TIFF(1),electrode(k).coord_TIFF(2));
% end
end

