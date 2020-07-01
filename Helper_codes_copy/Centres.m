%% Permet de calculer les coordon�es de tous les centres de pixels de l'implant
function [CentresPixels, nbTotalPixel, rayonPixelTIFF] = Centres(nom_image)
%% INITIALISATION
I=imread(nom_image);
imshow(I); % permet � l'utilisateur de visualiser la photo TIF
title('\bf Cliquez au centre du pixel central puis au centre d''un pixel aligne le plus en peripherie');
X1=ginput(1);%r�cup�re les coordon�es du pixel central [colonne ligne]
Xn=ginput(1);%r�cup�re les coordon�es du pixel p�riph�rique [colonne ligne]
close figure 1;
n=input('Combien de pixels separent les deux selectionnes ? (en comptant les deux) : ');
hex=struct('coord',{},'precedent_ID',{},'fait',{});
hex(1).coord=X1;
hex(1).precedent_ID=2; %utile pour initialisation
hex(1).fait=false;

%d�duire le voisin de l'hexagone central..
Xvoisin=X1-(X1-Xn)/(n-1);
hex(2).coord=Xvoisin;
hex(2).precedent_ID=1;
hex(2).fait=false;

% Nombre total de pixels � trouver
nbTotalPixel=1+6*((n-1)*(n)/2);
%Rayon du pixel
rayonPixelTIFF=(sqrt((X1(1)-Xn(1))^2+(X1(2)-Xn(2))^2))/(2*(n-1));

%% BOUCLE 
while (size(hex,2)<nbTotalPixel)
    hex=creation(hex,rayonPixelTIFF);
end

%% FINALISATION : cr�er la matrice contenant toutes les coordon�es des centres en colonne
CentresPixels=zeros(2,nbTotalPixel);
for i=1:nbTotalPixel
    CentresPixels(:,i)=hex(i).coord;
end
end

