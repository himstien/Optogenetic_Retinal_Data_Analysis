%% Retourne les 6 coord des centres autour du pixel indiqué par "ID_hexa_centre"
function [M_6_hex] = autour(ID_hexa_centre, hex)
%M_6_hex est une matrice
% ligne 1: abscices des centres
% ligne 2: ordonnée des centres
% 6 colonnes pour 6 vecteurs
angle=60*pi/180;
M=[cos(angle),-sin(angle);sin(angle),cos(angle)];
centre=hex(ID_hexa_centre).coord;
if size(centre,1)==1 %si centre n'est pas en colonne, le remettre en colonne
    centre=centre';
end
peripherique=hex(hex(ID_hexa_centre).precedent_ID).coord;
if size(peripherique,1)==1
    peripherique=peripherique';
end
for i=1:6
    M_6_hex(:,i)=M^i*(peripherique-centre)+centre;
end
end


