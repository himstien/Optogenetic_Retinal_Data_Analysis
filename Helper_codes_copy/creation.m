%% On créé la structure hex qui contiendra toutes les coordonnées de tous les pixels
function [hex] = creation(hex,rayonPixel)
%prendre plus petit ID avec fait=0 --> hexa_centre
i=1;
while (hex(i).fait==true)
    i=i+1;
end
ID_hexa_centre=i;
[M_6_hex]=autour(ID_hexa_centre,hex);
%si M_6_hex(1->6) pas déjà dans structure hex --> rajouter dans structure
for k=1:6
    j=1;
    trouve=false;
    while (j<=size(hex,2) && trouve==false)
        distance=sqrt((M_6_hex(1,k)-hex(j).coord(1))^2+(M_6_hex(2,k)-hex(j).coord(2))^2); %calcul distance entre le kième point et pixels existant
        if distance<rayonPixel
            trouve=true;
        end
        j=j+1;
    end
    if trouve==false %Si ce centre n'existe pas déjà, on l'ajoute à la structure
        hex(end+1).coord=M_6_hex(:,k);
        hex(end).precedent_ID=ID_hexa_centre;
        hex(end).fait=false;
    end
end
hex(ID_hexa_centre).fait=true; %on a bien testé/créé tous les centres autour de ce pixel 
end


