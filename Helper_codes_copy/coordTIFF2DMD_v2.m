%% Calcul les coordon�es dans rep�re DMD d'un point en coordon�es TIFF
function [coordDMD] = coordTIFF2DMD_v2(ligne_R_TIFF,colonne_R_TIFF,coordTIFF, dmdType)
%correction de la translation par rapport aux origines :

if (dmdType == 2)
    coordDMD=coordTIFF-[colonne_R_TIFF(1);ligne_R_TIFF(1)];

    d12=sqrt((ligne_R_TIFF(1)-ligne_R_TIFF(2))^2+(colonne_R_TIFF(1)-colonne_R_TIFF(2))^2);
    d13=sqrt((ligne_R_TIFF(1)-ligne_R_TIFF(3))^2+(colonne_R_TIFF(1)-colonne_R_TIFF(3))^2);

    %Matrice de rotation :
    teta=asin((colonne_R_TIFF(3)-colonne_R_TIFF(1))/(d13));
    teta=-teta; 
    R=[cos(teta) -sin(teta); sin(teta) cos(teta)]; 
    %ATTENTION : 
    %M rotation valable pour origine = (0;0)... or ici origine = (1;1)... 
    % � utiliser : P'=R*(P-[1;1])+[1;1]
    coordDMD=round(R*(coordDMD-[1;1])+[1;1]);

    % proportion des distances
    coordDMD(1)=(920*coordDMD(1))/(d13*cos(teta)); %nouveau test � faire !
    coordDMD(2)=(480*coordDMD(2))/(d12*cos(teta)); %nouveau test � faire !

    %IL FAUT RAJOUTER une translation selon la position des croix sur le BIN
    %qui est de 300px par 500px.
    coordDMD=coordDMD+[300;500];
    
else
    coordDMD=coordTIFF-[ligne_R_TIFF(1);colonne_R_TIFF(1)];

    d12=sqrt((ligne_R_TIFF(1)-ligne_R_TIFF(2))^2+(colonne_R_TIFF(1)-colonne_R_TIFF(2))^2);
    d13=sqrt((ligne_R_TIFF(1)-ligne_R_TIFF(3))^2+(colonne_R_TIFF(1)-colonne_R_TIFF(3))^2);
    d23=sqrt((ligne_R_TIFF(2)-ligne_R_TIFF(3))^2+(colonne_R_TIFF(2)-colonne_R_TIFF(3))^2);


    %Matrice de rotation :
    teta=asin((ligne_R_TIFF(2)-ligne_R_TIFF(3))/(d23));
    teta=-teta; 
    R=[cos(teta) -sin(teta); sin(teta) cos(teta)]; 
    %ATTENTION : 
    %M rotation valable pour origine = (0;0)... or ici origine = (1;1)... 
    % � utiliser : P'=R*(P-[1;1])+[1;1]
    coordDMD=round(R*(coordDMD-[1;1])+[1;1]);

    % proportion des distances
    coordDMD(1)=(445*coordDMD(1))/(d12*cos(teta)); %nouveau test � faire !
    coordDMD(2)=(235*coordDMD(2))/(d23*cos(teta)); %nouveau test � faire !

    %IL FAUT RAJOUTER une translation selon la position des croix sur le BIN
    %qui est de 300px par 500px.
    coordDMD=coordDMD + [320;150];
end
end

