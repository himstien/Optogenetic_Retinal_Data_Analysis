%% Calcul les coordon�es dans rep�re DMD d'un point en coordon�es TIFF
function [coordTIFF] = coordDMD2TIFF_rotBar(ligne_R_TIFF,colonne_R_TIFF, coordDMD, dmdType)
%correction de la translation par rapport aux origines :

    if (dmdType == 1)
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
        
       
        d12=sqrt((ligne_R_TIFF(1)-ligne_R_TIFF(2))^2+(colonne_R_TIFF(1)-colonne_R_TIFF(2))^2);
        d23=sqrt((ligne_R_TIFF(2)-ligne_R_TIFF(3))^2+(colonne_R_TIFF(2)-colonne_R_TIFF(3))^2);


        %Matrice de rotation :
        teta=asin((ligne_R_TIFF(2)-ligne_R_TIFF(3))/(d23));
%         teta=-teta; 
        R=[cos(teta) -sin(teta); sin(teta) cos(teta)]; 
        
        for i=1:size(coordDMD, 2)

            %IL FAUT RAJOUTER une translation selon la position des croix sur le BIN
            %qui est de 300px par 500px.
            coordDMD(:,i) = coordDMD(:,i) - [320;150];

            % proportion des distances
            coordDMD(1, i)=(coordDMD(1, i)/445)*(d12*cos(teta)); %nouveau test � faire !
            coordDMD(2, i)=(coordDMD(2, i)/235)*(d23*cos(teta)); %nouveau test � faire !

        %ATTENTION : 
        %M rotation valable pour origine = (0;0)... or ici origine = (1;1)... 
        % � utiliser : P'=R*(P-[1;1])+[1;1]

            coordDMD(:, i) = round(R*(coordDMD(:,i)-[1;1])+[1;1]);            
        end
        
        coordTIFF = coordDMD + repmat([ligne_R_TIFF(1);colonne_R_TIFF(1)], 1, size(coordDMD,2));


    end
end

