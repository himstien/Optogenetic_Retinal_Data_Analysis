%% Permet � l'utilisateur de calibrer le rep�re du DMD avec 3 points cliqu�s
function [ligne_R_TIFF,colonne_R_TIFF] = calibrationRepere(pictureNameCalibration)
%X_R_DMD : abscisses des points du rep�re du DMD : Origine, ligne (bas
%gauche), colonne (haut droite)
%Y_R_DMD : ordon�es des points du rep�re du DMD : Origine, ligne (bas
%gauche), colonne (haut droite)

I=imread(pictureNameCalibration);
imshow(I); % permet � l'utilisateur de visualiser la photo TIF
for i=1:3
    switch i
        case 1
            title('\bf Click the haut-gauche/top-left mark on the DMD.');
            
[colonne_R_TIFF, ligne_R_TIFF]=ginput(3); %r�cup�re les trois points cliqu�s
close figure 1;
end