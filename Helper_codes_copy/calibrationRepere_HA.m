%% Permet � l'utilisateur de calibrer le rep�re du DMD avec 3 points cliqu�s
function [X, Y, I] = calibrationRepere_HA(pictureNameCalibration)
%X_R_DMD : abscisses des points du rep�re du DMD : Origine, ligne (bas
%gauche), colonne (haut droite)
%Y_R_DMD : ordon�es des points du rep�re du DMD : Origine, ligne (bas
%gauche), colonne (haut droite)

X = [];
Y = [];

I=imread(pictureNameCalibration);


I = squeeze(I(:,:,2)); % ROmain's setup
I = 2^16 - I; % ROmain's setup

I = I'; % for Romain's DMD

imshow(I); % permet � l'utilisateur de visualiser la photo TIF
hold on;


for i=1:3
    switch i
        case 1
            title('\bf Click the haut-gauche/top-left mark on the DMD.');
            [x_TL, y_TL, but] = ginput(1);
            drawCircle([x_TL, y_TL], 3, 'r');
            if(~isequal(but, 1))
                return;
            end
        case 2
            title('\bf Click the basse-gauche/bottom-left mark on the DMD.');
            [x_BL, y_BL, but] = ginput(1); 
            drawCircle([x_BL, y_BL], 3, 'r');
            if(~isequal(but, 1))
                return;
            end    
        case 3
            title('\bf Click the haut-droit/bottom-right mark on the DMD.');
            [x_TR, y_TR, but] = ginput(1);
            drawCircle([x_TR, y_TR], 3, 'r');
            if(~isequal(but, 1))
                return;
            end
    end            
end            
            
X = [x_TL x_BL x_TR];
Y = [y_TL y_BL y_TR];

close figure 1;
end