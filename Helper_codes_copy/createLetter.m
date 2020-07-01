
function image = createLetter(letter, width)


switch letter
    case 'E'
%         disp('EE');
        image = zeros(width*5, width*3);
        image(:, 1:width) = 255; %Vertical bar
        image(1:width, :) = 255; % Up Hor bar
        image(width*2:width*2+width, :) = 255; % Middle Hor bar
        image(width*4:width*4+width, :) = 255; % Bott Hor bar    
    case 'X'
%         disp('XX');
        image1 = zeros(width*5, width*3);
        image1(:, size(image1,2)/2-width/2:size(image1,2)/2+width/2) = 255; % stadard vertical bar at the middle of the image
        image = imrotate(image1, 30, 'crop');  % rotate the bar by 30 deg
        image = image + imrotate(image1, -30, 'crop'); % % rotate the bar by -30 deg and add it to the previous bar
        image(image>=255) = 255;
    case 'T'
        image = zeros(width*5, width*3);
        image(:, size(image,2)/2-width/2:size(image,2)/2+width/2) = 255; % stadard vertical bar at the middle of the image  
        image(1:width, :) = 255; % Top Hor bar      
    case 'P'
        [x, y] = meshgrid(1:3*width/2, 1:3*width);
        image = zeros(width*5, width*3);
        image(:, 1:width) = 255; %Vertical bar
        image(1:width, 1:3*width/2) = 255; % Up Hor bar
        image(width*2:width*2+width, 1:3*width/2) = 255; % Middle Hor bar
        temp = 255* (hypot(x , y - 3*width/2) < 3*width/2 & hypot(x , y - 3*width/2) > width/2);
        image(1:3*width, 2*width-width/2:3*width-1) = temp;
    otherwise
        image = [];
        return;
end
        
        