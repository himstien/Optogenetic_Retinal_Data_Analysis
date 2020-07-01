% function createNeighborStimMatrix()

% numNeighbors = input('Enter number of pixels to create: ');
clear temp;

matrix = [];
finalMatrix = [];

for eachPixel = implant.pixel.number
    temp = zeros(1, 91);
    temp(1, eachPixel) = 1;
    matrix = [matrix; [temp my_bin2dec(num2str(temp))]];
end

for eachPixel = implant.pixel.number
    for numNeighbors = 0
        vectorCombinations = combnk(1:6, numNeighbors);
        for comb = 1:size(vectorCombinations, 1)
            for dir = 1:size(vectorCombinations, 2)
                temp = zeros(1, 91);
                temp(1, eachPixel) = 1;
                tempNum = implant.pixel.neighborsHex{eachPixel, vectorCombinations(comb, dir)};
                index = find(tempNum(1) == implant.pixel.coordinates(:,1) & tempNum(2) == implant.pixel.coordinates(:,2) & tempNum(3) == implant.pixel.coordinates(:,3));
                temp(1, index) = 1;
                matrix = [matrix; [temp my_bin2dec(num2str(temp))]];                
            end
        end
    end
end
[~, ia] = unique(matrix(:, 92));
finalMatrix = [finalMatrix; matrix(ia, 1:end-1)]; 
