function result = ring_hex(center, scalingFactor, radius, directions)
    result = [];    
    
    temp =  center + directions{5}*radius*scalingFactor;    
    
    for dir = 1:6        
        for r = 1:radius
            result = [result; temp];
            temp = temp + directions{dir}*scalingFactor;        
        end
    end      
%     result = result;
end

