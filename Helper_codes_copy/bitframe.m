classdef bitframe < handle
    %SPIKES_TREATED Summary of this class goes here
    %   Detailed explanation goes here
    properties
        DoubleFrames
        BitFrames
        fid
    end
    
    methods
        function obj = bitframe(fid)   
            obj.fid = fid;
        end
        function [] = add(obj,frame)
                pasta = reshape(frame,[],1);
                bytes_nb = 1:8:length(pasta)-7;
                    bytes = pasta(bytes_nb)*2^7 +  pasta(bytes_nb+1)*2^6 + ...
                         pasta(bytes_nb+2)*2^5 +  pasta(bytes_nb+3)*2^4 + ...
                        pasta(bytes_nb+4)*2^3 +  pasta(bytes_nb+5)*2^2 + ...
                       pasta(bytes_nb+6)*2^1 +  pasta(bytes_nb+7)*2^0;
                    fwrite(obj.fid,uint8(bytes),'uint8');
                
        end
    end
    
end


