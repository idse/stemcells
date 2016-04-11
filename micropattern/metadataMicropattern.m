classdef metadataMicropattern < metadata
    % metadata with additional properties for Andor IQ output
    
    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------

    properties
        
        colRadiiMicron
        colRadiiPixel
        colMargin
    end
    
    methods
        
        function this = metadataMicropattern(filename)

            this = this@metadata(filename);
        end
    end
end