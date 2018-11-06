classdef MetadataMicropatternAndor < MetadataAndor
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
        
        function this = MetadataMicropatternAndor(filename)

            this = this@MetadataAndor(filename);
            
            % default margin outside colony to process, in pixels
            this.colMargin = 10; 
        end
    end
end