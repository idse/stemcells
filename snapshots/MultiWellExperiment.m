classdef MultiWellExperiment < handle
    
    properties
        
        metadata 
        
        conditionLabels;    % cell array of condition names
        conditionIndex;     % cell array of vectors of position indices
        positions;          % nLocation x nTime array
    end
    
    methods
       
        function this = multiwellExperiment(metadata, conditionLabels)

            this.metadata = metadata;
            this.conditionLabels = conditionLabels;
        end
        
        % data processing
        %---------------------
        
        function extractData(this)
            
            for pi = 1:numel(this.positions)
                
                this.positions(pi).extractData();
            end
        end
    end
end