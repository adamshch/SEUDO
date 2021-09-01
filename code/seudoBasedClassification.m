function classOut = seudoBasedClassification(se, varargin)

if nargin > 1
    threshVal = varargin{1};
else
    threshVal = 0.5;
end

classOut = cell(se.nCells,1);                                              % Initialize the class output

for cc = 1:se.nCells
    
    classOut{cc}.profLevels =  se.tcSeudo.perTran.tc{cc};
    classOut{cc}.blobLevels = sum(se.tcSeudo.perTran.extras.tcBlobs{cc},2);
    classOut{cc}.ratios     = classOut{cc}.profLevels./classOut{cc}.blobLevels;
    if (~isempty(threshVal))&&(~isnan(threshVal))
        classOut{cc}.autoClass = classOut{cc}.ratios>threshVal;
    end
end

end



