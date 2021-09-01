function [x,y] = getFalconPlotPoints(se,computeNewValues,userStructChoice,varargin)
% compute values for FaLCon (Frequency and Level of Contamination) plots
% 
%
% USAGE
%
%   [x,y] = seudoObject.getFalconPlotPoints(computeNewValues,whichTcStruct,...)
%
%
% INPUTS
%
%   computeNewValues - if true,  use stored values (e.g. tcDefault.transientInfo().autoClass.cfrac)
%                      if false, compute new values using autoClassifyTransients
%
%   userStructChoice - which time course struct to operate on, e.g. 'default'
%                      see documentation for pickTcStruct in seudo.m
%   
%
% PARAMS
%
%   same parameters as for autoClassifyTransients
%
%
% EXAMPLES
%
%   [x,y] = seudoObject.getFalconPlotPoints;
%   [x,y] = seudoObject.getFalconPlotPoints(true); 
%   [x,y] = seudoObject.getFalconPlotPoints([],'default');
%   [x,y] = seudoObject.getFalconPlotPoints([],[],'alpha',0.03,'radius',[2 7]);  
%
%
% 2018 Jeff Gauthier





% PARSE & VALIDATE INPUTS

% select time course struct
if ~exist('userStructChoice','var') || isempty(userStructChoice), userStructChoice = 'default'; end
structSpec = pickTcStruct(se,userStructChoice);


% identify whether all fields exist and have enough entries
hasFields = false(se.nCells,1);
for cc = 1:se.nCells
    if isfield(se.(structSpec{1})(structSpec{2}).transientInfo(cc),'autoClass') && ...
            isfield(se.(structSpec{1})(structSpec{2}).transientInfo(cc).autoClass,'cfrac') && ...
            length(se.(structSpec{1})(structSpec{2}).transientInfo(cc).autoClass.cfrac) == 1 && ...
            isfield(se.(structSpec{1})(structSpec{2}).transientInfo(cc).autoClass,'rfrac') && ...
            length(se.(structSpec{1})(structSpec{2}).transientInfo(cc).autoClass.rfrac) == 1
        hasFields(cc) = true;
    else
        break
    end
end
hasFields = all(hasFields);


% decide whether to compute new values
if exist('computeNewValues','var') && ~isempty(computeNewValues)
    % check for error
    if ~computeNewValues && ~all(hasFields)
        error(sprintf(['Auto classification does not exist. Use autoClassify, or try:\n'...
            '   seudoObject.getFaLConPlotPoints(true,...)']))
    end
else
    % default to using existing parameters, unless they're not there or the user specified classification parameters
    if hasFields && nargin < 3
        computeNewValues = false;
    else
        computeNewValues = true;
    end
end




% GET POINTS

if computeNewValues
    % compute statistics from scratch
    
    fprintf('computing parameters for FaLCon plot points...\n')
    
    % automatically classify transients
    autoClass = se.autoClassifyTransients('default','saveResults',false,'plotComparison',false,varargin{:});
    
    switch autoClass(1).params.method
        case 'lbq'
            x = [autoClass.cfrac];
            y = [autoClass.rfrac];
        case 'corr'
            x = [autoClass.cfrac];
            y = [autoClass.rfrac];
        otherwise
            error('auto classification method ''%s'' not recognized, likely a bug in the code, so please report it!',autoClass.method)
    end
    
    
else
    % use existing auto classification
    
    fprintf('using existing parameters for FaLCon plot points\n')
    
    x = nan(se.nCells,1);
    y = nan(se.nCells,1);
    for cc = 1:se.nCells
        x(cc) = se.(structSpec{1})(structSpec{2}).transientInfo(cc).autoClass.cfrac;
        y(cc) = se.(structSpec{1})(structSpec{2}).transientInfo(cc).autoClass.rfrac;
    end
end




