function [meanKeptFractionTrue,meanKeptFractionFalse,keptFractionsTrueByCell,keptFractionsFalseByCell] = evaluateNewTimeCourses(tcDefault,tcNew)
% evaluateNewTimeCourses(tcDefault,tcNew)
%
% INPUTS
%
%   tcDefault - struct specifying original time courses, needs fields:
%                   transientInfo.classification
%                   tc
%
%   tcNew - F x C matrix of new time courses
%
%


keptFractionsByCell = cell(size(tcDefault.tc,2),1);

keptFractionsTrueByCell = cell(size(tcDefault.tc,2),1);
keptFractionsFalseByCell = cell(size(tcDefault.tc,2),1);

keptFractionsTrue = [];
keptFractionFalse = [];

for cc = 1:size(tcDefault.tc,2)
    
    % get info on each transient
    
    %ti = tcDefault.transientInfo(cc);
    keptFractionsByCell{cc} = nan(size(tcDefault.transientInfo(cc).times,1),1);
    for tt = 1:size(tcDefault.transientInfo(cc).times,1)
        % get original and new time courses
        origTrace = tcDefault.tc(tcDefault.transientInfo(cc).times(tt,1):tcDefault.transientInfo(cc).times(tt,2),cc);
        newTrace = tcNew(tcDefault.transientInfo(cc).times(tt,1):tcDefault.transientInfo(cc).times(tt,2),cc);
        
        % compute kept fraction
        keptFractionsByCell{cc}(tt) = sum(newTrace)/sum(origTrace);
    end
    
    
    % analyze transients, grouped by transient type
    
    keptFractionsTrue = [keptFractionsTrue; keptFractionsByCell{cc}(tcDefault.transientInfo(cc).classification == seudo.valTrue)];
    keptFractionFalse = [keptFractionFalse; keptFractionsByCell{cc}(tcDefault.transientInfo(cc).classification == seudo.valFalse)];
    
    keptFractionsTrueByCell{cc} = keptFractionsByCell{cc}(tcDefault.transientInfo(cc).classification == seudo.valTrue);
    keptFractionsFalseByCell{cc} = keptFractionsByCell{cc}(tcDefault.transientInfo(cc).classification == seudo.valFalse);
end

% get overall mean
%meanKeptFractionTrue = mean(keptFractionsTrue);
%meanKeptFractionFalse = mean(keptFractionFalse);

meanKeptFractionTrue = keptFractionsTrue;
meanKeptFractionFalse = keptFractionFalse;


