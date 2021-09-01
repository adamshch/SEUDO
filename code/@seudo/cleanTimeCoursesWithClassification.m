function tcClean = cleanTimeCoursesWithClassification(se,whichTimeCourses,varargin)
% tcClean = seudoObject.cleanTimeCoursesWithClassification(se,whichTimeCourses,...)
%
% use current classification to zero unwanted transients
%
% 
% EXAMPLES:
%
%     tcClean = se.cleanTimeCoursesWithClassification('default');
%     tcClean = se.cleanTimeCoursesWithClassification('default','keepUnclassified',0);
%
%
% OUTPUT: 
%
%     tcClean - T x C matrix, the cleaned time courses 
%
%
% PARAMETERS:
%   
%   keepTrue - if true (default), leave transients classified as true intact. otherwise, set to zero.
%   keepFalse - default = false
%   keepAmbiguous - default = false
%   keepUnclassified - default = true
%
%
% 2018 Jeff Gauthier



p = inputParser;
p.addParamValue('keepTrue',true);
p.addParamValue('keepFalse',false);
p.addParamValue('keepMixed',false);
p.addParamValue('keepUnclassified',true);
parse(p,varargin{:});



% get time course struct of
if nargin < 2, whichTimeCourses = 'default'; end
structSpec = pickTcStruct(se,whichTimeCourses);
S = se.(structSpec{1})(structSpec{2});




% start with standard time courses
tcClean = S.tc;

% go through each cell and classified transient
for cc = 1:size(S.tc,2)
    T = S.transientInfo(cc);
    for tt = 1:length(T.classification)
        
        % zero this transient if it should not be kept
        
        if T.classification(tt) == seudo.valTrue
            % true
            if ~p.Results.keepTrue
                tcClean(T.times(tt,1):T.times(tt,2),cc) = 0;
            end
            
        elseif T.classification(tt) == seudo.valMix
            % Mixed
            if ~p.Results.keepMixed
                tcClean(T.times(tt,1):T.times(tt,2),cc) = 0;
            end
            
        elseif T.classification(tt) == seudo.valFalse
            % false
            if ~p.Results.keepFalse
                tcClean(T.times(tt,1):T.times(tt,2),cc) = 0;
            end
            
        elseif isnan(T.classification(tt))
            % unclassified
            if ~ p.Results.keepUnclassified
                tcClean(T.times(tt,1):T.times(tt,2),cc) = 0;
            end
            
        else
            error('classification for cell %d, transient %d not recognized',cc,tt)
        end
        
        
    end
end



