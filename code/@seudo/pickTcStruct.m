function structSpec = pickTcStruct(se,userChoice)
% parse user selection of a particular struct of time courses to do analysis on
%
%
%       specify which time courses to inspect with one argument:
%           'default' - time courses provided by user, or LSQ
%           'contam' - high contamination time courses
%           'seudo' - seudo time courses computed most recently, i.e. the last set
%           {'seudo',#} - seudo time courses, # is which set of seudo results
%               if # is negative, it counts from the end
%
%       examples:
%           se.classifyTransients('default')
%           se.classifyTransients('contam')
%           se.classifyTransients('seudo')
%           se.classifyTransients({'seudo',3})
%           se.classifyTransients({'seudo',-1})
%
%
% OUTPUT:    { tcFieldName , whichElementInStruct }
%
%


% get user input for field name
if ischar(userChoice)
    whichTC = userChoice;
elseif iscell(userChoice) && length(userChoice) == 2
    whichTC = userChoice{1};
else
    disp(userChoice)
    error('time course struct specification unrecognized')
end


% get actual field name
switch lower(whichTC)
    case {'default','tcdefault'}
        structName = 'tcDefault';
    case {'contam','tccontam'}
        structName = 'tcContam';
    case {'seudo','tcseudo'}
        structName = 'tcSeudo';
    otherwise
        error('time course struct specification unrecognized: %s',userChoice)
end


% default to last one
if iscell(userChoice)
    structNum = userChoice{2};
else
    structNum = length(se.(structName));
end




% ensure number is in range
if structNum > length(se.(structName)) || structNum + length(se.(structName)) < 1
    error('specified set %d, but field %s only has %d sets',structNum,structName,length(se.(structName)))
end

% if negative, take from end
if structNum < 0
    structNum = length(se.(structName)) + structNum; 
end


% combine in to a cell array
structSpec = {structName,structNum};



