function [xStair,yStair] = lineToStairCase(x,y)
% make points for staircase (like the top edge of a bar plot)
%
% [xStair,yStair] = lineToStairCase(x,y)
%






% ensure vectors
assert(numel(x) == length(x))
assert(numel(y) == length(y))

% ensure same length
assert(length(x) == length(y))

if isempty(x)
    xStair = [];
    yStair = [];
    return
end

% make column
x = x(:);
y = y(:);

if length(x) == 1
    xStair = x;
    yStair = y;
    return
end


leftOverhang = (x(2)-x(1))/2;
rightOverhang = (x(end)-x(end-1))/2;


xMid = nanmean([x((1:end-1)') x((2:end)')],2);

xStair = [x(1)-leftOverhang; reshape([xMid xMid]',[],1); x(end)+rightOverhang];
yStair = reshape([y(:)'; y(:)'],[],1);


% add points to base at 0
xStair = xStair([1 1:end end]);
yStair = [0; yStair; 0];

