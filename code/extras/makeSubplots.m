function [h,hMat] = makeSubplots(figNum,nX,nY,gapX,gapY,bounds)
% [h,hMat] = makeSubplots(figNum,nX,nY,gapX,gapY,bounds)
%
% figNum        which figure to use
% nX, nY        number of subplots
% gapX, gapY    gap between plots, expressed as fraction of plot width
% bounds        4-element vector of bounds within the figure.  [0 0 1 1] is whole figure
%
% h             vector of all subplots
% hMat          matrix of subplots, arranged like the subplots themselves
%


% check inputs
%if ~isPositiveIntegerValuedNumeric([nX nY])
if ~(nX>0&&nY>0&&round(nY)==nY&&round(nX)==nX)
    error('nX (%0.1f) and nY (%0.1f) must be positive integers',nX,nY)
end
if gapX<0 || gapY<0
    error('gapX (%0.1f) and gapY (%0.1f) must be positive',gapX,gapY)
end



% get figure number
if figNum==0;figure;figNum=gcf;end

% be sure figure exists
if ~ishandle(figNum)
    figure(figNum)
end


% identify subplot size
spacingX = bounds(3)/(nX+gapX/(1+gapX));
spacingY = bounds(4)/(nY+gapY/(1+gapY));
sizeX = spacingX/(1+gapX);
sizeY = spacingY/(1+gapY);


% initialize hMat
hMat = zeros(nY,nX);

% start counter
ss = 1;
% cycle through all locations, beginning in the upper left and moving right
for yy = 1:nY
    for xx = 1:nX
        % get handle to this subplot
        %hS = subplot('position',[bounds(1)+spacingX*(xx-1) bounds(2)+bounds(4)-spacingY*yy sizeX sizeY]);
        hS = axes('parent',figNum,'position',[bounds(1)+spacingX*(gapX + xx-1) bounds(2)+bounds(4)-spacingY*yy sizeX sizeY]);
        % store in the matrix
        hMat(yy,xx) = hS;
        % increment counter
        ss = ss+1;
    end
end

% reshape hMat as a vector, being sure to start in the upper left and moving right
h = reshape(hMat',[],1);
