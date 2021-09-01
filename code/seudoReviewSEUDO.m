function varargout = seudoReviewSEUDO(varargin)
% SEUDOREVIEWSEUDO MATLAB code for seudoReviewSEUDO.fig
%      SEUDOREVIEWSEUDO, by itself, creates a new SEUDOREVIEWSEUDO or raises the existing
%      singleton*.
%
%      H = SEUDOREVIEWSEUDO returns the handle to a new SEUDOREVIEWSEUDO or the handle to
%      the existing singleton*.
%
%      SEUDOREVIEWSEUDO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEUDOREVIEWSEUDO.M with the given input arguments.
%
%      SEUDOREVIEWSEUDO('Property','Value',...) creates a new SEUDOREVIEWSEUDO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before seudoReviewSEUDO_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to seudoReviewSEUDO_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help seudoReviewSEUDO

% Last Modified by GUIDE v2.5 19-Oct-2018 01:03:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @seudoReviewSEUDO_OpeningFcn, ...
                   'gui_OutputFcn',  @seudoReviewSEUDO_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before seudoReviewSEUDO is made visible.
function seudoReviewSEUDO_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to seudoReviewSEUDO (see VARARGIN)

% Choose default command line output for seudoReviewSEUDO
handles.output = hObject;








% store inputs in handles struct
handles.se = varargin{1};
handles.tcStructSpecDefault = varargin{2};
handles.tcStructSpecSeudo = varargin{3};


% note # cells
handles.nCells = handles.se.nCells;


% hide axesTransientsArea, since it is only used to set bounds for subplots
set(handles.axesPlotAreaTransients,'visible','off')

% initialize set of transient subplots to be empty
%handles.axesTransSubplots = [];



% get maximum of all profiles as an image
profMax = max(bsxfun(@rdivide,handles.se.profiles,max(max(handles.se.profiles))),[],3);
profMax = repmat(profMax / max(profMax(:)),[1 1 3]);
handles.cellProfilesMax = profMax;

% get COMs
handles.profileCOMs = seudo.computeRoiCOMs(handles.se.profiles);



% set slider range for cell
set(handles.sliderCellNumber,'Max',handles.nCells)
set(handles.sliderCellNumber,'Min',1)
if handles.nCells > 1
    set(handles.sliderCellNumber,'SliderStep',[1 10]/(handles.nCells-1))
else
    
end


% SELECT INITIAL VALUES

% cell number
handles.thisCell = 1;

% transient count
handles.nTransients = 3;

% frame count
handles.nFrames = 10;

% seudo number
handles.thisSeudoNumber = handles.tcStructSpecSeudo{2};
set(handles.popupmenuSeudoNumber,'value',handles.thisSeudoNumber)
%set(handles.popupmenuSeudoNumber,'string',num2cellstr(1:length(handles.se.tcSeudo)))
set(handles.popupmenuSeudoNumber,'string',cellstr(num2str((1:length(handles.se.tcSeudo))')))

%  cell selection
handles.thisCell = 1;
set(handles.sliderCellNumber,'Value',1)
set(handles.editTextCellNumber,'String','1')

% transient selection


% frames selection



axis(handles.axesPatch,'off')
patch([0 1 1 0 0],[0 0 1 1 0],'k','facecolor',[1 1 1]*.5,'parent',handles.axesPatch)



% call functions to update all plots
handles = updateOverviewPlots(hObject,handles);

% scatter contamination
handles = computeScatterContam(hObject,handles);

handles = makeNewTransientSubplots(hObject,handles);
handles = setNewTransientCount(hObject,handles);
handles = setNewSeudoNumber(hObject,handles);






% set colors of text labels
set(handles.text_cell1,'foregroundcolor',seudo.pickColor('SEcell'))
set(handles.text_blobs1,'foregroundcolor',seudo.pickColor('SEblobs'))
set(handles.text_LSQ,'foregroundcolor',seudo.pickColor('LSQ'))
set(handles.text_cell2,'foregroundcolor',seudo.pickColor('SEcell'))
set(handles.text_blobs2,'foregroundcolor',seudo.pickColor('SEblobs'))
set(handles.text_true,'foregroundcolor',seudo.pickColor('true'))
set(handles.text_false,'foregroundcolor',seudo.pickColor('false'))

set(handles.text_LSQ2,'foregroundcolor',seudo.pickColor('LSQ'))
set(handles.text_cell3,'foregroundcolor',seudo.pickColor('SEcell'))
set(handles.text_blobs3,'foregroundcolor',seudo.pickColor('SEblobs'))




colormap(handles.axesPlotAreaMovie,'gray')
    
    
% Update handles structure
guidata(hObject, handles);





% --- Outputs from this function are returned to the command line.
function varargout = seudoReviewSEUDO_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;














% MAKING A NEW SELECTION








% NEW TRANSIENT COUNT

function handles = setNewTransientCount(figHandle,handles)

% make plots
handles = makeNewTransientSubplots(figHandle,handles);

% act as if new cell was chosen to update transient slider and transient numbers to plot
handles = setNewCell(figHandle,handles);





% NEW SET OF SEUDO RESULTS

function handles = setNewSeudoNumber(hObject,handles)

handles.thisSeudoNumber = get(handles.popupmenuSeudoNumber,'value');

if DEBUG_TEXT, fprintf('SEUDO number: %d\n',handles.thisSeudoNumber), end


% set title
set(hObject,'name',sprintf('%s: SEUDO results (%d), [%d x %d] pixels, %d frames, %d cells',...
    handles.se.name,handles.tcStructSpecSeudo{2},handles.se.movY,handles.se.movX,handles.se.movF,handles.se.nCells))

% delete current frames
handles.transientMovie = struct;

% update all plots
handles = updateOverviewPlots(hObject,handles);
handles = updateCellPlots(hObject,handles);
handles = updateTransientPlots(hObject,handles);
handles = updateFramePlots(hObject,handles);

guidata(hObject,handles)




% NEW CELL

function handles = setNewCell(hObject,handles)

if DEBUG_TEXT, fprintf('selected cell %d\n',handles.thisCell), end

% update cell plots (this sets sort order for transients)
handles = updateCellPlots(hObject,handles);


% note number of transients in this cell
transCount = size(handles.se.(handles.tcStructSpecDefault{1})(handles.tcStructSpecDefault{2}).transientInfo(handles.thisCell).times,1);
set(handles.textTransientCount,'string',sprintf('%d transients',transCount))


% set up slider

% select new transients
whichTransients = nan(1,handles.nTransients);
nonNanTrans = min(handles.nTransients,transCount);
whichTransients(1:nonNanTrans) = handles.transientSortOrder(1:nonNanTrans);

handles.whichTransients = whichTransients;


% note new selection
set(handles.editTransientsToPlot,'string',makeCommaDelimitedString(handles.whichTransients))

% set slider
set(handles.sliderTransients,'enable','on')
sliderValueCount = transCount - handles.nTransients + 1;
if sliderValueCount > 1
    set(handles.sliderTransients,'Max',sliderValueCount)
    set(handles.sliderTransients,'Min',1)
    set(handles.sliderTransients,'SliderStep',[1 handles.nTransients]/(sliderValueCount-1))  
    set(handles.sliderTransients,'Value',1)
else
    set(handles.sliderTransients,'enable','off')
end


% initialize transients to display
handles = setNewTransients(hObject,handles);







% NEW SET OF TRANSIENTS

function handles = setNewTransients(hObject,handles)

% plot the transients
handles = updateTransientPlots(hObject,handles);


% choose set of movie frames

handles.transientForFrames = handles.whichTransients(1);


% set slider
set(handles.sliderTransientForFrames,'enable','on')
transCount = length(handles.whichTransients);
if transCount > 1
    set(handles.sliderTransientForFrames,'Max',transCount)
    set(handles.sliderTransientForFrames,'Min',1)
    set(handles.sliderTransientForFrames,'SliderStep',[1 10]/(transCount-1))  
    set(handles.sliderTransientForFrames,'Value',1)
else
    set(handles.sliderTransientForFrames,'enable','off')
end


% initialize movie frames to display
handles = setNewTransientForFrames(hObject,handles);





% NEW TRANSIENT FOR MOVIE FRAMES

function handles = setNewTransientForFrames(hObject,handles)

% update display of which transient is being shown
set(handles.textTransientForFrames,'string',sprintf('frames for\ntransient %d',handles.transientForFrames))

% get frame count
timeWindow = handles.se.(handles.tcStructSpecDefault{1})(handles.tcStructSpecDefault{2}).transientInfo(handles.thisCell).times(handles.transientForFrames,:);
frameCount = diff(timeWindow)+1;


% update frame count
set(handles.textFrameCount,'string',sprintf('%d frames',frameCount))


% set slider
sliderCount = frameCount - handles.nFrames + 1;
set(handles.sliderMovieFrames,'enable','on')
if sliderCount > 1
    set(handles.sliderMovieFrames,'Max',sliderCount)
    set(handles.sliderMovieFrames,'Min',1)
    set(handles.sliderMovieFrames,'SliderStep',[1 handles.nFrames]/(sliderCount-1))  
    set(handles.sliderMovieFrames,'Value',1)
else
    set(handles.sliderMovieFrames,'enable','off')
end

% select frames
handles.whichFrames = 1:handles.nFrames;

% update plot
handles = updateFramePlots(hObject,handles);



% NEW NUMBER OF TRANSIENTS


function handles = makeNewTransientSubplots(figHandle,handles)


if DEBUG_TEXT, fprintf('create %d transient plots\n',handles.nTransients), end


% make subplots
%posFig = get(figHandle,'position'); % use if fixed size
posFig = [1 1 1 1]; % use if figure is resizable


% make subplots
nX = handles.nTransients;
nY = 5;

% create one big subplot, just to be sure previous transient plots
subplot('position',get(handles.axesPlotAreaTransients,'position')./posFig([3 4 3 4]),...
    'parent',figHandle)


[~,handles.axesTransSubplots] = makeSubplots(get(handles.axesPlotAreaTransients,'parent'),...
    nX,nY,.1,.2,get(handles.axesPlotAreaTransients,'position')./posFig([3 4 3 4]));


% for h = handles.axesTransSubplots(:)'
%     axis(h,'image')
%     axis(h,'ij')
%     axis(h,'off')
%     hold(h,'on')
%     imagesc(rand(10),'parent',h)
% end















% PLOTTING











% UPDATE OVERVIEW PLOTS

function handles = updateOverviewPlots(hObject,handles)


if DEBUG_TEXT, fprintf('updating overview plots for SEUDO number: %d\n',handles.thisSeudoNumber), end






% get kept fraction for all transients



% pull out relevant variables
se = handles.se;
tcDefault =  handles.se.(handles.tcStructSpecDefault{1})(handles.tcStructSpecDefault{2}); %.transientInfo(handles.thisCell).times
tcStructSeudo = handles.se.(handles.tcStructSpecSeudo{1})(handles.thisSeudoNumber);

allFracTrue = nan(se.nCells,1);
allFracFalse = nan(se.nCells,1);

for cc = 1:se.nCells
    
    extras = tcStructSeudo.extras{cc};
    params = tcStructSeudo.params{cc};
    tcLSQ = tcStructSeudo.tcLSQ(:,cc);
    tcSEUDO = tcStructSeudo.tc(:,cc);
    
    ti = tcDefault.transientInfo(cc);
    
    nTrans = size(ti.times,1);
    
    
    % transient summary
    
    fractionCell = nan(nTrans,1);
    transArea = nan(nTrans,1);
    
    % for each transient, get total fraction kept
    for tt = 1:nTrans
        inTrans = ti.times(tt,1):ti.times(tt,2);
        fractionCell(tt) = evaluateOneTimeCourse(tcLSQ(inTrans),tcSEUDO(inTrans),'method','area');
        transArea(tt) = sum(tcLSQ(inTrans));
    end
    
    % keep in [0 1]
    fractionCell(fractionCell<0) = 0;
    fractionCell(fractionCell>1) = 1;
    
    % get overall fraction
    isTrue = ti.classification == seudo.valTrue;
    isFalse = ti.classification == seudo.valFalse;
    switch 2
        case 1 % mean of fractions
            fracTrue = mean(fractionCell(isTrue));
            fracFalse = mean(fractionCell(isFalse));
        case 2 % total fractional area
            fracTrue = sum(fractionCell(isTrue).*transArea(isTrue))./sum(transArea(isTrue));
            fracFalse = sum(fractionCell(isFalse).*transArea(isFalse))./sum(transArea(isFalse));
    end
    
    % if no true or false transients, give negative value
    if sum(isTrue) == 0, fracTrue = -.1; end
    if sum(isFalse) == 0, fracFalse = -.1; end
    
    allFracTrue(cc) = fracTrue;
    allFracFalse(cc) = fracFalse;
    
    
end



% plot
h = handles.axesScatterEval;
cla(h)

% lines
set(h,'box','off','tickdir','out','xtick',0:.2:1,'ytick',0:.2:1)
grid(h,'on')
hold(h,'on')
plot(h,[0 1],[0 1],'-','color',[1 1 1]*.5)
axis(h,'image')
plot(h,[0 0 1 1 0],[0 1 1 0 0],'k-')

% points
plot(h,allFracFalse,allFracTrue,'.','color',[1 1 1]*.2)


% labels
switch 2
    case 1 % mean of fractions
        xlabel(h,'false transients kept')
        ylabel(h,'true transients kept')
    case 2 % total fractional area
        xlabel(h,'false transient area kept')
        ylabel(h,'true transient area kept')
end


xlim(h,[-.2 1.1])
ylim(h,[-.2 1.1])


% add numbers
if isfield(handles,'axesScatterEvalNumbers')
    delete(handles.axesScatterEvalNumbers(isgraphics(handles.axesScatterEvalNumbers)))
end

handles.axesScatterEvalNumbers = nan(se.nCells,1);
for cc = 1:se.nCells
    handles.axesScatterEvalNumbers(cc) = text(allFracFalse(cc),allFracTrue(cc),[' ' num2str(cc)],'color',[1 1 1]*.2,'fontsize',8,'parent',h);
end

if get(handles.checkboxShowCellNumbers,'value')
    set(handles.axesScatterEvalNumbers,'visible','on')
else
    set(handles.axesScatterEvalNumbers,'visible','off')
end
    


% store result
handles.allFracTrue = allFracTrue;
handles.allFracFalse = allFracFalse;
guidata(hObject,handles)




% UPDATE CELL PLOTS

function handles = updateCellPlots(hObject,handles)


if DEBUG_TEXT, fprintf('updating cell plots for cell: %d\n',handles.thisCell), end


% pull out relevant variables
se = handles.se;
tcDefault =  se.(handles.tcStructSpecDefault{1})(handles.tcStructSpecDefault{2}); %.transientInfo(handles.thisCell).times
tcStructSeudo = se.(handles.tcStructSpecSeudo{1})(handles.thisSeudoNumber);
extras = tcStructSeudo.extras{handles.thisCell};
params = tcStructSeudo.params{handles.thisCell};
tcLSQ = tcStructSeudo.tcLSQ(:,handles.thisCell);
tcSEUDO = tcStructSeudo.tc(:,handles.thisCell);

ti = tcDefault.transientInfo(handles.thisCell);

nTrans = size(ti.times,1);



% get sum of blobs in the ROI
prof = se.profiles(:,:,handles.thisCell) > 0;
tcBLOBS_ROI = sum(extras.tcBlobs(:,prof(extras.whichPixels)),2)/sum(extras.blobFits(:,1));





% identify window for this cell

% get pixels used in SEUDO analysis
whichPixIm = reshape(extras.whichPixels,se.movY,se.movX);
winY = sum(whichPixIm,2)>0;
winX = sum(whichPixIm,1)>0;

% store window coordinates
handles.thisCellWindowGlobal = [find(winY,1,'first') find(winY,1,'last') find(winX,1,'first') find(winX,1,'last')];

% get pixels used for transient info shape analysis
win = ti.window;
whichPixInTIwindow = whichPixIm(win(1):win(2),win(3):win(4));

% note which part of this region was analyzed in SEUDO
winY = sum(whichPixInTIwindow,2)>0;
winX = sum(whichPixInTIwindow,1)>0;
handles.thisCellWindowInTIwindow = [find(winY,1,'first') find(winY,1,'last') find(winX,1,'first') find(winX,1,'last')];





% highlight profile in image of all profiles
thisProfile = handles.se.profiles(:,:,handles.thisCell);
thisProfile = thisProfile / max(thisProfile(:));
profMax = handles.cellProfilesMax;
col = CELL_HIGHLIGHT_COLOR;
thisProfile = cat(3,thisProfile*col(1),thisProfile*col(2),thisProfile*col(3));
imToPlot = 1-(1-profMax*.7).*(1-thisProfile);
%imToPlot = imToPlot .^ .5;

h = handles.axesAllProfiles;
cla(h)
imagesc(imToPlot,'parent',h)
axis(h,'image')
axis(h,'off')
hold(h,'on')
cutWin = handles.thisCellWindowGlobal + [-.5 .5 -.5 .5];
plot(h,cutWin([3 4 4 3 3]),cutWin([1 1 2 2 1]),'-','linewidth',1,'color',[0 1 1])
title(h,'all profiles','color','k','fontsize',12,'fontweight','normal','fontname','helvetica')

% add cell numbers

% add cell number labels
x = handles.profileCOMs(:,1);
y = handles.profileCOMs(:,2);
handles.axesProfileCOMs = nan(handles.se.nCells,1);
for cc = 1:handles.se.nCells
    handles.axesProfileCOMs(cc) = text(x(cc),y(cc),[' ' num2str(cc)],'color',[1 1 1]*1,'fontsize',8,'parent',h);
end

if get(handles.checkboxShowCellNumbers,'value')
    set(handles.axesProfileCOMs,'visible','on')
else
    set(handles.axesProfileCOMs,'visible','off')
end






% get image of this profile
win = handles.thisCellWindowGlobal;
thisProfile = se.profiles(win(1):win(2),win(3):win(4),handles.thisCell);


% get contour line for this cell
fym = contourc(thisProfile,[1 1]*max(thisProfile(:))/10);
C = {};
currInd = 1;
while currInd < size(fym,2)
    nPts = fym(2,currInd);
    lineInds = currInd+(1:nPts);
    C{end+1} = fym(:,lineInds);
    currInd = currInd + 1 + nPts;
end
handles.thisCellContours = C;


% plot cell shape
h = handles.axesProfile;
cla(h)
imagesc(thisProfile,'parent',h)
axis(h,'image')
axis(h,'off')
colormap(h,'gray')


% add contour lines
hold(h,'on')
for cc = 1:length(C)
    plot(h,C{cc}(1,:),C{cc}(2,:),'-','color',CELL_HIGHLIGHT_COLOR)
end




% transient summary

fractionCell = nan(nTrans,1);
fractionBlobs = nan(nTrans,1);
transArea = nan(nTrans,1);

% for each transient, get total fraction kept
for tt = 1:nTrans
    inTrans = ti.times(tt,1):ti.times(tt,2);
    fractionCell(tt) = evaluateOneTimeCourse(tcLSQ(inTrans),tcSEUDO(inTrans),'method','area');
    fractionBlobs(tt) = evaluateOneTimeCourse(tcLSQ(inTrans),tcBLOBS_ROI(inTrans),'method','area');
    transArea(tt) = sum(tcLSQ(inTrans));
end

% choose order for sorting transients
switch get(handles.popupmenuSortTransients,'value')
    
    case 1 % order by time
        sortOrder = 1:nTrans;
        
    case 2 % order by area
        [~,sortOrder] = sort(transArea,'descend');
        
    case 3 % order by classification
        isTrue = find(ti.classification == seudo.valTrue);
        isFalse = find(ti.classification == seudo.valFalse);
        isUnclassified = find(isnan(ti.classification));
        %isMixed = find(ti.classification == seudo.valMix);
        isMixed = find(~isnan(ti.classification) & ~ismember(ti.classification,[seudo.valTrue seudo.valFalse]));
        
        sortOrder = [isTrue; isFalse; isMixed; isUnclassified];
        
end


% plot fraction cells and blobs
h = handles.axesTransientFractions;
cla(h)
colormap(h,[seudo.pickColor('SEcell');seudo.pickColor('SEblobs')]);

b = bar([fractionCell(sortOrder) fractionBlobs(sortOrder)],'stacked','parent',h,'linestyle','none','facecolor','flat');
ylim(h,[-.1 1.1])
xlim(h,[0 nTrans+1])
box(h,'off')
set(h,'tickdir','out')
set(h,'xtick',[])
axis(h,'off')
%set(h,'visible','off'); set(get(h,'child'),'visible','on')

% manually set colormap depending on matlab version, because isn't matlab just the best?
if any(isprop(b,'CData'))
    for bb = 1:length(b), b(bb).CData = bb;end
end



% plot area
h = handles.axesTransientAreas;
cla(h)
colormap(h,[seudo.pickColor('true');seudo.pickColor('false');seudo.pickColor('mixed');[1 1 1]*.3]);


% add offset to actual area be sure transient color is visible
transArea = transArea + max(transArea)/10;

% put into format for stacked bar plot
T = repmat(transArea,1,4);
T(ti.classification ~= seudo.valTrue,1) = 0;
T(ti.classification ~= seudo.valFalse,2) = 0;
T(~isnan(ti.classification),4) = 0;
T(any(T(:,[1 2 4])~=0,2),3) = 0;

T = T(sortOrder,:);


b = bar(T,'stacked','parent',h,'linestyle','none','facecolor','flat');
ylim(h,[0 max(transArea)])
xlim(h,[0 nTrans+1])
box(h,'off')
set(h,'tickdir','out')
set(h,'xtick',[])
axis(h,'off')

% manually set colormap depending on matlab version, because isn't matlab just the best?
if any(isprop(b,'CData'))
    for bb = 1:length(b), b(bb).CData = bb;end
end


% store sort order
handles.transientSortOrder = sortOrder;

% highlight location in overview plots
if isfield(handles,'overviewHighlightHandles')% && isgraphics(handles.overviewHighlightHandles)
    delete(handles.overviewHighlightHandles(isgraphics(handles.overviewHighlightHandles)))
end

handles.overviewHighlightHandles = [];
handles.overviewHighlightHandles(1) = ...
    plot(handles.axesScatterEval,handles.allFracFalse(handles.thisCell),handles.allFracTrue(handles.thisCell),...
    'o','color',CELL_HIGHLIGHT_COLOR);

handles.overviewHighlightHandles(2) = ...
    plot(handles.axesScatterContam,handles.scatterContamX(handles.thisCell),handles.scatterContamY(handles.thisCell),...
    'o','color',CELL_HIGHLIGHT_COLOR);


guidata(hObject,handles)





% UPDATE TRANSIENT PLOTS

function handles = updateTransientPlots(hObject,handles)


if DEBUG_TEXT, fprintf('plotting transients: %s\n',makeCommaDelimitedString(handles.whichTransients)), end


% if first time for this cell, do lots of axes adjustments
% if not, avoid them to reduce plotting time
firstTime = 1;

tic


% clear axes
if firstTime
    for h = handles.axesTransSubplots(:)'
        cla(h)
    end
end


transToPlot = handles.whichTransients;
transToPlot = transToPlot(~isnan(transToPlot));

for tt = 1:length(transToPlot)
    
    thisTransient = transToPlot(tt);
    
    
    % pull out relevant variables
    se = handles.se;
    tcDefault =  handles.se.(handles.tcStructSpecDefault{1})(handles.tcStructSpecDefault{2}); %.transientInfo(handles.thisCell).times
    tcStructSeudo = handles.se.(handles.tcStructSpecSeudo{1})(handles.thisSeudoNumber);
    extras = tcStructSeudo.extras{handles.thisCell};
    params = tcStructSeudo.params{handles.thisCell};
    tcLSQ = tcStructSeudo.tcLSQ(:,handles.thisCell);
    tcSEUDO = tcStructSeudo.tc(:,handles.thisCell);
    
    % which frames
    framesToPlot = tcDefault.transientInfo(handles.thisCell).times(thisTransient,:);
    framesToPlot = framesToPlot(1):framesToPlot(2);
    
    
    % which spatial subset
    
    %win = tcDefault.transientInfo(handles.thisCell).window;
    
    %wY = diff(win(1:2))+1;
    %wX = diff(win(3:4))+1;

            
    
    win = handles.thisCellWindowGlobal;
    %handles.thisCellWindowInTIwindow;
    
    wY = diff(win([1 2]))+1;
    wX = diff(win([3 4]))+1;

    
%     whichPixIm = reshape(extras.whichPixels,se.movY,se.movX);
%     winY = sum(whichPixIm,2)>0;
%     winX = sum(whichPixIm,1)>0;
%     
%     wY = sum(winY)
%     wX = sum(winX)
    
    
    
    
    apply_blobs = @(z) reshape(convn(reshape(z, wY,wX), extras.one_blob, 'same'),wY*wX, 1);
    apply_blobs2D = @(z) convn(z, extras.one_blob, 'same');
    
    
    % get blob time course
    if params.saveBlobTimeCourse
        % if individual blob time courses were saved, combine across space to get the total
        
        % all blobs
        tcBlobsAll = nansum(extras.tcBlobs(framesToPlot,:),2)/sum(extras.blobFits(:,1));
        
        % blobs in the profile
        roi = reshape(se.profiles(win(1):win(2),win(3):win(4),handles.thisCell),[],1);
        tcBlobsInProfile = nansum(extras.tcBlobs(framesToPlot,roi>0),2)/sum(extras.blobFits(:,1));
        
        % image of the blob sum
        B = reshape(extras.tcBlobs(framesToPlot,:)',wY,wX,[]);
        imBlobs = max(convn(B, extras.one_blob, 'same'),[],3);
        
        % combine with cell
        imCell = max(tcSEUDO(framesToPlot)) * se.profiles(win(1):win(2),win(3):win(4),handles.thisCell);
        vals = [imBlobs(:); imCell(:)];
        colCell = seudo.pickColor('SEcell');
        colBlobs = seudo.pickColor('SEblobs');
        %imBlobs = cat(3,imBlobs,imCell,imBlobs) / max(vals);
        M = [];
        for cc = 1:3
            M(:,:,cc) = 1 - (1 - colCell(cc)*imCell) .* (1 - colBlobs(cc)*imBlobs);
        end
        imBlobs = M / max(vals(:));
        
    else
        % if individual blob time courses were not saved, combine across space to get the total
        tcBlobsAll = extras.tcBlobs(:,1)/sum(extras.blobFits(:,1));
        tcBlobsInProfile = extras.tcBlobs(framesToPlot,2)/sum(extras.blobFits(:,1));
        imBlobs = [];
    end
    
    

    % plot results
    

    % get color based on classification
    thisClass = tcDefault.transientInfo(handles.thisCell).classification(thisTransient);
    if ~isempty(thisClass)
        if      isnan(thisClass),            col = 'k';
        elseif  thisClass == seudo.valTrue,  col = seudo.pickColor('true');
        elseif  thisClass == seudo.valFalse, col = seudo.pickColor('false');
        else                                 col = seudo.pickColor('mixed');
        end
    else
        col = 'k';
    end
    
    
    % active shape
    h = handles.axesTransSubplots(1,tt);
    
    transProfile = nan(se.movY,se.movX);
    winTI = tcDefault.transientInfo(handles.thisCell).window;
    transProfile(winTI(1):winTI(2),winTI(3):winTI(4)) = tcDefault.transientInfo(handles.thisCell).shapes(:,:,thisTransient);
    transProfile = transProfile(win(1):win(2),win(3):win(4));
    
    
    if firstTime
        imagesc(transProfile,'parent',h)
        axis(h,'image')
        axis(h,'off')
        colormap(h,'gray')
        hold(h,'on')
        C = handles.thisCellContours;
        for cc = 1:length(C)
            plot(h,C{cc}(1,:),C{cc}(2,:),'-','color',CELL_HIGHLIGHT_COLOR)
        end
    else
        kids = get(h,'child');
        set(kids(end),'CData',transProfile)
    end
    title(h,sprintf('transient %d',thisTransient),'fontsize',12,'color',col)
    
    
    % fit of cell and blob max
    h = handles.axesTransSubplots(2,tt);
    %apply_blobs = @(z) reshape(convn(reshape(z, se.movY,se.movX), extras.one_blob, 'same'),se.movX*se.movY, 1);
    if firstTime
        imagesc(imBlobs,'parent',h)
        axis(h,'image')
        set(h,'xtick',[],'ytick',[])
        colormap(h,'gray')
        axis(h,'on')
    else
        set(get(h,'child'),'CData',imBlobs)
    end
    %axis(h,'off')
    
    
    
    % found cell time course
    h = handles.axesTransSubplots(3,tt);
    cla(h)
    plot(h,[1 length(framesToPlot)],[0 0],'-','color',[1 1 1]*.5)
    hold(h,'on')
    plot(h,tcLSQ(framesToPlot),'color',seudo.pickColor('LSQ'),'linewidth',3)
    plot(h,tcSEUDO(framesToPlot),'color',seudo.pickColor('SEcell'))
    plot(h,tcBlobsInProfile,'--','color',seudo.pickColor('SEblobs'))
    axis(h,'tight')
    set(h,'xtick',[])
    axis(h,'on')
    
    
    % blob TC
    h = handles.axesTransSubplots(4,tt);
    cla(h)
    plot(h,[1 length(framesToPlot)],[0 0],'-','color',[1 1 1]*.5)
    hold(h,'on')
    plot(h,tcBlobsAll,'-','color',seudo.pickColor('SEblobs'))
    plot(h,tcBlobsInProfile,'--','color',seudo.pickColor('SEblobs'))
    axis(h,'tight')
    set(h,'xtick',[])
    if all(tcBlobsAll==0), ylim(h,[0 1]),end
    axis(h,'on')
    
    
    % costs
    h = handles.axesTransSubplots(5,tt);
    cla(h)
    tK = extras.lsqCosts(framesToPlot);
    tR = extras.blobCosts(framesToPlot);
    plotK = tK <= tR;
    iK = bwlabel(plotK);
    iR = bwlabel(~plotK);
    %     for bb = 1:max(iK)
    %         %patch([find(iK==bb); flipud(find(iK==bb))],[tK(iK==bb); flipud(tR(iK==bb))],'k','facealpha',.3,'lineStyle','none','parent',h)
    %         %patch([find(iK==bb); find(iK==bb,1,'last')+.5; find(iK==bb,1,'first')-.5],[tK(iK==bb); 0; 0],'k','facealpha',.3,'lineStyle','none','parent',h)
    %         hold(h,'on')
    %     end
    %     for bb = 1:max(iR)
    %         %patch([find(iR==bb); flipud(find(iR==bb))],[tR(iR==bb); flipud(tK(iR==bb))],'r','facealpha',.5,'lineStyle','none','parent',h)
    %         %patch([find(iR==bb); find(iR==bb,1,'last')+.5; find(iR==bb,1,'first')-.5],[tR(iR==bb); 0; 0],'r','facealpha',.5,'lineStyle','none','parent',h)
    %     end
    plot(h,[1 length(framesToPlot)],[0 0],'-','color',[1 1 1]*.5)
    hold(h,'on')
    plot(h,find(~iK),tK(iK==0),'o','markersize',5,'color',seudo.pickColor('true'))
    plot(h,find(~iR),tR(iR==0),'o','markersize',5,'color',seudo.pickColor('false'))
    plot(h,find(iK),tK(iK>0),'.','markersize',15,'color',seudo.pickColor('true'))
    plot(h,find(iR),tR(iR>0),'.','markersize',15,'color',seudo.pickColor('false'))
    axis(h,'on')
    ylim(h,[min([tK; tR]) max([tK; tR])])
    
    axis(h,'tight')
    
    
    
    % axis limits
    %linkaxes(handles.axesTransSubplots(3:5,tt),'x')
    for pp = 3:5
        %xlim(handles.axesTransSubplots(pp,tt),[0 length(framesToPlot)+1])
    end
    
end




for aa = reshape(handles.axesTransSubplots(:,tt+1:end),1,[])
    cla(aa)
    axis(aa,'off')
    title(aa,'')
end


    
for tt = 1:length(transToPlot)
    
    % space plots
    %linkaxes(handles.axesTransSubplots([1 2 ],tt))
end


% highlight which transients are plotted
h = handles.axesWhichTransients;
delete(get(h,'child'))
plot(h,find(ismember(handles.transientSortOrder,handles.whichTransients)),0,'.','color','k','markerSize',10)
ylim(h,[-1 1])
xlim(h,xlim(handles.axesTransientFractions))
%axis(h,'off')


if DEBUG_TEXT, toc, end








% UPDATE FRAME PLOTS

function handles = updateFramePlots(hObject,handles)


if DEBUG_TEXT, fprintf('frames %s for transient %d from cell %d\n',makeCommaDelimitedString(handles.whichFrames),handles.transientForFrames,handles.thisCell); end


se = handles.se;
tcDefault =  handles.se.(handles.tcStructSpecDefault{1})(handles.tcStructSpecDefault{2}); %.transientInfo(handles.thisCell).times
tcStructSeudo = handles.se.(handles.tcStructSpecSeudo{1})(handles.thisSeudoNumber);
extras = tcStructSeudo.extras{handles.thisCell};
params = tcStructSeudo.params{handles.thisCell};


ti = tcDefault.transientInfo(handles.thisCell);




% get transient movie

% load from handles struct if already created
if isfield(handles,'transientMovie') && ...
        isfield(handles.transientMovie,'whichCell') && handles.transientMovie.whichCell == handles.thisCell && ...
        isfield(handles.transientMovie,'whichTransient') && handles.transientMovie.whichTransient == handles.transientForFrames
    
    % retrieve data from struct
    M = handles.transientMovie.data;
    addY = handles.transientMovie.addY;
    nY = handles.transientMovie.nY;
    nX = handles.transientMovie.nX;
    
    
    % ensure not too many frames are requested
    handles.whichFrames = intersect(handles.whichFrames,1:size(M,3));
    
    % get relevant farmes
    Mplot = seudo.makeThumbnailMatrix(reshape(M(:,:,handles.whichFrames,:),nY+addY,nX,[]),'nY',4);
    
    % swap in image data (much faster than using imagesc and setting all properties, e.g. colormap)
    h = handles.axesPlotAreaMovie;
    set(get(h,'child'),'CData',Mplot)
    
else
    
    % otherwise create it from scratch
    
    
    % identify time subset
    transFrames = ti.times(handles.transientForFrames,1):ti.times(handles.transientForFrames,2);
    
    % ensure not too many frames are requested
    handles.whichFrames = intersect(handles.whichFrames,1:length(transFrames));
    
    
    % identify space subset that was analyzed
    win = handles.thisCellWindowGlobal;
    nY = diff(win(1:2))+1;
    nX = diff(win(3:4))+1;

    
    % movie frames
    F = se.getFrames(transFrames,win(1):win(2),win(3):win(4));
    
    % apply downsampling
    dsT = params.dsTime;
    F = convn(F,cat(3,ones(1,1,dsT),zeros(1,1,dsT-1))/dsT,'same');
    %F = convn(F,ones(1,1,dsT)/dsT,'same');
    
    
    % least squares fit
    whichCells = handles.thisCell;
    L = reshape(reshape(se.profiles(win(1):win(2),win(3):win(4),whichCells),[],length(whichCells))*tcStructSeudo.tcLSQ(transFrames,whichCells)',size(F));
    %L = reshape(P(:)*tcLSQ(transFrames)',size(theFrames));
    
    % seudo fit
    S = reshape(reshape(se.profiles(win(1):win(2),win(3):win(4),whichCells),[],length(whichCells))*tcStructSeudo.tc(transFrames,whichCells)',size(F));
    
    % blobs
    if params.saveBlobTimeCourse
        % make image that shows where each blob was located
        B = reshape(extras.tcBlobs(transFrames,:)',nY,nX,[]);
        B = convn(B, extras.one_blob, 'same');
    else
        % make solid color block that shows overall block amplitude
        B = repmat(permute(extras.tcBlobs(transFrames,1),[2 3 1]),nY,nX,1)/nY/nX*100;
    end
    
    
    
    
    % assemble into large matrix
    M = cat(4,F,L,S,B);
    
    % ensure it's as tall as it is wide
    addY = size(M,2) - size(M,1);
    if addY > 0
        addTop = round(addY/2);
        addBottom = addY-addTop;
        M = padarray(M,[addTop 0 0],nan,'pre');
        M = padarray(M,[addBottom 0 0],nan,'post');
    else
        addY = 0;
    end
    
    % reshape
    Mplot = seudo.makeThumbnailMatrix(reshape(M(:,:,handles.whichFrames,:),nY+addY,nX,[]),'nY',4);

    % both in one line (might be faster?)
    %M = seudo.makeThumbnailMatrix(cat(3,F(:,:,handles.whichFrames),L(:,:,handles.whichFrames),S(:,:,handles.whichFrames),B(:,:,handles.whichFrames)),'nY',4);
    
    
    % get max and min pixels values (smooth movie first to prevent spurious values from skewing the color scale)
    allVals = cat(3,convn(F,ones(3)/9,'same'),L,S,B);
    allVals = allVals(:);
    colorLim = [min(allVals) max(allVals)];
    
    
    % store in handles struct
    handles.transientMovie.data = M;
    handles.transientMovie.addY = addY;
    handles.transientMovie.nY = nY;
    handles.transientMovie.nX = nX;
    handles.transientMovie.whichCell = handles.thisCell;
    handles.transientMovie.whichTransient = handles.transientForFrames;
    
    guidata(hObject,handles)
    
    
    
    % creat plot
    h = handles.axesPlotAreaMovie;
    cla(h)
    imagesc(Mplot,'parent',h)
    axis(h,'image')
    axis(h,'off')
    %set(h,'xtick',[],'ytick',[])
    caxis(h,colorLim)

    
end



% update the patch object
h = handles.axesPatch;

% position of transient being plotted
posAbove = get(handles.axesTransSubplots(end,handles.transientForFrames == handles.whichTransients),'position');
posBelow = get(handles.axesPlotAreaMovie);
posHere = get(h,'position');

% identify position of above relative to here
aboveStart = (posAbove(1)-posHere(1))/posHere(3);
aboveEnd = (posAbove(1)+posAbove(3)-posHere(1))/posHere(3);
aboveWidth = aboveEnd - aboveStart;

% identify position of frames above 
transPlotFraction = (handles.whichFrames([1 end])-1)/(size(handles.transientMovie.data,3)-1);

% identify aspect ratio (width / height) or image and axes
imAspect = get(handles.axesPlotAreaMovie,'plotboxaspectratio');
imAspect = imAspect(1)/imAspect(2);
set(handles.axesPlotAreaMovie,'units','pixels')
axesAspect = get(handles.axesPlotAreaMovie,'position');
set(handles.axesPlotAreaMovie,'units','normalized')
axesAspect = axesAspect(3)/axesAspect(4);

if imAspect > axesAspect
    % if image is longer, it takes up full plot, so use whole thing
    bottomLeft = 0;
    bottomRight = 1;
else
    % if axes is longer, calculate subset
    plotFraction = imAspect / axesAspect;
    bottomLeft = (1-plotFraction)/2;
    bottomRight = 1-(1-plotFraction)/2;
end
    


transStartPoint = transPlotFraction(1)*aboveWidth + aboveStart;
transEndPoint = transPlotFraction(2)*aboveWidth + aboveStart;


topKink = 0.92;
bottomKink = 0.85;

cla(h)
%plot(h,[aboveStart aboveEnd],[0 0],'*-')
%hold(h,'on')
%plot(h,[transStartPoint transEndPoint],[0 0],'*-')
patch('xdata',[transStartPoint transEndPoint transEndPoint  bottomRight bottomRight bottomLeft bottomLeft transStartPoint transStartPoint],...
    'ydata',[1 1 topKink bottomKink 0 0 bottomKink topKink 1],'parent',h,'facecolor',[1 1 1]*.5,'linestyle','none')
%hold(h,'on')
%patch('xdata',[transStartPoint transEndPoint  bottomRight bottomLeft transStartPoint ],...
%    'ydata',[topKink topKink bottomKink bottomKink topKink],'parent',h,'facecolor',[1 1 1]*.95,'linestyle','none')
xlim(h,[0 1])
ylim(h,[0 1])
axis(h,'off')

%[posAbove(1) posHere(3)]



% USER INTERFACE ELEMENTS













% MOVE CELL SLIDER

function sliderCellNumber_Callback(hObject, eventdata, handles)

% get slider value
handles.thisCell = round(get(hObject,'Value'));

% set editable text box
set(handles.editTextCellNumber,'String',num2str(handles.thisCell))

% change cell
handles = setNewCell(hObject,handles);

% store result
guidata(hObject,handles)


function sliderCellNumber_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end





% EDIT CELL NUMBER TEXT

function editTextCellNumber_Callback(hObject, eventdata, handles)

% get user input
newCellNumber = round(str2num(get(hObject,'String')));

% if not valid, ignore
if isempty(newCellNumber)
    set(hObject,'string',num2str(handles.thisCell))
    return
end

% ensure in range
newCellNumber = min(max(1,newCellNumber),handles.se.nCells);

% store value
handles.thisCell = newCellNumber;

% update slider
set(handles.sliderCellNumber,'Value',handles.thisCell)

% ensure editable text box is clear
set(hObject,'String',num2str(handles.thisCell))

% change cell
handles = setNewCell(hObject,handles);

% store result
guidata(hObject,handles)


function editTextCellNumber_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% SELECT SEUDO NUMBER WITH POPUP

function popupmenuSeudoNumber_Callback(hObject, eventdata, handles)

setNewSeudoNumber(get(hObject,'parent'),handles);

function popupmenuSeudoNumber_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% EDIT TRANSIENT COUNT TEXT

function editNtrans_Callback(hObject, eventdata, handles)

% get value
newTransCount = str2double(get(handles.editNtrans,'string'));

% ensure reasonable value
if ~isnumeric(newTransCount) || isempty(newTransCount)
    set(handles.editNtrans,'string',num2str(handles.nTransients))
    return
end
if length(newTransCount) > 1, newTransCount = newTransCount(1); end
newTransCount = min(newTransCount,MAX_TRANSIENT_PLOTS);
newTransCount = max(newTransCount,1);

% update value
handles.nTransients = newTransCount;

% display value
set(handles.editNtrans,'string',num2str(handles.nTransients))

% plot new transients
handles = setNewTransientCount(get(hObject,'parent'),handles);

% store handles
guidata(hObject,handles)


function editNtrans_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% MOVE TRANSIENT SLIDER

function sliderTransients_Callback(hObject, eventdata, handles)

% get value from slider
firstTransient = round(get(handles.sliderTransients,'value'));

% update slider
set(handles.sliderTransients,'Value',firstTransient)

% set transient count
handles.whichTransients = handles.transientSortOrder(firstTransient + (0:handles.nTransients-1));

% update text field
set(handles.editTransientsToPlot,'string',makeCommaDelimitedString(handles.whichTransients))

% plot new transients
handles = setNewTransients(hObject,handles);

% store in handles
guidata(hObject,handles)

function sliderTransients_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% EDIT TRANSIENT NUMBERS

function editTransientsToPlot_Callback(hObject, eventdata, handles)

% parse input
wt = textscan(get(handles.editTransientsToPlot,'string'),'%f',inf,'delimiter',',');
wt = wt{1};

% if nothing parsable, return to previous value
if isempty(wt)
    set(handles.editTransientsToPlot,'string',makeCommaDelimitedString(handles.whichTransients))
    return
end

% make sure it's only natural numbers in the range, and no duplicates
wt = round(wt);
transCount = size(handles.se.(handles.tcStructSpecDefault{1})(handles.tcStructSpecDefault{2}).transientInfo(handles.thisCell).times,1);
wt = wt(wt>0 & wt <= transCount);
[~,toKeep] = unique(wt);
wt = wt(sort(toKeep));


% if more transients than allowed, discard excess
wt = wt(1:min(length(wt),MAX_TRANSIENT_PLOTS));


% if transient count changed, update transient count
if length(wt) ~= handles.nTransients
    handles.nTransients = length(wt);
    handles = setNewTransientCount(get(hObject,'parent'),handles);
    set(handles.editNtrans,'string',num2str(handles.nTransients))
end


% note new transients
handles.whichTransients = wt;

% set string
set(handles.editTransientsToPlot,'string',makeCommaDelimitedString(wt))

% plot new transients
handles = setNewTransients(hObject,handles);

% store in handles
guidata(hObject,handles)


function editTransientsToPlot_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% MOVE SLIDER FOR TRANSIENT TO SHOW FRAMES

function sliderTransientForFrames_Callback(hObject, eventdata, handles)

handles.transientForFrames = handles.whichTransients(get(hObject,'Value'));

handles = setNewTransientForFrames(hObject,handles);

guidata(hObject,handles)


function sliderTransientForFrames_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end







% EDIT FRAME COUNT TEXT

function editNframes_Callback(hObject, eventdata, handles)


% get value
newFramesCount = str2double(get(handles.editNframes,'string'));

% ensure reasonable value
if ~isnumeric(newFramesCount) || isempty(newFramesCount)
    set(handles.editNframes,'string',num2str(handles.nFrames))
    return
end
if length(newFramesCount) > 1, newFramesCount = newFramesCount(1); end
newFramesCount = min(newFramesCount,MAX_FRAME_COUNT);
newFramesCount = max(newFramesCount,1);

% update value
handles.nFrames = newFramesCount;

% insert value
set(handles.editNframes,'string',num2str(handles.nFrames))


% plot new frames
handles = setNewTransientForFrames(hObject,handles);

% store handles
guidata(hObject,handles)



function editNframes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% MOVE FRAMES SLIDER

function sliderMovieFrames_Callback(hObject, eventdata, handles)

% get value from slider
firstFrame = round(get(handles.sliderMovieFrames,'value'));

% set transient count
handles.whichFrames = firstFrame + (0:handles.nFrames-1);

% plot new transients
handles = updateFramePlots(hObject,handles);

% store in handles
guidata(hObject,handles)



function sliderMovieFrames_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end





% CHANGE TRANSIENT SORT ORDER

function popupmenuSortTransients_Callback(hObject, eventdata, handles)

handles = setNewCell(hObject,handles);

guidata(hObject,handles)

function popupmenuSortTransients_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% RECOMPUTE CONTAMINATION SCATTER

function editScatterContamParams_Callback(hObject, eventdata, handles)
handles = computeScatterContam(hObject,handles);

% update to be sure cell highlighting is preserved
handles = updateCellPlots(hObject,handles);

guidata(hObject,handles)


function editScatterContamParams_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% SHOW PARAMETERS

function pushbuttonPrintParams_Callback(hObject, eventdata, handles)

% get parameters
paramsHere = handles.se.(handles.tcStructSpecSeudo{1})(handles.thisSeudoNumber).params{handles.thisCell};
paramsDefault = handles.se.estimateTimeCoursesWithSEUDO('showDefaults',1);

% remove fields from both that don't affect the result
fldsToRemove = {'whichCells','whichFrames','verbose','validateParams','saveBlobTimeCourse','saveOtherCellTimeCourses','showDefaults'};
for ff = 1:length(fldsToRemove)
    if isfield(paramsHere,fldsToRemove{ff})
        paramsHere = rmfield(paramsHere,fldsToRemove{ff});
    end
    if isfield(paramsDefault,fldsToRemove{ff})
        paramsDefault = rmfield(paramsDefault,fldsToRemove{ff});
    end
end

% show all parameters used to generate the plotted result
fprintf('parameters for seudo set %d, cell %d (code version %s):\n\n',...
    handles.thisSeudoNumber,handles.thisCell,num2str(paramsHere.version))
paramsHere = rmfield(paramsHere,'version');
disp(paramsHere)


% show parameters that differ from from the current default


fprintf('parameters that differ from current default (code version %s):\n\n',num2str(paramsDefault.version))
paramsDefault = rmfield(paramsDefault,'version');

% identify differences
fn = fieldnames(paramsHere);
paramsDiff = struct;
for ff = 1:length(fn)
    if ~isfield(paramsDefault,fn{ff})
        fprintf('STRANGE: field %s is not found in current default parameters. Different SEUDO versions?\n',fn{ff})
        continue
    end
    
    if ~isequaln(paramsHere.(fn{ff}),paramsDefault.(fn{ff}))
        paramsDiff.(fn{ff}) = paramsHere.(fn{ff});
        %diffString = sprintf('%s''%s'',
    end
end

disp(paramsDiff)




% RELOAD ALL GRAPHICS

function pushbuttonReload_Callback(hObject, eventdata, handles) %#ok<*DEFNU,*INUSL>

figHandle = get(hObject,'parent');

handles.transientMovie.whichCell = 0;

handles = updateOverviewPlots(hObject,handles);

handles = makeNewTransientSubplots(figHandle,handles);
handles = setNewCell(figHandle,handles);
%handles = setNewTransientCount(figHandle,handles);
%handles = setNewSeudoNumber(figHandle,handles);


set(handles.editTextCellNumber,'String',num2str(handles.thisCell))
set(handles.sliderCellNumber,'Value',handles.thisCell)

guidata(hObject,handles)




% SHOW/HIDE CELL NUMBERS

function checkboxShowCellNumbers_Callback(hObject, eventdata, handles)


if get(handles.checkboxShowCellNumbers,'value')
    set(handles.axesScatterEvalNumbers,'visible','on')
    set(handles.axesScatterContamNumbers,'visible','on')
    set(handles.axesProfileCOMs,'visible','on')
else
    set(handles.axesScatterEvalNumbers,'visible','off')
    set(handles.axesScatterContamNumbers,'visible','off')
    set(handles.axesProfileCOMs,'visible','off')
end








% COMPUTE CONTAMINATION SCATTER

function handles = computeScatterContam(hObject,handles)

% show wait indicator
set(handles.textWaitScatter,'visible','on')
drawnow


% extract user specified parameters
eval(sprintf('userParams = struct(%s);',get(handles.editScatterContamParams,'string')));


% run code with user params
transientMetricsByCell = handles.se.autoClassifyTransients(...
    {handles.tcStructSpecDefault{1},handles.tcStructSpecDefault{2}},...
    'saveResults',false,userParams);


% hide wait indicator
set(handles.textWaitScatter,'visible','off')



% plot
h = handles.axesScatterContam;
cla(h)
x = [transientMetricsByCell.cfrac];
y = [transientMetricsByCell.rfrac];
plot(h,x,y,'.','color','k')
xlabel(h,'contamination rate')
%ylabel(h,'unexplained structure residual')
%ylabel(h,'unexplained transient profile variance')
ylabel(h,'contamination severity')
xlim(h,[0 1])
ylim(h,[0 1])
set(h,'xtick',0:.2:1,'ytick',0:.2:1)
hold(h,'on')


% add cell number labels
handles.axesScatterContamNumbers = nan(handles.se.nCells,1);
for cc = 1:handles.se.nCells
    handles.axesScatterContamNumbers(cc) = text(x(cc),double(y(cc)),[' ' num2str(cc)],'color',[1 1 1]*.2,'fontsize',8,'parent',h);
end

if get(handles.checkboxShowCellNumbers,'value')
    set(handles.axesScatterContamNumbers,'visible','on')
else
    set(handles.axesScatterContamNumbers,'visible','off')
end


% store results
handles.scatterContamX = x;
handles.scatterContamY = y;







% HELPER FUNCTIONS

function theStr = makeCommaDelimitedString(vals)

vals = vals(~isnan(vals));
theStr = sprintf('%d,',vals);
theStr = theStr(1:end-1);








% CONSTANTS

function maxCount = MAX_TRANSIENT_PLOTS()
maxCount = 5;

function maxCount = MAX_FRAME_COUNT()
maxCount = 40;

function tf = DEBUG_TEXT()
%tf = true;
tf = false;

function col = CELL_HIGHLIGHT_COLOR()
col = [0 .8 .8];
