function varargout = seudoReviewParamSearch(varargin)
% SEUDOREVIEWPARAMSEARCH MATLAB code for seudoReviewParamSearch.fig
%      SEUDOREVIEWPARAMSEARCH, by itself, creates a new SEUDOREVIEWPARAMSEARCH or raises the existing
%      singleton*.
%
%      H = SEUDOREVIEWPARAMSEARCH returns the handle to a new SEUDOREVIEWPARAMSEARCH or the handle to
%      the existing singleton*.
%
%      SEUDOREVIEWPARAMSEARCH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEUDOREVIEWPARAMSEARCH.M with the given input arguments.
%
%      SEUDOREVIEWPARAMSEARCH('Property','Value',...) creates a new SEUDOREVIEWPARAMSEARCH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before seudoReviewParamSearch_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to seudoReviewParamSearch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help seudoReviewParamSearch

% Last Modified by GUIDE v2.5 18-Jul-2018 13:40:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @seudoReviewParamSearch_OpeningFcn, ...
    'gui_OutputFcn',  @seudoReviewParamSearch_OutputFcn, ...
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


% --- Executes just before seudoReviewParamSearch is made visible.
function seudoReviewParamSearch_OpeningFcn(hObject, eventdata, handles, varargin)
% launch GUI to examine results of calling seudoObject.paramSearch
%
% USAGE: seudoReviewParamSearch(resultsStruct)
%
% PARAMETERS:
%
%   whichTransients - which transients to use
%                       T-length logical vector specifying which to include
%                       vector of which elements in results.pickCellTransients
%                       2-column matrix, each row [ cellID transientID ]
%


% clean up from previous if figure was not closed
if isfield(handles,'hLists')
    for hh = 1:length(handles.hLists)
        delete(handles.hLists{hh})
    end
end
if isfield(handles,'toDelete')
    delete(handles.toDelete)
end


handles.results = varargin{1};


% parse inputs
p = inputParser;
p.addParamValue('whichTransients',[]);
parse(p,varargin{2:end});




if isempty(p.Results.whichTransients)
    % default to all transients
    handles.whichTransients = 1:size(handles.results.pickedCellTransients,1);
else
    
    % specify subset
    if islogical(p.Results.whichTransients)
        % logical vector
        if length(p.Results.whichTransients) ~= size(handles.results.pickedCellTransients,1)
            error('incorrect logical transient specification: needed %d elements, only %d provided',...
                size(handles.results.pickedCellTransients,1),length(p.Results.whichTransients))
        end
        handles.whichTransients = find(p.Results.whichTransients);
        
    elseif isvector(p.Results.whichTransients) && length(p.Results.whichTransients) <= size(handles.results.pickedCellTransients,1)
        % specified vector
        handles.whichTransients = p.Results.whichTransients;
        
    elseif size(p.Results.whichTransients,2) == 2
        % list of transients
        handles.whichTransients = find(ismember(handles.results.pickedCellTransients(:,1:2),p.Results.whichTransients,'rows'));
    else
        error('transient specification not recognized')
    end
    
    
end




% set name
% name, time, # param sets, # variables, # transients
set(hObject,'name',sprintf('%d parameter sets applied to %d true & %d false transients from %s (computed %s)',...
    length(handles.results.paramSets),...
    sum(handles.results.pickedCellTransients(handles.whichTransients,3)==1),sum(handles.results.pickedCellTransients(handles.whichTransients,3)==0),...
    handles.results.name,datestr(handles.results.startTime,'mmm dd,HH:MMpm')))


% validate input
%  COMING SOON!


% populate popup buttons
set(handles.popupROC,'string',[handles.results.nonSingletonParamNames(:); {'none'}])
set(handles.popupX,'string',handles.results.nonSingletonParamNames)
set(handles.popupY,'string',handles.results.nonSingletonParamNames)
set(handles.popupZ,'string',handles.results.nonSingletonParamNames)




% get name and possible values of each parameter
paramCodex = handles.results.paramCodex;
paramNames = handles.results.paramNames;
for pp=1:length(paramNames)
    thisParamIndex = find(strcmp(paramNames{pp},handles.results.paramNames));
    % identify unique values, and an example parameter set for each
    [paramTableValues,exampleSets] = unique(paramCodex(:,thisParamIndex));
    % exclude cases where default value was used
    exampleSets = exampleSets(paramTableValues~=0);
    paramVals{pp} = cellfun(@(x)valueToString(x.(paramNames{pp})),handles.results.paramSets(exampleSets),'uniformoutput',0);
end



% show singleton params and their values
singletonParamString = '';
for qq = 1:length(handles.results.singletonParamNames)
    pp = find(strcmp(handles.results.singletonParamNames{qq},paramNames));
    singletonParamString = [singletonParamString sprintf('%s: %s\n',paramNames{pp},valueToString(paramVals{pp}{1}))];
end


otherParamsString = '';
for qq = 1:length(handles.results.nonSingletonParamNames)
    pp = find(strcmp(handles.results.nonSingletonParamNames{qq},paramNames));
    fym = cellfun(@valueToString,paramVals{pp},'uniformoutput',0);
    otherParamsString = [otherParamsString paramNames{pp} ': ' sprintf('%s, ',fym{:}) sprintf('\n')];
end

%set(handles.textShowParams,'String',otherParamsString)


% first list shows the singleton parameters
set(handles.listboxSingletons,'string',singletonParamString(1:end-1),'max',0,'min',1)




% create lists of parameter values

posFig = getpixelposition(hObject);

%set(hObject,'units','pixels')

set(handles.axesListboxes,'units','pixels')
pos = get(handles.axesListboxes,'position');
set(handles.axesListboxes,'visible','off')

H = makeSubplots(hObject,length(handles.results.nonSingletonParamNames),1,.1,0,pos./posFig([3 4 3 4]));
handles.toDelete = H;


handles.hLists = {};
for hh = 1:length(handles.results.nonSingletonParamNames)
    
    h = H(hh);
    axis(h,'off')
    
    % needs to be specified in pixels
    set(h,'units','pixels')
    
    % create listbox
    handles.hLists{hh} = uicontrol('Parent',hObject,'Style','listbox',...
        'Position',get(H(hh),'position'),'Callback', @list_callback,...
        'fontsize',12);
    
    % subsequent lists are for non-singleton parameters
    
    % get number (in the list of parameters)
    pp = find(strcmp(handles.results.nonSingletonParamNames{hh},paramNames));
    
    
    if any(paramCodex(:,pp)==0)
        thisString = sprintf('[default]\n');
    else
        thisString = '';
    end
    
    fym = cellfun(@valueToString,paramVals{pp},'uniformoutput',0);
    thisString = [thisString sprintf('%s\n',fym{:})];
    codexValues = unique(paramCodex(:,pp));
    
    set(handles.hLists{hh},'string',thisString(1:end-1),'min',0,'max',2,'value',1:length(codexValues))
    title(h,paramNames{pp},'fontsize',12)
    
    handles.codexValues{hh} = codexValues;
    
end



set(handles.listboxSingletons,'value',[],'min',0,'max',2)



% Choose default command line output for seudoReviewParamSearch
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% make plot
updateROCplot(handles)







% --- Outputs from this function are returned to the command line.
function varargout = seudoReviewParamSearch_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;







% CHOOSE PARAMTER TO COLOR IN POPUP MENU

function popupROC_Callback(hObject, eventdata, handles)
updateROCplot(handles)

function popupROC_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% OPEN MATRIX VIEWER FIGURE

function pushbuttonGenerateMatrices_Callback(hObject, eventdata, handles)
generateMatrices(handles)

function popupX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function popupZ_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function generateMatrices(handles)



% identify columns of vars to plot
paramNameX = handles.results.nonSingletonParamNames{get(handles.popupX,'value')};
paramNameY = handles.results.nonSingletonParamNames{get(handles.popupY,'value')};
paramNameZ = handles.results.nonSingletonParamNames{get(handles.popupZ,'value')};
paramIndexX = find(strcmp(paramNameX,handles.results.paramNames));
paramIndexY = find(strcmp(paramNameY,handles.results.paramNames));
paramIndexZ = find(strcmp(paramNameZ,handles.results.paramNames));



paramCodex = handles.results.paramCodex;


[paramValsX,paramExampleX] = unique(paramCodex(:,paramIndexX));
[paramValsY,paramExampleY] = unique(paramCodex(:,paramIndexY));
[paramValsZ,paramExampleZ] = unique(paramCodex(:,paramIndexZ));


npX = length(paramValsX);
npY = length(paramValsY);
npZ = length(paramValsZ);



valNamesX = cellfun(@(x)x.(paramNameX),handles.results.paramSets(paramExampleX),'uniformoutput',0);
valNamesY = cellfun(@(x)x.(paramNameY),handles.results.paramSets(paramExampleY),'uniformoutput',0);
valNamesZ = cellfun(@(x)x.(paramNameZ),handles.results.paramSets(paramExampleZ),'uniformoutput',0);



[meanKeptFractionTrue,meanKeptFractionFalse] = computeKeptFraction(handles);


% color in the points

dataColored = nan(npY,npX,3,npZ);

%ptToCol = @(x,y)hsv2rgb(cat(3,(1-(x+y)/2)*2/3+1/3,((y-x+1)/2).^3*.7+.3,ones(size(x))));
ptToCol = @(x,y)hsv2rgb(cat(3,(1-(x+y)/2)*2/3+1/3,max(1-x,y).^3*.8+.2,ones(size(x))));


for xx = 1:npX
    for yy = 1:npY
        for zz=1:npZ
            
            theLines = paramCodex(:,paramIndexX)==paramValsX(xx) & ...
                paramCodex(:,paramIndexY)==paramValsY(yy) & ...
                paramCodex(:,paramIndexZ)==paramValsZ(zz);
            
            dataColored(yy,xx,:,zz) = ptToCol(...
                nanmean(meanKeptFractionTrue(theLines)),nanmean(meanKeptFractionFalse(theLines)));
            
        end
    end
end



figure(103);clf

H = makeSubplots(gcf,npZ,1,.1,.3,[.2 0 .8 1]);



% plot legend showing colormap

pos = get(H(1),'position');
h = subplot('position',[.03 pos(2) .14 pos(4)]);

bs = 0.05;
x=bs/2:bs:(1-bs/2);
[X,Y] = meshgrid(x,x);
imagesc(ptToCol(X,Y),'xdata',[bs/2 1-bs/2],'ydata',[bs/2 1-bs/2],'parent',h)
axis(h,'xy'); axis(h,'image')
hold(h,'on')
plot(h,meanKeptFractionFalse,meanKeptFractionTrue,'ok','markersize',5)
xlabel(h,'false transients kept')
ylabel(h,'true transients kept')


% plot each matrix
for zz=1:npZ
    h = H(zz);
    imagesc(dataColored(:,:,:,zz),'parent',h)
    title(h,sprintf('%s = %s',paramNameZ,valueToString(valNamesZ{zz})))
    
    set(h,'xtick',1:npX,'ytick',1:npY,'xticklabel',cellfun(@valueToString,valNamesX,'uniformoutput',0))
    if zz == 1
        xlabel(h,paramNameX)
        set(h,'yticklabel',cellfun(@valueToString,valNamesY,'uniformoutput',0))
        ylabel(h,paramNameY)
    else
        set(h,'yticklabel','')
    end
end


% --- Executes on button press in pushbuttonSingleCells.
function pushbuttonSingleCells_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSingleCells (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% 
[meanKeptFractionTrue,meanKeptFractionFalse,keptFractions] = computeKeptFraction(handles);



figure
trueTrans = find(handles.results.pickedCellTransients(handles.whichTransients,3)==1)';
falseTrans = find(handles.results.pickedCellTransients(handles.whichTransients,3)==0)';
transToPlot = {trueTrans,falseTrans};
for pp=1:2
    subplot(1,2,pp)
    switch pp
        case 1
            meanCurve = meanKeptFractionTrue;
            tl = 'true transients';
        case 2
            meanCurve = meanKeptFractionFalse;
            tl = 'false transients';
    end
    
    % sort points according to mean curve
    [~,po] = sort(meanCurve);
    
    % pick colors
    cols = distinguishable_colors(length(transToPlot{pp}),'w');
    %cols = linspecer(length(transToPlot{pp}),'w')
    
    % plot points and fit curve from each cell
    for cc = 1:length(transToPlot{pp})
        
        % points
        x = meanCurve(po);
        y = keptFractions(transToPlot{pp}(cc),po);
        plot(x,y,'.','color',cols(cc,:))
        hold on
        
        % fit curve
        switch 1
            case 1 % polynomial
                p = polyfit(x,y,5);
                x = linspace(min(x),max(x),50);
                plot(x,polyval(p,x),'-','color',cols(cc,:),'linewidth',2)
                
            case 2 % powers
                f = fit(x(:),y(:),'power1');
                plot(x,f(x),'-','color',cols(cc,:))
                
            case 3 % powers 
                B = lsqcurvefit(@(b,X)X.^abs(1/b),1,x,y);
                plot(x,x.^abs(1/B),'-','color',cols(cc,:))
        end
        
    end
    
    % add guide lines
    plot([0 1],[0 1],'k-')
    plot([0 1 1 0 0],[0 0 1 1 0],'k-')
    axis image
    title(tl)
    xlabel('mean fraction kept')
    ylabel('fraction kept')
end










% UPDATE PLOT WHEN USER INTERFACE ELEMENTS ARE CHANGED

function list_callback(src,~)
updateROCplot(guidata(get(src,'parent')))

function checkboxShowNum_Callback(hObject, eventdata, handles)
updateROCplot(handles)

function pushbuttonSelectAll_Callback(hObject, eventdata, handles)
% select all rows
for hh=1:length(handles.hLists)
    set(handles.hLists{hh},'value',1:length(handles.codexValues{hh}));
end
% update plot
updateROCplot(handles)




% DISALLOW SELECTIONS IN THE LIST OF SINGLETON PARAMETERS

function listboxSingletons_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'min',0,'max',2)



function listboxSingletons_Callback(hObject, eventdata, handles)
set(hObject,'value',[],'min',0,'max',2)
%uicontrol(handles.popupROC)

jPeer = get(handles.figure1, 'JavaFrame');
jPeer.getAxisComponent.requestFocus







% EVALUATE SEUDO TIME COURSES

function [meanKeptFractionTrue,meanKeptFractionFalse,keptFractions] = computeKeptFraction(handles)


wt = handles.whichTransients;

switch 1
    case 1
        %% old way
        
        
        nParamSets = length(handles.results.paramSets);
        nTransients = length(wt);
        
        
        % measure how well each transient preserved after running SEUDO
        
        keptFractions = nan(nTransients,nParamSets);
        transAreas = nan(nTransients,nParamSets);
        
        switch get(handles.popupTransientMeasure,'value')
            case 1 
                args = {'method','area'};
            case 2 
                args = {'method','tall span','threshSigTrans',10,'threshTall',0.2};
            case 3
                args = {'method','tall span threshold','threshSigTrans',10,'threshTall',0.2,'minDuration',0.5};
        end
        
        
        evalFcn = seudo.makeTimeCourseEvalFunction(args{:});
        
        
        
        for tt = 1:nTransients
            
            thisTransient = wt(tt);
            
            % get original trace
            origTrace = handles.results.tc(thisTransient).tcDefault.tc(:,1);
            
            for pp = 1:nParamSets
                % get seudo trace
                seudoTrace = handles.results.tc(thisTransient).tcSeudo(pp).tc;
                
                % compute measure
                keptFractions(tt,pp) = evalFcn(origTrace,seudoTrace); 
            end
                
            % compute area
            transAreas(tt,:) = sum(origTrace(origTrace>0));
            
        end
        
        
        
        % compute overall measure for each group of transients (true and false from each parameter set)
        
        isTrue = handles.results.pickedCellTransients(wt,3)==1;
        isFalse = handles.results.pickedCellTransients(wt,3)==0;
        
        switch get(handles.popupCombineTransients,'value')
            
            case 1 % average
                meanKeptFractionTrue = nanmean(keptFractions(isTrue,:),1);
                meanKeptFractionFalse = nanmean(keptFractions(isFalse,:),1);
                
            case 2 % average weighted by total area
                
                for pp = 1:nParamSets
                    meanKeptFractionTrue(1,pp) = nansum(keptFractions(isTrue,pp).*transAreas(isTrue,pp)) / sum(transAreas(isTrue,pp));
                    meanKeptFractionFalse(1,pp) = nansum(keptFractions(isFalse,pp).*transAreas(isFalse,pp)) / sum(transAreas(isFalse,pp));
                end
                
        end
        
        

        % total area preserved
        
        
    case 2
        
        %% new way
        
        
        
        meanKeptFractionTrue = nan(length(handles.results.paramSets),1);
        meanKeptFractionFalse = nan(length(handles.results.paramSets),1);
        
        % cycle through transients and param sets
        for pp = 1:length(handles.results.paramSets)
            
            allFracTrue = [];
            allFracFalse = [];
            
            for tt = 1:length(handles.results.tc)
                
                thisTransient = wt(tt);
            
                % make tcDefault
                tcDefault = handles.results.tc(thisTransient).tcDefault;
                % store classification
                handles.results.pickedCellTransients(thisTransient,3) = tcDefault.transientInfo.classification;
                
                % make tcNew (from SEUDO)
                tcNew = handles.results.tc(thisTransient).tcSeudo(pp).tc;
                
                % get evaluation
                [T,F] = evaluateNewTimeCourses(tcDefault,tcNew);
                
                allFracTrue = [allFracTrue; T];
                allFracFalse = [allFracFalse; F];
            end
            
            meanKeptFractionTrue(pp,1) = nanmean(allFracTrue);
            meanKeptFractionFalse(pp,1) = nanmean(allFracFalse);
        end
        
        
end





% MAKE A STRONG FROM A GIVEN PARAMETER VALUE
function str = valueToString(theVal)

if isnumeric(theVal)
    str = num2str(theVal);
elseif ischar(theVal)
    str = theVal;
elseif islogical(theVal)
    if theVal, str = 'true'; else str = false; end
else
    error('unrecognized type: %s',class(theVal))
end






% UPDATE PLOT

function updateROCplot(handles)


nParamSets = length(handles.results.paramSets);


[meanKeptFractionTrue,meanKeptFractionFalse] = computeKeptFraction(handles);



% identify points to plot based on which items are selected in the lists

paramCodex = handles.results.paramCodex;
plotSel = true(size(paramCodex,1),1);
for hh=1:length(handles.hLists)
    
    % identify which non-singleton param this is
    pp = find(strcmp(handles.results.nonSingletonParamNames{hh},handles.results.paramNames));

    % identify rows matching this selection
    sel = get(handles.hLists{hh},'value');
    rowsFromThisSelection = ismember(paramCodex(:,pp),handles.codexValues{hh}(sel));
    
    % update list of rows
    plotSel = plotSel & rowsFromThisSelection;
end



dimColor = @(c)1-0.15*(1-c);


% set up plot
h = handles.axesROC;
cla(h)
set(h,'box','on','tickdir','out','xtick',0:.2:1,'ytick',0:.2:1)
grid(h,'on')
hold(h,'on')
plot(h,[0 1],[0 1],'-','color',[1 1 1]*.5)
xlabel(h,'false transients kept')
ylabel(h,'true transients kept')
axis(h,'image')



% plot points
if get(handles.popupROC,'value') == length(get(handles.popupROC,'string'))
    % all the same color
    plot(h,meanKeptFractionFalse(~plotSel),meanKeptFractionTrue(~plotSel),'.','markersize',15,'color',dimColor([0 0 0]))
    plot(h,meanKeptFractionFalse(plotSel),meanKeptFractionTrue(plotSel),'.k','markersize',15)
    
    cla(handles.axesLegend)
    set(handles.axesLegend,'visible','off')
else
    % plot each set in a different color
    
    % identify # of plots
    thisParamName = handles.results.nonSingletonParamNames{get(handles.popupROC,'value')};
    thisParamIndex = find(strcmp(thisParamName,handles.results.paramNames));
    paramCodex = handles.results.paramCodex;
    
    possibleParamValues = unique(paramCodex(:,thisParamIndex));
    nPlots = length(possibleParamValues);
    
    % set colors
    %cols = parula(nPlots);
    cols = hsv2rgb([linspace(.1,1,nPlots)' ones(nPlots,1) linspace(.5,1,nPlots)']);
    %cols = distinguishable_colors(nPlots,'w');
    %cols = linspecer(nPlots,'w');
    
    % plot points with each value in a different color
    
    % plot points that don't match selection in a dimmed color
    for nn = 1:nPlots
        toPlot = paramCodex(:,thisParamIndex) == possibleParamValues(nn) & ~plotSel;
        plot(h,meanKeptFractionFalse(toPlot),meanKeptFractionTrue(toPlot),'.','markersize',15,'color',dimColor(cols(nn,:)))
    end
    
    % plot points matching selection in a bright color
    for nn = 1:nPlots
        toPlot = paramCodex(:,thisParamIndex) == possibleParamValues(nn) & plotSel;
        plot(h,meanKeptFractionFalse(toPlot),meanKeptFractionTrue(toPlot),'.','markersize',15,'color',cols(nn,:))
    end
    
    
    % show legend
    h = handles.axesLegend;
    cla(h)
    
    for nn = 1:nPlots
        exampleParamSet = find(paramCodex(:,thisParamIndex) == possibleParamValues(nn),1,'first');
        if isfield(handles.results.paramSets{exampleParamSet},handles.results.paramNames{thisParamIndex})
            theVal = handles.results.paramSets{exampleParamSet}.(handles.results.paramNames{thisParamIndex});
            theStr = valueToString(theVal);
        else
            theStr = '[default]';
        end
        plot(h,1,nn,'.','color',cols(nn,:),'markersize',15)
        hold(h,'on')
        text(2,nn,theStr,'color',cols(nn,:),'parent',h,'fontsize',18)
    end
    ylim(h,[0 nPlots+1])
    xlim(h,[0 6])
    set(h,'xtick',[],'ytick',[],'box','on')
    
    set(handles.axesLegend,'visible','on')
end



% add text labels
if get(handles.checkboxShowNum,'value')
    for ss=find(plotSel')
        text(meanKeptFractionFalse(ss),meanKeptFractionTrue(ss),sprintf('\n%d',ss),...
            'verticalalign','baseline','horizontalAlign','center','color','k','fontsize',18,'parent',handles.axesROC,'clipping','on')
    end
end




% indicate if none match
if ~any(plotSel)
    text(.5,.5,'no paramter sets match selection','horizontalAlign','center','parent',handles.axesROC,'fontsize',18)
end

























% previous verion
%
% function reviewParamSearch(results)
% 
% if seudo time course has new transients, basically ignore them when computing stats.
% though technically possible, they should be super rare
%
% seudo transients should be subset of original transients
% blob transients should be about equal to the difference
%
% for each parameter set:
%   for each transient in the original time course, compute for the SEUDO time course:
%       - integrated area
%       - peak
%       - mean
%       - correlation with original time course
%
%   get histogram of these stats (maybe just integrated area divided by orig time course area) for true and false transients
%
% later in visualization, show user the histograms for each param set
%
%
%
% 
% nParamSets = length(results.paramSets);
% nTransients = length(results.tc);
% 
% % validate input
% 
% 
% % get fraction of each transient preserved after running SEUDO
% 
% keptFractions = nan(nTransients,nParamSets);
% 
% for tt = 1:nTransients
%     
%     % original trace
%     origTrace = results.tc(tt).tcDefault.tc;
%     origSum = sum(origTrace);
%     
%     for ss = 1:nParamSets
%         % seudo trace
%         seudoTrace = results.tc(tt).tcSeudo(ss).tc;
%         keptFractions(tt,ss) = sum(seudoTrace)/origSum;
%     end
% end
% 
% 
% % compute mean transients k
% meanKeptFractionTrue = nanmean(keptFractions(results.pickedCellTransients(:,3)==1,:),1);
% meanKeptFractionFalse = nanmean(keptFractions(results.pickedCellTransients(:,3)==0,:),1);
% 
% % display results
% 
% 
% 
% %% plot trajectory for each cell
% 
% 
% figure;
% trueTrans = find(results.pickedCellTransients(:,3)==1)';
% falseTrans = find(results.pickedCellTransients(:,3)==0)';
% transToPlot = {trueTrans,falseTrans};
% for pp=1:2
%     subplot(1,2,pp)
%     switch pp
%         case 1
%             meanThing = meanKeptFractionTrue;
%             tl = 'true transients';
%         case 2
%             meanThing = meanKeptFractionFalse;
%             tl = 'false transients';
%     end
%     [~,po] = sort(meanThing);
%     for cc = transToPlot{pp}
%         %plot(meanThing(po),keptFractions(cc,po)','.-')
%         x = meanThing(po);
%         y = keptFractions(cc,po);
%         plot(x,y,'.')
%         hold on
%         p = polyfit(x,y,6);
%         plot(x,polyval(p,x),'-')
%     end
%     plot([0 1],[0 1],'k-')
%     plot([0 1 1 0 0],[0 0 1 1 0],'k-')
%     axis image
%     title(tl)
%     xlabel('mean fraction kept')
%     ylabel('fraction kept')
% end
% 
% 
% 
% %%
% 
% 
% figure
% %plot([0 1],[0 1],'k',[1 1],[0 1],'k')
% plot([0 1],[0 1],'-','color',[1 1 1]*.5)
% hold on
% plot(meanKeptFractionFalse,meanKeptFractionTrue,'.k','markersize',15)
% xlabel('false transients kept')
% ylabel('true transients kept')
% axis image
% 
% for ss=1:nParamSets
%     
%     %text(meanKeptFractionFalse(ss),meanKeptFractionTrue(ss),sprintf('\n%d',ss),'verticalalign','baseline','horizontalAlign','center','color','k','fontsize',18)
%     
% end
% 
% 
% 
% 
% % show what each parameter set was
% 
% fprintf('\nCOMMON PARAMETERS:\n\n')
% disp(rmfield(results.paramSets{1},results.nonSingletonParamNames))
% 
% 
% fprintf('PARAMETERS THAT VARIED:\n\n')
% for ss = 1:nParamSets
%     fprintf(' set %d\n',ss)
%     disp(rmfield(results.paramSets{ss},results.singletonParamNames))
% end
% 
% 
% 
% 
% ptToCol = @(x,y)hsv2rgb(cat(3,(1-(x+y)/2)*2/3+1/3,((y-x+1)/2).^3*.7+.3,ones(size(x))));
% 
% 
% 
% if 0||0
% %%    
%     figure(11);clf
%     
%     x=0:.05:1; y=0:.05:1;
%     [X,Y] = meshgrid(x,y);
%     
%     H = (X+Y)/2;
%     %H = H-min(H(:));
%     %H = H/max(H(:));
%     %H = H/2;
%     
%     %S = (Y-X + 1)/2;
%     S = max(1-X,Y);
%     %S = S+1;
%     %S = S/2;
%     
%     H = mod((1-H)*1/2+5/12+eps,1);
%     %H = mod((1-H)*2/3+1/3+eps,1);
%     
%     S = S.^3*.8+.2-eps;
%     
%     C = hsv2rgb(cat(3,H,S,ones(size(H))));
%     
%     subplot(1,4,1)
%     imagesc(H);axis image; axis xy;colormap gray;colorbar
%     
%     subplot(1,4,2)
%     imagesc(S);axis image; axis xy;colormap gray;colorbar
%     
%     subplot(1,4,3)
%     imagesc(C,'xdata',[0 1],'ydata',[0 1]);axis image; axis xy;colormap gray;colorbar
%     
%     %ptToCol = @(x,y)hsv2rgb(cat(3,(1-(x+y)/2)*2/3+1/3,((y-x+1)/2).^3*.8+.2,ones(size(x))));
%     ptToCol = @(x,y)hsv2rgb(cat(3,(1-(x+y)/2)*2/3+1/3,max(1-x,y).^3*.8+.2,ones(size(x))));
%     
%     subplot(1,4,4)
%     
%     for xx=x
%         for yy=y
%             plot(xx,yy,'.','color',ptToCol(xx,yy),'markersize',30)
%             hold on
%         end
%     end
%     xlim([0 1]);ylim([0 1])
%     axis image
%     colorbar
% end
% 
% 



% PLOT ALL TRANSIENTS

function pushbuttonPlotTransients_Callback(hObject, eventdata, handles)
plotTransients(handles,str2double(get(handles.editPlotTransients,'string')))

function editPlotTransients_Callback(hObject, eventdata, handles)
plotTransients(handles,str2double(get(handles.editPlotTransients,'string')))



% --- Executes during object creation, after setting all properties.
function editPlotTransients_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editPlotTransients (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function plotTransients(handles,whichParamSet)

tF = find(handles.results.pickedCellTransients(handles.whichTransients,3)==0);
tT = find(handles.results.pickedCellTransients(handles.whichTransients,3)==1);
nFalse = length(tF);
nTrue = length(tT);

nX = 2;
nY = max(nFalse,nTrue);

figure(1);clf
[~,hMat] = makeSubplots(gcf,nX,nY,.1,.2,[0 0 1 1]);
set(gcf,'color','w')

for tt=1:length(handles.whichTransients)

    if ismember(tt,tF)
        h = hMat(find(tF==tt),2);
    else
        h = hMat(find(tT==tt),1);
    end
    
    
    xL = handles.results.tc(handles.whichTransients(tt)).tcDefault.tc(:,1); % default time course
    %xL = handles.results.tc(handles.whichTransients(tt)).tcSeudo(whichParamSet).tcLSQ; % within-SEUDO time course
    xS = handles.results.tc(handles.whichTransients(tt)).tcSeudo(whichParamSet).tc;
    
    switch 1
        case 1 % patch
            patch([1 1:length(xL) length(xL)],[0; xL; 0],seudo.pickColor('LSQ'),...
                'faceAlpha',.3,'parent',h,'edgecolor',seudo.pickColor('LSQ'),'linewidth',2)
            hold(h,'on')
            patch([1 1:length(xS) length(xS)],[0; xS; 0],seudo.pickColor('SEcell'),...
                'faceAlpha',.3,'parent',h,'edgecolor',seudo.pickColor('SEcell'),'linewidth',2)
        case 2 % lines
            plot(h,xL,'color',seudo.pickColor('LSQ'),'linewidth',3);
            hold(h,'on')
            plot(h,xS,'color',seudo.pickColor('SEcell'));
            
            
    end
    axis(h,'tight')
    box(h,'off')
    yt = get(h,'ytick');
    set(h,'xtick',[],'ytick',[0 yt(end)],'xcolor','w','tickdir','out')
    ylim(h,[0 100])
    
    hL = plot(h,xlim(h),[0 0],'-','color',[1 1 1]*.6);
    uistack(hL,'bottom')
end






% USER INTERFACE ELEMENTS






% COMBINE TRANSIENTS

function popupCombineTransients_Callback(hObject, eventdata, handles)

updateROCplot(handles)

function popupCombineTransients_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% MEASURE OF EACH TRANSIENT

function popupTransientMeasure_Callback(hObject, eventdata, handles)

updateROCplot(handles)


function popupTransientMeasure_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
