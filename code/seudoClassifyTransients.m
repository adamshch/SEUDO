function varargout = seudoClassifyTransients(varargin)
% SEUDOCLASSIFYTRANSIENTS MATLAB code for seudoClassifyTransients.fig
%      SEUDOCLASSIFYTRANSIENTS, by itself, creates a new SEUDOCLASSIFYTRANSIENTS or raises the existing
%      singleton*.
%
%      H = SEUDOCLASSIFYTRANSIENTS returns the handle to a new SEUDOCLASSIFYTRANSIENTS or the handle to
%      the existing singleton*.
%
%      SEUDOCLASSIFYTRANSIENTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEUDOCLASSIFYTRANSIENTS.M with the given input arguments.
%
%      SEUDOCLASSIFYTRANSIENTS('Property','Value',...) creates a new SEUDOCLASSIFYTRANSIENTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before seudoClassifyTransients_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to seudoClassifyTransients_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help seudoClassifyTransients

% Last Modified by GUIDE v2.5 12-Dec-2018 17:50:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @seudoClassifyTransients_OpeningFcn, ...
                   'gui_OutputFcn',  @seudoClassifyTransients_OutputFcn, ...
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










% --- Executes just before seudoClassifyTransients is made visible.
function seudoClassifyTransients_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to seudoClassifyTransients (see VARARGIN)



% store inputs in handles struct
handles.se = varargin{1};
handles.tcStructSpec = varargin{2};


% set title
set(hObject,'name',sprintf('%s: transient classification for %s(%d), [%d x %d] pixels, %d frames, %d cells',...
    handles.se.name,handles.tcStructSpec{1},handles.tcStructSpec{2},handles.se.movY,handles.se.movX,handles.se.movF,size(handles.se.profiles,3)))


% set tool tips

set(handles.popupExhaustiveBehavior,'tooltipstring',...
    ['<html>class of transients after first click:<br />' ...
    '<b>&nbsp;&nbsp;T-F - </b>true up to clicked transient, then false<br />' ...
    '<b>&nbsp;&nbsp;T-U - </b>true up to clicked transient, then unclassified</html>']);

set(handles.popupColormap,'tooltipstring',...
    ['<html>colormap for displaying transient profiles:<br />' ...
    '<b>&nbsp;&nbsp;gray, all pix - </b>black-white colormap, normalized to darkest and brightest pixels<br />' ...
    '<b>&nbsp;&nbsp;gray, profile pix - </b>black-white colormap, normalized to darkest and brightest pixels in the profile<br />' ...
    '<b>&nbsp;&nbsp;magenta/green, all pix - </b>magenta for negative pixels, green for positive</html>']);

set(handles.popupSortOrder,'tooltipstring',...
    ['<html>how to sort transient profiles:<br />' ...
    '<b>&nbsp;&nbsp;correlation - </b>correlation with source profile, high to low<br />' ...
    '<b>&nbsp;&nbsp;time - </b>in the order they occured<br />' ...
    '<b>&nbsp;&nbsp;amplitude - </b>peak amplitude,large to small<br />' ...
    '<b>&nbsp;&nbsp;classification - </b>unclassified, true, mixed, false<br />' ...
    '<b>&nbsp;&nbsp;clustering - </b>rough attempt to cluster, sorted by cluster number<br />' ...
    '<b>&nbsp;&nbsp;contamination severity - </b>see manuscript for definition, low to high</html>']);

set(handles.popupTransientTitle,'tooltipstring',...
    ['<html>title to display for each transient profile:<br />' ...
    '<b>&nbsp;&nbsp;transient ID</b><br />' ...
    '<b>&nbsp;&nbsp;correlation - </b>correlation with source profile<br />' ...
    '<b>&nbsp;&nbsp;best match source - </b>source ID of source with highest correlation<br />' ...
    '<b>&nbsp;&nbsp;peak frame - </b>frame with peak amplitude<br />' ...
    '<b>&nbsp;&nbsp;contamination severity - </b>see manuscript for definition<br />' ...
    '<b>&nbsp;&nbsp;none</b></html>']);





% note # cells
handles.nCells = size(handles.se.profiles,3);




% hide axesTransientsArea, since it is only used to set bounds for subplots
set(handles.axesTransientsArea,'visible','off')

% initialize set of transient subplots to be empty
%handles.axesTransSubplots = [];



% compute residual ratios
handles.transStats = handles.se.autoClassifyTransients('default','saveResults',0);


% set slider range for cell
set(handles.sliderCellNumber,'Max',handles.nCells)
set(handles.sliderCellNumber,'Min',1)
if handles.nCells > 1
    set(handles.sliderCellNumber,'SliderStep',[1 10]/(handles.nCells-1))
else
    
end


% set slider range for unclassified transients slider
% NOTE: this just moves one unit left or right, it's not like a real slider
set(handles.sliderUncl,'Min',-1,'Max',1,'SliderStep',[1 1],'Value',0)



% identify which profiles are at which locations
pMap = convn(handles.se.profiles>0,ones(3,3,1),'same');
for yy = 1:handles.se.movY
    for xx = 1:handles.se.movX
        handles.pLookup{yy}{xx} = sort(find(pMap(yy,xx,:)));
    end
end
handles.lastClick = nan;
handles.lastClick = 0;


% initialize display parameters

handles.thisCell = 1;
handles.lastCell = 1;
set(handles.sliderCellNumber,'Value',1)
set(handles.editTextCellNumber,'String','1')
set(handles.textNewPlots,'visible','off')

handles.oldBlurRadius = 0;
handles.convKernel = [];

handles.nTransY = 3;
handles.nTransX = 5;
set(handles.editNY,'string',num2str(3));
set(handles.editNX,'string',num2str(5));

handles = makeNewSubplots(hObject,handles);

%colormap(hObject,'gray')



set(handles.textTrue,'ForegroundColor',getPlotColor(seudo.valTrue))
set(handles.textFalse,'ForegroundColor',getPlotColor(seudo.valFalse))
set(handles.textMixed,'ForegroundColor',getPlotColor(seudo.valMix))
set(handles.textUnclassified,'ForegroundColor',getPlotColor(seudo.valUnc))
set(handles.textTrueCount,'ForegroundColor',getPlotColor(seudo.valTrue))
set(handles.textFalseCount,'ForegroundColor',getPlotColor(seudo.valFalse))
set(handles.textMixCount,'ForegroundColor',getPlotColor(seudo.valMix))
set(handles.textUncCount,'ForegroundColor',getPlotColor(seudo.valUnc))


set(handles.textNonArt,'ForegroundColor',[1 1 1]*0)
set(handles.textNonArtCount,'ForegroundColor',[1 1 1]*0)
set(handles.textArt,'ForegroundColor',seudo.pickColor('artifact'))
set(handles.textArtCount,'ForegroundColor',seudo.pickColor('artifact'))


handles.overlayProfile = 0;



% catch keyboard inputs
set(hObject,'WindowKeyPressFcn',@(o,e)seudoClassifyTransients_WindowKeyPressFcn(o,e));




% function for switching classificaiton
set(handles.uibuttongroupClass,'SelectionChangedFcn',@(o,e)uibuttongroupClassCallback(o,e))



% Choose default command line output for seudoClassifyTransients
handles.output = hObject;


% Update handles structure
guidata(hObject, handles);


% Set default focus UI element
handles.defaultFocus = handles.sliderCellNumber;
%handles.defaultFocus = handles.sliderUncl;


handles = setNewCell(handles);



updatePlot(hObject, handles)



% --- Outputs from this function are returned to the command line.
function varargout = seudoClassifyTransients_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;





% SWITCH CLASSIFICATION TYPE

function uibuttongroupClassCallback(hObject,~)

handles = guidata(hObject);

% update instruction text
updateInstructionsText(handles)

% update to show change
updatePlot(hObject,handles)



function updateInstructionsText(handles)

% identify # transients
ti = handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo(handles.thisCell);
nTrans = size(ti.times,1);

% reset the instructions text
if nTrans == 0
    set(handles.textInstructions,'string','This cell has no transients')
elseif all(isnan(ti.classification)) && get(handles.radiobuttonExhaustive,'value')
    % exhaustive
    switch get(handles.popupExhaustiveBehavior,'value')
        case 1 % T -> F
            set(handles.textInstructions,'string','Click on the last true transient (others will be marked false)')
        case 2 % T -> U
            set(handles.textInstructions,'string','Click on the last true transient (others will remain unclassified)')
    end
else
    % selective
    set(handles.textInstructions,'string','Click on a transient to cycle classification')
end


switch get(handles.radiobuttonExhaustive,'value')
    case 0
        set(handles.popupExhaustiveBehavior,'enable','off')
    case 1
        set(handles.popupExhaustiveBehavior,'enable','on')
end



function popupExhaustiveBehavior_Callback(hObject, eventdata, handles)

% just update instruction text
updateInstructionsText(handles)

function popupExhaustiveBehavior_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% SHOW CLASSIFICATION

function pushbuttonShowClass_Callback(hObject, eventdata, handles)

% get true and false classified transients
ti = handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo;
transSpec = [];
for cc = 1:length(ti)
    ft = find(ti(cc).classification == seudo.valFalse);
    tt = find(ti(cc).classification == seudo.valTrue);
    if ~isempty(ft)
        transSpec = [transSpec; cc*ones(length(ft),1) ft zeros(length(ft),1) ];
    end
    if ~isempty(tt)
        transSpec = [transSpec; cc*ones(length(tt),1) tt ones(length(tt),1) ];
    end
end

% print result to command window
if isempty(transSpec)
    fprintf('no transients were classified as true or false\n')
else
    str = '[';
    for tt = 1:size(transSpec,1)
        if transSpec(tt,3) == 1, transStr = 'true'; else transStr = 'false'; end
        str = [str sprintf('%d %d %s; ',transSpec(tt,1:2),transStr)];
    end
    str = [str(1:end-2) ']'];
    
    fprintf('transients classified as true or false in %s:\n%s\n',handles.se.name,str)
end


% return focus to slider
uicontrol(handles.defaultFocus)





% SAVE TO VARIABLE

function pushbuttonClassToVar_Callback(hObject, eventdata, handles)

% get true and false classified transients
ti = handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo;
transSpec = [];
for cc = 1:length(ti)
    ft = find(ti(cc).classification == seudo.valFalse);
    tt = find(ti(cc).classification == seudo.valTrue);
    if ~isempty(ft)
        transSpec = [transSpec; cc*ones(length(ft),1) ft zeros(length(ft),1) ];
    end
    if ~isempty(tt)
        transSpec = [transSpec; cc*ones(length(tt),1) tt ones(length(tt),1) ];
    end
end

% print result to command window
if isempty(transSpec)
    fprintf('no transients were classified as true or false\n')
else
    % save result to variable in base workspace
    varName = ['transients_' datestr(now,'YYYY_MM_DD__hh_mm_ss')];
    assignin('base',varName,transSpec)
    fprintf('saved to variable ''%s''\n\n',varName)
end


% return focus to slider
uicontrol(handles.defaultFocus)








% KEYBOARD SHORTCUTS

function seudoClassifyTransients_WindowKeyPressFcn(hObject, eventdata)

% skip if modifier key
if ~isempty(eventdata.Modifier), return, end

% skip if entering text
G = gco;
if ~isempty(isprop(G,'Style')) && isprop(G,'Style') && isequal(G.Style,'edit'), return, end

% get handles object
handles = guidata(hObject);

% identify which key was pressed
switch eventdata.Key
    
    case {'period','comma'}
        % move vertical slider
        
        % comma moves up, period moves down
        if strcmp(eventdata.Key,'period'), mul = -1; else mul = 1; end
        
        % get current value and slider properties
        currVal = get(handles.sliderTransients,'value');
        stepSize = get(handles.sliderTransients,'SliderStep'); 
        maxVal = get(handles.sliderTransients,'max');
        minVal = get(handles.sliderTransients,'min');
        
        % identify what the new value would be
        newVal = currVal + mul * (maxVal - minVal) * stepSize(1);
        
        % if it's valid, execute
        if maxVal > minVal && newVal >= minVal && newVal <= maxVal
            set(handles.sliderTransients,'value',newVal)
        end

        % update to show change
        updatePlot(hObject,handles)
        
    case {'o','v'}
        % cycle overlay
        cycleProfileOverlay(hObject,[])
        
    case {'a'}
        % toggle artifact checkbox
        set(handles.checkboxIsArtifact,'value',1-get(handles.checkboxIsArtifact,'value'));
        checkboxIsArtifact_Callback(handles.checkboxIsArtifact, [], handles)
        
    case {'c'}
        % cycle colormap
        set(handles.popupColormap,'value',mod(get(handles.popupColormap,'value'),3)+1)
        popupColormap_Callback(handles.popupColormap, [], handles)
        
        
    case {'0','1','2','3','4','5','6','7','8','9'}
        set(handles.editBlur,'string',eventdata.Key)
        editBlur_Callback(handles.editBlur, [], handles)
        
    case {'backquote'}
        set(handles.editBlur,'string','0')
        editBlur_Callback(handles.editBlur, [], handles)
end

 
 
 








% CHANGE CELL NUMBER BY MOVING SLIDER

function sliderCellNumber_Callback(hObject, eventdata, handles) %#ok<*INUSL,*DEFNU>

% note which cell this is (it will be the next "last cell")
handles.lastCell = handles.thisCell;

% get slider value
handles.thisCell = round(get(hObject,'Value'));

% set editable text box
set(handles.editTextCellNumber,'String',num2str(handles.thisCell))

% change cell
handles = setNewCell(handles);

% update plot
updatePlot(hObject,handles)



function sliderCellNumber_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end





% CHANGE CELL NUMBER BY EDITING TEXT

function editTextCellNumber_Callback(hObject, eventdata, handles)

% note which cell this is (it will be the next "last cell")
handles.lastCell = handles.thisCell;

% get value
handles.thisCell = round(str2num(get(hObject,'String')));

% ensure in range
handles.thisCell = min(max(1,handles.thisCell),size(handles.se.profiles,3));

% update slider
set(handles.sliderCellNumber,'Value',handles.thisCell)

% ensure editable text box is clear
set(hObject,'String',num2str(handles.thisCell))

% change cell
handles = setNewCell(handles);

% update plot
updatePlot(hObject,handles)


function editTextCellNumber_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% SLIDER TO CHANGE BETWEEN SCREEN OF TRANSIENTS


function sliderTransients_Callback(hObject, eventdata, handles)

set(handles.sliderTransients,'value',round(get(handles.sliderTransients,'value')))
updatePlot(hObject,handles)


function sliderTransients_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end







% FIELDS CONTROLING # OF SUBPLOTS IN X AND Y

function editNX_Callback(hObject, eventdata, handles)

% disable interface while plotting
objToDisable=findobj(get(hObject,'parent'),'Enable','on');
set(objToDisable,'Enable','off');

% get value
val = round(str2double(get(hObject,'String')));

% ensure in range
handles.nTransX = max(min(20,val),2);

% ensure editable text box is clear
set(hObject,'String',num2str(handles.nTransX))

% update subplots
handles = makeNewSubplots(get(hObject,'parent'),handles);

% change cell
handles = setNewCell(handles);

% update plot
updatePlot(hObject,handles)

% turn interface back on
set(objToDisable,'Enable','on');



function editNX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function editNY_Callback(hObject, eventdata, handles)

% disable interface while plotting
objToDisable=findobj(get(hObject,'parent'),'Enable','on');
set(objToDisable,'Enable','off');

% get value
val = round(str2double(get(hObject,'String')));

% ensure in range
handles.nTransY = max(min(20,val),2);

% ensure editable text box is clear
set(hObject,'String',num2str(handles.nTransY))

% update subplots
handles = makeNewSubplots(get(hObject,'parent'),handles);

% change cell
handles = setNewCell(handles);

% update plot
updatePlot(hObject,handles)

% turn interface back on
set(objToDisable,'Enable','on');




function editNY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% SET SORT ORDER



function popupSortOrder_Callback(hObject, eventdata, handles)

% update plot
updatePlot(hObject,handles)


function popupSortOrder_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% CLEAR CLASSIFICATION


function pushbuttonClearClassificaiton_Callback(hObject, eventdata, handles)

% reset classification
handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo(handles.thisCell).classification(:) = nan;

% change cell
handles = setNewCell(handles);

% update plot
updatePlot(hObject,handles)






% CHANGE COLORMAP

function popupColormap_Callback(hObject, eventdata, handles)
updatePlot(hObject,handles)

function popupColormap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% SWITCH TO LAST CELL

function pushbuttonBackCell_Callback(hObject, eventdata, handles)

% note which cell this is (it will be the next last cell)
switchingFrom = handles.thisCell;

% set new cell
handles.thisCell = handles.lastCell;

% note the last cell
handles.lastCell = switchingFrom;

% set editable text box
set(handles.editTextCellNumber,'String',num2str(handles.thisCell))

% change cell
handles = setNewCell(handles);

% update plot
updatePlot(hObject,handles)






% SAVE CLASSIFICATION TO DISK

function pushbuttonSaveClassification_Callback(hObject, eventdata, handles)
% save time courses and classification to a file

% get save destination
[FileName,PathName] = uiputfile('','Choose save location and name',[handles.se.name '_class.mat']);

% abort if cancelled
if ~ischar(FileName), return, end

saveName = [PathName FileName];

% save as struct
fprintf('saving classification of %s to %s...',handles.se.name,saveName)
tcStruct = handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2});
save(saveName,'-struct','tcStruct','-v7.3')

% help text
tcStructType = handles.tcStructSpec{1};
fprintf('done\nto load, create seudo object with parameter pair:\n    seudo(...,''%s'',''%s'')\n',...
    tcStructType,saveName)










% MAKE NEW SUBPLOTS

function handles = makeNewSubplots(figHandle,handles)

% make subplots
%posFig = get(figHandle,'position');
posFig = [1 1 1 1];

% create one big subplot, just to be sure previous transient plots
subplot('position',get(handles.axesTransientsArea,'position')./posFig([3 4 3 4]),...
    'parent',get(handles.axesProfile,'parent'))


set(handles.textNewPlots,'visible','on')
drawnow

handles.axesTransSubplots = makeSubplots(get(handles.axesTransientsArea,'parent'),...
    handles.nTransX,handles.nTransY,.1,.2,get(handles.axesTransientsArea,'position')./posFig([3 4 3 4]));

handles.axesTransSubplotTitles = nan(size(handles.axesTransSubplots));
for hh = 1:length(handles.axesTransSubplots);
    %h = handles.axesTransSubplots(:)'
    h = handles.axesTransSubplots(hh);
    axis(h,'image')
    axis(h,'ij')
    axis(h,'off')
    hold(h,'on')
    handles.axesTransSubplotTitles(hh) = title(h,'');
end

set(handles.textNewPlots,'visible','off')







% CHANGE TRANSIENT PLOT TITLES

function popupTransientTitle_Callback(hObject, eventdata, handles)

updatePlot(hObject,handles)


function popupTransientTitle_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end









% CHANGE BLUR RADIUS

function editBlur_Callback(hObject, eventdata, handles)

% validate/constrain input
newInput = str2num(get(handles.editBlur,'string'));
if newInput < 0
    newInput = 0;
elseif length(newInput) ~= 1 || ~isfinite(newInput) || ~isreal(newInput)
    newInput = handles.oldBlurRadius;
end
set(handles.editBlur,'string',num2str(newInput));
handles.oldBlurRadius = newInput;

% set convolution kernel
if newInput > 0
    handles.convKernel = fspecial('gauss',ceil([1 1]*(newInput+1)*2),newInput);
else
    handles.convKernel = [];
end

updatePlot(hObject,handles)

function editBlur_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end










% SET TRANSIENT PLOT COLORS

function [plotCol,plotSym] = getPlotColor(cl)
plotSym = '.';

if     isnan(cl),               plotCol = [1 1 1]*0.5;
elseif cl == seudo.valTrue,     plotCol = seudo.pickColor('true');
elseif cl == seudo.valFalse,    plotCol = seudo.pickColor('false');
else                            plotCol = seudo.pickColor('mixed'); plotSym = 'o';
end








% MOVE BETWEEN UNCLASSIFIED SOURCES

function sliderUncl_Callback(hObject, eventdata, handles)

% identify unclassified sources (i.e. at least 1 transient, all transients unclassified, not an artifact)
isUncl = cellfun(@(x)all(isnan(x))&~isempty(x),{handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo.classification});
isArt = [handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo.isArtifact] ~= 0;
isUncl = find(isUncl & ~isArt);

% get ID of next unclassified source (based on whether slider was moved up or down)
switch get(hObject,'value')
    case 1
        nextCell = find(isUncl > handles.thisCell,1,'first');
    case -1
        nextCell = find(isUncl < handles.thisCell,1,'last');
    otherwise
        disp(get(hObject,'value'))
        error('BUG: weird value for unclassified cell slider')
end

% reset value of this slider
set(hObject,'value',0)

% update plot
if ~isempty(nextCell)
    % store value of previous cell
    handles.lastCell = handles.thisCell;
    % set new cell
    handles.thisCell = isUncl(nextCell);
    % update data and plots
    handles = setNewCell(handles);
    updatePlot(hObject,handles)
end

% return focus to main slider
uicontrol(handles.defaultFocus)


function sliderUncl_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end







% IS ARTIFACT CHECK BOX

function checkboxIsArtifact_Callback(hObject, eventdata, handles)

% store value in seudo struct
handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo(handles.thisCell).isArtifact = ...
    get(hObject,'value');

% update plot
updatePlot(hObject,handles)

% return focus to slider
uicontrol(handles.defaultFocus)








% HANDLE CLICKING ON TRANSIENT SHAPE TO CHANGE CLASSIFICATION

function updateClassification(axesObj,event,handles,whichCell,trueTransNumber,allTransToThisOne)

% get handles to axes and figure
%axesObj = get(imageObj,'parent');
hFigure = get(axesObj,'parent');


switch get(hFigure,'selectiontype')
    case 'alt'
        % generate GUI to examine transient
        
        % get movie
        transInfo = handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo(whichCell);
        whichFrames = transInfo.times(trueTransNumber,1):transInfo.times(trueTransNumber,2);
        theMov = handles.se.getFrames(whichFrames,...
            transInfo.window(1):transInfo.window(2),transInfo.window(3):transInfo.window(4));
        
        % set name
        guiName = sprintf('%s: source %d, transient %d (frames %d to %d)',...
            handles.se.name,whichCell,trueTransNumber,transInfo.times(trueTransNumber,:));
        
        % launch GUI
        theProfile = handles.se.profiles(transInfo.window(1):transInfo.window(2),transInfo.window(3):transInfo.window(4),whichCell);
        seudoViewTransient(theMov,whichFrames,theProfile,guiName)
        return
end



isExhaustive = get(handles.radiobuttonExhaustive,'value');

% if no classification was performed yet...

if isExhaustive && handles.noClicksThisCell && ...
        all(isnan(handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo(whichCell).classification))
    
    % first click makes all good from first trans to this one, all subsequent bad
    
    
    
    switch get(handles.popupExhaustiveBehavior,'value')
        case 1
            transClassFirst = seudo.valTrue;
            transClassLast = seudo.valFalse;
        case 2
            transClassFirst = seudo.valTrue;
            transClassLast = seudo.valUnc;
    end

    % set classification
    handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo(whichCell).classification(:) = transClassLast;
    handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo(whichCell).classification(allTransToThisOne) = transClassFirst;
    
    % change color of time courses and scatter plot points
    for tt = allTransToThisOne
        set(handles.trPlots(tt),'color',getPlotColor(transClassFirst))
        set(handles.trScatterPlots(tt),'color',getPlotColor(transClassFirst))
    end
    for tt = setdiff(1:length(handles.trPlots),allTransToThisOne)
        set(handles.trPlots(tt),'color',getPlotColor(transClassLast))
        set(handles.trScatterPlots(tt),'color',getPlotColor(transClassLast))
    end
    
    handles.noClicksThisCell = false;
    
    % update plot
    updatePlot(hFigure,handles)
    
else
    % otherwise, cycle through classifications
    
    
    % identify new classification.
    % udpate rule:  unclassified -> false -> true -> mixed -> unclassified
    thisTransClass = handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo(whichCell).classification(trueTransNumber);
    if      isnan(thisTransClass),            thisTransClass = seudo.valFalse;
    elseif  thisTransClass == seudo.valFalse, thisTransClass = seudo.valTrue;
    elseif  thisTransClass == seudo.valTrue,  thisTransClass = seudo.valMix;
    else                                      thisTransClass = seudo.valUnc;
    end
        
    % store in seudo object
    handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo(whichCell).classification(trueTransNumber) = thisTransClass;
    
    % set color of time course
    set(handles.trPlots(trueTransNumber),'color',getPlotColor(thisTransClass))
    
    % set color of scatter plot point
    set(handles.trScatterPlots(trueTransNumber),'color',getPlotColor(thisTransClass))
    
    % set color of transient image frame
    set(findall(get(axesObj,'child'),'type','line','linewidth',borderLineWidth),...
        'color',getPlotColor(thisTransClass));
    
end

% update instructions
updateInstructionsText(handles)


% note # true, false transients
updateCellCount(handles)


% return focus to slider
uicontrol(handles.defaultFocus)




% CYCLE PROFILE OVERLAY

function cycleProfileOverlay(o,~)

% increment one step
handles = guidata(o);
handles.overlayProfile = mod(handles.overlayProfile + 1,3);

% hide if needed
if handles.overlayProfile ~= 2
    set(handles.contourPlotHandle,'visible','off')
else
    set(handles.contourPlotHandle,'visible','on')
end

% store value and plot
guidata(o,handles)
updatePlot(o,handles)









% UPDATE CELL NUMBER


function handles = setNewCell(handles)
% perform steps required when a new cell is selected


handles.noClicksThisCell = true;

% get transientInfo struct for this cell
ti = handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo(handles.thisCell);


% set up slider for moving between sets of transients

% get subplot size
nTransY = handles.nTransY;
nTransX = handles.nTransX;

% identify # transients
nTrans = size(ti.times,1);

% identify # of sets of subplots
handles.nScreens = ceil((nTrans-nTransX*(nTransY-1))/nTransX);

% set slider range
if handles.nScreens > 1
    set(handles.sliderTransients,'visible','on')
    set(handles.sliderTransients,'Max',handles.nScreens)
    set(handles.sliderTransients,'Min',1)
    set(handles.sliderTransients,'SliderStep',[1 handles.nTransY]/(handles.nScreens-1))
    set(handles.textUpDown,'visible','on')
else
    % make slider disappear
    set(handles.sliderTransients,'Max',1)
    set(handles.sliderTransients,'Min',1)
    set(handles.sliderTransients,'visible','off')
    set(handles.textUpDown,'visible','off')
    drawnow
end

% start with topmost screen
set(handles.sliderTransients,'Value',handles.nScreens)



% reset the instructions text
updateInstructionsText(handles)





% get correlation between cell profile and each transient
prof = handles.se.profiles(ti.window(1):ti.window(2),ti.window(3):ti.window(4),handles.thisCell);
if ~isempty(ti.shapes)
    %handles.corrVals = correlationVectorMatrix(prof(:),reshape(ti.shapes,[],size(ti.shapes,3)));
    handles.corrVals = handles.transStats(handles.thisCell).corrs;
else
    handles.corrVals = [];
end

% get amplitude of each transient
handles.tAmp = nan(size(ti.times,1),1);
for tt = 1:size(ti.times,1)
    handles.tAmp(tt) = max(handles.se.(handles.tcStructSpec{1}).tc(ti.times(tt,1):ti.times(tt,2),handles.thisCell));
end


% get order for plotting via clustering
if size(ti.shapes,3) > 1
    Z = linkage([reshape(ti.shapes,[],size(ti.shapes,3))]');
    handles.transientsClusterOrder = optimalleaforder(Z,pdist([ reshape(ti.shapes,[],size(ti.shapes,3))]'));
else
    handles.transientsClusterOrder = 1;
end
%C = corr([prof(:) reshape(ti.shapes,[],size(ti.shapes,3))]);
%Z = linkage(C);
%plotOrder = optimalleaforder(Z,C);

%plotOrder = setdiff(plotOrder,length(plotOrder));
%plotOrder = plotOrder - 1;

%[~,~,plotOrder] = dendrogram(Z);

%O = cluster(Z,'maxclust',2);
%[~,plotOrder] = sort(O);
%plotOrder = Z(:,1);
%plotOrder = [1 29 10 17 18 9 16 26 4 14 23 25 24 22 28 19 15 2 3 6 8 7 20 11 13 5 12 21 27 30];
        
        


% plot profile
cla(handles.axesProfile)
cutWin = ti.window;
handles.cellProfile = handles.se.profiles(cutWin(1):cutWin(2),cutWin(3):cutWin(4),handles.thisCell);
imagesc(handles.cellProfile,'parent',handles.axesProfile,'ButtonDownFcn',@(o,e)cycleProfileOverlay(o,e))
axis(handles.axesProfile,'image')
axis(handles.axesProfile,'off')
colormap(handles.axesProfile,'gray')

% add contour
handles.thisContour = contourc(double(handles.cellProfile>0),[1 1]);%,'color','r','linewidth',3);
hold(handles.axesProfile,'on')
handles.contourPlotHandle = plot(...
    handles.axesProfile,handles.thisContour(1,2:end),handles.thisContour(2,2:end),'-',...
    'color',cellHighlight,'linewidth',2);

% hide if needed
if handles.overlayProfile ~= 2
    set(handles.contourPlotHandle,'visible','off')
end


% plot time course
tc = handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).tc(:,handles.thisCell);
cla(handles.axesTimeCourse)
plot(handles.axesTimeCourse,tc,'color',[1 1 1]*.9)
hold(handles.axesTimeCourse,'on')
set(handles.axesTimeCourse,'box','off')

% plot each transient, colored according to its classification
tTimes = ti.times;
handles.trPlots = nan(size(tTimes,1),1);
for tt=1:size(tTimes,1)
    handles.trPlots(tt) = plot(handles.axesTimeCourse,tTimes(tt,1):tTimes(tt,2),tc(tTimes(tt,1):tTimes(tt,2)),'-',...
        'color',getPlotColor(ti.classification(tt)));
end
set(handles.axesTimeCourse,'tickdir','out')
axis(handles.axesTimeCourse,'tight')



% highlight profile in image of all profiles
thisProfile = handles.se.profiles(:,:,handles.thisCell);
thisProfile = thisProfile / max(thisProfile(:));
profMax = max(bsxfun(@rdivide,handles.se.profiles,max(max(handles.se.profiles))),[],3);
profMax = repmat(profMax / max(profMax(:)),[1 1 3]);
thisProfile = cat(3,thisProfile*0,thisProfile*1,thisProfile*.5);
imToPlot = 1-(1-profMax*.7).*(1-thisProfile);
%imToPlot = imToPlot .^ .5;

% plot image
cla(handles.axesAllProfiles)
imagesc(imToPlot,'parent',handles.axesAllProfiles,'hittest','off')
axis(handles.axesAllProfiles,'image')
set(handles.axesAllProfiles,'xtick',[],'ytick',[])
hold(handles.axesAllProfiles,'on')

% add line of window of cells that were cut out
cutWin = cutWin + [-.5 .5 -.5 .5];
plot(handles.axesAllProfiles,cutWin([3 4 4 3 3]),cutWin([1 1 2 2 1]),'-','linewidth',1,'color',cellHighlight,'hittest','off')
title(handles.axesAllProfiles,sprintf('all %d profiles',handles.se.nCells),'color','k','fontsize',12,'fontweight','normal','fontname','helvetica')


% set click callback function
set(handles.axesAllProfiles,'ButtonDownFcn',@(o,e)chooseCellFromClick(o,e))


% choose what to plot on Y axis
yVal = 2;

% scatter plot of transient parameters
handles.trScatterPlots = nan(size(tTimes,1),1);
cla(handles.axesCorrelation)
hold(handles.axesCorrelation,'on')
for tt=1:size(tTimes,1)
    X = handles.corrVals(tt);
    switch yVal
        case 1 % amplitude
            Y = handles.tAmp(tt);
        case 2 % residual ratio
            Y = handles.transStats(handles.thisCell).resRatios(tt);
    end
    [plotCol,plotSym] = getPlotColor(ti.classification(tt));
    handles.trScatterPlots(tt) = plot(handles.axesCorrelation,X,Y,plotSym,...
        'color',plotCol);
end
xlim(handles.axesCorrelation,[-1 1])
xlabel(handles.axesCorrelation,'correlation')
switch yVal
    case 1 % amplitude
        ylabel(handles.axesCorrelation,'amplitude')
    case 2 % residual ratio
        ylabel(handles.axesCorrelation,'contam. severity')
end





% CHOOSE NEW CELL BASED ON CLICKING THE IMAGE OF ALL CELLS

function chooseCellFromClick(axesObject,eventData)

figObject = get(axesObject,'parent');
handles = guidata(figObject);

% get click location
%xy = get(handles.axesAllProfiles,'CurrentPoint');
xy = round(eventData.IntersectionPoint(1:2));

% constrain
x = min(max(1,xy(1)),handles.se.movX);
y = min(max(1,xy(2)),handles.se.movY);

% identify which cells are at this location
cellList = handles.pLookup{y}{x};

% if no cell here, do nothing
if isempty(cellList), return, end

% identify how many times this location has been clicked
if isequal(handles.lastClick,[y x])
    % if the same as last click, increment click count
    handles.lastClickCount = mod(handles.lastClickCount,length(cellList))+1;
else
    % otherwise, just note this click location and reset the count
    handles.lastClick = [y x];
    handles.lastClickCount = 1;
end

% if different, switch to selected cell
if handles.thisCell ~= cellList(handles.lastClickCount)
    
    % note previous cell
    handles.lastCell = handles.thisCell;
    
    % set new cell and plot
    handles.thisCell = cellList(handles.lastClickCount);
    handles = setNewCell(handles);
    updatePlot(figObject,handles)
    
    % return focus to main slider
    uicontrol(handles.defaultFocus)
    
end





% UPDATE CELL COUNT

function updateCellCount(handles)
% show # of true, false transients classified

isArtSources = [handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo.isArtifact];

% get counts
classVals = cat(1,handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo(isArtSources==0).classification);
trueCount = sum(classVals == seudo.valTrue);
falseCount = sum(classVals == seudo.valFalse);
uncCount = sum(isnan(classVals));
mixCount = sum(~isnan(classVals) & ~ismember(classVals,[seudo.valTrue seudo.valFalse]));

artCount = sum(isArtSources~=0);
nonArtCount = sum(isArtSources==0);

% update text fields
set(handles.textTrueCount,'string',sprintf('%d',trueCount))
set(handles.textFalseCount,'string',sprintf('%d',falseCount))
set(handles.textMixCount,'string',sprintf('%d',mixCount))
set(handles.textUncCount,'string',sprintf('%d',uncCount))

set(handles.textArtCount,'string',sprintf('%d',artCount))
set(handles.textNonArtCount,'string',sprintf('%d',nonArtCount))











% UPDATE PLOT


function updatePlot(hObject,handles)


% identify current cell
thisCell = handles.thisCell;


% update cell-specific plots
%handles = setNewCell(handles);
 


% UPDATE TRANSIENT SUBPLOTS
%
% create many small subplots, one for each transient
% label each by how it was classified

nTransY = handles.nTransY;
nTransX = handles.nTransX;






% clear previous plots



% plot transients


% get current set of active shapes
ti = handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo(thisCell);

% identify which transients to plot
nTrans = size(ti.times,1);
if handles.nScreens > 1
    firstTransToPlot = ((get(handles.sliderTransients,'Max')-get(handles.sliderTransients,'Value'))) * nTransX + 1;
else
    firstTransToPlot = 1;
end
lastTransToPlot = min(nTrans,firstTransToPlot+nTransY*nTransX-1);
transToPlot = firstTransToPlot:lastTransToPlot;








% determine plotting order
switch get(handles.popupSortOrder,'value')
    
    case 1 % correlation
        [~,plotOrder] = sort(handles.corrVals,'descend');
        
    case 2 % time
        plotOrder = 1:nTrans;
        
    case 3 % amplitude
        [~,plotOrder] = sort(handles.tAmp,'descend');
        
    case 4 % classification
        [~,plotOrder] = sortrows([ti.classification handles.tAmp]);
        plotOrder = flipud(plotOrder);
        
    case 5 % clustering
        plotOrder = handles.transientsClusterOrder;
        
    case 6 % residual ratio
        [~,plotOrder] = sort(handles.transStats(handles.thisCell).resRatios);
        
    otherwise
        error('bug in the code: cannot handle sort order argument')
end






% plot each transient
for tt=1:length(transToPlot)

    h = handles.axesTransSubplots(tt);
    
    set(h,'visible','on')
    cla(h)
    
    trueTransNumber = plotOrder(transToPlot(tt));
    
    % get transient info for this cell
    ti = handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo(handles.thisCell);
    
    % get image of this transient
    transIm = ti.shapes(:,:,trueTransNumber);
    %     transIm = [ti.shapes(:,:,trueTransNumber),ti.shapesMean(:,:,trueTransNumber)];

    % convolve
    if ~isempty(handles.convKernel)
        transIm = conv2(transIm,handles.convKernel,'same');
    end
    
    % overlay profile, if descired
    addContour = false;
    switch handles.overlayProfile
        case 0 % no overlay
            
        case 1 % colored ghost
            
            % get colored version of cell profile
            thisProfile = handles.cellProfile;
            thisProfile = thisProfile / max(thisProfile(:));
            %thisProfile = cat(3,zeros([size(thisProfile) 2]),thisProfile); % profile in blue channel
            thisProfile = cat(3,zeros(size(thisProfile)),thisProfile,zeros(size(thisProfile))); % profile in green channel
            
            % combine transient shape and cell profile
            transIm = transIm - min(transIm(:));
            transIm = repmat(transIm / max(transIm(:)),[1 1 3]);
            transIm = 1-(1-transIm).*(1-thisProfile);
            
        case 2 % profile contour
            
            % it gets plotted after the image
            addContour = true;
            
        otherwise
            disp(handles.overlayProfile)
            error('BUG IN CODE: strange value of handles.overlayProfile')
    end
    
    % plot transient profiles
    imagesc(transIm ,'parent',h,'hittest','off')
    
    % set click behavior
    set(h,'ButtonDownFcn',@(o,e)updateClassification(o,e,handles,handles.thisCell,trueTransNumber,plotOrder(1:(tt+firstTransToPlot-1))))
    
    % colormap
    switch get(handles.popupColormap,'value')
        case 1 % gray
            colormap(h,'gray')
            caxis(h,[min(transIm(:)) max(transIm(:))]);
            
        case 2 % gray
            colormap(h,'gray')
            prof = handles.se.profiles(ti.window(1):ti.window(2),ti.window(3):ti.window(4),handles.thisCell);
            pixInProf = prof>0;
            caxis(h,[min(transIm(pixInProf)) max(transIm(pixInProf))]);
            
        case 3 % color to show polarity
            colormap(h,seudo.colormapGreenMagenta([],.5))
            caxis(h,[-1 1]*max(abs(transIm(:))));
        otherwise
            error('# of expected colormaps does not match selection')
    end
    
    % delete existing frame overlay
    delete(findall(get(h,'child'),'type','line'))
    
    % choose overlay color
    if ti.isArtifact, col = seudo.pickColor('artifact');
    else col = getPlotColor(ti.classification(trueTransNumber)); end
    
    % add overlay
    xl = [.5 size(transIm,2)+.5];
    yl = [.5 size(transIm,1)+.5];
    r = plot(h,xl([1 2 2 1 1 2]),yl([1 1 2 2 1 1]),'-','linewidth',borderLineWidth,'color',col);
    uistack(r,'bottom')
    
    
    % add contour, if applicable
    if addContour
        plot(h,handles.thisContour(1,2:end),handles.thisContour(2,2:end),'-','linewidth',2,'color',cellHighlight,'hittest','off');
    end
    
    
    % remove axis labels
    %axis(h,'off') % takes a LOOOOONG time
    set(h,'xtick',[],'ytick',[]) % faster
    
    % add title
    col = 'k';
    switch get(handles.popupTransientTitle,'value')
        case 1 % transient ID
            str = num2str(trueTransNumber);
        case 2 % correlation
            str = sprintf('%0.3f',handles.corrVals(trueTransNumber));
        case 3 % best match cell
            win = ti.window;
            corrPerCell = correlationVectorMatrix(reshape(ti.shapes(:,:,trueTransNumber),[],1),...
                reshape(handles.se.profiles(win(1):win(2),win(3):win(4),:),[],handles.se.nCells));
            [~,bestMatch] = max(corrPerCell);
            str = num2str(bestMatch);
            if bestMatch == handles.thisCell, col = seudo.pickColor('true');
            else  col = seudo.pickColor('false'); end
        case 4 % peak frame
            T = ti.times(trueTransNumber,:);
            [~,peakFrame] = max(handles.se.tcDefault.tc(T(1):T(2),handles.thisCell));
            peakFrame = peakFrame + T(1) - 1;
            str = num2str(peakFrame);
        case 5 % residual ratio
            resRatios = handles.transStats(handles.thisCell).resRatios;
            thisresRatio = resRatios(trueTransNumber);% / max(resRatios);
            str = sprintf('%0.3f',thisresRatio);
        case 6 % none
            str = '';
    end
    set(handles.axesTransSubplotTitles(tt),'string',str,'color',col)
end


% delete remaining plots
if isempty(tt), tt = 0; end
for pp = tt+1:length(handles.axesTransSubplots)
    %delete(handles.axesTransSubplots(pp))
    cla(handles.axesTransSubplots(pp))
    set(handles.axesTransSubplots(pp),'visible','off')
    set(handles.axesTransSubplotTitles(pp),'string','')
end



% just be sure all the display items are up to date
set(handles.editTextCellNumber,'String',num2str(handles.thisCell))
set(handles.sliderCellNumber,'Value',handles.thisCell)


% set value of artifact checkbox
set(handles.checkboxIsArtifact,'value',...
    handles.se.(handles.tcStructSpec{1})(handles.tcStructSpec{2}).transientInfo(handles.thisCell).isArtifact~=0);



% note # true, false transients
updateCellCount(handles)


% save values
guidata(hObject, handles);


% return focus to slider
uicontrol(handles.defaultFocus)






% CONSTANTS

function col = cellHighlight
col = [0 1 1];

function lw = borderLineWidth
lw = 5;


% --- Executes on button press in pushbuttonPlotTimeCourse.
function pushbuttonPlotTimeCourse_Callback(hObject, eventdata, handles)

handles.se.plotTransients(handles.thisCell)
