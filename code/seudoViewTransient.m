function varargout = seudoViewTransient(varargin)
% SEUDOVIEWTRANSIENT MATLAB code for seudoViewTransient.fig
%      SEUDOVIEWTRANSIENT, by itself, creates a new SEUDOVIEWTRANSIENT or raises the existing
%      singleton*.
%
%      H = SEUDOVIEWTRANSIENT returns the handle to a new SEUDOVIEWTRANSIENT or the handle to
%      the existing singleton*.
%
%      SEUDOVIEWTRANSIENT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEUDOVIEWTRANSIENT.M with the given input arguments.
%
%      SEUDOVIEWTRANSIENT('Property','Value',...) creates a new SEUDOVIEWTRANSIENT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before seudoViewTransient_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to seudoViewTransient_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help seudoViewTransient

% Last Modified by GUIDE v2.5 04-Jan-2019 22:21:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @seudoViewTransient_OpeningFcn, ...
                   'gui_OutputFcn',  @seudoViewTransient_OutputFcn, ...
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


% --- Executes just before seudoViewTransient is made visible.
function seudoViewTransient_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to seudoViewTransient (see VARARGIN)

% Choose default command line output for seudoViewTransient
handles.output = hObject;


% store inputs
handles.movieOrig = varargin{1};
handles.whichFrames = varargin{2};
handles.profile = varargin{3};
handles.nameString = varargin{4};

if size(handles.movieOrig,4) > 1
    handles.sourceFitLSQ = handles.movieOrig(:,:,:,2);
    handles.sourceFitSEUDO = handles.movieOrig(:,:,:,3);
    handles.kernelsFitSEUDO = handles.movieOrig(:,:,:,4);
    handles.movieOrig = handles.movieOrig(:,:,:,1);
    handles.isSEUDO = true;
else
    handles.isSEUDO = false;
end

set(hObject,'name',handles.nameString)

% note # frames
handles.nFrames = size(handles.movieOrig,3);


% make subplots
posFig = [1 1 1 1];
if handles.isSEUDO
    nX = 3;
    nY = 2;
else
    nX = 3;
    nY = 1;
end

% create one big subplot, just to be sure previous transient plots
%subplot('position',get(handles.axesPlotArea,'position')./posFig([3 4 3 4]),...
%    'parent',get(handles.axesPlotArea,'parent'))

% make subplots
[~,handles.axesSubplots] = makeSubplots(get(handles.axesPlotArea,'parent'),...
    nX,nY,.1,.2,get(handles.axesPlotArea,'position')./posFig([3 4 3 4]));
delete(handles.axesPlotArea)

% set up for image plotting
for h = [handles.axesSubplots(:); handles.axesFrames]'
    axis(h,'image')
    axis(h,'ij')
    axis(h,'off')
    hold(h,'on')
    colormap(h,'gray')
end


% assign axes
if handles.isSEUDO
    handles.axesSourceProfile = handles.axesSubplots(1,1);
    handles.axesMovie = handles.axesSubplots(1,3);
    handles.axesMovieOverlay = handles.axesSubplots(2,3);
    handles.axesSourceFitLSQ = handles.axesSubplots(1,2);
    handles.axesSourceFitSEUDO = handles.axesSubplots(2,2);
    handles.axeskernelsFitSEUDO = handles.axesSubplots(2,1);
    handles.axesToClear = handles.axesSubplots(2:6);
else
    handles.axesSourceProfile = handles.axesSubplots(1);
    handles.axesMovie = handles.axesSubplots(2);
    handles.axesMovieOverlay = handles.axesSubplots(3);
    handles.axesToClear = handles.axesSubplots(2:3);
end

% add titles
title(handles.axesSourceProfile,'source profile')
title(handles.axesMovie,'movie')
title(handles.axesMovieOverlay,'movie with overlay')
if handles.isSEUDO
    title(handles.axesSourceFitLSQ,'source fit LSQ')
    title(handles.axesSourceFitSEUDO,'source fit SEUDO')
    title(handles.axeskernelsFitSEUDO,'kernels')
    title(handles.axesFrames,'movie frames')
end

% show profile
imagesc(handles.profile,'parent',handles.axesSourceProfile)

% get contour
handles.contour = contourc(double(handles.profile>0),[1 1]);

% set slider range for frames
set(handles.sliderFrame,'Max',handles.nFrames)
set(handles.sliderFrame,'Min',1)
if handles.nFrames > 1
    set(handles.sliderFrame,'SliderStep',[1 10]/(handles.nFrames-1))
else 
end
set(handles.sliderFrame,'value',1)


% initialize user-entered values
handles.convKernel = [];
handles.oldBlurRadius = 0;
handles.movie = handles.movieOrig;
handles.currFrameAvg = 1;
set(handles.editFrameAvg,'string',num2str(handles.currFrameAvg))


% initialize min and max display levels
handles.minValueLine = [];
handles.maxValueLine = [];

% show thumbnails
handles = updateMovieValues(hObject,handles);

% show current frame
updateFramePlot(hObject,handles);

% UIWAIT makes seudoViewTransient wait for user response (see UIRESUME)
% uiwait(handles.figure1);




% --- Outputs from this function are returned to the command line.
function varargout = seudoViewTransient_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% SLIDER TO MOVE BETWEEN FRAMES

function sliderFrame_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function sliderFrame_Callback(hObject, eventdata, handles)
updateFramePlot(hObject,handles);




% SET FRAME AVERAGE

function editFrameAvg_Callback(hObject, eventdata, handles)

% validate/constrain input
newInput = str2double(get(handles.editFrameAvg,'string'));
if newInput < 1
    newInput = 1;
elseif length(newInput) ~= 1 || ~isfinite(newInput) || ~isreal(newInput)
    newInput = handles.currFrameAvg;
end
set(handles.editFrameAvg,'string',num2str(newInput));
handles.currFrameAvg = newInput;

% update movie and display
updateMovieValues(hObject,handles);

function editFrameAvg_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% SET SPATIAL BLUR

function editBlur_Callback(hObject, eventdata, handles)

% validate/constrain input
newInput = str2double(get(handles.editBlur,'string'));
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

% store changes in figure
guidata(hObject, handles);

% update movie and display
updateMovieValues(hObject,handles);


function editBlur_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% UPDATE HISTOGRAM OF PIXEL VALUES

function handles = updatePixelValuesHistogram(hObject,handles)

h = handles.axesPixelValues;

cla(h)
handles.pixelValues = handles.movie(:);
hist(h,handles.pixelValues,100)

% set min and max values
handles.minDisplayValue = min(handles.movie(:));
handles.maxDisplayValue = max(handles.movie(:));

set(h,'tickdir','out','ytick',[])
ylabel(h,'# pixels')
xlabel(h,'fluorescence value')
axis(h,'tight')
box(h,'off')
hold(h,'on')






% SET MIN AND MAX VALUES FOR DISPLAY

function pushbuttonMinVal_Callback(hObject, eventdata, handles)

% get value for minimum
axes(handles.axesPixelValues)
[x,~] = ginput(1);

% update value in 
handles.minDisplayValue = max(min(handles.pixelValues),x);

% save changes to figure
%guidata(hObject, handles);

% update display of current frame
updateFramePlot(hObject,handles)


function pushbuttonMaxVal_Callback(hObject, eventdata, handles)

% get value for minimum
axes(handles.axesPixelValues)
[x,~] = ginput(1);

% update value in 
handles.maxDisplayValue = min(max(handles.pixelValues),x);

% save changes to figure
%guidata(hObject, handles);

% update display of current frame
updateFramePlot(hObject,handles)






% UPDATE MOVIE VALUES

function handles = updateMovieValues(hObject,handles)

% do frame averaging
if handles.currFrameAvg > 1
    handles.movie = convn(handles.movieOrig,ones(1,1,handles.currFrameAvg)/handles.currFrameAvg,'same');
else
    handles.movie = handles.movieOrig;
end


% go gaussian blurring
if ~isempty(handles.convKernel)
    handles.movie = convn(handles.movie,handles.convKernel,'same');
end


% update thumbnails
nT = 5;
whichFrames = round(linspace(1,handles.nFrames,nT));
thumbIms = seudo.makeThumbnailMatrix(handles.movie(:,:,whichFrames),'nY',1,'padVal',inf);
cla(handles.axesFrames)
imagesc(thumbIms,'parent',handles.axesFrames)


% update histogram of values
handles = updatePixelValuesHistogram(hObject,handles);

% update display of current frame
handles = updateFramePlot(hObject,handles);

% save changes to figure
guidata(hObject, handles);





% UPDATE PLOT TO SHOW CURRENT FRAME

function handles = updateFramePlot(hObject,handles)


% note current frame
currFrame = round(get(handles.sliderFrame,'value'));

% set color range
%colorRange = [min(handles.movie(:)) max(handles.movie(:))];
colorRange = [handles.minDisplayValue handles.maxDisplayValue];

% show min and max value

% remove old min and max lines
delete(handles.minValueLine)
delete(handles.maxValueLine)

h = handles.axesPixelValues;
handles.minValueLine = plot(h,[1 1]*handles.minDisplayValue,ylim(h),'r-');
handles.maxValueLine = plot(h,[1 1]*handles.maxDisplayValue,ylim(h),'g-');
%handles.minDisplayValue

% get image of transient
transIm = handles.movie(:,:,currFrame);

% clear axes
for h = handles.axesToClear, cla(h), end


% plot current frame
imagesc(transIm,'parent',handles.axesMovie,colorRange)
imagesc(transIm,'parent',handles.axesMovieOverlay,colorRange)

% add contour
col = [0 1 1];
hold(handles.axesMovieOverlay,'on')
plot(handles.axesMovieOverlay,handles.contour(1,2:end),handles.contour(2,2:end),...
    '-','linewidth',2,'color',col);
%title(handles.axesSubplots(2),sprintf('movie\nframe %d',currFrame))

% add frame number label
set(handles.textFrameNumber,'string',sprintf('frame %d (%d of %d)',...
    handles.whichFrames(currFrame),currFrame,length(handles.whichFrames)))

if handles.isSEUDO
    imagesc(handles.sourceFitLSQ(:,:,currFrame),'parent',handles.axesSourceFitLSQ,colorRange)
    imagesc(handles.sourceFitSEUDO(:,:,currFrame),'parent',handles.axesSourceFitSEUDO,colorRange)
    imagesc(handles.kernelsFitSEUDO(:,:,currFrame),'parent',handles.axeskernelsFitSEUDO,colorRange)
end


% save changes to figure
guidata(hObject, handles);

% return focus to slider
uicontrol(handles.sliderFrame)


