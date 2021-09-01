function varargout = seudoQuickTest(varargin)
% SEUDOQUICKTEST MATLAB code for seudoQuickTest.fig
%      SEUDOQUICKTEST, by itself, creates a new SEUDOQUICKTEST or raises the existing
%      singleton*.
%
%      H = SEUDOQUICKTEST returns the handle to a new SEUDOQUICKTEST or the handle to
%      the existing singleton*.
%
%      SEUDOQUICKTEST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEUDOQUICKTEST.M with the given input arguments.
%
%      SEUDOQUICKTEST('Property','Value',...) creates a new SEUDOQUICKTEST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before seudoQuickTest_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to seudoQuickTest_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help seudoQuickTest

% Last Modified by GUIDE v2.5 14-Dec-2018 00:23:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @seudoQuickTest_OpeningFcn, ...
                   'gui_OutputFcn',  @seudoQuickTest_OutputFcn, ...
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


% --- Executes just before seudoQuickTest is made visible.
function seudoQuickTest_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to seudoQuickTest (see VARARGIN)



% To launch, the GUI is provided with a set of miniature seudo objects.
% Each object has a movie that  only includes one transient and the pixels near the cell



% store inputs in handles struct
handles.minis = varargin{1};
handles.cellSpec = varargin{2};
handles.name = varargin{3};

% set figure name
set(hObject,'name',handles.name)



% hide axesTransientsArea, since it is only used to set bounds for subplots
set(handles.axesPlotArea,'visible','off')

set(handles.textCalculating,'visible','off')

% make subplots
nX = length(handles.minis);
nY = 7;

%posFig = get(hObject,'position');
posFig = [0 0 1 1];

% create one big subplot, just to be sure previous transient plots
subplot('position',get(handles.axesPlotArea,'position')./posFig([3 4 3 4]),...
    'parent',hObject)

% create subplots for each cell/graphic
[~,handles.axesTransSubplots] = makeSubplots(hObject,nX,nY,.1,.3,get(handles.axesPlotArea,'position')./posFig([3 4 3 4]));


% delete any existing buttons
delete(findall(get(hObject,'child'),'style','pushbutton','string','movie'))

% make buttons to show movie of each transient
for xx=1:nX
    
    h = handles.axesTransSubplots(7,xx);
    
    set(h,'units','pixels')
    pos = get(h,'position');
    delete(handles.axesTransSubplots(7,xx))
    
    % set size
    bWidth = min(pos(3),120);
    bHeight = min(pos(4),40);
    
    pos(1) = pos(1) + pos(3)/2 - bWidth/2;
    pos(3) = bWidth;
    
    pos(2) = pos(2) + pos(4)/2 - bHeight/2;
    pos(4) = bHeight;
    
    h = uicontrol('Style','pushbutton','String','movie','fontsize',14,'fontweight','bold',...
        'Position',pos,'Callback',@(o,e)showMovie(o,xx));
    
    
    set(h,'units','normalized')
end



% update positions of text labels
% NEVERMIND MATLAB MAKES THIS TOO DIFFICULT
%setVertPos = @(x,p)[x(1) p(2) x(3:4)];
%set(handles.text1,'position',setVertPos(get(handles.text1,'position'),get(handles.axesTransSubplots(1,1),'position')))
%get(handles.text1,'position')
%get(handles.axesTransSubplots(2,1),'position')



% set colors of text labels
set(handles.text3b,'foregroundcolor',seudo.pickColor('LSQ'))
set(handles.text4a,'foregroundcolor',seudo.pickColor('SEcell'))
set(handles.text4b,'foregroundcolor',seudo.pickColor('SEblobs'))
set(handles.text11,'foregroundcolor',seudo.pickColor('true'))
set(handles.text12,'foregroundcolor',seudo.pickColor('false'))
set(handles.text6,'foregroundcolor',seudo.pickColor('SEblobs'))

%set(hObject,'color','w')

% color transient labels, if classification was provided
if size(handles.cellSpec,2) >= 3
    handles.cellClass = handles.cellSpec(:,3);
else
    handles.cellClass = [];
end


% choose default command line output for seudoQuickTest
handles.output = hObject;


% initialize list of previous parameter sets
handles.prevParams = {};
set(handles.popupmenuPrevParams,'string',{})


% compute seudo
handles = computeNewResults(hObject,handles);

% plot results
updatePlot(hObject,handles)

% choose horizontal zoom tool
set(zoom,'Enable','on','motion','horizontal');

% update handles structure
guidata(hObject, handles);






% DISPLAY MOVIE

function showMovie(buttonObject,whichTrans)

% get seudo object holding this cell
handles = guidata(get(buttonObject,'parent'));
se = handles.minis{whichTrans};



% get struct of SEUDO results
popupValue = get(handles.popupmenuPrevParams,'value');
tcStruct = handles.prevParams{popupValue,whichTrans}.tcStruct;




% assemble matrix of data

% profile
P = se.profiles(:,:,1);

% movie
M = se.movieSpec;
%maxVal = prctile(M(:),99.9);
%M(M>maxVal) = maxVal;

% account for dsTime
D = tcStruct.params.dsTime;
M = convn(M,ones(1,1,D));
M = M(:,:,1:se.movF);
movNorm      = ones(1,1,se.movF)*D;
movNorm(1:D) = 1:D;
M            = bsxfun(@rdivide,M,movNorm);

% least squares fit
L = bsxfun(@times,P,permute(se.tcSeudo(end).tcLSQ,[2 3 1]));

% seudo blobs
B  = nan(size(M));
nY = sum(sum(reshape(tcStruct.extras.whichPixels(:,1),se.movY,se.movX),2)>0);
nX = sum(sum(reshape(tcStruct.extras.whichPixels(:,1),se.movY,se.movX),1)>0);
apply_blobs = @(z) reshape(convn(reshape(z, nY,nX), tcStruct.extras.one_blob, 'same'),nX*nY, 1);
for ff = 1:se.movF
    B(:,:,ff) = reshape(apply_blobs(tcStruct.extras.tcBlobs{1}(ff,:)),nY,nX);
end
%B = B/sum(tcStruct.extras.blobFits(:,1));

% seudo cell
S = bsxfun(@times,P,permute(se.tcSeudo(end).tc,[2 3 1]));


% normalize profile
P = P/max(P(:))*max(M(:));


% collect data for showing in a new window
theMov = cat(4,M,L,S,B);
whichFrames = 1:size(theMov,3);
%guiName = 'test';
guiName = sprintf('transient %d from quickSEUDO',whichTrans);

% display
seudoViewTransient(theMov,whichFrames,P,guiName)



% if 0
%     
%     % single blob
%     b = zeros(size(P));
%     obY = size(tcStruct.extras.one_blob,1);
%     obX = size(tcStruct.extras.one_blob,2);
%     yStart = round(size(b,1)/2-obY/2);
%     xStart = round(size(b,2)/2-obX/2);
%     b(yStart+(0:obY-1),xStart+(0:obX-1)) = tcStruct.extras.one_blob;
%     b = b/max(b(:))*max(M(:));
%     
%     
%     % concatenate components
%     mov = [ repmat(P,1,1,se.movF) M L ; repmat(b,1,1,se.movF) B S ];
%     
%     % old viewer, deprecated
%     plotImagesSequentially(imageSeries(mov))
%     
% end


    
    
    
    
    
    
% --- Outputs from this function are returned to the command line.
function varargout = seudoQuickTest_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;





% SWITCH TO PREVIOUS PARAMETER SETS


function popupmenuPrevParams_Callback(hObject, eventdata, handles)

updatePlot(hObject,handles)


function popupmenuPrevParams_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end







% PERFORM SEUDO

function handles = computeNewResults(hObject,handles)


set(handles.textCalculating,'visible','on')
set(handles.textError,'visible','off')
drawnow


paramString = get(handles.editParams,'string');

    
% for each transient...
for ss = 1:length(handles.minis)
    
    
    % get the seudo object with this data
    se = handles.minis{ss};
    
    
    % perform SEUDO using the specified parameters
    
    tic
    try
        switch 2
            case 0 % version before windowed SEUDO was introduced
                eval(sprintf('se.estimateTimeCoursesWithSEUDO_apr27(''verbose'',0,''whichCells'',%d,%s);',1,paramString))
            case 1 % version before windowed SEUDO was introduced
                eval(sprintf('se.estimateTimeCoursesWithSEUDO_may17(''verbose'',0,''whichCells'',%d,%s);',1,paramString))
            case 2 % version after windowed SEUDO was introduced
                eval(sprintf('se.estimateTimeCoursesWithSEUDO(''verbose'',0,''saveBlobTimeCourse'',1,''padSpace'',inf,''whichCells'',%d,%s);',1,paramString))
        end
    catch err
        set(handles.textCalculating,'visible','off')
        set(handles.textError,'visible','on')
        err_ = struct;
        err_.message = err.message;
        err_.stack = err.stack(1);
        error(err_)
        %rethrow(err)
        %throwAsCaller(err)
    end
    compTime = toc;
    
    
    
    
    % save result in list of previous params
    resStruct = struct;
    resStruct.tcStruct = se.tcSeudo(end);
    resStruct.paramString = paramString;
    resStruct.compTime = compTime;
    if ss == 1
        handles.prevParams{size(handles.prevParams,1)+1,ss} = resStruct;
    else
        handles.prevParams{size(handles.prevParams,1),ss} = resStruct;
    end
    
    
end


% show parameters (if desired)
if get(handles.checkboxPrint,'value')
    fprintf('parameters: %s\n',handles.prevParams{end,ss}.paramString)
end


set(handles.textVersion,'string',sprintf('SEUDO code\nversion %s',num2str(se.tcSeudo(end).params.version)))

set(handles.textCalculating,'visible','off')


% add to list of accumulating param strings

% start with current list
newList = get(handles.popupmenuPrevParams,'string');

% add newline if needed
%if ~isempty(newList), newList = [newList sprintf('\n')]; end

% append param string
%newList = [newList paramString];
newList{end+1} = paramString;
set(handles.popupmenuPrevParams,'string',newList,'value',length(newList))



% save values
guidata(hObject, handles);




%  DISPLAY RESULT

function updatePlot(hObject,handles)

% clear previous plots
for h = reshape(handles.axesTransSubplots(1:6,:),1,[]), cla(h), end


popupValue = get(handles.popupmenuPrevParams,'value');
    

% for each transient...
for ss = 1:length(handles.minis)
    
    
    % get the seudo object with this data
    se = handles.minis{ss};
    
    % get struct of SEUDO results
    tcStruct = handles.prevParams{popupValue,ss}.tcStruct;
   
    % get dimensions
    nY = sum(sum(reshape(tcStruct.extras.whichPixels(:,1),se.movY,se.movX),2)>0);
    nX = sum(sum(reshape(tcStruct.extras.whichPixels(:,1),se.movY,se.movX),1)>0);
    
    apply_blobs = @(z) reshape(convn(reshape(z, nY,nX), tcStruct.extras.one_blob, 'same'),nY,nX, 1);
    
    if isfield(tcStruct.params,'saveBlobTimeCourse')
        if tcStruct.params.saveBlobTimeCourse
            tcBlobs = tcStruct.extras.tcBlobs{1};
            
            tcBlobsAll = nansum(tcBlobs,2)/sum(tcStruct.extras.blobFits(:,1));
            roi = reshape(se.profiles(:,:,1),[],1);
            roi = roi(tcStruct.extras.whichPixels(:,1));
            tcBlobsROI = nansum(tcBlobs(:,roi>0),2)/sum(tcStruct.extras.blobFits(:,1));
            imBlobSum = reshape(apply_blobs(max(tcBlobs,[],1)),nY,nX);
        else
            tcBlobsAll = tcStruct.extras.tcBlobs(:,1)/sum(tcStruct.extras.blobFits(:,1));
            tcBlobsROI = tcStruct.extras.tcBlobs(:,2)/sum(tcStruct.extras.blobFits(:,1));
            imBlobSum = [];
        end
    else
        roi = se.profiles(:,:,1);
        if isfield(tcStruct.extras,'blobFits')
            tcBlobsAll = sum(tcStruct.extras.tcBlobs,1) / sum(tcStruct.extras.blobFits);
            tcBlobsROI = sum(tcStruct.extras.tcBlobs(roi(:)>0,:),1) / sum(tcStruct.extras.blobFits);
        else
            tcBlobsAll = zeros(size(tcStruct.tcLSQ));
            tcBlobsROI = zeros(size(tcStruct.tcLSQ));
        end
        imBlobSum = reshape(apply_blobs(max(tcStruct.extras.tcBlobs,[],2)),se.movY,se.movX);
    end
    
    
    roi = se.profiles(:,:,1);
    roi = reshape(roi(tcStruct.extras.whichPixels),nY,nX);
    
    
    % get blob time course
    if tcStruct.params.saveBlobTimeCourse
        % if individual blob time courses were saved, combine across space to get the total
        
        extras = tcStruct.extras;
        
        % all blobs
        tcBlobsAll = nansum(extras.tcBlobs{1}(:,:),2)/sum(extras.blobFits(:,1));
        
        % blobs in the profile
        %roi = reshape(se.profiles(win(1):win(2),win(3):win(4),1),[],1);
        %tcBlobsInProfile = nansum(extras.tcBlobs(framesToPlot,roi>0),2)/sum(extras.blobFits(:,1));
        
        % image of the blob sum
        B = reshape(extras.tcBlobs{1}(:,:)',nY,nX,[]);
        imBlobs = max(convn(B, extras.one_blob, 'same'),[],3);
        
        % combine with cell
        imCell   = max(tcStruct.tcLSQ) * roi;
        vals     = [imBlobs(:); imCell(:)];
        colCell  = seudo.pickColor('SEcell');
        colBlobs = seudo.pickColor('SEblobs');
        %imBlobs = cat(3,imBlobs,imCell,imBlobs) / max(vals);
        M = [];
        for cc = 1:3
            %M(:,:,cc) = 1 - (1 - colCell(cc)*imCell) .* (1 - colBlobs(cc)*imBlobs);
            M = cat(3,M,1 - (1 - colCell(cc)*imCell) .* (1 - colBlobs(cc)*imBlobs));
        end
        imBlobs = M / max(vals(:));

    else
        % if individual blob time courses were not saved, combine across space to get the total
        tcBlobsAll = extras.tcBlobs(:,1)/sum(extras.blobFits(:,1));
        tcBlobsInProfile = extras.tcBlobs(framesToPlot,2)/sum(extras.blobFits(:,1));
        imBlobs = [];
    end
    
    imBlobSum = imBlobs;
    
    
    
    
    
    
    % plot results
    
    % cell shape
    h = handles.axesTransSubplots(1,ss);
    imagesc(roi,'parent',h)
    axis(h,'image')
    axis(h,'off')
    if ~isempty(handles.cellClass)
        switch handles.cellClass(ss)
            case seudo.valTrue,  col = seudo.pickColor('true');
            case 0, col = seudo.pickColor('false');
            otherwise,           col = 'k';
        end
    else
        col = 'k';
    end
    title(h,sprintf('cell %d, transient %d',handles.cellSpec(ss,1),handles.cellSpec(ss,2)),'fontsize',12,'color',col)
    colormap(h,'gray')
    
    % active shape
    h = handles.axesTransSubplots(2,ss);
    imagesc(se.tcDefault.transientInfo(1).shapes,'parent',h)
    axis(h,'image')
    axis(h,'off')
    colormap(h,'gray')
    
    
    % found cell time course
    h = handles.axesTransSubplots(3,ss);
    plot(h,tcStruct.tcLSQ,'color',seudo.pickColor('LSQ'),'linewidth',3)
    hold(h,'on')
    plot(h,tcStruct.tc,'color',seudo.pickColor('SEcell'),'linewidth',1)
    plot(h,tcBlobsROI,'--','color',seudo.pickColor('SEblobs'),'linewidth',1)
    axis(h,'tight')
    
    % note % preserved
    evalFcn_area = seudo.makeTimeCourseEvalFunction('method','area');
    title(h,sprintf('%d%% preserved',round(100*evalFcn_area(tcStruct.tcLSQ,tcStruct.tc))))
    
    % blob TC
    h = handles.axesTransSubplots(4,ss);
    plot(h,tcBlobsAll,'-','color',seudo.pickColor('SEblobs'),'linewidth',1)
    hold(h,'on')
    plot(h,tcBlobsROI,'--','color',seudo.pickColor('SEblobs'),'linewidth',1)
    axis(h,'tight')
    if isequal(ylim(h),[-1 1]), ylim(h,[0 1]),end
    
    % costs
    h = handles.axesTransSubplots(5,ss);
    tK = tcStruct.extras.lsqCosts;
    tR = tcStruct.extras.blobCosts;
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
    hold(h,'on')
    plot(h,find(~iK),tK(iK==0),'o','markersize',5,'color',seudo.pickColor('true'))
    plot(h,find(~iR),tR(iR==0),'o','markersize',5,'color',seudo.pickColor('false'))
    plot(h,find(iK),tK(iK>0),'.','markersize',15,'color',seudo.pickColor('true'))
    plot(h,find(iR),tR(iR>0),'.','markersize',15,'color',seudo.pickColor('false'))
    
    
    %axis(h,'tight')
    
    % blob max
    h = handles.axesTransSubplots(6,ss);
    imagesc(imBlobSum,'parent',h)
    axis(h,'image')
    set(h,'xtick',[],'ytick',[])
    colormap(h,'gray')
    %axis(h,'off')
    
    xlabel(h,sprintf('computed in %0.01f sec',handles.prevParams{popupValue,ss}.compTime))
    
end


% link axes
for ss = 1:length(handles.minis)
    linkaxes(handles.axesTransSubplots(3:5,ss),'x')
    linkaxes(handles.axesTransSubplots([1 2 6],ss))
end




% save values
guidata(hObject, handles);




% SHOW DEFAULT PARAMETERS

function pushbuttonShowDef_Callback(hObject, eventdata, handles)
handles.minis{1}.estimateTimeCoursesWithSEUDO('showDefault',true);


% CHECKBOX FOR PRINTING OUTPUT

function checkboxPrint_Callback(hObject, eventdata, handles)


% RUN SEUDO AND RE-PLOT

function pushbuttonGo_Callback(hObject, eventdata, handles)
handles = computeNewResults(hObject,handles);
set(handles.popupmenuPrevParams,'value',numel(get(handles.popupMenu, 'string')))
updatePlot(hObject,handles)

function editParams_Callback(hObject, eventdata, handles)
handles = computeNewResults(hObject,handles);
updatePlot(hObject,handles)

function editParams_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
