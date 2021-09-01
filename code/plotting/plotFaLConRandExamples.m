function varargout = plotFaLConRandExamples(se, m_cfrac, m_rfrac, varargin)

% plotFaLConRandExamples(se_cell, m_cfrac, m_rfrac, varargin)
%
% Code to plot example source and transient profiles from different the
% four quadrants of the FaLCon plot.
%
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;
p.addParameter('ex_cells'     , []);                                       % Optional list of cells to plot
p.addParameter('cell_cts'     , []);                                       % Optional likse of cell counts
p.addParameter('trans_list'   , []);                                       % Optional list of transients to plot
p.addParameter('plotClassOnly', false);                                    % True/false - whether to only plot classified transients
p.addParameter('figNum'       , []);                                        % Optional figure number to plot to

parse(p,varargin{:});

params = p.Results;
clear p

ex_cells      = params.ex_cells;
cell_cts      = params.cell_cts;
trans_list    = params.trans_list;
plotClassOnly = params.plotClassOnly;

if ~iscell(se)                                                             % This block makes sure that se is a cell
    TMP   = se;
    se    = cell(1);
    se{1} = TMP;
    clear TMP
end

if plotClassOnly                                                           % This block 
    classCells = [];
    for kk = 1:numel(se)
        for ll = 1:se{kk}.nCells
            classCells = cat(2,classCells,...
              any(~isnan(se{kk}.tcDefault.transientInfo(ll).classification))); % 
        end
    end
else
    classCells = true(size(m_cfrac));
end

% max_src = numel(m_rfrac);
m_rfrac = m_rfrac(classCells==1);                                          % Get only the residual fractions of the classified cells (if requested)
m_cfrac = m_cfrac(classCells==1);                                          % Get only the contamination fractions of the classified cells (if requested)

if ~exist('ex_cells','var')||isempty(ex_cells)
%     UL = find(m_rfrac>0.5&m_cfrac<0.3);
%     UR = find(m_rfrac>0.5&m_cfrac>0.5);
%     LL = find(m_rfrac<0.3&m_rfrac>0&m_cfrac<0.2);
%     LR = find(m_rfrac<0.4&m_cfrac>0.6);
    
    UL = find(m_rfrac>0.5&m_cfrac<0.3);
    UR = find(m_rfrac>0.5&m_cfrac>0.5);
    LL = find(m_rfrac<0.3&m_rfrac>0&m_cfrac<0.2);
    LR = find(m_rfrac<0.4&m_cfrac>0.6);
    
    if ~isempty(UL)
        ex_cells(1) = randsample(UL,1);
    else
        ex_cells(1) = randsample(1:numel(m_rfrac),1);
    end    
    if ~isempty(UR)
        ex_cells(2) = randsample(UR,1);
    else
        ex_cells(2) = randsample(1:numel(m_rfrac),1);
    end
    if ~isempty(LL)
        ex_cells(3) = randsample(LL,1);
    else
        ex_cells(3) = randsample(1:numel(m_rfrac),1);
    end
    if ~isempty(LR)
        ex_cells(4) = randsample(LR,1);
    else
        ex_cells(4) = randsample(1:numel(m_rfrac),1);
    end
    clear UL UR LL LR
else
%     if plotClassOnly
%         TMP_IDX     = 1:max_src;
%         TMP_IDX     = TMP_IDX(classCells==1);
%         for ll = 1:numel(ex_cells)
%             ex_cells(ll) = find(TMP_IDX == ex_cells(ll));
%         end
%         
%     end
end

classColors = [0  , 8  , 228;
               80 , 200, 80 ;
               255, 4  , 20 ;
               200, 200, 200]/256;                                         % Set up array of colors corresponding to different classifications

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial calculations, find the cells in the SEUDO object

if isempty(cell_cts)
    cell_cts = zeros(numel(se),1);
    for kk = 1:numel(se)
        cell_cts(kk) = numel(se{kk}.tcDefault.transientInfo);
    end
end
ec_batch = 1+sum(bsxfun(@gt, ex_cells, cumsum(cell_cts(:))),1);            % Get the batch number for each cell
in_batch = ex_cells - sum(bsxfun(@times, cell_cts(:), bsxfun(@gt, ...
                                      ex_cells, cumsum(cell_cts(:)))),1);  % Figure out which frame inside that file to load


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the various examples

if ~isempty(params.figNum)
    h = figure(params.figNum);
else
    h = figure();
end

for ll = 1:4
    if ~isempty(trans_list)
        tran_plot = trans_list(:,ll);
    else
        tran_plot = 1:size(se{ec_batch(ll)}.tcDefault.transientInfo(in_batch(ll)).shapes,3);
        if plotClassOnly
            tran_plot = tran_plot(~isnan(se{ec_batch(ll)}.tcDefault.transientInfo(in_batch(ll)).classification));
        end
        % tran_plot = randperm(size(se{ec_batch(ll)}.tcDefault.transientInfo(in_batch(ll)).shapes,3));
        tran_plot = tran_plot(randperm(numel(tran_plot)));
    end
    
    profSz    = size(se{ec_batch(ll)}.tcDefault.transientInfo(in_batch(ll)).shapes(:,:,1));
    tran_plot = sort(tran_plot(1:min(numel(tran_plot),12)),'ascend');
    class_inf = se{ec_batch(ll)}.tcDefault.transientInfo(in_batch(ll)).classification(tran_plot);
    
    [class_inf, IX] = sort(class_inf,'descend');
    IXnan           = isnan(class_inf);
    class_inf       = cat(1,vec(class_inf(~IXnan)), vec(class_inf(IXnan)));
    IX              = cat(1,vec(IX(~IXnan)), vec(IX(IXnan)));
    tran_plot       = tran_plot(IX);
    
    I_tmp     = basis2img2(reshape(se{ec_batch(ll)}.tcDefault.transientInfo(in_batch(ll)).shapes(:,:,tran_plot),[],numel(tran_plot)),...
        profSz,[3,4], true, class_inf);

    W        = se{ec_batch(ll)}.tcDefault.transientInfo(in_batch(ll)).window;
    prof_tmp = se{ec_batch(ll)}.profiles(W(1):W(2),W(3):W(4),in_batch(ll));
    prof_tmp = basis2img2([ones(numel(prof_tmp),1),vec(prof_tmp),ones(numel(prof_tmp),1)],...
        size(se{ec_batch(ll)}.tcDefault.transientInfo(in_batch(ll)).shapes(:,:,1)),[3,1],true);
    prof_tmp(1:floor(size(prof_tmp,1)/3),:,:)    = 1;
    prof_tmp(ceil(2*size(prof_tmp,1)/3):end,:,:) = 1;
    I_tmp    = cat(2, prof_tmp, ones(size(prof_tmp,1),10,3), I_tmp);
    
    dotOffset = [round(0.2*profSz(1)),round(0.1*profSz(2))];
    dotOffset = mean(dotOffset)*[1,1];
    [X,Y] = meshgrid(13+profSz(2)*(2:5)-dotOffset(2), 1+dotOffset(1)+profSz(1)*(0:2));
    X = vec(X.');
    X = X(1:numel(class_inf));
    Y = vec(Y.');
    Y = Y(1:numel(class_inf));

    ha = subplot(2,2,ll);
    imagesc(I_tmp)
    subplot(2,2,ll), hold on;
    
    scatter(X(class_inf== 1)   , Y(class_inf    ==1),'o','filled','MarkerEdgeColor',classColors(1,:),'MarkerFaceColor',classColors(1,:));
    scatter(X(class_inf== 0)   , Y(class_inf    ==0),'o','filled','MarkerEdgeColor',classColors(2,:),'MarkerFaceColor',classColors(2,:));
    scatter(X(class_inf==-1)   , Y(class_inf   ==-1),'o','filled','MarkerEdgeColor',classColors(3,:),'MarkerFaceColor',classColors(3,:));
    scatter(X(isnan(class_inf)), Y(isnan(class_inf)),'o','filled','MarkerEdgeColor',classColors(4,:),'MarkerFaceColor',classColors(4,:));
    subplot(2,2,ll), hold off;
    title(sprintf('Cell number %d (batch %d/source %d)', ex_cells(ll),...
                                              ec_batch(ll), in_batch(ll)))
    axis image; colormap gray
    set(gca, 'XTick',[],'YTick',[],'xcolor','none','ycolor','none')
%     xlabel(sprintf('FaLCon coordinates: (%1.3f,%1.3f)', ...
%                             m_cfrac(ex_cells(ll)), m_rfrac(ex_cells(ll))))
    ha.XAxis.Label.Color=[0 0 0];
    ha.XAxis.Label.Visible='on';
    box off
end

if nargout == 1
    varargout{1} = h;
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extra function for plotting multiple images with the same baseline

function I = basis2img2(A, dimI, dims, varargin)

% function I = basis2img2(A, dims)
%
% A   - matrix of image patches
% dim - a 2x1 vector specifying number of rows and cols
%       for reshaping A
% dimI - a 2x1 vector specifying number of rows and cols
%       for reshaping each of the columns of A
% outputs I, an RGB image.  red dots indicate the relative max 
% pixel intensity level.


if nargin > 3
    nonneg = varargin{1};
else
    nonneg = false;
end

% if nargin > 4
%     classInfo = varargin{2};
% else
%     classInfo = [];
% end
% 
% if nargin > 5
%     classColors = varargin{3};
% else
%     classColors = [0  , 8  , 228;
%                255, 4  , 20 ;
%                80 , 200, 80 ;
%                200, 200, 200]/256;                                         % Set up array of colors corresponding to different classifications
% end


N = size(A,1);
% M = size(A,2);

  
spacing = 1;
if nargin < 3
  a = floor(sqrt(size(A,2)));
  if size(A,2) > a*(a+1),  nrows = a+1;
  else                  ,  nrows = a; end
  ncols = ceil(sqrt(size(A,2)));
else
  nrows = dims(1); ncols = dims(2); 
end

if nargin < 2
    wx = sqrt(N);
    wy = sqrt(N);
else
    wx = dimI(1);
    wy = dimI(2);
end
  
I = zeros(wx*nrows+(nrows+1)*spacing, wy*ncols+(ncols+1)*spacing,3);


for i = 1:size(A,2)
  patch = reshape(A(:,i), wx, wy);
  sx    = (wy + spacing) * rem(i-1,ncols) + 1 + 1;
  sy    = (wx + spacing) * floor((i-1)/ncols) + 1 + 1;

  prange = max(abs(patch(:)));
  if prange > .00001
      if nonneg
          patch(patch<0) = 0;
          patch = patch / prange; 
      else
          patch = .5 + .5*patch / prange; 
      end
  else
      if nonneg
          patch = patch + .5;
      else
      end
  end 
%  
%   if ~isempty(classInfo)
%     if classInfo(i) == seudo.valTrue
%        I((sy-1):(sy+wx),(sx-1):(sx+wy),:) = repmat(reshape(classColors(1,:) ,[1,1,3]),[wx+2,wy+2,1]); % true
%     elseif isnan(classInfo(i))
%        I((sy-1):(sy+wx),(sx-1):(sx+wy),:) = repmat(reshape(classColors(4,:) ,[1,1,3]),[wx+2,wy+2,1]); % unclassified
%     elseif classInfo(i) == seudo.valFalse
%        I((sy-1):(sy+wx),(sx-1):(sx+wy),:) = repmat(reshape(classColors(2,:) ,[1,1,3]),[wx+2,wy+2,1]); % false
%     else
%        I((sy-1):(sy+wx),(sx-1):(sx+wy),:) = repmat(reshape(classColors(3,:) ,[1,1,3]),[wx+2,wy+2,1]); % mixed
%     end 
%   end


  I(sy:(sy+wx-1), sx:(sx+wy-1),:) = repmat(patch,[1 1 3]);
end

end

% ex_cells = [475, 202, 1708, 446];                                          % Example cells in: UL,UR,LL, LR,
% trans_list = [04,10,13,20,28,26,02,21;
%               12,13,04,05,03,02,10,07;
%               26,38,12,11,20,01,28,27;
%               01,02,03,04,05,06,07,08].'; % Transients for [475, 202, 1708, 446]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
