function se = parallelSEUDO(se, whichCells, varargin)

% function se = parallelSEUDO(se, whichCells, varargin)
%
% Function to run SEUDO on all cells in parallel, restricted to the
% transient times (to speed up the altorithm). 
%
% 2020 - Adam Charles & Jeff Gauthier

if isempty(whichCells); whichCells = 1:se.nCells; end                      % In the case of an empty whichCells, run all cells
nCells = numel(whichCells);                                                % Record the number of cells to run

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop through each source
results   = cell(nCells,1);
startTime = now;
seMini    = cell(nCells,1);

for cc = 1:nCells 
    keepFrames = false(se.movF,1);                                         % identify transient times, keep only those time points
    curCell    = whichCells(cc);                                           % Set the current cell ID
    ti         = se.tcDefault.transientInfo(curCell);                      % Get the time-course for the current cell
    win        = se.tcDefault.transientInfo(curCell).window;               % Get the spatial window for the current cell
    for tt = 1:size(ti.times,1)                                            % Loop through the transient IDs
        keepFrames(ti.times(tt,1):ti.times(tt,2)) = true;                  % Set the selection to pick out all of the frames with transients
    end
    if ~any(keepFrames), continue, end                                     % skip if no transients
    
    seMini{cc} = seudo(se.movieSpec(win(1):win(2), win(3):win(4), keepFrames),...
              se.profiles(win(1):win(2), win(3):win(4),:),...
                          'timeCourses',se.tcDefault.tc(keepFrames,:),...
                                            'computeTransientInfo',false); % package into seudo object
end
parfor cc = 1:nCells 
    seMini{cc}.estimateTransientsOnlyWithSEUDO('whichCells', whichCells(cc), varargin{:});   % estimate time courses with SEUDO
    results{cc}.tcSeudo   = seMini{cc}.tcSeudo(end);                       % store result
    results{cc}.tcDefault = seMini{cc}.tcDefault;
    fprintf('finished %d of %d cells in %0.1f min\n',cc,nCells,...
                                                    24*60*(now-startTime)) % show progress
end

clear seMini

nFrames           = se.movF;
thisTcSeudo       = struct;
thisTcSeudo.tc    = nan(se.movF,se.nCells);
thisTcSeudo.tcLSQ = nan(se.movF,se.nCells);

for cc = 1:nCells
    curCell    = whichCells(cc);
    analyzedFrames = false(nFrames,1);
    transTimes     = se.tcDefault.transientInfo(curCell).times;
    for rr = 1:size(transTimes,1)
        analyzedFrames(transTimes(rr,1):transTimes(rr,2)) = true;
    end

    assert(sum(analyzedFrames) == length(results{cc}.tcSeudo.tc))

    thisTcSeudo.tc(analyzedFrames,curCell)    = results{cc}.tcSeudo.tc;
    thisTcSeudo.tcLSQ(analyzedFrames,curCell) = results{cc}.tcSeudo.tcLSQ;

    thisTcSeudo.extras{curCell} = results{cc}.tcSeudo.extras;
    thisTcSeudo.params{curCell} = results{cc}.tcSeudo.params;

    thisTcSeudo.extras{curCell}.lsqCosts = nan(se.movF,1);
    thisTcSeudo.extras{curCell}.lsqCosts(analyzedFrames) = results{cc}.tcSeudo.extras.lsqCosts;

    thisTcSeudo.extras{curCell}.blobCosts = nan(se.movF,1);
    thisTcSeudo.extras{curCell}.blobCosts(analyzedFrames) = results{cc}.tcSeudo.extras.blobCosts;

    thisTcSeudo.extras{curCell}.tcCellsWithBlobs = nan(se.movF,1);
    thisTcSeudo.extras{curCell}.tcCellsWithBlobs(analyzedFrames) = results{cc}.tcSeudo.extras.tcCellsWithBlobs;

%     thisTcSeudo.extras{curCell}.tcBlobs = nan(se.movF,size(results{cc}.tcSeudo.extras.tcBlobs{1},2));
%     thisTcSeudo.extras{curCell}.tcBlobs(analyzedFrames,:) = results{cc}.tcSeudo.extras.tcBlobs{1};
end
 
% load empty transient info
thisTcSeudo.transientInfo = struct;

% put this seudo struct in first, append any exisitng structs afterwards
se.tcSeudo = [thisTcSeudo; se.tcSeudo];   


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%