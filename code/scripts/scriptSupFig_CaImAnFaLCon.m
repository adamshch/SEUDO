

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get all the filenames

path_to_caiman = '/home/adam/GITrepos/CaImAn-MATLAB/';                     % Point to where CaImAn is located
addpath(path_to_caiman);                                                   % Add utilities path
addpath(genpath([path_to_caiman,'utilities']));                            % Add utilities path
addpath(genpath([path_to_caiman,'deconvolution']));                        % Add deconvolution path

base_dir  = '/mnt/bucket/adam/SEUDOdata/caiman/';                          % Path to point to data directory
listing   = dir(base_dir);                                                 % Get list of files in the directory
filename  = cell(numel(listing),1);
num_files = 0;
for kk = 1:numel(filename)
    TMP_name = listing(kk).name;                                           % filename to be processed
    if (~listing(kk).isdir)&&(strcmp(TMP_name(end-3:end),'.mat'))          % If not a directory and ends in .mat, add it to the list
        num_files           = num_files + 1;                               % Increase the number of found mat files
        filename{num_files} = TMP_name;                                    % Save the filename
    end
end
filename = filename(1:num_files);                                          % Store only the relevant file-names

clear listing TMP_name kk num_files

% Load up the CAIMAN-processed data 
m_cfrac = [];
m_rfrac = [];
se_cell = cell(numel(filename),1);

for kk = 1:numel(filename)
    fprintf('Loading block %d...',kk)
    seudoData   = load([base_dir,filename{kk}]);                             % Load data
    se_cell{kk} = seudo(seudoData.CNM.Y, seudoData.caiman_out.rois,'timeCourses',...
                        seudoData.CNM.C','name',['caiman_',filename{kk}]); % Make seudo object
%     CA = se_cell{kk}.autoClassifyTransients('default','saveResults',false,...
%                       'plotComparison',false,'radius',[3 7],'alpha',0.01,'tauBounds',[0.05,0.1]); % automatically classify transients
    CA = se_cell{kk}.autoClassifyTransients('default','saveResults',false,'plotComparison',false,'method', 'corr', 'corrThresh', 0.5);
    for ll = 1:numel(CA)
        m_cfrac = cat(1,m_cfrac,CA(ll).cfrac);
        m_rfrac = cat(1,m_rfrac,CA(ll).rfrac);
    end
    clear seudoData se CA
    fprintf('done.\n')
end

clear kk ll
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get cell counts

cell_cts = zeros(numel(se_cell),1);
for kk = 1:numel(se_cell)
    cell_cts(kk) = numel(se_cell{kk}.tcDefault.transientInfo);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Recalculate FaLCon plot points

m_cfrac       = [];
m_rfrac       = [];
countPerCell  = [];
for kk = 1:numel(filename)
    fprintf('Loading block %d...',kk)
%     CA        = se_cell{kk}.autoClassifyTransients('default','saveResults',false,...
%                       'plotComparison',false,'radius',[3 7],'alpha',0.01,'tauBounds',[0.05,0.1]); % automatically classify transients
    CA = se_cell{kk}.autoClassifyTransients('default','saveResults',...
                     false, 'plotComparison',false,'method', 'corr',...
                                'ignoreZeros',false, 'corrThresh', 0.5);  % Claculate automatic classification
    
    se_cell{kk}.autoClassifyTransients('default','saveResults',...
                     true, 'overwrite', true, 'plotComparison',false,'method', 'corr',...
                                'ignoreZeros',false, 'corrThresh', 0.5);  % Claculate automatic classification
    
    m_cfrac = cat(1,m_cfrac,[CA.cfrac].');
    m_rfrac = cat(1,m_rfrac,[CA.rfrac].');
    
    for cc = 1:se_cell{kk}.nCells
        C             = se_cell{kk}.tcDefault.transientInfo(cc).classification;
        countPerCell  = cat(1, countPerCell, [sum(C==1) sum(C==0) sum(C==-1), sum(isnan(C))]);
    end
    clear CA
    fprintf('done.\n')
end

clear kk ll

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

plotFaLCon(m_cfrac, m_rfrac, 'counts', sum(countPerCell,2))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

