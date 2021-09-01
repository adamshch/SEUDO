function [tracesOut, params, varargout] = globalSEUDOWrapper(mov, spatialProfiles, timeTraces, varargin)

% [tracesOut, params, varargout] = globalSEUDOWrapper(mov, spatialProfiles, timeTraces, varargin)
%
% Wrapper to run SEUDO easily on a single dataset, incorporating speedups
% from parallelization and restricting computation to transient-active
% frames. The inputs are:
%   mov             - 3D array containing the fluorescence movie (3rd
%                        dimension is time)
%   spatialProfiles - 3D array containing the gray-scale spatial profiles
%                       (3rd dimension is Neuron ID)
%   timeTraces      - 2D array (number frames X number neurons) of the
%                       time-traces corresponding the spatial profiles
%   Extra inputs:
%   Extra inputs can be included in the format (..., 'name', value, ...)
%   Setting these properties over-rides both the default AND the automatic
%   settings based on data. These properties incude: 
%       lambdaBlob  - lambda value for SEUDO (default = 10)   
%       dsTime      - number of frames for running average(default = 3)   
%       sigma2      - noise variance for SEUDO (default = 0.0002)
%       padSpace    - number of pixels to pad with (default = 5)
%       blobRadius  - size of blob in SEUDO cell construction (default = 3)   
%       seName      - Optional name for SEUDO object (default = 'Unknown')
%       scaleFactor - Optional scale factor: auto computed otherwise
% 
% The outputs of this function are:
%
% tracesOut - 2D array (number frames X number neurons) of the
%             decontaminated time-traces corresponding the spatial profiles
% params    - Struct containing the parameters that ended up being used
%             (final result of automatic setting, user over-riding and
%             default parameters)
% extras    - OPTIONAL cell array containing structs that store extra infro
%             for the SEUDO results (e.g., when SEUDO cells were used etc.)
%
% 2020 - Adam Charles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input parsing

p = inputParser;
p.addParameter('padSpace'    , []);      
p.addParameter('blobRadius'  , []);   
p.addParameter('lambdaBlob'  , []);   
p.addParameter('dsTime'      , []);   
p.addParameter('sigma2'      , []);
p.addParameter('scaleFactor' , []);
p.addParameter('whichCells'  , []);                                        % Passing an empty array of cells will automatically run on all cells
p.addParameter('seName'      , 'Unknown');   
parse(p,varargin{:});
p = p.Results;

defaultPars.padSpace   = 5;                                                % We find that 5 pixels of padding is sufficient for SEUDO
defaultPars.blobRadius = 3;                                                % This value ensures that the blobs are not too big (looks too much like entire cells) or too small (looks like single-pixel noise)
defaultPars.lambdaBlob = 10;                                               % A reasonable parameter
defaultPars.dsTime     = 3;                                                % 3 frame averages are sufficient to get a decent estimate
defaultPars.sigma2     = 0.0020;                                           % A fairly low-noise default setting


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start a parallel pool if one does not exist

if isempty(gcp('nocreate')); parpool(); end                                % Check if parallel pool exists and make one if it doesn't

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start the process by making the SEUDO object

se = seudo(mov, spatialProfiles,'timeCourses',timeTraces,'name',p.seName); % Create the SEUDO object
clear mov spatialProfiles timeTraces

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up the parameters to run SEUDO

% whichCells = [];                                                           % Passing an empty array of cells will automatically run on all cells
parsAuto   = se.suggestParameters();                                       % Get a set of suggested parameters based on movie statistics
paNames    = fieldnames(parsAuto);                                         % Get the names of the automatically set parameters
pNames     = fieldnames(defaultPars);                                      % Get the names of the pre-set parameters

for ll = 1:numel(pNames)                                                   % Loop through all the parameter names...
    if any(strcmp(paNames, pNames{ll}))                                    % ... If a parameter is automatically set...
        defaultPars.(pNames{ll}) = parsAuto.(pNames{ll});                  % ... Replace that parameter in the parameter set to use.
    end
    if isempty(p.(pNames{ll}))                                             % Check if the parameter was unspecified by the user 
        p.(pNames{ll}) = defaultPars.(pNames{ll});                         % If so use ther default/suggested values
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run SEUDO in parallel mode

se.parallelSEUDO(p.whichCells, 'dsTime', p.dsTime, 'sigma2', p.sigma2, ...
              'lambdaBlob', p.lambdaBlob, 'blobRadius', p.blobRadius, ...
                         'padSpace', p.padSpace, 'saveBlobTimeCourse', 1); % Run SEUDO


baseOrigEst = nanmedian(se.tcDefault.tc,2);
baseSEUDO   = nanmedian(se.tcDefault.tc,2);

nonZeros = vec(sum(se.profiles,3))~=0;
tmpProf  = reshape(se.profiles,[],size(se.profiles,3));
tmpProf  = tmpProf(nonZeros,:);

if isempty(p.scaleFactor)
    reconMovie = tmpProf*(se.tcDefault.tc)';
    scaleFactor = mean(se.movieSpec(nonZeros,:))./mean(reconMovie);      % There may be a scale difference between the movie and time traces: compute that!
else
    scaleFactor = p.scaleFactor;
end



tracesOut = se.tcDefault.tc;                                               % Start with the old time traces
tracesOut(~isnan(se.tcSeudo.tc)) = (1./scaleFactor)*se.tcSeudo.tc(...
                                                   ~isnan(se.tcSeudo.tc)); % Replace transient times with the [scaled] SEUDO outputs
params    = se.tcSeudo.params;                                             % Output the parameter struct
if nargout > 2;    varargout{1} = se.tcSeudo.extras;    end                % If requested, output all the extra computation outputs

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
