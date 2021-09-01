function curFrame = getFrameMC(nf, keepY,keepX, meanIm, zeroLevel, allMatFiles, fileFrameCnt)

% curFrame = getFrameMC(nf, keepY,keepX, meanIm, zeroLevel, allMatFiles, fileFrameCnt)
% 
% Function to load one frame from a set of files
% 
% 2018 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get a single image frame

meanIm    = meanIm(keepY,keepX);                                           % Crop mean image
skip_blks = nf>cumsum(fileFrameCnt);
currFile  = sum(skip_blks)+1;                                              % Find which file to load from
inFrame   = nf - sum(fileFrameCnt(skip_blks));                             % Figure out which frame inside that file to load
curFrame  = allMatFiles{currFile}.movie(:,:,inFrame);                      % Load image
curFrame(isnan(curFrame)) = zeroLevel;                                     % Change the current frame
curFrame  = (curFrame(keepY,keepX) - meanIm)./(meanIm - zeroLevel);        % Normalize the frame appropriately

curFrame(isnan(curFrame)) = 0;
curFrame(isinf(curFrame)) = 0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% curFrame  = (curFrame(keepY,keepX) - meanIm)./(meanIm - zeroLevel + ...
%                                                      (meanIm==zeroLevel)); % Normalize the frame appropriately