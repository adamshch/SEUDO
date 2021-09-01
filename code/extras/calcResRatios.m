function r = calcResRatios(profs, T, cellID)
    
% r = calcResRatios(profs, T, cellID)
%
% Calculate for a given cell the residuals 
%
% 2018 - Adam Charles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

x_prof  = profs(T.window(1):T.window(2), T.window(3):T.window(4), cellID); % get profile for this cell
x_other = profs(T.window(1):T.window(2),...                      
                       T.window(3):T.window(4),[1:cellID-1,cellID+1:end]); % also get the profiles for other cells (but keep them separate)

switch 1
    case 1                                                                 % cut to the smallest bounding rectangle
        winSubY = any(x_prof>0,2);
        winSubX = any(x_prof>0,1);

    case 2                                                                 % use all pixels in the shape calculation
        winSubY = true(size(x_prof,1),1);
        winSubX = true(size(x_prof,2),1);
end

% keep only the pixels in the window
x_prof  = x_prof(winSubY,winSubX);                                         % Restrict the cell profile to the given window...
Ffits   = T.shapes(winSubY,winSubX,:);                                     %    ... as well as the data...
x_other = x_other(winSubY,winSubX,:);                                      %    ... as well as the other cell profiles

Ffits        = imgaussfilt(Ffits,1);                                       % Low-pass filter Fluorescence data with a Gaussian filter

% note size of final window
wY = sum(winSubY);                                                         % 
wX = sum(winSubX);                                                         % 

% reshape
x_prof  = x_prof(:);
Ffits   = reshape(Ffits,wX*wY,[]);
x_other = reshape(x_other,wX*wY,[]);
x_other = x_other(:,sum(x_other)~=0);

% keepPix = x_prof~=0;
% x_prof  = x_prof(keepPix,:);
% Ffits   = Ffits(keepPix,:);
% x_other = x_other(keepPix,:);

fit_vals = [x_prof,x_other]\Ffits;                                         % Calculate the fits to all ROIs in the image
% X = [x_prof,x_other];
% fit_vals = zeros(size(X,2), size(Ffits,2));
% for ll = 1:size(Ffits,2)
%     fit_vals = quadprog(2*(X'*X), -2*X'*Ffits(:,ll), [], [], [], [], ...
%                       zeros(size(X,2),1),[],[],optimset('Display','off'));
% end
res_evt  = Ffits - x_other*fit_vals(2:end,:);                              % Calculate event-triggering profile: residual images removing all but ROI in question
res      = Ffits - [x_prof,x_other]*fit_vals;                              % Calculate the residual images - based on all ROIs

r.resNorms  = vec(sum(res.^2,1));                                          % Get residual norms using all profiles
r.resNoCell = vec(sum(res_evt.^2,1));                                      % Get residual norms using only other profiles 
r.resRatios = vec(sum(res.^2,1))./vec(sum(res_evt.^2,1));                  % Get residual ratios

r.fits = Ffits;
r.zeros.numA = vec(sum(fit_vals<0,1));
r.zeros.numP = vec(sum(fit_vals(1,:)<0,1));
r.zeros.numO = vec(sum(fit_vals<0,1));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%