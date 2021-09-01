function op = prox_l1pos( q , N ) 
%PROX_L1POS    L1 norm, restricted to x >= 0
%    OP = PROX_L1( q ) implements the nonsmooth function
%        OP(X) = norm(q.*X,1) + indicator_{ X >= 0 }
%    Q is optional; if omitted, Q=1 is assumed. But if Q is supplied,
%    then it must be a positive real scalar (or must be same size as X).

% New in v1.0d

if nargin == 0,
        q = 1;
elseif ~isnumeric( q ) || ~isreal( q ) ||  any( q < 0 ) || all(q==0) %|| numel( q ) ~= 1
    error( 'Argument must be positive.' );
end

op = tfocs_prox( @(x) f(x,q,N), @(x,t) prox_f(x,t,q), 'vector');       

end

function v = f(x,q)
    if any( x((N+1):end) < 0 ) 
        v = Inf;
    elseif isscalar(q)
        v = q*sum( x(:) );
    %elseif numel(q) = numel(x((N+1):end))
    %    v = sum( q(:).*x((N+1):end));
    elseif numel(q) == numel(x)
        v = sum( q(:).*x(:) );
    end 
end

% The proximity operator is a simplified version of shrinkage:

%function x = prox_f(x,t,q)  
%    x   = max(0, x - t*q);
%end

function x = prox_f(q,x,t)
    tq = t .* q; % March 2012, allowing vectorized stepsizes
    s  = 1 - min( tq./abs(x), 1 );
    x  = x .* s;
end

