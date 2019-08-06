function [X,info] = TVinpaint(B,M,delta,eps_rel)
%TVINPAINT  Total variation image inpainting
% 
% X = TVinpaint(B,M,delta)
% [X,info] = TVinpaint(B,M,delta,eps_rel)
%
% This function solves the TV inpainting problem
%
%    min  TV(X)  subject to   || X(Ic) - B(Ic) ||_F <= delta
%
% where B is a noisy image with missing pixels, Ic are the indices to
% the intact pixels, X is the reconstruction, and delta is an upper
% bound for the residual norm.  The TV function is the 1-norm
% of the gradient magnitude, computed via neighbor pixel differences.
% At the image borders, we imposed reflexive boundary conditions for
% the gradient computations.
%
% The information about the intact and missing pixels is given in the
% form of the mask M, which is a matrix of the same size as B, and whose
% nonzero elements indicate missing pixels.
%
% The parameter delta should be of the same size as the norm of the
% image noise.  If the image is m-times-n, and sigma is the standard
% deviation of the image noise in a pixel, then we recommend to use
% delta = tau*sqrt(m*n)*sigma, where tau is slightly smaller than one,
% say, tau = 0.85.
%
% The function returns an epsilon-optimal solution X, meaning that
% if X* is the exact solution, then our solution X satisfies
%
%     TV(X) - TV(X*) <= epsilon = max(B(Ic))*m*n*eps_rel,
%
% where eps_rel is the specified relative accuracy; the default is
% eps_rel = 1e-3.
% 
% The solution status is returned in the stuct info, with info.STATUS
% having one of the settings
%  'EPSILON-OPTIMAL-SOLUTION': X is an epsilon-optimal solution
%  'NOT-EPSILON-OPTIMAL-SOLUTION': X is not an epsilon-optimal solution
%  'MAXIMUM-NUMBER-OF-ITERATIONS-EXCEEDED': X is not an epsilon-optimal
%   solution when the maximum number of iterations was reached. 
% Other fields of info:
%  info.NDENOISE      Iteration bound for the problem given, as
%                     computed by the function.
%  info.ITERATIONS_K  The number of iterations used.	
%  info.EPS_REL_K     The relative accuracy reached.
%  info.TIME          Time in seconds for the program to run. 
%
% See also: TVdeblur, TVdenoise.

% J. Dahl^1, P.C. Hansen^2, S.H. Jensen^1 & T.L. Jensen^1
% CSI project, (1) Aalborg University, (2)Technical University of Denmark
% April 28, 2009.

% Check input parameters.
if nargin < 3	
	error('Too few input parameters');
elseif nargin == 3
	eps_rel = 1e-3;
end
if size(M) ~= size(B)
    error('The mask M must have the same dimensions as B');
end

tic

% Indices to intact and missing pixels.
Ic = int32( find(M(:) == 0) );
I  = int32( find(M(:) ~= 0) );

%Check for the constant solution and print out if we might have numerical problems
alpha = sum(B(Ic))/numel(Ic);
X = alpha*ones(size(B));
mdelta = norm(X-B,'fro');
if  mdelta < delta
    info = info_type_inpaint(1,0,0,0,toc);
    return;
elseif mdelta < 1.1*delta
    disp('The algorithm might experience numerical problem, because the');
    disp('solution X is almost alpha*ones(size(B)). Proceed with care');
end

% Dynamical range in input image.
R = max(B(Ic));
S = min(B(Ic));

% Set parameter gamma for dealing with the missing pixels.
gamma = 0.5*(R-S)*sqrt(length(I));
d = (R-S)/2 + S;  % Center for the prox-function.

% More parameters for the inpainting algorithm.
mn = numel(B);
epsilon = R*mn*eps_rel;
mu = epsilon/mn;
Lmu = 8/mu ;
N = int32( ceil(4*sqrt(2*mn)*sqrt(delta^2+gamma^2)/epsilon) );

% Compute TV solution via C function.
[X,k,epsilon_k] = tv_inpaint(B,I,Ic,delta,gamma,d,epsilon,Lmu,mu,N);

if nargout == 2
    
	if k >= N
	   info = info_type_inpaint(3,N,k,epsilon_k/(R*mn),toc);
    elseif epsilon_k > epsilon
	    info = info_type_inpaint(2,N,k,epsilon_k/(R*mn),toc);
	else
		info = info_type_inpaint(1,N,k,epsilon_k/(R*mn),toc);	
    end

end
function info = info_type_inpaint(type,N,k,eps_rel_k,time)
%
% Type-like definitions for TVinpaint to be used with mxTV software
%
    
    info_type = struct('STATUS',{'EPSILON-OPTIMAL-SOLUTION',...
        'NOT-EPSILON-OPTIMAL-SOLUTION',...
        'MAXIMUM-NUMBER-OF-ITERATIONS-EXCEEDED'},...
        'NINPAINT',{N},'ITERATIONS_K',{k},...
        'EPS_REL_K',{eps_rel_k},...
        'TIME',{time});

    info = info_type(type);