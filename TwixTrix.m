function [score,Z,C] = TwixTrix(varargin)
% TWIXTRIX - Double 2-way t-test network inference
% TWIXTRIX implements a double 2-way t-test network inference method
% using unmoderated t-statistics
%
% INPUT:
%   - D   : data matrix, rows are genes, columns are samples
%   - tfs : list of regulators in the form of row numbers (optional)
%
% OUTPUT:
%   - score : matrix of size #tfs x #genes containing t-statistic for each
%             TF in the critical contrast of each gene
%   - Z     : Z-transformed score
%
% Copyright (C) 2012, Jianlong Qi and Tom Michoel
%
%   Jianlong Qi <jianlong.qi@frias.uni-freiburg.de>
%   Tom Michoel <tom.michoel@roslin.ed.ac.uk>
%
%   http://www.roslin.ed.ac.uk/tom-michoel
%

D = varargin{1}; % data
N = size(D,1); % number of genes
K = size(D,2); % number of samples

if nnz(isnan(D))>0
    error('This version of TwixTrix does not yet handle missing values, consider replacing them by zeros');
end

switch nargin
    case 1 % all-against-all testing
        tfs = 1:N;
    case 2
        tfs = varargin{2};
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Critical contrasts %%%
%%%%%%%%%%%%%%%%%%%%%%%%

% sort data rows
[D2,P] = sort(D,2);
% cumulative row sums
R = cumsum(D2,2);
% cumulative sum of squares
S = cumsum(D2.^2,2);
% number of elements in C1
M = 1:K-1;

% mean & std C1; note that sig = (n-1)*sigma^2
n1 = repmat(M, [N 1]);
mu1 = R(:,1:end-1)./n1;
sig1 = S(:,1:end-1) - n1.*mu1.^2;
% mean & std C2; note that sig = (n-1)*sigma^2
n2 = repmat(K-M, [N 1]);
mu2 = (repmat(R(:,end), [1 K-1]) - R(:,1:end-1))./n2;
sig2 = (repmat(S(:,end), [1 K-1]) - S(:,1:end-1)) - n2.*mu2.^2;

% all t-statistics
T = abs(mu1 - mu2)./(sqrt((sig1+sig2)./(n1+n2-2)).*sqrt(1./n1 + 1./n2));

% critical contrasts
[~,I] = max(T,[],2);
C = true(N,K);
for n=1:N
    C(n,P(n,1:I(n)))=false;
end
C = C'; % transpose for matrix multiplication later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DE of TFs in critical contrast of targets %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TF data
D3 = D(tfs,:);
L = length(tfs);

% mean & std C1; note that sig = (n-1)*sigma^2
n1 = repmat(sum(C==0), [L 1]);
mu1 = (D3*~C)./n1;
sig1 = (D3.^2*~C) - n1.*mu1.^2;
% mean & std C2; note that sig = (n-1)*sigma^2
n2 = repmat(sum(C==1), [L 1]);
mu2 = (D3*C)./n2;
sig2 = (D3.^2*C) - n2.*mu2.^2;

% all t-statistics
score = abs(mu1 - mu2)./(sqrt((sig1+sig2)./(n1+n2-2)).*sqrt(1./n1 + 1./n2));

% remove self-interactions
for k=1:length(tfs)
    score(k,tfs(k)) = 0;
end

% normalize over genes
Z = (score - repmat(mean(score), [L 1]))./repmat(std(score),[L 1]);

% return original C if requested
if nargin>3
    C = C';
end