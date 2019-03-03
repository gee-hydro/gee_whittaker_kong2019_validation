function [yr, amp, phi] = HANTS(y, ts, HiLo, nf, ylu, nptperyear, fet, noutmax, delta)
% [yr, amp,phi] = HANTS(nptperyear,nf,y,ts,HiLo,low,high,fet,dod,delta)
% HANTS processing
% 
% NOTE: This version is tested in MATLAB V2010b. In some older version you
% might get an error on line 117. Refer to the solution provided on that
% line.
%
% Modified:
%   Apply suppression of high amplitudes for near-singular case by
%	adding a number delta to the diagonal elements of matrix A,
%	except element (1,1), because the average should not be affected
%
%	Output of reconstructed time series in array yr June 2005
%
%   Change call and input arguments to accommodate a base period length (nptperyear)
%   All frequencies from 1 (base period) until nf are included
% 
% Parameters
%
% Inputs:
%   y     = colvec, array of input sample values (e.g. NDVI values)
%   ts    = colvec, array of size ni of time sample indicators
%           (indicates virtual sample number relative to the base period);
%           numbers in array ts maybe greater than nptperyear
%           If no aux file is used (no time samples), we assume ts(i)= i,
%           where i=1, ..., ni
%   nptperyear = length of the base period, measured in virtual samples
%           (days, dekads, months, etc.)
%   nf    = number of frequencies to be considered above the zero frequency
%   HiLo  = indicating rejection of -1: low, 1: high outliers. 
%   ylu   = [low, high] of time-series y (values outside the valid range are rejeced
%           right away)
%   fet   = fit error tolerance (points deviating more than fet from curve
%           fit are rejected)
%   %%dod   = degree of overdeterminedness (iteration stops if number of
%           points reaches the minimum required for curve fitting, plus
%           dod). This is a safety measure. min n_remain = dod + nr
%   noutmax = maximum throw-away points
%   delta = small positive number (e.g. 0.1) to suppress high amplitudes
%
% Outputs:
%
% amp   = returned array of amplitudes, first element is the average of
%         the curve
% phi   = returned array of phases, first element is zero
% yr	= array holding reconstructed time series
%
% Authors:
% Dongdong Kong (2017-10-27)
% Mohammad Abouali (2011), Converted to MATLAB
%
% References
% Wout Verhoef, NLR, Remote Sensing Dept. June 1998

% coder.varsize('y', 'ts', [100, 1], [1, 0])
ni  = length(y); % total number of actual samples of the time series

mat = zeros(min(2*nf+1,ni),ni,'single');
amp = zeros(nf+1,1, 'single');
phi = zeros(nf+1,1, 'single');
yr  = zeros(ni,1  , 'single');
zr  = zeros(ni*2+1, 'single');

% ra = zeros(nf, 1, 'single');
% rb = zeros(nf, 1, 'single');
% if (Opt.FirstRun==true)
HiLo_neg = -HiLo;

nr      = min(2*nf+1, ni);
% noutmax = ni - nr - dod;
dg      = 180.0/pi;

ang = 2.*pi*(0:nptperyear-1)/nptperyear;
cs  = cos(ang);
sn  = sin(ang);
%     Opt.FirstRun=false;
% end
mat(1,:)=1.0;

% f*2*pi*[0:nptperyear-1]/nptperyear; mod replace it
for i = 1:nf
    I = 1 + mod(i * (ts - 1), nptperyear);
    mat(2*i  , :) = cs(I);
    mat(2*i+1, :) = sn(I);
end
% i=1:nf;
% for j=1:ni
%     index = 1 + mod(i*(ts(j)-1), nptperyear);
%     mat(2*i  ,j) = cs(index);
%     mat(2*i+1,j) = sn(index);
% end

% p: weights of every points
p = ones(ni,1);
if ~isempty(ylu)
    p(y < ylu(1) | y > ylu(2) | isnan(y) ) = 0;
end
nout = sum(p==0);

if (nout>noutmax)
    %     disp('Not enough data points')
    %     disp(['nout:' num2str(nout) ' ,noutmax:' num2str(noutmax)])
    %     error('nout > noutmax')
    return
end

ready    = false;
nloop    = 0;
nloopmax = ni; %ni, control d to break loops

while ((~ready)&&(nloop < nloopmax))
    nloop=nloop+1;
    za = mat*(p.*y);
    
    A = mat * diag(p) * mat'; %how to know A was amplitude
    A = A + diag(ones(nr,1))*delta;
    A(1,1) = A(1,1) - delta;
    zr = A\za;
    
    yr = mat'*zr;
    if HiLo_neg == 0
        diffVec = abs(yr - y); %20171026, add top and low outliers removing method
    else
        diffVec = HiLo_neg * (yr - y);
    end
    err = p.*diffVec;
    
    [~, rankVec] = sort(err,'ascend');
    % The above line may not be recognized on some older MATLAB versions.
    % Simply comment the above line and uncomment the line below.
    %    [tmp, rankVec]=sort(err,'ascend');
    
    maxerr = diffVec(rankVec(ni));
    ready  = (maxerr <= fet) || (nout==noutmax);
    
    if (~ready)
        i=ni;
        j=rankVec(i);
        while ( (p(j)*diffVec(j) > maxerr*0.5)&&(nout<noutmax) )
            p(j) = 0;
            nout = nout+1;
            i = i-1;
            j = rankVec(i); %fixerror, rank, 20171026
        end
    end
end

if nargout > 1
    amp(1)   = zr(1);
    phi(1)   = 0.0;
    
    % zr(ni+1) = 0.0;
    i        = (2:2:nr)';
    ifr      = (i+2)/2;
    ra       = zr(i);
    rb       = zr(i+1);
    amp(ifr) = sqrt(ra.*ra+rb.*rb);
    phase    = atan2(rb, ra)*dg;
    phase(phase<0) = phase(phase<0) + 360;
    phi(ifr) = phase;
end
end