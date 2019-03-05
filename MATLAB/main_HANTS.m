clear, clc
close all
%%
% matlab -nosplash -nodesktop -minimize -r main_HANTS

%% PARAMETERS
ylu        = [-0.1, 1.0];
HiLo = -1;

% HANTS
nf_HANTS   = 4; % The number of frequency
nptperyear = 23;

noutmax    = 9/nptperyear*nptperyear; % DOD=5
delta      = 0.5;  % damping factor, 0.5
fet        = 0.02; % fit error tolerance, 0.05 (Zhou 2015)

% WMHA
nf_WMHA    = 1;
dm         = 3;   % D = (2m+1)*16 = 100,
thr        = 0.1; % MWHA, threshold for low value

t = (1:nptperyear);

%% 

file_y = 'TSM_CANS6_y.txt';
file_w = 'TSM_CANS6_w.txt';

y = importdata(file_y)';
w = importdata(file_w)';

n    = length(y);
x    = (1:n);

[y_1, y_2, y_3, y_or]=MWHA(y, x, x, ylu, nf_WMHA, dm, HiLo, thr, fet);

figure
subplot(211)
plot([y, y_1, y_2, y_3, y_or])
legend({'yi', 'y-1', 'y-2', 'y-3', 'y-or'})
ys = table(y, y_1, y_2, y_3);

subplot(212)
[yr, amp, phi] = HANTS(y, x', ylu, nf_HANTS, HiLo, nptperyear, fet, noutmax, delta);
plot([y, yr]);
legend({'y', 'yr'})

%% 
indir = 'F:/whit_eval';
percs = [0.1, 0.3, 0.5];
read_input = @(file) importdata(file, '\t', 1);

% parpool(6)

for i = 1:length(percs)
    perc = percs(i);
    file_y = sprintf('%s/perc_%d/TSF_whit_eval_y_noise%d.txt', indir, perc*100, perc*100);
    file_w = sprintf('%s/perc_%d/TSF_whit_eval_w_noise%d.txt', indir, perc*100, perc*100);
    
    outfile_MWHA = sprintf('fitting_MWHA_noises%d.csv', perc*100);
    outfile_HANTS = sprintf('fitting_HANTS_noises%d.csv', perc*100);
    disp(outfile_MWHA)
    
    mat_y = read_input(file_y); mat_y = mat_y.data';
    mat_w = read_input(file_w); mat_w = mat_w.data';
    
    [NTIME, NGRID] = size(mat_y); 
    res_MWHA  = nan(NGRID, NTIME);
    res_HANTS = nan(NGRID, NTIME);
    
    parfor j = 1:NGRID%/1e5
        if mod(j, 10000) == 0, fprintf('running %.1f%%: j = %d\n', j/NGRID*100, j); end
        y = mat_y(:, j); % rowvec
        w = mat_w(:, j); % rowvec
        
%         [y_1, y_2, y_3, y_or]=MWHA(y, t, t, ylu, nf_WMHA, dm, HiLo, thr, fet);
%         res_MWHA(j, :) = y_or';
        
        [yr, amp, phi] = HANTS(y, t', ylu, nf_HANTS, HiLo, nptperyear, fet, noutmax, delta);
        res_HANTS(j, :) = yr';
%         plot([y, y_1, y_2, y_3, y_or])
%         legend({'yi', 'y-1', 'y-2', 'y-3', 'y-or'})
    end
%     csvwrite(outfile_MWHA, res_MWHA)
    csvwrite(outfile_HANTS, res_HANTS)
end                          

% test about parallel
