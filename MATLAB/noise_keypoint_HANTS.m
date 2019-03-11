clear, clc
close all
%%
% matlab -nosplash -nodesktop -minimize -r main_HANTS
parameter_HANTS

%% 
indir = 'F:/whit_eval';
types = {'real', 'maxK', 'maxDer'};
percs = [0.1, 0.3, 0.5];
read_input = @(file) importdata(file, '\t', 1);

% parpool(6)
perc = 0.1;

for i = 2:length(types)
    typeI = types{i};
%     perc = percs(i);
    file_y = sprintf('%s/%s/perc_%d/TSF_whit_eval_y_noise%d.txt', indir, typeI, perc*100, perc*100);
    file_w = sprintf('%s/%s/perc_%d/TSF_whit_eval_w_noise%d.txt', indir, typeI, perc*100, perc*100);
    
    outfile_MWHA  = sprintf('fitting_MWHA_%s_noises%d.csv', typeI, perc*100);
    outfile_HANTS = sprintf('fitting_HANTS_%s_noises%d.csv', typeI, perc*100);
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
        
        [y_1, y_2, y_3, y_or]=MWHA(y, t, t, ylu, nf_WMHA, dm, HiLo, thr, fet);
        res_MWHA(j, :) = y_or';
        
        [yr, amp, phi] = HANTS(y, t', ylu, nf_HANTS, HiLo, nptperyear, fet, noutmax, delta);
        res_HANTS(j, :) = yr';
%         plot([y, y_1, y_2, y_3, y_or])
%         legend({'yi', 'y-1', 'y-2', 'y-3', 'y-or'})
    end
    csvwrite(outfile_MWHA, res_MWHA)
    csvwrite(outfile_HANTS, res_HANTS)
end                          

% test about parallel

run HANTS year by year
