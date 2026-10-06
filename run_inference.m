%% run_inference.m
close all
clear

% A step-by-step comprehensive script for inferring population spiking rate from
% wide-field fluorescence using any of the five methods in this repository.

% This script is designed to be flexible, clear, and simple to use with your
% data.
% To ease initial engagement, this script is pre-configured to run with a
% real dataset containing multiple trials from a single pixel of the dorsal
% cortex (GCaMP6s). The dataset, Musall_G6sData.mat pixel 1 (all trials), is
% included in the repository.
% Every location that requires changes in information to fit specifically
% your data is marked with: <<< UPDATE HERE.

% Script Structure:
% 1) Header  -- load your data, specify the calcium decay rate (gamma),
% choose the inference method, and set its parameter or a range for parameter search.
% 2) Truncate data if needed - split long traces into overlapping chunks (chunk_trace.m)
% 3) Find best parameter  -- identify the optimal method parameter if you
% do not already know it via odd/even cross-validation (search_best_param_oddeven.m).
% 4) Infer  -- run inference (run_method.m)
% 5) Re-build the full traces  -- stitch chunks back together if needed (stitch_chunks.m)
% 6) Example plots  -- display an example trial or pixel, and the
% trial-averaged result if applicable

%% ============================== HEADER ==============================

% 1.1) YOUR DATA GOES HERE
%    Load a Time x Pixels/Trials matrix of delta-F/F fluorescence and
%    name it dffed_fluor.
load('Musall_G6sData.mat');                             % <<< UPDATE HERE: your data file
dffed_fluor = squeeze(G6sData(1, 1:204, :));            % <<< UPDATE HERE: your T x P matrix, named dffed_fluor
n_timepoints = size(dffed_fluor,1);                     % useful to have
n_traces = size(dffed_fluor,2);                         
fprintf('data includes %d pixels/trials, each with %d sample time points.\n', n_traces, n_timepoints);
recording_rate = 30; %hz                                % <<< UPDATE HERE: your recording rate

% rescale for faster inference, fluorescence units are arbitrary 
dffed_fluor  = 10* dffed_fluor / mean(abs(dffed_fluor (:)));

% 1.2) YOUR GAMMA GOES HERE
%    The fraction your calcium indicator decays per time bin. 
%    If you need to calculate your calcium decay, use calcium_decay_finder.m 
%    The example below uses GCaMP6s's known decay for this dataset's 30 Hz sampling rate.
gamma = 0.983;                                          % <<< UPDATE HERE: your gamma

% 3) CHOOSE A METHOD  (pick one number)
%    1 = Continuously-Varying (Convar)   -- recommended default
%    2 = Dynamically-Binning  (Dynbin)
%    3 = First-Differences    (Firdif)
%    4 = Wiener-Filter        (Wiener)
%    5 = Lucy-Richardson      (Lucy)
method = 5;                                             % <<< UPDATE HERE: 1, 2, 3, 4, or 5

% 4) METHOD PARAMETER
%    If you already know the right value for your data, set
%    param_known = true and put it in param_value. Otherwise leave
%    param_known = false and Part 2 below will search for the best
%    value automatically, over param_search_range (a sensible default
%    range is filled in per method below if you leave this empty).
param_known = false;                                    % <<< UPDATE HERE: true or false
param_value = [];                                       % <<< UPDATE HERE, only used if param_known = true
param_search_range = [];                                % <<< UPDATE HERE to override the default search range

switch method
    case 1  % Continuously-Varying
        param_name = 'lambda';
        default_search_range = logspace(-4, 3, 50);
    case 2  % Dynamically-Binning
        param_name = 'lambda';
        default_search_range = logspace(-4, 3, 50);
    case 3  % First-Differences
        param_name = 'smt';
        default_search_range = unique(round(logspace(0, 2, 50)));
    case 4  % Wiener-Filter
        param_name = 'k';
        default_search_range = logspace(-3, 4, 50);
    case 5  % Lucy-Richardson
        param_name = 'iter';
        default_search_range = 1:1:29;
        param_known = true;
        param_value = 10;   % well-supported default number of iterations, if you just want to use it
        
    otherwise
        error('run_inference:unknownMethod', 'method must be 1, 2, 3, 4, or 5.');
end
if isempty(param_search_range)
    param_search_range = default_search_range;
end

% =======================================================================

%% ---- Part 1: split into chunks if the trace is long ----

% Deconvolution time grows quickly with trace length because it involves inverting a matrix whose size depends on the number of time points.
% Hence, traces longer than ~600 samples are split into overlapping chunks here and are stitched back together after inference (Part 4). Traces shorter than chunk_size are left as is. 
% Most trial-structured datasets (like the example data here) would pass through this step unchanged, 
% with data_size_chunked as the original data size and chunk_starts = 1.
% When chunking does happen, chunk_starts lists where each chunk begins in the
% original traces (e.g. [1 451 901]), the same for every pixel/trial.

chunk_size = 600;                                       % <<< UPDATE HERE if you want a different chunk length. A reasonable upper size limit for quick deconvolution is around 2500. 
[dffed_fluor_chunked, chunk_starts] = chunk_trace(dffed_fluor, chunk_size);
data_size_chuncked = size(dffed_fluor_chunked);    
fprintf('%d chunk(s) per trace, of up to %d samples each.\n', numel(chunk_starts), chunk_size);

%% ---- Part 2: find the best parameter, if not already known ----
if param_known
    best_param = param_value;
    fprintf('Using known %s = %g.\n', param_name, best_param);
else
    search_maxp = min([n_traces 100]); % maximum number of trials to use for the parameter search. 
    [best_param,~,~] = search_best_param_oddeven(dffed_fluor_chunked(:,randperm(n_traces,search_maxp)), method, gamma, param_search_range);
    fprintf('Best %s found: %g\n', param_name, best_param);
end

%% ---- Part 3: infer ----

% The following code can also be found inside run_method(y, gamma, method, param). 
% You can replace it with the following single-line code:
% [r, r_idx_start, r0, beta0] = run_method(dffed_fluor_chunked, gamma, method, best_param);
% It is included here explicitly to demonstrate how to use the different methods.


% for inferred calcium if used later
beta0 = []; r0 = [];

switch method
    case 1  % Continuously-Varying
        [r, r0, beta0] = convar(dffed_fluor_chunked, gamma, best_param);
        r_idx_start = 2;

    case 2  % Dynamically-Binning
        [r, r0, beta0, iter] = dynbin_wstop(dffed_fluor_chunked, gamma, best_param, 0.0001);
        r_idx_start = 2;

    case 3  % First-Differences
        r = firdif(dffed_fluor_chunked, gamma, best_param);
        r_idx_start = 2;
        % for inferred calcium if used later 
        r0 = dffed_fluor_chunked(1,:);

    case 4  % Wiener-Filter
        r = fft_wiener(dffed_fluor_chunked, gamma, best_param);
        r_idx_start = 1;

    case 5  % Lucy-Richardson
        % p_num sets the length of the calcium-decay filter kernel. 
        % It is restricted by the data size, up to half of its length. 
        % It is set here long enough for gamma^p_num to decay below 0.1%, 
        % or third of the trial timepoints, whichever is shorter.
        p_num = min([ceil(log(0.001) / log(gamma)) floor(data_size_chuncked(1)/3)]);
        r = lucric(dffed_fluor_chunked, gamma, p_num, best_param);
        r_idx_start = 1;

    otherwise
        error('run_method:unknownMethod', 'method must be an integer 1-5.');
end

%% ---- Part 4: stitch chunks back together, if needed ----

% This step concatenates r traces to retrieve r_full with the original data
% structure.
% If the original traces in the data were short and not chunked, 
% this step would also not change the result and would retrieve r_full = r;

r_full = stitch_chunks(r, r_idx_start, chunk_starts);

clear dffed_fluor_chunked r
%% ---- Part 5: plot ----

% plot all traces and their mean 
dt = 1/recording_rate;  
t = (1:n_timepoints) * dt;
figure;
subplot(2,1,1);
title('full trials/pixels and mean-subtracted average')
for i_trace = 1:n_traces
    plot(t, dffed_fluor(:, i_trace),'LineWidth',0.5);
    hold on;
end
plot(t, mean(dffed_fluor-repmat(mean(dffed_fluor,1),n_timepoints,1), 2),'k-','LineWidth',2);
ylabel('fluorescence')
subplot(2,1,2);
title('Average across all trials/pixels');
for i_trace = 1:n_traces
    plot(t(r_idx_start:end), r_full(:, i_trace),'LineWidth',0.5);
    hold on;
end
plot(t(r_idx_start:end), mean(r_full, 2),'k-','LineWidth',2);
ylabel('inferred spiking rate')
xlabel('time (s)');


% plot an example
timepoints_show = min(n_timepoints, 200);                        % <<< UPDATE HERE to choose a time span window to show
t = (1:timepoints_show) * dt;
i_example = 19;                                                  % <<< UPDATE HERE to choose which trial/pixel to show individually

figure;
subplot(2,1,1);
plot(t(1:timepoints_show), dffed_fluor(1:timepoints_show, i_example));
legend('fluorescence')
hold on;
subplot(2,1,2);
plot(t(r_idx_start:timepoints_show), r_full(1:timepoints_show-r_idx_start+1, i_example));
legend('inferred rate');
xlabel('time (s)');
title(sprintf('Example trial/pixel %d', i_example));

% if the original data set was not chuncked for inference it is possible to
% build the inferred calcium for comparison 
% The following code can also be found inside 
% rebuild_calcium(r, r0, r_idx_start, shift, gamma) with parameters set according to the method
% You can replace it with the following two-line code:
% c = rebuild_calcium(r_full, r0, r_idx_start, beta0, gamma);
% c_example = c(1:timepoints_show,i_example);
% It is included here explicitly to demonstrate how to use the different methods.


if isscalar(chunk_starts)
    Dinv = zeros(timepoints_show);
    insert_vec = 1;
    for i_t = 1:length(Dinv)
        Dinv(i_t,1:i_t) = insert_vec;
        insert_vec = [gamma^i_t, insert_vec];
    end
    
    switch method
        case 1  % Continuously-Varying
            c_example = Dinv*[r0(i_example);r_full(1:timepoints_show-1,i_example)]+beta0(i_example);
        case 2  % Dynamically-Binning
            c_example = Dinv*[r0(i_example);r_full(1:timepoints_show-1,i_example)]+beta0(i_example);
        case 3  % First-Differences
            c_example = Dinv*[dffed_fluor(1,i_example);r_full(1:timepoints_show-1,i_example)];
        case 4  % Wiener-Filter
            c_example = Dinv*r_full(1:timepoints_show,i_example);
        case 5  % Lucy-Richardson
            c_example = Dinv*r_full(1:timepoints_show,i_example);
        otherwise
            error('run_method:unknownMethod', 'method must be an integer 1-5.');
    end
   
    subplot(2,1,1);
    plot(t(1:timepoints_show), c_example);
    legend('fluorescence', 'inferred calcium')
end
