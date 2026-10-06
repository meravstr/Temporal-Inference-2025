function [chunked, chunk_starts] = chunk_trace(dffed_fluor, chunk_size)
% CHUNK_TRACE splits long fluorescence traces into overlapping traces (chunks) along time, 
% so each new trace is short enough for quick inference. Deconvolution cost grows quickly with trace length, so multi-thousand-sample traces should be broken into shorter pieces before inference, and stitched back together afterwards (see stitch_chunks.m).
% Consecutive chunks overlap by one quarter of chunk_size. That overlap allows stitch_chunks.m to splice the pieces back together.

% INPUTS dffed_fluor - original traces, TxN (time x pixels/trials);
% chunk_size - the maximal trace length after splitting to chunks.
% OUTPUTS chunked - fluorescence traces after splitting, t x (n_chunks * N) matrix. The columns are ordered trace by trace: all chunks of trace 1 (in time order), then all chunks of trace 2, and so on.
% chunk_starts - a 1 x n_chunks vector: the row of the original traces where each chunk begins (the same for every trace), e.g. [1 451 901].
% If T <= chunk_size, no chunking is needed: chunked = dffed_fluor and chunk_starts = 1.

[T, N] = size(dffed_fluor);

if T <= chunk_size
    chunked = dffed_fluor;
    chunk_starts = 1;
    return;
end

if chunk_size < 8
    error('chunk_trace:chunkTooSmall', 'chunk_size must be at least 8 samples.');
end

overlap = round(chunk_size / 4);
step = chunk_size - overlap;

chunk_starts = 1:step:(T - chunk_size + 1);
if chunk_starts(end) + chunk_size - 1 < T
    chunk_starts(end+1) = T - chunk_size + 1;   % make sure the last chunk reaches the end of the trace
end
n_chunks = numel(chunk_starts);

chunked = zeros(chunk_size, n_chunks * N);
for i = 1:n_chunks
    % chunk i of every trace: columns i, i + n_chunks, i + 2*n_chunks, ...
    chunked(:, i:n_chunks:end) = dffed_fluor(chunk_starts(i):chunk_starts(i) + chunk_size - 1, :);
end

end
