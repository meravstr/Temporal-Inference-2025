function r_full = stitch_chunks(r, r_idx_start, chunk_starts)
% STITCH_CHUNKS concatenates the inferred rate of every short trace (chunk) (produced by chunk_trace.m + inference) back into one continuous trace per original pixel/trial. 
% Since the inferred rate lacks baseline (constant shift) information, adjacent chunks overlap (see chunk_trace.m), and the overlapping parts are used to align the constant shift. 

% INPUTS  r - the inferred rate (chunk_size - r_idx_start + 1) x (n_chunks * N)
% r_idx_start - the row of the fluorescence that the first row of r corresponds to ( 1 or 2)
% chunk_starts - vector of numbers marking the rows in which the original fluorescence was cut into the next chunk (1 if no chunking).
% OUTPUT r_full - (T - r_idx_start + 1) x N, where T is the original trace length. 
% If data were never chunked (chunk_starts = 1), r_full = r.

chunk_starts = chunk_starts(:)';
n_chunks = numel(chunk_starts);

if n_chunks == 1
    r_full = r;
    return;
end

n_cols = size(r, 2);
if mod(n_cols, n_chunks) ~= 0
    error('stitch_chunks:sizeMismatch', ...
        'The number of columns of r (%d) must be a multiple of numel(chunk_starts) (%d).', n_cols, n_chunks);
end
n_traces = n_cols / n_chunks;

chunk_size = size(r, 1) + r_idx_start - 1;   % length of each fluorescence chunk
T = chunk_starts(end) + chunk_size - 1;

% last original index kept from chunk i = midpoint of its overlap with chunk i+1
split_after = zeros(1, n_chunks - 1);
for i = 1:n_chunks - 1
    actual_overlap = chunk_starts(i) + chunk_size - chunk_starts(i+1);
    split_after(i) = chunk_starts(i+1) + floor(actual_overlap / 2) - 1;
end

r_full = nan(T - r_idx_start + 1, n_traces);

for i = 1:n_chunks
    if i == 1
        keep_start = 1;
    else
        keep_start = split_after(i-1) + 1;
    end
    if i == n_chunks
        keep_end = T;
    else
        keep_end = split_after(i);
    end

    abs_idx   = chunk_starts(i) - 1 + (r_idx_start:chunk_size)';   % original time index of each row of r
    keep_mask = abs_idx >= keep_start & abs_idx <= keep_end;

    % chunk i of every trace sits in columns i, i + n_chunks, i + 2*n_chunks, ...
    r_full(abs_idx(keep_mask) - r_idx_start + 1, :) = r(keep_mask, i:n_chunks:end);
end

end
