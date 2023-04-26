function P = assemble_space_projector(N,I)

P = repmat({I},1,N);
P = blkdiag(P{:});
