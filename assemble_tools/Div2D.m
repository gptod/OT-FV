function [div] = Div2D(ncells,nsig_in,ind,sigma)

% divergence matrix

div = sparse(ncells,nsig_in);

for sig=1:nsig_in
    K = sigma(ind.internal(sig),1);
    L = sigma(ind.internal(sig),2);
    div(K,sig) = sigma(ind.internal(sig),4);
    div(L,sig) = -sigma(ind.internal(sig),4);
end