function ind = indices(ncell,sigma)

% INPUT:
% ncell: number of cells
% edges: edges structure
% OUPUT:
% ind: indices structure

% this function creates indices for a simplified access to the matrices
% these indices enable to select the entries of submatrices of the FV
% matrices

% indices for internal edges:

internal = find(sigma(:,2)~=0);

ind.internal = internal;
ind.bound = find(sigma(:,2)==0);

ind.nsig_in = length(internal);
ind.nsig_b = length(ind.bound);

% OLD

% ind.i_KK = sub2ind([ncell ncell],edges(internal,3),edges(internal,3));
% ind.i_LL = sub2ind([ncell ncell],edges(internal,4),edges(internal,4));
% ind.i_KL = sub2ind([ncell ncell],edges(internal,3),edges(internal,4));
% ind.i_LK = sub2ind([ncell ncell],edges(internal,4),edges(internal,3));

% [ind.u_i,~,ind.o_i] = unique(edges(internal,3));
% [ind.u_e,~,ind.o_e] = unique(edges(internal,4));

% [ind.u_KK,~,ind.o_KK] = unique(ind.i_KK);
% [ind.u_LL,~,ind.o_LL] = unique(ind.i_LL);
% [ind.u_KL,~,ind.o_KL] = unique(ind.i_KL);
% [ind.u_LK,~,ind.o_LK] = unique(ind.i_LK);

% ind.et_K = sub2ind([ind.nei ncell],[1:ind.nei]',edges(internal,3));
% ind.et_L = sub2ind([ind.nei ncell],[1:ind.nei]',edges(internal,4));

%[ind.etu_K,~,ind.eto_K] = unique(ind.et_K);
%[ind.etu_L,~,ind.eto_L] = unique(ind.et_L);

% indices for all edges:

% [ind.allu_i,~,ind.allo_i] = unique(edges(:,3));
% [ind.allu_e,~,ind.allo_e] = unique(edges(internal,4));
% 
% ind.alli_ii = sub2ind([ncell ncell],edges(:,3),edges(:,3));
% ind.alli_ee = sub2ind([ncell ncell],edges(internal,4),edges(internal,4));
% ind.alli_ie = sub2ind([ncell ncell],edges(internal,3),edges(internal,4));
% ind.alli_ei = sub2ind([ncell ncell],edges(internal,4),edges(internal,3));
% 
% [ind.allu_ii,~,ind.allo_ii] = unique(ind.alli_ii);
% [ind.allu_ee,~,ind.allo_ee] = unique(ind.alli_ee);
% [ind.allu_ie,~,ind.allo_ie] = unique(ind.alli_ie);
% [ind.allu_ei,~,ind.allo_ei] = unique(ind.alli_ei);
% 
% nedge = size(edges,1);
% ind.allet_i = sub2ind([nedge ncell],[1:nedge]',edges(:,3));
% ind.allet_e = sub2ind([nedge ncell],internal,edges(internal,4));
% [ind.alletu_i,~,ind.alleto_i] = unique(ind.allet_i);
% [ind.alletu_e,~,ind.alleto_e] = unique(ind.allet_e);

end

