function [Fvector] = F(x, w, n_ipt_neurs, n_hdn_neurs, n_opt_neurs)
  [nn_out, hdn_out] = nn(x,w, n_ipt_neurs, n_hdn_neurs, n_opt_neurs);
  Fvector = zeros(size(x));
  Fvector(1) = nn_out;
  for i=2:length(Fvector)
    Fvector(i) = x(i-1);
  end