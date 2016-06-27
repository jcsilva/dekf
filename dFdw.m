function [d] = dFdw(x, w0, n_inp, h_hdn, n_out)
  perturb = 1e-8;
  d=zeros(n_inp,length(w0));

  y_no_perturb = nn(x, w0, n_inp, h_hdn, n_out);
  for i=1:length(w0)
    w_pert = w0;
    w_pert(i) = w0(i) + perturb;
    y_pert = nn(x, w_pert, n_inp, h_hdn, n_out);
    d(1, i) = (y_pert - y_no_perturb) / perturb;
  end