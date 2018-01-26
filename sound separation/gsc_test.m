mu0s = 0.03:0.0025:0.13;
mu_maxs = 0.2:0.005:0.4;

err = zeros(length(mu0s),length(mu_maxs));

mu0_num = 0;
for mu0 = mu0s
  mu0_num = mu0_num + 1
  mu_max_num = 0;
  for mu_max = mu_maxs
    mu_max_num = mu_max_num+1;
    gsc
    err(mu0_num,mu_max_num) = sqrt(sum((s2-o).^2));
  end
end

