function kldiv = build_KL_divergence(m,Kmats,q_mu,q_sigma);

kldiv = 0;
for kk = 1:m.dx
    kldiv = kldiv + KL_divergence(Kmats.Kzzi{kk},q_mu{kk},q_sigma{kk});
end
