function [mu_h,var_h,mu_k,var_k] = predict_MultiOutputGP_fromSpikes(KmatsSpike,Kzzi,q_mu,q_sigma,C,b,trEval,index)

[mu_k,var_k] = predict_posteriorGP_fromSpikes(KmatsSpike,Kzzi,q_mu,q_sigma,trEval);

for nn = 1:length(trEval)
    mu_h{nn} = bsxfun(@plus,sum(bsxfun(@times,mu_k{nn},C(index{trEval(nn)},:)),2), b(index{trEval(nn)}));
    var_h{nn} = sum(bsxfun(@times,var_k{nn},C(index{trEval(nn)},:).^2),2);
end
