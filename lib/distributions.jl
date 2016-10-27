module Distributions

#Binomial Theorem
function probOfKSuccess(n, k, pSuccess)
    nCk = Binomial(n,k)

    probability = nCk * pSuccess^n * (1 - pSuccess) ^ k

    return probability

end

function probOfKorMoreSuccess(n, k, pSuccess)

    totalProb = 0
    for i in [k:n]
        totalProb += proOfSuccess(n, i, psuccess)
    end

end


end
