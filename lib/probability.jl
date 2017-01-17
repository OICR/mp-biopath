module Probability

#Binomial Theorem
function KSuccess(n, k, pSuccess)
    possibilities = Binomial(n,k)
    nCk =  n == k? 1 : possibilities
    probability = nCk * pSuccess ^ k * (1 - pSuccess) ^ (n-k)

    return probability
end

function KorMoreSuccess(n, k, pSuccess)
    totalProb = Float64(0)
    for i in [k:n;]
        prob = KSuccess(n, i, pSuccess)
        totalProb = totalProb + prob
    end

    return totalProb
end

function Binomial(n,k)
    return factorial(BigInt(n))/(factorial(BigInt(n-k)) * factorial(BigInt(k)))
end

end
