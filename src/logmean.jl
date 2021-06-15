"""
    logmean(aL,aR)

Logarithmic mean (aL-aR)/(log(aL)-log(aR)) for two positive states aL,aR
"""

@inline function logmean(aL,aR)
    return logmean(aL,aR,log(aL),log(aR))
end

"""
    logmean(aL,aR,logL,logR)

Logarithmic mean (aL-aR)/(log(aL)-log(aR)) for two positive states aL,aR.
Uses pre-computed log values logL = log(aL),logR = log(aR)
"""
@inline function logmean(aL,aR,logL,logR)

    # From "Entropy stable num. approx. for the isothermal and polytropic Euler"
    # by Winters et al 2019.
    da = aR-aL
    aavg = .5*(aR+aL)
    f = da/aavg
    v = f*f
    if abs(f)<1e-4
        # numbers assume γ = 1.4 (Winters provides formulas for general γ)
        return aavg*(1 + v*(-.2-v*(.0512 - v*0.026038857142857)))
    else
        return -da/(logL-logR)
    end
end
