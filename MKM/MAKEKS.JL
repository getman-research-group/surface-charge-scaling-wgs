function makeks(qreact, qts, qprod, erxn, eact, consts)
    kb   = consts.kb
    T    = consts.Te
    h    = consts.h
    kbeV = consts.kbeV

    # Preserve the original implementation's zero floor for reverse barriers.
    # The SI prints eact-erxn without this floor; the neutral TOF is sensitive
    # to this distinction for the selected nonlinear four-descriptor model.
    erev = max(eact - erxn, zero(eact - erxn))

    kf  = (kb * T / h) * (qts / qreact) * exp(-eact / (kbeV * T))
    keq = (qprod / qreact) * exp(-erxn / (kbeV * T))
    kr  = (kb * T / h) * (qts / qprod) * exp(-erev / (kbeV * T))

    return kf, kr, keq
end
