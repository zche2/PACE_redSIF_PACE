# α mapping for surrogate_meas SVD forward model: T₂ = T₁^α_coeff
const ALPHA_COEFF_MIN = 1.0
const ALPHA_COEFF_MAX = 11.0

"""Map unconstrained α_raw to physical α_coeff ∈ [ALPHA_COEFF_MIN, ALPHA_COEFF_MAX]."""
alpha_coeff_from_raw(α_raw::Real) =
    (ALPHA_COEFF_MAX - ALPHA_COEFF_MIN) / (1.0 + exp(-α_raw)) + ALPHA_COEFF_MIN
