##################################################################
#                      TOI-451-2gp-rvs-ngts
#                    Oscar Barragan, March 2026
#
# This input file reproduces the pyaneti 2D model of RV + NGTS 
# model of TOI-451 as in Barragan et al., 2026, MNRAS, 556, ag087
# If you use this file to create the setup of your analysis
# please cite the paper.
##################################################################

# -------------------------------------------------------------------------
# IMPORTANT: Special formatting of the "RV" input file for multi-dimensional GP
# For contemporaneous time series (i.e. they do not have the exact time stamps)
# -------------------------------------------------------------------------
#
# In this setup we are using a multi-dimensional Gaussian Process (MQX).
# Therefore, the fname_rv file MUST follow a 6-column structure:
#
#   time   value   value_error   offset_label   jitter_label   dimension_label
#
# where:
#
# time:
#   Time of observation (typically BJD or BJD_TDB, in days).
#
# value:
#   Measured value of the observable (e.g., RV in km/s, FWHM in km/s, BIS in km/s, photometry, etc.).
#
# value_error:
#   Uncertainty associated with "value", in the same units.
#
# offset_label:
#   Identifier used to define independent zero-point offsets (v0 terms).
#   Measurements sharing the same offset_label will share the same systemic offset.
#   This is typically used to separate instruments (e.g., HARPS, ESPRESSO, etc.).
#
# jitter_label:
#   Identifier used to define independent jitter terms.
#   Measurements sharing the same jitter_label will share the same additional
#   white-noise term (sigma_jitter).
#   In many cases this matches the offset_label, but it can be different if needed.
#
# dimension_label:
#   Identifier that defines which GP dimension the data point belongs to.
#   This tells pyaneti which time-series (e.g., RV, FWHM, BIS, photometry)
#   the point corresponds to within the multi-dimensional GP framework.
#
# IMPORTANT:
#   - Points with the same dimension_label are treated as belonging to the same
#     GP component (same row/column of the covariance matrix).
#   - The number of UNIQUE dimension_label entries must match the MQX kernel
#     definition (e.g., MQ2 -> 2 dimensions, MQ3 -> 3 dimensions, etc.).
#
# This structure allows full flexibility:
#   - Multiple instruments per dimension (via offset_label)
#   - Independent jitter terms (via jitter_label)
#   - Multi-dimensional GP correlations (via dimension_label)
#
# Make sure that the labels used here are consistent with the "telescopes"
# vector defined below in this input file.
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
# I.  Data files & general model selection
# -------------------------------------------------------------------------
# RV / timeseries file(s). This should contain all spectroscopic time-series
# and activity indicators following the logic described before in this file
fname_rv = ['mgp_all_toi451-final.dat']

# Number of planets to model (scalar). 3 planets.
nplanets = 3

# Which model pieces to fit for each planet:
# - fit_rv: fit radial-velocity model for each planet (True/False per planet)
# - fit_tr: fit transit model for each planet (True/False per planet)
# Here: we fit RV for all three planets but DO NOT fit transits (transits disabled).
fit_rv = [True, True, True]
fit_tr = [False, False, False]

# -------------------------------------------------------------------------
# II. MCMC / run controls
# -------------------------------------------------------------------------
# Thin factor for chains: only every `thin_factor` sample is kept.
thin_factor = 1

# niter: number of iterations taken into account (the burn-in total = thin_factor * niter).
niter = 500

# Number of independent Markov chains for the ensemble sampler.
nchains = 100

# Method selection:
# - 'mcmc' -> run the MCMC fit
# - 'plot'  -> create plots assuming a previous run exists
method = 'mcmc'
#method = 'plot'

# Plotting / run behaviour flags
is_plot_correlations = True   # produce correlation plots
is_plot_chains = True         # produce trace plots
is_plot_posteriors = True     # produce posterior distributions


# -------------------------------------------------------------------------
# IV. Stellar parameters (used by pyaneti for gaussian priors / derived params)
# -------------------------------------------------------------------------
# Stellar mass and radius priors (mean, sigma). Values from Barragan et al., 2026
mstar_mean  = 0.93
mstar_sigma = 0.04
rstar_mean  = 0.85
rstar_sigma = 0.03

# Stellar effective temperature (K) and uncertainty.
tstar_mean  = 5490
tstar_sigma = 115

# J-band magnitude (used to compute TSM following Kempton+2018).
mag_j = 9.636

# Labels for plotting y-axis on spectroscopic panels (example strings).
rv_ylabel = ['RV (m/s)', 'NGTS (ppt)']

# Output units for derived planet masses: 'earth', 'jupiter' or 'solar'.
unit_mass = 'earth'

# -------------------------------------------------------------------------
# V. Priors: flags and interpretation
# -------------------------------------------------------------------------
# Definition of fit flag meanings (used throughout `fit_*` arrays):
#  - 'f' : fixed parameter (value taken from corresponding min_* entry)
#  - 'u' : uniform prior between min_* and max_*
#  - 'g' : gaussian prior with mean = min_* and sigma = max_*
#  - 'b' : beta distribution with shape parameters Beta(min_*, max_*)
#
# NOTE: The code below follows the original file's choices (Table 3 of paper).

# Which parameters to sample / fix for each planet (lists indexed by planet).
fit_t0 = ['g', 'g', 'g']   # time of minimal conjuction; gaussian priors (mean, sigma)
fit_P  = ['g', 'g', 'g']   # orbital period; gaussian priors (mean, sigma)
fit_e  = ['f', 'f', 'f']   # eccentricity fixed (works only if is_ew False)
fit_w  = ['f', 'f', 'f']   # argument of periastron fixed
fit_b  = ['f', 'f', 'f']   # impact parameter; uniform priors
fit_a  = ['f', 'f', 'f']   # scaled semi-major axis or stellar density prior (explained below)
fit_rp = ['f', 'f', 'f']   # scaled radius Rp/Rs; uniform priors
fit_k  = ['u', 'u', 'u']   # Doppler semi-amplitude; uniform priors
fit_v0 = 'u'               # systemic velocities / offsets sampled with uniform prior


# -------------------------------------------------------------------------
# VI. Priors: numeric values per parameter
# -------------------------------------------------------------------------

# Prior ranges for all planets follow the indexing: element i -> planet i.
# IMPORTANT: For gaussian priors ('g'), `min_X` is the mean and `max_X` the sigma.

# Time of mid-transit (days)
# NOTE: original file appears to have an unusual pair: min_t0 contains three large numbers,
# while max_t0 contains three small numbers. We PRESERVE the original values (no logic change),
# but flag below that this looks like the original file's direct copy/format and will be used
# by pyaneti according to the flag interpretation above (here 'g' means mean=sig).
min_t0  = [10312.4421585, 10314.6581586, 10314.9661451]
max_t0  = [0.0084355, 0.0026255, 0.0033971]

# Orbital period (days) -- note same caution as for t0: min_P / max_P correspond to
# gaussian mean and sigma because fit_P = 'g' above.
min_P   = [1.8587022, 9.1925799, 16.3649323]
max_P   = [1e-5, 0.0000148, 0.0000340]

# Impact parameter b range (0..1)
min_b   = [0.0] * nplanets
max_b   = [1.0] * nplanets

# Doppler semi-amplitude K (units consistent with pyaneti setup; here small values assumed)
min_k   = [0.0] * nplanets
max_k   = [0.05] * nplanets

# Scaled planetary radii (Rp/Rs)
min_rp  = [0.0] * nplanets
max_rp  = [0.1] * nplanets

#Semi-major xis
min_a  = [1.1] * nplanets
max_a  = [10.1] * nplanets

# Eccentricity priors (unused because fit_e='f' for all planets, but declared for completeness)
min_e   = [0, 0, 0]
max_e   = [1, 1, 1]

# Argument of periastron (radians) ranges
min_w = [0.0] * nplanets
max_w = [2 * np.pi] * nplanets

# -------------------------------------------------------------------------
# VIII. Jitter terms
# -------------------------------------------------------------------------
# Choose whether to sample a jitter term for RV and transit time-series.
is_jitter_rv = True   # sample an additional RV jitter term
is_jitter_tr = False  # do not sample transit jitter (transits not modelled here)

# -------------------------------------------------------------------------
# IX. Multi-dimensional GP configuration (multi-GP)
# -------------------------------------------------------------------------
# Details for the multidimensional GP component applied to the spectroscopic time-series.
# The timeseries input file must label each time-series as if it were an 'instrument'
# in the fourth column (this allows activity indicators to be treated like extra instruments).

# The name tokens that correspond to the instruments/in-dataset labels.
telescopes = ['rv_serval', 'ngts']    # final assignment in the original file
telescopes_labels = telescopes        # set plotting labels to match final tokens

# Kernel selection for the multi-GP:
# - 'MQ1', 'MQ2', 'MQ3', ... indicate the quasi-periodic multi-dimensional GP
#   with X time-series (MQ2 => 2D GP, MQ3 => 3D GP, etc.). Here final selection is MQ2.
kernel_rv = 'MQ2'

# fit_krv: list of flags controlling sampling for the multi-GP amplitudes + kernel hyperparameters.
# The ordering convention used by pyaneti: [A0, A1, A2, A3, A4, A5, ..., kernel_hyperparams...]
fit_krv = ['f'] * 7
fit_krv[0] = 'u'   # allow sampling for first amplitude
fit_krv[1] = 'u'   # second amplitude
fit_krv[2] = 'u'   # third amplitude
fit_krv[3] = 'f'   # fixed (example: force A_3 = 0)
# We then set remaining amplitude flags to 'u' (sampled)
fit_krv[4] = 'u'
fit_krv[5] = 'u'
fit_krv[6] = 'u'

# Numeric prior ranges for the A_i amplitudes and the kernel hyper-parameters.
# These follow the fit_krv order: the first elements correspond to amplitudes ranges,
# then the kernel hyper-parameters are appended.
krv_priors = [
    -0.1, 0.1,    # A0 range
    -0.5, 0.5,    # A1 range
    -0.5, 0.0,    # A2 range
    0.0, 0.5      # A3 range (here fixed so range is effectively unused)
]

# Quasi-periodic kernel hyper-parameter priors:
# Convention used in examples:
#  [lambda_e_min, lambda_e_max, lambda_p_min, lambda_p_max, Pgp_min, Pgp_max]
QP_priors = [5., 100.0, 0.1, 2., 4.9, 5.3]

# Concatenate amplitude priors and QP hyper-parameter priors.
krv_priors = np.concatenate([krv_priors, QP_priors])

##################################################################
# END
##################################################################
