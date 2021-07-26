import numpy as np
# from math import comb # No longer needed because we use faster Pascal method to calculate n_choose_k
from scipy.special import gamma

callCount = 0

outArray = []

# This should be set by caller, to handle computational progress
slider_progress_callback = None

# This can be set by caller, to halt ongoing computation (only works with Riemann)
quit_computation_flag = False

RIEMANN_ITER_LIMIT = 150  # Default. Can be overridden
NK2_array = []
NK1_array = []

# Precompute table of coefficients for lookup. This reduces Riemann computation time from O(n^2) to O(n)
def precompute_coeffs():
    global NK2_array, NK1_array
    NK2_array = np.zeros(RIEMANN_ITER_LIMIT)
    print("Precomputing coefficients for Riemann/Dirichlet sum using Euler's transformation...")

    # Precompute N_choose_k / 2^(N+1) coefficients.
    NK1_array = np.zeros(shape=(RIEMANN_ITER_LIMIT, RIEMANN_ITER_LIMIT))
    NK1_array[0, 0] = 0.5  # This will be (n_choose_k) / 2^(n+1)
    for n in range(1, RIEMANN_ITER_LIMIT):
        NK1_array[n, 0] = NK1_array[n - 1, 0] / 2
        for k in range(1, n + 1):
            # Pascal's triangle, but with an additional divide by 2 at each row.
            NK1_array[n, k] = (NK1_array[n - 1, k - 1] + NK1_array[n - 1, k]) / 2

    # Precompute sum of above coefficients for each value of k. These will be used in Euler transform
    for k in range(0, RIEMANN_ITER_LIMIT):
        tmp_sum = 0
        for n in range(k, RIEMANN_ITER_LIMIT):
            tmp_sum += NK1_array[n, k]  # comb(n,k) / (2 ** (n+1))
        NK2_array[k] = ((-1) ** k) * tmp_sum

    print("Done precomputing coefficients\n")


def EtaToZetaScale(v):
    # Scale factor converts Dirichlet eta function to Riemann zeta function
    return 1 / (1 - 2 ** (1 - v))


def print_status():
    global callCount
    callCount = callCount + 1
    if slider_progress_callback is not None:
        if np.mod(callCount, 100) == 0:
            slider_progress_callback(callCount)
    else:
        if (np.mod(callCount, 20000) == 0) & (slider_progress_callback is not None):
            print("   \r" + str(int(callCount / 1000)) + "k ", end='')  # Print status every 10k samples


#
# Calculate Riemann zeta function for complex input s.
#
# If second argument is true, then also return array size of intermediate sum calculation (used to draw arrows)
#
# Approximate accuracy for real s is (very) roughly proportional to magnitude of final summation
# term, which is about 1/(2^ITER * ITER^s). Hence, number of digits of precision is roughly
# log_10(2^ITER * ITER^s) = .301 * ITER + s * log_10(ITER). For example:
#
# If ITER = 100, digits of precision is roughly 30 + 2*Re(s)
# If ITER = 1000, digits of precision is roughly 300 + 3*Re(s)
#
# Adding n to the ITER count gives roughly (0.3 + Re(s)/ITER) * n additional digits of
# precision, assuming ITER >> n. Precision improves slightly as Im(s) increases from 0, I think.
#
# Note that convergence is very poor for Re(s) large and negative. So we use the functional equation
# for Re(s) < 0.
#
def Riemann(s, get_array_size=False, do_eta=False):
    global callCount

    if np.size(s) > 1:
        return [0.0 if quit_computation_flag else Riemann(x, do_eta=do_eta) for x in s]

    if s == 1.0:
        # Calculation blows up at 1.0, so return nan
        print_status()
        return np.nan

    if np.real(s) < 0:
        if do_eta:
            # Use functional equation. Don't call printstatus yet
            return -s * Riemann(1 - s, do_eta=do_eta) * gamma(- s) * np.sin(-s * np.pi / 2) * (np.pi ** (s - 1)) * (
                        1 - 2 ** (s - 1)) / (1 - 2 ** s)
        else:
            # Use functional equation. Don't call printstatus yet
            return Riemann(1 - s, do_eta=do_eta) * gamma(1 - s) * np.sin(s * np.pi / 2) * (np.pi ** (s - 1)) * (2 ** s)

    cum_sum = 0 + 0j
    # Need first element zero so line segment will draw correctly
    store_intermediates = len(outArray) > 0
    if store_intermediates:
        outArray[0] = cum_sum

    if do_eta:
        scale1 = 1
    else:
        # Scale factor converts Dirichlet eta function to Riemann zeta function
        scale1 = 1 / (1 - 2 ** (1 - s))

    if store_intermediates or get_array_size:

        partial_sum_index = 1

        # Calculate terms of Dirichlet eta function, then apply scale1 factor to get Riemann zeta
        for k in range(0, RIEMANN_ITER_LIMIT):
            cum_sum = cum_sum + scale1 * NK2_array[k] / ((k + 1) ** s)
            if k < 600 or np.mod(k, 2) == 0:
                # After first 600 points, only plot every other point, to speed up graphics.
                # Should get rid of this, now that graphics have been sped up with draw_artist
                if store_intermediates:
                    outArray[partial_sum_index] = cum_sum
                partial_sum_index = partial_sum_index + 1
    else:

        # Calculate terms of Dirichlet eta function, then apply scale1 factor to get Riemann zeta
        for k in range(0, RIEMANN_ITER_LIMIT):
            cum_sum = cum_sum + NK2_array[k] / ((k + 1) ** s)

        if not do_eta:
            # Scale factor converts Dirichlet eta function to Riemann zeta function
            cum_sum = cum_sum * scale1

    # Make sure final point is included
    if np.mod(RIEMANN_ITER_LIMIT - 1, 10) != 0:
        if store_intermediates:
            outArray[partial_sum_index] = cum_sum

        if store_intermediates or get_array_size:
            partial_sum_index = partial_sum_index + 1

    print_status()

    if get_array_size:
        return cum_sum, partial_sum_index
    else:
        return cum_sum


#
# When multiplied by gamma(s/2) * pi ^(-s/2), the result has 180-deg rotational symmetry around s = 0.5 + 0j
#
def RiemannSymmetric(s):
    return Riemann(s) * gamma(s / 2) * (np.pi ** (-s / 2))
#    return (np.pi ** (-s/2))
