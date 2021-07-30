import typing
import bisect
import numpy as np
from scipy.special import gamma
import time

outArray = []

# This should be set by caller, to handle computational progress
computation_progress_callback: typing.Callable = None

# This can be set by caller, to halt ongoing computation (only works with Riemann)
quit_computation_flag = False

# True if we compute Riemann zeta using cached exponentials (gives ~10x speedup)
USE_CACHED_FUNC = True

RIEMANN_ITER_LIMIT = 40  # Default. Can be overridden
NK2_array = []
NK1_array = []

# Used to determine when to update progress bar, and also tracks which precomputed denominator coefficients to use.
row_count = 0


# Precompute table of coefficients for lookup using Euler's transformation.
# This reduces Riemann computation time from O(n^2) to O(n) where n is number of terms
def precompute_coeffs():
    global NK2_array, NK1_array
    NK2_array = np.zeros(RIEMANN_ITER_LIMIT)
    #    print("Precomputing " + str(RIEMANN_ITER_LIMIT) + " coefficients for Riemann/Dirichlet sum", end="")

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


def eta_zeta_scale(v):
    # Scale factor converts Dirichlet eta function to Riemann zeta function
    return 1 / (1 - 2 ** (1 - v))


# This tracks where Re(s) switches from negative to positive. We store this when first row is calculated, so that
# we can ensure split occurs at exactly the same place for every row. Otherwise rounding errors might cause slight
# differences
array_zero_split = 0
elapsed_time = 0
entry_count = 0


#
# Calculate Riemann zeta function for complex input s.
#
# If second argument is true, then also return array size of intermediate sum calculation (used to draw arrows)
#
# Note that s can be a 1D or 2D array, in which it will return output in the same
# format as the input

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
def riemann(s, get_array_size=False, do_eta=False, use_zero_for_nan=True):
    global row_count, array_zero_split, elapsed_time, entry_count

    if s.ndim == 2:
        if entry_count > 1:
            # Prevent this from running more than once concurrently
            return 0

        entry_count += 1

        t1 = time.time()
        precompute_denom(s)
        delay = time.time() - t1
        print("Precomputed denominators in " + "{:1.4f}".format(delay) + " seconds")
        row_count = 0
        t1 = time.time()
        # When doing heatmaps, s is initially a 2d ndarray, which behaves like a list of 1d ndarrays.
        # The following will call itself recursively for each item (row), and build a list of 1d arrays
        out = [0.0 if quit_computation_flag else riemann(x, False, do_eta, use_zero_for_nan) for x in s]
        # Convert list of 1d arrays into a single 2d ndarray
        out = np.stack(out)
        elapsed_time = time.time() - t1
        entry_count -= 1
        return out

    # Now we have 1D vector or single value
    if np.size(s) > 1:
        # We have 1D vector. Update progress bar, then call Riemann for individual values
        row_count += 1

        if np.mod(row_count, 50) == 0:
            # Update progress bar every 50 rows
            if computation_progress_callback is not None:
                computation_progress_callback(np.size(s) * 50)

        # RiemannArray requires elements to either be all negative or all non-negative. So we split s into two arrays.
        neg_array = s[0:array_zero_split]
        pos_array = s[array_zero_split:len(s)]

        #            a1 = riemann_array(neg_array, row_count - 1, do_eta=do_eta)
        #            a2 = riemann_array(pos_array, row_count - 1, do_eta=do_eta)
        return np.concatenate((riemann_array(neg_array, row_count - 1, do_eta=do_eta),
                               riemann_array(pos_array, row_count - 1, do_eta=do_eta)))

    #
    # Single value. Heatmap calculations never come here, since we used vectorized version instead.
    #
    if s == 1.0:
        # Calculation blows up at 1.0, so don't bother to calculate
        if use_zero_for_nan:
            # Do this to avoid warnings.
            return 0
        else:
            return np.nan

    if np.real(s) < 0:
        # Use functional equation
        if do_eta:
            return -s * riemann(1 - s, do_eta=do_eta) \
                   * gamma(- s) * np.sin(-s * np.pi / 2) * (np.pi ** (s - 1)) \
                   * (1 - 2 ** (s - 1)) / (1 - 2 ** s)
        else:
            return riemann(1 - s, do_eta=do_eta) \
                   * gamma(1 - s) * np.sin(s * np.pi / 2) * (np.pi ** (s - 1)) * (2 ** s)

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
        have_final_point = False

        # Calculate terms of Dirichlet eta function, then apply scale1 factor to get Riemann zeta
        for k in range(0, RIEMANN_ITER_LIMIT):
            cum_sum = cum_sum + scale1 * NK2_array[k] / ((k + 1) ** s)
            if k < 600 or np.mod(k, 2) == 0:
                # After first 600 points, only plot every other point, to speed up graphics.
                # Should get rid of this, now that graphics have been sped up with draw_artist
                if store_intermediates:
                    outArray[partial_sum_index] = cum_sum
                partial_sum_index = partial_sum_index + 1
                if k == RIEMANN_ITER_LIMIT - 1:
                    have_final_point = True
    else:

        # Calculate terms of Dirichlet eta function, then apply scale1 factor to get Riemann zeta
        for k in range(0, RIEMANN_ITER_LIMIT):
            cum_sum = cum_sum + NK2_array[k] / ((k + 1) ** s)

        if not do_eta:
            # Scale factor converts Dirichlet eta function to Riemann zeta function
            cum_sum = cum_sum * scale1

    if store_intermediates or get_array_size:

        # Make sure final point is included
        if not have_final_point:
            partial_sum_index += 1
            outArray[partial_sum_index] = cum_sum

        if get_array_size:
            return cum_sum, partial_sum_index

    return cum_sum


pre_computed_denom_right_mag: np.ndarray  # (k+1)^(Re(-s))
pre_computed_denom_left_mag: np.ndarray  # (k+1)^(1j*Im(-s))
pre_computed_denom_right_phase = []  # (k+1)^(Re(-(1-s)))
pre_computed_denom_left_phase = []  # (k+1)^(1j*Im(-(1-s)))
pre_computed_func1_mag: np.ndarray  # (2*pi)^Re(s) / pi for s < 0
pre_computed_func1_phase: np.ndarray  # phase for above, for s < 0
pre_computed_func2_mag: np.ndarray  # n/a
pre_computed_func2_phase: np.ndarray  # n/a


# Precomputes denominator coefficients, resulting in a roughly 4x speedup
# Assume mesh is a 2d ndarray, i.e. a list of lists
def precompute_denom(mesh):
    global array_zero_split, \
        pre_computed_func1_mag, pre_computed_func1_phase, \
        pre_computed_func2_mag, pre_computed_func2_phase, \
        pre_computed_denom_left_mag, pre_computed_denom_right_mag

    # Computes (k + 1) ^ -s for all s in grid.
    # Works by decomposing s = a + bi, then precomputing the following:
    #   (k + 1) ^ a   = varies by x position but not y
    #   (k + 1) ^ bi  = varies by y position but not x

    # Look at bottom row
    row1 = mesh[0]
    # Find where Re(s) switches to positive
    array_zero_split = bisect.bisect_left(row1, 0)
    # Put all negative values (< 0) in one array
    row1_neg_vals = row1[0:array_zero_split]
    # Put all non-negative (>= 0) values in a second array
    row1_pos_vals = row1[array_zero_split:len(row1)]

    # Delete any previous coefficients
    pre_computed_denom_left_phase.clear()
    pre_computed_denom_right_phase.clear()
    pre_computed_denom_left_mag = np.empty((RIEMANN_ITER_LIMIT, len(row1_neg_vals)), dtype=np.complex)
    pre_computed_denom_right_mag = np.empty((RIEMANN_ITER_LIMIT, len(row1_pos_vals)), dtype=np.complex)

    for k in range(0, RIEMANN_ITER_LIMIT):
        # Denominator magnitudes as a function of x. This array is sorted by k, then column
        pre_computed_denom_right_mag[k] = np.power(k + 1, -np.real(row1_pos_vals))
        pre_computed_denom_left_mag[k] = np.power(k + 1, -np.real(1 - row1_neg_vals))

    col1 = np.empty(len(mesh), dtype=np.complex)
    for row in range(0, len(mesh)):
        # Denominator phase is a function of y (row). This array is sorted by row, then k
        col1[row] = mesh[row][0]
        s = col1[row]
        k_list = np.linspace(1, RIEMANN_ITER_LIMIT, RIEMANN_ITER_LIMIT)
        pre_computed_denom_right_phase.append(np.power(k_list, -1j * np.imag(s)))
        pre_computed_denom_left_phase.append(np.power(k_list, -1j * np.imag(1 - s)))

    # Left half plane [(2*pi) ^ s] / pi
    pre_computed_func1_mag = np.power(2 * np.pi, np.real(row1_neg_vals)) / np.pi  # Magnitude is function of x
    pre_computed_func1_phase = np.power(2 * np.pi, 1j * np.imag(col1))  # Phase is function of y


#
# Calculate Riemann zeta function for a 1D ndarray (i.e. vector)
#
# Assumes elements are all non-negative, or all negative (in the latter case, will use functional equation
# to assure convergence) and are in the same order as in the pre_computed denominator arrays, i.e. that they
# are consecutive points in a single row of a mesh_grid
def riemann_array(s,
                  row_num,  # Which row of mesh grid? 0 is bottom
                  do_eta=False,  # Do Dirichlet eta instead of Riemann zeta
                  left_half_plane=False):  # True if called RECURSIVELY for Re(s) < 0. False otherwise, even if Re(s)<0

    global quit_computation_flag

    if len(s) == 0:
        return []

    if np.real(s[0]) < 0:
        # Use functional equation if Re[s] is negative, or else sum won't converge
        # Only checks the first element - we assume all are negative, or all are non-negative
        if do_eta:
            # Must set left_half_plane=True because we have replaced s with 1 - s, so arg values are now
            # all positive
            return -2 * s * riemann_array(1 - s, row_num, do_eta=do_eta, left_half_plane=True) \
                   * gamma(- s) * np.sin(-s * np.pi / 2) * np.power(np.pi, s - 1) \
                   * (1 - np.power(2, s - 1)) / (1 - np.power(2, s))
        else:
            if USE_CACHED_FUNC:
                # Combine and cached exponential function 2^s * pi ^ (s-1) to speed up.
                # Gives another 30-40% improvement when added on top of all earlier optimizations!
                # Now have ~2M samples/second on home office computer
                # Theoretically we could also cache sin since sin(a+bi) = sinacoshb + icosasinhb, but
                # that requires generating 4 more pre-calculated arrays, which feels like a pain
                return riemann_array(1 - s, row_num, do_eta=do_eta, left_half_plane=True) \
                       * gamma(1 - s) * np.sin(s * np.pi / 2) \
                       * pre_computed_func1_mag * pre_computed_func1_phase[row_num]
            else:
                # Old method is foolproof, and only slightly slower
                return riemann_array(1 - s, row_num, do_eta=do_eta, left_half_plane=True) \
                       * gamma(1 - s) * np.sin(s * np.pi / 2) \
                       * np.power(np.pi, s - 1) \
                       * np.power(2, s)  # This can be combined with above exponential. Keep this way for readability

    # We are in zone where sum will converge.

    if USE_CACHED_FUNC:
        if left_half_plane:
            try:
                cum_sum = np.dot(np.multiply(NK2_array, pre_computed_denom_left_phase[row_num]),
                                 pre_computed_denom_left_mag)  # Dot product is much faster than loop
            except:
                print("Error during riemann_arrays. Do you have two concurrent operations underway?")
                quit_computation_flag = True
                return []

            #        phase = pre_computed_denom_left_phase[row_num]
            #        for k in range(0, RIEMANN_ITER_LIMIT):
            #            cum_sum += NK2_array[k] * pre_computed_denom_left_mag[k] * phase[k]
            if do_eta:
                return cum_sum
            else:
                return cum_sum / (1 - 2 * pre_computed_denom_left_mag[1] * pre_computed_denom_left_phase[row_num][1])
        else:
            cum_sum = np.dot(np.multiply(NK2_array, pre_computed_denom_right_phase[row_num]),
                             pre_computed_denom_right_mag)  # Dot product is much faster than loop
            #             phase = pre_computed_denom_right_phase[row_num]
            #             for k in range(0, RIEMANN_ITER_LIMIT):
            #                 cum_sum += NK2_array[k] * pre_computed_denom_right_mag[k] * phase[k]
            if do_eta:
                return cum_sum
            else:
                return cum_sum / (1 - 2 * pre_computed_denom_right_mag[1] * pre_computed_denom_right_phase[row_num][1])

    else:
        # Old version: about 10x slower since we recompute exponential every time
        # We still benefit from vectorized operations, which gave 100x improvement.
        cum_sum = [0 + 0j] * len(s)
        for k in range(0, RIEMANN_ITER_LIMIT):
            cum_sum += NK2_array[k] * np.power(k + 1, - s)  # Calculate the denominator the traditional (4x slower) way.

        if do_eta:
            return cum_sum
        else:
            # Scale factor converts Dirichlet eta function to Riemann zeta function
            return cum_sum / (1 - np.power(2, 1 - s))


#
# When multiplied by gamma(s/2) * pi ^(-s/2), the result has 180-deg rotational symmetry around s = 0.5 + 0j
#
def riemann_symmetric(s):
    r = riemann(s)
    if quit_computation_flag:
        # If we were interrupted, we will have a "ragged" list, in which some elements are a single "0",
        # while others are a full length. Don't try
        # to multiple this by gamma stuff, or else will get an error that we can't multiply a sequence
        # by a single number
        return r

    return r * gamma(s / 2) * (np.pi ** (-s / 2))


# Not currently used
def riemann_gamma(s):
    return riemann(s) * gamma(s / 2)


# Calculate gamma(s) while also updating progress bar
def gamma_with_progress(s):
    global row_count

    if s.ndim == 2:
        # This will convert 2d ndarray into list of 1d ndarrays
        out = [0.0 if quit_computation_flag else gamma_with_progress(x) for x in s]
        # Convert back to single 2d ndarray
        return np.stack(out)

    # Now we have 1D vector
    row_count += 1

    if row_count == 50:
        # Update progress bar every 50 rows
        row_count = 0
        if computation_progress_callback is not None:
            computation_progress_callback(np.size(s) * 50)

    # Call vectorized version
    return gamma(s)
