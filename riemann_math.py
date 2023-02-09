import typing
import bisect
import numpy as np
from scipy.special import gamma, loggamma
import time

outArray = []

# This should be set by caller, to handle computational progress
computation_progress_callback: typing.Optional[typing.Callable] = None

# This can be set by caller, to halt ongoing computation (only works with Riemann)
quit_computation_flag = False

# True if we compute Riemann zeta using cached exponentials (gives ~10x speedup)
USE_CACHED_FUNC = True

RIEMANN_ITER_LIMIT = 40  # Default. Can be overridden
NK2_array = []  # 1d array of partial sums of NK1_array
D_array = []  # 1d array of coefficients for Borwein method
B_array = []  # 1d array of coefficients for Borwein method

# Counts how many rows of riemann() have been computed so far.
# Used to determine when to update progress bar, and also tracks which precomputed denominator coefficients to use.
row_count = 0

# Used by riemann() to track where Re(s) switches from negative to positive. We store this when first row is calculated,
# so that we can ensure split occurs at exactly the same place for every row. Otherwise rounding errors might cause
# slight differences between rows
array_zero_split = 0
elapsed_time = 0
entry_count = 0
use_complex256 = False  # This supports higher precision calculations, but only seems to work on MacOS or Linux


# Precompute table of coefficients for lookup using Euler's transformation.
# This reduces Riemann computation time from O(n^2) to O(n) where n is number of terms
def precompute_coeffs():
    global NK2_array

    if len(NK2_array) == RIEMANN_ITER_LIMIT:
        # No need to recalculate
        return

    t1 = time.time()

    if use_complex256:
        val_type = np.complex256
    else:
        val_type = complex

    NK2_array = np.zeros(RIEMANN_ITER_LIMIT, dtype=val_type)
    #    print("Precomputing " + str(RIEMANN_ITER_LIMIT) + " coefficients for Riemann/Dirichlet sum", end="")

    # Precompute a "cumulative", "normalized" Pascal's triangle,
    # in which each row is half the sum of the two values above it,
    # and the first item in each row is 1.
    #
    # The values in each row are the cumulative sum of values in the traditional Pascal's triangle,
    # normalized to 1/2^n.
    #
    # The Euler coefficients will be the values in the last row, excluding the first item
    # 1
    # 1   1/2
    # 1   3/4   1/4
    # 1   7/8   4/8   1/8
    # 1   15/16 11/16 5/16   1/16
    # 1   31/32 26/32 16/32  6/32  1/32
    # 1   63/64 57/64 42/64 22/64  7/64  1/64
    # ...
    NK1_array = np.zeros(RIEMANN_ITER_LIMIT + 1, dtype=val_type)
    NK1_array_tmp = np.zeros(RIEMANN_ITER_LIMIT + 1, dtype=val_type)
    NK1_array[0] = 1  #

    for n in range(1, RIEMANN_ITER_LIMIT + 1):
        NK1_array_tmp[0] = 1
        NK1_array_tmp[1:(n + 1)] = (NK1_array[0:n] + NK1_array[1:(n + 1)]) / 2
        NK1_array = NK1_array_tmp

    # Exclude first value in the above cumulative sum, and alternate sign
    for k in range(0, RIEMANN_ITER_LIMIT):
        NK2_array[k] = ((-1) ** k) * NK1_array[k + 1]

    delay = time.time() - t1
    print("Precomputed eta coefficients for %d iterations in %1.4f seconds " % (RIEMANN_ITER_LIMIT, delay))


# Precompute table of coefficients for lookup using Euler's transformation.
# This reduces Riemann computation time from O(n^2) to O(n) where n is number of terms
def precompute_borwein():
    global B_array, D_array

    if len(B_array) == RIEMANN_ITER_LIMIT:
        # No need to recalculate
        return

    t1 = time.time()

    B_array = np.zeros(RIEMANN_ITER_LIMIT, dtype=np.clongdouble)
    D_array = np.zeros(RIEMANN_ITER_LIMIT + 1, dtype=np.clongdouble)
    #    print("Precomputing " + str(RIEMANN_ITER_LIMIT) + " coefficients for Riemann/Dirichlet sum", end="")

    # Precompute N_choose_k / 2^(N+1) coefficients.
    NK_array = np.zeros(shape=(RIEMANN_ITER_LIMIT * 2 + 1, RIEMANN_ITER_LIMIT * 2 + 1), dtype=np.clongdouble)
    NK_array[0, 0] = 1  # This will be (n_choose_k)

    # Compute n choose k. This is faster because it uses addition instead of multiplication/factorials
    for n in range(1, RIEMANN_ITER_LIMIT * 2 + 1):
        NK_array[n, 0] = NK_array[n - 1, 0]
        NK_array[n, 1:(n + 1)] = NK_array[n - 1, 0:n] + NK_array[n - 1, 1:(n + 1)]

    # Compute borwein coefficients
    for k in range(0, RIEMANN_ITER_LIMIT + 1):
        D_array[k] = 0
        for l1 in range(0, k + 1):
            D_array[k] += NK_array[RIEMANN_ITER_LIMIT + l1,
                                   RIEMANN_ITER_LIMIT - l1] / (RIEMANN_ITER_LIMIT + l1) * (4.0 ** l1)

        D_array[k] *= RIEMANN_ITER_LIMIT

    for k in range(0, RIEMANN_ITER_LIMIT):
        B_array[k] = (-1 ** k) * (D_array[k] / D_array[RIEMANN_ITER_LIMIT] - 1)

    delay = time.time() - t1
    print("Precomputed Borwein coefficients for %d iterations in %1.4f seconds " % (RIEMANN_ITER_LIMIT, delay))


def eta_zeta_scale(v):
    # Scale factor converts Dirichlet eta function to Riemann zeta function
    return 1 / (1 - 2 ** (1 - v))


def riemann_borwein(s):
    if len(B_array) != RIEMANN_ITER_LIMIT:
        precompute_borwein()

    cum_sum = np.clongdouble(0)
    s = complex(s)  # Because np.power will reject negative integer power (but float or complex is OK)
    for k in range(0, RIEMANN_ITER_LIMIT):
        cum_sum += B_array[k] * np.power(k + 1, -s)

    scale1 = 1 / (1 - 2 ** (1 - s))
    return cum_sum * scale1


#
# Calculate Riemann zeta function for complex input s.
#
# If get_array_size is True, then also return array size of intermediate sum calculation (used to draw arrows)
# This option requires Re[s] >= 0.
#
# If s is a 2D array, values must have mesh arrangement, in which Re[s] is constant for all values in any column,
# and Imag[s] is constant for all values in any row. Also, all values having Re[s] < 0 must be in columns to left
# values having Re[s] >= 0
#
# If s is 1D array, then is_vertical option can be set True, to specify that values all have same Re[s], and evenly
# spaced Imag[s]. This allows faster computation of k-powers.
#
# If return_log is true, then return log[Riemann[s]], instead of Riemann[s]. This avoids problem where Riemann[s]
# rounds to zero for large abs[s]
#
# Approximate accuracy for real s is (very) roughly proportional to magnitude of final summation
# term, which is about 1/(2^ITER * ITER^s). Hence, number of digits of precision is roughly
# log_10(2^ITER * ITER^s) = .301 * ITER + s * log_10(ITER), e.g. 30 + 2*Re[s] digits for ITER = 100,
# and roughly 300 + 3*Re(s) digits for ITER = 1000.
#
# Accuracy goes down when Imag[s] is large, so we need ITER roughly equal to Imag[s]
#
# Note that convergence is very poor for Re(s) large and negative. So we use the functional equation
# for Re(s) < 0.
#
def riemann(s, get_array_size=False,
            do_eta=False,  # Calculate dirichlet eta function instead of Riemann zeta
            use_zero_for_nan=True,  # This is now ignored
            is_vertical=False,  # Applies to 1D array only. Specifies values with fixed Re[s] and evenly spaced imag[s]
            return_log=False):  # If true, return log(Riemann(s)) instead of Riemann(s). Useful for large Abs[s]
    global row_count, array_zero_split, elapsed_time, entry_count

    if np.size(s) > 1:
        if s.ndim == 2:  # Somehow, I get error if I merge the comparisons with "and" and if s is single value or list

            # Handle 2D array. Must have mesh structure: imag part depends on row, real part depends on column
            if np.max(np.abs(s)) >= 400 and not return_log:
                # Switch to log algorithm
                r = riemann(s, get_array_size=get_array_size, do_eta=do_eta,
                            use_zero_for_nan=use_zero_for_nan,
                            return_log=True)

                # For any complex c, we have: log(abs(c)) = Re(log(c))
                # Hence, we can reduce magnitude of s by reducing the real portion of log(s)
                # Compress magnitude range by 1000-fold, or else will round to zero for high values
                r = r - np.real(r) * 0.999
                # Convert log back to "normal"
                return np.exp(r)

            if entry_count >= 1:
                # Prevent this from running more than once concurrently
                print("Warning: can't have two concurrent riemann calculations at the same time.")
                return 0

            entry_count += 1

            t1 = time.time()
            precompute_denom(s)
            delay = time.time() - t1
            print("Precomputed powers in " + "{:1.4f}".format(delay) + " seconds")
            row_count = 0
            t1 = time.time()
            # When doing heatmaps, s is initially a 2d ndarray, which the following converts to list of 1d ndarrays.
            out = [0.0 if quit_computation_flag else riemann_row(x, do_eta=do_eta, return_log=return_log) for x in s]

            entry_count -= 1
            if not quit_computation_flag:
                elapsed_time = time.time() - t1
                # Convert list of 1d arrays back to single 2d ndarray
                return np.stack(out)
            else:
                # User has quit, just return
                elapsed_time = 0
                return None

    #
    # Either have single value, or 1d array.

    # Heatmap calculations never come here, since their input is 2D array, which is handled above.
    if np.size(s) > 1.0:

        # Make sure we don't retain any negative real components
        s[np.real(s) < 0] = 0

        # Handle 1d array using vectorized operations (but without cached denominators)
        out_vals = riemann_row_non_negative(s, row_num=-1,
                                            do_eta=do_eta,
                                            is_vertical=is_vertical,
                                            USE_CACHED=False)

        if return_log:
            return np.log(out_vals)
        else:
            return out_vals

    #        return [riemann(x, do_eta=do_eta, use_zero_for_nan=use_zero_for_nan) for x in s]

    # If we are here, we are computing for a single value. Use slow method, without vectorization
    # or cached k powers
    if np.size(s) == 1 and s == 1.0:
        # Calculation blows up at 1.0, so don't bother to calculate
        if use_zero_for_nan:
            # Do this to avoid warnings.
            return 0
        else:
            return np.nan

    if np.real(s) < 0:
        # Use functional equation
        if do_eta:
            if return_log:
                return np.log(-2 * s) + riemann(1 - s, do_eta=True, return_log=True) \
                    + loggamma(- s) + logsin(-s * np.pi / 2) + (s - 1) * np.log(np.pi) \
                    + np.log((1 - 2 ** (s - 1)) / (1 - 2 ** s))
            else:
                return -2 * s * riemann(1 - s, do_eta=True, return_log=False) \
                    * gamma(- s) * np.sin(-s * np.pi / 2) * (np.pi ** (s - 1)) \
                    * (1 - 2 ** (s - 1)) / (1 - 2 ** s)
        else:
            if return_log:
                return riemann(1 - s, do_eta=False, return_log=True) \
                    + loggamma(1 - s) + logsin(s * np.pi / 2) + (s - 1) * np.log(np.pi) + s * np.log(2)
            else:
                return riemann(1 - s, do_eta=False, return_log=False) \
                    * gamma(1 - s) * np.sin(s * np.pi / 2) * (np.pi ** (s - 1)) * (2 ** s)

    cum_sum = 0 + 0j

    if return_log:
        # Have single value
        s = complex(s)  # np.power can't handle negative integer exponents, so make sure it is complex and float
        for k in range(0, RIEMANN_ITER_LIMIT):
            cum_sum += NK2_array[k] * np.power(k + 1, -s)

        if do_eta:
            return np.log(cum_sum)
        else:
            # Scale factor converts Dirichlet eta function to Riemann zeta function
            return np.log(cum_sum) - np.log(1 - np.power(2, 1 - s))

    # Need first element zero so line segment will draw correctly
    store_intermediates = len(outArray) > 0
    if store_intermediates:
        outArray[0] = cum_sum

    if do_eta:
        scale1 = 1
    else:
        # Scale factor converts Dirichlet eta function to Riemann zeta function
        scale1 = 1 / (1 - 2 ** (1 - s))

    if len(NK2_array) != RIEMANN_ITER_LIMIT:
        precompute_coeffs()

    if store_intermediates or get_array_size:

        # These two options always calculates riemann, and are incompatible with do_eta option
        # get_array_size is incompatible with negative real s, since it requires two return values whereas
        # functional equation won't return the second return value
        partial_sum_index = 1
        have_final_point = False

        # Calculate terms of Dirichlet eta function, while storing intermediate values.
        # Apply scale1 factor to get Riemann zeta
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

        # Make sure final point is included
        if not have_final_point:
            partial_sum_index += 1
            outArray[partial_sum_index] = cum_sum

        if get_array_size:
            return cum_sum, partial_sum_index
        else:
            return cum_sum

    else:
        # If not saving intermediate values, use faster simple loop
        for k in range(0, RIEMANN_ITER_LIMIT):
            cum_sum = cum_sum + NK2_array[k] / ((k + 1) ** s)

        if return_log:
            # Convert to log, if requested
            cum_sum = np.log(cum_sum)

        if not do_eta:
            # Scale factor converts Dirichlet eta function to Riemann zeta function
            if return_log:
                cum_sum = cum_sum + np.log(scale1)
            else:
                cum_sum = cum_sum * scale1

        return cum_sum


#
# For c = a + bi, calculates sin(c)/exp(abs(b))
#
# This generally has magnitude ~1.0, avoiding the problem where sin(c) blows up for large b.
#
# This is based on the identity:
#
# sin(a+bi) = sin(a)cosh(y) + icos(x)sinh(y)
#
# which expands to:
#
# sin(a+bi) = e^abs(b)[sin(a)[1+e^(-2abs(b))]/2 + icos(a)sign(b)[1-e^(-2abs(b)]/2
#
def sin2(c):
    a = np.real(c)
    b = np.imag(c)
    return np.sin(a) * (1 + np.exp(-2 * np.abs(b))) / 2 + 1j * np.cos(a) * np.sign(b) * (1 - np.exp(-2 * np.abs(b))) / 2


#
# For c = a + bi, calculate log(sin(c))
#
# Because sin(c) blows up for large b, we use the sin2 function above, which pulls out the large e^abs(b)
# term. We can then calculate log(sin(c)) using the identity:
#
# log(sin(c)) = log[e^abs(b) * sin2(c)] = abs(b) + log(sin2(c))
def logsin(c):
    return np.abs(np.imag(c)) + np.log(sin2(c))


# Handles one row of 2d matrix using optional cached powers.
# If inputs have both positive and negative real portions, will call riemann_row_non_negative with
# negative s replaced by 1-s, and then apply functional equation to correct those values
def riemann_row(s, do_eta=False, USE_CACHED_DENOM=True, return_log=False):
    global row_count, array_zero_split

    # We have 1D vector. Update progress bar, then call Riemann for individual values
    row_count += 1

    if np.mod(row_count, 50) == 0:
        # Update progress bar every 50 rows
        if computation_progress_callback is not None:
            computation_progress_callback(np.size(s) * 50)

    # RiemannArray requires elements to be non-negative. So we split s into two arrays.
    s0 = s[0:array_zero_split]
    s1 = s[array_zero_split:len(s)]

    # Transform negative values into 1-s, and make new array
    in_vals = np.concatenate(((1 - s0), s1))

    out_vals = riemann_row_non_negative(in_vals, row_count - 1, do_eta=do_eta, USE_CACHED=USE_CACHED_DENOM)

    if array_zero_split > 0:
        out_vals[0:array_zero_split] = np.conj(out_vals[0:array_zero_split])

    if return_log:
        # If reporting log(Riemann) then also use loggamma and logsin
        out_vals = np.log(out_vals)

        if array_zero_split > 0:
            # Use functional equation to correct results for Re[s] < 0
            if do_eta:
                out_vals[0:array_zero_split] += np.log(-2 * s0) \
                                                + loggamma(- s0) + logsin(-s0 * np.pi / 2) \
                                                + np.log(np.pi) * (s0 - 1) \
                                                + np.log((1 - np.power(2, s0 - 1)) / (1 - np.power(2, s0)))
            else:
                out_vals[0:array_zero_split] += loggamma(1 - s0) + logsin(s0 * np.pi / 2) \
                                                + np.log(np.pi) * (s0 - 1) + np.log(2) * s0

        out_vals = np.real(out_vals) + 1j * (np.remainder(np.imag(out_vals) + np.pi, 2 * np.pi) - np.pi)
        return out_vals

    if array_zero_split > 0:
        # Use functional equation to correct results for Re[s] < 0
        if do_eta:
            out_vals[0:array_zero_split] *= -2 * s0 \
                                            * gamma(- s0) * np.sin(-s0 * np.pi / 2) \
                                            * np.power(np.pi, s0 - 1) \
                                            * (1 - np.power(2, s0 - 1)) / (1 - np.power(2, s0))
        else:
            if USE_CACHED_FUNC:
                try:
                    # Use combined+cached exponential function 2^s * pi ^ (s-1) to speed up.
                    # Gives another 30-40% improvement when added on top of all earlier optimizations!
                    # Now have ~2M samples/second on home office computer
                    # Theoretically we could also cache sin since sin(a+bi) = sin(a)cosh(b) + i*cos(a)sinh(b), but
                    # that requires generating 4 more pre-calculated arrays, which feels like a pain
                    out_vals[0:array_zero_split] *= gamma(1 - s0) * np.sin(s0 * np.pi / 2) \
                                                    * cached_pi_power_mag * cached_pi_power_phase[row_count - 1]
                except:
                    print("error")
            else:
                # Old method is foolproof, and only slightly slower
                out_vals[0:array_zero_split] *= gamma(1 - s0) * np.sin(s0 * np.pi / 2) \
                                                * np.power(np.pi, s0 - 1) \
                                                * np.power(2, s0)  # Can be combined with above exponential.

    if return_log:
        return np.log(out_vals)

    #            a1 = riemann_array(neg_array, row_count - 1, do_eta=do_eta)
    #            a2 = riemann_array(pos_array, row_count - 1, do_eta=do_eta)
    return out_vals


cached_powers_mag: np.ndarray  # 2d array. First arg is k, second is mesh col #. (k+1)^(Re(-s))
cached_powers_phase = []  # 2d array. First arg is mesh row #, second is k. (k+1)^(Re(-(1-s)))
cached_pi_power_mag: np.ndarray  # (2*pi)^Re(s) / pi for s < 0
cached_pi_power_phase: np.ndarray  # phase for above, for s < 0


# Precomputes denominator coefficients, resulting in a roughly 4x speedup
# Assume mesh is a 2d ndarray, i.e. a list of lists
def precompute_denom(mesh):
    global array_zero_split, \
        cached_pi_power_mag, cached_pi_power_phase, \
        cached_powers_mag, cached_powers_phase

    # Computes (k + 1) ^ -s for all s in grid.
    # Works by decomposing s = a + bi, then precomputing the following:
    #   (k + 1) ^ a   = varies by x position but not y
    #   (k + 1) ^ bi  = varies by y position but not x

    # Look at bottom row
    row1 = mesh[0]
    # Find where Re(s) switches to positive
    array_zero_split = bisect.bisect_left(np.real(row1), 0)
    # Put all negative values (< 0) in one array
    s0 = row1[0:array_zero_split]
    # Put all non-negative (>= 0) values in a second array
    s1 = row1[array_zero_split:len(row1)]

    # Delete any previous coefficients
    cached_powers_phase.clear()
    cached_powers_mag = np.empty((RIEMANN_ITER_LIMIT, len(row1)), dtype=complex)

    for k in range(0, RIEMANN_ITER_LIMIT):
        # Denominator magnitudes as a function of x. This array is sorted by k, then column
        cached_powers_mag[k] = np.concatenate((np.power(k + 1, -np.real(1 - s0)),
                                               np.power(k + 1, -np.real(s1))))

    col1 = np.empty(len(mesh), dtype=complex)
    for row in range(0, len(mesh)):
        # Denominator phase is a function of y (row). This array is sorted by row, then k
        col1[row] = mesh[row][0]
        s = col1[row]
        k_list = np.linspace(1, RIEMANN_ITER_LIMIT, RIEMANN_ITER_LIMIT)
        cached_powers_phase.append(np.power(k_list, -1j * np.imag(s)))

    # Left half plane [(2*pi) ^ s] / pi
    cached_pi_power_mag = np.power(2 * np.pi, np.real(s0)) / np.pi  # Magnitude is function of x
    cached_pi_power_phase = np.power(2 * np.pi, 1j * np.imag(col1))  # Phase is function of y


USE_MATRIX_MULT = True
bases = []  # List of numbers k=1, 2, 3, ...
powers = []  # List of values k^s_imag, where s_imag = Imag(s)
imag_part = []
MAX_MEMORY = 200000000  # max # of complex elements in array


def make_powers(s):
    global powers, bases, imag_part

    len_s = len(s)
    imag_part = np.imag(s)

    if len(bases) != RIEMANN_ITER_LIMIT:
        bases = np.array([np.arange(1, RIEMANN_ITER_LIMIT + 1, 1, dtype=complex)]).T

    if len_s * RIEMANN_ITER_LIMIT > MAX_MEMORY:
        len_s = int(MAX_MEMORY / RIEMANN_ITER_LIMIT)

    if len(powers) == 0:
        powers = np.empty((RIEMANN_ITER_LIMIT, len_s), dtype=complex)
    elif len_s > powers.shape[1]:
        powers = np.empty((RIEMANN_ITER_LIMIT, len_s), dtype=complex)

    make_powers2(len_s)


#
# Recursively pre-computes k^imag(s) for k = 1, 2, 3, ... ITER_LIMIT
#
def make_powers2(limit):
    global bases, powers

    if limit == 1:
        powers[:, 0] = np.power(bases, -1j * imag_part[0]).flatten()
        return

    half_s = int(np.ceil(limit / 2))
    is_odd = np.mod(limit, 2) == 1

    make_powers2(half_s)  # np.power(bases2, -im_part * 1j)

    # Now calculate second half
    imag_diff = imag_part[half_s] - imag_part[0]  # Imaginary diff
    ratio = np.power(bases, -1j * imag_diff)

    if is_odd:
        limit1 = half_s - 1
    else:
        limit1 = half_s

    # Second half of range is obtained by scalar multiplication from the first half
    powers[:, half_s:limit] = np.multiply(ratio, powers[:, 0:limit1])


#
# Calculate Riemann zeta function for a 1D ndarray (i.e. vector) of non-negative values.
#
# Uses vectorized operations and (optional) pre-computed (cached) power array.
# Note that this does not call gamma(), hence will work for very large input values
#
def riemann_row_non_negative(s,
                             row_num,  # Which row of mesh grid? 0 is bottom
                             do_eta=False,  # Do Dirichlet eta instead of Riemann zeta
                             USE_CACHED=True,
                             is_vertical=False):  # True if evenly spaced along vertical line
    global quit_computation_flag
    global bases, powers

    if len(s) == 0:
        return []

    if USE_CACHED:

        if USE_MATRIX_MULT:
            # Calculate k ^ (-s) for k = 1, 2, 3, ... ITER
            #
            # We can pre-compute the above for the real and imaginary parts of s, then multiply them together.
            #
            # NK2_array has shape (ITER, ), and holds the Euler weights
            #
            # cached_powers_phase is N x ITER matrix, where N is # of plot mesh rows
            # cached_powers_phase[row_num] has shape (ITER, ), with elements representing k (= 1, 2, ... ITER)
            # raised to power of the imaginary component of current row
            #
            # np.multiply() performs element by element multiplication with phase. Output is then multiplied
            # (again element by element) with magnitude, then summed (hence we use dot product instead of np.multiply)
            # to give final output
            #
            cum_sum = np.dot(np.multiply(NK2_array, cached_powers_phase[row_num]),
                             cached_powers_mag)
        else:
            phase = cached_powers_phase[row_num]
            cum_sum = 0
            for k in range(0, RIEMANN_ITER_LIMIT):
                cum_sum += NK2_array[k] * cached_powers_mag[k] * phase[k]
        if do_eta:
            return cum_sum
        else:
            return cum_sum / (1 - 2 * cached_powers_mag[1] * cached_powers_phase[row_num][1])

    else:

        if is_vertical:  #
            # Use pre-computed k-power table, allowing 2D vectorization, and further speedup

            # In vertical mode, all values of s have the same real part. So calculate for the first
            # one and assume it is same for all other points
            real_part = np.real(s[0])
            bases_real = np.power(bases, -real_part)
            NK2_real_power = np.multiply(NK2_array, bases_real.flatten())

            if USE_MATRIX_MULT:
                len_s = len(s)
                if powers.shape[1] == len_s:
                    cum_sum = np.dot(np.multiply(NK2_array, bases_real.flatten()), powers)
                else:
                    startPoint = 0
                    width = powers.shape[1]
                    cum_sum = np.empty(len_s, dtype=complex)
                    while startPoint <= len_s:
                        ratios = np.power(bases, -(imag_part[startPoint] - imag_part[0]) * 1j)
                        NK2_power = np.multiply(NK2_real_power, ratios.flatten())
                        endPoint = startPoint + width
                        if endPoint <= len_s:
                            cum_sum[startPoint:endPoint] = np.dot(NK2_power, powers)
                        else:
                            cum_sum[startPoint:len_s] = np.dot(NK2_power, powers[:, 0:(len_s - startPoint)])
                        startPoint += width
            else:
                # Use slower method that requires lots of complex exponential calculations
                # Replicate this column for each input value, creating a 2D array of size RIEMANN_ITER_LIMT x len(s)
                bases2 = np.repeat(bases, len(s), axis=1)
                # Raise all numbers in 2D array to each value in s. In matlab, we would do bsxfun(@power, bases, s)
                # on the 1D arrays, and would not need extra step of creating 2D array.
                powers = np.power(bases2, -np.imag(s) * 1j)
                # Now do matrix product
                cum_sum = np.dot(np.multiply(NK2_array, bases_real), powers)

            # Strangely, although np.dot() is much faster than loop, the 2D vectorized np.power()
            # is actually very slow, and takes LONGER than the loop below. To really speed things up, need
            # to pre-compute k-powers.
        else:
            # Use 1D vectorization, which is MUCH faster than nothing at all.
            # There is still a loop across k's. I tried vectorizing that dimension also, but it actually made things
            # 10-50% slower. By comparison matlab did the 2D vectorized operation in about 5 seconds on home office
            # computer, about 5x faster.
            cum_sum = [0 + 0j] * len(s)

            for k in range(0, RIEMANN_ITER_LIMIT):
                cum_sum += NK2_array[k] * np.power(k + 1, -s)  # Calculate power the traditional (4x slower) way.

        if do_eta:
            return cum_sum
        else:
            # Scale factor converts Dirichlet eta function to Riemann zeta function
            return cum_sum / (1 - np.power(2, 1 - s))


#
# Calculate Riemann Xi function, which is Riemann(s) * gamma(s/2) * pi ^(-s/2) * (1/2)s(s-1)
#
def riemann_symmetric(s, return_log=False, is_vertical=False):
    if not return_log:

        if np.max(abs(s)) > 400:
            # Switch to log algorithm
            r = riemann_symmetric(s, return_log=True, is_vertical=is_vertical)

            # Compress magnitude range by 100-fold
            r = r - np.real(r) * 0.99
            return np.exp(r)

        r = riemann(s, is_vertical=is_vertical)
        if quit_computation_flag:
            # If we were interrupted, we will have a "ragged" list, in which some elements are a single "0",
            # while others are a full length vector. Don't try to multiply by anything, as array sizes won't match.
            return r
        return 0.5 * s * (s - 1) * r * gamma(s / 2) * (np.pi ** (-s / 2))

    # For very large abs(s), we want to calculate log

    # For any complex c = a+bi
    #     log(c) = log(abs(c)) + i*Arg(c)
    # For any complex c = r * e^(it)
    #     log(c) = log(r) + it
    return np.log(0.5) + np.log(s) + np.log(s - 1) \
        + riemann(s, return_log=True, is_vertical=is_vertical) \
        + loggamma(s / 2) - s * np.log(np.pi) / 2


#
# Calculate a modified Riemann function that is always real on the critical line
#
def riemann_real(s, is_vertical=False):
    r = riemann(s, is_vertical=is_vertical)
    if quit_computation_flag:
        # If we were interrupted, we will have a "ragged" list, in which some elements are a single "0",
        # while others are a full length vector. Don't try to multiply by anything, as array sizes won't match.
        return r

    # Use gamma function to "unwind" Riemann zeta phase
    gp = loggamma(s / 2)
    # Force gamma magnitude to be one, since it becomes extremely small for large Imag(s)
    gp = gp - np.real(gp)

    # Unwind phase using gamma phase and pi ^ (-s/2)
    return r * (np.pi ** (-s / 2)) * np.exp(gp)


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

    if np.max(np.abs(s)) >= 400:
        r = loggamma(s)
        r = r - np.real(r) * 0.9
        return np.exp(r)

    # Call vectorized version
    return gamma(s)
