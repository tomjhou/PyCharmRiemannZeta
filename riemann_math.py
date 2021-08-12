import typing
import bisect
import numpy as np
from scipy.special import gamma, loggamma
import time

outArray = []

# This should be set by caller, to handle computational progress
computation_progress_callback: typing.Callable = None

# This can be set by caller, to halt ongoing computation (only works with Riemann)
quit_computation_flag = False

# True if we compute Riemann zeta using cached exponentials (gives ~10x speedup)
USE_CACHED_FUNC = True

RIEMANN_ITER_LIMIT = 40  # Default. Can be overridden
NK1_array = []  # 2d array of values (n_choose_k) / 2^(n+1)
NK2_array = []  # 1d array of partial sums of NK1_array

# Counts how many rows of riemann() have been computed so far.
# Used to determine when to update progress bar, and also tracks which precomputed denominator coefficients to use.
row_count = 0

# Used by riemann() to track where Re(s) switches from negative to positive. We store this when first row is calculated,
# so that we can ensure split occurs at exactly the same place for every row. Otherwise rounding errors might cause
# slight differences between rows
array_zero_split = 0
elapsed_time = 0
entry_count = 0


# Precompute table of coefficients for lookup using Euler's transformation.
# This reduces Riemann computation time from O(n^2) to O(n) where n is number of terms
def precompute_coeffs():
    global NK2_array, NK1_array

    if len(NK2_array) == RIEMANN_ITER_LIMIT:
        # No need to recalculate
        return

    VECTORIZED = True

    t1 = time.time()

    NK2_array = np.zeros(RIEMANN_ITER_LIMIT)
    #    print("Precomputing " + str(RIEMANN_ITER_LIMIT) + " coefficients for Riemann/Dirichlet sum", end="")

    # Precompute N_choose_k / 2^(N+1) coefficients.
    NK1_array = np.zeros(shape=(RIEMANN_ITER_LIMIT, RIEMANN_ITER_LIMIT))
    NK1_array[0, 0] = 0.5  # This will be (n_choose_k) / 2^(n+1)

    if VECTORIZED:

        for n in range(1, RIEMANN_ITER_LIMIT):
            NK1_array[n, 0] = NK1_array[n - 1, 0] / 2
            NK1_array[n, 1:(n+1)] = (NK1_array[n-1, 0:n] + NK1_array[n-1, 1:(n+1)])/2

        # Precompute sum of above coefficients for each value of k. These will be used in Euler transform
        for k in range(0, RIEMANN_ITER_LIMIT):
            NK2_array[k] = ((-1) ** k) * np.sum(NK1_array[k:RIEMANN_ITER_LIMIT, k])

    else:
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

    delay = time.time() - t1
    print("Precomputed eta coefficients for %d iterations in %1.4f seconds " % (RIEMANN_ITER_LIMIT, delay))


def eta_zeta_scale(v):
    # Scale factor converts Dirichlet eta function to Riemann zeta function
    return 1 / (1 - 2 ** (1 - v))


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
# Can this replace the explicit log_riemann method we wrote at end of this file?
#
def riemann(s, get_array_size=False, do_eta=False, use_zero_for_nan=True, use_log=False):
    global row_count, array_zero_split, elapsed_time, entry_count

    if np.size(s) > 1:
        if s.ndim == 2:  # Somehow, I get error if I merge the comparisons with "and" and if s is single value or list

            if np.max(np.abs(s)) >= 400 and not use_log:
                # Switch to log algorithm
                r = riemann(s, get_array_size=get_array_size, do_eta=do_eta,
                                   use_zero_for_nan=use_zero_for_nan,
                                   use_log=True)

                # Compress magnitude range by 100-fold
                r = r - np.real(r) * 0.99
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
            # When doing heatmaps, s is initially a 2d ndarray, which behaves like a list of 1d ndarrays.
            # The following will call itself recursively for each item (row), and build a list of 1d arrays
            out = [0.0 if quit_computation_flag else riemann_row(x, do_eta=do_eta, use_log=use_log) for x in s]

            entry_count -= 1
            if not quit_computation_flag:
                elapsed_time = time.time() - t1
                # Convert list of 1d arrays into a single 2d ndarray
                return np.stack(out)
            else:
                elapsed_time = 0
                return None

    #
    # Either have single value, or 1d array.
    # Heatmap calculations never come here, since they have have dedicated code that takes
    # advantage of the predictable structure of the 2d array
    #
    if np.size(s) > 1.0:
        # Handle 1d array using vectorized operations (but without cached denominators)
        out_vals = riemann_row_non_negative(s, row_num=-1,
                                        do_eta=do_eta, USE_CACHED_FUNC=False)

        if use_log:
            return np.log(out_vals)

        return out_vals

    #        return [riemann(x, do_eta=do_eta, use_zero_for_nan=use_zero_for_nan) for x in s]

    # If we are here, we are computing for a single value. Use slow method, without vectorization
    # or cached denominators
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
            return -s * riemann(1 - s, do_eta=do_eta, use_log=use_log) \
                   * gamma(- s) * np.sin(-s * np.pi / 2) * (np.pi ** (s - 1)) \
                   * (1 - 2 ** (s - 1)) / (1 - 2 ** s)
        else:
            return riemann(1 - s, do_eta=do_eta, use_log=use_log) \
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

    if len(NK2_array) != RIEMANN_ITER_LIMIT:
        precompute_coeffs()

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
            if use_log:
                return np.log(cum_sum), partial_sum_index

            return cum_sum, partial_sum_index

    if use_log:
        return np.log(cum_sum)

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
    return np.sin(a)*(1+np.exp(-2*np.abs(b)))/2 + 1j*np.cos(a)*np.sign(b)*(1-np.exp(-2*np.abs(b)))/2


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
def riemann_row(s, do_eta=False, USE_CACHED_DENOM=True, use_log=False):
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

    out_vals = riemann_row_non_negative(in_vals, row_count - 1, do_eta=do_eta, USE_CACHED_FUNC=USE_CACHED_DENOM)

    if array_zero_split > 0:
        out_vals[0:array_zero_split] = np.conj(out_vals[0:array_zero_split])

    if use_log:
        # If reporting log(Riemann) then also use loggamma and logsin
        out_vals = np.log(out_vals)

        if array_zero_split > 0:
            # Use functional equation to correct results for Re[s] < 0
            if do_eta:
                out_vals[0:array_zero_split] += np.log(-2) + np.log(s0) \
                                                + loggamma(- s0) + logsin(-s0 * np.pi / 2) \
                                                + np.log(np.pi) * (s0 - 1) \
                                                + np.log(1 - np.power(2, s0 - 1)) \
                                                - np.log(1 - np.power(2, s0))
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
                                                * np.power(2,
                                                           s0)  # This can be combined with above exponential. Keep this way for readability

    if use_log:
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
    cached_powers_mag = np.empty((RIEMANN_ITER_LIMIT, len(row1)), dtype=np.complex)

    for k in range(0, RIEMANN_ITER_LIMIT):
        # Denominator magnitudes as a function of x. This array is sorted by k, then column
        cached_powers_mag[k] = np.concatenate((np.power(k + 1, -np.real(1 - s0)),
                                               np.power(k + 1, -np.real(s1))))

    col1 = np.empty(len(mesh), dtype=np.complex)
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
bases = []
powers = []
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

    if powers == []:
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
    ratios = np.power(bases, -1j * imag_diff)

    if is_odd:
        limit1 = half_s - 1
    else:
        limit1 = half_s

    # Second half of range is obtained by scalar multiplication from the first half
    powers[:, half_s:limit] = np.multiply(ratios, powers[:, 0:limit1])


#
# Calculate Riemann zeta function for a 1D ndarray (i.e. vector) of non-negative values.
#
# Uses vectorized operations and (optional) pre-computed (cached) power array.
# Note that this does not call gamma(), hence will work for very large input values
#
def riemann_row_non_negative(s,
                             row_num,  # Which row of mesh grid? 0 is bottom
                             do_eta=False,  # Do Dirichlet eta instead of Riemann zeta
                             USE_CACHED_FUNC=True,
                             is_vertical=False):  # True if evenly spaced along vertical line
    global quit_computation_flag
    global bases, powers

    if len(s) == 0:
        return []

    if USE_CACHED_FUNC:

        if USE_MATRIX_MULT:
            # Dot product is much faster than loop
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
            # Use pre-computed power table, allowing 2D vectorization, and further speedup
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
                            cum_sum[startPoint:len_s] = np.dot(NK2_power, powers[:,0:(len_s - startPoint)])
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
            # is actually very slow, and takes LONGER than the loop below.
        else:
            # Use 1D vectorization, which is MUCH faster than nothing at all.
            # There is still a loop across the k's. I tried vectorizing that dimension also, but it actually made things
            # slower - about 29 sec instead of 26 on home office computer, and 12 sec vs 8 sec on mac mini. By comparison
            # matlab did the 2D vectorized operation in about 5 seconds on home office computer, or about 5x faster.
            cum_sum = [0 + 0j] * len(s)

            for k in range(0, RIEMANN_ITER_LIMIT):
                cum_sum += NK2_array[k] * np.power(k + 1, -s)  # Calculate power the traditional (4x slower) way.

        if do_eta:
            return cum_sum
        else:
            # Scale factor converts Dirichlet eta function to Riemann zeta function
            return cum_sum / (1 - np.power(2, 1 - s))


# Return log of riemann(s). Input must be 1d array or single value, with Re[s] >= 0 for all s
# Uses vectorized operations, and for 2d arrays also uses cached denominators
def log_riemann(s, do_eta=False, is_vertical=False):
    global row_count, elapsed_time

    if type(s) == np.ndarray:
        if s.ndim == 2:
            t1 = time.time()
            precompute_denom(s)
            delay = time.time() - t1
            print("Precomputed powers in " + "{:1.4f}".format(delay) + " seconds")
            row_count = 0
            t1 = time.time()

            row_count = 0
            out = [0 if quit_computation_flag else riemann_row(x, do_eta=False, use_log=True) for x in s]

            elapsed_time = time.time() - t1
            return np.stack(out)

    # Have 1D array or single value
    if np.size(s) > 1:
        # Remove values with Re[s] < 0 for vector
        s[np.real(s) < 0] = 0
        cum_sum = riemann_row_non_negative(s, 0, do_eta=False, USE_CACHED_FUNC=False, is_vertical=is_vertical)
        return np.log(cum_sum)
    else:
        # Have single value
        s = float(s)  # np.power can't handle negative integer exponents, so make sure it is float
        cum_sum = 0 + 0j
        for k in range(0, RIEMANN_ITER_LIMIT):
            cum_sum += NK2_array[k] * np.power(k + 1, -s)

        if do_eta:
            return np.log(cum_sum)
        else:
            # Scale factor converts Dirichlet eta function to Riemann zeta function
            return np.log(cum_sum) - np.log(1 - np.power(2, 1 - s))


#
# When multiplied by gamma(s/2) * pi ^(-s/2), the result has 180-deg rotational symmetry around s = 0.5 + 0j
#
# When multiplied by an additional (1/2)s(s-1), we get the riemann xi function
#
def riemann_symmetric(s, use_log=False, is_vertical=False):
    if not use_log:

        if np.max(abs(s)) > 400:
            # Switch to log algorithm
            r = riemann_symmetric(s, use_log=True, is_vertical=is_vertical)

            # Compress magnitude range by 100-fold
            r = r - np.real(r) * 0.99
            return np.exp(r)

        r = riemann(s)
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
           + log_riemann(s, is_vertical=is_vertical) \
           + loggamma(s / 2) - s * np.log(np.pi) / 2


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

        return

    # Call vectorized version
    return gamma(s)
