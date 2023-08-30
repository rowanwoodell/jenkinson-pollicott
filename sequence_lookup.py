import contfrac as cf
import math
import itertools
from scipy.optimize import root_scalar
from mpmath import mp
import time

mp.dps = 52


class SequenceLookup:
    def __init__(self, parent, power):
        self.parent = parent
        self.power = power
        self.weight = None

    def setParent(self, p):
        self.parent = p

    def setWeight(self, w):
        self.weight = w


def rotate_list(lst, n):
    return lst[n:] + lst[:n]


def find_all_sequences(family, n):
    if n <= 0:
        return ()
    if n == 1:
        return [tuple([x]) for x in family]
    else:
        seq = []
        remainder = find_all_sequences(family, n - 1)
        for s in remainder:
            for x in family:
                r = list(s)
                r.append(x)
                seq.append(tuple(r))
        return seq


def weight(s, acc):
    # seq = []
    # for x in s:
    #     seq.append(int(x))
    seq = list(s)
    w = mp.mpf(1)
    for i in range(len(seq)):
        # w *= cf.evaluate([0] + (rotate_list(seq, i) * acc))
        # w *= my_cf_eval(rotate_list(seq, i) * acc)
        w = mp.fmul(w, precise_cf_evaluate([0] + rotate_list(seq, i) * acc))
    # mp.nprint(w)
    return w


def precise_cf_evaluate(cont_frac):
    """Computes the floating point value of a finite continued fraction
    representation.

    That is the value of ``c[0] + 1/(c[1] + 1/(c[2] + 1/(c[3] + ...)))``
    for an input ``c``.

    Example:
        For an input of ``[2,3,4]`` is ``2 + 1/(3 + 1/4) = 30/13`` expressed as
        2.3076923076923075.

    Args:
        cont_frac (Iterable[Union[int, float]]): representation of a continued
            fraction as iterable of numbers.

    Returns:
        float: the value of the continued fraction.
    """
    if len(cont_frac) == 0:
        return 0
    fraction = 0
    for i, coefficient in enumerate(reversed(cont_frac[1:])):
        fraction = mp.fdiv(1, (mp.fadd(coefficient, fraction)))
    return mp.fadd(cont_frac[0], fraction)


def my_cf_eval(cont_frac):
    p_continuants = [0, 1, 0]
    q_continuants = [1, 0, 1]
    n = len(cont_frac)

    p_n = p(n, cont_frac, p_continuants)
    p_n_minus_one = p_continuants[n + 1]

    q_n = q(n, cont_frac, q_continuants)
    q_n_minus_one = q_continuants[n + 1]



    # print(p_n, p_n_minus_one)
    # print(q_n, q_n_minus_one)

    # print(p_n - q_n_minus_one)
    # print(math.sqrt(math.pow(q_n_minus_one - p_n, 2) + 4 * q_n * p_n_minus_one))
    # print(2 * q_n)

    numerator = p_n - q_n_minus_one + math.sqrt(math.pow(q_n_minus_one - p_n, 2) + 4 * q_n * p_n_minus_one)
    denominator = 2 * q_n

    return(numerator / denominator)


def p(n, cont_frac, p_continuants):#, continuants=[0, 1, 0]):
    m = n + 2

    if m <= len(p_continuants) - 1:
        return p_continuants[m]
    
    a_n = cont_frac[n - 1]
    p_n = a_n * p(n - 1, cont_frac, p_continuants) + p(n - 2, cont_frac, p_continuants)

    p_continuants.append(p_n)

    return p_n


def q(n, cont_frac, q_continuants):
    m = n + 2

    if m <= len(q_continuants) - 1:
        return q_continuants[m]
    
    a_n = cont_frac[n - 1]
    q_n = a_n * q(n - 1, cont_frac, q_continuants) + q(n - 2, cont_frac, q_continuants)

    q_continuants.append(q_n)

    return q_n

# for i in range(5):
#     print(find_all_sequences(alph, i))

# for i in range(4):
#     print(rotate_list(find_all_sequences(alph, 4)[1], i))

def get_cyclic_perm_classes(seq, n):
    cyclic_perm_classes = []

    for s in seq:
        break_flag = False
        for c in cyclic_perm_classes:
            for i in range(n):
                r = rotate_list(s, i)
                if r == c[0]:
                    c.append(s)
                    break_flag = True
                    break
                
            if break_flag:
                break
            
        if break_flag:
            continue
        
        cyclic_perm_classes.append([s])

    return cyclic_perm_classes


def initialise_lookup(alph, n):

    seq = find_all_sequences(alph, n)

    seq_parents = []
    for s in seq:
        for i in range(1, n + 1):
            r = rotate_list(s, i)
            if s == r:
                # print(s, r)
                seq_parents.append(SequenceLookup(s[:i], int(n / i)))
                break

    seq_lookup = dict(zip(seq, seq_parents))

    cyclic_perm_classes = get_cyclic_perm_classes(seq, n)

    for s in seq_lookup:
        for c in cyclic_perm_classes:
                if s in c:
                    slice = n // seq_lookup[s].power
                    seq_lookup[s].setParent(c[0][:slice])
                    break

    return seq_lookup


def get_full_lookup(alph, N, acc):
    lookup =  [initialise_lookup(alph, n) for n in range(N + 1)]

    for l in lookup:
        for s in l:
            if l[s].power == 1:
                w = weight(l[s].parent, acc)
            else:
                par = l[s].parent
                pow = l[s].power
                n = len(par)
                
                w = mp.power(lookup[n][par].weight, pow)
            
            l[s].setWeight(w)
    
    return lookup


# From https://jeromekelleher.net/generating-integer-partitions.html
def accel_asc(n):
    a = [0 for i in range(n + 1)]
    k = 1
    y = n - 1
    while k != 0:
        x = a[k - 1] + 1
        k -= 1
        while 2 * x <= y:
            a[k] = x
            y -= x
            k += 1
        l = k + 1
        while x <= y:
            a[k] = x
            a[l] = y
            yield a[:k + 2]
            x += 1
            y -= 1
        a[k] = x + y
        y = x + y - 1
        yield a[:k + 1]


# Returns a list containing k lists of lists of integers.
# The n-th list consists of all partitions of k of length n.
def r_tuples(k):
    if k == 0:
        return None
    
    partitions = list(accel_asc(k))
    partitions.sort(key=lambda p:len(p))

    permuted_partitions = []
    for p in partitions:
        permuted_partitions.append(list(itertools.permutations(p)))

    consolidated_partitions = []
    for p in permuted_partitions:
        consolidated_partitions.append(list(dict.fromkeys(p)))

    grouped_partitions = [[] for i in range(k)]

    i = 0
    for p in consolidated_partitions:
        for q in p:
            if len(q) != i+1:
                i += 1
            grouped_partitions[i].append(list(q))

    return grouped_partitions


# Returns C_n(s) for a given n and s
def C(s, n, weights):
    if n % 2 == 1 or n <= 0:
        return 0
    
    sum = mp.mpf(0)
    for w in weights:
        # print(w)
        frac = mp.fdiv(mp.power(w, mp.fmul(2 ,s)), (1 - mp.power(w, 2)))
        sum = mp.fadd(sum, frac)
    # print (2 / n)
    return mp.fmul(sum, mp.fdiv(2, n))


# Returns b_2k(s) for a given k and s
def b(s, k, partitions, C_n):
    if k <= 0:
        return 0
    
    sum = mp.mpf(0)
    for j in range(k):
        for p in partitions[k][j]:
            # print(p)
            r = len(p)
            factor = mp.fdiv(mp.power(-1, r), mp.factorial(r))
            prod = mp.mpf(1)
            for n in p:
                # print(n)
                prod = mp.fmul(prod, C_n[2*n])
            sum = mp.fadd(sum, (mp.fmul(factor, prod)))
    # print(sum)
    return sum




# print(p(4, [1, 5, 2, 2]))
# print(q(4, [1, 5, 2, 2]))

# print(my_cf_eval([1, 1] * 20))


alph = [mp.mpf(2), mp.mpf(5)]
# N = 10
acc = 200
# s_N = 0.5306685908375786932401360892819672638338139332442


print("family=", alph)
print("acc =", acc)

for N in range(2, 17, 2):
    time_start = time.time()

    lookup = get_full_lookup(alph, N, acc)

    partitions = [r_tuples(k) for k in range(0, (N // 2) + 1)]

    weights = [[l[s].weight for s in l] for l in lookup]
    
    def est(s, N=N, partitions=partitions, weights=weights):
        C_n =  [C(s, n, weights[n]) for n in range(N + 1)]
        sum = mp.mpf(1)
        for k in range(1, (N // 2) + 1):
            sum = mp.fadd(sum, b(s, k, partitions, C_n))
        return sum
    
    print("N =", N)
    mp.nprint(mp.findroot(est, (0.5, 1), solver="ridder"), mp.dps)

    time_end = time.time()

    print("time:", time_end - time_start, "seconds")