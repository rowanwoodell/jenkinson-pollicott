import math
import mpmath as mp
import itertools
import contfrac as cf
from scipy.optimize import root_scalar
# import matplotlib.pyplot as plt

# Helper function
def rotate_list(lst, n):
    return lst[n:] + lst[:n]


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


# Returns all sequences of length n with entries in the array 'family'.
def find_all_sequences(family, n):
    if n <= 0:
        return []
    if n == 1:
        return [[x] for x in family]
    else:
        seq = []
        remainder = find_all_sequences(family, n - 1)
        for s in remainder:
            for x in family:
                r = s.copy()
                r.append(x)
                seq.append(r)
        return seq


# Takes an array of sequences and returns one which contains only one 
# sequence of each cyclic permutation class, and also only those 
# which cannot be constructed by repeating shorter sequences
def consolidate(seq):
    con_seq = seq.copy()
    for s in seq:
        # if s[:-1] != s[1:]:
            if s in con_seq:
                for i in range(1, len(s)):
                    try:
                        r = rotate_list(s, i)
                        # print(s, i, r)
                        con_seq.remove(r)
                    except:
                        pass
        
    return con_seq


# Returns an array of pairs consisting of the weights w_i along with 
# the multiplicity of that weight
def find_weights(all_seq, n, acc, all_weights):
    weights = []
    for s in all_seq[n]:
        weights.append([weight(s, acc), n])

    for m in range(1, n):
        if n % m == 0:
            divisor_weights = all_weights[m]
            for w in divisor_weights:
                weights.append([math.pow(w[0], n/m), m])

    return weights


# Returs weights with multiplicities for all sequences
def find_all_weights(all_seq, acc):
    all_weights = []
    for n in range(len(all_seq)):
        all_weights.append(find_weights(all_seq, n, acc, all_weights))
    return all_weights


# Returns the weight of a given sequence (acc determines the accuracy)
def weight(s, acc):
    w = 1
    for i in range(len(s)):
        w *= cf.evaluate([0] + (rotate_list(s, i) * acc))
    return w


# Returns C_n(s) for a given n and s
def C(s, n, weights):
    sum = 0
    for w in weights[n]:
        # print(w)
        frac = math.pow(w[0], 2*s) / (1 - math.pow(w[0], 2))
        frac *= w[1]
        sum += frac
    # print (2 / n)
    return sum * (2 / n)


# Returns b_2k(s) for a given k and s
def b(s, k, weights, partitions):
    sum = 0
    for j in range(k):
        for p in partitions[k][j]:
            # print(p)
            r = len(p)
            factor = math.pow(-1, r) / mp.fac(r)
            prod = 1
            for n in p:
                # print(n)
                prod *= C(s, 2*n, weights)
            sum += (factor * prod)
    # print(sum)
    return sum


# Returns 1 + sum_n=1^N b_n(s). This is the function we
# aim to find the root of.
def pressure_est(s, N, weights, partitions):
    sum = 1
    for i in range(1, N+1):
        sum += b(s, i, weights, partitions)
    return sum


alph = [1, 2]
N = 10

all_cs = [consolidate(find_all_sequences(alph, x)) for x in range(0, N + 1)]

# s = find_all_sequences([1, 2], 2)
# print(s)

# cs = consolidate_sequences(s)
# print(cs)

# for the same level of accuracy, this seems to need to be larger when n (sequence length) is smaller
# accuracy_factor = 100
# for s in all_cs[4]:
#     w_i = 1
#     for i in range(len(s)):
#         w_i *= cf.evaluate([0] + (rotate_list(s, i) * accuracy_factor))
#     print(w_i)

a_w = find_all_weights(all_cs, 100)

p = [r_tuples(k) for k in range(0, (N // 2) + 1)]

# for x in p:
#     print(x)

# for i in range(len(all_cs)):
#     # print(all_cs[i])
#     print(a_w[i])

# print(pressure_est(0.1, 2, a_w))

def pressure(s, N=N // 2, weights=a_w, partitions=p):
    return pressure_est(s, N, weights, partitions)

mp.plot(pressure, [0, 1])

print(root_scalar(pressure_est, args=(N // 2, a_w, p), bracket=[0.5, 1]))

# print(b(0.5, 1, a_w, p))