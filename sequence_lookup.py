import contfrac as cf
import math
import itertools
from scipy.optimize import root_scalar


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
        return ''
    if n == 1:
        return family
    else:
        seq = []
        remainder = find_all_sequences(family, n - 1)
        for s in remainder:
            for x in family:
                r = s
                r += x
                seq.append(r)
        return seq
    

def weight(s, acc):
    seq = []
    for x in s:
        seq.append(int(x))
    w = 1
    for i in range(len(seq)):
        w *= cf.evaluate([0] + (rotate_list(seq, i) * acc))
    return w




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

                w = math.pow(lookup[n][par].weight, pow)

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
    
    sum = 0
    for w in weights[n]:
        # print(w)
        frac = math.pow(w, 2*s) / (1 - math.pow(w, 2))
        sum += frac
    # print (2 / n)
    return sum * (2 / n)


# Returns b_2k(s) for a given k and s
def b(s, k, partitions, C_n):
    if k <= 0:
        return 0
    
    sum = 0
    for j in range(k):
        for p in partitions[k][j]:
            # print(p)
            r = len(p)
            factor = math.pow(-1, r) / math.factorial(r)
            prod = 1
            for n in p:
                # print(n)
                prod *= C_n[2 * n]
            sum += (factor * prod)
    # print(sum)
    return sum


def est(s, N, partitions, weights):
    C_n =  [C(s, n, weights) for n in range(N + 1)]
    sum = 1
    for k in range(1, (N // 2) + 1):
        sum += b(s, k, partitions, C_n)
    return sum


alph = ['1', '2', '3', '4', '5', '6', '7', '8', '9']
N = 6
acc = 100
# s_N = 0.5306685908375786932401360892819672638338139332442


lookup = get_full_lookup(alph, N, acc)
# print(lookup[2])


partitions = [r_tuples(k) for k in range(0, (N // 2) + 1)]

# weights = []
# for n in range(N + 1):
#     l = lookup[n]
#     weights.append([l[s].weight for s in l])

weights = [[l[s].weight for s in l] for l in lookup]

# C_n =  [C(s_N, n, weights) for n in range(N + 1)]
# b_2k = [b(s_N, k, C_n) for k in range((N // 2) + 1)]

print(root_scalar(est, args=(N, partitions, weights), bracket=[0, 1]))

# 21:05 was working now isn't lol