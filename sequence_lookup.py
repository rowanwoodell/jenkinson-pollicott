import itertools
from mpmath import mp
import time


class SequenceLookup:
    '''
    Class to store information about sequences in order to find all weights 
    without computing them all.
    '''
    def __init__(self, atom, freq):
        self.atom = atom
        self.freq = freq
        self.weight = None

    def setAtom(self, a):
        self.atom = a

    def setWeight(self, w):
        self.weight = w

    def isAtomic(self):
        if self.freq == 1:
            return True
        else:
            return False


def rotate_list(lst, n):
    '''
    Returns a list with the same entries as lst, cycled n places to the left.
    '''
    return lst[n:] + lst[:n]


def find_all_sequences(family, n):
    '''
    Returns all sequences of length n with entries in the list family.
    '''
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
    '''
    Returns the weight of the sequence s.
    '''
    seq = list(s)
    w = mp.mpf(1)

    # Cycle seq, find the value of the purely periodic cf whose coefficients 
    # are the cycled seq, and multiply w by this value. Repeat for each cycle.
    for i in range(len(seq)):
        w = mp.fmul(w, precise_cf_evaluate(rotate_list(seq, i) * acc))
    return w


def precise_cf_evaluate(cont_frac):
    '''
    Returns the real number with continued fraction coefficients given by the 
    list cont_frac. Adapted from the cf_evalute function in the contfrac 
    library in order to support mpmath operations. Also assumes the leading
    coefficient is a zero and has been omitted, i.e. the value is between 0 
    and 1 (as is always the case for us).
    '''
    if len(cont_frac) == 0:
        return 0
    
    fraction = 0
    for coefficient in reversed(cont_frac):
        fraction = mp.fdiv(1, (mp.fadd(coefficient, fraction)))
    return fraction


def get_cyclic_perm_classes(seq, n):
    '''
    Returns the cyclic permutation classes for the sequences of length n in the
    list seq.
    '''
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


def initialise_lookup(family, n):
    '''
    Returns a list of SequenceLookup objects, one for each sequence of length 
    n. Each SequenceLookup has its atom and freq fields correctly assigned,
    but no weight assigned, as this is handled by the function get_full_lookup.
    '''
    seq = find_all_sequences(family, n)

    # For each sequence of length n, find its atom and create a SequenceLookup
    # object with the atom and the frequency of the sequence.
    # A potential speedup here would be to break when i > n/2, since the 
    # frequency must divide n.
    seq_atoms = []
    for s in seq:
        for i in range(1, n + 1):
            r = rotate_list(s, i)
            # We cycle the sequence until the cycled version is equal to the
            # uncycled version; the sequence is atomic iff we cycle all the way
            # through the sequence, else we know the atomic part and create a
            # SequenceLookup object with the relevant information.
            if s == r:
                atom = s[:i]
                freq = n // i
                seq_atoms.append(SequenceLookup(atom, freq))
                break
    
    # Create a dictionary which assigns SequenceLookup objects to their
    # original sequences.
    seq_lookup = dict(zip(seq, seq_atoms))

    # Assign the same atom to all sequences in the same cyclic permutation 
    # class.
    cyclic_perm_classes = get_cyclic_perm_classes(seq, n)
    for s in seq_lookup:
        for c in cyclic_perm_classes:
            if s in c:
                period = n // seq_lookup[s].freq
                seq_lookup[s].setAtom(c[0][:period])
                break

    return seq_lookup


def get_full_lookup(family, N, acc):
    '''
    Returns a list of SequenceLookup objects, one for each sequence of length 
    m with entries in family, with 1 <= m <= n. Each SequenceLookup has all of
    its fields correctly assigned.
    '''
    lookup =  [initialise_lookup(family, n) for n in range(N + 1)]

    # Assign weights to elements of lookup by the following logic:
    # if the sequence is atomic, calculate its weight from scratch;
    # else, use its freq and the weight of its atom to find its weight. 
    for l in lookup:
        for s in l:
            if l[s].isAtomic():
                w = weight(l[s].atom, acc)
            else:
                atm = l[s].atom
                frq = l[s].freq
                n = len(atm)
                w = mp.power(lookup[n][atm].weight, frq)
            
            l[s].setWeight(w)
    
    return lookup


def accel_asc(n):
    '''
    Returns a list containing all partitions of n.
    Code from https://jeromekelleher.net/generating-integer-partitions.html
    '''    
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


def find_all_compositions(k):
    '''
    Returns a list containing k lists of lists of integers. The n-th list 
    consists of all compositions of k of length n.
    '''
    if k == 0:
        return None
    
    # Find all partitions of k and sort them by length.
    partitions = list(accel_asc(k))
    partitions.sort(key=lambda p:len(p))

    # Permute all partitions to get all compositions.
    compositions = []
    for p in partitions:
        compositions.append(list(itertools.permutations(p)))

    # Remove duplicate compositions (weird Python way).
    consolidated_compositions = []
    for c in compositions:
        consolidated_compositions.append(list(dict.fromkeys(c)))

    # Group partitions of the same length into lists.
    grouped_partitions = [[] for i in range(k)]
    i = 0
    for c in consolidated_compositions:
        for q in c:
            if len(q) != i+1:
                i += 1
            grouped_partitions[i].append(list(q))

    return grouped_partitions


def C(s, n, weights):
    '''
    Returns C_n(s) for a given n and s.
    '''
    if n % 2 == 1 or n <= 0:
        return 0
    
    sum = mp.mpf(0)
    for w in weights:
        frac = mp.fdiv(mp.power(w, mp.fmul(2 ,s)), (1 - mp.power(w, 2)))
        sum = mp.fadd(sum, frac)
    return mp.fmul(sum, mp.fdiv(2, n))


def b(k, compositions, C_n):
    '''
    Returns b_2k(s) for a given k and precomputed partitions and C_n(s).
    '''
    if k <= 0:
        return 0
    
    sum = mp.mpf(0)
    for j in range(k):
        for p in compositions[k][j]:
            r = len(p)
            factor = mp.fdiv(mp.power(-1, r), mp.factorial(r))
            prod = mp.mpf(1)
            for n in p:
                prod = mp.fmul(prod, C_n[2*n])
            sum = mp.fadd(sum, (mp.fmul(factor, prod)))
    return sum


# Set the decimal precision for the mpmath package.
mp.dps = 52

# The subset of the natural numbers which we restict our coefficients to.
family = [mp.mpf(2), mp.mpf(5)]

# Accuracy parameter - this controls how accurate our continued fraction
# calculations are. Higher is more accurate but slower. A value of 200 gives
# no error (compared to JP2001) when family = [1, 2], so we take this to be a
# safe value for all sets.
acc = 200


print("family=", family)
print("acc =", acc)

for N in range(2, 17, 2):
    time_start = time.time()

    lookup = get_full_lookup(family, N, acc)

    compositions = [find_all_compositions(k) for k in range(0, (N // 2) + 1)]

    weights = [[l[s].weight for s in l] for l in lookup]

    def est(s, N=N, compositions=compositions, weights=weights):
        '''
        Returns an estimate of the function D(s, 1). In particular, it is the
        degree-N truncation of the sum.
        '''
        C_n =  [C(s, n, weights[n]) for n in range(N + 1)]
        sum = mp.mpf(1)
        for k in range(1, (N // 2) + 1):
            sum = mp.fadd(sum, b(k, compositions, C_n))
        return sum

    print("N =", N)

    # mpmath function to find the root of est in the interval (0, 1).
    mp.nprint(mp.findroot(est, (0, 1), solver="ridder", verbose=False), mp.dps)

    time_end = time.time()

    print("time:", time_end - time_start, "seconds")