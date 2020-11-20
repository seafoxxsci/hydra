#!/usr/bin/python

import sys
import os
from math import log, sqrt
import gmpy

def total_diversity(multi_species_count,outfile=sys.stdout):
    print >> outfile, \
         "Habitat\tShannon\tmax_Shannon\tEvenness\tSimpson\t"+ \
         "Chao1\tstdv(Chao1)\tVariance"
    for i in multi_species_count:
        chao1_r, chao1_stdv = chao1(multi_species_count[i])
        shannon_r = shannon(multi_species_count[i])
        max_shannon_r = max_shannon(multi_species_count[i])
        print >> outfile, \
           "%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f" % (
           i,shannon_r,
             max_shannon_r,
             shannon_r/max_shannon_r,
             simpson(multi_species_count[i]),             
             chao1_r,
             chao1_stdv,
             species_variance(multi_species_count[i],
                              len(multi_species_count[i]))
             )


def sorensen(species_1, species_2):
    """
    beta-diversity using the Sorensen formula.

    """
    # species_n is a dictionary {'species_name': number}
    n_common = 0.
    S1 = len(species_1)
    S2 = len(species_2)
    for i in species_1:
        if i in species_2:
            n_common += 1.
    return 2.*n_common/(S1+S2)

def whittaker(species_1, species_2):
    # Warning! incomplete
    S = len(species_1) + len(species_2)
    
def chao1(species_count):
    a = 0.
    b = 0.
    S1 = -1
    stdv  = 0.
    S_obs = len(species_count)
    for i in species_count:
        if i == 1: 
            a += 1
        elif i == 2:
            b += 1
    if b > 0:
        S1 = S_obs + a**2/(2*b)
        stdv = b * ( (a/(4*b))**4 + (a/b)**3 + (a/(2*b))**2 )
    return S1, stdv

def C(n,k):
    """
    Binomial n choose k
    """
    return gmpy.bincoef(n,k)

def expected_species_number(species_count,n):

    # n is the simulated sample size
    N = 0.          # number of individuals
    S = 0          # number of species
    insum = 0.
    for i in species_count:
        N += i
        if i > 0:
            S += 1.
    for i in species_count:
        insum += C((N-i),n)
#    print S, insum

    return S - 1./C(N,n) * insum

def species_variance(species_count,n):
    insum_1 = 0.
    insum_2 = 0.
    N = 0.
    S = 0.
    for i in species_count:
        N += i
        if i > 0:
            S += 1.

    N_on_n = C(N,n)
    for Ni in species_count:
        insum_1 += C(N-Ni, n) * (1. - (C(N-Ni, n) / N_on_n))
    
    for i, Ni in enumerate(species_count):
        for Nj in species_count[:i]:
            insum_2 += C(N-Ni-Nj,n) - (C(N-Ni,n)*C(N-Nj, n))/N_on_n

    return 1./N_on_n*(insum_1 + 2*insum_2)

def rarefaction(species_count):
    N = 0
    rf_vector = []
    for i in species_count:
        N += i
    for n in range(1,N+1):
        sp_var = species_variance(species_count,n)
#        print n, sp_var
        if abs(sp_var) < 1e-7:
            sp_var = 0.0
        try:
            rf_vector.append((n, 
                          expected_species_number(species_count, n),
                          sqrt(sp_var)))
        except:
            print "foobar",sp_var; raise
#        print "%d\t%.2f\t%.2f" % rf_vector[-1]
    return rf_vector

def plot_rarefaction(rarefaction_vector,outfile=sys.stdout):
    for i in rarefaction_vector:
        print >> outfile, "%d\t%.2f\t%.2f" % i

def all_rarefactions(multi_species_count):
    # multi_species_count is a dictionary containing multiple species counts
    # This function returns a dictionary of rarefaction vectors, using the
    # same keys as the ones of multi_species_count

    rarefaction_vectors = {}
    for i in multi_species_count:
        dummy_vec = [int(round(j)) for j in multi_species_count[i]]
#        rarefaction_vectors[i] = rarefaction(multi_species_count[i])
        rarefaction_vectors[i] = rarefaction(dummy_vec)
    return rarefaction_vectors

def table_rarefactions(rarefaction_vectors, outfile, sepchar="\t"):
    # Print 
    # print the titles first
    rvs = rarefaction_vectors.keys()
    rvs.sort()
    for title in rvs:
        outfile.write("N%s" % sepchar)
        for i in title[:-1]:
            outfile.write("%s" % str(i)) 
            outfile.write("-")
        outfile.write("%s%sstdv%s" % (title[-1], sepchar, sepchar))
    outfile.write(os.linesep)
    outfile.flush()

    # Now print the rows
    max_rows = max([len(rarefaction_vectors[i]) for i in rarefaction_vectors])
    for row in range(max_rows):
        for i in rvs:
            if row < len(rarefaction_vectors[i]):
                print rarefaction_vectors[i][row]
                outfile.write("%d\t%.2f\t%.2f\t" % 
                rarefaction_vectors[i][row])
            else:
                outfile.write(3*sepchar)
        outfile.write(os.linesep)
        outfile.flush()


if __name__ == '__main__':
    pass

