from numpy import array, zeros, log, sqrt
from math import factorial

sam_1 = [[23, 64, 14, 0, 0, 3, 1]]  # This is sample 'A' where each element = one observation of species at index i

# These samples are worked examples from https://www.ohio.edu/plantbio/staff/mccarthy/dendro/LEC5.pdf
sam_2 = array([100, 50, 30, 20, 1])  # Expected shannon(sam_2) ~ 1.201, simpson(sam_2) ~ 0.338 and inv_simpson(sam_2)
# ~2.96
sam_3 = array([5, 5, 5, 5, 5])  # Expected brillouin(sam_3) ~ 1.362

### Read in these samples from a single file

def counts(indices, result=None):
    """
    :param indices: Given a vector of feature indices,
    :param result: and an empty 1-dimensional array to fill,
    :return: an array of counts per species i
    """
    if result is None:
        result = zeros(max(indices)+1)
    for i in indices:
        result[i] += 1
    return result


def singletons(counts):
    """
    :param counts: Given an array of counts per species i,
    :return: the number of features that only occur once (singletons)
    """
    return sum(counts == 1)


def doubletons(counts):
    """
    :param counts: Given an array of counts per species i,
    :return: the number of features that only occur twice (doubletons)
    """
    return sum(counts == 2)


def observed_species(counts):
    """
    :param counts: Given an array of counts per species i,
    :return: the number of distinct features
    """
    return sum(counts != 0)


def margalef(counts):
    """
    Formula: Magurran AE. (2004) "Measuring Biological Diversity". Blackwell Publishing, Ltd.
    Original citation: Clifford and Stephenson (1975)

    :param counts: Given an array of counts per species i,
    :return: Margalef's index D = (S - 1) / ln(N), where S = number of observed species and N = total number of
    observations (features)
    """
    o = observed_species(counts)
    n = sum(counts)
    return (o-1) / log(n)


def menhinick(counts):
    """
    Formula: Magurran AE. (2004) "Measuring Biological Diversity". Blackwell Publishing, Ltd.
    Original citation: Whittaker (1977)

    :param counts: Given an array of counts per species i,
    :return: Menhinick's index D = S / sqrt(N), where S = number of observed species and N = total number of
    observations (features)
    """
    o = observed_species(counts)
    n = sum(counts)
    return o/sqrt(n)


def simpson(counts):
    """
    Formula: Simpson EH. (1949) "Measurement of Diversity"
    Gives the approximate probability that any two species sampled are actually the same (sum of squares of
    probabilities)

    :param counts: Given an array of counts per species i,
    :return: Simpson's lambda L = sigma S (Pi**2), where Pi = proportional abundance of species i and N = total number
    of observations (features)
    """
    n = sum(counts)
    p = counts/n
    return sum(p**2)


def gini_simpson(counts):
    """
    Gives the probability that any two species sampled are actually of different types (1 - sum of squares of
    probabilities)

    :param counts: Given an array of counts per species i,
    :return: Gini-Simpson index = 1 - Simpson's lambda L (also known as the probability of an interspecific encounter,
    or PIE)
    """
    return 1 - simpson(counts)


def inv_simpson(counts):
    """
    :param counts: Given an array of counts per species i,
    :return: The inverse of Simpson's lambda L (also known as the effective number of species, or ENS)
    """
    return 1/simpson(counts)


def simpson_e(counts):
    """
    Formula: Morris EK et al. (2014) "Choosing and using diversity indices: insights for ecological applications from
    the German Biodiversity Exploratories". Ecol Evol, 4(18), 3514–3524.

    :param counts: Given an array of counts per species i,
    :return: Simpson's evenness E = ENS / S, where ENS = inverse of Simpson's lambda L and S = number of observed
    species
    """
    o = observed_species(counts)
    return inv_simpson(counts) / o


def shannon(counts):
    """
    Formula: Magurran AE. (2004) "Measuring Biological Diversity". Blackwell Publishing, Ltd.
    Original citation: Shannon CE. (1948) "A Mathematical Theory of Communication".

    :param counts: Given an array of counts per species i,
    :return: Shannon's entropy H' = -sigma S (Pi * ln(Pi))) in bits, where S = number of species
    """
    s = sum(counts)
    p = counts/s
    nonzero_p = p[p.nonzero()]  # sum(nonzero_p) should approximate 1.
    return -sum(nonzero_p * log(nonzero_p))


def pielou(counts):
    """
    Formula: Magurran AE. (2004) "Measuring Biological Diversity". Blackwell Publishing, Ltd.
    Original citation: Pielou (1969)

    :param counts: Given an array of counts per species i,
    :return: Pielou's evenness J' = H' / H'(max), where H' = Shannon's entropy in bits and H'(max) = log(s), where all
    species s have equal abundances.
    """
    s = sum(counts)
    h_max = log(s)
    return shannon(counts) / h_max


def berger_parker(counts):
    """
    Original citation: Berger WH and Parker FL. (1970) "Diversity of Planktonic Foraminifera in Deep-Sea Sediments".

    :param counts: Given an array of counts per species i,
    :return: Berger-Parker's dominance Dd = Nmax / N, where Nmax = counts of the most abundant species i and N = total
    number of observations (features)
    """
    n_max = max(counts)
    n = sum(counts)
    return n_max/n


def mcintosh(counts):
    """
    Original citation: McIntosh (1967)

    :param counts: Given an array of counts per species i,
    :return: McIntosh's dominance D = (N - U) / N - sqrt(N), where N = total number of observations (features) and U =
    sqrt(sigma counts**2), which represents the Euclidean distance of the community from the origin
    """
    u = sqrt(sum(counts**2))
    n = sum(counts)
    return (n-u) / (n-sqrt(n))


def brillouin(counts):
    """
    Formula: Zar JH (2010) "Biostatistical analysis, 5th ed." Pearson Publishing.
    Original citation: Pielou (1975)

    :param counts: Given an array of counts per species i,
    :return: Brillouin's diversity index Hb = (log(N!) - sigma k (log s!)) / N, where N = total number of observations
    (features) and s = counts of species k, not including species with zero counts because log(0) = inf...
    """
    nonzero_c = counts[counts.nonzero()]
    n = sum(counts)
    fact_n = float(factorial(n))

    c = []
    for i in nonzero_c:
        c.append(factorial(i))
    fact_s = array(c)
    return (log(fact_n) - sum(log(fact_s))) / n


# TODO: ace(), chao_1(), chao_2(), faith(), fisher(), goods(), heip(), lladser(), renyi(), strong(), kempton_taylor(),
#  robbins(), michaelis_menten(), rarefaction()