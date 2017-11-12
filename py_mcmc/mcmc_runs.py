from random import random
from operator import mul
from collections import Counter

# hum_tr = [
#         # 0->0, 0->1
#         [2020606, 35101],
#         # 1->0, 1->1
#         [35101, 73093]]

# chi_tr = [
#         # 0->0, 0->1
#         [2016134, 39573],
#         # 1->0, 1->1
#         [39573, 68621]]

# gor_tr = [
#         # 0->0, 0->1
#         [2018784, 36923],
#         # 1->0, 1->1
#         [36923, 71271]]


### noHLA
# hum_tr = [
#         # 0->0, 0->1
#         [2018907, 35016],
#         # 1->0, 1->1
#         [35016, 72479]]

# chi_tr = [
#         # 0->0, 0->1
#         [2014030, 39497],
#         # 1->0, 1->1
#         [39497, 68394]]

# gor_tr = [
#         # 0->0, 0->1
#         [2016583, 36875],
#         # 1->0, 1->1
#         [36875, 71085]]

# ## these have no MHC
# hum_tr = [
#         # 0->0, 0->1
#         [2016670, 34975],
#         # 1->0, 1->1
#         [34975, 72361]]

# chi_tr = [
#         # 0->0, 0->1
#         [2011758, 39460],
#         # 1->0, 1->1
#         [39460, 68303]]

# transp_tr = [
#         # 0->0, 0->1
#         [2157020, 117],
#         # 1->0, 1->1
#         [117, 1727]]

humchi_tr = [
        # 0->0, 0->1
        [2149731, 3326],
        # 1->0, 1->1
        [3326, 2598]]

transp_tr = [
        # 0->0, 0->1
        [2157020, 117],
        # 1->0, 1->1
        [117, 1727]]


def get_tr_probs(tr_counts):
    tr_matrix = [[
            tr_counts[0][0]/float(sum(tr_counts[0])),
            tr_counts[0][1]/float(sum(tr_counts[0]))],[
            tr_counts[1][0]/float(sum(tr_counts[1])),
            tr_counts[1][1]/float(sum(tr_counts[1]))]]
    return tr_matrix

def get_st_prob(tr_counts):
    on_count = sum(tr_counts[1])
    total_count = sum(tr_counts[0])+sum(tr_counts[1])
    return on_count/float(total_count)

def get_nsites(tr_counts):
    return sum(tr_counts[0] + tr_counts[1]) + 1

#------

class MCMCTester(object):

    def __init__(self, tr_matrix, st_prob):
        self.tr_matrix = tr_matrix
        self.st_prob = st_prob

        # counter for transitions
        self.tr_counter = 0

        # randomize start state
        self.state = (random() < self.st_prob)

        # overlapped flag - for counting consecutive
        self.overlapped = 0

    def next(self):
        self.tr_counter += 1
        r = random()
        if self.tr_matrix[self.state][self.state] < r:
            # flip the state.
            self.state = not self.state
        # reset overlap if we are now off, else leave alone
        self.overlapped = self.overlapped * self.state

    def sum_states(self, n):
        on_count = 0
        for i in range(n):
            self.next()
            on_count += self.state
        return on_count

    def concur_states(self, other, n):
        on_concur_count = 0
        for i in range(n):
            self.next()
            other.next()
            on_concur_count += self.state * other.state
        return on_concur_count

    def concur_runs(self, other, n):
        on_concur_count = 0
        on_concur_runs = 0
        for i in range(n):
            self.next()
            other.next()
            concur = self.state * other.state
            on_concur_count += concur
            # only count a run if we've not previously concurred this run
            on_concur_runs += concur * (not self.overlapped)
            # if we are concurring now, or we were overlapped previously, 
            # then set overlap flag
            self.overlapped = concur or self.overlapped

        return (on_concur_count, on_concur_runs)

    def concur_many_others(self, others, n):
        on_concur_count = 0
        for i in range(n):

            # get next state for self and all others
            all_objs = [self] + others
            [obj.next() for obj in all_objs]

            # take the product of other.state for all others
            # use 'reduce() and mul()' to successively multiply each value
            # with the last, cumulatively
            # https://docs.python.org/2/library/functions.html#reduce
            on_concur_count += reduce(mul, [obj.state for obj in all_objs], 1)

        return on_concur_count

if __name__ == '__main__':

    #change these - FIRST is the one where RUNS counts
    trs = [transp_tr, humchi_tr]

    trps = [get_tr_probs(tr) for tr in trs]
    stps = [get_st_prob(tr) for tr in trs]
    mcmcs = [MCMCTester(trp, stp) for trp, stp in zip(trps, stps)]

    n_chrom = 22
    n_tests = 1000
    expected_sites = get_nsites(trs[0])

    #number of sites across all species:
    assert all([get_nsites(tr) == expected_sites for tr in trs]), \
            "Number of Sites are not equal across all species."

    #simulate for all separate chromosomes
    # run_all_chrs = lambda: sum([
    #         mcmcs[0].concur_many_others(
    #                 others= mcmcs[1:],
    #                 n= expected_sites/n_chrom) for i in range(n_chrom)])
    
    # concur runs, not states (ORDER OF MCMCS matters here)
    def run_all_chrs(n_chrom, expected_sites):
        one_chr = lambda: mcmcs[0].concur_runs(
                other= mcmcs[1],
                n= expected_sites/n_chrom)
        mcmc_chrs = [one_chr() for i in range(n_chrom)]
        (counts, runs) = zip(*mcmc_chrs)
        return "\t".join([str(sum(counts)), str(sum(runs))])

    print "\n".join(
            [run_all_chrs(n_chrom, expected_sites) for i in range(n_tests)])



#to put in background
#python mcmc_runs.py >> comb_chimphuman_transspec_runs.txt &




