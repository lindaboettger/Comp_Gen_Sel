from random import random
from operator import mul


hum_tr = [
        # 0->0, 0->1
        [986121, 17409],
        # 1->0, 1->1
        [17409, 154765]]


gor_tr = [
        # 0->0, 0->1
        [987486, 18842],
        # 1->0, 1->1
        [18842, 150534]]

chi_tr = [
        # 0->0, 0->1
        [989134, 17506],
        # 1->0, 1->1
        [17506, 151558]]

mus_tr = [
        # 0->0, 0->1
        [985813, 18856],
        # 1->0, 1->1
        [18856, 152179]]


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
        self.tr_counter = 0
        self.state = (random() < self.st_prob)

    def next(self):
        self.tr_counter += 1
        r = random()
        last_state = self.state
        if self.tr_matrix[self.state][self.state] < r:
            self.state = not self.state

    def sum_states(self, n):
        on_count = 0
        for i in range(n):
            self.next()
            on_count += obj1.state
        return on_count

    def concur_states(self, other, n):
        on_concur_count = 0
        for i in range(n):
            self.next()
            other.next()
            on_concur_count += self.state * other.state
        return on_concur_count

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

    #change these
    trs = [hum_tr, gor_tr, chi_tr]

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
    run_all_chrs = lambda: sum([
            mcmcs[0].concur_many_others(
                    others= mcmcs[1:],
                    n= expected_sites/n_chrom) for i in range(n_chrom)])

    print [run_all_chrs() for i in range(n_tests)]



#to put in background
#python mcmc_4species.py >> four_species_binwiden.txt &




