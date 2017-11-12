from random import random


################# human-chimp comp ###################
hum_tr = [
        # 0->0, 0->1
        [2419372, 40084],
        # 1->0, 1->1
        [40084, 89645]]

chi_tr = [
        # 0->0, 0->1
        [2503663, 42761],
        # 1->0, 1->1
        [42761, 82193]]


################# human-mus comp ###################
# hum_tr = [
#         # 0->0, 0->1
#         [1716784, 29177],
#         # 1->0, 1->1
#         [29177, 55469]]

# chi_tr = [
#         # 0->0, 0->1
#         [1719148, 24476],
#         # 1->0, 1->1
#         [24476, 62506]]

############## chimp-mus comp ##############

# ## chimp
# hum_tr = [
#         # 0->0, 0->1
#         [1725594, 32662],
#         # 1->0, 1->1
#         [32662, 52760]]

# #mus
# chi_tr = [
#         # 0->0, 0->1
#         [1731952, 24282],
#         # 1->0, 1->1
#         [24282, 63162]]





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

hum_trp = get_tr_probs(hum_tr)
chi_trp = get_tr_probs(chi_tr)

hum_st_prob = get_st_prob(hum_tr)
chi_st_prob = get_st_prob(chi_tr)


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
            other.next()%
            on_concur_count += self.state * other.state
        return on_concur_count

# quick test that it is indeed ~5%
# obj1 = MCMCTester(hum_trp, hum_st_prob)
# sum([obj1.sum_states(1000) for i in range(1000)])/float(1000)

# Testing Two Objects for Random Concurrence


if __name__ == '__main__':
    hum_mcmc = MCMCTester(hum_trp, hum_st_prob)
    chi_mcmc = MCMCTester(chi_trp, chi_st_prob)

#simulate for 22 separate chromosomes
    run_22_chrs = lambda: sum([
        hum_mcmc.concur_states(chi_mcmc, 1843679/22) for i in range(22)])

    print [run_22_chrs() for i in range(10)]


#to put in background
#python mcmc_test.py >> human_mus_full.txt &




