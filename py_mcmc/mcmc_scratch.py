import csv

#load tajD states file
fn = '/Users/boettger/postdoc/bal_sel/TajD/py_mcmc/hum_chi_tajD_states.csv'

def pop_transition_states(state_list):

    # initialize counters:
    states_total = len(state_list)
    states_on = sum(state_list)
    states_tr = [[0,0],[0,0]]

    # 1. Initialize the first state
    prev_state = state_list[0]
    states_on += state_list[0]
    states_total += 1

    # 3. Increment on_state counts for first state
    for state in state_list[2:]:
        


with open(fn, 'r') as csvfile:
    csvreader = csv.reader(csvfile)

    #throw away the header
    csvreader.next()

    prev_chr_num = 1
    a_states = []
    b_states = []

    for line in csvreader:
        chr_num, a, b = line


        if prev_chr_num == int(chr_num):
            a_states.append(int(a))
            b_states.append(int(b))
        else:
            a_tr_data = pop_transition_states(a_states)
            b_tr_data = pop_transition_states(b_states)
            prev_chr_num = int(chr_num)





    # Deal with the first line separately.
    # 1. Keep track of chromosome
    prev_chr_num = first_line.pop(0)



    # 4. Keep track of previous state for transitions
    prev_state = first_line


        chr_num = line.pop(0)

        # First, see if chr_num is the same.
        # If it has changed, let's stop for now.
        if prev_chr_num != chr_num:
            print 'Chr 1 done.'
            break

        states_total += 1
        states_on[0] += first_line[0]
        states_on[1] += first_line[1]






