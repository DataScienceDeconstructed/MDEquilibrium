
from ComputationalEquilibriums import ReferenceDistribution
from scipy.stats import chi2_contingency
import numpy as np


if __name__ == "__main__":
    # start out with empty distributions, a negative index so that it goes to 0 at the start, the testing significance level,
    # and a tracking variable for how many times in a row the test has failed.
    dist = [ReferenceDistribution(_type="Binary", _reference=0.0, _dist=[0, 0]), ReferenceDistribution(_type="Binary", _reference=0.0, _dist=[0, 0])]
    index = -1
    significance_level = 0.05
    sim_track = 0;
    info_lag = []

    #grab the number of particles and the name of the experiment from the save file
    parts = None
    name = None
    with open(r"1st10frames.txt", 'r') as fp:
        for i, line in enumerate(fp):
            if i == 0:
                parts = int(line.strip())
            if i == 1:
                name = line.strip()
            if i > 1:
                break
    #process the file
    with open(r"1st10frames.txt", 'r') as fp:
        for i, line in enumerate(fp):
            split_line = line.strip().split("\t") # split the file line into its components

            # if we have 2 distributions saved, compare them
            # there is a line for each particle plus a line for the number of particles and name of experiment.
            # we have to go through 2 saved steps to run our test. that's why we have 2 *(parts+2))
            if i % (2*(parts + 2)) == 0:
                if i > 0:
                    #process distributions for Chi squared metric
                    #run this for the first and last info sets.
                    # iterative sets just show that the density growth is not exceptional.
                    # probably need to handle it in 2^N steps.
                    # code below needs refactor to handle multiple lag deltas.
                    info_lag.append([dist[0].Distribution, dist[1].Distribution])

                    # reset the distributions so that processing can continue.
                    dist = [ReferenceDistribution(_type="Binary", _reference=0.0, _dist=[0, 0]),
                        ReferenceDistribution(_type="Binary", _reference=0.0, _dist=[0, 0])]

            #swap the index each time a new saved frame is available
            if split_line[0] == name:
                index = (index + 1) % 2

            #grab and record the highest z value
            if split_line[0] == '1': # this value will always exist even on break lines with the number of particles and name of exp.
                if float(split_line[3]) > dist[index].ReferenceValue: # check to see if z value is higher than current max
                    dist[index].update_reference(float(split_line[3]))

            #identify the NP inside and outside the brush
            if split_line[0] == '2': # looking at NPs now
                dist[index].update_distribution(float(split_line[3]))


    for i in range(0, len(info_lag)):
        for j in [ 2**c for c in range( np.log2(len(info_lag)/2) )]:
            stat, p, dof, arr = chi2_contingency([info_lag[i][0], info_lag[i][1]])
            if p <= significance_level:
                print(p, 'Reject NULL HYPOTHESIS. distributions too different')
                sim_track = 0
            else:
                print(p, 'ACCEPT NULL HYPOTHESIS. distributions very similar',info_lag[i][0], info_lag[i][1])
                sim_track += 1
                if sim_track > 3:
                # the 3 is in a way arbitrary, but it means the distributions have been fairly similar for 3 times in a row.
                    break