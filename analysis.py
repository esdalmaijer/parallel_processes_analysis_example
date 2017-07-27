import os
from multiprocessing import cpu_count, Process, Queue

import numpy


# # # # #
# SINGLE PARTICIPANT ANALYSIS

def analyse_participants(nrs, datadir, queue):
    
    """Computes the median response time for a list of single participant.
    
    Arguments
    
    nrs         -   A list of participant numbers, for example [0, 1, 5].
    
    datadir     -   A string that point to the data directory, for example
                    "C:\\Users\\Example\\experiment\\data"
    
    queue       -   A multiprocessing.Queue instance that will be used to
                    store the computed data in.
    """
    
    # Loop through all participant numbers.
    for i in nrs:
    
        # Construct the path to this participants file.
        fpath = os.path.join(datadir, 'pp%d.txt' % (i))
        
        # Check if the path exists, and skip to the next participant if it
        # doesn't.
        if not os.path.isfile(fpath):
            continue
        
        # Load the data. This will provide us with a NumPy array with shape
        # (K, N) where K is the number of variables, and N the number of
        # trials+1. Each vector in raw, e.g. raw[0,:] will contain the data
        # for one variable (e.g. RT), for all trials, and at the first index
        # it will have the variable name. For example, raw[0,:3] could be:
        # array(['RT', '540.717892422', '540.196351231'], dtype='|S13')
        raw = numpy.loadtxt(fpath, dtype=str, unpack=True)
        
        # Next up is unpacking the raw data to a dictionary with one key for
        # every column (=variable) in the data file.
        data = {}
        # Loop through all variables in the raw data.
        for colnr in range(raw.shape[0]):
            # Get the variable name from the raw data.
            var = raw[colnr, 0]
            # Get the values for each trial.
            val = raw[colnr, 1:]
            # Save the values in the data dict. First, we attempt to convert
            # the data to a number data type. If that doesn't work, we save
            # them as string data.
            try:
                data[var] = val.astype(float)
            except:
                data[var] = val
        
        # Remove all the NaNs from the data.
        not_nan = numpy.isnan(data['RT']) == False
        all_rt = data['RT'][not_nan]

        # Compute the median response time for this participant.
        m = numpy.median(all_rt)
        sd = numpy.std(all_rt)
        
        # Put the data in the queue.
        queue.put([i, m, sd])
        
        # REMOVE THE PI COMPUTATION IF YOU WANT TO USE THIS EXAMPLE SERIOUSLY.
        # Do some silly additional computations to add to the computational
        # load per participant. In this case, let's compute pi.
        n_computations = 1000000
        s = sum(1.0/k**2 for k in range(1,n_computations+1))
        pi = (6*s)**0.5
    
    return


# # # # #
# RUN ANALYSIS

# This if statement makes sure that we only run the full analysis script when
# we attempt to run it, but not in the parallel Processes that we will spawn.
if __name__ == "__main__":
    
    import time
    from matplotlib import pyplot
    

    # FILES AND FOLDERS
    # Get the path to the directory that this file is in.
    DIR = os.path.dirname(os.path.abspath(__file__))
    # Construct the path to the data directory.
    DATADIR = os.path.join(DIR, 'data')
    # Check if the data directory exists, and throw an error if it doesn't.
    if not os.path.isdir(DATADIR):
        raise Exception("ERROR: Data directory not found!")
    # Construct a path to the output directory.
    OUTDIR = os.path.join(DIR, 'output')
    # Check if the output directory exists, and create it if it doesn't.
    if not os.path.isdir(OUTDIR):
        os.mkdir(OUTDIR)

    # ANALYSIS SETTINGS    
    # Define the number of participants.
    n_participants = 10
    # Define the number of parallel processes.
    n_processes = 5
    # Check whether we can actually use this many processes, and throw an
    # error if we don't have enough CPU cores.
    if cpu_count() < n_processes + 1:
        raise Exception("ERROR: You requested %d processes; only %d available" \
            % (n_processes, cpu_count()))
    # Compute how many participants will be run per process.
    p_per_p = n_participants // n_processes
    
    # RUN ANALYSIS.
    # Get the starting time.
    t0 = time.time()
    # Start with a Queue to hold the data in, and an empty list to hold all
    # the proceses that we will spawn.
    data_queue = Queue()
    processes = []
    # Start the requested number of Processes.
    print("\nProcessing data for %d participants in %d Processes." % \
        (n_participants, n_processes))
    for i in range(n_processes):
        # Compute the first (si) and last (ei-1) participant number.
        si = i * p_per_p
        ei = (i+1) * p_per_p
        # If this is the last Process, use it to analyse all of the remaining
        # participants.
        if i == n_processes - 1:
            ei = max(ei, n_participants)
        # Create a list of all participant numbers for this Process.
        nrs = range(si, ei)
        # Create a new Process.
        p = Process(target=analyse_participants, args=[nrs,DATADIR,data_queue])
        p.name = "analyser_%d" % (i)
        print("\tStarting Process for participants %s" % (nrs))
        p.start()
        processes.append(p)
    # Wait until all processes have finished.
    print("\nWaiting for all Processes to finish")
    for p in processes:
        p.join()
        print("Process %s is done!" % (p.name))
    # Get the ending time.
    t1 = time.time()
    # Report how quick the data processing happened.
    print("\nProcessed %d participants in %d Processes in %.3f seconds." % \
        (n_participants, n_processes, t1-t0))
    
    # COMBINE DATA.
    print("\nGathering the data.")
    # Store the data in three lists: One for the subject names, one for the
    # median RT, and one for the standard deviation.
    participants = []
    medians = []
    sds = []
    while not data_queue.empty():
        p, m, sd = data_queue.get()
        participants.append(p)
        medians.append(m)
        sds.append(sd)
    
    # PLOT THE DATA.
    print("Plotting the data.")
    # Create a new pyplot Figure with an Axis.
    fig, ax = pyplot.subplots(nrows=1, ncols=1)
    # Loop through all participants to plot all individual data.
    x = 1.05
    for i, pnr in enumerate(participants):
        # Draw this participants' data.
        ax.errorbar(x, medians[i], yerr=sds[i], fmt='o', color='green', \
            ecolor='black', alpha=0.3)
        # Annotate this participant's number.
        ax.annotate('pp%d' % (pnr), (x, medians[i]), fontsize=8, alpha=0.3)
        # Update the x-value.
        x += 0.1
    # Compute the group average and standard deviation.
    M = numpy.mean(medians)
    SD = numpy.std(medians)
    # Compute where to draw the group mean.
    x = 1.0 + (n_participants//2) * 0.1
    # Plot the group average and standard deviation.
    ax.errorbar(x, M, yerr=SD, fmt='o', color='green', ecolor='black', \
        label="group average")
    # Set a limit on the x-axis.
    ax.set_xlim([0, 2+n_participants*0.1])
    # Write the axis labels.
    ax.set_ylabel("Response time (ms)")
    # Draw the legend.
    ax.legend()
    # Save the figure.
    fig.savefig(os.path.join(OUTDIR, "group_average_RT.png"))
    # Close the figure.
    pyplot.close(fig)
    
    print("All done!")
