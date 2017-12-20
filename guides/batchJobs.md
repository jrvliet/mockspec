
# Guide to Swinburne Queue

The best way to run *mockspec* is by submitting jobs to the g2 queue. This guide
walks through the basics of submitting to the queue. 

## Different Machines
The g2 computer at Swinburne has four queues that the jobs can be submitted to:

1. sstar
2. gstar
3. largemem
4. manygpu

The details of each of these different machines can be found at
https://supercomputing.swin.edu.au/job-queue/. In summary:

| Machine  | Cores | Mem  | Max walltime |
|--------- |-------|------|--------------|
| sstar    | 16    | 64   | 7 days       |
| gstar    | 12    | 48   | 7 days       |
| largemem | 32    | 512  | 2 days       |
| manygpu  | 12    | 48   | 7 days       |

For *mockspec* it is best to use largemem. If rates, cellfinder, and TPCF are
note being run, it is safe to run *mockspec* with either gstar or sstar. In
general sstar is slightly more efficient than gstar. There is no reason to ever
use manygpu since *mockspec* does not work with GPUs. 


## PBS 

The job to be submitted to the queue is described in a batch file in the
Portable Batch Script (PBS) format. This is standard for supercomputers the
world over and dozens of guides can be found online for it. A typical script for
*mockspec* is shown below:

```
#!/bin/csh

# queue name
#PBS -q largemem

# resource requests
#PBS -l nodes=1:ppn=4
#PBS -l walltime=00:48:00:00
#PBS -l pmem=30g

# Name the job
#PBS -N jobName

# Output files
#PBS -m ae
#PBS -j oe

# load modules
module load intel

module list
echo Deploying job to CPUs...
cat $PBS_NODEFILE

date
cd /path/to/galaxy/expansion
/path/to/python /path/to/mockspec/
date
```

All PBS commands start with #PBS, followed by a flag then the options. The
options used here are summarized below.

| Flag | Description |
|-----------------|-------------|
| -q   | Sets the queue to use |
| -l   | Defines the resources to request from the system |
| -N   | Name of the job |
| -m   | Mail alert at (b)eginning, (e)nd, and (a)bortion of job |
| -j   | Job log files (o)utput and (e)rror |


The most important part is the resource requests. The nodes arguement details
how many nodes are needed for the code to run, which is just one. The ppn is the
number of processors per node and should be set equal to the number of
processors outlined in mockspec.config. The default is four. The walltime
arguement specifies how long the code will run as measured by a clock on the
wall. If the code takes longer than the specified wall time, the system will
terminate the run. If the code takes significantly less time than the specified 
wall time, nothing bad happen, but it wastes computer resources. Try to keep the
amount of time requested accurate to how long the code will actually take to
run. The format of the arguement is *days:hours:minutes:seconds*.

Everything after the # load module line is code that is executed when the job
runs. Any standard linux command can be put here. The module load intel line
loads common intel packages to make sure the code executes properly. The module
list and cat $PBS_NODEFILE lines output the state of the system, which is useful
for debugging if something goes wrong. 

The output of the run is written to a log file formatted as <jobName>.o<jobID>
in the directory the job was submitted. 


## Submitting

To actually submit the job to the queue, use the qsub command as:

```
qsub <name_of_batch_file>
```

If successful, the output will be:
```
<jobID>.pbs.hpc.swin.edu.au
```

You can submit up to 100 jobs at a time.


## Status

To check the status of your queue, use the qstat command. If run without any
arguements, it will output the status of all jobs currently running, regardless
of who submitted it or which queue is being used. A more useful command is:

```
qstat -u <userID>
```

where <userID> is your username on the system. This will only output the status
of jobs that you submitted. An example output of this is:

| Field | Example Value | Description|
|-------|---------------|------------|
| Job ID | 6483016.pbs.hpc.swin.edu.au | ID of the Job |
| Username | jvander | Name of the user that submitted the job |
| Queue | largemem | Name of the machine being used |
| Jobname | dwarf-mock | Name of the job in the batch file by the -N flag |
| SessID | 4063 | Session ID |
| NDS | 1 | Number of nodes being used |
| TSK | 4 | Number of processors being used |
| Req'd Memory | -- | Amount of memory being used |
| Req'd Time | 30:00 | Walltime requested in batch file |
| S | R | Status (Q=queued, R=running) |
| Elap Time | 01:48 | Time job has been running


## Logs
The standard output of the job is stored in a file named <jobName>.o<jobID>
placed in the directory that the job was submitted from. 







