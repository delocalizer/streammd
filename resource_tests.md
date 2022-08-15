## Testing mem, cpu and walltime usage

### picard

```
Job Id: 17824602.hpcpbs02
    Job_Name = md
    resources_used.cpupercent = 239
    resources_used.cput = 01:22:27
    resources_used.mem = 30948696kb
    resources_used.walltime = 01:17:34
```

### streammd

For all runs the following streammd configuration was used:
```
DEFAULT_FPRATE = 1e-6
DEFAULT_NITEMS = 1e9
DEFAULT_NWORKERS = 10
```

Note that the memory usage for runs that weren't on cgroup enabled node are
incorrect as they sum the child process shared memory values. That's all
runs below except #1.

#### run 1

4 samtools input reader threads; 4 samtools output writer threads.
```
samtools view -@ 4 -h scratch/140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.qname.bam | python3 src/streammd/markdups.py |samtools view -@ 4 -O bam > scratch/140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.qname.streammd_samtools_threads.bam
```

```
Job Id: 17940747.hpcpbs02
    Job_Name = streammd_samtools_threads
    resources_used.cpupercent = 1036
    resources_used.cput = 04:54:58
    resources_used.mem = 8388608kb
    resources_used.walltime = 00:27:27
```

#### run 2

One samtools input reader thread; 4 samtools output writer threads.
```
samtools view -h scratch/140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.qname.bam | python3 src/streammd/markdups.py |samtools view -@ 4 -O bam > scratch/140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.qname.streammd_samtools_threads_2.bam
```

```
Job Id: 17940997.hpcpbs02
    Job_Name = streammd_samtools_threads_2
    resources_used.cpupercent = 937
    resources_used.cput = 03:51:35
    resources_used.mem = 38849744kb
    resources_used.walltime = 00:25:22
```

#### run 3

Eight samtools input reader threads; 4 samtools output writer threads.
```
samtools view -@ 8 -h scratch/140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.qname.bam | python3 src/streammd/markdups.py |samtools view -@ 4 -O bam > scratch/140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.qname.streammd_samtools_threads_3.bam
```

```
Job Id: 17940999.hpcpbs02
    Job_Name = streammd_samtools_threads_3
    resources_used.cpupercent = 1044
    resources_used.cput = 04:08:56
    resources_used.mem = 38871888kb
    resources_used.walltime = 00:25:35
```

#### run 4

Eight samtools input reader threads; one output writer thread

```
samtools view -@ 8 -h scratch/140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.qname.bam | python3 src/streammd/markdups.py |samtools view -O bam > scratch/140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.qname.streammd_samtools_threads_4.bam
```

```
Job Id: 17940998.hpcpbs02
    Job_Name = streammd_samtools_threads_4
    resources_used.cpupercent = 288
    resources_used.cput = 04:12:14
    resources_used.mem = 38861064kb
    resources_used.walltime = 01:43:58

```

#### run 5

Eight samtools input reader threads; no output writer (direct to SAM)

```
samtools view -@ 8 -h scratch/140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.qname.bam | python3 src/streammd/markdups.py > scratch/140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.qname.streammd_samtools_threads_5.sam
```

```
Job Id: 17941320.hpcpbs02
    Job_Name = streammd_samtools_threads_5
    resources_used.cpupercent = 1208
    resources_used.cput = 02:46:13
    resources_used.mem = 38853776kb
    resources_used.walltime = 00:14:16
```

#### run 6

One samtools input reader thread; no output writer (direct to SAM)

```
samtools view -h scratch/140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.qname.bam | python3 src/streammd/markdups.py > scratch/140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.qname.streammd_samtools_threads_6.sam
```

```
Job Id: 17941644.hpcpbs02
    Job_Name = streammd_samtools_threads_6
    resources_used.cpupercent = 1123
    resources_used.cput = 02:51:42
    resources_used.mem = 38832620kb
    resources_used.walltime = 00:16:07
```

Aside: `top` during this one looked like this —

```
113356 conradL   20   0 7328928   3.4g   3.3g R  94.9  0.9   6:45.31 python3                     
113358 conradL   20   0 7328928   3.4g   3.3g R  94.9  0.9   6:45.97 python3                     
113359 conradL   20   0 7328928   3.4g   3.3g R  94.9  0.9   6:47.22 python3                     
113362 conradL   20   0 7328928   3.4g   3.3g R  94.9  0.9   6:46.67 python3                     
113363 conradL   20   0 7328928   3.4g   3.3g R  94.9  0.9   6:48.27 python3                     
113364 conradL   20   0 7328928   3.4g   3.3g R  94.9  0.9   6:47.91 python3                     
113355 conradL   20   0 7328928   3.4g   3.3g R  94.2  0.9   6:45.39 python3                     
113357 conradL   20   0 7328928   3.4g   3.3g R  93.5  0.9   6:45.70 python3                     
113360 conradL   20   0 7328928   3.4g   3.3g R  93.5  0.9   6:46.23 python3                     
113361 conradL   20   0 7328928   3.4g   3.3g R  93.5  0.9   6:45.46 python3                     
113346 conradL   20   0   88380   3428   2456 R  91.3  0.0   5:54.28 samtools                    
113348 conradL   20   0  488848  57428   1512 R  75.4  0.0   5:57.16 python3   
```
so all the worker threads pegged at ~ 100% CPU; yays.

#### run 7

One samtools input reader thread; 8 output writer threads

```
samtools view -h scratch/140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.qname.bam | python3 src/streammd/markdups.py |samtools view -@ 8 -O bam > scratch/140924_EXTERN_0588_AC1DTBACXX.lane_2.nobc.qname.streammd.bam
```

```
Job Id: 17941786.hpcpbs02
    Job_Name = streammd_samtools_threads
    resources_used.cpupercent = 1708
    resources_used.cput = 03:45:00
    resources_used.mem = 38854188kb
    resources_used.walltime = 00:13:28
```

### streammd summary

All the following used 10 worker threads and `n=1e9, p=1e-6`.

| input resources | output resources | cpupercent | walltime | run    |
| --------------- | ---------------- | ---------- | -------- | ---    |
| 1 thread        | 8 threads to bam | 1708       | 00:13:28 |  7     |
| 8 threads       | 0 threads to sam | 1208       | 00:14:16 |  5     |
| 1 thread        | 0 threads to sam | 1123       | 00:16:07 |  6     |
| 1 thread        | 4 threads to bam | 937        | 00:25:22 |  2     |
| 4 threads       | 4 threads to bam | 1036       | 00:27:27 |  1     |
| 8 threads       | 4 threads to bam | 1044       | 00:25:35 |  3     |
| [picard]        | [picard]         | 239        | 01:17:34 |[picard]|
| 8 threads       | 1 threads to bam | 288        | 01:43:58 |  4     |





Conclusions:
 * With 10 worker threads, speed is limited only by how fast we can process the
   outputs. Fastest results are writing directly to SAM or using many threads
   in `samtools view -O bam`. Note that 8 output threads is essentially twice
   as fast as 4.
 * Additional threads for the input reader makes essentially no difference.
   Compare runs 2,1,3 that all used 4 output threads but 1,4 and 8 input threads
   — they all came in around the 26 minute mark and indeed the 1 input thread run
   (run 2) was marginally the fastest.
 * With single output thread the workers are very under-used (cpupercent) and
   walltime is 8 times longer than with 8 output threads.
 * Summary: streammd itself is _fast_.
