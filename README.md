RQA for very long data series
=============================

General
-------

Programme to calculate the main measures of the recurrence quantification
analysis (RQA). The calculation is performed in an efficient way directly
on the time series, without preceeding calculation of the recurrence plot
and using multithreading (based on openMP).

It allows to calculate the RQA measures from single-column files (like
from time series) and multi-column files (like multi-dimensional data), 
and additionally time delay embedding can be applied.

Available recurrence definitions: Euclidean, maximum, and minimum norm, 
as well as order patterns.

Files
-----

- `rqa_omp.cpp`: C++ code for efficient, multithreaded RQA calculation
- `lorenz.dat`: example data file representing the three coordinates of
the Lorenz system

Compilation
-----------

- Compile without multithread support:  `g++ -O3 -o rqa_omp rqa_omp.cpp`
- Compile for multithread support:  `g++ -O3 -fopenmp -o rqa_omp rqa_omp.cpp`

Usage
-----

`rqa [options]`

|Option|Explanation|
|-|-|
|`-i <string>`  |data filename (input)|
|`-o <string>`  |filename RQA measures (output)|
|`-p <string>`  |filename histogramme diagonal line lengths (output)|
|`-q <string>`  |filename histogramme vertical line lengths (output)|
|`-n <string>`  |distance norm (`EUCLIDEAN`, `MAX`, `MIN`, `OP`), default=`EUCLIDEAN`|
|`-m <number>`  |  embedding dimension, default=1|
|`-t <number>`  |embedding delay, default=1|
|`-e <number>`  |threshold, default=1|
|`-l <number>`  |l_min, default=2|
|`-v <number>`  |v_min, default=2|
|`-w <number>`  |Theiler window, default=1|
|`-s`           |silent (no messages displayed)|
|`-h`           |print this help text|

__Examples:__

1. Calculate RQA of 3-columns Lorenz using threshold 10 and Euclidean norm:
```
./rqa_omp -i lorenz.dat -e 10 -n EUCLIDEAN
```

2. Calculate RQA of a sine signal using order patterns with 3 dimensions
and delay 10, and store results in file rqa_results.dat
```
./rqa_omp -i sin.dat -m 3 -t 10 -n OP -o rqa_results.dat
```

RQA Measures
------------

- RR
- DET
- DET/RR (only available in stdout, not in results file)
- L_max
- DIV (only available in stdout, not in results file)
- L_mean
- L_entr
- LAM
- LAM/DET (only available in stdout, not in results file)
- TT
- V_max
- V_entr (only available in stdout, not in results file)

Input file
----------

Use a simple ASCII file, structured as single-column or multi-column data. Each row represents one time point; if multiple columns are present, they correspond to different variables measured at the same time.

Output file
-----------

Results can be exported as RQA measures (option `-o`) and as histogrammes of
the diagonal and vertical line length distributions (options `-p` and `-q`).
The RQA results file has eight columns corresponding to the following RQA
measures:

1. RR
2. DET
3. L_max
4. L_mean
5. L_entr
6. LAM
7. TT
8. V_max


Copyright
---------

Norbert Marwan\
Potsdam Institute for Climate Impact Research\
4/2009

License: GPLv3
