- [Compilation](#compilation)
  * [Personal computer](#personal-computer)
  * [Shared cluster](#shared-cluster)
- [Execution](#execution)
  * [Personal computer](#personal-computer-1)
  * [Shared cluster](#shared-cluster-1)
- [Help](#help)
- [Arguments](#arguments)

# simulator_metapop_hermaphrodites
  
This is a simulator of metapopulations. Demes are composed of hermaphrodites, with optional super-males (androdioecy) or super-females (gynodioecy).  
  
# Compilation
## Personal computer  
```
gcc quantiSex.c -L/usr/local/lib -lgsl -lgslcblas -lm -Wall -Wextra -Wshadow -Werror -O3 -o quantiSex
```
  
## Shared cluster  
```
LD_LIBRARY_PATH=/shared/home/croux/gsl/lib
gcc -I/shared/home/croux/gsl/include -L/shared/home/croux/gsl/lib quantiSex.c -lgsl -lgslcblas -lm -Wall -Wextra -Wshadow -Werror -O3 -o quantiSex
```

# Execution  
## Personal computer  
```
./quantiSex 100 100   100   10 0.0001 1 0.00001   100   1 0.05   0.1 1 1   0 0.1   20   123
```

## Shared cluster
```
# before the first execution
module load gcc
LD_LIBRARY_PATH=/home/roux/gsl/lib

# then simply do everytimes
./quantiSex 100 100   100   10 0.0001 1 0.00001   100   1 0.05   0.1 1 1   0 0.1   20   123
```

# Help  
```
./quantiSex
```
  
# Arguments  
1.  Number of demes (>0)  
2.  Max number of individuals per deme (>0)  
  
3.  Number of generations (>0)  
  
4.  Number of neutral loci (>=0)  
5.  Neutral mutation rate (in [0-1]))  
6.  Number of quantitative loci (>0)  
7.  Quantitative mutation rate (in [0-1])  
  
8.  Max number of offsprings per hermaphrodite (>0)  
  
9. Immigration rate (Poisson distributed; >=0)  
10. Rate of pollen dispersal  
  
11. Extinction rate (Binomialy distributed; in [0-1])  
12. Number of individuals recolonizing an extincted deme (>0)  
13. Colonization model, 'migration pool' (=0) or 'propagule pool' (=1) models  
  
14. sexualSystem is equal to 0 if autosomal, equal to 1 if XY and equal to 2 if ZW  
15. Selfing rate of hermaphrodites, fixed over time (in [0-1])  
  
16. Frequency at which statistics are written in the outfile (>0)  
  
17. Seed for the random generator (>0)  

