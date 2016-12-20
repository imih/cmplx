for i in {1..160}
do
 mpiexec -quiet -n 11 ./benchmark_main -n $i -d
 mpiexec -quiet -n 11 ./benchmark_main -n $i -s
 mpiexec -quiet -n 11 ./benchmark_main -n $i -m
 mpiexec -quiet -n 11 ./benchmark_main -n $i -f

 mpiexec -quiet -n 11 ./benchmark_main -n $i -d -b
 mpiexec -quiet -n 11 ./benchmark_main -n $i -s -b 
 mpiexec -quiet -n 11 ./benchmark_main -n $i -m -b
 mpiexec -quiet -n 11 ./benchmark_main -n $i -f -b
done
