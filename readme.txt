Our code can also be seen on this github repo:

https://github.com/LihengGong/parallel_proj


We made great effort to make sure that the code can compile and run on both MacOS and Linux(Ubuntu)


====================================================================================================
If TBB code is to be tested, please make sure that TBB is correctly installed on your system:

https://software.intel.com/en-us/intel-tbb?_ga=2.47479911.1144381025.1556767450-348508815.1556490817&elq_cid=4863694&erpm_id=7727868
====================================================================================================



Then compiling and running the program would be easy:

====================================================================================================
- To compile and run OpenMP version of the code 
please execute openmp_macos.sh for MacOS and 
execute openmp_linux.sh for Linux
====================================================================================================


====================================================================================================
- To compile and run TBB version fo the code
please execute inteltbb_macos.sh for MacOS and
execute inteltbb_linux.sh for Linux
====================================================================================================

Excerpt of a sample output would be:

.....
.....
[100%] Linking CXX executable parallel_proj_bin
make[2]: warning:  Clock skew detected.  Your build may be incomplete.
[100%] Built target parallel_proj_bin


>>>>>>Now use 1 thread as baseline...<<<<<<


scale bunny
num of thread 1
num of thread 1

real    0m5.962s
user    0m5.916s
sys 0m0.028s


>>>>>>Now use 8 threads to parallelize...<<<<<<


scale bunny
num of thread 8

real    0m1.603s
user    0m11.636s
sys 0m0.136s


