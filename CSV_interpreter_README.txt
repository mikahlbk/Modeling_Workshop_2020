The CSV_interpreter program is made to work with 
the batchmaker.cpp program, so compile by

g++ batchmaker.cpp -o batchGenerator.out
g++ CSV_interpreter.cpp -o CSV_interpreter.out

and run it with

./CSV_interpreter.out <csv file name>.csv


In: CSV file organized as follows:

-batchmaker_flag_1,flag_value_1, ...,-batchmaker_flag_M,flag_value_M,-test,NAME
-parameter_flag_1,...,-parameter_flag_P
x11,x12,...,x1P
...
xN1,XN2,...,XNP

Out: N simulations submitted with parameter vectors (Xk1,...,XkP)
submitted via slurm as NAME_1, ... , NAME_N. 

List of batchmaker_flags (Optional marked with *)
-test <NAME>
-hours <##> 
-cores <##>
-nodes <##>
-partition <short/batch/intel>
*-mem <#> (# is in GB, default 2)
*-nodes <#> (defaults to 1).
