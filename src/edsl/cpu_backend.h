#define NODE n
#define VERTICAL_LEVEL nz

#define DO_NODES_AND_VERTICAL_LEVELS do NODE = 1, myDim_nod2d
#define END_NODES_AND_VERTICAL_LEVELS end do

#define OPERATE_AT(dimension, target) dimension = target
#define END_OPERATE_AT

#define OPERATE_ON_RANGE(dimension, lower_bound, upper_bound) do dimension = lower_bound, upper_bound
#define END_OPERATE_RANGE end do
