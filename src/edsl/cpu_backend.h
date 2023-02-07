#define DO_NODES_AND_VERTICAL_LEVELS(NODE_ITERATOR, NODE_LB, NODE_UB, VERTICAL_LEVEL_ITERATOR, VERTICAL_LEVEL_LB, VERTICAL_LEVEL_UB) \
    do NODE_ITERATOR = NODE_LB, NODE_UB
#define END_NODES_AND_VERTICAL_LEVELS \
    end do

#define OPERATE_AT(dimension, target) \
    dimension = target
#define END_OPERATE_AT

#define OPERATE_ON_RANGE(dimension, lower_bound, upper_bound) \
    do dimension = lower_bound, upper_bound
#define END_OPERATE_RANGE \
    end do
