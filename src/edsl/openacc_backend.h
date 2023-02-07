#define DO_NODES_AND_VERTICAL_LEVELS(NODE_ITERATOR, NODE_LB, NODE_UB, VERTICAL_LEVEL_ITERATOR, VERTICAL_LEVEL_LB, VERTICAL_LEVEL_UB) \
    do NODE_ITERATOR = NODE_LB, NODE_UB; \
        do VERTICAL_LEVEL_ITERATOR = VERTICAL_LEVEL_LB, VERTICAL_LEVEL_UB

#define END_NODES_AND_VERTICAL_LEVELS \
        end do; \
    end do

#define OPERATE_AT(dimension, target) \
    if (dimension == target) then

#define END_OPERATE_AT \
    end if

#define OPERATE_ON_RANGE(dimension, lower_bound, upper_bound) \
    if(lower_bound <= dimension .and. dimension <= upper_bound) then

#define END_OPERATE_RANGE \
    end if
