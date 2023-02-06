#define NODE n
#define VERTICAL_LEVEL nz

#define DO_NODES_AND_VERTICAL_LEVELS \
    do NODE = 1, myDim_nod2d; \
        do VERTICAL_LEVEL = 1, nl

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
