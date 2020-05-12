//
// Created by Gonzalo Lera Romero on 12/05/2020.
//

#ifndef TDTSPTW_RELAXATION_ALGORITHM_H
#define TDTSPTW_RELAXATION_ALGORITHM_H

#include "goc/goc.h"

namespace tdtsptw
{
class RelaxationAlgorithm
{
public:
	enum Type { NGLTI, NGLTD };
	enum Direction { Forward, Backward, Bidirectional };
	Type type;
	Direction direction;


};
}

#endif //TDTSPTW_RELAXATION_ALGORITHM_H
