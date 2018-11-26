/*
 * Ray.h
 *
 *  Created on: Nov 24, 2018
 *      Author: Matt
 */

#ifndef _RAY_H_
#define _RAY_H_

#include "VectorMath.h"

class Ray
{
public:

    /**
     * Constructor
     */
    Ray(VectorMath::point origin, VectorMath::point direction);

    /**
     * Destructor
     */
    ~Ray();

    /**
     * getOrigin - Gets the Origin of the Ray
     */
    VectorMath::point getOrigin();

    /**
     * getDirection - Gets the Direction of the RAy
     */
    VectorMath::point getDirection();

private:
    /** Origin of the Ray */
    VectorMath::point _origin;

    /** Direction of the Ray */
    VectorMath::point _direction;
};

#endif /* _RAY_H_ */
