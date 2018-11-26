/*
 * Ray.cpp
 *
 *  Created on: Nov 24, 2018
 *      Author: Matt
 */

#include "Ray.h"

/**
 * Constructor
 */
Ray::Ray(VectorMath::point origin, VectorMath::point direction)
{
   // Initialize Ray
    _origin = origin;
    _direction = direction;
}

/**
 * Destructor
 */
Ray::~Ray()
{

}

/**
 * getOrigin - Gets the Origin of the Ray
 */
VectorMath::point Ray::getOrigin()
{
   return _origin;
}

/**
 * getDirection - Gets the Direction of the Ray
 */
VectorMath::point Ray::getDirection()
{
    return _direction;
}
