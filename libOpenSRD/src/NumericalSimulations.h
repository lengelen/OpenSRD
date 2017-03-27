/*
	OpenSRD 1.0.0	 NumericalSimulations.h: for numerical validation, explained in accompanying paper
    Copyright (C) 2017  Lukas Engelen

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef NUMERICALSIMULATIONS_H_
#define NUMERICALSIMULATIONS_H_


#include "objectTypes.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <random>
#include <cmath>
vector<Corner> ShuffleCorners(vector<Corner> CorrectCorners, float std);
vector<Point3f> ShufflePoints(vector<Point3f> CorrectCorners, float stdd);

#endif /* NUMERICALSIMULATIONS_H_ */
