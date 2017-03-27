/*
	OpenSRD 1.0.0	 NumericalSimulations.cpp: for numerical validation, explained in accompanying paper
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
#include "NumericalSimulations.h"

template <class Type> class Shuffle
{ // Gives a new corner, (u,v) coordinate shifted with std
    private:
        float stdd;
    public:
            // Constructor
        Shuffle (float i_std) {

        	 stdd=i_std;
            }

            // The function call
            Type operator ( ) ( Corner& initial ) const
            {
            	Point3f coords=initial.getCoords();
            	std::normal_distribution<float> dist(0,stdd);

            	std::random_device rd;
            	std::mt19937 gen(rd());
            	float shiftx=dist(gen);
            	float shifty=dist(gen);
            	Point3f shuffled=Point3f(coords.x+shiftx, coords.y+shifty, 1);
            	return Corner(shuffled, initial.getiD());
            }
        };
template <class Type> class ShuffleFeatures
{ // Gives a new Point3f, (u,v) coordinate shifted with std
    private:
        float stdd;
    public:
            // Constructor
        ShuffleFeatures (float i_std) {

        	 stdd=i_std;
            }

            // The function call
            Type operator ( ) ( Point3f& initial ) const
            {
            	std::normal_distribution<float> dist(0,stdd);

            	std::random_device rd;
            	std::mt19937 gen(rd());
            	float shiftx=dist(gen);
            	float shifty=dist(gen);
            	return Point3f(initial.x+shiftx, initial.y+shifty, 0);
            }
        };
vector<Corner> ShuffleCorners(vector<Corner> CorrectCorners, float stdd)
{ // Shuffle list of corner-objects based on given standard deviation stdd -> N(0,stdd)

	vector<Corner> ShuffledCorners(CorrectCorners);
	transform(CorrectCorners.begin(), CorrectCorners.end(),ShuffledCorners.begin(), Shuffle<Corner> (stdd));
	return ShuffledCorners;
}
vector<Point3f> ShufflePoints(vector<Point3f> CorrectCorners, float stdd)
{ // Shuffle list of corner-objects based on given standard deviation stdd -> N(0,stdd)

	vector<Point3f> ShuffledCorners(CorrectCorners);
	transform(CorrectCorners.begin(), CorrectCorners.end(),ShuffledCorners.begin(), ShuffleFeatures<Point3f> (stdd));
	return ShuffledCorners;
}

