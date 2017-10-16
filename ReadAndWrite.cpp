/*
	OpenSRD 1.0.0	
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
#include "ReadAndWrite.h"

//Reading functions
bool readStringList(vector<string>& l, string dir )
{///finds all images in directory and stores their names in vector of strings
    l.clear();
    string filepath;
    DIR *dp;
    struct dirent *dirp;
    struct stat filestat;

    dp = opendir( dir.c_str() );
    if (dp == NULL)
    {
        cout << "Error opening " << dir << endl;
        return false;
    }

    while ((dirp = readdir( dp )))
    {
        filepath = dir + "/" + dirp->d_name;

        // If the file is a directory (or is in some way invalid) we'll skip it
        if (stat( filepath.c_str(), &filestat )) continue;
        if (S_ISDIR( filestat.st_mode ))         continue;
        l.push_back((filepath));

    }

    closedir( dp );
    std::sort( l.begin(), l.end() );
    cout <<"Number of images loaded: "<<l.size()<<endl;
    return true;
}

// Write functions
void writeVecToFile(vector<Corner> features, string  path)
{///write a Float-Vector of corner points to file with path name path
	std::ofstream fout(path);
    if(!fout)
    {
        cout<<"File Not Opened: "<<path<<endl;  return;
    }

    for(size_t i=0; i<features.size(); i++)
    {
              fout<< features[i].getiD()<<"\t" <<features[i].printCoordinates();
              fout <<endl;
    }
    fout.close();
}
void writeCentersToFile(vector<vector<Point2f>> features, string  path)
{///write a Float-Vector of corner points to file with path name path
	std::ofstream fout(path);
    if(!fout)
    {
        cout<<"File Not Opened: "<<path<<endl;  return;
    }

    for(size_t t=0; t<features.size(); t++)
    {for(size_t i=0; i<features[t].size(); i++)
	{
              fout<< features[t][i].x<<"\t" <<features[t][i].y;
              fout <<endl;
    }}
    fout.close();
}
void writePointsToFile(vector<Point3f> features,  string  file, string outputDirectory)
{///write a Float-Vector of 3D feature points to file with path name out
	std::ofstream fout(outputDirectory+file);
    if(!fout)
    {
        cout<<"File Not Opened: "<<file<<endl;  return;
    }

    for(size_t i=0; i<features.size(); i++)
    {
              fout<< features[i].x<<"\t" <<features[i].y<<"\t" <<features[i].z;
              fout <<endl;
    }
    fout.close();

}

void writeArray(vector<vector<real_1d_array> >planes,  string  file, string outputDirectory)
{///write an Array of surface coefficients to file with path name out
    ofstream fout(outputDirectory+file);

	if(!fout)
	{
		cout<<"File Not Opened"<<endl;  return;
	}

	//Concatenate vector of vectors
	for( size_t t=0; t<planes.size();t++){
		for(size_t i=0; i<planes[t].size(); i++)
			{
    		fout <<planes[t][i].tostring(5).c_str()<<endl;
    		}
	}

    fout.close();
}
void saveCameraParams( const string& file, string outputDirectory, Mat Rotationmatrix, Mat TranslationVector)
{ ///Saves the results into file with filename
    FileStorage fs( outputDirectory+file, FileStorage::WRITE ); //Open filestorage file

    //Writing out the results
    fs << "Rotationmatrix" << Rotationmatrix;
    fs << "TranslationVector" << TranslationVector;
    fs.release();
	}



