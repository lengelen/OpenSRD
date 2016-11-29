/*
 * ReadAndWrite.cpp
 *
 *  Created on: Nov 27, 2016
 *      Author: lengelen
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
void writeMatToFile(cv::Mat& m, const char* filename, string  outputDirectory)
{ ///write a Float-matrix to outputdirectory, using filename
    std::ofstream fout(outputDirectory+filename);

    if(!fout)
    {
        cout<<"File Not Opened"<<endl;  return;
    }

    for(int i=0; i<m.rows; i++)
    {
        for(int j=0; j<m.cols; j++)
        {
            fout<<m.at<float>(i,j)<<"\t";
        }
        fout<<endl;
    }

    fout.close();
}
void writeVecToFile(vector<Corner> features, string  file)
{///write a Float-Vector of 3D feature points to file with path name out
	std::ofstream fout(file);
    if(!fout)
    {
        cout<<"File Not Opened: "<<file<<endl;  return;
    }

    for(size_t i=0; i<features.size(); i++)
    {
              fout<< features[i].getiD()<<"\t" <<features[i].printCoordinates();
              fout <<endl;
    }
    fout.close();
}
void writeArray(vector<vector<real_1d_array> >planes, string  out)
{///write an Array of surface coefficients to file with path name out
    ofstream fout(out);

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
void writeArrayErrors(vector<vector<double> > errors, string  out)
{///write an Array of surface coefficients to file with path name out
    ofstream fout(out);

    if(!fout)
    {
        cout<<"File Not Opened"<<endl;  return;
    }


    //Concatenate vector of vectors
    	for( size_t t=0; t<errors.size();t++){
    		for(size_t i=0; i<errors[t].size(); i++)
    		fout << errors[t][i]<<endl;

    }

    fout.close();
}
void saveCameraParams( const string& filename, Mat Rotationmatrix, Mat TranslationVector)
{ ///Saves the results into file with filename
    FileStorage fs( filename, FileStorage::WRITE ); //Open filestorage file

    //Writing out the results
    fs << "Rotationmatrix" << Rotationmatrix;
    fs << "TranslationVector" << TranslationVector;
    fs.release();
	}



