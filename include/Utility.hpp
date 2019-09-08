#ifndef __UTILITY_HPP_
#define __UTILITY_HPP_

#include <vector>
#include <string>
#include <fstream>
#include <iostream>

#include <sys/types.h>  // directory check
#include <sys/stat.h>
#include <dirent.h>

#include <boost/algorithm/string.hpp>

/*
std::vector<double> BreakCSVLine(std::string line, std::string delimeter = " ") {
    std::vector<double> d_words;
    int ind = 0;

}*/

std::vector<std::vector<double>> ReadCsv(std::string fileName, std::string delimeter = " ,\t")
{
	std::ifstream file(fileName);
	assert(file.is_open());

	std::vector<std::vector<double>> dataList;

	std::string line = "";
	// Iterate through each line and split the content using delimeter
	std::size_t n_row = 0;
	std::size_t n_col = 0;
	while (getline(file, line))
	{
	    if(line[0] == '#') {continue;}
		std::vector<std::string> vec;
		//line = "1111 1.000 2.0 3.0 -4.0";
		boost::split(vec, line, boost::is_any_of(delimeter));
		//std::cout << line << std::endl;
		//for(int i = 0; i < vec.size(); ++i) {
        //    std::cout << vec[i] << std::endl;
		//}
		if( n_row == 0 ) {
            n_col = vec.size();
            dataList.resize(n_col);
		} else {
            assert(vec.size() == n_col);
            for(int i = 0; i < n_col; ++i) {
                dataList[i].push_back(std::stod(vec[i]));
            }
		}
		n_row++;
	}
	// Close the File
	file.close();

	return dataList;
}

void GetIncidentField_Experiment(std::vector<double>& t, std::vector<double>& e_inc, bool normalized = false) {
    std::vector<std::vector<double>> data = ReadCsv("out/EO_sampling_data.csv");
    auto& t_samples = data[0];
    auto& e_inc_samples = data[1];

    std::size_t n_pts = t_samples.size();
    t.resize(n_pts);
    e_inc.resize(n_pts);
    double e_max = 0.0;
    for(std::size_t i = 0; i < n_pts; ++i) {
        t[i] = t_samples[i] * 1.0e-12;
        e_inc[i] = e_inc_samples[i];

        if(std::abs(e_inc[i]) > e_max) {
            e_max = e_inc[i];
        }
        //std::cout << "t: " << t_samples[i] << "  , E: " << e_inc_samples[i] << std::endl;
    }

    if(normalized) {
        for(std::size_t i = 0; i < n_pts; ++i) {
            e_inc[i] /= e_max;
        }
    }
}

void GetFourierTransform(std::vector<double>& t, std::vector<double>& e_t,
                         double f_0, double f_1, std::size_t n_f,
                         std::vector<double>& f, std::vector<std::complex<double>>& e_f) {
    std::size_t n_t = t.size();
    assert(e_t.size() == n_t);
    f.resize(n_f);
    e_f.resize(n_f);
    std::complex<double> _j(0.0, 1.0);
    for(std::size_t i = 0; i < n_f; ++i) {
        f[i] = f_0 + (double)i/(n_f - 1) * (f_1 - f_0);
        e_f[i] = 0.0;

        double w_i = 2.0*M_PI*f[i];
        for(std::size_t j = 0; j < n_t; ++j) {
            e_f[i] += e_t[j]*std::exp(-_j*w_i*t[j]);
        }
        e_f[i] *= (t[1] - t[0]);
    }
}

void CreateFolderIfItDoesNotExists(std::string folder, bool delete_content = false) {
    struct stat info;

    if( stat( folder.c_str(), &info ) != 0 ) {
        //printf( "cannot access %s\n", folder.c_str() ); // no such directory

        int dir_err = mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (-1 == dir_err) {
            printf("Error creating directory %s!\n", folder.c_str());
        } else {
            printf("Created directory %s!\n", folder.c_str());
        }
    } else if( info.st_mode & S_IFDIR ) {
        printf( "%s is a directory\n", folder.c_str() );
        //delete
        if(delete_content) {
            std::string syscommand = std::string("rm ") + folder + "/*";
            std::system( syscommand.c_str() );
            std::cout << "content delete!" << std::endl;
        }
    }
}


#endif // __UTILITY_HPP_

