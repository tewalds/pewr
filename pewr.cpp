
#define _USE_MATH_DEFINES
#include <cmath>

#include <complex>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdint.h>

#include <fftw3.h>

#define NDEBUG
#include "Array.h"
#include "time.h"

// Compile with g++ -lfftw3 -lm -fopenmp pewr.cpp -o pewr

using namespace std;
using namespace Array;

typedef std::complex<double> Complex;

typedef array2<Complex> ArrayComplex;
typedef array2<double>  ArrayDouble;
typedef array2<bool>    ArrayBool;

template <class T> std::string to_str(T a){
	std::stringstream out;
	out << a;
	return out.str();
}

struct Plane {
	int size, padding;
	double       fval;      // focal plane
	ArrayDouble  amplitude; // amplitudes, stores the initial image until converted to amplitudes
	ArrayComplex prop;      // propagation value to defocus the exit wave
	ArrayComplex ew;        // exit wave plane in the real domain
	fftw_plan    fftfwd;    // fast fourier transform in the forward direction space -> freq
	fftw_plan    fftbwd;    // fast fourier transform in the reverse direction freq -> space
	
	Plane(int _size, int _padding) : 
		size(_size), padding(_padding),
		amplitude(size, size, sizeof(double)),
		prop(     padding, padding, sizeof(Complex)),
		ew(       padding, padding, sizeof(Complex))
		{
		fftfwd = fftw_plan_dft_2d(padding, padding, reinterpret_cast<fftw_complex*>(ew()), reinterpret_cast<fftw_complex*>(ew()), FFTW_FORWARD, FFTW_MEASURE);
		fftbwd = fftw_plan_dft_2d(padding, padding, reinterpret_cast<fftw_complex*>(ew()), reinterpret_cast<fftw_complex*>(ew()), FFTW_BACKWARD, FFTW_MEASURE);
	}

	template <class T> void import(const string & name){
		ifstream ifs(name.c_str(), ios::in | ios::binary);
		T in;
		for(int x = 0; x < size; x++){
			for(int y = 0; y < size; y++){
				ifs.read( (char *) & in, sizeof(T));
				amplitude[x][y] = in;
			}
		}
		ifs.close();
	}

	void dump(const string & name){
		ofstream ofs(name.c_str(), ios::out | ios::binary);
		for(int x = 0; x < padding; x++)
			for(int y = 0; y < padding; y++)
				ofs.write((const char*)& ew[x][y], sizeof(Complex));
		ofs.close();
	}

	double mean(){
		double sum = 0;
		for(int x = 0; x < size; x++)
			for(int y = 0; y < size; y++)
				sum += amplitude[x][y];
		return sum/(size*size);
	}

	void compute_amplitudes(){
		for(int x = 0; x < size; x++)
			for(int y = 0; y < size; y++)
				amplitude[x][y] = sqrt(abs(amplitude[x][y]));
	}
};

void die(int code, const string & str){
	printf("%s\n", str.c_str());
	exit(code);
}

class PEWR {
	bool            verbose; // output the timings of each part of the algorithm
	int             size;    // size of the planes without padding
	int             padding; // size of the planes with padding
	int             nplanes; // number of planes
	int             iters;   // number of iterations
	double          lambda;  // wavelength of electrons
	double          psize;   // pixel size
	double          qmax;    // aperature size for the top hat filter
	vector<Plane *> planes;  // stack of planes
	ArrayBool       tophat;  // precompute q2 boundary
	ArrayComplex    ew;      // current best guess of the exit wave in space domain
	ArrayComplex    ewfft;   // current best guess of the exit wave in frequency domain
	string          output;  // prefix for the name of the output file
	int             outputfreq;  // output on a linear scale, ie if x=10, output 10, 20, 30, 40, 50, etc
	double          outputpower; // output on an exponential scale, ie if x=2, output 2, 4, 8, 16, 32, etc
	int             outputlast;  // output the last few iterations

	double q2(int x, int y){
		double qx = (((x + padding/2) % padding) - padding/2) / ( padding * psize);
		double qy = (((y + padding/2) % padding) - padding/2) / ( padding * psize);
		return qx*qx + qy*qy;
	}

public:
	PEWR(const string & config){
		verbose = false;
		size    = 0;
		padding = 0;
		nplanes = 0;
		qmax    = 0;
		lambda  = 0;
		psize   = 0;
		iters   = 0;
		outputfreq = 0;
		outputpower = 0;
		outputlast = 1;

		string type;
		int startiter = 1;

		ifstream ifs(config.c_str(), ifstream::in);

		string::size_type dirpos = config.find("/");
		if(dirpos != string::npos){
			string dir = config.substr(0, dirpos);
			if(chdir(dir.c_str()) == -1)
				die(1, "Cannot change directories to " + dir);
		}

		Time start;
		cout << "Parsing config file, loading data ... ";
		cout.flush();

		while(ifs.good()){
			string cmd;
			ifs >> cmd;

			if(cmd == "")
				break;
			
//			cout << cmd << endl;
			
			if(cmd == "size"){
				ifs >> size;
			}else if(cmd == "padding"){
				ifs >> padding;

				ew.Allocate(   padding, padding, sizeof(Complex));
				ewfft.Allocate(padding, padding, sizeof(Complex));
			}else if(cmd == "verbose"){
				verbose = true;
			}else if(cmd == "nplanes"){
				if(size == 0 || padding == 0)
					die(1, "padding and size must be set before nplanes");

				ifs >> nplanes;
				for(int i = 0; i < nplanes; i++)
					planes.push_back(new Plane(size, padding));
			}else if(cmd == "qmax"){
				ifs >> qmax;
			}else if(cmd == "lambda"){
				ifs >> lambda;
			}else if(cmd == "psize"){
				ifs >> psize;
			}else if(cmd == "iters"){
				ifs >> iters;
			}else if(cmd == "type"){
				ifs >> type;
			}else if(cmd == "output"){
				ifs >> output;
			}else if(cmd == "outputfreq"){
				ifs >> outputfreq;
			}else if(cmd == "outputpower"){
				ifs >> outputpower;
			}else if(cmd == "outputlast"){
				ifs >> outputlast;
			}else if(cmd == "planes"){
				if(nplanes == 0)
					die(1, "nplanes size must be set before planes");

				if(type == "")
					die(1, "type must be set before planes");

				for(int i = 0; i < nplanes; i++){
					string name;
					ifs >> name;

					if(     type == "uint8" ) planes[i]->import<uint8_t>(name);
					else if(type ==  "int8" ) planes[i]->import< int8_t>(name);
					else if(type == "uint16") planes[i]->import<uint16_t>(name);
					else if(type ==  "int16") planes[i]->import< int16_t>(name);
					else if(type == "uint32") planes[i]->import<uint32_t>(name);
					else if(type ==  "int32") planes[i]->import< int32_t>(name);
					else if(type == "float" ) planes[i]->import<float>(name);
					else if(type == "double") planes[i]->import<double>(name);
					else die(1, "Unknown type " + type);
				}
			}else if(cmd == "fvals"){
				if(nplanes == 0)
					die(1, "nplanes size must be set before fvals");

				for(int i = 0; i < nplanes; i++)
					ifs >> planes[i]->fval;
			}else if(cmd == "frange"){
				if(nplanes == 0)
					die(1, "nplanes size must be set before frange");

				int start, incr;
				ifs >> start >> incr;

				for(int i = 0; i < nplanes; i++){
					planes[i]->fval = start;
					start += incr;
				}
			}else if(cmd == "guess"){
				string name;
				ifs >> name >> startiter;

				ifstream ifguess(name.c_str(), ios::in | ios::binary);
				for(int x = 0; x < size; x++)
					for(int y = 0; y < size; y++)
						ifguess.read( (char *) & ew[x][y], sizeof(Complex));
				ifguess.close();
			}else{
				die(1, "Unknown command " + cmd);
			}
		}

		ifs.close();

		cout << "precomputing data ... ";
		cout.flush();

		//normalize
		double mean = 0;
		#pragma omp parallel for schedule(guided) reduction(+:mean)
		for(int i = 0; i < nplanes; i++)
			mean += planes[i]->mean();
		mean /= nplanes;
		#pragma omp parallel for schedule(guided)
		for(int i = 0; i < nplanes; i++)
			planes[i]->amplitude *= 1.0/mean;

		//compute amplitudes	
		#pragma omp parallel for schedule(guided)
		for(int i = 0; i < nplanes; i++)
			planes[i]->compute_amplitudes();

		//precompute q2
		double qmax2 = qmax*qmax;
		tophat.Allocate(padding, padding, sizeof(double));
		#pragma omp parallel for schedule(guided)
		for(int x = 0; x < padding; x++)
			for(int y = 0; y < padding; y++)
				tophat[x][y] = (q2(x, y) <= qmax2);

		//compute propagation arrays
		#pragma omp parallel for schedule(guided)
		for(int i = 0; i < nplanes; i++){
			for(int x = 0; x < padding; x++){
				for(int y = 0; y < padding; y++){
					double chi = M_PI * lambda * planes[i]->fval * q2(x, y);
					planes[i]->prop[x][y] = polar(1.0, -chi);
				}
			}
		}

		//setup the base fftw plans
		fftw_plan fftfwd = fftw_plan_dft_2d(padding, padding, reinterpret_cast<fftw_complex*>(ew()), reinterpret_cast<fftw_complex*>(ewfft()), FFTW_FORWARD, FFTW_MEASURE);
		fftw_plan fftbwd = fftw_plan_dft_2d(padding, padding, reinterpret_cast<fftw_complex*>(ewfft()), reinterpret_cast<fftw_complex*>(ew()), FFTW_BACKWARD, FFTW_MEASURE);

		//setup the initial approximations
		if(startiter > 1)
			#pragma omp parallel for schedule(guided)
			for(int x = 0; x < padding; x++)
				for(int y = 0; y < padding; y++)
					ew[x][y] = Complex(1, 0);

		fftw_execute(fftfwd);

		cout << "done in " << (int)((Time() - start)*1000) << " msec\n";
		cout.flush();

		double nextpoweroutput = 1;
		if(outputpower > 0)
			while(nextpoweroutput > startiter)
				nextpoweroutput *= outputpower;

		// Run iterations
		for(int iter = startiter; iter <= iters; iter++){

			Time startiter;
			cout << "Iter " << iter << " ...";
			cout.flush();

			double timedelta[8];
			for(int i = 0; i < 8; i++)
				timedelta[i] = 0;

			#pragma omp parallel for schedule(dynamic)
			for(int p = 0; p < nplanes; p++){
				Plane * plane = planes[p];

				Time time1, time2;

				// Propagate EW to each plane
				for(int x = 0; x < padding; x++){
					for(int y = 0; y < padding; y++){
						if(tophat[x][y]){
							plane->ew[x][y] = ewfft[x][y]*plane->prop[x][y];
						}else{
							plane->ew[x][y] = 0;
						}
					}
				}

				if(verbose){
					time2 = Time();
					timedelta[0] += time2 - time1;
				}

				fftw_execute(plane->fftbwd); //plane->ewfft => plane->ew

				if(verbose){
					time1 = Time();
					timedelta[1] += time1 - time2;
				}

				plane->ew *= 1.0/(padding*padding);

				if(verbose){
					time2 = Time();
					timedelta[2] += time2 - time1;
				}

				// Replace EW amplitudes
				for(int x = 0; x < size; x++)
					for(int y = 0; y < size; y++)
						plane->ew[x][y] = polar(plane->amplitude[x][y], arg(plane->ew[x][y]));

				if(verbose){
					time1 = Time();
					timedelta[3] += time1 - time2;
				}

				// Back propagate EW to zero plane, backpropagation is merged with mean
				fftw_execute(plane->fftfwd); //plane->ew => plane->ewfft

				if(verbose){
					time2 = Time();
					timedelta[4] += time2 - time1;
				}
			}
		
			Time time1, time2;

			// Backpropagate and find mean EW
			#pragma omp parallel for schedule(guided)
			for(int x = 0; x < padding; x++){
				for(int y = 0; y < padding; y++){
					if(tophat[x][y]){
						Complex mean = 0;
						for(int p = 0; p < nplanes; p++)
							mean += planes[p]->ew[x][y] * conj(planes[p]->prop[x][y]);
						ewfft[x][y] = mean / (double)nplanes;
					}else{
						ewfft[x][y] = 0;
					}
				}
			}

			if(verbose){
				time2 = Time();
				timedelta[6] += time2 - time1;
			}

			//output exit wave
			if(((outputfreq > 0 && iter % outputfreq == 0) ||
			    (outputpower > 0 && nextpoweroutput <= iter) ||
			    (iters - iter < outputlast)) && output.size() > 0){

				if(outputpower > 0)
					nextpoweroutput *= outputpower;

				fftw_execute(fftbwd); //ewfft -> ew

				ew *= 1.0/(padding*padding);

				ofstream ofs((output + "." + to_str(iter)).c_str(), ios::out | ios::binary);
				for(int x = 0; x < padding; x++)
					for(int y = 0; y < padding; y++)
						ofs.write((const char*)& ew[x][y], sizeof(Complex));
				ofs.close();
			}

			if(verbose){
				time1 = Time();
				timedelta[7] += time1 - time2;
			}

			if(verbose)
				for(int i = 0; i < 8; i++)
					cout << " " << (int)(timedelta[i]*1000);

			cout << " done in " << (int)((Time() - startiter)*1000) << " msec\n";
		}
	}
};

int main(int argc, char **argv){
	if(argc < 2)
		die(1, "Must pass the name of a config file");

	PEWR pewr(argv[1]);
}


