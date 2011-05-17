
#define _USE_MATH_DEFINES
#include <cmath>

#include <complex>
#include <string>
#include <vector>
#include <iostream>
#include <stdint.h>

#include "Array.h"
#include "fftw++.h"
#include "time.h"

// Compile with g++ pewr.cpp fftw++.cc -lfftw3

using namespace std;
using namespace Array;
using namespace fftwpp;

typedef array2<Complex> ArrayComplex;
typedef array2<double>  ArrayDouble;


struct Plane {
	int size, padding;
	double       fval;      // focal plane
	ArrayDouble  image;     // initial image
	ArrayDouble  amplitude; // initial amplitude
	ArrayComplex prop;      // propagation value to defocus the exit wave
	ArrayComplex backprop;  // back propagation value to refocus the exit wave
	ArrayComplex ew;        // exit wave plane in the real domain
	ArrayComplex ewfft;     // exit wave plane in the frequency domain
	fft2d        fftfwd;    // fast fourier transform in the forward direction space -> freq
	fft2d        fftbwd;    // fast fourier transform in the reverse direction freq -> space
	
	Plane(int _size, int _padding) : 
		size(_size), padding(_padding),
		image(    padding, padding, sizeof(double)),
		amplitude(padding, padding, sizeof(double)),
		prop(     padding, padding, sizeof(Complex)),
		backprop( padding, padding, sizeof(Complex)),
		ew(       padding, padding, sizeof(Complex)),
		ewfft(    padding, padding, sizeof(Complex)),
		fftfwd(-1, ew, ewfft), 
		fftbwd( 1, ewfft, ew) 
		{
	}

	template <class T> void import(const string & name){
		ifstream ifs(name.c_str(), ios::in | ios::binary);
		T in;
		for(int x = 0; x < size; x++){
			for(int y = 0; y < size; y++){
				ifs >> in;
//				ifs.read( &in, sizeof(T));

				image[x][y] = in;
			}
		}
		ifs.close();
	}

	double mean(){
		double sum = 0;
		for(int x = 0; x < size; x++)
			for(int y = 0; y < size; y++)
				sum += image[x][y];
		return sum/(size*size);
	}

	void normalize(double divisor){
		double factor = 1.0/divisor;
		for(int x = 0; x < size; x++)
			for(int y = 0; y < size; y++)
				image[x][y] *= factor;
	}

	void compute_amplitudes(){
		for(int x = 0; x < size; x++)
			for(int y = 0; y < size; y++)
				amplitude[x][y] = sqrt(abs(image[x][y]));
	}
};

void die(int code, const string & str){
	printf("%s\n", str.c_str());
	exit(code);
}

class PEWR {
	int             size;    // size of the image without padding
	int             padding; // size of the planes with padding
	int             nplanes; // number of planes
	int             iters;   // number of iterations
	double          lambda;  // wavelength of electrons
	double          psize;   // pixel size
	double          qmax;    // aperature size for the top hat filter
	vector<Plane *> planes;  // stack of planes
	ArrayComplex    ew;      // current best guess of the exit wave in space domain
	ArrayComplex    ewfft;   // current best guess of the exit wave in frequency domain

	double q2(int x, int y){
		double qx = -0.5/psize + x/(padding*psize);
		double qy = -0.5/psize + y/(padding*psize);
		return qx*qx + qy*qy;
	}

public:
	PEWR(const string & config){
		size    = 0;
		padding = 0;
		nplanes = 0;
		qmax    = 0;
		lambda  = 0;
		psize   = 0;
		iters   = 0;

		string type;

		ifstream ifs(config.c_str(), ifstream::in);

		Time start;
		cout << "Loading config file ... ";

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
					die(1, "nplanes size must be set before planes");

				for(int i = 0; i < nplanes; i++)
					ifs >> planes[i]->fval;
			}else{
				die(1, "Unknown command " + cmd);
			}
		}

		ifs.close();

		cout << "normalizing data ... ";

		//normalize
		double mean = 0;
		for(int i = 0; i < nplanes; i++)
			mean += planes[i]->mean();
		mean /= nplanes;
		for(int i = 0; i < nplanes; i++)
			planes[i]->normalize(mean);

		//compute amplitudes	
		for(int i = 0; i < nplanes; i++)
			planes[i]->compute_amplitudes();

		//compute propagation arrays
		for(int i = 0; i < nplanes; i++){
			for(int x = 0; x < size; x++){
				for(int y = 0; y < size; y++){
					double chi = M_PI * lambda * planes[i]->fval * q2(x, y);
					planes[i]->prop[x][y] = polar(1.0, -chi);
					planes[i]->backprop[x][y] = polar(1.0, chi);
				}
			}
		}

		cout << "initializing approx ... ";

		//setup the initial approximations
		ew.Allocate(   padding, padding, sizeof(Complex));
		ewfft.Allocate(padding, padding, sizeof(Complex));
		
		for(int x = 0; x < size; x++)
			for(int y = 0; y < size; y++)
				ew[x][y] = Complex(1, 0);

		fft2d fftfwd(-1, ew, ewfft);
		fftfwd.fft(ew, ewfft);

		cout << "done in " << (int)((Time() - start)*1000) << " msec\n";

		// Run iterations
		for(int iter = 0; iter < iters; iter++){

			Time startiter;
			cout << "Iter " << iter << " ... ";

			#pragma omp parallel
			for(int p = 0; p < nplanes; p++){
				Plane * plane = planes[p];

				// Propagate EW to each plane
				// EWplanes(p,:) = ifft(EWfft.*prop(p,:).*A);
				for(int x = 0; x < padding; x++){
					for(int y = 0; y < padding; y++){
						if(q2(x,y) <= qmax*qmax){
							plane->ewfft[x][y] *= plane->prop[x][y];
						}else{
							plane->ewfft[x][y] = 0;
						}
					}
				}
				plane->fftbwd.fft(plane->ewfft, plane->ew);
				
				
		
				// Replace EW amplitudes
				//EWplanes(p,1:(N/2)) = Astack(p,:).*exp(1i*angle(EWplanes(p,1:(N/2))));
				for(int x = 0; x < size; x++)
					for(int y = 0; y < size; y++)
						plane->ew[x][y] = polar(plane->amplitude[x][y], arg(plane->ew[x][y]));

				// Back propagate EW to zero plane
				// EWplanesfft(a2,:) = fft(EWplanes(a2,:)).*backprop(a2,:).*A;
				for(int x = 0; x < padding; x++){
					for(int y = 0; y < padding; y++){
						if(q2(x,y) <= qmax*qmax){
							plane->ew[x][y] *= plane->backprop[x][y];
						}else{
							plane->ew[x][y] = 0;
						}
					}
				}
				plane->fftfwd.fft(plane->ew, plane->ewfft);
			}
		
			// Find mean EW and output old phase
			//EWfft = mean(EWplanesfft,1);	
			#pragma omp parallel
			for(int x = 0; x < padding; x++){
				for(int y = 0; y < padding; y++){
					Complex mean = 0;
					for(int p = 0; p < nplanes; p++)
						mean += planes[p]->ewfft[x][y];
					ewfft[x][y] = mean / (double)nplanes;
				}
			}
			cout << " done in " << (int)((Time() - startiter)*1000) << " msec\n";
		}
	}
};

int main(int argc, char **argv){
	if(argc < 2)
		die(1, "Must pass the name of a config file");

	PEWR pewr(argv[1]);
}


