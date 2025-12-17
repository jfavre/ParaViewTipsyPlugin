#ifndef TIPSY_FILE_H
#define TIPSY_FILE_H
//#define USE_VTK_SWAP 1
/*
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 * 		Tim Dykes 
 *

	File: tipsy_file.h
	Purpose: Encapsulates a Tipsy file, with option to swap endian, and have 1 int
			 padding on the header. Both default to yes.

			Structures gas/dark/star_particle and header all sourced from tipsydefs.h 
			from tipsy tools found here http://www-hpcc.astro.washington.edu/tools/tipsy/tipsy.html

 */

#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#ifdef USE_VTK_SWAP
#include "vtkByteSwap.h"
#endif
#define MAXDIM 3

/* 
	Structures for tipsy particle types & header
*/

struct gas_particle {
    float mass;
    float pos[MAXDIM];
    float vel[MAXDIM];
    float rho;
    float temp;
    float hsmooth;
    float metals;
    float phi;
};

struct dark_particle {
    float mass;
    float pos[MAXDIM];
    float vel[MAXDIM];
    float eps;
    float phi;
};

struct star_particle {
    float mass;
    float pos[MAXDIM];
    float vel[MAXDIM];
    float metals;
    float tform;
    float eps;
    float phi;
};

struct header {
    double time;
    int nbodies;
    int ndim;
    int nsph;
    int ndark;
    int nstar;
    int pad;
};

class TipsyFile{
public:
  //  Data
  header h;
  gas_particle*   sph;
  dark_particle* dark;
  star_particle* star;

  // Extra
  const char* name;
  std::ifstream src;
  bool swap_endian;
  bool header_read;

  const size_t static stride_of_gas_particle = sizeof(gas_particle)/sizeof(float);
  const size_t static  stride_of_star_particle = sizeof(star_particle)/sizeof(float);
  const size_t static  stride_of_dark_particle = sizeof(dark_particle)/sizeof(float);

  TipsyFile() { header_read = false; sph = nullptr; dark = nullptr; star = nullptr;}

  TipsyFile(const char* filename, bool swap = false)
    {
    header_read = false; 
    swap_endian = false;
    sph = nullptr; 
    dark = nullptr; 
    star = nullptr;
    FileOpen(filename, swap);
    }

  const float* gas_ptr(){ return (const float*)sph;};
  const float* dark_ptr(){ return (const float*)dark;};
  const float* star_ptr(){ return (const float*)star;};
  
  // Create a new tipsy file
  //void create(){}

  ~TipsyFile()
    {
    FileClose();
    std::cout << "TipsyReader: releasing all resources\n";
    }
  void FileOpen(const char* filename, bool swap = false)
    {
    if(src.is_open())
      std::cerr << "TipsyFile: File " << filename << " is already open!\n";

    src.open(filename, std::ios::binary);
    if(!src.is_open())
      std::cerr << "TipsyFile: Cannot open file: " <<  filename << "\n";

    swap_endian = swap;		
    name = filename;
    }

  void Free_Sph_Buffer()
    {
    if(sph)
      {
      free(sph); sph = nullptr;
      std::cout << "free(sph)...\n";
      }
    }
  void Free_Dark_Buffer()
    {
    if(dark)
      {
      free(dark); dark = nullptr;
      std::cout << "free(dark)...\n";
      }
    }
  void Free_Star_Buffer()
    {
    if(star)
      {
      free(star); star = nullptr;
      std::cout << "free(star)...\n";
      }
    }

  void FileClose()
    {
    Free_Star_Buffer();
    Free_Dark_Buffer();
    Free_Sph_Buffer();
    if(src.is_open())
      src.close();
    }

  bool read_header(bool hasPad = true)
    {
    if(header_read)
      return true;

    int pad = 0;
    if(!hasPad)
      {
      pad = sizeof(int);
      }
// read first and do not swap. Check for length
    src.read((char*)&h, sizeof(header)-pad);

    src.seekg(0, src.end);
    long length = src.tellg();
    src.seekg(sizeof(header), src.beg);
    long size = sizeof(header) + 
                h.ndark * sizeof(dark_particle) +
                h.nsph * sizeof(gas_particle) +
                h.nstar * sizeof(star_particle);

    if(size != length)
      {
      //std::cerr << " detected different endianness. All data will be byte-swapped\n";
#ifdef USE_VTK_SWAP
      vtkByteSwap::Swap8BE(&h.time);
      vtkByteSwap::Swap4BE(&h.nbodies);
      vtkByteSwap::Swap4BE(&h.ndim);
      vtkByteSwap::Swap4BE(&h.nsph);
      vtkByteSwap::Swap4BE(&h.ndark);
      vtkByteSwap::Swap4BE(&h.nstar);
#else
      byteswap(&h.time);
      byteswap(&h.nbodies);
      byteswap(&h.ndim);
      byteswap(&h.nsph);
      byteswap(&h.ndark);
      byteswap(&h.nstar);
#endif
      swap_endian = true;
      }
    header_read = h.nbodies == (h.nsph + h.ndark + h.nstar);
    return header_read;
    }

  bool report_header()
    {
    if(header_read)
      {
/*
      std::cerr << __LINE__ << ": TipsyFile: "<< name << "\ntime=" << h.time 
                << " nbodies=" << h.nbodies << " ndim=" << h.ndim 
                << " nsph="<< h.nsph << " ndark=" << h.ndark 
                << " nstar=" << h.nstar 
                << ". SwapEndian = " << swap_endian << std::endl;
*/
      return h.nbodies == (h.nsph + h.ndark + h.nstar);
      }
    else
      std::cerr << "TipsyFile: report_header(): havent read header yet\n";
    return false;
    }

void read_all(bool hasPad = true)
  {
    if(!src.is_open())
      std::cerr << "TipsyFile: read_all(): file is not open\n";

    if(!header_read)
      read_header(hasPad);
    int base_offset = src.tellg();
    //std::cerr << __LINE__ << ": current_offset = " << src.tellg() << std::endl;
    if(h.nsph > 0){
      sph = (gas_particle*)malloc(h.nsph*sizeof(gas_particle));
      std::cout << " allocating " << h.nsph  << " SPH  particles of size " << sizeof(gas_particle)  << " = " << h.nsph*sizeof(gas_particle)/1024 << " Kbytes\n";
      }

    if(h.ndark > 0)
      dark = (dark_particle*)malloc(h.ndark*sizeof(dark_particle));

    if(h.nstar > 0)
      star = (star_particle*)malloc(h.nstar*sizeof(star_particle));

    if( (h.nsph > 0 && sph == nullptr) || (h.ndark > 0 && dark == nullptr) || (h.nstar > 0 && star == nullptr))
     {
     // Couldnt allocate, cleanup and quit
     FileClose();
     std::cerr << "Could not allocate memory for particles...\n";
     }

// by adding explicit seekg(), we make it possible to skip reading particles
// of a particular type
  // Read sph
  if(h.nsph)
    {
    src.seekg(base_offset, src.beg);
    src.read((char*)sph, h.nsph * sizeof(gas_particle));
    if(swap_endian)
      {
#ifdef USE_VTK_SWAP
       vtkByteSwap::Swap4BERange(sph, h.nsph * stride_of_gas_particle);
#else
       gas_particle* pp = sph;
       for(auto i = 0; i < h.nsph; i++, pp++)
         {
         for(size_t j = 0; j < stride_of_gas_particle; j++)
           byteswap(&((float*)pp)[j]);
         }
#endif
      }
    }

// Read dark
  if(h.ndark)
    {
    src.seekg(h.nsph * sizeof(gas_particle) + base_offset, src.beg);
    src.read((char*)dark, h.ndark * sizeof(dark_particle));
    if(swap_endian)
      {
#ifdef USE_VTK_SWAP
      vtkByteSwap::Swap4BERange(dark, h.ndark * stride_of_dark_particle);
#else
      dark_particle* pp = dark;
      for(auto i = 0; i < h.ndark; i++, pp++)
        {
        for(size_t j = 0; j < stride_of_dark_particle; j++)
          byteswap(&((float*)pp)[j]);
        }
#endif
      }
    }

// Read star
  if(h.nstar)
    {
    src.seekg(h.nsph * sizeof(gas_particle) + h.ndark * sizeof(dark_particle) + base_offset, src.beg);
    src.read((char*)star, h.nstar * sizeof(star_particle));
    if(swap_endian)
      {
#ifdef USE_VTK_SWAP
      vtkByteSwap::Swap4BERange(star, h.nstar * stride_of_star_particle);
#else
      star_particle* pp = star;
      for(auto i = 0; i < h.nstar; i++, pp++)
        {
        for(size_t j = 0; j < stride_of_star_particle; j++)
          byteswap(&((float*)pp)[j]);
        }
#endif
      }
    }

  std::cout << __LINE__ << ": TipsyFile: read file " << name
            << "\nnbodies: " << h.nbodies
            << "\nnsph: " << h.nsph
            << "\nndark: " << h.ndark
            << "\nnstar: " << h.nstar
            << "\nswapped endian: " << swap_endian << std::endl;
  }

int split_particlesSet(int N, int piece, int numPieces, int& standard_load)
  {
// I call standard_load, the size of all the first pieces except the last one
// I call load, the size of the last piece
  int load;
  if(numPieces == 1)
    standard_load = load = N;
  else
    {
    standard_load = load = N / numPieces;
    if (piece == (numPieces-1))
      {
      load = N - (numPieces-1) * standard_load;
      }
    }
  return load;
  }

void read_gas_piece(int piece, int numPieces, int &n1, [[maybe_unused]] bool hasPad = true)
  {
// if n1 or n2 or n3 = 0, it means we do not read the particular type.
  int standard_load[3]={0,0,0};

  if(!src.is_open())
    std::cerr << "TipsyFile: read_all(): file is not open\n";

  long base_offset = sizeof(header);

  if(h.nsph > 0 and n1)
    {
    n1 = split_particlesSet(h.nsph, piece, numPieces, standard_load[0]);
    sph = (gas_particle*)malloc(n1 * sizeof(gas_particle));
    std::cout << "CPU " << piece << " allocating " << n1  << " SPH  particles of size " << sizeof(gas_particle)  << " = " << n1*sizeof(gas_particle)/1024 << " Kbytes\n";
    }
  else n1=0;


  if( (n1 > 0 && sph == nullptr))
    {
    // Couldnt allocate, cleanup and quit
    FileClose();
    std::cerr << __LINE__ << "Could not allocate memory for particles...\n";
    }

// Read sph
  if(n1)
    {
    long offset = piece * standard_load[0] * sizeof(gas_particle) + base_offset;
    src.seekg(offset, src.beg);

    src.read((char*)sph, n1 * sizeof(gas_particle));

    if(swap_endian)
      {
#ifdef USE_VTK_SWAP
       vtkByteSwap::Swap4BERange(sph, n1 * sizeof(gas_particle)/sizeof(float));
#else
	  gas_particle* pp = sph;
	  for(auto i = 0; i < n1; i++, pp++)
        {
	    for(size_t j = 0; j < sizeof(gas_particle)/sizeof(float); j++)
		  byteswap(&((float*)pp)[j]);
	    }
#endif
	  }
/*
    float min= 1e30;
    float max=-min;
    float *fp = (float*)sph + 7;
    for(unsigned i = 0; i < n1; i++)
      {
      if(*fp > max) max = *fp; if(*fp < min) min = *fp;
      fp += sizeof(gas_particle)/sizeof(float);
      }
    errs << "Density SPH: min=" << min << ", max="<< max << endl;
*/
/*
    fp = (float*)sph; min= 1e30; max=-min;
    for(unsigned i = 0; i < n1; i++)
      {

      if(*fp > max) max = *fp; if(*fp < min) min = *fp;
      fp += sizeof(gas_particle)/sizeof(float);
      }
    errs << "Mass SPH: min=" << min << ", max="<< max << endl;
*/
    }

  std::cout << __LINE__ << ": TipsyFile: read file " << name
            << "\ntime   : " << h.time
            << "\nnbodies: " << h.nbodies
            << "\nnsph   : " << h.nsph
            << "\nndark  : " << h.ndark
            << "\nnstar  : " << h.nstar
            << "\nswapped endian: " << swap_endian << std::endl;
} // read_gas_piece()

void read_dark_matter_piece(int piece, int numPieces, int &n2, [[maybe_unused]] bool hasPad = true)
  {
// if n1 or n2 or n3 = 0, it means we do not read the particular type.
  int standard_load[3]={0,0,0};

  if(!src.is_open())
    std::cerr << "TipsyFile: read_all(): file is not open\n";

  long base_offset = sizeof(header);

  if(h.ndark > 0  and n2)
    {
    n2 = split_particlesSet(h.ndark, piece, numPieces, standard_load[1]);
    dark = (dark_particle*)malloc(n2*sizeof(dark_particle));
    //std::cerr << "CPU " << piece << " allocating " << n2  << " DARK particles of size " << sizeof(dark_particle) << "= " << n2*sizeof(dark_particle)<< "\n";
    }
  else n2=0;

  if( (n2 > 0 && dark == nullptr))
    {
    // Couldnt allocate, cleanup and quit
    FileClose();
    std::cerr << __LINE__ << "Could not allocate memory for particles...\n";
    }

//Read dark
  if(n2)
    {
    long offset = piece * standard_load[1] * sizeof(dark_particle) +
                  h.nsph * sizeof(gas_particle) +
                  base_offset;
    src.seekg(offset, src.beg);

    src.read((char*)dark, n2 * sizeof(dark_particle));

    if(n2 * sizeof(dark_particle) != src.gcount())
      std::cerr << "ERROR reading DARK: only " << src.gcount() << " read. instead of "<< n2 * sizeof(dark_particle) << std::endl;

	if(swap_endian)
	  {
#ifdef USE_VTK_SWAP
     vtkByteSwap::Swap4BERange(dark, n2 * sizeof(dark_particle)/sizeof(float));
#else
	  dark_particle* pp = dark;
	  for(auto i = 0; i < n2; i++, pp++)
	    {
	    for(size_t j = 0; j < sizeof(dark_particle)/sizeof(float); j++)
	   byteswap(&((float*)pp)[j]);
	    }
#endif
	  }
    }
} // read_dark_matter_piece()

void read_star_piece(int piece, int numPieces, int &n3, [[maybe_unused]] bool hasPad = true)
  {
// if n1 or n2 or n3 = 0, it means we do not read the particular type.
  int standard_load[3]={0,0,0};

  if(!src.is_open())
    std::cerr << "TipsyFile: read_all(): file is not open\n";

  long base_offset = sizeof(header);

  if(h.nstar > 0  and n3)
    {
    n3 = split_particlesSet(h.nstar, piece, numPieces, standard_load[2]);
    star = (star_particle*)malloc(n3*sizeof(star_particle));
    //std::cerr << "CPU " << piece << " allocating " << n3  << " STAR particles of size " << sizeof(star_particle) << "= " << n3*sizeof(star_particle)<< "\n";
    }
  else n3=0;

  if((n3 > 0 && star == nullptr))
    {
    // Couldnt allocate, cleanup and quit
    FileClose();
    std::cerr << __LINE__ << "Could not allocate memory for particles...\n";
    }

// Read star
  if(n3)
    {
    long offset = piece * standard_load[2] * sizeof(star_particle) +
                  h.nsph * sizeof(gas_particle) +
                  h.ndark * sizeof(dark_particle) +
                  base_offset;
    src.seekg(offset, src.beg);

    src.read((char*)star, n3 * sizeof(star_particle));

    if(n3 * sizeof(star_particle) != src.gcount())
      std::cerr << "ERROR reading STAR: only " << src.gcount() << " read. instead of "<< n3 * sizeof(star_particle) << std::endl;

	if(swap_endian)
	  {
#ifdef USE_VTK_SWAP
     vtkByteSwap::Swap4BERange(star, n3 * sizeof(star_particle)/sizeof(float));
#else
	  star_particle* pp = star;
	  for(auto i = 0; i < n3; i++, pp++)
	    {
	    for(size_t j = 0; j < sizeof(star_particle)/sizeof(float); j++)
	   byteswap(&((float*)pp)[j]);
	    }
#endif
	  }
    }
} // read_star_piece()

void write(std::string fileout, bool hasPad = true)
  {
// Output file
  std::ofstream out(fileout.c_str(), std::ios::binary);

  if(!out.is_open())
    {
    std::cerr << "TipsyFile: write() could not open file " << fileout.c_str() << " for output.\n";
    }

// Header struct includes padding, if we dont want padding dont write the final int.
// We always swap endianness of header - but if we didnt also swap the endianness 
// of the file, then we should swap back header

  int nsph = h.nsph;
  int ndark = h.ndark;
  int nstar = h.nstar;
/*
  if(!swap_endian)
    {
#ifdef USE_VTK_SWAP
          vtkByteSwap::Swap8BE(&h.time);
          vtkByteSwap::Swap4BE(&h.nbodies);
          vtkByteSwap::Swap4BE(&h.ndim);
          vtkByteSwap::Swap4BE(&h.nsph);
          vtkByteSwap::Swap4BE(&h.ndark);
          vtkByteSwap::Swap4BE(&h.nstar);
          vtkByteSwap::Swap4BE(&h.nstar);
#else
			byteswap(&h.time);
	        byteswap(&h.nbodies);
	        byteswap(&h.ndim);
	        byteswap(&h.nsph);
	        byteswap(&h.ndark);
	        byteswap(&h.nstar);
#endif
    }
*/
  out.write((char*)&h, sizeof(header) - (!hasPad ? sizeof(int) : 0) );

  // Write gas
  if(nsph > 0)
    out.write((char*)sph, sizeof(gas_particle) * nsph);

  // Write dark
  if(ndark > 0)
    out.write((char*)dark, sizeof(dark_particle) * ndark);

  // Write star
  if(nstar > 0)
    out.write((char*)star, sizeof(star_particle) * nstar);

  out.close();
  };
  
private:
#ifndef USE_VTK_SWAP
template <typename T>
void byteswap(T* in)
{
    unsigned size = sizeof(T);
    for(unsigned i = 0; i < size/2; i++)
    {
        // Swap bytes
        ((char*)in)[i]      ^= ((char*)in)[size-i-1];
        ((char*)in)[size-i-1] ^= ((char*)in)[i];
        ((char*)in)[i]      ^= ((char*)in)[size-i-1];
    }
}
#endif
};

#endif
