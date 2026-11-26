#ifndef TIPSY_FILE_H
#define TIPSY_FILE_H

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

	Modernized:
	  - particle buffers are std::vector (allocation failure is reported, not fatal
	    fall-through into src.read() with a null pointer)
	  - 64-bit clean: file offsets and sizes use std::streamoff, so files larger
	    than 2GB are handled correctly
	  - the header size check honors the hasPad option (previously an un-padded
	    little-endian file was falsely detected as big-endian and byte-swapped)
	  - read_*_piece() honor the hasPad option
	  - endianness detection is followed by a structural validation of the header
	  - vtkByteSwap is always used (vectorized), the hand-rolled byteswap is gone
	  - all debug chatter to std::cout removed; only errors go to std::cerr

 */

#include "vtkByteSwap.h"

#include <cstddef>
#include <fstream>
#include <iostream>
#include <new>
#include <string>
#include <vector>

/*
	Structures for tipsy particle types & header
*/

struct gas_particle {
    float mass;
    float pos[3];
    float vel[3];
    float rho;
    float temp;
    float hsmooth;
    float metals;
    float phi;
};

struct dark_particle {
    float mass;
    float pos[3];
    float vel[3];
    float eps;
    float phi;
};

struct star_particle {
    float mass;
    float pos[3];
    float vel[3];
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

// the whole class depends on these structs having no compiler-added padding
static_assert(sizeof(gas_particle) == 12 * sizeof(float),
              "gas_particle must be 12 floats with no padding");
static_assert(sizeof(dark_particle) == 9 * sizeof(float),
              "dark_particle must be 9 floats with no padding");
static_assert(sizeof(star_particle) == 11 * sizeof(float),
              "star_particle must be 11 floats with no padding");
static_assert(sizeof(header) == sizeof(double) + 6 * sizeof(int),
              "header must be a double plus 6 ints with no padding");

class TipsyFile{
public:
  //  Data
  header h{};
  std::vector<gas_particle>  sph;
  std::vector<dark_particle> dark;
  std::vector<star_particle> star;

  // Extra
  std::string name;
  bool swap_endian = false;
  bool header_read = false;

  static constexpr std::size_t stride_of_gas_particle  = sizeof(gas_particle)  / sizeof(float);
  static constexpr std::size_t stride_of_dark_particle = sizeof(dark_particle) / sizeof(float);
  static constexpr std::size_t stride_of_star_particle = sizeof(star_particle) / sizeof(float);

  TipsyFile() = default;

  TipsyFile(const char* filename, bool swap = false)
    {
    this->FileOpen(filename, swap);
    }

  // the ifstream member makes the class non-copyable anyway, be explicit
  TipsyFile(const TipsyFile&)            = delete;
  TipsyFile& operator=(const TipsyFile&) = delete;

  ~TipsyFile()
    {
    this->FileClose();
    }

  const float* gas_ptr()  const { return reinterpret_cast<const float*>(this->sph.data());  }
  const float* dark_ptr() const { return reinterpret_cast<const float*>(this->dark.data()); }
  const float* star_ptr() const { return reinterpret_cast<const float*>(this->star.data()); }

  // returns false if the file could not be opened ; the caller reports
  // (CanReadFile probes arbitrary files and must stay silent)
  bool FileOpen(const char* filename, bool swap = false)
    {
    if (src.is_open())
      {
      src.close();
      }
    src.clear();
    header_read = false;
    swap_endian = swap;
    name        = filename ? filename : "";
    src.open(name, std::ios::binary);
    return src.is_open();
    }

  // release a particle buffer (and its capacity) once its content has been
  // transferred to VTK arrays
  void Free_Sph_Buffer()  { sph   = {}; }
  void Free_Dark_Buffer() { dark  = {}; }
  void Free_Star_Buffer() { star  = {}; }

  void FileClose()
    {
    Free_Star_Buffer();
    Free_Dark_Buffer();
    Free_Sph_Buffer();
    if (src.is_open())
      {
      src.close();
      }
    }

  // Reads and validates the header. On a size mismatch (the file length does
  // not correspond to what the header announces) the header is byte-swapped
  // and swap_endian is set for the particle sections.
  bool read_header(bool hasPad = true)
    {
    if (header_read)
      {
      return true;
      }
    if (!src.is_open())
      {
      return false;
      }

    // a header without padding lacks the trailing int
    const std::streamoff pad        = hasPad ? 0 : static_cast<std::streamoff>(sizeof(int));
    const std::streamoff dataOffset = static_cast<std::streamoff>(sizeof(header)) - pad;

    h = header{};
    src.seekg(0, src.beg);
    src.read(reinterpret_cast<char*>(&h), dataOffset);
    if (!src)
      {
      return false;
      }

    src.seekg(0, src.end);
    const std::streamoff length = src.tellg();

    // detected wrong endianness : all floats/integers will be byte-swapped
    if (this->announced_size(hasPad) != length)
      {
      swap_header();
      swap_endian = true;
      }

    header_read = this->validate_header();
    if (header_read)
      {
      src.seekg(dataOffset, src.beg);
      }
    return header_read;
    }

  bool report_header() const
    {
    if (header_read)
      {
      return this->validate_header();
      }
    std::cerr << "TipsyFile: report_header(): havent read header yet\n";
    return false;
    }

  // read all particle types of the (whole) file
  bool read_all(bool hasPad = true)
    {
    if (!header_read && !read_header(hasPad))
      {
      return false;
      }
    int n1 = h.nsph, n2 = h.ndark, n3 = h.nstar;
    bool res = read_gas_piece(0, 1, n1, hasPad);
    res = read_dark_matter_piece(0, 1, n2, hasPad) && res;
    res = read_star_piece(0, 1, n3, hasPad) && res;
    return res;
    }

  // compute the number of particles ('load') of 'piece' out of numPieces,
  // 'standard_load' receives the size of all the pieces except the last one,
  // the last piece gets the remainder
  int split_particlesSet(int N, int piece, int numPieces, int& standard_load) const
    {
    int load;
    if (numPieces == 1)
      {
      standard_load = load = N;
      }
    else
      {
      standard_load = load = N / numPieces;
      if (piece == (numPieces - 1))
        {
        load = N - (numPieces - 1) * standard_load;
        }
      }
    return load;
    }

  // read a contiguous piece of the particle type, n is in/out: caller sets it
  // to a nonzero value to request a read, on return it holds the number of
  // particles read (0 on error or empty type)
  bool read_gas_piece(int piece, int numPieces, int& n, bool hasPad = true)
    {
    return read_section(sph, h.nsph, this->data_offset(hasPad), piece, numPieces, n);
    }

  bool read_dark_matter_piece(int piece, int numPieces, int& n, bool hasPad = true)
    {
    return read_section(dark, h.ndark,
                        this->data_offset(hasPad)
                          + to_size(h.nsph) * static_cast<std::streamoff>(sizeof(gas_particle)),
                        piece, numPieces, n);
    }

  bool read_star_piece(int piece, int numPieces, int& n, bool hasPad = true)
    {
    return read_section(star, h.nstar,
                        this->data_offset(hasPad)
                          + to_size(h.nsph)  * static_cast<std::streamoff>(sizeof(gas_particle))
                          + to_size(h.ndark) * static_cast<std::streamoff>(sizeof(dark_particle)),
                        piece, numPieces, n);
    }

  bool write(std::string fileout, bool hasPad = true)
    {
    std::ofstream out(fileout, std::ios::binary);
    if (!out.is_open())
      {
      std::cerr << "TipsyFile: write() could not open file " << fileout << " for output.\n";
      return false;
      }

    // Header struct includes padding, if we dont want padding dont write the final int.
    // The header is written in its current (host) endianness, as are the particles.
    out.write(reinterpret_cast<const char*>(&h),
              sizeof(header) - (hasPad ? 0 : sizeof(int)));

    if (!sph.empty())
      {
      out.write(reinterpret_cast<const char*>(sph.data()),
                static_cast<std::streamsize>(sph.size() * sizeof(gas_particle)));
      }
    if (!dark.empty())
      {
      out.write(reinterpret_cast<const char*>(dark.data()),
                static_cast<std::streamsize>(dark.size() * sizeof(dark_particle)));
      }
    if (!star.empty())
      {
      out.write(reinterpret_cast<const char*>(star.data()),
                static_cast<std::streamsize>(star.size() * sizeof(star_particle)));
      }

    const bool res = out.good();
    out.close();
    if (!res)
      {
      std::cerr << "TipsyFile: write() error while writing " << fileout << "\n";
      }
    return res;
    }

protected:
  std::ifstream src;

  // offset of the first particle section, honoring the header padding option
  static std::streamoff data_offset(bool hasPad)
    {
    return static_cast<std::streamoff>(sizeof(header))
         - (hasPad ? 0 : static_cast<std::streamoff>(sizeof(int)));
    }

  static std::streamoff to_size(int n)
    {
    return static_cast<std::streamoff>(static_cast<long long>(n));
    }

  // total file size announced by the (current) header
  std::streamoff announced_size(bool hasPad) const
    {
    return this->data_offset(hasPad)
         + to_size(h.ndark) * static_cast<std::streamoff>(sizeof(dark_particle))
         + to_size(h.nsph)  * static_cast<std::streamoff>(sizeof(gas_particle))
         + to_size(h.nstar) * static_cast<std::streamoff>(sizeof(star_particle));
    }

  // structural sanity of the header
  bool validate_header() const
    {
    if (h.ndim < 1 || h.ndim > 3)
      {
      return false;
      }
    if (h.nbodies < 0 || h.nsph < 0 || h.ndark < 0 || h.nstar < 0)
      {
      return false;
      }
    return h.nbodies == (h.nsph + h.ndark + h.nstar);
    }

  void swap_header()
    {
    vtkByteSwap::Swap8BE(&h.time);
    vtkByteSwap::Swap4BE(&h.nbodies);
    vtkByteSwap::Swap4BE(&h.ndim);
    vtkByteSwap::Swap4BE(&h.nsph);
    vtkByteSwap::Swap4BE(&h.ndark);
    vtkByteSwap::Swap4BE(&h.nstar);
    }

  // common implementation of the read_*_piece() methods : compute this piece's
  // particle count, resize the buffer, seek to the 64-bit clean offset, read,
  // verify the read length and byte-swap in place if needed
  template <typename T>
  bool read_section(std::vector<T>& buffer, int count, std::streamoff sectionStart,
                    int piece, int numPieces, int& n)
    {
    if (!src.is_open())
      {
      std::cerr << "TipsyFile: " << name << ": file is not open\n";
      return false;
      }
    if (count <= 0 || n <= 0)
      {
      n = 0;
      return true;
      }
    if (piece < 0 || numPieces <= 0 || piece >= numPieces)
      {
      std::cerr << "TipsyFile: " << name << ": bad piece request "
                << piece << "/" << numPieces << "\n";
      n = 0;
      return false;
      }

    int standard_load = 0;
    n = this->split_particlesSet(count, piece, numPieces, standard_load);

    buffer.clear();
    try
      {
      buffer.resize(static_cast<std::size_t>(n));
      }
    catch (const std::bad_alloc&)
      {
      std::cerr << "TipsyFile: " << name << ": could not allocate memory for "
                << n << " particles\n";
      n = 0;
      return false;
      }

    const std::streamoff offset =
          sectionStart
        + static_cast<std::streamoff>(static_cast<long long>(piece) * standard_load)
              * static_cast<std::streamoff>(sizeof(T));
    src.seekg(offset, src.beg);
    src.read(reinterpret_cast<char*>(buffer.data()),
             static_cast<std::streamsize>(n) * static_cast<std::streamsize>(sizeof(T)));

    const std::streamsize wanted =
          static_cast<std::streamsize>(n) * static_cast<std::streamsize>(sizeof(T));
    if (src.gcount() != wanted)
      {
      std::cerr << "TipsyFile: " << name << ": ERROR reading particles at offset " << offset
                << ": only " << src.gcount() << " of " << wanted << " bytes read\n";
      buffer = {};
      n = 0;
      return false;
      }

    if (swap_endian)
      {
      vtkByteSwap::Swap4BERange(buffer.data(),
                                static_cast<vtkIdType>(static_cast<std::size_t>(n) * (sizeof(T) / sizeof(float))));
      }
    return true;
    }
};

#endif
