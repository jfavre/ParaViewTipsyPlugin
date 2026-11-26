
#include <cstdlib>
#include <iostream>
#include "tipsy_file.h"

// reads a Tipsy file and writes it back out (a no-op round-trip when both
// the input and the output use the same header padding)
int main(int argc, char * argv[])
{
  if (argc < 3)
    {
    std::cerr << "usage: TipsyReadAndSwap fileIn fileOut\n";
    return EXIT_FAILURE;
    }
  TipsyFile filein;
  bool res = filein.FileOpen(argv[1]);
  if (!res)
    {
    std::cerr << "Cannot open " << argv[1] << "\n";
    }
  res = res && filein.read_all();
  res = filein.write(argv[2], true) && res;
  return res ? EXIT_SUCCESS : EXIT_FAILURE;
}
