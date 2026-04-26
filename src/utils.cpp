#include "utils.hpp"

void print_progress(int step, int total, double eta) {
  const int bar_width = 40;

  double progress = (double)step / total;
  int pos = bar_width * progress;

  std::cout << "\r[";
  for (int i = 0; i < bar_width; ++i) {
    if (i < pos)
      std::cout << "#";
    else if (i == pos)
      std::cout << ">";
    else
      std::cout << "-";
  }

  std::cout << "] " << int(progress * 100) << "% "
            << "(" << step << "/" << total << ") "
            << "ETA: " << int(eta) << "s" << std::flush;
}
