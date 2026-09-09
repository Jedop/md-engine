#include "utils.hpp"

void print_progress(int step, int total, double elapsed) {
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

    std::cout << "] " << int(progress * 100) << "% (" << step << "/" << total << ") ";

    if (step == total) {
        // Print final elapsed time when finished
        std::cout << "Elapsed Time: " << std::fixed << std::setprecision(2) << elapsed << "s    \n";
    } else {
        // Calculate and print ETA during the run
        double eta = elapsed * (1.0 / progress - 1.0);
        std::cout << "ETA: " << int(eta) << "s    ";
    }
    std::cout << std::flush;
}
