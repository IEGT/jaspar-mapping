#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <filesystem> // test and creating for directories
#include <algorithm> // For removing brackets, replace
#include <iomanip>   // For formatted output / setting precision
#include <chrono>    // For progress timing
#include <getopt.h>  // For GNU Getopt
#include <math.h>

/** \brief Function to display the progress bar
 */
void displayProgressBar(const float& progress) {
    int barWidth = 75;  // Width of the progress bar
    std::cerr << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cerr << "=";
        else if (i == pos) std::cerr << ">";
        else std::cerr << " ";
    }
    std::cerr << "] " << std::fixed << std::setprecision(3) << progress * 100.0 << " %\r";
    std::cerr.flush();
}
