#include <iostream>
#include "timer.hpp"

int main()
{
    // create a timer object
    timer stopwatch;
    // start timer
    stopwatch.start();
    // Compute sum of first ten million positive integers.
    unsigned long sum = 0;
    for (unsigned long k = 1; k < 10000001; k++)
    {
        sum += k;
    }
    // end timer
    stopwatch.stop();

    std::cout << "Elapsed time for computing the sum of first\n"
              << "ten million positive integers is " << std::scientific
              << stopwatch.get_elapsed_time() << " seconds.\n"
              << std::endl;
    return 0;
}