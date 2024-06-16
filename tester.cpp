#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <array>
#include <memory>
#include <algorithm>
#include <numeric>

// Structure for results
struct Result
{
    std::string benchmark;
    int threads;
    double minTime;
    double avgTime;
    double maxTime;
};

// Function to execute system command and capture output
std::string exec(const char *cmd)
{
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe)
    {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
    {
        result += buffer.data();
    }
    return result;
}

// Function to log to file
void logToFile(const std::string &filename, const std::string &message)
{
    std::ofstream logFile;
    logFile.open(filename, std::ios_base::app); // Append mode
    if (logFile.is_open())
    {
        logFile << message << std::endl;
        logFile.close();
    }
    else
    {
        std::cerr << "Unable to open log file: " << filename << std::endl;
    }
}

int main()
{
    std::vector<std::string> benchmarks = {
        "alpine", "ackley", "sphere", "schwefel", "rastrigin", "griewank", "csendes", "colville",
        "easom", "michalewicz", "shekel", "schwefel2_4", "schwefel2_6", "schaffer",
        "sumSquares", "step2", "quartic", "powell", "rosenbrock", "dixonPrice", "schwefel1_2", "schwefel2_21", "schwefel2_22"};
    std::vector<int> threads = {1, 2, 4, 8, 16, 32};

    std::string summaryLogFilename = "cuckoo_search_summary_log.txt";
    std::string detailedLogFilename = "cuckoo_search_detailed_log.txt";
    std::string csvFilename = "cuckoo_search_results.csv";

    // Remove old log files
    std::remove(summaryLogFilename.c_str());
    std::remove(detailedLogFilename.c_str());

    std::ofstream csvFile(csvFilename);
    csvFile << "Benchmark Function,Number of Threads,Min Time,Avg Time,Max Time\n";

    for (const auto &benchmark : benchmarks)
    {
        for (const auto &numThreads : threads)
        {
            std::vector<double> execTimes;
            for (int i = 0; i < 5; ++i) // Repeat 5 times
            {
                std::string command = "cmd /c \"set OMP_NUM_THREADS=" + std::to_string(numThreads) + " && cuckoo_search " + std::to_string(numThreads) + " " + benchmark + " " + detailedLogFilename + "\"";
                std::cout << "Running: " << command << std::endl;
                auto output = exec(command.c_str());

                // Parse execution time
                std::istringstream output_stream(output);
                std::string line;
                while (std::getline(output_stream, line))
                {
                    if (line.find("Vreme izvršavanja") != std::string::npos)
                    {
                        double execTime = std::stod(line.substr(line.find(":") + 1));
                        execTimes.push_back(execTime);
                    }
                }

                // Log detailed output
                logToFile(detailedLogFilename, output);
            }

            // Calculate min, avg, and max times
            if (!execTimes.empty())
            {
                double minTime = *std::min_element(execTimes.begin(), execTimes.end());
                double avgTime = std::accumulate(execTimes.begin(), execTimes.end(), 0.0) / execTimes.size();
                double maxTime = *std::max_element(execTimes.begin(), execTimes.end());

                // Write to CSV
                csvFile << benchmark << "," << numThreads << "," << minTime << "," << avgTime << "," << maxTime << "\n";

                // Log summary
                logToFile(summaryLogFilename, "Benchmark: " + benchmark + ", Threads: " + std::to_string(numThreads) + ", Min Time: " + std::to_string(minTime) + ", Avg Time: " + std::to_string(avgTime) + ", Max Time: " + std::to_string(maxTime));
            }
        }
    }

    csvFile.close();
    std::cout << "Testing completed." << std::endl;

    return 0;
}