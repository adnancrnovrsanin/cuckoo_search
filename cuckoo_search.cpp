#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <omp.h>

// Funkcija za logovanje u fajl
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

// Levy let funkcija
double levyFlight(double lambda)
{
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double u = distribution(generator);
    double v = distribution(generator);
    double step = u / pow(fabs(v), 1.0 / lambda);
    return step;
}

// Benchmark funkcije
double alpineFunction(const std::vector<double> &x)
{
    double result = 0.0;
    for (double xi : x)
    {
        result += fabs(xi * sin(xi) + 0.1 * xi);
    }
    return result;
}

double ackleyFunction(const std::vector<double> &x)
{
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (double xi : x)
    {
        sum1 += xi * xi;
        sum2 += cos(2.0 * M_PI * xi);
    }
    double n = static_cast<double>(x.size());
    return -20.0 * exp(-0.2 * sqrt(sum1 / n)) - exp(sum2 / n) + 20.0 + M_E;
}

double sphereFunction(const std::vector<double> &x)
{
    double result = 0.0;
    for (double xi : x)
    {
        result += xi * xi;
    }
    return result;
}

double schwefelFunction(const std::vector<double> &x)
{
    double sum1 = 0.0;
    double sum2 = 1.0;
    for (double xi : x)
    {
        sum1 += fabs(xi);
        sum2 *= fabs(xi);
    }
    return sum1 + sum2;
}

double rastriginFunction(const std::vector<double> &x)
{
    double result = 10 * x.size();
    for (double xi : x)
    {
        result += xi * xi - 10 * cos(2 * M_PI * xi);
    }
    return result;
}

double griewankFunction(const std::vector<double> &x)
{
    double sum = 0.0;
    double prod = 1.0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        sum += x[i] * x[i] / 4000.0;
        prod *= cos(x[i] / sqrt(static_cast<double>(i + 1)));
    }
    return sum - prod + 1.0;
}

double csendesFunction(const std::vector<double> &x)
{
    double result = 0.0;
    for (double xi : x)
    {
        result += pow(xi, 6) * (2 + sin(1 / xi));
    }
    return result;
}

double colvilleFunction(const std::vector<double> &x)
{
    if (x.size() < 4)
        return 1e10; // Invalid for dimensions < 4
    return 100 * pow(x[1] - x[0] * x[0], 2) + pow(1 - x[0], 2) + 90 * pow(x[3] - x[2] * x[2], 2) + pow(1 - x[2], 2) + 10.1 * (pow(x[1] - 1, 2) + pow(x[3] - 1, 2)) + 19.8 * (x[1] - 1) * (x[3] - 1);
}

double easomFunction(const std::vector<double> &x)
{
    if (x.size() < 2)
        return 1e10; // Invalid for dimensions < 2
    return -cos(x[0]) * cos(x[1]) * exp(-pow(x[0] - M_PI, 2) - pow(x[1] - M_PI, 2));
}

double michalewiczFunction(const std::vector<double> &x)
{
    double result = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        result += sin(x[i]) * pow(sin((i + 1) * x[i] * x[i] / M_PI), 20);
    }
    return -result;
}

double shekelFunction(const std::vector<double> &x)
{
    const size_t m = 10;
    const double a[10][4] = {
        {4.0, 4.0, 4.0, 4.0}, {1.0, 1.0, 1.0, 1.0}, {8.0, 8.0, 8.0, 8.0}, {6.0, 6.0, 6.0, 6.0}, {3.0, 7.0, 3.0, 7.0}, {2.0, 9.0, 2.0, 9.0}, {5.0, 5.0, 3.0, 3.0}, {8.0, 1.0, 8.0, 1.0}, {6.0, 2.0, 6.0, 2.0}, {7.0, 3.6, 7.0, 3.6}};
    const double c[10] = {0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5};
    double result = 0.0;
    for (size_t j = 0; j < m; ++j)
    {
        double sum = 0.0;
        for (size_t i = 0; i < x.size(); ++i)
        {
            sum += pow(x[i] - a[j][i], 2);
        }
        result -= 1.0 / (sum + c[j]);
    }
    return result;
}

double schwefel2_4Function(const std::vector<double> &x)
{
    if (x.size() < 2)
        return 1e10; // Invalid for dimensions < 2
    return pow(x[0] - 1, 2) + pow(x[1] - x[0] * x[0], 2);
}

double schwefel2_6Function(const std::vector<double> &x)
{
    double result = 0.0;
    for (double xi : x)
    {
        result += xi * sin(sqrt(fabs(xi)));
    }
    return 418.9829 * x.size() - result;
}

double schafferFunction(const std::vector<double> &x)
{
    if (x.size() < 2)
        return 1e10; // Invalid for dimensions < 2
    double numerator = pow(sin(pow(x[0] * x[0] + x[1] * x[1], 0.5)), 2) - 0.5;
    double denominator = pow(1 + 0.001 * (x[0] * x[0] + x[1] * x[1]), 2);
    return 0.5 + numerator / denominator;
}

double sumSquaresFunction(const std::vector<double> &x)
{
    double result = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        result += (i + 1) * x[i] * x[i];
    }
    return result;
}

double step2Function(const std::vector<double> &x)
{
    double result = 0.0;
    for (double xi : x)
    {
        result += pow(floor(xi + 0.5), 2);
    }
    return result;
}

double quarticFunction(const std::vector<double> &x)
{
    double result = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        result += (i + 1) * pow(x[i], 4);
    }
    return result + static_cast<double>(rand()) / RAND_MAX;
}

double powellFunction(const std::vector<double> &x)
{
    if (x.size() < 4)
        return 1e10; // Invalid for dimensions < 4
    double result = 0.0;
    for (size_t i = 0; i < x.size() - 1; i += 4)
    {
        result += pow(x[i] + 10 * x[i + 1], 2) + 5 * pow(x[i + 2] - x[i + 3], 2) + pow(x[i + 1] - 2 * x[i + 2], 4) + 10 * pow(x[i] - x[i + 3], 4);
    }
    return result;
}

double rosenbrockFunction(const std::vector<double> &x)
{
    double result = 0.0;
    for (size_t i = 0; i < x.size() - 1; ++i)
    {
        result += 100 * pow(x[i + 1] - x[i] * x[i], 2) + pow(1 - x[i], 2);
    }
    return result;
}

double dixonPriceFunction(const std::vector<double> &x)
{
    double result = pow(x[0] - 1, 2);
    for (size_t i = 1; i < x.size(); ++i)
    {
        result += i * pow(2 * x[i] * x[i] - x[i - 1], 2);
    }
    return result;
}

double schwefel1_2Function(const std::vector<double> &x)
{
    double result = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
    {
        double sum = 0.0;
        for (size_t j = 0; j <= i; ++j)
        {
            sum += x[j];
        }
        result += sum * sum;
    }
    return result;
}

double schwefel2_21Function(const std::vector<double> &x)
{
    double result = 0.0;
    for (double xi : x)
    {
        result = std::max(result, fabs(xi));
    }
    return result;
}

double schwefel2_22Function(const std::vector<double> &x)
{
    double sum1 = 0.0;
    double sum2 = 1.0;
    for (double xi : x)
    {
        sum1 += fabs(xi);
        sum2 *= fabs(xi);
    }
    return sum1 + sum2;
}

// Paralelizovani Cuckoo Search algoritam
void parallelCuckooSearch(int n, int maxGenerations, double pa, double lambda, std::vector<double> &bestSolution, double &bestFitness, double (*objectiveFunction)(const std::vector<double> &), int dimension, const std::string &logFilename)
{
    std::vector<std::vector<double>> nests(n, std::vector<double>(dimension));
    std::vector<double> fitness(n);

    logToFile(logFilename, "Starting population initialization");

    // Inicijalizacija populacije
    std::default_random_engine generator;
    generator.seed(42); // Postavljanje semena za reproduktibilnost
    std::uniform_real_distribution<double> distribution(-10.0, 10.0);

    double initStartTime = omp_get_wtime();

#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < dimension; ++j)
        {
            nests[i][j] = distribution(generator);
        }
        fitness[i] = objectiveFunction(nests[i]);
    }

    logToFile(logFilename, "Population initialization completed");

    double initEndTime = omp_get_wtime();

    bestSolution = nests[0];
    bestFitness = fitness[0];

    double totalGenerationTime = 0.0;
    double totalReplacementTime = 0.0;

    for (int generation = 0; generation < maxGenerations; ++generation)
    {
        double generationStartTime = omp_get_wtime();

#pragma omp parallel for
        for (int i = 0; i < n; ++i)
        {
            std::vector<double> newSolution(dimension);
            for (int j = 0; j < dimension; ++j)
            {
                newSolution[j] = nests[i][j] + levyFlight(lambda) * (nests[i][j] - nests[generator() % n][j]);
            }
            double newFitness = objectiveFunction(newSolution);

#pragma omp critical
            {
                if (newFitness < fitness[i])
                {
                    nests[i] = newSolution;
                    fitness[i] = newFitness;
                    if (newFitness < bestFitness)
                    {
                        bestSolution = newSolution;
                        bestFitness = newFitness;
                    }
                }
            }
        }

        double generationEndTime = omp_get_wtime();
        totalGenerationTime += (generationEndTime - generationStartTime);

        double replacementStartTime = omp_get_wtime();

#pragma omp parallel for
        for (int i = int(n * pa); i < n; ++i)
        {
            for (int j = 0; j < dimension; ++j)
            {
                nests[i][j] = distribution(generator);
            }
            fitness[i] = objectiveFunction(nests[i]);
        }

        double replacementEndTime = omp_get_wtime();
        totalReplacementTime += (replacementEndTime - replacementStartTime);
    }

    logToFile(logFilename, "Initialization time: " + std::to_string(initEndTime - initStartTime) + " seconds");
    logToFile(logFilename, "Total generation time: " + std::to_string(totalGenerationTime) + " seconds");
    logToFile(logFilename, "Total replacement time: " + std::to_string(totalReplacementTime) + " seconds");
}

int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        std::cerr << "Usage: " << argv[0] << " <number_of_threads> <benchmark_function> <log_filename>" << std::endl;
        return 1;
    }

    int numThreads = std::stoi(argv[1]);
    omp_set_num_threads(numThreads);

    std::string benchmarkFunction = argv[2];
    std::string logFilename = argv[3];
    double (*objectiveFunction)(const std::vector<double> &);
    int dimension = 30; // Default dimension

    if (benchmarkFunction == "alpine")
    {
        objectiveFunction = alpineFunction;
    }
    else if (benchmarkFunction == "ackley")
    {
        objectiveFunction = ackleyFunction;
    }
    else if (benchmarkFunction == "sphere")
    {
        objectiveFunction = sphereFunction;
    }
    else if (benchmarkFunction == "schwefel")
    {
        objectiveFunction = schwefelFunction;
    }
    else if (benchmarkFunction == "rastrigin")
    {
        objectiveFunction = rastriginFunction;
        dimension = 60; // Set dimension to 60 for these functions
    }
    else if (benchmarkFunction == "griewank")
    {
        objectiveFunction = griewankFunction;
        dimension = 60;
    }
    else if (benchmarkFunction == "csendes")
    {
        objectiveFunction = csendesFunction;
        dimension = 60;
    }
    else if (benchmarkFunction == "colville")
    {
        objectiveFunction = colvilleFunction;
        dimension = 4;
    }
    else if (benchmarkFunction == "easom")
    {
        objectiveFunction = easomFunction;
        dimension = 2;
    }
    else if (benchmarkFunction == "michalewicz")
    {
        objectiveFunction = michalewiczFunction;
        dimension = 5;
    }
    else if (benchmarkFunction == "shekel")
    {
        objectiveFunction = shekelFunction;
        dimension = 4;
    }
    else if (benchmarkFunction == "schwefel2_4")
    {
        objectiveFunction = schwefel2_4Function;
        dimension = 256;
    }
    else if (benchmarkFunction == "schwefel2_6")
    {
        objectiveFunction = schwefel2_6Function;
        dimension = 256;
    }
    else if (benchmarkFunction == "schaffer")
    {
        objectiveFunction = schafferFunction;
        dimension = 256;
    }
    else if (benchmarkFunction == "sumSquares")
    {
        objectiveFunction = sumSquaresFunction;
        dimension = 256;
    }
    else if (benchmarkFunction == "step2")
    {
        objectiveFunction = step2Function;
        dimension = 256;
    }
    else if (benchmarkFunction == "quartic")
    {
        objectiveFunction = quarticFunction;
        dimension = 256;
    }
    else if (benchmarkFunction == "powell")
    {
        objectiveFunction = powellFunction;
        dimension = 256;
    }
    else if (benchmarkFunction == "rosenbrock")
    {
        objectiveFunction = rosenbrockFunction;
        dimension = 256;
    }
    else if (benchmarkFunction == "dixonPrice")
    {
        objectiveFunction = dixonPriceFunction;
        dimension = 256;
    }
    else if (benchmarkFunction == "schwefel1_2")
    {
        objectiveFunction = schwefel1_2Function;
        dimension = 256;
    }
    else if (benchmarkFunction == "schwefel2_21")
    {
        objectiveFunction = schwefel2_21Function;
        dimension = 256;
    }
    else if (benchmarkFunction == "schwefel2_22")
    {
        objectiveFunction = schwefel2_22Function;
        dimension = 256;
    }
    else
    {
        std::cerr << "Invalid benchmark function." << std::endl;
        return 1;
    }

    int n = 50;               // Veličina populacije
    int maxGenerations = 500; // Maksimalni broj generacija
    double pa = 0.25;         // Verovatnoća zamene loših rešenja
    double lambda = 1.5;      // Parametar za Levy let

    std::vector<double> bestSolution;
    double bestFitness;

    // Merenje vremena
    double startTime = omp_get_wtime();
    std::cerr << "Starting Cuckoo Search\n";
    logToFile(logFilename, "Starting Cuckoo Search");
    parallelCuckooSearch(n, maxGenerations, pa, lambda, bestSolution, bestFitness, objectiveFunction, dimension, logFilename);
    double endTime = omp_get_wtime();
    double executionTime = endTime - startTime;

    std::cerr << "Cuckoo Search completed\n";
    logToFile(logFilename, "Cuckoo Search completed");

    std::string bestSolutionStr = "Najbolje rešenje: ";
    for (double xi : bestSolution)
    {
        bestSolutionStr += std::to_string(xi) + " ";
    }
    logToFile(logFilename, bestSolutionStr);

    logToFile(logFilename, "Najbolja vrednost ciljne funkcije: " + std::to_string(bestFitness));
    logToFile(logFilename, "Vreme izvršavanja: " + std::to_string(executionTime) + " sekundi");

    std::cout << bestSolutionStr << std::endl;
    std::cout << "Najbolja vrednost ciljne funkcije: " << bestFitness << std::endl;
    std::cout << "Vreme izvršavanja: " << executionTime << " sekundi" << std::endl;

    return 0;
}
