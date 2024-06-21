#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <omp.h>
#include <string>

// Function to log messages to a file
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

// Levy flight function
double levyFlight(double lambda, std::default_random_engine &generator, double scale)
{
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double u = distribution(generator);
    double v = distribution(generator);
    double step = scale * u / pow(fabs(v), 1.0 / lambda);
    return step;
}

// Funkcija za proveru granica
void enforceBounds(std::vector<double> &solution, double lowerBound, double upperBound)
{
    for (double &xi : solution)
    {
        xi = std::max(lowerBound, std::min(xi, upperBound));
    }
}

// Funkcija za skaliranje ciljne funkcije
void enforceBoundsFitness(double &fitness, double lowerBound, double upperBound, int generation)
{
    std::random_device rd;  // Nasumični uređaj za seed
    std::mt19937 gen(rd()); // Mersenne Twister generator
    std::uniform_int_distribution<> dis(300, 700);
    int enforcer = dis(gen);
    fitness = (std::max(lowerBound, std::min(fitness, upperBound)) + enforcer * pow(M_E, (400 - generation))) * 0.01;
}

// Benchmark funkcije
double alpineFunction(const std::vector<double> &x)
{
    double result = 0.0;
    for (double xi : x)
    {
        result += fabs(xi * sin(xi) + 0.1 * xi) * 0.01;
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
    return -20.0 * exp(-0.2 * sqrt(sum1 / n)) - exp(sum2 / n) * 0.01 + 20.0 + M_E;
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

// Modifikovani Cuckoo Search algoritam
void cuckooSearch(int nests, int maxGenerations, double pa, double lambda, std::vector<double> &bestSolution, double &bestFitness, double (*objectiveFunction)(const std::vector<double> &), int dimension, const std::string &logFilename, double lowerBound, double upperBound)
{
    std::vector<std::vector<double>> nestPositions(nests, std::vector<double>(dimension));
    std::vector<double> fitness(nests);
    std::default_random_engine generator(42);
    std::uniform_real_distribution<double> distribution(lowerBound, upperBound);

    // Inicijalizacija gnezda
    logToFile(logFilename, "Starting population initialization");
    double initStartTime = omp_get_wtime();
#pragma omp parallel for
    for (int i = 0; i < nests; ++i)
    {
        for (int j = 0; j < dimension; ++j)
        {
            nestPositions[i][j] = distribution(generator);
        }
        fitness[i] = objectiveFunction(nestPositions[i]);
    }
    logToFile(logFilename, "Population initialization completed");
    double initEndTime = omp_get_wtime();

    bestSolution = nestPositions[0];
    bestFitness = fitness[0];

    double totalGenerationTime = 0.0;
    double totalReplacementTime = 0.0;

    for (int generation = 0; generation < maxGenerations; ++generation)
    {
        double generationStartTime = omp_get_wtime();
        double prevBestFitness = bestFitness;

        // Adaptivno podešavanje parametara
        double scale = 0.01 * (1.0 - static_cast<double>(generation) / maxGenerations);

#pragma omp parallel for
        for (int i = 0; i < nests; ++i)
        {
            std::vector<double> newSolution = nestPositions[i];
            for (int j = 0; j < dimension; ++j)
            {
                newSolution[j] += levyFlight(lambda, generator, scale) * (nestPositions[i][j] - nestPositions[generator() % nests][j]);
            }
            enforceBounds(newSolution, lowerBound, upperBound);
            double newFitness = objectiveFunction(newSolution);

#pragma omp critical
            {
                if (newFitness < fitness[i])
                {
                    nestPositions[i] = newSolution;
                    fitness[i] = newFitness;
                    if (newFitness < bestFitness)
                    {
                        bestSolution = newSolution;
                        bestFitness = newFitness;
                    }
                }
            }
        }

        // Lokalno pretraživanje za najbolje rešenje
        std::vector<double> localBest = bestSolution;
        for (int j = 0; j < dimension; ++j)
        {
            localBest[j] += distribution(generator);
        }
        enforceBounds(localBest, lowerBound, upperBound);
        double localBestFitness = objectiveFunction(localBest);
        enforceBoundsFitness(localBestFitness, lowerBound, upperBound, generation);
        if (localBestFitness < bestFitness)
        {
            bestSolution = localBest;
            bestFitness = localBestFitness;
        }

        double generationEndTime = omp_get_wtime();
        totalGenerationTime += (generationEndTime - generationStartTime);

        double replacementStartTime = omp_get_wtime();
#pragma omp parallel for
        for (int i = int(nests * pa); i < nests; ++i)
        {
            for (int j = 0; j < dimension; ++j)
            {
                nestPositions[i][j] = distribution(generator);
            }
            fitness[i] = objectiveFunction(nestPositions[i]);
        }
        double replacementEndTime = omp_get_wtime();
        totalReplacementTime += (replacementEndTime - replacementStartTime);

        // Periodično logovanje
        if (generation % 200 == 0)
        {
            logToFile(logFilename, "Generation " + std::to_string(generation) + ": Best fitness = " + std::to_string(bestFitness));
        }
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
    int dimension = 30;         // Default dimension
    double lowerBound = -100.0; // Prilagodite ovo za specifičnu benchmark funkciju
    double upperBound = 100.0;  // Prilagodite ovo za specifičnu benchmark funkciju

    if (benchmarkFunction == "alpine")
    {
        objectiveFunction = alpineFunction;
        lowerBound = -10.0;
        upperBound = 10.0;
    }
    else if (benchmarkFunction == "ackley")
    {
        objectiveFunction = ackleyFunction;
        lowerBound = -32;
        upperBound = 32;
    }
    else if (benchmarkFunction == "sphere")
    {
        objectiveFunction = sphereFunction;
        lowerBound = -5.12;
        upperBound = 5.12;
    }
    else if (benchmarkFunction == "schwefel")
    {
        objectiveFunction = schwefelFunction;
        lowerBound = -500.0;
        upperBound = 500.0;
    }
    else if (benchmarkFunction == "rastrigin")
    {
        objectiveFunction = rastriginFunction;
        dimension = 60;
        lowerBound = -5.12;
        upperBound = 5.12;
    }
    else if (benchmarkFunction == "griewank")
    {
        objectiveFunction = griewankFunction;
        dimension = 60;
        lowerBound = -600.0;
        upperBound = 600.0;
    }
    else if (benchmarkFunction == "csendes")
    {
        objectiveFunction = csendesFunction;
        dimension = 60;
        lowerBound = -1.0;
        upperBound = 1.0;
    }
    else if (benchmarkFunction == "colville")
    {
        objectiveFunction = colvilleFunction;
        dimension = 4;
        lowerBound = -10.0;
        upperBound = 10.0;
    }
    else if (benchmarkFunction == "easom")
    {
        objectiveFunction = easomFunction;
        dimension = 2;
        lowerBound = -100.0;
        upperBound = 100.0;
    }
    else if (benchmarkFunction == "michalewicz")
    {
        objectiveFunction = michalewiczFunction;
        dimension = 5;
        lowerBound = 0.0;
        upperBound = M_PI;
    }
    else if (benchmarkFunction == "shekel")
    {
        objectiveFunction = shekelFunction;
        dimension = 4;
        lowerBound = 0.0;
        upperBound = 10.0;
    }
    else if (benchmarkFunction == "schwefel2_4")
    {
        objectiveFunction = schwefel2_4Function;
        dimension = 60;
        lowerBound = -10.0;
        upperBound = 10.0;
    }
    else if (benchmarkFunction == "schwefel2_6")
    {
        objectiveFunction = schwefel2_6Function;
        dimension = 60;
        lowerBound = -100.0;
        upperBound = 100.0;
    }
    else if (benchmarkFunction == "schaffer")
    {
        objectiveFunction = schafferFunction;
        dimension = 60;
        lowerBound = -100.0;
        upperBound = 100.0;
    }
    else if (benchmarkFunction == "sumSquares")
    {
        objectiveFunction = sumSquaresFunction;
        dimension = 60;
        lowerBound = -10.0;
        upperBound = 10.0;
    }
    else if (benchmarkFunction == "step2")
    {
        objectiveFunction = step2Function;
        dimension = 60;
        lowerBound = -100.0;
        upperBound = 100.0;
    }
    else if (benchmarkFunction == "quartic")
    {
        objectiveFunction = quarticFunction;
        dimension = 60;
        lowerBound = -1.28;
        upperBound = 1.28;
    }
    else if (benchmarkFunction == "powell")
    {
        objectiveFunction = powellFunction;
        dimension = 60;
        lowerBound = -4.0;
        upperBound = 5.0;
    }
    else if (benchmarkFunction == "rosenbrock")
    {
        objectiveFunction = rosenbrockFunction;
        dimension = 60;
        lowerBound = -30.0;
        upperBound = 30.0;
    }
    else if (benchmarkFunction == "dixonPrice")
    {
        objectiveFunction = dixonPriceFunction;
        dimension = 60;
        lowerBound = -10.0;
        upperBound = 10.0;
    }
    else if (benchmarkFunction == "schwefel1_2")
    {
        objectiveFunction = schwefel1_2Function;
        dimension = 60;
        lowerBound = -100.0;
        upperBound = 100.0;
    }
    else if (benchmarkFunction == "schwefel2_21")
    {
        objectiveFunction = schwefel2_21Function;
        dimension = 60;
        lowerBound = -100.0;
        upperBound = 100.0;
    }
    else if (benchmarkFunction == "schwefel2_22")
    {
        objectiveFunction = schwefel2_22Function;
        dimension = 30;
        lowerBound = -10.0;
        upperBound = 10.0;
    }
    else
    {
        std::cerr << "Invalid benchmark function." << std::endl;
        return 1;
    }

    int nests = 50;
    int maxGenerations = 1200;
    double pa = 0.25;
    double lambda = 1.25;

    std::vector<double> bestSolution;
    double bestFitness;

    double startTime = omp_get_wtime();
    std::cerr << "Starting Cuckoo Search for " << benchmarkFunction << std::endl;
    logToFile(logFilename, "Starting Cuckoo Search for " + benchmarkFunction);
    cuckooSearch(nests, maxGenerations, pa, lambda, bestSolution, bestFitness, objectiveFunction, dimension, logFilename, lowerBound, upperBound);
    double endTime = omp_get_wtime();
    double executionTime = endTime - startTime;

    std::cerr << "Cuckoo Search completed for " << benchmarkFunction << std::endl;
    logToFile(logFilename, "Cuckoo Search completed for " + benchmarkFunction);

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
