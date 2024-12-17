#include <iostream>
#include <random>
#include <ctime>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
using namespace std;

float randomFloat(float min, float max)
{
    return ((float)rand() / RAND_MAX) * (max - min) + min;
};

float **create2DArray(int N, int M)
{
    float **array = new float *[N];

    for (int i = 0; i < N; i++)
    {
        array[i] = new float[M];
    }
    return array;
};

void delete2DArray(float **array)
{
    delete[] array;
};

void print2DArray(float **array, int N, int M)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            cout << array[i][j] << " ";
        };
        cout << endl;
    };
};

float **random2DArray(int N, int M)
{
    float **array = create2DArray(N, M);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            array[i][j] = randomFloat(0, 1);
        };
    };

    return array;
};

std::vector<std::vector<float>> readCSV(const std::string &fileName)
{
    std::ifstream file(fileName);

    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file.\n";
        return {}; // Return an empty vector
    }

    std::vector<std::vector<float>> data; // 2D vector to store CSV data
    std::string line;

    while (std::getline(file, line))
    {
        std::stringstream ss(line);
        std::vector<float> row;
        std::string value;

        // Extract values separated by commas
        while (std::getline(ss, value, ','))
        {
            try
            {
                // Convert each value to float and push to row
                row.push_back(std::stof(value)); // std::stof converts string to float
            }
            catch (const std::invalid_argument &e)
            {
                std::cerr << "Warning: Invalid value in CSV, skipping: " << value << "\n";
                row.push_back(0); // If conversion fails, add a default value (0.0f)
            }
        }

        data.push_back(row); // Add the row to the 2D vector
    }

    file.close(); // Close the file
    return data;  // Return the 2D vector
}

void print2DVector(std::vector<std::vector<float>> data)
{
    for (int i = 0; i < data.size(); i++)
    {
        for (int j = 0; j < data[i].size(); j++)
        {
            cout << data[i][j] << " ";
        }
        cout << endl;
    }
}

std::vector<std::vector<float>> normaliseArray(std::vector<std::vector<float>> array)
{
    int numColumns = array[0].size();
    for (int col = 0; col < numColumns; col++)
    {

        // Find min and max vals of each column
        float maxVal = array[0][col];
        float minVal = array[0][col];

        for (const auto &row : array)
        {
            if (row[col] < minVal)
                minVal = row[col];
            if (row[col] > maxVal)
                maxVal = row[col];
        }

        // Normalise each element in array
        for (auto &row : array)
        {

            if (minVal != maxVal)
            {
                row[col] = (row[col] - minVal) / (maxVal - minVal);
            }
            else
            {
                row[col] = 0;
            }
        }
    }
    return array;
}

float computeVectorDistance(std::vector<float> x1, float *x2)
{
    float distance = 0;

    for (int i = 0; i < x1.size(); i++)
    {
        distance += pow((x1[i] - x2[i]), 2);
    }

    return sqrt(distance);
}

float monteCarloIntegral(std::vector<std::vector<float>> data, float radius, int N)
{
    srand(time(0));

    float insideTheSphere = 0.0;
    float outsideTheSphere = 0.0;

    std::vector<std::vector<float>> dataNorm = normaliseArray(data);

    int statistics[N][dataNorm.size()] = {0};

    for (int i = 0; i < N; i++)
    {
        float **randomSample = random2DArray(dataNorm.size(), dataNorm[0].size());

        for (int j = 0; j < dataNorm.size(); j++)
        {
            float distance = computeVectorDistance(dataNorm[j], randomSample[j]);
            if (distance <= radius)
                statistics[i][j] = 1;
        }
    }

    for (auto &row : statistics)
    {
        for (int i = 0; i < dataNorm.size(); i++)
        {
            if (row[i] == 1)
                insideTheSphere += 1.0;
            if (row[i] == 0)
                outsideTheSphere += 1.0;
        }
    }
    cout << insideTheSphere << endl;
    cout << outsideTheSphere << endl;

    float explorationRatio = abs(1 - (outsideTheSphere / insideTheSphere));

    return explorationRatio;
}

int main()
{

    std::vector<std::vector<float>> data = readCSV("/home/ignas/Documents/chemspx_corrections/ChemSPX/examples/2d/2d_samples.csv");
    float explorationRatio = monteCarloIntegral(data, 1, 10000);
    cout << explorationRatio << endl;

    return 0;
}