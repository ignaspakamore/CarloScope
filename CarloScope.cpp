#include <iostream>
#include <random>
#include <ctime>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>
#include <array>
using namespace std;

float randomFloat(float min, float max)
{
    return ((float)rand() / RAND_MAX) * (max - min) + min;
};

double roundToTwoDecimalPlaces(double value) { 
    return std::round(value * 100.0) / 100.0; 
} 

double **create2DArray(int N, int M)
{
    double **array = new double *[N];

    for (int i = 0; i < N; i++)
    {
        array[i] = new double[M];
    }
    return array;
};

void delete2DArray(double **array)
{
    delete[] array;
};

void print2DArray(double **array)
{   int N = sizeof(array);
    int M = sizeof(array[0]);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            cout << array[i][j] << " ";
        };
        cout << endl;
    };
};

double **random2DArray(int N, int M)
{
    double **array = create2DArray(N, M);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            array[i][j] = randomFloat(0, 1);
        };
    };

    return array;
}

double** readCSV(string *&filename) {
    std::ifstream file(*filename);
    int cols = 0;
    int rows = 0;

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return nullptr;
    }

    std::vector<std::vector<double>> data;
    std::string line;

    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream lineStream(line);
        std::string cell;

        while (std::getline(lineStream, cell, ',')) {
            try {row.push_back(std::stod(cell));} // Convert cell to double and add to row
            catch(const std::invalid_argument& e){}
        }

        if (row.size() > cols) {
            cols = row.size(); // Update the number of columns
        }

        data.push_back(row);
    }

    file.close();

    rows = data.size();

    // Allocate 2D array
    double** array = new double*[rows];

    for (size_t i = 0; i < rows; ++i) {
        array[i] = new double[cols](); // Initialize with zeros
        for (size_t j = 0; j < data[i].size(); ++j) {
            array[i][j] = data[i][j];
        }
    }

    return array;
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

double **normaliseArray(double **array)
{   
    int numColumns = sizeof(array[0]);
    int numRows = sizeof(array);

    for (int col = 0; col < numColumns; col++)
    {

        // Find min and max vals of each column
        float maxVal = array[0][col];
        float minVal = array[0][col];

        for (int row=0; row<numRows; row++)
        {
            if (array[row][col] < minVal)
                minVal = array[row][col];
            if (array[row][col] > maxVal)
                maxVal = array[row][col];
        }

        // Normalise each element in array
        for (int row=0; row<numRows; row++)
        {

            if (minVal != maxVal)
            {
                array[row][col] = (array[row][col] - minVal) / (maxVal - minVal);
            }
            else
            {
                array[row][col] = 0;
            }
        }
    }
    return array;
}

double computeVectorDistance(double *x1, double *x2)
{
    double distance = 0;

    for (int i = 0; i < sizeof(x1); i++)
    {
        distance += pow((x1[i] - x2[i]), 2);
    }

    return sqrt(distance);
}

double monteCarloIntegral(double **data, double *&radius, long long *&N)
{
    srand(time(0));

    float insideTheSphere = 0.0;
    float outsideTheSphere = 0.0;
    
    double **dataNorm = normaliseArray(data);

    int sampleSize = sizeof(dataNorm[0]);

    double **randomArray = random2DArray(*N, sampleSize);

    for (int i=0; i < *N; i++)
    {
        for (int j=0; j < sampleSize; j++)
        {
            double distance = computeVectorDistance(randomArray[i], dataNorm[j]);

            if (distance < *radius)
            {
                insideTheSphere += 1.0;
            }

            else if (distance > *radius)
            {
                outsideTheSphere += 1.0;
            }
            
        }

    }

    double explorationRatio = insideTheSphere / outsideTheSphere;

    explorationRatio = roundToTwoDecimalPlaces(explorationRatio);

    return explorationRatio;
}

void printHelp(){
    cout<<"Usage program [options]\n"
        <<"Options:"
        <<" -help   Shows this help message\n"
        <<" -file   Input CSV file\n"
        <<" -r      Radius\n"
        <<" -n      Number of Monte Carlo samples\n";
}


void getArguments(int argc, char *argv[], string *&path, double *&radius, long long *&N){

    for (int i=1; i<argc; i++){
        string arg = argv[i];

        if (arg == "-help"){
            printHelp();
        }

        else if (arg == "-file"){
            if (i+1 < argc){
                path = new string(argv[i+1]);
            }
            else{
                cout<<"Error: -file requires a CSV file path argument.\n";
                exit(1);
            }
        }

        else if (arg == "-r"){
            if(i+1 < argc){
               radius = new double(std::stod(argv[i+1]));
            }
            else{
                cout<<"Error: -r requires a float radius value.\n";
                exit(1);
            }
        }

        else if (arg == "-n"){
            if(i+1 < argc){
                N = new long long(std::stoll(argv[i+1]));
            }
            else{
                cout<<"Error: -n requires integer N Monte Carlo sample number.\n";
                exit(1);
            }
        }
    }

    if (*path == "" && *radius ==0 && *N == 0)
    {
        cout<<"Error: One or more input parameters are NULL!\n";
        exit(1);
    }
    return;
}

int main(int argc, char *argv[])
{   
    string* path=nullptr;
    double* radius=nullptr;
    long long* N=nullptr;

    cout << R"(
_________             ______     ________                         
__  ____/_____ __________  /_______  ___/________________________ 
_  /    _  __ `/_  ___/_  /_  __ \____ \_  ___/  __ \__  __ \  _ \
/ /___  / /_/ /_  /   _  / / /_/ /___/ // /__ / /_/ /_  /_/ /  __/
\____/  \__,_/ /_/    /_/  \____//____/ \___/ \____/_  .___/\___/ 
                                                    /_/                                                     
    )" <<endl;

    getArguments(argc, argv, path, radius, N);
 
    double **data = readCSV(path);
    double explorationRatio = monteCarloIntegral(data, radius, N);

    cout<<"Design Space:\n"
        <<"Occupied: "<<explorationRatio<<endl
        <<"Free:     "<<1-explorationRatio<<endl;

    return 0;
}