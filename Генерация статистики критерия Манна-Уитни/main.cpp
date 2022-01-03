#include <iostream>
#include <ctime>
#include <random>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string>
#include <chrono>

using namespace std;

// Для генерации распределений внутри функций
static std::random_device rd; // random device engine, usually based on /dev/random on UNIX-like systems
// initialize Mersennes' twister using rd to generate the seed
static std::mt19937 rng{rd()};

double statisticsForNLittle(vector<double> &sample1, vector<double> &sample2)
{
    // Заполнение общего вектора для дальнейшего нахождения рангов
    vector<vector<double>> R;
    R.resize(sample1.size() + sample2.size());
    for(int i = 0; i < R.size(); i++)
        R[i].resize(2);

    for(int i = 0; i < sample1.size(); i++)
        R[i][0] = sample1[i];

    for(int i = sample1.size(), j = 0; i < R.size(); i++, j++)
        R[i][0] = sample2[j];

    sort(R.begin(), R.end());

    int count = 0;
    double sum = 0;

    // Расчёт ранга последних элементов при условии, что они равны друг другу
    for(int i = R.size()-1; R[i][0] == R[i-1][0]; i--)
        count++;

    // Если в конце вектора есть элементы, равные друг другу, находится их ранг
    if(count > 0)
    {
        // Находится сумма рангов этих элементов
        for(int i = R.size()-1, j = 0; j <= count; i--, j++)
            sum += i+1;

        // Среднее арифмитическое (count + 1, т.к. если одинаковых элемента 2, то count, согласно коду выше, будет равен 1, и из-за этого предпоследний элемент не инициализируется)
        for(int i = R.size()-1, j = 0; j <= count; i--, j++)
            R[i][1] = sum / (count + 1);
    }

    // Иначе последний элемент получает самое высокое значение ранга
    else R[R.size()-1][1] = R.size();

    // Временная переменная для определения индекса элемента, до которого нужно инициализировать оставшиеся элементы
    int tmp = count;

    count = 0;
    sum = 0;

    // Инициализация рангов всех элементов кроме последнего (последних)
    for(int i = 0; i < R.size() - 1 - tmp; i++)
    {
        // Если текущий элемент не равен предыдущему (о чём говорит count) и не равен следующему элементам, то присваивается ранг по порядку
        if(count == 0 && R[i][0] != R[i+1][0])
            R[i][1] = i+1;
        // Иначе
        else
        {
            // Если текущий элемент равен следующему, увеличивается счётчик
            if(R[i][0] == R[i+1][0])
                count++;
            // Иначе...
            else
            {
                // Счётчик увеличивается в последний раз для определения всех одинаковых элементов
                count++;
                for(int j = 0; j < count; j++)
                    sum += i-j;

                for(int j = 0; j < count; j++)
                    R[i-j][1] = sum / count;

                sum = 0;
                count = 0;
            }
        }
    }

    vector<double> vTmp;
    vTmp.resize(R.size());
    for(int i = 0; i < vTmp.size(); i++)
    {
        vTmp[i] = R[i][0];
    }

    // Ранг 1-й выборки
    int R1 = 0;
    // Ранг 2-й выборки
    int R2 = 0;
    // Переменные для нахождения индекса элемента по значению
    vector<double>::iterator it;
    int index = 0;

    // Расчёт 1-го ранга
    for(int i = 0; i < sample1.size(); i++)
    {
        // Поиск элемента по значению
        it = find(vTmp.begin(), vTmp.end(), sample1[i]);
        // Нахождение индекса элемента
        index = distance(vTmp.begin(), it);

        // Если элемент есть в общем векторе, ранг увеличивается
        if(index < vTmp.size())
            R1 += R[index][1];
    }

    //cout << "R1 = " << R1 << endl;

    // Расчёт 2-го ранга
    for(int i = 0; i < sample2.size(); i++)
    {
        it = find(vTmp.begin(), vTmp.end(), sample2[i]);
        index = distance(vTmp.begin(), it);

        if(index < vTmp.size())
            R2 += R[index][1];
    }

    //cout << R1 + R2 << endl;

    //cout << "R2 = " << R2 << endl;

    int m = sample1.size();
    int n = sample2.size();

    double U1 = m * n + ((m * (m - 1)) / 2) - R1;
    //cout << "U1 = " << U1 << endl;
    double U2 = m * n + ((n * (n - 1)) / 2) - R2;
    //cout << "U2 = " << U2 << endl;

    return min(U1, U2);
}

double MannWhitneyUtest(vector<double> &sample1, vector<double> &sample2)
{
    int m = sample1.size();
    int n = sample2.size();

    if(m + n > 2 && m >= 1 && n >= 1)
    {
        double sum = 0;
        for(int i = 0; i < sample1.size(); i++)
        {
            for(int j = 0; j < sample2.size(); j++)
            {
                if(sample2[j] < sample1[i])
                    sum += 1;
                else if(sample2[j] == sample1[i])
                    sum += 0.5;
            }
        }
        return ((abs(sum - ((m * n) / 2))) / sqrt((m * n * (m + n + 1)) / 12));
    }

    else return statisticsForNLittle(sample1, sample2);
}

vector<double> normRGenerate(int n, double D = 1)
{
    normal_distribution<double> normR(0, D);

    vector<double> sample;
    sample.resize(n);
    for(int i = 0; i < sample.size(); i++)
        sample[i] = normR(rng);

    return sample;
}

vector<double> cauchyRGenerate(int n, double D = 1)
{
    cauchy_distribution<double> cauchyR(0, D);

    vector<double> sample;
    sample.resize(n);
    for(int i = 0; i < sample.size(); i++)
        sample[i] = cauchyR(rng);

    return sample;
}

vector<double> expRGenerate(int n, double D = 1)
{
    exponential_distribution<double> expR(D);

    vector<double> sample;
    sample.resize(n);
    for(int i = 0; i < sample.size(); i++)
        sample[i] = expR(rng);

    return sample;
}

vector<double> MonteKarloMethod(int N, int m, int n, string typeOfRaspr = "norm", int varDif = 1)
{
    vector<double> statSample;
    vector<double> sample1;
    vector<double> sample2;

    if(typeOfRaspr == "norm")
    {
        for(int i = 0; i < N; i++)
        {
            sample1 = normRGenerate(m);
            sample2 = normRGenerate(n, varDif);
            statSample.push_back(MannWhitneyUtest(sample1, sample2));
        }
    }

    else if(typeOfRaspr == "cauchy")
    {
        for(int i = 0; i < N; i++)
        {
            sample1 = cauchyRGenerate(m);
            sample2 = cauchyRGenerate(n, varDif);
            statSample.push_back(MannWhitneyUtest(sample1, sample2));
        }
    }

    else
    {
        for(int i = 0; i < N; i++)
        {
            sample1 = expRGenerate(m);
            sample2 = expRGenerate(n, varDif);
            statSample.push_back(MannWhitneyUtest(sample1, sample2));
        }
    }

    return statSample;
}

void samplesGenerate(int N, int n, int m, string typeOfRaspr = "norm", int varDif = 1)
{
    vector<double> sample;
    sample = MonteKarloMethod(N, m, n, typeOfRaspr, varDif);

    string sampleName = typeOfRaspr + to_string(n) + "x" + to_string(m) + "_var" + to_string(varDif) + ".dat";
    ofstream sampleStream(sampleName);
    sampleStream << sampleName << endl;
    sampleStream << "0 " << N << endl;
    for (int i = 0; i < N; i++)
    {
        sampleStream << sample[i] << endl;
    }
}

int main()
{
    int n, m, N;

    N = 100;

    cout << "Enter n: ";
    cin >> n;
    m = n;
    cout << "Enter N: ";
    cin >> N;
    string dist;
    cout << "Enter the distribution (norm, cauchy or exp): ";
    cin >> dist;
    int varDif;
    cout << "Enter the variance difference: ";
    cin >> varDif;
    samplesGenerate(N, n, m, dist, varDif);

    return 0;
}
