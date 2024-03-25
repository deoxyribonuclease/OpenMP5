#include <iostream>
#include <omp.h>
#include <climits>
#include <cstdlib>
#include <ctime>
#include <thread>
#include <chrono>

const int rows = 10000;
const int cols = 10000;

int arr[rows][cols];

void init_arr() {
    srand(time(0));
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            arr[i][j] = rand() % 100;
        }
    }
}

long long get_sum(int num_threads, double& execution_time) {
    long long sum = 0;
    double t1 = omp_get_wtime();

#pragma omp parallel num_threads(num_threads) reduction(+:sum)
    {
#pragma omp for
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                sum += arr[i][j];
            }
        }
    }
    double t2 = omp_get_wtime();
    execution_time = t2 - t1;
    return sum;
}

int get_min_row(int num_threads, int& min_row, int& min_sum, double& execution_time) {
    min_sum = INT_MAX;
    min_row = -1;
    double t1 = omp_get_wtime();

#pragma omp parallel for num_threads(num_threads)
    for (int i = 0; i < rows; i++) {
        long long row_sum = 0;
        for (int j = 0; j < cols; j++) {
            row_sum += arr[i][j];
        }
        if (row_sum < min_sum) {
#pragma omp critical
            {
                if (row_sum < min_sum) {
                    min_sum = row_sum;
                    min_row = i;
                }
            }
        }
    }
    double t2 = omp_get_wtime();
    execution_time = t2 - t1;
}

int main() {
    int const max_thread = 8;
    long long sum[max_thread];
    int min_sum[max_thread];
    int min_row[max_thread];
    double executed_time[2][max_thread];
    init_arr();


    omp_set_nested(1);
    double t1 = omp_get_wtime();
#pragma omp parallel sections
    {
#pragma omp section
        {
            for (int i = 1; i <= max_thread; i += 1)
            {
               sum[i-1] = get_sum(i,executed_time[0][i - 1]);
            }
        }
#pragma omp section
        {
            for (int i = 1; i <= max_thread; i += 1)
            {
                get_min_row(i, min_row[i - 1], min_sum[i - 1], executed_time[1][i - 1]);
            }
        }
    }
    double t2 = omp_get_wtime();

    for (int i = 1; i <= max_thread; i += 1)
    {
        std::cout << "Threads used: " << i << " all sum - " << sum[i - 1] << " and executed time - " 
            << executed_time[0][i - 1] << " seconds" << std::endl;
    }
      
    for (int i = 1; i <= max_thread; i += 1)
    {
        std::cout << "Threads used: " << i << " Min row - " << min_row[i - 1] << " with min sum " << min_sum[i - 1]
            << " and executed time - " << executed_time[1][i - 1] << " seconds" << std::endl;
    }
    std::cout << "Total time - " << t2 - t1 << " seconds" << std::endl;

    return 0;
}
