#include <stdio.h>
#include <algorithm>
#include <pthread.h>
#include <math.h>

#include <time.h>
#include "square_root_ispc.h"
#include "square_root_serial.h"
using namespace ispc;

extern void square_root_serial(float* input, float* output, int count, int maxIter);

static void verifyResult(int N, float* result, float* gold) {
    for (int i=0; i<N; i++) {
        if (fabs(result[i] - gold[i]) > 1e-4) {
            printf("Error: [%d] Got %f expected %f\n", i, result[i], gold[i]);
        }
    }
}

int main() {

    const unsigned int count = 20 * 1000 * 1000;
    const float initialGuess = 1.0f;
    clock_t startTime,endTime;
    float* input = new float[count];
    float* output = new float[count];
    float* gold = new float[count];
    const unsigned int maxIter = 32;
    for (unsigned int i=0; i<count; i++)
    {
        // random input values
        input[i] = .001f + 8.998f * static_cast<float>(rand()) / RAND_MAX;
        // TODO: Try different input values here.
        output[i] = 0.f;
    }

    // generate a gold version to check results
    for (unsigned int i=0; i<count; i++)
        gold[i] = sqrt(input[i]);

    //
    // And run the serial implementation 3 times, again reporting the
    // minimum time.
    //
    double minSerial = 1e30;
    for (int i = 0; i < 5; ++i) {
        double startTime = clock();
        square_root_serial(input, output, count, maxIter);
        double endTime = clock();
        minSerial = std::min(minSerial, (endTime - startTime)/CLOCKS_PER_SEC);
    }

    printf("[sqrt serial]:\t\t[%.3f] ms\n", minSerial * 1000);

    verifyResult(count, output, gold);

    //
    // Compute the image using the ispc implementation; report the minimum
    // time of three runs.
    //
    double minISPC = 1e30;
    for (int i = 0; i < 5; ++i) {
        double startTime = clock();
        square_root_serial(input, output, count, maxIter);
        double endTime = clock();
        minISPC = std::min(minISPC, (endTime - startTime)/CLOCKS_PER_SEC);
    }

    printf("[sqrt ispc]:\t\t[%.3f] ms\n", minISPC * 1000);

    verifyResult(count, output, gold);

    // Clear out the buffer
    for (unsigned int i = 0; i < count; ++i)
        output[i] = 0;

    //
    // Tasking version of the ISPC code
    //
    double minTaskISPC = 1e30;
    for (int i = 0; i < 3; ++i) {
        double startTime = clock();
        square_root_serial(input, output, count, maxIter);
        double endTime = clock();
        minTaskISPC = std::min(minTaskISPC, (endTime - startTime)/CLOCKS_PER_SEC);
    }

    printf("[sqrt task ispc]:\t[%.3f] ms\n", minTaskISPC * 1000);

    verifyResult(count, output, gold);

    printf("\t\t\t\t(%.2fx speedup from ISPC)\n", minSerial/minISPC);
    printf("\t\t\t\t(%.2fx speedup from task ISPC)\n", minSerial/minTaskISPC);

    delete[] input;
    delete[] output;
    delete[] gold;

    return 0;
}
