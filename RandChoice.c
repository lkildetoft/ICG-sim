#include <stdlib.h>
#include "RandChoice.h"

float randChoice(float *randlist, float *weights, int len) {
    /*
    Generate weighted random numbers by calculating the cumulative
    sum of the weights and finding the weight which is greater or
    equal to the generated random number.

    Arguments:
    *randlist: array of floats containing numbers to choose from.
    *weigths: array of floats containing the weight of each number
    in the distribution. Must be indexed and equal in length to 
    randlist. 
    len: the length of the arrays. 

    Returns:
    A weigthed random pick.
    */
    int low = 0;
    int mid = 0;
    float wsum = 0;
    float *cweights = calloc(len, sizeof(float));
    //Calculate cumulative sums of weights
    if (cweights == NULL) {
        return -1;
    } else {
        for (int i = 0; i < len; i++) {
            wsum += weights[i];
            cweights[i] = wsum;
        }
    };
    //Generate a random number between 0 and wsum
    float rnd = ((float)rand()/(RAND_MAX))*wsum;
    //Find the weighted choice by binary search
    while (low < len) {
        mid = (low + len - 1)/2;
        if (rnd == cweights[mid]) {
            break;
        } else if (rnd > cweights[mid]) {
            low = mid + 1;
        } else {
            len = mid - 1;
        }
    };

    free(cweights);

    return randlist[mid];
};