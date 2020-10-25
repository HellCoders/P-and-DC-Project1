
export void square_root_ispc(uniform float input[], uniform float output[], uniform int count, uniform int maxIter) {
    
    //X is te variable which will output the square root and S is the input real number
    float X = 1, S;

    //Using foreach for SPMD implementation within the core
    foreach (i = 0 ... count) {
        S = input[i];
        for (uniform int j = 0; j < maxIter; j++) {
            X = 0.5 * (X + (S/X));
        }
        output[i] = X;
    }

}


task void square_root_ispc_task(uniform int span, uniform float input[], uniform float output[], uniform int count, uniform int maxIter) {

    //X is te variable which will output the square root and S is the input real number
    float X = 1, S;

    uniform int indexStart = taskIndex * span;
    uniform int indexEnd = min(count, indexStart + span);

    foreach (i = indexStart ... indexEnd) {
        S = input[i];
        for (uniform int j = 0; j < maxIter; j++) {
            X = 0.5 * (X + (S/X));
        }
        output[i] = X;
    }

}


export void square_root_ispctask(uniform float input[], uniform float output[], uniform int count, uniform int maxIter)
{

    uniform int span = count / 64;  // 64 tasks

    launch[count/span] square_root_ispc_task(span, input, output, count, maxIter);
}
