By completing the programming assignment, students will demonstrate the following skills:
- Profile a given program to determine the  execution times of the different functions and their time fractions of the program overall execution time
- Use the profiling results to identify the functions that represent significant fractions of the program execution time (more than 20%)
- Apply different low-level optimization techniques to improve the performance of the functions representing a significant fraction
- Determine the speedup of the fractions after applying the optimization techniques on the fractions
- Test the correctness of the optimized functions against the original functions
- Discuss and compare the results of the different optimization techniques
- Determine the cache performance of a given code

# Description
You are provided with a C++ program that generates four fractal images and then processes each image in three stages: Gamma Application phase, Blue Tinting phase, and Image Blurring phase. Each phase is performed by a separate function (applyGamma, applyTint, and applyBlur).
The program generates an output image at the end of the three phases for each fractal image.

In this assignment, you are asked to :
1. Profile the code in prog3.cpp to determine the execution time of the three functions applyGamma, applyTint, applyBlur. Determine the fraction of time that each function represents from the execution time to process one image (execution time of one iteration of the function process_images).
2. Identify which of the three functions represent important fractions of the  process_images execution time (> 20%) 
3. Apply different optimization techniques to the functions identified in 2. You must create  six (6) different versions of the functions using a combination of one or more technique  from the following list: 
..1. code motion
..2. reducing the number of function calls
..3. eliminating memory references when possible
..4. loop unrolling with different factors
..5. loop unrolling with different factors and multiple accumulators
..6.  loop unrolling with different factors, multiple accumulators, and re-association
4. Evaluate the effect of each optimized function by determining its speedup and the overall speedup on the time to process one image. 
5. Check the correctness of each optimized function against the original function by comparing the output images of the two functions. The program should display the result of the correctness check for each optimized function.
6. Discuss the results obtained by comparing the performance results obtained for the different optimized functions and the limitations if any. 
 
You are provided with the following files:
- `image.h`: header file where the classes Image and Fractal are declared
- `image.cpp`:  file where the classes Image and Fractal are implemented
- `prog3.cpp`: file where the functions applyGamma, applyTint, applyBlur, writeImage, process_images, and the main function are defined
- `makefile`: file to compile and generate the executable prog3 using the files image.cpp and prog3.cpp

You should modify `prog3.cpp` to complete the following:
1. Add the definitions of six optimized functions
2. Call the added functions and measure their execution time
3. Make sure that any added function passes the correctness test by comparing its output to the output of the original function. Use the function `checkCorrectness`

Your submission should consist of the updated `prog3.cpp` and a write up of two pages max that includes the following:
1. A description of  the optimization techniques you used for each function, a presentation of the results obtained, and a discussion of the results
2. The cache performance (miss rate) of the functions applyGamma, applyTint, and applyBlur if the cache is direct-mapped, the size of the cache line is 64 bytes, and the size of the cache is 32KB.

There is no autograding for this assignment since you have the choice for the type of optimization you will use.

