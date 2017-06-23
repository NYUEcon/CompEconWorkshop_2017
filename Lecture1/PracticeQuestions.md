# Computational Economics: Tools 1 Practice

## Getting to know your computer

1. Evaluate `"1101011101010000"` as a floating point number
2. Given a matrix `A = [5 3; 0 5]` (semi-colon denotes next row). Compute the eigenvalues of `A`. What is the condition number of `A`?
3. Create a large matrix of normally distributed random numbers (for example, `1000 x 1000`). Write two for loops (first/outer loop iterating over the elements of the column and second loop iterating over the elements of the row) and time how long it takes to sum up all of the elements in this matrix. Now reverse the order of the loops and time it. Which is faster? What does this mean about whether the programming language is column or row major? If you do this in Python then you will need to speed up the loops using `numba`.

## Root finding and optimization

1. Write a function that computes the root of a one-dimensional function.
2. Write a function that computes the optimal value of a one-dimensional function.
3. There are two countries. One country, A, produces apples and the other country, B, produces bananas. In each country there is a single representative consumer who likes to consume bundles of apples and bananas. The consumer in country A has a utility function `u_A = (w*a^x + (1-w)*b^x)^(1/x)` whereas the consumer in country B has utility function `u_B = (w*b^x + (1-w)*a^x)^(1/x)`. Country A produces 5 apples and country B produces 10 bananas. Both countries have a home bias parameter (`w`) of `0.9` and an elasticity of substitution parameter (`x`) of `-0.5`. How many apples/bananas does each country demand? What is the relative price of bananas to apples? What happens to the price if the number of bananas produced drops to 5?

