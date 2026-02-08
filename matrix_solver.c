/*** Yiğit Baki ***/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*** MATRIX STRUCTURE ***/
/* A matrix structure that is suitable for dynamic matrix operations */
typedef struct{
    float **M;  // Pointer for a 2-dimensional (2D) float type array
    int n, m;   // Numbers of rows and columns for the 2D array
} Matrix;


/*** FUNCTION PROTOTYPES ***/
Matrix construct(Matrix x, int n, int m);   // Returns a Matrix variable by allocating memory for n-by-m sized 2D array
Matrix destruct(Matrix x);                  // Returns a Matrix variable by deallocating the memory allocated for the 2D array
Matrix configure(Matrix x);                 // Returns a Matrix variable by initializing its members with the user-specified parameters
void print(Matrix x);                       // Prints the 2D array elements held in a Matrix variable
float determinant(Matrix x);                // Returns the calculated determinant value for a Matrix variable
Matrix cofactor(Matrix x);                  // Returns a Matrix variable that holds the cofactors of the 2D array of a Matrix variable
Matrix transpose(Matrix x);                 // Returns a Matrix variable that holds the transpose of the 2D array of a Matrix variable
Matrix inverse(Matrix x);                   // Returns a Matrix variable that holds the inverse of the 2D array of a Matrix variable
Matrix multiply(Matrix x, Matrix y);        // Returns a Matrix variable that holds the multiplication of 2D array of 2 Matrix variables
Matrix solve(Matrix A, Matrix b);           // Returns the solution in Matrix type for a linear equations system specified by 2 Matrix variables



/*** MAIN FUNCTION ***/
// Define and solve a linear equations system (Ax=b)
int main(){
    Matrix A = configure(A);                                        // Define and configure the coefficient matrix (n-by-n)
    Matrix b = configure(b);                                        // Define and configure the constant (right-hand side) vector (n-by-1)
    Matrix x = solve(A, b);                                         // Solve the system and store the solution in the variable vector (n-by-1)
    getchar();                                                      // Press enter to see the solution
    for(int i=0; i<x.n; i++) printf("x%d = %g\t", i+1, x.M[i][0]);  // Print the solution
    putchar('\n');                                                  // Give a new line after printing the solution
    getchar();                                                      // Press enter to terminate the program
    return 0;
}

/*** FUNCTION IMPLEMENTATIONS ***/

/* Returns a Matrix variable by allocating memory for n-by-m sized 2D array */
Matrix construct(Matrix x, int n, int m){
    // Initialize the numbers of rows and columns for array
    x.n = n;
    x.m = m;
    // Allocate memory for rows of array
    x.M = (float **)malloc(n * sizeof(float *));
    // Allocate memory for array elements in each row
    int i;
    for ( i = 0; i < n; i++)
    {
        x.M[i] = (float *)malloc(m * sizeof(float));
    }

    return x;       // Return the constructed Matrix argument
}


/* Returns a Matrix variable by deallocating the memory allocated for the 2D array */
Matrix destruct(Matrix x){
    // Deallocate the memory allocated for columns of array for each row
    // Deallocate the memory allocated for rows of array
    int i;
    for ( i = 0; i < x.n; i++)
    {
        free(x.M[i]);
    }
    free(x.M);
    x.M = NULL;
    x.n=0;
    x.m=0;
    
    return x;                               // Return the destructed Matrix argument
}


/* Returns a Matrix variable by initializing its members with the user-specified parameters */
Matrix configure(Matrix x){
    // Echo instructions that "Enter the matrix size as numbers of rows and columns: "
    printf("Enter the matrix size as numbers of rows and columns: ");
    // Read numbers of rows and columns into members n an m respectively
    scanf("%d %d", &x.n, &x.m);
    // Call construct() with numbers of rows and columns read for Matrix argument
    x = construct(x, x.n, x.m);
    // Echo instructions that "Enter the matrix elements row by row:\n"
    printf("Enter the matrix elements row by row:\n");
    // Read matrix element values one-by-one based iterating through rows and columns
    int i, j;
    for (i = 0; i < x.n; i++)
    {
        for ( j = 0; j < x.m; j++)
        {
            scanf("%f", &x.M[i][j]);
        }
        
    }
    // Put a new line character after reading matrix
    printf("\n");

    return x;                                                                   // Return the configured Matrix argument
}


/* Prints the 2D array elements held in a Matrix variable */
void print(Matrix x){
    // Put a new line character at beginning
    printf("\n");
    // Print each element row-by-row and separating columns by tab characters and put a new line character after each row
    int i, j;
    for ( i = 0; i < x.n; i++)
    {
        for ( j = 0; j < x.m; j++)
        {
            printf("%f\t", x.M[i][j]);
        }
        printf("\n");
    }
    
    return;  // Return nothing
}


/* Returns the calculated determinant value for a Matrix variable */
float determinant(Matrix x){
    float det;                                           // Define a float variable for the determinant value
    int n = x.n;                                         // Define matrix size with a single variable by assuming a square matrix
    
    
    if (x.n == x.m){    // If the matrix is square then start to calculate the determinant value
        if (n <= 0){    // If the matrix size is not positive then echo an error that "\nThe matrix does not exist!\n" and immediately return with the uninitialized determinant value
        printf("\nThe matrix does not exist!\n");
        return det;
        }

        if (n == 1){ // If the matrix has a single element then the determinant value will be equal to its value        
        return x.M[0][0];
        }

        if (n == 2){ // If the matrix size is 2-by-2 then the determinant value will be calculated by immediate calculation       
            return x.M[0][0] * x.M[1][1] - x.M[0][1] * x.M[1][0]; 
        }
    
        if (n > 2) // If the matrix is larger then calculate the determinant via triple nested loops (O(n^3))
        {
            int i; 
            for ( i = 0; i < n; i++){
            Matrix minor = construct(minor , n-1, n-1); // Define a minor matrix by reducing the the original matrix size by 1
        
            int j;
            for ( j = 1; j < n; j++)
            {
                int col = 0;    // Define a manual index to iterate the columns of the minor matrix
                int k;
                for (k = 0; k < n; k++)
                {
                    if (k != i)     // Iterate through rows and then through columns for each row (except the first one) of the original matrix and also skip the column if it corresponds the element that the minor is looked for
                    {
                        minor.M[j-1][col] = x.M[j][k];  // Assign the original matrix elements to the minor matrix by handling all indices by incrementing and resetting the manuel column iterator of the minor matrix as well as the regular loop indicies
                        col++;
                    }
                    
                }    
            }
            det += pow(-1, i) * x.M[0][i] * determinant(minor);     // Accumulate the minor determinants calculated via recursive call of determinant()
            minor = destruct(minor);    
        }
        return det;   
        }
    }
        
    if (x.n != x.m) // If the matrix is not square then echo an error that "\nThe matrix is non-square and has no determinant!\n"
    {
        printf("\nThe matrix is non-square and has no determinant!\n");
        return det;
    }    
    return det;          // Return the calculated determinant value
}


/* Returns a Matrix variable that holds the cofactors of the 2D array of a Matrix variable*/
Matrix cofactor(Matrix x){
    int n = x.n;                                                        // Define the matrix size with a single variable by assuming a square matrix
    Matrix cfc = construct(cfc, n, n);                                  // Define and construct a square Matrix variable for the cofactor matrix
    // If the matrix is square then start to calculate the cofactor matrix
    if (x.n == x.m)
    {
        if (n <= 0){     // If the matrix size is not positive then echo an error that "\nThe matrix does not exist!\n" and immediately return with the uninitialized cofactor matrix        
            printf("\nThe matrix does not exist!\n");
            return cfc;
        }

        if (n == 1){     // If the matrix has single element then its cofactor becomes 1        
            cfc.M[0][0] = 1;
            return cfc;
        }

        if (n > 1)      // If the matrix is larger, then obtain cofactors via quadruple nested loops (O(n^4))
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Matrix minor = construct(minor, n-1, n-1);      // Define a minor matrix by reducing the the original matrix size by 1
                    int row = 0;                                    // Define two manual indices to iterate the rows and columns of the minor matrix
                    for (int k = 0; k < n; k++)
                    {
                        if (k != i)
                        {
                            int col = 0;
                            for (int l = 0; l < n; l++)
                            {
                                if (l != j)     // Then reiterate through rows and then through columns for each row of the original matrix by skipping the row or the column if it corresponds the element that the minor is looked for
                                {
                                    minor.M[row][col] = x.M[k][l];  // Calculate minor matrices for each row and then each column for each row of the original matrix
                                    col++;
                                }   
                            }  
                            row++; 
                        }   
                    }
                    cfc.M[i][j] = pow(-1, i+ j) * determinant(minor);        // Calculate values of elements of the cofactor matrix by obtaining determinant values for obtained minor matrices (SEE (2) IN HOMEWORK DOCUMENT.)
                    destruct(minor);  
                }   
            }   
        }    
        return cfc;
    }
    
    if (x.n != x.m) // If the matrix is not square then echo an error that "\nThe matrix is non-square and has no cofactor!\n"
    {
        printf("\nThe matrix is non-square and has no cofactor!\n");
        return cfc;
    }
    return cfc;             // Return the Matrix variable that holds the cofactors
}


 /* Returns a Matrix variable that holds the transpose of the 2D array of a Matrix variable */
Matrix transpose(Matrix x){
    Matrix tr = construct(tr, x.n, x.m); // REWRITE THIS LINE and both define and construct (in a single line) a reversed-sized (columns-by-rows) Matrix variable for transpose
    for (int i = 0; i < x.n; i++)
    {
        for (int j = 0; j < x.m; j++)   // Assign rows of the original matrix to columns of the transpose matrix
        {
            tr.M[j][i] = x.M[i][j];
        }   
    }
    
    return tr;                                         // Return the transposed Matrix variable
}


/* Returns a Matrix variable that holds the inverse of the 2D array of a Matrix variable */
Matrix inverse(Matrix x){
    int n = x.n;                                                                        // Define the matrix size with a single variable by assuming a square matrix
    float det = determinant(x);                                                         // Calculate the determinant value of the matrix
    Matrix inv = construct(inv, n, n); // both define and construct a square Matrix variable for the inverse matrix

    if (det == 0){    // If the the matrix determinant is 0, then echo an error that "\nThe matrix is non-square or singular and has no inverse!\n"
        printf("\nThe matrix is non-square or singular and has no inverse!\n");
        return inv;
    }else{      // If the matrix determinant is nonzero then calculate the inverse matrix
        Matrix cfc = cofactor(x);
        Matrix adjoint = transpose(cfc);    // Define and obtain the adjoint matrix as the transpose of the cofactor matrix
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                inv.M[i][j] = adjoint.M[i][j] / det;    // Calculate the inverse matrix elements by dividing adjoints with the determinant value by iterating through rows and columns
            }   
        }
        destruct(cfc);
        destruct(adjoint);      // Adjoint matrix may be destructed after all the inverse matrix elements calculated
        return inv;              // Return the inverted Matrix variable
    }                                                        
}


/* Returns a Matrix variable that holds the multiplication of 2D arrays of 2 Matrix variables */
Matrix multiply(Matrix x, Matrix y){
    Matrix mlt = construct(mlt, x.n, y.m); // both define and construct a proper-sized Matrix variable for the matrix multiplication
    if (x.m != y.n)     // If sizes of matrices does not match then echo an error that "\nSizes of matrices do not match for the matrix multiplication!\n"
    {
        printf("\nSizes of matrices do not match for the matrix multiplication!\n");
        return mlt;
    }

    if (x.m == y.n)     // If sizes of the matrices match for the matrix multiplication
    {
        int n = x.n;     // Iterate through the rows and columns of the multiplication matrix
        int m = y.m;

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                mlt.M[i][j] = 0;         // Initialize each element of the multiplication matrix as 0 by default and then assign sum of products to elements of the multiplication matrix by iterating both matrix arguments matching the rows of first multiplicand matrix to the columns of second multiplicand matrix
                for (int k = 0; k < x.m; k++)
                {
                    mlt.M[i][j] += x.M[i][k] * y.M[k][j];
                }
                
            }
            
        }
        return mlt;         // Return the Matrix variable that holds the multiplication of matrices
    }
    return mlt;
}  
/* Returns the solution in Matrix type for a linear equations system specified by 2 Matrix variables */
Matrix solve(Matrix A, Matrix b){
    Matrix invA = inverse(A);
    Matrix x = multiply(invA, b);
    destruct(invA);
    return x;       // Return the calculated solution for the linear equations system in Matrix type via multiply() and inverse()
}
