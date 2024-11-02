#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <strings.h>
#include <string.h>

#define MAX_ROWS 100
#define MAX_COLS 100

int rows, cols;                           // dimensions of the grid.
 
// structure to hold the minimum and maximum thresholds
typedef struct {
    double min;
    double max;
} Threshold;

// structure for a 2D segment tree to store statistical information of the data grid
typedef struct {
    int rows, cols;
    double ***minMatrix; // 3D matrix to hold minimum values for each parameter
    double ***maxMatrix; // 3D matrix to hold maximum values for each parameter
    double ***sumMatrix; // 3D matrix to hold sum values for each parameter
} SegmentTree2D;

// structure for holding the primary data matrices
typedef struct {
    int rows, cols;
    double **temp;       // matrix for temperature values
    double **humidity;   // matrix for humidity values
    double **rain;       // matrix for rain values
} Data;

// structure to define a query with various parameters
typedef struct {
    int x1, y1, x2, y2;       // coordinates for the query range
    int row, col;      // specific row, column, and value for point queries
    double value;
    char rangeType[10];       // type of range for the query (e.g., min, max, sum)
    int valid;                // flag to indicate if the query is valid
    double result;            // result of the query
    double timeTaken;         // time taken to execute the query
    int attribute;            // specifies the attribute (e.g., temp, humidity, rain) for the query
} Query;


// helper function to allocate 2D Segment Tree matrices
SegmentTree2D *create2DSegmentTree(int rows, int cols) {
    // allocate memory for the SegmentTree2D structure
    SegmentTree2D *tree = malloc(sizeof(SegmentTree2D));
    tree->rows = rows; // set number of rows
    tree->cols = cols; // set number of columns

    // allocate memory for min, max, and sum matrices, each having 4*rows
    tree->minMatrix = (double***)malloc(((4 * rows) + 1) * sizeof(double**));
    tree->maxMatrix = (double***)malloc(((4 * rows) + 1) * sizeof(double**));
    tree->sumMatrix = (double***)malloc(((4 * rows) + 1) * sizeof(double**));

    // loop through each row of the 3D matrices
    for (int k = 0; k < (4 * rows) + 1; ++k) {
        // allocate memory for each row's columns in the matrices
        tree->minMatrix[k] = (double**)malloc(((4 * cols) + 1) * sizeof(double*));
        tree->maxMatrix[k] = (double**)malloc(((4 * cols) + 1) * sizeof(double*));
        tree->sumMatrix[k] = (double**)malloc(((4 * cols) + 1) * sizeof(double*));

        // loop through each column in the current row
        for (int i = 0; i < (4 * cols) + 1; ++i) {
            // allocate space for 3 parameters (e.g., temperature, humidity, rain) in each cell
            tree->minMatrix[k][i] = (double*)malloc(3 * sizeof(double));
            tree->maxMatrix[k][i] = (double*)malloc(3 * sizeof(double));
            tree->sumMatrix[k][i] = (double*)malloc(3 * sizeof(double));

            // initialize the values for each parameter
            for (int j = 0; j < 3; j++) {
                // set initial minimum values to the largest possible double
                tree->minMatrix[k][i][j] = DBL_MAX;
                // set initial maximum values to the smallest possible double
                tree->maxMatrix[k][i][j] = -DBL_MAX;
                // initialize sum values to zero
                tree->sumMatrix[k][i][j] = 0.0;
            }
        }
    }

    return tree; // return the pointer to the created segment tree
}

Data *createDataMatrix(int rows, int cols) {
    // allocate memory for the Data structure
    Data *data = (Data *)malloc(sizeof(Data));
    data->rows = rows; // set the number of rows
    data->cols = cols; // set the number of columns

    // allocate memory for temperature, humidity, and rain matrices
    data->temp = (double **)malloc(rows * sizeof(double *));
    data->humidity = (double **)malloc(rows * sizeof(double *));
    data->rain = (double **)malloc(rows * sizeof(double *));

    // loop through each row to allocate columns and initialize values
    for (int i = 0; i < rows; i++) {
        // allocate memory for each column in the matrices
        data->temp[i] = (double *)malloc(cols * sizeof(double));
        data->humidity[i] = (double *)malloc(cols * sizeof(double));
        data->rain[i] = (double *)malloc(cols * sizeof(double));

        // initialize each cell to 0.0
        for (int j = 0; j < cols; j++) {
            data->temp[i][j] = 0.0;       // set initial temperature to 0.0
            data->humidity[i][j] = 0.0;   // set initial humidity to 0.0
            data->rain[i][j] = 0.0;       // set initial rain value to 0.0
        }
    }

    return data; // return the pointer to the initialized data matrix
}

// helper function to free the Segment Tree matrices
void free2DSegmentTree(SegmentTree2D *tree) {
    // free each allocated level in the segment tree matrices
    for (int k = 0; k < (4 * tree->rows) + 1; ++k) {
        for (int i = 0; i < (4 * tree->cols) + 1; ++i) {
            free(tree->minMatrix[k][i]);
            free(tree->maxMatrix[k][i]);
            free(tree->sumMatrix[k][i]);
        }
        free(tree->minMatrix[k]);
        free(tree->maxMatrix[k]);
        free(tree->sumMatrix[k]);
    }

    // free the top-level pointers for minMatrix, maxMatrix, and sumMatrix
    free(tree->minMatrix);
    free(tree->maxMatrix);
    free(tree->sumMatrix);

    // finally, free the SegmentTree2D structure itself
    free(tree);
}

// helper function to free the allocated memory in the Data struct
void free_data(Data *data, int rows) {
    for (int i = 0; i < rows; i++) {
        free(data->temp[i]);
        free(data->humidity[i]);
        free(data->rain[i]);
    }
    free(data->temp);
    free(data->humidity);
    free(data->rain);
}

// utility function to get maximum of two numbers.
int max(int a, int b)
{
    return (a > b) ? a : b;
}

// utility function to get minimum of two numbers.
int min(int a, int b)
{
    return (a < b) ? a : b;
}

// row: index in the segment tree for rows.
// lx, rx: range of rows being covered by the current node.
// col: index in the segment tree for columns.
// ly, ry: range of columns being covered by the current node.

// function to build the column segment tree for a particular row.
void buildRow(Data* data, SegmentTree2D *segTree, int row, int lx, int rx, int col, int ly, int ry, Threshold* thresholds , int** marked) {
    // case for a single column.
    if (ly == ry) {
        // case for a single row.
        if (lx == rx) {
            // Initialize sum and max for each parameter with the cell value.
            segTree->sumMatrix[row][col][0] = data->temp[lx][ly];
            segTree->maxMatrix[row][col][0] = data->temp[lx][ly];
            segTree->minMatrix[row][col][0] = data->temp[lx][ly];
            segTree->sumMatrix[row][col][1] = data->humidity[lx][ly];
            segTree->maxMatrix[row][col][1] = data->humidity[lx][ly];
            segTree->minMatrix[row][col][1] = data->humidity[lx][ly];
            segTree->sumMatrix[row][col][2] = data->rain[lx][ly];
            segTree->maxMatrix[row][col][2] = data->rain[lx][ly];
            segTree->minMatrix[row][col][2] = data->rain[lx][ly];

            // set marked array based on the temperature threshold.
            if((data->temp[lx][ly] > thresholds[0].min && data->temp[lx][ly] < thresholds[0].max) && (data->humidity[lx][ly] > thresholds[1].min && data->humidity[lx][ly] < thresholds[1].max)&&(data->rain[lx][ly] > thresholds[2].min && data->rain[lx][ly] < thresholds[2].max)) marked[lx][ly] = 1;
        }
        // case when not a single row.
        else {
            // merge the max and sum values for the row's children for each parameter.
            for (int i = 0; i < 3; i++) {
                segTree->sumMatrix[row][col][i] = segTree->sumMatrix[2 * row][col][i] + segTree->sumMatrix[2 * row + 1][col][i];
                segTree->maxMatrix[row][col][i] = fmax(segTree->maxMatrix[2 * row][col][i], segTree->maxMatrix[2 * row + 1][col][i]);
                segTree->minMatrix[row][col][i] = fmin(segTree->minMatrix[2 * row][col][i], segTree->minMatrix[2 * row + 1][col][i]);
            }
        }
    }
    // case when not a single column.
    else {
        // recursively build for each half.
        int my = (ly + ry) / 2; // middle column index for splitting.
        buildRow(data, segTree, row, lx, rx, 2 * col, ly, my, thresholds , marked);        // build left child.
        buildRow(data, segTree, row, lx, rx, 2 * col + 1, my + 1, ry, thresholds , marked); // build right child.

        // merge values from left and right children for each parameter.
        for (int i = 0; i < 3; i++) {
            segTree->sumMatrix[row][col][i] = segTree->sumMatrix[row][2 * col][i] + segTree->sumMatrix[row][2 * col + 1][i];
            segTree->maxMatrix[row][col][i] = fmax(segTree->maxMatrix[row][2 * col][i], segTree->maxMatrix[row][2 * col + 1][i]);
            segTree->minMatrix[row][col][i] = fmin(segTree->minMatrix[row][2 * col][i], segTree->minMatrix[row][2 * col + 1][i]);
        }
    }
}

// row: index in the segment tree for rows
// lx, rx: range of rows being covered by the current node

// function to build the entire 2D segment tree.
void build(Data* data, SegmentTree2D *segTree, int row, int lx, int rx, int cols, Threshold* thresholds , int **marked) {
    // Case when not a single row.
    if (lx != rx) {
        int mid = (lx + rx) / 2;                        // Calculate the middle row index.
        build(data, segTree, 2 * row, lx, mid, cols, thresholds , marked);         // Build left child for rows.
        build(data, segTree, 2 * row + 1, mid + 1, rx, cols, thresholds , marked); // Build right child for rows.
    }

    // Build the column segment for this row segment.
    buildRow(data, segTree, row, lx, rx, 1, 0, cols - 1, thresholds , marked);
}


// function to query the sum of a rectangular sub-region for a specific row segment
// parameters:
// - row: the index representing the current row segment in the segment tree
// - col: the index representing the current column segment in the segment tree
// - t_ly, t_ry: the range of columns that the current node covers
// - ly, ry: the range of columns we want to query for the sum
double sumQueryRow(SegmentTree2D *segTree, int row, int col, int t_ly, int t_ry, int ly, int ry , int attribute)
{
    // if the queried range is invalid (e.g., start is greater than end), return 0 as there's no sum
    if (ly > ry)
        return 0.0;

    // if the current segment matches exactly with the query range, return the sum stored in this node
    if (ly == t_ly && ry == t_ry)
    {
        return segTree->sumMatrix[row][col][attribute];
    }

    // find the middle column of the current segment
    int tmy = (t_ly + t_ry) / 2;

    // recursively query the left child (first half of the column range)
    // combine it with the result of querying the right child (second half of the column range)
    return sumQueryRow(segTree, row, 2 * col, t_ly, tmy, ly, min(ry, tmy) , attribute) +            // left child query
           sumQueryRow(segTree, row, 2 * col + 1, tmy + 1, t_ry, max(ly, tmy + 1), ry , attribute); // right child query
}

// function to query the sum of a rectangular sub-region in the 2D segment tree
// parameters:
// - row: the index representing the current row segment in the segment tree
// - t_lx, t_rx: the range of rows that the current node covers
// - lx, rx: the range of rows we want to query for the sum
// - ly, ry: the range of columns we want to query for the sum
double sumQuery(Data* data, SegmentTree2D *segTree, int row, int t_lx, int t_rx, int lx, int rx, int ly, int ry , int attribute)
{
    // if the queried range is invalid (e.g., start is greater than end), return 0 as there's no sum
    if (lx > rx)
        return 0.0;

    // if the current segment matches exactly with the row query range,
    // delegate the column query to sumQueryRow function
    if (lx == t_lx && rx == t_rx)
    {
        return sumQueryRow(segTree, row, 1, 0, cols - 1, ly, ry , attribute);
    }

    // find the middle row of the current segment
    int t_mid = (t_lx + t_rx) / 2;

    // recursively query the left child (first half of the row range)
    // combine it with the result of querying the right child (second half of the row range)
    return sumQuery(data, segTree, 2 * row, t_lx, t_mid, lx, min(rx, t_mid), ly, ry , attribute) +            // left child query
           sumQuery(data, segTree, 2 * row + 1, t_mid + 1, t_rx, max(lx, t_mid + 1), rx, ly, ry , attribute); // right child query
}

double avgQuery(Data* data, SegmentTree2D *segTree, int row, int t_lx, int t_rx, int lx, int rx, int ly, int ry, int attribute) {
    // call sumQuery to get the total sum of the attribute in the specified range
    double totalSum = sumQuery(data, segTree, row, t_lx, t_rx, lx, rx, ly, ry, attribute);

    // calculate the number of cells in the subgrid (rows * columns)
    int numCells = (rx - lx + 1) * (ry - ly + 1);

    // calculate and return the average value
    return totalSum / (double)numCells;
}

// function to query the maximum value of a rectangular sub-region for a specific row segment
// parameters:
// - row: the index representing the current row segment in the segment tree
// - col: the index representing the current column segment in the segment tree
// - t_ly, t_ry: the range of columns that the current node covers
// - ly, ry: the range of columns we want to query for the maximum value
double maxQueryRow(SegmentTree2D *segTree, int row, int col, int t_ly, int t_ry, int ly, int ry , int attribute)
{
    // if the queried range is invalid (e.g., start is greater than end), return INT_MIN as there's no valid value
    if (ly > ry)
        return DBL_MIN;

    // if the current segment matches exactly with the query range, return the maximum value stored in this node
    if (ly == t_ly && ry == t_ry)
    {
        return segTree->maxMatrix[row][col][attribute];
    }

    // find the middle column of the current segment
    int t_mid = (t_ly + t_ry) / 2;

    // recursively query the left child (first half of the column range)
    // combine it with the result of querying the right child (second half of the column range) to find the maximum
    return fmax(maxQueryRow(segTree, row, 2 * col, t_ly, t_mid, ly, min(ry, t_mid) , attribute),              // left child query
                maxQueryRow(segTree, row, 2 * col + 1, t_mid + 1, t_ry, max(ly, t_mid + 1), ry , attribute)); // right child query
}

// function to query the maximum value of a rectangular sub-region in the 2d segment tree
// parameters:
// - row: the index representing the current row segment in the segment tree
// - t_lx, t_rx: the range of rows that the current node covers
// - lx, rx: the range of rows we want to query for the maximum value
// - ly, ry: the range of columns we want to query for the maximum value
double maxQuery(SegmentTree2D *segTree, int row, int t_lx, int t_rx, int lx, int rx, int ly, int ry , int attribute)
{
    // if the queried range is invalid (e.g., start is greater than end), return INT_MIN as there's no valid value
    if (lx > rx)
        return DBL_MIN;

    // if the current segment matches exactly with the row query range,
    // delegate the column query to maxQueryRow function
    if (lx == t_lx && rx == t_rx)
    {
        return maxQueryRow(segTree, row, 1, 0, cols - 1, ly, ry , attribute);
    }

    // find the middle row of the current segment
    int t_mid = (t_lx + t_rx) / 2;

    // recursively query the left child (first half of the row range)
    // combine it with the result of querying the right child (second half of the row range) to find the maximum
    return fmax(maxQuery(segTree, 2 * row, t_lx, t_mid, lx, min(rx, t_mid), ly, ry , attribute),              // left child query
                maxQuery(segTree, 2 * row + 1, t_mid + 1, t_rx, max(lx, t_mid + 1), rx, ly, ry , attribute)); // right child query
}

double minQueryRow(SegmentTree2D *segTree, int row, int col, int t_ly, int t_ry, int ly, int ry , int attribute)
{
    // if the queried range is invalid (e.g., start is greater than end), return INT_MIN as there's no valid value
    if (ly > ry)
        return DBL_MAX;

    // if the current segment matches exactly with the query range, return the maximum value stored in this node
    if (ly == t_ly && ry == t_ry)
    {
        return segTree->minMatrix[row][col][attribute];
    }

    // find the middle column of the current segment
    int t_mid = (t_ly + t_ry) / 2;

    // recursively query the left child (first half of the column range)
    // combine it with the result of querying the right child (second half of the column range) to find the maximum
    return fmin(minQueryRow(segTree, row, 2 * col, t_ly, t_mid, ly, min(ry, t_mid) , attribute),              // left child query
                minQueryRow(segTree, row, 2 * col + 1, t_mid + 1, t_ry, max(ly, t_mid + 1), ry , attribute)); // right child query
}

// function to query the maximum value of a rectangular sub-region in the 2d segment tree
// parameters:
// - row: the index representing the current row segment in the segment tree
// - t_lx, t_rx: the range of rows that the current node covers
// - lx, rx: the range of rows we want to query for the maximum value
// - ly, ry: the range of columns we want to query for the maximum value
double minQuery(SegmentTree2D *segTree, int row, int t_lx, int t_rx, int lx, int rx, int ly, int ry , int attribute)
{
    // if the queried range is invalid (e.g., start is greater than end), return INT_MIN as there's no valid value
    if (lx > rx)
        return DBL_MAX;

    // if the current segment matches exactly with the row query range,
    // delegate the column query to maxQueryRow function
    if (lx == t_lx && rx == t_rx)
    {
        return minQueryRow(segTree, row, 1, 0, cols - 1, ly, ry , attribute);
    }

    // find the middle row of the current segment
    int t_mid = (t_lx + t_rx) / 2;

    // recursively query the left child (first half of the row range)
    // combine it with the result of querying the right child (second half of the row range) to find the maximum
    return fmin(minQuery(segTree, 2 * row, t_lx, t_mid, lx, min(rx, t_mid), ly, ry , attribute),              // left child query
                minQuery(segTree, 2 * row + 1, t_mid + 1, t_rx, max(lx, t_mid + 1), rx, ly, ry , attribute)); // right child query
}

// function to update a single point in the segment tree for a specific row segment
// parameters:
// - row: the index representing the current row segment in the segment tree
// - col: the index representing the current column segment in the segment tree
// - lx, rx: the range of rows being covered by the current node
// - ly, ry: the range of columns being covered by the current node
// - a, b: the coordinates (row, column) of the point we want to update
// - value: the new value to update the point to
void updateRow(Data* data, SegmentTree2D *segTree, int row, int lx, int rx, int col, int ly, int ry, int a, int b, double value , Threshold* thresholds , int** marked , int attribute)
{
    // if we're at a single column (a leaf node in the column segment tree)
    if (ly == ry)
    {
        // check if we're at a single row (a leaf node in the row segment tree)
        if (lx == rx)
        {
            // if both the row and column are leaf nodes, update the point directly
                segTree->sumMatrix[row][col][attribute] = value; // update the sum value at the point
                segTree->maxMatrix[row][col][attribute] = value;           // update the max value at the point
                segTree->minMatrix[row][col][attribute] = value;           // update the max value at the point
            if(value > thresholds[attribute].min && value < thresholds[attribute].max) marked[lx][ly] = 1;
            else marked[lx][ly] = 0;
        }
        else
        {
            // if it's not a single row, merge the updates from both child row nodes
                segTree->sumMatrix[row][col][attribute] = segTree->sumMatrix[2 * row][col][attribute] + segTree->sumMatrix[2 * row + 1][col][attribute]; // update the sum
                segTree->maxMatrix[row][col][attribute] = fmax(segTree->maxMatrix[2 * row][col][attribute], segTree->maxMatrix[2 * row + 1][col][attribute]);                          // update the max
                segTree->sumMatrix[row][col][attribute] = segTree->sumMatrix[2 * row][col][attribute] + segTree->sumMatrix[2 * row + 1][col][attribute]; // update the sum
                segTree->minMatrix[row][col][attribute] = fmin(segTree->minMatrix[2 * row][col][attribute], segTree->minMatrix[2 * row + 1][col][attribute]);                          // update the max
        }
    }
    else
    {
        // if it's not a single column, split the columns and update the relevant child
        int mid = (ly + ry) / 2;

        // if the column index is within the left half, update the left child
        if (b <= mid)
        {
            updateRow(data, segTree, row, lx, rx, 2 * col, ly, mid, a, b, value , thresholds , marked , attribute);
        }
        else
        {
            // otherwise, update the right child
            updateRow(data, segTree, row, lx, rx, 2 * col + 1, mid + 1, ry, a, b, value ,thresholds, marked , attribute);
        }

        // after updating the children, merge their sums and max values
            segTree->sumMatrix[row][col][attribute] = segTree->sumMatrix[row][2 * col][attribute] + segTree->sumMatrix[row][2 * col + 1][attribute]; // merge sums
            segTree->maxMatrix[row][col][attribute] = fmax(segTree->maxMatrix[row][2 * col][attribute], segTree->maxMatrix[row][2 * col + 1][attribute]);                          // merge maxes
            segTree->minMatrix[row][col][attribute] = fmin(segTree->minMatrix[row][2 * col][attribute], segTree->minMatrix[row][2 * col + 1][attribute]);                          // merge maxes
    }
}

// function to update a single point in the 2d segment tree
// parameters:
// - row: the index representing the current row segment in the segment tree
// - lx, rx: the range of rows being covered by the current node
// - a, b: the coordinates (row, column) of the point we want to update
// - value: the new value to update the point to
void update(Data* data, SegmentTree2D *segTree, int row, int lx, int rx, int a, int b, double value , Threshold *thresholds, int** marked , int attribute)
{
    // if it's not a single row (non-leaf node in the row segment tree)
    if (lx != rx)
    {
        int mid = (lx + rx) / 2;

        // if the row index is within the left half, update the left child
        if (a <= mid)
        {
            update(data, segTree, 2 * row, lx, mid, a, b, value , thresholds , marked , attribute);
        }
        else
        {
            // otherwise, update the right child
            update(data, segTree, 2 * row + 1, mid + 1, rx, a, b, value , thresholds ,  marked , attribute);
        }
    }

    // after updating the rows, update the column segment tree for this row segment
    updateRow(data, segTree, row, lx, rx, 1, 0, cols - 1, a, b, value ,thresholds, marked , attribute);
}

int *nearestSmallerToLeft(int *heights, int size, int pseudo)
{
    int *v = (int *)malloc(size * sizeof(int));
    int stack[size][2]; // stack to store elements and their indices
    int top = -1;

    for (int i = 0; i < size; i++)
    {
        // pop elements from stack that are greater than or equal to the current element
        while (top >= 0 && stack[top][0] >= heights[i])
        {
            top--;
        }
        // store the index of the nearest smaller element to the left, or -1 if none
        v[i] = (top == -1) ? pseudo : stack[top][1];

        // push current element and index to the stack
        stack[++top][0] = heights[i];
        stack[top][1] = i;
    }
    return v;
}

// function to find the nearest smaller elements to the right for each element in `heights`
int *nearestSmallerToRight(int *heights, int size, int pseudo)
{
    int *v1 = (int *)malloc(size * sizeof(int));
    int stack[size][2]; // Stack to store elements and their indices
    int top = -1;

    for (int i = size - 1; i >= 0; i--)
    {
        // Pop elements from stack that are greater than or equal to the current element
        while (top >= 0 && stack[top][0] >= heights[i])
        {
            top--;
        }
        // Store the index of the nearest smaller element to the right, or `pseudo` if none
        v1[i] = (top == -1) ? pseudo : stack[top][1];

        // Push current element and index to the stack
        stack[++top][0] = heights[i];
        stack[top][1] = i;
    }
    return v1;
}

// function to calculate the maximum area in a histogram-like row
int subregionArea(int *subRegions, int size, int *leftCoord, int *rightCoord, int row)
{
    int *left = nearestSmallerToLeft(subRegions, size, -1);     // Indices of nearest smaller elements to the left
    int *right = nearestSmallerToRight(subRegions, size, size); // Indices of nearest smaller elements to the right
    int max_area = 0;

    // Iterate through each element in the row to calculate potential area
    for (int i = 0; i < size; i++)
    {
        int width = right[i] - left[i] - 1; // Calculate width of the rectangle
        int area = subRegions[i] * width;   // Calculate area for this rectangle

        // Update max_area and coordinates if a new maximum area is found
        if (area > max_area)
        {
            max_area = area;
            *leftCoord = left[i] + 1;   // Left column of the rectangle
            *rightCoord = right[i] - 1; // Right column of the rectangle
        }
    }

    // Free allocated memory for left and right arrays
    free(left);
    free(right);
    return max_area;
}

// Main function to find the maximal rectangular region in a binary matrix
int maximalRegion(int **region, int x, int y, int *topLeftX, int *topLeftY, int *bottomRightX, int *bottomRightY)
{
    int *temp = (int *)malloc(y * sizeof(int)); // Temp array to store heights of histogram columns

    // Initialize temp with the first row of region
    for (int j = 0; j < y; j++)
    {
        temp[j] = region[0][j];
    }

    // Initialize coordinates for the maximal region
    *topLeftX = *topLeftY = *bottomRightX = *bottomRightY = 0;

    int maxArea = subregionArea(temp, y, topLeftY, bottomRightY, 0); // Calculate initial max area for the first row

    // Process each subsequent row to calculate the maximal rectangular region
    for (int i = 1; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            // Update temp array based on the current row; reset height if 0, otherwise accumulate height
            temp[j] = (region[i][j] == 0) ? 0 : temp[j] + region[i][j];
        }

        // Track the left-most and right-most columns of the max area subregion
        int leftCoord, rightCoord;
        int area = subregionArea(temp, y, &leftCoord, &rightCoord, i);

        // Update max area and coordinates if a new maximum area is found
        if (area > maxArea)
        {
            maxArea = area;
            *topLeftX = i - temp[leftCoord] + 1; // Top row of the rectangle
            *topLeftY = leftCoord;               // Left column
            *bottomRightX = i;                   // Bottom row
            *bottomRightY = rightCoord;          // Right column
        }
    }

    free(temp); // Free allocated memory for temp array
    return maxArea;
}

// function to check if a given cell (r, c) can be included in BFS.
bool isSafe(int **M, int r, int c, int ROW, int COL)
{
    return (r >= 0) && (r < ROW) && (c >= 0) &&
           (c < COL) && (M[r][c] == 1);
}

// breadth-first search to visit all cells in the current connected region (island).
void findRegion(int **M, int sr, int sc, int ROW, int COL)
{
    // defining row and column offsets for the 8 neighboring cells.
    int dr[] = {-1, -1, -1, 0, 0, 1, 1, 1};
    int dc[] = {-1, 0, 1, -1, 1, -1, 0, 1};

    // initializing a queue to store coordinates for BFS.
    int queue[ROW * COL][2];
    int front = 0, rear = 0; // setting up front and rear pointers for the queue.
    queue[rear][0] = sr;     // enqueue the starting row.
    queue[rear][1] = sc;     // enqueue the starting column.
    rear++;                  // moving the rear pointer forward.
    M[sr][sc] = 0;           // marking the starting cell as visited by setting it to 0.

    // continue until there are no more cells to process in the queue.
    while (front < rear)
    {
        int r = queue[front][0]; // get the current row index from the front of the queue.
        int c = queue[front][1]; // get the current column index from the front of the queue.
        front++;                 // move to the next cell in the queue.

        // exploring all 8 neighboring cells.
        for (int k = 0; k < 8; k++)
        {
            int row = r + dr[k]; // calculate the new row index.
            int col = c + dc[k]; // calculate the new column index.
            if (isSafe(M, row, col, ROW, COL))
            {                         // check if the new cell is safe to visit.
                queue[rear][0] = row; // enqueue the new row index.
                queue[rear][1] = col; // enqueue the new column index.
                rear++;               // move the rear pointer forward.
                M[row][col] = 0;      // mark the cell as visited.
            }
        }
    }
}

// this function counts the number of connected regions (islands) in the binary matrix.
int countRegions(int **regions, int numRow, int numCol)
{
    int numRegions = 0; // starting the count for regions.
    // loop through each cell in the matrix to find unvisited regions.
    for (int r = 0; r < numRow; r++)
    {
        for (int c = 0; c < numCol; c++)
        {
            if (regions[r][c] == 1)
            {                                              // if the cell is part of an unvisited region.
                findRegion(regions, r, c, numRow, numCol); // call BFS to mark all connected cells.
                numRegions++;                              // increment the count of regions found.
            }
        }
    }
    return numRegions; // return the total number of connected regions.
}

int read_data(const char *file_path, Data* data, int rows, int cols) {
    FILE *file = fopen(file_path, "r");
    if (!file) {
        perror("Unable to open file");
        return -1;
    }

    char line[256];
    int n = 0; // Counter for filling the grid

    // Skip the header line
    fgets(line, sizeof(line), file);

    // Process each line in the CSV
    while (fgets(line, sizeof(line), file)) {
        double rain_val, temp_val, humidity_val;
        // Parse rain, temp, humidity, and aqi values from CSV line
        if (sscanf(line, "%*lf,%*lf,%lf,%lf,%lf", &temp_val, &humidity_val, &rain_val) == 3) {
            int row = n / cols;
            int col = n % cols;

            // Check if indices are within the matrix bounds
            if (row < rows && col < cols) {
                // Assign parsed values to the respective matrices in the struct
                data->rain[row][col] = rain_val;
                data->temp[row][col] = temp_val;
                data->humidity[row][col] = humidity_val;
                n++; // Increment the counter
            } else {
                fprintf(stderr, "Index out of bounds for n: %d\n", n);
                break;
            }
        } else {
            fprintf(stderr, "Failed to parse line: %s\n", line);
        }
    }

    fclose(file);
    return 0;
}

void updateMarked(int **marked, Data *data, Threshold *thresholds, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if ((data->temp[i][j] > thresholds[0].min && data->temp[i][j] < thresholds[0].max) &&
                (data->humidity[i][j] > thresholds[1].min && data->humidity[i][j] < thresholds[1].max) &&
                (data->rain[i][j] > thresholds[2].min && data->rain[i][j] < thresholds[2].max)) {
                marked[i][j] = 1;
            } else {
                marked[i][j] = 0;
            }
        }
    }
}

void processQuery(FILE *outputFile, Data *data, SegmentTree2D *segTree, char *type, int** marked , Threshold *thresholds) {
    clock_t start, end;
    Query query;
    
    query.valid = 1;  // Assume the query is valid initially
    query.result = 0.0;
    double lat_1 , long_1 , lat_2 , long_2;
    if (strcasecmp(type, "Range") == 0) {
        printf("Enter range query type (Sum, Avg, Max or Min): ");
        scanf("%s", query.rangeType);
        printf("Enter coordinates (x1, y1, x2, y2): ");
        scanf("%lf %lf %lf %lf", &lat_1, &long_1, &lat_2, &long_2);
        query.x1 = (int)round((lat_1 - 9.00) / 0.20);
        query.x2 = (int)round((lat_2 - 9.00) / 0.20);
        query.y1 = (int)round((long_1 - 69.00) / 0.20);
        query.y2 = (int)round((long_2 - 69.00) / 0.20);
        printf("Enter Temperature-0, Humidity-1 or Rain-2: ");
        scanf("%d", &query.attribute);

        start = clock();
        if (strcasecmp(query.rangeType, "Avg") == 0) {
            // Call avg query function with data and other parameters
            query.result = avgQuery(data, segTree, 1, 0, rows, query.x1, query.x2, query.y1, query.y2, query.attribute);
            
        } 
        else if (strcasecmp(query.rangeType, "Sum") == 0) {
            // Call sum query function with data and other parameters
            query.result = sumQuery(data, segTree, 1, 0, rows, query.x1, query.x2, query.y1, query.y2, query.attribute);
            
        } 
        else if (strcasecmp(query.rangeType, "Max") == 0) {
            // Call max query function with data and other parameters
            query.result = maxQuery( segTree, 1, 0, rows, query.x1, query.x2, query.y1, query.y2, query.attribute);
            
        } 
        else if (strcasecmp(query.rangeType, "Min") == 0) {
            // Call min query function with data and other parameters
            query.result = minQuery( segTree, 1, 0, rows, query.x1, query.x2, query.y1, query.y2, query.attribute);
            
        } 
        else {
            fprintf(outputFile, "Query failed: Unsupported range query type '%s'\n", query.rangeType);
            query.valid = 0;
            return;
        }
        end = clock();
    } 
    else if (strcasecmp(type, "Update") == 0) {
        printf("Enter coordinates to update (row, column): ");
        scanf("%lf %lf", &lat_1, &long_1);
        query.row = (int)round((lat_1 - 9.0) / 0.2);
        query.col = (int)round((long_1 - 69.0) / 0.2);

        printf("Enter Temperature-0, Humidity-1 or Rain-2: ");
        scanf("%d", &query.attribute);

        printf("Enter value to update at (%.2f, %.2f): ", lat_1, long_1);
        scanf("%lf", &query.value);

        start = clock();
        // Call update function with row, column, and value
        update(data,segTree, 1, 0, rows, query.row, query.col, query.value, thresholds, marked , query.attribute);
        end = clock();
        printf("Updated successfully!");
        fprintf(outputFile, "Update performed at coordinates (%.2f, %.2f) with value %.2f\n", lat_1, long_1, query.value);
    }
    else if (strcasecmp(type, "Combined") == 0){
        printf("Enter min and max values for temperature: ");
        scanf("%lf %lf", &thresholds[0].min, &thresholds[0].max);
        printf("Enter min and max values for humidity: ");
        scanf("%lf %lf", &thresholds[1].min, &thresholds[1].max);
        printf("Enter min and max values for rain: ");
        scanf("%lf %lf", &thresholds[2].min, &thresholds[2].max);

        updateMarked(marked , data , thresholds , rows , cols);
        int topLeftX, topLeftY, bottomRightX, bottomRightY;
        int maxArea = maximalRegion(marked, rows, cols, &topLeftX, &topLeftY, &bottomRightX, &bottomRightY);
        int numRegions = countRegions(marked , rows , cols);
        double tl_lat = 9.0+(topLeftX * 0.2) , br_lat = 9.0+(bottomRightX * 0.2);
        double tl_long = 69.0+(topLeftY * 0.2) , br_long = 69.0+(bottomRightY * 0.2);
        fprintf(outputFile ,"Temperature thresholds: min = %.2f, max = %.2f\n", thresholds[0].min, thresholds[0].max);
        fprintf(outputFile ,"Humidity thresholds: min = %.2f, max = %.2f\n", thresholds[1].min, thresholds[1].max);
        fprintf(outputFile ,"Rain thresholds: min = %.2f, max = %.2f\n", thresholds[2].min, thresholds[2].max);
        fprintf(outputFile ,"Maximal Rectangle Area: %d\n", maxArea);
        fprintf(outputFile ,"Top-left corner: (%.2f, %.2f)\n", tl_lat, tl_long);
        fprintf(outputFile ,"Bottom-right corner: (%.2f, %.2f)\n", br_lat, br_long); 
        fprintf(outputFile ,"No. of Regions are: %d\n" , numRegions);

    }
    else {
        fprintf(outputFile, "Query failed: Unsupported query type '%s'\n", type);
        query.valid = 0;
        return;
    }

    query.timeTaken = ((double)(end - start)) / CLOCKS_PER_SEC;

    if (query.valid && strcasecmp(type, "Range") == 0) {
        
        fprintf(outputFile, "Query Type: %s\n", type);
        fprintf(outputFile, "Range Query Type: %s\n", query.rangeType);
        fprintf(outputFile, "Coordinates: (%.2f, %.2f) to (%.2f, %.2f)\n", lat_1, long_1, lat_2, long_2);
        fprintf(outputFile, "Result: %.2f\n", query.result);
        fprintf(outputFile, "Time Taken: %.3f seconds\n\n", query.timeTaken);
        if (query.attribute==0){
            fprintf(outputFile, "Attribute: Temperature\n" );
        }
        else if (query.attribute==1){
            fprintf(outputFile, "Attribute: Humidity\n" );
        }
        else if (query.attribute==2){
            fprintf(outputFile, "Attribute: Rain\n" );
        }
        
    } 
    else if (strcasecmp(type, "Update") == 0) {
        fprintf(outputFile, "Time Taken for Update: %.3f seconds\n\n", query.timeTaken);
    }
}

int main()
{
    // input sample grid size.
    rows = MAX_ROWS;
    cols = MAX_COLS;

    //initializing required arrays.
    SegmentTree2D *segTree = create2DSegmentTree(rows, cols);
    Data *data = createDataMatrix(rows , cols);
    int** markedArray = (int**)malloc(rows * sizeof(int*));
    for (int i = 0; i < rows; i++) {
        markedArray[i] = (int*)malloc(cols * sizeof(int));
        for (int j = 0; j < cols; j++) {
            markedArray[i][j] = 0;
        }
    }

    //defining default thresholds
    Threshold* thresholds = malloc(3 * sizeof(Threshold));
    thresholds[0] = (Threshold){0.0, 100.0};  // Temperature threshold
    thresholds[1] = (Threshold){0.0, 100.0};   // Rain threshold
    thresholds[2] = (Threshold){0.0, 100.0};  

    // Read data from the file and populate the Data struct
    if (read_data("environmental_data_2022.csv", data, rows, cols) != 0) {
        fprintf(stderr, "Error reading data from file.\n");
        return EXIT_FAILURE;
    }

    build(data, segTree, 1, 0, rows - 1 , cols , thresholds , markedArray); //building the segment tree.
    //taking the output.
    FILE *outputFile = fopen("output_summary.txt", "w");
    if (outputFile == NULL) {
        fprintf(stderr, "Error opening file for writing\n");
        return 1;
    }
    time_t now;
    time(&now);
    fprintf(outputFile, "Output Metadata:\n");
    fprintf(outputFile, "Dataset: environmental_data_2022.csv\n");
    fprintf(outputFile, "Generated on: %s\n\n", ctime(&now));
    int processedQueries = 0;
    char type[256];
    while (1) {
        printf("\nEnter query type (Range, Update, Combined): or type 'exit' to finish: ");
        scanf("%s", type);

        if (strcasecmp(type, "exit") == 0) {
            break;
        }
        else{
            processQuery(outputFile , data,segTree , type, markedArray , thresholds);
            processedQueries++;
        }
    }
    fprintf(outputFile, "Total Queries Processed: %d\n", processedQueries);
    fclose(outputFile);


    //Free allocated memory.
    free_data(data, rows);
    free2DSegmentTree(segTree);
    for (int i = 0; i < rows; i++) {
        free(markedArray[i]);
    }
    free(markedArray);

    return 0;
}
