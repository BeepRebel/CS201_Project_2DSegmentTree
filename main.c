#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>

#define MAX_ROWS 1000 // defined the maximum no. of rows.
#define MAX_COLS 1000 // defined the maximum no. of columns.
#define GRID_SIZE 100

int rows, cols;                           // dimensions of the grid.


int marked[MAX_ROWS][MAX_COLS]; // Matrix to track cells above the threshold
int threshold;                  // Define a threshold for marking cells

typedef struct {
    int rows, cols;
    double ***minMatrix; // 3D matrix to hold min values for each parameter
    double ***maxMatrix; // 3D matrix to hold max values for each parameter
    double ***sumMatrix; // 3D matrix to hold sum values for each parameter
} SegmentTree2D;

typedef struct {
    int rows , cols;
    double **temp;
    double **humidity;
    double **rain;
} Data;



// Helper function to allocate 2D Segment Tree matrices
SegmentTree2D *create2DSegmentTree(int rows, int cols) {
    SegmentTree2D *tree = malloc(sizeof(SegmentTree2D));
    tree->rows = rows;
    tree->cols = cols;

    // Allocate memory for each parameter (temperature, humidity, rain)
    tree->minMatrix = (double***)malloc(4 * rows * sizeof(double**));
    tree->maxMatrix = (double***)malloc(4 * rows * sizeof(double**));
    tree->sumMatrix = (double***)malloc(4 * rows * sizeof(double**));

    for (int k = 0; k < 4 * rows; ++k) {
        tree->minMatrix[k] = (double**)malloc(4 * cols * sizeof(double*));
        tree->maxMatrix[k] = (double**)malloc(4 * cols * sizeof(double*));
        tree->sumMatrix[k] = (double**)malloc(4 * cols * sizeof(double*));

        for (int i = 0; i < 4*cols; ++i) {
            tree->minMatrix[k][i] = (double*)malloc(3 * sizeof(double));
            tree->maxMatrix[k][i] = (double*)malloc(3 * sizeof(double));
            tree->sumMatrix[k][i] = (double*)malloc(3 * sizeof(double));
            for(int j = 0 ; j < 3 ; j++){
            tree->minMatrix[k][i][j] = DBL_MAX;
            tree->maxMatrix[k][i][j] = -DBL_MAX;
            tree->sumMatrix[k][i][j] = 0.0;
            }
        }
    }

    return tree;
}

Data *createDataMatrix(int rows, int cols)
{
    Data *data = (Data *)malloc(sizeof(Data));
    data->rows = rows;
    data->cols = cols;

    // Allocate 4x larger matrices for segment tree data
    data->temp = (double **)malloc(rows * sizeof(double *));
    data->humidity = (double **)malloc(rows * sizeof(double *));
    data->rain = (double **)malloc(rows * sizeof(double *));

    for (int i = 0; i < rows; i++)
    {
        data->temp[i] = (double *)malloc(cols * sizeof(double));
        data->humidity[i] = (double *)malloc(cols * sizeof(double));
        data->rain[i] = (double *)malloc(cols * sizeof(double));
        for (int j = 0; j < cols; j++)
        {
            data->temp[i][j] = 0.0;
            data->humidity[i][j] = 0.0;
            data->rain[i][j] = 0.0;
        }
    }
    return data;
}

// Helper function to free the Segment Tree matrices
void freeSegmentTree(SegmentTree2D *tree, int parameters) {
    // Free each parameter's matrices (minMatrix, maxMatrix, sumMatrix)
    for (int k = 0; k < parameters; ++k) {
        for (int i = 0; i < tree->rows; ++i) {
            free(tree->minMatrix[k][i]);
            free(tree->maxMatrix[k][i]);
            free(tree->sumMatrix[k][i]);
        }
        free(tree->minMatrix[k]);
        free(tree->maxMatrix[k]);
        free(tree->sumMatrix[k]);
    }

    // Free the outermost arrays
    free(tree->minMatrix);
    free(tree->maxMatrix);
    free(tree->sumMatrix);

    // Free the SegmentTree2D structure itself
    free(tree);
}

// Helper function to free the Segment Tree matrices
void freeData(Data *data)
{
    for (int i = 0; i < data->rows; i++)
    {
        free(data->temp[i]);
        free(data->humidity[i]);
        free(data->rain[i]);
    }
        free(data->temp);
        free(data->humidity);
        free(data->rain);
    free(data);
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

void buildRow(Data* data, SegmentTree2D *segTree, int row, int lx, int rx, int col, int ly, int ry, int threshold) {
    // Case for a single column.
    if (ly == ry) {
        // Case for a single row.
        if (lx == rx) {
            // Initialize sum and max for each parameter with the cell value.
            segTree->sumMatrix[row][col][0] = data->temp[lx][ly];
            segTree->maxMatrix[row][col][0] = data->temp[lx][ly];
            segTree->sumMatrix[row][col][1] = data->humidity[lx][ly];
            segTree->maxMatrix[row][col][1] = data->humidity[lx][ly];
            segTree->sumMatrix[row][col][2] = data->rain[lx][ly];
            segTree->maxMatrix[row][col][2] = data->rain[lx][ly];

            // Set marked array based on the temperature threshold.
            marked[lx][ly] = (data->temp[lx][ly] > threshold) ? 1 : 0;
        }
        // Case when not a single row.
        else {
            // Merge the max and sum values for the row's children for each parameter.
            for (int i = 0; i < 3; i++) {
                segTree->sumMatrix[row][col][i] = segTree->sumMatrix[2 * row][col][i] + segTree->sumMatrix[2 * row + 1][col][i];
                segTree->maxMatrix[row][col][i] = fmax(segTree->maxMatrix[2 * row][col][i], segTree->maxMatrix[2 * row + 1][col][i]);
            }
        }
    }
    // Case when not a single column.
    else {
        // Recursively build for each half.
        int my = (ly + ry) / 2; // Middle column index for splitting.
        buildRow(data, segTree, row, lx, rx, 2 * col, ly, my, threshold);        // Build left child.
        buildRow(data, segTree, row, lx, rx, 2 * col + 1, my + 1, ry, threshold); // Build right child.

        // Merge values from left and right children for each parameter.
        for (int i = 0; i < 3; i++) {
            segTree->sumMatrix[row][col][i] = segTree->sumMatrix[row][2 * col][i] + segTree->sumMatrix[row][2 * col + 1][i];
            segTree->maxMatrix[row][col][i] = fmax(segTree->maxMatrix[row][2 * col][i], segTree->maxMatrix[row][2 * col + 1][i]);
        }
    }
}

// row: index in the segment tree for rows
// lx, rx: range of rows being covered by the current node

// Function to build the entire 2D segment tree.
void build(Data* data, SegmentTree2D *segTree, int row, int lx, int rx, int cols, int threshold) {
    // Case when not a single row.
    if (lx != rx) {
        int mid = (lx + rx) / 2;                        // Calculate the middle row index.
        build(data, segTree, 2 * row, lx, mid, cols, threshold);         // Build left child for rows.
        build(data, segTree, 2 * row + 1, mid + 1, rx, cols, threshold); // Build right child for rows.
    }

    // Build the column segment for this row segment.
    buildRow(data, segTree, row, lx, rx, 1, 0, cols - 1, threshold);
}


// function to query the sum of a rectangular sub-region for a specific row segment
// parameters:
// - row: the index representing the current row segment in the segment tree
// - col: the index representing the current column segment in the segment tree
// - t_ly, t_ry: the range of columns that the current node covers
// - ly, ry: the range of columns we want to query for the sum
double sumQueryRow(Data* data, SegmentTree2D *segTree, int row, int col, int t_ly, int t_ry, int ly, int ry , int attribute)
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
    return sumQueryRow(data, segTree, row, 2 * col, t_ly, tmy, ly, min(ry, tmy) , attribute) +            // left child query
           sumQueryRow(data, segTree, row, 2 * col + 1, tmy + 1, t_ry, max(ly, tmy + 1), ry , attribute); // right child query
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
        return sumQueryRow(data, segTree, row, 1, 0, cols - 1, ly, ry , attribute);
    }

    // find the middle row of the current segment
    int t_mid = (t_lx + t_rx) / 2;

    // recursively query the left child (first half of the row range)
    // combine it with the result of querying the right child (second half of the row range)
    return sumQuery(data, segTree, 2 * row, t_lx, t_mid, lx, min(rx, t_mid), ly, ry , attribute) +            // left child query
           sumQuery(data, segTree, 2 * row + 1, t_mid + 1, t_rx, max(lx, t_mid + 1), rx, ly, ry , attribute); // right child query
}

// function to query the maximum value of a rectangular sub-region for a specific row segment
// parameters:
// - row: the index representing the current row segment in the segment tree
// - col: the index representing the current column segment in the segment tree
// - t_ly, t_ry: the range of columns that the current node covers
// - ly, ry: the range of columns we want to query for the maximum value
double maxQueryRow(Data* data, SegmentTree2D *segTree, int row, int col, int t_ly, int t_ry, int ly, int ry , int attribute)
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
    return fmax(maxQueryRow(data, segTree, row, 2 * col, t_ly, t_mid, ly, min(ry, t_mid) , attribute),              // left child query
                maxQueryRow(data, segTree, row, 2 * col + 1, t_mid + 1, t_ry, max(ly, t_mid + 1), ry , attribute)); // right child query
}

// function to query the maximum value of a rectangular sub-region in the 2d segment tree
// parameters:
// - row: the index representing the current row segment in the segment tree
// - t_lx, t_rx: the range of rows that the current node covers
// - lx, rx: the range of rows we want to query for the maximum value
// - ly, ry: the range of columns we want to query for the maximum value
double maxQuery(Data* data, SegmentTree2D *segTree, int row, int t_lx, int t_rx, int lx, int rx, int ly, int ry , int attribute)
{
    // if the queried range is invalid (e.g., start is greater than end), return INT_MIN as there's no valid value
    if (lx > rx)
        return DBL_MIN;

    // if the current segment matches exactly with the row query range,
    // delegate the column query to maxQueryRow function
    if (lx == t_lx && rx == t_rx)
    {
        return maxQueryRow(data, segTree, row, 1, 0, cols - 1, ly, ry , attribute);
    }

    // find the middle row of the current segment
    int t_mid = (t_lx + t_rx) / 2;

    // recursively query the left child (first half of the row range)
    // combine it with the result of querying the right child (second half of the row range) to find the maximum
    return fmax(maxQuery(data, segTree, 2 * row, t_lx, t_mid, lx, min(rx, t_mid), ly, ry , attribute),              // left child query
                maxQuery(data, segTree, 2 * row + 1, t_mid + 1, t_rx, max(lx, t_mid + 1), rx, ly, ry , attribute)); // right child query
}

// function to update a single point in the segment tree for a specific row segment
// parameters:
// - row: the index representing the current row segment in the segment tree
// - col: the index representing the current column segment in the segment tree
// - lx, rx: the range of rows being covered by the current node
// - ly, ry: the range of columns being covered by the current node
// - a, b: the coordinates (row, column) of the point we want to update
// - value: the new value to update the point to
void updateRow(Data* data, SegmentTree2D *segTree, int row, int lx, int rx, int col, int ly, int ry, int a, int b, int value)
{
    // if we're at a single column (a leaf node in the column segment tree)
    if (ly == ry)
    {
        // check if we're at a single row (a leaf node in the row segment tree)
        if (lx == rx)
        {
            // if both the row and column are leaf nodes, update the point directly
            for(int i = 0 ; i < 3 ; i++){
                segTree->sumMatrix[row][col][i] = value; // update the sum value at the point
                segTree->maxMatrix[row][col][i] = value;           // update the max value at the point
            }
            marked[lx][ly] = (value > threshold) ? 1 : 0;
        }
        else
        {
            // if it's not a single row, merge the updates from both child row nodes
            for(int i = 0 ; i < 3 ; i++){
                segTree->sumMatrix[row][col][i] = segTree->sumMatrix[2 * row][col][i] + segTree->sumMatrix[2 * row + 1][col][i]; // update the sum
                segTree->maxMatrix[row][col][i] = fmax(segTree->maxMatrix[2 * row][col][i], segTree->maxMatrix[2 * row + 1][col][i]);                          // update the max
            }
        }
    }
    else
    {
        // if it's not a single column, split the columns and update the relevant child
        int mid = (ly + ry) / 2;

        // if the column index is within the left half, update the left child
        if (b <= mid)
        {
            updateRow(data, segTree, row, lx, rx, 2 * col, ly, mid, a, b, value);
        }
        else
        {
            // otherwise, update the right child
            updateRow(data, segTree, row, lx, rx, 2 * col + 1, mid + 1, ry, a, b, value);
        }

        // after updating the children, merge their sums and max values
        for(int i = 0 ; i < 3 ; i++){
            segTree->sumMatrix[row][col][i] = segTree->sumMatrix[row][2 * col][i] + segTree->sumMatrix[row][2 * col + 1][i]; // merge sums
            segTree->maxMatrix[row][col][i] = fmax(segTree->maxMatrix[row][2 * col][i], segTree->maxMatrix[row][2 * col + 1][i]);                          // merge maxes
        }
    }
}

// function to update a single point in the 2d segment tree
// parameters:
// - row: the index representing the current row segment in the segment tree
// - lx, rx: the range of rows being covered by the current node
// - a, b: the coordinates (row, column) of the point we want to update
// - value: the new value to update the point to
void update(Data* data, SegmentTree2D *segTree, int row, int lx, int rx, int a, int b, int value)
{
    // if it's not a single row (non-leaf node in the row segment tree)
    if (lx != rx)
    {
        int mid = (lx + rx) / 2;

        // if the row index is within the left half, update the left child
        if (a <= mid)
        {
            update(data, segTree, 2 * row, lx, mid, a, b, value);
        }
        else
        {
            // otherwise, update the right child
            update(data, segTree, 2 * row + 1, mid + 1, rx, a, b, value);
        }
    }

    // after updating the rows, update the column segment tree for this row segment
    updateRow(data, segTree, row, lx, rx, 1, 0, cols - 1, a, b, value);
}

int *nearestSmallerToLeft(int *heights, int size, int pseudo)
{
    int *v = (int *)malloc(size * sizeof(int));
    int stack[size][2]; // Stack to store elements and their indices
    int top = -1;

    for (int i = 0; i < size; i++)
    {
        // Pop elements from stack that are greater than or equal to the current element
        while (top >= 0 && stack[top][0] >= heights[i])
        {
            top--;
        }
        // Store the index of the nearest smaller element to the left, or -1 if none
        v[i] = (top == -1) ? pseudo : stack[top][1];

        // Push current element and index to the stack
        stack[++top][0] = heights[i];
        stack[top][1] = i;
    }
    return v;
}

// Function to find the nearest smaller elements to the right for each element in `heights`
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

// Function to calculate the maximum area in a histogram-like row
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
// Function to free the allocated memory in the Data struct
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
        double rain_val, temp_val, humidity_val, aqi_val;
        // Parse rain, temp, humidity, and aqi values from CSV line
        if (sscanf(line, "%lf,%*lf,%lf,%*lf,%lf,%lf,%*lf,%*lf", &rain_val, &temp_val, &humidity_val, &aqi_val) == 4) {
            int row = n / GRID_SIZE;
            int col = n % GRID_SIZE;

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

int main()
{
    // input sample grid size.
    rows = GRID_SIZE;
    cols = GRID_SIZE;

    SegmentTree2D *segTree = create2DSegmentTree(rows, cols);
    Data *data = createDataMatrix(rows , cols);

    // sample grid
    // Read data from the file and populate the Data struct
    if (read_data("C:/Users/Krishna/OneDrive/Desktop/CS201 Project/combined_data_2019.csv", data, rows, cols) != 0) {
        fprintf(stderr, "Error reading data from file.\n");
        return EXIT_FAILURE;
    }

    for(int i = 0 ; i < 10 ; i++){
        for(int j = 0 ; j < 10 ; j++){
            printf("%.2f " , data->temp[i][j]);
        }
        printf("\n");
    }

    build(data, segTree, 1, 0, rows - 1 , cols , 10.0);

    printf("\n");
    for(int i = 0 ; i < 10 ; i++){
        for(int j = 0 ; j < 10 ; j++){
            printf("%.2f " , segTree->sumMatrix[i][j][0]);
        }
        printf("\n");
    }
    printf("\nSum of subgrid (1, 1) to (3, 3): %.2f\n", maxQuery(data, segTree, 1, 0, rows - 1, 0, 80, 1, 80 ,0));

    threshold = 10;
    free_data(data, rows);
    freeSegmentTree(segTree , 3);
    // Copy example grid to global grid
    // for (int i = 0; i < n; i++)
    // {
    //     for (int j = 0; j < m; j++)
    //     {
    //         grid[i][j] = temp[i][j];
    //     }
    // }

    // //building the segment tree for the input grid.
    // build(1, 0, n - 1);

    // // performing sum query
    // printf("Sum of subgrid (1, 1) to (3, 3): %d\n", sumQuery(1, 0, n - 1, 1, 3, 1, 3));

    // // performing max query
    // printf("Max value in subgrid (1, 1) to (3, 3): %d\n", maxQuery(1, 0, n - 1, 1, 3, 1, 3));

    // // update a value and query again
    // update(1, 0, n - 1, 2, 2, 20); // update grid[2][2] to 20
    // // sum
    // printf("Sum of subgrid (1, 1) to (3, 3) after update: %d\n", sumQuery(1, 0, n - 1, 1, 3, 1, 3));
    // // max
    // printf("Max value in subgrid (1, 1) to (3, 3) after update: %d\n", maxQuery(1, 0, n - 1, 1, 3, 1, 3));
    // int** markedArray = (int**)malloc(n * sizeof(int*));
    // for (int i = 0; i < n; i++) {
    //     markedArray[i] = (int*)malloc(m * sizeof(int));
    //     for (int j = 0; j < m; j++) {
    //         markedArray[i][j] = marked[i][j];
    //     }
    // }
    //     printf("markedArray:\n");
    // for (int i = 0; i < n; i++) {
    //     for (int j = 0; j < m; j++) {
    //         printf("%d ", markedArray[i][j]);
    //     }
    //     printf("\n");
    // }

    // int topLeftX, topLeftY, bottomRightX, bottomRightY;
    // int maxArea = maximalRegion(markedArray, n, m, &topLeftX, &topLeftY, &bottomRightX, &bottomRightY);

    // printf("Maximal Rectangle Area: %d\n", maxArea);
    // printf("Top-left corner: (%d, %d)\n", topLeftX, topLeftY);
    // printf("Bottom-right corner: (%d, %d)\n", bottomRightX, bottomRightY);

    // int numRegions = countRegions(markedArray , n , m);
    // printf("%d\n" , numRegions);

    // // Free allocated memory for markedArray
    // for (int i = 0; i < n; i++) {
    //     free(markedArray[i]);
    // }
    // free(markedArray);

    return 0;
}
