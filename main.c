#include <stdio.h>
#include <limits.h>

#define MAX_ROWS 1000 //defined the maximum no. of rows.
#define MAX_COLS 1000 // defined the maximum no. of columns.

int n, m;  //dimensions of the grid.
int RangeSum[4 * MAX_ROWS][4 * MAX_COLS]; // for sum queries.
int RangeMax[4 * MAX_ROWS][4 * MAX_COLS]; // for max queries.
int grid[MAX_ROWS][MAX_COLS];  //input grid defined globally.

int marked[MAX_ROWS][MAX_COLS]; // Matrix to track cells above the threshold
int threshold; // Define a threshold for marking cells


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

void buildRow(int row, int lx, int rx, int col, int ly, int ry)
{
    // case for a single column.
    if (ly == ry)
    {
        // case for single row.
        if (lx == rx)
        {
            // initializing the sum and column with the cell value.
            RangeSum[row][col] = grid[lx][ly];
            RangeMax[row][col] = grid[lx][ly];
            marked[lx][ly] = (grid[lx][ly] > threshold) ? 1 : 0;
        }
        //case when not a single row.
        else
        {
            // merge the max and sum value for the row's children.
            RangeSum[row][col] = RangeSum[2 * row][col] + RangeSum[2 * row + 1][col];
            RangeMax[row][col] = max(RangeMax[2 * row][col], RangeMax[2 * row + 1][col]);
        }
    }
    // case when not a single column.
    else
    {
        // recursively building for each half.
        int my = (ly + ry) / 2; // middle column index for splitting.
        buildRow(row, lx, rx, 2 * col, ly, my);           // building the left child.
        buildRow(row, lx, rx, 2 * col + 1, my + 1, ry);   // building the right child.
        // merging the left and right child values.
        RangeSum[row][col] = RangeSum[row][2 * col] + RangeSum[row][2 * col + 1];
        RangeMax[row][col] = max(RangeMax[row][2 * col], RangeMax[row][2 * col + 1]);
    }
}

// row: index in the segment tree for rows
// lx, rx: range of rows being covered by the current node

// function to build the entire 2D segment tree for answering queries.

void build(int row, int lx, int rx)
{
    // case when not a single row.
    if (lx != rx)
    {
        int mid = (lx + rx) / 2;       // calculating the middle row index.
        build(2 * row, lx, mid);       // building the left child for the row segment.
        build(2 * row + 1, mid + 1, rx); // building the right child for the row segment.
    }

    // building the column segment for this row segment.
    buildRow(row, lx, rx, 1, 0, m - 1);
}

// function to query the sum of a rectangular sub-region for a specific row segment
// parameters:
// - row: the index representing the current row segment in the segment tree
// - col: the index representing the current column segment in the segment tree
// - t_ly, t_ry: the range of columns that the current node covers
// - ly, ry: the range of columns we want to query for the sum
int sumQueryRow(int row, int col, int t_ly, int t_ry, int ly, int ry)
{
    // if the queried range is invalid (e.g., start is greater than end), return 0 as there's no sum
    if (ly > ry)
        return 0;

    // if the current segment matches exactly with the query range, return the sum stored in this node
    if (ly == t_ly && ry == t_ry)
    {
        return RangeSum[row][col]; 
    }

    // find the middle column of the current segment
    int tmy = (t_ly + t_ry) / 2;

    // recursively query the left child (first half of the column range)
    // combine it with the result of querying the right child (second half of the column range)
    return sumQueryRow(row, 2 * col, t_ly, tmy, ly, min(ry, tmy)) +             // left child query
           sumQueryRow(row, 2 * col + 1, tmy + 1, t_ry, max(ly, tmy + 1), ry);  // right child query
}

// function to query the sum of a rectangular sub-region in the 2D segment tree
// parameters:
// - row: the index representing the current row segment in the segment tree
// - t_lx, t_rx: the range of rows that the current node covers
// - lx, rx: the range of rows we want to query for the sum
// - ly, ry: the range of columns we want to query for the sum
int sumQuery(int row, int t_lx, int t_rx, int lx, int rx, int ly, int ry)
{
    // if the queried range is invalid (e.g., start is greater than end), return 0 as there's no sum
    if (lx > rx)
        return 0;

    // if the current segment matches exactly with the row query range,
    // delegate the column query to sumQueryRow function
    if (lx == t_lx && rx == t_rx)
    {
        return sumQueryRow(row, 1, 0, m - 1, ly, ry); 
    }

    // find the middle row of the current segment
    int t_mid = (t_lx + t_rx) / 2;

    // recursively query the left child (first half of the row range)
    // combine it with the result of querying the right child (second half of the row range)
    return sumQuery(2 * row, t_lx, t_mid, lx, min(rx, t_mid), ly, ry) +             // left child query
           sumQuery(2 * row + 1, t_mid + 1, t_rx, max(lx, t_mid + 1), rx, ly, ry);  // right child query
}

// function to query the maximum value of a rectangular sub-region for a specific row segment
// parameters:
// - row: the index representing the current row segment in the segment tree
// - col: the index representing the current column segment in the segment tree
// - t_ly, t_ry: the range of columns that the current node covers
// - ly, ry: the range of columns we want to query for the maximum value
int maxQueryRow(int row, int col, int t_ly, int t_ry, int ly, int ry)
{
    // if the queried range is invalid (e.g., start is greater than end), return INT_MIN as there's no valid value
    if (ly > ry)
        return INT_MIN;

    // if the current segment matches exactly with the query range, return the maximum value stored in this node
    if (ly == t_ly && ry == t_ry)
    {
        return RangeMax[row][col]; 
    }

    // find the middle column of the current segment
    int t_mid = (t_ly + t_ry) / 2;

    // recursively query the left child (first half of the column range)
    // combine it with the result of querying the right child (second half of the column range) to find the maximum
    return max(maxQueryRow(row, 2 * col, t_ly, t_mid, ly, min(ry, t_mid)),               // left child query
               maxQueryRow(row, 2 * col + 1, t_mid + 1, t_ry, max(ly, t_mid + 1), ry));  // right child query
}

// function to query the maximum value of a rectangular sub-region in the 2d segment tree
// parameters:
// - row: the index representing the current row segment in the segment tree
// - t_lx, t_rx: the range of rows that the current node covers
// - lx, rx: the range of rows we want to query for the maximum value
// - ly, ry: the range of columns we want to query for the maximum value
int maxQuery(int row, int t_lx, int t_rx, int lx, int rx, int ly, int ry)
{
    // if the queried range is invalid (e.g., start is greater than end), return INT_MIN as there's no valid value
    if (lx > rx)
        return INT_MIN;

    // if the current segment matches exactly with the row query range,
    // delegate the column query to maxQueryRow function
    if (lx == t_lx && rx == t_rx)
    {
        return maxQueryRow(row, 1, 0, m - 1, ly, ry); 
    }

    // find the middle row of the current segment
    int t_mid = (t_lx + t_rx) / 2;

    // recursively query the left child (first half of the row range)
    // combine it with the result of querying the right child (second half of the row range) to find the maximum
    return max(maxQuery(2 * row, t_lx, t_mid, lx, min(rx, t_mid), ly, ry),               // left child query
               maxQuery(2 * row + 1, t_mid + 1, t_rx, max(lx, t_mid + 1), rx, ly, ry));  // right child query
}

// function to update a single point in the segment tree for a specific row segment
// parameters:
// - row: the index representing the current row segment in the segment tree
// - col: the index representing the current column segment in the segment tree
// - lx, rx: the range of rows being covered by the current node
// - ly, ry: the range of columns being covered by the current node
// - a, b: the coordinates (row, column) of the point we want to update
// - value: the new value to update the point to
void updateRow(int row, int lx, int rx, int col, int ly, int ry, int a, int b, int value)
{
    // if we're at a single column (a leaf node in the column segment tree)
    if (ly == ry)
    { 
        // check if we're at a single row (a leaf node in the row segment tree)
        if (lx == rx)
        {                               
            // if both the row and column are leaf nodes, update the point directly
            RangeSum[row][col] = value; // update the sum value at the point
            RangeMax[row][col] = value; // update the max value at the point
            marked[lx][ly] = (value > threshold) ? 1 : 0;
        }
        else
        { 
            // if it's not a single row, merge the updates from both child row nodes
            RangeSum[row][col] = RangeSum[2 * row][col] + RangeSum[2 * row + 1][col]; // update the sum
            RangeMax[row][col] = max(RangeMax[2 * row][col], RangeMax[2 * row + 1][col]); // update the max
        }
    }
    else
    { 
        // if it's not a single column, split the columns and update the relevant child
        int mid = (ly + ry) / 2;
        
        // if the column index is within the left half, update the left child
        if (b <= mid)
        {
            updateRow(row, lx, rx, 2 * col, ly, mid, a, b, value);
        }
        else
        {
            // otherwise, update the right child
            updateRow(row, lx, rx, 2 * col + 1, mid + 1, ry, a, b, value);
        }

        // after updating the children, merge their sums and max values
        RangeSum[row][col] = RangeSum[row][2 * col] + RangeSum[row][2 * col + 1]; // merge sums
        RangeMax[row][col] = max(RangeMax[row][2 * col], RangeMax[row][2 * col + 1]); // merge maxes
    }
}

// function to update a single point in the 2d segment tree
// parameters:
// - row: the index representing the current row segment in the segment tree
// - lx, rx: the range of rows being covered by the current node
// - a, b: the coordinates (row, column) of the point we want to update
// - value: the new value to update the point to
void update(int row, int lx, int rx, int a, int b, int value)
{
    // if it's not a single row (non-leaf node in the row segment tree)
    if (lx != rx)
    {
        int mid = (lx + rx) / 2;
        
        // if the row index is within the left half, update the left child
        if (a <= mid)
        {
            update(2 * row, lx, mid, a, b, value); 
        }
        else
        {
            // otherwise, update the right child
            update(2 * row + 1, mid + 1, rx, a, b, value);
        }
    }
    
    // after updating the rows, update the column segment tree for this row segment
    updateRow(row, lx, rx, 1, 0, m - 1, a, b, value); 
}

int* nearestSmallerToLeft(int* heights, int size, int pseudo) {
    int* v = (int*)malloc(size * sizeof(int));
    int stack[size][2];  // Stack to store elements and their indices
    int top = -1;

    for (int i = 0; i < size; i++) {
        // Pop elements from stack that are greater than or equal to the current element
        while (top >= 0 && stack[top][0] >= heights[i]) {
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
int* nearestSmallerToRight(int* heights, int size, int pseudo) {
    int* v1 = (int*)malloc(size * sizeof(int));
    int stack[size][2];  // Stack to store elements and their indices
    int top = -1;

    for (int i = size - 1; i >= 0; i--) {
        // Pop elements from stack that are greater than or equal to the current element
        while (top >= 0 && stack[top][0] >= heights[i]) {
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
int subregionArea(int* subRegions, int size, int* leftCoord, int* rightCoord, int row) {
    int* left = nearestSmallerToLeft(subRegions, size, -1);     // Indices of nearest smaller elements to the left
    int* right = nearestSmallerToRight(subRegions, size, size); // Indices of nearest smaller elements to the right
    int max_area = 0;

    // Iterate through each element in the row to calculate potential area
    for (int i = 0; i < size; i++) {
        int width = right[i] - left[i] - 1;  // Calculate width of the rectangle
        int area = subRegions[i] * width;    // Calculate area for this rectangle

        // Update max_area and coordinates if a new maximum area is found
        if (area > max_area) {
            max_area = area;
            *leftCoord = left[i] + 1;        // Left column of the rectangle
            *rightCoord = right[i] - 1;      // Right column of the rectangle
        }
    }

    // Free allocated memory for left and right arrays
    free(left);
    free(right);
    return max_area;
}

// Main function to find the maximal rectangular region in a binary matrix
int maximalRegion(int** region, int x, int y, int* topLeftX, int* topLeftY, int* bottomRightX, int* bottomRightY) {
    int* temp = (int*)malloc(y * sizeof(int));  // Temp array to store heights of histogram columns
    
    // Initialize temp with the first row of region
    for (int j = 0; j < y; j++) {
        temp[j] = region[0][j];
    }

    // Initialize coordinates for the maximal region
    *topLeftX = *topLeftY = *bottomRightX = *bottomRightY = 0;
    
    int maxArea = subregionArea(temp, y, topLeftY, bottomRightY, 0);  // Calculate initial max area for the first row

    // Process each subsequent row to calculate the maximal rectangular region
    for (int i = 1; i < x; i++) {
        for (int j = 0; j < y; j++) {
            // Update temp array based on the current row; reset height if 0, otherwise accumulate height
            temp[j] = (region[i][j] == 0) ? 0 : temp[j] + region[i][j];
        }

        // Track the left-most and right-most columns of the max area subregion
        int leftCoord, rightCoord;
        int area = subregionArea(temp, y, &leftCoord, &rightCoord, i);
        
        // Update max area and coordinates if a new maximum area is found
        if (area > maxArea) {
            maxArea = area;
            *topLeftX = i - temp[leftCoord] + 1; // Top row of the rectangle
            *topLeftY = leftCoord;               // Left column
            *bottomRightX = i;                   // Bottom row
            *bottomRightY = rightCoord;          // Right column
        }
    }

    free(temp);  // Free allocated memory for temp array
    return maxArea;
}

int main()
{
    // input sample grid size.
    n = 4;
    m = 4;

    // sample grid

    int temp[4][4] = {
        {1, 2, 3, 4},
        {5, 6, 7, 8},
        {9, 10, 11, 12},
        {13, 14, 15, 16}
    };

    // copying the grid to the declared global grid.
    threshold = 10; // Example threshold
    // Copy example grid to global grid
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            grid[i][j] = temp[i][j];
        }
    }

    //building the segment tree for the input grid.
    build(1, 0, n - 1);

    // performing sum query 
    printf("Sum of subgrid (1, 1) to (3, 3): %d\n", sumQuery(1, 0, n - 1, 1, 3, 1, 3));

    // performing max query 
    printf("Max value in subgrid (1, 1) to (3, 3): %d\n", maxQuery(1, 0, n - 1, 1, 3, 1, 3));

    // update a value and query again
    update(1, 0, n - 1, 2, 2, 20); // update grid[2][2] to 20
    // sum
    printf("Sum of subgrid (1, 1) to (3, 3) after update: %d\n", sumQuery(1, 0, n - 1, 1, 3, 1, 3));
    // max
    printf("Max value in subgrid (1, 1) to (3, 3) after update: %d\n", maxQuery(1, 0, n - 1, 1, 3, 1, 3));
    int** markedArray = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++) {
        markedArray[i] = (int*)malloc(m * sizeof(int));
        for (int j = 0; j < m; j++) {
            markedArray[i][j] = marked[i][j];
        }
    }
        printf("markedArray:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            printf("%d ", markedArray[i][j]);
        }
        printf("\n");
    }


    int topLeftX, topLeftY, bottomRightX, bottomRightY;
    int maxArea = maximalRegion(markedArray, n, m, &topLeftX, &topLeftY, &bottomRightX, &bottomRightY);

    printf("Maximal Rectangle Area: %d\n", maxArea);
    printf("Top-left corner: (%d, %d)\n", topLeftX, topLeftY);
    printf("Bottom-right corner: (%d, %d)\n", bottomRightX, bottomRightY);

    int numRegions = countRegions(markedArray , n , m);
    printf("%d\n" , numRegions);

    // Free allocated memory for markedArray
    for (int i = 0; i < n; i++) {
        free(markedArray[i]);
    }
    free(markedArray);

    return 0;
}
