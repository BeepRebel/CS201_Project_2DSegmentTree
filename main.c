#include <stdio.h>
#include <limits.h>

#define MAX_ROWS 1000 //defined the maximum no. of rows.
#define MAX_COLS 1000 // defined the maximum no. of columns.

int n, m;  //dimensions of the grid.
int RangeSum[4 * MAX_ROWS][4 * MAX_COLS]; // for sum queries.
int RangeMax[4 * MAX_ROWS][4 * MAX_COLS]; // for max queries.
int grid[MAX_ROWS][MAX_COLS];  //input grid defined globally.

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
    
    return 0;
}
