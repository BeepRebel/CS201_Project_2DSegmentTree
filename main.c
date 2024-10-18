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

    return 0;
}
