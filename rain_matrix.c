#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define GRID_SIZE 100
#define MAX_LINE_LENGTH 256
#define MISSING_VALUE -999

// Function to read CSV file and populate rain matrix
void read_csv(const char *filename, double matrix[GRID_SIZE][GRID_SIZE]) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Unable to open file");
        exit(EXIT_FAILURE);
    }

    char line[MAX_LINE_LENGTH];
    double rain;
    int n = 0; // Counter for filling the grid

    // Initialize the rain_matrix with missing values
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            matrix[i][j] = MISSING_VALUE; // Initialize with MISSING_VALUE
        }
    }

    // Read each line in the CSV file, skipping the header
    fgets(line, sizeof(line), file); // Skip header line
    while (fgets(line, sizeof(line), file)) {
        // Parse the rain value from the CSV line
        if (sscanf(line, "%*lf,%*lf,%lf", &rain) == 1) {
            // Calculate row and column based on counter n
            int row = n / GRID_SIZE;
            int col = n % GRID_SIZE;

            // Check if indices are within the matrix bounds
            if (row < GRID_SIZE && col < GRID_SIZE) {
                matrix[row][col] = rain; // Store rain value
                n++; // Increment the counter
            } else {
                fprintf(stderr, "Index out of bounds for n: %d\n", n);
                break; // Break the loop if out of bounds
            }
        } else {
            fprintf(stderr, "Failed to parse line: %s\n", line);
        }
    }

    fclose(file);
}

// Function to print the rain matrix
void print_rain_matrix(double rain_matrix[GRID_SIZE][GRID_SIZE]) {
    for (int i = 0; i < GRID_SIZE; i++) {
        for (int j = 0; j < GRID_SIZE; j++) {
            if (rain_matrix[i][j] == MISSING_VALUE) {
                printf("  -999 "); // Indicate missing value
            } else {
                printf("%7.2f ", rain_matrix[i][j]);
            }
        }
        printf("\n");
    }
}

int main() {
    double rain_matrix[GRID_SIZE][GRID_SIZE];

    // Read the CSV data and populate the rain matrix
    read_csv("combined_data_2019.csv", rain_matrix);

    // Print the rain matrix
    print_rain_matrix(rain_matrix);

    return 0;
}
