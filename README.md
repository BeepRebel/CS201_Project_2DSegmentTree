# CS201 Project (Group 13)
# Geographical Region Statistics for Environmental Monitoring Using 2D Segment Tree

## Project Overview

This project uses a **2D Segment Tree** for efficient storage, querying, and updating of environmental data such as **temperature, humidity,** and **rainfall** across a grid-based geographical region. Designed to support real-time environmental monitoring, it provides flexible tools for dynamic data analysis. Key features include range queries, update operations, and complex threshold-based queries that make it an ideal solution for weather forecasting and disaster risk prediction.

## Authors

- **Charvi Pahuja** - (2023CSB1114)
- **Krishna Agarwal** - (2023CSB1131)
- **Tamanna** - (2023CSB1169)

Mentoring TA
- **Shubham Thawait**


## Files Included

The project includes two data files:
- **`data_2022.txt`**: Environmental data for the year 2022
- **`data_2023.txt`**: Environmental data for the year 2023
- **`main.c`**; contains the main code

The program reads data from the selected year’s file and builds a 2D Segment Tree to support operations based on this structured data.

---

## How to Run 
To clone the repository to your local machine, run:
```bash
git clone https://github.com/your-username/your-repository.git
cd your-repository
```
After navigating into the cloned directory, compile the C program main.c using:
```bash
gcc main.c -o main
```
Once compiled, you can run the program with the following command:
```
./main
```
---
## Usage

After selecting the data year and starting the program, you can perform various query operations to interact with the environmental data. Below are the query types supported by the program, along with their inputs and expected outputs.

### 1. Range Query
Retrieve statistics for a specified environmental attribute within a selected region.

- **Input**:
  - Query Type: `Range Query`
  - Coordinates: Define the region with starting (x1, y1) and ending (x2, y2) points with latitude and longitudes.
  - Statistic: Choose a statistic (e.g., `sum`, `min`, `max`, or `average`).
  - Attribute: Select an attribute (e.g., `temperature`, `humidity`, or `rainfall`).
- **Output**: Returns the calculated statistic for the specified attribute within the selected region.

**Example**:
```plaintext
Query Type: Range
Range Query Type: Sum
Coordinates: (9.20, 69.20) to (9.80, 69.60)
Result: 342.70
Time Taken: 0.000 seconds
Attribute: Temperature
```
### 2. Update Query

Modify values for specific attributes at a given coordinate.

- **Input**:
  - Query Type: `Update`
  - Coordinates: Specify the point (x, y) where the update is to be applied.
  - New Value: Provide a new value for any attribute (e.g., `temperature`, `humidity`, or `rainfall`).
- **Output**: Confirms that the attribute was updated in the 2D grid and the segment tree, allowing subsequent queries to reflect the change.

### 3. Combined Threshold Query

Identify regions meeting specified threshold conditions for temperature, humidity, and rainfall.

- **Input**:
  - Query Type: `Combined Query`
  - Thresholds: Set threshold ranges for each attribute:
- **Output**: Provides:
  - The maximum area of regions that satisfy the thresholds.
  - The count of regions meeting the criteria.
  - Coordinates of one such region.
    
### Output Summary

The final output for the queries is saved in a file named `output_summary.txt`. This file provides a summary of the query results, including any calculated statistics, updated values, and details of regions that meet specified criteria. 

Be sure to check `output_summary.txt` after running your queries to view the complete results.
