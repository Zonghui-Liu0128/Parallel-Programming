#include <mpi.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>

using namespace std;
//*************************************************************************************
//		MPI Programming Coursework-Solving the Wave Equation
//*************************************************************************************
/*
	This program uses a finite difference method to solve
	Wave's equation using parallel programming technology(MPI).
//
//  Author:
//  Zonghui Liu
//*************************************************************************************
*/

//	Global Parameters:
//
//    Global, int id, the MPI process id.
//    Global, int p, the number of MPI processes.
//    Global, clock_t mpi_start, mpi_end, to get the current time.
//    Global, int imax = 100, jmax = 100, the total number of points in x direction and y direction.
//	  Global, double x_max,y_max, maximum size of x direction and y direction.
//	  Global, double t_max, maximum of time.
//    Global, vector <vector<double>> grid, new_grid, old_grid, the status of the current node(3 time steps)
//    Global, double dt, the time step for iteration.
//    Global, double dt_out, the time step for writing to file.
//    Global, double t, t_out, the current time for iteration and writing to file.
//    Global, double dx, dy, the step in x direction and y direction.
//	  Global, int rows, the number of rows in MPI processes
//	  Global, int columns, the number of columns in MPI processes
//	  Global, double c, the parameter in the Wave's equation
//	  Global, int id_row, the MPI process id in row.
//	  Global, int id_column, the MPI process id in column.
//	  MPI_Datatype for Send data to neighbors : Datatype_left, Datatype_right, Datatype_top, Datatype_bottom;
//	  MPI_Datatype for Send data to neighbors : Datatype_left_recv, Datatype_right_recv, Datatype_top_recv, Datatype_bottom_recv;
//    bool is_Neumann, is_periodic, the flag to choose boundary conditions
//--------------------------------------------------------------------------------------------------------------------------------
//
int id, p;
clock_t mpi_start, mpi_end;
int imax = 100, jmax = 100;
double x_max = 10.0, y_max = 10.0;
double t_max = 30.0;
vector <vector<double>> grid, new_grid, old_grid;
double t, dt, t_out = 0.0, dt_out = 0.04;
double dx, dy;
double c = 1;
int rows, columns;
int id_row, id_column;
MPI_Datatype Datatype_left, Datatype_right, Datatype_top, Datatype_bottom;
MPI_Datatype Datatype_left_recv, Datatype_right_recv, Datatype_top_recv, Datatype_bottom_recv;
bool is_Neumann = false;
bool is_periodic = false;

/**
 * Todo: Divide the whole domain into several nodes and
 *       get the number of rows and columns of the MPI processors
 * @param p: the number of MPI processors
 * @param rows: the number of rows for the MPI processors
 * @param columns: the number of columns for the MPI processors
 */
void find_dimensions(int p, int &rows, int &columns) {
    int min_gap = p;
    int top = sqrt(p) + 1;
    for (int i = 1; i <= top; i++) {
        if (p % i == 0) {
            int gap = abs(p / i - i);
            if (gap < min_gap) {
                min_gap = gap;
                rows = i;
                columns = p / i;
            }
        }
    }
//    if (id == 0) {
//        cout << "The total number of processor is: " << p << "." << endl;
//        cout << "Divide " << p << " into " << rows << " by " << columns << " grid" << endl;
//    }
}

/**
 * Todo: Assign each node to each processor and
 *       get the id_row and id_column according to the id
 * @param id: the MPI processor rank
 * @param id_row: the MPI processor id in row
 * @param id_column: the MPI processor id in column
 */
void id_to_index(int id, int &id_row, int &id_column) {
    id_column = id % columns;
    id_row = id / columns;
}

/**
 * Todo: Get the id of MPI processor according to the MPI processor id in row and that in column
 * @param id_row: the MPI processor id in row
 * @param id_column: the MPI processor id in column
 * @return
 */
int id_from_index(int id_row, int id_column) {
    if (id_row >= rows || id_row < 0)
        return -1;
    if (id_column >= columns || id_column < 0)
        return -1;
    return id_row * columns + id_column;
}

/**
 * Todo: Allocate each node with their own data
 * @param imax_current: the number of rows of the current node with ghost layer
 * @param jmax_current: the number of columns of the current node with ghost layer
 */
void Allocate_to_node(int &imax_current, int &jmax_current) {
    imax_current = imax / rows + 2;
    jmax_current = jmax / columns + 2;
    int rows_remain = imax % rows;
    int columns_remain = jmax % columns;
    int imax_current_remain = imax_current + rows_remain;
    int jmax_current_remain = jmax_current + columns_remain;

    // The last row and last column node has the biggest size,
    // all the remaining points will be added to this node
    if (id_row == rows - 1 && id_column == columns - 1) {
        old_grid.resize(imax_current_remain, vector<double>(jmax_current_remain));
        grid.resize(imax_current_remain, vector<double>(jmax_current_remain));
        new_grid.resize(imax_current_remain, vector<double>(jmax_current_remain));

    }
        // The nodes in the last row,
        // the remaining points in y-direction will be added into these nodes
    else if (id_row == rows - 1) {
        old_grid.resize(imax_current_remain, vector<double>(jmax_current));
        grid.resize(imax_current_remain, vector<double>(jmax_current));
        new_grid.resize(imax_current_remain, vector<double>(jmax_current));

    }
        // The nodes in the last column,
        // the remaining points in x-direction will be added into these nodes
    else if (id_column == columns - 1) {
        old_grid.resize(imax_current, vector<double>(jmax_current_remain));
        grid.resize(imax_current, vector<double>(jmax_current_remain));
        new_grid.resize(imax_current, vector<double>(jmax_current_remain));
    }
        // The internal nodes
    else {
        old_grid.resize(imax_current, vector<double>(jmax_current));
        grid.resize(imax_current, vector<double>(jmax_current));
        new_grid.resize(imax_current, vector<double>(jmax_current));
    }
}

/**
 * Todo: Initialize the system
 * @param dx: the step in x direction
 * @param dy: the step in y direction
 * @param imax_current: the number of rows of the current node with ghost layer
 * @param jmax_current: the number of columns of the current node with ghost layer
 */
void Initialize(double dx, double dy, int &imax_current, int &jmax_current) {
    int imax_nude = imax_current - 2;
    int jmax_nude = jmax_current - 2;

    double r_splash = 1.0;
    double x_splash = 3.0;
    double y_splash = 3.0;

    // Note: we don't want to initialize the points in the ghost layer
    for (int i = 1; i < imax_current - 1; i++)
        for (int j = 1; j < jmax_current - 1; j++) {
            double x = dx * (i - 1) + dx * imax_nude * id_row;
            double y = dy * (j - 1) + dy * jmax_nude * id_column;
            double dist = sqrt(pow(x - x_splash, 2.0) + pow(y - y_splash, 2.0));
            if (dist < r_splash) {
                double h = 5.0 * (cos(dist / r_splash * M_PI) + 1.0);
                // Initialize those points that receive the disturbance to h
                grid[i][j] = h;
                old_grid[i][j] = h;
            } else {
                // Initialize those points that did not receive the disturbance to zero
                grid[i][j] = 0;
                old_grid[i][j] = 0;
            }
        }
}

/**
 * Todo: Write the status of each points to files
 * @param id: the MPI processor rank
 * @param out: the count of output
 * @param grid: the grid we want to write
 */
void grid_to_file(int id, int out, vector <vector<double>> &grid) {
    int imax_current = grid.size();
    int jmax_current = grid[0].size();
    stringstream fname;
    fstream f1;
    // Note: the path could be changed into an absolute path
    fname << "/Users/lzh/CLionProjects/WaveEquation/out/output" << "_" << out << "_p" << id << ".dat";
    f1.open(fname.str().c_str(), ios_base::out);
    // Note: we only need the data without ghost points
    for (int i = 1; i < imax_current - 1; i++) {
        for (int j = 1; j < jmax_current - 1; j++)
            f1 << grid[i][j] << "\t";
        f1 << endl;
    }
    f1.close();
}

/**
 * Todo: Create MPI data types for sending data
 * @param send_data: the buffer for sending data
 * @param m: the number of rows of grid
 * @param n: the number of columns of grid
 */
void createdatatypes_send(double **send_data, int m, int n) {
    vector<int> block_lengths;
    vector <MPI_Datatype> typelist;
    vector <MPI_Aint> addresses;
    MPI_Aint add_start;

    // The left boundary(the id of column is zero)
    for (int i = 0; i < m; i++) {
        block_lengths.push_back(1);
        typelist.push_back(MPI_DOUBLE);
        MPI_Aint temp_address;
        MPI_Get_address(&send_data[i][0], &temp_address);
        addresses.push_back(temp_address);
    }
    MPI_Get_address(send_data, &add_start);
    for (int i = 0; i < m; i++) addresses[i] = addresses[i] - add_start;
    MPI_Type_create_struct(m, block_lengths.data(), addresses.data(), typelist.data(), &Datatype_left);
    MPI_Type_commit(&Datatype_left);

    // The right boundary(the id of column is n - 1)
    block_lengths.resize(0);
    typelist.resize(0);
    addresses.resize(0);
    for (int i = 0; i < m; i++) {
        block_lengths.push_back(1);
        typelist.push_back(MPI_DOUBLE);
        MPI_Aint temp_address;
        MPI_Get_address(&send_data[i][n - 1], &temp_address);
        addresses.push_back(temp_address);
    }
    for (int i = 0; i < m; i++) addresses[i] = addresses[i] - add_start;
    MPI_Type_create_struct(m, block_lengths.data(), addresses.data(), typelist.data(), &Datatype_right);
    MPI_Type_commit(&Datatype_right);

    // The top boundary(the id of row is zero) - only need one value
    int block_length = n;
    MPI_Datatype typeval = MPI_DOUBLE;
    MPI_Aint address;
    MPI_Get_address(send_data[0], &address);
    address = address - add_start;
    MPI_Type_create_struct(1, &block_length, &address, &typeval, &Datatype_top);
    MPI_Type_commit(&Datatype_top);

    // The bottom boundary(the id of row is m - 1) - only need one value
    MPI_Get_address(send_data[m - 1], &address);
    address = address - add_start;
    MPI_Type_create_struct(1, &block_length, &address, &typeval, &Datatype_bottom);
    MPI_Type_commit(&Datatype_bottom);
}

/**
 * Todo: Create MPI data types for sending data
 * @param recv_data: the buffer for receiving data
 * @param m: the number of rows of grid
 * @param n: the number of columns of grid
 */
void createdatatypes_recv(double **recv_data, int m, int n) {
    vector<int> block_lengths;
    vector <MPI_Datatype> typelist;
    vector <MPI_Aint> addresses;
    MPI_Aint add_start;

    // The left boundary(the id of column is zero)
    // Note: there are ghost points in the recv_data buffer, need to pay attention to the index
    for (int i = 1; i < m - 1; i++) {
        block_lengths.push_back(1);
        typelist.push_back(MPI_DOUBLE);
        MPI_Aint temp_address;
        MPI_Get_address(&recv_data[i][0], &temp_address);
        addresses.push_back(temp_address);
    }
    MPI_Get_address(recv_data, &add_start);
    for (int i = 0; i < m - 2; i++) addresses[i] = addresses[i] - add_start;
    MPI_Type_create_struct(m - 2, block_lengths.data(), addresses.data(), typelist.data(), &Datatype_left_recv);
    MPI_Type_commit(&Datatype_left_recv);

    // The right boundary(the id of column is n - 1)
    block_lengths.resize(0);
    typelist.resize(0);
    addresses.resize(0);
    for (int i = 1; i < m - 1; i++) {
        block_lengths.push_back(1);
        typelist.push_back(MPI_DOUBLE);
        MPI_Aint temp_address;
        MPI_Get_address(&recv_data[i][n - 1], &temp_address);
        addresses.push_back(temp_address);
    }
    for (int i = 0; i < m - 2; i++) addresses[i] = addresses[i] - add_start;
    MPI_Type_create_struct(m - 2, block_lengths.data(), addresses.data(), typelist.data(), &Datatype_right_recv);
    MPI_Type_commit(&Datatype_right_recv);

    // The top boundary - only need one value
    // Note: we only need 4 neighbours to communicate, so the index should begin with 1 instead of 0
    int block_length = n - 2;
    MPI_Datatype typeval = MPI_DOUBLE;
    MPI_Aint address;
    MPI_Get_address(&recv_data[0][1], &address);
    address = address - add_start;
    MPI_Type_create_struct(1, &block_length, &address, &typeval, &Datatype_top_recv);
    MPI_Type_commit(&Datatype_top_recv);

    // The bottom boundary - only need one value
    // Note: we only need 4 neighbours to communicate, so the index should begin with m-1 and 1 instead of m and 0
    MPI_Get_address(&recv_data[m - 1][1], &address);
    address = address - add_start;
    MPI_Type_create_struct(1, &block_length, &address, &typeval, &Datatype_bottom_recv);
    MPI_Type_commit(&Datatype_bottom_recv);
}

/**
 * Todo: Copy the data from the receiving buffer to the current grid
 * @param grid: the current grid we want to update
 * @param rows_of_recv: the number of rows of the receiving array
 * @param cols_of_recv: the number of columns of the receiving array
 * @param recv_data: the receiving buffer(2-D array)
 * @param recv_left: the flag to check if there are some data from the left node
 * @param recv_top: the flag to check if there are some data from the top node
 * @param recv_right: the flag to check if there are some data from the right node
 * @param recv_bottom: the flag to check if there are some data from the bottom node
 */
void copy_boundry(vector <vector<double>> &grid, int rows_of_recv, int cols_of_recv, double **recv_data,
                  bool recv_left, bool recv_top, bool recv_right, bool recv_bottom) {
    // Get the size of the current node
    int imax_current = grid.size();
    int jmax_current = grid[0].size();
    // Left
    if (recv_left) {
        // Copy data from the buffer(recv_data) to the left ghost layer
        // Note: we don't need the ghost layer of neighbouring node
        for (int i = 1; i < imax_current - 1; i++) {
            grid[i][0] = recv_data[i][0];
        }
    }
    // Right
    if (recv_right) {
        // Copy data from the buffer(recv_data) to the right ghost layer
        // Note: we don't need the ghost layer of neighbouring node
        for (int i = 1; i < imax_current - 1; i++) {
            grid[i][jmax_current - 1] = recv_data[i][cols_of_recv - 1];
        }
    }
    // Top
    if (recv_top) {
        // Copy data from the buffer(recv_data) to the top ghost layer
        // Note: we don't need the ghost layer of neighbouring node
        for (int i = 1; i < jmax_current - 1; i++) {
            grid[0][i] = recv_data[0][i];
        }
    }
    // Bottom
    if (recv_bottom) {
        // Copy data from the buffer(recv_data) to the bottom ghost layer
        // Note: we don't need the ghost layer of neighbouring node
        for (int i = 1; i < jmax_current - 1; i++) {
            grid[imax_current - 1][i] = recv_data[rows_of_recv - 1][i];
        }
    }
}

/**
 * Todo: Do iteration, update the status of each point
 * @param is_Neumann: the flag to check whether the boundary condition is Neumann boundary conditions or not
 */
void do_iteration(bool is_Neumann) {
    //gird size (include ghost layer)
    int imax_current = grid.size();
    int jmax_current = grid[0].size();

    // The current node is an internal node(don't need to consider the boundary)
    if (id_row != 0 && id_row != rows - 1 && id_column != 0 && id_column != columns - 1) {
        // Note: the imax_current and jmax_current include the ghost layer, but we don't want to consider the ghost layer
        for (int i = 1; i < imax_current - 1; i++)
            for (int j = 1; j < jmax_current - 1; j++)
                // Update the status of the current grid
                new_grid[i][j] = pow(dt * c, 2.0) *
                                 ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) +
                                  (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) +
                                 2.0 * grid[i][j] - old_grid[i][j];
    }

    // The current node is at the top
    if (id_row == 0) {
        // The current node is at the top and the left
        if (id_column == 0) {
            // Note: we don't need to update the boundary points,
            // they will be sent from the neighbouring nodes
            for (int i = 2; i < imax_current - 1; i++)
                for (int j = 2; j < jmax_current - 1; j++) {
                    new_grid[i][j] = pow(dt * c, 2.0) *
                                     ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) +
                                      (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) +
                                     2.0 * grid[i][j] - old_grid[i][j];
                }
        }
            // The current node is at the top and the right
        else if (id_column == columns - 1) {
            // Note: we don't need to update the boundary points,
            // they will be sent from the neighbouring nodes
            for (int i = 2; i < imax_current - 1; i++)
                for (int j = 1; j < jmax_current - 2; j++) {
                    new_grid[i][j] = pow(dt * c, 2.0) *
                                     ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) +
                                      (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) +
                                     2.0 * grid[i][j] - old_grid[i][j];
                }
        } else {
            for (int i = 2; i < imax_current - 1; i++)
                for (int j = 1; j < jmax_current - 1; j++) {
                    new_grid[i][j] = pow(dt * c, 2.0) *
                                     ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) +
                                      (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) +
                                     2.0 * grid[i][j] - old_grid[i][j];
                }
        }
    }

    // The current node is at the bottom
    if (id_row == rows - 1) {
        // The current node is at the bottom and the left
        if (id_column == 0) {
            // Note: we don't need to update the boundary points,
            // they will be sent from the neighbouring nodes
            for (int i = 1; i < imax_current - 2; i++)
                for (int j = 2; j < jmax_current - 1; j++) {
                    new_grid[i][j] = pow(dt * c, 2.0) *
                                     ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) +
                                      (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) +
                                     2.0 * grid[i][j] - old_grid[i][j];
                }
        }
            // The current node is at the bottom and the right
        else if (id_column == columns - 1) {
            // Note: we don't need to update the boundary points,
            // they will be sent from the neighbouring nodes
            for (int i = 1; i < imax_current - 2; i++)
                for (int j = 1; j < jmax_current - 2; j++) {
                    new_grid[i][j] = pow(dt * c, 2.0) *
                                     ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) +
                                      (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) +
                                     2.0 * grid[i][j] - old_grid[i][j];
                }
        } else {
            // Note: we don't need to update the boundary points,
            // they will be sent from the neighbouring nodes
            for (int i = 1; i < imax_current - 2; i++)
                for (int j = 1; j < jmax_current - 1; j++) {
                    new_grid[i][j] = pow(dt * c, 2.0) *
                                     ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) +
                                      (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) +
                                     2.0 * grid[i][j] - old_grid[i][j];
                }
        }
    }

    // The current node is at the left
    if (id_column == 0) {
        for (int i = 1; i < imax_current - 1; ++i) {
            for (int j = 2; j < jmax_current - 1; ++j) {
                new_grid[i][j] = pow(dt * c, 2.0) *
                                 ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) +
                                  (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) +
                                 2.0 * grid[i][j] - old_grid[i][j];
            }
        }
    }

    // The current node is at the right
    if (id_column == columns - 1) {
        for (int i = 1; i < imax_current - 1; ++i) {
            for (int j = 1; j < jmax_current - 2; ++j) {
                new_grid[i][j] = pow(dt * c, 2.0) *
                                 ((grid[i + 1][j] - 2.0 * grid[i][j] + grid[i - 1][j]) / pow(dx, 2.0) +
                                  (grid[i][j + 1] - 2.0 * grid[i][j] + grid[i][j - 1]) / pow(dy, 2.0)) +
                                 2.0 * grid[i][j] - old_grid[i][j];
            }
        }
    }

    /** Set Neumann boundaries according to process position **/
    if (is_Neumann) {
        // if it is in top boundary
        if (id_row == 0) {
            // Set the top Neumann boundaries
            for (int j = 1; j < jmax_current - 1; j++) {
                new_grid[1][j] = new_grid[2][j];
            }
        }
        // if it is in bottom boundary
        if (id_row == rows - 1) {
            // Set the bottom Neumann boundaries
            for (int j = 1; j < jmax_current - 1; j++) {
                new_grid[imax_current - 2][j] = new_grid[imax_current - 3][j];
            }
        }
        // if it is in left boundary
        if (id_column == 0) {
            // Set the left Neumann boundaries
            for (int i = 1; i < imax_current - 1; i++) {
                new_grid[i][1] = new_grid[i][2];
            }
        }
        // if it is in right boundary
        if (id_column == columns - 1) {
            // Set the right Neumann boundaries
            for (int i = 1; i < imax_current - 1; i++) {
                new_grid[i][jmax_current - 2] = new_grid[i][jmax_current - 3];
            }
        }
    }
        /** Set Dirichlet boundaries according to process position **/
    else {
        // if it is in top boundary
        if (id_row == 0) {
            // Set the top Neumann boundaries
            for (int j = 1; j < jmax_current - 1; j++) {
                new_grid[1][j] = 0;
            }
        }
        // if it is in bottom boundary
        if (id_row == rows - 1) {
            // Set the bottom Neumann boundaries
            for (int j = 1; j < jmax_current - 1; j++) {
                new_grid[imax_current - 2][j] = 0;
            }
        }
        // if it is in left boundary
        if (id_column == 0) {
            // Set the left Neumann boundaries
            for (int i = 1; i < imax_current - 1; i++) {
                new_grid[i][1] = 0;
            }
        }
        // if it is in right boundary
        if (id_column == columns - 1) {
            // Set the right Neumann boundaries
            for (int i = 1; i < imax_current - 1; i++) {
                new_grid[i][jmax_current - 2] = 0;
            }
        }
    }
    // Update the time step
    t += dt;

    // Update the status of grid
    old_grid.swap(new_grid);
    old_grid.swap(grid);
}

/**
 * Todo: Write the status of each points to the terminal(Test)
 * @param grid: the grid we want to print
 */
void print_grid(vector <vector<double>> &grid) {
    for (int i = 0; i < grid.size(); ++i) {
        for (int j = 0; j < grid[0].size(); ++j) {
            cout << grid[i][j] << "\t";
            cout.flush();
        }
        cout << endl;
        cout.flush();
    }
}

/**
 * Todo: Find the id of neighboring nodes
 * @param id: the id of the current node
 * @param id_row: the id of row in the MPI processor
 * @param id_column: the id of column in the MPI processor
 * @param id_left: the id of the left neighbour
 * @param id_right: the id of the right neighbour
 * @param id_top: the id of the top neighbour
 * @param id_bottom: the id of the bottom neighbour
 */
void get_neighbour_ids(int id, int id_row, int id_column, int &id_left, int &id_right, int &id_top, int &id_bottom) {

    // The internal node
    if (id_row > 0 && id_row < rows - 1 && id_column > 0 && id_column < columns - 1) {
        id_left = id - 1;
        id_right = id + 1;
        id_top = id + columns;
        id_bottom = id - columns;
    }
        // The most bottom nodes
    else if (id_row == 0 && id_column < columns - 1 && id_column > 0) {
        id_left = id - 1;
        id_right = id + 1;
        id_bottom = id + (rows - 1) * columns;
        id_top = id + columns;
    }

        // The most top nodes
    else if (id_row == rows - 1 && id_column < columns - 1 && id_column > 0) {
        id_left = id - 1;
        id_right = id + 1;
        id_bottom = id - columns;
        id_top = id - (rows - 1) * columns;
    }

        //	The most left nodes
    else if (id_column == 0 && id_row < rows - 1 && id_row > 0) {
        id_left = id + columns - 1;
        id_right = id + 1;
        id_bottom = id - columns;
        id_top = id + columns;
    }

        // The most right nodes
    else if (id_column == columns - 1 && id_row < rows - 1 && id_row > 0) {
        id_left = id - 1;
        id_right = id - columns + 1;
        id_bottom = id - columns;
        id_top = id + columns;
    }

        //	The most top and most right node
    else if (id_column == columns - 1 && id_row == rows - 1) {
        id_left = id - 1;
        id_right = id - columns + 1;
        id_bottom = id - columns;
        id_top = id - (rows - 1) * columns;
    }

        //	The most top and most left node
    else if (id_column == 0 && id_row == rows - 1) {
        id_left = id + columns - 1;
        id_right = id + 1;
        id_bottom = id - columns;
        id_top = id - (rows - 1) * columns;
    }

        //  The most bottom and most right node
    else if (id_row == 0 && id_column == columns - 1) {
        id_left = id - 1;
        id_right = id - columns + 1;
        id_bottom = id + (rows - 1) * columns;
        id_top = id + columns;
    }

        //  The most bottom and most left node
    else if (id_row == 0 && id_column == 0) {
        id_left = id + columns - 1;
        id_right = id + 1;
        id_bottom = id + (rows - 1) * columns;
        id_top = id + columns;
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    mpi_start = clock();
//    cout << "Processor " << id << " in " << p << endl;
//    cout.flush();

    /** Divide the domain into several nodes, the number of nodes = the number of processors **/
    find_dimensions(p, rows, columns);
//    if (id == 0) {
//        cout << "The number of processors(" << p << ") can be divided into " << rows << " * " << columns << endl;
//    }

    /** Get the row-id and column-id based on the process id for each node/processor **/
    id_to_index(id, id_row, id_column);
//    cout << id << ": " << "(" << id_row << "," << id_column << ")" << endl;
//    cout.flush();

    // the size of the current node with ghost layer
    int imax_current, jmax_current;

    /** Allocate each node to each process based on id **/
    Allocate_to_node(imax_current, jmax_current);
//    cout << id << ": " << "The size of current node is " << grid.size() << " * " << grid[0].size() << endl;
//    cout.flush();

    dx = x_max / ((double) imax - 1);
    dy = y_max / ((double) jmax - 1);

    t = 0.0;

    dt = 0.1 * min(dx, dy) / c;

    //output time index，iterate time index in iteration
    int out_cnt = 0, it = 0;

    /** Find the neighbours **/
    int id_left = 0;
    int id_right = 0;
    int id_top = 0;
    int id_bottom = 0;
    get_neighbour_ids(id, id_row, id_column, id_left, id_right, id_top, id_bottom);
//    cout << id << ": " << id_top << "-" << id_bottom << endl;
//    cout << id << ": " << id_left << "-" << id_right << "-" << id_top << "-" << id_bottom << endl;
//    cout.flush();

    /** Initialize the system: set half sinusoidal disturbance in the whole domain **/
    Initialize(dx, dy, imax_current, jmax_current);
    // Wrtie the initial output to the file
    grid_to_file(id, out_cnt, grid);
    out_cnt++;
    t_out += dt_out;

    /** Initialize the sending and receiving buffer **/
    double **send_data;
    double **recv_data;
    // !Note: should use grid.size() instead of imax_current, for the current node's size
    int rows_of_send_data = grid.size() - 2;
    int cols_of_send_data = grid[0].size() - 2;
    send_data = new double *[rows_of_send_data];
    for (int i = 0; i < rows_of_send_data; i++) {
        send_data[i] = new double[cols_of_send_data];
    }

    // To avoid data overflow errors,
    // the maximum size of the nodes will be used for initializing the 2-D array(recv_data)
    int rows_remain = imax % rows;
    int columns_remain = jmax % columns;
    int rows_of_recv_data = imax_current + rows_remain;
    int cols_of_recv_data = jmax_current + columns_remain;
    recv_data = new double *[rows_of_recv_data];
    for (int i = 0; i < rows_of_recv_data; i++) {
        recv_data[i] = new double[cols_of_recv_data];
    }

//    cout << id << ": " << endl;
//    cout.flush();
//    cout << "The size of the send_data buffer: " << rows_of_send_data << " * " << cols_of_send_data << endl;
//    cout.flush();
//    cout << "The size of the recv_data buffer: " << rows_of_recv_data << " * " << cols_of_recv_data << endl;
//    cout.flush();

    /** Create the MPI data types for sending and receiving data from neighbouring nodes **/
    // the MPI data type for sending data to neighbouring nodes
    createdatatypes_send(send_data, rows_of_send_data, cols_of_send_data);
    // the MPI data type for receiving data from neighbouring nodes
    createdatatypes_recv(recv_data, rows_of_recv_data, cols_of_recv_data);
    if (is_periodic) {
        while (t < t_max) {
            MPI_Request *requests = nullptr;
            requests = new MPI_Request[4 * 2];
            // communicate for sending data to ghost cells
            int cnt = 0;
            MPI_Isend(send_data, 1, Datatype_left, id_left, it, MPI_COMM_WORLD, &requests[cnt]);
            cnt++;
            MPI_Isend(send_data, 1, Datatype_right, id_right, it, MPI_COMM_WORLD, &requests[cnt]);
            cnt++;
            MPI_Isend(send_data, 1, Datatype_top, id_top, it, MPI_COMM_WORLD, &requests[cnt]);
            cnt++;
            MPI_Isend(send_data, 1, Datatype_bottom, id_bottom, it, MPI_COMM_WORLD, &requests[cnt]);
            cnt++;

            MPI_Irecv(recv_data, 1, Datatype_top_recv, id_top, 2, MPI_COMM_WORLD, &requests[cnt]);
            cnt++;
            MPI_Irecv(recv_data, 1, Datatype_left_recv, id_left, 2, MPI_COMM_WORLD, &requests[cnt]);
            cnt++;
            MPI_Irecv(recv_data, 1, Datatype_right_recv, id_right, 2, MPI_COMM_WORLD, &requests[cnt]);
            cnt++;
            MPI_Irecv(recv_data, 1, Datatype_bottom_recv, id_bottom, 2, MPI_COMM_WORLD, &requests[cnt]);
            cnt++;

            // Copy the data of edges in recv_array to gird(only copy the edge it received, not all four edges).
            copy_boundry(grid, rows_of_recv_data, cols_of_recv_data, recv_data, true, true, true,
                         true);
            // Perform the calculation for new_gird using grid and old_grid, and set the boundry value(Neumann boundaries) if process is in boundry.
            do_iteration(is_Neumann);
            // Write the output if it should
            if (t_out <= t) {
                grid_to_file(id, out_cnt, grid);
                out_cnt++;
                t_out += dt_out;
            }
            it++;
            delete[] requests;
        }
    } else {
        while (t < t_max) {
            // Ready for send data
            int imax_nude = grid.size() - 2;
            int jmax_nude = grid[0].size() - 2;
            for (int i = 0; i < imax_nude; i++) {
                for (int j = 0; j < jmax_nude; j++) {
                    // Load the data from grid to the sending buffer
                    send_data[i][j] = (double) grid[i + 1][j + 1];
                }
            }
            // Initialize the recv_array
            for (int i = 0; i < rows_of_recv_data; i++)
                for (int j = 0; j < cols_of_recv_data; j++) {
                    // Initialize the recv buffer with -1
                    recv_data[i][j] = -1;
                }
            int cnt = 0;
            MPI_Request *request = new MPI_Request[4 * 2];
            // The flags to mark the position of communication
            bool recv_top = false;
            bool recv_left = false;
            bool recv_right = false;
            bool recv_bottom = false;

            // Iterate the neighbor processes
            for (int i = -1; i <= 1; i++)
                for (int j = -1; j <= 1; j++) {
                    // The process sended row index and column index
                    int com_i = id_row + i;
                    int com_j = id_column + j;
                    // Get the id of the sending destination
                    int com_id = id_from_index(com_i, com_j);
                    // If process sended is not current process id, and it is an existed id
                    if (com_id != id && com_id >= 0 && com_id < p) {
                        // Communication for the top boundary
                        if (com_i == id_row - 1 && com_j == id_column) {
                            MPI_Isend(send_data, 1, Datatype_top, com_id, it, MPI_COMM_WORLD, &request[cnt * 2]);
                            MPI_Irecv(recv_data, 1, Datatype_top_recv, com_id, it, MPI_COMM_WORLD,
                                      &request[cnt * 2 + 1]);
                            recv_top = true;
                            cnt++;
                        }
                            // Communication for the bottom boundary
                        else if (com_i == id_row + 1 && com_j == id_column) {
                            MPI_Isend(send_data, 1, Datatype_bottom, com_id, it, MPI_COMM_WORLD, &request[cnt * 2]);
                            MPI_Irecv(recv_data, 1, Datatype_bottom_recv, com_id, it, MPI_COMM_WORLD,
                                      &request[cnt * 2 + 1]);
                            recv_bottom = true;
                            cnt++;
                        }
                            // Communication for the left boundary
                        else if (com_i == id_row && com_j == id_column - 1) {
                            MPI_Isend(send_data, 1, Datatype_left, com_id, it, MPI_COMM_WORLD, &request[cnt * 2]);
                            MPI_Irecv(recv_data, 1, Datatype_left_recv, com_id, it, MPI_COMM_WORLD,
                                      &request[cnt * 2 + 1]);
                            recv_left = true;
                            cnt++;
                        }
                            // Communication for the right boundary
                        else if (com_i == id_row && com_j == id_column + 1) {
                            MPI_Isend(send_data, 1, Datatype_right, com_id, it, MPI_COMM_WORLD, &request[cnt * 2]);
                            MPI_Irecv(recv_data, 1, Datatype_right_recv, com_id, it, MPI_COMM_WORLD,
                                      &request[cnt * 2 + 1]);
                            recv_right = true;
                            cnt++;
                        }
                    }
                }
            // Wait for all sending and receiving for this current process finished
            MPI_Waitall(cnt * 2, request, MPI_STATUS_IGNORE);
            // Copy the data of edges in recv_array to gird(only copy the edge it received, not all four edges).
            copy_boundry(grid, rows_of_recv_data, cols_of_recv_data, recv_data, recv_left, recv_top, recv_right,
                         recv_bottom);
            // Perform the calculation for new_gird using grid and old_grid, and set the boundry value(Neumann boundaries) if process is in boundry.
            do_iteration(is_Neumann);
            // Write the output if it should
            if (t_out <= t) {
                grid_to_file(id, out_cnt, grid);
                out_cnt++;
                t_out += dt_out;
            }
            it++;
            delete[] request;
        }
    }
    // Calculate the program run time
    if (id == 0) {
        mpi_end = clock();
        double time = double(mpi_end - mpi_start) / CLOCKS_PER_SEC;
        cout << "Total time consumption: " << time << " s" << endl;
    }

    MPI_Type_free(&Datatype_left);
    MPI_Type_free(&Datatype_right);
    MPI_Type_free(&Datatype_top);
    MPI_Type_free(&Datatype_bottom);
    MPI_Type_free(&Datatype_left_recv);
    MPI_Type_free(&Datatype_right_recv);
    MPI_Type_free(&Datatype_top_recv);
    MPI_Type_free(&Datatype_bottom_recv);

    MPI_Finalize();
    
    delete[] send_data;
    delete[] recv_data;
}
