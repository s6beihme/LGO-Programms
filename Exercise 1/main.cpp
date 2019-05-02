#include <iostream>
#include <memory>
#include <vector>
#include <string>
#include <string>
#include <fstream>

//basic implementation of the Fourier-Motzkin Elimination
//this is extremely inefficient, especially regarding memory!!!

//We put everything into 1 file because the the program isnt so big

typedef std::unique_ptr<std::unique_ptr<double[]>[]> Matrix;

void print_matrix(Matrix& A, int m, int n) {
	std::cout << "\n";
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			std::cout << A[i][j] << " ";
		}
		std::cout << "\n";
	}
}


//class representing an LP of the form of Ax<=b
class LP {
public:
	//constructor allocating memory
	LP(int _m, int _n);

	//constructor reading LP from file filename
	LP(std::string filename);

	void print();
	//dimensions of matrix
	int m,n;

	std::unique_ptr<double[]> c;
	Matrix A;
	std::unique_ptr<double[]> b;
};

LP::LP(int _m, int _n)  {
	m = _m;
	n = _n;
	A = std::make_unique<std::unique_ptr<double[]>[]>(m);
	for (int i = 0; i < m; i++) {
		A[i] = std::make_unique<double[]>(n);
	}
	b = std::make_unique<double[]>(m);
	c = std::make_unique<double[]>(n);
}

LP::LP(std::string filename) {
	//open file
	std::ifstream myfile;
	myfile.open(filename);

	if (!myfile.is_open()) { 
		std::cout << "FILE DIDNT OPEN!\n";
		exit(0);
	}

	//read in dimensions
	myfile >> m >> n;

	//allocate memory for lp
	A = std::make_unique<std::unique_ptr<double[]>[]>(m);
	for (int i = 0; i < m; i++) {
		A[i] = std::make_unique<double[]>(n);
	}
	b = std::make_unique<double[]>(m);
	c = std::make_unique<double[]>(n);

	//read in LP
	for (int i = 0; i < n; i++) myfile >> c[i];
	for (int i = 0; i < m; i++) myfile >> b[i];
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			myfile >> A[i][j];
		}
	}
}

void LP::print() {
	std::cout << "\nc:\n";
	for (int i = 0; i < n; i++) std::cout << c[i] << ", ";
	std::cout << "\nA:\n";
	print_matrix(A, m, n);
	std::cout << "\nb:\n";
	for (int i = 0; i < m; i++) std::cout << b[i] << ", ";
	std::cout << std::endl;
}


//the solution vectors have to have memory allocated when passing them to this function
//this function is to be called on the original LP with History=Im and computes a solution
//or a certificate recursively
bool Fourier_Motzkin_Elimination_step(std::vector<double> &solution_infeasible,
									  std::vector<double> &solution_feasible,
									  LP& old_LP,
									  Matrix &History,
									  int original_row_count,
									  int original_column_count) {

	//make vectors U, L and N like in the script to store the indices
	//of positive, negative and zero-rows respectively
	double temp;
	std::vector<int> U;
	std::vector<int> L;
	std::vector<int> N;
	int positive_row_count = 0;
	int negative_row_count = 0;
	int zero_row_count = 0;

	for (int i = 0; i < old_LP.m; i++) {
		if (old_LP.A[i][0] > 0) {
			temp = old_LP.A[i][0];
			positive_row_count++;
			U.push_back(i);
			//divide row and history array by A[i][0] 
			for (int j = 0; j < old_LP.n; j++) {
				old_LP.A[i][j] /= temp;
			}
			old_LP.b[i]/= temp;
			for (int j = 0; j < original_row_count; j++) {
				History[i][j] /= temp;
			}
		}
		if (old_LP.A[i][0] < 0) {
			temp=old_LP.A[i][0];
			negative_row_count++;
			L.push_back(i);
			//divide row and history array by A[i][0] 
			for (int j = 0; j < old_LP.n; j++) {
				old_LP.A[i][j] = old_LP.A[i][j] /(-temp);
			}
			old_LP.b[i]= old_LP.b[i] /(-temp);
			for (int j = 0; j < original_row_count; j++) {
				History[i][j] = History[i][j] /(-temp);
			}
		}
		if (old_LP.A[i][0] == 0) {
			zero_row_count++;
			N.push_back(i);
		}
	}

	//create new LP that will be equivalent (regarding feasibility) to the old one
	LP new_LP(positive_row_count*negative_row_count + zero_row_count, old_LP.n - 1); 

	Matrix new_History = std::make_unique<std::unique_ptr<double[]>[]> (new_LP.m);
	for (int i = 0; i < new_LP.m; i++) {
		new_History[i] = std::make_unique<double[]>(original_row_count);
	}

	//create new inequalities
	int current_row = 0;
	for (int u : U) {
		for (int l : L) {
			for (int i = 0; i < new_LP.n; i++) {
				new_LP.A[current_row][i] = old_LP.A[u][i + 1] + old_LP.A[l][i + 1];
			}
			new_LP.b[current_row] = old_LP.b[u] + old_LP.b[l];
			for (int i = 0; i < original_row_count; i++) {
				new_History[current_row][i] = History[u][i] + History[l][i];
			}
			current_row++;
		}
	}
	for (int n : N) {
		for (int i = 0; i < new_LP.n; i++) {
			new_LP.A[current_row][i] = old_LP.A[n][i + 1];
		}
		new_LP.b[current_row] = old_LP.b[n];
		for (int i = 0; i < original_row_count; i++) {
			new_History[current_row][i] = History[n][i];
		}
		current_row++;
	}

	if (new_LP.n == 0) {
		//check if LP is infeasible. if so, retur the coefficients
		//of the rows of the original LP that led to a contradiction
		for (int i = 0; i < new_LP.m; i++) {
			if (new_LP.b[i] < 0) {
				for (int j = 0; j < original_row_count; j++) solution_infeasible[j]=(new_History[i][j]);
				return false;
			}
		}
		return true;
	}

	if (Fourier_Motzkin_Elimination_step(solution_infeasible, solution_feasible, new_LP, new_History, original_row_count, original_column_count) == false) {
		return false;
	}
	else {
		//THIS DOESNT WORK PROPERLY FOR SOME REASON
		//find min{a_l(^T)*x-b_l | l in L}
		double lower_bound;
		bool lower_bound_initialized=false;
		double upper_bound;
		bool upper_bound_initialized=false;
		for (int l : L) {
			double sum = 0;
			//compute vector product (a_l)^T*solution_feasible
			for (int i = 0; i < new_LP.n - 1; i++) {
				sum += new_LP.A[l][new_LP.n-i-1] * solution_feasible[original_column_count - i-1];
			}

			if (lower_bound_initialized == false) {
				lower_bound = sum - new_LP.b[l];
				lower_bound_initialized = true;
			}
			else {
				if (sum - new_LP.b[l] > lower_bound) lower_bound = sum - new_LP.b[l];
			}
		}
		for (int u : U) {
			double sum = 0;
			//compute vector product (a_l)^T*solution_feasible
			for (int i = 0; i < new_LP.n - 1; i++) {
				sum += new_LP.A[u][new_LP.n - i-1] * solution_feasible[original_column_count - i-1];
			}

			if (upper_bound_initialized == false) {
				upper_bound = new_LP.b[u] - sum;
				upper_bound_initialized = true;
			}
			else {
				if (new_LP.b[u] - sum < lower_bound) lower_bound = new_LP.b[u] - sum;
			}
		}

		if (lower_bound_initialized == true) solution_feasible[new_LP.n - 1] = lower_bound;
		else if (upper_bound_initialized == true) solution_feasible[new_LP.n - 1] = upper_bound;
		else solution_feasible[new_LP.n - 1] = 0;
		return true;
	}
}

int main(int argc, char* argv[]) {
	//read in the LP
	std::string filename(argv[1]);
	LP test(filename);

	test.print();

	//set up solutions and history
	std::vector<double> solution_feas (test.n);
	std::vector<double> solution_infeas (test.m);
	Matrix history = std::make_unique<std::unique_ptr<double[]>[]>(test.m);
	for (int i = 0; i < test.m; i++) {
		history[i] = std::make_unique<double[]>(test.m);
	}

	//initialize history
	for (int i = 0; i < test.m; i++) {
		history[i][i] = 1;
	}

	if (Fourier_Motzkin_Elimination_step(solution_infeas, solution_feas, test, history, test.m, test.n) == 0) {
		std::cout << "\nEmpty!\n";
		for (int i = 0; i < test.m; i++) std::cout << solution_infeas[i] << ", ";
	}
	else {
		std::cout << "\n";
		for (int i = 0; i < test.n; i++) std::cout << solution_feas[i] << ", ";
	}
	return 0;
}

