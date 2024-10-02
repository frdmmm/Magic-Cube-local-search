#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <cmath>
#include <limits>
#include <chrono>

using namespace std;


int magic_number(int n) {
    return (n * (n * n * n + 1)) / 2;
}

const int n = 5;
const int MAGIC_NUMBER = magic_number(n);


int sum_of_error(const vector<int>& state) {
    int total_error = 0;
    
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            // Garis sepanjang sumbu X (panjang)
            int line_sum = 0;
            for (int k = 0; k < n; ++k) {
                line_sum += state[i * n * n + j * n + k];
            }
            total_error += abs(line_sum - MAGIC_NUMBER);
            
            // Garis sepanjang sumbu Y (lebar)
            line_sum = 0;
            for (int k = 0; k < n; ++k) {
                line_sum += state[i * n * n + k * n + j];
            }
            total_error += abs(line_sum - MAGIC_NUMBER);
            
            // Garis sepanjang sumbu Z (tinggi)
            line_sum = 0;
            for (int k = 0; k < n; ++k) {
                line_sum += state[k * n * n + i * n + j];
            }
            total_error += abs(line_sum - MAGIC_NUMBER);
        }
    }

    // Garis diagonal di bidang XY
    for (int i = 0; i < n; ++i) {
        int line_sum = 0;
        for (int j = 0; j < n; ++j) {
            line_sum += state[i * n * n + j * n + j];
        }
        total_error += abs(line_sum - MAGIC_NUMBER);
        
        line_sum = 0;
        for (int j = 0; j < n; ++j) {
            line_sum += state[i * n * n + j * n + (n - 1 - j)];
        }
        total_error += abs(line_sum - MAGIC_NUMBER);
    }

    // Garis diagonal di bidang XZ
    for (int i = 0; i < n; ++i) {
        int line_sum = 0;
        for (int j = 0; j < n; ++j) {
            line_sum += state[j * n * n + i * n + j];
        }
        total_error += abs(line_sum - MAGIC_NUMBER);
        
        line_sum = 0;
        for (int j = 0; j < n; ++j) {
            line_sum += state[(n - 1 - j) * n * n + i * n + j];
        }
        total_error += abs(line_sum - MAGIC_NUMBER);
    }

    // Garis diagonal di bidang YZ
    for (int i = 0; i < n; ++i) {
        int line_sum = 0;
        for (int j = 0; j < n; ++j) {
            line_sum += state[j * n * n + j * n + i];
        }
        total_error += abs(line_sum - MAGIC_NUMBER);
        
        line_sum = 0;
        for (int j = 0; j < n; ++j) {
            line_sum += state[j * n * n + (n - 1 - j) * n + i];
        }
        total_error += abs(line_sum - MAGIC_NUMBER);
    }

    // Garis diagonal di ruang
    int line_sum = 0;
    for (int i = 0; i < n; ++i) {
        line_sum += state[i * n * n + i * n + i];
    }
    total_error += abs(line_sum - MAGIC_NUMBER);
    
    line_sum = 0;
    for (int i = 0; i < n; ++i) {
        line_sum += state[i * n * n + i * n + (n - 1 - i)];
    }
    total_error += abs(line_sum - MAGIC_NUMBER);
    
    line_sum = 0;
    for (int i = 0; i < n; ++i) {
        line_sum += state[(n - 1 - i) * n * n + i * n + i];
    }
    total_error += abs(line_sum - MAGIC_NUMBER);
    
    line_sum = 0;
    for (int i = 0; i < n; ++i) {
        line_sum += state[i * n * n + (n - 1 - i) * n + i];
    }
    total_error += abs(line_sum - MAGIC_NUMBER);

    return -total_error; // Total error, kita balik nilainya agar semakin mendekati 0 semakin baik
}

// Pertukaran dua elemen acak dalam state
void swap_random(vector<int>& state, mt19937& gen) {
    uniform_int_distribution<> dis(0, state.size() - 1);
    
    int i = dis(gen);
    int j = dis(gen);
    swap(state[i], state[j]);
}

// Simulated Annealing
vector<int> simulated_annealing(vector<int> state, double initial_temp, double cooling_rate, int max_iter) {
    vector<int> current_state = state;
    vector<int> best_state = current_state;
    double current_temp = initial_temp;

    // Gunakan seed berdasarkan waktu saat ini
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 gen(seed); 

    for (int i = 0; i < max_iter; ++i) {
        vector<int> new_state = current_state;
        swap_random(new_state, gen);
        int current_energy = sum_of_error(current_state);
        int new_energy = sum_of_error(new_state);
        
        if (new_energy > current_energy) {
            current_state = new_state;
            if (new_energy > sum_of_error(best_state)) {
                best_state = new_state;
            }
        } else {
            double acceptance_prob = exp((new_energy - current_energy) / current_temp);
            if ((double)rand() / RAND_MAX < acceptance_prob) {
                current_state = new_state;
            }
        }
        current_temp *= cooling_rate;
    }

    return best_state;
}

// Mengambil semua tetangga dari state
vector<vector<int>> get_neighbors(const vector<int>& state) {
    vector<vector<int>> neighbors;
    for (size_t i = 0; i < state.size(); ++i) {
        for (size_t j = i + 1; j < state.size(); ++j) {
            vector<int> new_state = state;
            swap(new_state[i], new_state[j]);
            neighbors.push_back(new_state);
        }
    }
    return neighbors;
}

// Steepest Ascent Hill Climbing
vector<int> steepest_ascent_hill_climbing(vector<int> state, int max_iter) {
    vector<int> current_state = state;
    int current_energy = sum_of_error(current_state);

    for (int i = 0; i < max_iter; ++i) {
        vector<vector<int>> neighbors = get_neighbors(current_state);
        vector<int> best_neighbor;
        int best_energy = current_energy;

        // Cek semua tetangga untuk menemukan yang terbaik
        for (const auto& neighbor : neighbors) {
            int neighbor_energy = sum_of_error(neighbor);
            if (neighbor_energy >= best_energy) {
                best_energy = neighbor_energy;
                best_neighbor = neighbor;
            }
        }

        // Jika tidak ada tetangga yang lebih baik, berhenti
        if (best_neighbor.empty()) {
            break;
        }

        // Update state saat ini ke tetangga terbaik
        current_state = best_neighbor;
        current_energy = best_energy;
    }

    return current_state;
}

// Random Restart Hill Climbing
vector<int> random_restart_hill_climbing(int max_restarts, int max_iter) {
    vector<int> best_state;
    int best_energy = numeric_limits<int>::min();

    // Gunakan seed berdasarkan waktu saat ini
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 gen(seed); 

    for (int i = 0; i < max_restarts; ++i) {
        vector<int> state(125);
        iota(state.begin(), state.end(), 1);
        shuffle(state.begin(), state.end(), gen); // Mengacak state awal

        vector<int> candidate_state = steepest_ascent_hill_climbing(state, max_iter);
        int candidate_energy = sum_of_error(candidate_state);

        // Jika state baru lebih baik, update best state
        if (candidate_energy > best_energy) {
            best_state = candidate_state;
            best_energy = candidate_energy;
        }
    }

    return best_state;
}

// Contoh eksekusi
int main() {
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 gen(seed); // Mengatur seed generator angka acak

    // Mengacak state awal
    vector<int> initial_state(125);
    iota(initial_state.begin(), initial_state.end(), 1);
    shuffle(initial_state.begin(), initial_state.end(), gen);

    vector<int> sa_solution = simulated_annealing(initial_state, 100000.0, 0.09, 1000000);
    cout << "SA Solution: ";
    for (int num : sa_solution) {
        cout << num << " ";
    }
    cout << "Energy: " << sum_of_error(sa_solution) << endl;

    vector<int> steepest = steepest_ascent_hill_climbing(initial_state, 100);
    cout << "Steepest Ascent Solution: ";
    for (int num : steepest) {
        cout << num << " ";
    }
    cout << "Energy: " << sum_of_error(steepest) << endl;

    return 0;
}
