#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <chrono>
#include <boost/multiprecision/cpp_int.hpp>

using namespace boost::multiprecision;

// 素数生成 (エラトステネスの篩)
std::vector<int> generate_primes(int max_val) {
    std::vector<bool> is_p(max_val + 1, true);
    std::vector<int> primes;
    for (int p = 2; p <= max_val; p++) {
        if (is_p[p]) {
            primes.push_back(p);
            for (int i = p * 2; i <= max_val; i += p) is_p[i] = false;
        }
    }
    return primes;
}

cpp_int gcd(cpp_int a, cpp_int b) {
    while (b != 0) { cpp_int t = b; b = a % b; a = t; }
    return a;
}

// 取得した関係式 (Relation) を保持する構造体
struct Relation {
    cpp_int z;
    cpp_int v;
    std::vector<int> exponents;     // 素因数分解の指数
    std::vector<uint8_t> gf2_vector; // 指数のパリティ (0 or 1)
};

// 複数多項式ではなく単一の基本二次篩法 (Basic Quadratic Sieve)
void run_qs(cpp_int N) {
    std::cout << "[*] Starting Quadratic Sieve for N = " << N << std::endl;
    
    // 【ステップ1】素数ベース (Factor Base) の生成
    // 巨大な数の場合はここを数万〜数百万にします
    int B = 10000; 
    std::vector<int> primes = generate_primes(B);
    int C = primes.size() + 1; // +1は負号 (-1) のため
    
    std::cout << "[Phase 1] Factor base size: " << C - 1 << " primes." << std::endl;
    
    cpp_int S = sqrt(N);
    std::vector<Relation> relations;
    
    long long i = 0;
    long long sign = 1;
    
    int target_relations = C + 20; // 未知数より少し多めに関係式を集める
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "[Phase 2] Sieving for relations (target " << target_relations << ")..." << std::endl;
    
    // 【ステップ2】関係式 (Relations) の収集 (シーブ処理の代替)
    while (relations.size() < target_relations) {
        cpp_int z = S + (i * sign);
        if (sign > 0) i++;
        sign = -sign;
        
        if (z <= 0) continue;
        
        cpp_int v = z * z - N;
        if (v == 0) {
            std::cout << "Perfect square found! " << z << std::endl;
            return;
        }
        
        cpp_int abs_v = v < 0 ? -v : v;
        
        std::vector<int> exps(C, 0);
        std::vector<uint8_t> gf2(C, 0);
        
        // 負の場合はインデックス0に記録
        if (v < 0) {
            exps[0] = 1;
            gf2[0] = 1;
        }
        
        // Trial division over Factor Base
        for (size_t p_idx = 0; p_idx < primes.size(); p_idx++) {
            int p = primes[p_idx];
            while (abs_v % p == 0) {
                abs_v /= p;
                exps[p_idx + 1]++;
                gf2[p_idx + 1] ^= 1;
            }
        }
        
        // 完全にFactor Base内の素数だけで分解できた場合 (B-smooth)
        if (abs_v == 1) {
            Relation r;
            r.z = z;
            r.v = v;
            r.exponents = exps;
            r.gf2_vector = gf2;
            relations.push_back(r);
            if (relations.size() % 100 == 0) {
                std::cout << "\r  Found " << relations.size() << " relations" << std::flush;
            }
        }
    }
    std::cout << std::endl;
    
    // 【ステップ3】GF(2) 上での線形代数 (ガウスの消去法)
    std::cout << "[Phase 3] Building GF(2) matrix and finding null space..." << std::endl;
    
    int R = relations.size();
    std::vector<std::vector<uint8_t>> matrix(R, std::vector<uint8_t>(C + R, 0));
    for (int r = 0; r < R; r++) {
        for (int c = 0; c < C; c++) matrix[r][c] = relations[r].gf2_vector[c];
        matrix[r][C + r] = 1; // どの関係式を使ったか追跡するための単位行列
    }
    
    int pivot_row = 0;
    for (int c = 0; c < C && pivot_row < R; c++) {
        int r = pivot_row;
        while (r < R && matrix[r][c] == 0) r++;
        if (r == R) continue;
        
        std::swap(matrix[pivot_row], matrix[r]);
        
        for (int j = 0; j < R; j++) {
            if (j != pivot_row && matrix[j][c] == 1) {
                for (int k = c; k < C + R; k++) {
                    matrix[j][k] ^= matrix[pivot_row][k]; // GF(2)での足し算(XOR)
                }
            }
        }
        pivot_row++;
    }
    
    // 【ステップ4】平方根の計算とGCDによる因数発見
    std::cout << "[Phase 4 & 5] Computing Square Roots and GCD..." << std::endl;
    
    bool found = false;
    for (int r = pivot_row; r < R; r++) {
        std::vector<int> used_relations;
        for (int j = 0; j < R; j++) {
            if (matrix[r][C + j] == 1) used_relations.push_back(j);
        }
        if (used_relations.empty()) continue;
        
        cpp_int X = 1;
        std::vector<int> total_exps(C, 0);
        
        for (int idx : used_relations) {
            X = (X * relations[idx].z) % N;
            for (int c = 0; c < C; c++) {
                total_exps[c] += relations[idx].exponents[c];
            }
        }
        
        cpp_int Y = 1;
        for (int c = 1; c < C; c++) {
            int half_exp = total_exps[c] / 2;
            for (int k = 0; k < half_exp; k++) {
                Y = (Y * primes[c - 1]) % N;
            }
        }
        
        cpp_int diff = X > Y ? X - Y : Y - X;
        cpp_int factor = gcd(diff, N);
        
        if (factor > 1 && factor < N) {
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end_time - start_time;
            
            std::cout << "\n================================================================" << std::endl;
            std::cout << "[SUCCESS] Non-trivial factor found using Quadratic Sieve!" << std::endl;
            std::cout << "  Factor 1: " << factor << std::endl;
            std::cout << "  Factor 2: " << N / factor << std::endl;
            std::cout << "  Time    : " << elapsed.count() << " seconds." << std::endl;
            std::cout << "================================================================" << std::endl;
            found = true;
            break;
        }
    }
    
    if (!found) {
        std::cout << "[FAILURE] Could not find a non-trivial factor. Need a larger factor base." << std::endl;
    }
}

int main() {
    std::cout << "================================================================" << std::endl;
    std::cout << "  Quadratic Sieve (Predecessor & Core Architecture of GNFS)     " << std::endl;
    std::cout << "================================================================" << std::endl;
    std::cout << "素因数分解したい数値を入力してください: ";
    
    std::string input_str;
    if (!(std::cin >> input_str)) return 0;
    
    cpp_int N;
    try {
        N = cpp_int(input_str);
    } catch (...) {
        std::cerr << "無効な数値形式です。" << std::endl;
        return 1;
    }
    
    run_qs(N);
    return 0;
}
