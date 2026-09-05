#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/miller_rabin.hpp>
#include <boost/random.hpp>

using namespace boost::multiprecision;
using namespace boost::random;

// 最大公約数を計算する関数
cpp_int gcd(cpp_int a, cpp_int b) {
    while (b != 0) {
        cpp_int temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

// ポラード・ロー素因数分解アルゴリズム (Pollard's rho algorithm)
// 数十桁（60〜100bit程度）までの合成数の素因数を高速に見つけます。
cpp_int pollard_rho(const cpp_int& n) {
    if (n % 2 == 0) return 2;
    
    cpp_int x, y, d, c;
    c = 1;
    
    while (true) {
        x = 2;
        y = 2;
        d = 1;
        
        // 疑似乱数関数 f(x) = (x^2 + c) mod n
        auto f = [&](const cpp_int& val) {
            return (val * val + c) % n;
        };
        
        uint32_t step = 0;
        while (d == 1) {
            x = f(x);
            y = f(f(y));
            cpp_int diff = (x > y) ? (x - y) : (y - x);
            d = gcd(diff, n);
            
            if (d == n) {
                break; // 失敗：非自明な約数が見つかる前にループに入ってしまった場合
            }
            
            step++;
            if (step % 100000 == 0) {
                std::cout << "  ... searching (step " << step << ")" << std::endl;
            }
        }
        
        if (d != 1 && d != n) {
            return d; // 非自明な素因数を発見
        }
        
        c++; // 定数cを変更して再挑戦
    }
}

// 再帰的に完全に素因数分解を行う関数
void factorize(const cpp_int& n, std::vector<cpp_int>& factors) {
    if (n == 1) return;
    
    // Miller-Rabin 素数判定テストを用いて、すでに素数かどうかを確認
    mt11213b base_gen(clock());
    independent_bits_engine<mt11213b, 256, cpp_int> gen(base_gen);
    
    if (miller_rabin_test(n, 25, gen)) {
        factors.push_back(n);
        return;
    }
    
    // 合成数の場合は素因数を一つ見つける
    cpp_int divisor = pollard_rho(n);
    
    // 見つかった約数と、残りの部分について再帰的に素因数分解を行う
    factorize(divisor, factors);
    factorize(n / divisor, factors);
}

int main() {
    std::cout << "================================================================" << std::endl;
    std::cout << "  Practical Integer Factorizer (Pollard's rho + Miller-Rabin)   " << std::endl;
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
    
    std::cout << "\n[*] Target N = " << N << std::endl;
    std::cout << "[*] Starting factorization...\n" << std::endl;
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::vector<cpp_int> factors;
    factorize(N, factors);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;
    
    std::cout << "\n[SUCCESS] Factorization completed in " << elapsed.count() << " seconds." << std::endl;
    std::cout << "Factors of " << N << " are:" << std::endl;
    
    for (size_t i = 0; i < factors.size(); ++i) {
        std::cout << "  Prime " << i+1 << ": " << factors[i] << std::endl;
    }
    
    // 検算
    cpp_int verify = 1;
    for (const auto& f : factors) verify *= f;
    
    if (verify == N) {
        std::cout << "\n[Verify] OK: The product of factors exactly matches N." << std::endl;
    } else {
        std::cout << "\n[Verify] FAILED." << std::endl;
    }
    
    std::cout << "================================================================" << std::endl;
    return 0;
}
