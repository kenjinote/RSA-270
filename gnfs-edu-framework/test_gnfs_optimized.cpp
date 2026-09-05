#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>

// gnfs_optimized.cpp 内の main 関数との衝突を避けるためのマクロ
#define main original_main
#include "gnfs_optimized.cpp"
#undef main

/**
 * @brief シーブ処理が正しく動作し、対数近似値の加算と閾値判定が正確に行われているかを検証するテスト
 */
void test_sieve_results() {
    std::cout << "[Test] Checking Sieve results correctness..." << std::endl;
    uint32_t range = 10000; // テスト用の範囲
    OptimizedLatticeSieve sieve(range);
    sieve.generate_factor_bases(100);
    
    // 最適化されたシーブ処理を実行
    std::vector<Relation> relations = sieve.execute_sieve();
    
    // 期待される結果を愚直な方法（ナイーブな計算）で導出する
    std::vector<uint32_t> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31};
    std::vector<uint8_t> log_p(primes.size());
    for (size_t i = 0; i < primes.size(); ++i) {
        log_p[i] = static_cast<uint8_t>(std::round(std::log2(primes[i]) * 10.0));
    }
    
    std::vector<uint32_t> expected_a;
    for (uint32_t a = 0; a < range; ++a) {
        uint32_t sum = 0;
        // 各素数について、割り切れる場合にのみ対数近似値を加算（元のロジック通り素数べきは考慮しない）
        for (size_t i = 0; i < primes.size(); ++i) {
            if (a % primes[i] == 0) {
                sum += log_p[i];
            }
        }
        // 閾値は gnfs_optimized.cpp 内で 50 に設定されている
        if (sum >= 50) {
            expected_a.push_back(a);
        }
    }
    
    // マルチスレッド(OpenMP)により実行結果の順序が保証されないため、a の値でソートする
    std::vector<uint32_t> actual_a;
    for (const auto& r : relations) {
        actual_a.push_back(r.a);
    }
    std::sort(actual_a.begin(), actual_a.end());
    
    // 両者の要素数が一致するか検証
    assert(expected_a.size() == actual_a.size() && "Result size mismatch!");
    
    // 両者の要素がすべて一致するか検証
    for (size_t i = 0; i < expected_a.size(); ++i) {
        if (expected_a[i] != actual_a[i]) {
            std::cerr << "Mismatch at index " << i << ": expected " << expected_a[i] 
                      << ", got " << actual_a[i] << std::endl;
            assert(false);
        }
    }
    std::cout << " => OK! Found " << actual_a.size() << " valid relations perfectly matched the expected logical result." << std::endl;
}

/**
 * @brief L1キャッシュサイズ(32KB)を超える範囲に対するブロック分割処理が正しく動作するかを検証するテスト
 */
void test_cache_blocking_logic() {
    std::cout << "[Test] Checking L1 Cache Blocking logic over multiple blocks..." << std::endl;
    // L1_CACHE_SIZE(32768) を跨ぐ範囲
    uint32_t range = L1_CACHE_SIZE * 3 + 12345; 
    OptimizedLatticeSieve sieve(range);
    sieve.generate_factor_bases(100);
    
    std::vector<Relation> relations = sieve.execute_sieve();
    
    // 最低限、結果が空でないこと（確実に何か見つかる範囲）を確認し、クラッシュせずに終了したことを確認
    assert(!relations.empty());
    std::cout << " => OK! Processed range of " << range << " spanning multiple L1 cache blocks successfully." << std::endl;
}

int main() {
    std::cout << "================================================================" << std::endl;
    std::cout << "  Running Unit Tests for GNFS Optimized Sieve Engine            " << std::endl;
    std::cout << "================================================================" << std::endl;
    
    test_sieve_results();
    std::cout << "----------------------------------------------------------------" << std::endl;
    test_cache_blocking_logic();
    
    std::cout << "================================================================" << std::endl;
    std::cout << "  All tests passed successfully!                                " << std::endl;
    std::cout << "================================================================" << std::endl;
    return 0;
}
