#include <map>
#include <boost/optional/optional_io.hpp>
#include <chrono>

#include "btree.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

void fail()
{
    std::cout << "FAILED\n";
    exit(-1);
}

class TimeCounter {
public:
    TimeCounter()
        : begin(high_resolution_clock::now())
    {}

    ~TimeCounter()
    {
        auto duration = duration_cast<milliseconds>(
                             high_resolution_clock::now() - begin).count();
        std::cout << "Execution wall time, ms: " << duration << '\n';
    }

private:
    high_resolution_clock::time_point begin;
};

template <class Key, class Value>
void print_tree_info(const MyBTree::BTree<Key, Value>& tree)
{
    std::cout << "io count: " << tree.get_io_count() << '\t';
    std::cout << "node binary size: " << tree.get_node_binary_size() << '\t';
}

int calc_value(int key)
{
    return key * 10;
}

// #define TREE_DUMP
template <class T>
void tree_dump(T& tree)
{
#ifdef TREE_DUMP
    std::cout << '\n';
    tree.dump();
#endif
}

std::vector<int> test_create(
        const char* filename,
        size_t min_degree,
        const std::vector<int>& keys_to_add,
        const std::vector<int>& keys_to_remove)
{
    TimeCounter time_counter;

    std::cout << "CREATE ";
    std::cout << std::setw(13) << std::left
              << "min_degree: " << min_degree << '\t';

    MyBTree::BTree<int, int> tree(min_degree, filename);

    for (const auto& key : keys_to_add) {
        tree.insert(key, calc_value(key));
    }

    tree_dump(tree);

    auto residual_keys(keys_to_add);
    for (const auto& key : keys_to_remove) {
#ifdef TREE_DUMP
        std::cout << "REMOVING " << key << '\n';
#endif
        tree.remove(key);

        tree_dump(tree);

        residual_keys.erase(
                std::remove(residual_keys.begin(), residual_keys.end(), key),
                residual_keys.end()
        );

    }

    print_tree_info(tree);

    return residual_keys;
}

void test_read(
        const char* filename,
        const std::vector<int>& existing_keys,
        const std::vector<int>& not_existing_keys)
{
    TimeCounter time_counter;

    std::cout << std::setw(21) << std::left << "READ\t";

    MyBTree::BTree<int, int> tree(filename);

    tree_dump(tree);

    if (tree.size() != existing_keys.size()) {
        std::cout << "\ntree size: " << tree.size() << '\n'
                  << "expected: " << existing_keys.size() << '\n';
        fail();
    }

    for (const auto& key : existing_keys) {
        const auto value = tree.find(key);
        auto expected_value = calc_value(key);
        if (value != expected_value) {
            std::cout << "\nfind() failed. key: " << key
                      << ", value = " << value
                      << ", expected value = " << expected_value
                      << '\n';
            fail();
        }
    }

    for (const auto& key : not_existing_keys) {
        const auto value = tree.find(key);
        if (value) {
            std::cout << "found removed key: " << key
                      << ". Found value: " << *value << '\n';
            fail();
        }
    }

    print_tree_info(tree);
}

void test(
        size_t min_degree,
        const std::vector<int>& keys_to_add,
        const std::vector<int>& keys_to_remove = {})
{
    const char* filename = "db.dat";

    auto residual_keys = test_create(
            filename, min_degree, keys_to_add, keys_to_remove);
    test_read(filename, residual_keys, keys_to_remove);

    std::cout << '\n';
}

std::vector<int> gen_vector(int size)
{
    std::vector<int> v(size);
    std::iota(v.begin(), v.end(), 0);
    return v;
}

std::vector<int> gen_vector(int begin, int size)
{
    std::vector<int> v(size);
    std::iota(v.begin(), v.end(), begin);
    return v;
}

std::vector<int> reverse(std::vector<int> v) {
    std::reverse(v.begin(), v.end());
    return v;
}

void performance_test_impl(
        const std::vector<int>& keys_to_add,
        const std::vector<int>& keys_to_delete)
{
    test(2, keys_to_add, keys_to_delete);
    test(16, keys_to_add, keys_to_delete);
    test(32, keys_to_add, keys_to_delete);
    test(512, keys_to_add, keys_to_delete);
    test(2048, keys_to_add, keys_to_delete);
}

void performance_test()
{
    std::cout << "PERFORMANCE TEST\n";

    auto keys_to_add = gen_vector(30000);
    performance_test_impl(keys_to_add, {});

    std::cout << "PERFORMANCE TEST with deletion\n";
    performance_test_impl(keys_to_add, gen_vector(13000, 10000));
}

void test_addition(int min_degree)
{
    test(min_degree, {});
    test(min_degree, {1});
    test(min_degree, {1,2});
    test(min_degree, {1, 2, 3});
    test(min_degree, {1, 2, 3, 4});
    test(min_degree, {2, 1});
    test(min_degree, {4, 3, 2, 1});
    test(min_degree, {4, 3, 1, 2, 0, 100, 150, 99});

    test(min_degree, gen_vector(6));
    test(min_degree, gen_vector(1000));

    test(min_degree, reverse(gen_vector(10)));
    test(min_degree, reverse(gen_vector(1000)));
}

void test_deletion_2_case(int min_degree)
{
    // 2-a
    test(min_degree, {7,8,9,1}, {});
    test(min_degree, {7,8,9,1}, {8});

    // 2-b
    test(min_degree, {7,8,9,10}, {});
    test(min_degree, {7,8,9,10}, {8});

    // 2-c without root elimination
    test(min_degree, {5,6,9,10,7,8}, {});
    test(min_degree, {5,6,9,10,7,8}, {8});
    test(min_degree, {5,6,9,10,7,8}, {8,6});

    test(min_degree, {5,6,9,10,7,8}, {8,9});

    // 2-c with root elimination
    test(min_degree, {7,8,9,10}, {});
    test(min_degree, {7,8,9,10}, {10, 8});
}

void test_deletion_3_case(int min_degree)
{
    test(min_degree, {7,8,9,1}, {});
    test(min_degree, {7,8,9,1}, {1});
    test(min_degree, {7,8,9,1}, {7});

    test(min_degree, {7,8,9,10}, {});
    test(min_degree, {7,8,9,10}, {9});
    test(min_degree, {7,8,9,10}, {10});

    // 3-a

    test(min_degree, {7,8,9,1}, {});
    test(min_degree, {7,8,9,1}, {9});

    test(min_degree, {7,8,9,10}, {});
    test(min_degree, {7,8,9,10}, {7});

    // 3-b

    test(min_degree, {7,8,9,10}, {7,10});
    test(min_degree, {7,8,9,10}, {7,8});
}

void test_deletion(int min_degree)
{
    test(min_degree, {}, {1});
    test(min_degree, {1}, {1});
    test(min_degree, {1, 2}, {1});
    test(min_degree, {1, 2, 3}, {1, 1, 0});

    test_deletion_2_case(min_degree);
    test_deletion_3_case(min_degree);

    test(min_degree, gen_vector(10), gen_vector(1));
    test(min_degree, gen_vector(10), gen_vector(2));
    test(min_degree, gen_vector(10), gen_vector(9));
    test(min_degree, gen_vector(10), gen_vector(10));
    test(min_degree, gen_vector(100), reverse(gen_vector(100)));
    test(min_degree, reverse(gen_vector(10)), reverse(gen_vector(11)));
    test(min_degree, reverse(gen_vector(10)), reverse(gen_vector(10)));
    test(min_degree, reverse(gen_vector(12)), reverse(gen_vector(12)));
    test(min_degree, gen_vector(1000), reverse(gen_vector(500)));
    test(min_degree, reverse(gen_vector(1000)), gen_vector(500));
    test(min_degree, reverse(gen_vector(1000)), reverse(gen_vector(500)));
    test(min_degree, gen_vector(1000), gen_vector(100, 300));
}

int main()
{
    for (int i = 2; i < 20; ++i) {
        test_addition(i);
        test_deletion(i);
    }
    performance_test();

    return 0;
}
