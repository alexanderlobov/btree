#pragma once

#include "log.h"
#include "serialization.h"

#include <boost/optional.hpp>

#include <fstream>
#include <vector>
#include <iostream>
#include <iomanip>

// Implementation according to
// "Introduction to Algorithms" by Cormen, Leiserson, Rivest, Stein

namespace MyBTree
{

// tested only for Key=int, Value=int
// Need some additional support for other types, e.g
// custom comparison functions and complex types serialization

template <class Key, class Value>
class BTree {
public:
    // creates new file
    BTree(size_t minimum_degree, const char* filename);

    // reads existing file
    explicit BTree(const char* filename);

    BTree(const BTree&) = delete;
    BTree& operator=(const BTree&) = delete;

    ~BTree();

    boost::optional<Value> find(const Key& key);
    void insert(const Key& key, const Value& value);
    void remove(const Key& key);
    size_t size();

    // methods for testing

    void dump()
    {
        std::cout << "total_number_of_keys = " << total_number_of_keys << '\n';
        dump_node(root_node, root_index, 0);
    }

    int get_io_count() const
    {
        return io_operations_count;
    }

    int get_node_binary_size() const
    {
        return node_binary_size();
    }

private:

    using NodeIndex = size_t;

    static size_t max_keys_size(size_t minimum_degree)
    {
        return max_children_size(minimum_degree) - 1;
    }

    static size_t max_children_size(size_t minimum_degree)
    {
        return 2 * minimum_degree;
    }

    struct Node {
        Node()
        {}

        Node(size_t minimum_degree)
            : keys(max_keys_size(minimum_degree))
            , values(max_keys_size(minimum_degree))
            , children(max_children_size(minimum_degree))
        {}

        void resize(size_t minimum_degree)
        {
            keys.resize(max_keys_size(minimum_degree));
            values.resize(max_keys_size(minimum_degree));
            children.resize(max_children_size(minimum_degree));
        }

        bool is_full() const
        {
            return number_of_keys == keys.size();
        }

        void dump_node() const
        {
            LOG_FIELD(number_of_keys);
            LOG_FIELD(is_leaf);
            print_vector(keys);
            print_vector(values);
            print_vector(children);
        }

        int number_of_keys = 0;
        bool is_leaf = true;
        std::vector<Key> keys;
        std::vector<Value> values;
        std::vector<NodeIndex> children;
    };

    size_t node_offset(NodeIndex index) const
    {
        return metadata_size() + index * node_binary_size();
    }

    size_t metadata_size() const
    {
        return sizeof(minimum_degree) +
               sizeof(num_nodes) +
               sizeof(root_index) +
               sizeof(total_number_of_keys);
    }

    size_t node_binary_size() const
    {
        return sizeof(root_node.number_of_keys) +
               sizeof(root_node.is_leaf) +
               binary_size(root_node.keys) +
               binary_size(root_node.values) +
               binary_size(root_node.children);
    }

    void write_node(const Node& node, NodeIndex index)
    {
        LOG_DEBUG("writing node with index " << index);

        file.seekp(node_offset(index));
        write(file, node.number_of_keys);
        write(file, node.is_leaf);
        write(file, node.keys);
        write(file, node.values);
        write(file, node.children);

        ++io_operations_count;

#ifndef NDEBUG
        node.dump_node();
#endif
    }

    void read_node(Node& node, NodeIndex index)
    {
        LOG_DEBUG("reading node with index " << index);

        file.seekg(node_offset(index));
        read(file, node.number_of_keys);
        read(file, node.is_leaf);
        read(file, node.keys);
        read(file, node.values);
        read(file, node.children);

        ++io_operations_count;

#ifndef NDEBUG
        node.dump_node();
#endif
    }

    void write_metadata()
    {
        file.seekp(0);
        LOG_DEBUG("write_metadata()");

        write(file, minimum_degree);
        write(file, num_nodes);
        write(file, root_index);
        write(file, total_number_of_keys);

        dump_metadata();
    }

    void read_metadata()
    {
        LOG_DEBUG("read_metadata()")

        read(file, minimum_degree);
        read(file, num_nodes);
        read(file, root_index);
        read(file, total_number_of_keys);

        dump_metadata();
    }

    void dump_metadata() const
    {
        LOG_FIELD(minimum_degree);
        LOG_FIELD(num_nodes);
        LOG_FIELD(root_index);
        LOG_FIELD(total_number_of_keys);
    }

    void split_child(Node& node, NodeIndex node_index,
                     size_t child_pos, Node& child, NodeIndex child_index,
                     Node& new_child)

    {
        LOG_DEBUG("split_child()")
        assert(!node.is_full());

        // Node new_child(minimum_degree);
        new_child.is_leaf = child.is_leaf;
        new_child.number_of_keys = minimum_degree - 1;

        NodeIndex new_child_index = num_nodes++;

        // copy second half of keys to new child
        for (size_t j = 0; j < new_child.number_of_keys; ++j) {
            new_child.keys[j] = child.keys[j + minimum_degree];
            new_child.values[j] = child.values[j + minimum_degree];
        }

        // copy second half of children to new child
        if (! child.is_leaf) {
            for (size_t j = 0; j < new_child.number_of_keys + 1; ++j) {
                new_child.children[j] = child.children[j + minimum_degree];
            }
        }

        child.number_of_keys = minimum_degree - 1;

        for (size_t j = node.number_of_keys; j > child_pos; --j) {
            node.children[j+1] = node.children[j];
        }
        node.children[child_pos + 1] = new_child_index;

        for (int j = node.number_of_keys - 1; j >= static_cast<int>(child_pos); --j) {
            node.keys[j + 1] = node.keys[j];
            node.values[j + 1] = node.values[j];
        }
        node.keys[child_pos] = child.keys[minimum_degree - 1];
        node.values[child_pos] = child.values[minimum_degree - 1];
        ++node.number_of_keys;

        write_node(node, node_index);
        write_node(child, child_index);
        write_node(new_child, new_child_index);
    }

    void insert_nonfull(Node& node, NodeIndex node_index,
                   const Key& key, const Value& value)
    {
        LOG_DEBUG("insert_nonfull()")
        assert(!node.is_full());
        int i = node.number_of_keys - 1;
        if (node.is_leaf) {
            while (i >= 0 && key < node.keys[i]) {
                node.keys[i+1] = node.keys[i];
                node.values[i+1] = node.values[i];
                --i;
            }
            node.keys[i+1] = key;
            node.values[i+1] = value;
            ++node.number_of_keys;
            ++total_number_of_keys;
            write_node(node, node_index);
        } else {
            while (i >= 0 && key < node.keys[i]) {
                --i;
            }
            ++i;
            Node child(minimum_degree);
            NodeIndex child_index = node.children[i];
            read_node(child, child_index);

            if (child.is_full()) {
                Node new_child(minimum_degree);
                split_child(node, node_index, i, child, child_index, new_child);
                if (key > node.keys[i]) {
                    insert_nonfull(new_child, node.children[i + 1], key, value);
                } else {
                    insert_nonfull(child, node.children[i], key, value);
                }
            } else {
                insert_nonfull(child, node.children[i], key, value);
            }
        }
    }

    boost::optional<Value> find(const Node& node, const Key& key)
    {
        LOG_DEBUG("find(" << key << ")");
        size_t i = 0;
        while (i < node.number_of_keys && key > node.keys[i]) {
            ++i;
        }
        if (i < node.number_of_keys && key == node.keys[i]) {
            return node.values[i];
        } else if (node.is_leaf) {
            return boost::none;
        } else {
            Node child(minimum_degree);
            read_node(child, node.children[i]);
            return find(child, key);
        }
    }


    // test functions

    void dump_child(const Node& node, int child_pos, int margin)
    {
        if (!node.is_leaf) {
            Node child(minimum_degree);
            const auto index = node.children[child_pos];
            read_node(child, index);
            dump_node(child, index, margin);
        }
    }

    void dump_node(const Node& node, NodeIndex index, int margin)
    {
        static int shift = 10;
        for (int i = 0; i < node.number_of_keys; ++i) {
            dump_child(node, i, margin + shift);
            std::cout << std::setw(margin) << std::right
                      << "i" << index << " "
                      << node.number_of_keys
                      << "(" << node.keys[i]
                      << ", " << node.values[i]
                      << ")"
                      << '\n';
        }
        dump_child(node, node.number_of_keys, margin + shift);
    }

    void merge_remove(
            Node& x, NodeIndex x_index,
            Node& left_child, NodeIndex left_child_index,
            Node& right_child, NodeIndex right_child_index,
            int pos_in_parent,
            const Key& key)
    {
        LOG_DEBUG("merge_remove()");
        left_child.keys[minimum_degree - 1] = x.keys[pos_in_parent];
        left_child.values[minimum_degree - 1] = x.values[pos_in_parent];
        for (int k = 0; k < minimum_degree - 1; ++k) {
            left_child.keys[k + minimum_degree] = right_child.keys[k];
            left_child.values[k + minimum_degree] = right_child.values[k];
            left_child.children[k + minimum_degree] = right_child.children[k];
        }
        left_child.children[2 * minimum_degree - 1] =
            right_child.children[minimum_degree - 1];
        left_child.number_of_keys = 2 * minimum_degree - 1;
        right_child.number_of_keys = 0;

        for (int k = pos_in_parent; k < x.number_of_keys - 1; ++k) {
            x.keys[k] = x.keys[k + 1];
            x.values[k] = x.values[k + 1];
            x.children[k + 1] = x.children[k + 2];
        }
        --x.number_of_keys;

        write_node(x, x_index);
        write_node(left_child, left_child_index);
        write_node(right_child, right_child_index);

        if (x.number_of_keys == 0) {
            assert(x_index == root_index);

            root_index = left_child_index;
            // change to swap
            root_node = left_child;
            remove(root_node, root_index, key);
        } else {
            remove(left_child, left_child_index, key);
        }
    }

    void remove(Node& x, NodeIndex x_index, const Key& key)
    {
        LOG_DEBUG("remove impl from node " << x_index);
        if (x.is_leaf) {
            remove_from_leaf(x, x_index, key);
        } else {
            remove_from_internal_node(x, x_index, key);
        }
    }

    void remove_from_leaf(Node& x, NodeIndex x_index, const Key& key)
    {
        LOG_DEBUG("remove from leaf");
        int i = 0;
        while (i < x.number_of_keys && key > x.keys[i]) {
            ++i;
        }
        if (i < x.number_of_keys && key == x.keys[i]) {
            for (int j = i + 1; j < x.number_of_keys; ++j) {
                x.keys[j - 1] = x.keys[j];
                x.values[j - 1] = x.values[j];
            }
            --x.number_of_keys;
            --total_number_of_keys;
            write_node(x, x_index);
        }
    }

    void remove_from_internal_node(Node& x, NodeIndex x_index, const Key& key)
    {
        int i = 0;
        while (i < x.number_of_keys && key > x.keys[i]) {
            ++i;
        }
        if (i < x.number_of_keys && key == x.keys[i]) {
            remove_found_key_in_internal_node(x, x_index, key, i);
        } else {
            remove_not_found_key_in_internal_node(x, x_index, key, i);
        }
    }

    void remove_found_key_in_internal_node(
            Node& x, NodeIndex x_index, const Key& key, int key_pos)
    {
        LOG_DEBUG("key in the current node");
        auto left_child_index = x.children[key_pos];
        Node left_child(minimum_degree);
        read_node(left_child, left_child_index);
        if (left_child.number_of_keys >= minimum_degree) {
            Key pred_key;
            Value pred_value;
            std::tie(pred_key, pred_value) = remove_max(
                    left_child, left_child_index);
            x.keys[key_pos] = pred_key;
            x.values[key_pos] = pred_value;
            write_node(x, x_index);
        } else {
            auto right_child_index = x.children[key_pos+1];
            Node right_child(minimum_degree);
            read_node(right_child, right_child_index);
            if (right_child.number_of_keys >= minimum_degree) {
                Key succ_key;
                Value succ_value;
                std::tie(succ_key, succ_value) = remove_min(
                        right_child, right_child_index);
                x.keys[key_pos] = succ_key;
                x.values[key_pos] = succ_value;
                write_node(x, x_index);
            } else {
                merge_remove(
                        x, x_index,
                        left_child, left_child_index,
                        right_child, right_child_index,
                        key_pos, key);
            }
        }
    }

    void remove_not_found_key_in_internal_node(
            Node& x, NodeIndex x_index, const Key& key, int child_pos)
    {
        LOG_DEBUG("key is not in the current node");
        Node child(minimum_degree);
        NodeIndex child_index = x.children[child_pos];
        read_node(child, child_index);
        if (child.number_of_keys > minimum_degree - 1) {
            LOG_DEBUG("remove from child");
            remove(child, child_index, key);
        } else {
            Node left_sibling(minimum_degree);
            NodeIndex left_sibling_index;

            LOG_DEBUG("try to read left sibling");
            bool has_left = has_child(
                    x, child_pos-1, left_sibling, left_sibling_index);
            LOG_FIELD(has_left);

            if (has_left && left_sibling.number_of_keys > minimum_degree - 1) {
                borrow_from_left(
                        x, x_index,
                        child, child_index,
                        left_sibling, left_sibling_index,
                        child_pos, key);
            } else {
                Node right_sibling(minimum_degree);
                NodeIndex right_sibling_index;

                LOG_DEBUG("try to read right sibling");
                LOG_DEBUG(child_pos);
                bool has_right = has_child(
                        x, child_pos + 1, right_sibling, right_sibling_index);
                LOG_FIELD(has_right);

                if (has_right && right_sibling.number_of_keys > minimum_degree - 1) {
                    borrow_from_right(
                            x, x_index,
                            child, child_index,
                            right_sibling, right_sibling_index,
                            child_pos, key);
                } else if (has_left) {
                    LOG_DEBUG("merge with left sibling");
                    merge_remove(
                            x, x_index,
                            left_sibling, left_sibling_index,
                            child, child_index,
                            child_pos - 1, key);
                } else {
                    LOG_DEBUG("merge with right sibling");
                    merge_remove(
                            x, x_index,
                            child, child_index,
                            right_sibling, right_sibling_index,
                            child_pos, key);
                }
            }
        }
    }

    void borrow_from_left(
            Node& x, NodeIndex x_index,
            Node& child, NodeIndex child_index,
            Node& left_sibling, NodeIndex left_sibling_index,
            int child_pos, const Key& key
    ) {
        LOG_DEBUG("\n3-a borrow from left");
        auto left_size = left_sibling.number_of_keys;

        insert_to_beginning(child, x.keys[child_pos-1], x.values[child_pos-1],
                            left_sibling.children[left_size]);

        // moving key from left sibling to x
        x.keys[child_pos-1] = left_sibling.keys[left_size - 1];
        x.values[child_pos-1] = left_sibling.values[left_size - 1];
        --left_sibling.number_of_keys;

        write_node(left_sibling, left_sibling_index);
        write_node(x, x_index);
        write_node(child, child_index);

        remove(child, child_index, key);
    }

    void borrow_from_right(
            Node& x, NodeIndex x_index,
            Node& child, NodeIndex child_index,
            Node& right_sibling, NodeIndex right_sibling_index,
            int child_pos, const Key& key
    ) {
        LOG_DEBUG("3-a borrow from right");
        auto right_size = right_sibling.number_of_keys;

        // insert key from x to the end of child
        child.keys[child.number_of_keys] = x.keys[child_pos];
        child.values[child.number_of_keys] = x.values[child_pos];
        child.children[child.number_of_keys + 1] = right_sibling.children[0];
        ++child.number_of_keys;

        // moving key from right sibling to x
        x.keys[child_pos] = right_sibling.keys[0];
        x.values[child_pos] = right_sibling.values[0];

        for (int j = 0; j < right_size - 1; ++j) {
            right_sibling.keys[j] = right_sibling.keys[j+1];
            right_sibling.values[j] = right_sibling.values[j+1];
            right_sibling.children[j] = right_sibling.children[j+1];
        }
        right_sibling.children[right_size - 1] = right_sibling.children[right_size];
        --right_sibling.number_of_keys;

        write_node(right_sibling, right_sibling_index);
        write_node(x, x_index);
        write_node(child, child_index);

        remove(child, child_index, key);
    }

    void insert_to_beginning(Node& x, const Key& key, const Value& value,
                             NodeIndex child_index)
    {
        assert(x.number_of_keys < 2 * minimum_degree - 1);
        for (int j = x.number_of_keys - 1; j >= 0; --j) {
            x.keys[j + 1] = x.keys[j];
            x.values[j + 1] = x.values[j];
            x.children[j + 2] = x.children[j + 1];
        }
        x.children[1] = x.children[0];

        x.keys[0] = key;
        x.values[0] = value;
        x.children[0] = child_index;

        ++x.number_of_keys;
    }

    bool has_child(const Node& x, int pos, Node& child, NodeIndex& child_index)
    {
        if (pos >= 0 && pos <= x.number_of_keys) {
            child_index = x.children[pos];
            read_node(child, child_index);
            return true;
        }
        return false;

    }

    std::tuple<Key, Value> remove_max(Node& x, NodeIndex x_index)
    {
        if (x.is_leaf) {
            auto index = x.number_of_keys - 1;
            auto key = x.keys[index];
            auto value = x.values[index];
            remove(x, x_index, key);
            return std::make_tuple(key, value);
        } else {
            Node n(minimum_degree);
            NodeIndex n_index = x.children[x.number_of_keys];
            read_node(n, n_index);
            return remove_max(n, n_index);
        }
    }

    std::tuple<Key, Value> remove_min(Node& x, NodeIndex x_index)
    {
        if (x.is_leaf) {
            auto index = 0;
            auto key = x.keys[index];
            auto value = x.values[index];
            remove(x, x_index, key);
            return std::make_tuple(key, value);
        } else {
            Node n(minimum_degree);
            NodeIndex n_index = x.children[0];
            read_node(n, n_index);
            return remove_min(n, n_index);
        }
    }

    size_t minimum_degree = 2;
    NodeIndex num_nodes = 1;
    NodeIndex root_index = 0;
    size_t total_number_of_keys = 0;

    NodeIndex current_node_index = 0;

    Node root_node;

    int io_operations_count = 0;

    std::fstream file;
};

template <class Key, class Value>
BTree<Key, Value>::BTree(size_t minimum_degree, const char* filename)
    : minimum_degree(minimum_degree)
    , root_node(minimum_degree)
    , file(filename, std::ios_base::out | std::ios_base::in |
                     std::ios_base::trunc | std::ios_base::binary)
{
    LOG_DEBUG("creating new btree. Filename is " << filename);
    if (minimum_degree < 2) {
        throw std::invalid_argument("minimum degree should be greater that 1");
    }
}

template <class Key, class Value>
BTree<Key, Value>::BTree(const char* filename)
    : file(filename, std::ios_base::out | std::ios_base::in | std::ios_base::binary)
{
    if (!file) {
        throw std::invalid_argument(std::string("can't open file ") + filename);
    }

    LOG_DEBUG("Reading existing btree. Filename is " << filename);

    read_metadata();

    root_node.resize(minimum_degree);
    read_node(root_node, root_index);
}

template <class Key, class Value>
BTree<Key, Value>::~BTree()
{
    write_metadata();
}

template <class Key, class Value>
boost::optional<Value> BTree<Key, Value>::find(const Key& key)
{
    return find(root_node, key);
}

template <class Key, class Value>
void BTree<Key, Value>::insert(const Key& key, const Value& value)
{
    LOG_DEBUG("INSERT key: " << key << ", value: " << value);
    if (root_node.is_full()) {
        LOG_DEBUG("root is full, splitting");
        Node new_root(minimum_degree);
        new_root.is_leaf = false;
        new_root.children[0] = root_index;

        NodeIndex new_root_index = num_nodes++;

        Node new_child(minimum_degree);
        split_child(new_root, new_root_index, 0, root_node, root_index,
                    new_child);
        insert_nonfull(new_root, root_index, key, value);

        LOG_DEBUG("after calling insert_nonfull")

        root_index = new_root_index;
        // TODO: replace by swap()
        root_node = new_root;
    } else {
        insert_nonfull(root_node, root_index, key, value);
    }
}

template <class Key, class Value>
void BTree<Key, Value>::remove(const Key& key)
{
    LOG_DEBUG("REMOVE " << key);
    remove(root_node, root_index, key);
}

template <class Key, class Value>
size_t BTree<Key, Value>::size()
{
    return total_number_of_keys;
}

} // namespace BTree
