#ifndef HEAP_H
#define HEAP_H
#include <cmath>
#include <vector>

template<typename T> 
class Heap {
public:
    // the nodes of the heap
    struct heap_reference {
        int pos;
        T* item;
        int key;

        heap_reference() : pos(-1), item(0), key(0) {}

        heap_reference(int p, T* i, int k) :
            pos(p), item(i), key(k) {}

        bool is_valid() { return pos >= 0; }
        void invalidate() { pos = -1; }
    };

    Heap() : heap() {};

    ~Heap() {};

    int size() { return heap.size(); }

    // remove and return the item with the min key
    T* pop() {
        auto top = heap[0].item;
        auto last = heap.back();
        heap.pop_back();

        if (heap.size() > 0) {
            heap[0] = last;
            siftdown(0);
        }
        return top;
    }

    // insert a new item with the given key
    heap_reference insert(T* n, int key) {
        heap.push_back(heap_reference(heap.size()-1, n, key));
        return siftup(heap.size()-1);
    }

    void increase_key(heap_reference n, int new_key) {
        heap[n.pos].key = new_key;
        siftup(n.pos);
    }

private:
    std::vector<heap_reference> heap;

    heap_reference siftup(int hole) {
        auto moving = heap[hole];

        auto p = heap_parent(hole);
        while (p != 0 && heap[hole].key > moving.key) {
            heap[hole] = heap[p];
            heap[hole].pos = hole;
            hole = p;
            p = heap_parent(hole);
        }
        heap[hole] = moving;
        heap[hole].pos = hole;
        return heap[hole];
    }

    void siftdown(int hole) {
        auto moving = heap[hole];
        auto c = minchild(hole);
        while (c > 0 && heap[c].key < moving.key) {
            heap[hole] = heap[c];
            heap[hole].pos = hole;
            hole = c;
            c = minchild(c);
        }
        heap[hole] = moving;
        moving.pos = hole;
    }

    int heap_parent(int pos) {
        if (pos == 0) return -1;
        return (pos - 1) / 2;
    }

    int minchild(int pos) {
        auto left = 2*pos + 1;
        auto right = 2*pos + 2;
        // both children exist
        if (right < heap.size()) {
            // return the smallest
            return (heap[right].key < heap[left].key) ? right : left;

        // left child is the only child, so it's the smallest
        } else if (left < heap.size()) {
            return left;

        // no children
        } else return -1;
    }

};
#endif

