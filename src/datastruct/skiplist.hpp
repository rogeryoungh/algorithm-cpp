#ifndef ALGO_DATASTRUCT_SKIPLIST
#define ALGO_DATASTRUCT_SKIPLIST

#include "../base.hpp"

#include <array>
#include <optional>
#include <queue>
#include <random>

template <typename K, typename V>
struct SkipList {
  enum : u64 {
    MAXL = 32,                                 // 最大层数
    P = 4,                                     // 概率（1 / 4)
    S = 0xFFFF,                                // 限制随机值
    PS = S / P,                                // 概率
    INVALID = std::numeric_limits<u64>::max(), // 最大值（哨兵）
  };
  /*
      跳表结点（层数，值，指针）
  */
  struct SkipListNode {
    K key;
    V value;
    i32 level{};
    SkipListNode **forward;
    SkipListNode() = default;
    SkipListNode(K k, V v, i32 l, SkipListNode *nxt = nullptr) {
      key = k;
      value = v;
      level = l;
      forward = new SkipListNode *[l + 1];
      for (i32 i = 0; i <= l; ++i)
        forward[i] = nxt;
    }
    ~SkipListNode() {
      delete[] forward;
    }
  };
  using Node = SkipListNode;
  Node *head, *tail; // 两排哨兵
  i32 length;        // 链表长度 L0 层
  i32 level;         // 层数
  std::mt19937 rng;
  SkipList() : rng(std::random_device{}()) {
    level = length = 0;
    tail = new Node(INVALID, 0, 0);
    head = new Node(INVALID, 0, MAXL, tail);
  }
  ~SkipList() {
    delete head;
    delete tail;
  }
  /*
      第一层 50 % ， 第二层 25 % ，第三层 12.5% 依此类推， P = 1 / 4
  */
  i32 randomLevel() {
    i32 lv = 1;
    while ((rng() & S) < PS)
      ++lv;
    return MAXL > lv ? lv : MAXL;
  }
  /*
      插入新结点，
  */
  V &insert(const K &key, const V &value) {
    std::array<Node *, MAXL + 1> update{};
    Node *p = head;
    for (i32 i = level; i >= 0; --i) {
      while (p->forward[i]->key < key) {
        p = p->forward[i];
      }
      update[i] = p;
    }
    p = p->forward[0];
    if (p->key == key) {
      return p->value = value;
    }
    i32 lv = randomLevel();
    if (lv > level) { // 一次只会比当前层数高一层
      lv = ++level;
      update[lv] = head;
    }
    Node *newNode = new Node(key, value, lv);
    for (i32 i = lv; i >= 0; --i) {
      p = update[i];
      newNode->forward[i] = p->forward[i];
      p->forward[i] = newNode;
    }
    ++length;
    return newNode->value;
  }
  bool erase(const K &key) {
    std::array<Node *, MAXL + 1> update;
    Node *p = head;
    for (i32 i = level; i >= 0; --i) {
      while (p->forward[i]->key < key) {
        p = p->forward[i];
      }
      update[i] = p;
    }
    p = p->forward[0];
    if (p->key != key)
      return false;
    for (i32 i = 0; i <= level; ++i) {
      if (update[i]->forward[i] != p) {
        break;
      }
      update[i]->forward[i] = p->forward[i];
    }
    delete p;
    while (level > 0 && head->forward[level] == tail)
      --level;
    --length;
    return true;
  }
  std::optional<V> find(const K &key) {
    Node *p = head;
    for (i32 i = level; i >= 0; --i) {
      while (p->forward[i]->key < key) {
        p = p->forward[i];
      }
    }
    p = p->forward[0];
    if (p->key == key)
      return p->value;
    return std::nullopt;
  }
  V &operator[](const K &key) {
    Node *p = head;
    for (i32 i = level; i >= 0; --i) {
      while (p->forward[i]->key < key) {
        p = p->forward[i];
      }
    }
    p = p->forward[0];
    if (p->key == key) {
      return p->value;
    } else {
      return insert(key, 0);
    }
  }
  bool count(const K &key) {
    return find(key) != tail->value;
  }
};

#endif
