// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2026 Rasmus Munk Larsen <rmlarsen@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
// SPDX-License-Identifier: MPL-2.0

#ifndef EIGEN_LRU_CACHE_H
#define EIGEN_LRU_CACHE_H

#include <cstddef>
#include <functional>
#include <list>
#include <unordered_map>
#include <utility>

#include "../InternalHeaderCheck.h"

namespace Eigen {
namespace internal {

// Bounded least-recently-used cache.
//
// Keyed on Key (hashed by Hash, compared by KeyEqual) and storing movable
// Value objects. find() and insert() are O(1) average: an unordered_map
// lookup plus a constant-time std::list splice. On hit, the touched entry
// is promoted to the front of the list. On insert into a full cache, the
// back of the list (least-recently-used) is destroyed before the new entry
// is added — Value's destructor runs, so an RAII Value handles eviction
// cleanup without any additional callback machinery.
//
// Thread safety: none. Callers must serialize.
//
// Not intended for hot-loop O(n) workloads where n is large; for the small
// caches Eigen uses (handful of entries, low shape cardinality), the
// node-allocation overhead of std::list / std::unordered_map is amortized
// across many hits per insert.
template <typename Key, typename Value, typename Hash = std::hash<Key>, typename KeyEqual = std::equal_to<Key>>
class LruCache {
 public:
  using key_type = Key;
  using value_type = Value;

  explicit LruCache(std::size_t capacity) : capacity_(capacity) {
    eigen_assert(capacity_ > 0 && "LruCache capacity must be positive");
    // Pre-size the bucket array so the map never rehashes while filling. A
    // rehash is otherwise guaranteed once load factor crosses ~1.0.
    index_.reserve(capacity_);
  }

  LruCache(const LruCache&) = delete;
  LruCache& operator=(const LruCache&) = delete;

  LruCache(LruCache&& o) noexcept : capacity_(o.capacity_), items_(std::move(o.items_)), index_(std::move(o.index_)) {}

  LruCache& operator=(LruCache&& o) noexcept {
    if (this != &o) {
      capacity_ = o.capacity_;
      items_ = std::move(o.items_);
      index_ = std::move(o.index_);
    }
    return *this;
  }

  ~LruCache() = default;

  // Returns a pointer to the cached value and marks it most-recently-used,
  // or nullptr if the key is not present.
  Value* find(const Key& key) {
    auto map_it = index_.find(key);
    if (map_it == index_.end()) return nullptr;
    items_.splice(items_.begin(), items_, map_it->second);
    return &map_it->second->second;
  }

  // Inserts (key, value), evicting the least-recently-used entry if the cache
  // is at capacity. If the key already exists, the existing entry is destroyed
  // and a new entry is inserted as most-recently-used.
  // Returns a pointer to the inserted value, or nullptr if capacity is 0
  // (degenerate case — assert-firing in debug; the early return keeps release
  // builds from dereferencing items_.back() on an empty list).
  Value* insert(const Key& key, Value value) {
    if (capacity_ == 0) return nullptr;
    auto map_it = index_.find(key);
    if (map_it != index_.end()) {
      auto old_it = map_it->second;
      items_.emplace_front(key, std::move(value));
      map_it->second = items_.begin();
      items_.erase(old_it);
      return &map_it->second->second;
    }
    if (items_.size() >= capacity_) {
      index_.erase(items_.back().first);
      items_.pop_back();
    }
    items_.emplace_front(key, std::move(value));
    index_.emplace(key, items_.begin());
    return &items_.front().second;
  }

  void clear() {
    items_.clear();
    index_.clear();
  }

  std::size_t size() const { return items_.size(); }
  std::size_t capacity() const { return capacity_; }
  bool empty() const { return items_.empty(); }

 private:
  using list_type = std::list<std::pair<Key, Value>>;

  std::size_t capacity_;
  list_type items_;
  std::unordered_map<Key, typename list_type::iterator, Hash, KeyEqual> index_;
};

}  // namespace internal
}  // namespace Eigen

#endif  // EIGEN_LRU_CACHE_H
