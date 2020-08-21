#include <iostream>
#include <algorithm>
#include <fstream>
#include <cassert>
#include <vector>
#include <random>
#include <cmath>
#include <string>
#include <chrono>
#include <unordered_map>

using namespace std;
using namespace std::chrono;

auto randomNumberBetween = [](uint64_t low, uint64_t high) {
  auto randomFunc = [distribution_ = std::uniform_int_distribution<uint64_t>(low, high), random_engine_ = std::mt19937{std::random_device{}()}]() mutable {
    return distribution_(random_engine_);
  };
  return randomFunc;
};

static std::vector<uint64_t> generator(unsigned size, uint64_t maxValue) {
  vector<uint64_t> numbers;
  std::generate_n(std::back_inserter(numbers), size, randomNumberBetween(1, maxValue));
  sort(std::begin(numbers), std::end(numbers));
  return numbers;
}

unsigned computeLog(unsigned n) {
  auto lg = 0;
  while ((1u << lg) <= n) ++lg;
  return lg - 1;
}

class Builder {
public:
  size_t h;
  std::vector<unsigned> tree;
  Builder(size_t h) : h(h) {
    // Create a tree of size h
    tree.resize(1u << h);
    
    // And fill the tree
    fill(1, 1, h);
  }
  
  void fill(unsigned u, unsigned val, size_t h) {
    tree[u] = val;
    if (h == 1)
      return;
    // Reduce the height and continue with the same node
    auto currHeight = h / 2;
    auto restHeight = h - currHeight;
    fill(u, val, currHeight);

    // Then iterate over the remaining sons and also fill them
    auto sons = 1u << currHeight;
    auto firstNode = sons * u;
    auto firstVal = val + sons - 1;
    auto sonSize = (1u << restHeight) - 1;
    for (unsigned index = 0; index != sons; ++index)
      fill(firstNode + index, firstVal + index * sonSize, restHeight);
  }
  
  std::vector<std::vector<std::pair<uint32_t, uint32_t>>> extract() {
    // Extract the information about the periodical structure of the tree
    std::vector<std::vector<std::pair<uint32_t, uint32_t>>> sigma;
    sigma.resize(h + 1);
    sigma[1].resize(1);
    sigma[1][0] = make_pair(1, 0);
    for (unsigned index = 1; index != h; ++index) {
      // Iterate through each level and find the periods
      std::unordered_map<unsigned, std::pair<unsigned, unsigned>> info;
      unsigned start = (1u << index);
      
      // TODO: is it correct that we take 'diff' as the key? Maybe there are multiple 'diffs'?
      for (unsigned ptr = 0; ptr != index; ++ptr) {
        auto diff = tree[start + (1u << ptr)] - tree[start + (1u << ptr) - 1];
        auto it = info.find(diff);
        if (it == info.end())
          info[diff] = make_pair(ptr, tree[start + (1u << ptr)] - tree[start]);
      }
      
      // Gather and sort them by the period (first=period, second=partial sum of jumps)
      std::vector<pair<unsigned, unsigned>> store;
      for (auto elem : info)
        store.push_back(make_pair(elem.second.first, elem.second.second));
      std::sort(store.begin(), store.end());
      
      // And build the sigma table
      sigma[index + 1].resize(store.size() + 1);
      sigma[index + 1][0] = make_pair(tree[1u << index], store.size());
      for (unsigned ptr = 0; ptr != store.size(); ++ptr)
        sigma[index + 1][ptr + 1] = make_pair(store[store.size() - ptr - 1].first, store[store.size() - ptr - 1].second);
#if 0
      cout << "Level: " << (index + 1) << " starts with " << tree[1u << index] << ": ";
      for (auto elem : store)
        cout << "(period=" << elem.first << ", jump=" << elem.second << "), ";
      cout << endl;
#endif
    }
    return sigma;
  }
  
  void pretty() {
    // Print the tree onto the console
    auto buildTab = [&](size_t size) -> string {
      string ret;
      for (unsigned index = 0; index != size; ++index)
        ret += ' ';
      return ret;
    };
    auto tab = (1u << h) * std::to_string(tree[(1u << h) - 1]).size() / 2;
    for (unsigned index = 0; index != h; ++index) {
      auto tmp = buildTab(tab);
      for (unsigned u = (1u << index); u != (1u << (index + 1)); ++u) {
        cout << tmp;
        cout << tree[u];
        auto aux = buildTab(tab - std::to_string(tree[u]).size());
        cout << aux;
      }
      cout << endl;
      tab >>= 1;
    }
  }
};

class Solver {
public:
  size_t height, finalSize;
  std::vector<uint64_t> aux, tree;
  std::vector<std::vector<std::pair<unsigned, unsigned>>> sigma;
  Solver(std::vector<uint64_t>& val) {
    auto vectorSize = val.size(), totalSize = 2 * vectorSize;
    aux.resize(totalSize);
    tree.resize(totalSize);
    
    // Build the artificial BST
    for (unsigned index = vectorSize; index != totalSize; ++index)
      aux[index] = val[index - vectorSize];
    for (unsigned index = vectorSize - 1; index; --index)
      aux[index] = std::max(aux[2 * index], aux[2 * index + 1]);
    
    // And build the tree
    height = computeLog(totalSize);
    Builder builder(height);
    sigma = builder.extract();
    
    // Create the vam Emde Boas-layout
    for (unsigned h = height; h; --h)
      for (unsigned tmp = (1u << (h - 1)), u = tmp; u != (tmp << 1); ++u)
        tree[builder.tree[u]] = aux[u];
    
#if 0
    std::cerr << "DEBUG LAYOUT" << std::endl;
    for (unsigned index = 1; index != (1u << height); ++index) {
      std::cerr << index << " -> " << tree[index] << std::endl;
    }
#endif
    finalSize = totalSize;
#if 0
    for (unsigned index = 0; index != size; ++index)
      aux[index] = aux[index + size];
    aux.resize(size);
    std::cerr << "AUX: " << std::endl;
    for (auto elem : aux) {
      std::cerr << elem << ", ";
    }
    std::cerr << std::endl;
#else
    aux.clear();
#endif
  }
  
  uint64_t search(uint64_t key) {
    unsigned delta = 1, pi = 1;
    while (delta != height) {
      unsigned relative = 2 * pi - (1u << delta);
      unsigned left = sigma[delta + 1][0].first, capacity = sigma[delta + 1][0].second;
      for (unsigned index = 0; index != capacity && relative; ++index) {
        unsigned q = (relative >> sigma[delta + 1][index + 1].first);
        left += q * sigma[delta + 1][index + 1].second;
        relative -= (q << sigma[delta + 1][index + 1].first);
      }
      pi = (pi << 1) + (key > tree[left]);
      delta++;
    }
    return pi - finalSize; 
  }
};

// Final version of the cache-oblivious static search-tree
class BoostedSolver {
public:
  unsigned height, vectorSize, totalSize, width;
  std::vector<uint64_t> aux, tree;
  std::vector<std::vector<std::pair<unsigned, unsigned>>> sigma;
  std::vector<unsigned> compressed;
  BoostedSolver(std::vector<uint64_t>& val) {
    vectorSize = val.size(), totalSize = 2 * vectorSize;
    aux.resize(totalSize);
    tree.resize(totalSize);
    
    // Build the artificial BST
    for (unsigned index = vectorSize; index != totalSize; ++index)
      aux[index] = val[index - vectorSize];
    for (unsigned index = vectorSize - 1; index; --index)
      aux[index] = std::max(aux[2 * index], aux[2 * index + 1]);
    
    // And build the tree
    height = computeLog(totalSize);
    Builder builder(height);
    sigma = builder.extract();
    
    // Compute the maximal width of each level
    width = 0;
    for (unsigned level = 1; level <= height; ++level)
      width = std::max(width, sigma[level][0].second);
    
    // And compress the table, creating a single vector
    compressed.resize((1 + 2 * width) * (height - 1));
    auto currPtr = 0;
    for (unsigned level = 2; level <= height; ++level) {
      compressed[currPtr++] = sigma[level][0].first;
      for (unsigned index = 0; index != sigma[level][0].second; ++index) {
        compressed[currPtr + 2 * index] = sigma[level][index + 1].first;
        compressed[currPtr + 2 * index + 1] = sigma[level][index + 1].second;
      }
      currPtr += 2 * width;
    }
    assert(currPtr == compressed.size());
    
    // Create the vam Emde Boas-layout
    for (unsigned h = height; h; --h)
      for (unsigned tmp = (1u << (h - 1)), u = tmp; u != (tmp << 1); ++u)
        tree[builder.tree[u]] = aux[u];
    aux.clear();
  }
  
  uint64_t search(uint64_t key) {
    unsigned delta = 1, pi = 1, ptr = 0;
    // Iterate over all levels
    while (delta != height) {
      // Relative is the relative position within the current level
      unsigned relative = 2 * pi - (1u << delta);
      
      // Build 'left', by taking a look at each period, until the relative position becomes zero
      // Note that the periods are stored in descending order
      unsigned left = compressed[ptr++];
      for (unsigned index = 0; relative; ++index) {
        // On the even positions: the periods, on the odd ones: their partial sums 
        unsigned q = (relative >> compressed[ptr + 2 * index]);
        left += q * compressed[ptr + 2 * index + 1];
        relative -= (q << compressed[ptr + 2 * index]);
      }
      pi = (pi << 1) + (key > tree[left]);
      ptr += 2 * width;
      delta++;
    }
    
    // And return the position of the succesor within the array
    return pi - vectorSize; 
  }
};

class BinarySearch {
public:
  std::vector<uint64_t> v;
  BinarySearch(std::vector<uint64_t> val) : v(move(val)) {}
  
  uint64_t search(uint64_t key) {
    return std::lower_bound(v.begin(), v.end(), key) - v.begin();
  }
};

int main(void) {
  // Run: ./a.out < 1.in
  ifstream input("1.in");
  unsigned n, type;
  input >> n >> type;
  
  static constexpr unsigned runs = 1e7;
  static constexpr uint64_t inf = 1e8;
  auto analyze = [&]() -> std::vector<uint64_t> {
    if (!type) {
      vector<uint64_t> val = generator(n, true ? std::numeric_limits<uint64_t>::max() : inf);
      assert(std::is_sorted(val.begin(), val.end()));
      auto it = std::unique(val.begin(), val.end());
      size_t dist = it - val.begin(), complete = val.size();
      while (dist != complete) {
        val.pop_back();
        ++dist;
      }
      auto lg = computeLog(val.size());
      val.resize(1u << lg);
      assert(std::is_sorted(val.begin(), val.end()));
      return val;
    }
    
    std::vector<uint64_t> val(n);
    for (unsigned index = 0; index != n; ++index)
      input >> val[index];
    return val;
  };
  
  vector<uint64_t> val = analyze();
  std::cerr << "size=" << val.size() << std::endl;
  std::cerr << "min=" << val.front() << " max=" << val.back() << std::endl;
  Solver solver(val);
  BoostedSolver boostedSolver(val);
  BinarySearch bs(val);
#if 0
  for (unsigned index = val.front(); index <= val.back(); ++index) {
    auto query = index;
    if (boostedSolver.search(query) != bs.search(query)) {
      std::cerr << "key=" << query << ": boosted=" << boostedSolver.search(query) << " vs solver=" << solver.search(query) << " vs bs=" << bs.search(query) << std::endl;
      assert(0);
    }
  }
#elif 1
  vector<uint64_t> queries;
  std::generate_n(std::back_inserter(queries), runs, randomNumberBetween(val.front(), val.back()));
  auto measureTime = [&](unsigned mode) -> double {
    auto start = high_resolution_clock::now();
    if (mode == 2) {
      for (auto query : queries)
        boostedSolver.search(query);
    } else if (mode == 1) {
      for (auto query : queries)
        solver.search(query);
    } else {
      for (auto query : queries)
        bs.search(query);
    }
    auto stop = high_resolution_clock::now();
    return duration_cast<nanoseconds>(stop - start).count() / queries.size();
  };
  
  double boostedTime = measureTime(2), solverTime = measureTime(1), bsTime = measureTime(0);
  std::cerr << "BoostedSolver: " << boostedTime <<  ", Solver: " << solverTime << ", BS: " << bsTime << std::endl;
  
  auto checkForCorrectness = [&]() -> bool {
    std::sort(queries.begin(), queries.end());
    for (auto query : queries) {
      if (boostedSolver.search(query) != bs.search(query)) {
        std::cerr << "key=" << query << ": boosted=" << boostedSolver.search(query) << " vs solver=" << solver.search(query) << " vs bs=" << bs.search(query) << std::endl;
        return false;
      }
    }
    return true;
  };
  assert(checkForCorrectness());
#else 
  Builder builder(6);
  builder.extract();
#endif
  return 0;
}
