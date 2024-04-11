#include <vector>
#include <map>
#include <set>
#include <iostream>

void dbg_out() {
    std::cerr << '\n';
}

template<typename H, typename... T>
void dbg_out(const H h, const T... t) {
    std::cerr << ' ' << h; dbg_out(t...);
}

template<typename T, typename U>
void dbg_out(const std::pair<T, U>& a) {
    std::cerr << " {" << a.first << ", " << a.second << "}\n";
}

template<typename T, std::size_t N>
void dbg_out(const std::array<T, N>& a) {
    for (const T& b : a)
        std::cerr << ' ' << b;
    std::cerr << '\n';
}

template<typename T>
void dbg_out(const std::vector<T>& a) {
    for (const T& b : a)
        std::cerr << ' ' << b;
    std::cerr << '\n';
}

template<typename T, typename U>
void dbg_out(const std::vector<std::pair<T, U>>& a) {
    for (const std::pair<T, U>& b : a)
        std::cerr << " {" << b.first << ", " << b.second << '}';
    std::cerr << '\n';
}

template<typename T, std::size_t N> 
void dbg_out(const std::vector<std::array<T, N>>& a) {
    std::cerr << '\n';
    for (int i = 0; i < a.size(); i++) {
        std::cerr << '{' << i << "}:";
        dbg_out(a[i]);
    }
}

template<typename T>
void dbg_out(const std::vector<std::vector<T>>& a) {
    std::cerr << '\n';
    for (int i = 0; i < a.size(); i++) {
        std::cerr << '{' << i << "}:";
        dbg_out(a[i]);
    }
}

template<typename T, typename U>
void dbg_out(const std::map<T, U>& a) {
    for (const std::pair<T, U>& b : a)
        std::cerr << " {" << b.first << ", " << b.second << '}';
    std::cerr << '\n';
}

template<typename T>
void dbg_out(const std::set<T>& a) {
    for (const T& b : a)
        std::cerr << ' ' << b;
    std::cerr << '\n';
}

#define dbg(...) std::cerr << '[' << __LINE__ << "] (" << #__VA_ARGS__ << "):", dbg_out(__VA_ARGS__), std::cerr << '\n';