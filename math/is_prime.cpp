// O(sqrt(n))
bool is_prime(int n) {
    if (n <= 1)
        return 0;

    if (n == 2 or n == 3)
        return 1;

    if (n % 2 == 0 or n % 3 == 0)
        return 0;

    for (int i = 5; i * i <= n; i += 6)
        if (n % i == 0 or n % (i + 2) == 0)
            return 0;

    return 1;
}