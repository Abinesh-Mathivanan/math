#include <iostream>
#include <vector>

long long partition(int n) {
    if (n < 0) return 0;

    // Initialdddize a vector to store partition counts
    std::vector<long long> partitions(n + 1, 0);
    partitions[0] = 1;

    // Calculate partition counts using the Hardy-Ramanujan-Rademacher formula
    for (int i = 1; i <= n; ++i) {
        int j = 1;
        int k = 1;
        partitions[i] = 0;

        while (j * (3 * j - 1) / 2 <= i) {
            partitions[i] += partitions[i - j * (3 * j - 1) / 2] * k;
            j += 1;

            if (j >= 0) k *= -1;
        }
   
        j = 1;
        k = 1;

        while (j * (3 * j + 1) / 2 <= i) {
            partitions[i] += partitions[i - j * (3 * j + 1) / 2] * k;
            j += 1;

            if (j >= 0) k *= -1;
        }
    }

    return partitions[n];
}

int main() {
    int n;
    std:: cout << "Enter a number: ";
    std:: cin >> n;
    long long partition_count = partition(n);
    std::cout << "The number of partitions of " << n << " is " << partition_count << std::endl;
    return 0;
}
